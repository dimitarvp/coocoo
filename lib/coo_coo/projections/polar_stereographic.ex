defmodule CooCoo.Projections.PolarStereographic do
  @moduledoc """
  Provides Polar Stereographic projection conversions.
  This module implements the core mathematical formulas for the projection.
  UPS is a specific application of Polar Stereographic.
  """

  alias CooCoo.Constants
  alias CooCoo.Coordinates.GeodeticCoordinates
  alias CooCoo.Coordinates.MapProjectionCoordinates
  import CooCoo.MathHelpers, only: [equal?: 3, zero?: 2, finite_float?: 1]

  @boundary_epsilon 1.0e-15

  def polar_pow(eccentricity, sin_latitude) do
    if not (is_number(eccentricity) and is_number(sin_latitude) and
              abs(sin_latitude) <= 1.0 and eccentricity >= 0 and eccentricity < 1.0) do
      {:error, :invalid_input_to_polar_pow}
    else
      es_sin = eccentricity * sin_latitude

      if abs(es_sin) > 1.0 or equal?(abs(es_sin), 1.0, @boundary_epsilon) do
        {:error, :es_sin_out_of_bounds_for_polar_pow}
      else
        denominator = 1.0 + es_sin
        numerator = 1.0 - es_sin

        # Numerator can be 0 if es_sin is 1 (caught by abs(es_sin) < 1)
        if numerator < 0.0 or equal?(denominator, 0.0, @boundary_epsilon) do
          {:error, :invalid_base_for_polar_pow}
        else
          base = numerator / denominator
          exponent = eccentricity / 2.0
          {:ok, :math.pow(base, exponent)}
        end
      end
    end
  end

  defp calculate_rho_forward(a, es, t_val, projection_params) do
    standard_parallel_param = Keyword.get(projection_params, :standard_parallel)
    scale_factor_param = Keyword.get(projection_params, :scale_factor)

    cond do
      sp = standard_parallel_param ->
        eff_sp = abs(sp)

        if equal?(abs(eff_sp), Constants.pi_over_2(), @boundary_epsilon) do
          k90 = :math.sqrt(:math.pow(1.0 + es, 1.0 + es) * :math.pow(1.0 - es, 1.0 - es))

          if zero?(k90, @boundary_epsilon),
            do: {:error, :k90_is_zero_rho_fwd_sp_pole},
            else: {:ok, 2.0 * a * t_val / k90}
        else
          sinolat_eff_sp = :math.sin(eff_sp)
          essin_eff_sp = es * sinolat_eff_sp
          sqrt_arg_mc = 1.0 - essin_eff_sp * essin_eff_sp

          if sqrt_arg_mc <= 0.0 do
            {:error, :sqrt_domain_error_for_mc_rho_fwd}
          else
            mc_val = :math.cos(eff_sp) / :math.sqrt(sqrt_arg_mc)
            a_mc_val = a * mc_val

            # Nested `with` for polar_pow result
            with {:ok, pow_es_eff_sp} <- polar_pow(es, sinolat_eff_sp),
                 true <-
                   (finite_float?(pow_es_eff_sp) and pow_es_eff_sp > 0.0) or
                     {:error, {:polar_pow_eff_sp_failed_rho_fwd, pow_es_eff_sp}} do
              tc_val = :math.tan(Constants.pi_over_4() - eff_sp / 2.0) / pow_es_eff_sp

              if zero?(tc_val, @boundary_epsilon),
                do: {:error, :tc_val_is_zero_rho_fwd},
                else: {:ok, a_mc_val * t_val / tc_val}
            else
              # Handles errors from the nested `with` (polar_pow or its validation)
              error_reason ->
                {:error, {:rho_fwd_std_parallel_calc_failed_inner_with, error_reason}}
            end
          end
        end

      sf = scale_factor_param ->
        k90 = :math.sqrt(:math.pow(1.0 + es, 1.0 + es) * :math.pow(1.0 - es, 1.0 - es))

        if zero?(k90, @boundary_epsilon),
          do: {:error, :k90_is_zero_rho_fwd_sf_mode},
          else: {:ok, 2.0 * a * sf * t_val / k90}

      true ->
        {:error, :missing_standard_parallel_or_scale_factor_for_rho_fwd}
    end
  end

  defp calculate_t_val_inverse(a, es, rho, projection_params) do
    standard_parallel_param = Keyword.get(projection_params, :standard_parallel)
    scale_factor_param = Keyword.get(projection_params, :scale_factor)

    cond do
      sp = standard_parallel_param ->
        eff_sp = abs(sp)

        if equal?(abs(eff_sp), Constants.pi_over_2(), @boundary_epsilon) do
          k90 = :math.sqrt(:math.pow(1.0 + es, 1.0 + es) * :math.pow(1.0 - es, 1.0 - es))
          divisor = 2.0 * a

          if zero?(k90, @boundary_epsilon) or zero?(divisor, @boundary_epsilon) do
            {:error, :zero_divisor_for_t_val_inv_sp_pole}
          else
            {:ok, rho * k90 / divisor}
          end
        else
          sinolat_eff_sp = :math.sin(eff_sp)
          essin_eff_sp = es * sinolat_eff_sp
          sqrt_arg_mc = 1.0 - essin_eff_sp * essin_eff_sp

          if sqrt_arg_mc <= 0.0 do
            {:error, :sqrt_domain_error_for_mc_t_inv}
          else
            mc_val = :math.cos(eff_sp) / :math.sqrt(sqrt_arg_mc)
            a_mc_val = a * mc_val

            if zero?(a_mc_val, @boundary_epsilon) do
              {:error, :a_mc_val_is_zero_t_inv}
            else
              with {:ok, pow_es_eff_sp} <- polar_pow(es, sinolat_eff_sp),
                   true <-
                     (finite_float?(pow_es_eff_sp) and pow_es_eff_sp > 0.0) or
                       {:error, {:polar_pow_eff_sp_failed_t_inv, pow_es_eff_sp}} do
                tc_val = :math.tan(Constants.pi_over_4() - eff_sp / 2.0) / pow_es_eff_sp

                # tc_val could be zero if eff_sp is pi/2, but that case is handled by the outer `if`
                {:ok, rho * tc_val / a_mc_val}
              else
                error_reason ->
                  {:error, {:t_inv_std_parallel_calc_failed_inner_with, error_reason}}
              end
            end
          end
        end

      sf = scale_factor_param ->
        k90 = :math.sqrt(:math.pow(1.0 + es, 1.0 + es) * :math.pow(1.0 - es, 1.0 - es))
        divisor_sf = 2.0 * a * sf

        if zero?(k90, @boundary_epsilon) or zero?(divisor_sf, @boundary_epsilon) do
          {:error, :zero_divisor_for_t_val_inv_sf_mode}
        else
          {:ok, rho * k90 / divisor_sf}
        end

      true ->
        {:error, :missing_standard_parallel_or_scale_factor_for_t_inv}
    end
  end

  @doc """
  Converts Geodetic coordinates to Polar Stereographic coordinates.
  Returns `{:ok, %MapProjectionCoordinates{}}` or `{:error, reason}`.
  """
  def convert_from_geodetic(
        %GeodeticCoordinates{} = geodetic_coords,
        projection_params,
        ellipsoid_params
      ) do
    # Extract parameters first
    a = Keyword.get(ellipsoid_params, :a, Constants.wgs84_a())
    es = Keyword.get(ellipsoid_params, :es, Constants.wgs84_es())
    # Renamed to avoid conflict
    central_meridian_param = Keyword.fetch!(projection_params, :central_meridian)
    false_easting = Keyword.fetch!(projection_params, :false_easting)
    false_northing = Keyword.fetch!(projection_params, :false_northing)

    # Determine operational hemisphere
    is_southern_hemisphere =
      if hem = Keyword.get(projection_params, :hemisphere) do
        hem == ?S
      else
        sp_param = Keyword.get(projection_params, :standard_parallel, geodetic_coords.latitude)
        sp_param < 0.0
      end

    # Effective latitude, longitude, and central_meridian for N-Pole centric math
    latitude_eff =
      if is_southern_hemisphere, do: -geodetic_coords.latitude, else: geodetic_coords.latitude

    longitude_eff =
      if is_southern_hemisphere, do: -geodetic_coords.longitude, else: geodetic_coords.longitude

    op_central_meridian =
      if is_southern_hemisphere, do: -central_meridian_param, else: central_meridian_param

    # Now proceed with the `with` block using these pre-calculated values
    # Ensure a is available if not directly used but needed for helpers
    with {:ok, a_val} <- {:ok, a},
         # Ensure es is available
         {:ok, es_val} <- {:ok, es},
         # Check if at the pole, which is a simple calculation not needing `with` for its condition
         # Special tuple to signal early successful exit
         true <-
           abs(abs(latitude_eff) - Constants.pi_over_2()) >= 1.0e-10 or
             {:bypass_with_pole_coords,
              {:ok,
               %MapProjectionCoordinates{
                 easting: false_easting,
                 northing: false_northing,
                 warning_message: ""
               }}},
         {:ok, dlam} <-
           {:ok,
            (
              dlam_raw = longitude_eff - op_central_meridian

              cond do
                dlam_raw > Constants.pi() -> dlam_raw - Constants.two_pi()
                dlam_raw < -Constants.pi() -> dlam_raw + Constants.two_pi()
                true -> dlam_raw
              end
            )},
         {:ok, sin_lat_eff} <- {:ok, :math.sin(latitude_eff)},
         # Pass es_val
         {:ok, pow_es_val} <- polar_pow(es_val, sin_lat_eff),
         true <-
           (finite_float?(pow_es_val) and pow_es_val > 0.0) or
             {:error, {:polar_pow_failed_or_non_positive, pow_es_val}},
         # Added to prevent division by zero for t_val
         true <-
           not zero?(pow_es_val, @boundary_epsilon) or
             {:error, :polar_pow_result_is_zero_division_error},
         {:ok, t_val} <-
           {:ok, :math.tan(Constants.pi_over_4() - latitude_eff / 2.0) / pow_es_val},
         # Pass a_val, es_val
         {:ok, rho} <- calculate_rho_forward(a_val, es_val, t_val, projection_params) do
      temp_easting = rho * :math.sin(dlam)
      temp_northing = -rho * :math.cos(dlam)

      {easting_final, northing_final} =
        if is_southern_hemisphere do
          {-temp_easting + false_easting, -temp_northing + false_northing}
        else
          {temp_easting + false_easting, temp_northing + false_northing}
        end

      {:ok,
       %MapProjectionCoordinates{
         easting: easting_final,
         northing: northing_final,
         warning_message: ""
       }}
    else
      # This block handles errors from the `with` clause OR the special bypass tuple
      # Pass through the success
      {:bypass_with_pole_coords, ok_pole_coords_tuple} -> ok_pole_coords_tuple
      {:error, error_details} -> {:error, {:convert_from_geodetic_ps_failed, error_details}}
    end
  end

  @doc """
  Converts Polar Stereographic coordinates to Geodetic coordinates.
  Returns `{:ok, %GeodeticCoordinates{}}` or `{:error, reason}`.
  """
  def convert_to_geodetic(
        %MapProjectionCoordinates{} = proj_coords,
        projection_params,
        ellipsoid_params
      ) do
    a = Keyword.fetch!(ellipsoid_params, :a)
    es = Keyword.fetch!(ellipsoid_params, :es)

    central_meridian = Keyword.fetch!(projection_params, :central_meridian)
    false_easting = Keyword.fetch!(projection_params, :false_easting)
    false_northing = Keyword.fetch!(projection_params, :false_northing)

    # Determine hemisphere_char first
    hemisphere_char =
      Keyword.get_lazy(projection_params, :hemisphere, fn ->
        if proj_coords.northing >= false_northing, do: ?N, else: ?S
      end)

    is_southern_hemisphere = hemisphere_char == ?S
    op_central_meridian = if is_southern_hemisphere, do: -central_meridian, else: central_meridian

    dx = proj_coords.easting - false_easting
    dy = proj_coords.northing - false_northing

    {dx_calc, dy_calc} =
      if is_southern_hemisphere, do: {-dx, -dy}, else: {dx, dy}

    rho = :math.sqrt(dx_calc * dx_calc + dy_calc * dy_calc)

    if zero?(rho, @boundary_epsilon) do
      latitude_at_pole =
        if is_southern_hemisphere, do: -Constants.pi_over_2(), else: Constants.pi_over_2()

      {:ok,
       %GeodeticCoordinates{
         latitude: latitude_at_pole,
         longitude: central_meridian,
         warning_message: ""
       }}
    else
      # Pass the explicitly determined hemisphere_char to calculate_t_val_inverse if it needs it,
      # or ensure projection_params for it are set up correctly.
      # The current calculate_t_val_inverse uses standard_parallel or scale_factor from projection_params.
      # We might need to adjust how it infers hemisphere if :hemisphere key isn't present in projection_params for it.
      # For now, assuming projection_params passed to calculate_t_val_inverse is sufficient.
      # The `hemisphere_char` variable is used to determine `is_southern_hemisphere` which influences `op_central_meridian`.
      # The call to `calculate_t_val_inverse` does not directly take `hemisphere_char` but uses `projection_params`.

      with {:ok, t_val} <- calculate_t_val_inverse(a, es, rho, projection_params) do
        phi_initial = Constants.pi_over_2() - 2.0 * :math.atan(t_val)

        latitude_eff_result =
          Enum.reduce_while(1..10, phi_initial, fn _iter, current_phi ->
            with {:ok, pow_es_val} <- polar_pow(es, :math.sin(current_phi)),
                 true <-
                   finite_float?(pow_es_val) or
                     {:halt, {:error, {:polar_pow_failed_in_phi_iteration, pow_es_val}}} do
              next_phi = Constants.pi_over_2() - 2.0 * :math.atan(t_val * pow_es_val)

              if abs(next_phi - current_phi) < 1.0e-10 do
                {:halt, {:ok, next_phi}}
              else
                {:cont, next_phi}
              end
            else
              error_reason -> {:halt, {:error, {:phi_iteration_polar_pow_failed, error_reason}}}
            end
          end)

        case latitude_eff_result do
          {:ok, latitude_eff} ->
            longitude_eff = op_central_meridian + :math.atan2(dx_calc, -dy_calc)

            latitude = if is_southern_hemisphere, do: -latitude_eff, else: latitude_eff
            longitude = if is_southern_hemisphere, do: -longitude_eff, else: longitude_eff

            longitude =
              cond do
                longitude > Constants.pi() -> longitude - Constants.two_pi()
                longitude < -Constants.pi() -> longitude + Constants.two_pi()
                true -> longitude
              end

            latitude =
              cond do
                latitude > Constants.pi_over_2() -> Constants.pi_over_2()
                latitude < -Constants.pi_over_2() -> -Constants.pi_over_2()
                true -> latitude
              end

            longitude =
              cond do
                longitude > Constants.pi() -> Constants.pi()
                longitude < -Constants.pi() -> -Constants.pi()
                true -> longitude
              end

            {:ok,
             %GeodeticCoordinates{
               longitude: longitude,
               latitude: latitude,
               height: 0.0,
               warning_message: ""
             }}

          {:error, reason} ->
            {:error, {:latitude_calculation_failed_in_to_geodetic, reason}}
        end
      else
        # This else matches errors from calculate_t_val_inverse
        {:error, error_details} -> {:error, {:convert_to_geodetic_ps_failed, error_details}}
      end
    end
  end
end
