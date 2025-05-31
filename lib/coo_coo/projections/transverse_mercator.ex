defmodule CooCoo.Projections.TransverseMercator do
  @moduledoc """
  Provides Transverse Mercator projection conversions.
  This module implements the core mathematical formulas for the projection.
  UTM is a specific application of Transverse Mercator.

  Based on NGA.SIG.0012_2.0.0_UTMUPS "The Universal Grids and the Transverse
  Mercator and Polar Stereographic Map Projections".
  """

  alias CooCoo.Constants
  alias CooCoo.Coordinates.GeodeticCoordinates
  alias CooCoo.Coordinates.MapProjectionCoordinates
  alias CooCoo.Projections.TMHelpers

  import CooCoo.MathHelpers, only: [equal?: 3, zero?: 2, finite_float?: 1]

  # Value taken from TMHelpers for consistency in this module
  @n_terms 6
  @max_delta_long Constants.pi_over_180() * 70.0

  # Coefficients are pre-calculated for WGS84 for potential optimization if needed,
  # but the conversion functions will use get_tm_coefficients for flexibility.
  @wgs84_coefficients_tuple TMHelpers.generate_coefficients(Constants.wgs84_f())
  @wgs84_a_coeffs elem(@wgs84_coefficients_tuple, 0)
  @wgs84_b_coeffs elem(@wgs84_coefficients_tuple, 1)
  @wgs84_r4oa elem(@wgs84_coefficients_tuple, 2)

  @boundary_epsilon 1.0e-15

  def n_terms(), do: @n_terms

  # --- Helper for ArcHyperbolicTangent ---
  def atanh(x) when is_number(x) do
    cond do
      equal?(x, 1.0, @boundary_epsilon) ->
        {:ok, :infinity}

      equal?(x, -1.0, @boundary_epsilon) ->
        {:ok, :"-infinity"}

      # abs(x) < 1.0 is generally safe from epsilon issues for this check
      abs(x) < 1.0 ->
        denominator = 1.0 - x
        # This zero? check might be redundant if the equal?(x, 1.0) above catches it.
        if zero?(denominator, @boundary_epsilon) do
          # This implies x is extremely close to 1.0, should be caught by the first clause.
          # If it gets here, it's a very nuanced floating point state.
          {:error, :division_by_zero_in_atanh_log}
        else
          numerator = 1.0 + x
          # If numerator or denominator evaluate to non-positive due to x being near -1 or 1
          # Denominator strictly < 0 if x > 1
          if numerator <= 0.0 or denominator < 0.0 do
            {:error, :log_domain_error_in_atanh}
          else
            value = numerator / denominator
            # Check result of division for log's domain
            if value <= 0.0 do
              {:error, :log_domain_error_in_atanh_after_division}
            else
              {:ok, 0.5 * :math.log(value)}
            end
          end
        end

      # This means abs(x) > 1.0 (and not caught by the epsilon checks for +/- 1.0)
      true ->
        {:error, :atanh_domain_error_abs_x_gt_1}
    end
  end

  # --- Helper to compute geodetic latitude from conformal latitude ---
  def geodetic_latitude_from_conformal(sin_chi, eccentricity) do
    s_initial = sin_chi

    case do_geodetic_lat_iteration(sin_chi, eccentricity, s_initial, s_initial + 1.0e99, 30) do
      {:ok, final_s} ->
        clamped_final_s = max(-1.0, min(1.0, final_s))
        {:ok, :math.asin(clamped_final_s)}

      {:error, reason} ->
        {:error, {:geodetic_lat_iteration_failed, reason}}
    end
  end

  def do_geodetic_lat_iteration(_sin_chi, _eccentricity, _current_s, _s_old, 0) do
    {:error, :max_iterations_reached_in_geodetic_lat}
  end

  def do_geodetic_lat_iteration(sin_chi, eccentricity, current_s, s_old, iterations_left) do
    # Use a specific, small epsilon for convergence check of the iteration itself
    iteration_epsilon = 1.0e-12
    # Using helper
    if equal?(current_s, s_old, iteration_epsilon) do
      {:ok, current_s}
    else
      val_for_atanh = eccentricity * current_s

      p_arg =
        cond do
          val_for_atanh >= 1.0 -> 1.0 - @boundary_epsilon
          val_for_atanh <= -1.0 -> -1.0 + @boundary_epsilon
          true -> val_for_atanh
        end

      with {:ok, atanh_val_p_arg} <- atanh(p_arg),
           true <-
             is_number(atanh_val_p_arg) or
               {:error, {:atanh_returned_non_numeric_for_p_exp_iteration, atanh_val_p_arg}},
           {:ok, p_exp_arg} <- {:ok, eccentricity * atanh_val_p_arg},
           {:ok, p_val} <-
             {:ok,
              cond do
                p_exp_arg > 709.0 -> 1.0e300
                p_exp_arg < -709.0 -> 1.0e-300
                true -> :math.exp(p_exp_arg)
              end} do
        # p_val is now a number
        p_sq = p_val * p_val
        one_plus_sin_chi = 1.0 + sin_chi
        one_minus_sin_chi = 1.0 - sin_chi
        denominator = one_plus_sin_chi * p_sq + one_minus_sin_chi

        # Using helper
        if zero?(denominator, @boundary_epsilon) do
          {:error, :zero_denominator_in_geodetic_lat_iteration}
        else
          next_s_val_calc = (one_plus_sin_chi * p_sq - one_minus_sin_chi) / denominator

          if finite_float?(next_s_val_calc) do
            do_geodetic_lat_iteration(
              sin_chi,
              eccentricity,
              next_s_val_calc,
              current_s,
              iterations_left - 1
            )
          else
            {:error, {:invalid_next_s_in_iteration, next_s_val_calc}}
          end
        end
      else
        {:error, reason} -> {:error, {:iteration_step_failed_in_with, reason}}
      end
    end
  end

  # --- Helper to compute hyperbolic series terms ---
  def compute_hyperbolic_series(two_x, num_terms) do
    c2kx0 = :math.cosh(two_x)
    s2kx0 = :math.sinh(two_x)

    {c_list_rev, s_list_rev} =
      Enum.reduce(1..(num_terms - 1), {[c2kx0], [s2kx0]}, fn _k_idx, {acc_c, acc_s} ->
        prev_c = List.first(acc_c)
        prev_s = List.first(acc_s)

        next_c = prev_c * c2kx0 + prev_s * s2kx0
        next_s = prev_s * c2kx0 + prev_c * s2kx0
        {[next_c | acc_c], [next_s | acc_s]}
      end)

    {List.to_tuple(Enum.reverse(c_list_rev)), List.to_tuple(Enum.reverse(s_list_rev))}
  end

  # --- Helper to compute trigonometric series terms ---
  def compute_trig_series(two_y, num_terms) do
    c2ky0 = :math.cos(two_y)
    s2ky0 = :math.sin(two_y)

    {c_list_rev, s_list_rev} =
      Enum.reduce(1..(num_terms - 1), {[c2ky0], [s2ky0]}, fn _k_idx, {acc_c, acc_s} ->
        prev_c = List.first(acc_c)
        prev_s = List.first(acc_s)

        next_c = prev_c * c2ky0 - prev_s * s2ky0
        next_s = prev_s * c2ky0 + prev_c * s2ky0
        {[next_c | acc_c], [next_s | acc_s]}
      end)

    {List.to_tuple(Enum.reverse(c_list_rev)), List.to_tuple(Enum.reverse(s_list_rev))}
  end

  defp get_tm_coefficients(a, f) do
    if a == Constants.wgs84_a() and f == Constants.wgs84_f() do
      {:ok, {@wgs84_a_coeffs, @wgs84_b_coeffs, @wgs84_r4oa}}
    else
      # TMHelpers.generate_coefficients returns a raw tuple, wrap it in {:ok, ...}
      {:ok, TMHelpers.generate_coefficients(f)}
    end
  end

  defp nudge_value_if_at_boundary(value) do
    if abs(value) >= 1.0, do: value / (1.0 + 1.0e-15), else: value
  end

  @doc """
  Converts Geodetic coordinates (latitude, longitude) to Transverse Mercator
  easting and northing.
  Returns `{:ok, %MapProjectionCoordinates{}}` or `{:error, reason}`.
  """
  def convert_from_geodetic(
        %GeodeticCoordinates{} = geodetic_coords,
        projection_params,
        ellipsoid_params
      ) do
    a = Keyword.get(ellipsoid_params, :a, Constants.wgs84_a())
    f = Keyword.get(ellipsoid_params, :f, Constants.wgs84_f())
    es = Keyword.get(ellipsoid_params, :es, Constants.wgs84_es())

    origin_lon = Keyword.fetch!(projection_params, :central_meridian)
    false_easting = Keyword.fetch!(projection_params, :false_easting)
    false_northing = Keyword.fetch!(projection_params, :false_northing)
    k0 = Keyword.fetch!(projection_params, :scale_factor)

    with {:ok, {a_coeffs_tuple, _b_coeffs_tuple_ignored, r4oa_val}} <- get_tm_coefficients(a, f),
         {:ok, k0r4} <- {:ok, r4oa_val * k0 * a},
         {:ok, lambda} <-
           {:ok,
            (
              lambda_raw = geodetic_coords.longitude - origin_lon

              cond do
                lambda_raw > Constants.pi() -> lambda_raw - Constants.two_pi()
                lambda_raw < -Constants.pi() -> lambda_raw + Constants.two_pi()
                true -> lambda_raw
              end
            )},
         {:ok, sin_phi} <- {:ok, :math.sin(geodetic_coords.latitude)},
         {:ok, cos_phi} <- {:ok, :math.cos(geodetic_coords.latitude)},
         {:ok, val_for_atanh1} <- {:ok, es * sin_phi},
         {:ok, p_arg1} <- {:ok, nudge_value_if_at_boundary(val_for_atanh1)},
         {:ok, atanh_val1} <- atanh(p_arg1),
         true <-
           is_number(atanh_val1) or {:error, {:atanh_returned_non_numeric_for_p_exp, atanh_val1}},
         {:ok, p_exp_arg_fwd} <- {:ok, es * atanh_val1},
         {:ok, p_val} <-
           {:ok,
            cond do
              p_exp_arg_fwd > 709.0 -> 1.0e300
              p_exp_arg_fwd < -709.0 -> 1.0e-300
              true -> :math.exp(p_exp_arg_fwd)
            end},
         true <- not zero?(p_val, @boundary_epsilon) or {:error, :p_val_is_zero_from_exp},
         {:ok, part1} <- {:ok, (1.0 + sin_phi) / p_val},
         {:ok, part2} <- {:ok, (1.0 - sin_phi) * p_val},
         {:ok, denom_chi} <- {:ok, part1 + part2},
         true <- not zero?(denom_chi, @boundary_epsilon) or {:error, :denom_chi_is_zero},
         {:ok, cos_chi} <- {:ok, 2.0 * cos_phi / denom_chi},
         {:ok, sin_chi} <- {:ok, (part1 - part2) / denom_chi},
         {:ok, val_for_atanh2} <- {:ok, cos_chi * :math.sin(lambda)},
         {:ok, p_arg2} <- {:ok, nudge_value_if_at_boundary(val_for_atanh2)},
         {:ok, u_prime} <- atanh(p_arg2),
         true <-
           is_number(u_prime) or {:error, {:atanh_returned_non_numeric_for_u_prime, u_prime}} do
      warning_message =
        if abs(lambda) > @max_delta_long,
          do: "Longitude is far from Central Meridian, distortion may be significant",
          else: ""

      v_prime = :math.atan2(sin_chi, cos_chi * :math.cos(lambda))

      {c2ku_list_tup, s2ku_list_tup} = compute_hyperbolic_series(2.0 * u_prime, @n_terms)
      {c2kv_list_tup, s2kv_list_tup} = compute_trig_series(2.0 * v_prime, @n_terms)

      x_star =
        Enum.reduce(0..(@n_terms - 1), 0.0, fn k, acc ->
          acc + elem(a_coeffs_tuple, k) * elem(s2ku_list_tup, k) * elem(c2kv_list_tup, k)
        end) + u_prime

      y_star =
        Enum.reduce(0..(@n_terms - 1), 0.0, fn k, acc ->
          acc + elem(a_coeffs_tuple, k) * elem(c2ku_list_tup, k) * elem(s2kv_list_tup, k)
        end) + v_prime

      easting = k0r4 * x_star + false_easting
      northing = k0r4 * y_star + false_northing

      {:ok,
       %MapProjectionCoordinates{
         easting: easting,
         northing: northing,
         warning_message: warning_message
       }}
    else
      # This block handles errors from the `with` clause
      # The first clause that fails will result in its error tuple being passed here
      {:error, error_details} -> {:error, {:convert_from_geodetic_failed, error_details}}
    end
  end

  @doc """
  Converts Transverse Mercator easting and northing to Geodetic coordinates
  (latitude, longitude).
  Returns `{:ok, %GeodeticCoordinates{}}` or `{:error, reason}`.
  """
  def convert_to_geodetic(
        %MapProjectionCoordinates{} = proj_coords,
        projection_params,
        ellipsoid_params
      ) do
    a = Keyword.get(ellipsoid_params, :a, Constants.wgs84_a())
    f = Keyword.get(ellipsoid_params, :f, Constants.wgs84_f())
    es = Keyword.get(ellipsoid_params, :es, Constants.wgs84_es())

    origin_lon = Keyword.fetch!(projection_params, :central_meridian)
    false_easting = Keyword.fetch!(projection_params, :false_easting)
    false_northing = Keyword.fetch!(projection_params, :false_northing)
    k0 = Keyword.fetch!(projection_params, :scale_factor)

    case get_tm_coefficients(a, f) do
      {:ok, {_a_coeffs_tuple_ignored, b_coeffs_tuple, r4oa_val}} ->
        with {:ok, k0r4_product} <- {:ok, r4oa_val * k0 * a},
             true <-
               not zero?(k0r4_product, @boundary_epsilon) or {:error, :k0r4_product_is_zero},
             {:ok, k0r4_inv} <- {:ok, 1.0 / k0r4_product},
             {:ok, x_star} <- {:ok, (proj_coords.easting - false_easting) * k0r4_inv},
             {:ok, y_star} <- {:ok, (proj_coords.northing - false_northing) * k0r4_inv},
             {:ok, {c2kx_list_tup, s2kx_list_tup}} <-
               {:ok, compute_hyperbolic_series(2.0 * x_star, @n_terms)},
             {:ok, {c2ky_list_tup, s2ky_list_tup}} <-
               {:ok, compute_trig_series(2.0 * y_star, @n_terms)},
             {:ok, u_prime} <-
               {:ok,
                Enum.reduce(0..(@n_terms - 1), 0.0, fn k, acc ->
                  acc + elem(b_coeffs_tuple, k) * elem(s2kx_list_tup, k) * elem(c2ky_list_tup, k)
                end) + x_star},
             {:ok, v_prime} <-
               {:ok,
                Enum.reduce(0..(@n_terms - 1), 0.0, fn k, acc ->
                  acc + elem(b_coeffs_tuple, k) * elem(c2kx_list_tup, k) * elem(s2ky_list_tup, k)
                end) + y_star},
             {:ok, lambda_from_cm} <- {:ok, :math.atan2(:math.sinh(u_prime), :math.cos(v_prime))},
             {:ok, cosh_u_prime} <- {:ok, :math.cosh(u_prime)},
             true <-
               not zero?(cosh_u_prime, @boundary_epsilon) or
                 {:error, :cosh_u_prime_is_zero_for_sin_chi},
             {:ok, sin_chi_raw} <- {:ok, :math.sin(v_prime) / cosh_u_prime},
             {:ok, clamped_sin_chi} <- {:ok, max(-1.0, min(1.0, sin_chi_raw))},
             {:ok, latitude} <- geodetic_latitude_from_conformal(clamped_sin_chi, es) do
          longitude_raw = origin_lon + lambda_from_cm

          longitude =
            cond do
              longitude_raw > Constants.pi() -> longitude_raw - Constants.two_pi()
              longitude_raw < -Constants.pi() -> longitude_raw + Constants.two_pi()
              true -> longitude_raw
            end

          warning_message =
            if abs(lambda_from_cm) > @max_delta_long,
              do: "Longitude is far from Central Meridian, distortion may be significant",
              else: ""

          {:ok,
           %GeodeticCoordinates{
             longitude: longitude,
             latitude: latitude,
             height: 0.0,
             warning_message: warning_message
           }}
        else
          # This block handles errors from the inner `with` clause
          {:error, error_details_inner} ->
            {:error, {:convert_to_geodetic_failed_inner_with, error_details_inner}}
        end

        # This else matches the `case get_tm_coefficients(a,f)` if it could return an error
        # Since get_tm_coefficients currently always returns {:ok, ...} or TMHelpers raises,
        # this error case for get_tm_coefficients is not strictly reachable from its current form.
        # However, it's good practice if get_tm_coefficients were to be refactored to return errors.

        # {:error, reason_coeffs} ->
        #   {:error, {:get_tm_coefficients_failed, reason_coeffs}}
    end
  end
end
