defmodule CooCoo.Projections.PolarStereographic do
  @moduledoc """
  Provides Polar Stereographic projection conversions.
  This module implements the core mathematical formulas for the projection.
  UPS is a specific application of Polar Stereographic.
  """

  alias CooCoo.Constants
  alias CooCoo.Coordinates.GeodeticCoordinates
  alias CooCoo.Coordinates.MapProjectionCoordinates

  def polar_pow(eccentricity, sin_latitude) do
    es_sin = eccentricity * sin_latitude
    denominator = 1.0 + es_sin

    if denominator == 0.0 do
      # This is an edge case, should not occur with valid ellipsoid and lat.
      # If (1-es_sin) is also 0, then it's 0^positive_power = 0.
      # If (1-es_sin) is non-zero, it's effectively division by zero in the base.
      # For now, let's return something that might indicate an issue or a large value.
      # A robust solution might raise an error or return a specific atom.
      # If es_sin = -1, then base is (1 - (-1)) / (1 + (-1)) = 2 / 0 -> problematic.
      # Given es < 1 and |sin_latitude| <= 1, es_sin will be < 1.
      # So 1+es_sin will not be 0 unless es = -1 (not possible).
      # And 1-es_sin will not be 0 unless es = 1 (not possible).
      # The base (1.0 - es_sin) / (1.0 + es_sin) should always be well-defined and positive.
      # Fallback, but should be reviewed if triggered.
      0.0
    else
      base = (1.0 - es_sin) / denominator
      exponent = eccentricity / 2.0
      :math.pow(base, exponent)
    end
  end

  @doc """
  Converts Geodetic coordinates to Polar Stereographic coordinates.
  """
  def convert_from_geodetic(
        %GeodeticCoordinates{} = geodetic_coords,
        projection_params,
        ellipsoid_params
      ) do
    a = Keyword.fetch!(ellipsoid_params, :a)
    es = Keyword.fetch!(ellipsoid_params, :es)

    central_meridian = Keyword.fetch!(projection_params, :central_meridian)
    false_easting = Keyword.fetch!(projection_params, :false_easting)
    false_northing = Keyword.fetch!(projection_params, :false_northing)

    input_latitude = geodetic_coords.latitude
    input_longitude = geodetic_coords.longitude

    # Determine if standard parallel or scale factor mode
    standard_parallel_param = Keyword.get(projection_params, :standard_parallel)
    scale_factor_param = Keyword.get(projection_params, :scale_factor)

    # Determine operational hemisphere
    is_southern_hemisphere =
      if hem = Keyword.get(projection_params, :hemisphere) do
        hem == ?S
      else
        # Infer from standard_parallel if provided, otherwise from input_latitude
        if standard_parallel_param, do: standard_parallel_param < 0.0, else: input_latitude < 0.0
      end

    # Adjust input lat/lon and central_meridian for internal N-Pole centric math
    latitude_eff = if is_southern_hemisphere, do: -input_latitude, else: input_latitude
    longitude_eff = if is_southern_hemisphere, do: -input_longitude, else: input_longitude
    op_central_meridian = if is_southern_hemisphere, do: -central_meridian, else: central_meridian

    # At the pole itself
    if abs(abs(latitude_eff) - Constants.pi_over_2()) < 1.0e-10 do
      %MapProjectionCoordinates{
        easting: false_easting,
        northing: false_northing,
        warning_message: ""
      }
    else
      dlam = longitude_eff - op_central_meridian

      dlam =
        cond do
          dlam > Constants.pi() -> dlam - Constants.two_pi()
          dlam < -Constants.pi() -> dlam + Constants.two_pi()
          true -> dlam
        end

      sin_lat_eff = :math.sin(latitude_eff)
      # essin is used in polar_pow
      # essin = es * sin_lat_eff
      pow_es_val = polar_pow(es, sin_lat_eff)
      t = :math.tan(Constants.pi_over_4() - latitude_eff / 2.0) / pow_es_val

      rho =
        cond do
          # Standard Parallel mode
          sp = standard_parallel_param ->
            # Effective standard parallel for calculation (always positive for internal math)
            eff_sp = abs(sp)

            if abs(abs(eff_sp) - Constants.pi_over_2()) <= 1.0e-10 do
              # Standard parallel is at the pole, implies k0=1 at pole for this formula path
              k90 = :math.sqrt(:math.pow(1.0 + es, 1.0 + es) * :math.pow(1.0 - es, 1.0 - es))
              2.0 * a * t / k90
            else
              # Standard parallel not at the pole
              sinolat_eff_sp = :math.sin(eff_sp)
              essin_eff_sp = es * sinolat_eff_sp
              pow_es_eff_sp = polar_pow(es, sinolat_eff_sp)
              cosolat_eff_sp = :math.cos(eff_sp)

              mc_val = cosolat_eff_sp / :math.sqrt(1.0 - essin_eff_sp * essin_eff_sp)
              a_mc_val = a * mc_val
              tc_val = :math.tan(Constants.pi_over_4() - eff_sp / 2.0) / pow_es_eff_sp

              if tc_val == 0.0, do: 0.0, else: a_mc_val * t / tc_val
            end

          # Scale Factor mode (scale_factor_param is k0 at the pole)
          sf = scale_factor_param ->
            k90 = :math.sqrt(:math.pow(1.0 + es, 1.0 + es) * :math.pow(1.0 - es, 1.0 - es))
            2.0 * a * sf * t / k90

          # Fallback if parameters are missing (should be caught by earlier validation)
          true ->
            throw(
              "PolarStereographic: Missing standard_parallel or scale_factor in projection_params"
            )
        end

      temp_easting = rho * :math.sin(dlam)
      # For North pole base projection
      temp_northing = -rho * :math.cos(dlam)

      {easting, northing} =
        if is_southern_hemisphere do
          {-temp_easting + false_easting, -temp_northing + false_northing}
        else
          {temp_easting + false_easting, temp_northing + false_northing}
        end

      %MapProjectionCoordinates{easting: easting, northing: northing, warning_message: ""}
    end
  end

  @doc """
  Converts Polar Stereographic coordinates to Geodetic coordinates.
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

    standard_parallel_param = Keyword.get(projection_params, :standard_parallel)
    scale_factor_param = Keyword.get(projection_params, :scale_factor)

    is_southern_hemisphere =
      if hem = Keyword.get(projection_params, :hemisphere) do
        hem == ?S
      else
        # Heuristic if not specified
        if standard_parallel_param,
          do: standard_parallel_param < 0.0,
          else: proj_coords.northing < false_northing
      end

    op_central_meridian = if is_southern_hemisphere, do: -central_meridian, else: central_meridian

    dx = proj_coords.easting - false_easting
    dy = proj_coords.northing - false_northing

    {dx_calc, dy_calc} =
      if is_southern_hemisphere, do: {-dx, -dy}, else: {dx, dy}

    rho = :math.sqrt(dx_calc * dx_calc + dy_calc * dy_calc)

    latitude_at_pole =
      if is_southern_hemisphere, do: -Constants.pi_over_2(), else: Constants.pi_over_2()

    if abs(rho) < 1.0e-10 do
      %GeodeticCoordinates{
        latitude: latitude_at_pole,
        longitude: central_meridian,
        warning_message: ""
      }
    else
      t_val =
        cond do
          sp = standard_parallel_param ->
            eff_sp = abs(sp)

            if abs(abs(eff_sp) - Constants.pi_over_2()) <= 1.0e-10 do
              k90 = :math.sqrt(:math.pow(1.0 + es, 1.0 + es) * :math.pow(1.0 - es, 1.0 - es))
              rho * k90 / (2.0 * a)
            else
              sinolat_eff_sp = :math.sin(eff_sp)
              essin_eff_sp = es * sinolat_eff_sp
              pow_es_eff_sp = polar_pow(es, sinolat_eff_sp)
              cosolat_eff_sp = :math.cos(eff_sp)
              mc_val = cosolat_eff_sp / :math.sqrt(1.0 - essin_eff_sp * essin_eff_sp)
              a_mc_val = a * mc_val
              tc_val = :math.tan(Constants.pi_over_4() - eff_sp / 2.0) / pow_es_eff_sp
              rho * tc_val / a_mc_val
            end

          sf = scale_factor_param ->
            k90 = :math.sqrt(:math.pow(1.0 + es, 1.0 + es) * :math.pow(1.0 - es, 1.0 - es))
            rho * k90 / (2.0 * a * sf)

          true ->
            throw("PolarStereographic Inverse: Missing standard_parallel or scale_factor")
        end

      phi_initial = Constants.pi_over_2() - 2.0 * :math.atan(t_val)

      latitude_eff =
        Enum.reduce(1..5, phi_initial, fn _iter, current_phi ->
          # essin_iter = es * :math.sin(current_phi) # This essin was unused in C++ loop as well
          pow_es_val = polar_pow(es, :math.sin(current_phi))
          Constants.pi_over_2() - 2.0 * :math.atan(t_val * pow_es_val)
        end)

      longitude_eff = op_central_meridian + :math.atan2(dx_calc, -dy_calc)

      latitude = if is_southern_hemisphere, do: -latitude_eff, else: latitude_eff
      longitude = if is_southern_hemisphere, do: -longitude_eff, else: longitude_eff

      longitude =
        cond do
          longitude > Constants.pi() -> longitude - Constants.two_pi()
          longitude < -Constants.pi() -> longitude + Constants.two_pi()
          true -> longitude
        end

      %GeodeticCoordinates{
        longitude: longitude,
        latitude: latitude,
        height: 0.0,
        warning_message: ""
      }
    end
  end
end
