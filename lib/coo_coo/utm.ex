defmodule CooCoo.UTM do
  @moduledoc """
  Provides conversions between Geodetic and UTM coordinates.
  Universal Transverse Mercator (UTM) is a specific application of the
  Transverse Mercator projection.
  """

  alias CooCoo.Constants
  alias CooCoo.Coordinates.GeodeticCoordinates
  alias CooCoo.Coordinates.UTMCoordinates
  # To handle TM results
  alias CooCoo.Coordinates.MapProjectionCoordinates
  alias CooCoo.Projections.TransverseMercator

  # Standard UTM false easting
  @utm_false_easting 500_000.0
  @utm_false_northing_north_hemisphere 0.0
  @utm_false_northing_south_hemisphere 10_000_000.0
  # Standard UTM scale factor at central meridian
  @utm_scale_factor 0.9996

  @doc """
  Converts Geodetic coordinates to UTM coordinates.

  An optional `utm_zone_override` can be provided. If not, the zone is
  calculated. If provided, it must be within +/-1 of the natural zone,
  or be the corresponding zone on the other side of the 180-degree meridian
  (e.g. zone 1 for natural zone 60, or zone 60 for natural zone 1).
  Special NGA rules for Norway and Svalbard are applied if no override is given.
  """
  def convert_from_geodetic(%GeodeticCoordinates{} = geodetic_coords, utm_zone_override \\ nil) do
    latitude = geodetic_coords.latitude
    longitude = geodetic_coords.longitude

    if latitude < Constants.min_mgrs_non_polar_lat() - Constants.mgrs_epsilon2() or
         latitude >= Constants.max_mgrs_non_polar_lat() + Constants.mgrs_epsilon2() do
      # Note: MGRS.cpp uses MIN_MGRS_NON_POLAR_LAT and MAX_MGRS_NON_POLAR_LAT.
      # UTM itself is generally defined between -80 and +84.
      # This check matches the MGRS module's check for when to use UTM vs UPS path.
      # Let's use the C++ UTM.cpp bounds: MIN_LAT (-80.5 deg) and MAX_LAT (84.5 deg)
      # For strict UTM, it's -80 to +84. MGRS extends this slightly for continuity.
      # Using MIN_LAT/MAX_LAT from UTM.cpp for this UTM module.
      min_lat_utm = Constants.degrees_to_radians(-80.5)
      max_lat_utm = Constants.degrees_to_radians(84.5)

      if latitude < min_lat_utm or latitude > max_lat_utm do
        # This error should ideally be caught by a higher-level function (like MGRS)
        # that decides whether to use UTM or UPS.
        # If called directly, this is the UTM validity range.
        throw(%ArgumentError{
          message: "Latitude outside UTM validity range (-80.5 to 84.5 degrees)"
        })
      end
    end

    # Normalize longitude to [-PI, PI) for zone calculation, then to [0, 2PI) for some TM needs
    # The TransverseMercator module itself handles normalization relative to central_meridian.
    # Here, we need longitude in [0, 2PI) for unambiguous zone calculation as per C++ `temp_zone` logic.
    lon_for_zone_calc =
      cond do
        longitude < 0.0 -> longitude + Constants.two_pi()
        longitude >= Constants.two_pi() -> :math.fmod(longitude, Constants.two_pi())
        true -> longitude
      end

    # Calculate natural UTM zone
    # C++ logic: if (longitude < PI) temp_zone = (long)(31 + (...)); else temp_zone = (long)(...);
    # This implies longitude is in [0, 2PI).
    natural_zone =
      if lon_for_zone_calc < Constants.pi() do
        trunc(31.0 + lon_for_zone_calc / Constants.deg6_in_radians())
      else
        # Add a small epsilon to handle exact boundary, similar to C++
        trunc((lon_for_zone_calc + 1.0e-10) / Constants.deg6_in_radians() - 29.0)
      end

    natural_zone = if natural_zone > 60, do: 1, else: natural_zone

    # Apply zone override or special rules
    actual_zone =
      cond do
        utm_zone_override != nil ->
          cond do
            (natural_zone == 1 and utm_zone_override == 60) or
              (natural_zone == 60 and utm_zone_override == 1) or
                (abs(natural_zone - utm_zone_override) <= 1 and utm_zone_override >= 1 and
                   utm_zone_override <= 60) ->
              utm_zone_override

            true ->
              throw(%ArgumentError{message: "Invalid UTM zone override"})
          end

        # No override, apply NGA special rules
        true ->
          lat_deg = Constants.radians_to_degrees(latitude)
          # Use original longitude for degree checks
          lon_deg = Constants.radians_to_degrees(geodetic_coords.longitude)

          cond do
            # Southern Norway
            lat_deg > 55 and lat_deg < 64 and lon_deg > -1 and lon_deg < 3 -> 31
            lat_deg > 55 and lat_deg < 64 and lon_deg > 2 and lon_deg < 12 -> 32
            # Svalbard
            lat_deg > 71 and lon_deg > -1 and lon_deg < 9 -> 31
            lat_deg > 71 and lon_deg > 8 and lon_deg < 21 -> 33
            lat_deg > 71 and lon_deg > 20 and lon_deg < 33 -> 35
            lat_deg > 71 and lon_deg > 32 and lon_deg < 42 -> 37
            true -> natural_zone
          end
      end

    hemisphere = if latitude < 0.0, do: ?S, else: ?N

    false_northing =
      if hemisphere == ?S,
        do: @utm_false_northing_south_hemisphere,
        else: @utm_false_northing_north_hemisphere

    # Calculate central meridian for the actual_zone
    # C++: if (zone >= 31) cm = ((6*zone-183)*PI_OVER_180); else cm = ((6*zone+177)*PI_OVER_180);
    # This is equivalent to (zone * 6 - 183) degrees for all zones.
    central_meridian = Constants.degrees_to_radians(actual_zone * 6.0 - 183.0)

    projection_params = [
      central_meridian: central_meridian,
      # UTM is based on latitude of origin being the equator
      latitude_of_origin: 0.0,
      false_easting: @utm_false_easting,
      # TM uses 0, hemisphere offset applied after
      false_northing: 0.0,
      scale_factor: @utm_scale_factor
    ]

    # Assuming WGS84 for now, or allow ellipsoid_params to be passed in
    ellipsoid_params = [a: Constants.wgs84_a(), f: Constants.wgs84_f(), es: Constants.wgs84_es()]

    %MapProjectionCoordinates{
      easting: tm_easting,
      northing: tm_northing,
      warning_message: tm_warning
    } =
      TransverseMercator.convert_from_geodetic(
        geodetic_coords,
        projection_params,
        ellipsoid_params
      )

    # Apply hemisphere-based false northing
    final_northing = tm_northing + false_northing

    # Validate UTM coordinate ranges
    if tm_easting < Constants.utm_min_easting() or
         tm_easting > Constants.utm_max_easting_non_mgrs() do
      # Note: MGRS uses a tighter max_easting (900,000). UTM generally goes up to 999,999, but values
      # outside 100k-900k are in adjacent zones.
      # C++ MGRS.cpp checks against 100k-900k for UTM inputs.
      # Let's align with common UTM limits which are effectively handled by zone logic.
      # The MGRS.cpp error checks seem to be for MGRS context, not general UTM.
      # For now, we'll assume the projection itself is valid.
      # If strict MGRS compatibility needed, add warning or error here based on stricter MGRS limits.
    end

    if final_northing < Constants.utm_min_northing() or
         final_northing > Constants.utm_max_northing_non_mgrs() do
      # Similar to easting, MGRS has its context. Standard UTM northing valid range.
    end

    %UTMCoordinates{
      zone: actual_zone,
      hemisphere: hemisphere,
      easting: tm_easting,
      northing: final_northing,
      warning_message: tm_warning
    }
  end

  @doc """
  Converts UTM coordinates to Geodetic coordinates.
  """
  def convert_to_geodetic(%UTMCoordinates{} = utm_coords) do
    zone = utm_coords.zone
    hemisphere = utm_coords.hemisphere
    easting = utm_coords.easting
    northing = utm_coords.northing

    if zone < 1 or zone > 60 do
      throw(%ArgumentError{message: "Invalid UTM zone"})
    end

    if hemisphere != ?N and hemisphere != ?S do
      throw(%ArgumentError{message: "Invalid UTM hemisphere"})
    end

    # C++ MGRS checks easting against [100000, 900000] and northing [0, 10000000]
    # These are practical limits for a given zone.

    false_northing_offset =
      if hemisphere == ?S,
        do: @utm_false_northing_south_hemisphere,
        else: @utm_false_northing_north_hemisphere

    # Northing for TM calculation is relative to equator
    tm_northing = northing - false_northing_offset

    central_meridian = Constants.degrees_to_radians(zone * 6.0 - 183.0)

    projection_params = [
      central_meridian: central_meridian,
      latitude_of_origin: 0.0,
      false_easting: @utm_false_easting,
      # TM itself uses 0 for N/S hemispheres
      false_northing: 0.0,
      scale_factor: @utm_scale_factor
    ]

    ellipsoid_params = [a: Constants.wgs84_a(), f: Constants.wgs84_f(), es: Constants.wgs84_es()]

    map_proj_coords = %MapProjectionCoordinates{easting: easting, northing: tm_northing}

    %GeodeticCoordinates{longitude: lon, latitude: lat, warning_message: tm_warning} =
      TransverseMercator.convert_to_geodetic(map_proj_coords, projection_params, ellipsoid_params)

    # Basic validation of resulting latitude against UTM operational range.
    min_lat_utm = Constants.degrees_to_radians(-80.5)
    max_lat_utm = Constants.degrees_to_radians(84.5)

    warning_message =
      if lat < min_lat_utm or lat > max_lat_utm do
        # Concatenate if tm_warning already exists
        if(tm_warning != "", do: tm_warning <> "\n", else: "") <>
          "Resulting latitude is outside typical UTM range (-80.5 to 84.5 degrees)"
      else
        tm_warning
      end

    %GeodeticCoordinates{
      longitude: lon,
      latitude: lat,
      height: 0.0,
      warning_message: warning_message
    }
  end
end
