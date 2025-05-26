defmodule CooCoo.UPS do
  @moduledoc """
  Provides conversions between Geodetic and Universal Polar Stereographic (UPS) coordinates.
  UPS is a specific application of the Polar Stereographic projection.
  """

  alias CooCoo.Constants
  alias CooCoo.Coordinates.GeodeticCoordinates
  alias CooCoo.Coordinates.UPSCoordinates
  # To handle PS results
  alias CooCoo.Coordinates.MapProjectionCoordinates

  # UPS specific parameters
  @ups_false_easting 2_000_000.0
  @ups_false_northing 2_000_000.0
  # Central meridian for UPS is 0 degrees (by definition)
  @ups_origin_longitude 0.0
  # Scale factor at the pole for UPS
  @ups_scale_factor 0.994

  # UPS operational latitude limits
  # North Pole: 83.5 to 90 ( Technically, TM extends to 84, so UPS is > 84 for MGRS context)
  # South Pole: -79.5 to -90 (Technically, TM extends to -80, so UPS is < -80 for MGRS context)
  # Using the MGRS context boundaries for applicability here.
  # Approx 84 deg N
  @min_north_latitude_rad Constants.max_mgrs_non_polar_lat()
  # Approx 80 deg S
  @max_south_latitude_rad Constants.min_mgrs_non_polar_lat()

  @doc """
  Converts Geodetic coordinates to UPS coordinates.
  """
  def convert_from_geodetic(%GeodeticCoordinates{} = geodetic_coords) do
    latitude = geodetic_coords.latitude
    # longitude = geodetic_coords.longitude # UPS central meridian is fixed at 0

    hemisphere =
      cond do
        # Includes North Pole exactly at 0.0 lat for this check
        latitude >= 0.0 ->
          if latitude < @min_north_latitude_rad do
            throw(%ArgumentError{
              message: "Latitude not in Northern UPS zone (must be >= ~84 degrees N)"
            })
          end

          ?N

        latitude < 0.0 ->
          if latitude > @max_south_latitude_rad do
            throw(%ArgumentError{
              message: "Latitude not in Southern UPS zone (must be <= ~80 degrees S)"
            })
          end

          ?S
      end

    # For Polar Stereographic, the 'standard_parallel' or 'scale_factor' + 'hemisphere'
    # determines the projection math. UPS is defined with a scale factor at the pole.
    projection_params = [
      # UPS central meridian is 0.0
      central_meridian: @ups_origin_longitude,
      scale_factor: @ups_scale_factor,
      # Pass hemisphere to PS to guide its internal math adjustments
      hemisphere: hemisphere,
      false_easting: @ups_false_easting,
      false_northing: @ups_false_northing
    ]

    ellipsoid_params = [a: Constants.wgs84_a(), es: Constants.wgs84_es()]

    %MapProjectionCoordinates{
      easting: ps_easting,
      northing: ps_northing,
      warning_message: ps_warning
    } =
      CooCoo.Projections.PolarStereographic.convert_from_geodetic(
        geodetic_coords,
        projection_params,
        ellipsoid_params
      )

    # UPS coordinates are defined within a certain range from the pole center
    # Min/Max East/North for UPS: 0 to 4,000,000 (MGRS.cpp defines)
    # This validation can also be done here if needed, or by the caller (e.g. MGRS)

    %UPSCoordinates{
      hemisphere: hemisphere,
      easting: ps_easting,
      northing: ps_northing,
      warning_message: ps_warning
    }
  end

  @doc """
  Converts UPS coordinates to Geodetic coordinates.
  """
  def convert_to_geodetic(%UPSCoordinates{} = ups_coords) do
    hemisphere = ups_coords.hemisphere
    easting = ups_coords.easting
    northing = ups_coords.northing

    if hemisphere != ?N and hemisphere != ?S do
      throw(%ArgumentError{message: "Invalid UPS hemisphere"})
    end

    # UPS coordinates range check (as per MGRS.cpp for UPS derived values)
    if easting < Constants.ups_min_east_north() or easting > Constants.ups_max_east_north_mgrs() or
         northing < Constants.ups_min_east_north() or
         northing > Constants.ups_max_east_north_mgrs() do
      # MGRS.cpp throws an error for this if converting MGRS polar to UPS.
      # Here, we can return a warning or error. For now, proceed.
      # A warning might be appropriate.
    end

    projection_params = [
      central_meridian: @ups_origin_longitude,
      scale_factor: @ups_scale_factor,
      hemisphere: hemisphere,
      false_easting: @ups_false_easting,
      false_northing: @ups_false_northing
    ]

    ellipsoid_params = [a: Constants.wgs84_a(), es: Constants.wgs84_es()]

    map_proj_coords = %MapProjectionCoordinates{easting: easting, northing: northing}

    %GeodeticCoordinates{longitude: lon, latitude: lat, warning_message: ps_warning} =
      CooCoo.Projections.PolarStereographic.convert_to_geodetic(
        map_proj_coords,
        projection_params,
        ellipsoid_params
      )

    # Validate resulting latitude is in the correct polar region
    warning_message =
      cond do
        hemisphere == ?N and lat < @min_north_latitude_rad ->
          if(ps_warning != "", do: ps_warning <> "\n", else: "") <>
            "Resulting latitude outside Northern UPS zone"

        hemisphere == ?S and lat > @max_south_latitude_rad ->
          if(ps_warning != "", do: ps_warning <> "\n", else: "") <>
            "Resulting latitude outside Southern UPS zone"

        true ->
          ps_warning
      end

    %GeodeticCoordinates{
      longitude: lon,
      latitude: lat,
      height: 0.0,
      warning_message: warning_message
    }
  end
end
