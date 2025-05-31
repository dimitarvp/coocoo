defmodule CooCoo do
  @moduledoc """
  Main API for CooCoo library, providing high-level coordinate conversions.
  """

  alias CooCoo.Coordinates.GeodeticCoordinates
  alias CooCoo.Coordinates.MGRSCoordinates
  # The module doing the actual work
  alias CooCoo.MGRS
  # For degree/radian conversion
  alias CooCoo.Constants

  @doc """
  Converts GPS (latitude, longitude) coordinates to an MGRS string.

  Latitude and longitude are expected in decimal degrees.
  Height is in meters (optional, defaults to 0.0, not used in MGRS string itself).

  Options:
    - `:precision` (integer 0-5): Desired MGRS precision. Defaults to 5 (1-meter).
      - 0: 100km
      - 1: 10km
      - 2: 1km
      - 3: 100m
      - 4: 10m
      - 5: 1m

  Returns:
    - `{:ok, mgrs_string}` on success.
    - `{:ok, {mgrs_string, warning_message}}` on success with a warning.
    - `{:error, reason}` on failure (e.g., invalid input).
  """
  def gps_to_mgrs(latitude_deg, longitude_deg, opts \\ []) do
    precision = Keyword.get(opts, :precision, 5)
    height = Keyword.get(opts, :height, 0.0)

    try do
      geodetic_coords = %GeodeticCoordinates{
        latitude: latitude_deg * Constants.pi_over_180(),
        longitude: longitude_deg * Constants.pi_over_180(),
        height: height
      }

      %MGRSCoordinates{mgrs_string: mgrs_str, warning_message: warning} =
        MGRS.convert_from_geodetic(geodetic_coords, precision)

      if warning == "" do
        {:ok, mgrs_str}
      else
        {:ok, {mgrs_str, warning}}
      end
    rescue
      # Catch specific errors from underlying modules
      e in [ArgumentError, RuntimeError] ->
        {:error, Exception.message(e)}
        # Potentially other specific exceptions if MGRS.convert raises them
    end
  end

  @doc """
  Converts an MGRS string to GPS (latitude, longitude) coordinates.

  Returns:
    - `{:ok, %{latitude_deg: float, longitude_deg: float, height_m: float}}` on success.
    - `{:ok, {%{latitude_deg: float, longitude_deg: float, height_m: float}, warning_message}}` on success with a warning.
    - `{:error, reason}` on failure.

  The returned latitude and longitude are in decimal degrees. Height is in meters (will be 0.0 as MGRS is 2D).
  """
  def mgrs_to_gps(mgrs_string) when is_binary(mgrs_string) do
    try do
      # Provide a default precision; it will be correctly parsed from the string
      # by the underlying MGRS.break_mgrs_string function.
      # Default precision, will be overridden
      mgrs_coords_input = %MGRSCoordinates{mgrs_string: mgrs_string, precision: 0}

      %GeodeticCoordinates{
        latitude: lat_rad,
        longitude: lon_rad,
        height: h,
        warning_message: warning
      } =
        MGRS.convert_to_geodetic(mgrs_coords_input)

      output_map = %{
        latitude_deg: Constants.radians_to_degrees(lat_rad),
        longitude_deg: Constants.radians_to_degrees(lon_rad),
        height_m: h
      }

      if warning == "" do
        {:ok, output_map}
      else
        {:ok, {output_map, warning}}
      end
    rescue
      e in [ArgumentError, RuntimeError] ->
        {:error, Exception.message(e)}
    end
  end
end
