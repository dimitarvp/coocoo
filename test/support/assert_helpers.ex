defmodule CooCoo.Test.AssertHelpers do
  @moduledoc """
  Custom assertion helpers for CooCoo tests.
  """

  import ExUnit.Assertions
  # Alias CooCoo.Constants
  alias CooCoo.Constants

  @default_float_delta 1.0e-9
  @default_meter_delta 0.001
  # @default_degree_delta 1.0e-7 # This was unused, removing for now.
  # We can re-add if a specific helper needs it.

  @doc """
  Asserts that two floating-point numbers are close within a given delta.
  Includes a custom message.
  """
  def assert_floats_close(actual, expected, delta \\ @default_float_delta, message \\ nil) do
    message = message || "Expected #{expected} to be close to #{actual} (delta: #{delta})"
    assert_in_delta actual, expected, delta, message
  end

  @doc """
  Asserts that two GeodeticCoordinates structs are approximately equal.
  Compares latitude and longitude in radians, and optionally height in meters.
  """
  def assert_geodetic_coords_approx_equal(
        %CooCoo.Coordinates.GeodeticCoordinates{} = actual,
        %CooCoo.Coordinates.GeodeticCoordinates{} = expected,
        opts \\ []
      ) do
    # Renamed for clarity
    lat_lon_delta_rad = Keyword.get(opts, :lat_lon_delta_rad, @default_float_delta)
    # Renamed for clarity
    height_delta_m = Keyword.get(opts, :height_delta_m, @default_meter_delta)
    check_height = Keyword.get(opts, :check_height, true)

    assert_in_delta actual.latitude,
                    expected.latitude,
                    lat_lon_delta_rad,
                    "Latitude mismatch. Expected: #{expected.latitude}, Got: #{actual.latitude}"

    assert_in_delta actual.longitude,
                    expected.longitude,
                    lat_lon_delta_rad,
                    "Longitude mismatch. Expected: #{expected.longitude}, Got: #{actual.longitude}"

    if check_height do
      assert_in_delta actual.height,
                      expected.height,
                      height_delta_m,
                      "Height mismatch. Expected: #{expected.height}, Got: #{actual.height}"
    end
  end

  @doc """
  Asserts that two MapProjectionCoordinates (or similar like UTM/UPS) structs
  are approximately equal. Compares easting and northing in meters.
  """
  def assert_projected_coords_approx_equal(
        actual,
        expected,
        delta_meters \\ @default_meter_delta
      ) do
    # Check if actual and expected have :easting and :northing fields
    # This is a bit of dynamic checking, but makes the helper more flexible.
    unless Map.has_key?(actual, :easting) and Map.has_key?(actual, :northing) and
             Map.has_key?(expected, :easting) and Map.has_key?(expected, :northing) do
      flunk(
        "Actual or expected struct does not have :easting and :northing fields for comparison. Actual: #{inspect(actual)}, Expected: #{inspect(expected)}"
      )
    end

    assert_in_delta actual.easting,
                    expected.easting,
                    delta_meters,
                    "Easting mismatch. Expected: #{expected.easting}, Got: #{actual.easting} for #{inspect(actual)}"

    assert_in_delta actual.northing,
                    expected.northing,
                    delta_meters,
                    "Northing mismatch. Expected: #{expected.northing}, Got: #{actual.northing} for #{inspect(actual)}"
  end

  @doc """
  Asserts that two MGRSCoordinate structs are equal.
  Compares the MGRS string and precision.
  """
  def assert_mgrs_coords_equal(
        %CooCoo.Coordinates.MGRSCoordinates{} = actual,
        %CooCoo.Coordinates.MGRSCoordinates{} = expected
      ) do
    assert actual.mgrs_string == expected.mgrs_string,
           "MGRS string mismatch. Expected: \"#{expected.mgrs_string}\", Got: \"#{actual.mgrs_string}\""

    assert actual.precision == expected.precision,
           "MGRS precision mismatch. Expected: #{expected.precision}, Got: #{actual.precision}"
  end

  @mgrs_precision_deltas %{
    5 => 0.5,
    4 => 5.0,
    3 => 50.0,
    2 => 500.0,
    1 => 5000.0,
    0 => 50000.0
  }

  @doc """
  Returns an appropriate delta in meters for comparing projected coordinates
  derived from MGRS strings of a given precision.
  """
  def mgrs_precision_to_meter_delta(precision_level) when precision_level in 0..5 do
    Map.get(@mgrs_precision_deltas, precision_level, @default_meter_delta)
  end

  @doc """
  Returns an appropriate delta in radians for comparing geodetic coordinates
  derived from MGRS strings of a given precision.
  """
  def mgrs_precision_to_radian_delta(precision_level) when precision_level in 0..5 do
    meter_delta = mgrs_precision_to_meter_delta(precision_level)
    degree_delta_approx = meter_delta / 111_000.0
    # Correctly aliased call
    degree_delta_approx * Constants.pi_over_180()
  end
end
