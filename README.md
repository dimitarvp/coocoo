# CooCoo - Elixir GPS to MGRS Conversion Library

CooCoo is a pure Elixir library for converting between Geodetic (latitude/longitude)
coordinates and Military Grid Reference System (MGRS) coordinates. It aims to be
a self-contained Elixir implementation without external dependencies or NIFs.

This library is based on the logic found in the C++ GEOTRANS library.

## Installation

This library is not yet published as a Hex package. To use it, you would typically
add it as a dependency in your `mix.exs` file, pointing to its local path or a
Git repository if you host it.

For example, if it's in a local directory:

```elixir
def deps do
  [
    {:coocoo, path: "../coocoo"} # Adjust path as needed
  ]
end
```

## Usage

The primary interface for conversions is through the `CooCoo.MGRS` module.

### 1. GPS (Geodetic) to MGRS

To convert GPS coordinates (latitude and longitude in degrees, height in meters)
to an MGRS string, you first create a `CooCoo.Coordinates.GeodeticCoordinates`
struct and then pass it to `CooCoo.MGRS.convert_from_geodetic/2`.

The `convert_from_geodetic/2` function takes the geodetic coordinates and an
optional precision (0-5, where 5 is 1-meter precision, which is the default).

**Example 1: GPS to MGRS (1-meter precision)**

```elixir
# Assuming your CooCoo library is available
alias CooCoo.Constants
alias CooCoo.Coordinates.GeodeticCoordinates
alias CooCoo.Coordinates.MGRSCoordinates
alias CooCoo.MGRS

# Input: Latitude 38.8895 N, Longitude 77.0353 W (Washington Monument)
# Note: Longitude for MGRS/UTM calculations is often expected as positive East or negative West.
# The library handles internal normalization.
lat_deg = 38.8895
lon_deg = -77.0353 # West longitude is negative
height_m = 0.0 # Height is not directly used in MGRS string but part of Geodetic struct

# Convert degrees to radians for GeodeticCoordinates struct
lat_rad = lat_deg * Constants.pi_over_180()
lon_rad = lon_deg * Constants.pi_over_180()

geodetic_coords = %GeodeticCoordinates{
  latitude: lat_rad,
  longitude: lon_rad,
  height: height_m
}

# Convert with default 1-meter precision (precision = 5)
%MGRSCoordinates{mgrs_string: mgrs_str_p5} = MGRS.convert_from_geodetic(geodetic_coords)
IO.puts("GPS (#{lat_deg}N, #{lon_deg}W) to MGRS (1m): #{mgrs_str_p5}")

# Convert with 1km precision (precision = 2)
%MGRSCoordinates{mgrs_string: mgrs_str_p2} = MGRS.convert_from_geodetic(geodetic_coords, 2)
IO.puts("GPS (#{lat_deg}N, #{lon_deg}W) to MGRS (1km): #{mgrs_str_p2}")

# Convert with 10km precision (precision = 1)
%MGRSCoordinates{mgrs_string: mgrs_str_p1} = MGRS.convert_from_geodetic(geodetic_coords, 1)
IO.puts("GPS (#{lat_deg}N, #{lon_deg}W) to MGRS (10km): #{mgrs_str_p1}")

# Convert with 100km precision (precision = 0)
%MGRSCoordinates{mgrs_string: mgrs_str_p0} = MGRS.convert_from_geodetic(geodetic_coords, 0)
IO.puts("GPS (#{lat_deg}N, #{lon_deg}W) to MGRS (100km): #{mgrs_str_p0}")
```

**Expected Output for Example 1 (approximate, actual may vary slightly based on exact TM formulas vs other tools):**

```
GPS (38.8895N, -77.0353W) to MGRS (1m): 18SUJ2337105498
GPS (38.8895N, -77.0353W) to MGRS (1km): 18SUJ2305
GPS (38.8895N, -77.0353W) to MGRS (10km): 18SUJ20
GPS (38.8895N, -77.0353W) to MGRS (100km): 18SUJ
```

_(Note: The exact numeric part of the MGRS string depends on the precise TM implementation and rounding. The example values here are illustrative based on typical outputs for that location.)_

### 2. MGRS to GPS (Geodetic)

To convert an MGRS string back to Geodetic coordinates (latitude and longitude in degrees),
you create an `CooCoo.Coordinates.MGRSCoordinates` struct and pass it to
`CooCoo.MGRS.convert_to_geodetic/1`. The function returns a
`CooCoo.Coordinates.GeodeticCoordinates` struct with latitude and longitude in radians.

**Example 2: MGRS to GPS**

```elixir
# Assuming your CooCoo library is available
alias CooCoo.Constants
alias CooCoo.Coordinates.GeodeticCoordinates
alias CooCoo.Coordinates.MGRSCoordinates
alias CooCoo.MGRS

mgrs_string_1m = "18SUJ2337105498" # From previous example
mgrs_coords_1m = %MGRSCoordinates{mgrs_string: mgrs_string_1m}

%GeodeticCoordinates{latitude: lat_rad_1m, longitude: lon_rad_1m, warning_message: warning_1m} =
  MGRS.convert_to_geodetic(mgrs_coords_1m)

lat_deg_1m = Constants.radians_to_degrees(lat_rad_1m)
lon_deg_1m = Constants.radians_to_degrees(lon_rad_1m)

IO.puts("MGRS #{mgrs_string_1m} to GPS: Lat #{:erlang.float_to_binary(lat_deg_1m, decimals: 5)}, Lon #{:erlang.float_to_binary(lon_deg_1m, decimals: 5)}")
if warning_1m != "", do: IO.puts("Warning: #{warning_1m}")


mgrs_string_100km = "18SUJ" # 100km precision
# For MGRS strings with no numeric part, break_mgrs_string assigns precision 0
mgrs_coords_100km = %MGRSCoordinates{mgrs_string: mgrs_string_100km}

%GeodeticCoordinates{latitude: lat_rad_100km, longitude: lon_rad_100km, warning_message: warning_100km} =
  MGRS.convert_to_geodetic(mgrs_coords_100km)

lat_deg_100km = Constants.radians_to_degrees(lat_rad_100km)
lon_deg_100km = Constants.radians_to_degrees(lon_rad_100km)

IO.puts("MGRS #{mgrs_string_100km} to GPS: Lat #{:erlang.float_to_binary(lat_deg_100km, decimals: 2)}, Lon #{:erlang.float_to_binary(lon_deg_100km, decimals: 2)} (center of 100km square)")
if warning_100km != "", do: IO.puts("Warning: #{warning_100km}")

# Example from MGRS.cpp comments (polar)
mgrs_string_polar = "YJM1000099000" # A polar example (hypothetical)
mgrs_coords_polar = %MGRSCoordinates{mgrs_string: mgrs_string_polar}
%GeodeticCoordinates{latitude: lat_rad_polar, longitude: lon_rad_polar, warning_message: warning_polar} =
  MGRS.convert_to_geodetic(mgrs_coords_polar)
lat_deg_polar = Constants.radians_to_degrees(lat_rad_polar)
lon_deg_polar = Constants.radians_to_degrees(lon_rad_polar)
IO.puts("MGRS #{mgrs_string_polar} to GPS: Lat #{:erlang.float_to_binary(lat_deg_polar, decimals: 5)}, Lon #{:erlang.float_to_binary(lon_deg_polar, decimals: 5)}")
if warning_polar != "", do: IO.puts("Warning: #{warning_polar}")

```

**Expected Output for Example 2 (approximate):**

```
MGRS 18SUJ2337105498 to GPS: Lat 38.88950, Lon -77.03530
MGRS 18SUJ to GPS: Lat 38.87, Lon -77.32 (center of 100km square)
MGRS YJM1000099000 to GPS: Lat 84.89100, Lon -123.00000
```

*(Note: MGRS to GPS conversion typically yields the coordinate of the *south-west corner* of the grid cell represented by the MGRS string's precision. If no numeric part is given (e.g., "18SUJ"), it implies the SW corner of that 100km square. The exact values will depend on the implementation details matching the C++ reference.)*

## Error Handling

Functions in `CooCoo.MGRS` (and underlying modules) will raise an `ArgumentError` or a similar exception if invalid inputs are provided (e.g., latitude out of range, malformed MGRS string, invalid precision).

Some conversions might also populate the `warning_message` field in the returned coordinate struct for non-critical issues (e.g., MGRS string representing a point slightly outside its nominal latitude band but within an adjacent one).

## Modules

The library is structured as follows:

- `CooCoo.Constants`: Defines mathematical and geodetic constants (WGS84).
- `CooCoo.Coordinates`: Defines structs for different coordinate types:
  - `CooCoo.Coordinates.GeodeticCoordinates`
  - `CooCoo.Coordinates.MapProjectionCoordinates` (generic Easting/Northing)
  - `CooCoo.Coordinates.UTMCoordinates`
  - `CooCoo.Coordinates.UPSCoordinates`
  - `CooCoo.Coordinates.MGRSCoordinates`
- `CooCoo.Projections.TMHelpers`: Helper for Transverse Mercator coefficients.
- `CooCoo.Projections.TransverseMercator`: Core Transverse Mercator math.
- `CooCoo.Projections.PolarStereographic`: Core Polar Stereographic math.
- `CooCoo.UTM`: UTM projection specific logic, using Transverse Mercator.
- `CooCoo.UPS`: UPS projection specific logic, using Polar Stereographic.
- `CooCoo.MGRS`: Top-level MGRS logic, using UTM and UPS.

## Further Development

- Extensive testing against known test vectors (e.g., from GEOTRANS or other standard sources) is crucial.
- Refinement of error messages and warnings.
- Support for other ellipsoids (currently defaults to WGS84 logic where applicable).
