defmodule CooCoo.Constants do
  @moduledoc """
  Module holding various mathematical and geodetic constants.
  """

  # Mathematical Constants
  @pi :math.pi()
  @pi_over_2 @pi / 2.0
  @pi_over_4 @pi / 4.0
  @pi_over_180 @pi / 180.0
  @two_pi 2.0 * @pi
  # For clarity if preferred over @pi_over_180
  @one_degree_in_radians @pi_over_180

  # WGS84 Ellipsoid Parameters (Most commonly used for GPS and often MGRS)
  # Source: NGA.STND.0036_1.0.0_WGS84 (defines a and 1/f)
  # Semi-major axis in meters
  @wgs84_a 6_378_137.0
  # Inverse flattening
  @wgs84_inv_f 298.257223563
  # Flattening
  @wgs84_f 1.0 / @wgs84_inv_f
  # Eccentricity squared (es2) = 2f - f^2
  @wgs84_es2 2.0 * @wgs84_f - @wgs84_f * @wgs84_f
  # Eccentricity (es) = sqrt(es2)
  @wgs84_es :math.sqrt(@wgs84_es2)

  # MGRS/USNG Specific Numerical Constants from MGRS.cpp
  @one_hundred_thousand 100_000.0
  @two_million 2_000_000.0

  # Maximum precision of easting & northing for MGRS/USNG
  @max_mgrs_precision 5

  # Latitude boundaries for MGRS (non-polar vs polar) in radians
  # MIN_MGRS_NON_POLAR_LAT  (-80.0 * ( PI / 180.0 )) /* -80 deg in rad */
  # MAX_MGRS_NON_POLAR_LAT  ( 84.0 * ( PI / 180.0 )) /*  84 deg in rad */
  @min_mgrs_non_polar_lat -80.0 * @pi_over_180
  @max_mgrs_non_polar_lat 84.0 * @pi_over_180

  # Specific degree to radian conversions often used in MGRS logic
  # In Elixir, it's often cleaner to calculate these on the fly or via a helper,
  # but for direct translation:
  @deg6_in_radians 6.0 * @pi_over_180
  @deg8_in_radians 8.0 * @pi_over_180
  @deg72_in_radians 72.0 * @pi_over_180
  @deg80_in_radians 80.0 * @pi_over_180
  # Used for latitude band C southern boundary
  @deg80_5_in_radians 80.5 * @pi_over_180
  # Used for latitude band X northern boundary
  @deg84_5_in_radians 84.5 * @pi_over_180

  # MGRS Related constants from MGRS.cpp for makeMGRSString / breakMGRSString
  # Used for rounding eastings/northings
  @mgrs_epsilon2 4.99e-4

  # UTM constants that MGRS uses (these are also somewhat general but appear in MGRS context)
  @utm_min_easting 100_000.0
  # Max easting for UTM (MGRS uses a slightly different check internally sometimes)
  @utm_max_easting_non_mgrs 900_000.0
  @utm_min_northing 0.0
  # Max northing for UTM
  @utm_max_northing_non_mgrs 10_000_000.0

  # UPS constants for MGRS (these are for polar regions)
  @ups_min_east_north 0.0
  # MAX_EAST_NORTH in MGRS.cpp seems to be this for UPS based coords
  @ups_max_east_north_mgrs 3_999_999.0

  # Letter representations (0-indexed, A=0) - these are heavily used in MGRS logic
  # These will likely be used in pattern matching or case statements rather than direct array indexing.
  @letter_a 0
  @letter_b 1
  @letter_c 2
  @letter_d 3
  @letter_e 4
  @letter_f 5
  @letter_g 6
  @letter_h 7
  # @letter_i 8  (I is skipped in MGRS)
  @letter_j 9
  @letter_k 10
  @letter_l 11
  @letter_m 12
  @letter_n 13
  # @letter_o 14 (O is skipped in MGRS)
  @letter_p 15
  @letter_q 16
  @letter_r 17
  @letter_s 18
  @letter_t 19
  @letter_u 20
  @letter_v 21
  @letter_w 22
  @letter_x 23
  @letter_y 24
  @letter_z 25

  # Alphabet for MGRS string construction
  @alphabet "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

  # Public API for accessing constants
  def pi, do: @pi
  def pi_over_2, do: @pi_over_2
  def pi_over_4, do: @pi_over_4
  def pi_over_180, do: @pi_over_180
  def two_pi, do: @two_pi
  def one_degree_in_radians, do: @one_degree_in_radians

  def wgs84_a, do: @wgs84_a
  def wgs84_f, do: @wgs84_f
  def wgs84_inv_f, do: @wgs84_inv_f
  def wgs84_es, do: @wgs84_es
  def wgs84_es2, do: @wgs84_es2

  def one_hundred_thousand, do: @one_hundred_thousand
  def two_million, do: @two_million
  def max_mgrs_precision, do: @max_mgrs_precision
  def min_mgrs_non_polar_lat, do: @min_mgrs_non_polar_lat
  def max_mgrs_non_polar_lat, do: @max_mgrs_non_polar_lat

  def deg6_in_radians, do: @deg6_in_radians
  def deg8_in_radians, do: @deg8_in_radians
  def deg72_in_radians, do: @deg72_in_radians
  def deg80_in_radians, do: @deg80_in_radians
  def deg80_5_in_radians, do: @deg80_5_in_radians
  def deg84_5_in_radians, do: @deg84_5_in_radians

  def mgrs_epsilon2, do: @mgrs_epsilon2

  def utm_min_easting, do: @utm_min_easting
  def utm_max_easting_non_mgrs, do: @utm_max_easting_non_mgrs
  def utm_min_northing, do: @utm_min_northing
  def utm_max_northing_non_mgrs, do: @utm_max_northing_non_mgrs

  def ups_min_east_north, do: @ups_min_east_north
  def ups_max_east_north_mgrs, do: @ups_max_east_north_mgrs

  # Alphabet related
  def alphabet_string, do: @alphabet
  def letter_value(char) when is_integer(char) and char >= 0 and char <= 25, do: char

  def letter_value(char_atom) when is_atom(char_atom) do
    case Atom.to_string(char_atom) do
      "A" -> @letter_a
      "B" -> @letter_b
      "C" -> @letter_c
      "D" -> @letter_d
      "E" -> @letter_e
      "F" -> @letter_f
      "G" -> @letter_g
      "H" -> @letter_h
      "J" -> @letter_j
      "K" -> @letter_k
      "L" -> @letter_l
      "M" -> @letter_m
      "N" -> @letter_n
      "P" -> @letter_p
      "Q" -> @letter_q
      "R" -> @letter_r
      "S" -> @letter_s
      "T" -> @letter_t
      "U" -> @letter_u
      "V" -> @letter_v
      "W" -> @letter_w
      "X" -> @letter_x
      "Y" -> @letter_y
      "Z" -> @letter_z
      # Or raise error
      _ -> nil
    end
  end

  def letter_char(value) when is_integer(value) and value >= 0 and value <= 25 do
    String.at(@alphabet, value)
  end

  @doc """
  Converts degrees to radians.
  """
  def degrees_to_radians(degrees) when is_number(degrees) do
    degrees * @pi_over_180
  end

  @doc """
  Converts radians to degrees.
  """
  def radians_to_degrees(radians) when is_number(radians) do
    radians / @pi_over_180
  end
end
