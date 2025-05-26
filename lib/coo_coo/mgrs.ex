defmodule CooCoo.MGRS do
  @moduledoc """
  Provides conversions between Geodetic (latitude/longitude) coordinates
  and Military Grid Reference System (MGRS) coordinates.
  """

  alias CooCoo.Constants
  alias CooCoo.Coordinates.GeodeticCoordinates
  alias CooCoo.Coordinates.UTMCoordinates
  alias CooCoo.Coordinates.UPSCoordinates
  alias CooCoo.Coordinates.MGRSCoordinates
  alias CooCoo.UTM
  alias CooCoo.UPS

  @one_hundred_thousand Constants.one_hundred_thousand()
  @two_million Constants.two_million()
  @mgrs_epsilon2 Constants.mgrs_epsilon2()

  # Define @max_mgrs_precision using the constant from CooCoo.Constants
  @max_mgrs_precision Constants.max_mgrs_precision()

  @latitude_band_table [
    {:C, 1_100_000.0, -72.0, -80.5, 0.0},
    {:D, 2_000_000.0, -64.0, -72.0, 2_000_000.0},
    {:E, 2_800_000.0, -56.0, -64.0, 2_000_000.0},
    {:F, 3_700_000.0, -48.0, -56.0, 2_000_000.0},
    {:G, 4_600_000.0, -40.0, -48.0, 4_000_000.0},
    {:H, 5_500_000.0, -32.0, -40.0, 4_000_000.0},
    {:J, 6_400_000.0, -24.0, -32.0, 6_000_000.0},
    {:K, 7_300_000.0, -16.0, -24.0, 6_000_000.0},
    {:L, 8_200_000.0, -8.0, -16.0, 8_000_000.0},
    {:M, 9_100_000.0, 0.0, -8.0, 8_000_000.0},
    {:N, 0.0, 8.0, 0.0, 0.0},
    {:P, 800_000.0, 16.0, 8.0, 0.0},
    {:Q, 1_700_000.0, 24.0, 16.0, 0.0},
    {:R, 2_600_000.0, 32.0, 24.0, 2_000_000.0},
    {:S, 3_500_000.0, 40.0, 32.0, 2_000_000.0},
    {:T, 4_400_000.0, 48.0, 40.0, 4_000_000.0},
    {:U, 5_300_000.0, 56.0, 48.0, 4_000_000.0},
    {:V, 6_200_000.0, 64.0, 56.0, 6_000_000.0},
    {:W, 7_000_000.0, 72.0, 64.0, 6_000_000.0},
    {:X, 7_900_000.0, 84.5, 72.0, 6_000_000.0}
  ]

  @ups_constant_table [
    {:A, :J, :Z, :Z, 800_000.0, 800_000.0},
    {:B, :A, :R, :Z, 2_000_000.0, 800_000.0},
    {:Y, :J, :Z, :P, 800_000.0, 1_300_000.0},
    {:Z, :A, :J, :P, 2_000_000.0, 1_300_000.0}
  ]

  defp get_latitude_letter_code(latitude_rad) do
    cond do
      latitude_rad >= Constants.deg72_in_radians() and
          latitude_rad < Constants.deg84_5_in_radians() ->
        Constants.letter_value(:X)

      # Changed > to >= for -80.5
      latitude_rad >= -Constants.deg80_5_in_radians() and
          latitude_rad < Constants.deg72_in_radians() ->
        band_float =
          (latitude_rad + Constants.deg80_in_radians()) / Constants.deg8_in_radians() + 1.0e-12

        band_index = trunc(band_float)
        band_index = if band_index < 0, do: 0, else: band_index

        if band_index >= 0 and band_index < length(@latitude_band_table) do
          # The table is ordered C..X. The index maps to this.
          # E.g. index 0 = C, 1 = D ... 7 = L, 8 = M ... 18 = W, 19 = X
          elem(Enum.at(@latitude_band_table, band_index), 0) |> Constants.letter_value()
        else
          throw(%ArgumentError{
            message: "Latitude out of calculable band range for MGRS (index: #{band_index})"
          })
        end

      true ->
        throw(%ArgumentError{message: "Latitude out of MGRS range (lat_rad: #{latitude_rad})"})
    end
  end

  defp compute_scale(precision) when precision >= 0 and precision <= @max_mgrs_precision do
    :math.pow(10.0, 5 - precision)
  end

  defp compute_scale(_precision) do
    throw(%ArgumentError{message: "Invalid MGRS precision"})
  end

  defp make_mgrs_string(zone_num, letters_indices, easting, northing, precision) do
    zone_str =
      if zone_num > 0 do
        :io_lib.format("~2.2.0w", [zone_num]) |> List.to_string()
      else
        "  "
      end

    letter_chars = Enum.map(letters_indices, &Constants.letter_char(&1)) |> Enum.join("")

    divisor = compute_scale(precision)

    easting_within_100k = :math.fmod(easting, @one_hundred_thousand)

    easting_within_100k =
      if easting_within_100k >= 99999.5, do: 99999.0, else: easting_within_100k

    east_val = trunc((easting_within_100k + @mgrs_epsilon2) / divisor)

    northing_within_100k = :math.fmod(northing, @one_hundred_thousand)

    northing_within_100k =
      if northing_within_100k >= 99999.5, do: 99999.0, else: northing_within_100k

    north_val = trunc((northing_within_100k + @mgrs_epsilon2) / divisor)

    east_str = Integer.to_string(east_val) |> String.pad_leading(precision, "0")
    north_str = Integer.to_string(north_val) |> String.pad_leading(precision, "0")

    zone_str <> letter_chars <> east_str <> north_str
  end

  defp get_grid_values(zone) do
    set_number = rem(zone, 6)
    set_number = if set_number == 0, do: 6, else: set_number

    pattern_offset =
      if rem(set_number, 2) == 0, do: 500_000.0, else: 0.0

    {ltr2_low_value, ltr2_high_value} =
      cond do
        set_number == 1 or set_number == 4 ->
          {Constants.letter_value(:A), Constants.letter_value(:H)}

        set_number == 2 or set_number == 5 ->
          {Constants.letter_value(:J), Constants.letter_value(:R)}

        set_number == 3 or set_number == 6 ->
          {Constants.letter_value(:S), Constants.letter_value(:Z)}
      end

    {ltr2_low_value, ltr2_high_value, pattern_offset}
  end

  # --- Recursive helper for while loop like behavior ---
  defp normalize_northing_loop(grid_northing, limit) when grid_northing >= limit do
    normalize_northing_loop(grid_northing - limit, limit)
  end

  defp normalize_northing_loop(grid_northing, _limit), do: grid_northing

  @doc """
  Converts Geodetic coordinates to an MGRS coordinate string.
  """
  def convert_from_geodetic(%GeodeticCoordinates{} = geodetic_coords, precision \\ 5) do
    latitude_rad = geodetic_coords.latitude
    # Retain for passing to from_utm_internal
    longitude_rad = geodetic_coords.longitude

    if latitude_rad < -Constants.pi_over_2() or latitude_rad > Constants.pi_over_2() do
      throw(%ArgumentError{message: "Latitude out of range [-90, 90] degrees"})
    end

    if precision < 0 or precision > @max_mgrs_precision do
      throw(%ArgumentError{message: "Precision must be between 0 and 5"})
    end

    if latitude_rad >= Constants.min_mgrs_non_polar_lat() - Constants.mgrs_epsilon2() and
         latitude_rad < Constants.max_mgrs_non_polar_lat() + Constants.mgrs_epsilon2() do
      utm_coords = UTM.convert_from_geodetic(geodetic_coords)
      from_utm_internal(utm_coords, longitude_rad, latitude_rad, precision)
    else
      ups_coords = UPS.convert_from_geodetic(geodetic_coords)
      from_ups_internal(ups_coords, precision)
    end
  end

  # longitude_rad parameter is kept for potential future use if natural zone calculation/override needed here
  defp from_utm_internal(%UTMCoordinates{} = utm_coords, _longitude_rad, latitude_rad, precision) do
    letters = List.duplicate(0, 3)

    letters = List.replace_at(letters, 0, get_latitude_letter_code(latitude_rad))

    zone_num = utm_coords.zone
    easting = utm_coords.easting
    northing = utm_coords.northing

    {ltr2_low_value, _ltr2_high_value, pattern_offset} = get_grid_values(zone_num)

    easting_100k_index_raw = trunc(easting / @one_hundred_thousand) - 1
    letter1_val = ltr2_low_value + easting_100k_index_raw

    letter1_val =
      if ltr2_low_value == Constants.letter_value(:J) and
           letter1_val > Constants.letter_value(:N),
         do: letter1_val + 1,
         else: letter1_val

    letters = List.replace_at(letters, 1, letter1_val)

    grid_northing_intermediate = normalize_northing_loop(northing, @two_million)
    grid_northing_final = grid_northing_intermediate + pattern_offset

    grid_northing_final =
      if grid_northing_final >= @two_million,
        do: grid_northing_final - @two_million,
        else: grid_northing_final

    letter2_val = trunc(grid_northing_final / @one_hundred_thousand)

    letter2_val =
      if letter2_val > Constants.letter_value(:H), do: letter2_val + 1, else: letter2_val

    letter2_val =
      if letter2_val > Constants.letter_value(:N), do: letter2_val + 1, else: letter2_val

    letters = List.replace_at(letters, 2, letter2_val)

    mgrs_str = make_mgrs_string(zone_num, letters, easting, northing, precision)
    %MGRSCoordinates{mgrs_string: mgrs_str, precision: precision}
  end

  defp from_ups_internal(%UPSCoordinates{} = ups_coords, precision) do
    letters = List.duplicate(0, 3)
    easting = ups_coords.easting
    northing = ups_coords.northing
    hemisphere = ups_coords.hemisphere

    {first_letter_atom, ltr2_low_atom, _ltr2_high_atom, _ltr3_high_atom, false_e, false_n} =
      cond do
        hemisphere == ?N and easting >= @two_million -> List.keyfind!(@ups_constant_table, :Z, 0)
        hemisphere == ?N -> List.keyfind!(@ups_constant_table, :Y, 0)
        hemisphere == ?S and easting >= @two_million -> List.keyfind!(@ups_constant_table, :B, 0)
        hemisphere == ?S -> List.keyfind!(@ups_constant_table, :A, 0)
      end

    letters = List.replace_at(letters, 0, Constants.letter_value(first_letter_atom))
    ltr2_low_value = Constants.letter_value(ltr2_low_atom)

    grid_easting = easting - false_e
    letter1_val = ltr2_low_value + trunc(grid_easting / @one_hundred_thousand)

    letter1_val =
      if easting < @two_million do
        temp_val =
          if letter1_val > Constants.letter_value(:L), do: letter1_val + 3, else: letter1_val

        if temp_val > Constants.letter_value(:U), do: temp_val + 2, else: temp_val
      else
        temp_val =
          if letter1_val > Constants.letter_value(:C), do: letter1_val + 2, else: letter1_val

        temp_val = if temp_val > Constants.letter_value(:H), do: temp_val + 1, else: temp_val
        if temp_val > Constants.letter_value(:L), do: temp_val + 3, else: temp_val
      end

    letters = List.replace_at(letters, 1, letter1_val)

    grid_northing = northing - false_n
    letter2_val = trunc(grid_northing / @one_hundred_thousand)

    letter2_val =
      if letter2_val > Constants.letter_value(:H), do: letter2_val + 1, else: letter2_val

    letter2_val =
      if letter2_val > Constants.letter_value(:N), do: letter2_val + 1, else: letter2_val

    letters = List.replace_at(letters, 2, letter2_val)

    mgrs_str = make_mgrs_string(0, letters, easting, northing, precision)
    %MGRSCoordinates{mgrs_string: mgrs_str, precision: precision}
  end
end
