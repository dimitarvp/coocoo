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

  @max_mgrs_precision Constants.max_mgrs_precision()

  @latitude_band_table [
    {:C, Constants.letter_value(:C), 1_100_000.0, -72.0, -80.5, 0.0},
    {:D, Constants.letter_value(:D), 2_000_000.0, -64.0, -72.0, 2_000_000.0},
    {:E, Constants.letter_value(:E), 2_800_000.0, -56.0, -64.0, 2_000_000.0},
    {:F, Constants.letter_value(:F), 3_700_000.0, -48.0, -56.0, 2_000_000.0},
    {:G, Constants.letter_value(:G), 4_600_000.0, -40.0, -48.0, 4_000_000.0},
    {:H, Constants.letter_value(:H), 5_500_000.0, -32.0, -40.0, 4_000_000.0},
    {:J, Constants.letter_value(:J), 6_400_000.0, -24.0, -32.0, 6_000_000.0},
    {:K, Constants.letter_value(:K), 7_300_000.0, -16.0, -24.0, 6_000_000.0},
    {:L, Constants.letter_value(:L), 8_200_000.0, -8.0, -16.0, 8_000_000.0},
    {:M, Constants.letter_value(:M), 9_100_000.0, 0.0, -8.0, 8_000_000.0},
    {:N, Constants.letter_value(:N), 0.0, 8.0, 0.0, 0.0},
    {:P, Constants.letter_value(:P), 800_000.0, 16.0, 8.0, 0.0},
    {:Q, Constants.letter_value(:Q), 1_700_000.0, 24.0, 16.0, 0.0},
    {:R, Constants.letter_value(:R), 2_600_000.0, 32.0, 24.0, 2_000_000.0},
    {:S, Constants.letter_value(:S), 3_500_000.0, 40.0, 32.0, 2_000_000.0},
    {:T, Constants.letter_value(:T), 4_400_000.0, 48.0, 40.0, 4_000_000.0},
    {:U, Constants.letter_value(:U), 5_300_000.0, 56.0, 48.0, 4_000_000.0},
    {:V, Constants.letter_value(:V), 6_200_000.0, 64.0, 56.0, 6_000_000.0},
    {:W, Constants.letter_value(:W), 7_000_000.0, 72.0, 64.0, 6_000_000.0},
    {:X, Constants.letter_value(:X), 7_900_000.0, 84.5, 72.0, 6_000_000.0}
  ]

  @ups_constant_table [
    {:A, Constants.letter_value(:A), :J, :Z, :Z, 800_000.0, 800_000.0},
    {:B, Constants.letter_value(:B), :A, :R, :Z, 2_000_000.0, 800_000.0},
    {:Y, Constants.letter_value(:Y), :J, :Z, :P, 800_000.0, 1_300_000.0},
    {:Z, Constants.letter_value(:Z), :A, :J, :P, 2_000_000.0, 1_300_000.0}
  ]

  # --- Helper Functions ---

  defp get_latitude_letter_code(latitude_rad) do
    cond do
      latitude_rad >= Constants.deg72_in_radians() and
          latitude_rad < Constants.deg84_5_in_radians() ->
        Constants.letter_value(:X)

      latitude_rad >= -Constants.deg80_5_in_radians() and
          latitude_rad < Constants.deg72_in_radians() ->
        band_float =
          (latitude_rad + Constants.deg80_in_radians()) / Constants.deg8_in_radians() + 1.0e-12

        band_index = trunc(band_float)
        band_index = if band_index < 0, do: 0, else: band_index

        if band_index >= 0 and band_index < length(@latitude_band_table) do
          elem(Enum.at(@latitude_band_table, band_index), 1)
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

  defp normalize_northing_loop(grid_northing, limit) when grid_northing >= limit do
    normalize_northing_loop(grid_northing - limit, limit)
  end

  defp normalize_northing_loop(grid_northing, _limit), do: grid_northing

  defp get_latitude_band_min_northing(letter_val) do
    found_band =
      Enum.find(@latitude_band_table, fn {_l_atom, l_val, _, _, _, _} -> l_val == letter_val end)

    if found_band do
      {_l_atom, _l_val, min_n, _n_deg, _s_deg, n_offset} = found_band
      {min_n, n_offset}
    else
      throw(%ArgumentError{
        message: "Invalid MGRS latitude band letter value: #{Constants.letter_char(letter_val)}"
      })
    end
  end

  defp in_latitude_range(letter_val, latitude_rad, border_rad) do
    # Ensure letter_val is valid before trying to find it
    # (C..H, J..N, P..X)
    valid_letter_val =
      (letter_val >= Constants.letter_value(:C) and letter_val <= Constants.letter_value(:H)) or
        (letter_val >= Constants.letter_value(:J) and letter_val <= Constants.letter_value(:N)) or
        (letter_val >= Constants.letter_value(:P) and letter_val <= Constants.letter_value(:X))

    if valid_letter_val do
      found_band =
        Enum.find(@latitude_band_table, fn {_l_atom, l_val, _, _, _, _} -> l_val == letter_val end)

      # This should always find a band if valid_letter_val is true, as our table covers these.
      if found_band do
        {_l_atom, _l_val, _min_n, north_deg, south_deg, _n_offset} = found_band
        north_boundary_rad = Constants.degrees_to_radians(north_deg)
        south_boundary_rad = Constants.degrees_to_radians(south_deg)

        south_boundary_rad - border_rad <= latitude_rad and
          latitude_rad <= north_boundary_rad + border_rad
      else
        # Should not be reached if valid_letter_val logic is correct and table is complete
        false
      end
    else
      # letter_val was outside the valid MGRS range (e.g. after decrementing past C)
      false
    end
  end

  defp break_mgrs_string(mgrs_string) do
    temp_mgrs_string = mgrs_string |> String.replace(" ", "") |> String.upcase()

    unless Regex.match?(~r/^[A-Z0-9]+$/i, temp_mgrs_string) do
      throw(%ArgumentError{message: "Invalid characters in MGRS string: " <> mgrs_string})
    end

    {zone_str_match, remaining_after_zone} =
      case Regex.run(~r/^(\d{1,2})/, temp_mgrs_string, return: :index) do
        [{start_idx, len} | _] when start_idx == 0 ->
          zone_s = String.slice(temp_mgrs_string, start_idx, len)
          remaining = String.slice(temp_mgrs_string, len, String.length(temp_mgrs_string) - len)
          {zone_s, remaining}

        _ ->
          {"", temp_mgrs_string}
      end

    zone = if zone_str_match == "", do: 0, else: String.to_integer(zone_str_match)

    if zone < 0 or zone > 60,
      do: throw(%ArgumentError{message: "Invalid zone in MGRS string: " <> zone_str_match})

    if String.length(remaining_after_zone) < 3 do
      throw(%ArgumentError{message: "MGRS string too short after zone: " <> remaining_after_zone})
    end

    ltr1_char = String.at(remaining_after_zone, 0)
    ltr2_char = String.at(remaining_after_zone, 1)
    ltr3_char = String.at(remaining_after_zone, 2)

    unless Regex.match?(~r/^[A-HJ-NP-X]$/, ltr1_char) and
             Regex.match?(~r/^[A-HJ-NP-Z]$/, ltr2_char) and
             Regex.match?(~r/^[A-HJ-NP-V]$/, ltr3_char) do
      throw(%ArgumentError{
        message: "Invalid MGRS letters: " <> ltr1_char <> ltr2_char <> ltr3_char
      })
    end

    letters = [
      Constants.letter_value(String.to_atom(ltr1_char)),
      Constants.letter_value(String.to_atom(ltr2_char)),
      Constants.letter_value(String.to_atom(ltr3_char))
    ]

    if zone == 0 do
      first_letter_as_char_atom = Constants.letter_char(List.first(letters)) |> String.to_atom()

      unless Enum.member?([:A, :B, :Y, :Z], first_letter_as_char_atom) do
        throw(%ArgumentError{message: "Invalid polar band letter for MGRS string without zone."})
      end
    end

    digits_str = String.slice(remaining_after_zone, 3, String.length(remaining_after_zone) - 3)
    num_digits = String.length(digits_str)

    if rem(num_digits, 2) != 0 or num_digits > 10 do
      throw(%ArgumentError{
        message: "Invalid number of digits for easting/northing in MGRS string: " <> digits_str
      })
    end

    precision = div(num_digits, 2)

    {easting, northing} =
      if num_digits > 0 do
        n_half = precision
        easting_str = String.slice(digits_str, 0, n_half)
        northing_str = String.slice(digits_str, n_half, n_half)
        multiplier = compute_scale(precision)

        {String.to_integer(easting_str) * multiplier,
         String.to_integer(northing_str) * multiplier}
      else
        {0.0, 0.0}
      end

    {zone, letters, easting, northing, precision}
  end

  # --- Main Public Conversion Functions ---

  @doc """
  Converts Geodetic coordinates to an MGRS coordinate string.
  Precision defaults to 5 (1-meter).
  """
  def convert_from_geodetic(%GeodeticCoordinates{} = geodetic_coords, precision \\ 5) do
    latitude_rad = geodetic_coords.latitude
    longitude_rad = geodetic_coords.longitude

    if latitude_rad < -Constants.pi_over_2() or latitude_rad > Constants.pi_over_2() do
      throw(%ArgumentError{message: "Latitude out of range [-90, 90] degrees"})
    end

    if precision < 0 or precision > @max_mgrs_precision do
      throw(%ArgumentError{message: "Precision must be between 0 and #{@max_mgrs_precision}"})
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

    grid_northing_temp = normalize_northing_loop(northing, @two_million)
    grid_northing_temp = grid_northing_temp + pattern_offset

    grid_northing_final =
      if grid_northing_temp >= @two_million,
        do: grid_northing_temp - @two_million,
        else: grid_northing_temp

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

  # --- MGRS to Geodetic Conversion ---
  @doc """
  Converts an MGRS coordinate string to Geodetic coordinates.
  """
  def convert_to_geodetic(%MGRSCoordinates{mgrs_string: mgrs_str_in}) do
    {zone, letters_indices, mgrs_easting, mgrs_northing, precision} =
      break_mgrs_string(mgrs_str_in)

    if zone > 0 do
      utm_coords_intermediate =
        to_utm_internal(zone, letters_indices, mgrs_easting, mgrs_northing, precision)

      UTM.convert_to_geodetic(utm_coords_intermediate)
    else
      ups_coords_intermediate = to_ups_internal(letters_indices, mgrs_easting, mgrs_northing)
      UPS.convert_to_geodetic(ups_coords_intermediate)
    end
  end

  defp to_utm_internal(zone, letters_indices, mgrs_easting, mgrs_northing, precision) do
    first_letter_val = Enum.at(letters_indices, 0)
    second_letter_val = Enum.at(letters_indices, 1)
    third_letter_val = Enum.at(letters_indices, 2)

    if first_letter_val == Constants.letter_value(:X) and Enum.member?([32, 34, 36], zone) do
      throw(%ArgumentError{message: "Invalid MGRS string (X band rule for zone #{zone})"})
    end

    if first_letter_val == Constants.letter_value(:V) and zone == 31 and
         second_letter_val > Constants.letter_value(:D) do
      throw(%ArgumentError{message: "Invalid MGRS string (V band rule for zone 31)"})
    end

    hemisphere = if first_letter_val < Constants.letter_value(:N), do: ?S, else: ?N

    {ltr2_low_value, ltr2_high_value, pattern_offset} = get_grid_values(zone)

    if second_letter_val < ltr2_low_value or
         second_letter_val > ltr2_high_value or
         third_letter_val > Constants.letter_value(:V) do
      throw(%ArgumentError{message: "Invalid MGRS grid letters for zone #{zone}"})
    end

    grid_easting = (second_letter_val - ltr2_low_value + 1) * @one_hundred_thousand

    grid_easting =
      if ltr2_low_value == Constants.letter_value(:J) and
           second_letter_val > Constants.letter_value(:O) do
        grid_easting - @one_hundred_thousand
      else
        grid_easting
      end

    row_letter_northing = third_letter_val * @one_hundred_thousand

    row_letter_northing =
      if third_letter_val > Constants.letter_value(:O),
        do: row_letter_northing - @one_hundred_thousand,
        else: row_letter_northing

    row_letter_northing =
      if third_letter_val > Constants.letter_value(:I),
        do: row_letter_northing - @one_hundred_thousand,
        else: row_letter_northing

    row_letter_northing =
      if row_letter_northing >= @two_million,
        do: row_letter_northing - @two_million,
        else: row_letter_northing

    {min_northing_for_band, northing_offset_for_band} =
      get_latitude_band_min_northing(first_letter_val)

    grid_northing = row_letter_northing - pattern_offset
    grid_northing = if grid_northing < 0.0, do: grid_northing + @two_million, else: grid_northing
    grid_northing = grid_northing + northing_offset_for_band

    grid_northing =
      if grid_northing < min_northing_for_band,
        do: grid_northing + @two_million,
        else: grid_northing

    easting = grid_easting + mgrs_easting
    northing = grid_northing + mgrs_northing

    temp_utm_for_validation = %UTMCoordinates{
      zone: zone,
      hemisphere: hemisphere,
      easting: easting,
      northing: northing
    }

    %GeodeticCoordinates{latitude: calculated_latitude_rad} =
      UTM.convert_to_geodetic(temp_utm_for_validation)

    divisor_for_border = :math.pow(10.0, precision)

    border_rad =
      if divisor_for_border == 0.0, do: 0.0, else: Constants.pi_over_180() / divisor_for_border

    warning_message =
      if in_latitude_range(first_letter_val, calculated_latitude_rad, border_rad) do
        ""
      else
        prev_band_val_initial = first_letter_val - 1

        prev_band_val =
          if prev_band_val_initial == Constants.letter_value(:I) or
               prev_band_val_initial == Constants.letter_value(:O) do
            prev_band_val_initial - 1
          else
            prev_band_val_initial
          end

        prev_band_val =
          if first_letter_val == Constants.letter_value(:C),
            do: first_letter_val,
            else: prev_band_val

        next_band_val_initial = first_letter_val + 1

        next_band_val =
          if next_band_val_initial == Constants.letter_value(:I) or
               next_band_val_initial == Constants.letter_value(:O) do
            next_band_val_initial + 1
          else
            next_band_val_initial
          end

        next_band_val =
          if first_letter_val == Constants.letter_value(:X),
            do: first_letter_val,
            else: next_band_val

        if in_latitude_range(prev_band_val, calculated_latitude_rad, border_rad) or
             in_latitude_range(next_band_val, calculated_latitude_rad, border_rad) do
          "MGRS coordinates are in an adjacent latitude band."
        else
          throw(%ArgumentError{message: "MGRS coordinates are outside valid latitude bands."})
        end
      end

    %UTMCoordinates{
      zone: zone,
      hemisphere: hemisphere,
      easting: easting,
      northing: northing,
      warning_message: warning_message
    }
  end

  defp to_ups_internal(letters_indices, mgrs_easting, mgrs_northing) do
    first_letter_val = Enum.at(letters_indices, 0)
    second_letter_val = Enum.at(letters_indices, 1)
    third_letter_val = Enum.at(letters_indices, 2)

    found_ups_const =
      Enum.find(@ups_constant_table, fn {_atom, val, _, _, _, _, _} -> val == first_letter_val end)

    if is_nil(found_ups_const) do
      throw(%ArgumentError{
        message: "Invalid MGRS polar band letter: #{Constants.letter_char(first_letter_val)}"
      })
    end

    {_band_atom, _band_val, ltr2_low_atom, ltr2_high_atom, ltr3_high_atom, false_e, false_n} =
      found_ups_const

    hemisphere =
      if first_letter_val == Constants.letter_value(:Y) or
           first_letter_val == Constants.letter_value(:Z),
         do: ?N,
         else: ?S

    ltr2_low_value = Constants.letter_value(ltr2_low_atom)
    ltr2_high_value = Constants.letter_value(ltr2_high_atom)
    ltr3_high_value = Constants.letter_value(ltr3_high_atom)

    skipped_grid_letters = [:D, :E, :M, :N, :V, :W] |> Enum.map(&Constants.letter_value(&1))

    if second_letter_val < ltr2_low_value or
         second_letter_val > ltr2_high_value or
         Enum.member?(skipped_grid_letters, second_letter_val) or
         third_letter_val > ltr3_high_value do
      throw(%ArgumentError{message: "Invalid MGRS UPS grid letters"})
    end

    grid_northing = third_letter_val * @one_hundred_thousand + false_n

    grid_northing =
      if third_letter_val > Constants.letter_value(:I),
        do: grid_northing - @one_hundred_thousand,
        else: grid_northing

    grid_northing =
      if third_letter_val > Constants.letter_value(:O),
        do: grid_northing - @one_hundred_thousand,
        else: grid_northing

    base_grid_easting = (second_letter_val - ltr2_low_value) * @one_hundred_thousand + false_e

    # Y, Z bands
    # A, B bands
    current_ge =
      if ltr2_low_value != Constants.letter_value(:A) do
        temp_ge = base_grid_easting

        temp_ge =
          if second_letter_val > Constants.letter_value(:L),
            do: temp_ge - 300_000.0,
            else: temp_ge

        if second_letter_val > Constants.letter_value(:U), do: temp_ge - 200_000.0, else: temp_ge
      else
        temp_ge = base_grid_easting

        temp_ge =
          if second_letter_val > Constants.letter_value(:C),
            do: temp_ge - 200_000.0,
            else: temp_ge

        temp_ge =
          if second_letter_val > Constants.letter_value(:I),
            do: temp_ge - @one_hundred_thousand,
            else: temp_ge

        if second_letter_val > Constants.letter_value(:L), do: temp_ge - 300_000.0, else: temp_ge
      end

    easting = current_ge + mgrs_easting
    northing = grid_northing + mgrs_northing

    %UPSCoordinates{
      hemisphere: hemisphere,
      easting: easting,
      northing: northing,
      warning_message: ""
    }
  end
end
