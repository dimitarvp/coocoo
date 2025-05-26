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

  @n_terms 6
  @max_delta_long Constants.pi_over_180() * 70.0

  @wgs84_coefficients_tuple TMHelpers.generate_coefficients(Constants.wgs84_f())
  @wgs84_a_coeffs elem(@wgs84_coefficients_tuple, 0)
  @wgs84_b_coeffs elem(@wgs84_coefficients_tuple, 1)
  @wgs84_r4oa elem(@wgs84_coefficients_tuple, 2)

  # --- Helper for ArcHyperbolicTangent ---
  defp atanh(x) when is_number(x) do
    cond do
      x == 1.0 ->
        :infinity

      x == -1.0 ->
        :"-infinity"

      abs(x) < 1.0 ->
        0.5 * :math.log((1.0 + x) / (1.0 - x))

      true ->
        # For abs(x) > 1.0, log input is negative, results in NaN or error.
        # This case should ideally not be reached with valid TM inputs.
        0.5 * :math.log((1.0 + x) / (1.0 - x))
    end
  end

  # --- Helper to compute geodetic latitude from conformal latitude ---
  defp geodetic_latitude_from_conformal(sin_chi, eccentricity) do
    s_initial = sin_chi

    final_s =
      do_geodetic_lat_iteration(sin_chi, eccentricity, s_initial, s_initial + 1.0e99, 30)

    clamped_final_s = max(-1.0, min(1.0, final_s))
    :math.asin(clamped_final_s)
  end

  defp do_geodetic_lat_iteration(_sin_chi, _eccentricity, current_s, _s_old, 0) do
    current_s
  end

  defp do_geodetic_lat_iteration(sin_chi, eccentricity, current_s, s_old, iterations_left) do
    if abs(current_s - s_old) < 1.0e-12 do
      current_s
    else
      val_for_atanh = eccentricity * current_s

      p_arg =
        cond do
          val_for_atanh >= 1.0 -> 1.0 - 1.0e-15
          val_for_atanh <= -1.0 -> -1.0 + 1.0e-15
          true -> val_for_atanh
        end

      p_exp_arg = eccentricity * atanh(p_arg)

      p =
        cond do
          # Approx limit for :math.exp
          p_exp_arg > 709.0 -> 1.0e300
          p_exp_arg < -709.0 -> 1.0e-300
          true -> :math.exp(p_exp_arg)
        end

      p_sq = p * p
      one_plus_sin_chi = 1.0 + sin_chi
      one_minus_sin_chi = 1.0 - sin_chi

      denominator = one_plus_sin_chi * p_sq + one_minus_sin_chi

      next_s =
        if denominator == 0.0 do
          current_s
        else
          (one_plus_sin_chi * p_sq - one_minus_sin_chi) / denominator
        end

      do_geodetic_lat_iteration(
        sin_chi,
        eccentricity,
        next_s,
        current_s,
        iterations_left - 1
      )
    end
  end

  # --- Helper to compute hyperbolic series terms ---
  defp compute_hyperbolic_series(two_x, num_terms) do
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
  defp compute_trig_series(two_y, num_terms) do
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

  @doc """
  Converts Geodetic coordinates (latitude, longitude) to Transverse Mercator
  easting and northing.
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

    {a_coeffs_tuple, _b_coeffs_tuple, r4oa_val} =
      if a == Constants.wgs84_a() and f == Constants.wgs84_f() do
        {@wgs84_a_coeffs, @wgs84_b_coeffs, @wgs84_r4oa}
      else
        TMHelpers.generate_coefficients(f)
      end

    k0r4 = r4oa_val * k0 * a

    lambda = geodetic_coords.longitude - origin_lon

    lambda =
      cond do
        lambda > Constants.pi() -> lambda - Constants.two_pi()
        lambda < -Constants.pi() -> lambda + Constants.two_pi()
        true -> lambda
      end

    # Warning for large lambda can be added here or in UTM module
    warning_message =
      if abs(lambda) > @max_delta_long,
        do: "Longitude is far from Central Meridian, distortion may be significant",
        else: ""

    sin_phi = :math.sin(geodetic_coords.latitude)
    cos_phi = :math.cos(geodetic_coords.latitude)

    val_for_atanh = es * sin_phi
    p_arg = if abs(val_for_atanh) >= 1.0, do: val_for_atanh / (1.0 + 1.0e-15), else: val_for_atanh

    p_exp_arg_fwd = es * atanh(p_arg)

    p_val =
      cond do
        p_exp_arg_fwd > 709.0 -> 1.0e300
        p_exp_arg_fwd < -709.0 -> 1.0e-300
        true -> :math.exp(p_exp_arg_fwd)
      end

    part1 = (1.0 + sin_phi) / p_val
    part2 = (1.0 - sin_phi) * p_val
    denom_chi = part1 + part2
    cos_chi = if denom_chi == 0.0, do: 0.0, else: 2.0 * cos_phi / denom_chi
    sin_chi = if denom_chi == 0.0, do: 0.0, else: (part1 - part2) / denom_chi

    u_prime = atanh(cos_chi * :math.sin(lambda))
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

    %MapProjectionCoordinates{
      easting: easting,
      northing: northing,
      warning_message: warning_message
    }
  end

  @doc """
  Converts Transverse Mercator easting and northing to Geodetic coordinates
  (latitude, longitude).
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

    {_a_coeffs_tuple, b_coeffs_tuple, r4oa_val} =
      if a == Constants.wgs84_a() and f == Constants.wgs84_f() do
        {@wgs84_a_coeffs, @wgs84_b_coeffs, @wgs84_r4oa}
      else
        TMHelpers.generate_coefficients(f)
      end

    k0r4_inv = if r4oa_val * k0 * a == 0.0, do: 0.0, else: 1.0 / (r4oa_val * k0 * a)

    x_star = (proj_coords.easting - false_easting) * k0r4_inv
    y_star = (proj_coords.northing - false_northing) * k0r4_inv

    {c2kx_list_tup, s2kx_list_tup} = compute_hyperbolic_series(2.0 * x_star, @n_terms)
    {c2ky_list_tup, s2ky_list_tup} = compute_trig_series(2.0 * y_star, @n_terms)

    u_prime =
      Enum.reduce(0..(@n_terms - 1), 0.0, fn k, acc ->
        acc + elem(b_coeffs_tuple, k) * elem(s2kx_list_tup, k) * elem(c2ky_list_tup, k)
      end) + x_star

    v_prime =
      Enum.reduce(0..(@n_terms - 1), 0.0, fn k, acc ->
        acc + elem(b_coeffs_tuple, k) * elem(c2kx_list_tup, k) * elem(s2ky_list_tup, k)
      end) + y_star

    lambda = :math.atan2(:math.sinh(u_prime), :math.cos(v_prime))

    cosh_u_prime = :math.cosh(u_prime)
    sin_chi = if cosh_u_prime == 0.0, do: 0.0, else: :math.sin(v_prime) / cosh_u_prime
    sin_chi = max(-1.0, min(1.0, sin_chi))

    latitude = geodetic_latitude_from_conformal(sin_chi, es)
    longitude = origin_lon + lambda

    longitude =
      cond do
        longitude > Constants.pi() -> longitude - Constants.two_pi()
        longitude < -Constants.pi() -> longitude + Constants.two_pi()
        true -> longitude
      end

    warning_message =
      if abs(lambda) > @max_delta_long,
        do: "Longitude is far from Central Meridian, distortion may be significant",
        else: ""

    %GeodeticCoordinates{
      longitude: longitude,
      latitude: latitude,
      height: 0.0,
      warning_message: warning_message
    }
  end
end
