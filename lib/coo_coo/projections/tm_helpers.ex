defmodule CooCoo.Projections.TMHelpers do
  @moduledoc """
  Helper functions for Transverse Mercator projection coefficient generation.
  These are based on Krueger's n-series expansions.
  """

  # Number of terms to use for alpha and beta series
  @n_terms 6

  # --- Private helpers for individual alpha coefficients (alpha_k) ---
  # Each function calculates one term of the series for x* and y*
  # Arguments are powers of n: n1=n, n2=n^2, etc.

  def calculate_alpha_1(n1, n2, n3, n4, n5, n6, n7, n8) do
    0.5 * n1 - 2.0 / 3.0 * n2 + 5.0 / 16.0 * n3 + 41.0 / 180.0 * n4 - 127.0 / 288.0 * n5 +
      7891.0 / 37800.0 * n6 - 72161.0 / 387_072.0 * n7 + 18_975_107.0 / 50_803_200.0 * n8
  end

  def calculate_alpha_2(_n1, n2, n3, n4, n5, n6, n7, n8) do
    13.0 / 48.0 * n2 - 3.0 / 5.0 * n3 + 557.0 / 1440.0 * n4 + 281.0 / 630.0 * n5 -
      1_983_433.0 / 1_935_360.0 * n6 + 13769.0 / 28800.0 * n7 - 148_003_883.0 / 174_182_400.0 * n8
  end

  def calculate_alpha_3(_n1, _n2, n3, n4, n5, n6, n7, n8) do
    61.0 / 240.0 * n3 - 103.0 / 140.0 * n4 + 15061.0 / 26880.0 * n5 + 167_603.0 / 181_440.0 * n6 -
      67_102_379.0 / 29_030_400.0 * n7 + 79_682_431.0 / 79_833_600.0 * n8
  end

  def calculate_alpha_4(_n1, _n2, _n3, n4, n5, n6, n7, n8) do
    49561.0 / 161_280.0 * n4 - 179.0 / 168.0 * n5 + 6_601_661.0 / 7_257_600.0 * n6 +
      97445.0 / 49896.0 * n7 - 40_176_129_013.0 / 7_664_025_600.0 * n8
  end

  def calculate_alpha_5(_n1, _n2, _n3, _n4, n5, n6, n7, n8) do
    34729.0 / 80640.0 * n5 - 3_418_889.0 / 1_995_840.0 * n6 + 14_644_087.0 / 9_123_840.0 * n7 +
      2_605_413_599.0 / 622_702_080.0 * n8
  end

  def calculate_alpha_6(_n1, _n2, _n3, _n4, _n5, n6, n7, n8) do
    212_378_941.0 / 319_334_400.0 * n6 - 30_705_481.0 / 10_378_368.0 * n7 +
      175_214_326_799.0 / 58_118_860_800.0 * n8
  end

  def calculate_alpha_7(_n1, _n2, _n3, _n4, _n5, _n6, n7, n8) do
    1_522_256_789.0 / 1_383_782_400.0 * n7 - 16_759_934_899.0 / 3_113_510_400.0 * n8
  end

  def calculate_alpha_8(_n1, _n2, _n3, _n4, _n5, _n6, _n7, n8) do
    1_424_729_850_961.0 / 743_921_418_240.0 * n8
  end

  # --- Private helpers for individual beta coefficients (beta_k) ---
  # Each function calculates one term of the series for inverse conversion

  def calculate_beta_1(n1, n2, n3, n4, n5, n6, n7, n8) do
    -0.5 * n1 + 2.0 / 3.0 * n2 - 37.0 / 96.0 * n3 + 1.0 / 360.0 * n4 + 81.0 / 512.0 * n5 -
      96199.0 / 604_800.0 * n6 + 5_406_467.0 / 38_707_200.0 * n7 - 7_944_359.0 / 67_737_600.0 * n8
  end

  def calculate_beta_2(_n1, n2, n3, n4, n5, n6, n7, n8) do
    -(1.0 / 48.0) * n2 - 1.0 / 15.0 * n3 + 437.0 / 1440.0 * n4 - 46.0 / 105.0 * n5 +
      1_118_711.0 / 3_870_720.0 * n6 - 51841.0 / 1_209_600.0 * n7 -
      24_749_483.0 / 348_364_800.0 * n8
  end

  def calculate_beta_3(_n1, _n2, n3, n4, n5, n6, n7, n8) do
    -(17.0 / 480.0) * n3 + 37.0 / 840.0 * n4 + 209.0 / 4480.0 * n5 - 5569.0 / 90720.0 * n6 -
      9_261_899.0 / 58_060_800.0 * n7 + 6_457_463.0 / 17_740_800.0 * n8
  end

  def calculate_beta_4(_n1, _n2, _n3, n4, n5, n6, n7, n8) do
    -(4397.0 / 161_280.0) * n4 + 11.0 / 504.0 * n5 + 830_251.0 / 7_257_600.0 * n6 -
      466_511.0 / 2_494_800.0 * n7 - 324_154_477.0 / 7_664_025_600.0 * n8
  end

  def calculate_beta_5(_n1, _n2, _n3, _n4, n5, n6, n7, n8) do
    -(4583.0 / 161_280.0) * n5 + 108_847.0 / 3_991_680.0 * n6 + 8_005_831.0 / 63_866_880.0 * n7 -
      22_894_433.0 / 124_540_416.0 * n8
  end

  def calculate_beta_6(_n1, _n2, _n3, _n4, _n5, n6, n7, n8) do
    -(20_648_693.0 / 638_668_800.0) * n6 + 16_363_163.0 / 518_918_400.0 * n7 +
      2_204_645_983.0 / 12_915_302_400.0 * n8
  end

  def calculate_beta_7(_n1, _n2, _n3, _n4, _n5, _n6, n7, n8) do
    -(219_941_297.0 / 5_535_129_600.0) * n7 - 497_323_811.0 / 12_454_041_600.0 * n8
  end

  def calculate_beta_8(_n1, _n2, _n3, _n4, _n5, _n6, _n7, n8) do
    -(191_773_887_257.0 / 3_719_607_091_200.0) * n8
  end

  # --- Public function to generate all coefficients ---
  @doc """
  Generates Krueger's n-series coefficients for Transverse Mercator.
  Returns a tuple: `{a_coeffs_tuple, b_coeffs_tuple, r4oa}`
  where `a_coeffs_tuple` and `b_coeffs_tuple` contain `@n_terms` elements.
  """
  def generate_coefficients(flattening) do
    inv_f = 1.0 / flattening
    # Helmert's n (denoted as n in NGA document)
    n1 = 1.0 / (2.0 * inv_f - 1.0)

    n2 = n1 * n1
    n3 = n2 * n1
    n4 = n3 * n1
    n5 = n4 * n1
    n6 = n5 * n1
    n7 = n6 * n1
    n8 = n7 * n1
    # n9  = n8 * n1 # Not used in the 8-term series shown
    # n^10
    n10 = n8 * n2

    all_a_coeffs = [
      calculate_alpha_1(n1, n2, n3, n4, n5, n6, n7, n8),
      calculate_alpha_2(n1, n2, n3, n4, n5, n6, n7, n8),
      calculate_alpha_3(n1, n2, n3, n4, n5, n6, n7, n8),
      calculate_alpha_4(n1, n2, n3, n4, n5, n6, n7, n8),
      calculate_alpha_5(n1, n2, n3, n4, n5, n6, n7, n8),
      calculate_alpha_6(n1, n2, n3, n4, n5, n6, n7, n8),
      calculate_alpha_7(n1, n2, n3, n4, n5, n6, n7, n8),
      calculate_alpha_8(n1, n2, n3, n4, n5, n6, n7, n8)
    ]

    all_b_coeffs = [
      calculate_beta_1(n1, n2, n3, n4, n5, n6, n7, n8),
      calculate_beta_2(n1, n2, n3, n4, n5, n6, n7, n8),
      calculate_beta_3(n1, n2, n3, n4, n5, n6, n7, n8),
      calculate_beta_4(n1, n2, n3, n4, n5, n6, n7, n8),
      calculate_beta_5(n1, n2, n3, n4, n5, n6, n7, n8),
      calculate_beta_6(n1, n2, n3, n4, n5, n6, n7, n8),
      calculate_beta_7(n1, n2, n3, n4, n5, n6, n7, n8),
      calculate_beta_8(n1, n2, n3, n4, n5, n6, n7, n8)
    ]

    a_coeffs_final =
      all_a_coeffs
      |> Enum.take(@n_terms)
      |> List.to_tuple()

    b_coeffs_final =
      all_b_coeffs
      |> Enum.take(@n_terms)
      |> List.to_tuple()

    # R_G/a or R4/a (Radius of rectifying sphere / semi-major axis)
    # NGA.SIG.0012_2.0.0_UTMUPS eq 3-30: (1/(1+n)) * (1 + 1/4 n^2 + 1/64 n^4 + 1/256 n^6 + 25/16384 n^8 + ...)
    # The C++ code uses n10 as well.
    r4oa_sum_terms =
      1.0 + n2 / 4.0 + n4 / 64.0 + n6 / 256.0 + 25.0 * n8 / 16384.0 + 49.0 * n10 / 65536.0

    r4oa = r4oa_sum_terms / (1.0 + n1)

    {a_coeffs_final, b_coeffs_final, r4oa}
  end

  # Public accessor for the module attribute @n_terms if needed by other modules (e.g., tests)
  def n_terms(), do: @n_terms
end
