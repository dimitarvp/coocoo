defmodule CooCoo.CoreHelpersTest do
  use ExUnit.Case, async: true

  alias CooCoo.Constants
  alias CooCoo.MathHelpers
  alias CooCoo.Projections.TMHelpers
  alias CooCoo.Projections.TransverseMercator
  alias CooCoo.Projections.PolarStereographic

  import CooCoo.Test.AssertHelpers

  @high_precision_delta 1.0e-14
  @moderate_precision_delta 1.0e-9

  describe "CooCoo.Constants" do
    test "mathematical constants are correct" do
      assert_floats_close(Constants.pi(), :math.pi(), @high_precision_delta, "Pi constant")

      assert_floats_close(
        Constants.pi_over_2(),
        :math.pi() / 2.0,
        @high_precision_delta,
        "Pi_over_2"
      )

      assert_floats_close(
        Constants.pi_over_4(),
        :math.pi() / 4.0,
        @high_precision_delta,
        "Pi_over_4"
      )

      assert_floats_close(
        Constants.pi_over_180(),
        :math.pi() / 180.0,
        @high_precision_delta,
        "Pi_over_180"
      )

      assert_floats_close(Constants.two_pi(), 2.0 * :math.pi(), @high_precision_delta, "Two_pi")
    end

    test "WGS84 ellipsoid parameters are defined and consistent" do
      a = Constants.wgs84_a()
      f = Constants.wgs84_f()
      inv_f = Constants.wgs84_inv_f()
      es = Constants.wgs84_es()
      es2 = Constants.wgs84_es2()

      assert a == 6_378_137.0, "WGS84 Semi-major axis"
      assert inv_f == 298.257223563, "WGS84 Inverse Flattening"
      assert_floats_close(f, 1.0 / inv_f, @high_precision_delta, "WGS84 Flattening derivation")

      assert_floats_close(
        es2,
        2.0 * f - f * f,
        @high_precision_delta,
        "WGS84 EccentricitySquared derivation"
      )

      assert_floats_close(
        es,
        :math.sqrt(es2),
        @high_precision_delta,
        "WGS84 Eccentricity derivation"
      )
    end

    test "degree/radian conversion helpers work" do
      assert_floats_close(
        Constants.degrees_to_radians(180.0),
        Constants.pi(),
        @high_precision_delta
      )

      assert_floats_close(
        Constants.radians_to_degrees(Constants.pi()),
        180.0,
        @high_precision_delta
      )

      assert_floats_close(
        Constants.degrees_to_radians(90.0),
        Constants.pi_over_2(),
        @high_precision_delta
      )

      assert_floats_close(
        Constants.radians_to_degrees(Constants.pi_over_2()),
        90.0,
        @high_precision_delta
      )

      assert_floats_close(Constants.degrees_to_radians(0.0), 0.0, @high_precision_delta)
    end

    test "MGRS specific numerical constants are correct" do
      assert Constants.one_hundred_thousand() == 100_000.0
      assert Constants.two_million() == 2_000_000.0
      assert Constants.max_mgrs_precision() == 5

      assert_floats_close(
        Constants.min_mgrs_non_polar_lat(),
        -80.0 * Constants.pi_over_180(),
        @high_precision_delta
      )

      assert_floats_close(
        Constants.max_mgrs_non_polar_lat(),
        84.0 * Constants.pi_over_180(),
        @high_precision_delta
      )
    end
  end

  describe "CooCoo.MathHelpers" do
    test "equal?/3 for floats" do
      assert MathHelpers.equal?(1.0, 1.00000000001, 1.0e-9)
      refute MathHelpers.equal?(1.0, 1.000001, 1.0e-9)
      assert MathHelpers.equal?(1.23456789, 1.23456788, 1.0e-7)
      assert MathHelpers.equal?(-5.0, -5.000000000000001, @high_precision_delta)
    end

    test "zero?/2 for floats" do
      assert MathHelpers.zero?(0.0)
      assert MathHelpers.zero?(1.0e-13)
      refute MathHelpers.zero?(1.0e-11)
      assert MathHelpers.zero?(-1.0e-14, 1.0e-13)
    end

    test "finite_float?/1" do
      assert MathHelpers.finite_float?(0.0)
      assert MathHelpers.finite_float?(123.456)
      assert MathHelpers.finite_float?(-1.0e300)
      refute MathHelpers.finite_float?(:infinity)
      refute MathHelpers.finite_float?(:"-infinity")
      # Elixir's `0.0/0.0` raises ArithmeticError, so we can't pass it directly.
      # Test with a known NaN producing operation if available from another source,
      # or trust that `value == value` handles actual NaN float values.
      # For now, ensuring it doesn't crash with standard floats and infinities.
      refute MathHelpers.finite_float?("not a float")
      refute MathHelpers.finite_float?(nil)
      # Test NaN if we can construct one safely, e.g., if a function is known to return it.
      # For now, assume `value == value` correctly handles actual NaN bit patterns.
      # nan = :math.acos(2.0) # This would produce NaN
      # refute MathHelpers.finite_float?(nan) # Currently :math.acos(2.0) errors.
    end
  end

  describe "CooCoo.Projections.TMHelpers.generate_coefficients/1" do
    test "calculate_beta_1/8 generates correct value for WGS84 n" do
      f_wgs84 = Constants.wgs84_f()
      inv_f_wgs84 = 1.0 / f_wgs84
      n1 = 1.0 / (2.0 * inv_f_wgs84 - 1.0)
      n2 = n1 * n1
      n3 = n2 * n1
      n4 = n3 * n1
      n5 = n4 * n1
      n6 = n5 * n1
      n7 = n6 * n1
      n8 = n7 * n1

      # This is the formula for beta_1 from TMHelpers.calculate_beta_1
      # beta_1 = 0.5*n1 - (2/3)*n2 + (37/96)*n3 - (1/360)*n4 - (81/512)*n5 + (96199/604800)*n6 - (5406467/38707200)*n7 + (7944359/67737600)*n8
      # The expected value from the C++ hardcoded "WE" coefficients for the first beta term (bCoeff[0])
      # which corresponds to our beta_1 in the series.
      expected_beta1_wgs84 = -8.3773216405794867707e-04

      # Note: TMHelpers.calculate_beta_1 is not exported by default from TMHelpers if it's defp.
      # You confirmed you made them public (def).
      actual_beta1 = TMHelpers.calculate_beta_1(n1, n2, n3, n4, n5, n6, n7, n8)

      assert_floats_close(
        actual_beta1,
        expected_beta1_wgs84,
        # Using a very small delta for this direct formula check
        @high_precision_delta
      )
    end

    test "generates WGS84 coefficients correctly" do
      # Access @n_terms from the module where it's defined for this context
      n_terms = CooCoo.Projections.TMHelpers.n_terms()
      {a_coeffs, b_coeffs, r4oa} = TMHelpers.generate_coefficients(Constants.wgs84_f())

      expected_a1_wgs84 = 8.3773182062446983032e-04
      expected_b1_wgs84 = -8.3773216405794867707e-04
      expected_r4oa_wgs84 = 0.999600000383918

      assert tuple_size(a_coeffs) == n_terms
      assert tuple_size(b_coeffs) == n_terms

      assert_floats_close(
        elem(a_coeffs, 0),
        expected_a1_wgs84,
        @high_precision_delta,
        "TM alpha1 coeff for WGS84"
      )

      assert_floats_close(
        elem(b_coeffs, 0),
        expected_b1_wgs84,
        @high_precision_delta
      )

      assert_floats_close(r4oa, expected_r4oa_wgs84, @high_precision_delta)
    end
  end

  describe "CooCoo.Projections.TransverseMercator public helpers" do
    # Ensure atanh/1, geodetic_latitude_from_conformal/2 etc. are `def` in TransverseMercator.ex
    test "atanh/1 works correctly" do
      assert TransverseMercator.atanh(0.0) == {:ok, 0.0}
      # The value 0.5493061443340549 is atanh(0.5)
      {:ok, atanh_0_5} = TransverseMercator.atanh(0.5)
      assert_floats_close(atanh_0_5, 0.5493061443340549, @high_precision_delta)
      assert TransverseMercator.atanh(1.0) == {:ok, :infinity}
      assert TransverseMercator.atanh(-1.0) == {:ok, :"-infinity"}
      assert TransverseMercator.atanh(2.0) == {:error, :atanh_domain_error_abs_x_gt_1}
      assert TransverseMercator.atanh(-2.0) == {:error, :atanh_domain_error_abs_x_gt_1}
    end

    test "geodetic_latitude_from_conformal/2 calculates correctly" do
      es_wgs84 = Constants.wgs84_es()
      sin_chi_val = 0.5

      # Expected value from external trusted source (e.g., Snyder or other geodetic libraries)
      # For WGS84, es = 0.0818191908426, sin(conformal_lat) = 0.5
      # Resulting geodetic latitude approx 0.522704908 radians (29.949901 degrees)
      # Increased precision for expected value
      expected_phi_rad = 0.5227049080181342

      case TransverseMercator.geodetic_latitude_from_conformal(sin_chi_val, es_wgs84) do
        {:ok, actual_phi_rad} ->
          assert_floats_close(actual_phi_rad, expected_phi_rad, @moderate_precision_delta)

        {:error, reason} ->
          flunk("geodetic_latitude_from_conformal failed: #{inspect(reason)}")
      end

      case TransverseMercator.geodetic_latitude_from_conformal(0.0, es_wgs84) do
        {:ok, actual_phi_rad} ->
          assert_floats_close(actual_phi_rad, 0.0, @high_precision_delta)

        {:error, reason} ->
          flunk("geodetic_latitude_from_conformal failed for sin_chi=0: #{inspect(reason)}")
      end
    end

    test "compute_hyperbolic_series/2 and compute_trig_series/2 return correct tuple sizes" do
      # Access @n_terms from where it is defined
      n_terms_tm = CooCoo.Projections.TransverseMercator.n_terms()

      assert elem(TransverseMercator.compute_hyperbolic_series(0.5, n_terms_tm), 0) |> tuple_size ==
               n_terms_tm

      assert elem(TransverseMercator.compute_hyperbolic_series(0.5, n_terms_tm), 1) |> tuple_size ==
               n_terms_tm

      assert elem(TransverseMercator.compute_trig_series(0.5, n_terms_tm), 0) |> tuple_size ==
               n_terms_tm

      assert elem(TransverseMercator.compute_trig_series(0.5, n_terms_tm), 1) |> tuple_size ==
               n_terms_tm
    end
  end

  describe "CooCoo.Projections.PolarStereographic public helpers" do
    # Ensure polar_pow/2 is `def` in PolarStereographic.ex
    test "polar_pow/2 calculates correctly" do
      es_wgs84 = Constants.wgs84_es()
      sin_60_deg = :math.sin(60.0 * Constants.pi_over_180())

      expected_polar_pow_val = 0.994209594256263

      case PolarStereographic.polar_pow(es_wgs84, sin_60_deg) do
        {:ok, actual_value} ->
          assert_floats_close(
            actual_value,
            expected_polar_pow_val,
            @high_precision_delta
          )

        {:error, reason} ->
          flunk("polar_pow/2 returned an error for valid inputs: #{inspect(reason)}")
      end

      # Test edge case: sin_latitude = 0 (equator)
      case PolarStereographic.polar_pow(es_wgs84, 0.0) do
        {:ok, actual_value} ->
          assert_floats_close(actual_value, 1.0, @high_precision_delta)

        {:error, reason} ->
          flunk("polar_pow/2 returned an error for sin_latitude = 0.0: #{inspect(reason)}")
      end

      # Test for an error case if polar_pow is designed to return one for certain inputs
      # For example, if sin_latitude > 1.0 (which shouldn't happen with :math.sin)
      # or if eccentricity was such that es_sin >= 1
      # This would be an invalid input for sin_latitude
      invalid_sin_lat = 1.1
      # polar_pow now has internal validation for this
      assert PolarStereographic.polar_pow(es_wgs84, invalid_sin_lat) ==
               {:error, :invalid_input_to_polar_pow}

      # Test with es_sin very close to 1 (eccentricity = 0.9999999999, sin_latitude = 1.0)
      # This should produce a very small positive number, not an error,
      # based on direct C++ equivalent math.
      expected_val_near_boundary = 7.071068112959566e-6

      case PolarStereographic.polar_pow(0.9999999999, 1.0) do
        {:ok, actual_val} ->
          assert_floats_close(actual_val, expected_val_near_boundary, @high_precision_delta)

        {:error, reason} ->
          flunk("polar_pow(0.9999999999, 1.0) returned an unexpected error: #{inspect(reason)}")
      end
    end
  end
end
