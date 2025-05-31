defmodule CooCoo.MathHelpers do
  @moduledoc """
  Provides helper functions for floating-point comparisons with epsilon tolerance.
  """

  # A general-purpose epsilon
  @default_epsilon 1.0e-12

  @doc """
  Checks if two floats are equal within a given epsilon.
  Default epsilon is #{@default_epsilon}.
  """
  def equal?(a, b, epsilon \\ @default_epsilon)
      when is_number(a) and is_number(b) and is_number(epsilon) do
    abs(a - b) < epsilon
  end

  @doc """
  Checks if a float is approximately zero within a given epsilon.
  Default epsilon is #{@default_epsilon}.
  """
  def zero?(a, epsilon \\ @default_epsilon) when is_number(a) and is_number(epsilon) do
    abs(a) < epsilon
  end
end
