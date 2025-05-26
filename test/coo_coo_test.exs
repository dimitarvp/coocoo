defmodule CooCooTest do
  use ExUnit.Case
  doctest CooCoo

  test "greets the world" do
    assert CooCoo.hello() == :world
  end
end
