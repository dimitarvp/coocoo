defmodule CooCoo.Coordinates.MapProjectionCoordinates do
  @moduledoc """
  Represents generic Map Projection coordinates (easting, northing).
  """
  @typedoc """
  - `:easting` - Easting value in meters.
  - `:northing` - Northing value in meters.
  - `:warning_message` - Any warning message associated with the coordinate.
  - `:error_message` - Any error message associated with the coordinate.
  """
  @enforce_keys [:easting, :northing]
  defstruct easting: 0.0,
            northing: 0.0,
            warning_message: "",
            error_message: ""

  @type t :: %__MODULE__{
          easting: float,
          northing: float,
          warning_message: String.t(),
          error_message: String.t()
        }
end
