defmodule CooCoo.Coordinates.UPSCoordinates do
  @moduledoc """
  Represents Universal Polar Stereographic (UPS) coordinates.
  """
  @typedoc """
  - `:hemisphere` - Hemisphere ('N' or 'S').
  - `:easting` - Easting value in meters.
  - `:northing` - Northing value in meters.
  - `:warning_message` - Any warning message associated with the coordinate.
  - `:error_message` - Any error message associated with the coordinate.
  """
  @enforce_keys [:hemisphere, :easting, :northing]
  # Using char ?N and ?S
  defstruct hemisphere: ?N,
            easting: 0.0,
            northing: 0.0,
            warning_message: "",
            error_message: ""

  @type t :: %__MODULE__{
          hemisphere: char,
          easting: float,
          northing: float,
          warning_message: String.t(),
          error_message: String.t()
        }
end
