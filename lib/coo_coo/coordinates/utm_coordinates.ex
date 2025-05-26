defmodule CooCoo.Coordinates.UTMCoordinates do
  @moduledoc """
  Represents Universal Transverse Mercator (UTM) coordinates.
  """
  @typedoc """
  - `:zone` - UTM zone number (1-60).
  - `:hemisphere` - Hemisphere ('N' or 'S').
  - `:easting` - Easting value in meters.
  - `:northing` - Northing value in meters.
  - `:warning_message` - Any warning message associated with the coordinate.
  - `:error_message` - Any error message associated with the coordinate.
  """
  @enforce_keys [:zone, :hemisphere, :easting, :northing]
  defstruct zone: 0,
            # Using char ?N and ?S
            hemisphere: ?N,
            easting: 0.0,
            northing: 0.0,
            warning_message: "",
            error_message: ""

  @type t :: %__MODULE__{
          zone: integer,
          hemisphere: char,
          easting: float,
          northing: float,
          warning_message: String.t(),
          error_message: String.t()
        }
end
