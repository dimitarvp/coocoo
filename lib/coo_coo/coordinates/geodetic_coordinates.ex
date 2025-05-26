defmodule CooCoo.Coordinates.GeodeticCoordinates do
  @moduledoc """
  Represents Geodetic coordinates.
  """
  @typedoc """
  - `:longitude` - Longitude in radians.
  - `:latitude` - Latitude in radians.
  - `:height` - Height above the ellipsoid in meters (defaults to 0.0).
  - `:warning_message` - Any warning message associated with the coordinate.
  - `:error_message` - Any error message associated with the coordinate.
  """
  @enforce_keys [:longitude, :latitude]
  defstruct longitude: 0.0,
            latitude: 0.0,
            height: 0.0,
            warning_message: "",
            error_message: ""

  @type t :: %__MODULE__{
          longitude: float,
          latitude: float,
          height: float,
          warning_message: String.t(),
          error_message: String.t()
        }
end
