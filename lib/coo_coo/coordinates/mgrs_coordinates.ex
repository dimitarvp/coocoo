defmodule CooCoo.Coordinates.MGRSCoordinates do
  @moduledoc """
  Represents Military Grid Reference System (MGRS) or U.S. National Grid (USNG) coordinates.
  """
  @typedoc """
  - `:mgrs_string` - The MGRS/USNG coordinate string.
  - `:precision` - Precision level (0-5).
  - `:warning_message` - Any warning message associated with the coordinate.
  - `:error_message` - Any error message associated with the coordinate.
  """
  @enforce_keys [:mgrs_string, :precision]
  defstruct mgrs_string: "",
            # Default to highest precision (1m)
            precision: 5,
            warning_message: "",
            error_message: ""

  @type t :: %__MODULE__{
          mgrs_string: String.t(),
          precision: 0..5,
          warning_message: String.t(),
          error_message: String.t()
        }
end
