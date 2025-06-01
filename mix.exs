defmodule CooCoo.MixProject do
  use Mix.Project

  def project do
    [
      app: :coocoo,
      version: "0.1.0",
      elixir: "~> 1.18",
      start_permanent: Mix.env() == :prod,
      deps: deps(),
      elixirc_paths: elixirc_paths(Mix.env()),
      dialyzer: dialyzer(Mix.env())
    ]
  end

  def application do
    [
      extra_applications: [:logger]
    ]
  end

  defp deps do
    [
      {:decimal, "~> 2.3"},
      {:dialyxir, "~> 1.4", only: :dev, runtime: false}
    ]
  end

  defp elixirc_paths(:test), do: ["lib", "test/support"]
  defp elixirc_paths(_), do: ["lib"]

  defp dialyzer(_env) do
    [
      # Specifies the directory where core PLTs (OTP, Elixir stdlib) are stored.
      plt_core_path: "priv/plts/",
      # Specifies the path to the final project PLT file, which includes dependencies.
      # Using {:no_warn, ...} suppresses warnings if the file doesn't exist initially.
      plt_file: {:no_warn, "priv/plts/core.plt"}
      # You might want to add other options here, for example:
      # applications: [:credo, :another_dep], # To add specific apps to the PLT
      # flags: ["-Wunmatched_returns", ...], # Custom Dialyzer flags
      # For most projects, the defaults for 'applications' (includes project + deps) are fine.
    ]
  end
end
