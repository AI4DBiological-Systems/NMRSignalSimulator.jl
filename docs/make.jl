using Documenter
using NMRSignalSimulator

makedocs(
    sitename = "NMRSignalSimulator",
    format = Documenter.HTML(),
    modules = [NMRSignalSimulator]
)

makedocs(
    sitename="NMRSignalSimulator.jl",
    modules=[NMRSignalSimulator],
    format=Documenter.HTML(prettyurls = get(ENV, "CI", nothing)=="true"),
    pages=[
        "Overview" => "index.md",
        "Public API" => "api.md",
        "Demo: code walk-through" => "demo_code.md",
    ],
)

deploydocs(
    repo = "github.com/AI4DBiological-Systems/NMRSignalSimulator.jl.git",
    target = "build",
    branch = "gh-pages",
    versions = ["stable" => "v^", "v#.#" ],
)