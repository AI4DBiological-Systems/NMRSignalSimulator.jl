using Documenter
using NMRSignalSimulator

makedocs(
    sitename = "NMRSignalSimulator",
    format = Documenter.HTML(),
    pages=[
        "Overview" => "index.md",
        "Public API" => "api.md",
        "Demo: resonance groups" => "generated/glutamine.md",
        #"Demo: frequency-domain surrogate" => "generated/surrogate.md",

    ],
)

makedocs(
    sitename="NMRSignalSimulator.jl",
    modules=[NMRSignalSimulator],
    format=Documenter.HTML(prettyurls = get(ENV, "CI", nothing)=="true"),
    pages=[
        "Overview" => "index.md",
        "Public API" => "api.md",
        "Demo: resonance groups" => "generated/glutamine.md",
    ],
)

deploydocs(
    repo = "github.com/AI4DBiological-Systems/NMRSignalSimulator.jl.git",
    target = "build",
    branch = "gh-pages",
    versions = ["stable" => "v^", "v#.#" ],
)