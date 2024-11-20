using Documenter
using JuTrack

makedocs(
    sitename = "JuTrack.jl",
    modules = [JuTrack],
    format = Documenter.HTML(),
    pages = [
        "Home" => "index.md",
        "API" => "api.md"
    ]
)

deploydocs(
    repo = "github.com/MSU-Beam-Dynamics/JuTrack.jl.git",
    target = "gh-pages",
)