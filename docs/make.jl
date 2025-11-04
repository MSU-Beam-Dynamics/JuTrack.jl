using Pkg
Pkg.develop(path=joinpath(@__DIR__, ".."))
Pkg.instantiate()

using Documenter
using JuTrack

makedocs(
    sitename = "JuTrack.jl",
    modules = [JuTrack],
    format = Documenter.HTML(),
    warnonly = [:missing_docs, :cross_references],
    pages = [
        "Home" => "index.md",
        "API"  => "api.md",
    ],
)

deploydocs(
    repo      = "github.com/MSU-Beam-Dynamics/JuTrack.jl.git",
    devbranch = "main",
    versions  = ["stable" => "v^", "dev" => "main"],
)
