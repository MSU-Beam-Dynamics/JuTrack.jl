using Pkg
Pkg.activate(".")                      # change "." to your path of JuTrack.jl
Pkg.add(["SnoopCompile"])              # SnoopCompile usually re-exports the macro
try
    Pkg.add(["SnoopCompileCore"])      # present on some versions
catch
end

using JuTrack
using SnoopCompile

tinf = SnoopCompile.@snoopi_deep begin
    include(joinpath(@__DIR__, "_agg_warmup.jl"))
end

pc  = SnoopCompile.parcel(tinf)
mkpath(@__DIR__)
out = joinpath(@__DIR__, "precompile_statements.jl")
SnoopCompile.write(out, pc)
println("Wrote precompile statements to: ", out)
