using Pkg
Pkg.activate("."); Pkg.instantiate() # change "." to your path of JuTrack.jl
Pkg.add("PackageCompiler")
using PackageCompiler

ext = Sys.iswindows() ? ".dll" : (Sys.isapple() ? ".dylib" : ".so")
img = joinpath(@__DIR__, "JuTrack_sysimage" * ext)

warmup = joinpath(@__DIR__, "_agg_warmup.jl")

traces = filter(isfile, [
    joinpath(@__DIR__, "trace_base.jl"),
    joinpath(@__DIR__, "trace_core.jl"),
    joinpath(@__DIR__, "trace_optics.jl"),
    joinpath(@__DIR__, "trace_tpsa.jl"),
])

kwargs = (
    sysimage_path = img,
    incremental = false,
    cpu_target = get(ENV, "JUTRACK_CPU", "generic"),
)

opt = kwargs
if isfile(warmup)
    opt = merge(opt, (precompile_execution_file = warmup,))
end
if !isempty(traces)
    opt = merge(opt, (precompile_statements_file = traces,))
end

create_sysimage([:JuTrack]; opt...)

println("Sysimage built at: ", img)
