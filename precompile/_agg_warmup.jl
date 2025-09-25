# keep data tiny and avoid external IO.

using Pkg
Pkg.activate("."); Pkg.instantiate() # change "." to your path of JuTrack.jl
using JuTrack

# Add the function calls you want Julia to precompile in this file, and remove any you don't need to reduce the size of the sysimage.
include("../test/alsu1.jl")
include("../test/phase_advance_tuning_demo.jl")
include("../test/jacobian_benchmark.jl")
include("../test/RDT_demo.jl")
include("../test/strongbb_demo.jl")
include("../test/wakefield_demo.jl")
include("../src/demo/ESR/esr-main.jl")
include("../src/demo/spear3/spear3.jl")
include("../src/demo/transverse_space_charge/emit_gradient_withSC.jl")
include("../src/demo/TwissTuning/TwissTuning_example.jl")
