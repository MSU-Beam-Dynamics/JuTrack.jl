include("../src/JuTrack.jl")
using .JuTrack
using Plots  
using LaTeXStrings
using Enzyme
Enzyme.API.runtimeActivity!(true)



function f(x1) # treat beam as a constant will lead to 0 gradient
    # particles = randn(50000, 6) / 1e3
    # ebeam = Beam(particles, znbin=2)
    RLCwake = LongitudinalRLCWake(x1, 5.5e3, 3.0)

    vbase=3.42*8.5e6
    ϕs=10.0
    vact=vbase/cos(ϕs*π/180.0)
    mainRFe=AccelCavity(591e6, vact, 7560.0, π-ϕs*π/180.0)
    tunex, tuney=50.08, 44.14
    alphac=3.42/tunex/tunex
    lmap=LongitudinalRFMap(alphac, mainRFe)

    particles = randn(50000, 6) / 1e3
    ebeam = Beam(particles, znbin=10)
    histogram1DinZ!(ebeam)

    linepass!([lmap, RLCwake], ebeam)
    get_emittance!(ebeam)

    return ebeam.emittance
end

grad = autodiff(Forward, f, Duplicated,  Duplicated(180e9, 1.0)) # 11.713 s (1750155 allocations: 74.56 GiB)
println(grad)

# @btime begin
#     grad = autodiff(Forward, fp, Duplicated,  Duplicated(180e9, 1.0))
# end
