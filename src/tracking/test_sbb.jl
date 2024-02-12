
include("../JuTrack.jl")
using. JuTrack
using Enzyme
Enzyme.API.runtimeActivity!(true)
function f(x, ebeam)

    opIPp=optics4DUC(0.8, 0.0, 0.072, 0.0)
    pstrong=StrongGaussianBeam(1.0, m_p, 1.0, 68800000000, x,  opIPp, [95e-6, 8.5e-6, 0.06], 11)
    initilize_zslice!(pstrong, :gaussian, :evennpar, 5.0)
    linepass!([pstrong], ebeam)
    # get_emittance!(ebeam)
    return ebeam.emittance[1]
end

ebeam = Beam(1e9, 172000000000, 100, charge=-1.0, emittance=[2e-8, 1.3e-9, 1.36e-4])
opIP=optics4DUC(0.45,0.0,0.056,0.0)

vbase = 3.42*8.5e6
phis = 10.0
vact = vbase/cos(phis*π/180)
mainRF=AccelCavity(591e6, vact, 7560.0, π-phis*π/180.0)
tunex, tuney=50.08, 44.14
alphac = 3.42/tunex/tuney
lmap = LongitudinalRFMap(alphac, mainRF)
initilize_6DGaussiandist!(ebeam, opIP, lmap, 1.8e-3)
get_emittance!(ebeam)

println(f(275e9, ebeam))
grad = autodiff(Forward, f, Duplicated(275e9, 1.0), Const(ebeam))
# not successful yet. Because of the interaction between SpecialFunctions and Complex numbers.
println(grad)
# scatter(ebeam.r[:,6], ebeam.r[:,1]; markersize=0.1, color=:oslo, xlabel="z [m]", ylabel="x [m]")

