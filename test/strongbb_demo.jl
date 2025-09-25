# This code simulate strong bb interaction and calculate the final emittance of the electron beam
# The derivatives of the final results w.r.t. the emittance of the proton beam at IP are automatically calculated via AD
# This script is just a test for the AD feature in JuTrack. The results may not be physically meaningful.
using Pkg
Pkg.activate("."); Pkg.instantiate() # change "." to your path of JuTrack.jl
using JuTrack
set_tps_dim(3)
function f1(x1::DTPSAD{N, T}, x2::DTPSAD{N, T}, x3::DTPSAD{N, T}) where {N, T}
    vbase=3.42*8.5e6
    ϕs=10.0
    vact=vbase/cos(ϕs*π/180.0)
    mainRFe=AccelCavity(freq=DTPSAD(591e6), volt=vact, h=7560.0, phis=π-ϕs*π/180.0)
    tunex, tuney=50.08, 44.14
    alphac=3.42/tunex/tunex

    lmap=LongitudinalRFMap(DTPSAD(alphac), mainRFe)
    opIPp=optics4DUC(0.8, 0.0, 0.072, 0.0)
    opIPe=optics4DUC(0.45,0.0,0.056,0.0)
    pstrong=StrongGaussianBeam(DTPSAD(1.0), DTPSAD(m_p), DTPSAD(1.0), Int(0.688e11), DTPSAD(275e9), opIPp, [x1, x2, x3], 9)
    initilize_zslice!(pstrong, :gaussian, :evennpar, 7.0)
    crab_ratio=0.33
    overcrab=1.0
    pcrab1st = easyCRABCAVITY(freq=197.0e6, halfthetac=overcrab*12.5e-3*(1+crab_ratio))
    pcrab2nd = easyCRABCAVITY(freq=197.0e6*2.0, halfthetac=-overcrab*12.5e-3*crab_ratio)
    crab_crossing_setup!(pstrong, 12.5e-3, pcrab1st, pcrab2nd)
    ebeam=Beam(zeros(DTPSAD{NVAR(), Float64}, 5000, 6), np = Int(1.72e11*3), energy = 10e9, emittance=[20e-9, 1.3e-9, 1.36e-4])
    initilize_6DGaussiandist!(ebeam, opIPe, lmap)
    linepass!([pstrong], ebeam)
    get_emittance!(ebeam)
    return ebeam.emittance
end

x1 = DTPSAD(95e-6, 1)
x2 = DTPSAD(8.5e-6, 2)
x3 = DTPSAD(0.06, 3)

println("x emittance: ", f1(x1, x2, x3)[1], "\ny emittance: ", f1(x1, x2, x3)[2], "\nz emittance: ", f1(x1, x2, x3)[3])
println("Jacobian: ", Jacobian(f1, [95e-6, 8.5e-6, 0.06])) # or Jacobian(f1, [x1, x2, x3])


###################################################################################################
# Legacy code of Using Enzyme to calculate the gradient of a strong beam-beam interaction function
# Uncomment the following lines to try Enzyme for automatic differentiation
###################################################################################################
# using JuTrack

# function f1(x)
#     vbase=3.42*8.5e6
#     ϕs=10.0
#     vact=vbase/cos(ϕs*π/180.0)
#     mainRFe=AccelCavity(freq=591e6, volt=vact, h=7560.0, phis=π-ϕs*π/180.0)
#     tunex, tuney=50.08, 44.14
#     alphac=3.42/tunex/tunex
#     lmap=LongitudinalRFMap(alphac, mainRFe)
#     opIPp=optics4DUC(0.8, 0.0, 0.072, 0.0)
#     opIPe=optics4DUC(0.45,0.0,0.056,0.0)
#     pstrong=StrongGaussianBeam(1.0, m_p, 1.0, Int(0.688e11), 275e9,  opIPp, [x, 8.5e-6, 0.06], 9)
#     initilize_zslice!(pstrong, :gaussian, :evennpar, 7.0)
#     crab_ratio=0.33
#     overcrab=1.0
#     pcrab1st = easyCRABCAVITY(freq=197.0e6, halfthetac=overcrab*12.5e-3*(1+crab_ratio))
#     pcrab2nd = easyCRABCAVITY(freq=197.0e6*2.0, halfthetac=-overcrab*12.5e-3*crab_ratio)
#     crab_crossing_setup!(pstrong, 12.5e-3, pcrab1st, pcrab2nd)
#     ebeam=Beam(zeros(5000, 6), np = Int(1.72e11*3), energy = 10e9, emittance=[20e-9, 1.3e-9, 1.36e-4])
#     initilize_6DGaussiandist!(ebeam, opIPe, lmap)
#     linepass!([pstrong], ebeam)
#     get_emittance!(ebeam)
#     return ebeam.emittance[1]
# end

# println(f1(95e-6))
# g = autodiff(ForwardWithPrimal, f1, Duplicated(95e-6, 1.0))
# print("grad: ", g)

