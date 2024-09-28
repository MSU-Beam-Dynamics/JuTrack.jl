using JuTrack

function f1(x)
    vbase=3.42*8.5e6
    ϕs=10.0
    vact=vbase/cos(ϕs*π/180.0)
    mainRFe=AccelCavity(freq=591e6, volt=vact, h=7560.0, phis=π-ϕs*π/180.0)
    tunex, tuney=50.08, 44.14
    αc=3.42/tunex/tunex
    lmap=LongitudinalRFMap(αc, mainRFe)
    opIPp=optics4DUC(0.8, 0.0, 0.072, 0.0)
    opIPe=optics4DUC(0.45,0.0,0.056,0.0)
    pstrong=StrongGaussianBeam(1.0, m_p, 1.0, Int(0.688e11), 275e9,  opIPp, [x, 8.5e-6, 0.06], 9)
    initilize_zslice!(pstrong, :gaussian, :evennpar, 7.0)
    crab_ratio=0.33
    overcrab=1.0
    pcrab1st = easyCRABCAVITY(freq=197.0e6, halfthetac=overcrab*12.5e-3*(1+crab_ratio))
    pcrab2nd = easyCRABCAVITY(freq=197.0e6*2.0, halfthetac=-overcrab*12.5e-3*crab_ratio)
    crab_crossing_setup!(pstrong, 12.5e-3, pcrab1st, pcrab2nd)
    ebeam=Beam(zeros(5000, 6), np = Int(1.72e11*3), energy = 10e9, emittance=[20e-9, 1.3e-9, 1.36e-4])
    initilize_6DGaussiandist!(ebeam, opIPe, lmap)
    linepass!([pstrong], ebeam)
    get_emittance!(ebeam)
    return ebeam.emittance[1]
end

println(f1(95e-6))
g = autodiff(Forward, f1, Duplicated(95e-6, 1.0))
print("grad: ", g)

