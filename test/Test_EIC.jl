using JuTrack
using Test
using BenchmarkTools
using Plots
using DelimitedFiles
using Serialization
using StructArrays
using FFTW
using LaTeXStrings
using ProgressMeter
using Distributions

theme(:ggplot2)
#theme(:vibrant)

struct record
    cx::Float64
    cy::Float64
    cz::Float64
    sx::Float64
    sy::Float64
    sz::Float64
    ex::Float64
    ey::Float64
    ez::Float64
    lumi::Float64
end


begin
    num_particles = 1000
    turns = 100

    pbeam = Beam(zeros(1000, 6), np=Int(0.688e11), energy=275e9, emittance = [11.3e-9, 1e-9, 3.7e-2])  #[11.3e-9, 1e-9, 3.7e-2] #3.7e-2 -> 6cm
    opIPp = optics4DUC(0.8, 0.0, 0.072, 0.0)
    mainRF = AccelCavity(freq=591e6, volt=15.8e6, h=7560.0, phis=Float64(pi))

    αc=1.5e-3
    lmap = LongitudinalRFMap(αc, mainRF)
    initilize_6DGaussiandist!(pbeam, opIPp, lmap)

    # define oneturnmap
    tunex=0.228
    tuney=0.21
    oneturn = TransferMap4DChrom(opIPp, tunex, tuney, 1.0, 1.0)
    
    #estrong beam
    opIPe = optics4DUC(0.45, 0.0, 0.056, 0.0)
    estrong = StrongGaussianBeam(1.0, m_e, 1.0, Int(1.72e11), 10e9,  opIPe, [95e-6, 8.5e-6, 0.007], 5)
    initilize_zslice!(estrong, :gaussian, :evennpar, 5.0)


    # define crab cavity
    crab_ratio=0.33
    overcrab=1.0

    pcrabu = easyCRABCAVITY(freq=197.0e6, halfthetac=overcrab*12.5e-3*(1+crab_ratio))
    pcrabu2nd = easyCRABCAVITY(freq=197.0e6*2.0, halfthetac=-overcrab*12.5e-3*crab_ratio)
    pcrabd = easyCRABCAVITY(freq=197.0e6, halfthetac=-overcrab*12.5e-3*(1+crab_ratio))
    pcrabd2nd = easyCRABCAVITY(freq=197.0e6*2.0, halfthetac=overcrab*12.5e-3*crab_ratio)


    """
    pcrab1st = easyCRABCAVITY(freq=197.0e6, halfthetac=overcrab*12.5e-3*(1+crab_ratio))
    pcrab2nd = easyCRABCAVITY(freq=197.0e6*2.0, halfthetac=-overcrab*12.5e-3*crab_ratio)
    crab_crossing_setup!(pstrong, 12.5e-3, pcrab1st, pcrab2nd) #ma pass! è x track

    """

    # define Lorentz boost
    lb = LorentzBoost(12.5e-3)
    invlb = InvLorentzBoost(12.5e-3)
 
    #define SC_lens for SC_kick BE
    nSC = 4
    ds = 3800 / (nSC + 1)
    opSC = optics4DUC(1.0, 0.0, 1.0, 0.0)
    sc_BE = SC_lens(opSC, ds, nSC)
    turns_rec = 0
    w = 0.0

    #phi_advx = LinRange(0, 29.228, nSC+1)
    #phi_advy = LinRange(0, 30.210, nSC+1)
    TM1 = Array{TransferMap4D, 1}(undef, nSC)
    TM1_inv = Array{Inverse_TransferMap4D, 1}(undef, nSC)
    
    # records
    particles_turns = zeros(Float64, num_particles, 6, turns)
    particles_turns_SC = zeros(Float64, num_particles, 6, turns)
    delta_px = zeros(Float64, num_particles)
    records = StructArray{record}(undef, turns)

    bb = true
    lumi=0.0

    """
    TM1 = TransferMap4D(opIPp, opSC, phi_advx[2], phi_advy[2])
    TM1_inv = Inverse_TransferMap4D(opIPp, opSC, phi_advx[2], phi_advy[2])
    """

    ring = [oneturn, sc_BE]
    #ring = [oneturn, TM1, sc_BE, TM1_inv, TM1, sc_BE, TM1_inv, TM1, sc_BE, TM1_inv, TM1, sc_BE, TM1_inv, TM1, sc_BE, TM1_inv]
    #ring = [oneturn, pcrabu, pcrabu2nd, lb, estrong, invlb, pcrabd, pcrabd2nd] #se TransferMap4D dentro SC function

end


begin
    N = 1000
    new_emit1 = zeros(N+1, 3)
    NLOST = zeros(N)
    println("Start tracking")
    prog = Progress(N)
    X1 = zeros(N, 4)
    for i in 1:N
        # println("Turn: ", i)
        if i == 1
            get_emittance!(pbeam)
          new_emit1[i, :] = pbeam.emittance
        end

        linepass!(ring, pbeam)
        NLOST[i] = sum(pbeam.lost_flag)
        get_2nd_moment!(pbeam)
        get_emittance!(pbeam)
        new_emit1[i+1, :] = pbeam.emittance
        # save
        if i % 100 == 0
            writedlm("emit1_SCBE_2.txt", new_emit1)
            writedlm("lost_SCBE_2.txt", NLOST)
        end
        X1[i, :] = pbeam.r[1, 1:4] #track first particle
        next!(prog)
    end
end


get_emittance!(pbeam)
println(pbeam.emittance)

plot(new_emit1[:,1])
plot(X1[:, 1], X1[:, 2])
plot(pbeam.r[:,1], pbeam.r[:,2], seriestype=:scatter, markersize=0.9, xlabel="x", ylabel= "px")
plot(pbeam.r[:,3], pbeam.r[:,4], seriestype=:scatter, markersize=0.9, xlabel="x", ylabel= "px")

