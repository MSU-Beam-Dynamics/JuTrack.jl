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
    num_particles = 5000
    turns = 10000

    pbeam = Beam(zeros(num_particles, 6), np=Int(0.688e11), energy=275e9, emittance = [11.3e-9, 1e-9, 3.7e-2], mass=938.272e6)  #[11.3e-9, 1e-9, 3.7e-2] #3.7e-2 -> 6cm
    opIPp = optics4DUC(0.8, 0.0, 0.072, 0.0)
    mainRF = AccelCavity(freq=591e6, volt=15.8e6, h=7560.0, phis=Float64(pi))

    αc=1.5e-3
    lmap = LongitudinalRFMap(αc, mainRF)
    initilize_6DGaussiandist!(pbeam, opIPp, lmap)

    # define oneturnmap
    tunex = 0.228
    tuney = 0.210
    oneturn = TransferMap4DChrom(opIPp, tunex, tuney, 1.0, 1.0)
    
    #estrong beam
    opIPe = optics4DUC(0.45, 0.0, 0.056, 0.0)
    estrong = StrongGaussianBeam(1.0, m_e, 1.0, Int(1.72e11), 10e9,  opIPe, [95e-6, 8.5e-6, 0.007], 5)
    initilize_zslice!(estrong, :gaussian, :evennpar, 5.0)

    # define crab cavity
    crab_ratio = 0.33
    overcrab = 1.0

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
    nSC = 5
    ds = 3800 /(nSC + 1)
    opSC = optics4DUC(1.0, 0.0, 1.0, 0.0)
    sc_BE = SC_lens(opSC, ds, nSC)
    w = 0.0

    phi_advx = LinRange(0, 0.228, nSC+1)
    phi_advy = LinRange(0, 0.210, nSC+1)

    #phi_advx = sort(0.228 .* rand(nSC+1))
    #phi_advy = sort(0.210 .* rand(nSC+1))
    
    TM1 = TransferMap4D(opIPp, opSC, phi_advx[1], phi_advy[1])
    TM1_inv = Inverse_TransferMap4D(opIPp, opSC, phi_advx[1], phi_advy[1])

    TM2 = TransferMap4D(opIPp, opSC, phi_advx[2], phi_advy[2])
    TM2_inv = Inverse_TransferMap4D(opIPp, opSC, phi_advx[2], phi_advy[2])

    TM3 = TransferMap4D(opIPp, opSC, phi_advx[3], phi_advy[3])
    TM3_inv = Inverse_TransferMap4D(opIPp, opSC, phi_advx[3], phi_advy[3])
    
    TM4 = TransferMap4D(opIPp, opSC, phi_advx[4], phi_advy[4])
    TM4_inv = Inverse_TransferMap4D(opIPp, opSC, phi_advx[4], phi_advy[4])
    
    TM5 = TransferMap4D(opIPp, opSC, phi_advx[5], phi_advy[5])
    TM5_inv = Inverse_TransferMap4D(opIPp, opSC, phi_advx[5], phi_advy[5])
    

    # records
    particles_turns = zeros(Float64, num_particles, 6, turns)
    particles_turns_SC = zeros(Float64, num_particles, 6, turns)
    delta_px = zeros(Float64, num_particles)
    records = StructArray{record}(undef, turns)

    bb = false
    lumi = 0.0

    ring = [oneturn, TM1, sc_BE, TM1_inv, TM2, sc_BE, TM2_inv, TM3, sc_BE, TM3_inv, TM4, sc_BE, TM4_inv, TM5, sc_BE, TM5_inv] 

    #ring = [oneturn, lmap]
    #ring = [oneturn, TM1, TM1_inv, TM2, TM2_inv, TM3, TM3_inv, TM4, TM4_inv, TM5, TM5_inv]
    #ring = [oneturn, pcrabu, pcrabu2nd, lb, estrong, invlb, pcrabd, pcrabd2nd] #se TransferMap4D dentro SC function

end


get_emittance!(pbeam)
pbeam.emittance
pbeam.beamsize


begin
    N = turns
    new_emit = zeros(N, 3)
    new_bs = zeros(N, 6)
    NLOST = zeros(N)
    println("Start tracking")
    prog = Progress(N)
    X1 = zeros(N, 4)
    
    for i in 1:N

        #linepass!(ring, pbeam)

        linepass!(ring, pbeam)

        particles_turns[:, 1, i] .= pbeam.r[:, 1]
        particles_turns[:, 2, i] .= pbeam.r[:, 2]
        particles_turns[:, 3, i] .= pbeam.r[:, 3]
        particles_turns[:, 4, i] .= pbeam.r[:, 4]
        particles_turns[:, 5, i] .= pbeam.r[:, 5]
        particles_turns[:, 6, i] .= pbeam.r[:, 6]
    

        #linepass!(ring2, pbeam)

        particles_turns_SC[:, 1, i] .= pbeam.r[:, 1]
        particles_turns_SC[:, 2, i] .= pbeam.r[:, 2]
        particles_turns_SC[:, 3, i] .= pbeam.r[:, 3]
        particles_turns_SC[:, 4, i] .= pbeam.r[:, 4]
        particles_turns_SC[:, 5, i] .= pbeam.r[:, 5]
        particles_turns_SC[:, 6, i] .= pbeam.r[:, 6]

        delta_px .=  particles_turns[:,2,i] .- particles_turns_SC[:,2,i]


        NLOST[i] = sum(pbeam.lost_flag)
        
        get_emittance!(pbeam)
        records[i]=record(pbeam.centroid[1], pbeam.centroid[3], pbeam.centroid[5],
                        pbeam.beamsize[1], pbeam.beamsize[3], pbeam.beamsize[5],
                        pbeam.emittance[1], pbeam.emittance[2], pbeam.emittance[3], lumi)



        # save
        """if i % 100 == 0
            writedlm("emit1_SCBE_2.txt", new_emit1)
            writedlm("lost_SCBE_2.txt", NLOST)
        end
        """
        X1[i, :] = pbeam.r[1, 1:4] #track first particle
        next!(prog)


    end
end


########################################
plot(records.sx, label = L"$σ_x$", xlabel="# Turns")
plot(records.ex, label = L"$ε_x$", xlabel="# Turns")
plot(records.sz, label = L"$σ_z$", xlabel="# Turns")
plot(records.cx, label = L"$σ_x$ centroid", xlabel="# Turns")

plot(records.ez, label = L"$ε_z$", xlabel="# Turns")


##########################################
scatter(X1[:, 1], X1[:, 2])
plot(pbeam.r[:,1], pbeam.r[:,2], seriestype=:scatter, markersize=0.9, xlabel="x", ylabel= "px")
plot(pbeam.r[:,3], pbeam.r[:,4], seriestype=:scatter, markersize=0.9, xlabel="y", ylabel= "py")
plot(pbeam.r[:,5], pbeam.r[:,6], seriestype=:scatter, markersize=0.9, xlabel="z", ylabel= "δ")


delta_px
particles_turns_SC[:,2,:].- particles_turns[:,2,:]
scatter(pbeam.r[:,1], delta_px)


using PyCall
@pyimport NAFFlib

#half=turns ÷ 2
half = 9000

tunex_1=NAFFlib.multiparticle_tunes(particles_turns_SC[:, 1, 8000:half])
tunex_2=NAFFlib.multiparticle_tunes(particles_turns_SC[:, 1, half:end])
tuney_1=NAFFlib.multiparticle_tunes(particles_turns_SC[:, 3, 8000:half])
tuney_2=NAFFlib.multiparticle_tunes(particles_turns_SC[:, 3, half:end])

diff_tunex = sqrt.((tunex_1 .- tunex_2) .^2  .+ (tuney_1 .- tuney_2) .^2 )

maximum(diff_tunex)

scatter(tunex_1, tuney_1, marker_z = log10.(diff_tunex .+ 1e-15), markersize = .9,  color = :jet, clim=(-10,-2), 
    aspect_ratio=:equal, legend=:topleft, xlabel="Horizontal Tune", ylabel= "Vertical Tune", label="nSC = 10", dpi=300, size=(600,500),
    xlim = (0.225, 0.230), ylim= (0.207, 0.212))
    

maximum(tunex_1)-minimum(tunex_1)
maximum(tuney_1)-minimum(tuney_1)

# plot some particles last turns
pp = rand(1:1000, 5)

plot(particles_turns_SC[pp, 1, 8000:10000]', particles_turns_SC[pp, 2, 8000:10000]', seriestype=:scatter, markersize=2, xlabel="x", ylabel= "px", labels=["particle 1" "particle 2" "particle 3" "particle 4" "particle 5"])
plot(particles_turns_SC[pp, 3, 8000:10000]', particles_turns_SC[pp, 4, 8000:10000]', seriestype=:scatter, markersize=2, xlabel="y", ylabel= "py", labels=["particle 1" "particle 2" "particle 3" "particle 4" "particle 5"])

open("test_EIC_jutrack_tune_.jls", "w") do f
    serialize(f, particles_turns_SC)
end


"""
factorSC = 2*pbeam.classrad0/pbeam.beta^2/pbeam.gamma^3*sc_BE.ds*pbeam.np
pbeam.classrad0
#in EPIC = 1.535e-18
pbeam.gamma
# in EPIC 300, qui 540e3 ??? - beta molto simile
"""
