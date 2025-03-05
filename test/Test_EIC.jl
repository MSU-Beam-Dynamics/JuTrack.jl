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

############## start tracking #################
for i in 1:turns
    #if i < turns_rec
    track!(pbeam, oneturn)

    if bb == true
        pass!(lmap, rin, num_particles, pbeam)
        pass!(pcrabu, rin, num_particles, pbeam) #questi cambiano mi sa
        pass!(pcrabu2nd, rin, num_particles, pbeam)
        pass!(lb, rin, num_particles, pbeam)
        
        lumi=linepass!([estrong], pbeam)
        
        pass!(invlb, rin, num_particles, pbeam)
        pass!(pcrabd, rin, num_particles, pbeam)
        pass!(pcrabd2nd, rin, num_particles, pbeam)
    end

    #poi togliere!!!
    particles_turns[:, 1, i] .= pbeam.r[:, 1]
    particles_turns[:, 2, i] .= pbeam.r[:, 2]
    particles_turns[:, 3, i] .= pbeam.r[:, 3]
    particles_turns[:, 4, i] .= pbeam.r[:, 4]

    for j in 1:nSC
    
        TM1[j] = TransferMap4D(opIPp, opSC, phi_advx[j], phi_advy[j])
        track!(pbeam, TM1[j])
        
        #sc_BE.sigma_xi_SC[i,j], sc_BE.sigma_yi_SC[i,j] = pass!(sc_BE, rin, num_particles, pbeam)
        pass!(sc_BE, rin, num_particles, pbeam)
        
        TM1_inv[j] = Inverse_TransferMap4D(opIPp, opSC, phi_advx[j], phi_advy[j])
        track!(pbeam, TM1_inv[j])

    end

        # record SC beam size in file
        #open("pbeam_records_sx_SC_Sympl_TEST.jls", "w") do f
        #    serialize(f, sc_BE.sigma_xi_SC)
        #end

        #open("pbeam_records_sy_SC_Sympl_TEST.jls", "w") do f
        #    serialize(f, sc_BE.sigma_yi_SC)
        #end
        
    """elseif i == turns_rec
        for j in 1:nSC
            sc_BE.sigma_x_SC[j] = mean(sc_BE.sigma_xi_SC[200:turns_rec,j])
            sc_BE.sigma_y_SC[j] = mean(sc_BE.sigma_yi_SC[200:turns_rec,j])
        end


    elseif i > turns_rec

        track!(pbeam, oneturn)
        
        if bb == true
            
            pass!(pbeam, lmap)
            pass!(pbeam, pcrabu)
            pass!(pbeam, pcrabu2nd)
            pass!(pbeam, lb)
            lumi=pass!(pbeam, estrong)
            pass!(pbeam, invlb)
            pass!(pbeam, pcrabd)
            pass!(pbeam, pcrabd2nd)
        end
    
        particles_turns_SC[:, 1, i] .= pbeam.r[:, 1]
        particles_turns_SC[:, 2, i] .= pbeam.r[:, 2]
        particles_turns_SC[:, 3, i] .= pbeam.r[:, 3]
        particles_turns_SC[:, 4, i] .= pbeam.r[:, 4]
        

        #SC BE model 
        for j in 1:(nSC)
            
            TM1[j] = TransferMap4D(opIPp, opSC, phi_advx[j], phi_advy[j])
            track!(pbeam, TM1[j])
            
            sc_BE.sigma_x_SC[j], sc_BE.sigma_y_SC[j] = SC_kick!(sc_BE, pbeam, w, sc_BE.sigma_x_SC[j], sc_BE.sigma_y_SC[j])
            

            TM1_inv[j] = Inverse_TransferMap4D(opIPp, opSC, phi_advx[j], phi_advy[j])
            track!(pbeam, TM1_inv[j])
        end
        """
    particles_turns_SC[:, 1, i] .= pbeam.r[:, 1]
    particles_turns_SC[:, 2, i] .= pbeam.r[:, 2]
    particles_turns_SC[:, 3, i] .= pbeam.r[:, 3]
    particles_turns_SC[:, 4, i] .= pbeam.r[:, 4]
    
    #if i > turns_rec

    """get_emittance!(pbeam)
    records[i]=record(pbeam.centroid[1], pbeam.centroid[3], pbeam.centroid[5],
                        pbeam.beamsize[1], pbeam.beamsize[3], pbeam.beamsize[5],
                        pbeam.emittance[1], pbeam.emittance[2], pbeam.emittance[3], lumi)

    delta_px .=  particles_turns[:,2,i] .- particles_turns_SC[:,2,i] """
    #end
        
    #end
 
end

get_emittance!(pbeam)
println(pbeam.emittance)


# plot distribution

plot(beam1.r[:,1], beam1.r[:,2], seriestype=:scatter, markersize=0.9, xlabel="x", ylabel= "px")
plot(beam1.r[:,3], beam1.r[:,4], seriestype=:scatter, markersize=0.9, xlabel="y", ylabel= "py")
plot(beam1.r[:,1], beam1.r[:,3], seriestype=:scatter, markersize=0.9, xlabel="y", ylabel= "py")
plot(beam1.r[:,5], beam1.r[:,6], seriestype=:scatter, markersize=0.9, xlabel="z", ylabel= "δ")


plot(beam1.r[:,1], delta_px, seriestype=:scatter, markersize=0.9, xlabel="x", ylabel= "Δpx")



#########################################

plot(records.sx, label = L"$σ_x$", xlabel="# Turns")
plot(records.ex, label = L"$ε_x$", xlabel="# Turns")
plot(records.sz, label = L"$σ_z$", xlabel="# Turns")
plot(records.cx, label = L"$σ_x$ centroid", xlabel="# Turns")

##########################################
using PyCall
@pyimport NAFFlib

#half = turns ÷ 2
half = (turns-turns_rec) ÷ 2

tunex_1=NAFFlib.multiparticle_tunes(particles_turns_SC[:, 1, turns_rec+1:half ])
tunex_2=NAFFlib.multiparticle_tunes(particles_turns_SC[:, 1, half:end])
tuney_1=NAFFlib.multiparticle_tunes(particles_turns_SC[:, 3, turns_rec+1:half])
tuney_2=NAFFlib.multiparticle_tunes(particles_turns_SC[:, 3, half:end])

diff_tunex = sqrt.((tunex_1 .- tunex_2) .^2  .+ (tuney_1 .- tuney_2) .^2 )

maximum(diff_tunex)

scatter(tunex_1, tuney_1, marker_z = log10.(diff_tunex .+ 1e-15), markersize = .7,  color = :jet, clim=(-10,-2), 
aspect_ratio=:equal, legend=:topleft, xlabel="Horizontal Tune", ylabel= "Vertical Tune", label="nSC = 10", dpi=300)#,  xlim=(0.22, 0.23), ylim=(0.185, 0.215))

maximum(tunex_1)-minimum(tunex_1)
maximum(tuney_1)-minimum(tuney_1)

####################################################################

fftx=fft(@view records_after.cx[5001:end]) # N/2 + 1 to N/2
ffty=fft(@view records_after.cy[5001:end])
fftz=fft(@view records_after.cz[5001:end])
fftf=fftfreq(5000, 1.0)

turns = 10000
plot(@view(fftf[1:turns÷2]), abs.(@view fftx[1:turns÷2]), xlim=(0.2,0.25), legendfontsize = 10, label = L"$x$ $centroid$", xlabel =L"$tune$")
plot(@view(fftf[1:turns÷2]), abs.(@view ffty[1:turns÷2]), xlim=(0.2,0.25), legendfontsize = 10, label = L"$y$ $centroid$", xlabel =L"$tune$")


records_after = open(deserialize,"pbeam_records_BBSC_10kturns_10kpart_w05.jls")
plot(records_after.sx, legendfontsize = 10, label = L"$σ_x$", xlabel="# Turns")
plot(records_after.ex, legendfontsize = 10, label = L"$ε_x$", xlabel="# Turns")