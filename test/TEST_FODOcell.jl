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

    data_fodo= readdlm("test/fodocell_parameter.txt")
    betax = data_fodo[1, :]
    betay = data_fodo[2,:]
    alphax = data_fodo[3,:]
    alphay = data_fodo[4,:]
    dmux = data_fodo[5,:]
    dmuy = data_fodo[6,:]
end

begin
    num_particles = 1000
    turns = 1000

    pbeam = Beam(zeros(1000, 6), np=Int(0.25e11), energy=1e9, emittance = [5.75424e-7, 4.82822e-7, 2.3e-5])
    opbeam = optics4DUC(0.746, -1.396, 0.893, 1.473)
    mainRF = AccelCavity(freq=10e6, volt=85e3, h=1.0, phis=Float64(0.0))
    αc=0.0
    lmap = LongitudinalRFMap(αc, mainRF)
    initilize_6DGaussiandist!(pbeam, opbeam, lmap)

    # define oneturnmap
    tunex=0.2495
    tuney=0.2030
    oneturn = TransferMap4DChrom(opbeam, tunex, tuney, 0.0, 0.0)

    # single element optics
    opD1 = optics4DUC(betax[1], alphax[1], betay[1], alphay[1])
    opQ1 = optics4DUC(betax[2], alphax[2], betay[2], alphay[2])
    opD2 = optics4DUC(betax[3], alphax[3], betay[3], alphay[3])
    opQ2 = optics4DUC(betax[4], alphax[4], betay[4], alphay[4])
    opD3 = optics4DUC(betax[5], alphax[5], betay[5], alphay[5]) #ritorna uguale a opbeam
    
    # single element phase advance
    Qx_D1 = dmux[1]
    Qx_Q1 = dmux[2]
    Qx_D2 = dmux[3]
    Qx_Q2 = dmux[4]
    Qx_D3 = dmux[5]

    Qy_D1 = dmuy[1]
    Qy_Q1 = dmuy[2]
    Qy_D2 = dmuy[3]
    Qy_Q2 = dmuy[4]
    Qy_D3 = dmuy[5]

    TM_D1= TransferMap4D(opbeam, opD1, Qx_D1, Qy_D1)
    TM_Q1= TransferMap4D(opbeam, opQ1, Qx_Q1, Qy_Q1)
    TM_D2= TransferMap4D(opbeam, opD2, Qx_D2, Qy_D2)
    TM_Q2= TransferMap4D(opbeam, opQ2, Qx_Q2, Qy_Q2)
    TM_D3= TransferMap4D(opbeam, opD3, Qx_D3, Qy_D3)
    
    TM_D1_inv= Inverse_TransferMap4D(opbeam, opD1, Qx_D1, Qy_D1)
    TM_Q1_inv= Inverse_TransferMap4D(opbeam, opQ1, Qx_Q1, Qy_Q1)
    TM_D2_inv= Inverse_TransferMap4D(opbeam, opD2, Qx_D2, Qy_D2)
    TM_Q2_inv= Inverse_TransferMap4D(opbeam, opQ2, Qx_Q2, Qy_Q2)
    TM_D3_inv= Inverse_TransferMap4D(opbeam, opD3, Qx_D3, Qy_D3)
    
    l_D1 = 0.2
    l_Q1 = 0.1
    l_D2 = 0.4
    l_Q2 = 0.1
    l_D3 = 0.2

    nsc_D1 = 1
    nsc_Q1 = 1
    nsc_D2 = 1
    nsc_Q2 = 1
    nsc_D3 = 1

    opSC = optics4DUC(1.0, 0.0, 1.0, 0.0) # non mi serve in single element

    sc_D1 = SC_lens(opSC, l_D1, nsc_D1)
    sc_Q1 = SC_lens(opSC, l_Q1, nsc_Q1)
    sc_D2 = SC_lens(opSC, l_D2, nsc_D2)
    sc_Q2 = SC_lens(opSC, l_Q2, nsc_Q2)
    sc_D3 = SC_lens(opSC, l_D3, nsc_D3)

    # records
    particles_turns = zeros(Float64, num_particles, 6, turns)
    particles_turns_SC = zeros(Float64, num_particles, 6, turns)
    delta_px = zeros(Float64, num_particles)
    records = StructArray{record}(undef, turns)

    ring= [TM_D1, sc_D1, TM_D1_inv, TM_Q1, sc_Q1, TM_Q1_inv, TM_D2, sc_D2, TM_D2_inv,
            TM_Q2, sc_Q2, TM_Q2_inv, TM_D3, sc_D3]

end

get_emittance!(pbeam)
pbeam.emittance
pbeam.beamsize

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
pbeam.beamsize
pbeam.emittance

plot(new_emit1[:,1])
plot(X1[:, 1], X1[:, 2])
plot(pbeam.r[:,1], pbeam.r[:,2], seriestype=:scatter, markersize=0.9, xlabel="x", ylabel= "px")
plot(pbeam.r[:,3], pbeam.r[:,4], seriestype=:scatter, markersize=0.9, xlabel="x", ylabel= "px")
plot(pbeam.r[:,5], pbeam.r[:,6], seriestype=:scatter, markersize=0.9, xlabel="x", ylabel= "px")

plot(particles_turns_SC[:, 1, :], delta_px, seriestype=:scatter, markersize=0.9, xlabel="x", ylabel= "Δpx")



using JuTrack
using Distributions
using ProgressMeter
using DelimitedFiles
using PyCall
np = pyimport("numpy")
plt = pyimport("matplotlib.pyplot")

begin

    D1L = 0.2
    D2L = 0.4
    D3L = 0.2
    Q1L = 0.1
    Q2L = 0.1
    Q1k = 29.6
    Q2k = -29.6

    a = 13e-3
    b = 13e-3
    nl = 15
    nm = 15

    D1 = DRIFT(len=D1L)
    D2 = DRIFT(len=D2L)
    D3 = DRIFT(len=D3L)
    Q1 = KQUAD(len=Q1L, k1=Q1k)
    Q2 = KQUAD(len=Q2L, k1=Q2k)

    opSC = optics4DUC(1.0, 0.0, 1.0, 0.0) # non mi serve in single element

    sc_D1 = SC_lens(opSC, D1L, 1)
    sc_Q1 = SC_lens(opSC, Q1L, 1)
    sc_D2 = SC_lens(opSC, D2L, 1)
    sc_Q2 = SC_lens(opSC, Q2L, 1)
    sc_D3 = SC_lens(opSC, D3L, 1)


    line_SC = [D1, sc_D1, Q1, sc_Q1, D2, sc_D2, Q2, sc_Q2, D3, sc_D3]

    distparam = [
        3.677529920673089E-004  ,   # sigx
        8.428925532276500E-004 ,   # sigpx
        -0.828277121044551 ,   # muxpx
        1.0,   # xscale
        1.0,   # pxscale
        0.0,   # xmu1 (mean x)
        0.0,   # xmu2 (mean px)
        3.677529304933903E-004  ,   # sigy
        8.428931246578997E-004  ,   # sigpy
        0.828276927537804 ,   # muypy
        1.0,   # yscale
        1.0,   # pyscale
        0.0,   # xmu3 (mean y)
        0.0,   # xmu4 (mean py)
        1.0,   # sigz
        0.1,   # sigpz
        0.5,   # muzpz
        0.0,   # zscale
        0.0,   # pzscale
        0.0,   # xmu5 (mean z)
        0.0    # xmu6 (mean pz)
    ]
    Npt = 1000
    Pts1 = Gauss3_Dist(distparam, Npt, seed=1234)

    beam = Beam(Pts1, energy=1.0e9, current=20.0, mass=m_p, charge=1.0)
    beam1 = Beam(beam)

    N = 1000
    new_emit1 = zeros(N+1, 3)
    NLOST = zeros(N)

end


println("Start tracking")
prog = Progress(N)
X1 = zeros(N, 4)

for i in 1:N
    
    # println("Turn: ", i)
    if i == 1
        get_emittance!(beam1)
        new_emit1[i, :] = beam1.emittance
    end

    linepass!(line_SC, beam1)

    NLOST[i] = sum(beam1.lost_flag)
    get_emittance!(beam1)
    new_emit1[i+1, :] = beam1.emittance
    # save
    if i % 100 == 0
        writedlm("emit1.txt", new_emit1)
        writedlm("lost.txt", NLOST)
    end
    X1[i, :] = beam1.r[1, 1:4]
    next!(prog)
end

using Plots
plot(new_emit1)


plt.figure(figsize=(9, 4))
plt.subplot(1, 2, 1)
plt.plot(np.arange(N+1), new_emit[1:N+1, 1].*1e6, label="without SC")
plt.plot(np.arange(N+1), new_emit1[1:N+1, 1].*1e6, label="with SC")
plt.xlabel("periods", fontsize=16, fontname="Times New Roman")
plt.ylabel("emit (mm*mrad)", fontsize=16, fontname="Times New Roman")
plt.title("x emittance", fontsize=16, fontname="Times New Roman")
plt.xticks(fontsize=14, fontname="Times New Roman")
plt.yticks(fontsize=14, fontname="Times New Roman")
# plt.yscale("log")
plt.legend(prop=Dict("family"=>"Times New Roman"), frameon=false)
# plt.ylim(0.172, 0.175)
plt.subplot(1, 2, 2)
plt.plot(np.arange(N+1), new_emit[1:N+1, 2].*1e6, label="without SC")
plt.plot(np.arange(N+1), new_emit1[1:N+1, 2].*1e6, label="with SC")
plt.xlabel("periods", fontsize=16, fontname="Times New Roman")
plt.ylabel("emit (mm*mrad)", fontsize=16, fontname="Times New Roman")
plt.title("y emittance", fontsize=16, fontname="Times New Roman")
plt.legend(prop=Dict("family"=>"Times New Roman"), frameon=false)
plt.xticks(fontsize=14, fontname="Times New Roman")
plt.yticks(fontsize=14, fontname="Times New Roman")
# plt.ylim(0.172, 0.175)
# plt.yscale("log")
plt.tight_layout()
plt.show()

emit_growth = zeros(N)
for i in 1:N
    emit_growth[i] = (new_emit1[i+1, 1] / new_emit1[1, 1] * new_emit1[i+1, 2] / new_emit1[1, 2] - 1) * 100
end
plt.plot(np.arange(10000), emit_growth[1:10000], label="emit growth 5000")
plt.xlabel("periods", fontsize=16, fontname="Times New Roman")
plt.ylabel("emit growth (%)", fontsize=16, fontname="Times New Roman")
plt.xticks(fontsize=14, fontname="Times New Roman")
plt.yticks(fontsize=14, fontname="Times New Roman")
plt.tight_layout()
plt.show()