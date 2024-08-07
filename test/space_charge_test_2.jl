# include("../src/JuTrack.jl")
using JuTrack
using Distributions, Plots
using ProgressMeter
using DelimitedFiles
using LinearAlgebra

# include("../src/demo/ssrf_ring.jl")
# ring = ssrf(-1.063770, 0)

Dr = DRIFT(name="Dr", len=0.5)
HalfDr = DRIFT(name="HalfDr", len=0.25)
p2Dr = DRIFT(name="p2Dr", len=0.2)
SF = KSEXT(name="SF", len=0.1, k2=1.0)
SD = KSEXT(name="SD", len=0.1, k2=-1.0)
B1 = SBEND(name="B", len= 1.0, angle=2*pi/40.0)
Q1 = QUAD(name="Q1", len=0.5, k1=1.2) # optimized k1 starting from -1.0
Q2 = QUAD(name="Q2", len=0.5, k1=-1.2)

cell = [HalfDr, B1, p2Dr, SF, p2Dr, Q1, Dr, B1, p2Dr, SD, p2Dr, Q2, HalfDr]
ring = [cell..., cell..., cell..., cell..., cell..., cell..., cell..., cell..., cell..., cell...,
        cell..., cell..., cell..., cell..., cell..., cell..., cell..., cell..., cell..., cell...]

function insert_space_charge(lattice, dphi, a, b, Nl, Nm)
    twi = twissring(lattice, 0.0, 1)
    phi0 = 0.0
    s0 = 0.0
    s = spos(lattice)
    insert_idx = [] 
    insert_len = []

    if twi[end].dmux < dphi
        println("The phase advance of the whole ring is less than ", dphi, ", only one space charge element is inserted.")
        sc = SPACECHARGE(effective_len=s[end], a=a, b=b, Nl=Nl, Nm=Nm)
        new_lattice = [lattice..., sc]
        return new_lattice
    end
    for i in 1:length(lattice)
        if twi[i].dmux - phi0 > dphi
            len = s[i] - s0
            push!(insert_idx, i)
            push!(insert_len, len)
            s0 = s[i]
            phi0 = twi[i].dmux
        end
    end
    # new list of AbstractElement
    new_lattice = Vector{AbstractElement}()
    for i in 1:length(lattice)
        push!(new_lattice, lattice[i])
        if i in insert_idx
            num = findfirst(x -> x == i, insert_idx)
            sc = SPACECHARGE(effective_len=insert_len[num], a=a, b=b, Nl=Nl, Nm=Nm)
            push!(new_lattice, sc)
        end
    end
    println("Space charge elements are inserted at: ", insert_idx)
    return new_lattice
end


new_RING = insert_space_charge(ring, pi/4, 10e-3, 10e-3, 15, 15)



beam = Beam(zeros(10001, 6), energy=1.0e9, current=50.0, mass=m_p, charge=1.0, emittance=[1e-6, 1e-6, 0.0])

beta = beam.beta
gamma = beam.gamma
# emit_norm = 1e-8
# emit_phys = emit_norm  / (beta * gamma)
# beam.emittance = [emit_phys, emit_phys, 0.0]


vbase=3.42*8.5e6
ϕs=10.0
vact=vbase/cos(ϕs*π/180.0)
freq = 591e6
mainRFe=AccelCavity(freq, vact, 7560.0, π-ϕs*π/180.0)
tunex, tuney=50.08, 44.14
αc=3.42/tunex/tunex
lmap=LongitudinalRFMap(αc, mainRFe)
opt=optics4DUC(10.0, 0.0, 10.0, 0.0)
initilize_6DGaussiandist!(beam, opt, lmap)
beam.r[:, 5] .= 0.0
beam.r[:, 6] .= 0.0
beam.r[1, :] .= 0.0 # reference particle
beam1 = Beam(beam)

init_r = copy(beam.r)
# tracking
N = 2000

rout0 = pringpass!(ring, beam, N, true)
rout1 = pringpass!(new_RING, beam1, N, true)


function get_tune(RING, rout, sur_idx, nturns)
    twi = twissring(RING, 0.0, 1)
    np = length(sur_idx)
    betax = twi[1].betax
    betay = twi[1].betay
    alphax = twi[1].alphax
    alphay = twi[1].alphay

    x = zeros(nturns, np)
    px = zeros(nturns, np)
    y = zeros(nturns, np)
    py = zeros(nturns, np)
    for i in 1:nturns
        x[i,:] = rout[i][sur_idx,1] ./ sqrt(betax)
        px[i,:] = -rout[i][sur_idx,2] .* sqrt(betax) .- alphax .* x[i]
        y[i,:] = rout[i][sur_idx,3] ./ sqrt(betay)
        py[i,:] = -rout[i][sur_idx,4] .* sqrt(betay) .- alphay .* y[i]
    end

    nux1 = zeros(np)
    nux2 = zeros(np)
    nuy1 = zeros(np)
    nuy2 = zeros(np)
    for i in 1:np
        nux1[i], nux2[i] = naff(nturns, x[:,i], px[:,i])
        nuy1[i], nuy2[i] = naff(nturns, y[:,i], py[:,i])
    end

    diff_nux = log10.((nux2 .- nux1).^2 .+ (nuy2 .- nuy1).^2)
    return diff_nux, nux1, nux2, nuy1, nuy2
end

survived = findall(x -> x == 0, beam1.lost_flag)

diff_nux, nux1, nux2, nuy1, nuy2 = get_tune(new_RING, rout1, survived, N)
diff_nux0, nux10, nux20, nuy10, nuy20 = get_tune(ring, rout0, survived, N)
scatter(init_r[survived, 1], init_r[survived, 3], zcolor=diff_nux, markersize=2)
scatter(nux1, nuy1, label="nux1-nuy1")

nux1_copy = copy(nux1)
nuy1_copy = copy(nuy1)
nux2_copy = copy(nux2)
nuy2_copy = copy(nuy2)

for i in 1:length(nux1_copy)
    if nux1_copy[i] > 0.5
        nux1_copy[i] = 1.0 - nux1_copy[i]
    end
    if nuy1_copy[i] < 0.5
        nuy1_copy[i] = 1.0 - nuy1_copy[i]
    end
    if nux2_copy[i] > 0.5
        nux2_copy[i] = 1.0 - nux2_copy[i]
    end
    if nuy2_copy[i] < 0.5
        nuy2_copy[i] = 1.0 - nuy2_copy[i]
    end
end

# scatter(nux1_copy, nuy1_copy, label="W/ SC", xlabel="nux", ylabel="nuy", 
#      xlim=(0.19, 0.25), ylim=(0.86, 0.95), markersize=2, figsize=(600, 400))
# scatter!([0.2199], [0.9178], marker=:p, label="W/O SC")

# scatter(init_r[survived, 1]*1e3, init_r[survived, 3]*1e3, markersize=2, legend=false, xlabel="x (m)", ylabel="y (m)")

using PyCall
np = pyimport("numpy")
plt = pyimport("matplotlib.pyplot")

plt.figure(figsize=(5, 3.5))
plt.scatter(nux1_copy, nuy1_copy, s=5, label="W/ SC")
plt.scatter([0.2199], [0.9178],marker="*", s=15, label="W/O SC")
plt.xlabel("nux", fontsize=16, fontname="Times New Roman")
plt.ylabel("nuy", fontsize=16, fontname="Times New Roman")
plt.xticks(fontsize=14, fontname="Times New Roman")
plt.yticks(fontsize=14, fontname="Times New Roman")
plt.xlim(0.19, 0.25)
plt.ylim(0.86, 0.95)
plt.legend(prop=Dict("family"=>"Times New Roman"))
plt.tight_layout()
plt.show()
# plt.figure(figsize=(6, 4))
# plt.scatter(init_r[survived, 1]*1e3, init_r[survived, 3]*1e3, c=diff_nux, s=2)
# plt.xlabel("x (mm)", fontsize=16, fontname="Times New Roman")
# plt.ylabel("y (mm)", fontsize=16, fontname="Times New Roman")
# plt.xticks(fontsize=14, fontname="Times New Roman")
# plt.yticks(fontsize=14, fontname="Times New Roman")
# plt.colorbar()
# plt.show()