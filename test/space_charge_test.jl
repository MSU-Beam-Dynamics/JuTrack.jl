# include("../src/JuTrack.jl")
using JuTrack
using Distributions, Plots
using ProgressMeter
using DelimitedFiles

include("../src/demo/ssrf_ring.jl")
# RING = ssrf(-1.063770, 0)

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


# new_RING = insert_space_charge(RING, pi/4, 10e-3, 10e-3, 15, 15)

D1 = DRIFT(len=0.2)
D2 = DRIFT(len=0.4)
D3 = DRIFT(len=0.2)
Q1 = QUAD(len=0.1, k1=29.037)
Q2 = QUAD(len=0.1, k1=-29.037)
line = [D1, Q1, D2, Q2, D3]
M66 = findm66(line, 0.0, 1)
phase = acos((M66[1, 1] + M66[2, 2]) / 2)
println("The phase advance of the whole ring is ", phase * 180 / pi, " degree.")

a = 10e-3
b = 10e-3
nl = 15
nm = 15
SC_D1 = SPACECHARGE(effective_len=0.2, a=a, b=b, Nl=nl, Nm=nm)
SC_D2 = SPACECHARGE(effective_len=0.4, a=a, b=b, Nl=nl, Nm=nm)
SC_D3 = SPACECHARGE(effective_len=0.2, a=a, b=b, Nl=nl, Nm=nm)
SC_Q1 = SPACECHARGE(effective_len=0.1, a=a, b=b, Nl=nl, Nm=nm)
SC_Q2 = SPACECHARGE(effective_len=0.1, a=a, b=b, Nl=nl, Nm=nm)





line_sc = [D1, SC_D1, Q1, SC_Q1, D2, SC_D2, Q2, SC_Q2, D3, SC_D3]

beam = Beam(zeros(50000, 6), energy=1.0e9, current=450.0, mass=m_p, charge=1.0, emittance=[1e-6, 1e-6, 0.0])

beta = beam.beta
gamma = beam.gamma
emit_norm = 1e-6
emit_phys = emit_norm  / (beta * gamma)
beam.emittance = [emit_phys, emit_phys, 0.0]


vbase=3.42*8.5e6
ϕs=10.0
vact=vbase/cos(ϕs*π/180.0)
freq = 591e6
mainRFe=AccelCavity(freq, vact, 7560.0, π-ϕs*π/180.0)
tunex, tuney=50.08, 44.14
αc=3.42/tunex/tunex
lmap=LongitudinalRFMap(αc, mainRFe)
opt=optics4DUC(5.0, 0.0, 5.0, 0.0)
initilize_6DGaussiandist!(beam, opt, lmap)
beam.r[:, 5] .= 0.0
beam.r[:, 6] .= 0.0
beam1 = Beam(beam)


N = 50000
new_emit = zeros(N+1, 3)
new_emit1 = zeros(N+1, 3)
NLOST = zeros(N)

println("Start tracking")
prog = Progress(N)
for i in 1:N
    # println("Turn: ", i)
    if i == 1
        get_emittance!(beam)
        get_emittance!(beam1)
        new_emit[i, :] = beam.emittance
        new_emit1[i, :] = beam1.emittance
    end
    linepass!(line, beam)
    linepass!(line_sc, beam1)
    NLOST[i] = sum(beam1.lost_flag)
    get_emittance!(beam)
    get_emittance!(beam1)
    new_emit[i+1, :] = beam.emittance
    new_emit1[i+1, :] = beam1.emittance
    next!(prog)
    if mod(i, 1000) == 0
        all_emit = [new_emit new_emit1]
        writedlm("emit_withoutSC_withSC_1GeV_50000n_450A.txt", all_emit)
    end
end

# using DelimitedFiles
# # writedlm("emit.txt", new_emit)
# new_emit = readdlm("emit.txt")

# p1=plot(0:200, new_emit[1:100:end, 1].*1e6, title = "x emit", xlabel = "turns*100", ylabel = "emit (mm*mrad)", label = "without SC")
# plot!(0:200, new_emit1[1:100:end, 1].*1e6, title = "x emit", xlabel = "turns*100", ylabel = "emit (mm*mrad)",  label = "with SC")
# p2=plot(0:200, new_emit1[1:100:end, 2].*1e6, title = "y emit", xlabel = "turns*100", ylabel = "emit (mm*mrad)", label = "with SC")
# plot!(0:200, new_emit[1:100:end, 2].*1e6, title = "y emit", xlabel = "turns*100", ylabel = "emit (mm*mrad)", label = "without SC")
# plot(p1, p2, layout = (1, 2))

using PyCall
np = pyimport("numpy")
plt = pyimport("matplotlib.pyplot")
N=25000
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
plt.legend()
plt.subplot(1, 2, 2)
plt.plot(np.arange(N+1), new_emit[1:N+1, 2].*1e6, label="without SC")
plt.plot(np.arange(N+1), new_emit1[1:N+1, 2].*1e6, label="with SC")
plt.xlabel("periods", fontsize=16, fontname="Times New Roman")
plt.ylabel("emit (mm*mrad)", fontsize=16, fontname="Times New Roman")
plt.title("y emittance", fontsize=16, fontname="Times New Roman")
plt.legend(prop=Dict("family"=>"Times New Roman"))
plt.xticks(fontsize=14, fontname="Times New Roman")
plt.yticks(fontsize=14, fontname="Times New Roman")
# plt.yscale("log")
plt.tight_layout()
plt.show()

emit_growth = zeros(N)
for i in 1:N
    emit_growth[i] = (new_emit1[i+1, 1] / new_emit[1, 1] * new_emit1[i+1, 2] / new_emit[1, 2] - 1) * 100
end
plt.plot(np.arange(N), emit_growth, label="emit growth")
plt.xlabel("periods", fontsize=16, fontname="Times New Roman")
plt.ylabel("emit growth (%)", fontsize=16, fontname="Times New Roman")
plt.xticks(fontsize=14, fontname="Times New Roman")
plt.yticks(fontsize=14, fontname="Times New Roman")
plt.show()

idx_sur = findall(x -> x == 0, beam1.lost_flag)
plt.figure(figsize=(11, 8))
plt.subplot(2, 2, 1)
plt.scatter(beam.r[:, 1], beam.r[:, 2], s=0.1)
plt.xlabel("x (m)", fontsize=16, fontname="Times New Roman")
plt.ylabel("px", fontsize=16, fontname="Times New Roman")
plt.title("Without space charge", fontsize=16, fontname="Times New Roman")
plt.xlim(-0.008, 0.008)
plt.ylim(-0.017, 0.017)
plt.subplot(2, 2, 2)
plt.scatter(beam1.r[idx_sur, 1], beam1.r[idx_sur, 2], s=0.1)
plt.xlabel("x (m)", fontsize=16, fontname="Times New Roman")
plt.ylabel("px", fontsize=16, fontname="Times New Roman")
plt.title("With space charge", fontsize=16, fontname="Times New Roman")
plt.xlim(-0.008, 0.008)
plt.ylim(-0.017, 0.017)
plt.subplot(2, 2, 3)
plt.scatter(beam.r[:, 3], beam.r[:, 4], s=0.1)
plt.xlabel("y (m)", fontsize=16, fontname="Times New Roman")
plt.ylabel("py", fontsize=16, fontname="Times New Roman")
plt.title("Without space charge", fontsize=16, fontname="Times New Roman")
plt.xlim(-0.008, 0.008)
plt.ylim(-0.017, 0.017)
plt.subplot(2, 2, 4)
plt.scatter(beam1.r[idx_sur, 3], beam1.r[idx_sur, 4], s=0.1)
plt.xlabel("y (m)", fontsize=16, fontname="Times New Roman")
plt.ylabel("py", fontsize=16, fontname="Times New Roman")
plt.title("With space charge", fontsize=16, fontname="Times New Roman")
plt.xlim(-0.008, 0.008)
plt.ylim(-0.017, 0.017)
plt.tight_layout()
plt.show()

