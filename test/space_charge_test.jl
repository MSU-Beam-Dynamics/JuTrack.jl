include("../src/JuTrack.jl")
using .JuTrack
using Distributions, Plots
include("../src/demo/ssrf_ring.jl")
RING = ssrf(-1.063770, 0)

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


new_RING = insert_space_charge(RING, pi/4, 10e-3, 10e-3, 15, 15)

D1 = DRIFT(len=0.2)
D2 = DRIFT(len=0.4)
D3 = DRIFT(len=0.2)
Q1 = KQUAD(len=0.1, k1=29.6)
Q2 = KQUAD(len=0.1, k1=-29.6)

a = 13e-3
b = 13e-3
nl = 12
nm = 12
SC_D1 = SPACECHARGE(effective_len=0.2, a=a, b=b, Nl=nl, Nm=nm)
SC_D2 = SPACECHARGE(effective_len=0.4, a=a, b=b, Nl=nl, Nm=nm)
SC_D3 = SPACECHARGE(effective_len=0.2, a=a, b=b, Nl=nl, Nm=nm)
SC_Q1 = SPACECHARGE(effective_len=0.1, a=a, b=b, Nl=nl, Nm=nm)
SC_Q2 = SPACECHARGE(effective_len=0.1, a=a, b=b, Nl=nl, Nm=nm)


line = [D1, Q1, D2, Q2, D3]
line_sc = [D1, SC_D1, Q1, SC_Q1, D2, SC_D2, Q2, SC_Q2, D3, SC_D3]

beam = Beam(zeros(5000, 6), energy=1.0e9, current=200.0, mass=m_p, charge=1.0, emittance=[1e-6, 1e-6, 0.0])

beta = beam.beta
gamma = beam.gamma
emit_norm = 1e-6
emit_phys = emit_norm / (beta * gamma)
beam.emittance = [emit_phys, emit_phys, 0.0]


vbase=3.42*8.5e6
ϕs=10.0
vact=vbase/cos(ϕs*π/180.0)
freq = 591e6
mainRFe=AccelCavity(freq, vact, 7560.0, π-ϕs*π/180.0)
tunex, tuney=50.08, 44.14
αc=3.42/tunex/tunex
lmap=LongitudinalRFMap(αc, mainRFe)
opt=optics4DUC(1.0, 0.0, 1.0, 0.0)
initilize_6DGaussiandist!(beam, opt, lmap)
beam1 = Beam(beam)


N = 1000
new_emit = zeros(N, 3)
new_emit1 = zeros(N, 3)
for i in 1:N
    linepass!(line, beam)
    linepass!(line_sc, beam1)
    get_emittance!(beam)
    get_emittance!(beam1)
    new_emit[i, :] = beam.emittance
    new_emit1[i, :] = beam1.emittance
end

plot(layout = (1, 2), legend = false, size = (1000, 400))
plot!(new_emit[:, 1].*1e6, title = "x emit", xlabel = "turns", ylabel = "emit (mm*mrad)", subplot = 1)
plot!(new_emit[:, 2].*1e6, title = "y emit", xlabel = "turns", ylabel = "emit (mm*mrad)", subplot = 2)
plot!(new_emit1[:, 1].*1e6, title = "x emit", xlabel = "turns", ylabel = "emit (mm*mrad)", subplot = 1)
plot!(new_emit1[:, 2].*1e6, title = "y emit", xlabel = "turns", ylabel = "emit (mm*mrad)", subplot = 2)

# # Extract x, px, y, py
# x = beam.r[:, 1]
# px = beam.r[:, 2]
# y = beam.r[:, 3]
# py = beam.r[:, 4]
# x1 = beam1.r[:, 1]
# px1 = beam1.r[:, 2]
# y1 = beam1.r[:, 3]
# py1 = beam1.r[:, 4]

# plot(layout = (2, 2), legend = false)

# plot!(x, px, seriestype = :scatter, subplot = 1, title = "x vs px", xlabel = "x", ylabel = "px", markersize = 1.0)
# plot!(y, py, seriestype = :scatter, subplot = 2, title = "y vs py", xlabel = "y", ylabel = "py", markersize = 1.0)
# plot!(x1, px1, seriestype = :scatter, subplot = 3, title = "x vs px", xlabel = "x", ylabel = "px", markersize = 1.0)
# plot!(y1, py1, seriestype = :scatter, subplot = 4, title = "y vs py", xlabel = "y", ylabel = "py", markersize = 1.0)