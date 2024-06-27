include("../src/JuTrack.jl")
using .JuTrack
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

    new_lattice = []
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

xlist = [i*0.001 - 0.02 for i in 1:40]
ylist = [i*0.001 for i in 1:20]
particles = zeros(800, 6)
for i in 1:40
    for j in 1:20
        particles[(i-1)*20+j, 1] = xlist[i]
        particles[(i-1)*20+j, 3] = ylist[j]
    end
end
beam = Beam(copy(particles), 3.5e9, current=3.0)
beam1 = Beam(copy(particles), 3.5e9, current=3.0)

diff_nux, nux1, nux2, nuy1, nuy2 = FMA(RING, beam, 1000)
diff_nux1, nux11, nux21, nuy11, nuy21 = FMA(new_RING, beam1, 1000)