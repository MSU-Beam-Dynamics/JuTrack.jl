include("../src/JuTrack.jl")
using .JuTrack
include("ssrf_ring.jl")
using Plots  
using Enzyme
using ProgressMeter
using LinearAlgebra
using LaTeXStrings

RING = ssrf(-1.063770, 0)
# x: -0.04 to 0.04
xlist = [-0.04 + 0.001 * i for i in 0:80]
ylist = [0.0 + 0.001 * i for i in 0:40]
N = length(xlist) * length(ylist)

particles = zeros(length(xlist) * length(ylist), 6)
for i in 1:length(xlist)
    for j in 1:length(ylist)
        particles[(i-1)*length(ylist) + j, 1] = xlist[i]
        particles[(i-1)*length(ylist) + j, 3] = ylist[j]
    end
end

beam = Beam(copy(particles), energy=3.5e9)
pringpass!(RING, beam, 1)
survived = findall(x -> x == 0, beam.lost_flag)


function eigen_TPSA(xx)
    global RING
    x = CTPS(xx[1], 1, 6, 1)
    px = CTPS(xx[2], 2, 6, 1)
    y = CTPS(xx[3], 3, 6, 1)
    py = CTPS(xx[4], 4, 6, 1)
    z = CTPS(0.0, 5, 6, 1)
    delta = CTPS(0.0, 6, 6, 1)
    rin = [x, px, y, py, z, delta]
    ringpass_TPSA!(RING, rin, 1)
    map = [rin[1].map[2] rin[1].map[3] rin[1].map[4] rin[1].map[5];
            rin[2].map[2] rin[2].map[3] rin[2].map[4] rin[2].map[5];
            rin[3].map[2] rin[3].map[3] rin[3].map[4] rin[3].map[5];
            rin[4].map[2] rin[4].map[3] rin[4].map[4] rin[4].map[5]]
    e, v = qr_eigen(map)
    # println(map)
    return e[1], e[2]
end

E = zeros(ComplexF64, N, 2)
Vectors = zeros(ComplexF64, N, 2, 2)
grad = []
p = Progress(N, 1)
for i in 1:N
    next!(p)
    g = jacobian(Forward, final_x, [particles[i, 1], 0.0, particles[i, 3], 0.0], Val(4))
    # e = eigen([g[1][1] g[1][2] g[1][3] g[1][4]; g[2][1] g[2][2] g[2][3] g[2][4];
    #             g[3][1] g[3][2] g[3][3] g[3][4]; g[4][1] g[4][2] g[4][3] g[4][4]]).values
    # vectors = eigen([g[1][1] g[1][2] g[1][3] g[1][4]; g[2][1] g[2][2] g[2][3] g[2][4];
    #             g[3][1] g[3][2] g[3][3] g[3][4]; g[4][1] g[4][2] g[4][3] g[4][4]]).vectors
    e = eigen([g[1][1] g[1][2]; g[2][1] g[2][2]]).values
    vectors = eigen([g[1][1] g[1][2]; g[2][1] g[2][2]]).vectors
    Vectors[i, :, :] = vectors
    push!(grad, g)
    E[i, :] = e
    sleep(0.01)
end

