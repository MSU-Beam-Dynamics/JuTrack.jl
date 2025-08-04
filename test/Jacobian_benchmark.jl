# Benchmarking the Jacobian calculation in JuTrack.jl
# Test case: single pass through the SPEAR3 lattice with one particle.
# Number of variables: 6 (x0, xp0, y0, yp0, z0, delta0)
# Number of output variables: 6 (xf, xpf, yf, ypf, zf, deltaf)
# Three automatic differentiation methods are compared:
# 1. JuTrack's first-order TPSA
# 2. Enzyme's Forward mode AD
# 3. JuTrack's high-order TPSA
# This benchmark test is done on 08/01/2025 with the latest version of JuTrack.jl
# with Julia version 1.10.4, Enzyme version v0.13.3
include("../src/demo/SPEAR3/spear3.jl")
using JuTrack
using BenchmarkTools

set_tps_dim(6)
RING = spear3()
RING_TPSAD = Number2TPSAD(RING)

###############################################################
# Calculate the Jacobian using JuTrack's first-order TPSA
###############################################################
function fast_TPSA_jacobian(x, xp, y, yp, z, delta)
    beam = TBeam([x xp y yp z delta])
    ringpass!(RING_TPSAD, beam, 1)
    return beam.r[1, :]
end
@btime Jacobian(fast_TPSA_jacobian, [0.001, 0.0, 0.001, 0.0, 0.0, 0.0])
# 2.568 ms (33856 allocations: 8.39 MiB)

###############################################################
# Calculate the Jacobian using Enzyme AD
###############################################################
function AD_jacobian(x)
    beam = Beam([x[1] x[2] x[3] x[4] x[5] x[6]])
    ringpass!(RING, beam, 1)
    return beam.r
end
@btime jacobian(Forward, AD_jacobian, [0.001, 0.0, 0.001, 0.0, 0.0, 0.0])
# 3.091 ms (21211 allocations: 1.65 MiB)

##############################################################
# Calculate the Jacobian using high-order TPSA
##############################################################
function TPSA_track_jacobian(x, y)
    X = CTPS(x, 1, 6, 1)
    PX = CTPS(0.0, 2, 6, 1)
    Y = CTPS(y, 3, 6, 1)
    PY = CTPS(0.0, 4, 6, 1)
    Z = CTPS(0.0, 5, 6, 1)
    DELTA = CTPS(0.0, 6, 6, 1)
    rin = [X, PX, Y, PY, Z, DELTA]
    ringpass_TPSA!(RING, rin, 1, E0=3.5e9, m0=m_e)
    jaco = zeros(6, 6)
    for i in 1:6
        jaco[i, :] = rin[i].map[2:7]
    end
    return jaco
end
@btime g = TPSA_track_jacobian(0.001, 0.001)
# 137.778 ms (5412862 allocations: 462.68 MiB)