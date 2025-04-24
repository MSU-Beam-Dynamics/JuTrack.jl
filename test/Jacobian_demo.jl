include("../src/demo/SPEAR3/spear3.jl")
using JuTrack
using Serialization
using BenchmarkTools

RING = spear3()


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
@btime begin # run 10000 times or reach max time
    g = TPSA_track_jacobian(0.001, 0.001)
end
# 132.694 ms (4920663 allocations: 425.79 MiB)

function AD_jacobian(x)
    beam = Beam([x[1] x[2] x[3] x[4] x[5] x[6]], energy=3.5e9)
    ringpass!(RING, beam, 1)
    return beam.r
end
@btime begin # run 10000 times or reach max time
    g = jacobian(Forward, AD_jacobian, [0.001, 0.0, 0.001, 0.0, 0.0, 0.0])
end
# 2.847 ms (20576 allocations: 1.61 MiB)