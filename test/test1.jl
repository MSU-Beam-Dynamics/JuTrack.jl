include("../src/JuTrack.jl")
using .JuTrack
# using Enzyme
using Zygote

function f(k1, k2)
    x = CTPS(0.0, 1, 6, 2)
    px = CTPS(0.0, 2, 6, 2)
    y = CTPS(0.0, 3, 6, 2)
    py = CTPS(0.0, 4, 6, 2)
    dp = 0.0

    delta = CTPS(dp, 5, 6, 2)
    z = CTPS(0.0, 6, 6, 2)

    D1 = Drift("D1", 0.5)
    Q1 = Quad("Q1", 1.0, k1, 0)
    D2 = Drift("D2", 0.5)
    Q2 = Quad("Q2", 1.0, k2, 0)
    seq = [Q1, D1, Q2, D2]

    # rin = [x, px, y, py, delta, z]
    x1, px1, y1, py1, delta1, z1 = track(seq, x, px, y, py, delta, z)
    return x1.map
end
# println(f(-3.0, 2.0))
# grad = Zygote.jacobian(f, -3.0, 2.0)
# println(grad)


function f2(k1, k2)
    x = CTPS(0.0, 1, 6, 2)
    px = CTPS(0.0, 2, 6, 2)
    y = CTPS(0.0, 3, 6, 2)
    py = CTPS(0.0, 4, 6, 2)
    dp = 0.0

    delta = CTPS(dp, 5, 6, 2)
    z = CTPS(0.0, 6, 6, 2)

    D1 = Drift("D1", 0.5)
    Q1 = Quad("Q1", 1.0, k1, 0)
    D2 = Drift("D2", 0.5)
    Q2 = Quad("Q2", 1.0, k2, 0)
    seq = [Q1, D1, Q2, D2]

    twissin = EdwardsTengTwiss(betx=1.0, bety=1.0, alfx=0.0, alfy=0.0, dx=0.0, dy=0.0, dpx=0.0, 
                                dpy=0.0, mux=0.0, muy=0.0, R11=0.0, R12=0.0, R21=0.0, R22=0.0, mode=1)
    ss,names,ret = twissPropagate(twissin, seq, 0.0, 4)
    println(ret.betax+ret.betay)
    return ret.betax+ret.betay
end
grad_twiss = Zygote.gradient(f2, -3.0, 2.0)
println(grad_twiss)