include("../src/JuTrack.jl")
using .JuTrack
# using Enzyme
using Zygote

x = CTPS(0.0, 1, 6, 2)
px = CTPS(0.0, 2, 6, 2)
y = CTPS(0.0, 3, 6, 2)
py = CTPS(0.0, 4, 6, 2)
dp = 0.0

delta = CTPS(dp, 5, 6, 2)
z = CTPS(0.0, 6, 6, 2)

D1 = Drift("D1", 0.5)
Q1 = Quad("Q1", 1.0, -3.0, 0)
D2 = Drift("D2", 0.5)
Q2 = Quad("Q2", 1.0, 2.0, 0)
seq = [Q1, D1, Q2, D2]

rin = [x, px, y, py, delta, z]
rout = track(seq, rin)
println(rout[1].map)
println(rout[2].map)
println(rout[3].map)
println(rout[4].map)