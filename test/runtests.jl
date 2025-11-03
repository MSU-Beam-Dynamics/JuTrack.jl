# in Julia REPL, type ] to enter package mode
# activate JuTrack
# test JuTrack

using JuTrack
using Test

function f(x)
    k1 = x[1]
    k2 = x[2]
    x1 = CTPS(1.0, 1, 6, 3)
    x2 = CTPS(2.0, 2, 6, 3)
    y = k1*x1^2 + k2*x2^2
    return y.map
end

@testset "JuTrack.jl" begin
ctps = CTPS(Float64, 6, 3)
ctps1 = CTPS(2.0, 6, 3)
ctps2 = CTPS(3.0, 1, 6, 3)
findindex(ctps2, [1, 0, 0, 0, 0, 0, 0])
ind = findindex(ctps2, [1, 0, 0, 0, 0, 0, 0])
ctps3 = ctps2*ctps2
ctps4 = 1.0 / ctps2
ctps5 = exp(ctps2)
ctps6 = log(ctps2)
ctps7 = sqrt(ctps2)
ctps8 = ctps2^ 2
ctps9 = sin(ctps2)
ctps10 = cos(ctps2)
ctps11 = asin(CTPS(0.5, 1, 6, 3))
ctps12 = acos(CTPS(0.5, 1, 6, 3))
ctps13 = sinh(ctps2)
ctps14 = cosh(ctps2)
println(ctps1.map)
println(ctps2.map)
println(ctps3.map)
println(ctps4.map)
println(ctps5.map)
println(ctps6.map)
println(ctps7.map)
println(ctps8.map)
println(ctps9.map)
println(ctps10.map)
println(ctps11.map)
println(ctps12.map)
println(ctps13.map)
println(ctps14.map)

k1 = 3.0
k2 = 1.5
y = f([k1, k2])
println(y)
grad = jacobian(Forward, f, [k1, k2])
println(grad)

# particle tracking
beam = Beam([0.001 0.0 0.0 0.0 0.0 0.0])
D1 = DRIFT(len=1.0)
Q1 = QUAD(k1=1.2, len=0.5)
Q2 = QUAD(k1=-1.2, len=0.5)
line = [D1, Q1, D1, Q2, D1]
linepass!(line, beam)
println(beam.r)

# fast tpsa
set_tps_dim(6)
x = DTPSAD(0.001, 1)
xp = DTPSAD(0.0, 2)
y = DTPSAD(0.0, 3)
yp = DTPSAD(0.0, 4)
z = DTPSAD(0.0, 5)
dp = DTPSAD(0.0, 6)
beam_tpsa = Beam([x xp y yp z dp])
D1 = DRIFT(len=DTPSAD(1.0))
Q1 = QUAD(k1=DTPSAD(1.2), len=DTPSAD(0.5))
Q2 = QUAD(k1=DTPSAD(-1.2), len=DTPSAD(0.5))
line_ = [D1, Q1, D1, Q2, D1]
linepass!(line_, beam_tpsa)
println(beam_tpsa.r)
end
