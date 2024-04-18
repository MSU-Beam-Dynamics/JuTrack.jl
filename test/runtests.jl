# # in Julia REPL, run ] to enter package mode
# # activate JuTrack
# # test JuTrack

# using JuTrack
# using Test
# using Enzyme

# function f(x)
#     k1 = x[1]
#     k2 = x[2]
#     x1 = CTPS(1.0, 1, 6, 3)
#     x2 = CTPS(2.0, 2, 6, 3)
#     y = k1*x1^2 + k2*x2^2
#     return y.map
# end

# @testset "JuTrack.jl" begin
# ctps = CTPS(Float64, 6, 3)
# ctps1 = CTPS(2.0, 6, 3)
# ctps2 = CTPS(3.0, 1, 6, 3)
# findindex(ctps2, [1, 0, 0, 0, 0, 0, 0])
# ind = findindex(ctps2, [1, 0, 0, 0, 0, 0, 0])
# ctps3 = ctps2*ctps2
# ctps4 = inv(ctps2)
# ctps5 = exp(ctps2)
# ctps6 = log(ctps2)
# ctps7 = sqrt(ctps2)
# ctps8 = ctps2^ 2
# ctps9 = sin(ctps2)
# ctps10 = cos(ctps2)
# ctps11 = asin(CTPS(0.5, 1, 6, 3))
# ctps12 = acos(CTPS(0.5, 1, 6, 3))
# ctps13 = sinh(ctps2)
# ctps14 = cosh(ctps2)
# println(ctps1.map)
# println(ctps2.map)
# println(ctps3.map)
# println(ctps4.map)
# println(ctps5.map)
# println(ctps6.map)
# println(ctps7.map)
# println(ctps8.map)
# println(ctps9.map)
# println(ctps10.map)
# println(ctps11.map)
# println(ctps12.map)
# println(ctps13.map)
# println(ctps14.map)

# k1 = 3.0
# k2 = 1.5
# y = f([k1, k2])
# println(y)
# # grad = jacobian(Forward, f, [k1, k2])
# # println(grad)

# end
