include("../src/JuTPSA.jl")
using .JuTPSA
using Test
# using Enzyme
using Zygote

function f(k1, k2)
    ctps1 = CTPS(k1, 1, 3, 3) 
    ctps2 = k2*CTPS(1.0, 2, 3, 3)
    ctps3 = ctps1 * sinh(ctps2)
    # result = k2*(pow(ctps2, k1))
    result = derivative(ctps3, 1, 1)
    map = result.map
    nterm = result.terms
    return map
end
k1 = 3.0
k2 = 1.0
println(f(k1 ,k2))
# grad_g = Zygote.jacobian(f, k1)
grad_g = Zygote.jacobian(f, k1, k2)

println(grad_g)

# function f(x, y)
#     k1 = x[1] 
#     ctps1 = CTPS(1.0, 1, 3, 3) 
#     result = integrate(ctps1, 1, k1)
#     y[1] = result.map[1]
#     return Nothing
# end
# x  = [3.0]
# bx = [0.0]
# y  = [0.0]
# by = [1.0]
# Enzyme.autodiff(Reverse, FODO_track_result, Duplicated(x, bx), Duplicated(y, by))
