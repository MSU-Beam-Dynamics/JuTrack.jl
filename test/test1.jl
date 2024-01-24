include("../src/JuTrack.jl")
include("ssrf_ring.jl")
using. JuTrack
using Enzyme
using BenchmarkTools

Enzyme.API.runtimeActivity!(true)
function f(x)
    SSRF = ssrf(x[1])
    particle = [0.001 0.0001 0.0005 0.0002 0.0 0.0; 0.001 0.0 0.0 0.0 0.0 0.0]
    Po = 6.849327668039155e+03
    # linepass!(SSRF, particle, Po, 0.0)
    ringpass!(SSRF, particle, 1, Po, 0.0)
    return particle[1,1]
end
# x = [-1.063770]
# println(f(x))
# @btime gradient(Forward, f, x)

# grad = gradient(Forward, f, x)
# println(grad)

function f_dup(x, y)
    SSRF = ssrf(x[1])
    particle = [0.001 0.0001 0.0005 0.0002 0.0 0.0; 0.001 0.0 0.0 0.0 0.0 0.0]
    Po = 6.849327668039155e+03
    # linepass!(SSRF, particle, Po, 0.0)
    ringpass!(SSRF, particle, 1, Po, 0.0)
    y[1] = particle[1,1]
    return nothing
end
# x = [-1.063770]
# y = [0.0]
# bx = [0.0]
# by = [1.0]
# Enzyme.autodiff(Reverse, f_dup, Duplicated(x, bx), Duplicated(y, by))
# println(y[1])

# function ff(xx)

#     # Sext = KSEXT(name="Q",len=xx[1], k2=xx[2], nSlices=10, synch_rad=0)
#     # Quad = KQUAD(name="Q",len=xx[1], k1=xx[2], nSlices=10, synch_rad=0)
#     dipole = CSBEND(name="Q",len=xx[1], angle=xx[2], nSlices=10, synch_rad=0)
#     particle = [0.001 0.0001 0.0005 0.0002 0.0 0.0; 0.001 0.0 0.0 0.0 0.0 0.0]
#     pass!(particle, dipole, 6.849327668039155e3, 0.0)
     
#     return particle[1,:]
# end
# x = [0.72, 0.07853981633974483]
# println(ff(x))
# println(ff([0.2, 1.555155/0.2*2]))
# grad = gradient(Forward, f, x)
# println(grad)