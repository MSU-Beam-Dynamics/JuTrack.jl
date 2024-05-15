include("../src/JuTrack.jl")
using. JuTrack
using Enzyme
using BenchmarkTools
Enzyme.API.runtimeActivity!(true)

# an simple example that optimize the transfer matrix of a ring

D1 = DRIFT(len=1.0)
B1 = SBEND(len= 1.0, angle=pi/6.0)
Q1 = KQUAD(len=1.0, k1=-0.9325169994516977) # optimized k1 starting from -1.0
Q2 = KQUAD(len=1.0, k1=0.3)

cell = [D1, B1, Q1, B1, Q2, B1, Q1, B1, D1]
ring = [D1, B1, Q1, B1, Q2, B1, Q1, B1, D1, 
        D1, B1, Q1, B1, Q2, B1, Q1, B1, D1, 
        D1, B1, Q1, B1, Q2, B1, Q1, B1, D1]
# particles = [0.001 0.0001 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0]
# beam = Beam(particles)
# linepass!(ring, beam)
# println(beam.r)


function f1(k, ring)
    changed_idx = [3, 7, 12, 16, 21, 25]
    changed_ele = [KQUAD(len=1.0, k1=k), KQUAD(len=1.0, k1=k), 
                    KQUAD(len=1.0, k1=k), KQUAD(len=1.0, k1=k), 
                    KQUAD(len=1.0, k1=k), KQUAD(len=1.0, k1=k)]
    m66 = ADfindm66(ring, 0.0, 3, changed_idx, changed_ele)
    trace = m66[1, 1] + m66[2, 2]
    return trace
end
function f2(k, ring)
    changed_idx = [5, 14, 23]
    changed_ele = [KQUAD(len=1.0, k1=k), KQUAD(len=1.0, k1=k), KQUAD(len=1.0, k1=k)]
    m66 = ADfindm66(ring, 0.0, 3, changed_idx, changed_ele)
    trace = m66[1, 1] + m66[2, 2]
    return trace
end

grad1 = autodiff(Forward, f1, DuplicatedNoNeed, Duplicated(-1.0, 1.0),  Const(ring))
grad2 = autodiff(Forward, f2, DuplicatedNoNeed, Duplicated(1.0, 1.0),  Const(ring))
println(grad1)
println(grad2)

k1 = -1.0
k2 = 1.0
# optimize k1 to make the |trace| of the transfer matrix less than 2
x0 = [k1, k2]
niter = 100
step = 0.0001

x0_vals = Float64[]
trace_vals = Float64[]
grad_vals = Float64[]
for i in 1:niter
    trace0 = f1(x0[1], ring)
    grad = autodiff(Forward, f1, DuplicatedNoNeed, Duplicated(x0[1], 1.0),  Const(ring))
    if trace0 > 2
        x0[1] -= step * grad[1]
    else
        x0[1] += step * grad[1]
    end
    trace1 = f1(x0[1], ring)
    println("trace0: ", trace0, " trace1: ", trace1, "grad:", grad, "at step ", i)
    push!(x0_vals, x0[1])
    push!(trace_vals, trace1)
    push!(grad_vals, grad[1])
    if trace1 < 2.0 && trace1 > -2.0
        println("tuning finished at step ", i, " with trace1: ", trace1, " and target: ", 2.0)
        break
    end
end
println("k1: ", x0[1], " trace: ", f1(x0[1], ring))