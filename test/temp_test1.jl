include("../src/JuTrack.jl")
using. JuTrack
# Enzyme.API.runtimeActivity!(true)

# an simple example that optimize the transfer matrix of a ring
Dr = DRIFT(name="Dr", len=0.5)
HalfDr = DRIFT(name="HalfDr", len=0.25)
p2Dr = DRIFT(name="p2Dr", len=0.2)
SF = KSEXT(name="SF", len=0.1, k2=1.0)
SD = KSEXT(name="SD", len=0.1, k2=-1.0)
B1 = SBEND(name="B", len= 1.0, angle=2*pi/40.0)
Q1 = KQUAD(name="Q1", len=0.5, k1=1.2) # optimized k1 starting from -1.0
Q2 = KQUAD(name="Q2", len=0.5, k1=-1.2)

cell = [HalfDr, B1, p2Dr, SF, p2Dr, Q1, Dr, B1, p2Dr, SD, p2Dr, Q2, HalfDr]
ring = [cell..., cell..., cell..., cell..., cell..., cell..., cell..., cell..., cell..., cell...,
        cell..., cell..., cell..., cell..., cell..., cell..., cell..., cell..., cell..., cell...]



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