include("ssrf_ring.jl")
using JuTrack
using BenchmarkTools
using Plots
# Enzyme.API.runtimeActivity!(true)


SSRF = ssrf(-1.063770, 0)

function twiss_test(xx, ring=SSRF)
    # we don't suggest to create a long lattice inside the function. Use twiss_test(ring, xx) instead.
    # the ring can be set as Const in autodiff function
    # or use global variable 
    changed_idx = findelem(SSRF, :name, "QL1")
    changed_elems = [KQUAD(len=0.32, k1=xx) for i in 1:length(changed_idx)]
    twi0 = ADperiodicEdwardsTengTwiss(SSRF, 0.0, 1, changed_idx, changed_elems)
    # twi0 = periodicEdwardsTengTwiss(ring, 0.0, 1)
    
    return twi0.betax
end

function tuning_test(target)
    x0 = -1.063770*1.0
    niter = 20
    step = 0.0001

    x0_vals = Float64[]
    beta_vals = Float64[]
    grad_vals = Float64[]
    for i in 1:niter
        beta0 = twiss_test(x0)
        grad = autodiff(Forward, twiss_test, DuplicatedNoNeed, Duplicated(x0, 1.0), Const(SSRF))
        x0 -= step * grad[1]
        beta1 = twiss_test(x0)
        println("beta0: ", beta0, " beta1: ", beta1, " grad:", grad, " at step ", i)
        push!(x0_vals, x0)
        push!(beta_vals, beta1)
        push!(grad_vals, grad[1])
        if beta1 < target
            println("tuning finished at step ", i, " with beta1: ", beta1, " and target: ", target)
            break
        end
    end
    return x0_vals, beta_vals, grad_vals
end
x0_vals, beta_vals, grad_vals = tuning_test(7.0)

using LaTeXStrings
p1 = plot(1:length(x0_vals), x0_vals, title = L"Evolution\ of\ k_1", xlabel = L"Iteration", ylabel = L"Strength (m^{-2})", legend = false, line=:dash, marker=:circle)
p2 = plot(1:length(beta_vals), beta_vals, title = L"Evolution\ of\ \beta_x", xlabel = L"Iteration", ylabel = L"\beta_x(m)", legend = false, line=:dash, marker=:circle)
p3 = plot(1:length(grad_vals), grad_vals, title = L"Evolution\ of\ \frac{d \beta_x}{d k_1}", xlabel = L"Iteration", 
    ylabel = L"d \beta_x /d k_1", legend = false, line=:dash, marker=:circle)
plot(p1, p2, p3, layout = (3, 1), size=(800, 650))


