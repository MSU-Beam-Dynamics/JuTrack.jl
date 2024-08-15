using JuTrack
using Serialization
using Plots
using Printf

RING = deserialize("src/demo/SPEAR3/spear3.jls")

# index = findelem(RING, MARKER)
# dlist, s = computeRDT(RING, index)
# g_terms = [:h21000, :h10110, :h30000, :h10200,  :h10020]

# function plot_RDT(dlist, s, name::Symbol)
#     y = [getfield(dlist[i], name)[1] for i in 1:length(dlist)]
#     plot!(s, y, label=String(name))
# end

# plot( title="RDT", xlabel="s (m)", ylabel="RDT", legend=:topleft)
# for i in 1:length(g_terms)
#     plot_RDT(dlist, s, g_terms[i])
# end
# display(plot!())


SDM_index = findelem(RING, :name, "SDM") # 0.21m, -17 /m^-3
SFM_index = findelem(RING, :name, "SFM") # 0.21m, 15 /m^-3

for i in 1:length(SDM_index)
    RING[SDM_index[i]].PolynomB[3] = -10.0 
end
for i in 1:length(SFM_index)
    RING[SFM_index[i]].PolynomB[3] = 10.0 
end

function obj(dlist)
    tot_21000 = 0.0
    tot_10110 = 0.0
    tot_30000 = 0.0
    tot_10200 = 0.0
    tot_10020 = 0.0
    for i in 1:length(dlist)
        tot_21000 += dlist[i].h21000[1]
        tot_10110 += dlist[i].h10110[1]
        tot_30000 += dlist[i].h30000[1]
        tot_10200 += dlist[i].h10200[1]
        tot_10020 += dlist[i].h10020[1]
    end
    return tot_21000 + tot_10110 + tot_30000 + tot_10200 + tot_10020
end
function f(x1, x2)
    changed_id1 = findelem(RING, :name, "SDM")
    changed_id2 = findelem(RING, :name, "SFM")
    changed_ids = vcat(changed_id1, changed_id2)
    changed_elems = [i in changed_id1 ? KSEXT(len=0.21, k2=x1) : KSEXT(len=0.21, k2=x2) for i in 1:length(RING)]

    index = findelem(RING, MARKER)
    dlist, s = ADcomputeRDT(RING, index, changed_ids, changed_elems)
    return obj(dlist)
end


# g1 = autodiff(Forward, f, Duplicated, Duplicated(-10.0, 1.0), Const(10.0))
# g2 = autodiff(Forward, f, Duplicated, Const(-10.0), Duplicated(10.0, 1.0))
# g = autodiff(Forward, f, BatchDuplicated, BatchDuplicated(-17.0, (1.0, 0.0)), BatchDuplicated(15.0, (0.0, 1.0))) # problematic

function gradient_descent(f, x1_init, x2_init; lr=0.001, tol=1e-6, max_iter=100)
    x1, x2 = x1_init, x2_init
    for iter in 1:max_iter
        # Calculate gradients
        g1 = autodiff(Forward, fSDM, Duplicated, Duplicated(x1, 1.0), Const(x2))
        g2 = autodiff(Forward, fSDM, Duplicated, Const(x1), Duplicated(x2, 1.0))
        grad_x1 = g1[2]
        grad_x2 = g2[2]

        # Update variables
        x1 -= lr * grad_x1
        x2 -= lr * grad_x2

        @printf("Iteration %d: x1=%.5f, x2=%.5f, grad1=%.5f, grad2=%.5f, f=%.5f\n", 
            iter, x1, x2, grad_x1, grad_x2, g1[1])
        # Check for convergence
        if abs(grad_x1) < tol && abs(grad_x2) < tol
            println("Converged after $iter iterations.")
            break
        end
    end
    return x1, x2
end

# Initial guesses
x1_init = -10.0
x2_init = 10.0

# Run gradient descent
optimal_x1, optimal_x2 = gradient_descent(f, x1_init, x2_init)

println("Optimized x1: $optimal_x1")
println("Optimized x2: $optimal_x2")

