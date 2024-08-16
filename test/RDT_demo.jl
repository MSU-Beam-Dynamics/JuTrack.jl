using JuTrack
using Serialization
using Plots
using Printf

# this example calculate the sum of third order resonance terms for the SPEAR3 lattice
# the objective is to minimize the sum of the following terms:  h21000, h10110, h30000, h10200, h10020
# the optimization variables are the sextupole strengths of SDM and SFM
# the optimization is done by gradient descent. The gradient is calculated by autodiff @Enzyme
RING = deserialize("src/demo/SPEAR3/spear3.jls")

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
    changed_elems = [changed_ids[i] in changed_id1 ? KSEXT(len=0.21, k2=x1) : KSEXT(len=0.21, k2=x2) for i in 1:length(changed_ids)]

    index = findelem(RING, MARKER)
    dlist, s = ADcomputeRDT(RING, index, changed_ids, changed_elems)
    return obj(dlist)
end


# g1 = autodiff(Forward, f, Duplicated, Duplicated(-10.0, 1.0), Const(10.0))
# g2 = autodiff(Forward, f, Duplicated, Const(-10.0), Duplicated(10.0, 1.0))
# g = autodiff(Forward, f, BatchDuplicated, BatchDuplicated(-17.0, (1.0, 0.0)), BatchDuplicated(15.0, (0.0, 1.0))) # problematic

function gradient_descent(x1_init, x2_init; lr=0.001, tol=1e-6, max_iter=100)
    x1_his = [x1_init]
    x2_his = [x2_init]
    g1_his = []
    g2_his = []
    f_his = []
    x1, x2 = x1_init, x2_init
    for iter in 1:max_iter
        # Calculate gradients
        g1 = autodiff(Forward, f, Duplicated, Duplicated(x1, 1.0), Const(x2))
        g2 = autodiff(Forward, f, Duplicated, Const(x1), Duplicated(x2, 1.0))
        grad_x1 = g1[2]
        grad_x2 = g2[2]

        # Update variables
        x1 -= lr * grad_x1
        x2 -= lr * grad_x2
        push!(x1_his, x1)
        push!(x2_his, x2)
        push!(g1_his, grad_x1)
        push!(g2_his, grad_x2)
        push!(f_his, g1[1])

        @printf("Iteration %d: x1=%.5f, x2=%.5f, grad1=%.5f, grad2=%.5f, f=%.5f\n", 
            iter, x1, x2, grad_x1, grad_x2, g1[1])
        # Check for convergence
        if abs(grad_x1) < tol && abs(grad_x2) < tol
            println("Converged after $iter iterations.")
            break
        end
    end
    return x1_his, x2_his, g1_his, g2_his, f_his
end

# Initial guesses
x1_init = -10.0
x2_init = 10.0

# Run gradient descent
x1_his, x2_his, g1_his, g2_his, f_his = gradient_descent(x1_init, x2_init, lr=0.01, max_iter=1000)

# Plot the results
p1 = plot(0:length(x1_his)-1, x1_his, label="x1", xlabel="Iteration", ylabel="Value", title="Gradient Descent")
plot!(0:length(x2_his)-1, x2_his, label="x2")
p2 = plot(g1_his, label="grad1", xlabel="Iteration", ylabel="Gradient")
plot!(g2_his, label="grad2")
p3 = plot(f_his, label="f", xlabel="Iteration", ylabel="Objective")
plot(p1, p2, p3, layout=(3,1), legend=:topleft)

# using PyCall
# @pyimport matplotlib.pyplot as plt
# plt.plot(figsize=(15, 10))
# plt.subplot(3, 1, 1)
# plt.plot(0:length(x1_his)-1, x1_his, label="x1")
# plt.plot(0:length(x2_his)-1, x2_his, label="x2")
# plt.xlabel("Iteration")
# plt.ylabel("Value")
# plt.legend()
# plt.subplot(3, 1, 2)
# plt.plot(g1_his, label="grad1")
# plt.plot(g2_his, label="grad2")
# plt.xlabel("Iteration")
# plt.ylabel("Gradient")
# plt.legend()
# plt.subplot(3, 1, 3)
# plt.plot(f_his, label="f")  
# plt.xlabel("Iteration")
# plt.ylabel("Objective")
# plt.legend()
# plt.tight_layout()
# plt.show()

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
