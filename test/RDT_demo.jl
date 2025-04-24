using JuTrack
using Serialization
# using Plots
using Printf

# this example calculate the sum of third order resonance terms for the SPEAR3 lattice
# the objective is to minimize the sum of the following terms:  h21000, h10110, h30000, h10200, h10020
# the optimization variables are the sextupole strengths of SDM and SFM. 
# The code may need further optimization for large memory allocations.
include("../src/demo/SPEAR3/spear3.jl")
RING = spear3()
SDM_index = findelem(RING, :name, "SDM") # 0.21m, -17 /m^-3
SFM_index = findelem(RING, :name, "SFM") # 0.21m, 15 /m^-3

for i in 1:length(SDM_index)
    RING[SDM_index[i]].PolynomB[3] = -100.0 
end
for i in 1:length(SFM_index)
    RING[SFM_index[i]].PolynomB[3] = 100.0 
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
    # Find elements with names "SDM" and "SFM" in the RING
    # changed_id1 = findelem(RING, :name, "SDM")
    # changed_id2 = findelem(RING, :name, "SFM")
    
    # Combine the found element indices
    changed_ids = vcat(SDM_index, SFM_index)
    # Create new elements based on the found indices
    changed_elems1 = [KSEXT(len=0.21, k2=x1) for id in SDM_index]
    changed_elems2 = [KSEXT(len=0.21, k2=x2) for id in SFM_index]
    
    changed_elems = [changed_elems1..., changed_elems2...] #vcat(changed_elems1, changed_elems2)
    
    # Find the index of the MARKER element in the RING
    index = findelem(RING, MARKER)
    
    # Compute the RDT
    dlist, s = ADcomputeRDT(RING, index, changed_ids, changed_elems, E0=3.5e9, m0=m_e)
    
    # Return the objective function value
    return obj(dlist)
end


function gradient_descent(x1_init, x2_init; lr=0.001, tol=1e-6, max_iter=100)
    x1_his = [x1_init]
    x2_his = [x2_init]
    g1_his = []
    g2_his = []
    f_his = []
    x1, x2 = x1_init, x2_init
    for iter in 1:max_iter
        # Calculate gradients
        g1 = autodiff(ForwardWithPrimal, f, Duplicated(x1, 1.0), Const(x2))
        g2 = autodiff(ForwardWithPrimal, f, Const(x1), Duplicated(x2, 1.0))
        grad_x1 = g1[1]
        grad_x2 = g2[1]

        # Update variables
        x1 -= lr * grad_x1
        x2 -= lr * grad_x2
        push!(x1_his, x1)
        push!(x2_his, x2)
        push!(g1_his, grad_x1)
        push!(g2_his, grad_x2)
        push!(f_his, g1[2])

        @printf("Iteration %d: x1=%.5f, x2=%.5f, grad1=%.5f, grad2=%.5f, f=%.5f\n", 
            iter, x1, x2, grad_x1, grad_x2, g1[2])
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
x1_his, x2_his, g1_his, g2_his, f_his = gradient_descent(x1_init, x2_init, lr=1e-7, max_iter=20)

# using PyCall
# @pyimport matplotlib.pyplot as plt
# using LaTeXStrings
# default(;fontfamily="Times Roman")
# plt.figure(figsize=(8, 6))

# plt.subplot(3, 1, 1)
# plt.plot(x1_his[1:600], label="KSDM", linestyle="-", linewidth=1.5)
# plt.plot(x2_his[1:600], label="KSFM", linestyle="-", linewidth=1.5)
# plt.xlabel("Iteration", fontsize=15, fontfamily="Times New Roman")
# plt.ylabel(L"Strength ($\mathrm{m}^{-3}$)", fontsize=15, fontfamily="Times New Roman")
# plt.legend(prop=Dict("family"=>"Times New Roman"), loc="best", fontsize=18, frameon=false)
# plt.xticks(fontsize=15, fontfamily="Times New Roman")
# plt.yticks(fontsize=15, fontfamily="Times New Roman")
# plt.text(-0.15, 1.05, "a", transform=plt.gca().transAxes, fontsize=18, fontweight="bold", fontfamily="Times New Roman")

# plt.subplot(3, 1, 2)
# plt.plot(g1_his[1:600], label="grad w.r.t. KSDM",  linestyle="-", linewidth=1.5)
# plt.plot(g2_his[1:600], label="grad w.r.t. KSFM",  linestyle="-", linewidth=1.5)
# plt.xlabel("Iteration", fontsize=15, fontfamily="Times New Roman")
# plt.ylabel("Gradient", fontsize=15, fontfamily="Times New Roman")
# plt.xticks(fontsize=15, fontfamily="Times New Roman")
# plt.yticks(fontsize=15, fontfamily="Times New Roman")
# plt.legend(prop=Dict("family"=>"Times New Roman"), loc="best", fontsize=18, frameon=false)
# plt.text(-0.15, 1.05, "b", transform=plt.gca().transAxes, fontsize=18, fontweight="bold", fontfamily="Times New Roman")

# plt.subplot(3, 1, 3)
# plt.plot(f_his[1:600], linestyle="-", linewidth=1.5)
# plt.xlabel("Iteration", fontsize=15, fontfamily="Times New Roman")
# plt.ylabel("Objective", fontsize=15, fontfamily="Times New Roman")
# plt.xticks(fontsize=15, fontfamily="Times New Roman")
# plt.yticks(fontsize=15, fontfamily="Times New Roman")
# plt.text(-0.15, 1.05, "c", transform=plt.gca().transAxes, fontsize=18, fontweight="bold", fontfamily="Times New Roman")

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
