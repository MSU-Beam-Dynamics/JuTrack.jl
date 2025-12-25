# Gradient descent optimization example using JuTrack.
# Objective:                minimize the sum of the following terms:  h21000, h10110, h30000, h10200, h10020 for the SPEAR3 lattice.
# Optimization variables:   SDM and SFM.
using Pkg
Pkg.activate("."); Pkg.instantiate() # change "." to your path of JuTrack.jl
using JuTrack
using Printf

set_tps_dim(2) # 2-D optimization

include("../src/demo/SPEAR3/spear3.jl")
RING = spear3()
RING = Number2TPSAD(RING)

SDM_index = findelem(RING, :name, "SDM") # 0.21m, -17 /m^-3
SFM_index = findelem(RING, :name, "SFM") # 0.21m, 15 /m^-3
idx_marker       = findelem(RING, MARKER)

for i in 1:length(SDM_index)
    RING[SDM_index[i]].k2 = -100.0 
end
for i in 1:length(SFM_index)
    RING[SFM_index[i]].k2 = 100.0 
end

function obj(dlist)
    tot_21000 = DTPSAD(0.0)
    tot_10110 = DTPSAD(0.0)
    tot_30000 = DTPSAD(0.0)
    tot_10200 = DTPSAD(0.0)
    tot_10020 = DTPSAD(0.0)
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
    for i in 1:length(SDM_index)
        RING[SDM_index[i]].k2 = x1
    end
    for i in 1:length(SFM_index)
        RING[SFM_index[i]].k2 = x2
    end

    dlist, _         = computeRDT(RING, idx_marker)
    return obj(dlist)               # <— scalar
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
        g, y = Gradient(f, [x1, x2], true)
        grad_x1 = g[1]
        grad_x2 = g[2]

        # Update variables
        x1 -= lr * grad_x1
        x2 -= lr * grad_x2
        push!(x1_his, x1)
        push!(x2_his, x2)
        push!(g1_his, grad_x1)
        push!(g2_his, grad_x2)
        push!(f_his, y)

        @printf("Iteration %d: x1=%.5f, x2=%.5f, grad1=%.5f, grad2=%.5f, f=%.5f\n", 
            iter, x1, x2, grad_x1, grad_x2, y)
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
x1_his, x2_his, g1_his, g2_his, f_his = gradient_descent(x1_init, x2_init, lr=1e-3, max_iter=200)

# # using Python's plotting library
# using PyCall
# @pyimport matplotlib.pyplot as plt
# Pkg.add("LaTeXStrings")
# using LaTeXStrings
# plt.figure(figsize=(8, 6))

# plt.subplot(3, 1, 1)
# plt.plot(x1_his, label="KSDM", linestyle="-", linewidth=1.5)
# plt.plot(x2_his, label="KSFM", linestyle="-", linewidth=1.5)
# plt.xlabel("Iteration", fontsize=15, fontfamily="Times New Roman")
# plt.ylabel(L"Strength ($\mathrm{m}^{-3}$)", fontsize=15, fontfamily="Times New Roman")
# plt.legend(prop=Dict("family"=>"Times New Roman"), loc="best", fontsize=18, frameon=false)
# plt.xticks(fontsize=15, fontfamily="Times New Roman")
# plt.yticks(fontsize=15, fontfamily="Times New Roman")
# plt.text(-0.15, 1.05, "a", transform=plt.gca().transAxes, fontsize=18, fontweight="bold", fontfamily="Times New Roman")

# plt.subplot(3, 1, 2)
# plt.plot(g1_his, label="grad w.r.t. KSDM",  linestyle="-", linewidth=1.5)
# plt.plot(g2_his, label="grad w.r.t. KSFM",  linestyle="-", linewidth=1.5)
# plt.xlabel("Iteration", fontsize=15, fontfamily="Times New Roman")
# plt.ylabel("Gradient", fontsize=15, fontfamily="Times New Roman")
# plt.xticks(fontsize=15, fontfamily="Times New Roman")
# plt.yticks(fontsize=15, fontfamily="Times New Roman")
# plt.legend(prop=Dict("family"=>"Times New Roman"), loc="best", fontsize=18, frameon=false)
# plt.text(-0.15, 1.05, "b", transform=plt.gca().transAxes, fontsize=18, fontweight="bold", fontfamily="Times New Roman")

# plt.subplot(3, 1, 3)
# plt.plot(f_his, linestyle="-", linewidth=1.5)
# plt.xlabel("Iteration", fontsize=15, fontfamily="Times New Roman")
# plt.ylabel("Objective", fontsize=15, fontfamily="Times New Roman")
# plt.xticks(fontsize=15, fontfamily="Times New Roman")
# plt.yticks(fontsize=15, fontfamily="Times New Roman")
# plt.text(-0.15, 1.05, "c", transform=plt.gca().transAxes, fontsize=18, fontweight="bold", fontfamily="Times New Roman")

# plt.tight_layout()
# plt.savefig("gradient_descent_optimization_SPEAR3.png", dpi=300)
# plt.show()


###############################################
# Legacy code. Using Enzyme for AD
###############################################

# include("../src/demo/SPEAR3/spear3.jl")
# RING = spear3()
# SDM_index = findelem(RING, :name, "SDM") # 0.21m, -17 /m^-3
# SFM_index = findelem(RING, :name, "SFM") # 0.21m, 15 /m^-3

# for i in 1:length(SDM_index)
#     RING[SDM_index[i]].k2 = -100.0 
# end
# for i in 1:length(SFM_index)
#     RING[SFM_index[i]].k2 = 100.0 
# end

# function obj(dlist)
#     tot_21000 = 0.0
#     tot_10110 = 0.0
#     tot_30000 = 0.0
#     tot_10200 = 0.0
#     tot_10020 = 0.0
#     for i in 1:length(dlist)
#         tot_21000 += dlist[i].h21000[1]
#         tot_10110 += dlist[i].h10110[1]
#         tot_30000 += dlist[i].h30000[1]
#         tot_10200 += dlist[i].h10200[1]
#         tot_10020 += dlist[i].h10020[1]
#     end
#     return tot_21000 + tot_10110 + tot_30000 + tot_10200 + tot_10020
# end

# function f(x1, x2, ring, sdm_id, sfm_id)
#     changed_ids      = vcat(sdm_id, sfm_id)
#     changed_elems    = vcat([KSEXT(len=0.21, k2=x1) for _ in sdm_id],
#                             [KSEXT(len=0.21, k2=x2) for _ in sfm_id])

#     idx_marker       = findelem(ring, MARKER)
#     dlist, _         = ADcomputeRDT(ring, idx_marker, changed_ids,
#                                     changed_elems; E0=3.5e9, m0=m_e)
#     return obj(dlist)               # <— scalar
# end

# function gradient_descent(x1_init, x2_init; lr=0.001, tol=1e-6, max_iter=100)
#     x1_his = [x1_init]
#     x2_his = [x2_init]
#     g1_his = []
#     g2_his = []
#     f_his = []
#     x1, x2 = x1_init, x2_init
#     for iter in 1:max_iter
#         # Calculate gradients
#         g1 = autodiff(ForwardWithPrimal, f, Duplicated(x1, 1.0), Const(x2), Const(RING), Const(SDM_index), Const(SFM_index))
#         g2 = autodiff(ForwardWithPrimal, f, Const(x1), Duplicated(x2, 1.0), Const(RING), Const(SDM_index), Const(SFM_index))
#         grad_x1 = g1[1]
#         grad_x2 = g2[1]

#         # Update variables
#         x1 -= lr * grad_x1
#         x2 -= lr * grad_x2
#         push!(x1_his, x1)
#         push!(x2_his, x2)
#         push!(g1_his, grad_x1)
#         push!(g2_his, grad_x2)
#         push!(f_his, g1[2])

#         @printf("Iteration %d: x1=%.5f, x2=%.5f, grad1=%.5f, grad2=%.5f, f=%.5f\n", 
#             iter, x1, x2, grad_x1, grad_x2, g1[2])
#         # Check for convergence
#         if abs(grad_x1) < tol && abs(grad_x2) < tol
#             println("Converged after $iter iterations.")
#             break
#         end
#     end
#     return x1_his, x2_his, g1_his, g2_his, f_his
# end

# # Initial guesses
# x1_init = -10.0
# x2_init = 10.0

# # Run gradient descent
# x1_his, x2_his, g1_his, g2_his, f_his = gradient_descent(x1_init, x2_init, lr=1e-4, max_iter=20)