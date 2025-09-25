# Example of using JuTrack to tune the ESR lattice for a specific phase advance between two crab cavities
# Warning: This code is for demonstration purposes. Simply tuning the quadrupoles may result in unstable solutions.
using Pkg
Pkg.activate("."); Pkg.instantiate() # change "." to your path of JuTrack.jl
using JuTrack
using Serialization
using Printf
set_tps_dim(7)  # 7 variables for the TPSA
E0 = 17.846262619763e9

# Load the ESR main lattice 
ESR_crab = deserialize("src/demo/ESR/esr_main.jls")
ESR_crab = Number2TPSAD(ESR_crab) # Convert to TPSA format

function Q_perturb(ESR_crab)
    for i in eachindex(ESR_crab)
        if ESR_crab[i] isa KQUAD
            k1 = ESR_crab[i].k1
            k1 = k1 * (1 + 0.001 * randn_approx())
            new_KQ = KQUAD(name=ESR_crab[i].name, len=ESR_crab[i].len, k1=k1)
            ESR_crab[i] = new_KQ
        end
    end
    return ESR_crab
end

ESR_perturb = Q_perturb(ESR_crab) # Perturbed ESR lattice
idx = findelem(ESR_crab, CRABCAVITY)

function get_phase14_zero(x1, x2, x3, x4, x5, x6, x7)
    X = [x1, x2, x3, x4, x5, x6, x7]
    changed_idx = [9,13,17,23,27,31,5537] # indices of the quadrupoles to be changed
    E0 = 17.846262619763e9

    refpts = [i for i in 1:length(ESR_crab)]
    for i in 1:length(changed_idx)
        ESR_perturb[changed_idx[i]].k1 = X[i]
    end
    # refpts = [35, 5533, length(ESR_crab)]
    twi = twissring(ESR_perturb, 0.0, 0, refpts, E0=E0, m0=m_e)

    phase41 = twi[35].mux + twi[end].mux - twi[5533].mux
    return phase41 - 2*pi
end

function multi_val_op(x0, niter, step, RING)
    target = 0.015
    x0_vals = zeros(7, niter)
    goal_vals = []
    grad_vals = zeros(7, niter)
    g0 = get_phase14_zero(x0[1], x0[2], x0[3], x0[4], x0[5], x0[6], x0[7])
    for i in 1:niter
        grads, new_phase = Gradient(get_phase14_zero, x0, true)
        grad9 = grads[1]
        grad13 = grads[2]
        grad17 = grads[3]
        grad23 = grads[4]
        grad27 = grads[5]
        grad31 = grads[6]
        grad5537 = grads[7]

        x0[1] -= step * grad9[1]
        x0[2] -= step * grad13[1]
        x0[3] -= step * grad17[1]
        x0[4] -= step * grad23[1]
        x0[5] -= step * grad27[1]
        x0[6] -= step * grad31[1]
        x0[7] -= step * grad5537[1]

        println(@sprintf("init: %.6f now: %.6f at step %d, 
        k1 : %.6f, k2 : %.6f, k3 : %.6f, k4 : %.6f, k5 : %.6f, k6 : %.6f, k7 : %.6f", g0.val, new_phase, i,
            x0[1], x0[2], x0[3], x0[4], x0[5], x0[6], x0[7]))
        x0_vals[:, i] = x0
        push!(goal_vals, new_phase)
        grad_vals[:, i] = [grad9[1], grad13[1], grad17[1], grad23[1], grad27[1], grad31[1], grad5537[1]]
        if new_phase < target 
            println("tuning finished at step ", i)
            break
        end
    end
    return x0_vals, goal_vals, grad_vals
end

xinit = [-1e-6, 1e-6, -1e-6, 1e-6, -1e-6, 1e-6, -1e-6] # Initial guess for the quadrupole strengths
x0_vals, goal_vals, grad_vals = multi_val_op(xinit, 100, 1e-4, ESR_crab)

# using LaTeXStrings
# plot_steps = length(goal_vals)
# p1 = plot(1:plot_steps, x0_vals[1, 1:plot_steps], title = L"Evolution\ of\ k", xlabel = L"Iterations", ylabel = L"Strength (m^{-1})", label=L"k_1", line=:dash, marker=:circle,framestyle=:box)
# for i in 2:7
#     label_str = "k_{$i}"  
#     full_label = latexstring(label_str)
#     plot!(1:plot_steps, x0_vals[i, 1:plot_steps], label=full_label, line=:dash, marker=:circle)
# end
# p2 = plot(1:plot_steps, goal_vals[1:plot_steps], title = L"Evolution\ of\ \Delta \phi", xlabel = L"Iterations", 
#     ylabel = L"phase\ advance(rad)", legend = false, line=:dash, marker=:circle,framestyle=:box)
# p3 = plot(1:plot_steps, grad_vals[1, 1:plot_steps], title = L"Evolution\ of\ gradient", xlabel = L"Iterations", 
#     ylabel = L"\partial \frac{\Delta \phi}{\partial k}", label=L"k_1", line=:dash, marker=:circle,framestyle=:box)
# for i in 2:7
#     label_str = "k_{$i}"  
#     full_label = latexstring(label_str)
#     plot!(1:plot_steps, grad_vals[i, 1:plot_steps], label=full_label, line=:dash, marker=:circle)
# end
# plot(p1, p2, p3, layout = (3, 1), size=(800, 650))
# # savefig("plot1.png")


##########################################################
# Legacy code of using Enzyme for gradient calculation
##########################################################
# E0 = 17.846262619763e9

# ESR_crab = deserialize("src/demo/ESR/esr_main_linearquad.jls")

# function Q_perturb(ESR_crab)
#     for i in eachindex(ESR_crab)
#         if ESR_crab[i] isa KQUAD
#             k1 = ESR_crab[i].k1
#             k1 = k1 * (1 + 0.001 * randn_approx())
#             new_KQ = KQUAD(name=ESR_crab[i].name, len=ESR_crab[i].len, k1=k1)
#             ESR_crab[i] = new_KQ
#         end
#     end
#     return ESR_crab
# end

# ESR_perturb = Q_perturb(ESR_crab) # Perturbed ESR lattice

# idx = findelem(ESR_crab,CRABCAVITY)

# xinit = [-1e-6, 1e-6, -1e-6, 1e-6, -1e-6, 1e-6, -1e-6]

# zero_idx = [9,13,17,23,27,31,5537]
# function get_phase14_zero(x, ESR_perturb)
#     changed_idx = [9,13,17,23,27,31,5537]
#     L = 0.25
#     E0 = 17.846262619763e9
#     new_Q1 = QUAD(len=L, k1=x[1])
#     new_Q2 = QUAD(len=L, k1=x[2])
#     new_Q3 = QUAD(len=L, k1=x[3])
#     new_Q4 = QUAD(len=L, k1=x[4])
#     new_Q5 = QUAD(len=L, k1=x[5])
#     new_Q6 = QUAD(len=L, k1=x[6])
#     new_Q7 = QUAD(len=L, k1=x[7])
#     changed_ele = [new_Q1, new_Q2, new_Q3, new_Q4, new_Q5, new_Q6, new_Q7]
#     refpts = [i for i in 1:5550]
#     # refpts = [35, 5533, length(ESR_crab)]
#     twi = ADtwissring(ESR_perturb, 0.0, 0, refpts, changed_idx, changed_ele, E0=E0, m0=m_e)

#     phase41 = twi[35].mux + twi[end].mux - twi[5533].mux
#     return phase41 - 2*pi
# end

# function multi_val_op(x0, niter, step, RING)
#     target = 0.01
#     x0_vals = zeros(7, niter)
#     goal_vals = []
#     grad_vals = zeros(7, niter)
#     zero_idx = [9,13,17,23,27,31,5537]
#     g0 = get_phase14_zero(x0, ESR_perturb)
#     for i in 1:niter
#         grads = jacobian(Forward, get_phase14_zero, x0, Const(ESR_perturb))

#         x0[1] -= step * grads[1][1]
#         x0[2] -= step * grads[1][2]
#         x0[3] -= step * grads[1][3]
#         x0[4] -= step * grads[1][4]
#         x0[5] -= step * grads[1][5]
#         x0[6] -= step * grads[1][6]
#         x0[7] -= step * grads[1][7]

#         new_phase = get_phase14_zero(x0, ESR_perturb)
#         println("init: ", g0, " now: ", new_phase, "at step ", i)
#         x0_vals[:, i] = x0
#         push!(goal_vals, new_phase)
#         grad_vals[:, i] = [grads[1][1], grads[1][2], grads[1][3], grads[1][4], grads[1][5], grads[1][6], grads[1][7]]
#         if new_phase < target
#             println("tuning finished at step ", i)
#             break
#         end
#     end
#     return x0_vals, goal_vals, grad_vals
# end

# x0_vals, goal_vals, grad_vals = multi_val_op(xinit, 10, 5e-6, ESR_crab)