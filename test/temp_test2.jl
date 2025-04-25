using JuTrack
using Serialization
# using LaTeXStrings
# using Plots

E0 = 17.846262619763e9

ESR_crab = deserialize("src/demo/ESR/esr_main_linearquad.jls")

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

ESR_perturb = Q_perturb(ESR_crab)

idx = findelem(ESR_crab,CRABCAVITY)

xinit = [ESR_perturb[9].k1; ESR_perturb[13].k1; ESR_perturb[17].k1; ESR_perturb[23].k1; ESR_perturb[27].k1; ESR_perturb[31].k1; ESR_perturb[5537].k1]

zero_idx = [9,13,17,23,27,31,5537]
function get_phase14_zero(x)
    changed_idx = [9,13,17,23,27,31,5537]
    L = 0.25
    E0 = 17.846262619763e9
    new_Q1 = KQUAD(len=L, k1=x[1])
    new_Q2 = KQUAD(len=L, k1=x[2])
    new_Q3 = KQUAD(len=L, k1=x[3])
    new_Q4 = KQUAD(len=L, k1=x[4])
    new_Q5 = KQUAD(len=L, k1=x[5])
    new_Q6 = KQUAD(len=L, k1=x[6])
    new_Q7 = KQUAD(len=L, k1=x[7])
    changed_ele = [new_Q1, new_Q2, new_Q3, new_Q4, new_Q5, new_Q6, new_Q7]
    refpts = [i for i in 1:5550]
    # refpts = [35, 5533, length(ESR_crab)]
    twi = ADtwissring(ESR_perturb, 0.0, 0, refpts, changed_idx, changed_ele, E0=E0, m0=m_e)

    phase41 = twi[35].mux + twi[end].mux - twi[5533].mux
    return phase41 - 2*pi
end

function multi_val_op(x0, niter, step, RING)
    target = 0.01
    x0_vals = zeros(7, niter)
    goal_vals = []
    grad_vals = zeros(7, niter)
    zero_idx = [9,13,17,23,27,31,5537]
    g0 = get_phase14_zero(x0)
    for i in 1:niter
        # grads = gradient(Forward, get_phase14_zero, x0)
        # grad9 = grads[1]
        # grad13 = grads[2]
        # grad17 = grads[3]
        # grad23 = grads[4]
        # grad27 = grads[5]
        # grad31 = grads[6]
        # grad5537 = grads[7]
        grad9 = autodiff(Forward, get_phase14_zero, Duplicated(x0, [1.0,0.0,0.0,0.0,0.0,0.0,0.0]))
        grad13 = autodiff(Forward, get_phase14_zero, Duplicated(x0, [0.0,1.0,0.0,0.0,0.0,0.0,0.0]))
        grad17 = autodiff(Forward, get_phase14_zero, Duplicated(x0, [0.0,0.0,1.0,0.0,0.0,0.0,0.0]))
        grad23 = autodiff(Forward, get_phase14_zero, Duplicated(x0, [0.0,0.0,0.0,1.0,0.0,0.0,0.0]))
        grad27 = autodiff(Forward, get_phase14_zero, Duplicated(x0, [0.0,0.0,0.0,0.0,1.0,0.0,0.0]))
        grad31 = autodiff(Forward, get_phase14_zero, Duplicated(x0, [0.0,0.0,0.0,0.0,0.0,1.0,0.0]))
        grad5537 = autodiff(Forward, get_phase14_zero, Duplicated(x0, [0.0,0.0,0.0,0.0,0.0,0.0,1.0]))

        x0[1] -= step * grad9[1]
        x0[2] -= step * grad13[1]
        x0[3] -= step * grad17[1]
        x0[4] -= step * grad23[1]
        x0[5] -= step * grad27[1]
        x0[6] -= step * grad31[1]
        x0[7] -= step * grad5537[1]

        new_phase = get_phase14_zero(x0)
        println("init: ", g0, " now: ", new_phase, "at step ", i)
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

x0_vals, goal_vals, grad_vals = multi_val_op(xinit, 10, 5e-6, ESR_crab)

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
