include("../src/JuTrack.jl")
using. JuTrack
using Serialization
using Enzyme
using LaTeXStrings
using Plots

Enzyme.API.runtimeActivity!(true)


E0 = 17.846262619763e9


ESR_crab = deserialize("test/esr_main_linearquad.jls")
# ESR_nocrab = deserialize("test/esr_main_rad_craboff.jls")

function get_phase14(x3, RING)
    # change the 3rd quad, optimize phase advance between CC1-35 and CC4-5533
    changed_idx = [3]
    new_D1 = KQUAD(len=RING[3].len, k1=x3)
    changed_ele = [new_D1]
    refpts = [i for i in 1:length(RING)]
    twi = ADtwissring(RING, 0.0, 1, refpts, changed_idx, changed_ele)

    phase41 = twi[35].dmux + twi[end].dmux - twi[5533].dmux

    return phase41 - 2*pi
end

function get_phase23(x, RING)
    # change the 3rd quad, optimize phase advance between CC1-35 and CC4-5533
    changed_idx = [3]
    new_D1 = KQUAD(len=RING[3].len, k1=x)
    changed_ele = [new_D1]
    refpts = [i for i in 1:length(RING)]
    twi = ADtwissring(RING, 0.0, 1, refpts, changed_idx, changed_ele)

    phase23 = twi[952].dmux - twi[916].dmux

    return phase23 - 2*pi
end

phi14 = get_phase14(ESR_crab[3].k1, ESR_crab)
phi23 = get_phase23(ESR_crab[3].k1, ESR_crab)

function Q_perturb(ESR_crab)
    for i in eachindex(ESR_crab)
        if ESR_crab[i] isa KQUAD
            k1 = ESR_crab[i].k1
            k1 = k1 * (1 + 0.001 * randn())
            new_KQ = KQUAD(name=ESR_crab[i].name, len=ESR_crab[i].len, k1=k1)
            ESR_crab[i] = new_KQ
        end
    end
    return ESR_crab
end
function optics(Twi)
    beta = zeros(length(Twi), 2)
    beta[:, 1] = [Twi[i].betax for i in eachindex(Twi)]
    beta[:, 2] = [Twi[i].betay for i in eachindex(Twi)]
    alpha = zeros(length(Twi), 2)
    alpha[:, 1] = [Twi[i].alphax for i in eachindex(Twi)]
        alpha[:, 2] = [Twi[i].alphay for i in eachindex(Twi)]
    gamma = zeros(length(Twi), 2)
    gamma[:, 1] = [Twi[i].gammax for i in eachindex(Twi)]
    gamma[:, 2] = [Twi[i].gammay for i in eachindex(Twi)]
    mu = zeros(length(Twi), 2)
    mu[:, 1] = [Twi[i].dmux for i in eachindex(Twi)]
    mu[:, 2] = [Twi[i].dmuy for i in eachindex(Twi)]
    dp = zeros(length(Twi), 4)
    dp[:, 1] = [Twi[i].dx for i in eachindex(Twi)]
    dp[:, 2] = [Twi[i].dy for i in eachindex(Twi)]
    dp[:, 3] = [Twi[i].dpx for i in eachindex(Twi)]
    dp[:, 4] = [Twi[i].dpy for i in eachindex(Twi)]
    return beta, alpha, gamma, mu, dp
end
ESR_perturb = Q_perturb(ESR_crab)

idx = findelem(ESR_crab,CRABCAVITY)

xinit = [ESR_perturb[9].k1; ESR_perturb[13].k1; ESR_perturb[17].k1; ESR_perturb[23].k1; ESR_perturb[27].k1; ESR_perturb[31].k1; ESR_perturb[5537].k1]

zero_idx = [9,13,17,23,27,31,5537]
function get_phase14_zero(x, RING)
    changed_idx = [9,13,17,23,27,31,5537]
    new_D1 = KQUAD(len=RING[9].len, k1=x[1])
    new_D2 = KQUAD(len=RING[13].len, k1=x[2])
    new_D3 = KQUAD(len=RING[17].len, k1=x[3])
    new_D4 = KQUAD(len=RING[23].len, k1=x[4])
    new_D5 = KQUAD(len=RING[27].len, k1=x[5])
    new_D6 = KQUAD(len=RING[31].len, k1=x[6])
    new_D7 = KQUAD(len=RING[5537].len, k1=x[7])
    changed_ele = [new_D1, new_D2, new_D3, new_D4, new_D5, new_D6, new_D7]
    refpts = [i for i in 1:length(RING)]
    twi = ADtwissring(RING, 0.0, 1, refpts, changed_idx, changed_ele)

    phase41 = twi[35].dmux + twi[end].dmux - twi[5533].dmux

    return phase41 - 2*pi
end

function get_phase23_zero(x, RING)
    changed_idx = [9,13,17,23,27,31,5537]
    new_D1 = KQUAD(len=RING[9].len, k1=x[1])
    new_D2 = KQUAD(len=RING[13].len, k1=x[2])
    new_D3 = KQUAD(len=RING[17].len, k1=x[3])
    new_D4 = KQUAD(len=RING[23].len, k1=x[4])
    new_D5 = KQUAD(len=RING[27].len, k1=x[5])
    new_D6 = KQUAD(len=RING[31].len, k1=x[6])
    new_D7 = KQUAD(len=RING[5537].len, k1=x[7])
    changed_ele = [new_D1, new_D2, new_D3, new_D4, new_D5, new_D6, new_D7]
    refpts = [i for i in 1:length(RING)]
    twi = ADtwissring(RING, 0.0, 1, refpts, changed_idx, changed_ele)

    phase23 = twi[952].dmux - twi[916].dmux

    return phase23 - 2*pi
end

function multi_val_op(x0, niter, step, RING)
    target = 0.01
    x0_vals = zeros(7, niter)
    goal_vals = []
    grad_vals = zeros(7, niter)
    zero_idx = [9,13,17,23,27,31,5537]
    g0 = get_phase14_zero(x0, RING)
    for i in 1:niter
        grad9 = autodiff(Forward, get_phase14_zero, Duplicated, Duplicated(x0, [1.0,0.0,0.0,0.0,0.0,0.0,0.0]), Const(ESR_perturb))
        grad13 = autodiff(Forward, get_phase14_zero, Duplicated, Duplicated(x0, [0.0,1.0,0.0,0.0,0.0,0.0,0.0]), Const(ESR_perturb))
        grad17 = autodiff(Forward, get_phase14_zero, Duplicated, Duplicated(x0, [0.0,0.0,1.0,0.0,0.0,0.0,0.0]), Const(ESR_perturb))
        grad23 = autodiff(Forward, get_phase14_zero, Duplicated, Duplicated(x0, [0.0,0.0,0.0,1.0,0.0,0.0,0.0]), Const(ESR_perturb))
        grad27 = autodiff(Forward, get_phase14_zero, Duplicated, Duplicated(x0, [0.0,0.0,0.0,0.0,1.0,0.0,0.0]), Const(ESR_perturb))
        grad31 = autodiff(Forward, get_phase14_zero, Duplicated, Duplicated(x0, [0.0,0.0,0.0,0.0,0.0,1.0,0.0]), Const(ESR_perturb))
        grad5537 = autodiff(Forward, get_phase14_zero, Duplicated, Duplicated(x0, [0.0,0.0,0.0,0.0,0.0,0.0,1.0]), Const(ESR_perturb))

        x0[1] -= step * grad9[2]
        x0[2] -= step * grad13[2]
        x0[3] -= step * grad17[2]
        x0[4] -= step * grad23[2]
        x0[5] -= step * grad27[2]
        x0[6] -= step * grad31[2]
        x0[7] -= step * grad5537[2]

        new_phase = get_phase14_zero(x0, RING)
        println("init: ", g0, " now: ", new_phase, "at step ", i)
        x0_vals[:, i] = x0
        push!(goal_vals, new_phase)
        grad_vals[:, i] = [grad9[2], grad13[2], grad17[2], grad23[2], grad27[2], grad31[2], grad5537[2]]
        if abs(new_phase) < target 
            println("tuning finished at step ", i)
            break
        end
    end
    return x0_vals, goal_vals, grad_vals
end

x0_vals, goal_vals, grad_vals = multi_val_op(xinit, 10, 1e-5, ESR_crab)

plot_steps = 5
p1 = plot(1:plot_steps, x0_vals[1, 1:plot_steps], title = L"Evolution\ of\ k", xlabel = L"Iterations", ylabel = L"Strength (m^{-1})", label=L"k_1", line=:dash, marker=:circle)
for i in 2:7
    label_str = "k_{$i}"  
    full_label = latexstring(label_str)
    plot!(1:plot_steps, x0_vals[i, 1:plot_steps], label=full_label, line=:dash, marker=:circle)
end
p2 = plot(1:plot_steps, goal_vals[1:plot_steps], title = L"Evolution\ of\ \Delta \phi", xlabel = L"Iterations", ylabel = L"phase\ advance(rad)", legend = false, line=:dash, marker=:circle)
p3 = plot(1:plot_steps, grad_vals[1, 1:plot_steps], title = L"Evolution\ of\ gradient", xlabel = L"Iterations", 
    ylabel = L"\partial \frac{\Delta \phi}{\partial k}", label=L"k_1", line=:dash, marker=:circle)
for i in 2:7
    label_str = "k_{$i}"  
    full_label = latexstring(label_str)
    plot!(1:plot_steps, grad_vals[i, 1:plot_steps], label=full_label, line=:dash, marker=:circle)
end
plot(p1, p2, p3, layout = (3, 1), size=(800, 650))
