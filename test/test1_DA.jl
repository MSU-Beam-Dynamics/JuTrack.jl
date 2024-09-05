using JuTrack
using Serialization
using Plots
using Printf

RING =spear()
function sigmoid(x, k=1.0, limit=1.0)
    return 1.0 / (1.0 + exp(-k*(abs(x)-limit)))
end
function htan(x, k=1.0, limit=1.0)
    return (exp(k*(abs(x)-limit)) - exp(-k*(abs(x)-limit))) / (exp(k*(abs(x)-limit)) + exp(-k*(abs(x)-limit)))+1.0
end
function softplus(x)
    return log(1.0 + exp(abs(x)))
end

function f(x1, x2)
    # RING = [DRIFT(len=0.1), KSEXT(name="SDM", len=0.21, k2=0.0), DRIFT(len=0.1), KSEXT(name="SFM",len=0.21, k2=0.0), DRIFT(len=0.1)]
    changed_id1 = [19,33,384,398,433,447,793,807]
    changed_id2 = [23,29,388,394,437,443,797,803]
    changed_ids = vcat(changed_id1, changed_id2)
    changed_elems = [changed_ids[i] in changed_id1 ? KSEXT(len=0.21, k2=x1) : KSEXT(len=0.21, k2=x2) for i in 1:length(changed_ids)]
    # changed_ids = changed_id1
    # changed_elems = [KSEXT(len=0.21, k2=x1) for i in 1:length(changed_ids)]
    nturns = 1000
    amp_nstep = 50
    amp_step = 0.001
    angle_steps = 11
    E = 3.0e9

    # println("number of threads used for parallel comptuing: ", Threads.nthreads())
    angle_list = [pi/(angle_steps-1) * i for i in 0:angle_steps-1]
    amp_list = [amp_step * i for i in 1:amp_nstep]

    particles = zeros(length(angle_list) * length(amp_list), 6)
    for i in 1:length(angle_list)
        for j in 1:length(amp_list)
            particles[(i-1)*length(amp_list) + j, 1] = amp_list[j] * cos(angle_list[i])
            particles[(i-1)*length(amp_list) + j, 3] = amp_list[j] * sin(angle_list[i])
        end
    end

    beam = Beam(particles, energy=E)
    ADringpass!(RING, beam, nturns, changed_ids, changed_elems)
    loss = 0.0
    for i in 1:length(angle_list) * length(amp_list)
        loss += htan(beam.r[i, 1], 10.0, 1.0)
    end
    return loss
end


@time g1 = autodiff(Forward, f, Duplicated, Duplicated(-10.0, 1.0), Const(10.0))
@time g2 = autodiff(Forward, f, Duplicated, Const(-100.0), Duplicated(50.0, 1.0))


xlist = -1:0.01:1
ylist = sigmoid.(xlist, 15.0)
plot(xlist, ylist, xlabel="x", ylabel="y", title="Sigmoid function")


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

        lr1 = lr *abs( g1[1] / grad_x1)
        lr2 = lr *abs( g2[1] / grad_x2)
        # Update variables
        x1 -= lr1 * grad_x1
        x2 -= lr2 * grad_x2
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
x1_init = -30.0
x2_init = 30.0

# Run gradient descent
x1_his, x2_his, g1_his, g2_his, f_his = gradient_descent(x1_init, x2_init, lr=0.01, max_iter=1000)

# Plot the results
p1 = plot(0:length(x1_his)-1, x1_his, label="x1", xlabel="Iteration", ylabel="Value", title="Gradient Descent")
plot!(0:length(x2_his)-1, x2_his, label="x2")
p2 = plot(g1_his, label="grad1", xlabel="Iteration", ylabel="Gradient")
plot!(g2_his, label="grad2")
p3 = plot(f_his, label="f", xlabel="Iteration", ylabel="Objective")
plot(p1, p2, p3, layout=(3,1), legend=:topleft)

amp_nstep = 50
amp_step = 0.0001
angle_steps = 20
E = 3.0e9

# println("number of threads used for parallel comptuing: ", Threads.nthreads())
angle_list = [pi/(angle_steps-1) * i for i in 0:angle_steps-1]
amp_list = [amp_step * i for i in 1:amp_nstep]

particles = zeros(length(angle_list) * length(amp_list), 6)
for i in 1:length(angle_list)
    for j in 1:length(amp_list)
        particles[(i-1)*length(amp_list) + j, 1] = amp_list[j] * cos(angle_list[i])
        particles[(i-1)*length(amp_list) + j, 3] = amp_list[j] * sin(angle_list[i])
    end
end
beam = Beam(particles, energy=E)
beam1 = Beam(particles, energy=E)
@btime ringpass!(RING, beam, 100)
@btime ringpass!(RING1, beam1, 100)