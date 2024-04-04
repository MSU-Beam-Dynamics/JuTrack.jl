using Random
using Plots
using Enzyme

function kernel(x)
    return 1 / sqrt(2 * pi) * exp(-x^2 / 2)
end

# Define the KDE function
function kde(x, x_drawn, h)
    pdf = 0.0
    for xi in x_drawn
        pdf += kernel((x - xi) / h)
    end
    return pdf / (length(x_drawn) * h)
end

gamma = 1.0  # Damping coefficient
D = 2.0  # Diffusion coefficient
Δt = 0.01  # Time step
N = 1000  # Number of steps
n = 100000  # Number of particles

# Simulate the motion
function simulate(v0, gamma, D, Δt, N)
    v = zeros(N)
    v[1] = v0
    for n in 1:(N-1)
        ξ = randn() 
        v[n+1] = v[n] - gamma * v[n] * Δt + sqrt(2 * D * Δt) * ξ
    end
    return v
end

v_init = randn(n)  # Initial velocity
v_final = zeros(n)

for i in 1:n
    v_final[i] = simulate(v_init[i], gamma, D, Δt, N)[end]
end

x_grad = zeros(1000, 2)
x = range(min(v_init...), stop=max(v_init...), length=1000)

for i in eachindex(x)
    grad = autodiff(Forward, kde, Duplicated, Duplicated(x[i], 1.0), Const(v_final), Const(0.1))
    x_grad[i, 1] = grad[1]
    x_grad[i, 2] = grad[2]
end
p1 = plot(x, x_grad[:, 1], label="KDE density", xlabel="x",ylabel="density")
plot!(x, x_grad[:, 2], label="gradient", xlabel="x", xlims=(-4.0, 4.0))
# histogram normalized to 1
histogram!(v_final, bins=100, alpha=0.3,xlabel="x",xlims=(-4.0, 4.0), normed=true,label=false)