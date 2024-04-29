include("../src/JuTrack.jl")
using .JuTrack
using Enzyme
using LinearAlgebra
using Plots

function HenonMap(x, p)
    v = 0.195
    lamda = 1.0
    x_new = cos(2*pi*v) * x + sin(2*pi*v) * (p - lamda*x^2)
    p_new = -sin(2*pi*v) * x + cos(2*pi*v) * (p - lamda*x^2)
    return x_new, p_new
end

xlist = -0.6:0.01:0.6
plist = -0.6:0.01:0.6

data = zeros(length(xlist) * length(plist), 2)
for i in 1:length(xlist)
    for j in 1:length(xlist)
        data[(i-1)*length(xlist) + j, 1] = xlist[i]
        data[(i-1)*length(xlist) + j, 2] = plist[j]
    end
end

function henon_nturn(xin)
    x = xin[1]
    p = xin[2]
    for j in 1:100
        x, p = HenonMap(x, p)
    end

    return x, p
end

grad_n = []
# complex E
E = zeros(ComplexF64, size(data, 1), 2)
idx = []
for i in 1:size(data, 1)
    global idx
    grad = jacobian(Forward, henon_nturn, data[i, :], Val(2))
    if isnan(grad[1][1]) || isnan(grad[1][2]) || isnan(grad[2][1]) || isnan(grad[2][2]) || isinf(grad[1][1]) || isinf(grad[1][2]) || isinf(grad[2][1]) || isinf(grad[2][2]) || henon_nturn(data[i, :])[1] > 2.0 || henon_nturn(data[i, :])[2] > 2.0
        e = [0.0, 0.0]
    else
        e = eigen([grad[1][1] grad[1][2]; grad[2][1] grad[2][2]]).values
        idx = push!(idx, i)
    end
    E[i, :] = e
    push!(grad_n, grad)
end

imag_list = zeros(size(data, 1), 2)
real_list = zeros(size(data, 1), 2)
for i in 1:size(data, 1)
    imag_list[i, :] = imag(E[i, :])
    real_list[i, :] = real(E[i, :])
end

log_real_list1 = log10.(abs.(real_list[idx, 1]).+1e-5)
log_real_list2 = log10.(abs.(real_list[idx, 2]).+1e-5)
# find the index of log_real_list1 less than 5
idx1 = findall(x -> x < 5, log_real_list1)
idx2 = findall(x -> x < 5, log_real_list2)
p1 = scatter(data[idx[idx1], 1], data[idx[idx1], 2], zcolor=log_real_list1[idx1], label="log(|Real Eigenvalue1|+1e-5)", markersize=2, colormap=:jet)
p2 = scatter(data[idx[idx2], 1], data[idx[idx2], 2], zcolor=log_real_list2[idx2], label="log(|Real Eigenvalue2|+1e-5)", markersize=2, colormap=:jet)
p3 = scatter(data[idx, 1], data[idx, 2], zcolor=imag_list[idx, 1], label="Imag Eigenvalue1", markersize=2, colormap=:jet)
p4 = scatter(data[idx, 1], data[idx, 2], zcolor=imag_list[idx, 2], label="Imag Eigenvalue2", markersize=2, colormap=:jet)
plot(p1, p2, p3, p4, layout=(2, 2), size=(1200, 800))

savefig("eigen_100turns.png")