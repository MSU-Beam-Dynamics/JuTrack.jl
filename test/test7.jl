include("../src/JuTrack.jl")
using .JuTrack
include("ssrf_ring.jl")
using Plots  
using Enzyme
using ProgressMeter
using LinearAlgebra
using LaTeXStrings

RING = ssrf(-1.063770, 0)
xlist = randn(1000) * 0.003 
ylist = randn(1000) * 0.003 



function final_x(x)
    global RING
    beam = Beam([x[1] x[2] x[3] x[4] 0.0 0.0], energy=3.5e9)
    ringpass!(RING, beam, 1)
    return beam.r
end

grad = []
p = Progress(1000, 1)
for i in 1:1000
    next!(p)
    g = jacobian(Forward, final_x, [xlist[i], 0.0, ylist[i], 0.0], Val(4))
    push!(grad, g)
    sleep(0.1)
end

grad_m = zeros(1000, 24)
for i in 1:1000
    grad_m[i, :] = grad[i]
end

# 
function eigen_jacobian(J)
    # J is jacobian obtained from Enzyme. J is 24 * 1 Vector. 
    A = zeros(6, 4)

    # derivative wrt x
    A[1, 1] = J[1]
    A[2, 1] = J[2]
    A[3, 1] = J[3]
    A[4, 1] = J[4]
    A[5, 1] = J[5]
    A[6, 1] = J[6]
    # derivative wrt px
    A[1, 2] = J[7]
    A[2, 2] = J[8]
    A[3, 2] = J[9]
    A[4, 2] = J[10]
    A[5, 2] = J[11]
    A[6, 2] = J[12]
    # derivative wrt y
    A[1, 3] = J[13]
    A[2, 3] = J[14]
    A[3, 3] = J[15]
    A[4, 3] = J[16]
    A[5, 3] = J[17]
    A[6, 3] = J[18]
    # derivative wrt py
    A[1, 4] = J[19]
    A[2, 4] = J[20]
    A[3, 4] = J[21]
    A[4, 4] = J[22]
    A[5, 4] = J[23]
    A[6, 4] = J[24]
    eig_result = eigen(A[1:4, 1:4])
end

eigens = zeros(ComplexF64, 1000, 4)
for i in 1:1000
    eigens[i, :] = eigen_jacobian(grad_m[i, :]).values
end

vector_len = zeros(1000, 4)
for i in 1:1000
    e1 = eigen_jacobian(grad_m[i, :]).vectors
    vector_len[i, 1] = norm(e1[1, :])
    vector_len[i, 2] = norm(e1[2, :])
    vector_len[i, 3] = norm(e1[3, :])
    vector_len[i, 4] = norm(e1[4, :])
end

imag_list = zeros(1000, 4)
real_list = zeros(1000, 4)
for i in 1:1000
    imag_list[i, :] = imag(eigens[i, :])
    real_list[i, :] = real(eigens[i, :])
end
scatter(real_list[:,1], imag_list[:,1], label="1", xlabel="Real", ylabel="Imaginary")
scatter!(real_list[:,2], imag_list[:,2], label="2")
