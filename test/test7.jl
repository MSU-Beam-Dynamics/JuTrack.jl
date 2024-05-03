include("../src/JuTrack.jl")
using .JuTrack
include("ssrf_ring.jl")
using Plots  
using Enzyme
using ProgressMeter
using LinearAlgebra
using LaTeXStrings

RING = ssrf(-1.063770, 0)
# x: -0.04 to 0.04
xlist = [-0.04 + 0.001 * i for i in 0:80]
ylist = [0.0 + 0.001 * i for i in 0:40]
N = length(xlist) * length(ylist)

particles = zeros(length(xlist) * length(ylist), 6)
for i in 1:length(xlist)
    for j in 1:length(ylist)
        particles[(i-1)*length(ylist) + j, 1] = xlist[i]
        particles[(i-1)*length(ylist) + j, 3] = ylist[j]
    end
end

beam = Beam(copy(particles), energy=3.5e9)
pringpass!(RING, beam, 1000)

# function final_x(x)
#     global RING
#     beam = Beam([x[1] x[2] x[3] x[4] 0.0 0.0], energy=3.5e9)
#     ringpass!(RING, beam, 1000)
#     return beam.r
# end

# grad = []
# p = Progress(N, 1)
# for i in 1:N
#     next!(p)
#     g = jacobian(Forward, final_x, [particles[i, 1], 0.0, particles[i, 3], 0.0], Val(4))
#     push!(grad, g)
#     sleep(0.1)
# end

# grad_m = zeros(N, 24)
# for i in 1:N
#     grad_m[i, :] = grad[i]
# end

# # 
# function eigen_jacobian(J)
#     # J is jacobian obtained from Enzyme. J is 24 * 1 Vector. 
#     A = zeros(6, 4)

#     # derivative wrt x
#     A[1, 1] = J[1]
#     A[2, 1] = J[2]
#     A[3, 1] = J[3]
#     A[4, 1] = J[4]
#     A[5, 1] = J[5]
#     A[6, 1] = J[6]
#     # derivative wrt px
#     A[1, 2] = J[7]
#     A[2, 2] = J[8]
#     A[3, 2] = J[9]
#     A[4, 2] = J[10]
#     A[5, 2] = J[11]
#     A[6, 2] = J[12]
#     # derivative wrt y
#     A[1, 3] = J[13]
#     A[2, 3] = J[14]
#     A[3, 3] = J[15]
#     A[4, 3] = J[16]
#     A[5, 3] = J[17]
#     A[6, 3] = J[18]
#     # derivative wrt py
#     A[1, 4] = J[19]
#     A[2, 4] = J[20]
#     A[3, 4] = J[21]
#     A[4, 4] = J[22]
#     A[5, 4] = J[23]
#     A[6, 4] = J[24]
#     eig_result = eigen(A[1:4, 1:4])
# end

# eigens = zeros(ComplexF64, N, 4)
# for i in 1:N
#     eigens[i, :] = eigen_jacobian(grad_m[i, :]).values
# end

# vector_len = zeros(N, 4)
# for i in 1:N
#     e1 = eigen_jacobian(grad_m[i, :]).vectors
#     vector_len[i, 1] = norm(e1[1, :])
#     vector_len[i, 2] = norm(e1[2, :])
#     vector_len[i, 3] = norm(e1[3, :])
#     vector_len[i, 4] = norm(e1[4, :])
# end

# imag_list = zeros(N, 4)
# real_list = zeros(N, 4)
# for i in 1:N
#     imag_list[i, :] = imag(eigens[i, :])
#     real_list[i, :] = real(eigens[i, :])
# end

# beam = Beam(particles, energy=3.5e9)
# pringpass!(RING, beam, 1000)

# idx_num = findall(x -> !isnan(x), real_list[:, 1])
# scatter(real_list[:,1], imag_list[:,1], label="1", xlabel="Real", ylabel="Imaginary")
# # scatter!(real_list[:,2], imag_list[:,2], label="2")

# norm_real = sqrt.(real_list[:, 1].^2 + real_list[:, 2].^2 + real_list[:, 3].^2 + real_list[:, 4].^2)
# norm_imag = sqrt.(imag_list[:, 1].^2 + imag_list[:, 2].^2 + imag_list[:, 3].^2 + imag_list[:, 4].^2)
# idx1 = findall(x -> x < 100, log10.(norm_real))
# idx2 = findall(x -> x < 100, log10.(norm_imag))
# p1 = scatter(particles[idx1, 1], particles[idx1, 3], zcolor=log10.(norm_real[idx1]), label="log(Norm of Real Eigenvalue)", markersize=5, colormap=:jet)
# p2 = scatter(particles[idx2, 1], particles[idx2, 3], zcolor=log10.(norm_imag[idx2].+1e-5), label="log(Norm of imag Eigenvalue + 1e-5)", markersize=5, colormap=:jet)
# plot(p1, p2, layout=(1, 2), size=(1200, 400))