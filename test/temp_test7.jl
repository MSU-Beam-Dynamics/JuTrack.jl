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
pringpass!(RING, beam, 100)
survived = findall(x -> x == 0, beam.lost_flag)

function final_x(x)
    global RING
    beam = Beam([x[1] x[2] x[3] x[4] 0.0 0.0], energy=3.5e9)
    ringpass!(RING, beam, 1)
    return beam.r[1, 1], beam.r[1, 2], beam.r[1, 3], beam.r[1, 4]
end


E = zeros(ComplexF64, N, 2)
Vectors = zeros(ComplexF64, N, 2, 2)
grad = []
p = Progress(N, 1)
for i in 1:N
    next!(p)
    g = jacobian(Forward, final_x, [particles[i, 1], 0.0, particles[i, 3], 0.0], Val(4))
    # e = eigen([g[1][1] g[1][2] g[1][3] g[1][4]; g[2][1] g[2][2] g[2][3] g[2][4];
    #             g[3][1] g[3][2] g[3][3] g[3][4]; g[4][1] g[4][2] g[4][3] g[4][4]]).values
    # vectors = eigen([g[1][1] g[1][2] g[1][3] g[1][4]; g[2][1] g[2][2] g[2][3] g[2][4];
    #             g[3][1] g[3][2] g[3][3] g[3][4]; g[4][1] g[4][2] g[4][3] g[4][4]]).vectors
    e = eigen([g[1][1] g[1][2]; g[2][1] g[2][2]]).values
    vectors = eigen([g[1][1] g[1][2]; g[2][1] g[2][2]]).vectors
    Vectors[i, :, :] = vectors
    push!(grad, g)
    E[i, :] = e
    sleep(0.01)
end


imag_list = zeros(N, 2)
real_list = zeros(N, 2)
for i in 1:N
    imag_list[i, 1] = imag(E[i, 1])
    imag_list[i, 2] = imag(E[i, 2])
    # imag_list[i, 3] = imag(E[i, 3])
    # imag_list[i, 4] = imag(E[i, 4])
    real_list[i, 1] = real(E[i, 1])
    real_list[i, 2] = real(E[i, 2])
    # real_list[i, 3] = real(E[i, 3])
    # real_list[i, 4] = real(E[i, 4])
end

# calculate the length of eigenvectors
function complex_vector_norm(v)
    return sqrt(sum(abs(c)^2 for c in v))
end

idx1 = findall(x -> x < 10, abs.(real_list[survived, 1]))
idx2 = findall(x -> x < 10, abs.(real_list[survived, 2]))
idx3 = findall(x -> x < 10, abs.(imag_list[survived, 1]))
idx4 = findall(x -> x < 10, abs.(imag_list[survived, 2]))
# finx the index of both real1 and real12 less than 0

p1 = scatter(particles[survived, 1], particles[survived, 3], zcolor=log10.(abs.(real_list[survived, 1])),
    label="Re 1", markersize=3, colormap=:jet, markershape=:rect, xlabel="x(m)", 
    ylabel="y(m)", ylimits=(0.0, 0.04), legend=false, title="log(|Re 1|)", colorbar=true)
p2 = scatter(particles[survived, 1], particles[survived, 3], zcolor=log10.(abs.(real_list[survived, 2])),
    label="Re 2", markersize=3, colormap=:jet, markershape=:rect, xlabel="x(m)", 
    ylabel="y(m)", ylimits=(0.0, 0.04), legend=false, title="log(|Re 2|)", colorbar=true)
p3 = scatter(particles[survived, 1], particles[survived, 3], zcolor=log10.(abs.(imag_list[survived, 1]).+1e-5),
    label="Im 1", markersize=3, colormap=:jet, markershape=:rect, xlabel="x(m)", 
    ylabel="y(m)", ylimits=(0.0, 0.04), legend=false, title="log(|Im 2|+1e-5)", colorbar=true)
p4 = scatter(particles[survived, 1], particles[survived, 3], zcolor=log10.(abs.(imag_list[survived, 2]).+1e-5),
    label="Im 2", markersize=3, colormap=:jet, markershape=:rect, xlabel="x(m)", 
    ylabel="y(m)", ylimits=(0.0, 0.04), legend=false, title="log(|Re 2|+1e-5)", colorbar=true)
plot(p1, p2, p3, p4, layout=(2, 2), size=(900, 600))

# norm_vecs = zeros(N, 4)
# for i in 1:N
#     for j in 1:4
#         norm_vecs[i, j] = complex_vector_norm(Vectors[i, :, j])
#     end
# end

# max_over_min = zeros(N, 1)
# for i in 1:N
#     max_over_min[i] = maximum(norm_vecs[i, :]) / minimum(norm_vecs[i, :])
# end
# scatter(particles[survived, 1], particles[survived, 3], zcolor=max_over_min[survived], 
#     label="log(Norm of imag Eigenvalue + 1e-5)", markersize=3, colormap=:jet, markershape=:rect, xlabel="x(m)", 
#     ylabel="y(m)", ylimits=(0.0, 0.04), legend=false, colorbar=true, title="Max/Min of Eigenvectors")

# max_over_min = zeros(N, 4)
# for i in 1:N
#     max_over_min[i, 1] = maximum(abs.(Vectors[i, :, 1])) / minimum(abs.(Vectors[i, :, 1]))
#     max_over_min[i, 2] = maximum(abs.(Vectors[i, :, 2])) / minimum(abs.(Vectors[i, :, 2]))
#     max_over_min[i, 3] = maximum(abs.(Vectors[i, :, 3])) / minimum(abs.(Vectors[i, :, 3]))
#     max_over_min[i, 4] = maximum(abs.(Vectors[i, :, 4])) / minimum(abs.(Vectors[i, :, 4]))
#     if max_over_min[i, 1] == Inf
#         max_over_min[i, 1] = 1e10
#     end
#     if max_over_min[i, 2] == Inf
#         max_over_min[i, 2] = 1e10
#     end
#     if max_over_min[i, 3] == Inf
#         max_over_min[i, 3] = 1e10
#     end
#     if max_over_min[i, 4] == Inf
#         max_over_min[i, 4] = 1e10
#     end
# end
# idx1 = findall(x -> x < 1e10, max_over_min[survived, 1])
# p1 = scatter(particles[survived[idx1], 1], particles[survived[idx1], 3], zcolor=log10.(max_over_min[survived[idx1], 1]), 
#     label="log(Norm of imag Eigenvalue + 1e-5)", markersize=3, colormap=:jet, markershape=:rect, xlabel="x(m)", 
#     ylabel="y(m)", ylimits=(0.0, 0.04), legend=false, colorbar=true, title="Max/Min of Eigenvector1")
# p2 = scatter(particles[survived, 1], particles[survived, 3], zcolor=log10.(max_over_min[survived, 2]),
#     label="log(Norm of imag Eigenvalue + 1e-5)", markersize=3, colormap=:jet, markershape=:rect, xlabel="x(m)", 
#     ylabel="y(m)", ylimits=(0.0, 0.04), legend=false, colorbar=true, title="Max/Min of Eigenvector2")
# p3 = scatter(particles[survived, 1], particles[survived, 3], zcolor=log10.(max_over_min[survived, 3]),
#     label="log(Norm of imag Eigenvalue + 1e-5)", markersize=3, colormap=:jet, markershape=:rect, xlabel="x(m)", 
#     ylabel="y(m)", ylimits=(0.0, 0.04), legend=false, colorbar=true, title="Max/Min of Eigenvector3")
# p4 = scatter(particles[survived, 1], particles[survived, 3], zcolor=log10.(max_over_min[survived, 4]),
#     label="log(Norm of imag Eigenvalue + 1e-5)", markersize=3, colormap=:jet, markershape=:rect, xlabel="x(m)", 
#     ylabel="y(m)", ylimits=(0.0, 0.04), legend=false, colorbar=true, title="Max/Min of Eigenvector4")
# plot(p1, p2, p3, p4, layout=(2, 2), size=(900, 900))

# norm_real = sqrt.(real_list[survived, 1].^2 + real_list[survived, 2].^2 + real_list[survived, 3].^2 + real_list[survived, 4].^2)
# norm_imag = sqrt.(imag_list[survived, 1].^2 + imag_list[survived, 2].^2 + imag_list[survived, 3].^2 + imag_list[survived, 4].^2)
# idx1 = findall(x -> x < 5, log10.(norm_real))
# idx2 = findall(x -> x < 5, log10.(norm_imag))
# p1 = scatter(particles[survived[idx1], 1], particles[survived[idx1], 3], zcolor=log10.(norm_real[idx1]), 
#     label="log(Norm of Real Eigenvalue)", markersize=3, colormap=:jet, markershape=:rect, xlabel="x(m)", 
#     ylabel="y(m)", ylimits=(0.0, 0.04), legend=false)
# p2 = scatter(particles[survived[idx2], 1], particles[survived[idx2], 3], zcolor=log10.(norm_imag[idx2].+1e-5), 
#     label="log(Norm of imag Eigenvalue + 1e-5)", markersize=3, colormap=:jet, markershape=:rect, xlabel="x(m)", 
#     ylabel="y(m)", ylimits=(0.0, 0.04), legend=false)
# plot(p1, p2, layout=(1, 2), size=(900, 330))

# data_all = hcat(particles[:,1], particles[:,3], real_list, imag_list)
# data_survived = data_all[survived, :]


# DA calculation
# na = 37
# nl = 80
# angle_list = range(0, stop=pi, length=na)
# llist = [(0.0005*i) for i in 1:nl]
# particle = zeros(na*nl, 6)
# DA = zeros(na, 2)
# p = Progress(na, 1)
# for i in 1:na
#     next!(p)
#     for j in 1:nl
#         particle = [llist[j]*cos(angle_list[i]) 0.0 llist[j]*sin(angle_list[i]) 0.0 0.0 0.0]
#         beam = Beam(particle, energy=3.5e9)
#         ringpass!(RING, beam, 1000)
#         lost_flag = beam.lost_flag
#         if lost_flag[1] == 1
#             DA[i, 1] = llist[j] * cos(angle_list[i])
#             DA[i, 2] = llist[j] * sin(angle_list[i])
#             break
#         end
#     end
# end
