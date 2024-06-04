include("../src/JuTrack.jl")
using .JuTrack
include("../src/demo/ssrf_ring.jl")
using Plots  
using Enzyme
using LinearAlgebra
using LaTeXStrings
using DelimitedFiles

RING = ssrf(-1.063770, 0)
# x: -0.04 to 0.04
angle_list = [ 0 + pi/36 * i for i in 0:36]
amp_list = [0.0 + 0.001 * i for i in 0:40]

N = length(angle_list) * length(amp_list)
particles = zeros(length(angle_list) * length(amp_list), 6)
for i in 1:length(angle_list)
    for j in 1:length(amp_list)
        particles[(i-1)*length(amp_list) + j, 1] = amp_list[j] * cos(angle_list[i])
        particles[(i-1)*length(amp_list) + j, 3] = amp_list[j] * sin(angle_list[i])
    end
end

beam = Beam(copy(particles), energy=3.5e9)
pringpass!(RING, beam, 1)
survived = findall(x -> x == 0, beam.lost_flag)

N_survived = length(survived)
survived_particles = zeros(N_survived, 6)
for i in 1:N_survived
    survived_particles[i, :] = particles[survived[i], :]
end

function TPSA_track_jacobian(amp, angle)
    x = amp * cos(angle)
    y = amp * sin(angle)
    X = CTPS(x, 1, 6, 1)
    PX = CTPS(0.0, 2, 6, 1)
    Y = CTPS(y, 3, 6, 1)
    PY = CTPS(0.0, 4, 6, 1)
    Z = CTPS(0.0, 5, 6, 1)
    DELTA = CTPS(0.0, 6, 6, 1)
    rin = [X, PX, Y, PY, Z, DELTA]
    ringpass_TPSA!(RING, rin, 1)
    jaco = zeros(6, 6)
    for i in 1:6
        jaco[i, :] = rin[i].map[2:7]
    end
    eigenvalues, eigenvectors = qr_eigen(jaco)
    return eigenvalues
end

eigens = zeros(N_survived, 6)
for i in 1:N_survived
    angle = atan(survived_particles[i, 3], survived_particles[i, 1])
    amplitude = sqrt(survived_particles[i, 1]^2 + survived_particles[i, 3]^2)
    eigenvalues = TPSA_track_jacobian(amplitude, angle)
    eigens[i, :] = eigenvalues
    println("Progress: ", i, "/", N_survived)
end
idx = findall(x -> x < 10, abs.(eigens[:,1]))
# disable legend
p1 = scatter(survived_particles[idx,1],survived_particles[idx,3],title=L"eigenvalue_1", colormap=:jet, markersize=3, legend=false, zcolor=eigens[idx,1], colorbar=true)
idx = findall(x -> x < 10, abs.(eigens[:,1]))
p2 = scatter(survived_particles[idx,1],survived_particles[idx,3], title=L"eigenvalue_2", colormap=:jet, markersize=3, legend=false ,zcolor=eigens[idx,2], colorbar=true)
idx = findall(x -> x < 10, abs.(eigens[:,1]))
p3 = scatter(survived_particles[idx,1],survived_particles[idx,3], title=L"eigenvalue_3", colormap=:jet, markersize=3, legend=false,zcolor=eigens[idx,3], colorbar=true)
idx = findall(x -> x < 10, abs.(eigens[:,1]))
p4 = scatter(survived_particles[idx,1],survived_particles[idx,3], title=L"eigenvalue_4", colormap=:jet, markersize=3, legend=false,zcolor=eigens[idx,4], colorbar=true)
p5 = scatter(survived_particles[:,1],survived_particles[:,3], title=L"eigenvalue_5", colormap=:jet, markersize=3, legend=false,zcolor=eigens[:,5], colorbar=true)
p6 = scatter(survived_particles[:,1],survived_particles[:,3], title=L"eigenvalue_6", colormap=:jet, markersize=3, legend=false,zcolor=eigens[:,6], colorbar=true)
plot(p1, p2, p3, p4, p5, p6, layout=(2,3), size=(1200,450))
# savefig("eigenvalues_abslessthan10.png")
# for i in 1046:N_survived
#     angle = atan(survived_particles[i, 3], survived_particles[i, 1])
#     amplitude = sqrt(survived_particles[i, 1]^2 + survived_particles[i, 3]^2)
#     gradx = autodiff(Forward, TPSA_track_jacobian, DuplicatedNoNeed, Duplicated(amplitude, 1.0), Const(angle))
#     open("outputx.csv", "a") do file
#         # Append a row of data
#         writedlm(file, reshape(gradx[1], 1, :), ',')
#     end
#     println("Progress: ", i, "/", N_survived)
# end

# gradx = readdlm("outputx.csv", ',', Float64)
# idx = findall(x -> x != 0, gradx[:,1])
# p1 = scatter(survived_particles[idx,1],survived_particles[idx,3],zcolor=log10.(abs.(gradx[idx,1])), title=L"de_1/dx", colormap=:jet)
# idx = findall(x -> x != 0, gradx[:,2])
# p2 = scatter(survived_particles[idx,1],survived_particles[idx,3],zcolor=log10.(abs.(gradx[idx,2])), title=L"de_2/dx", colormap=:jet)
# idx = findall(x -> x != 0, gradx[:,3])
# p3 = scatter(survived_particles[idx,1],survived_particles[idx,3],zcolor=log10.(abs.(gradx[idx,3])), title=L"de_3/dx", colormap=:jet)
# idx = findall(x -> x != 0, gradx[:,4])
# p4 = scatter(survived_particles[idx,1],survived_particles[idx,3],zcolor=log10.(abs.(gradx[idx,4])), title=L"de_4/dx", colormap=:jet)
# idx = findall(x -> x != 0, gradx[:,5])
# p5 = scatter(survived_particles[idx,1],survived_particles[idx,3],zcolor=log10.(abs.(gradx[idx,5])), title=L"de_5/dx", colormap=:jet)
# plot(p1, p2, p3, p4, p5, layout=(2,3), size=(1200,600))



# idx = findall(x -> x != 0, gradx[:,1].^2 .+ gradx[:,2].^2 .+ gradx[:,3].^2 .+ gradx[:,4].^2 .+ gradx[:,5].^2 + grady[:,1].^2 .+ grady[:,2].^2 .+ grady[:,3].^2 .+ grady[:,4].^2 .+ grady[:,5].^2)
# p5 = scatter(survived_particles[idx,1],survived_particles[idx,3],zcolor=log10.(gradx[idx,1].^2 .+ gradx[idx,2].^2 .+ gradx[idx,3].^2 
#     .+ gradx[idx,4].^2 .+ gradx[idx,5].^2 + grady[idx,1].^2 .+ grady[idx,2].^2 .+ grady[idx,3].^2 .+ grady[idx,4].^2 .+ gradx[idx,5].^2), title=L"sum(eigenvalues^2)", colormap=:jet)

# idx = findall(x -> x != 0, grady[:,1])
# p6 = scatter(survived_particles[idx,1],survived_particles[idx,3],zcolor=log10.(abs.(grady[idx,1])), title=L"de_1/dy", colormap=:jet)
# idx = findall(x -> x != 0, grady[:,2])
# p7 = scatter(survived_particles[idx,1],survived_particles[idx,3],zcolor=log10.(abs.(grady[idx,2])), title=L"de_2/dy", colormap=:jet)
# idx = findall(x -> x != 0, grady[:,3])
# p8 = scatter(survived_particles[idx,1],survived_particles[idx,3],zcolor=log10.(abs.(grady[idx,3])), title=L"de_3/dy", colormap=:jet)
# idx = findall(x -> x != 0, grady[:,4])
# p9 = scatter(survived_particles[idx,1],survived_particles[idx,3],zcolor=log10.(abs.(grady[idx,4])), title=L"de_4/dy", colormap=:jet)
# plot(p6, p7, p8, p9, layout=(2,2), size=(800,600))

# eigens = zeros(N_survived, 6)
# for i in 1:N_survived
#     eigens[i, :] = TPSA_track_jacobian(survived_particles[i, 1], survived_particles[i, 3])
#     println("Progress: ", i, "/", N_survived)
# end

# scatter(survived_particles[:,1],survived_particles[:,3],zcolor=log10.(eigens[:,1].^2 .+ eigens[:,2].^2 .+ eigens[:,3].^2 .+ eigens[:,4].^2), title=L"sum(eigenvalues^2)", colormap=:jet)