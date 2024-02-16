include("../JuTrack.jl")
using .JuTrack
using BenchmarkTools
using Distributed
include("drift_P.jl")

function matrix_to_array(matrix::Matrix{Float64})
    particles = zeros(Float64, size(matrix, 1)*size(matrix, 2))
    for i in 1:size(matrix, 1)
        for j in 1:size(matrix, 2)
            particles[(i-1)*size(matrix, 2)+j] = matrix[i, j]
        end
    end
    return particles
end

function array_to_matrix(array::Vector{Float64}, n::Int)
    particles = zeros(Float64, n, 6)
    for i in 1:n
        for j in 1:6
            particles[i, j] = array[(i-1)*6+j]
        end
    end
    return particles
end

np = 1000000
particles = zeros(Float64, np, 6)
particles[:, 2] .= 0.1
p_arr = matrix_to_array(particles)
D = DRIFT(len=3.0)
Q = KQUAD(k1=0.1, len=0.1)
line = [D]
beam = Beam(particles)
# # pass!(D, p_arr, np, beam, zeros(Float64, 6), zeros(Float64, 6, 6))
# # println(p_arr[1:6])
# using BenchmarkTools
# @btime begin
#     pass!(D, p_arr, np, beam, zeros(Float64, 6), zeros(Float64, 6, 6))
# end

function linepass_P!(line::Vector, particles::Beam)
    np = particles.nmacro
    particles6 = matrix_to_array(particles.r)
    if length(particles6) != np*6
        error("The number of particles does not match the length of the particle array")
    end
    noTarray = zeros(6)
    noRmatrix = [1.0 0.0 0.0 0.0 0.0 0.0; 
                 0.0 1.0 0.0 0.0 0.0 0.0; 
                 0.0 0.0 1.0 0.0 0.0 0.0; 
                 0.0 0.0 0.0 1.0 0.0 0.0; 
                 0.0 0.0 0.0 0.0 1.0 0.0; 
                 0.0 0.0 0.0 0.0 0.0 1.0]

    for i in eachindex(line)
        # ele = line[i]
        pass_P!(line[i], particles6, np, particles, noTarray, noRmatrix)        
    end
    rout = array_to_matrix(particles6, np)
    particles.r = rout
    return nothing
end

function linepass_distributed!(line::Vector, particles::Beam)
    np = particles.nmacro
    particles6 = matrix_to_array(particles.r)

    # Ensure correct size
    if length(particles6) != np*6
        error("The number of particles does not match the length of the particle array")
    end

    noTarray = zeros(6)
    noRmatrix = [1.0 0.0 0.0 0.0 0.0 0.0; 
                 0.0 1.0 0.0 0.0 0.0 0.0; 
                 0.0 0.0 1.0 0.0 0.0 0.0; 
                 0.0 0.0 0.0 1.0 0.0 0.0; 
                 0.0 0.0 0.0 0.0 1.0 0.0; 
                 0.0 0.0 0.0 0.0 0.0 1.0]

    # Wrap particle processing for pmap
    function process_wrapper(p_idx)
        particle_offset = (p_idx-1)*6 + 1
        particle_segment = @view particles6[particle_offset:particle_offset+5]
        # if !islost(particle_segment)
            for ele in line
                pass!(ele, particle_segment, 1, particles, noTarray, noRmatrix)
            end
        # end
        return particle_segment
    end

    # Use pmap for parallel execution across distributed workers
    updated_particles = pmap(process_wrapper, 1:np)

    # Reassemble the particles array
    for (i, particle_segment) in enumerate(updated_particles)
        particles6[((i-1)*6+1):(i*6)] = particle_segment
    end

    rout = array_to_matrix(particles6, np)
    particles.r = rout
end

@btime begin
    # @sync @distributed for i in 1:10000000
    #     # some complex computation
    #     i^2 + i - sin(i)
    # end
    linepass!(line, beam)
end
@btime begin
    # for i in 1:10000000
    #     # some complex computation
    #     i^2 + i - sin(i)
    # end
    linepass_P!(line, beam)
end

# @btime begin
#     # Threads.@threads for i in 1:10000000
#     #     # some complex computation
#     #     i^2 + i - sin(i)
#     # end
#     linepass_distributed!(line, beam)
# end