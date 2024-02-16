include("../lattice/canonical_elements.jl")
using BenchmarkTools
include("../lattice/beam.jl")
using Distributed

# Add workers
addprocs(14) 

@everywhere begin
    function ATmultmv!(r::AbstractVector{Float64}, A::Matrix{Float64})
        temp = zeros(6)
        for i in 1:6
            for j in 1:6
                temp[i] += A[i, j] * r[j]
            end
        end
        for i in 1:6
            r[i] = temp[i]
        end
        return nothing
    end

    function ATaddvv!(r::AbstractVector{Float64}, dr::Array{Float64,1})
        for i in 1:6
            r[i] += dr[i]
        end
        return nothing
    end
    
    function fastdrift!(r::AbstractVector{Float64}, NormL::Float64, le::Float64)
        # NormL is computed externally to speed up calculations
        # in the loop if momentum deviation (delta) does not change
        # such as in 4-th order symplectic integrator w/o radiation
        # AT uses small angle approximation pz = 1 + delta. 
        # Here we use pz = sqrt((1 + delta)^2 - px^2 - py^2) for precise calculation
        r[1] += NormL * r[2]
        r[3] += NormL * r[4]
        # r[6] += NormL * (r[2]^2 + r[4]^2) / (2*(1+r[5]))
        r[5] += NormL * (1.0 + r[6]) - le
        return nothing
    end
    
    function drift6!(r::AbstractVector{Float64}, le::Float64)
        # AT uses small angle approximation pz = 1 + delta. 
        # Here we use pz = sqrt((1 + delta)^2 - px^2 - py^2) for precise calculation
        NormL = le / sqrt(((1.0 + r[6])^2 - r[2]^2 - r[4]^2))
        r[1] += NormL * r[2]
        r[3] += NormL * r[4]
        # r[6] += NormL * (r[2]^2 + r[4]^2) / (2*(1+r[5])) # for linearized approximation
        r[5] += NormL * (1.0 + r[6]) - le
        return nothing
    end 
end

@everywhere function process_particle(r6::AbstractVector{Float64}, le::Float64, T1::Array{Float64,1}, T2::Array{Float64,1}, 
    R1::Array{Float64,2}, R2::Array{Float64,2}, noTarray::Array{Float64,1}, noRmatrix::Array{Float64,2})
    if T1 != noTarray
        ATaddvv!(r6, T1)
    end
    if R1 != noRmatrix
        ATmultmv!(r6, R1)
    end
    drift6!(r6, le)
    if R2 != noRmatrix
        ATmultmv!(r6, R2)
    end
    if T2 != noTarray
        ATaddvv!(r6, T2)
    end
    return r6
end

function DriftPass_parallel!(r_in::Array{Float64,1}, le::Float64, T1::Array{Float64,1}, T2::Array{Float64,1}, 
    R1::Array{Float64,2}, R2::Array{Float64,2}, RApertures::Array{Float64,1}, EApertures::Array{Float64,1}, 
    num_particles::Int, lost_flags::Array{Int64,1}, noTarray::Array{Float64,1}, noRmatrix::Array{Float64,2})
    
    futures = []
    for c in 1:num_particles
        if lost_flags[c] == 1
            continue
        end
        r6 = @view r_in[(c-1)*6+1:c*6]
        if !isnan(r6[1])
            push!(futures, @spawn process_particle(copy(r6), le, T1, T2, R1, R2, noTarray, noRmatrix))
        end
    end
    
    # Fetch results and update r_in
    for (index, future) in enumerate(futures)
        r6_updated = fetch(future)
        r_in[(index-1)*6+1:index*6] = r6_updated
    end

    return nothing
end
function pass_P!(ele::DRIFT, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam, noTarray::Array{Float64,1}, noRmatrix::Array{Float64,2})
    # ele: EDRIFT
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    lost_flags = particles.lost_flag
    DriftPass_parallel!(r_in, ele.len, ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, num_particles, lost_flags, noTarray, noRmatrix)
    return nothing
end

function matrix_to_array(matrix::Matrix{Float64})
    particles = zeros(Float64, size(matrix, 1)*size(matrix, 2))
    for i in 1:size(matrix, 1)
        for j in 1:size(matrix, 2)
            particles[(i-1)*size(matrix, 2)+j] = matrix[i, j]
        end
    end
    return particles
end
r = zeros(Float64, 100000, 6)
r[:, 2] .= 0.1
r_arr = matrix_to_array(r)
beam = Beam(r)
D = DRIFT(len=3.0)
pass_P!(D, r_arr, beam.nmacro, beam, zeros(Float64, 6), zeros(Float64, 6, 6))
println(r_arr[1:6])
println(Threads.nthreads())
@btime begin    
    pass_P!(D, r_arr, beam.nmacro, beam, zeros(Float64, 6), zeros(Float64, 6, 6))
end
# using Enzyme
# Enzyme.API.runtimeActivity!(true) 
# function fE(len)
#     # r = zeros(Float64, 1000000, 6)
#     # r[:, 2] .= 0.1
#     # r_arr = matrix_to_array(r)
#     # beam = Beam(r)
#     # D = DRIFT(len=len)
#     # pass_P!(D, r_arr, beam.nmacro, beam, zeros(Float64, 6), zeros(Float64, 6, 6))
#     s = 1.0
#     Threads.@threads for i in 1:10000
#         s += len[1]
#     end
#     return s
# end

# using Dagger

# function fE_dagger(len)
#     # Create an array to hold spawned tasks
#     tasks = [Dagger.@spawn len[1] for i in 1:10000]

#     # Collect the results of the tasks and sum them up
#     s = 1.0 + sum(collect(tasks))

#     return s
# end
# @btime begin
#     fE(3.0)
# end
# @btime begin
#     fE_dagger(3.0)
# end
# println(fE_dagger(3.0))
# grad = gradient(Forward, fE, [3.0])
# println(grad)