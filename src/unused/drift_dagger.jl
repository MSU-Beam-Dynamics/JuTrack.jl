using Distributed
addprocs(14; exeflags="--threads=2")
@everywhere using Dagger
@everywhere include("../lattice/canonical_elements.jl")
@everywhere using BenchmarkTools
@everywhere include("../lattice/beam.jl")



@everywhere begin
    const CoordLimit = 1.0
    const AngleLimit = 1.0
end

@everywhere begin
function ATmultmv!(r::AbstractVector{Float64}, A::Matrix{Float64})
    # multiplies 6-component column vector r by 6x6 matrix R: as in A*r
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
function DriftPass!(r_in::Array{Float64,1}, le::Float64, T1::Array{Float64,1}, T2::Array{Float64,1}, 
    R1::Array{Float64,2}, R2::Array{Float64, 2}, RApertures::Array{Float64,1}, EApertures::Array{Float64,1}, 
    num_particles::Int, lost_flags::Array{Int64,1}, noTarray::Array{Float64,1}, noRmatrix::Array{Float64,2})
    Threads.@threads for c in 1:num_particles
    # for c in 1:num_particles
        if lost_flags[c] == 1
            continue
        end
        r6 = @view r_in[(c-1)*6+1:c*6]
        if !isnan(r6[1])
            # Misalignment at entrance
            if T1 != noTarray
                ATaddvv!(r6, T1)
            end
            if R1 != noRmatrix
                ATmultmv!(r6, R1)
            end
            # Check physical apertures at the entrance of the magnet
            # if RApertures !== nothing
            #     checkiflostRectangularAp!(r6, RApertures)
            # end
            # if EApertures !== nothing
            #     checkiflostEllipticalAp!(r6, EApertures)
            # end
            drift6!(r6, le)
            # Check physical apertures at the exit of the magnet
            # if RApertures !== nothing
            #     checkiflostRectangularAp!(r6, RApertures)
            # end
            # if EApertures !== nothing
            #     checkiflostEllipticalAp!(r6, EApertures)
            # end
            # Misalignment at exit
            if R2 != noRmatrix
                ATmultmv!(r6, R2)
            end
            if T2 != noTarray
                ATaddvv!(r6, T2)
            end
            if r6[1] > CoordLimit || r6[2] > AngleLimit || r6[1] < -CoordLimit || r6[2] < -AngleLimit
                lost_flags[c] = 1
            end
        end
    end
    return nothing
end

function process_particle!(c, r_in, le, T1, T2, R1, R2, lost_flags, noTarray, noRmatrix)
    if lost_flags[c] != 1
        r6 = @view r_in[(c-1)*6+1:c*6]
        if !isnan(r6[1])
            # Misalignment at entrance
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
            # Assuming CoordLimit and AngleLimit are defined or passed to the function
            if r6[1] > CoordLimit || r6[2] > AngleLimit || r6[1] < -CoordLimit || r6[2] < -AngleLimit
                lost_flags[c] = 1
            end
        end
    end
end

end


function DriftPassDagger!(r_in, le, T1, T2, R1, R2, lost_flags, noTarray, noRmatrix, num_particles)
    tasks = [Dagger.@spawn process_particle!(c, r_in, le, T1, T2, R1, R2, lost_flags, noTarray, noRmatrix) for c in 1:num_particles]
    
    # Wait for all tasks to complete
    collect(tasks)
    
    return nothing
end

function pass_P!(ele::DRIFT, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam, noTarray::Array{Float64,1}, noRmatrix::Array{Float64,2})
    # ele: EDRIFT
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    lost_flags = particles.lost_flag
    # DriftPass!(r_in, ele.len, ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, num_particles, lost_flags, noTarray, noRmatrix)
    DriftPassDagger!(r_in, ele.len, ele.T1, ele.T2, ele.R1, ele.R2, lost_flags, noTarray, noRmatrix, num_particles)
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

# r = zeros(Float64, 100000, 6)
# r[:, 2] .= 0.1
# r_arr = matrix_to_array(r)
# beam = Beam(r)
# D = DRIFT(len=3.0)
# pass_P!(D, r_arr, beam.nmacro, beam, zeros(Float64, 6), zeros(Float64, 6, 6))
# println(r_arr[1:6])
# println(Threads.nthreads())
# @btime begin    
#     pass_P!(D, r_arr, beam.nmacro, beam, zeros(Float64, 6), zeros(Float64, 6, 6))
# end

# function fE(len)
#     r = zeros(Float64, 1000, 6)
#     r[:, 2] .= 0.1
#     r_arr = matrix_to_array(r)
#     beam = Beam(r)
#     D = DRIFT(len=len)
#     pass_P!(D, r_arr, beam.nmacro, beam, zeros(Float64, 6), zeros(Float64, 6, 6))
#     return r_arr[1]
# end
# println(fE(3.0))

# using Enzyme
# Enzyme.API.runtimeActivity!(true) 
# grad = gradient(Forward, fE, [3.0])
# println(grad)

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