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
    # in the loop if momentum deviation (delta) does not change
    # such as in 4-th order symplectic integrator w/o radiation

    if use_exact_Hamiltonian == 1
        r[1] += NormL * r[2]
        r[3] += NormL * r[4]
        r[5] += NormL * (1.0 + r[6]) - le
    else
        r[1] += NormL * r[2]
        r[3] += NormL * r[4]
        r[5] += NormL * (r[2]^2 + r[4]^2) / (2.0*(1.0+r[6]))
    end
    return nothing
end

function drift6!(r::AbstractVector{Float64}, le::Float64)
    # AT uses small angle approximation pz = 1 + delta. 
    # Here we use pz = sqrt((1 + delta)^2 - px^2 - py^2) for precise calculation
    if use_exact_Hamiltonian == 1
        NormL = le / sqrt(((1.0 + r[6])^2 - r[2]^2 - r[4]^2))
        r[5] += NormL * (1.0 + r[6]) - le
    else
        NormL = le / (1.0 + r[6])
        r[5] += NormL * (r[2]^2 + r[4]^2) / (2.0*(1.0+r[6])) # for linearized approximation
    end
    r[1] += NormL * r[2]
    r[3] += NormL * r[4]
    return nothing
end 
function DriftPass!(r_in::Array{Float64,1}, le::Float64, T1::Array{Float64,1}, T2::Array{Float64,1}, 
    R1::Array{Float64,2}, R2::Array{Float64, 2}, RApertures::Array{Float64,1}, EApertures::Array{Float64,1}, 
    num_particles::Int, lost_flags::Array{Int64,1})
    # Threads.@threads for c in 1:num_particles
    for c in 1:num_particles
        if lost_flags[c] == 1
            continue
        end
        r6 = @view r_in[(c-1)*6+1:c*6]
        if !isnan(r6[1])
            # Misalignment at entrance
            if T1 != zeros(6)
                ATaddvv!(r6, T1)
            end
            if R1 != zeros(6, 6)
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
            if R2 != zeros(6, 6)
                ATmultmv!(r6, R2)
            end
            if T2 != zeros(6)
                ATaddvv!(r6, T2)
            end
            if r6[1] > CoordLimit || r6[2] > AngleLimit || r6[1] < -CoordLimit || r6[2] < -AngleLimit
                lost_flags[c] = 1
            end
        end
    end
    return nothing
end

function pass!(ele::DRIFT, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam)
    # ele: EDRIFT
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    lost_flags = particles.lost_flag
    DriftPass!(r_in, ele.len, ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, num_particles, lost_flags)
    return nothing
end

function pass!(ele::MARKER, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam)
    return nothing
end

################################################################################
# multi-threading
function DriftPass_P!(r_in::Array{Float64,1}, le::Float64, T1::Array{Float64,1}, T2::Array{Float64,1}, 
    R1::Array{Float64,2}, R2::Array{Float64, 2}, RApertures::Array{Float64,1}, EApertures::Array{Float64,1}, 
    num_particles::Int, lost_flags::Array{Int64,1})
    Threads.@threads for c in 1:num_particles
    # for c in 1:num_particles
        if lost_flags[c] == 1
            continue
        end
        r6 = @view r_in[(c-1)*6+1:c*6]
        if !isnan(r6[1])
            # Misalignment at entrance
            if T1 != zeros(6)
                ATaddvv!(r6, T1)
            end
            if R1 != zeros(6, 6)
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
            if R2 != zeros(6, 6)
                ATmultmv!(r6, R2)
            end
            if T2 != zeros(6)
                ATaddvv!(r6, T2)
            end
            if r6[1] > CoordLimit || r6[2] > AngleLimit || r6[1] < -CoordLimit || r6[2] < -AngleLimit
                lost_flags[c] = 1
            end
        end
    end
    return nothing
end

function pass_P!(ele::DRIFT, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam)
    # ele: EDRIFT
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    lost_flags = particles.lost_flag
    DriftPass_P!(r_in, ele.len, ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, num_particles, lost_flags)
    return nothing
end

function pass_P!(ele::MARKER, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam)
    return nothing
end