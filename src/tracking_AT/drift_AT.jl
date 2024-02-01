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

function fastdrift!(r::AbstractVector{Float64}, NormL::Float64)
    # NormL=(Physical Length)/(1+delta)  is computed externally to speed up calculations
    # in the loop if momentum deviation (delta) does not change
    # such as in 4-th order symplectic integrator w/o radiation
    r[1] += NormL * r[2]
    r[3] += NormL * r[4]
    r[6] += NormL * (r[2]^2 + r[4]^2) / (2*(1+r[5]))
    return nothing
end

function drift6!(r::AbstractVector{Float64}, le::Float64)
    p_norm = 1.0 / (1.0 + r[5])
    NormL = le * p_norm
    r[1] += NormL * r[2]
    r[3] += NormL * r[4]
    r[6] += NormL * p_norm * (r[2] * r[2] + r[4] * r[4]) / 2.0
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
            if !isnothing(T1)
                ATaddvv!(r6, T1)
            end
            if !isnothing(R1)
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
            if !isnothing(R2)
                ATmultmv!(r6, R2)
            end
            if !isnothing(T2)
                ATaddvv!(r6, T2)
            end
            if r6[1] > CoordLimit || r6[2] > AngleLimit
                lost_flags[c] = 1
            end
        end
    end
    return nothing
end

function pass!(ele::DRIFT, r_in::Array{Float64,1}, num_particles::Int64, lost_flags::Array{Int64,1})
    # ele: EDRIFT
    # r_in: 6-by-num_particles array
    # num_particles: number of particles

    DriftPass!(r_in, ele.len, ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, num_particles, lost_flags)
    return nothing
end
