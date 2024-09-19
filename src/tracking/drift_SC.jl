function DriftPass_SC!(r_in::Array{Float64,1}, le::Float64, T1::Array{Float64,1}, T2::Array{Float64,1}, 
    R1::Array{Float64,2}, R2::Array{Float64, 2}, RApertures::Array{Float64,1}, EApertures::Array{Float64,1}, 
    num_particles::Int, lost_flags::Array{Int64,1}, a, b, Nl, Nm, K)
    # Ref [Qiang, Ji. "Differentiable self-consistent space-charge simulation for accelerator design." Physical Review Accelerators and Beams 26.2 (2023): 024601.]

    # half-kick-half
    for c in 1:num_particles
        if isone(lost_flags[c])
            continue
        end
        r6 = @view r_in[(c-1)*6+1:c*6]
        if !isnan(r6[1])
            # Misalignment at entrance
            if !iszero(T1)
                addvv!(r6, T1)
            end
            if !iszero(R1)
                multmv!(r6, R1)
            end
            drift6!(r6, le/2.0)
            # Misalignment at exit
            if !iszero(R2)
                multmv!(r6, R2)
            end
            if !iszero(T2)
                addvv!(r6, T2)
            end
            if check_lost(r6)
                lost_flags[c] = 1
            end
        end
    end

    space_charge!(r_in, K, Nl, Nm, a/Nl, b/Nm, a, b, num_particles, le, lost_flags)

    for c in 1:num_particles
        if isone(lost_flags[c])
            continue
        end
        r6 = @view r_in[(c-1)*6+1:c*6]
        if !isnan(r6[1])
            # Misalignment at entrance
            if !iszero(T1)
                addvv!(r6, T1)
            end
            if !iszero(R1)
                multmv!(r6, R1)
            end
            drift6!(r6, le/2.0)
            # Misalignment at exit
            if !iszero(R2)
                multmv!(r6, R2)
            end
            if !iszero(T2)
                addvv!(r6, T2)
            end
            if check_lost(r6)
                lost_flags[c] = 1
            end
        end
    end

    return nothing
end

function pass!(ele::DRIFT_SC, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam)
    # ele: EDRIFT
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    lost_flags = particles.lost_flag
    K = calculate_K(particles, particles.current)
    lstep = ele.len / ele.Nsteps
    for i in 1:ele.Nsteps
        DriftPass_SC!(r_in, lstep, ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, num_particles, lost_flags,
            ele.a, ele.b, ele.Nl, ele.Nm, K)
    end
    return nothing
end



################################################################################
# multi-threading
function DriftPass_SC_P!(r_in::Array{Float64,1}, le::Float64, T1::Array{Float64,1}, T2::Array{Float64,1}, 
    R1::Array{Float64,2}, R2::Array{Float64, 2}, RApertures::Array{Float64,1}, EApertures::Array{Float64,1}, 
    num_particles::Int, lost_flags::Array{Int64,1}, a, b, Nl, Nm, K)

    # half-kick-half
    Threads.@threads for c in 1:num_particles
        if isone(lost_flags[c])
            continue
        end
        r6 = @view r_in[(c-1)*6+1:c*6]
        if !isnan(r6[1])
            # Misalignment at entrance
            if !iszero(T1)
                addvv!(r6, T1)
            end
            if !iszero(R1)
                multmv!(r6, R1)
            end
            drift6!(r6, le/2.0)
            # Misalignment at exit
            if !iszero(R2)
                multmv!(r6, R2)
            end
            if !iszero(T2)
                addvv!(r6, T2)
            end
            if check_lost(r6)
                lost_flags[c] = 1
            end
        end
    end

    space_charge_P!(r_in, K, Nl, Nm, a/Nl, b/Nm, a, b, num_particles, le, lost_flags)

    Threads.@threads for c in 1:num_particles
        if isone(lost_flags[c])
            continue
        end
        r6 = @view r_in[(c-1)*6+1:c*6]
        if !isnan(r6[1])
            # Misalignment at entrance
            if !iszero(T1)
                addvv!(r6, T1)
            end
            if !iszero(R1)
                multmv!(r6, R1)
            end
            drift6!(r6, le/2.0)
            # Misalignment at exit
            if !iszero(R2)
                multmv!(r6, R2)
            end
            if !iszero(T2)
                addvv!(r6, T2)
            end
            if check_lost(r6)
                lost_flags[c] = 1
            end
        end
    end

    return nothing
end

function pass_P!(ele::DRIFT_SC, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam)
    # ele: EDRIFT
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    lost_flags = particles.lost_flag
    K = calculate_K(particles, particles.current)
    lstep = ele.len / ele.Nsteps
    for i in 1:ele.Nsteps
        DriftPass_SC_P!(r_in, lstep, ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, num_particles, lost_flags,
            ele.a, ele.b, ele.Nl, ele.Nm, K)
    end
    return nothing
end