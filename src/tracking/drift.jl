function multmv!(r::AbstractVector{Float64}, A::Matrix{Float64})
    # multiplies 6-component column vector r by 6x6 matrix R: as in A*r
    @inbounds begin
        r1 = A[1, 1] * r[1] + A[1, 2] * r[2] + A[1, 3] * r[3] + A[1, 4] * r[4] + A[1, 5] * r[5] + A[1, 6] * r[6]
        r2 = A[2, 1] * r[1] + A[2, 2] * r[2] + A[2, 3] * r[3] + A[2, 4] * r[4] + A[2, 5] * r[5] + A[2, 6] * r[6]
        r3 = A[3, 1] * r[1] + A[3, 2] * r[2] + A[3, 3] * r[3] + A[3, 4] * r[4] + A[3, 5] * r[5] + A[3, 6] * r[6]
        r4 = A[4, 1] * r[1] + A[4, 2] * r[2] + A[4, 3] * r[3] + A[4, 4] * r[4] + A[4, 5] * r[5] + A[4, 6] * r[6]
        r5 = A[5, 1] * r[1] + A[5, 2] * r[2] + A[5, 3] * r[3] + A[5, 4] * r[4] + A[5, 5] * r[5] + A[5, 6] * r[6]
        r6 = A[6, 1] * r[1] + A[6, 2] * r[2] + A[6, 3] * r[3] + A[6, 4] * r[4] + A[6, 5] * r[5] + A[6, 6] * r[6]
        r[1] = r1
        r[2] = r2
        r[3] = r3
        r[4] = r4
        r[5] = r5
        r[6] = r6
    end
    return nothing
end

function addvv!(r::AbstractVector{Float64}, dr::Array{Float64,1})
    # adds 6-component column vector dr to 6-component column vector r: as in r = r + dr
    for i in 1:6
        r[i] += dr[i]
    end
    return nothing
end

@inline function addvv_row!(r::Matrix{Float64}, c::Int, dr::Vector{Float64})
    @inbounds for i in 1:6
        r[c, i] += dr[i]
    end
    return nothing
end

@inline function multmv_row!(r::Matrix{Float64}, c::Int, A::Matrix{Float64})
    @inbounds begin
        r1 = A[1, 1] * r[c, 1] + A[1, 2] * r[c, 2] + A[1, 3] * r[c, 3] + A[1, 4] * r[c, 4] + A[1, 5] * r[c, 5] + A[1, 6] * r[c, 6]
        r2 = A[2, 1] * r[c, 1] + A[2, 2] * r[c, 2] + A[2, 3] * r[c, 3] + A[2, 4] * r[c, 4] + A[2, 5] * r[c, 5] + A[2, 6] * r[c, 6]
        r3 = A[3, 1] * r[c, 1] + A[3, 2] * r[c, 2] + A[3, 3] * r[c, 3] + A[3, 4] * r[c, 4] + A[3, 5] * r[c, 5] + A[3, 6] * r[c, 6]
        r4 = A[4, 1] * r[c, 1] + A[4, 2] * r[c, 2] + A[4, 3] * r[c, 3] + A[4, 4] * r[c, 4] + A[4, 5] * r[c, 5] + A[4, 6] * r[c, 6]
        r5 = A[5, 1] * r[c, 1] + A[5, 2] * r[c, 2] + A[5, 3] * r[c, 3] + A[5, 4] * r[c, 4] + A[5, 5] * r[c, 5] + A[5, 6] * r[c, 6]
        r6 = A[6, 1] * r[c, 1] + A[6, 2] * r[c, 2] + A[6, 3] * r[c, 3] + A[6, 4] * r[c, 4] + A[6, 5] * r[c, 5] + A[6, 6] * r[c, 6]
        r[c, 1] = r1
        r[c, 2] = r2
        r[c, 3] = r3
        r[c, 4] = r4
        r[c, 5] = r5
        r[c, 6] = r6
    end
    return nothing
end

@inline function _linearized_drift_phifac(px, py, dp_p, gamma2i::Float64)
    if use_exact_Hamiltonian == 2
        phifac = (px^2 + py^2 + dp_p^2 * gamma2i) / 2.0
        return (phifac / (1.0 + dp_p) - dp_p * gamma2i) / (1.0 + dp_p)
    end
    return (px^2 + py^2) / (2.0 * (1.0 + dp_p)^2)
end

function fastdrift!(r::AbstractVector{Float64}, NormL::Float64, le::Float64, beti::Float64, gamma2i::Float64 = 0.0)
    # Provide an option to use exact Hamiltonian or linearized approximation
    # AT uses small angle approximation pz = 1 + delta. 
    # MADX use pz = sqrt((1 + 2*delta/beta + delta^2 - px^2 - py^2).

    if use_exact_Hamiltonian == 1
        r[1] += NormL * r[2]
        r[3] += NormL * r[4]
        r[5] += NormL * (1.0*beti + r[6]) - le*beti
    else
        r[1] += NormL * r[2]
        r[3] += NormL * r[4]
        r[5] += le * _linearized_drift_phifac(r[2], r[4], r[6], gamma2i)
    end
    return nothing
end

function drift6!(r::AbstractVector{Float64}, le::Float64, beti::Float64, gamma2i::Float64 = 0.0)
    # Provide an option to use exact Hamiltonian or linearized approximation
    # AT uses small angle approximation pz = 1 + delta. 
    # MADX use pz = sqrt((1 + 2*delta/beta + delta^2 - px^2 - py^2).
    if use_exact_Hamiltonian == 1
        if 1.0 + 2.0*r[6]*beti + r[6]^2 - r[2]^2 - r[4]^2 <= 0.0
            # This is a special case when the particle is lost.
            # We set the particle to NaN to indicate it is lost.
            r[1] = NaN
            r[2] = NaN
            r[3] = NaN
            r[4] = NaN
            r[5] = NaN
            r[6] = NaN
            return nothing
        end
        NormL = le / sqrt(1.0 + 2.0*r[6]*beti + r[6]^2 - r[2]^2 - r[4]^2)
        r[5] += NormL * (1.0*beti + r[6]) - le*beti
    else
        NormL = le / (1.0 + r[6])
        r[5] += le * _linearized_drift_phifac(r[2], r[4], r[6], gamma2i)
    end
    r[1] += NormL * r[2]
    r[3] += NormL * r[4]
    return nothing
end 

@inline function drift6_row!(r::Matrix{Float64}, c::Int, le::Float64, beti::Float64, gamma2i::Float64 = 0.0)
    @inbounds begin
        px = r[c, 2]
        py = r[c, 4]
        dp_p = r[c, 6]
        NormL = 0.0
        if use_exact_Hamiltonian == 1
            discr = 1.0 + 2.0 * dp_p * beti + dp_p^2 - px^2 - py^2
            if discr <= 0.0
                r[c, 1] = NaN
                r[c, 2] = NaN
                r[c, 3] = NaN
                r[c, 4] = NaN
                r[c, 5] = NaN
                r[c, 6] = NaN
                return nothing
            end
            NormL = le / sqrt(discr)
            r[c, 5] += NormL * (beti + dp_p) - le * beti
        else
            NormL = le / (1.0 + dp_p)
            r[c, 5] += le * _linearized_drift_phifac(px, py, dp_p, gamma2i)
        end
        r[c, 1] += NormL * px
        r[c, 3] += NormL * py
    end
    return nothing
end

@inline function _check_lost_row(r::Matrix{Float64}, c::Int)
    @inbounds begin
        x = r[c, 1]
        px = r[c, 2]
        y = r[c, 3]
        py = r[c, 4]
        dp_p = r[c, 6]
        if isnan(x) || isinf(x)
            return true
        end
        if max(abs(x), abs(px), abs(y), abs(py)) > CoordLimit || abs(dp_p) > CoordLimit
            return true
        end
        if px^2 + py^2 > 1.0 + dp_p^2
            return true
        end
    end
    return false
end

@inline function _check_lost_aperture_row(r::Matrix{Float64}, c::Int, RApertures::Vector{Float64}, EApertures::Vector{Float64})
    @inbounds begin
        if !iszero(RApertures)
            x = r[c, 1]
            y = r[c, 3]
            if x < RApertures[1] || x > RApertures[2] || y < RApertures[3] || y > RApertures[4]
                return true
            end
        end
        if !iszero(EApertures)
            ax = EApertures[1]
            ay = EApertures[2]
            if ax > 0.0 && ay > 0.0
                x = r[c, 1]
                y = r[c, 3]
                if x^2 / ax^2 + y^2 / ay^2 > 1.0
                    return true
                end
            end
        end
    end
    return false
end

function DriftPass!(r_in::Matrix{Float64}, le::Float64, beti::Float64, gamma2i::Float64, T1::Array{Float64,1}, T2::Array{Float64,1}, 
    R1::Array{Float64,2}, R2::Array{Float64, 2}, RApertures::Array{Float64,1}, EApertures::Array{Float64,1}, 
    num_particles::Int, lost_flags::Array{Int64,1})
    # Modified based on AT function. Ref[Terebilo, Andrei. "Accelerator modeling with MATLAB accelerator toolbox." PACS2001 (2001)].
    for c in 1:num_particles
        if isone(lost_flags[c]) || isnan(r_in[c, 1])
            continue
        end
        if !iszero(T1)
            addvv_row!(r_in, c, T1)
        end
        if !iszero(R1)
            multmv_row!(r_in, c, R1)
        end

        drift6_row!(r_in, c, le, beti, gamma2i)

        if !iszero(R2)
            multmv_row!(r_in, c, R2)
        end
        if !iszero(T2)
            addvv_row!(r_in, c, T2)
        end
        if _check_lost_row(r_in, c) || _check_lost_aperture_row(r_in, c, RApertures, EApertures)
            lost_flags[c] = 1
        end
    end
    return nothing
end

"""
    pass!(ele::DRIFT, r_in::Matrix{Float64}, num_particles::Int64, particles::Beam{Float64})

This is a function to track particles through a drift element.

# Arguments
- ele::DRIFT: a drift element
- r_in::Matrix{Float64}: num_particles-by-6 matrix
- num_particles::Int64: number of particles
- particles::Beam{Float64}: beam object
"""
function pass!(ele::DRIFT, r_in::Matrix{Float64}, num_particles::Int64, particles::Beam{Float64})
    lost_flags = particles.lost_flag
    if use_exact_beti == 1
        beti = 1.0 / particles.beta
    else
        beti = 1.0 
    end
    gamma2i = use_exact_Hamiltonian == 2 ? 1.0 / (particles.gamma^2) : 0.0
    DriftPass!(r_in, ele.len, beti, gamma2i, ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, num_particles, lost_flags)
    return nothing
end

function pass!(ele::MARKER, r_in::Matrix{Float64}, num_particles::Int64, particles::Beam{Float64})
    return nothing
end

################################################################################
# multi-threading
function DriftPass_P!(r_in::Matrix{Float64}, le::Float64, beti::Float64, gamma2i::Float64, T1::Array{Float64,1}, T2::Array{Float64,1}, 
    R1::Array{Float64,2}, R2::Array{Float64, 2}, RApertures::Array{Float64,1}, EApertures::Array{Float64,1}, 
    num_particles::Int, lost_flags::Array{Int64,1})
    Threads.@threads for c in 1:num_particles
        if isone(lost_flags[c]) || isnan(r_in[c, 1])
            continue
        end
        if !iszero(T1)
            addvv_row!(r_in, c, T1)
        end
        if !iszero(R1)
            multmv_row!(r_in, c, R1)
        end

        drift6_row!(r_in, c, le, beti, gamma2i)

        if !iszero(R2)
            multmv_row!(r_in, c, R2)
        end
        if !iszero(T2)
            addvv_row!(r_in, c, T2)
        end
        if _check_lost_row(r_in, c) || _check_lost_aperture_row(r_in, c, RApertures, EApertures)
            lost_flags[c] = 1
        end
    end
    return nothing
end

function pass_P!(ele::DRIFT, r_in::Matrix{Float64}, num_particles::Int64, particles::Beam{Float64})
    # ele: EDRIFT
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    lost_flags = particles.lost_flag
    if use_exact_beti == 1
        beti = 1.0 / particles.beta
    else
        beti = 1.0 
    end
    gamma2i = use_exact_Hamiltonian == 2 ? 1.0 / (particles.gamma^2) : 0.0
    DriftPass_P!(r_in, ele.len, beti, gamma2i, ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, num_particles, lost_flags)
    return nothing
end

function pass_P!(ele::MARKER, r_in::Matrix{Float64}, num_particles::Int64, particles::Beam{Float64})
    return nothing
end

##################################
# Space charge
####################################
function DriftPass_SC!(r_in::Matrix{Float64}, le::Float64, beti::Float64, gamma2i::Float64, T1::Array{Float64,1}, T2::Array{Float64,1}, 
    R1::Array{Float64,2}, R2::Array{Float64, 2}, RApertures::Array{Float64,1}, EApertures::Array{Float64,1}, 
    num_particles::Int, lost_flags::Array{Int64,1}, a, b, Nl, Nm, K)
    # Ref [Qiang, Ji. "Differentiable self-consistent space-charge simulation for accelerator design." Physical Review Accelerators and Beams 26.2 (2023): 024601.]

    # half-kick-half
    for c in 1:num_particles
        if isone(lost_flags[c])
            continue
        end
        r6 = @view r_in[c, :]
        if !isnan(r6[1])
            # Misalignment at entrance
            if !iszero(T1)
                addvv!(r6, T1)
            end
            if !iszero(R1)
                multmv!(r6, R1)
            end
            drift6!(r6, le/2.0, beti, gamma2i)
            # Misalignment at exit
            if !iszero(R2)
                multmv!(r6, R2)
            end
            if !iszero(T2)
                addvv!(r6, T2)
            end
            if check_lost(r6) || check_lost_aperture(r6, RApertures, EApertures)
                lost_flags[c] = 1
            end
        end
    end

    space_charge!(r_in, K, Nl, Nm, a/Nl, b/Nm, a, b, num_particles, le, lost_flags)

    for c in 1:num_particles
        if isone(lost_flags[c])
            continue
        end
        r6 = @view r_in[c, :]
        if !isnan(r6[1])
            # Misalignment at entrance
            if !iszero(T1)
                addvv!(r6, T1)
            end
            if !iszero(R1)
                multmv!(r6, R1)
            end
            drift6!(r6, le/2.0, beti, gamma2i)
            # Misalignment at exit
            if !iszero(R2)
                multmv!(r6, R2)
            end
            if !iszero(T2)
                addvv!(r6, T2)
            end
            if check_lost(r6) || check_lost_aperture(r6, RApertures, EApertures)
                lost_flags[c] = 1
            end
        end
    end

    return nothing
end

function DriftPass_SC!(r_in::Matrix{DTPSAD{N, T}}, le::DTPSAD{N, T}, beti::Float64, gamma2i::Float64, T1::Array{DTPSAD{N, T},1}, T2::Array{DTPSAD{N, T},1}, 
    R1::Array{DTPSAD{N, T},2}, R2::Array{DTPSAD{N, T}, 2}, RApertures::Array{Float64,1}, EApertures::Array{Float64,1}, 
    num_particles::Int, lost_flags::Array{Int64,1}, a::DTPSAD{N, T}, b::DTPSAD{N, T}, Nl::Int64, Nm::Int64, K::DTPSAD{N, T}) where {N, T}
    # Ref [Qiang, Ji. "Differentiable self-consistent space-charge simulation for accelerator design." Physical Review Accelerators and Beams 26.2 (2023): 024601.]

    # half-kick-half
    @inbounds for c in 1:num_particles
        if isone(lost_flags[c])
            continue
        end
        r6 = @view r_in[c, :]
        if !isnan(r6[1])
            # Misalignment at entrance
            if !iszero(T1)
                addvv!(r6, T1)
            end
            if !iszero(R1)
                multmv!(r6, R1)
            end
            drift6!(r6, le/2.0, gamma2i)
            # Misalignment at exit
            if !iszero(R2)
                multmv!(r6, R2)
            end
            if !iszero(T2)
                addvv!(r6, T2)
            end
            if check_lost(r6) || check_lost_aperture(r6, RApertures, EApertures)
                lost_flags[c] = 1
            end
        end
    end

    space_charge!(r_in, K, Nl, Nm, a/Nl, b/Nm, a, b, num_particles, le, lost_flags)

    for c in 1:num_particles
        if isone(lost_flags[c])
            continue
        end
        r6 = @view r_in[c, :]
        if !isnan(r6[1])
            # Misalignment at entrance
            if !iszero(T1)
                addvv!(r6, T1)
            end
            if !iszero(R1)
                multmv!(r6, R1)
            end
            drift6!(r6, le/2.0, gamma2i)
            # Misalignment at exit
            if !iszero(R2)
                multmv!(r6, R2)
            end
            if !iszero(T2)
                addvv!(r6, T2)
            end
            if check_lost(r6) || check_lost_aperture(r6, RApertures, EApertures)
                lost_flags[c] = 1
            end
        end
    end

    return nothing
end

function pass!(ele::DRIFT_SC{Float64}, r_in::Matrix{Float64}, num_particles::Int64, particles::Beam{Float64})
    # ele: EDRIFT
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    lost_flags = particles.lost_flag
    K = calculate_K(particles, particles.current)
    lstep = ele.len / ele.Nsteps
    if use_exact_beti == 1
        beti = 1.0 / particles.beta
    else
        beti = 1.0 
    end
    for i in 1:ele.Nsteps
        gamma2i = use_exact_Hamiltonian == 2 ? 1.0 / (particles.gamma^2) : 0.0
        DriftPass_SC!(r_in, lstep, beti, gamma2i, ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, num_particles, lost_flags,
            ele.a, ele.b, ele.Nl, ele.Nm, K)
    end
    return nothing
end

function pass!(ele::DRIFT_SC{DTPSAD{N, T}}, r_in::Matrix{DTPSAD{N, T}}, num_particles::Int64, particles::Beam{DTPSAD{N, T}}) where {N, T}
    # ele: EDRIFT
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    lost_flags = particles.lost_flag
    K = calculate_K(particles, particles.current)
    lstep = ele.len / ele.Nsteps
    if use_exact_beti == 1
        beti = 1.0 / particles.beta.val
    else
        beti = 1.0 
    end
    for i in 1:ele.Nsteps
        gamma2i = use_exact_Hamiltonian == 2 ? 1.0 / particles.gamma.val^2 : 0.0
        DriftPass_SC!(r_in, lstep, beti, gamma2i, ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, num_particles, lost_flags,
            ele.a, ele.b, ele.Nl, ele.Nm, K)
    end
    return nothing
end

function DriftPass_SC2P5D!(ele::DRIFT_SC2P5D{Float64}, r_in::Matrix{Float64}, le::Float64, beti::Float64,
    gamma2i::Float64,
    num_particles::Int, lost_flags::Array{Int64,1}, particles::Beam{Float64})
    lstep = le / ele.Nsteps

    for c in 1:num_particles
        if isone(lost_flags[c])
            continue
        end
        if !isnan(r_in[c, 1])
            if !iszero(ele.T1)
                addvv_row!(r_in, c, ele.T1)
            end
            if !iszero(ele.R1)
                multmv_row!(r_in, c, ele.R1)
            end
            if _check_lost_row(r_in, c) || _check_lost_aperture_row(r_in, c, ele.RApertures, ele.EApertures)
                lost_flags[c] = 1
            end
        end
    end

    for _ in 1:ele.Nsteps
        _sc2p5d_track_step!(ele, lstep, r_in, particles)

        for c in 1:num_particles
            if isone(lost_flags[c])
                continue
            end
            if !isnan(r_in[c, 1])
                drift6_row!(r_in, c, lstep, beti, gamma2i)
                if _check_lost_row(r_in, c) || _check_lost_aperture_row(r_in, c, ele.RApertures, ele.EApertures)
                    lost_flags[c] = 1
                end
            end
        end
    end

    for c in 1:num_particles
        if isone(lost_flags[c])
            continue
        end
        if !isnan(r_in[c, 1])
            if !iszero(ele.R2)
                multmv_row!(r_in, c, ele.R2)
            end
            if !iszero(ele.T2)
                addvv_row!(r_in, c, ele.T2)
            end
            if _check_lost_row(r_in, c) || _check_lost_aperture_row(r_in, c, ele.RApertures, ele.EApertures)
                lost_flags[c] = 1
            end
        end
    end

    return nothing
end

function DriftPass_SC2P5D!(ele::DRIFT_SC2P5D{DTPSAD{N, T}}, r_in::Matrix{DTPSAD{N, T}}, le::DTPSAD{N, T},
    beti::Float64, gamma2i::Float64, num_particles::Int, lost_flags::Array{Int64,1}, particles::Beam{DTPSAD{N, T}}) where {N, T}
    lstep = le / ele.Nsteps

    @inbounds for c in 1:num_particles
        if isone(lost_flags[c])
            continue
        end
        r6 = @view r_in[c, :]
        if !isnan(r6[1])
            if !iszero(ele.T1)
                addvv!(r6, ele.T1)
            end
            if !iszero(ele.R1)
                multmv!(r6, ele.R1)
            end
            if check_lost(r6) || check_lost_aperture(r6, ele.RApertures, ele.EApertures)
                lost_flags[c] = 1
            end
        end
    end

    for _ in 1:ele.Nsteps
        _sc2p5d_track_step!(ele, lstep, r_in, particles)

        @inbounds for c in 1:num_particles
            if isone(lost_flags[c])
                continue
            end
            r6 = @view r_in[c, :]
            if !isnan(r6[1])
                drift6!(r6, lstep, gamma2i)
                if check_lost(r6) || check_lost_aperture(r6, ele.RApertures, ele.EApertures)
                    lost_flags[c] = 1
                end
            end
        end
    end

    @inbounds for c in 1:num_particles
        if isone(lost_flags[c])
            continue
        end
        r6 = @view r_in[c, :]
        if !isnan(r6[1])
            if !iszero(ele.R2)
                multmv!(r6, ele.R2)
            end
            if !iszero(ele.T2)
                addvv!(r6, ele.T2)
            end
            if check_lost(r6) || check_lost_aperture(r6, ele.RApertures, ele.EApertures)
                lost_flags[c] = 1
            end
        end
    end

    return nothing
end

function pass!(ele::DRIFT_SC2P5D{Float64}, r_in::Matrix{Float64}, num_particles::Int64, particles::Beam{Float64})
    lost_flags = particles.lost_flag
    if use_exact_beti == 1
        beti = 1.0 / particles.beta
    else
        beti = 1.0
    end
    gamma2i = use_exact_Hamiltonian == 2 ? 1.0 / (particles.gamma^2) : 0.0
    DriftPass_SC2P5D!(ele, r_in, ele.len, beti, gamma2i, num_particles, lost_flags, particles)
    return nothing
end

function pass!(ele::DRIFT_SC2P5D{DTPSAD{N, T}}, r_in::Matrix{DTPSAD{N, T}}, num_particles::Int64,
    particles::Beam{DTPSAD{N, T}}) where {N, T}
    lost_flags = particles.lost_flag
    if use_exact_beti == 1
        beti = 1.0 / particles.beta.val
    else
        beti = 1.0
    end
    gamma2i = use_exact_Hamiltonian == 2 ? 1.0 / particles.gamma.val^2 : 0.0
    DriftPass_SC2P5D!(ele, r_in, ele.len, beti, gamma2i, num_particles, lost_flags, particles)
    return nothing
end


################################################################################
# multi-threading
function DriftPass_SC_P!(r_in::Matrix{Float64}, le::Float64, beti::Float64, gamma2i::Float64, T1::Array{Float64,1}, T2::Array{Float64,1}, 
    R1::Array{Float64,2}, R2::Array{Float64, 2}, RApertures::Array{Float64,1}, EApertures::Array{Float64,1}, 
    num_particles::Int, lost_flags::Array{Int64,1}, a, b, Nl, Nm, K)

    # half-kick-half
    Threads.@threads for c in 1:num_particles
        if isone(lost_flags[c])
            continue
        end
        r6 = @view r_in[c, :]
        if !isnan(r6[1])
            # Misalignment at entrance
            if !iszero(T1)
                addvv!(r6, T1)
            end
            if !iszero(R1)
                multmv!(r6, R1)
            end
            drift6!(r6, le/2.0, beti, gamma2i)
            # Misalignment at exit
            if !iszero(R2)
                multmv!(r6, R2)
            end
            if !iszero(T2)
                addvv!(r6, T2)
            end
            if check_lost(r6) || check_lost_aperture(r6, RApertures, EApertures)
                lost_flags[c] = 1
            end
        end
    end

    space_charge_P!(r_in, K, Nl, Nm, a/Nl, b/Nm, a, b, num_particles, le, lost_flags)

    Threads.@threads for c in 1:num_particles
        if isone(lost_flags[c])
            continue
        end
        r6 = @view r_in[c, :]
        if !isnan(r6[1])
            # Misalignment at entrance
            if !iszero(T1)
                addvv!(r6, T1)
            end
            if !iszero(R1)
                multmv!(r6, R1)
            end
            drift6!(r6, le/2.0, beti, gamma2i)
            # Misalignment at exit
            if !iszero(R2)
                multmv!(r6, R2)
            end
            if !iszero(T2)
                addvv!(r6, T2)
            end
            if check_lost(r6) || check_lost_aperture(r6, RApertures, EApertures)
                lost_flags[c] = 1
            end
        end
    end

    return nothing
end

function pass_P!(ele::DRIFT_SC, r_in::Matrix{Float64}, num_particles::Int64, particles::Beam{Float64})
    # ele: EDRIFT
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    lost_flags = particles.lost_flag
    K = calculate_K(particles, particles.current)
    lstep = ele.len / ele.Nsteps
    if use_exact_beti == 1
        beti = 1.0 / particles.beta
    else
        beti = 1.0 
    end
    gamma2i = use_exact_Hamiltonian == 2 ? 1.0 / (particles.gamma^2) : 0.0
    for i in 1:ele.Nsteps
        DriftPass_SC_P!(r_in, lstep, beti, gamma2i, ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, num_particles, lost_flags,
            ele.a, ele.b, ele.Nl, ele.Nm, K)
    end
    return nothing
end

function pass_P!(ele::DRIFT_SC2P5D, r_in::Matrix{Float64}, num_particles::Int64, particles::Beam{Float64})
    pass!(ele, r_in, num_particles, particles)
    return nothing
end

################################################################################
# High-order TPSA
################################################################################
function multmv!(r::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}, A::Matrix{Float64}) where {T, TPS_Dim, Max_TPS_Degree}
    # multiplies 6-component column vector r by 6x6 matrix R: as in A*r
    temp = Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}(undef, 6)

    for i in 1:6
        temp[i] = CTPS(0.0, TPS_Dim, Max_TPS_Degree)
        for j in 1:6
            temp[i] += A[i, j] * r[j]
        end
    end
    r .= temp
    return nothing
end

function addvv!(r::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}, dr::Array{Float64, 1}) where {T, TPS_Dim, Max_TPS_Degree}
    for i in 1:6
        r[i] += dr[i]
    end
    return nothing
end

function fastdrift!(r::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}, NormL::CTPS{T, TPS_Dim, Max_TPS_Degree}, 
    le::Float64, beti::Float64, gamma2i::Float64 = 0.0) where {T, TPS_Dim, Max_TPS_Degree}
    # Provide an option to use exact Hamiltonian or linearized approximation
    # AT uses small angle approximation pz = 1 + delta. 
    # MADX use pz = sqrt((1 + 2*delta/beta + delta^2 - px^2 - py^2).
    if use_exact_Hamiltonian == 1
        r[1] += NormL * r[2]
        r[3] += NormL * r[4]
        r[5] += NormL * (1.0*beti + r[6]) - le*beti
    else
        r[1] += NormL * r[2]
        r[3] += NormL * r[4]
        r[5] += le * _linearized_drift_phifac(r[2], r[4], r[6], gamma2i)
    end
    return nothing 
end

function drift6!(r::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}, le::Float64, beti::Float64,
    gamma2i::Float64 = 0.0) where {T, TPS_Dim, Max_TPS_Degree}
    # Provide an option to use exact Hamiltonian or linearized approximation
    # AT uses small angle approximation pz = 1 + delta. 
    # MADX use pz = sqrt((1 + 2*delta/beta + delta^2 - px^2 - py^2).
    if use_exact_Hamiltonian == 1
        NormL = le / sqrt(1.0 + 2.0*r[6]*beti + r[6]^2 - r[2]^2 - r[4]^2)
        r[5] += NormL * (1.0*beti + r[6]) - le*beti
    else
        NormL = le / (1.0 + r[6])
        r[5] += le * _linearized_drift_phifac(r[2], r[4], r[6], gamma2i)
    end
    r[1] += NormL * r[2]
    r[3] += NormL * r[4]
    return nothing
end 
function DriftPass_TPSA!(r_in::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}, le::Float64, beti::Float64, gamma2i::Float64,
    T1::Array{Float64,1}, T2::Array{Float64,1}, R1::Array{Float64,2}, R2::Array{Float64,2}, 
    RApertures::Array{Float64,1}, EApertures::Array{Float64,1}) where {T, TPS_Dim, Max_TPS_Degree}

    if !iszero(T1)
        addvv!(r_in, T1)
    end
    if !iszero(R1)
        multmv!(r_in, R1)
    end

    drift6!(r_in, le, beti, gamma2i)

    # Misalignment at exit
    if !iszero(R2)
        multmv!(r_in, R2)
    end
    if !iszero(T2)
        addvv!(r_in, T2)
    end

    return nothing
end

function pass_TPSA!(ele::DRIFT, r_in::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}; E0::Float64=0.0, m0::Float64=m_e) where {T, TPS_Dim, Max_TPS_Degree}
    # ele: EDRIFT
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    gamma = (E0 + m0) / m0
    beta = sqrt(1.0 - 1.0 / gamma^2)
    if use_exact_beti == 1
        beti = 1.0 / beta
    else
        beti = 1.0 
    end
    gamma2i = use_exact_Hamiltonian == 2 ? 1.0 / gamma^2 : 0.0
    DriftPass_TPSA!(r_in, ele.len, beti, gamma2i, ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures)
    return nothing
end

function pass_TPSA!(ele::MARKER, r_in::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}; E0::Float64=0.0, m0::Float64=m_e) where {T, TPS_Dim, Max_TPS_Degree}
    return nothing
end
