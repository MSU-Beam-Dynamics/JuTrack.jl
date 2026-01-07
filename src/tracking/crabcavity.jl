function CRABCAVITYPass!(r::Matrix{Float64}, cavity::CRABCAVITY, beta::Float64, E::Float64, num_particles::Int64, lost_flags::Array{Int,1})
    if use_exact_beti == 1
        beti = 1.0 / beta
    else
        beti = 1.0 
    end
    for c in 1:num_particles
        if isone(lost_flags[c])
            continue
        end
        r6 = @view r[c, :] 
        ang = cavity.k * r6[5] + cavity.phi
        if cavity.len == 0.0
            # r6[2] += (cavity.volt/beta2E) * sin(ang)
            # r6[6] += (-cavity.k * cavity.volt/beta2E) * r6[1] * cos(ang)
            # this is symplectic form
            r6[2] += (cavity.volt/E) * sin(ang/beta)
            r6[6] += (-cavity.k * cavity.volt/E/beta) * r6[1] * cos(ang/beta)
        else
            drift6!(r6, cavity.len / 2.0, beti)
            # r6[2] += (cavity.volt/beta2E) * sin(ang)
            # r6[6] += (-cavity.k * cavity.volt/beta2E) * r6[1] * cos(ang)
            # this is symplectic form
            r6[2] += (cavity.volt/E) * sin(ang/beta)
            r6[6] += (-cavity.k * cavity.volt/E/beta) * r6[1] * cos(ang/beta)
            drift6!(r6, cavity.len / 2.0, beti)
        end
    end
    return nothing
end

function CRABCAVITYK2Pass!(r::Matrix{Float64}, cavity::CRABCAVITY_K2, beta::Float64, E::Float64, num_particles::Int64, lost_flags::Array{Int,1})
    if use_exact_beti == 1
        beti = 1.0 / beta
    else
        beti = 1.0 
    end
    for c in 1:num_particles
        if isone(lost_flags[c])
            continue
        end
        r6 = @view r[c, :] 
        ang = cavity.k * r6[5] + cavity.phi
        if cavity.len == 0.0
            # r6[2] += (cavity.volt/beta2E) * sin(ang)
            # r6[6] += (-cavity.k * cavity.volt/beta2E) * r6[1] * cos(ang)
            # this is symplectic form
            r6[2] += (cavity.volt/E) * sin(ang/beta)
            r6[6] += (-cavity.k * cavity.volt/E/beta) * r6[1] * cos(ang/beta)
            # sextupole kick
            r6[2] -= cavity.k2 * (r6[1]^2 - r6[3]^2) * sin(ang)
            r6[4] += 2.0 * cavity.k2 * r6[1] * r6[3] * sin(ang)
            r6[6] -= (cavity.k2 * cavity.k / 3.0) * (r6[1]^3 - 3.0 * r6[1] * r6[3]^2) * cos(ang)
        else
            drift6!(r6, cavity.len / 2.0, beti)
            # r6[2] += (cavity.volt/beta2E) * sin(ang)
            # r6[6] += (-cavity.k * cavity.volt/beta2E) * r6[1] * cos(ang)            
            # this is symplectic form
            r6[2] += (cavity.volt/E) * sin(ang/beta)
            r6[6] += (-cavity.k * cavity.volt/E/beta) * r6[1] * cos(ang/beta)
            # sextupole kick
            r6[2] -= cavity.k2 * (r6[1]^2 - r6[3]^2) * sin(ang)
            r6[4] += 2.0 * cavity.k2 * r6[1] * r6[3] * sin(ang)
            r6[6] -= (cavity.k2 * cavity.k / 3.0) * (r6[1]^3 - 3.0 * r6[1] * r6[3]^2) * cos(ang)
            drift6!(r6, cavity.len / 2.0, beti)
        end
    end
    return nothing
end


function easyCRABCAVITYPass!(r::Matrix{Float64}, cavity::easyCRABCAVITY, num_particles::Int64, lost_flags::Array{Int,1})
    for c in 1:num_particles
        if isone(lost_flags[c])
            continue
        end
        r6 = @view r[c, :] 
        ang = cavity.k * r6[5] + cavity.phi
        r6[1] += (cavity.halfthetac/cavity.k) * sin(ang)
        r6[6] += cavity.halfthetac * r6[2] * cos(ang)
    end
    return nothing
end

function AccelCavityPass!(r::Matrix{Float64}, cavity::AccelCavity, beta2E::Float64, num_particles::Int64, lost_flags::Array{Int,1})
    v_beta2E = cavity.volt/beta2E
    # Threads.@threads for c in 1:num_particles
    for c in 1:num_particles
        if lost_flags[c] == 1
            continue
        end
        r6 = @view r[c, :] 
        sv = sin(cavity.k * r6[5] + cavity.phis) - sin(cavity.phis)
        r6[6] += v_beta2E * sv
    end
    return nothing
end

function pass!(ele::CRABCAVITY, r_in::Matrix{Float64}, num_particles::Int64, particles::Beam{Float64})
    lost_flags = particles.lost_flag
    CRABCAVITYPass!(r_in, ele, particles.beta, particles.energy, num_particles,lost_flags)
    return nothing
end

function pass!(ele::CRABCAVITY_K2, r_in::Matrix{Float64}, num_particles::Int64, particles::Beam{Float64})
    lost_flags = particles.lost_flag
    CRABCAVITYK2Pass!(r_in, ele, particles.beta, particles.energy, num_particles,lost_flags)
    return nothing
end

function pass!(ele::easyCRABCAVITY, r_in::Matrix{Float64}, num_particles::Int64, particles::Beam{Float64})
    lost_flags = particles.lost_flag
    easyCRABCAVITYPass!(r_in, ele, num_particles,lost_flags)
    return nothing
end

function pass!(ele::AccelCavity, r_in::Matrix{Float64}, num_particles::Int64, particles::Beam{Float64})
    lost_flags = particles.lost_flag
    beta2E = particles.beta^2 * particles.energy
    AccelCavityPass!(r_in, ele, beta2E, num_particles,lost_flags)
    return nothing
end

##########################################################################################
# multi-threading
function CRABCAVITYPass_P!(r::Matrix{Float64}, cavity::CRABCAVITY, beta::Float64, E::Float64, num_particles::Int64, lost_flags::Array{Int,1})
    if use_exact_beti == 1
        beti = 1.0 / beta
    else
        beti = 1.0 
    end
    Threads.@threads for c in 1:num_particles
    # for c in 1:num_particles
        if isone(lost_flags[c])
            continue
        end
        r6 = @view r[c, :] 
        ang = cavity.k * r6[5] + cavity.phi
        if cavity.len == 0.0
            # r6[2] += (cavity.volt/beta2E) * sin(ang)
            # r6[6] += (-cavity.k * cavity.volt/beta2E) * r6[1] * cos(ang)
            # this is symplectic form
            r6[2] += (cavity.volt/E) * sin(ang/beta)
            r6[6] += (-cavity.k * cavity.volt/E/beta) * r6[1] * cos(ang/beta)
        else
            drift6!(r6, cavity.len / 2.0, beti)
            # r6[2] += (cavity.volt/beta2E) * sin(ang)
            # r6[6] += (-cavity.k * cavity.volt/beta2E) * r6[1] * cos(ang)
            # this is symplectic form
            r6[2] += (cavity.volt/E) * sin(ang/beta)
            r6[6] += (-cavity.k * cavity.volt/E/beta) * r6[1] * cos(ang/beta)
            drift6!(r6, cavity.len / 2.0, beti)
        end
    end
    return nothing
end

function CRABCAVITYK2Pass_P!(r::Matrix{Float64}, cavity::CRABCAVITY_K2, beta::Float64, E::Float64, num_particles::Int64, lost_flags::Array{Int,1})
    if use_exact_beti == 1
        beti = 1.0 / beta
    else
        beti = 1.0 
    end
    Threads.@threads for c in 1:num_particles
        if isone(lost_flags[c])
            continue
        end
        r6 = @view r[c, :] 
        ang = cavity.k * r6[5] + cavity.phi
        if cavity.len == 0.0
            # r6[2] += (cavity.volt/beta2E) * sin(ang)
            # r6[6] += (-cavity.k * cavity.volt/beta2E) * r6[1] * cos(ang)
            # this is symplectic form
            r6[2] += (cavity.volt/E) * sin(ang/beta)
            r6[6] += (-cavity.k * cavity.volt/E/beta) * r6[1] * cos(ang/beta)
            # sextupole kick
            r6[2] -= cavity.k2 * (r6[1]^2 - r6[3]^2) * sin(ang)
            r6[4] += 2.0 * cavity.k2 * r6[1] * r6[3] * sin(ang)
            r6[6] -= (cavity.k2 * cavity.k / 3.0) * (r6[1]^3 - 3.0 * r6[1] * r6[3]^2) * cos(ang)
        else
            drift6!(r6, cavity.len / 2.0, beti)
            # r6[2] += (cavity.volt/beta2E) * sin(ang)
            # r6[6] += (-cavity.k * cavity.volt/beta2E) * r6[1] * cos(ang)
            # this is symplectic form
            r6[2] += (cavity.volt/E) * sin(ang/beta)
            r6[6] += (-cavity.k * cavity.volt/E/beta) * r6[1] * cos(ang/beta)
            # sextupole kick
            r6[2] -= cavity.k2 * (r6[1]^2 - r6[3]^2) * sin(ang)
            r6[4] += 2.0 * cavity.k2 * r6[1] * r6[3] * sin(ang)
            r6[6] -= (cavity.k2 * cavity.k / 3.0) * (r6[1]^3 - 3.0 * r6[1] * r6[3]^2) * cos(ang)
            drift6!(r6, cavity.len / 2.0, beti)
        end
    end
    return nothing
end

function easyCRABCAVITYPass_P!(r::Matrix{Float64}, cavity::easyCRABCAVITY, num_particles::Int64, lost_flags::Array{Int,1})
    Threads.@threads for c in 1:num_particles
    # for c in 1:num_particles
        if isone(lost_flags[c])
            continue
        end
        r6 = @view r[c, :] 
        ang = cavity.k * r6[5] + cavity.phi
        r6[1] += (cavity.halfthetac/cavity.k) * sin(ang)
        r6[6] += cavity.halfthetac * r6[2] * cos(ang)
    end
    return nothing
end

function AccelCavityPass_P!(r::Matrix{Float64}, cavity::AccelCavity, beta2E::Float64, num_particles::Int64, lost_flags::Array{Int,1})
    v_beta2E = cavity.volt/beta2E
    Threads.@threads for c in 1:num_particles
    # for c in 1:num_particles
        if lost_flags[c] == 1
            continue
        end
        r6 = @view r[c, :] 
        sv = sin(cavity.k * r6[5] + cavity.phis) - sin(cavity.phis)
        r6[6] += v_beta2E * sv
    end
    return nothing
end

function pass_P!(ele::CRABCAVITY, r_in::Matrix{Float64}, num_particles::Int64, particles::Beam{Float64})
    lost_flags = particles.lost_flag
    CRABCAVITYPass_P!(r_in, ele, particles.beta, particles.energy, num_particles,lost_flags)
    return nothing
end

function pass_P!(ele::CRABCAVITY_K2, r_in::Matrix{Float64}, num_particles::Int64, particles::Beam{Float64})
    lost_flags = particles.lost_flag
    CRABCAVITYK2Pass_P!(r_in, ele, particles.beta, particles.energy, num_particles,lost_flags)
    return nothing
end

function pass_P!(ele::easyCRABCAVITY, r_in::Matrix{Float64}, num_particles::Int64, particles::Beam{Float64})
    lost_flags = particles.lost_flag
    easyCRABCAVITYPass_P!(r_in, ele, num_particles,lost_flags)
    return nothing
end

function pass_P!(ele::AccelCavity, r_in::Matrix{Float64}, num_particles::Int64, particles::Beam{Float64})
    lost_flags = particles.lost_flag
    beta2E = particles.beta^2 * particles.energy
    AccelCavityPass_P!(r_in, ele, beta2E, num_particles,lost_flags)
    return nothing
end


##########################################################################################
# High-order TPSA
function pass_TPSA!(ele::CRABCAVITY, r_in::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}; E0::Float64=0.0, m0::Float64=m_e) where {T, TPS_Dim, Max_TPS_Degree}
    if E0 == 0.0
        println("Warning: beam energy is not defined")
    end
    gamma = (E0 + m0) / m0
    beta = sqrt(1.0 - 1.0 / gamma^2)
    if use_exact_beti == 1
        beti = 1.0 / beta
    else
        beti = 1.0 
    end
    beta2E = beta^2 * E0
    ang = ele.k * r_in[5] + ele.phi
    if ele.len == 0.0
        # r_in[2] += ele.volt/beta2E * sin(ang)
        # r_in[6] += (-ele.k * ele.volt/beta2E) * r_in[1] * cos(ang)
        # this is symplectic form
        r_in[2] += (ele.volt/E0) * sin(ang/beta)
        r_in[6] += (-ele.k * ele.volt/E0/beta) * r_in[1] * cos(ang/beta)
        return nothing
    else
        drift6!(r_in, ele.len / 2.0, beti)
        # r_in[2] += ele.volt/beta2E * sin(ang)
        # r_in[6] += (-ele.k * ele.volt/beta2E) * r_in[1] * cos(ang)
        # this is symplectic form
        r_in[2] += (ele.volt/E0) * sin(ang/beta)
        r_in[6] += (-ele.k * ele.volt/E0/beta) * r_in[1] * cos(ang/beta)
        # println("crabcavity is not implemented in TPSA")
        drift6!(r_in, ele.len / 2.0, beti)
    end
    return nothing
end

function pass_TPSA!(ele::CRABCAVITY_K2, r_in::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}; E0::Float64=0.0, m0::Float64=m_e) where {T, TPS_Dim, Max_TPS_Degree}
    if E0 == 0.0
        println("Warning: beam energy is not defined")
    end
    gamma = (E0 + m0) / m0
    beta = sqrt(1.0 - 1.0 / gamma^2)
    ang = ele.k * r_in[5] + ele.phi
    if use_exact_beti == 1
        beti = 1.0 / beta
    else
        beti = 1.0 
    end
    if ele.len == 0.0
        # r_in[2] += ele.volt/beta2E * sin(ang)
        # r_in[6] += (-ele.k * ele.volt/beta2E) * r_in[1] * cos(ang)
        # this is symplectic form
        r_in[2] += (ele.volt/E0) * sin(ang/beta)
        r_in[6] += (-ele.k * ele.volt/E0/beta) * r_in[1] * cos(ang/beta)
        # sextupole kick
        r_in[2] -= ele.k2 * (r_in[1]^2 - r_in[3]^2) * sin(ang)
        r_in[4] += 2.0 * ele.k2 * r_in[1] * r_in[3] * sin(ang)
        r_in[6] -= (-ele.k2 * ele.k / 3.0) * (r_in[1]^3 - 3.0 * r_in[1] * r_in[3]^2) * cos(ang)
        return nothing
    else
        drift6!(r_in, ele.len / 2.0, beti)
        # r_in[2] += ele.volt/beta2E * sin(ang)
        # r_in[6] += (-ele.k * ele.volt/beta2E) * r_in[1] * cos(ang)
        # this is symplectic form
        r_in[2] += (ele.volt/E0) * sin(ang/beta)
        r_in[6] += (-ele.k * ele.volt/E0/beta) * r_in[1] * cos(ang/beta)
        # sextupole kick
        r_in[2] -= ele.k2 * (r_in[1]^2 - r_in[3]^2) * sin(ang)
        r_in[4] += 2.0 * ele.k2 * r_in[1] * r_in[3] * sin(ang)
        r_in[6] -= (-ele.k2 * ele.k / 3.0) * (r_in[1]^3 - 3.0 * r_in[1] * r_in[3]^2) * cos(ang)
        drift6!(r_in, ele.len / 2.0, beti)
    end
    return nothing
end

function pass_TPSA!(ele::easyCRABCAVITY, r_in::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}; E0::Float64=0.0, m0::Float64=m_e) where {T, TPS_Dim, Max_TPS_Degree}
    gamma = (E0 + m0) / m0
    beta = sqrt(1.0 - 1.0 / gamma^2)
    if use_exact_beti == 1
        beti = 1.0 / beta
    else
        beti = 1.0 
    end
    if ele.len == 0.0
        return nothing
    else
        drift6!(r_in, ele.len, beti)
    end
    println("warning: easyCRABCAVITY is not implemented in TPSA")
    return nothing
end

function pass_TPSA!(ele::AccelCavity, r_in::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}; E0::Float64=0.0, m0::Float64=m_e) where {T, TPS_Dim, Max_TPS_Degree}
    gamma = (E0 + m0) / m0
    beta = sqrt(1.0 - 1.0 / gamma^2)
    if use_exact_beti == 1
        beti = 1.0 / beta
    else
        beti = 1.0 
    end
    if ele.len == 0.0
        return nothing
    else
        drift6!(r_in, ele.len, beti)
    end
    println("warning: AccelCavity is not implemented in TPSA")
    return nothing
end