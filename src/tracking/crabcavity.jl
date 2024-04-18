function CRABCAVITYPass!(r::Array{Float64,1}, cavity::CRABCAVITY, beta2E::Float64, num_particles::Int64, lost_flags::Array{Int,1})
    for c in 1:num_particles
        if isone(lost_flags[c])
            continue
        end
        r6 = @view r[(c-1)*6+1:c*6] 
        ang = -cavity.k * r6[5] + cavity.phi
        if cavity.len == 0.0
            r6[2] += (cavity.volt/beta2E) * sin(ang)
            r6[6] += (-cavity.k * cavity.volt/beta2E) * r6[1] * cos(ang)
        else
            drift6!(r6, cavity.len / 2.0)
            r6[2] += (cavity.volt/beta2E) * sin(ang)
            r6[6] += (-cavity.k * cavity.volt/beta2E) * r6[1] * cos(ang)
            drift6!(r6, cavity.len / 2.0)
        end
    end
    return nothing
end

function easyCRABCAVITYPass!(r::Array{Float64,1}, cavity::easyCRABCAVITY, num_particles::Int64, lost_flags::Array{Int,1})
    for c in 1:num_particles
        if isone(lost_flags[c])
            continue
        end
        r6 = @view r[(c-1)*6+1:c*6] 
        ang = -cavity.k * r6[5] + cavity.phi
        r6[1] += (cavity.halfthetac/cavity.k) * sin(ang)
        r6[6] += cavity.halfthetac * r6[2] * cos(ang)
    end
    return nothing
end

function AccelCavityPass!(r::Array{Float64,1}, cavity::AccelCavity, beta2E::Float64, num_particles::Int64, lost_flags::Array{Int,1})
    v_beta2E = cavity.volt/beta2E
    # Threads.@threads for c in 1:num_particles
    for c in 1:num_particles
        if lost_flags[c] == 1
            continue
        end
        r6 = @view r[(c-1)*6+1:c*6] 
        sv = sin(-cavity.k * r6[5] + cavity.phis) - sin(cavity.phis)
        r6[6] += v_beta2E * sv
    end
    return nothing
end

function pass!(ele::CRABCAVITY, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam)
    beta2E = particles.beta^2 * particles.energy
    lost_flags = particles.lost_flag
    CRABCAVITYPass!(r_in, ele, beta2E, num_particles,lost_flags)
    return nothing
end

function pass!(ele::easyCRABCAVITY, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam)
    lost_flags = particles.lost_flag
    easyCRABCAVITYPass!(r_in, ele, num_particles,lost_flags)
    return nothing
end

function pass!(ele::AccelCavity, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam)
    lost_flags = particles.lost_flag
    beta2E = particles.beta^2 * particles.energy
    AccelCavityPass!(r_in, ele, beta2E, num_particles,lost_flags)
    return nothing
end

##########################################################################################
# multi-threading
function CRABCAVITYPass_P!(r::Array{Float64,1}, cavity::CRABCAVITY, beta2E::Float64, num_particles::Int64, lost_flags::Array{Int,1})
    Threads.@threads for c in 1:num_particles
    # for c in 1:num_particles
        if isone(lost_flags[c])
            continue
        end
        r6 = @view r[(c-1)*6+1:c*6] 
        ang = -cavity.k * r6[5] + cavity.phi
        if cavity.len == 0.0
            r6[2] += (cavity.volt/beta2E) * sin(ang)
            r6[6] += (-cavity.k * cavity.volt/beta2E) * r6[1] * cos(ang)
        else
            drift6!(r6, cavity.len / 2.0)
            r6[2] += (cavity.volt/beta2E) * sin(ang)
            r6[6] += (-cavity.k * cavity.volt/beta2E) * r6[1] * cos(ang)
            drift6!(r6, cavity.len / 2.0)
        end
    end
    return nothing
end

function easyCRABCAVITYPass_P!(r::Array{Float64,1}, cavity::easyCRABCAVITY, num_particles::Int64, lost_flags::Array{Int,1})
    Threads.@threads for c in 1:num_particles
    # for c in 1:num_particles
        if isone(lost_flags[c])
            continue
        end
        r6 = @view r[(c-1)*6+1:c*6] 
        ang = -cavity.k * r6[5] + cavity.phi
        r6[1] += (cavity.halfthetac/cavity.k) * sin(ang)
        r6[6] += cavity.halfthetac * r6[2] * cos(ang)
    end
    return nothing
end

function AccelCavityPass_P!(r::Array{Float64,1}, cavity::AccelCavity, beta2E::Float64, num_particles::Int64, lost_flags::Array{Int,1})
    v_beta2E = cavity.volt/beta2E
    Threads.@threads for c in 1:num_particles
    # for c in 1:num_particles
        if lost_flags[c] == 1
            continue
        end
        r6 = @view r[(c-1)*6+1:c*6] 
        sv = sin(-cavity.k * r6[5] + cavity.phis) - sin(cavity.phis)
        r6[6] += v_beta2E * sv
    end
    return nothing
end

function pass_P!(ele::CRABCAVITY, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam)
    beta2E = particles.beta^2 * particles.energy
    lost_flags = particles.lost_flag
    CRABCAVITYPass_P!(r_in, ele, beta2E, num_particles,lost_flags)
    return nothing
end

function pass_P!(ele::easyCRABCAVITY, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam)
    lost_flags = particles.lost_flag
    easyCRABCAVITYPass_P!(r_in, ele, num_particles,lost_flags)
    return nothing
end

function pass_P!(ele::AccelCavity, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam)
    lost_flags = particles.lost_flag
    beta2E = particles.beta^2 * particles.energy
    AccelCavityPass_P!(r_in, ele, beta2E, num_particles,lost_flags)
    return nothing
end