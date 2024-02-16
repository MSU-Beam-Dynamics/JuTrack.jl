function CrabCavityPass_P!(r::Array{Float64,1}, cavity::CrabCavity, beta2E::Float64, num_particles::Int64, lost_flags::Array{Int,1})
    Threads.@threads for c in 1:num_particles
    # for c in 1:num_particles
        if lost_flags[c] == 1
            continue
        end
        r6 = @view r[(c-1)*6+1:c*6] 
        ang = -cavity.k * r6[5] + cavity.phi
        r6[2] += cavity.volt/beta2E * sin(ang)
        r6[6] += (-cavity.k * cavity.volt/beta2E) * r[1] * cos(ang)
    end
    return nothing
end

function easyCrabCavityPass_P!(r::Array{Float64,1}, cavity::easyCrabCavity, num_particles::Int64, lost_flags::Array{Int,1})
    Threads.@threads for c in 1:num_particles
    # for c in 1:num_particles
        if lost_flags[c] == 1
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

function pass_P!(ele::CrabCavity, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam, noTarray::Array{Float64,1}, noRmatrix::Array{Float64,2})
    beta2E = particles.beta^2 * particles.energy
    lost_flags = particles.lost_flag
    CrabCavityPass_P!(r_in, ele, beta2E, num_particles,lost_flags)
    return nothing
end

function pass_P!(ele::easyCrabCavity, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam, noTarray::Array{Float64,1}, noRmatrix::Array{Float64,2})
    lost_flags = particles.lost_flag
    easyCrabCavityPass_P!(r_in, ele, num_particles,lost_flags)
    return nothing
end

function pass_P!(ele::AccelCavity, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam, noTarray::Array{Float64,1}, noRmatrix::Array{Float64,2})
    lost_flags = particles.lost_flag
    beta2E = particles.beta^2 * particles.energy
    AccelCavityPass_P!(r_in, ele, beta2E, num_particles,lost_flags)
    return nothing
end