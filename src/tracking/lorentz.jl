"""
    function pass!(ele::LorentzBoost, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam)

This is a function to apply linearized Lorentz Boost at IP for converting from lab frame so that the beam is colliding with the opposing beam head-on.

# Arguments
- ele::LorentzBoost: a Lorentz boost element
- r_in::Array{Float64,1}: 6-by-num_particles array
- num_particles::Int64: number of particles
- particles::Beam: beam object
"""
function pass!(ele::LorentzBoost, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam)
    lost_flags = particles.lost_flag
    if ele.mode == 0
        invcosang = 1.0 / ele.cosang
        for c in 1:num_particles
            if isone(lost_flags[c])
                continue
            end
            r6 = @view r_in[(c-1)*6+1:c*6] 
            r6[1] += ele.tanang * r6[5]
            r6[6] -= ele.tanang * r6[2]
            r6[2] *= invcosang
            r6[4] *= invcosang
            r6[5] *= invcosang
        end
    end
    return nothing
end


"""
    function pass!(ele::InvLorentzBoost, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam)

This is a function to apply linearized inverse Lorentz Boost at IP for converting from boosted frame to lab frame

# Arguments
- ele::InvLorentzBoost: an inverse Lorentz boost element
- r_in::Array{Float64,1}: 6-by-num_particles array
- num_particles::Int64: number of particles
- particles::Beam: beam object
"""
function pass!(ele::InvLorentzBoost, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam)
    lost_flags = particles.lost_flag
    if ele.mode == 0
        invcosang = 1.0 / ele.cosang
        for c in 1:num_particles
            if isone(lost_flags[c])
                continue
            end
            r6 = @view r_in[(c-1)*6+1:c*6] 
            r6[1] -= ele.sinang * r6[5]
            r6[6] += ele.sinang * r6[2]
            r6[2] *= ele.cosang
            r6[4] *= ele.cosang
            r6[5] *= ele.cosang
        end
    end
    return nothing
end



function pass_P!(ele::LorentzBoost, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam)
    lost_flags = particles.lost_flag
    if ele.mode == 0
        invcosang = 1.0 / ele.cosang
        Threads.@threads for c in 1:num_particles
        # for c in 1:num_particles
            if isone(lost_flags[c])
                continue
            end
            r6 = @view r_in[(c-1)*6+1:c*6] 
            r6[1] += ele.tanang * r6[5]
            r6[6] -= ele.tanang * r6[2]
            r6[2] *= invcosang
            r6[4] *= invcosang
            r6[5] *= invcosang
        end
    end
    return nothing
end

function pass_P!(ele::InvLorentzBoost, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam)
    lost_flags = particles.lost_flag
    if ele.mode == 0
        invcosang = 1.0 / ele.cosang
        Threads.@threads for c in 1:num_particles
        # for c in 1:num_particles
            if isone(lost_flags[c])
                continue
            end
            r6 = @view r_in[(c-1)*6+1:c*6] 
            r6[1] -= ele.sinang * r6[5]
            r6[6] += ele.sinang * r6[2]
            r6[2] *= ele.cosang
            r6[4] *= ele.cosang
            r6[5] *= ele.cosang
        end
    end
    return nothing
end