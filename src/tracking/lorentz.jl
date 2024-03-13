function pass!(ele::LorentzBoost, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam)
    lost_flags = particles.lost_flag
    if ele.mode == 0
        invcosang = 1.0 / ele.cosang
        for c in 1:num_particles
            if lost_flags[c] == 1
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

function pass!(ele::InvLorentzBoost, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam)
    lost_flags = particles.lost_flag
    if ele.mode == 0
        invcosang = 1.0 / ele.cosang
        for c in 1:num_particles
            if lost_flags[c] == 1
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
            if lost_flags[c] == 1
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
            if lost_flags[c] == 1
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