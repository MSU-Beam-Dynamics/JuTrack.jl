function pass!(ele::SOLENOID, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam{Float64})
    lost_flags = particles.lost_flag
    T1 = ele.T1
    R1 = ele.R1
    T2 = ele.T2
    R2 = ele.R2
    if use_exact_beti == 1
        beti = 1.0 / particles.beta
    else
        beti = 1.0 
    end
    # Threads.@threads for c in 1:num_particles
    if ele.ks != 0.0
        for c in 1:num_particles
            if lost_flags[c] == 1
                continue
            end
            r6 = @view r_in[(c-1)*6+1:c*6]
            if !isnan(r6[1]) 
                p_norm = 1.0 / (1.0 + r6[6]) # use linearized p_norm
    
                # Misalignment at entrance
                if !iszero(T1)
                    addvv!(r6, T1)
                end
                if !iszero(R1)
                    multmv!(r6, R1)
                end
    
                x = r6[1]
                xpr = r6[2]*p_norm
                y = r6[3]
                ypr = r6[4]*p_norm
                H = ele.ks * p_norm / 2.0
                S = sin(ele.len * H)
                C = cos(ele.len * H)
                r6[1] = x*C*C + xpr*C*S/H + y*C*S + ypr*S*S/H
                r6[2] = (-x*H*C*S + xpr*C*C - y*H*S*S + ypr*C*S) / p_norm
                r6[3] = -x*C*S - xpr*S*S/H + y*C*C + ypr*C*S/H
                r6[4] = (x*H*S*S - xpr*C*S - y*C*S*H + ypr*C*C) / p_norm
                r6[5] += ele.len*(H*H*(x*x+y*y) + 2.0*H*(xpr*y-ypr*x) +xpr*xpr+ypr*ypr)/2.0
    
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
    else
        for c in 1:num_particles
            if lost_flags[c] == 1
                continue
            end
            r6 = @view r_in[(c-1)*6+1:c*6]
            drift6!(r6, ele.len, beti)
            if check_lost(r6)
                lost_flags[c] = 1
            end
        end
    end

    return nothing

end

##########################################################################################
# multi-threading
function pass_P!(ele::SOLENOID, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam{Float64})
    lost_flags = particles.lost_flag
    T1 = ele.T1
    R1 = ele.R1
    T2 = ele.T2
    R2 = ele.R2
    beti = 1.0 / particles.beta
    if ele.ks != 0.0
        Threads.@threads for c in 1:num_particles
            if lost_flags[c] == 1
                continue
            end
            r6 = @view r_in[(c-1)*6+1:c*6]
            if !isnan(r6[1]) 
                p_norm = 1.0 / (1.0 + r6[6])
    
                # Misalignment at entrance
                if !iszero(T1)
                    addvv!(r6, T1)
                end
                if !iszero(R1)
                    multmv!(r6, R1)
                end
    
                x = r6[1]
                xpr = r6[2]*p_norm
                y = r6[3]
                ypr = r6[4]*p_norm
                H = ele.ks * p_norm / 2.0
                S = sin(ele.len * H)
                C = cos(ele.len * H)
                r6[1] = x*C*C + xpr*C*S/H + y*C*S + ypr*S*S/H
                r6[2] = (-x*H*C*S + xpr*C*C - y*H*S*S + ypr*C*S) / p_norm
                r6[3] = -x*C*S - xpr*S*S/H + y*C*C + ypr*C*S/H
                r6[4] = (x*H*S*S - xpr*C*S - y*C*S*H + ypr*C*C) / p_norm
                r6[5] += ele.len*(H*H*(x*x+y*y) + 2.0*H*(xpr*y-ypr*x) +xpr*xpr+ypr*ypr)/2.0
    
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
    else
        Threads.@threads for c in 1:num_particles
            if lost_flags[c] == 1
                continue
            end
            r6 = @view r_in[(c-1)*6+1:c*6]
            if !iszero(T1)
                addvv!(r_in, T1)
            end
            if !iszero(R1)
                multmv!(r_in, R1)
            end
            drift6!(r6, ele.len, beti)
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