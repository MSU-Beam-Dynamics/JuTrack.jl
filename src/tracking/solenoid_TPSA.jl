include("drift_TPSA.jl")
function pass_TPSA!(ele::SOLENOID, r_in::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}) where {T, TPS_Dim, Max_TPS_Degree}
    T1 = ele.T1
    R1 = ele.R1
    T2 = ele.T2
    R2 = ele.R2
    # Threads.@threads for c in 1:num_particles
    if ele.ks != 0.0
        # for c in 1:num_particles
        #     if lost_flags[c] == 1
        #         continue
        #     end
            # r6 = @view r_in[(c-1)*6+1:c*6]
            # if !isnan(r6[1]) 
            if use_exact_Hamiltonian == 1
                p_norm = 1.0 / sqrt((1.0 + r_in[6])^2 - r_in[2]^2 - r_in[4]^2) # use exact p_norm
            else
                p_norm = 1.0 / sqrt(1.0 + r_in[6]) # use linearized p_norm
            end
                # Misalignment at entrance
                if T1 != zeros(6)
                    ATaddvv!(r_in, T1)
                end
                if R1 != zeros(6, 6)
                    ATmultmv!(r_in, R1)
                end
    
                x = r_in[1]
                xpr = r_in[2]*p_norm
                y = r_in[3]
                ypr = r_in[4]*p_norm
                H = ele.ks * p_norm / 2.0
                S = sin(ele.len * H)
                C = cos(ele.len * H)
                r_in[1] = x*C*C + xpr*C*S/H + y*C*S + ypr*S*S/H
                r_in[2] = (-x*H*C*S + xpr*C*C - y*H*S*S + ypr*C*S) / p_norm
                r_in[3] = -x*C*S - xpr*S*S/H + y*C*C + ypr*C*S/H
                r_in[4] = (x*H*S*S - xpr*C*S - y*C*S*H + ypr*C*C) / p_norm
                r_in[5] += ele.len*(H*H*(x*x+y*y) + 2.0*H*(xpr*y-ypr*x) +xpr*xpr+ypr*ypr)/2.0
    
                if R2 != zeros(6, 6)
                    ATmultmv!(r_in, R2)
                end
                if T2 != zeros(6) 
                    ATaddvv!(r_in, T2)
                end
            # end
        # end
    else
        # for c in 1:num_particles
        #     if lost_flags[c] == 1
        #         continue
        #     end
        #     r6 = @view r_in[(c-1)*6+1:c*6]
            drift6!(r_in, ele.len)
        # end
    end

    return nothing

end