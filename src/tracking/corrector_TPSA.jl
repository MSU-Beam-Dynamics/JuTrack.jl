function CorrectorPass_TPSA!(r::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}, le::Float64, xkick::Float64, ykick::Float64,
    T1::Array{Float64,1}, T2::Array{Float64,1}, R1::Array{Float64,2}, R2::Array{Float64,2}) where {T, TPS_Dim, Max_TPS_Degree}

    # Threads.@threads for c in 1:num_particles
    # for c in 1:num_particles
    #     if lost_flags[c] == 1
    #         continue
    #     end
        # r6 = @view r[(c-1)*6+1:c*6]
        # if !isnan(r6[1])
            # p_norm = 1.0 / sqrt((1.0 + r6[6])^2 - r6[2]^2 - r6[4]^2)
            if false
                p_norm = 1.0 / sqrt((1.0 + r[6])^2 - r[2]^2 - r[4]^2)
            else
                p_norm = 1.0 / sqrt(1.0 + r[6])
            end
            NormL = le * p_norm
            # Misalignment at entrance
            if T1 != zeros(6)
                ATaddvv!(r, T1)
            end
            if R1 != zeros(6, 6)
                ATmultmv!(r, R1)
            end
            
            r[5] += NormL*p_norm*(xkick*xkick/3.0 + ykick*ykick/3.0 +
   		            r[2]*r[2] + r[4]*r[4] +
   		            r[2]*xkick + r[4]*ykick)/2.0
            r[1] += NormL*(r[2]+xkick/2.0)
		    r[2] += xkick
		    r[3] += NormL*(r[4]+ykick/2.0)
   		    r[4] += ykick

            # Misalignment at exit
            if R2 != zeros(6, 6)
                ATmultmv!(r, R2)
            end
            if T2 != zeros(6)
                ATaddvv!(r, T2)
            end

        # end
    # end

    return nothing
end

function pass_TPSA!(ele::CORRECTOR, r_in::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}) where {T, TPS_Dim, Max_TPS_Degree}
    # ele: CORRECTOR
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    CorrectorPass_TPSA!(r_in, ele.len, ele.xkick, ele.ykick, ele.T1, ele.T2, ele.R1, ele.R2,)
    return nothing
end

