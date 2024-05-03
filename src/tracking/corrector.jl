function CorrectorPass!(r::Array{Float64,1}, le::Float64, xkick::Float64, ykick::Float64,
    T1::Array{Float64,1}, T2::Array{Float64,1}, R1::Array{Float64,2}, R2::Array{Float64,2}, 
    num_particles::Int, lost_flags::Array{Int64,1})
    # Modified based on AT function. Ref[Terebilo, Andrei. "Accelerator modeling with MATLAB accelerator toolbox." PACS2001 (2001)].

    for c in 1:num_particles
        if lost_flags[c] == 1
            continue
        end
        r6 = @view r[(c-1)*6+1:c*6]
        if !isnan(r6[1])
            p_norm = 1.0 / (1.0 + r6[6])
            
            NormL = le * p_norm
            # Misalignment at entrance
            if !iszero(T1)
                addvv!(r6, T1)
            end
            if !iszero(R1)
                multmv!(r6, R1)
            end

            
            r6[5] += NormL*p_norm*(xkick*xkick/3 + ykick*ykick/3 +
   		            r6[2]*r6[2] + r6[4]*r6[4] +
   		            r6[2]*xkick + r6[4]*ykick)/2
            r6[1] += NormL*(r6[2]+xkick/2)
		    r6[2] += xkick
		    r6[3] += NormL*(r6[4]+ykick/2)
   		    r6[4] += ykick


            # Misalignment at exit
            if !iszero(R2)
                multmv!(r6, R2)
            end
            if !iszero(T2)
                addvv!(r6, T2)
            end
            if abs(r6[1]) > CoordLimit || abs(r6[2]) > AngleLimit || isnan(r6[1]) || isinf(r6[1])
                lost_flags[c] = 1
            end
        end
    end

    return nothing
end

function pass!(ele::CORRECTOR, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam)
    # ele: CORRECTOR
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    lost_flags = particles.lost_flag
    CorrectorPass!(r_in, ele.len, ele.xkick, ele.ykick, ele.T1, ele.T2, ele.R1, ele.R2, num_particles, lost_flags)
    return nothing
end

function CorrectorPass_P!(r::Array{Float64,1}, le::Float64, xkick::Float64, ykick::Float64,
    T1::Array{Float64,1}, T2::Array{Float64,1}, R1::Array{Float64,2}, R2::Array{Float64,2}, 
    num_particles::Int, lost_flags::Array{Int64,1})

    Threads.@threads for c in 1:num_particles
    # for c in 1:num_particles
        if lost_flags[c] == 1
            continue
        end
        r6 = @view r[(c-1)*6+1:c*6]
        if !isnan(r6[1])
            if 1.0 + r6[6] >= 0
                p_norm = 1.0 / sqrt(1.0 + r6[6])
            else
                lost_flags[c] = 1
                continue
            end
            NormL = le * p_norm
            # Misalignment at entrance
            if !iszero(T1)
                addvv!(r6, T1)
            end
            if !iszero(R1)
                multmv!(r6, R1)
            end


            r6[5] += NormL*p_norm*(xkick*xkick/3 + ykick*ykick/3 +
   		            r6[2]*r6[2] + r6[4]*r6[4] +
   		            r6[2]*xkick + r6[4]*ykick)/2
            r6[1] += NormL*(r6[2]+xkick/2)
		    r6[2] += xkick
		    r6[3] += NormL*(r6[4]+ykick/2)
   		    r6[4] += ykick

            # Misalignment at exit
            if !iszero(R2)
                multmv!(r6, R2)
            end
            if !iszero(T2)
                addvv!(r6, T2)
            end
            if abs(r6[1]) > CoordLimit || abs(r6[2]) > AngleLimit || isnan(r6[1]) || isinf(r6[1])
                lost_flags[c] = 1
            end
        end
    end

    return nothing
end

function pass_P!(ele::CORRECTOR, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam)
    # ele: CORRECTOR
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    lost_flags = particles.lost_flag
    CorrectorPass_P!(r_in, ele.len, ele.xkick, ele.ykick, ele.T1, ele.T2, ele.R1, ele.R2, num_particles, lost_flags)
    return nothing
end