function CorrectorPass!(r::Array{Float64,1}, le::Float64, xkick::Float64, ykick::Float64,
    T1::Array{Float64,1}, T2::Array{Float64,1}, R1::Array{Float64,2}, R2::Array{Float64,2}, 
    num_particles::Int, lost_flags::Array{Int64,1})

    # Threads.@threads for c in 1:num_particles
    for c in 1:num_particles
        if lost_flags[c] == 1
            continue
        end
        r6 = @view r[(c-1)*6+1:c*6]
        if !isnan(r6[1])
            if false
                p_norm = 1.0 / sqrt((1.0 + r6[6])^2 - r6[2]^2 - r6[4]^2)
            else
                p_norm = 1.0 / sqrt(1.0 + r6[6])
            end
            NormL = le * p_norm
            # Misalignment at entrance
            if T1 != zeros(6)
                ATaddvv!(r6, T1)
            end
            if R1 != zeros(6, 6)
                ATmultmv!(r6, R1)
            end

            # Check physical apertures at the entrance of the magnet
            # if RApertures != nothing
            #     checkiflostRectangularAp(r6, RApertures)
            # end
            # if EApertures != nothing
            #     checkiflostEllipticalAp(r6, EApertures)
            # end
            
            r6[5] += NormL*p_norm*(xkick*xkick/3 + ykick*ykick/3 +
   		            r6[2]*r6[2] + r6[4]*r6[4] +
   		            r6[2]*xkick + r6[4]*ykick)/2
            r6[1] += NormL*(r6[2]+xkick/2)
		    r6[2] += xkick
		    r6[3] += NormL*(r6[4]+ykick/2)
   		    r6[4] += ykick
            # Check physical apertures at the exit of the magnet
            # if RApertures != nothing
            #     checkiflostRectangularAp(r6, RApertures)
            # end
            # if EApertures != nothing
            #     checkiflostEllipticalAp(r6, EApertures)
            # end

            # Misalignment at exit
            if R2 != zeros(6, 6)
                ATmultmv!(r6, R2)
            end
            if T2 != zeros(6)
                ATaddvv!(r6, T2)
            end
            if r6[1] > CoordLimit || r6[2] > AngleLimit || r6[1] < -CoordLimit || r6[2] < -AngleLimit
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
            if false
                p_norm = 1.0 / sqrt((1.0 + r6[6])^2 - r6[2]^2 - r6[4]^2)
            else
                p_norm = 1.0 / sqrt(1.0 + r6[6])
            end
            NormL = le * p_norm
            # Misalignment at entrance
            if T1 != zeros(6)
                ATaddvv!(r6, T1)
            end
            if R1 != zeros(6, 6)
                ATmultmv!(r6, R1)
            end

            # Check physical apertures at the entrance of the magnet
            # if RApertures != nothing
            #     checkiflostRectangularAp(r6, RApertures)
            # end
            # if EApertures != nothing
            #     checkiflostEllipticalAp(r6, EApertures)
            # end
            
            r6[5] += NormL*p_norm*(xkick*xkick/3 + ykick*ykick/3 +
   		            r6[2]*r6[2] + r6[4]*r6[4] +
   		            r6[2]*xkick + r6[4]*ykick)/2
            r6[1] += NormL*(r6[2]+xkick/2)
		    r6[2] += xkick
		    r6[3] += NormL*(r6[4]+ykick/2)
   		    r6[4] += ykick
            # Check physical apertures at the exit of the magnet
            # if RApertures != nothing
            #     checkiflostRectangularAp(r6, RApertures)
            # end
            # if EApertures != nothing
            #     checkiflostEllipticalAp(r6, EApertures)
            # end

            # Misalignment at exit
            if R2 != zeros(6, 6)
                ATmultmv!(r6, R2)
            end
            if T2 != zeros(6)
                ATaddvv!(r6, T2)
            end
            if r6[1] > CoordLimit || r6[2] > AngleLimit || r6[1] < -CoordLimit || r6[2] < -AngleLimit
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