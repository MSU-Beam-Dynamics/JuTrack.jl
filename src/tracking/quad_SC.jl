function QuadLinearPass_SC!(r::Array{Float64,1}, le::Float64, k1::Float64, 
    T1::Array{Float64,1}, T2::Array{Float64,1}, R1::Array{Float64,2}, R2::Array{Float64,2}, 
    RApertures::Array{Float64,1}, EApertures::Array{Float64,1},
    num_particles::Int, lost_flags::Array{Int64,1}, a, b, Nl, Nm, K)
    # This is linear quadrupole pass function. Not used for symplectic tracking.

    # half-kick-half
    for c in 1:num_particles
        if lost_flags[c] == 1
            continue
        end
        r6 = @view r[(c-1)*6+1:c*6]
        if !isnan(r6[1])
            p_norm = 1.0 / (1.0 + r6[6])

            if iszero(k1)
                drift6!(r6, le/2.0)
            else
                # Misalignment at entrance
                if !iszero(T1)
                    addvv!(r6, T1)
                end
                if !iszero(R1)
                    multmv!(r6, R1)
                end

                x   = r6[1]
                xpr = r6[2]*p_norm
                y   = r6[3]
                ypr = r6[4]*p_norm
                
                g  = abs(k1)/(1.0 + r6[6])
                t  = sqrt(g)
                lt = le/2.0 * t
                
                if k1>0
                    MHD = cos(lt)
                    M12 = sin(lt)/t
                    M21 = -M12*g
                    MVD = cosh(lt)
                    M34 = sinh(lt)/t
                    M43 = M34*g
                else
                    MHD = cosh(lt)
                    M12 = sinh(lt)/t
                    M21 = M12*g
                    MVD = cos(lt)
                    M34 = sin(lt)/t
                    M43 = -M34*g
                end

                r6[1] =  MHD*x + M12*xpr
                r6[2] = (M21*x + MHD*xpr)/p_norm
                r6[3] =  MVD*y + M34*ypr
                r6[4] = (M43*y + MVD*ypr)/p_norm

                r6[5]+= g*(x*x*(le/2.0-MHD*M12)-y*y*(le/2.0-MVD*M34))/4.0
                r6[5]+= (xpr*xpr*(le/2.0+MHD*M12)+ypr*ypr*(le/2.0+MVD*M34))/4.0
                r6[5]+= (x*xpr*M12*M21 + y*ypr*M34*M43)/2.0

                # Misalignment at exit
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
    end

    space_charge!(r, K, Nl, Nm, a/Nl, b/Nm, a, b, num_particles, le, lost_flags)

    for c in 1:num_particles
        if lost_flags[c] == 1
            continue
        end
        r6 = @view r[(c-1)*6+1:c*6]
        if !isnan(r6[1])
            p_norm = 1.0 / (1.0 + r6[6])

            if iszero(k1)
                drift6!(r6, le/2.0)
            else
                # Misalignment at entrance
                if !iszero(T1)
                    addvv!(r6, T1)
                end
                if !iszero(R1)
                    multmv!(r6, R1)
                end

                x   = r6[1]
                xpr = r6[2]*p_norm
                y   = r6[3]
                ypr = r6[4]*p_norm
                
                g  = abs(k1)/(1.0 + r6[6])
                t  = sqrt(g)
                lt = le/2.0 * t
                
                if k1>0
                    MHD = cos(lt)
                    M12 = sin(lt)/t
                    M21 = -M12*g
                    MVD = cosh(lt)
                    M34 = sinh(lt)/t
                    M43 = M34*g
                else
                    MHD = cosh(lt)
                    M12 = sinh(lt)/t
                    M21 = M12*g
                    MVD = cos(lt)
                    M34 = sin(lt)/t
                    M43 = -M34*g
                end

                r6[1] =  MHD*x + M12*xpr
                r6[2] = (M21*x + MHD*xpr)/p_norm
                r6[3] =  MVD*y + M34*ypr
                r6[4] = (M43*y + MVD*ypr)/p_norm

                r6[5]+= g*(x*x*(le/2.0-MHD*M12)-y*y*(le/2.0-MVD*M34))/4.0
                r6[5]+= (xpr*xpr*(le/2.0+MHD*M12)+ypr*ypr*(le/2.0+MVD*M34))/4.0
                r6[5]+= (x*xpr*M12*M21 + y*ypr*M34*M43)/2.0

                # Misalignment at exit
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
    end
    return nothing
end

function pass!(ele::QUAD_SC, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam)
    # ele: KQUAD
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    lost_flags = particles.lost_flag
    K = calculate_K(particles, particles.current)
    lstep = ele.len/ele.Nsteps
    for i in 1:ele.Nsteps
        QuadLinearPass_SC!(r_in, lstep, ele.k1, ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, num_particles, lost_flags,
            particles.a, particles.b, particles.Nl, particles.Nm, K)
    end
    return nothing
end

function QuadLinearPass_SC_P!(r::Array{Float64,1}, le::Float64, k1::Float64, 
    T1::Array{Float64,1}, T2::Array{Float64,1}, R1::Array{Float64,2}, R2::Array{Float64,2}, 
    RApertures::Array{Float64,1}, EApertures::Array{Float64,1},
    num_particles::Int, lost_flags::Array{Int64,1}, a, b, Nl, Nm, K)
    # This is linear quadrupole pass function. Not used for symplectic tracking.

    # half-kick-half
    Threads.@threads for c in 1:num_particles
        if lost_flags[c] == 1
            continue
        end
        r6 = @view r[(c-1)*6+1:c*6]
        if !isnan(r6[1])
            p_norm = 1.0 / (1.0 + r6[6])

            if iszero(k1)
                drift6!(r6, le/2.0)
            else
                # Misalignment at entrance
                if !iszero(T1)
                    addvv!(r6, T1)
                end
                if !iszero(R1)
                    multmv!(r6, R1)
                end

                x   = r6[1]
                xpr = r6[2]*p_norm
                y   = r6[3]
                ypr = r6[4]*p_norm
                
                g  = abs(k1)/(1.0 + r6[6])
                t  = sqrt(g)
                lt = le/2.0 * t
                
                if k1>0
                    MHD = cos(lt)
                    M12 = sin(lt)/t
                    M21 = -M12*g
                    MVD = cosh(lt)
                    M34 = sinh(lt)/t
                    M43 = M34*g
                else
                    MHD = cosh(lt)
                    M12 = sinh(lt)/t
                    M21 = M12*g
                    MVD = cos(lt)
                    M34 = sin(lt)/t
                    M43 = -M34*g
                end

                r6[1] =  MHD*x + M12*xpr
                r6[2] = (M21*x + MHD*xpr)/p_norm
                r6[3] =  MVD*y + M34*ypr
                r6[4] = (M43*y + MVD*ypr)/p_norm

                r6[5]+= g*(x*x*(le/2.0-MHD*M12)-y*y*(le/2.0-MVD*M34))/4.0
                r6[5]+= (xpr*xpr*(le/2.0+MHD*M12)+ypr*ypr*(le/2.0+MVD*M34))/4.0
                r6[5]+= (x*xpr*M12*M21 + y*ypr*M34*M43)/2.0

                # Misalignment at exit
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
    end

    space_charge_P!(r, K, Nl, Nm, a/Nl, b/Nm, a, b, num_particles, le, lost_flags)

    Threads.@threads for c in 1:num_particles
        if lost_flags[c] == 1
            continue
        end
        r6 = @view r[(c-1)*6+1:c*6]
        if !isnan(r6[1])
            p_norm = 1.0 / (1.0 + r6[6])

            if iszero(k1)
                drift6!(r6, le/2.0)
            else
                # Misalignment at entrance
                if !iszero(T1)
                    addvv!(r6, T1)
                end
                if !iszero(R1)
                    multmv!(r6, R1)
                end

                x   = r6[1]
                xpr = r6[2]*p_norm
                y   = r6[3]
                ypr = r6[4]*p_norm
                
                g  = abs(k1)/(1.0 + r6[6])
                t  = sqrt(g)
                lt = le/2.0 * t
                
                if k1>0
                    MHD = cos(lt)
                    M12 = sin(lt)/t
                    M21 = -M12*g
                    MVD = cosh(lt)
                    M34 = sinh(lt)/t
                    M43 = M34*g
                else
                    MHD = cosh(lt)
                    M12 = sinh(lt)/t
                    M21 = M12*g
                    MVD = cos(lt)
                    M34 = sin(lt)/t
                    M43 = -M34*g
                end

                r6[1] =  MHD*x + M12*xpr
                r6[2] = (M21*x + MHD*xpr)/p_norm
                r6[3] =  MVD*y + M34*ypr
                r6[4] = (M43*y + MVD*ypr)/p_norm

                r6[5]+= g*(x*x*(le/2.0-MHD*M12)-y*y*(le/2.0-MVD*M34))/4.0
                r6[5]+= (xpr*xpr*(le/2.0+MHD*M12)+ypr*ypr*(le/2.0+MVD*M34))/4.0
                r6[5]+= (x*xpr*M12*M21 + y*ypr*M34*M43)/2.0

                # Misalignment at exit
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
    end
    return nothing
end

function pass_P!(ele::QUAD_SC, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam)
    # ele: KQUAD
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    lost_flags = particles.lost_flag
    K = calculate_K(particles, particles.current)
    lstep = ele.len/ele.Nsteps
    for i in 1:ele.Nsteps
        QuadLinearPass_SC!(r_in, lstep, ele.k1, ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, num_particles, lost_flags,
            particles.a, particles.b, particles.Nl, particles.Nm, K)
    end
    return nothing
end