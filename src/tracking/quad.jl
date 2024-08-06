function QuadLinearPass!(r::Array{Float64,1}, le::Float64, k1::Float64, 
    T1::Array{Float64,1}, T2::Array{Float64,1}, R1::Array{Float64,2}, R2::Array{Float64,2}, 
    RApertures::Array{Float64,1}, EApertures::Array{Float64,1},
    num_particles::Int, lost_flags::Array{Int64,1})
    # This is linear quadrupole pass function. Not used for symplectic tracking.

    
    for c in 1:num_particles
        if lost_flags[c] == 1
            continue
        end
        r6 = @view r[(c-1)*6+1:c*6]
        if !isnan(r6[1])
            p_norm = 1.0 / (1.0 + r6[6])

            if iszero(k1)
                drift6!(r6, le)
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
                lt = le * t
                
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

                r6[5]+= g*(x*x*(le-MHD*M12)-y*y*(le-MVD*M34))/4.0
                r6[5]+= (xpr*xpr*(le+MHD*M12)+ypr*ypr*(le+MVD*M34))/4.0
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

function pass!(ele::QUAD, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam)
    # ele: KQUAD
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    lost_flags = particles.lost_flag
    QuadLinearPass!(r_in, ele.len, ele.k1, ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, num_particles, lost_flags)
    
    return nothing
end

function QuadLinearPass_P!(r::Array{Float64,1}, le::Float64, k1::Float64, 
    T1::Array{Float64,1}, T2::Array{Float64,1}, R1::Array{Float64,2}, R2::Array{Float64,2}, 
    RApertures::Array{Float64,1}, EApertures::Array{Float64,1},
    num_particles::Int, lost_flags::Array{Int64,1})
    # Modified based on AT function. Ref[Terebilo, Andrei. "Accelerator modeling with MATLAB accelerator toolbox." PACS2001 (2001)].

    
    Threads.@threads for c in 1:num_particles
        if lost_flags[c] == 1
            continue
        end
        r6 = @view r[(c-1)*6+1:c*6]
        if !isnan(r6[1])
            p_norm = 1.0 / (1.0 + r6[6])

            if iszero(k1)
                drift6!(r6, le)
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
                lt = le * t
                
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

                r6[5]+= g*(x*x*(le-MHD*M12)-y*y*(le-MVD*M34))/4.0
                r6[5]+= (xpr*xpr*(le+MHD*M12)+ypr*ypr*(le+MVD*M34))/4.0
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

function pass_P!(ele::QUAD, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam)
    # ele: KQUAD
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    lost_flags = particles.lost_flag
    QuadLinearPass_P!(r_in, ele.len, ele.k1, ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, num_particles, lost_flags)
    
    return nothing
end

###############################
# TPSA
###############################

function QuadLinearPass!(r::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}, le::Float64, k1::Float64, 
    T1::Array{Float64,1}, T2::Array{Float64,1}, R1::Array{Float64,2}, R2::Array{Float64,2}) where {T, TPS_Dim, Max_TPS_Degree}

    p_norm = 1.0 / (1.0 + r[6])

    if iszero(k1)
        drift6!(r, le)
    else
        # Misalignment at entrance
        if !iszero(T1)
            addvv!(r, T1)
        end
        if !iszero(R1)
            multmv!(r, R1)
        end

        x   = r[1]
        xpr = r[2]*p_norm
        y   = r[3]
        ypr = r[4]*p_norm
        
        g  = abs(k1)/(1.0 + r[6])
        t  = sqrt(g)
        lt = le * t
        
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

        r[1] =  MHD*x + M12*xpr
        r[2] = (M21*x + MHD*xpr)/p_norm
        r[3] =  MVD*y + M34*ypr
        r[4] = (M43*y + MVD*ypr)/p_norm

        r[5]+= g*(x*x*(le-MHD*M12)-y*y*(le-MVD*M34))/4.0
        r[5]+= (xpr*xpr*(le+MHD*M12)+ypr*ypr*(le+MVD*M34))/4.0
        r[5]+= (x*xpr*M12*M21 + y*ypr*M34*M43)/2.0

        # Misalignment at exit
        if !iszero(R2)
            multmv!(r6, R2)
        end
        if !iszero(T2)
            addvv!(r6, T2)
        end
    end
        
    return nothing
end

function pass_TPSA!(ele::QUAD, r_in::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}) where {T, TPS_Dim, Max_TPS_Degree}
    # ele: KQUAD
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    QuadLinearPass!(r_in, ele.len, ele.k1, ele.T1, ele.T2, ele.R1, ele.R2)
    
    return nothing
end