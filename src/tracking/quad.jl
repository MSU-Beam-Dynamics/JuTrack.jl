@inline function _quad_linear_transport_row!(r::Matrix{Float64}, c::Int, le::Float64, k1::Float64)
    @inbounds begin
        p_norm = 1.0 / (1.0 + r[c, 6])
        x = r[c, 1]
        xpr = r[c, 2] * p_norm
        y = r[c, 3]
        ypr = r[c, 4] * p_norm

        g = abs(k1) / (1.0 + r[c, 6])
        t = sqrt(g)
        lt = le * t

        if k1 > 0
            MHD = cos(lt)
            M12 = sin(lt) / t
            M21 = -M12 * g
            MVD = cosh(lt)
            M34 = sinh(lt) / t
            M43 = M34 * g
        else
            MHD = cosh(lt)
            M12 = sinh(lt) / t
            M21 = M12 * g
            MVD = cos(lt)
            M34 = sin(lt) / t
            M43 = -M34 * g
        end

        r[c, 1] = MHD * x + M12 * xpr
        r[c, 2] = (M21 * x + MHD * xpr) / p_norm
        r[c, 3] = MVD * y + M34 * ypr
        r[c, 4] = (M43 * y + MVD * ypr) / p_norm

        r[c, 5] += g * (x * x * (le - MHD * M12) - y * y * (le - MVD * M34)) / 4.0
        r[c, 5] += (xpr * xpr * (le + MHD * M12) + ypr * ypr * (le + MVD * M34)) / 4.0
        r[c, 5] += (x * xpr * M12 * M21 + y * ypr * M34 * M43) / 2.0
    end
    return nothing
end

function QuadLinearPass!(r::Matrix{Float64}, le::Float64, k1::Float64, beti::Float64,
    T1::Array{Float64,1}, T2::Array{Float64,1}, R1::Array{Float64,2}, R2::Array{Float64,2}, 
    RApertures::Array{Float64,1}, EApertures::Array{Float64,1},
    num_particles::Int, lost_flags::Array{Int64,1})
    # This is linear quadrupole pass function. Not used for symplectic tracking.

    
    for c in 1:num_particles
        if isone(lost_flags[c]) || isnan(r[c, 1])
            continue
        end
        if !iszero(T1)
            addvv_row!(r, c, T1)
        end
        if !iszero(R1)
            multmv_row!(r, c, R1)
        end

        if iszero(k1)
            drift6_row!(r, c, le, beti)
        else
            _quad_linear_transport_row!(r, c, le, k1)
        end

        if !iszero(R2)
            multmv_row!(r, c, R2)
        end
        if !iszero(T2)
            addvv_row!(r, c, T2)
        end
        if _check_lost_row(r, c) || _check_lost_aperture_row(r, c, RApertures, EApertures)
            lost_flags[c] = 1
        end
    end
    return nothing
end

function pass!(ele::QUAD, r_in::Matrix{Float64}, num_particles::Int64, particles::Beam{Float64})
    # ele: KQUAD
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    lost_flags = particles.lost_flag
    if use_exact_beti == 1
        beti = 1.0 / particles.beta
    else
        beti = 1.0 
    end
    QuadLinearPass!(r_in, ele.len, ele.k1, beti, ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, num_particles, lost_flags)
    
    return nothing
end

function QuadLinearPass_P!(r::Matrix{Float64}, le::Float64, k1::Float64, beti::Float64,
    T1::Array{Float64,1}, T2::Array{Float64,1}, R1::Array{Float64,2}, R2::Array{Float64,2}, 
    RApertures::Array{Float64,1}, EApertures::Array{Float64,1},
    num_particles::Int, lost_flags::Array{Int64,1})
    # Modified based on AT function. Ref[Terebilo, Andrei. "Accelerator modeling with MATLAB accelerator toolbox." PACS2001 (2001)].

    
    Threads.@threads for c in 1:num_particles
        if isone(lost_flags[c]) || isnan(r[c, 1])
            continue
        end
        if !iszero(T1)
            addvv_row!(r, c, T1)
        end
        if !iszero(R1)
            multmv_row!(r, c, R1)
        end

        if iszero(k1)
            drift6_row!(r, c, le, beti)
        else
            _quad_linear_transport_row!(r, c, le, k1)
        end

        if !iszero(R2)
            multmv_row!(r, c, R2)
        end
        if !iszero(T2)
            addvv_row!(r, c, T2)
        end
        if _check_lost_row(r, c) || _check_lost_aperture_row(r, c, RApertures, EApertures)
            lost_flags[c] = 1
        end
    end
    return nothing
end

function pass_P!(ele::QUAD, r_in::Matrix{Float64}, num_particles::Int64, particles::Beam{Float64})
    # ele: KQUAD
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    lost_flags = particles.lost_flag
    if use_exact_beti == 1
        beti = 1.0 / particles.beta
    else
        beti = 1.0 
    end
    QuadLinearPass_P!(r_in, ele.len, ele.k1, beti, ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, num_particles, lost_flags)
    
    return nothing
end

###############################
# TPSA
###############################

function QuadLinearPass!(r::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}, le::Float64, k1::Float64, beti::Float64,
    T1::Array{Float64,1}, T2::Array{Float64,1}, R1::Array{Float64,2}, R2::Array{Float64,2}) where {T, TPS_Dim, Max_TPS_Degree}

    p_norm = 1.0 / (1.0 + r[6])

    if iszero(k1)
        drift6!(r, le, beti)
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

function pass_TPSA!(ele::QUAD, r_in::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}; E0::Float64=0.0, m0::Float64=m_e) where {T, TPS_Dim, Max_TPS_Degree}
    # ele: KQUAD
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    gamma = E0 / m0
    beta = sqrt(1 - 1/gamma^2)
    if use_exact_beti == 1
        beti = 1.0 / beta
    else
        beti = 1.0 
    end
    QuadLinearPass!(r_in, ele.len, ele.k1, beti, ele.T1, ele.T2, ele.R1, ele.R2)
    
    return nothing
end

######################################
# space-charge
######################################
function QuadLinearPass_SC!(r::Matrix{Float64}, le::Float64, k1::Float64, beti::Float64, 
    T1::Array{Float64,1}, T2::Array{Float64,1}, R1::Array{Float64,2}, R2::Array{Float64,2}, 
    RApertures::Array{Float64,1}, EApertures::Array{Float64,1},
    num_particles::Int, lost_flags::Array{Int64,1}, a::Float64, b::Float64, Nl::Int, Nm::Int, K::Float64, Nsteps::Int)
    # This is linear quadrupole pass function. Not used for symplectic tracking.

    lstep = le/Nsteps
    for step in 1:Nsteps
        # half-kick-half
        for c in 1:num_particles
            if lost_flags[c] == 1
                continue
            end
            r6 = @view r[c, :]
            p_norm = 1.0 / (1.0 + r6[6])

            if step == 1
                # Misalignment at entrance
                if !iszero(T1)
                    addvv!(r6, T1)
                end
                if !iszero(R1)
                    multmv!(r6, R1)
                end
            end

            if iszero(k1)
                drift6!(r6, lstep/2.0, beti)
            else
                x   = r6[1]
                xpr = r6[2]*p_norm
                y   = r6[3]
                ypr = r6[4]*p_norm
                
                g  = abs(k1)/(1.0 + r6[6])
                t  = sqrt(g)
                lt = lstep/2.0 * t
                
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

                r6[5]+= g*(x*x*(lstep/2.0-MHD*M12)-y*y*(lstep/2.0-MVD*M34))/4.0
                r6[5]+= (xpr*xpr*(lstep/2.0+MHD*M12)+ypr*ypr*(lstep/2.0+MVD*M34))/4.0
                r6[5]+= (x*xpr*M12*M21 + y*ypr*M34*M43)/2.0

                if check_lost(r6) || check_lost_aperture(r6, RApertures, EApertures)
                    lost_flags[c] = 1
                end
            end
        end
    
        space_charge!(r, K, Nl, Nm, a/Nl, b/Nm, a, b, num_particles, lstep, lost_flags)
    
        for c in 1:num_particles
            if lost_flags[c] == 1
                continue
            end
            r6 = @view r[c, :]
            p_norm = 1.0 / (1.0 + r6[6])
    
            if iszero(k1)
                drift6!(r6, lstep/2.0, beti)
            else
                x   = r6[1]
                xpr = r6[2]*p_norm
                y   = r6[3]
                ypr = r6[4]*p_norm
                
                g  = abs(k1)/(1.0 + r6[6])
                t  = sqrt(g)
                lt = lstep/2.0 * t
                
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

                r6[5]+= g*(x*x*(lstep/2.0-MHD*M12)-y*y*(lstep/2.0-MVD*M34))/4.0
                r6[5]+= (xpr*xpr*(lstep/2.0+MHD*M12)+ypr*ypr*(lstep/2.0+MVD*M34))/4.0
                r6[5]+= (x*xpr*M12*M21 + y*ypr*M34*M43)/2.0
            end

            if step == Nsteps
                # Misalignment at exit
                if !iszero(R2)
                    multmv!(r6, R2)
                end
                if !iszero(T2)
                    addvv!(r6, T2)
                end
            end
            if check_lost(r6) || check_lost_aperture(r6, RApertures, EApertures)
                lost_flags[c] = 1
            end
        end        
    end
    return nothing
end

function pass!(ele::QUAD_SC, r_in::Matrix{Float64}, num_particles::Int64, particles::Beam{Float64})
    # ele: KQUAD
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    lost_flags = particles.lost_flag
    K = calculate_K(particles, particles.current)
    if use_exact_beti == 1
        beti = 1.0 / particles.beta
    else
        beti = 1.0 
    end
    QuadLinearPass_SC!(r_in, ele.len, ele.k1, beti, ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, num_particles, lost_flags,
            ele.a, ele.b, ele.Nl, ele.Nm, K, ele.Nsteps)
    return nothing
end

function QuadLinearPass_SC!(r::Matrix{DTPSAD{N, T}}, le::DTPSAD{N, T}, k1::DTPSAD{N, T}, beti::Float64, 
    T1::Array{DTPSAD{N, T},1}, T2::Array{DTPSAD{N, T},1}, R1::Array{DTPSAD{N, T},2}, R2::Array{DTPSAD{N, T},2}, 
    RApertures::Array{Float64,1}, EApertures::Array{Float64,1},
    num_particles::Int, lost_flags::Array{Int64,1}, a::DTPSAD{N, T}, b::DTPSAD{N, T}, Nl::Int, Nm::Int, K::DTPSAD{N, T}, Nsteps::Int) where {N, T}
    # This is linear quadrupole pass function. Not used for symplectic tracking.

    lstep = le/Nsteps
    for step in 1:Nsteps
        # half-kick-half
        @inbounds for c in 1:num_particles
            if lost_flags[c] == 1
                continue
            end
            r6 = @view r[c, :]
            p_norm = 1.0 / (1.0 + r6[6])

            if step == 1
                # Misalignment at entrance
                if !iszero(T1)
                    addvv!(r6, T1)
                end
                if !iszero(R1)
                    multmv!(r6, R1)
                end
            end

            if iszero(k1)
                drift6!(r6, lstep/2.0)
            else
                x   = r6[1]
                xpr = r6[2]*p_norm
                y   = r6[3]
                ypr = r6[4]*p_norm
                
                g  = abs(k1)/(1.0 + r6[6])
                t  = sqrt(g)
                lt = lstep/2.0 * t
                
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

                r6[5]+= g*(x*x*(lstep/2.0-MHD*M12)-y*y*(lstep/2.0-MVD*M34))/4.0
                r6[5]+= (xpr*xpr*(lstep/2.0+MHD*M12)+ypr*ypr*(lstep/2.0+MVD*M34))/4.0
                r6[5]+= (x*xpr*M12*M21 + y*ypr*M34*M43)/2.0

                if check_lost(r6) || check_lost_aperture(r6, RApertures, EApertures)
                    lost_flags[c] = 1
                end
            end
        end
    
        space_charge!(r, K, Nl, Nm, a/Nl, b/Nm, a, b, num_particles, lstep, lost_flags)

        @inbounds for c in 1:num_particles
            if lost_flags[c] == 1
                continue
            end
            r6 = @view r[c, :]
            p_norm = 1.0 / (1.0 + r6[6])
    
            if iszero(k1)
                drift6!(r6, lstep/2.0)
            else
                x   = r6[1]
                xpr = r6[2]*p_norm
                y   = r6[3]
                ypr = r6[4]*p_norm
                
                g  = abs(k1)/(1.0 + r6[6])
                t  = sqrt(g)
                lt = lstep/2.0 * t
                
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

                r6[5]+= g*(x*x*(lstep/2.0-MHD*M12)-y*y*(lstep/2.0-MVD*M34))/4.0
                r6[5]+= (xpr*xpr*(lstep/2.0+MHD*M12)+ypr*ypr*(lstep/2.0+MVD*M34))/4.0
                r6[5]+= (x*xpr*M12*M21 + y*ypr*M34*M43)/2.0
            end

            if step == Nsteps
                # Misalignment at exit
                if !iszero(R2)
                    multmv!(r6, R2)
                end
                if !iszero(T2)
                    addvv!(r6, T2)
                end
            end
            if check_lost(r6) || check_lost_aperture(r6, RApertures, EApertures)
                lost_flags[c] = 1
            end
        end        
    end
    return nothing
end

function pass!(ele::QUAD_SC{DTPSAD{N, T}}, r_in::Matrix{DTPSAD{N, T}}, num_particles::Int64, particles::Beam{DTPSAD{N, T}}) where {N, T}
    # ele: KQUAD
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    lost_flags = particles.lost_flag
    K = calculate_K(particles, particles.current)
    if use_exact_beti == 1
        beti = 1.0 / particles.beta.val
    else
        beti = 1.0 
    end
    QuadLinearPass_SC!(r_in, ele.len, ele.k1, beti, ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, num_particles, lost_flags,
            ele.a, ele.b, ele.Nl, ele.Nm, K, ele.Nsteps)
    return nothing
end

@inline function _quad_sc2p5d_quad1!(r6::AbstractVector{Float64}, length::Float64, k1::Float64,
    beti::Float64, gamma2i::Float64)
    if iszero(k1)
        drift6!(r6, length, beti, gamma2i)
        return nothing
    end

    dp_p = r6[6]
    x = r6[1]
    px = r6[2]
    y = r6[3]
    py = r6[4]

    sqrt_k1 = sqrt(abs(k1))
    lt = sqrt_k1 * length
    if k1 > 0.0
        cx = cos(lt)
        sx = sin(lt)
        cy = cosh(lt)
        sy = sinh(lt)
        m11 = cx
        m12 = sx / sqrt_k1
        m21 = -sx * sqrt_k1
        m22 = cx
        m33 = cy
        m34 = sy / sqrt_k1
        m43 = sy * sqrt_k1
        m44 = cy
    else
        cx = cosh(lt)
        sx = sinh(lt)
        cy = cos(lt)
        sy = sin(lt)
        m11 = cx
        m12 = sx / sqrt_k1
        m21 = sx * sqrt_k1
        m22 = cx
        m33 = cy
        m34 = sy / sqrt_k1
        m43 = -sy * sqrt_k1
        m44 = cy
    end

    r6[1] = x * m11 + px * m12
    r6[2] = x * m21 + px * m22
    r6[3] = y * m33 + py * m34
    r6[4] = y * m43 + py * m44
    # JuTrack's internal z convention is opposite to PyORBIT's transport sign.
    r6[5] -= dp_p * gamma2i * length
    return nothing
end

@inline function _quad_sc2p5d_quad1_row!(r::Matrix{Float64}, c::Int, length::Float64, k1::Float64,
    beti::Float64, gamma2i::Float64)
    if iszero(k1)
        drift6_row!(r, c, length, beti, gamma2i)
        return nothing
    end

    @inbounds begin
        dp_p = r[c, 6]
        x = r[c, 1]
        px = r[c, 2]
        y = r[c, 3]
        py = r[c, 4]

        sqrt_k1 = sqrt(abs(k1))
        lt = sqrt_k1 * length
        if k1 > 0.0
            cx = cos(lt)
            sx = sin(lt)
            cy = cosh(lt)
            sy = sinh(lt)
            m11 = cx
            m12 = sx / sqrt_k1
            m21 = -sx * sqrt_k1
            m22 = cx
            m33 = cy
            m34 = sy / sqrt_k1
            m43 = sy * sqrt_k1
            m44 = cy
        else
            cx = cosh(lt)
            sx = sinh(lt)
            cy = cos(lt)
            sy = sin(lt)
            m11 = cx
            m12 = sx / sqrt_k1
            m21 = sx * sqrt_k1
            m22 = cx
            m33 = cy
            m34 = sy / sqrt_k1
            m43 = -sy * sqrt_k1
            m44 = cy
        end

        r[c, 1] = x * m11 + px * m12
        r[c, 2] = x * m21 + px * m22
        r[c, 3] = y * m33 + py * m34
        r[c, 4] = y * m43 + py * m44
        r[c, 5] -= dp_p * gamma2i * length
    end
    return nothing
end

@inline function _quad_sc2p5d_quad2!(r6::AbstractVector{Float64}, length::Float64, gamma2i::Float64)
    dp_p = r6[6]
    knl = 1.0 / (1.0 + dp_p)
    r6[1] -= knl * length * dp_p * r6[2]
    r6[3] -= knl * length * dp_p * r6[4]
    phifac = (r6[2] * r6[2] + r6[4] * r6[4] + dp_p * dp_p * gamma2i) / 2.0
    phifac = (phifac * knl + dp_p * dp_p * gamma2i) * knl
    r6[5] += length * phifac
    return nothing
end

@inline function _quad_sc2p5d_quad2_row!(r::Matrix{Float64}, c::Int, length::Float64, gamma2i::Float64)
    @inbounds begin
        dp_p = r[c, 6]
        knl = 1.0 / (1.0 + dp_p)
        r[c, 1] -= knl * length * dp_p * r[c, 2]
        r[c, 3] -= knl * length * dp_p * r[c, 4]
        phifac = (r[c, 2] * r[c, 2] + r[c, 4] * r[c, 4] + dp_p * dp_p * gamma2i) / 2.0
        phifac = (phifac * knl + dp_p * dp_p * gamma2i) * knl
        r[c, 5] += length * phifac
    end
    return nothing
end

@inline function _quad_sc2p5d_quad1!(r6::AbstractVector{DTPSAD{N, T}}, length::DTPSAD{N, T},
    k1::DTPSAD{N, T}, gamma2i::Float64) where {N, T}
    if iszero(k1)
        drift6!(r6, length, gamma2i)
        return nothing
    end

    dp_p = r6[6]
    x = r6[1]
    px = r6[2]
    y = r6[3]
    py = r6[4]

    sqrt_k1 = sqrt(abs(k1))
    lt = sqrt_k1 * length
    if k1 > 0.0
        cx = cos(lt)
        sx = sin(lt)
        cy = cosh(lt)
        sy = sinh(lt)
        m11 = cx
        m12 = sx / sqrt_k1
        m21 = -sx * sqrt_k1
        m22 = cx
        m33 = cy
        m34 = sy / sqrt_k1
        m43 = sy * sqrt_k1
        m44 = cy
    else
        cx = cosh(lt)
        sx = sinh(lt)
        cy = cos(lt)
        sy = sin(lt)
        m11 = cx
        m12 = sx / sqrt_k1
        m21 = sx * sqrt_k1
        m22 = cx
        m33 = cy
        m34 = sy / sqrt_k1
        m43 = -sy * sqrt_k1
        m44 = cy
    end

    r6[1] = x * m11 + px * m12
    r6[2] = x * m21 + px * m22
    r6[3] = y * m33 + py * m34
    r6[4] = y * m43 + py * m44
    r6[5] -= dp_p * gamma2i * length
    return nothing
end

@inline function _quad_sc2p5d_quad2!(r6::AbstractVector{DTPSAD{N, T}}, length::DTPSAD{N, T},
    gamma2i::Float64) where {N, T}
    dp_p = r6[6]
    knl = 1.0 / (1.0 + dp_p)
    r6[1] -= knl * length * dp_p * r6[2]
    r6[3] -= knl * length * dp_p * r6[4]
    phifac = (r6[2] * r6[2] + r6[4] * r6[4] + dp_p * dp_p * gamma2i) / 2.0
    phifac = (phifac * knl + dp_p * dp_p * gamma2i) * knl
    r6[5] += length * phifac
    return nothing
end

function QuadLinearPass_SC2P5D!(ele::QUAD_SC2P5D{Float64}, r::Matrix{Float64}, le::Float64, k1::Float64,
    beti::Float64, gamma2i::Float64, T1::Array{Float64,1}, T2::Array{Float64,1}, R1::Array{Float64,2},
    R2::Array{Float64,2}, RApertures::Array{Float64,1}, EApertures::Array{Float64,1},
    num_particles::Int, lost_flags::Array{Int64,1}, particles::Beam{Float64})

    lstep = le / ele.Nsteps

    for c in 1:num_particles
        if isone(lost_flags[c])
            continue
        end
        if !isnan(r[c, 1])
            if !iszero(T1)
                addvv_row!(r, c, T1)
            end
            if !iszero(R1)
                multmv_row!(r, c, R1)
            end
            if _check_lost_row(r, c) || _check_lost_aperture_row(r, c, RApertures, EApertures)
                lost_flags[c] = 1
            end
        else
            lost_flags[c] = 1
        end
    end

    for step in 1:ele.Nsteps
        _sc2p5d_track_step!(ele, lstep, r, particles)

        for c in 1:num_particles
            if isone(lost_flags[c])
                continue
            end
            if !isnan(r[c, 1])
                half_step = lstep / 2.0
                _quad_sc2p5d_quad1_row!(r, c, half_step, k1, beti, gamma2i)
                _quad_sc2p5d_quad2_row!(r, c, half_step, gamma2i)
                _quad_sc2p5d_quad2_row!(r, c, half_step, gamma2i)
                _quad_sc2p5d_quad1_row!(r, c, half_step, k1, beti, gamma2i)

                if step == ele.Nsteps
                    if !iszero(R2)
                        multmv_row!(r, c, R2)
                    end
                    if !iszero(T2)
                        addvv_row!(r, c, T2)
                    end
                end
                if _check_lost_row(r, c) || _check_lost_aperture_row(r, c, RApertures, EApertures)
                    lost_flags[c] = 1
                end
            else
                lost_flags[c] = 1
            end
        end
    end
    return nothing
end

function pass!(ele::QUAD_SC2P5D{Float64}, r_in::Matrix{Float64}, num_particles::Int64, particles::Beam{Float64})
    lost_flags = particles.lost_flag
    if use_exact_beti == 1
        beti = 1.0 / particles.beta
    else
        beti = 1.0
    end
    gamma2i = use_exact_Hamiltonian == 2 ? 1.0 / (particles.gamma^2) : 0.0
    QuadLinearPass_SC2P5D!(ele, r_in, ele.len, ele.k1, beti, gamma2i, ele.T1, ele.T2, ele.R1, ele.R2,
        ele.RApertures, ele.EApertures, num_particles, lost_flags, particles)
    return nothing
end

function QuadLinearPass_SC2P5D!(ele::QUAD_SC2P5D{DTPSAD{N, T}}, r::Matrix{DTPSAD{N, T}}, le::DTPSAD{N, T},
    k1::DTPSAD{N, T}, gamma2i::Float64, T1::Array{DTPSAD{N, T},1}, T2::Array{DTPSAD{N, T},1},
    R1::Array{DTPSAD{N, T},2}, R2::Array{DTPSAD{N, T},2}, RApertures::Array{Float64,1},
    EApertures::Array{Float64,1}, num_particles::Int, lost_flags::Array{Int64,1},
    particles::Beam{DTPSAD{N, T}}) where {N, T}

    lstep = le / ele.Nsteps

    @inbounds for c in 1:num_particles
        if isone(lost_flags[c])
            continue
        end
        r6 = @view r[c, :]
        if !isnan(r6[1])
            if !iszero(T1)
                addvv!(r6, T1)
            end
            if !iszero(R1)
                multmv!(r6, R1)
            end
            if check_lost(r6) || check_lost_aperture(r6, RApertures, EApertures)
                lost_flags[c] = 1
            end
        else
            lost_flags[c] = 1
        end
    end

    for step in 1:ele.Nsteps
        _sc2p5d_track_step!(ele, lstep, r, particles)

        @inbounds for c in 1:num_particles
            if isone(lost_flags[c])
                continue
            end
            r6 = @view r[c, :]
            if !isnan(r6[1])
                half_step = lstep / 2.0
                _quad_sc2p5d_quad1!(r6, half_step, k1, gamma2i)
                _quad_sc2p5d_quad2!(r6, half_step, gamma2i)
                _quad_sc2p5d_quad2!(r6, half_step, gamma2i)
                _quad_sc2p5d_quad1!(r6, half_step, k1, gamma2i)

                if step == ele.Nsteps
                    if !iszero(R2)
                        multmv!(r6, R2)
                    end
                    if !iszero(T2)
                        addvv!(r6, T2)
                    end
                end
                if check_lost(r6) || check_lost_aperture(r6, RApertures, EApertures)
                    lost_flags[c] = 1
                end
            else
                lost_flags[c] = 1
            end
        end
    end
    return nothing
end

function pass!(ele::QUAD_SC2P5D{DTPSAD{N, T}}, r_in::Matrix{DTPSAD{N, T}}, num_particles::Int64,
    particles::Beam{DTPSAD{N, T}}) where {N, T}
    lost_flags = particles.lost_flag
    gamma2i = use_exact_Hamiltonian == 2 ? 1.0 / (particles.gamma.val^2) : 0.0
    QuadLinearPass_SC2P5D!(ele, r_in, ele.len, ele.k1, gamma2i, ele.T1, ele.T2, ele.R1, ele.R2,
        ele.RApertures, ele.EApertures, num_particles, lost_flags, particles)
    return nothing
end

function QuadLinearPass_SC_P!(r::Matrix{Float64}, le::Float64, k1::Float64, beti::Float64,
    T1::Array{Float64,1}, T2::Array{Float64,1}, R1::Array{Float64,2}, R2::Array{Float64,2}, 
    RApertures::Array{Float64,1}, EApertures::Array{Float64,1},
    num_particles::Int, lost_flags::Array{Int64,1}, a::Float64, b::Float64, Nl::Int, Nm::Int, K::Float64, Nsteps::Int)
    # This is linear quadrupole pass function. Not used for symplectic tracking.

    lstep = le/Nsteps
    for step in 1:Nsteps
        # half-kick-half
        Threads.@threads for c in 1:num_particles
            if lost_flags[c] == 1
                continue
            end
            r6 = @view r[c, :]
            p_norm = 1.0 / (1.0 + r6[6])

            if step == 1
                # Misalignment at entrance
                if !iszero(T1)
                    addvv!(r6, T1)
                end
                if !iszero(R1)
                    multmv!(r6, R1)
                end
            end

            if iszero(k1)
                drift6!(r6, lstep/2.0, beti)
            else
                x   = r6[1]
                xpr = r6[2]*p_norm
                y   = r6[3]
                ypr = r6[4]*p_norm
                
                g  = abs(k1)/(1.0 + r6[6])
                t  = sqrt(g)
                lt = lstep/2.0 * t
                
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

                r6[5]+= g*(x*x*(lstep/2.0-MHD*M12)-y*y*(lstep/2.0-MVD*M34))/4.0
                r6[5]+= (xpr*xpr*(lstep/2.0+MHD*M12)+ypr*ypr*(lstep/2.0+MVD*M34))/4.0
                r6[5]+= (x*xpr*M12*M21 + y*ypr*M34*M43)/2.0

                if check_lost(r6) || check_lost_aperture(r6, RApertures, EApertures)
                    lost_flags[c] = 1
                end
            end
        end
    
        space_charge_P!(r, K, Nl, Nm, a/Nl, b/Nm, a, b, num_particles, lstep, lost_flags)
    
        Threads.@threads for c in 1:num_particles
            if lost_flags[c] == 1
                continue
            end
            r6 = @view r[c, :]
            p_norm = 1.0 / (1.0 + r6[6])
    
            if iszero(k1)
                drift6!(r6, lstep/2.0, beti)
            else
                x   = r6[1]
                xpr = r6[2]*p_norm
                y   = r6[3]
                ypr = r6[4]*p_norm
                
                g  = abs(k1)/(1.0 + r6[6])
                t  = sqrt(g)
                lt = lstep/2.0 * t
                
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

                r6[5]+= g*(x*x*(lstep/2.0-MHD*M12)-y*y*(lstep/2.0-MVD*M34))/4.0
                r6[5]+= (xpr*xpr*(lstep/2.0+MHD*M12)+ypr*ypr*(lstep/2.0+MVD*M34))/4.0
                r6[5]+= (x*xpr*M12*M21 + y*ypr*M34*M43)/2.0
            end

            if step == Nsteps
                # Misalignment at exit
                if !iszero(R2)
                    multmv!(r6, R2)
                end
                if !iszero(T2)
                    addvv!(r6, T2)
                end
            end
            if check_lost(r6) || check_lost_aperture(r6, RApertures, EApertures)
                lost_flags[c] = 1
            end
        end        
    end
    return nothing
end

function pass_P!(ele::QUAD_SC, r_in::Matrix{Float64}, num_particles::Int64, particles::Beam{Float64})
    # ele: KQUAD
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    lost_flags = particles.lost_flag
    K = calculate_K(particles, particles.current)
    if use_exact_beti == 1
        beti = 1.0 / particles.beta
    else
        beti = 1.0 
    end
    QuadLinearPass_SC!(r_in, ele.len, ele.k1, beti, ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, num_particles, lost_flags,
            ele.a, ele.b, ele.Nl, ele.Nm, K, ele.Nsteps)
    return nothing
end

function pass_P!(ele::QUAD_SC2P5D, r_in::Matrix{Float64}, num_particles::Int64, particles::Beam{Float64})
    pass!(ele, r_in, num_particles, particles)
    return nothing
end
