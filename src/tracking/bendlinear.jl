function bend6!(r::AbstractVector{Float64}, L::Float64, b_angle::Float64, grd::Float64, ByError::Float64)
    # AT function. Ref[Terebilo, Andrei. "Accelerator modeling with MATLAB accelerator toolbox." PACS2001 (2001)].
    Kx = b_angle/L
    p_norm = 1/(1+r[6])
    G1 = (Kx*Kx+grd)*p_norm
    G2 = -grd*p_norm

    # predefine the matrix elements. maintain type stability for auto-differentiation.
    MHD = 0.0
    M12 = 0.0
    M21 = 0.0
    MVD = 0.0
    M34 = 0.0
    M43 = 0.0
    arg1 = 0.0
    arg2 = 0.0
    sqrtG1 = 0.0
    sqrtG2 = 0.0
    if G1 == 0.0
        MHD = 1.0
        M12 = L
        M21 = 0.0
    elseif G1 > 0.0
        sqrtG1 = sqrt(G1)
        arg1 = L*sqrtG1
        MHD = cos(arg1)
        M12 = sin(arg1)/sqrtG1
        M21 = -sin(arg1)*sqrtG1
    elseif G1 < 0.0
        sqrtG1 = sqrt(-G1)
        arg1 = L*sqrtG1
        MHD = cosh(arg1)
        M12 = sinh(arg1)/sqrtG1
        M21 = sinh(arg1)*sqrtG1
    end

    if G2 == 0.0
        MVD = 1.0
        M34 = L
        M43 = 0.0
    elseif G2 > 0.0	        
        sqrtG2 = sqrt(G2)
        arg2 = L*sqrtG2
        MVD = cos(arg2)
        M34 = sin(arg2)/sqrtG2
        M43 = -sin(arg2)*sqrtG2
    elseif G2 < 0.0	
        sqrtG2 = sqrt(-G2)
        arg2 = L*sqrtG2
        MVD = cosh(arg2)
        M34 = sinh(arg2)/sqrtG2
        M43 = sinh(arg2)*sqrtG2
    end

    x = r[1]
    xpr = r[2]*p_norm
    y = r[3]
    ypr = r[4]*p_norm
    delta = r[6]

    r[1] = MHD*x + M12*xpr
    r[2] = (M21*x + MHD*xpr)/p_norm

    if G1 == 0.0
        r[1] += (delta*p_norm-ByError)*L*L*Kx/2.0
        r[2] += (delta*p_norm-ByError)*L*Kx/p_norm 
    elseif G1 > 0.0
        r[1] += (delta*p_norm-ByError)*(1.0-cos(arg1))*Kx/G1
        r[2] += (delta*p_norm-ByError)*sin(arg1)*Kx/(sqrtG1*p_norm)
    elseif G1 < 0.0
        r[1] += (delta*p_norm-ByError)*(1.0-cosh(arg1))*Kx/G1
        r[2] += (delta*p_norm-ByError)*sinh(arg1)*Kx/(sqrtG1*p_norm) 
    end

    r[3]=  MVD*y + M34*ypr
    r[4]= (M43*y + MVD*ypr)/p_norm 
    r[5]+= xpr*xpr*(L+MHD*M12)/4.0

    if G1 != 0.0
        r[5]+= (L-MHD*M12)*(x*x*G1+(delta*p_norm-ByError)*(delta*p_norm-ByError)*Kx*Kx/G1-2.0*x*Kx*(delta*p_norm-ByError))/4.0
        r[5]+= M12*M21*(x*xpr - xpr*(delta*p_norm-ByError)*Kx/G1)/2.0
        r[5]+= Kx*x*M12+xpr*(1.0-MHD)*Kx/G1+(delta*p_norm-ByError)*(L-M12)*Kx*Kx/G1
    end
    r[5]+= ((L-MVD*M34)*y*y*G2 + ypr*ypr*(L+MVD*M34))/4.0
    r[5]+= M34*M43*x*xpr/2.0
end

function BendLinearPass!(r::Array{Float64,1}, le::Float64, grd::Float64, ba::Float64, bye::Float64,
    entrance_angle::Float64, exit_angle::Float64,
    fint1::Float64, fint2::Float64, gap::Float64, T1::Array{Float64,1}, T2::Array{Float64,1}, 
    R1::Array{Float64,2}, R2::Array{Float64,2}, RApertures::Array{Float64,1}, EApertures::Array{Float64,1}, 
    num_particles::Int, lost_flags::Array{Int64,1})
    # AT function. Ref[Terebilo, Andrei. "Accelerator modeling with MATLAB accelerator toolbox." PACS2001 (2001)].

    irho = ba/le

    for c in 1:num_particles
        if isone(lost_flags[c])
            continue
        end
        r6 = @view r[(c-1)*6+1:c*6]
        # Misalignment at entrance
        if !iszero(T1)
            addvv!(r6, T1)
        end
        if !iszero(R1)
            multmv!(r6, R1)
        end
        # Edge focus at entrance
        edge_fringe_entrance!(r6, irho, entrance_angle, fint1, gap, 1)

        bend6!(r6, le, ba, grd, bye)

        # Edge focus at exit
        edge_fringe_exit!(r6, irho, exit_angle, fint2, gap, 1)

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
    return nothing
end

function BendLinearPass_P!(r::Array{Float64,1}, le::Float64, grd::Float64, ba::Float64, bye::Float64,
    entrance_angle::Float64, exit_angle::Float64,
    fint1::Float64, fint2::Float64, gap::Float64, T1::Array{Float64,1}, T2::Array{Float64,1}, 
    R1::Array{Float64,2}, R2::Array{Float64,2}, RApertures::Array{Float64,1}, EApertures::Array{Float64,1}, 
    num_particles::Int, lost_flags::Array{Int64,1})
    # AT function. Ref[Terebilo, Andrei. "Accelerator modeling with MATLAB accelerator toolbox." PACS2001 (2001)].

    irho = ba/le

    Threads.@threads for c in 1:num_particles
        if isone(lost_flags[c])
            continue
        end
        r6 = @view r[(c-1)*6+1:c*6]
        # Misalignment at entrance
        if !iszero(T1)
            addvv!(r6, T1)
        end
        if !iszero(R1)
            multmv!(r6, R1)
        end
        # Edge focus at entrance
        edge_fringe_entrance!(r6, irho, entrance_angle, fint1, gap, 1)

        # Edge focus at exit
        edge_fringe_exit!(r6, irho, exit_angle, fint2, gap, 1)

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
    return nothing
end



function pass!(ele::LBEND, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam{Float64})
    # ele: LBEND
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    lost_flags = particles.lost_flag
    irho = ele.angle / ele.len
    BendLinearPass!(r_in, ele.len, ele.K, ele.angle, ele.ByError,
        ele.e1, ele.e2,
        ele.fint1, ele.fint2, ele.FullGap, ele.T1, ele.T2,
        ele.R1, ele.R2, ele.RApertures, ele.EApertures,
        num_particles, lost_flags)
    return nothing
end

function pass_P!(ele::LBEND, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam{Float64})
    # ele: LBEND
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    lost_flags = particles.lost_flag
    irho = ele.angle / ele.len
    BendLinearPass_P!(r_in, ele.len, ele.K, ele.angle, ele.ByError,
        ele.e1, ele.e2,
        ele.fint1, ele.fint2, ele.FullGap, ele.T1, ele.T2,
        ele.R1, ele.R2, ele.RApertures, ele.EApertures,
        num_particles, lost_flags)
    return nothing
end
