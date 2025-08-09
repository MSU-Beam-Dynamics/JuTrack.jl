# include("fringe.jl")
function B2perp(bx, by, irho, x, xpr, y, ypr)
    v_norm2 = 1.0 / ((1.0 + x*irho)^2 + xpr^2 + ypr^2)
    return ((by * (1.0 + x*irho))^2 + (bx *(1.0 + x*irho))^2 + (bx*ypr - by*xpr)^2) * v_norm2
end
function bndthinkick!(r::AbstractVector{Float64}, A::Array{Float64,1}, B::Array{Float64,1}, 
    L::Float64, irho::Float64, max_order::Int, beti::Float64)
    # AT function. Ref[Terebilo, Andrei. "Accelerator modeling with MATLAB accelerator toolbox." PACS2001 (2001)].
    ReSum = B[max_order + 1]
    ImSum = A[max_order + 1]
    ReSumTemp = 0.0

    for i in max_order-1:-1:0
        ReSumTemp = ReSum * r[1] - ImSum * r[3] + B[i+1]
        ImSum = ImSum * r[1] + ReSum * r[3] + A[i+1]
        ReSum = ReSumTemp
    end

    r[2] -= L * (ReSum - (r[6] - r[1] * irho) * irho)
    r[4] += L * ImSum
    r[5] += L * irho * r[1] * beti # Path length
    return nothing
end

function bndthinkickrad!(r::AbstractVector{Float64}, A::Array{Float64,1}, B::Array{Float64,1}, 
    L::Float64, irho::Float64, E0::Float64, max_order::Int, beti::Float64)
    # Modified based on AT function. Ref[Terebilo, Andrei. "Accelerator modeling with MATLAB accelerator toolbox." PACS2001 (2001)].
    ReSum = B[max_order + 1]
    ImSum = A[max_order + 1]
    ReSumTemp = 0.0
    CRAD = CGAMMA * E0^3 / (2.0*pi*1e27) # [m]/[GeV^3] M.Sands (4.1)

    for i in max_order-1:-1:0
        ReSumTemp = ReSum * r[1] - ImSum * r[3] + B[i+1]
        ImSum = ImSum * r[1] + ReSum * r[3] + A[i+1]
        ReSum = ReSumTemp
    end

    # angles from momentums
    p_norm = 1.0 / (1.0*beti + r[6])
    x = r[1]
    xpr = r[2] * p_norm
    y = r[3]
    ypr = r[4] * p_norm
    B2P = B2perp(ImSum, ReSum + irho, irho, x, xpr, y, ypr)

    dp_0 = r[6]
    r[6] = r[6] - CRAD * (1.0*beti+r[6])^2 * B2P * (1.0 + x*irho + (xpr^2 + ypr^2) / 2.0) * L
    
    # momentums after losing energy
    p_norm = 1.0 / (1.0*beti + r[6])
    r[2] = xpr / p_norm
    r[4] = ypr / p_norm

    r[2] -= L * (ReSum - (dp_0 - r[1]*irho)*irho)
    r[4] += L * ImSum
    r[5] += L * irho * r[1] * beti # Path length
    return nothing
end

function BendSymplecticPassRad!(r::Array{Float64,1}, le::Float64, beti::Float64, irho::Float64, A::Array{Float64,1}, B::Array{Float64,1}, 
    max_order::Int, num_int_steps::Int, entrance_angle::Float64, exit_angle::Float64, FringeBendEntrance::Int, FringeBendExit::Int,
    fint1::Float64, fint2::Float64, gap::Float64, FringeQuadEntrance::Int, FringeQuadExit::Int,
    fringeIntM0::Array{Float64,1}, fringeIntP0::Array{Float64,1}, T1::Array{Float64,1}, T2::Array{Float64,1}, 
    R1::Array{Float64,2}, R2::Array{Float64,2}, RApertures::Array{Float64,1}, EApertures::Array{Float64,1},
    KickAngle::Array{Float64,1}, E0::Float64, num_particles::Int, lost_flags::Array{Int64,1})
    # AT function. Ref[Terebilo, Andrei. "Accelerator modeling with MATLAB accelerator toolbox." PACS2001 (2001)].

    DRIFT1 = 0.6756035959798286638
    DRIFT2 = -0.1756035959798286639
    KICK1 = 1.351207191959657328
    KICK2 = -1.702414383919314656
    SL = le / num_int_steps
    L1 = SL * DRIFT1
    L2 = SL * DRIFT2
    K1 = SL * KICK1
    K2 = SL * KICK2

    if FringeQuadEntrance==2# && !isnothing(fringeIntM0) && !isnothing(fringeIntP0)
        useLinFrEleEntrance = 1
    else
        useLinFrEleEntrance = 0
    end
    if FringeQuadExit==2# && !isnothing(fringeIntM0) && !isnothing(fringeIntP0)
        useLinFrEleExit = 1
    else
        useLinFrEleExit = 0
    end

    B[1] -= sin(KickAngle[1]) / le
    A[1] += sin(KickAngle[2]) / le


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
        edge_fringe_entrance!(r6, irho, entrance_angle, fint1, gap, FringeBendEntrance)

        # Quadrupole gradient fringe entrance
        if !iszero(FringeQuadEntrance) && !iszero(B[2])
            if useLinFrEleEntrance == 1
                linearQuadFringeElegantEntrance!(r6, B[2], fringeIntM0, fringeIntP0)
            else
                QuadFringePassP!(r6, B[2])
            end
        end

        # Integrator
        for m in 1:num_int_steps
            drift6!(r6, L1, beti)
            bndthinkickrad!(r6, A, B, K1, irho, E0, max_order, beti)
            drift6!(r6, L2, beti)
            bndthinkickrad!(r6, A, B, K2, irho, E0, max_order, beti)
            drift6!(r6, L2, beti)
            bndthinkickrad!(r6, A, B, K1, irho, E0, max_order, beti)
            drift6!(r6, L1, beti)
        end

        # Quadrupole gradient fringe exit
        if !iszero(FringeQuadExit) && !iszero(B[2])
            if useLinFrEleExit == 1
                linearQuadFringeElegantExit!(r6, B[2], fringeIntM0, fringeIntP0)
            else
                QuadFringePassN!(r6, B[2])
            end
        end

        # Edge focus at exit
        edge_fringe_exit!(r6, irho, exit_angle, fint2, gap, FringeBendExit)

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

    
    B[1] += sin(KickAngle[1]) / le
    A[1] -= sin(KickAngle[2]) / le
    return nothing
end

function BendSymplecticPass!(r::Array{Float64,1}, le::Float64, beti::Float64, irho::Float64, A::Array{Float64,1}, B::Array{Float64,1}, 
    max_order::Int, num_int_steps::Int, entrance_angle::Float64, exit_angle::Float64, FringeBendEntrance::Int, FringeBendExit::Int,
    fint1::Float64, fint2::Float64, gap::Float64, FringeQuadEntrance::Int, FringeQuadExit::Int,
    fringeIntM0::Array{Float64,1}, fringeIntP0::Array{Float64,1}, T1::Array{Float64,1}, T2::Array{Float64,1}, 
    R1::Array{Float64,2}, R2::Array{Float64,2}, RApertures::Array{Float64,1}, EApertures::Array{Float64,1},
    KickAngle::Array{Float64,1}, num_particles::Int, lost_flags::Array{Int64,1})
    # Modified based on AT function. Ref[Terebilo, Andrei. "Accelerator modeling with MATLAB accelerator toolbox." PACS2001 (2001)].

    DRIFT1 = 0.6756035959798286638
    DRIFT2 = -0.1756035959798286639
    KICK1 = 1.351207191959657328
    KICK2 = -1.702414383919314656
    SL = le / num_int_steps
    L1 = SL * DRIFT1
    L2 = SL * DRIFT2
    K1 = SL * KICK1
    K2 = SL * KICK2

    if FringeQuadEntrance==2# && !isnothing(fringeIntM0) && !isnothing(fringeIntP0)
        useLinFrEleEntrance = 1
    else
        useLinFrEleEntrance = 0
    end
    if FringeQuadExit==2 #&& !isnothing(fringeIntM0) && !isnothing(fringeIntP0)
        useLinFrEleExit = 1
    else
        useLinFrEleExit = 0
    end

    B[1] -= sin(KickAngle[1]) / le
    A[1] += sin(KickAngle[2]) / le


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
        edge_fringe_entrance!(r6, irho, entrance_angle, fint1, gap, FringeBendEntrance)

        # Quadrupole gradient fringe entrance
        if !iszero(FringeQuadEntrance) && !iszero(B[2])
            if isone(useLinFrEleEntrance)
                linearQuadFringeElegantEntrance!(r6, B[2], fringeIntM0, fringeIntP0)
            else
                QuadFringePassP!(r6, B[2])
            end
        end

        # Integrator
        for m in 1:num_int_steps
            drift6!(r6, L1, beti)
            bndthinkick!(r6, A, B, K1, irho, max_order, beti)
            drift6!(r6, L2, beti)
            bndthinkick!(r6, A, B, K2, irho, max_order, beti)
            drift6!(r6, L2, beti)
            bndthinkick!(r6, A, B, K1, irho, max_order, beti)
            drift6!(r6, L1, beti)
        end

        # Quadrupole gradient fringe exit
        if !iszero(FringeQuadExit) && !iszero(B[2])
            if isone(useLinFrEleExit)
                linearQuadFringeElegantExit!(r6, B[2], fringeIntM0, fringeIntP0)
            else
                QuadFringePassN!(r6, B[2])
            end
        end

        # Edge focus at exit
        edge_fringe_exit!(r6, irho, exit_angle, fint2, gap, FringeBendExit)

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

    B[1] += sin(KickAngle[1]) / le
    A[1] -= sin(KickAngle[2]) / le
    return nothing
end


function pass!(ele::SBEND, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam{Float64})
    # ele: SBEND
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    # particles: Beam{Float64} object
    lost_flags = particles.lost_flag
    irho = ele.angle / ele.len
    E0 = particles.energy
    if use_exact_beti == 1
        beti = 1.0 / particles.beta
    else
        beti = 1.0 
    end
    if ele.rad == 0
        BendSymplecticPass!(r_in, ele.len, beti, irho, ele.PolynomA, ele.PolynomB, ele.MaxOrder, ele.NumIntSteps,
            ele.e1, ele.e2,
            ele.FringeBendEntrance, ele.FringeBendExit,
            ele.fint1, ele.fint2, ele.gap,
            ele.FringeQuadEntrance, ele.FringeQuadExit,
            ele.FringeIntM0, ele.FringeIntP0,
            ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures,
            ele.KickAngle, num_particles, lost_flags)
    else
        BendSymplecticPassRad!(r_in, ele.len, beti, irho, ele.PolynomA, ele.PolynomB, ele.MaxOrder, ele.NumIntSteps,
            ele.e1, ele.e2,
            ele.FringeBendEntrance, ele.FringeBendExit,
            ele.fint1, ele.fint2, ele.gap,
            ele.FringeQuadEntrance, ele.FringeQuadExit,
            ele.FringeIntM0, ele.FringeIntP0,
            ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures,
            ele.KickAngle, E0, num_particles, lost_flags)
    end
    return nothing
end

######################
# multi-threading
function BendSymplecticPassRad_P!(r::Array{Float64,1}, le::Float64, beti::Float64, irho::Float64, A::Array{Float64,1}, B::Array{Float64,1}, 
    max_order::Int, num_int_steps::Int, entrance_angle::Float64, exit_angle::Float64, FringeBendEntrance::Int, FringeBendExit::Int,
    fint1::Float64, fint2::Float64, gap::Float64, FringeQuadEntrance::Int, FringeQuadExit::Int,
    fringeIntM0::Array{Float64,1}, fringeIntP0::Array{Float64,1}, T1::Array{Float64,1}, T2::Array{Float64,1}, 
    R1::Array{Float64,2}, R2::Array{Float64,2}, RApertures::Array{Float64,1}, EApertures::Array{Float64,1},
    KickAngle::Array{Float64,1}, E0::Float64, num_particles::Int, lost_flags::Array{Int64,1})
    
    DRIFT1 = 0.6756035959798286638
    DRIFT2 = -0.1756035959798286639
    KICK1 = 1.351207191959657328
    KICK2 = -1.702414383919314656
    SL = le / num_int_steps
    L1 = SL * DRIFT1
    L2 = SL * DRIFT2
    K1 = SL * KICK1
    K2 = SL * KICK2

    if FringeQuadEntrance==2 #&& !isnothing(fringeIntM0) && !isnothing(fringeIntP0)
        useLinFrEleEntrance = 1
    else
        useLinFrEleEntrance = 0
    end
    if FringeQuadExit==2 #&& !isnothing(fringeIntM0) && !isnothing(fringeIntP0)
        useLinFrEleExit = 1
    else
        useLinFrEleExit = 0
    end

    B[1] -= sin(KickAngle[1]) / le
    A[1] += sin(KickAngle[2]) / le

    Threads.@threads for c in 1:num_particles
    # for c in 1:num_particles
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
        edge_fringe_entrance!(r6, irho, entrance_angle, fint1, gap, FringeBendEntrance)

        # Quadrupole gradient fringe entrance
        if !iszero(FringeQuadEntrance) && !iszero(B[2])
            if isone(useLinFrEleEntrance)
                linearQuadFringeElegantEntrance!(r6, B[2], fringeIntM0, fringeIntP0)
            else
                QuadFringePassP!(r6, B[2])
            end
        end

        # Integrator
        for m in 1:num_int_steps
            drift6!(r6, L1, beti)
            bndthinkickrad!(r6, A, B, K1, irho, E0, max_order, beti)
            drift6!(r6, L2, beti)
            bndthinkickrad!(r6, A, B, K2, irho, E0, max_order, beti)
            drift6!(r6, L2, beti)
            bndthinkickrad!(r6, A, B, K1, irho, E0, max_order, beti)
            drift6!(r6, L1, beti)
        end

        # Quadrupole gradient fringe exit
        if !iszero(FringeQuadExit) && !iszero(B[2])
            if isone(useLinFrEleExit)
                linearQuadFringeElegantExit!(r6, B[2], fringeIntM0, fringeIntP0)
            else
                QuadFringePassN!(r6, B[2])
            end
        end

        # Edge focus at exit
        edge_fringe_exit!(r6, irho, exit_angle, fint2, gap, FringeBendExit)

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

    B[1] += sin(KickAngle[1]) / le
    A[1] -= sin(KickAngle[2]) / le
    return nothing
end

function BendSymplecticPass_P!(r::Array{Float64,1}, le::Float64, beti::Float64, irho::Float64, A::Array{Float64,1}, B::Array{Float64,1}, 
    max_order::Int, num_int_steps::Int, entrance_angle::Float64, exit_angle::Float64, FringeBendEntrance::Int, FringeBendExit::Int,
    fint1::Float64, fint2::Float64, gap::Float64, FringeQuadEntrance::Int, FringeQuadExit::Int,
    fringeIntM0::Array{Float64,1}, fringeIntP0::Array{Float64,1}, T1::Array{Float64,1}, T2::Array{Float64,1}, 
    R1::Array{Float64,2}, R2::Array{Float64,2}, RApertures::Array{Float64,1}, EApertures::Array{Float64,1},
    KickAngle::Array{Float64,1}, num_particles::Int, lost_flags::Array{Int64,1})
    
    DRIFT1 = 0.6756035959798286638
    DRIFT2 = -0.1756035959798286639
    KICK1 = 1.351207191959657328
    KICK2 = -1.702414383919314656
    SL = le / num_int_steps
    L1 = SL * DRIFT1
    L2 = SL * DRIFT2
    K1 = SL * KICK1
    K2 = SL * KICK2

    if FringeQuadEntrance==2#&& !isnothing(fringeIntM0) && !isnothing(fringeIntP0)
        useLinFrEleEntrance = 1
    else
        useLinFrEleEntrance = 0
    end
    if FringeQuadExit==2# && !isnothing(fringeIntM0) && !isnothing(fringeIntP0)
        useLinFrEleExit = 1
    else
        useLinFrEleExit = 0
    end

    B[1] -= sin(KickAngle[1]) / le
    A[1] += sin(KickAngle[2]) / le


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
        edge_fringe_entrance!(r6, irho, entrance_angle, fint1, gap, FringeBendEntrance)

        # Quadrupole gradient fringe entrance
        if !iszero(FringeQuadEntrance) && !iszero(B[2])
            if isone(useLinFrEleEntrance)
                linearQuadFringeElegantEntrance!(r6, B[2], fringeIntM0, fringeIntP0)
            else
                QuadFringePassP!(r6, B[2])
            end
        end

        # Integrator
        for m in 1:num_int_steps
            drift6!(r6, L1, beti)
            bndthinkick!(r6, A, B, K1, irho, max_order, beti)
            drift6!(r6, L2, beti)
            bndthinkick!(r6, A, B, K2, irho, max_order, beti)
            drift6!(r6, L2, beti)
            bndthinkick!(r6, A, B, K1, irho, max_order, beti)
            drift6!(r6, L1, beti)
        end

        # Quadrupole gradient fringe exit
        if !iszero(FringeQuadExit) && !iszero(B[2])
            if isone(useLinFrEleExit)
                linearQuadFringeElegantExit!(r6, B[2], fringeIntM0, fringeIntP0)
            else
                QuadFringePassN!(r6, B[2])
            end
        end

        # Edge focus at exit
        edge_fringe_exit!(r6, irho, exit_angle, fint2, gap, FringeBendExit)

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
    
    B[1] += sin(KickAngle[1]) / le
    A[1] -= sin(KickAngle[2]) / le
    return nothing
end


function pass_P!(ele::SBEND, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam{Float64})
    # ele: SBEND
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    lost_flags = particles.lost_flag
    irho = ele.angle / ele.len
    E0 = particles.energy
    if use_exact_beti == 1
        beti = 1.0 / particles.beta
    else
        beti = 1.0 
    end
    if ele.rad == 0
        BendSymplecticPass_P!(r_in, ele.len, irho, beti, ele.PolynomA, ele.PolynomB, ele.MaxOrder, ele.NumIntSteps,
            ele.e1, ele.e2,
            ele.FringeBendEntrance, ele.FringeBendExit,
            ele.fint1, ele.fint2, ele.gap,
            ele.FringeQuadEntrance, ele.FringeQuadExit,
            ele.FringeIntM0, ele.FringeIntP0,
            ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures,
            ele.KickAngle, num_particles, lost_flags)
    else
        BendSymplecticPassRad_P!(r_in, ele.len, beti, irho, ele.PolynomA, ele.PolynomB, ele.MaxOrder, ele.NumIntSteps,
            ele.e1, ele.e2,
            ele.FringeBendEntrance, ele.FringeBendExit,
            ele.fint1, ele.fint2, ele.gap,
            ele.FringeQuadEntrance, ele.FringeQuadExit,
            ele.FringeIntM0, ele.FringeIntP0,
            ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures,
            ele.KickAngle, E0, num_particles, lost_flags)
    end
    return nothing
end



# Tracking functions for exact sector bend
function pxyz(dp1, px, py)
    return sqrt(dp1*dp1 - px*px - py*py)
end

function Yrot!(r::AbstractVector{Float64}, phi::Float64, beti::Float64)
    # Forest 10.26, rotation in free space
    if phi != 0.0
        dp1 = 1.0*beti + r[6]
        c = cos(phi)
        s = sin(phi)
        pz = pxyz(dp1, r[2], r[4])
        p = c * pz - s * r[2]
        px = s * pz + c * r[2]
        x = r[1] * pz / p
        dy = r[1] * r[4] * s / p
        dct = dp1 * r[1] * s / p
        r[1] = x
        r[2] = px
        r[3] += dy
        r[5] += dct
    end
    return nothing
end

function Sec(x::Float64)
    return 1.0 / cos(x)
end

function bend_fringe!(r::AbstractVector{Float64}, irho::Float64, gK::Float64, beti::Float64)
    # Forest 13.13, bend fringe in the hard-edge limit
    b0 = irho 
    pz = pxyz(1.0*beti+r[6], r[2], r[4])
    px = r[2]
    py = r[4]
    d = r[6]
    xp = px / pz
    yp = py / pz
    phi = -b0 * tan( b0 * gK * (1.0 + xp*xp*(2.0 + yp*yp))*pz - atan(xp / (1.0 + yp*yp)))
    px2 = px*px
    px4 = px2*px2
    py2 = py*py
    py4 = py2*py2
    py6 = py4*py2
    pz2 = pz*pz
    pz3 = pz2*pz
    pz4 = pz2*pz2
    pz5 = pz4*pz
    pz6 = pz4*pz2
    py2z2 = (py2 + pz2) * (py2 + pz2)
    # powsec = pow(Sec((b0*gK*(pz4 + px2*(py2 + 2*pz2)))/pz3 - atan((px*pz)/(py2 + pz2))),2)
    powsec = Sec((b0*gK*(pz4 + px2*(py2 + 2.0*pz2)))/pz3 - atan((px*pz)/(py2 + pz2)))^2
    denom = (pz5*(py4 + px2*pz2 + 2.0*py2*pz2 + pz4))

    dpx = -(b0*(px2*pz4*(py2 - pz2) - pz6*(py2 + pz2) +
            b0*gK*px*(pz2*py2z2*(2.0*py2 + 3.0*pz2) + px4*(3.0*py2*pz2 + 2.0*pz4) +
            px2*(3.0*py6 + 8.0*py4*pz2 + 9.0*py2*pz4 + 5.0*pz6)))*powsec)/denom
    dpy = -(b0*py*(px*pz4*(py2 + pz2) +
            b0*gK*(-(pz4*py2z2) + px4*(3.0*py2*pz2 + 4.0*pz4) +
                   px2*(3.0*py6 + 10.0*py4*pz2 + 11.0*py2*pz4 + 3.0*pz6)))*powsec)/denom
    dd = (b0*(1.0*beti + d)*(px*pz4*(py2 - pz2) + b0*gK*
                      (-(pz4*py2z2) + px4*(3.0*py2*pz2 + 2.0*pz4) +
                       px2*(3.0*py6 + 8.0*py4*pz2 + 7.0*py2*pz4 + pz6)))*powsec)/denom

    yf = (2.0 * r[3]) / (1.0 + sqrt(1.0 - 2.0 * dpy * r[3]))
    dxf = 0.5 * dpx * yf * yf
    dct = 0.5 * dd * yf * yf
    dpyf = phi * yf

    r[3] = yf
    r[1] += dxf
    r[4] -= dpyf
    r[5] -= dct
    return nothing
end

function multipole_fringe!(r6::AbstractVector{Float64}, le::Float64, polya::Array{Float64,1}, polyb::Array{Float64,1}, 
        max_order::Int, edge::Float64, skip_b0::Int, beti::Float64)
    # Forest 13.29
    FX = 0.0
    FY = 0.0
    FX_X = 0.0
    FX_Y = 0.0
    FY_X = 0.0
    FY_Y = 0.0
  
    RX = 1.0
    IX = 0.0

    for n in 0:max_order
        B = polyb[n + 1]  
        A = polya[n + 1] 
    
        j = n + 1.0
    
        DRX = RX
        DIX = IX
    
        # Complex multiplications
        RX = DRX * r6[1] - DIX * r6[3]
        IX = DRX * r6[3] + DIX * r6[1]
        
        U, V, DU, DV = 0.0, 0.0, 0.0, 0.0
        if n == 0 && skip_b0 != 0
            U -= A * IX
            V += A * RX
            DU -= A * DIX
            DV += A * DRX
        else
            U += B * RX - A * IX
            V += B * IX + A * RX
            DU += B * DRX - A * DIX
            DV += B * DIX + A * DRX
        end
    
        f1 = -edge / 4.0 / (j + 1.0)
    
        U *= f1
        V *= f1
        DU *= f1
        DV *= f1
    
        DUX = j * DU
        DVX = j * DV
        DUY = -j * DV
        DVY = j * DU
    
        nf = 1.0 * (j + 2.0) / j
    
        FX += U * r6[1] + nf * V * r6[3]
        FY += U * r6[3] - nf * V * r6[1]
    
        FX_X += DUX * r6[1] + U + nf * r6[3] * DVX
        FX_Y += DUY * r6[1] + nf * V + nf * r6[3] * DVY
    
        FY_X += DUX * r6[3] - nf * V - nf * r6[1] * DVX
        FY_Y += DUY * r6[3] + U - nf * r6[1] * DVY
    end

    DEL = 1.0 / (1.0*beti + r6[6])
    A = 1.0 - FX_X * DEL
    B = -FY_X * DEL
    D = 1.0 - FY_Y * DEL
    C = -FX_Y * DEL

    r6[1] -= FX * DEL
    r6[3] -= FY * DEL

    pxf = (D * r6[2] - B * r6[4]) / (A * D - B * C)
    pyf = (A * r6[4] - C * r6[2]) / (A * D - B * C)
    r6[4] = pyf
    r6[2] = pxf
    r6[5] -= (r6[2] * FX + r6[4] * FY) * DEL * DEL
    return nothing
end

function exact_bend!(r6::AbstractVector{Float64}, irho::Float64, L::Float64, beti::Float64)
    # Forest 12.18, bend-kick split, map W(L, irho)

    dp1 = 1.0*beti + r6[6]  # r6[delta_]
    pz = pxyz(dp1, r6[2], r6[4])  # r6[px_], r6[py_]

    if abs(irho) < 1e-6
        NormL = L / pz
        r6[1] += r6[2] * NormL  # r6[x_]
        r6[3] += r6[4] * NormL  # r6[y_]
        r6[5] += NormL * dp1    # r6[ct_], absolute path length
    else
        pzmx = pz - (1.0 + r6[1] * irho)  # r6[x_]
        cs = cos(irho * L)
        sn = sin(irho * L)
        px = r6[2] * cs + pzmx * sn  # r6[px_]
        d2 = pxyz(dp1, 0.0, r6[4])  # r6[py_]
        dasin = L + (asin(r6[2] / d2) - asin(px / d2)) / irho
        x = (pxyz(dp1, px, r6[4]) - pzmx * cs + r6[2] * sn - 1.0) / irho  # r6[x_]
        dy = r6[4] * dasin  # r6[py_]
        dct = dp1 * dasin   # r6[ct_], absolute path length

        r6[1] = x
        r6[2] = px
        r6[3] += dy
        r6[5] += dct
    end
    return nothing
end

function B2perp_exact_bnd(bx, by, irho, x, xpr, y, ypr)
    nrm = (1.0 + x *irho)^2
    v_norm2 = nrm + xpr^2*(1.0-nrm) + ypr^2*(1.0-nrm)
    return bx*bx + by*by - (bx*xpr + by*ypr)^2/v_norm2
end
function exactbndthinkick_rad!(r::AbstractVector{Float64}, A::Array{Float64,1}, B::Array{Float64,1}, 
    L::Float64, irho::Float64, max_order::Int, beti::Float64, rad_const::Float64)
    # AT function. Ref[Terebilo, Andrei. "Accelerator modeling with MATLAB accelerator toolbox." PACS2001 (2001)].
    ReSum = B[max_order + 1]
    ImSum = A[max_order + 1]
    ReSumTemp = 0.0

    for i in max_order-1:-1:0
        ReSumTemp = ReSum * r[1] - ImSum * r[3] + B[i+1]
        ImSum = ImSum * r[1] + ReSum * r[3] + A[i+1]
        ReSum = ReSumTemp
    end

    p_norm = 1.0 / (1.0 + r[6])
    x = r[1]
    xpr = r[2] * p_norm
    y = r[3]
    ypr = r[4] * p_norm

    B2P = B2perp_exact_bnd(ImSum, ReSum+irho, irho, x, xpr, y, ypr)

    r[6] -= rad_const * (1.0 + r[6])^2 * B2P * (1.0 + x*irho) * L / sqrt(1.0 - xpr*xpr - ypr*ypr)

    p_norm = 1.0 / (1.0 + r[6])
    r[2] = xpr/p_norm
    r[4] = ypr/p_norm

    r[2] -= L * ReSum
    r[4] += L * ImSum
    return nothing
end

function bend_edge!(r6::AbstractVector{Float64}, rhoinv::Float64, theta::Float64, beti::Float64)
    # Forest 12.41, ideal wedge, map U(theta, rhoinv)

    if abs(rhoinv) >= 1e-6
        dp1 = 1.0*beti + r6[6]  # r6[delta_]
        c = cos(theta)
        s = sin(theta)
        pz = pxyz(dp1, r6[2], r6[4])  # r6[px_], r6[py_]
        d2 = pxyz(dp1, 0.0, r6[4])    # r6[py_]
        px = r6[2] * c + (pz - rhoinv * r6[1]) * s  # r6[px_]
        dasin = asin(r6[2] / d2) - asin(px / d2)
        num = r6[1] * (r6[2] * sin(2.0 * theta) + s * s * (2.0 * pz - rhoinv * r6[1]))
        den = pxyz(dp1, px, r6[4]) + pxyz(dp1, r6[2], r6[4]) * c - r6[2] * s
        x = r6[1] * c + num / den  # r6[x_]
        dy = r6[4] * theta / rhoinv + r6[4] / rhoinv * dasin  # r6[py_]
        dct = dp1 / rhoinv * (theta + dasin)  # r6[ct_]

        r6[1] = x
        r6[2] = px
        r6[3] += dy
        r6[5] += dct
    end
    return nothing
end

function ExactSectorBend!(r::Array{Float64,1}, le::Float64, beti::Float64, angle::Float64, A::Array{Float64,1}, B::Array{Float64,1}, 
    max_order::Int, num_int_steps::Int, entrance_angle::Float64, exit_angle::Float64, FringeBendEntrance::Int, FringeBendExit::Int,
    FringeQuadEntrance::Int, FringeQuadExit::Int, gk::Float64,
    T1::Array{Float64,1}, T2::Array{Float64,1}, 
    R1::Array{Float64,2}, R2::Array{Float64,2}, RApertures::Array{Float64,1}, EApertures::Array{Float64,1},
    KickAngle::Array{Float64,1}, num_particles::Int, lost_flags::Array{Int64,1})
    
    irho = angle / le
    DRIFT1 = 0.6756035959798286638
    DRIFT2 = -0.1756035959798286639
    KICK1 = 1.351207191959657328
    KICK2 = -1.702414383919314656
    SL = le / num_int_steps
    L1 = SL * DRIFT1
    L2 = SL * DRIFT2
    K1 = SL * KICK1
    K2 = SL * KICK2

    B0 = B[1]
    A0 = A[1]

    if !iszero(KickAngle[1])
        B[1] -= sin(KickAngle[1]) / le
    end
    if !iszero(KickAngle[2])
        A[1] += sin(KickAngle[2]) / le
    end

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

        Yrot!(r6, entrance_angle, beti)

        if FringeBendEntrance != 0
            bend_fringe!(r6, irho, gk, beti)
        end

        if FringeQuadEntrance != 0
            multipole_fringe!(r6, le, A, B, max_order, 1.0, 1, beti)
        end

        bend_edge!(r6, irho, -entrance_angle, beti)

        # Integrator
        if num_int_steps == 0
            exact_bend!(r6, irho, le, beti)
        else
            for m in 1:num_int_steps
                exact_bend!(r6, irho, L1, beti)
                strthinkick!(r6, A, B, K1, max_order)
                exact_bend!(r6, irho, L2, beti)
                strthinkick!(r6, A, B, K2, max_order)
                exact_bend!(r6, irho, L2, beti)
                strthinkick!(r6, A, B, K1, max_order)
                exact_bend!(r6, irho, L1, beti)
            end
        end

        r6[5] -= le*beti  # r6[ct_], absolute path length

        bend_edge!(r6, irho, -exit_angle, beti)
        if FringeQuadExit != 0
            multipole_fringe!(r6, le, A, B, max_order, -1.0, 1, beti)
        end
        if FringeBendExit != 0
            bend_fringe!(r6, -irho, gk, beti)
        end
        Yrot!(r6, exit_angle, beti)

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
    
    if !iszero(KickAngle[1])
        B[1] += sin(KickAngle[1]) / le
    end
    if !iszero(KickAngle[2])
        A[1] -= sin(KickAngle[2]) / le
    end
    return nothing
end

function ExactSectorBend_rad!(r::Array{Float64,1}, le::Float64, rad_const::Float64, beti::Float64, angle::Float64, A::Array{Float64,1}, B::Array{Float64,1}, 
    max_order::Int, num_int_steps::Int, entrance_angle::Float64, exit_angle::Float64, FringeBendEntrance::Int, FringeBendExit::Int,
    FringeQuadEntrance::Int, FringeQuadExit::Int, gk::Float64,
    T1::Array{Float64,1}, T2::Array{Float64,1}, 
    R1::Array{Float64,2}, R2::Array{Float64,2}, RApertures::Array{Float64,1}, EApertures::Array{Float64,1},
    KickAngle::Array{Float64,1}, num_particles::Int, lost_flags::Array{Int64,1})
    
    irho = angle / le
    DRIFT1 = 0.6756035959798286638
    DRIFT2 = -0.1756035959798286639
    KICK1 = 1.351207191959657328
    KICK2 = -1.702414383919314656
    SL = le / num_int_steps
    L1 = SL * DRIFT1
    L2 = SL * DRIFT2
    K1 = SL * KICK1
    K2 = SL * KICK2

    B0 = B[1]
    A0 = A[1]

    if !iszero(KickAngle[1])
        B[1] -= sin(KickAngle[1]) / le
    end
    if !iszero(KickAngle[2])
        A[1] += sin(KickAngle[2]) / le
    end

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

        Yrot!(r6, entrance_angle, beti)

        if FringeBendEntrance != 0
            bend_fringe!(r6, irho, gk, beti)
        end

        if FringeQuadEntrance != 0
            multipole_fringe!(r6, le, A, B, max_order, 1.0, 1, beti)
        end

        bend_edge!(r6, irho, -entrance_angle, beti)

        # Integrator
        if num_int_steps == 0
            exact_bend!(r6, irho, le, beti)
        else
            for m in 1:num_int_steps
                exact_bend!(r6, irho, L1, beti)
                exactbndthinkick_rad!(r6, A, B, K1, irho, max_order, beti, rad_const)
                exact_bend!(r6, irho, L2, beti)
                exactbndthinkick_rad!(r6, A, B, K2, irho, max_order, beti, rad_const)
                exact_bend!(r6, irho, L2, beti)
                exactbndthinkick_rad!(r6, A, B, K1, irho, max_order, beti, rad_const)
                exact_bend!(r6, irho, L1, beti)
            end
        end

        r6[5] -= le*beti  # r6[ct_], absolute path length

        bend_edge!(r6, irho, -exit_angle, beti)
        if FringeQuadExit != 0
            multipole_fringe!(r6, le, A, B, max_order, -1.0, 1, beti)
        end
        if FringeBendExit != 0
            bend_fringe!(r6, -irho, gk, beti)
        end
        Yrot!(r6, exit_angle, beti)

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
    
    if !iszero(KickAngle[1])
        B[1] += sin(KickAngle[1]) / le
    end
    if !iszero(KickAngle[2])
        A[1] -= sin(KickAngle[2]) / le
    end
    return nothing
end

function ExactSectorBend_rad_P!(r::Array{Float64,1}, le::Float64, rad_const::Float64, beti::Float64, angle::Float64, A::Array{Float64,1}, B::Array{Float64,1}, 
    max_order::Int, num_int_steps::Int, entrance_angle::Float64, exit_angle::Float64, FringeBendEntrance::Int, FringeBendExit::Int,
    FringeQuadEntrance::Int, FringeQuadExit::Int, gk::Float64,
    T1::Array{Float64,1}, T2::Array{Float64,1}, 
    R1::Array{Float64,2}, R2::Array{Float64,2}, RApertures::Array{Float64,1}, EApertures::Array{Float64,1},
    KickAngle::Array{Float64,1}, num_particles::Int, lost_flags::Array{Int64,1})
    
    irho = angle / le
    DRIFT1 = 0.6756035959798286638
    DRIFT2 = -0.1756035959798286639
    KICK1 = 1.351207191959657328
    KICK2 = -1.702414383919314656
    SL = le / num_int_steps
    L1 = SL * DRIFT1
    L2 = SL * DRIFT2
    K1 = SL * KICK1
    K2 = SL * KICK2

    B0 = B[1]
    A0 = A[1]

    if !iszero(KickAngle[1])
        B[1] -= sin(KickAngle[1]) / le
    end
    if !iszero(KickAngle[2])
        A[1] += sin(KickAngle[2]) / le
    end

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

        Yrot!(r6, entrance_angle, beti)

        if FringeBendEntrance != 0
            bend_fringe!(r6, irho, gk, beti)
        end

        if FringeQuadEntrance != 0
            multipole_fringe!(r6, le, A, B, max_order, 1.0, 1, beti)
        end

        bend_edge!(r6, irho, -entrance_angle, beti)

        # Integrator
        if num_int_steps == 0
            exact_bend!(r6, irho, le, beti)
        else
            for m in 1:num_int_steps
                exact_bend!(r6, irho, L1, beti)
                exactbndthinkick_rad!(r6, A, B, K1, irho, max_order, beti, rad_const)
                exact_bend!(r6, irho, L2, beti)
                exactbndthinkick_rad!(r6, A, B, K2, irho, max_order, beti, rad_const)
                exact_bend!(r6, irho, L2, beti)
                exactbndthinkick_rad!(r6, A, B, K1, irho, max_order, beti, rad_const)
                exact_bend!(r6, irho, L1, beti)
            end
        end

        r6[5] -= le*beti  # r6[ct_], absolute path length

        bend_edge!(r6, irho, -exit_angle, beti)
        if FringeQuadExit != 0
            multipole_fringe!(r6, le, A, B, max_order, -1.0, 1, beti)
        end
        if FringeBendExit != 0
            bend_fringe!(r6, -irho, gk, beti)
        end
        Yrot!(r6, exit_angle, beti)

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
    
    if !iszero(KickAngle[1])
        B[1] += sin(KickAngle[1]) / le
    end
    if !iszero(KickAngle[2])
        A[1] -= sin(KickAngle[2]) / le
    end
    return nothing
end

function ExactSectorBend_P!(r::Array{Float64,1}, le::Float64, beti::Float64, angle::Float64, A::Array{Float64,1}, B::Array{Float64,1}, 
    max_order::Int, num_int_steps::Int, entrance_angle::Float64, exit_angle::Float64, FringeBendEntrance::Int, FringeBendExit::Int,
    FringeQuadEntrance::Int, FringeQuadExit::Int, gk::Float64,
    T1::Array{Float64,1}, T2::Array{Float64,1}, 
    R1::Array{Float64,2}, R2::Array{Float64,2}, RApertures::Array{Float64,1}, EApertures::Array{Float64,1},
    KickAngle::Array{Float64,1}, num_particles::Int, lost_flags::Array{Int64,1})
    
    irho = angle / le
    DRIFT1 = 0.6756035959798286638
    DRIFT2 = -0.1756035959798286639
    KICK1 = 1.351207191959657328
    KICK2 = -1.702414383919314656
    SL = le / num_int_steps
    L1 = SL * DRIFT1
    L2 = SL * DRIFT2
    K1 = SL * KICK1
    K2 = SL * KICK2

    B0 = B[1]
    A0 = A[1]

    if !iszero(KickAngle[1])
        B[1] -= sin(KickAngle[1]) / le
    end
    if !iszero(KickAngle[2])
        A[1] += sin(KickAngle[2]) / le
    end

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

        Yrot!(r6, entrance_angle, beti)

        if FringeBendEntrance != 0
            bend_fringe!(r6, irho, gk, beti)
        end

        if FringeQuadEntrance != 0
            multipole_fringe!(r6, le, A, B, max_order, 1.0, 1, beti)
        end

        bend_edge!(r6, irho, -entrance_angle, beti)

        # Integrator
        if num_int_steps == 0
            exact_bend!(r6, irho, le, beti)
        else
            for m in 1:num_int_steps
                exact_bend!(r6, irho, L1, beti)
                strthinkick!(r6, A, B, K1, max_order)
                exact_bend!(r6, irho, L2, beti)
                strthinkick!(r6, A, B, K2, max_order)
                exact_bend!(r6, irho, L2, beti)
                strthinkick!(r6, A, B, K1, max_order)
                exact_bend!(r6, irho, L1, beti)
            end
        end

        r6[5] -= le*beti  # r6[ct_], absolute path length

        bend_edge!(r6, irho, -exit_angle, beti)
        if FringeQuadExit != 0
            multipole_fringe!(r6, le, A, B, max_order, -1.0, 1, beti)
        end
        if FringeBendExit != 0
            bend_fringe!(r6, -irho, gk, beti)
        end
        Yrot!(r6, exit_angle, beti)

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
    
    if !iszero(KickAngle[1])
        B[1] += sin(KickAngle[1]) / le
    end
    if !iszero(KickAngle[2])
        A[1] -= sin(KickAngle[2]) / le
    end
    return nothing
end

function pass!(ele::ESBEND, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam{Float64})
    # ele: ESBEND
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    lost_flags = particles.lost_flag
    E0 = particles.energy
    rad_const = 0.0
    if use_exact_beti == 1
        beti = 1.0 / particles.beta
    else
        beti = 1.0 
    end
    if ele.rad == 0
        ExactSectorBend!(r_in, ele.len, beti, ele.angle, ele.PolynomA, ele.PolynomB, ele.MaxOrder, ele.NumIntSteps,
            ele.e1, ele.e2,
            ele.FringeBendEntrance, ele.FringeBendExit,
            ele.FringeQuadEntrance, ele.FringeQuadExit,
            ele.gK,
            ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures,
            ele.KickAngle, num_particles, lost_flags)
    else
        if particles.mass == m_e
            rad_const = RAD_CONST_E * particles.gamma^3
        elseif particles.mass == m_p
            rad_const = RAD_CONST_P * particles.gamma^3
        else
            rad_const = 0.0
            println("SR is not implemented for this particle mass.")
        end
        ExactSectorBend_rad!(r_in, ele.len, rad_const, beti, ele.angle, ele.PolynomA, ele.PolynomB, ele.MaxOrder, ele.NumIntSteps,
            ele.e1, ele.e2,
            ele.FringeBendEntrance, ele.FringeBendExit,
            ele.FringeQuadEntrance, ele.FringeQuadExit,
            ele.gK,
            ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures,
            ele.KickAngle, num_particles, lost_flags)
    end
    return nothing
end

function pass_P!(ele::ESBEND, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam{Float64})
    # ele: ESBEND
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    lost_flags = particles.lost_flag
    E0 = particles.energy
    rad_const = 0.0
    if use_exact_beti == 1
        beti = 1.0 / particles.beta
    else
        beti = 1.0 
    end
    if ele.rad == 0
        ExactSectorBend_P!(r_in, ele.len, beti, ele.angle, ele.PolynomA, ele.PolynomB, ele.MaxOrder, ele.NumIntSteps,
            ele.e1, ele.e2,
            ele.FringeBendEntrance, ele.FringeBendExit,
            ele.FringeQuadEntrance, ele.FringeQuadExit,
            ele.gK,
            ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures,
            ele.KickAngle, num_particles, lost_flags)
    else
        if particles.mass == m_e
            rad_const = RAD_CONST_E * particles.gamma^3
        elseif particles.mass == m_p
            rad_const = RAD_CONST_P * particles.gamma^3
        else
            rad_const = 0.0
            println("SR is not implemented for this particle mass.")
        end
        ExactSectorBend_rad_P!(r_in, ele.len, rad_const, beti, ele.angle, ele.PolynomA, ele.PolynomB, ele.MaxOrder, ele.NumIntSteps,
            ele.e1, ele.e2,
            ele.FringeBendEntrance, ele.FringeBendExit,
            ele.FringeQuadEntrance, ele.FringeQuadExit,
            ele.gK,
            ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures,
            ele.KickAngle, num_particles, lost_flags)
    end
    return nothing
end

# TPSA tracking function for exact sector bend
function Yrot!(r::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}, phi::Float64, beti::Float64) where {T, TPS_Dim, Max_TPS_Degree}
    # Forest 10.26, rotation in free space
    if phi != 0.0
        dp1 = 1.0*beti + r[6]
        c = cos(phi)
        s = sin(phi)
        pz = pxyz(dp1, r[2], r[4])
        p = c * pz - s * r[2]
        px = s * pz + c * r[2]
        x = r[1] * pz / p
        dy = r[1] * r[4] * s / p
        dct = dp1 * r[1] * s / p
        r[1] = x
        r[2] = px
        r[3] += dy
        r[5] += dct
    end
    return nothing
end

function Sec(x::CTPS{T, TPS_Dim, Max_TPS_Degree}) where {T, TPS_Dim, Max_TPS_Degree}
    return 1.0 / cos(x)
end

function bend_fringe!(r::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}, irho::Float64, gK::Float64, beti::Float64) where {T, TPS_Dim, Max_TPS_Degree}
    # Forest 13.13, bend fringe in the hard-edge limit
    b0 = irho 
    pz = pxyz(1.0*beti+r[6], r[2], r[4])
    px = r[2]
    py = r[4]
    d = r[6]
    xp = px / pz
    yp = py / pz
    phi = -b0 * tan( b0 * gK * (1.0 + xp*xp*(2.0 + yp*yp))*pz - atan(xp / (1.0 + yp*yp)))
    px2 = px*px
    px4 = px2*px2
    py2 = py*py
    py4 = py2*py2
    py6 = py4*py2
    pz2 = pz*pz
    pz3 = pz2*pz
    pz4 = pz2*pz2
    pz5 = pz4*pz
    pz6 = pz4*pz2
    py2z2 = (py2 + pz2) * (py2 + pz2)
    # powsec = pow(Sec((b0*gK*(pz4 + px2*(py2 + 2*pz2)))/pz3 - atan((px*pz)/(py2 + pz2))),2)
    powsec = Sec((b0*gK*(pz4 + px2*(py2 + 2.0*pz2)))/pz3 - atan((px*pz)/(py2 + pz2)))^2
    denom = (pz5*(py4 + px2*pz2 + 2.0*py2*pz2 + pz4))

    dpx = -(b0*(px2*pz4*(py2 - pz2) - pz6*(py2 + pz2) +
            b0*gK*px*(pz2*py2z2*(2.0*py2 + 3.0*pz2) + px4*(3.0*py2*pz2 + 2.0*pz4) +
            px2*(3.0*py6 + 8.0*py4*pz2 + 9.0*py2*pz4 + 5.0*pz6)))*powsec)/denom
    dpy = -(b0*py*(px*pz4*(py2 + pz2) +
            b0*gK*(-(pz4*py2z2) + px4*(3.0*py2*pz2 + 4.0*pz4) +
                   px2*(3.0*py6 + 10.0*py4*pz2 + 11.0*py2*pz4 + 3.0*pz6)))*powsec)/denom
    dd = (b0*(1.0*beti + d)*(px*pz4*(py2 - pz2) + b0*gK*
                      (-(pz4*py2z2) + px4*(3.0*py2*pz2 + 2.0*pz4) +
                       px2*(3.0*py6 + 8.0*py4*pz2 + 7.0*py2*pz4 + pz6)))*powsec)/denom

    yf = (2.0 * r[3]) / (1.0 + sqrt(1.0 - 2.0 * dpy * r[3]))
    dxf = 0.5 * dpx * yf * yf
    dct = 0.5 * dd * yf * yf
    dpyf = phi * yf

    r[3] = yf
    r[1] += dxf
    r[4] -= dpyf
    r[5] -= dct
    return nothing
end

function multipole_fringe!(r6::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}, le::Float64, polya::Array{Float64,1}, polyb::Array{Float64,1}, 
    max_order::Int, edge::Float64, skip_b0::Int, beti::Float64) where {T, TPS_Dim, Max_TPS_Degree}
    # Forest 13.29
    FX = CTPS(T(0.0), TPS_Dim, Max_TPS_Degree)
    FY = CTPS(T(0.0), TPS_Dim, Max_TPS_Degree)
    FX_X = CTPS(T(0.0), TPS_Dim, Max_TPS_Degree)
    FX_Y = CTPS(T(0.0), TPS_Dim, Max_TPS_Degree)
    FY_X = CTPS(T(0.0), TPS_Dim, Max_TPS_Degree)
    FY_Y = CTPS(T(0.0), TPS_Dim, Max_TPS_Degree)
  
    RX = CTPS(T(1.0), TPS_Dim, Max_TPS_Degree)
    IX = CTPS(T(0.0), TPS_Dim, Max_TPS_Degree)

    for n in 0:max_order
        B = polyb[n + 1]  
        A = polya[n + 1] 
    
        j = n + 1.0
    
        DRX = CTPS(RX)
        DIX = CTPS(IX)
    
        # Complex multiplications
        RX = DRX * r6[1] - DIX * r6[3]
        IX = DRX * r6[3] + DIX * r6[1]
        
        U, V, DU, DV = CTPS(T(0.0), TPS_Dim, Max_TPS_Degree), CTPS(T(0.0), TPS_Dim, Max_TPS_Degree), CTPS(T(0.0), TPS_Dim, Max_TPS_Degree), CTPS(T(0.0), TPS_Dim, Max_TPS_Degree)
        if n == 0 && skip_b0 != 0
            U -= A * IX
            V += A * RX
            DU -= A * DIX
            DV += A * DRX
        else
            U += B * RX - A * IX
            V += B * IX + A * RX
            DU += B * DRX - A * DIX
            DV += B * DIX + A * DRX
        end
    
        f1 = -edge / 4.0 / (j + 1.0)
    
        U *= f1
        V *= f1
        DU *= f1
        DV *= f1
    
        DUX = j * DU
        DVX = j * DV
        DUY = -j * DV
        DVY = j * DU
    
        nf = 1.0 * (j + 2.0) / j
    
        FX += U * r6[1] + nf * V * r6[3]
        FY += U * r6[3] - nf * V * r6[1]
    
        FX_X += DUX * r6[1] + U + nf * r6[3] * DVX
        FX_Y += DUY * r6[1] + nf * V + nf * r6[3] * DVY
    
        FY_X += DUX * r6[3] - nf * V - nf * r6[1] * DVX
        FY_Y += DUY * r6[3] + U - nf * r6[1] * DVY
    end

    DEL = 1.0 / (1.0*beti + r6[6])
    A = 1.0 - FX_X * DEL
    B = -FY_X * DEL
    D = 1.0 - FY_Y * DEL
    C = -FX_Y * DEL

    r6[1] -= FX * DEL
    r6[3] -= FY * DEL

    pxf = (D * r6[2] - B * r6[4]) / (A * D - B * C)
    pyf = (A * r6[4] - C * r6[2]) / (A * D - B * C)
    r6[4] = pyf
    r6[2] = pxf
    r6[5] -= (r6[2] * FX + r6[4] * FY) * DEL * DEL
    return nothing
end

function exact_bend!(r6::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}, irho::Float64, L::Float64, beti::Float64) where {T, TPS_Dim, Max_TPS_Degree}
    # Forest 12.18, bend-kick split, map W(L, irho)

    dp1 = 1.0*beti + r6[6]  # r6[delta_]
    pz = pxyz(dp1, r6[2], r6[4])  # r6[px_], r6[py_]

    if abs(irho) < 1e-6
        NormL = L / pz
        r6[1] += r6[2] * NormL  # r6[x_]
        r6[3] += r6[4] * NormL  # r6[y_]
        r6[5] += NormL * dp1    # r6[ct_], absolute path length
    else
        pzmx = pz - (1.0 + r6[1] * irho)  # r6[x_]
        cs = cos(irho * L)
        sn = sin(irho * L)
        px = r6[2] * cs + pzmx * sn  # r6[px_]
        d2 = pxyz(dp1, 0.0, r6[4])  # r6[py_]
        dasin = L + (asin(r6[2] / d2) - asin(px / d2)) / irho
        x = (pxyz(dp1, px, r6[4]) - pzmx * cs + r6[2] * sn - 1.0) / irho  # r6[x_]
        dy = r6[4] * dasin  # r6[py_]
        dct = dp1 * dasin   # r6[ct_], absolute path length

        r6[1] = x
        r6[2] = px
        r6[3] += dy
        r6[5] += dct
    end
    return nothing
end

function bend_edge!(r6::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}, rhoinv::Float64, theta::Float64, beti::Float64) where {T, TPS_Dim, Max_TPS_Degree}
    # Forest 12.41, ideal wedge, map U(theta, rhoinv)

    if abs(rhoinv) >= 1e-6
        dp1 = 1.0*beti + r6[6]  # r6[delta_]
        c = cos(theta)
        s = sin(theta)
        pz = pxyz(dp1, r6[2], r6[4])  # r6[px_], r6[py_]
        d2 = pxyz(dp1, 0.0, r6[4])    # r6[py_]
        px = r6[2] * c + (pz - rhoinv * r6[1]) * s  # r6[px_]
        dasin = asin(r6[2] / d2) - asin(px / d2)
        num = r6[1] * (r6[2] * sin(2.0 * theta) + s * s * (2.0 * pz - rhoinv * r6[1]))
        den = pxyz(dp1, px, r6[4]) + pxyz(dp1, r6[2], r6[4]) * c - r6[2] * s
        x = r6[1] * c + num / den  # r6[x_]
        dy = r6[4] * theta / rhoinv + r6[4] / rhoinv * dasin  # r6[py_]
        dct = dp1 / rhoinv * (theta + dasin)  # r6[ct_]

        r6[1] = x
        r6[2] = px
        r6[3] += dy
        r6[5] += dct
    end
    return nothing
end

function exactbndthinkick_rad!(r::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}, A::Array{Float64,1}, B::Array{Float64,1}, 
    L::Float64, irho::Float64, max_order::Int, beti::Float64, rad_const::Float64) where {T, TPS_Dim, Max_TPS_Degree}
    # AT function. Ref[Terebilo, Andrei. "Accelerator modeling with MATLAB accelerator toolbox." PACS2001 (2001)].
    ReSum = CTPS(T(B[max_order + 1]), TPS_Dim, Max_TPS_Degree)
    ImSum = CTPS(T(A[max_order + 1]), TPS_Dim, Max_TPS_Degree)
    ReSumTemp = CTPS(T(0.0), TPS_Dim, Max_TPS_Degree)

    for i in max_order-1:-1:0
        ReSumTemp = ReSum * r[1] - ImSum * r[3] + B[i+1]
        ImSum = ImSum * r[1] + ReSum * r[3] + A[i+1]
        ReSum = ReSumTemp
    end

    p_norm = 1.0 / (1.0 + r[6])
    x = r[1]
    xpr = r[2] * p_norm
    y = r[3]
    ypr = r[4] * p_norm

    B2P = B2perp_exact_bnd(ImSum, ReSum+irho, irho, x, xpr, y, ypr)

    r[6] -= rad_const * (1.0 + r[6])^2 * B2P * (1.0 + x*irho) * L / sqrt(1.0 - xpr*xpr - ypr*ypr)

    p_norm = 1.0 / (1.0 + r[6])
    r[2] = xpr/p_norm
    r[4] = ypr/p_norm

    r[2] -= L * ReSum
    r[4] += L * ImSum
    return nothing
end

function ExactSectorBend!(r6::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}, le::Float64, beti::Float64, angle::Float64, A::Array{Float64,1}, B::Array{Float64,1}, 
    max_order::Int, num_int_steps::Int, entrance_angle::Float64, exit_angle::Float64, FringeBendEntrance::Int, FringeBendExit::Int,
    FringeQuadEntrance::Int, FringeQuadExit::Int, gk::Float64,
    T1::Array{Float64,1}, T2::Array{Float64,1}, 
    R1::Array{Float64,2}, R2::Array{Float64,2}, RApertures::Array{Float64,1}, EApertures::Array{Float64,1},
    KickAngle::Array{Float64,1}) where {T, TPS_Dim, Max_TPS_Degree}
    
    irho = angle / le
    DRIFT1 = 0.6756035959798286638
    DRIFT2 = -0.1756035959798286639
    KICK1 = 1.351207191959657328
    KICK2 = -1.702414383919314656
    SL = le / num_int_steps
    L1 = SL * DRIFT1
    L2 = SL * DRIFT2
    K1 = SL * KICK1
    K2 = SL * KICK2

    B0 = B[1]
    A0 = A[1]

    if !iszero(KickAngle[1])
        B[1] -= sin(KickAngle[1]) / le
    end
    if !iszero(KickAngle[2])
        A[1] += sin(KickAngle[2]) / le
    end


    # Misalignment at entrance
    if !iszero(T1)
        addvv!(r6, T1)
    end
    if !iszero(R1)
        multmv!(r6, R1)
    end

    Yrot!(r6, entrance_angle, beti)

    if FringeBendEntrance != 0
        bend_fringe!(r6, irho, gk, beti)
    end

    if FringeQuadEntrance != 0
        multipole_fringe!(r6, le, A, B, max_order, 1.0, 1, beti)
    end

    bend_edge!(r6, irho, -entrance_angle, beti)

    # Integrator
    if num_int_steps == 0
        exact_bend!(r6, irho, le, beti)
    else
        for m in 1:num_int_steps
            exact_bend!(r6, irho, L1, beti)
            strthinkick!(r6, A, B, K1, max_order)
            exact_bend!(r6, irho, L2, beti)
            strthinkick!(r6, A, B, K2, max_order)
            exact_bend!(r6, irho, L2, beti)
            strthinkick!(r6, A, B, K1, max_order)
            exact_bend!(r6, irho, L1, beti)
        end
    end

    r6[5] -= le*beti  # r6[ct_], absolute path length

    bend_edge!(r6, irho, -exit_angle, beti)
    if FringeQuadExit != 0
        multipole_fringe!(r6, le, A, B, max_order, -1.0, 1, beti)
    end
    if FringeBendExit != 0
        bend_fringe!(r6, -irho, gk, beti)
    end
    Yrot!(r6, exit_angle, beti)

    # Misalignment at exit
    if !iszero(R2)
        multmv!(r6, R2)
    end
    if !iszero(T2)
        addvv!(r6, T2)
    end

    if !iszero(KickAngle[1])
        B[1] += sin(KickAngle[1]) / le
    end
    if !iszero(KickAngle[2])
        A[1] -= sin(KickAngle[2]) / le
    end
    return nothing
end

function ExactSectorBend_rad!(r6::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}, le::Float64, rad_const::Float64, beti::Float64, angle::Float64, A::Array{Float64,1}, B::Array{Float64,1}, 
    max_order::Int, num_int_steps::Int, entrance_angle::Float64, exit_angle::Float64, FringeBendEntrance::Int, FringeBendExit::Int,
    FringeQuadEntrance::Int, FringeQuadExit::Int, gk::Float64,
    T1::Array{Float64,1}, T2::Array{Float64,1}, 
    R1::Array{Float64,2}, R2::Array{Float64,2}, RApertures::Array{Float64,1}, EApertures::Array{Float64,1},
    KickAngle::Array{Float64,1}) where {T, TPS_Dim, Max_TPS_Degree}
    
    irho = angle / le
    DRIFT1 = 0.6756035959798286638
    DRIFT2 = -0.1756035959798286639
    KICK1 = 1.351207191959657328
    KICK2 = -1.702414383919314656
    SL = le / num_int_steps
    L1 = SL * DRIFT1
    L2 = SL * DRIFT2
    K1 = SL * KICK1
    K2 = SL * KICK2

    B0 = B[1]
    A0 = A[1]

    if !iszero(KickAngle[1])
        B[1] -= sin(KickAngle[1]) / le
    end
    if !iszero(KickAngle[2])
        A[1] += sin(KickAngle[2]) / le
    end


    # Misalignment at entrance
    if !iszero(T1)
        addvv!(r6, T1)
    end
    if !iszero(R1)
        multmv!(r6, R1)
    end

    Yrot!(r6, entrance_angle, beti)

    if FringeBendEntrance != 0
        bend_fringe!(r6, irho, gk, beti)
    end

    if FringeQuadEntrance != 0
        multipole_fringe!(r6, le, A, B, max_order, 1.0, 1, beti)
    end

    bend_edge!(r6, irho, -entrance_angle, beti)

    # Integrator
    if num_int_steps == 0
        exact_bend!(r6, irho, le, beti)
    else
        for m in 1:num_int_steps
            exact_bend!(r6, irho, L1, beti)
            exactbndthinkick_rad!(r6, A, B, K1, irho, max_order, beti, rad_const)
            exact_bend!(r6, irho, L2, beti)
            exactbndthinkick_rad!(r6, A, B, K2, irho, max_order, beti, rad_const)
            exact_bend!(r6, irho, L2, beti)
            exactbndthinkick_rad!(r6, A, B, K1, irho, max_order, beti, rad_const)
            exact_bend!(r6, irho, L1, beti)
        end
    end

    r6[5] -= le*beti  # r6[ct_], absolute path length

    bend_edge!(r6, irho, -exit_angle, beti)
    if FringeQuadExit != 0
        multipole_fringe!(r6, le, A, B, max_order, -1.0, 1, beti)
    end
    if FringeBendExit != 0
        bend_fringe!(r6, -irho, gk, beti)
    end
    Yrot!(r6, exit_angle, beti)

    # Misalignment at exit
    if !iszero(R2)
        multmv!(r6, R2)
    end
    if !iszero(T2)
        addvv!(r6, T2)
    end
    
    if !iszero(KickAngle[1])
        B[1] += sin(KickAngle[1]) / le
    end
    if !iszero(KickAngle[2])
        A[1] -= sin(KickAngle[2]) / le
    end
    return nothing
end

function pass_TPSA!(ele::ESBEND, r_in::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}; E0::Float64=0.0, m0::Float64=m_e) where {T, TPS_Dim, Max_TPS_Degree}
    gamma = E0 / m0
    beta = sqrt(1 - 1 / (gamma * gamma))
    rad_const = 0.0
    if use_exact_beti == 1
        beti = 1.0 / beta
    else
        beti = 1.0 
    end
    if ele.rad == 0
        ExactSectorBend!(r_in, ele.len, beti, ele.angle, ele.PolynomA, ele.PolynomB, ele.MaxOrder, ele.NumIntSteps,
            ele.e1, ele.e2,
            ele.FringeBendEntrance, ele.FringeBendExit,
            ele.FringeQuadEntrance, ele.FringeQuadExit,
            ele.gK,
            ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures,
            ele.KickAngle)
    else
        if m0 == m_e
            rad_const = RAD_CONST_E * gamma^3
        elseif m0 == m_p
            rad_const = RAD_CONST_P * gamma^3
        else
            rad_const = 0.0
            println("SR is not implemented for this particle mass.")
        end
        ExactSectorBend_rad!(r_in, ele.len, rad_const, beti, ele.angle, ele.PolynomA, ele.PolynomB, ele.MaxOrder, ele.NumIntSteps,
            ele.e1, ele.e2,
            ele.FringeBendEntrance, ele.FringeBendExit,
            ele.FringeQuadEntrance, ele.FringeQuadExit,
            ele.gK,
            ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures,
            ele.KickAngle)
    end
    return nothing
end