function bndthinkickrad!(r::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}, A::Array{Float64,1}, B::Array{Float64,1}, 
    L::Float64, irho::Float64, E0::Float64, max_order::Int, beti::Float64) where {T, TPS_Dim, Max_TPS_Degree}
    ReSum = CTPS(T(B[max_order + 1]), TPS_Dim, Max_TPS_Degree)
    ImSum = CTPS(T(A[max_order + 1]), TPS_Dim, Max_TPS_Degree)
    ReSumTemp = CTPS(zero(T), TPS_Dim, Max_TPS_Degree)
    CRAD = CGAMMA * E0^3 / (2.0*pi*1e27) # [m]/[GeV^3] M.Sands (4.1)

    for i in reverse(1:max_order)
        ReSumTemp = ReSum * r[1] - ImSum * r[3] + B[i]
        ImSum = ImSum * r[1] + ReSum * r[3] + A[i]
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
    r[5] += L * irho * r[1] * beti
    return nothing
end

function bndthinkick!(
    r::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}},
    A::Vector{Float64}, B::Vector{Float64},
    L::Float64, irho::Float64,
    max_order::Int, beti::Float64) where {T, TPS_Dim, Max_TPS_Degree}
  
    ReSum = CTPS(T, TPS_Dim, Max_TPS_Degree)
    ImSum = CTPS(T, TPS_Dim, Max_TPS_Degree)
    ReSum.map[1] = B[max_order + 1]
    ImSum.map[1] = A[max_order + 1]

    for i in max_order:-1:1
        # always stay in CTPS land:
        if i == max_order
            ReSumTemp = B[max_order + 1] * r[1] - A[max_order + 1] * r[3] + B[i]
            ImSum     = A[max_order + 1] * r[1] + B[max_order + 1] * r[3] + A[i]
            ReSum      = CTPS(ReSumTemp)
        else
            ReSumTemp = ReSum * r[1] - ImSum * r[3] + B[i]
            ImSum     = ImSum * r[1] + ReSum * r[3] + A[i]
            ReSum      = ReSumTemp
        end
    end

    r[2] -= L * (ReSum - (r[6] - r[1]*irho) * irho)
    r[4] += L *  ImSum
    r[5] += L * (irho * r[1]) * beti
    return nothing
end
  

function BendSymplecticPassRad!(r::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}, le::Float64, beti::Float64, irho::Float64, A::Array{Float64,1}, B::Array{Float64,1}, 
    max_order::Int, num_int_steps::Int, entrance_angle::Float64, exit_angle::Float64, FringeBendEntrance::Int, FringeBendExit::Int,
    fint1::Float64, fint2::Float64, gap::Float64, FringeQuadEntrance::Int, FringeQuadExit::Int,
    fringeIntM0::Array{Float64,1}, fringeIntP0::Array{Float64,1}, T1::Array{Float64,1}, T2::Array{Float64,1}, 
    R1::Array{Float64,2}, R2::Array{Float64,2}, RApertures::Array{Float64,1}, EApertures::Array{Float64,1},
    KickAngle::Array{Float64,1}, E0::Float64) where {T, TPS_Dim, Max_TPS_Degree}

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


    # Misalignment at entrance
    if !iszero(T1)
        addvv!(r, T1)
    end
    if !iszero(R1)
        multmv!(r, R1)
    end
    # Edge focus at entrance
    edge_fringe_entrance!(r, irho, entrance_angle, fint1, gap, FringeBendEntrance)

    # Quadrupole gradient fringe entrance
    if !iszero(FringeQuadEntrance) && !iszero(B[2])
        if useLinFrEleEntrance == 1
            linearQuadFringeElegantEntrance!(r, B[2], fringeIntM0, fringeIntP0)
        else
            QuadFringePassP!(r, B[2])
        end
    end

    # Integrator
    for m in 1:num_int_steps
        drift6!(r, L1, beti)
        bndthinkickrad!(r, A, B, K1, irho, E0, max_order, beti)
        drift6!(r, L2, beti)
        bndthinkickrad!(r, A, B, K2, irho, E0, max_order, beti)
        drift6!(r, L2, beti)
        bndthinkickrad!(r, A, B, K1, irho, E0, max_order, beti)
        drift6!(r, L1, beti)
    end

    # Quadrupole gradient fringe exit
    if !iszero(FringeQuadExit) && !iszero(B[2])
        if useLinFrEleExit == 1
            linearQuadFringeElegantExit!(r, B[2], fringeIntM0, fringeIntP0)
        else
            QuadFringePassN!(r, B[2])
        end
    end

    # Edge focus at exit
    edge_fringe_exit!(r, irho, exit_angle, fint2, gap, FringeBendExit)

    # Misalignment at exit
    if !iszero(R2)
        multmv!(r, R2)
    end
    if !iszero(T2)
        addvv!(r, T2)
    end
    
    B[1] += sin(KickAngle[1]) / le
    A[1] -= sin(KickAngle[2]) / le
    return nothing
end

function BendSymplecticPass!(r::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}, le::Float64, beti::Float64, irho::Float64, 
    A::Array{Float64,1}, B::Array{Float64,1}, max_order::Int, num_int_steps::Int,
    entrance_angle::Float64, exit_angle::Float64,
    FringeBendEntrance::Int, FringeBendExit::Int,
    fint1::Float64, fint2::Float64, gap::Float64,
    FringeQuadEntrance::Int, FringeQuadExit::Int,
    fringeIntM0::Array{Float64,1}, fringeIntP0::Array{Float64,1},
    T1::Array{Float64,1}, T2::Array{Float64,1}, R1::Array{Float64,2}, R2::Array{Float64,2}, 
    RApertures::Array{Float64,1}, EApertures::Array{Float64,1}, KickAngle::Array{Float64,1}) where {T, TPS_Dim, Max_TPS_Degree}

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

    # Misalignment at entrance
    if !iszero(T1)
        addvv!(r, T1)
    end
    if !iszero(R1)
        multmv!(r, R1)
    end

    # Edge focus at entrance
    edge_fringe_entrance!(r, irho, entrance_angle, fint1, gap, FringeBendEntrance)

    # Quadrupole gradient fringe entrance
    if !iszero(FringeQuadEntrance) && !iszero(B[2])
        if useLinFrEleEntrance == 1
            linearQuadFringeElegantEntrance!(r, B[2], fringeIntM0, fringeIntP0)
        else
            QuadFringePassP!(r, B[2])
        end
    end

    # Integrator
    for m in 1:num_int_steps
        drift6!(r, L1, beti)
        bndthinkick!(r, A, B, K1, irho, max_order, beti)
        drift6!(r, L2, beti)
        bndthinkick!(r, A, B, K2, irho, max_order, beti)
        drift6!(r, L2, beti)
        bndthinkick!(r, A, B, K1, irho, max_order, beti)
        drift6!(r, L1, beti)
    end

    # Quadrupole gradient fringe exit
    if !iszero(FringeQuadExit) && !iszero(B[2])
        if useLinFrEleExit == 1
            linearQuadFringeElegantExit!(r, B[2], fringeIntM0, fringeIntP0)
        else
            QuadFringePassN!(r, B[2])
        end
    end

    # Edge focus at exit
    edge_fringe_exit!(r, irho, exit_angle, fint2, gap, FringeBendExit)

    # Misalignment at exit
    if !iszero(R2)
        multmv!(r, R2)
    end
    if !iszero(T2)
        addvv!(r, T2)
    end

    B[1] += sin(KickAngle[1]) / le
    A[1] -= sin(KickAngle[2]) / le
    return nothing
end
    
    
function pass_TPSA!(ele::SBEND, r_in::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}; E0::Float64=0.0, m0::Float64=m_e) where {T, TPS_Dim, Max_TPS_Degree}
    # ele: SBEND
    # r_in: 6-by-1 TPSA array
    irho = ele.angle / ele.len
    gamma = (E0 + m0) / m0
    beta = sqrt(1.0 - 1.0 / gamma^2)
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
                ele.KickAngle)
    else
        BendSymplecticPassRad!(r_in, ele.len, beti, irho, ele.PolynomA, ele.PolynomB, ele.MaxOrder, ele.NumIntSteps,
                ele.e1, ele.e2,
                ele.FringeBendEntrance, ele.FringeBendExit,
                ele.fint1, ele.fint2, ele.gap,
                ele.FringeQuadEntrance, ele.FringeQuadExit,
                ele.FringeIntM0, ele.FringeIntP0,
                ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures,
                ele.KickAngle, E0)
    end
    return nothing
end
