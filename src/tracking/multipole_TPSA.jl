function StrB2perp(bx::CTPS{T, TPS_Dim, Max_TPS_Degree}, by::CTPS{T, TPS_Dim, Max_TPS_Degree}, 
                x::CTPS{T, TPS_Dim, Max_TPS_Degree}, xpr::CTPS{T, TPS_Dim, Max_TPS_Degree}, 
                y::CTPS{T, TPS_Dim, Max_TPS_Degree}, ypr::CTPS{T, TPS_Dim, Max_TPS_Degree}) where {T, TPS_Dim, Max_TPS_Degree}
    # Calculates sqr(|B x e|) , where e is a unit vector in the direction of velocity
    # v_norm2 = 1.0 / (1.0 + xpr^2 + ypr^2)
    # return (by^2 + bx^2 + (bx*ypr - by*xpr)^2) * v_norm2
    return bx*bx + by*by + (bx*xpr - by*ypr)^2
end
function strthinkickrad!(r::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}, A::AbstractVector{Float64}, B::AbstractVector{Float64},
                        L::Float64, E0::Float64, max_order::Int, rad_const::Float64) where {T, TPS_Dim, Max_TPS_Degree}
    # Modified based on AT function. Ref[Terebilo, Andrei. "Accelerator modeling with MATLAB accelerator toolbox." PACS2001 (2001)].
    ReSum = CTPS(T(B[max_order + 1]), TPS_Dim, Max_TPS_Degree)
    ImSum = CTPS(T(A[max_order + 1]), TPS_Dim, Max_TPS_Degree)
    ReSumTemp = CTPS(zero(T), TPS_Dim, Max_TPS_Degree)

    for i in reverse(1:max_order)
        ReSumTemp = ReSum * r[1] - ImSum * r[3] + B[i]
        ImSum = ImSum * r[1] + ReSum * r[3] + A[i]
        ReSum = ReSumTemp
    end

    # angles for momentum
    p_norm = 1.0 / (1.0 + r[6])

    x = r[1]
    xpr = r[2] * p_norm
    y = r[3]
    ypr = r[4] * p_norm
    B2P = StrB2perp(ImSum, ReSum , x , xpr, y ,ypr)
    factor = L / (p_norm)^2 / sqrt(1.0 - xpr^2 - ypr^2)

    r[6] -= rad_const * B2P * factor

    # momentums after losing energy
    p_norm = 1.0 / (1.0 + r[6])

    r[2] = xpr / p_norm
    r[4] = ypr / p_norm

    r[2] -= L * ReSum
    r[4] += L * ImSum
    return nothing
end

function StrMPoleSymplectic4RadPass(r::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}, le::Float64, rad_const::Float64, beti::Float64,
    A::Array{Float64,1}, B::Array{Float64,1}, 
    max_order::Int, num_int_step::Int, 
    FringeQuadEntrance::Int, FringeQuadExit::Int, #(no fringe), 1 (lee-whiting) or 2 (lee-whiting+elegant-like) 
    T1::Array{Float64,1}, T2::Array{Float64,1}, R1::Array{Float64,2}, R2::Array{Float64,2}, 
    RApertures::Array{Float64,1}, EApertures::Array{Float64,1}, KickAngle::Array{Float64,1}, E0::Float64) where {T, TPS_Dim, Max_TPS_Degree}

    DRIFT1  =  0.6756035959798286638
    DRIFT2 = -0.1756035959798286639
    KICK1 = 1.351207191959657328
    KICK2 = -1.702414383919314656

    SL = le/num_int_step
    L1 = SL*DRIFT1
    L2 = SL*DRIFT2
    K1 = SL*KICK1
    K2 = SL*KICK2

    # if FringeQuadEntrance==2 && !isnothing(fringeIntM0) && !isnothing(fringeIntP0)
    #     useLinFrEleEntrance = 1
    # else
    #     useLinFrEleEntrance = 0
    # end
    # if FringeQuadExit==2 && !isnothing(fringeIntM0) && !isnothing(fringeIntP0)
    #     useLinFrEleExit = 1
    # else
    #     useLinFrEleExit = 0
    # end

    if le > 0
        B[1] -= sin(KickAngle[1])/le
        A[1] += sin(KickAngle[2])/le
    end
    # Misalignment at entrance
    if !iszero(T1)
        addvv!(r, T1)
    end
    if !iszero(R1)
        multmv!(r, R1)
    end

    # Check physical apertures at the entrance of the magnet
    # if RApertures != nothing
    #     checkiflostRectangularAp(r, RApertures)
    # end
    # if EApertures != nothing
    #     checkiflostEllipticalAp(r, EApertures)
    # end

    # if FringeQuadEntrance != 0 && B[2] != 0
    #     if useLinFrEleEntrance == 1
    #         linearQuadFringeElegantEntrance!(r, B[2], fringeIntM0, fringeIntP0)
    #     else
    #         QuadFringePassP!(r, B[2])
    #     end
    # end
    if FringeQuadEntrance != 0 
        multipole_fringe!(r6, le, A, B, max_order, 1.0, 1, beti)
    end


    # Integrator
    for m in 1:num_int_step
        drift6!(r, L1, beti)
        strthinkickrad!(r, A, B, K1, E0, max_order, rad_const)
        drift6!(r, L2, beti)
        strthinkickrad!(r, A, B, K2, E0, max_order, rad_const)
        drift6!(r, L2, beti)
        strthinkickrad!(r, A, B, K1, E0, max_order, rad_const)
        drift6!(r, L1, beti)
    end
    # if FringeQuadExit != 0 && B[2] != 0
    #     if useLinFrEleExit == 1
    #         linearQuadFringeElegantExit!(r, B[2], fringeIntM0, fringeIntP0)
    #     else
    #         QuadFringePassN!(r, B[2])
    #     end
    # end
    if FringeQuadExit != 0 
        multipole_fringe!(r6, le, A, B, max_order, -1.0, 1, beti)
    end

    # Check physical apertures at the exit of the magnet
    # if RApertures != nothing
    #     checkiflostRectangularAp(r, RApertures)
    # end
    # if EApertures != nothing
    #     checkiflostEllipticalAp(r, EApertures)
    # end

    # Misalignment at exit
    if !iszero(R2)
        multmv!(r, R2)
    end
    if !iszero(T2)
        addvv!(r, T2)
    end

    if le > 0
        B[1] += sin(KickAngle[1]) / le
        A[1] -= sin(KickAngle[2]) / le
    end
    return nothing
end

function strthinkick!(r::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}, A::AbstractVector{Float64}, B::AbstractVector{Float64}, 
                        L::Float64, max_order::Int) where {T, TPS_Dim, Max_TPS_Degree}
    ReSum = CTPS(T(B[max_order + 1]), TPS_Dim, Max_TPS_Degree)
    ImSum = CTPS(T(A[max_order + 1]), TPS_Dim, Max_TPS_Degree)
    ReSumTemp = CTPS(zero(T), TPS_Dim, Max_TPS_Degree)

    for i in max_order-1: -1: 0
        ReSumTemp = ReSum * r[1] - ImSum * r[3] + B[i+1]
        ImSum = ImSum * r[1] + ReSum * r[3] + A[i+1]
        ReSum = ReSumTemp
    end

    r[2] -= L * ReSum
    r[4] += L * ImSum
    return nothing
end


function StrMPoleSymplectic4Pass!(r::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}, le::Float64, beti::Float64, 
    A::AbstractVector{Float64}, B::AbstractVector{Float64}, max_order::Int, num_int_step::Int, 
    FringeQuadEntrance::Int, FringeQuadExit::Int, #(no fringe), 1 (lee-whiting) or 2 (lee-whiting+elegant-like) 
    T1::AbstractVector{Float64}, T2::AbstractVector{Float64}, R1::Array{Float64,2}, R2::Array{Float64,2}, 
    RApertures::AbstractVector{Float64}, EApertures::AbstractVector{Float64}, KickAngle::AbstractVector{Float64}) where {T, TPS_Dim, Max_TPS_Degree}

    DRIFT1  =  0.6756035959798286638
    DRIFT2 = -0.1756035959798286639
    KICK1 = 1.351207191959657328
    KICK2 = -1.702414383919314656

    SL = le/num_int_step
    L1 = SL*DRIFT1
    L2 = SL*DRIFT2
    K1 = SL*KICK1
    K2 = SL*KICK2

    # if FringeQuadEntrance==2 && !isnothing(fringeIntM0) && !isnothing(fringeIntP0)
    #     useLinFrEleEntrance = 1
    # else
    #     useLinFrEleEntrance = 0
    # end
    # if FringeQuadExit==2 && !isnothing(fringeIntM0) && !isnothing(fringeIntP0)
    #     useLinFrEleExit = 1
    # else
    #     useLinFrEleExit = 0
    # end

    if le > 0
        B[1] -= sin(KickAngle[1])/le
        A[1] += sin(KickAngle[2])/le
    end

    # Misalignment at entrance
    if !iszero(T1)
        addvv!(r, T1)
    end
    if !iszero(R1)
        multmv!(r, R1)
    end
    if FringeQuadEntrance != 0 
        multipole_fringe!(r6, le, A, B, max_order, 1.0, 1, beti)
    end

    # Integrator
    for m in 1:num_int_step
        drift6!(r, L1, beti)
        strthinkick!(r, A, B, K1, max_order)
        drift6!(r, L2, beti)
        strthinkick!(r, A, B, K2, max_order)
        drift6!(r, L2, beti)
        strthinkick!(r, A, B, K1, max_order)
        drift6!(r, L1, beti)
    end

    if FringeQuadExit != 0 
        multipole_fringe!(r6, le, A, B, max_order, -1.0, 1, beti)
    end
    # Misalignment at exit
    if !iszero(R2)
        multmv!(r, R2)
    end
    if !iszero(T2)
        addvv!(r, T2)
    end            
        # end
    # end
    if le > 0
        B[1] += sin(KickAngle[1]) / le
        A[1] -= sin(KickAngle[2]) / le
    end
    return nothing
end

function pass_TPSA!(ele::KQUAD, r_in::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}; E0::Float64=0.0, m0::Float64=m_e) where {T, TPS_Dim, Max_TPS_Degree}
    # ele: KQUAD
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    gamma = E0 / m0
    beta = sqrt(1.0 - 1.0 / (gamma^2))
    if use_exact_beti == 1
        beti = 1.0 / beta
    else
        beti = 1.0 
    end
    rad_const = 0.0
    PolynomB = zeros(4)
    if ele.PolynomB[1] == 0.0 && ele.PolynomB[2] == 0.0 && ele.PolynomB[3] == 0.0 && ele.PolynomB[4] == 0.0
        PolynomB[2] = ele.k1
        if ele.rad == 0
            StrMPoleSymplectic4Pass!(r_in, ele.len, beti, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit, 
                ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle)
        else
            if m0 == m_e
                rad_const = RAD_CONST_E * gamma^3
            elseif m0 == m_p
                rad_const = RAD_CONST_P * gamma^3
            else
                rad_const = 0.0
                println("SR is not implemented for this particle mass.")
            end
            StrMPoleSymplectic4RadPass(r_in, ele.len,rad_const, beti, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit, 
                ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, E0)
        end
    else
        PolynomB[1] = ele.PolynomB[1]
        PolynomB[2] = ele.PolynomB[2] 
        PolynomB[3] = ele.PolynomB[3] / 2.0
        PolynomB[4] = ele.PolynomB[4] / 6.0
        if ele.rad == 0
            StrMPoleSymplectic4Pass!(r_in, ele.len, beti, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit, 
                ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle)
        else
            if m0 == m_e
                rad_const = RAD_CONST_E * gamma^3
            elseif m0 == m_p
                rad_const = RAD_CONST_P * gamma^3
            else
                rad_const = 0.0
                println("SR is not implemented for this particle mass.")
            end
            StrMPoleSymplectic4RadPass(r_in, ele.len, rad_const, beti, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit, 
                ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, E0)
        end
    end
    return nothing
end

function pass_TPSA!(ele::KSEXT, r_in::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}; E0::Float64=0.0, m0::Float64=m_e) where {T, TPS_Dim, Max_TPS_Degree}
    # ele: KSEXT
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    gamma = E0 / m0
    beta = sqrt(1.0 - 1.0 / (gamma^2))
    if use_exact_beti == 1
        beti = 1.0 /beta
    else
        beti = 1.0 
    end
    PolynomB = zeros(4)
    rad_const = 0.0

    if ele.PolynomB[1] == 0.0 && ele.PolynomB[2] == 0.0 && ele.PolynomB[3] == 0.0 && ele.PolynomB[4] == 0.0
        PolynomB[3] = ele.k2 / 2.0
        if ele.rad == 0
            StrMPoleSymplectic4Pass!(r_in, ele.len, beti, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit,
                ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle)
        else
            if m0 == m_e
                rad_const = RAD_CONST_E * gamma^3
            elseif m0 == m_p
                rad_const = RAD_CONST_P * gamma^3
            else
                rad_const = 0.0
                println("SR is not implemented for this particle mass.")
            end
            StrMPoleSymplectic4RadPass(r_in, ele.len, rad_const, beti, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit, 
                ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, E0)
        end
    else
        PolynomB[1] = ele.PolynomB[1]
        PolynomB[2] = ele.PolynomB[2] 
        PolynomB[3] = ele.PolynomB[3] / 2.0
        PolynomB[4] = ele.PolynomB[4] / 6.0
        if ele.rad == 0
            StrMPoleSymplectic4Pass!(r_in, ele.len, beti, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit, 
                ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle)
        else
            if m0 == m_e
                rad_const = RAD_CONST_E * gamma^3
            elseif m0 == m_p
                rad_const = RAD_CONST_P * gamma^3
            else
                rad_const = 0.0
                println("SR is not implemented for this particle mass.")
            end
            StrMPoleSymplectic4RadPass(r_in, ele.len, rad_const, beti, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit, 
                ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, E0)
        end
    end
    return nothing 
end

function pass_TPSA!(ele::KOCT, r_in::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}; E0::Float64=0.0, m0::Float64=m_e) where {T, TPS_Dim, Max_TPS_Degree}
    # ele: KOCT
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    gamma = E0 / m0
    beta = sqrt(1.0 - 1.0 / (gamma^2))
    if use_exact_beti == 1
        beti = 1.0 / beta
    else
        beti = 1.0 
    end
    rad_const = 0.0
    PolynomB = zeros(4)
    if ele.PolynomB[1] == 0.0 && ele.PolynomB[2] == 0.0 && ele.PolynomB[3] == 0.0 && ele.PolynomB[4] == 0.0
        PolynomB[4] = ele.k3 / 6.0
        if ele.rad == 0
            StrMPoleSymplectic4Pass!(r_in, ele.len, beti, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit, 
                ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle)
        else
            if m0 == m_e
                rad_const = RAD_CONST_E * gamma^3
            elseif m0 == m_p
                rad_const = RAD_CONST_P * gamma^3
            else
                rad_const = 0.0
                println("SR is not implemented for this particle mass.")
            end
            StrMPoleSymplectic4RadPass(r_in, ele.len, rad_const, beti, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit, 
                ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, E0)
        end
    else
        PolynomB[1] = ele.PolynomB[1]
        PolynomB[2] = ele.PolynomB[2] 
        PolynomB[3] = ele.PolynomB[3] / 2.0
        PolynomB[4] = ele.PolynomB[4] / 6.0
        if ele.rad == 0
            StrMPoleSymplectic4Pass!(r_in, ele.len, beti, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit, 
                ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle)
        else
            if m0 == m_e
                rad_const = RAD_CONST_E * gamma^3
            elseif m0 == m_p
                rad_const = RAD_CONST_P * gamma^3
            else
                rad_const = 0.0
                println("SR is not implemented for this particle mass.")
            end
            StrMPoleSymplectic4RadPass(r_in, ele.len, rad_const, beti, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit, 
                ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, E0)
        end
    end
    return nothing
end
