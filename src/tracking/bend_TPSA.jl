# include("fringe_TPSA.jl")

function bndthinkick!(r::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}, A, B, L, irho, max_order) where {T, TPS_Dim, Max_TPS_Degree}
    ### this section is type unstable
    # ReSum = B[max_order + 1]
    # ImSum = A[max_order + 1]
    # ReSumTemp = 0.0
    
    # i = 1
    ReSumTemp = B[max_order + 1] * r[1] - A[max_order + 1] * r[3] + B[1]
    ImSum = A[max_order + 1] * r[1] + B[max_order + 1] * r[3] + A[1]
    ReSum = CTPS(ReSumTemp)
    for i in reverse(2:max_order)
        ReSumTemp = ReSum * r[1] - ImSum * r[3] + B[i]
        ImSum = ImSum * r[1] + ReSum * r[3] + A[i]
        ReSum = CTPS(ReSumTemp)
    end

    r[2] -= L * (ReSum - (r[6] - r[1] * irho) * irho)
    r[4] += L * ImSum
    r[5] += L * irho * r[1]  # Path length
    return nothing
end
    
function BendSymplecticPass!(r::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}, le, irho, A, B, max_order, num_int_steps,
    entrance_angle, exit_angle,
    FringeBendEntrance, FringeBendExit,
    fint1, fint2, gap,
    FringeQuadEntrance, FringeQuadExit,
    fringeIntM0, fringeIntP0,
    T1, T2, R1, R2, RApertures, EApertures,
    KickAngle) where {T, TPS_Dim, Max_TPS_Degree}

    DRIFT1 = 0.6756035959798286638
    DRIFT2 = -0.1756035959798286639
    KICK1 = 1.351207191959657328
    KICK2 = -1.702414383919314656
    SL = le / num_int_steps
    L1 = SL * DRIFT1
    L2 = SL * DRIFT2
    K1 = SL * KICK1
    K2 = SL * KICK2

    # if FringeQuadEntrance==2 #&& !isnothing(fringeIntM0) && !isnothing(fringeIntP0)
    #     useLinFrEleEntrance = 1
    # else
    #     useLinFrEleEntrance = 0
    # end
    # if FringeQuadExit==2 #&& !isnothing(fringeIntM0) && !isnothing(fringeIntP0)
    #     useLinFrEleExit = 1
    # else
    #     useLinFrEleExit = 0
    # end

    B[1] -= sin(KickAngle[1]) / le
    A[1] += sin(KickAngle[2]) / le

    if isone(use_exact_Hamiltonian)
        NormL1 = L1 / sqrt((1.0 + r[6])^2 - r[2]^2 - r[4]^2)
        NormL2 = L2 / sqrt((1.0 + r[6])^2 - r[2]^2 - r[4]^2)
    else
        NormL1 = L1 / (1.0 + r[6])
        NormL2 = L2 / (1.0 + r[6])
    end
    # Misalignment at entrance
    if !iszero(T1)
        addvv!(r, T1)
    end
    if !iszero(R1)
        multmv!(r, R1)
    end

    # Edge focus at entrance
    edge_fringe_entrance!(r, irho, entrance_angle, fint1, gap, FringeBendEntrance)

    # Integrator
    for m in 1:num_int_steps
        fastdrift!(r, NormL1, L1)
        bndthinkick!(r, A, B, K1, irho, max_order)
        fastdrift!(r, NormL2, L2)
        bndthinkick!(r, A, B, K2, irho, max_order)
        fastdrift!(r, NormL2, L2)
        bndthinkick!(r, A, B, K1, irho, max_order)
        fastdrift!(r, NormL1, L1)
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
    
    
function pass_TPSA!(ele::SBEND, r_in::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}) where {T, TPS_Dim, Max_TPS_Degree}
    # ele: SBEND
    # r_in: 6-by-1 TPSA array

    irho = ele.angle / ele.len
    BendSymplecticPass!(r_in, ele.len, irho, ele.PolynomA, ele.PolynomB, ele.MaxOrder, ele.NumIntSteps,
            ele.e1, ele.e2,
            ele.FringeBendEntrance, ele.FringeBendExit,
            ele.fint1, ele.fint2, ele.gap,
            ele.FringeQuadEntrance, ele.FringeQuadExit,
            ele.FringeIntM0, ele.FringeIntP0,
            ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures,
            ele.KickAngle)
    return nothing
end
# function pass_TPSA!(ele::RBEND, r_in::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}) where {T, TPS_Dim, Max_TPS_Degree}
#     # ele: RBEND
#     # r_in: 6-by-1 TPSA array

#     irho = ele.angle / ele.len
#     BendSymplecticPass!(r_in, ele.len, irho, ele.PolynomA, ele.PolynomB, ele.MaxOrder, ele.NumIntSteps,
#             ele.e1, ele.e2,
#             ele.FringeBendEntrance, ele.FringeBendExit,
#             ele.fint1, ele.fint2, ele.gap,
#             ele.FringeQuadEntrance, ele.FringeQuadExit,
#             ele.FringeIntM0, ele.FringeIntP0,
#             ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures,
#             ele.KickAngle)
#     return nothing
# end