# include("fringe_TPSA.jl")

function strthinkick!(r::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}, A, B, L, max_order) where {T, TPS_Dim, Max_TPS_Degree}
    ### this section is type unstable
    # ReSum = B[max_order + 1]
    # ImSum = A[max_order + 1]
    # ReSumTemp = zeros(length(r[1]))

    # i = 1
    ReSumTemp = B[max_order + 1] * r[1] - A[max_order + 1] * r[3] + B[1]
    ImSum = A[max_order + 1] * r[1] + B[max_order + 1] * r[3] + A[1]
    ReSum = CTPS(ReSumTemp)

    for i in reverse(2:max_order)
        ReSumTemp = ReSum * r[1] - ImSum * r[3] + B[i]
        ImSum = ImSum * r[1] + ReSum * r[3] + A[i]
        ReSum = CTPS(ReSumTemp)
    end

    r[2] -= L * ReSum
    r[4] += L * ImSum
    return nothing
end


function StrMPoleSymplectic4Pass!(r::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}, le, A, B, max_order, num_int_step, 
    FringeQuadEntrance, FringeQuadExit, #(no fringe), 1 (lee-whiting) or 2 (lee-whiting+elegant-like) 
    fringeIntM0,  # I0m/K1, I1m/K1, I2m/K1, I3m/K1, Lambda2m/K1 
    fringeIntP0,  # I0p/K1, I1p/K1, I2p/K1, I3p/K1, Lambda2p/K1
    T1, T2, R1, R2, RApertures, EApertures, KickAngle) where {T, TPS_Dim, Max_TPS_Degree}

    DRIFT1  =  0.6756035959798286638
    DRIFT2 = -0.1756035959798286639
    KICK1 = 1.351207191959657328
    KICK2 = -1.702414383919314656

    SL = le/num_int_step
    L1 = SL*DRIFT1
    L2 = SL*DRIFT2
    K1 = SL*KICK1
    K2 = SL*KICK2

    if FringeQuadEntrance==2 && !isnothing(fringeIntM0) && !isnothing(fringeIntP0)
        useLinFrEleEntrance = 1
    else
        useLinFrEleEntrance = 0
    end
    if FringeQuadExit==2 && !isnothing(fringeIntM0) && !isnothing(fringeIntP0)
        useLinFrEleExit = 1
    else
        useLinFrEleExit = 0
    end

    if le > 0
        B[1] -= sin(KickAngle[1])/le
        A[1] += sin(KickAngle[2])/le
    end
    # Threads.@threads for c in 1:num_particles
    # for c in 1:num_particles
    if use_exact_Hamiltonian == 1
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

    # Integrator
    for m in 1:num_int_step
        fastdrift!(r, NormL1, L1)
        strthinkick!(r, A, B, K1, max_order)
        fastdrift!(r, NormL2, L2)
        strthinkick!(r, A, B, K2, max_order)
        fastdrift!(r, NormL2, L2)
        strthinkick!(r, A, B, K1, max_order)
        fastdrift!(r, NormL1, L1)
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

function pass_TPSA!(ele::KQUAD, r_in::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}) where {T, TPS_Dim, Max_TPS_Degree}
    # ele: KQUAD
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    PolynomB = zeros(4)
    # if ele.rad != 0
    #     println("Synchrtron radiation is not implemented for Quad in TPSA")
    # end
    if ele.PolynomB[1] == 0.0 && ele.PolynomB[2] == 0.0 && ele.PolynomB[3] == 0.0 && ele.PolynomB[4] == 0.0
        PolynomB[2] = ele.k1
        StrMPoleSymplectic4Pass!(r_in, ele.len, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
            ele.FringeQuadEntrance, ele.FringeQuadExit, ele.FringeIntM0, ele.FringeIntP0, 
            ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle)
    else
        PolynomB[1] = ele.PolynomB[1]
        PolynomB[2] = ele.PolynomB[2] 
        PolynomB[3] = ele.PolynomB[3] / 2.0
        PolynomB[4] = ele.PolynomB[4] / 6.0
        StrMPoleSymplectic4Pass!(r_in, ele.len, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
            ele.FringeQuadEntrance, ele.FringeQuadExit, ele.FringeIntM0, ele.FringeIntP0, 
            ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle)
    end
    return nothing
end

function pass_TPSA!(ele::KSEXT, r_in::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}) where {T, TPS_Dim, Max_TPS_Degree}
    # ele: KSEXT
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    PolynomB = zeros(4)
    # if ele.rad != 0
    #     println("Synchrtron radiation is not implemented for Sext in TPSA")
    # end
    if ele.PolynomB[1] == 0.0 && ele.PolynomB[2] == 0.0 && ele.PolynomB[3] == 0.0 && ele.PolynomB[4] == 0.0
        PolynomB[3] = ele.k2 / 2.0
        StrMPoleSymplectic4Pass!(r_in, ele.len, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
            ele.FringeQuadEntrance, ele.FringeQuadExit, ele.FringeIntM0, ele.FringeIntP0, 
            ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle)
    else
        PolynomB[1] = ele.PolynomB[1]
        PolynomB[2] = ele.PolynomB[2] 
        PolynomB[3] = ele.PolynomB[3] / 2.0
        PolynomB[4] = ele.PolynomB[4] / 6.0
        StrMPoleSymplectic4Pass!(r_in, ele.len, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
            ele.FringeQuadEntrance, ele.FringeQuadExit, ele.FringeIntM0, ele.FringeIntP0, 
            ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle)
    end
    return nothing 
end

function pass_TPSA!(ele::KOCT, r_in::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}) where {T, TPS_Dim, Max_TPS_Degree}
    # ele: KOCT
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    PolynomB = zeros(4)
    # if ele.rad != 0
    #     println("Synchrtron radiation is not implemented for OCT in TPSA")
    # end
    if ele.PolynomB[1] == 0.0 && ele.PolynomB[2] == 0.0 && ele.PolynomB[3] == 0.0 && ele.PolynomB[4] == 0.0
        PolynomB[4] = ele.k3 / 6.0
        StrMPoleSymplectic4Pass!(r_in, ele.len, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
            ele.FringeQuadEntrance, ele.FringeQuadExit, ele.FringeIntM0, ele.FringeIntP0, 
            ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle)
    else
        PolynomB[1] = ele.PolynomB[1]
        PolynomB[2] = ele.PolynomB[2] 
        PolynomB[3] = ele.PolynomB[3] / 2.0
        PolynomB[4] = ele.PolynomB[4] / 6.0
        StrMPoleSymplectic4Pass!(r_in, ele.len, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
            ele.FringeQuadEntrance, ele.FringeQuadExit, ele.FringeIntM0, ele.FringeIntP0, 
            ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle)
    end
    return nothing
end
