# include("drift.jl")

function StrB2perp(bx::Float64, by::Float64, x::Float64, xpr::Float64, y::Float64, ypr::Float64)
    # Calculates sqr(|B x e|) , where e is a unit vector in the direction of velocity
    # v_norm2 = 1.0 / (1.0 + xpr^2 + ypr^2)
    # return (by^2 + bx^2 + (bx*ypr - by*xpr)^2) * v_norm2
    return bx*bx + by*by + (bx*xpr - by*ypr)^2
end

function strthinkickrad!(r::AbstractVector{Float64}, A::AbstractVector{Float64}, B::AbstractVector{Float64}, 
    L::Float64, E0::Float64, max_order::Int, rad_const::Float64)
    # Modified based on AT function. Ref[Terebilo, Andrei. "Accelerator modeling with MATLAB accelerator toolbox." PACS2001 (2001)].
    ReSum = B[max_order + 1]
    ImSum = A[max_order + 1]
    ReSumTemp = 0.0

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

    dp_0 = r[6]
    r[6] -= rad_const * B2P * factor

    # momentums after losing energy
    p_norm = 1.0 / (1.0 + r[6])

    r[2] = xpr / p_norm
    r[4] = ypr / p_norm

    r[2] -= L * ReSum
    r[4] += L * ImSum
    return nothing
end

function strthinkick!(r::AbstractVector{Float64}, A::AbstractVector{Float64}, B::AbstractVector{Float64}, L::Float64, max_order::Int)
    # Modified based on AT function. Ref[Terebilo, Andrei. "Accelerator modeling with MATLAB accelerator toolbox." PACS2001 (2001)].

    ReSum = B[max_order + 1]
    ImSum = A[max_order + 1]
    ReSumTemp = 0.0

    for i in max_order-1: -1: 0
        ReSumTemp = ReSum * r[1] - ImSum * r[3] + B[i+1]
        ImSum = ImSum * r[1] + ReSum * r[3] + A[i+1]
        ReSum = ReSumTemp
    end

    r[2] -= L * ReSum
    r[4] += L * ImSum
    return nothing
end


function StrMPoleSymplectic4Pass!(r::Array{Float64,1}, le::Float64, beti::Float64, A::Array{Float64,1}, B::Array{Float64,1}, 
    max_order::Int, num_int_step::Int, 
    FringeQuadEntrance::Int, FringeQuadExit::Int, #(no fringe), 1 (lee-whiting) or 2 (lee-whiting+elegant-like) 
    T1::Array{Float64,1}, T2::Array{Float64,1}, R1::Array{Float64,2}, R2::Array{Float64,2}, 
    RApertures::Array{Float64,1}, EApertures::Array{Float64,1}, KickAngle::Array{Float64,1}, 
    num_particles::Int, lost_flags::Array{Int64,1})
    # Modified based on AT function. Ref[Terebilo, Andrei. "Accelerator modeling with MATLAB accelerator toolbox." PACS2001 (2001)].

    DRIFT1  =  0.6756035959798286638
    DRIFT2 = -0.1756035959798286639
    KICK1 = 1.351207191959657328
    KICK2 = -1.702414383919314656

    SL = le/num_int_step
    L1 = SL*DRIFT1
    L2 = SL*DRIFT2
    K1 = SL*KICK1
    K2 = SL*KICK2

    if le > 0
        B[1] -= sin(KickAngle[1])/le
        A[1] += sin(KickAngle[2])/le
    end
    for c in 1:num_particles
        if lost_flags[c] == 1
            continue
        end
        r6 = @view r[(c-1)*6+1:c*6]
        if !isnan(r6[1])
            # Misalignment at entrance
            if !iszero(T1)
                addvv!(r6, T1)
            end
            if !iszero(R1)
                multmv!(r6, R1)
            end

            if FringeQuadEntrance != 0 
                multipole_fringe!(r6, le, A, B, max_order, 1.0, 1, beti)
            end

            # Integrator
            for m in 1:num_int_step
                drift6!(r6, L1, beti)
                strthinkick!(r6, A, B, K1, max_order)
                drift6!(r6, L2, beti)
                strthinkick!(r6, A, B, K2, max_order)
                drift6!(r6, L2, beti)
                strthinkick!(r6, A, B, K1, max_order)
                drift6!(r6, L1, beti)
            end

            if FringeQuadExit != 0
                multipole_fringe!(r6, le, A, B, max_order, -1.0, 1, beti)
            end

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
    if le > 0
        B[1] += sin(KickAngle[1]) / le
        A[1] -= sin(KickAngle[2]) / le
    end

    return nothing
end

function StrMPoleSymplectic4RadPass!(r::Array{Float64,1}, le::Float64, rad_const::Float64, beti::Float64, A::Array{Float64,1}, B::Array{Float64,1}, 
    max_order::Int, num_int_step::Int, 
    FringeQuadEntrance::Int, FringeQuadExit::Int, #(no fringe), 1 (lee-whiting) or 2 (lee-whiting+elegant-like) 
    T1::Array{Float64,1}, T2::Array{Float64,1}, R1::Array{Float64,2}, R2::Array{Float64,2}, 
    RApertures::Array{Float64,1}, EApertures::Array{Float64,1}, KickAngle::Array{Float64,1}, E0::Float64,
    num_particles::Int, lost_flags::Array{Int64,1})
    # Modified based on AT function. Ref[Terebilo, Andrei. "Accelerator modeling with MATLAB accelerator toolbox." PACS2001 (2001)].

    DRIFT1  =  0.6756035959798286638
    DRIFT2 = -0.1756035959798286639
    KICK1 = 1.351207191959657328
    KICK2 = -1.702414383919314656

    SL = le/num_int_step
    L1 = SL*DRIFT1
    L2 = SL*DRIFT2
    K1 = SL*KICK1
    K2 = SL*KICK2

    if le > 0
        B[1] -= sin(KickAngle[1])/le
        A[1] += sin(KickAngle[2])/le
    end
    for c in 1:num_particles
        if lost_flags[c] == 1
            continue
        end
        r6 = @view r[(c-1)*6+1:c*6]
        if !isnan(r6[1])
            # Misalignment at entrance
            if !iszero(T1)
                addvv!(r6, T1)
            end
            if !iszero(R1)
                multmv!(r6, R1)
            end

            if FringeQuadEntrance != 0 
                multipole_fringe!(r6, le, A, B, max_order, 1.0, 1, beti)
            end

            # Integrator
            for m in 1:num_int_step
                drift6!(r6, L1, beti)
                strthinkickrad!(r6, A, B, K1, E0, max_order, rad_const)
                drift6!(r6, L2, beti)
                strthinkickrad!(r6, A, B, K2, E0, max_order, rad_const)
                drift6!(r6, L2, beti)
                strthinkickrad!(r6, A, B, K1, E0, max_order, rad_const)
                drift6!(r6, L1, beti)
            end

            if FringeQuadExit != 0
                multipole_fringe!(r6, le, A, B, max_order, -1.0, 1, beti)
            end

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
    if le > 0
        B[1] += sin(KickAngle[1]) / le
        A[1] -= sin(KickAngle[2]) / le
    end

    return nothing
end

function pass!(ele::KQUAD, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam{Float64})
    # ele: KQUAD
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    lost_flags = particles.lost_flag
    PolynomB = zeros(4)
    E0 = particles.energy
    rad_const = 0.0

    if use_exact_beti == 1
        beti = 1.0 / particles.beta
    else
        beti = 1.0 
    end
    if ele.PolynomB[1] == 0.0 && ele.PolynomB[2] == 0.0 && ele.PolynomB[3] == 0.0 && ele.PolynomB[4] == 0.0
        PolynomB[2] = ele.k1
        if ele.rad == 0
            StrMPoleSymplectic4Pass!(r_in, ele.len, beti, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit, ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, num_particles, lost_flags)
        else
            if particles.mass == m_e
                rad_const = RAD_CONST_E * particles.gamma^3
            elseif particles.mass == m_p
                rad_const = RAD_CONST_P * particles.gamma^3
            else
                rad_const = 0.0
                println("SR is not implemented for this particle mass.")
            end
            StrMPoleSymplectic4RadPass!(r_in, ele.len, rad_const, beti, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit, ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, E0, num_particles, lost_flags)
        end
    else
        PolynomB[1] = ele.PolynomB[1]
        PolynomB[2] = ele.PolynomB[2] 
        PolynomB[3] = ele.PolynomB[3] / 2.0
        PolynomB[4] = ele.PolynomB[4] / 6.0
        if ele.rad == 0
            StrMPoleSymplectic4Pass!(r_in, ele.len, beti, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit, 
                ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, num_particles, lost_flags)
        else
            if particles.mass == m_e
                rad_const = RAD_CONST_E * particles.gamma^3
            elseif particles.mass == m_p
                rad_const = RAD_CONST_P * particles.gamma^3
            else
                rad_const = 0.0
                println("SR is not implemented for this particle mass.")
            end
            StrMPoleSymplectic4RadPass!(r_in, ele.len, rad_const, beti, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit, 
                ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, E0, num_particles, lost_flags)
        end
    end
    return nothing
end

function pass!(ele::KSEXT, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam{Float64})
    # ele: KSEXT
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    rad_const = 0.0

    lost_flags = particles.lost_flag
    PolynomB = zeros(4)
    E0 = particles.energy
    if use_exact_beti == 1
        beti = 1.0 / particles.beta
    else
        beti = 1.0 
    end
    if ele.PolynomB[1] == 0.0 && ele.PolynomB[2] == 0.0 && ele.PolynomB[3] == 0.0 && ele.PolynomB[4] == 0.0
        PolynomB[3] = ele.k2 / 2.0
        if ele.rad == 0
            StrMPoleSymplectic4Pass!(r_in, ele.len, beti,ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit, ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, num_particles, lost_flags)
        else
            if particles.mass == m_e
                rad_const = RAD_CONST_E * particles.gamma^3
            elseif particles.mass == m_p
                rad_const = RAD_CONST_P * particles.gamma^3
            else
                rad_const = 0.0
                println("SR is not implemented for this particle mass.")
            end
            StrMPoleSymplectic4RadPass!(r_in, ele.len, rad_const, beti, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit, ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, E0, num_particles, lost_flags)
        end
    else
        PolynomB[1] = ele.PolynomB[1]
        PolynomB[2] = ele.PolynomB[2] 
        PolynomB[3] = ele.PolynomB[3] / 2.0
        PolynomB[4] = ele.PolynomB[4] / 6.0
        if ele.rad == 0
            StrMPoleSymplectic4Pass!(r_in, ele.len, beti, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit, ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, num_particles, lost_flags)
        else
            if particles.mass == m_e
                rad_const = RAD_CONST_E * particles.gamma^3
            elseif particles.mass == m_p
                rad_const = RAD_CONST_P * particles.gamma^3
            else
                rad_const = 0.0
                println("SR is not implemented for this particle mass.")
            end
            StrMPoleSymplectic4RadPass!(r_in, ele.len, rad_const, beti, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit, ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, E0, num_particles, lost_flags)
        end
    end
    return nothing
end

function pass!(ele::KOCT, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam{Float64})
    # ele: KOCT
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    rad_const = 0.0

    lost_flags = particles.lost_flag
    PolynomB = zeros(4)
    E0 = particles.energy
    if use_exact_beti == 1
        beti = 1.0 / particles.beta
    else
        beti = 1.0 
    end
    if ele.PolynomB[1] == 0.0 && ele.PolynomB[2] == 0.0 && ele.PolynomB[3] == 0.0 && ele.PolynomB[4] == 0.0
        PolynomB[4] = ele.k3 / 6.0
        if ele.rad == 0
            StrMPoleSymplectic4Pass!(r_in, ele.len, beti, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit, ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, num_particles, lost_flags)
        else
            if particles.mass == m_e
                rad_const = RAD_CONST_E * particles.gamma^3
            elseif particles.mass == m_p
                rad_const = RAD_CONST_P * particles.gamma^3
            else
                rad_const = 0.0
                println("SR is not implemented for this particle mass.")
            end
            StrMPoleSymplectic4RadPass!(r_in, ele.len, beti, particles.mass, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit, ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, E0, num_particles, lost_flags)
        end
    else
        PolynomB[1] = ele.PolynomB[1]
        PolynomB[2] = ele.PolynomB[2] 
        PolynomB[3] = ele.PolynomB[3] / 2.0
        PolynomB[4] = ele.PolynomB[4] / 6.0
        if ele.rad == 0
            StrMPoleSymplectic4Pass!(r_in, ele.len, beti, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit, ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, num_particles, lost_flags)
        else
            if particles.mass == m_e
                rad_const = RAD_CONST_E * particles.gamma^3
            elseif particles.mass == m_p
                rad_const = RAD_CONST_P * particles.gamma^3
            else
                rad_const = 0.0
                println("SR is not implemented for this particle mass.")
            end
            StrMPoleSymplectic4RadPass!(r_in, ele.len, rad_const, beti, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit, ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, E0, num_particles, lost_flags)
        end
    end
    return nothing
end


###################
# multi-threading
function StrMPoleSymplectic4Pass_P!(r::Array{Float64,1}, le::Float64, beti::Float64, A::Array{Float64,1}, B::Array{Float64,1}, 
    max_order::Int, num_int_step::Int, 
    FringeQuadEntrance::Int, FringeQuadExit::Int, #(no fringe), 1 (lee-whiting) or 2 (lee-whiting+elegant-like) 
    T1::Array{Float64,1}, T2::Array{Float64,1}, R1::Array{Float64,2}, R2::Array{Float64,2}, 
    RApertures::Array{Float64,1}, EApertures::Array{Float64,1}, KickAngle::Array{Float64,1}, 
    num_particles::Int, lost_flags::Array{Int64,1})

    DRIFT1  =  0.6756035959798286638
    DRIFT2 = -0.1756035959798286639
    KICK1 = 1.351207191959657328
    KICK2 = -1.702414383919314656

    SL = le/num_int_step
    L1 = SL*DRIFT1
    L2 = SL*DRIFT2
    K1 = SL*KICK1
    K2 = SL*KICK2

    B[1] -= sin(KickAngle[1])/le
    A[1] += sin(KickAngle[2])/le

    Threads.@threads for c in 1:num_particles
        if lost_flags[c] == 1
            continue
        end
        r6 = @view r[(c-1)*6+1:c*6]
        if !isnan(r6[1])
            # Misalignment at entrance
            if !iszero(T1)
                addvv!(r6, T1)
            end
            if !iszero(R1)
                multmv!(r6, R1)
            end


            if FringeQuadEntrance != 0 
                multipole_fringe!(r6, le, A, B, max_order, 1.0, 1, beti)
            end


            # Integrator
            for m in 1:num_int_step
                drift6!(r6, L1, beti)
                strthinkick!(r6, A, B, K1, max_order)
                drift6!(r6, L2, beti)
                strthinkick!(r6, A, B, K2, max_order)
                drift6!(r6, L2, beti)
                strthinkick!(r6, A, B, K1, max_order)
                drift6!(r6, L1, beti)
            end

            if FringeQuadExit != 0
                multipole_fringe!(r6, le, A, B, max_order, -1.0, 1, beti)
            end

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

    B[1] += sin(KickAngle[1]) / le
    A[1] -= sin(KickAngle[2]) / le
    return nothing
end

function StrMPoleSymplectic4RadPass_P!(r::Array{Float64,1}, le::Float64, rad_const::Float64, beti::Float64, A::Array{Float64,1}, B::Array{Float64,1}, 
    max_order::Int, num_int_step::Int, 
    FringeQuadEntrance::Int, FringeQuadExit::Int, #(no fringe), 1 (lee-whiting) or 2 (lee-whiting+elegant-like) 
    T1::Array{Float64,1}, T2::Array{Float64,1}, R1::Array{Float64,2}, R2::Array{Float64,2}, 
    RApertures::Array{Float64,1}, EApertures::Array{Float64,1}, KickAngle::Array{Float64,1}, E0::Float64,
    num_particles::Int, lost_flags::Array{Int64,1})

    DRIFT1  =  0.6756035959798286638
    DRIFT2 = -0.1756035959798286639
    KICK1 = 1.351207191959657328
    KICK2 = -1.702414383919314656

    SL = le/num_int_step
    L1 = SL*DRIFT1
    L2 = SL*DRIFT2
    K1 = SL*KICK1
    K2 = SL*KICK2

    if le > 0
        B[1] -= sin(KickAngle[1])/le
        A[1] += sin(KickAngle[2])/le
    end

    Threads.@threads for c in 1:num_particles
    # for c in 1:num_particles
        if lost_flags[c] == 1
            continue
        end
        r6 = @view r[(c-1)*6+1:c*6]
        if !isnan(r6[1])
            # Misalignment at entrance
            if !iszero(T1)
                addvv!(r6, T1)
            end
            if !iszero(R1)
                multmv!(r6, R1)
            end

            # Check physical apertures at the entrance of the magnet
            # if RApertures != nothing
            #     checkiflostRectangularAp(r6, RApertures)
            # end
            # if EApertures != nothing
            #     checkiflostEllipticalAp(r6, EApertures)
            # end


            if FringeQuadEntrance != 0 
                multipole_fringe!(r6, le, A, B, max_order, 1.0, 1, beti)
            end

            # Integrator
            for m in 1:num_int_step
                drift6!(r6, L1, beti)
                strthinkickrad!(r6, A, B, K1, E0, max_order, rad_const)
                drift6!(r6, L2, beti)
                strthinkickrad!(r6, A, B, K2, E0, max_order, rad_const)
                drift6!(r6, L2, beti)
                strthinkickrad!(r6, A, B, K1, E0, max_order, rad_const)
                drift6!(r6, L1, beti)
            end

            if FringeQuadExit != 0
                multipole_fringe!(r6, le, A, B, max_order, -1.0, 1, beti)
            end

            # Check physical apertures at the exit of the magnet
            # if RApertures != nothing
            #     checkiflostRectangularAp(r6, RApertures)
            # end
            # if EApertures != nothing
            #     checkiflostEllipticalAp(r6, EApertures)
            # end

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

    if le > 0
        B[1] += sin(KickAngle[1]) / le
        A[1] -= sin(KickAngle[2]) / le
    end
    return nothing
end

function pass_P!(ele::KQUAD, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam{Float64})
    # ele: KQUAD
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    rad_const = 0.0

    lost_flags = particles.lost_flag
    PolynomB = zeros(4)
    E0 = particles.energy
    if use_exact_beti == 1
        beti = 1.0 / particles.beta
    else
        beti = 1.0 
    end
    if ele.PolynomB[1] == 0.0 && ele.PolynomB[2] == 0.0 && ele.PolynomB[3] == 0.0 && ele.PolynomB[4] == 0.0
        PolynomB[2] = ele.k1
        if ele.rad == 0
            StrMPoleSymplectic4Pass_P!(r_in, ele.len, beti, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit, ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, num_particles, lost_flags)
        else
            if particles.mass == m_e
                rad_const = RAD_CONST_E * particles.gamma^3
            elseif particles.mass == m_p
                rad_const = RAD_CONST_P * particles.gamma^3
            else
                rad_const = 0.0
                println("SR is not implemented for this particle mass.")
            end
            StrMPoleSymplectic4RadPass_P!(r_in, ele.len, rad_const, beti, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit, ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, E0, num_particles, lost_flags)
        end
    else
        PolynomB[1] = ele.PolynomB[1]
        PolynomB[2] = ele.PolynomB[2] 
        PolynomB[3] = ele.PolynomB[3] / 2.0
        PolynomB[4] = ele.PolynomB[4] / 6.0
        if ele.rad == 0
            StrMPoleSymplectic4Pass_P!(r_in, ele.len, beti, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit, 
                ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, num_particles, lost_flags)
        else
            if particles.mass == m_e
                rad_const = RAD_CONST_E * particles.gamma^3
            elseif particles.mass == m_p
                rad_const = RAD_CONST_P * particles.gamma^3
            else
                rad_const = 0.0
                println("SR is not implemented for this particle mass.")
            end
            StrMPoleSymplectic4RadPass_P!(r_in, ele.len, rad_const, beti, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit, 
                ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, E0, num_particles, lost_flags)
        end
    end
    return nothing
end

function pass_P!(ele::KSEXT, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam{Float64})
    # ele: KSEXT
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    rad_const = 0.0

    lost_flags = particles.lost_flag
    PolynomB = zeros(4)
    E0 = particles.energy
    if use_exact_beti == 1
        beti = 1.0 / particles.beta
    else
        beti = 1.0 
    end
    if ele.PolynomB[1] == 0.0 && ele.PolynomB[2] == 0.0 && ele.PolynomB[3] == 0.0 && ele.PolynomB[4] == 0.0
        PolynomB[3] = ele.k2 / 2.0
        if ele.rad == 0
            StrMPoleSymplectic4Pass_P!(r_in, ele.len, beti, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit, 
                ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, num_particles, lost_flags)
        else
            if particles.mass == m_e
                rad_const = RAD_CONST_E * particles.gamma^3
            elseif particles.mass == m_p
                rad_const = RAD_CONST_P * particles.gamma^3
            else
                rad_const = 0.0
                println("SR is not implemented for this particle mass.")
            end
            StrMPoleSymplectic4RadPass_P!(r_in, ele.len, rad_const, beti, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit, 
                ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, E0, num_particles, lost_flags)
        end
    else
        PolynomB[1] = ele.PolynomB[1]
        PolynomB[2] = ele.PolynomB[2]
        PolynomB[3] = ele.PolynomB[3] / 2.0
        PolynomB[4] = ele.PolynomB[4] / 6.0
        if ele.rad == 0
            StrMPoleSymplectic4Pass_P!(r_in, ele.len, beti, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit, 
                ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, num_particles, lost_flags)
        else
            if particles.mass == m_e
                rad_const = RAD_CONST_E * particles.gamma^3
            elseif particles.mass == m_p
                rad_const = RAD_CONST_P * particles.gamma^3
            else
                rad_const = 0.0
                println("SR is not implemented for this particle mass.")
            end
            StrMPoleSymplectic4RadPass_P!(r_in, ele.len, rad_const, beti, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit, 
                ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, E0, num_particles, lost_flags)
        end
    end
    return nothing
end

function pass_P!(ele::KOCT, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam{Float64})
    # ele: KOCT
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    rad_const = 0.0

    lost_flags = particles.lost_flag
    PolynomB = zeros(4)
    E0 = particles.energy
    if use_exact_beti == 1
        beti = 1.0 / particles.beta
    else
        beti = 1.0 
    end
    if ele.PolynomB[1] == 0.0 && ele.PolynomB[2] == 0.0 && ele.PolynomB[3] == 0.0 && ele.PolynomB[4] == 0.0
        PolynomB[4] = ele.k3 / 6.0
        if ele.rad == 0
            StrMPoleSymplectic4Pass_P!(r_in, ele.len, beti, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit, 
                ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, num_particles, lost_flags)
        else
            if particles.mass == m_e
                rad_const = RAD_CONST_E * particles.gamma^3
            elseif particles.mass == m_p
                rad_const = RAD_CONST_P * particles.gamma^3
            else
                rad_const = 0.0
                println("SR is not implemented for this particle mass.")
            end
            StrMPoleSymplectic4RadPass_P!(r_in, ele.len, rad_const, beti, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit, 
                ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, E0, num_particles, lost_flags)
        end
    else
        PolynomB[1] = ele.PolynomB[1]
        PolynomB[2] = ele.PolynomB[2] 
        PolynomB[3] = ele.PolynomB[3] / 2.0
        PolynomB[4] = ele.PolynomB[4] / 6.0
        if ele.rad == 0
            StrMPoleSymplectic4Pass_P!(r_in, ele.len, beti, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit, ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, num_particles, lost_flags)
        else
            if particles.mass == m_e
                rad_const = RAD_CONST_E * particles.gamma^3
            elseif particles.mass == m_p
                rad_const = RAD_CONST_P * particles.gamma^3
            else
                rad_const = 0.0
                println("SR is not implemented for this particle mass.")
            end
            StrMPoleSymplectic4RadPass_P!(r_in, ele.len, rad_const, beti,ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit, ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, E0, num_particles, lost_flags)
        end
    end
    return nothing
end