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


function StrMPoleSymplectic4Pass!(r::Matrix{Float64}, le::Float64, beti::Float64, A::Array{Float64,1}, B::Array{Float64,1}, 
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
        r6 = @view r[c, :]
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

function StrMPoleSymplectic4RadPass!(r::Matrix{Float64}, le::Float64, rad_const::Float64, beti::Float64, A::Array{Float64,1}, B::Array{Float64,1}, 
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
        r6 = @view r[c, :]
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

function pass!(ele::KQUAD, r_in::Matrix{Float64}, num_particles::Int64, particles::Beam{Float64})
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

    PolynomB[1] = ele.k0
    PolynomB[2] = ele.k1 
    PolynomB[3] = ele.k2 / 2.0
    PolynomB[4] = ele.k3 / 6.0
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
    return nothing
end

function pass!(ele::KSEXT, r_in::Matrix{Float64}, num_particles::Int64, particles::Beam{Float64})
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

    PolynomB[1] = ele.k0
    PolynomB[2] = ele.k1 
    PolynomB[3] = ele.k2 / 2.0
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
        StrMPoleSymplectic4RadPass!(r_in, ele.len, rad_const, beti, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
            ele.FringeQuadEntrance, ele.FringeQuadExit, ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, E0, num_particles, lost_flags)
    end
    return nothing
end

function pass!(ele::KOCT, r_in::Matrix{Float64}, num_particles::Int64, particles::Beam{Float64})
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

    PolynomB[1] = ele.k0
    PolynomB[2] = ele.k1 
    PolynomB[3] = ele.k2 / 2.0
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
        StrMPoleSymplectic4RadPass!(r_in, ele.len, rad_const, beti, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
            ele.FringeQuadEntrance, ele.FringeQuadExit, ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, E0, num_particles, lost_flags)
    end
    return nothing
end


###################
# multi-threading
function StrMPoleSymplectic4Pass_P!(r::Matrix{Float64}, le::Float64, beti::Float64, A::Array{Float64,1}, B::Array{Float64,1}, 
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
        r6 = @view r[c, :]
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

function StrMPoleSymplectic4RadPass_P!(r::Matrix{Float64}, le::Float64, rad_const::Float64, beti::Float64, A::Array{Float64,1}, B::Array{Float64,1}, 
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
        r6 = @view r[c, :]
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

function pass_P!(ele::KQUAD, r_in::Matrix{Float64}, num_particles::Int64, particles::Beam{Float64})
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

    PolynomB[1] = ele.k0
    PolynomB[2] = ele.k1
    PolynomB[3] = ele.k2 / 2.0
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
    return nothing
end

function pass_P!(ele::KSEXT, r_in::Matrix{Float64}, num_particles::Int64, particles::Beam{Float64})
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

    PolynomB[1] = ele.k0
    PolynomB[2] = ele.k1
    PolynomB[3] = ele.k2 / 2.0
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
    return nothing
end

function pass_P!(ele::KOCT, r_in::Matrix{Float64}, num_particles::Int64, particles::Beam{Float64})
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

    PolynomB[1] = ele.k0
    PolynomB[2] = ele.k1 
    PolynomB[3] = ele.k2 / 2.0
    PolynomB[4] = ele.k3 / 6.0
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
    return nothing
end


########################################################
# Space charge
function StrMPoleSymplectic4Pass_SC!(r::Matrix{Float64}, le::Float64, beti::Float64, A::Array{Float64,1}, B::Array{Float64,1}, 
    max_order::Int, num_int_step::Int, 
    FringeQuadEntrance::Int, FringeQuadExit::Int, #(no fringe), 1 (lee-whiting) or 2 (lee-whiting+elegant-like) 
    T1::Array{Float64,1}, T2::Array{Float64,1}, R1::Array{Float64,2}, R2::Array{Float64,2}, 
    RApertures::Array{Float64,1}, EApertures::Array{Float64,1}, KickAngle::Array{Float64,1}, 
    num_particles::Int, lost_flags::Array{Int64,1}, a::Float64, b::Float64, Nl::Int, Nm::Int, K::Float64, Nsteps::Int)
    # Ref[Terebilo, Andrei. "Accelerator modeling with MATLAB accelerator toolbox." PACS2001 (2001)].
    # and [Qiang, Ji. "Differentiable self-consistent space-charge simulation for accelerator design." Physical Review Accelerators and Beams 26.2 (2023): 024601.]

    DRIFT1  =  0.6756035959798286638
    DRIFT2 = -0.1756035959798286639
    KICK1 = 1.351207191959657328
    KICK2 = -1.702414383919314656

    lstep = le / Nsteps
    SL = lstep / 2.0 /num_int_step
    L1 = SL*DRIFT1
    L2 = SL*DRIFT2
    K1 = SL*KICK1
    K2 = SL*KICK2

    if le > 0
        B[1] -= sin(KickAngle[1])/le
        A[1] += sin(KickAngle[2])/le 
    end

    for step in 1:Nsteps
        for c in 1:num_particles
            if lost_flags[c] == 1
                continue
            end
            r6 = @view r[c, :]
            
            if step == 1
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

            if check_lost(r6)
                lost_flags[c] = 1
            end
        end

        space_charge!(r, K, Nl, Nm, a/Nl, b/Nm, a, b, num_particles, lstep, lost_flags)

        for c in 1:num_particles
            if lost_flags[c] == 1
                continue
            end
            r6 = @view r[c, :]

            NormL1 = L1 / (1.0 + r6[6])
            NormL2 = L2 / (1.0 + r6[6])
    
            # Integrator
            for m in 1:num_int_step
                drift6!(r6, L1,beti)
                strthinkick!(r6, A, B, K1, max_order)
                drift6!(r6, L2, beti)
                strthinkick!(r6, A, B, K2, max_order)
                drift6!(r6, L2, beti)
                strthinkick!(r6, A, B, K1, max_order)
                drift6!(r6, L1, beti)
            end
            
            if step == Nsteps
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

function StrMPoleSymplectic4RadPass_SC!(r::Matrix{Float64}, le::Float64, rad_const::Float64, beti::Float64, A::Array{Float64,1}, B::Array{Float64,1}, 
    max_order::Int, num_int_step::Int, 
    FringeQuadEntrance::Int, FringeQuadExit::Int, #(no fringe), 1 (lee-whiting) or 2 (lee-whiting+elegant-like) 
    T1::Array{Float64,1}, T2::Array{Float64,1}, R1::Array{Float64,2}, R2::Array{Float64,2}, 
    RApertures::Array{Float64,1}, EApertures::Array{Float64,1}, KickAngle::Array{Float64,1}, E0::Float64,
    num_particles::Int, lost_flags::Array{Int64,1}, a::Float64, b::Float64, Nl::Int, Nm::Int, K::Float64, Nsteps::Int)
    # Ref[Terebilo, Andrei. "Accelerator modeling with MATLAB accelerator toolbox." PACS2001 (2001)].
    # and [Qiang, Ji. "Differentiable self-consistent space-charge simulation for accelerator design." Physical Review Accelerators and Beams 26.2 (2023): 024601.]

    DRIFT1  =  0.6756035959798286638
    DRIFT2 = -0.1756035959798286639
    KICK1 = 1.351207191959657328
    KICK2 = -1.702414383919314656

    lstep = le / Nsteps
    SL = lstep / 2.0 /num_int_step
    L1 = SL*DRIFT1
    L2 = SL*DRIFT2
    K1 = SL*KICK1
    K2 = SL*KICK2

    if le > 0
        B[1] -= sin(KickAngle[1])/ le 
        A[1] += sin(KickAngle[2])/ le
    end

    for step in 1:Nsteps
        for c in 1:num_particles
            if lost_flags[c] == 1
                continue
            end
            r6 = @view r[c, :]
            if step == 1
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
            if check_lost(r6)
                lost_flags[c] = 1
            end
        end

        space_charge!(r, K, Nl, Nm, a/Nl, b/Nm, a, b, num_particles, le, lost_flags)

        for c in 1:num_particles
            if lost_flags[c] == 1
                continue
            end
            r6 = @view r[c, :]
    
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
            
            if step == Nsteps
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

function pass!(ele::KQUAD_SC, r_in::Matrix{Float64}, num_particles::Int64, particles::Beam{Float64})
    # ele: KQUAD_SC
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    lost_flags = particles.lost_flag
    K = calculate_K(particles, particles.current)
    PolynomB = zeros(4)
    E0 = particles.energy
    if use_exact_beti == 1
        beti = 1.0 / particles.beta
    else
        beti = 1.0 
    end

    PolynomB[1] = ele.k0
    PolynomB[2] = ele.k1 
    PolynomB[3] = ele.k2 / 2.0
    PolynomB[4] = ele.k3 / 6.0
    if ele.rad == 0
        StrMPoleSymplectic4Pass_SC!(r_in,  ele.len, beti, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit, 
                ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, num_particles, lost_flags,
                ele.a, ele.b, ele.Nl, ele.Nm, K, ele.Nsteps)
    else
        if particles.mass == m_e
            rad_const = RAD_CONST_E * particles.gamma^3
        elseif particles.mass == m_p
            rad_const = RAD_CONST_P * particles.gamma^3
        else
            rad_const = 0.0
            println("SR is not implemented for this particle mass.")
        end
        StrMPoleSymplectic4RadPass_SC!(r_in,  ele.len, rad_const, beti, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit, 
                ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, E0, num_particles, lost_flags,
                ele.a, ele.b, ele.Nl, ele.Nm, K, ele.Nsteps)
    end
    return nothing
end

function pass!(ele::KSEXT_SC, r_in::Matrix{Float64}, num_particles::Int64, particles::Beam{Float64})
    # ele: KSEXT
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    lost_flags = particles.lost_flag
    K = calculate_K(particles, particles.current)
    PolynomB = zeros(4)
    E0 = particles.energy
    if use_exact_beti == 1
        beti = 1.0 / particles.beta
    else
        beti = 1.0 
    end

    PolynomB[1] = ele.k0
    PolynomB[2] = ele.k1 
    PolynomB[3] = ele.k2 / 2.0
    PolynomB[4] = ele.k3 / 6.0
    if ele.rad == 0
        StrMPoleSymplectic4Pass_SC!(r_in, ele.len, beti, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit, 
                ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, num_particles, lost_flags,
                ele.a, ele.b, ele.Nl, ele.Nm, K, ele.Nsteps)
    else
        if particles.mass == m_e
            rad_const = RAD_CONST_E * particles.gamma^3
        elseif particles.mass == m_p
            rad_const = RAD_CONST_P * particles.gamma^3
        else
            rad_const = 0.0
            println("SR is not implemented for this particle mass.")
        end
        StrMPoleSymplectic4RadPass_SC!(r_in, ele.len, rad_const, beti, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit, 
                ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, E0, num_particles, lost_flags,
                ele.a, ele.b, ele.Nl, ele.Nm, K, ele.Nsteps)
    end
    
    return nothing
end

function pass!(ele::KOCT_SC, r_in::Matrix{Float64}, num_particles::Int64, particles::Beam{Float64})
    # ele: KOCT_SC
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    lost_flags = particles.lost_flag
    PolynomB = zeros(4)
    K = calculate_K(particles, particles.current)
    E0 = particles.energy
    if use_exact_beti == 1
        beti = 1.0 / particles.beta
    else
        beti = 1.0 
    end

    PolynomB[1] = ele.k0
    PolynomB[2] = ele.k1 
    PolynomB[3] = ele.k2 / 2.0
    PolynomB[4] = ele.k3 / 6.0
    if ele.rad == 0
        StrMPoleSymplectic4Pass_SC!(r_in, ele.len, beti, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit, 
                ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, num_particles, lost_flags,
                ele.a, ele.b, ele.Nl, ele.Nm, K, ele.Nsteps)
    else
        if particles.mass == m_e
            rad_const = RAD_CONST_E * particles.gamma^3
        elseif particles.mass == m_p
            rad_const = RAD_CONST_P * particles.gamma^3
        else
            rad_const = 0.0
            println("SR is not implemented for this particle mass.")
        end
        StrMPoleSymplectic4RadPass_SC!(r_in, ele.len, rad_const, beti, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit,
                ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, E0, num_particles, lost_flags,
                ele.a, ele.b, ele.Nl, ele.Nm, K, ele.Nsteps)
    end
    return nothing
end


function StrMPoleSymplectic4Pass_SC!(r::Matrix{DTPSAD{N, T}}, le::DTPSAD{N, T}, beti::Float64, A::Array{DTPSAD{N, T},1}, B::Array{DTPSAD{N, T},1}, 
    max_order::Int, num_int_step::Int, 
    FringeQuadEntrance::Int, FringeQuadExit::Int, #(no fringe), 1 (lee-whiting) or 2 (lee-whiting+elegant-like) 
    T1::Array{DTPSAD{N, T},1}, T2::Array{DTPSAD{N, T},1}, R1::Array{DTPSAD{N, T},2}, R2::Array{DTPSAD{N, T},2}, 
    RApertures::Array{Float64,1}, EApertures::Array{Float64,1}, KickAngle::Array{DTPSAD{N, T},1}, 
    num_particles::Int, lost_flags::Array{Int64,1}, a::DTPSAD{N, T}, b::DTPSAD{N, T}, Nl::Int, Nm::Int, K::DTPSAD{N, T}, Nsteps::Int) where {N, T}
    # Ref[Terebilo, Andrei. "Accelerator modeling with MATLAB accelerator toolbox." PACS2001 (2001)].
    # and [Qiang, Ji. "Differentiable self-consistent space-charge simulation for accelerator design." Physical Review Accelerators and Beams 26.2 (2023): 024601.]

    DRIFT1  =  0.6756035959798286638
    DRIFT2 = -0.1756035959798286639
    KICK1 = 1.351207191959657328
    KICK2 = -1.702414383919314656

    lstep = le / Nsteps
    SL = lstep / 2.0 /num_int_step
    L1 = SL*DRIFT1
    L2 = SL*DRIFT2
    K1 = SL*KICK1
    K2 = SL*KICK2

    if le > 0
        B[1] -= sin(KickAngle[1])/le
        A[1] += sin(KickAngle[2])/le 
    end

    for step in 1:Nsteps
        @inbounds for c in 1:num_particles
            if lost_flags[c] == 1
                continue
            end
            r6 = @view r[c, :]
            
            if step == 1
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
            end

            # Integrator
            for m in 1:num_int_step
                drift6!(r6, L1)
                strthinkick!(r6, A, B, K1, max_order)
                drift6!(r6, L2)
                strthinkick!(r6, A, B, K2, max_order)
                drift6!(r6, L2)
                strthinkick!(r6, A, B, K1, max_order)
                drift6!(r6, L1)
            end

            if check_lost(r6)
                lost_flags[c] = 1
            end
        end

        space_charge!(r, K, Nl, Nm, a/Nl, b/Nm, a, b, num_particles, lstep, lost_flags)

        for c in 1:num_particles
            if lost_flags[c] == 1
                continue
            end
            r6 = @view r[c, :]

            NormL1 = L1 / (1.0 + r6[6])
            NormL2 = L2 / (1.0 + r6[6])
    
            # Integrator
            for m in 1:num_int_step
                drift6!(r6, L1)
                strthinkick!(r6, A, B, K1, max_order)
                drift6!(r6, L2)
                strthinkick!(r6, A, B, K2, max_order)
                drift6!(r6, L2)
                strthinkick!(r6, A, B, K1, max_order)
                drift6!(r6, L1)
            end
            
            if step == Nsteps
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

function StrMPoleSymplectic4RadPass_SC!(r::Matrix{DTPSAD{N, T}}, le::DTPSAD{N, T}, rad_const::DTPSAD{N, T}, beti::Float64, A::Array{DTPSAD{N, T},1}, B::Array{DTPSAD{N, T},1}, 
    max_order::Int, num_int_step::Int, 
    FringeQuadEntrance::Int, FringeQuadExit::Int, #(no fringe), 1 (lee-whiting) or 2 (lee-whiting+elegant-like) 
    T1::Array{DTPSAD{N, T},1}, T2::Array{DTPSAD{N, T},1}, R1::Array{DTPSAD{N, T},2}, R2::Array{DTPSAD{N, T},2}, 
    RApertures::Array{Float64,1}, EApertures::Array{Float64,1}, KickAngle::Array{DTPSAD{N, T},1}, E0::DTPSAD{N, T},
    num_particles::Int, lost_flags::Array{Int64,1}, a::DTPSAD{N, T}, b::DTPSAD{N, T}, Nl::Int, Nm::Int, K::DTPSAD{N, T}, Nsteps::Int) where {N, T}
    # Ref[Terebilo, Andrei. "Accelerator modeling with MATLAB accelerator toolbox." PACS2001 (2001)].
    # and [Qiang, Ji. "Differentiable self-consistent space-charge simulation for accelerator design." Physical Review Accelerators and Beams 26.2 (2023): 024601.]

    DRIFT1  =  0.6756035959798286638
    DRIFT2 = -0.1756035959798286639
    KICK1 = 1.351207191959657328
    KICK2 = -1.702414383919314656

    lstep = le / Nsteps
    SL = lstep / 2.0 /num_int_step
    L1 = SL*DRIFT1
    L2 = SL*DRIFT2
    K1 = SL*KICK1
    K2 = SL*KICK2

    if le > 0
        B[1] -= sin(KickAngle[1])/ le 
        A[1] += sin(KickAngle[2])/ le
    end

    for step in 1:Nsteps
        @inbounds for c in 1:num_particles
            if lost_flags[c] == 1
                continue
            end
            r6 = @view r[c, :]
            if step == 1
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
            end
            # Integrator
            for m in 1:num_int_step
                drift6!(r6, L1)
                strthinkickrad!(r6, A, B, K1, E0, max_order, rad_const)
                drift6!(r6, L2)
                strthinkickrad!(r6, A, B, K2, E0, max_order, rad_const)
                drift6!(r6, L2)
                strthinkickrad!(r6, A, B, K1, E0, max_order, rad_const)
                drift6!(r6, L1)
            end
            if check_lost(r6)
                lost_flags[c] = 1
            end
        end

        space_charge!(r, K, Nl, Nm, a/Nl, b/Nm, a, b, num_particles, le, lost_flags)

        for c in 1:num_particles
            if lost_flags[c] == 1
                continue
            end
            r6 = @view r[c, :]
    
            # Integrator
            for m in 1:num_int_step
                drift6!(r6, L1)
                strthinkickrad!(r6, A, B, K1, E0, max_order, rad_const)
                drift6!(r6, L2)
                strthinkickrad!(r6, A, B, K2, E0, max_order, rad_const)
                drift6!(r6, L2)
                strthinkickrad!(r6, A, B, K1, E0, max_order, rad_const)
                drift6!(r6, L1)
            end
            
            if step == Nsteps
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

function pass!(ele::KQUAD_SC{DTPSAD{N, T}}, r_in::Matrix{DTPSAD{N, T}}, num_particles::Int64, particles::Beam{DTPSAD{N, T}}) where {N, T}
    # ele: KQUAD_SC
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    lost_flags = particles.lost_flag
    K = calculate_K(particles, particles.current)
    PolynomB = zeros(DTPSAD{N, T}, 4)
    E0 = particles.energy
    if use_exact_beti == 1
        beti = 1.0 / particles.beta.val
    else
        beti = 1.0 
    end

    PolynomB[1] = ele.k0
    PolynomB[2] = ele.k1 
    PolynomB[3] = ele.k2 / 2.0
    PolynomB[4] = ele.k3 / 6.0
    if ele.rad == 0
        StrMPoleSymplectic4Pass_SC!(r_in,  ele.len, beti, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit, 
                ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, num_particles, lost_flags,
                ele.a, ele.b, ele.Nl, ele.Nm, K, ele.Nsteps)
    else
        if particles.mass == m_e
            rad_const = RAD_CONST_E * particles.gamma^3
        elseif particles.mass == m_p
            rad_const = RAD_CONST_P * particles.gamma^3
        else
            rad_const = 0.0
            println("SR is not implemented for this particle mass.")
        end
        StrMPoleSymplectic4RadPass_SC!(r_in,  ele.len, rad_const, beti, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit, 
                ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, E0, num_particles, lost_flags,
                ele.a, ele.b, ele.Nl, ele.Nm, K, ele.Nsteps)
    end
    return nothing
end

function pass!(ele::KSEXT_SC{DTPSAD{N, T}}, r_in::Matrix{DTPSAD{N, T}}, num_particles::Int64, particles::Beam{DTPSAD{N, T}}) where {N, T}
    # ele: KSEXT
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    lost_flags = particles.lost_flag
    K = calculate_K(particles, particles.current)
    PolynomB = zeros(DTPSAD{N, T}, 4)
    E0 = particles.energy
    if use_exact_beti == 1
        beti = 1.0 / particles.beta.val
    else
        beti = 1.0 
    end

    PolynomB[1] = ele.k0
    PolynomB[2] = ele.k1 
    PolynomB[3] = ele.k2 / 2.0
    PolynomB[4] = ele.k3 / 6.0
    if ele.rad == 0
        StrMPoleSymplectic4Pass_SC!(r_in, ele.len, beti, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit, 
                ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, num_particles, lost_flags,
                ele.a, ele.b, ele.Nl, ele.Nm, K, ele.Nsteps)
    else
        if particles.mass == m_e
            rad_const = RAD_CONST_E * particles.gamma^3
        elseif particles.mass == m_p
            rad_const = RAD_CONST_P * particles.gamma^3
        else
            rad_const = 0.0
            println("SR is not implemented for this particle mass.")
        end
        StrMPoleSymplectic4RadPass_SC!(r_in, ele.len, rad_const, beti, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit, 
                ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, E0, num_particles, lost_flags,
                ele.a, ele.b, ele.Nl, ele.Nm, K, ele.Nsteps)
    end
    return nothing
end

function pass!(ele::KOCT_SC{DTPSAD{N, T}}, r_in::Matrix{DTPSAD{N, T}}, num_particles::Int64, particles::Beam{DTPSAD{N, T}}) where {N, T}
    # ele: KOCT_SC
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    lost_flags = particles.lost_flag
    PolynomB = zeros(DTPSAD{N, T}, 4)
    K = calculate_K(particles, particles.current)
    E0 = particles.energy
    if use_exact_beti == 1
        beti = 1.0 / particles.beta.val
    else
        beti = 1.0 
    end

    PolynomB[1] = ele.k0
    PolynomB[2] = ele.k1 
    PolynomB[3] = ele.k2 / 2.0
    PolynomB[4] = ele.k3 / 6.0
    if ele.rad == 0
        StrMPoleSymplectic4Pass_SC!(r_in, ele.len, beti, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit, 
                ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, num_particles, lost_flags,
                ele.a, ele.b, ele.Nl, ele.Nm, K, ele.Nsteps)
    else
        if particles.mass == m_e
            rad_const = RAD_CONST_E * particles.gamma^3
        elseif particles.mass == m_p
            rad_const = RAD_CONST_P * particles.gamma^3
        else
            rad_const = 0.0
            println("SR is not implemented for this particle mass.")
        end
        StrMPoleSymplectic4RadPass_SC!(r_in, ele.len, rad_const, beti, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit,
                ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, E0, num_particles, lost_flags,
                ele.a, ele.b, ele.Nl, ele.Nm, K, ele.Nsteps)
    end
    return nothing
end
###################
# multi-threading
function StrMPoleSymplectic4Pass_P_SC!(r::Matrix{Float64}, le::Float64, beti::Float64, A::Array{Float64,1}, B::Array{Float64,1}, 
    max_order::Int, num_int_step::Int, 
    FringeQuadEntrance::Int, FringeQuadExit::Int, #(no fringe), 1 (lee-whiting) or 2 (lee-whiting+elegant-like) 
    T1::Array{Float64,1}, T2::Array{Float64,1}, R1::Array{Float64,2}, R2::Array{Float64,2}, 
    RApertures::Array{Float64,1}, EApertures::Array{Float64,1}, KickAngle::Array{Float64,1}, 
    num_particles::Int, lost_flags::Array{Int64,1}, a::Float64, b::Float64, Nl::Int, Nm::Int, K::Float64, Nsteps::Int)
    # Modified based on AT function. Ref[Terebilo, Andrei. "Accelerator modeling with MATLAB accelerator toolbox." PACS2001 (2001)].
    # and [Qiang, Ji. "Differentiable self-consistent space-charge simulation for accelerator design." Physical Review Accelerators and Beams 26.2 (2023): 024601.]

    DRIFT1  =  0.6756035959798286638
    DRIFT2 = -0.1756035959798286639
    KICK1 = 1.351207191959657328
    KICK2 = -1.702414383919314656

    lstep = le / Nsteps
    SL = lstep / 2.0 /num_int_step
    L1 = SL*DRIFT1
    L2 = SL*DRIFT2
    K1 = SL*KICK1
    K2 = SL*KICK2

    if le > 0
        B[1] -= sin(KickAngle[1])/le
        A[1] += sin(KickAngle[2])/le 
    end

    for step in 1:Nsteps
        Threads.@threads for c in 1:num_particles
            if lost_flags[c] == 1
                continue
            end
            r6 = @view r[c, :]
            
            if step == 1
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

            if check_lost(r6)
                lost_flags[c] = 1
            end
        end

        space_charge_P!(r, K, Nl, Nm, a/Nl, b/Nm, a, b, num_particles, lstep, lost_flags)

        Threads.@threads for c in 1:num_particles
            if lost_flags[c] == 1
                continue
            end
            r6 = @view r[c, :]
    
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
            
            if step == Nsteps
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

function StrMPoleSymplectic4RadPass_P_SC!(r::Matrix{Float64}, len::Float64, rad_const::Float64, beti::Float64, A::Array{Float64,1}, B::Array{Float64,1}, 
    max_order::Int, num_int_step::Int, 
    FringeQuadEntrance::Int, FringeQuadExit::Int, #(no fringe), 1 (lee-whiting) or 2 (lee-whiting+elegant-like) 
    T1::Array{Float64,1}, T2::Array{Float64,1}, R1::Array{Float64,2}, R2::Array{Float64,2}, 
    RApertures::Array{Float64,1}, EApertures::Array{Float64,1}, KickAngle::Array{Float64,1}, E0::Float64,
    num_particles::Int, lost_flags::Array{Int64,1}, a::Float64, b::Float64, Nl::Int, Nm::Int, K::Float64, Nsteps::Int)
    # Ref[Terebilo, Andrei. "Accelerator modeling with MATLAB accelerator toolbox." PACS2001 (2001)].
    # and [Qiang, Ji. "Differentiable self-consistent space-charge simulation for accelerator design." Physical Review Accelerators and Beams 26.2 (2023): 024601.]

    DRIFT1  =  0.6756035959798286638
    DRIFT2 = -0.1756035959798286639
    KICK1 = 1.351207191959657328
    KICK2 = -1.702414383919314656

    lstep = le / Nsteps
    SL = lstep / 2.0 /num_int_step
    L1 = SL*DRIFT1
    L2 = SL*DRIFT2
    K1 = SL*KICK1
    K2 = SL*KICK2
    
    if le > 0
        B[1] -= sin(KickAngle[1])/ le 
        A[1] += sin(KickAngle[2])/ le
    end

    for step in 1:Nsteps
        Threads.@threads for c in 1:num_particles
            if lost_flags[c] == 1
                continue
            end
            r6 = @view r[c, :]
            if step == 1
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
            if check_lost(r6)
                lost_flags[c] = 1
            end
        end

        space_charge_P!(r, K, Nl, Nm, a/Nl, b/Nm, a, b, num_particles, le, lost_flags)

        Threads.@threads for c in 1:num_particles
            if lost_flags[c] == 1
                continue
            end
            r6 = @view r[c, :]
    
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
            
            if step == Nsteps
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

function pass_P!(ele::KQUAD_SC, r_in::Matrix{Float64}, num_particles::Int64, particles::Beam{Float64})
    # ele: KQUAD
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    lost_flags = particles.lost_flag
    K = calculate_K(particles, particles.current)
    PolynomB = zeros(4)
    E0 = particles.energy
    if use_exact_beti == 1
        beti = 1.0 / particles.beta
    else
        beti = 1.0 
    end

    PolynomB[1] = ele.k0
    PolynomB[2] = ele.k1 
    PolynomB[3] = ele.k2 / 2.0
    PolynomB[4] = ele.k3 / 6.0
    if ele.rad == 0
        StrMPoleSymplectic4Pass_P_SC!(r_in, ele.len, beti, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit, 
                ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, num_particles, lost_flags,
                ele.a, ele.b, ele.Nl, ele.Nm, K, ele.Nsteps)
    else
        rad_const = 0.0
        if particles.mass == m_e
            rad_const = RAD_CONST_E * particles.gamma^3
        elseif particles.mass == m_p
            rad_const = RAD_CONST_P * particles.gamma^3
        else
            rad_const = 0.0
            println("SR is not implemented for this particle mass.")
        end
        StrMPoleSymplectic4RadPass_P_SC!(r_in, ele.len, rad_const, beti, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit,
                ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, E0, num_particles, lost_flags,
                ele.a, ele.b, ele.Nl, ele.Nm, K, ele.Nsteps)
    end
    return nothing
end

function pass_P!(ele::KSEXT_SC, r_in::Matrix{Float64}, num_particles::Int64, particles::Beam{Float64})
    # ele: KSEXT_SC
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    lost_flags = particles.lost_flag
    PolynomB = zeros(4)
    K = calculate_K(particles, particles.current)
    E0 = particles.energy
    if use_exact_beti == 1
        beti = 1.0 / particles.beta
    else
        beti = 1.0 
    end

    PolynomB[1] = ele.k0
    PolynomB[2] = ele.k1
    PolynomB[3] = ele.k2 / 2.0
    PolynomB[4] = ele.k3 / 6.0
    if ele.rad == 0
            StrMPoleSymplectic4Pass_P_SC!(r_in, ele.len, beti, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit, 
                ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, num_particles, lost_flags,
                ele.a, ele.b, ele.Nl, ele.Nm, K, ele.Nsteps)
    else
        StrMPoleSymplectic4RadPass_P_SC!(r_in, ele.len, rad_const, beti, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps,
            ele.FringeQuadEntrance, ele.FringeQuadExit,
            ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, E0, num_particles, lost_flags,
            ele.a, ele.b, ele.Nl, ele.Nm, K, ele.Nsteps)
    end
    return nothing
end

function pass_P!(ele::KOCT_SC, r_in::Matrix{Float64}, num_particles::Int64, particles::Beam{Float64})
    # ele: KOCT_SC
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    lost_flags = particles.lost_flag
    PolynomB = zeros(4)
    K = calculate_K(particles, particles.current)
    E0 = particles.energy
    if use_exact_beti == 1
        beti = 1.0 / particles.beta
    else
        beti = 1.0 
    end

    PolynomB[1] = ele.k0
    PolynomB[2] = ele.k1
    PolynomB[3] = ele.k2 / 2.0
    PolynomB[4] = ele.k3 / 6.0
    if ele.rad == 0
        StrMPoleSymplectic4Pass_P_SC!(r_in, ele.len, beti, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
            ele.FringeQuadEntrance, ele.FringeQuadExit, 
            ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, num_particles, lost_flags,
            ele.a, ele.b, ele.Nl, ele.Nm, K, ele.Nsteps)
    else
        rad_const = 0.0
        if particles.mass == m_e
            rad_const = RAD_CONST_E * particles.gamma^3
        elseif particles.mass == m_p
            rad_const = RAD_CONST_P * particles.gamma^3
        else
            rad_const = 0.0
            println("SR is not implemented for this particle mass.")
        end
        StrMPoleSymplectic4RadPass_P_SC!(r_in, ele.len, rad_const, beti, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps,
            ele.FringeQuadEntrance, ele.FringeQuadExit,
            ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, E0, num_particles, lost_flags,
            ele.a, ele.b, ele.Nl, ele.Nm, K, ele.Nsteps)
    end
    return nothing
end

#############################################
# High-order TPSA
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
    gamma = (E0 + m0) / m0
    beta = sqrt(1.0 - 1.0 / (gamma^2))
    if use_exact_beti == 1
        beti = 1.0 / beta
    else
        beti = 1.0 
    end
    rad_const = 0.0
    PolynomB = zeros(4)

    PolynomB[1] = ele.k0
    PolynomB[2] = ele.k1
    PolynomB[3] = ele.k2 / 2.0
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
    return nothing
end

function pass_TPSA!(ele::KSEXT, r_in::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}; E0::Float64=0.0, m0::Float64=m_e) where {T, TPS_Dim, Max_TPS_Degree}
    # ele: KSEXT
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    gamma = (E0 + m0) / m0
    beta = sqrt(1.0 - 1.0 / (gamma^2))
    if use_exact_beti == 1
        beti = 1.0 /beta
    else
        beti = 1.0 
    end
    PolynomB = zeros(4)
    rad_const = 0.0

    PolynomB[1] = ele.k0
    PolynomB[2] = ele.k1 
    PolynomB[3] = ele.k2 / 2.0
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
    return nothing 
end

function pass_TPSA!(ele::KOCT, r_in::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}; E0::Float64=0.0, m0::Float64=m_e) where {T, TPS_Dim, Max_TPS_Degree}
    # ele: KOCT
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    gamma = (E0 + m0) / m0
    beta = sqrt(1.0 - 1.0 / (gamma^2))
    if use_exact_beti == 1
        beti = 1.0 / beta
    else
        beti = 1.0 
    end
    rad_const = 0.0
    PolynomB = zeros(4)

    PolynomB[1] = ele.k0
    PolynomB[2] = ele.k1 
    PolynomB[3] = ele.k2 / 2.0
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
    return nothing
end
