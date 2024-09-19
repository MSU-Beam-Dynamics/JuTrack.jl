function StrMPoleSymplectic4Pass_SC!(r::Array{Float64,1}, le::Float64, A::Array{Float64,1}, B::Array{Float64,1}, 
    max_order::Int, num_int_step::Int, 
    FringeQuadEntrance::Int, FringeQuadExit::Int, #(no fringe), 1 (lee-whiting) or 2 (lee-whiting+elegant-like) 
    fringeIntM0::Array{Float64,1},  # I0m/K1, I1m/K1, I2m/K1, I3m/K1, Lambda2m/K1 
    fringeIntP0::Array{Float64,1},  # I0p/K1, I1p/K1, I2p/K1, I3p/K1, Lambda2p/K1
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

    if le > 0
        B[1] -= sin(KickAngle[1])/le
        A[1] += sin(KickAngle[2])/le 
    end

    for step in 1:Nsteps
        for c in 1:num_particles
            if lost_flags[c] == 1
                continue
            end
            r6 = @view r[(c-1)*6+1:c*6]
            NormL1 = L1 / (1.0 + r6[6])
            NormL2 = L2 / (1.0 + r6[6])
            
            if step == 1
                # Misalignment at entrance
                if !iszero(T1)
                    addvv!(r6, T1)
                end
                if !iszero(R1)
                    multmv!(r6, R1)
                end

                if FringeQuadEntrance != 0 && B[2] != 0
                    if useLinFrEleEntrance == 1
                        linearQuadFringeElegantEntrance!(r6, B[2], fringeIntM0, fringeIntP0)
                    else
                        QuadFringePassP!(r6, B[2])
                    end
                end
            end

            # Integrator
            for m in 1:num_int_step
                fastdrift!(r6, NormL1, L1)
                strthinkick!(r6, A, B, K1, max_order)
                fastdrift!(r6, NormL2, L2)
                strthinkick!(r6, A, B, K2, max_order)
                fastdrift!(r6, NormL2, L2)
                strthinkick!(r6, A, B, K1, max_order)
                fastdrift!(r6, NormL1, L1)
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
            r6 = @view r[(c-1)*6+1:c*6]

            NormL1 = L1 / (1.0 + r6[6])
            NormL2 = L2 / (1.0 + r6[6])
    
            # Integrator
            for m in 1:num_int_step
                fastdrift!(r6, NormL1, L1)
                strthinkick!(r6, A, B, K1, max_order)
                fastdrift!(r6, NormL2, L2)
                strthinkick!(r6, A, B, K2, max_order)
                fastdrift!(r6, NormL2, L2)
                strthinkick!(r6, A, B, K1, max_order)
                fastdrift!(r6, NormL1, L1)
            end
            
            if step == Nsteps
                if FringeQuadExit != 0 && B[2] != 0
                    if useLinFrEleExit == 1
                        linearQuadFringeElegantExit!(r6, B[2], fringeIntM0, fringeIntP0)
                    else
                        QuadFringePassN!(r6, B[2])
                    end
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

function StrMPoleSymplectic4RadPass_SC!(r::Array{Float64,1}, le::Float64, A::Array{Float64,1}, B::Array{Float64,1}, 
    max_order::Int, num_int_step::Int, 
    FringeQuadEntrance::Int, FringeQuadExit::Int, #(no fringe), 1 (lee-whiting) or 2 (lee-whiting+elegant-like) 
    fringeIntM0::Array{Float64,1},  # I0m/K1, I1m/K1, I2m/K1, I3m/K1, Lambda2m/K1 
    fringeIntP0::Array{Float64,1},  # I0p/K1, I1p/K1, I2p/K1, I3p/K1, Lambda2p/K1
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

    if FringeQuadEntrance==2# && !isnothing(fringeIntM0) && !isnothing(fringeIntP0)
        useLinFrEleEntrance = 1
    else
        useLinFrEleEntrance = 0
    end
    if FringeQuadExit==2#&& !isnothing(fringeIntM0) && !isnothing(fringeIntP0)
        useLinFrEleExit = 1
    else
        useLinFrEleExit = 0
    end

    if le > 0
        B[1] -= sin(KickAngle[1])/ le 
        A[1] += sin(KickAngle[2])/ le
    end

    for step in 1:Nsteps
        for c in 1:num_particles
            if lost_flags[c] == 1
                continue
            end
            r6 = @view r[(c-1)*6+1:c*6]
            if step == 1
                # Misalignment at entrance
                if !iszero(T1)
                    addvv!(r6, T1)
                end
                if !iszero(R1)
                    multmv!(r6, R1)
                end
    
                if FringeQuadEntrance != 0 && B[2] != 0
                    if useLinFrEleEntrance == 1
                        linearQuadFringeElegantEntrance!(r6, B[2], fringeIntM0, fringeIntP0)
                    else
                        QuadFringePassP!(r6, B[2])
                    end
                end
            end
            # Integrator
            for m in 1:num_int_step
                drift6!(r6, L1)
                strthinkickrad!(r6, A, B, K1, E0, max_order)
                drift6!(r6, L2)
                strthinkickrad!(r6, A, B, K2, E0, max_order)
                drift6!(r6, L2)
                strthinkickrad!(r6, A, B, K1, E0, max_order)
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
            r6 = @view r[(c-1)*6+1:c*6]
    
            # Integrator
            for m in 1:num_int_step
                drift6!(r6, L1)
                strthinkickrad!(r6, A, B, K1, E0, max_order)
                drift6!(r6, L2)
                strthinkickrad!(r6, A, B, K2, E0, max_order)
                drift6!(r6, L2)
                strthinkickrad!(r6, A, B, K1, E0, max_order)
                drift6!(r6, L1)
            end
            
            if step == Nsteps
                if FringeQuadExit != 0 && B[2] != 0
                    if useLinFrEleExit == 1
                        linearQuadFringeElegantExit!(r6, B[2], fringeIntM0, fringeIntP0)
                    else
                        QuadFringePassN!(r6, B[2])
                    end
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

function pass!(ele::KQUAD_SC, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam)
    # ele: KQUAD_SC
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    lost_flags = particles.lost_flag
    K = calculate_K(particles, particles.current)
    PolynomB = zeros(4)
    E0 = particles.energy
    if ele.PolynomB[1] == 0.0 && ele.PolynomB[2] == 0.0 && ele.PolynomB[3] == 0.0 && ele.PolynomB[4] == 0.0
        PolynomB[2] = ele.k1
        if ele.rad == 0
            StrMPoleSymplectic4Pass_SC!(r_in,  ele.len, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                    ele.FringeQuadEntrance, ele.FringeQuadExit, ele.FringeIntM0, ele.FringeIntP0, 
                    ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, num_particles, lost_flags,
                    ele.a, ele.b, ele.Nl, ele.Nm, K, ele.Nsteps)
        else
            StrMPoleSymplectic4RadPass_SC!(r_in,  ele.len, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                    ele.FringeQuadEntrance, ele.FringeQuadExit, ele.FringeIntM0, ele.FringeIntP0, 
                    ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, E0, num_particles, lost_flags,
                    ele.a, ele.b, ele.Nl, ele.Nm, K, ele.Nsteps)
        end
    else
        PolynomB[1] = ele.PolynomB[1]
        PolynomB[2] = ele.PolynomB[2] 
        PolynomB[3] = ele.PolynomB[3] / 2.0
        PolynomB[4] = ele.PolynomB[4] / 6.0
        if ele.rad == 0
            StrMPoleSymplectic4Pass_SC!(r_in,  ele.len, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                    ele.FringeQuadEntrance, ele.FringeQuadExit, ele.FringeIntM0, ele.FringeIntP0, 
                    ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, num_particles, lost_flags,
                    ele.a, ele.b, ele.Nl, ele.Nm, K, ele.Nsteps)
        else
            StrMPoleSymplectic4RadPass_SC!(r_in,  ele.len, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                    ele.FringeQuadEntrance, ele.FringeQuadExit, ele.FringeIntM0, ele.FringeIntP0, 
                    ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, E0, num_particles, lost_flags,
                    ele.a, ele.b, ele.Nl, ele.Nm, K, ele.Nsteps)
        end
    end
    return nothing
end

function pass!(ele::KSEXT_SC, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam)
    # ele: KSEXT
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    lost_flags = particles.lost_flag
    K = calculate_K(particles, particles.current)
    PolynomB = zeros(4)
    E0 = particles.energy
    if ele.PolynomB[1] == 0.0 && ele.PolynomB[2] == 0.0 && ele.PolynomB[3] == 0.0 && ele.PolynomB[4] == 0.0
        PolynomB[3] = ele.k2 / 2.0
        if ele.rad == 0
            StrMPoleSymplectic4Pass_SC!(r_in, ele.len, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                    ele.FringeQuadEntrance, ele.FringeQuadExit, ele.FringeIntM0, ele.FringeIntP0, 
                    ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, num_particles, lost_flags,
                    ele.a, ele.b, ele.Nl, ele.Nm, K, ele.Nsteps)
        else
            StrMPoleSymplectic4RadPass_SC!(r_in, ele.len, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                    ele.FringeQuadEntrance, ele.FringeQuadExit, ele.FringeIntM0, ele.FringeIntP0, 
                    ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, E0, num_particles, lost_flags,
                    ele.a, ele.b, ele.Nl, ele.Nm, K, ele.Nsteps)
        end
    else
        PolynomB[1] = ele.PolynomB[1]
        PolynomB[2] = ele.PolynomB[2] 
        PolynomB[3] = ele.PolynomB[3] / 2.0
        PolynomB[4] = ele.PolynomB[4] / 6.0
        if ele.rad == 0
            StrMPoleSymplectic4Pass_SC!(r_in, ele.len, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                    ele.FringeQuadEntrance, ele.FringeQuadExit, ele.FringeIntM0, ele.FringeIntP0, 
                    ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, num_particles, lost_flags,
                    ele.a, ele.b, ele.Nl, ele.Nm, K, ele.Nsteps)
        else
            StrMPoleSymplectic4RadPass_SC!(r_in, ele.len, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                    ele.FringeQuadEntrance, ele.FringeQuadExit, ele.FringeIntM0, ele.FringeIntP0, 
                    ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, E0, num_particles, lost_flags,
                    ele.a, ele.b, ele.Nl, ele.Nm, K, ele.Nsteps)
        end
    end
    return nothing
end

function pass!(ele::KOCT_SC, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam)
    # ele: KOCT_SC
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    lost_flags = particles.lost_flag
    PolynomB = zeros(4)
    K = calculate_K(particles, particles.current)
    E0 = particles.energy
    if ele.PolynomB[1] == 0.0 && ele.PolynomB[2] == 0.0 && ele.PolynomB[3] == 0.0 && ele.PolynomB[4] == 0.0
        PolynomB[4] = ele.k3 / 6.0
        if ele.rad == 0
            StrMPoleSymplectic4Pass_SC!(r_in, ele.len, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                    ele.FringeQuadEntrance, ele.FringeQuadExit, ele.FringeIntM0, ele.FringeIntP0, 
                    ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, num_particles, lost_flags,
                    ele.a, ele.b, ele.Nl, ele.Nm, K, ele.Nsteps)
        else
            StrMPoleSymplectic4RadPass_SC!(r_in, ele.len, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                    ele.FringeQuadEntrance, ele.FringeQuadExit, ele.FringeIntM0, ele.FringeIntP0, 
                    ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, E0, num_particles, lost_flags,
                    ele.a, ele.b, ele.Nl, ele.Nm, K, ele.Nsteps)
        end
    else
        PolynomB[1] = ele.PolynomB[1]
        PolynomB[2] = ele.PolynomB[2] 
        PolynomB[3] = ele.PolynomB[3] / 2.0
        PolynomB[4] = ele.PolynomB[4] / 6.0
        if ele.rad == 0
            StrMPoleSymplectic4Pass_SC!(r_in, ele.len, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                    ele.FringeQuadEntrance, ele.FringeQuadExit, ele.FringeIntM0, ele.FringeIntP0, 
                    ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, num_particles, lost_flags,
                    ele.a, ele.b, ele.Nl, ele.Nm, K, ele.Nsteps)
        else
            StrMPoleSymplectic4RadPass_SC!(r_in, ele.len, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                    ele.FringeQuadEntrance, ele.FringeQuadExit, ele.FringeIntM0, ele.FringeIntP0, 
                    ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, E0, num_particles, lost_flags,
                    ele.a, ele.b, ele.Nl, ele.Nm, K, ele.Nsteps)
        end
    end
    return nothing
end


###################
# multi-threading
function StrMPoleSymplectic4Pass_P_SC!(r::Array{Float64,1}, le::Float64, A::Array{Float64,1}, B::Array{Float64,1}, 
    max_order::Int, num_int_step::Int, 
    FringeQuadEntrance::Int, FringeQuadExit::Int, #(no fringe), 1 (lee-whiting) or 2 (lee-whiting+elegant-like) 
    fringeIntM0::Array{Float64,1},  # I0m/K1, I1m/K1, I2m/K1, I3m/K1, Lambda2m/K1 
    fringeIntP0::Array{Float64,1},  # I0p/K1, I1p/K1, I2p/K1, I3p/K1, Lambda2p/K1
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

    if le > 0
        B[1] -= sin(KickAngle[1])/le
        A[1] += sin(KickAngle[2])/le 
    end

    for step in 1:Nsteps
        Threads.@threads for c in 1:num_particles
            if lost_flags[c] == 1
                continue
            end
            r6 = @view r[(c-1)*6+1:c*6]
            NormL1 = L1 / (1.0 + r6[6])
            NormL2 = L2 / (1.0 + r6[6])
            
            if step == 1
                # Misalignment at entrance
                if !iszero(T1)
                    addvv!(r6, T1)
                end
                if !iszero(R1)
                    multmv!(r6, R1)
                end

                if FringeQuadEntrance != 0 && B[2] != 0
                    if useLinFrEleEntrance == 1
                        linearQuadFringeElegantEntrance!(r6, B[2], fringeIntM0, fringeIntP0)
                    else
                        QuadFringePassP!(r6, B[2])
                    end
                end
            end

            # Integrator
            for m in 1:num_int_step
                fastdrift!(r6, NormL1, L1)
                strthinkick!(r6, A, B, K1, max_order)
                fastdrift!(r6, NormL2, L2)
                strthinkick!(r6, A, B, K2, max_order)
                fastdrift!(r6, NormL2, L2)
                strthinkick!(r6, A, B, K1, max_order)
                fastdrift!(r6, NormL1, L1)
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
            r6 = @view r[(c-1)*6+1:c*6]

            NormL1 = L1 / (1.0 + r6[6])
            NormL2 = L2 / (1.0 + r6[6])
    
            # Integrator
            for m in 1:num_int_step
                fastdrift!(r6, NormL1, L1)
                strthinkick!(r6, A, B, K1, max_order)
                fastdrift!(r6, NormL2, L2)
                strthinkick!(r6, A, B, K2, max_order)
                fastdrift!(r6, NormL2, L2)
                strthinkick!(r6, A, B, K1, max_order)
                fastdrift!(r6, NormL1, L1)
            end
            
            if step == Nsteps
                if FringeQuadExit != 0 && B[2] != 0
                    if useLinFrEleExit == 1
                        linearQuadFringeElegantExit!(r6, B[2], fringeIntM0, fringeIntP0)
                    else
                        QuadFringePassN!(r6, B[2])
                    end
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

function StrMPoleSymplectic4RadPass_P_SC!(r::Array{Float64,1}, le::Float64, A::Array{Float64,1}, B::Array{Float64,1}, 
    max_order::Int, num_int_step::Int, 
    FringeQuadEntrance::Int, FringeQuadExit::Int, #(no fringe), 1 (lee-whiting) or 2 (lee-whiting+elegant-like) 
    fringeIntM0::Array{Float64,1},  # I0m/K1, I1m/K1, I2m/K1, I3m/K1, Lambda2m/K1 
    fringeIntP0::Array{Float64,1},  # I0p/K1, I1p/K1, I2p/K1, I3p/K1, Lambda2p/K1
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

    if FringeQuadEntrance==2# && !isnothing(fringeIntM0) && !isnothing(fringeIntP0)
        useLinFrEleEntrance = 1
    else
        useLinFrEleEntrance = 0
    end
    if FringeQuadExit==2#&& !isnothing(fringeIntM0) && !isnothing(fringeIntP0)
        useLinFrEleExit = 1
    else
        useLinFrEleExit = 0
    end

    if le > 0
        B[1] -= sin(KickAngle[1])/ le 
        A[1] += sin(KickAngle[2])/ le
    end

    for step in 1:Nsteps
        Threads.@threads for c in 1:num_particles
            if lost_flags[c] == 1
                continue
            end
            r6 = @view r[(c-1)*6+1:c*6]
            if step == 1
                # Misalignment at entrance
                if !iszero(T1)
                    addvv!(r6, T1)
                end
                if !iszero(R1)
                    multmv!(r6, R1)
                end
    
                if FringeQuadEntrance != 0 && B[2] != 0
                    if useLinFrEleEntrance == 1
                        linearQuadFringeElegantEntrance!(r6, B[2], fringeIntM0, fringeIntP0)
                    else
                        QuadFringePassP!(r6, B[2])
                    end
                end
            end
            # Integrator
            for m in 1:num_int_step
                drift6!(r6, L1)
                strthinkickrad!(r6, A, B, K1, E0, max_order)
                drift6!(r6, L2)
                strthinkickrad!(r6, A, B, K2, E0, max_order)
                drift6!(r6, L2)
                strthinkickrad!(r6, A, B, K1, E0, max_order)
                drift6!(r6, L1)
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
            r6 = @view r[(c-1)*6+1:c*6]
    
            # Integrator
            for m in 1:num_int_step
                drift6!(r6, L1)
                strthinkickrad!(r6, A, B, K1, E0, max_order)
                drift6!(r6, L2)
                strthinkickrad!(r6, A, B, K2, E0, max_order)
                drift6!(r6, L2)
                strthinkickrad!(r6, A, B, K1, E0, max_order)
                drift6!(r6, L1)
            end
            
            if step == Nsteps
                if FringeQuadExit != 0 && B[2] != 0
                    if useLinFrEleExit == 1
                        linearQuadFringeElegantExit!(r6, B[2], fringeIntM0, fringeIntP0)
                    else
                        QuadFringePassN!(r6, B[2])
                    end
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

function pass_P!(ele::KQUAD_SC, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam)
    # ele: KQUAD
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    lost_flags = particles.lost_flag
    K = calculate_K(particles, particles.current)
    PolynomB = zeros(4)
    E0 = particles.energy
    if ele.PolynomB[1] == 0.0 && ele.PolynomB[2] == 0.0 && ele.PolynomB[3] == 0.0 && ele.PolynomB[4] == 0.0
        PolynomB[2] = ele.k1
        if ele.rad == 0
            StrMPoleSymplectic4Pass_P_SC!(r_in, ele.len, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                    ele.FringeQuadEntrance, ele.FringeQuadExit, ele.FringeIntM0, ele.FringeIntP0, 
                    ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, num_particles, lost_flags,
                    ele.a, ele.b, ele.Nl, ele.Nm, K, ele.Nsteps)
        else
            StrMPoleSymplectic4RadPass_P_SC!(r_in, ele.len, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                    ele.FringeQuadEntrance, ele.FringeQuadExit, ele.FringeIntM0, ele.FringeIntP0, 
                    ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, E0, num_particles, lost_flags,
                    ele.a, ele.b, ele.Nl, ele.Nm, K, ele.Nsteps)
        end
    else
        PolynomB[1] = ele.PolynomB[1]
        PolynomB[2] = ele.PolynomB[2] 
        PolynomB[3] = ele.PolynomB[3] / 2.0
        PolynomB[4] = ele.PolynomB[4] / 6.0
        if ele.rad == 0
            StrMPoleSymplectic4Pass_P_SC!(r_in, ele.len, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                    ele.FringeQuadEntrance, ele.FringeQuadExit, ele.FringeIntM0, ele.FringeIntP0, 
                    ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, num_particles, lost_flags,
                    ele.a, ele.b, ele.Nl, ele.Nm, K, ele.Nsteps)
        else
            StrMPoleSymplectic4RadPass_P_SC!(r_in, ele.len, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                    ele.FringeQuadEntrance, ele.FringeQuadExit, ele.FringeIntM0, ele.FringeIntP0, 
                    ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, E0, num_particles, lost_flags,
                    ele.a, ele.b, ele.Nl, ele.Nm, K, ele.Nsteps)
        end
    end
    return nothing
end

function pass_P!(ele::KSEXT_SC, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam)
    # ele: KSEXT_SC
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    lost_flags = particles.lost_flag
    PolynomB = zeros(4)
    K = calculate_K(particles, particles.current)
    E0 = particles.energy
    if ele.PolynomB[1] == 0.0 && ele.PolynomB[2] == 0.0 && ele.PolynomB[3] == 0.0 && ele.PolynomB[4] == 0.0
        PolynomB[3] = ele.k2 / 2.0
        if ele.rad == 0
                StrMPoleSymplectic4Pass_P_SC!(r_in, ele.len, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                    ele.FringeQuadEntrance, ele.FringeQuadExit, ele.FringeIntM0, ele.FringeIntP0, 
                    ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, num_particles, lost_flags,
                    ele.a, ele.b, ele.Nl, ele.Nm, K, ele.Nsteps)
        else
                StrMPoleSymplectic4RadPass_P_SC!(r_in, ele.len, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                    ele.FringeQuadEntrance, ele.FringeQuadExit, ele.FringeIntM0, ele.FringeIntP0, 
                    ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, E0, num_particles, lost_flags,
                    ele.a, ele.b, ele.Nl, ele.Nm, K, ele.Nsteps)
        end
    else
        PolynomB[1] = ele.PolynomB[1]
        PolynomB[2] = ele.PolynomB[2]
        PolynomB[3] = ele.PolynomB[3] / 2.0
        PolynomB[4] = ele.PolynomB[4] / 6.0
        if ele.rad == 0
                StrMPoleSymplectic4Pass_P_SC!(r_in, ele.len, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                    ele.FringeQuadEntrance, ele.FringeQuadExit, ele.FringeIntM0, ele.FringeIntP0, 
                    ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, num_particles, lost_flags,
                    ele.a, ele.b, ele.Nl, ele.Nm, K, ele.Nsteps)
        else
                StrMPoleSymplectic4RadPass_P_SC!(r_in, ele.len, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                    ele.FringeQuadEntrance, ele.FringeQuadExit, ele.FringeIntM0, ele.FringeIntP0, 
                    ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, E0, num_particles, lost_flags,
                    ele.a, ele.b, ele.Nl, ele.Nm, K, ele.Nsteps)
        end
    end
    return nothing
end

function pass_P!(ele::KOCT_SC, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam)
    # ele: KOCT_SC
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    lost_flags = particles.lost_flag
    PolynomB = zeros(4)
    K = calculate_K(particles, particles.current)
    E0 = particles.energy
    if ele.PolynomB[1] == 0.0 && ele.PolynomB[2] == 0.0 && ele.PolynomB[3] == 0.0 && ele.PolynomB[4] == 0.0
        PolynomB[4] = ele.k3 / 6.0
        if ele.rad == 0
                StrMPoleSymplectic4Pass_P_SC!(r_in, ele.len, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                    ele.FringeQuadEntrance, ele.FringeQuadExit, ele.FringeIntM0, ele.FringeIntP0, 
                    ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, num_particles, lost_flags,
                    ele.a, ele.b, ele.Nl, ele.Nm, K, ele.Nsteps)
        else
                StrMPoleSymplectic4RadPass_P_SC!(r_in, ele.len, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                    ele.FringeQuadEntrance, ele.FringeQuadExit, ele.FringeIntM0, ele.FringeIntP0, 
                    ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, E0, num_particles, lost_flags,
                    ele.a, ele.b, ele.Nl, ele.Nm, K, ele.Nsteps)

        end
    else
        PolynomB[1] = ele.PolynomB[1]
        PolynomB[2] = ele.PolynomB[2] 
        PolynomB[3] = ele.PolynomB[3] / 2.0
        PolynomB[4] = ele.PolynomB[4] / 6.0
        if ele.rad == 0
                StrMPoleSymplectic4Pass_P_SC!(r_in, ele.len, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                    ele.FringeQuadEntrance, ele.FringeQuadExit, ele.FringeIntM0, ele.FringeIntP0, 
                    ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, num_particles, lost_flags,
                    ele.a, ele.b, ele.Nl, ele.Nm, K, ele.Nsteps)
        else
                StrMPoleSymplectic4RadPass_P_SC!(r_in, ele.len, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                    ele.FringeQuadEntrance, ele.FringeQuadExit, ele.FringeIntM0, ele.FringeIntP0, 
                    ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, E0, num_particles, lost_flags,
                    ele.a, ele.b, ele.Nl, ele.Nm, K, ele.Nsteps)
        end
    end
    return nothing
end