function BendSymplecticPassRad!(r::Array{Float64,1}, le::Float64, irho::Float64, A::Array{Float64,1}, B::Array{Float64,1}, 
    max_order::Int, num_int_steps::Int, entrance_angle::Float64, exit_angle::Float64, FringeBendEntrance::Int, FringeBendExit::Int,
    fint1::Float64, fint2::Float64, gap::Float64, FringeQuadEntrance::Int, FringeQuadExit::Int,
    fringeIntM0::Array{Float64,1}, fringeIntP0::Array{Float64,1}, T1::Array{Float64,1}, T2::Array{Float64,1}, 
    R1::Array{Float64,2}, R2::Array{Float64,2}, RApertures::Array{Float64,1}, EApertures::Array{Float64,1},
    KickAngle::Array{Float64,1}, E0::Float64, num_particles::Int, lost_flags::Array{Int64,1},
    a::Float64, b::Float64, Nl::Int, Nm::Int, K::Float64, Nsteps::Int)
    # Ref[Terebilo, Andrei. "Accelerator modeling with MATLAB accelerator toolbox." PACS2001 (2001)].
    # and [Qiang, Ji. "Differentiable self-consistent space-charge simulation for accelerator design." Physical Review Accelerators and Beams 26.2 (2023): 024601.]

    DRIFT1 = 0.6756035959798286638
    DRIFT2 = -0.1756035959798286639
    KICK1 = 1.351207191959657328
    KICK2 = -1.702414383919314656
    lstep = le / Nsteps
    SL = lstep / 2.0 / num_int_steps
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

    for step in 1:Nsteps
        for c in 1:num_particles
            if isone(lost_flags[c])
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
                if step == 1
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
                end

                # Integrator
                for m in 1:num_int_steps
                    drift6!(r6, L1)
                    bndthinkickrad!(r6, A, B, K1, irho, E0, max_order)
                    drift6!(r6, L2)
                    bndthinkickrad!(r6, A, B, K2, irho, E0, max_order)
                    drift6!(r6, L2)
                    bndthinkickrad!(r6, A, B, K1, irho, E0, max_order)
                    drift6!(r6, L1)
                end

                if check_lost(r6)
                    lost_flags[c] = 1
                end
            end
        end

        space_charge!(r, K, Nl, Nm, a/Nl, b/Nm, a, b, num_particles, lstep, lost_flags)

        for c in 1:num_particles
            if isone(lost_flags[c])
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
                if step == 1
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
                end

                # Integrator
                for m in 1:num_int_steps
                    drift6!(r6, L1)
                    bndthinkickrad!(r6, A, B, K1, irho, E0, max_order)
                    drift6!(r6, L2)
                    bndthinkickrad!(r6, A, B, K2, irho, E0, max_order)
                    drift6!(r6, L2)
                    bndthinkickrad!(r6, A, B, K1, irho, E0, max_order)
                    drift6!(r6, L1)
                end
                
                if step == Nsteps
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
                end

                if check_lost(r6)
                    lost_flags[c] = 1
                end
            end
        end
    end

    
    B[1] += sin(KickAngle[1]) / le
    A[1] -= sin(KickAngle[2]) / le
    return nothing
end

function BendSymplecticPass_SC!(r::Array{Float64,1}, le::Float64, irho::Float64, A::Array{Float64,1}, B::Array{Float64,1}, 
    max_order::Int, num_int_steps::Int, entrance_angle::Float64, exit_angle::Float64, FringeBendEntrance::Int, FringeBendExit::Int,
    fint1::Float64, fint2::Float64, gap::Float64, FringeQuadEntrance::Int, FringeQuadExit::Int,
    fringeIntM0::Array{Float64,1}, fringeIntP0::Array{Float64,1}, T1::Array{Float64,1}, T2::Array{Float64,1}, 
    R1::Array{Float64,2}, R2::Array{Float64,2}, RApertures::Array{Float64,1}, EApertures::Array{Float64,1},
    KickAngle::Array{Float64,1}, num_particles::Int, lost_flags::Array{Int64,1}, 
    a::Float64, b::Float64, Nl::Int, Nm::Int, K::Float64, Nsteps::Int)
    # Modified based on AT function. Ref[Terebilo, Andrei. "Accelerator modeling with MATLAB accelerator toolbox." PACS2001 (2001)].
    # and [Qiang, Ji. "Differentiable self-consistent space-charge simulation for accelerator design." Physical Review Accelerators and Beams 26.2 (2023): 024601.]

    DRIFT1 = 0.6756035959798286638
    DRIFT2 = -0.1756035959798286639
    KICK1 = 1.351207191959657328
    KICK2 = -1.702414383919314656

    lstep = le / Nsteps
    SL = lstep / 2.0 / num_int_steps
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


    for step in 1:Nsteps
        for c in 1:num_particles
            if isone(lost_flags[c])
                continue
            end
            r6 = @view r[(c-1)*6+1:c*6]
            if !isnan(r6[1])
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
                end
    
                # Integrator
                for m in 1:num_int_steps
                    fastdrift!(r6, NormL1, L1)
                    bndthinkick!(r6, A, B, K1, irho, max_order)
                    fastdrift!(r6, NormL2, L2)
                    bndthinkick!(r6, A, B, K2, irho, max_order)
                    fastdrift!(r6, NormL2, L2)
                    bndthinkick!(r6, A, B, K1, irho, max_order)
                    fastdrift!(r6, NormL1, L1)
                end

                if check_lost(r6)
                    lost_flags[c] = 1
                end
            end
        end

        space_charge!(r, K, Nl, Nm, a/Nl, b/Nm, a, b, num_particles, lstep, lost_flags)

        for c in 1:num_particles
            if isone(lost_flags[c])
                continue
            end
            r6 = @view r[(c-1)*6+1:c*6]
            NormL1 = L1 / (1.0 + r6[6])
            NormL2 = L2 / (1.0 + r6[6])
            
            # Integrator
            for m in 1:num_int_steps
                fastdrift!(r6, NormL1, L1)
                bndthinkick!(r6, A, B, K1, irho, max_order)
                fastdrift!(r6, NormL2, L2)
                bndthinkick!(r6, A, B, K2, irho, max_order)
                fastdrift!(r6, NormL2, L2)
                bndthinkick!(r6, A, B, K1, irho, max_order)
                fastdrift!(r6, NormL1, L1)
            end

            if step == Nsteps
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

function pass!(ele::SBEND_SC, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam)
    # ele: SBEND_SC
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    # particles: Beam object
    lost_flags = particles.lost_flag
    irho = ele.angle / ele.len
    K = calculate_K(particles, particles.current)
    E0 = particles.energy
    if ele.rad == 0
        BendSymplecticPass_SC!(r_in, ele.len, irho, ele.PolynomA, ele.PolynomB, ele.MaxOrder, ele.NumIntSteps,
            ele.e1, ele.e2,
            ele.FringeBendEntrance, ele.FringeBendExit,
            ele.fint1, ele.fint2, ele.gap,
            ele.FringeQuadEntrance, ele.FringeQuadExit,
            ele.FringeIntM0, ele.FringeIntP0,
            ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures,
            ele.KickAngle, num_particles, lost_flags, ele.a, ele.b, ele.Nl, ele.Nm, K, ele.Nsteps)
    else
        BendSymplecticPassRad_SC!(r_in, ele.len, irho, ele.PolynomA, ele.PolynomB, ele.MaxOrder, ele.NumIntSteps,
            ele.e1, ele.e2,
            ele.FringeBendEntrance, ele.FringeBendExit,
            ele.fint1, ele.fint2, ele.gap,
            ele.FringeQuadEntrance, ele.FringeQuadExit,
            ele.FringeIntM0, ele.FringeIntP0,
            ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures,
            ele.KickAngle, E0, num_particles, lost_flags, ele.a, ele.b, ele.Nl, ele.Nm, K, ele.Nsteps)
    end
    return nothing
end

######################
# multi-threading
function BendSymplecticPassRad_P_SC!(r::Array{Float64,1}, le::Float64, irho::Float64, A::Array{Float64,1}, B::Array{Float64,1}, 
    max_order::Int, num_int_steps::Int, entrance_angle::Float64, exit_angle::Float64, FringeBendEntrance::Int, FringeBendExit::Int,
    fint1::Float64, fint2::Float64, gap::Float64, FringeQuadEntrance::Int, FringeQuadExit::Int,
    fringeIntM0::Array{Float64,1}, fringeIntP0::Array{Float64,1}, T1::Array{Float64,1}, T2::Array{Float64,1}, 
    R1::Array{Float64,2}, R2::Array{Float64,2}, RApertures::Array{Float64,1}, EApertures::Array{Float64,1},
    KickAngle::Array{Float64,1}, E0::Float64, num_particles::Int, lost_flags::Array{Int64,1},
    a::Float64, b::Float64, Nl::Int, Nm::Int, K::Float64, Nsteps::Int)
    # Ref[Terebilo, Andrei. "Accelerator modeling with MATLAB accelerator toolbox." PACS2001 (2001)].
    # and [Qiang, Ji. "Differentiable self-consistent space-charge simulation for accelerator design." Physical Review Accelerators and Beams 26.2 (2023): 024601.]

    DRIFT1 = 0.6756035959798286638
    DRIFT2 = -0.1756035959798286639
    KICK1 = 1.351207191959657328
    KICK2 = -1.702414383919314656
    lstep = le / Nsteps
    SL = lstep / 2.0 / num_int_steps
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

    for step in 1:Nsteps
        Threads.@threads for c in 1:num_particles
            if isone(lost_flags[c])
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
                if step == 1
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
                end

                # Integrator
                for m in 1:num_int_steps
                    drift6!(r6, L1)
                    bndthinkickrad!(r6, A, B, K1, irho, E0, max_order)
                    drift6!(r6, L2)
                    bndthinkickrad!(r6, A, B, K2, irho, E0, max_order)
                    drift6!(r6, L2)
                    bndthinkickrad!(r6, A, B, K1, irho, E0, max_order)
                    drift6!(r6, L1)
                end

                if check_lost(r6)
                    lost_flags[c] = 1
                end
            end
        end

        space_charge_P!(r, K, Nl, Nm, a/Nl, b/Nm, a, b, num_particles, lstep, lost_flags)

        Threads.@threads for c in 1:num_particles
            if isone(lost_flags[c])
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
                if step == 1
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
                end

                # Integrator
                for m in 1:num_int_steps
                    drift6!(r6, L1)
                    bndthinkickrad!(r6, A, B, K1, irho, E0, max_order)
                    drift6!(r6, L2)
                    bndthinkickrad!(r6, A, B, K2, irho, E0, max_order)
                    drift6!(r6, L2)
                    bndthinkickrad!(r6, A, B, K1, irho, E0, max_order)
                    drift6!(r6, L1)
                end
                
                if step == Nsteps
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
                end

                if check_lost(r6)
                    lost_flags[c] = 1
                end
            end
        end
    end

    
    B[1] += sin(KickAngle[1]) / le
    A[1] -= sin(KickAngle[2]) / le
    return nothing
end

function BendSymplecticPass_P_SC!(r::Array{Float64,1}, le::Float64, irho::Float64, A::Array{Float64,1}, B::Array{Float64,1}, 
    max_order::Int, num_int_steps::Int, entrance_angle::Float64, exit_angle::Float64, FringeBendEntrance::Int, FringeBendExit::Int,
    fint1::Float64, fint2::Float64, gap::Float64, FringeQuadEntrance::Int, FringeQuadExit::Int,
    fringeIntM0::Array{Float64,1}, fringeIntP0::Array{Float64,1}, T1::Array{Float64,1}, T2::Array{Float64,1}, 
    R1::Array{Float64,2}, R2::Array{Float64,2}, RApertures::Array{Float64,1}, EApertures::Array{Float64,1},
    KickAngle::Array{Float64,1}, num_particles::Int, lost_flags::Array{Int64,1}, 
    a::Float64, b::Float64, Nl::Int, Nm::Int, K::Float64, Nsteps::Int)
    # Modified based on AT function. Ref[Terebilo, Andrei. "Accelerator modeling with MATLAB accelerator toolbox." PACS2001 (2001)].
    # and [Qiang, Ji. "Differentiable self-consistent space-charge simulation for accelerator design." Physical Review Accelerators and Beams 26.2 (2023): 024601.]

    DRIFT1 = 0.6756035959798286638
    DRIFT2 = -0.1756035959798286639
    KICK1 = 1.351207191959657328
    KICK2 = -1.702414383919314656

    lstep = le / Nsteps
    SL = lstep / 2.0 / num_int_steps
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


    for step in 1:Nsteps
        Threads.@threads for c in 1:num_particles
            if isone(lost_flags[c])
                continue
            end
            r6 = @view r[(c-1)*6+1:c*6]
            if !isnan(r6[1])
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
                end
    
                # Integrator
                for m in 1:num_int_steps
                    fastdrift!(r6, NormL1, L1)
                    bndthinkick!(r6, A, B, K1, irho, max_order)
                    fastdrift!(r6, NormL2, L2)
                    bndthinkick!(r6, A, B, K2, irho, max_order)
                    fastdrift!(r6, NormL2, L2)
                    bndthinkick!(r6, A, B, K1, irho, max_order)
                    fastdrift!(r6, NormL1, L1)
                end

                if check_lost(r6)
                    lost_flags[c] = 1
                end
            end
        end

        space_charge_P!(r, K, Nl, Nm, a/Nl, b/Nm, a, b, num_particles, lstep, lost_flags)

        Threads.@threads for c in 1:num_particles
            if isone(lost_flags[c])
                continue
            end
            r6 = @view r[(c-1)*6+1:c*6]
            NormL1 = L1 / (1.0 + r6[6])
            NormL2 = L2 / (1.0 + r6[6])
            
            # Integrator
            for m in 1:num_int_steps
                fastdrift!(r6, NormL1, L1)
                bndthinkick!(r6, A, B, K1, irho, max_order)
                fastdrift!(r6, NormL2, L2)
                bndthinkick!(r6, A, B, K2, irho, max_order)
                fastdrift!(r6, NormL2, L2)
                bndthinkick!(r6, A, B, K1, irho, max_order)
                fastdrift!(r6, NormL1, L1)
            end

            if step == Nsteps
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

function pass_P!(ele::SBEND_SC, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam)
    # ele: SBEND_SC
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    lost_flags = particles.lost_flag
    irho = ele.angle / ele.len
    E0 = particles.energy
    K = calculate_K(particles, particles.current)
    if ele.rad == 0
        BendSymplecticPass_P_SC!(r_in, ele.len, irho, ele.PolynomA, ele.PolynomB, ele.MaxOrder, ele.NumIntSteps,
            ele.e1, ele.e2,
            ele.FringeBendEntrance, ele.FringeBendExit,
            ele.fint1, ele.fint2, ele.gap,
            ele.FringeQuadEntrance, ele.FringeQuadExit,
            ele.FringeIntM0, ele.FringeIntP0,
            ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures,
            ele.KickAngle, num_particles, lost_flags, ele.a, ele.b, ele.Nl, ele.Nm, K, ele.Nsteps)
    else
        BendSymplecticPassRad_P_SC!(r_in, ele.len, irho, ele.PolynomA, ele.PolynomB, ele.MaxOrder, ele.NumIntSteps,
            ele.e1, ele.e2,
            ele.FringeBendEntrance, ele.FringeBendExit,
            ele.fint1, ele.fint2, ele.gap,
            ele.FringeQuadEntrance, ele.FringeQuadExit,
            ele.FringeIntM0, ele.FringeIntP0,
            ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures,
            ele.KickAngle, E0, num_particles, lost_flags, ele.a, ele.b, ele.Nl, ele.Nm, K, ele.Nsteps)
    end
    return nothing
end