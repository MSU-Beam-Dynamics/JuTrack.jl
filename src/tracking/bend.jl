include("drift.jl")
include("fringe.jl")

function B2perp(bx, by, irho, x, xpr, y, ypr)
    v_norm2 = 1.0 / ((1.0 + x*irho)^2 + xpr^2 + ypr^2)
    return ((by * (1.0 + x*irho))^2 + (bx *(1.0 + x*irho))^2 + (bx*ypr - by*xpr)^2) * v_norm2
end
function bndthinkick!(r::AbstractVector{Float64}, A, B, L, irho, max_order)
# Calculate multipole kick in a curved elemrnt (bending magnet)
# The reference coordinate system  has the curvature given by the inverse 
# (design) radius irho.
# IMPORTANT !!!
# The magnetic field Bo that provides this curvature MUST NOT be included in the dipole term
# PolynomB[1](MATLAB notation)(C: B[0] in this function) of the By field expansion

# The kick is given by

#            e L      L delta      L x
# theta  = - --- B  + -------  -  -----  , 
#      x     p    y     rho           2
#             0                    rho

#          e L
# theta  = --- B
#      y    p   x
#            0

    ReSum = B[max_order + 1]
    ImSum = A[max_order + 1]
    ReSumTemp = 0.0

    for i in reverse(1:max_order)
        ReSumTemp = ReSum * r[1] - ImSum * r[3] + B[i]
        ImSum = ImSum * r[1] + ReSum * r[3] + A[i]
        ReSum = ReSumTemp
    end

    r[2] -= L * (ReSum - (r[6] - r[1] * irho) * irho)
    r[4] += L * ImSum
    r[5] += L * irho * r[1]  # Path length
    return nothing
end

function bndthinkickrad!(r::AbstractVector{Float64}, A, B, L, irho, E0, max_order)
        ReSum = B[max_order + 1]
        ImSum = A[max_order + 1]
        ReSumTemp = 0.0
        CRAD = CGAMMA * E0^3 / (2.0*pi*1e27) # [m]/[GeV^3] M.Sands (4.1)

        for i in reverse(1:max_order)
            ReSumTemp = ReSum * r[1] - ImSum * r[3] + B[i]
            ImSum = ImSum * r[1] + ReSum * r[3] + A[i]
            ReSum = ReSumTemp
        end

        # angles from momentums
        if use_exact_Hamiltonian == 1
            p_norm = 1.0 / sqrt((1.0 + r[6])^2 - r[2]^2 - r[4]^2)
        else
            p_norm = 1.0 / (1.0 + r[6])
        end
        x = r[1]
        xpr = r[2] * p_norm
        y = r[3]
        ypr = r[4] * p_norm
        B2P = B2perp(ImSum, ReSum + irho, irho, x, xpr, y, ypr)

        dp_0 = r[6]
        r[6] = r[6] - CRAD * (1.0+r[6])^2 * B2P * (1.0 + x*irho + (xpr^2 + ypr^2) / 2.0) * L
    
        # momentums after losing energy
        if use_exact_Hamiltonian == 1
            p_norm = 1.0 / sqrt((1.0 + r[6])^2 - r[2]^2 - r[4]^2)
        else
            p_norm = 1.0 / (1.0 + r[6])
        end
        r[2] = xpr / p_norm
        r[4] = ypr / p_norm

        r[2] -= L * (ReSum - (dp_0 - r[1]*irho)*irho)
        r[4] += L * ImSum
        r[5] += L * irho * r[1]
        return nothing
    end

function BendSymplecticPassRad!(r::Array{Float64,1}, le::Float64, irho::Float64, A::Array{Float64,1}, B::Array{Float64,1}, 
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

    B[1] -= sin(KickAngle[1]) / le
    A[1] += sin(KickAngle[2]) / le


    # Threads.@threads for c in 1:num_particles
    for c in 1:num_particles
        if lost_flags[c] == 1
            continue
        end
        r6 = @view r[(c-1)*6+1:c*6]
        if !isnan(r6[1])
            # Misalignment at entrance
            if T1 != zeros(6)
                ATaddvv!(r6, T1)
            end
            if R1 != zeros(6, 6)
                ATmultmv!(r6, R1)
            end

            # Aperture checks at the entrance
            # if RApertures !== nothing
            #     checkiflostRectangularAp!(r6, RApertures)
            # end
            # if EApertures !== nothing
            #     checkiflostEllipticalAp!(r6, EApertures)
            # end

            # Edge focus at entrance
            edge_fringe_entrance!(r6, irho, entrance_angle, fint1, gap, FringeBendEntrance)

            # Quadrupole gradient fringe entrance
            if FringeQuadEntrance != 0 && B[2] != 0
                if useLinFrEleEntrance == 1
                    linearQuadFringeElegantEntrance!(r6, B[2], fringeIntM0, fringeIntP0)
                else
                    QuadFringePassP!(r6, B[2])
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

            # Quadrupole gradient fringe exit
            if FringeQuadExit != 0 && B[2] != 0
                if useLinFrEleExit == 1
                    linearQuadFringeElegantExit!(r6, B[2], fringeIntM0, fringeIntP0)
                else
                    QuadFringePassN!(r6, B[2])
                end
            end

            # Edge focus at exit
            edge_fringe_exit!(r6, irho, exit_angle, fint2, gap, FringeBendExit)

            # Aperture checks at the exit
            # if RApertures !== nothing
            #     checkiflostRectangularAp!(r6, RApertures)
            # end
            # if EApertures !== nothing
            #     checkiflostEllipticalAp!(r6, EApertures)
            # end

            # Misalignment at exit
            if R2 != zeros(6, 6)
                ATmultmv!(r6, R2)
            end
            if T2 != zeros(6) 
                ATaddvv!(r6, T2)
            end
            if r6[1] > CoordLimit || r6[2] > AngleLimit || r6[1] < -CoordLimit || r6[2] < -AngleLimit
                lost_flags[c] = 1
            end
        end
    end

    
    B[1] += sin(KickAngle[1]) / le
    A[1] -= sin(KickAngle[2]) / le
    return nothing
end

function BendSymplecticPass!(r::Array{Float64,1}, le::Float64, irho::Float64, A::Array{Float64,1}, B::Array{Float64,1}, 
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

    B[1] -= sin(KickAngle[1]) / le
    A[1] += sin(KickAngle[2]) / le


    # Threads.@threads for c in 1:num_particles
    for c in 1:num_particles
        if lost_flags[c] == 1
            continue
        end
        r6 = @view r[(c-1)*6+1:c*6]
        if !isnan(r6[1])
            if use_exact_Hamiltonian == 1
                NormL1 = L1 / sqrt((1.0 + r6[6])^2 - r6[2]^2 - r6[4]^2)
                NormL2 = L2 / sqrt((1.0 + r6[6])^2 - r6[2]^2 - r6[4]^2)
            else
                NormL1 = L1 / (1.0 + r6[6])
                NormL2 = L2 / (1.0 + r6[6])
            end

            # Misalignment at entrance
            if T1 != zeros(6)
                ATaddvv!(r6, T1)
            end
            if R1 != zeros(6, 6)
                ATmultmv!(r6, R1)
            end

            # Aperture checks at the entrance
            # if RApertures !== nothing
            #     checkiflostRectangularAp!(r6, RApertures)
            # end
            # if EApertures !== nothing
            #     checkiflostEllipticalAp!(r6, EApertures)
            # end

            # Edge focus at entrance
            edge_fringe_entrance!(r6, irho, entrance_angle, fint1, gap, FringeBendEntrance)

            # Quadrupole gradient fringe entrance
            if FringeQuadEntrance != 0 && B[2] != 0
                if useLinFrEleEntrance == 1
                    linearQuadFringeElegantEntrance!(r6, B[2], fringeIntM0, fringeIntP0)
                else
                    QuadFringePassP!(r6, B[2])
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

            # Quadrupole gradient fringe exit
            if FringeQuadExit != 0 && B[2] != 0
                if useLinFrEleExit == 1
                    linearQuadFringeElegantExit!(r6, B[2], fringeIntM0, fringeIntP0)
                else
                    QuadFringePassN!(r6, B[2])
                end
            end

            # Edge focus at exit
            edge_fringe_exit!(r6, irho, exit_angle, fint2, gap, FringeBendExit)

            # Aperture checks at the exit
            # if RApertures !== nothing
            #     checkiflostRectangularAp!(r6, RApertures)
            # end
            # if EApertures !== nothing
            #     checkiflostEllipticalAp!(r6, EApertures)
            # end

            # Misalignment at exit
            if R2 != zeros(6, 6)
                ATmultmv!(r6, R2)
            end
            if T2 != zeros(6) 
                ATaddvv!(r6, T2)
            end
            if r6[1] > CoordLimit || r6[2] > AngleLimit || r6[1] < -CoordLimit || r6[2] < -AngleLimit
                lost_flags[c] = 1
            end
        end
    end

    
    B[1] += sin(KickAngle[1]) / le
    A[1] -= sin(KickAngle[2]) / le
    return nothing
end


function pass!(ele::SBEND, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam)
    # ele: SBEND
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    lost_flags = particles.lost_flag
    irho = ele.angle / ele.len
    E0 = particles.energy
    if ele.rad == 0
        BendSymplecticPass!(r_in, ele.len, irho, ele.PolynomA, ele.PolynomB, ele.MaxOrder, ele.NumIntSteps,
            ele.e1, ele.e2,
            ele.FringeBendEntrance, ele.FringeBendExit,
            ele.fint1, ele.fint2, ele.gap,
            ele.FringeQuadEntrance, ele.FringeQuadExit,
            ele.FringeIntM0, ele.FringeIntP0,
            ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures,
            ele.KickAngle, num_particles, lost_flags)
    else
        BendSymplecticPassRad!(r_in, ele.len, irho, ele.PolynomA, ele.PolynomB, ele.MaxOrder, ele.NumIntSteps,
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
function pass!(ele::RBEND, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam)
    # ele: RBEND
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    lost_flags = particles.lost_flag
    irho = ele.angle / ele.len
    E0 = particles.energy
    if ele.rad == 0
        BendSymplecticPass!(r_in, ele.len, irho, ele.PolynomA, ele.PolynomB, ele.MaxOrder, ele.NumIntSteps,
            ele.e1, ele.e2,
            ele.FringeBendEntrance, ele.FringeBendExit,
            ele.fint1, ele.fint2, ele.gap,
            ele.FringeQuadEntrance, ele.FringeQuadExit,
            ele.FringeIntM0, ele.FringeIntP0,
            ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures,
            ele.KickAngle, num_particles, lost_flags)
    else
        BendSymplecticPassRad!(r_in, ele.len, irho, ele.PolynomA, ele.PolynomB, ele.MaxOrder, ele.NumIntSteps,
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
function BendSymplecticPassRad_P!(r::Array{Float64,1}, le::Float64, irho::Float64, A::Array{Float64,1}, B::Array{Float64,1}, 
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

    B[1] -= sin(KickAngle[1]) / le
    A[1] += sin(KickAngle[2]) / le


    # Threads.@threads for c in 1:num_particles
    for c in 1:num_particles
        if lost_flags[c] == 1
            continue
        end
        r6 = @view r[(c-1)*6+1:c*6]
        if !isnan(r6[1])
            # Misalignment at entrance
            if T1 != zeros(6)
                ATaddvv!(r6, T1)
            end
            if R1 != zeros(6, 6)
                ATmultmv!(r6, R1)
            end

            # Aperture checks at the entrance
            # if RApertures !== nothing
            #     checkiflostRectangularAp!(r6, RApertures)
            # end
            # if EApertures !== nothing
            #     checkiflostEllipticalAp!(r6, EApertures)
            # end

            # Edge focus at entrance
            edge_fringe_entrance!(r6, irho, entrance_angle, fint1, gap, FringeBendEntrance)

            # Quadrupole gradient fringe entrance
            if FringeQuadEntrance != 0 && B[2] != 0
                if useLinFrEleEntrance == 1
                    linearQuadFringeElegantEntrance!(r6, B[2], fringeIntM0, fringeIntP0)
                else
                    QuadFringePassP!(r6, B[2])
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

            # Quadrupole gradient fringe exit
            if FringeQuadExit != 0 && B[2] != 0
                if useLinFrEleExit == 1
                    linearQuadFringeElegantExit!(r6, B[2], fringeIntM0, fringeIntP0)
                else
                    QuadFringePassN!(r6, B[2])
                end
            end

            # Edge focus at exit
            edge_fringe_exit!(r6, irho, exit_angle, fint2, gap, FringeBendExit)

            # Aperture checks at the exit
            # if RApertures !== nothing
            #     checkiflostRectangularAp!(r6, RApertures)
            # end
            # if EApertures !== nothing
            #     checkiflostEllipticalAp!(r6, EApertures)
            # end

            # Misalignment at exit
            if R2 != zeros(6, 6)
                ATmultmv!(r6, R2)
            end
            if T2 != zeros(6) 
                ATaddvv!(r6, T2)
            end
            if r6[1] > CoordLimit || r6[2] > AngleLimit || r6[1] < -CoordLimit || r6[2] < -AngleLimit
                lost_flags[c] = 1
            end
        end
    end

    
    B[1] += sin(KickAngle[1]) / le
    A[1] -= sin(KickAngle[2]) / le
    return nothing
end

function BendSymplecticPass_P!(r::Array{Float64,1}, le::Float64, irho::Float64, A::Array{Float64,1}, B::Array{Float64,1}, 
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

    B[1] -= sin(KickAngle[1]) / le
    A[1] += sin(KickAngle[2]) / le


    # Threads.@threads for c in 1:num_particles
    for c in 1:num_particles
        if lost_flags[c] == 1
            continue
        end
        r6 = @view r[(c-1)*6+1:c*6]
        if !isnan(r6[1])
            if use_exact_Hamiltonian == 1
                NormL1 = L1 / sqrt((1.0 + r6[6])^2 - r6[2]^2 - r6[4]^2)
                NormL2 = L2 / sqrt((1.0 + r6[6])^2 - r6[2]^2 - r6[4]^2)
            else
                NormL1 = L1 / (1.0 + r6[6])
                NormL2 = L2 / (1.0 + r6[6])
            end

            # Misalignment at entrance
            if T1 != zeros(6)
                ATaddvv!(r6, T1)
            end
            if R1 != zeros(6, 6)
                ATmultmv!(r6, R1)
            end

            # Aperture checks at the entrance
            # if RApertures !== nothing
            #     checkiflostRectangularAp!(r6, RApertures)
            # end
            # if EApertures !== nothing
            #     checkiflostEllipticalAp!(r6, EApertures)
            # end

            # Edge focus at entrance
            edge_fringe_entrance!(r6, irho, entrance_angle, fint1, gap, FringeBendEntrance)

            # Quadrupole gradient fringe entrance
            if FringeQuadEntrance != 0 && B[2] != 0
                if useLinFrEleEntrance == 1
                    linearQuadFringeElegantEntrance!(r6, B[2], fringeIntM0, fringeIntP0)
                else
                    QuadFringePassP!(r6, B[2])
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

            # Quadrupole gradient fringe exit
            if FringeQuadExit != 0 && B[2] != 0
                if useLinFrEleExit == 1
                    linearQuadFringeElegantExit!(r6, B[2], fringeIntM0, fringeIntP0)
                else
                    QuadFringePassN!(r6, B[2])
                end
            end

            # Edge focus at exit
            edge_fringe_exit!(r6, irho, exit_angle, fint2, gap, FringeBendExit)

            # Aperture checks at the exit
            # if RApertures !== nothing
            #     checkiflostRectangularAp!(r6, RApertures)
            # end
            # if EApertures !== nothing
            #     checkiflostEllipticalAp!(r6, EApertures)
            # end

            # Misalignment at exit
            if R2 != zeros(6, 6)
                ATmultmv!(r6, R2)
            end
            if T2 != zeros(6) 
                ATaddvv!(r6, T2)
            end
            if r6[1] > CoordLimit || r6[2] > AngleLimit || r6[1] < -CoordLimit || r6[2] < -AngleLimit
                lost_flags[c] = 1
            end
        end
    end

    
    B[1] += sin(KickAngle[1]) / le
    A[1] -= sin(KickAngle[2]) / le
    return nothing
end


function pass_P!(ele::SBEND, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam)
    # ele: SBEND
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    lost_flags = particles.lost_flag
    irho = ele.angle / ele.len
    E0 = particles.energy
    if ele.rad == 0
        BendSymplecticPass_P!(r_in, ele.len, irho, ele.PolynomA, ele.PolynomB, ele.MaxOrder, ele.NumIntSteps,
            ele.e1, ele.e2,
            ele.FringeBendEntrance, ele.FringeBendExit,
            ele.fint1, ele.fint2, ele.gap,
            ele.FringeQuadEntrance, ele.FringeQuadExit,
            ele.FringeIntM0, ele.FringeIntP0,
            ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures,
            ele.KickAngle, num_particles, lost_flags)
    else
        BendSymplecticPassRad_P!(r_in, ele.len, irho, ele.PolynomA, ele.PolynomB, ele.MaxOrder, ele.NumIntSteps,
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
function pass_P!(ele::RBEND, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam)
    # ele: RBEND
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    lost_flags = particles.lost_flag
    irho = ele.angle / ele.len
    E0 = particles.energy
    if ele.rad == 0
        BendSymplecticPass_P!(r_in, ele.len, irho, ele.PolynomA, ele.PolynomB, ele.MaxOrder, ele.NumIntSteps,
            ele.e1, ele.e2,
            ele.FringeBendEntrance, ele.FringeBendExit,
            ele.fint1, ele.fint2, ele.gap,
            ele.FringeQuadEntrance, ele.FringeQuadExit,
            ele.FringeIntM0, ele.FringeIntP0,
            ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures,
            ele.KickAngle, num_particles, lost_flags)
    else
        BendSymplecticPassRad_P!(r_in, ele.len, irho, ele.PolynomA, ele.PolynomB, ele.MaxOrder, ele.NumIntSteps,
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
