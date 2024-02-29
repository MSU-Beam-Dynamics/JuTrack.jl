include("drift.jl")
include("fringe.jl")

function StrB2perp(bx::Float64, by::Float64, x::Float64, xpr::Float64, y::Float64, ypr::Float64)
    # Calculates sqr(|B x e|) , where e is a unit vector in the direction of velocity
    v_norm2 = 1.0 / (1.0 + xpr^2 + ypr^2)
    return (by^2 + bx^2 + (bx*ypr - by*xpr)^2) * v_norm2
end

function strthinkickrad!(r::AbstractVector{Float64}, A, B, L, E0, max_order)
    ReSum = B[max_order + 1]
    ImSum = A[max_order + 1]
    ReSumTemp = 0.0
    irho = 0.0 # straight elements
    CRAD = CGAMMA * E0^3 / (2.0*pi*1e27) # [m]/[GeV^3] M.Sands (4.1)

    for i in reverse(1:max_order)
        ReSumTemp = ReSum * r[1] - ImSum * r[3] + B[i]
        ImSum = ImSum * r[1] + ReSum * r[3] + A[i]
        ReSum = ReSumTemp
    end

    # angles for momentum
    # p_norm = 1.0 / (1.0 + r[6])
    p_norm = 1.0 / sqrt((1.0 + r[6])^2 - r[2]^2 - r[4]^2)
    x = r[1]
    xpr = r[2] * p_norm
    y = r[3]
    ypr = r[4] * p_norm
    B2P = StrB2perp(ImSum, ReSum , x , xpr, y ,ypr)

    dp_0 = r[6]
    r[6] = r[6] - CRAD * (1.0+r[6])^2 * B2P * (1.0 + x*irho + (xpr^2 + ypr^2) / 2.0) * L

    # momentums after losing energy
    p_norm = 1.0 / sqrt((1.0 + r[6])^2 - r[2]^2 - r[4]^2)
    r[2] = xpr / p_norm
    r[4] = ypr / p_norm

    r[2] -= L * (ReSum - (dp_0 - r[1]*irho)*irho)
    r[4] += L * ImSum
    r[5] += L * irho * r[1]
    return nothing
end

function strthinkick!(r::AbstractVector{Float64}, A, B, L, max_order)
    # Calculate and apply a multipole kick to a 6-dimentional
    # phase space vector in a straight element (quadrupole)
    
    # IMPORTANT !!!
    # The reference coordinate system is straight but the field expansion may still
    # contain dipole terms A[1], B[1]

    ReSum = B[max_order + 1]
    ImSum = A[max_order + 1]
    ReSumTemp = 0.0

    for i in reverse(1:max_order)
        ReSumTemp = ReSum * r[1] - ImSum * r[3] + B[i]
        ImSum = ImSum * r[1] + ReSum * r[3] + A[i]
        ReSum = ReSumTemp
    end

    r[2] -= L * ReSum
    r[4] += L * ImSum
    return nothing
end


function StrMPoleSymplectic4Pass!(r::Array{Float64,1}, le::Float64, A::Array{Float64,1}, B::Array{Float64,1}, 
    max_order::Int, num_int_step::Int, 
    FringeQuadEntrance::Int, FringeQuadExit::Int, #(no fringe), 1 (lee-whiting) or 2 (lee-whiting+elegant-like) 
    fringeIntM0::Array{Float64,1},  # I0m/K1, I1m/K1, I2m/K1, I3m/K1, Lambda2m/K1 
    fringeIntP0::Array{Float64,1},  # I0p/K1, I1p/K1, I2p/K1, I3p/K1, Lambda2p/K1
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

    B[1] -= sin(KickAngle[1])/le
    A[1] += sin(KickAngle[2])/le

    # Threads.@threads for c in 1:num_particles
    for c in 1:num_particles
        if lost_flags[c] == 1
            continue
        end
        r6 = @view r[(c-1)*6+1:c*6]
        if !isnan(r6[1])
            NormL1 = L1 / sqrt((1.0 + r6[6])^2 - r6[2]^2 - r6[4]^2)
            NormL2 = L2 / sqrt((1.0 + r6[6])^2 - r6[2]^2 - r6[4]^2)
            # NormL1 = L1 / (1.0 + r6[6])
            # NormL2 = L2 / (1.0 + r6[6])

            # Misalignment at entrance
            if T1 != zeros(6)
                ATaddvv!(r6, T1)
            end
            if R1 != zeros(6, 6)
                ATmultmv!(r6, R1)
            end

            # Check physical apertures at the entrance of the magnet
            # if RApertures != nothing
            #     checkiflostRectangularAp(r6, RApertures)
            # end
            # if EApertures != nothing
            #     checkiflostEllipticalAp(r6, EApertures)
            # end

            if FringeQuadEntrance != 0 && B[2] != 0
                if useLinFrEleEntrance == 1
                    linearQuadFringeElegantEntrance!(r6, B[2], fringeIntM0, fringeIntP0)
                else
                    QuadFringePassP!(r6, B[2])
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

            if FringeQuadExit != 0 && B[2] != 0
                if useLinFrEleExit == 1
                    linearQuadFringeElegantExit!(r6, B[2], fringeIntM0, fringeIntP0)
                else
                    QuadFringePassN!(r6, B[2])
                end
            end

            # Check physical apertures at the exit of the magnet
            # if RApertures != nothing
            #     checkiflostRectangularAp(r6, RApertures)
            # end
            # if EApertures != nothing
            #     checkiflostEllipticalAp(r6, EApertures)
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

function StrMPoleSymplectic4RadPass!(r::Array{Float64,1}, le::Float64, A::Array{Float64,1}, B::Array{Float64,1}, 
    max_order::Int, num_int_step::Int, 
    FringeQuadEntrance::Int, FringeQuadExit::Int, #(no fringe), 1 (lee-whiting) or 2 (lee-whiting+elegant-like) 
    fringeIntM0::Array{Float64,1},  # I0m/K1, I1m/K1, I2m/K1, I3m/K1, Lambda2m/K1 
    fringeIntP0::Array{Float64,1},  # I0p/K1, I1p/K1, I2p/K1, I3p/K1, Lambda2p/K1
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

    B[1] -= sin(KickAngle[1])/le
    A[1] += sin(KickAngle[2])/le

    # Threads.@threads for c in 1:num_particles
    for c in 1:num_particles
        if lost_flags[c] == 1
            continue
        end
        r6 = @view r[(c-1)*6+1:c*6]
        if !isnan(r6[1])
            # NormL1 = L1 / (1.0 + r6[6])
            # NormL2 = L2 / (1.0 + r6[6])
            NormL1 = L1 / sqrt((1.0 + r6[6])^2 - r6[2]^2 - r6[4]^2)
            NormL2 = L2 / sqrt((1.0 + r6[6])^2 - r6[2]^2 - r6[4]^2)

            # Misalignment at entrance
            if T1 != zeros(6)
                ATaddvv!(r6, T1)
            end
            if R1 != zeros(6, 6)
                ATmultmv!(r6, R1)
            end

            # Check physical apertures at the entrance of the magnet
            # if RApertures != nothing
            #     checkiflostRectangularAp(r6, RApertures)
            # end
            # if EApertures != nothing
            #     checkiflostEllipticalAp(r6, EApertures)
            # end

            if FringeQuadEntrance != 0 && B[2] != 0
                if useLinFrEleEntrance == 1
                    linearQuadFringeElegantEntrance!(r6, B[2], fringeIntM0, fringeIntP0)
                else
                    QuadFringePassP!(r6, B[2])
                end
            end

            # Integrator
            for m in 1:num_int_step
                fastdrift!(r6, NormL1, L1)
                strthinkickrad!(r6, A, B, K1, E0, max_order)
                fastdrift!(r6, NormL2, L2)
                strthinkickrad!(r6, A, B, K2, E0, max_order)
                fastdrift!(r6, NormL2, L2)
                strthinkickrad!(r6, A, B, K1, E0, max_order)
                fastdrift!(r6, NormL1, L1)
            end

            if FringeQuadExit != 0 && B[2] != 0
                if useLinFrEleExit == 1
                    linearQuadFringeElegantExit!(r6, B[2], fringeIntM0, fringeIntP0)
                else
                    QuadFringePassN!(r6, B[2])
                end
            end

            # Check physical apertures at the exit of the magnet
            # if RApertures != nothing
            #     checkiflostRectangularAp(r6, RApertures)
            # end
            # if EApertures != nothing
            #     checkiflostEllipticalAp(r6, EApertures)
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

function pass!(ele::KQUAD, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam)
    # ele: KQUAD
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    lost_flags = particles.lost_flag
    PolynomB = zeros(4)
    E0 = particles.energy
    if ele.PolynomB[1] == 0.0 && ele.PolynomB[2] == 0.0 && ele.PolynomB[3] == 0.0 && ele.PolynomB[4] == 0.0
        PolynomB[2] = ele.k1
        if ele.rad == 0
            StrMPoleSymplectic4Pass!(r_in, ele.len, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit, ele.FringeIntM0, ele.FringeIntP0, 
                ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, num_particles, lost_flags)
        else
            StrMPoleSymplectic4RadPass!(r_in, ele.len, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit, ele.FringeIntM0, ele.FringeIntP0, 
                ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, E0, num_particles, lost_flags)
        end
    else
        PolynomB[1] = ele.PolynomB[1]
        PolynomB[2] = ele.PolynomB[2]
        PolynomB[3] = ele.PolynomB[3]
        PolynomB[4] = ele.PolynomB[4]
        if ele.rad == 0
            StrMPoleSymplectic4Pass!(r_in, ele.len, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit, ele.FringeIntM0, ele.FringeIntP0, 
                ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, num_particles, lost_flags)
        else
            StrMPoleSymplectic4RadPass!(r_in, ele.len, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit, ele.FringeIntM0, ele.FringeIntP0, 
                ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, E0, num_particles, lost_flags)
        end
    end
    return nothing
end

function pass!(ele::KSEXT, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam)
    # ele: KSEXT
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    lost_flags = particles.lost_flag
    PolynomB = zeros(4)
    E0 = particles.energy
    if ele.PolynomB[1] == 0.0 && ele.PolynomB[2] == 0.0 && ele.PolynomB[3] == 0.0 && ele.PolynomB[4] == 0.0
        PolynomB[3] = ele.k2
        if ele.rad == 0
            StrMPoleSymplectic4Pass!(r_in, ele.len, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit, ele.FringeIntM0, ele.FringeIntP0, 
                ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, num_particles, lost_flags)
        else
            StrMPoleSymplectic4RadPass!(r_in, ele.len, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit, ele.FringeIntM0, ele.FringeIntP0, 
                ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, E0, num_particles, lost_flags)
        end
    else
        PolynomB[1] = ele.PolynomB[1]
        PolynomB[2] = ele.PolynomB[2]
        PolynomB[3] = ele.PolynomB[3]
        PolynomB[4] = ele.PolynomB[4]
        if ele.rad == 0
            StrMPoleSymplectic4Pass!(r_in, ele.len, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit, ele.FringeIntM0, ele.FringeIntP0, 
                ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, num_particles, lost_flags)
        else
            StrMPoleSymplectic4RadPass!(r_in, ele.len, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit, ele.FringeIntM0, ele.FringeIntP0, 
                ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, E0, num_particles, lost_flags)
        end
    end
    return nothing
end

function pass!(ele::KOCT, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam)
    # ele: KOCT
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    lost_flags = particles.lost_flag
    PolynomB = zeros(4)
    E0 = particles.energy
    if ele.PolynomB[1] == 0.0 && ele.PolynomB[2] == 0.0 && ele.PolynomB[3] == 0.0 && ele.PolynomB[4] == 0.0
        PolynomB[4] = ele.k3
        if ele.rad == 0
            StrMPoleSymplectic4Pass!(r_in, ele.len, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit, ele.FringeIntM0, ele.FringeIntP0, 
                ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, num_particles, lost_flags)
        else
            StrMPoleSymplectic4RadPass!(r_in, ele.len, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit, ele.FringeIntM0, ele.FringeIntP0, 
                ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, E0, num_particles, lost_flags)
        end
    else
        PolynomB[1] = ele.PolynomB[1]
        PolynomB[2] = ele.PolynomB[2]
        PolynomB[3] = ele.PolynomB[3]
        PolynomB[4] = ele.PolynomB[4]
        if ele.rad == 0
            StrMPoleSymplectic4Pass!(r_in, ele.len, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit, ele.FringeIntM0, ele.FringeIntP0, 
                ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, num_particles, lost_flags)
        else
            StrMPoleSymplectic4RadPass!(r_in, ele.len, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit, ele.FringeIntM0, ele.FringeIntP0, 
                ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, E0, num_particles, lost_flags)
        end
    end
    return nothing
end

# using Enzyme
# # q = KQUAD(PolynomialB=[0.0, 1.0, 0.0, 0.0])
# include("../lattice/canonical_elements_AT.jl")
# function f(k)
#     particles = [0.001, 0.0001, 0.0005, 0.0002, 0.0, 0.0, 0.001, 0.0, 0.0, 0.0, 0.0, 0.0]
#     Q = KQUAD(k1=k, len=1.0)
#     pass!(Q, particles, 2)
#     return particles
# end
# using BenchmarkTools
# @btime f(0.1)
