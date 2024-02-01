include("drift.jl")
include("fringe.jl")

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

    r[2] -= L * (ReSum - (r[5] - r[1] * irho) * irho)
    r[4] += L * ImSum
    r[6] += L * irho * r[1]  # Path length
    return nothing
end

function BndMPoleSymplectic4Pass!(r::Array{Float64,1}, le::Float64, irho::Float64, A::Array{Float64,1}, B::Array{Float64,1}, 
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
            p_norm = 1.0 / (1.0 + r6[5])
            NormL1 = L1 * p_norm
            NormL2 = L2 * p_norm

            # Misalignment at entrance
            if isnothing(T1)
                ATaddvv!(r6, T1)
            end
            if isnothing(R1)
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
                fastdrift!(r6, NormL1)
                bndthinkick!(r6, A, B, K1, irho, max_order)
                fastdrift!(r6, NormL2)
                bndthinkick!(r6, A, B, K2, irho, max_order)
                fastdrift!(r6, NormL2)
                bndthinkick!(r6, A, B, K1, irho, max_order)
                fastdrift!(r6, NormL1)
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
            if R2 !== nothing
                ATmultmv!(r6, R2)
            end
            if T2 !== nothing
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


function pass!(ele::SBEND, r_in::Array{Float64,1}, num_particles::Int64, lost_flags::Array{Int64,1})
    # ele: SBEND
    # r_in: 6-by-num_particles array
    # num_particles: number of particles

    irho = ele.angle / ele.len
    BndMPoleSymplectic4Pass!(r_in, ele.len, irho, ele.PolynomA, ele.PolynomB, ele.MaxOrder, ele.NumIntSteps,
        ele.e1, ele.e2,
        ele.FringeBendEntrance, ele.FringeBendExit,
        ele.fint1, ele.fint2, ele.gap,
        ele.FringeQuadEntrance, ele.FringeQuadExit,
        ele.FringeIntM0, ele.FringeIntP0,
        ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures,
        ele.KickAngle, num_particles, lost_flags)
    return nothing
end
# function f(L, angle)
#     particles = [0.001 0.0001 0.0005 0.0002 0.0 0.0 0.001 0.0 0.0 0.0 0.0 0.0]
#     T1 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
#     T2 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
#     R1 = [1.0 0.0 0.0 0.0 0.0 0.0; 
#           0.0 1.0 0.0 0.0 0.0 0.0; 
#           0.0 0.0 1.0 0.0 0.0 0.0; 
#           0.0 0.0 0.0 1.0 0.0 0.0; 
#           0.0 0.0 0.0 0.0 1.0 0.0; 
#           0.0 0.0 0.0 0.0 0.0 1.0]
#     R2 = [1.0 0.0 0.0 0.0 0.0 0.0;
#             0.0 1.0 0.0 0.0 0.0 0.0;
#             0.0 0.0 1.0 0.0 0.0 0.0;
#             0.0 0.0 0.0 1.0 0.0 0.0;
#             0.0 0.0 0.0 0.0 1.0 0.0;
#             0.0 0.0 0.0 0.0 0.0 1.0]
#     irho = angle/L 
    
#     BndMPoleSymplectic4Pass!(particles, L, irho, [0.0], [0.0], 0, 10, 0.01, 0.02, 0, 0, 
#         0.0, 0.0, 0.0, 0, 0, nothing, nothing, T1, T2, R1, R2, nothing, nothing, 
#         [0.0, 0.0], 1)
#     # println(particles)
#     return particles
#     end
#     print(f(1.0,0.01))