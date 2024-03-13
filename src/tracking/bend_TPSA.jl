# include("drift_TPSA.jl")
include("drift_TPSA.jl")
include("fringe_TPSA.jl")

function bndthinkick!(r::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}, A, B, L, irho, max_order) where {T, TPS_Dim, Max_TPS_Degree}
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
        # for c in 1:num_particles
        #     r6 = @view r[(c-1)*6+1:c*6]
            # if !isnan(r6[1])
            if use_exact_Hamiltonian == 1
                NormL1 = L1 / sqrt((1.0 + r[6])^2 - r[2]^2 - r[4]^2)
                NormL2 = L2 / sqrt((1.0 + r[6])^2 - r[2]^2 - r[4]^2)
            else
                NormL1 = L1 / (1.0 + r[6])
                NormL2 = L2 / (1.0 + r[6])
            end
                # Misalignment at entrance
                if T1 != zeros(6)
                    ATaddvv!(r, T1)
                end
                if R1 != zeros(6, 6)
                    ATmultmv!(r, R1)
                end
    
                # Edge focus at entrance
                edge_fringe_entrance!(r, irho, entrance_angle, fint1, gap, FringeBendEntrance)
    
                # Quadrupole gradient fringe entrance
                # if FringeQuadEntrance != 0 && B[2] != 0
                #     if useLinFrEleEntrance == 1
                #         linearQuadFringeElegantEntrance!(r6, B[2], fringeIntM0, fringeIntP0)
                #     else
                #         QuadFringePassP!(r6, B[2])
                #     end
                # end
    
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
    
                # Quadrupole gradient fringe exit
                # if FringeQuadExit != 0 && B[2] != 0
                #     if useLinFrEleExit == 1
                #         linearQuadFringeElegantExit!(r6, B[2], fringeIntM0, fringeIntP0)
                #     else
                #         QuadFringePassN!(r6, B[2])
                #     end
                # end
    
                # Edge focus at exit
                edge_fringe_exit!(r, irho, exit_angle, fint2, gap, FringeBendExit)
    
    
                # Misalignment at exit
                if R2 != zeros(6, 6)
                    ATmultmv!(r, R2)
                end
                if T2 != zeros(6)
                    ATaddvv!(r, T2)
                end
            # end
        # end
    
        
        B[1] += sin(KickAngle[1]) / le
        A[1] -= sin(KickAngle[2]) / le
        return nothing
    end
    
    
function pass_TPSA!(ele::SBEND, r_in::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}) where {T, TPS_Dim, Max_TPS_Degree}
    # ele: SBEND
    # r_in: 6-by-1 TPSA array
    if ele.rad != 0
        println("Synchrotron radiation is not included in the SBEND element for TPSA tracking.")
    end
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
    function pass_TPSA!(ele::RBEND, r_in::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}) where {T, TPS_Dim, Max_TPS_Degree}
        # ele: RBEND
        # r_in: 6-by-1 TPSA array
        if ele.rad != 0
            println("Synchrotron radiation is not included in the SBEND element for TPSA tracking.")
        end
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
# using BenchmarkTools
# using Enzyme
# include("../lattice/canonical_elements_AT.jl")
# # include("../TPSA_Enzyme/arrayTPSA_fixedmap.jl")
# include("../TPSA_Enzyme/TPSA_fixedmap.jl")
# # # # q = KQUAD(PolynomialB=[0.0, 1.0, 0.0, 0.0])
# function f(xx)
#     q = SBEND(angle=xx[1], len=1.0)
#     x = CTPS(0.0, 1, 6, 3)
#     xp = CTPS(0.0, 2, 6, 3)
#     y = CTPS(0.0, 3, 6, 3)
#     yp = CTPS(0.0, 4, 6, 3)
#     z = CTPS(0.0, 5, 6, 3)
#     delta = CTPS(0.0, 6, 6, 3)
#     rin = [x, xp, y, yp, z, delta]
#     pass_TPSA!(q, rin, 1)
# return rin[1].map[2]
# end
# x = [0.1573]
# println(f(x))
# @btime grad = gradient(Forward, f, x)
# println(grad)