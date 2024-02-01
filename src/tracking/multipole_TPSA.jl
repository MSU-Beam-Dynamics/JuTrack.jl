include("drift_TPSA.jl")
# include("fringe_AT.jl")

function strthinkick!(r::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}, A, B, L, max_order) where {T, TPS_Dim, Max_TPS_Degree}
    # Calculate and apply a multipole kick to a 6-dimentional
    # phase space vector in a straight element (quadrupole)
    
    # IMPORTANT !!!
    # The reference coordinate system is straight but the field expansion may still
    # contain dipole terms A[1], B[1]

    ### this section is type unstable
    # ReSum = B[max_order + 1]
    # ImSum = A[max_order + 1]
    # ReSumTemp = zeros(length(r[1]))

    # i = 1
    ReSumTemp = B[max_order + 1] * r[1] - A[max_order + 1] * r[3] + B[1]
    ImSum = A[max_order + 1] * r[1] + B[max_order + 1] * r[3] + A[1]
    ReSum = CTPS(ReSumTemp)
    # ReSumTemp = tadd(tminus(tmult(B[max_order + 1], r[1]), tmult(A[max_order + 1], r[3])), B[1])
    # ImSum = tadd(tadd(tmult(A[max_order + 1], r[1]), tmult(B[max_order + 1], r[3])), A[1])
    # ReSum = CTPS(ReSumTemp)
    for i in reverse(2:max_order)
        # ReSumTemp = tadd(tminus(tmult(ReSum, r[1]), tmult(ImSum, r[3])), B[i])
        # ImSum = tadd(tadd(tmult(ImSum, r[1]), tmult(ReSum, r[3])), A[i])
        # ReSum = CTPS(ReSumTemp)
        ReSumTemp = ReSum * r[1] - ImSum * r[3] + B[i]
        ImSum = ImSum * r[1] + ReSum * r[3] + A[i]
        ReSum = CTPS(ReSumTemp)
    end

    r[2] -= L * ReSum
    r[4] += L * ImSum
    # r[2] = tminus(r[2], tmult(L, ReSum))
    # r[4] = tadd(r[4], tmult(L, ImSum))

    return nothing
end


function StrMPoleSymplectic4Pass!(r::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}, le, A, B, max_order, num_int_step, 
    FringeQuadEntrance, FringeQuadExit, #(no fringe), 1 (lee-whiting) or 2 (lee-whiting+elegant-like) 
    fringeIntM0,  # I0m/K1, I1m/K1, I2m/K1, I3m/K1, Lambda2m/K1 
    fringeIntP0,  # I0p/K1, I1p/K1, I2p/K1, I3p/K1, Lambda2p/K1
    T1, T2, R1, R2, RApertures, EApertures, KickAngle, num_particles) where {T, TPS_Dim, Max_TPS_Degree}

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
            norm  = 1.0 / (1.0 + r[5])
            NormL1 = L1 * norm
            NormL2 = L2 * norm
            # norm = tdiv(1.0, tadd(1.0, r[5]))
            # NormL1 = tmult(L1, norm)
            # NormL2 = tmult(L2, norm)

            # Misalignment at entrance
            if !isnothing(T1)
                ATaddvv!(r, T1)
            end
            if !isnothing(R1)
                ATmultmv!(r, R1)
            end

            # fringe effect is not implemented for TPSA
            # if FringeQuadEntrance != 0 && B[2] != 0
            #     if useLinFrEleEntrance == 1
            #         linearQuadFringeElegantEntrance!(r6, B[2], fringeIntM0, fringeIntP0)
            #     else
            #         QuadFringePassP!(r6, B[2])
            #     end
            # end

            # Integrator
            for m in 1:num_int_step
                fastdrift!(r, NormL1)
                strthinkick!(r, A, B, K1, max_order)
                fastdrift!(r, NormL2)
                strthinkick!(r, A, B, K2, max_order)
                fastdrift!(r, NormL2)
                strthinkick!(r, A, B, K1, max_order)
                fastdrift!(r, NormL1)
            end

            # if FringeQuadExit != 0 && B[2] != 0
            #     if useLinFrEleExit == 1
            #         linearQuadFringeElegantExit!(r6, B[2], fringeIntM0, fringeIntP0)
            #     else
            #         QuadFringePassN!(r6, B[2])
            #     end
            # end

            # Misalignment at exit
            if !isnothing(R2)
                ATmultmv!(r, R2)
            end
            if !isnothing(T2)
                ATaddvv!(r, T2)
            end
        # end
    end

    B[1] += sin(KickAngle[1]) / le
    A[1] -= sin(KickAngle[2]) / le
    return nothing
end

function pass_TPSA!(ele::KQUAD, r_in::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}, num_particles::Int64) where {T, TPS_Dim, Max_TPS_Degree}
    # ele: KQUAD
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    PolynomB = zeros(4)
    if ele.PolynomB[1] == 0.0 && ele.PolynomB[2] == 0.0 && ele.PolynomB[3] == 0.0 && ele.PolynomB[4] == 0.0
        PolynomB[2] = ele.k1
        StrMPoleSymplectic4Pass!(r_in, ele.len, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
            ele.FringeQuadEntrance, ele.FringeQuadExit, ele.FringeIntM0, ele.FringeIntP0, 
            ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, num_particles)
    else
        PolynomB[1] = ele.PolynomB[1]
        PolynomB[2] = ele.PolynomB[2]
        PolynomB[3] = ele.PolynomB[3]
        PolynomB[4] = ele.PolynomB[4]
        StrMPoleSymplectic4Pass!(r_in, ele.len, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
            ele.FringeQuadEntrance, ele.FringeQuadExit, ele.FringeIntM0, ele.FringeIntP0, 
            ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, num_particles)
    end
    return nothing
end

function pass_TPSA!(ele::KSEXT, r_in::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}, num_particles::Int64) where {T, TPS_Dim, Max_TPS_Degree}
    # ele: KSEXT
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    PolynomB = zeros(4)
    if ele.PolynomB[1] == 0.0 && ele.PolynomB[2] == 0.0 && ele.PolynomB[3] == 0.0 && ele.PolynomB[4] == 0.0
        PolynomB[3] = ele.k2
        StrMPoleSymplectic4Pass!(r_in, ele.len, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
            ele.FringeQuadEntrance, ele.FringeQuadExit, ele.FringeIntM0, ele.FringeIntP0, 
            ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, num_particles)
    else
        PolynomB[1] = ele.PolynomB[1]
        PolynomB[2] = ele.PolynomB[2]
        PolynomB[3] = ele.PolynomB[3]
        PolynomB[4] = ele.PolynomB[4]
        StrMPoleSymplectic4Pass!(r_in, ele.len, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
            ele.FringeQuadEntrance, ele.FringeQuadExit, ele.FringeIntM0, ele.FringeIntP0, 
            ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, num_particles)
    end
    return nothing 
end

function pass_TPSA!(ele::KOCT, r_in::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}, num_particles::Int64) where {T, TPS_Dim, Max_TPS_Degree}
    # ele: KOCT
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    PolynomB = zeros(4)
    if ele.PolynomB[1] == 0.0 && ele.PolynomB[2] == 0.0 && ele.PolynomB[3] == 0.0 && ele.PolynomB[4] == 0.0
        PolynomB[4] = ele.k3
        StrMPoleSymplectic4Pass!(r_in, ele.len, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
            ele.FringeQuadEntrance, ele.FringeQuadExit, ele.FringeIntM0, ele.FringeIntP0, 
            ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, num_particles)
    else
        PolynomB[1] = ele.PolynomB[1]
        PolynomB[2] = ele.PolynomB[2]
        PolynomB[3] = ele.PolynomB[3]
        PolynomB[4] = ele.PolynomB[4]
        StrMPoleSymplectic4Pass!(r_in, ele.len, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
            ele.FringeQuadEntrance, ele.FringeQuadExit, ele.FringeIntM0, ele.FringeIntP0, 
            ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, num_particles)
    end
    return nothing
end

# using Enzyme
# include("../lattice/canonical_elements.jl")
# include("../TPSA_Enzyme/TPSA_fixedmap.jl")
# # # # q = KQUAD(PolynomialB=[0.0, 1.0, 0.0, 0.0])
# function f(xx)
#     q = KQUAD(len=1.0, k1=xx[1])
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
# x = [-1.063770]
# println(f(x))
# # grad = gradient(Forward, f, x)
# using BenchmarkTools
# # # @btime f([1.0])
# @btime grad = gradient(Forward, f, [1.0])
# # println(grad)