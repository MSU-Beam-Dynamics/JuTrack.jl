include("drift_AT.jl")
include("fringe_AT.jl")

function strthinkick!(r, A, B, L, max_order)
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
end


function StrMPoleSymplectic4Pass!(r, le, A, B, max_order, num_int_step, 
    FringeQuadEntrance, FringeQuadExit, #(no fringe), 1 (lee-whiting) or 2 (lee-whiting+elegant-like) 
    fringeIntM0,  # I0m/K1, I1m/K1, I2m/K1, I3m/K1, Lambda2m/K1 
    fringeIntP0,  # I0p/K1, I1p/K1, I2p/K1, I3p/K1, Lambda2p/K1
    T1, T2, R1, R2, RApertures, EApertures, KickAngle, num_particles)

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

    Threads.@threads for c in 1:num_particles
        r6 = @view r[(c-1)*6+1:c*6]
        if !isnan(r6[1])
            norm = 1.0 / (1.0 + r6[5])
            NormL1 = L1 * norm
            NormL2 = L2 * norm

            # Misalignment at entrance
            if !isnothing(T1)
                ATaddvv!(r6, T1)
            end
            if !isnothing(R1)
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
                fastdrift!(r6, NormL1)
                strthinkick!(r6, A, B, K1, max_order)
                fastdrift!(r6, NormL2)
                strthinkick!(r6, A, B, K2, max_order)
                fastdrift!(r6, NormL2)
                strthinkick!(r6, A, B, K1, max_order)
                fastdrift!(r6, NormL1)
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
            if !isnothing(R2)
                ATmultmv!(r6, R2)
            end
            if !isnothing(T2)
                ATaddvv!(r6, T2)
            end
        end
    end

    B[1] += sin(KickAngle[1]) / le
    A[1] -= sin(KickAngle[2]) / le
end

# using Base: @kwdef
# abstract type AbstractElement end

# @kwdef struct KQUAD <: AbstractElement
#     name::String  = "Quad"                              # element name  
#     PolynomA::Array{Float64,1} = [0.0, 0.0, 0.0, 0.0]    
#     PolynomB::Array{Float64,1} = [0.0, 0.0, 0.0, 0.0]
#     MaxOrder::Int64 = 1
#     NumIntSteps::Int64 = 10
#     FringeQuadEntrance::Int64 = 0
#     FringeQuadExit::Int64 = 0
#     FringeIntM0::Array{Float64,1} = [0.0, 0.0, 0.0, 0.0, 0.0]
#     FringeIntP0::Array{Float64,1} = [0.0, 0.0, 0.0, 0.0, 0.0]
#     T1::Array{Float64,1} = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
#     T2::Array{Float64,1} = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
#     R1::Array{Float64,2} = [1.0 0.0 0.0 0.0 0.0 0.0; 
#                             0.0 1.0 0.0 0.0 0.0 0.0; 
#                             0.0 0.0 1.0 0.0 0.0 0.0; 
#                             0.0 0.0 0.0 1.0 0.0 0.0; 
#                             0.0 0.0 0.0 0.0 1.0 0.0; 
#                             0.0 0.0 0.0 0.0 0.0 1.0]
#     R2::Array{Float64,2} = [1.0 0.0 0.0 0.0 0.0 0.0;
#                             0.0 1.0 0.0 0.0 0.0 0.0;
#                             0.0 0.0 1.0 0.0 0.0 0.0;
#                             0.0 0.0 0.0 1.0 0.0 0.0;
#                             0.0 0.0 0.0 0.0 1.0 0.0;
#                             0.0 0.0 0.0 0.0 0.0 1.0]         
#     RApertures::Array{Float64,1} = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
#     EApertures::Array{Float64,1} = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
#     KickAngle::Array{Float64,1} = [0.0, 0.0]
# end

# using Enzyme
# # q = KQUAD(PolynomialB=[0.0, 1.0, 0.0, 0.0])
# function f(k)
# particles = [0.001 0.0001 0.0005 0.0002 0.0 0.0 0.001 0.0 0.0 0.0 0.0 0.0]
# T1 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
# T2 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
# R1 = [1.0 0.0 0.0 0.0 0.0 0.0; 
#       0.0 1.0 0.0 0.0 0.0 0.0; 
#       0.0 0.0 1.0 0.0 0.0 0.0; 
#       0.0 0.0 0.0 1.0 0.0 0.0; 
#       0.0 0.0 0.0 0.0 1.0 0.0; 
#       0.0 0.0 0.0 0.0 0.0 1.0]
# R2 = [1.0 0.0 0.0 0.0 0.0 0.0;
#         0.0 1.0 0.0 0.0 0.0 0.0;
#         0.0 0.0 1.0 0.0 0.0 0.0;
#         0.0 0.0 0.0 1.0 0.0 0.0;
#         0.0 0.0 0.0 0.0 1.0 0.0;
#         0.0 0.0 0.0 0.0 0.0 1.0]
# PolynomB = [0.0, k[1], 0.0, 0.0]
# StrMPoleSymplectic4Pass!(particles, 1.0, [0.0, 0.0, 0.0, 0.0], PolynomB, 1, 10, 0, 0, 
# [0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0], T1, T2, R1, R2, [0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 
# [0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0], 2)
# # println(particles)
# return particles[1]
# end
# print(f([1.0]))
# grad = gradient(Forward, f, [1.0])
# println(grad)