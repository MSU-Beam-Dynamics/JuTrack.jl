include("drift.jl")
include("fringe.jl")
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


function StrMPoleSymplectic4Pass_P!(r::Array{Float64,1}, le::Float64, A::Array{Float64,1}, B::Array{Float64,1}, 
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

    Threads.@threads for c in 1:num_particles
    # for c in 1:num_particles
        if lost_flags[c] == 1
            continue
        end
        r6 = @view r[(c-1)*6+1:c*6]
        if !isnan(r6[1])
            # norm = 1.0 / (1.0 + r6[5])
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

function pass_P!(ele::KQUAD, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam)
    # ele: KQUAD
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    lost_flags = particles.lost_flag
    PolynomB = zeros(4)
    if ele.PolynomB[1] == 0.0 && ele.PolynomB[2] == 0.0 && ele.PolynomB[3] == 0.0 && ele.PolynomB[4] == 0.0
        PolynomB[2] = ele.k1
        StrMPoleSymplectic4Pass_P!(r_in, ele.len, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
            ele.FringeQuadEntrance, ele.FringeQuadExit, ele.FringeIntM0, ele.FringeIntP0, 
            ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, num_particles, lost_flags)
    else
        PolynomB[1] = ele.PolynomB[1]
        PolynomB[2] = ele.PolynomB[2]
        PolynomB[3] = ele.PolynomB[3]
        PolynomB[4] = ele.PolynomB[4]
        StrMPoleSymplectic4Pass_P!(r_in, ele.len, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
            ele.FringeQuadEntrance, ele.FringeQuadExit, ele.FringeIntM0, ele.FringeIntP0, 
            ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, num_particles, lost_flags)
    end
    return nothing
end

function pass_P!(ele::KSEXT, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam)
    # ele: KSEXT
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    lost_flags = particles.lost_flag
    PolynomB = zeros(4)
    if ele.PolynomB[1] == 0.0 && ele.PolynomB[2] == 0.0 && ele.PolynomB[3] == 0.0 && ele.PolynomB[4] == 0.0
        PolynomB[3] = ele.k2
        StrMPoleSymplectic4Pass_P!(r_in, ele.len, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
            ele.FringeQuadEntrance, ele.FringeQuadExit, ele.FringeIntM0, ele.FringeIntP0, 
            ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, num_particles, lost_flags)
    else
        PolynomB[1] = ele.PolynomB[1]
        PolynomB[2] = ele.PolynomB[2]
        PolynomB[3] = ele.PolynomB[3]
        PolynomB[4] = ele.PolynomB[4]
        StrMPoleSymplectic4Pass_P!(r_in, ele.len, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
            ele.FringeQuadEntrance, ele.FringeQuadExit, ele.FringeIntM0, ele.FringeIntP0, 
            ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, num_particles, lost_flags)
    end
    return nothing
end

function pass_P!(ele::KOCT, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam)
    # ele: KOCT
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    lost_flags = particles.lost_flag
    PolynomB = zeros(4)
    if ele.PolynomB[1] == 0.0 && ele.PolynomB[2] == 0.0 && ele.PolynomB[3] == 0.0 && ele.PolynomB[4] == 0.0
        PolynomB[4] = ele.k3
        StrMPoleSymplectic4Pass_P!(r_in, ele.len, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
            ele.FringeQuadEntrance, ele.FringeQuadExit, ele.FringeIntM0, ele.FringeIntP0, 
            ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, num_particles, lost_flags)
    else
        PolynomB[1] = ele.PolynomB[1]
        PolynomB[2] = ele.PolynomB[2]
        PolynomB[3] = ele.PolynomB[3]
        PolynomB[4] = ele.PolynomB[4]
        StrMPoleSymplectic4Pass_P!(r_in, ele.len, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
            ele.FringeQuadEntrance, ele.FringeQuadExit, ele.FringeIntM0, ele.FringeIntP0, 
            ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, num_particles, lost_flags)
    end
    return nothing
end

# function matrix_to_array(matrix::Matrix{Float64})
#     particles = zeros(Float64, size(matrix, 1)*size(matrix, 2))
#     for i in 1:size(matrix, 1)
#         for j in 1:size(matrix, 2)
#             particles[(i-1)*size(matrix, 2)+j] = matrix[i, j]
#         end
#     end
#     return particles
# end
# r = zeros(Float64, 1000000, 6)
# r[:, 2] .= 0.1
# r_arr = matrix_to_array(r)
# beam = Beam(r)
# Q = KQUAD(k1=0.1, len=1.0)
# # pass_P!(Q, r_arr, beam.nmacro, beam, zeros(Float64, 6), zeros(Float64, 6, 6))
# # println(r_arr[1:6])
# # println(Threads.nthreads())
# # @btime begin    
# #     pass_P!(Q, r_arr, beam.nmacro, beam, zeros(Float64, 6), zeros(Float64, 6, 6))
# # end

# function f_multithreading(x1, x2)
#     r = zeros(Float64, 1000000, 6)
#     r[:, 2] .= 0.1
#     r_arr = matrix_to_array(r)
#     beam = Beam(r)
#     Q1 = KQUAD(k1=x1, len=1.0)
#     Q2 = KQUAD(k1=x2, len=1.0)
#     pass_P!(Q1, r_arr, beam.nmacro, beam, zeros(Float64, 6), zeros(Float64, 6, 6))
#     pass_P!(Q2, r_arr, beam.nmacro, beam, zeros(Float64, 6), zeros(Float64, 6, 6))
#     return r_arr[1]
# end
# println(f_multithreading(3.0, 1.0))
# # grad = gradient(Forward, f_multithreading, [3.0, 1.0])
# grad = autodiff(Forward, f_multithreading,  Duplicated, Duplicated(3.0, 1.0), Duplicated(1.0, 1.0))
# @time grad = autodiff(Forward, f_multithreading, Duplicated, Duplicated(3.0, 1.0), Duplicated(1.0, 1.0))
# println(grad)