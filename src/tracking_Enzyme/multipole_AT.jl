function fastdrift!(r, NormL)
# NormL=(Physical Length)/(1+delta)  is computed externally to speed up calculations
# in the loop if momentum deviation (delta) does not change
# such as in 4-th order symplectic integrator w/o radiation
    r[1] += NormL * r[2]
    r[3] += NormL * r[4]
    r[6] += NormL * (r[2]^2 + r[4]^2) / (2*(1+r[5]))
end

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

function ATaddvv!(r, dr)
    for i in 1:6
        r[i] += dr[i]
    end
end

function ATmultmv!(r, A)
    # multiplies 6-component column vector r by 6x6 matrix R: as in A*r
    temp = zeros(6)
    for i in 1:6
        for j in 1:6
            temp[i] += A[i, j] * r[j]
        end
    end
    for i in 1:6
        r[i] = temp[i]
    end
end

function QuadFringePassP!(r, b2)
    u = b2 / (12.0 * (1.0 + r[5]))
    x2 = r[1]^2
    z2 = r[3]^2
    xz = r[1] * r[3]
    gx = u * (x2 + 3 * z2) * r[1]
    gz = u * (z2 + 3 * x2) * r[3]
    r1tmp = 0.0
    r3tmp = 0.0

    r[1] += gx
    r1tmp = 3 * u * (2 * xz * r[4] - (x2 + z2) * r[2])
    
    r[3] -= gz
    
    r3tmp = 3 * u * (2 * xz * r[2] - (x2 + z2) * r[4])
    r[6] -= (gz * r[4] - gx * r[2]) / (1 + r[5])
    
    r[2] += r1tmp
    r[4] -= r3tmp
end

function QuadFringePassN!(r, b2)
    u = b2 / (12.0 * (1.0 + r[5]))
    x2 = r[1]^2
    z2 = r[3]^2
    xz = r[1] * r[3]
    gx = u * (x2 + 3 * z2) * r[1]
    gz = u * (z2 + 3 * x2) * r[3]
    r1tmp = 0.0
    r3tmp = 0.0

    r[1] -= gx
    r1tmp = 3 * u * (2 * xz * r[4] - (x2 + z2) * r[2])
    
    r[3] += gz
    
    r3tmp = 3 * u * (2 * xz * r[2] - (x2 + z2) * r[4])
    r[6] += (gz * r[4] - gx * r[2]) / (1 + r[5])
    
    r[2] -= r1tmp
    r[4] += r3tmp
end

function quadPartialFringeMatrix!(R, K1, inFringe, fringeInt, part)
    J1x, J2x, J3x, J1y, J2y, J3y = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    K1sqr = K1 * K1

    R[5, 5] = R[6, 6] = 1.0

    if part == 1
        J1x = inFringe * (K1 * fringeInt[2] - 2 * K1sqr * fringeInt[4] / 3)
        J2x = inFringe * (K1 * fringeInt[3])
        J3x = inFringe * (K1sqr * (fringeInt[3] + fringeInt[5]))

        K1 = -K1
        J1y = inFringe * (K1 * fringeInt[2] - 2 * K1sqr * fringeInt[4] / 3)
        J2y = -J2x
        J3y = J3x
    else
        J1x = inFringe * (K1 * fringeInt[2] + K1sqr * fringeInt[1] * fringeInt[3] / 2)
        J2x = inFringe * (K1 * fringeInt[3])
        J3x = inFringe * (K1sqr * (fringeInt[5] - fringeInt[1] * fringeInt[2]))

        K1 = -K1
        J1y = inFringe * (K1 * fringeInt[2] + K1sqr * fringeInt[1] * fringeInt[3])
        J2y = -J2x
        J3y = J3x
    end
    expJ1x = exp(J1x)
    R[1, 1] = expJ1x
    R[1, 2] = J2x / expJ1x
    R[2, 1] = expJ1x * J3x
    R[2, 2] = (1 + J2x * J3x) / expJ1x

    expJ1y = exp(J1y)
    R[3, 3] = expJ1y
    R[3, 4] = J2y / expJ1y
    R[4, 3] = expJ1y * J3y
    R[4, 4] = (1 + J2y * J3y) / expJ1y
end

function linearQuadFringeElegantEntrance!(r6, b2, fringeIntM0, fringeIntP0)
    R = zeros(6, 6)
    inFringe = -1.0
    fringeIntM = fringeIntP0
    fringeIntP = fringeIntM0
    delta = r6[5]

    # determine first linear matrix for this delta
    quadPartialFringeMatrix!(R, b2 / (1 + delta), inFringe, fringeIntM, 1)
    r6[1] = R[1, 1] * r6[1] + R[1, 2] * r6[2]
    r6[2] = R[2, 1] * r6[1] + R[2, 2] * r6[2]
    r6[3] = R[3, 3] * r6[3] + R[3, 4] * r6[4]
    r6[4] = R[4, 3] * r6[3] + R[4, 4] * r6[4]

    # nonlinear fringe field
    QuadFringePassP!(r6, b2)  

    # determine and apply second linear matrix, from elegant code
    quadPartialFringeMatrix!(R, b2 / (1 + delta), inFringe, fringeIntP, 2)
    r6[1] = R[1, 1] * r6[1] + R[1, 2] * r6[2]
    r6[2] = R[2, 1] * r6[1] + R[2, 2] * r6[2]
    r6[3] = R[3, 3] * r6[3] + R[3, 4] * r6[4]
    r6[4] = R[4, 3] * r6[3] + R[4, 4] * r6[4]
end

function linearQuadFringeElegantExit!(r6, b2, fringeIntM0, fringeIntP0)
    R = zeros(6, 6)
    inFringe = 1.0
    fringeIntM = fringeIntM0
    fringeIntP = fringeIntP0
    delta = r6[5]

    # Determine first linear matrix for this delta
    quadPartialFringeMatrix!(R, b2 / (1 + delta), inFringe, fringeIntM, 1)
    r6[1] = R[1, 1] * r6[1] + R[1, 2] * r6[2]
    r6[2] = R[2, 1] * r6[1] + R[2, 2] * r6[2]
    r6[3] = R[3, 3] * r6[3] + R[3, 4] * r6[4]
    r6[4] = R[4, 3] * r6[3] + R[4, 4] * r6[4]

    # Nonlinear fringe field 
    QuadFringePassN!(r6, b2)

    # Determine and apply second linear matrix
    quadPartialFringeMatrix!(R, b2 / (1 + delta), inFringe, fringeIntP, 2)
    r6[1] = R[1, 1] * r6[1] + R[1, 2] * r6[2]
    r6[2] = R[2, 1] * r6[1] + R[2, 2] * r6[2]
    r6[3] = R[3, 3] * r6[3] + R[3, 4] * r6[4]
    r6[4] = R[4, 3] * r6[3] + R[4, 4] * r6[4]
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

using Base: @kwdef
abstract type AbstractElement end

@kwdef struct KQUAD <: AbstractElement
    name::String  = "Quad"                              # element name  
    PolynomA::Array{Float64,1} = [0.0, 0.0, 0.0, 0.0]    
    PolynomB::Array{Float64,1} = [0.0, 0.0, 0.0, 0.0]
    MaxOrder::Int64 = 1
    NumIntSteps::Int64 = 10
    FringeQuadEntrance::Int64 = 0
    FringeQuadExit::Int64 = 0
    FringeIntM0::Array{Float64,1} = [0.0, 0.0, 0.0, 0.0, 0.0]
    FringeIntP0::Array{Float64,1} = [0.0, 0.0, 0.0, 0.0, 0.0]
    T1::Array{Float64,1} = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    T2::Array{Float64,1} = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    R1::Array{Float64,2} = [1.0 0.0 0.0 0.0 0.0 0.0; 
                            0.0 1.0 0.0 0.0 0.0 0.0; 
                            0.0 0.0 1.0 0.0 0.0 0.0; 
                            0.0 0.0 0.0 1.0 0.0 0.0; 
                            0.0 0.0 0.0 0.0 1.0 0.0; 
                            0.0 0.0 0.0 0.0 0.0 1.0]
    R2::Array{Float64,2} = [1.0 0.0 0.0 0.0 0.0 0.0;
                            0.0 1.0 0.0 0.0 0.0 0.0;
                            0.0 0.0 1.0 0.0 0.0 0.0;
                            0.0 0.0 0.0 1.0 0.0 0.0;
                            0.0 0.0 0.0 0.0 1.0 0.0;
                            0.0 0.0 0.0 0.0 0.0 1.0]         
    RApertures::Array{Float64,1} = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    EApertures::Array{Float64,1} = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    KickAngle::Array{Float64,1} = [0.0, 0.0]
end

using Enzyme
# q = KQUAD(PolynomialB=[0.0, 1.0, 0.0, 0.0])
function f(k)
particles = [0.001 0.0001 0.0005 0.0002 0.0 0.0 0.001 0.0 0.0 0.0 0.0 0.0]
T1 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
T2 = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
R1 = [1.0 0.0 0.0 0.0 0.0 0.0; 
      0.0 1.0 0.0 0.0 0.0 0.0; 
      0.0 0.0 1.0 0.0 0.0 0.0; 
      0.0 0.0 0.0 1.0 0.0 0.0; 
      0.0 0.0 0.0 0.0 1.0 0.0; 
      0.0 0.0 0.0 0.0 0.0 1.0]
R2 = [1.0 0.0 0.0 0.0 0.0 0.0;
        0.0 1.0 0.0 0.0 0.0 0.0;
        0.0 0.0 1.0 0.0 0.0 0.0;
        0.0 0.0 0.0 1.0 0.0 0.0;
        0.0 0.0 0.0 0.0 1.0 0.0;
        0.0 0.0 0.0 0.0 0.0 1.0]
PolynomB = [0.0, k[1], 0.0, 0.0]
StrMPoleSymplectic4Pass!(particles, 1.0, [0.0, 0.0, 0.0, 0.0], PolynomB, 1, 10, 0, 0, 
[0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0], T1, T2, R1, R2, [0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 
[0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0], 2)
# println(particles)
return particles[1]
end
grad = gradient(Forward, f, [1.0])
println(grad)