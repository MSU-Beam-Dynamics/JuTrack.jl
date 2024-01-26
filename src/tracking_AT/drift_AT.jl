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

function ATaddvv!(r, dr)
    for i in 1:6
        r[i] += dr[i]
    end
end

function fastdrift!(r, NormL)
    # NormL=(Physical Length)/(1+delta)  is computed externally to speed up calculations
    # in the loop if momentum deviation (delta) does not change
    # such as in 4-th order symplectic integrator w/o radiation
        r[1] += NormL * r[2]
        r[3] += NormL * r[4]
        r[6] += NormL * (r[2]^2 + r[4]^2) / (2*(1+r[5]))
end

function drift6!(r, le)
    p_norm = 1.0 / (1.0 + r[5])
    NormL = le * p_norm
    r[1] += NormL * r[2]
    r[3] += NormL * r[4]
    r[6] += NormL * p_norm * (r[2] * r[2] + r[4] * r[4]) / 2.0
end 
function DriftPass!(r_in, le, T1, T2, R1, R2, RApertures, EApertures, num_particles)
    Threads.@threads for c in 1:num_particles
        r6 = @view r_in[(c-1)*6+1:c*6]
        if !isnan(r6[1])
            # Misalignment at entrance
            if !isnothing(T1)
                ATaddvv!(r6, T1)
            end
            if !isnothing(R1)
                ATmultmv!(r6, R1)
            end
            # Check physical apertures at the entrance of the magnet
            # if RApertures !== nothing
            #     checkiflostRectangularAp!(r6, RApertures)
            # end
            # if EApertures !== nothing
            #     checkiflostEllipticalAp!(r6, EApertures)
            # end
            drift6!(r6, le)
            # Check physical apertures at the exit of the magnet
            # if RApertures !== nothing
            #     checkiflostRectangularAp!(r6, RApertures)
            # end
            # if EApertures !== nothing
            #     checkiflostEllipticalAp!(r6, EApertures)
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
end


# function f_test(L)
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
#     DriftPass!(particles, L[1], T1, T2, R1, R2, nothing, nothing, 2)
#     # println(particles)
#     return particles
#     end
#     print(f_test([1.0]))