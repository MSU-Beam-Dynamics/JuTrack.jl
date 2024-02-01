function ATmultmv!(r::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}, A::Matrix{Float64}) where {T, TPS_Dim, Max_TPS_Degree}
    # multiplies 6-component column vector r by 6x6 matrix R: as in A*r
    temp = Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}(undef, 6)

    for i in 1:6
        temp[i] = CTPS(0.0, TPS_Dim, Max_TPS_Degree)
        for j in 1:6
            temp[i] += A[i, j] * r[j]
        end
    end
    for i in 1:6
        r[i] = temp[i]
    end
    return nothing
end

function ATaddvv!(r::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}, dr::Array{Float64, 1}) where {T, TPS_Dim, Max_TPS_Degree}
    for i in 1:6
        r[i] += dr[i]
    end
    return nothing
end

function fastdrift!(r::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}, NormL::CTPS{T, TPS_Dim, Max_TPS_Degree}) where {T, TPS_Dim, Max_TPS_Degree}
    # NormL=(Physical Length)/(1+delta)  is computed externally to speed up calculations
    # in the loop if momentum deviation (delta) does not change
    # such as in 4-th order symplectic integrator w/o radiation
    # r[1] = tadd(r[1], tmult(NormL, r[2]))
    # r[3] = tadd(r[3], tmult(NormL, r[4]))
    # r[6] = tadd(r[6], tmult(NormL, tdiv(tadd(tmult(r[2], r[2]), tmult(r[4], r[4])), tmult(2.0, tadd(1.0, r[5])))))    
    r[1] += NormL * r[2]
    r[3] += NormL * r[4]
    r[6] += NormL * (r[2]^2 + r[4]^2) / (2.0*(1.0+r[5]))
    return nothing 
end

function drift6!(r::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}, le::Float64) where {T, TPS_Dim, Max_TPS_Degree}
    p_norm = 1.0 / (1.0 + r[5])
    NormL = le * p_norm
    # r[1] = tadd(r[1], tmult(NormL, r[2]))
    # r[3] = tadd(r[3], tmult(NormL, r[4]))
    # r[6] = tadd(r[6], tmult(NormL, tdiv(tadd(tmult(r[2], r[2]), tmult(r[4], r[4])), tmult(2.0, tadd(1.0, r[5])))))
    r[1] += NormL * r[2]
    r[3] += NormL * r[4]
    r[6] += NormL * p_norm * (r[2] * r[2] + r[4] * r[4]) / 2.0
    return nothing
end 
function DriftPass_TPSA!(r_in::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}, le, T1, T2, R1, R2, 
    RApertures, EApertures, num_particles) where {T, TPS_Dim, Max_TPS_Degree}
    # Threads.@threads for c in 1:num_particles
    for c in 1:num_particles
        # if !isnan(r_in[1].map[1])
            # Misalignment at entrance
            if !isnothing(T1)
                ATaddvv!(r_in, T1)
            end
            if !isnothing(R1)
                ATmultmv!(r_in, R1)
            end

            drift6!(r_in, le)

            # Misalignment at exit
            if !isnothing(R2)
                ATmultmv!(r_in, R2)
            end
            if !isnothing(T2)
                ATaddvv!(r_in, T2)
            end
        # end
    end
    return nothing
end

function pass_TPSA!(ele::DRIFT, r_in::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}, num_particles::Int64) where {T, TPS_Dim, Max_TPS_Degree}
    # ele: EDRIFT
    # r_in: 6-by-num_particles array
    # num_particles: number of particles

    DriftPass_TPSA!(r_in, ele.len, ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, num_particles)
    return nothing
end

# function f(xx)
# x = CTPS(0.0, 1, 6, 3)
# xp = CTPS(0.0, 2, 6, 3)
# y = CTPS(0.0, 3, 6, 3)
# yp = CTPS(0.0, 4, 6, 3)
# z = CTPS(0.0, 5, 6, 3)
# delta = CTPS(0.0, 6, 6, 3)
# rin = [x, xp, y, yp, z, delta]
# D1 = DRIFT(name="D1",len=xx[1])
# pass_TPSA!(D1, rin, 1)
# println(rin[1])
# return rin[1][3]
# end
# using Enzyme
# grad = gradient(Forward, f, [2.0])
# println(grad)