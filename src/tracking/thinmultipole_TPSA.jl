include("drift_TPSA.jl")
include("multipole_TPSA.jl")

function ThinMPolePass_TPSA!(r::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}, le::Float64, 
    A::Array{Float64,1}, B::Array{Float64,1}, 
    # bax::Float64, bay::Float64,
    max_order::Int, 
    T1::Array{Float64,1}, T2::Array{Float64,1}, R1::Array{Float64,2}, R2::Array{Float64,2}, 
    RApertures::Array{Float64,1}, EApertures::Array{Float64,1}, KickAngle::Array{Float64,1}) where {T, TPS_Dim, Max_TPS_Degree}
    # no bending
    bax = 0.0 
    bay = 0.0

    B[1] -= KickAngle[1]
    A[1] += KickAngle[2]

    # Threads.@threads for c in 1:num_particles
    # for c in 1:num_particles
    #     if lost_flags[c] == 1
    #         continue
    #     end
    #     r = @view r[(c-1)*6+1:c*6]
        # if !isnan(r[1])
            # Misalignment at entrance
            if T1 != zeros(6)
                ATaddvv!(r, T1)
            end
            if R1 != zeros(6, 6)
                ATmultmv!(r, R1)
            end

            strthinkick!(r, A, B, 1.0, max_order)
            r[2] += bax * r[6]
            r[4] -= bay * r[6]
            r[6] -= bax * r[1] - bay * r[3]  # Path lenghtening


            # Misalignment at exit
            if R2 != zeros(6, 6)
                ATmultmv!(r, R2)
            end
            if T2 != zeros(6)
                ATaddvv!(r, T2)
            end
            # if r[1] > CoordLimit || r[2] > AngleLimit || r[1] < -CoordLimit || r[2] < -AngleLimit
            #     lost_flags[c] = 1
            # end
        # end
    # end

    B[1] += KickAngle[1]
    A[1] -= KickAngle[2]
    return nothing
end

function pass_TPSA!(ele::thinMULTIPOLE, r_in::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}) where {T, TPS_Dim, Max_TPS_Degree}
    # ele: KQUAD
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    ThinMPolePass_TPSA!(r_in, ele.len, ele.PolynomA, ele.PolynomB, ele.MaxOrder, ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle)
    return nothing
end

