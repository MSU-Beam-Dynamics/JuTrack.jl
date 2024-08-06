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

    # Misalignment at entrance
    if !iszero(T1)
        addvv!(r, T1)
    end
    if !iszero(R1)
        multmv!(r, R1)
    end

    strthinkick!(r, A, B, 1.0, max_order)
    r[2] += bax * r[6]
    r[4] -= bay * r[6]
    r[6] -= bax * r[1] - bay * r[3]  # Path lenghtening


    # Misalignment at exit
    if !iszero(R2)
        multmv!(r, R2)
    end
    if !iszero(T2)
        addvv!(r, T2)
    end

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

