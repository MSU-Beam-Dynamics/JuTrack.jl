function CorrectorPass_TPSA!(r::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}, le::Float64, xkick::Float64, ykick::Float64,
    T1::Array{Float64,1}, T2::Array{Float64,1}, R1::Array{Float64,2}, R2::Array{Float64,2}) where {T, TPS_Dim, Max_TPS_Degree}
    
    p_norm = 1.0 / (1.0 + r[6])
    NormL = le * p_norm
    # Misalignment at entrance
    if !iszero(T1)
        addvv!(r, T1)
    end
    if !iszero(R1)
        multmv!(r, R1)
    end
    
    r[5] += NormL*p_norm*(xkick*xkick/3.0 + ykick*ykick/3.0 +
               r[2]*r[2] + r[4]*r[4] +
               r[2]*xkick + r[4]*ykick)/2.0
    r[1] += NormL*(r[2]+xkick/2.0)
    r[2] += xkick
    r[3] += NormL*(r[4]+ykick/2.0)
       r[4] += ykick

    # Misalignment at exit
    if !iszero(R2)
        multmv!(r, R2)
    end
    if !iszero(T2)
        addvv!(r, T2)
    end    


    return nothing
end

function pass_TPSA!(ele::CORRECTOR, r_in::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}) where {T, TPS_Dim, Max_TPS_Degree}
    # ele: CORRECTOR
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    CorrectorPass_TPSA!(r_in, ele.len, ele.xkick, ele.ykick, ele.T1, ele.T2, ele.R1, ele.R2,)
    return nothing
end

