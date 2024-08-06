function pass_TPSA!(ele::CRABCAVITY, r_in::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}) where {T, TPS_Dim, Max_TPS_Degree}
    E = ele.energy
    gamma = E / m_e
    beta = sqrt(1.0 - 1.0 / gamma^2)
    beta2E = beta^2 * E
    if ele.len == 0.0
        # ang = -ele.k * r_in[5] + ele.phi
        # r_in[2] += ele.volt/beta2E * sin(ang)
        # r_in[6] += (-ele.k * ele.volt/beta2E) * r_in[1] * cos(ang)
        # println("crabcavity is not implemented in TPSA")
        return nothing
    else
        drift6!(r_in, ele.len / 2.0)
        # ang = -ele.k * r_in[5] + ele.phi
        # r_in[2] += ele.volt/beta2E * sin(ang)
        # r_in[6] += (-ele.k * ele.volt/beta2E) * r_in[1] * cos(ang)
        # println("crabcavity is not implemented in TPSA")
        drift6!(r_in, ele.len / 2.0)
    end
    return nothing
end

function pass_TPSA!(ele::easyCRABCAVITY, r_in::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}) where {T, TPS_Dim, Max_TPS_Degree}
    if ele.len == 0.0
        return nothing
    else
        drift6!(r_in, ele.len)
    end
    println("warning: easyCRABCAVITY is not implemented in TPSA")
    return nothing
end

function pass_TPSA!(ele::AccelCavity, r_in::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}) where {T, TPS_Dim, Max_TPS_Degree}
    if ele.len == 0.0
        return nothing
    else
        drift6!(r_in, ele.len)
    end
    println("warning: AccelCavity is not implemented in TPSA")
    return nothing
end