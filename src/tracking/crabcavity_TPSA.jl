function pass_TPSA!(ele::CRABCAVITY, r_in::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}; E0::Float64=0.0, m0::Float64=m_e) where {T, TPS_Dim, Max_TPS_Degree}
    if E0 == 0.0
        println("Warning: beam energy is not defined")
    end
    gamma = (E0 + m0) / m0
    beta = sqrt(1.0 - 1.0 / gamma^2)
    if use_exact_beti == 1
        beti = 1.0 / beta
    else
        beti = 1.0 
    end
    beta2E = beta^2 * E0
    ang = ele.k * r_in[5] + ele.phi
    if ele.len == 0.0
        # r_in[2] += ele.volt/beta2E * sin(ang)
        # r_in[6] += (-ele.k * ele.volt/beta2E) * r_in[1] * cos(ang)
        # this is symplectic form
        r_in[2] += (ele.volt/E0) * sin(ang/beta)
        r_in[6] += (-ele.k * ele.volt/E0/beta) * r_in[1] * cos(ang/beta)
        return nothing
    else
        drift6!(r_in, ele.len / 2.0, beti)
        # r_in[2] += ele.volt/beta2E * sin(ang)
        # r_in[6] += (-ele.k * ele.volt/beta2E) * r_in[1] * cos(ang)
        # this is symplectic form
        r_in[2] += (ele.volt/E0) * sin(ang/beta)
        r_in[6] += (-ele.k * ele.volt/E0/beta) * r_in[1] * cos(ang/beta)
        # println("crabcavity is not implemented in TPSA")
        drift6!(r_in, ele.len / 2.0, beti)
    end
    return nothing
end

function pass_TPSA!(ele::CRABCAVITY_K2, r_in::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}; E0::Float64=0.0, m0::Float64=m_e) where {T, TPS_Dim, Max_TPS_Degree}
    if E0 == 0.0
        println("Warning: beam energy is not defined")
    end
    gamma = (E0 + m0) / m0
    beta = sqrt(1.0 - 1.0 / gamma^2)
    ang = ele.k * r_in[5] + ele.phi
    if use_exact_beti == 1
        beti = 1.0 / beta
    else
        beti = 1.0 
    end
    if ele.len == 0.0
        # r_in[2] += ele.volt/beta2E * sin(ang)
        # r_in[6] += (-ele.k * ele.volt/beta2E) * r_in[1] * cos(ang)
        # this is symplectic form
        r_in[2] += (ele.volt/E0) * sin(ang/beta)
        r_in[6] += (-ele.k * ele.volt/E0/beta) * r_in[1] * cos(ang/beta)
        # sextupole kick
        r_in[2] -= ele.k2 * (r_in[1]^2 - r_in[3]^2) * sin(ang)
        r_in[4] += 2.0 * ele.k2 * r_in[1] * r_in[3] * sin(ang)
        r_in[6] -= (-ele.k2 * ele.k / 3.0) * (r_in[1]^3 - 3.0 * r_in[1] * r_in[3]^2) * cos(ang)
        return nothing
    else
        drift6!(r_in, ele.len / 2.0, beti)
        # r_in[2] += ele.volt/beta2E * sin(ang)
        # r_in[6] += (-ele.k * ele.volt/beta2E) * r_in[1] * cos(ang)
        # this is symplectic form
        r_in[2] += (ele.volt/E0) * sin(ang/beta)
        r_in[6] += (-ele.k * ele.volt/E0/beta) * r_in[1] * cos(ang/beta)
        # sextupole kick
        r_in[2] -= ele.k2 * (r_in[1]^2 - r_in[3]^2) * sin(ang)
        r_in[4] += 2.0 * ele.k2 * r_in[1] * r_in[3] * sin(ang)
        r_in[6] -= (-ele.k2 * ele.k / 3.0) * (r_in[1]^3 - 3.0 * r_in[1] * r_in[3]^2) * cos(ang)
        drift6!(r_in, ele.len / 2.0, beti)
    end
    return nothing
end

function pass_TPSA!(ele::easyCRABCAVITY, r_in::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}; E0::Float64=0.0, m0::Float64=m_e) where {T, TPS_Dim, Max_TPS_Degree}
    gamma = (E0 + m0) / m0
    beta = sqrt(1.0 - 1.0 / gamma^2)
    if use_exact_beti == 1
        beti = 1.0 / beta
    else
        beti = 1.0 
    end
    if ele.len == 0.0
        return nothing
    else
        drift6!(r_in, ele.len, beti)
    end
    println("warning: easyCRABCAVITY is not implemented in TPSA")
    return nothing
end

function pass_TPSA!(ele::AccelCavity, r_in::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}; E0::Float64=0.0, m0::Float64=m_e) where {T, TPS_Dim, Max_TPS_Degree}
    gamma = (E0 + m0) / m0
    beta = sqrt(1.0 - 1.0 / gamma^2)
    if use_exact_beti == 1
        beti = 1.0 / beta
    else
        beti = 1.0 
    end
    if ele.len == 0.0
        return nothing
    else
        drift6!(r_in, ele.len, beti)
    end
    println("warning: AccelCavity is not implemented in TPSA")
    return nothing
end