include("drift_TPSA.jl")
function pass_TPSA!(ele::CrabCavity, r_in::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}) where {T, TPS_Dim, Max_TPS_Degree}
    if ele.len == 0.0
        return nothing
    else
        drift6!(r_in, ele.len)
    end
    println("warning: CrabCavity is not implemented in TPSA")
    return nothing
end

function pass_TPSA!(ele::easyCrabCavity, r_in::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}) where {T, TPS_Dim, Max_TPS_Degree}
    if ele.len == 0.0
        return nothing
    else
        drift6!(r_in, ele.len)
    end
    println("warning: easyCrabCavity is not implemented in TPSA")
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