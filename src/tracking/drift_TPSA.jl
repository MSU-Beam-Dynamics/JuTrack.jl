function multmv!(r::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}, A::Matrix{Float64}) where {T, TPS_Dim, Max_TPS_Degree}
    # multiplies 6-component column vector r by 6x6 matrix R: as in A*r
    temp = Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}(undef, 6)

    for i in 1:6
        temp[i] = CTPS(0.0, TPS_Dim, Max_TPS_Degree)
        for j in 1:6
            temp[i] += A[i, j] * r[j]
        end
    end
    r .= temp
    return nothing
end

function addvv!(r::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}, dr::Array{Float64, 1}) where {T, TPS_Dim, Max_TPS_Degree}
    for i in 1:6
        r[i] += dr[i]
    end
    return nothing
end

function fastdrift!(r::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}, NormL::CTPS{T, TPS_Dim, Max_TPS_Degree}, 
    le::Float64) where {T, TPS_Dim, Max_TPS_Degree}
    # similar to AT fastdrift. Provide an option to use exact Hamiltonian or linearized approximation
    # in the loop if momentum deviation (delta) does not change
    # such as in 4-th order symplectic integrator w/o radiation

    if isone(use_exact_Hamiltonian)
        r[1] += NormL * r[2]
        r[3] += NormL * r[4]
        r[5] += NormL * (1.0 + r[6]) - le
    else
        r[1] += NormL * r[2]
        r[3] += NormL * r[4]
        r[5] += NormL * (r[2]^2 + r[4]^2) / (2.0*(1.0+r[5]))
    end
    return nothing 
end

function drift6!(r::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}, le::Float64) where {T, TPS_Dim, Max_TPS_Degree}
    # similar to AT drift6. Provide an option to use exact Hamiltonian or linearized approximation
    if isone(use_exact_Hamiltonian)
        NormL = le / sqrt(((1.0 + r[6])^2 - r[2]^2 - r[4]^2))
        r[5] += NormL * (1.0 + r[6]) - le
    else
        NormL = le / (1.0 + r[6])
        r[5] += NormL * (r[2]^2 + r[4]^2) / (2.0*(1.0+r[5])) # for linearized approximation
    end
    r[1] += NormL * r[2]
    r[3] += NormL * r[4]
    return nothing
end 
function DriftPass_TPSA!(r_in::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}, le, T1, T2, R1, R2, 
    RApertures, EApertures) where {T, TPS_Dim, Max_TPS_Degree}

    if !iszero(T1)
        addvv!(r_in, T1)
    end
    if !iszero(R1)
        multmv!(r_in, R1)
    end

    drift6!(r_in, le)

    # Misalignment at exit
    if !iszero(R2)
        multmv!(r_in, R2)
    end
    if !iszero(T2)
        addvv!(r_in, T2)
    end

    return nothing
end

function pass_TPSA!(ele::DRIFT, r_in::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}) where {T, TPS_Dim, Max_TPS_Degree}
    # ele: EDRIFT
    # r_in: 6-by-num_particles array
    # num_particles: number of particles

    DriftPass_TPSA!(r_in, ele.len, ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures)
    return nothing
end

function pass_TPSA!(ele::MARKER, r_in::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}) where {T, TPS_Dim, Max_TPS_Degree}
    return nothing
end
