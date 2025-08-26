function strthinkick1!(r::AbstractVector{Float64}, A, B, L, max_order)
    # Modified based on AT function. Ref[Terebilo, Andrei. "Accelerator modeling with MATLAB accelerator toolbox." PACS2001 (2001)].

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
    return nothing
end
function ThinMPolePass!(r::Array{Float64,1}, le::Float64, A::Array{Float64,1}, B::Array{Float64,1}, 
    # bax::Float64, bay::Float64,
    max_order::Int, 
    T1::Array{Float64,1}, T2::Array{Float64,1}, R1::Array{Float64,2}, R2::Array{Float64,2}, 
    RApertures::Array{Float64,1}, EApertures::Array{Float64,1}, KickAngle::Array{Float64,1}, 
    num_particles::Int, lost_flags::Array{Int64,1})
    # Modified based on AT function. Ref[Terebilo, Andrei. "Accelerator modeling with MATLAB accelerator toolbox." PACS2001 (2001)].

    # no bending
    bax = 0.0 
    bay = 0.0

    B[1] -= KickAngle[1]
    A[1] += KickAngle[2]

    # Threads.@threads for c in 1:num_particles
    for c in 1:num_particles
        if lost_flags[c] == 1
            continue
        end
        r6 = @view r[(c-1)*6+1:c*6]
        if !isnan(r6[1])
            # Misalignment at entrance
            if !iszero(T1)
                addvv!(r6, T1)
            end
            if !iszero(R1)
                multmv!(r6, R1)
            end

            strthinkick1!(r6, A, B, 1.0, max_order)
            r6[2] += bax * r6[6]
            r6[4] -= bay * r6[6]
            r6[6] -= bax * r6[1] - bay * r6[3]  # Path lenghtening

            # Misalignment at exit
            if !iszero(R2)
                multmv!(r6, R2)
            end
            if !iszero(T2)
                addvv!(r6, T2)
            end
            if check_lost(r6)
                lost_flags[c] = 1
            end
        end
    end

    B[1] += KickAngle[1]
    A[1] -= KickAngle[2]
    return nothing
end

function pass!(ele::thinMULTIPOLE, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam{Float64})
    # ele: KQUAD
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    lost_flags = particles.lost_flag
    ThinMPolePass!(r_in, ele.len, ele.PolynomA, ele.PolynomB, ele.MaxOrder, ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, num_particles, lost_flags)
    return nothing
end

function ThinMPolePass_P!(r::Array{Float64,1}, le::Float64, A::Array{Float64,1}, B::Array{Float64,1}, 
    # bax::Float64, bay::Float64,
    max_order::Int, 
    T1::Array{Float64,1}, T2::Array{Float64,1}, R1::Array{Float64,2}, R2::Array{Float64,2}, 
    RApertures::Array{Float64,1}, EApertures::Array{Float64,1}, KickAngle::Array{Float64,1}, 
    num_particles::Int, lost_flags::Array{Int64,1})
    # no bending
    bax = 0.0 
    bay = 0.0

    B[1] -= KickAngle[1]
    A[1] += KickAngle[2]

    Threads.@threads for c in 1:num_particles
    # for c in 1:num_particles
        if lost_flags[c] == 1
            continue
        end
        r6 = @view r[(c-1)*6+1:c*6]
        if !isnan(r6[1])
            # Misalignment at entrance
            if !iszero(T1)
                addvv!(r6, T1)
            end
            if !iszero(R1)
                multmv!(r6, R1)
            end
            
            strthinkick1!(r6, A, B, 1.0, max_order)
            r6[2] += bax * r6[6]
            r6[4] -= bay * r6[6]
            r6[6] -= bax * r6[1] - bay * r6[3]  # Path lenghtening


            # Misalignment at exit
            if !iszero(R2)
                multmv!(r6, R2)
            end
            if !iszero(T2)
                addvv!(r6, T2)
            end
            if check_lost(r6)
                lost_flags[c] = 1
            end
        end
    end

    B[1] += KickAngle[1]
    A[1] -= KickAngle[2]
    return nothing
end

function pass_P!(ele::thinMULTIPOLE, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam{Float64})
    # ele: KQUAD
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    lost_flags = particles.lost_flag
    ThinMPolePass_P!(r_in, ele.len, ele.PolynomA, ele.PolynomB, ele.MaxOrder, ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, num_particles, lost_flags)
    return nothing
end

###########################################################################################
# high-order TPSA
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

function pass_TPSA!(ele::thinMULTIPOLE, r_in::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}; E0::Float64=0.0, m0::Float64=m_e) where {T, TPS_Dim, Max_TPS_Degree}
    # ele: KQUAD
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    ThinMPolePass_TPSA!(r_in, ele.len, ele.PolynomA, ele.PolynomB, ele.MaxOrder, ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle)
    return nothing
end

