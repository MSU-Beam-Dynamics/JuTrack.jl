# TRANSLATION and YROTATION are used to convert the MAD-X lattice files
# The transer map is the same as the one in MAD-X.
function pass!(elem::TRANSLATION, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam)
    for c in 1:num_particles
        if isone(particles.lost_flag[c])
            continue
        end
        r6 = @view r_in[(c-1)*6+1:c*6]
        pz = sqrt(1.0 + 2.0 * r6[6] / particles.beta + r6[6]^2 - r6[2]^2 - r6[4]^2)
        r6[1] -= elem.dx + elem.ds * r6[2] / pz
        r6[3] -= elem.dy + elem.ds * r6[4] / pz
        r6[5] += elem.ds * (1.0/particles.beta + r6[6]) / pz
        if check_lost(r6)
            particles.lost_flag[c] = 1
        end
    end
    return nothing
end

function pass_P!(elem::TRANSLATION, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam)
    Threads.@threads for c in 1:num_particles
        if isone(particles.lost_flag[c])
            continue
        end
        r6 = @view r_in[(c-1)*6+1:c*6]
        pz = sqrt(1.0 + 2.0 * r6[6] / particles.beta + r6[6]^2 - r6[2]^2 - r6[4]^2)
        r6[1] -= elem.dx + elem.ds * r6[2] / pz
        r6[3] -= elem.dy + elem.ds * r6[4] / pz
        r6[5] += elem.ds * (1.0 / particles.beta + r6[6]) / pz
        if check_lost(r6)
            particles.lost_flag[c] = 1
        end
    end
    return nothing
end

function pass_TPSA!(elem::TRANSLATION, r_in::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}; E0::Float64=0.0, m0::Float64=0.0) where {T, TPS_Dim, Max_TPS_Degree}
    gamma = (E0+m0) / m0
    beta = sqrt(1.0 - 1.0 / (gamma^2))
    pz = sqrt(1.0 + 2.0 * r_in[6] / beta + r_in[6]^2 - r_in[2]^2 - r_in[4]^2)
    r_in[1] -= elem.dx + elem.ds * r_in[2] / pz
    r_in[3] -= elem.dy + elem.ds * r_in[4] / pz
    r_in[5] += elem.ds * (1.0 / beta + r_in[6]) / pz
    return nothing
end

# YROTATION
function pass!(elem::YROTATION, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam)
    angle = -elem.angle
    if angle == 0.0
        return nothing
    end
    ca = cos(angle)
    sa = sin(angle)
    ta = tan(angle)
    for c in 1:num_particles
        if isone(particles.lost_flag[c])
            continue
        end
        r6 = @view r_in[(c-1)*6+1:c*6]
        x, px, y, py, t, pt = r6[1], r6[2], r6[3], r6[4], r6[5], r6[6]
        pz = sqrt(1.0 + 2.0 * pt / particles.beta + pt^2 - px^2 - py^2)
        ptt = 1.0 - ta*px/pz
        r6[1] = x/(ca*ptt)
        r6[2] = ca*px + sa*pz
        r6[3] = y + ta*x*py/(pz*ptt)
        r6[5] = t + ta*x*(1.0 / particles.beta+pt)/(pz*ptt)

        if check_lost(r6)
            particles.lost_flag[c] = 1
        end
    end
    return nothing
end

function pass_P!(elem::YROTATION, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam)
    angle = -elem.angle
    if angle == 0.0
        return nothing
    end
    ca = cos(angle)
    sa = sin(angle)
    ta = tan(angle)
    Threads.@threads for c in 1:num_particles
        if isone(particles.lost_flag[c])
            continue
        end
        r6 = @view r_in[(c-1)*6+1:c*6]
        x, px, y, py, t, pt = r6[1], r6[2], r6[3], r6[4], r6[5], r6[6]
        pz = sqrt(1.0 + 2.0 * pt / particles.beta + pt^2 - px^2 - py^2)
        ptt = 1.0 - ta*px/pz
        r6[1] = x/(ca*ptt)
        r6[2] = ca*px + sa*pz
        r6[3] = y + ta*x*py/(pz*ptt)
        r6[5] = t + ta*x*(1.0 / particles.beta+pt)/(pz*ptt)

        if check_lost(r6)
            particles.lost_flag[c] = 1
        end
    end
    return nothing
end

function pass_TPSA!(elem::YROTATION, r_in::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}; E0::Float64=0.0, m0::Float64=0.0) where {T, TPS_Dim, Max_TPS_Degree}
    angle = -elem.angle
    if angle == 0.0
        return nothing
    end
    ca = cos(angle)
    sa = sin(angle)
    ta = tan(angle)
    gamma = (E0+m0) / m0
    beta = sqrt(1.0 - 1.0 / (gamma^2))
    x, px, y, py, t, pt = r_in[1], r_in[2], r_in[3], r_in[4], r_in[5], r_in[6]
    pz = sqrt(1.0 + 2.0 * pt / beta + pt^2 - px^2 - py^2)
    ptt = 1.0 - ta*px/pz
    r_in[1] = x/(ca*ptt)
    r_in[2] = ca*px + sa*pz
    r_in[3] = y + ta*x*py/(pz*ptt)
    r_in[5] = t + ta*x*(1.0 / beta+pt)/(pz*ptt)
    return nothing
end
