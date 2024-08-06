function calculate_K(beam, I)
    m0 = beam.mass * 1.782662e-36 # kg
    charge = abs(beam.charge * charge_e)
    K = charge * I / (2.0 * pi * epsilon_0 * m0 * speed_of_light^3 * beam.beta^3 * beam.gamma^3)
    return K
end

function shape_function(x, deltax)
    if abs(x) < deltax / 2.0
        return 3.0 / 4.0 - (x / deltax)^2
    elseif deltax / 2.0 <= abs(x) < 3.0 * deltax / 2.0
        return 0.5 * (1.5 - abs(x / deltax))^2
    else
        return 0.0
    end
end

function d_shape_function(x, deltax)
    if abs(x) < deltax / 2.0
        return -2.0 * x / deltax^2
    elseif deltax / 2.0 <= abs(x) < 3.0 * deltax / 2.0 && x > 0
        return (-1.5 + x / deltax) / deltax
    elseif deltax / 2.0 <= abs(x) < 3.0 * deltax / 2.0 && x <= 0
        return (1.5 + x / deltax) / deltax
    else
        return 0.0
    end
end

# function calculate_philm(rin, Nl, Nm, dx, dy, a, b, Np)
#     rholm = zeros(Nl, Nm)
#     gamma2lm = zeros(Nl, Nm)
#     philm = zeros(Nl, Nm)
#     for i in 1:Nl
#         for j in 1:Nm
#             for k in 1:Np
#                 grid_x = (i - 1) * dx
#                 grid_y = (j - 1) * dy
#                 rholm[i, j] += shape_function(rin[(k-1)*6 + 1] - grid_x, dx) * shape_function(rin[(k-1)*6 + 3] - grid_y, dy)
#             end
#         end
#     end
#     rholm *= 1.0 / (dx * dy * Np)

#     for i in 1:Nl
#         for j in 1:Nm
#             al = i * pi / a
#             bm = j * pi / b
#             gamma2lm[i, j] = al^2 + bm^2
#             philm[i, j] = 4.0 * pi * rholm[i, j] / gamma2lm[i, j]
#         end
#     end

#     return philm, gamma2lm, rholm
# end
function calculate_philm(rin, Nl, Nm, dx, dy, a, b, Np, lost_flags)
    philm = zeros(Nl, Nm)
    gamma2lm = zeros(Nl, Nm)
    for i in 1:Nl
        for j in 1:Nm
            al = i * pi / a
            bm = j * pi / b
            gamma2lm[i, j] = al^2 + bm^2
            for k in 1:Np
                if lost_flags[k] == 1
                    continue
                end
                philm[i, j] += sin(al * rin[(k-1)*6 + 1]) * sin(bm * rin[(k-1)*6 + 3]) / gamma2lm[i, j]
            end
        end
    end
    philm .*= 4.0 * pi  * 4.0 / (a * b * Np)
    return philm, gamma2lm
end

function calculate_philm_P(rin, Nl, Nm, dx, dy, a, b, Np, lost_flags)
    philm = zeros(Nl, Nm)
    gamma2lm = zeros(Nl, Nm)
    for i in 1:Nl
        for j in 1:Nm
            al = i * pi / a
            bm = j * pi / b
            gamma2lm[i, j] = al^2 + bm^2

            local_sum = zeros(Threads.nthreads())  
        
            Threads.@threads for k in 1:Np
                tid = Threads.threadid()
                if lost_flags[k] == 1
                    continue
                end
                local_sum[tid] += sin(al * rin[(k-1)*6 + 1]) * sin(bm * rin[(k-1)*6 + 3]) / gamma2lm[i, j]
            end
            
            # Sum up all partial results from each thread
            philm[i, j] = sum(local_sum)
        end
    end
    philm .*= 4.0 * pi  * 4.0 / (a * b * Np)
    return philm, gamma2lm
end

function space_charge!(r_in, K, Nl, Nm, dx, dy, a, b, Np, dt, lost_flags)
    philm, gamma2lm = calculate_philm(r_in, Nl, Nm, dx, dy, a, b, Np, lost_flags)
    term1 = zeros(Np)
    term2 = zeros(Np)
    for i in 1:Nl
        for j in 1:Nm
            al = i * pi / a
            bm = j * pi / b
            term1 .+= (philm[i, j] * al .* cos.(al .* r_in[1:6:end]) .* sin.(bm .* r_in[3:6:end]))
            term2 .+= (philm[i, j] * bm .* sin.(al .* r_in[1:6:end]) .* cos.(bm .* r_in[3:6:end]))
        end
    end
    r_in[2:6:end] .-= (dt * K / 2.0) * term1
    r_in[4:6:end] .-= (dt * K / 2.0) * term2
end

function space_charge_P!(r_in, K, Nl, Nm, dx, dy, a, b, Np, dt, lost_flags)
    philm, gamma2lm = calculate_philm_P(r_in, Nl, Nm, dx, dy, a, b, Np, lost_flags)
    term1 = zeros(Np)
    term2 = zeros(Np)

    nthreads = Threads.nthreads()
    term1_thread = [zeros(Np) for _ in 1:nthreads]  
    term2_thread = [zeros(Np) for _ in 1:nthreads]

    Threads.@threads for i in 1:Nl
        tid = Threads.threadid()
        for j in 1:Nm
            al = i * pi / a
            bm = j * pi / b
            term1_thread[tid] .+= (philm[i, j] * al .* cos.(al .* r_in[1:6:end]) .* sin.(bm .* r_in[3:6:end]))
            term2_thread[tid] .+= (philm[i, j] * bm .* sin.(al .* r_in[1:6:end]) .* cos.(bm .* r_in[3:6:end]))
        end
    end

    # Combine results from all threads
    for t in 1:nthreads
        term1 .+= term1_thread[t]
        term2 .+= term2_thread[t]
    end
    r_in[2:6:end] .-= (dt * K / 2.0) * term1
    r_in[4:6:end] .-= (dt * K / 2.0) * term2
end

function pass!(ele::SPACECHARGE, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam)
    # ele: SPACECHARGE
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    lost_flags = particles.lost_flag
    I = particles.current
    dx = ele.a / ele.Nl
    dy = ele.b / ele.Nm
    # v = speed_of_light * sqrt(1.0 - 1.0 / particles.gamma^2)
    dt = ele.effective_len 
    K = calculate_K(particles, I)
    m0 = particles.mass * 1.782662e-36 # kg
    # p0 = particles.gamma * particles.beta * m0 * speed_of_light
    space_charge!(r_in, K, ele.Nl, ele.Nm, dx, dy, ele.a, ele.b, num_particles, dt, lost_flags)
    return nothing
end

function pass_TPSA!(ele::SPACECHARGE, r_in::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}) where {T, TPS_Dim, Max_TPS_Degree}
    println("TPSA is not implemented for space charge yet.")
    return nothing
end

function pass_P!(ele::SPACECHARGE, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam)
    # println("Parallel computing is not implemented for space charge yet.")
    lost_flags = particles.lost_flag
    I = particles.current
    dx = ele.a / ele.Nl
    dy = ele.b / ele.Nm
    # v = speed_of_light * sqrt(1.0 - 1.0 / particles.gamma^2)
    dt = ele.effective_len
    K = calculate_K(particles, I)
    # convert the mass from eV to kg
    m0 = particles.mass * 1.782662e-36 # kg
    # p0 = particles.gamma * particles.beta * m0 * speed_of_light
    space_charge_P!(r_in, K, ele.Nl, ele.Nm, dx, dy, ele.a, ele.b, num_particles, dt, lost_flags)
    return nothing
end
