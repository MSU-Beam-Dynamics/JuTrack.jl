# Ref [Qiang, Ji. "Differentiable self-consistent space-charge simulation for accelerator design." Physical Review Accelerators and Beams 26.2 (2023): 024601.]

function calculate_K(beam, I)
    m0 = beam.mass * 1.78266192e-36 # kg, may result in a slight round-off error
    charge = abs(beam.charge * charge_e)
    K = charge * I / (2.0 * pi * epsilon_0 * m0 * speed_of_light^3 * beam.beta^3 * beam.gamma^3)
    return K
end

# function shape_function(x, deltax)
#     if abs(x) < deltax / 2.0
#         return 3.0 / 4.0 - (x / deltax)^2
#     elseif deltax / 2.0 <= abs(x) < 3.0 * deltax / 2.0
#         return 0.5 * (1.5 - abs(x / deltax))^2
#     else
#         return 0.0
#     end
# end

# function d_shape_function(x, deltax)
#     if abs(x) < deltax / 2.0
#         return -2.0 * x / deltax^2
#     elseif deltax / 2.0 <= abs(x) < 3.0 * deltax / 2.0 && x > 0
#         return (-1.5 + x / deltax) / deltax
#     elseif deltax / 2.0 <= abs(x) < 3.0 * deltax / 2.0 && x <= 0
#         return (1.5 + x / deltax) / deltax
#     else
#         return 0.0
#     end
# end

function compute_delta_pxy!(r_in::Vector{T}, K, Nl, Nm, a, b, Np, dt, lost_flags) where T
    x = r_in[1:6:end]
    y = r_in[3:6:end]
    
    mask = 1 .- lost_flags          # mask[k] = 1 if particle k is active, 0 if lost
    
    xmin = minimum(x)
    ymin = minimum(y)
    x_shifted = x .- xmin           # Shift x positions
    y_shifted = y .- ymin           # Shift y positions
    
    al = (1:Nl) * π ./ a             
    bm = (1:Nm) * π ./ b            
    gamma2lm = al.^2 .+ (bm'.^2)    
    sin_al_xk = sin.(x_shifted * al')      # Np x Nl
    sin_bm_yk = sin.(y_shifted * bm')      # Np x Nm
    
    # Zero out lost particles
    sin_al_xk .= sin_al_xk .* mask      
    sin_bm_yk .= sin_bm_yk .* mask
    
    # φ_lm
    philm = (sin_al_xk' * sin_bm_yk) ./ gamma2lm 
    philm .*= (4.0 * π * 4.0) / (a * b * Np)
    
    cos_al_xi = cos.(x_shifted * al')      # Np x Nl
    cos_bm_yi = cos.(y_shifted * bm')      # Np x Nm
    
    # Apply mask to cos terms
    cos_al_xi .= cos_al_xi .* mask
    cos_bm_yi .= cos_bm_yi .* mask
    
    # Initialize term1 and term2
    term1 = zeros(T, Np)  # For Δp_x
    term2 = zeros(T, Np)  # For Δp_y
    
    # Compute Δp_x and Δp_y
    for l in 1:Nl
        for m in 1:Nm
            # Retrieve α_l and β_m
            α_l = al[l]
            β_m = bm[m]
            φ_lm = philm[l, m]
            # Compute terms for all particles
            term1 .+= φ_lm * α_l .* (cos_al_xi[:, l] .* sin_bm_yk[:, m])
            term2 .+= φ_lm * β_m .* (sin_al_xk[:, l] .* cos_bm_yi[:, m])
        end
    end
    
    # Apply the constant factor
    delta_px = - (dt * K / 2.0) * term1
    delta_py = - (dt * K / 2.0) * term2

    # Update p_x and p_y
    r_in[2:6:end] .+= delta_px
    r_in[4:6:end] .+= delta_py
    return nothing
end

function space_charge_P!(r_in, K, Nl, Nm, dx, dy, a, b, Np, dt, lost_flags)
    # call compute_delta_pxy! to compute Δp_x and Δp_y
    compute_delta_pxy!(r_in, K, Nl, Nm, a, b, Np, dt, lost_flags)
    return nothing
end
function space_charge!(r_in, K, Nl, Nm, dx, dy, a, b, Np, dt, lost_flags)
    # call compute_delta_pxy! to compute Δp_x and Δp_y
    compute_delta_pxy!(r_in, K, Nl, Nm, a, b, Np, dt, lost_flags)
    return nothing
end

function pass!(ele::SPACECHARGE, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam{Float64})
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
    # space_charge!(r_in, K, ele.Nl, ele.Nm, dx, dy, ele.a, ele.b, num_particles, dt, lost_flags)
    compute_delta_pxy!(r_in, K, ele.Nl, ele.Nm, ele.a, ele.b, num_particles, dt, lost_flags)
    return nothing
end

function pass_TPSA!(ele::SPACECHARGE, r_in::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}) where {T, TPS_Dim, Max_TPS_Degree}
    println("TPSA is not implemented for space charge yet.")
    return nothing
end

function pass_P!(ele::SPACECHARGE, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam{Float64})
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
    # space_charge_P!(r_in, K, ele.Nl, ele.Nm, dx, dy, ele.a, ele.b, num_particles, dt, lost_flags)
    compute_delta_pxy!(r_in, K, ele.Nl, ele.Nm, ele.a, ele.b, num_particles, dt, lost_flags)
    return nothing
end
