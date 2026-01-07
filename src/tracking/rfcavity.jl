
function RFCavityPass!(r_in::Matrix{Float64}, le::Float64, nv::Float64, freq::Float64, h::Float64, 
    lag::Float64, philag::Float64, nturn::Int, T0::Float64, beta::Float64, num_particles::Int, lost_flags::Array{Int64,1})
    # Modified based on AT function. Ref[Terebilo, Andrei. "Accelerator modeling with MATLAB accelerator toolbox." PACS2001 (2001)].
    # le - physical length
    # nv - peak voltage (V) normalized to the design enegy (eV)
    # freq - frequency (Hz)
    # h - harmonic number = round(freq * T0)
    # T0 - revolution period (s) 

    C0 = 2.99792458e8  # Speed of light in vacuum
    if use_exact_beti == 1
        beti = 1.0 / beta
    else
        beti = 1.0 
    end
    if le == 0
        for c in 1:num_particles
            if lost_flags[c] == 1
                continue
            end
            r6 = @view r_in[c, :]
            if !isnan(r6[1])
                r6[6] += -nv * sin(2 * pi * freq * ((r6[5] - lag) / C0 - (h / freq - T0) * nturn) - philag) / beta^2
            end
        end
    else
        halflength = le / 2
        for c in 1:num_particles
            if lost_flags[c] == 1
                continue
            end
            r6 = @view r_in[c, :]
            if !isnan(r6[1])
                # drift-kick-drift
                drift6!(r6, halflength, beti)
                r6[6] += -nv * sin(2 * pi * freq * ((r6[5] - lag) / C0 - (h / freq - T0) * nturn) - philag) / beta^2
                # r6[6] += -nv * sin(2 * pi * freq * ((-r6[5] - lag)  - (h / freq - T0) * nturn) - philag)
                drift6!(r6, halflength, beti)
            end
            if check_lost(r6)
                lost_flags[c] = 1
            end
        end
    end
    return nothing
end

function pass!(ele::RFCA, r_in::Matrix{Float64}, num_particles::Int64, particles::Beam{Float64})
    # ele: RFCA
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    lost_flags = particles.lost_flag
    T0 = particles.T0
    beta = particles.beta
    nturn = 0 #particles.nturn
    if ele.energy == 0
        println("Energy is not defined for RFCA ", ele.name)
    end
    nv = ele.volt / ele.energy
    RFCavityPass!(r_in, ele.len, nv, ele.freq, ele.h, ele.lag, ele.philag, nturn, T0, beta, num_particles, lost_flags)
    return nothing
end


##########################################################################################
# multi-threading
function RFCavityPass_P!(r_in::Matrix{Float64}, le::Float64, nv::Float64, freq::Float64, h::Float64, 
    lag::Float64, philag::Float64, nturn::Int, T0::Float64, beta::Float64, num_particles::Int, lost_flags::Array{Int64,1})
    # le - physical length
    # nv - peak voltage (V) normalized to the design enegy (eV)
    # freq - frequency (Hz)
    # h - harmonic number = round(freq * T0)
    # T0 - revolution period (s) 

    C0 = 2.99792458e8  # Speed of light in vacuum
    if use_exact_beti == 1
        beti = 1.0 / beta
    else
        beti = 1.0 
    end
    if le == 0
        Threads.@threads for c in 1:num_particles
        # for c in 1:num_particles
            if lost_flags[c] == 1
                continue
            end
            r6 = @view r_in[c, :]
            if !isnan(r6[1])
                r6[6] += -nv * sin(2 * pi * freq * ((r6[5] - lag) / C0 - (h / freq - T0) * nturn) - philag) / beta^2
            end
        end
    else
        halflength = le / 2
        Threads.@threads for c in 1:num_particles
        # for c in 1:num_particles
            if lost_flags[c] == 1
                continue
            end
            r6 = @view r_in[c, :]
            if !isnan(r6[1])
                # drift-kick-drift
                drift6!(r6, halflength, beti)
                r6[6] += -nv * sin(2 * pi * freq * ((r6[5] - lag) / C0 - (h / freq - T0) * nturn) - philag) / beta^2
                drift6!(r6, halflength, beti)
            end
            if check_lost(r6)
                lost_flags[c] = 1
            end
        end
    end
    return nothing
end

function pass_P!(ele::RFCA, r_in::Matrix{Float64}, num_particles::Int64, particles::Beam{Float64})
    # ele: RFCA
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    if ele.energy == 0
        println("Energy is not defined for RFCA ", ele.name)
    end
    lost_flags = particles.lost_flag
    T0=1.0/ele.freq      # Does not matter since nturns == 0
    beta = particles.beta
    nv = ele.volt / ele.energy
    RFCavityPass_P!(r_in, ele.len, nv, ele.freq, ele.h, ele.lag, ele.philag, 0, T0, beta, num_particles, lost_flags)
    return nothing
end


# TPSA
function RFCavityPass!(r_in::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}, le, nv, freq, h, lag, 
    philag, nturn, T0, beta) where {T, TPS_Dim, Max_TPS_Degree}
    # le - physical length
    # nv - peak voltage (V) normalized to the design enegy (eV)
    # freq - frequency (Hz)
    # h - harmonic number = round(freq * T0)
    # T0 - revolution period (s) 

    # if num_particles != 1
    #     error("RFCavityPass: num_particles must be 1 for TPSA")
    # end
    C0 = 2.99792458e8  # Speed of light in vacuum
    halflength = le / 2.0
    if use_exact_beti == 1
        beti = 1.0 / beta
    else
        beti = 1.0
    end

    if le == 0
        # for c in 1:num_particles
            # r6 = @view r_in[c, :]
            # if !isnan(r6[1])
            # r_in[5] = tminus(r_in[5], tmult(nv, tsin(tminus(tmult(2*pi*freq, tdiv(tminus(r_in[6], lag), C0)), philag))))
            r_in[6] += -nv * sin(2.0 * pi * freq * ((r_in[5] - lag) / C0 - (h / freq - T0) * nturn) - philag) / beta^2
            # end
        # end
        return nothing
    else
        # for c in 1:num_particles
            # r6 = @view r_in[c, :]
            # if !isnan(r6[1])
                # drift-kick-drift
                drift6!(r_in, halflength, beti)

                # nturn = 0
                # r_in[5] = tminus(r_in[5], tmult(nv, tsin(tminus(tmult(2*pi*freq, tdiv(tminus(r_in[6], lag), C0)), philag))))
                r_in[6] += -nv * sin(2 * pi * freq * ((r_in[5] - lag) / C0 - (h / freq - T0) * nturn) - philag) / beta^2
                drift6!(r_in, halflength, beti)
                # println("rfcavity is not implemented in TPSA")
            # end
        # end
    end
    return nothing
end

function pass_TPSA!(ele::RFCA, r_in::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}; E0::Float64=0.0, m0::Float64=m_e) where {T, TPS_Dim, Max_TPS_Degree}
    # ele: RFCA
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    if ele.energy == 0
        println("Energy is not defined for RFCA ", ele.name)
    end
    if E0 == 0.0
        println("Warning: beam energy is not defined")
    end
    # println("Beam energy is assumed ", E0, " eV in the RF cavity", ", mass is ", m0, " eV")
    T0=1.0/ele.freq      # Does not matter since nturns == 0
    nv = ele.volt / ele.energy
    gamma = (E0 + m0) / m0
    beta = sqrt(1.0 - 1.0 / gamma^2)
    RFCavityPass!(r_in, ele.len, nv, ele.freq, ele.h, ele.lag, ele.philag, 0, T0, beta)
    return nothing
end
