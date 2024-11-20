
function RFCavityPass!(r_in::Array{Float64,1}, le::Float64, nv::Float64, freq::Float64, h::Float64, 
    lag::Float64, philag::Float64, nturn::Int, T0::Float64, beta::Float64, num_particles::Int, lost_flags::Array{Int64,1})
    # Modified based on AT function. Ref[Terebilo, Andrei. "Accelerator modeling with MATLAB accelerator toolbox." PACS2001 (2001)].
    # le - physical length
    # nv - peak voltage (V) normalized to the design enegy (eV)
    # freq - frequency (Hz)
    # h - harmonic number = round(freq * T0)
    # T0 - revolution period (s) 

    C0 = 2.99792458e8  # Speed of light in vacuum

    if le == 0
        for c in 1:num_particles
            if lost_flags[c] == 1
                continue
            end
            r6 = @view r_in[(c-1)*6+1:c*6]
            if !isnan(r6[1])
                if isinf(2 * pi * freq * ((r6[5] - lag) / C0 - (h / freq - T0) * nturn) - philag)
                    println(c)
                    println(r6)
                end
                r6[6] += -nv * sin(2 * pi * freq * ((r6[5] - lag) / C0 - (h / freq - T0) * nturn) - philag) / beta^2
            end
        end
    else
        halflength = le / 2
        for c in 1:num_particles
            if lost_flags[c] == 1
                continue
            end
            r6 = @view r_in[(c-1)*6+1:c*6]
            if !isnan(r6[1])
                # drift-kick-drift
                drift6!(r6, halflength)
                r6[6] += -nv * sin(2 * pi * freq * ((r6[5] - lag) / C0 - (h / freq - T0) * nturn) - philag) / beta^2
                # r6[6] += -nv * sin(2 * pi * freq * ((-r6[5] - lag)  - (h / freq - T0) * nturn) - philag)
                drift6!(r6, halflength)
            end
            if check_lost(r6)
                lost_flags[c] = 1
            end
        end
    end
    return nothing
end

function pass!(ele::RFCA, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam)
    # ele: RFCA
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    lost_flags = particles.lost_flag
    T0 = particles.T0
    beta = particles.beta
    nturn = 0 #particles.nturn
    if ele.energy == 0
        error("Energy is not defined for RFCA ", ele.name)
    end
    nv = ele.volt / ele.energy
    RFCavityPass!(r_in, ele.len, nv, ele.freq, ele.h, ele.lag, ele.philag, nturn, T0, beta, num_particles, lost_flags)
    return nothing
end


##########################################################################################
# multi-threading
function RFCavityPass_P!(r_in::Array{Float64,1}, le::Float64, nv::Float64, freq::Float64, h::Float64, 
    lag::Float64, philag::Float64, nturn::Int, T0::Float64, beta::Float64, num_particles::Int, lost_flags::Array{Int64,1})
    # le - physical length
    # nv - peak voltage (V) normalized to the design enegy (eV)
    # freq - frequency (Hz)
    # h - harmonic number = round(freq * T0)
    # T0 - revolution period (s) 

    C0 = 2.99792458e8  # Speed of light in vacuum

    if le == 0
        Threads.@threads for c in 1:num_particles
        # for c in 1:num_particles
            if lost_flags[c] == 1
                continue
            end
            r6 = @view r_in[(c-1)*6+1:c*6]
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
            r6 = @view r_in[(c-1)*6+1:c*6]
            if !isnan(r6[1])
                # drift-kick-drift
                drift6!(r6, halflength)
                r6[6] += -nv * sin(2 * pi * freq * ((r6[5] - lag) / C0 - (h / freq - T0) * nturn) - philag) / beta^2
                drift6!(r6, halflength)
            end
            if check_lost(r6)
                lost_flags[c] = 1
            end
        end
    end
    return nothing
end

function pass_P!(ele::RFCA, r_in::Array{Float64,1}, num_particles::Int64, particles::Beam)
    # ele: RFCA
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    if ele.energy == 0
        error("Energy is not defined for RFCA ", ele.name)
    end
    lost_flags = particles.lost_flag
    T0=1.0/ele.freq      # Does not matter since nturns == 0
    beta = particles.beta
    nv = ele.volt / ele.energy
    RFCavityPass_P!(r_in, ele.len, nv, ele.freq, ele.h, ele.lag, ele.philag, 0, T0, beta, num_particles, lost_flags)
    return nothing
end
