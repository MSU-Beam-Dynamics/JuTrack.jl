include("drift_AT.jl")

function RFCavityPass!(r_in, le, nv, freq, h, lag, philag, nturn, T0, num_particles)
    # le - physical length
    # nv - peak voltage (V) normalized to the design enegy (eV)
    # freq - frequency (Hz)
    # h - harmonic number = round(freq * T0)
    # T0 - revolution period (s) = 1.0 / frequency

    C0 = 2.99792458e8  # Speed of light in vacuum

    if le == 0
        for c in 1:num_particles
            r6 = @view r_in[(c-1)*6+1:c*6]
            if !isnan(r6[1])
                r6[5] += -nv * sin(2 * pi * freq * ((r6[6] - lag) / C0 - (h / freq - T0) * nturn) - philag)
            end
        end
    else
        halflength = le / 2
        for c in 1:num_particles
            r6 = @view r_in[(c-1)*6+1:c*6]
            if !isnan(r6[1])
                # drift-kick-drift
                drift6!(r6, halflength)
                r6[5] += -nv * sin(2 * pi * freq * ((r6[6] - lag) / C0 - (h / freq - T0) * nturn) - philag)
                drift6!(r6, halflength)
            end
        end
    end
end


# function f(L, vol, h)
#     ring_length = 432.0
#     freq = h * 299792458.0 / ring_length
#     T0 = ring_length / 299792458.0
#     E = 1e9
#     nv = vol / E
#     particles = [0.001 0.0001 0.0005 0.0002 0.0 0.0 0.001 0.0 0.0 0.0 0.0 0.0]    
#     RFCavityPass!(particles, L, nv, freq, h, 0.0, 0.0, 0, T0, 2)
#     # println(particles)
#     return particles
# end

# print(f(1.0, 4.0e+6, 720))