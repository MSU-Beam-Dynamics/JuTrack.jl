include("drift.jl")

function RFCavityPass!(r_in::Array{Float64,1}, le::Float64, nv::Float64, freq::Float64, h::Float64, 
    lag::Float64, philag::Float64, nturn::Int, T0::Int, num_particles::Int, lost_flags::Array{Int64,1})
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
                r6[5] += -nv * sin(2 * pi * freq * ((r6[6] - lag) / C0 - (h / freq - T0) * nturn) - philag)
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
                r6[5] += -nv * sin(2 * pi * freq * ((r6[6] - lag) / C0 - (h / freq - T0) * nturn) - philag)
                drift6!(r6, halflength)
            end
            if r6[1] > CoordLimit || r6[2] > AngleLimit || r6[1] < -CoordLimit || r6[2] < -AngleLimit
                lost_flags[c] = 1
            end
        end
    end
    return nothing
end

function pass!(ele::RFCA, r_in::Array{Float64,1}, num_particles::Int64, lost_flags::Array{Int64,1})
    # ele: RFCA
    # r_in: 6-by-num_particles array
    # num_particles: number of particles

    T0=1.0/ele.freq      # Does not matter since nturns == 0
    nv = ele.volt / ele.energy
    RFCavityPass!(r_in, ele.len, nv, ele.freq, ele.h, ele.lag, ele.philag, 0, T0, num_particles, lost_flags)
    return nothing
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