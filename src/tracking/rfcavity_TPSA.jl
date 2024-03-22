include("drift_TPSA.jl")

function RFCavityPass!(r_in::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}, le, nv, freq, h, lag, 
    philag, nturn, T0) where {T, TPS_Dim, Max_TPS_Degree}
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

    if le == 0
        # for c in 1:num_particles
            # r6 = @view r_in[(c-1)*6+1:c*6]
            # if !isnan(r6[1])
            # r_in[5] = tminus(r_in[5], tmult(nv, tsin(tminus(tmult(2*pi*freq, tdiv(tminus(r_in[6], lag), C0)), philag))))
            # r_in[6] += -nv * sin(2.0 * pi * freq * ((r_in[5] - lag) / C0 - (h / freq - T0) * nturn) - philag)
            # end
        # end
    else
        # for c in 1:num_particles
            # r6 = @view r_in[(c-1)*6+1:c*6]
            # if !isnan(r6[1])
                # drift-kick-drift
                drift6!(r_in, halflength)

                # nturn = 0
                # r_in[5] = tminus(r_in[5], tmult(nv, tsin(tminus(tmult(2*pi*freq, tdiv(tminus(r_in[6], lag), C0)), philag))))
                # r_in[6] += -nv * sin(2 * pi * freq * ((r_in[5] - lag) / C0 - (h / freq - T0) * nturn) - philag)
                drift6!(r_in, halflength)
                println("rfcavity is not implemented in TPSA")
            # end
        # end
    end
    return nothing
end

function pass_TPSA!(ele::RFCA, r_in::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}) where {T, TPS_Dim, Max_TPS_Degree}
    # ele: RFCA
    # r_in: 6-by-num_particles array
    # num_particles: number of particles
    if ele.energy == 0
        error("Energy is not defined for RFCA ", ele.name)
    end
    T0=1.0/ele.freq      # Does not matter since nturns == 0
    nv = ele.volt / ele.energy
    RFCavityPass!(r_in, ele.len, nv, ele.freq, ele.h, ele.lag, ele.philag, 0, T0)
    return nothing
end

# using Enzyme
# include("../lattice/canonical_elements.jl")
# include("../TPSA_Enzyme/TPSA_fixedmap.jl")
# # # # q = KQUAD(PolynomialB=[0.0, 1.0, 0.0, 0.0])
# function f(xx)
#     q = RFCA(freq=xx[1], volt=xx[2], len=1.0, energy=3.5e9)
#     x = CTPS(0.0, 1, 6, 3)
#     xp = CTPS(0.0, 2, 6, 3)
#     y = CTPS(0.0, 3, 6, 3)
#     yp = CTPS(0.0, 4, 6, 3)
#     z = CTPS(0.0, 5, 6, 3)
#     delta = CTPS(0.0, 6, 6, 3)
#     rin = [x, xp, y, yp, z, delta]
#     pass_TPSA!(q, rin, 1)
# return rin[5].map[6]
# end
# x = [500e6, 1000000.0]
# println(f(x))
# using BenchmarkTools
# @btime grad = gradient(Forward, f, x, Val(2))
# grad = gradient(Reverse, f, x)
# println(grad)