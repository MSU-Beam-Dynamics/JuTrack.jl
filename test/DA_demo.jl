# include("../src/JuTrack.jl")

using JuTrack
using Serialization
using Plots

ESR = deserialize("src/demo/ESR/esr_main_rad.jls")
ESR_norad = deserialize("src/demo/ESR/esr_main.jls")
E0 = 17.846262619763e9

na = 13
nturns = 1000
DA1, survived1 = dynamic_aperture(ESR, nturns, 0.003, 0.0001, na, E0)
DA2, survived2 = dynamic_aperture(ESR_norad, nturns, 0.003, 0.0001, na, E0)


plot(DA1[:,1], DA1[:,2], label="With Radiation", xlabel="x(m)", ylabel="y(m)")
plot!(DA2[:,1], DA2[:,2], label="Without Radiation", xlabel="x(m)", ylabel="y(m)")