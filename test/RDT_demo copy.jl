using JuTrack
using Serialization
using Plots
using Printf

# this example calculate the sum of third order resonance terms for the SPEAR3 lattice
# the objective is to minimize the sum of the following terms:  h21000, h10110, h30000, h10200, h10020
# the optimization variables are the sextupole strengths of SDM and SFM
# the optimization is done by gradient descent. The gradient is calculated by autodiff @Enzyme
RING1 = deserialize("src/demo/SPEAR3/spear3.jls")

twi = twissring(RING, 0.0, 1)
beta, alpha, gamma, mu, dp = array_optics(twi)
s = spos(RING)
Sindex = findelem(RING, KSEXT)
plot(s, dp[:, 1], label="dx", xlabel="s (m)", ylabel="Dispersion (m)", title="Dispersion")
# plot!(s, dp[:, 2], label="dpx")
plot!(s[Sindex], dp[Sindex, 1], seriestype="scatter", label="dx (Sextupoles)")

SDM_index = findelem(RING, :name, "SDM") # 0.21m, -17 /m^-3
SFM_index = findelem(RING, :name, "SFM") # 0.21m, 15 /m^-3


function mcf(RING, dp=1e-8, E=3e9)
    orbit = findorbit(RING, dp)
    beam = Beam(orbit, energy=E)
    ringpass!(RING, beam, 1)
    pos = spos(RING)
    L = pos[end]
    z = beam.r[1, 5]
    alpha = z / L / dp
    # gammat = 1.0 / sqrt(alpha)
    return alpha
end

alphaM = mcf(RING)
twi = twissring(RING, 0.0, 1)
beta, alpha, gamma, mu, disp = array_optics(twi)
s = spos(RING)
C = s[end]
betap = -2.0 * alpha
FX = disp[1, 2] .* beta[1, 1] .- disp[1, 1] .* betap[1, 1] ./ 2.0
A2 = -disp[1,1] .* (1 .- cos.(mu[end, 1])) - FX .* sin.(mu[end, 1])
A1 = -1.0 ./ beta[1, 1] .* (disp[1, 1] .* sin.(mu[end, 1]) .- FX .* (1.0 .- cos.(mu[end, 1])) + A2 .* betap[1, 1] ./ 2.0)
mu_sych = mu[:, 1] + (1-cos(mu[end, 1])) * sqrt(())

delta = 1e-3
orbit = findorbit(RING, delta)
x = orbit[1]
xp = orbit[2]
ds = A1.*x + A2.*xp .- alphaM * C * delta

function get_ds(x, xp, delta)
    beam = Beam([x xp 0.0 0.0 0.0 delta], energy=3e9)
    ringpass!(RING, beam, 1)
    return beam.r[1, 5]
end
ds_track = get_ds(x, xp, delta)
# RING[1].volt = 3.2e6

# changednumber = 0
# for i in eachindex(RING)
#     if RING[i] isa QUAD || RING[i] isa KQUAD || RING[i] isa KSEXT || RING[i] isa KOCT || RING[i] isa SBEND
#         RING[i].rad = 1
#         changednumber += 1
#     end
# end

