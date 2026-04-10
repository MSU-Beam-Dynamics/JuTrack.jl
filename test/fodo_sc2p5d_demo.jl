# Example of using JuTrack 2.5-D space charge in a periodic FODO cell
using Pkg
Pkg.activate("."); Pkg.instantiate() # change "." to your path of JuTrack.jl
using JuTrack
using Random

use_exact_drift(2)

energy_eV = 1.0e9
mass_eV = m_p
charge_state = 1.0
current_A = 20.0
bunch_length_m = 1.0e-3
nmacro = 20_000
ncells = 1000
sample_every = 100
seed = 1234

D1 = DRIFT_SC2P5D(len = 0.2, xsize = 32, ysize = 32, zsize = 32, pipe_radius = 13e-3, Nsteps = 4)
D2 = DRIFT_SC2P5D(len = 0.4, xsize = 32, ysize = 32, zsize = 32, pipe_radius = 13e-3, Nsteps = 8)
D3 = DRIFT_SC2P5D(len = 0.2, xsize = 32, ysize = 32, zsize = 32, pipe_radius = 13e-3, Nsteps = 4)
Q1 = QUAD_SC2P5D(len = 0.1, k1 = 29.6, xsize = 32, ysize = 32, zsize = 32, pipe_radius = 13e-3, Nsteps = 4)
Q2 = QUAD_SC2P5D(len = 0.1, k1 = -29.6, xsize = 32, ysize = 32, zsize = 32, pipe_radius = 13e-3, Nsteps = 4)
line_sc = [D1, Q1, D2, Q2, D3]

# transverse distribution
distparam = [
    3.677529920673089e-4,
    8.428925532276500e-4,
    -0.828277121044551,
    1.0,
    1.0,
    0.0,
    0.0,
    3.677529304933903e-4,
    8.428931246578997e-4,
    0.828276927537804,
    1.0,
    1.0,
    0.0,
    0.0,
    1.0,
    0.1,
    0.5,
    0.0,
    0.0,
    0.0,
    0.0,
]

particles = Gauss3_Dist(distparam, nmacro; seed = seed)
dz = bunch_length_m / nmacro
z_centers = collect(range(-0.5 * bunch_length_m + 0.5 * dz, step = dz, length = nmacro))
# add longitudinal distribution (uniform in z with zero energy spread)
particles[:, 5] .= z_centers[randperm(MersenneTwister(seed + 1), nmacro)]
particles[:, 6] .= 0.0

gamma = (energy_eV + mass_eV) / mass_eV
beta = sqrt(1.0 - 1.0 / gamma^2)
total_particles = max(1, round(Int, current_A * bunch_length_m / (abs(charge_state) * charge_e * speed_of_light * beta)))

beam = Beam(
    particles,
    energy_eV;
    np = total_particles,
    current = current_A,
    charge = charge_state,
    mass = mass_eV,
)

get_emittance!(beam)
println("Initial RMS emittances [x, y, z] = ", beam.emittance)
println("Tracking JuTrack 2.5-D FODO demo for ", ncells, " cells ...")

for cell in 1:ncells
    linepass!(line_sc, beam)
    if (cell % sample_every == 0) || (cell == ncells)
        get_emittance!(beam)
        println(
            "cell = ", cell,
            ", emit = ", beam.emittance,
            ", lost = ", sum(beam.lost_flag), " / ", beam.nmacro,
        )
    end
end
