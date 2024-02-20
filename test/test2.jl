include("../src/JuTrack.jl")
using. JuTrack

solen = SOLENOID(len=0.5, ks=0.3)
particles = [0.001 0.0001 0.0005 0.0002 0.0 0.0; 0.001 0.0 0.0 0.0 0.0 0.0]
beam = Beam(particles)
line = [solen]
linepass!(line, beam)
println(beam.r)

x = CTPS(0.0, 1, 6, 3)
xp = CTPS(0.0, 2, 6, 3)
y = CTPS(0.0, 3, 6, 3)
yp = CTPS(0.0, 4, 6, 3)
z = CTPS(0.0, 5, 6, 3)
delta = CTPS(0.0, 6, 6, 3)
rin = [x, xp, y, yp, z, delta]
pass_TPSA!(solen, rin)
println(rin[1].map)