include("../src/JuTrack.jl")
using. JuTrack
using Serialization
using Enzyme
using BenchmarkTools
Enzyme.API.runtimeActivity!(true)

# circumference = 3834.0018419157
# solen = SOLENOID(len=0.5, ks=0.3)
# Q1 = KQUAD(len=0.5, k1=0.3)
# S1 = KSEXT(len=0.5, k2=0.3)
# K1 = HKICKER(len=0.5, xkick=0.3)
# R1 = RBEND(len=0.2, angle=0.01)
# RF0 = RFCA(name="rf0", len=4.01667, volt=0.0, freq=(299792458*1.0/circumference)*7560.0, h=7560.0, lag=0.5, energy=17.846262619763e9)
# crab = CrabCavity(name="rf_crab", len=4.0, volt=0.0, freq=3.94e8)
# particles = zeros(1, 6)
# particles[1, 1] = 0.001
# particles[1, 2] = 0.0001
# particles_100 = zeros(100, 6)
# particles_100[:, 1] .= 0.001
# particles_100[:, 2] .= 0.0001
# beam_100 = Beam(particles_100)
# beam = Beam(particles)
# line = [crab]
# @btime begin
#     linepass!(line, beam)
# end

# @btime begin
#     linepass!(line, beam_100)
# end

# Load from the file
@time esr = deserialize("esr_main_mutable_vector.jls")
function f(x, ring)
    particles = [0.001 0.0001 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0]
    beam = Beam(particles)
    nameofD1 = ring[2].name
    new_D1 = DRIFT(name=nameofD1, len=x)
    # the ring is pretty long
    new_ring = [ring[1], new_D1, ring[3:end]...]
    linepass!(new_ring, beam)
    return beam.r[1, 1]
end

println(f(0.1, esr))
println(f(0.11, esr))
grad = autodiff(Forward, f, Duplicated(0.1, 1.0), Const(esr))
println(grad)