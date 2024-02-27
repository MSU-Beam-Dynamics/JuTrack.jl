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
@time esr = deserialize("esr_main_vector.jls")
function f(x,k, ring)
    particles = [0.001 0.0001 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0]
    beam = Beam(particles)
    changed_idx = [2, 3]
    new_D1 = DRIFT(len=x)
    lenQ1 = ring[3].len
    new_Q1 = KQUAD(len=lenQ1, k1=k)
    changed_ele = [new_D1, new_Q1]
    ADlinepass!(ring, beam, changed_idx, changed_ele)
    return beam.r[1, 1]
end

println((f(0.1001, -0.22, esr)-f(0.1, -0.22, esr))/0.0001)
println((f(0.1, -0.22, esr)-f(0.1, -0.22001, esr))/0.00001)

dx = [1.0, 1.0]
# grad = autodiff(Forward, f, DuplicatedNoNeed, Duplicated(0.1, 1.0), Duplicated(-0.22, 1.0), Const(esr))
grad = autodiff(Forward, f, BatchDuplicated, BatchDuplicated(0.1, (1.0,0.0)), BatchDuplicated(-0.22, (0.0,1.0)), Const(esr))

println(grad)

# function f(x, y, z)
#     return x^2 + y^2 + z^2
# end
# grad = autodiff(Forward, f, BatchDuplicated, BatchDuplicated(2.0, (1.0,0.0)), BatchDuplicated(3.0, (0.0,1.0)), Const(4.0))
# println(grad)