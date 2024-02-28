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


function f1(x, ring)
    particles = [0.001 0.0001 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0]
    beam = Beam(particles)
    changed_idx = [2]
    new_D1 = DRIFT(len=x)
    # lenQ1 = ring[3].len
    # new_Q1 = KQUAD(len=lenQ1, k1=k)
    changed_ele = [new_D1]
    ADlinepass!(ring, beam, changed_idx, changed_ele)
    return beam.r[1, 1]
end
function f2(x, ring)
    particles = [0.001 0.0001 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0]
    beam = Beam(particles)
    changed_idx = [3]
    # new_D1 = DRIFT(len=x)
    lenQ1 = ring[3].len
    new_Q1 = KQUAD(len=lenQ1, k1=x)
    changed_ele = [new_Q1]
    ADlinepass!(ring, beam, changed_idx, changed_ele)
    return beam.r[1, 1]
end
# println((f1(0.1001, esr)-f1(0.1,  esr))/0.0001)
# println((f2(-0.22+ 1e-15, esr)-f2(-0.22,  esr))/1e-15)

# dx = [1.0, 1.0]
# @time grad1 = autodiff(Forward, f1, Duplicated, Duplicated(0.1, 1.0),  Const(esr))
# @time grad2 = autodiff(Forward, f2, Duplicated, Duplicated(-0.22, 1.0),  Const(esr))
# println(grad1, grad2)

function f(xx, ring)
    x = CTPS(0.0, 1, 6, 3)
    xp = CTPS(0.0, 2, 6, 3)
    y = CTPS(0.0, 3, 6, 3)
    yp = CTPS(0.0, 4, 6, 3)
    delta = CTPS(0.0, 5, 6, 3)
    z = CTPS(0.0, 6, 6, 3)
    rin = [x, xp, y, yp, delta, z]

    changed_idx = [2]
    new_D1 = DRIFT(len=xx)
    changed_ele = [new_D1]

    ADlinepass_TPSA!(ring, rin, changed_idx, changed_ele)
    return rin[1].map[3]
end
# println(f(0.1, esr))
# println(f(0.1001, esr))
# println((f(0.1001, esr)-f(0.1,  esr))/0.0001)
# grad = autodiff(Forward, f, Duplicated, Duplicated(0.1, 1.0),  Const(esr))
# println(grad)
# Load from the file
@time esr = deserialize("esr_main_vector.jls")
twi, pos = twissring(esr, 0.0, 3)

twi_matrix = zeros(length(esr), 7)
for i in eachindex(esr)
    twi_matrix[i, 1] = pos[i]
    twi_matrix[i, 2] = twi[i].betax
    twi_matrix[i, 3] = twi[i].betay
    twi_matrix[i, 4] = twi[i].alphax
    twi_matrix[i, 5] = twi[i].alphay
    twi_matrix[i, 6] = twi[i].dx
    twi_matrix[i, 7] = twi[i].dy
end
# save the matrix as a text file
using DelimitedFiles
writedlm("twiss_matrix.txt", twi_matrix)