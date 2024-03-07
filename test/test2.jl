include("../src/JuTrack.jl")
using. JuTrack
using Serialization
using Enzyme
using BenchmarkTools
Enzyme.API.runtimeActivity!(true)


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
    z = CTPS(0.0, 5, 6, 3)
    delta = CTPS(0.0, 6, 6, 3)
    rin = [x, xp, y, yp, z, delta]

    changed_idx = [2]
    new_D1 = DRIFT(len=xx)
    changed_ele = [new_D1]

    ADlinepass_TPSA!(ring, rin, changed_idx, changed_ele)
    return rin[1].map[3]
end

particle = [0.0005 0.0001 0.0 0.0 0.0 0.0; 0.0 0.0 0.0 0.0 0.0 0.0]
beam = Beam(particle, energy=18e9)

esr = deserialize("test/esr_main_vector.jls")
rout = ringpass!(esr, beam, 100, true)

using CSV, DataFrames
path = "C:/Users/WAN/Downloads/ESR/Version-6.2/tracking_results.tfs"
data = CSV.File(path, delim=' ', ignorerepeated=true, header=8, comment="#") |> DataFrame
MAD_result = zeros(100, 6)
for i in 1:100
    MAD_result[i, 1] = data[i+1, 3]
    MAD_result[i, 2] = data[i+1, 4]
    MAD_result[i, 3] = data[i+1, 5]
    MAD_result[i, 4] = data[i+1, 6]
    MAD_result[i, 5] = data[i+1, 7]
    MAD_result[i, 6] = data[i+1, 8]
end

using Plots
rout_m = zeros(100, 6)
for i in 1:100
    rout_m[i, 1] = rout[i][1,1]
    rout_m[i, 2] = rout[i][1,2]
    rout_m[i, 3] = rout[i][1,3]
    rout_m[i, 4] = rout[i][1,4]
    rout_m[i, 5] = rout[i][1,5]
    rout_m[i, 6] = rout[i][1,6]
end
p1 = plot(1:100, MAD_result[:, 1], label="MADX x", xlabel="Turn", ylabel="Position (m)", title="Tracking Results", legend=:topleft)
plot!(1:100, rout_m[:, 1], label="JuTrack x")
p2 = plot(1:100, MAD_result[:, 2], label="MADX xp", xlabel="Turn", ylabel="Position (m)", title="Tracking Results", legend=:topleft)
plot!(1:100, rout_m[:, 2], label="JuTrack xp")
p3 = plot(1:100, MAD_result[:, 3], label="MADX y", xlabel="Turn", ylabel="Position (m)", title="Tracking Results", legend=:topleft)
plot!(1:100, rout_m[:, 3], label="JuTrack y")
p4 = plot(1:100, MAD_result[:, 4], label="MADX yp", xlabel="Turn", ylabel="Position (m)", title="Tracking Results", legend=:topleft)
plot!(1:100, rout_m[:, 4], label="JuTrack yp")
p_all = plot(p1, p2, p3, p4, layout=(2, 2), size=(1000, 600))
# savefig(p_all, "tracking_results.png")

# p5 = plot(MAD_result[:, 1], MAD_result[:, 2], label="MADX", xlabel="x", ylabel="xp", title="H Phase Space", legend=:topleft)
# plot!(rout_m[:, 1], rout_m[:, 2], label="JuTrack")
# p6 = plot(MAD_result[:, 3], MAD_result[:, 4], label="MADX", xlabel="y", ylabel="yp", title="V Phase Space", legend=:topleft)
# plot!(rout_m[:, 3], rout_m[:, 4], label="JuTrack")
# p_all2 = plot(p5, p6, layout=(1, 2), size=(800, 360))