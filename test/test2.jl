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

function linepass2!(line, line2, particles::Beam, particles2::Beam)
    # Note!!! A lost particle's coordinate will not be marked as NaN or Inf like other softwares 
    # Check if the particle is lost by checking the lost_flag
    np = particles.nmacro
    particles6 = matrix_to_array(particles.r)
    particles6_2 = matrix_to_array(particles2.r)

    for i in eachindex(line)
        pass!(line[i], particles6, np, particles)  
        pass!(line2[i], particles6_2, np, particles2)     
        if abs(particles6[1] - particles6_2[1]) > 1e-4 
            println(i)
            println(particles6[1]," ", particles6_2[1])
            println(particles6[2], " ", particles6_2[2]) 
        end
    end
    rout = array_to_matrix(particles6, np)
    particles.r = rout
    return nothing
end

function linepass3!(line, particles::Beam)
    # Note!!! A lost particle's coordinate will not be marked as NaN or Inf like other softwares 
    # Check if the particle is lost by checking the lost_flag
    np = particles.nmacro
    particles6 = matrix_to_array(particles.r)
    x = zeros(length(line))
    px = zeros(length(line))
    y = zeros(length(line))
    py = zeros(length(line))
    z = zeros(length(line))
    delta = zeros(length(line))
    s = zeros(length(line))
    for i in eachindex(line)
        pass!(line[i], particles6, np, particles)  
        x[i] = particles6[1]
        px[i] = particles6[2]
        y[i] = particles6[3]
        py[i] = particles6[4]
        z[i] = particles6[5]
        delta[i] = particles6[6]
        s[i] = total_length(line[1:i])
    end
    # for i in eachindex(line)
    #     pass!(line[i], particles6, np, particles)  
    #     x[i+length(line)] = particles6[1]
    #     px[i+length(line)] = particles6[2]
    #     s[i+length(line)] = total_length(line[1:i])
    # end
    rout = array_to_matrix(particles6, np)
    particles.r = rout
    return x, px, y, py, z, delta, s
end
# particle = [0.0005 0.0 0.0001 0.0 0.0 0.0]

# L0 = 3833.99838673867
# C0 = 299792458.0
# T0 = L0 / C0
# beam = Beam(particle, energy=17.846262619763e9, T0=T0)

# esr = deserialize("test/esr_main.jls")
# esr_rad = deserialize("test/esr_main_rad.jls")
# rout = ringpass!(esr_rad, beam, 100, true)
# # println(rout[1, :])
# # println(use_exact_Hamiltonian)
# # use_exact_drift(1)
# # beam = Beam(particle, energy=18e9)
# # rout = ringpass!(esr, beam, 100, true)
# # println(rout[1, :])
# # println(use_exact_Hamiltonian)
# using CSV, DataFrames
# path = "C:/Users/WAN/Downloads/ESR/Version-6.2/tracking_results_1turns_rad.tfs"
# data = CSV.File(path, delim=' ', ignorerepeated=true, header=54, comment="#") |> DataFrame
# MAD_result = zeros(838, 7)
# for i in 1:838
#     MAD_result[i, 1] = data[i, 3]
#     MAD_result[i, 2] = data[i, 4]
#     MAD_result[i, 3] = data[i, 5]
#     MAD_result[i, 4] = data[i, 6]
#     MAD_result[i, 5] = data[i, 7]
#     MAD_result[i, 6] = data[i, 8]
#     MAD_result[i, 7] = data[i, 9]
# end

# using Plots
# rout_m = zeros(100, 6)
# for i in 1:100
#     rout_m[i, 1] = rout[i][1,1]
#     rout_m[i, 2] = rout[i][1,2]
#     rout_m[i, 3] = rout[i][1,3]
#     rout_m[i, 4] = rout[i][1,4]
#     rout_m[i, 5] = rout[i][1,5]
#     rout_m[i, 6] = rout[i][1,6]
# end
# p1 = plot(1:100, MAD_result[1:100, 1], label="MADX x", xlabel="Turn", ylabel="Position (m)", title="Tracking Results", legend=:topleft)
# plot!(1:100, rout_m[:, 1], label="JuTrack x")
# p2 = plot(1:100, MAD_result[1:100, 2], label="MADX xp", xlabel="Turn", ylabel="Position (m)", title="Tracking Results", legend=:topleft)
# plot!(1:100, rout_m[:, 2], label="JuTrack xp")
# p3 = plot(1:100, MAD_result[1:100, 3], label="MADX y", xlabel="Turn", ylabel="Position (m)", title="Tracking Results", legend=:topleft)
# plot!(1:100, rout_m[:, 3], label="JuTrack y")
# p4 = plot(1:100, MAD_result[1:100, 4], label="MADX yp", xlabel="Turn", ylabel="Position (m)", title="Tracking Results", legend=:topleft)
# plot!(1:100, rout_m[:, 4], label="JuTrack yp")
# p_all = plot(p1, p2, p3, p4, layout=(2, 2), size=(1000, 600))
# savefig(p_all, "tracking_results.png")
# open("output_file.txt", "w") do file
#     for i = 1:5550
#         write(file, esr[i].name * "\n")
#     end
# end
# p5 = plot(MAD_result[:, 1], MAD_result[:, 2], label="MADX", xlabel="x", ylabel="xp", title="H Phase Space", legend=:topleft)
# plot!(rout_m[:, 1], rout_m[:, 2], label="JuTrack")
# p6 = plot(MAD_result[:, 3], MAD_result[:, 4], label="MADX", xlabel="y", ylabel="yp", title="V Phase Space", legend=:topleft)
# plot!(rout_m[:, 3], rout_m[:, 4], label="JuTrack")
# p_all2 = plot(p5, p6, layout=(1, 2), size=(800, 360))

