include("../src/JuTrack.jl")
include("ssrf_ring.jl")
using. JuTrack
using Enzyme
# using BenchmarkTools
Enzyme.API.runtimeActivity!(true)

function test_track(xx)
    # we don't suggest to create lattice in the function, it's better to create lattice outside the function
    particles = [0.001 0.0001 0.0005 0.0002 0.0 0.0; 0.001 0.0 0.0 0.0 0.0 0.0]
    beam = Beam(particles)
    line = ssrf(xx[1])
    linepass!(line, beam)
    return beam.r[1,1]
end
# x = [-1.063770]
# println(test_track(x))
# grad = gradient(Forward, test_track, x, Val(1))
# println(grad)
# @time grad = gradient(Forward, test_track, x, Val(1))


function f_TPS(xx)
    # we don't suggest to create lattice in the function, it's better to create lattice outside the function
    SSRF = ssrf(xx[1])
    x = CTPS(0.0, 1, 6, 3)
    xp = CTPS(0.0, 2, 6, 3)
    y = CTPS(0.0, 3, 6, 3)
    yp = CTPS(0.0, 4, 6, 3)
    delta = CTPS(0.0, 5, 6, 3)
    z = CTPS(0.0, 6, 6, 3)
    rin = [x, xp, y, yp, delta, z]
    linepass_TPSA!(SSRF, rin)
    return rin[1].map[2]
end
# println(f_TPS([-1.063770]))
# x = [-1.063770]
# println(f_TPS(x))
# @btime f_TPS(x) 
# @btime grad = gradient(Forward, f_TPS, x)


function twiss_test(xx)
    # we don't suggest to create lattice in the function, it's better to create lattice outside the function
    SSRF = ssrf(xx[1])
    twiss_in = EdwardsTengTwiss(betax=1.0,betay=2.0)
    ss, name, twiss_out = Twissline(twiss_in, SSRF, 0.0, 1, length(SSRF))
    return twiss_out.betax
end
println(twiss_test([-1.063770]))
# grad = gradient(Forward, twiss_test, [-1.063770])
# println(grad)
# @btime begin 
#     twiss_test([-1.063770])
# end
# @btime begin 
#     grad = gradient(Forward, twiss_test, [-1.063770])
# end

function twiss_test_refpts(xx)
    # we don't suggest to create lattice in the function, it's better to create lattice outside the function
    SSRF = ssrf(xx[1])
    twiss_in = EdwardsTengTwiss(betax=1.0,betay=2.0)
    idx = collect(1:length(SSRF))
    ss, name, twiss_out = Twissline(twiss_in, SSRF, 0.0, 1, idx)
    return ss, name, twiss_out
end
# ss, name, twiss_out = twiss_test_refpts([-1.063770])

# beta = zeros(length(twiss_out), 2)
# for i in eachindex(twiss_out)
#     beta[i,1] = twiss_out[i].betax
#     beta[i,2] = twiss_out[i].betay
# end

# println("beta_x: ", beta[:,1])
# println("beta_y: ", beta[:,2])

# run the plot code in REPL
# using Plots
# gr()
# p = plot(ss, beta[:,1], label="betax", xlabel="s(m)", ylabel="beta(m)")
# plot!(ss, beta[:,2], label="betay")
# display(p)
# using DelimitedFiles
# writedlm("ss.txt", ss)
# writedlm("beta.txt", beta)
