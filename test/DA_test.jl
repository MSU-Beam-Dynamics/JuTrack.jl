include("../src/JuTrack.jl")
include("ssrf_ring.jl")
using. JuTrack
using Enzyme
using BenchmarkTools
using Plots
gr()

Enzyme.API.runtimeActivity!(true)

# function dynamic_aperture(x)
#     SSRF = ssrf(x[1])
#     nangle = 18
#     nsteps = 50
#     nturn = 1000

#     angle = LinRange(0, pi, nangle)
#     xl = LinRange(0.001, 0.050, nsteps)
#     DA_bound = zeros(nangle)
#     rin = zeros(nangle,6)
#     for i in 1:nsteps
#         for j in 1:nangle
#             rin[j,1] = xl[i] * cos(angle[j])
#             rin[j,3] = xl[i] * sin(angle[j])
#         end
#         beam = Beam(rin)
#         ringpass!(SSRF, beam, nturn)
#         lost_flag = beam.lost_flag
#         for j in 1:nangle
#             if lost_flag[j] == 1 && DA_bound[j] == 0.0
#                 DA_bound[j] = xl[i]
#             end
#         end
#         println("step no. ", i, " is done")
#     end
#     return DA_bound[9]
# end
# x = [-1.063770]
# # the function of calculating DA is not smooth and continuous, so it is not differentiable
# grad = gradient(Forward, dynamic_aperture, x)
# println(grad)

# DA = dynamic_aperture(x)
# println(DA)
# DA_bound = zeros(18, 2)
# angle = LinRange(0, pi, 18)
# for i in 1:18
#     DA_bound[i,1] = DA[i] * cos(angle[i])
#     DA_bound[i,2] = DA[i] * sin(angle[i])
# end
# plot(DA_bound[:,1], DA_bound[:,2], seriestype = :scatter, label = "DA bound", 
#     xlabel = "x(m)", ylabel = "y(m)", title = "Dynamic Aperture", legend = :topleft, 
#     xlims = (-0.035, 0.035), ylims = (-0.001, 0.035))
# # save
# savefig("DA_bound.png")


function stability_score(amplitude, k)
    # k = 0.6 
    if amplitude < 0
        amplitude = -amplitude
    end
    return 1 / (1 + exp(-10 * (amplitude - k)))
end

function calculate_DA(x)
    SSRF = ssrf(x[1])
    nangle = 18
    nsteps = 50
    nturn = 1000

    angle = LinRange(0, pi, nangle)
    xl = LinRange(0.001, 0.050, nsteps)
    rin = zeros(nsteps,6)
    DA_bound = zeros(nangle)

    for i in 1:nangle
        for j in 1:nsteps
            rin[j,1] = xl[j] * cos(angle[i])
            rin[j,3] = xl[j] * sin(angle[i])
        end
        beam = Beam(rin)
        ringpass!(SSRF, beam, nturn)
        scores = [stability_score(beam.r[i,1], 0.6) for i in 1:nsteps]
        cumulative_score = cumsum(scores)
        # println(cumulative_score, "at angle ", i)
        for j in 1:nsteps
            if cumulative_score[j] >= 0.5
                DA_bound[i] = cumulative_score[j] * xl[j]
                break
            end
        end
        # println(DA_bound[i], "at angle ", i)

        println("step no. ", i, " is done")
    end

    return DA_bound
end
# DA = calculate_DA([-1.063770])
# println(DA)
grad = gradient(Forward, calculate_DA, [-1.063770])
println(grad)
# ([-1.7458596704903375e6, -2.929149727424697e27, -9.87161306985279e10, -1.0855096299906201e27, -4.05758220223344e19, -5.872441061104684e9, 1.5475047969855115e12, -1.8899567826262295e10, 1.5212146823181528e14, -2.925262176288771e30, -2.7770622661384356e16, 1.5155426669651458e36, 1.3268510154830861e8, 2.756127521817982e12, 2.417128170029123e12, 626695.8142906541, -7.350280751321233e17, -8.431664871189388e27],)
# the results are not correct, the function is not differentiable