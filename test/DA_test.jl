include("../src/JuTrack.jl")
include("ssrf_ring.jl")
using. JuTrack
using Enzyme
using BenchmarkTools
using Plots
gr()

Enzyme.API.runtimeActivity!(true)

function dynamic_aperture(x)
    SSRF = ssrf(x[1])
    nangle = 18
    nsteps = 50
    nturn = 1000

    angle = LinRange(0, pi, nangle)
    xl = LinRange(0.001, 0.050, nsteps)
    DA_bound = zeros(nangle)
    rin = zeros(nangle,6)
    for i in 1:nsteps
        for j in 1:nangle
            rin[j,1] = xl[i] * cos(angle[j])
            rin[j,3] = xl[i] * sin(angle[j])
        end
        beam = Beam(rin)
        ringpass!(SSRF, beam, nturn)
        lost_flag = beam.lost_flag
        for j in 1:nangle
            if lost_flag[j] == 1 && DA_bound[j] == 0.0
                DA_bound[j] = xl[i]
            end
        end
        println("step no. ", i, " is done")
    end
    return DA_bound[9]
end
x = [-1.063770]
# the function of calculating DA is not smooth and continuous, so it is not differentiable
grad = gradient(Forward, dynamic_aperture, x)
println(grad)

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


