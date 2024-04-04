include("../src/JuTrack.jl")
using .JuTrack
include("ssrf_ring.jl")
using Plots  
using Enzyme
using ProgressMeter

RING = ssrf(-1.063770, 0)
xlist = randn(1000) * 0.003 
ylist = randn(1000) * 0.003 
# particle = zeros(1000, 6)
# particle[:, 1] = xlist
# particle[:, 3] = ylist

# xlist = range(-0.005, stop=0.005, length=1000)


function final_x(x)
    global RING
    beam = Beam([x[1] x[2] x[3] x[4] 0.0 0.0], energy=3.5e9)
    ringpass!(RING, beam, 1000)
    return beam.r
end

grad = []
p = Progress(1000, 1)
for i in 1:1000
    next!(p)
    g = jacobian(Forward, final_x, [xlist[i], 0.0, ylist[i], 0.0], Val(4))
    push!(grad, g)
    sleep(0.1)
end

grad_m = zeros(1000, 12)
for i in 1:1000
    grad_m[i, :] = grad[i]
end

# plot x-y with grad_m as color
using LaTeXStrings
p1=scatter(xlist, ylist, zcolor=log10.(abs.(grad_m[:,1])), color=:viridis, markerstrokewidth=0, label=L"log(|dx_{final}/dx_{init}|)", xlabel="x(m)", ylabel="y(m)", colorbar=true, markersize=2, colormap=:jet)
p2=scatter(xlist, ylist, zcolor=log10.(abs.(grad_m[:,2])), color=:viridis, markerstrokewidth=0, label=L"log(|dp_x_{final}/dx_{init}|)", xlabel="x(m)", ylabel="y(m)", colorbar=true, markersize=2, colormap=:jet)
p3=scatter(xlist, ylist, zcolor=log10.(abs.(grad_m[:,3])), color=:viridis, markerstrokewidth=0, label=L"log(|dy_{final}/dx_{init}|)", xlabel="x(m)", ylabel="y(m)", colorbar=true, markersize=2, colormap=:jet)
p4=scatter(xlist, ylist, zcolor=log10.(abs.(grad_m[:,4])), color=:viridis, markerstrokewidth=0, label=L"log(|dp_y_{final}/dx_{init}|)", xlabel="x(m)", ylabel="y(m)", colorbar=true, markersize=2, colormap=:jet)
plot(p1, p2, p3, p4, layout=(2,2), size=(850, 600))

p1=scatter(xlist, ylist, zcolor=log10.(abs.(grad_m[:,7])), color=:viridis, markerstrokewidth=0, label=L"log(|dx_{final}/dy_{init}|)", xlabel="x(m)", ylabel="y(m)", colorbar=true, markersize=2, colormap=:jet)
p2=scatter(xlist, ylist, zcolor=log10.(abs.(grad_m[:,8])), color=:viridis, markerstrokewidth=0, label=L"log(|dp_x_{final}/dy_{init}|)", xlabel="x(m)", ylabel="y(m)", colorbar=true, markersize=2, colormap=:jet)
p3=scatter(xlist, ylist, zcolor=log10.(abs.(grad_m[:,9])), color=:viridis, markerstrokewidth=0, label=L"log(|dy_{final}/dy_{init}|)", xlabel="x(m)", ylabel="y(m)", colorbar=true, markersize=2, colormap=:jet)
p4=scatter(xlist, ylist, zcolor=log10.(abs.(grad_m[:,10])), color=:viridis, markerstrokewidth=0, label=L"log(|dp_y_{final}/dy_{init}|)", xlabel="x(m)", ylabel="y(m)", colorbar=true, markersize=2, colormap=:jet)
plot(p1, p2, p3, p4, layout=(2,2), size=(850, 600))


# na = 13
# nl = 20
# angle_list = range(0, stop=pi, length=na)
# xl = range(0.001, stop=0.02, length=nl)
# particle = zeros(na*nl, 6)
# for i in 1:na
#     for j in 1:nl
#         particle[(i-1)*nl+j, 1] = xl[j]*cos(angle_list[i])
#         particle[(i-1)*nl+j, 3] = xl[j]*sin(angle_list[i])
#     end
# end
# plot(particle[:,1], particle[:,3], seriestype=:scatter, label="Initial particle distribution")
# beam = Beam(particle, energy=3.5e9)
# pringpass!(RING, beam, 1000)
# lost_flag = beam.lost_flag
# index = findall(lost_flag .== 0)
# index1 = findall(lost_flag .== 1)
# plot(particle[index,1], particle[index,3], seriestype=:scatter, label="Survived")
# plot!(particle[index1,1], particle[index1,3], seriestype=:scatter, label="Lost", xlabel="x(m)", ylabel="y(m)")