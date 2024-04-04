include("../src/JuTrack.jl")
include("ssrf_ring.jl")
using. JuTrack
using Enzyme
using Serialization
using Plots

ESR = deserialize("test/esr_main_rad.jls")
ESR_nocrab = deserialize("test/esr_main_rad_craboff.jls")
ESR_norad = deserialize("test/esr_main.jls")
E0 = 17.846262619763e9

na = 13
nl = 20
angle_list = range(0, stop=pi, length=na)
xlist = range(0.0001, stop=0.002, length=nl)
ylist = range(0.00005, stop=0.001, length=nl)
particle = zeros(na*nl, 6)
for i in 1:na
    for j in 1:nl
        particle[(i-1)*nl+j, 1] = xlist[j]*cos(angle_list[i])
        particle[(i-1)*nl+j, 3] = ylist[j]*sin(angle_list[i])
    end
end
plot(particle[:,1], particle[:,3], seriestype=:scatter, label="Initial particle distribution")

beam = Beam(particle, energy=E0)
beam1 = Beam(particle, energy=E0)

@time ringpass!(ESR_norad, beam, 1000)
@time pringpass!(ESR_norad, beam1, 1000)

# plot whose lost_flag == 0
lost_flag = beam1.lost_flag
index = findall(lost_flag .== 0)
index1 = findall(lost_flag .== 1)
plot(particle[index,1], particle[index,3], seriestype=:scatter, label="Survived")
plot!(particle[index1,1], particle[index1,3], seriestype=:scatter, label="Lost", xlabel="x(m)", ylabel="y(m)")
savefig("DA_CCoff_radoff.png")
