# include("../src/JuTrack.jl")
# include("ssrf_ring.jl")
# using. JuTrack
# using DelimitedFiles


# RF0 = RFCA(name="RF0", len=4.01667, volt=3.78e6, h=7560.0, freq=591.1397738e6, energy=17.846262619763e9, philag=0.0)
# # particle = [4.698798040560031e-5 -4.880534093529371e-6 -2.2043075293955083e-8 -6.922060539225577e-9 -0.00028773568241259454 -0.0007275489333155954]
# particle = [0.0 0.0 0.0 0.0 0.0001 0.0]
# beam = Beam(particle, energy=17.846262619763e9)
# linepass!([RF0],beam)

# ring = ssrf(-1.063770, 1)
# particle = [0.0005 0.0001 0.0002 0.00005 0.0 0.0]
# beam = Beam(particle, energy=3500.0e6)

# function linepass1!(line, particles::Beam, ref)
#     # Note!!! A lost particle's coordinate will not be marked as NaN or Inf like other softwares 
#     # Check if the particle is lost by checking the lost_flag
#     np = particles.nmacro
#     particles6 = matrix_to_array(particles.r)
#     if length(particles6) != np*6
#         error("The number of particles does not match the length of the particle array")
#     end
#     save_beam = []
#     for i in eachindex(line)
#         # ele = line[i]
#         pass!(line[i], particles6, np, particles)        
#         if isnan(particles6[1]) || isinf(particles6[1])
#             println("The particle is lost at element ", i, "element name is ", line[i].name)
#             rout = array_to_matrix(particles6, np)
#             particles.r = rout
            
#             return nothing
#         end
#         push!(save_beam, copy(particles6))
#     end
#     rout = array_to_matrix(particles6, np)
#     particles.r = rout
#     return save_beam
# end

# rout = ringpass!(ring, beam, 100, true)
# rout_m = zeros(length(rout), 6)
# for i in 1:length(rout)
#     rout_m[i, :] = rout[i]
# end
# data = readdlm("C:/Users/WAN/Desktop/SSRF/temp.txt", ',')

# using Plots
# s = spos(ring)
# plot(1:100, data[5,:], label="AT")
# plot!(1:100, rout_m[:,6], label="jutrack")