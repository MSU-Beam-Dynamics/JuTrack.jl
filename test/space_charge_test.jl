# include("../src/JuTrack.jl")
using JuTrack
using Distributions, Plots
using ProgressMeter
using DelimitedFiles


function calculate_twiss_parameters(x, px)
    # Number of particles
    N = length(x)
    
    # Mean values
    mean_x = mean(x)
    mean_px = mean(px)
    
    # Second moments
    mean_x2 = mean(x.^2)
    mean_px2 = mean(px.^2)
    mean_xpx = mean(x .* px)
    
    # Covariances
    cov_x2 = mean_x2 - mean_x^2
    cov_px2 = mean_px2 - mean_px^2
    cov_xpx = mean_xpx - mean_x * mean_px
    
    # Emittance
    epsilon = sqrt(cov_x2 * cov_px2 - cov_xpx^2)
    
    # Twiss parameters
    beta = cov_x2 / epsilon
    gamma = cov_px2 / epsilon
    alpha = -cov_xpx / epsilon
    
    return alpha, beta, gamma, epsilon
end
# alpha, beta, gamma, epsilon = calculate_twiss_parameters(particles[:, 1], particles[:, 2])

# σx = 3.677529920673089e-4
# σpx = 8.428925532276500e-4
# μxpx = -0.828277121044551 
# σy = 3.677529304933903e-4
# σpy = 8.428931246578997e-4
# μypy = 0.828276927537804 
# Σ = [
#     σx^2  μxpx*σx*σpx  0             0;
#     μxpx*σx*σpx  σpx^2  0             0;
#     0       0      σy^2  μypy*σy*σpy;
#     0       0      μypy*σy*σpy  σpy^2
# ]

# # Define the multivariate normal distribution
# dist = MvNormal(zeros(4), Σ)  # Using zeros(4) as there are no explicit means given

# # Generate samples
# num_samples = 5000
# samples = rand(dist, num_samples)
# # genrate 6D beam
# particles = [samples[1, :] samples[2, :] samples[3, :] samples[4, :] zeros(num_samples) zeros(num_samples)]
# beam = Beam(particles, energy=1.0e9, current=200.0, mass=m_p, charge=1.0)


function f(D3L)
    beam = Beam(zeros(5000,6), energy=1.0e9, current=200.0, mass=m_p, charge=1.0, 
        emittance=[1.7309472352590595e-07, 1.7346194711492812e-07, 0.0])
    
    vbase=3.42*8.5e6
    ϕs=10.0
    vact=vbase/cos(ϕs*π/180.0)
    freq = 591e6
    mainRFe=AccelCavity(freq=freq, volt=vact, h=7560.0, phis=π-ϕs*π/180.0)
    tunex, tuney=50.08, 44.14
    αc=3.42/tunex/tunex
    lmap=LongitudinalRFMap(αc, mainRFe)
    opt=optics4DUC(0.7772150933108849, -1.4747086861162038, 0.7800845105254005, 1.4768023077802253)
    initilize_6DGaussiandist!(beam, opt, lmap)
    beam.r[:, 5] .= 0.0
    beam.r[:, 6] .= 0.0
    beam.r[1, :] .= 0.0

    D1L = 0.2
    D2L = 0.4
    # D3L = 0.2
    Q1L = 0.1
    Q2L = 0.1
    Q1k = 29.6
    Q2k = -29.6
    
    D1 = DRIFT(len=D1L/4)
    D2 = DRIFT(len=D2L/8)
    D3 = DRIFT(len=D3L/4)
    Q1 = QUAD(len=Q1L/4, k1=Q1k)
    Q2 = QUAD(len=Q2L/4, k1=Q2k)
    a = 13e-3
    b = 13e-3
    nl = 12
    nm = 12
    SC_D1 = SPACECHARGE(effective_len=0.2/4, a=a, b=b, Nl=nl, Nm=nm)
    SC_D2 = SPACECHARGE(effective_len=0.4/8, a=a, b=b, Nl=nl, Nm=nm)
    SC_D3 = SPACECHARGE(effective_len=0.2/4, a=a, b=b, Nl=nl, Nm=nm)
    SC_Q1 = SPACECHARGE(effective_len=0.1/4, a=a, b=b, Nl=nl, Nm=nm)
    SC_Q2 = SPACECHARGE(effective_len=0.1/4, a=a, b=b, Nl=nl, Nm=nm)
    # line = [D1,D1,D1,D1, Q1,Q1,Q1,Q1, D2,D2,D2,D2,D2,D2,D2,D2, Q2,Q2,Q2,Q2, D3,D3,D3,D3]
    line_SC = [D1,SC_D1,D1,SC_D1,D1,SC_D1,D1,SC_D1, 
                Q1,SC_Q1,Q1,SC_Q1,Q1,SC_Q1,Q1,SC_Q1, 
                D2,SC_D2,D2,SC_D2,D2,SC_D2,D2,SC_D2, D2,SC_D2,D2,SC_D2,D2,SC_D2,D2,SC_D2,
                Q2,SC_Q2,Q2,SC_Q2,Q2,SC_Q2,Q2,SC_Q2, 
                D3,SC_D3,D3,SC_D3,D3,SC_D3,D3,SC_D3]
    
    N = 8
    new_emit= zeros(N+1, 3)
    # new_emit1= zeros(N+1, 3)
    # beam1 = Beam(beam)
    for i in 1:N
        # println("Turn: ", i)
        if i == 1
            get_emittance!(beam)
            new_emit[i, :] = beam.emittance
            # new_emit1[i, :] = beam1.emittance
        end
        linepass!(line_SC, beam)
        get_emittance!(beam)
        new_emit[i+1, :] = beam.emittance
        # get_emittance!(beam1)
        # new_emit1[i+1, :] = beam1.emittance
    end
    return new_emit[end,:] ./ [1.7309472352590595e-07, 1.7346194711492812e-07, 1.0]
end

grad_Q1k = autodiff(Forward, f, Duplicated, Duplicated(29.6, 1.0))
grad_Q1k_fd = (f(29.6+1e-6) .- f(29.6-1e-6)) / 2e-6
grad_Q1L = autodiff(Forward, f, Duplicated, Duplicated(0.1, 1.0))
grad_Q1L_fd = (f(0.1+1e-6) .- f(0.1-1e-6)) / 2e-6
grad_Q2k = autodiff(Forward, f, Duplicated, Duplicated(-29.6, 1.0))
grad_Q2k_fd = (f(-29.6+1e-6) .- f(-29.6-1e-6)) / 2e-6
grad_Q2L = autodiff(Forward, f, Duplicated, Duplicated(0.1, 1.0))
grad_Q2L_fd = (f(0.1+1e-6) .- f(0.1-1e-6)) / 2e-6
grad_D1L = autodiff(Forward, f, Duplicated, Duplicated(0.2, 1.0))
grad_D1L_fd = (f(0.2+1e-6) .- f(0.2-1e-6)) / 2e-6
grad_D2L = autodiff(Forward, f, Duplicated, Duplicated(0.4, 1.0))
grad_D2L_fd = (f(0.4+1e-6) .- f(0.4-1e-6)) / 2e-6
grad_D3L = autodiff(Forward, f, Duplicated, Duplicated(0.2, 1.0))
grad_D3L_fd = (f(0.2+1e-6) .- f(0.2-1e-6)) / 2e-6

final = zeros(7, 4)
final[1, 1:2] = grad_D1L[2][1:2]
final[2, 1:2] = grad_Q1L[2][1:2]
final[3, 1:2] = grad_Q1k[2][1:2]
final[4, 1:2] = grad_D2L[2][1:2]
final[5, 1:2] = grad_Q2L[2][1:2]
final[6, 1:2] = grad_Q2k[2][1:2]
final[7, 1:2] = grad_D3L[2][1:2]
final[1, 3:4] = grad_D1L_fd[1:2]
final[2, 3:4] = grad_Q1L_fd[1:2]
final[3, 3:4] = grad_Q1k_fd[1:2]
final[4, 3:4] = grad_D2L_fd[1:2]
final[5, 3:4] = grad_Q2L_fd[1:2]
final[6, 3:4] = grad_Q2k_fd[1:2]
final[7, 3:4] = grad_D3L_fd[1:2]

using DelimitedFiles
writedlm("final.txt", final)

beam = Beam(zeros(5000,6), energy=1.0e9, current=200.0, mass=m_p, charge=1.0, 
emittance=[1.7309472352590595e-07, 1.7346194711492812e-07, 0.0])

vbase=3.42*8.5e6
ϕs=10.0
vact=vbase/cos(ϕs*π/180.0)
freq = 591e6
mainRFe=AccelCavity(freq=freq, volt=vact, h=7560.0, phis=π-ϕs*π/180.0)
tunex, tuney=50.08, 44.14
αc=3.42/tunex/tunex
lmap=LongitudinalRFMap(αc, mainRFe)
opt=optics4DUC(0.7772150933108849, -1.4747086861162038, 0.7800845105254005, 1.4768023077802253)
initilize_6DGaussiandist!(beam, opt, lmap)
beam1 = Beam(beam)
# ini = copy(beam.r)

# N = 800
# new_emit = zeros(N+1, 3)
# new_emit1 = zeros(N+1, 3)
# NLOST = zeros(N)

# println("Start tracking")
# prog = Progress(N)
# for i in 1:N
#     # println("Turn: ", i)
#     if i == 1
#         get_emittance!(beam)
#         get_emittance!(beam1)
#         new_emit[i, :] = beam.emittance
#         new_emit1[i, :] = beam1.emittance
#     end
#     linepass!(line, beam)
#     linepass!(line_SC, beam1)
#     NLOST[i] = sum(beam1.lost_flag)
#     get_emittance!(beam)
#     get_emittance!(beam1)
#     new_emit[i+1, :] = beam.emittance
#     new_emit1[i+1, :] = beam1.emittance
#     next!(prog)
# end

# using DelimitedFiles
# # writedlm("emit.txt", new_emit)
# new_emit = readdlm("emit.txt")

# p1=plot(0:200, new_emit[1:100:end, 1].*1e6, title = "x emit", xlabel = "turns*100", ylabel = "emit (mm*mrad)", label = "without SC")
# plot!(0:200, new_emit1[1:100:end, 1].*1e6, title = "x emit", xlabel = "turns*100", ylabel = "emit (mm*mrad)",  label = "with SC")
# p2=plot(0:200, new_emit1[1:100:end, 2].*1e6, title = "y emit", xlabel = "turns*100", ylabel = "emit (mm*mrad)", label = "with SC")
# plot!(0:200, new_emit[1:100:end, 2].*1e6, title = "y emit", xlabel = "turns*100", ylabel = "emit (mm*mrad)", label = "without SC")
# plot(p1, p2, layout = (1, 2))

# using PyCall
# np = pyimport("numpy")
# plt = pyimport("matplotlib.pyplot")
# N=25000
# plt.figure(figsize=(9, 4))
# plt.subplot(1, 2, 1)
# plt.plot(np.arange(N+1), new_emit[1:N+1, 1].*1e6, label="without SC")
# plt.plot(np.arange(N+1), new_emit1[1:N+1, 1].*1e6, label="with SC")
# plt.xlabel("periods", fontsize=16, fontname="Times New Roman")
# plt.ylabel("emit (mm*mrad)", fontsize=16, fontname="Times New Roman")
# plt.title("x emittance", fontsize=16, fontname="Times New Roman")
# plt.xticks(fontsize=14, fontname="Times New Roman")
# plt.yticks(fontsize=14, fontname="Times New Roman")
# # plt.yscale("log")
# plt.legend()
# plt.subplot(1, 2, 2)
# plt.plot(np.arange(N+1), new_emit[1:N+1, 2].*1e6, label="without SC")
# plt.plot(np.arange(N+1), new_emit1[1:N+1, 2].*1e6, label="with SC")
# plt.xlabel("periods", fontsize=16, fontname="Times New Roman")
# plt.ylabel("emit (mm*mrad)", fontsize=16, fontname="Times New Roman")
# plt.title("y emittance", fontsize=16, fontname="Times New Roman")
# plt.legend(prop=Dict("family"=>"Times New Roman"))
# plt.xticks(fontsize=14, fontname="Times New Roman")
# plt.yticks(fontsize=14, fontname="Times New Roman")
# # plt.yscale("log")
# plt.tight_layout()
# plt.show()

# emit_growth = zeros(N)
# for i in 1:N
#     emit_growth[i] = (new_emit1[i+1, 1] / new_emit[1, 1] * new_emit1[i+1, 2] / new_emit[1, 2] - 1) * 100
# end
# plt.plot(np.arange(N), emit_growth, label="emit growth")
# plt.xlabel("periods", fontsize=16, fontname="Times New Roman")
# plt.ylabel("emit growth (%)", fontsize=16, fontname="Times New Roman")
# plt.xticks(fontsize=14, fontname="Times New Roman")
# plt.yticks(fontsize=14, fontname="Times New Roman")
# plt.show()

# idx_sur = findall(x -> x == 0, beam1.lost_flag)
# plt.figure(figsize=(11, 8))
# plt.subplot(2, 2, 1)
# plt.scatter(beam.r[:, 1], beam.r[:, 2], s=0.1)
# plt.xlabel("x (m)", fontsize=16, fontname="Times New Roman")
# plt.ylabel("px", fontsize=16, fontname="Times New Roman")
# plt.title("Without space charge", fontsize=16, fontname="Times New Roman")
# plt.xlim(-0.008, 0.008)
# plt.ylim(-0.017, 0.017)
# plt.subplot(2, 2, 2)
# plt.scatter(beam1.r[idx_sur, 1], beam1.r[idx_sur, 2], s=0.1)
# plt.xlabel("x (m)", fontsize=16, fontname="Times New Roman")
# plt.ylabel("px", fontsize=16, fontname="Times New Roman")
# plt.title("With space charge", fontsize=16, fontname="Times New Roman")
# plt.xlim(-0.008, 0.008)
# plt.ylim(-0.017, 0.017)
# plt.subplot(2, 2, 3)
# plt.scatter(beam.r[:, 3], beam.r[:, 4], s=0.1)
# plt.xlabel("y (m)", fontsize=16, fontname="Times New Roman")
# plt.ylabel("py", fontsize=16, fontname="Times New Roman")
# plt.title("Without space charge", fontsize=16, fontname="Times New Roman")
# plt.xlim(-0.008, 0.008)
# plt.ylim(-0.017, 0.017)
# plt.subplot(2, 2, 4)
# plt.scatter(beam1.r[idx_sur, 3], beam1.r[idx_sur, 4], s=0.1)
# plt.xlabel("y (m)", fontsize=16, fontname="Times New Roman")
# plt.ylabel("py", fontsize=16, fontname="Times New Roman")
# plt.title("With space charge", fontsize=16, fontname="Times New Roman")
# plt.xlim(-0.008, 0.008)
# plt.ylim(-0.017, 0.017)
# plt.tight_layout()
# plt.show()

