# Calculate the gradient of the emittance with respect to the lattice parameters of a FODO cell under transverse space charge effects
using JuTrack
# using PyCall
# np = pyimport("numpy")
# plt = pyimport("matplotlib.pyplot")

function f(D3L)
    distparam = [
        3.677529920673089E-004  ,   # sigx
        8.428925532276500E-004 ,   # sigpx
        -0.828277121044551 ,   # muxpx
        1.0,   # xscale
        1.0,   # pxscale
        0.0,   # xmu1 (mean x)
        0.0,   # xmu2 (mean px)
        3.677529304933903E-004  ,   # sigy
        8.428931246578997E-004  ,   # sigpy
        0.828276927537804 ,   # muypy
        1.0,   # yscale
        1.0,   # pyscale
        0.0,   # xmu3 (mean y)
        0.0,   # xmu4 (mean py)
        1.0,   # sigz
        0.1,   # sigpz
        0.5,   # muzpz
        0.0,   # zscale
        0.0,   # pzscale
        0.0,   # xmu5 (mean z)
        0.0    # xmu6 (mean pz)
    ]
    Npt = 5000
    Pts1 = Gauss3_Dist(distparam, Npt)
    beam = Beam(Pts1, energy=1.0e9, current=200.0, mass=m_p, charge=1.0)

    D1L = 0.2
    D2L = 0.4
    # D3L = 0.2
    Q1L = 0.1
    Q2L = 0.1
    Q1k = 29.6
    Q2k = -29.6
    
    a = 13e-3
    b = 13e-3
    nl = 12
    nm = 12

    D1 = DRIFT_SC(len=D1L, a=a, b=b, Nl=nl, Nm=nm, Nsteps=4)
    D2 = DRIFT_SC(len=D2L, a=a, b=b, Nl=nl, Nm=nm, Nsteps=8)
    D3 = DRIFT_SC(len=D3L, a=a, b=b, Nl=nl, Nm=nm, Nsteps=4)
    Q1 = KQUAD_SC(len=Q1L, k1=Q1k, NumIntSteps=20, a=a, b=b, Nl=nl, Nm=nm, Nsteps=4)
    Q2 = KQUAD_SC(len=Q2L, k1=Q2k, NumIntSteps=20, a=a, b=b, Nl=nl, Nm=nm, Nsteps=4)

    line_SC = [D1, Q1, D2, Q2, D3]
    N = 1
    new_emit= zeros(N+1, 3)
    for i in 1:N
        # println("Turn: ", i)
        if i == 1
            get_emittance!(beam)
            new_emit[i, :] = beam.emittance
        end
        linepass!(line_SC, beam)
        get_emittance!(beam)
        new_emit[i+1, :] = beam.emittance
    end
    return new_emit[end,:] ./ [new_emit[1,1], new_emit[1,2], 1.0]
end

# grad_Q1k = autodiff(Forward, f, Duplicated, Duplicated(29.6, 1.0))
# grad_Q1k_fd = (f(29.6+1e-6) .- f(29.6-1e-6)) / 2e-6
# grad_Q1L = autodiff(Forward, f, Duplicated, Duplicated(0.1, 1.0))
# grad_Q1L_fd = (f(0.1+1e-6) .- f(0.1-1e-6)) / 2e-6
# grad_Q2k = autodiff(Forward, f, Duplicated, Duplicated(-29.6, 1.0))
# grad_Q2k_fd = (f(-29.6+1e-6) .- f(-29.6-1e-6)) / 2e-6
# grad_Q2L = autodiff(Forward, f, Duplicated, Duplicated(0.1, 1.0))
# grad_Q2L_fd = (f(0.1+1e-6) .- f(0.1-1e-6)) / 2e-6
# grad_D1L = autodiff(Forward, f, Duplicated, Duplicated(0.2, 1.0))
# grad_D1L_fd = (f(0.2+1e-6) .- f(0.2-1e-6)) / 2e-6
# grad_D2L = autodiff(Forward, f, Duplicated, Duplicated(0.4, 1.0))
# grad_D2L_fd = (f(0.4+1e-6) .- f(0.4-1e-6)) / 2e-6
grad_D3L = autodiff(Forward, f, Duplicated(0.2, 1.0))
grad_D3L_fd = (f(0.2+1e-6) .- f(0.2-1e-6)) / 2e-6


# final = zeros(7, 4)
# final[1, 1:2] = grad_D1L[1][1:2]
# final[2, 1:2] = grad_Q1L[1][1:2]
# final[3, 1:2] = grad_Q1k[1][1:2]
# final[4, 1:2] = grad_D2L[1][1:2]
# final[5, 1:2] = grad_Q2L[1][1:2]
# final[6, 1:2] = grad_Q2k[1][1:2]
# final[7, 1:2] = grad_D3L[1][1:2]
# final[1, 3:4] = grad_D1L_fd[1:2]
# final[2, 3:4] = grad_Q1L_fd[1:2]
# final[3, 3:4] = grad_Q1k_fd[1:2]
# final[4, 3:4] = grad_D2L_fd[1:2]
# final[5, 3:4] = grad_Q2L_fd[1:2]
# final[6, 3:4] = grad_Q2k_fd[1:2]
# final[7, 3:4] = grad_D3L_fd[1:2]

# using DelimitedFiles
# writedlm("final.txt", final)

# plt.figure(figsize=(7, 4))
# plt.plot(np.arange(1, 8), final[:, 3], "+", label="FD, X", markersize=8)
# plt.plot(np.arange(1, 8), final[:, 1], "x", label="AD, X", markersize=8)
# plt.plot(np.arange(1, 8), final[:, 4], "*", label="FD, Y", markersize=8)
# plt.plot(np.arange(1, 8), final[:, 2], "s", label="AD, Y", markersize=8, markerfacecolor="none")
# plt.xticks(np.arange(1, 8), ["D1L", "Q1L", "Q1k", "D2L", "Q2L", "Q2k", "D3L"], fontsize=14, fontname="Times New Roman")
# plt.yticks(np.arange(-0.5, 0.6, 0.1), fontsize=14, fontname="Times New Roman")
# plt.legend(prop=Dict("family"=>"Times New Roman"), frameon=false)
# plt.ylabel("Gradient", fontsize=16, fontname="Times New Roman")
# plt.ylim(-0.5, 0.5)
# plt.show()
