# This script shows how to use the machine learning-based space charge model to track the particles through the beamline, 
# and calculate the emittance growth.

# To run the ML model, you must have PyCall, StatsBase, and ONNXRuntime installed in your Julia environment. 
# The Python path of PyCall should be correctly set. The Python environment should have the following packages: Pytorch, numpy and scipy.

# The machine learning model is over 100MB, so it is not included in the repository. It will be downloaded in this script.

using JuTrack
using Downloads
using ProgressMeter
using DelimitedFiles
import ONNXRunTime as ORT
include("ML_SC.jl")
pushfirst!(PyVector(pyimport("sys")."path"), "src/MLmodel")
GAN_model = pyimport("GAN_model")
url = "https://github.com/MSU-Beam-Dynamics/JuTrack.jl/releases/download/v0.3.0/generator_TVregularization_epoch_5000.onnx"
dest_path = "src/MLmodel/model.onnx"
Downloads.download(url, dest_path)
model = ORT.load_inference(dest_path)


x_mean = 6.102458832397462e-05
x_std = 0.0004449417922594988
y_mean = 0.02980841330718319
y_std = 0.46651152626100306
xedges = range(-0.01, 0.01; length=128+1)  
yedges = range(-0.01, 0.01; length=128+1)
delta = 0.02 / 128
xaxis = np.arange(-0.01+delta/2, 0.01+delta/2, delta)

D1L = 0.2
D2L = 0.4
D3L = 0.2
Q1L = 0.1
Q2L = 0.1
Q1k = 29.6
Q2k = -29.6

a = 13e-3
b = 13e-3
nl = 15
nm = 15

D1 = DRIFT_SC_ML(len=D1L, a=a, b=b, Nl=nl, Nm=nm, Nsteps=4)
D2 = DRIFT_SC_ML(len=D2L, a=a, b=b, Nl=nl, Nm=nm, Nsteps=8)
D3 = DRIFT_SC_ML(len=D3L, a=a, b=b, Nl=nl, Nm=nm, Nsteps=4)
Q1 = KQUAD_SC_ML(len=Q1L, k1=Q1k, NumIntSteps=20, a=a, b=b, Nl=nl, Nm=nm, Nsteps=4)
Q2 = KQUAD_SC_ML(len=Q2L, k1=Q2k, NumIntSteps=20, a=a, b=b, Nl=nl, Nm=nm, Nsteps=4)

line_SC = [D1,Q1,D2,Q2,D3]

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
Npt = 50000
Pts1 = Gauss3_Dist(distparam, Npt, seed=1234)

beam = Beam(Pts1, energy=1.0e9, current=100.0, mass=m_p, charge=1.0)
beam1 = Beam(beam)

N = 10000
new_emit1 = zeros(N+1, 3)
NLOST = zeros(N)

println("Start tracking")
prog = Progress(N)
X1 = zeros(N, 4)
for i in 1:N
    # println("Turn: ", i)
    if i == 1
        get_emittance!(beam1)
        new_emit1[i, :] = beam1.emittance
    end
    linepass_ML!(line_SC, beam1, model, x_mean, x_std, y_mean, y_std, xedges, yedges, xaxis, delta)
    NLOST[i] = sum(beam1.lost_flag)
    get_emittance!(beam1)
    new_emit1[i+1, :] = beam1.emittance
    # save
    if i % 1000 == 0
        writedlm("emit_100A.txt", new_emit1)
        writedlm("lost_100A.txt", NLOST)
    end
    X1[i, :] = beam1.r[1, 1:4]
    next!(prog)
end


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
# plt.legend(prop=Dict("family"=>"Times New Roman"), frameon=false)
# # plt.ylim(0.172, 0.175)
# plt.subplot(1, 2, 2)
# plt.plot(np.arange(N+1), new_emit[1:N+1, 2].*1e6, label="without SC")
# plt.plot(np.arange(N+1), new_emit1[1:N+1, 2].*1e6, label="with SC")
# plt.xlabel("periods", fontsize=16, fontname="Times New Roman")
# plt.ylabel("emit (mm*mrad)", fontsize=16, fontname="Times New Roman")
# plt.title("y emittance", fontsize=16, fontname="Times New Roman")
# plt.legend(prop=Dict("family"=>"Times New Roman"), frameon=false)
# plt.xticks(fontsize=14, fontname="Times New Roman")
# plt.yticks(fontsize=14, fontname="Times New Roman")
# # plt.ylim(0.172, 0.175)
# # plt.yscale("log")
# plt.tight_layout()
# plt.show()

# emit_growth = zeros(N)
# for i in 1:N
#     emit_growth[i] = (new_emit1[i+1, 1] / new_emit1[1, 1] * new_emit1[i+1, 2] / new_emit1[1, 2] - 1) * 100
# end
# plt.plot(np.arange(10000), emit_growth[1:10000], label="emit growth 5000")
# plt.xlabel("periods", fontsize=16, fontname="Times New Roman")
# plt.ylabel("emit growth (%)", fontsize=16, fontname="Times New Roman")
# plt.xticks(fontsize=14, fontname="Times New Roman")
# plt.yticks(fontsize=14, fontname="Times New Roman")
# plt.tight_layout()
# plt.show()