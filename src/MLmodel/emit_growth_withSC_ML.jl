"""
This script demonstrates how to use a machine learning-based space charge model (transverse space charge only) 
to track particles through a beamline and calculate emittance growth.
Related paper: "A symplectic machine learning model for fast simulation of space-charge effects" by J. Wan, Y. Hao and J. Qiang.

The model is based on a neural network trained on a large dataset of particle tracking simulations.
The beam is assumed coasting beam within a 13mm x 13mm perfectly conducting pipe.
The model is trained using PyTorch and saved in ONNX format.
The pre-trained model is over 100MB and is not included in the repository. It will be automatically downloaded in this script.

Essential packagess in Julia environment: 
- JuTrack
- PyCall
- StatsBase
- ONNXRuntime

The Python functionality is provided by the PyCall package, which allows Julia to call Python code and use Python libraries.
The Python path of PyCall should be correctly set up to use the required packages.
The Python environment linked to PyCall must have the following packages: 
- Pytorch
- numpy 
- scipy
"""
using JuTrack
using Downloads
using ProgressMeter
using DelimitedFiles
import ONNXRunTime as ORT

# Julia functions for particle tracking with space charge
# Adjust the path to the file if you are using a different directory. 
# Or download it at https://github.com/MSU-Beam-Dynamics/JuTrack.jl/blob/main/src/MLmodel/ML_SC.jl
include("ML_SC.jl") 

# Load Python functions
# correctly set the path to the python file if you are using a different directory.
# Or download it at https://github.com/MSU-Beam-Dynamics/JuTrack.jl/blob/main/src/MLmodel/GAN_model.py
pushfirst!(PyVector(pyimport("sys")."path"), "src/MLmodel") 
GAN_model = pyimport("GAN_model")

# Download the pre-trained model
url = "https://github.com/MSU-Beam-Dynamics/JuTrack.jl/releases/download/v0.3.0/generator_TVregularization_epoch_5000.onnx"
dest_path = "src/MLmodel/model.onnx"
Downloads.download(url, dest_path)
model = ORT.load_inference(dest_path)

# Parameters for the pre-trained model
x_mean = 6.102458832397462e-05
x_std = 0.0004449417922594988
y_mean = 0.02980841330718319
y_std = 0.46651152626100306
xedges = range(-0.01, 0.01; length=128+1)  
yedges = range(-0.01, 0.01; length=128+1)
delta = 0.02 / 128
xaxis = np.arange(-0.01+delta/2, 0.01+delta/2, delta)

# Beam line
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