# JuTrack use matplotlib in Python to plot the lattice
# PyCall is used to call the Python function from Julia
# matplotlib and numpy are required Python packages in the Python environment associated with PyCall
using Pkg
Pkg.activate("."); Pkg.instantiate() # change "." to your path of JuTrack.jl
using JuTrack
include("../src/utils/lattice_plot.jl")

D1 = DRIFT(len=0.2)
Q1 = QUAD(len=0.1, k1=29.6)
Q2 = QUAD(len=0.3, k1=-29.6)
D2 = DRIFT(len=0.4)
D3 = DRIFT(len=0.2)
S1 = KSEXT(len=0.2, k2=0.0)
O1 = KOCT(len=0.1, k3=0.0)
B1 = SBEND(len=0.6, angle=π/6)
B2 = RBEND(len=0.4, angle=-π/6)
B3 = RBEND(len=0.2, angle=π/6)
B4 = SBEND(len=0.2, angle=-π/6)
RF = RFCA(len=0.2, volt=3.42*8.5e6, freq=591e6)
line = [RF, D1, Q1, D2, Q2, D3, B1, D1, Q1, D2, Q2, D3,
        B2, D1, Q1, D2, Q2, D3, B3, D1, Q1, D2, Q2, D3,
        B4, D1, Q1, D2, Q2, D3, S1, D2, O1, D2, S1, D2]

# Call the plot_lattice function to plot the lattice
plot_lattice(line, 0.3, true)