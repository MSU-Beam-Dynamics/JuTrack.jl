# JuTrack

A Julia-based package that enables advanced auto differentiation (AD) for symplectic 6-D particle tracking in particle accelerators.
A mannual can be found at [here](docs/JuTrack_manual.pdf).

# Citation
```
@article{WAN2025109497,
title = {JuTrack: A Julia package for auto-differentiable accelerator modeling and particle tracking},
journal = {Computer Physics Communications},
volume = {309},
pages = {109497},
year = {2025},
issn = {0010-4655},
doi = {https://doi.org/10.1016/j.cpc.2024.109497},
url = {https://www.sciencedirect.com/science/article/pii/S001046552400420X},
author = {Jinyu Wan and Helena Alamprese and Christian Ratcliff and Ji Qiang and Yue Hao}
}
```

# Installation

* Install Julia at [here](https://julialang.org/downloads/oldreleases/) (1.10.4 is required! Use other Julia versions may result in an error).

* Download the package:
```
git clone https://github.com/MSU-Beam-Dynamics/JuTrack.jl
```

* Since Julia and its packages are under very active development, it is strongly recommended to create a separate environment for each piece of work.
To create a JuTrack environment and install its dependencies (errors may occur if not using this method to install all dependencies):
```
cd JuTrack.jl
julia --project=. -e "using Pkg; Pkg.instantiate()"
```

# Use the package in Julia
Activate the JuTrack environment in the Julia code (change the path to your JuTrack directory accordingly):
```
using Pkg
Pkg.activate("/path/to/JuTrack.jl")
Pkg.instantiate()
```

Import JuTrack:
```
using JuTrack
```

# Lattice definition
```
D1 = DRIFT(name="D1", len=1.0)
D2 = DRIFT(name="D2", len=1.0)
D3 = DRIFT(name="D3", len=1.0)
D4 = DRIFT(name="D4", len=1.0)
Q1 = KQUAD(name="Q1", len=1.0, k1=-0.9) 
Q2 = KQUAD(name="Q2", len=1.0, k1=0.3)
B1 = SBEND(name="B1", len=0.6, angle=pi/15.0)
```
Create the lattice as a Julia vector
```
LINE = [D1, Q1, D2, B1, D3, Q2, D4]
```

# Particle tracking
Particles' coordinates are represented as a N * 6 matrix, saved in beam.r
```
particles = rand(10, 6) / 1000
beam = Beam(particles, energy=3.5e9)
```

Tracking
```
linepass!(LINE, beam) # or ringpass!(RING, beam, nturns) for multi-turn tracking
println(beam.r) 
```

# Optics calculation
Obtain periodic Twiss parameters of a ring accelerator
```
twi = periodicEdwardsTengTwiss(RING, 0.0, 0)
```

# Automatic differentiation
We provide two options for fast AD calculation, fastTPSA and Enzyme. fastTPSA is a first-order dual-number module, and Enzyme is a third-party package that runs AD at compiler level. Both methods provide comparably fast AD. The fastTPSA is written in pure Julia that is supposed to be more stable. While Enzyme supports all kinds of differentiable functions, which is not limited by the tracking code in JuTrack.

fastTPSA example:
Obtain derivatives of tracking result w.r.t the quadrupole strength k1
```
set_tps_dim(2) # 2 variables
function tracking_wrt_k1(x1::DTPSAD{NVAR(), Float64}, x2::DTPSAD{NVAR(), Float64})
    D1 = TDRIFT(len=1.0)        # in fastTPSA, we use T+ElementName.  
    D2 = TDRIFT(len=1.0)        # The use of TELEMENT is the same as standard ELEMENT 
    Q1 = TKQUAD(len=1.0)        # All parameters in TELEMENT are described as TPS variables
    Q2 = TKQUAD(len=1.0)

    Q1.k1 = x1
    Q2.k1 = x2

    beam = TBeam([0.1 0.0 0.0 0.0 0.0 0.0])

    LINE = [D1, Q1, D2, Q2] 
    linepass!(LINE, beam)
    return beam.r
end
k1 = -0.9
k2 = 0.3
g, r = Jacobian(tracking_wrt_k1, [k1, k2], true)
```

Enzyme example:
Obtain derivatives of tracking result w.r.t the quadrupole strength k1
```
function tracking_wrt_k1(X)
    D1 = DRIFT(len=1.0)
    D2 = DRIFT(len=1.0)
    Q1 = KQUAD(len=1.0, k1=X[1]) 
    Q2 = KQUAD(len=1.0, k1=X[2])

    beam = Beam([0.1 0.0 0.0 0.0 0.0 0.0])

    # !!! avoid creating large lattice in the differentiable function.
    # !!! load or create your lattice outside the function if it is large.
    LINE = [D1, Q1, D2, Q2] 
    linepass!(LINE, beam)
    return beam.r
end
k1 = -0.9
k2 = 0.3
derivatives, results = jacobian(ForwardWithPrimal, tracking_wrt_k1, [k1, k2])
```

# Parallel computation setting
Multi-threading is available for multi-particle tracking. 
Before using it, please ensure the Julia multi-threading is set up correctly by typing:
```
println("Number of threads in use: ", Threads.nthreads())
```

Parallel computing is available for multi-particle tracking using:
```
plinepass!(beamline, beam)
```
or 
```
pringpass!(beamline, beam, nturns)
```

# Known issues
* JuTrack is actively under development. If you encounter any issues, please open an issue on GitHub or email wan@frib.msu.edu.
* This package currently supports forward AD. Backward AD is still under development.
* Downgrade/upgrade the Julia environment to 1.10.4 for any conflicts.
* Downgrade/upgrade Enzyme version to 0.13.3 with ```Pkg.add(name="Enzyme", version="0.13.3")``` for any conflicts. 
