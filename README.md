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

* Install Julia at [here](https://julialang.org/downloads/oldreleases/) (1.10.4 is preferred).

* Download the package for offline isntallation.
* Or using online installation via
```
using Pkg
Pkg.add(url="https://github.com/MSU-Beam-Dynamics/JuTrack.jl")
```

* Installation of Enzyme. v0.13.3 is perferred. 
```
Pkg.add(name="Enzyme", version="0.13.3")
```

# Import the package in Julia
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
Obtain derivatives of tracking result w.r.t the quadrupole strength k1
```
function tracking_wrt_k1(x)
    D1 = DRIFT(len=1.0)
    D2 = DRIFT(len=1.0)
    Q1 = KQUAD(len=1.0, k1=x) 
    Q2 = KQUAD(len=1.0, k1=0.3)

    beam = Beam([0.1 0.0 0.0 0.0 0.0 0.0], energy=3.5e9)

    # !!! avoid creating large lattice in the differentiable function.
    # !!! load or create your lattice outside the function if it is large.
    LINE = [D1, Q1, D2, Q2] 
    linepass!(LINE, beam)
    return beam.r
end
k1 = -0.9
derivatives, results = autodiff(ForwardWithPrimal, tracking_wrt_k1, Duplicated(k1, 1.0))
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
