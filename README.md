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

* Install Julia at [here](https://julialang.org/downloads/oldreleases/).

* Download the package:
```
git clone https://github.com/MSU-Beam-Dynamics/JuTrack.jl
```

* Since Julia and its packages are under very active development, it is strongly recommended to create a separate environment for each piece of work.
To create a "JuTrack" environment and install all dependencies (errors may occur if not using this method to install all dependencies):
```
cd JuTrack.jl
julia --project=. -e "using Pkg; Pkg.instantiate()"
```

# Use the package in Julia
Activate the JuTrack environment we just created in the Julia code (change below path to your JuTrack directory accordingly):
```
using Pkg
Pkg.activate("/path/to/JuTrack.jl") # use the correct Julia environment
Pkg.instantiate()                   # check if all dependencies are correctly installed
```

Import JuTrack in Julia:
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
set_tps_dim(2) # 2 variables for fastTPSA
function tracking_wrt_k1(x1::DTPSAD{NVAR(), Float64}, x2::DTPSAD{NVAR(), Float64})
    D1 = DRIFT(len=DTPSAD(1.0))        # We use parametric elements in JuTrack.  
    D2 = DRIFT(len=DTPSAD(1.0))        # When any of the parameters is a DTPSAD type, the element will be ELEMENT{DTPSAD}
    Q1 = KQUAD(len=DTPSAD(1.0))        # Otherwise, the element will be ELEMENT{Float64}
    Q2 = KQUAD(len=DTPSAD(1.0))

    Q1.k1 = x1
    Q2.k1 = x2

    beam = Beam(DTPSAD.([0.1 0.0 0.0 0.0 0.0 0.0]))

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

# Precompile
Julia compiles functions the first time they are called (Just-In-Time, JIT). For complex code, this first call can take from seconds to minutes, and the cost repeats every time a new Julia process starts.

A custom sysimage precompiles those methods ahead of time, so loading and running functions has little to no latency, making Julia feel much closer to a traditional compiled language.

We provide ready-made helper scripts under the precompile/ directory. The file _agg_warmup.jl includes some sample functions. Add the function calls you want Julia to precompile in this file, and remove any you don't need to reduce the size of the sysimage.

1) Test precompile statements. Move to the JuTrack.jl folder, and type

Windows (PowerShell/CMD):
```
julia --project=. --trace-compile=precompile\trace_base.jl precompile\run_warmup.jl
```

macOS/Linux:
```
julia --project=. --trace-compile=precompile/trace_base.jl precompile/run_warmup.jl
```


2) Build the sysimage

Windows:
```
julia --project=. precompile\build_sysimage.jl
```

macOS/Linux:
```
julia --project=. precompile/build_sysimage.jl
```

This creates precompile/JuTrack_sysimage.(dll|dylib|so) for your platform.

3) Use it

Windows:
```
julia --sysimage=precompile\JuTrack_sysimage.dll your_code.jl
```

macOS:
```
julia --sysimage=precompile/JuTrack_sysimage.dylib your_code.jl
```

Linux:
```
julia --sysimage=precompile/JuTrack_sysimage.so your_code.jl
```

# Known issues
* JuTrack is actively under development. If you encounter any issues, please open an issue on GitHub or email wan@frib.msu.edu.
* Please ensure using the same OS and Julia verion when using a precompiled sysimage.