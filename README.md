# JuTrack

A Julia-based package that enables advanced auto differentiation (AD) for symplectic 6-D particle tracking in particle accelerators.
A mannual can be found at [here](https://msu-beam-dynamics.github.io/JuTrack.jl/).

# Installation

* Install Julia at [here](https://julialang.org/downloads/oldreleases/) (1.10.4 is preferred).

* Online installation of the package via
```
using Pkg
Pkg.add(url="https://github.com/MSU-Beam-Dynamics/JuTrack.jl")
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
D5 = DRIFT(name="D5", len=1.0)
Q1 = KQUAD(name="Q1", len=1.0, k1=-0.9) 
Q2 = KQUAD(name="Q2", len=1.0, k1=0.3)
B1 = SBEND(name="B1", len= 0.6, angle=pi/15.0)
B2 = SBEND(name="B2", len= 0.6, angle=-pi/15.0)
```
Create the lattice as a Julia vector
```
LINE = [D1, Q1, D2, B1, D3, Q2, D4, B2, D5, B2, D4, Q2, D3, B1, D2, Q1, D1]
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
twi = twissring(RING, 0.0, 1)
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

    # !!! Creating a large lattice in the function you try to differentiate 
    # !!! will slow the computation and may result in an error.
    # !!! Create/load the lattice outside of the function if it is large.
    LINE = [D1, Q1, D2, Q2] 
    linepass!(LINE, beam)
    return beam.r
end
k1 = -0.9
derivatives, results = autodiff(ForwardWithPrimal, tracking_wrt_k1, Duplicated(k1, 1.0))
```

# Parallel computation setting
Multi-threading is available for multi-particle tracking. Before using it, one has to set up the number of threads.

On Linux or macOS:
```
export JULIA_NUM_THREADS=N
``` 
or add the above line to your .bashrc file. Change N to the desired number of threads.

On Windows:
Maually add JULIA_NUM_THREADS as the variable name and N (the number of threads) as the variable value in system variable.

It is recommended to use Visual Studio Code to permanently set up the Julia environment. 
Open the command palette (Ctrl+Shift+P or Cmd+Shift+P on macOS). 
Type "Preferences: Open Settings (JSON)" and select it to open the settings file. 
Add or modify the Julia settings to include the environment variable like this:
```
"julia.executablePath": "path/to/julia", // change it to your Julia path
"julia.environmentVariables": {
    "JULIA_NUM_THREADS": "48" // Adjust the number to the desired number of threads
}
```

To check if the multi-threading is set up correctly, open the Julia REPL, and type:
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
* This package currently supports forward AD. Backward AD is still under development.