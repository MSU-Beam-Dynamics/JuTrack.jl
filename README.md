# JuTrack

A Julia-based package that enables advanced auto differentiation (AD) for symplectic 6-D particle tracking in particle accelerators.

# Installation

* Install Julia 1.9.4 at [here](https://julialang.org/downloads/oldreleases/). Current version is developed based on Julia 1.9.4. Using an unexpected version may result in an error.

* Download the package
```
git clone https://github.com/MSU-Beam-Dynamics/JuTrack.jl.git
```

* Open Julia REPL. Move to the package folder.
```
cd("path-to-the-package")
```

* Install the required dependencies
```
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```

# Import the package in Julia
```
using Pkg # optional if you are using another Julia environment
Pkg.activate(".") # optional if you are using another Julia environment
include("path-to-the-package/src/JuTrack.jl")
using .JuTrack
```
# Lattice definition
```
D1 = DRIFT(len=1.0)
D2 = DRIFT(len=1.0)
D3 = DRIFT(len=1.0)
D4 = DRIFT(len=1.0)
D5 = DRIFT(len=1.0)
Q1 = KQUAD(len=1.0, k1=-0.9) 
Q2 = KQUAD(len=1.0, k1=0.3)
B1 = SBEND(len= 0.6, angle=pi/15.0)
B2 = SBEND(len= 0.6, angle=-pi/15.0)

# Create the lattice as a Julia vector
LINE = [D1, Q1, D2, B1, D3, Q2, D4, B2, D5, B2, D4, Q2, D3, B1, D2, Q1, D1]
```

# Particle tracking
```
# Particles' coordinates are represent as a N * 6 matrix, saved in beam.r
particles = rand(10, 6) / 1000
# Create the beam
beam = Beam(particles, energy=3.5e9)
# Tracking
linepass!(LINE, beam) # or ringpass!(RING, beam, nturns) for multi-turn tracking
println(beam.r) 
```

# Optics calculation
```
# Obtain periodic Twiss parameters of a ring accelerator
twi = twissring(RING, 0.0, 1)
```

# Aotumatic differentiation
```
# Obtain derivatives of tracking result w.r.t the initial coordinate x0
# The warm-up time may be long due to Julia's JIT feature. The computation will be fast after the compilation is done. 
function tracking(x)
    beam = Beam([x 0.0 0.0 0.0 0.0 0.0], energy=3.5e9)
    ringpass!(LINE, beam, 1)
    return beam.r
end
x0 = 0.01
result = autodiff(Forward, tracking, Duplicated, Duplicated(x0, 1.0))
```
```
# obtain 4*4 Jacobian matrix of a lattice of specific initial condition [0.01 0.0 0.01 0.0]
function obtain_jacobian(x)
    beam = Beam([x[1] x[2] x[3] x[4] 0.0 0.0], energy=3.5e9)
    ringpass!(LINE, beam, 1)
    return beam.r[1:4]
end
g = jacobian(Forward, obtain_jacobian, [0.01, 0.0, 0.01, 0.0], Val(4))
```

# Parallel computing setting
This package uses Julia's multi-threading for parallel computing. Before using it, one has to set up the number of threads.
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
    "JULIA_NUM_THREADS": "4" // Adjust the number to the desired number of threads
}
```

To check if the multi-threading is set up correctly, open the Julia REPL, and type:
```
println(Threads.nthreads())
```

Parallel computing is available for multi-particle tracking using:
```
plinepass(beamline, beam)
```
or 
```
pringpass(beamline, beam, nturns)
```

# Known issues
* This package currently supports forward AD. Backward AD is still under development.

* Creating long lattice vectors in the differentiated function may result in an error. To avoid it, please create/load the lattice before the differentiation, and then take it as a constant variable or global variable for the differentiated function. 

* Terminated with error code (-1073741819). This is a Windows-specific error code, indicating possible conflicts of the package with the memomry management of Windows system. It is always recommended to run it on Linux platform for better performance.

* Current stable version is on Julia 1.9.4. Please up/downgrade the Julia version if there is a issue.

