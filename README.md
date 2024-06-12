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

# Parallel computing setting
This package uses Julia's multi-threading for parallel computing. Before using it, one has to set up the number of threads.
On Linux or macOS:
```
export JULIA_NUM_THREADS=N
``` 
or add the above line to your .bashrc file.

On Windows:
Maually add JULIA_NUM_THREADS as the variable name and N (the number of threads you wish to use) as the variable value in system variable.

It is recommended to use VS code to permanently set up the Julia environment. 
Open the Command Palette (Ctrl+Shift+P or Cmd+Shift+P on macOS). 
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

* Current stable version is on Julia 1.9.4. Please up/downgrade the Julia version if there is a issue.

