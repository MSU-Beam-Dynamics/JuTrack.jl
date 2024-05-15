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
include("path-to-the-package/src/JuTrack.jl")
using .JuTrack
```

# Known issues
* This package currently supports forward AD. Backward AD is still under development.

* For large lattice files, Julia vector array will be much more efficiency than Julia tuple. 

* Creating long lattice arrays in the differentiated function may result in an error. To avoid it, please create/load the lattice before the differentiation, and then take it as a constant variable or global variable for the differentiated function. 

* Current stable version is on Julia 1.9.4. Please up/downgrade the Julia version if there is a issue.

