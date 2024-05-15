# JuTrack

A Julia-based package that enables advanced auto differentiation (AD) for symplectic 6-D particle tracking in particle accelerators.

# Installation

* Install Julia 1.9.4 at [here](https://julialang.org/downloads/oldreleases/). Current version is developed based on Julia 1.9.4. Using an unexpected version may result in an error.

* Download the package
'''
git clone https://github.com/MSU-Beam-Dynamics/JuTrack.jl.git

* Install dependencies
Open Julia REPL. Move the the package folder.
'''
cd("path-to-the-package")
'''
Install the required dependencies
'''
using Pkg
Pkg.activate(".")
Pkg.instantiate()
'''

# To use this module, import it in Julia
'''
include("path-to-the-package/src/JuTrack.jl")
using .JuTrack
'''

# Known issues
1. Backward AD is not enable for most elements. Use Forward AD instead.
2. To deal with large lattice file, use Julia vector array instead of Julia tuple for better efficiency. 
3. Create long lattice arrays in the differentiated function may result in an error. To avoid it, please create/load the lattice before the differentiation, and then take it as a constant variable or global variable for the differentiated function. 
4. Current stable version is on Julia 1.9.4. Please up/downgrade the Julia version if there is a issue.

