# JuTrack

A Julia-based package that enables advanced auto-differentiation for symplectic 6-D particle tracking in particle accelerators.

Known issues
1. Backward AD is not enable for most elements. Use Forward AD instead.
2. To deal with large lattice file, use Julia vector array instead of Julia tuple for better effciency. 
3. Do not include the code of creating a long lattice array in the differentiated function. A reasonable method is creating/loading the lattice vector before runing the differntiated function, and then take it as a constant variable for the differentiated function. 
E.g., for a function f(k, lattice) = betax, call AD with autodiff(Forward, Duplicated, Duplicated(k, 1.0), Const(lattice_vector)).

[![Build Status](https://github.com/Jinyu95/JuTrack.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Jinyu95/JuTrack.jl/actions/workflows/CI.yml?query=branch%3Amain)
