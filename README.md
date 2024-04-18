# JuTrack

A Julia-based package that enables advanced auto differentiation (AD) for symplectic 6-D particle tracking in particle accelerators.

Known issues
1. Backward AD is not enable for most elements. Use Forward AD instead.
2. To deal with large lattice file, use Julia vector array instead of Julia tuple for better efficiency. 
3. Create long lattice arrays in the differentiated function may result in an error. To avoid it, please create/load the lattice before the differentiation, and then take it as a constant variable or global variable for the differentiated function. 
4. Current stable version is on Julia 1.9.4. Please up/downgrade the Julia version if there is a issue.

