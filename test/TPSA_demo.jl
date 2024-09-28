# This demo shows how to use multivariate TPSA in JuTrack.
# The 6-D coordinates are represented by TPSA objects. 
# The TPSA objects are passed to the linepass_TPSA! function to track the particles through the lattice.
# The derivatives of the TPSA coefficients with respect to the quadrupole strength k is calculated using autodiff.

using JuTrack

function f(k)
    x = CTPS(0.0, 1, 6, 2) # create a TPSA object with 6 dimensions and 2nd order. 1 represents the index of the TPSA object.
    px = CTPS(0.0, 2, 6, 2)
    y = CTPS(0.0, 3, 6, 2)
    py = CTPS(0.0, 4, 6, 2)
    z = CTPS(0.0, 5, 6, 2)
    delta = CTPS(0.0, 6, 6, 2)

    Q1 = QUAD(len=0.2, k1=-k)
    D1 = DRIFT(len=0.2)
    Q2 = QUAD(len=0.2, k1=k)
    D2 = DRIFT(len=0.2)

    line = [Q1, D1, Q2, D2]
    rin = [x, px, y, py, z, delta]

    linepass_TPSA!(line, rin)
    return rin[1].map
end

grad_k = autodiff(Forward, f, Duplicated(10.0, 1.0))
println("Derivatives of the TPSA coefficients with respect to the quadrupole strength k: ", grad_k[1])