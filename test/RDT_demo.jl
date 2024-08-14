using JuTrack
using Serialization
RING1 = deserialize("src/demo/SPEAR3/spear3.jls")

k = -1.3864672452
L = 0.3533895

function f(x::Float64, RING)
    changed_ids = [5]
    changed_elems = [KQUAD(len=0.3533895, k1=x)]
    index = findelem(RING, MARKER)
    dlist, s = ADcomputeRDT(RING, index, changed_ids, changed_elems)
    return dlist[end].h21000[1]
end
g = autodiff(Forward, f, Duplicated, Duplicated(k, 1.0), Const(RING1))
@time g1 = autodiff(Forward, f, Duplicated, Duplicated(k+1e-6, 1.0), Const(RING1))
diff = (g1[1] - g[1]) / 1e-6
println("numerical gradient: ", diff)
println("autodiff gradient: ", g[2])
