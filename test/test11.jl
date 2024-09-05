using JuTrack
    RING = Lattice(nelems=3)
    D1 = DRIFT(len=-1.38)
    Q1 = KQUAD(len=1.0, k1=-1.3)
    S1 = SBEND(len=1.0, angle=0.5)
    
    add!(RING, Q1)
    add!(RING, D1)
    add!(RING, S1)
function f(k)
    # RING = Lattice(nelems=3)
    # D1 = DRIFT(len=k)
    # Q1 = KQUAD(len=1.0, k1=-1.3)
    # S1 = SBEND(len=1.0, angle=0.5)
    
    # add!(RING, Q1)
    # add!(RING, D1)
    # add!(RING, S1)
    rin = [0.0 0.0001 0.0 0.0 0.0 0.0; 0.01 0.001 0.0 0.0 0.0 0.0]
    beam = Beam(rin)
    # RING1 = [Q1, D1, S1]
    # changed_id = [5]
    # changed_elem = [KQUAD(len=0.3533895, k1=k)]
    RING = spear()
    RING.kquads[1].k1 = k
    plinepass!(RING, beam)
    return beam.r[1,:]
end

# using BenchmarkTools
# @btime g1 = Enzyme.autodiff(Forward, f, Duplicated, Duplicated(-1.38, 1.0))

for i in 1:1
    g1 = Enzyme.autodiff(Forward, f, Duplicated, Duplicated(-1.38, 1.0))
    println("iteration $i: ", g1[1])
end

# beam = Beam([0.001 0.0 0.0 0.0 0.0 0.0; 0.001 0.001 0.0001 0.0001 0.0 0.0], energy=3e9)
# linepass!(RING, beam)
# beam1 = Beam([0.001 0.0 0.0 0.0 0.0 0.0; 0.001 0.001 0.0001 0.0001 0.0 0.0], energy=3e9)
# linepass!(RING1, beam1)
# println(beam1.r.-beam.r)