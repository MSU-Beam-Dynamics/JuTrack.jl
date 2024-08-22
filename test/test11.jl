using JuTrack

function f(k)
    RING = buildlattice()
    D1 = DRIFT(len=k)
    Q1 = QUAD(len=1.0, k1=-1.3)
    S1 = SBEND(len=1.0, angle=pi/4)
    
    add!(RING, Q1)
    add!(RING, D1)
    add!(RING, S1)
    rin = [0.0 0.01 0.0 0.0 0.0 0.0; 0.01 0.001 0.0 0.0 0.0 0.0]
    beam = Beam(rin)

    changed_id = [3]
    changed_elem = [SBEND(len=1.0, angle=k)]
    ADplinepass!(RING, beam, changed_id, changed_elem)
    return beam.r[1,:]
end

using BenchmarkTools
@btime g1 = Enzyme.autodiff(Forward, f, Duplicated, Duplicated(0.5, 1.0))

for i in 1:50
    g1 = Enzyme.autodiff(Forward, f, Duplicated, Duplicated(0.5, 1.0))
    println(i)
end

# beam = Beam([0.001 0.0 0.0 0.0 0.0 0.0; 0.001 0.001 0.0001 0.0001 0.0 0.0], energy=3e9)
# linepass!(RING, beam)
# beam1 = Beam([0.001 0.0 0.0 0.0 0.0 0.0; 0.001 0.001 0.0001 0.0001 0.0 0.0], energy=3e9)
# linepass!(RING1, beam1)
# println(beam1.r.-beam.r)