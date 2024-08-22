mutable struct Lattice
    rfcas::Vector{RFCA}
    drifts::Vector{DRIFT}
    markers::Vector{MARKER}
    quads::Vector{QUAD}
    kquads::Vector{KQUAD}
    ksexts::Vector{KSEXT}
    kocts::Vector{KOCT}
    sbends::Vector{SBEND}
    thinmultipoles::Vector{thinMULTIPOLE}
    solenoids::Vector{SOLENOID}
    correctors::Vector{CORRECTOR}
    crabcavities::Vector{CRABCAVITY}
    spacecharges::Vector{SPACECHARGE}
    LongitudinalRLCWakes::Vector{LongitudinalRLCWake}
    # LongitudinalWakes::Vector{LongitudinalWake}
    # StrongGaussianBeams::Vector{StrongGaussianBeam}
    element_order::Vector{Tuple{Int, Int}}  # Stores (type, index) tuples
    nelems::Int64
end

function add!(lattice::Lattice, rfca::RFCA)
    push!(lattice.rfcas, rfca)
    push!(lattice.element_order, (1, length(lattice.rfcas)))
    lattice.nelems += 1
end
function add!(lattice::Lattice, drift::DRIFT)
    push!(lattice.drifts, drift)
    push!(lattice.element_order, (2, length(lattice.drifts)))
    lattice.nelems += 1
end
function add!(lattice::Lattice, marker::MARKER)
    push!(lattice.markers, marker)
    push!(lattice.element_order, (3, length(lattice.markers)))
    lattice.nelems += 1
end
function add!(lattice::Lattice, quad::QUAD)
    push!(lattice.quads, quad)
    push!(lattice.element_order, (4, length(lattice.quads)))
    lattice.nelems += 1
end
function add!(lattice::Lattice, kquad::KQUAD)
    push!(lattice.kquads, kquad)
    push!(lattice.element_order, (5, length(lattice.kquads)))
    lattice.nelems += 1
end
function add!(lattice::Lattice, ksext::KSEXT)
    push!(lattice.ksexts, ksext)
    push!(lattice.element_order, (6, length(lattice.ksexts)))
    lattice.nelems += 1
end
function add!(lattice::Lattice, koct::KOCT)
    push!(lattice.kocts, koct)
    push!(lattice.element_order, (7, length(lattice.kocts)))
    lattice.nelems += 1
end
function add!(lattice::Lattice, sbend::SBEND)
    push!(lattice.sbends, sbend)
    push!(lattice.element_order, (8, length(lattice.sbends)))
    lattice.nelems += 1
end
function add!(lattice::Lattice, thinmultipole::thinMULTIPOLE)
    push!(lattice.thinmultipoles, thinmultipole)
    push!(lattice.element_order, (9, length(lattice.thinmultipoles)))
    lattice.nelems += 1
end
function add!(lattice::Lattice, solenoid::SOLENOID)
    push!(lattice.solenoids, solenoid)
    push!(lattice.element_order, (10, length(lattice.solenoids)))
    lattice.nelems += 1
end
function add!(lattice::Lattice, corrector::CORRECTOR)
    push!(lattice.correctors, corrector)
    push!(lattice.element_order, (11, length(lattice.correctors)))
    lattice.nelems += 1
end
function add!(lattice::Lattice, crabcavity::CRABCAVITY)
    push!(lattice.crabcavities, crabcavity)
    push!(lattice.element_order, (12, length(lattice.crabcavities)))
    lattice.nelems += 1
end
function add!(lattice::Lattice, spacecharge::SPACECHARGE)
    push!(lattice.spacecharges, spacecharge)
    push!(lattice.element_order, (13, length(lattice.spacecharges)))
    lattice.nelems += 1
end
function add!(lattice::Lattice, longitudinalRLCwake::LongitudinalRLCWake)
    push!(lattice.LongitudinalRLCWakes, longitudinalRLCwake)
    push!(lattice.element_order, (14, length(lattice.LongitudinalRLCWakes)))
    lattice.nelems += 1
end
function add!(lattice::Lattice, longitudinalwake::LongitudinalWake)
    push!(lattice.LongitudinalWakes, longitudinalwake)
    push!(lattice.element_order, (15, length(lattice.LongitudinalWakes)))
    lattice.nelems += 1
end
function add!(lattice::Lattice, stronggaussianbeam::StrongGaussianBeam)
    push!(lattice.StrongGaussianBeams, stronggaussianbeam)
    push!(lattice.element_order, (16, length(lattice.StrongGaussianBeams)))
    lattice.nelems += 1
end

# function buildlattice()
#     return Lattice(RFCA[], DRIFT[], MARKER[], QUAD[], KQUAD[], KSEXT[], KOCT[], SBEND[], 
#     thinMULTIPOLE[], SOLENOID[], CORRECTOR[], CRABCAVITY[], SPACECHARGE[], LongitudinalRLCWake[], 
#     LongitudinalWake[], StrongGaussianBeam[], Tuple{Int64, Int64}[], 0)
# end
function buildlattice()
    return Lattice([RFCA()], [DRIFT()], [MARKER()], [QUAD()], [KQUAD()], [KSEXT()], [KOCT()], [SBEND()], 
    [thinMULTIPOLE()], [SOLENOID()], [CORRECTOR()], [CRABCAVITY()], [SPACECHARGE()], [LongitudinalRLCWake()], [(0, 0)], 0)
end