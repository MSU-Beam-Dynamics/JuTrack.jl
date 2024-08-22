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
    LongitudinalWakes::Vector{LongitudinalWake}
    StrongGaussianBeams::Vector{StrongGaussianBeam}
    element_order::Vector{Tuple{Symbol, Int}}  # Stores (type, index) tuples
end

function add!(lattice::Lattice, rfca::RFCA)
    push!(lattice.rfcas, rfca)
    push!(lattice.element_order, (:RFCA, length(lattice.rfcas)))
end
function add!(lattice::Lattice, drift::DRIFT)
    push!(lattice.drifts, drift)
    push!(lattice.element_order, (:DRIFT, length(lattice.drifts)))
end
function add!(lattice::Lattice, marker::MARKER)
    push!(lattice.markers, marker)
    push!(lattice.element_order, (:MARKER, length(lattice.markers)))
end
function add!(lattice::Lattice, quad::QUAD)
    push!(lattice.quads, quad)
    push!(lattice.element_order, (:QUAD, length(lattice.quads)))
end
function add!(lattice::Lattice, kquad::KQUAD)
    push!(lattice.kquads, kquad)
    push!(lattice.element_order, (:KQUAD, length(lattice.kquads)))
end
function add!(lattice::Lattice, ksext::KSEXT)
    push!(lattice.ksexts, ksext)
    push!(lattice.element_order, (:KSEXT, length(lattice.ksexts)))
end
function add!(lattice::Lattice, koct::KOCT)
    push!(lattice.kocts, koct)
    push!(lattice.element_order, (:KOCT, length(lattice.kocts)))
end
function add!(lattice::Lattice, sbend::SBEND)
    push!(lattice.sbends, sbend)
    push!(lattice.element_order, (:SBEND, length(lattice.sbends)))
end
function add!(lattice::Lattice, thinmultipole::thinMULTIPOLE)
    push!(lattice.thinmultipoles, thinmultipole)
    push!(lattice.element_order, (:thinMULTIPOLE, length(lattice.thinmultipoles)))
end
function add!(lattice::Lattice, solenoid::SOLENOID)
    push!(lattice.solenoids, solenoid)
    push!(lattice.element_order, (:SOLENOID, length(lattice.solenoids)))
end
function add!(lattice::Lattice, corrector::CORRECTOR)
    push!(lattice.correctors, corrector)
    push!(lattice.element_order, (:CORRECTOR, length(lattice.correctors)))
end
function add!(lattice::Lattice, crabcavity::CRABCAVITY)
    push!(lattice.crabcavities, crabcavity)
    push!(lattice.element_order, (:CRABCAVITY, length(lattice.crabcavities)))
end
function add!(lattice::Lattice, spacecharge::SPACECHARGE)
    push!(lattice.spacecharges, spacecharge)
    push!(lattice.element_order, (:SPACECHARGE, length(lattice.spacecharges)))
end
function add!(lattice::Lattice, longitudinalRLCwake::LongitudinalRLCWake)
    push!(lattice.LongitudinalRLCWakes, longitudinalRLCwake)
    push!(lattice.element_order, (:LongitudinalRLCWake, length(lattice.LongitudinalRLCWakes)))
end
function add!(lattice::Lattice, longitudinalwake::LongitudinalWake)
    push!(lattice.LongitudinalWakes, longitudinalwake)
    push!(lattice.element_order, (:LongitudinalWake, length(lattice.LongitudinalWakes)))
end
function add!(lattice::Lattice, stronggaussianbeam::StrongGaussianBeam)
    push!(lattice.StrongGaussianBeams, stronggaussianbeam)
    push!(lattice.element_order, (:StrongGaussianBeam, length(lattice.StrongGaussianBeams)))
end

function buildlattice()
    return Lattice(RFCA[], DRIFT[], MARKER[], QUAD[], KQUAD[], KSEXT[], KOCT[], SBEND[], 
    thinMULTIPOLE[], SOLENOID[], CORRECTOR[], CRABCAVITY[], SPACECHARGE[], LongitudinalRLCWake[], 
    LongitudinalWake[], StrongGaussianBeam[], Tuple{Symbol, Int}[])
end
