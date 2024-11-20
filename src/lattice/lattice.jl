# temporary lattice definition. Not used in the current version
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
    rfca_index::Int64
    drift_index::Int64
    marker_index::Int64
    quad_index::Int64
    kquad_index::Int64
    ksext_index::Int64
    koct_index::Int64
    sbend_index::Int64
    thinmultipole_index::Int64
    solenoid_index::Int64
    corrector_index::Int64
    crabcavity_index::Int64
    spacecharge_index::Int64
    LongitudinalRLCWake_index::Int64
    element_order::Vector{Tuple{Int64, Int64}}  # Stores (type, index) tuples
    nelems::Int64
    function Lattice(;nelems::Int64=0)
        rfcas=Vector{RFCA}(undef, nelems)
        drifts=Vector{DRIFT}(undef, nelems)
        markers=Vector{MARKER}(undef, nelems)
        quads=Vector{QUAD}(undef, nelems)
        kquads=Vector{KQUAD}(undef, nelems)
        ksexts=Vector{KSEXT}(undef, nelems)
        kocts=Vector{KOCT}(undef, nelems)
        sbends=Vector{SBEND}(undef, nelems)
        thinmultipoles=Vector{thinMULTIPOLE}(undef, nelems)
        solenoids=Vector{SOLENOID}(undef, nelems)
        correctors=Vector{CORRECTOR}(undef, nelems)
        crabcavities=Vector{CRABCAVITY}(undef, nelems)
        spacecharges=Vector{SPACECHARGE}(undef, nelems)
        LongitudinalRLCWakes=Vector{LongitudinalRLCWake}(undef, nelems)
        element_order=Vector{Tuple{Int64, Int64}}(undef, nelems)
        # nelems::Int64=0
        new(rfcas, drifts, markers, quads, kquads, ksexts, kocts, sbends, thinmultipoles, solenoids, correctors, crabcavities, 
            spacecharges, LongitudinalRLCWakes, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, element_order, 0)
    end
end

function add!(lattice::Lattice, elem::RFCA)
    index = lattice.rfca_index
    lattice.rfcas[index] = elem
    lattice.element_order[lattice.nelems + 1] = (1, index)
    lattice.nelems += 1
    lattice.rfca_index += 1
end
function add!(lattice::Lattice, elem::DRIFT)
    index = lattice.drift_index
    lattice.drifts[index] = elem
    lattice.element_order[lattice.nelems + 1] = (2, index)
    lattice.nelems += 1
    lattice.drift_index += 1
end
function add!(lattice::Lattice, elem::MARKER)
    index = lattice.marker_index
    lattice.markers[index] = elem
    lattice.element_order[lattice.nelems + 1] = (3, index)
    lattice.nelems += 1
    lattice.marker_index += 1
end
function add!(lattice::Lattice, elem::QUAD)
    index = lattice.quad_index
    lattice.quads[index] = elem
    lattice.element_order[lattice.nelems + 1] = (4, index)
    lattice.nelems += 1
    lattice.quad_index += 1
end
function add!(lattice::Lattice, elem::KQUAD)
    index = lattice.kquad_index
    lattice.kquads[index] = elem
    lattice.element_order[lattice.nelems + 1] = (5, index)
    lattice.nelems += 1
    lattice.kquad_index += 1
end
function add!(lattice::Lattice, elem::KSEXT)
    index = lattice.ksext_index
    lattice.ksexts[index] = elem
    lattice.element_order[lattice.nelems + 1] = (6, index)
    lattice.nelems += 1
    lattice.ksext_index += 1
end
function add!(lattice::Lattice, elem::KOCT)
    index = lattice.koct_index
    lattice.kocts[index] = elem
    lattice.element_order[lattice.nelems + 1] = (7, index)
    lattice.nelems += 1
    lattice.koct_index += 1
end
function add!(lattice::Lattice, elem::SBEND)
    index = lattice.sbend_index
    lattice.sbends[index] = elem
    lattice.element_order[lattice.nelems + 1] = (8, index)
    lattice.nelems += 1
    lattice.sbend_index += 1
end
function add!(lattice::Lattice, elem::thinMULTIPOLE)
    index = lattice.thinmultipole_index
    lattice.thinmultipoles[index] = elem
    lattice.element_order[lattice.nelems + 1] = (9, index)
    lattice.nelems += 1
    lattice.thinmultipole_index += 1
end
function add!(lattice::Lattice, elem::SOLENOID)
    index = lattice.solenoid_index
    lattice.solenoids[index] = elem
    lattice.element_order[lattice.nelems + 1] = (10, index)
    lattice.nelems += 1
    lattice.solenoid_index += 1
end
function add!(lattice::Lattice, elem::CORRECTOR)
    index = lattice.corrector_index
    lattice.correctors[index] = elem
    lattice.element_order[lattice.nelems + 1] = (11, index)
    lattice.nelems += 1
    lattice.corrector_index += 1
end
function add!(lattice::Lattice, elem::CRABCAVITY)
    index = lattice.crabcavity_index
    lattice.crabcavities[index] = elem
    lattice.element_order[lattice.nelems + 1] = (12, index)
    lattice.nelems += 1
    lattice.crabcavity_index += 1
end
function add!(lattice::Lattice, elem::SPACECHARGE)
    index = lattice.spacecharge_index
    lattice.spacecharges[index] = elem
    lattice.element_order[lattice.nelems + 1] = (13, index)
    lattice.nelems += 1
    lattice.spacecharge_index += 1
end
function add!(lattice::Lattice, elem::LongitudinalRLCWake)
    index = lattice.LongitudinalRLCWake_index
    lattice.LongitudinalRLCWakes[index] = elem
    lattice.element_order[lattice.nelems + 1] = (14, index)
    lattice.nelems += 1
    lattice.LongitudinalRLCWake_index += 1
end
