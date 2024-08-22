function linepass!(lattice::Lattice, particles::Beam)
    np = particles.nmacro
    particles6 = matrix_to_array(particles.r)
    for (typ, idx) in lattice.element_order
        if typ == :RFCA
            pass!(lattice.rfcas[idx], particles6, np, particles)
        elseif typ == :DRIFT
            pass!(lattice.drifts[idx], particles6, np, particles)
        elseif typ == :MARKER
            pass!(lattice.markers[idx], particles6, np, particles)
        elseif typ == :QUAD
            pass!(lattice.quads[idx], particles6, np, particles)
        elseif typ == :KQUAD
            pass!(lattice.kquads[idx], particles6, np, particles)
        elseif typ == :KSEXT
            pass!(lattice.ksexts[idx], particles6, np, particles)
        elseif typ == :KOCT
            pass!(lattice.kocts[idx], particles6, np, particles)
        elseif typ == :SBEND
            pass!(lattice.sbends[idx], particles6, np, particles)
        elseif typ == :thinMULTIPOLE
            pass!(lattice.thinmultipoles[idx], particles6, np, particles)
        elseif typ == :SOLENOID
            pass!(lattice.solenoids[idx], particles6, np, particles)
        elseif typ == :CORRECTOR
            pass!(lattice.correctors[idx], particles6, np, particles)
        elseif typ == :CRABCAVITY
            pass!(lattice.crabcavities[idx], particles6, np, particles)
        elseif typ == :SPACECHARGE
            pass!(lattice.spacecharges[idx], particles6, np, particles)
        elseif typ == :LongitudinalRLCWake
            pass!(lattice.LongitudinalRLCWakes[idx], particles6, np, particles)
        elseif typ == :LongitudinalWake
            pass!(lattice.LongitudinalWakes[idx], particles6, np, particles)
        else
            println("Unknown element type: $typ")
        end
    end
    rout = array_to_matrix(particles6, np)
    particles.r = rout
    return nothing
end

function ADlinepass!(lattice::Lattice, particles::Beam, id::Vector{Int}, elems::Vector)
    np = particles.nmacro
    particles6 = matrix_to_array(particles.r)
    c = 1
    for i in 1:length(lattice.element_order)
        if i in id
            pass!(elems[c], particles6, np, particles)
            c += 1
        else
            if lattice.element_order[i][1] == :RFCA
                pass!(lattice.rfcas[lattice.element_order[i][2]], particles6, np, particles)
            elseif lattice.element_order[i][1] == :DRIFT
                pass!(lattice.drifts[lattice.element_order[i][2]], particles6, np, particles)
            elseif lattice.element_order[i][1] == :MARKER
                pass!(lattice.markers[lattice.element_order[i][2]], particles6, np, particles)
            elseif lattice.element_order[i][1] == :QUAD
                pass!(lattice.quads[lattice.element_order[i][2]], particles6, np, particles)
            elseif lattice.element_order[i][1] == :KQUAD
                pass!(lattice.kquads[lattice.element_order[i][2]], particles6, np, particles)
            elseif lattice.element_order[i][1] == :KSEXT
                pass!(lattice.ksexts[lattice.element_order[i][2]], particles6, np, particles)
            elseif lattice.element_order[i][1] == :KOCT
                pass!(lattice.kocts[lattice.element_order[i][2]], particles6, np, particles)
            elseif lattice.element_order[i][1] == :SBEND
                pass!(lattice.sbends[lattice.element_order[i][2]], particles6, np, particles)
            elseif lattice.element_order[i][1] == :thinMULTIPOLE
                pass!(lattice.thinmultipoles[lattice.element_order[i][2]], particles6, np, particles)
            elseif lattice.element_order[i][1] == :SOLENOID
                pass!(lattice.solenoids[lattice.element_order[i][2]], particles6, np, particles)
            elseif lattice.element_order[i][1] == :CORRECTOR
                pass!(lattice.correctors[lattice.element_order[i][2]], particles6, np, particles)
            elseif lattice.element_order[i][1] == :CRABCAVITY
                pass!(lattice.crabcavities[lattice.element_order[i][2]], particles6, np, particles)
            elseif lattice.element_order[i][1] == :SPACECHARGE
                pass!(lattice.spacecharges[lattice.element_order[i][2]], particles6, np, particles)
            elseif lattice.element_order[i][1] == :LongitudinalRLCWake
                pass!(lattice.LongitudinalRLCWakes[lattice.element_order[i][2]], particles6, np, particles)
            elseif lattice.element_order[i][1] == :LongitudinalWake
                pass!(lattice.LongitudinalWakes[lattice.element_order[i][2]], particles6, np, particles)
            else
                println("Unknown element type: $(lattice.element_order[i][1])")
            end
        end
    end
    rout = array_to_matrix(particles6, np)
    particles.r = rout
    return c
end

function ringpass!(lattice::Lattice, particles::Beam, nturn::Int)
    for i in 1:nturn
        linepass!(lattice, particles)    
    end
    return nothing
end
function ADringpass!(lattice::Lattice, particles::Beam, nturn::Int, id::Vector{Int}, elems::Vector)
    for i in 1:nturn
        ADlinepass!(lattice, particles, id, elems)    
    end
    return nothing
end

# multi-threading
function plinepass!(lattice::Lattice, particles::Beam)
    np = particles.nmacro
    particles6 = matrix_to_array(particles.r)
    for (typ, idx) in lattice.element_order
        if typ == :RFCA
            pass_P!(lattice.rfcas[idx], particles6, np, particles)
        elseif typ == :DRIFT
            pass_P!(lattice.drifts[idx], particles6, np, particles)
        elseif typ == :MARKER
            pass_P!(lattice.markers[idx], particles6, np, particles)
        elseif typ == :QUAD
            pass_P!(lattice.quads[idx], particles6, np, particles)
        elseif typ == :KQUAD
            pass_P!(lattice.kquads[idx], particles6, np, particles)
        elseif typ == :KSEXT
            pass_P!(lattice.ksexts[idx], particles6, np, particles)
        elseif typ == :KOCT
            pass_P!(lattice.kocts[idx], particles6, np, particles)
        elseif typ == :SBEND
            pass_P!(lattice.sbends[idx], particles6, np, particles)
        elseif typ == :thinMULTIPOLE
            pass_P!(lattice.thinmultipoles[idx], particles6, np, particles)
        elseif typ == :SOLENOID
            pass_P!(lattice.solenoids[idx], particles6, np, particles)
        elseif typ == :CORRECTOR
            pass_P!(lattice.correctors[idx], particles6, np, particles)
        elseif typ == :CRABCAVITY
            pass_P!(lattice.crabcavities[idx], particles6, np, particles)
        elseif typ == :SPACECHARGE
            pass_P!(lattice.spacecharges[idx], particles6, np, particles)
        elseif typ == :LongitudinalRLCWake
            pass_P!(lattice.LongitudinalRLCWakes[idx], particles6, np, particles)
        elseif typ == :LongitudinalWake
            pass_P!(lattice.LongitudinalWakes[idx], particles6, np, particles)
        else
            println("Unknown element type: $typ")
        end
    end
    rout = array_to_matrix(particles6, np)
    particles.r = rout
    return nothing
end

function ADplinepass!(lattice::Lattice, particles::Beam, id::Vector{Int}, elems::Vector)
    np = particles.nmacro
    particles6 = matrix_to_array(particles.r)
    c = 1
    for i in 1:length(lattice.element_order)
        if i in id
            pass_P!(elems[c], particles6, np, particles)
            c += 1
        else
            if lattice.element_order[i][1] == :RFCA
                pass_P!(lattice.rfcas[lattice.element_order[i][2]], particles6, np, particles)
            elseif lattice.element_order[i][1] == :DRIFT
                pass_P!(lattice.drifts[lattice.element_order[i][2]], particles6, np, particles)
            elseif lattice.element_order[i][1] == :MARKER
                pass_P!(lattice.markers[lattice.element_order[i][2]], particles6, np, particles)
            elseif lattice.element_order[i][1] == :QUAD
                pass_P!(lattice.quads[lattice.element_order[i][2]], particles6, np, particles)
            elseif lattice.element_order[i][1] == :KQUAD
                pass_P!(lattice.kquads[lattice.element_order[i][2]], particles6, np, particles)
            elseif lattice.element_order[i][1] == :KSEXT
                pass_P!(lattice.ksexts[lattice.element_order[i][2]], particles6, np, particles)
            elseif lattice.element_order[i][1] == :KOCT
                pass_P!(lattice.kocts[lattice.element_order[i][2]], particles6, np, particles)
            elseif lattice.element_order[i][1] == :SBEND
                pass_P!(lattice.sbends[lattice.element_order[i][2]], particles6, np, particles)
            elseif lattice.element_order[i][1] == :thinMULTIPOLE
                pass_P!(lattice.thinmultipoles[lattice.element_order[i][2]], particles6, np, particles)
            elseif lattice.element_order[i][1] == :SOLENOID
                pass_P!(lattice.solenoids[lattice.element_order[i][2]], particles6, np, particles)
            elseif lattice.element_order[i][1] == :CORRECTOR
                pass_P!(lattice.correctors[lattice.element_order[i][2]], particles6, np, particles)
            elseif lattice.element_order[i][1] == :CRABCAVITY
                pass_P!(lattice.crabcavities[lattice.element_order[i][2]], particles6, np, particles)
            elseif lattice.element_order[i][1] == :SPACECHARGE
                pass_P!(lattice.spacecharges[lattice.element_order[i][2]], particles6, np, particles)
            elseif lattice.element_order[i][1] == :LongitudinalRLCWake
                pass_P!(lattice.LongitudinalRLCWakes[lattice.element_order[i][2]], particles6, np, particles)
            elseif lattice.element_order[i][1] == :LongitudinalWake
                pass_P!(lattice.LongitudinalWakes[lattice.element_order[i][2]], particles6, np, particles)
            else
                println("Unknown element type: $(lattice.element_order[i][1])")
            end
        end
    end
    rout = array_to_matrix(particles6, np)
    particles.r = rout
    return c
end

function pringpass!(lattice::Lattice, particles::Beam, nturn::Int)
    for i in 1:nturn
        plinepass!(lattice, particles)    
    end
    return nothing
end
function ADpringpass!(lattice::Lattice, particles::Beam, nturn::Int, id::Vector{Int}, elems::Vector)
    for i in 1:nturn
        ADplinepass!(lattice, particles, id, elems)    
    end
    return nothing
end

# tracking function for TPSA
function linepass_TPSA!(line::Lattice, rin::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}) where {T, TPS_Dim, Max_TPS_Degree}
    if length(rin) != 6
        error("The length of TPSA must be 6")
    end
    for (typ, idx) in line.element_order
        if typ == :RFCA
            pass_TPSA!(line.rfcas[idx], rin)
        elseif typ == :DRIFT
            pass_TPSA!(line.drifts[idx], rin)
        elseif typ == :MARKER
            pass_TPSA!(line.markers[idx], rin)
        elseif typ == :QUAD
            pass_TPSA!(line.quads[idx], rin)
        elseif typ == :KQUAD
            pass_TPSA!(line.kquads[idx], rin)
        elseif typ == :KSEXT
            pass_TPSA!(line.ksexts[idx], rin)
        elseif typ == :KOCT
            pass_TPSA!(line.kocts[idx], rin)
        elseif typ == :SBEND
            pass_TPSA!(line.sbends[idx], rin)
        elseif typ == :thinMULTIPOLE
            pass_TPSA!(line.thinmultipoles[idx], rin)
        elseif typ == :SOLENOID
            pass_TPSA!(line.solenoids[idx], rin)
        elseif typ == :CORRECTOR
            pass_TPSA!(line.correctors[idx], rin)
        elseif typ == :CRABCAVITY
            pass_TPSA!(line.crabcavities[idx], rin)
        elseif typ == :SPACECHARGE
            pass_TPSA!(line.spacecharges[idx], rin)
        elseif typ == :LongitudinalRLCWake
            pass_TPSA!(line.LongitudinalRLCWakes[idx], rin)
        elseif typ == :LongitudinalWake
            pass_TPSA!(line.LongitudinalWakes[idx], rin)
        else
            println("Unknown element type: $typ")
        end
    end
    return nothing
end

function ringpass_TPSA!(line::Lattice, rin::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}, nturn::Int) where {T, TPS_Dim, Max_TPS_Degree}
    if length(rin) != 6
        error("The length of TPSA must be 6")
    end
    for i in 1:nturn
        linepass_TPSA!(line, rin)    
    end
    return nothing
end

function ADlinepass_TPSA!(line::Lattice, rin::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}, changed_idx::Vector, changed_ele::Vector) where {T, TPS_Dim, Max_TPS_Degree}
    if length(rin) != 6
        error("The length of TPSA must be 6")
    end
    count = 1
    for i in eachindex(line.element_order)
        if i in changed_idx
            pass_TPSA!(changed_ele[count], rin)
            count += 1
        else
            if line.element_order[i][1] == :RFCA
                pass_TPSA!(line.rfcas[line.element_order[i][2]], rin)
            elseif line.element_order[i][1] == :DRIFT
                pass_TPSA!(line.drifts[line.element_order[i][2]], rin)
            elseif line.element_order[i][1] == :MARKER
                pass_TPSA!(line.markers[line.element_order[i][2]], rin)
            elseif line.element_order[i][1] == :QUAD
                pass_TPSA!(line.quads[line.element_order[i][2]], rin)
            elseif line.element_order[i][1] == :KQUAD
                pass_TPSA!(line.kquads[line.element_order[i][2]], rin)
            elseif line.element_order[i][1] == :KSEXT
                pass_TPSA!(line.ksexts[line.element_order[i][2]], rin)
            elseif line.element_order[i][1] == :KOCT
                pass_TPSA!(line.kocts[line.element_order[i][2]], rin)
            elseif line.element_order[i][1] == :SBEND
                pass_TPSA!(line.sbends[line.element_order[i][2]], rin)
            elseif line.element_order[i][1] == :thinMULTIPOLE
                pass_TPSA!(line.thinmultipoles[line.element_order[i][2]], rin)
            elseif line.element_order[i][1] == :SOLENOID
                pass_TPSA!(line.solenoids[line.element_order[i][2]], rin)
            elseif line.element_order[i][1] == :CORRECTOR
                pass_TPSA!(line.correctors[line.element_order[i][2]], rin)
            elseif line.element_order[i][1] == :CRABCAVITY
                pass_TPSA!(line.crabcavities[line.element_order[i][2]], rin)
            elseif line.element_order[i][1] == :SPACECHARGE
                pass_TPSA!(line.spacecharges[line.element_order[i][2]], rin)
            elseif line.element_order[i][1] == :LongitudinalRLCWake
                pass_TPSA!(line.LongitudinalRLCWakes[line.element_order[i][2]], rin)
            elseif line.element_order[i][1] == :LongitudinalWake
                pass_TPSA!(line.LongitudinalWakes[line.element_order[i][2]], rin)
            else
                println("Unknown element type: $(line.element_order[i][1])")
            end
        end
    end
    return count
end

function ADringpass_TPSA!(line::Lattice, rin::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}, nturn::Int, changed_idx::Vector, changed_ele::Vector) where {T, TPS_Dim, Max_TPS_Degree}
    if length(rin) != 6
        error("The length of TPSA must be 6")
    end
    for i in 1:nturn
        ADlinepass_TPSA!(line, rin, changed_idx, changed_ele)    
    end
    return nothing
end
