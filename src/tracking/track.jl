function linepass!(lattice::Lattice, particles::Beam)
    np = particles.nmacro
    particles6 = matrix_to_array(particles.r)
    for i in 1:lattice.nelems
        if lattice.element_order[i][1] == 1
            pass!(lattice.rfcas[lattice.element_order[i][2]], particles6, np, particles)
        elseif lattice.element_order[i][1] == 2
            pass!(lattice.drifts[lattice.element_order[i][2]], particles6, np, particles)
        elseif lattice.element_order[i][1] == 3
            pass!(lattice.markers[lattice.element_order[i][2]], particles6, np, particles)
        elseif lattice.element_order[i][1] == 4
            pass!(lattice.quads[lattice.element_order[i][2]], particles6, np, particles)
        elseif lattice.element_order[i][1] == 5
            pass!(lattice.kquads[lattice.element_order[i][2]], particles6, np, particles)
        elseif lattice.element_order[i][1] == 6
            pass!(lattice.ksexts[lattice.element_order[i][2]], particles6, np, particles)
        elseif lattice.element_order[i][1] == 7
            pass!(lattice.kocts[lattice.element_order[i][2]], particles6, np, particles)
        elseif lattice.element_order[i][1] == 8
            pass!(lattice.sbends[lattice.element_order[i][2]], particles6, np, particles)
        elseif lattice.element_order[i][1] == 9
            pass!(lattice.thinmultipoles[lattice.element_order[i][2]], particles6, np, particles)
        elseif lattice.element_order[i][1] == 10
            pass!(lattice.solenoids[lattice.element_order[i][2]], particles6, np, particles)
        elseif lattice.element_order[i][1] == 11
            pass!(lattice.correctors[lattice.element_order[i][2]], particles6, np, particles)
        elseif lattice.element_order[i][1] == 12
            pass!(lattice.crabcavities[lattice.element_order[i][2]], particles6, np, particles)
        elseif lattice.element_order[i][1] == 13
            pass!(lattice.spacecharges[lattice.element_order[i][2]], particles6, np, particles)
        elseif lattice.element_order[i][1] == 14
            pass!(lattice.LongitudinalRLCWakes[lattice.element_order[i][2]], particles6, np, particles)
        elseif lattice.element_order[i][1] == 15
            pass!(lattice.LongitudinalWakes[lattice.element_order[i][2]], particles6, np, particles)
        # else
        #     println("Unknown element type: $lattice.element_order[i][1]")
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
    for i in 1:lattice.nelems
        # println("i = $(i), ADlinepass! element: $(lattice.element_order[i])")
        if i in id
            pass!(elems[c], particles6, np, particles)
            # println("ADlinepass! element: $(elems[c])")
            c += 1
        else
            if lattice.element_order[i][1] == 1
                pass!(lattice.rfcas[lattice.element_order[i][2]], particles6, np, particles)
            elseif lattice.element_order[i][1] == 2
                pass!(lattice.drifts[lattice.element_order[i][2]], particles6, np, particles)
            elseif lattice.element_order[i][1] == 3
                pass!(lattice.markers[lattice.element_order[i][2]], particles6, np, particles)
            elseif lattice.element_order[i][1] == 4
                pass!(lattice.quads[lattice.element_order[i][2]], particles6, np, particles)
            elseif lattice.element_order[i][1] == 5
                pass!(lattice.kquads[lattice.element_order[i][2]], particles6, np, particles)
            elseif lattice.element_order[i][1] == 6
                pass!(lattice.ksexts[lattice.element_order[i][2]], particles6, np, particles)
            elseif lattice.element_order[i][1] == 7
                pass!(lattice.kocts[lattice.element_order[i][2]], particles6, np, particles)
            elseif lattice.element_order[i][1] == 8
                pass!(lattice.sbends[lattice.element_order[i][2]], particles6, np, particles)
            elseif lattice.element_order[i][1] == 9
                pass!(lattice.thinmultipoles[lattice.element_order[i][2]], particles6, np, particles)
            elseif lattice.element_order[i][1] == 10
                pass!(lattice.solenoids[lattice.element_order[i][2]], particles6, np, particles)
            elseif lattice.element_order[i][1] == 11
                pass!(lattice.correctors[lattice.element_order[i][2]], particles6, np, particles)
            elseif lattice.element_order[i][1] == 12
                pass!(lattice.crabcavities[lattice.element_order[i][2]], particles6, np, particles)
            elseif lattice.element_order[i][1] == 13
                pass!(lattice.spacecharges[lattice.element_order[i][2]], particles6, np, particles)
            elseif lattice.element_order[i][1] == 14
                pass!(lattice.LongitudinalRLCWakes[lattice.element_order[i][2]], particles6, np, particles)
            elseif lattice.element_order[i][1] == 15
                pass!(lattice.LongitudinalWakes[lattice.element_order[i][2]], particles6, np, particles)
            # else
            #     println("Unknown element type: $(lattice.element_order[i][1])")
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
    for i in 1:lattice.nelems
        if lattice.element_order[i][1] == 1
            pass_P!(lattice.rfcas[lattice.element_order[i][2]], particles6, np, particles)
        elseif lattice.element_order[i][1] == 2
            pass_P!(lattice.drifts[lattice.element_order[i][2]], particles6, np, particles)
        elseif lattice.element_order[i][1] == 3
            pass_P!(lattice.markers[lattice.element_order[i][2]], particles6, np, particles)
        elseif lattice.element_order[i][1] == 4
            pass_P!(lattice.quads[lattice.element_order[i][2]], particles6, np, particles)
        elseif lattice.element_order[i][1] == 5
            pass_P!(lattice.kquads[lattice.element_order[i][2]], particles6, np, particles)
        elseif lattice.element_order[i][1] == 6
            pass_P!(lattice.ksexts[lattice.element_order[i][2]], particles6, np, particles)
        elseif lattice.element_order[i][1] == 7
            pass_P!(lattice.kocts[lattice.element_order[i][2]], particles6, np, particles)
        elseif lattice.element_order[i][1] == 8
            pass_P!(lattice.sbends[lattice.element_order[i][2]], particles6, np, particles)
        elseif lattice.element_order[i][1] == 9
            pass_P!(lattice.thinmultipoles[lattice.element_order[i][2]], particles6, np, particles)
        elseif lattice.element_order[i][1] == 10
            pass_P!(lattice.solenoids[lattice.element_order[i][2]], particles6, np, particles)
        elseif lattice.element_order[i][1] == 11
            pass_P!(lattice.correctors[lattice.element_order[i][2]], particles6, np, particles)
        elseif lattice.element_order[i][1] == 12
            pass_P!(lattice.crabcavities[lattice.element_order[i][2]], particles6, np, particles)
        elseif lattice.element_order[i][1] == 13
            pass_P!(lattice.spacecharges[lattice.element_order[i][2]], particles6, np, particles)
        elseif lattice.element_order[i][1] == 14
            pass_P!(lattice.LongitudinalRLCWakes[lattice.element_order[i][2]], particles6, np, particles)
        elseif lattice.element_order[i][1] == 15
            pass_P!(lattice.LongitudinalWakes[lattice.element_order[i][2]], particles6, np, particles)
        # else
        #     println("Unknown element type: $lattice.element_order[i][1]")
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
    for i in 1:lattice.nelems
        if i in id
            pass_P!(elems[c], particles6, np, particles)
            c += 1
        else
            if lattice.element_order[i][1] == 1
                pass_P!(lattice.rfcas[lattice.element_order[i][2]], particles6, np, particles)
            elseif lattice.element_order[i][1] == 2
                pass_P!(lattice.drifts[lattice.element_order[i][2]], particles6, np, particles)
            elseif lattice.element_order[i][1] == 3
                pass_P!(lattice.markers[lattice.element_order[i][2]], particles6, np, particles)
            elseif lattice.element_order[i][1] == 4
                pass_P!(lattice.quads[lattice.element_order[i][2]], particles6, np, particles)
            elseif lattice.element_order[i][1] == 5
                pass_P!(lattice.kquads[lattice.element_order[i][2]], particles6, np, particles)
            elseif lattice.element_order[i][1] == 6
                pass_P!(lattice.ksexts[lattice.element_order[i][2]], particles6, np, particles)
            elseif lattice.element_order[i][1] == 7
                pass_P!(lattice.kocts[lattice.element_order[i][2]], particles6, np, particles)
            elseif lattice.element_order[i][1] == 8
                pass_P!(lattice.sbends[lattice.element_order[i][2]], particles6, np, particles)
            elseif lattice.element_order[i][1] == 9
                pass_P!(lattice.thinmultipoles[lattice.element_order[i][2]], particles6, np, particles)
            elseif lattice.element_order[i][1] == 10
                pass_P!(lattice.solenoids[lattice.element_order[i][2]], particles6, np, particles)
            elseif lattice.element_order[i][1] == 11
                pass_P!(lattice.correctors[lattice.element_order[i][2]], particles6, np, particles)
            elseif lattice.element_order[i][1] == 12
                pass_P!(lattice.crabcavities[lattice.element_order[i][2]], particles6, np, particles)
            elseif lattice.element_order[i][1] == 13
                pass_P!(lattice.spacecharges[lattice.element_order[i][2]], particles6, np, particles)
            elseif lattice.element_order[i][1] == 14
                pass_P!(lattice.LongitudinalRLCWakes[lattice.element_order[i][2]], particles6, np, particles)
            elseif lattice.element_order[i][1] == 15
                pass_P!(lattice.LongitudinalWakes[lattice.element_order[i][2]], particles6, np, particles)
            # else
            #     println("Unknown element type: $(lattice.element_order[i][1])")
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
    for i in 1:line.nelems
        if lattice.element_order[i][1] == 1
            pass_TPSA!(line.rfcas[lattice.element_order[i][2]], rin)
        elseif lattice.element_order[i][1] == 2
            pass_TPSA!(line.drifts[lattice.element_order[i][2]], rin)
        elseif lattice.element_order[i][1] == 3
            pass_TPSA!(line.markers[lattice.element_order[i][2]], rin)
        elseif lattice.element_order[i][1] == 4
            pass_TPSA!(line.quads[lattice.element_order[i][2]], rin)
        elseif lattice.element_order[i][1] == 5
            pass_TPSA!(line.kquads[lattice.element_order[i][2]], rin)
        elseif lattice.element_order[i][1] == 6
            pass_TPSA!(line.ksexts[lattice.element_order[i][2]], rin)
        elseif lattice.element_order[i][1] == 7
            pass_TPSA!(line.kocts[lattice.element_order[i][2]], rin)
        elseif lattice.element_order[i][1] == 8
            pass_TPSA!(line.sbends[lattice.element_order[i][2]], rin)
        elseif lattice.element_order[i][1] == 9
            pass_TPSA!(line.thinmultipoles[lattice.element_order[i][2]], rin)
        elseif lattice.element_order[i][1] == 10
            pass_TPSA!(line.solenoids[lattice.element_order[i][2]], rin)
        elseif lattice.element_order[i][1] == 11
            pass_TPSA!(line.correctors[lattice.element_order[i][2]], rin)
        elseif lattice.element_order[i][1] == 12
            pass_TPSA!(line.crabcavities[lattice.element_order[i][2]], rin)
        elseif lattice.element_order[i][1] == 13
            pass_TPSA!(line.spacecharges[lattice.element_order[i][2]], rin)
        elseif lattice.element_order[i][1] == 14
            pass_TPSA!(line.LongitudinalRLCWakes[lattice.element_order[i][2]], rin)
        elseif lattice.element_order[i][1] == 15
            pass_TPSA!(line.LongitudinalWakes[lattice.element_order[i][2]], rin)
        # else
        #     println("Unknown element type: $typ")
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
    for i in 1:line.nelems
        if i in changed_idx
            pass_TPSA!(changed_ele[count], rin)
            count += 1
        else
            if line.element_order[i][1] == 1
                pass_TPSA!(line.rfcas[line.element_order[i][2]], rin)
            elseif line.element_order[i][1] == 2
                pass_TPSA!(line.drifts[line.element_order[i][2]], rin)
            elseif line.element_order[i][1] == 3
                pass_TPSA!(line.markers[line.element_order[i][2]], rin)
            elseif line.element_order[i][1] == 4
                pass_TPSA!(line.quads[line.element_order[i][2]], rin)
            elseif line.element_order[i][1] == 5
                pass_TPSA!(line.kquads[line.element_order[i][2]], rin)
            elseif line.element_order[i][1] == 6
                pass_TPSA!(line.ksexts[line.element_order[i][2]], rin)
            elseif line.element_order[i][1] == 7
                pass_TPSA!(line.kocts[line.element_order[i][2]], rin)
            elseif line.element_order[i][1] == 8
                pass_TPSA!(line.sbends[line.element_order[i][2]], rin)
            elseif line.element_order[i][1] == 9
                pass_TPSA!(line.thinmultipoles[line.element_order[i][2]], rin)
            elseif line.element_order[i][1] == 10
                pass_TPSA!(line.solenoids[line.element_order[i][2]], rin)
            elseif line.element_order[i][1] == 11
                pass_TPSA!(line.correctors[line.element_order[i][2]], rin)
            elseif line.element_order[i][1] == 12
                pass_TPSA!(line.crabcavities[line.element_order[i][2]], rin)
            elseif line.element_order[i][1] == 13
                pass_TPSA!(line.spacecharges[line.element_order[i][2]], rin)
            elseif line.element_order[i][1] == 14
                pass_TPSA!(line.LongitudinalRLCWakes[line.element_order[i][2]], rin)
            elseif line.element_order[i][1] == 15
                pass_TPSA!(line.LongitudinalWakes[line.element_order[i][2]], rin)
            # else
            #     println("Unknown element type: $(line.element_order[i][1])")
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
