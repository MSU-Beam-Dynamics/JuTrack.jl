function total_length(ring)
    length = 0.0
    for i in eachindex(ring)
        length += ring[i].len
    end
    return length
end

function spos(ring)
    pos = zeros(length(ring))
    for i in eachindex(ring)
        pos[i] = total_length(ring[1:i])
    end
    return pos
end

function findelem(ring, field::Symbol, value)
    # Example: findelem(ring, :name, "QF")
    ele_index = []
    for i in eachindex(ring)
        if field in fieldnames(typeof(ring[i])) && getfield(ring[i], field) == value
            push!(ele_index, i)
        end
    end
    return ele_index
end

function findelem(ring, type::Type)
    # Example: findelem(ring, DRFIT)
    ele_index = []
    for i in eachindex(ring)
        if typeof(ring[i]) == type
            push!(ele_index, i)
        end
    end
    return ele_index
end

function use_exact_drift(flag)
    if flag == 1
        global use_exact_Hamiltonian = 1
    else
        global use_exact_Hamiltonian = 0
    end
end