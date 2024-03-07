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