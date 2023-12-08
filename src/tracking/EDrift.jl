# function EDrift(r::Array{T}, L::T) where T
#     x_new = r[1] + L * r[2]
#     y_new = r[3] + L * r[4]
#     s_new = r[5] + L * sqrt(1 + r[2]^2 + r[4]^2)
#     return [x_new, r[2], y_new, r[4], s_new, r[6]]
# end

function EDrift(part, np, length)
    return [
        [
            coord[1] + coord[2] * length,
            coord[2],
            coord[3] + coord[4] * length,
            coord[4],
            coord[5] + length * sqrt(1 + coord[2]^2 + coord[4]^2),
            coord[6]
        ]
        for coord in part
    ]
end