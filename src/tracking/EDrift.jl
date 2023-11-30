function EDrift(r::Array{T}, L::T) where T
    x_new = r[1] + L * r[2]
    y_new = r[3] + L * r[4]
    s_new = r[5] + L * sqrt(1 + r[2]^2 + r[4]^2)
    return [x_new, r[2], y_new, r[4], s_new, r[6]]
end