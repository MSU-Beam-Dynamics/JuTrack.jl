# include("../lattice/canonical_elements.jl")

# function EDrift_pass(part, np, OO)
#     len = OO.len
#     for i in 1:np
#         part[i,1] += part[i,2] * len
#         part[i,3] += part[i,4] * len
#         part[i,5] += len * sqrt(1.0 + part[i,2]^2 + part[i,4]^2)
#     end
# end

function pass!(part::Matrix{Float64}, OO::EDRIFT, Po::Float64, sigmaDelta2::Float64)
    len = OO.len
    for i in 1:length(part[:,1])
        part[i,1] += part[i,2] * len
        part[i,3] += part[i,4] * len
        part[i,5] += len * sqrt(1.0 + part[i,2]^2 + part[i,4]^2)
    end
end

# include("../TPSA_Enzyme/TPSA_fixedmap.jl")
# function pass_TPSA(x::CTPS, xp::CTPS, y::CTPS, yp::CTPS, z::CTPS, delta::CTPS, 
#                     OO::EDRIFT, Po::Float64, sigmaDelta2::Float64)
#     length = OO.len
#     x += xp * length
#     y += yp * length
#     z += length * sqrt(1.0 + xp^2 + yp^2)
#     return x, xp, y, yp, z, delta
# end



# using Enzyme
# function ff(x, y)
#     k = x[1]
#     particle = [Float64[0.001, 0.0001, 0.0005, 0.0002, 0.0, 0.0], Float64[0.001, 0.0, 0.0, 0.0, 1.0, 0.0]]
#     EDrift(particle, 2, k)    
#     y[1] = particle[1][1]
#     return nothing
# end

# function test_E(k)
#     x = [k, 0.0]
#     y = [0.0]
#     bx = [0.0, 0.0]
#     by = [1.0] 
#     Enzyme.autodiff(Reverse, ff, Duplicated(x, bx), Duplicated(y, by))
# end
# k = 1.0

# x = [k, 0.0]
# y = [0.0]
# bx = [0.0, 0.0]
# by = [1.0] 
# ff(x, y)
# @time ff(x, y)
# test_E(k)
# @time test_E(k)
