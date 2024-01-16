# ==============================================================================
# This file is part of the TPSA (Truncated Power Series Algebra) Julia package.
#
# Author: Jinyu Wan
# Email: wan@frib.msu.edu
# Version: 1.0
# Created Date: 11-01-2023
# Modified Date: 11-06-2023

include("mathfunc.jl")
# using Zygote

struct PolyMap
    dim::Int
    max_order::Int
    map::Matrix{Int}

    function PolyMap(dim::Int, order::Int)
        new(dim, order, setindexmap(dim, order))
    end
end
function Base.copy(pm::PolyMap)
    new_pm = PolyMap(pm.dim, pm.max_order)
    return new_pm
end

function decomposite(n::Int, dim::Int)
    result = zeros(dim + 1)
    itemp = n + 1
    for i in dim:-1:1
        k = i - 1
        while binomial(k, i) < itemp
            k += 1
        end
        itemp -= binomial(k - 1, i)
        result[dim - i + 1] = k - i  
    end
    for i in dim:-1:2  
        result[i] = result[i - 1] - result[i]
    end
    return result
end


function setindexmap(dim::Int, max_order::Int)
    totallength = binomial(max_order + dim, dim)
    # map = [decomposite(i, dim) for i in 0:totallength-1]
    # create map as a matrix
    map = zeros(Int, totallength, dim + 1)
    for i in 0:totallength-1
        map[i + 1, :] = decomposite(i, dim)
    end
    return map
end

# function getindexmap(p::PolyMap, i::Int)
#     if i < 1 || i > length(p.map)
#         error("index out of range")
#     end
#     return p.map[i]
# end
function getindexmap(p::PolyMap, i::Int)
    if i < 1 || i > length(p.map)
        error("index out of range")
    end
    return p.map[i,:]
end
# z = PolyMap(6, 2)
# println(getindexmap(z, 1))
