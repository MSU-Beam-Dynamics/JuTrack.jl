# ==============================================================================
# This file is part of the TPSA (Truncated Power Series Algebra) Julia package.
#
# Author: Jinyu Wan
# Email: wan@frib.msu.edu
# Version: 1.0
# Created Date: 11-01-2023
# Modified Date: 11-13-2023
include("polymap.jl")
# using StaticArrays

struct CTPS{T, TPS_Dim, Max_TPS_Degree}
    degree::Int
    terms::Int
    map::Vector{T}
    polymap::Ref{PolyMap}
end

global polyMapCache = Dict{Tuple{Int, Int}, PolyMap}()
# Get or create a PolyMap 
function getOrCreatePolyMap(TPS_Dim::Int, Max_TPS_Degree::Int)
    key = (TPS_Dim, Max_TPS_Degree)
    if haskey(polyMapCache, key)
        return polyMapCache[key]
    else
        newPolyMap = PolyMap(TPS_Dim, Max_TPS_Degree)
        polyMapCache[key] = newPolyMap
        return newPolyMap
    end
end

function CTPS(T::Type, TPS_Dim::Int, Max_TPS_Degree::Int) 
    polymap = getOrCreatePolyMap(TPS_Dim, Max_TPS_Degree)
    terms = binomial(TPS_Dim + Max_TPS_Degree, Max_TPS_Degree)
    return CTPS{T, TPS_Dim, Max_TPS_Degree}(0, terms, zeros(T, terms), Ref(polymap)) 
end

function CTPS(a::T, TPS_Dim::Int, Max_TPS_Degree::Int) where T
    polymap = getOrCreatePolyMap(TPS_Dim, Max_TPS_Degree)
    terms = binomial(TPS_Dim + Max_TPS_Degree, Max_TPS_Degree)
    map = zeros(T, terms)
    map[1] = a
    return CTPS{T, TPS_Dim, Max_TPS_Degree}(0, terms, map, Ref(polymap)) 
end

function CTPS(a::T, n::Int, TPS_Dim::Int, Max_TPS_Degree::Int) where T
    if n <= TPS_Dim && n > 0
        terms = binomial(TPS_Dim + Max_TPS_Degree, Max_TPS_Degree)
        map  = zeros(T, terms)
        # map = zeros(SVector{terms, T})
        map[n+1] = one(T)
        map[1] = a
        polymap = getOrCreatePolyMap(TPS_Dim, Max_TPS_Degree)
        return CTPS{T, TPS_Dim, Max_TPS_Degree}(1, terms, map, Ref(polymap)) 
    else
        error("Num of var out of range in CTPS")
    end
end

function CTPS(M::CTPS{T, TPS_Dim, Max_TPS_Degree}) where {T, TPS_Dim, Max_TPS_Degree}
    map = zeros(T, length(M.map))
    for i in eachindex(map)
        map[i] = M.map[i]
    end
    return CTPS{T, TPS_Dim, Max_TPS_Degree}(Max_TPS_Degree, length(map), map, M.polymap)
end

function cst(ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}) where {T, TPS_Dim, Max_TPS_Degree}
    a0 = ctps.map[1]
    return a0
end

function findindex(ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}, indexmap::Vector{Int}) where {T, TPS_Dim, Max_TPS_Degree}
    # find the index of the indexmap in the map vector
    # indexmap is a vector of length TPS_Dim + 1, e.g. [0, 1, 1] for x1^1 * x2^1
    dim = TPS_Dim
    if length(indexmap) == dim
        total = sum(indexmap)
        newindexmap = [total; indexmap]
        return findindex(ctps, newindexmap)
    end
    if length(indexmap) != (dim + 1)
        error("Index map does not have correct length")
    end
    sum = zeros(Int, length(indexmap))
    for i in 1:length(indexmap)
        sum[i] = indexmap[i]
    end
    # sum = copy(indexmap)
    # sum_buffer = Zygote.Buffer(sum) # Buffer for Zygote
    # for i in 1:length(sum_buffer)
    #     sum_buffer[i] = sum[i]
    # end
    for i in 2:dim+1
        if indexmap[i] < 0
            error("The index map has invalid component")
        end
        sum[i] = sum[i-1] - indexmap[i]
    end
    # sum = copy(sum_buffer)
    result = Int(1)
    for i in dim:-1:1
        if sum[dim - i + 1] == 0
            break
        end
        result += binomial(sum[dim - i + 1] - 1 + i, i)
    end
    return result
end

# function findpower(ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}, n::Int) where {T, TPS_Dim, Max_TPS_Degree}
#     if n < ctps.terms
#         return getindexmap(ctps.polymap[], n)
#     else
#         error("The index is out of range")
#     end
# end

# function redegree!(ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}, degree::Int) where {T, TPS_Dim, Max_TPS_Degree}
#     ctps.degree = min(degree, Max_TPS_Degree)
#     ctps.terms = binomial(TPS_Dim + ctps.degree, ctps.degree)
#     new_map = [i <= length(ctps.map) ? ctps.map[i] : zero(T) for i in 1:ctps.terms]
#     ctps.map = new_map
# end
# function redegree(ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}, degree::Int) where {T, TPS_Dim, Max_TPS_Degree}
#     degree = min(degree, Max_TPS_Degree)
#     terms = binomial(TPS_Dim + degree, degree)
#     new_map = zeros(T, terms)
#     for i in 1:ctps.terms
#         new_map[i] = ctps.map[i]
#     end
#     # new_map = [i <= length(ctps.map) ? ctps.map[i] : 0.0 for i in 1:terms]
#     # polymap = getOrCreatePolyMap(TPS_Dim, Max_TPS_Degree)
#     ctps_new = CTPS{T, TPS_Dim, Max_TPS_Degree}(degree, terms, new_map, ctps.polymap)
#     return ctps_new
# end
# function redegree(ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}, degree::Int) where {T, TPS_Dim, Max_TPS_Degree}
#     degree = min(degree, Max_TPS_Degree)
#     terms = binomial(TPS_Dim + degree, degree)
#     new_map = zeros(T, terms)
#     new_map_buffer = Zygote.Buffer(new_map)
#     for i in 1:ctps.terms
#         new_map_buffer[i] = ctps.map[i]
#     end
#     for i in ctps.terms+1:terms
#         new_map_buffer[i] = zero(T)
#     end
#     new_map = copy(new_map_buffer)
#     ctps_new = CTPS{T, TPS_Dim, Max_TPS_Degree}(degree, terms, new_map, PolyMap(TPS_Dim, Max_TPS_Degree))
#     return ctps_new
# end
function assign!(ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}, a::T, n_var::Int) where {T, TPS_Dim, Max_TPS_Degree}
    if n_var <= TPS_Dim && n_var > 0
        map = ctps.map
        map[n_var+1] = one(T)
        map[1] = a
        return nothing
    else
        error("Num of var out of range in CTPS")
    end
end

function assign!(ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}, a::T) where {T, TPS_Dim, Max_TPS_Degree}
    ctps.map[1] = a
    return nothing
end

function reassign!(ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}, a::T, n_var::Int) where {T, TPS_Dim, Max_TPS_Degree}
    if n_var <= TPS_Dim && n_var > 0
        map = ctps.map
        for i in eachindex(map)
            map[i] = zero(T)
        end
        map[n_var+1] = one(T)
        map[1] = a
        return nothing
    else
        error("Num of var out of range in CTPS")
    end
end

# function element(ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}, ind::Int) where {T, TPS_Dim, Max_TPS_Degree}
#     if ind < 1 || ind > ctps.terms
#         error("Element index out of range in CTPS")
#     end
#     return ctps.map[ind]
# end

function element(ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}, ind::Vector{Int}) where {T, TPS_Dim, Max_TPS_Degree}
    result = findindex(ctps, ind)
    return ctps.map[result]
end

# function evaluate(ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}, inivalue::Vector{U}) where {T, TPS_Dim, Max_TPS_Degree, U}
#     if length(inivalue) != TPS_Dim
#         error("Inconsistent dimension to evaluate CTPS")
#     end
#     sum = U(ctps.map[1])
#     for i in 2:ctps.terms
#         temp = getindexmap(ctps.polymap[], i)
#         product = U(1)
#         for j in 1:TPS_Dim
#             dimpower = U(inivalue[j])^temp[j+1]
#             product *= dimpower
#         end
#         sum += product * U(ctps.map[i])
#     end
#     return sum
# end

# function derivative(ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}, ndim::Int, order::Int) where {T, TPS_Dim, Max_TPS_Degree}
#     if order <= 0
#         return CTPS(ctps)
#     end
#     if order > ctps.degree
#         return CTPS(T, TPS_Dim, Max_TPS_Degree)
#     end
#     if ndim <= TPS_Dim && ndim > 0
#         derivative_ctps = CTPS(T, TPS_Dim, Max_TPS_Degree)
#         derivative_ctps = redegree(derivative_ctps, ctps.degree - order)
#         for i in 2:ctps.terms
#             temp = getindexmap(ctps.polymap[], i)
#             if temp[ndim + 1] >= order
#                 thisdim = temp[ndim + 1]
#                 # buf = Zygote.Buffer(temp)  # Buffer for Zygote
#                 # for j in 1:length(buf)
#                 #     buf[j] = temp[j]
#                 # end
#                 temp[ndim + 1] -= order
#                 temp[1] -= order
#                 # new_temp = copy(buf)  
#                 index = findindex(ctps, temp)
#                 derivative_ctps_map_buffer = copy(derivative_ctps.map)
#                 # for j in 1:length(derivative_ctps_map_buffer)
#                 #     derivative_ctps_map_buffer[j] = derivative_ctps.map[j]
#                 # end
#                 derivative_ctps_map_buffer[index] = factorial(new_temp[ndim + 1] + order) / factorial(new_temp[ndim + 1]) * ctps.map[i]
#                 derivative_ctps = CTPS{T, TPS_Dim, Max_TPS_Degree}(derivative_ctps.degree, derivative_ctps.terms, derivative_ctps_map_buffer, derivative_ctps.polymap)
#                 # derivative_ctps.map = copy(derivative_ctps_map_buffer)
#                 # derivative_ctps.map[index] = factorial(new_temp[ndim + 1] + order) / factorial(new_temp[ndim + 1]) * ctps.map[i]
#             end
#         end
#         return derivative_ctps
#     else
#         error("The dimension is out of range")
#     end
# end

# function integrate(ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}, ndim::Int, a0::T) where {T, TPS_Dim, Max_TPS_Degree}
#     if ndim <= TPS_Dim && ndim > 0
#         temp = CTPS(a0, TPS_Dim, Max_TPS_Degree)
#         temp = redegree(temp, ctps.degree + 1)
#         map_buffer = copy(temp.map)
#         # map_buffer = Zygote.Buffer(temp.map)
#         # for i in 1:temp.terms
#         #     map_buffer[i] = temp.map[i]
#         # end
#         for i in 1:ctps.terms
#             indexlist = getindexmap(ctps.polymap[], i)
#             thisdim = indexlist[ndim+1]
#             # indexlist_buffer = Zygote.Buffer(indexlist) 
#             # for j in 1:length(indexlist_buffer)
#             #     indexlist_buffer[j] = indexlist[j]
#             # end
#             indexlist[ndim+1] += 1
#             indexlist[1] += 1
#             # indexlist = copy(indexlist_buffer)
#             new_i = findindex(ctps, indexlist)
#             map_buffer[new_i] = ctps.map[i] / (thisdim + 1)
#         end
#         temp = CTPS{T, TPS_Dim, Max_TPS_Degree}(temp.degree, temp.terms, map_buffer, temp.polymap)
#         # temp.map = copy(map_buffer)
#         return temp
#     else
#         error("Inconsistent dimension to integrate")
#     end
# end

# Overloaded operations
import Base: +, -, *, /, sin, cos, tan, sinh, cosh, asin, acos, sqrt, ^, inv, exp, log

# +
function +(ctps1::CTPS{T, TPS_Dim, Max_TPS_Degree}, ctps2::CTPS{T, TPS_Dim, Max_TPS_Degree}) where {T, TPS_Dim, Max_TPS_Degree}
    ctps_new = CTPS(ctps1)
    for i in eachindex(ctps_new.map)
        ctps_new.map[i] += ctps2.map[i]
    end
    return ctps_new
    # ctps_map_buffer = zeros(T, ctps1.terms)
    # for i in 1:ctps2.terms
    #     ctps_map_buffer[i] = ctps1.map[i] + ctps2.map[i]
    # end
    # return CTPS{T, TPS_Dim, Max_TPS_Degree}(ctps1.degree, ctps1.terms, ctps_map_buffer, ctps1.polymap)
end



function +(ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}, a::T) where {T, TPS_Dim, Max_TPS_Degree}
    ctps_new = CTPS(ctps)
    ctps_new.map[1] += a
    return ctps_new
end
function +(a::T, ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}) where {T, TPS_Dim, Max_TPS_Degree}
    return ctps + a
end
# -
function -(ctps1::CTPS{T, TPS_Dim, Max_TPS_Degree}, ctps2::CTPS{T, TPS_Dim, Max_TPS_Degree}) where {T, TPS_Dim, Max_TPS_Degree}
    ctps_new = CTPS(ctps1)
    for i in eachindex(ctps_new.map)
        ctps_new.map[i] -= ctps2.map[i]
    end
    return ctps_new
    # ctps_map_buffer = zeros(T, ctps1.terms)
    # for i in 1:ctps2.terms
    #     ctps_map_buffer[i] = ctps1.map[i] - ctps2.map[i]
    # end
    # return CTPS{T, TPS_Dim, Max_TPS_Degree}(ctps1.degree, ctps1.terms, ctps_map_buffer, ctps1.polymap)
end
function -(ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}, a::T) where {T, TPS_Dim, Max_TPS_Degree}
    ctps_new = CTPS(ctps)
    ctps_new.map[1] -= a
    return ctps_new
    # ctps_map_buffer = zeros(T, ctps.terms)
    # ctps_map_buffer[1] = ctps.map[1] - a
    # for i in 2:ctps.terms
    #     ctps_map_buffer[i] = ctps.map[i]
    # end
    # ctps_new = CTPS{T, TPS_Dim, Max_TPS_Degree}(ctps.degree, ctps.terms, ctps_map_buffer, ctps.polymap)
    # return ctps_new
end
function -(a::T, ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}) where {T, TPS_Dim, Max_TPS_Degree}
    ctps_new = CTPS(ctps)
    ctps_new.map[1] = a - ctps_new.map[1]
    for i in eachindex(ctps_new.map)[2:end]
        ctps_new.map[i] = -ctps_new.map[i]
    end
    return ctps_new
    # ctps_map_buffer = zeros(T, ctps.terms)
    # ctps_map_buffer[1] = a - ctps.map[1]
    # for i in 2:ctps.terms
    #     ctps_map_buffer[i] = -ctps.map[i]
    # end
    # ctps_new = CTPS{T, TPS_Dim, Max_TPS_Degree}(ctps.degree, ctps.terms, ctps_map_buffer, ctps.polymap)
    # return ctps_new
end
function -(ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}) where {T, TPS_Dim, Max_TPS_Degree}
    ctps_new = CTPS(ctps)
    for i in eachindex(ctps_new.map)
        ctps_new.map[i] = -ctps_new.map[i]
    end
    # ctps_map_buffer = zeros(T, ctps.terms)
    # for i in 1:ctps.terms
    #     ctps_map_buffer[i] = -ctps.map[i]
    # end
    # ctps_new = CTPS{T, TPS_Dim, Max_TPS_Degree}(ctps.degree, ctps.terms, ctps_map_buffer, ctps.polymap)
    return ctps_new
end

# *
function *(ctps1::CTPS{T, TPS_Dim, Max_TPS_Degree}, ctps2::CTPS{T, TPS_Dim, Max_TPS_Degree}) where {T, TPS_Dim, Max_TPS_Degree}
    # if ctps2.degree == 0
    #     ctps_map_buffer = zeros(T, ctps1.terms)
    #     for i in 1:ctps1.terms
    #         ctps_map_buffer[i] = ctps1.map[i] * ctps2.map[1]
    #     end
    #     # map = ctps1.map * ctps2.map[1]
    #     ctps_new = CTPS{T, TPS_Dim, Max_TPS_Degree}(ctps1.degree, ctps1.terms, ctps_map_buffer, ctps1.polymap)
    #     return ctps_new
    # elseif ctps1.degree == 0
    #     ctps_map_buffer = zeros(T, ctps2.terms)
    #     for i in 1:ctps2.terms
    #         ctps_map_buffer[i] = ctps2.map[i] * ctps1.map[1]
    #     end
    #     # map = ctps2.map * ctps1.map[1]
    #     ctps_new = CTPS{T, TPS_Dim, Max_TPS_Degree}(ctps2.degree, ctps2.terms, ctps_map_buffer, ctps2.polymap)
    #     return ctps_new
    # end

    ctps_new = CTPS(ctps1)
    # ctps_new = redegree(ctps1, ctps1.degree + ctps2.degree)
    ctps_map_buffer = ctps_new.map
    # ctps_map_buffer = zeros(T, length(ctps_new.map))
    for i in eachindex(ctps_map_buffer)
        ctps_map_buffer[i] = 0.0
    end

    for i in 1:length(ctps1.map)
        if ctps1.map[i] == 0
            continue
        end
        temp1 = getindexmap(ctps1.polymap[], i)
        j_max = min(length(ctps2.map), binomial(TPS_Dim + Max_TPS_Degree - temp1[1], TPS_Dim))
        for j in 1:j_max
            if ctps2.map[j] == 0
                continue
            end 
            temp2 = getindexmap(ctps2.polymap[], j)
            # no array operation
            temp = zeros(Int, TPS_Dim + 1)
            for k in 1:TPS_Dim+1
                temp[k] = temp1[k] + temp2[k]
            end
            # temp = temp1 + temp2
            index = findindex(ctps_new, temp)
            ctps_map_buffer[index] += ctps1.map[i] * ctps2.map[j]
        end
    end
    return ctps_new
    # return CTPS{T, TPS_Dim, Max_TPS_Degree}(ctps_new.degree, ctps_new.terms, ctps_map_buffer, ctps_new.polymap)
end
function *(ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}, a::T) where {T, TPS_Dim, Max_TPS_Degree}
    ctps_new = CTPS(ctps)
    for i in eachindex(ctps_new.map)
        ctps_new.map[i] *= a
    end
    # map = zeros(T, length(ctps.map))
    # for i in eachindex(map)
    #     map[i] = ctps.map[i] * a
    # end
    # ctps_new = CTPS{T, TPS_Dim, Max_TPS_Degree}(ctps.degree, ctps.terms, map, ctps.polymap)
    return ctps_new
end
function *(a::T, ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}) where {T, TPS_Dim, Max_TPS_Degree}
    return ctps * a
end

# /
function inv(ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}) where {T, TPS_Dim, Max_TPS_Degree}
    if cst(ctps) == zero(T)
        error("Divide by zero in CTPS")
    end
    temp = CTPS(ctps) - cst(ctps)
    # sum = CTPS(zero(T), TPS_Dim, Max_TPS_Degree)
    term_by_oder = CTPS(1.0 / ctps.map[1], TPS_Dim, Max_TPS_Degree)

    # sum = sum + term_by_oder
    sum = CTPS(term_by_oder)
    for i in 1:Max_TPS_Degree
        term_by_oder *= -temp/cst(ctps)
        sum += term_by_oder
    end
    return sum
end

function /(ctps1::CTPS{T, TPS_Dim, Max_TPS_Degree}, ctps2::CTPS{T, TPS_Dim, Max_TPS_Degree}) where {T, TPS_Dim, Max_TPS_Degree}
    if cst(ctps2) == zero(T)
        error("Divide by zero in CTPS")
    end

    # if ctps2.degree == 0
    #     map = ctps1.map
    #     for i in eachindex(map)
    #         map[i] /= ctps2.map[1]
    #     end
    #     return ctps1
    # end

    return ctps1 * inv(ctps2)
end
function /(ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}, a::T) where {T, TPS_Dim, Max_TPS_Degree}
    if a == zero(T)
        error("Divide by zero in CTPS")
    end
    ctps_new = CTPS(ctps)
    for i in eachindex(ctps_new.map)
        ctps_new.map[i] /= a
    end
    # ctps_map_buffer = zeros(T, ctps.terms)
    # for i in eachindex(ctps_map_buffer)
    #     ctps_map_buffer[i] = ctps.map[i] / a
    # end
    # ctps_new = CTPS{T, TPS_Dim, Max_TPS_Degree}(ctps.degree, ctps.terms, ctps_map_buffer, ctps.polymap)
    return ctps_new
end
function /(a::T, ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}) where {T, TPS_Dim, Max_TPS_Degree}
    if cst(ctps) == zero(T)
        error("Divide by zero in CTPS")
    end
    return a * inv(ctps)
end

# exponential
function exp(ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}) where {T, TPS_Dim, Max_TPS_Degree}
    if cst(ctps) == zero(T)
        return CTPS(one(T), TPS_Dim, Max_TPS_Degree)
    end
    temp = CTPS(ctps)
    temp -= cst(temp)
    term_by_oder = CTPS(one(T), TPS_Dim, Max_TPS_Degree)
    sum = CTPS(one(T), TPS_Dim, Max_TPS_Degree)
    for i in 1:Max_TPS_Degree
        index = 1.0 / factorial(i)
        term_by_oder = term_by_oder * temp
        sum = sum + (term_by_oder * T(index))
    end
    sum = sum * T(Base.exp(cst(ctps)))
    return sum
end

# logarithm
function log(ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}) where {T, TPS_Dim, Max_TPS_Degree}
    if cst(ctps) == zero(T)
        error("Log of zero in CTPS")
    end
    temp = CTPS(ctps)
    temp -= cst(temp)
    term_by_oder = temp / cst(ctps)
    sum = CTPS(zero(T), TPS_Dim, Max_TPS_Degree) + term_by_oder

    for i in 2:Max_TPS_Degree
        term_by_oder = term_by_oder * (-temp / cst(ctps))
        sum = sum + (term_by_oder / T(i))
    end
    sum = sum + Base.log(cst(ctps))
    return sum
end

# square root
function sqrt(ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}) where {T, TPS_Dim, Max_TPS_Degree}
    if cst(ctps) < zero(T)
        error("Square root of negative number in CTPS")
    end
    a0 = sqrt(cst(ctps))
    temp = CTPS(ctps) - cst(ctps)
    # temp = temp 
    term_by_oder = temp / a0
    sum = term_by_oder/2.0

    for i in 2:Max_TPS_Degree
        index = 1.0 * doublefactorial(2 * i - 3) / doublefactorial(2 * i)
        term_by_oder = (-temp) * term_by_oder / cst(ctps)
        sum += term_by_oder * index
    end
    sum = sum + a0
    return sum
end

# power
function pow(ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}, b::Int) where {T, TPS_Dim, Max_TPS_Degree}
    if b == 1
        return ctps
    elseif b == 0
        return CTPS(one(T), TPS_Dim, Max_TPS_Degree)
    end

    temp = CTPS(ctps)
    index = T(b)
    if cst(ctps) == zero(T)
        # if mod(b, 1.0) == 0 && b > 1
        if b > 1
            sum = CTPS(ctps)
            for i in 2:b
                sum = sum * ctps
            end
            return sum
        else
            error("Divide by zero, in CTPS::pow")
        end
    end
    temp = temp - cst(temp)
    term_by_oder = CTPS(one(T), TPS_Dim, Max_TPS_Degree)
    factor = cst(ctps) ^ b
    sum = CTPS(factor, TPS_Dim, Max_TPS_Degree)

    for i in 1:Max_TPS_Degree
        factor = factor / cst(ctps) * index / i
        index -= 1.0
        term_by_oder = term_by_oder * temp
        sum = sum + (term_by_oder * factor)
        if index == 0.0
            break
        end
    end
    return sum
end
function ^(ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}, b::Int) where {T, TPS_Dim, Max_TPS_Degree}
    return pow(ctps, b)
end

# sin
function sin(ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}) where {T, TPS_Dim, Max_TPS_Degree}
    temp = CTPS(ctps)
    a0 = cst(ctps)
    sin_a0 = sin(a0)
    cos_a0 = cos(a0)
    temp = temp - a0
    term_by_oder = CTPS(one(T), TPS_Dim, Max_TPS_Degree)
    sum = CTPS(zero(T), TPS_Dim, Max_TPS_Degree)

    is_odd_iteration = 1  
    for i in 1:Max_TPS_Degree
        if is_odd_iteration == 1
            index = cos_a0 * (-1) ^ ((i - 1) / 2) / doublefactorial(i)
        else
            index = sin_a0 * (-1) ^ (i / 2) / doublefactorial(i)
        end
        term_by_oder = term_by_oder * temp
        sum = sum + (term_by_oder * index)

        is_odd_iteration = -is_odd_iteration 
    end

    # for i in 1:Max_TPS_Degree
    #     if mod(i, 2) == 1
    #         index = cos_a0 * (-1) ^ ((i - 1) / 2) / doublefactorial(i)
    #     else
    #         index = sin_a0 * (-1) ^ (i / 2) / doublefactorial(i)
    #     end
    #     term_by_oder = term_by_oder * temp
    #     sum = sum + (term_by_oder * index)
    # end
    sum = sum + sin_a0
    return sum
end

# cos
function cos(ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}) where {T, TPS_Dim, Max_TPS_Degree}
    temp = CTPS(ctps)
    a0 = cst(ctps)
    sin_a0 = sin(a0)
    cos_a0 = cos(a0)
    temp = temp - a0
    term_by_oder = CTPS(one(T), TPS_Dim, Max_TPS_Degree)
    sum = CTPS(zero(T), TPS_Dim, Max_TPS_Degree)

    is_odd_iteration = 1  
    for i in 1:Max_TPS_Degree
        if is_odd_iteration == 1
            index = sin_a0 * (-1) ^ ((i + 1) / 2) / doublefactorial(i)
        else
            index = cos_a0 * (-1) ^ (i / 2) / doublefactorial(i)
        end
        term_by_oder = term_by_oder * temp
        sum = sum + (term_by_oder * index)
    
        is_odd_iteration = -is_odd_iteration  
    end
    
    # for i in 1:Max_TPS_Degree
    #     if mod(i, 2) == 1
    #         index = sin_a0 * (-1) ^ ((i + 1) / 2) / factorial(i)
    #     else
    #         index = cos_a0 * (-1) ^ (i / 2) / factorial(i)
    #     end
    #     term_by_oder = term_by_oder * temp
    #     sum = sum + (term_by_oder * index)
    # end
    sum = sum + cos_a0
    return sum
end

# arcsin
function asin(ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}) where {T, TPS_Dim, Max_TPS_Degree}
    temp = CTPS(ctps)
    a0 = cst(ctps)
    arcsin_a0 = asin(a0)
    cos_y0 = sqrt(one(T) - a0 ^ 2)
    temp = temp - a0
    temp1 = CTPS(temp)
    for i in 1:Max_TPS_Degree + 10
        temp1 = temp1 - sin(temp1) + (temp + a0*(1.0-cos(temp1)))/cos_y0
    end
    return temp1 + arcsin_a0
end

# arccos
function acos(ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}) where {T, TPS_Dim, Max_TPS_Degree}
    return T(pi / 2) - asin(ctps)
end

# tangent
function tan(ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}) where {T, TPS_Dim, Max_TPS_Degree}
    return sin(ctps) / cos(ctps)
end

# hyperbolic sin
function sinh(ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}) where {T, TPS_Dim, Max_TPS_Degree}
    temp = CTPS(ctps)
    a0 = cst(ctps)
    sinh_a0 = sinh(a0)
    cosh_a0 = cosh(a0)
    temp = temp - a0
    term_by_oder = CTPS(one(T), TPS_Dim, Max_TPS_Degree)
    sum = CTPS(zero(T), TPS_Dim, Max_TPS_Degree)

    is_odd_iteration = 1  
    for i in 1:Max_TPS_Degree
        if is_odd_iteration == 1
            index = cosh_a0 / doublefactorial(i)
        else
            index = sinh_a0 / doublefactorial(i)
        end
        term_by_oder = term_by_oder * temp
        sum = sum + (term_by_oder * index)
        is_odd_iteration = -is_odd_iteration
    end
    sum = sum + sinh_a0
    return sum
end

# hyperbolic cos
function cosh(ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}) where {T, TPS_Dim, Max_TPS_Degree}
    temp = CTPS(ctps)
    a0 = cst(ctps)
    sinh_a0 = sinh(a0)
    cosh_a0 = cosh(a0)
    temp = temp - a0
    term_by_oder = CTPS(one(T), TPS_Dim, Max_TPS_Degree)
    sum = CTPS(zero(T), TPS_Dim, Max_TPS_Degree)

    is_odd_iteration = 1
    for i in 1:Max_TPS_Degree
        if is_odd_iteration == 1
            index = sinh_a0 / doublefactorial(i)
        else
            index = cosh_a0 / doublefactorial(i)
        end
        term_by_oder = term_by_oder * temp
        sum = sum + (term_by_oder * index)
        is_odd_iteration = -is_odd_iteration
    end
    sum = sum + cosh_a0
    return sum
end