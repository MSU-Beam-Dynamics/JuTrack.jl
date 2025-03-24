# ==============================================================================
# This file is part of the TPSA (Truncated Power Series Algebra) Julia package.
# Translated from Yue Hao's C++ code
# Author: Jinyu Wan
# Email: wan@frib.msu.edu
# Created Date: 11-01-2023
# Modified Date: 03-21-2025
include("polymap.jl")

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

"""
    CTPS(T::Type, TPS_Dim::Int, Max_TPS_Degree::Int)

Create a truncated power series (TPS) object with type `T`, dimension `TPS_Dim`, and maximum degree `Max_TPS_Degree`.
"""
function CTPS(T::Type, TPS_Dim::Int, Max_TPS_Degree::Int) 
    polymap = getOrCreatePolyMap(TPS_Dim, Max_TPS_Degree)
    terms = binomial(TPS_Dim + Max_TPS_Degree, Max_TPS_Degree)
    return CTPS{T, TPS_Dim, Max_TPS_Degree}(0, terms, zeros(T, terms), Ref(polymap)) 
end

"""
    CTPS(a::T, TPS_Dim::Int, Max_TPS_Degree::Int)

Create a truncated power series (TPS) object with type `T`, dimension `TPS_Dim`, and maximum degree `Max_TPS_Degree`, and set the constant term to `a`.
"""
function CTPS(a::T, TPS_Dim::Int, Max_TPS_Degree::Int) where T
    polymap = getOrCreatePolyMap(TPS_Dim, Max_TPS_Degree)
    terms = binomial(TPS_Dim + Max_TPS_Degree, Max_TPS_Degree)
    map = zeros(T, terms)
    map[1] = a
    return CTPS{T, TPS_Dim, Max_TPS_Degree}(0, terms, map, Ref(polymap)) 
end

"""
    CTPS(a::T, n::Int, TPS_Dim::Int, Max_TPS_Degree::Int)

Create a truncated power series (TPS) object with type `T`, dimension `TPS_Dim`, and maximum degree `Max_TPS_Degree`, and set the `n`-th variable to `a`.
"""
function CTPS(a::T, n::Int, TPS_Dim::Int, Max_TPS_Degree::Int) where T
    if n <= TPS_Dim && n > 0
        terms = binomial(TPS_Dim + Max_TPS_Degree, Max_TPS_Degree)
        map  = zeros(T, terms)
        map[n+1] = one(T)
        map[1] = a
        polymap = getOrCreatePolyMap(TPS_Dim, Max_TPS_Degree)
        return CTPS{T, TPS_Dim, Max_TPS_Degree}(1, terms, map, Ref(polymap)) 
    else
        error("Num of var out of range in CTPS")
    end
end

"""
    CTPS(M::CTPS{T, TPS_Dim, Max_TPS_Degree})

Copy a truncated power series (TPS) object `M`.
"""
function CTPS(M::CTPS{T, TPS_Dim, Max_TPS_Degree}) where {T, TPS_Dim, Max_TPS_Degree}
    map = zeros(T, length(M.map))
    for i in eachindex(map)
        map[i] = M.map[i]
    end
    return CTPS{T, TPS_Dim, Max_TPS_Degree}(M.degree, length(map), map, M.polymap)
end

"""
    cst(ctps::CTPS{T, TPS_Dim, Max_TPS_Degree})

Return the constant term of a truncated power series (TPS) object `ctps`.
"""
function cst(ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}) where {T, TPS_Dim, Max_TPS_Degree}
    a0 = ctps.map[1]
    return a0
end

"""
    findindex(ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}, indexmap::Vector{Int})

Find the index of the indexmap in the map vector.
"""
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
    SUM = zeros(Int, length(indexmap))
    for i in 1:length(indexmap)
        SUM[i] = indexmap[i]
    end

    for i in 2:dim+1
        if indexmap[i] < 0
            error("The index map has invalid component")
        end
        SUM[i] = SUM[i-1] - indexmap[i]
    end
    # sum = copy(sum_buffer)
    result = Int(1)
    for i in dim:-1:1
        if SUM[dim - i + 1] == 0
            break
        end
        result += binomial(SUM[dim - i + 1] - 1 + i, i)
    end
    return result
end

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
end

function +(ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}, a::T) where {T, TPS_Dim, Max_TPS_Degree}
    ctps_new = CTPS(ctps)
    ctps_new.map[1] += a
    return ctps_new
end
function +(a::T, ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}) where {T, TPS_Dim, Max_TPS_Degree}
    return ctps + a
end

# Addition with complex and float64
function +(ctps::CTPS{ComplexF64, TPS_Dim, Max_TPS_Degree}, a::Float64) where {TPS_Dim, Max_TPS_Degree}
    ctps_new = CTPS(ctps)
    ctps_new.map[1] += a
    return ctps_new
end

function +(a::Float64, ctps::CTPS{ComplexF64, TPS_Dim, Max_TPS_Degree}) where {TPS_Dim, Max_TPS_Degree}
    return ctps + a
end
# -
function -(ctps1::CTPS{T, TPS_Dim, Max_TPS_Degree}, ctps2::CTPS{T, TPS_Dim, Max_TPS_Degree}) where {T, TPS_Dim, Max_TPS_Degree}
    ctps_new = CTPS(ctps1)
    for i in eachindex(ctps_new.map)
        ctps_new.map[i] -= ctps2.map[i]
    end
    return ctps_new
end
function -(ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}, a::T) where {T, TPS_Dim, Max_TPS_Degree}
    ctps_new = CTPS(ctps)
    ctps_new.map[1] -= a
    return ctps_new
end
function -(a::T, ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}) where {T, TPS_Dim, Max_TPS_Degree}
    ctps_new = CTPS(ctps)
    ctps_new.map[1] = a - ctps_new.map[1]
    for i in eachindex(ctps_new.map)[2:end]
        ctps_new.map[i] = -ctps_new.map[i]
    end
    return ctps_new
end
function -(ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}) where {T, TPS_Dim, Max_TPS_Degree}
    ctps_new = CTPS(ctps)
    for i in eachindex(ctps_new.map)
        ctps_new.map[i] = -ctps_new.map[i]
    end
    return ctps_new
end
# operate with complex number
function -(ctps::CTPS{ComplexF64, TPS_Dim, Max_TPS_Degree}, a::Float64) where {TPS_Dim, Max_TPS_Degree}
    ctps_new = CTPS(ctps)
    ctps_new.map[1] -= a
    return ctps_new
end
function -(a::Float64, ctps::CTPS{ComplexF64, TPS_Dim, Max_TPS_Degree}) where {TPS_Dim, Max_TPS_Degree}
    ctps_new = CTPS(ctps)
    ctps_new.map[1] = a - ctps_new.map[1]
    for i in eachindex(ctps_new.map)[2:end]
        ctps_new.map[i] = -ctps_new.map[i]
    end
    return ctps_new
end

# *
function *(ctps1::CTPS{T, TPS_Dim, Max_TPS_Degree}, ctps2::CTPS{T, TPS_Dim, Max_TPS_Degree}) where {T, TPS_Dim, Max_TPS_Degree}
    ctps_new = CTPS(ctps1)
    ctps_map_buffer = ctps_new.map
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
            index = findindex(ctps_new, temp)
            ctps_map_buffer[index] += ctps1.map[i] * ctps2.map[j]
        end
    end
    return ctps_new
end
function *(ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}, a::T) where {T, TPS_Dim, Max_TPS_Degree}
    ctps_new = CTPS(ctps)
    for i in eachindex(ctps_new.map)
        ctps_new.map[i] *= a
    end
    return ctps_new
end
function *(a::T, ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}) where {T, TPS_Dim, Max_TPS_Degree}
    return ctps * a
end

# operate with complex number
function *(ctps::CTPS{ComplexF64, TPS_Dim, Max_TPS_Degree}, a::Float64) where {TPS_Dim, Max_TPS_Degree}
    ctps_new = CTPS(ctps)
    for i in eachindex(ctps_new.map)
        ctps_new.map[i] *= a
    end
    return ctps_new
end
function *(a::Float64, ctps::CTPS{ComplexF64, TPS_Dim, Max_TPS_Degree}) where {TPS_Dim, Max_TPS_Degree}
    return ctps * a
end

# /
function inv(ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}) where {T, TPS_Dim, Max_TPS_Degree}
    if cst(ctps) == zero(T)
        error("Divide by zero in CTPS")
    end
    temp = CTPS(ctps) - cst(ctps)
    term_by_oder = CTPS(1.0 / ctps.map[1], TPS_Dim, Max_TPS_Degree)

    SUM = CTPS(term_by_oder)
    for i in 1:Max_TPS_Degree
        term_by_oder *= -temp/cst(ctps)
        SUM += term_by_oder
    end
    return SUM
end

function /(ctps1::CTPS{T, TPS_Dim, Max_TPS_Degree}, ctps2::CTPS{T, TPS_Dim, Max_TPS_Degree}) where {T, TPS_Dim, Max_TPS_Degree}
    if cst(ctps2) == zero(T)
        error("Divide by zero in CTPS")
    end

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
    return ctps_new
end
function /(a::T, ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}) where {T, TPS_Dim, Max_TPS_Degree}
    if cst(ctps) == zero(T)
        error("Divide by zero in CTPS")
    end
    return a * inv(ctps)
end

# operate with complex number
function /(ctps::CTPS{ComplexF64, TPS_Dim, Max_TPS_Degree}, a::Float64) where {TPS_Dim, Max_TPS_Degree}
    if a == zero(Float64)
        error("Divide by zero in CTPS")
    end
    ctps_new = CTPS(ctps)
    for i in eachindex(ctps_new.map)
        ctps_new.map[i] /= a
    end
    return ctps_new
end
function /(a::Float64, ctps::CTPS{ComplexF64, TPS_Dim, Max_TPS_Degree}) where {TPS_Dim, Max_TPS_Degree}
    if cst(ctps) == zero(ComplexF64)
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
    SUM = CTPS(one(T), TPS_Dim, Max_TPS_Degree)
    for i in 1:Max_TPS_Degree
        index = 1.0 / factorial_double(i)
        term_by_oder = term_by_oder * temp
        SUM = SUM + (term_by_oder * T(index))
    end
    SUM = SUM * T(Base.exp(cst(ctps)))
    return SUM
end

# logarithm
function log(ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}) where {T, TPS_Dim, Max_TPS_Degree}
    if cst(ctps) == zero(T)
        error("Log of zero in CTPS")
    end
    temp = CTPS(ctps)
    temp -= cst(temp)
    term_by_oder = temp / cst(ctps)
    SUM = CTPS(zero(T), TPS_Dim, Max_TPS_Degree) + term_by_oder

    for i in 2:Max_TPS_Degree
        term_by_oder = term_by_oder * (-temp / cst(ctps))
        SUM = SUM + (term_by_oder / T(i))
    end
    SUM = SUM + Base.log(cst(ctps))
    return SUM
end

# square root
function sqrt(ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}) where {T, TPS_Dim, Max_TPS_Degree}
    if real(cst(ctps)) < 0.0
        error("Square root of negative number in CTPS")
    end
    a0 = sqrt(cst(ctps))
    temp = CTPS(ctps) - cst(ctps)
    # temp = temp 
    term_by_oder = temp / a0
    SUM = term_by_oder/2.0

    for i in 2:Max_TPS_Degree
        index = 1.0 * doublefactorial_double(2 * i - 3) / doublefactorial_double(2 * i)
        term_by_oder = (-temp) * term_by_oder / cst(ctps)
        SUM += term_by_oder * index
    end
    SUM = SUM + a0
    return SUM
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
            SUM = CTPS(ctps)
            for i in 2:b
                SUM = SUM * ctps
            end
            return SUM
        else
            error("Divide by zero, in CTPS::pow")
        end
    end
    temp = temp - cst(temp)
    term_by_oder = CTPS(one(T), TPS_Dim, Max_TPS_Degree)
    factor = cst(ctps) ^ b
    SUM = CTPS(factor, TPS_Dim, Max_TPS_Degree)

    for i in 1:Max_TPS_Degree
        factor = factor / cst(ctps) * index / i
        index -= 1.0
        term_by_oder = term_by_oder * temp
        SUM = SUM + (term_by_oder * factor)
        if index == 0.0
            break
        end
    end
    return SUM
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
    SUM = CTPS(zero(T), TPS_Dim, Max_TPS_Degree)

    is_odd_iteration = 1  
    for i in 1:Max_TPS_Degree
        if is_odd_iteration == 1
            index = cos_a0 * (-1) ^ ((i - 1) / 2) / factorial_double(i)
        else
            index = sin_a0 * (-1) ^ (i / 2) / factorial_double(i)
        end
        term_by_oder = term_by_oder * temp
        SUM = SUM + (term_by_oder * index)

        is_odd_iteration = -is_odd_iteration 
    end
    SUM = SUM + sin_a0
    return SUM
end

# cos
function cos(ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}) where {T, TPS_Dim, Max_TPS_Degree}
    temp = CTPS(ctps)
    a0 = cst(ctps)
    sin_a0 = sin(a0)
    cos_a0 = cos(a0)
    temp = temp - a0
    term_by_oder = CTPS(one(T), TPS_Dim, Max_TPS_Degree)
    SUM = CTPS(zero(T), TPS_Dim, Max_TPS_Degree)

    is_odd_iteration = 1  
    for i in 1:Max_TPS_Degree
        if is_odd_iteration == 1
            index = sin_a0 * (-1) ^ ((i + 1) / 2) / factorial_double(i)
        else
            index = cos_a0 * (-1) ^ (i / 2) / factorial_double(i)
        end
        term_by_oder = term_by_oder * temp
        SUM = SUM + (term_by_oder * index)
    
        is_odd_iteration = -is_odd_iteration  
    end
    SUM = SUM + cos_a0
    return SUM
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
    SUM = CTPS(zero(T), TPS_Dim, Max_TPS_Degree)

    is_odd_iteration = 1  
    for i in 1:Max_TPS_Degree
        if is_odd_iteration == 1
            index = cosh_a0 / factorial_double(i)
        else
            index = sinh_a0 / factorial_double(i)
        end
        term_by_oder = term_by_oder * temp
        SUM = SUM + (term_by_oder * index)
        is_odd_iteration = -is_odd_iteration
    end
    SUM = SUM + sinh_a0
    return SUM
end

# hyperbolic cos
function cosh(ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}) where {T, TPS_Dim, Max_TPS_Degree}
    temp = CTPS(ctps)
    a0 = cst(ctps)
    sinh_a0 = sinh(a0)
    cosh_a0 = cosh(a0)
    temp = temp - a0
    term_by_oder = CTPS(one(T), TPS_Dim, Max_TPS_Degree)
    SUM = CTPS(zero(T), TPS_Dim, Max_TPS_Degree)

    is_odd_iteration = 1
    for i in 1:Max_TPS_Degree
        if is_odd_iteration == 1
            index = sinh_a0 / factorial_double(i)
        else
            index = cosh_a0 / factorial_double(i)
        end
        term_by_oder = term_by_oder * temp
        SUM = SUM + (term_by_oder * index)
        is_odd_iteration = -is_odd_iteration
    end
    SUM = SUM + cosh_a0
    return SUM
end