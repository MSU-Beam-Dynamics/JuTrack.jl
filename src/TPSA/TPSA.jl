# # This file is part of the multivariate high-order TPSA (Truncated Power Series Algebra) package.
# # Translated from Yue Hao's C++ code
# # Author: Jinyu Wan
# # Email: wan@frib.msu.edu
# # Created Date: 11-01-2023
# # Modified Date: 08-25-2025

module HighOrderTPS

include("polymap.jl")

export CTPS, cst, assign!, reassign!

"""
    struct CTPS{T, TPS_Dim, Max_TPS_Degree}
A truncated power series in `TPS_Dim` variables with coefficients
of type `T` and maximum total degree `Max_TPS_Degree`.
"""
struct CTPS{T, TPS_Dim, Max_TPS_Degree}
    degree::Int
    terms::Int
    map::Vector{T}
    polymap::Ref{PolyMap}
end

"""
Global cache for storing and reusing `PolyMap` objects.  
"""
const polyMapCache = Dict{Tuple{Int, Int}, PolyMap}()

"""
    getOrCreatePolyMap(TPS_Dim::Int, Max_TPS_Degree::Int) -> PolyMap
Return an existing `PolyMap` from the global cache or construct a new
one on demand. 
"""
function getOrCreatePolyMap(TPS_Dim::Int, Max_TPS_Degree::Int)
    key = (TPS_Dim, Max_TPS_Degree)
    if @inbounds haskey(polyMapCache, key)
        return polyMapCache[key]
    else
        newPolyMap = PolyMap(TPS_Dim, Max_TPS_Degree)
        polyMapCache[key] = newPolyMap
        return newPolyMap
    end
end

"""
    CTPS(T::Type, TPS_Dim::Int, Max_TPS_Degree::Int)

Construct a zero power series of element type `T`, dimension
`TPS_Dim` and maximum total degree `Max_TPS_Degree`. 
"""
function CTPS(T::Type, TPS_Dim::Int, Max_TPS_Degree::Int)
    polymap = getOrCreatePolyMap(TPS_Dim, Max_TPS_Degree)
    terms  = binomial(TPS_Dim + Max_TPS_Degree, Max_TPS_Degree)
    map    = Vector{T}(undef, terms)
    fill!(map, zero(T))
    return CTPS{T, TPS_Dim, Max_TPS_Degree}(0, terms, map, Ref(polymap))
end

"""
    CTPS(a::T, TPS_Dim::Int, Max_TPS_Degree::Int)

Construct a power series whose constant term is `a` and all
higher-order terms are zero.
"""
function CTPS(a::T, TPS_Dim::Int, Max_TPS_Degree::Int) where T
    polymap = getOrCreatePolyMap(TPS_Dim, Max_TPS_Degree)
    terms  = binomial(TPS_Dim + Max_TPS_Degree, Max_TPS_Degree)
    map    = Vector{T}(undef, terms)
    fill!(map, zero(T))
    map[1] = a
    return CTPS{T, TPS_Dim, Max_TPS_Degree}(0, terms, map, Ref(polymap))
end

"""
    CTPS(a::T, n::Int, TPS_Dim::Int, Max_TPS_Degree::Int)

Construct a series representing `a + x_n`.  The constant term is set
to `a` and the coefficient of the n-th variable is one. 
"""
function CTPS(a::T, n::Int, TPS_Dim::Int, Max_TPS_Degree::Int) where T
    n > 0 && n ≤ TPS_Dim || throw(ArgumentError("variable index out of range"))
    polymap = getOrCreatePolyMap(TPS_Dim, Max_TPS_Degree)
    terms  = binomial(TPS_Dim + Max_TPS_Degree, Max_TPS_Degree)
    map    = Vector{T}(undef, terms)
    fill!(map, zero(T))
    map[1]     = a
    map[n + 1] = one(T)
    return CTPS{T, TPS_Dim, Max_TPS_Degree}(1, terms, map, Ref(polymap))
end

"""
    CTPS(M::CTPS)

Copy constructor for truncated power series. 
"""
function CTPS(M::CTPS{T, TPS_Dim, Max_TPS_Degree}) where {T, TPS_Dim, Max_TPS_Degree}
    map = copy(M.map)
    return CTPS{T, TPS_Dim, Max_TPS_Degree}(M.degree, length(map), map, M.polymap)
end

"""
    cst(ctps::CTPS)

Return the constant term (zeroth order coefficient) of the truncated
power series `ctps`.
"""
@inline function cst(ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}) where {T, TPS_Dim, Max_TPS_Degree}
    return ctps.map[1]
end

"""
    findindex(ctps, indexmap)

Given an exponent vector `indexmap`, return the corresponding linear
index into the coefficient array.  
"""
function findindex(ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}, indexmap::AbstractVector{Int}) where {T, TPS_Dim, Max_TPS_Degree}
    dim = TPS_Dim
    len = length(indexmap)
    if len == dim
        # Compute total degree and construct a new vector of length dim+1
        total = zero(Int)
        newmap = Vector{Int}(undef, dim + 1)
        @inbounds for i in 1:dim
            EXP = indexmap[i]
            EXP ≥ 0 || throw(ArgumentError("negative exponent in indexmap"))
            total += EXP
            newmap[i + 1] = EXP
        end
        newmap[1] = total
        return findindex(ctps, newmap)
    elseif len == dim + 1
        # First entry is the total degree
        total = indexmap[1]
        total ≥ 0 || throw(ArgumentError("negative total degree in indexmap"))
        result = 1
        # Walk from the highest dimension downwards, updating the cumulative sum
        # `current_sum` represents the remaining degree after accounting for
        # exponents of the lower dimensions.
        current_sum = total
        @inbounds for i in dim:-1:1
            # If no degree remains, all higher terms vanish
            if current_sum == 0
                break
            end
            result += binomial(current_sum - 1 + i, i)
            # Subtract the exponent of the current variable (indexed by dim - i + 2)
            current_sum -= indexmap[dim - i + 2]
        end
        return result
    else
        throw(ArgumentError("indexmap must have length $(dim) or $(dim + 1)"))
    end
end

"""
    assign!(ctps, a)

Assign the constant term of the truncated power series `ctps` to `a`. 
"""
@inline function assign!(ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}, a::T) where {T, TPS_Dim, Max_TPS_Degree}
    ctps.map[1] = a
    return nothing
end

"""
    assign!(ctps, a, n_var)

Assign `ctps` to represent the constant `a` plus a unit coefficient in
dimension `n_var`.  
"""
function assign!(ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}, a::T, n_var::Int) where {T, TPS_Dim, Max_TPS_Degree}
    n_var > 0 && n_var ≤ TPS_Dim || throw(ArgumentError("variable index out of range"))
    ctps.map[n_var + 1] = one(T)
    ctps.map[1]        = a
    return nothing
end

"""
    reassign!(ctps, a, n_var)

Reset all coefficients of `ctps` to zero, then assign the constant
term to `a` and the coefficient of the `n_var`-th variable to one.
"""
function reassign!(ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}, a::T, n_var::Int) where {T, TPS_Dim, Max_TPS_Degree}
    n_var > 0 && n_var ≤ TPS_Dim || throw(ArgumentError("variable index out of range"))
    fill!(ctps.map, zero(T))
    ctps.map[1]        = a
    ctps.map[n_var + 1] = one(T)
    return nothing
end

"""
    element(ctps, ind)

Return the coefficient corresponding to the exponent vector `ind`.
The exponent vector may be of length `TPS_Dim` or `TPS_Dim + 1`.
"""
@inline function element(ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}, ind::AbstractVector{Int}) where {T, TPS_Dim, Max_TPS_Degree}
    idx = findindex(ctps, ind)
    return ctps.map[idx]
end

#############################################################################################################################
##  Overloaded arithmetic and transcendental operations
#############################################################################################################################

import Base: +, -, *, /, ^, sin, cos, tan, exp, log, sinh, cosh, tanh, asin, acos, atan, asinh, acosh, atanh, sqrt

function +(ctps1::CTPS{T, TPS_Dim, Max_TPS_Degree}, ctps2::CTPS{T, TPS_Dim, Max_TPS_Degree}) where {T, TPS_Dim, Max_TPS_Degree}
    ctps_new = CTPS(ctps1)
    @inbounds for i in eachindex(ctps_new.map)
        ctps_new.map[i] += ctps2.map[i]
    end
    return ctps_new
end

function +(ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}, a::T) where {T, TPS_Dim, Max_TPS_Degree}
    ctps_new = CTPS(ctps)
    ctps_new.map[1] += a
    return ctps_new
end

@inline +(a::T, ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}) where {T, TPS_Dim, Max_TPS_Degree} = ctps + a

function -(ctps1::CTPS{T, TPS_Dim, Max_TPS_Degree}, ctps2::CTPS{T, TPS_Dim, Max_TPS_Degree}) where {T, TPS_Dim, Max_TPS_Degree}
    ctps_new = CTPS(ctps1)
    @inbounds for i in eachindex(ctps_new.map)
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
    # Reverse the sign of all coefficients
    @inbounds for i in eachindex(ctps_new.map)
        ctps_new.map[i] = -ctps_new.map[i]
    end
    ctps_new.map[1] += a
    return ctps_new
end

function -(ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}) where {T, TPS_Dim, Max_TPS_Degree}
    ctps_new = CTPS(ctps)
    @inbounds for i in eachindex(ctps_new.map)
        ctps_new.map[i] = -ctps_new.map[i]
    end
    return ctps_new
end

function *(ctps1::CTPS{T, TPS_Dim, Max_TPS_Degree}, ctps2::CTPS{T, TPS_Dim, Max_TPS_Degree}) where {T, TPS_Dim, Max_TPS_Degree}
    # Initialise an empty result with zero coefficients
    ctps_new = CTPS{T, TPS_Dim, Max_TPS_Degree}(0, ctps1.terms, Vector{T}(undef, ctps1.terms), ctps1.polymap)
    fill!(ctps_new.map, zero(T))
    # Temporary buffer for exponent addition
    temp = Vector{Int}(undef, TPS_Dim + 1)
    # Loop over all non-zero terms of ctps1 and ctps2
    @inbounds for i in eachindex(ctps1.map)
        c1 = ctps1.map[i]
        if c1 == zero(T)
            continue
        end
        exp1 = getindexmap(ctps1.polymap[], i)
        # The maximum degree that can be achieved by multiplying exp1 with exp2
        deg1 = exp1[1]
        max_j_terms = binomial(TPS_Dim + (Max_TPS_Degree - deg1), TPS_Dim)
        # Restrict iteration to avoid computing terms beyond the truncation
        j_lim = min(ctps2.terms, max_j_terms)
        for j in 1:j_lim
            c2 = ctps2.map[j]
            if c2 == zero(T)
                continue
            end
            exp2 = getindexmap(ctps2.polymap[], j)
            # Compute the sum of exponent vectors exp1 and exp2 into temp
            @inbounds for k in 1:(TPS_Dim + 1)
                temp[k] = exp1[k] + exp2[k]
            end
            # Insert the product into the appropriate position
            idx = findindex(ctps_new, temp)
            ctps_new.map[idx] += c1 * c2
        end
    end
    return ctps_new
end

function *(ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}, a::T) where {T, TPS_Dim, Max_TPS_Degree}
    ctps_new = CTPS(ctps)
    @inbounds for i in eachindex(ctps_new.map)
        ctps_new.map[i] *= a
    end
    return ctps_new
end

@inline *(a::T, ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}) where {T, TPS_Dim, Max_TPS_Degree} = ctps * a

function /(ctps1::CTPS{T, TPS_Dim, Max_TPS_Degree}, ctps2::CTPS{T, TPS_Dim, Max_TPS_Degree}) where {T, TPS_Dim, Max_TPS_Degree}
    c0 = cst(ctps2)
    c0 == zero(T) && throw(DivideError())
    return ctps1 * inv(ctps2)
end

function /(ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}, a::T) where {T, TPS_Dim, Max_TPS_Degree}
    a == zero(T) && throw(DivideError())
    ctps_new = CTPS(ctps)
    @inbounds for i in eachindex(ctps_new.map)
        ctps_new.map[i] /= a
    end
    return ctps_new
end

function /(a::T, ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}) where {T, TPS_Dim, Max_TPS_Degree}
    c0 = cst(ctps)
    c0 == zero(T) && throw(DivideError())
    return a * inv(ctps)
end

###############################################################
##  Utility functions for reciprocals and elementary functions
###############################################################

function inv(ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}) where {T, TPS_Dim, Max_TPS_Degree}
    c0 = cst(ctps)
    c0 == zero(T) && throw(DivideError())
    # Separate out the non-constant part
    temp = ctps - c0
    inv_c0 = one(T) / c0
    term_by_order = CTPS(inv_c0, TPS_Dim, Max_TPS_Degree)
    SUM = CTPS(term_by_order)
    for i in 1:Max_TPS_Degree
        term_by_order *= -(temp / c0)
        SUM += term_by_order
    end
    return SUM
end

function exp(ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}) where {T, TPS_Dim, Max_TPS_Degree}
    a0 = cst(ctps)
    # When the constant term is zero, exp(ctps) starts at one
    if a0 == zero(T)
        return CTPS(one(T), TPS_Dim, Max_TPS_Degree)
    end
    temp = ctps - a0
    term_by_order = CTPS(one(T), TPS_Dim, Max_TPS_Degree)
    SUM = CTPS(one(T), TPS_Dim, Max_TPS_Degree)
    for i in 1:Max_TPS_Degree
        term_by_order *= temp
        coeff = one(T) / T(factorial_double(i))
        SUM += term_by_order * coeff
    end
    return SUM * T(exp(a0))
end

function log(ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}) where {T, TPS_Dim, Max_TPS_Degree}
    a0 = cst(ctps)
    a0 == zero(T) && throw(DomainError(ctps, "log of zero in CTPS"))
    temp = ctps - a0
    term_by_order = temp / a0
    SUM = CTPS(zero(T), TPS_Dim, Max_TPS_Degree) + term_by_order
    for i in 2:Max_TPS_Degree
        term_by_order *= -(temp / a0)
        SUM += term_by_order / T(i)
    end
    return SUM + T(log(a0))
end

function sqrt(ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}) where {T<:Real, TPS_Dim, Max_TPS_Degree}
    a0 = cst(ctps)
    a0 < zero(T) && throw(DomainError(ctps, "square root of negative number in CTPS"))
    root_a0 = sqrt(a0)
    temp = ctps - a0
    term_by_order = temp / root_a0
    SUM = term_by_order / T(2)
    for i in 2:Max_TPS_Degree
        coeff = doublefactorial_double(2 * i - 3) / doublefactorial_double(2 * i)
        term_by_order *= -temp / a0
        SUM += term_by_order * T(coeff)
    end
    return SUM + root_a0
end

function pow(ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}, b::Int) where {T, TPS_Dim, Max_TPS_Degree}
    b == 1 && return ctps
    b == 0 && return CTPS(one(T), TPS_Dim, Max_TPS_Degree)
    if b > 1
        # Handle simple positive integer powers by repeated multiplication
        result = CTPS(ctps)
        for i in 2:b
            result *= ctps
        end
        return result
    end
    # Negative or zero constant term requires inversion
    a0 = cst(ctps)
    a0 == zero(T) && throw(DomainError(ctps, "power with zero constant term"))
    temp = ctps - a0
    index = T(b)
    factor = a0^b
    term_by_order = CTPS(one(T), TPS_Dim, Max_TPS_Degree)
    SUM = CTPS(factor, TPS_Dim, Max_TPS_Degree)
    for i in 1:Max_TPS_Degree
        factor = factor / a0 * (index / T(i))
        index -= one(T)
        term_by_order *= temp
        SUM += term_by_order * factor
        index == zero(T) && break
    end
    return SUM
end

@inline ^(ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}, b::Int) where {T, TPS_Dim, Max_TPS_Degree} = pow(ctps, b)

###############################################################
##  Trigonometric and hyperbolic functions
###############################################################

function sin(ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}) where {T, TPS_Dim, Max_TPS_Degree}
    a0 = cst(ctps)
    sin_a0 = sin(a0)
    cos_a0 = cos(a0)
    temp = ctps - a0
    term_by_order = CTPS(one(T), TPS_Dim, Max_TPS_Degree)
    SUM = CTPS(zero(T), TPS_Dim, Max_TPS_Degree)
    odd = true
    for i in 1:Max_TPS_Degree
        coeff = if odd
            T(cos_a0 * ((-1) ^ ((i - 1) ÷ 2)) / factorial_double(i))
        else
            T(sin_a0 * ((-1) ^ (i ÷ 2)) / factorial_double(i))
        end
        term_by_order *= temp
        SUM += term_by_order * coeff
        odd = !odd
    end
    return SUM + sin_a0
end

function cos(ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}) where {T, TPS_Dim, Max_TPS_Degree}
    a0 = cst(ctps)
    sin_a0 = sin(a0)
    cos_a0 = cos(a0)
    temp = ctps - a0
    term_by_order = CTPS(one(T), TPS_Dim, Max_TPS_Degree)
    SUM = CTPS(zero(T), TPS_Dim, Max_TPS_Degree)
    odd = true
    for i in 1:Max_TPS_Degree
        coeff = if odd
            T(sin_a0 * ((-1) ^ ((i + 1) ÷ 2)) / factorial_double(i))
        else
            T(cos_a0 * ((-1) ^ (i ÷ 2)) / factorial_double(i))
        end
        term_by_order *= temp
        SUM += term_by_order * coeff
        odd = !odd
    end
    return SUM + cos_a0
end

@inline function tan(ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}) where {T, TPS_Dim, Max_TPS_Degree}
    return sin(ctps) / cos(ctps)
end

function asin(ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}) where {T, TPS_Dim, Max_TPS_Degree}
    a0 = cst(ctps)
    arcsin_a0 = asin(a0)
    cos_y0 = sqrt(one(T) - a0^2)
    temp = ctps - a0
    y = CTPS(temp)
    # Newton iteration to solve y = asin(x): y + sin(y) - (x - a0) = a0*(1 - cos(y))/cos_y0
    for _ in 1:(Max_TPS_Degree + 10)
        y -= sin(y) - (temp + a0 * (one(T) - cos(y))) / cos_y0
    end
    return y + arcsin_a0
end

@inline function acos(ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}) where {T, TPS_Dim, Max_TPS_Degree}
    return T(pi/2) - asin(ctps)
end

@inline function atan(ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}) where {T, TPS_Dim, Max_TPS_Degree}
    return asin(ctps / sqrt(one(T) + ctps*ctps))
end

function sinh(ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}) where {T, TPS_Dim, Max_TPS_Degree}
    a0 = cst(ctps)
    sinh_a0 = sinh(a0)
    cosh_a0 = cosh(a0)
    temp = ctps - a0
    term_by_order = CTPS(one(T), TPS_Dim, Max_TPS_Degree)
    SUM = CTPS(zero(T), TPS_Dim, Max_TPS_Degree)
    odd = true
    for i in 1:Max_TPS_Degree
        coeff = if odd
            T(cosh_a0 / factorial_double(i))
        else
            T(sinh_a0 / factorial_double(i))
        end
        term_by_order *= temp
        SUM += term_by_order * coeff
        odd = !odd
    end
    return SUM + sinh_a0
end

function cosh(ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}) where {T, TPS_Dim, Max_TPS_Degree}
    a0 = cst(ctps)
    sinh_a0 = sinh(a0)
    cosh_a0 = cosh(a0)
    temp = ctps - a0
    term_by_order = CTPS(one(T), TPS_Dim, Max_TPS_Degree)
    SUM = CTPS(zero(T), TPS_Dim, Max_TPS_Degree)
    odd = true
    for i in 1:Max_TPS_Degree
        coeff = if odd
            T(sinh_a0 / factorial_double(i))
        else
            T(cosh_a0 / factorial_double(i))
        end
        term_by_order *= temp
        SUM += term_by_order * coeff
        odd = !odd
    end
    return SUM + cosh_a0
end

end 