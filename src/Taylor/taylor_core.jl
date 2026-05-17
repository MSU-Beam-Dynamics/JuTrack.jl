module TaylorCore

using StaticArrays: SVector
using ..TPSAadStatic: DTPSAD
using ..HighOrderTPS: CTPS, cst, findindex

abstract type AbstractTaylorBackend end
struct DTPSADBackend <: AbstractTaylorBackend end
struct CTPSBackend <: AbstractTaylorBackend end
struct HOTPSABackend <: AbstractTaylorBackend end

struct MonomialTable{N,O}
    exponents::Vector{NTuple{N,Int}}
    degree_offsets::Vector{Int}
    index_map::Dict{NTuple{N,Int},Int}
end

Base.length(table::MonomialTable) = length(table.exponents)

const _MONOMIAL_TABLE_CACHE = Dict{Tuple{Int,Int},Any}()

@inline jet_backend(::Type{DTPSAD{N,T}}) where {N,T} = DTPSADBackend()
@inline jet_backend(::DTPSAD{N,T}) where {N,T} = DTPSADBackend()
@inline jet_backend(::Type{CTPS{T,N,O}}) where {T,N,O} = CTPSBackend()
@inline jet_backend(::CTPS{T,N,O}) where {T,N,O} = CTPSBackend()

@inline jet_nvars(::Type{DTPSAD{N,T}}) where {N,T} = N
@inline jet_nvars(::DTPSAD{N,T}) where {N,T} = N
@inline jet_nvars(::Type{CTPS{T,N,O}}) where {T,N,O} = N
@inline jet_nvars(::CTPS{T,N,O}) where {T,N,O} = N

@inline jet_order(::Type{DTPSAD{N,T}}) where {N,T} = 1
@inline jet_order(::DTPSAD{N,T}) where {N,T} = 1
@inline jet_order(::Type{CTPS{T,N,O}}) where {T,N,O} = O
@inline jet_order(::CTPS{T,N,O}) where {T,N,O} = O

@inline jet_primal(x::DTPSAD{N,T}) where {N,T} = x.val
@inline jet_primal(x::CTPS{T,N,O}) where {T,N,O} = cst(x)

@inline function _unit_exponents(::Val{N}, i::Int) where {N}
    return ntuple(j -> j == i ? 1 : 0, N)
end

function _generate_degree_exponents(::Val{1}, degree::Int)
    return NTuple{1,Int}[(degree,)]
end

function _generate_degree_exponents(::Val{N}, degree::Int) where {N}
    exponents = NTuple{N,Int}[]
    for head in 0:degree
        for tail in _generate_degree_exponents(Val(N - 1), degree - head)
            push!(exponents, (head, tail...))
        end
    end
    return exponents
end

function _ctps_sort_index(::Val{N}, ::Val{O}, exponent::NTuple{N,Int}) where {N,O}
    proto = CTPS(Float64, N, O)
    return findindex(proto, collect(exponent))
end

function _build_monomial_table(::Val{N}, ::Val{O}) where {N,O}
    exponents = NTuple{N,Int}[]
    degree_offsets = Vector{Int}(undef, O + 2)
    next_index = 1
    for degree in 0:O
        degree_offsets[degree + 1] = next_index
        degree_terms = _generate_degree_exponents(Val(N), degree)
        sort!(degree_terms, by = exponent -> _ctps_sort_index(Val(N), Val(O), exponent))
        append!(exponents, degree_terms)
        next_index += length(degree_terms)
    end
    degree_offsets[O + 2] = next_index
    index_map = Dict{NTuple{N,Int},Int}()
    for (i, exponent) in enumerate(exponents)
        index_map[exponent] = i
    end
    return MonomialTable{N,O}(exponents, degree_offsets, index_map)
end

function monomial_table(::Val{N}, ::Val{O}) where {N,O}
    key = (N, O)
    if !haskey(_MONOMIAL_TABLE_CACHE, key)
        _MONOMIAL_TABLE_CACHE[key] = _build_monomial_table(Val(N), Val(O))
    end
    return _MONOMIAL_TABLE_CACHE[key]::MonomialTable{N,O}
end

@inline monomial_table(::Type{DTPSAD{N,T}}) where {N,T} = monomial_table(Val(N), Val(1))
@inline monomial_table(::Type{CTPS{T,N,O}}) where {T,N,O} = monomial_table(Val(N), Val(O))

@inline function monomial_index(table::MonomialTable{N,O}, exponent::NTuple{N,Int}) where {N,O}
    return get(table.index_map, exponent, 0)
end

@inline monomial_exponents(table::MonomialTable, index::Int) = table.exponents[index]

function degree_range(table::MonomialTable{N,O}, degree::Int) where {N,O}
    0 <= degree <= O || throw(ArgumentError("degree $degree is outside 0:$O"))
    return table.degree_offsets[degree + 1]:(table.degree_offsets[degree + 2] - 1)
end

function monomial_product_index(table::MonomialTable{N,O}, left::NTuple{N,Int}, right::NTuple{N,Int}) where {N,O}
    total = ntuple(i -> left[i] + right[i], N)
    return sum(total) > O ? 0 : monomial_index(table, total)
end

function monomial_product_index(table::MonomialTable{N,O}, left::Int, right::Int) where {N,O}
    return monomial_product_index(table, table.exponents[left], table.exponents[right])
end

function jet_type(::DTPSADBackend, nvars::Integer; order::Integer=1, T::Type{<:Number}=Float64)
    order == 1 || throw(ArgumentError("DTPSAD backend supports only first order"))
    nvars > 0 || throw(ArgumentError("nvars must be positive"))
    return DTPSAD{Int(nvars),T}
end

function jet_type(::CTPSBackend, nvars::Integer; order::Integer, T::Type{<:Number}=Float64)
    nvars > 0 || throw(ArgumentError("nvars must be positive"))
    order >= 0 || throw(ArgumentError("order must be non-negative"))
    return CTPS{T,Int(nvars),Int(order)}
end

function jet_seed(::Type{DTPSAD{N,T}}, value::Real, i::Integer) where {N,T}
    1 <= i <= N || throw(ArgumentError("seed index $i is outside 1:$N"))
    deriv = SVector{N,T}(ntuple(j -> j == i ? one(T) : zero(T), N))
    return DTPSAD{N,T}(convert(T, value), deriv)
end

function jet_seed(::Type{CTPS{T,N,O}}, value::Real, i::Integer) where {T,N,O}
    1 <= i <= N || throw(ArgumentError("seed index $i is outside 1:$N"))
    return CTPS(convert(T, value), i, N, O)
end

jet_constant(x) = jet_primal(x)

function jet_coefficient(x::DTPSAD{N,T}, exponent::NTuple{N,Int}) where {N,T}
    degree = sum(exponent)
    if degree == 0
        return x.val
    elseif degree == 1
        found = 0
        @inbounds for i in 1:N
            if exponent[i] == 1
                found == 0 || return zero(T)
                found = i
            elseif exponent[i] != 0
                return zero(T)
            end
        end
        return found == 0 ? zero(T) : x.deriv[found]
    end
    return zero(T)
end

function jet_coefficient(x::CTPS{T,N,O}, exponent::NTuple{N,Int}) where {T,N,O}
    sum(exponent) <= O || return zero(T)
    table = monomial_table(Val(N), Val(O))
    index = monomial_index(table, exponent)
    return index == 0 ? zero(T) : x.map[index]
end

function jet_gradient(x)
    N = jet_nvars(x)
    return [jet_coefficient(x, _unit_exponents(Val(N), i)) for i in 1:N]
end

function jet_hessian(x)
    N = jet_nvars(x)
    primal_type = typeof(jet_primal(x))
    H = zeros(primal_type, N, N)
    for i in 1:N
        for j in i:N
            exponent = ntuple(k -> k == i ? 1 : 0, N)
            if i == j
                exponent = ntuple(k -> k == i ? 2 : 0, N)
                value = 2 * jet_coefficient(x, exponent)
            else
                exponent = ntuple(k -> (k == i || k == j) ? 1 : 0, N)
                value = jet_coefficient(x, exponent)
            end
            H[i, j] = value
            H[j, i] = value
        end
    end
    return H
end

function seeded_jets(::Type{S}, x::AbstractVector{<:Real}) where {S}
    return [jet_seed(S, x[i], i) for i in eachindex(x)]
end

function _unwrap_scalar_output(y)
    if y isa AbstractVector
        length(y) == 1 || throw(ArgumentError("function must return a single value"))
        return y[1]
    elseif y isa Tuple
        throw(ArgumentError("function must return a single value"))
    end
    return y
end

function _unwrap_vector_output(y)
    return y isa AbstractVector ? y : collect(y)
end

function TaylorGradient(::Type{S}, f, x::AbstractVector{<:Real}; primal::Bool=false) where {S}
    X = seeded_jets(S, x)
    y = _unwrap_scalar_output(f(X...))
    gradient = jet_gradient(y)
    return primal ? (gradient, jet_primal(y)) : gradient
end

function TaylorJacobian(::Type{S}, F, x::AbstractVector{<:Real}; primal::Bool=false) where {S}
    X = seeded_jets(S, x)
    outputs = _unwrap_vector_output(F(X...))
    m = length(outputs)
    n = length(x)
    J = Matrix{Float64}(undef, m, n)
    values = Vector{Float64}(undef, m)
    @inbounds for row in 1:m
        J[row, :] .= Float64.(jet_gradient(outputs[row]))
        values[row] = Float64(jet_primal(outputs[row]))
    end
    return primal ? (J, values) : J
end

export AbstractTaylorBackend, DTPSADBackend, CTPSBackend, HOTPSABackend
export MonomialTable, monomial_table, monomial_index, monomial_exponents, degree_range, monomial_product_index
export jet_backend, jet_nvars, jet_order, jet_type, jet_primal, jet_constant, jet_seed
export jet_coefficient, jet_gradient, jet_hessian, seeded_jets
export TaylorGradient, TaylorJacobian

end
