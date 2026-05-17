using StaticArrays: SVector, MVector

abstract type AbstractHOTPSA <: Number end

@inline _hotpsa_nterms(nvars::Integer, order::Integer) = Int(binomial(nvars + order, order))

"""
    HOTPSA{N,O,T,L} <: AbstractHOTPSA

Fast packed high-order Taylor scalar with `N` variables, truncation order `O`,
coefficient type `T`, and `L = binomial(N + O, O)` stored coefficients.
"""
struct HOTPSA{N,O,T<:Number,L} <: AbstractHOTPSA
    coeffs::SVector{L,T}
end

@inline function hotpsa_type(nvars::Integer, order::Integer; T::Type{<:Number}=Float64)
    nvars > 0 || throw(ArgumentError("nvars must be positive"))
    order >= 0 || throw(ArgumentError("order must be non-negative"))
    return HOTPSA{Int(nvars), Int(order), T, _hotpsa_nterms(nvars, order)}
end

@inline _hotpsa_coefftype(::Type{HOTPSA{N,O,T,L}}) where {N,O,T,L} = T
@inline _hotpsa_nvars(::Type{HOTPSA{N,O,T,L}}) where {N,O,T,L} = N
@inline _hotpsa_order(::Type{HOTPSA{N,O,T,L}}) where {N,O,T,L} = O
@inline _hotpsa_length(::Type{HOTPSA{N,O,T,L}}) where {N,O,T,L} = L

@inline function _hotpsa_constant_coeffs(::Type{T}, ::Val{L}, value::Number) where {T,L}
    return SVector{L,T}(ntuple(i -> i == 1 ? convert(T, value) : zero(T), L))
end

@inline function _hotpsa_zero_mvector(::Type{T}, ::Val{L}) where {T,L}
    return MVector{L,T}(ntuple(_ -> zero(T), L))
end

@inline function _hotpsa_promote_coeffs(coeffs::SVector{L,S}, ::Type{T}) where {L,S,T}
    return map(T, coeffs)
end

@inline function HOTPSA{N,O,T,L}(value::Number) where {N,O,T<:Number,L}
    return HOTPSA{N,O,T,L}(_hotpsa_constant_coeffs(T, Val(L), value))
end

@inline function HOTPSA{N,O,T,L}(value::Number, i::Integer) where {N,O,T<:Number,L}
    return TaylorCore.jet_seed(HOTPSA{N,O,T,L}, value, i)
end

Base.convert(::Type{HOTPSA{N,O,T,L}}, x::Number) where {N,O,T<:Number,L} = HOTPSA{N,O,T,L}(x)
Base.convert(::Type{HOTPSA{N,O,T,L}}, x::HOTPSA{N,O,S,L}) where {N,O,T<:Number,S<:Number,L} =
    HOTPSA{N,O,T,L}(_hotpsa_promote_coeffs(x.coeffs, T))

Base.broadcastable(x::HOTPSA{N,O,T,L}) where {N,O,T,L} = Ref(x)
Base.zero(::Type{HOTPSA{N,O,T,L}}) where {N,O,T,L} = HOTPSA{N,O,T,L}(zero(T))
Base.one(::Type{HOTPSA{N,O,T,L}}) where {N,O,T,L} = HOTPSA{N,O,T,L}(one(T))
Base.iszero(x::HOTPSA{N,O,T,L}) where {N,O,T,L} = iszero(x.coeffs[1])
Base.isequal(a::HOTPSA{N,O,T,L}, b::HOTPSA{N,O,T,L}) where {N,O,T,L} = a.coeffs == b.coeffs
Base.:(==)(a::HOTPSA{N,O,T,L}, b::HOTPSA{N,O,T,L}) where {N,O,T,L} = a.coeffs == b.coeffs
Base.show(io::IO, x::HOTPSA{N,O,T,L}) where {N,O,T,L} = print(io, "HOTPSA{$N,$O}(", x.coeffs, ")")

Base.promote_rule(::Type{HOTPSA{N,O,T1,L}}, ::Type{HOTPSA{N,O,T2,L}}) where {N,O,T1<:Number,T2<:Number,L} =
    HOTPSA{N,O,promote_type(T1, T2),L}
Base.promote_rule(::Type{HOTPSA{N,O,T1,L}}, ::Type{T2}) where {N,O,T1<:Number,T2<:Number,L} =
    HOTPSA{N,O,promote_type(T1, T2),L}

@inline function _hotpsa_product_table(::Val{N}, ::Val{O}) where {N,O}
    key = (N, O)
    cache = get!(_HOTPSA_PRODUCT_CACHE, key) do
        table = monomial_table(Val(N), Val(O))
        L = length(table)
        prod = Matrix{Int}(undef, L, L)
        @inbounds for i in 1:L
            for j in 1:L
                prod[i, j] = monomial_product_index(table, i, j)
            end
        end
        prod
    end
    return cache::Matrix{Int}
end

const _HOTPSA_PRODUCT_CACHE = Dict{Tuple{Int,Int},Any}()

@inline function _hotpsa_mul_coeffs(a::SVector{L,T1}, b::SVector{L,T2}, ::Val{N}, ::Val{O}) where {L,T1<:Number,T2<:Number,N,O}
    T = promote_type(T1, T2)
    prod = _hotpsa_product_table(Val(N), Val(O))
    acc = _hotpsa_zero_mvector(T, Val(L))
    @inbounds for i in 1:L
        ai = a[i]
        iszero(ai) && continue
        for j in 1:L
            bj = b[j]
            iszero(bj) && continue
            idx = prod[i, j]
            idx == 0 && continue
            acc[idx] += ai * bj
        end
    end
    return SVector{L,T}(acc)
end

@inline function _hotpsa_nilpotent(x::HOTPSA{N,O,T,L}) where {N,O,T,L}
    coeffs = MVector{L,T}(x.coeffs)
    coeffs[1] = zero(T)
    return HOTPSA{N,O,T,L}(SVector{L,T}(coeffs))
end

function _hotpsa_compose(x::HOTPSA{N,O,T,L}, coeffs::AbstractVector{S}) where {N,O,T<:Number,L,S<:Number}
    RT = promote_type(T, S)
    RType = HOTPSA{N,O,RT,L}
    u = convert(RType, _hotpsa_nilpotent(x))
    result = RType(coeffs[end])
    @inbounds for k in (length(coeffs) - 1):-1:1
        result = coeffs[k] + u * result
    end
    return result
end

Base.:-(x::HOTPSA{N,O,T,L}) where {N,O,T,L} = HOTPSA{N,O,T,L}(-x.coeffs)

function Base.:+(a::HOTPSA{N,O,T1,L}, b::HOTPSA{N,O,T2,L}) where {N,O,T1<:Number,T2<:Number,L}
    T = promote_type(T1, T2)
    return HOTPSA{N,O,T,L}(_hotpsa_promote_coeffs(a.coeffs, T) + _hotpsa_promote_coeffs(b.coeffs, T))
end

function Base.:+(a::HOTPSA{N,O,T,L}, b::Number) where {N,O,T<:Number,L}
    coeffs = MVector{L,T}(a.coeffs)
    coeffs[1] += convert(T, b)
    return HOTPSA{N,O,T,L}(SVector{L,T}(coeffs))
end

Base.:+(a::Number, b::HOTPSA{N,O,T,L}) where {N,O,T<:Number,L} = b + a

function Base.:-(a::HOTPSA{N,O,T1,L}, b::HOTPSA{N,O,T2,L}) where {N,O,T1<:Number,T2<:Number,L}
    T = promote_type(T1, T2)
    return HOTPSA{N,O,T,L}(_hotpsa_promote_coeffs(a.coeffs, T) - _hotpsa_promote_coeffs(b.coeffs, T))
end

function Base.:-(a::HOTPSA{N,O,T,L}, b::Number) where {N,O,T<:Number,L}
    coeffs = MVector{L,T}(a.coeffs)
    coeffs[1] -= convert(T, b)
    return HOTPSA{N,O,T,L}(SVector{L,T}(coeffs))
end

Base.:-(a::Number, b::HOTPSA{N,O,T,L}) where {N,O,T<:Number,L} = convert(HOTPSA{N,O,T,L}, a) - b

function Base.:*(a::HOTPSA{N,O,T1,L}, b::HOTPSA{N,O,T2,L}) where {N,O,T1<:Number,T2<:Number,L}
    coeffs = _hotpsa_mul_coeffs(a.coeffs, b.coeffs, Val(N), Val(O))
    T = eltype(coeffs)
    return HOTPSA{N,O,T,L}(coeffs)
end

function Base.:*(a::HOTPSA{N,O,T,L}, b::Number) where {N,O,T<:Number,L}
    U = promote_type(T, typeof(b))
    return HOTPSA{N,O,U,L}(_hotpsa_promote_coeffs(a.coeffs, U) .* convert(U, b))
end

Base.:*(a::Number, b::HOTPSA{N,O,T,L}) where {N,O,T<:Number,L} = b * a

Base.inv(x::HOTPSA{N,O,T,L}) where {N,O,T<:Number,L} = begin
    x0 = x.coeffs[1]
    iszero(x0) && throw(DomainError(x0, "inverse is undefined for zero primal value"))
    coeffs = Vector{typeof(inv(x0))}(undef, O + 1)
    coeffs[1] = inv(x0)
    @inbounds for k in 1:O
        coeffs[k + 1] = -coeffs[k] / x0
    end
    _hotpsa_compose(x, coeffs)
end

Base.:/(a::HOTPSA{N,O,T1,L}, b::HOTPSA{N,O,T2,L}) where {N,O,T1<:Number,T2<:Number,L} = a * inv(b)
Base.:/(a::HOTPSA{N,O,T,L}, b::Number) where {N,O,T<:Number,L} = a * inv(convert(promote_type(T, typeof(b)), b))
Base.:/(a::Number, b::HOTPSA{N,O,T,L}) where {N,O,T<:Number,L} = convert(HOTPSA{N,O,promote_type(typeof(a), T),L}, a) / b

function Base.:^(x::HOTPSA{N,O,T,L}, p::Integer) where {N,O,T<:Number,L}
    p < 0 && return inv(x^(-p))
    p == 0 && return one(HOTPSA{N,O,T,L})
    p == 1 && return x
    result = one(HOTPSA{N,O,T,L})
    base = x
    n = p
    while n > 0
        if isodd(n)
            result *= base
        end
        n >>= 1
        n == 0 && break
        base *= base
    end
    return result
end

Base.:^(x::HOTPSA{N,O,T,L}, p::Real) where {N,O,T<:Number,L} = exp(p * log(x))

Base.exp(x::HOTPSA{N,O,T,L}) where {N,O,T<:Number,L} = begin
    x0 = x.coeffs[1]
    coeffs = Vector{typeof(exp(x0))}(undef, O + 1)
    coeffs[1] = exp(x0)
    @inbounds for k in 1:O
        coeffs[k + 1] = coeffs[k] / k
    end
    _hotpsa_compose(x, coeffs)
end

Base.log(x::HOTPSA{N,O,T,L}) where {N,O,T<:Number,L} = begin
    x0 = x.coeffs[1]
    if abs(x0) < eps(real(x0))
        throw(DomainError(x0, "log argument too small"))
    end
    invx = inv(x0)
    coeffs = Vector{typeof(log(x0))}(undef, O + 1)
    coeffs[1] = log(x0)
    invxpow = invx
    sgn = one(invx)
    @inbounds for k in 1:O
        coeffs[k + 1] = sgn * invxpow / k
        invxpow *= invx
        sgn = -sgn
    end
    _hotpsa_compose(x, coeffs)
end

Base.log1p(x::HOTPSA{N,O,T,L}) where {N,O,T<:Number,L} = begin
    x0 = x.coeffs[1]
    base = one(x0) + x0
    if abs(base) < eps(real(base))
        throw(DomainError(x0, "log1p argument too small"))
    end
    invbase = inv(base)
    coeffs = Vector{typeof(log1p(x0))}(undef, O + 1)
    coeffs[1] = log1p(x0)
    invpow = invbase
    sgn = one(invbase)
    @inbounds for k in 1:O
        coeffs[k + 1] = sgn * invpow / k
        invpow *= invbase
        sgn = -sgn
    end
    _hotpsa_compose(x, coeffs)
end

Base.sqrt(x::HOTPSA{N,O,T,L}) where {N,O,T<:Number,L} = begin
    x0 = x.coeffs[1]
    if (x0 isa Real) && x0 < zero(x0)
        throw(DomainError(x0, "sqrt of negative number"))
    end
    if iszero(x0)
        nil = _hotpsa_nilpotent(x)
        any(c -> !iszero(c), nil.coeffs) && throw(DomainError(x0, "sqrt derivative undefined at zero with non-zero higher-order coefficients"))
        return HOTPSA{N,O,T,L}(zero(x0))
    end
    coeffs = Vector{typeof(sqrt(x0))}(undef, O + 1)
    coeffs[1] = sqrt(x0)
    invx = inv(x0)
    @inbounds for k in 1:O
        coeffs[k + 1] = coeffs[k] * ((3 / 2 - k) / k) * invx
    end
    _hotpsa_compose(x, coeffs)
end

function _hotpsa_sin_cos_coeffs(x0, order::Int)
    s_coeffs = Vector{typeof(sin(x0))}(undef, order + 1)
    c_coeffs = Vector{typeof(cos(x0))}(undef, order + 1)
    s_coeffs[1] = sin(x0)
    c_coeffs[1] = cos(x0)
    @inbounds for k in 1:order
        s_coeffs[k + 1] = c_coeffs[k] / k
        c_coeffs[k + 1] = -s_coeffs[k] / k
    end
    return s_coeffs, c_coeffs
end

function _hotpsa_sinh_cosh_coeffs(x0, order::Int)
    s_coeffs = Vector{typeof(sinh(x0))}(undef, order + 1)
    c_coeffs = Vector{typeof(cosh(x0))}(undef, order + 1)
    s_coeffs[1] = sinh(x0)
    c_coeffs[1] = cosh(x0)
    @inbounds for k in 1:order
        s_coeffs[k + 1] = c_coeffs[k] / k
        c_coeffs[k + 1] = s_coeffs[k] / k
    end
    return s_coeffs, c_coeffs
end

function _hotpsa_atan_coeffs(x0, order::Int)
    coeffs = Vector{typeof(atan(x0))}(undef, order + 1)
    coeffs[1] = atan(x0)
    order == 0 && return coeffs

    q0 = one(x0) + x0 * x0
    iszero(q0) && throw(DomainError(x0, "atan derivative is singular at the primal value"))
    q1 = 2 * x0
    q2 = one(x0)

    inv_coeffs = Vector{typeof(inv(q0))}(undef, order)
    inv_coeffs[1] = inv(q0)
    @inbounds for n in 1:(order - 1)
        acc = q1 * inv_coeffs[n]
        if n > 1
            acc += q2 * inv_coeffs[n - 1]
        end
        inv_coeffs[n + 1] = -acc / q0
    end
    @inbounds for n in 1:order
        coeffs[n + 1] = inv_coeffs[n] / n
    end
    return coeffs
end

Base.sin(x::HOTPSA{N,O,T,L}) where {N,O,T<:Number,L} = _hotpsa_compose(x, first(_hotpsa_sin_cos_coeffs(x.coeffs[1], O)))
Base.cos(x::HOTPSA{N,O,T,L}) where {N,O,T<:Number,L} = _hotpsa_compose(x, last(_hotpsa_sin_cos_coeffs(x.coeffs[1], O)))
Base.sinh(x::HOTPSA{N,O,T,L}) where {N,O,T<:Number,L} = _hotpsa_compose(x, first(_hotpsa_sinh_cosh_coeffs(x.coeffs[1], O)))
Base.cosh(x::HOTPSA{N,O,T,L}) where {N,O,T<:Number,L} = _hotpsa_compose(x, last(_hotpsa_sinh_cosh_coeffs(x.coeffs[1], O)))
Base.tan(x::HOTPSA{N,O,T,L}) where {N,O,T<:Number,L} = sin(x) / cos(x)
Base.asin(x::HOTPSA{N,O,T,L}) where {N,O,T<:Number,L} = atan(x / sqrt(one(HOTPSA{N,O,T,L}) - x * x))
Base.acos(x::HOTPSA{N,O,T,L}) where {N,O,T<:Number,L} = (pi / 2) - asin(x)
Base.atan(x::HOTPSA{N,O,T,L}) where {N,O,T<:Number,L} = _hotpsa_compose(x, _hotpsa_atan_coeffs(x.coeffs[1], O))
Base.tanh(x::HOTPSA{N,O,T,L}) where {N,O,T<:Number,L} = sinh(x) / cosh(x)

Base.sign(x::HOTPSA{N,O,T,L}) where {N,O,T<:Real,L} = sign(x.coeffs[1])
Base.abs(x::HOTPSA{N,O,T,L}) where {N,O,T<:Real,L} = begin
    s = sign(x.coeffs[1])
    HOTPSA{N,O,T,L}(s .* x.coeffs)
end

Base.:<(a::HOTPSA{N,O,T,L}, b::HOTPSA{N,O,T,L}) where {N,O,T<:Real,L} = a.coeffs[1] < b.coeffs[1]
Base.:>(a::HOTPSA{N,O,T,L}, b::HOTPSA{N,O,T,L}) where {N,O,T<:Real,L} = a.coeffs[1] > b.coeffs[1]
Base.:<=(a::HOTPSA{N,O,T,L}, b::HOTPSA{N,O,T,L}) where {N,O,T<:Real,L} = a.coeffs[1] <= b.coeffs[1]
Base.:>=(a::HOTPSA{N,O,T,L}, b::HOTPSA{N,O,T,L}) where {N,O,T<:Real,L} = a.coeffs[1] >= b.coeffs[1]
Base.isless(a::HOTPSA{N,O,T,L}, b::HOTPSA{N,O,T,L}) where {N,O,T<:Real,L} = a.coeffs[1] < b.coeffs[1]
Base.:<(a::HOTPSA{N,O,T,L}, b::Real) where {N,O,T<:Real,L} = a.coeffs[1] < b
Base.:>(a::HOTPSA{N,O,T,L}, b::Real) where {N,O,T<:Real,L} = a.coeffs[1] > b
Base.:<=(a::HOTPSA{N,O,T,L}, b::Real) where {N,O,T<:Real,L} = a.coeffs[1] <= b
Base.:>=(a::HOTPSA{N,O,T,L}, b::Real) where {N,O,T<:Real,L} = a.coeffs[1] >= b
Base.isless(a::HOTPSA{N,O,T,L}, b::Real) where {N,O,T<:Real,L} = a.coeffs[1] < b
Base.:<(a::Real, b::HOTPSA{N,O,T,L}) where {N,O,T<:Real,L} = a < b.coeffs[1]
Base.:>(a::Real, b::HOTPSA{N,O,T,L}) where {N,O,T<:Real,L} = a > b.coeffs[1]
Base.:<=(a::Real, b::HOTPSA{N,O,T,L}) where {N,O,T<:Real,L} = a <= b.coeffs[1]
Base.:>=(a::Real, b::HOTPSA{N,O,T,L}) where {N,O,T<:Real,L} = a >= b.coeffs[1]
Base.isless(a::Real, b::HOTPSA{N,O,T,L}) where {N,O,T<:Real,L} = a < b.coeffs[1]

TaylorCore.jet_backend(::Type{HOTPSA{N,O,T,L}}) where {N,O,T<:Number,L} = TaylorCore.HOTPSABackend()
TaylorCore.jet_backend(::HOTPSA{N,O,T,L}) where {N,O,T<:Number,L} = TaylorCore.HOTPSABackend()
TaylorCore.jet_nvars(::Type{HOTPSA{N,O,T,L}}) where {N,O,T<:Number,L} = N
TaylorCore.jet_nvars(::HOTPSA{N,O,T,L}) where {N,O,T<:Number,L} = N
TaylorCore.jet_order(::Type{HOTPSA{N,O,T,L}}) where {N,O,T<:Number,L} = O
TaylorCore.jet_order(::HOTPSA{N,O,T,L}) where {N,O,T<:Number,L} = O
TaylorCore.jet_primal(x::HOTPSA{N,O,T,L}) where {N,O,T<:Number,L} = x.coeffs[1]
TaylorCore.monomial_table(::Type{HOTPSA{N,O,T,L}}) where {N,O,T<:Number,L} = monomial_table(Val(N), Val(O))

function TaylorCore.jet_type(::TaylorCore.HOTPSABackend, nvars::Integer; order::Integer, T::Type{<:Number}=Float64)
    return hotpsa_type(nvars, order; T=T)
end

function TaylorCore.jet_seed(::Type{HOTPSA{N,O,T,L}}, value::Real, i::Integer) where {N,O,T<:Number,L}
    O >= 1 || throw(ArgumentError("HOTPSA seed requires order >= 1"))
    1 <= i <= N || throw(ArgumentError("seed index $i is outside 1:$N"))
    coeffs = _hotpsa_zero_mvector(T, Val(L))
    coeffs[1] = convert(T, value)
    table = monomial_table(Val(N), Val(O))
    exponent = ntuple(j -> j == i ? 1 : 0, N)
    idx = monomial_index(table, exponent)
    idx == 0 && throw(ArgumentError("monomial index for seed variable $i is missing"))
    coeffs[idx] = one(T)
    return HOTPSA{N,O,T,L}(SVector{L,T}(coeffs))
end

function TaylorCore.jet_coefficient(x::HOTPSA{N,O,T,L}, exponent::NTuple{N,Int}) where {N,O,T<:Number,L}
    sum(exponent) <= O || return zero(T)
    table = monomial_table(Val(N), Val(O))
    idx = monomial_index(table, exponent)
    return idx == 0 ? zero(T) : x.coeffs[idx]
end
