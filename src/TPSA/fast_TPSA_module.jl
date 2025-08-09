module TPSAadStatic
using StaticArrays

const NVAR_REF = Ref(6)

"""Return the current number of differentiation variables."""
NVAR() = NVAR_REF[]

"""Set the number of differentiation variables and return the new value."""
function set_tps_dim(n::Integer)
    NVAR_REF[] = n
    return n
end

struct DTPSAD{N,T<:Number}
    val::T
    deriv::SVector{N,T}
end

function DTPSAD(a::Number)
    T = typeof(float(a))
    z = zero(SVector{NVAR(),T})
    return DTPSAD{NVAR(),T}(convert(T, a), z)
end

function DTPSAD(a::Number, i::Integer)
    T = typeof(float(a))
    z = zero(SVector{NVAR(),T})
    # set the i-th derivative to one(T)
    z = setindex(z, one(T), i)
    return DTPSAD{NVAR(),T}(convert(T, a), z)
end

function DTPSAD(a::DTPSAD{N,T}) where {N,T}
    return DTPSAD{N,T}(a.val, a.deriv)
end
function Base.convert(::Type{DTPSAD{N,T}}, b::Number) where {N,T<:Number}
    DTPSAD{N,T}(convert(T, b), zero(SVector{N,T}))
end

Base.:-(x::DTPSAD{N,T}) where {N,T} = DTPSAD{N,T}(-x.val, -x.deriv)

Base.:+(a::DTPSAD{N,T}, b::DTPSAD{N,T}) where {N,T} =
    DTPSAD{N,T}(a.val + b.val, a.deriv + b.deriv)

Base.:+(a::DTPSAD{N,T}, b::Real) where {N,T} =
    DTPSAD{N,T}(a.val + convert(T, b), a.deriv)

Base.:+(a::Real, b::DTPSAD{N,T}) where {N,T} = b + a

Base.:-(a::DTPSAD{N,T}, b::DTPSAD{N,T}) where {N,T} =
    DTPSAD{N,T}(a.val - b.val, a.deriv - b.deriv)

Base.:-(a::DTPSAD{N,T}, b::Real) where {N,T} =
    DTPSAD{N,T}(a.val - convert(T, b), a.deriv)

Base.:-(a::Real, b::DTPSAD{N,T}) where {N,T} =
    DTPSAD{N,T}(convert(T, a) - b.val, -b.deriv)

Base.:*(a::DTPSAD{N,T}, b::DTPSAD{N,T}) where {N,T} = begin
    v = a.val * b.val
    dv = b.val .* a.deriv + a.val .* b.deriv
    return DTPSAD{N,T}(v, dv)
end

Base.:*(a::DTPSAD{N,T}, b::Real) where {N,T} =
    DTPSAD{N,T}(a.val * convert(T, b), a.deriv .* convert(T, b))

Base.:*(a::Real, b::DTPSAD{N,T}) where {N,T} = b * a

Base.:/(a::DTPSAD{N,T}, b::DTPSAD{N,T}) where {N,T} = begin
    # derivative of a/b: (a' b - a b') / b^2
    invb = one(T) / b.val
    v = a.val * invb
    dv = (a.deriv .* b.val .- a.val .* b.deriv) .* (invb * invb)
    return DTPSAD{N,T}(v, dv)
end

Base.:/(a::DTPSAD{N,T}, b::Real) where {N,T} = begin
    invb = one(T) / convert(T, b)
    return DTPSAD{N,T}(a.val * invb, a.deriv .* invb)
end

Base.:/(a::Real, b::DTPSAD{N,T}) where {N,T} = begin
    invb0 = one(T) / b.val
    v = convert(T, a) * invb0
    dv = (-convert(T, a) .* b.deriv) .* (invb0 * invb0)
    return DTPSAD{N,T}(v, dv)
end

function Base.:^(x::DTPSAD{N,T}, b::Real) where {N,T}
    b == 0 && return DTPSAD{N,T}(one(T), zero(SVector{N,T}))
    b == 1 && return x
    x0 = x.val
    val = x0^b
    # derivative: d(x^b) = b * x^(b-1) * x'
    factor = b * x0^(b - 1)
    return DTPSAD{N,T}(val, factor .* x.deriv)
end

Base.exp(x::DTPSAD{N,T}) where {N,T} = begin
    e0 = exp(x.val)
    DTPSAD{N,T}(e0, e0 .* x.deriv)
end

Base.log(x::DTPSAD{N,T}) where {N,T} = begin
    # guard against zero or near zero; works for real and complex
    if abs(x.val) < eps(real(x.val))
        throw(DomainError(x.val, "log argument too small"))
    end
    v = log(x.val)
    dv = x.deriv ./ x.val
    DTPSAD{N,T}(v, dv)
end

Base.sqrt(x::DTPSAD{N,T}) where {N,T} = begin
    # for real negative values error; complex values accepted
    if (x.val isa Real) && x.val < zero(x.val)
        throw(DomainError(x.val, "sqrt of negative number"))
    end
    s0 = sqrt(x.val)

    if iszero(x.val)
        if all(iszero, x.deriv)
            return DTPSAD{N,T}(s0, zero(SVector{N,T}))
        else
            # If any derivative is non-zero, sqrt'(0) would be infinite
            throw(DomainError(x.val, "sqrt derivative undefined at zero with non-zero input derivatives"))
        end
    end
    
    factor = one(T) / (2 * s0)
    DTPSAD{N,T}(s0, factor .* x.deriv)
end

Base.sin(x::DTPSAD{N,T}) where {N,T} = begin
    s0 = sin(x.val)
    c0 = cos(x.val)
    DTPSAD{N,T}(s0, c0 .* x.deriv)
end

Base.cos(x::DTPSAD{N,T}) where {N,T} = begin
    s0 = sin(x.val)
    c0 = cos(x.val)
    DTPSAD{N,T}(c0, (-s0) .* x.deriv)
end

Base.tan(x::DTPSAD{N,T}) where {N,T} = begin
    t0 = tan(x.val)
    sec2 = one(T) + t0 * t0
    DTPSAD{N,T}(t0, sec2 .* x.deriv)
end

Base.asin(x::DTPSAD{N,T}) where {N,T} = begin
    # for real values require |x|<1, complex values allowed
    if (x.val isa Real) && abs(x.val) ≥ one(T)
        throw(DomainError(x.val, "asin domain"))
    end
    v = asin(x.val)
    factor = one(T) / sqrt(one(T) - x.val * x.val)
    DTPSAD{N,T}(v, factor .* x.deriv)
end

Base.acos(x::DTPSAD{N,T}) where {N,T} = begin
    if (x.val isa Real) && abs(x.val) ≥ one(T)
        throw(DomainError(x.val, "acos domain"))
    end
    v = acos(x.val)
    factor = -(one(T) / sqrt(one(T) - x.val * x.val))
    DTPSAD{N,T}(v, factor .* x.deriv)
end

Base.atan(x::DTPSAD{N,T}) where {N,T} = begin
    v = atan(x.val)
    factor = one(T) / (one(T) + x.val * x.val)
    DTPSAD{N,T}(v, factor .* x.deriv)
end

Base.atan(y::DTPSAD{N,T}, x::DTPSAD{N,T}) where {N,T} = begin
    r2 = x.val^2 + y.val^2  
    if r2 == 0
        throw(DomainError((y.val, x.val), "atan2(0,0) is undefined"))
    end
    result_val = atan(y.val, x.val)  
    deriv = (x.val .* y.deriv .- y.val .* x.deriv) ./ r2
    DTPSAD{N,T}(result_val, deriv)
end


Base.sinh(x::DTPSAD{N,T}) where {N,T} = begin
    s0 = sinh(x.val)
    c0 = cosh(x.val)
    DTPSAD{N,T}(s0, c0 .* x.deriv)
end

Base.cosh(x::DTPSAD{N,T}) where {N,T} = begin
    s0 = sinh(x.val)
    c0 = cosh(x.val)
    DTPSAD{N,T}(c0, s0 .* x.deriv)
end

Base.tanh(x::DTPSAD{N,T}) where {N,T} = begin
    t0 = tanh(x.val)
    factor = one(T) - t0 * t0
    DTPSAD{N,T}(t0, factor .* x.deriv)
end

Base.asinh(x::DTPSAD{N,T}) where {N,T} = begin
    v = asinh(x.val)
    factor = one(T) / sqrt(one(T) + x.val * x.val)
    DTPSAD{N,T}(v, factor .* x.deriv)
end

Base.acosh(x::DTPSAD{N,T}) where {N,T} = begin
    if (x.val isa Real) && x.val ≤ one(T)
        throw(DomainError(x.val, "acosh domain"))
    end
    v = acosh(x.val)
    factor = one(T) / sqrt(x.val * x.val - one(T))
    DTPSAD{N,T}(v, factor .* x.deriv)
end

Base.atanh(x::DTPSAD{N,T}) where {N,T} = begin
    if (x.val isa Real) && abs(x.val) ≥ one(T)
        throw(DomainError(x.val, "atanh domain"))
    end
    v = atanh(x.val)
    factor = one(T) / (one(T) - x.val * x.val)
    DTPSAD{N,T}(v, factor .* x.deriv)
end
# matrix with scalar
Base.:*(a::DTPSAD{N,T}, b::Matrix{DTPSAD{N,T}}) where {N,T} = begin
    return [a * x for x in b]
end
Base.:*(a::Matrix{DTPSAD{N,T}}, b::DTPSAD{N,T}) where {N,T} = b*a
Base.:/(a::DTPSAD{N,T}, b::Matrix{DTPSAD{N,T}}) where {N,T} = begin
    return [a / x for x in b]
end
Base.:/(a::Matrix{DTPSAD{N,T}}, b::DTPSAD{N,T}) where {N,T} = begin
    return [x / b for x in a]
end
Base.:+(a::DTPSAD{N,T}, b::Matrix{DTPSAD{N,T}}) where {N,T} = begin
    return [a + x for x in b]
end
Base.:+(a::Matrix{DTPSAD{N,T}}, b::DTPSAD{N,T}) where {N,T} = b + a
Base.:-(a::DTPSAD{N,T}, b::Matrix{DTPSAD{N,T}}) where {N,T} = begin
    return [a - x for x in b]
end
Base.:-(a::Matrix{DTPSAD{N,T}}, b::DTPSAD{N,T}) where {N,T} = begin
    return [x - b for x in a]
end
# Vector with scalar
Base.:*(a::DTPSAD{N,T}, b::Vector{DTPSAD{N,T}}) where {N,T} = begin
    return [a * x for x in b]
end
Base.:*(a::Vector{DTPSAD{N,T}}, b::DTPSAD{N,T}) where {N,T} = b*a
Base.:/(a::DTPSAD{N,T}, b::Vector{DTPSAD{N,T}}) where {N,T} = begin
    return [a / x for x in b]
end
Base.:/(a::Vector{DTPSAD{N,T}}, b::DTPSAD{N,T}) where {N,T} = begin
    return [x / b for x in a]
end
Base.:+(a::DTPSAD{N,T}, b::Vector{DTPSAD{N,T}}) where {N,T} = begin
    return [a + x for x in b]
end
Base.:+(a::Vector{DTPSAD{N,T}}, b::DTPSAD{N,T}) where {N,T} = b + a
Base.:-(a::DTPSAD{N,T}, b::Vector{DTPSAD{N,T}}) where {N,T} = begin
    return [a - x for x in b]
end
Base.:-(a::Vector{DTPSAD{N,T}}, b::DTPSAD{N,T}) where {N,T} = begin
    return [x - b for x in a]
end

function Base.isequal(a::DTPSAD{N,T}, b::DTPSAD{N,T}) where {N,T}
    return (a.val == b.val) && (a.deriv == b.deriv)
end

Base.broadcastable(x::DTPSAD{N,T}) where {N,T} = Ref(x)

Base.iszero(x::DTPSAD{N,T}) where {N,T} = iszero(x.val) #&& all(iszero, x.deriv)

Base.sign(x::DTPSAD{N,T}) where {N,T} = sign(x.val)

Base.abs(x::DTPSAD{N,T}) where {N,T} = begin
    s = sign(x.val)
    DTPSAD{N,T}(convert(T, abs(x.val)), s .* x.deriv)
end

# Promotion rules for mixed real/complex operations
Base.promote_rule(::Type{DTPSAD{N,T1}}, ::Type{DTPSAD{N,T2}}) where {N,T1,T2} = 
    DTPSAD{N,promote_type(T1,T2)}

Base.promote_rule(::Type{DTPSAD{N,T1}}, ::Type{T2}) where {N,T1,T2<:Number} = 
    DTPSAD{N,promote_type(T1,T2)}

# Convert between different DTPSAD types
Base.convert(::Type{DTPSAD{N,T}}, x::DTPSAD{N,S}) where {N,T,S} = 
    DTPSAD{N,T}(convert(T, x.val), convert(SVector{N,T}, x.deriv))

# Mixed-type arithmetic operations
Base.:+(a::DTPSAD{N,T1}, b::DTPSAD{N,T2}) where {N,T1,T2} = begin
    T = promote_type(T1, T2)
    DTPSAD{N,T}(a.val + b.val, a.deriv + b.deriv)
end

Base.:-(a::DTPSAD{N,T1}, b::DTPSAD{N,T2}) where {N,T1,T2} = begin
    T = promote_type(T1, T2)
    DTPSAD{N,T}(a.val - b.val, a.deriv - b.deriv)
end

Base.:*(a::DTPSAD{N,T1}, b::DTPSAD{N,T2}) where {N,T1,T2} = begin
    T = promote_type(T1, T2)
    v = a.val * b.val
    dv = b.val .* a.deriv + a.val .* b.deriv
    DTPSAD{N,T}(v, dv)
end

Base.:/(a::DTPSAD{N,T1}, b::DTPSAD{N,T2}) where {N,T1,T2} = begin
    T = promote_type(T1, T2)
    invb = one(T) / b.val
    v = a.val * invb
    dv = (a.deriv .* b.val .- a.val .* b.deriv) .* (invb * invb)
    DTPSAD{N,T}(v, dv)
end

# Complex conjugate
Base.conj(x::DTPSAD{N,T}) where {N,T} = DTPSAD{N,T}(conj(x.val), conj.(x.deriv))

# Real and imaginary parts
Base.real(x::DTPSAD{N,T}) where {N,T} = DTPSAD{N,real(T)}(real(x.val), real.(x.deriv))
Base.imag(x::DTPSAD{N,T}) where {N,T} = DTPSAD{N,real(T)}(imag(x.val), imag.(x.deriv))

# Absolute value for complex numbers
Base.abs(x::DTPSAD{N,T}) where {N,T<:Complex} = begin
    r = abs(x.val)
    if r == 0
        # Handle zero case to avoid division by zero
        DTPSAD{N,real(T)}(zero(real(T)), zero(SVector{N,real(T)}))
    else
        # |z|' = Re(z*conj(z'))/|z| = Re(z*conj(z'))/|z|
        factor = real(conj(x.val)) / r
        DTPSAD{N,real(T)}(r, factor .* real.(x.deriv) + imag(x.val) / r .* imag.(x.deriv))
    end
end

# Angle/phase for complex numbers
Base.angle(x::DTPSAD{N,T}) where {N,T<:Complex} = begin
    r2 = abs2(x.val)  # |x|^2
    if r2 == 0
        throw(DomainError(x.val, "angle of zero"))
    end
    phase = angle(x.val)
    # d/dz arg(z) = -i/(2*|z|^2) * (conj(z) - z) = Im(conj(z))/|z|^2
    factor_real = -imag(x.val) / r2
    factor_imag = real(x.val) / r2
    deriv_real = factor_real .* real.(x.deriv) + factor_imag .* imag.(x.deriv)
    DTPSAD{N,real(T)}(phase, deriv_real)
end

# Utility function to create complex DTPSAD from real and imaginary parts
function Base.complex(re::DTPSAD{N,T1}, im::DTPSAD{N,T2}) where {N,T1<:Real,T2<:Real}
    T = Complex{promote_type(T1,T2)}
    val = complex(re.val, im.val)
    deriv = complex.(re.deriv, im.deriv)
    DTPSAD{N,T}(val, deriv)
end

# Operations between DTPSAD{N,Float64} and Complex numbers
# Addition
Base.:+(a::DTPSAD{N,T}, b::ComplexF64) where {N,T<:Real} = begin
    DTPSAD{N,ComplexF64}(a.val + b, complex.(a.deriv, zero(T)))
end
Base.:+(a::DTPSAD{N,ComplexF64}, b::ComplexF64) where {N} = begin
    DTPSAD{N,ComplexF64}(a.val + b, a.deriv)
end

Base.:+(a::ComplexF64, b::DTPSAD{N,T}) where {N,T<:Real} = b + a
Base.:+(a::ComplexF64, b::DTPSAD{N,ComplexF64}) where {N} = b + a

# Subtraction
Base.:-(a::DTPSAD{N,T}, b::ComplexF64) where {N,T<:Real} = begin
    DTPSAD{N,ComplexF64}(a.val - b, complex.(a.deriv, zero(T)))
end
Base.:-(a::DTPSAD{N,ComplexF64}, b::ComplexF64) where {N} = begin
    DTPSAD{N,ComplexF64}(a.val - b, a.deriv)
end

Base.:-(a::ComplexF64, b::DTPSAD{N,T}) where {N,T<:Real} = begin
    DTPSAD{N,ComplexF64}(a - b.val, complex.(-b.deriv, zero(T)))
end
Base.:-(a::ComplexF64, b::DTPSAD{N,ComplexF64}) where {N} = begin
    DTPSAD{N,ComplexF64}(a - b.val, -b.deriv)
end

# Multiplication
Base.:*(a::DTPSAD{N,T}, b::ComplexF64) where {N,T<:Real} = begin
    DTPSAD{N,ComplexF64}(a.val * b, b .* complex.(a.deriv, zero(T)))
end
Base.:*(a::DTPSAD{N,ComplexF64}, b::ComplexF64) where {N} = begin
    DTPSAD{N,ComplexF64}(a.val * b, a.deriv * b)
end

Base.:*(a::ComplexF64, b::DTPSAD{N,T}) where {N,T<:Real} = b * a
Base.:*(a::ComplexF64, b::DTPSAD{N,ComplexF64}) where {N} = b * a

# Division
Base.:/(a::DTPSAD{N,T}, b::ComplexF64) where {N,T<:Real} = begin
    invb = one(ComplexF64) / b
    DTPSAD{N,ComplexF64}(a.val * invb, invb .* complex.(a.deriv, zero(T)))
end
Base.:/(a::DTPSAD{N,ComplexF64}, b::ComplexF64) where {N} = begin
    invb = one(ComplexF64) / b
    DTPSAD{N,ComplexF64}(a.val * invb, a.deriv * invb)
end

Base.:/(a::ComplexF64, b::DTPSAD{N,T}) where {N,T<:Real} = begin
    invb_val = one(ComplexF64) / b.val
    val = a * invb_val
    # d/dx (a/x) = -a/x^2 * dx
    factor = -a * invb_val * invb_val
    deriv = factor .* complex.(b.deriv, zero(T))
    DTPSAD{N,ComplexF64}(val, deriv)
end
Base.:/(a::ComplexF64, b::DTPSAD{N,ComplexF64}) where {N} = begin
    invb_val = one(ComplexF64) / b.val
    val = a * invb_val
    # d/dx (a/x) = -a/x^2 * dx
    factor = -a * invb_val * invb_val
    deriv = factor .* b.deriv
    DTPSAD{N,ComplexF64}(val, deriv)
end


# Matrix/Vector operations with complex scalars
Base.:*(a::DTPSAD{N,T}, b::Matrix{ComplexF64}) where {N,T<:Real} = begin
    dtpsad_complex = DTPSAD{N,ComplexF64}(complex(a.val, zero(T)), complex.(a.deriv, zero(T)))
    return [dtpsad_complex * x for x in b]
end

Base.:*(a::Matrix{ComplexF64}, b::DTPSAD{N,T}) where {N,T<:Real} = b * a

Base.:+(a::DTPSAD{N,T}, b::Matrix{ComplexF64}) where {N,T<:Real} = begin
    dtpsad_complex = DTPSAD{N,ComplexF64}(complex(a.val, zero(T)), complex.(a.deriv, zero(T)))
    return [dtpsad_complex + x for x in b]
end

Base.:+(a::Matrix{ComplexF64}, b::DTPSAD{N,T}) where {N,T<:Real} = b + a

Base.:-(a::DTPSAD{N,T}, b::Matrix{ComplexF64}) where {N,T<:Real} = begin
    dtpsad_complex = DTPSAD{N,ComplexF64}(complex(a.val, zero(T)), complex.(a.deriv, zero(T)))
    return [dtpsad_complex - x for x in b]
end

Base.:-(a::Matrix{ComplexF64}, b::DTPSAD{N,T}) where {N,T<:Real} = begin
    dtpsad_complex = DTPSAD{N,ComplexF64}(complex(b.val, zero(T)), complex.(b.deriv, zero(T)))
    return [x - dtpsad_complex for x in a]
end

Base.:/(a::DTPSAD{N,T}, b::Matrix{ComplexF64}) where {N,T<:Real} = begin
    dtpsad_complex = DTPSAD{N,ComplexF64}(complex(a.val, zero(T)), complex.(a.deriv, zero(T)))
    return [dtpsad_complex / x for x in b]
end

Base.:/(a::Matrix{ComplexF64}, b::DTPSAD{N,T}) where {N,T<:Real} = begin
    dtpsad_complex = DTPSAD{N,ComplexF64}(complex(b.val, zero(T)), complex.(b.deriv, zero(T)))
    return [x / dtpsad_complex for x in a]
end

# Vector operations with complex scalars
Base.:*(a::DTPSAD{N,T}, b::Vector{ComplexF64}) where {N,T<:Real} = begin
    dtpsad_complex = DTPSAD{N,ComplexF64}(complex(a.val, zero(T)), complex.(a.deriv, zero(T)))
    return [dtpsad_complex * x for x in b]
end

Base.:*(a::Vector{ComplexF64}, b::DTPSAD{N,T}) where {N,T<:Real} = b * a

Base.:+(a::DTPSAD{N,T}, b::Vector{ComplexF64}) where {N,T<:Real} = begin
    dtpsad_complex = DTPSAD{N,ComplexF64}(complex(a.val, zero(T)), complex.(a.deriv, zero(T)))
    return [dtpsad_complex + x for x in b]
end

Base.:+(a::Vector{ComplexF64}, b::DTPSAD{N,T}) where {N,T<:Real} = b + a

Base.:-(a::DTPSAD{N,T}, b::Vector{ComplexF64}) where {N,T<:Real} = begin
    dtpsad_complex = DTPSAD{N,ComplexF64}(complex(a.val, zero(T)), complex.(a.deriv, zero(T)))
    return [dtpsad_complex - x for x in b]
end

Base.:-(a::Vector{ComplexF64}, b::DTPSAD{N,T}) where {N,T<:Real} = begin
    dtpsad_complex = DTPSAD{N,ComplexF64}(complex(b.val, zero(T)), complex.(b.deriv, zero(T)))
    return [x - dtpsad_complex for x in a]
end

Base.:/(a::DTPSAD{N,T}, b::Vector{ComplexF64}) where {N,T<:Real} = begin
    dtpsad_complex = DTPSAD{N,ComplexF64}(complex(a.val, zero(T)), complex.(a.deriv, zero(T)))
    return [dtpsad_complex / x for x in b]
end

Base.:/(a::Vector{ComplexF64}, b::DTPSAD{N,T}) where {N,T<:Real} = begin
    dtpsad_complex = DTPSAD{N,ComplexF64}(complex(b.val, zero(T)), complex.(b.deriv, zero(T)))
    return [x / dtpsad_complex for x in a]
end

# Allow Complex{DTPSAD} construction
Base.Complex(re::DTPSAD{N,T}, im::DTPSAD{N,T}) where {N,T<:Real} = complex(re, im)
Base.Complex(re::DTPSAD{N,T}) where {N,T<:Real} = complex(re, zero(re))

Base.isnan(x::DTPSAD{N,T}) where {N,T} = isnan(x.val) || any(isnan, x.deriv)

Base.isinf(x::DTPSAD{N,T}) where {N,T} = isinf(x.val) || any(isinf, x.deriv)

Base.length(x::DTPSAD{N,T}) where {N,T} = 1

Base.:<(a::DTPSAD{N,T}, b::DTPSAD{N,T}) where {N,T<:Real} = a.val < b.val
Base.:>(a::DTPSAD{N,T}, b::DTPSAD{N,T}) where {N,T<:Real} = a.val > b.val
Base.:<(a::DTPSAD{N,T}, b::Real) where {N,T<:Real} = a.val < b
Base.:>(a::DTPSAD{N,T}, b::Real) where {N,T<:Real} = a.val > b
Base.:<(a::Real, b::DTPSAD{N,T}) where {N,T<:Real} = a < b.val
Base.:>(a::Real, b::DTPSAD{N,T}) where {N,T<:Real} = a > b.val
Base.isless(a::DTPSAD{N,T}, b::DTPSAD{N,T}) where {N,T<:Real} = a.val < b.val
Base.isless(a::DTPSAD{N,T}, b::Real) where {N,T<:Real} = a.val < b
Base.isless(a::Real, b::DTPSAD{N,T}) where {N,T<:Real} = a < b.val
Base.:<=(a::DTPSAD{N,T}, b::DTPSAD{N,T}) where {N,T<:Real} = a.val ≤ b.val
Base.:>=(a::DTPSAD{N,T}, b::DTPSAD{N,T}) where {N,T<:Real} = a.val ≥ b.val
Base.:<=(a::DTPSAD{N,T}, b::Real) where {N,T<:Real} = a.val ≤ b
Base.:>=(a::DTPSAD{N,T}, b::Real) where {N,T<:Real} = a.val ≥ b
Base.:<=(a::Real, b::DTPSAD{N,T}) where {N,T<:Real} = a ≤ b.val
Base.:>=(a::Real, b::DTPSAD{N,T}) where {N,T<:Real} = a ≥ b.val
Base.:(==)(a::DTPSAD{N,T}, b::DTPSAD{N,T}) where {N,T} = isequal(a, b)
Base.:(==)(a::DTPSAD{N,T}, b::Real) where {N,T} = a.val == b

Base.zero(::Type{DTPSAD{N,T}}) where {N,T} = DTPSAD{N,T}(zero(T), zero(SVector{N,T}))
Base.one(::Type{DTPSAD{N,T}}) where {N,T} = DTPSAD{N,T}(one(T), zero(SVector{N,T}))
Base.zero(x::DTPSAD{N,T}) where {N,T} = DTPSAD{N,T}(zero(T), zero(SVector{N,T}))
Base.one(x::DTPSAD{N,T}) where {N,T}  = DTPSAD{N,T}(one(T), zero(SVector{N,T}))

# 1-D array of DTPSAD like zeros(n) or ones(n)
Base.zeros(::Type{DTPSAD{N,T}}, n::Integer) where {N,T} =
    [DTPSAD{N,T}(zero(T), zero(SVector{N,T})) for _ in 1:n]

Base.ones(::Type{DTPSAD{N,T}}, n::Integer) where {N,T} =
    [DTPSAD{N,T}(one(T), zero(SVector{N,T})) for _ in 1:n]

# 2-D array of DTPSAD like zeros(m, n) or ones(m, n)
Base.zeros(::Type{DTPSAD{N,T}}, m::Integer, n::Integer) where {N,T} =
    [DTPSAD{N,T}(zero(T), zero(SVector{N,T})) for _ in 1:m, _ in 1:n]

Base.ones(::Type{DTPSAD{N,T}}, m::Integer, n::Integer) where {N,T} =
    [DTPSAD{N,T}(one(T), zero(SVector{N,T})) for _ in 1:m, _ in 1:n]

# 3-D array of DTPSAD like zeros(m, n, p) or ones(m, n, p)
Base.zeros(::Type{DTPSAD{N,T}}, m::Integer, n::Integer, p::Integer) where {N,T} =
    [DTPSAD{N,T}(zero(T), zero(SVector{N,T})) for _ in 1:m, _ in 1:n, _ in 1:p]
Base.ones(::Type{DTPSAD{N,T}}, m::Integer, n::Integer, p::Integer) where {N,T} =
    [DTPSAD{N,T}(one(T), zero(SVector{N,T})) for _ in 1:m, _ in 1:n, _ in 1:p]

Base.similar(x::DTPSAD{N,T}) where {N,T} = DTPSAD{N,T}(zero(T), zero(SVector{N,T}))
Base.transpose(x::DTPSAD{N,T}) where {N,T} = x  
# copy
Base.copy(x::DTPSAD{N,T}) where {N,T} = DTPSAD{N,T}(x.val, x.deriv)
Base.copy(x::AbstractArray{<:DTPSAD{N,T}}) where {N,T} = 
    [copy(x[i]) for i in eachindex(x)]  # copy each element
export DTPSAD, set_tps_dim, NVAR


const _DefaultT = Float64
"""
    DTPSAD{N}(a::Number)
"""
function DTPSAD{N}(a::Number) where {N}
    DTPSAD{N,_DefaultT}(convert(_DefaultT, a), zero(SVector{N,_DefaultT}))
end

"""
    DTPSAD{N}(a::Number, i::Integer)

Construct a `DTPSAD` instance of dimension `N` with default element
type `Float64` and a derivative of one at index `i`.
"""
function DTPSAD{N}(a::Number, i::Integer) where {N}
    z = zero(SVector{N,_DefaultT})
    z = setindex(z, one(_DefaultT), i)
    DTPSAD{N,_DefaultT}(convert(_DefaultT, a), z)
end


function Base.convert(::Type{DTPSAD{N}}, b::Number) where {N}
    DTPSAD{N,_DefaultT}(convert(_DefaultT, b), zero(SVector{N,_DefaultT}))
end


Base.zero(::Type{DTPSAD{N}}) where {N} = DTPSAD{N,_DefaultT}(zero(_DefaultT), zero(SVector{N,_DefaultT}))
Base.one(::Type{DTPSAD{N}}) where {N} = DTPSAD{N,_DefaultT}(one(_DefaultT), zero(SVector{N,_DefaultT}))
Base.zero(x::DTPSAD{N}) where {N} = DTPSAD{N,_DefaultT}(zero(_DefaultT), zero(SVector{N,_DefaultT}))
Base.one(x::DTPSAD{N}) where {N}  = DTPSAD{N,_DefaultT}(one(_DefaultT), zero(SVector{N,_DefaultT}))

# Arrays of one-parameter DTPSAD
Base.zeros(::Type{DTPSAD{N}}, n::Integer) where {N} =
    [DTPSAD{N,_DefaultT}(zero(_DefaultT), zero(SVector{N,_DefaultT})) for _ in 1:n]
Base.ones(::Type{DTPSAD{N}}, n::Integer) where {N} =
    [DTPSAD{N,_DefaultT}(one(_DefaultT), zero(SVector{N,_DefaultT})) for _ in 1:n]
Base.zeros(::Type{DTPSAD{N}}, m::Integer, n::Integer) where {N} =
    [DTPSAD{N,_DefaultT}(zero(_DefaultT), zero(SVector{N,_DefaultT})) for _ in 1:m, _ in 1:n]
Base.ones(::Type{DTPSAD{N}}, m::Integer, n::Integer) where {N} =
    [DTPSAD{N,_DefaultT}(one(_DefaultT), zero(SVector{N,_DefaultT})) for _ in 1:m, _ in 1:n]


Base.similar(x::DTPSAD{N}) where {N} = DTPSAD{N,_DefaultT}(zero(_DefaultT), zero(SVector{N,_DefaultT}))

function Gradient(f, x::AbstractVector{<:Number}, Primal::Bool=false)
    n = length(x)
    NVAR() != n && set_tps_dim(n)           # make sure the dual length matches
    # promote input to duals with unit derivative
    X = [DTPSAD(x[i], i) for i in 1:n]
    y = f(X...)                            # splat so users can write f(x1,x2,…) or f(X)
    y_der = y.deriv                       # StaticVector
    if Primal
        return collect(y_der), y.val  # return both derivative and value
    else
        return collect(y_der)           # return only the derivative
    end
end

function Jacobian(F, x::AbstractVector{<:Number}, Primal::Bool=false)
    n = length(x)
    NVAR() != n && set_tps_dim(n)
    X = [DTPSAD(x[i], i) for i in 1:n]

    # call F and coerce result to an AbstractVector of DTPSADs
    Y = F(X...) isa AbstractVector ? F(X...) : collect(F(X...))
    m = length(Y)
    J = Matrix{Float64}(undef, m, n)
    @inbounds for k in 1:m
        J[k, :] .= Y[k].deriv     # each output’s derivative row
    end
    if Primal
        # if Primal is true, return both Jacobian and function value
        return J, [y.val for y in Y]
    else
        return J
    end
end
export Gradient, Jacobian
end # module TPSAadStatic