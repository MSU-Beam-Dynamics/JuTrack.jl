module PureErrorFunctions

export erfcx, erf, erfinv
using ..TPSAadStatic  # Note the .. to access parent module
using ..TPSAadStatic: DTPSAD

# ---------------- Common helpers ----------------
@inline invsqrtpi(::Type{T}) where {T<:AbstractFloat} = one(T) / sqrt(T(pi))
const SQRTPI = sqrt(pi)
const INV_SQRTPI = 1 / SQRTPI

# ---- Cody constants for real erfcx ----
struct ErfConsts{T<:AbstractFloat}
    SQRPI::T; THRESH::T; SIXTEN::T; XINF::T; XNEG::T; XSMALL::T; XBIG::T; XHUGE::T; XMAX::T
end
@inline function consts(::Type{Float64})
    T = Float64
    ErfConsts{T}(invsqrtpi(T), T(0.46875), T(16), floatmax(T),
                 -sqrt(log(floatmax(T)/T(2))),
                 eps(T)/T(2),
                 T(26.543),
                 inv(T(2)*sqrt(eps(T)/T(2))),
                 min(floatmax(T), inv(sqrt(T(pi))*floatmin(T))))
end
@inline function consts(::Type{Float32})
    T = Float32
    ErfConsts{T}(invsqrtpi(T), T(0.46875), T(16), floatmax(T),
                 -sqrt(log(floatmax(T)/T(2))),
                 eps(T)/T(2),
                 T(9.194),
                 inv(T(2)*sqrt(eps(T)/T(2))),
                 min(floatmax(T), inv(sqrt(T(pi))*floatmin(T))))
end

# Cody coefficients (stored as Float64; cast to T on use)
const A_dp = (3.16112374387056560e+00,
              1.13864154151050156e+02,
              3.77485237685302021e+02,
              3.20937758913846947e+03,
              1.85777706184603153e-01)
const B_dp = (2.36012909523441209e+01,
              2.44024637934444173e+02,
              1.28261652607737228e+03,
              2.84423683343917062e+03)
const C_dp = (5.64188496988670089e-01,
              8.88314979438837594e+00,
              6.61191906371416295e+01,
              2.98635138197400131e+02,
              8.81952221241769090e+02,
              1.71204761263407058e+03,
              2.05107837782607147e+03,
              1.23033935479799725e+03,
              2.15311535474403846e-08)
const D_dp = (1.57449261107098347e+01,
              1.17693950891312499e+02,
              5.37181101862009858e+02,
              1.62138957456669019e+03,
              3.29079923573345963e+03,
              4.36261909014324716e+03,
              3.43936767414372164e+03,
              1.23033935480374942e+03)
const P_dp = (3.05326634961232344e-01,
              3.60344899949804439e-01,
              1.25781726111229246e-01,
              1.60837851487422766e-02,
              6.58749161529837803e-04,
              1.63153871373020978e-02)
const Q_dp = (2.56852019228982242e+00,
              1.87295284992346047e+00,
              5.27905102951428412e-01,
              6.05183413124413191e-02,
              2.33520497626869185e-03)

@inline toT(T, tpl) = ntuple(i -> T(tpl[i]), length(tpl))

# =========================================================
# Part A) Real erfcx (Cody/NETLIB branching; Float32/64)
# =========================================================
function erfcx(x::T) where {T<:AbstractFloat}
    C = consts(T)
    A = toT(T, A_dp); B = toT(T, B_dp)
    Cc = toT(T, C_dp); D = toT(T, D_dp)
    P = toT(T, P_dp); Q = toT(T, Q_dp)

    X = x
    Y = abs(X)

    # Region 1: small |x|
    if Y ≤ C.THRESH
        ysq = Y > C.XSMALL ? Y*Y : zero(T)
        xnum = A[5]*ysq
        xden = ysq
        @inbounds for i in 1:3
            xnum = (xnum + A[i]) * ysq
            xden = (xden + B[i]) * ysq
        end
        erf_approx = X * (xnum + A[4]) / (xden + B[4])
        return exp(ysq) * (one(T) - erf_approx)
    end

    # Region 2: moderate |x|
    if Y ≤ T(4)
        xnum = Cc[9]*Y
        xden = Y
        @inbounds for i in 1:7
            xnum = (xnum + Cc[i]) * Y
            xden = (xden + D[i])  * Y
        end
        res = (xnum + Cc[8]) / (xden + D[8])
        if X ≥ 0
            return res
        else
            if X < C.XNEG
                return C.XINF
            end
            ysq = floor(Y*C.SIXTEN)/C.SIXTEN
            del = (X - ysq)*(X + ysq)
            ey2 = exp(ysq*ysq) * exp(del)
            return (ey2 + ey2) - res
        end
    end

    # Region 3: large |x|
    if Y ≥ C.XBIG
        if (X ≥ 0) && (Y ≥ C.XMAX)
            return zero(T)
        elseif (X ≥ 0) && (Y ≥ C.XHUGE)
            return C.SQRPI / Y
        elseif (X < 0) && (X < C.XNEG)
            return C.XINF
        end
    end

    ysq = inv(Y*Y)
    xnum = P[6]*ysq
    xden = ysq
    @inbounds for i in 1:4
        xnum = (xnum + P[i]) * ysq
        xden = (xden + Q[i]) * ysq
    end
    res = ysq * (xnum + P[5]) / (xden + Q[5])
    res = (C.SQRPI - res) / Y
    if X ≥ 0
        return res
    else
        if X < C.XNEG
            return C.XINF
        end
        ysq2 = floor(Y*C.SIXTEN)/C.SIXTEN
        del = (X - ysq2)*(X + ysq2)
        ey2 = exp(ysq2*ysq2) * exp(del)
        return (ey2 + ey2) - res
    end
end

function consts(::Type{DTPSAD{N, T}}) where {N, T<:AbstractFloat}
    return consts(T)
end
function erfcx(x::DTPSAD{N, T}) where {N, T<:AbstractFloat}
    C = consts(T)
    A = toT(T, A_dp); B = toT(T, B_dp)
    Cc = toT(T, C_dp); D = toT(T, D_dp)
    P = toT(T, P_dp); Q = toT(T, Q_dp)

    X = x
    Y = abs(X)

    # Region 1: small |x|
    if Y ≤ C.THRESH
        ysq = Y > C.XSMALL ? Y*Y : zero(T)
        xnum = A[5]*ysq
        xden = ysq
        @inbounds for i in 1:3
            xnum = (xnum + A[i]) * ysq
            xden = (xden + B[i]) * ysq
        end
        erf_approx = X * (xnum + A[4]) / (xden + B[4])
        return exp(ysq) * (one(T) - erf_approx)
    end

    # Region 2: moderate |x|
    if Y ≤ T(4)
        xnum = Cc[9]*Y
        xden = Y
        @inbounds for i in 1:7
            xnum = (xnum + Cc[i]) * Y
            xden = (xden + D[i])  * Y
        end
        res = (xnum + Cc[8]) / (xden + D[8])
        if X ≥ 0
            return res
        else
            if X < C.XNEG
                return C.XINF
            end
            ysq = floor(Y*C.SIXTEN)/C.SIXTEN
            del = (X - ysq)*(X + ysq)
            ey2 = exp(ysq*ysq) * exp(del)
            return (ey2 + ey2) - res
        end
    end

    # Region 3: large |x|
    if Y ≥ C.XBIG
        if (X ≥ 0) && (Y ≥ C.XMAX)
            return zero(T)
        elseif (X ≥ 0) && (Y ≥ C.XHUGE)
            return C.SQRPI / Y
        elseif (X < 0) && (X < C.XNEG)
            return C.XINF
        end
    end

    ysq = inv(Y*Y)
    xnum = P[6]*ysq
    xden = ysq
    @inbounds for i in 1:4
        xnum = (xnum + P[i]) * ysq
        xden = (xden + Q[i]) * ysq
    end
    res = ysq * (xnum + P[5]) / (xden + Q[5])
    res = (C.SQRPI - res) / Y
    if X ≥ 0
        return res
    else
        if X < C.XNEG
            return C.XINF
        end
        ysq2 = floor(Y*C.SIXTEN)/C.SIXTEN
        del = (X - ysq2)*(X + ysq2)
        ey2 = exp(ysq2*ysq2) * exp(del)
        return (ey2 + ey2) - res
    end
end
# =========================================================
# Part B) Complex erfcx via Faddeeva w(z)
#          w(z) = exp(-z^2) erfc(-i z),  erfcx(z) = w(i z)
# =========================================================

# (1) Power series for erf(z) (entire) — used only inside w(z) for small |z|
function _erf_series(z::Complex{T}; tol=eps(T)*4, maxiter::Int=800) where {T<:Real}
    z2 = z*z
    term = (2/SQRTPI) * z
    s = term
    n = 1
    while n ≤ maxiter
        term *= -z2 / n          # step recurrence
        add = term / (2n + 1)
        s += add
        if abs(add) ≤ tol*max(T(1), abs(s))
            break
        end
        n += 1
    end
    return s
end

# (2) Laplace continued fraction for w(z) in Im z ≥ 0 (large |z| / sizeable Im z)
function _wofz_cf(z::Complex{T}; n::Int=40) where {T<:Real}
    f = z
    @inbounds for k in 1:n
        f = z - (k/2) / f
    end
    return (im * INV_SQRTPI) / f
end

# (3) w(z) in upper half-plane: choose between CF and series route
function _wofz_upper(z::Complex{T}) where {T<:Real}
    y = imag(z); x = real(z)
    if abs(z) ≥ 4 || y ≥ 0.1*abs(x)
        return _wofz_cf(z)
    else
        return exp(-z*z) * (one(z) - _erf_series(-im*z))
    end
end

# (4) w(z) everywhere via symmetry: w(z) = 2exp(-z^2) - w(-z) for Im z < 0
function _wofz(z::Complex{T}) where {T<:Real}
    if imag(z) < 0
        return 2*exp(-z*z) - _wofz_upper(-z)
    else
        return _wofz_upper(z)
    end
end

# Public complex erfcx: erfcx(z) = w(i z)
erfcx(z::Complex{T}) where {T<:Real} = _wofz(im*z)

function _erf_series(z::DTPSAD{N, T}) where {N, T<:Complex{<:Real}}
    z2 = z*z
    term = (2/SQRTPI) * z
    s = term
    n = 1
    tol = eps(real(T))*4
    maxiter = 800
    while n ≤ maxiter
        term *= -z2 / n          # step recurrence
        add = term / (2n + 1)
        s += add
        if abs(add) ≤ tol*max(real(T)(1), abs(s))
            break
        end
        n += 1
    end
    return s
end

function _wofz_cf(z::DTPSAD{N, T}) where {N, T<:Complex{<:Real}}
    f = z
    @inbounds for k in 1:40
        f = z - (k/2) / f
    end
    return (1.0im * INV_SQRTPI) / f
end

function _wofz_upper(z::DTPSAD{N, T}) where {N, T<:Complex{<:Real}}
    y = imag(z); x = real(z)
    if abs(z) ≥ 4 || y ≥ 0.1*abs(x)
        return _wofz_cf(z)
    else
        return exp(-z*z) * (one(z) - _erf_series(-1.0im*z))
    end
end

function _wofz(z::DTPSAD{N, T}) where {N, T<:Complex{<:Real}}
    if imag(z) < 0
        return 2*exp(-z*z) - _wofz_upper(-z)
    else
        return _wofz_upper(z)
    end
end
erfcx(z::DTPSAD{N, T}) where {N, T<:Complex{<:Real}} = _wofz(1.0im*z)

# =========================================================
# Part C) erf (via erfcx identity) for Real & Complex
#          erf(z) = 1 - exp(-z^2) * erfcx(z)
# =========================================================
@inline erf(x::Real) = begin
    T = float(x)
    one(T) - exp(-T*T) * erfcx(T)
end
@inline erf(z::Complex) = one(z) - exp(-z*z) * erfcx(z)

# DTPSAD version
@inline erf(x::DTPSAD{N, T}) where {N, T<:AbstractFloat} = begin
    one(x) - exp(-x*x) * erfcx(x)
end
@inline erf(z::DTPSAD{N, T}) where {N, T<:Complex{<:Real}} = begin
    one(z) - exp(-z*z) * erfcx(z)
end

# =========================================================
# Part D) erfinv for Real (Winitzki seed + Halley refinement)
# =========================================================
@inline function _winitzki_erfinv(ax::T) where {T<:AbstractFloat}
    a = T(0.147)
    l = log1p(-ax*ax)            # log(1 - x^2) with accuracy near 0
    m = T(2)/(T(pi)*a) + l/T(2)
    sqrt(max(zero(T), sqrt(m*m - l/a) - m))
end

@inline erfinv(x::Real) = erfinv(float(x))
function erfinv(x::T) where {T<:AbstractFloat}
    if isnan(x); return T(NaN); end
    if x ==  one(T); return  T(Inf); end
    if x == -one(T); return -T(Inf); end
    if x < -one(T) || x > one(T); return T(NaN); end
    if x == zero(T); return zero(T); end

    s = copysign(one(T), x)
    ax = abs(x)
    y  = _winitzki_erfinv(ax) * s

    invsp = invsqrtpi(T)
    niter = T === Float64 ? 3 : 2
    @inbounds for _ in 1:niter
        # f(y) = erf(y) - x; f'(y) = 2/√π * exp(-y^2)
        e = (one(T) - exp(-y*y) * erfcx(y)) - x
        u = T(2) * invsp * exp(-y*y)
        r = e / u
        y -= r / (one(T) + y*r)   # Halley step
    end
    return y
end

@inline function _winitzki_erfinv(ax::DTPSAD{N, T}) where {N, T<:AbstractFloat}
    a = T(0.147)
    l = log1p(-ax*ax)            # log(1 - x^2) with accuracy near 0
    m = T(2)/(T(pi)*a) + l/T(2)
    sqrt(max(zero(T), sqrt(m*m - l/a) - m))
end
function erfinv(x::DTPSAD{N, T}) where {N, T<:AbstractFloat}
    if isnan(x); return DTPSAD{N, T}(NaN); end
    if x ==  one(T); return  DTPSAD{N, T}(Inf); end
    if x == -one(T); return DTPSAD{N, T}(-Inf); end
    if x < -one(T) || x > one(T); return DTPSAD{N, T}(NaN); end
    if x == zero(T); return DTPSAD{N, T}(zero(T)); end

    s = copysign(one(T), x)
    ax = abs(x)
    y  = _winitzki_erfinv(ax) * s

    invsp = invsqrtpi(T)
    niter = T === Float64 ? 3 : 2
    @inbounds for _ in 1:niter
        # f(y) = erf(y) - x; f'(y) = 2/√π * exp(-y^2)
        e = (one(T) - exp(-y*y) * erfcx(y)) - x
        u = T(2) * invsp * exp(-y*y)
        r = e / u
        y -= r / (one(T) + y*r)   # Halley step
    end
    return y
end 

end
