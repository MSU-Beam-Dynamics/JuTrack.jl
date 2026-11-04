"""
    WigglerRadiation(ring, lindata; energy=3e9, nstep=60, details=false)

Radiation integrals for aligned wigglers. `lindata` contains element-exit
Twiss data. Native `WIGGLER.lw` and flattened six-entry harmonic blocks are
supported, as are legacy `Lw` and 6-by-N harmonic matrices. `By` bends in x;
`Bx` bends in y. Harmonic entries are [id, amplitude, kx/kw, ky/kw, kz/kw, phase].

The default return remains `(I1, I2, I3, I4, I5)` with horizontal I5. With
`details=true`, a named tuple additionally supplies `I5y`, `cone_x`, `cone_y`
and the model description. Cone integrals must be multiplied by 13*Cq/55
and divided by the appropriate damping integral to obtain emittances.

I2/I3 use the squared/absolute-cubed on-axis *curvature*, including harmonic
interference. I1 uses the periodic oscillatory dispersion. I4 is neglected.
As in the earlier AT-style helper, I5 is the maximum of the transported
entrance-H contribution and intrinsic wiggler excitation; their cross terms
are not modeled. Intrinsic excitation and photon-cone integrals use uncoupled
averaged planar focusing. This is not a coupled radiation-diffusion solution.
Offsets, rolls and non-planar transverse harmonic shapes are rejected.

`nstep` is the minimum quadrature intervals per shortest harmonic period.
Integration uses the declared length; native tracking currently rounds a
noninteger number of wiggler periods, so use integer-period native devices
when comparing with tracking. Radiation flags do not affect these integrals.
"""
function WigglerRadiation(ring::Vector{<:AbstractElement{T}},
                         lindata::Vector{<:EdwardsTengTwiss{T}};
                         energy::Float64=3e9, nstep::Int=60,
                         details::Bool=false) where T
    length(lindata) == length(ring) || error("lindata must have the same length as ring")
    isfinite(energy) && energy > 0 || throw(ArgumentError("energy must be finite and positive"))
    nstep > 0 || throw(ArgumentError("nstep must be positive"))
    totals = zeros(T, 8)
    for i in eachindex(ring)
        elem = ring[i]
        iswig = elem isa WIGGLER || (hasproperty(elem, :eletype) &&
            occursin("WIGGLER", uppercase(elem.eletype)))
        iswig || continue
        dini = lindata[i == 1 ? end : i-1]
        vals = _wiggler_integrals(elem, dini, energy / speed_of_light, nstep)
        totals .+= vals
    end
    basic = Tuple(totals[1:5])
    details || return basic
    return (I1=totals[1], I2=totals[2], I3=totals[3], I4=totals[4],
            I5=totals[5], I5y=totals[6], cone_x=totals[7], cone_y=totals[8],
            model="Aligned on-axis curvature; uncoupled averaged planar focusing; I4 and external/intrinsic dispersion cross terms neglected")
end

function _wiggler_period(elem)
    period = hasproperty(elem, :lw) ? elem.lw :
             hasproperty(elem, :Lw) ? elem.Lw :
             throw(ArgumentError("Wiggler has neither lw nor Lw"))
    isfinite(period) && period > 0 || throw(ArgumentError("Wiggler period must be finite and positive"))
    return period
end

function _wiggler_harmonics(elem, field::Symbol)
    hasproperty(elem, field) || return zeros(6, 0)
    raw = getproperty(elem, field)
    if raw isa AbstractVector
        length(raw) % 6 == 0 || throw(ArgumentError("$field must contain complete six-entry harmonic blocks"))
        harmonics = reshape(raw, 6, :)
    elseif raw isa AbstractMatrix && size(raw, 1) == 6
        harmonics = raw
    else
        throw(ArgumentError("$field must be a flattened vector or a 6-by-N matrix"))
    end
    all(isfinite, harmonics) || throw(ArgumentError("Nonfinite wiggler harmonic data"))
    all(>(0), harmonics[5, :]) || throw(ArgumentError("Longitudinal harmonic numbers must be positive"))
    return harmonics
end

# Combine coincident frequencies before squaring: opposite harmonics cancel.
function _wiggler_modes(harmonics, rhoinv)
    modes = Dict{Float64,ComplexF64}()
    for j in axes(harmonics, 2)
        n = Float64(harmonics[5, j])
        modes[n] = get(modes, n, 0.0im) + rhoinv * harmonics[2, j] * cis(harmonics[6, j])
    end
    return [(n, c) for (n, c) in sort!(collect(modes); by=first) if abs(c) > 0]
end

"""On-axis transverse curvatures `(Bx/Brho, By/Brho)`, in inverse metres."""
function Baxis(elem::AbstractElement, s::AbstractVector, kw::Float64, rhoinv::Float64)
    bx = zeros(length(s)); by = zeros(length(s))
    for (n, c) in _wiggler_modes(_wiggler_harmonics(elem, :Bx), rhoinv)
        bx .+= real.(c .* cis.(n * kw .* s))
    end
    for (n, c) in _wiggler_modes(_wiggler_harmonics(elem, :By), rhoinv)
        by .-= real.(c .* cis.(n * kw .* s))
    end
    return bx, by
end

function _wiggler_simpson(values, length)
    intervals = size(values, 1)-1
    result = values[1, :] + values[end, :]
    for j in 2:intervals
        result += (iseven(j) ? 4 : 2) * values[j, :]
    end
    return result * (length / (3 * intervals))
end

function _wiggler_twiss(beta, alpha, focusing, s)
    beta > 0 || throw(ArgumentError("Wiggler entrance beta must be positive"))
    gamma = (1+alpha^2)/beta
    if focusing > 0
        q = sqrt(focusing)
        a, b, c, d = cos(q*s), sin(q*s)/q, -q*sin(q*s), cos(q*s)
    else
        a, b, c, d = 1.0, s, 0.0, 1.0
    end
    betas = a*a*beta - 2*a*b*alpha + b*b*gamma
    alphas = -a*c*beta + (a*d+b*c)*alpha - b*d*gamma
    return betas, alphas, (1+alphas^2)/betas
end

function _wiggler_intrinsic(modes, kw, s)
    dispersion = 0.0; slope = 0.0; periodic = 0.0; curvature = 0.0
    for (n, c) in modes
        k = n*kw
        wave = c*cis(k*s)
        curvature += real(wave)
        periodic -= real(wave)/(k*k)
        # Zero entrance dispersion/slope. Its overall sign does not affect
        # intrinsic H; no interference with external dispersion is claimed.
        dispersion += (real(c)-real(wave)-k*s*imag(c))/(k*k)
        slope += (imag(wave)-imag(c))/k
    end
    return dispersion, slope, periodic, curvature
end

function _wiggler_integrals(elem, dini, Brho, nstep)
    isfinite(Brho) && Brho > 0 || throw(ArgumentError("Brho must be finite and positive"))
    nstep > 0 || throw(ArgumentError("nstep must be positive"))
    period = _wiggler_period(elem)
    le = elem.len
    isfinite(le) && le >= 0 || throw(ArgumentError("Wiggler length must be finite and nonnegative"))
    isfinite(elem.Bmax) || throw(ArgumentError("Wiggler field must be finite"))
    by = _wiggler_harmonics(elem, :By); bx = _wiggler_harmonics(elem, :Bx)
    for (harmonics, flat_coordinate, varying_coordinate) in ((by, 3, 4), (bx, 4, 3))
        for j in axes(harmonics, 2)
            harmonics[2, j] == 0 && continue
            (harmonics[flat_coordinate, j] == 0 &&
             abs(harmonics[varying_coordinate, j]) == harmonics[5, j]) ||
                throw(ArgumentError("Wiggler radiation approximation requires standard planar transverse harmonics"))
        end
    end
    for field in (:T1, :T2)
        hasproperty(elem, field) && any(!iszero, getproperty(elem, field)) &&
            throw(ArgumentError("Offset wigglers require an orbit-dependent radiation calculation"))
    end
    for field in (:R1, :R2)
        if hasproperty(elem, field)
            transform = getproperty(elem, field)
            (all(iszero, transform) || transform == Matrix{eltype(transform)}(I, 6, 6)) ||
                throw(ArgumentError("Rotated wigglers require a coupled radiation calculation"))
        end
    end
    mx = _wiggler_modes(by, elem.Bmax/Brho) # horizontal bending
    my = _wiggler_modes(bx, elem.Bmax/Brho) # vertical bending
    (le == 0 || (isempty(mx) && isempty(my))) && return zeros(8)
    kw = 2*pi/period
    highest = maximum(first, vcat(mx, my))
    intervals = max(2, ceil(Int, nstep * highest * le/period))
    intervals += isodd(intervals)
    # Leading planar focusing, averaged over complete harmonic periods.
    focusing_x = sum(abs2(c) for (_, c) in my; init=0.0)/2
    focusing_y = sum(abs2(c) for (_, c) in mx; init=0.0)/2
    integrands = zeros(intervals+1, 8)
    for j in 0:intervals
        s = le*j/intervals
        dx, dpx, px, hx = _wiggler_intrinsic(mx, kw, s)
        dy, dpy, py, hy = _wiggler_intrinsic(my, kw, s)
        curvature2 = hx*hx+hy*hy
        curvature = sqrt(curvature2); abs3 = curvature2*curvature
        betax, alphax, gammax = _wiggler_twiss(dini.betax, dini.alphax, focusing_x, s)
        betay, alphay, gammay = _wiggler_twiss(dini.betay, dini.alphay, focusing_y, s)
        hxdisp = gammax*dx*dx + 2*alphax*dx*dpx + betax*dpx*dpx
        hydisp = gammay*dy*dy + 2*alphay*dy*dpy + betay*dpy*dpy
        integrands[j+1, :] = [px*hx+py*hy, curvature2, abs3, 0.0,
                              hxdisp*abs3, hydisp*abs3,
                              betax*hy*hy*curvature, betay*hx*hx*curvature]
    end
    result = _wiggler_simpson(integrands, le)
    # Exact phase-independent I2/I3 for a single mode over complete periods.
    modes = vcat(mx, my)
    if length(modes) == 1
        n, c = only(modes)
        periods = n*le/period
        if isapprox(periods, round(periods); rtol=1e-12, atol=1e-12)
            result[1] = -le*abs2(c)/(2*(n*kw)^2)
            result[2] = le*abs2(c)/2
            result[3] = 4*le*abs(c)^3/(3*pi)
        end
    end
    Hx = (dini.dx^2+(dini.betax*dini.dpx+dini.alphax*dini.dx)^2)/dini.betax
    Hy = (dini.dy^2+(dini.betay*dini.dpy+dini.alphay*dini.dy)^2)/dini.betay
    result[5] = max(Hx*result[3], result[5])
    result[6] = max(Hy*result[3], result[6])
    return result
end

"""Single-device counterpart of `WigglerRadiation`; preserves the five-value API."""
function wigrad(elem::AbstractElement, dini::EdwardsTengTwiss, Brho::Float64, nstep::Int)
    return Tuple(_wiggler_integrals(elem, dini, Brho, nstep)[1:5])
end
