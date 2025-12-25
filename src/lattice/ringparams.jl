"""
    ElementRadiation(ring::Vector{<:AbstractElement{T}}, lindata::Vector{<:EdwardsTengTwiss{T}}; 
                     UseQuadrupoles::Bool=true)

Compute the radiation integrals in dipoles and quadrupoles.

Analytical integration from:
EVALUATION OF SYNCHROTRON RADIATION INTEGRALS
R.H. Helm, M.J. Lee, P.L. Morton and M. Sands
SLAC-PUB-1193, March 1973

# Arguments
- `ring::Vector{<:AbstractElement}`: Vector of lattice elements
- `lindata::Vector{<:EdwardsTengTwiss}`: Vector of Twiss parameters at element exits (length = length(ring))
- `UseQuadrupoles::Bool=true`: Include quadrupoles in radiation integrals calculation

# Returns
- `I1, I2, I3, I4, I5, I6, Iv`: Seven radiation integrals

# Note
Unlike MATLAB AT, JuTrack's `twissring` returns Twiss parameters at element exits only.
This function handles the indexing difference automatically.

# Example
```julia
twiss = twissring(ring, 0.0, 0)
I1, I2, I3, I4, I5, I6, Iv = ElementRadiation(ring, twiss)
```
"""
function ElementRadiation(ring::Vector{<:AbstractElement{T}}, lindata::Vector{<:EdwardsTengTwiss{T}}; 
                         UseQuadrupoles::Bool=true) where T
    
    # Verify that lindata has the correct length
    if length(lindata) != length(ring)
        error("lindata must have the same length as ring. Got length(lindata)=$(length(lindata)), length(ring)=$(length(ring))")
    end
    
    # Identify dipole elements with non-zero bending angles
    isdipole = Vector{Bool}(undef, length(ring))
    for i in eachindex(ring)
        elem = ring[i]
        if hasproperty(elem, :angle)
            isdipole[i] = abs(elem.angle) > 0.0
        else
            isdipole[i] = false
        end
    end
    
    # Identify quadrupole elements with non-zero focusing
    isquadrupole = fill(false, length(ring))
    if UseQuadrupoles
        for i in eachindex(ring)
            elem = ring[i]
            # Check if element has PolynomB field with quadrupole component (index 2)
            if hasproperty(elem, :PolynomB) && length(elem.PolynomB) >= 2
                isquadrupole[i] = abs(elem.PolynomB[2]) > 0.0
            elseif hasproperty(elem, :k1)
                isquadrupole[i] = abs(elem.k1) > 0.0
            end
        end
    end
    
    # Combine dipole and quadrupole elements
    iselement = isdipole .| isquadrupole
    
    # Get initial and final Twiss parameters for each element
    # In JuTrack, lindata[i] contains Twiss at exit of element i
    # For element i, entrance Twiss is exit of element i-1 (or last element for i=1)
    element_indices = findall(iselement)
    
    vini = Vector{EdwardsTengTwiss{T}}(undef, length(element_indices))
    vend = Vector{EdwardsTengTwiss{T}}(undef, length(element_indices))
    
    for (idx, elem_idx) in enumerate(element_indices)
        # Exit Twiss is straightforward
        vend[idx] = lindata[elem_idx]
        
        # Entrance Twiss is exit of previous element (periodic boundary)
        if elem_idx == 1
            vini[idx] = lindata[end]  # For first element, use last element's exit
        else
            vini[idx] = lindata[elem_idx - 1]
        end
    end
    
    # Compute radiation integrals for each element
    di1_vec = T[]
    di2_vec = T[]
    di3_vec = T[]
    di4_vec = T[]
    di5_vec = T[]
    di6_vec = T[]
    div_vec = T[]
    
    for (idx, elem_idx) in enumerate(element_indices)
        di1, di2, di3, di4, di5, di6, div = elrad(ring[elem_idx], vini[idx], vend[idx])
        push!(di1_vec, di1)
        push!(di2_vec, di2)
        push!(di3_vec, di3)
        push!(di4_vec, di4)
        push!(di5_vec, di5)
        push!(di6_vec, di6)
        push!(div_vec, div)
    end
    
    # Sum up contributions from all elements
    I1 = sum(di1_vec)
    I2 = sum(di2_vec)
    I3 = sum(di3_vec)
    I4 = sum(di4_vec)
    I5 = sum(di5_vec)
    I6 = sum(di6_vec)
    Iv = sum(div_vec)
    
    return I1, I2, I3, I4, I5, I6, Iv
end

"""
    elrad(elem::AbstractElement{T}, dini::EdwardsTengTwiss{T}, dend::EdwardsTengTwiss{T}) where T

Compute radiation integrals for a single element.

# Arguments
- `elem::AbstractElement`: Lattice element
- `dini::EdwardsTengTwiss`: Twiss parameters at element entrance
- `dend::EdwardsTengTwiss`: Twiss parameters at element exit

# Returns
- Tuple of (di1, di2, di3, di4, di5, di6, div): Radiation integral contributions
"""
function elrad(elem::AbstractElement{T}, dini::EdwardsTengTwiss{T}, dend::EdwardsTengTwiss{T}) where T
    betax0 = dini.betax
    betaz0 = dini.betay  # Note: betay in JuTrack corresponds to betaz in MATLAB
    alphax0 = dini.alphax
    eta0 = dini.dx
    etap0 = dini.dpx
    
    # Get bending angle
    if hasproperty(elem, :angle)
        theta = elem.angle
    else
        # For quadrupoles without bending angle, compute from orbit change
        # Note: Assuming closed orbit is stored in R matrix or separate field
        # For now, set theta to 0 for non-bending elements
        theta = T(0.0)
    end
    
    # Skip if angle is too small
    if abs(theta) < T(1.0e-7)
        return T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0), T(0.0)
    end
    
    # Get entrance and exit angles
    ange = hasproperty(elem, :e1) ? elem.e1 : T(0.0)
    ango = hasproperty(elem, :e2) ? elem.e2 : T(0.0)
    
    # Element parameters
    ll = elem.len
    rho = ll / theta
    rho2 = rho * rho
    K = getfoc(elem)
    kx2 = K + 1.0 / rho2
    kz2 = -K
    eps1 = tan(ange) / rho
    eps2 = tan(ango) / rho
    
    # Exit dispersion and alpha
    eta3 = dend.dx
    alphax1 = alphax0 - betax0 * eps1
    gammax1 = (1.0 + alphax1 * alphax1) / betax0
    alphaz1 = dini.alphay + betaz0 * eps1
    alphaz2 = dend.alphay - dend.betay * eps2
    gammaz1 = (1.0 + alphaz1 * alphaz1) / betaz0
    etap1 = etap0 + eta0 * eps1
    etap2 = dend.dpx - eta3 * eps2
    
    # H function at entrance
    h0 = gammax1 * eta0 * eta0 + 2.0 * alphax1 * eta0 * etap1 + betax0 * etap1 * etap1
    
    # Compute averages depending on focusing strength
    if kx2 != 0.0
        if kx2 > 0.0  # Focusing
            kl = ll * sqrt(kx2)
            ss = sin(kl) / kl
            cc = cos(kl)
        else  # Defocusing
            kl = ll * sqrt(-kx2)
            ss = sinh(kl) / kl
            cc = cosh(kl)
        end
        
        eta_ave = (theta - (etap2 - etap1)) / kx2 / ll
        bb = 2.0 * (alphax1 * eta0 + betax0 * etap1) * rho
        aa = -2.0 * (alphax1 * etap1 + gammax1 * eta0) * rho
        h_ave = h0 + (aa * (1.0 - ss) + bb * (1.0 - cc) / ll + 
                     gammax1 * (3.0 - 4.0 * ss + ss * cc) / 2.0 / kx2 - 
                     alphax1 * (1.0 - cc)^2 / kx2 / ll + 
                     betax0 * (1.0 - ss * cc) / 2.0) / kx2 / rho2
    else
        eta_ave = 0.5 * (eta0 + eta3) - ll * ll / 12.0 / rho
        hp0 = 2.0 * (alphax1 * eta0 + betax0 * etap1) / rho
        h2p0 = 2.0 * (-alphax1 * etap1 + betax0 / rho - gammax1 * eta0) / rho
        h_ave = h0 + hp0 * ll / 2.0 + h2p0 * ll * ll / 6.0 - 
                alphax1 * ll^3 / 4.0 / rho2 + 
                gammax1 * ll^4 / 20.0 / rho2
    end
    
    # Average vertical beta function
    if kz2 != 0.0
        bz_ave = (gammaz1 + kz2 * betaz0 + (alphaz2 - alphaz1) / ll) / 2.0 / kz2
    else
        bz_ave = betaz0 - alphaz1 * ll + gammaz1 * ll * ll / 3.0
    end
    
    # Compute radiation integrals
    di1 = eta_ave * ll / rho
    di2 = ll / rho2
    di3 = ll / abs(rho) / rho2
    di4 = eta_ave * ll * (2.0 * K + 1.0 / rho2) / rho - 
          (eta0 * eps1 + eta3 * eps2) / rho
    di5 = h_ave * ll / abs(rho) / rho2
    di6 = kz2^2 * eta_ave^2 * ll
    div = abs(theta) / rho2 * bz_ave
    
    return di1, di2, di3, di4, di5, di6, div
end

"""
    getfoc(elem::AbstractElement{T}) where T

Get the focusing strength (K value) from an element.

# Arguments
- `elem::AbstractElement`: Lattice element

# Returns
- `K::Float64`: Quadrupole focusing strength (K = ∂²By/∂x²)
"""
function getfoc(elem::AbstractElement{T}) where T
    K = T(0.0)
    
    # Check for PolynomB field (normalized multipole component)
    if hasproperty(elem, :PolynomB) && length(elem.PolynomB) >= 2
        K = elem.PolynomB[2]
        
        # Check for consistency if both K and PolynomB are defined
        if hasproperty(elem, :k1) && abs(elem.k1 - K) > 1e-10
            @warn "Values in k1 and PolynomB[2] are different. Using PolynomB[2]"
        end
    elseif hasproperty(elem, :k1)
        K = elem.k1
    end
    
    return K
end

"""
    WigglerRadiation(ring::Vector{<:AbstractElement{T}}, lindata::Vector{<:EdwardsTengTwiss{T}}; 
                     energy::Float64=3e9, nstep::Int=60) where T

Compute the radiation integrals in wigglers.

WigglerRadiation computes the radiation integrals for all wigglers with
the following approximations:
- The self-induced dispersion is neglected in I4 and I5, but is used as
  a lower limit for the I5 contribution
- I1, I2 are integrated analytically
- I3 is integrated analytically for a single harmonic, numerically otherwise

# Arguments
- `ring::Vector{<:AbstractElement}`: Vector of lattice elements
- `lindata::Vector{<:EdwardsTengTwiss}`: Vector of Twiss parameters at element exits
- `energy::Float64=3e9`: Beam energy in eV
- `nstep::Int=60`: Number of integration steps for numerical I3 calculation

# Returns
- `I1, I2, I3, I4, I5`: Five radiation integrals for wigglers

# Note
This function assumes wiggler elements have the following fields:
- `Lw`: Wiggler period length
- `Bmax`: Maximum magnetic field
- `By`: Horizontal field harmonics (matrix where each column is [_, amplitude, _, _, harmonic, phase])
- `Bx`: Vertical field harmonics (matrix where each column is [_, amplitude, _, _, harmonic, phase])

# Example
```julia
twiss = twissring(ring, 0.0, 0)
I1, I2, I3, I4, I5 = WigglerRadiation(ring, twiss, energy=3e9)
```
"""
function WigglerRadiation(ring::Vector{<:AbstractElement{T}}, lindata::Vector{<:EdwardsTengTwiss{T}}; 
                         energy::Float64=3e9, nstep::Int=60) where T
    
    # Verify that lindata has the correct length
    if length(lindata) != length(ring)
        error("lindata must have the same length as ring. Got length(lindata)=$(length(lindata)), length(ring)=$(length(ring))")
    end
    
    # Identify wiggler elements
    # Note: JuTrack may not have a standard Wiggler class yet
    # This checks for elements with "WIGGLER" in eletype
    iswiggler = falses(length(ring))
    for i in eachindex(ring)
        elem = ring[i]
        if hasproperty(elem, :eletype) && occursin("WIGGLER", uppercase(elem.eletype))
            # Additional check: must have wiggler-specific fields
            if hasproperty(elem, :Lw) && hasproperty(elem, :Bmax)
                iswiggler[i] = true
            end
        end
    end
    
    if !any(iswiggler)
        # No wigglers found, return zeros
        return T(0.0), T(0.0), T(0.0), T(0.0), T(0.0)
    end
    
    # Compute magnetic rigidity Brho = p/q = (E/c)/q for relativistic particles
    # Brho [T⋅m] = γβmc / q = √(E² - (mc²)²) / (qc)
    # For highly relativistic particles: Brho ≈ E / (qc) 
    # E in eV, c in m/s, q = e (elementary charge)
    # Brho [T⋅m] = E[eV] / (c[m/s] * e[C]) = E[eV] / (c[m/s] * 1.602e-19[C]) * 1.602e-19[C]
    # Simplified: Brho [T⋅m] = E[eV] / (c[m/s])
    Brho = energy / speed_of_light  # This gives Brho in T⋅m for E in eV
    
    # Get entrance Twiss for each wiggler
    wiggler_indices = findall(iswiggler)
    
    di1_vec = T[]
    di2_vec = T[]
    di3_vec = T[]
    di4_vec = T[]
    di5_vec = T[]
    
    for elem_idx in wiggler_indices
        # Get entrance Twiss (exit of previous element)
        if elem_idx == 1
            dini = lindata[end]
        else
            dini = lindata[elem_idx - 1]
        end
        
        di1, di2, di3, di4, di5 = wigrad(ring[elem_idx], dini, Brho, nstep)
        push!(di1_vec, di1)
        push!(di2_vec, di2)
        push!(di3_vec, di3)
        push!(di4_vec, di4)
        push!(di5_vec, di5)
    end
    
    I1 = sum(di1_vec)
    I2 = sum(di2_vec)
    I3 = sum(di3_vec)
    I4 = sum(di4_vec)
    I5 = sum(di5_vec)
    
    return I1, I2, I3, I4, I5
end

"""
    wigrad(elem::AbstractElement, dini::EdwardsTengTwiss, Brho::Float64, nstep::Int)

Compute radiation integrals for a single wiggler element.

# Arguments
- `elem::AbstractElement`: Wiggler element
- `dini::EdwardsTengTwiss`: Twiss parameters at wiggler entrance
- `Brho::Float64`: Magnetic rigidity (T⋅m)
- `nstep::Int`: Number of integration steps

# Returns
- Tuple of (di1, di2, di3, di4, di5): Radiation integral contributions
"""
function wigrad(elem::AbstractElement, dini::EdwardsTengTwiss, Brho::Float64, nstep::Int)
    le = elem.len
    alphax0 = dini.alphax
    betax0 = dini.betax
    gammax0 = (alphax0 * alphax0 + 1.0) / betax0
    eta0 = dini.dx
    etap0 = dini.dpx
    H0 = gammax0 * eta0 * eta0 + 2.0 * alphax0 * eta0 * etap0 + betax0 * etap0 * etap0
    avebetax = betax0 + alphax0 * le + gammax0 * le * le / 3.0
    
    kw = 2.0 * π / elem.Lw
    rhoinv = elem.Bmax / Brho
    
    # Extract Fourier coefficients
    # By and Bx are matrices where each column represents a harmonic
    # Row 2 contains the amplitude coefficients
    coefh = hasproperty(elem, :By) && size(elem.By, 2) > 0 ? elem.By[2, :] * rhoinv : Float64[]
    coefv = hasproperty(elem, :Bx) && size(elem.Bx, 2) > 0 ? elem.Bx[2, :] * rhoinv : Float64[]
    coef2 = vcat(coefh, coefv)
    
    # Compute I3
    if length(coef2) == 1
        # Analytical I3 for single harmonic
        di3 = coef2[1]^3 * 4.0 * le / 3.0 / π
    else
        # Numerical I3 for multiple harmonics
        s_vals = LinRange(0.0, elem.Lw, nstep + 1)
        bx, bz = Baxis(elem, s_vals, kw, rhoinv)
        B2 = bx .^ 2 .+ bz .^ 2
        rinv = sqrt.(B2)
        # Trapezoidal integration
        di3 = trapz(s_vals, rinv .^ 3) * le / elem.Lw
    end
    
    # Compute I2
    di2 = le * (dot(coefh, coefh) + dot(coefv, coefv)) / 2.0
    
    # Compute I1
    di1 = -di2 / kw / kw
    
    # I4 is zero (dispersion neglected)
    di4 = 0.0
    
    # Compute I5 with lower limit
    if !isempty(coefh)
        d5lim = 4.0 * avebetax * le * coefh[1]^5 / 15.0 / π / kw / kw
    else
        d5lim = 0.0
    end
    di5 = max(H0 * di3, d5lim)
    
    return di1, di2, di3, di4, di5
end

"""
    Baxis(elem::AbstractElement, s::AbstractVector, kw::Float64, rhoinv::Float64)

Compute the magnetic field on the axis of a generic wiggler.

# Arguments
- `elem::AbstractElement`: Wiggler element
- `s::AbstractVector`: Position along the wiggler
- `kw::Float64`: Wiggler wave number (2π/Lw)
- `rhoinv::Float64`: Inverse of magnetic rigidity scaled by Bmax

# Returns
- `bx, bz`: Horizontal and vertical field components (normalized by Bmax)
"""
function Baxis(elem::AbstractElement, s::AbstractVector, kw::Float64, rhoinv::Float64)
    Bmax = elem.Bmax
    kws = kw .* s
    
    # Initialize field arrays
    bx = zeros(length(s))
    bz = zeros(length(s))
    
    # Process horizontal wiggler harmonics (By field)
    if hasproperty(elem, :By) && size(elem.By, 2) > 0
        for i in 1:size(elem.By, 2)
            pb = elem.By[:, i]
            # By field: Bz/Bmax = -By2 * cos(By5*kw*s + By6)
            # pb[2] = amplitude, pb[5] = harmonic number, pb[6] = phase
            bz .-= Bmax * pb[2] * cos.(pb[5] * kws .+ pb[6])
        end
    end
    
    # Process vertical wiggler harmonics (Bx field)
    if hasproperty(elem, :Bx) && size(elem.Bx, 2) > 0
        for i in 1:size(elem.Bx, 2)
            pb = elem.Bx[:, i]
            # Bx field: Bx/Bmax = Bx2 * cos(Bx5*kw*s + Bx6)
            bx .+= Bmax * pb[2] * cos.(pb[5] * kws .+ pb[6])
        end
    end
    
    return bx, bz
end

"""
    trapz(x::AbstractVector, y::AbstractVector)

Trapezoidal numerical integration.

# Arguments
- `x::AbstractVector`: Independent variable (must be sorted)
- `y::AbstractVector`: Dependent variable

# Returns
- Integral of y with respect to x using trapezoidal rule
"""
function trapz(x::AbstractVector, y::AbstractVector)
    if length(x) != length(y)
        error("x and y must have the same length")
    end
    if length(x) < 2
        return 0.0
    end
    
    integral = 0.0
    for i in 1:(length(x)-1)
        integral += (x[i+1] - x[i]) * (y[i+1] + y[i]) / 2.0
    end
    
    return integral
end

"""
    ElossRadiation(ring::Vector{<:AbstractElement{T}}, lindata::Vector{<:EdwardsTengTwiss{T}}) where T

Compute the radiation integrals for energy loss elements.

This function handles special energy loss elements like `SimpleQuantumDiffusion` 
or `EnergyLoss` that model radiation effects without actual bending magnets.

# Arguments
- `ring::Vector{<:AbstractElement}`: Vector of lattice elements
- `lindata::Vector{<:EdwardsTengTwiss}`: Vector of Twiss parameters at element exits

# Returns
- `I1, I2, I3, I4, I5`: Five radiation integrals for energy loss elements

# Note
Currently returns zeros as JuTrack does not yet have energy loss element types.
This function is a placeholder for future implementation.

# Example
```julia
twiss = twissring(ring, 0.0, 0)
I1e, I2e, I3e, I4e, I5e = ElossRadiation(ring, twiss)
```
"""
function ElossRadiation(ring::Vector{<:AbstractElement{T}}, lindata::Vector{<:EdwardsTengTwiss{T}}) where T
    # Verify that lindata has the correct length
    if length(lindata) != length(ring)
        error("lindata must have the same length as ring. Got length(lindata)=$(length(lindata)), length(ring)=$(length(ring))")
    end
    
    # TODO: Implement energy loss element detection when such element types are added to JuTrack
    # For now, return zeros as no energy loss elements are defined
    return T(0.0), T(0.0), T(0.0), T(0.0), T(0.0)
end

"""
    ringpara(ring::Vector{<:AbstractElement}; 
             energy::Float64=3e9, 
             Vrf::Float64=0.0,
             harm::Int=1,
             freq_rf::Float64=0.0,
             dp::Float64=0.0,
             print_summary::Bool=true)

Calculate and optionally print ring parameters including equilibrium emittance, damping times, 
energy spread, and RF-dependent parameters.

This function computes radiation integrals from dipoles, quadrupoles, and wigglers, then 
calculates equilibrium beam parameters based on synchrotron radiation theory.

# Arguments
- `ring::Vector{<:AbstractElement}`: Vector of lattice elements
- `energy::Float64=3e9`: Beam energy in eV
- `Vrf::Float64=0.0`: RF voltage per cell [V]
- `harm::Int=1`: Harmonic number
- `freq_rf::Float64=0.0`: RF frequency in Hz
- `dp::Float64=0.0`: Momentum deviation (for off-momentum calculations)
- `print_summary::Bool=true`: If true, print parameter summary

# Returns
- NamedTuple with all calculated parameters

# Example
```julia
using JuTrack
include("src/demo/SPEAR3/spear3.jl")
ring = spear3()

# Print summary
ringpara(ring, energy=3e9, Vrf=3.2e6, harm=372, freq_rf=476e6)

# Get parameters as structure
params = ringpara(ring, energy=3e9, Vrf=3e6, harm=372, freq_rf=476e6, print_summary=false)
println("Natural emittance: ", params.emittx * 1e9, " nm⋅rad")
```

# Reference
Based on AT's ringpara function. See also:
- H. Wiedemann, "Particle Accelerator Physics"
- M. Sands, "The Physics of Electron Storage Rings"
"""
function ringpara(ring::Vector{<:AbstractElement}; 
                  energy::Float64=3e9, 
                  Vrf::Float64=0.0,
                  harm::Int=1,
                  freq_rf::Float64=0.0,
                  dp::Float64=0.0,
                  print_summary::Bool=false)
    
    energy_GeV = energy * 1.0e-9  # eV to GeV
    gamma = energy_GeV / __E0
    beta = sqrt(1.0 - 1.0/gamma^2)
    
    twiss = twissring(ring, dp, 0)
    circ = sum([elem.len for elem in ring])
    
    # Compute radiation integrals
    I1d, I2d, I3d, I4d, I5d, I6d, Iv = ElementRadiation(ring, twiss)
    I1w, I2w, I3w, I4w, I5w = WigglerRadiation(ring, twiss, energy=energy)
    I1e, I2e, I3e, I4e, I5e = ElossRadiation(ring, twiss)
    
    I1 = I1d + I1w + I1e
    I2 = I2d + I2w + I2e
    I3 = I3d + I3w + I3e
    I4 = I4d + I4w + I4e
    I5 = I5d + I5w + I5e
    
    # Compute basic parameters
    T0 = circ / speed_of_light / beta  # Revolution period [s]
    frev = 1.0 / T0  # Revolution frequency [Hz]
    R = circ / (2.0 * pi)  # Average radius [m]
    
    # Momentum compaction factor
    alphac = I1 / circ
    etac = 1.0 / gamma^2 - alphac  # Slip factor
    
    # Energy loss per turn [eV]
    # U0 = 1.0e9 * CGAMMA / (2.0 * pi) * energy_GeV^4 * I2
    U0 = tracking_U0(ring, energy, m_e)  
    
    # Damping partition numbers
    Jx = 1.0 - I4 / I2
    Jy = 1.0
    Je = 2.0 + I4 / I2
    
    # Natural emittance [m⋅rad]
    emittx = Cq * gamma^2 * I5 / (I2 - I4)
    
    # Natural energy spread (relative)
    sigma_E = gamma * sqrt(Cq * I3 / (2.0 * I2 + I4))
    
    # Damping rates [1/s]
    alpha0 = U0 / 2.0 / T0 / energy
    alphax = Jx * alpha0  # Horizontal damping rate
    alphay = Jy * alpha0  # Vertical damping rate
    alphaE = Je * alpha0  # Energy damping rate
    
    # Damping times [s]
    dampingtime_x = 1.0 / alphax
    dampingtime_y = 1.0 / alphay
    dampingtime_E = 1.0 / alphaE
    
    # Tunes and chromaticity
    # twissring returns twiss parameters at each element exit
    # The last element gives the phase advance after one complete turn
    nux = twiss[end].mux / (2.0 * pi)
    nuy = twiss[end].muy / (2.0 * pi)
    
    chroms = getchrom(ring, dp=dp, energy=energy)
    chromx = chroms[1]
    chromy = chroms[2]
    
    # Vertical emittance from vertical dispersion
    dipindex = findall(i -> hasproperty(ring[i], :angle) && abs(ring[i].angle) > 0.0, 1:length(ring))
    if !isempty(dipindex)
        beta_y = [twiss[i].betay for i in dipindex]
        alpha_y = [twiss[i].alphay for i in dipindex]
        Dy = [twiss[i].dy for i in dipindex]
        Dpy = [twiss[i].dpy for i in dipindex]
        theta = [ring[i].angle for i in dipindex]
        len = [ring[i].len for i in dipindex]
        
        curVavg1 = (1.0 ./ beta_y) .* (Dy.^2 .+ (beta_y .* Dpy .+ alpha_y .* Dy).^2)
        curVavg = sum(curVavg1 .* abs.(theta)) / sum(len)
        emitty_d = Cq * gamma^2 * curVavg / Jy
    else
        emitty_d = 0.0
    end
    
    # Vertical emittance limit from 1/gamma cone opening angle
    emitty_lim = 13.0 / 55.0 * Cq / Jy / I2 * Iv
    
    # RF-dependent parameters
    if Vrf > 0.0 && U0 <= Vrf
        phi_s = pi - asin(U0 / Vrf)  # Synchronous phase [rad]
        nus = sqrt(harm * Vrf * abs(etac * cos(phi_s)) / (2.0 * pi * energy)) / beta  # Synchrotron tune
        
        # RF bucket height (linear approximation)
        delta_max = sqrt(2.0 * U0 / (pi * alphac * harm * energy)) * 
                    sqrt(sqrt((Vrf / U0)^2 - 1.0) - acos(U0 / Vrf))
        
        # Bunch length [m]
        bunchtime = sigma_E * harm * abs(etac) / nus / (2.0 * pi * freq_rf)  # [s]
        bunchlength = beta * speed_of_light * bunchtime  # [m]
    else
        phi_s = NaN
        nus = NaN
        delta_max = NaN
        bunchlength = NaN
    end
    
    # Create results structure
    results = (
        E0 = energy,
        E0_GeV = energy_GeV,
        gamma = gamma,
        beta = beta,
        R = R,
        circumference = circ,
        T0 = T0,
        frev = frev,
        alphac = alphac,
        etac = etac,
        nux = nux,
        nuy = nuy,
        chromx = chromx,
        chromy = chromy,
        I1 = I1,
        I2 = I2,
        I3 = I3,
        I4 = I4,
        I5 = I5,
        I6 = I6d,
        Iv = Iv,
        Jx = Jx,
        Jy = Jy,
        Je = Je,
        U0 = U0,
        sigma_E = sigma_E,
        emittx = emittx,
        emitty_d = emitty_d,
        emitty_lim = emitty_lim,
        dampingtime_x = dampingtime_x,
        dampingtime_y = dampingtime_y,
        dampingtime_E = dampingtime_E,
        dampingalpha_x = alphax,
        dampingalpha_y = alphay,
        dampingalpha_E = alphaE,
        Vrf = Vrf,
        harm = harm,
        freq_rf = freq_rf,
        phi_s = phi_s,
        nus = nus,
        delta_max = delta_max,
        bunchlength = bunchlength
    )
    
    # Print summary if requested
    if print_summary
        println()
        println("   ******** JuTrack Ring Parameter Summary ********")
        println("   Energy:                     $(energy_GeV) GeV")
        println("   Circumference:              $(circ) m")
        println("   Revolution time:            $(T0*1e9) ns ($(frev*1e-6) MHz)")
        println("   Betatron tune H:            $(nux) ($(((nux - floor(nux))/T0)*1e-3) kHz)")
        println("                 V:            $(nuy) ($(((nuy - floor(nuy))/T0)*1e-3) kHz)")
        if chromx != 0.0 || chromy != 0.0
            println("   Chromaticity H:             $(chromx)")
            println("                V:             $(chromy)")
        end
        println("   Synchrotron Integral 1:     $(I1) m")
        println("                        2:     $(I2) m^-1")
        println("                        3:     $(I3) m^-2")
        println("                        4:     $(I4) m^-1")
        println("                        5:     $(I5) m^-1")
        println("   Damping Partition H:        $(Jx)")
        println("                     V:        $(Jy)")
        println("                     E:        $(Je)")
        println("   Radiation Loss:             $(U0/1000.0) keV")
        println("   Natural Energy Spread:      $(sigma_E)")
        println("   Natural Emittance:          $(emittx*1e9) nm⋅rad")
        println("   Radiation Damping H:        $(dampingtime_x*1e3) ms")
        println("                     V:        $(dampingtime_y*1e3) ms")
        println("                     E:        $(dampingtime_E*1e3) ms")
        println("   Momentum Compaction Factor: $(alphac)")
        println("   Slip factor:                $(etac)")
        println()
        if Vrf > 0.0
            println("   Assuming cavities Voltage:  $(Vrf/1e3) kV")
            println("                   Frequency:  $(freq_rf/1e6) MHz")
            println("             Harmonic Number:  $(harm)")
            if !isnan(phi_s)
                println("   Synchronous Phase:          $(phi_s) rad ($(phi_s*180/pi) deg)")
                println("   Linear Energy Acceptance:   $(delta_max*100) %")
                println("   Synchrotron Tune:           $(nus) ($(nus/T0*1e-3) kHz or $(1/nus) turns)")
                println("   Bunch Length:               $(bunchlength*1e3) mm")
            else
                println("   WARNING: U0 > Vrf - RF voltage insufficient!")
            end
            println()
        end
        println("   Emitty from Dy:             $(emitty_d*1e9) nm⋅rad")
        println("   Emitty 1/gamma cone limit:  $(emitty_lim*1e12) pm⋅rad")
        println()
        return results
    else
        return results
    end
end
