using LinearAlgebra
function get_len(ele::AbstractElement)
    return get_len_value(ele.len)
end
function get_len_value(L::Float64)
    return L[1]
end
function get_len_value(L::DTPSAD{N, T}) where {N, T}
    return L
end

"""
    total_length(ring::Vector)
Calculate the total length of the lattice.
# Arguments
- ring::Vector: a vector of beam line elements
# Return
- leng::Float64: total length of the lattice
"""
function total_length(ring::Vector)
    leng = 0.0
    for i in eachindex(ring)
        leng += get_len(ring[i])
    end
    return leng
end

"""
    spos(ring::Vector)

Calculate the s position of each element in the lattice.

# Arguments
- ring::Vector: a vector of beam line elements

# Return
- pos::Vector{Float64}: a vector of s positions
"""
function spos(ring::Vector{AbstractElement{Float64}})
    pos = zeros(length(ring))
    len = 0.0
    for i in eachindex(ring)
        len += get_len(ring[i])
        pos[i] = len
    end
    return pos
end
function spos(ring::Vector{AbstractElement{DTPSAD{N, T}}}) where {N, T}
    pos = zeros(DTPSAD{N, T}, length(ring))
    len = zero(DTPSAD{N, T})
    for i in eachindex(ring)
        len += get_len(ring[i])
        pos[i] = len
    end
    return pos
end

"""
    spos(ring::Vector, idx::Vector)

Calculate the s position of the specified elements in the lattice.

# Arguments
- ring::Vector: a vector of beam line elements
- idx::Vector: a vector of indices of the elements

# Return
- pos::Vector{Float64}: a vector of s positions of the specified elements
"""
function spos(ring::Vector{AbstractElement{Float64}}, idx::Vector{Int})
    pos_all = spos(ring)
    pos = zeros(length(idx))
    for i in 1:length(idx)
        pos[i] = pos_all[idx[i]]
    end
    return pos
end
function spos(ring::Vector{AbstractElement{DTPSAD{N, T}}}, idx::Vector{Int}) where {N, T}
    pos_all = spos(ring)
    pos = zeros(DTPSAD{N, T}, length(idx))
    for i in 1:length(idx)
        pos[i] = pos_all[idx[i]]
    end
    return pos
end
# function findelem(ring::Vector, field::Symbol, value)
#     # Example: findelem(ring, :name, "QF")
#     # warning: may not compatible with autodifferentiation
#     ele_index = []
#     for i in eachindex(ring)
#         if field in fieldnames(typeof(ring[i])) && getfield(ring[i], field) == value
#             push!(ele_index, i)
#         end
#     end
#     return ele_index
# end
""" 
    findelem(ring::Vector, field::Symbol, value)

Find the index of elements with the specified field value in the ring.

# Arguments
- ring::Vector: a vector of beam line elements
- field::Symbol: the field name
- value: the value of the field

# Return
- ele_index::Vector{Int}: a vector of indices of the elements with the specified field value

# Example
```julia
ele_index = findelem(ring, :name, "QF")
```
"""
function findelem(ring::Vector, field::Symbol, value)
    c = 0
    for i in eachindex(ring)
        if field in fieldnames(typeof(ring[i])) && getfield(ring[i], field) == value
            c += 1
        end
    end
    ele_index = zeros(Int, c)
    c = 0
    for i in eachindex(ring)
        if field in fieldnames(typeof(ring[i])) && getfield(ring[i], field) == value
            c += 1
            ele_index[c] = i
        end
    end
    return ele_index
end

# function findelem(ring::Vector, type::Type)
#     # Example: findelem(ring, DRFIT)    
#     # warning: may not compatible with autodifferentiation
#     ele_index = []
#     for i in eachindex(ring)
#         if typeof(ring[i]) == type
#             push!(ele_index, i)
#         end
#     end
#     return ele_index
# end
"""
    findelem(ring::Vector, type::Type)

Find the index of elements with the specified type in the ring.

# Arguments
- ring::Vector: a vector of beam line elements
- type::Type: the type of the element

# Return
- ele_index::Vector{Int}: a vector of indices of the elements with the specified type

# Example
```julia
ele_index = findelem(ring, DRIFT)
```
"""
function findelem(ring::Vector, type::Type)
    c = 0
    for i in eachindex(ring)
        if isa(ring[i], type)
            c += 1
        end
    end
    ele_index = zeros(Int, c)
    c = 0
    for i in eachindex(ring)
        if isa(ring[i], type)
            c += 1
            ele_index[c] = i
        end
    end
    return ele_index
end

function use_exact_drift(flag)
    if flag == 1
        global use_exact_Hamiltonian = 1
    else
        global use_exact_Hamiltonian = 0
    end
end

function array_optics(Twi)
    # use type of Twi[1].betax to determine the type of beta, alpha, gamma, mu, dp
    T = typeof(Twi[1].betax)
    beta = zeros(T, length(Twi), 2)
    beta[:, 1] = [Twi[i].betax for i in eachindex(Twi)]
    beta[:, 2] = [Twi[i].betay for i in eachindex(Twi)]
    alpha = zeros(T, length(Twi), 2)
    alpha[:, 1] = [Twi[i].alphax for i in eachindex(Twi)]
    alpha[:, 2] = [Twi[i].alphay for i in eachindex(Twi)]
    gamma = zeros(T, length(Twi), 2)
    gamma[:, 1] = [Twi[i].gammax for i in eachindex(Twi)]
    gamma[:, 2] = [Twi[i].gammay for i in eachindex(Twi)]
    mu = zeros(T, length(Twi), 2)
    mu[:, 1] = [Twi[i].mux for i in eachindex(Twi)]
    mu[:, 2] = [Twi[i].muy for i in eachindex(Twi)]
    dp = zeros(T, length(Twi), 4)
    dp[:, 1] = [Twi[i].dx for i in eachindex(Twi)]
    dp[:, 2] = [Twi[i].dpx for i in eachindex(Twi)]
    dp[:, 3] = [Twi[i].dy for i in eachindex(Twi)]
    dp[:, 4] = [Twi[i].dpy for i in eachindex(Twi)]
    return beta, alpha, gamma, mu, dp
end


function symplectic(M66::Array{Float64,2})
    # check if a transfer map is symplectic
    # the canonical coordinates are (x, px/p0, y, py/p0, tau, -dE/beta0/c/p0)
    # z = tau, dp/p0 ~ dE/beta0/c/p0. Therefore, the last block of J is -J2
    J2 = [0.0 1.0; -1.0 0.0]
    J = [J2 zeros((2,2)) zeros((2,2)); zeros((2,2)) J2 zeros((2,2)); zeros((2,2)) zeros((2,2)) -J2]
    delta = transpose(M66) * J * M66 .- J
    println("max deviation: ", maximum(abs.(delta)))
    return maximum(abs.(delta))
end

function rad_on!(ring)
    for ele in ring
        if ele isa SBEND || ele isa KQUAD || ele isa KSEXT || ele isa ESBEND
            ele.rad = 1
        end
    end
end
function rad_off!(ring)
    for ele in ring
        if ele isa SBEND || ele isa KQUAD || ele isa KSEXT || ele isa ESBEND
            ele.rad = 0
        end
    end
end

function tracking_U0(ring::Vector{AbstractElement{T}}, energy, mass) where {T}
    rin = zeros(T, 1, 6)
    beam = Beam(rin, energy=energy, mass=mass)
    linepass!(ring, beam)   
    U0 = -beam.r[6] * energy
    return U0
end

function integral_U0(ring, energy, mass)
    cspeed   = speed_of_light
    e_radius = 2.8179403205e-15
    Cgamma   = 4π * e_radius / (3 * m_e^3)      # [m / eV^3]
    Cgamma = Cgamma * (m_e / mass)^4
    Brho     = sqrt(energy^2 - mass^2) / cspeed
    coef     = Cgamma / (2π) * energy^4            

    # 1) Dipole contribution
    dip_idx = findelem(ring, ESBEND)
    θ       = [ring[i].angle for i in dip_idx]
    Ld      = [ring[i].len for i in dip_idx]
    I2d     = sum(abs.(θ .^ 2 ./ Ld))

    # The following code is commented out because no wigglers or energy loss elements are available in JuTrack.
    # 2) Wiggler contribution
    # wig_idx = findall(e -> e.eletype == "Wiggler", ring)
    # I2w     = sum(wiggler_i2(e, Brho) for e in ring[wig_idx])

    # # 3) EnergyLoss‐element contribution
    # el_idx = findall(e -> e.Class == "EnergyLoss" && e.PassMethod != "IdentityPass", ring)
    # I2e    = sum(eloss_i2(e, coef) for e in ring[el_idx])

    # 4) Any user-supplied I2 fields
    # i2_idx = findelem(ring, "I2")
    # I2x    = sum(getfield.(ring[i2_idx], :I2))

    I2_total = I2d #+ I2w + I2e + I2x
    return coef * I2_total
end

wiggler_i2(elem, Brho) = begin
    rhoinv = elem.Bmax / Brho
    hc     = elem.By[2,:] * rhoinv      # horizontal Fourier-B term
    vc     = elem.Bx[2,:] * rhoinv      # vertical Fourier-B term
    return elem.len * (dot(hc,hc) + dot(vc,vc)) / 2
end

"""
    find_closed_orbit_4d(ring::Vector; dp::Float64=0.0, x0=zeros(4), energy::Float64=3.5e9, mass::Float64=m_e,
    tol::Float64=1e-8, maxiter::Int=20)

Find the 4-D closed orbit of a ring using autodiff.
!!! Don't use this function for automatic differentiation because it already uses AD. 
!!! Use this function for AD will result in second order differentiation (slow or crash).
"""
function find_closed_orbit_4d(ring::Vector; dp::Float64=0.0, x0=zeros(4), energy::Float64=3.5e9, mass::Float64=m_e,
    tol::Float64=1e-8, maxiter::Int=20)
    x = copy(x0)

    function fmap(θ::Vector{Float64})
        rin = zeros(Float64, 6)
        rin[1] = θ[1]
        rin[2] = θ[2]
        rin[3] = θ[3]
        rin[4] = θ[4]
        rin[5] = 0.0
        rin[6] = dp
        b = Beam(reshape(rin, 1, 6), energy=energy, mass=mass)          
        ringpass!(ring, b, 1)                               
        return b.r[1, 1:4]                              
    end

    for iter in 1:maxiter
        J, x_out = jacobian(ForwardWithPrimal, Const(fmap), x)     # value and Jacobian
        Δ = x_out .- x                                      # fixed-point error

        if norm(Δ) < tol
            return x                                       # found closed orbit
        end

        # Newton update: solve (I - J) Δx = Δ
        Δx = (I - J[1] + 1e-6 * I) \ Δ                                  # requires invertibility
        x  += Δx
    end
    println("Closed orbit did not converge after $maxiter iterations")
    return x
end

"""
    find_closed_orbit_6d(ring::Vector; x0=zeros(6), energy::Float64=3.5e9, mass::Float64=m_e,
    tol::Float64=1e-8, maxiter::Int=20)

Find the 6-D closed orbit of a ring using autodiff.
To find a 6-D closed orbit, the ring must have active RF cavities.
!!! Don't use this function for automatic differentiation because it already uses AD. 
!!! Use this function for AD will result in second order differentiation (slow or crash).
"""
function find_closed_orbit_6d(ring::Vector; x0=zeros(6), energy::Float64=3.5e9, mass::Float64=m_e,
    tol::Float64=1e-8, maxiter::Int=20)
    x = copy(x0)

    function fmap(θ::Vector{Float64})
        b = Beam(reshape(θ, 1, 6), energy=energy, mass=mass)          
        ringpass!(ring, b, 1)                               
        return b.r[1, :]                              
    end

    for iter in 1:maxiter
        J, x_out = jacobian(ForwardWithPrimal, Const(fmap), x)     # value and Jacobian
        Δ = x_out .- x                                      # fixed-point error

        if norm(Δ) < tol
            return x                                       # found closed orbit
        end

        # Newton update: solve (I - J) Δx = Δ
        Δx = (I - J[1] + 1e-6 * I) \ Δ                                  # requires invertibility
        x  += Δx
    end
    println("Closed orbit did not converge after $maxiter iterations")
    return x
end

function numerical_jacobian(f::Function, x::Vector{Float64}; h::Float64=1e-6)
    n = length(x)
    y0 = f(x)
    m = length(y0)
    J = zeros(m, n)
    for j in 1:n
        dx = zeros(n)
        dx[j] = h
        y_plus  = f(x .+ dx)
        y_minus = f(x .- dx)
        @inbounds J[:, j] = (y_plus .- y_minus) ./ (2h)
    end
    return J, y0
end

# find_closed_orbit using finite differences
function fast_closed_orbit_6d(ring::Vector; x0=zeros(6), energy::Float64=3.5e9,
                                 mass::Float64=m_e, tol::Float64=1e-8, maxiter::Int=20)
    x = copy(x0)

    function fmap(θ::Vector{Float64})
        b = Beam(reshape(θ, 1, 6), energy=energy, mass=mass)
        ringpass!(ring, b, 1)
        return b.r[1, :]
    end

    for iter in 1:maxiter
        J, x_out = numerical_jacobian(fmap, x)    # replace AD with FD
        Δ = x_out .- x

        if norm(Δ) < tol
            return x
        end

        # Newton update: solve (I - J) Δx = Δ
        Δx = (I - J + 1e-6*I) \ Δ
        x  += Δx
    end

    println("Closed orbit did not converge after $maxiter iterations")
    return x
end

function fast_closed_orbit_4d(ring::Vector; x0=zeros(4), energy::Float64=3.5e9,
    mass::Float64=m_e, tol::Float64=1e-8, maxiter::Int=20)
    x = copy(x0)

    function fmap(θ::Vector{Float64})
        rin = [θ[1] θ[2] θ[3] θ[4] 0.0 0.0]
        b = Beam(rin, energy=energy, mass=mass)
        ringpass!(ring, b, 1)
        return b.r[1, 1:4]
    end

    for iter in 1:maxiter
        J, x_out = numerical_jacobian(fmap, x)    # replace AD with FD
        Δ = x_out .- x

        if norm(Δ) < tol
        return x
    end

    # Newton update: solve (I - J) Δx = Δ
    Δx = (I - J + 1e-6*I) \ Δ
    x  += Δx
    end

    println("Closed orbit did not converge after $maxiter iterations")
    return x
end

import Base: mod
function mod(x::DTPSAD{N, T}, y::Float64) where {N, T}
    return x - floor(x / y) * y
end
function getchrom(RING; dp=0.0, order=0, energy=1e9, mass=m_e)
    M66 = findm66(RING,dp, order, E0=energy, m0=mass)
    M44 = M66[1:4, 1:4]
    cos_mu_x = (M44[1,1]+M44[2,2])/2
    cos_mu_y = (M44[3,3]+M44[4,4])/2

    sin_mu_x = sign(M44[1,2])*sqrt(abs(-M44[1,2]*M44[2,1]-(M44[1,1]-M44[2,2])^2/4))
    sin_mu_y = sign(M44[3,4])*sqrt(abs(-M44[3,4]*M44[4,3]-(M44[3,3]-M44[4,4])^2/4))
    tune = mod.(atan.([sin_mu_x, sin_mu_y], [cos_mu_x, cos_mu_y]) / (2π), 1.0)

    M66_dp = findm66(RING,dp+1e-8, order, E0=energy, m0=mass)
    M44_dp = M66_dp[1:4, 1:4]
    cos_mu_x_dp = (M44_dp[1,1]+M44_dp[2,2])/2
    cos_mu_y_dp = (M44_dp[3,3]+M44_dp[4,4])/2

    sin_mu_x_dp = sign(M44_dp[1,2])*sqrt(abs(-M44_dp[1,2]*M44_dp[2,1]-(M44_dp[1,1]-M44_dp[2,2])^2/4))
    sin_mu_y_dp = sign(M44_dp[3,4])*sqrt(abs(-M44_dp[3,4]*M44_dp[4,3]-(M44_dp[3,3]-M44_dp[4,4])^2/4))
    tune_dp = mod.(atan.([sin_mu_x_dp, sin_mu_y_dp], [cos_mu_x_dp, cos_mu_y_dp]) / (2π), 1.0)

    chrom = (tune_dp .- tune) ./ 1e-8
    return chrom
end

function gettune(RING; dp=0.0, order=0, energy=1e9, mass=m_e)
    M66 = findm66(RING,dp, order, E0=energy, m0=mass)
    M44 = M66[1:4, 1:4]
    cos_mu_x = (M44[1,1]+M44[2,2])/2
    cos_mu_y = (M44[3,3]+M44[4,4])/2

    sin_mu_x = sign(M44[1,2])*sqrt(abs(-M44[1,2]*M44[2,1]-(M44[1,1]-M44[2,2])^2/4))
    sin_mu_y = sign(M44[3,4])*sqrt(abs(-M44[3,4]*M44[4,3]-(M44[3,3]-M44[4,4])^2/4))
    tune = mod.(atan.([sin_mu_x, sin_mu_y], [cos_mu_x, cos_mu_y]) / (2π), 1.0)

    return tune
end