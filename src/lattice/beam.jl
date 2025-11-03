"""
    Beam

A struct that contains the information of a beam. The struct contains the following fields:

- `r::Matrix{Float64}`: Nx6 matrix of the 6D phase space coordinates of the beam particles.
- `np::Int`: Number of particles.
- `nmacro::Int`: Number of macro particles.
- `energy::Float64`: Energy of the beam in eV.
- `lost_flag::Vector{Int}`: Vector of length `nmacro` that contains the lost flag of each particle.
- `charge::Float64`: Charge of the beam particles. -1.0 for electrons.
- `mass::Float64`: Mass of the beam in eV.
- `gamma::Float64`: Lorentz factor of the beam particles.
- `beta::Float64`: Velocity of the beam particles in units of the speed of light.
- `atomnum::Float64`: Atomic number of the beam particles.
- `classrad0::Float64`: Classical radiation constant. (not used)
- `radconst::Float64`: Radiation constant. (not used)
- `T0::Float64`: Revolution period of the beam particles. (not used)
- `nturn::Int`: Number of turns in the simulation. (not used)
- `znbin::Int`: Number of bins in the z direction.
- `inzindex::Vector{Int}`: Index of the z bin.
- `zhist::Vector{Float64}`: Histogram of the z direction.
- `zhist_edges::Vector{Float64}`: Edges of the z histogram.
- `temp1::Vector{Float64}`: Temporary variable.
- `temp2::Vector{Float64}`: Temporary variable.
- `temp3::Vector{Float64}`: Temporary variable.
- `temp4::Vector{Float64}`: Temporary variable.
- `temp5::Vector{Float64}`: Temporary variable.
- `emittance::Vector{Float64}`: Emittance in x, y, z directions.
- `centroid::Vector{Float64}`: Centroid in x, px, y, py, z, pz directions.
- `moment2nd::Matrix{Float64}`: 2nd momentum matrix.
- `beamsize::Vector{Float64}`: Beam size in x, px, y, py, z, pz directions.
- `current::Float64`: Beam current for space charge calculation in A.
"""
mutable struct Beam{T}
    r::Matrix{T} # x, px, y, py, z, dp
    np::Int # default is the same as nmacro 
    nmacro::Int 
    energy::T
    lost_flag::Vector{Int}
    charge::T 
    mass::T 
    gamma::T
    beta::T
    atomnum::T
    classrad0::T
    radconst::T
    T0::T # revolution period
    nturn::Int # number of turns
    znbin::Int # number of bins in z
    inzindex::Vector{Int} # index of z bin
    zhist::Vector{T}
    zhist_edges::Vector{T}
    temp1::Vector{T}
    temp2::Vector{T}
    temp3::Vector{T}
    temp4::Vector{T}
    temp5::Vector{T}
    emittance::Vector{T} # emittance in x, y, z
    centroid::Vector{T} # centroid in x, px, y, py, z, pz
    moment2nd::Matrix{T} # 2nd momentum matrix
    beamsize::Vector{T} # Beam size in x, px, y, py, z, pz
    current::T # beam current for space charge calculation
end
"""
    Beam(r::Matrix{Float64}, energy::Float64; np::Int=size(r, 1), charge::Float64=-1.0, mass::Float64=m_e, atn::Float64=1.0,
        emittance::Vector{Float64}=zeros(Float64, 3), centroid::Vector{Float64}=zeros(Float64, 6), T0::Float64=0.0, znbin::Int=99, 
        current::Float64=0.0)

Construct a `Beam` object with coordinates of particles and beam energy. Other parameters are optional.
"""
function Beam(r::Matrix{Float64}, energy::Float64; np::Int=size(r, 1), charge::Float64=-1.0, mass::Float64=m_e, atn::Float64=1.0,
    emittance::Vector{Float64}=zeros(Float64, 3), centroid::Vector{Float64}=zeros(Float64, 6), T0::Float64=0.0, znbin::Int=99, 
    current::Float64=0.0)  

    nmacro = size(r, 1)
    lost_flag = zeros(Int, nmacro)
    gamma = (energy + mass) / mass
    beta = sqrt(1.0 - 1.0 / gamma^2)
    classrad0=charge*charge/(atn*mass)/4/pi/55.26349406*1e-6
    radconst=4*pi/3*classrad0/mass/mass/mass
    inzindex = zeros(Int, nmacro)
    zhist = zeros(Float64, znbin)
    zhist_edges = zeros(Float64, znbin+1)
    return Beam{Float64}(r, np, nmacro, energy, lost_flag, charge, mass, gamma, beta, atn, classrad0, radconst, 
    T0, 1, znbin, inzindex, zhist, zhist_edges, zeros(Float64, nmacro), zeros(Float64, nmacro), zeros(Float64, nmacro), 
    zeros(Float64, nmacro), zeros(Float64, nmacro), emittance, centroid, zeros(Float64, 6, 6), zeros(Float64, 6), current)
end
"""
Beam(r::Matrix{Float64}; energy::Float64=1e9,np::Int=size(r, 1), charge::Float64=-1.0, mass::Float64=m_e, atn::Float64=1.0,
        emittance::Vector{Float64}=zeros(Float64, 3), centroid::Vector{Float64}=zeros(Float64, 6), T0::Float64=0.0, znbin::Int=99, 
        current::Float64=0.0)

Construct a `Beam` object with particle coordinates. Other parameters are optional.
"""
function Beam(r::Matrix{Float64}; energy::Float64=1e9,np::Int=size(r, 1), charge::Float64=-1.0, mass::Float64=m_e, atn::Float64=1.0,
    emittance::Vector{Float64}=zeros(Float64, 3), centroid::Vector{Float64}=zeros(Float64, 6), T0::Float64=0.0, znbin::Int=99, 
    current::Float64=0.0)

    nmacro = size(r, 1)
    lost_flag = zeros(Int, nmacro)
    gamma = (energy + mass) / mass
    beta = sqrt(1.0 - 1.0 / gamma^2)
    classrad0=charge*charge/(atn*mass)/4/pi/55.26349406*1e-6
    radconst=4*pi/3*classrad0/mass/mass/mass
    inzindex = zeros(Int, nmacro)
    zhist = zeros(Float64, znbin)
    zhist_edges = zeros(Float64, znbin+1)
    return Beam{Float64}(r, np, nmacro, energy, lost_flag, charge, mass, gamma, beta, atn, classrad0, radconst, T0, 1, znbin, 
    inzindex, zhist, zhist_edges, zeros(Float64, nmacro), zeros(Float64, nmacro), zeros(Float64, nmacro), zeros(Float64, nmacro), 
    zeros(Float64, nmacro), emittance, centroid, zeros(Float64, 6, 6), zeros(Float64, 6), current)
end
"""
    Beam(;r::Matrix{Float64}=zeros(Float64, 1,6), energy::Float64=1e9, np::Int=size(r, 1), charge::Float64=-1.0, mass::Float64=m_e, 
    atn::Float64=1.0, emittance::Vector{Float64}=zeros(Float64, 3), centroid::Vector{Float64}=zeros(Float64, 6), T0::Float64=0.0, 
    znbin::Int=99, current::Float64=0.0)

Construct a `Beam` object with default parameters. All parameters are optional.
"""
function Beam(;r::Matrix{Float64}=zeros(Float64, 1,6), energy::Float64=1e9, np::Int=size(r, 1), charge::Float64=-1.0, mass::Float64=m_e,
    atn::Float64=1.0, emittance::Vector{Float64}=zeros(Float64, 3), centroid::Vector{Float64}=zeros(Float64, 6), T0::Float64=0.0, 
    znbin::Int=99, current::Float64=0.0)

    nmacro = size(r, 1)
    lost_flag = zeros(Int, nmacro)
    gamma = (energy + mass) / mass
    beta = sqrt(1.0 - 1.0 / gamma^2)
    classrad0=charge*charge/(atn*mass)/4/pi/55.26349406*1e-6
    radconst=4*pi/3*classrad0/mass/mass/mass
    inzindex = zeros(Int, nmacro)
    zhist = zeros(Float64, znbin)
    zhist_edges = zeros(Float64, znbin+1)
    return Beam{Float64}(r, np, nmacro, energy, lost_flag, charge, mass, gamma, beta, atn, classrad0, radconst, T0, 1, znbin, inzindex, zhist, zhist_edges, zeros(Float64, nmacro), zeros(Float64, nmacro), zeros(Float64, nmacro), zeros(Float64, nmacro), zeros(Float64, nmacro), emittance, centroid, zeros(Float64, 6, 6), zeros(Float64, 6), current)
end

function Beam(energy::Float64, np::Int, nmacro::Int; charge::Float64=-1.0, mass::Float64=m_e, atn::Float64=1.0,
    emittance::Vector{Float64}=zeros(Float64, 3), centroid::Vector{Float64}=zeros(Float64, 6), T0::Float64=0.0, 
    znbin::Int=99, current::Float64=0.0)
    r = zeros(Float64, nmacro, 6)
    lost_flag = zeros(Int, nmacro)
    gamma = (energy + mass) / mass
    beta = sqrt(1.0 - 1.0 / gamma^2)
    classrad0=charge*charge/(atn*mass)/4/pi/55.26349406*1e-6
    radconst=4*pi/3*classrad0/mass/mass/mass
    inzindex = zeros(Int, nmacro)
    zhist = zeros(Float64, znbin)
    zhist_edges = zeros(Float64, znbin+1)
    return Beam{Float64}(r, np, nmacro, energy, lost_flag, charge, mass, gamma, beta, atn, classrad0, radconst, T0, 1, znbin, inzindex, zhist, zhist_edges, zeros(Float64, nmacro), zeros(Float64, nmacro), zeros(Float64, nmacro), zeros(Float64, nmacro), zeros(Float64, nmacro), emittance, centroid, zeros(Float64, 6, 6), zeros(Float64, 6), current)
end

"""
    Beam(beam::Beam)

Copy a `Beam` object.
"""
function Beam(beam::Beam{T}) where T
    r = copy(beam.r)
    lost_flag = copy(beam.lost_flag)
    inzindex = copy(beam.inzindex)
    zhist = copy(beam.zhist)
    zhist_edges = copy(beam.zhist_edges)
    emittance = copy(beam.emittance)
    centroid = copy(beam.centroid)
    moment2nd = copy(beam.moment2nd)
    beamsize = copy(beam.beamsize)
    return Beam{T}(r, beam.np, beam.nmacro, beam.energy, lost_flag, beam.charge, beam.mass, beam.gamma, beam.beta, beam.atomnum, beam.classrad0, beam.radconst, beam.T0, beam.nturn, beam.znbin, inzindex, zhist, zhist_edges, beam.temp1, beam.temp2, beam.temp3, beam.temp4, beam.temp5, emittance, centroid, moment2nd, beamsize, beam.current)
end
"""
    Beam(r_GTPS::Matrix{DTPSAD{N, T}}; energy::Float64=1e9, np::Int=size(r_GTPS, 1), charge::Float64=-1.0, mass::Float64=m_e, 
    atn::Float64=1.0, emittance::Vector{Float64}=zeros(Float64, 3), centroid::Vector{Float64}=zeros(Float64, 6), T0::Float64=0.0, 
    znbin::Int=99, current::Float64=0.0) where {N, T <: Number}
Construct a `Beam` object with coordinates of particles in TPSA (DTPSAD type) format. Other parameters are optional.
"""
function Beam(r_GTPS::Matrix{DTPSAD{N, T}}; energy::Float64=1e9, np::Int=size(r_GTPS, 1), charge::Float64=-1.0, mass::Float64=m_e, 
    atn::Float64=1.0, emittance::Vector{Float64}=zeros(Float64, 3), centroid::Vector{Float64}=zeros(Float64, 6), T0::Float64=0.0, 
    znbin::Int=99, current::Float64=0.0) where {N, T <: Number}

    nmacro = size(r_GTPS, 1)
    lost_flag = zeros(Int, nmacro)
    gamma = (DTPSAD(energy) + DTPSAD(mass)) / DTPSAD(mass)
    beta = sqrt(1.0 - 1.0 / gamma^2)
    classrad0 = DTPSAD(charge) * DTPSAD(charge) / (DTPSAD(atn) * DTPSAD(mass)) / 4.0 / Float64(pi) / 55.26349406 * 1e-6
    radconst = 4 * Float64(pi) / 3 * classrad0 / DTPSAD(mass)^3
    inzindex = zeros(Int, nmacro)
    zhist = zeros(DTPSAD{N, T}, znbin)
    zhist_edges = zeros(DTPSAD{N, T}, znbin + 1)
    return Beam{DTPSAD{N, T}}(r_GTPS, np, nmacro, DTPSAD(energy), lost_flag, DTPSAD(charge), DTPSAD(mass), gamma, beta, DTPSAD(atn), classrad0, radconst,
                T0, 1, znbin, inzindex, zhist, zhist_edges,
                zeros(DTPSAD{N, T}, nmacro), zeros(DTPSAD{N, T}, nmacro), zeros(DTPSAD{N, T}, nmacro),
                zeros(DTPSAD{N, T}, nmacro), zeros(DTPSAD{N, T}, nmacro), emittance,
                centroid, zeros(DTPSAD{N, T}, 6, 6), zeros(DTPSAD{N, T}, 6), DTPSAD(current))
end


