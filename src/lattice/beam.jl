mutable struct Beam
    r::Matrix{Float64} # x, px, y, py, z, dp
    np::Int # default is the same as nmacro 
    nmacro::Int 
    energy::Float64
    lost_flag::Vector{Int}
    charge::Float64 
    mass::Float64 
    gamma::Float64
    beta::Float64
    atomnum::Float64
    classrad0::Float64
    radconst::Float64
    T0::Float64 # revolution period
    nturn::Int # number of turns
    znbin::Int # number of bins in z
    inzindex::Vector{Int} # index of z bin
    zhist::Vector{Float64}
    zhist_edges::Vector{Float64}
    temp1::Vector{Float64}
    temp2::Vector{Float64}
    temp3::Vector{Float64}
    temp4::Vector{Float64}
    temp5::Vector{Float64}
    emittance::Vector{Float64} # emittance in x, y, z
    centroid::Vector{Float64} # centroid in x, px, y, py, z, pz
    moment2nd::Matrix{Float64} # 2nd momentum matrix
    beamsize::Vector{Float64} # Beam size in x, px, y, py, z, pz
end
function Beam(r::Matrix{Float64}, energy::Float64; np::Int=size(r, 1), charge::Float64=-1.0, mass::Float64=m_e, atn::Float64=1.0,
        emittance::Vector{Float64}=zeros(Float64, 3), centroid::Vector{Float64}=zeros(Float64, 6), T0::Float64=0.0, znbin::Int=100)  
    nmacro = size(r, 1)
    lost_flag = zeros(Int, nmacro)
    gamma = energy / mass
    beta = sqrt(1.0 - 1.0 / gamma^2)
    classrad0=charge*charge/(atn*mass)/4/pi/55.26349406*1e-6
    radconst=4*pi/3*classrad0/mass/mass/mass
    inzindex = zeros(Int, nmacro)
    zhist = zeros(Float64, znbin)
    zhist_edges = zeros(Float64, znbin+1)
    return Beam(r, np, nmacro, energy, lost_flag, charge, mass, gamma, beta, atn, classrad0, radconst, T0, 1, znbin, inzindex, zhist, zhist_edges, zeros(Float64, nmacro), zeros(Float64, nmacro), zeros(Float64, nmacro), zeros(Float64, nmacro), zeros(Float64, nmacro), emittance, centroid, zeros(Float64, 6, 6), zeros(Float64, 6))
end
function Beam(r::Matrix{Float64}; energy::Float64=1e9,np::Int=size(r, 1), charge::Float64=-1.0, mass::Float64=m_e, atn::Float64=1.0,
        emittance::Vector{Float64}=zeros(Float64, 3), centroid::Vector{Float64}=zeros(Float64, 6), T0::Float64=0.0, znbin::Int=100)
    nmacro = size(r, 1)
    lost_flag = zeros(Int, nmacro)
    gamma = energy / mass
    beta = sqrt(1.0 - 1.0 / gamma^2)
    classrad0=charge*charge/(atn*mass)/4/pi/55.26349406*1e-6
    radconst=4*pi/3*classrad0/mass/mass/mass
    inzindex = zeros(Int, nmacro)
    zhist = zeros(Float64, znbin)
    zhist_edges = zeros(Float64, znbin+1)
    return Beam(r, np, nmacro, energy, lost_flag, charge, mass, gamma, beta, atn, classrad0, radconst, T0, 1, znbin, inzindex, zhist, zhist_edges, zeros(Float64, nmacro), zeros(Float64, nmacro), zeros(Float64, nmacro), zeros(Float64, nmacro), zeros(Float64, nmacro), emittance, centroid, zeros(Float64, 6, 6), zeros(Float64, 6))
end
function Beam(;r::Matrix{Float64}=zeros(Float64, 1,6), energy::Float64=1e9, np::Int=size(r, 1), charge::Float64=-1.0, mass::Float64=m_e, atn::Float64=1.0,
        emittance::Vector{Float64}=zeros(Float64, 3), centroid::Vector{Float64}=zeros(Float64, 6), T0::Float64=0.0, znbin::Int=100)
    nmacro = size(r, 1)
    lost_flag = zeros(Int, nmacro)
    gamma = energy / mass
    beta = sqrt(1.0 - 1.0 / gamma^2)
    classrad0=charge*charge/(atn*mass)/4/pi/55.26349406*1e-6
    radconst=4*pi/3*classrad0/mass/mass/mass
    inzindex = zeros(Int, nmacro)
    zhist = zeros(Float64, znbin)
    zhist_edges = zeros(Float64, znbin+1)
    return Beam(r, np, nmacro, energy, lost_flag, charge, mass, gamma, beta, atn, classrad0, radconst, T0, 1, znbin, inzindex, zhist, zhist_edges, zeros(Float64, nmacro), zeros(Float64, nmacro), zeros(Float64, nmacro), zeros(Float64, nmacro), zeros(Float64, nmacro), emittance, centroid, zeros(Float64, 6, 6), zeros(Float64, 6))
end

function Beam(energy::Float64, np::Int, nmacro::Int; charge::Float64=-1.0, mass::Float64=m_e, atn::Float64=1.0,
    emittance::Vector{Float64}=zeros(Float64, 3), centroid::Vector{Float64}=zeros(Float64, 6), T0::Float64=0.0, znbin::Int=100)
    r = zeros(Float64, nmacro, 6)
    lost_flag = zeros(Int, nmacro)
    gamma = energy / mass
    beta = sqrt(1.0 - 1.0 / gamma^2)
    classrad0=charge*charge/(atn*mass)/4/pi/55.26349406*1e-6
    radconst=4*pi/3*classrad0/mass/mass/mass
    inzindex = zeros(Int, nmacro)
    zhist = zeros(Float64, znbin)
    zhist_edges = zeros(Float64, znbin+1)
    return Beam(r, np, nmacro, energy, lost_flag, charge, mass, gamma, beta, atn, classrad0, radconst, T0, 1, znbin, inzindex, zhist, zhist_edges, zeros(Float64, nmacro), zeros(Float64, nmacro), zeros(Float64, nmacro), zeros(Float64, nmacro), zeros(Float64, nmacro), emittance, centroid, zeros(Float64, 6, 6), zeros(Float64, 6))
end

function Beam(beam::Beam)
    return Beam(beam.r, beam.np, beam.nmacro, beam.energy, beam.lost_flag, beam.charge, beam.mass, beam.gamma, beam.beta, beam.atomnum, beam.classrad0, beam.radconst, beam.T0, beam.nturn, beam.znbin, beam.inzindex, beam.zhist, beam.zhist_edges, beam.temp1, beam.temp2, beam.temp3, beam.temp4, beam.temp5, beam.emittance, beam.centroid, beam.moment2nd, beam.beamsize)
end

