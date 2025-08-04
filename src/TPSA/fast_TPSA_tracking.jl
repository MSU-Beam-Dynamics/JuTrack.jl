# use first-order TPSA. It is faster than other TPSA modules.
using .TPSAadStatic
mutable struct TBeam{N, T <: Number}
    r::Matrix{DTPSAD{N, T}}
    np::Int # default is the same as nmacro
    nmacro::Int
    energy::DTPSAD{N, T} # energy in eV
    lost_flag::Vector{Int}
    charge::DTPSAD{N, T}
    mass::DTPSAD{N, T}
    gamma::DTPSAD{N, T}
    beta::DTPSAD{N, T}
    atomnum::DTPSAD{N, T}
    classrad0::DTPSAD{N, T}
    radconst::DTPSAD{N, T}
    T0::Float64 # revolution period
    nturn::Int # number of turns
    znbin::Int # number of bins in z
    inzindex::Vector{Int} # index of z bin
    zhist::Vector
    zhist_edges::Vector
    temp1::Vector
    temp2::Vector
    temp3::Vector
    temp4::Vector
    temp5::Vector
    emittance::Vector # emittance in x, y, z
    centroid::Vector # centroid in x, px, y, py, z, pz
    moment2nd::Matrix # 2nd momentum matrix
    beamsize::Vector # Beam size in x, px, y, py, z, pz
    current::DTPSAD{N, T} # beam current for space charge calculation
end

function TBeam(r::Matrix{Float64}; energy::Float64=1e9, np::Int=size(r, 1), charge::Float64=-1.0, mass::Float64=m_e, atn::Float64=1.0,
        emittance::Vector{Float64}=zeros(Float64, 3), centroid::Vector{Float64}=zeros(Float64, 6), T0::Float64=0.0, znbin::Int=99, current::Float64=0.0)
    nmacro = size(r, 1)
    r_GTPS = DTPSAD.(r)
    lost_flag = zeros(Int, nmacro)
    gamma = DTPSAD(energy) / DTPSAD(mass)
    beta = sqrt(1.0 - 1.0 / gamma^2)
    classrad0 = DTPSAD(charge) * DTPSAD(charge) / (DTPSAD(atn) * DTPSAD(mass)) / 4.0 / Float64(pi) / 55.26349406 * 1e-6
    radconst = 4 * Float64(pi) / 3 * classrad0 / DTPSAD(mass)^3
    inzindex = zeros(Int, nmacro)
    zhist = zeros(DTPSAD{NVAR(), Float64}, znbin)
    zhist_edges = zeros(DTPSAD{NVAR(), Float64}, znbin + 1)
    return TBeam{NVAR(), Float64}(r_GTPS, np, nmacro, DTPSAD(energy), lost_flag, DTPSAD(charge), DTPSAD(mass), gamma, beta, DTPSAD(atn), classrad0, radconst,
                T0, 1, znbin, inzindex, zhist, zhist_edges,
                zeros(DTPSAD{NVAR(), Float64}, nmacro), zeros(DTPSAD{NVAR(), Float64}, nmacro), zeros(DTPSAD{NVAR(), Float64}, nmacro),
                zeros(DTPSAD{NVAR(), Float64}, nmacro), zeros(DTPSAD{NVAR(), Float64}, nmacro), emittance,
                centroid, zeros(DTPSAD{NVAR(), Float64}, 6, 6), zeros(DTPSAD{NVAR(), Float64}, 6), DTPSAD(current))
end

function TBeam(r_GTPS::Matrix{DTPSAD{N, T}}; energy::Float64=1e9, np::Int=size(r_GTPS, 1), charge::Float64=-1.0, mass::Float64=m_e, atn::Float64=1.0,
        emittance::Vector{Float64}=zeros(Float64, 3), centroid::Vector{Float64}=zeros(Float64, 6), T0::Float64=0.0, znbin::Int=99, current::Float64=0.0) where {N, T <: Number}
    nmacro = size(r_GTPS, 1)
    lost_flag = zeros(Int, nmacro)
    gamma = DTPSAD(energy) / DTPSAD(mass)
    beta = sqrt(1.0 - 1.0 / gamma^2)
    classrad0 = DTPSAD(charge) * DTPSAD(charge) / (DTPSAD(atn) * DTPSAD(mass)) / 4.0 / Float64(pi) / 55.26349406 * 1e-6
    radconst = 4 * Float64(pi) / 3 * classrad0 / DTPSAD(mass)^3
    inzindex = zeros(Int, nmacro)
    zhist = zeros(DTPSAD{N, T}, znbin)
    zhist_edges = zeros(DTPSAD{N, T}, znbin + 1)
    return TBeam{N, T}(r_GTPS, np, nmacro, DTPSAD(energy), lost_flag, DTPSAD(charge), DTPSAD(mass), gamma, beta, DTPSAD(atn), classrad0, radconst,
                T0, 1, znbin, inzindex, zhist, zhist_edges,
                zeros(DTPSAD{N, T}, nmacro), zeros(DTPSAD{N, T}, nmacro), zeros(DTPSAD{N, T}, nmacro),
                zeros(DTPSAD{N, T}, nmacro), zeros(DTPSAD{N, T}, nmacro), emittance,
                centroid, zeros(DTPSAD{N, T}, 6, 6), zeros(DTPSAD{N, T}, 6), DTPSAD(current))
end

abstract type AbstractTPSAElement <: AbstractElement end

mutable struct TMARKER{N, T} <: AbstractTPSAElement
    name::String
    len::DTPSAD{N, T}
    eletype::String
end
function TMARKER(; name::String = "MARKER", len::T = 0.0) where T <: Number
    TMARKER{NVAR(), T}(name, DTPSAD(len), "MARKER")
end
function TMARKER(name::String, len::DTPSAD{N, T}, eletype::String) where {N, T <: Number}
    TMARKER{N, T}(name, len, eletype)
end

mutable struct TDRIFT{N, T} <: AbstractTPSAElement
    name::String
    len::DTPSAD{N, T}
    T1::Array{DTPSAD{N, T},1}
    T2::Array{DTPSAD{N, T},1}
    R1::Array{DTPSAD{N, T},2}
    R2::Array{DTPSAD{N, T},2}
    RApertures::Array{Float64,1}
    EApertures::Array{Float64,1}
    eletype::String
end
function TDRIFT(;name::String = "DRIFT", len::T = 0.0, T1::Array{T,1} = zeros(6), 
                T2::Array{T,1} = zeros(6), R1::Array{T,2} = zeros(6, 6), R2::Array{T,2} = zeros(6, 6), 
                RApertures::Array{Float64,1} = zeros(6), EApertures::Array{Float64,1} = zeros(6)) where T <: Number
    TDRIFT{NVAR(), T}(name, DTPSAD(len), DTPSAD.(T1), DTPSAD.(T2), DTPSAD.(R1), DTPSAD.(R2), RApertures, EApertures, "DRIFT")
end
TDRIFT(name::String, len::DTPSAD{N, T}, T1::Array{DTPSAD{N, T},1}, 
        T2::Array{DTPSAD{N, T},1}, R1::Array{DTPSAD{N, T},2}, R2::Array{DTPSAD{N, T},2}, 
        RApertures::Array{Float64,1}, EApertures::Array{Float64,1}, eletype::String) where {N, T <: Number} =
    TDRIFT{N, T}(name, len, T1, T2, R1, R2, RApertures, EApertures, eletype)

mutable struct TQUAD{N, T} <: AbstractTPSAElement
    name::String
    len::DTPSAD{N, T}
    k1::DTPSAD{N, T}
    rad::Int64
    T1::Array{DTPSAD{N, T},1}
    T2::Array{DTPSAD{N, T},1}
    R1::Array{DTPSAD{N, T},2}
    R2::Array{DTPSAD{N, T},2}
    RApertures::Array{Float64,1}
    EApertures::Array{Float64,1}
    eletype::String
end
function TQUAD(;name::String = "Quad", len::T = 0.0, k1::T = 0.0, rad::Int64 = 0, 
                T1::Array{T,1} = zeros(6), T2::Array{T,1} = zeros(6), 
                R1::Array{T,2} = zeros(6, 6), R2::Array{T,2} = zeros(6, 6), 
                RApertures::Array{Float64,1} = zeros(6), EApertures::Array{Float64,1} = zeros(6)) where T <: Number
    TQUAD{NVAR(), T}(name, DTPSAD(len), DTPSAD(k1), rad, DTPSAD.(T1), DTPSAD.(T2), DTPSAD.(R1), DTPSAD.(R2), RApertures, EApertures, "QUAD")
end
TQUAD(name::String, len::DTPSAD{N, T}, k1::DTPSAD{N, T}, rad::Int64, 
        T1::Array{DTPSAD{N, T},1}, T2::Array{DTPSAD{N, T},1}, 
        R1::Array{DTPSAD{N, T},2}, R2::Array{DTPSAD{N, T},2}, 
        RApertures::Array{Float64,1}, EApertures::Array{Float64,1}, eletype::String) where {N, T <: Number} =
    TQUAD{N, T}(name, len, k1, rad, T1, T2, R1, R2, RApertures, EApertures, eletype)

mutable struct TKQUAD{N, T} <: AbstractTPSAElement
    name::String
    len::DTPSAD{N, T}
    k1::DTPSAD{N, T}
    PolynomA::Array{DTPSAD{N, T},1}
    PolynomB::Array{DTPSAD{N, T},1}
    MaxOrder::Int64
    NumIntSteps::Int64
    rad::Int64
    FringeQuadEntrance::Int64
    FringeQuadExit::Int64
    T1::Array{DTPSAD{N, T},1}
    T2::Array{DTPSAD{N, T},1}
    R1::Array{DTPSAD{N, T},2}
    R2::Array{DTPSAD{N, T},2}
    RApertures::Array{Float64,1}
    EApertures::Array{Float64,1}
    KickAngle::Array{DTPSAD{N, T},1}
    eletype::String
end
function TKQUAD(;name::String = "Quad", len::T = 0.0, k1::T = 0.0, 
            PolynomA::Array{T,1} = zeros(4), 
            PolynomB::Array{T,1} = zeros(4), MaxOrder::Int64=1, 
            NumIntSteps::Int64 = 10, rad::Int64=0, FringeQuadEntrance::Int64 = 0, 
            FringeQuadExit::Int64 = 0, T1::Array{T,1} = zeros(6), 
            T2::Array{T,1} = zeros(6), R1::Array{T,2} = zeros(6, 6), 
            R2::Array{T,2} = zeros(6, 6), RApertures::Array{Float64,1} = zeros(Float64, 6), 
            EApertures::Array{Float64,1} = zeros(Float64, 6), KickAngle::Array{T,1} = zeros(2)) where T <: Number
    if k1 != 0.0 && PolynomB[2] == 0.0
        PolynomB[2] = k1
    end
    TKQUAD{NVAR(), T}(name, DTPSAD(len), DTPSAD(k1), DTPSAD.(PolynomA), DTPSAD.(PolynomB), MaxOrder, NumIntSteps, rad, FringeQuadEntrance, FringeQuadExit,
        DTPSAD.(T1), DTPSAD.(T2), DTPSAD.(R1), DTPSAD.(R2), RApertures, EApertures, DTPSAD.(KickAngle), "KQUAD")
end
TKQUAD(name::String, len::DTPSAD{N, T}, k1::DTPSAD{N, T}, 
        PolynomA::Array{DTPSAD{N, T},1}, PolynomB::Array{DTPSAD{N, T},1}, MaxOrder::Int64, 
        NumIntSteps::Int64, rad::Int64, FringeQuadEntrance::Int64, FringeQuadExit::Int64, 
        T1::Array{DTPSAD{N, T},1}, T2::Array{DTPSAD{N, T},1}, 
        R1::Array{DTPSAD{N, T},2}, R2::Array{DTPSAD{N, T},2}, 
        RApertures::Array{Float64,1}, EApertures::Array{Float64,1}, KickAngle::Array{DTPSAD{N, T},1}, eletype::String) where {N, T <: Number} =
    TKQUAD{N, T}(name, len, k1, PolynomA, PolynomB, MaxOrder, NumIntSteps, rad, FringeQuadEntrance, FringeQuadExit,
        T1, T2, R1, R2, RApertures, EApertures, KickAngle, eletype)

mutable struct TKSEXT{N, T} <: AbstractTPSAElement
    name::String
    len::DTPSAD{N, T}
    k2::DTPSAD{N, T}
    PolynomA::Array{DTPSAD{N, T},1}
    PolynomB::Array{DTPSAD{N, T},1}
    MaxOrder::Int64
    NumIntSteps::Int64
    rad::Int64
    FringeQuadEntrance::Int64
    FringeQuadExit::Int64
    T1::Array{DTPSAD{N, T},1}
    T2::Array{DTPSAD{N, T},1}
    R1::Array{DTPSAD{N, T},2}
    R2::Array{DTPSAD{N, T},2}
    RApertures::Array{Float64,1}
    EApertures::Array{Float64,1}
    KickAngle::Array{DTPSAD{N, T},1}
    eletype::String
end
function TKSEXT(;name::String = "Sext", len::T = 0.0, k2::T = 0.0, 
                PolynomA::Array{T,1} = zeros(4), 
                PolynomB::Array{T,1} = zeros(4), MaxOrder::Int64=2, 
                NumIntSteps::Int64 = 10, rad::Int64=0, FringeQuadEntrance::Int64 = 0, 
                FringeQuadExit::Int64 = 0, T1::Array{T,1} = zeros(6), 
                T2::Array{T,1} = zeros(6), R1::Array{T,2} = zeros(6, 6), 
                R2::Array{T,2} = zeros(6, 6), RApertures::Array{Float64,1} = zeros(Float64, 6), 
                EApertures::Array{Float64,1} = zeros(Float64, 6), KickAngle::Array{T,1} = zeros(2)) where T <: Number
    if k2 != 0.0 && PolynomB[3] == 0.0
        PolynomB[3] = k2
    end
    TKSEXT{NVAR(), T}(name, DTPSAD(len), DTPSAD(k2), DTPSAD.(PolynomA), DTPSAD.(PolynomB), MaxOrder, NumIntSteps, rad, FringeQuadEntrance, FringeQuadExit,
        DTPSAD.(T1), DTPSAD.(T2), DTPSAD.(R1), DTPSAD.(R2), RApertures, EApertures, DTPSAD.(KickAngle), "KSEXT")
end
TKSEXT(name::String, len::DTPSAD{N, T}, k2::DTPSAD{N, T}, 
        PolynomA::Array{DTPSAD{N, T},1}, PolynomB::Array{DTPSAD{N, T},1}, MaxOrder::Int64, 
        NumIntSteps::Int64, rad::Int64, FringeQuadEntrance::Int64, FringeQuadExit::Int64, 
        T1::Array{DTPSAD{N, T},1}, T2::Array{DTPSAD{N, T},1}, 
        R1::Array{DTPSAD{N, T},2}, R2::Array{DTPSAD{N, T},2}, 
        RApertures::Array{Float64,1}, EApertures::Array{Float64,1}, KickAngle::Array{DTPSAD{N, T},1}, eletype::String) where {N, T <: Number} =
    TKSEXT{N, T}(name, len, k2, PolynomA, PolynomB, MaxOrder, NumIntSteps, rad, FringeQuadEntrance, FringeQuadExit,
        T1, T2, R1, R2, RApertures, EApertures, KickAngle, eletype)

mutable struct TKOCT{N, T} <: AbstractTPSAElement
    name::String
    len::DTPSAD{N, T}
    k3::DTPSAD{N, T}
    PolynomA::Array{DTPSAD{N, T},1}
    PolynomB::Array{DTPSAD{N, T},1}
    MaxOrder::Int64
    NumIntSteps::Int64
    rad::Int64
    FringeQuadEntrance::Int64
    FringeQuadExit::Int64
    T1::Array{DTPSAD{N, T},1}
    T2::Array{DTPSAD{N, T},1}
    R1::Array{DTPSAD{N, T},2}
    R2::Array{DTPSAD{N, T},2}
    RApertures::Array{Float64,1}
    EApertures::Array{Float64,1}
    KickAngle::Array{DTPSAD{N, T},1}
    eletype::String
end
function TKOCT(;name::String = "OCT", len::T = 0.0, k3::T = 0.0, 
                PolynomA::Array{T,1} = zeros(4), 
                PolynomB::Array{T,1} = zeros(4), MaxOrder::Int64=3, 
                NumIntSteps::Int64 = 10, rad::Int64=0, FringeQuadEntrance::Int64 = 0, 
                FringeQuadExit::Int64 = 0, T1::Array{T,1} = zeros(6), 
                T2::Array{T,1} = zeros(6), R1::Array{T,2} = zeros(6, 6), 
                R2::Array{T,2} = zeros(6, 6), RApertures::Array{Float64,1} = zeros(Float64, 6), 
                EApertures::Array{Float64,1} = zeros(Float64, 6), KickAngle::Array{T,1} = zeros(2)) where T <: Number
    if k3 != 0.0 && PolynomB[4] == 0.0
        PolynomB[4] = k3
    end
    TKOCT{NVAR(), T}(name, DTPSAD(len), DTPSAD(k3), DTPSAD.(PolynomA), DTPSAD.(PolynomB), MaxOrder, NumIntSteps, rad, FringeQuadEntrance, FringeQuadExit,
        DTPSAD.(T1), DTPSAD.(T2), DTPSAD.(R1), DTPSAD.(R2), RApertures, EApertures, DTPSAD.(KickAngle), "KOCT")
end
TKOCT(name::String, len::DTPSAD{N, T}, k3::DTPSAD{N, T}, 
        PolynomA::Array{DTPSAD{N, T},1}, PolynomB::Array{DTPSAD{N, T},1}, MaxOrder::Int64, 
        NumIntSteps::Int64, rad::Int64, FringeQuadEntrance::Int64, FringeQuadExit::Int64, 
        T1::Array{DTPSAD{N, T},1}, T2::Array{DTPSAD{N, T},1}, 
        R1::Array{DTPSAD{N, T},2}, R2::Array{DTPSAD{N, T},2}, 
        RApertures::Array{Float64,1}, EApertures::Array{Float64,1}, KickAngle::Array{DTPSAD{N, T},1}, eletype::String) where {N, T <: Number} =
    TKOCT{N, T}(name, len, k3, PolynomA, PolynomB, MaxOrder, NumIntSteps, rad, FringeQuadEntrance, FringeQuadExit,
        T1, T2, R1, R2, RApertures, EApertures, KickAngle, eletype)

mutable struct TCORRECTOR{N, T} <: AbstractTPSAElement
    name::String
    len::DTPSAD{N, T}
    xkick::DTPSAD{N, T}
    ykick::DTPSAD{N, T}
    T1::Array{DTPSAD{N, T},1}
    T2::Array{DTPSAD{N, T},1}
    R1::Array{DTPSAD{N, T},2}
    R2::Array{DTPSAD{N, T},2}
    eletype::String
end
function TCORRECTOR(;name::String = "CORRECTOR", len::T = 0.0, xkick::T = 0.0, ykick::T = 0.0, 
                    T1::Array{T,1} = zeros(6), T2::Array{T,1} = zeros(6), R1::Array{T,2} = zeros(6, 6), 
                    R2::Array{T,2} = zeros(6, 6)) where T <: Number
    TCORRECTOR{NVAR(), T}(name, DTPSAD(len), DTPSAD(xkick), DTPSAD(ykick), DTPSAD.(T1), DTPSAD.(T2), DTPSAD.(R1), DTPSAD.(R2), "CORRECTOR")
end
TCORRECTOR(name::String, len::DTPSAD{N, T}, xkick::DTPSAD{N, T}, ykick::DTPSAD{N, T}, 
            T1::Array{DTPSAD{N, T},1}, T2::Array{DTPSAD{N, T},1}, 
            R1::Array{DTPSAD{N, T},2}, R2::Array{DTPSAD{N, T},2}, eletype::String) where {N, T <: Number} =
    TCORRECTOR{N, T}(name, len, xkick, ykick, T1, T2, R1, R2, eletype)

mutable struct TSOLENOID{N, T} <: AbstractTPSAElement
    name::String
    len::DTPSAD{N, T}
    ks::DTPSAD{N, T} # rad/m
    T1::Array{DTPSAD{N, T},1}
    T2::Array{DTPSAD{N, T},1}
    R1::Array{DTPSAD{N, T},2}
    R2::Array{DTPSAD{N, T},2}
    eletype::String
end
function TSOLENOID(;name::String = "Solenoid", len::T = 0.0, ks::T = 0.0, T1::Array{T,1} = zeros(6), 
                T2::Array{T,1} = zeros(6), R1::Array{T,2} = zeros(6, 6), R2::Array{T,2} = zeros(6, 6)) where T <: Number
    TSOLENOID{NVAR(), T}(name, DTPSAD(len), DTPSAD(ks), DTPSAD.(T1), DTPSAD.(T2), DTPSAD.(R1), DTPSAD.(R2), "SOLENOID")
end
TSOLENOID(name::String, len::DTPSAD{N, T}, ks::DTPSAD{N, T}, 
            T1::Array{DTPSAD{N, T},1}, T2::Array{DTPSAD{N, T},1}, 
            R1::Array{DTPSAD{N, T},2}, R2::Array{DTPSAD{N, T},2}, eletype::String) where {N, T <: Number} =
    TSOLENOID{N, T}(name, len, ks, T1, T2, R1, R2, eletype)

mutable struct TSBEND{N, T} <: AbstractTPSAElement
    name::String
    len::DTPSAD{N, T}
    angle::DTPSAD{N, T}
    e1::DTPSAD{N, T}
    e2::DTPSAD{N, T}
    PolynomA::Array{DTPSAD{N, T},1}
    PolynomB::Array{DTPSAD{N, T},1}
    MaxOrder::Int64
    NumIntSteps::Int64
    rad::Int64
    fint1::DTPSAD{N, T}
    fint2::DTPSAD{N, T}
    gap::DTPSAD{N, T}
    FringeBendEntrance::Int64
    FringeBendExit::Int64
    FringeQuadEntrance::Int64
    FringeQuadExit::Int64
    FringeIntM0::Array{DTPSAD{N, T},1}
    FringeIntP0::Array{DTPSAD{N, T},1}
    T1::Array{DTPSAD{N, T},1}
    T2::Array{DTPSAD{N, T},1}
    R1::Array{DTPSAD{N, T},2}
    R2::Array{DTPSAD{N, T},2}
    RApertures::Array{Float64,1}
    EApertures::Array{Float64,1}
    KickAngle::Array{DTPSAD{N, T},1}
    eletype::String
end
function TSBEND(;name::String = "SBend", len::T = 0.0, angle::T = 0.0, e1::T = 0.0, e2::T = 0.0, 
                PolynomA::Array{T,1} = zeros(4), PolynomB::Array{T,1} = zeros(4), 
                MaxOrder::Int64=0, NumIntSteps::Int64 = 10, rad::Int64=0, fint1::T = 0.0, fint2::T = 0.0, 
                gap::T = 0.0, FringeBendEntrance::Int64 = 1, FringeBendExit::Int64 = 1, 
                FringeQuadEntrance::Int64 = 0, FringeQuadExit::Int64 = 0, FringeIntM0::Array{T,1} = zeros(5), 
                FringeIntP0::Array{T,1} = zeros(5), T1::Array{T,1} = zeros(6), 
                T2::Array{T,1} = zeros(6), R1::Array{T,2} = zeros(6, 6), 
                R2::Array{T,2} = zeros(6, 6), RApertures::Array{Float64,1} = zeros(Float64, 6), 
                EApertures::Array{Float64,1} = zeros(Float64, 6), KickAngle::Array{T,1} = zeros(2)) where T <: Number
    if PolynomB[2] != 0.0
        MaxOrder = 1
    end
    if PolynomB[3] != 0.0
        MaxOrder = 2
    end
    if PolynomB[4] != 0.0
        MaxOrder = 3
    end
    TSBEND{NVAR(), T}(name, DTPSAD(len), DTPSAD(angle), DTPSAD(e1), DTPSAD(e2), DTPSAD.(PolynomA), DTPSAD.(PolynomB), MaxOrder, NumIntSteps, rad, DTPSAD(fint1), DTPSAD(fint2), DTPSAD(gap),
        FringeBendEntrance, FringeBendExit, FringeQuadEntrance, FringeQuadExit, DTPSAD.(FringeIntM0), DTPSAD.(FringeIntP0),
        DTPSAD.(T1), DTPSAD.(T2), DTPSAD.(R1), DTPSAD.(R2), RApertures, EApertures, DTPSAD.(KickAngle), "SBEND")
end
TSBEND(name::String, len::DTPSAD{N, T}, angle::DTPSAD{N, T}, e1::DTPSAD{N, T}, e2::DTPSAD{N, T}, 
        PolynomA::Array{DTPSAD{N, T},1}, PolynomB::Array{DTPSAD{N, T},1}, MaxOrder::Int64, 
        NumIntSteps::Int64, rad::Int64, fint1::DTPSAD{N, T}, fint2::DTPSAD{N, T}, 
        gap::DTPSAD{N, T}, FringeBendEntrance::Int64, FringeBendExit::Int64, 
        FringeQuadEntrance::Int64, FringeQuadExit::Int64, FringeIntM0::Array{DTPSAD{N, T},1}, 
        FringeIntP0::Array{DTPSAD{N, T},1}, T1::Array{DTPSAD{N, T},1}, 
        T2::Array{DTPSAD{N, T},1}, R1::Array{DTPSAD{N, T},2}, R2::Array{DTPSAD{N, T},2},
        RApertures::Array{Float64,1}, EApertures::Array{Float64,1}, KickAngle::Array{DTPSAD{N, T},1},
        eletype::String) where {N, T <: Number} =
    TSBEND{N, T}(name, len, angle, e1, e2, PolynomA, PolynomB, MaxOrder, NumIntSteps, rad, fint1, fint2, gap,
        FringeBendEntrance, FringeBendExit, FringeQuadEntrance, FringeQuadExit,
        FringeIntM0, FringeIntP0,
        T1, T2, R1, R2, RApertures, EApertures, KickAngle, eletype)

function TRBEND(;name::String = "RBend", len::T = 0.0, angle::T = 0.0, PolynomA::Array{T,1} = zeros(4), 
                PolynomB::Array{T,1} = zeros(4), MaxOrder::Int64=0, NumIntSteps::Int64 = 10, rad::Int64=0, fint1::T = 0.0, 
                fint2::T = 0.0, gap::T = 0.0, FringeBendEntrance::Int64 = 1, FringeBendExit::Int64 = 1, 
                FringeQuadEntrance::Int64 = 0, FringeQuadExit::Int64 = 0, FringeIntM0::Array{T,1} = zeros(5), 
                FringeIntP0::Array{T,1} = zeros(5), T1::Array{T,1} = zeros(6), T2::Array{T,1} = zeros(6), 
                R1::Array{T,2} = zeros(6, 6), R2::Array{T,2} = zeros(6, 6), RApertures::Array{Float64,1} = zeros(Float64, 6), 
                EApertures::Array{Float64,1} = zeros(Float64, 6), KickAngle::Array{T,1} = zeros(2)) where T <: Number
    e1 = angle/2.0
    e2 = angle/2.0
    if PolynomB[2] != 0.0
        MaxOrder = 1
    end
    if PolynomB[3] != 0.0
        MaxOrder = 2
    end
    if PolynomB[4] != 0.0
        MaxOrder = 3
    end
    return TSBEND(name=name, len=len, angle=angle, e1=e1, e2=e2, PolynomA=PolynomA, PolynomB=PolynomB, MaxOrder=MaxOrder, 
                NumIntSteps=NumIntSteps, rad=rad, fint1=fint1, fint2=fint2, gap=gap, FringeBendEntrance=FringeBendEntrance, 
                FringeBendExit=FringeBendExit, FringeQuadEntrance=FringeQuadEntrance, FringeQuadExit=FringeQuadExit, 
                FringeIntM0=FringeIntM0, FringeIntP0=FringeIntP0, T1=T1, T2=T2, R1=R1, R2=R2, RApertures=RApertures, 
                EApertures=EApertures, KickAngle=KickAngle)
end

mutable struct TESBEND{N, T} <: AbstractTPSAElement
    name::String
    len::DTPSAD{N, T}
    angle::DTPSAD{N, T}
    e1::DTPSAD{N, T}
    e2::DTPSAD{N, T}
    PolynomA::Array{DTPSAD{N, T},1}
    PolynomB::Array{DTPSAD{N, T},1}
    MaxOrder::Int64
    NumIntSteps::Int64
    rad::Int64
    gK::DTPSAD{N, T}
    FringeBendEntrance::Int64
    FringeBendExit::Int64
    FringeQuadEntrance::Int64
    FringeQuadExit::Int64
    T1::Array{DTPSAD{N, T},1}
    T2::Array{DTPSAD{N, T},1}
    R1::Array{DTPSAD{N, T},2}
    R2::Array{DTPSAD{N, T},2}
    RApertures::Array{Float64,1}
    EApertures::Array{Float64,1}
    KickAngle::Array{DTPSAD{N, T},1}
    eletype::String
end
function TESBEND(;name::String = "ESBend", len::T = 0.0, angle::T = 0.0, e1::T = 0.0, e2::T = 0.0, 
                PolynomA::Array{T,1} = zeros(4), PolynomB::Array{T,1} = zeros(4), 
                MaxOrder::Int64=0, NumIntSteps::Int64 = 10, rad::Int64=0,
                gK::T = 0.0, FringeBendEntrance::Int64 = 1, FringeBendExit::Int64 = 1, 
                FringeQuadEntrance::Int64 = 0, FringeQuadExit::Int64 = 0, T1::Array{T,1} = zeros(6), 
                T2::Array{T,1} = zeros(6), R1::Array{T,2} = zeros(6, 6), 
                R2::Array{T,2} = zeros(6, 6), RApertures::Array{Float64,1} = zeros(Float64, 6), 
                EApertures::Array{Float64,1} = zeros(Float64, 6), KickAngle::Array{T,1} = zeros(2)) where T <: Number
    if PolynomB[2] != 0.0
        MaxOrder = 1
    end
    if PolynomB[3] != 0.0
        MaxOrder = 2
    end
    if PolynomB[4] != 0.0
        MaxOrder = 3
    end
    TESBEND{NVAR(), T}(name, DTPSAD(len), DTPSAD(angle), DTPSAD(e1), DTPSAD(e2), DTPSAD.(PolynomA), DTPSAD.(PolynomB), MaxOrder, NumIntSteps, rad, DTPSAD(gK), 
        FringeBendEntrance, FringeBendExit, FringeQuadEntrance, FringeQuadExit, 
        DTPSAD.(T1), DTPSAD.(T2), DTPSAD.(R1), DTPSAD.(R2), RApertures, EApertures, DTPSAD.(KickAngle), "ESBEND")
end
TESBEND(name::String, len::DTPSAD{N, T}, angle::DTPSAD{N, T}, e1::DTPSAD{N, T}, e2::DTPSAD{N, T}, 
        PolynomA::Array{DTPSAD{N, T},1}, PolynomB::Array{DTPSAD{N, T},1}, MaxOrder::Int64, 
        NumIntSteps::Int64, rad::Int64, gK::DTPSAD{N, T}, FringeBendEntrance::Int64, FringeBendExit::Int64, 
        FringeQuadEntrance::Int64, FringeQuadExit::Int64, T1::Array{DTPSAD{N, T},1}, 
        T2::Array{DTPSAD{N, T},1}, R1::Array{DTPSAD{N, T},2}, R2::Array{DTPSAD{N, T},2},
        RApertures::Array{Float64,1}, EApertures::Array{Float64,1}, KickAngle::Array{DTPSAD{N, T},1},
        eletype::String) where {N, T <: Number} =
    TESBEND{N, T}(name, len, angle, e1, e2, PolynomA, PolynomB, MaxOrder, NumIntSteps, rad, gK,
        FringeBendEntrance, FringeBendExit, FringeQuadEntrance, FringeQuadExit,
        T1, T2, R1, R2, RApertures, EApertures, KickAngle, eletype)

function TERBEND(;name::String = "ESBend", len::T = 0.0, angle::T = 0.0, 
    PolynomA::Array{T,1} = zeros(4), PolynomB::Array{T,1} = zeros(4), 
    MaxOrder::Int64=0, NumIntSteps::Int64 = 10, rad::Int64=0,
    gK::T = 0.0, FringeBendEntrance::Int64 = 1, FringeBendExit::Int64 = 1, 
    FringeQuadEntrance::Int64 = 0, FringeQuadExit::Int64 = 0, T1::Array{T,1} = zeros(6), 
    T2::Array{T,1} = zeros(6), R1::Array{T,2} = zeros(6, 6), 
    R2::Array{T,2} = zeros(6, 6), RApertures::Array{Float64,1} = zeros(Float64, 6), 
    EApertures::Array{Float64,1} = zeros(Float64, 6), KickAngle::Array{T,1} = zeros(2)) where T <: Number

    e1 = angle/2.0
    e2 = angle/2.0
    if PolynomB[2] != 0.0
        MaxOrder = 1
    end
    if PolynomB[3] != 0.0
        MaxOrder = 2
    end
    if PolynomB[4] != 0.0
        MaxOrder = 3
    end
    return TESBEND(name=name, len=len, angle=angle, e1=e1, e2=e2, PolynomA=PolynomA, PolynomB=PolynomB, MaxOrder=MaxOrder, 
                NumIntSteps=NumIntSteps, rad=rad, gK=gK, FringeBendEntrance=FringeBendEntrance, 
                FringeBendExit=FringeBendExit, FringeQuadEntrance=FringeQuadEntrance, FringeQuadExit=FringeQuadExit, 
                T1=T1, T2=T2, R1=R1, R2=R2, RApertures=RApertures, EApertures=EApertures, KickAngle=KickAngle)
end

mutable struct TRFCA{N, T} <: AbstractTPSAElement
    name::String
    len::DTPSAD{N, T}
    volt::DTPSAD{N, T}
    freq::DTPSAD{N, T}
    h::DTPSAD{N, T}
    lag::DTPSAD{N, T}
    philag::DTPSAD{N, T}
    energy::DTPSAD{N, T}
    eletype::String
end
function TRFCA(;name::String = "RFCA", len::T = 0.0, volt::T = 0.0, freq::T = 0.0, h::T = 1.0, 
                lag::T = 0.0, philag::T = 0.0, energy::T = 0.0) where T <: Number
    TRFCA{NVAR(), T}(name, DTPSAD(len), DTPSAD(volt), DTPSAD(freq), DTPSAD(h), DTPSAD(lag), DTPSAD(philag), DTPSAD(energy), "RFCA")
end
TRFCA(name::String, len::DTPSAD{N, T}, volt::DTPSAD{N, T}, freq::DTPSAD{N, T}, 
        h::DTPSAD{N, T}, lag::DTPSAD{N, T}, philag::DTPSAD{N, T}, energy::DTPSAD{N, T}, eletype::String) where {N, T <: Number} =
    TRFCA{N, T}(name, len, volt, freq, h, lag, philag, energy, eletype)

mutable struct TCRABCAVITY{N, T} <: AbstractTPSAElement
    name::String 
    len::DTPSAD{N, T}
    volt::DTPSAD{N, T}  # voltage
    freq::DTPSAD{N, T}  # frequency
    k::DTPSAD{N, T}  # wave number
    phi::DTPSAD{N, T}  # phase
    errors::Array{DTPSAD{N, T},1} # 1: Voltage error, 2: Phase error
    energy::DTPSAD{N, T}
    eletype::String 
end
function TCRABCAVITY(;name::String = "CRABCAVITY", len::T = 0.0, volt::T = 0.0, 
    freq::T = 0.0, phi::T = 0.0, errors::Array{T,1} = zeros(2), energy::T = 1e9) where T <: Number
    k = 2*π*freq/2.99792458e8
    return TCRABCAVITY{NVAR(), T}(name, DTPSAD(len), DTPSAD(volt), DTPSAD(freq), DTPSAD(k), DTPSAD(phi), DTPSAD.(errors), DTPSAD(energy), "CRABCAVITY")
end
TCRABCAVITY(name::String, len::DTPSAD{N, T}, volt::DTPSAD{N, T}, freq::DTPSAD{N, T}, k::DTPSAD{N, T},
            phi::DTPSAD{N, T}, errors::Array{DTPSAD{N, T},1}, energy::DTPSAD{N, T}, eletype::String) where {N, T <: Number} =
    TCRABCAVITY{N, T}(name, len, volt, freq, k, phi, errors, energy, eletype)

mutable struct TCRABCAVITYK2{N, T} <: AbstractTPSAElement
    name::String 
    len::DTPSAD{N, T}
    volt::DTPSAD{N, T}  # voltage
    freq::DTPSAD{N, T}  # frequency
    k::DTPSAD{N, T}  # wave number
    phi::DTPSAD{N, T}  # phase
    errors::Array{DTPSAD{N, T},1} # 1: Voltage error, 2: Phase error
    k2::DTPSAD{N, T}
    energy::DTPSAD{N, T}
    eletype::String 
end
function TCRABCAVITYK2(;name::String = "CRABCAVITY", len::T = 0.0, volt::T = 0.0, 
    freq::T = 0.0, phi::T = 0.0, k2::T = 0.0, errors::Array{T,1} = zeros(2), energy::T = 1e9) where T <: Number
    k = 2*π*freq/2.99792458e8
    return TCRABCAVITYK2{NVAR(), T}(name, DTPSAD(len), DTPSAD(volt), DTPSAD(freq), DTPSAD(k), DTPSAD(phi), DTPSAD(k2), DTPSAD.(errors), DTPSAD(energy), "CRABCAVITY")
end
TCRABCAVITYK2(name::String, len::DTPSAD{N, T}, volt::DTPSAD{N, T}, freq::DTPSAD{N, T}, k::DTPSAD{N, T},
            phi::DTPSAD{N, T}, k2::DTPSAD{N, T}, errors::Array{DTPSAD{N, T},1}, energy::DTPSAD{N, T}, eletype::String) where {N, T <: Number} =
    TCRABCAVITYK2{N, T}(name, len, volt, freq, k, phi, k2, errors, energy, eletype)

mutable struct TthinMULTIPOLE{N, T} <: AbstractTPSAElement
    name::String
    len::DTPSAD{N, T}
    PolynomA::Array{DTPSAD{N, T},1}
    PolynomB::Array{DTPSAD{N, T},1}
    MaxOrder::Int64
    NumIntSteps::Int64
    rad::Int64
    FringeQuadEntrance::Int64
    FringeQuadExit::Int64
    T1::Array{DTPSAD{N, T},1}
    T2::Array{DTPSAD{N, T},1}
    R1::Array{DTPSAD{N, T},2}
    R2::Array{DTPSAD{N, T},2}
    RApertures::Array{Float64,1}
    EApertures::Array{Float64,1}
    KickAngle::Array{DTPSAD{N, T},1}
    eletype::String
end
function TthinMULTIPOLE(;name::String = "thinMULTIPOLE", len::T = 0.0, PolynomA::Array{T,1} = zeros(4), 
                PolynomB::Array{T,1} = zeros(4), MaxOrder::Int64=1, NumIntSteps::Int64 = 1, rad::Int64=0, 
                FringeQuadEntrance::Int64 = 0, FringeQuadExit::Int64 = 0, T1::Array{T,1} = zeros(6), 
                T2::Array{T,1} = zeros(6), R1::Array{T,2} = zeros(6, 6), 
                R2::Array{T,2} = zeros(6, 6), RApertures::Array{Float64,1} = zeros(Float64, 6), 
                EApertures::Array{Float64,1} = zeros(Float64, 6), KickAngle::Array{T,1} = zeros(2)) where T <: Number

    if PolynomB[3] != 0.0
        MaxOrder = 2
    end
    if PolynomB[4] != 0.0
        MaxOrder = 3
    end
    TthinMULTIPOLE{NVAR(), T}(name, DTPSAD(len), DTPSAD.(PolynomA), DTPSAD.(PolynomB), MaxOrder, NumIntSteps, rad, FringeQuadEntrance, FringeQuadExit,
        DTPSAD.(T1), DTPSAD.(T2), DTPSAD.(R1), DTPSAD.(R2), RApertures, EApertures, DTPSAD.(KickAngle), "thinMULTIPOLE")
end
TthinMULTIPOLE(name::String, len::DTPSAD{N, T}, PolynomA::Array{DTPSAD{N, T},1}, 
                PolynomB::Array{DTPSAD{N, T},1}, MaxOrder::Int64, NumIntSteps::Int64, rad::Int64, 
                FringeQuadEntrance::Int64, FringeQuadExit::Int64, T1::Array{DTPSAD{N, T},1}, 
                T2::Array{DTPSAD{N, T},1}, R1::Array{DTPSAD{N, T},2}, R2::Array{DTPSAD{N, T},2},
                RApertures::Array{Float64,1}, EApertures::Array{Float64,1}, KickAngle::Array{DTPSAD{N, T},1},
                eletype::String) where {N, T <: Number} =
    TthinMULTIPOLE{N, T}(name, len, PolynomA, PolynomB, MaxOrder, NumIntSteps, rad, FringeQuadEntrance, FringeQuadExit,
        T1, T2, R1, R2, RApertures, EApertures, KickAngle, eletype)

# TRANSLATION and YROTATION are used to convert the MAD-X lattice files
mutable struct TTRANSLATION{N, T} <: AbstractTPSAElement
    name::String
    len::DTPSAD{N, T}
    dx::DTPSAD{N, T}
    dy::DTPSAD{N, T}
    ds::DTPSAD{N, T}
    eletype::String
end
function TTRANSLATION(;name::String = "TRANSLATION", len::T = 0.0, dx::T = 0.0, dy::T = 0.0, ds::T = 0.0) where T <: Number
    TTRANSLATION{NVAR(), T}(name, DTPSAD(len), DTPSAD(dx), DTPSAD(dy), DTPSAD(ds), "TRANSLATION")
end
TTRANSLATION(name::String, len::DTPSAD{N, T}, dx::DTPSAD{N, T}, dy::DTPSAD{N, T}, ds::DTPSAD{N, T}, eletype::String) where {N, T <: Number} =
    TTRANSLATION{N, T}(name, len, dx, dy, ds, eletype)

mutable struct TYROTATION{N, T} <: AbstractTPSAElement
    name::String
    len::DTPSAD{N, T}
    angle::DTPSAD{N, T}
    eletype::String
end
function TYROTATION(;name::String = "YROTATION", len::T = 0.0, angle::T = 0.0) where T <: Number
    TYROTATION{NVAR(), T}(name, DTPSAD(len), DTPSAD(angle), "YROTATION")
end
TYROTATION(name::String, len::DTPSAD{N, T}, angle::DTPSAD{N, T}, eletype::String) where {N, T <: Number} =
    TYROTATION{N, T}(name, len, angle, eletype)

function strthinkick!(r::SubArray, A::Vector{DTPSAD{N, T}}, B::Vector{DTPSAD{N, T}}, L::DTPSAD{N, T}, max_order::Int) where {N, T <: Number}
    # Modified based on AT function. Ref[Terebilo, Andrei. "Accelerator modeling with MATLAB accelerator toolbox." PACS2001 (2001)].

    ReSum = B[max_order + 1]
    ImSum = A[max_order + 1]
    ReSumTemp = 0.0

    for i in max_order-1: -1: 0
        ReSumTemp = ReSum * r[1] - ImSum * r[3] + B[i+1]
        ImSum = ImSum * r[1] + ReSum * r[3] + A[i+1]
        ReSum = ReSumTemp
    end

    r[2] -= L * ReSum
    r[4] += L * ImSum
    return nothing
end

function strthinkick1!(r::SubArray, A::Vector{DTPSAD{N, T}}, B::Vector{DTPSAD{N, T}}, L::DTPSAD{N, T}, max_order::Int) where {N, T <: Number}
    # Modified based on AT function. Ref[Terebilo, Andrei. "Accelerator modeling with MATLAB accelerator toolbox." PACS2001 (2001)].

    # Calculate and apply a multipole kick to a 6-dimentional
    # phase space vector in a straight element (quadrupole)
    
    # IMPORTANT !!!
    # The reference coordinate system is straight but the field expansion may still
    # contain dipole terms A[1], B[1]

    ReSum = B[max_order + 1]
    ImSum = A[max_order + 1]
    ReSumTemp = 0.0

    for i in reverse(1:max_order)
        ReSumTemp = ReSum * r[1] - ImSum * r[3] + B[i]
        ImSum = ImSum * r[1] + ReSum * r[3] + A[i]
        ReSum = ReSumTemp
    end

    r[2] -= L * ReSum
    r[4] += L * ImSum
    return nothing
end

function bndthinkick!(r::SubArray, A::Array{DTPSAD{N, T},1}, B::Array{DTPSAD{N, T},1}, 
    L::DTPSAD{N, T}, irho::DTPSAD{N, T}, max_order::Int, beti::Float64) where {N, T <: Number}
    # AT function. Ref[Terebilo, Andrei. "Accelerator modeling with MATLAB accelerator toolbox." PACS2001 (2001)].
    ReSum = B[max_order + 1]
    ImSum = A[max_order + 1]
    ReSumTemp = 0.0

    for i in max_order-1:-1:0
        ReSumTemp = ReSum * r[1] - ImSum * r[3] + B[i+1]
        ImSum = ImSum * r[1] + ReSum * r[3] + A[i+1]
        ReSum = ReSumTemp
    end

    r[2] -= L * (ReSum - (r[6] - r[1] * irho) * irho)
    r[4] += L * ImSum
    r[5] += L * irho * r[1] * beti # Path length
    return nothing
end

function bndthinkickrad!(r::SubArray, A::Array{DTPSAD{N, T},1}, B::Array{DTPSAD{N, T},1}, 
    L::DTPSAD{N, T}, irho::DTPSAD{N, T}, E0::DTPSAD{N, T}, max_order::Int, beti::Float64) where {N, T <: Number}
    # Modified based on AT function. Ref[Terebilo, Andrei. "Accelerator modeling with MATLAB accelerator toolbox." PACS2001 (2001)].
    ReSum = B[max_order + 1]
    ImSum = A[max_order + 1]
    ReSumTemp = 0.0
    CRAD = CGAMMA * E0^3 / (2.0*pi*1e27) # [m]/[GeV^3] M.Sands (4.1)

    for i in max_order-1:-1:0
        ReSumTemp = ReSum * r[1] - ImSum * r[3] + B[i+1]
        ImSum = ImSum * r[1] + ReSum * r[3] + A[i+1]
        ReSum = ReSumTemp
    end

    # angles from momentums
    p_norm = 1.0 / (1.0*beti + r[6])
    x = r[1]
    xpr = r[2] * p_norm
    y = r[3]
    ypr = r[4] * p_norm
    B2P = B2perp(ImSum, ReSum + irho, irho, x, xpr, y, ypr)

    dp_0 = r[6]
    r[6] = r[6] - CRAD * (1.0*beti+r[6])^2 * B2P * (1.0 + x*irho + (xpr^2 + ypr^2) / 2.0) * L
    
    # momentums after losing energy
    p_norm = 1.0 / (1.0*beti + r[6])
    r[2] = xpr / p_norm
    r[4] = ypr / p_norm

    r[2] -= L * (ReSum - (dp_0 - r[1]*irho)*irho)
    r[4] += L * ImSum
    r[5] += L * irho * r[1] * beti # Path length
    return nothing
end

function Yrot!(r::SubArray, phi::DTPSAD{N, T}, beti::Float64) where {N, T <: Number}
    # Forest 10.26, rotation in free space
    if phi != 0.0
        dp1 = 1.0*beti + r[6]
        c = cos(phi)
        s = sin(phi)
        pz = pxyz(dp1, r[2], r[4])
        p = c * pz - s * r[2]
        px = s * pz + c * r[2]
        x = r[1] * pz / p
        dy = r[1] * r[4] * s / p
        dct = dp1 * r[1] * s / p
        r[1] = x
        r[2] = px
        r[3] += dy
        r[5] += dct
    end
    return nothing
end

function Sec(x::DTPSAD{N, T}) where {N, T <: Number}
    return 1.0 / cos(x)
end

function bend_fringe!(r::SubArray, irho::DTPSAD{N, T}, gK::DTPSAD{N, T}, beti::Float64) where {N, T <: Number}
    # Forest 13.13, bend fringe in the hard-edge limit
    b0 = irho 
    pz = pxyz(1.0*beti+r[6], r[2], r[4])
    px = r[2]
    py = r[4]
    d = r[6]
    xp = px / pz
    yp = py / pz
    phi = -b0 * tan( b0 * gK * (1.0 + xp*xp*(2.0 + yp*yp))*pz - atan(xp / (1.0 + yp*yp)))
    px2 = px*px
    px4 = px2*px2
    py2 = py*py
    py4 = py2*py2
    py6 = py4*py2
    pz2 = pz*pz
    pz3 = pz2*pz
    pz4 = pz2*pz2
    pz5 = pz4*pz
    pz6 = pz4*pz2
    py2z2 = (py2 + pz2) * (py2 + pz2)
    # powsec = pow(Sec((b0*gK*(pz4 + px2*(py2 + 2*pz2)))/pz3 - atan((px*pz)/(py2 + pz2))),2)
    powsec = Sec((b0*gK*(pz4 + px2*(py2 + 2.0*pz2)))/pz3 - atan((px*pz)/(py2 + pz2)))^2
    denom = (pz5*(py4 + px2*pz2 + 2.0*py2*pz2 + pz4))

    dpx = -(b0*(px2*pz4*(py2 - pz2) - pz6*(py2 + pz2) +
            b0*gK*px*(pz2*py2z2*(2.0*py2 + 3.0*pz2) + px4*(3.0*py2*pz2 + 2.0*pz4) +
            px2*(3.0*py6 + 8.0*py4*pz2 + 9.0*py2*pz4 + 5.0*pz6)))*powsec)/denom
    dpy = -(b0*py*(px*pz4*(py2 + pz2) +
            b0*gK*(-(pz4*py2z2) + px4*(3.0*py2*pz2 + 4.0*pz4) +
                   px2*(3.0*py6 + 10.0*py4*pz2 + 11.0*py2*pz4 + 3.0*pz6)))*powsec)/denom
    dd = (b0*(1.0*beti + d)*(px*pz4*(py2 - pz2) + b0*gK*
                      (-(pz4*py2z2) + px4*(3.0*py2*pz2 + 2.0*pz4) +
                       px2*(3.0*py6 + 8.0*py4*pz2 + 7.0*py2*pz4 + pz6)))*powsec)/denom

    yf = (2.0 * r[3]) / (1.0 + sqrt(1.0 - 2.0 * dpy * r[3]))
    dxf = 0.5 * dpx * yf * yf
    dct = 0.5 * dd * yf * yf
    dpyf = phi * yf

    r[3] = yf
    r[1] += dxf
    r[4] -= dpyf
    r[5] -= dct
    return nothing
end

function multipole_fringe!(r6::SubArray, le::DTPSAD{N, T}, polya::Array{DTPSAD{N, T},1}, polyb::Array{DTPSAD{N, T},1}, 
        max_order::Int, edge::DTPSAD{N, T}, skip_b0::Int, beti::Float64) where {N, T <: Number}
    # Forest 13.29
    FX = 0.0
    FY = 0.0
    FX_X = 0.0
    FX_Y = 0.0
    FY_X = 0.0
    FY_Y = 0.0
  
    RX = 1.0
    IX = 0.0

    for n in 0:max_order
        B = polyb[n + 1]  
        A = polya[n + 1] 
    
        j = n + 1.0
    
        DRX = RX
        DIX = IX
    
        # Complex multiplications
        RX = DRX * r6[1] - DIX * r6[3]
        IX = DRX * r6[3] + DIX * r6[1]
        
        U, V, DU, DV = 0.0, 0.0, 0.0, 0.0
        if n == 0 && skip_b0 != 0
            U -= A * IX
            V += A * RX
            DU -= A * DIX
            DV += A * DRX
        else
            U += B * RX - A * IX
            V += B * IX + A * RX
            DU += B * DRX - A * DIX
            DV += B * DIX + A * DRX
        end
    
        f1 = -edge / 4.0 / (j + 1.0)
    
        U *= f1
        V *= f1
        DU *= f1
        DV *= f1
    
        DUX = j * DU
        DVX = j * DV
        DUY = -j * DV
        DVY = j * DU
    
        nf = 1.0 * (j + 2.0) / j
    
        FX += U * r6[1] + nf * V * r6[3]
        FY += U * r6[3] - nf * V * r6[1]
    
        FX_X += DUX * r6[1] + U + nf * r6[3] * DVX
        FX_Y += DUY * r6[1] + nf * V + nf * r6[3] * DVY
    
        FY_X += DUX * r6[3] - nf * V - nf * r6[1] * DVX
        FY_Y += DUY * r6[3] + U - nf * r6[1] * DVY
    end

    DEL = 1.0 / (1.0*beti + r6[6])
    A = 1.0 - FX_X * DEL
    B = -FY_X * DEL
    D = 1.0 - FY_Y * DEL
    C = -FX_Y * DEL

    r6[1] -= FX * DEL
    r6[3] -= FY * DEL

    pxf = (D * r6[2] - B * r6[4]) / (A * D - B * C)
    pyf = (A * r6[4] - C * r6[2]) / (A * D - B * C)
    r6[4] = pyf
    r6[2] = pxf
    r6[5] -= (r6[2] * FX + r6[4] * FY) * DEL * DEL
    return nothing
end

function exact_bend!(r6::SubArray, irho::DTPSAD{N, T}, L::DTPSAD{N, T}, beti::Float64) where {N, T <: Number}
    # Forest 12.18, bend-kick split, map W(L, irho)

    dp1 = 1.0*beti + r6[6]  # r6[delta_]
    pz = pxyz(dp1, r6[2], r6[4])  # r6[px_], r6[py_]

    if abs(irho) < 1e-6
        NormL = L / pz
        r6[1] += r6[2] * NormL  # r6[x_]
        r6[3] += r6[4] * NormL  # r6[y_]
        r6[5] += NormL * dp1    # r6[ct_], absolute path length
    else
        pzmx = pz - (1.0 + r6[1] * irho)  # r6[x_]
        cs = cos(irho * L)
        sn = sin(irho * L)
        px = r6[2] * cs + pzmx * sn  # r6[px_]
        d2 = pxyz(dp1, 0.0, r6[4])  # r6[py_]
        dasin = L + (asin(r6[2] / d2) - asin(px / d2)) / irho
        x = (pxyz(dp1, px, r6[4]) - pzmx * cs + r6[2] * sn - 1.0) / irho  # r6[x_]
        dy = r6[4] * dasin  # r6[py_]
        dct = dp1 * dasin   # r6[ct_], absolute path length

        r6[1] = x
        r6[2] = px
        r6[3] += dy
        r6[5] += dct
    end
    return nothing
end

function exactbndthinkick_rad!(r::SubArray, A::Array{DTPSAD{N, T},1}, B::Array{DTPSAD{N, T},1}, 
    L::DTPSAD{N, T}, irho::DTPSAD{N, T}, max_order::Int, beti::Float64, rad_const::DTPSAD{N, T}) where {N, T <: Number}
    # AT function. Ref[Terebilo, Andrei. "Accelerator modeling with MATLAB accelerator toolbox." PACS2001 (2001)].
    ReSum = B[max_order + 1]
    ImSum = A[max_order + 1]
    ReSumTemp = 0.0

    for i in max_order-1:-1:0
        ReSumTemp = ReSum * r[1] - ImSum * r[3] + B[i+1]
        ImSum = ImSum * r[1] + ReSum * r[3] + A[i+1]
        ReSum = ReSumTemp
    end

    p_norm = 1.0 / (1.0 + r[6])
    x = r[1]
    xpr = r[2] * p_norm
    y = r[3]
    ypr = r[4] * p_norm

    B2P = B2perp_exact_bnd(ImSum, ReSum+irho, irho, x, xpr, y, ypr)

    r[6] -= rad_const * (1.0 + r[6])^2 * B2P * (1.0 + x*irho) * L / sqrt(1.0 - xpr*xpr - ypr*ypr)

    p_norm = 1.0 / (1.0 + r[6])
    r[2] = xpr/p_norm
    r[4] = ypr/p_norm

    r[2] -= L * ReSum
    r[4] += L * ImSum
    return nothing
end

function bend_edge!(r6::SubArray, rhoinv::DTPSAD{N, T}, theta::DTPSAD{N, T}, beti::Float64) where {N, T <: Number}
    # Forest 12.41, ideal wedge, map U(theta, rhoinv)

    if abs(rhoinv) >= 1e-6
        dp1 = 1.0*beti + r6[6]  # r6[delta_]
        c = cos(theta)
        s = sin(theta)
        pz = pxyz(dp1, r6[2], r6[4])  # r6[px_], r6[py_]
        d2 = pxyz(dp1, 0.0, r6[4])    # r6[py_]
        px = r6[2] * c + (pz - rhoinv * r6[1]) * s  # r6[px_]
        dasin = asin(r6[2] / d2) - asin(px / d2)
        num = r6[1] * (r6[2] * sin(2.0 * theta) + s * s * (2.0 * pz - rhoinv * r6[1]))
        den = pxyz(dp1, px, r6[4]) + pxyz(dp1, r6[2], r6[4]) * c - r6[2] * s
        x = r6[1] * c + num / den  # r6[x_]
        dy = r6[4] * theta / rhoinv + r6[4] / rhoinv * dasin  # r6[py_]
        dct = dp1 / rhoinv * (theta + dasin)  # r6[ct_]

        r6[1] = x
        r6[2] = px
        r6[3] += dy
        r6[5] += dct
    end
    return nothing
end

function ExactSectorBend!(r::Array{DTPSAD{N, T},1}, le::DTPSAD{N, T}, beti::Float64, angle::DTPSAD{N, T}, A::Array{DTPSAD{N, T},1}, B::Array{DTPSAD{N, T},1},
    max_order::Int, num_int_steps::Int, entrance_angle::DTPSAD{N, T}, exit_angle::DTPSAD{N, T}, FringeBendEntrance::Int, FringeBendExit::Int,
    FringeQuadEntrance::Int, FringeQuadExit::Int, gk::DTPSAD{N, T},
    T1::Array{DTPSAD{N, T},1}, T2::Array{DTPSAD{N, T},1},
    R1::Array{DTPSAD{N, T},2}, R2::Array{DTPSAD{N, T},2}, RApertures::Array{Float64,1}, EApertures::Array{Float64,1},
    KickAngle::Array{DTPSAD{N, T},1}, num_particles::Int, lost_flags::Array{Int64,1}) where {N, T <: Number}

    irho = angle / le
    DRIFT1 = 0.6756035959798286638
    DRIFT2 = -0.1756035959798286639
    KICK1 = 1.351207191959657328
    KICK2 = -1.702414383919314656
    SL = le / num_int_steps
    L1 = SL * DRIFT1
    L2 = SL * DRIFT2
    K1 = SL * KICK1
    K2 = SL * KICK2

    B0 = B[1]
    A0 = A[1]

    if !iszero(KickAngle[1])
        B[1] -= sin(KickAngle[1]) / le
    end
    if !iszero(KickAngle[2])
        A[1] += sin(KickAngle[2]) / le
    end

    @inbounds for c in 1:num_particles
        if isone(lost_flags[c])
            continue
        end
        r6 = @view r[(c-1)*6+1:c*6]
        # Misalignment at entrance
        # Use `all(iszero, ...)` to test if arrays are all zero.
        if !all(iszero, T1)
            addvv!(r6, T1)
        end
        if !all(iszero, R1)
            multmv!(r6, R1)
        end

        Yrot!(r6, entrance_angle, beti)

        if FringeBendEntrance != 0
            bend_fringe!(r6, irho, gk, beti)
        end

        if FringeQuadEntrance != 0
            multipole_fringe!(r6, le, A, B, max_order, 1.0, 1, beti)
        end

        bend_edge!(r6, irho, -entrance_angle, beti)

        # Integrator
        if num_int_steps == 0
            exact_bend!(r6, irho, le, beti)
        else
            for m in 1:num_int_steps
                exact_bend!(r6, irho, L1, beti)
                strthinkick!(r6, A, B, K1, max_order)
                exact_bend!(r6, irho, L2, beti)
                strthinkick!(r6, A, B, K2, max_order)
                exact_bend!(r6, irho, L2, beti)
                strthinkick!(r6, A, B, K1, max_order)
                exact_bend!(r6, irho, L1, beti)
            end
        end

        r6[5] -= le*beti  # r6[ct_], absolute path length

        bend_edge!(r6, irho, -exit_angle, beti)
        if FringeQuadExit != 0
            multipole_fringe!(r6, le, A, B, max_order, -1.0, 1, beti)
        end
        if FringeBendExit != 0
            bend_fringe!(r6, -irho, gk, beti)
        end
        Yrot!(r6, exit_angle, beti)

        # Misalignment at exit
        if !all(iszero, R2)
            multmv!(r6, R2)
        end
        if !all(iszero, T2)
            addvv!(r6, T2)
        end
        if check_lost_GTPSA(r6)
            lost_flags[c] = 1
        end
    end
    
    if !iszero(KickAngle[1])
        B[1] += sin(KickAngle[1]) / le
    end
    if !iszero(KickAngle[2])
        A[1] -= sin(KickAngle[2]) / le
    end
    return nothing
end

function ExactSectorBend_rad!(r::Array{DTPSAD{N, T},1}, le::DTPSAD{N, T}, rad_const::DTPSAD{N, T}, beti::Float64, angle::DTPSAD{N, T}, A::Array{DTPSAD{N, T},1}, B::Array{DTPSAD{N, T},1}, 
    max_order::Int, num_int_steps::Int, entrance_angle::DTPSAD{N, T}, exit_angle::DTPSAD{N, T}, FringeBendEntrance::Int, FringeBendExit::Int,
    FringeQuadEntrance::Int, FringeQuadExit::Int, gk::DTPSAD{N, T},
    T1::Array{DTPSAD{N, T},1}, T2::Array{DTPSAD{N, T},1}, 
    R1::Array{DTPSAD{N, T},2}, R2::Array{DTPSAD{N, T},2}, RApertures::Array{Float64,1}, EApertures::Array{Float64,1},
    KickAngle::Array{DTPSAD{N, T},1}, num_particles::Int, lost_flags::Array{Int64,1}) where {N, T <: Number}
    
    irho = angle / le
    DRIFT1 = 0.6756035959798286638
    DRIFT2 = -0.1756035959798286639
    KICK1 = 1.351207191959657328
    KICK2 = -1.702414383919314656
    SL = le / num_int_steps
    L1 = SL * DRIFT1
    L2 = SL * DRIFT2
    K1 = SL * KICK1
    K2 = SL * KICK2

    B0 = B[1]
    A0 = A[1]

    if !iszero(KickAngle[1])
        B[1] -= sin(KickAngle[1]) / le
    end
    if !iszero(KickAngle[2])
        A[1] += sin(KickAngle[2]) / le
    end

    @inbounds for c in 1:num_particles
        if isone(lost_flags[c])
            continue
        end
        r6 = @view r[(c-1)*6+1:c*6]
        # Misalignment at entrance
        if !all(iszero, T1)
            addvv!(r6, T1)
        end
        if !all(iszero, R1)
            multmv!(r6, R1)
        end

        Yrot!(r6, entrance_angle, beti)

        if FringeBendEntrance != 0
            bend_fringe!(r6, irho, gk, beti)
        end

        if FringeQuadEntrance != 0
            multipole_fringe!(r6, le, A, B, max_order, 1.0, 1, beti)
        end

        bend_edge!(r6, irho, -entrance_angle, beti)

        # Integrator
        if num_int_steps == 0
            exact_bend!(r6, irho, le, beti)
        else
            for m in 1:num_int_steps
                exact_bend!(r6, irho, L1, beti)
                exactbndthinkick_rad!(r6, A, B, K1, irho, max_order, beti, rad_const)
                exact_bend!(r6, irho, L2, beti)
                exactbndthinkick_rad!(r6, A, B, K2, irho, max_order, beti, rad_const)
                exact_bend!(r6, irho, L2, beti)
                exactbndthinkick_rad!(r6, A, B, K1, irho, max_order, beti, rad_const)
                exact_bend!(r6, irho, L1, beti)
            end
        end

        r6[5] -= le*beti  # r6[ct_], absolute path length

        bend_edge!(r6, irho, -exit_angle, beti)
        if FringeQuadExit != 0
            multipole_fringe!(r6, le, A, B, max_order, -1.0, 1, beti)
        end
        if FringeBendExit != 0
            bend_fringe!(r6, -irho, gk, beti)
        end
        Yrot!(r6, exit_angle, beti)

        # Misalignment at exit
        if !all(iszero, R2)
            multmv!(r6, R2)
        end
        if !all(iszero, T2)
            addvv!(r6, T2)
        end
        if check_lost_GTPSA(r6)
            lost_flags[c] = 1
        end
    end
    
    if !iszero(KickAngle[1])
        B[1] += sin(KickAngle[1]) / le
    end
    if !iszero(KickAngle[2])
        A[1] -= sin(KickAngle[2]) / le
    end
    return nothing
end

function multmv!(r::SubArray, A::Matrix{DTPSAD{N, T}}) where {N, T <: Number}
    temp = zeros(DTPSAD{N, T}, 6)
    for i in 1:6
        for j in 1:6
            temp[i] += A[i, j] * r[j]
        end
    end
    r .= temp
    return nothing
end

function addvv!(r::SubArray, dr::Array{DTPSAD{N, T},1}) where {N, T <: Number}
    r .= r .+ dr
    return nothing
end

function check_lost_GTPSA(r6)
    if isnan(r6[1]) || isinf(r6[1])
        return true
    end
    if maximum(abs.(r6[1:4])) > CoordLimit || abs(r6[6]) > CoordLimit
        return true
    end
    # sqrt(1.0 + 2.0*r[6]*beti + r[6]^2 - r[2]^2 - r[4]^2) must be real
    if r6[2]^2 + r6[4]^2 > 1.0  + r6[6]^2
        return true
    end
    return false
end

function drift6!(r::SubArray, le::DTPSAD{N, T}) where {N, T <: Number}
    if isone(use_exact_Hamiltonian)
        NormL = le / sqrt(1.0 + 2.0*r[6] + r[6]^2 - r[2]^2 - r[4]^2)
        r[5] += NormL * (1.0 + r[6]) - le
    else
        NormL = le / (1.0 + r[6])
        r[5] += NormL * (r[2]^2 + r[4]^2) / (2.0*(1.0+r[6])) # for linearized approximation
    end
    r[1] += NormL * r[2]
    r[3] += NormL * r[4]
    return nothing
end 

function pass!(ele::TDRIFT, rin::Vector{DTPSAD{N, T}}, npart::Int, particles::TBeam) where {N, T <: Number}
    @inbounds for c in 1:npart
        lost_flag = particles.lost_flag[c]
        if lost_flag == 1
            continue
        end
        r6 = @view rin[(c-1)*6+1:c*6]  
        if !all(iszero, ele.T1)
            addvv!(r6, ele.T1)
        end
        if !all(iszero, ele.R1)
            multmv!(r6, ele.R1)
        end
        drift6!(r6, ele.len)
        if !all(iszero, ele.R2)
            multmv!(r6, ele.R2)
        end
        if !all(iszero, ele.T2)
            addvv!(r6, ele.T2)
        end
        lost_flag = check_lost_GTPSA(r6) ? 1 : lost_flag
        particles.lost_flag[c] = lost_flag
    end
    return nothing
end

function pass!(ele::TQUAD, rin::Vector{DTPSAD{N, T}}, npart::Int, particles::TBeam) where {N, T <: Number}
    @inbounds for c in 1:npart
        lost_flag = particles.lost_flag[c]
        if lost_flag == 1
            continue
        end
        r6 = @view rin[(c-1)*6+1:c*6]  
        p_norm = 1.0 / (1.0 + r6[6])
        
        if !all(iszero, ele.T1)
            addvv!(r6, ele.T1)
        end
        if !all(iszero, ele.R1)
            multmv!(r6, ele.R1)
        end
        
        if iszero(ele.k1)
            drift6!(r6, ele.len)
        else
            x   = r6[1]
            xpr = r6[2] * p_norm
            y   = r6[3]
            ypr = r6[4] * p_norm
            
            g  = abs(ele.k1) / (1.0 + r6[6])
            t  = sqrt(g)
            lt = ele.len * t
            
            if ele.k1 > 0
                MHD = cos(lt)
                M12 = sin(lt)/t
                M21 = -M12*g
                MVD = cosh(lt)
                M34 = sinh(lt)/t
                M43 = M34*g
            else
                MHD = cosh(lt)
                M12 = sinh(lt)/t
                M21 = M12*g
                MVD = cos(lt)
                M34 = sin(lt)/t
                M43 = -M34*g
            end
            
            r6[1] =  MHD*x + M12*xpr
            r6[2] = (M21*x + MHD*xpr)/p_norm
            r6[3] =  MVD*y + M34*ypr
            r6[4] = (M43*y + MVD*ypr)/p_norm

            r6[5]+= g*(x*x*(ele.len-MHD*M12)-y*y*(ele.len-MVD*M34))/4.0
            r6[5]+= (xpr*xpr*(ele.len+MHD*M12)+ypr*ypr*(ele.len+MVD*M34))/4.0
            r6[5]+= (x*xpr*M12*M21 + y*ypr*M34*M43)/2.0
        end

        if !all(iszero, ele.R2)
            multmv!(r6, ele.R2)
        end
        if !all(iszero, ele.T2)
            addvv!(r6, ele.T2)
        end
        lost_flag = check_lost_GTPSA(r6) ? 1 : lost_flag
        particles.lost_flag[c] = lost_flag
    end
    return nothing
end

function pass!(ele::TCORRECTOR, rin::Vector{DTPSAD{N, T}}, npart::Int, particles::TBeam) where {N, T <: Number}
    @inbounds for c in 1:npart
        lost_flag = particles.lost_flag[c]
        if lost_flag == 1
            continue
        end
        r6 = @view rin[(c-1)*6+1:c*6]  
        if !all(iszero, ele.T1)
            addvv!(r6, ele.T1)
        end
        if !all(iszero, ele.R1)
            multmv!(r6, ele.R1)
        end

        p_norm = 1.0 / (1.0 + r6[6])
        NormL = ele.len * p_norm

        r6[5] += NormL*p_norm*(ele.xkick*ele.xkick/3.0 + ele.ykick*ele.ykick/3.0 +
   		            r6[2]*r6[2] + r6[4]*r6[4] +
   		            r6[2]*ele.xkick + r6[4]*ele.ykick)/2.0
        r6[1] += NormL*(r6[2]+ele.xkick/2.0)
        r6[2] += ele.xkick
        r6[3] += NormL*(r6[4]+ele.ykick/2.0)
        r6[4] += ele.ykick

        if !all(iszero, ele.R2)
            multmv!(r6, ele.R2)
        end
        if !all(iszero, ele.T2)
            addvv!(r6, ele.T2)
        end
        lost_flag = check_lost_GTPSA(r6) ? 1 : lost_flag
        particles.lost_flag[c] = lost_flag
    end
    return nothing
end

function pass!(ele::TMARKER, rin::Vector{DTPSAD{N, T}}, npart::Int, particles::TBeam) where {N, T <: Number}
    return nothing
end

function pass!(ele::TSOLENOID, rin::Vector{DTPSAD{N, T}}, npart::Int, particles::TBeam) where {N, T <: Number}
    @inbounds for c in 1:npart
        lost_flag = particles.lost_flag[c]
        if lost_flag == 1
            continue
        end
        r6 = @view rin[(c-1)*6+1:c*6]  

        if !all(iszero, ele.T1)
            addvv!(r6, ele.T1)
        end
        if !all(iszero, ele.R1)
            multmv!(r6, ele.R1)
        end

        if ele.ks != 0.0
            p_norm = 1.0 / (1.0 + r6[6])
            x = r6[1]
            xpr = r6[2]*p_norm
            y = r6[3]
            ypr = r6[4]*p_norm

            H = ele.ks * p_norm / 2.0
            S = sin(ele.len * H)
            C = cos(ele.len * H)
            r6[1] = x*C*C + xpr*C*S/H + y*C*S + ypr*S*S/H
            r6[2] = (-x*H*C*S + xpr*C*C - y*H*S*S + ypr*C*S) / p_norm
            r6[3] = -x*C*S - xpr*S*S/H + y*C*C + ypr*C*S/H
            r6[4] = (x*H*S*S - xpr*C*S - y*C*S*H + ypr*C*C) / p_norm
            r6[5] += ele.len*(H*H*(x*x+y*y) + 2.0*H*(xpr*y-ypr*x) +xpr*xpr+ypr*ypr)/2.0
        else
            drift6!(r6, ele.len)
        end
        if !all(iszero, ele.R2)
            multmv!(r6, ele.R2)
        end
        if !all(iszero, ele.T2)
            addvv!(r6, ele.T2)
        end
        lost_flag = check_lost_GTPSA(r6) ? 1 : lost_flag
        particles.lost_flag[c] = lost_flag
    end
    return nothing
end

function multipole_fringe!(r6::SubArray, le::DTPSAD{N, T}, polya::Array{DTPSAD{N, T},1}, polyb::Array{DTPSAD{N, T},1}, 
        max_order::Int, edge::Float64, skip_b0::Int, beti::Float64) where {N, T <: Number}
    # Forest 13.29
    FX = 0.0
    FY = 0.0
    FX_X = 0.0
    FX_Y = 0.0
    FY_X = 0.0
    FY_Y = 0.0
  
    RX = 1.0
    IX = 0.0

    for n in 0:max_order
        B = polyb[n + 1]  
        A = polya[n + 1] 
    
        j = n + 1.0
    
        DRX = RX
        DIX = IX
    
        # Complex multiplications
        RX = DRX * r6[1] - DIX * r6[3]
        IX = DRX * r6[3] + DIX * r6[1]
        
        U, V, DU, DV = 0.0, 0.0, 0.0, 0.0
        if n == 0 && skip_b0 != 0
            U -= A * IX
            V += A * RX
            DU -= A * DIX
            DV += A * DRX
        else
            U += B * RX - A * IX
            V += B * IX + A * RX
            DU += B * DRX - A * DIX
            DV += B * DIX + A * DRX
        end
    
        f1 = -edge / 4.0 / (j + 1.0)
    
        U *= f1
        V *= f1
        DU *= f1
        DV *= f1
    
        DUX = j * DU
        DVX = j * DV
        DUY = -j * DV
        DVY = j * DU
    
        nf = 1.0 * (j + 2.0) / j
    
        FX += U * r6[1] + nf * V * r6[3]
        FY += U * r6[3] - nf * V * r6[1]
    
        FX_X += DUX * r6[1] + U + nf * r6[3] * DVX
        FX_Y += DUY * r6[1] + nf * V + nf * r6[3] * DVY
    
        FY_X += DUX * r6[3] - nf * V - nf * r6[1] * DVX
        FY_Y += DUY * r6[3] + U - nf * r6[1] * DVY
    end

    DEL = 1.0 / (1.0*beti + r6[6])
    A = 1.0 - FX_X * DEL
    B = -FY_X * DEL
    D = 1.0 - FY_Y * DEL
    C = -FX_Y * DEL

    r6[1] -= FX * DEL
    r6[3] -= FY * DEL

    pxf = (D * r6[2] - B * r6[4]) / (A * D - B * C)
    pyf = (A * r6[4] - C * r6[2]) / (A * D - B * C)
    r6[4] = pyf
    r6[2] = pxf
    r6[5] -= (r6[2] * FX + r6[4] * FY) * DEL * DEL
    return nothing
end

function StrB2perp(bx::DTPSAD{N, T}, by::DTPSAD{N, T}, x::DTPSAD{N, T}, xpr::DTPSAD{N, T}, y::DTPSAD{N, T}, ypr::DTPSAD{N, T}) where {N, T <: Number}
    # Calculates sqr(|B x e|) , where e is a unit vector in the direction of velocity
    # v_norm2 = 1.0 / (1.0 + xpr^2 + ypr^2)
    # return (by^2 + bx^2 + (bx*ypr - by*xpr)^2) * v_norm2
    return bx*bx + by*by + (bx*xpr - by*ypr)^2
end

function strthinkickrad!(r::SubArray, A::Vector{DTPSAD{N, T}}, B::Vector{DTPSAD{N, T}}, 
    L::DTPSAD{N, T}, E0::DTPSAD{N, T}, max_order::Int, rad_const::DTPSAD{N, T}) where {N, T <: Number}
    # Modified based on AT function. Ref[Terebilo, Andrei. "Accelerator modeling with MATLAB accelerator toolbox." PACS2001 (2001)].
    ReSum = B[max_order + 1]
    ImSum = A[max_order + 1]
    ReSumTemp = 0.0

    for i in reverse(1:max_order)
        ReSumTemp = ReSum * r[1] - ImSum * r[3] + B[i]
        ImSum = ImSum * r[1] + ReSum * r[3] + A[i]
        ReSum = ReSumTemp
    end

    # angles for momentum
    p_norm = 1.0 / (1.0 + r[6])
    x = r[1]
    xpr = r[2] * p_norm
    y = r[3]
    ypr = r[4] * p_norm
    B2P = StrB2perp(ImSum, ReSum , x , xpr, y ,ypr)
    factor = L / (p_norm)^2 / sqrt(1.0 - xpr^2 - ypr^2)

    dp_0 = r[6]
    r[6] -= rad_const * B2P * factor

    # momentums after losing energy
    p_norm = 1.0 / (1.0 + r[6])

    r[2] = xpr / p_norm
    r[4] = ypr / p_norm

    r[2] -= L * ReSum
    r[4] += L * ImSum
    return nothing
end

function StrMPoleSymplectic4Pass!(r::Vector{DTPSAD{N, T}}, le::DTPSAD{N, T}, beti::Float64, A::Array{DTPSAD{N, T},1}, B::Array{DTPSAD{N, T},1},
    max_order::Int, num_int_step::Int,
    FringeQuadEntrance::Int, FringeQuadExit::Int, # (no fringe), 1 (lee-whiting) or 2 (lee-whiting+elegant-like)
    T1::Array{DTPSAD{N, T},1}, T2::Array{DTPSAD{N, T},1}, R1::Array{DTPSAD{N, T},2}, R2::Array{DTPSAD{N, T},2},
    RApertures::Array{Float64,1}, EApertures::Array{Float64,1}, KickAngle::Array{DTPSAD{N, T},1},
    num_particles::Int, lost_flags::Array{Int64,1}) where {N, T <: Number}
    # Modified based on AT function. Ref[Terebilo, Andrei. "Accelerator modeling with MATLAB accelerator toolbox." PACS2001 (2001)].

    DRIFT1  =  0.6756035959798286638
    DRIFT2 = -0.1756035959798286639
    KICK1 = 1.351207191959657328
    KICK2 = -1.702414383919314656

    SL = le/num_int_step
    L1 = SL*DRIFT1
    L2 = SL*DRIFT2
    K1 = SL*KICK1
    K2 = SL*KICK2

    if le > 0
        B[1] -= sin(KickAngle[1])/le
        A[1] += sin(KickAngle[2])/le
    end
    @inbounds for c in 1:num_particles
        if lost_flags[c] == 1
            continue
        end
        r6 = @view r[(c-1)*6+1:c*6]
        # Misalignment at entrance
        if !all(iszero, T1)
            addvv!(r6, T1)
        end
        if !all(iszero, R1)
            multmv!(r6, R1)
        end

        if FringeQuadEntrance != 0 
            multipole_fringe!(r6, le, A, B, max_order, 1.0, 1, beti)
        end

        # Integrator
        for m in 1:num_int_step
            if check_lost_GTPSA(r6)
                lost_flags[c] = 1
                break
            end
            drift6!(r6, L1)
            strthinkick!(r6, A, B, K1, max_order)
            if check_lost_GTPSA(r6)
                lost_flags[c] = 1
                break
            end
            drift6!(r6, L2)
            strthinkick!(r6, A, B, K2, max_order)
            if check_lost_GTPSA(r6)
                lost_flags[c] = 1
                break
            end
            drift6!(r6, L2)
            strthinkick!(r6, A, B, K1, max_order)
            if check_lost_GTPSA(r6)
                lost_flags[c] = 1
                break
            end
            drift6!(r6, L1)
        end

        if FringeQuadExit != 0
            multipole_fringe!(r6, le, A, B, max_order, -1.0, 1, beti)
        end

        # Misalignment at exit
        if !all(iszero, R2)
            multmv!(r6, R2)
        end
        if !all(iszero, T2)
            addvv!(r6, T2)
        end
        if check_lost_GTPSA(r6)
            lost_flags[c] = 1
        end
    end
    if le > 0
        B[1] += sin(KickAngle[1]) / le
        A[1] -= sin(KickAngle[2]) / le
    end
    return nothing
end

function StrMPoleSymplectic4RadPass!(r::Vector{DTPSAD{N, T}}, le::DTPSAD{N, T}, rad_const::DTPSAD{N, T}, beti::Float64, A::Array{DTPSAD{N, T},1}, B::Array{DTPSAD{N, T},1},
    max_order::Int, num_int_step::Int,
    FringeQuadEntrance::Int, FringeQuadExit::Int, # (no fringe), 1 (lee-whiting) or 2 (lee-whiting+elegant-like)
    T1::Array{DTPSAD{N, T},1}, T2::Array{DTPSAD{N, T},1}, R1::Array{DTPSAD{N, T},2}, R2::Array{DTPSAD{N, T},2},
    RApertures::Array{Float64,1}, EApertures::Array{Float64,1}, KickAngle::Array{DTPSAD{N, T},1}, E0::DTPSAD{N, T},
    num_particles::Int, lost_flags::Array{Int64,1}) where {N, T <: Number}
    # Modified based on AT function. Ref[Terebilo, Andrei. "Accelerator modeling with MATLAB accelerator toolbox." PACS2001 (2001)].

    DRIFT1  =  0.6756035959798286638
    DRIFT2 = -0.1756035959798286639
    KICK1 = 1.351207191959657328
    KICK2 = -1.702414383919314656

    SL = le/num_int_step
    L1 = SL*DRIFT1
    L2 = SL*DRIFT2
    K1 = SL*KICK1
    K2 = SL*KICK2

    if le > 0
        B[1] -= sin(KickAngle[1])/le
        A[1] += sin(KickAngle[2])/le
    end
    @inbounds for c in 1:num_particles
        if lost_flags[c] == 1
            continue
        end
        r6 = @view r[(c-1)*6+1:c*6]
        # Misalignment at entrance
        if !all(iszero, T1)
            addvv!(r6, T1)
        end
        if !all(iszero, R1)
            multmv!(r6, R1)
        end

        if FringeQuadEntrance != 0 
            multipole_fringe!(r6, le, A, B, max_order, 1.0, 1, beti)
        end

        # Integrator
        for m in 1:num_int_step
            if check_lost_GTPSA(r6)
                lost_flags[c] = 1
                break
            end
            drift6!(r6, L1)
            strthinkickrad!(r6, A, B, K1, E0, max_order, rad_const)
            if check_lost_GTPSA(r6)
                lost_flags[c] = 1
                break
            end
            drift6!(r6, L2)
            strthinkickrad!(r6, A, B, K2, E0, max_order, rad_const)
            if check_lost_GTPSA(r6)
                lost_flags[c] = 1
                break
            end
            drift6!(r6, L2)
            strthinkickrad!(r6, A, B, K1, E0, max_order, rad_const)
            if check_lost_GTPSA(r6)
                lost_flags[c] = 1
                break
            end
            drift6!(r6, L1)
        end

        if FringeQuadExit != 0
            multipole_fringe!(r6, le, A, B, max_order, -1.0, 1, beti)
        end

        # Misalignment at exit
        if !all(iszero, R2)
            multmv!(r6, R2)
        end
        if !all(iszero, T2)
            addvv!(r6, T2)
        end
        if check_lost_GTPSA(r6)
            lost_flags[c] = 1
        end
    end
    if le > 0
        B[1] += sin(KickAngle[1]) / le
        A[1] -= sin(KickAngle[2]) / le
    end

    return nothing
end

function pass!(ele::TKQUAD, r_in::Vector{DTPSAD{N, T}}, num_particles::Int64, particles::TBeam) where {N, T <: Number}
    lost_flags = particles.lost_flag
    PolynomB = zeros(DTPSAD{N, T}, 4)
    E0 = particles.energy
    rad_const = DTPSAD(0.0)

    if ele.PolynomB[1] == 0.0 && ele.PolynomB[2] == 0.0 && ele.PolynomB[3] == 0.0 && ele.PolynomB[4] == 0.0
        PolynomB[2] = ele.k1
        if ele.rad == 0
            StrMPoleSymplectic4Pass!(r_in, ele.len, 1.0, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit, ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, num_particles, lost_flags)
        else
            if particles.mass == m_e
                rad_const = RAD_CONST_E * particles.gamma^3
            elseif particles.mass == m_p
                rad_const = RAD_CONST_P * particles.gamma^3
            else
                rad_const = DTPSAD(0.0)
                println("SR is not implemented for this particle mass.")
            end
            StrMPoleSymplectic4RadPass!(r_in, ele.len, rad_const, 1.0, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit, ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, E0, num_particles, lost_flags)
        end
    else
        PolynomB[1] = ele.PolynomB[1]
        PolynomB[2] = ele.PolynomB[2] 
        PolynomB[3] = ele.PolynomB[3] / 2.0
        PolynomB[4] = ele.PolynomB[4] / 6.0
        if ele.rad == 0
            StrMPoleSymplectic4Pass!(r_in, ele.len, 1.0, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit, 
                ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, num_particles, lost_flags)
        else
            if particles.mass == m_e
                rad_const = RAD_CONST_E * particles.gamma^3
            elseif particles.mass == m_p
                rad_const = RAD_CONST_P * particles.gamma^3
            else
                rad_const = DTPSAD(0.0)
                println("SR is not implemented for this particle mass.")
            end
            StrMPoleSymplectic4RadPass!(r_in, ele.len, rad_const, 1.0, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit, 
                ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, E0, num_particles, lost_flags)
        end
    end
    return nothing
end

function pass!(ele::TKSEXT, r_in::Vector{DTPSAD{N, T}}, num_particles::Int64, particles::TBeam) where {N, T <: Number}
    rad_const = DTPSAD(0.0)
    lost_flags = particles.lost_flag
    PolynomB = zeros(DTPSAD{N, T}, 4)
    E0 = particles.energy

    if ele.PolynomB[1] == 0.0 && ele.PolynomB[2] == 0.0 && ele.PolynomB[3] == 0.0 && ele.PolynomB[4] == 0.0
        PolynomB[3] = ele.k2 / 2.0
        if ele.rad == 0
            StrMPoleSymplectic4Pass!(r_in, ele.len, 1.0,ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit, ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, num_particles, lost_flags)
        else
            if particles.mass == m_e
                rad_const = RAD_CONST_E * particles.gamma^3
            elseif particles.mass == m_p
                rad_const = RAD_CONST_P * particles.gamma^3
            else
                rad_const = DTPSAD(0.0)
                println("SR is not implemented for this particle mass.")
            end
            StrMPoleSymplectic4RadPass!(r_in, ele.len, rad_const, 1.0, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit, ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, E0, num_particles, lost_flags)
        end
    else
        PolynomB[1] = ele.PolynomB[1]
        PolynomB[2] = ele.PolynomB[2] 
        PolynomB[3] = ele.PolynomB[3] / 2.0
        PolynomB[4] = ele.PolynomB[4] / 6.0
        if ele.rad == 0
            StrMPoleSymplectic4Pass!(r_in, ele.len, 1.0, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit, ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, num_particles, lost_flags)
        else
            if particles.mass == m_e
                rad_const = RAD_CONST_E * particles.gamma^3
            elseif particles.mass == m_p
                rad_const = RAD_CONST_P * particles.gamma^3
            else
                rad_const = DTPSAD(0.0)
                println("SR is not implemented for this particle mass.")
            end
            StrMPoleSymplectic4RadPass!(r_in, ele.len, rad_const, 1.0, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit, ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, E0, num_particles, lost_flags)
        end
    end
    return nothing
end

function pass!(ele::TKOCT, r_in::Vector{DTPSAD{N, T}}, num_particles::Int64, particles::TBeam) where {N, T <: Number}
    rad_const = DTPSAD(0.0)

    lost_flags = particles.lost_flag
    PolynomB = zeros(DTPSAD{N, T}, 4)
    E0 = particles.energy

    if ele.PolynomB[1] == 0.0 && ele.PolynomB[2] == 0.0 && ele.PolynomB[3] == 0.0 && ele.PolynomB[4] == 0.0
        PolynomB[4] = ele.k3 / 6.0
        if ele.rad == 0
            StrMPoleSymplectic4Pass!(r_in, ele.len, 1.0, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit, ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, num_particles, lost_flags)
        else
            if particles.mass == m_e
                rad_const = RAD_CONST_E * particles.gamma^3
            elseif particles.mass == m_p
                rad_const = RAD_CONST_P * particles.gamma^3
            else
                rad_const = DTPSAD(0.0)
                println("SR is not implemented for this particle mass.")
            end
            StrMPoleSymplectic4RadPass!(r_in, ele.len, 1.0, particles.mass, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit, ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, E0, num_particles, lost_flags)
        end
    else
        PolynomB[1] = ele.PolynomB[1]
        PolynomB[2] = ele.PolynomB[2] 
        PolynomB[3] = ele.PolynomB[3] / 2.0
        PolynomB[4] = ele.PolynomB[4] / 6.0
        if ele.rad == 0
            StrMPoleSymplectic4Pass!(r_in, ele.len, 1.0, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit, ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, num_particles, lost_flags)
        else
            if particles.mass == m_e
                rad_const = RAD_CONST_E * particles.gamma^3
            elseif particles.mass == m_p
                rad_const = RAD_CONST_P * particles.gamma^3
            else
                rad_const = DTPSAD(0.0)
                println("SR is not implemented for this particle mass.")
            end
            StrMPoleSymplectic4RadPass!(r_in, ele.len, rad_const, 1.0, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
                ele.FringeQuadEntrance, ele.FringeQuadExit, ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, E0, num_particles, lost_flags)
        end
    end
    return nothing
end

function BendSymplecticPass!(r::Vector{DTPSAD{N, T}}, le::DTPSAD{N, T}, beti::Float64, irho::DTPSAD{N, T}, A::Array{DTPSAD{N, T},1}, B::Array{DTPSAD{N, T},1}, 
    max_order::Int, num_int_steps::Int, entrance_angle::DTPSAD{N, T}, exit_angle::DTPSAD{N, T}, FringeBendEntrance::Int, FringeBendExit::Int,
    fint1::DTPSAD{N, T}, fint2::DTPSAD{N, T}, gap::DTPSAD{N, T}, FringeQuadEntrance::Int, FringeQuadExit::Int,
    fringeIntM0::Array{DTPSAD{N, T},1}, fringeIntP0::Array{DTPSAD{N, T},1}, T1::Array{DTPSAD{N, T},1}, T2::Array{DTPSAD{N, T},1}, 
    R1::Array{DTPSAD{N, T},2}, R2::Array{DTPSAD{N, T},2}, RApertures::Array{Float64,1}, EApertures::Array{Float64,1},
    KickAngle::Array{DTPSAD{N, T},1}, num_particles::Int, lost_flags::Array{Int64,1}) where {N, T <: Number}
    # Modified based on AT function. Ref[Terebilo, Andrei. "Accelerator modeling with MATLAB accelerator toolbox." PACS2001 (2001)].

    DRIFT1 = 0.6756035959798286638
    DRIFT2 = -0.1756035959798286639
    KICK1 = 1.351207191959657328
    KICK2 = -1.702414383919314656
    SL = le / num_int_steps
    L1 = SL * DRIFT1
    L2 = SL * DRIFT2
    K1 = SL * KICK1
    K2 = SL * KICK2

    if FringeQuadEntrance==2# && !isnothing(fringeIntM0) && !isnothing(fringeIntP0)
        useLinFrEleEntrance = 1
    else
        useLinFrEleEntrance = 0
    end
    if FringeQuadExit==2 #&& !isnothing(fringeIntM0) && !isnothing(fringeIntP0)
        useLinFrEleExit = 1
    else
        useLinFrEleExit = 0
    end

    B[1] -= sin(KickAngle[1]) / le
    A[1] += sin(KickAngle[2]) / le


    @inbounds for c in 1:num_particles
        if isone(lost_flags[c])
            continue
        end
        r6 = @view r[(c-1)*6+1:c*6]
        # Misalignment at entrance
        if !all(iszero, T1)
            addvv!(r6, T1)
        end
        if !all(iszero, R1)
            multmv!(r6, R1)
        end

        # Edge focus at entrance
        edge_fringe_entrance!(r6, irho, entrance_angle, fint1, gap, FringeBendEntrance)

        # Quadrupole gradient fringe entrance
        if !iszero(FringeQuadEntrance) && !iszero(B[2])
            if isone(useLinFrEleEntrance)
                linearQuadFringeElegantEntrance!(r6, B[2], fringeIntM0, fringeIntP0)
            else
                QuadFringePassP!(r6, B[2])
            end
        end

        # Integrator
        for m in 1:num_int_steps
            drift6!(r6, L1)
            bndthinkick!(r6, A, B, K1, irho, max_order, beti)
            drift6!(r6, L2)
            bndthinkick!(r6, A, B, K2, irho, max_order, beti)
            drift6!(r6, L2)
            bndthinkick!(r6, A, B, K1, irho, max_order, beti)
            drift6!(r6, L1)
        end

        # Quadrupole gradient fringe exit
        if !iszero(FringeQuadExit) && !iszero(B[2])
            if isone(useLinFrEleExit)
                linearQuadFringeElegantExit!(r6, B[2], fringeIntM0, fringeIntP0)
            else
                QuadFringePassN!(r6, B[2])
            end
        end

        # Edge focus at exit
        edge_fringe_exit!(r6, irho, exit_angle, fint2, gap, FringeBendExit)

        # Misalignment at exit
        if !all(iszero, R2)
            multmv!(r6, R2)
        end
        if !all(iszero, T2)
            addvv!(r6, T2)
        end
        if check_lost_GTPSA(r6)
            lost_flags[c] = 1
        end
    end

    B[1] += sin(KickAngle[1]) / le
    A[1] -= sin(KickAngle[2]) / le
    return nothing
end

function BendSymplecticPassRad!(r::Vector{DTPSAD{N, T}}, le::DTPSAD{N, T}, beti::Float64, irho::DTPSAD{N, T}, A::Array{DTPSAD{N, T},1}, B::Array{DTPSAD{N, T},1}, 
    max_order::Int, num_int_steps::Int, entrance_angle::DTPSAD{N, T}, exit_angle::DTPSAD{N, T}, FringeBendEntrance::Int, FringeBendExit::Int,
    fint1::DTPSAD{N, T}, fint2::DTPSAD{N, T}, gap::DTPSAD{N, T}, FringeQuadEntrance::Int, FringeQuadExit::Int,
    fringeIntM0::Array{DTPSAD{N, T},1}, fringeIntP0::Array{DTPSAD{N, T},1}, T1::Array{DTPSAD{N, T},1}, T2::Array{DTPSAD{N, T},1}, 
    R1::Array{DTPSAD{N, T},2}, R2::Array{DTPSAD{N, T},2}, RApertures::Array{Float64,1}, EApertures::Array{Float64,1},
    KickAngle::Array{DTPSAD{N, T},1}, E0::DTPSAD{N, T}, num_particles::Int, lost_flags::Array{Int64,1}) where {N, T <: Number}
    # AT function. Ref[Terebilo, Andrei. "Accelerator modeling with MATLAB accelerator toolbox." PACS2001 (2001)].

    DRIFT1 = 0.6756035959798286638
    DRIFT2 = -0.1756035959798286639
    KICK1 = 1.351207191959657328
    KICK2 = -1.702414383919314656
    SL = le / num_int_steps
    L1 = SL * DRIFT1
    L2 = SL * DRIFT2
    K1 = SL * KICK1
    K2 = SL * KICK2

    if FringeQuadEntrance==2# && !isnothing(fringeIntM0) && !isnothing(fringeIntP0)
        useLinFrEleEntrance = 1
    else
        useLinFrEleEntrance = 0
    end
    if FringeQuadExit==2# && !isnothing(fringeIntM0) && !isnothing(fringeIntP0)
        useLinFrEleExit = 1
    else
        useLinFrEleExit = 0
    end

    B[1] -= sin(KickAngle[1]) / le
    A[1] += sin(KickAngle[2]) / le


    @inbounds for c in 1:num_particles
        if isone(lost_flags[c])
            continue
        end
        r6 = @view r[(c-1)*6+1:c*6]
        # Misalignment at entrance
        if !all(iszero, T1)
            addvv!(r6, T1)
        end
        if !all(iszero, R1)
            multmv!(r6, R1)
        end
        # Edge focus at entrance
        edge_fringe_entrance!(r6, irho, entrance_angle, fint1, gap, FringeBendEntrance)

        # Quadrupole gradient fringe entrance
        if !iszero(FringeQuadEntrance) && !iszero(B[2])
            if useLinFrEleEntrance == 1
                linearQuadFringeElegantEntrance!(r6, B[2], fringeIntM0, fringeIntP0)
            else
                QuadFringePassP!(r6, B[2])
            end
        end

        # Integrator
        for m in 1:num_int_steps
            drift6!(r6, L1)
            bndthinkickrad!(r6, A, B, K1, irho, E0, max_order, beti)
            drift6!(r6, L2)
            bndthinkickrad!(r6, A, B, K2, irho, E0, max_order, beti)
            drift6!(r6, L2)
            bndthinkickrad!(r6, A, B, K1, irho, E0, max_order, beti)
            drift6!(r6, L1)
        end

        # Quadrupole gradient fringe exit
        if !iszero(FringeQuadExit) && !iszero(B[2])
            if useLinFrEleExit == 1
                linearQuadFringeElegantExit!(r6, B[2], fringeIntM0, fringeIntP0)
            else
                QuadFringePassN!(r6, B[2])
            end
        end

        # Edge focus at exit
        edge_fringe_exit!(r6, irho, exit_angle, fint2, gap, FringeBendExit)

        # Misalignment at exit
        if !all(iszero, R2)
            multmv!(r6, R2)
        end
        if !all(iszero, T2)
            addvv!(r6, T2)
        end
        if check_lost_GTPSA(r6)
            lost_flags[c] = 1
        end
    end

    
    B[1] += sin(KickAngle[1]) / le
    A[1] -= sin(KickAngle[2]) / le
    return nothing
end

function pass!(ele::TSBEND, r_in::Vector{DTPSAD{N, T}}, num_particles::Int64, particles::TBeam) where {N, T <: Number}
    lost_flags = particles.lost_flag
    irho = ele.angle / ele.len
    E0 = particles.energy
    if ele.rad == 0
        BendSymplecticPass!(r_in, ele.len, 1.0, irho, ele.PolynomA, ele.PolynomB, ele.MaxOrder, ele.NumIntSteps,
            ele.e1, ele.e2,
            ele.FringeBendEntrance, ele.FringeBendExit,
            ele.fint1, ele.fint2, ele.gap,
            ele.FringeQuadEntrance, ele.FringeQuadExit,
            ele.FringeIntM0, ele.FringeIntP0,
            ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures,
            ele.KickAngle, num_particles, lost_flags)
    else
        BendSymplecticPassRad!(r_in, ele.len, 1.0, irho, ele.PolynomA, ele.PolynomB, ele.MaxOrder, ele.NumIntSteps,
            ele.e1, ele.e2,
            ele.FringeBendEntrance, ele.FringeBendExit,
            ele.fint1, ele.fint2, ele.gap,
            ele.FringeQuadEntrance, ele.FringeQuadExit,
            ele.FringeIntM0, ele.FringeIntP0,
            ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures,
            ele.KickAngle, E0, num_particles, lost_flags)
    end
    return nothing
end

function pass!(ele::TESBEND, r_in::Array{DTPSAD{N, T},1}, num_particles::Int64, particles::TBeam) where {N, T <: Number}
    lost_flags = particles.lost_flag
    E0 = particles.energy
    rad_const = DTPSAD(0.0)
    if ele.rad == 0
        ExactSectorBend!(r_in, ele.len, 1.0, ele.angle, ele.PolynomA, ele.PolynomB, ele.MaxOrder, ele.NumIntSteps,
            ele.e1, ele.e2,
            ele.FringeBendEntrance, ele.FringeBendExit,
            ele.FringeQuadEntrance, ele.FringeQuadExit,
            ele.gK,
            ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures,
            ele.KickAngle, num_particles, lost_flags)
    else
        if particles.mass == m_e
            rad_const = RAD_CONST_E * particles.gamma^3
        elseif particles.mass == m_p
            rad_const = RAD_CONST_P * particles.gamma^3
        else
            rad_const = DTPSAD(0.0)
            println("SR is not implemented for this particle mass.")
        end
        ExactSectorBend_rad!(r_in, ele.len, rad_const, 1.0, ele.angle, ele.PolynomA, ele.PolynomB, ele.MaxOrder, ele.NumIntSteps,
            ele.e1, ele.e2,
            ele.FringeBendEntrance, ele.FringeBendExit,
            ele.FringeQuadEntrance, ele.FringeQuadExit,
            ele.gK,
            ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures,
            ele.KickAngle, num_particles, lost_flags)
    end
    return nothing
end

function pass!(ele::TRFCA, rin::Vector{DTPSAD{N, T}}, npart::Int64, particles::TBeam) where {N, T <: Number}
    if ele.energy == 0
        println("Energy is not defined for RFCA ", ele.name)
    end
    nv = ele.volt / ele.energy
    if ele.len == 0
        @inbounds for c in 1:npart
            if particles.lost_flag[c] == 1
                continue
            end
            r6 = @view rin[(c-1)*6+1:c*6]
            r6[6] += -nv * sin(2 * pi * ele.freq * ((r6[5] - ele.lag) / speed_of_light - 
                (ele.h / ele.freq - particles.T0) * 0.0) - ele.philag) / particles.beta^2
            if check_lost_GTPSA(r6)
                particles.lost_flag[c] = 1
            end
        end
    else
        halflength = ele.len / 2.0
        @inbounds for c in 1:npart
            if particles.lost_flag[c] == 1
                continue
            end
            r6 = @view rin[(c-1)*6+1:c*6]
            drift6!(r6, halflength)
            r6[6] += -nv * sin(2 * pi * ele.freq * ((r6[5] - ele.lag) / speed_of_light - 
                (ele.h / ele.freq - particles.T0) * 0.0) - ele.philag) / particles.beta^2
            drift6!(r6, halflength)
            if check_lost_GTPSA(r6)
                particles.lost_flag[c] = 1
            end
        end
    end
    return nothing
end

function pass!(cavity::TCRABCAVITY, rin::Vector{DTPSAD{N, T}}, npart::Int64, particles::TBeam) where {N, T <: Number}
    beta = particles.beta
    @inbounds for c in 1:npart
        if isone(particles.lost_flag[c])
            continue
        end
        r6 = @view rin[(c-1)*6+1:c*6]
        ang = cavity.k * r6[5] + cavity.phi
        if cavity.len == 0.0
            # r6[2] += (cavity.volt/beta2E) * sin(ang)
            # r6[6] += (-cavity.k * cavity.volt/beta2E) * r6[1] * cos(ang)
            # this is symplectic form
            r6[2] += (cavity.volt/particles.energy) * sin(ang/beta)
            r6[6] += (-cavity.k * cavity.volt/particles.energy/beta) * r6[1] * cos(ang/beta)
        else
            drift6!(r6, cavity.len / 2.0)
            # r6[2] += (cavity.volt/beta2E) * sin(ang)
            # r6[6] += (-cavity.k * cavity.volt/beta2E) * r6[1] * cos(ang)
            # this is symplectic form
            r6[2] += (cavity.volt/particles.energy) * sin(ang/beta)
            r6[6] += (-cavity.k * cavity.volt/particles.energy/beta) * r6[1] * cos(ang/beta)
            drift6!(r6, cavity.len / 2.0)
        end
    end
    return nothing
end

function pass!(cavity::TCRABCAVITYK2, rin::Vector{DTPSAD{N, T}}, npart::Int64, particles::TBeam) where {N, T <: Number}
    beta = particles.beta
    @inbounds for c in 1:npart
        if isone(particles.lost_flag[c])
            continue
        end
        r6 = @view rin[(c-1)*6+1:c*6]
        ang = cavity.k * r6[5] + cavity.phi
        if cavity.len == 0.0
            # r6[2] += (cavity.volt/beta2E) * sin(ang)
            # r6[6] += (-cavity.k * cavity.volt/beta2E) * r6[1] * cos(ang)
            # this is symplectic form
            r6[2] += (cavity.volt/particles.energy) * sin(ang/beta)
            r6[6] += (-cavity.k * cavity.volt/particles.energy/beta) * r6[1] * cos(ang/beta)
            # sextupole kick
            r6[2] -= cavity.k2 * (r6[1]^2 - r6[3]^2) * sin(ang)
            r6[4] += 2.0 * cavity.k2 * r6[1] * r6[3] * sin(ang)
            r6[6] -= (cavity.k2 * cavity.k / 3.0) * (r6[1]^3 - 3.0 * r6[1] * r6[3]^2) * cos(ang)
        else
            drift6!(r6, cavity.len / 2.0)
            # r6[2] += (cavity.volt/beta2E) * sin(ang)
            # r6[6] += (-cavity.k * cavity.volt/beta2E) * r6[1] * cos(ang)            
            # this is symplectic form
            r6[2] += (cavity.volt/particles.energy) * sin(ang/beta)
            r6[6] += (-cavity.k * cavity.volt/particles.energy/beta) * r6[1] * cos(ang/beta)
            # sextupole kick
            r6[2] -= cavity.k2 * (r6[1]^2 - r6[3]^2) * sin(ang)
            r6[4] += 2.0 * cavity.k2 * r6[1] * r6[3] * sin(ang)
            r6[6] -= (cavity.k2 * cavity.k / 3.0) * (r6[1]^3 - 3.0 * r6[1] * r6[3]^2) * cos(ang)
            drift6!(r6, cavity.len / 2.0)
        end
    end
    return nothing
end

function pass!(ele::TthinMULTIPOLE, rin::Vector{DTPSAD{N, T}}, npart::Int64, particles::TBeam) where {N, T <: Number}
    # Modified based on AT function. Ref[Terebilo, Andrei. "Accelerator modeling with MATLAB accelerator toolbox." PACS2001 (2001)].

    # no bending
    bax = 0.0 
    bay = 0.0

    A = [ele.PolynomA[i] for i in 1:4]
    B = [ele.PolynomB[i] for i in 1:4]
    B[1] -= ele.KickAngle[1]
    A[1] += ele.KickAngle[2]

    # for c in 1:npart
    @inbounds for c in 1:npart
        if particles.lost_flag[c] == 1
            continue
        end
        r6 = @view rin[(c-1)*6+1:c*6]
        # Misalignment at entrance
        if !all(iszero, ele.T1)
            addvv!(r6, ele.T1)
        end
        if !all(iszero, ele.R1)
            multmv!(r6, ele.R1)
        end

        strthinkick1!(r6, A, B, DTPSAD(1.0), ele.MaxOrder)
        r6[2] += bax * r6[6]
        r6[4] -= bay * r6[6]
        r6[6] -= bax * r6[1] - bay * r6[3]  # Path lenghtening

        # Misalignment at exit
        if !all(iszero, ele.R2)
            multmv!(r6, ele.R2)
        end
        if !all(iszero, ele.T2)
            addvv!(r6, ele.T2)
        end
        if check_lost(r6)
            particles.lost_flag[c] = 1
        end
    end

    B[1] += ele.KickAngle[1]
    A[1] -= ele.KickAngle[2]
    return nothing
end

function pass!(elem::TTRANSLATION, r_in::Vector{DTPSAD{N, T}}, num_particles::Int64, particles::TBeam) where {N, T <: Number}
    @inbounds for c in 1:num_particles
        if isone(particles.lost_flag[c])
            continue
        end
        r6 = @view r_in[(c-1)*6+1:c*6]
        if use_exact_beti == 1
            pz = sqrt(1.0 + 2.0 * r6[6] / particles.beta + r6[6]^2 - r6[2]^2 - r6[4]^2)
            r6[1] -= elem.dx + elem.ds * r6[2] / pz
            r6[3] -= elem.dy + elem.ds * r6[4] / pz
            r6[5] += elem.ds * (1.0/particles.beta + r6[6]) / pz
        else
            pz = sqrt(1.0 + 2.0 * r6[6] + r6[6]^2 - r6[2]^2 - r6[4]^2)
            r6[1] -= elem.dx + elem.ds * r6[2] / pz
            r6[3] -= elem.dy + elem.ds * r6[4] / pz
            r6[5] += elem.ds * (1.0 + r6[6]) / pz
        end
        if check_lost_GTPSA(r6)
            particles.lost_flag[c] = 1
        end
    end
    return nothing
end

function pass!(elem::TYROTATION, r_in::Vector{DTPSAD{N, T}}, num_particles::Int64, particles::TBeam) where {N, T <: Number}
    angle = -elem.angle
    if angle == 0.0
        return nothing
    end
    ca = cos(angle)
    sa = sin(angle)
    ta = tan(angle)
    if use_exact_beti == 1
        beta = particles.beta
    else
        beta = 1.0
    end
    @inbounds for c in 1:num_particles
        if isone(particles.lost_flag[c])
            continue
        end
        r6 = @view r_in[(c-1)*6+1:c*6]
        x, px, y, py, t, pt = r6[1], r6[2], r6[3], r6[4], r6[5], r6[6]
        pz = sqrt(1.0 + 2.0 * pt / beta + pt^2 - px^2 - py^2)
        ptt = 1.0 - ta*px/pz
        r6[1] = x/(ca*ptt)
        r6[2] = ca*px + sa*pz
        r6[3] = y + ta*x*py/(pz*ptt)
        r6[5] = t + ta*x*(1.0 / beta+pt)/(pz*ptt)

        if check_lost_GTPSA(r6)
            particles.lost_flag[c] = 1
        end
    end
    return nothing
end

function linepass!(line::Vector{<:AbstractElement}, particles::TBeam)
    np = particles.nmacro
    particles6 = matrix_to_array(particles.r)
    if length(particles6) != np * 6
        error("The number of particles does not match the length of the particle array")
    end
    for ele in line
        pass!(ele, particles6, np, particles)
        # pass!(ele, particles.r, np, particles)
    end
    rout = array_to_matrix(particles6, np)
    particles.r = rout
    return nothing
end

function linepass!(line::Vector{<:AbstractElement}, particles::TBeam, refpts::Vector)
    # Note!!! A lost particle's coordinate will not be marked as NaN or Inf like other softwares 
    # Check if the particle is lost by checking the lost_flag
    np = particles.nmacro
    particles6 = matrix_to_array(particles.r)
    if length(particles6) != np * 6
        error("The number of particles does not match the length of the particle array")
    end
    saved_particles = []
    for i in eachindex(line)
        # ele = line[i]
        pass!(line[i], particles6, np, particles)        
        if i in refpts
            push!(saved_particles, copy(array_to_matrix(particles6, np)))
        end
    end
    rout = array_to_matrix(particles6, np)
    particles.r = rout
    return saved_particles
end

function ringpass!(line::Vector{<:AbstractElement}, particles::TBeam, nturns::Int)
    for turn in 1:nturns
        linepass!(line, particles)
    end
    return nothing
end
function matrix_to_array(r::Matrix{DTPSAD{N, T}}) where {N, T <: Number}
    np = size(r, 1)
    if size(r, 2) != 6
        error("Input matrix must have 6 columns")
    end
    arr = zeros(DTPSAD{N, T}, np * 6)
    for i in 1:np
        for j in 1:6
            arr[(i-1)*6 + j] = r[i, j]
        end
    end
    return arr
end
function array_to_matrix(r::Vector{DTPSAD{N, T}}, np::Int) where {N, T <: Number}
    if length(r) != np * 6
        error("Input vector length must match the number of particles times 6")
    end
    mat = zeros(DTPSAD{N, T}, np, 6)
    for i in 1:np
        for j in 1:6
            mat[i, j] = r[(i-1)*6 + j]
        end
    end
    return mat
end

function _strip(param)
    if isa(param, String)
        return param
    elseif isa(param, Float64)
        return param
    elseif isa(param, Int64) || isa(param, Int32)
        return param
    elseif isa(param, DTPSAD{NVAR(), Float64}) 
        return param.val
    elseif isa(param, Vector{DTPSAD{NVAR(), Float64}}) 
        return [_strip(p) for p in param]
    elseif isa(param, Matrix{DTPSAD{NVAR(), Float64}})
        Mat = zeros(Float64, size(param))
        for i in eachindex(param)
            Mat[i] = _strip(param[i])
        end
        return Mat
    elseif isa(param, Vector{Float64}) || isa(param, Vector{Int64}) || isa(param, Vector{Int32})
        return param
    elseif isa(param, Matrix{Float64}) || isa(param, Matrix{Int64}) || isa(param, Matrix{Int32})
        return Param
    end
end
function to_Number(te::AbstractTPSAElement)
    real_sym = Symbol(startswith(string(nameof(typeof(te))), "T") ? 
                string(nameof(typeof(te)))[2:end] : 
                error("Type $(typeof(te)) does not start with T"))
    real_type = getfield(parentmodule(typeof(te)), real_sym)

    vals = map(f -> _strip(getfield(te, f)), fieldnames(typeof(te)))

    # 2-c. call inner constructor directly
    return real_type(vals...)
end

# convert real elements to their corresponding types
function TPSAD2Number(line::Vector{<:AbstractTPSAElement})
    return [to_Number(ele) for ele in line]
end

## reverse the conversion
function _strip_inverse(param, fieldname::Symbol)
    if isa(param, String)
        return param
    elseif isa(param, Float64)
        return DTPSAD(param)
    elseif isa(param, Int64) || isa(param, Int32) || 
        isa(param, Vector{Int64}) || isa(param, Vector{Int32}) || 
        isa(param, Matrix{Int64}) || isa(param, Matrix{Int32})
        return param
    elseif isa(param, Vector{Float64}) 
        if fieldname == :RApertures || fieldname == :EApertures
            return param
        else
            return [DTPSAD(p) for p in param]
        end
    elseif isa(param, Matrix{Float64}) 
        Mat = zeros(DTPSAD{NVAR(), Float64}, size(param, 1), size(param, 2))
        for i in eachindex(param)
            Mat[i] = DTPSAD(param[i])
        end
        return Mat
    end
end
function to_TPSAD(ele::AbstractElement)
    sym = Symbol("T" * string(nameof(typeof(ele))))
    vals = map(f -> _strip_inverse(getfield(ele, f), f), fieldnames(typeof(ele)))
    TPSAD_type = getfield(parentmodule(typeof(ele)), sym)
    return TPSAD_type(vals...)
end
function Number2TPSAD(line::Vector{<:AbstractElement})
    return [to_TPSAD(ele) for ele in line]
end