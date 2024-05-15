include("optics.jl")
using Base: @kwdef
abstract type AbstractElement end

@kwdef struct MARKER <: AbstractElement
    name::String = "MARKER"
    len::Float64 = 0.0
    eletype::String = "MARKER"
end

@kwdef struct DRIFT <: AbstractElement
    name::String = "EDrift"
    len::Float64 = 0.0
    T1::Array{Float64,1} = zeros(6)
    T2::Array{Float64,1} = zeros(6)
    R1::Array{Float64,2} = zeros(6,6)
    R2::Array{Float64,2} = zeros(6,6)        
    RApertures::Array{Float64,1} = zeros(6)
    EApertures::Array{Float64,1} = zeros(6)
    eletype::String = "DRIFT"
end


@kwdef struct KQUAD <: AbstractElement
    name::String  = "Quad"                                      # element name  
    len::Float64 = 0.0
    k1::Float64 = 0.0                                           # use k1 if PolynomB is not given
    PolynomA::Array{Float64,1} = zeros(4)    
    PolynomB::Array{Float64,1} = zeros(4)          # PolynomB has higher priority than k1, k2, k3
    MaxOrder::Int64 = 1
    NumIntSteps::Int64 = 10
    rad::Int64 = 0
    FringeQuadEntrance::Int64 = 0
    FringeQuadExit::Int64 = 0
    FringeIntM0::Array{Float64,1} = zeros(5)
    FringeIntP0::Array{Float64,1} = zeros(5)
    T1::Array{Float64,1} = zeros(6)
    T2::Array{Float64,1} = zeros(6)
    R1::Array{Float64,2} = zeros(6,6)
    R2::Array{Float64,2} = zeros(6,6)         
    RApertures::Array{Float64,1} = zeros(6)
    EApertures::Array{Float64,1} = zeros(6)
    KickAngle::Array{Float64,1} = zeros(2)
    eletype::String = "KQUAD"
end

@kwdef struct KSEXT <: AbstractElement
    name::String  = "Sext"                                      # element name  
    len::Float64 = 0.0
    k2::Float64 = 0.0                                           # use k2 if PolynomB is not given
    PolynomA::Array{Float64,1} = zeros(4)    
    PolynomB::Array{Float64,1} = zeros(4)           # PolynomB has higher priority than k1, k2, k3
    MaxOrder::Int64 = 2
    NumIntSteps::Int64 = 10
    rad::Int64 = 0
    FringeQuadEntrance::Int64 = 0
    FringeQuadExit::Int64 = 0
    FringeIntM0::Array{Float64,1} = zeros(5)
    FringeIntP0::Array{Float64,1} = zeros(5)
    T1::Array{Float64,1} = zeros(6)
    T2::Array{Float64,1} = zeros(6)
    R1::Array{Float64,2} = zeros(6, 6)
    R2::Array{Float64,2} = zeros(6, 6)         
    RApertures::Array{Float64,1} = zeros(6)
    EApertures::Array{Float64,1} = zeros(6)
    KickAngle::Array{Float64,1} = zeros(2)
    eletype::String = "KSEXT"
end

@kwdef struct KOCT <: AbstractElement
    name::String  = "OCT"                                       # element name  
    len::Float64 = 0.0
    k3::Float64 = 0.0                                           # use k3 if PolynomB is not given
    PolynomA::Array{Float64,1} = zeros(4)    
    PolynomB::Array{Float64,1} = zeros(4)           # PolynomB has higher priority than k1, k2, k3
    MaxOrder::Int64 = 3
    NumIntSteps::Int64 = 10
    rad::Int64 = 0
    FringeQuadEntrance::Int64 = 0
    FringeQuadExit::Int64 = 0
    FringeIntM0::Array{Float64,1} = zeros(5)
    FringeIntP0::Array{Float64,1} = zeros(5)
    T1::Array{Float64,1} = zeros(6)
    T2::Array{Float64,1} = zeros(6)
    R1::Array{Float64,2} = zeros(6, 6)
    R2::Array{Float64,2} = zeros(6, 6)        
    RApertures::Array{Float64,1} = zeros(6)
    EApertures::Array{Float64,1} = zeros(6)
    KickAngle::Array{Float64,1} = zeros(2)
    eletype::String = "KOCT"
end

@kwdef struct thinMULTIPOLE <: AbstractElement
    name::String  = "thinMULTIPOLE"                                       # element name  
    len::Float64 = 0.0
    PolynomA::Array{Float64,1} = zeros(4)    
    PolynomB::Array{Float64,1} = zeros(4)           # PolynomB has higher priority than k1, k2, k3
    MaxOrder::Int64 = 1
    NumIntSteps::Int64 = 1
    rad::Int64 = 0
    FringeQuadEntrance::Int64 = 0
    FringeQuadExit::Int64 = 0
    FringeIntM0::Array{Float64,1} = zeros(5)
    FringeIntP0::Array{Float64,1} = zeros(5)
    T1::Array{Float64,1} = zeros(6)
    T2::Array{Float64,1} = zeros(6)
    R1::Array{Float64,2} = zeros(6, 6)
    R2::Array{Float64,2} = zeros(6, 6)        
    RApertures::Array{Float64,1} = zeros(6)
    EApertures::Array{Float64,1} = zeros(6)
    KickAngle::Array{Float64,1} = zeros(2)
    eletype::String = "thinMULTIPOLE"
end

@kwdef struct SBEND <: AbstractElement
    name::String = "SBend"
    len::Float64 = 0.0
    angle::Float64 = 0.0
    e1::Float64 = 0.0
    e2::Float64 = 0.0
    PolynomA::Array{Float64,1} = zeros(4)    
    PolynomB::Array{Float64,1} = zeros(4)
    MaxOrder::Int64 = 0
    NumIntSteps::Int64 = 10
    rad::Int64 = 0
    fint1::Float64 = 0.0
    fint2::Float64 = 0.0
    gap::Float64 = 0.0
    FringeBendEntrance::Int64 = 1
    FringeBendExit::Int64 = 1
    FringeQuadEntrance::Int64 = 0
    FringeQuadExit::Int64 = 0
    FringeIntM0::Array{Float64,1} = zeros(5)
    FringeIntP0::Array{Float64,1} = zeros(5)
    T1::Array{Float64,1} = zeros(6)
    T2::Array{Float64,1} = zeros(6)
    R1::Array{Float64,2} = zeros(6,6)
    R2::Array{Float64,2} = zeros(6,6)         
    RApertures::Array{Float64,1} = zeros(6)
    EApertures::Array{Float64,1} = zeros(6)
    KickAngle::Array{Float64,1} = zeros(2)
    eletype::String = "SBEND"
end

struct RBEND <: AbstractElement
    name::String
    len::Float64 
    angle::Float64 
    e1::Float64 
    e2::Float64 
    PolynomA::Array{Float64,1}
    PolynomB::Array{Float64,1} 
    MaxOrder::Int64 
    NumIntSteps::Int64
    rad::Int64
    fint1::Float64
    fint2::Float64
    gap::Float64
    FringeBendEntrance::Int64
    FringeBendExit::Int64
    FringeQuadEntrance::Int64
    FringeQuadExit::Int64
    FringeIntM0::Array{Float64,1}
    FringeIntP0::Array{Float64,1}
    T1::Array{Float64,1}
    T2::Array{Float64,1}
    R1::Array{Float64,2}
    R2::Array{Float64,2}     
    RApertures::Array{Float64,1}
    EApertures::Array{Float64,1}
    KickAngle::Array{Float64,1}
    eletype::String
end
function RBEND(;name::String = "RBend", len::Float64 = 0.0, angle::Float64 = 0.0, PolynomA::Array{Float64,1} = zeros(4), 
                PolynomB::Array{Float64,1} = zeros(4), MaxOrder::Int64 = 0, NumIntSteps::Int64 = 10, rad::Int64=0, fint1::Float64 = 0.0, 
                fint2::Float64 = 0.0, gap::Float64 = 0.0, FringeBendEntrance::Int64 = 1, FringeBendExit::Int64 = 1, 
                FringeQuadEntrance::Int64 = 0, FringeQuadExit::Int64 = 0, FringeIntM0::Array{Float64,1} = zeros(5), 
                FringeIntP0::Array{Float64,1} = zeros(5), T1::Array{Float64,1} = zeros(6), T2::Array{Float64,1} = zeros(6), 
                R1::Array{Float64,2} = zeros(6,6), R2::Array{Float64,2} = zeros(6,6), RApertures::Array{Float64,1} = zeros(6), 
                EApertures::Array{Float64,1} = zeros(6), KickAngle::Array{Float64,1} = zeros(2))
    e1 = angle/2.0
    e2 = angle/2.0
    return RBEND(name, len, angle, e1, e2, PolynomA, PolynomB, MaxOrder, NumIntSteps, rad, fint1, fint2, gap, FringeBendEntrance, FringeBendExit, FringeQuadEntrance, FringeQuadExit, FringeIntM0, FringeIntP0, T1, T2, R1, R2, RApertures, EApertures, KickAngle, "RBEND")
end


@kwdef struct RFCA <: AbstractElement
    name::String = "RFCA"
    len::Float64 = 0.0
    volt::Float64 = 0.0
    freq::Float64 = 0.0
    h::Float64 = 1.0
    lag::Float64 = 0.0
    philag::Float64 = 0.0
    energy::Float64 = 0.0 # eV
    eletype::String = "RFCA"
end

@kwdef struct SOLENOID <: AbstractElement
    name::String = "Solenoid"
    len::Float64 = 0.0
    ks::Float64 = 0.0
    T1::Array{Float64,1} = zeros(6)
    T2::Array{Float64,1} = zeros(6)
    R1::Array{Float64,2} = zeros(6,6)
    R2::Array{Float64,2} = zeros(6,6)
    eletype::String = "SOLENOID"
end

@kwdef struct CORRECTOR <: AbstractElement
    name::String = "HKicker"
    len::Float64 = 0.0
    xkick::Float64 = 0.0
    ykick::Float64 = 0.0
    T1::Array{Float64,1} = zeros(6)
    T2::Array{Float64,1} = zeros(6)
    R1::Array{Float64,2} = zeros(6,6)
    R2::Array{Float64,2} = zeros(6,6)
    eletype::String = "CORRECTOR"
end
function HKICKER(;name="HKicker", len=0.0, xkick=0.0)
    return CORRECTOR(name, len, xkick, 0.0, zeros(6), zeros(6), zeros(6,6), zeros(6,6), "HKICKER")
end
function VKICKER(;name="VKicker", len=0.0, ykick=0.0)
    return CORRECTOR(name, len, 0.0, ykick, zeros(6), zeros(6), zeros(6,6), zeros(6,6), "VKICKER")
end

# non-canonical elements
@kwdef struct QUAD <: AbstractElement
    name::String  = "Quad"                                      # element name  
    len::Float64 = 0.0
    k1::Float64 = 0.0                                           # use k1 if PolynomB is not given
    rad::Int64 = 0
    T1::Array{Float64,1} = zeros(6)
    T2::Array{Float64,1} = zeros(6)
    R1::Array{Float64,2} = zeros(6,6)
    R2::Array{Float64,2} = zeros(6,6)         
    RApertures::Array{Float64,1} = zeros(6)
    EApertures::Array{Float64,1} = zeros(6)
    eletype::String = "QUAD"
end

###########################################
# the following elements may not be symplectic and may not work with Enzyme
struct CRABCAVITY <: AbstractElement
    name::String 
    len::Float64 
    volt::Float64  # voltage
    freq::Float64  # frequency
    k::Float64  # wave number
    phi::Float64  # phase
    errors::Array{Float64,1} # 1: Voltage error, 2: Phase error
    energy::Float64
    eletype::String 
end
function CRABCAVITY(;name::String = "CRABCAVITY", len::Float64 = 0.0, volt::Float64 = 0.0, freq::Float64 = 0.0, phi::Float64 = 0.0, errors::Array{Float64,1} = zeros(2), energy::Float64 = 1e9)
    k = 2*π*freq/2.99792458e8
    return CRABCAVITY(name, len, volt, freq, k, phi, errors, energy, "CRABCAVITY")
end

struct easyCRABCAVITY <: AbstractElement
    name::String 
    len::Float64 
    halfthetac::Float64 
    freq::Float64 
    k::Float64 
    phi::Float64 
    errors::Array{Float64,1}  # 1: Voltage error, 2: Phase error
    eletype::String
end
function easyCRABCAVITY(;name::String = "easyCRABCAVITY", len::Float64 = 0.0, halfthetac::Float64 = 0.0, freq::Float64 = 0.0, phi::Float64 = 0.0, errors::Array{Float64,1} = [0.0, 0.0])
    k = 2*π*freq/2.99792458e8
    return easyCRABCAVITY(name, len, halfthetac, freq, k, phi, errors, "easyCRABCAVITY")
end

@kwdef struct AccelCavity <: AbstractElement
    name::String = "AccelCavity"
    len::Float64 = 0.0
    volt::Float64 = 0.0 # voltage
    freq::Float64 = 0.0 # frequency
    k::Float64 = 0.0 # wave number
    h::Float64 = 1.0 # harmonic number
    phis::Float64 = 0.0 # synchronous phase π/2 for accelerating on crest
    eletype::String = "AccelCavity"
end
function AccelCavity(freq, nv, h, phis; name="AccelCavity", len=0.0)
    k = 2*π*freq/2.99792458e8
    return AccelCavity(name, len, nv, freq, k, h, phis, "AccelCavity")
end

abstract type AbstractTransferMap <:AbstractElement end
abstract type AbstractTransverseMap <:AbstractTransferMap end
abstract type AbstractLongitudinalMap <:AbstractTransferMap end
struct LongitudinalRFMap <: AbstractLongitudinalMap
    alphac::Float64
    RF::AbstractElement
    LongitudinalRFMap(alphac::Float64, RF::AbstractElement)=new(alphac, RF)
end

struct LorentzBoost <: AbstractElement
    angle::Float64
    cosang::Float64
    tanang::Float64
    mode::Int
    LorentzBoost(angle)=new(angle, cos(angle), tan(angle), 0)
end

struct InvLorentzBoost <: AbstractElement
    angle::Float64
    sinang::Float64
    cosang::Float64
    mode::Int
    InvLorentzBoost(angle)=new(angle, sin(angle), cos(angle), 0)
end


######### strong beam-beam
abstract type AbstractStrongBeamBeam <:AbstractElement end

struct StrongThinGaussianBeam <: AbstractStrongBeamBeam
    amplitude::Float64
    rmssizex::Float64
    rmssizey::Float64
    zloc::Float64
    xoffset::Float64
    yoffset::Float64
    StrongThinGaussianBeam(amp::Float64, rx::Float64, ry::Float64, zloc::Float64=0.0, xoff::Float64=0.0, yoff::Float64=0.0)=new(amp,rx,ry,zloc,xoff,yoff)
end

struct StrongGaussianBeam <: AbstractStrongBeamBeam  # Strong Beam with transverse Gaussian distribution
    # particle::ParticleType
    charge::Float64  
    mass::Float64  
    atomnum::Float64  
    classrad0::Float64  
    radconst::Float64  
    num_particle::Int  # Number of particles
    total_energy::Float64 # Total energy of the beam
    momentum::Float64  # Design Momentum of the beam
    gamma::Float64  # Relativistic gamma
    beta::Float64  # Relativistic beta v/c
    optics::AbstractOptics4D # optics @IP
    beamsize::Vector{Float64} # Beam size at IP
    nzslice::Int # Number of slices in z direction
    zslice_center::Vector{Float64} # z center of each slice
    zslice_npar::Vector{Float64} # amplitude of each slice
    xoffsets::Vector{Float64} # x offset of each slice
    yoffsets::Vector{Float64} # y offset of each slice
    function StrongGaussianBeam(charge::Float64, mass::Float64, atomnum::Float64, 
            np::Int, energy::Float64, op::AbstractOptics4D, bs::Vector{Float64}, nz::Int)
        momentum=sqrt(energy*energy-mass*mass)  
        gamma=energy/mass
        beta=momentum/energy
        classrad0=charge*charge/(atomnum*mass)/4/pi/55.26349406*1e-6
        radconst=4*pi/3*classrad0/mass/mass/mass
        new(charge,mass,atomnum,classrad0,radconst,np,energy,momentum,gamma,beta, op, bs, nz, zeros(nz), zeros(nz), zeros(nz), zeros(nz))
    end  
end
