abstract type AbstractElement end
struct Drift <: AbstractElement
    name::String
    len::Float64
end

struct Quad <: AbstractElement
    name::String
    len::Float64
    k1::Float64 
    k1s::Float64
end

struct ThinQuad <: AbstractElement
    name::String
    len::Float64
    k1l::Float64 
    k1sl::Float64
    function ThinQuad(name::String, k1l::Float64, k1sl::Float64)
        return new(name, 0.0, k1l, k1sl)
    end
end

struct SBend <: AbstractElement
    name::String
    len::Float64
    angle::Float64
    k1::Float64
    e1::Float64
    e2::Float64
    f_int1::Float64
    f_int2::Float64
end

struct RBend <: AbstractElement
    name::String
    len::Float64
    angle::Float64
    k1::Float64
    e1::Float64
    e2::Float64
    f_int1::Float64
    f_int2::Float64
end

struct Bend <: AbstractElement
    name::String
    len::Float64
    angle::Float64
    e1::Float64
    e2::Float64
    grad::Float64
end

struct DipEdge <: AbstractElement
    name::String
    len::Float64
    h::Float64
    e1::Float64
    f_int::Float64
    function DipEdge(name::String, h::Float64, e1::Float64, f_int::Float64)
        return new(name, 0.0, h, e1, f_int)
    end
end

struct DipBody <: AbstractElement
    name::String
    len::Float64
    angle::Float64
    k1::Float64
end

struct ThinCrabCavity <: AbstractElement
    name::String
    len::Float64
    strength::Tuple{Float64,Float64}
    frequency::Float64
    phase::Float64
    kcc::Float64
end
function ThinCrabCavity(;name::String="",strength::Tuple{Float64,Float64},frequency::Float64,phase::Float64=0.0)
	kcc::Float64=2pi*Frequency/Float64(299792458)
	ThinCrabCavity(name,0,strength,frequency,phase,kcc)
end

struct Marker <: AbstractElement
    name::String
    len::Float64
    function Marker(name::String)
        return new(name, 0.0)
    end
end

struct Solenoid <: AbstractElement
    name::String
    len::Float64
    ks::Float64
end

struct LorentzBoost <: AbstractElement
    name::String
    len::Float64
    thetac::Float64
    thetas::Float64
    function LorentzBoost(name::String,thetac::Float64,thetas::Float64)
        return new(name, 0.0, thetac, thetas)
    end
end

struct RevLorentzBoost <: AbstractElement
    name::String
    len::Float64
    thetac::Float64
    thetas::Float64
    function RevLorentzBoost(name::String,thetac::Float64,thetas::Float64)
        return new(name, 0.0, thetac, thetas)
    end
end