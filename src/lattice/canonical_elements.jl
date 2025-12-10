"""
    AbstractElement{T}

Abstract type for all elements in the lattice.
"""
abstract type AbstractElement{T} end

"""
    MARKER(;name::String = "MARKER", len::Float64 = 0.0)

A marker element.
Example:
```julia
marker = MARKER(name="MARKER1")
```
"""
mutable struct MARKER{T} <: AbstractElement{T}
    name::String
    len::T
    eletype::String

    MARKER(name::String, len::T, eletype::String) where {T} = new{T}(name, len, eletype)
end
function MARKER(; name::String = "MARKER", len = 0.0)
    return MARKER(name, len, "MARKER")
end

"""
    DRIFT(;name::String = "DRIFT", len::Float64 = 0.0, T1::Array{Float64,1} = zeros(6), 
        T2::Array{Float64,1} = zeros(6), R1::Array{Float64,2} = zeros(6,6), R2::Array{Float64,2} = zeros(6,6), 
        RApertures::Array{Float64,1} = zeros(6), EApertures::Array{Float64,1} = zeros(6))

A drift element.

# Arguments
- name::String: element name
- len::Float64: element length
- T1::Array{Float64,1}: misalignment at entrance
- T2::Array{Float64,1}: misalignment at exit
- R1::Array{Float64,2}: rotation at entrance
- R2::Array{Float64,2}: rotation at exit
- RApertures::Array{Float64,1}: rectangular apertures. Not implemented yet.
- EApertures::Array{Float64,1}: elliptical apertures. Not implemented yet.
Example:
```julia
drift = DRIFT(name="D1", len=1.0)
```
"""
mutable struct DRIFT{T} <: AbstractElement{T}
    name::String
    len::T
    T1::Array{T,1}
    T2::Array{T,1}
    R1::Array{T,2}
    R2::Array{T,2}        
    RApertures::Array{Float64,1}
    EApertures::Array{Float64,1}
    eletype::String

    DRIFT(name::String, len::T, T1::Array{T,1}, T2::Array{T,1}, 
        R1::Array{T,2}, R2::Array{T,2}, RApertures::Array{Float64,1}, EApertures::Array{Float64,1}, eletype::String) where {T} = 
        new{T}(name, len, T1, T2, R1, R2, RApertures, EApertures, eletype)
end
function DRIFT(;name::String = "DRIFT", len = 0.0, T1 = zeros(6), 
        T2 = zeros(6), R1 = zeros(6,6), R2 = zeros(6,6), 
        RApertures = zeros(6), EApertures = zeros(6))
        if len isa DTPSAD || T1[1] isa DTPSAD || T2[1] isa DTPSAD || R1[1,1] isa DTPSAD || R2[1,1] isa DTPSAD
            return DRIFT(name, DTPSAD(len), DTPSAD.(T1), DTPSAD.(T2), 
                DTPSAD.(R1), DTPSAD.(R2), RApertures, EApertures, "DRIFT")
        end
    return DRIFT(name, len, T1, T2, R1, R2, RApertures, EApertures, "DRIFT")
end

"""
    DRIFT_SC(;name::String = "DRIFT_SC", len::Float64 = 0.0, T1::Array{Float64,1} = zeros(6), 
        T2::Array{Float64,1} = zeros(6), R1::Array{Float64,2} = zeros(6,6), R2::Array{Float64,2} = zeros(6,6), 
        RApertures::Array{Float64,1} = zeros(6), EApertures::Array{Float64,1} = zeros(6), a::Float64 = 1.0, b::Float64 = 1.0,
        Nl::Int64 = 10, Nm::Int64 = 10, Nsteps::Int64=1)

A drift element with space charge.

# Arguments
- name::String: element name
- len::Float64: element length
- T1::Array{Float64,1}: misalignment at entrance
- T2::Array{Float64,1}: misalignment at exit
- R1::Array{Float64,2}: rotation at entrance
- R2::Array{Float64,2}: rotation at exit
- RApertures::Array{Float64,1}: rectangular apertures. Not implemented yet.
- EApertures::Array{Float64,1}: elliptical apertures. Not implemented yet.
- a::Float64: horizontal size of the perfectly conducting pipe
- b::Float64: vertical size of the perfectly conducting pipe
- Nl::Int64: number of mode in the horizontal direction
- Nm::Int64: number of mode in the vertical direction
- Nsteps::Int64: number of steps for space charge calculation. One step represents a half-kick-half.
Example:
```julia
drift = DRIFT_SC(name="D1_SC", len=0.5, a=13e-3, b=13e-3, Nl=15, Nm=15)
```
"""
mutable struct DRIFT_SC{T} <: AbstractElement{T}
    name::String
    len::T
    T1::Array{T,1}
    T2::Array{T,1}
    R1::Array{T,2}
    R2::Array{T,2}        
    RApertures::Array{Float64,1}
    EApertures::Array{Float64,1}
    a::T
    b::T
    Nl::Int64
    Nm::Int64
    Nsteps::Int64 # Number of steps for space charge calculation. One step represents a half-kick-half.
    eletype::String

    # constructor with all parameters
    DRIFT_SC(name::String, len::T, T1::Array{T,1}, T2::Array{T,1}, 
        R1::Array{T,2}, R2::Array{T,2}, RApertures::Array{Float64,1}, EApertures::Array{Float64,1}, 
        a::T, b::T, Nl::Int64, Nm::Int64, Nsteps::Int64, eletype::String) where T = 
        new{T}(name, len, T1, T2, R1, R2, RApertures, EApertures, a, b, Nl, Nm, Nsteps, eletype)
end
function DRIFT_SC(;name::String = "DRIFT_SC", len = 0.0, T1 = zeros(6), 
    T2 = zeros(6), R1 = zeros(6,6), R2 = zeros(6,6), 
    RApertures = zeros(6), EApertures = zeros(6), a = 1.0, b = 1.0,
    Nl = 10, Nm = 10, Nsteps = 1) 
    if len isa DTPSAD || T1[1] isa DTPSAD || T2[1] isa DTPSAD || R1[1,1] isa DTPSAD || R2[1,1] isa DTPSAD ||
        a isa DTPSAD || b isa DTPSAD
        return DRIFT_SC(name, DTPSAD(len), DTPSAD.(T1), DTPSAD.(T2), 
            DTPSAD.(R1), DTPSAD.(R2), RApertures, EApertures, 
            DTPSAD(a), DTPSAD(b), Nl, Nm, Nsteps, "DRIFT_SC")
    end
    return DRIFT_SC(name, len, T1, T2, R1, R2, RApertures, EApertures, a, b, Nl, Nm, Nsteps, "DRIFT_SC")
end


"""
    KQUAD(;name::String = "Quad", len::Float64 = 0.0, k1::Float64 = 0.0, 
        PolynomA::Array{Float64,1} = zeros(Float64, 4), PolynomB::Array{Float64,1} = zeros(Float64, 4), 
        MaxOrder::Int64=1, NumIntSteps::Int64 = 10, rad::Int64=0, FringeQuadEntrance::Int64 = 0, 
        FringeQuadExit::Int64 = 0, T1::Array{Float64,1} = zeros(Float64, 6), 
        T2::Array{Float64,1} = zeros(Float64, 6), R1::Array{Float64,2} = zeros(Float64, 6, 6), 
        R2::Array{Float64,2} = zeros(Float64, 6, 6), RApertures::Array{Float64,1} = zeros(Float64, 6), 
        EApertures::Array{Float64,1} = zeros(Float64, 6), KickAngle::Array{Float64,1} = zeros(Float64, 2))

A canonical quadrupole element.

Example:
```julia
quad = KQUAD(name="Q1", len=0.5, k1=0.5)
```
"""
mutable struct KQUAD{T} <: AbstractElement{T}
    name::String
    len::T
    k1::T
    PolynomA::Array{T,1}
    PolynomB::Array{T,1}
    MaxOrder::Int64
    NumIntSteps::Int64
    rad::Int64
    FringeQuadEntrance::Int64
    FringeQuadExit::Int64
    T1::Array{T,1}
    T2::Array{T,1}
    R1::Array{T,2}
    R2::Array{T,2}
    RApertures::Array{Float64,1}
    EApertures::Array{Float64,1}
    KickAngle::Array{T,1}
    eletype::String

    # constructor with all parameters
    KQUAD(name::String, len::T, k1::T, PolynomA::Array{T,1}, 
        PolynomB::Array{T,1}, MaxOrder::Int64, NumIntSteps::Int64, rad::Int64, 
        FringeQuadEntrance::Int64, FringeQuadExit::Int64, T1::Array{T,1}, 
        T2::Array{T,1}, R1::Array{T,2}, R2::Array{T,2}, 
        RApertures::Array{Float64,1}, EApertures::Array{Float64,1}, KickAngle::Array{T,1}, eletype::String) where T =
        new{T}(name, len, k1, PolynomA, PolynomB, MaxOrder, NumIntSteps, rad, FringeQuadEntrance, FringeQuadExit,
            T1, T2, R1, R2, RApertures, EApertures, KickAngle, eletype)
end
function KQUAD(;name::String = "Quad", len = 0.0, k1 = 0.0, 
                PolynomA = zeros(4), 
                PolynomB = zeros(4), MaxOrder::Int64=1, 
                NumIntSteps::Int64 = 10, rad::Int64=0, FringeQuadEntrance::Int64 = 0, 
                FringeQuadExit::Int64 = 0, T1 = zeros(6), 
                T2 = zeros(6), R1 = zeros(6, 6), 
                R2 = zeros(6, 6), RApertures = zeros(6), 
                EApertures = zeros(6), KickAngle = zeros(2))
    
    if len isa DTPSAD || PolynomA[1] isa DTPSAD || PolynomB[1] isa DTPSAD || k1 isa DTPSAD ||
        T1[1] isa DTPSAD || T2[1] isa DTPSAD || R1[1,1] isa DTPSAD || R2[1,1] isa DTPSAD || KickAngle[1] isa DTPSAD
        if k1 != 0.0 && PolynomB[2] == 0.0
            PolynomB_T = [DTPSAD(PolynomB[1]), DTPSAD(k1), DTPSAD(PolynomB[3]), DTPSAD(PolynomB[4])]
            return KQUAD(name, DTPSAD(len), DTPSAD(k1), 
                DTPSAD.(PolynomA), PolynomB_T, MaxOrder, NumIntSteps, rad, FringeQuadEntrance, FringeQuadExit,
                DTPSAD.(T1), DTPSAD.(T2), DTPSAD.(R1), DTPSAD.(R2), RApertures, EApertures, DTPSAD.(KickAngle), "KQUAD")
        else
            return KQUAD(name, DTPSAD(len), DTPSAD(k1), 
                DTPSAD.(PolynomA), DTPSAD.(PolynomB), MaxOrder, NumIntSteps, rad, FringeQuadEntrance, FringeQuadExit,
                DTPSAD.(T1), DTPSAD.(T2), DTPSAD.(R1), DTPSAD.(R2), RApertures, EApertures, DTPSAD.(KickAngle), "KQUAD")
        end
    end

    if k1 != 0.0 && PolynomB[2] == 0.0
        PolynomB_T = [PolynomB[1], k1, PolynomB[3], PolynomB[4]]
        return KQUAD(name, len, k1, 
            PolynomA, PolynomB_T, MaxOrder, NumIntSteps, rad, FringeQuadEntrance, FringeQuadExit,
            T1, T2, R1, R2, RApertures, EApertures, KickAngle, "KQUAD")
    end
    return KQUAD(name, len, k1, 
        PolynomA, PolynomB, MaxOrder, NumIntSteps, rad, FringeQuadEntrance, FringeQuadExit,
        T1, T2, R1, R2, RApertures, EApertures, KickAngle, "KQUAD")
end

"""
    KQUAD_SC(;name::String = "Quad", len::Float64 = 0.0, k1::Float64 = 0.0, 
        PolynomA::Array{Float64,1} = zeros(Float64, 4), PolynomB::Array{Float64,1} = zeros(Float64, 4), 
        MaxOrder::Int64=1, NumIntSteps::Int64 = 10, rad::Int64=0, FringeQuadEntrance::Int64 = 0, 
        FringeQuadExit::Int64 = 0, T1::Array{Float64,1} = zeros(Float64, 6), 
        T2::Array{Float64,1} = zeros(Float64, 6), R1::Array{Float64,2} = zeros(Float64, 6, 6), 
        R2::Array{Float64,2} = zeros(Float64, 6, 6), RApertures::Array{Float64,1} = zeros(Float64, 6), 
        EApertures::Array{Float64,1} = zeros(Float64, 6), KickAngle::Array{Float64,1} = zeros(Float64, 2),
        a::Float64 = 1.0, b::Float64 = 1.0, Nl::Int64 = 10, Nm::Int64 = 10, Nsteps::Int64=1)

A canonical quadrupole element with space charge.
Example:
```julia
quad = KQUAD_SC(name="Q1_SC", len=0.5, k1=0.5, a=13e-3, b=13e-3, Nl=15, Nm=15)
```
"""
mutable struct KQUAD_SC{T} <: AbstractElement{T}
    name::String
    len::T
    k1::T
    PolynomA::Array{T,1}
    PolynomB::Array{T,1}
    MaxOrder::Int64
    NumIntSteps::Int64
    rad::Int64
    FringeQuadEntrance::Int64
    FringeQuadExit::Int64
    T1::Array{T,1}
    T2::Array{T,1}
    R1::Array{T,2}
    R2::Array{T,2}
    RApertures::Array{Float64,1}
    EApertures::Array{Float64,1}
    KickAngle::Array{T,1}
    a::T
    b::T
    Nl::Int64
    Nm::Int64
    Nsteps::Int64
    eletype::String

    # constructor with all parameters
    KQUAD_SC(name::String, len::T, k1::T, PolynomA::Array{T,1}, 
        PolynomB::Array{T,1}, MaxOrder::Int64, NumIntSteps::Int64, rad::Int64, 
        FringeQuadEntrance::Int64, FringeQuadExit::Int64, T1::Array{T,1}, 
        T2::Array{T,1}, R1::Array{T,2}, R2::Array{T,2}, 
        RApertures::Array{Float64,1}, EApertures::Array{Float64,1}, KickAngle::Array{T,1},
        a::T, b::T, Nl::Int64, Nm::Int64, Nsteps::Int64, eletype::String) where T =
        new{T}(name, len, k1, PolynomA, PolynomB, MaxOrder, NumIntSteps, rad, FringeQuadEntrance,
            FringeQuadExit, T1, T2, R1, R2, RApertures, EApertures, KickAngle,
            a, b, Nl, Nm, Nsteps, eletype)
end

function KQUAD_SC(;name::String = "Quad", len = 0.0, k1 = 0.0, 
                PolynomA = zeros(4), 
                PolynomB = zeros(4), MaxOrder::Int64=1, 
                NumIntSteps::Int64 = 10, rad::Int64=0, FringeQuadEntrance::Int64 = 0, 
                FringeQuadExit::Int64 = 0, T1 = zeros(6), 
                T2 = zeros(6), R1 = zeros(6, 6), 
                R2 = zeros(6, 6), RApertures = zeros(6), 
                EApertures = zeros(6), KickAngle = zeros(2),
                a::Float64 = 1.0, b::Float64 = 1.0, Nl::Int64 = 10, Nm::Int64 = 10, Nsteps::Int64=1)
    if len isa DTPSAD || PolynomA[1] isa DTPSAD || PolynomB[1] isa DTPSAD || k1 isa DTPSAD ||
        T1[1] isa DTPSAD || T2[1] isa DTPSAD || R1[1,1] isa DTPSAD || R2[1,1] isa DTPSAD || KickAngle[1] isa DTPSAD ||
        a isa DTPSAD || b isa DTPSAD
        if k1 != 0.0 && PolynomB[2] == 0.0
            PolynomB_T = [DTPSAD(PolynomB[1]), DTPSAD(k1), DTPSAD(PolynomB[3]), DTPSAD(PolynomB[4])]
            return KQUAD_SC(name, DTPSAD(len), DTPSAD(k1), 
                DTPSAD.(PolynomA), PolynomB_T, MaxOrder, NumIntSteps, rad, FringeQuadEntrance, FringeQuadExit,
                DTPSAD.(T1), DTPSAD.(T2), DTPSAD.(R1), DTPSAD.(R2), RApertures, EApertures, 
                DTPSAD.(KickAngle), DTPSAD(a), DTPSAD(b), Nl, Nm, Nsteps, "KQUAD_SC")
        end
        return KQUAD_SC(name, DTPSAD(len), DTPSAD(k1), 
            DTPSAD.(PolynomA), DTPSAD.(PolynomB), MaxOrder, NumIntSteps, rad, FringeQuadEntrance, FringeQuadExit,
            DTPSAD.(T1), DTPSAD.(T2), DTPSAD.(R1), DTPSAD.(R2), RApertures, EApertures, 
            DTPSAD.(KickAngle), DTPSAD(a), DTPSAD(b), Nl, Nm, Nsteps, "KQUAD_SC")
    end
    if k1 != 0.0 && PolynomB[2] == 0.0
        PolynomB_T = [PolynomB[1], k1, PolynomB[3], PolynomB[4]]
        return KQUAD_SC(name, len, k1, 
            PolynomA, PolynomB_T, MaxOrder, NumIntSteps, rad, FringeQuadEntrance, FringeQuadExit,
            T1, T2, R1, R2, RApertures, EApertures, KickAngle, a, b, Nl, Nm, Nsteps, "KQUAD_SC")
    end
    return KQUAD_SC(name, len, k1, 
        PolynomA, PolynomB, MaxOrder, NumIntSteps, rad, FringeQuadEntrance, FringeQuadExit,
        T1, T2, R1, R2, RApertures, EApertures, KickAngle, a, b, Nl, Nm, Nsteps, "KQUAD_SC")
end

"""
    KSEXT(;name::String = "Sext", len::Float64 = 0.0, k2::Float64 = 0.0, 
        PolynomA::Array{Float64,1} = zeros(Float64, 4), PolynomB::Array{Float64,1} = zeros(Float64, 4), 
        MaxOrder::Int64=2, NumIntSteps::Int64 = 10, rad::Int64=0, FringeQuadEntrance::Int64 = 0, 
        FringeQuadExit::Int64 = 0, T1::Array{Float64,1} = zeros(Float64, 6), 
        T2::Array{Float64,1} = zeros(Float64, 6), R1::Array{Float64,2} = zeros(Float64, 6, 6), 
        R2::Array{Float64,2} = zeros(Float64, 6, 6), RApertures::Array{Float64,1} = zeros(Float64, 6), 
        EApertures::Array{Float64,1} = zeros(Float64, 6), KickAngle::Array{Float64,1} = zeros(Float64, 2))

A canonical sextupole element.
Example:
```julia
sext = KSEXT(name="S1", len=0.5, k2=0.5)
```
"""
mutable struct KSEXT{T} <: AbstractElement{T}
    name::String
    len::T
    k2::T
    PolynomA::Array{T,1}
    PolynomB::Array{T,1}
    MaxOrder::Int64
    NumIntSteps::Int64
    rad::Int64
    FringeQuadEntrance::Int64
    FringeQuadExit::Int64
    T1::Array{T,1}
    T2::Array{T,1}
    R1::Array{T,2}
    R2::Array{T,2}
    RApertures::Array{Float64,1}
    EApertures::Array{Float64,1}
    KickAngle::Array{T,1}
    eletype::String
    # constructor with all parameters
    KSEXT(name::String, len::T, k2::T, PolynomA::Array{T,1}, 
        PolynomB::Array{T,1}, MaxOrder::Int64, NumIntSteps::Int64, rad::Int64, 
        FringeQuadEntrance::Int64, FringeQuadExit::Int64, T1::Array{T,1}, 
        T2::Array{T,1}, R1::Array{T,2}, R2::Array{T,2}, 
        RApertures::Array{Float64,1}, EApertures::Array{Float64,1}, KickAngle::Array{T,1}, eletype::String) where T= 
        new{T}(name, len, k2, PolynomA, PolynomB, MaxOrder, NumIntSteps, rad, FringeQuadEntrance,
            FringeQuadExit, T1, T2, R1, R2, RApertures, EApertures, KickAngle, eletype)
end
function KSEXT(;name::String = "Sext", len = 0.0, k2 = 0.0, 
                PolynomA = zeros(4), 
                PolynomB = zeros(4), MaxOrder::Int64=2, 
                NumIntSteps::Int64 = 10, rad::Int64=0, FringeQuadEntrance::Int64 = 0, 
                FringeQuadExit::Int64 = 0, T1 = zeros(6), 
                T2 = zeros(6), R1 = zeros(6, 6), 
                R2 = zeros(6, 6), RApertures = zeros(6), 
                EApertures = zeros(6), KickAngle = zeros(2))

    if len isa DTPSAD || PolynomA[1] isa DTPSAD || PolynomB[1] isa DTPSAD || k2 isa DTPSAD ||
        T1[1] isa DTPSAD || T2[1] isa DTPSAD || R1[1,1] isa DTPSAD || R2[1,1] isa DTPSAD || KickAngle[1] isa DTPSAD
        if k2 != 0.0 && PolynomB[3] == 0.0
            PolynomB_T = [DTPSAD(PolynomB[1]), DTPSAD(PolynomB[2]), DTPSAD(k2), DTPSAD(PolynomB[4])]
            return KSEXT(name, DTPSAD(len), DTPSAD(k2), 
                DTPSAD.(PolynomA), PolynomB_T, MaxOrder, NumIntSteps, rad, FringeQuadEntrance, FringeQuadExit,
                DTPSAD.(T1), DTPSAD.(T2), DTPSAD.(R1), DTPSAD.(R2), RApertures, EApertures, DTPSAD.(KickAngle), "KSEXT")
        else
            return KSEXT(name, DTPSAD(len), DTPSAD(k2), 
                DTPSAD.(PolynomA), DTPSAD.(PolynomB), MaxOrder, NumIntSteps, rad, FringeQuadEntrance, FringeQuadExit,
                DTPSAD.(T1), DTPSAD.(T2), DTPSAD.(R1), DTPSAD.(R2), RApertures, EApertures, DTPSAD.(KickAngle), "KSEXT")
        end
    end
    if k2 != 0.0 && PolynomB[3] == 0.0
        PolynomB_T = [PolynomB[1], PolynomB[2], k2, PolynomB[4]]
        return KSEXT(name, len, k2, 
            PolynomA, PolynomB_T, MaxOrder, NumIntSteps, rad, FringeQuadEntrance, FringeQuadExit,
            T1, T2, R1, R2, RApertures, EApertures, KickAngle, "KSEXT")
    end
    return KSEXT(name, len, k2, 
        PolynomA, PolynomB, MaxOrder, NumIntSteps, rad, FringeQuadEntrance, FringeQuadExit,
        T1, T2, R1, R2, RApertures, EApertures, KickAngle, "KSEXT")
end

"""
    KSEXT_SC(;name::String = "Sext", len::Float64 = 0.0, k2::Float64 = 0.0, 
        PolynomA::Array{Float64,1} = zeros(Float64, 4), PolynomB::Array{Float64,1} = zeros(Float64, 4), 
        MaxOrder::Int64=2, NumIntSteps::Int64 = 10, rad::Int64=0, FringeQuadEntrance::Int64 = 0, 
        FringeQuadExit::Int64 = 0, T1::Array{Float64,1} = zeros(Float64, 6), 
        T2::Array{Float64,1} = zeros(Float64, 6), R1::Array{Float64,2} = zeros(Float64, 6, 6), 
        R2::Array{Float64,2} = zeros(Float64, 6, 6), RApertures::Array{Float64,1} = zeros(Float64, 6),
        EApertures::Array{Float64,1} = zeros(Float64, 6), KickAngle::Array{Float64,1} = zeros(Float64, 2),
        a::Float64 = 1.0, b::Float64 = 1.0, Nl::Int64 = 10, Nm::Int64 = 10, Nsteps::Int64=1)

A canonical sextupole element with space charge.
Example:
```julia
sext = KSEXT_SC(name="S1_SC", len=0.5, k2=0.5, a=13e-3, b=13e-3, Nl=15, Nm=15)
```
"""
mutable struct KSEXT_SC{T} <: AbstractElement{T}
    name::String
    len::T
    k2::T
    PolynomA::Array{T,1}
    PolynomB::Array{T,1}
    MaxOrder::Int64
    NumIntSteps::Int64
    rad::Int64
    FringeQuadEntrance::Int64
    FringeQuadExit::Int64
    T1::Array{T,1}
    T2::Array{T,1}
    R1::Array{T,2}
    R2::Array{T,2}
    RApertures::Array{Float64,1}
    EApertures::Array{Float64,1}
    KickAngle::Array{T,1}
    a::T
    b::T
    Nl::Int64
    Nm::Int64
    Nsteps::Int64
    eletype::String

    # constructor with all parameters
    KSEXT_SC(name::String, len::T, k2::T, PolynomA::Array{T,1}, 
        PolynomB::Array{T,1}, MaxOrder::Int64, NumIntSteps::Int64, rad::Int64, 
        FringeQuadEntrance::Int64, FringeQuadExit::Int64, T1::Array{T,1}, 
        T2::Array{T,1}, R1::Array{T,2}, R2::Array{T,2}, 
        RApertures::Array{Float64,1}, EApertures::Array{Float64,1}, KickAngle::Array{T,1}, 
        a::T, b::T, Nl::Int64, Nm::Int64, Nsteps::Int64, eletype::String) where T = 
        new{T}(name, len, k2, PolynomA, PolynomB, MaxOrder, NumIntSteps, rad, FringeQuadEntrance,
            FringeQuadExit, T1, T2, R1, R2, RApertures, EApertures, KickAngle,
            a, b, Nl, Nm, Nsteps, eletype)
end
function KSEXT_SC(;name::String = "Sext", len = 0.0, k2 = 0.0, 
                PolynomA = zeros(4), 
                PolynomB = zeros(4), MaxOrder::Int64=2, 
                NumIntSteps::Int64 = 10, rad::Int64=0, FringeQuadEntrance::Int64 = 0, 
                FringeQuadExit::Int64 = 0, T1 = zeros(6), 
                T2 = zeros(6), R1 = zeros(6, 6), 
                R2 = zeros(6, 6), RApertures = zeros(6),
                EApertures = zeros(6), KickAngle = zeros(2),
                a::Float64 = 1.0, b::Float64 = 1.0, Nl::Int64 = 10, Nm::Int64 = 10, Nsteps::Int64=1)
    if len isa DTPSAD || PolynomA[1] isa DTPSAD || PolynomB[1] isa DTPSAD || k2 isa DTPSAD ||
        T1[1] isa DTPSAD || T2[1] isa DTPSAD || R1[1,1] isa DTPSAD || R2[1,1] isa DTPSAD || KickAngle[1] isa DTPSAD ||
        a isa DTPSAD || b isa DTPSAD
        if k2 != 0.0 && PolynomB[3] == 0.0
            PolynomB_T = [DTPSAD(PolynomB[1]), DTPSAD(PolynomB[2]), DTPSAD(k2), DTPSAD(PolynomB[4])]
            return KSEXT_SC(name, DTPSAD(len), DTPSAD(k2), 
                DTPSAD.(PolynomA), PolynomB_T, MaxOrder, NumIntSteps, rad, FringeQuadEntrance, FringeQuadExit,
                DTPSAD.(T1), DTPSAD.(T2), DTPSAD.(R1), DTPSAD.(R2), RApertures, EApertures, 
                DTPSAD.(KickAngle), DTPSAD(a), DTPSAD(b), Nl, Nm, Nsteps, "KSEXT_SC")
        end
        return KSEXT_SC(name, DTPSAD(len), DTPSAD(k2), 
            DTPSAD.(PolynomA), DTPSAD.(PolynomB), MaxOrder, NumIntSteps, rad, FringeQuadEntrance, FringeQuadExit,
            DTPSAD.(T1), DTPSAD.(T2), DTPSAD.(R1), DTPSAD.(R2), RApertures, EApertures, 
            DTPSAD.(KickAngle), DTPSAD(a), DTPSAD(b), Nl, Nm, Nsteps, "KSEXT_SC")
    end
    if k2 != 0.0 && PolynomB[3] == 0.0
        PolynomB_T = [PolynomB[1], PolynomB[2], k2, PolynomB[4]]
        return KSEXT_SC(name, len, k2, 
            PolynomA, PolynomB_T, MaxOrder, NumIntSteps, rad, FringeQuadEntrance, FringeQuadExit,
            T1, T2, R1, R2, RApertures, EApertures, KickAngle, a, b, Nl, Nm, Nsteps, "KSEXT_SC")
    end
    return KSEXT_SC(name, len, k2, 
        PolynomA, PolynomB, MaxOrder, NumIntSteps, rad, FringeQuadEntrance, FringeQuadExit,
        T1, T2, R1, R2, RApertures, EApertures, KickAngle, a, b, Nl, Nm, Nsteps, "KSEXT_SC")
end
"""
    KOCT(;name::String = "OCT", len::Float64 = 0.0, k3::Float64 = 0.0, 
        PolynomA::Array{Float64,1} = zeros(Float64, 4), PolynomB::Array{Float64,1} = zeros(Float64, 4), 
        MaxOrder::Int64=3, NumIntSteps::Int64 = 10, rad::Int64=0, FringeQuadEntrance::Int64 = 0, 
        FringeQuadExit::Int64 = 0, T1::Array{Float64,1} = zeros(Float64, 6), 
        T2::Array{Float64,1} = zeros(Float64, 6), R1::Array{Float64,2} = zeros(Float64, 6, 6), 
        R2::Array{Float64,2} = zeros(Float64, 6, 6), RApertures::Array{Float64,1} = zeros(Float64, 6), 
        EApertures::Array{Float64,1} = zeros(Float64, 6), KickAngle::Array{Float64,1} = zeros(Float64, 2))

A canonical octupole element.
Example:
```julia
oct = KOCT(name="O1", len=0.5, k3=0.5)
```
"""
mutable struct KOCT{T} <: AbstractElement{T}
    name::String
    len::T
    k3::T
    PolynomA::Array{T,1}
    PolynomB::Array{T,1}
    MaxOrder::Int64
    NumIntSteps::Int64
    rad::Int64
    FringeQuadEntrance::Int64
    FringeQuadExit::Int64
    T1::Array{T,1}
    T2::Array{T,1}
    R1::Array{T,2}
    R2::Array{T,2}
    RApertures::Array{Float64,1}
    EApertures::Array{Float64,1}
    KickAngle::Array{T,1}
    eletype::String

    # constructor with all parameters
    KOCT(name::String, len::T, k3::T, PolynomA::Array{T,1}, 
        PolynomB::Array{T,1}, MaxOrder::Int64, NumIntSteps::Int64, rad::Int64, 
        FringeQuadEntrance::Int64, FringeQuadExit::Int64, T1::Array{T,1}, 
        T2::Array{T,1}, R1::Array{T,2}, R2::Array{T,2}, 
        RApertures::Array{Float64,1}, EApertures::Array{Float64,1}, KickAngle::Array{T,1}, eletype::String) where T = 
        new{T}(name, len, k3, PolynomA, PolynomB, MaxOrder, NumIntSteps, rad, FringeQuadEntrance,
            FringeQuadExit, T1, T2, R1, R2, RApertures, EApertures, KickAngle, eletype)
end
function KOCT(;name::String = "OCT", len=0.0, k3=0.0, 
                PolynomA=zeros(4), 
                PolynomB=zeros(4), MaxOrder::Int64=3, 
                NumIntSteps::Int64 = 10, rad::Int64=0, FringeQuadEntrance::Int64 = 0, 
                FringeQuadExit::Int64 = 0, T1 = zeros(6), 
                T2 = zeros(6), R1 = zeros(6, 6), 
                R2 = zeros(6, 6), RApertures = zeros(6), 
                EApertures = zeros(6), KickAngle = zeros(2))
    if len isa DTPSAD || PolynomA[1] isa DTPSAD || PolynomB[1] isa DTPSAD || k3 isa DTPSAD ||
        T1[1] isa DTPSAD || T2[1] isa DTPSAD || R1[1,1] isa DTPSAD || R2[1,1] isa DTPSAD || KickAngle[1] isa DTPSAD
        if k3 != 0.0 && PolynomB[4] == 0
            PolynomB_T = [DTPSAD(PolynomB[1]), DTPSAD(PolynomB[2]), DTPSAD(PolynomB[3]), DTPSAD(k3)]
            return KOCT(name, DTPSAD(len), DTPSAD(k3), 
                DTPSAD.(PolynomA), PolynomB_T, MaxOrder, NumIntSteps, rad, FringeQuadEntrance, FringeQuadExit,
                DTPSAD.(T1), DTPSAD.(T2), DTPSAD.(R1), DTPSAD.(R2), RApertures, EApertures, DTPSAD.(KickAngle), "KOCT")
        else
            return KOCT(name, DTPSAD(len), DTPSAD(k3), 
                DTPSAD.(PolynomA), DTPSAD.(PolynomB), MaxOrder, NumIntSteps, rad, FringeQuadEntrance, FringeQuadExit,
                DTPSAD.(T1), DTPSAD.(T2), DTPSAD.(R1), DTPSAD.(R2), RApertures, EApertures, DTPSAD.(KickAngle), "KOCT")
        end
    end
    if k3 != 0.0 && PolynomB[4] == 0
        PolynomB_T = [PolynomB[1], PolynomB[2], PolynomB[3], k3]
        return KOCT(name, len, k3, 
            PolynomA, PolynomB_T, MaxOrder, NumIntSteps, rad, FringeQuadEntrance, FringeQuadExit,
            T1, T2, R1, R2, RApertures, EApertures, KickAngle, "KOCT")
    end
    return KOCT(name, len, k3, 
        PolynomA, PolynomB, MaxOrder, NumIntSteps, rad, FringeQuadEntrance, FringeQuadExit,
        T1, T2, R1, R2, RApertures, EApertures, KickAngle, "KOCT")
end

"""
    KOCT_SC(;name::String = "OCT", len::Float64 = 0.0, k3::Float64 = 0.0, 
        PolynomA::Array{Float64,1} = zeros(Float64, 4), PolynomB::Array{Float64,1} = zeros(Float64, 4), 
        MaxOrder::Int64=3, NumIntSteps::Int64 = 10, rad::Int64=0, FringeQuadEntrance::Int64 = 0, 
        FringeQuadExit::Int64 = 0, T1::Array{Float64,1} = zeros(Float64, 6), 
        T2::Array{Float64,1} = zeros(Float64, 6), R1::Array{Float64,2} = zeros(Float64, 6, 6), 
        R2::Array{Float64,2} = zeros(Float64, 6, 6), RApertures::Array{Float64,1} = zeros(Float64, 6),
        EApertures::Array{Float64,1} = zeros(Float64, 6), KickAngle::Array{Float64,1} = zeros(Float64, 2),
        a::Float64 = 1.0, b::Float64 = 1.0, Nl::Int64 = 10, Nm::Int64 = 10, Nsteps::Int64=1)

A canonical octupole element with space charge.
Example:
```julia
oct = KOCT_SC(name="O1_SC", len=0.5, k3=0.5, a=13e-3, b=13e-3, Nl=15, Nm=15)
```
"""
mutable struct KOCT_SC{T} <: AbstractElement{T}
    name::String
    len::T
    k3::T
    PolynomA::Array{T,1}
    PolynomB::Array{T,1}
    MaxOrder::Int64
    NumIntSteps::Int64
    rad::Int64
    FringeQuadEntrance::Int64
    FringeQuadExit::Int64
    T1::Array{T,1}
    T2::Array{T,1}
    R1::Array{T,2}
    R2::Array{T,2}
    RApertures::Array{Float64,1}
    EApertures::Array{Float64,1}
    KickAngle::Array{T,1}
    a::T
    b::T
    Nl::Int64
    Nm::Int64
    Nsteps::Int64
    eletype::String

    # constructor with all parameters
    KOCT_SC(name::String, len::T, k3::T, PolynomA::Array{T,1}, 
        PolynomB::Array{T,1}, MaxOrder::Int64, NumIntSteps::Int64, rad::Int64, 
        FringeQuadEntrance::Int64, FringeQuadExit::Int64, T1::Array{T,1}, 
        T2::Array{T,1}, R1::Array{T,2}, R2::Array{T,2}, 
        RApertures::Array{Float64,1}, EApertures::Array{Float64,1}, KickAngle::Array{T,1}, 
        a::T, b::T, Nl::Int64, Nm::Int64, Nsteps::Int64, eletype::String) where T = 
        new{T}(name, len, k3, PolynomA, PolynomB, MaxOrder, NumIntSteps, rad, FringeQuadEntrance,
            FringeQuadExit, T1, T2, R1, R2, RApertures, EApertures, KickAngle, a, b, Nl, Nm, Nsteps, eletype)
end
function KOCT_SC(;name::String = "OCT", len = 0.0, k3 = 0.0, 
                PolynomA = zeros(4), 
                PolynomB = zeros(4), MaxOrder::Int64=3, 
                NumIntSteps::Int64 = 10, rad::Int64=0, FringeQuadEntrance::Int64 = 0, 
                FringeQuadExit::Int64 = 0, T1 = zeros(6), 
                T2 = zeros(6), R1 = zeros(6, 6), 
                R2 = zeros(6, 6), RApertures = zeros(6),
                EApertures = zeros(6), KickAngle = zeros(2),
                a::Float64 = 1.0, b::Float64 = 1.0, Nl::Int64 = 10, Nm::Int64 = 10, Nsteps::Int64=1)
    if len isa DTPSAD || PolynomA[1] isa DTPSAD || PolynomB[1] isa DTPSAD || k3 isa DTPSAD ||
        T1[1] isa DTPSAD || T2[1] isa DTPSAD || R1[1,1] isa DTPSAD || R2[1,1] isa DTPSAD || KickAngle[1] isa DTPSAD ||
        a isa DTPSAD || b isa DTPSAD
        if k3 != 0.0 && PolynomB[4] == 0
            PolynomB_T = [DTPSAD(PolynomB[1]), DTPSAD(PolynomB[2]), DTPSAD(PolynomB[3]), DTPSAD(k3)]
            return KOCT_SC(name, DTPSAD(len), DTPSAD(k3), 
                DTPSAD.(PolynomA), PolynomB_T, MaxOrder, NumIntSteps, rad, FringeQuadEntrance, FringeQuadExit,
                DTPSAD.(T1), DTPSAD.(T2), DTPSAD.(R1), DTPSAD.(R2), RApertures, EApertures, 
                DTPSAD.(KickAngle), DTPSAD(a), DTPSAD(b), Nl, Nm, Nsteps, "KOCT_SC")
        end
        return KOCT_SC(name, DTPSAD(len), DTPSAD(k3), 
            DTPSAD.(PolynomA), DTPSAD.(PolynomB), MaxOrder, NumIntSteps, rad, FringeQuadEntrance, FringeQuadExit,
            DTPSAD.(T1), DTPSAD.(T2), DTPSAD.(R1), DTPSAD.(R2), RApertures, EApertures, 
            DTPSAD.(KickAngle), DTPSAD(a), DTPSAD(b), Nl, Nm, Nsteps, "KOCT_SC")
    end
    if k3 != 0.0 && PolynomB[4] == 0.0
        PolynomB_T = [PolynomB[1], PolynomB[2], PolynomB[3], k3]
        return KOCT_SC(name, len, k3, 
            PolynomA, PolynomB_T, MaxOrder, NumIntSteps, rad, FringeQuadEntrance, FringeQuadExit,
            T1, T2, R1, R2, RApertures, EApertures, KickAngle, a, b, Nl, Nm, Nsteps, "KOCT_SC")
    end
    return KOCT_SC(name, len, k3, 
        PolynomA, PolynomB, MaxOrder, NumIntSteps, rad, FringeQuadEntrance, FringeQuadExit,
        T1, T2, R1, R2, RApertures, EApertures, KickAngle, a, b, Nl, Nm, Nsteps, "KOCT_SC")
end

"""
    thinMULTIPOLE(;name::String = "thinMULTIPOLE", len::Float64 = 0.0, PolynomA::Array{Float64,1} = zeros(Float64, 4), 
        PolynomB::Array{Float64,1} = zeros(Float64, 4), MaxOrder::Int64=1, NumIntSteps::Int64 = 1, rad::Int64=0, 
        FringeQuadEntrance::Int64 = 0, FringeQuadExit::Int64 = 0, T1::Array{Float64,1} = zeros(Float64, 6), 
        T2::Array{Float64,1} = zeros(Float64, 6), R1::Array{Float64,2} = zeros(Float64, 6, 6), 
        R2::Array{Float64,2} = zeros(Float64, 6, 6), RApertures::Array{Float64,1} = zeros(Float64, 6), 
        EApertures::Array{Float64,1} = zeros(Float64, 6), KickAngle::Array{Float64,1} = zeros(Float64, 2))

A thin multipole element.
PolynomA and PolynomB are the skew and normal components of the multipole.

Example:
```julia
multipole = thinMULTIPOLE(name="M1", len=0.5, PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
```
"""
mutable struct thinMULTIPOLE{T} <: AbstractElement{T}
    name::String
    len::T
    PolynomA::Array{T,1}
    PolynomB::Array{T,1}
    MaxOrder::Int64
    NumIntSteps::Int64
    rad::Int64
    FringeQuadEntrance::Int64
    FringeQuadExit::Int64
    T1::Array{T,1}
    T2::Array{T,1}
    R1::Array{T,2}
    R2::Array{T,2}
    RApertures::Array{Float64,1}
    EApertures::Array{Float64,1}
    KickAngle::Array{T,1}
    eletype::String
    # constructor with all parameters
    thinMULTIPOLE(name::String, len::T, PolynomA::Array{T,1}, 
    PolynomB::Array{T,1}, MaxOrder::Int64, NumIntSteps::Int64, rad::Int64, FringeQuadEntrance::Int64, 
    FringeQuadExit::Int64, T1::Array{T,1}, T2::Array{T,1}, R1::Array{T,2}, R2::Array{T,2}, 
    RApertures::Array{Float64,1}, EApertures::Array{Float64,1}, KickAngle::Array{T,1}, eletype::String) where T = new{T}(name, len, 
    PolynomA, PolynomB, MaxOrder, NumIntSteps, rad, FringeQuadEntrance, FringeQuadExit,
            T1, T2, R1, R2, RApertures, EApertures, KickAngle, eletype)
end
function thinMULTIPOLE(;name::String = "thinMULTIPOLE", len = 0.0, PolynomA = zeros(4), 
                PolynomB = zeros(4), MaxOrder::Int64=1, NumIntSteps::Int64 = 1, rad::Int64=0, 
                FringeQuadEntrance::Int64 = 0, FringeQuadExit::Int64 = 0, T1= zeros(6), 
                T2= zeros(6), R1= zeros(6, 6), R2= zeros(6, 6), RApertures= zeros(6), 
                EApertures= zeros(6), KickAngle= zeros(2))
    if PolynomB[3] != 0.0
        MaxOrder = 2
    end
    if PolynomB[4] != 0.0
        MaxOrder = 3
    end
    if len isa DTPSAD || PolynomA[1] isa DTPSAD || PolynomB[1] isa DTPSAD ||
        T1[1] isa DTPSAD || T2[1] isa DTPSAD || R1[1,1] isa DTPSAD || R2[1,1] isa DTPSAD || KickAngle[1] isa DTPSAD
        return thinMULTIPOLE(name, DTPSAD(len), 
            DTPSAD.(PolynomA), DTPSAD.(PolynomB), MaxOrder, NumIntSteps, rad, FringeQuadEntrance, FringeQuadExit,
            DTPSAD.(T1), DTPSAD.(T2), DTPSAD.(R1), DTPSAD.(R2), RApertures, EApertures, DTPSAD.(KickAngle), "thinMULTIPOLE")
    end
    return thinMULTIPOLE(name, len, 
        PolynomA, PolynomB, MaxOrder, NumIntSteps, rad, FringeQuadEntrance, FringeQuadExit,
        T1, T2, R1, R2, RApertures, EApertures, KickAngle, "thinMULTIPOLE")
end

"""
    SBEND(;name::String = "SBend", len::Float64 = 0.0, angle::Float64 = 0.0, e1::Float64 = 0.0, e2::Float64 = 0.0, 
        PolynomA::Array{Float64,1} = zeros(Float64, 4), PolynomB::Array{Float64,1} = zeros(Float64, 4), 
        MaxOrder::Int64=0, NumIntSteps::Int64 = 10, rad::Int64=0, fint1::Float64 = 0.0, fint2::Float64 = 0.0, 
        gap::Float64 = 0.0, FringeBendEntrance::Int64 = 1, FringeBendExit::Int64 = 1, FringeQuadEntrance::Int64 = 0, 
        FringeQuadExit::Int64 = 0, FringeIntM0::Array{Float64,1} = zeros(Float64, 5), FringeIntP0::Array{Float64,1} = zeros(Float64, 5), 
        T1::Array{Float64,1} = zeros(Float64, 6), T2::Array{Float64,1} = zeros(Float64, 6), R1::Array{Float64,2} = zeros(Float64, 6, 6), 
        R2::Array{Float64,2} = zeros(Float64, 6, 6), RApertures::Array{Float64,1} = zeros(Float64, 6), EApertures::Array{Float64,1} = zeros(Float64, 6), 
        KickAngle::Array{Float64,1} = zeros(Float64, 2))

A sector bending magnet.
Example:
```julia
bend = SBEND(name="B1", len=0.5, angle=0.5)
```
"""
mutable struct SBEND{T} <: AbstractElement{T}
    name::String
    len::T
    angle::T
    e1::T
    e2::T
    PolynomA::Array{T,1}
    PolynomB::Array{T,1}
    MaxOrder::Int64
    NumIntSteps::Int64
    rad::Int64
    fint1::T
    fint2::T
    gap::T
    FringeBendEntrance::Int64
    FringeBendExit::Int64
    FringeQuadEntrance::Int64
    FringeQuadExit::Int64
    FringeIntM0::Array{T,1}
    FringeIntP0::Array{T,1}
    T1::Array{T,1}
    T2::Array{T,1}
    R1::Array{T,2}
    R2::Array{T,2}
    RApertures::Array{Float64,1}
    EApertures::Array{Float64,1}
    KickAngle::Array{T,1}
    eletype::String
    # constructor with all parameters
    SBEND(name::String, len::T, angle::T, e1::T, e2::T, 
        PolynomA::Array{T,1}, PolynomB::Array{T,1}, MaxOrder::Int64, NumIntSteps::Int64, rad::Int64, 
        fint1::T, fint2::T, gap::T, FringeBendEntrance::Int64, FringeBendExit::Int64, 
        FringeQuadEntrance::Int64, FringeQuadExit::Int64, FringeIntM0::Array{T,1}, 
        FringeIntP0::Array{T,1}, T1::Array{T,1}, T2::Array{T,1}, R1::Array{T,2}, 
        R2::Array{T,2}, RApertures::Array{Float64,1}, EApertures::Array{Float64,1}, KickAngle::Array{T,1}, 
        eletype::String) where T = new{T}(name, len, angle, e1, e2, PolynomA, PolynomB, MaxOrder, NumIntSteps, rad, fint1, fint2,
            gap, FringeBendEntrance, FringeBendExit, FringeQuadEntrance,
            FringeQuadExit, FringeIntM0, FringeIntP0,
            T1, T2, R1, R2, RApertures,
            EApertures,
            KickAngle,
            eletype)
end
function SBEND(;name::String = "SBend", len = 0.0, angle = 0.0, e1 = 0.0, e2 = 0.0, 
                PolynomA = zeros(4), PolynomB = zeros(4), 
                MaxOrder::Int64=0, NumIntSteps::Int64 = 10, rad::Int64=0, fint1 = 0.0, fint2 = 0.0, 
                gap = 0.0, FringeBendEntrance::Int64 = 1, FringeBendExit::Int64 = 1, 
                FringeQuadEntrance::Int64 = 0, FringeQuadExit::Int64 = 0, FringeIntM0 = zeros(5), 
                FringeIntP0 = zeros(5), T1 = zeros(6), 
                T2 = zeros(6), R1 = zeros(6, 6), 
                R2 = zeros(6, 6), RApertures = zeros(6), 
                EApertures = zeros(6), KickAngle = zeros(2))
    if PolynomB[2] != 0.0
        MaxOrder = 1
    end
    if PolynomB[3] != 0.0
        MaxOrder = 2
    end
    if PolynomB[4] != 0.0
        MaxOrder = 3
    end
    if len isa DTPSAD || angle isa DTPSAD || e1 isa DTPSAD || e2 isa DTPSAD || PolynomA[1] isa DTPSAD || PolynomB[1] isa DTPSAD ||
        fint1 isa DTPSAD || fint2 isa DTPSAD || gap isa DTPSAD || FringeIntM0[1] isa DTPSAD || FringeIntP0[1] isa DTPSAD ||
        T1[1] isa DTPSAD || T2[1] isa DTPSAD || R1[1,1] isa DTPSAD || R2[1,1] isa DTPSAD || KickAngle[1] isa DTPSAD
        return SBEND(name, DTPSAD(len), DTPSAD(angle), DTPSAD(e1), DTPSAD(e2), 
            DTPSAD.(PolynomA), DTPSAD.(PolynomB), MaxOrder, NumIntSteps, rad, DTPSAD(fint1), DTPSAD(fint2), 
            DTPSAD(gap), FringeBendEntrance, FringeBendExit, FringeQuadEntrance, FringeQuadExit,
            DTPSAD.(FringeIntM0), DTPSAD.(FringeIntP0),
            DTPSAD.(T1), DTPSAD.(T2), DTPSAD.(R1), DTPSAD.(R2), RApertures,
            EApertures,
            DTPSAD.(KickAngle),
            "SBEND")
    end
    return SBEND(name, len, angle, e1, e2, 
        PolynomA, PolynomB, MaxOrder, NumIntSteps, rad, fint1, fint2, 
        gap, FringeBendEntrance, FringeBendExit, FringeQuadEntrance,
        FringeQuadExit, FringeIntM0, FringeIntP0,
        T1, T2, R1, R2, RApertures,
        EApertures,
        KickAngle,
        "SBEND")
end

"""
    SBEND_SC(;name::String = "SBend", len::Float64 = 0.0, angle::Float64 = 0.0, e1::Float64 = 0.0, e2::Float64 = 0.0, 
        PolynomA::Array{Float64,1} = zeros(Float64, 4), PolynomB::Array{Float64,1} = zeros(Float64, 4), 
        MaxOrder::Int64=0, NumIntSteps::Int64 = 10, rad::Int64=0, fint1::Float64 = 0.0, fint2::Float64 = 0.0, 
        gap::Float64 = 0.0, FringeBendEntrance::Int64 = 1, FringeBendExit::Int64 = 1, 
        FringeQuadEntrance::Int64 = 0, FringeQuadExit::Int64 = 0, FringeIntM0::Array{Float64,1} = zeros(Float64, 5), 
        FringeIntP0::Array{Float64,1} = zeros(Float64, 5), T1::Array{Float64,1} = zeros(Float64, 6), 
        T2::Array{Float64,1} = zeros(Float64, 6), R1::Array{Float64,2} = zeros(Float64, 6, 6), 
        R2::Array{Float64,2} = zeros(Float64, 6, 6), RApertures::Array{Float64,1} = zeros(Float64, 6), 
        EApertures::Array{Float64,1} = zeros(Float64, 6), KickAngle::Array{Float64,1} = zeros(Float64, 2),
        a::Float64 = 1.0, b::Float64 = 1.0, Nl::Int64 = 10, Nm::Int64 = 10, Nsteps::Int64=1)

A sector bending magnet with space charge.
Example:
```julia
bend = SBEND_SC(name="B1_SC", len=0.5, angle=0.5, a=13e-3, b=13e-3, Nl=15, Nm=15)
```
"""
mutable struct SBEND_SC{T} <: AbstractElement{T}
    name::String
    len::T
    angle::T
    e1::T
    e2::T
    PolynomA::Array{T,1}
    PolynomB::Array{T,1}
    MaxOrder::Int64
    NumIntSteps::Int64
    rad::Int64
    fint1::T
    fint2::T
    gap::T
    FringeBendEntrance::Int64
    FringeBendExit::Int64
    FringeQuadEntrance::Int64
    FringeQuadExit::Int64
    FringeIntM0::Array{T,1}
    FringeIntP0::Array{T,1}
    T1::Array{T,1}
    T2::Array{T,1}
    R1::Array{T,2}
    R2::Array{T,2}
    RApertures::Array{Float64,1}
    EApertures::Array{Float64,1}
    KickAngle::Array{T,1}
    a::T
    b::T
    Nl::Int64
    Nm::Int64
    Nsteps::Int64
    eletype::String
    # constructor with all parameters
    SBEND_SC(name::String, len::T, angle::T, e1::T, e2::T, PolynomA::Array{T,1}, 
    PolynomB::Array{T,1}, MaxOrder::Int64, NumIntSteps::Int64, rad::Int64, fint1::T, fint2::T, 
    gap::T, FringeBendEntrance::Int64, FringeBendExit::Int64, FringeQuadEntrance::Int64, FringeQuadExit::Int64, 
    FringeIntM0::Array{T,1}, FringeIntP0::Array{T,1}, T1::Array{T,1}, T2::Array{T,1}, 
    R1::Array{T,2}, R2::Array{T,2}, RApertures::Array{Float64,1}, EApertures::Array{Float64,1}, 
    KickAngle::Array{T,1}, a::T, b::T, Nl::Int64, Nm::Int64, Nsteps::Int64, eletype::String) where T = new{T}(name, len, 
    angle, e1, e2, PolynomA, PolynomB, MaxOrder, NumIntSteps, rad, fint1, fint2, gap, 
            FringeBendEntrance, FringeBendExit, FringeQuadEntrance, FringeQuadExit, FringeIntM0, FringeIntP0, 
            T1, T2, R1, R2, RApertures, EApertures, KickAngle, a, b, Nl, Nm, Nsteps, eletype)
end
function SBEND_SC(;name::String = "SBend", len = 0.0, angle = 0.0, e1 = 0.0, e2 = 0.0, 
                PolynomA = zeros(4), PolynomB = zeros(4), 
                MaxOrder::Int64=0, NumIntSteps::Int64 = 10, rad::Int64=0, fint1 = 0.0, fint2 = 0.0, 
                gap = 0.0, FringeBendEntrance::Int64 = 1, FringeBendExit::Int64 = 1, 
                FringeQuadEntrance::Int64 = 0, FringeQuadExit::Int64 = 0, FringeIntM0 = zeros(5), 
                FringeIntP0 = zeros(5), T1 = zeros(6), 
                T2 = zeros(6), R1 = zeros(6, 6), 
                R2 = zeros(6, 6), RApertures = zeros(6), 
                EApertures = zeros(6), KickAngle = zeros(2),
                a::Float64 = 1.0, b::Float64 = 1.0, Nl::Int64 = 10, Nm::Int64 = 10, Nsteps::Int64=1)
    if PolynomB[2] != 0.0
        MaxOrder = 1
    end
    if PolynomB[3] != 0.0
        MaxOrder = 2
    end
    if PolynomB[4] != 0.0
        MaxOrder = 3
    end
    if len isa DTPSAD || angle isa DTPSAD || e1 isa DTPSAD || e2 isa DTPSAD || PolynomA[1] isa DTPSAD || PolynomB[1] isa DTPSAD ||
        fint1 isa DTPSAD || fint2 isa DTPSAD || gap isa DTPSAD || FringeIntM0[1] isa DTPSAD || FringeIntP0[1] isa DTPSAD ||
        T1[1] isa DTPSAD || T2[1] isa DTPSAD || R1[1,1] isa DTPSAD || R2[1,1] isa DTPSAD || KickAngle[1] isa DTPSAD ||
        a isa DTPSAD || b isa DTPSAD
        return SBEND_SC(name, DTPSAD(len), DTPSAD(angle), DTPSAD(e1), DTPSAD(e2), 
            DTPSAD.(PolynomA), DTPSAD.(PolynomB), MaxOrder, NumIntSteps, rad, DTPSAD(fint1), DTPSAD(fint2), 
            DTPSAD(gap), FringeBendEntrance, FringeBendExit, FringeQuadEntrance, FringeQuadExit,
            DTPSAD.(FringeIntM0), DTPSAD.(FringeIntP0),
            DTPSAD.(T1), DTPSAD.(T2), DTPSAD.(R1), DTPSAD.(R2), RApertures,
            EApertures,
            DTPSAD.(KickAngle),
            DTPSAD(a), DTPSAD(b), Nl, Nm, Nsteps, "SBEND_SC")
    end
    return SBEND_SC(name, len, angle, e1, e2,
        PolynomA, PolynomB, MaxOrder, NumIntSteps, rad, fint1, fint2, 
        gap, FringeBendEntrance, FringeBendExit, FringeQuadEntrance,
        FringeQuadExit, FringeIntM0, FringeIntP0,
        T1, T2, R1, R2, RApertures,
        EApertures,
        KickAngle,
        a, b, Nl, Nm, Nsteps, "SBEND_SC")   
end

"""
    RBEND(;name::String = "RBend", len::Float64 = 0.0, angle::Float64 = 0.0, PolynomA::Array{Float64,1} = zeros(Float64, 4), 
        PolynomB::Array{Float64,1} = zeros(Float64, 4), MaxOrder::Int64=0, NumIntSteps::Int64 = 10, rad::Int64=0, 
        fint1::Float64 = 0.0, fint2::Float64 = 0.0, gap::Float64 = 0.0, FringeBendEntrance::Int64 = 1, 
        FringeBendExit::Int64 = 1, FringeQuadEntrance::Int64 = 0, FringeQuadExit::Int64 = 0, 
        FringeIntM0::Array{Float64,1} = zeros(Float64, 5), FringeIntP0::Array{Float64,1} = zeros(Float64, 5), 
        T1::Array{Float64,1} = zeros(Float64, 6), T2::Array{Float64,1} = zeros(Float64, 6), 
        R1::Array{Float64,2} = zeros(Float64, 6, 6), R2::Array{Float64,2} = zeros(Float64, 6, 6), 
        RApertures::Array{Float64,1} = zeros(Float64, 6), EApertures::Array{Float64,1} = zeros(Float64, 6), 
        KickAngle::Array{Float64,1} = zeros(Float64, 2))

A rectangular bending magnet.
Example:
```julia
bend = RBEND(name="B1", len=0.5, angle=0.5)
```
"""
function RBEND(;name::String = "RBend", len = 0.0, angle = 0.0, PolynomA = zeros(4), 
                PolynomB = zeros(4), MaxOrder::Int64=0, NumIntSteps::Int64 = 10, rad::Int64=0, fint1 = 0.0, 
                fint2 = 0.0, gap = 0.0, FringeBendEntrance::Int64 = 1, FringeBendExit::Int64 = 1, 
                FringeQuadEntrance::Int64 = 0, FringeQuadExit::Int64 = 0, FringeIntM0 = zeros(5), 
                FringeIntP0 = zeros(5), T1 = zeros(6), T2 = zeros(6), 
                R1 = zeros(6,6), R2 = zeros(6,6), RApertures = zeros(6), 
                EApertures = zeros(6), KickAngle = zeros(2))
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
    return SBEND(name=name, len=len, angle=angle, e1=e1, e2=e2, PolynomA=PolynomA, PolynomB=PolynomB, MaxOrder=MaxOrder, 
                NumIntSteps=NumIntSteps, rad=rad, fint1=fint1, fint2=fint2, gap=gap, FringeBendEntrance=FringeBendEntrance, 
                FringeBendExit=FringeBendExit, FringeQuadEntrance=FringeQuadEntrance, FringeQuadExit=FringeQuadExit, 
                FringeIntM0=FringeIntM0, FringeIntP0=FringeIntP0, T1=T1, T2=T2, R1=R1, R2=R2, RApertures=RApertures, 
                EApertures=EApertures, KickAngle=KickAngle)
end

function RBEND_SC(;name::String = "RBend", len::Float64 = 0.0, angle::Float64 = 0.0, PolynomA::Array{Float64,1} = zeros(4), 
                PolynomB::Array{Float64,1} = zeros(4), MaxOrder::Int64=0, NumIntSteps::Int64 = 10, rad::Int64=0, fint1::Float64 = 0.0, 
                fint2::Float64 = 0.0, gap::Float64 = 0.0, FringeBendEntrance::Int64 = 1, FringeBendExit::Int64 = 1, 
                FringeQuadEntrance::Int64 = 0, FringeQuadExit::Int64 = 0, FringeIntM0::Array{Float64,1} = zeros(5), 
                FringeIntP0::Array{Float64,1} = zeros(5), T1::Array{Float64,1} = zeros(6), T2::Array{Float64,1} = zeros(6), 
                R1::Array{Float64,2} = zeros(6,6), R2::Array{Float64,2} = zeros(6,6), RApertures::Array{Float64,1} = zeros(6), 
                EApertures::Array{Float64,1} = zeros(6), KickAngle::Array{Float64,1} = zeros(2), a::Float64 = 1.0, b::Float64 = 1.0, 
                Nl::Int64 = 10, Nm::Int64 = 10, Nsteps::Int64=1)
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
    return SBEND_SC(name=name, len=len, angle=angle, e1=e1, e2=e2, PolynomA=PolynomA, PolynomB=PolynomB, MaxOrder=MaxOrder, 
                NumIntSteps=NumIntSteps, rad=rad, fint1=fint1, fint2=fint2, gap=gap, FringeBendEntrance=FringeBendEntrance, 
                FringeBendExit=FringeBendExit, FringeQuadEntrance=FringeQuadEntrance, FringeQuadExit=FringeQuadExit, 
                FringeIntM0=FringeIntM0, FringeIntP0=FringeIntP0, T1=T1, T2=T2, R1=R1, R2=R2, RApertures=RApertures, 
                EApertures=EApertures, KickAngle=KickAngle, a=a, b=b, Nl=Nl, Nm=Nm, Nsteps=Nsteps)
end

"""
    LBEND(;name::String = "Bend", len::Float64 = 0.0, angle::Float64 = 0.0, e1::Float64 = 0.0, e2::Float64 = 0.0, K::Float64 = 0.0,
        fint1::Float64 = 0.0, fint2::Float64 = 0.0, FullGap::Float64 = 0.0,
        T1::Array{Float64,1} = zeros(Float64, 6), T2::Array{Float64,1} = zeros(Float64, 6), R1::Array{Float64,2} = zeros(Float64, 6, 6), 
        R2::Array{Float64,2} = zeros(Float64, 6, 6), RApertures::Array{Float64,1} = zeros(Float64, 6), EApertures::Array{Float64,1} = zeros(Float64, 6))

A sector bending magnet with linear map.
Example:
```julia
bend = LBEND(name="B1", len=0.5, angle=0.5)
```
"""
mutable struct LBEND{T} <: AbstractElement{T}
    name::String
    len::T
    angle::T
    e1::T
    e2::T
    K::T
    ByError::T
    fint1::T
    fint2::T
    FullGap::T
    T1::Array{T,1}
    T2::Array{T,1}
    R1::Array{T,2}
    R2::Array{T,2}
    RApertures::Array{Float64,1}
    EApertures::Array{Float64,1}
    eletype::String

    # constructor with all parameters
    LBEND(name::String, len::T, angle::T, e1::T, e2::T, 
    K::T, ByError::T, fint1::T, fint2::T, FullGap::T, 
    T1::Array{T,1}, T2::Array{T,1}, R1::Array{T,2}, R2::Array{T,2}, 
    RApertures::Array{Float64,1}, EApertures::Array{Float64,1}, eletype::String) where T = new{T}(name, len, angle, e1, e2, 
    K, ByError, fint1, fint2, FullGap, T1, T2, R1, R2, RApertures, EApertures, eletype)
end
function LBEND(;name::String = "LBEND", len = 0.0, angle = 0.0, e1 = 0.0, e2 = 0.0, 
    K = 0.0, ByError = 0.0, fint1 = 0.0, fint2 = 0.0, FullGap = 0.0, 
    T1 = zeros(Float64, 6), T2 = zeros(Float64, 6), 
    R1 = zeros(Float64, 6, 6), R2 = zeros(Float64, 6, 6), 
    RApertures = zeros(Float64, 6), EApertures = zeros(Float64, 6))
    if len isa DTPSAD || angle isa DTPSAD || e1 isa DTPSAD || e2 isa DTPSAD || K isa DTPSAD || ByError isa DTPSAD ||
        fint1 isa DTPSAD || fint2 isa DTPSAD || FullGap isa DTPSAD || T1[1] isa DTPSAD || T2[1] isa DTPSAD ||
        R1[1,1] isa DTPSAD || R2[1,1] isa DTPSAD
        return LBEND(name, DTPSAD(len), DTPSAD(angle), DTPSAD(e1), DTPSAD(e2), 
            DTPSAD(K), DTPSAD(ByError), DTPSAD(fint1), DTPSAD(fint2), DTPSAD(FullGap), 
            DTPSAD.(T1), DTPSAD.(T2), DTPSAD.(R1), DTPSAD.(R2), RApertures, EApertures, "LBEND")
    end
    return LBEND(name, len, angle, e1, e2, K, ByError, fint1, fint2, FullGap, T1, T2, R1, R2, RApertures, EApertures, "LBEND")
end


mutable struct ESBEND{T} <: AbstractElement{T}
    name::String
    len::T
    angle::T
    e1::T
    e2::T
    PolynomA::Array{T,1}
    PolynomB::Array{T,1}
    MaxOrder::Int64
    NumIntSteps::Int64
    rad::Int64
    gK::T
    FringeBendEntrance::Int64
    FringeBendExit::Int64
    FringeQuadEntrance::Int64
    FringeQuadExit::Int64
    T1::Array{T,1}
    T2::Array{T,1}
    R1::Array{T,2}
    R2::Array{T,2}
    RApertures::Array{Float64,1}
    EApertures::Array{Float64,1}
    KickAngle::Array{T,1}
    eletype::String
    # constructor with all parameters
    ESBEND(name::String, len::T, angle::T, e1::T, e2::T, 
        PolynomA::Array{T,1}, PolynomB::Array{T,1}, MaxOrder::Int64, 
        NumIntSteps::Int64, rad::Int64, gK::T, FringeBendEntrance::Int64, 
        FringeBendExit::Int64, FringeQuadEntrance::Int64, FringeQuadExit::Int64, 
        T1::Array{T,1}, T2::Array{T,1}, R1::Array{T,2}, R2::Array{T,2}, 
        RApertures::Array{Float64,1}, EApertures::Array{Float64,1}, KickAngle::Array{T,1}, eletype::String) where T = new{T}(name, len, angle, 
        e1, e2, PolynomA, PolynomB, MaxOrder, NumIntSteps, rad, gK, FringeBendEntrance, FringeBendExit, 
        FringeQuadEntrance, FringeQuadExit, T1, T2, R1, R2, RApertures, EApertures, KickAngle, eletype)
end
function ESBEND(;name::String = "ESBend", len = 0.0, angle = 0.0, e1 = 0.0, e2 = 0.0, 
                PolynomA = zeros(4), PolynomB = zeros(4), 
                MaxOrder::Int64=0, NumIntSteps::Int64 = 10, rad::Int64=0,
                gK = 0.0, FringeBendEntrance::Int64 = 1, FringeBendExit::Int64 = 1,
                FringeQuadEntrance::Int64 = 0, FringeQuadExit::Int64 = 0, T1 = zeros(6),
                T2 = zeros(6), R1 = zeros(6, 6), R2 = zeros(6, 6), RApertures = zeros(6),
                EApertures = zeros(6), KickAngle = zeros(2))
    if PolynomB[2] != 0.0
        MaxOrder = 1
    end
    if PolynomB[3] != 0.0
        MaxOrder = 2
    end
    if PolynomB[4] != 0.0
        MaxOrder = 3
    end
    if len isa DTPSAD || angle isa DTPSAD || e1 isa DTPSAD || e2 isa DTPSAD || PolynomA[1] isa DTPSAD || PolynomB[1] isa DTPSAD ||
        gK isa DTPSAD || T1[1] isa DTPSAD || T2[1] isa DTPSAD || R1[1,1] isa DTPSAD || R2[1,1] isa DTPSAD || KickAngle[1] isa DTPSAD
        return ESBEND(name, DTPSAD(len), DTPSAD(angle), DTPSAD(e1), DTPSAD(e2), 
            DTPSAD.(PolynomA), DTPSAD.(PolynomB), MaxOrder, NumIntSteps, rad, 
            DTPSAD(gK), FringeBendEntrance, FringeBendExit, FringeQuadEntrance, FringeQuadExit,
            DTPSAD.(T1), DTPSAD.(T2), DTPSAD.(R1), DTPSAD.(R2), RApertures,
            EApertures,
            DTPSAD.(KickAngle),
            "ESBEND")
    end
    return ESBEND(name, len, angle, e1, e2, 
        PolynomA, PolynomB, MaxOrder, NumIntSteps, rad, gK, FringeBendEntrance, FringeBendExit, 
        FringeQuadEntrance, FringeQuadExit, T1, T2, R1, R2, RApertures,
        EApertures,
        KickAngle,
        "ESBEND")
end

function ERBEND(;name::String = "ESBend", len = 0.0, angle = 0.0, 
    PolynomA = zeros(4), PolynomB = zeros(4), 
    MaxOrder::Int64=0, NumIntSteps::Int64 = 10, rad::Int64=0,
    gK = 0.0, FringeBendEntrance::Int64 = 1, FringeBendExit::Int64 = 1, 
    FringeQuadEntrance::Int64 = 0, FringeQuadExit::Int64 = 0, T1 = zeros(6), 
    T2 = zeros(6), R1 = zeros(6, 6), R2 = zeros(6, 6), RApertures = zeros(6), 
    EApertures = zeros(6), KickAngle = zeros(2))

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
    if len isa DTPSAD || angle isa DTPSAD || e1 isa DTPSAD || e2 isa DTPSAD || PolynomA[1] isa DTPSAD || PolynomB[1] isa DTPSAD ||
        gK isa DTPSAD || T1[1] isa DTPSAD || T2[1] isa DTPSAD || R1[1,1] isa DTPSAD || R2[1,1] isa DTPSAD || KickAngle[1] isa DTPSAD
        return ESBEND(name, DTPSAD(len), DTPSAD(angle), DTPSAD(e1), DTPSAD(e2), 
            DTPSAD.(PolynomA), DTPSAD.(PolynomB), MaxOrder, NumIntSteps, rad, 
            DTPSAD(gK), FringeBendEntrance, FringeBendExit, FringeQuadEntrance, FringeQuadExit,
            DTPSAD.(T1), DTPSAD.(T2), DTPSAD.(R1), DTPSAD.(R2), RApertures,
            EApertures,
            DTPSAD.(KickAngle),
            "ESBEND")
    end
    return ESBEND(name, len, angle, e1, e2, 
        PolynomA, PolynomB, MaxOrder, NumIntSteps, rad, gK, FringeBendEntrance, FringeBendExit, 
        FringeQuadEntrance, FringeQuadExit, T1, T2, R1, R2, RApertures,
        EApertures,
        KickAngle,
        "ESBEND")
end

"""
    RFCA(;name::String = "RFCA", len::Float64 = 0.0, volt::Float64 = 0.0, freq::Float64 = 0.0, h::Float64 = 1.0, 
        lag::Float64 = 0.0, philag::Float64 = 0.0, energy::Float64 = 0.0)

A RF cavity element.
Example:
```julia
rf = RFCA(name="RF1", len=0.5, volt=1e6, freq=1e6)
```
"""
mutable struct RFCA{T} <: AbstractElement{T}
    name::String
    len::T
    volt::T
    freq::T
    h::T
    lag::T
    philag::T
    energy::T
    eletype::String
    # constructor with all parameters
    RFCA(name::String, len::T, volt::T, freq::T, h::T, lag::T, philag::T, energy::T, eletype::String) where T = 
        new{T}(name, len, volt, freq, h, lag, philag, energy, eletype)
end
function RFCA(;name::String = "RFCA", len = 0.0, volt = 0.0, freq = 0.0, h = 1.0, 
                lag = 0.0, philag = 0.0, energy = 0.0)
    if len isa DTPSAD || volt isa DTPSAD || freq isa DTPSAD || h isa DTPSAD || lag isa DTPSAD || philag isa DTPSAD || energy isa DTPSAD
        return RFCA(name, DTPSAD(len), DTPSAD(volt), DTPSAD(freq), DTPSAD(h), 
                    DTPSAD(lag), DTPSAD(philag), DTPSAD(energy), "RFCA")
    end
    return RFCA(name, len, volt, freq, h, lag, philag, energy, "RFCA")
end

"""
    SOLENOID(;name::String = "Solenoid", len::Float64 = 0.0, ks::Float64 = 0.0, T1::Array{Float64,1} = zeros(6), 
        T2::Array{Float64,1} = zeros(6), R1::Array{Float64,2} = zeros(6,6), R2::Array{Float64,2} = zeros(6,6))

A solenoid element.
Example:
```julia
solenoid = SOLENOID(name="S1", len=0.5, ks=1.0)
```
"""
mutable struct SOLENOID{T} <: AbstractElement{T}
    name::String
    len::T
    ks::T # rad/m
    T1::Array{T,1}
    T2::Array{T,1}
    R1::Array{T,2}
    R2::Array{T,2}
    eletype::String
    # constructor with all parameters
    SOLENOID(name::String, len::T, ks::T, T1::Array{T,1}, T2::Array{T,1}, 
    R1::Array{T,2}, R2::Array{T,2}, eletype::String) where T = new{T}(name, len, ks, T1, T2, R1, R2, eletype)
end
function SOLENOID(;name::String = "Solenoid", len = 0.0, ks = 0.0, T1 = zeros(6), 
                T2 = zeros(6), R1 = zeros(6,6), R2 = zeros(6,6))
    if len isa DTPSAD || ks isa DTPSAD || T1[1] isa DTPSAD || T2[1] isa DTPSAD || R1[1,1] isa DTPSAD || R2[1,1] isa DTPSAD
        return SOLENOID(name, DTPSAD(len), DTPSAD(ks), DTPSAD.(T1), DTPSAD.(T2), DTPSAD.(R1), DTPSAD.(R2), "SOLENOID")
    end
    return SOLENOID(name, len, ks, T1, T2, R1, R2, "SOLENOID")
end

"""
    CORRECTOR(;name::String = "CORRECTOR", len::Float64 = 0.0, xkick::Float64 = 0.0, ykick::Float64 = 0.0, 
        T1::Array{Float64,1} = zeros(6), T2::Array{Float64,1} = zeros(6), R1::Array{Float64,2} = zeros(6,6), 
        R2::Array{Float64,2} = zeros(6,6))

A corrector element.
Example:
```julia
corrector = CORRECTOR(name="C1", len=0.5, xkick=1e-3)
```
"""
mutable struct CORRECTOR{T} <: AbstractElement{T}
    name::String
    len::T
    xkick::T
    ykick::T
    T1::Array{T,1}
    T2::Array{T,1}
    R1::Array{T,2}
    R2::Array{T,2}
    eletype::String

    # constructor with all parameters
    CORRECTOR(name::String, len::T, xkick::T, ykick::T, 
    T1::Array{T,1}, T2::Array{T,1}, R1::Array{T,2}, R2::Array{T,2}, eletype::String) where T = 
        new{T}(name, len, xkick, ykick, T1, T2, R1, R2, eletype)
end
function CORRECTOR(;name::String = "CORRECTOR", len = 0.0, xkick = 0.0, ykick = 0.0, 
                    T1 = zeros(6), T2 = zeros(6), R1 = zeros(6,6), 
                    R2 = zeros(6,6))
    if len isa DTPSAD || xkick isa DTPSAD || ykick isa DTPSAD || T1[1] isa DTPSAD || T2[1] isa DTPSAD || R1[1,1] isa DTPSAD || R2[1,1] isa DTPSAD
        return CORRECTOR(name, DTPSAD(len), DTPSAD(xkick), DTPSAD(ykick), 
                        DTPSAD.(T1), DTPSAD.(T2), DTPSAD.(R1), DTPSAD.(R2), "CORRECTOR")
    end
    return CORRECTOR(name, len, xkick, ykick, T1, T2, R1, R2, "CORRECTOR")
end

function HKICKER(;name::String = "HKicker", len = 0.0, xkick = 0.0)
    return CORRECTOR(name=name, len=len, xkick=xkick, ykick=0.0)
end
function VKICKER(;name::String = "VKicker", len = 0.0, ykick = 0.0)
    return CORRECTOR(name=name, len=len, xkick=0.0, ykick=ykick)
end

mutable struct SPACECHARGE{T} <: AbstractElement{T}
    # spectral space charge
    # this element is treated as an integrated effect of space charge over a length of effective_len
    name::String
    len::T
    effective_len::T
    Nl::Int64
    Nm::Int64
    a::T
    b::T
    eletype::String

    # constructor with all parameters
    SPACECHARGE(name::String, len::T, effective_len::T, Nl::Int64, Nm::Int64, a::T, b::T, eletype::String) where T = 
        new{T}(name, len, effective_len, Nl, Nm, a, b, eletype)
end
function SPACECHARGE(;name::String = "SPACECHARGE", len = 0.0, effective_len = 0.0, Nl::Int64 = 10, Nm::Int64 = 10, a = 1.0, b = 1.0)
    if len isa DTPSAD || effective_len isa DTPSAD || a isa DTPSAD || b isa DTPSAD
        return SPACECHARGE(name, DTPSAD(len), DTPSAD(effective_len), Nl, Nm, DTPSAD(a), DTPSAD(b), "SPACECHARGE")
    end
    return SPACECHARGE(name, len, effective_len, Nl, Nm, a, b, "SPACECHARGE")
end

# TRANSLATION and YROTATION are used to convert the MAD-X lattice files
mutable struct TRANSLATION{T} <: AbstractElement{T}
    name::String
    len::T
    dx::T
    dy::T
    ds::T
    eletype::String
    # constructor with all parameters
    TRANSLATION(name::String, len::T, dx::T, dy::T, ds::T, eletype::String) where T = 
        new{T}(name, len, dx, dy, ds, eletype)
end
function TRANSLATION(;name::String = "TRANSLATION", len = 0.0, dx = 0.0, dy = 0.0, ds = 0.0)
    if len isa DTPSAD || dx isa DTPSAD || dy isa DTPSAD || ds isa DTPSAD
        return TRANSLATION(name, DTPSAD(len), DTPSAD(dx), DTPSAD(dy), DTPSAD(ds), "TRANSLATION")
    end
    return TRANSLATION(name, len, dx, dy, ds, "TRANSLATION")
end

mutable struct YROTATION{T} <: AbstractElement{T}
    name::String
    len::T
    angle::T
    eletype::String

    # constructor with all parameters
    YROTATION(name::String, len::T, angle::T, eletype::String) where T = 
        new{T}(name, len, angle, eletype)
end
function YROTATION(;name::String = "YROTATION", len = 0.0, angle = 0.0)
    if len isa DTPSAD || angle isa DTPSAD
        return YROTATION(name, DTPSAD(len), DTPSAD(angle), "YROTATION")
    end
    return YROTATION(name, len, angle, "YROTATION")
end

# non-canonical elements
"""
    QUAD(;name::String = "Quad", len::Float64 = 0.0, k1::Float64 = 0.0, rad::Int64 = 0, 
        T1::Array{Float64,1} = zeros(6), T2::Array{Float64,1} = zeros(6), R1::Array{Float64,2} = zeros(6,6), 
        R2::Array{Float64,2} = zeros(6,6), RApertures::Array{Float64,1} = zeros(6), EApertures::Array{Float64,1} = zeros(6))

A quadrupole element using matrix formalism.
Example:
```julia
quad = QUAD(name="Q1", len=0.5, k1=1.0)
```
"""
mutable struct QUAD{T} <: AbstractElement{T}
    name::String
    len::T
    k1::T
    rad::Int64
    T1::Array{T,1}
    T2::Array{T,1}
    R1::Array{T,2}
    R2::Array{T,2}
    RApertures::Array{Float64,1}
    EApertures::Array{Float64,1}
    eletype::String
    # constructor with all parameters
    QUAD(name::String, len::T, k1::T, rad::Int64, 
    T1::Array{T,1}, T2::Array{T,1}, R1::Array{T,2}, R2::Array{T,2}, 
    RApertures::Array{Float64,1}, EApertures::Array{Float64,1}, eletype::String) where T = 
        new{T}(name, len, k1, rad, T1, T2, R1, R2, RApertures, EApertures, eletype)
end
function QUAD(;name::String = "Quad", len=0.0, k1 = 0.0, rad::Int64 = 0, 
                T1 = zeros(6), T2 = zeros(6), 
                R1 = zeros(6,6), R2 = zeros(6,6), 
                RApertures = zeros(6), EApertures = zeros(6))
    if len isa DTPSAD || k1 isa DTPSAD || T1[1] isa DTPSAD || T2[1] isa DTPSAD || R1[1,1] isa DTPSAD || R2[1,1] isa DTPSAD
        return QUAD(name, DTPSAD(len), DTPSAD(k1), rad, 
                    DTPSAD.(T1), DTPSAD.(T2), DTPSAD.(R1), DTPSAD.(R2), 
                    RApertures, EApertures, "QUAD")
    end
    return QUAD(name, len, k1, rad, T1, T2, R1, R2, RApertures, EApertures, "QUAD")
end

"""
    QUAD_SC(;name::String = "Quad", len::Float64 = 0.0, k1::Float64 = 0.0, rad::Int64 = 0, 
        T1::Array{Float64,1} = zeros(6), T2::Array{Float64,1} = zeros(6), R1::Array{Float64,2} = zeros(6,6), 
        R2::Array{Float64,2} = zeros(6,6), RApertures::Array{Float64,1} = zeros(6), EApertures::Array{Float64,1} = zeros(6), 
        a::Float64 = 1.0, b::Float64 = 1.0, Nl::Int64 = 10, Nm::Int64 = 10, Nsteps::Int64=1)

A quadrupole element with space charge.
Example:
```julia
quad = QUAD_SC(name="Q1_SC", len=0.5, k1=1.0, a=13e-3, b=13e-3, Nl=15, Nm=15)
```
"""
mutable struct QUAD_SC{T} <: AbstractElement{T}
    name::String
    len::T
    k1::T
    rad::Int64
    T1::Array{T,1}
    T2::Array{T,1}
    R1::Array{T,2}
    R2::Array{T,2}
    RApertures::Array{Float64,1}
    EApertures::Array{Float64,1}
    a::T
    b::T
    Nl::Int64
    Nm::Int64
    Nsteps::Int64
    eletype::String
    
    # constructor with all parameters
    QUAD_SC(name::String, len::T, k1::T, rad::Int64,
    T1::Array{T,1}, T2::Array{T,1}, R1::Array{T,2}, R2::Array{T,2},
    RApertures::Array{Float64,1}, EApertures::Array{Float64,1}, a::T, b::T,
    Nl::Int64, Nm::Int64, Nsteps::Int64, eletype::String) where T =
        new{T}(name, len, k1, rad, T1, T2, R1, R2, RApertures, EApertures, a, b, Nl, Nm, Nsteps, eletype)
end
function QUAD_SC(;name::String = "Quad", len = 0.0, k1 = 0.0, rad::Int64 = 0, 
                T1 = zeros(6), T2 = zeros(6), 
                R1 = zeros(6,6), R2 = zeros(6,6), 
                RApertures = zeros(6), EApertures = zeros(6), 
                a = 1.0, b = 1.0, Nl::Int64 = 10, Nm::Int64 = 10, Nsteps::Int64=1)
    if len isa DTPSAD || k1 isa DTPSAD || T1[1] isa DTPSAD || T2[1] isa DTPSAD || R1[1,1] isa DTPSAD || R2[1,1] isa DTPSAD
        return QUAD_SC(name, DTPSAD(len), DTPSAD(k1), rad,
                    DTPSAD.(T1), DTPSAD.(T2), DTPSAD.(R1), DTPSAD.(R2),
                    RApertures, EApertures, DTPSAD(a), DTPSAD(b),
                    Nl, Nm, Nsteps, "QUAD_SC")
    end
    return QUAD_SC(name, len, k1, rad, T1, T2, R1, R2, RApertures, EApertures, a, b,
                    Nl, Nm, Nsteps, "QUAD_SC")
end


"""
    WIGGLER(;name::String = "WIGGLER", len::Float64 = 0.0, lw::Float64 = 0.0, Bmax::Float64 = 0.0, 
            Nsteps::Int64 = 10, By::Array{Int} = [1;1;0;1;1;0], Bx::Array{Int} = Int[], energy::Float64 = 1e9,
            R1::Array{Float64,2} = zeros(6,6), R2::Array{Float64,2} = zeros(6,6), 
            T1::Array{Float64,1} = zeros(6), T2::Array{Float64,1} = zeros(6))
A wiggler element.
- arameters:
    - len: total length of the wiggler (m)
    - lw: period length of the wiggler (m)
    - Bmax: Peak magnetic field (T)
    - Nsteps: number of integration steps
    - By: wiggler harmonics for horizontal wigglers. Default [1;1;0;1;1;0]
    - Bx: wiggler harmonics for vertical wigglers. Default []
    - energy: reference energy (eV)
Example:
```julia
wiggler = WIGGLER(name="W1", len=1.0, lw=0.1, Bmax=1.0)
```
"""
mutable struct WIGGLER{T} <: AbstractElement{T}
    name::String 
    len::T 
    lw::T 
    Bmax::T 
    Nsteps::Int64 
    By::Array{Int,1}
    Bx::Array{Int,1}
    energy::T
    NHharm::Int64
    NVharm::Int64
    rad::Int64
    R1::Array{T,2}
    R2::Array{T,2}
    T1::Array{T,1}
    T2::Array{T,1}
    eletype::String 

    # constructor with all parameters
    WIGGLER(name::String, len::T, lw::T, Bmax::T, Nsteps::Int64, By::Array{Int,1}, Bx::Array{Int,1},
    energy::T, NHharm::Int64, NVharm::Int64, rad::Int64, R1::Array{T,2}, R2::Array{T,2}, T1::Array{T,1}, T2::Array{T,1}, eletype::String) where T = 
        new{T}(name, len, lw, Bmax, Nsteps, By, Bx, energy, NHharm, NVharm, rad, R1, R2, T1, T2, eletype)
end
function WIGGLER(;name::String = "WIGGLER", len = 0.0, lw = 0.0, Bmax = 0.0, 
    Nsteps::Int64 = 10, By = [1;1;0;1;1;0], Bx = Int[], energy = 1e9, rad::Int64 = 0,
    R1 = zeros(6,6), R2 = zeros(6,6), T1 = zeros(6), T2 = zeros(6))

    NHharm = length(By) ÷ 6
    NVharm = length(Bx) ÷ 6
    if len isa DTPSAD || lw isa DTPSAD || Bmax isa DTPSAD || By isa DTPSAD || Bx isa DTPSAD || energy isa DTPSAD ||
        R1[1,1] isa DTPSAD || R2[1,1] isa DTPSAD || T1[1] isa DTPSAD || T2[1] isa DTPSAD
        return WIGGLER(name, DTPSAD(len), DTPSAD(lw), DTPSAD(Bmax), Nsteps, By, Bx,
            DTPSAD(energy), NHharm, NVharm, rad, DTPSAD.(R1), DTPSAD.(R2), DTPSAD.(T1), DTPSAD.(T2), "WIGGLER")
    end
    return WIGGLER(name, len, lw, Bmax, Nsteps, By, Bx,
        energy, NHharm, NVharm, rad, R1, R2, T1, T2, "WIGGLER")
end

###########################################
"""
    CRABCAVITY(;name::String = "CRABCAVITY", len::Float64 = 0.0, volt::Float64 = 0.0, freq::Float64 = 0.0, 
        phi::Float64 = 0.0, errors::Array{Float64,1} = zeros(2), energy::Float64 = 1e9)

A crab cavity element.
Example:
```julia
crab = CRABCAVITY(name="CRAB1", len=0.5, volt=1e6, freq=1e6)
```
"""
mutable struct CRABCAVITY{T} <: AbstractElement{T}
    name::String 
    len::T 
    volt::T  # voltage
    freq::T  # frequency
    k::T  # wave number
    phi::T  # phase
    errors::Array{T,1} # 1: Voltage error, 2: Phase error
    energy::T
    eletype::String 

    # constructor with all parameters
    CRABCAVITY(name::String, len::T, volt::T, freq::T, k::T, phi::T,
    errors::Array{T,1}, energy::T, eletype::String) where T = 
        new{T}(name, len, volt, freq, k, phi, errors, energy, eletype)
end
function CRABCAVITY(;name::String = "CRABCAVITY", len = 0.0, volt = 0.0, 
    freq = 0.0, phi = 0.0, errors = zeros(2), energy = 1e9)
    k = 2*π*freq/2.99792458e8
    if len isa DTPSAD || volt isa DTPSAD || freq isa DTPSAD || phi isa DTPSAD || errors[1] isa DTPSAD || energy isa DTPSAD
        return CRABCAVITY(name, DTPSAD(len), DTPSAD(volt), DTPSAD(freq), 
            DTPSAD(k), DTPSAD(phi), DTPSAD.(errors), DTPSAD(energy), "CRABCAVITY")
    end
    return CRABCAVITY(name, len, volt, freq, k, phi, errors, energy, "CRABCAVITY")
end

mutable struct CRABCAVITY_K2{T} <: AbstractElement{T}
    name::String 
    len::T 
    volt::T  # voltage
    freq::T  # frequency
    k::T  # wave number
    phi::T  # phase
    k2::T  # k2
    errors::Array{T,1} # 1: Voltage error, 2: Phase error
    energy::T
    eletype::String 

    # constructor with all parameters
    CRABCAVITY_K2(name::String, len::T, volt::T, freq::T, k::T, phi::T,
    k2::T, errors::Array{T,1}, energy::T, eletype::String) where T = 
        new{T}(name, len, volt, freq, k, phi, k2, errors, energy, eletype)
end
function CRABCAVITY_K2(;name::String = "CRABCAVITY_K2", len = 0.0, volt = 0.0, 
    freq = 0.0, phi = 0.0, k2 = 0.0, errors = zeros(2), energy = 1e9)
    k = 2*π*freq/2.99792458e8
    if len isa DTPSAD || volt isa DTPSAD || freq isa DTPSAD || phi isa DTPSAD || k2 isa DTPSAD || errors[1] isa DTPSAD || energy isa DTPSAD
        return CRABCAVITY_K2(name, DTPSAD(len), DTPSAD(volt), DTPSAD(freq), 
            DTPSAD(k), DTPSAD(phi), DTPSAD(k2), DTPSAD.(errors), DTPSAD(energy), "CRABCAVITY_K2")
    end
    return CRABCAVITY_K2(name, len, volt, freq, k, phi, k2, errors, energy, "CRABCAVITY_K2")
end


mutable struct easyCRABCAVITY{T} <: AbstractElement{T}
    name::String 
    len::T 
    halfthetac::T 
    freq::T 
    k::T 
    phi::T 
    errors::Array{T,1}  # 1: Voltage error, 2: Phase error
    eletype::String

    # constructor with all parameters
    easyCRABCAVITY(name::String, len::T, halfthetac::T, freq::T,
    k::T, phi::T, errors::Array{T,1}, eletype::String) where T =
        new{T}(name, len, halfthetac, freq, k, phi, errors, eletype)
end
function easyCRABCAVITY(;name::String = "easyCRABCAVITY", len = 0.0, halfthetac = 0.0, 
    freq = 0.0, k = 0.0, phi = 0.0, errors = zeros(2), eletype::String = "easyCRABCAVITY")
    if len isa DTPSAD || halfthetac isa DTPSAD || freq isa DTPSAD || k isa DTPSAD || phi isa DTPSAD || errors[1] isa DTPSAD
        return easyCRABCAVITY(name, DTPSAD(len), DTPSAD(halfthetac), DTPSAD(freq), 
            DTPSAD(k), DTPSAD(phi), DTPSAD.(errors), eletype)
    end
    return easyCRABCAVITY(name, len, halfthetac, freq, k, phi, errors, eletype)
end

mutable struct AccelCavity{T} <: AbstractElement{T}
    name::String 
    len::T 
    volt::T  # voltage
    freq::T  # frequency
    k::T  # wave number
    h::T  # harmonic number
    phis::T  # synchronous phase π/2 for accelerating on crest
    eletype::String 

    # constructor with all parameters
    AccelCavity(name::String, len::T, volt::T, freq::T, k::T, h::T, phis::T, eletype::String) where T = 
        new{T}(name, len, volt, freq, k, h, phis, eletype)
end
function AccelCavity(;name::String = "AccelCavity", len = 0.0, volt = 0.0, 
    freq = 0.0, h = 1.0, phis = 0.0)
    k = 2*π*freq/2.99792458e8
    if len isa DTPSAD || volt isa DTPSAD || freq isa DTPSAD || h isa DTPSAD || phis isa DTPSAD
        return AccelCavity(name, DTPSAD(len), DTPSAD(volt), DTPSAD(freq), 
            DTPSAD(k), DTPSAD(h), DTPSAD(phis), "AccelCavity")
    end
    return AccelCavity(name, len, volt, freq, k, h, phis, "AccelCavity")
end


abstract type AbstractTransferMap{T} <: AbstractElement{T} end
abstract type AbstractTransverseMap{T} <: AbstractTransferMap{T} end
abstract type AbstractLongitudinalMap{T} <:AbstractTransferMap{T} end
mutable struct LongitudinalRFMap{T} <: AbstractLongitudinalMap{T}
    alphac::T
    RF::AbstractElement{T}
    LongitudinalRFMap(alphac::T, RF::AbstractElement{T}) where T = new{T}(alphac, RF)
end

mutable struct LorentzBoost{T} <: AbstractElement{T}
    angle::T
    cosang::T
    tanang::T
    mode::Int
    LorentzBoost(angle::T) where T = new{T}(angle, cos(angle), tan(angle), 0)
end

mutable struct InvLorentzBoost{T} <: AbstractElement{T}
    angle::T
    sinang::T
    cosang::T
    mode::Int
    InvLorentzBoost(angle::T) where T = new{T}(angle, sin(angle), cos(angle), 0)
end


# optics
abstract type AbstractOptics end
abstract type AbstractOptics2D <:AbstractOptics end
abstract type AbstractOptics4D <:AbstractOptics end
mutable struct optics2D <: AbstractOptics2D
    beta::Float64
    alpha::Float64
    gamma::Float64
    phase::Float64
    eta::Float64
    etap::Float64
end
optics2D(b::Float64, a::Float64, ph::Float64, eta::Float64, etap::Float64)=optics2D(b,a,(1.0+a^2)/b,ph,eta,etap)
optics2D(b::Float64, a::Float64)=optics2D(b,a,(1.0+a^2)/b,0.0,0.0,0.0)

mutable struct optics4DUC <: AbstractOptics4D # 4D linear transformation uncoupled
    optics_x::optics2D
    optics_y::optics2D
end
"""
    optics4DUC(bx::Float64, ax::Float64, by::Float64, ay::Float64)

Construct a 4D optics element with uncoupled optics.

# Arguments
- `bx::Float64`: beta function in x direction
- `ax::Float64`: alpha function in x direction
- `by::Float64`: beta function in y direction
- `ay::Float64`: alpha function in y direction

# Returns
- `optics4DUC`: 4D optics element with uncoupled optics
"""
optics4DUC(bx::Float64, ax::Float64, by::Float64, ay::Float64)=optics4DUC(optics2D(bx,ax),optics2D(by,ay))

######### strong beam-beam
abstract type AbstractStrongBeamBeam{T} <:AbstractElement{T} end

mutable struct StrongThinGaussianBeam{T} <: AbstractStrongBeamBeam{T}
    amplitude::T
    rmssizex::T
    rmssizey::T
    zloc::T
    xoffset::T
    yoffset::T
    
    function StrongThinGaussianBeam(amp::T, rx::T, ry::T, zloc::T=zero(T), 
                                    xoff::T=zero(T), yoff::T=zero(T)) where {T}
        new{T}(amp, rx, ry, zloc, xoff, yoff)
    end
end

"""
    StrongGaussianBeam(charge::Float64, mass::Float64, atomnum::Float64, np::Int, energy::Float64, 
    op::AbstractOptics4D, bs::Vector{Float64}, nz::Int)

Construct a strong beam-beam element with Gaussian distribution.
# Arguments
- `charge::Float64`: charge of the particle
- `mass::Float64`: mass of the particle
- `atomnum::Float64`: atomic number of the particle
- `np::Int`: number of particles in the beam
- `energy::Float64`: total energy of the beam
- `op::AbstractOptics4D`: optics at the interaction point
- `bs::Vector{Float64}`: beam size at the interaction point
- `nz::Int`: number of slices in z direction
# Returns
- `StrongGaussianBeam`: strong beam-beam element with Gaussian distribution
"""
mutable struct StrongGaussianBeam{T} <: AbstractStrongBeamBeam{T}  # Strong Beam with transverse Gaussian distribution
    # particle::ParticleType
    charge::T
    mass::T
    atomnum::T
    num_particle::Int  # Number of particles
    total_energy::T # Total energy of the beam
    momentum::T  # Design Momentum of the beam
    gamma::T  # Relativistic gamma
    beta::T  # Relativistic beta v/c
    optics::AbstractOptics4D # optics @IP
    beamsize::Vector{T} # Beam size at IP
    nzslice::Int # Number of slices in z direction
    zslice_center::Vector{T} # z center of each slice
    zslice_npar::Vector{T} # amplitude of each slice
    xoffsets::Vector{T} # x offset of each slice
    yoffsets::Vector{T} # y offset of each slice
    # constructor with all parameters
    StrongGaussianBeam(charge::T, mass::T, atomnum::T, 
            num_particle::Int, total_energy::T, momentum::T, gamma::T, beta::T, 
            optics::AbstractOptics4D, beamsize::Vector{T}, nzslice::Int, zslice_center::Vector{T}, 
            zslice_npar::Vector{T}, xoffsets::Vector{T}, yoffsets::Vector{T}) where {T} = 
        new{T}(charge, mass, atomnum, num_particle, total_energy, momentum, gamma, beta,
            optics, beamsize, nzslice, zslice_center, zslice_npar, xoffsets, yoffsets)
end
function StrongGaussianBeam(charge::Float64, mass::Float64, atomum::Float64, np::Int, 
        energy::Float64, op::AbstractOptics4D, bs::Vector{Float64}, nz::Int)
    momentum=sqrt(energy*energy-mass*mass)
    gamma = energy / mass
    beta = momentum / energy
    return StrongGaussianBeam(charge, mass, atomum, np, energy, 
        momentum, gamma, beta, op, bs, nz, 
        zeros(Float64, nz), zeros(Float64, nz), zeros(Float64, nz), zeros(Float64, nz))
end
function StrongGaussianBeam(charge::DTPSAD{N, T}, mass::DTPSAD{N, T}, atomum::DTPSAD{N, T}, np::Int, 
        energy::DTPSAD{N, T}, op::AbstractOptics4D, bs::Vector{DTPSAD{N, T}}, nz::Int) where {N, T}
    momentum=sqrt(energy*energy-mass*mass)
    gamma = energy / mass
    beta = momentum / energy
    return StrongGaussianBeam(charge, mass, atomum, np, energy, 
        momentum, gamma, beta, op, bs, nz,
        zeros(DTPSAD{N, T}, nz), zeros(DTPSAD{N, T}, nz), zeros(DTPSAD{N, T}, nz), zeros(DTPSAD{N, T}, nz))
end 

# wake field
function linear_interpolate(x, x_points, y_points)
    """
    Interpolates or extrapolates a value using linear interpolation.

    x: The point to interpolate or extrapolate.
    x_points: The x-coordinates of the data points.
    y_points: The y-coordinates of the data points.

    Returns the interpolated or extrapolated value at x.
    """
    if x <= x_points[1]
        slope = (y_points[2] - y_points[1]) / (x_points[2] - x_points[1])
        return y_points[1] + slope * (x - x_points[1])
    elseif x >= x_points[end]
        slope = (y_points[end] - y_points[end - 1]) / (x_points[end] - x_points[end - 1])
        return y_points[end] + slope * (x - x_points[end])
    else
        for i in 2:length(x_points)
            if x < x_points[i]
                slope = (y_points[i] - y_points[i - 1]) / (x_points[i] - x_points[i - 1])
                return y_points[i - 1] + slope * (x - x_points[i - 1])
            end
        end
    end
end

"""
    LongitudinalRLCWake(;freq::Float64=1.0e9, Rshunt::Float64=1.0e6, Q0::Float64=1.0)

A longitudinal RLC wake element.
"""
mutable struct LongitudinalRLCWake{T} <: AbstractElement{T}
    name::String
    freq::T
    Rshunt::T
    Q0::T
    # constructor with all parameters
    LongitudinalRLCWake(name::String, freq::T, Rshunt::T, Q0::T) where T = new{T}(name, freq, Rshunt, Q0)
end
function LongitudinalRLCWake(;name::String="RLCWake", freq::T=1.0e9, Rshunt::T=1.0e6, Q0::T=1.0) where {T}
    if freq isa DTPSAD || Rshunt isa DTPSAD || Q0 isa DTPSAD
        return LongitudinalRLCWake(name, DTPSAD(freq), DTPSAD(Rshunt), DTPSAD(Q0))
    end
    return LongitudinalRLCWake(name, freq, Rshunt, Q0)
end

function wakefieldfunc_RLCWake(rlcwake::LongitudinalRLCWake, t::T) where {T}
    Q0p=sqrt(rlcwake.Q0^2 - 1.0/4.0)
    w0 = 2*pi*rlcwake.freq
    w0p= w0/rlcwake.Q0*Q0p
    if t>0
        return zero(T)
    else
        return rlcwake.Rshunt * w0 /rlcwake.Q0 * (cos(w0p * t) +  sin(w0p * t) / 2.0 / Q0p) * exp(w0 * t / 2.0 / rlcwake.Q0)
    end
end

mutable struct LongitudinalWake{T} <: AbstractElement{T}
    name::String
    times::AbstractVector{T}
    wakefields::AbstractVector{T}
    wakefield::Function
    # constructor with all parameters
    LongitudinalWake(name::String, times::AbstractVector{T}, wakefields::AbstractVector{T}, wakefield::Function) where {T} =
        new{T}(name, times, wakefields, wakefield)
end
"""
    LongitudinalWake(times::AbstractVector, wakefields::AbstractVector, wakefield::Function)

Create longitudinal wake element.
# Arguments
- `times::AbstractVector`: time points
- `wakefields::AbstractVector`: wakefield values at the time points
- `fliphalf::Float64=-1.0`: flip the wakefield for positrons if needed
# Returns
- `LongitudinalWake`: the created longitudinal wake element
"""
function LongitudinalWake(times::AbstractVector, wakefields::AbstractVector, fliphalf::Float64=-1.0; name::String="LongitudinalWake")
    wakefield = function (t::Float64)
        t>times[1]*fliphalf && return 0.0
        return linear_interpolate(t*fliphalf, times, wakefields)
    end
    return LongitudinalWake(name, times, wakefields, wakefield)
end