"""
    AbstractElement

Abstract type for all elements in the lattice.
"""
abstract type AbstractElement end

"""
    MARKER(;name::String = "MARKER", len::Float64 = 0.0)

A marker element.
Example:
```julia
marker = MARKER(name="MARKER1")
```
"""
mutable struct MARKER <: AbstractElement
    name::String
    len::Float64
    eletype::String

    MARKER(;name::String = "MARKER", len::Float64 = 0.0) = new(name, len, "MARKER")
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
mutable struct DRIFT <: AbstractElement
    name::String
    len::Float64
    T1::Array{Float64,1}
    T2::Array{Float64,1}
    R1::Array{Float64,2}
    R2::Array{Float64,2}        
    RApertures::Array{Float64,1}
    EApertures::Array{Float64,1}
    eletype::String

    DRIFT(;name::String = "DRIFT", len::Float64 = 0.0, T1::Array{Float64,1} = zeros(6), 
        T2::Array{Float64,1} = zeros(6), R1::Array{Float64,2} = zeros(6,6), R2::Array{Float64,2} = zeros(6,6), 
        RApertures::Array{Float64,1} = zeros(6), EApertures::Array{Float64,1} = zeros(6)) = new(name, len, T1, T2, R1, R2, RApertures, EApertures, "DRIFT")
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
mutable struct DRIFT_SC <: AbstractElement
    name::String
    len::Float64
    T1::Array{Float64,1}
    T2::Array{Float64,1}
    R1::Array{Float64,2}
    R2::Array{Float64,2}        
    RApertures::Array{Float64,1}
    EApertures::Array{Float64,1}
    a::Float64
    b::Float64
    Nl::Int64
    Nm::Int64
    Nsteps::Int64 # Number of steps for space charge calculation. One step represents a half-kick-half.
    eletype::String

    DRIFT_SC(;name::String = "DRIFT_SC", len::Float64 = 0.0, T1::Array{Float64,1} = zeros(6), 
        T2::Array{Float64,1} = zeros(6), R1::Array{Float64,2} = zeros(6,6), R2::Array{Float64,2} = zeros(6,6), 
        RApertures::Array{Float64,1} = zeros(6), EApertures::Array{Float64,1} = zeros(6), a::Float64 = 1.0, b::Float64 = 1.0,
        Nl::Int64 = 10, Nm::Int64 = 10, Nsteps::Int64=1) = new(name, len, T1, T2, R1, R2, RApertures, EApertures, a, b, Nl, Nm, Nsteps, "DRIFT_SC")
end

"""
    KQUAD(;name::String = "Quad", len::Float64 = 0.0, k1::Float64 = 0.0, 
        PolynomA::Array{Float64,1} = zeros(Float64, 4), PolynomB::Array{Float64,1} = zeros(Float64, 4), 
        MaxOrder::Int64=1, NumIntSteps::Int64 = 10, rad::Int64=0, FringeQuadEntrance::Int64 = 0, 
        FringeQuadExit::Int64 = 0, FringeIntM0::Array{Float64,1} = zeros(Float64, 5), 
        FringeIntP0::Array{Float64,1} = zeros(Float64, 5), T1::Array{Float64,1} = zeros(Float64, 6), 
        T2::Array{Float64,1} = zeros(Float64, 6), R1::Array{Float64,2} = zeros(Float64, 6, 6), 
        R2::Array{Float64,2} = zeros(Float64, 6, 6), RApertures::Array{Float64,1} = zeros(Float64, 6), 
        EApertures::Array{Float64,1} = zeros(Float64, 6), KickAngle::Array{Float64,1} = zeros(Float64, 2))

A canonical quadrupole element.

Example:
```julia
quad = KQUAD(name="Q1", len=0.5, k1=0.5)
```
"""
mutable struct KQUAD <: AbstractElement
    name::String
    len::Float64
    k1::Float64
    PolynomA::Array{Float64,1}
    PolynomB::Array{Float64,1}
    MaxOrder::Int64
    NumIntSteps::Int64
    rad::Int64
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

    function KQUAD(;name::String = "Quad", len::Float64 = 0.0, k1::Float64 = 0.0, 
                    PolynomA::Array{Float64,1} = zeros(Float64, 4), 
                    PolynomB::Array{Float64,1} = zeros(Float64, 4), MaxOrder::Int64=1, 
                    NumIntSteps::Int64 = 10, rad::Int64=0, FringeQuadEntrance::Int64 = 0, 
                    FringeQuadExit::Int64 = 0, FringeIntM0::Array{Float64,1} = zeros(Float64, 5), 
                    FringeIntP0::Array{Float64,1} = zeros(Float64, 5), T1::Array{Float64,1} = zeros(Float64, 6), 
                    T2::Array{Float64,1} = zeros(Float64, 6), R1::Array{Float64,2} = zeros(Float64, 6, 6), 
                    R2::Array{Float64,2} = zeros(Float64, 6, 6), RApertures::Array{Float64,1} = zeros(Float64, 6), 
                    EApertures::Array{Float64,1} = zeros(Float64, 6), KickAngle::Array{Float64,1} = zeros(Float64, 2))
        if k1 != 0.0 && PolynomB[2] == 0.0
            PolynomB[2] = k1
        end
        new(name, len, k1, PolynomA, PolynomB, MaxOrder, NumIntSteps, rad, FringeQuadEntrance, FringeQuadExit, 
            FringeIntM0, FringeIntP0, T1, T2, R1, R2, RApertures, EApertures, KickAngle, "KQUAD")
    end
end

"""
    KQUAD_SC(;name::String = "Quad", len::Float64 = 0.0, k1::Float64 = 0.0, 
        PolynomA::Array{Float64,1} = zeros(Float64, 4), PolynomB::Array{Float64,1} = zeros(Float64, 4), 
        MaxOrder::Int64=1, NumIntSteps::Int64 = 10, rad::Int64=0, FringeQuadEntrance::Int64 = 0, 
        FringeQuadExit::Int64 = 0, FringeIntM0::Array{Float64,1} = zeros(Float64, 5), 
        FringeIntP0::Array{Float64,1} = zeros(Float64, 5), T1::Array{Float64,1} = zeros(Float64, 6), 
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
mutable struct KQUAD_SC <: AbstractElement
    name::String
    len::Float64
    k1::Float64
    PolynomA::Array{Float64,1}
    PolynomB::Array{Float64,1}
    MaxOrder::Int64
    NumIntSteps::Int64
    rad::Int64
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
    a::Float64
    b::Float64
    Nl::Int64
    Nm::Int64
    Nsteps::Int64
    eletype::String

    function KQUAD_SC(;name::String = "Quad", len::Float64 = 0.0, k1::Float64 = 0.0, 
                    PolynomA::Array{Float64,1} = zeros(Float64, 4), 
                    PolynomB::Array{Float64,1} = zeros(Float64, 4), MaxOrder::Int64=1, 
                    NumIntSteps::Int64 = 10, rad::Int64=0, FringeQuadEntrance::Int64 = 0, 
                    FringeQuadExit::Int64 = 0, FringeIntM0::Array{Float64,1} = zeros(Float64, 5), 
                    FringeIntP0::Array{Float64,1} = zeros(Float64, 5), T1::Array{Float64,1} = zeros(Float64, 6), 
                    T2::Array{Float64,1} = zeros(Float64, 6), R1::Array{Float64,2} = zeros(Float64, 6, 6), 
                    R2::Array{Float64,2} = zeros(Float64, 6, 6), RApertures::Array{Float64,1} = zeros(Float64, 6), 
                    EApertures::Array{Float64,1} = zeros(Float64, 6), KickAngle::Array{Float64,1} = zeros(Float64, 2),
                    a::Float64 = 1.0, b::Float64 = 1.0, Nl::Int64 = 10, Nm::Int64 = 10, Nsteps::Int64=1)
        if k1 != 0.0 && PolynomB[2] == 0.0
            PolynomB[2] = k1
        end
        new(name, len, k1, PolynomA, PolynomB, MaxOrder, NumIntSteps, rad, FringeQuadEntrance, FringeQuadExit, 
            FringeIntM0, FringeIntP0, T1, T2, R1, R2, RApertures, EApertures, KickAngle, a, b, Nl, Nm, Nsteps, "KQUAD_SC")
    end
end

"""
    KSEXT(;name::String = "Sext", len::Float64 = 0.0, k2::Float64 = 0.0, 
        PolynomA::Array{Float64,1} = zeros(Float64, 4), PolynomB::Array{Float64,1} = zeros(Float64, 4), 
        MaxOrder::Int64=2, NumIntSteps::Int64 = 10, rad::Int64=0, FringeQuadEntrance::Int64 = 0, 
        FringeQuadExit::Int64 = 0, FringeIntM0::Array{Float64,1} = zeros(Float64, 5), 
        FringeIntP0::Array{Float64,1} = zeros(Float64, 5), T1::Array{Float64,1} = zeros(Float64, 6), 
        T2::Array{Float64,1} = zeros(Float64, 6), R1::Array{Float64,2} = zeros(Float64, 6, 6), 
        R2::Array{Float64,2} = zeros(Float64, 6, 6), RApertures::Array{Float64,1} = zeros(Float64, 6), 
        EApertures::Array{Float64,1} = zeros(Float64, 6), KickAngle::Array{Float64,1} = zeros(Float64, 2))

A canonical sextupole element.
Example:
```julia
sext = KSEXT(name="S1", len=0.5, k2=0.5)
```
"""
mutable struct KSEXT <: AbstractElement
    name::String
    len::Float64
    k2::Float64
    PolynomA::Array{Float64,1}
    PolynomB::Array{Float64,1}
    MaxOrder::Int64
    NumIntSteps::Int64
    rad::Int64
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

    function KSEXT(;name::String = "Sext", len::Float64 = 0.0, k2::Float64 = 0.0, 
                    PolynomA::Array{Float64,1} = zeros(Float64, 4), 
                    PolynomB::Array{Float64,1} = zeros(Float64, 4), MaxOrder::Int64=2, 
                    NumIntSteps::Int64 = 10, rad::Int64=0, FringeQuadEntrance::Int64 = 0, 
                    FringeQuadExit::Int64 = 0, FringeIntM0::Array{Float64,1} = zeros(Float64, 5), 
                    FringeIntP0::Array{Float64,1} = zeros(Float64, 5), T1::Array{Float64,1} = zeros(Float64, 6), 
                    T2::Array{Float64,1} = zeros(Float64, 6), R1::Array{Float64,2} = zeros(Float64, 6, 6), 
                    R2::Array{Float64,2} = zeros(Float64, 6, 6), RApertures::Array{Float64,1} = zeros(Float64, 6), 
                    EApertures::Array{Float64,1} = zeros(Float64, 6), KickAngle::Array{Float64,1} = zeros(Float64, 2))
        if k2 != 0.0 && PolynomB[3] == 0.0
            PolynomB[3] = k2
        end
        new(name, len, k2, PolynomA, PolynomB, MaxOrder, NumIntSteps, rad, FringeQuadEntrance, FringeQuadExit, 
            FringeIntM0, FringeIntP0, T1, T2, R1, R2, RApertures, EApertures, KickAngle, "KSEXT")
    end
end

"""
    KSEXT_SC(;name::String = "Sext", len::Float64 = 0.0, k2::Float64 = 0.0, 
        PolynomA::Array{Float64,1} = zeros(Float64, 4), PolynomB::Array{Float64,1} = zeros(Float64, 4), 
        MaxOrder::Int64=2, NumIntSteps::Int64 = 10, rad::Int64=0, FringeQuadEntrance::Int64 = 0, 
        FringeQuadExit::Int64 = 0, FringeIntM0::Array{Float64,1} = zeros(Float64, 5), 
        FringeIntP0::Array{Float64,1} = zeros(Float64, 5), T1::Array{Float64,1} = zeros(Float64, 6), 
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
mutable struct KSEXT_SC <: AbstractElement
    name::String
    len::Float64
    k2::Float64
    PolynomA::Array{Float64,1}
    PolynomB::Array{Float64,1}
    MaxOrder::Int64
    NumIntSteps::Int64
    rad::Int64
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
    a::Float64
    b::Float64
    Nl::Int64
    Nm::Int64
    Nsteps::Int64
    eletype::String

    function KSEXT_SC(;name::String = "Sext", len::Float64 = 0.0, k2::Float64 = 0.0, 
                    PolynomA::Array{Float64,1} = zeros(Float64, 4), 
                    PolynomB::Array{Float64,1} = zeros(Float64, 4), MaxOrder::Int64=2, 
                    NumIntSteps::Int64 = 10, rad::Int64=0, FringeQuadEntrance::Int64 = 0, 
                    FringeQuadExit::Int64 = 0, FringeIntM0::Array{Float64,1} = zeros(Float64, 5), 
                    FringeIntP0::Array{Float64,1} = zeros(Float64, 5), T1::Array{Float64,1} = zeros(Float64, 6), 
                    T2::Array{Float64,1} = zeros(Float64, 6), R1::Array{Float64,2} = zeros(Float64, 6, 6), 
                    R2::Array{Float64,2} = zeros(Float64, 6, 6), RApertures::Array{Float64,1} = zeros(Float64, 6),
                    EApertures::Array{Float64,1} = zeros(Float64, 6), KickAngle::Array{Float64,1} = zeros(Float64, 2),
                    a::Float64 = 1.0, b::Float64 = 1.0, Nl::Int64 = 10, Nm::Int64 = 10, Nsteps::Int64=1)
        if k2 != 0.0 && PolynomB[3] == 0.0
            PolynomB[3] = k2
        end
        new(name, len, k2, PolynomA, PolynomB, MaxOrder, NumIntSteps, rad, FringeQuadEntrance, FringeQuadExit, 
            FringeIntM0, FringeIntP0, T1, T2, R1, R2, RApertures, EApertures, KickAngle, a, b, Nl, Nm, Nsteps, "KSEXT_SC")
    end
end

"""
    KOCT(;name::String = "OCT", len::Float64 = 0.0, k3::Float64 = 0.0, 
        PolynomA::Array{Float64,1} = zeros(Float64, 4), PolynomB::Array{Float64,1} = zeros(Float64, 4), 
        MaxOrder::Int64=3, NumIntSteps::Int64 = 10, rad::Int64=0, FringeQuadEntrance::Int64 = 0, 
        FringeQuadExit::Int64 = 0, FringeIntM0::Array{Float64,1} = zeros(Float64, 5), 
        FringeIntP0::Array{Float64,1} = zeros(Float64, 5), T1::Array{Float64,1} = zeros(Float64, 6), 
        T2::Array{Float64,1} = zeros(Float64, 6), R1::Array{Float64,2} = zeros(Float64, 6, 6), 
        R2::Array{Float64,2} = zeros(Float64, 6, 6), RApertures::Array{Float64,1} = zeros(Float64, 6), 
        EApertures::Array{Float64,1} = zeros(Float64, 6), KickAngle::Array{Float64,1} = zeros(Float64, 2))

A canonical octupole element.
Example:
```julia
oct = KOCT(name="O1", len=0.5, k3=0.5)
```
"""
mutable struct KOCT <: AbstractElement
    name::String
    len::Float64
    k3::Float64
    PolynomA::Array{Float64,1}
    PolynomB::Array{Float64,1}
    MaxOrder::Int64
    NumIntSteps::Int64
    rad::Int64
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

    function KOCT(;name::String = "OCT", len::Float64 = 0.0, k3::Float64 = 0.0, 
                    PolynomA::Array{Float64,1} = zeros(Float64, 4), 
                    PolynomB::Array{Float64,1} = zeros(Float64, 4), MaxOrder::Int64=3, 
                    NumIntSteps::Int64 = 10, rad::Int64=0, FringeQuadEntrance::Int64 = 0, 
                    FringeQuadExit::Int64 = 0, FringeIntM0::Array{Float64,1} = zeros(Float64, 5), 
                    FringeIntP0::Array{Float64,1} = zeros(Float64, 5), T1::Array{Float64,1} = zeros(Float64, 6), 
                    T2::Array{Float64,1} = zeros(Float64, 6), R1::Array{Float64,2} = zeros(Float64, 6, 6), 
                    R2::Array{Float64,2} = zeros(Float64, 6, 6), RApertures::Array{Float64,1} = zeros(Float64, 6), 
                    EApertures::Array{Float64,1} = zeros(Float64, 6), KickAngle::Array{Float64,1} = zeros(Float64, 2))
        if k3 != 0.0 && PolynomB[4] == 0.0
            PolynomB[4] = k3
        end
        new(name, len, k3, PolynomA, PolynomB, MaxOrder, NumIntSteps, rad, FringeQuadEntrance, FringeQuadExit, 
            FringeIntM0, FringeIntP0, T1, T2, R1, R2, RApertures, EApertures, KickAngle, "KOCT")
    end
end

"""
    KOCT_SC(;name::String = "OCT", len::Float64 = 0.0, k3::Float64 = 0.0, 
        PolynomA::Array{Float64,1} = zeros(Float64, 4), PolynomB::Array{Float64,1} = zeros(Float64, 4), 
        MaxOrder::Int64=3, NumIntSteps::Int64 = 10, rad::Int64=0, FringeQuadEntrance::Int64 = 0, 
        FringeQuadExit::Int64 = 0, FringeIntM0::Array{Float64,1} = zeros(Float64, 5), 
        FringeIntP0::Array{Float64,1} = zeros(Float64, 5), T1::Array{Float64,1} = zeros(Float64, 6), 
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
mutable struct KOCT_SC <: AbstractElement
    name::String
    len::Float64
    k3::Float64
    PolynomA::Array{Float64,1}
    PolynomB::Array{Float64,1}
    MaxOrder::Int64
    NumIntSteps::Int64
    rad::Int64
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
    a::Float64
    b::Float64
    Nl::Int64
    Nm::Int64
    Nsteps::Int64
    eletype::String

    function KOCT_SC(;name::String = "OCT", len::Float64 = 0.0, k3::Float64 = 0.0, 
                    PolynomA::Array{Float64,1} = zeros(Float64, 4), 
                    PolynomB::Array{Float64,1} = zeros(Float64, 4), MaxOrder::Int64=3, 
                    NumIntSteps::Int64 = 10, rad::Int64=0, FringeQuadEntrance::Int64 = 0, 
                    FringeQuadExit::Int64 = 0, FringeIntM0::Array{Float64,1} = zeros(Float64, 5), 
                    FringeIntP0::Array{Float64,1} = zeros(Float64, 5), T1::Array{Float64,1} = zeros(Float64, 6), 
                    T2::Array{Float64,1} = zeros(Float64, 6), R1::Array{Float64,2} = zeros(Float64, 6, 6), 
                    R2::Array{Float64,2} = zeros(Float64, 6, 6), RApertures::Array{Float64,1} = zeros(Float64, 6),
                    EApertures::Array{Float64,1} = zeros(Float64, 6), KickAngle::Array{Float64,1} = zeros(Float64, 2),
                    a::Float64 = 1.0, b::Float64 = 1.0, Nl::Int64 = 10, Nm::Int64 = 10, Nsteps::Int64=1)
        if k3 != 0.0 && PolynomB[4] == 0.0
            PolynomB[4] = k3
        end
        new(name, len, k3, PolynomA, PolynomB, MaxOrder, NumIntSteps, rad, FringeQuadEntrance, FringeQuadExit, 
            FringeIntM0, FringeIntP0, T1, T2, R1, R2, RApertures, EApertures, KickAngle, a, b, Nl, Nm, Nsteps, "KOCT_SC")
    end
end

"""
    thinMULTIPOLE(;name::String = "thinMULTIPOLE", len::Float64 = 0.0, PolynomA::Array{Float64,1} = zeros(Float64, 4), 
        PolynomB::Array{Float64,1} = zeros(Float64, 4), MaxOrder::Int64=1, NumIntSteps::Int64 = 1, rad::Int64=0, 
        FringeQuadEntrance::Int64 = 0, FringeQuadExit::Int64 = 0, FringeIntM0::Array{Float64,1} = zeros(Float64, 5), 
        FringeIntP0::Array{Float64,1} = zeros(Float64, 5), T1::Array{Float64,1} = zeros(Float64, 6), 
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
mutable struct thinMULTIPOLE <: AbstractElement
    name::String
    len::Float64
    PolynomA::Array{Float64,1}
    PolynomB::Array{Float64,1}
    MaxOrder::Int64
    NumIntSteps::Int64
    rad::Int64
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

    function thinMULTIPOLE(;name::String = "thinMULTIPOLE", len::Float64 = 0.0, PolynomA::Array{Float64,1} = zeros(Float64, 4), 
                    PolynomB::Array{Float64,1} = zeros(Float64, 4), MaxOrder::Int64=1, NumIntSteps::Int64 = 1, rad::Int64=0, 
                    FringeQuadEntrance::Int64 = 0, FringeQuadExit::Int64 = 0, FringeIntM0::Array{Float64,1} = zeros(Float64, 5), 
                    FringeIntP0::Array{Float64,1} = zeros(Float64, 5), T1::Array{Float64,1} = zeros(Float64, 6), 
                    T2::Array{Float64,1} = zeros(Float64, 6), R1::Array{Float64,2} = zeros(Float64, 6, 6), 
                    R2::Array{Float64,2} = zeros(Float64, 6, 6), RApertures::Array{Float64,1} = zeros(Float64, 6), 
                    EApertures::Array{Float64,1} = zeros(Float64, 6), KickAngle::Array{Float64,1} = zeros(Float64, 2))

        if PolynomB[3] != 0.0
            MaxOrder = 2
        end
        if PolynomB[4] != 0.0
            MaxOrder = 3
        end
        new(name, len, PolynomA, PolynomB, MaxOrder, NumIntSteps, rad, FringeQuadEntrance, FringeQuadExit,
            FringeIntM0, FringeIntP0, T1, T2, R1, R2, RApertures, EApertures, KickAngle, "thinMULTIPOLE")
    end
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
mutable struct SBEND <: AbstractElement
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

    function SBEND(;name::String = "SBend", len::Float64 = 0.0, angle::Float64 = 0.0, e1::Float64 = 0.0, e2::Float64 = 0.0, 
                    PolynomA::Array{Float64,1} = zeros(Float64, 4), PolynomB::Array{Float64,1} = zeros(Float64, 4), 
                    MaxOrder::Int64=0, NumIntSteps::Int64 = 10, rad::Int64=0, fint1::Float64 = 0.0, fint2::Float64 = 0.0, 
                    gap::Float64 = 0.0, FringeBendEntrance::Int64 = 1, FringeBendExit::Int64 = 1, 
                    FringeQuadEntrance::Int64 = 0, FringeQuadExit::Int64 = 0, FringeIntM0::Array{Float64,1} = zeros(Float64, 5), 
                    FringeIntP0::Array{Float64,1} = zeros(Float64, 5), T1::Array{Float64,1} = zeros(Float64, 6), 
                    T2::Array{Float64,1} = zeros(Float64, 6), R1::Array{Float64,2} = zeros(Float64, 6, 6), 
                    R2::Array{Float64,2} = zeros(Float64, 6, 6), RApertures::Array{Float64,1} = zeros(Float64, 6), 
                    EApertures::Array{Float64,1} = zeros(Float64, 6), KickAngle::Array{Float64,1} = zeros(Float64, 2))
        if PolynomB[2] != 0.0
            MaxOrder = 1
        end
        if PolynomB[3] != 0.0
            MaxOrder = 2
        end
        if PolynomB[4] != 0.0
            MaxOrder = 3
        end
        new(name, len, angle, e1, e2, PolynomA, PolynomB, MaxOrder, NumIntSteps, rad, fint1, fint2, gap, 
            FringeBendEntrance, FringeBendExit, FringeQuadEntrance, FringeQuadExit, FringeIntM0, FringeIntP0, 
            T1, T2, R1, R2, RApertures, EApertures, KickAngle, "SBEND")
    end
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
mutable struct SBEND_SC <: AbstractElement
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
    a::Float64
    b::Float64
    Nl::Int64
    Nm::Int64
    Nsteps::Int64
    eletype::String

    function SBEND_SC(;name::String = "SBend", len::Float64 = 0.0, angle::Float64 = 0.0, e1::Float64 = 0.0, e2::Float64 = 0.0, 
                    PolynomA::Array{Float64,1} = zeros(Float64, 4), PolynomB::Array{Float64,1} = zeros(Float64, 4), 
                    MaxOrder::Int64=0, NumIntSteps::Int64 = 10, rad::Int64=0, fint1::Float64 = 0.0, fint2::Float64 = 0.0, 
                    gap::Float64 = 0.0, FringeBendEntrance::Int64 = 1, FringeBendExit::Int64 = 1, 
                    FringeQuadEntrance::Int64 = 0, FringeQuadExit::Int64 = 0, FringeIntM0::Array{Float64,1} = zeros(Float64, 5),
                    FringeIntP0::Array{Float64,1} = zeros(Float64, 5), T1::Array{Float64,1} = zeros(Float64, 6),
                    T2::Array{Float64,1} = zeros(Float64, 6), R1::Array{Float64,2} = zeros(Float64, 6, 6),
                    R2::Array{Float64,2} = zeros(Float64, 6, 6), RApertures::Array{Float64,1} = zeros(Float64, 6),
                    EApertures::Array{Float64,1} = zeros(Float64, 6), KickAngle::Array{Float64,1} = zeros(Float64, 2),
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
        new(name, len, angle, e1, e2, PolynomA, PolynomB, MaxOrder, NumIntSteps, rad, fint1, fint2, gap, 
            FringeBendEntrance, FringeBendExit, FringeQuadEntrance, FringeQuadExit, FringeIntM0, FringeIntP0, 
            T1, T2, R1, R2, RApertures, EApertures, KickAngle, a, b, Nl, Nm, Nsteps, "SBEND_SC")
    end
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
function RBEND(;name::String = "RBend", len::Float64 = 0.0, angle::Float64 = 0.0, PolynomA::Array{Float64,1} = zeros(4), 
                PolynomB::Array{Float64,1} = zeros(4), MaxOrder::Int64=0, NumIntSteps::Int64 = 10, rad::Int64=0, fint1::Float64 = 0.0, 
                fint2::Float64 = 0.0, gap::Float64 = 0.0, FringeBendEntrance::Int64 = 1, FringeBendExit::Int64 = 1, 
                FringeQuadEntrance::Int64 = 0, FringeQuadExit::Int64 = 0, FringeIntM0::Array{Float64,1} = zeros(5), 
                FringeIntP0::Array{Float64,1} = zeros(5), T1::Array{Float64,1} = zeros(6), T2::Array{Float64,1} = zeros(6), 
                R1::Array{Float64,2} = zeros(6,6), R2::Array{Float64,2} = zeros(6,6), RApertures::Array{Float64,1} = zeros(6), 
                EApertures::Array{Float64,1} = zeros(6), KickAngle::Array{Float64,1} = zeros(2))
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
    RFCA(;name::String = "RFCA", len::Float64 = 0.0, volt::Float64 = 0.0, freq::Float64 = 0.0, h::Float64 = 1.0, 
        lag::Float64 = 0.0, philag::Float64 = 0.0, energy::Float64 = 0.0)

A RF cavity element.
Example:
```julia
rf = RFCA(name="RF1", len=0.5, volt=1e6, freq=1e6)
```
"""
mutable struct RFCA <: AbstractElement
    name::String
    len::Float64
    volt::Float64
    freq::Float64
    h::Float64
    lag::Float64
    philag::Float64
    energy::Float64
    eletype::String

    function RFCA(;name::String = "RFCA", len::Float64 = 0.0, volt::Float64 = 0.0, freq::Float64 = 0.0, h::Float64 = 1.0, 
                    lag::Float64 = 0.0, philag::Float64 = 0.0, energy::Float64 = 0.0)
        new(name, len, volt, freq, h, lag, philag, energy, "RFCA")
    end
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
mutable struct SOLENOID <: AbstractElement
    name::String
    len::Float64
    ks::Float64 # rad/m
    T1::Array{Float64,1}
    T2::Array{Float64,1}
    R1::Array{Float64,2}
    R2::Array{Float64,2}
    eletype::String

    function SOLENOID(;name::String = "Solenoid", len::Float64 = 0.0, ks::Float64 = 0.0, T1::Array{Float64,1} = zeros(6), 
                    T2::Array{Float64,1} = zeros(6), R1::Array{Float64,2} = zeros(6,6), R2::Array{Float64,2} = zeros(6,6))
        new(name, len, ks, T1, T2, R1, R2, "SOLENOID")
    end
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
mutable struct CORRECTOR <: AbstractElement
    name::String
    len::Float64
    xkick::Float64
    ykick::Float64
    T1::Array{Float64,1}
    T2::Array{Float64,1}
    R1::Array{Float64,2}
    R2::Array{Float64,2}
    eletype::String

    function CORRECTOR(;name::String = "CORRECTOR", len::Float64 = 0.0, xkick::Float64 = 0.0, ykick::Float64 = 0.0, 
                        T1::Array{Float64,1} = zeros(6), T2::Array{Float64,1} = zeros(6), R1::Array{Float64,2} = zeros(6,6), 
                        R2::Array{Float64,2} = zeros(6,6))
        new(name, len, xkick, ykick, T1, T2, R1, R2, "CORRECTOR")
    end
end

function HKICKER(;name::String = "HKicker", len::Float64 = 0.0, xkick::Float64 = 0.0)
    return CORRECTOR(name=name, len=len, xkick=xkick, ykick=0.0)
end
function VKICKER(;name::String = "VKicker", len::Float64 = 0.0, ykick::Float64 = 0.0)
    return CORRECTOR(name=name, len=len, xkick=0.0, ykick=ykick)
end

mutable struct SPACECHARGE <: AbstractElement
    # spectral space charge
    # this element is treated as an integrated effect of space charge over a length of effective_len
    name::String
    len::Float64
    effective_len::Float64
    Nl::Int64
    Nm::Int64
    a::Float64
    b::Float64
    eletype::String

    function SPACECHARGE(;name::String = "SPACECHARGE", len::Float64 = 0.0, effective_len::Float64 = 0.0, Nl::Int64 = 15, 
                        Nm::Int64 = 15, a::Float64 = 10e-3, b::Float64 = 10e-3)
        new(name, len, effective_len, Nl, Nm, a, b, "SPACECHARGE")
    end
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
mutable struct QUAD <: AbstractElement
    name::String
    len::Float64
    k1::Float64
    rad::Int64
    T1::Array{Float64,1}
    T2::Array{Float64,1}
    R1::Array{Float64,2}
    R2::Array{Float64,2}
    RApertures::Array{Float64,1}
    EApertures::Array{Float64,1}
    eletype::String

    function QUAD(;name::String = "Quad", len::Float64 = 0.0, k1::Float64 = 0.0, rad::Int64 = 0, 
                    T1::Array{Float64,1} = zeros(6), T2::Array{Float64,1} = zeros(6), 
                    R1::Array{Float64,2} = zeros(6,6), R2::Array{Float64,2} = zeros(6,6), 
                    RApertures::Array{Float64,1} = zeros(6), EApertures::Array{Float64,1} = zeros(6))
        new(name, len, k1, rad, T1, T2, R1, R2, RApertures, EApertures, "QUAD")
    end
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
mutable struct QUAD_SC <: AbstractElement
    name::String
    len::Float64
    k1::Float64
    rad::Int64
    T1::Array{Float64,1}
    T2::Array{Float64,1}
    R1::Array{Float64,2}
    R2::Array{Float64,2}
    RApertures::Array{Float64,1}
    EApertures::Array{Float64,1}
    a::Float64
    b::Float64
    Nl::Int64
    Nm::Int64
    Nsteps::Int64
    eletype::String
    
    function QUAD_SC(;name::String = "Quad", len::Float64 = 0.0, k1::Float64 = 0.0, rad::Int64 = 0, 
                    T1::Array{Float64,1} = zeros(6), T2::Array{Float64,1} = zeros(6), 
                    R1::Array{Float64,2} = zeros(6,6), R2::Array{Float64,2} = zeros(6,6), 
                    RApertures::Array{Float64,1} = zeros(6), EApertures::Array{Float64,1} = zeros(6), 
                    a::Float64 = 1.0, b::Float64 = 1.0, Nl::Int64 = 10, Nm::Int64 = 10, Nsteps::Int64=1)
        new(name, len, k1, rad, T1, T2, R1, R2, RApertures, EApertures, a, b, Nl, Nm, Nsteps, "QUAD_SC")
    end
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
mutable struct CRABCAVITY <: AbstractElement
    name::String 
    len::Float64 
    volt::Float64  # voltage
    freq::Float64  # frequency
    k::Float64  # wave number
    phi::Float64  # phase
    errors::Array{Float64,1} # 1: Voltage error, 2: Phase error
    energy::Float64
    eletype::String 
    function CRABCAVITY(;name::String = "CRABCAVITY", len::Float64 = 0.0, volt::Float64 = 0.0, 
        freq::Float64 = 0.0, phi::Float64 = 0.0, errors::Array{Float64,1} = zeros(2), energy::Float64 = 1e9)
        k = 2*π*freq/2.99792458e8
        return new(name, len, volt, freq, k, phi, errors, energy, "CRABCAVITY")
    end
end


mutable struct easyCRABCAVITY <: AbstractElement
    name::String 
    len::Float64 
    halfthetac::Float64 
    freq::Float64 
    k::Float64 
    phi::Float64 
    errors::Array{Float64,1}  # 1: Voltage error, 2: Phase error
    eletype::String
    function easyCRABCAVITY(;name::String = "easyCRABCAVITY", len::Float64 = 0.0, halfthetac::Float64 = 0.0, 
        freq::Float64 = 0.0, k::Float64 = 0.0, phi::Float64 = 0.0, errors::Array{Float64,1} = zeros(2), eletype::String = "easyCRABCAVITY")
        return new(name, len, halfthetac, freq, k, phi, errors, eletype)
    end
end

mutable struct AccelCavity <: AbstractElement
    name::String 
    len::Float64 
    volt::Float64  # voltage
    freq::Float64  # frequency
    k::Float64  # wave number
    h::Float64  # harmonic number
    phis::Float64  # synchronous phase π/2 for accelerating on crest
    eletype::String 
    function AccelCavity(;name::String = "AccelCavity", len::Float64 = 0.0, volt::Float64 = 0.0, 
        freq::Float64 = 0.0, h::Float64 = 1.0, phis::Float64 = 0.0)
        k = 2*π*freq/2.99792458e8
        return new(name, len, volt, freq, k, h, phis, "AccelCavity")
    end
end


abstract type AbstractTransferMap <:AbstractElement end
abstract type AbstractTransverseMap <:AbstractTransferMap end
abstract type AbstractLongitudinalMap <:AbstractTransferMap end
mutable struct LongitudinalRFMap <: AbstractLongitudinalMap
    alphac::Float64
    RF::AbstractElement
    LongitudinalRFMap(alphac::Float64, RF::AbstractElement)=new(alphac, RF)
end

mutable struct LorentzBoost <: AbstractElement
    angle::Float64
    cosang::Float64
    tanang::Float64
    mode::Int
    LorentzBoost(angle)=new(angle, cos(angle), tan(angle), 0)
end

mutable struct InvLorentzBoost <: AbstractElement
    angle::Float64
    sinang::Float64
    cosang::Float64
    mode::Int
    InvLorentzBoost(angle)=new(angle, sin(angle), cos(angle), 0)
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
abstract type AbstractStrongBeamBeam <:AbstractElement end

mutable struct StrongThinGaussianBeam <: AbstractStrongBeamBeam
    amplitude::Float64
    rmssizex::Float64
    rmssizey::Float64
    zloc::Float64
    xoffset::Float64
    yoffset::Float64
    StrongThinGaussianBeam(amp::Float64, rx::Float64, ry::Float64, zloc::Float64=0.0, xoff::Float64=0.0, yoff::Float64=0.0)=new(amp,rx,ry,zloc,xoff,yoff)
end

"""
    StrongGaussianBeam(charge::Float64, mass::Float64, atomnum::Float64, np::Int, energy::Float64, op::AbstractOptics4D, bs::Vector{Float64}, nz::Int)

Construct a strong beam-beam element with Gaussian distribution.
"""
mutable struct StrongGaussianBeam <: AbstractStrongBeamBeam  # Strong Beam with transverse Gaussian distribution
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
mutable struct LongitudinalRLCWake <: AbstractElement
    freq::Float64
    Rshunt::Float64
    Q0::Float64
    LongitudinalRLCWake(;freq::Float64=1.0e9, Rshunt::Float64=1.0e6, Q0::Float64=1.0)=new(freq, Rshunt, Q0)
end
function wakefieldfunc_RLCWake(rlcwake::LongitudinalRLCWake, t::Float64)
    Q0p=sqrt(rlcwake.Q0^2 - 1.0/4.0)
    w0 = 2*pi*rlcwake.freq
    w0p= w0/rlcwake.Q0*Q0p
    if t>0
        return 0.0
    else
        return rlcwake.Rshunt * w0 /rlcwake.Q0 * (cos(w0p * t) +  sin(w0p * t) / 2 / Q0p) * exp(w0 * t / 2 / rlcwake.Q0)
    end
end

mutable struct LongitudinalWake <: AbstractElement
    times::AbstractVector
    wakefields::AbstractVector
    wakefield::Function
end
"""
    LongitudinalWake(times::AbstractVector, wakefields::AbstractVector, wakefield::Function)

Create longitudinal wake element.
"""
function LongitudinalWake(times::AbstractVector, wakefields::AbstractVector, fliphalf::Float64=-1.0)
    wakefield = function (t::Float64)
        t>times[1]*fliphalf && return 0.0
        return linear_interpolate(t*fliphalf, times, wakefields)
    end
    return LongitudinalWake(times, wakefields, wakefield)
end