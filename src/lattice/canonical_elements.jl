using Base: @kwdef
abstract type AbstractElement end

@kwdef struct DRIFT <: AbstractElement
    name::String = "EDrift"
    len::Float64 = 0.0
    T1::Array{Float64,1} = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    T2::Array{Float64,1} = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    R1::Array{Float64,2} = [1.0 0.0 0.0 0.0 0.0 0.0; 
                            0.0 1.0 0.0 0.0 0.0 0.0; 
                            0.0 0.0 1.0 0.0 0.0 0.0; 
                            0.0 0.0 0.0 1.0 0.0 0.0; 
                            0.0 0.0 0.0 0.0 1.0 0.0; 
                            0.0 0.0 0.0 0.0 0.0 1.0]
    R2::Array{Float64,2} = [1.0 0.0 0.0 0.0 0.0 0.0;
                            0.0 1.0 0.0 0.0 0.0 0.0;
                            0.0 0.0 1.0 0.0 0.0 0.0;
                            0.0 0.0 0.0 1.0 0.0 0.0;
                            0.0 0.0 0.0 0.0 1.0 0.0;
                            0.0 0.0 0.0 0.0 0.0 1.0]         
    RApertures::Array{Float64,1} = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    EApertures::Array{Float64,1} = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
end


@kwdef struct KQUAD <: AbstractElement
    name::String  = "Quad"                                      # element name  
    len::Float64 = 0.0
    k1::Float64 = 0.0                                           # use k1 if PolynomB is not given
    PolynomA::Array{Float64,1} = [0.0, 0.0, 0.0, 0.0]    
    PolynomB::Array{Float64,1} = [0.0, 0.0, 0.0, 0.0]           # PolynomB has higher priority than k1, k2, k3
    MaxOrder::Int64 = 1
    NumIntSteps::Int64 = 10
    FringeQuadEntrance::Int64 = 0
    FringeQuadExit::Int64 = 0
    FringeIntM0::Array{Float64,1} = [0.0, 0.0, 0.0, 0.0, 0.0]
    FringeIntP0::Array{Float64,1} = [0.0, 0.0, 0.0, 0.0, 0.0]
    T1::Array{Float64,1} = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    T2::Array{Float64,1} = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    R1::Array{Float64,2} = [1.0 0.0 0.0 0.0 0.0 0.0; 
                            0.0 1.0 0.0 0.0 0.0 0.0; 
                            0.0 0.0 1.0 0.0 0.0 0.0; 
                            0.0 0.0 0.0 1.0 0.0 0.0; 
                            0.0 0.0 0.0 0.0 1.0 0.0; 
                            0.0 0.0 0.0 0.0 0.0 1.0]
    R2::Array{Float64,2} = [1.0 0.0 0.0 0.0 0.0 0.0;
                            0.0 1.0 0.0 0.0 0.0 0.0;
                            0.0 0.0 1.0 0.0 0.0 0.0;
                            0.0 0.0 0.0 1.0 0.0 0.0;
                            0.0 0.0 0.0 0.0 1.0 0.0;
                            0.0 0.0 0.0 0.0 0.0 1.0]         
    RApertures::Array{Float64,1} = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    EApertures::Array{Float64,1} = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    KickAngle::Array{Float64,1} = [0.0, 0.0]
end

@kwdef struct KSEXT <: AbstractElement
    name::String  = "Sext"                                      # element name  
    len::Float64 = 0.0
    k2::Float64 = 0.0                                           # use k2 if PolynomB is not given
    PolynomA::Array{Float64,1} = [0.0, 0.0, 0.0, 0.0]    
    PolynomB::Array{Float64,1} = [0.0, 0.0, 0.0, 0.0]           # PolynomB has higher priority than k1, k2, k3
    MaxOrder::Int64 = 2
    NumIntSteps::Int64 = 10
    FringeQuadEntrance::Int64 = 0
    FringeQuadExit::Int64 = 0
    FringeIntM0::Array{Float64,1} = [0.0, 0.0, 0.0, 0.0, 0.0]
    FringeIntP0::Array{Float64,1} = [0.0, 0.0, 0.0, 0.0, 0.0]
    T1::Array{Float64,1} = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    T2::Array{Float64,1} = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    R1::Array{Float64,2} = [1.0 0.0 0.0 0.0 0.0 0.0; 
                            0.0 1.0 0.0 0.0 0.0 0.0; 
                            0.0 0.0 1.0 0.0 0.0 0.0; 
                            0.0 0.0 0.0 1.0 0.0 0.0; 
                            0.0 0.0 0.0 0.0 1.0 0.0; 
                            0.0 0.0 0.0 0.0 0.0 1.0]
    R2::Array{Float64,2} = [1.0 0.0 0.0 0.0 0.0 0.0;
                            0.0 1.0 0.0 0.0 0.0 0.0;
                            0.0 0.0 1.0 0.0 0.0 0.0;
                            0.0 0.0 0.0 1.0 0.0 0.0;
                            0.0 0.0 0.0 0.0 1.0 0.0;
                            0.0 0.0 0.0 0.0 0.0 1.0]         
    RApertures::Array{Float64,1} = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    EApertures::Array{Float64,1} = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    KickAngle::Array{Float64,1} = [0.0, 0.0]
end

@kwdef struct KOCT <: AbstractElement
    name::String  = "OCT"                                       # element name  
    len::Float64 = 0.0
    k3::Float64 = 0.0                                           # use k3 if PolynomB is not given
    PolynomA::Array{Float64,1} = [0.0, 0.0, 0.0, 0.0]    
    PolynomB::Array{Float64,1} = [0.0, 0.0, 0.0, 0.0]           # PolynomB has higher priority than k1, k2, k3
    MaxOrder::Int64 = 3
    NumIntSteps::Int64 = 10
    FringeQuadEntrance::Int64 = 0
    FringeQuadExit::Int64 = 0
    FringeIntM0::Array{Float64,1} = [0.0, 0.0, 0.0, 0.0, 0.0]
    FringeIntP0::Array{Float64,1} = [0.0, 0.0, 0.0, 0.0, 0.0]
    T1::Array{Float64,1} = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    T2::Array{Float64,1} = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    R1::Array{Float64,2} = [1.0 0.0 0.0 0.0 0.0 0.0; 
                            0.0 1.0 0.0 0.0 0.0 0.0; 
                            0.0 0.0 1.0 0.0 0.0 0.0; 
                            0.0 0.0 0.0 1.0 0.0 0.0; 
                            0.0 0.0 0.0 0.0 1.0 0.0; 
                            0.0 0.0 0.0 0.0 0.0 1.0]
    R2::Array{Float64,2} = [1.0 0.0 0.0 0.0 0.0 0.0;
                            0.0 1.0 0.0 0.0 0.0 0.0;
                            0.0 0.0 1.0 0.0 0.0 0.0;
                            0.0 0.0 0.0 1.0 0.0 0.0;
                            0.0 0.0 0.0 0.0 1.0 0.0;
                            0.0 0.0 0.0 0.0 0.0 1.0]         
    RApertures::Array{Float64,1} = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    EApertures::Array{Float64,1} = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    KickAngle::Array{Float64,1} = [0.0, 0.0]
end

@kwdef struct SBEND <: AbstractElement
    name::String = "SBend"
    len::Float64 = 0.0
    angle::Float64 = 0.0
    e1::Float64 = 0.0
    e2::Float64 = 0.0
    PolynomA::Array{Float64,1} = [0.0, 0.0, 0.0, 0.0]    
    PolynomB::Array{Float64,1} = [0.0, 0.0, 0.0, 0.0]
    MaxOrder::Int64 = 0
    NumIntSteps::Int64 = 10
    fint1::Float64 = 0.0
    fint2::Float64 = 0.0
    gap::Float64 = 0.0
    FringeBendEntrance::Int64 = 0
    FringeBendExit::Int64 = 0
    FringeQuadEntrance::Int64 = 0
    FringeQuadExit::Int64 = 0
    FringeIntM0::Array{Float64,1} = [0.0, 0.0, 0.0, 0.0, 0.0]
    FringeIntP0::Array{Float64,1} = [0.0, 0.0, 0.0, 0.0, 0.0]
    T1::Array{Float64,1} = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    T2::Array{Float64,1} = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    R1::Array{Float64,2} = [1.0 0.0 0.0 0.0 0.0 0.0; 
                            0.0 1.0 0.0 0.0 0.0 0.0; 
                            0.0 0.0 1.0 0.0 0.0 0.0; 
                            0.0 0.0 0.0 1.0 0.0 0.0; 
                            0.0 0.0 0.0 0.0 1.0 0.0; 
                            0.0 0.0 0.0 0.0 0.0 1.0]
    R2::Array{Float64,2} = [1.0 0.0 0.0 0.0 0.0 0.0;
                            0.0 1.0 0.0 0.0 0.0 0.0;
                            0.0 0.0 1.0 0.0 0.0 0.0;
                            0.0 0.0 0.0 1.0 0.0 0.0;
                            0.0 0.0 0.0 0.0 1.0 0.0;
                            0.0 0.0 0.0 0.0 0.0 1.0]         
    RApertures::Array{Float64,1} = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    EApertures::Array{Float64,1} = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    KickAngle::Array{Float64,1} = [0.0, 0.0]
end


@kwdef struct RFCA <: AbstractElement
    name::String = "RFCA"
    len::Float64 = 0.0
    volt::Float64 = 0.0
    freq::Float64 = 0.0
    h::Float64 = 1.0
    lag::Float64 = 0.0
    philag::Float64 = 0.0
    energy::Float64 = 0.0
end
