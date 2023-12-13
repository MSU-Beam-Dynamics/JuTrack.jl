using Base: @kwdef
abstract type AbstractElement end
struct EDRIFT <: AbstractElement
    name::String
    len::Float64
end
function EDRIFT(; name="EDrift", len=0.0)
    EDRIFT(name, len)
end     

@kwdef struct KQUAD <: AbstractElement
    name::String = "KQuad"                              # element name  
    k1::Float64 = 0.0                                   # quadrupole strength  
    len::Float64 = 0.0                                  # length of quadrupole
    bore::Float64 = 0.0                                 # bore radius
    fse::Float64 = 0.0                       
    tilt::Float64 = 0.0                      
    dx::Float64 = 0.0                        
    dy::Float64 = 0.0                        
    dz::Float64 = 0.0                        
    edge1_effects::Int64 = 0                
    edge2_effects::Int64 = 0                
    edge1Linear::Int64 = 1                
    edge1NonlinearFactor::Float64 = 1.0      
    edge2Linear::Int64 = 1                
    edge2NonlinearFactor::Float64 = 1.0      
    fringeIntM::Array{Float64, 1} = zeros(Float64, 5)      
    fringeIntP::Array{Float64, 1} = zeros(Float64, 5)      
    radial::Int64 = 0                     
    integration_order::Int64 = 4          
    nSlices::Int64 = 4                    
    xkick::Float64 = 0.0                    
    ykick::Float64 = 0.0                    
    synch_rad::Int64 = 0                  
    isr::Int64 = 0                   
    isr1Particle::Int64 = 1       
end
function KQUAD(name, len, k1; bore=0.0, fse=0.0, tilt=0.0, dx=0.0, dy=0.0, dz=0.0, 
    edge1_effects=0, edge2_effects=0, edge1Linear=1, edge1NonlinearFactor=1.0, 
    edge2Linear=1, edge2NonlinearFactor=1.0, fringeIntM=zeros(Float64, 5), fringeIntP=zeros(Float64, 5), 
    radial=0, integration_order=4, nSlices=4, xkick=0.0, ykick=0.0, synch_rad=0, isr=0)
    KQUAD(name, k1, len, bore, fse, tilt, dx, dy, dz, edge1_effects, edge2_effects, edge1Linear, edge1NonlinearFactor, edge2Linear, edge2NonlinearFactor, fringeIntM, fringeIntP, radial, integration_order, nSlices, xkick, ykick, synch_rad, isr)
end
function KQUAD(name, len, k1, dx, dy, dz, edge1_effects, edge2_effects, nSlices; bore=0.0, fse=0.0, tilt=0.0, 
    edge1Linear=1, edge1NonlinearFactor=1.0, 
    edge2Linear=1, edge2NonlinearFactor=1.0, 
    fringeIntM=zeros(Float64, 5), fringeIntP=zeros(Float64, 5), 
    radial=0, integration_order=4, xkick=0.0, ykick=0.0, synch_rad=0, isr=0)
    KQUAD(name, k1, len, bore, fse, tilt, dx, dy, dz, edge1_effects, edge2_effects, edge1Linear, edge1NonlinearFactor, edge2Linear, edge2NonlinearFactor, fringeIntM, fringeIntP, radial, integration_order, nSlices, xkick, ykick, synch_rad, isr)
end

@kwdef struct KSEXT <: AbstractElement
    name::String                 = "KSext"      # element name
    k2::Float64                  = 0.0          # sextupole strength
    len::Float64                 = 0.0          # length of sextupole
    bore::Float64                = 0.0          # bore radius
    fse::Float64                 = 0.0          # fractional strength error
    tilt::Float64                = 0.0          # rotation about longitudinal axis
    dx::Float64                  = 0.0          # misalignment
    dy::Float64                  = 0.0          # misalignment
    dz::Float64                  = 0.0          # misalignment
    radial::Int64                = 0.0          # include radial dependence? not yet implemented
    integration_order::Int64     = 4            # integration order
    nSlices::Int64               = 4            # number of integrator steps
    xkick::Float64               = 0.0          # kick in x
    ykick::Float64               = 0.0          # kick in y
    synch_rad::Int64             = 0            # include synchrotron radiation?
    isr::Int64                   = 0            # include incoherent synchrotron radiation (quantum excitation)? not implemented yet
end
function KSEXT(name, len, k2; bore=0.0, fse=0.0, tilt=0.0, dx=0.0, dy=0.0, dz=0.0, 
    radial=0, integration_order=4, nSlices=4, xkick=0.0, ykick=0.0, synch_rad=0, isr=0)
    KSEXT(name, k2, len, bore, fse, tilt, dx, dy, dz, radial, integration_order, nSlices, xkick, ykick, synch_rad, isr)
end
function KSEXT(name, len, k2, dx, dy, dz, nSlices; bore=0.0, fse=0.0, tilt=0.0, 
    radial=0, integration_order=4, xkick=0.0, ykick=0.0, synch_rad=0, isr=0)
    KSEXT(name, k2, len, bore, fse, tilt, dx, dy, dz, radial, integration_order, nSlices, xkick, ykick, synch_rad, isr)
end

@kwdef struct KOCT <: AbstractElement 
    name::String    = "KOct"               # element name
    k3::Float64     = 0.0                  # octupole strength
    len::Float64    =0.0                   # length of octupole
    bore::Float64   =0.0                   # bore radius
    fse::Float64    =0.0                   # fractional strength error
    tilt::Float64   =0.0                   # rotation about longitudinal axis
    dx::Float64     =0.0                   # misalignment
    dy::Float64     =0.0                   # misalignment
    dz::Float64     =0.0                   # misalignment
    radial::Int64   =0.0                   # include radial dependence? not yet implemented
    integration_order::Int64   = 4         # integration order
    nSlices::Int64             = 4         # number of integrator steps
    xkick::Float64             = 0.0       # kick in x
    ykick::Float64             = 0.0       # kick in y
    synch_rad::Int64           = 0         # include synchrotron radiation?
    isr::Int64                 = 0         # include incoherent synchrotron radiation (quantum excitation)? not implemented yet
end
function KOCT(name, len, k3; bore=0.0, fse=0.0, tilt=0.0, dx=0.0, dy=0.0, dz=0.0, 
    radial=0, integration_order=4, nSlices=4, xkick=0.0, ykick=0.0, synch_rad=0, isr=0)
    KOCT(name, k3, len, bore, fse, tilt, dx, dy, dz, radial, integration_order, nSlices, xkick, ykick, synch_rad, isr)
end
function KOCT(name, len, k3, dx, dy, dz, nSlices; bore=0.0, fse=0.0, tilt=0.0,
    radial=0, integration_order=4, xkick=0.0, ykick=0.0, synch_rad=0, isr=0)
    KOCT(name, k3, len, bore, fse, tilt, dx, dy, dz, radial, integration_order, nSlices, xkick, ykick, synch_rad, isr)
end

@kwdef struct CSBEND <: AbstractElement
    name::String = "CSBend"
    len::Float64 = 0.0
    angle::Float64 = 0.0
    e1::Float64 = 0.0
    e2::Float64 = 0.0
    fint::Float64 = 0.5
    hgap::Float64 = 0.0
    k1::Float64 = 0.0
    k2::Float64 = 0.0
    k3::Float64 = 0.0
    k4::Float64 = 0.0
    k5::Float64 = 0.0
    k6::Float64 = 0.0
    k7::Float64 = 0.0
    k8::Float64 = 0.0
    b1::Float64 = 0.0
    b2::Float64 = 0.0
    b3::Float64 = 0.0
    b4::Float64 = 0.0
    b5::Float64 = 0.0
    b6::Float64 = 0.0
    b7::Float64 = 0.0
    b8::Float64 = 0.0
    dx::Float64 = 0.0
    dy::Float64 = 0.0
    dz::Float64 = 0.0
    edge1_effects::Int64 = 1
    edge2_effects::Int64 = 1
    edge_order::Int64 = 1
    e1_kick_limit::Float64 = -1.0
    e2_kick_limit::Float64 = -1.0
    kick_limit_scaling::Float64 = 0.0
    nonlinear::Int64 = 1
    use_bn::Int64 = 0
    synch_rad::Int64 = 0
    isr::Int64 = 0
    isr1Particle::Int64 = 1
    distributionBasedRadiation::Int64 = 0
    includeOpeningAngle::Int64 = 1
    tilt::Float64 = 0.0
    etilt::Float64 = 0.0
    fse::Float64 = 0.0
    h1::Float64 = 0.0
    h2::Float64 = 0.0
    integration_order::Int64 = 4
    nSlices::Int64 = 4
    expansionOrder::Int64 = 0
end
function CSBEND(name, length, angle; e1=0.0, e2=0.0, fint=0.5, hgap=0.0, 
    k1=0.0, k2=0.0, k3=0.0, k4=0.0, k5=0.0, k6=0.0, k7=0.0, k8=0.0, 
    b1=0.0, b2=0.0, b3=0.0, b4=0.0, b5=0.0, b6=0.0, b7=0.0, b8=0.0, 
    dx=0.0, dy=0.0, dz=0.0, edge1_effects=1, edge2_effects=1, edge_order=1, 
    e1_kick_limit=-1.0, e2_kick_limit=-1.0, kick_limit_scaling=0.0, nonlinear=1, 
    use_bn=0, synch_rad=0, isr=0, isr1Particle=1, distributionBasedRadiation=0, 
    includeOpeningAngle=1, tilt=0.0, etilt=0.0, fse=0.0, h1=0.0, h2=0.0, integration_order=4, 
    nSlices=4, expansionOrder=0)
    CSBEND(name, length, angle, e1, e2, fint, hgap, k1, k2, k3, k4, k5, k6, k7, k8, b1, b2, b3, b4, b5, b6, b7, b8, dx, dy, dz, edge1_effects, edge2_effects, edge_order, e1_kick_limit, e2_kick_limit, kick_limit_scaling, nonlinear, use_bn, synch_rad, isr, isr1Particle, distributionBasedRadiation, includeOpeningAngle, tilt, etilt, fse, h1, h2, integration_order, nSlices, expansionOrder)
end
function CSBEND(name, length, angle, e1, e2, nSlices; fint=0.5, hgap=0.0, 
    k1=0.0, k2=0.0, k3=0.0, k4=0.0, k5=0.0, k6=0.0, k7=0.0, k8=0.0, 
    b1=0.0, b2=0.0, b3=0.0, b4=0.0, b5=0.0, b6=0.0, b7=0.0, b8=0.0, 
    dx=0.0, dy=0.0, dz=0.0, edge1_effects=1, edge2_effects=1, edge_order=1, 
    e1_kick_limit=-1.0, e2_kick_limit=-1.0, kick_limit_scaling=0.0, nonlinear=1, 
    use_bn=0, synch_rad=0, isr=0, isr1Particle=1, distributionBasedRadiation=0, 
    includeOpeningAngle=1, tilt=0.0, etilt=0.0, fse=0.0, h1=0.0, h2=0.0, integration_order=4, expansionOrder=0)
    CSBEND(name, length, angle, e1, e2, fint, hgap, k1, k2, k3, k4, k5, k6, k7, k8, b1, b2, b3, b4, b5, b6, b7, b8, dx, dy, dz, edge1_effects, edge2_effects, edge_order, e1_kick_limit, e2_kick_limit, kick_limit_scaling, nonlinear, use_bn, synch_rad, isr, isr1Particle, distributionBasedRadiation, includeOpeningAngle, tilt, etilt, fse, h1, h2, integration_order, nSlices, expansionOrder)
end

@kwdef struct RFCA <: AbstractElement
    name::String = "RFCA"
    freq::Float64 = 0.0
    volt::Float64 = 0.0
    phase::Float64 = 0.0
    phase_reference::Int64 = 0
    phase_fiducial::Float64 = 0.0
    fiducial_mode = nothing  
    tReference::Float64 = -1.0
    Q::Float64 = 0.0
    len::Float64 = 0.0
    nSlices::Int64 = 1
    dx::Float64 = 0.0
    dy::Float64 = 0.0
    change_p0::Int64 = 0
    change_t::Int64 = 0
    linearize::Int64 = 0
    lockPhase::Int64 = 0
    end1Focus::Int64 = 0
    end2Focus::Int64 = 0
    bodyFocusModel::String = "none"
end
function RFCA(name, len, freq, volt, phase, nSlices; phase_reference=0, phase_fiducial=0.0, fiducial_mode=nothing, 
                tReference=-1.0, Q=0.0, dx=0.0, dy=0.0, change_p0=0, change_t=0, 
                linearize=0, lockPhase=0, end1Focus=0, end2Focus=0, bodyFocusModel="none")
    RFCA(name, freq, volt, phase, phase_reference, phase_fiducial, fiducial_mode, tReference, Q, len, nSlices, dx, dy, change_p0, change_t, linearize, lockPhase, end1Focus, end2Focus, bodyFocusModel)
end
