abstract type AbstractElement end
struct EDRIFT <: AbstractElement
    name::String
    len::Float64
end
function EDRIFT(; name="EDrift", len=0.0)
    EDRIFT(name, len)
end     

struct KQUAD <: AbstractElement
    name::String                       # element name
    k1::Float64                        # quadrupole strength
    len::Float64                       # length of quadrupole
    bore::Float64                      # bore radius
    fse::Float64                       # fractional strength error
    tilt::Float64                      # rotation about longitudinal axis
    dx::Float64                        # misalignment
    dy::Float64                        # misalignment
    dz::Float64                        # misalignment
    # pitch::Float64
    # yaw::Float64
    edge1_effects::Int16               # include entrance edge effects?
    edge2_effects::Int16               # include exit edge effects?
    edge1Linear::Int16                 # selectively turn off linear part if EDGE1_EFFECTS nonzero.
    edge1NonlinearFactor::Float64      # selectively scale nonlinear entrance edge effects if EDGE1_EFFECTS>1
    edge2Linear::Int16                 # selectively turn off linear part if EDGE2_EFFECTS nonzero.
    edge2NonlinearFactor::Float64      # selectively scale nonlinear exit edge effects if EDGE2_EFFECTS>1 
    fringeIntM::Array{Float64, 1}      # + fringe integral
    fringeIntP::Array{Float64, 1}      # - fringe integral
    radial::Int16                      # include radial dependence? not yet implemented
    integration_order::Int16           # integration order
    nSlices::Int16                     # number of integrator steps
    xkick::Float64                     # kick in x
    ykick::Float64                     # kick in y
    synch_rad::Int16                   # include synchrotron radiation?
    isr::Int16                         # include incoherent synchrotron radiation (quantum excitation)?
end
function KQUAD(; name="KQuad", k1=0.0, len=0.0, bore=0.0, fse=0.0, tilt=0.0, dx=0.0, dy=0.0, dz=0.0, 
                edge1_effects=0, edge2_effects=0, edge1Linear=1, edge1NonlinearFactor=1.0, 
                edge2Linear=1, edge2NonlinearFactor=1.0, fringeIntM=zeros(Float64, 5), fringeIntP=zeros(Float64, 5), 
                radial=0, integration_order=4, nSlices=1, xkick=0.0, ykick=0.0, synch_rad=0, isr=0)
    KQUAD(name, k1, len, bore, fse, tilt, dx, dy, dz, edge1_effects, edge2_effects, edge1Linear, edge1NonlinearFactor, edge2Linear, edge2NonlinearFactor, fringeIntM, fringeIntP, radial, integration_order, nSlices, xkick, ykick, synch_rad, isr)
end
function KQUAD(name, k1, len; bore=0.0, fse=0.0, tilt=0.0, dx=0.0, dy=0.0, dz=0.0, 
    edge1_effects=0, edge2_effects=0, edge1Linear=1, edge1NonlinearFactor=1.0, 
    edge2Linear=1, edge2NonlinearFactor=1.0, fringeIntM=zeros(Float64, 5), fringeIntP=zeros(Float64, 5), 
    radial=0, integration_order=4, nSlices=1, xkick=0.0, ykick=0.0, synch_rad=0, isr=0)
    KQUAD(name, k1, len, bore, fse, tilt, dx, dy, dz, edge1_effects, edge2_effects, edge1Linear, edge1NonlinearFactor, edge2Linear, edge2NonlinearFactor, fringeIntM, fringeIntP, radial, integration_order, nSlices, xkick, ykick, synch_rad, isr)
end

struct KSEXT <: AbstractElement
    name::String                       # element name
    k2::Float64                        # sextupole strength
    len::Float64                       # length of sextupole
    bore::Float64                      # bore radius
    fse::Float64                       # fractional strength error
    tilt::Float64                      # rotation about longitudinal axis
    dx::Float64                        # misalignment
    dy::Float64                        # misalignment
    dz::Float64                        # misalignment
    radial::Int16                      # include radial dependence? not yet implemented
    integration_order::Int16           # integration order
    nSlices::Int16                     # number of integrator steps
    xkick::Float64                     # kick in x
    ykick::Float64                     # kick in y
    synch_rad::Int16                   # include synchrotron radiation?
    isr::Int16                         # include incoherent synchrotron radiation (quantum excitation)? not implemented yet
end
function KSEXT(; name="KSext", k2=0.0, len=0.0, bore=0.0, fse=0.0, tilt=0.0, dx=0.0, dy=0.0, dz=0.0, 
                radial=0, integration_order=4, nSlices=1, xkick=0.0, ykick=0.0, synch_rad=0, isr=0)
    KSEXT(name, k2, len, bore, fse, tilt, dx, dy, dz, radial, integration_order, nSlices, xkick, ykick, synch_rad, isr)
end
function KSEXT(name, k2, len; bore=0.0, fse=0.0, tilt=0.0, dx=0.0, dy=0.0, dz=0.0, 
    radial=0, integration_order=4, nSlices=1, xkick=0.0, ykick=0.0, synch_rad=0, isr=0)
    KSEXT(name, k2, len, bore, fse, tilt, dx, dy, dz, radial, integration_order, nSlices, xkick, ykick, synch_rad, isr)
end

struct KOCT <: AbstractElement 
    name::String                       # element name
    k3::Float64                        # octupole strength
    len::Float64                       # length of octupole
    bore::Float64                      # bore radius
    fse::Float64                       # fractional strength error
    tilt::Float64                      # rotation about longitudinal axis
    dx::Float64                        # misalignment
    dy::Float64                        # misalignment
    dz::Float64                        # misalignment
    radial::Int16                      # include radial dependence? not yet implemented
    integration_order::Int16           # integration order
    nSlices::Int16                     # number of integrator steps
    xkick::Float64                     # kick in x
    ykick::Float64                     # kick in y
    synch_rad::Int16                   # include synchrotron radiation?
    isr::Int16                         # include incoherent synchrotron radiation (quantum excitation)? not implemented yet
end
function KOCT(; name="KOct", k3=0.0, len=0.0, bore=0.0, fse=0.0, tilt=0.0, dx=0.0, dy=0.0, dz=0.0, 
                radial=0, integration_order=4, nSlices=1, xkick=0.0, ykick=0.0, synch_rad=0, isr=0)
    KOCT(name, k3, len, bore, fse, tilt, dx, dy, dz, radial, integration_order, nSlices, xkick, ykick, synch_rad, isr)
end
function KOCT(name, k3, len; bore=0.0, fse=0.0, tilt=0.0, dx=0.0, dy=0.0, dz=0.0, 
    radial=0, integration_order=4, nSlices=1, xkick=0.0, ykick=0.0, synch_rad=0, isr=0)
    KOCT(name, k3, len, bore, fse, tilt, dx, dy, dz, radial, integration_order, nSlices, xkick, ykick, synch_rad, isr)
end

struct CSBEND <: AbstractElement
    name::String
    length::Float64
    angle::Float64
    e1::Float64
    e2::Float64
    fint::Float64
    hgap::Float64
    k1::Float64
    k2::Float64
    k3::Float64
    k4::Float64
    k5::Float64
    k6::Float64
    k7::Float64
    k8::Float64
    b1::Float64
    b2::Float64
    b3::Float64
    b4::Float64
    b5::Float64
    b6::Float64
    b7::Float64
    b8::Float64
    dx::Float64
    dy::Float64
    dz::Float64
    edge1_effects::Int64
    edge2_effects::Int64
    edge_order::Int64
    e1_kick_limit::Float64
    e2_kick_limit::Float64
    kick_limit_scaling::Float64
    nonlinear::Int64
    use_bn::Int64
    synch_rad::Int64
    isr::Int64
    isr1Particle::Int64
    distributionBasedRadiation::Int64
    includeOpeningAngle::Int64
    tilt::Float64
    etilt::Float64
    fse::Float64
    h1::Float64
    h2::Float64
    integration_order::Int64
    nSlice::Int64
    expansionOrder::Int64
end
function CSBEND(; name="CSBend", length=0.0, angle=0.0, e1=0.0, e2=0.0, fint=0.5, hgap=0.0, 
                k1=0.0, k2=0.0, k3=0.0, k4=0.0, k5=0.0, k6=0.0, k7=0.0, k8=0.0, 
                b1=0.0, b2=0.0, b3=0.0, b4=0.0, b5=0.0, b6=0.0, b7=0.0, b8=0.0, 
                dx=0.0, dy=0.0, dz=0.0, edge1_effects=1, edge2_effects=1, edge_order=1, 
                e1_kick_limit=-1.0, e2_kick_limit=-1.0, kick_limit_scaling=0.0, nonlinear=1, 
                use_bn=0, synch_rad=0, isr=0, isr1Particle=1, distributionBasedRadiation=0, 
                includeOpeningAngle=1, tilt=0.0, etilt=0.0, fse=0.0, h1=0.0, h2=0.0, integration_order=4, 
                nSlice=1, expansionOrder=0)
    CSBEND(name, length, angle, e1, e2, fint, hgap, k1, k2, k3, k4, k5, k6, k7, k8, b1, b2, b3, b4, b5, b6, b7, b8, dx, dy, dz, edge1_effects, edge2_effects, edge_order, e1_kick_limit, e2_kick_limit, kick_limit_scaling, nonlinear, use_bn, synch_rad, isr, isr1Particle, distributionBasedRadiation, includeOpeningAngle, tilt, etilt, fse, h1, h2, integration_order, nSlice, expansionOrder)
end
function CSBEND(name, length, angle; e1=0.0, e2=0.0, fint=0.5, hgap=0.0, 
    k1=0.0, k2=0.0, k3=0.0, k4=0.0, k5=0.0, k6=0.0, k7=0.0, k8=0.0, 
    b1=0.0, b2=0.0, b3=0.0, b4=0.0, b5=0.0, b6=0.0, b7=0.0, b8=0.0, 
    dx=0.0, dy=0.0, dz=0.0, edge1_effects=1, edge2_effects=1, edge_order=1, 
    e1_kick_limit=-1.0, e2_kick_limit=-1.0, kick_limit_scaling=0.0, nonlinear=1, 
    use_bn=0, synch_rad=0, isr=0, isr1Particle=1, distributionBasedRadiation=0, 
    includeOpeningAngle=1, tilt=0.0, etilt=0.0, fse=0.0, h1=0.0, h2=0.0, integration_order=4, 
    nSlice=1, expansionOrder=0)
    CSBEND(name, length, angle, e1, e2, fint, hgap, k1, k2, k3, k4, k5, k6, k7, k8, b1, b2, b3, b4, b5, b6, b7, b8, dx, dy, dz, edge1_effects, edge2_effects, edge_order, e1_kick_limit, e2_kick_limit, kick_limit_scaling, nonlinear, use_bn, synch_rad, isr, isr1Particle, distributionBasedRadiation, includeOpeningAngle, tilt, etilt, fse, h1, h2, integration_order, nSlice, expansionOrder)
end