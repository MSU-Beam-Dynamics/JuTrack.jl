abstract type AbstractElement end
struct EDrift <: AbstractElement
    name::String
    len::Float64
end

struct KQUAD <: AbstractElement
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


struct KSEXT <: AbstractElement
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

struct KOCT <: AbstractElement 
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

struct CSBEND <: AbstractElement
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