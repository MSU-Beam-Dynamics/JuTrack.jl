# Native-coordinate compatibility elements for measured-machine optics models.
# This model is only used for the HEPS model calibration and should not be 
# directly transferred to other machines.
# No external application, machine configuration, or lattice is loaded here.
using LinearAlgebra: I

function _fringe_array(value, shape, label)
    a = Float64.(value)
    size(a) == shape || throw(ArgumentError("$label must have shape $shape"))
    all(isfinite, a) || throw(ArgumentError("$label must be finite"))
    return copy(a)
end

function _fringe_scalars(values...)
    all(isfinite, values) || throw(ArgumentError("Fringe-element fields must be finite"))
end

"""
    QuadLinearFintFudgePass(; name="QFF", len=0., k1=0., fudge=1.,
                             fin=I(6), fout=I(6), T1, T2, R1, R2, ...)

Finite analytic quadrupole with effective strength `k1*fudge`, preceded and
followed by explicit 6x6 fringe matrices. All vectors and matrices use native
JuTrack coordinates [x, px, y, py, z, delta], NOT AT coordinate order.
`fin`/`fout` default to identity, not zero. External rotations follow JuTrack's
zero-matrix-means-disabled convention. k1 is in m^-2; fudge is dimensionless.

This compatibility map retains the legacy paraxial path-length convention even
when the global drift mode is exact. Float64 tracking and finite-difference
optics/LOCO are supported; radiation and Taylor/parameter-AD tracking are not.
"""
mutable struct QuadLinearFintFudgePass{T} <: AbstractElement{T}
    name::String
    len::T
    k1::T
    fudge::T
    fin::Matrix{T}
    fout::Matrix{T}
    T1::Vector{T}
    T2::Vector{T}
    R1::Matrix{T}
    R2::Matrix{T}
    RApertures::Vector{Float64}
    EApertures::Vector{Float64}
    rad::Int64
    eletype::String
end

function QuadLinearFintFudgePass(;name::String="QFF", len=0., k1=0., fudge=1.,
    fin=Matrix{Float64}(I,6,6), fout=Matrix{Float64}(I,6,6),
    T1=zeros(6), T2=zeros(6), R1=zeros(6,6), R2=zeros(6,6),
    RApertures=zeros(6), EApertures=zeros(6), rad::Int=0)
    _fringe_scalars(len,k1,fudge)
    len >= 0 || throw(ArgumentError("Quadrupole length must be nonnegative"))
    rad == 0 || throw(ArgumentError("QuadLinearFintFudgePass does not implement radiation"))
    return QuadLinearFintFudgePass{Float64}(name,Float64(len),Float64(k1),Float64(fudge),
        _fringe_array(fin,(6,6),"fin"),_fringe_array(fout,(6,6),"fout"),
        _fringe_array(T1,(6,),"T1"),_fringe_array(T2,(6,),"T2"),
        _fringe_array(R1,(6,6),"R1"),_fringe_array(R2,(6,6),"R2"),
        _fringe_array(RApertures,(6,),"RApertures"),_fringe_array(EApertures,(6,),"EApertures"),
        0,"QuadLinearFintFudgePass")
end

"""
    BendLinearFriTrapPass(;name="BFT",len,angle,k1=0.,kq=0.,i1a=0.,i1b=0.,...)

Analytic combined-function bend with legacy trapezoidal quadrupole fringe maps.
`angle`, entrance_angle and exit_angle are radians; len/gap are metres; k1/kq
are m^-2, and kq*i1a / kq*i1b are dimensionless fringe factors. `byerror` is
the legacy normalized dipole-field error. Entrance sequence: nonlinear
quadrupole fringe (when i1a!=0), linear fringe, dipole edge; exit is the reverse
sequence with exit signs. The nonlinear fringe is NOT QuadFringePassP/N.

All coordinates, translations and rotations are native [x,px,y,py,z,delta].
This is a compatibility map for Float64 tracking and finite-difference optics;
radiation integrals/tracking and Taylor/parameter-AD tracking are not supported.
"""
mutable struct BendLinearFriTrapPass{T} <: AbstractElement{T}
    name::String
    len::T
    angle::T
    entrance_angle::T
    exit_angle::T
    k1::T
    kq::T
    i1a::T
    i1b::T
    fint1::T
    fint2::T
    gap::T
    byerror::T
    T1::Vector{T}
    T2::Vector{T}
    R1::Matrix{T}
    R2::Matrix{T}
    RApertures::Vector{Float64}
    EApertures::Vector{Float64}
    rad::Int64
    eletype::String
end

function BendLinearFriTrapPass(;name::String="BFT",len,angle,k1=0.,kq=0.,i1a=0.,i1b=0.,
    entrance_angle=0.,exit_angle=0.,fint1=0.,fint2=0.,gap=0.,byerror=0.,
    T1=zeros(6),T2=zeros(6),R1=zeros(6,6),R2=zeros(6,6),
    RApertures=zeros(6),EApertures=zeros(6),rad::Int=0)
    _fringe_scalars(len,angle,k1,kq,i1a,i1b,entrance_angle,exit_angle,fint1,fint2,gap,byerror)
    len > 0 || throw(ArgumentError("Bend length must be positive"))
    gap >= 0 || throw(ArgumentError("Full gap must be nonnegative"))
    rad == 0 || throw(ArgumentError("BendLinearFriTrapPass does not implement radiation"))
    return BendLinearFriTrapPass{Float64}(name,Float64(len),Float64(angle),
        Float64(entrance_angle),Float64(exit_angle),Float64(k1),Float64(kq),
        Float64(i1a),Float64(i1b),Float64(fint1),Float64(fint2),Float64(gap),Float64(byerror),
        _fringe_array(T1,(6,),"T1"),_fringe_array(T2,(6,),"T2"),
        _fringe_array(R1,(6,6),"R1"),_fringe_array(R2,(6,6),"R2"),
        _fringe_array(RApertures,(6,),"RApertures"),_fringe_array(EApertures,(6,),"EApertures"),
        0,"BendLinearFriTrapPass")
end
