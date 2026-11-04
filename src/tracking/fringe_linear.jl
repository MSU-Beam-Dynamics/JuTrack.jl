# The analytic body maps reuse JuTrack's QUAD/LBEND kernels. The additional
# fringe operations implement the legacy FintFudge/FriTrap map convention in
# native coordinates. Do not replace this sequence by Elegant fringe maps.
const _MachineFringeElement = Union{QuadLinearFintFudgePass,BendLinearFriTrapPass}

function _machine_fringe_core!(e::QuadLinearFintFudgePass{Float64},r,c)
    multmv_row!(r,c,e.fin)
    if !all(isfinite,@view(r[c,:])) || r[c,6]<=-1
        r[c,:].=NaN
        return
    end
    strength=e.k1*e.fudge
    if iszero(strength)
        # The legacy zero-strength limit is paraxial, independent of the global
        # exact-Hamiltonian drift switch (essential for compatibility).
        norm=1/(1+r[c,6])
        r[c,1]+=e.len*norm*r[c,2]
        r[c,3]+=e.len*norm*r[c,4]
        r[c,5]+=e.len*norm^2*(r[c,2]^2+r[c,4]^2)/2
    else
        _quad_linear_transport_row!(r,c,e.len,strength)
    end
    multmv_row!(r,c,e.fout)
end

function _fritrap_linear_fringe!(r,integral)
    factor=1+integral/(1+r[6])
    if !isfinite(factor) || iszero(factor)
        r .= NaN
        return
    end
    r[1]*=factor; r[2]/=factor
    r[3]/=factor; r[4]*=factor
end

function _fritrap_nonlinear_fringe!(r,k)
    # Rotate by 45 degrees, apply the two ordered shear maps, rotate back.
    a=inv(sqrt(2.))
    x,px,y,py=r[1],r[2],r[3],r[4]
    r[1]=(x+y)*a; r[2]=(px+py)*a
    r[3]=(-x+y)*a; r[4]=(-px+py)*a
    norm=1/(1+r[6])
    x,px,y,py=r[1],r[2],r[3],r[4]
    r[1]=x-k*y^3*norm/6
    r[4]=py+k*y^2*px*norm/2
    r[2]=px+k*r[1]^2*r[4]*norm/2
    r[3]=y-k*r[1]^3*norm/6
    x,px,y,py=r[1],r[2],r[3],r[4]
    r[1]=(x-y)*a; r[2]=(px-py)*a
    r[3]=(x+y)*a; r[4]=(px+py)*a
end

function _machine_fringe_core!(e::BendLinearFriTrapPass{Float64},r,c)
    v=@view r[c,:]
    !iszero(e.i1a) && _fritrap_nonlinear_fringe!(v,e.kq)
    _fritrap_linear_fringe!(v,-e.kq*e.i1a)
    all(isfinite,v) || return
    edge_fringe_entrance!(v,e.angle/e.len,e.entrance_angle,e.fint1,e.gap,1)
    bend6!(v,e.len,e.angle,e.k1,e.byerror)
    edge_fringe_exit!(v,e.angle/e.len,e.exit_angle,e.fint2,e.gap,1)
    _fritrap_linear_fringe!(v,e.kq*e.i1b)
    all(isfinite,v) || return
    !iszero(e.i1b) && _fritrap_nonlinear_fringe!(v,-e.kq)
end

function _machine_fringe_row!(e,r,c,lost)
    lost[c]==1 && return
    if !all(isfinite,@view(r[c,:])) || r[c,6]<=-1
        lost[c]=1; return
    end
    !iszero(e.T1) && addvv_row!(r,c,e.T1)
    !iszero(e.R1) && multmv_row!(r,c,e.R1)
    if !all(isfinite,@view(r[c,:])) || r[c,6]<=-1
        lost[c]=1; return
    end
    _machine_fringe_core!(e,r,c)
    !iszero(e.R2) && multmv_row!(r,c,e.R2)
    !iszero(e.T2) && addvv_row!(r,c,e.T2)
    if !all(isfinite,@view(r[c,:])) || _check_lost_row(r,c) ||
       _check_lost_aperture_row(r,c,e.RApertures,e.EApertures)
        lost[c]=1
    end
end

function _check_machine_fringe(e)
    e.rad==0 || throw(ArgumentError("Radiation is unsupported for $(e.eletype)"))
    use_exact_beti==0 || throw(ArgumentError("$(e.eletype) requires delta=p/p0-1 coordinates"))
    isfinite(e.len) && e.len>=0 || throw(ArgumentError("Invalid fringe-element length"))
    if e isa BendLinearFriTrapPass
        e.len>0 || throw(ArgumentError("Bend length must be positive"))
    end
end

function pass!(e::_MachineFringeElement,r::Matrix{Float64},n::Int64,b::Beam{Float64})
    _check_machine_fringe(e)
    for c in 1:n
        _machine_fringe_row!(e,r,c,b.lost_flag)
    end
    nothing
end

function pass_P!(e::_MachineFringeElement,r::Matrix{Float64},n::Int64,b::Beam{Float64})
    _check_machine_fringe(e)
    Threads.@threads for c in 1:n
        _machine_fringe_row!(e,r,c,b.lost_flag)
    end
    nothing
end

function pass_TPSA!(e::_MachineFringeElement,r::AbstractVector;kwargs...)
    throw(ArgumentError("$(e.eletype) currently supports Float64 tracking and finite-difference optics, not Taylor tracking"))
end

function elrad(e::Union{QuadLinearFintFudgePass{T},BendLinearFriTrapPass{T}},
               initial::EdwardsTengTwiss{T},final::EdwardsTengTwiss{T}) where T
    throw(ArgumentError("Radiation integrals are not validated for $(e.eletype); do not silently substitute an ordinary element"))
end
