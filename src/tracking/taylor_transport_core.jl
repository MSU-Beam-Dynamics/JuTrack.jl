@inline _taylor_branch_value(x) = x
@inline _taylor_branch_value(x::AbstractHOTPSA) = jet_primal(x)
@inline _taylor_iszero(x) = _taylor_branch_value(x) == 0
@inline _taylor_allzero(xs) = all(_taylor_iszero, xs)

@inline function _taylor_addvv!(r::AbstractVector{S}, dr::AbstractVector) where {S}
    @inbounds begin
        r[1] += dr[1]
        r[2] += dr[2]
        r[3] += dr[3]
        r[4] += dr[4]
        r[5] += dr[5]
        r[6] += dr[6]
    end
    return nothing
end

@inline function _taylor_multmv!(r::AbstractVector{S}, A::AbstractMatrix) where {S}
    @inbounds begin
        r1_in = r[1]
        r2_in = r[2]
        r3_in = r[3]
        r4_in = r[4]
        r5_in = r[5]
        r6_in = r[6]

        r1 = A[1, 1] * r1_in + A[1, 2] * r2_in + A[1, 3] * r3_in + A[1, 4] * r4_in + A[1, 5] * r5_in + A[1, 6] * r6_in
        r2 = A[2, 1] * r1_in + A[2, 2] * r2_in + A[2, 3] * r3_in + A[2, 4] * r4_in + A[2, 5] * r5_in + A[2, 6] * r6_in
        r3 = A[3, 1] * r1_in + A[3, 2] * r2_in + A[3, 3] * r3_in + A[3, 4] * r4_in + A[3, 5] * r5_in + A[3, 6] * r6_in
        r4 = A[4, 1] * r1_in + A[4, 2] * r2_in + A[4, 3] * r3_in + A[4, 4] * r4_in + A[4, 5] * r5_in + A[4, 6] * r6_in
        r5 = A[5, 1] * r1_in + A[5, 2] * r2_in + A[5, 3] * r3_in + A[5, 4] * r4_in + A[5, 5] * r5_in + A[5, 6] * r6_in
        r6 = A[6, 1] * r1_in + A[6, 2] * r2_in + A[6, 3] * r3_in + A[6, 4] * r4_in + A[6, 5] * r5_in + A[6, 6] * r6_in

        r[1] = r1
        r[2] = r2
        r[3] = r3
        r[4] = r4
        r[5] = r5
        r[6] = r6
    end
    return nothing
end

@inline function _taylor_linearized_drift_phifac(px, py, dp_p, gamma2i)
    denom = 1.0 + dp_p
    if use_exact_Hamiltonian == 2
        phifac = (px^2 + py^2 + dp_p^2 * gamma2i) / 2
        return (phifac / denom - dp_p * gamma2i) / denom
    end
    return (px^2 + py^2) / (2 * denom^2)
end

@inline function _taylor_fastdrift!(r::AbstractVector{S}, NormL, le, beti, gamma2i=0.0) where {S}
    if use_exact_Hamiltonian == 1
        r[1] += NormL * r[2]
        r[3] += NormL * r[4]
        r[5] += NormL * (beti + r[6]) - le * beti
    else
        r[1] += NormL * r[2]
        r[3] += NormL * r[4]
        r[5] += le * _taylor_linearized_drift_phifac(r[2], r[4], r[6], gamma2i)
    end
    return nothing
end

@inline function _taylor_drift6!(r::AbstractVector{S}, le, beti, gamma2i=0.0) where {S}
    if use_exact_Hamiltonian == 1
        NormL = le / sqrt(1.0 + 2.0 * r[6] * beti + r[6]^2 - r[2]^2 - r[4]^2)
        r[5] += NormL * (beti + r[6]) - le * beti
    else
        NormL = le / (1.0 + r[6])
        r[5] += le * _taylor_linearized_drift_phifac(r[2], r[4], r[6], gamma2i)
    end
    r[1] += NormL * r[2]
    r[3] += NormL * r[4]
    return nothing
end

@inline function _taylor_drift_pass!(r::AbstractVector{S}, le, beti, gamma2i,
    T1::AbstractVector, T2::AbstractVector,
    R1::AbstractMatrix, R2::AbstractMatrix) where {S}
    hasT1 = !_taylor_allzero(T1)
    hasR1 = !_taylor_allzero(R1)
    hasR2 = !_taylor_allzero(R2)
    hasT2 = !_taylor_allzero(T2)

    hasT1 && _taylor_addvv!(r, T1)
    hasR1 && _taylor_multmv!(r, R1)
    _taylor_drift6!(r, le, beti, gamma2i)
    hasR2 && _taylor_multmv!(r, R2)
    hasT2 && _taylor_addvv!(r, T2)
    return nothing
end

@inline function _taylor_corrector_pass!(r::AbstractVector{S}, le, xkick, ykick,
    T1::AbstractVector, T2::AbstractVector,
    R1::AbstractMatrix, R2::AbstractMatrix) where {S}
    hasT1 = !_taylor_allzero(T1)
    hasR1 = !_taylor_allzero(R1)
    hasR2 = !_taylor_allzero(R2)
    hasT2 = !_taylor_allzero(T2)

    hasT1 && _taylor_addvv!(r, T1)
    hasR1 && _taylor_multmv!(r, R1)

    p_norm = 1.0 / (1.0 + r[6])
    NormL = le * p_norm
    r[5] += NormL * p_norm * (xkick * xkick / 3 + ykick * ykick / 3 + r[2] * r[2] + r[4] * r[4] + r[2] * xkick + r[4] * ykick) / 2
    r[1] += NormL * (r[2] + xkick / 2)
    r[2] += xkick
    r[3] += NormL * (r[4] + ykick / 2)
    r[4] += ykick

    hasR2 && _taylor_multmv!(r, R2)
    hasT2 && _taylor_addvv!(r, T2)
    return nothing
end

@inline function _taylor_quad_linear_transport!(r::AbstractVector{S}, le, k1) where {S}
    p_norm = 1.0 / (1.0 + r[6])
    x = r[1]
    xpr = r[2] * p_norm
    y = r[3]
    ypr = r[4] * p_norm

    g = abs(k1) / (1.0 + r[6])
    t = sqrt(g)
    lt = le * t

    if k1 > 0
        MHD = cos(lt)
        M12 = sin(lt) / t
        M21 = -M12 * g
        MVD = cosh(lt)
        M34 = sinh(lt) / t
        M43 = M34 * g
    else
        MHD = cosh(lt)
        M12 = sinh(lt) / t
        M21 = M12 * g
        MVD = cos(lt)
        M34 = sin(lt) / t
        M43 = -M34 * g
    end

    r[1] = MHD * x + M12 * xpr
    r[2] = (M21 * x + MHD * xpr) / p_norm
    r[3] = MVD * y + M34 * ypr
    r[4] = (M43 * y + MVD * ypr) / p_norm
    r[5] += g * (x * x * (le - MHD * M12) - y * y * (le - MVD * M34)) / 4
    r[5] += (xpr * xpr * (le + MHD * M12) + ypr * ypr * (le + MVD * M34)) / 4
    r[5] += (x * xpr * M12 * M21 + y * ypr * M34 * M43) / 2
    return nothing
end

@inline function _taylor_quad_linear_pass!(r::AbstractVector{S}, le, k1, beti, gamma2i,
    T1::AbstractVector, T2::AbstractVector,
    R1::AbstractMatrix, R2::AbstractMatrix) where {S}
    hasT1 = !_taylor_allzero(T1)
    hasR1 = !_taylor_allzero(R1)
    hasR2 = !_taylor_allzero(R2)
    hasT2 = !_taylor_allzero(T2)

    hasT1 && _taylor_addvv!(r, T1)
    hasR1 && _taylor_multmv!(r, R1)

    if _taylor_iszero(k1)
        _taylor_drift6!(r, le, beti, gamma2i)
    else
        _taylor_quad_linear_transport!(r, le, k1)
    end

    hasR2 && _taylor_multmv!(r, R2)
    hasT2 && _taylor_addvv!(r, T2)
    return nothing
end

@inline function _taylor_translation_pass!(r::AbstractVector{S}, dx, dy, ds, beta) where {S}
    pz = sqrt(1.0 + 2.0 * r[6] / beta + r[6]^2 - r[2]^2 - r[4]^2)
    r[1] -= dx + ds * r[2] / pz
    r[3] -= dy + ds * r[4] / pz
    r[5] += ds * (1.0 / beta + r[6]) / pz
    return nothing
end

@inline function _taylor_yrotation_pass!(r::AbstractVector{S}, angle, beta) where {S}
    if _taylor_iszero(angle)
        return nothing
    end

    ca = cos(-angle)
    sa = sin(-angle)
    ta = tan(-angle)
    x, px, y, py, t, pt = r[1], r[2], r[3], r[4], r[5], r[6]
    pz = sqrt(1.0 + 2.0 * pt / beta + pt^2 - px^2 - py^2)
    ptt = 1.0 - ta * px / pz
    r[1] = x / (ca * ptt)
    r[2] = ca * px + sa * pz
    r[3] = y + ta * x * py / (pz * ptt)
    r[5] = t + ta * x * (1.0 / beta + pt) / (pz * ptt)
    return nothing
end

@inline function _taylor_strthinkick!(r::AbstractVector{S}, A::AbstractVector, B::AbstractVector,
    L, max_order::Int, a1_shift=0.0, b1_shift=0.0) where {S}
    ReSum = B[max_order + 1]
    ImSum = A[max_order + 1]
    ReSumTemp = zero(ReSum)

    for i in (max_order - 1):-1:0
        bcoeff = i == 0 ? B[1] + b1_shift : B[i + 1]
        acoeff = i == 0 ? A[1] + a1_shift : A[i + 1]
        ReSumTemp = ReSum * r[1] - ImSum * r[3] + bcoeff
        ImSum = ImSum * r[1] + ReSum * r[3] + acoeff
        ReSum = ReSumTemp
    end

    r[2] -= L * ReSum
    r[4] += L * ImSum
    return nothing
end

@inline function _taylor_thinmultipole_pass!(r::AbstractVector{S}, A::AbstractVector, B::AbstractVector,
    max_order::Int, T1::AbstractVector, T2::AbstractVector,
    R1::AbstractMatrix, R2::AbstractMatrix, KickAngle::AbstractVector) where {S}
    hasT1 = !_taylor_allzero(T1)
    hasR1 = !_taylor_allzero(R1)
    hasR2 = !_taylor_allzero(R2)
    hasT2 = !_taylor_allzero(T2)

    hasT1 && _taylor_addvv!(r, T1)
    hasR1 && _taylor_multmv!(r, R1)
    _taylor_strthinkick!(r, A, B, 1.0, max_order, KickAngle[2], -KickAngle[1])
    hasR2 && _taylor_multmv!(r, R2)
    hasT2 && _taylor_addvv!(r, T2)
    return nothing
end

@inline function _taylor_multipole_integrator_pass!(r::AbstractVector{S}, le, beti, gamma2i,
    A::AbstractVector, B::AbstractVector, max_order::Int, num_int_step::Int,
    FringeQuadEntrance::Int, FringeQuadExit::Int,
    T1::AbstractVector, T2::AbstractVector,
    R1::AbstractMatrix, R2::AbstractMatrix, KickAngle::AbstractVector) where {S}
    if FringeQuadEntrance != 0 || FringeQuadExit != 0
        throw(ArgumentError("HOTPSA multipole fringe fields are not implemented yet"))
    end

    DRIFT1 = 0.6756035959798286638
    DRIFT2 = -0.1756035959798286639
    KICK1 = 1.351207191959657328
    KICK2 = -1.702414383919314656

    SL = le / num_int_step
    L1 = SL * DRIFT1
    L2 = SL * DRIFT2
    K1 = SL * KICK1
    K2 = SL * KICK2

    a1_shift = zero(eltype(A))
    b1_shift = zero(eltype(B))
    if !_taylor_iszero(le)
        a1_shift = sin(KickAngle[2]) / le
        b1_shift = -sin(KickAngle[1]) / le
    end

    hasT1 = !_taylor_allzero(T1)
    hasR1 = !_taylor_allzero(R1)
    hasR2 = !_taylor_allzero(R2)
    hasT2 = !_taylor_allzero(T2)

    hasT1 && _taylor_addvv!(r, T1)
    hasR1 && _taylor_multmv!(r, R1)
    for _ in 1:num_int_step
        _taylor_drift6!(r, L1, beti, gamma2i)
        _taylor_strthinkick!(r, A, B, K1, max_order, a1_shift, b1_shift)
        _taylor_drift6!(r, L2, beti, gamma2i)
        _taylor_strthinkick!(r, A, B, K2, max_order, a1_shift, b1_shift)
        _taylor_drift6!(r, L2, beti, gamma2i)
        _taylor_strthinkick!(r, A, B, K1, max_order, a1_shift, b1_shift)
        _taylor_drift6!(r, L1, beti, gamma2i)
    end
    hasR2 && _taylor_multmv!(r, R2)
    hasT2 && _taylor_addvv!(r, T2)
    return nothing
end

@inline function _taylor_solenoid_pass!(r::AbstractVector{S}, le, ks, beti,
    T1::AbstractVector, T2::AbstractVector,
    R1::AbstractMatrix, R2::AbstractMatrix) where {S}
    hasT1 = !_taylor_allzero(T1)
    hasR1 = !_taylor_allzero(R1)
    hasR2 = !_taylor_allzero(R2)
    hasT2 = !_taylor_allzero(T2)

    hasT1 && _taylor_addvv!(r, T1)
    hasR1 && _taylor_multmv!(r, R1)

    if _taylor_iszero(ks)
        _taylor_drift6!(r, le, beti)
    else
        p_norm = 1.0 / (1.0 + r[6])
        x = r[1]
        xpr = r[2] * p_norm
        y = r[3]
        ypr = r[4] * p_norm
        H = ks * p_norm / 2.0
        S1 = sin(le * H)
        C1 = cos(le * H)
        Hinv = inv(H)

        r[1] = x * C1 * C1 + xpr * C1 * S1 * Hinv + y * C1 * S1 + ypr * S1 * S1 * Hinv
        r[2] = (-x * H * C1 * S1 + xpr * C1 * C1 - y * H * S1 * S1 + ypr * C1 * S1) / p_norm
        r[3] = -x * C1 * S1 - xpr * S1 * S1 * Hinv + y * C1 * C1 + ypr * C1 * S1 * Hinv
        r[4] = (x * H * S1 * S1 - xpr * C1 * S1 - y * C1 * S1 * H + ypr * C1 * C1) / p_norm
        r[5] += le * (H * H * (x * x + y * y) + 2.0 * H * (xpr * y - ypr * x) + xpr * xpr + ypr * ypr) / 2.0
    end

    hasR2 && _taylor_multmv!(r, R2)
    hasT2 && _taylor_addvv!(r, T2)
    return nothing
end

@inline function _taylor_rfca_pass!(r::AbstractVector{S}, le, nv, freq, lag, philag, beta, beti) where {S}
    C0 = 2.99792458e8
    arg = 2.0 * pi * freq * ((r[5] - lag) / C0) - philag
    if _taylor_iszero(le)
        r[6] += -nv * sin(arg) / beta^2
    else
        halflength = le / 2.0
        _taylor_drift6!(r, halflength, beti)
        r[6] += -nv * sin(arg) / beta^2
        _taylor_drift6!(r, halflength, beti)
    end
    return nothing
end

@inline function _taylor_crabcavity_pass!(r::AbstractVector{S}, le, volt, k, phi, E0, beta, beti) where {S}
    ang = k * r[5] + phi
    if _taylor_iszero(le)
        r[2] += (volt / E0) * sin(ang / beta)
        r[6] += (-k * volt / E0 / beta) * r[1] * cos(ang / beta)
    else
        halflength = le / 2.0
        _taylor_drift6!(r, halflength, beti)
        r[2] += (volt / E0) * sin(ang / beta)
        r[6] += (-k * volt / E0 / beta) * r[1] * cos(ang / beta)
        _taylor_drift6!(r, halflength, beti)
    end
    return nothing
end

@inline function _taylor_crabcavity_k2_pass!(r::AbstractVector{S}, le, volt, k, phi, k2, E0, beta, beti) where {S}
    ang = k * r[5] + phi
    if _taylor_iszero(le)
        r[2] += (volt / E0) * sin(ang / beta)
        r[6] += (-k * volt / E0 / beta) * r[1] * cos(ang / beta)
        r[2] -= k2 * (r[1]^2 - r[3]^2) * sin(ang)
        r[4] += 2.0 * k2 * r[1] * r[3] * sin(ang)
        r[6] -= (k2 * k / 3.0) * (r[1]^3 - 3.0 * r[1] * r[3]^2) * cos(ang)
    else
        halflength = le / 2.0
        _taylor_drift6!(r, halflength, beti)
        r[2] += (volt / E0) * sin(ang / beta)
        r[6] += (-k * volt / E0 / beta) * r[1] * cos(ang / beta)
        r[2] -= k2 * (r[1]^2 - r[3]^2) * sin(ang)
        r[4] += 2.0 * k2 * r[1] * r[3] * sin(ang)
        r[6] -= (k2 * k / 3.0) * (r[1]^3 - 3.0 * r[1] * r[3]^2) * cos(ang)
        _taylor_drift6!(r, halflength, beti)
    end
    return nothing
end

@inline function _taylor_easycrabcavity_pass!(r::AbstractVector{S}, halfthetac, k, phi) where {S}
    ang = k * r[5] + phi
    r[1] += (halfthetac / k) * sin(ang)
    r[6] += halfthetac * r[2] * cos(ang)
    return nothing
end

@inline function _taylor_accelcavity_pass!(r::AbstractVector{S}, volt, k, phis, beta2E) where {S}
    sv = sin(k * r[5] + phis) - sin(phis)
    r[6] += volt / beta2E * sv
    return nothing
end

@inline function _taylor_lorentz_boost_pass!(r::AbstractVector{S}, cosang, tanang, mode::Int) where {S}
    if mode == 0
        invcosang = 1.0 / cosang
        r[1] += tanang * r[5]
        r[6] -= tanang * r[2]
        r[2] *= invcosang
        r[4] *= invcosang
        r[5] *= invcosang
    end
    return nothing
end

@inline function _taylor_inv_lorentz_boost_pass!(r::AbstractVector{S}, sinang, cosang, mode::Int) where {S}
    if mode == 0
        r[1] -= sinang * r[5]
        r[6] += sinang * r[2]
        r[2] *= cosang
        r[4] *= cosang
        r[5] *= cosang
    end
    return nothing
end

@inline function _taylor_bend6!(r::AbstractVector{S}, L, b_angle, grd, ByError) where {S}
    Kx = b_angle / L
    p_norm = 1.0 / (1.0 + r[6])
    G1 = (Kx * Kx + grd) * p_norm
    G2 = -grd * p_norm

    MHD = zero(r[1])
    M12 = zero(r[1])
    M21 = zero(r[1])
    MVD = zero(r[1])
    M34 = zero(r[1])
    M43 = zero(r[1])
    arg1 = zero(r[1])
    sqrtG1 = zero(r[1])
    arg2 = zero(r[1])
    sqrtG2 = zero(r[1])

    G1p = _taylor_branch_value(G1)
    if G1p == 0.0
        MHD = one(r[1])
        M12 = L
        M21 = zero(r[1])
    elseif G1p > 0.0
        sqrtG1 = sqrt(G1)
        arg1 = L * sqrtG1
        MHD = cos(arg1)
        M12 = sin(arg1) / sqrtG1
        M21 = -sin(arg1) * sqrtG1
    else
        sqrtG1 = sqrt(-G1)
        arg1 = L * sqrtG1
        MHD = cosh(arg1)
        M12 = sinh(arg1) / sqrtG1
        M21 = sinh(arg1) * sqrtG1
    end

    G2p = _taylor_branch_value(G2)
    if G2p == 0.0
        MVD = one(r[1])
        M34 = L
        M43 = zero(r[1])
    elseif G2p > 0.0
        sqrtG2 = sqrt(G2)
        arg2 = L * sqrtG2
        MVD = cos(arg2)
        M34 = sin(arg2) / sqrtG2
        M43 = -sin(arg2) * sqrtG2
    else
        sqrtG2 = sqrt(-G2)
        arg2 = L * sqrtG2
        MVD = cosh(arg2)
        M34 = sinh(arg2) / sqrtG2
        M43 = sinh(arg2) * sqrtG2
    end

    x = r[1]
    xpr = r[2] * p_norm
    y = r[3]
    ypr = r[4] * p_norm
    delta = r[6]
    dispersion_error = delta * p_norm - ByError

    r[1] = MHD * x + M12 * xpr
    r[2] = (M21 * x + MHD * xpr) / p_norm

    if G1p == 0.0
        r[1] += dispersion_error * L * L * Kx / 2.0
        r[2] += dispersion_error * L * Kx / p_norm
    elseif G1p > 0.0
        r[1] += dispersion_error * (1.0 - cos(arg1)) * Kx / G1
        r[2] += dispersion_error * sin(arg1) * Kx / (sqrtG1 * p_norm)
    else
        r[1] += dispersion_error * (1.0 - cosh(arg1)) * Kx / G1
        r[2] += dispersion_error * sinh(arg1) * Kx / (sqrtG1 * p_norm)
    end

    r[3] = MVD * y + M34 * ypr
    r[4] = (M43 * y + MVD * ypr) / p_norm
    r[5] += xpr * xpr * (L + MHD * M12) / 4.0

    if G1p != 0.0
        r[5] += (L - MHD * M12) * (x * x * G1 + dispersion_error * dispersion_error * Kx * Kx / G1 -
            2.0 * x * Kx * dispersion_error) / 4.0
        r[5] += M12 * M21 * (x * xpr - xpr * dispersion_error * Kx / G1) / 2.0
        r[5] += Kx * x * M12 + xpr * (1.0 - MHD) * Kx / G1 + dispersion_error * (L - M12) * Kx * Kx / G1
    end
    r[5] += ((L - MVD * M34) * y * y * G2 + ypr * ypr * (L + MVD * M34)) / 4.0
    r[5] += M34 * M43 * x * xpr / 2.0
    return nothing
end

@inline function _taylor_lbend_pass!(r::AbstractVector{S}, le, grd, ba, bye,
    entrance_angle, exit_angle, fint1, fint2, gap,
    T1::AbstractVector, T2::AbstractVector,
    R1::AbstractMatrix, R2::AbstractMatrix) where {S}

    irho = ba / le
    !_taylor_allzero(T1) && _taylor_addvv!(r, T1)
    !_taylor_allzero(R1) && _taylor_multmv!(r, R1)
    edge_fringe_entrance!(r, irho, entrance_angle, fint1, gap, 1)
    _taylor_bend6!(r, le, ba, grd, bye)
    edge_fringe_exit!(r, irho, exit_angle, fint2, gap, 1)
    !_taylor_allzero(R2) && _taylor_multmv!(r, R2)
    !_taylor_allzero(T2) && _taylor_addvv!(r, T2)
    return nothing
end

function QuadFringePassP!(r::AbstractVector{S}, b2) where {S<:AbstractHOTPSA}
    u = b2 / (12.0 * (1.0 + r[6]))
    x2 = r[1]^2
    z2 = r[3]^2
    xz = r[1] * r[3]
    gx = u * (x2 + 3.0 * z2) * r[1]
    gz = u * (z2 + 3.0 * x2) * r[3]
    r1tmp = zero(r[1])
    r3tmp = zero(r[3])

    r[1] += gx
    r1tmp = 3.0 * u * (2.0 * xz * r[4] - (x2 + z2) * r[2])
    r[3] -= gz
    r3tmp = 3.0 * u * (2.0 * xz * r[2] - (x2 + z2) * r[4])
    r[5] -= (gz * r[4] - gx * r[2]) / (1.0 + r[6])
    r[2] += r1tmp
    r[4] -= r3tmp
    return nothing
end

function QuadFringePassN!(r::AbstractVector{S}, b2) where {S<:AbstractHOTPSA}
    u = b2 / (12.0 * (1.0 + r[6]))
    x2 = r[1]^2
    z2 = r[3]^2
    xz = r[1] * r[3]
    gx = u * (x2 + 3.0 * z2) * r[1]
    gz = u * (z2 + 3.0 * x2) * r[3]
    r1tmp = zero(r[1])
    r3tmp = zero(r[3])

    r[1] -= gx
    r1tmp = 3.0 * u * (2.0 * xz * r[4] - (x2 + z2) * r[2])
    r[3] += gz
    r3tmp = 3.0 * u * (2.0 * xz * r[2] - (x2 + z2) * r[4])
    r[5] += (gz * r[4] - gx * r[2]) / (1.0 + r[6])
    r[2] -= r1tmp
    r[4] += r3tmp
    return nothing
end

function linearQuadFringeElegantEntrance!(r6::AbstractVector{S}, b2, fringeIntM0, fringeIntP0) where {S<:AbstractHOTPSA}
    RT = promote_type(S, typeof(b2))
    R = zeros(RT, 6, 6)
    inFringe = -1.0
    fringeIntM = fringeIntP0
    fringeIntP = fringeIntM0
    delta = r6[6]

    quadPartialFringeMatrix!(R, b2 / (1 + delta), inFringe, fringeIntM, 1)
    r6[1] = R[1, 1] * r6[1] + R[1, 2] * r6[2]
    r6[2] = R[2, 1] * r6[1] + R[2, 2] * r6[2]
    r6[3] = R[3, 3] * r6[3] + R[3, 4] * r6[4]
    r6[4] = R[4, 3] * r6[3] + R[4, 4] * r6[4]

    QuadFringePassP!(r6, b2)

    quadPartialFringeMatrix!(R, b2 / (1 + delta), inFringe, fringeIntP, 2)
    r6[1] = R[1, 1] * r6[1] + R[1, 2] * r6[2]
    r6[2] = R[2, 1] * r6[1] + R[2, 2] * r6[2]
    r6[3] = R[3, 3] * r6[3] + R[3, 4] * r6[4]
    r6[4] = R[4, 3] * r6[3] + R[4, 4] * r6[4]
    return nothing
end

function linearQuadFringeElegantExit!(r6::AbstractVector{S}, b2, fringeIntM0, fringeIntP0) where {S<:AbstractHOTPSA}
    RT = promote_type(S, typeof(b2))
    R = zeros(RT, 6, 6)
    inFringe = 1.0
    fringeIntM = fringeIntM0
    fringeIntP = fringeIntP0
    delta = r6[6]

    quadPartialFringeMatrix!(R, b2 / (1 + delta), inFringe, fringeIntM, 1)
    r6[1] = R[1, 1] * r6[1] + R[1, 2] * r6[2]
    r6[2] = R[2, 1] * r6[1] + R[2, 2] * r6[2]
    r6[3] = R[3, 3] * r6[3] + R[3, 4] * r6[4]
    r6[4] = R[4, 3] * r6[3] + R[4, 4] * r6[4]

    QuadFringePassN!(r6, b2)

    quadPartialFringeMatrix!(R, b2 / (1 + delta), inFringe, fringeIntP, 2)
    r6[1] = R[1, 1] * r6[1] + R[1, 2] * r6[2]
    r6[2] = R[2, 1] * r6[1] + R[2, 2] * r6[2]
    r6[3] = R[3, 3] * r6[3] + R[3, 4] * r6[4]
    r6[4] = R[4, 3] * r6[3] + R[4, 4] * r6[4]
    return nothing
end

function Yrot!(r::AbstractVector{S}, phi::Number, beti::Float64) where {S<:AbstractHOTPSA}
    if phi != 0.0
        dp1 = 1.0 * beti + r[6]
        c = cos(phi)
        s = sin(phi)
        pz = pxyz(dp1, r[2], r[4])
        p = c * pz - s * r[2]
        px = s * pz + c * r[2]
        x = r[1] * pz / p
        dy = r[1] * r[4] * s / p
        dct = dp1 * r[1] * s / p
        r[1] = x
        r[2] = px
        r[3] += dy
        r[5] += dct
    end
    return nothing
end

function bend_fringe!(r::AbstractVector{S}, irho::Float64, gK, beti::Float64) where {S<:AbstractHOTPSA}
    b0 = irho
    pz = pxyz(1.0 * beti + r[6], r[2], r[4])
    px = r[2]
    py = r[4]
    d = r[6]
    xp = px / pz
    yp = py / pz
    phi = -b0 * tan(b0 * gK * (1.0 + xp * xp * (2.0 + yp * yp)) * pz - atan(xp / (1.0 + yp * yp)))
    px2 = px * px
    px4 = px2 * px2
    py2 = py * py
    py4 = py2 * py2
    py6 = py4 * py2
    pz2 = pz * pz
    pz3 = pz2 * pz
    pz4 = pz2 * pz2
    pz5 = pz4 * pz
    pz6 = pz4 * pz2
    py2z2 = (py2 + pz2) * (py2 + pz2)
    powsec = (1.0 / cos((b0 * gK * (pz4 + px2 * (py2 + 2.0 * pz2))) / pz3 - atan((px * pz) / (py2 + pz2))))^2
    denom = pz5 * (py4 + px2 * pz2 + 2.0 * py2 * pz2 + pz4)

    dpx = -(b0 * (px2 * pz4 * (py2 - pz2) - pz6 * (py2 + pz2) +
        b0 * gK * px * (pz2 * py2z2 * (2.0 * py2 + 3.0 * pz2) + px4 * (3.0 * py2 * pz2 + 2.0 * pz4) +
        px2 * (3.0 * py6 + 8.0 * py4 * pz2 + 9.0 * py2 * pz4 + 5.0 * pz6))) * powsec) / denom
    dpy = -(b0 * py * (px * pz4 * (py2 + pz2) +
        b0 * gK * (-(pz4 * py2z2) + px4 * (3.0 * py2 * pz2 + 4.0 * pz4) +
        px2 * (3.0 * py6 + 10.0 * py4 * pz2 + 11.0 * py2 * pz4 + 3.0 * pz6))) * powsec) / denom
    dd = (b0 * (1.0 * beti + d) * (px * pz4 * (py2 - pz2) + b0 * gK *
        (-(pz4 * py2z2) + px4 * (3.0 * py2 * pz2 + 2.0 * pz4) +
        px2 * (3.0 * py6 + 8.0 * py4 * pz2 + 7.0 * py2 * pz4 + pz6))) * powsec) / denom

    yf = (2.0 * r[3]) / (1.0 + sqrt(1.0 - 2.0 * dpy * r[3]))
    dxf = 0.5 * dpx * yf * yf
    dct = 0.5 * dd * yf * yf
    dpyf = phi * yf

    r[3] = yf
    r[1] += dxf
    r[4] -= dpyf
    r[5] -= dct
    return nothing
end

function multipole_fringe!(r6::AbstractVector{S}, le::Float64, polya::AbstractVector, polyb::AbstractVector,
    max_order::Int, edge::Float64, skip_b0::Int, beti::Float64) where {S<:AbstractHOTPSA}
    FX = zero(r6[1])
    FY = zero(r6[1])
    FX_X = zero(r6[1])
    FX_Y = zero(r6[1])
    FY_X = zero(r6[1])
    FY_Y = zero(r6[1])
    RX = one(r6[1])
    IX = zero(r6[1])

    for n in 0:max_order
        B = polyb[n + 1]
        A = polya[n + 1]
        j = n + 1.0
        DRX = RX
        DIX = IX

        RX = DRX * r6[1] - DIX * r6[3]
        IX = DRX * r6[3] + DIX * r6[1]

        U = zero(r6[1])
        V = zero(r6[1])
        DU = zero(r6[1])
        DV = zero(r6[1])
        if n == 0 && skip_b0 != 0
            U -= A * IX
            V += A * RX
            DU -= A * DIX
            DV += A * DRX
        else
            U += B * RX - A * IX
            V += B * IX + A * RX
            DU += B * DRX - A * DIX
            DV += B * DIX + A * DRX
        end

        f1 = -edge / 4.0 / (j + 1.0)
        U *= f1
        V *= f1
        DU *= f1
        DV *= f1

        DUX = j * DU
        DVX = j * DV
        DUY = -j * DV
        DVY = j * DU
        nf = (j + 2.0) / j

        FX += U * r6[1] + nf * V * r6[3]
        FY += U * r6[3] - nf * V * r6[1]
        FX_X += DUX * r6[1] + U + nf * r6[3] * DVX
        FX_Y += DUY * r6[1] + nf * V + nf * r6[3] * DVY
        FY_X += DUX * r6[3] - nf * V - nf * r6[1] * DVX
        FY_Y += DUY * r6[3] + U - nf * r6[1] * DVY
    end

    DEL = 1.0 / (1.0 * beti + r6[6])
    A = 1.0 - FX_X * DEL
    B = -FY_X * DEL
    D = 1.0 - FY_Y * DEL
    C = -FX_Y * DEL

    r6[1] -= FX * DEL
    r6[3] -= FY * DEL
    pxf = (D * r6[2] - B * r6[4]) / (A * D - B * C)
    pyf = (A * r6[4] - C * r6[2]) / (A * D - B * C)
    r6[4] = pyf
    r6[2] = pxf
    r6[5] -= (r6[2] * FX + r6[4] * FY) * DEL * DEL
    return nothing
end

function exact_bend!(r6::AbstractVector{S}, irho::Float64, L::Float64, beti::Float64) where {S<:AbstractHOTPSA}
    dp1 = 1.0 * beti + r6[6]
    pz = pxyz(dp1, r6[2], r6[4])

    if abs(irho) < 1e-6
        NormL = L / pz
        r6[1] += r6[2] * NormL
        r6[3] += r6[4] * NormL
        r6[5] += NormL * dp1
    else
        pzmx = pz - (1.0 + r6[1] * irho)
        cs = cos(irho * L)
        sn = sin(irho * L)
        px = r6[2] * cs + pzmx * sn
        d2 = pxyz(dp1, 0.0, r6[4])
        dasin = L + (asin(r6[2] / d2) - asin(px / d2)) / irho
        x = (pxyz(dp1, px, r6[4]) - pzmx * cs + r6[2] * sn - 1.0) / irho
        dy = r6[4] * dasin
        dct = dp1 * dasin

        r6[1] = x
        r6[2] = px
        r6[3] += dy
        r6[5] += dct
    end
    return nothing
end

function bend_edge!(r6::AbstractVector{S}, rhoinv::Float64, theta::Float64, beti::Float64) where {S<:AbstractHOTPSA}
    if abs(rhoinv) >= 1e-6
        dp1 = 1.0 * beti + r6[6]
        c = cos(theta)
        s = sin(theta)
        pz = pxyz(dp1, r6[2], r6[4])
        d2 = pxyz(dp1, 0.0, r6[4])
        px = r6[2] * c + (pz - rhoinv * r6[1]) * s
        dasin = asin(r6[2] / d2) - asin(px / d2)
        num = r6[1] * (r6[2] * sin(2.0 * theta) + s * s * (2.0 * pz - rhoinv * r6[1]))
        den = pxyz(dp1, px, r6[4]) + pxyz(dp1, r6[2], r6[4]) * c - r6[2] * s
        x = r6[1] * c + num / den
        dy = r6[4] * theta / rhoinv + r6[4] / rhoinv * dasin
        dct = dp1 / rhoinv * (theta + dasin)

        r6[1] = x
        r6[2] = px
        r6[3] += dy
        r6[5] += dct
    end
    return nothing
end

function bndthinkick!(r::AbstractVector{S}, A::AbstractVector, B::AbstractVector,
    L::Float64, irho::Float64, max_order::Int, beti::Float64) where {S<:AbstractHOTPSA}
    ReSum = B[max_order + 1]
    ImSum = A[max_order + 1]
    ReSumTemp = zero(ReSum)

    for i in max_order-1:-1:0
        ReSumTemp = ReSum * r[1] - ImSum * r[3] + B[i + 1]
        ImSum = ImSum * r[1] + ReSum * r[3] + A[i + 1]
        ReSum = ReSumTemp
    end

    r[2] -= L * (ReSum - (r[6] - r[1] * irho) * irho)
    r[4] += L * ImSum
    r[5] += L * (irho * r[1]) * beti
    return nothing
end

function BendSymplecticPass!(r::AbstractVector{S}, le::Float64, beti::Float64, irho::Float64,
    A::AbstractVector, B::AbstractVector, max_order::Int, num_int_steps::Int,
    entrance_angle::Float64, exit_angle::Float64,
    FringeBendEntrance::Int, FringeBendExit::Int,
    fint1, fint2, gap,
    FringeQuadEntrance::Int, FringeQuadExit::Int,
    fringeIntM0, fringeIntP0,
    T1::AbstractVector, T2::AbstractVector, R1::AbstractMatrix, R2::AbstractMatrix,
    RApertures, EApertures, KickAngle::AbstractVector) where {S<:AbstractHOTPSA}
    DRIFT1 = 0.6756035959798286638
    DRIFT2 = -0.1756035959798286639
    KICK1 = 1.351207191959657328
    KICK2 = -1.702414383919314656
    SL = le / num_int_steps
    L1 = SL * DRIFT1
    L2 = SL * DRIFT2
    K1 = SL * KICK1
    K2 = SL * KICK2

    useLinFrEleEntrance = FringeQuadEntrance == 2 ? 1 : 0
    useLinFrEleExit = FringeQuadExit == 2 ? 1 : 0

    B[1] -= sin(KickAngle[1]) / le
    A[1] += sin(KickAngle[2]) / le

    if !iszero(T1)
        _taylor_addvv!(r, T1)
    end
    if !iszero(R1)
        _taylor_multmv!(r, R1)
    end

    edge_fringe_entrance!(r, irho, entrance_angle, fint1, gap, FringeBendEntrance)
    if !iszero(FringeQuadEntrance) && !iszero(B[2])
        if useLinFrEleEntrance == 1
            linearQuadFringeElegantEntrance!(r, B[2], fringeIntM0, fringeIntP0)
        else
            QuadFringePassP!(r, B[2])
        end
    end

    for _ in 1:num_int_steps
        _taylor_drift6!(r, L1, beti)
        bndthinkick!(r, A, B, K1, irho, max_order, beti)
        _taylor_drift6!(r, L2, beti)
        bndthinkick!(r, A, B, K2, irho, max_order, beti)
        _taylor_drift6!(r, L2, beti)
        bndthinkick!(r, A, B, K1, irho, max_order, beti)
        _taylor_drift6!(r, L1, beti)
    end

    if !iszero(FringeQuadExit) && !iszero(B[2])
        if useLinFrEleExit == 1
            linearQuadFringeElegantExit!(r, B[2], fringeIntM0, fringeIntP0)
        else
            QuadFringePassN!(r, B[2])
        end
    end
    edge_fringe_exit!(r, irho, exit_angle, fint2, gap, FringeBendExit)

    if !iszero(R2)
        _taylor_multmv!(r, R2)
    end
    if !iszero(T2)
        _taylor_addvv!(r, T2)
    end

    B[1] += sin(KickAngle[1]) / le
    A[1] -= sin(KickAngle[2]) / le
    return nothing
end

function ExactSectorBend!(r6::AbstractVector{S}, le::Float64, beti::Float64, angle::Float64,
    A::AbstractVector, B::AbstractVector, max_order::Int, num_int_steps::Int,
    entrance_angle::Float64, exit_angle::Float64, FringeBendEntrance::Int, FringeBendExit::Int,
    FringeQuadEntrance::Int, FringeQuadExit::Int, gk,
    T1::AbstractVector, T2::AbstractVector, R1::AbstractMatrix, R2::AbstractMatrix,
    RApertures, EApertures, KickAngle::AbstractVector) where {S<:AbstractHOTPSA}
    irho = angle / le
    DRIFT1 = 0.6756035959798286638
    DRIFT2 = -0.1756035959798286639
    KICK1 = 1.351207191959657328
    KICK2 = -1.702414383919314656
    SL = le / num_int_steps
    L1 = SL * DRIFT1
    L2 = SL * DRIFT2
    K1 = SL * KICK1
    K2 = SL * KICK2

    if !iszero(KickAngle[1])
        B[1] -= sin(KickAngle[1]) / le
    end
    if !iszero(KickAngle[2])
        A[1] += sin(KickAngle[2]) / le
    end

    if !iszero(T1)
        _taylor_addvv!(r6, T1)
    end
    if !iszero(R1)
        _taylor_multmv!(r6, R1)
    end

    Yrot!(r6, entrance_angle, beti)
    if FringeBendEntrance != 0
        bend_fringe!(r6, irho, gk, beti)
    end
    if FringeQuadEntrance != 0
        multipole_fringe!(r6, le, A, B, max_order, 1.0, 1, beti)
    end
    bend_edge!(r6, irho, -entrance_angle, beti)

    if num_int_steps == 0
        exact_bend!(r6, irho, le, beti)
    else
        for _ in 1:num_int_steps
            exact_bend!(r6, irho, L1, beti)
            _taylor_strthinkick!(r6, A, B, K1, max_order)
            exact_bend!(r6, irho, L2, beti)
            _taylor_strthinkick!(r6, A, B, K2, max_order)
            exact_bend!(r6, irho, L2, beti)
            _taylor_strthinkick!(r6, A, B, K1, max_order)
            exact_bend!(r6, irho, L1, beti)
        end
    end

    r6[5] -= le * beti
    bend_edge!(r6, irho, -exit_angle, beti)
    if FringeQuadExit != 0
        multipole_fringe!(r6, le, A, B, max_order, -1.0, 1, beti)
    end
    if FringeBendExit != 0
        bend_fringe!(r6, -irho, gk, beti)
    end
    Yrot!(r6, exit_angle, beti)

    if !iszero(R2)
        _taylor_multmv!(r6, R2)
    end
    if !iszero(T2)
        _taylor_addvv!(r6, T2)
    end

    if !iszero(KickAngle[1])
        B[1] += sin(KickAngle[1]) / le
    end
    if !iszero(KickAngle[2])
        A[1] -= sin(KickAngle[2]) / le
    end
    return nothing
end

@inline function _taylor_reference_kinematics(E0::Number, m0::Number)
    gamma = (E0 + m0) / m0
    beta = sqrt(1.0 - 1.0 / gamma^2)
    beti = use_exact_beti == 1 ? 1.0 / beta : 1.0
    gamma2i = use_exact_Hamiltonian == 2 ? 1.0 / gamma^2 : 0.0
    return beti, gamma2i
end

function pass_Taylor!(ele::DRIFT{ET}, r_in::AbstractVector{S}; E0::Number=0.0, m0::Number=m_e) where {ET,S}
    beti, gamma2i = _taylor_reference_kinematics(E0, m0)
    _taylor_drift_pass!(r_in, ele.len, beti, gamma2i, ele.T1, ele.T2, ele.R1, ele.R2)
    return nothing
end

function pass_Taylor!(::MARKER{ET}, r_in::AbstractVector{S}; E0::Number=0.0, m0::Number=m_e) where {ET,S}
    return nothing
end

function pass_Taylor!(ele::CORRECTOR{ET}, r_in::AbstractVector{S}; E0::Number=0.0, m0::Number=m_e) where {ET,S}
    _taylor_corrector_pass!(r_in, ele.len, ele.xkick, ele.ykick, ele.T1, ele.T2, ele.R1, ele.R2)
    return nothing
end

function pass_Taylor!(ele::QUAD{ET}, r_in::AbstractVector{S}; E0::Number=0.0, m0::Number=m_e) where {ET,S}
    beti, gamma2i = _taylor_reference_kinematics(E0, m0)
    _taylor_quad_linear_pass!(r_in, ele.len, ele.k1, beti, gamma2i, ele.T1, ele.T2, ele.R1, ele.R2)
    return nothing
end

function pass_Taylor!(ele::TRANSLATION{ET}, r_in::AbstractVector{S}; E0::Number=0.0, m0::Number=0.0) where {ET,S}
    beta = 1.0
    if use_exact_beti == 1
        gamma = (E0 + m0) / m0
        beta = sqrt(1.0 - 1.0 / gamma^2)
    end
    _taylor_translation_pass!(r_in, ele.dx, ele.dy, ele.ds, beta)
    return nothing
end

function pass_Taylor!(ele::YROTATION{ET}, r_in::AbstractVector{S}; E0::Number=0.0, m0::Number=0.0) where {ET,S}
    beta = 1.0
    if use_exact_beti == 1
        gamma = (E0 + m0) / m0
        beta = sqrt(1.0 - 1.0 / gamma^2)
    end
    _taylor_yrotation_pass!(r_in, ele.angle, beta)
    return nothing
end

function pass_Taylor!(ele::thinMULTIPOLE{ET}, r_in::AbstractVector{S}; E0::Number=0.0, m0::Number=m_e) where {ET,S}
    _taylor_thinmultipole_pass!(r_in, ele.PolynomA, ele.PolynomB, ele.MaxOrder,
        ele.T1, ele.T2, ele.R1, ele.R2, ele.KickAngle)
    return nothing
end

function pass_Taylor!(ele::KQUAD{ET}, r_in::AbstractVector{S}; E0::Number=0.0, m0::Number=m_e) where {ET,S}
    ele.rad == 0 || throw(ArgumentError("HOTPSA KQUAD radiation is not implemented yet"))
    beti, gamma2i = _taylor_reference_kinematics(E0, m0)
    PolynomB = SVector(ele.k0, ele.k1, ele.k2 / 2, ele.k3 / 6)
    _taylor_multipole_integrator_pass!(r_in, ele.len, beti, gamma2i, ele.PolynomA, PolynomB,
        ele.MaxOrder, ele.NumIntSteps, ele.FringeQuadEntrance, ele.FringeQuadExit,
        ele.T1, ele.T2, ele.R1, ele.R2, ele.KickAngle)
    return nothing
end

function pass_Taylor!(ele::KSEXT{ET}, r_in::AbstractVector{S}; E0::Number=0.0, m0::Number=m_e) where {ET,S}
    ele.rad == 0 || throw(ArgumentError("HOTPSA KSEXT radiation is not implemented yet"))
    beti, gamma2i = _taylor_reference_kinematics(E0, m0)
    PolynomB = SVector(ele.k0, ele.k1, ele.k2 / 2, ele.k3 / 6)
    _taylor_multipole_integrator_pass!(r_in, ele.len, beti, gamma2i, ele.PolynomA, PolynomB,
        ele.MaxOrder, ele.NumIntSteps, ele.FringeQuadEntrance, ele.FringeQuadExit,
        ele.T1, ele.T2, ele.R1, ele.R2, ele.KickAngle)
    return nothing
end

function pass_Taylor!(ele::KOCT{ET}, r_in::AbstractVector{S}; E0::Number=0.0, m0::Number=m_e) where {ET,S}
    ele.rad == 0 || throw(ArgumentError("HOTPSA KOCT radiation is not implemented yet"))
    beti, gamma2i = _taylor_reference_kinematics(E0, m0)
    PolynomB = SVector(ele.k0, ele.k1, ele.k2 / 2, ele.k3 / 6)
    _taylor_multipole_integrator_pass!(r_in, ele.len, beti, gamma2i, ele.PolynomA, PolynomB,
        ele.MaxOrder, ele.NumIntSteps, ele.FringeQuadEntrance, ele.FringeQuadExit,
        ele.T1, ele.T2, ele.R1, ele.R2, ele.KickAngle)
    return nothing
end

function pass_Taylor!(ele::SOLENOID{ET}, r_in::AbstractVector{S}; E0::Number=0.0, m0::Number=m_e) where {ET,S}
    beti, _ = _taylor_reference_kinematics(E0, m0)
    _taylor_solenoid_pass!(r_in, ele.len, ele.ks, beti, ele.T1, ele.T2, ele.R1, ele.R2)
    return nothing
end

function pass_Taylor!(ele::RFCA{ET}, r_in::AbstractVector{S}; E0::Number=0.0, m0::Number=m_e) where {ET,S}
    gamma = (E0 + m0) / m0
    beta = sqrt(1.0 - 1.0 / gamma^2)
    beti = use_exact_beti == 1 ? 1.0 / beta : 1.0
    nv = ele.volt / ele.energy
    _taylor_rfca_pass!(r_in, ele.len, nv, ele.freq, ele.lag, ele.philag, beta, beti)
    return nothing
end

function pass_Taylor!(ele::CRABCAVITY{ET}, r_in::AbstractVector{S}; E0::Number=0.0, m0::Number=m_e) where {ET,S}
    gamma = (E0 + m0) / m0
    beta = sqrt(1.0 - 1.0 / gamma^2)
    beti = use_exact_beti == 1 ? 1.0 / beta : 1.0
    _taylor_crabcavity_pass!(r_in, ele.len, ele.volt, ele.k, ele.phi, E0, beta, beti)
    return nothing
end

function pass_Taylor!(ele::CRABCAVITY_K2{ET}, r_in::AbstractVector{S}; E0::Number=0.0, m0::Number=m_e) where {ET,S}
    gamma = (E0 + m0) / m0
    beta = sqrt(1.0 - 1.0 / gamma^2)
    beti = use_exact_beti == 1 ? 1.0 / beta : 1.0
    _taylor_crabcavity_k2_pass!(r_in, ele.len, ele.volt, ele.k, ele.phi, ele.k2, E0, beta, beti)
    return nothing
end

function pass_Taylor!(ele::easyCRABCAVITY{ET}, r_in::AbstractVector{S}; E0::Number=0.0, m0::Number=m_e) where {ET,S}
    _taylor_easycrabcavity_pass!(r_in, ele.halfthetac, ele.k, ele.phi)
    return nothing
end

function pass_Taylor!(ele::AccelCavity{ET}, r_in::AbstractVector{S}; E0::Number=0.0, m0::Number=m_e) where {ET,S}
    gamma = (E0 + m0) / m0
    beta = sqrt(1.0 - 1.0 / gamma^2)
    beta2E = beta^2 * E0
    _taylor_accelcavity_pass!(r_in, ele.volt, ele.k, ele.phis, beta2E)
    return nothing
end

function pass_Taylor!(ele::LorentzBoost{ET}, r_in::AbstractVector{S}; E0::Number=0.0, m0::Number=m_e) where {ET,S}
    _taylor_lorentz_boost_pass!(r_in, ele.cosang, ele.tanang, ele.mode)
    return nothing
end

function pass_Taylor!(ele::InvLorentzBoost{ET}, r_in::AbstractVector{S}; E0::Number=0.0, m0::Number=m_e) where {ET,S}
    _taylor_inv_lorentz_boost_pass!(r_in, ele.sinang, ele.cosang, ele.mode)
    return nothing
end

function pass_Taylor!(ele::LBEND{ET}, r_in::AbstractVector{S}; E0::Number=0.0, m0::Number=m_e) where {ET,S}
    _taylor_lbend_pass!(r_in, ele.len, ele.K, ele.angle, ele.ByError,
        ele.e1, ele.e2, ele.fint1, ele.fint2, ele.FullGap,
        ele.T1, ele.T2, ele.R1, ele.R2)
    return nothing
end

function pass_Taylor!(ele::ESBEND{ET}, r_in::AbstractVector{S}; E0::Number=0.0, m0::Number=m_e) where {ET,S}
    ele.rad == 0 || throw(ArgumentError("HOTPSA ESBEND radiation is not implemented yet"))
    beti, _ = _taylor_reference_kinematics(E0, m0)
    ExactSectorBend!(r_in, ele.len, beti, ele.angle, ele.PolynomA, ele.PolynomB, ele.MaxOrder, ele.NumIntSteps,
        ele.e1, ele.e2, ele.FringeBendEntrance, ele.FringeBendExit,
        ele.FringeQuadEntrance, ele.FringeQuadExit, ele.gK,
        ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle)
    return nothing
end

function pass_Taylor!(ele::SBEND{ET}, r_in::AbstractVector{S}; E0::Number=0.0, m0::Number=m_e) where {ET,S}
    ele.rad == 0 || throw(ArgumentError("HOTPSA SBEND radiation is not implemented yet"))
    beti, _ = _taylor_reference_kinematics(E0, m0)
    irho = ele.angle / ele.len
    BendSymplecticPass!(r_in, ele.len, beti, irho, ele.PolynomA, ele.PolynomB, ele.MaxOrder, ele.NumIntSteps,
        ele.e1, ele.e2, ele.FringeBendEntrance, ele.FringeBendExit,
        ele.fint1, ele.fint2, ele.gap, ele.FringeQuadEntrance, ele.FringeQuadExit,
        ele.FringeIntM0, ele.FringeIntP0,
        ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle)
    return nothing
end

function linepass_Taylor!(line::Vector{<:AbstractElement}, rin::AbstractVector{S}; E0::Number=3e9, m0::Number=m_e) where {S}
    length(rin) == 6 || error("Taylor state vectors must have length 6")
    for ele in line
        pass_Taylor!(ele, rin, E0=E0, m0=m0)
    end
    return nothing
end

function ADlinepass_Taylor!(line::Vector{<:AbstractElement}, rin::AbstractVector{S},
    changed_idx::Vector, changed_ele::Vector; E0::Number=3e9, m0::Number=m_e) where {S}
    length(rin) == 6 || error("Taylor state vectors must have length 6")
    for i in eachindex(line)
        idx_in_changed = findfirst(==(i), changed_idx)
        if idx_in_changed !== nothing
            pass_Taylor!(changed_ele[idx_in_changed], rin, E0=E0, m0=m0)
        else
            pass_Taylor!(line[i], rin, E0=E0, m0=m0)
        end
    end
    return nothing
end

function ringpass_Taylor!(line::Vector{<:AbstractElement}, rin::AbstractVector{S}, nturn::Int;
    E0::Number=3e9, m0::Number=m_e) where {S}
    for _ in 1:nturn
        linepass_Taylor!(line, rin, E0=E0, m0=m0)
    end
    return nothing
end

function pass!(ele::AbstractElement, r_in::AbstractVector{S}; E0::Number=0.0, m0::Number=m_e) where {S<:AbstractHOTPSA}
    pass_Taylor!(ele, r_in, E0=E0, m0=m0)
    return nothing
end

function pass!(ele::AbstractElement, r_in::AbstractMatrix{S}, num_particles::Integer, particles::Beam{S}) where {S<:AbstractHOTPSA}
    E0 = Float64(jet_primal(particles.energy))
    m0 = Float64(jet_primal(particles.mass))
    @inbounds for c in 1:num_particles
        particles.lost_flag[c] == 1 && continue
        pass_Taylor!(ele, @view(r_in[c, :]), E0=E0, m0=m0)
    end
    return nothing
end

function linepass!(line::Vector{<:AbstractElement}, rin::AbstractVector{S};
    E0::Number=3e9, m0::Number=m_e) where {S<:AbstractHOTPSA}
    linepass_Taylor!(line, rin, E0=E0, m0=m0)
    return nothing
end

function linepass!(line::Vector{<:AbstractElement}, particles::Beam{S}) where {S<:AbstractHOTPSA}
    np = particles.nmacro
    for i in eachindex(line)
        pass!(line[i], particles.r, np, particles)
    end
    return nothing
end

function linepass!(line::Vector{<:AbstractElement}, particles::Beam{S}, refpts::Vector) where {S<:AbstractHOTPSA}
    np = particles.nmacro
    saved_particles = []
    for i in eachindex(line)
        pass!(line[i], particles.r, np, particles)
        if i in refpts
            push!(saved_particles, copy(particles.r))
        end
    end
    return saved_particles
end

function ADlinepass!(line::Vector{<:AbstractElement}, rin::AbstractVector{S},
    changed_idx::Vector, changed_ele::Vector; E0::Number=3e9, m0::Number=m_e) where {S<:AbstractHOTPSA}
    ADlinepass_Taylor!(line, rin, changed_idx, changed_ele, E0=E0, m0=m0)
    return nothing
end

function ADlinepass!(line::Vector{<:AbstractElement}, particles::Beam{S},
    changed_idx::Vector, changed_ele::Vector) where {S<:AbstractHOTPSA}
    np = particles.nmacro
    for i in eachindex(line)
        idx_in_changed = findfirst(==(i), changed_idx)
        if idx_in_changed !== nothing
            pass!(changed_ele[idx_in_changed], particles.r, np, particles)
        else
            pass!(line[i], particles.r, np, particles)
        end
    end
    return nothing
end

function ADlinepass!(line::Vector{<:AbstractElement}, particles::Beam{S}, refpts::Vector,
    changed_idx::Vector, changed_ele::Vector) where {S<:AbstractHOTPSA}
    np = particles.nmacro
    saved_particles = []
    for i in eachindex(line)
        idx_in_changed = findfirst(==(i), changed_idx)
        if idx_in_changed !== nothing
            pass!(changed_ele[idx_in_changed], particles.r, np, particles)
        else
            pass!(line[i], particles.r, np, particles)
        end
        if i in refpts
            push!(saved_particles, copy(particles.r))
        end
    end
    return saved_particles
end

function ringpass!(line::Vector{<:AbstractElement}, rin::AbstractVector{S}, nturn::Int;
    E0::Number=3e9, m0::Number=m_e) where {S<:AbstractHOTPSA}
    ringpass_Taylor!(line, rin, nturn, E0=E0, m0=m0)
    return nothing
end

function ringpass!(line::Vector{<:AbstractElement}, particles::Beam{S}, nturn::Int) where {S<:AbstractHOTPSA}
    for _ in 1:nturn
        linepass!(line, particles)
    end
    return nothing
end

function ringpass!(line::Vector{<:AbstractElement}, particles::Beam{S}, nturn::Int, save::Bool) where {S<:AbstractHOTPSA}
    save_beam = []
    for _ in 1:nturn
        linepass!(line, particles)
        if save
            push!(save_beam, copy(particles.r))
        end
    end
    return save_beam
end
