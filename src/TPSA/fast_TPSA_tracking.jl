# use first-order TPSA. It is faster than other TPSA modules.
using .TPSAadStatic

function strthinkick!(r::SubArray, A::Vector{DTPSAD{N, T}}, B::Vector{DTPSAD{N, T}}, L::DTPSAD{N, T}, max_order::Int) where {N, T <: Number}
    # Modified based on AT function. Ref[Terebilo, Andrei. "Accelerator modeling with MATLAB accelerator toolbox." PACS2001 (2001)].

    ReSum = B[max_order + 1]
    ImSum = A[max_order + 1]
    ReSumTemp = 0.0

    for i in max_order-1: -1: 0
        ReSumTemp = ReSum * r[1] - ImSum * r[3] + B[i+1]
        ImSum = ImSum * r[1] + ReSum * r[3] + A[i+1]
        ReSum = ReSumTemp
    end

    r[2] -= L * ReSum
    r[4] += L * ImSum
    return nothing
end

function strthinkick1!(r::SubArray, A::Vector{DTPSAD{N, T}}, B::Vector{DTPSAD{N, T}}, L::DTPSAD{N, T}, max_order::Int) where {N, T <: Number}
    # Modified based on AT function. Ref[Terebilo, Andrei. "Accelerator modeling with MATLAB accelerator toolbox." PACS2001 (2001)].

    # Calculate and apply a multipole kick to a 6-dimentional
    # phase space vector in a straight element (quadrupole)
    
    # IMPORTANT !!!
    # The reference coordinate system is straight but the field expansion may still
    # contain dipole terms A[1], B[1]

    ReSum = B[max_order + 1]
    ImSum = A[max_order + 1]
    ReSumTemp = 0.0

    for i in reverse(1:max_order)
        ReSumTemp = ReSum * r[1] - ImSum * r[3] + B[i]
        ImSum = ImSum * r[1] + ReSum * r[3] + A[i]
        ReSum = ReSumTemp
    end

    r[2] -= L * ReSum
    r[4] += L * ImSum
    return nothing
end

function bndthinkick!(r::SubArray, A::Array{DTPSAD{N, T},1}, B::Array{DTPSAD{N, T},1}, 
    L::DTPSAD{N, T}, irho::DTPSAD{N, T}, max_order::Int, beti::Float64) where {N, T <: Number}
    # AT function. Ref[Terebilo, Andrei. "Accelerator modeling with MATLAB accelerator toolbox." PACS2001 (2001)].
    ReSum = B[max_order + 1]
    ImSum = A[max_order + 1]
    ReSumTemp = 0.0

    for i in max_order-1:-1:0
        ReSumTemp = ReSum * r[1] - ImSum * r[3] + B[i+1]
        ImSum = ImSum * r[1] + ReSum * r[3] + A[i+1]
        ReSum = ReSumTemp
    end

    r[2] -= L * (ReSum - (r[6] - r[1] * irho) * irho)
    r[4] += L * ImSum
    r[5] += L * irho * r[1] * beti # Path length
    return nothing
end

function bndthinkickrad!(r::SubArray, A::Array{DTPSAD{N, T},1}, B::Array{DTPSAD{N, T},1}, 
    L::DTPSAD{N, T}, irho::DTPSAD{N, T}, E0::DTPSAD{N, T}, max_order::Int, beti::Float64) where {N, T <: Number}
    # Modified based on AT function. Ref[Terebilo, Andrei. "Accelerator modeling with MATLAB accelerator toolbox." PACS2001 (2001)].
    ReSum = B[max_order + 1]
    ImSum = A[max_order + 1]
    ReSumTemp = 0.0
    CRAD = CGAMMA * E0^3 / (2.0*pi*1e27) # [m]/[GeV^3] M.Sands (4.1)

    for i in max_order-1:-1:0
        ReSumTemp = ReSum * r[1] - ImSum * r[3] + B[i+1]
        ImSum = ImSum * r[1] + ReSum * r[3] + A[i+1]
        ReSum = ReSumTemp
    end

    # angles from momentums
    p_norm = 1.0 / (1.0*beti + r[6])
    x = r[1]
    xpr = r[2] * p_norm
    y = r[3]
    ypr = r[4] * p_norm
    B2P = B2perp(ImSum, ReSum + irho, irho, x, xpr, y, ypr)

    dp_0 = r[6]
    r[6] = r[6] - CRAD * (1.0*beti+r[6])^2 * B2P * (1.0 + x*irho + (xpr^2 + ypr^2) / 2.0) * L
    
    # momentums after losing energy
    p_norm = 1.0 / (1.0*beti + r[6])
    r[2] = xpr / p_norm
    r[4] = ypr / p_norm

    r[2] -= L * (ReSum - (dp_0 - r[1]*irho)*irho)
    r[4] += L * ImSum
    r[5] += L * irho * r[1] * beti # Path length
    return nothing
end

function Yrot!(r::SubArray, phi::DTPSAD{N, T}, beti::Float64) where {N, T <: Number}
    # Forest 10.26, rotation in free space
    if phi != 0.0
        dp1 = 1.0*beti + r[6]
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

function Sec(x::DTPSAD{N, T}) where {N, T <: Number}
    return 1.0 / cos(x)
end

function bend_fringe!(r::SubArray, irho::DTPSAD{N, T}, gK::DTPSAD{N, T}, beti::Float64) where {N, T <: Number}
    # Forest 13.13, bend fringe in the hard-edge limit
    b0 = irho 
    pz = pxyz(1.0*beti+r[6], r[2], r[4])
    px = r[2]
    py = r[4]
    d = r[6]
    xp = px / pz
    yp = py / pz
    phi = -b0 * tan( b0 * gK * (1.0 + xp*xp*(2.0 + yp*yp))*pz - atan(xp / (1.0 + yp*yp)))
    px2 = px*px
    px4 = px2*px2
    py2 = py*py
    py4 = py2*py2
    py6 = py4*py2
    pz2 = pz*pz
    pz3 = pz2*pz
    pz4 = pz2*pz2
    pz5 = pz4*pz
    pz6 = pz4*pz2
    py2z2 = (py2 + pz2) * (py2 + pz2)
    # powsec = pow(Sec((b0*gK*(pz4 + px2*(py2 + 2*pz2)))/pz3 - atan((px*pz)/(py2 + pz2))),2)
    powsec = Sec((b0*gK*(pz4 + px2*(py2 + 2.0*pz2)))/pz3 - atan((px*pz)/(py2 + pz2)))^2
    denom = (pz5*(py4 + px2*pz2 + 2.0*py2*pz2 + pz4))

    dpx = -(b0*(px2*pz4*(py2 - pz2) - pz6*(py2 + pz2) +
            b0*gK*px*(pz2*py2z2*(2.0*py2 + 3.0*pz2) + px4*(3.0*py2*pz2 + 2.0*pz4) +
            px2*(3.0*py6 + 8.0*py4*pz2 + 9.0*py2*pz4 + 5.0*pz6)))*powsec)/denom
    dpy = -(b0*py*(px*pz4*(py2 + pz2) +
            b0*gK*(-(pz4*py2z2) + px4*(3.0*py2*pz2 + 4.0*pz4) +
                   px2*(3.0*py6 + 10.0*py4*pz2 + 11.0*py2*pz4 + 3.0*pz6)))*powsec)/denom
    dd = (b0*(1.0*beti + d)*(px*pz4*(py2 - pz2) + b0*gK*
                      (-(pz4*py2z2) + px4*(3.0*py2*pz2 + 2.0*pz4) +
                       px2*(3.0*py6 + 8.0*py4*pz2 + 7.0*py2*pz4 + pz6)))*powsec)/denom

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

function multipole_fringe!(r6::SubArray, le::DTPSAD{N, T}, polya::Array{DTPSAD{N, T},1}, polyb::Array{DTPSAD{N, T},1}, 
        max_order::Int, edge::DTPSAD{N, T}, skip_b0::Int, beti::Float64) where {N, T <: Number}
    # Forest 13.29
    FX = 0.0
    FY = 0.0
    FX_X = 0.0
    FX_Y = 0.0
    FY_X = 0.0
    FY_Y = 0.0
  
    RX = 1.0
    IX = 0.0

    for n in 0:max_order
        B = polyb[n + 1]  
        A = polya[n + 1] 
    
        j = n + 1.0
    
        DRX = RX
        DIX = IX
    
        # Complex multiplications
        RX = DRX * r6[1] - DIX * r6[3]
        IX = DRX * r6[3] + DIX * r6[1]
        
        U, V, DU, DV = 0.0, 0.0, 0.0, 0.0
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
    
        nf = 1.0 * (j + 2.0) / j
    
        FX += U * r6[1] + nf * V * r6[3]
        FY += U * r6[3] - nf * V * r6[1]
    
        FX_X += DUX * r6[1] + U + nf * r6[3] * DVX
        FX_Y += DUY * r6[1] + nf * V + nf * r6[3] * DVY
    
        FY_X += DUX * r6[3] - nf * V - nf * r6[1] * DVX
        FY_Y += DUY * r6[3] + U - nf * r6[1] * DVY
    end

    DEL = 1.0 / (1.0*beti + r6[6])
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

function exact_bend!(r6::SubArray, irho::DTPSAD{N, T}, L::DTPSAD{N, T}, beti::Float64) where {N, T <: Number}
    # Forest 12.18, bend-kick split, map W(L, irho)

    dp1 = 1.0*beti + r6[6]  # r6[delta_]
    pz = pxyz(dp1, r6[2], r6[4])  # r6[px_], r6[py_]

    if abs(irho) < 1e-6
        NormL = L / pz
        r6[1] += r6[2] * NormL  # r6[x_]
        r6[3] += r6[4] * NormL  # r6[y_]
        r6[5] += NormL * dp1    # r6[ct_], absolute path length
    else
        pzmx = pz - (1.0 + r6[1] * irho)  # r6[x_]
        cs = cos(irho * L)
        sn = sin(irho * L)
        px = r6[2] * cs + pzmx * sn  # r6[px_]
        d2 = pxyz(dp1, 0.0, r6[4])  # r6[py_]
        dasin = L + (asin(r6[2] / d2) - asin(px / d2)) / irho
        x = (pxyz(dp1, px, r6[4]) - pzmx * cs + r6[2] * sn - 1.0) / irho  # r6[x_]
        dy = r6[4] * dasin  # r6[py_]
        dct = dp1 * dasin   # r6[ct_], absolute path length

        r6[1] = x
        r6[2] = px
        r6[3] += dy
        r6[5] += dct
    end
    return nothing
end

function exactbndthinkick_rad!(r::SubArray, A::Array{DTPSAD{N, T},1}, B::Array{DTPSAD{N, T},1}, 
    L::DTPSAD{N, T}, irho::DTPSAD{N, T}, max_order::Int, beti::Float64, rad_const::DTPSAD{N, T}) where {N, T <: Number}
    # AT function. Ref[Terebilo, Andrei. "Accelerator modeling with MATLAB accelerator toolbox." PACS2001 (2001)].
    ReSum = B[max_order + 1]
    ImSum = A[max_order + 1]
    ReSumTemp = 0.0

    for i in max_order-1:-1:0
        ReSumTemp = ReSum * r[1] - ImSum * r[3] + B[i+1]
        ImSum = ImSum * r[1] + ReSum * r[3] + A[i+1]
        ReSum = ReSumTemp
    end

    p_norm = 1.0 / (1.0 + r[6])
    x = r[1]
    xpr = r[2] * p_norm
    y = r[3]
    ypr = r[4] * p_norm

    B2P = B2perp_exact_bnd(ImSum, ReSum+irho, irho, x, xpr, y, ypr)

    r[6] -= rad_const * (1.0 + r[6])^2 * B2P * (1.0 + x*irho) * L / sqrt(1.0 - xpr*xpr - ypr*ypr)

    p_norm = 1.0 / (1.0 + r[6])
    r[2] = xpr/p_norm
    r[4] = ypr/p_norm

    r[2] -= L * ReSum
    r[4] += L * ImSum
    return nothing
end

@inline function exactbndthinkick_rad_row!(r::Matrix{DTPSAD{N, T}}, c::Int, A::Array{DTPSAD{N, T},1}, B::Array{DTPSAD{N, T},1},
    L::DTPSAD{N, T}, irho::DTPSAD{N, T}, max_order::Int, beti::Float64, rad_const::DTPSAD{N, T}) where {N, T <: Number}
    ReSum = B[max_order + 1]
    ImSum = A[max_order + 1]
    ReSumTemp = 0.0
    x = r[c, 1]
    y = r[c, 3]

    for i in max_order-1:-1:0
        ReSumTemp = ReSum * x - ImSum * y + B[i+1]
        ImSum = ImSum * x + ReSum * y + A[i+1]
        ReSum = ReSumTemp
    end

    p_norm = 1.0 / (1.0 + r[c, 6])
    xpr = r[c, 2] * p_norm
    ypr = r[c, 4] * p_norm
    B2P = B2perp_exact_bnd(ImSum, ReSum + irho, irho, x, xpr, y, ypr)

    r[c, 6] -= rad_const * (1.0 + r[c, 6])^2 * B2P * (1.0 + x*irho) * L / sqrt(1.0 - xpr*xpr - ypr*ypr)

    p_norm = 1.0 / (1.0 + r[c, 6])
    r[c, 2] = xpr / p_norm
    r[c, 4] = ypr / p_norm

    r[c, 2] -= L * ReSum
    r[c, 4] += L * ImSum
    return nothing
end

function bend_edge!(r6::SubArray, rhoinv::DTPSAD{N, T}, theta::DTPSAD{N, T}, beti::Float64) where {N, T <: Number}
    # Forest 12.41, ideal wedge, map U(theta, rhoinv)

    if abs(rhoinv) >= 1e-6
        dp1 = 1.0*beti + r6[6]  # r6[delta_]
        c = cos(theta)
        s = sin(theta)
        pz = pxyz(dp1, r6[2], r6[4])  # r6[px_], r6[py_]
        d2 = pxyz(dp1, 0.0, r6[4])    # r6[py_]
        px = r6[2] * c + (pz - rhoinv * r6[1]) * s  # r6[px_]
        dasin = asin(r6[2] / d2) - asin(px / d2)
        num = r6[1] * (r6[2] * sin(2.0 * theta) + s * s * (2.0 * pz - rhoinv * r6[1]))
        den = pxyz(dp1, px, r6[4]) + pxyz(dp1, r6[2], r6[4]) * c - r6[2] * s
        x = r6[1] * c + num / den  # r6[x_]
        dy = r6[4] * theta / rhoinv + r6[4] / rhoinv * dasin  # r6[py_]
        dct = dp1 / rhoinv * (theta + dasin)  # r6[ct_]

        r6[1] = x
        r6[2] = px
        r6[3] += dy
        r6[5] += dct
    end
    return nothing
end

function ExactSectorBend!(r::Matrix{DTPSAD{N, T}}, le::DTPSAD{N, T}, beti::Float64, angle::DTPSAD{N, T}, A::Array{DTPSAD{N, T},1}, B::Array{DTPSAD{N, T},1},
    max_order::Int, num_int_steps::Int, entrance_angle::DTPSAD{N, T}, exit_angle::DTPSAD{N, T}, FringeBendEntrance::Int, FringeBendExit::Int,
    FringeQuadEntrance::Int, FringeQuadExit::Int, gk::DTPSAD{N, T},
    T1::Array{DTPSAD{N, T},1}, T2::Array{DTPSAD{N, T},1},
    R1::Array{DTPSAD{N, T},2}, R2::Array{DTPSAD{N, T},2}, RApertures::Array{Float64,1}, EApertures::Array{Float64,1},
    KickAngle::Array{DTPSAD{N, T},1}, num_particles::Int, lost_flags::Array{Int64,1}) where {N, T <: Number}

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

    B0 = B[1]
    A0 = A[1]

    if !iszero(KickAngle[1])
        B[1] -= sin(KickAngle[1]) / le
    end
    if !iszero(KickAngle[2])
        A[1] += sin(KickAngle[2]) / le
    end

    hasT1 = !all(iszero, T1)
    hasR1 = !all(iszero, R1)
    hasR2 = !all(iszero, R2)
    hasT2 = !all(iszero, T2)

    @inbounds for c in 1:num_particles
        if isone(lost_flags[c])
            continue
        end
        # Misalignment at entrance
        hasT1 && addvv_row!(r, c, T1)
        hasR1 && multmv_row!(r, c, R1)

        Yrot_row!(r, c, entrance_angle, beti)

        if FringeBendEntrance != 0
            r6 = @view r[c, :]
            bend_fringe!(r6, irho, gk, beti)
        end

        if FringeQuadEntrance != 0
            r6 = @view r[c, :]
            multipole_fringe!(r6, le, A, B, max_order, 1.0, 1, beti)
        end

        bend_edge_row!(r, c, irho, -entrance_angle, beti)

        # Integrator
        if num_int_steps == 0
            exact_bend_row!(r, c, irho, le, beti)
        else
            for m in 1:num_int_steps
                exact_bend_row!(r, c, irho, L1, beti)
                strthinkick_row!(r, c, A, B, K1, max_order)
                exact_bend_row!(r, c, irho, L2, beti)
                strthinkick_row!(r, c, A, B, K2, max_order)
                exact_bend_row!(r, c, irho, L2, beti)
                strthinkick_row!(r, c, A, B, K1, max_order)
                exact_bend_row!(r, c, irho, L1, beti)
            end
        end

        r[c, 5] -= le*beti

        bend_edge_row!(r, c, irho, -exit_angle, beti)
        if FringeQuadExit != 0
            r6 = @view r[c, :]
            multipole_fringe!(r6, le, A, B, max_order, -1.0, 1, beti)
        end
        if FringeBendExit != 0
            r6 = @view r[c, :]
            bend_fringe!(r6, -irho, gk, beti)
        end
        Yrot_row!(r, c, exit_angle, beti)

        # Misalignment at exit
        hasR2 && multmv_row!(r, c, R2)
        hasT2 && addvv_row!(r, c, T2)
        if check_lost_GTPSA_row(r, c)
            lost_flags[c] = 1
        end
    end
    
    if !iszero(KickAngle[1])
        B[1] += sin(KickAngle[1]) / le
    end
    if !iszero(KickAngle[2])
        A[1] -= sin(KickAngle[2]) / le
    end
    return nothing
end

function ExactSectorBend_rad!(r::Matrix{DTPSAD{N, T}}, le::DTPSAD{N, T}, rad_const::DTPSAD{N, T}, beti::Float64, angle::DTPSAD{N, T}, A::Array{DTPSAD{N, T},1}, B::Array{DTPSAD{N, T},1}, 
    max_order::Int, num_int_steps::Int, entrance_angle::DTPSAD{N, T}, exit_angle::DTPSAD{N, T}, FringeBendEntrance::Int, FringeBendExit::Int,
    FringeQuadEntrance::Int, FringeQuadExit::Int, gk::DTPSAD{N, T},
    T1::Array{DTPSAD{N, T},1}, T2::Array{DTPSAD{N, T},1}, 
    R1::Array{DTPSAD{N, T},2}, R2::Array{DTPSAD{N, T},2}, RApertures::Array{Float64,1}, EApertures::Array{Float64,1},
    KickAngle::Array{DTPSAD{N, T},1}, num_particles::Int, lost_flags::Array{Int64,1}) where {N, T <: Number}
    
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

    B0 = B[1]
    A0 = A[1]

    if !iszero(KickAngle[1])
        B[1] -= sin(KickAngle[1]) / le
    end
    if !iszero(KickAngle[2])
        A[1] += sin(KickAngle[2]) / le
    end

    hasT1 = !all(iszero, T1)
    hasR1 = !all(iszero, R1)
    hasR2 = !all(iszero, R2)
    hasT2 = !all(iszero, T2)

    @inbounds for c in 1:num_particles
        if isone(lost_flags[c])
            continue
        end
        # Misalignment at entrance
        hasT1 && addvv_row!(r, c, T1)
        hasR1 && multmv_row!(r, c, R1)

        Yrot_row!(r, c, entrance_angle, beti)

        if FringeBendEntrance != 0
            r6 = @view r[c, :]
            bend_fringe!(r6, irho, gk, beti)
        end

        if FringeQuadEntrance != 0
            r6 = @view r[c, :]
            multipole_fringe!(r6, le, A, B, max_order, 1.0, 1, beti)
        end

        bend_edge_row!(r, c, irho, -entrance_angle, beti)

        # Integrator
        if num_int_steps == 0
            exact_bend_row!(r, c, irho, le, beti)
        else
            for m in 1:num_int_steps
                exact_bend_row!(r, c, irho, L1, beti)
                exactbndthinkick_rad_row!(r, c, A, B, K1, irho, max_order, beti, rad_const)
                exact_bend_row!(r, c, irho, L2, beti)
                exactbndthinkick_rad_row!(r, c, A, B, K2, irho, max_order, beti, rad_const)
                exact_bend_row!(r, c, irho, L2, beti)
                exactbndthinkick_rad_row!(r, c, A, B, K1, irho, max_order, beti, rad_const)
                exact_bend_row!(r, c, irho, L1, beti)
            end
        end

        r[c, 5] -= le*beti

        bend_edge_row!(r, c, irho, -exit_angle, beti)
        if FringeQuadExit != 0
            r6 = @view r[c, :]
            multipole_fringe!(r6, le, A, B, max_order, -1.0, 1, beti)
        end
        if FringeBendExit != 0
            r6 = @view r[c, :]
            bend_fringe!(r6, -irho, gk, beti)
        end
        Yrot_row!(r, c, exit_angle, beti)

        # Misalignment at exit
        hasR2 && multmv_row!(r, c, R2)
        hasT2 && addvv_row!(r, c, T2)
        if check_lost_GTPSA_row(r, c)
            lost_flags[c] = 1
        end
    end
    
    if !iszero(KickAngle[1])
        B[1] += sin(KickAngle[1]) / le
    end
    if !iszero(KickAngle[2])
        A[1] -= sin(KickAngle[2]) / le
    end
    return nothing
end

function multmv!(r::SubArray, A::Matrix{DTPSAD{N, T}}) where {N, T <: Number}
    @inbounds begin
        r1 = A[1, 1] * r[1] + A[1, 2] * r[2] + A[1, 3] * r[3] + A[1, 4] * r[4] + A[1, 5] * r[5] + A[1, 6] * r[6]
        r2 = A[2, 1] * r[1] + A[2, 2] * r[2] + A[2, 3] * r[3] + A[2, 4] * r[4] + A[2, 5] * r[5] + A[2, 6] * r[6]
        r3 = A[3, 1] * r[1] + A[3, 2] * r[2] + A[3, 3] * r[3] + A[3, 4] * r[4] + A[3, 5] * r[5] + A[3, 6] * r[6]
        r4 = A[4, 1] * r[1] + A[4, 2] * r[2] + A[4, 3] * r[3] + A[4, 4] * r[4] + A[4, 5] * r[5] + A[4, 6] * r[6]
        r5 = A[5, 1] * r[1] + A[5, 2] * r[2] + A[5, 3] * r[3] + A[5, 4] * r[4] + A[5, 5] * r[5] + A[5, 6] * r[6]
        r6 = A[6, 1] * r[1] + A[6, 2] * r[2] + A[6, 3] * r[3] + A[6, 4] * r[4] + A[6, 5] * r[5] + A[6, 6] * r[6]
        r[1] = r1
        r[2] = r2
        r[3] = r3
        r[4] = r4
        r[5] = r5
        r[6] = r6
    end
    return nothing
end

@inline function bend_edge_row!(r::Matrix{DTPSAD{N, T}}, c::Int, rhoinv::DTPSAD{N, T}, theta::DTPSAD{N, T}, beti::Float64) where {N, T <: Number}
    if abs(rhoinv) >= 1e-6
        dp1 = 1.0*beti + r[c, 6]
        ctheta = cos(theta)
        stheta = sin(theta)
        pz = pxyz(dp1, r[c, 2], r[c, 4])
        d2 = pxyz(dp1, 0.0, r[c, 4])
        px = r[c, 2] * ctheta + (pz - rhoinv * r[c, 1]) * stheta
        dasin = asin(r[c, 2] / d2) - asin(px / d2)
        num = r[c, 1] * (r[c, 2] * sin(2.0 * theta) + stheta * stheta * (2.0 * pz - rhoinv * r[c, 1]))
        den = pxyz(dp1, px, r[c, 4]) + pxyz(dp1, r[c, 2], r[c, 4]) * ctheta - r[c, 2] * stheta
        x = r[c, 1] * ctheta + num / den
        dy = r[c, 4] * theta / rhoinv + r[c, 4] / rhoinv * dasin
        dct = dp1 / rhoinv * (theta + dasin)

        r[c, 1] = x
        r[c, 2] = px
        r[c, 3] += dy
        r[c, 5] += dct
    end
    return nothing
end

@inline function exact_bend_row!(r::Matrix{DTPSAD{N, T}}, c::Int, irho::DTPSAD{N, T}, L::DTPSAD{N, T}, beti::Float64) where {N, T <: Number}
    dp1 = 1.0*beti + r[c, 6]
    pz = pxyz(dp1, r[c, 2], r[c, 4])

    if abs(irho) < 1e-6
        NormL = L / pz
        r[c, 1] += r[c, 2] * NormL
        r[c, 3] += r[c, 4] * NormL
        r[c, 5] += NormL * dp1
    else
        pzmx = pz - (1.0 + r[c, 1] * irho)
        cs = cos(irho * L)
        sn = sin(irho * L)
        px = r[c, 2] * cs + pzmx * sn
        d2 = pxyz(dp1, 0.0, r[c, 4])
        dasin = L + (asin(r[c, 2] / d2) - asin(px / d2)) / irho
        x = (pxyz(dp1, px, r[c, 4]) - pzmx * cs + r[c, 2] * sn - 1.0) / irho
        dy = r[c, 4] * dasin
        dct = dp1 * dasin

        r[c, 1] = x
        r[c, 2] = px
        r[c, 3] += dy
        r[c, 5] += dct
    end
    return nothing
end

@inline function Yrot_row!(r::Matrix{DTPSAD{N, T}}, c::Int, phi::DTPSAD{N, T}, beti::Float64) where {N, T <: Number}
    if phi != 0.0
        dp1 = 1.0*beti + r[c, 6]
        cphi = cos(phi)
        sphi = sin(phi)
        pz = pxyz(dp1, r[c, 2], r[c, 4])
        p = cphi * pz - sphi * r[c, 2]
        px = sphi * pz + cphi * r[c, 2]
        x = r[c, 1] * pz / p
        dy = r[c, 1] * r[c, 4] * sphi / p
        dct = dp1 * r[c, 1] * sphi / p
        r[c, 1] = x
        r[c, 2] = px
        r[c, 3] += dy
        r[c, 5] += dct
    end
    return nothing
end

function addvv!(r::SubArray, dr::Array{DTPSAD{N, T},1}) where {N, T <: Number}
    @inbounds for i in 1:6
        r[i] += dr[i]
    end
    return nothing
end

@inline function addvv_row!(r::Matrix{DTPSAD{N, T}}, c::Int, dr::AbstractVector{DTPSAD{N, T}}) where {N, T <: Number}
    @inbounds for i in 1:6
        r[c, i] += dr[i]
    end
    return nothing
end

@inline function multmv_row!(r::Matrix{DTPSAD{N, T}}, c::Int, A::AbstractMatrix{DTPSAD{N, T}}) where {N, T <: Number}
    @inbounds begin
        r1 = A[1, 1] * r[c, 1] + A[1, 2] * r[c, 2] + A[1, 3] * r[c, 3] + A[1, 4] * r[c, 4] + A[1, 5] * r[c, 5] + A[1, 6] * r[c, 6]
        r2 = A[2, 1] * r[c, 1] + A[2, 2] * r[c, 2] + A[2, 3] * r[c, 3] + A[2, 4] * r[c, 4] + A[2, 5] * r[c, 5] + A[2, 6] * r[c, 6]
        r3 = A[3, 1] * r[c, 1] + A[3, 2] * r[c, 2] + A[3, 3] * r[c, 3] + A[3, 4] * r[c, 4] + A[3, 5] * r[c, 5] + A[3, 6] * r[c, 6]
        r4 = A[4, 1] * r[c, 1] + A[4, 2] * r[c, 2] + A[4, 3] * r[c, 3] + A[4, 4] * r[c, 4] + A[4, 5] * r[c, 5] + A[4, 6] * r[c, 6]
        r5 = A[5, 1] * r[c, 1] + A[5, 2] * r[c, 2] + A[5, 3] * r[c, 3] + A[5, 4] * r[c, 4] + A[5, 5] * r[c, 5] + A[5, 6] * r[c, 6]
        r6 = A[6, 1] * r[c, 1] + A[6, 2] * r[c, 2] + A[6, 3] * r[c, 3] + A[6, 4] * r[c, 4] + A[6, 5] * r[c, 5] + A[6, 6] * r[c, 6]
        r[c, 1] = r1
        r[c, 2] = r2
        r[c, 3] = r3
        r[c, 4] = r4
        r[c, 5] = r5
        r[c, 6] = r6
    end
    return nothing
end

@inline function _nan_like_dtpsa(x::DTPSAD{N, T}) where {N, T <: Number}
    return DTPSAD{N, T}(T(NaN), x.deriv * T(NaN))
end

function check_lost_GTPSA(r6)
    if isnan(r6[1]) || isinf(r6[1])
        return true
    end
    if max(abs(r6[1]), abs(r6[2]), abs(r6[3]), abs(r6[4])) > CoordLimit || abs(r6[6]) > CoordLimit
        return true
    end
    # sqrt(1.0 + 2.0*r[6]*beti + r[6]^2 - r[2]^2 - r[4]^2) must be real
    if r6[2]^2 + r6[4]^2 > 1.0  + r6[6]^2
        return true
    end
    return false
end

@inline function check_lost_GTPSA_row(r::Matrix{DTPSAD{N, T}}, c::Int) where {N, T <: Number}
    @inbounds begin
        x = r[c, 1]
        px = r[c, 2]
        y = r[c, 3]
        py = r[c, 4]
        dp_p = r[c, 6]
        if isnan(x) || isinf(x)
            return true
        end
        if max(abs(x), abs(px), abs(y), abs(py)) > CoordLimit || abs(dp_p) > CoordLimit
            return true
        end
        if px^2 + py^2 > 1.0 + dp_p^2
            return true
        end
    end
    return false
end

function drift6!(r::SubArray, le::DTPSAD{N, T}, gamma2i::Float64 = 0.0) where {N, T <: Number}
    if use_exact_Hamiltonian == 1
        # check if sqrt(1.0 + 2.0*r[6] + r[6]^2 - r[2]^2 - r[4]^2) is real
        if 1.0 + 2.0*r[6] + r[6]^2 - r[2]^2 - r[4]^2 < 0
            r .= NaN
            return nothing
        end
        NormL = le / sqrt(1.0 + 2.0*r[6] + r[6]^2 - r[2]^2 - r[4]^2)
        r[5] += NormL * (1.0 + r[6]) - le
    else
        NormL = le / (1.0 + r[6])
        r[5] += le * _linearized_drift_phifac(r[2], r[4], r[6], gamma2i)
    end
    r[1] += NormL * r[2]
    r[3] += NormL * r[4]
    return nothing
end 

@inline function drift6_row!(r::Matrix{DTPSAD{N, T}}, c::Int, le::DTPSAD{N, T}, gamma2i::Float64 = 0.0) where {N, T <: Number}
    @inbounds begin
        px = r[c, 2]
        py = r[c, 4]
        dp_p = r[c, 6]
        if use_exact_Hamiltonian == 1
            discr = 1.0 + 2.0 * dp_p + dp_p^2 - px^2 - py^2
            if discr < 0
                nanv = _nan_like_dtpsa(r[c, 1])
                for j in 1:6
                    r[c, j] = nanv
                end
                return nothing
            end
            NormL = le / sqrt(discr)
            r[c, 5] += NormL * (1.0 + dp_p) - le
        else
            NormL = le / (1.0 + dp_p)
            r[c, 5] += le * _linearized_drift_phifac(px, py, dp_p, gamma2i)
        end
        r[c, 1] += NormL * px
        r[c, 3] += NormL * py
    end
    return nothing
end

@inline function strthinkick_row!(r::Matrix{DTPSAD{N, T}}, c::Int, A::Vector{DTPSAD{N, T}}, B::Vector{DTPSAD{N, T}},
    L::DTPSAD{N, T}, max_order::Int) where {N, T <: Number}
    ReSum = B[max_order + 1]
    ImSum = A[max_order + 1]
    ReSumTemp = 0.0
    x = r[c, 1]
    y = r[c, 3]
    for i in max_order-1:-1:0
        ReSumTemp = ReSum * x - ImSum * y + B[i + 1]
        ImSum = ImSum * x + ReSum * y + A[i + 1]
        ReSum = ReSumTemp
    end
    r[c, 2] -= L * ReSum
    r[c, 4] += L * ImSum
    return nothing
end

@inline function strthinkick1_row!(r::Matrix{DTPSAD{N, T}}, c::Int, A::Vector{DTPSAD{N, T}}, B::Vector{DTPSAD{N, T}},
    L::DTPSAD{N, T}, max_order::Int) where {N, T <: Number}
    ReSum = B[max_order + 1]
    ImSum = A[max_order + 1]
    ReSumTemp = 0.0
    x = r[c, 1]
    y = r[c, 3]
    for i in reverse(1:max_order)
        ReSumTemp = ReSum * x - ImSum * y + B[i]
        ImSum = ImSum * x + ReSum * y + A[i]
        ReSum = ReSumTemp
    end
    r[c, 2] -= L * ReSum
    r[c, 4] += L * ImSum
    return nothing
end

@inline function bndthinkick_row!(r::Matrix{DTPSAD{N, T}}, c::Int, A::Array{DTPSAD{N, T}, 1}, B::Array{DTPSAD{N, T}, 1},
    L::DTPSAD{N, T}, irho::DTPSAD{N, T}, max_order::Int, beti::Float64) where {N, T <: Number}
    ReSum = B[max_order + 1]
    ImSum = A[max_order + 1]
    ReSumTemp = 0.0
    x = r[c, 1]
    y = r[c, 3]
    dp_0 = r[c, 6]
    for i in max_order-1:-1:0
        ReSumTemp = ReSum * x - ImSum * y + B[i + 1]
        ImSum = ImSum * x + ReSum * y + A[i + 1]
        ReSum = ReSumTemp
    end
    r[c, 2] -= L * (ReSum - (dp_0 - x * irho) * irho)
    r[c, 4] += L * ImSum
    r[c, 5] += L * irho * x * beti
    return nothing
end

@inline function bndthinkickrad_row!(r::Matrix{DTPSAD{N, T}}, c::Int, A::Array{DTPSAD{N, T}, 1}, B::Array{DTPSAD{N, T}, 1},
    L::DTPSAD{N, T}, irho::DTPSAD{N, T}, E0::DTPSAD{N, T}, max_order::Int, beti::Float64) where {N, T <: Number}
    ReSum = B[max_order + 1]
    ImSum = A[max_order + 1]
    ReSumTemp = 0.0
    CRAD = CGAMMA * E0^3 / (2.0*pi*1e27)

    x = r[c, 1]
    y = r[c, 3]
    for i in max_order-1:-1:0
        ReSumTemp = ReSum * x - ImSum * y + B[i+1]
        ImSum = ImSum * x + ReSum * y + A[i+1]
        ReSum = ReSumTemp
    end

    p_norm = 1.0 / (1.0*beti + r[c, 6])
    xpr = r[c, 2] * p_norm
    ypr = r[c, 4] * p_norm
    B2P = B2perp(ImSum, ReSum + irho, irho, x, xpr, y, ypr)

    dp_0 = r[c, 6]
    r[c, 6] = r[c, 6] - CRAD * (1.0*beti+r[c, 6])^2 * B2P * (1.0 + x*irho + (xpr^2 + ypr^2) / 2.0) * L

    p_norm = 1.0 / (1.0*beti + r[c, 6])
    r[c, 2] = xpr / p_norm
    r[c, 4] = ypr / p_norm

    r[c, 2] -= L * (ReSum - (dp_0 - x*irho)*irho)
    r[c, 4] += L * ImSum
    r[c, 5] += L * irho * x * beti
    return nothing
end

function pass!(ele::DRIFT{DTPSAD{N, T}}, rin::Matrix{DTPSAD{N, T}}, npart::Int, particles::Beam{DTPSAD{N, T}}) where {N, T <: Number}
    gamma2i = use_exact_Hamiltonian == 2 ? 1.0 / particles.gamma.val^2 : 0.0
    hasT1 = !all(iszero, ele.T1)
    hasR1 = !all(iszero, ele.R1)
    hasR2 = !all(iszero, ele.R2)
    hasT2 = !all(iszero, ele.T2)
    @inbounds for c in 1:npart
        lost_flag = particles.lost_flag[c]
        if lost_flag == 1
            continue
        end
        hasT1 && addvv_row!(rin, c, ele.T1)
        hasR1 && multmv_row!(rin, c, ele.R1)
        drift6_row!(rin, c, ele.len, gamma2i)
        hasR2 && multmv_row!(rin, c, ele.R2)
        hasT2 && addvv_row!(rin, c, ele.T2)
        lost_flag = check_lost_GTPSA_row(rin, c) ? 1 : lost_flag
        particles.lost_flag[c] = lost_flag
    end
    return nothing
end

function pass!(ele::QUAD{DTPSAD{N, T}}, rin::Matrix{DTPSAD{N, T}}, npart::Int, particles::Beam{DTPSAD{N, T}}) where {N, T <: Number}
    gamma2i = use_exact_Hamiltonian == 2 ? 1.0 / particles.gamma.val^2 : 0.0
    hasT1 = !all(iszero, ele.T1)
    hasR1 = !all(iszero, ele.R1)
    hasR2 = !all(iszero, ele.R2)
    hasT2 = !all(iszero, ele.T2)
    @inbounds for c in 1:npart
        lost_flag = particles.lost_flag[c]
        if lost_flag == 1
            continue
        end
        p_norm = 1.0 / (1.0 + rin[c, 6])
        
        hasT1 && addvv_row!(rin, c, ele.T1)
        hasR1 && multmv_row!(rin, c, ele.R1)
        
        if iszero(ele.k1)
            drift6_row!(rin, c, ele.len, gamma2i)
        else
            x   = rin[c, 1]
            xpr = rin[c, 2] * p_norm
            y   = rin[c, 3]
            ypr = rin[c, 4] * p_norm
            
            g  = abs(ele.k1) / (1.0 + rin[c, 6])
            t  = sqrt(g)
            lt = ele.len * t
            
            if ele.k1 > 0
                MHD = cos(lt)
                M12 = sin(lt)/t
                M21 = -M12*g
                MVD = cosh(lt)
                M34 = sinh(lt)/t
                M43 = M34*g
            else
                MHD = cosh(lt)
                M12 = sinh(lt)/t
                M21 = M12*g
                MVD = cos(lt)
                M34 = sin(lt)/t
                M43 = -M34*g
            end
            
            rin[c, 1] =  MHD*x + M12*xpr
            rin[c, 2] = (M21*x + MHD*xpr)/p_norm
            rin[c, 3] =  MVD*y + M34*ypr
            rin[c, 4] = (M43*y + MVD*ypr)/p_norm

            rin[c, 5] += g*(x*x*(ele.len-MHD*M12)-y*y*(ele.len-MVD*M34))/4.0
            rin[c, 5] += (xpr*xpr*(ele.len+MHD*M12)+ypr*ypr*(ele.len+MVD*M34))/4.0
            rin[c, 5] += (x*xpr*M12*M21 + y*ypr*M34*M43)/2.0
        end

        hasR2 && multmv_row!(rin, c, ele.R2)
        hasT2 && addvv_row!(rin, c, ele.T2)
        lost_flag = check_lost_GTPSA_row(rin, c) ? 1 : lost_flag
        particles.lost_flag[c] = lost_flag
    end
    return nothing
end

function pass!(ele::CORRECTOR{DTPSAD{N, T}}, rin::Matrix{DTPSAD{N, T}}, npart::Int, particles::Beam{DTPSAD{N, T}}) where {N, T <: Number}
    hasT1 = !all(iszero, ele.T1)
    hasR1 = !all(iszero, ele.R1)
    hasR2 = !all(iszero, ele.R2)
    hasT2 = !all(iszero, ele.T2)
    @inbounds for c in 1:npart
        lost_flag = particles.lost_flag[c]
        if lost_flag == 1
            continue
        end
        hasT1 && addvv_row!(rin, c, ele.T1)
        hasR1 && multmv_row!(rin, c, ele.R1)

        p_norm = 1.0 / (1.0 + rin[c, 6])
        NormL = ele.len * p_norm

        rin[c, 5] += NormL*p_norm*(ele.xkick*ele.xkick/3.0 + ele.ykick*ele.ykick/3.0 +
   		            rin[c, 2]*rin[c, 2] + rin[c, 4]*rin[c, 4] +
   		            rin[c, 2]*ele.xkick + rin[c, 4]*ele.ykick)/2.0
        rin[c, 1] += NormL*(rin[c, 2]+ele.xkick/2.0)
        rin[c, 2] += ele.xkick
        rin[c, 3] += NormL*(rin[c, 4]+ele.ykick/2.0)
        rin[c, 4] += ele.ykick

        hasR2 && multmv_row!(rin, c, ele.R2)
        hasT2 && addvv_row!(rin, c, ele.T2)
        lost_flag = check_lost_GTPSA_row(rin, c) ? 1 : lost_flag
        particles.lost_flag[c] = lost_flag
    end
    return nothing
end

function pass!(ele::MARKER{DTPSAD{N, T}}, rin::Matrix{DTPSAD{N, T}}, npart::Int, particles::Beam{DTPSAD{N, T}}) where {N, T <: Number}
    return nothing
end

function pass!(ele::SOLENOID{DTPSAD{N, T}}, rin::Matrix{DTPSAD{N, T}}, npart::Int, particles::Beam{DTPSAD{N, T}}) where {N, T <: Number}
    gamma2i = use_exact_Hamiltonian == 2 ? 1.0 / particles.gamma.val^2 : 0.0
    hasT1 = !all(iszero, ele.T1)
    hasR1 = !all(iszero, ele.R1)
    hasR2 = !all(iszero, ele.R2)
    hasT2 = !all(iszero, ele.T2)
    @inbounds for c in 1:npart
        lost_flag = particles.lost_flag[c]
        if lost_flag == 1
            continue
        end
        hasT1 && addvv_row!(rin, c, ele.T1)
        hasR1 && multmv_row!(rin, c, ele.R1)

        if ele.ks != 0.0
            p_norm = 1.0 / (1.0 + rin[c, 6])
            x = rin[c, 1]
            xpr = rin[c, 2]*p_norm
            y = rin[c, 3]
            ypr = rin[c, 4]*p_norm

            H = ele.ks * p_norm / 2.0
            S = sin(ele.len * H)
            C = cos(ele.len * H)
            rin[c, 1] = x*C*C + xpr*C*S/H + y*C*S + ypr*S*S/H
            rin[c, 2] = (-x*H*C*S + xpr*C*C - y*H*S*S + ypr*C*S) / p_norm
            rin[c, 3] = -x*C*S - xpr*S*S/H + y*C*C + ypr*C*S/H
            rin[c, 4] = (x*H*S*S - xpr*C*S - y*C*S*H + ypr*C*C) / p_norm
            rin[c, 5] += ele.len*(H*H*(x*x+y*y) + 2.0*H*(xpr*y-ypr*x) +xpr*xpr+ypr*ypr)/2.0
        else
            drift6_row!(rin, c, ele.len, gamma2i)
        end
        hasR2 && multmv_row!(rin, c, ele.R2)
        hasT2 && addvv_row!(rin, c, ele.T2)
        lost_flag = check_lost_GTPSA_row(rin, c) ? 1 : lost_flag
        particles.lost_flag[c] = lost_flag
    end
    return nothing
end

function multipole_fringe!(r6::SubArray, le::DTPSAD{N, T}, polya::Array{DTPSAD{N, T},1}, polyb::Array{DTPSAD{N, T},1}, 
        max_order::Int, edge::Float64, skip_b0::Int, beti::Float64) where {N, T <: Number}
    # Forest 13.29
    FX = 0.0
    FY = 0.0
    FX_X = 0.0
    FX_Y = 0.0
    FY_X = 0.0
    FY_Y = 0.0
  
    RX = 1.0
    IX = 0.0

    for n in 0:max_order
        B = polyb[n + 1]  
        A = polya[n + 1] 
    
        j = n + 1.0
    
        DRX = RX
        DIX = IX
    
        # Complex multiplications
        RX = DRX * r6[1] - DIX * r6[3]
        IX = DRX * r6[3] + DIX * r6[1]
        
        U, V, DU, DV = 0.0, 0.0, 0.0, 0.0
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
    
        nf = 1.0 * (j + 2.0) / j
    
        FX += U * r6[1] + nf * V * r6[3]
        FY += U * r6[3] - nf * V * r6[1]
    
        FX_X += DUX * r6[1] + U + nf * r6[3] * DVX
        FX_Y += DUY * r6[1] + nf * V + nf * r6[3] * DVY
    
        FY_X += DUX * r6[3] - nf * V - nf * r6[1] * DVX
        FY_Y += DUY * r6[3] + U - nf * r6[1] * DVY
    end

    DEL = 1.0 / (1.0*beti + r6[6])
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

function StrB2perp(bx::DTPSAD{N, T}, by::DTPSAD{N, T}, x::DTPSAD{N, T}, xpr::DTPSAD{N, T}, y::DTPSAD{N, T}, ypr::DTPSAD{N, T}) where {N, T <: Number}
    # Calculates sqr(|B x e|) , where e is a unit vector in the direction of velocity
    # v_norm2 = 1.0 / (1.0 + xpr^2 + ypr^2)
    # return (by^2 + bx^2 + (bx*ypr - by*xpr)^2) * v_norm2
    return bx*bx + by*by + (bx*xpr - by*ypr)^2
end

function strthinkickrad!(r::SubArray, A::Vector{DTPSAD{N, T}}, B::Vector{DTPSAD{N, T}}, 
    L::DTPSAD{N, T}, E0::DTPSAD{N, T}, max_order::Int, rad_const::DTPSAD{N, T}) where {N, T <: Number}
    # Modified based on AT function. Ref[Terebilo, Andrei. "Accelerator modeling with MATLAB accelerator toolbox." PACS2001 (2001)].
    ReSum = B[max_order + 1]
    ImSum = A[max_order + 1]
    ReSumTemp = 0.0

    for i in reverse(1:max_order)
        ReSumTemp = ReSum * r[1] - ImSum * r[3] + B[i]
        ImSum = ImSum * r[1] + ReSum * r[3] + A[i]
        ReSum = ReSumTemp
    end

    # angles for momentum
    p_norm = 1.0 / (1.0 + r[6])
    x = r[1]
    xpr = r[2] * p_norm
    y = r[3]
    ypr = r[4] * p_norm
    B2P = StrB2perp(ImSum, ReSum , x , xpr, y ,ypr)
    factor = L / (p_norm)^2 / sqrt(1.0 - xpr^2 - ypr^2)

    dp_0 = r[6]
    r[6] -= rad_const * B2P * factor

    # momentums after losing energy
    p_norm = 1.0 / (1.0 + r[6])

    r[2] = xpr / p_norm
    r[4] = ypr / p_norm

    r[2] -= L * ReSum
    r[4] += L * ImSum
    return nothing
end

function StrMPoleSymplectic4Pass!(r::Matrix{DTPSAD{N, T}}, le::DTPSAD{N, T}, beti::Float64, A::Array{DTPSAD{N, T},1}, B::Array{DTPSAD{N, T},1},
    max_order::Int, num_int_step::Int,
    FringeQuadEntrance::Int, FringeQuadExit::Int, # (no fringe), 1 (lee-whiting) or 2 (lee-whiting+elegant-like)
    T1::Array{DTPSAD{N, T},1}, T2::Array{DTPSAD{N, T},1}, R1::Array{DTPSAD{N, T},2}, R2::Array{DTPSAD{N, T},2},
    RApertures::Array{Float64,1}, EApertures::Array{Float64,1}, KickAngle::Array{DTPSAD{N, T},1},
    num_particles::Int, lost_flags::Array{Int64,1}, gamma2i::Float64 = 0.0) where {N, T <: Number}
    # Modified based on AT function. Ref[Terebilo, Andrei. "Accelerator modeling with MATLAB accelerator toolbox." PACS2001 (2001)].

    DRIFT1  =  0.6756035959798286638
    DRIFT2 = -0.1756035959798286639
    KICK1 = 1.351207191959657328
    KICK2 = -1.702414383919314656

    SL = le/num_int_step
    L1 = SL*DRIFT1
    L2 = SL*DRIFT2
    K1 = SL*KICK1
    K2 = SL*KICK2

    if le > 0
        B[1] -= sin(KickAngle[1])/le
        A[1] += sin(KickAngle[2])/le
    end
    hasT1 = !all(iszero, T1)
    hasR1 = !all(iszero, R1)
    hasR2 = !all(iszero, R2)
    hasT2 = !all(iszero, T2)
    @inbounds for c in 1:num_particles
        if lost_flags[c] == 1
            continue
        end
        # Misalignment at entrance
        hasT1 && addvv_row!(r, c, T1)
        hasR1 && multmv_row!(r, c, R1)

        if FringeQuadEntrance != 0 
            r6 = @view r[c, :]
            multipole_fringe!(r6, le, A, B, max_order, 1.0, 1, beti)
        end

        # Integrator
        for m in 1:num_int_step
            if check_lost_GTPSA_row(r, c)
                lost_flags[c] = 1
                break
            end
            drift6_row!(r, c, L1, gamma2i)
            strthinkick_row!(r, c, A, B, K1, max_order)
            if check_lost_GTPSA_row(r, c)
                lost_flags[c] = 1
                break
            end
            drift6_row!(r, c, L2, gamma2i)
            strthinkick_row!(r, c, A, B, K2, max_order)
            if check_lost_GTPSA_row(r, c)
                lost_flags[c] = 1
                break
            end
            drift6_row!(r, c, L2, gamma2i)
            strthinkick_row!(r, c, A, B, K1, max_order)
            if check_lost_GTPSA_row(r, c)
                lost_flags[c] = 1
                break
            end
            drift6_row!(r, c, L1, gamma2i)
        end

        if FringeQuadExit != 0
            r6 = @view r[c, :]
            multipole_fringe!(r6, le, A, B, max_order, -1.0, 1, beti)
        end

        # Misalignment at exit
        hasR2 && multmv_row!(r, c, R2)
        hasT2 && addvv_row!(r, c, T2)
        if check_lost_GTPSA_row(r, c)
            lost_flags[c] = 1
        end
    end
    if le > 0
        B[1] += sin(KickAngle[1]) / le
        A[1] -= sin(KickAngle[2]) / le
    end
    return nothing
end

function StrMPoleSymplectic4RadPass!(r::Matrix{DTPSAD{N, T}}, le::DTPSAD{N, T}, rad_const::DTPSAD{N, T}, beti::Float64, A::Array{DTPSAD{N, T},1}, B::Array{DTPSAD{N, T},1},
    max_order::Int, num_int_step::Int,
    FringeQuadEntrance::Int, FringeQuadExit::Int, # (no fringe), 1 (lee-whiting) or 2 (lee-whiting+elegant-like)
    T1::Array{DTPSAD{N, T},1}, T2::Array{DTPSAD{N, T},1}, R1::Array{DTPSAD{N, T},2}, R2::Array{DTPSAD{N, T},2},
    RApertures::Array{Float64,1}, EApertures::Array{Float64,1}, KickAngle::Array{DTPSAD{N, T},1}, E0::DTPSAD{N, T},
    num_particles::Int, lost_flags::Array{Int64,1}, gamma2i::Float64 = 0.0) where {N, T <: Number}
    # Modified based on AT function. Ref[Terebilo, Andrei. "Accelerator modeling with MATLAB accelerator toolbox." PACS2001 (2001)].

    DRIFT1  =  0.6756035959798286638
    DRIFT2 = -0.1756035959798286639
    KICK1 = 1.351207191959657328
    KICK2 = -1.702414383919314656

    SL = le/num_int_step
    L1 = SL*DRIFT1
    L2 = SL*DRIFT2
    K1 = SL*KICK1
    K2 = SL*KICK2

    if le > 0
        B[1] -= sin(KickAngle[1])/le
        A[1] += sin(KickAngle[2])/le
    end
    hasT1 = !all(iszero, T1)
    hasR1 = !all(iszero, R1)
    hasR2 = !all(iszero, R2)
    hasT2 = !all(iszero, T2)
    @inbounds for c in 1:num_particles
        if lost_flags[c] == 1
            continue
        end
        # Misalignment at entrance
        hasT1 && addvv_row!(r, c, T1)
        hasR1 && multmv_row!(r, c, R1)

        if FringeQuadEntrance != 0 
            r6 = @view r[c, :]
            multipole_fringe!(r6, le, A, B, max_order, 1.0, 1, beti)
        end

        # Integrator
        for m in 1:num_int_step
            if check_lost_GTPSA_row(r, c)
                lost_flags[c] = 1
                break
            end
            r6 = @view r[c, :]
            drift6_row!(r, c, L1, gamma2i)
            strthinkickrad!(r6, A, B, K1, E0, max_order, rad_const)
            if check_lost_GTPSA_row(r, c)
                lost_flags[c] = 1
                break
            end
            r6 = @view r[c, :]
            drift6_row!(r, c, L2, gamma2i)
            strthinkickrad!(r6, A, B, K2, E0, max_order, rad_const)
            if check_lost_GTPSA_row(r, c)
                lost_flags[c] = 1
                break
            end
            r6 = @view r[c, :]
            drift6_row!(r, c, L2, gamma2i)
            strthinkickrad!(r6, A, B, K1, E0, max_order, rad_const)
            if check_lost_GTPSA_row(r, c)
                lost_flags[c] = 1
                break
            end
            drift6_row!(r, c, L1, gamma2i)
        end

        if FringeQuadExit != 0
            r6 = @view r[c, :]
            multipole_fringe!(r6, le, A, B, max_order, -1.0, 1, beti)
        end

        # Misalignment at exit
        hasR2 && multmv_row!(r, c, R2)
        hasT2 && addvv_row!(r, c, T2)
        if check_lost_GTPSA_row(r, c)
            lost_flags[c] = 1
        end
    end
    if le > 0
        B[1] += sin(KickAngle[1]) / le
        A[1] -= sin(KickAngle[2]) / le
    end

    return nothing
end

function pass!(ele::KQUAD{DTPSAD{N, T}}, r_in::Matrix{DTPSAD{N, T}}, num_particles::Int64, particles::Beam{DTPSAD{N, T}}) where {N, T <: Number}
    lost_flags = particles.lost_flag
    PolynomB = zeros(DTPSAD{N, T}, 4)
    E0 = particles.energy
    rad_const = DTPSAD(0.0)
    gamma2i = use_exact_Hamiltonian == 2 ? 1.0 / particles.gamma.val^2 : 0.0

    PolynomB[1] = ele.k0
    PolynomB[2] = ele.k1
    PolynomB[3] = ele.k2 / 2.0
    PolynomB[4] = ele.k3 / 6.0
    if ele.rad == 0
        StrMPoleSymplectic4Pass!(r_in, ele.len, 1.0, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
            ele.FringeQuadEntrance, ele.FringeQuadExit, 
            ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, num_particles, lost_flags, gamma2i)
    else
        if particles.mass == m_e
            rad_const = RAD_CONST_E * particles.gamma^3
        elseif particles.mass == m_p
            rad_const = RAD_CONST_P * particles.gamma^3
        else
            rad_const = DTPSAD(0.0)
            println("SR is not implemented for this particle mass.")
        end
        StrMPoleSymplectic4RadPass!(r_in, ele.len, rad_const, 1.0, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
            ele.FringeQuadEntrance, ele.FringeQuadExit, 
            ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, E0, num_particles, lost_flags, gamma2i)
    end
    return nothing
end

function pass!(ele::KSEXT{DTPSAD{N, T}}, r_in::Matrix{DTPSAD{N, T}}, num_particles::Int64, particles::Beam{DTPSAD{N, T}}) where {N, T <: Number}
    rad_const = DTPSAD(0.0)
    lost_flags = particles.lost_flag
    PolynomB = zeros(DTPSAD{N, T}, 4)
    E0 = particles.energy
    gamma2i = use_exact_Hamiltonian == 2 ? 1.0 / particles.gamma.val^2 : 0.0

    PolynomB[1] = ele.k0
    PolynomB[2] = ele.k1 
    PolynomB[3] = ele.k2 / 2.0
    PolynomB[4] = ele.k3 / 6.0
    if ele.rad == 0
        StrMPoleSymplectic4Pass!(r_in, ele.len, 1.0, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
            ele.FringeQuadEntrance, ele.FringeQuadExit, ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, num_particles, lost_flags, gamma2i)
    else
        if particles.mass == m_e
            rad_const = RAD_CONST_E * particles.gamma^3
        elseif particles.mass == m_p
            rad_const = RAD_CONST_P * particles.gamma^3
        else
            rad_const = DTPSAD(0.0)
            println("SR is not implemented for this particle mass.")
        end
        StrMPoleSymplectic4RadPass!(r_in, ele.len, rad_const, 1.0, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
            ele.FringeQuadEntrance, ele.FringeQuadExit, ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, E0, num_particles, lost_flags, gamma2i)
    end
    return nothing
end

function pass!(ele::KOCT{DTPSAD{N, T}}, r_in::Matrix{DTPSAD{N, T}}, num_particles::Int64, particles::Beam{DTPSAD{N, T}}) where {N, T <: Number}
    rad_const = DTPSAD(0.0)

    lost_flags = particles.lost_flag
    PolynomB = zeros(DTPSAD{N, T}, 4)
    E0 = particles.energy
    gamma2i = use_exact_Hamiltonian == 2 ? 1.0 / particles.gamma.val^2 : 0.0

    PolynomB[1] = ele.k0
    PolynomB[2] = ele.k1 
    PolynomB[3] = ele.k2 / 2.0
    PolynomB[4] = ele.k3 / 6.0
    if ele.rad == 0
        StrMPoleSymplectic4Pass!(r_in, ele.len, 1.0, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
            ele.FringeQuadEntrance, ele.FringeQuadExit, ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, num_particles, lost_flags, gamma2i)
    else
        if particles.mass == m_e
            rad_const = RAD_CONST_E * particles.gamma^3
        elseif particles.mass == m_p
            rad_const = RAD_CONST_P * particles.gamma^3
        else
            rad_const = DTPSAD(0.0)
            println("SR is not implemented for this particle mass.")
        end
        StrMPoleSymplectic4RadPass!(r_in, ele.len, rad_const, 1.0, ele.PolynomA, PolynomB, ele.MaxOrder, ele.NumIntSteps, 
            ele.FringeQuadEntrance, ele.FringeQuadExit, ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures, ele.KickAngle, E0, num_particles, lost_flags, gamma2i)
    end
    return nothing
end

function BendSymplecticPass!(r::Matrix{DTPSAD{N, T}}, le::DTPSAD{N, T}, beti::Float64, irho::DTPSAD{N, T}, A::Array{DTPSAD{N, T},1}, B::Array{DTPSAD{N, T},1}, 
    max_order::Int, num_int_steps::Int, entrance_angle::DTPSAD{N, T}, exit_angle::DTPSAD{N, T}, FringeBendEntrance::Int, FringeBendExit::Int,
    fint1::DTPSAD{N, T}, fint2::DTPSAD{N, T}, gap::DTPSAD{N, T}, FringeQuadEntrance::Int, FringeQuadExit::Int,
    fringeIntM0::Array{DTPSAD{N, T},1}, fringeIntP0::Array{DTPSAD{N, T},1}, T1::Array{DTPSAD{N, T},1}, T2::Array{DTPSAD{N, T},1}, 
    R1::Array{DTPSAD{N, T},2}, R2::Array{DTPSAD{N, T},2}, RApertures::Array{Float64,1}, EApertures::Array{Float64,1},
    KickAngle::Array{DTPSAD{N, T},1}, num_particles::Int, lost_flags::Array{Int64,1}) where {N, T <: Number}
    # Modified based on AT function. Ref[Terebilo, Andrei. "Accelerator modeling with MATLAB accelerator toolbox." PACS2001 (2001)].

    DRIFT1 = 0.6756035959798286638
    DRIFT2 = -0.1756035959798286639
    KICK1 = 1.351207191959657328
    KICK2 = -1.702414383919314656
    SL = le / num_int_steps
    L1 = SL * DRIFT1
    L2 = SL * DRIFT2
    K1 = SL * KICK1
    K2 = SL * KICK2

    if FringeQuadEntrance==2# && !isnothing(fringeIntM0) && !isnothing(fringeIntP0)
        useLinFrEleEntrance = 1
    else
        useLinFrEleEntrance = 0
    end
    if FringeQuadExit==2 #&& !isnothing(fringeIntM0) && !isnothing(fringeIntP0)
        useLinFrEleExit = 1
    else
        useLinFrEleExit = 0
    end

    B[1] -= sin(KickAngle[1]) / le
    A[1] += sin(KickAngle[2]) / le

    hasT1 = !all(iszero, T1)
    hasR1 = !all(iszero, R1)
    hasR2 = !all(iszero, R2)
    hasT2 = !all(iszero, T2)

    @inbounds for c in 1:num_particles
        if isone(lost_flags[c])
            continue
        end
        # Misalignment at entrance
        hasT1 && addvv_row!(r, c, T1)
        hasR1 && multmv_row!(r, c, R1)

        # Edge focus at entrance
        r6 = @view r[c, :]
        edge_fringe_entrance!(r6, irho, entrance_angle, fint1, gap, FringeBendEntrance)

        # Quadrupole gradient fringe entrance
        if !iszero(FringeQuadEntrance) && !iszero(B[2])
            if isone(useLinFrEleEntrance)
                linearQuadFringeElegantEntrance!(r6, B[2], fringeIntM0, fringeIntP0)
            else
                QuadFringePassP!(r6, B[2])
            end
        end

        # Integrator
        for m in 1:num_int_steps
            drift6_row!(r, c, L1)
            bndthinkick_row!(r, c, A, B, K1, irho, max_order, beti)
            drift6_row!(r, c, L2)
            bndthinkick_row!(r, c, A, B, K2, irho, max_order, beti)
            drift6_row!(r, c, L2)
            bndthinkick_row!(r, c, A, B, K1, irho, max_order, beti)
            drift6_row!(r, c, L1)
        end

        # Quadrupole gradient fringe exit
        if !iszero(FringeQuadExit) && !iszero(B[2])
            r6 = @view r[c, :]
            if isone(useLinFrEleExit)
                linearQuadFringeElegantExit!(r6, B[2], fringeIntM0, fringeIntP0)
            else
                QuadFringePassN!(r6, B[2])
            end
        end

        # Edge focus at exit
        edge_fringe_exit!(r6, irho, exit_angle, fint2, gap, FringeBendExit)

        # Misalignment at exit
        hasR2 && multmv_row!(r, c, R2)
        hasT2 && addvv_row!(r, c, T2)
        if check_lost_GTPSA_row(r, c)
            lost_flags[c] = 1
        end
    end

    B[1] += sin(KickAngle[1]) / le
    A[1] -= sin(KickAngle[2]) / le
    return nothing
end

function BendSymplecticPassRad!(r::Matrix{DTPSAD{N, T}}, le::DTPSAD{N, T}, beti::Float64, irho::DTPSAD{N, T}, A::Array{DTPSAD{N, T},1}, B::Array{DTPSAD{N, T},1}, 
    max_order::Int, num_int_steps::Int, entrance_angle::DTPSAD{N, T}, exit_angle::DTPSAD{N, T}, FringeBendEntrance::Int, FringeBendExit::Int,
    fint1::DTPSAD{N, T}, fint2::DTPSAD{N, T}, gap::DTPSAD{N, T}, FringeQuadEntrance::Int, FringeQuadExit::Int,
    fringeIntM0::Array{DTPSAD{N, T},1}, fringeIntP0::Array{DTPSAD{N, T},1}, T1::Array{DTPSAD{N, T},1}, T2::Array{DTPSAD{N, T},1}, 
    R1::Array{DTPSAD{N, T},2}, R2::Array{DTPSAD{N, T},2}, RApertures::Array{Float64,1}, EApertures::Array{Float64,1},
    KickAngle::Array{DTPSAD{N, T},1}, E0::DTPSAD{N, T}, num_particles::Int, lost_flags::Array{Int64,1}) where {N, T <: Number}
    # AT function. Ref[Terebilo, Andrei. "Accelerator modeling with MATLAB accelerator toolbox." PACS2001 (2001)].

    DRIFT1 = 0.6756035959798286638
    DRIFT2 = -0.1756035959798286639
    KICK1 = 1.351207191959657328
    KICK2 = -1.702414383919314656
    SL = le / num_int_steps
    L1 = SL * DRIFT1
    L2 = SL * DRIFT2
    K1 = SL * KICK1
    K2 = SL * KICK2

    if FringeQuadEntrance==2# && !isnothing(fringeIntM0) && !isnothing(fringeIntP0)
        useLinFrEleEntrance = 1
    else
        useLinFrEleEntrance = 0
    end
    if FringeQuadExit==2# && !isnothing(fringeIntM0) && !isnothing(fringeIntP0)
        useLinFrEleExit = 1
    else
        useLinFrEleExit = 0
    end

    B[1] -= sin(KickAngle[1]) / le
    A[1] += sin(KickAngle[2]) / le

    hasT1 = !all(iszero, T1)
    hasR1 = !all(iszero, R1)
    hasR2 = !all(iszero, R2)
    hasT2 = !all(iszero, T2)

    @inbounds for c in 1:num_particles
        if isone(lost_flags[c])
            continue
        end
        # Misalignment at entrance
        hasT1 && addvv_row!(r, c, T1)
        hasR1 && multmv_row!(r, c, R1)
        # Edge focus at entrance
        r6 = @view r[c, :]
        edge_fringe_entrance!(r6, irho, entrance_angle, fint1, gap, FringeBendEntrance)

        # Quadrupole gradient fringe entrance
        if !iszero(FringeQuadEntrance) && !iszero(B[2])
            if useLinFrEleEntrance == 1
                linearQuadFringeElegantEntrance!(r6, B[2], fringeIntM0, fringeIntP0)
            else
                QuadFringePassP!(r6, B[2])
            end
        end

        # Integrator
        for m in 1:num_int_steps
            drift6_row!(r, c, L1)
            bndthinkickrad_row!(r, c, A, B, K1, irho, E0, max_order, beti)
            drift6_row!(r, c, L2)
            bndthinkickrad_row!(r, c, A, B, K2, irho, E0, max_order, beti)
            drift6_row!(r, c, L2)
            bndthinkickrad_row!(r, c, A, B, K1, irho, E0, max_order, beti)
            drift6_row!(r, c, L1)
        end

        # Quadrupole gradient fringe exit
        if !iszero(FringeQuadExit) && !iszero(B[2])
            r6 = @view r[c, :]
            if useLinFrEleExit == 1
                linearQuadFringeElegantExit!(r6, B[2], fringeIntM0, fringeIntP0)
            else
                QuadFringePassN!(r6, B[2])
            end
        end

        # Edge focus at exit
        edge_fringe_exit!(r6, irho, exit_angle, fint2, gap, FringeBendExit)

        # Misalignment at exit
        hasR2 && multmv_row!(r, c, R2)
        hasT2 && addvv_row!(r, c, T2)
        if check_lost_GTPSA_row(r, c)
            lost_flags[c] = 1
        end
    end

    
    B[1] += sin(KickAngle[1]) / le
    A[1] -= sin(KickAngle[2]) / le
    return nothing
end

function pass!(ele::SBEND{DTPSAD{N, T}}, r_in::Matrix{DTPSAD{N, T}}, num_particles::Int64, particles::Beam{DTPSAD{N, T}}) where {N, T <: Number}
    lost_flags = particles.lost_flag
    irho = ele.angle / ele.len
    E0 = particles.energy
    if ele.rad == 0
        BendSymplecticPass!(r_in, ele.len, 1.0, irho, ele.PolynomA, ele.PolynomB, ele.MaxOrder, ele.NumIntSteps,
            ele.e1, ele.e2,
            ele.FringeBendEntrance, ele.FringeBendExit,
            ele.fint1, ele.fint2, ele.gap,
            ele.FringeQuadEntrance, ele.FringeQuadExit,
            ele.FringeIntM0, ele.FringeIntP0,
            ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures,
            ele.KickAngle, num_particles, lost_flags)
    else
        BendSymplecticPassRad!(r_in, ele.len, 1.0, irho, ele.PolynomA, ele.PolynomB, ele.MaxOrder, ele.NumIntSteps,
            ele.e1, ele.e2,
            ele.FringeBendEntrance, ele.FringeBendExit,
            ele.fint1, ele.fint2, ele.gap,
            ele.FringeQuadEntrance, ele.FringeQuadExit,
            ele.FringeIntM0, ele.FringeIntP0,
            ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures,
            ele.KickAngle, E0, num_particles, lost_flags)
    end
    return nothing
end

function pass!(ele::ESBEND{DTPSAD{N, T}}, r_in::Matrix{DTPSAD{N, T}}, num_particles::Int64, particles::Beam{DTPSAD{N, T}}) where {N, T <: Number}
    lost_flags = particles.lost_flag
    E0 = particles.energy
    rad_const = DTPSAD(0.0)
    if ele.rad == 0
        ExactSectorBend!(r_in, ele.len, 1.0, ele.angle, ele.PolynomA, ele.PolynomB, ele.MaxOrder, ele.NumIntSteps,
            ele.e1, ele.e2,
            ele.FringeBendEntrance, ele.FringeBendExit,
            ele.FringeQuadEntrance, ele.FringeQuadExit,
            ele.gK,
            ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures,
            ele.KickAngle, num_particles, lost_flags)
    else
        if particles.mass == m_e
            rad_const = RAD_CONST_E * particles.gamma^3
        elseif particles.mass == m_p
            rad_const = RAD_CONST_P * particles.gamma^3
        else
            rad_const = DTPSAD(0.0)
            println("SR is not implemented for this particle mass.")
        end
        ExactSectorBend_rad!(r_in, ele.len, rad_const, 1.0, ele.angle, ele.PolynomA, ele.PolynomB, ele.MaxOrder, ele.NumIntSteps,
            ele.e1, ele.e2,
            ele.FringeBendEntrance, ele.FringeBendExit,
            ele.FringeQuadEntrance, ele.FringeQuadExit,
            ele.gK,
            ele.T1, ele.T2, ele.R1, ele.R2, ele.RApertures, ele.EApertures,
            ele.KickAngle, num_particles, lost_flags)
    end
    return nothing
end

function pass!(ele::RFCA{DTPSAD{N, T}}, rin::Matrix{DTPSAD{N, T}}, npart::Int64, particles::Beam{DTPSAD{N, T}}) where {N, T <: Number}
    if ele.energy == 0
        println("Energy is not defined for RFCA ", ele.name)
    end
    nv = ele.volt / ele.energy
    gamma2i = use_exact_Hamiltonian == 2 ? 1.0 / particles.gamma.val^2 : 0.0
    if ele.len == 0
        @inbounds for c in 1:npart
            if particles.lost_flag[c] == 1
                continue
            end
            rin[c, 6] += -nv * sin(2 * pi * ele.freq * ((rin[c, 5] - ele.lag) / speed_of_light - 
                (ele.h / ele.freq - particles.T0) * 0.0) - ele.philag) / particles.beta^2
            if check_lost_GTPSA_row(rin, c)
                particles.lost_flag[c] = 1
            end
        end
    else
        halflength = ele.len / 2.0
        @inbounds for c in 1:npart
            if particles.lost_flag[c] == 1
                continue
            end
            drift6_row!(rin, c, halflength, gamma2i)
            rin[c, 6] += -nv * sin(2 * pi * ele.freq * ((rin[c, 5] - ele.lag) / speed_of_light - 
                (ele.h / ele.freq - particles.T0) * 0.0) - ele.philag) / particles.beta^2
            drift6_row!(rin, c, halflength, gamma2i)
            if check_lost_GTPSA_row(rin, c)
                particles.lost_flag[c] = 1
            end
        end
    end
    return nothing
end

function pass!(cavity::CRABCAVITY{DTPSAD{N, T}}, rin::Matrix{DTPSAD{N, T}}, npart::Int64, particles::Beam{DTPSAD{N, T}}) where {N, T <: Number}
    beta = particles.beta
    gamma2i = use_exact_Hamiltonian == 2 ? 1.0 / particles.gamma.val^2 : 0.0
    @inbounds for c in 1:npart
        if isone(particles.lost_flag[c])
            continue
        end
        ang = cavity.k * rin[c, 5] + cavity.phi
        if cavity.len == 0.0
            rin[c, 2] += (cavity.volt/particles.energy) * sin(ang/beta)
            rin[c, 6] += (-cavity.k * cavity.volt/particles.energy/beta) * rin[c, 1] * cos(ang/beta)
        else
            drift6_row!(rin, c, cavity.len / 2.0, gamma2i)
            rin[c, 2] += (cavity.volt/particles.energy) * sin(ang/beta)
            rin[c, 6] += (-cavity.k * cavity.volt/particles.energy/beta) * rin[c, 1] * cos(ang/beta)
            drift6_row!(rin, c, cavity.len / 2.0, gamma2i)
        end
    end
    return nothing
end

function pass!(cavity::CRABCAVITY_K2{DTPSAD{N, T}}, rin::Matrix{DTPSAD{N, T}}, npart::Int64, particles::Beam{DTPSAD{N, T}}) where {N, T <: Number}
    beta = particles.beta
    gamma2i = use_exact_Hamiltonian == 2 ? 1.0 / particles.gamma.val^2 : 0.0
    @inbounds for c in 1:npart
        if isone(particles.lost_flag[c])
            continue
        end
        ang = cavity.k * rin[c, 5] + cavity.phi
        if cavity.len == 0.0
            rin[c, 2] += (cavity.volt/particles.energy) * sin(ang/beta)
            rin[c, 6] += (-cavity.k * cavity.volt/particles.energy/beta) * rin[c, 1] * cos(ang/beta)
            rin[c, 2] -= cavity.k2 * (rin[c, 1]^2 - rin[c, 3]^2) * sin(ang)
            rin[c, 4] += 2.0 * cavity.k2 * rin[c, 1] * rin[c, 3] * sin(ang)
            rin[c, 6] -= (cavity.k2 * cavity.k / 3.0) * (rin[c, 1]^3 - 3.0 * rin[c, 1] * rin[c, 3]^2) * cos(ang)
        else
            drift6_row!(rin, c, cavity.len / 2.0, gamma2i)
            rin[c, 2] += (cavity.volt/particles.energy) * sin(ang/beta)
            rin[c, 6] += (-cavity.k * cavity.volt/particles.energy/beta) * rin[c, 1] * cos(ang/beta)
            rin[c, 2] -= cavity.k2 * (rin[c, 1]^2 - rin[c, 3]^2) * sin(ang)
            rin[c, 4] += 2.0 * cavity.k2 * rin[c, 1] * rin[c, 3] * sin(ang)
            rin[c, 6] -= (cavity.k2 * cavity.k / 3.0) * (rin[c, 1]^3 - 3.0 * rin[c, 1] * rin[c, 3]^2) * cos(ang)
            drift6_row!(rin, c, cavity.len / 2.0, gamma2i)
        end
    end
    return nothing
end

function pass!(ele::thinMULTIPOLE{DTPSAD{N, T}}, rin::Matrix{DTPSAD{N, T}}, npart::Int64, particles::Beam{DTPSAD{N, T}}) where {N, T <: Number}
    # Modified based on AT function. Ref[Terebilo, Andrei. "Accelerator modeling with MATLAB accelerator toolbox." PACS2001 (2001)].

    # no bending
    bax = 0.0 
    bay = 0.0

    A = [ele.PolynomA[i] for i in 1:4]
    B = [ele.PolynomB[i] for i in 1:4]
    B[1] -= ele.KickAngle[1]
    A[1] += ele.KickAngle[2]

    # for c in 1:npart
    @inbounds for c in 1:npart
        if particles.lost_flag[c] == 1
            continue
        end
        # Misalignment at entrance
        if !all(iszero, ele.T1)
            addvv_row!(rin, c, ele.T1)
        end
        if !all(iszero, ele.R1)
            multmv_row!(rin, c, ele.R1)
        end

        strthinkick1_row!(rin, c, A, B, DTPSAD(1.0), ele.MaxOrder)
        rin[c, 2] += bax * rin[c, 6]
        rin[c, 4] -= bay * rin[c, 6]
        rin[c, 6] -= bax * rin[c, 1] - bay * rin[c, 3]  # Path lenghtening

        # Misalignment at exit
        if !all(iszero, ele.R2)
            multmv_row!(rin, c, ele.R2)
        end
        if !all(iszero, ele.T2)
            addvv_row!(rin, c, ele.T2)
        end
        if check_lost_GTPSA_row(rin, c)
            particles.lost_flag[c] = 1
        end
    end

    B[1] += ele.KickAngle[1]
    A[1] -= ele.KickAngle[2]
    return nothing
end

function pass!(elem::TRANSLATION{DTPSAD{N, T}}, r_in::Matrix{DTPSAD{N, T}}, num_particles::Int64, particles::Beam{DTPSAD{N, T}}) where {N, T <: Number}
    @inbounds for c in 1:num_particles
        if isone(particles.lost_flag[c])
            continue
        end
        if use_exact_beti == 1
            pz = sqrt(1.0 + 2.0 * r_in[c, 6] / particles.beta + r_in[c, 6]^2 - r_in[c, 2]^2 - r_in[c, 4]^2)
            r_in[c, 1] -= elem.dx + elem.ds * r_in[c, 2] / pz
            r_in[c, 3] -= elem.dy + elem.ds * r_in[c, 4] / pz
            r_in[c, 5] += elem.ds * (1.0/particles.beta + r_in[c, 6]) / pz
        else
            pz = sqrt(1.0 + 2.0 * r_in[c, 6] + r_in[c, 6]^2 - r_in[c, 2]^2 - r_in[c, 4]^2)
            r_in[c, 1] -= elem.dx + elem.ds * r_in[c, 2] / pz
            r_in[c, 3] -= elem.dy + elem.ds * r_in[c, 4] / pz
            r_in[c, 5] += elem.ds * (1.0 + r_in[c, 6]) / pz
        end
        if check_lost_GTPSA_row(r_in, c)
            particles.lost_flag[c] = 1
        end
    end
    return nothing
end

function pass!(elem::YROTATION{DTPSAD{N, T}}, r_in::Matrix{DTPSAD{N, T}}, num_particles::Int64, particles::Beam{DTPSAD{N, T}}) where {N, T <: Number}
    angle = -elem.angle
    if angle == 0.0
        return nothing
    end
    ca = cos(angle)
    sa = sin(angle)
    ta = tan(angle)
    if use_exact_beti == 1
        beta = particles.beta
    else
        beta = 1.0
    end
    @inbounds for c in 1:num_particles
        if isone(particles.lost_flag[c])
            continue
        end
        x, px, y, py, t, pt = r_in[c, 1], r_in[c, 2], r_in[c, 3], r_in[c, 4], r_in[c, 5], r_in[c, 6]
        pz = sqrt(1.0 + 2.0 * pt / beta + pt^2 - px^2 - py^2)
        ptt = 1.0 - ta*px/pz
        r_in[c, 1] = x/(ca*ptt)
        r_in[c, 2] = ca*px + sa*pz
        r_in[c, 3] = y + ta*x*py/(pz*ptt)
        r_in[c, 5] = t + ta*x*(1.0 / beta+pt)/(pz*ptt)

        if check_lost_GTPSA_row(r_in, c)
            particles.lost_flag[c] = 1
        end
    end
    return nothing
end

function pass!(elem::LongitudinalRLCWake{DTPSAD{N, T}}, r_in::Matrix{DTPSAD{N, T}}, num_particles::Int64, particles::Beam{DTPSAD{N, T}}) where {N, T <: Number}
    histogram1DinZ!(particles, particles.znbin, particles.inzindex, particles.zhist, particles.zhist_edges)
    eN_b2E=particles.np*1.6021766208e-19*particles.charge^2/particles.energy/particles.beta/particles.beta/particles.atomnum
    zhist_center = zeros(DTPSAD{NVAR(), Float64}, particles.znbin)
    zhist_center .= ((particles.zhist_edges[1:end-1]) .+ (particles.zhist_edges[2:end]))./2.0
    wakefield = zeros(DTPSAD{NVAR(), Float64}, particles.znbin)
    for i in 1:particles.znbin
        t = (zhist_center[i] -  0.0*zhist_center[end]) / 2.99792458e8
        wakefield[i] = wakefieldfunc_RLCWake(elem, t)
    end
    wakepotential = zeros(DTPSAD{NVAR(), Float64}, particles.znbin)

    halfzn = particles.znbin ÷ 2
    for i=1:particles.znbin
        for j=-halfzn:halfzn
            if i-j>0 && i-j<=particles.znbin
                wakepotential[i]+=wakefield[j+halfzn+1]*particles.zhist[i-j]/particles.nmacro
            end
        end
    end

    wakeatedge = zeros(DTPSAD{NVAR(), Float64}, particles.znbin+1)
    wakeatedge[2:end-1] .= ((wakepotential[1:end-1]) .+ (wakepotential[2:end])) ./ 2.0
    wakeatedge[1] = 2*wakeatedge[2]-wakeatedge[3]
    wakeatedge[end] = 2*wakeatedge[end-1]-wakeatedge[end-2]

    zsep = particles.zhist_edges[2]-particles.zhist_edges[1]
    @inbounds for c in 1:num_particles
        if isone(particles.lost_flag[c])
            continue
        end
        zloc = r_in[c, 5]
        zindex=particles.inzindex[c]
        wake1=wakeatedge[zindex]
        wake2=wakeatedge[zindex+1]
        wakezloc=wake1+(wake2-wake1)*(zloc-particles.zhist_edges[zindex])/zsep
        r_in[c, 6] -= wakezloc*eN_b2E

        if check_lost_GTPSA_row(r_in, c)
            particles.lost_flag[c] = 1
        end
    end
    return nothing
end

function Bassetti_Erskine_xgty!(res::Vector{DTPSAD{N, T}}, x::DTPSAD{N, T}, y::DTPSAD{N, T}, sigmax::DTPSAD{N, T}, sigmay::DTPSAD{N, T}) where {N, T <: Number}
    # Only positive y is valid for this function
    # for y<0, Ex = Ex, Ey = -Ey
    if y < 0.0
        Bassetti_Erskine_xgty!(res, x, -y, sigmax, sigmay)
        res[2] = -res[2]
        return nothing
    end
    termexp=exp(-x*x/2.0/sigmax/sigmax-y*y/2.0/sigmay/sigmay)
	sqrtδsigma2=sqrt(Complex(2*(sigmax*sigmax-sigmay*sigmay)))
	term1=erfcx(-1.0im*(x+1.0im*y)/sqrtδsigma2)
	term2=erfcx(-1.0im*(x*sigmay/sigmax+1.0im*y*sigmax/sigmay)/sqrtδsigma2)

	complex_e=-1.0im*2*sqrt(pi)/sqrtδsigma2*(term1-termexp*term2)
	res[1]=real(complex_e)
    res[2]=-imag(complex_e)
    res[3]=termexp/2.0/π/sigmax/sigmay
    return nothing
end
function Bassetti_Erskine_ygtx!(res::Vector{DTPSAD{N, T}}, x::DTPSAD{N, T}, y::DTPSAD{N, T}, sigmax::DTPSAD{N, T}, sigmay::DTPSAD{N, T}) where {N, T <: Number}
    # Only negative x is valid for this function
    # for x>0, Ex = -Ex, Ey = Ey
    if x > 0.0
        Bassetti_Erskine_ygtx!(res, -x, y, sigmax, sigmay)
        res[1] = -res[1]
        return nothing
    end
    termexp=exp(-x*x/2/sigmax/sigmax-y*y/2/sigmay/sigmay)
	sqrtδsigma2=sqrt(Complex(2*(sigmax*sigmax-sigmay*sigmay)))
	term1=erfcx(-1.0im*(x+1.0im*y)/sqrtδsigma2)
	term2=erfcx(-1.0im*(x*sigmay/sigmax+1.0im*y*sigmax/sigmay)/sqrtδsigma2)

	complex_e=-1.0im*2*sqrt(pi)/sqrtδsigma2*(term1-termexp*term2)
    res[1]=real(complex_e)
    res[2]=-imag(complex_e)
    res[3]=termexp/2.0/π/sigmax/sigmay
	return nothing
end
function Bassetti_Erskine!(res::Vector{DTPSAD{N, T}}, x::DTPSAD{N, T}, y::DTPSAD{N, T}, sigmax::DTPSAD{N, T}, sigmay::DTPSAD{N, T}) where {N, T <: Number}
    if sigmax > sigmay
        Bassetti_Erskine_xgty!(res, x, y, sigmax, sigmay)
        return nothing
    else
        Bassetti_Erskine_ygtx!(res, x, y, sigmax, sigmay) 
        return nothing
    end
end
function track_sbb!(rin::Matrix{DTPSAD{N, T}}, num_macro::Int64, temp1::Vector{DTPSAD{N, T}}, temp2::Vector{DTPSAD{N, T}}, 
    temp3::Vector{DTPSAD{N, T}}, temp4::Vector{DTPSAD{N, T}}, temp5::Vector{DTPSAD{N, T}}, sgb::StrongGaussianBeam, factor::DTPSAD{N, T}) where {N, T <: Number}
    #factor=wb.particle.classrad0/wb.gamma*wb.particle.charge*sgb.particle.charge
    
    lumi=DTPSAD{N, T}(0.0)
    fieldvec = zeros(DTPSAD{N, T}, 3)

    @inbounds for i in 1:sgb.nzslice
        slicelumi=DTPSAD{N, T}(0.0)
        @inbounds for j in 1:num_macro
            r6 = @view rin[j, :]
            # temp1: collision zlocation, temp2: beamsize x, temp3: beamsize y, temp4: beta x, temp5: beta y
            temp1[j] = (r6[5] .+ sgb.zslice_center[i])./2.0
            temp4[j] = sgb.optics.optics_x.beta .+ sgb.optics.optics_x.gamma .* temp1[j] .* temp1[j] .- 2.0 .* sgb.optics.optics_x.alpha .* temp1[j]
            temp2[j] = sgb.beamsize[1] .* sqrt.(temp4[j] ./ sgb.optics.optics_x.beta)
            temp5[j] = sgb.optics.optics_y.beta .+ sgb.optics.optics_y.gamma .* temp1[j] .* temp1[j] .- 2.0 .* sgb.optics.optics_y.alpha .* temp1[j]
            temp3[j] = sgb.beamsize[2] .* sqrt.(temp5[j] ./ sgb.optics.optics_y.beta)
        
            # temp4 and temp5 are free to change now.
            r6[1] += (r6[2] .* temp1[j])
            r6[3] += (r6[4] .* temp1[j])
            Bassetti_Erskine!(fieldvec, r6[1], r6[3], temp2[j], temp3[j])
            r6[2] += (sgb.zslice_npar[i]*factor) * fieldvec[1]
            r6[4] += (sgb.zslice_npar[i]*factor) * fieldvec[2]
            slicelumi += fieldvec[3]

            r6[1] -= (r6[2] .* temp1[j])
            r6[3] -= (r6[4] .* temp1[j])
        end
       
        lumi += slicelumi * sgb.zslice_npar[i] #  Will do it outside* wb.num_particle / wb.num_macro
    end

    return lumi

end
function pass!(sgb::StrongGaussianBeam{DTPSAD{N, T}}, r_in::Matrix{DTPSAD{N, T}}, num_macro::Int, wb::Beam{DTPSAD{N, T}}) where {N, T <: Number}
    factor=wb.classrad0/wb.gamma*wb.charge*sgb.charge
    lumi=track_sbb!(r_in, num_macro, wb.temp1, wb.temp2, wb.temp3, wb.temp4, wb.temp5, sgb, factor)
    lumi *= wb.np / wb.nmacro
    return nothing
end

function linepass!(line::Vector{<:AbstractElement{DTPSAD{N, T}}}, particles::Beam{DTPSAD{N, T}}) where {N, T <: Number}
    np = particles.nmacro
    for ele in line
        pass!(ele, particles.r, np, particles)
    end
    return nothing
end

function linepass!(line::Vector{<:AbstractElement{DTPSAD{N, T}}}, particles::Beam{DTPSAD{N, T}}, refpts::Vector) where {N, T <: Number}
    # Note!!! A lost particle's coordinate will not be marked as NaN or Inf like other softwares 
    # Check if the particle is lost by checking the lost_flag
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

function ringpass!(line::Vector{<:AbstractElement{DTPSAD{N, T}}}, particles::Beam{DTPSAD{N, T}}, nturns::Int) where {N, T <: Number}
    for turn in 1:nturns
        linepass!(line, particles)
    end
    return nothing
end
# Conversion functions for AbstractElement{Float64} ↔ AbstractElement{DTPSAD{N, T}}

function _strip(param::DTPSAD{N, T}) where {N, T}
    return param.val
end

function _strip(param::Vector{DTPSAD{N, T}}) where {N, T}
    return [_strip(p) for p in param]
end

function _strip(param::Matrix{DTPSAD{N, T}}) where {N, T}
    Mat = zeros(Float64, size(param))
    for i in eachindex(param)
        Mat[i] = _strip(param[i])
    end
    return Mat
end

function _strip(param::Union{String, Float64, Int64, Int32, Vector{Float64}, Vector{Int64}, Vector{Int32}, Matrix{Float64}, Matrix{Int64}, Matrix{Int32}})
    return param
end

function _strip(param)
    # Fallback for any other type
    return param
end

# Convert AbstractElement{DTPSAD{N, T}} to AbstractElement{Float64}
function to_Number(te::AbstractElement{DTPSAD{N, T}}) where {N, T}
    # Get the base element type by removing the DTPSAD parameterization
    base_type = typeof(te).name.wrapper
    
    # Extract field values and convert DTPSAD to Float64
    vals = map(f -> _strip(getfield(te, f)), fieldnames(typeof(te)))
    
    # Create the Float64 version
    return base_type(vals...)
end

# Convert vector of TPSA elements to Number elements
function TPSAD2Number(line::Vector{<:AbstractElement{DTPSAD{N, T}}}) where {N, T}
    return [to_Number(ele) for ele in line]
end

# Reverse conversion functions
function _strip_inverse(param::Float64, fieldname::Symbol, ::Type{DTPSAD{N, T}}) where {N, T}
    if fieldname === :green_dx || fieldname === :green_dy
        return param
    end
    return DTPSAD(param)
end

function _strip_inverse(param::Vector{Float64}, fieldname::Symbol, ::Type{DTPSAD{N, T}}) where {N, T}
    if fieldname == :RApertures || fieldname == :EApertures ||
       fieldname == :z_grid || fieldname == :z_deriv_grid
        return param  # Keep apertures as Float64 vectors
    else
        return [DTPSAD(p) for p in param]
    end
end

function _strip_inverse(param::Matrix{Float64}, fieldname::Symbol, ::Type{DTPSAD{N, T}}) where {N, T}
    if fieldname == :rho_grid || fieldname == :phi_grid
        return param
    end
    Mat = zeros(DTPSAD{N, T}, size(param)...)
    for i in eachindex(param)
        Mat[i] = DTPSAD(param[i])
    end
    return Mat
end

function _strip_inverse(param::Union{String, Int64, Int32, Vector{Int64}, Vector{Int32}, Matrix{Int64}, Matrix{Int32}}, fieldname::Symbol, ::Type{DTPSAD{N, T}}) where {N, T}
    return param  # Keep non-Float64 types as-is
end

function _strip_inverse(param, fieldname::Symbol, ::Type{DTPSAD{N, T}}) where {N, T}
    return param  # Fallback for any other type
end

# Convert AbstractElement{Float64} to AbstractElement{DTPSAD{N, T}}
function to_TPSAD(ele::AbstractElement{Float64}, ::Type{DTPSAD{N, T}}) where {N, T}
    # Get the base element type
    base_type = typeof(ele).name.wrapper
    
    # Extract field values and convert Float64 to DTPSAD{N, T}
    vals = map(f -> _strip_inverse(getfield(ele, f), f, DTPSAD{N, T}), fieldnames(typeof(ele)))
    
    # Create the DTPSAD version
    return base_type(vals...)
end

# Convenience function using current NVAR()
function to_TPSAD(ele::AbstractElement{Float64})
    return to_TPSAD(ele, DTPSAD{NVAR(), Float64})
end

# Convert vector of Number elements to TPSA elements
function Number2TPSAD(line::Vector{<:AbstractElement{Float64}}, ::Type{DTPSAD{N, T}}) where {N, T}
    return [to_TPSAD(ele, DTPSAD{N, T}) for ele in line]
end

# Convenience function using current NVAR()
function Number2TPSAD(line::Vector{<:AbstractElement{Float64}})
    return Number2TPSAD(line, DTPSAD{NVAR(), Float64})
end

# Helper functions for specific NVAR dimensions
function to_TPSAD_N(ele::AbstractElement{Float64}, N::Int)
    return to_TPSAD(ele, DTPSAD{N, Float64})
end

function Number2TPSAD_N(line::Vector{<:AbstractElement{Float64}}, N::Int)
    return Number2TPSAD(line, DTPSAD{N, Float64})
end

# Generic conversion between different DTPSAD dimensions
function convert_TPSAD_dimension(te::AbstractElement{DTPSAD{N1, T}}, ::Type{DTPSAD{N2, T}}) where {N1, N2, T}
    # First convert to Float64, then to target DTPSAD type
    number_ele = to_Number(te)
    return to_TPSAD(number_ele, DTPSAD{N2, T})
end

function convert_TPSAD_dimension(line::Vector{<:AbstractElement{DTPSAD{N1, T}}}, ::Type{DTPSAD{N2, T}}) where {N1, N2, T}
    return [convert_TPSAD_dimension(ele, DTPSAD{N2, T}) for ele in line]
end

# Convenience functions for dimension conversion
function change_TPSA_dimension(te::AbstractElement{DTPSAD{N1, T}}, N2::Int) where {N1, T}
    return convert_TPSAD_dimension(te, DTPSAD{N2, T})
end

function change_TPSA_dimension(line::Vector{<:AbstractElement{DTPSAD{N1, T}}}, N2::Int) where {N1, T}
    return convert_TPSAD_dimension(line, DTPSAD{N2, T})
end
