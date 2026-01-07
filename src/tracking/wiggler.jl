# WIGGLER tracking Modified based on AT function 'GWigSymplecticPass' and 'GWigSymplecticRadPass'. 
# Ref[Terebilo, Andrei. "Accelerator modeling with MATLAB accelerator toolbox." PACS2001 (2001)].

# When synchrtron radtiaon is enabled, this code may result in small differences compared to 
# 'GWigSymplecticRadPass' of AT because AT sometimes skips the last step of SR kicks due to 
# floating-point rounding, which is unphysical.

function _sinc_taylor(x::Float64)
    x2 = x*x
    return 1.0 - x2/6.0*(1.0 - x2/20.0*(1.0 - x2/42.0*(1.0 - x2/72.0)))
end

function _wig_ax_axpy!(ele::WIGGLER, r::AbstractVector{Float64}, Zw::Float64, Aw::Float64, Po::Float64)
    x = r[1]
    y = r[3]
    kw = 2π / ele.lw

    ax = 0.0
    axpy = 0.0

    By = ele.By
    NHharm = ele.NHharm
    
    for i in 1:NHharm
        if 6*i > length(By)
            break
        end
        base = (i-1)*6

        HCw = float(By[base+2]) * Aw / Po
        kx = float(By[base+3]) * kw
        ky = float(By[base+4]) * kw
        kz = float(By[base+5]) * kw
        tz = float(By[base+6])

        cx = cos(kx*x)
        chy = cosh(ky*y)
        sz = sin(kz*Zw + tz)
        ax += HCw * (kw/kz) * cx * chy * sz

        shy = sinh(ky*y)
        sxkx = abs(kx/kw) > 1e-6 ? sin(kx*x)/kx : x * _sinc_taylor(kx*x)
        axpy += HCw * (kw/kz) * ky * sxkx * shy * sz
    end

    Bx = ele.Bx
    NVharm = ele.NVharm
    
    for i in 1:NVharm
        if 6*i > length(Bx)
            break
        end
        base = (i-1)*6
        VCw = float(Bx[base+2]) * Aw / Po
        kx = float(Bx[base+3]) * kw
        ky = float(Bx[base+4]) * kw
        kz = float(Bx[base+5]) * kw
        tz = float(Bx[base+6])

        shx = sinh(kx*x)
        sy = sin(ky*y)
        sz = sin(kz*Zw + tz)
        ax += VCw * (kw/kz) * (ky/kx) * shx * sy * sz

        chx = cosh(kx*x)
        cy = cos(ky*y)
        axpy += VCw * (kw/kz) * (ky/kx)^2 * chx * cy * sz
    end

    return ax, axpy
end

function _wig_ay_aypx!(ele::WIGGLER, r::AbstractVector{Float64}, Zw::Float64, Aw::Float64, Po::Float64)
    x = r[1]
    y = r[3]
    kw = 2π / ele.lw

    ay = 0.0
    aypx = 0.0

    By = ele.By
    NHharm = ele.NHharm
    
    for i in 1:NHharm
        if 6*i > length(By)
            break
        end
        base = (i-1)*6
        HCw = float(By[base+2]) * Aw / Po
        kx = float(By[base+3]) * kw
        ky = float(By[base+4]) * kw
        kz = float(By[base+5]) * kw
        tz = float(By[base+6])

        sx = sin(kx*x)
        shy = sinh(ky*y)
        sz = sin(kz*Zw + tz)
        ay += HCw * (kw/kz) * (kx/ky) * sx * shy * sz

        cx = cos(kx*x)
        chy = cosh(ky*y)
        aypx += HCw * (kw/kz) * (kx/ky)^2 * cx * chy * sz
    end

    Bx = ele.Bx
    NVharm = ele.NVharm
    
    for i in 1:NVharm
        if 6*i > length(Bx)
            break
        end
        base = (i-1)*6
        VCw = float(Bx[base+2]) * Aw / Po
        kx = float(Bx[base+3]) * kw
        ky = float(Bx[base+4]) * kw
        kz = float(Bx[base+5]) * kw
        tz = float(Bx[base+6])

        chx = cosh(kx*x)
        cy = cos(ky*y)
        sz = sin(kz*Zw + tz)
        ay += VCw * (kw/kz) * chx * cy * sz

        shx = sinh(kx*x)
        syky = abs(ky/kw) > 1e-6 ? sin(ky*y)/ky : y * _sinc_taylor(ky*y)
        aypx += VCw * (kw/kz) * kx * shx * syky * sz
    end

    return ay, aypx
end

function _wig_map_2nd!(ele::WIGGLER, r::AbstractVector{Float64}, dl::Float64, Zw_ref::Ref{Float64}, Aw::Float64, Po::Float64)
    delta = r[6]
    dld = dl / (1.0 + delta)
    dl2 = 0.5 * dl
    dl2d = dl2 / (1.0 + delta)

    Zw_ref[] += dl2

    ay, aypx = _wig_ay_aypx!(ele, r, Zw_ref[], Aw, Po)
    r[2] -= aypx
    r[4] -= ay

    r[3] += dl2d * r[4]
    r[5] += 0.5 * dl2d * (r[4]^2) / (1.0 + delta)

    ay, aypx = _wig_ay_aypx!(ele, r, Zw_ref[], Aw, Po)
    r[2] += aypx
    r[4] += ay

    ax, axpy = _wig_ax_axpy!(ele, r, Zw_ref[], Aw, Po)
    r[2] -= ax
    r[4] -= axpy

    r[1] += dld * r[2]
    r[5] += 0.5 * dld * (r[2]^2) / (1.0 + delta)

    ax, axpy = _wig_ax_axpy!(ele, r, Zw_ref[], Aw, Po)
    r[2] += ax
    r[4] += axpy

    ay, aypx = _wig_ay_aypx!(ele, r, Zw_ref[], Aw, Po)
    r[2] -= aypx
    r[4] -= ay

    r[3] += dl2d * r[4]
    r[5] += 0.5 * dl2d * (r[4]^2) / (1.0 + delta)

    ay, aypx = _wig_ay_aypx!(ele, r, Zw_ref[], Aw, Po)
    r[2] += aypx
    r[4] += ay

    Zw_ref[] += dl2

    return nothing
end

# 4th-order pass: composition of three 2nd-order maps
function _wig_pass_4th!(ele::WIGGLER, r::AbstractVector{Float64})
    PN = ele.Nsteps
    Nw = Int(round(ele.len / ele.lw))
    Nstep = PN * max(1, Nw)
    dl = ele.lw / PN
    dl1 = dl * 1.3512071919596573
    dl0 = dl * -1.7024143839193146

    Po = sqrt((ele.energy / m_e)^2 - 1.0)
    C0 = 299792458.0  # m/s
    __E0_GeV = m_e * 1.0e-9  # Convert eV to GeV: 510998.95 eV = 0.00051099895 GeV
    TWOPI = 2.0 * π
    Aw = 1.0e-9 * C0 / __E0_GeV / TWOPI * ele.lw * ele.Bmax
    Zw_ref = Ref(0.0)

    for i in 1:Nstep
        _wig_map_2nd!(ele, r, dl1, Zw_ref, Aw, Po)
        _wig_map_2nd!(ele, r, dl0, Zw_ref, Aw, Po)
        _wig_map_2nd!(ele, r, dl1, Zw_ref, Aw, Po)
    end
    return nothing
end

# Synchrotron radiation
function _wig_B!(ele::WIGGLER, r::AbstractVector{Float64}, Zw::Float64, B::Vector{Float64})
    """Compute magnetic field [Bx, By, Bz] at particle location from analytic formula"""
    x  = r[1]
    y  = r[3]
    kw = 2π / ele.lw
    PB0 = ele.Bmax

    B[1] = 0.0  # Bx
    B[2] = 0.0  # By
    B[3] = 0.0  # Bz

    # Horizontal wiggler harmonics (By array)
    By = ele.By
    NHharm = ele.NHharm
    for i in 1:NHharm
        if 6i > length(By)
            break
        end
        base = (i - 1) * 6
        HCw_raw = float(By[base + 2])
        kx = float(By[base + 3]) * kw
        ky = float(By[base + 4]) * kw
        kz = float(By[base + 5]) * kw
        tz = float(By[base + 6])

        sx  = sin(kx * x)
        cx  = cos(kx * x)
        chy = cosh(ky * y)
        shy = sinh(ky * y)
        cz  = cos(kz * Zw + tz)
        sz  = sin(kz * Zw + tz)

        B[1] += PB0 * HCw_raw * (kx/ky) * sx * shy * cz   # Bx
        B[2] -= PB0 * HCw_raw * cx * chy * cz             # By
        B[3] += PB0 * HCw_raw * (kz/ky) * cx * shy * sz   # Bz
    end

    # Vertical wiggler harmonics (Bx array)
    Bx = ele.Bx
    NVharm = ele.NVharm
    for i in 1:NVharm
        if 6i > length(Bx)
            break
        end
        base = (i - 1) * 6
        VCw_raw = float(Bx[base + 2])
        kx = float(Bx[base + 3]) * kw
        ky = float(Bx[base + 4]) * kw
        kz = float(Bx[base + 5]) * kw
        tz = float(Bx[base + 6])

        shx = sinh(kx * x)
        chx = cosh(kx * x)
        sy  = sin(ky * y)
        cy  = cos(ky * y)
        cz  = cos(kz * Zw + tz)
        sz  = sin(kz * Zw + tz)

        B[1] += PB0 * VCw_raw * chx * cy * cz             # Bx
        B[2] -= PB0 * VCw_raw * (ky/kx) * shx * sy * cz   # By
        B[3] -= PB0 * VCw_raw * (kz/kx) * cy * shx * sz   # Bz
    end

    return nothing
end

# Apply synchrotron radiation energy-loss kicks (classical loss only).
function _wig_radiation_kicks!(r::AbstractVector{Float64},
                               B::Vector{Float64},
                               Po::Float64,
                               srCoef::Float64,
                               dl::Float64)
    # B^2 in T^2, using only transverse components (Bx, By) exactly like GWigRadiationKicks
    B2 = B[1]*B[1] + B[2]*B[2]
    if B2 == 0.0
        return nothing
    end

    # Beam rigidity H in T·m
    # C: H = 1.0e9 * pWig->Po * __E0 / C0
    # with __E0 (GeV) = m_e(eV)*1e-9  ⇒ 1e9*__E0 = m_e (eV)
    # => H = Po * m_e / C0
    C0 = 2.99792458e8
    H  = Po * m_e / C0

    # 1/rho^2
    irho2 = B2 / (H * H)

    # (1 + delta)^2 (delta is r[6] in JuTrack, X[4] in AT)
    delta   = r[6]
    dFactor = (1.0 + delta) * (1.0 + delta)

    # dDelta = -srCoef * (1+delta)^2 * (1/rho^2) * dl
    dDelta = -srCoef * dFactor * irho2 * dl

    # Update δ, px, py (X[4], X[1], X[3] in C)
    r[6] += dDelta
    scale = 1.0 + dDelta
    r[2] *= scale
    r[4] *= scale

    return nothing
end

# 4th-order integrator with synchrotron radiation (rad == 1)
function _wig_pass_4th_rad!(ele::WIGGLER, r::AbstractVector{Float64})
    # --- Integration parameters (match _wig_pass_4th! and GWigSymplecticRadPass) ---
    PN   = ele.Nsteps                          # Nstep per period
    Nw   = Int(round(ele.len / ele.lw))        # number of periods (≈ le/Lw)
    Nstep = PN * max(1, Nw)                    # Niter in C: Nstep*(le/Lw)
    SL   = ele.lw / PN                         # SL = Lw/Nstep
    dl1  = SL * 1.3512071919596573             # KICK1
    dl0  = SL * -1.7024143839193146            # KICK2

    # --- Precompute Po, Aw (same as in _wig_pass_4th!) ---
    # Po = sqrt(gamma^2 - 1), gamma = energy / rest_energy (in eV)
    gamma = ele.energy / m_e
    Po = sqrt(gamma^2 - 1.0)

    C0 = 2.99792458e8
    __E0_GeV = m_e * 1.0e-9
    TWOPI = 2.0 * π
    Aw = 1.0e-9 * C0 / __E0_GeV / TWOPI * ele.lw * ele.Bmax

    # --- SR coefficient (from GWigInit2) ---
    # Wig->srCoef = 2/3 * __RE * gamma^3
    srCoef = (2.0 / 3.0) * __RE * gamma^3

    # Internal longitudinal coordinate Zw (same as pWig->Zw)
    Zw_ref = Ref(0.0)

    # Magnetic field workspace
    B = [0.0, 0.0, 0.0]

    # ===== Entrance SR kick (at Zw = 0) =====
    # We mimic the sequence:
    #   GWigAx, GWigAy, GWigB;
    #   X[1] -= ax; X[3] -= ay;  // mechanical
    #   GWigRadiationKicks;
    #   X[1] += ax; X[3] += ay;  // back to canonical
    ax, axpy = _wig_ax_axpy!(ele, r, Zw_ref[], Aw, Po)
    ay, aypx = _wig_ay_aypx!(ele, r, Zw_ref[], Aw, Po)

    r[2] -= ax
    r[4] -= ay
    _wig_B!(ele, r, Zw_ref[], B)
    _wig_radiation_kicks!(r, B, Po, srCoef, SL)
    r[2] += ax
    r[4] += ay

    # ===== Main 4th-order integrator with SR kicks after each slice =====
    for i in 1:Nstep
        if i == 50
            nothing
        end
        _wig_map_2nd!(ele, r, dl1, Zw_ref, Aw, Po)
        _wig_map_2nd!(ele, r, dl0, Zw_ref, Aw, Po)
        _wig_map_2nd!(ele, r, dl1, Zw_ref, Aw, Po)

        ax, axpy = _wig_ax_axpy!(ele, r, Zw_ref[], Aw, Po)
        ay, aypx = _wig_ay_aypx!(ele, r, Zw_ref[], Aw, Po)

        r[2] -= ax
        r[4] -= ay
        _wig_B!(ele, r, Zw_ref[], B)
        _wig_radiation_kicks!(r, B, Po, srCoef, SL)
        r[2] += ax
        r[4] += ay
    end

    # Important: we **do not** call a GWigGauge-equivalent here.
    # Your non-SR _wig_pass_4th! already operates in the external
    # (canonical) coordinates that match AT at the element boundaries.
    # This SR version reduces exactly to _wig_pass_4th! if srCoef = 0.

    return nothing
end


# Main pass implementation using the WIGGLER struct; only 4th-order used
function pass!(ele::WIGGLER, r_in::Matrix{Float64}, num_particles::Int64, particles::Beam{Float64})
    T1 = ele.T1
    T2 = ele.T2
    R1 = ele.R1
    R2 = ele.R2
    lost_flags = particles.lost_flag

    for c in 1:num_particles
        r6 = @view r_in[c, :]
        if lost_flags[c] == 1
            continue
        end

        # entrance misalignment/rotation
        if !iszero(T1)
            addvv!(r6, T1)
        end
        if !iszero(R1)
            multmv!(r6, R1)
        end

        # 4th-order integrator (non-radiation only)
        if ele.rad == 1
            _wig_pass_4th_rad!(ele, r6)
        else
            _wig_pass_4th!(ele, r6)
        end

        # exit misalignment/rotation
        if !iszero(R2)
            multmv!(r6, R2)
        end
        if !iszero(T2)
            addvv!(r6, T2)
        end

        # check lost
        if check_lost(r6)
            lost_flags[c] = 1
        end
    end
    return nothing
end

function pass_P!(ele::WIGGLER, r_in::Matrix{Float64}, num_particles::Int64, particles::Beam{Float64})
    T1 = ele.T1
    T2 = ele.T2
    R1 = ele.R1
    R2 = ele.R2
    lost_flags = particles.lost_flag

    Threads.@threads for c in 1:num_particles
        r6 = @view r_in[c, :]
        if lost_flags[c] == 1
            continue
        end

        # entrance misalignment/rotation
        if !iszero(T1)
            addvv!(r6, T1)
        end
        if !iszero(R1)
            multmv!(r6, R1)
        end

        # 4th-order integrator (non-radiation only)
        _wig_pass_4th!(ele, r6)

        # exit misalignment/rotation
        if !iszero(R2)
            multmv!(r6, R2)
        end
        if !iszero(T2)
            addvv!(r6, T2)
        end

        # check lost
        if check_lost(r6)
            lost_flags[c] = 1
        end
    end
    return nothing
end

##########################################################################################
# High-order TPSA
# TPSA versions of helper functions for single TPSA vector (one particle)
function _wig_ax_axpy_TPSA!(ele::WIGGLER, r::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}, Zw::Float64, Aw::Float64, Po::Float64) where {T, TPS_Dim, Max_TPS_Degree}
    x = r[1]
    y = r[3]
    kw = 2π / ele.lw

    ax = 0.0
    axpy = 0.0

    By = ele.By
    NHharm = ele.NHharm
    
    for i in 1:NHharm
        if 6*i > length(By)
            break
        end
        base = (i-1)*6
        HCw = float(By[base+2]) * Aw / Po
        kx = float(By[base+3]) * kw
        ky = float(By[base+4]) * kw
        kz = float(By[base+5]) * kw
        tz = float(By[base+6])

        cx = cos(kx*x)
        chy = cosh(ky*y)
        sz = sin(kz*Zw + tz)
        ax += HCw * (kw/kz) * cx * chy * sz

        shy = sinh(ky*y)
        sxkx = abs(kx/kw) > 1e-6 ? sin(kx*x)/kx : x * _sinc_taylor(kx)
        axpy += HCw * (kw/kz) * ky * sxkx * shy * sz
    end

    Bx = ele.Bx
    NVharm = ele.NVharm
    
    for i in 1:NVharm
        if 6*i > length(Bx)
            break
        end
        base = (i-1)*6
        VCw = float(Bx[base+2]) * Aw / Po
        kx = float(Bx[base+3]) * kw
        ky = float(Bx[base+4]) * kw
        kz = float(Bx[base+5]) * kw
        tz = float(Bx[base+6])

        shx = sinh(kx*x)
        sy = sin(ky*y)
        sz = sin(kz*Zw + tz)
        ax += VCw * (kw/kz) * (ky/kx) * shx * sy * sz

        chx = cosh(kx*x)
        cy = cos(ky*y)
        axpy += VCw * (kw/kz) * (ky/kx)^2 * chx * cy * sz
    end

    return ax, axpy
end

function _wig_ay_aypx_TPSA!(ele::WIGGLER, r::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}, Zw::Float64, Aw::Float64, Po::Float64) where {T, TPS_Dim, Max_TPS_Degree}
    x = r[1]
    y = r[3]
    kw = 2π / ele.lw

    ay = 0.0
    aypx = 0.0

    By = ele.By
    NHharm = ele.NHharm
    
    for i in 1:NHharm
        if 6*i > length(By)
            break
        end
        base = (i-1)*6
        HCw = float(By[base+2]) * Aw / Po
        kx = float(By[base+3]) * kw
        ky = float(By[base+4]) * kw
        kz = float(By[base+5]) * kw
        tz = float(By[base+6])

        sx = sin(kx*x)
        shy = sinh(ky*y)
        sz = sin(kz*Zw + tz)
        ay += HCw * (kw/kz) * (kx/ky) * sx * shy * sz

        cx = cos(kx*x)
        chy = cosh(ky*y)
        aypx += HCw * (kw/kz) * (kx/ky)^2 * cx * chy * sz
    end

    Bx = ele.Bx
    NVharm = ele.NVharm
    
    for i in 1:NVharm
        if 6*i > length(Bx)
            break
        end
        base = (i-1)*6
        VCw = float(Bx[base+2]) * Aw / Po
        kx = float(Bx[base+3]) * kw
        ky = float(Bx[base+4]) * kw
        kz = float(Bx[base+5]) * kw
        tz = float(Bx[base+6])

        chx = cosh(kx*x)
        cy = cos(ky*y)
        sz = sin(kz*Zw + tz)
        ay += VCw * (kw/kz) * chx * cy * sz

        shx = sinh(kx*x)
        syky = abs(ky/kw) > 1e-6 ? sin(ky*y)/ky : y * _sinc_taylor(ky)
        aypx += VCw * (kw/kz) * kx * shx * syky * sz
    end

    return ay, aypx
end

function _wig_map_2nd_TPSA!(ele::WIGGLER, r::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}, dl::Float64, Zw_ref::Ref{Float64}, Aw::Float64, Po::Float64) where {T, TPS_Dim, Max_TPS_Degree}
    delta = r[6]
    dld = dl / (1.0 + delta)
    dl2 = 0.5 * dl
    dl2d = dl2 / (1.0 + delta)

    # Step1: half step in internal z
    Zw_ref[] += dl2

    # Step2: half drift in y
    ay, aypx = _wig_ay_aypx_TPSA!(ele, r, Zw_ref[], Aw, Po)
    r[2] -= aypx
    r[4] -= ay

    r[3] += dl2d * r[4]
    r[5] += 0.5 * dl2d * (r[4]^2) / (1.0 + delta)

    ay, aypx = _wig_ay_aypx_TPSA!(ele, r, Zw_ref[], Aw, Po)
    r[2] += aypx
    r[4] += ay

    # Step3: full drift in x
    ax, axpy = _wig_ax_axpy_TPSA!(ele, r, Zw_ref[], Aw, Po)
    r[2] -= ax
    r[4] -= axpy

    r[1] += dld * r[2]
    r[5] += 0.5 * dld * (r[2]^2) / (1.0 + delta)

    ax, axpy = _wig_ax_axpy_TPSA!(ele, r, Zw_ref[], Aw, Po)
    r[2] += ax
    r[4] += axpy

    # Step4: half drift in y
    ay, aypx = _wig_ay_aypx_TPSA!(ele, r, Zw_ref[], Aw, Po)
    r[2] -= aypx
    r[4] -= ay

    r[3] += dl2d * r[4]
    r[5] += 0.5 * dl2d * (r[4]^2) / (1.0 + delta)

    ay, aypx = _wig_ay_aypx_TPSA!(ele, r, Zw_ref[], Aw, Po)
    r[2] += aypx
    r[4] += ay

    # Step5: half step in internal z
    Zw_ref[] += dl2

    return nothing
end

function _wig_pass_4th_TPSA!(ele::WIGGLER, r::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}) where {T, TPS_Dim, Max_TPS_Degree}
    PN = ele.Nsteps
    Nw = Int(round(ele.len / ele.lw))
    Nstep = PN * max(1, Nw)
    dl = ele.lw / PN
    dl1 = dl * 1.3512071919596573
    dl0 = dl * -1.7024143839193146

    # precompute Aw and Po
    Po = sqrt((ele.energy / m_e)^2 - 1.0)
    C0 = 299792458.0
    __E0_GeV = m_e * 1.0e-9
    TWOPI = 2.0 * π
    Aw = 1.0e-9 * C0 / __E0_GeV / TWOPI * ele.lw * ele.Bmax
    Zw_ref = Ref(0.0)

    for i in 1:Nstep
        _wig_map_2nd_TPSA!(ele, r, dl1, Zw_ref, Aw, Po)
        _wig_map_2nd_TPSA!(ele, r, dl0, Zw_ref, Aw, Po)
        _wig_map_2nd_TPSA!(ele, r, dl1, Zw_ref, Aw, Po)
    end
    return nothing
end

function _wig_B_TPSA!(ele::WIGGLER, r::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}, Zw::Float64, B::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}) where {T, TPS_Dim, Max_TPS_Degree}
    """Compute magnetic field [Bx, By, Bz] at particle location for TPSA"""
    x  = r[1]
    y  = r[3]
    kw = 2π / ele.lw
    PB0 = ele.Bmax

    B[1] = CTPS(0.0, TPS_Dim, Max_TPS_Degree)  # Bx
    B[2] = CTPS(0.0, TPS_Dim, Max_TPS_Degree)  # By
    B[3] = CTPS(0.0, TPS_Dim, Max_TPS_Degree)  # Bz

    # Horizontal wiggler harmonics
    By = ele.By
    NHharm = ele.NHharm
    for i in 1:NHharm
        if 6i > length(By)
            break
        end
        base = (i - 1) * 6
        HCw_raw = float(By[base + 2])
        kx = float(By[base + 3]) * kw
        ky = float(By[base + 4]) * kw
        kz = float(By[base + 5]) * kw
        tz = float(By[base + 6])

        sx  = sin(kx * x)
        cx  = cos(kx * x)
        chy = cosh(ky * y)
        shy = sinh(ky * y)
        cz  = cos(kz * Zw + tz)
        sz  = sin(kz * Zw + tz)

        B[1] += PB0 * HCw_raw * (kx/ky) * sx * shy * cz
        B[2] -= PB0 * HCw_raw * cx * chy * cz
        B[3] += PB0 * HCw_raw * (kz/ky) * cx * shy * sz
    end

    # Vertical wiggler harmonics
    Bx = ele.Bx
    NVharm = ele.NVharm
    for i in 1:NVharm
        if 6i > length(Bx)
            break
        end
        base = (i - 1) * 6
        VCw_raw = float(Bx[base + 2])
        kx = float(Bx[base + 3]) * kw
        ky = float(Bx[base + 4]) * kw
        kz = float(Bx[base + 5]) * kw
        tz = float(Bx[base + 6])

        shx = sinh(kx * x)
        chx = cosh(kx * x)
        sy  = sin(ky * y)
        cy  = cos(ky * y)
        cz  = cos(kz * Zw + tz)
        sz  = sin(kz * Zw + tz)

        B[1] += PB0 * VCw_raw * chx * cy * cz
        B[2] -= PB0 * VCw_raw * (ky/kx) * shx * sy * cz
        B[3] -= PB0 * VCw_raw * (kz/kx) * cy * shx * sz
    end

    return nothing
end

function _wig_radiation_kicks_TPSA!(r::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}, 
                                    B::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}},
                                    Po::Float64,
                                    srCoef::Float64,
                                    dl::Float64) where {T, TPS_Dim, Max_TPS_Degree}
    """Apply synchrotron radiation kicks for TPSA"""
    B2 = B[1]*B[1] + B[2]*B[2]
    if B2 == 0.0
        return nothing
    end

    C0 = 2.99792458e8
    H  = Po * m_e / C0

    irho2 = B2 / (H * H)
    delta = r[6]
    dFactor = (1.0 + delta) * (1.0 + delta)
    dDelta = -srCoef * dFactor * irho2 * dl

    r[6] += dDelta
    scale = 1.0 + dDelta
    r[2] *= scale
    r[4] *= scale

    return nothing
end

function _wig_pass_4th_rad_TPSA!(ele::WIGGLER, r::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}) where {T, TPS_Dim, Max_TPS_Degree}
    """4th-order integrator with synchrotron radiation for TPSA"""
    PN   = ele.Nsteps
    Nw   = Int(round(ele.len / ele.lw))
    Nstep = PN * max(1, Nw)
    SL   = ele.lw / PN
    dl1  = SL * 1.3512071919596573
    dl0  = SL * -1.7024143839193146

    gamma = ele.energy / m_e
    Po = sqrt(gamma^2 - 1.0)
    C0 = 2.99792458e8
    __E0_GeV = m_e * 1.0e-9
    TWOPI = 2.0 * π
    Aw = 1.0e-9 * C0 / __E0_GeV / TWOPI * ele.lw * ele.Bmax
    srCoef = (2.0 / 3.0) * __RE * gamma^3

    Zw_ref = Ref(0.0)
    B = [CTPS(0.0, TPS_Dim, Max_TPS_Degree) for _ in 1:3]

    # Entrance SR kick
    ax, axpy = _wig_ax_axpy_TPSA!(ele, r, Zw_ref[], Aw, Po)
    ay, aypx = _wig_ay_aypx_TPSA!(ele, r, Zw_ref[], Aw, Po)
    r[2] -= ax
    r[4] -= ay
    _wig_B_TPSA!(ele, r, Zw_ref[], B)
    _wig_radiation_kicks_TPSA!(r, B, Po, srCoef, SL)
    r[2] += ax
    r[4] += ay

    # Main integration loop
    for i in 1:Nstep
        _wig_map_2nd_TPSA!(ele, r, dl1, Zw_ref, Aw, Po)
        _wig_map_2nd_TPSA!(ele, r, dl0, Zw_ref, Aw, Po)
        _wig_map_2nd_TPSA!(ele, r, dl1, Zw_ref, Aw, Po)

        ax, axpy = _wig_ax_axpy_TPSA!(ele, r, Zw_ref[], Aw, Po)
        ay, aypx = _wig_ay_aypx_TPSA!(ele, r, Zw_ref[], Aw, Po)
        r[2] -= ax
        r[4] -= ay
        _wig_B_TPSA!(ele, r, Zw_ref[], B)
        _wig_radiation_kicks_TPSA!(r, B, Po, srCoef, SL)
        r[2] += ax
        r[4] += ay
    end

    return nothing
end

function pass_TPSA!(ele::WIGGLER, r_in::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}; E0::Float64=0.0, m0::Float64=m_e) where {T, TPS_Dim, Max_TPS_Degree}
    T1 = ele.T1
    T2 = ele.T2
    R1 = ele.R1
    R2 = ele.R2

    # entrance misalignment/rotation
    if !iszero(T1)
        addvv!(r_in, T1)
    end
    if !iszero(R1)
        multmv!(r_in, R1)
    end

    # 4th-order integrator with or without radiation
    if ele.rad == 1
        _wig_pass_4th_rad_TPSA!(ele, r_in)
    else
        _wig_pass_4th_TPSA!(ele, r_in)
    end

    # exit misalignment/rotation
    if !iszero(R2)
        multmv!(r_in, R2)
    end
    if !iszero(T2)
        addvv!(r_in, T2)
    end

    return nothing
end

##########################################################################################
# DTPSAD 
function _sinc_taylor(x::DTPSAD{N, T}) where {N, T <: Number}
    x2 = x*x
    return 1.0 - x2/6.0*(1.0 - x2/20.0*(1.0 - x2/42.0*(1.0 - x2/72.0)))
end

function _wig_ax_axpy_DTPSAD!(ele::WIGGLER{DTPSAD{N, T}}, r::SubArray, Zw::DTPSAD{N, T}, Aw::DTPSAD{N, T}, Po::DTPSAD{N, T}) where {N, T <: Number}
    x = r[1]
    y = r[3]
    kw = 2π / ele.lw

    ax = zero(DTPSAD{N, T})
    axpy = zero(DTPSAD{N, T})

    By = ele.By
    NHharm = ele.NHharm
    
    for i in 1:NHharm
        if 6*i > length(By)
            break
        end
        base = (i-1)*6
        HCw = By[base+2] * Aw / Po
        kx = By[base+3] * kw
        ky = By[base+4] * kw
        kz = By[base+5] * kw
        tz = By[base+6]

        cx = cos(kx*x)
        chy = cosh(ky*y)
        sz = sin(kz*Zw + tz)
        ax += HCw * (kw/kz) * cx * chy * sz

        shy = sinh(ky*y)
        sxkx = abs(kx/kw) > 1e-6 ? sin(kx*x)/kx : x * _sinc_taylor(kx)
        axpy += HCw * (kw/kz) * ky * sxkx * shy * sz
    end

    Bx = ele.Bx
    NVharm = ele.NVharm
    
    for i in 1:NVharm
        if 6*i > length(Bx)
            break
        end
        base = (i-1)*6
        VCw = Bx[base+2] * Aw / Po
        kx = Bx[base+3] * kw
        ky = Bx[base+4] * kw
        kz = Bx[base+5] * kw
        tz = Bx[base+6]

        shx = sinh(kx*x)
        sy = sin(ky*y)
        sz = sin(kz*Zw + tz)
        ax += VCw * (kw/kz) * (ky/kx) * shx * sy * sz

        chx = cosh(kx*x)
        cy = cos(ky*y)
        axpy += VCw * (kw/kz) * (ky/kx)^2 * chx * cy * sz
    end

    return ax, axpy
end

function _wig_ay_aypx_DTPSAD!(ele::WIGGLER{DTPSAD{N, T}}, r::SubArray, Zw::DTPSAD{N, T}, Aw::DTPSAD{N, T}, Po::DTPSAD{N, T}) where {N, T <: Number}
    x = r[1]
    y = r[3]
    kw = 2π / ele.lw

    ay = zero(DTPSAD{N, T})
    aypx = zero(DTPSAD{N, T})

    By = ele.By
    NHharm = ele.NHharm
    
    for i in 1:NHharm
        if 6*i > length(By)
            break
        end
        base = (i-1)*6
        HCw = By[base+2] * Aw / Po
        kx = By[base+3] * kw
        ky = By[base+4] * kw
        kz = By[base+5] * kw
        tz = By[base+6]

        sx = sin(kx*x)
        shy = sinh(ky*y)
        sz = sin(kz*Zw + tz)
        ay += HCw * (kw/kz) * (kx/ky) * sx * shy * sz

        cx = cos(kx*x)
        chy = cosh(ky*y)
        aypx += HCw * (kw/kz) * (kx/ky)^2 * cx * chy * sz
    end

    Bx = ele.Bx
    NVharm = ele.NVharm
    
    for i in 1:NVharm
        if 6*i > length(Bx)
            break
        end
        base = (i-1)*6
        VCw = Bx[base+2] * Aw / Po
        kx = Bx[base+3] * kw
        ky = Bx[base+4] * kw
        kz = Bx[base+5] * kw
        tz = Bx[base+6]

        chx = cosh(kx*x)
        cy = cos(ky*y)
        sz = sin(kz*Zw + tz)
        ay += VCw * (kw/kz) * chx * cy * sz

        shx = sinh(kx*x)
        syky = abs(ky/kw) > 1e-6 ? sin(ky*y)/ky : y * _sinc_taylor(ky)
        aypx += VCw * (kw/kz) * kx * shx * syky * sz
    end

    return ay, aypx
end

function _wig_map_2nd_DTPSAD!(ele::WIGGLER{DTPSAD{N, T}}, r::SubArray, dl::DTPSAD{N, T}, Zw_ref::Ref{DTPSAD{N, T}}, Aw::DTPSAD{N, T}, Po::DTPSAD{N, T}) where {N, T <: Number}
    delta = r[6]
    dld = dl / (1.0 + delta)
    dl2 = 0.5 * dl
    dl2d = dl2 / (1.0 + delta)

    # Step1: half step in internal z
    Zw_ref[] += dl2

    # Step2: half drift in y
    ay, aypx = _wig_ay_aypx_DTPSAD!(ele, r, Zw_ref[], Aw, Po)
    r[2] -= aypx
    r[4] -= ay

    r[3] += dl2d * r[4]
    r[5] += 0.5 * dl2d * (r[4]^2) / (1.0 + delta)

    ay, aypx = _wig_ay_aypx_DTPSAD!(ele, r, Zw_ref[], Aw, Po)
    r[2] += aypx
    r[4] += ay

    # Step3: full drift in x
    ax, axpy = _wig_ax_axpy_DTPSAD!(ele, r, Zw_ref[], Aw, Po)
    r[2] -= ax
    r[4] -= axpy

    r[1] += dld * r[2]
    r[5] += 0.5 * dld * (r[2]^2) / (1.0 + delta)

    ax, axpy = _wig_ax_axpy_DTPSAD!(ele, r, Zw_ref[], Aw, Po)
    r[2] += ax
    r[4] += axpy

    # Step4: half drift in y
    ay, aypx = _wig_ay_aypx_DTPSAD!(ele, r, Zw_ref[], Aw, Po)
    r[2] -= aypx
    r[4] -= ay

    r[3] += dl2d * r[4]
    r[5] += 0.5 * dl2d * (r[4]^2) / (1.0 + delta)

    ay, aypx = _wig_ay_aypx_DTPSAD!(ele, r, Zw_ref[], Aw, Po)
    r[2] += aypx
    r[4] += ay

    # Step5: half step in internal z
    Zw_ref[] += dl2

    return nothing
end

function _wig_pass_4th_DTPSAD!(ele::WIGGLER{DTPSAD{N, T}}, r::SubArray) where {N, T <: Number}
    PN = ele.Nsteps
    Nw = Int(round(ele.len / ele.lw))
    Nstep = PN * max(1, Nw)
    dl = ele.lw / PN
    dl1 = dl * 1.3512071919596573
    dl0 = dl * -1.7024143839193146

    # precompute Aw and Po (in DTPSAD type)
    Po = sqrt((ele.energy / m_e)^2 - 1.0)
    C0 = 299792458.0
    __E0_GeV = m_e * 1.0e-9
    TWOPI = 2.0 * π
    Aw = 1.0e-9 * C0 / __E0_GeV / TWOPI * ele.lw * ele.Bmax
    Zw_ref = Ref(zero(r[1]))

    for i in 1:Nstep
        _wig_map_2nd_DTPSAD!(ele, r, dl1, Zw_ref, Aw, Po)
        _wig_map_2nd_DTPSAD!(ele, r, dl0, Zw_ref, Aw, Po)
        _wig_map_2nd_DTPSAD!(ele, r, dl1, Zw_ref, Aw, Po)
    end
    return nothing
end

function _wig_B_DTPSAD!(ele::WIGGLER{DTPSAD{N, T}}, r::SubArray, Zw::DTPSAD{N, T}, B::Vector{DTPSAD{N, T}}) where {N, T <: Number}
    """Compute magnetic field [Bx, By, Bz] at particle location for DTPSAD"""
    x  = r[1]
    y  = r[3]
    kw = 2π / ele.lw
    PB0 = ele.Bmax

    B[1] = zero(x)
    B[2] = zero(x)
    B[3] = zero(x)

    # Horizontal wiggler harmonics
    By = ele.By
    NHharm = ele.NHharm
    for i in 1:NHharm
        if 6i > length(By)
            break
        end
        base = (i - 1) * 6
        HCw_raw = By[base + 2]
        kx = By[base + 3] * kw
        ky = By[base + 4] * kw
        kz = By[base + 5] * kw
        tz = By[base + 6]

        sx  = sin(kx * x)
        cx  = cos(kx * x)
        chy = cosh(ky * y)
        shy = sinh(ky * y)
        cz  = cos(kz * Zw + tz)
        sz  = sin(kz * Zw + tz)

        B[1] += PB0 * HCw_raw * (kx/ky) * sx * shy * cz
        B[2] -= PB0 * HCw_raw * cx * chy * cz
        B[3] += PB0 * HCw_raw * (kz/ky) * cx * shy * sz
    end

    # Vertical wiggler harmonics
    Bx = ele.Bx
    NVharm = ele.NVharm
    for i in 1:NVharm
        if 6i > length(Bx)
            break
        end
        base = (i - 1) * 6
        VCw_raw = Bx[base + 2]
        kx = Bx[base + 3] * kw
        ky = Bx[base + 4] * kw
        kz = Bx[base + 5] * kw
        tz = Bx[base + 6]

        shx = sinh(kx * x)
        chx = cosh(kx * x)
        sy  = sin(ky * y)
        cy  = cos(ky * y)
        cz  = cos(kz * Zw + tz)
        sz  = sin(kz * Zw + tz)

        B[1] += PB0 * VCw_raw * chx * cy * cz
        B[2] -= PB0 * VCw_raw * (ky/kx) * shx * sy * cz
        B[3] -= PB0 * VCw_raw * (kz/kx) * cy * shx * sz
    end

    return nothing
end

function _wig_radiation_kicks_DTPSAD!(r::SubArray, 
                                      B::Vector{DTPSAD{N, T}},
                                      Po::DTPSAD{N, T},
                                      srCoef::DTPSAD{N, T},
                                      dl::DTPSAD{N, T}) where {N, T <: Number}
    """Apply synchrotron radiation kicks for DTPSAD"""
    B2 = B[1]*B[1] + B[2]*B[2]
    if B2 == 0.0
        return nothing
    end

    C0 = 2.99792458e8
    m_e_val = 510998.95  # eV
    H  = Po * m_e_val / C0

    irho2 = B2 / (H * H)
    delta = r[6]
    dFactor = (1.0 + delta) * (1.0 + delta)
    dDelta = -srCoef * dFactor * irho2 * dl

    r[6] += dDelta
    scale = 1.0 + dDelta
    r[2] *= scale
    r[4] *= scale

    return nothing
end

function _wig_pass_4th_rad_DTPSAD!(ele::WIGGLER{DTPSAD{N, T}}, r::SubArray) where {N, T <: Number}
    """4th-order integrator with synchrotron radiation for DTPSAD"""
    PN   = ele.Nsteps
    Nw   = Int(round(ele.len / ele.lw))
    Nstep = PN * max(1, Nw)
    SL   = ele.lw / PN
    dl1  = SL * 1.3512071919596573
    dl0  = SL * -1.7024143839193146

    gamma = ele.energy / m_e
    Po = sqrt(gamma^2 - 1.0)
    C0 = 299792458.0
    __E0_GeV = m_e * 1.0e-9
    TWOPI = 2.0 * π
    Aw = 1.0e-9 * C0 / __E0_GeV / TWOPI * ele.lw * ele.Bmax
    srCoef = (2.0 / 3.0) * __RE * gamma^3

    Zw_ref = Ref(zero(r[1]))
    B = [zero(r[1]), zero(r[1]), zero(r[1])]

    # Entrance SR kick
    ax, axpy = _wig_ax_axpy_DTPSAD!(ele, r, Zw_ref[], Aw, Po)
    ay, aypx = _wig_ay_aypx_DTPSAD!(ele, r, Zw_ref[], Aw, Po)
    r[2] -= ax
    r[4] -= ay
    _wig_B_DTPSAD!(ele, r, Zw_ref[], B)
    _wig_radiation_kicks_DTPSAD!(r, B, Po, srCoef, SL)
    r[2] += ax
    r[4] += ay

    # Main integration loop
    for i in 1:Nstep
        _wig_map_2nd_DTPSAD!(ele, r, dl1, Zw_ref, Aw, Po)
        _wig_map_2nd_DTPSAD!(ele, r, dl0, Zw_ref, Aw, Po)
        _wig_map_2nd_DTPSAD!(ele, r, dl1, Zw_ref, Aw, Po)

        ax, axpy = _wig_ax_axpy_DTPSAD!(ele, r, Zw_ref[], Aw, Po)
        ay, aypx = _wig_ay_aypx_DTPSAD!(ele, r, Zw_ref[], Aw, Po)
        r[2] -= ax
        r[4] -= ay
        _wig_B_DTPSAD!(ele, r, Zw_ref[], B)
        _wig_radiation_kicks_DTPSAD!(r, B, Po, srCoef, SL)
        r[2] += ax
        r[4] += ay
    end

    return nothing
end

function pass!(ele::WIGGLER{DTPSAD{N, T}}, rin::Matrix{DTPSAD{N, T}}, npart::Int, particles::Beam{DTPSAD{N, T}}) where {N, T <: Number}
    @inbounds for c in 1:npart
        lost_flag = particles.lost_flag[c]
        if lost_flag == 1
            continue
        end
        r6 = @view rin[c, :]

        # entrance misalignment/rotation
        if !all(iszero, ele.T1)
            addvv!(r6, ele.T1)
        end
        if !all(iszero, ele.R1)
            multmv!(r6, ele.R1)
        end

        # 4th-order integrator with or without radiation
        if ele.rad == 1
            _wig_pass_4th_rad_DTPSAD!(ele, r6)
        else
            _wig_pass_4th_DTPSAD!(ele, r6)
        end

        # exit misalignment/rotation
        if !all(iszero, ele.R2)
            multmv!(r6, ele.R2)
        end
        if !all(iszero, ele.T2)
            addvv!(r6, ele.T2)
        end

        lost_flag = check_lost_GTPSA(r6) ? 1 : lost_flag
        particles.lost_flag[c] = lost_flag
    end
    return nothing
end

