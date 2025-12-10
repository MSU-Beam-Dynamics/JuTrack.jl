using FFTW
using PyCall
using Statistics

function dominant_frequency_phase_complex(x::AbstractVector{<:Real}, p::AbstractVector{<:Real})
    N = length(x)
    @assert N == length(p) "x and p must have same length"
    if any(!isfinite, x) || any(!isfinite, p)
        return (NaN, NaN, NaN)
    end
    # Remove means to suppress DC
    xz = x .- mean(x)
    pz = p .- mean(p)
    z = ComplexF64.(xz, .-pz) # x - i p

    w = @. 0.5 - 0.5*cos(2π*(0:(N-1))/(N-1))
    zw = z .* w
    L = 1 << (ceil(Int, log2(N)) + 3)
    buf = (L == N) ? zw : vcat(zw, zeros(ComplexF64, L - N))
    F = fft(buf)

    mag = abs.(F)
    upper = L >>> 1
    kmax = argmax(@view mag[2:upper]) + 1

    k = kmax
    if 2 ≤ k ≤ (upper - 1)
        ALPHA, BETA, GAMMA = mag[k-1], mag[k], mag[k+1]
        denom = (ALPHA - 2BETA + GAMMA)
        δ = iszero(denom) ? 0.0 : 0.5*(ALPHA - GAMMA)/denom
        kHAT = k + δ
    else
        kHAT = k
    end


    nu = kHAT / L
    ϕ = angle(F[k])
    A = 2*mag[k]/sum(w)
    return nu, A, ϕ
end

"""
    compute_tunes_from_tracking(series; include_x=true, include_y=true)

Compute tunes (nux, nuy) and amplitudes from time series. 
Returns `(nux, nuy, Ax, Ay, φx, φxp, φy, φyp)` where phases come from
x vs xp and y vs yp and are used to adjust to the correct half-plane.
"""
function compute_tunes_from_tracking(series; include_x::Bool=true, include_y::Bool=true, use_complex::Bool=true, fold_to_half::Bool=true)
    x, xp, y, yp = series.x, series.xp, series.y, series.yp
    N = length(x)
    @assert N ≥ 16 "Need at least 16 turns for tune extraction"

    nux = 0.0; nuy = 0.0; Ax = 0.0; Ay = 0.0
    φx = 0.0; φxp = 0.0; φy = 0.0; φyp = 0.0

    if include_x
        if use_complex
            nux, Ax, φx = dominant_frequency_phase_complex(x, xp)
        else
            nux, Ax, φx = dominant_frequency_phase(x)
            _,  _,  φxp = dominant_frequency_phase(xp)
            nux = adjust_tune_half_plane(nux, φx, φxp)
        end
        if fold_to_half
            nux = min(nux, 1 - nux)
        end
    end
    if include_y
        if use_complex
            nuy, Ay, φy = dominant_frequency_phase_complex(y, yp)
        else
            nuy, Ay, φy = dominant_frequency_phase(y)
            _,  _,  φyp = dominant_frequency_phase(yp)
            nuy = adjust_tune_half_plane(nuy, φy, φyp)
        end
        if fold_to_half
            nuy = min(nuy, 1 - nuy)
        end
    end
    return nux, nuy, Ax, Ay, φx, φxp, φy, φyp
end

function adjust_tune_half_plane(nu::Real, ϕ0::Real, ϕ1::Real)
    Δ = ϕ0 - ϕ1
    # unwrap across π
    if abs(Δ) > π
        if ϕ0 < ϕ1
            ϕ0 += 2π
        else
            ϕ1 += 2π
        end
    end
    return (ϕ0 < ϕ1) ? nu : (1 - nu)
end

pack_segment(x, xp, y, yp) = (x=x, xp=xp, y=y, yp=yp)

function dominant_frequency_phase(sig::AbstractVector{<:Real})
    N = length(sig)
    @assert N ≥ 16 "Need at least 16 samples"
    if any(!isfinite, sig)
        return (NaN, NaN, NaN)
    end
    # Hann window
    w = @. 0.5 - 0.5*cos(2π*(0:(N-1))/(N-1))
    sw = sig .* w

    # Zero-pad by 8x
    L = 1 << (ceil(Int, log2(N)) + 3)
    buf = (L == N) ? sw : vcat(sw, zeros(Float64, L - N))
    F = fft(buf)

    mag = abs.(F)
    # skip DC, search up to Nyquist
    upper = L >>> 1
    kmax = argmax(@view mag[2:upper]) + 1

    k = kmax
    if 2 ≤ k ≤ (upper - 1)
        ALF, BET, GAM = mag[k-1], mag[k], mag[k+1]
        denom = (ALF - 2BET + GAM)
        δ = iszero(denom) ? 0.0 : 0.5*(ALF - GAM)/denom
        KHAT = k + δ
    else
        KHAT = k
    end

    nu = KHAT / L
    ϕ = angle(F[k])
    A = 2*mag[k]/sum(w)
    return nu, A, ϕ
end

"""
    fma_map_from_segments(seg1, seg2=nothing; include_x=true, include_y=true)

Frequency-map analysis for *many particles* using user-supplied trajectories.

`seg1` and `seg2` are NamedTuples of the form `(x, xp, y, yp)` where each field
is a matrix of size `(turns, N)` — rows are turns, columns are particles.

If `seg2` is provided, diffusion metrics are computed between the two segments
(similar to the C code: Δnux, Δnuy, √(Δnux²+Δnuy²), ΔAx, ΔAy, and log10(Δnux²+Δnuy²)).

Returns a Vector of NamedTuples with fields:
  `:nux, :nuy, :Ax, :Ay, :dnu_x, :dnu_y, :dnu, :dA_x, :dA_y, :diffusion`.
"""
function fma_map_from_segments(seg1::NamedTuple{(:x,:xp,:y,:yp)},
                               seg2::Union{Nothing,NamedTuple{(:x,:xp,:y,:yp)}}=nothing;
                               include_x::Bool=true, include_y::Bool=true)
    x1, xp1, y1, yp1 = seg1.x, seg1.xp, seg1.y, seg1.yp
    @assert size(x1) == size(xp1) == size(y1) == size(yp1) "All seg1 arrays must have same size (turns, N)"
    T, N = size(x1)

    has2 = seg2 !== nothing
    if has2
        x2, xp2, y2, yp2 = seg2.x, seg2.xp, seg2.y, seg2.yp
        @assert size(x2) == (T, N) && size(x2)==size(xp2)==size(y2)==size(yp2) "seg2 must match seg1 in size"
    end

    out = Vector{NamedTuple}(undef, N)
    @inbounds for j in 1:N
        s1_ok = all(isfinite, view(x1, :, j)) && all(isfinite, view(xp1, :, j)) &&
                all(isfinite, view(y1, :, j)) && all(isfinite, view(yp1, :, j))

        nux1 = NaN; nuy1 = NaN; Ax1 = NaN; Ay1 = NaN
        dnux = 0.0; dnuy = 0.0; dnu = 0.0; dAx = 0.0; dAy = 0.0; diffusion = 0.0

        if s1_ok
            nux1, nuy1, Ax1, Ay1, _, _, _, _ = compute_tunes_from_tracking((
                x  = view(x1, :, j),
                xp = view(xp1, :, j),
                y  = view(y1, :, j),
                yp = view(yp1, :, j),
            ); include_x, include_y)
        end

        if has2
            x2, xp2, y2, yp2 = seg2.x, seg2.xp, seg2.y, seg2.yp
            s2_ok = all(isfinite, view(x2, :, j)) && all(isfinite, view(xp2, :, j)) &&
                    all(isfinite, view(y2, :, j)) && all(isfinite, view(yp2, :, j))
            if s1_ok && s2_ok
                nux2, nuy2, Ax2, Ay2, _, _, _, _ = compute_tunes_from_tracking((
                    x  = view(x2, :, j),
                    xp = view(xp2, :, j),
                    y  = view(y2, :, j),
                    yp = view(yp2, :, j),
                ); include_x, include_y)
                dnux = abs(nux2 - nux1)
                dnuy = abs(nuy2 - nuy1)
                dnu  = hypot(nux2 - nux1, nuy2 - nuy1)
                dAx = abs(Ax1 - Ax2)
                dAy = abs(Ay1 - Ay2)
                val = (nux2 - nux1)^2 + (nuy2 - nuy1)^2
                diffusion = val > 0 ? log10(val) : NaN
            else
                nux1 = NaN; nuy1 = NaN; Ax1 = NaN; Ay1 = NaN; dnux = NaN; dnuy = NaN; dnu = NaN; dAx = NaN; dAy = NaN; diffusion = NaN
            end
        end

        out[j] = (
            nux = nux1, nuy = nuy1,
            Ax = Ax1, Ay = Ay1,
            dnu_x = dnux, dnu_y = dnuy, dnu = dnu,
            dA_x = dAx, dA_y = dAy,
            diffusion = diffusion,
        )
    end
    return out
end

"""
    FMA(RING, nturns; 
    xmin=-3e-3+1e-6, xmax=3e-3+1e-6, ymin=1e-6, ymax=3e-3+1e-6,
    nx=61, ny=31,
    na=21, ns=31, amax=3e-3+1e-6, 
    energy=2e9, sampling_method="grid", normalize_coordinates=false)
Perform Frequency Map Analysis (FMA) on a given ring lattice.
# Arguments
- RING::Vector{<:AbstractElement{Float64}}: a ring lattice
- nturns::Int: number of turns to track
- xmin, xmax, ymin, ymax::Float64: bounds for grid sampling
- nx, ny::Int: number of points in x and y for grid sampling
- na, ns::Int: number of angles and steps for radial sampling
- amax::Float64: maximum amplitude for radial sampling
- energy::Float64: beam energy [eV]
- sampling_method::String: "grid" or "radial"
- normalize_coordinates::Bool: whether to normalize coordinates using Twiss parameters
# Returns
- rows::Vector{NamedTuple}: FMA results for each particle
"""
function FMA(RING, nturns; 
    xmin=-3e-3+1e-6, xmax=3e-3+1e-6, ymin=1e-6, ymax=3e-3+1e-6,
    nx=61, ny=31,
    na=21, ns=31, amax=3e-3+1e-6, 
    energy=2e9, sampling_method="grid", normalize_coordinates=false)

    if sampling_method == "grid"
        # sample particles uniformly in x-y plane
        particles = zeros(nx * ny, 6)
        for i in 1:nx
            for j in 1:ny
                particle_idx = (i - 1) * ny + j
                x = xmin + (xmax - xmin) * (i - 1) / (nx - 1)
                y = ymin + (ymax - ymin) * (j - 1) / (ny - 1)
                particles[particle_idx, 1] = x  # x position
                particles[particle_idx, 3] = y  # y position
            end
        end
    elseif sampling_method == "radial"
        # sample particles in a radial pattern
        particles = zeros(na * ns, 6)
        particle_idx = 0
        for i in 1:na
            for j in 1:ns
                particle_idx += 1
                angle = π * (i - 1) / na
                amplitude = amax * (j - 1) / (ns - 1)
                x = amplitude * cos(angle)
                y = amplitude * sin(angle)
                particles[particle_idx, 1] = x  # x position
                particles[particle_idx, 3] = y  # y position
            end
        end
    else
        error("Unknown sampling method: $sampling_method")
    end

    beam = Beam(particles, energy=energy)
    rout = pringpass!(RING, beam, nturns, true)  # (turns, N) array

    if normalize_coordinates
        twi = periodicEdwardsTengTwiss(RING, 0.0, 0, E0=energy)
        betax = twi.betax
        betay = twi.betay
        alphax = twi.alphax
        alphay = twi.alphay
        n_particles = size(particles, 1)
        X = zeros(nturns, n_particles)
        XP = zeros(nturns, n_particles)
        Y = zeros(nturns, n_particles)
        YP = zeros(nturns, n_particles)
        for i in 1:nturns
            for j in 1:n_particles
                X[i, j] = rout[i][j, 1] *sqrt(betax) + rout[i][j, 2]*alphax/sqrt(betax)
                XP[i, j] = rout[i][j, 2] /sqrt(betax)
                Y[i, j] = rout[i][j, 3] *sqrt(betay) + rout[i][j, 4]*alphay/sqrt(betay)
                YP[i, j] = rout[i][j, 4] /sqrt(betay)
            end
        end
    else
        n_particles = size(particles, 1)
        X = zeros(nturns, n_particles)
        XP = zeros(nturns, n_particles)
        Y = zeros(nturns, n_particles)
        YP = zeros(nturns, n_particles)
        for i in 1:nturns
            for j in 1:n_particles
                X[i, j] = rout[i][j, 1]
                XP[i, j] = rout[i][j, 2]
                Y[i, j] = rout[i][j, 3]
                YP[i, j] = rout[i][j, 4]
            end
        end
    end
    half_turns = div(nturns, 2)
    seg1 = pack_segment(X[1:half_turns, :], XP[1:half_turns, :], Y[1:half_turns, :], YP[1:half_turns, :])
    seg2 = pack_segment(X[half_turns+1:end, :], XP[half_turns+1:end, :], Y[half_turns+1:end, :], YP[half_turns+1:end, :])
    rows = fma_map_from_segments(seg1, seg2; include_x=true, include_y=true)
    # add initial conditions to each row
    for j in 1:n_particles
        rows[j] = merge(rows[j], (x0=particles[j, 1], px0=particles[j, 2], y0=particles[j, 3], py0=particles[j, 4]))
    end
    return rows
end

function plot_resonance_line_extended(n, m, k, color, style, alpha, linewidth)
    plt = try
        pyimport("matplotlib.pyplot")
    catch _
        error("PyCall and Matplotlib are required for plotting. Please install them first.")
    end
    nu_range_extended = range(0.0, 1.0, length=100)
    if m == 0 && n != 0  # Vertical lines: nux = k/n
        if 0 <= k/n <= 1.0
            plt.axvline(x=k/n, color=color, linestyle=style, alpha=alpha, linewidth=linewidth)
        end
    elseif n == 0 && m != 0  # Horizontal lines: nuy = k/m  
        if 0 <= k/m <= 1.0
            plt.axhline(y=k/m, color=color, linestyle=style, alpha=alpha, linewidth=linewidth)
        end
    else  # nuy = (k - n*nux)/m
        nuy_line = (k .- n .* nu_range_extended) ./ m
        valid_idx = (nuy_line .>= 0.0) .& (nuy_line .<= 1.0) .& 
                   (nu_range_extended .>= 0.0) .& (nu_range_extended .<= 1.0)
        if any(valid_idx)
            plt.plot(nu_range_extended[valid_idx], nuy_line[valid_idx], color=color, 
                    linestyle=style, alpha=alpha, linewidth=linewidth)
        end
    end
end

"""
    plot_fma(rows; 
    figsize=(10,4), s=10, 
    x_min=-0.003, x_max=0.003, y_min=0.0, y_max=0.003,
    resonance_lines=true, resonance_orders=[1,2,3,4],
    filepath="fma_plot.png")
Plot Frequency Map Analysis (FMA) results. 
The plot function requires PyCall and Matplotlib installed in the associated Python environment.
# Arguments
- rows::Vector{NamedTuple}: FMA results from `FMA` function
- figsize::Tuple{Int,Int}: figure size
- s::Int: marker size
- x_min, x_max, y_min, y_max::Float64: axis limits
- resonance_lines::Bool: whether to plot resonance lines
- resonance_orders::Vector{Int}: resonance orders to plot
- filepath::String: output file path for saving the plot
"""
function plot_fma(rows; 
                    figsize=(10,4), s=10, 
                    x_min=-0.003, x_max=0.003, y_min=0.0, y_max=0.003,
                    resonance_lines=true, resonance_orders=[1,2,3,4],
                    filepath="fma_plot.png")
    nux = [r.nux for r in rows]
    nuy = [r.nuy for r in rows]
    diffusion = [r.diffusion for r in rows]
    x = [r.x0 for r in rows]
    y = [r.y0 for r in rows]

    # First order resonances 
    resonances_1st = [
        # n*nux + m*nuy = 0
        (1, 0, 0),   # nux = 0
        (0, 1, 0),   # nuy = 0
        (1, 1, 0),   # nux + nuy = 0
        (1, -1, 0),  # nux - nuy = 0
        (-1, 1, 0),  # -nux + nuy = 0

        # n*nux + m*nuy = 1
        (1, 0, 1),   # nux = 1
        (0, 1, 1),   # nuy = 1  
        (1, 1, 1),   # nux + nuy = 1
        (1, -1, 1),  # nux - nuy = 1
        (-1, 1, 1),  # -nux + nuy = 1
        (2, 0, 1),   # 2nux = 1 (half-integer)
        (0, 2, 1),   # 2nuy = 1 (half-integer)
    ]

    # Second order resonances
    resonances_2nd = [
        # n*nux + m*nuy = 0
        (2, 0, 0),   # 2nux = 0
        (0, 2, 0),   # 2nuy = 0
        (2, 1, 0),   # 2nux + nuy = 0
        (1, 2, 0),   # nux + 2nuy = 0
        (2, -1, 0),  # 2nux - nuy = 0
        (-1, 2, 0),  # -nux + 2nuy = 0
        
        # n*nux + m*nuy = 1
        (2, 1, 1),   # 2nux + nuy = 1
        (1, 2, 1),   # nux + 2nuy = 1
        (2, -1, 1),  # 2nux - nuy = 1
        (-1, 2, 1),  # -nux + 2nuy = 1
        
        # n*nux + m*nuy = 2
        (2, 0, 2),   # 2nux = 2 (integer)
        (0, 2, 2),   # 2nuy = 2 (integer)
        (1, 1, 2),   # nux + nuy = 2
    ]

    # Third order resonances  
    resonances_3rd = [
        # n*nux + m*nuy = 0
        (3, 0, 0),   # 3nux = 0
        (0, 3, 0),   # 3nuy = 0
        (3, 1, 0), (1, 3, 0), (3, -1, 0), (-1, 3, 0),
        (3, 2, 0), (2, 3, 0), (3, -2, 0), (-2, 3, 0),
        
        # n*nux + m*nuy = 1  
        (3, 0, 1),   # 3nux = 1 (third-integer)
        (0, 3, 1),   # 3nuy = 1 (third-integer)
        (3, 1, 1), (1, 3, 1), (3, -1, 1), (-1, 3, 1),
        (3, 2, 1), (2, 3, 1), (3, -2, 1), (-2, 3, 1),
    ]

    # Fourth order
    resonances_4th = [
        (4, 0, 1), (0, 4, 1),  # 4nux = 1, 4nuy = 1 (quarter-integer)
        (2, 2, 1),             # 2nux + 2nuy = 1
    ]

    plt = try
        pyimport("matplotlib.pyplot")
    catch _
        error("PyCall and Matplotlib are required for plotting. Please install them first.")
    end
    plt.figure(figsize=figsize)
    plt.subplot(1,2,1)
    plt.scatter(x, y, c=diffusion, s=s, cmap="jet")
    plt.colorbar(label="log10(Δnux² + Δnuy²)")
    plt.xlabel("x [m]")
    plt.ylabel("y [m]")
    plt.xlim(x_min, x_max)
    plt.ylim(y_min, y_max)
    plt.subplot(1,2,2)
    sc = plt.scatter(nux, nuy, c=diffusion, cmap="jet", s=s)
    plt.colorbar(sc, label="log10(Δnux² + Δnuy²)")
    plt.xlabel("nux")
    plt.ylabel("nuy")
    if resonance_lines
        if 1 in resonance_orders
            for (n, m, k) in resonances_1st
                plot_resonance_line_extended(n, m, k, "black", "-", 0.8, 1.5)
            end
        end
        if 2 in resonance_orders
            for (n, m, k) in resonances_2nd
                plot_resonance_line_extended(n, m, k, "red", "--", 0.5, 1.0)
            end
        end
        if 3 in resonance_orders
            for (n, m, k) in resonances_3rd
                plot_resonance_line_extended(n, m, k, "blue", "-.", 0.4, 0.8)
            end
        end
        if 4 in resonance_orders
            for (n, m, k) in resonances_4th
                plot_resonance_line_extended(n, m, k, "purple", ":", 0.3, 0.6)
            end
        end
    end
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.tight_layout()
    plt.savefig(filepath, dpi=300)
    plt.show()
end