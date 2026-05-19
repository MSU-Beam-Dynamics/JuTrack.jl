# Frequency Map Analysis (FMA) for JuTrack.
# This code uses JuTrack for tracking and PyCall to access the 
# nafflib Python library for tune analysis.
# To use, nafflib and matplotlib must be installed in the Python environment of PyCall.

using PyCall
using LinearAlgebra

const np = pyimport("numpy")
const nafflib = pyimport("nafflib")

# Particle generation
function fma_grid(;
    xmin=-5e-3,
    xmax=5e-3,
    ymin=0.0,
    ymax=3e-3,
    nx::Int=101,
    ny::Int=31,
    px0=0.0,
    py0=0.0,
    dp0=0.0,
    ct0=0.0,
)
    particles = zeros(nx * ny, 6)

    k = 0
    for ix in 1:nx
        x = nx == 1 ? 0.5 * (xmin + xmax) :
            xmin + (xmax - xmin) * (ix - 1) / (nx - 1)

        for iy in 1:ny
            y = ny == 1 ? 0.5 * (ymin + ymax) :
                ymin + (ymax - ymin) * (iy - 1) / (ny - 1)

            k += 1
            particles[k, 1] = x
            particles[k, 2] = px0
            particles[k, 3] = y
            particles[k, 4] = py0
            particles[k, 5] = ct0
            particles[k, 6] = dp0
        end
    end

    return particles
end


# Tracking function
function tracking_output_to_arrays(rout, nturns::Int, npart::Int)
    X  = fill(NaN, nturns, npart)
    PX = fill(NaN, nturns, npart)
    Y  = fill(NaN, nturns, npart)
    PY = fill(NaN, nturns, npart)

    @inbounds for it in 1:nturns
        r = rout[it]

        for ip in 1:npart
            X[it, ip]  = r[ip, 1]
            PX[it, ip] = r[ip, 2]
            Y[it, ip]  = r[ip, 3]
            PY[it, ip] = r[ip, 4]
        end
    end

    return X, PX, Y, PY
end


function track_fma(RING, particles::Matrix{Float64}, nturns::Int; energy=6e9)
    particles0 = copy(particles)
    npart = size(particles0, 1)

    beam = Beam(copy(particles0), energy=energy)

    # Expected output:
    # rout[turn][particle, coordinate]
    rout = pringpass!(RING, beam, nturns, true)

    X, PX, Y, PY = tracking_output_to_arrays(rout, nturns, npart)

    return (
        particles = particles0,
        X = X,
        PX = PX,
        Y = Y,
        PY = PY,
        nturns = nturns,
        npart = npart,
        energy = energy,
    )
end


# Optional Courant-Snyder normalization
function cs_normalize_trajectory(traj, twi)
    sx = sqrt(Float64(twi.betax))
    sy = sqrt(Float64(twi.betay))

    X  = traj.X
    PX = traj.PX
    Y  = traj.Y
    PY = traj.PY

    Xn  = X ./ sx
    PXn = Float64(twi.alphax) .* X ./ sx .+ sx .* PX

    Yn  = Y ./ sy
    PYn = Float64(twi.alphay) .* Y ./ sy .+ sy .* PY

    return Xn, PXn, Yn, PYn
end

# NAFF tune analysis
function demean_columns!(A::Matrix{Float64})
    T, N = size(A)

    @inbounds for j in 1:N
        μ = 0.0

        for i in 1:T
            μ += A[i, j]
        end

        μ /= T

        for i in 1:T
            A[i, j] -= μ
        end
    end

    return A
end


function good_naff_columns(X::AbstractMatrix, P::AbstractMatrix; min_rms=1e-14)
    @assert size(X) == size(P)

    T, N = size(X)
    good = trues(N)

    @inbounds for j in 1:N
        mx = 0.0
        mp = 0.0

        for i in 1:T
            x = X[i, j]
            p = P[i, j]

            if !isfinite(x) || !isfinite(p)
                good[j] = false
                break
            end

            mx += x
            mp += p
        end

        if !good[j]
            continue
        end

        mx /= T
        mp /= T

        s = 0.0
        for i in 1:T
            dx = X[i, j] - mx
            dp = P[i, j] - mp
            s += dx * dx + dp * dp
        end

        good[j] = sqrt(s / T) > min_rms
    end

    return good
end


function tune_mod(q; fold_to_half=false)
    if !isfinite(q)
        return NaN
    end

    ν = mod(Float64(q), 1.0)

    if fold_to_half
        ν = min(ν, 1.0 - ν)
    end

    return ν
end


function tune_distance(ν1, ν2; fold_to_half=false)
    if !isfinite(ν1) || !isfinite(ν2)
        return NaN
    end

    d = abs(Float64(ν2) - Float64(ν1))

    if fold_to_half
        return d
    else
        return min(d, 1.0 - d)
    end
end


function to_py_particles_by_turns(A::Matrix{Float64})
    # Julia: turns × particles
    # nafflib.multiparticle_tunes: particles × turns
    return np.array(permutedims(A); dtype=np.float64, order="C")
end


function naff_tunes(
    X::AbstractMatrix,
    P::AbstractMatrix;
    window_order::Int=2,
    window_type::String="hann",
    fold_to_half::Bool=false,
    subtract_mean::Bool=true,
    min_rms::Float64=1e-14,
    processes::Union{Nothing,Int}=nothing,
)
    @assert size(X) == size(P)

    T, N = size(X)
    @assert T >= 16 "Need at least 16 turns for NAFF."

    tunes = fill(NaN, N)

    good = good_naff_columns(X, P; min_rms=min_rms)
    goodidx = findall(good)

    if isempty(goodidx)
        return tunes
    end

    Xg = Matrix{Float64}(X[:, goodidx])
    Pg = Matrix{Float64}(P[:, goodidx])

    if subtract_mean
        demean_columns!(Xg)
        demean_columns!(Pg)
    end

    x_py = to_py_particles_by_turns(Xg)
    p_py = to_py_particles_by_turns(Pg)

    q_py = if processes === nothing
        nafflib.multiparticle_tunes(
            x_py,
            p_py;
            window_order=window_order,
            window_type=window_type,
        )
    else
        nafflib.multiparticle_tunes(
            x_py,
            p_py;
            window_order=window_order,
            window_type=window_type,
            processes=processes,
        )
    end

    q = Vector{Float64}(q_py)

    for (ii, j) in enumerate(goodidx)
        tunes[j] = tune_mod(q[ii]; fold_to_half=fold_to_half)
    end

    return tunes
end


# FMA
function compute_fma_from_tracking(
    traj;
    normalize::Bool=true,
    twiss=nothing,
    window_order::Int=2,
    window_type::String="hann",
    fold_to_half::Bool=false,
    subtract_mean::Bool=true,
    min_rms::Float64=1e-14,
    diffusion_floor::Float64=-40.0,
    processes::Union{Nothing,Int}=nothing,
)
    nturns = traj.nturns
    npart = traj.npart

    @assert iseven(nturns) "Use an even number of turns."
    @assert nturns >= 32 "Use at least 32 turns."

    if normalize
        @assert twiss !== nothing "Pass twiss=periodicEdwardsTengTwiss(...) when normalize=true."
        X, PX, Y, PY = cs_normalize_trajectory(traj, twiss)
    else
        X, PX, Y, PY = traj.X, traj.PX, traj.Y, traj.PY
    end

    h = div(nturns, 2)

    X1  = @view X[1:h, :]
    PX1 = @view PX[1:h, :]
    Y1  = @view Y[1:h, :]
    PY1 = @view PY[1:h, :]

    X2  = @view X[h+1:end, :]
    PX2 = @view PX[h+1:end, :]
    Y2  = @view Y[h+1:end, :]
    PY2 = @view PY[h+1:end, :]

    nux1 = naff_tunes(
        X1, PX1;
        window_order=window_order,
        window_type=window_type,
        fold_to_half=fold_to_half,
        subtract_mean=subtract_mean,
        min_rms=min_rms,
        processes=processes,
    )

    nuy1 = naff_tunes(
        Y1, PY1;
        window_order=window_order,
        window_type=window_type,
        fold_to_half=fold_to_half,
        subtract_mean=subtract_mean,
        min_rms=min_rms,
        processes=processes,
    )

    nux2 = naff_tunes(
        X2, PX2;
        window_order=window_order,
        window_type=window_type,
        fold_to_half=fold_to_half,
        subtract_mean=subtract_mean,
        min_rms=min_rms,
        processes=processes,
    )

    nuy2 = naff_tunes(
        Y2, PY2;
        window_order=window_order,
        window_type=window_type,
        fold_to_half=fold_to_half,
        subtract_mean=subtract_mean,
        min_rms=min_rms,
        processes=processes,
    )

    dnu_x = fill(NaN, npart)
    dnu_y = fill(NaN, npart)
    dnu   = fill(NaN, npart)
    diffusion = fill(NaN, npart)

    @inbounds for j in 1:npart
        dx = tune_distance(nux1[j], nux2[j]; fold_to_half=fold_to_half)
        dy = tune_distance(nuy1[j], nuy2[j]; fold_to_half=fold_to_half)

        if isfinite(dx) && isfinite(dy)
            dnu_x[j] = dx
            dnu_y[j] = dy
            dnu[j] = hypot(dx, dy)

            val = dx^2 + dy^2
            diffusion[j] = log10(max(val, 10.0^diffusion_floor))
        end
    end

    rows = Vector{NamedTuple}(undef, npart)

    @inbounds for j in 1:npart
        rows[j] = (
            x0 = traj.particles[j, 1],
            px0 = traj.particles[j, 2],
            y0 = traj.particles[j, 3],
            py0 = traj.particles[j, 4],

            nux1 = nux1[j],
            nuy1 = nuy1[j],
            nux2 = nux2[j],
            nuy2 = nuy2[j],

            # Use first-half tune as the plotted tune.
            nux = nux1[j],
            nuy = nuy1[j],

            dnu_x = dnu_x[j],
            dnu_y = dnu_y[j],
            dnu = dnu[j],
            diffusion = diffusion[j],
        )
    end

    return (
        rows = rows,
        nux1 = nux1,
        nuy1 = nuy1,
        nux2 = nux2,
        nuy2 = nuy2,
        dnu_x = dnu_x,
        dnu_y = dnu_y,
        dnu = dnu,
        diffusion = diffusion,
    )
end

function rowfield(rows, name::Symbol)
    v = fill(NaN, length(rows))

    for i in eachindex(rows)
        if name in propertynames(rows[i])
            v[i] = Float64(getproperty(rows[i], name))
        end
    end

    return v
end


function get_rows(obj)
    if obj isa NamedTuple && (:rows in keys(obj))
        return obj.rows
    else
        return obj
    end
end


function finite_mask(vecs...)
    n = length(vecs[1])
    mask = trues(n)

    for v in vecs
        @assert length(v) == n
        mask .&= isfinite.(v)
    end

    return mask
end


function resonance_style(order)
    if order == 1
        return ("black", "-", 0.8, 1.1)
    elseif order == 2
        return ("red", "--", 0.55, 0.9)
    elseif order == 3
        return ("blue", "-.", 0.45, 0.75)
    elseif order == 4
        return ("purple", ":", 0.4, 0.7)
    else
        return ("gray", ":", 0.25, 0.5)
    end
end


function gcd3(a, b, c)
    return gcd(gcd(abs(a), abs(b)), abs(c))
end


function plot_resonance_lines!(
    ax;
    xlim=(0.0, 1.0),
    ylim=(0.0, 1.0),
    orders=1:4,
    primitive=true,
)
    plt = pyimport("matplotlib.pyplot")

    max_order = maximum(collect(orders))
    xmin, xmax = xlim
    ymin, ymax = ylim

    for n in -max_order:max_order
        for m in -max_order:max_order
            if n == 0 && m == 0
                continue
            end

            order = abs(n) + abs(m)

            if !(order in orders)
                continue
            end

            # Avoid duplicate sign copies.
            if n < 0 || (n == 0 && m < 0)
                continue
            end

            vals = (
                n * xmin + m * ymin,
                n * xmin + m * ymax,
                n * xmax + m * ymin,
                n * xmax + m * ymax,
            )

            kmin = ceil(Int, minimum(vals) - 1e-12)
            kmax = floor(Int, maximum(vals) + 1e-12)

            for k in kmin:kmax
                if primitive && gcd3(n, m, k) != 1
                    continue
                end

                color, ls, alpha, lw = resonance_style(order)

                if m == 0
                    x = k / n
                    if xmin <= x <= xmax
                        ax.plot(
                            [x, x],
                            [ymin, ymax],
                            color=color,
                            linestyle=ls,
                            alpha=alpha,
                            linewidth=lw,
                        )
                    end

                elseif n == 0
                    y = k / m
                    if ymin <= y <= ymax
                        ax.plot(
                            [xmin, xmax],
                            [y, y],
                            color=color,
                            linestyle=ls,
                            alpha=alpha,
                            linewidth=lw,
                        )
                    end

                else
                    xs = collect(range(xmin, xmax, length=500))
                    ys = (k .- n .* xs) ./ m
                    good = (ys .>= ymin) .& (ys .<= ymax)

                    if any(good)
                        ax.plot(
                            xs[good],
                            ys[good],
                            color=color,
                            linestyle=ls,
                            alpha=alpha,
                            linewidth=lw,
                        )
                    end
                end
            end
        end
    end

    return nothing
end

function plot_fma(
    obj;
    figsize=(11, 5),
    s=8,
    cmap="jet",
    initial_unit=:mm,
    xlim_initial=nothing,
    ylim_initial=nothing,
    tune_xlim=(0.0, 1.0),
    tune_ylim=(0.0, 1.0),
    diffusion_clim=nothing,
    resonance_lines=true,
    resonance_orders=1:4,
    show_all_initial_particles=true,
    filepath=nothing,
    dpi=300,
    show=true,
)
    plt = pyimport("matplotlib.pyplot")

    rows = get_rows(obj)

    if initial_unit == :mm
        scale = 1e3
        unit = "mm"
    elseif initial_unit == :m
        scale = 1.0
        unit = "m"
    else
        error("initial_unit must be :mm or :m")
    end

    x0 = rowfield(rows, :x0) .* scale
    y0 = rowfield(rows, :y0) .* scale
    nux = rowfield(rows, :nux)
    nuy = rowfield(rows, :nuy)
    diffusion = rowfield(rows, :diffusion)

    fig = plt.figure(figsize=figsize)

    # --------------------------------------------------------
    # Left: initial coordinate map
    # --------------------------------------------------------
    ax1 = plt.subplot(1, 2, 1)

    mask_xy = finite_mask(x0, y0)
    mask_diff = finite_mask(x0, y0, diffusion)

    if show_all_initial_particles && any(mask_xy)
        ax1.scatter(
            x0[mask_xy],
            y0[mask_xy];
            s=max(1.0, 0.35 * s),
            color="0.82",
            marker="s",
            label="tracked initial points",
            zorder=1,
        )
    end

    if any(mask_diff)
        if diffusion_clim === nothing
            sc1 = ax1.scatter(
                x0[mask_diff],
                y0[mask_diff];
                c=diffusion[mask_diff],
                s=s,
                marker="s",
                cmap=cmap,
                zorder=2,
            )
        else
            sc1 = ax1.scatter(
                x0[mask_diff],
                y0[mask_diff];
                c=diffusion[mask_diff],
                s=s,
                cmap=cmap,
                marker="s",
                vmin=diffusion_clim[1],
                vmax=diffusion_clim[2],
                zorder=2,
            )
        end

        plt.colorbar(sc1, ax=ax1, label="log10(Δνx² + Δνy²)")
    else
        ax1.text(
            0.5,
            0.5,
            "No finite diffusion data",
            ha="center",
            va="center",
            transform=ax1.transAxes,
        )
    end

    ax1.set_xlabel("x₀ [$unit]")
    ax1.set_ylabel("y₀ [$unit]")
    # ax1.set_title("Initial coordinate map")

    if xlim_initial !== nothing
        ax1.set_xlim(xlim_initial[1], xlim_initial[2])
    end

    if ylim_initial !== nothing
        ax1.set_ylim(ylim_initial[1], ylim_initial[2])
    end

    # --------------------------------------------------------
    # Right: frequency map
    # --------------------------------------------------------
    ax2 = plt.subplot(1, 2, 2)

    mask_tune = finite_mask(nux, nuy, diffusion)

    if any(mask_tune)
        if diffusion_clim === nothing
            sc2 = ax2.scatter(
                nux[mask_tune],
                nuy[mask_tune];
                c=diffusion[mask_tune],
                s=s,
                marker="s",
                cmap=cmap,
                zorder=2,
            )
        else
            sc2 = ax2.scatter(
                nux[mask_tune],
                nuy[mask_tune];
                c=diffusion[mask_tune],
                s=s,
                marker="s",
                cmap=cmap,
                vmin=diffusion_clim[1],
                vmax=diffusion_clim[2],
                zorder=2,
            )
        end

        plt.colorbar(sc2, ax=ax2, label="log10(Δνx² + Δνy²)")
    else
        ax2.text(
            0.5,
            0.5,
            "No finite tune data",
            ha="center",
            va="center",
            transform=ax2.transAxes,
        )
    end

    if resonance_lines
        plot_resonance_lines!(
            ax2;
            xlim=tune_xlim,
            ylim=tune_ylim,
            orders=resonance_orders,
        )
    end

    ax2.set_xlabel("νx")
    ax2.set_ylabel("νy")
    # ax2.set_title("Frequency map")
    ax2.set_xlim(tune_xlim[1], tune_xlim[2])
    ax2.set_ylim(tune_ylim[1], tune_ylim[2])
    # ax2.set_aspect("equal", adjustable="box")

    plt.tight_layout()

    if filepath !== nothing
        plt.savefig(filepath, dpi=dpi, bbox_inches="tight")
    end

    if show
        plt.show()
    end

    return fig
end