using FFTW
using StaticArrays: SVector

@inline function _sc2p5d_macro_size(particles::Beam{T}) where {T}
    if particles.nmacro <= 0
        return zero(T)
    end
    return convert(T, particles.np / particles.nmacro)
end

@inline _sc2p5d_gridx(xmin, dx, ix::Int) = xmin + (ix - 1) * dx
@inline _sc2p5d_gridy(ymin, dy, iy::Int) = ymin + (iy - 1) * dy
@inline _sc2p5d_gridz(zmin, dz, iz::Int) = zmin + (iz - 0.5) * dz
@inline _sc2p5d_primal(x::Float64) = x
@inline _sc2p5d_primal(x::DTPSAD{N, T}) where {N, T <: Number} = Float64(x.val)

function _sc2p5d_extrema(r_in::Matrix{Float64}, lost_flags::Vector{Int})
    x_min = Inf
    x_max = -Inf
    y_min = Inf
    y_max = -Inf
    z_min = Inf
    z_max = -Inf
    n_active = 0
    @inbounds for c in 1:size(r_in, 1)
        if isone(lost_flags[c]) || isnan(r_in[c, 1])
            continue
        end
        x = r_in[c, 1]
        y = r_in[c, 3]
        z = r_in[c, 5]
        x_min = min(x_min, x)
        x_max = max(x_max, x)
        y_min = min(y_min, y)
        y_max = max(y_max, y)
        z_min = min(z_min, z)
        z_max = max(z_max, z)
        n_active += 1
    end
    return x_min, x_max, y_min, y_max, z_min, z_max, n_active
end

function _sc2p5d_adjust_bounds(x_min::Float64, x_max::Float64, y_min::Float64, y_max::Float64, xy_ratio::Float64)
    xy_ratio_beam = (x_max - x_min) / (y_max - y_min)
    if xy_ratio_beam > xy_ratio
        center = (y_max + y_min) / 2.0
        width = ((y_max - y_min) * (xy_ratio_beam / xy_ratio)) / 2.0
        y_min = center - width
        y_max = center + width
    else
        center = (x_max + x_min) / 2.0
        width = ((x_max - x_min) / (xy_ratio_beam / xy_ratio)) / 2.0
        x_min = center - width
        x_max = center + width
    end
    return x_min, x_max, y_min, y_max
end

function _sc2p5d_adjust_bounds(x_min::DTPSAD{N, T}, x_max::DTPSAD{N, T},
    y_min::DTPSAD{N, T}, y_max::DTPSAD{N, T}, xy_ratio::DTPSAD{N, T}) where {N, T <: Number}
    xy_ratio_beam = (x_max - x_min) / (y_max - y_min)
    if xy_ratio_beam > xy_ratio
        center = (y_max + y_min) / 2.0
        width = ((y_max - y_min) * (xy_ratio_beam / xy_ratio)) / 2.0
        y_min = center - width
        y_max = center + width
    else
        center = (x_max + x_min) / 2.0
        width = ((x_max - x_min) / (xy_ratio_beam / xy_ratio)) / 2.0
        x_min = center - width
        x_max = center + width
    end
    return x_min, x_max, y_min, y_max
end

function _sc2p5d_grid2d_ind_frac(coord::Float64, cmin::Float64, dc::Float64, n::Int)
    ind = trunc(Int, (coord - cmin) / dc + 0.5) + 1
    ind = clamp(ind, 2, n - 1)
    frac = (coord - (cmin + (ind - 1) * dc)) / dc
    return ind, frac
end

function _sc2p5d_bin_2d!(grid::Matrix{Float64}, x::Float64, y::Float64, weight::Float64,
    x_min::Float64, x_max::Float64, y_min::Float64, y_max::Float64)
    if x < x_min || x > x_max || y < y_min || y > y_max
        return nothing
    end
    nx, ny = size(grid)
    dx = (x_max - x_min) / (nx - 1)
    dy = (y_max - y_min) / (ny - 1)
    ix, x_frac = _sc2p5d_grid2d_ind_frac(x, x_min, dx, nx)
    iy, y_frac = _sc2p5d_grid2d_ind_frac(y, y_min, dy, ny)
    x_frac2 = x_frac * x_frac
    y_frac2 = y_frac * y_frac
    wxm = 0.5 * (0.25 - x_frac + x_frac2)
    wx0 = 0.75 - x_frac2
    wxp = 0.5 * (0.25 + x_frac + x_frac2)
    wym = 0.5 * (0.25 - y_frac + y_frac2)
    wy0 = 0.75 - y_frac2
    wyp = 0.5 * (0.25 + y_frac + y_frac2)
    @inbounds begin
        grid[ix - 1, iy - 1] += wxm * wym * weight
        grid[ix - 1, iy    ] += wxm * wy0 * weight
        grid[ix - 1, iy + 1] += wxm * wyp * weight
        grid[ix    , iy - 1] += wx0 * wym * weight
        grid[ix    , iy    ] += wx0 * wy0 * weight
        grid[ix    , iy + 1] += wx0 * wyp * weight
        grid[ix + 1, iy - 1] += wxp * wym * weight
        grid[ix + 1, iy    ] += wxp * wy0 * weight
        grid[ix + 1, iy + 1] += wxp * wyp * weight
    end
    return nothing
end

function _sc2p5d_grid2d_gradient(grid::Matrix{Float64}, x::Float64, y::Float64,
    x_min::Float64, x_max::Float64, y_min::Float64, y_max::Float64)
    nx, ny = size(grid)
    dx = (x_max - x_min) / (nx - 1)
    dy = (y_max - y_min) / (ny - 1)
    ix, x_frac = _sc2p5d_grid2d_ind_frac(x, x_min, dx, nx)
    iy, y_frac = _sc2p5d_grid2d_ind_frac(y, y_min, dy, ny)
    x_frac2 = x_frac * x_frac
    y_frac2 = y_frac * y_frac
    wxm = 0.5 * (0.25 - x_frac + x_frac2)
    wx0 = 0.75 - x_frac2
    wxp = 0.5 * (0.25 + x_frac + x_frac2)
    wym = 0.5 * (0.25 - y_frac + y_frac2)
    wy0 = 0.75 - y_frac2
    wyp = 0.5 * (0.25 + y_frac + y_frac2)
    dwxm = -0.5 + x_frac
    dwx0 = -2.0 * x_frac
    dwxp = 0.5 + x_frac
    dwym = -0.5 + y_frac
    dwy0 = -2.0 * y_frac
    dwyp = 0.5 + y_frac
    @inbounds begin
        ex = dwxm * wym * grid[ix - 1, iy - 1] +
             dwxm * wy0 * grid[ix - 1, iy    ] +
             dwxm * wyp * grid[ix - 1, iy + 1] +
             dwx0 * wym * grid[ix    , iy - 1] +
             dwx0 * wy0 * grid[ix    , iy    ] +
             dwx0 * wyp * grid[ix    , iy + 1] +
             dwxp * wym * grid[ix + 1, iy - 1] +
             dwxp * wy0 * grid[ix + 1, iy    ] +
             dwxp * wyp * grid[ix + 1, iy + 1]
        ey = wxm * dwym * grid[ix - 1, iy - 1] +
             wxm * dwy0 * grid[ix - 1, iy    ] +
             wxm * dwyp * grid[ix - 1, iy + 1] +
             wx0 * dwym * grid[ix    , iy - 1] +
             wx0 * dwy0 * grid[ix    , iy    ] +
             wx0 * dwyp * grid[ix    , iy + 1] +
             wxp * dwym * grid[ix + 1, iy - 1] +
             wxp * dwy0 * grid[ix + 1, iy    ] +
             wxp * dwyp * grid[ix + 1, iy + 1]
    end
    return ex / dx, ey / dy
end

function _sc2p5d_grid1d_get_ind_and_wz(z::Float64, z_min::Float64, dz::Float64, nz::Int)
    if nz <= 1
        return 1, 1, 1.0, 0.0
    end
    ind0 = trunc(Int, ((z - z_min) / dz) - 0.5)
    z_frac = (z - _sc2p5d_gridz(z_min, dz, ind0 + 1)) / dz
    if ind0 == 0 && z_frac < 0.0
        return nz, 1, -z_frac, 1.0 + z_frac
    end
    if ind0 == nz - 1
        if z_frac > 0.0
            return nz, 1, 1.0 - z_frac, z_frac
        end
        return nz - 1, nz, -z_frac, 1.0 + z_frac
    end
    return ind0 + 1, ind0 + 2, 1.0 - z_frac, z_frac
end

function _sc2p5d_bin_1d!(grid::Vector{Float64}, z::Float64, weight::Float64, z_min::Float64, z_max::Float64)
    if z < z_min || z > z_max
        return nothing
    end
    nz = length(grid)
    if nz == 0
        return nothing
    end
    dz = (z_max - z_min) / nz
    iz0, izp, wz0, wzp = _sc2p5d_grid1d_get_ind_and_wz(z, z_min, dz, nz)
    if nz > 1
        @inbounds begin
            grid[iz0] += wz0 * weight
            grid[izp] += wzp * weight
        end
    else
        @inbounds grid[1] += weight
    end
    return nothing
end

function _sc2p5d_grid1d_value(grid::Vector{Float64}, z::Float64, z_min::Float64, z_max::Float64)
    if z < z_min || z > z_max
        return 0.0
    end
    nz = length(grid)
    if nz == 0
        return 0.0
    end
    if nz == 1
        return grid[1]
    end
    dz = (z_max - z_min) / nz
    iz0, izp, wz0, wzp = _sc2p5d_grid1d_get_ind_and_wz(z, z_min, dz, nz)
    @inbounds return wz0 * grid[iz0] + wzp * grid[izp]
end

function _sc2p5d_build_green_fft!(ele::SPACECHARGE2P5D, dx::Float64, dy::Float64)
    if isapprox(ele.green_dx, dx; atol=1e-15, rtol=1e-12) &&
       isapprox(ele.green_dy, dy; atol=1e-15, rtol=1e-12)
        return nothing
    end
    nx2 = 2 * ele.xsize
    ny2 = 2 * ele.ysize
    greens = zeros(ComplexF64, nx2, ny2)
    half_x = div(nx2, 2)
    half_y = div(ny2, 2)
    @inbounds for iy in 0:half_y
        r_y = iy * dy
        for ix in 0:half_x
            r_x = ix * dx
            r2 = r_x * r_x + r_y * r_y
            val = (ix == 0 && iy == 0) ? 0.0 : -0.5 * log(r2)
            greens[ix + 1, iy + 1] = ComplexF64(val, 0.0)
        end
        for ix in half_x + 1:nx2 - 1
            greens[ix + 1, iy + 1] = greens[nx2 - ix + 1, iy + 1]
        end
    end
    @inbounds for iy in half_y + 1:ny2 - 1
        for ix in 0:nx2 - 1
            greens[ix + 1, iy + 1] = greens[ix + 1, ny2 - iy + 1]
        end
    end
    ele.green_fft .= fft(greens)
    ele.green_dx = dx
    ele.green_dy = dy
    return nothing
end

function _sc2p5d_find_potential!(ele::SPACECHARGE2P5D, dx::Float64, dy::Float64)
    _sc2p5d_build_green_fft!(ele, dx, dy)
    nx = ele.xsize
    ny = ele.ysize
    rho_pad = zeros(ComplexF64, 2 * nx, 2 * ny)
    @views rho_pad[1:nx, 1:ny] .= complex.(ele.rho_grid, 0.0)
    phi_pad = ifft(fft(rho_pad) .* ele.green_fft)
    @views ele.phi_grid .= real.(phi_pad[1:nx, 1:ny])
    return nothing
end

function _sc2p5d_calculate_long_derivative!(ele::SPACECHARGE2P5D, z_min::Float64, z_max::Float64)
    nz = ele.zsize
    dz = (z_max - z_min) / nz
    fill!(ele.z_deriv_grid, 0.0)
    if nz == 1
        return nothing
    elseif nz == 2
        z_mid = 0.5 * (_sc2p5d_gridz(z_min, dz, 1) + _sc2p5d_gridz(z_min, dz, 2))
        z_deriv = (_sc2p5d_grid1d_value(ele.z_grid, z_mid + 0.5 * dz, z_min, z_max) -
                   _sc2p5d_grid1d_value(ele.z_grid, z_mid - 0.5 * dz, z_min, z_max)) / dz
        ele.z_deriv_grid[1] = z_deriv / dz
        ele.z_deriv_grid[2] = z_deriv / dz
        return nothing
    end

    n_avg = min(ele.long_avg_n, nz)
    s = zeros(Float64, 5, 2)
    @inbounds for iz in 1:nz
        z = _sc2p5d_gridz(z_min, dz, iz)
        i_start = iz - div(n_avg, 2)
        if i_start < 1
            i_start = 1
        end
        if (i_start + n_avg - 1) > nz
            i_start = nz - n_avg + 1
        end
        i_stop = i_start + n_avg - 1
        fill!(s, 0.0)
        for i in i_start:i_stop
            x = i - i_start + 1
            y = ele.z_grid[i]
            xpow = 1.0
            for j in 0:4
                s[j + 1, 1] += xpow
                s[j + 1, 2] += xpow * y
                xpow *= x
            end
        end
        det = s[1, 1] * s[3, 1] * s[5, 1] -
              s[2, 1] * s[2, 1] * s[5, 1] -
              s[1, 1] * s[4, 1] * s[4, 1] +
              2.0 * s[2, 1] * s[3, 1] * s[4, 1] -
              s[3, 1] * s[3, 1] * s[3, 1]
        if abs(det) <= eps(Float64)
            ele.z_deriv_grid[iz] = 0.0
            continue
        end
        a = (s[1, 2] * s[2, 1] * s[4, 1] - s[2, 2] * s[1, 1] * s[4, 1] - s[1, 2] * s[3, 1] * s[3, 1] +
             s[2, 2] * s[2, 1] * s[3, 1] + s[3, 2] * s[1, 1] * s[3, 1] - s[3, 2] * s[2, 1] * s[2, 1]) / det
        b = (s[2, 2] * s[1, 1] * s[5, 1] - s[1, 2] * s[2, 1] * s[5, 1] + s[1, 2] * s[3, 1] * s[4, 1] -
             s[3, 2] * s[1, 1] * s[4, 1] - s[2, 2] * s[3, 1] * s[3, 1] + s[3, 2] * s[2, 1] * s[3, 1]) / det
        deriv = 2.0 * a * (z - _sc2p5d_gridz(z_min, dz, i_start) + dz) / dz + b
        deriv /= dz
        ele.z_deriv_grid[iz] = deriv / dz
    end
    return nothing
end

function _sc2p5d_bunch_analysis!(ele::SPACECHARGE2P5D, r_in::Matrix{Float64}, lost_flags::Vector{Int}, macro_size::Float64)
    x_min, x_max, y_min, y_max, z_min, z_max, n_active = _sc2p5d_extrema(r_in, lost_flags)
    if n_active < 2
        return nothing
    end
    if !(x_min < x_max && y_min < y_max && z_min < z_max)
        return nothing
    end

    x_min, x_max, y_min, y_max = _sc2p5d_adjust_bounds(x_min, x_max, y_min, y_max, Float64(ele.xy_ratio))
    fill!(ele.rho_grid, 0.0)
    fill!(ele.phi_grid, 0.0)
    fill!(ele.z_grid, 0.0)
    fill!(ele.z_deriv_grid, 0.0)

    @inbounds for c in 1:size(r_in, 1)
        if isone(lost_flags[c]) || isnan(r_in[c, 1])
            continue
        end
        _sc2p5d_bin_2d!(ele.rho_grid, r_in[c, 1], r_in[c, 3], macro_size, x_min, x_max, y_min, y_max)
        _sc2p5d_bin_1d!(ele.z_grid, r_in[c, 5], macro_size, z_min, z_max)
    end

    _sc2p5d_calculate_long_derivative!(ele, z_min, z_max)

    nx = ele.xsize
    ny = ele.ysize
    dx = (x_max - x_min) / (nx - 1)
    dy = (y_max - y_min) / (ny - 1)
    total_macrosize = 0.0
    x_avg = 0.0
    x2_avg = 0.0
    y_avg = 0.0
    y2_avg = 0.0
    @inbounds for ix in 1:nx
        x = _sc2p5d_gridx(x_min, dx, ix)
        for iy in 1:ny
            val = ele.rho_grid[ix, iy]
            # Mirror PyORBIT's getGridX(iy) usage here so the 2.5-D
            # longitudinal kick stays benchmark-compatible with SpaceChargeCalc2p5Drb.
            y = _sc2p5d_gridx(x_min, dx, iy)
            total_macrosize += val
            x_avg += x * val
            x2_avg += x * x * val
            y_avg += y * val
            y2_avg += y * y * val
        end
    end
    if total_macrosize <= 0.0
        return nothing
    end
    x_avg /= total_macrosize
    x2_avg /= total_macrosize
    y_avg /= total_macrosize
    y2_avg /= total_macrosize
    x2_avg = abs(x2_avg - x_avg * x_avg)
    y2_avg = abs(y2_avg - y_avg * y_avg)
    a_bunch = sqrt(max(0.0, 2.0 * (x2_avg + y2_avg)))
    return x_min, x_max, y_min, y_max, z_min, z_max, total_macrosize, x_avg, y_avg, a_bunch
end

function _sc2p5d_track!(ele::SPACECHARGE2P5D, r_in::Matrix{Float64}, particles::Beam{Float64})
    if ele.effective_len == 0.0 || size(r_in, 1) < 2 || particles.beta <= 0.0
        return nothing
    end
    if Float64(ele.pipe_radius) <= 0.0
        error("SPACECHARGE2P5D requires pipe_radius > 0.")
    end
    macro_size = _sc2p5d_macro_size(particles)
    if macro_size <= 0.0
        return nothing
    end
    analysis = _sc2p5d_bunch_analysis!(ele, r_in, particles.lost_flag, macro_size)
    analysis === nothing && return nothing
    x_min, x_max, y_min, y_max, z_min, z_max, total_macrosize, x_center, y_center, a_bunch = analysis
    _sc2p5d_find_potential!(ele, (x_max - x_min) / (ele.xsize - 1), (y_max - y_min) / (ele.ysize - 1))
    z_step = (z_max - z_min) / ele.zsize

    factor = 2.0 * ele.effective_len * particles.classrad0 / (particles.beta^2 * particles.gamma^3)
    factor /= (z_step * total_macrosize)

    long_sc_factor_in = a_bunch > 0.0 ? 1.0 + 2.0 * log(Float64(ele.pipe_radius) / a_bunch) : 0.0
    long_sc_factor_out = 2.0 * log(Float64(ele.pipe_radius))
    a_bunch_2 = a_bunch * a_bunch
    long_sc_factor = -ele.effective_len * particles.classrad0 * particles.mass / (particles.gamma^2)
    total_energy = particles.energy + particles.mass
    delta_scale = particles.beta^2 * total_energy
    if delta_scale == 0.0
        return nothing
    end

    @inbounds for c in 1:size(r_in, 1)
        if isone(particles.lost_flag[c]) || isnan(r_in[c, 1])
            continue
        end
        x = r_in[c, 1]
        y = r_in[c, 3]
        z = r_in[c, 5]
        ex, ey = _sc2p5d_grid2d_gradient(ele.phi_grid, x, y, x_min, x_max, y_min, y_max)
        ez = _sc2p5d_grid1d_value(ele.z_deriv_grid, z, z_min, z_max)
        lfactor = -_sc2p5d_grid1d_value(ele.z_grid, z, z_min, z_max) * factor

        r_in[c, 2] += ex * lfactor
        r_in[c, 4] += ey * lfactor

        r2 = (x - x_center)^2 + (y - y_center)^2
        if a_bunch_2 > 0.0
            if r2 <= a_bunch_2
                long_sc_coeff = long_sc_factor_in - r2 / a_bunch_2
            else
                long_sc_coeff = long_sc_factor_out - log(r2)
            end
            d_e = ez * long_sc_factor * long_sc_coeff
            r_in[c, 6] += d_e / delta_scale
        end

        r6 = @view r_in[c, :]
        if check_lost(r6)
            particles.lost_flag[c] = 1
        end
    end
    return nothing
end

function _sc2p5d_extrema(r_in::Matrix{DTPSAD{N, T}}, lost_flags::Vector{Int}) where {N, T <: Number}
    x_min = zero(DTPSAD{N, T})
    x_max = zero(DTPSAD{N, T})
    y_min = zero(DTPSAD{N, T})
    y_max = zero(DTPSAD{N, T})
    z_min = zero(DTPSAD{N, T})
    z_max = zero(DTPSAD{N, T})
    n_active = 0
    found_first = false
    @inbounds for c in 1:size(r_in, 1)
        if isone(lost_flags[c]) || isnan(r_in[c, 1])
            continue
        end
        x = r_in[c, 1]
        y = r_in[c, 3]
        z = r_in[c, 5]
        if !found_first
            x_min = x_max = x
            y_min = y_max = y
            z_min = z_max = z
            found_first = true
        else
            if x < x_min
                x_min = x
            end
            if x > x_max
                x_max = x
            end
            if y < y_min
                y_min = y
            end
            if y > y_max
                y_max = y
            end
            if z < z_min
                z_min = z
            end
            if z > z_max
                z_max = z
            end
        end
        n_active += 1
    end
    return x_min, x_max, y_min, y_max, z_min, z_max, n_active
end

function _sc2p5d_grid2d_ind_frac(coord::DTPSAD{N, T}, cmin::DTPSAD{N, T}, dc::DTPSAD{N, T}, n::Int) where {N, T <: Number}
    ind = trunc(Int, ((_sc2p5d_primal(coord) - _sc2p5d_primal(cmin)) / _sc2p5d_primal(dc)) + 0.5) + 1
    ind = clamp(ind, 2, n - 1)
    frac = (coord - (cmin + (ind - 1) * dc)) / dc
    return ind, frac
end

function _sc2p5d_bin_2d!(grid::Matrix{DTPSAD{N, T}}, x::DTPSAD{N, T}, y::DTPSAD{N, T}, weight::DTPSAD{N, T},
    x_min::DTPSAD{N, T}, x_max::DTPSAD{N, T}, y_min::DTPSAD{N, T}, y_max::DTPSAD{N, T}) where {N, T <: Number}
    if x < x_min || x > x_max || y < y_min || y > y_max
        return nothing
    end
    nx, ny = size(grid)
    dx = (x_max - x_min) / (nx - 1)
    dy = (y_max - y_min) / (ny - 1)
    ix, x_frac = _sc2p5d_grid2d_ind_frac(x, x_min, dx, nx)
    iy, y_frac = _sc2p5d_grid2d_ind_frac(y, y_min, dy, ny)
    x_frac2 = x_frac * x_frac
    y_frac2 = y_frac * y_frac
    wxm = 0.5 * (0.25 - x_frac + x_frac2)
    wx0 = 0.75 - x_frac2
    wxp = 0.5 * (0.25 + x_frac + x_frac2)
    wym = 0.5 * (0.25 - y_frac + y_frac2)
    wy0 = 0.75 - y_frac2
    wyp = 0.5 * (0.25 + y_frac + y_frac2)
    @inbounds begin
        grid[ix - 1, iy - 1] += wxm * wym * weight
        grid[ix - 1, iy    ] += wxm * wy0 * weight
        grid[ix - 1, iy + 1] += wxm * wyp * weight
        grid[ix    , iy - 1] += wx0 * wym * weight
        grid[ix    , iy    ] += wx0 * wy0 * weight
        grid[ix    , iy + 1] += wx0 * wyp * weight
        grid[ix + 1, iy - 1] += wxp * wym * weight
        grid[ix + 1, iy    ] += wxp * wy0 * weight
        grid[ix + 1, iy + 1] += wxp * wyp * weight
    end
    return nothing
end

function _sc2p5d_grid2d_gradient(grid::Matrix{DTPSAD{N, T}}, x::DTPSAD{N, T}, y::DTPSAD{N, T},
    x_min::DTPSAD{N, T}, x_max::DTPSAD{N, T}, y_min::DTPSAD{N, T}, y_max::DTPSAD{N, T}) where {N, T <: Number}
    nx, ny = size(grid)
    dx = (x_max - x_min) / (nx - 1)
    dy = (y_max - y_min) / (ny - 1)
    ix, x_frac = _sc2p5d_grid2d_ind_frac(x, x_min, dx, nx)
    iy, y_frac = _sc2p5d_grid2d_ind_frac(y, y_min, dy, ny)
    x_frac2 = x_frac * x_frac
    y_frac2 = y_frac * y_frac
    wxm = 0.5 * (0.25 - x_frac + x_frac2)
    wx0 = 0.75 - x_frac2
    wxp = 0.5 * (0.25 + x_frac + x_frac2)
    wym = 0.5 * (0.25 - y_frac + y_frac2)
    wy0 = 0.75 - y_frac2
    wyp = 0.5 * (0.25 + y_frac + y_frac2)
    dwxm = -0.5 + x_frac
    dwx0 = -2.0 * x_frac
    dwxp = 0.5 + x_frac
    dwym = -0.5 + y_frac
    dwy0 = -2.0 * y_frac
    dwyp = 0.5 + y_frac
    @inbounds begin
        ex = dwxm * wym * grid[ix - 1, iy - 1] +
             dwxm * wy0 * grid[ix - 1, iy    ] +
             dwxm * wyp * grid[ix - 1, iy + 1] +
             dwx0 * wym * grid[ix    , iy - 1] +
             dwx0 * wy0 * grid[ix    , iy    ] +
             dwx0 * wyp * grid[ix    , iy + 1] +
             dwxp * wym * grid[ix + 1, iy - 1] +
             dwxp * wy0 * grid[ix + 1, iy    ] +
             dwxp * wyp * grid[ix + 1, iy + 1]
        ey = wxm * dwym * grid[ix - 1, iy - 1] +
             wxm * dwy0 * grid[ix - 1, iy    ] +
             wxm * dwyp * grid[ix - 1, iy + 1] +
             wx0 * dwym * grid[ix    , iy - 1] +
             wx0 * dwy0 * grid[ix    , iy    ] +
             wx0 * dwyp * grid[ix    , iy + 1] +
             wxp * dwym * grid[ix + 1, iy - 1] +
             wxp * dwy0 * grid[ix + 1, iy    ] +
             wxp * dwyp * grid[ix + 1, iy + 1]
    end
    return ex / dx, ey / dy
end

function _sc2p5d_grid1d_get_ind_and_wz(z::DTPSAD{N, T}, z_min::DTPSAD{N, T}, dz::DTPSAD{N, T}, nz::Int) where {N, T <: Number}
    if nz <= 1
        return 1, 1, one(z_min), zero(z_min)
    end
    ind0 = trunc(Int, ((_sc2p5d_primal(z) - _sc2p5d_primal(z_min)) / _sc2p5d_primal(dz)) - 0.5)
    z_frac = (z - _sc2p5d_gridz(z_min, dz, ind0 + 1)) / dz
    if ind0 == 0 && z_frac < 0.0
        return nz, 1, -z_frac, 1.0 + z_frac
    end
    if ind0 == nz - 1
        if z_frac > 0.0
            return nz, 1, 1.0 - z_frac, z_frac
        end
        return nz - 1, nz, -z_frac, 1.0 + z_frac
    end
    return ind0 + 1, ind0 + 2, 1.0 - z_frac, z_frac
end

function _sc2p5d_bin_1d!(grid::Vector{DTPSAD{N, T}}, z::DTPSAD{N, T}, weight::DTPSAD{N, T},
    z_min::DTPSAD{N, T}, z_max::DTPSAD{N, T}) where {N, T <: Number}
    if z < z_min || z > z_max
        return nothing
    end
    nz = length(grid)
    if nz == 0
        return nothing
    end
    dz = (z_max - z_min) / nz
    iz0, izp, wz0, wzp = _sc2p5d_grid1d_get_ind_and_wz(z, z_min, dz, nz)
    if nz > 1
        @inbounds begin
            grid[iz0] += wz0 * weight
            grid[izp] += wzp * weight
        end
    else
        @inbounds grid[1] += weight
    end
    return nothing
end

function _sc2p5d_grid1d_value(grid::Vector{DTPSAD{N, T}}, z::DTPSAD{N, T}, z_min::DTPSAD{N, T}, z_max::DTPSAD{N, T}) where {N, T <: Number}
    if z < z_min || z > z_max
        return zero(DTPSAD{N, T})
    end
    nz = length(grid)
    if nz == 0
        return zero(DTPSAD{N, T})
    end
    if nz == 1
        return grid[1]
    end
    dz = (z_max - z_min) / nz
    iz0, izp, wz0, wzp = _sc2p5d_grid1d_get_ind_and_wz(z, z_min, dz, nz)
    @inbounds return wz0 * grid[iz0] + wzp * grid[izp]
end

function _sc2p5d_fft2(arr::Matrix{DTPSAD{N, T}}) where {N, T <: Number}
    vals = Array{ComplexF64}(undef, size(arr))
    derivs = ntuple(_ -> Array{ComplexF64}(undef, size(arr)), N)
    @inbounds for i in eachindex(arr)
        vals[i] = ComplexF64(arr[i].val)
        for k in 1:N
            derivs[k][i] = ComplexF64(arr[i].deriv[k])
        end
    end
    fft_vals = fft(vals)
    fft_derivs = map(fft, derivs)
    out = zeros(DTPSAD{N, ComplexF64}, size(arr)...)
    @inbounds for i in eachindex(out)
        out[i] = DTPSAD{N, ComplexF64}(fft_vals[i], SVector{N, ComplexF64}(ntuple(k -> fft_derivs[k][i], N)))
    end
    return out
end

function _sc2p5d_ifft2(arr::Matrix{DTPSAD{N, ComplexF64}}) where {N}
    vals = Array{ComplexF64}(undef, size(arr))
    derivs = ntuple(_ -> Array{ComplexF64}(undef, size(arr)), N)
    @inbounds for i in eachindex(arr)
        vals[i] = arr[i].val
        for k in 1:N
            derivs[k][i] = arr[i].deriv[k]
        end
    end
    ifft_vals = ifft(vals)
    ifft_derivs = map(ifft, derivs)
    out = zeros(DTPSAD{N, ComplexF64}, size(arr)...)
    @inbounds for i in eachindex(out)
        out[i] = DTPSAD{N, ComplexF64}(ifft_vals[i], SVector{N, ComplexF64}(ntuple(k -> ifft_derivs[k][i], N)))
    end
    return out
end

function _sc2p5d_find_potential_dtpsa(rho_grid::Matrix{DTPSAD{N, T}}, dx::DTPSAD{N, T}, dy::DTPSAD{N, T}) where {N, T <: Number}
    nx, ny = size(rho_grid)
    nx2 = 2 * nx
    ny2 = 2 * ny
    greens = zeros(DTPSAD{N, T}, nx2, ny2)
    half_x = div(nx2, 2)
    half_y = div(ny2, 2)
    @inbounds for iy in 0:half_y
        r_y = iy * dy
        for ix in 0:half_x
            r_x = ix * dx
            if ix == 0 && iy == 0
                greens[ix + 1, iy + 1] = zero(DTPSAD{N, T})
            else
                greens[ix + 1, iy + 1] = -0.5 * log(r_x * r_x + r_y * r_y)
            end
        end
        for ix in half_x + 1:nx2 - 1
            greens[ix + 1, iy + 1] = greens[nx2 - ix + 1, iy + 1]
        end
    end
    @inbounds for iy in half_y + 1:ny2 - 1
        for ix in 0:nx2 - 1
            greens[ix + 1, iy + 1] = greens[ix + 1, ny2 - iy + 1]
        end
    end
    rho_pad = zeros(DTPSAD{N, T}, nx2, ny2)
    @views rho_pad[1:nx, 1:ny] .= rho_grid
    phi_pad = _sc2p5d_ifft2(_sc2p5d_fft2(rho_pad) .* _sc2p5d_fft2(greens))
    phi_grid = zeros(DTPSAD{N, T}, nx, ny)
    @views phi_grid .= real.(phi_pad[1:nx, 1:ny])
    return phi_grid
end

function _sc2p5d_calculate_long_derivative!(z_grid::Vector{DTPSAD{N, T}}, z_deriv_grid::Vector{DTPSAD{N, T}},
    long_avg_n::Int, z_min::DTPSAD{N, T}, z_max::DTPSAD{N, T}) where {N, T <: Number}
    nz = length(z_grid)
    dz = (z_max - z_min) / nz
    fill!(z_deriv_grid, zero(DTPSAD{N, T}))
    if nz == 1
        return nothing
    elseif nz == 2
        z_mid = 0.5 * (_sc2p5d_gridz(z_min, dz, 1) + _sc2p5d_gridz(z_min, dz, 2))
        z_deriv = (_sc2p5d_grid1d_value(z_grid, z_mid + 0.5 * dz, z_min, z_max) -
                   _sc2p5d_grid1d_value(z_grid, z_mid - 0.5 * dz, z_min, z_max)) / dz
        z_deriv_grid[1] = z_deriv / dz
        z_deriv_grid[2] = z_deriv / dz
        return nothing
    end

    n_avg = min(long_avg_n, nz)
    s = zeros(DTPSAD{N, T}, 5, 2)
    @inbounds for iz in 1:nz
        z = _sc2p5d_gridz(z_min, dz, iz)
        i_start = iz - div(n_avg, 2)
        if i_start < 1
            i_start = 1
        end
        if (i_start + n_avg - 1) > nz
            i_start = nz - n_avg + 1
        end
        i_stop = i_start + n_avg - 1
        fill!(s, zero(DTPSAD{N, T}))
        for i in i_start:i_stop
            x = DTPSAD{N, T}(i - i_start + 1)
            y = z_grid[i]
            xpow = one(DTPSAD{N, T})
            for j in 0:4
                s[j + 1, 1] += xpow
                s[j + 1, 2] += xpow * y
                xpow *= x
            end
        end
        det = s[1, 1] * s[3, 1] * s[5, 1] -
              s[2, 1] * s[2, 1] * s[5, 1] -
              s[1, 1] * s[4, 1] * s[4, 1] +
              2.0 * s[2, 1] * s[3, 1] * s[4, 1] -
              s[3, 1] * s[3, 1] * s[3, 1]
        if abs(det) <= eps(Float64)
            z_deriv_grid[iz] = zero(DTPSAD{N, T})
            continue
        end
        a = (s[1, 2] * s[2, 1] * s[4, 1] - s[2, 2] * s[1, 1] * s[4, 1] - s[1, 2] * s[3, 1] * s[3, 1] +
             s[2, 2] * s[2, 1] * s[3, 1] + s[3, 2] * s[1, 1] * s[3, 1] - s[3, 2] * s[2, 1] * s[2, 1]) / det
        b = (s[2, 2] * s[1, 1] * s[5, 1] - s[1, 2] * s[2, 1] * s[5, 1] + s[1, 2] * s[3, 1] * s[4, 1] -
             s[3, 2] * s[1, 1] * s[4, 1] - s[2, 2] * s[3, 1] * s[3, 1] + s[3, 2] * s[2, 1] * s[3, 1]) / det
        deriv = 2.0 * a * (z - _sc2p5d_gridz(z_min, dz, i_start) + dz) / dz + b
        deriv /= dz
        z_deriv_grid[iz] = deriv / dz
    end
    return nothing
end

function _sc2p5d_bunch_analysis_dtpsa(ele::SPACECHARGE2P5D{DTPSAD{N, T}}, r_in::Matrix{DTPSAD{N, T}},
    lost_flags::Vector{Int}, macro_size::DTPSAD{N, T}) where {N, T <: Number}
    x_min, x_max, y_min, y_max, z_min, z_max, n_active = _sc2p5d_extrema(r_in, lost_flags)
    if n_active < 2
        return nothing
    end
    if !(x_min < x_max && y_min < y_max && z_min < z_max)
        return nothing
    end

    x_min, x_max, y_min, y_max = _sc2p5d_adjust_bounds(x_min, x_max, y_min, y_max, ele.xy_ratio)
    rho_grid = zeros(DTPSAD{N, T}, ele.xsize, ele.ysize)
    z_grid = zeros(DTPSAD{N, T}, ele.zsize)
    z_deriv_grid = zeros(DTPSAD{N, T}, ele.zsize)

    @inbounds for c in 1:size(r_in, 1)
        if isone(lost_flags[c]) || isnan(r_in[c, 1])
            continue
        end
        _sc2p5d_bin_2d!(rho_grid, r_in[c, 1], r_in[c, 3], macro_size, x_min, x_max, y_min, y_max)
        _sc2p5d_bin_1d!(z_grid, r_in[c, 5], macro_size, z_min, z_max)
    end

    _sc2p5d_calculate_long_derivative!(z_grid, z_deriv_grid, ele.long_avg_n, z_min, z_max)

    nx = ele.xsize
    ny = ele.ysize
    dx = (x_max - x_min) / (nx - 1)
    total_macrosize = zero(DTPSAD{N, T})
    x_avg = zero(DTPSAD{N, T})
    x2_avg = zero(DTPSAD{N, T})
    y_avg = zero(DTPSAD{N, T})
    y2_avg = zero(DTPSAD{N, T})
    @inbounds for ix in 1:nx
        x = _sc2p5d_gridx(x_min, dx, ix)
        for iy in 1:ny
            val = rho_grid[ix, iy]
            y = _sc2p5d_gridx(x_min, dx, iy)
            total_macrosize += val
            x_avg += x * val
            x2_avg += x * x * val
            y_avg += y * val
            y2_avg += y * y * val
        end
    end
    if total_macrosize <= 0.0
        return nothing
    end
    x_avg /= total_macrosize
    x2_avg /= total_macrosize
    y_avg /= total_macrosize
    y2_avg /= total_macrosize
    x2_avg = abs(x2_avg - x_avg * x_avg)
    y2_avg = abs(y2_avg - y_avg * y_avg)
    radial_second_moment = 2.0 * (x2_avg + y2_avg)
    a_bunch = radial_second_moment < 0.0 ? zero(DTPSAD{N, T}) : sqrt(radial_second_moment)
    return x_min, x_max, y_min, y_max, z_min, z_max, total_macrosize, x_avg, y_avg, a_bunch, rho_grid, z_grid, z_deriv_grid
end

function _sc2p5d_track!(ele::SPACECHARGE2P5D{DTPSAD{N, T}}, r_in::Matrix{DTPSAD{N, T}}, particles::Beam{DTPSAD{N, T}}) where {N, T <: Number}
    if _sc2p5d_primal(ele.effective_len) == 0.0 || size(r_in, 1) < 2 || particles.beta <= 0.0
        return nothing
    end
    if _sc2p5d_primal(ele.pipe_radius) <= 0.0
        error("SPACECHARGE2P5D requires pipe_radius > 0.")
    end
    macro_size = _sc2p5d_macro_size(particles)
    if macro_size <= 0.0
        return nothing
    end
    analysis = _sc2p5d_bunch_analysis_dtpsa(ele, r_in, particles.lost_flag, macro_size)
    analysis === nothing && return nothing
    x_min, x_max, y_min, y_max, z_min, z_max, total_macrosize, x_center, y_center, a_bunch, rho_grid, z_grid, z_deriv_grid = analysis
    phi_grid = _sc2p5d_find_potential_dtpsa(rho_grid, (x_max - x_min) / (ele.xsize - 1), (y_max - y_min) / (ele.ysize - 1))
    z_step = (z_max - z_min) / ele.zsize

    factor = 2.0 * ele.effective_len * particles.classrad0 / (particles.beta^2 * particles.gamma^3)
    factor /= (z_step * total_macrosize)

    long_sc_factor_in = a_bunch > 0.0 ? 1.0 + 2.0 * log(ele.pipe_radius / a_bunch) : zero(DTPSAD{N, T})
    long_sc_factor_out = 2.0 * log(ele.pipe_radius)
    a_bunch_2 = a_bunch * a_bunch
    long_sc_factor = -ele.effective_len * particles.classrad0 * particles.mass / (particles.gamma^2)
    total_energy = particles.energy + particles.mass
    delta_scale = particles.beta^2 * total_energy
    if delta_scale == 0.0
        return nothing
    end

    @inbounds for c in 1:size(r_in, 1)
        if isone(particles.lost_flag[c]) || isnan(r_in[c, 1])
            continue
        end
        x = r_in[c, 1]
        y = r_in[c, 3]
        z = r_in[c, 5]
        ex, ey = _sc2p5d_grid2d_gradient(phi_grid, x, y, x_min, x_max, y_min, y_max)
        ez = _sc2p5d_grid1d_value(z_deriv_grid, z, z_min, z_max)
        lfactor = -_sc2p5d_grid1d_value(z_grid, z, z_min, z_max) * factor

        r_in[c, 2] += ex * lfactor
        r_in[c, 4] += ey * lfactor

        if a_bunch_2 > 0.0
            r2 = (x - x_center)^2 + (y - y_center)^2
            if r2 <= a_bunch_2
                long_sc_coeff = long_sc_factor_in - r2 / a_bunch_2
            else
                long_sc_coeff = long_sc_factor_out - log(r2)
            end
            d_e = ez * long_sc_factor * long_sc_coeff
            r_in[c, 6] += d_e / delta_scale
        end

        r6 = @view r_in[c, :]
        if check_lost(r6)
            particles.lost_flag[c] = 1
        end
    end
    return nothing
end

function pass!(ele::SPACECHARGE2P5D{DTPSAD{N, T}}, r_in::Matrix{DTPSAD{N, T}}, num_particles::Int64,
    particles::Beam{DTPSAD{N, T}}) where {N, T <: Number}
    _sc2p5d_track!(ele, r_in, particles)
    return nothing
end

function pass!(ele::SPACECHARGE2P5D, r_in::Matrix{Float64}, num_particles::Int64, particles::Beam{Float64})
    _sc2p5d_track!(ele, r_in, particles)
    return nothing
end

function pass_P!(ele::SPACECHARGE2P5D, r_in::Matrix{Float64}, num_particles::Int64, particles::Beam{Float64})
    _sc2p5d_track!(ele, r_in, particles)
    return nothing
end

function pass_TPSA!(ele::SPACECHARGE2P5D, r_in::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}) where {T, TPS_Dim, Max_TPS_Degree}
    println("TPSA is not implemented for SPACECHARGE2P5D.")
    return nothing
end

function _insert_space_charge_2p5d_periodic(line::Vector{<:AbstractElement{Float64}}, sc_path_length_min::Float64; kwargs...)
    markers = Tuple{Int, Float64, Float64}[]
    length_total = 0.0
    running_path = 0.0
    for i in eachindex(line)
        if running_path > sc_path_length_min
            push!(markers, (i, length_total, running_path))
            running_path = 0.0
        end
        running_path += line[i].len
        length_total += line[i].len
    end
    rest_length = isempty(markers) ? length_total : (length_total - markers[end][2])
    pushfirst!(markers, (1, 0.0, rest_length))

    sc_specs = Tuple{Int, Float64}[]
    if length(markers) > 1
        for i in 1:length(markers) - 1
            push!(sc_specs, (markers[i][1], markers[i + 1][3]))
        end
    end
    push!(sc_specs, (markers[end][1], rest_length))

    new_line = AbstractElement{Float64}[]
    spec_index = 1
    for i in eachindex(line)
        while spec_index <= length(sc_specs) && sc_specs[spec_index][1] == i
            sc = SPACECHARGE2P5D(name="SC2P5D_$(spec_index)", len=0.0,
                effective_len=sc_specs[spec_index][2]; kwargs...)
            push!(new_line, sc)
            spec_index += 1
        end
        push!(new_line, line[i])
    end
    return new_line
end

function insert_space_charge_2p5d(line::Vector{<:AbstractElement{Float64}}, sc_path_length_min::Real;
    periodic::Bool = false, kwargs...)
    if isempty(line)
        return AbstractElement{Float64}[]
    end
    path_min = Float64(sc_path_length_min)
    if path_min <= 0.0
        error("insert_space_charge_2p5d requires sc_path_length_min > 0.")
    end
    if periodic
        return _insert_space_charge_2p5d_periodic(line, path_min; kwargs...)
    end
    new_line = AbstractElement{Float64}[]
    running_path = 0.0
    sc_index = 1
    for ele in line
        push!(new_line, ele)
        running_path += ele.len
        if running_path >= path_min
            sc = SPACECHARGE2P5D(name="SC2P5D_$(sc_index)", len=0.0, effective_len=running_path; kwargs...)
            push!(new_line, sc)
            running_path = 0.0
            sc_index += 1
        end
    end
    if running_path > 0.0
        sc = SPACECHARGE2P5D(name="SC2P5D_$(sc_index)", len=0.0, effective_len=running_path; kwargs...)
        push!(new_line, sc)
    end
    return new_line
end
