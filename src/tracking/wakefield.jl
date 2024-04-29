

# note: 
# replace collect function
# round, ceil, floor, Int64(x), the derivatives are 0.0
# use AD-friendly interpolation function



function linear_interpolate(x, x_points, y_points)
    """
    Interpolates or extrapolates a value using linear interpolation.

    x: The point to interpolate or extrapolate.
    x_points: The x-coordinates of the data points.
    y_points: The y-coordinates of the data points.

    Returns the interpolated or extrapolated value at x.
    """
    if x <= x_points[1]
        slope = (y_points[2] - y_points[1]) / (x_points[2] - x_points[1])
        return y_points[1] + slope * (x - x_points[1])
    elseif x >= x_points[end]
        slope = (y_points[end] - y_points[end - 1]) / (x_points[end] - x_points[end - 1])
        return y_points[end] + slope * (x - x_points[end])
    else
        for i in 2:length(x_points)
            if x < x_points[i]
                slope = (y_points[i] - y_points[i - 1]) / (x_points[i] - x_points[i - 1])
                return y_points[i - 1] + slope * (x - x_points[i - 1])
            end
        end
    end
end

struct LongitudinalRLCWake#<: AbstractLongiWakefield
    freq::Float64
    Rshunt::Float64
    Q0::Float64
    wakefield::Function
end
function LongitudinalRLCWake(freq::Float64, Rshunt::Float64, Q0::Float64)
    Q0p=sqrt(Q0^2 - 1.0/4.0)
    ω0 = 2*pi*freq
    ω0p= ω0/Q0*Q0p
    wakefield = function (t::Float64)
        t>0 && return 0.0
        return Rshunt * ω0 /Q0 * (cos(ω0p * t) +  sin(ω0p * t) / 2 / Q0p) * exp(ω0 * t / 2 / Q0)
    end
    return LongitudinalRLCWake(freq, Rshunt, Q0, wakefield)
end

struct LongitudinalWake #<: AbstractLongiWakefield
    times::AbstractVector
    wakefields::AbstractVector
    wakefield::Function
end
function LongitudinalWake(times::AbstractVector, wakefields::AbstractVector, fliptime::Float64=-1.0)
    # wf = linear_interpolation(times, wakefields, extrapolation_bc=Line()) 
    wakefield_function = function (t::Float64)
        t>times[1]*fliptime && return 0.0
        return linear_interpolate(t*fliptime, times, wakefields)
    end
    return LongitudinalWake(times, wakefields, wakefield_function)
end




function LongiWakefieldPass!(r, num_macro, rlcwake, inzindex, eN_b2E, nbins, zhist, zhist_edges)
    zhist_center = zeros(nbins)
    zhist_center_end = (zhist_edges[nbins] + zhist_edges[nbins+1]) / 2.0
    wakefield = zeros(nbins)

    for i in 1:nbins
        zhist_center[i] = (zhist_edges[i] + zhist_edges[i+1]) / 2.0
        wakefield[i] = rlcwake.wakefield((zhist_center[i] - zhist_center_end) / 2.99792458e8)
    end
    
    wakepotential = zeros(nbins)
    
    halfzn = nbins ÷ 2
    @inbounds for i=1:nbins
        for j=-halfzn:halfzn
            if i-j>0 && i-j<=nbins
                wakepotential[i]+=wakefield[j+halfzn+1]*zhist[i-j]/num_macro
            end
        end
    end

    wakeatedge = zeros(nbins+1)
    wakeatedge[2:end-1] .= ((wakepotential[1:end-1]) .+ (wakepotential[2:end])) ./ 2.0
    wakeatedge[1] = 2*wakeatedge[2]-wakeatedge[3]
    wakeatedge[end] = 2*wakeatedge[end-1]-wakeatedge[end-2]

    zsep=(zhist_edges[2]-zhist_edges[1])
    @inbounds for i in 1:num_macro
        r6 = @view r[(i-1)*6+1:i*6]
        zloc=r6[5]
        zindex=inzindex[i]
        wake1=wakeatedge[zindex]
        wake2=wakeatedge[zindex+1]
        wakezloc=wake1+(wake2-wake1)*(zloc-zhist_edges[zindex])/zsep
        r6[6]-=wakezloc*eN_b2E
    end
    return nothing
end
function LongiWakefieldPass_P!(r, num_macro, rlcwake, inzindex, eN_b2E, nbins, zhist, zhist_edges)
    zhist_center = zeros(nbins)
    zhist_center_end = (zhist_edges[nbins] + zhist_edges[nbins+1]) / 2.0
    wakefield = zeros(nbins)

    for i in 1:nbins
        zhist_center[i] = (zhist_edges[i] + zhist_edges[i+1]) / 2.0
        wakefield[i] = rlcwake.wakefield((zhist_center[i] - zhist_center_end) / 2.99792458e8)
    end
    
    wakepotential = zeros(nbins)

    halfzn = nbins ÷ 2
    @inbounds @Threads.threads for i=1:nbins
        for j=-halfzn:halfzn
            if i-j>0 && i-j<=nbins
                wakepotential[i]+=wakefield[j+halfzn+1]*zhist[i-j]/num_macro
            end
        end
    end
    
    wakeatedge = zeros(nbins+1)
    wakeatedge[2:end-1] .= ((wakepotential[1:end-1]) .+ (wakepotential[2:end])) ./ 2.0
    wakeatedge[1] = 2*wakeatedge[2]-wakeatedge[3]
    wakeatedge[end] = 2*wakeatedge[end-1]-wakeatedge[end-2]

    zsep=(zhist_edges[2]-zhist_edges[1])
    @Threads.threads for i in 1:num_macro
        r6 = @view r[(i-1)*6+1:i*6]
        zloc=r6[5]
        zindex=inzindex[i]
        wake1=wakeatedge[zindex]
        wake2=wakeatedge[zindex+1]
        wakezloc=wake1+(wake2-wake1)*(zloc-zhist_edges[zindex])/zsep
        r6[6]-=wakezloc*eN_b2E
    end
    return nothing
end


function pass!(rlcwake::LongitudinalRLCWake, r, np, beam)
    histogram1DinZ!(beam, beam.znbin, beam.inzindex, beam.zhist, beam.zhist_edges)
    eN_b2E=beam.np*1.6021766208e-19*beam.charge^2/beam.energy/beam.beta/beam.beta/beam.atomnum
    LongiWakefieldPass!(r, np, rlcwake, beam.inzindex, eN_b2E, beam.znbin, beam.zhist, beam.zhist_edges)

end
function pass_P!(rlcwake::LongitudinalRLCWake, r, np, beam)
    histogram1DinZ!(beam, beam.znbin, beam.inzindex, beam.zhist, beam.zhist_edges)
    eN_b2E=beam.np*1.6021766208e-19*beam.charge^2/beam.energy/beam.beta/beam.beta/beam.atomnum
    LongiWakefieldPass_P!(r, np, rlcwake, beam.inzindex, eN_b2E, beam.znbin, beam.zhist, beam.zhist_edges)

end

function pass!(lm::LongitudinalRFMap, r, np, beam) 
    gamma = beam.gamma
    eta = lm.alphac - 1.0 / gamma / gamma
    for c in 1:np
        if beam.lost_flag[c] == 1
            continue
        end
        r6 = @view r[(c-1)*6+1:c*6]
        if !isnan(r6[1])
            r6[5] -= (2π * lm.RF.h * eta / lm.RF.k) * r6[6]
        end
    end
    return nothing
end
function pass_P!(lm::LongitudinalRFMap, r, np, beam) 
    gamma = beam.gamma
    eta = lm.alphac - 1.0 / gamma / gamma
    @Threads.threads for c in 1:np
        if beam.lost_flag[c] == 1
            continue
        end
        r6 = @view r[(c-1)*6+1:c*6]
        if !isnan(r6[1])
            r6[5] -= (2π * lm.RF.h * eta / lm.RF.k) * r6[6]
        end
    end
    return nothing
end
