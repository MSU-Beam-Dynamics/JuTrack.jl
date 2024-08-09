

# note: 
# replace collect function
# round, ceil, floor, Int64(x), the derivatives are 0.0
# use AD-friendly interpolation function

function LongiWakefieldPass!(r, num_macro, rlcwake, inzindex, eN_b2E, nbins, zhist, zhist_edges)
    zhist_center = zeros(nbins)
    zhist_center .= ((zhist_edges[1:end-1]) .+ (zhist_edges[2:end]))./2.0
    wakefield = zeros(nbins)
    for i in 1:nbins
        t = (zhist_center[i] -  0.0*zhist_center[end]) / 2.99792458e8
        wakefield[i] = wakefieldfunc_RLCWake(rlcwake, t)
    end
    wakepotential = zeros(nbins)
    
    halfzn = nbins ÷ 2
    for i=1:nbins
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

    zsep = zhist_edges[2]-zhist_edges[1]
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
    zhist_center .= ((zhist_edges[1:end-1]) .+ (zhist_edges[2:end]))./2.0
    wakefield = zeros(nbins)
    for i in 1:nbins
        t = (zhist_center[i] -  0.0*zhist_center[end]) / 2.99792458e8
        wakefield[i] = wakefieldfunc_RLCWake(rlcwake, t)
    end
    wakepotential = zeros(nbins)
    
    halfzn = nbins ÷ 2
    for i=1:nbins
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

    zsep = zhist_edges[2]-zhist_edges[1]
    @inbounds Threads.@threads for i in 1:num_macro
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
