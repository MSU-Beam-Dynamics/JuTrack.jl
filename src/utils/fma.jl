using LinearAlgebra

function naff(nturn::Int, x::Vector{Float64}, px::Vector{Float64})
    qmin = 0.001
    qmax = 0.999
    afind = qmax
    bfind = qmin
    step_find = 0.002
    nux1 = 0.0
    nux2 = 0.0

    while step_find > 1.0e-12
        find = 0.0
        findnu = 0.0
        for xnu in afind:-step_find:bfind
            np1 = div(nturn, 2) + 1
            np = np1 - 1
            sumr = 0.0
            sumi = 0.0
            for ik in 0:np
                tw = 2.0 * ik / np - 1.0
                rr = 2.0 * pi * ik * xnu
                sumr += (x[ik + 1] * cos(rr) + px[ik + 1] * sin(rr)) * (1.0 + cos(pi * tw))
                sumi += (-x[ik + 1] * sin(rr) + px[ik + 1] * cos(rr)) * (1.0 + cos(pi * tw))
            end
            sumnorm = sqrt(sumr * sumr + sumi * sumi) / sqrt(np)
            if find < sumnorm
                find = sumnorm
                findnu = xnu
            end
        end
        afind = findnu + step_find
        bfind = findnu - step_find
        step_find /= 10.0
    end

    nux1 = findnu

    afind = qmax
    bfind = qmin
    step_find = 0.002
    findnu = 0.0
    while step_find > 1.0e-12
        find = 0.0
        findnu = 0.0
        for xnu in afind:-step_find:bfind
            np1 = div(nturn, 2) + 1
            np = np1 - 1
            sumr = 0.0
            sumi = 0.0
            for ik in 0:np
                ik1 = ik + div(nturn, 2)
                if ik1 >= length(x)
                    continue
                end
                tw = 2.0 * ik / np - 1.0
                rr = 2.0 * pi * ik * xnu
                sumr += (x[ik1 + 1] * cos(rr) + px[ik1 + 1] * sin(rr)) * (1.0 + cos(pi * tw))
                sumi += (-x[ik1 + 1] * sin(rr) + px[ik1 + 1] * cos(rr)) * (1.0 + cos(pi * tw))
            end
            sumnorm = sqrt(sumr * sumr + sumi * sumi) / sqrt(np)
            if find < sumnorm
                find = sumnorm
                findnu = xnu
            end
        end
        afind = findnu + step_find
        bfind = findnu - step_find
        step_find /= 10.0
    end

    nux2 = findnu

    return nux1, nux2
end

"""
    FMA(RING, beam, nturns)

do frequency map analysis (FMA) with NAFF

# Arguments
- RING: lattice
- beam: Beam object. Avoid zero initial coordinates.
- nturns: number of turns

# Return
- diff_nux: difference between nux1 and nux2
- nux1: horizontal tune 1
- nux2: horizontal tune 2
- nuy1: vertical tune 1
- nuy2: vertical tune 2
"""
function FMA(RING, beam, nturns)
    # do fma and return the tune
    # RING: lattice
    # beam: Beam object. Avoid zero initial coordinates.
    # nturns: number of turns
    twi = twissring(RING, 0.0, 1)
    rout = ringpass!(RING, beam, nturns, true)
    
    betax = twi[1].betax
    betay = twi[1].betay
    alphax = twi[1].alphax
    alphay = twi[1].alphay

    x = zeros(nturns, beam.np)
    px = zeros(nturns, beam.np)
    y = zeros(nturns, beam.np)
    py = zeros(nturns, beam.np)
    for i in 1:nturns
        x[i,:] = rout[i][:,1] ./ sqrt(betax)
        px[i,:] = -rout[i][:,2] .* sqrt(betax) .- alphax .* x[i]
        y[i,:] = rout[i][:,3] ./ sqrt(betay)
        py[i,:] = -rout[i][:,4] .* sqrt(betay) .- alphay .* y[i]
    end

    nux1 = zeros(beam.np)
    nux2 = zeros(beam.np)
    nuy1 = zeros(beam.np)
    nuy2 = zeros(beam.np)
    for i in 1:beam.np
        nux1[i], nux2[i] = naff(nturns, x[:,i], px[:,i])
        nuy1[i], nuy2[i] = naff(nturns, y[:,i], py[:,i])
    end

    diff_nux = log10.((nux2 .- nux1).^2 .+ (nuy2 .- nuy1).^2)
    return diff_nux, nux1, nux2, nuy1, nuy2
end
