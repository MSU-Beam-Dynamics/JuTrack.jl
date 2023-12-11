include("quadFringe.jl") 
include("../lattice/canonical_elements.jl")

# function fillPowerArray(x::Float64, order::Int)
#     xpow = Vector{Float64}(undef, order + 1)
#     xpow[1] = 1.0
#     for i in 2:order + 1
#         xpow[i] = xpow[i - 1] * x
#     end
#     return xpow
# end
function fillPowerArray(x::Float64, order::Int)
    # avoid muatating arrays for the implementation of Zygote.jl
    return [i == 1 ? 1.0 : x^(i-1) for i in 1:order + 1]
end


# function expansion_coefficients(n::Int)
#     # Calculate expansion coefficients with signs for (x+iy)^n/n!
#     expansion_coef = Vector{Float64}(undef, n + 1)
#     for i in 0:n
#         sign = isodd(div(i, 2)) ? -1.0 : 1.0
#         expansion_coef[i + 1] = sign / (factorial(i) * factorial(n - i))
#     end

#     return expansion_coef
# end
function expansion_coefficients(n::Int)
    # Calculate expansion coefficients with signs for (x+iy)^n/n!
    # avoid muatating arrays for the implementation of Zygote.jl
    return [(isodd(div(i, 2)) ? -1.0 : 1.0) / (factorial(i) * factorial(n - i)) for i in 0:n]
end



function apply_canonical_multipole_kicks(qx, qy, xpow, ypow, order, KnL, skew)
    sum_Fx = sum_Fy = 0.0
    coef = expansion_coefficients(order)  

    for i in 0:order
        if isodd(i)
            sum_Fx += coef[i + 1] * xpow[order - i + 1] * ypow[i + 1]
        else
            sum_Fy += coef[i + 1] * xpow[order - i + 1] * ypow[i + 1]
        end
    end

    if skew != 0
        sum_Fx, sum_Fy = -sum_Fy, sum_Fx
    end

    delta_qx = -KnL * sum_Fy
    delta_qy = KnL * sum_Fx

    return qx - KnL * sum_Fy, qy + KnL * sum_Fx, delta_qx, delta_qy
end

function convertSlopesToMomenta(xp, yp, delta)
    expandHamiltonian = 0
    if expandHamiltonian == 1
        qx = (1 + delta) * xp
        qy = (1 + delta) * yp
    else
        denom = sqrt(1 + xp^2 + yp^2)
        qx = (1 + delta) * xp / denom
        qy = (1 + delta) * yp / denom
    end
    return qx, qy
end

function convertMomentaToSlopes(qx, qy, delta)
    expandHamiltonian = 0
    if expandHamiltonian == 1
        xp = qx / (1 + delta)
        yp = qy / (1 + delta)
    else
        denom = (1 + delta)^2 - qx^2 - qy^2
        if denom < 0
            warn("Warning: particle acquired undefined slopes when integrating through kick multipole")
        end
        denom = sqrt(denom)
        xp = qx / denom
        yp = qy / denom
    end
    return xp, yp
end

# function offsetBeamCoordinates(coord, np, dx, dy, dz)
#     for ip in np:-1:1
#         part = coord[ip]
#         part[5] += dz * sqrt(1 + part[2]^2 + part[4]^2)
#         part[1] = part[1] - dx + dz * part[2]
#         part[3] = part[3] - dy + dz * part[4]
#         coord[ip] = part
#     end
#     return coord
# end
function offsetBeamCoordinates(coord, np, dx, dy, dz)
    return [
        [
            coord[ip][1] - dx + dz * coord[ip][2],
            coord[ip][2],
            coord[ip][3] - dy + dz * coord[ip][4],
            coord[ip][4],
            coord[ip][5] + dz * sqrt(1 + coord[ip][2]^2 + coord[ip][4]^2),
            coord[ip][6]
        ] 
        for ip in 1:np
    ]
end

function rotateBeamCoordinates(part, np, angle)
    if angle == 0 || abs(abs(angle) - 2 * π) < 1e-12
        return part
    end

    if abs(abs(angle) - π) < 1e-12
        cos_a, sin_a = -1, 0
    elseif abs(angle - π / 2) < 1e-12
        cos_a, sin_a = 0, 1
    elseif abs(angle + π / 2) < 1e-12
        cos_a, sin_a = 0, -1
    else
        cos_a, sin_a = cos(angle), sin(angle)
    end

    # for i in 1:np
    #     x, xp, y, yp = part[i][1], part[i][2], part[i][3], part[i][4]
    #     part[i][1] = x * cos_a + y * sin_a
    #     part[i][3] = -x * sin_a + y * cos_a
    #     part[i][2] = xp * cos_a + yp * sin_a
    #     part[i][4] = -xp * sin_a + yp * cos_a
    # end
    part_new = [
        [
            part[i][1] * cos_a + part[i][3] * sin_a,  # x
            part[i][2] * cos_a + part[i][4] * sin_a,  # xp
            -part[i][1] * sin_a + part[i][3] * cos_a, # y
            -part[i][2] * sin_a + part[i][4] * cos_a  # yp
        ]
        for i in 1:np
    ]
    return part_new
end


function integrate_kick_multipole_ordn(coord, dx, dy, xkick, ykick,
    Po, rad_coef, isr_coef,
    order, KnL, skew,
    n_parts, drift,
    integration_order,
    multData, steeringMultData,
    apData, dzLoss,
    radial, refTilt)

    BETA = 1.25992104989487316477
    COORD_LIMIT = 10
    SLOPE_LIMIT = 1

    driftFrac2 = [0.5, 0.5]
    kickFrac2 = [1.0, 0.0]

    driftFrac4 = [0.5 / (2 - BETA), (1 - BETA) / (2 - BETA) / 2, (1 - BETA) / (2 - BETA) / 2, 0.5 / (2 - BETA)]
    kickFrac4 = [1.0 / (2 - BETA), -BETA / (2 - BETA), 1.0 / (2 - BETA), 0.0]

    # From AOP-TN-2020-064
    driftFrac6 = [0.39225680523878, 0.5100434119184585, -0.47105338540975655, 0.0687531682525181,
                    0.0687531682525181, -0.47105338540975655, 0.5100434119184585, 0.39225680523878]
    kickFrac6 = [0.784513610477560, 0.235573213359357, -1.17767998417887, 1.3151863206839063,
                    -1.17767998417887,  0.235573213359357, 0.784513610477560, 0]
    nSubsteps = 0
    driftFrac = nothing
    kickFrac = nothing
    if integration_order == 2
        nSubsteps = 2
        driftFrac = driftFrac2
        kickFrac = kickFrac2
    elseif integration_order == 4
        nSubsteps = 4
        driftFrac = driftFrac4
        kickFrac = kickFrac4
    elseif integration_order == 6
        nSubsteps = 8
        driftFrac = driftFrac6
        kickFrac = kickFrac6
    else
        error("Error: invalid integration order: $integration_order")
    end
    drift = drift/n_parts
    xkick = xkick/n_parts
    ykick = ykick/n_parts

    x = coord[1]
    xp = coord[2]
    y = coord[3]
    yp = coord[4]
    s  = 0
    dp = coord[6]
    p = Po*(1+dp)
    beta0 = p/sqrt(p^2+1)

    if isnan(x) || isnan(xp) || isnan(y) || isnan(yp) || isnan(dp)
        return 0
    end

    if abs(x) > COORD_LIMIT || abs(y) > COORD_LIMIT || abs(xp) > SLOPE_LIMIT || abs(yp) > SLOPE_LIMIT
        return 0
    end

    qx, qy = convertSlopesToMomenta(xp, yp, dp)
    maxOrder = order
    xp, yp = convertMomentaToSlopes(qx, qy, dp)
    dzLoss = 0
    expandHamiltonian = 0
    for i_kick in 0:n_parts-1
        delta_qx = delta_qy = 0.0
        for step in 1:nSubsteps
            if drift != 0
                dsh = drift * driftFrac[step]
                x = x + xp * dsh
                y = y + yp * dsh
                if expandHamiltonian != 0
                    s = s + dsh * (1 + (xp^2 + yp^2) / 2)
                else
                    s = s + dsh *sqrt(1 + xp^2 + yp^2)
                end
                dzLoss = dzLoss + dsh
            end
            
            if kickFrac[step] == 0
                break
            end
            xpow = fillPowerArray(x, maxOrder)
            ypow = fillPowerArray(y, maxOrder)
            delta_qx = delta_qy = 0.0

            if radial == 0 ## no radial now
                for iOrder in 1:3
                    if KnL[iOrder] != 0
                        qx, qy, delta_qx, delta_qy = apply_canonical_multipole_kicks(qx, qy, xpow, ypow, 
                                                    order, KnL[iOrder]/n_parts*kickFrac[step], skew[iOrder])
                    end
                end
            # else
            #     applyRadialCanonicalMultipoleKicks
            end

            if xkick != 0
                qx, qy, delta_qx, delta_qy = apply_canonical_multipole_kicks(qx, qy, xpow, ypow, 
                                            0, -xkick*kickFrac[step], 0)
            end
            if ykick != 0
                qx, qy, delta_qx, delta_qy = apply_canonical_multipole_kicks(qx, qy, xpow, ypow, 
                                            0, ykick*kickFrac[step], 1)
            end

            if steeringMultData != nothing && steeringMultData.orders != 0
                for imult in 1:steeringMultData.orders
                    if steeringMultData.KnL[imult] != 0
                        qx, qy, delta_qx, delta_qy = apply_canonical_multipole_kicks(qx, qy, xpow, ypow, 
                                                    steeringMultData.order[imult], steeringMultData.KnL[imult]*xkick*kickFrac[step], 0)
                    end
                    if steeringMultData.JnL[imult] != 0
                        qx, qy, delta_qx, delta_qy = apply_canonical_multipole_kicks(qx, qy, xpow, ypow, 
                                                    steeringMultData.order[imult], steeringMultData.JnL[imult]*ykick*kickFrac[step], 1)
                    end
                end
            end

            if multData != nothing 
                # do kicks for spurious multipole
                for imult in 1:multData.orders
                    if !isnothing(multData.KnL) && multData.KnL[imult] != 0
                        qx, qy, delta_qx, delta_qy = apply_canonical_multipole_kicks(qx, qy, xpow, ypow, 
                                                    multData.order[imult], multData.KnL[imult]*kickFrac[step]/n_parts, 0)
                    end
                    if !isnothing(multData.JnL) && multData.JnL[imult] != 0
                        qx, qy, delta_qx, delta_qy = apply_canonical_multipole_kicks(qx, qy, xpow, ypow, 
                                                    multData.order[imult], multData.JnL[imult]*kickFrac[step]/n_parts, 1)
                    end
                end
            end

            xp, yp = convertMomentaToSlopes(qx, qy, dp)

            if (rad_coef != 0 || isr_coef != 0) && (drift !=0)
                qx = qx / (1 + dp)
                qy = qy / (1 + dp)
                deltaFactor = (1 + dp)^2

                delta_qx = delta_qx / kickFrac[step]
                delta_qy = delta_qy / kickFrac[step]
                F2 = (delta_qx/drift-xkick/drift)^2 + (delta_qy/drift+ykick/drift)^2

                delta_qx = 0.0
                delta_qy = 0.0
                dsFactor = sqrt(1 + xp^2 + yp^2)
                dsISRFactor = dsFactor*drift/(nSubsteps-1)
                dsFactor = dsFactor * kickFrac[step]
                if rad_coef != 0
                    dp = dp - rad_coef*deltaFactor*F2*dsFactor
                end
                if isr_coef != 0
                    srGaussianLimit = 3.0
                    dp = dp - isr_coef * deltaFactor * F2^0.75 *
                        sqrt(dsISRFactor)*gauss_rn_lim(0.0, 1.0, srGaussianLimit, random_2)
                end
                # if sigmaDelta2 != 0
                #     sigmaDelta2 = sigmaDelta2 + (isr_coef * deltaFactor)^2 * F2^1.5 * dsISRFactor
                # end
                qx = qx * (1 + dp)
                qy = qy * (1 + dp)
                xp, yp = convertMomentaToSlopes(qx, qy, dp)
                
            end
            
        end

    end

    # if i_part < 0 || i_part == n_parts - 1
    #     if edgeMultData != nothing && edgeMultData.orders != 0
    #         xpow = fillPowerArray(x, maxOrder)
    #         ypow = fillPowerArray(y, maxOrder)
    #         for imult in 1:edgeMultData.orders
    #             qx, qy, delta_qx, delta_qy = apply_canonical_multipole_kicks(qx, qy, xpow, ypow, 
    #                                         edgeMultData.order[imult], edgeMultData.KnL[imult], 0)
    #             qx, qy, delta_qx, delta_qy = apply_canonical_multipole_kicks(qx, qy, xpow, ypow, 
    #                                         edgeMultData.order[imult], edgeMultData.JnL[imult], 1)
    #         end
    #     end
    # end

    xp, yp = convertMomentaToSlopes(qx, qy, dp)

    if rad_coef != 0
        p = Po * (1 + dp)
        beta1 = p / sqrt(p^2 + 1)
        coord5 = beta1 * (coord[5]/beta0 + 2*s/(beta0 + beta1))
    else
        coord5 = coord[5] + s
    end
    coord_new = [x, xp, y, yp, coord5, dp]
    return coord_new
end


function multipole_tracking2(particle, n_part, elem, Po)
    # constants
    particleMass = 9.1093897e-31 # Electron mass (kg)
    particleCharge = 1.60217733e-19 # Electron charge (C)
    c_mks = 299792458.0 # Speed of light (m/s)
    epsilon_o = 8.854187817e-12 # Permittivity of vacuum (F/m)
    me_mev = 0.51099906 # Electron mass (MeV)
    particleRadius = particleCharge^2 / (4 * pi * epsilon_o * particleMass * c_mks^2) # Classical electron radius (m)
    
    skew = zeros(Int16, 3)
    dx = dy = dz = 0.0 
    nSlices = integ_order = 0
    i_part = i_top = 0
    coord = []
    drift = 0.0
    tilt = pitch = yaw = rad_coef = isr_coef = xkick = ykick = dzLoss = 0.0
    sextWarning = quadWarning = octWarning = false
    multData = steeringMultData = apData = nothing

    if elem isa KQUAD
        nSlices = elem.nSlices
        order = 1
        if elem.bore != 0
            KnL1 = elem.B / elem.bore * (particleCharge / (particleMass * c_mks * Po)) * elem.len * (1 + elem.fse)
            KnL = [KnL1, 0, 0]
        else
            KnL1 = elem.k1 * elem.len * (1 + elem.fse)
            KnL = [KnL1, 0, 0]
        end
    
        drift = elem.len
        tilt = elem.tilt
        # pitch = elem.pitch
        # yaw = elem.yaw
        dx = elem.dx
        dy = elem.dy
        dz = elem.dz
        xkick = elem.xkick
        ykick = elem.ykick
        integ_order = elem.integration_order

        if elem.synch_rad != 0
            rad_coef = particleCharge^2 * Po^3 / (6 * pi * epsilon_o * sqrt(c_mks) * particleMass)
        end
        # isr is not implemented yet
        # isr_coef = particleRadius *sqrt(55.0 / (24 * sqrt(3)) * Po^5 * 137.0359895)
        if elem.isr == 0 || (elem.isr1Particle == 0 && n_part == 1 )
            # Minus sign indicates we accumulate into sigmaDelta^2 only, don't perturb particles
            isr_coef *= -1
        end
        if elem.len < 1e-6 && (elem.isr != 0 || elem.synch_rad != 0)
            rad_coef = 0.0
            isr_coef = 0.0 # avoid unphysical results
            if quadWarning == false
                println("Warning: one or more quadrupoles with length < 1e-6 have had SYNCH_RAD=0 and ISR=0 forced to avoid unphysical results")
                quadWarning = true
            end 
        end

    elseif elem isa KSEXT
        nSlices = elem.nSlices
        order = 2
        if elem.bore != 0
            KnL2 = 2 * elem.B / elem.bore^2 * (particleCharge / (particleMass * c_mks * Po)) * elem.len * (1 + elem.fse)
            KnL = [0, KnL2, 0]
        else
            KnL2 = elem.k2 * elem.len * (1 + elem.fse)
            KnL = [0, KnL2, 0]
        end
        drift = elem.len
        tilt = elem.tilt
        # pitch = elem.pitch
        # yaw = elem.yaw
        dx = elem.dx
        dy = elem.dy
        dz = elem.dz
        integ_order = elem.integration_order

        if elem.synch_rad != 0
            rad_coef = particleCharge^2 * Po^3 / (6 * pi * epsilon_o * sqrt(c_mks) * particleMass)
        end

        # isr is not implemented
        # isr_coef = particleRadius *sqrt(55.0 / (24 * sqrt(3)) * Po^5 * 137.0359895)
        if elem.isr == 0 || (elem.isr1Particle == 0 && n_part == 1 )
            # Minus sign indicates we accumulate into sigmaDelta^2 only, don't perturb particles
            isr_coef *= -1
        end
        if elem.len < 1e-6 && (elem.isr != 0 || elem.synch_rad != 0)
            rad_coef = 0.0
            isr_coef = 0.0 # avoid unphysical results
            if sextWarning == false
                println("Warning: one or more sextupoles with length < 1e-6 have had SYNCH_RAD=0 and ISR=0 forced to avoid unphysical results")
                sextWarning = true
            end 
        end

    elseif elem isa KOCT
        nSlices = elem.nSlices
        order = 3
        if elem.bore != 0
            KnL3 = 6 * elem.B / elem.bore^3 * (particleCharge / (particleMass * c_mks * Po)) * elem.len * (1 + elem.fse)
            KnL = [0, 0, KnL3]
        else
            KnL3 = elem.k3 * elem.len * (1 + elem.fse)
            KnL = [0, 0, KnL3]
        end

        drift = elem.len
        tilt = elem.tilt
        # pitch = elem.pitch
        # yaw = elem.yaw
        dx = elem.dx
        dy = elem.dy
        dz = elem.dz
        integ_order = elem.integration_order

        if elem.synch_rad != 0
            rad_coef = particleCharge^2 * Po^3 / (6 * pi * epsilon_o * sqrt(c_mks) * particleMass)
        end

        # isr is not implemented
        # isr_coef = particleRadius *sqrt(55.0 / (24 * sqrt(3)) * Po^5 * 137.0359895)
        if elem.isr == 0 || (elem.isr1Particle == 0 && n_part == 1 )
            # Minus sign indicates we accumulate into sigmaDelta^2 only, don't perturb particles
            isr_coef *= -1
        end
        if elem.len < 1e-6 && (elem.isr != 0 || elem.synch_rad != 0)
            rad_coef = 0.0
            isr_coef = 0.0 # avoid unphysical results
            if octWarning == false
                println("Warning: one or more octupoles with length < 1e-6 have had SYNCH_RAD=0 and ISR=0 forced to avoid unphysical results")
                octWarning = true
            end 
        end

    end

    if order < 0
        error("Error: invalid order: $order")
    end
    if integ_order != 2 && integ_order != 4 && integ_order != 6
        error("Error: invalid integration order: $integ_order")
    end

    i_top = n_part -1


    if dx != 0 || dy != 0 || dz != 0 
        particle_new1 = offsetBeamCoordinates(particle, n_part, dx, dy, dz)
    else
        particle_new1 = particle
    end
    if tilt != 0
        particle_new2 = rotateBeamCoordinates(particle_new1, n_part, tilt)
    else
        particle_new2 = particle_new1
    end
    


    if elem isa KQUAD
        if elem.edge1_effects !=0
            particle_new3 = quadFringe(particle_new2, n_part, elem.k1, elem.fringeIntM, elem.fringeIntP, 
                                    false, -1, elem.edge1_effects, elem.edge1Linear, elem.edge1NonlinearFactor)
        else
            particle_new3 = particle_new2
        end
    else
        particle_new3 = particle_new2
    end

    # for i_part in 1:i_top+1
    #     coord = particle[i_part]
    #     if elem isa KQUAD
    #         coord = integrate_kick_multipole_ordn(coord, dx, dy, xkick, ykick, Po, rad_coef, 
    #                                             isr_coef, order, KnL, skew, nSlices, 
    #                                             drift, integ_order, multData, 
    #                                             steeringMultData, apData, dzLoss, elem.radial, tilt)
    #     else
    #         coord = integrate_kick_multipole_ordn(coord, dx, dy, xkick, ykick, Po, rad_coef, 
    #                                             isr_coef, order, KnL, skew, nSlices, 
    #                                             drift, integ_order, multData, 
    #                                             steeringMultData, apData, dzLoss, 0, tilt)
    #     end
    # end
    particles_new4 = [
    if elem isa KQUAD
        integrate_kick_multipole_ordn(particle_new3[i_part], dx, dy, xkick, ykick, Po, rad_coef, 
                                      isr_coef, order, KnL, skew, nSlices, 
                                      drift, integ_order, multData, 
                                      steeringMultData, apData, dzLoss, elem.radial, tilt)
    else
        integrate_kick_multipole_ordn(particle_new3[i_part], dx, dy, xkick, ykick, Po, rad_coef, 
                                      isr_coef, order, KnL, skew, nSlices, 
                                      drift, integ_order, multData, 
                                      steeringMultData, apData, dzLoss, 0, tilt)
    end
    for i_part in 1:i_top+1]


    if elem isa KQUAD
        if elem.edge2_effects !=0
            particles_new5 = quadFringe(particles_new4, n_part, elem.k1, elem.fringeIntM, elem.fringeIntP, 
                                    false, 1, elem.edge2_effects, elem.edge2Linear, elem.edge2NonlinearFactor)
        else
            particles_new5 = particles_new4
        end
    else
        particles_new5 = particles_new4
    end
    expandHamiltonian = 0
    return particles_new5
end

# double precision
# particle = [Float64[0.001, 0.0001, 0.0005, 0.0002, 0.0, 0.0], Float64[0.001, 0.0, 0.0, 0.0, 1.0, 0.0]]
# Quad = KQUAD("Q",1.0, 1.0, 0.0, 
#                 0.0, 0.0, 0.0, 0.0, 0.0, 1, 
#                 1, 1, 1.0, 1, 1.0, [1.0, 0.1, 0.0, 0.0, 0.0], [1.0, 0.1, 0.0, 0.0, 0.0], 0, 4, 4, 0.0, 0.0, 0, 0)

# Sext = KSEXT("S",1.0, 1.0, 0.0, 
#                 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4, 4, 0.0, 0.0, 0.0, 0.0)

# Oct = KOCT("O",1.0, 1.0, 0.0, 
#                 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4, 4, 0.0, 0.0, 0.0, 0.0)

# n_part = 2
# rout = multipole_tracking2(particle, n_part, Sext, 1.0)
# println(rout)

# function f(k)
#     Sext = KSEXT(k, 1.0, 0.0, 
#                 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4, 4, 0.0, 0.0, 0.0, 0.0)
#     Oct = KOCT(k, 1.0, 0.0, 
#                 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4, 4, 0.0, 0.0, 0.0, 0.0)
#     Quad = KQUAD(k, 1.0, 0.0, 
#                 0.0, 0.0, 0.0, 0.0, 0.0, 1, 
#                 1, 1, 1.0, 1, 1.0, [1.0, 0.1, 0.0, 0.0, 0.0], [1.0, 0.1, 0.0, 0.0, 0.0], 0, 4, 4, 0.0, 0.0, 0, 0)
#     particle = [Float64[0.001, 0.0001, 0.0005, 0.0002, 0.0, 0.0], Float64[0.001, 0.0, 0.0, 0.0, 0.0, 0.0]]
#     rout = multipole_tracking2(particle, 2, Quad, 1.0)
#     return rout[1]
# end
# using Zygote
# println(f(1.0))
# grad = jacobian(f, 1.0)
# println(grad)