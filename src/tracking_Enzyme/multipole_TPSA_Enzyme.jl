include("quadFringe_TPSA_Enzyme.jl") 
include("../lattice/canonical_elements.jl")
include("../TPSA_Enzyme/TPSA_fixedmap.jl")

function expansion_coefficients(n::Int)
    # is_even = 1  # Start with even (0 is even)
    # coefficients = zeros(Float64, n + 1)

    # for i in 0:n
    #     if is_even == 1
    #         coefficient = 1.0 / (doublefactorial(i) * doublefactorial(n - i))
    #     else
    #         coefficient = -1.0 / (doublefactorial(i) * doublefactorial(n - i))
    #     end
    #     coefficients[i + 1] = coefficient
    #     is_even = -is_even  
    # end

    # return coefficients
    if n == 1
        return [1.0, -1.0]
    elseif n == 2
        return [0.5, -1.0 , 0.5]
    elseif n == 3
        return [1.0/3.0, -0.5, 0.5, 1.0/3.0]
    elseif n == 4
        return [0.125, -1.0/3.0, 0.25, -1.0/3.0, 0.125]
    end
end



function apply_canonical_multipole_kicks(qx::CTPS{T, TPS_Dim, Max_TPS_Degree}, qy::CTPS{T, TPS_Dim, Max_TPS_Degree}, 
                x::CTPS{T, TPS_Dim, Max_TPS_Degree}, y::CTPS{T, TPS_Dim, Max_TPS_Degree}, 
                order::Int, KnL::Float64, skew::Int) where {T, TPS_Dim, Max_TPS_Degree}
    # different from multi-particle tracking
    sum_Fx = CTPS(0.0, TPS_Dim, Max_TPS_Degree)
    sum_Fy = CTPS(0.0, TPS_Dim, Max_TPS_Degree)
    coef = expansion_coefficients(order)  
    
    for i in 0:order
        if isodd(i)
            sum_Fx += coef[i + 1] * x^(order - i) * y^i #xpow[order - i + 1] * ypow[i + 1]
        else
            sum_Fy += coef[i + 1] * x^(order - i) * y^i #xpow[order - i + 1] * ypow[i + 1]
        end
    end


    # sum_Fy = coef[1] * x^order
    # sum_Fx = coef[2] * x^(order - 1) * y

    # odd_para = -1
    # for i in 2:order
    #     if odd_para == 1
    #         sum_Fx += coef[i + 1] * x^(order - i) * y^i
    #     else
    #         sum_Fy += coef[i + 1] * x^(order - i) * y^i
    #     end
    #     odd_para = -odd_para
    # end

    if skew != 0
        sum_Fx, sum_Fy = -sum_Fy, sum_Fx
    end

    delta_qx = -KnL * sum_Fy
    delta_qy = KnL * sum_Fx

    qx -= KnL * sum_Fy
    qy += KnL * sum_Fx
    return qx, qy, delta_qx, delta_qy
end

function convertSlopesToMomenta(xp::CTPS, yp::CTPS, delta::CTPS)
    expandHamiltonian = 0
    # if expandHamiltonian == 1
    #     qx = (1.0 + delta) * xp
    #     qy = (1.0 + delta) * yp
    # else
        denom = sqrt(1.0 + xp^2 + yp^2)
        qx = (1.0 + delta) * xp / denom
        qy = (1.0 + delta) * yp / denom
    # end
    return qx, qy
end

function convertMomentaToSlopes(qx::CTPS, qy::CTPS, delta::CTPS)
    expandHamiltonian = 0
    # if expandHamiltonian == 1
    #     xp = qx / (1.0 + delta)
    #     yp = qy / (1.0 + delta)
    # else
        denom = (1.0 + delta)^2 - qx^2 - qy^2
        # if denom < 0
        #     warn("Warning: particle acquired undefined slopes when integrating through kick multipole")
        # end
        # denom = sqrt(denom)
        xp = qx / sqrt(denom)
        yp = qy / sqrt(denom)
    # end
    return xp, yp
end


function offsetBeamCoordinates(x, xp, y, yp, z, delta, np, dx, dy, dz)
    return x - dx + dz * xp, xp, y - dy + dz * yp, yp, z + dz * sqrt(1.0 + xp^2 + yp^2), delta
end

function rotateBeamCoordinates(x, xp, y, yp, z, delta, np, angle)
    if angle == 0 || abs(abs(angle) - 2.0 * π) < 1e-12
        return x, xp, y, yp, z, delta
    end

    if abs(abs(angle) - π) < 1e-12
        cos_a, sin_a = -1.0, 0.0
    elseif abs(angle - π / 2) < 1e-12
        cos_a, sin_a = 0.0, 1.0
    elseif abs(angle + π / 2) < 1e-12
        cos_a, sin_a = 0.0, -1.0
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
    x_new = x * cos_a + y * sin_a
    xp_new = xp * cos_a + yp * sin_a
    y_new = -x * sin_a + y * cos_a
    yp_new = -xp * sin_a + yp * cos_a
    return x_new, xp_new, y_new, yp_new, z, delta
end


function integrate_kick_multipole_ordn(x, xp, y, yp, z, delta, dx, dy, xkick, ykick,
    Po, rad_coef, isr_coef, synch_rad, 
    order, KnL1, KnL2, KnL3, skew,
    n_parts, drift,
    integration_order,
    multData, steeringMultData,
    apData, dzLoss,
    radial, refTilt)

    BETA = 1.25992104989487316477
    # COORD_LIMIT = 10
    # SLOPE_LIMIT = 1

    driftFrac2 = [0.5, 0.5]
    kickFrac2 = [1.0, 0.0]

    driftFrac4 = [0.5 / (2.0 - BETA), (1.0 - BETA) / (2.0 - BETA) / 2.0, (1.0 - BETA) / (2.0 - BETA) / 2.0, 0.5 / (2.0 - BETA)]
    kickFrac4 = [1.0 / (2.0 - BETA), -BETA / (2.0 - BETA), 1.0 / (2.0 - BETA), 0.0]
 
    # From AOP-TN-2020-064
    driftFrac6 = [0.39225680523878, 0.5100434119184585, -0.47105338540975655, 0.0687531682525181,
                    0.0687531682525181, -0.47105338540975655, 0.5100434119184585, 0.39225680523878]
    kickFrac6 = [0.784513610477560, 0.235573213359357, -1.17767998417887, 1.3151863206839063,
                    -1.17767998417887,  0.235573213359357, 0.784513610477560, 0.0]

    if integration_order == 2
        nSubsteps = 2
        nSubsteps_float = 2.0
        driftFrac = driftFrac2
        kickFrac = kickFrac2
    elseif integration_order == 4
        nSubsteps = 4
        nSubsteps_float = 4.0
        driftFrac = driftFrac4
        kickFrac = kickFrac4
    elseif integration_order == 6
        nSubsteps = 8
        nSubsteps_float = 8.0
        driftFrac = driftFrac6
        kickFrac = kickFrac6
    else
        error("Error: invalid integration order: $integration_order")
    end
    drift /= n_parts
    xkick /= n_parts
    ykick /= n_parts


    s  = CTPS(z)
    s_map = s.map
    for s_map_i in 1:length(s_map)
        s_map[s_map_i] = 0.0
    end
    dp = CTPS(delta)
    p = (1.0+dp)
    p *= Po
    temp = p^2 + 1.0
    beta0 = p/sqrt(temp)

    qx, qy = convertSlopesToMomenta(xp, yp, dp)
    # xp, yp = convertMomentaToSlopes(qx, qy, dp)

    for i_kick in 0:n_parts-1
        for step in 1:nSubsteps
            if drift != 0.0
                dsh = drift * driftFrac[step]
                x += xp * dsh
                y += yp * dsh
                s += dsh *sqrt(1.0 + xp^2 + yp^2)
                dzLoss = dzLoss + dsh
            end
            
            if kickFrac[step] == 0.0
                break
            end

            # if radial == 0 ## no radial now
            #     for iOrder in 1:3
            #         if iOrder == 1 #&& KnL1 != 0.0
            #             qx, qy, delta_qx, delta_qy = apply_canonical_multipole_kicks(qx, qy, x, y, 
            #                                         order, KnL1/n_parts*kickFrac[step], skew[iOrder])
            #         elseif iOrder == 2 #&& KnL2 != 0.0
            #             qx, qy, delta_qx, delta_qy = apply_canonical_multipole_kicks(qx, qy, x, y, 
            #                                         order, KnL2/n_parts*kickFrac[step], skew[iOrder])
            #         elseif iOrder == 3 #&& KnL3 != 0.0
            #             qx, qy, delta_qx, delta_qy = apply_canonical_multipole_kicks(qx, qy, x, y, 
            #                                         order, KnL3/n_parts*kickFrac[step], skew[iOrder])
            #         end
            #     end
            # # else
            # #     applyRadialCanonicalMultipoleKicks
            # end

            if order == 1 #&& KnL1 != 0.0
                qx, qy, delta_qx, delta_qy = apply_canonical_multipole_kicks(qx, qy, x, y, 
                                                    order, KnL1/n_parts*kickFrac[step], skew[order])
            elseif order == 2 #&& KnL2 != 0.0
                qx, qy, delta_qx, delta_qy = apply_canonical_multipole_kicks(qx, qy, x, y, 
                                                    order, KnL2/n_parts*kickFrac[step], skew[order])
            elseif order == 3 #&& KnL3 != 0.0
                qx, qy, delta_qx, delta_qy = apply_canonical_multipole_kicks(qx, qy, x, y, 
                                                    order, KnL3/n_parts*kickFrac[step], skew[order])
            end

            if xkick != 0.0
                qx, qy, delta_qx, delta_qy = apply_canonical_multipole_kicks(qx, qy, x, y, 
                                            0, -xkick*kickFrac[step], 0)
            end
            if ykick != 0.0
                qx, qy, delta_qx, delta_qy = apply_canonical_multipole_kicks(qx, qy, x, y, 
                                            0, ykick*kickFrac[step], 1)
            end

            # xp, yp = convertMomentaToSlopes(qx, qy, dp)
            denom = (1.0 + dp)^2 - qx^2 - qy^2
            xp = qx / sqrt(denom)
            yp = qy / sqrt(denom)

            # if (rad_coef != 0.0 || isr_coef != 0.0) && (drift !=0.0)
            if synch_rad == 1  && drift !=0.0
                qx = qx / (1.0 + dp) 
                qy = qy / (1.0 + dp)
                deltaFactor = (1.0 + dp)^2

                delta_qx /= kickFrac[step]
                delta_qy /= kickFrac[step]
                F2 = (delta_qx/drift-xkick/drift)^2 + (delta_qy/drift+ykick/drift)^2

                dsFactor = sqrt(1.0 + xp^2 + yp^2) * drift * kickFrac[step]
                dsISRFactor = dsFactor*drift/(nSubsteps_float-1.0)
                dsFactor *= drift * kickFrac[step]

                # if rad_coef != 0.0
                    dp -= rad_coef*deltaFactor*F2*dsFactor
                # end

                # not implemented yet
                # if isr_coef > 0
                #     srGaussianLimit = 3.0
                #     dp -= isr_coef * deltaFactor * F2^0.75 *
                #         sqrt(dsISRFactor)*gauss_rn_lim(0.0, 1.0, srGaussianLimit, random_2)
                # end
                # if sigmaDelta2 != 0
                #     sigmaDelta2 = sigmaDelta2 + (isr_coef * deltaFactor)^2 * F2^1.5 * dsISRFactor
                # end

                qx = qx * (1.0 + dp)
                qy = qy * (1.0 + dp)

                denom1 = (1.0 + dp)^2 - qx^2 - qy^2
                xp = qx / sqrt(denom1)
                yp = qy / sqrt(denom1)
                # xp, yp = convertMomentaToSlopes(qx, qy, dp)
                
            end
            
        end

    end

    denom1 = (1.0 + dp)^2 - qx^2 - qy^2
    xp = qx / sqrt(denom1)
    yp = qy / sqrt(denom1)
    # xp, yp = convertMomentaToSlopes(qx, qy, dp)

    if synch_rad == 1
        p = Po * (1.0 + dp)
        beta1 = p / sqrt(p^2 + 1.0)
        coord5 = beta1 * (z/beta0 + 2.0*s/(beta0 + beta1))
    else
        coord5 = z + s
    end
    return x, xp, y, yp, coord5, dp
end


function pass_TPSA(x::CTPS, xp::CTPS, y::CTPS, yp::CTPS, z::CTPS, delta::CTPS, elem::KQUAD, Po::Float64, sigmaDelta2::Float64)
    # constants
    particleMass = 9.1093897e-31 # Electron mass (kg)
    particleCharge = 1.60217733e-19 # Electron charge (C)
    c_mks = 299792458.0 # Speed of light (m/s)
    epsilon_o = 8.854187817e-12 # Permittivity of vacuum (F/m)
    me_mev = 0.51099906 # Electron mass (MeV)
    particleRadius = particleCharge^2 / (4.0 * pi * epsilon_o * particleMass * c_mks^2) # Classical electron radius (m)
    
    skew = [0,0,0]
    dx = dy = dz = 0.0 
    nSlices = integ_order = 0
    i_part = i_top = 0
    coord = []
    drift = 0.0
    tilt = pitch = yaw = rad_coef = isr_coef = xkick = ykick = dzLoss = 0.0
    sextWarning = quadWarning = octWarning = false
    multData = steeringMultData = apData = nothing
    # fringeIntM = [elem.fringeIntM0, elem.fringeIntM1, elem.fringeIntM2, elem.fringeIntM3, elem.fringeIntM4]
    # fringeIntP = [elem.fringeIntP0, elem.fringeIntP1, elem.fringeIntP2, elem.fringeIntP3, elem.fringeIntP4]
    fringeIntM = elem.fringeIntM
    fringeIntP = elem.fringeIntP
    KnL1 = KnL2 = KnL3 = 0.0
    if elem isa KQUAD
        nSlices = elem.nSlices
        order = 1
        if elem.bore != 0.0
            KnL1 = elem.B / elem.bore * (particleCharge / (particleMass * c_mks * Po)) * elem.len * (1.0 + elem.fse)
            # KnL = [KnL1, 0.0, 0.0]
        else
            KnL1 = elem.k1 * elem.len * (1.0 + elem.fse)
            # KnL = [KnL1, 0.0, 0.0]
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
            rad_coef = particleCharge^2 * Po^3 / (6.0 * pi * epsilon_o * (c_mks^2) * particleMass)
        end
        # isr is not implemented yet
        isr_coef = particleRadius *sqrt(55.0 / (24.0 * sqrt(3.0)) * Po^5 * 137.0359895)
        if elem.isr == 0 || (elem.isr1Particle == 0 && n_part == 1 )
            # Minus sign indicates we accumulate into sigmaDelta^2 only, don't perturb particles
            isr_coef *= -1.0
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
        if elem.bore != 0.0
            KnL2 = 2.0 * elem.B / elem.bore^2 * (particleCharge / (particleMass * c_mks * Po)) * elem.len * (1.0 + elem.fse)
            # KnL = [0.0, KnL2, 0.0]
        else
            KnL2 = elem.k2 * elem.len * (1.0 + elem.fse)
            # KnL = [0.0, KnL2, 0.0]
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
            rad_coef = particleCharge^2 * Po^3 / (6.0 * pi * epsilon_o * (c_mks^2) * particleMass)
        end

        # isr is not implemented
        isr_coef = particleRadius *sqrt(55.0 / (24.0 * sqrt(3.0)) * Po^5 * 137.0359895)
        if elem.isr == 0 || (elem.isr1Particle == 0 && n_part == 1 )
            # Minus sign indicates we accumulate into sigmaDelta^2 only, don't perturb particles
            isr_coef *= -1.0
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
        if elem.bore != 0.0
            KnL3 = 6 * elem.B / elem.bore^3 * (particleCharge / (particleMass * c_mks * Po)) * elem.len * (1.0 + elem.fse)
            # KnL = [0.0, 0.0, KnL3]
        else
            KnL3 = elem.k3 * elem.len * (1.0 + elem.fse)
            # KnL = [0.0, 0.0, KnL3]
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
            rad_coef = particleCharge^2 * Po^3 / (6.0 * pi * epsilon_o * (c_mks^2) * particleMass)
        end

        # isr is not implemented
        isr_coef = particleRadius *sqrt(55.0 / (24.0 * sqrt(3.0)) * Po^5 * 137.0359895)
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

    # i_top = n_part -1


    if dx != 0 || dy != 0 || dz != 0 
        x, xp, y, yp, z, delta = offsetBeamCoordinates(x, xp, y, yp, z, delta, n_part, dx, dy, dz)
    else
        x, xp, y, yp, z, delta = x, xp, y, yp, z, delta
    end
    if tilt != 0
        x, xp, y, yp, z, delta = rotateBeamCoordinates(x, xp, y, yp, z, delta, n_part, tilt)
    else
        x, xp, y, yp, z, delta = x, xp, y, yp, z, delta
    end
    


    # if elem isa KQUAD
    #     if elem.edge1_effects !=0
    #         x, xp, y, yp, z, delta = quadFringe(x, xp, y, yp, z, delta, n_part, elem.k1, fringeIntM, fringeIntP, 
    #                                 false, -1, elem.edge1_effects, elem.edge1Linear, elem.edge1NonlinearFactor)
    #     else
    #         x, xp, y, yp, z, delta = x, xp, y, yp, z, delta
    #     end
    # end
  
    if elem isa KQUAD
        x, xp, y, yp, z, delta = integrate_kick_multipole_ordn(x, xp, y, yp, z, delta, dx, dy, xkick, ykick, Po, rad_coef, 
                                      isr_coef, elem.synch_rad, order, KnL1, KnL2, KnL3, skew, nSlices, 
                                      drift, integ_order, multData, 
                                      steeringMultData, apData, dzLoss, elem.radial, tilt)
    else
        x, xp, y, yp, z, delta = integrate_kick_multipole_ordn(x, xp, y, yp, z, delta, dx, dy, xkick, ykick, Po, rad_coef, 
                                      isr_coef, elem.synch_rad, order, KnL1, KnL2, KnL3, skew, nSlices, 
                                      drift, integ_order, multData, 
                                      steeringMultData, apData, dzLoss, 0, tilt)
    end
 


    # if elem isa KQUAD
    #     if elem.edge2_effects !=0
    #         x, xp, y, yp, z, delta = quadFringe(x, xp, y, yp, z, delta, n_part, elem.k1, fringeIntM, fringeIntP, 
    #                                 false, 1, elem.edge2_effects, elem.edge2Linear, elem.edge2NonlinearFactor)
    #     # else
    #     #     x, xp, y, yp, z, delta = x, xp, y, yp, z, delta
    #     end
    # end
    expandHamiltonian = 0
    return x, xp, y, yp, z, delta
end

function pass_TPSA(x::CTPS, xp::CTPS, y::CTPS, yp::CTPS, z::CTPS, delta::CTPS, elem::KSEXT, Po::Float64, sigmaDelta2::Float64)
    # constants
    particleMass = 9.1093897e-31 # Electron mass (kg)
    particleCharge = 1.60217733e-19 # Electron charge (C)
    c_mks = 299792458.0 # Speed of light (m/s)
    epsilon_o = 8.854187817e-12 # Permittivity of vacuum (F/m)
    me_mev = 0.51099906 # Electron mass (MeV)
    particleRadius = particleCharge^2 / (4.0 * pi * epsilon_o * particleMass * c_mks^2) # Classical electron radius (m)
    
    skew = [0,0,0]
    dx = dy = dz = 0.0 
    nSlices = integ_order = 0
    i_part = i_top = 0
    coord = []
    drift = 0.0
    tilt = pitch = yaw = rad_coef = isr_coef = xkick = ykick = dzLoss = 0.0
    sextWarning = quadWarning = octWarning = false
    multData = steeringMultData = apData = nothing
    # fringeIntM = [elem.fringeIntM0, elem.fringeIntM1, elem.fringeIntM2, elem.fringeIntM3, elem.fringeIntM4]
    # fringeIntP = [elem.fringeIntP0, elem.fringeIntP1, elem.fringeIntP2, elem.fringeIntP3, elem.fringeIntP4]
    fringeIntM = elem.fringeIntM
    fringeIntP = elem.fringeIntP
    KnL1 = KnL2 = KnL3 = 0.0
    if elem isa KQUAD
        nSlices = elem.nSlices
        order = 1
        if elem.bore != 0.0
            KnL1 = elem.B / elem.bore * (particleCharge / (particleMass * c_mks * Po)) * elem.len * (1.0 + elem.fse)
            # KnL = [KnL1, 0.0, 0.0]
        else
            KnL1 = elem.k1 * elem.len * (1.0 + elem.fse)
            # KnL = [KnL1, 0.0, 0.0]
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
            rad_coef = particleCharge^2 * Po^3 / (6.0 * pi * epsilon_o * (c_mks^2) * particleMass)
        end
        # isr is not implemented yet
        isr_coef = particleRadius *sqrt(55.0 / (24.0 * sqrt(3.0)) * Po^5 * 137.0359895)
        if elem.isr == 0 || (elem.isr1Particle == 0 && n_part == 1 )
            # Minus sign indicates we accumulate into sigmaDelta^2 only, don't perturb particles
            isr_coef *= -1.0
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
        if elem.bore != 0.0
            KnL2 = 2.0 * elem.B / elem.bore^2 * (particleCharge / (particleMass * c_mks * Po)) * elem.len * (1.0 + elem.fse)
            # KnL = [0.0, KnL2, 0.0]
        else
            KnL2 = elem.k2 * elem.len * (1.0 + elem.fse)
            # KnL = [0.0, KnL2, 0.0]
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
            rad_coef = particleCharge^2 * Po^3 / (6.0 * pi * epsilon_o * (c_mks^2) * particleMass)
        end

        # isr is not implemented
        isr_coef = particleRadius *sqrt(55.0 / (24.0 * sqrt(3.0)) * Po^5 * 137.0359895)
        if elem.isr == 0 || (elem.isr1Particle == 0 && n_part == 1 )
            # Minus sign indicates we accumulate into sigmaDelta^2 only, don't perturb particles
            isr_coef *= -1.0
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
        if elem.bore != 0.0
            KnL3 = 6 * elem.B / elem.bore^3 * (particleCharge / (particleMass * c_mks * Po)) * elem.len * (1.0 + elem.fse)
            # KnL = [0.0, 0.0, KnL3]
        else
            KnL3 = elem.k3 * elem.len * (1.0 + elem.fse)
            # KnL = [0.0, 0.0, KnL3]
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
            rad_coef = particleCharge^2 * Po^3 / (6.0 * pi * epsilon_o * (c_mks^2) * particleMass)
        end

        # isr is not implemented
        isr_coef = particleRadius *sqrt(55.0 / (24.0 * sqrt(3.0)) * Po^5 * 137.0359895)
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

    # i_top = n_part -1


    if dx != 0 || dy != 0 || dz != 0 
        x, xp, y, yp, z, delta = offsetBeamCoordinates(x, xp, y, yp, z, delta, n_part, dx, dy, dz)
    else
        x, xp, y, yp, z, delta = x, xp, y, yp, z, delta
    end
    if tilt != 0
        x, xp, y, yp, z, delta = rotateBeamCoordinates(x, xp, y, yp, z, delta, n_part, tilt)
    else
        x, xp, y, yp, z, delta = x, xp, y, yp, z, delta
    end
    


    # if elem isa KQUAD
    #     if elem.edge1_effects !=0
    #         x, xp, y, yp, z, delta = quadFringe(x, xp, y, yp, z, delta, n_part, elem.k1, fringeIntM, fringeIntP, 
    #                                 false, -1, elem.edge1_effects, elem.edge1Linear, elem.edge1NonlinearFactor)
    #     else
    #         x, xp, y, yp, z, delta = x, xp, y, yp, z, delta
    #     end
    # end
  
    if elem isa KQUAD
        x, xp, y, yp, z, delta = integrate_kick_multipole_ordn(x, xp, y, yp, z, delta, dx, dy, xkick, ykick, Po, rad_coef, 
                                      isr_coef, elem.synch_rad, order, KnL1, KnL2, KnL3, skew, nSlices, 
                                      drift, integ_order, multData, 
                                      steeringMultData, apData, dzLoss, elem.radial, tilt)
    else
        x, xp, y, yp, z, delta = integrate_kick_multipole_ordn(x, xp, y, yp, z, delta, dx, dy, xkick, ykick, Po, rad_coef, 
                                      isr_coef, elem.synch_rad, order, KnL1, KnL2, KnL3, skew, nSlices, 
                                      drift, integ_order, multData, 
                                      steeringMultData, apData, dzLoss, 0, tilt)
    end
 


    # if elem isa KQUAD
    #     if elem.edge2_effects !=0
    #         x, xp, y, yp, z, delta = quadFringe(x, xp, y, yp, z, delta, n_part, elem.k1, fringeIntM, fringeIntP, 
    #                                 false, 1, elem.edge2_effects, elem.edge2Linear, elem.edge2NonlinearFactor)
    #     # else
    #     #     x, xp, y, yp, z, delta = x, xp, y, yp, z, delta
    #     end
    # end
    expandHamiltonian = 0
    return x, xp, y, yp, z, delta
end

function pass_TPSA(x::CTPS, xp::CTPS, y::CTPS, yp::CTPS, z::CTPS, delta::CTPS, elem::KOCT, Po::Float64, sigmaDelta2::Float64)
    # constants
    particleMass = 9.1093897e-31 # Electron mass (kg)
    particleCharge = 1.60217733e-19 # Electron charge (C)
    c_mks = 299792458.0 # Speed of light (m/s)
    epsilon_o = 8.854187817e-12 # Permittivity of vacuum (F/m)
    me_mev = 0.51099906 # Electron mass (MeV)
    particleRadius = particleCharge^2 / (4.0 * pi * epsilon_o * particleMass * c_mks^2) # Classical electron radius (m)
    
    skew = [0,0,0]
    dx = dy = dz = 0.0 
    nSlices = integ_order = 0
    i_part = i_top = 0
    coord = []
    drift = 0.0
    tilt = pitch = yaw = rad_coef = isr_coef = xkick = ykick = dzLoss = 0.0
    sextWarning = quadWarning = octWarning = false
    multData = steeringMultData = apData = nothing
    # fringeIntM = [elem.fringeIntM0, elem.fringeIntM1, elem.fringeIntM2, elem.fringeIntM3, elem.fringeIntM4]
    # fringeIntP = [elem.fringeIntP0, elem.fringeIntP1, elem.fringeIntP2, elem.fringeIntP3, elem.fringeIntP4]
    fringeIntM = elem.fringeIntM
    fringeIntP = elem.fringeIntP
    KnL1 = KnL2 = KnL3 = 0.0
    if elem isa KQUAD
        nSlices = elem.nSlices
        order = 1
        if elem.bore != 0.0
            KnL1 = elem.B / elem.bore * (particleCharge / (particleMass * c_mks * Po)) * elem.len * (1.0 + elem.fse)
            # KnL = [KnL1, 0.0, 0.0]
        else
            KnL1 = elem.k1 * elem.len * (1.0 + elem.fse)
            # KnL = [KnL1, 0.0, 0.0]
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
            rad_coef = particleCharge^2 * Po^3 / (6.0 * pi * epsilon_o * (c_mks^2) * particleMass)
        end
        # isr is not implemented yet
        isr_coef = particleRadius *sqrt(55.0 / (24.0 * sqrt(3.0)) * Po^5 * 137.0359895)
        if elem.isr == 0 || (elem.isr1Particle == 0 && n_part == 1 )
            # Minus sign indicates we accumulate into sigmaDelta^2 only, don't perturb particles
            isr_coef *= -1.0
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
        if elem.bore != 0.0
            KnL2 = 2.0 * elem.B / elem.bore^2 * (particleCharge / (particleMass * c_mks * Po)) * elem.len * (1.0 + elem.fse)
            # KnL = [0.0, KnL2, 0.0]
        else
            KnL2 = elem.k2 * elem.len * (1.0 + elem.fse)
            # KnL = [0.0, KnL2, 0.0]
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
            rad_coef = particleCharge^2 * Po^3 / (6.0 * pi * epsilon_o * (c_mks^2) * particleMass)
        end

        # isr is not implemented
        isr_coef = particleRadius *sqrt(55.0 / (24.0 * sqrt(3.0)) * Po^5 * 137.0359895)
        if elem.isr == 0 || (elem.isr1Particle == 0 && n_part == 1 )
            # Minus sign indicates we accumulate into sigmaDelta^2 only, don't perturb particles
            isr_coef *= -1.0
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
        if elem.bore != 0.0
            KnL3 = 6 * elem.B / elem.bore^3 * (particleCharge / (particleMass * c_mks * Po)) * elem.len * (1.0 + elem.fse)
            # KnL = [0.0, 0.0, KnL3]
        else
            KnL3 = elem.k3 * elem.len * (1.0 + elem.fse)
            # KnL = [0.0, 0.0, KnL3]
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
            rad_coef = particleCharge^2 * Po^3 / (6.0 * pi * epsilon_o * (c_mks^2) * particleMass)
        end

        # isr is not implemented
        isr_coef = particleRadius *sqrt(55.0 / (24.0 * sqrt(3.0)) * Po^5 * 137.0359895)
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

    # i_top = n_part -1


    if dx != 0 || dy != 0 || dz != 0 
        x, xp, y, yp, z, delta = offsetBeamCoordinates(x, xp, y, yp, z, delta, n_part, dx, dy, dz)
    else
        x, xp, y, yp, z, delta = x, xp, y, yp, z, delta
    end
    if tilt != 0
        x, xp, y, yp, z, delta = rotateBeamCoordinates(x, xp, y, yp, z, delta, n_part, tilt)
    else
        x, xp, y, yp, z, delta = x, xp, y, yp, z, delta
    end
    


    # if elem isa KQUAD
    #     if elem.edge1_effects !=0
    #         x, xp, y, yp, z, delta = quadFringe(x, xp, y, yp, z, delta, n_part, elem.k1, fringeIntM, fringeIntP, 
    #                                 false, -1, elem.edge1_effects, elem.edge1Linear, elem.edge1NonlinearFactor)
    #     else
    #         x, xp, y, yp, z, delta = x, xp, y, yp, z, delta
    #     end
    # end
  
    if elem isa KQUAD
        x, xp, y, yp, z, delta = integrate_kick_multipole_ordn(x, xp, y, yp, z, delta, dx, dy, xkick, ykick, Po, rad_coef, 
                                      isr_coef, elem.synch_rad, order, KnL1, KnL2, KnL3, skew, nSlices, 
                                      drift, integ_order, multData, 
                                      steeringMultData, apData, dzLoss, elem.radial, tilt)
    else
        x, xp, y, yp, z, delta = integrate_kick_multipole_ordn(x, xp, y, yp, z, delta, dx, dy, xkick, ykick, Po, rad_coef, 
                                      isr_coef, elem.synch_rad, order, KnL1, KnL2, KnL3, skew, nSlices, 
                                      drift, integ_order, multData, 
                                      steeringMultData, apData, dzLoss, 0, tilt)
    end
 


    # if elem isa KQUAD
    #     if elem.edge2_effects !=0
    #         x, xp, y, yp, z, delta = quadFringe(x, xp, y, yp, z, delta, n_part, elem.k1, fringeIntM, fringeIntP, 
    #                                 false, 1, elem.edge2_effects, elem.edge2Linear, elem.edge2NonlinearFactor)
    #     # else
    #     #     x, xp, y, yp, z, delta = x, xp, y, yp, z, delta
    #     end
    # end
    expandHamiltonian = 0
    return x, xp, y, yp, z, delta
end



# x = CTPS(0.0, 1, 6, 2)
# xp = CTPS(0.0, 2, 6, 2)
# y = CTPS(0.0, 3, 6, 2)
# yp = CTPS(0.0, 4, 6, 2)
# delta = CTPS(0.0, 5, 6, 2)
# z = CTPS(0.0, 6, 6, 2)

# Quad = KQUAD(name="Q",k1=1.0,len=1.0)
# Sext = KSEXT(name="S",k2=2.0,len=1.0)
# Oct = KOCT(name="O",k3=6.0,len=1.0)
# CSB = CSBEND(name="CSB",angle=pi/20/2,len=0.72,e1=pi/20/2,e2=0.0)

# n_part = 1
# xout, xpout, yout, ypout, zout, dpout = multipole_tracking2(x, xp, y, yp, z, delta, n_part, Quad, 1000.0)
# println(xout)
# xvalue = evaluate(xout, [0.001, 0.0001, 0.0005, 0.0002, 0.0, 0.0])
# yvalue = evaluate(yout, [0.001, 0.0001, 0.0005, 0.0002, 0.0, 0.0])
# println(xvalue)
# println(yvalue)
# Map66 = [xout.map[2] xout.map[3] xout.map[4] xout.map[5] xout.map[6] xout.map[7];
# xpout.map[2] xpout.map[3] xpout.map[4] xpout.map[5] xpout.map[6] xpout.map[7];
# 			yout.map[2] yout.map[3] yout.map[4] yout.map[5] yout.map[6] yout.map[7];
# 			ypout.map[2] ypout.map[3] ypout.map[4] ypout.map[5] ypout.map[6] ypout.map[7];
# 			dpout.map[2] dpout.map[3] dpout.map[4] dpout.map[5] dpout.map[6] dpout.map[7];
# 			zout.map[2] zout.map[3] zout.map[4] zout.map[5] zout.map[6] zout.map[7]]
# println(Map66)

# using Enzyme
# # # Enzyme.API.runtimeActivity!(true)

# function f2(xx)
#     L = xx[1]
#     K = xx[2]
#     elem = KQUAD(len=L, k1=K, synch_rad=1)
#     # sext = KSEXT(len=L, k2=K, synch_rad=1)
#     # sext = KSEXT(len=1.0, k2=xx[1],nSlices=4,synch_rad=1)
#     x = CTPS(0.0, 1, 6, 3)
#     xp = CTPS(0.0, 2, 6, 3)
#     y = CTPS(0.0, 3, 6, 3)
#     yp = CTPS(0.0, 4, 6, 3)
#     delta = CTPS(0.0, 5, 6, 3)
#     z = CTPS(0.0, 6, 6, 3)

#     Po = 19569.50762296901

#     x, xp, y, yp, z, delta = pass_TPSA(x, xp, y, yp, z, delta, elem, Po, 0.0)


#     return x.map[2]
# end

# xx = [1.0, 1.0]
# using BenchmarkTools
# println(gradient(Forward, f2, xx))
# @btime f2(xx)