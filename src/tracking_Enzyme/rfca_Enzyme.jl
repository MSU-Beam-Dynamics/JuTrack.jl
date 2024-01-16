include("../lattice/canonical_elements.jl")
 
# 2.99792458e8 = 2.99792458e8

function add_to_particle_energy!(coord, timeOfFlight, Po, dgamma)
    P = Po * (1 + coord[6])
    gamma = sqrt(P^2 + 1)
    gamma1 = gamma + dgamma
    if gamma1<=1
        gamma1 = 1+1e-7
    end
    P1 = sqrt(gamma1^2 - 1)
    coord[6] = P1/Po - 1
    coord[5] = timeOfFlight * 2.99792458e8 * P1 / gamma1

    Pz = P/sqrt(1+coord[2]^2+coord[4]^2)
    Pz1 = sqrt(Pz^2 + gamma1^2 - gamma^2)
    PRatio = Pz/Pz1
    coord[2] *= PRatio
    coord[4] *= PRatio
    return nothing
end

function do_match_energy!(part, np, P_central, change_beam)
    active = 1
    if change_beam == 0
        P_average = 0.0
        if active == 1
            for ip in 1:np
                P_average += P_central*(1 + part[ip,6])
            end
        end
        P_average /= np
        if abs((P_average - P_central)/P_central) > 1e-14
            if active == 1
                # part_new =[[part[ip][1], part[ip][2], part[ip][3], part[ip][4], part[ip][5],
                # ((1 + part[ip][6]) * P_central - P_average) / P_average] for ip in 1:np]
                for ip in 1:np
                    part[ip,6] = ((1+part[ip,6])*P_central - P_average)/P_average
                end
            end
            P_central = P_average
        end
    else
        P_average = 0.0
        if active == 1
            for ip in 1:np
                P_average += P_central*(1 + part[ip,6])
            end
        end
        P_average /= np
        if active == 1
            dP_centroid = P_central - P_average
            for ip in 1:np
                P = (1+part[ip,6])*P_central
                t = part[ip,5]/(P/sqrt(P^2+1.0))
                P += dP_centroid
                # part_new = [[part[ip][1], part[ip][2], part[ip][3], part[ip][4], t*(P/sqrt(P^2+1)), 
                #                 (P-P_central)/P_central] for ip in 1:np]
                part[ip,6] = (P-P_central)/P_central
                part[ip,5] = t*(P/sqrt(P^2+1))
            end
        end
    end
    return P_central
end

function get_phase_reference(rfca_phase_reference, reference_ref_number, reference_phase, reference_flags, n_references)
    FL_REF_PHASE_SET = 1
    if rfca_phase_reference == 0
        return 0, 0.0
    end
    for i in 1:n_references
        if reference_ref_number[i] == rfca_phase_reference
            if reference_flags[i] & FL_REF_PHASE_SET != 0
                phase = reference_phase[i]
                return 1, phase
            end
            return 2, 0.0
        end
    end
    return 3,0.0
end

function parseFiducialMode(modeString)
    FID_MODE_TMEAN = 0x002
    if isnothing(modeString)
        return FID_MODE_TMEAN
    end

    FID_MODE_LIGHT = 0x001
    if modeString == 1 
        return FID_MODE_LIGHT<<0
    elseif modeString == 2
        return FID_MODE_LIGHT<<1
    elseif modeString == 3
        return FID_MODE_LIGHT<<2
    elseif modeString == 4
        return FID_MODE_LIGHT<<3
    else
        return 0
    end
end


function beta_from_delta(p, delta)
    p *= 1.0 + delta
    return p / sqrt(p * p + 1.0)
end

function findFiducialTime(part, np, s0, sOffset, p0, mode)
    FID_MODE_LIGHT = 0x001
    FID_MODE_TMEAN = 0x002
    FID_MODE_FIRST = 0x004
    FID_MODE_PMAX = 0x008
    tFid = 0.0

    if mode & FID_MODE_LIGHT != 0
        tFid = (s0 + sOffset) / 2.99792458e8
    elseif mode & FID_MODE_FIRST != 0
        tFid = (part[1,5] + sOffset) / (2.99792458e8 * beta_from_delta(p0, part[1,6]))

    elseif mode & FID_MODE_PMAX != 0
        best = part[1,6]
        ibest = 1
        for i in 2:np
            if best < part[i,6]
                best = part[i,6]
                ibest = i
            end
        end
        tFid = (part[ibest,5] + sOffset) / (2.99792458e8 * beta_from_delta(p0, part[ibest,6]))
    elseif mode & FID_MODE_TMEAN != 0
        tsum = 0.0
        for ip in 1:np
            tsum += (part[ip,5] + sOffset) / (2.99792458e8 * beta_from_delta(p0, part[ip,6]))
        end
        tFid = tsum / np
    else
        error("Invalid fiducial mode")
    end

    return tFid
end

function set_phase_reference(phase_ref_number, phase, reference_ref_number, reference_phase, reference_flags, n_references)
    if phase_ref_number == 0
        return 0, reference_ref_number, reference_phase, reference_flags
    end
    # reference_ref_number_Buffer = Zygote.Buffer(reference_ref_number)
    # reference_phase_Buffer = Zygote.Buffer(reference_phase)
    # reference_flags_Buffer = Zygote.Buffer(reference_flags)
    # for i in 1:n_references
    #     reference_ref_number_Buffer[i] = reference_ref_number[i]
    #     reference_phase_Buffer[i] = reference_phase[i]
    #     reference_flags_Buffer[i] = reference_flags[i]
    # end
    for i in 1:n_references
        if reference_ref_number[i] == phase_ref_number
            reference_phase[i] = phase
            reference_flags[i] = 1
            return 1, reference_ref_number, reference_phase, reference_flags, phase, n_references
        end
    end
    reference_ref_number[n_references] = phase_ref_number
    reference_phase[n_references] = phase
    reference_flags[n_references] = 1
    n_references += 1
    return 1, reference_ref_number, reference_phase, reference_flags, phase, n_references
end

function pass!(part::Matrix{Float64}, rfca::RFCA, Po::Float64, sigmaDelta2::Float64)
    np = size(part, 1)
    trackRfCavity!(part, np, rfca, Po)
end

function trackRfCavity!(part, np, rfca, P_central)
    particleMassMV = 0.51099906
    particleRelSign = 1 # relative to electron
    

    dgamma=0.0
    dgammaMax=0.0
    been_warned = 0
    dgammaOverGammaAve = 0.0
    dgammaOverGammaNp = 0.0
    lockPhase = 0
    if rfca.bodyFocusModel == "none"
        useSRSModel = 0
        matrixMethod = 0
    elseif  rfca.bodyFocusModel == "srs"
        useSRSModel = 1
        matrixMethod = 1
    else
        error("Unknown body focus model for RFCA: ", rfca.bodyFocusModel)
    end

    # if been_warned == 0
    #     if rfca.freq < 1e3
    #         println("Warning: RFCA frequency is less than 1 kHz.")
    #         been_warned = 1
    #     end
    #     if rfca.volt < 100
    #         println("Warning: RFCA voltage is less than 100 V.")
    #         been_warned = 1
    #     end
    # end

    length = rfca.len
    # if rfca.volt == 0.0 
    #     if length == 0.0
    #         return part, fiducial_seen, P_central, n_references, rfca_phase_reference, reference_ref_number, 
    #         reference_phase, reference_flags
    #     else
    #         for i in 1:np
    #             part[i,1] += part[i,2] * length
    #             part[i,3] += part[i,4] * length
    #             part[i,5] += length * sqrt(1.0 + part[i,2]^2 + part[i,4]^2)
    #         end
    #         return part, fiducial_seen, P_central, n_references, rfca_phase_reference, reference_ref_number, 
    #         reference_phase, reference_flags
    #     end
    # end 
    omega = 2*pi*rfca.freq
    volt = rfca.volt/(1e6 * particleMassMV * particleRelSign)
    if omega != 0 
        tau = rfca.Q/omega
    else
        tau = 0.0
    end

    if length == 0.0
        nKicks = 1
    else
        nKicks = rfca.nSlices
    end

    length /= nKicks
    volt /= nKicks
    dtLight = length/2.99792458e8

    # if rfca.phase_reference == 0
    #     reference_ref_number[1] = 2147483647
    #     reference_phase[1] = 0.0
    #     reference_flags[1] = 0
    #     n_references = n_references + 1
    #     rfca_phase_reference = 2147483647
    # end
    # get phase reference

    if rfca.fiducial_seen == 0
        if rfca.tReference != -1.0
            t0 = rfca.tReference
        else
            sOffset = length/2.0
            tsum = 0.0
            for ip in 1:np
               tsum += (part[ip,5] + sOffset) / (2.99792458e8 * beta_from_delta(P_central, part[ip,6]))
            end
            t0 = tsum / np
        #     t0 = findFiducialTime(part, np, zEnd-rfca.len, length/2, P_central, mode)
        end
        rfca.phase_fiducial = -omega*t0
        rfca.fiducial_seen = 1
    end

    if omega !=0.0
        t0 = -rfca.phase_fiducial/omega
        To = 2*pi/omega
    else
        t0 = 0.0
        To = 0.0
    end

    phase = rfca.phase_fiducial + rfca.phase * pi/180
    
    same_dgamma = 0
    if omega==0.0 && tau==0.0
        dgamma = volt*sin(phase)
        dgammaMax = volt
        same_dgamma = 1
    end
    timeOffset = 0.0

    if omega !=0.0 && rfca.change_t!=0
        coord = part[1,:]
        P = P_central*(1 + coord[6])
        gamma = sqrt(P^2 + 1)
        beta_i = P/gamma
        t = coord[5]/(beta_i*2.99792458e8)
        if omega!=0.0 && t>(0.9*To) && rfca.change_t!=0
            timeOffset = round(Int, t / To) * To
        end
    end

    # not use linearize and lockPhase
    # linearize = rfca.linearize
    # lockPhase = rfca.lockPhase
    # if linearize!=0 || lockPhase!=0
    #     tAve = 0
    #     if nKicks!=1
    #         error("Must use nKicks=1 for linearized RFCA")
    #     end
    #     for i in 1:np
    #         coord = part[i,:]
    #         P = P_central*(1 + coord[6])
    #         gamma = sqrt(P^2 + 1)
    #         beta_i = P/gamma
    #         t = (coord[5]+length/2)/(beta_i*2.99792458e8) - timeOffset
    #         tAve += t
    #     end
    #     tAve /= np
    #     if lockPhase!=0
    #         phase = pi/180*rfca.phase
    #         dgammaAve = volt*sin(phase)
    #     else
    #         dgammaAve = volt*sin(omega*(tAve-timeOffset)+phase)
    #     end
    # end
    
    for ip in 1:np
        part[ip,1] -= rfca.dx
        part[ip,3] -= rfca.dy
    end

    if matrixMethod==0
        inverseF = zeros(np)
        for ik in 0:nKicks-1
            dgammaOverGammaAve = 0.0
            dgammaOverGammaNp = 0.0
            for ip in 1:np
                coord = part[ip,:]
                # dx, dy
                if coord[6] == -1.0
                    continue
                end
                if length!=0
                    dc4 = length/2*sqrt(1+coord[2]^2+coord[4]^2)
                else
                    dc4 = 0.0
                end

                # energy kick
                P = P_central*(1.0 + coord[6])
                gamma = sqrt(P^2 + 1.0)
                beta_i = P/gamma
                t = (coord[5]+dc4)/(beta_i*2.99792458e8) - timeOffset
                if ik==0 && timeOffset!=0 && rfca.change_t!=0
                    coord[5] = t*beta_i*2.99792458e8 - dc4
                end
                dt = t-t0
                if dt < 0
                    dt = 0.0
                end
                if same_dgamma==0
                    # if linearize ==0 # not use linearize
                        # if lockPhase == 1
                        #     temp1 = tAve
                        # else
                            temp1 = 0.0
                        # end
                        if tau == 0
                            temp2 = 1.0
                        else
                            temp2 = sqrt(1.0 - exp(-dt / tau))
                        end
                        dgamma = volt * sin(omega * (t - temp1 - ik * dtLight) + phase) * temp2
                        # dgamma = volt * sin(omega * (t - (lockPhase==1 ? tAve : 0.0) - ik * dtLight) + phase) * 
                        #             (tau!=0 ? sqrt(1.0 - exp(-dt / tau)) : 1.0)
                    # else
                    #     dgamma = dgammaAve + volt*omega*(t-tAve)*cos(omega*(tAve-timeOffset)+phase)
                    # end
                end
                if gamma!=0
                    dgammaOverGammaNp += 1.0
                    dgammaOverGammaAve += dgamma/gamma
                end
                if length != 0
                    if rfca.end1Focus!=0 && ik==0
                        # focus kick
                        inverseF[ip] = dgamma/(2*gamma*length)
                        coord[2] -= coord[1]*inverseF[ip]
                        coord[4] -= coord[3]*inverseF[ip]
                    else
                        inverseF[ip] = 0
                    end
                    # initial drift
                    coord[1] += coord[2]*length/2
                    coord[3] += coord[4]*length/2
                    coord[5] += length/2*sqrt(1+coord[2]^2+coord[4]^2)
                else
                    inverseF[ip] = 0
                end
                # energy kick
                add_to_particle_energy!(coord, t, P_central, dgamma)
                gamma1 = gamma+dgamma
                if gamma1<=1
                    coord[6] = -1
                else
                    inverseF[ip] = -dgamma/(2*gamma1*length)
                end
                part[ip,:] = coord
            end

            if length!=0
                for ip in 1:np
                    coord = part[ip,:]
                    coord[1] += coord[2]*length/2
                    coord[3] += coord[4]*length/2
                    coord[5] += length/2*sqrt(1.0+coord[2]^2+coord[4]^2)

                    if rfca.end2Focus!=0 && ik==nKicks-1
                        # focus kick
                        coord[2] -= coord[1]*inverseF[ip]
                        coord[4] -= coord[3]*inverseF[ip]
                    end
                    part[ip,:] = coord
                end
            end
        end
    else
        # matrix method
        sin_phase = 0.0
        R11 = 1.0
        R21 = 0.0
        for ik in 0:nKicks-1
            dgammaOverGammaAve = 0.0
            dgammaOverGammaNp = 0.0
            for ip in 1:np
                coord = part[ip,:]
                # use matrix to propagate particles
                P = P_central*(1.0 + coord[6])
                gamma = sqrt(P^2 + 1.0)
                beta_i = P/gamma
                ds1 = length/2*sqrt(1.0+coord[2]^2+coord[4]^2)
                t = (coord[5]+ds1)/(beta_i*2.99792458e8) - timeOffset
                if timeOffset!=0 && rfca.change_t!=0
                    coord[5] = t*beta_i*2.99792458e8 - ds1
                end
                dt = t-t0
                if dt < 0
                    dt = 0.0
                end
                if same_dgamma==0
                    # if linearize ==0 # not use linearize
                        sin_phase = sin(omega * (t - 0.0 - ik * dtLight) + phase)
                        cos_phase = cos(omega * (t - 0.0 - ik * dtLight) + phase)
                        if tau == 0 
                            temp3 = 1.0
                        else
                            temp3 = sqrt(1.0 - exp(-dt / tau))
                        end
                        dgammaMax = volt * temp3
                        # dgammaMax = volt * (tau!=0 ? sqrt(1 - exp(-dt / tau)) : 1)
                        dgamma = dgammaMax * sin_phase
                    # else
                    #     cos_phase = cos(omega * tAvet + phase)
                    #     sin_phase = omega*(t-tAve)*cos_phase
                    #     if tau == 0
                    #         temp3 = 1.0
                    #     else
                    #         temp3 = sqrt(1.0 - exp(-dt / tau))
                    #     end
                    #     dgammaMax = volt * temp3
                    #     # dgammaMax = volt*(tau!=0 ? sqrt(1 - exp(-dt / tau)) : 1)
                    #     dgamma = dgammaAve + dgammaMax*sin_phase
                    # end
                end

                if rfca.end1Focus!=0 && ik==0 && length!=0
                    # focus kick
                    inverseF1 = dgamma/(2*gamma*length)
                    coord[2] -= coord[1]*inverseF1
                    coord[4] -= coord[3]*inverseF1
                end

                dP = sqrt((gamma+dgamma)^2 - 1.0) - P
                if gamma!=0
                    dgammaOverGammaNp += 1.0
                    dgammaOverGammaAve += dgamma/gamma
                end

                if useSRSModel != 0 
                    gammaf = gamma+dgamma
                    if abs(sin_phase) > 1e-6
                        alpha = log(gammaf/gamma)/(2*sqrt(2)*sin_phase)
                    else
                        alpha = dgammaMax/gamma/(2*sqrt(2))
                    end
                    R11 = cos(alpha)
                    R22 = R11*gamma/gammaf
                    R12 = 2*sqrt(2)*gamma*length/dgammaMax*sin(alpha)
                    R21 = -sin(alpha)*dgammaMax/(length*gammaf*2*sqrt(2))
                else
                    R22 = 1.0/(1.0+dP/P)
                    if abs(dp/p) > 1e-14
                        R12 = length*(P/dP*log(1.0+dP/P))
                    else
                        R12 = length
                    end
                end

                coord[5] += ds1
                x = coord[1]
                xp = coord[2]
                coord[1] = R11*x + R12*xp
                coord[2] = R21*x + R22*xp
                x = coord[3]
                xp = coord[4]
                coord[3] = R11*x + R12*xp
                coord[4] = R21*x + R22*xp
                coord[5] += length/2*sqrt(1+coord[2]^2+coord[4]^2)
                coord[6] = (P+dP-P_central)/P_central

                gamma += dgamma
                if gamma <= 1
                    coord[6] = -1
                end
                if rfca.end2Focus!=0 && ik==nKicks-1 && length!=0
                    # focus kick
                    inverseF1 = -dgamma/(2*gamma*length)
                    coord[2] -= coord[1]*inverseF1
                    coord[4] -= coord[3]*inverseF1
                end
                coord[5] = (P+dP)/gamma*coord[5]/beta_i
                part[ip,:] = coord
            end
        end
    end

    for ip in 1:np
        coord = part[ip,:]
        coord[1] += rfca.dx
        coord[3] += rfca.dy
        part[ip,:] = coord
    end
    if rfca.change_p0 != 0
        P_central = do_match_energy!(part, np, P_central, 0)
    end
    return nothing
    # return part, fiducial_seen, P_central, n_references, rfca_phase_reference, reference_ref_number, 
            #reference_phase, reference_flags
end


# function f(xx)
#     volt = xx[1]
#     freq = xx[2]
#     rfca = RFCA(name="rfca", freq=freq, volt=volt, phase=90.0, len=0.1, nSlices=10)
#     parts = [0.001 0.0001 0.0005 0.0002 0.0 0.0; 0.001 0.0 0.0 0.0 0.0 0.0]
#     Po = 19569.50762296901
#     pass!(parts, rfca, Po, 0.0)
#     # println((1+parts[1,6])*Po)
#     return parts[1,:]
# end
# volt = 1e6
# freq = 500e6
# println(f([volt, freq]))
# using Enzyme
# using BenchmarkTools
# @btime begin
#     # f([volt, freq])
#     # grad = jacobian(Reverse, f, [volt, freq], Val(6))
# end
# println(grad)