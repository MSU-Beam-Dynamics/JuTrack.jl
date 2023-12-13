include("EDrift.jl")
include("../lattice/canonical_elements.jl")
include("../TPSA/TPSA.jl")
using Zygote
 
const c_mks = 2.99792458e8

function add_to_particle_energy(x, xp, y, yp, z, delta, timeOfFlight, Po, dgamma)
    P = Po * (1 + delta)
    gamma = sqrt(P^2 + 1)
    gamma1 = gamma + dgamma
    # if gamma1<=1
    #     gamma1 = 1+1e-7
    # end
    P1 = sqrt(gamma1^2 - 1)
    coord_6 = P1/Po - 1
    coord_5 = timeOfFlight * c_mks * P1 / gamma1

    Pz = P/sqrt(1+xp^2+yp^2)
    Pz1 = sqrt(Pz^2 + gamma1^2 - gamma^2)
    PRatio = Pz/Pz1
    coord_2 = xp*PRatio
    coord_4 = yp*PRatio
    return x, coord_2, y, coord_4, coord_5, coord_6
end

function do_match_energy(z, delta, np, P_central, change_beam)
    active = 1
    if change_beam == 0
        P_average = 0
        if active == 1
            for ip in 1:np
                P_average += P_central*(1 + delta)
            end
        end
        P_average /= np
        if abs((P_average - P_central)/P_central) > 1e-14
            if active == 1
                delta = ((1 + delta) * P_central - P_average) / P_average
                # for ip in 1:np
                #     part[ip][6] = ((1+part[ip][6])*P_central - P_average)/P_average
                # end
            end
            P_central = P_average
        end
    else
        P_average = 0
        if active == 1
            for ip in 1:np
                P_average += P_central*(1 + delta)
            end
        end
        P_average /= np
        if active == 1
            dP_centroid = P_central - P_average
            for ip in 1:np
                P = (1+delta)*P_central
                t = z/(P/sqrt(P^2+1))
                P += dP_centroid
                z = t*(P/sqrt(P^2+1))
                delta = (P-P_central)/P_central
                # part[ip][6] = (P-P_central)/P_central
                # part[ip][5] = t*(P/sqrt(P^2+1))
            end
        end
    end
    return z, delta, P_central
end

function get_phase_reference(rfca_phase_reference, reference_ref_number, reference_phase, reference_flags, n_references)
    FL_REF_PHASE_SET = 1
    if rfca_phase_reference == 0
        return 0, 0
    end
    for i in 1:n_references
        if reference_ref_number[i] == rfca_phase_reference
            if reference_flags[i] & FL_REF_PHASE_SET != 0
                phase = reference_phase[i]
                return 1, phase
            end
            return 2, 0
        end
    end
    return 3,0
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
    p *= 1 + delta
    return p / sqrt(p * p + 1)
end

function findFiducialTime(x, xp, y, yp, z, delta, s0, sOffset, p0, mode)
    FID_MODE_LIGHT = 0x001
    FID_MODE_TMEAN = 0x002
    FID_MODE_FIRST = 0x004
    FID_MODE_PMAX = 0x008
    tFid = 0.0
    np = 1

    if mode & FID_MODE_LIGHT != 0
        tFid = (s0 + sOffset) / c_mks
    elseif mode & FID_MODE_FIRST != 0
        tFid = (z + sOffset) / (c_mks * beta_from_delta(p0, delta))

    elseif mode & FID_MODE_PMAX != 0
        best = delta
        ibest = 1
        for i in 2:np
            if best < delta
                best = delta
                ibest = i
            end
        end
        tFid = (z + sOffset) / (c_mks * beta_from_delta(p0,delta))
    elseif mode & FID_MODE_TMEAN != 0
        tsum = 0.0
        for ip in 1:np
            tsum += (z + sOffset) / (c_mks * beta_from_delta(p0, delta))
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
    reference_ref_number_Buffer = Zygote.Buffer(reference_ref_number)
    # reference_phase_Buffer = Zygote.Buffer(reference_phase)
    reference_flags_Buffer = Zygote.Buffer(reference_flags)
    for i in 1:n_references
        reference_ref_number_Buffer[i] = reference_ref_number[i]
        # reference_phase_Buffer[i] = reference_phase[i]
        reference_flags_Buffer[i] = reference_flags[i]
    end
    for i in 1:n_references
        if reference_ref_number[i] == phase_ref_number
            reference_phase_Buffer = [i == j ? phase : CTPS(reference_phase[j], 0, 6, 2) for j in 1:n_references]
            reference_flags_Buffer[i] = 1
            return 1, reference_ref_number, reference_phase_Buffer, copy(reference_flags_Buffer), phase, n_references
        end
    end
    reference_ref_number_Buffer[n_references] = phase_ref_number
    reference_phase_Buffer = [i==n_references ? phase : reference_phase[i] for i in 1:n_references]
    # reference_phase_Buffer[n_references] = phase
    reference_flags_Buffer[n_references] = 1
    # reference_ref_number[n_references] = phase_ref_number
    # reference_phase[n_references] = phase
    # reference_flags[n_references] = 1
    n_references += 1
    return 1, copy(reference_ref_number_Buffer), copy(reference_phase_Buffer), copy(reference_flags_Buffer), phase, n_references
end

function simple_rf_cavity(x, xp, y, yp, z, delta, rfca, Po, zEnd, fiducial_seen, rfca_phase_reference,
     reference_ref_number, reference_phase, reference_flags, n_references)
    return trackRfCavityWithWakes(x, xp, y, yp, z, delta, rfca, Po, zEnd, 
            nothing, nothing, nothing, nothing, 
            fiducial_seen,  rfca_phase_reference, reference_ref_number, reference_phase, reference_flags, n_references)
end

function trackRfCavityWithWakes(x, xp, y, yp, z, delta, rfca, P_central, zEnd, wake, trwake, LSCKick, wakesAtEnd,
    fiducial_seen, rfca_phase_reference, reference_ref_number, reference_phase, reference_flags, n_references)
    particleMassMV = 0.51099906
    particleRelSign = 1 # relative to electron
    np = 1

    dgamma=0.0
    dgammaMax=0.0
    been_warned = 0
    dgammaOverGammaAve = 0
    dgammaOverGammaNp = 0
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

    if been_warned == 0
        if rfca.freq < 1e3
            println("Warning: RFCA frequency is less than 1 kHz.")
            been_warned = 1
        end
        if rfca.volt < 100
            println("Warning: RFCA voltage is less than 100 V.")
            been_warned = 1
        end
    end

    length = rfca.len
    if rfca.volt == 0 && isnothing(wake) && isnothing(trwake) && isnothing(LSCKick)
        if length == 0
            return x, xp, y, yp, z, delta, fiducial_seen, P_central, n_references, rfca_phase_reference,
            reference_ref_number, reference_phase, reference_flags
        else
            return EDrift(x, xp, y, yp, z, delta, length), fiducial_seen, P_central, n_references, rfca_phase_reference,
            reference_ref_number, reference_phase, reference_flags
        end
    end 
    omega = 2*pi*rfca.freq
    volt = rfca.volt/(1e6 * particleMassMV * particleRelSign)
    if omega != 0 
        tau = rfca.Q/omega
    else
        tau = 0
    end

    if length == 0
        nKicks = 1
    else
        nKicks = rfca.nSlices
    end

    length /= nKicks
    volt /= nKicks
    dtLight = length/c_mks

    if rfca.phase_reference == 0
        reference_ref_number = [2147483647]
        reference_phase = [0.0]
        reference_flags = [0]
        n_references += 1
        rfca_phase_reference = 2147483647
    end
    # get phase reference
    ref_flag, phase = get_phase_reference(rfca_phase_reference, reference_ref_number, reference_phase, reference_flags, n_references)
    if ref_flag == 0
        error("Phase reference not found for RFCA: ", rfca.name)
    elseif ref_flag == 3 || ref_flag == 2
        if fiducial_seen == 0
            mode = parseFiducialMode(rfca.fiducial_mode)
            if mode == 0
                error("Invalid fiducial mode for RFCA: ", rfca.name)
            end
            if rfca.tReference != -1
                t0 = rfca.tReference
            else
                t0 = findFiducialTime(x, xp, y, yp, z, delta, zEnd-rfca.len, length/2, P_central, mode)
            end
            phase_fiducial = -omega*t0
            fiducial_seen = 1
        end
        set_phase_reference_flag, reference_ref_number, reference_phase, reference_flags, phase, n_references = set_phase_reference(rfca_phase_reference, 
                        phase_fiducial, reference_ref_number, reference_phase, reference_flags, n_references)
    end

    if omega !=0
        t0 = -phase_fiducial/omega
        To = 2*pi/omega
    else
        t0 = 0
        To = 0
    end

    phase = phase + rfca.phase * pi/180
    same_dgamma = 0
    if omega==0 && tau==0
        dgamma = volt*sin(phase)
        dgammaMax = volt
        same_dgamma = 1
    end
    timeOffset = 0

    if omega !=0 && rfca.change_t!=0
        # coord = part[1]
        P = P_central*(1 + delta)
        gamma = sqrt(P^2 + 1)
        beta_i = P/gamma
        t = z/(beta_i*c_mks)
        if omega!=0 && t>(0.9*To) && rfca.change_t!=0
            timeOffset = round(Int, t / To) * To
        end
    end
    linearize = rfca.linearize
    lockPhase = rfca.lockPhase
    if linearize!=0 || lockPhase!=0
        tAve = 0
        if nKicks!=1
            error("Must use nKicks=1 for linearized RFCA")
        end
        for i in 1:np
            # coord = part[i]
            P = P_central*(1 + delta)
            gamma = sqrt(P^2 + 1)
            beta_i = P/gamma
            t = (z+length/2)/(beta_i*c_mks) - timeOffset
            tAve += t
        end
        tAve /= np
        if lockPhase!=0
            phase = pi/180*rfca.phase
            dgammaAve = volt*sin(phase)
        else
            dgammaAve = volt*sin(omega*(tAve-timeOffset)+phase)
        end
    end

    # part_new = [[part[ip][1]-rfca.dx, part[ip][2], part[ip][3]-rfca.dy, part[ip][4], part[ip][5], part[ip][6]] for ip in 1:np]
    # for ip in 1:np
    #     coord = part[ip]
    #     coord[1] -= rfca.dx
    #     coord[3] -= rfca.dy
    # end

    if matrixMethod==0
        inverseF = 0.0
        # inverseF_Buffer = Zygote.Buffer(inverseF)
        for ik in 0:nKicks-1
            dgammaOverGammaAve = 0
            dgammaOverGammaNp = 0
            for ip in 1:np
                # dx, dy
                x -= rfca.dx
                y -= rfca.dy
                # if delta == -1
                #     continue
                # end
                if length!=0
                    dc4 = length/2*sqrt(1+xp^2+yp^2)
                else
                    dc4 = 0
                end

                # energy kick
                P = P_central*(1 + delta)
                gamma = sqrt(P^2 + 1)
                beta_i = P/gamma
                t = (z+dc4)/(beta_i*c_mks) - timeOffset
                if ik==0 && timeOffset!=0 && rfca.change_t!=0
                    z = t*beta_i*c_mks - dc4
                end
                dt = t-t0
                # if dt < 0
                #     dt = 0
                # end
                if same_dgamma==0
                    if linearize ==0
                        dgamma = volt * sin(omega * (t - (lockPhase==1 ? tAve : 0) - ik * dtLight) + phase) * 
                                    (tau!=0 ? sqrt(1 - exp(-dt / tau)) : 1)
                    else
                        dgamma = dgammaAve + volt*omega*(t-tAve)*cos(omega*(tAve-timeOffset)+phase)
                    end
                end
                if gamma!=0
                    dgammaOverGammaNp += 1
                    dgammaOverGammaAve += dgamma/gamma
                end
                if length != 0
                    if rfca.end1Focus!=0 && ik==0
                        # focus kick
                        inverseF = dgamma/(2*gamma*length)
                        xp -= x*inverseF
                        yp -= y*inverseF
                    else
                        inverseF = 0
                    end
                    # initial drift
                    x += xp*length/2
                    y += yp*length/2
                    z += length/2*sqrt(1+xp^2+yp^2)
                else
                    inverseF = 0
                end
                # energy kick
                x, xp, y, yp, z, delta = add_to_particle_energy(x, xp, y, yp, z, delta, t, P_central, dgamma)
                gamma1 = gamma+dgamma
                if cst(gamma1)<=1
                    delta = assian(delta, -1)
                else
                    if length !=0
                        inverseF = -dgamma/(2*gamma1*length)
                    else
                        inverseF = Inf
                    end
                end
            end

            if wakesAtEnd==0
                if !isnothing(wake)
                    #####################################################
                end
                if !isnothing(trwake)
                    #####################################################
                end
                if !isnothing(LSCKick)
                    #####################################################
                    if dgammaOverGammaNp!=0
                        dgammaOverGammaAve /= dgammaOverGammaNp
                    end
                ##########################addLSCKick
                end
            end
            if length!=0
                for ip in 1:np
                    x += xp*length/2
                    y += yp*length/2
                    z += length/2*sqrt(1+xp^2+yp^2)
                    if rfca.end2Focus!=0 && ik==nKicks-1
                        # focus kick
                        xp = xp - x*inverseF
                        yp = yp - y*inverseF
                        # coord[2] -= coord[1]*inverseF_Buffer[ip]
                        # coord[4] -= coord[3]*inverseF_Buffer[ip]
                    end
                end
            end
            if wakesAtEnd!=0
                if !isnothing(wake)
                    #####################################################
                end
                if !isnothing(trwake)
                    #####################################################
                end
                if !isnothing(LSCKick)
                    #####################################################
                    if dgammaOverGammaNp!=0
                        dgammaOverGammaAve /= dgammaOverGammaNp
                    end
                ##########################addLSCKick
                end
            end
        end
    else
        # matrix method
        sin_phase = 0.0
        R11 = 1
        R21 = 0
        for ik in 0:nKicks-1
            dgammaOverGammaAve = 0
            dgammaOverGammaNp = 0
            for ip in 1:np
                # use matrix to propagate particles
                P = P_central*(1 + delta)
                gamma = sqrt(P^2 + 1)
                beta_i = P/gamma
                ds1 = length/2*sqrt(1+xp^2+yp^2)
                t = (z+ds1)/(beta_i*c_mks) - timeOffset
                if timeOffset!=0 && rfca.change_t!=0
                    z = t*beta_i*c_mks - ds1
                end
                dt = t-t0
                if dt < 0
                    dt = 0
                end
                if same_dgamma==0
                    if linearize ==0
                        sin_phase = sin(omega * (t - (lockPhase==1 ? tAve : 0) - ik * dtLight) + phase)
                        cos_phase = cos(omega * (t - (lockPhase==1 ? tAve : 0) - ik * dtLight) + phase)
                        dgammaMax = volt * (tau!=0 ? sqrt(1 - exp(-dt / tau)) : 1)
                        dgamma = dgammaMax * sin_phase
                    else
                        cos_phase = cos(omega * tAvet + phase)
                        sin_phase = omega*(t-tAve)*cos_phase
                        dgammaMax = volt*(tau!=0 ? sqrt(1 - exp(-dt / tau)) : 1)
                        dgamma = dgammaAve + dgammaMax*sin_phase
                    end
                end

                if rfca.end1Focus!=0 && ik==0 && length!=0
                    # focus kick
                    inverseF = dgamma/(2*gamma*length)
                    xp -= x*inverseF
                    yp -= y*inverseF
                end

                dP = sqrt((gamma+dgamma)^2 - 1) - P
                if gamma!=0
                    dgammaOverGammaNp += 1
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
                    R22 = 1/(1+dP/P)
                    if abs(dp/p) > 1e-14
                        R12 = length*(P/dP*log(1+dP/P))
                    else
                        R12 = length
                    end
                end

                z += ds1
                x_temp = copy(x)
                xp_temp = copy(xp)
                x = R11*x_temp + R12*xp_temp
                xp = R21*x_temp + R22*xp_temp
                x_temp = copy(y)
                xp_temp = copy(yp)
                y = R11*x_temp + R12*xp_temp
                yp = R21*x_temp + R22*xp_temp
                z += length/2*sqrt(1+xp^2+xp^2)
                delta = (P+dP-P_central)/P_central

                gamma += dgamma
                if gamma <= 1
                    delta = assign(delta, -1)
                end
                if rfca.end2Focus!=0 && ik==nKicks-1 && length!=0
                    # focus kick
                    inverseF = -dgamma/(2*gamma*length)
                    xp -= x*inverseF
                    yp -= y*inverseF
                end
                z = (P+dP)/gamma*z/beta_i
            end
            if wakesAtEnd==0
                if !isnothing(wake)
                    #####################################################
                end
                if !isnothing(trwake)
                    #####################################################
                end
                if !isnothing(LSCKick)
                    #####################################################
                    if dgammaOverGammaNp!=0
                        dgammaOverGammaAve /= dgammaOverGammaNp
                    end
                ##########################addLSCKick
                end
            end
        end
    end
    for ip in 1:np
        x = x + rfca.dx
        y = y + rfca.dy
    end
    if rfca.change_p0 != 0
        z, delta, P_central = do_match_energy(z, delta, np, P_central, 0)
    end
    return x, xp, y, yp, z, delta, fiducial_seen, P_central, n_references, rfca_phase_reference,
            reference_ref_number, reference_phase, reference_flags
end


# test with Zygote

# rfca = RFCA("rfca", 500e6, 1e6, 90.0, 0, 0, nothing, -1.0, 0, 0.0, 10, 0, 0, 0, 0, 0, 0, 0, 0, "none")
# x = CTPS(0.0, 1, 6, 2)
# xp = CTPS(0.0, 2, 6, 2)
# y = CTPS(0.0, 3, 6, 2)
# yp = CTPS(0.0, 4, 6, 2)
# delta = CTPS(0.0, 5, 6, 2)
# z = CTPS(0.0, 6, 6, 2)
# x, xp, y, yp, z, delta, fiducial_seen, P_central, n_references, reference_ref_number, reference_phase, reference_flags= simple_rf_cavity(x, xp, y, yp, z, delta, 
#                          rfca, 19569.507622969009, 0.0, 0, 0, 0, [0.0], [0], 0)
# println(x)
# println(xp)
# println(y)
# println(yp)
# println(z)
# println(delta)

# xvalue = evaluate(x, [0.001, 0.0001, 0.0005, 0.0002, 0.0, 0.0])
# yvalue = evaluate(y, [0.001, 0.0001, 0.0005, 0.0002, 0.0, 0.0])
# println(xvalue)
# println(yvalue)
# Map66 = [x.map[2] x.map[3] x.map[4] x.map[5] x.map[6] x.map[7];
# xp.map[2] xp.map[3] xp.map[4] xp.map[5] xp.map[6] xp.map[7];
# 			y.map[2] y.map[3] y.map[4] y.map[5] y.map[6] y.map[7];
# 			yp.map[2] yp.map[3] yp.map[4] yp.map[5] yp.map[6] yp.map[7];
# 			delta.map[2] delta.map[3] delta.map[4] delta.map[5] delta.map[6] delta.map[7];
# 			z.map[2] z.map[3] z.map[4] z.map[5] z.map[6] z.map[7]]
# println(Map66)
