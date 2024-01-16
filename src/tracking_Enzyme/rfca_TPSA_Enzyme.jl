include("../lattice/canonical_elements.jl")
include("../TPSA_Enzyme/TPSA_fixedmap.jl")
# 2.99792458e8 = 2.99792458e8

function add_to_particle_energy(x::CTPS{T, TPS_Dim, Max_TPS_Degree}, xp::CTPS{T, TPS_Dim, Max_TPS_Degree}, y::CTPS{T, TPS_Dim, Max_TPS_Degree}, 
    yp::CTPS{T, TPS_Dim, Max_TPS_Degree}, z::CTPS{T, TPS_Dim, Max_TPS_Degree}, delta::CTPS{T, TPS_Dim, Max_TPS_Degree}, timeOfFlight, Po, dgamma) where {T, TPS_Dim, Max_TPS_Degree}
    P = Po * (1.0 + delta)
    gamma = sqrt(P^2 + 1.0)
    gamma1 = gamma + dgamma
    if gamma1.map[1]<=1.0
        gamma1.map[1] = 1.0+1e-7
    end
    P1 = sqrt(gamma1^2 - 1.0)
    delta_new = P1/Po - 1.0
    z_new = timeOfFlight * 2.99792458e8 * P1 / gamma1

    Pz = P/sqrt(1.0+xp^2+yp^2)
    Pz1 = sqrt(Pz^2 + gamma1^2 - gamma^2)
    PRatio = Pz/Pz1
    xp *= PRatio
    yp *= PRatio
    return xp, yp, z, delta
end

function do_match_energy(x::CTPS{T, TPS_Dim, Max_TPS_Degree}, xp::CTPS{T, TPS_Dim, Max_TPS_Degree}, y::CTPS{T, TPS_Dim, Max_TPS_Degree}, 
    yp::CTPS{T, TPS_Dim, Max_TPS_Degree}, z::CTPS{T, TPS_Dim, Max_TPS_Degree}, delta::CTPS{T, TPS_Dim, Max_TPS_Degree}, np::Int, P_central::Float64, change_beam::Int) where {T, TPS_Dim, Max_TPS_Degree}
    active = 1
    if change_beam == 0
        P_average = CTPS(0.0, TPS_Dim, Max_TPS_Degree)
        if active == 1
            # for ip in 1:np
                P_average += P_central*(1.0 + delta)
            # end
        end
        # P_average /= np
        if abs(((P_average - P_central)/P_central).map[1]) > 1e-14
            if active == 1
                # part_new =[[part[ip][1], part[ip][2], part[ip][3], part[ip][4], part[ip][5],
                # ((1 + part[ip][6]) * P_central - P_average) / P_average] for ip in 1:np]
                # for ip in 1:np
                    delta = ((1.0+delta)*P_central - P_average)/P_average
                # end
            end
            P_central = P_average
        end
    else
        P_average = 0.0
        if active == 1
            # for ip in 1:np
                P_average += P_central*(1.0 + delta)
            # end
        end
        # P_average /= np
        if active == 1
            dP_centroid = P_central - P_average
            # for ip in 1:np
                P = (1.0+delta)*P_central
                t = z/(P/sqrt(P^2+1.0))
                P += dP_centroid
                # part_new = [[part[ip][1], part[ip][2], part[ip][3], part[ip][4], t*(P/sqrt(P^2+1)), 
                #                 (P-P_central)/P_central] for ip in 1:np]
                delta = (P-P_central)/P_central
                z = t*(P/sqrt(P^2+1.0))
            # end
        end
    end
    return x, xp, y, yp, z, delta, P_central
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


function pass_TPSA(x::CTPS{T, TPS_Dim, Max_TPS_Degree}, xp::CTPS{T, TPS_Dim, Max_TPS_Degree}, y::CTPS{T, TPS_Dim, Max_TPS_Degree}, 
                yp::CTPS{T, TPS_Dim, Max_TPS_Degree}, z::CTPS{T, TPS_Dim, Max_TPS_Degree}, delta::CTPS{T, TPS_Dim, Max_TPS_Degree}, 
                rfca::RFCA, P_central::Float64, sigmaDelta2::Float64) where {T, TPS_Dim, Max_TPS_Degree}
    particleMassMV = 0.51099906
    particleRelSign = 1 # relative to electron
    np = 1

    dgamma=CTPS(0.0, TPS_Dim, Max_TPS_Degree)
    dgammaMax=CTPS(0.0, TPS_Dim, Max_TPS_Degree)
    been_warned = 0
    # dgammaOverGammaAve = 0.0
    # dgammaOverGammaNp = 0.0
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


    length = rfca.len

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

    # if rfca.fiducial_seen == 0
        if rfca.tReference != -1.0
            error("does not support tReference")#t0 = rfca.tReference
        end
        # else
            sOffset = length/2.0
            tsum = (z + sOffset) / (2.99792458e8 * beta_from_delta(P_central, delta))
            t0 = tsum #/ np # np is 1
        # rfca.phase_fiducial = -omega*t0
        phase_fiducial = -omega*t0
        rfca.fiducial_seen = 1
    # end

    if omega !=0.0
        t0 = -phase_fiducial/omega
        To = 2*pi/omega
    else
        t0 = CTPS(0.0, TPS_Dim, Max_TPS_Degree)
        To = 0.0
    end

    phase = phase_fiducial + rfca.phase * pi/180
    
    same_dgamma = 0
    if omega==0.0 && tau==0.0
        dgamma = sin(phase) * volt
        dgammaMax = CTPS(volt, TPS_Dim, Max_TPS_Degree)
        same_dgamma = 1
    end
    timeOffset = 0.0

    # if omega !=0.0 && rfca.change_t!=0 # does not support change_t
    #     # coord = part[1,:]
    #     P = P_central*(1 + delta)
    #     gamma = sqrt(P^2 + 1)
    #     beta_i = P/gamma
    #     t = z/(beta_i*2.99792458e8)
    #     if omega!=0.0 && t>(0.9*To) && rfca.change_t!=0
    #         timeOffset = round(Int, t / To) * To
    #     end
    # end

    
    for ip in 1:np
        x -= rfca.dx
        y -= rfca.dy
    end

    if matrixMethod==0
        for ik in 0:nKicks-1
            for ip in 1:np
                if delta.map[1] == -1.0
                    continue
                end
                if length!=0
                    dc4 = length/2.0*sqrt(1.0+xp^2+yp^2)
                else
                    dc4 = CTPS(0.0, TPS_Dim, Max_TPS_Degree)
                end

                # energy kick
                P = P_central*(1.0 + delta)
                gamma = sqrt(P^2 + 1.0)
                beta_i = P/gamma
                t = (z+dc4)/(beta_i*2.99792458e8) - timeOffset
                if ik==0 && timeOffset!=0 && rfca.change_t!=0
                    z = t*beta_i*2.99792458e8 - dc4
                end
                dt = t-t0
                if dt.map[1] < 0
                    dt.map[1] = 0.0
                end
                if same_dgamma==0
                        if tau == 0
                            temp2 = CTPS(1.0, TPS_Dim, Max_TPS_Degree)
                        else
                            temp2 = sqrt(1.0 - exp(-dt / tau))
                        end
                        dgamma = sin(omega * (t - ik * dtLight) + phase) * temp2 * volt
                end
                if length != 0
                    if rfca.end1Focus!=0 && ik==0
                        # focus kick
                        inverseF = dgamma/(gamma*length*2.0)
                        xp -= x*inverseF
                        yp -= y*inverseF
                    else
                        inverseF = CTPS(0.0, TPS_Dim, Max_TPS_Degree)
                    end
                    # initial drift
                    x += xp*length/2.0
                    y += yp*length/2.0
                    z += length/2.0*sqrt(1.0+xp^2+yp^2)
                else
                    inverseF = CTPS(0.0, TPS_Dim, Max_TPS_Degree)
                end
                # energy kick
                xp,yp,z,delta = add_to_particle_energy(x,xp,y,yp,z,delta, t, P_central, dgamma)
                gamma1 = gamma+dgamma
                if gamma1.map[1]<=1
                    delta.map[1] = -1.0
                else
                    inverseF = -dgamma/(2.0*gamma1*length)
                end
            end

            if length!=0
                for ip in 1:np
                    x += xp*length/2.0
                    y += yp*length/2.0
                    z += length/2.0*sqrt(1.0+xp^2+yp^2)

                    if rfca.end2Focus!=0 && ik==nKicks-1
                        # focus kick
                        xp -= x*inverseF
                        yp -= y*inverseF
                    end
                end
            end
        end
    else
        # matrix method
        sin_phase = CTPS(0.0, TPS_Dim, Max_TPS_Degree)
        R11 = CTPS(1.0, TPS_Dim, Max_TPS_Degree)
        R21 = CTPS(0.0, TPS_Dim, Max_TPS_Degree)
        for ik in 0:nKicks-1
            for ip in 1:np
                # use matrix to propagate particles
                P = P_central*(1.0 + delta)
                gamma = sqrt(P^2 + 1.0)
                beta_i = P/gamma
                ds1 = length/2.0*sqrt(1.0+xp^2+yp^2)
                t = (z+ds1)/(beta_i*2.99792458e8) - timeOffset
                if timeOffset!=0 && rfca.change_t!=0
                    z = t*beta_i*2.99792458e8 - ds1
                end
                dt = t-t0
                if dt.map[1] < 0.0
                    dt.map[1] = 0.0
                end
                if same_dgamma==0
                    # if linearize ==0 # not use linearize
                        sin_phase = sin(omega * (t - ik * dtLight) + phase)
                        cos_phase = cos(omega * (t - ik * dtLight) + phase)
                        if tau == 0 
                            temp3 = CTPS(1.0, TPS_Dim, Max_TPS_Degree)
                        else
                            temp3 = sqrt(1.0 - exp(-dt / tau))
                        end
                        dgammaMax = volt * temp3
                        dgamma = dgammaMax * sin_phase
                end

                if rfca.end1Focus!=0 && ik==0 && length!=0
                    # focus kick
                    inverseF1 = dgamma/(2.0*gamma*length)
                    xp -= x*inverseF1
                    yp -= y*inverseF1
                end

                dP = sqrt((gamma+dgamma)^2 - 1.0) - P

                if useSRSModel != 0 
                    gammaf = gamma+dgamma
                    if abs(sin_phase.map[1]) > 1e-6
                        alpha = log(gammaf/gamma)/(2.0*sqrt(2.0)*sin_phase)
                    else
                        alpha = dgammaMax/gamma/(2.0*sqrt(2.0))
                    end
                    R11 = cos(alpha)
                    R22 = R11*gamma/gammaf
                    R12 = 2.0*sqrt(2.0)*gamma*length/dgammaMax*sin(alpha)
                    R21 = -sin(alpha)*dgammaMax/(length*gammaf*2*sqrt(2.0))
                else
                    R22 = 1.0/(1.0+dP/P)
                    if abs((dp/p).map[1]) > 1e-14
                        R12 = length*(P/dP*log(1.0+dP/P))
                    else
                        R12 = CTPS(length, TPS_Dim, Max_TPS_Degree)
                    end
                end

                z += ds1
                x_temp = CTPS(x)
                xp_temp = CTPS(xp)
                x = R11*x_temp + R12*xp_temp
                xp = R21*x_temp + R22*xp_temp

                x_temp = CTPS(y)
                xp_temp = CTPS(yp)
                y = R11*x_temp + R12*xp_temp
                yp = R21*x_temp + R22*xp_temp
                z += length/2.0*sqrt(1.0+xp^2+yp^2)
                delta = (P+dP-P_central)/P_central

                gamma += dgamma
                if gamma.map[1] <= 1
                    delta.map[1] = -1
                end
                if rfca.end2Focus!=0 && ik==nKicks-1 && length!=0
                    # focus kick
                    inverseF1 = -dgamma/(2.0*gamma*length)
                    xp -= x*inverseF1
                    yp -= y*inverseF1
                end
                z = (P+dP)/gamma*z/beta_i
            end
        end
    end

    for ip in 1:np
        x += rfca.dx
        y += rfca.dy
    end
    if rfca.change_p0 != 0
        x,xp,y,yp,z,delta,P_central = do_match_energy(x,xp,y,yp,z,delta, np, P_central, 0)
    end
    return x, xp, y, yp, z, delta#, P_central
    # return part, fiducial_seen, P_central, n_references, rfca_phase_reference, reference_ref_number, 
            #reference_phase, reference_flags
end


# function f(xx)
#     volt = xx[1]
#     freq = xx[2]
#     rfca = RFCA(name="rfca", freq=freq, volt=volt, phase=90.0, len=0.1, nSlices=10)
#     # parts = [0.001 0.0001 0.0005 0.0002 0.0 0.0; 0.001 0.0 0.0 0.0 0.0 0.0]
#     x = CTPS(0.0, 1, 6, 3)
#     xp = CTPS(0.0, 2, 6, 3)
#     y = CTPS(0.0, 3, 6, 3)
#     yp = CTPS(0.0, 4, 6, 3)
#     z = CTPS(0.0, 5, 6, 3)
#     delta = CTPS(0.0, 6, 6, 3)
#     Po = 19569.50762296901
#     x, xp, y, yp, z, delta = pass_TPSA(x, xp, y, yp, z, delta, rfca, Po, 0.0)
#     # println((1+parts[1,6])*Po)
#     return x.map[3]
# end
# volt = 1e6
# freq = 500e6
# println(f([volt, freq]))
# using Enzyme
# using BenchmarkTools
# @btime grad = f([volt, freq])
# println(grad)
# 
# @btime begin
#     # f([volt, freq])
#     # grad = jacobian(Reverse, f, [volt, freq], Val(6))
# end
