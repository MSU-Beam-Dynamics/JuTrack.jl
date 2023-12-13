using Zygote

function qfringe_R_matrix(dk_dz, l)
    if l==0
        return 1.0, 0.0, 0.0, 1.0
    end
    if dk_dz == 0
        return 1.0, l, 0.0, 1.0
    end
    l3 = l^3

    R11 = 0
    R21 = 0
    term = 1
    n = 0
    while abs(term / R11) > 1e-16
        R11 += term
        R21 += n * term
        term *= -l3 * dk_dz / ((n + 3) * (n + 2))
        n += 3
    end
    R21 /= l

    R12 = 0
    R22 = 0
    term = l
    n = 1
    while abs(term / R12) > 1e-16
        R12 += term
        R22 += n * term
        term *= -l3 * dk_dz / ((n + 3) * (n + 2))
        n += 3
    end
    R22 /= l
    return R11, R12, R21, R22
end

function qfringe_T_matrix(dk_dz, l, reverse)
    if l==0 || dk_dz == 0
        return 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    end
    l3 = l^3
    T116 = T216 = 0
    term = 1
    n = 0
    while abs(term / (T116 != 0 ? T116 : 1)) > 1e-16
        T116 += -n * term / 3
        T216 += -(n^2) * term / 3
        term *= -l3 * dk_dz / ((n + 3) * (n + 2))
        n += 3
    end
    T216 /= l

    T126 = T226 = 0
    term = l
    n = 1
    while abs(term / (T126 != 0 ? T126 : 1)) > 1e-16
        T126 += -n * term / 3
        T226 += -n*(n-1) * term / 3
        term *= -l3 * dk_dz / ((n + 3) * (n + 2))
        n += 3
    end
    T226 /= l
    if reverse != 0
        T116, T226 = T226, T116
    end
    if reverse == 0
        ko = dk_dz * l
        term = ko^2*l^3
        mult = ko*l^2
        T511 = 1/20*term
        term *= mult
        T511 += -1/240*term
        term *= mult
        T511 += 13/79200*term
        term *= mult
        T511 += -19/4989600*term
        term *= mult
        T511 += 19/325721088*term

        T522 = l
        term = l
        T522 += -1/6*term
        term *= mult
        T522 += 5/252*term
        term *= mult
        T522 += -11/11340*term
        term *= mult
        T522 += 187/7076160*term
        term *= mult
        T522 += -391/849139200*term
        T522 /= 2

        term = mult
        T512 = -1/6*term
        term *= mult
        T512 += 1/30*term
        term *= mult
        T512 += -1/480*term
        term *= mult
        T512 += 1/14784*term
        term *= mult
        T512 += -1/739200*term
        term *= mult
        T512 += 1/54454400*term
    else
        ko = dk_dz * l

        # T511 terms
        term = ko^2 * l^3
        mult = ko * l^2
        T511 = 2 / 15 * term; term *= mult
        T511 += -1 / 80 * term; term *= mult
        T511 += 1 / 1848 * term; term *= mult
        T511 += -1 / 73920 * term; term *= mult
        T511 += 3 / 13613600 * term

        # T522 terms
        T522 = term = l; term *= mult
        T522 += -1 / 6 * term; term *= mult
        T522 += 1 / 70 * term; term *= mult
        T522 += -1 / 1680 * term; term *= mult
        T522 += 1 / 68640 * term; term *= mult
        T522 += -3 / 12812800 * term
        T522 /= 2

        # T512 terms
        term = mult
        T512 = -1 / 3 * term; term *= mult
        T512 += 1 / 20 * term; term *= mult
        T512 += -1 / 336 * term; term *= mult
        T512 += 1 / 10560 * term; term *= mult
        T512 += -3 / 1601600 * term; term *= mult
        T512 += 1 / 39603200 * term
    end
    return T116, T126, T216, T226, T511, T512, T522
end

function quad_fringe(l, ko, order, reverse, fse)
    ko *= (1 + fse)
    C = [0.0 0.0 0.0 0.0 l 0.0]
    R = zeros(Float64, 6, 6)
    T = zeros(Float64, 6, 6, 6)
    R_Buffer = Zygote.Buffer(R)
    T_Buffer = Zygote.Buffer(T)
    for i in 1:6
        for j in 1:6
            R_Buffer[i, j] = 0.0
        end
    end
    for i in 1:6
        for j in 1:6
            for k in 1:6
                T_Buffer[i, j, k] = 0.0
            end
        end
    end
    R_Buffer[5, 5] = 1.0
    R_Buffer[6, 6] = 1.0
    R11, R12, R21, R22 = qfringe_R_matrix(ko/l, l)
    R33, R34, R43, R44 = qfringe_R_matrix(-ko/l, l)
    R_Buffer[1, 1] = R11
    R_Buffer[1, 2] = R12
    R_Buffer[2, 1] = R21
    R_Buffer[2, 2] = R22
    R_Buffer[3, 3] = R33
    R_Buffer[3, 4] = R34
    R_Buffer[4, 3] = R43
    R_Buffer[4, 4] = R44
    if reverse != 0
        R_Buffer[1,1] = R22
        R_Buffer[2,2] = R11
        R_Buffer[3,3] = R44
        R_Buffer[4,4] = R33
    end
    if order >= 2
        T116, T126, T216, T226, T511, T512, T522 = qfringe_T_matrix(ko/l, l, reverse)
        T_Buffer[1, 6, 1] = T116
        T_Buffer[1, 6, 2] = T126
        T_Buffer[2, 6, 1] = T216
        T_Buffer[2, 6, 2] = T226
        T_Buffer[5, 1, 1] = T511
        T_Buffer[5, 2, 1] = T512
        T_Buffer[5, 2, 2] = T522

        T363, T364, T463, T464, T533, T543, T433 = qfringe_T_matrix(-ko/l, l, reverse)
        T_Buffer[3, 6, 3] = T363
        T_Buffer[3, 6, 4] = T364
        T_Buffer[4, 6, 3] = T463
        T_Buffer[4, 6, 4] = T464
        T_Buffer[5, 3, 3] = T533
        T_Buffer[5, 4, 3] = T543
        T_Buffer[5, 4, 4] = T433
    end
    return C, copy(R), copy(T)

end
function initialize_matrices(order::Int)

    if order == 3
        C = zeros(Float64, 6)
        R = zeros(Float64, 6, 6)
        T = zeros(Float64, 6, 6, 6)  # Adjust the dimensions as required
        Q = [zeros(Float64, min(j, k) + 1) for j in 1:6, k in 1:6]
    elseif order == 2
        C = zeros(Float64, 6)
        R = zeros(Float64, 6, 6)
        T = zeros(Float64, 6, 6, 6)  
        Q = nothing
    elseif order == 1
        C = zeros(Float64, 6)
        R = zeros(Float64, 6, 6)
        T = nothing
        Q = nothing
    else
        error("Invalid order: $order (initialize_matrices)")
    end
    return C, R, T, Q
end

function drift_matrix(length::Float64, order::Int)
    # C, R, T, Q = initialize_matrices(M.order)
    C = [0 0 0 0 length 0]
    R = [1.0 length 0.0 0.0 0.0 0.0;
         0.0 1.0 0.0 0.0 0.0 0.0;
         0.0 0.0 1.0 length 0.0 0.0;
         0.0 0.0 0.0 1.0 0.0 0.0;
         0.0 0.0 0.0 0.0 1.0 0.0;
         0.0 0.0 0.0 0.0 0.0 1.0]

    if order > 1
        T[5, 2, 2] = T[5, 4, 4] = length / 2
    end
    if order > 1
        T = [i == 5 && ((j == 2 && k == 2) || (j == 4 && k == 4)) ? length / 2 : 0.0 
            for i in 1:6, j in 1:6, k in 1:6]
    else
        T = zeros(Float64, 6, 6, 6)
    end
    return C, R, T
end

function quadrupole_matrix(k1, lHC, maximum_order, tilt, fse, xkick, ykick, 
                            edge1_effects, edge2_effects, fringeType, ffringe, fringeIntM, fringeIntP)
    # lHC is the "hard core" length (wherein K1 is constant)
    # lNominal is the effective length.
    # If fringe effects are off, these are the same.
    if maximum_order != 1 && maximum_order != 2 && maximum_order != 3
        error("Invalid maximum order: $maximum_order (quadrupole_matrix)")
    end
    lEdge = 0
    k1 *= (1 + fse)

    if fringeType != "none" || fringeType != "inset" || fringeType != "fixed-strength" || fringeType != "integrals"
        error("Unknown fringe type for QUAD")
    end

    if fringeType == "none"
        fringeCode = -1
    elseif fringeType == "inset"
        fringeCode = 0
    elseif fringeType == "fixed-strength"
        fringeCode = 1
    elseif fringeType == "integrals"
        fringeCode = 2
    end
    
    if k1 == 0 || lHC == 0
        M = drift_matrix(lHC, maximum_order)
    else
        R = zeros(Float64, 6, 6)
        R_Buffer = Zygote.Buffer(R)
        for i in 1:6
            for j in 1:6
                R_Buffer[i, j] = 0.0
            end
        end 
        lNominal = lHC
        if ffringe != 0 && fringeCode != 2
            if edge1_effects == 0 || edge2_effects == 0
                error("Edge effects must not be 0 for fringe effects")
            end
            lEdge = lNominal * ffringe/2
            if fringeCode == 1
                lHc = lNominal - lEdge
            elseif fringeCode == 0
                lHc = lNominal - 2 * lEdge
            end
        end

        k = sqrt(abs(k1))
        kl = k * lHC
        sin_kl = sin(kl)
        cos_kl = cos(kl)
        sinh_kl = sinh(kl)
        cosh_kl = cosh(kl)

        R_Buffer[5,5] = 1
        R_Buffer[6,6] = 1
        C = [0 0 0 lHC 0 0]
        if k1 > 0
            # focussing in horizontal plane
            R_Buffer[1,1] = cos_kl
            R_Buffer[1,2] = sin_kl / k
            R_Buffer[2,1] = -k * sin_kl
            R_Buffer[2,2] = cos_kl
            R_Buffer[3,3] = cosh_kl
            R_Buffer[3,4] = sinh_kl / k
            R_Buffer[4,3] = k * sinh_kl
            R_Buffer[4,4] = cosh_kl

            if maximum_order >= 2
                k2 = k^2
                k3 = k^3
                k4 = k^4
                l2 = lHC^2
                l3 = lHC^3
                l4 = lHC^4
                cos_kl = cos(k*lHC)
                cos_2kl = cos(2*k*lHC)
                cos_3kl = cos(3*k*lHC)
                sin_kl = sin(k*lHC)
                sin_2kl = sin(2*k*lHC)
                sin_3kl = sin(3*k*lHC)
                cosh_kl = cosh(k*lHC)
                cosh_2kl = cosh(2*k*lHC)
                cosh_3kl = cosh(3*k*lHC)
                sinh_kl = sinh(k*lHC)
                sinh_2kl = sinh(2*k*lHC)
                sinh_3kl = sinh(3*k*lHC)

                T = zeros(Float64, 6, 6, 6)
                T_Buffer = Zygote.Buffer(T)
                for i in 1:6
                    for j in 1:6
                        for k in 1:6
                            T_Buffer[i, j, k] = 0.0
                        end
                    end
                end
                T_Buffer[1, 6, 1] = kl * sin_kl / 2
                T_Buffer[2, 6, 2] = kl * sin_kl / 2
                T_Buffer[1, 6, 2] = sin_kl / (2 * k) - lHC * cos_kl / 2
                T_Buffer[2, 6, 1] = k/2*(kl*cos_kl + sin_kl)
                T_Buffer[3, 6, 3] = -kl/2*sinh_kl
                T_Buffer[3, 6, 4] = (sinh_kl/k - lHC*cosh_kl)/2
                T_Buffer[4, 6, 3] = -k / 2 * (kl * cosh_kl + sinh_kl)
                T_Buffer[4, 6, 4] = -kl / 2 * sinh_kl
                T_Buffer[5, 1, 1] = (k^2) * (lHC - sin_kl / k * cos_kl) / 4
                T_Buffer[5, 2, 1] = -(sin_kl^2) / 2
                T_Buffer[5, 2, 2] = (lHC + sin_kl / k * cos_kl) / 4
                T_Buffer[5, 3, 3] = -(k^2) * (lHC - sinh_kl / k * cosh_kl) / 4
                T_Buffer[5, 4, 3] = (sinh_kl^2) / 2
                T_Buffer[5, 4, 4] = (lHC + sinh_kl / k * cosh_kl) / 4
                if maximum_order >= 3
                    U = zeros(Float64, 6, 6, 6, 6)
                    U_Buffer = Zygote.Buffer(U)
                    for i in 1:6
                        for j in 1:6
                            for k in 1:6
                                for l in 1:6
                                    U_Buffer[i, j, k, l] = 0.0
                                end
                            end
                        end
                    end
                    U_Buffer[1, 1, 1, 1] = (-3 * k3 * lHC * sin_kl) / 16 + (3 * k2 * sin_kl * sin_2kl) / 32
                    U_Buffer[1, 2, 1, 1] = (-3 * k2 * lHC * cos_kl) / 16 + (21 * k * sin_kl) / 64 - (3 * k * sin_3kl) / 64
                    U_Buffer[1, 2, 2, 1] = (-9 * k * lHC * sin_kl) / 16 - (3 * sin_kl * sin_2kl) / 32
                    U_Buffer[1, 2, 2, 2] = (3 * lHC * cos_kl) / 16 - (21 * sin_kl) / (64 * k) + (3 * sin_3kl) / (64 * k)
                    U_Buffer[1, 3, 3, 1] = (k3*lHC*sin_kl)/8 + (k2*cos_kl*(sinh_kl^2))/16 - (3*k2*sin_kl*sinh_2kl)/32
                    U_Buffer[1, 3, 3, 2] = -(k2*lHC*cos_kl)/8 - (3*k*sin_kl)/32 + (k*cosh_2kl*sin_kl)/32 + (3*k*cos_kl*sinh_2kl)/32
                    U_Buffer[1, 4, 3, 1] = (k2*lHC*cos_kl)/4 - (7*k*sin_kl)/32 - (3*k*cosh_2kl*sin_kl)/32 + (k*cos_kl*sinh_2kl)/32
                    U_Buffer[1, 4, 3, 2] = (k*lHC*sin_kl)/4 + (3*cos_kl*(sinh_kl^2))/16 + (sin_kl*sinh_2kl)/32
                    U_Buffer[1, 4, 4, 1] = -(k*lHC*sin_kl)/8 + (cos_kl*(sinh_kl^2))/16 - (3*sin_kl*sinh_2kl)/32
                    U_Buffer[1, 4, 4, 2] = (lHC*cos_kl)/8 - (11*sin_kl)/(32*k) + (cosh_2kl*sin_kl)/(32*k) + (3*cos_kl*sinh_2kl)/(32*k)
                    U_Buffer[1, 6, 6, 1] = -(k2*l2*cos_kl)/8 - (3*k*lHC*sin_kl)/8
                    U_Buffer[1, 6, 6, 2] = (lHC*cos_kl)/8 - sin_kl/(8*k) - (k*l2*sin_kl)/8
                    U_Buffer[2, 1, 1, 1] = (-3*k4*lHC*cos_kl)/16 - (3*k3*sin_kl)/16 + (3*k3*cos_2kl*sin_kl)/16 + (3*k3*cos_kl*sin_2kl)/32
                    U_Buffer[2, 2, 1, 1] = (9*k2*cos_kl)/64 - (9*k2*cos_3kl)/64 + (3*k3*lHC*sin_kl)/16
                    U_Buffer[2, 2, 2, 1] = (-9*k2*lHC*cos_kl)/16 - (9*k*sin_kl)/16 - (3*k*cos_2kl*sin_kl)/16 - (3*k*cos_kl*sin_2kl)/32
                    U_Buffer[2, 2, 2, 2] = (-9*cos_kl)/64 + (9*cos_3kl)/64 - (3*k*lHC*sin_kl)/16
                    U_Buffer[2, 3, 3, 1] = (k4*lHC*cos_kl)/8 + (k3*sin_kl)/8 - (3*k3*cosh_2kl*sin_kl)/16 + (k3*cos_kl*cosh_kl*sinh_kl)/8 - (k3*sin_kl*(sinh_kl^2))/16 - (3*k3*cos_kl*sinh_2kl)/32
                    U_Buffer[2, 3, 3, 2] = (-7 * k2 * cos_kl) / 32 + (7 * k2 * cos_kl * cosh_2kl) / 32 + (k3 * lHC * sin_kl) / 8 - (k2 * sin_kl * sinh_2kl) / 32
                    U_Buffer[2, 4, 3, 1] = (k2 * cos_kl) / 32 - (k2 * cos_kl * cosh_2kl) / 32 - (k3 * lHC * sin_kl) / 4 - (7 * k2 * sin_kl * sinh_2kl) / 32
                    U_Buffer[2, 4, 3, 2] = (k2 * lHC * cos_kl) / 4 + (k * sin_kl) / 4 + (k * cosh_2kl * sin_kl) / 16 + (3 * k * cos_kl * cosh_kl * sinh_kl) / 8 - (3 * k * sin_kl * (sinh_kl^2)) / 16 + (k * cos_kl * sinh_2kl) / 32
                    U_Buffer[2, 4, 4, 1] = -(k2 * lHC * cos_kl) / 8 - (k * sin_kl) / 8 - (3 * k * cosh_2kl * sin_kl) / 16 + (k * cos_kl * cosh_kl * sinh_kl) / 8 - (k * sin_kl * (sinh_kl^2)) / 16 - (3 * k * cos_kl * sinh_2kl) / 32
                    U_Buffer[2, 4, 4, 2] = (-7 * cos_kl) / 32 + (7 * cos_kl * cosh_2kl) / 32 - (k * lHC * sin_kl) / 8 - (sin_kl * sinh_2kl) / 32
                    U_Buffer[2, 6, 6, 2] = (5 * k2 * lHC * cos_kl) / 8 + (3 * k * sin_kl) / 8 + (k3 * l2 * sin_kl) / 8
                    U_Buffer[2, 6, 6, 3] = -(k2 * l2 * cos_kl) / 8 - (3 * k * lHC * sin_kl) / 8
                    U_Buffer[3, 3, 1, 1] = (k2 * cosh_kl * (sin_kl^2)) / 16 + (k3 * lHC * sinh_kl) / 8 - (3 * k2 * sin_2kl * sinh_kl) / 32
                    U_Buffer[3, 3, 2, 1] = -(k2 * lHC * cosh_kl) / 4 - (k * cosh_kl * sin_2kl) / 32 + (7 * k * sinh_kl) / 32 + (3 * k * cos_2kl * sinh_kl) / 32
                    U_Buffer[3, 3, 2, 2] = -(cosh_kl * (sin_kl^2)) / 16 + (k * lHC * sinh_kl) / 8 + (3 * sin_2kl * sinh_kl) / 32
                    U_Buffer[3, 3, 3, 3] = (-3 * k3 * lHC * sinh_kl) / 16 + (3 * k2 * sinh_kl * sinh_2kl) / 32
                    U_Buffer[3, 4, 1, 1] = (k2 * lHC * cosh_kl) / 8 - (3 * k * cosh_kl * sin_2kl) / 32 + (3 * k * sinh_kl) / 32 - (k * cos_2kl * sinh_kl) / 32
                    U_Buffer[3, 4, 2, 1] = (-3 * cosh_kl) / 32 + (3 * cos_2kl * cosh_kl) / 32 - (k * lHC * sinh_kl) / 4 - (sin_2kl * sinh_kl) / 32
                    U_Buffer[3, 4, 2, 2] = (lHC * cosh_kl) / 8 + (3 * cosh_kl * sin_2kl) / (32 * k) - (11 * sinh_kl) / (32 * k) + (cos_2kl * sinh_kl) / (32 * k)
                    U_Buffer[3, 4, 3, 3] = (3 * k2 * lHC * cosh_kl) / 16 - (21 * k * sinh_kl) / 64 + (3 * k * sinh_3kl) / 64
                    U_Buffer[3, 4, 4, 3] = (9 * k * lHC * sinh_kl) / 16 + (3 * sinh_kl * sinh_2kl) / 32
                    U_Buffer[3, 4, 4, 4] = (3 * lHC * cosh_kl) / 16 - (21 * sinh_kl) / (64 * k) + (3 * sinh_3kl) / (64 * k)
                    U_Buffer[3, 6, 6, 3] = (k2 * l2 * cosh_kl) / 8 + (3 * k * lHC * sinh_kl) / 8
                    U_Buffer[3, 6, 6, 4] = (lHC * cosh_kl) / 8 - sinh_kl / (8 * k) + (k * l2 * sinh_kl) / 8
                    U_Buffer[4, 3, 1, 1] = (k4 * lHC * cosh_kl) / 8 + (k3 * cos_kl * cosh_kl * sin_kl) / 8 - (3 * k3 * cosh_kl * sin_2kl) / 32 + (k3 * sinh_kl) / 8 - (3 * k3 * cos_2kl * sinh_kl) / 16 + (k3 * (sin_kl^2) * sinh_kl) / 16
                    U_Buffer[4, 3, 2, 1] = -(k2 * cosh_kl) / 32 + (k2 * cos_2kl * cosh_kl) / 32 - (k3 * lHC * sinh_kl) / 4 - (7 * k2 * sin_2kl * sinh_kl) / 32
                    U_Buffer[4, 3, 2, 2] = (k2 * lHC * cosh_kl) / 8 - (k * cos_kl * cosh_kl * sin_kl) / 8 + (3 * k * cosh_kl * sin_2kl) / 32 + (k * sinh_kl) / 8 + (3 * k * cos_2kl * sinh_kl) / 16 - (k * (sin_kl^2) * sinh_kl) / 16
                    U_Buffer[4, 3, 3, 3] = (-3 * k4 * lHC * cosh_kl) / 16 - (3 * k3 * sinh_kl) / 16 + (3 * k3 * cosh_2kl * sinh_kl) / 16 + (3 * k3 * cosh_kl * sinh_2kl) / 32
                    U_Buffer[4, 4, 1, 1] = (7 * k2 * cosh_kl) / 32 - (7 * k2 * cos_2kl * cosh_kl) / 32 + (k3 * lHC * sinh_kl) / 8 - (k2 * sin_2kl * sinh_kl) / 32
                    U_Buffer[4, 4, 2, 1] = -(k2 * lHC * cosh_kl) / 4 - (7 * k * cosh_kl * sin_2kl) / 32 - (11 * k * sinh_kl) / 32 + (k * cos_2kl * sinh_kl) / 32
                    U_Buffer[4, 4, 2, 2] = (-7*cosh_kl)/32 + (7*cos_2kl*cosh_kl)/32 + (k*lHC*sinh_kl)/8 + (sin_2kl*sinh_kl)/32
                    U_Buffer[4, 4, 3, 3] = (-9*k2*cosh_kl)/64 + (9*k2*cosh_3kl)/64 + (3*k3*lHC*sinh_kl)/16
                    U_Buffer[4, 4, 4, 3] = (9*k2*lHC*cosh_kl)/16 + (9*k*sinh_kl)/16 + (3*k*cosh_2kl*sinh_kl)/16 + (3*k*cosh_kl*sinh_2kl)/32
                    U_Buffer[4, 4, 4, 4] = (-9*cosh_kl)/64 + (9*cosh_3kl)/64 + (3*k*lHC*sinh_kl)/16
                    U_Buffer[4, 6, 6, 3] = (5*k2*lHC*cosh_kl)/8 + (3*k*sinh_kl)/8 + (k3*l2*sinh_kl)/8
                    U_Buffer[4, 6, 6, 4] = (k2*l2*cosh_kl)/8 + (3*k*lHC*sinh_kl)/8
                end
            end
        else
            R_Buffer[3,3] = cos_kl
            R_Buffer[3,4] = sin_kl / k
            R_Buffer[4,3] = -k * sin_kl
            R_Buffer[4,4] = cos_kl
            R_Buffer[1,1] = cosh_kl
            R_Buffer[1,2] = sinh_kl / k
            R_Buffer[2,1] = k * sinh_kl
            R_Buffer[2,2] = cosh_kl

            if maximum_order >= 2
                k2 = k^2
                k3 = k^3
                k4 = k^4
                l2 = lHC^2
                l3 = lHC^3
                l4 = lHC^4
                cos_kl = cos(k*lHC)
                cos_2kl = cos(2*k*lHC)
                cos_3kl = cos(3*k*lHC)
                sin_kl = sin(k*lHC)
                sin_2kl = sin(2*k*lHC)
                sin_3kl = sin(3*k*lHC)
                cosh_kl = cosh(k*lHC)
                cosh_2kl = cosh(2*k*lHC)
                cosh_3kl = cosh(3*k*lHC)
                sinh_kl = sinh(k*lHC)
                sinh_2kl = sinh(2*k*lHC)
                sinh_3kl = sinh(3*k*lHC)

                T = zeros(Float64, 6, 6, 6)
                T_Buffer = Zygote.Buffer(T)
                for i in 1:6
                    for j in 1:6
                        for k in 1:6
                            T_Buffer[i, j, k] = 0.0
                        end
                    end
                end
                T_Buffer[3, 6, 3] = T_Buffer[4, 6, 4] = kl * sin_kl / 2
                T_Buffer[3, 6, 4] = sin_kl / (2 * k) - lHC * cos_kl / 2
                T_Buffer[4, 6, 3] = k / 2 * (kl * cos_kl + sin_kl)
                T_Buffer[1, 6, 1] = T_Buffer[2, 6, 2] = -kl / 2 * sinh_kl
                T_Buffer[1, 6, 2] = (sinh_kl / k - lHC * cosh_kl) / 2
                T_Buffer[2, 6, 1] = -k / 2 * (kl * cosh_kl + sinh_kl)
                T_Buffer[5, 1, 1] = -k^2 * (lHC - sinh_kl / k * cosh_kl) / 4
                T_Buffer[5, 2, 1] = sinh_kl^2 / 2
                T_Buffer[5, 2, 2] = (lHC + sinh_kl / k * cosh_kl) / 4
                T_Buffer[5, 3, 3] = k^2 * (lHC - sin_kl / k * cos_kl) / 4
                T_Buffer[5, 4, 3] = -sin_kl^2 / 2
                T_Buffer[5, 4, 4] = (lHC + sin_kl / k * cos_kl) / 4

                if maximum_order >= 3
                    U = zeros(Float64, 6, 6, 6, 6)
                    U_Buffer = Zygote.Buffer(U)
                    for i in 1:6
                        for j in 1:6
                            for k in 1:6
                                for l in 1:6
                                    U_Buffer[i, j, k, l] = 0.0
                                end
                            end
                        end
                    end
                    U_Buffer[1, 1, 1, 1] = (-3 * k3 * lHC * sinh_kl) / 16 + (3 * k2 * sinh_kl * sinh_2kl) / 32
                    U_Buffer[1, 2, 1, 1] = (3 * k2 * lHC * cosh_kl) / 16 - (21 * k * sinh_kl) / 64 + (3 * k * sinh_3kl) / 64
                    U_Buffer[1, 2, 2, 1] = (9 * k * lHC * sinh_kl) / 16 + (3 * sinh_kl * sinh_2kl) / 32
                    U_Buffer[1, 2, 2, 2] = (3 * lHC * cosh_kl) / 16 - (21 * sinh_kl) / (64 * k) + (3 * sinh_3kl) / (64 * k)
                    U_Buffer[1, 3, 3, 1] = (k2 * cosh_kl * sin_kl^2) / 16 + (k3 * lHC * sinh_kl) / 8 - (3 * k2 * sin_2kl * sinh_kl) / 32
                    U_Buffer[1, 3, 3, 2] = (k2 * lHC * cosh_kl) / 8 - (3 * k * cosh_kl * sin_2kl) / 32 + (3 * k * sinh_kl) / 32 - (k * cos_2kl * sinh_kl) / 32
                    U_Buffer[1, 4, 3, 1] = -(k2 * lHC * cosh_kl) / 4 - (k * cosh_kl * sin_2kl) / 32 + (7 * k * sinh_kl) / 32 + (3 * k * cos_2kl * sinh_kl) / 32
                    U_Buffer[1, 4, 3, 2] = (-3 * cosh_kl) / 32 + (3 * cos_2kl * cosh_kl) / 32 - (k * lHC * sinh_kl) / 4 - (sin_2kl * sinh_kl) / 32
                    U_Buffer[1, 4, 4, 1] = -(cosh_kl * sin_kl^2) / 16 + (k * lHC * sinh_kl) / 8 + (3 * sin_2kl * sinh_kl) / 32
                    U_Buffer[1, 4, 4, 2] = (lHC * cosh_kl) / 8 + (3 * cosh_kl * sin_2kl) / (32 * k) - (11 * sinh_kl) / (32 * k) + (cos_2kl * sinh_kl) / (32 * k)
                    U_Buffer[1, 6, 6, 1] = (k2 * l2 * cosh_kl) / 8 + (3 * k * lHC * sinh_kl) / 8
                    U_Buffer[1, 6, 6, 2] = (lHC * cosh_kl) / 8 - sinh_kl / (8 * k) + (k * l2 * sinh_kl) / 8
                    U_Buffer[2, 1, 1, 1] = (-3 * k4 * lHC * cosh_kl) / 16 - (3 * k3 * sinh_kl) / 16 + (3 * k3 * cosh_2kl * sinh_kl) / 16 + (3 * k3 * cosh_kl * sinh_2kl) / 32
                    U_Buffer[2, 2, 1, 1] = (-9 * k2 * cosh_kl) / 64 + (9 * k2 * cosh_3kl) / 64 + (3 * k3 * lHC * sinh_kl) / 16
                    U_Buffer[2, 2, 2, 1] = (9 * k2 * lHC * cosh_kl) / 16 + (9 * k * sinh_kl) / 16 + (3 * k * cosh_2kl * sinh_kl) / 16 + (3 * k * cosh_kl * sinh_2kl) / 32
                    U_Buffer[2, 2, 2, 2] = (-9 * cosh_kl) / 64 + (9 * cosh_3kl) / 64 + (3 * k * lHC * sinh_kl) / 16
                    U_Buffer[2, 3, 3, 1] = (k4 * lHC * cosh_kl) / 8 + (k3 * sin_kl) / 8 - (3 * k3 * cosh_2kl * sin_kl) / 16 + (k3 * cosh_kl * cosh_kl * sinh_kl) / 8 - (k3 * sin_kl * (sinh_kl^2)) / 16 - (3 * k3 * cosh_kl * sinh_2kl) / 32
                    U_Buffer[2, 3, 3, 2] = (7 * k2 * cosh_kl) / 32 - (7 * k2 * cos_2kl * cosh_kl) / 32 + (k3 * lHC * sinh_kl) / 8 - (k2 * sin_2kl * sinh_kl) / 32
                    U_Buffer[2, 4, 3, 1] = -(k2 * cosh_kl) / 32 + (k2 * cos_2kl * cosh_kl) / 32 - (k3 * lHC * sinh_kl) / 4 - (7 * k2 * sin_2kl * sinh_kl) / 32
                    U_Buffer[2, 4, 3, 2] = -(k2 * lHC * cosh_kl) / 4 - (7 * k * cosh_kl * sin_2kl) / 32 - (11 * k * sinh_kl) / 32 + (k * cos_2kl * sinh_kl) / 32
                    U_Buffer[2, 4, 4, 1] = (k2 * lHC * cosh_kl) / 8 - (k * sin_kl * cosh_kl * sin_kl) / 8 + (3 * k * cosh_kl * sin_2kl) / 32 + (k * sinh_kl) / 8 + (3 * k * cos_2kl * sinh_kl) / 16 - (k * (sin_kl^2) * sinh_kl) / 16
                    U_Buffer[2, 4, 4, 2] = (-7 * cosh_kl) / 32 + (7 * cos_2kl * cosh_kl) / 32 + (k * lHC * sinh_kl) / 8 + (sin_2kl * sinh_kl) / 32
                    U_Buffer[2, 6, 6, 1] = (5 * k2 * lHC * cosh_kl) / 8 + (3 * k * sinh_kl) / 8 + (k3 * l2 * sinh_kl) / 8
                    U_Buffer[2, 6, 6, 2] = (k2 * l2 * cosh_kl) / 8 - (3 * k * lHC * sinh_kl) / 8

                    U_Buffer[3, 3, 1, 1] = (k3 * lHC * sin_kl) / 8 + (k2 * cos_kl * sinh_kl^2) / 16 - (3 * k2 * sin_kl * sinh_2kl) / 32
                    U_Buffer[3, 3, 2, 1] = (k2 * lHC * cos_kl) / 4 - (7 * k * sin_kl) / 32 - (3 * k * cosh_2kl * sin_kl) / 32 + (k * cos_kl * sinh_2kl) / 32
                    U_Buffer[3, 3, 2, 2] = -(k * lHC * sin_kl) / 8 + (cos_kl * sinh_kl^2) / 16 - (3 * sin_kl * sinh_2kl) / 32
                    U_Buffer[3, 3, 3, 3] = (-3 * k3 * lHC * sin_kl) / 16 + (3 * k2 * sin_kl * sin_2kl) / 32
                    U_Buffer[3, 4, 1, 1] = -(k2 * lHC * cos_kl) / 8 - (3 * k * sin_kl) / 32 + (k * cosh_2kl * sin_kl) / 32 + (3 * k * cos_kl * sinh_2kl) / 32
                    U_Buffer[3, 4, 2, 1] = (k * lHC * sin_kl) / 4 + (3 * cos_kl * sinh_kl^2) / 16 + (sin_kl * sinh_2kl) / 32
                    U_Buffer[3, 4, 2, 2] = (lHC * cos_kl) / 8 - (11 * sin_kl) / (32 * k) + (cosh_2kl * sin_kl) / (32 * k) + (3 * cos_kl * sinh_2kl) / (32 * k)
                    U_Buffer[3, 4, 3, 3] = (-3 * k2 * lHC * cos_kl) / 16 + (21 * k * sin_kl) / 64 - (3 * k * sin_3kl) / 64
                    U_Buffer[3, 4, 4, 3] = (-9 * k * lHC * sin_kl) / 16 - (3 * sin_kl * sin_2kl) / 32
                    U_Buffer[3, 4, 4, 4] = (3 * lHC * cos_kl) / 16 - (21 * sin_kl) / (64 * k) + (3 * sin_3kl) / (64 * k)
                    U_Buffer[3, 6, 6, 3] = -(k2 * l2 * cos_kl) / 8 - (3 * k * lHC * sin_kl) / 8
                    U_Buffer[3, 6, 6, 4] = (lHC * cos_kl) / 8 - sin_kl / (8 * k) - (k * l2 * sin_kl) / 8
                    U_Buffer[4, 3, 1, 1] = (k4 * lHC * cos_kl) / 8 + (k3 * sin_kl) / 8 - (3 * k3 * cosh_2kl * sin_kl) / 16 + (k3 * cos_kl * cosh_kl * sinh_kl) / 8 - (k3 * sin_kl * sinh_kl^2) / 16 - (3 * k3 * cos_kl * sinh_2kl) / 32
                    U_Buffer[4, 3, 2, 1] = (k2 * cos_kl) / 32 - (k2 * cos_kl * cosh_2kl) / 32 - (k3 * lHC * sin_kl) / 4 - (7 * k2 * sin_kl * sinh_2kl) / 32
                    U_Buffer[4, 3, 2, 2] = -(k2 * lHC * cos_kl) / 8 - (k * sin_kl) / 8 - (3 * k * cosh_2kl * sin_kl) / 16 + (k * cos_kl * cosh_kl * sinh_kl) / 8 - (k * sin_kl * sinh_kl^2) / 16 - (3 * k * cos_kl * sinh_2kl) / 32
                    U_Buffer[4, 3, 3, 3] = (-3 * k4 * lHC * cos_kl) / 16 - (3 * k3 * sin_kl) / 16 + (3 * k3 * cos_2kl * sin_kl) / 16 + (3 * k3 * cos_kl * sin_2kl) / 32
                    U_Buffer[4, 4, 1, 1] = (-7 * k2 * cos_kl) / 32 + (7 * k2 * cos_kl * cosh_2kl) / 32 + (k3 * lHC * sin_kl) / 8 - (k2 * sin_kl * sinh_2kl) / 32
                    U_Buffer[4, 4, 2, 1] = (k2 * lHC * cos_kl) / 4 + (k * sin_kl) / 4 + (k * cosh_2kl * sin_kl) / 16 + (3 * k * cos_kl * cosh_kl * sinh_kl) / 8 - (3 * k * sin_kl * sinh_kl^2) / 16 + (k * cos_kl * sinh_2kl) / 32
                    U_Buffer[4, 4, 2, 2] = (-7 * cos_kl) / 32 + (7 * cos_kl * cosh_2kl) / 32 - (k * lHC * sin_kl) / 8 - (sin_kl * sinh_2kl) / 32
                    U_Buffer[4, 4, 3, 3] = (9 * k2 * cos_kl) / 64 - (9 * k2 * cos_3kl) / 64 + (3 * k3 * lHC * sin_kl) / 16
                    U_Buffer[4, 4, 4, 4] = (-9 * k2 * lHC * cos_kl) / 16 - (9 * k * sin_kl) / 16 - (3 * k * cos_2kl * sin_kl) / 16 - (3 * k * cos_kl * sin_2kl) / 32
                    U_Buffer[4, 4, 4, 5] = (-9 * cos_kl) / 64 + (9 * cos_3kl) / 64 - (3 * k * lHC * sin_kl) / 16
                    U_Buffer[4, 6, 6, 3] = (-5 * k2 * lHC * cos_kl) / 8 - (3 * k * sin_kl) / 8 + (k3 * l2 * sin_kl) / 8
                    U_Buffer[4, 6, 6, 4] = -(k2 * l2 * cos_kl) / 8 - (3 * k * lHC * sin_kl) / 8
                end
            end
        end
        
        if fringeCode == 0 || fringeCode == 1
            if lEdge != 0 && k1 !=0
                Mfringe = quad_fringe(lEdge, k1, maximum_order, 0, 0.0)
            end
            if fringeCode == 1
                Md = drift_matrix(-lEdge/2, maximum_order)
                
    end
end