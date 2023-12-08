using Zygote
include("../lattice/canonical_elements.jl")

function exactDrift(part, np, length)
    return [
        [
            coord[1] + coord[2] * length,
            coord[2],
            coord[3] + coord[4] * length,
            coord[4],
            coord[5] + length * sqrt(1 + coord[2]^2 + coord[4]^2),
            coord[6]
        ]
        for coord in part
    ]
end


function computeCSBENDFieldCoefficients(b, h, nonlinear, expansionOrder)
    if expansionOrder == 0
        # set the order to be <highestMultipole> + 2
        i = length(b)
        while i >= 1
            if b[i] != 0
                break
            end
            i -= 1
        end
        expansionOrder = i + 2  
        if expansionOrder < 4
            # make minimum value 4 for backward compatibility
            expansionOrder = 4
        end
    end
    
    expansionOrder1 = expansionOrder + 1
    if (expansionOrder1>11)
        error("expansion order >10 for CSBEND or CSRCSBEND")
    end

    Fx_xy = zeros(11, 11)
    Fy_xy = zeros(11, 11)
    Fx_xy_Buffer = Zygote.Buffer(Fx_xy)
    Fy_xy_Buffer = Zygote.Buffer(Fy_xy)
    for i in 1:11
        for j in 1:11
            Fx_xy_Buffer[i, j] = 0
            Fy_xy_Buffer[i, j] = 0
        end
    end
    
    h2 = h^2
    h3 = h^3
    h4 = h^4
    h5 = h^5

    Fx_xy_Buffer[1,2] = b[1]
    Fy_xy_Buffer[1,1] = 1
    Fy_xy_Buffer[2,1] = b[1]

    if nonlinear !=0 && any(b .!= 0)
        Fx_xy_Buffer[2,2] = b[2]
        Fx_xy_Buffer[3,2] = b[3]/2
        Fx_xy_Buffer[4,2] = b[4]/6
        Fx_xy_Buffer[5,2] = b[5]/24
        Fx_xy_Buffer[6,2] = b[6]/120
        Fx_xy_Buffer[7,2] = b[7]/720
        Fx_xy_Buffer[8,2] = b[8]/5040
        Fx_xy_Buffer[1,4] = (-b[3] - b[2]*h + b[1]*h2)/6
        Fx_xy_Buffer[2,4] = (-b[4] - b[3]*h + 2*b[2]*h2 - 2*b[1]*h3)/6
        Fx_xy_Buffer[3,4] = (-b[5] + h*(-b[4] + 3*h*(b[3] - 2*b[2]*h + 2*b[1]*h2)))/12
        Fx_xy_Buffer[4,4] = (-b[6] - h*(b[5] + 4*h*(-b[4] + 3*h*(b[3] - 2*b[2]*h + 2*b[1]*h2))))/36
        Fx_xy_Buffer[5,4] = (-b[7] + h*(-b[6] + 5*h*(b[5] + 4*h*(-b[4] + 3*h*(b[3] - 2*b[2]*h + 2*b[1]*h2)))))/ 144
        Fx_xy_Buffer[6,4] = (-b[8] - h*(b[7] + 6*h*(-b[6] + 5*h*(b[5] + 4*h*(-b[4] + 3*h*(b[3] - 2*b[2]*h + 2*b[1]*h2))))))/720
        Fx_xy_Buffer[7,4] = (h*(-b[8] + 7*h*(b[7] + 6*h*(-b[6] + 5*h*(b[5] + 4*h*(-b[4] + 3*h*(b[3] - 2*b[2]*h + 2*b[1]*h2)))))))/ 4320
        Fx_xy_Buffer[8,4] = (h2*(b[8] - 7*h*(b[7] + 6*h*(-b[6] + 5*h*(b[5] + 4*h*(-b[4] + 3*h*(b[3] - 2*b[2]*h + 2*b[1]*h2)))))))/3780
        Fx_xy_Buffer[1,6] = (b[5] + h*(2*b[4] - 3*h*(b[3] - b[2]*h + b[1]*h2)))/120
        Fx_xy_Buffer[2,6] = (b[6] + h*(2*b[5] + h*(-5*b[4] + 9*b[3]*h - 12*b[2]*h2 + 12*b[1]*h3)))/ 120
        Fx_xy_Buffer[3,6] = (b[7] + h*(2*b[6] + h*(-7*b[5] + 19*b[4]*h - 39*b[3]*h2 + 60*b[2]*h3 - 60*b[1]*h4)))/240
        Fx_xy_Buffer[4,6] = (b[8] + h*(2*b[7] + 3*h*(-3*b[6] + 11*b[5]*h - 32*b[4]*h2 + 72*b[3]*h3 - 120*b[2]*h4 + 120*b[1]*h5)))/720
        Fx_xy_Buffer[5,6] = (h*(2*b[8] + h*(-11*b[7] - 3*h*(-17*b[6] + 5*h*(13*b[5] + 8*h*(-5*b[4] + 3*h*(4*b[3] - 7*b[2]*h + 7*b[1]*h2)))))))/2880
        Fx_xy_Buffer[6,6] = (h2*(-13*b[8] + h*(73*b[7] + 12*h*(-29*b[6] + 5*h*(23*b[5] + 2*h*(-37*b[4] + 3*h*(31*b[3] + 56*h*(-b[2] + b[1]*h))))))))/14400
        Fx_xy_Buffer[1,8] = (-b[7] - 3*h*(b[6] + h*(-2*b[5] + h*(4*b[4] - 9*b[3]*h + 15*b[2]*h2 - 15*b[1]*h3))))/5040
        Fx_xy_Buffer[2,8] = (-b[8] - 3*h*(b[7] + h*(-3*b[6] + 8*b[5]*h - 21*b[4]*h2 + 51*b[3]*h3 - 90*b[2]*h4 + 90*b[1]*h5)))/5040
        Fx_xy_Buffer[3,8] = (h*(-b[8] + h*(4*b[7] + h*(-14*b[6] + 45*b[5]*h - 135*b[4]*h2 + 345*b[3]*h3 - 630*b[2]*h4 + 630*b[1]*h5))))/ 3360
        Fx_xy_Buffer[4,8] = (h2*(5*b[8] + h*(-22*b[7] + 3*h*(29*b[6] - 5*h*(21*b[5] + 4*h*(-17*b[4] + 45*b[3]*h - 84*b[2]*h2 + 84*b[1]*h3))))))/10080
        Fx_xy_Buffer[1,10] = (h*(4*b[8] - 5*h*(2*b[7] + 3*h*(-2*b[6] + 7*b[5]*h - 22*b[4]*h2 + 57*b[3]*h3 - 105*b[2]*h4 + 105*b[1]*h5))))/ 362880
        Fx_xy_Buffer[2,10] = (h2*(-14*b[8] + 5*h*(10*b[7] + 3*h*(-13*b[6] + 50*b[5]*h - 167*b[4]*h2 + 447*b[3]*h3 - 840*b[2]*h4 + 840*b[1]*h5))))/362880

        Fy_xy_Buffer[3,1] = b[2] / 2
        Fy_xy_Buffer[4,1] = b[3] / 6
        Fy_xy_Buffer[5,1] = b[4] / 24
        Fy_xy_Buffer[6,1] = b[5] / 120
        Fy_xy_Buffer[7,1] = b[6] / 720
        Fy_xy_Buffer[8,1] = b[7] / 5040
        Fy_xy_Buffer[9,1] = b[8] / 40320
        Fy_xy_Buffer[1,3] = (-b[2] - b[1] * h) / 2
        Fy_xy_Buffer[2,3] = (-b[3] - b[2] * h + b[1] * h2) / 2
        Fy_xy_Buffer[3,3] = (-b[4] - b[3] * h + 2 * b[2] * h2 - 2 * b[1] * h3) / 4
        Fy_xy_Buffer[4,3] = (-b[5] + h * (-b[4] + 3 * h * (b[3] - 2 * b[2] * h + 2 * b[1] * h2))) / 12
        Fy_xy_Buffer[5,3] = (-b[6] - h * (b[5] + 4 * h * (-b[4] + 3 * h * (b[3] - 2 * b[2] * h + 2 * b[1] * h2)))) / 48
        Fy_xy_Buffer[6,3] = (-b[7] + h * (-b[6] + 5 * h * (b[5] + 4 * h * (-b[4] + 3 * h * (b[3] - 2 * b[2] * h + 2 * b[1] * h2))))) / 240
        Fy_xy_Buffer[7,3] = (-b[8] - h * (b[7] + 6 * h * (-b[6] + 5 * h * (b[5] + 4 * h * (-b[4] + 3 * h * (b[3] - 2 * b[2] * h + 2 * b[1] * h2)))))) / 1440
        Fy_xy_Buffer[8,3] = (h * (-b[8] + 7 * h * (b[7] + 6 * h * (-b[6] + 5 * h * (b[5] + 4 * h * (-b[4] + 3 * h * (b[3] - 2 * b[2] * h + 2 * b[1] * h2))))))) / 10080
        Fy_xy_Buffer[9,3] = (h2 * (b[8] - 7 * h * (b[7] + 6 * h * (-b[6] + 5 * h * (b[5] + 4 * h * (-b[4] + 3 * h * (b[3] - 2 * b[2] * h + 2 * b[1] * h2))))))) / 10080
    
        Fy_xy_Buffer[1,5] = (b[4] + h*(2*b[3] + h*(-b[2] + b[1]*h)))/24
        Fy_xy_Buffer[2,5] = (b[5] + h*(2*b[4] - 3*h*(b[3] - b[2]*h + b[1]*h2)))/24
        Fy_xy_Buffer[3,5] = (b[6] + h*(2*b[5] + h*(-5*b[4] + 9*b[3]*h - 12*b[2]*h2 + 12*b[1]*h3)))/48
        Fy_xy_Buffer[4,5] = (b[7] + h*(2*b[6] + h*(-7*b[5] + 19*b[4]*h - 39*b[3]*h2 + 60*b[2]*h3 - 60*b[1]*h4)))/144
        Fy_xy_Buffer[5,5] = (b[8] + h*(2*b[7] + 3*h*(-3*b[6] + 11*b[5]*h - 32*b[4]*h2 + 72*b[3]*h3 - 120*b[2]*h4 + 120*b[1]*h5)))/576
        Fy_xy_Buffer[6,5] = (h*(2*b[8] + h*(-11*b[7] - 3*h*(-17*b[6] + 5*h*(13*b[5] + 8*h*(-5*b[4] + 3*h*(4*b[3] - 7*b[2]*h + 7*b[1]*h2)))))))/2880
        Fy_xy_Buffer[7,5] = (h2*(-13*b[8] + h*(73*b[7] + 12*h*(-29*b[6] + 5*h*(23*b[5] + 2*h*(-37*b[4] + 3*h*(31*b[3] + 56*h*(-b[2] + b[1]*h))))))))/17280
        Fy_xy_Buffer[1,7] = (-b[6] - 3*h*(b[5] + h*(-b[4] + 2*b[3]*h - 3*b[2]*h2 + 3*b[1]*h3)))/720
        Fy_xy_Buffer[2,7] = (-b[7] - 3*h*(b[6] + h*(-2*b[5] + h*(4*b[4] - 9*b[3]*h + 15*b[2]*h2 - 15*b[1]*h3))))/720
        Fy_xy_Buffer[3,7] = (-b[8] - 3*h*(b[7] + h*(-3*b[6] + 8*b[5]*h - 21*b[4]*h2 + 51*b[3]*h3 - 90*b[2]*h4 + 90*b[1]*h5)))/1440
        Fy_xy_Buffer[4,7] = (h*(-b[8] + h*(4*b[7] + h*(-14*b[6] + 45*b[5]*h - 135*b[4]*h2 + 345*b[3]*h3 - 630*b[2]*h4 + 630*b[1]*h5))))/ 1440
        Fy_xy_Buffer[5,7] = (h2*(5*b[8] + h*(-22*b[7] + 3*h*(29*b[6] - 5*h*(21*b[5] + 4*h*(-17*b[4] + 45*b[3]*h - 84*b[2]*h2 + 84*b[1]*h3))))))/5760
        Fy_xy_Buffer[1,9] = (b[8] + h*(4*b[7] + 3*h*(-2*b[6] + 6*b[5]*h - 17*b[4]*h2 + 42*b[3]*h3 - 75*b[2]*h4 + 75*b[1]*h5)))/40320
        Fy_xy_Buffer[2,9] = (h*(4*b[8] - 5*h*(2*b[7] + 3*h*(-2*b[6] + 7*b[5]*h - 22*b[4]*h2 + 57*b[3]*h3 - 105*b[2]*h4 + 105*b[1]*h5))))/ 40320.;
        Fy_xy_Buffer[3,9] = (h2*(-14*b[8] + 5*h*(10*b[7] + 3*h*(-13*b[6] + 50*b[5]*h - 167*b[4]*h2 + 447*b[3]*h3 - 840*b[2]*h4 + 840*b[1]*h5))))/80640
        Fy_xy_Buffer[1,11] = (h2*(2*b[8] + h*(-8*b[7] - 3*h*(-11*b[6] + h*(43*b[5] + 5*h*(-29*b[4] + 3*h*(26*b[3] + 49*h*(-b[2] + b[1]*h))))))))/725760
    end
    Fx_xy = copy(Fx_xy_Buffer)
    Fy_xy = copy(Fy_xy_Buffer)
    return Fx_xy, Fy_xy, expansionOrder
end

function rotate_coordinates(coord, angle)
    # Constants for comparing floating point numbers
    epsilon = 1e-12

    if angle == 0 || abs(abs(angle) - 2*pi) < epsilon
        return
    end

    if abs(abs(angle) - pi) < epsilon
        cos_a = -1
        sin_a = 0
    elseif abs(angle - pi/2) < epsilon
        cos_a = 0
        sin_a = 1
    elseif abs(angle + pi/2) < epsilon
        cos_a = 0
        sin_a = -1
    else
        cos_a = cos(angle)
        sin_a = sin(angle)
    end

    x, xp, y, yp = coord[1], coord[2], coord[3], coord[4]
    new_coord = [x * cos_a + y * sin_a, xp * cos_a + yp * sin_a, -x * sin_a + y * cos_a, -xp * sin_a + yp * cos_a, coord[5], coord[6]]
    return new_coord
end

function computeEtiltCentroidOffset(rho0, angle, etilt, tilt)
    # compute final offsets due to error-tilt of the magnet
    # see pages 90-93 of notebook 1 about this
    
    if etilt == 0
        dcoord_etilt = zeros(6)
        return dcoord_etilt
    end

    q1a = (1 - cos(angle)) * rho0 * (cos(etilt) - 1)
    q2a = 0
    q3a = (1 - cos(angle)) * rho0 * sin(etilt)
    qp1 = sin(angle) * cos(etilt)
    qp2 = cos(angle)
    k = sqrt(qp1^2 + qp2^2)
    qp1 /= k
    qp2 /= k
    qp3 = sin(angle) * sin(etilt) / k
    tan_alpha = 1 / tan(angle) / cos(etilt)
    q1b = q1a * tan_alpha / (tan(angle) + tan_alpha)
    q2b = -q1b * tan(angle)
    dz = sqrt((q1b - q1a)^2 + (q2b - q2a)^2)
    q3b = q3a + qp3 * dz

    dcoord_etilt = [sqrt(q1b^2 + q2b^2), tan(atan(tan_alpha) - (pi/2 - angle)), q3b, qp3, dz * sqrt(1 + qp3^2), 0]


    # Debugging information can be handled with Julia's logging or print statements
    # println("pre-tilt offsets due to ETILT=$etilt:  $(dcoord_etilt[1]) $(dcoord_etilt[2]) $(dcoord_etilt[3]) $(dcoord_etilt[4]) $(dcoord_etilt[5])")

    # rotate by tilt to get into same frame as bend equations.
    new_dcoord = rotate_coordinates(dcoord_etilt, tilt)
    return new_dcoord

    # println("offsets due to ETILT=$etilt:  $(dcoord_etilt[1]) $(dcoord_etilt[2]) $(dcoord_etilt[3]) $(dcoord_etilt[4]) $(dcoord_etilt[5])")
end

function determine_bend_flags(elem, edge1_effects::Int, edge2_effects::Int)
    bend_flags = 0

    # Traversing predecessors
    other = elem.pred
    while other !== nothing
        if other.type == elem.type && other.name == elem.name
            bend_flags |= SAME_BEND_PRECEDES
            break
        end
        if other.type == T_MARK || other.type == T_WATCH
            other = other.pred
        else
            break
        end
    end

    # Traversing successors
    other = elem.succ
    while other !== nothing
        if other.type == elem.type && other.name == elem.name
            bend_flags |= SAME_BEND_FOLLOWS
            break
        end
        if other.type == T_MARK || other.type == T_WATCH
            other = other.succ
        else
            break
        end
    end

    if edge1_effects != 0 && (bend_flags & SAME_BEND_PRECEDES) == 0
        bend_flags |= BEND_EDGE1_EFFECTS
    end
    if edge2_effects != 0 && (bend_flags & SAME_BEND_FOLLOWS) == 0
        bend_flags |= BEND_EDGE2_EFFECTS
    end

    bend_flags |= BEND_EDGE_DETERMINED
    return bend_flags
end

function apply_edge_effects(x, xp, y, yp, rho, n, beta, he, psi, which_edge)
    h = 1 / rho
    tan_beta = tan(beta)
    R21 = h * tan_beta
    R43 = -h * tan(beta - psi)

    h2 = h^2
    tan2_beta = tan_beta^2
    sec_beta = 1 / cos(beta)
    sec2_beta = sec_beta^2

    T111 = which_edge * h / 2 * tan2_beta
    T133 = -which_edge * h / 2 * sec2_beta
    T211 = which_edge == -1 ? -n * h2 * tan_beta : -h2 * (n + tan2_beta / 2) * tan_beta
    T331 = -which_edge * h * tan2_beta
    T221 = T331
    T441 = -T331
    T233 = which_edge == -1 ? h2 * (n + 0.5 + tan2_beta) * tan_beta : h2 * (n - tan2_beta / 2) * tan_beta
    T243 = which_edge * h * tan2_beta
    T431 = h2 * (2 * n + (which_edge == 1 ? sec2_beta : 0)) * tan_beta
    T432 = which_edge * h * sec2_beta

    if he != 0
        term = h / 2 * he * sec2_beta * sec_beta
        T211 += term
        T233 -= term
        T431 -= 2 * term
    end

    x0, xp0, y0, yp0 = copy(x), copy(xp), copy(y), copy(yp)
    x = x0 + T111 * x0^2 + T133 * y0^2
    xp = xp0 + R21 * x0 + T211 * x0^2 + T221 * x0 * xp0 + T233 * y0^2 + T243 * y0 * yp0
    y = y0 + T331 * x0 * y0
    yp = yp0 + R43 * y0 + T441 * yp0 * x0 + T431 * x0 * y0 + T432 * xp0 * y0
    return x, xp, y, yp
end

function computeCSBENDFields(x, y, Fx_xy, Fy_xy, expansionOrder1)
    xp = [x^(i - 1) for i in 1:expansionOrder1]
    yp = [y^(i - 1) for i in 1:expansionOrder1]

    sumFX = 0
    sumFY = 0
    # for i in 1:expansionOrder1
    #     for j in 2:2:(expansionOrder1 - i)
    #         sumFX += Fx_xy[i, j] * xp[i] * yp[j]
    #         sumFY += Fy_xy[i, j] * xp[i] * yp[j]
    #     end
    # end
    for i in 1:expansionOrder1
        for j in 2:2:(expansionOrder1 - i + 1)
            if Fx_xy[i, j] != 0
                sumFX += Fx_xy[i, j] * xp[i] * yp[j]
            end
        end
    end
    for i in 1:expansionOrder1
        for j in 1:2:(expansionOrder1 - i + 1)
            if Fy_xy[i, j] != 0
                sumFY += Fy_xy[i, j] * xp[i] * yp[j]
            end
        end
    end
    
    
    
    return sumFX, sumFY
end

function convertToDipoleCanonicalCoordinates(Qi, rho)
    f = (1 + Qi[6]) / sqrt((1 + Qi[1] / rho)^2 + Qi[2]^2 + Qi[4]^2)
    new_Qi = [Qi[1], Qi[2]*f, Qi[3], Qi[4]*f, Qi[5], Qi[6]]
    return new_Qi 
end
function convertFromDipoleCanonicalCoordinates(Qi, rho)
    f = (1 + Qi[1]/rho)/sqrt((1+Qi[6])^2-Qi[2]^2-Qi[4]^2);
    new_Qi = [Qi[1], Qi[2]*f, Qi[3], Qi[4]*f, Qi[5], Qi[6]]
    return new_Qi
end

function dipoleFringe(vec, h, inFringe, higherOrder)
    # vec = [x, qx, y, qy, s, delta]
    x = vec[1]
    px = vec[2]
    y = vec[3]
    py = vec[4]
    delta = vec[6]
    dx = dpx = dy = dpy = ds = 0

    a = inFringe * h / (8 * (1 + delta))

    if higherOrder != 0
        xr = [x^i for i in 0:10]
        yr = [y^i for i in 0:10]
        ar = [a^i for i in 0:10]

        dx = (a * (7 * ar[8] * xr[10] + 7 * ar[9] * xr[11] + ar[6] * xr[8] * (7 + 27 * ar[3] * yr[3]) +
        ar[7] * xr[9] * (7 + 30 * ar[3] * yr[3]) + ar[5] * xr[7] * (7 + 24 * ar[3] * yr[3] +
        30 * ar[5] * yr[5]) + ar[4] * xr[6] * (7 + 21 * ar[3] * yr[3] + 99 * ar[5] * yr[5]) +
        3 * a * x * yr[3] * (-7 + 35 * ar[3] * yr[3] - 91 * ar[5] * yr[5] + 246 * ar[7] * yr[7]) +
        ar[3] * xr[5] * (7 + 21 * ar[3] * yr[3] - 90 * ar[5] * yr[5] + 1140 * ar[7] * yr[7]) +
        xr[4] * (7 * a + 189 * ar[6] * yr[5] - 975 * ar[8] * yr[7]) +
        3 * yr[3] * (7 - 7 * ar[3] * yr[3] + 21 * ar[5] * yr[5] - 39 * ar[7] * yr[7] +
        82 * ar[9] * yr[9]) - 7 * xr[3] * (-1 - 6 * ar[3] * yr[3] + 21 * ar[5] * yr[5] -
        96 * ar[7] * yr[7] + 315 * ar[9] * yr[9])))/7

        dpx = (a * (2 * py * y * (7 - 7 * a * x - 28 * ar[3] * yr[3] + 70 * x * ar[4] * yr[3] + 56 * x * ar[6] * (xr[3] -
        6 * yr[3]) * yr[3] + 84 * ar[5] * yr[3] * (-xr[3] + yr[3]) +
        3 * x * ar[8] * yr[3] * (xr[5] - 262 * xr[3] * yr[3] + 409 * yr[5]) -
        4 * ar[7] * (5 * xr[5] * yr[3] - 162 * xr[3] * yr[5] + 57 * yr[7]) +
        20 * ar[9] * (33 * xr[5] * yr[5] - 162 * xr[3] * yr[7] + 29 * yr[9])) +
        px * (-3 * ar[8] * xr[7] * yr[3] - 24 * ar[7] * xr[6] * yr[3] * (-1 + 55 * ar[3] * yr[3]) +
        3 * ar[6] * xr[5] * yr[3] * (-28 + 655 * ar[3] * yr[3]) +
        24 * ar[5] * xr[4] * yr[3] * (7 - 90 * ar[3] * yr[3] + 630 * ar[5] * yr[5]) +
        3 * a * yr[3] * (-21 + 70 * ar[3] * yr[3] - 196 * ar[5] * yr[5] + 513 * ar[7] * yr[7]) -
        7 * a * xr[3] * (-1 + 30 * ar[3] * yr[3] - 240 * ar[5] * yr[5] + 1227 * ar[7] * yr[7]) -
        2 * x * (7 - 84 * ar[3] * yr[3] + 420 * ar[5] * yr[5] - 1596 * ar[7] * yr[7] +
        5220 * ar[9] * yr[9]))))/7

        dy = -(a * y * (ar[8] * xr[7] * yr[3] + 8 * ar[7] * xr[6] * yr[3] * (-1 + 33 * ar[3] * yr[3]) -
        8 * ar[5] * xr[4] * yr[3] * (7 - 54 * ar[3] * yr[3] + 270 * ar[5] * yr[5]) +
        xr[5] * (28 * ar[6] * yr[3] - 393 * ar[8] * yr[5]) + 3 * a * yr[3] * (7 - 14 * ar[3] * yr[3] +
        28 * ar[5] * yr[5] - 57 * ar[7] * yr[7]) + a * xr[3] * (-7 + 70 * ar[3] * yr[3] -
        336 * ar[5] * yr[5] + 1227 * ar[7] * yr[7]) + 2 * x * (7 - 28 * ar[3] * yr[3] +
        84 * ar[5] * yr[5] - 228 * ar[7] * yr[7] + 580 * ar[9] * yr[9])))/7

        dpy = (a * (-6 * px * y * (7 - 7 * a * x + 14 * ar[3] * (xr[3] - yr[3]) + 70 * x * ar[4] * yr[3] +
        7 * ar[5] * (xr[5] - 14 * xr[3] * yr[3] + 9 * yr[5]) + 7 * ar[6] * (xr[6] +
        18 * xr[4] * yr[3] - 39 * x * yr[5]) + 4 * ar[7] * (2 * xr[7] - 15 * xr[5] * yr[3] +
        168 * xr[3] * yr[5] - 39 * yr[7]) + ar[8] * (9 * xr[8] + 66 * xr[6] * yr[3] -
        975 * xr[4] * yr[5] + 984 * x * yr[7]) + 10 * ar[9] * (xr[9] + 2 * xr[7] * yr[3] +
        114 * xr[5] * yr[5] - 294 * xr[3] * yr[7] + 41 * yr[9])) + py * (63 * ar[8] * xr[9] +
        70 * ar[9] * xr[10] + 7 * ar[6] * xr[7] * (7 + 27 * ar[3] * yr[3]) +
        8 * ar[7] * xr[8] * (7 + 30 * ar[3] * yr[3]) + 6 * ar[5] * xr[6] * (7 + 24 * ar[3] * yr[3] +
        30 * ar[5] * yr[5]) + 5 * ar[4] * xr[5] * (7 + 21 * ar[3] * yr[3] + 99 * ar[5] * yr[5]) +
        3 * a * xr[3] * (7 + 189 * ar[5] * yr[5] - 975 * ar[7] * yr[7]) +
        3 * a * yr[3] * (-7 + 35 * ar[3] * yr[3] - 91 * ar[5] * yr[5] + 246 * ar[7] * yr[7]) +
        4 * ar[3] * xr[4] * (7 + 21 * ar[3] * yr[3] - 90 * ar[5] * yr[5] + 1140 * ar[7] * yr[7]) -
        14 * x * (-1 - 6 * ar[3] * yr[3] + 21 * ar[5] * yr[5] - 96 * ar[7] * yr[7] +
        315 * ar[9] * yr[9]))))/7

    else
        dx = a * (x^2 + 3 * y^2)
        dpx = 2 * a * (-px * x + py * y)
        dy = -2 * a * x * y
        dpy = 2 * a * (py * x - 3 * px * y)
    end

    ds = (a / (1 + delta)) * (2 * py * x * y - px * x^2 - 3 * px * y^2)

    new_vec = [vec[1]+dx, vec[2]+dpx, vec[3]+dy, vec[4]+dpy, vec[5]+ds, vec[6]]
    # vec[1] += dx
    # vec[2] += dpx
    # vec[3] += dy
    # vec[4] += dpy
    # vec[5] += ds
    return new_vec
end

function addRadiationKick(Qx, Qy, dPoP, sigmaDelta2,
    x, h0, Fx, Fy, ds, radCoef, dsISR, isrCoef,
    distributionBased, includeOpeningAngle, meanPhotonsPerMeter,
    normalizedCriticalEnergy0, Po)

    f = (1 + x * h0) / sqrt((1 + dPoP)^2 - Qx^2 - Qy^2)
    xp = Qx * f
    yp = Qy * f
    dsFactor = sqrt((1 + x * h0)^2 + xp^2 + yp^2)
    F2 = Fx^2 + Fy^2

    if distributionBased == 0 # now must be 0
        deltaFactor = (1 + dPoP)^2
        Qx /= (1 + dPoP)
        Qy /= (1 + dPoP)
        if radCoef != 0
            dPoP -= radCoef * deltaFactor * F2 * ds * dsFactor
        end
        # if isrCoef > 0
        # # gauss_rn_lim and random_2 are not defined yet
        #     dPoP -= isrCoef * deltaFactor * F2^0.75 * sqrt(dsISR * dsFactor) * gauss_rn_lim(0.0, 1.0, 3.0, random_2())
        # end
        if !isnothing(sigmaDelta2) 
            sigmaDelta2 += (isrCoef * deltaFactor)^2 * F2^1.5 * dsISR * dsFactor
        end
        Qx *= (1 + dPoP)
        Qy *= (1 + dPoP)
    else
        F = sqrt(F2)
        # incomplete code
        # nMean = meanPhotonsPerMeter * dsISR * dsFactor * F
        # nEmitted = inversePoissonCDF(nMean, random_2(1))
        # normalizedCriticalEnergy = normalizedCriticalEnergy0 * (1 + dPoP) * F

        # for i in 1:nEmitted
        #     y = pickNormalizedPhotonEnergy(random_2(1))
        #     dDelta = normalizedCriticalEnergy * (1 + dPoP) * y
        #     photonCount += 1
        #     energyCount += y
        #     dPoP -= dDelta
        #     if includeOpeningAngle
        #         logy = log10(y)
        #         thetaRms = dDelta * 10^(-0.2418673276661232 + logy * (-0.4472680955382907 + logy * (-0.045353504248236 - logy * 0.006181818621278201))) / Po
        #         xp += thetaRms * gauss_rn_lim(0.0, 1.0, 3.0, random_2())
        #         yp += thetaRms * gauss_rn_lim(0.0, 1.0, 3.0, random_2())
        #     end
        # end
        # f = (1 + dPoP) / EXSQRT((1 + x * h0)^2 + xp^2 + yp^2, sqrtOrder)
        # Qx = xp * f
        # Qy = yp * f
    end

    return Qx, Qy, dPoP, sigmaDelta2
end

function integrate_csbend_ord2(Qi, sigmaDelta2, s, n, rho0, p0, Fx_xy, Fy_xy, rho_actual, 
                                rad_coef, isrConstant, expansionOrder)
    particle_lost = false
    s_lost = 0.0

    if n < 1
        error("invalid number of steps (integrate_csbend_ord4)")
    end

    Qf = [Qi[i] for i in 1:6]
    Qf_Buffer = Zygote.Buffer(Qf)
    for i in 1:6
        Qf_Buffer[i] = Qf[i]
    end

    dist = 0.0
    ds = s/n
    dsh = ds/2
    for i in 1:n
        if i == 1
            # First drift
            f = (1 + Qf_Buffer[6])^2 - Qf_Buffer[4]^2
            if f <= 0
                particle_lost = true
                s_lost = dist
                Qf = copy(Qf_Buffer)
                return Qf, particle_lost, s_lost
            end
            f = sqrt(f)
            if abs(Qf_Buffer[2]/f) > 1
                particle_lost = true
                s_lost = dist
                Qf = copy(Qf_Buffer)
                return Qf, particle_lost, s_lost
            end
            sin_phi = Qf_Buffer[2] / f
            phi = asin(sin_phi)
            sine = sin(dsh / rho0 + phi)
            cosi = cos(dsh / rho0 + phi)
            if cosi == 0
                particle_lost = true
                s_lost = dist
                Qf = copy(Qf_Buffer)
                return Qf, particle_lost, s_lost
            end
            tang = sine / cosi
            cos_phi = cos(phi)
            Qf_Buffer[2] = f * sine
            factor = (rho0+Qf_Buffer[1])*cos_phi/f*(tang-sin_phi/cos_phi)
            Qf_Buffer[3] += Qf_Buffer[4] * factor
            dist += factor * (1 + Qf_Buffer[6])
            f = cos_phi / cosi
            Qf_Buffer[1] = rho0 * (f - 1) + f * Qf_Buffer[1]
        end

        # First kick
        x = Qf_Buffer[1]
        y = Qf_Buffer[3]
        Fx, Fy = computeCSBENDFields(x, y, Fx_xy, Fy_xy, expansionOrder+1)

        Qf_Buffer[2] += -ds * (1 + Qf_Buffer[1]/rho0) * Fy / rho_actual
        Qf_Buffer[4] += ds * (1 + Qf_Buffer[1]/rho0) * Fx / rho_actual
        if rad_coef != 0 || isrConstant != 0
            Qx, Qy, dPoP, sigmaDelta2 = addRadiationKick(Qf_Buffer[2], Qf_Buffer[4], Qf_Buffer[6], sigmaDelta2,
                                        Qf_Buffer[1], 1/rho0, Fx, Fy, ds, rad_coef, s/3, 0, 0, 0, 0, 0, p0)
            Qf_Buffer[2] = Qx
            Qf_Buffer[4] = Qy
            Qf_Buffer[6] = dPoP
        end

        # second drift
        if i == n
            f = (1 + Qf_Buffer[6])^2 - Qf_Buffer[4]^2
            if f <= 0
                particle_lost = true
                s_lost = dist
                Qf = copy(Qf_Buffer)
                return Qf, particle_lost, s_lost
            end
            f = sqrt(f)
            if abs(Qf_Buffer[2]/f) > 1
                particle_lost = true
                s_lost = dist
                Qf = copy(Qf_Buffer)
                return Qf, particle_lost, s_lost
            end
            sin_phi = Qf_Buffer[2] / f
            phi = asin(sin_phi)
            sine = sin(dsh / rho0 + phi)
            cosi = cos(dsh / rho0 + phi)
            if cosi == 0
                particle_lost = true
                s_lost = dist
                Qf = copy(Qf_Buffer)
                return Qf, particle_lost, s_lost
            end
            tang = sine / cosi
            cos_phi = cos(phi)
            Qf_Buffer[2] = f * sine
            factor = (rho0+Qf_Buffer[1])*cos_phi/f*(tang-sin_phi/cos_phi)
            Qf_Buffer[3] += Qf_Buffer[4] * factor
            dist += factor * (1 + Qf_Buffer[6])
            f = cos_phi / cosi
            Qf_Buffer[1] = rho0 * (f - 1) + f * Qf_Buffer[1]
        else
            f = (1 + Qf_Buffer[6])^2 - Qf_Buffer[4]^2
            if f <= 0
                particle_lost = true
                s_lost = dist
                Qf = copy(Qf_Buffer)
                return Qf, particle_lost, s_lost
            end
            f = sqrt(f)
            if abs(Qf_Buffer[2]/f) > 1
                particle_lost = true
                s_lost = dist
                Qf = copy(Qf_Buffer)
                return Qf, particle_lost, s_lost
            end
            sin_phi = Qf_Buffer[2] / f
            phi = asin(sin_phi)
            sine = sin(ds / rho0 + phi)
            cosi = cos(ds / rho0 + phi)
            if cosi == 0
                particle_lost = true
                s_lost = dist
                Qf = copy(Qf_Buffer)
                return Qf, particle_lost, s_lost
            end
            tang = sine / cosi
            cos_phi = cos(phi)
            Qf_Buffer[2] = f * sine
            factor = (rho0+Qf_Buffer[1])*cos_phi/f*(tang-sin_phi/cos_phi)
            Qf_Buffer[3] += Qf_Buffer[4] * factor
            dist += factor * (1 + Qf_Buffer[6])
            f = cos_phi / cosi
            Qf_Buffer[1] = rho0 * (f - 1) + f * Qf_Buffer[1]
        end
    end
    Qf_Buffer[5] += dist
    Qf = copy(Qf_Buffer)
    return Qf, particle_lost, s_lost, sigmaDelta2
end

function integrate_csbend_ord4(Qi, sigmaDelta2, s, n, rho0, p0, Fx_xy, Fy_xy, rho_actual, 
                                rad_coef, isrConstant, expansionOrder)
    BETA = 1.25992104989487316477
    particle_lost = false
    s_lost = 0.0

    if n < 1
        error("invalid number of steps (integrate_csbend_ord4)")
    end

    Qf = [Qi[i] for i in 1:6]
    Qf_Buffer = Zygote.Buffer(Qf)
    for i in 1:6
        Qf_Buffer[i] = Qf[i]
    end

    dist = 0.0
    s /= n
    for i in 1:n
        # First drift
        dsh = s / 2 / (2 - BETA)
        f = (1 + Qf_Buffer[6])^2 - Qf_Buffer[4]^2
        if f <= 0
            particle_lost = true
            s_lost = dist
            Qf = copy(Qf_Buffer)
            return Qf, particle_lost, s_lost
        end
        f = sqrt(f)
        if abs(Qf_Buffer[2]/f) > 1
            particle_lost = true
            s_lost = dist
            Qf = copy(Qf_Buffer)
            return Qf, particle_lost, s_lost
        end
        sin_phi = Qf_Buffer[2] / f
        phi = asin(sin_phi)
        sine = sin(dsh / rho0 + phi)
        cosi = cos(dsh / rho0 + phi)
        if cosi == 0
            particle_lost = true
            s_lost = dist
            Qf = copy(Qf_Buffer)
            return Qf, particle_lost, s_lost
        end
        tang = sine / cosi
        cos_phi = cos(phi)
        Qf_Buffer[2] = f * sine
        factor = (rho0+Qf_Buffer[1])*cos_phi/f*(tang-sin_phi/cos_phi)
        Qf_Buffer[3] += Qf_Buffer[4] * factor
        dist += factor * (1 + Qf_Buffer[6])
        f = cos_phi / cosi
        Qf_Buffer[1] = rho0 * (f - 1) + f * Qf_Buffer[1]

        # First kick
        ds = s / (2 - BETA)
        # Assuming computeCSBENDFields and addRadiationKick are defined in Julia
        x = Qf_Buffer[1]
        y = Qf_Buffer[3]
        Fx, Fy = computeCSBENDFields(x, y, Fx_xy, Fy_xy, expansionOrder+1)

        Qf_Buffer[2] += -ds * (1 + Qf_Buffer[1]/rho0) * Fy / rho_actual
        Qf_Buffer[4] += ds * (1 + Qf_Buffer[1]/rho0) * Fx / rho_actual
        if rad_coef != 0 || isrConstant != 0
            Qx, Qy, dPoP, sigmaDelta2 = addRadiationKick(Qf_Buffer[2], Qf_Buffer[4], Qf_Buffer[6], sigmaDelta2,
                                        Qf_Buffer[1], 1/rho0, Fx, Fy, ds, rad_coef, s/3, 0, 0, 0, 0, 0, p0)
            Qf_Buffer[2] = Qx
            Qf_Buffer[4] = Qy
            Qf_Buffer[6] = dPoP
        end

        # second drift
        dsh = s*(1 - BETA)/(2 - BETA)/2
        f = (1 + Qf_Buffer[6])^2 - Qf_Buffer[4]^2
        if f <= 0
            particle_lost = true
            s_lost = dist
            Qf = copy(Qf_Buffer)
            return Qf, particle_lost, s_lost
        end
        f = sqrt(f)
        if abs(Qf_Buffer[2]/f) > 1
            particle_lost = true
            s_lost = dist
            Qf = copy(Qf_Buffer)
            return Qf, particle_lost, s_lost
        end
        sin_phi = Qf_Buffer[2] / f
        phi = asin(sin_phi)
        sine = sin(dsh / rho0 + phi)
        cosi = cos(dsh / rho0 + phi)
        if cosi == 0
            particle_lost = true
            s_lost = dist
            Qf = copy(Qf_Buffer)
            return Qf, particle_lost, s_lost
        end
        tang = sine / cosi
        cos_phi = cos(phi)
        Qf_Buffer[2] = f * sine
        factor = (rho0+Qf_Buffer[1])*cos_phi/f*(tang-sin_phi/cos_phi)
        Qf_Buffer[3] += Qf_Buffer[4] * factor
        dist += factor * (1 + Qf_Buffer[6])
        f = cos_phi / cosi
        Qf_Buffer[1] = rho0 * (f - 1) + f * Qf_Buffer[1]

        # second kick
        ds = -s*BETA/(2 - BETA)
        x = Qf_Buffer[1]
        y = Qf_Buffer[3]
        Fx, Fy = computeCSBENDFields(x, y, Fx_xy, Fy_xy, expansionOrder+1)

        Qf_Buffer[2] += -ds * (1 + Qf_Buffer[1]/rho0) * Fy / rho_actual
        Qf_Buffer[4] += ds * (1 + Qf_Buffer[1]/rho0) * Fx / rho_actual
        if rad_coef != 0 || isrConstant != 0
            Qx, Qy, dPoP, sigmaDelta2 = addRadiationKick(Qf_Buffer[2], Qf_Buffer[4], Qf_Buffer[6], sigmaDelta2,
                                        Qf_Buffer[1], 1/rho0, Fx, Fy, ds, rad_coef, s/3, 0, 0, 0, 0, 0, p0)
            Qf_Buffer[2] = Qx
            Qf_Buffer[4] = Qy
            Qf_Buffer[6] = dPoP
        end

        # third drift
        dsh = s*(1 - BETA)/(2 - BETA)/2
        f = (1 + Qf_Buffer[6])^2 - Qf_Buffer[4]^2
        if f <= 0
            particle_lost = true
            s_lost = dist
            Qf = copy(Qf_Buffer)
            return Qf, particle_lost, s_lost
        end
        f = sqrt(f)
        if abs(Qf_Buffer[2]/f) > 1
            particle_lost = true
            s_lost = dist
            Qf = copy(Qf_Buffer)
            return Qf, particle_lost, s_lost
        end
        sin_phi = Qf_Buffer[2] / f
        phi = asin(sin_phi)
        sine = sin(dsh / rho0 + phi)
        cosi = cos(dsh / rho0 + phi)
        if cosi == 0
            particle_lost = true
            s_lost = dist
            Qf = copy(Qf_Buffer)
            return Qf, particle_lost, s_lost
        end
        tang = sine / cosi
        cos_phi = cos(phi)
        Qf_Buffer[2] = f * sine
        factor = (rho0+Qf_Buffer[1])*cos_phi/f*(tang-sin_phi/cos_phi)
        Qf_Buffer[3] += Qf_Buffer[4] * factor 
        dist += factor * (1 + Qf_Buffer[6])
        f = cos_phi / cosi
        Qf_Buffer[1] = rho0 * (f - 1) + f * Qf_Buffer[1]

        # third kick
        ds = s / (2 - BETA)
        Qf_Buffer[2] += -ds * (1 + Qf_Buffer[1]/rho0) * Fy / rho_actual
        Qf_Buffer[4] += ds * (1 + Qf_Buffer[1]/rho0) * Fx / rho_actual
        if rad_coef != 0 || isrConstant != 0
            Qx, Qy, dPoP, sigmaDelta2 = addRadiationKick(Qf_Buffer[2], Qf_Buffer[4], Qf_Buffer[6], sigmaDelta2,
                                                        Qf_Buffer[1], 1/rho0, Fx, Fy, ds, rad_coef, s/3, 0, 0, 0, 0, 0, p0)
            Qf_Buffer[2] = Qx
            Qf_Buffer[4] = Qy
            Qf_Buffer[6] = dPoP
        end

        # fourth drift
        dsh = s / 2 / (2 - BETA)
        f = (1 + Qf_Buffer[6])^2 - Qf_Buffer[4]^2
        if f <= 0
            particle_lost = true
            s_lost = dist
            Qf = copy(Qf_Buffer)
            return Qf, particle_lost, s_lost
        end
        f = sqrt(f)
        if abs(Qf_Buffer[2]/f) > 1
            particle_lost = true
            s_lost = dist
            Qf = copy(Qf_Buffer)
            return Qf, particle_lost, s_lost
        end
        sin_phi = Qf_Buffer[2] / f
        phi = asin(sin_phi)
        sine = sin(dsh / rho0 + phi)
        cosi = cos(dsh / rho0 + phi)
        if cosi == 0
            particle_lost = true
            s_lost = dist
            Qf = copy(Qf_Buffer)
            return Qf, particle_lost, s_lost
        end
        tang = sine / cosi
        cos_phi = cos(phi)
        Qf_Buffer[2] = f * sine
        factor = (rho0+Qf_Buffer[1])*cos_phi/f*(tang-sin_phi/cos_phi)
        Qf_Buffer[3] += Qf_Buffer[4] * factor
        dist += factor * (1 + Qf_Buffer[6])
        f = cos_phi / cosi
        Qf_Buffer[1] = rho0 * (f - 1) + f * Qf_Buffer[1]
    end
    Qf_Buffer[5] += dist
    Qf = copy(Qf_Buffer)
    return Qf, particle_lost, s_lost, sigmaDelta2
end

function track_one_part(coord, n, he1, he2, e1, e2, dxf, dyf, dzf, sin_ttilt, cos_ttilt, dcoord_etilt, psi1, psi2, csbend, Po, rho0, rho_actual, rad_coef, isrConstant, 
                        sigmaDelta2, Fx_xy, Fy_xy, e1_kick_limit, e2_kick_limit, expansionOrder)
    x = coord[1]*cos_ttilt + coord[3]*sin_ttilt
    y = coord[3]*cos_ttilt - coord[1]*sin_ttilt
    xp = coord[2]*cos_ttilt + coord[4]*sin_ttilt
    yp = coord[4]*cos_ttilt - coord[2]*sin_ttilt
    s = coord[5]
    dp = coord[6]
    dp0 = dp

    if csbend.edge1_effects != 0
        rho = (1+dp)*rho_actual
        if csbend.edge_order < 2 || csbend.edge1_effects > 1
            delta_xp = tan(e1)/rho*x
                if e1_kick_limit > 0 && abs(delta_xp) > e1_kick_limit
                    delta_xp = sign(delta_xp)*e1_kick_limit
                end
            xp += delta_xp
            yp -= tan(e1 - psi1/(1+dp))/rho*y
        else
            x, xp, y, yp = apply_edge_effects(x, xp, y, yp, rho, n, e1, he1, psi1*(1+dp), -1)
        end
    end

    xp *= (1+x/rho0)
    yp *= (1+x/rho0)

    Qi = [x, xp, y, yp, 0, dp]

    if csbend.edge1_effects !=0 && e1!=0 && rad_coef!=0
        Fx, Fy = computeCSBENDFields(x, y, Fx_xy, Fy_xy, expansionOrder+1)
        dp_prime = -rad_coef * (Fx^2 + Fy^2) * (1 + dp)^2 * sqrt((1+x/rho0)^2 + xp^2 + yp^2)
        Qi = [Qi[1], Qi[2], Qi[3], Qi[4], Qi[5], Qi[6] - dp_prime*x*tan(e1)]
        # Qi[6] -= dp_prime*x*tan(e1)
    end

    Qi = convertToDipoleCanonicalCoordinates(Qi, rho0)

    if csbend.edge1_effects > 1
        Qi = dipoleFringe(Qi, rho0, -1, csbend.edge1_effects-2)
    end

    particle_lost = 0
    if particle_lost == 0
        if csbend.integration_order == 4
            # Qf, particle_lost, s_lost, sigmaDelta2 = integrate_csbend_ord4(Qi, sigmaDelta2, csbend.length, csbend.nSlice, 
                                                # rho0, Po, Fx_xy, Fy_xy, rho_actual, rad_coef, isrConstant, expansionOrder)
            Qf, particle_lost, s_lost, sigmaDelta2 = integrate_csbend_ord4(Qi, sigmaDelta2, csbend.length, csbend.nSlice, 
                                                rho0, Po, Fx_xy, Fy_xy, rho_actual, rad_coef, isrConstant, expansionOrder)
        elseif csbend.integration_order == 2
            Qf, particle_lost, s_lost, sigmaDelta2 = integrate_csbend_ord2(Qi, sigmaDelta2, csbend.length, csbend.nSlice, 
                                                rho0, Po, Fx_xy, Fy_xy, rho_actual, rad_coef, isrConstant, expansionOrder)
        else
            error("invalid integration order (track_through_csbend)")
        end
    end

    if csbend.edge2_effects > 1
        Qf = dipoleFringe(Qf, rho0, 1, csbend.edge2_effects-2)
    end
    Qf = convertFromDipoleCanonicalCoordinates(Qf, rho0)

    ######
    # if particle_lost
    #####

    if csbend.edge2_effects !=0 && e2!=0 && rad_coef!=0
        x = Qf[1]
        xp = Qf[2]
        y = Qf[3]
        yp = Qf[4]
        dp = Qf[6]
        Fx, Fy = computeCSBENDFields(x, y, Fx_xy, Fy_xy, expansionOrder+1)
        dp_prime = -rad_coef * (Fx^2 + Fy^2) * (1 + dp)^2 * sqrt((1+x/rho0)^2 + xp^2 + yp^2)
        Qf = [Qf[1], Qf[2], Qf[3], Qf[4], Qf[5], Qf[6] - dp_prime*Qf[1]*tan(e2)]
        # Qf[6] -= dp_prime*Qf[1]*tan(e2)
    end

    # final coordinates
    if rad_coef !=0 || isrConstant != 0
        p0 = Po*(1+dp0)
        beta0 = p0/sqrt(p0^2+1)
        p1 = Po*(1+Qf[6])
        beta1 = p1/sqrt(p1^2+1)
        s = beta1*s/beta0 + Qf[5]
    else
        s += Qf[5]
    end
    x = Qf[1]
    xp = Qf[2]
    y = Qf[3]
    yp = Qf[4]
    dp = Qf[6]

    xp /= (1+x/rho0)
    yp /= (1+x/rho0)

    if csbend.edge2_effects != 0
        rho = (1+dp)*rho_actual
        if csbend.edge_order < 2 || csbend.edge2_effects > 1
            delta_xp = tan(e2)/rho*x
            if e2_kick_limit > 0 && abs(delta_xp) > e2_kick_limit
                delta_xp = sign(delta_xp)*e2_kick_limit
            end
            xp += delta_xp
            yp -= tan(e2 - psi2/(1+dp))/rho*y
        else
        x, xp, y, yp = apply_edge_effects(x, xp, y, yp, rho, n, e2, he2, psi2*(1+dp), 1)
        end
    end

    new_coord = [x*cos_ttilt - y*sin_ttilt + dcoord_etilt[1] + dxf + dzf*coord[2],
                xp*cos_ttilt - yp*sin_ttilt + dcoord_etilt[2],
                y*cos_ttilt + x*sin_ttilt + dcoord_etilt[3] + dyf + dzf*coord[4],
                yp*cos_ttilt + xp*sin_ttilt + dcoord_etilt[4],
                s + dzf*sqrt(1+coord[2]^2+coord[4]^2),
                dp]
    # coord[1] = x*cos_ttilt - y*sin_ttilt + dcoord_etilt[1]
    # coord[2] = xp*cos_ttilt - yp*sin_ttilt + dcoord_etilt[2]
    # coord[3] = y*cos_ttilt + x*sin_ttilt + dcoord_etilt[3]
    # coord[4] = yp*cos_ttilt + xp*sin_ttilt + dcoord_etilt[4]
    # coord[5] = s
    # coord[6] = dp

    # coord[1] += dxf + dzf*coord[2]
    # coord[3] += dyf + dzf*coord[4]
    # coord[5] += dzf*sqrt(1+coord[2]^2+coord[4]^2)
    return new_coord, particle_lost, s_lost, sigmaDelta2
end

function track_through_csbend(part, n_part, csbend, p_error, Po, sigmaDelta2)
    # constants
    particleMass = 9.1093897e-31 # Electron mass (kg)
    particleCharge = 1.60217733e-19 # Electron charge (C)
    c_mks = 299792458.0 # Speed of light (m/s)
    epsilon_o = 8.854187817e-12 # Permittivity of vacuum (F/m)
    me_mev = 0.51099906 # Electron mass (MeV)
    particleRadius = particleCharge^2 / (4 * pi * epsilon_o * particleMass * c_mks^2) # Classical electron radius (m)

    largeRhoWarning = 0
    if csbend.angle == 0
        part1 = exactDrift(part, n_part, csbend.length)
        return part1
    end

    rho0 = csbend.length / csbend.angle
    if csbend.use_bn != 0
        b = [csbend.b1, csbend.b2, csbend.b3, csbend.b4, csbend.b5, csbend.b6, csbend.b7, csbend.b8]
    else
        b = [csbend.k1*rho0, csbend.k2*rho0, csbend.k3*rho0, csbend.k4*rho0, 
                csbend.k5*rho0, csbend.k6*rho0, csbend.k7*rho0, csbend.k8*rho0]
    end

    he1 = csbend.h1
    he2 = csbend.h2

    if csbend.angle < 0
        angle = -csbend.angle
        e1 = -csbend.e1
        e2 = -csbend.e2
        etilt = csbend.etilt
        tilt = csbend.tilt + pi
        rho0 = csbend.length / angle
        b = [i % 2 == 1 ? -b[i] : b[i] for i in 1:8]
    else
        angle = csbend.angle
        e1 = csbend.e1
        e2 = csbend.e2
        etilt = csbend.etilt
        tilt = csbend.tilt
        rho0 = csbend.length / angle
    end

    if rho0 > 1e6
        largeRhoWarning = 1
        println("CSBEND Warning: large bend radius, rho0 = $rho0. Treated as drift.")
        part1 = exactDrift(part, n_part, csbend.length)
        return part1
    end

    fse = csbend.fse
    h = 1 / rho0
    n = -b[1]/h
    if fse > -1
        rho_actual = 1/((1+fse)*h)
    else
        rho_actual = 1e16/h
    end

    e1_kick_limit = csbend.e1_kick_limit
    e2_kick_limit = csbend.e2_kick_limit
    if csbend.kick_limit_scaling != 0
        e1_kick_limit *= rho0/rho_actual
        e2_kick_limit *= rho0/rho_actual
    end

    # angles for fringe field
    hgap = csbend.hgap
    fint = csbend.fint
    Kg = 2*hgap * fint
    psi1 = Kg/rho_actual/cos(e1)*(1+sin(e1)^2)
    psi2 = Kg/rho_actual/cos(e2)*(1+sin(e2)^2)
    # rad_coef is d((P-Po)/Po)/ds for the on-axis, on-momentum particle, where po is the momentum of the central particle.
    if csbend.synch_rad !=0
        rad_coef = particleCharge^2*Po^3*(1+fse)^2/(6*PI*epsilon_o*c_mks^2*particleMass*rho0^2)
    else
        rad_coef = 0
    end

    isrConstant = particleRadius*sqrt(55.0/(24*sqrt(3))*Po^5* 137.0359895/abs(rho_actual)^3)
    if csbend.isr == 0 || (csbend.isr1Particle == 0 && n_part == 1)
        isrConstant *= -1
    end

    # distributionBasedRadiation = csbend.distributionBasedRadiation
    # if csbend.distributionBasedRadiation != 0
    #     meanPhotonsPerRadian0 = 5.0/(2.0*sqrt(3))*Po/137.0359895
    #     meanPhotonsPerMeter0 = (5*c_mks*Po*particleMass*particleRadius)/(2*sqrt(3)*hbar_mks*rho_actual)
    #     normalizedCriticalEnergy0 = 3.0/2*hbar_mks*c_mks*Po^3/abs(rho_actual)/(Po*particleMass*c_mks^2)
    #     includeOpeningAngle = csbend.includeOpeningAngle
    # end

    Fx_xy, Fy_xy, expansionOrder = computeCSBENDFieldCoefficients(b, h, csbend.nonlinear, csbend.expansionOrder)

    ttilt = tilt + etilt
    if ttilt == 0 
        cos_ttilt = 1
        sin_ttilt = 0
    elseif abs(abs(ttilt)-pi)<1e-12
        cos_ttilt = -1
        sin_ttilt = 0
    elseif abs(abs(ttilt-pi/2))<1e-12
        cos_ttilt = 0
        sin_ttilt = 1
    elseif abs(abs(ttilt+pi/2))<1e-12
        cos_ttilt = 0
        sin_ttilt = -1
    else
        cos_ttilt = cos(ttilt)
        sin_ttilt = sin(ttilt)
    end
    dcoord_etilt = computeEtiltCentroidOffset(rho0, angle, etilt, tilt)

    dxi = -csbend.dx
    dzi =  csbend.dz
    dyi = -csbend.dy

    dxf = csbend.dx * cos(csbend.angle) + csbend.dz * sin(csbend.angle)
    dzf = csbend.dx * sin(csbend.angle) - csbend.dz * cos(csbend.angle)
    dyf = csbend.dy

    i_top = n_part-1
    if !isnothing(sigmaDelta2) 
        sigmaDelta2 = 0
    end

    if dxi != 0 || dyi != 0 || dzi != 0
        part1 = [
            [
                part[ip][1] + dxi + dzi * part[ip][2],
                part[ip][2],
                part[ip][3] + dyi + dzi * part[ip][4],
                part[ip][4],
                part[ip][5] + dzi * sqrt(1 + part[ip][2]^2 + part[ip][4]^2),
                part[ip][6]
            ] 
            for ip in 1:n_part
        ]
    else 
        part1 = copy(part)
    end 
    # for i in 1:i_top+1
    #     coord = [part[i][1]+dxi+dzi*part[i][2], part[i][2], part[i][3]+dyi+dzi*part[i][4], part[i][4], 
    #                 part[i][5]+dzi*sqrt(1+part[i][2]^2+part[i][4]^2), part[i][6]]
    #     coord, particle_lost, s_lost, sigmaDelta2 = track_one_part(coord, n, he1, he2, e1, e2, dxf, dyf, dzf, sin_ttilt, cos_ttilt, 
    #                                                                 dcoord_etilt, psi1, csbend, Po, rho0, rho_actual, rad_coef, 
    #                                                                 isrConstant, sigmaDelta2, Fx_xy, Fy_xy, e1_kick_limit, 
    #                                                                 e2_kick_limit, expansionOrder)
    #     part[i] = coord
    # end
    new_parts_Buffer = Zygote.Buffer(part1[1], n_part, 6)
    lost_list = zeros(Int64, n_part)
    lost_list_Buffer = Zygote.Buffer(lost_list)

    for i in 1:n_part
        coord, particle_lost, s_lost, sigmaDelta2 = track_one_part(part1[i], n, he1, he2, e1, e2, dxf, dyf, dzf, sin_ttilt, cos_ttilt, 
                                            dcoord_etilt, psi1, psi2, csbend, Po, rho0, rho_actual, rad_coef, 
                                            isrConstant, sigmaDelta2, Fx_xy, Fy_xy, e1_kick_limit, 
                                            e2_kick_limit, expansionOrder)
        lost_list_Buffer[i] = particle_lost
        for j in 1:6
            new_parts_Buffer[i, j] = coord[j]
        end
    end
    new_parts = copy(new_parts_Buffer)
    new_parts_list = [new_parts[i,:] for i in 1:n_part]
    lost_list = copy(lost_list_Buffer)

    if !isnothing(sigmaDelta2) 
        sigmaDelta2 /= i_top+1
    end
    return new_parts_list

end



# function f(L, Angle)
#     CSB = CSBEND(L, Angle, 0.01, 0.02, 0.5, 0.0, 0.0 ,0.0 ,0.0 ,0.0, 0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
#             0.0, 0.0, 0.0, 1, 1, 1, -1.0, -1.0, 0.0, 1, 0, 0, 0, 1, 0, 1, 0.0, 0.0, 0.0, 0.0, 0.0, 2, 4, 0)
#     particle = [Float64[0.001, 0.0001, 0.0005, 0.0002, 0.0, 0.0], Float64[0.001, 0.0, 0.0, 0.0, 0.0, 0.0]]

#     n_part = 2
#     rout = track_through_csbend(particle, n_part, CSB, 0.0, 1000.0, 0.0, 0.0, nothing)
#     # println(rout)
#     return rout[1]
# end
# println(f(0.72, -0.1571))

# grad = Zygote.jacobian(f, 0.72, 0.1571) 
# @time grad1 = Zygote.jacobian(f, 0.72, 0.1571) 

# println(grad)