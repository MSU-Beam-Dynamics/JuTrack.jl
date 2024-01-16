include("../lattice/canonical_elements.jl")
include("../TPSA_Enzyme/TPSA_fixedmap.jl")
function exactDrift(x::CTPS, xp::CTPS, y::CTPS, yp::CTPS, z::CTPS, delta::CTPS, np::Int, len::Float64)
    x = x + xp * len
    y = y + yp * len
    z = z + len * sqrt(1.0 + xp^2 + yp^2)
    return x, xp, y, yp, z, delta

end

function SIGN(x)
    if x >= 0
        return 1.0
    else
        return -1.0
    end
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
    
    h2 = h^2
    h3 = h^3
    h4 = h^4
    h5 = h^5

    Fx_xy[1,2] = b[1]
    Fy_xy[1,1] = 1.0
    Fy_xy[2,1] = b[1]

    if nonlinear !=0 && any(b .!= 0)
        Fx_xy[2,2] = b[2]
        Fx_xy[3,2] = b[3]/2
        Fx_xy[4,2] = b[4]/6
        Fx_xy[5,2] = b[5]/24
        Fx_xy[6,2] = b[6]/120
        Fx_xy[7,2] = b[7]/720
        Fx_xy[8,2] = b[8]/5040
        Fx_xy[1,4] = (-b[3] - b[2]*h + b[1]*h2)/6
        Fx_xy[2,4] = (-b[4] - b[3]*h + 2*b[2]*h2 - 2*b[1]*h3)/6
        Fx_xy[3,4] = (-b[5] + h*(-b[4] + 3*h*(b[3] - 2*b[2]*h + 2*b[1]*h2)))/12
        Fx_xy[4,4] = (-b[6] - h*(b[5] + 4*h*(-b[4] + 3*h*(b[3] - 2*b[2]*h + 2*b[1]*h2))))/36
        Fx_xy[5,4] = (-b[7] + h*(-b[6] + 5*h*(b[5] + 4*h*(-b[4] + 3*h*(b[3] - 2*b[2]*h + 2*b[1]*h2)))))/ 144
        Fx_xy[6,4] = (-b[8] - h*(b[7] + 6*h*(-b[6] + 5*h*(b[5] + 4*h*(-b[4] + 3*h*(b[3] - 2*b[2]*h + 2*b[1]*h2))))))/720
        Fx_xy[7,4] = (h*(-b[8] + 7*h*(b[7] + 6*h*(-b[6] + 5*h*(b[5] + 4*h*(-b[4] + 3*h*(b[3] - 2*b[2]*h + 2*b[1]*h2)))))))/ 4320
        Fx_xy[8,4] = (h2*(b[8] - 7*h*(b[7] + 6*h*(-b[6] + 5*h*(b[5] + 4*h*(-b[4] + 3*h*(b[3] - 2*b[2]*h + 2*b[1]*h2)))))))/3780
        Fx_xy[1,6] = (b[5] + h*(2*b[4] - 3*h*(b[3] - b[2]*h + b[1]*h2)))/120
        Fx_xy[2,6] = (b[6] + h*(2*b[5] + h*(-5*b[4] + 9*b[3]*h - 12*b[2]*h2 + 12*b[1]*h3)))/ 120
        Fx_xy[3,6] = (b[7] + h*(2*b[6] + h*(-7*b[5] + 19*b[4]*h - 39*b[3]*h2 + 60*b[2]*h3 - 60*b[1]*h4)))/240
        Fx_xy[4,6] = (b[8] + h*(2*b[7] + 3*h*(-3*b[6] + 11*b[5]*h - 32*b[4]*h2 + 72*b[3]*h3 - 120*b[2]*h4 + 120*b[1]*h5)))/720
        Fx_xy[5,6] = (h*(2*b[8] + h*(-11*b[7] - 3*h*(-17*b[6] + 5*h*(13*b[5] + 8*h*(-5*b[4] + 3*h*(4*b[3] - 7*b[2]*h + 7*b[1]*h2)))))))/2880
        Fx_xy[6,6] = (h2*(-13*b[8] + h*(73*b[7] + 12*h*(-29*b[6] + 5*h*(23*b[5] + 2*h*(-37*b[4] + 3*h*(31*b[3] + 56*h*(-b[2] + b[1]*h))))))))/14400
        Fx_xy[1,8] = (-b[7] - 3*h*(b[6] + h*(-2*b[5] + h*(4*b[4] - 9*b[3]*h + 15*b[2]*h2 - 15*b[1]*h3))))/5040
        Fx_xy[2,8] = (-b[8] - 3*h*(b[7] + h*(-3*b[6] + 8*b[5]*h - 21*b[4]*h2 + 51*b[3]*h3 - 90*b[2]*h4 + 90*b[1]*h5)))/5040
        Fx_xy[3,8] = (h*(-b[8] + h*(4*b[7] + h*(-14*b[6] + 45*b[5]*h - 135*b[4]*h2 + 345*b[3]*h3 - 630*b[2]*h4 + 630*b[1]*h5))))/ 3360
        Fx_xy[4,8] = (h2*(5*b[8] + h*(-22*b[7] + 3*h*(29*b[6] - 5*h*(21*b[5] + 4*h*(-17*b[4] + 45*b[3]*h - 84*b[2]*h2 + 84*b[1]*h3))))))/10080
        Fx_xy[1,10] = (h*(4*b[8] - 5*h*(2*b[7] + 3*h*(-2*b[6] + 7*b[5]*h - 22*b[4]*h2 + 57*b[3]*h3 - 105*b[2]*h4 + 105*b[1]*h5))))/ 362880
        Fx_xy[2,10] = (h2*(-14*b[8] + 5*h*(10*b[7] + 3*h*(-13*b[6] + 50*b[5]*h - 167*b[4]*h2 + 447*b[3]*h3 - 840*b[2]*h4 + 840*b[1]*h5))))/362880

        Fy_xy[3,1] = b[2] / 2
        Fy_xy[4,1] = b[3] / 6
        Fy_xy[5,1] = b[4] / 24
        Fy_xy[6,1] = b[5] / 120
        Fy_xy[7,1] = b[6] / 720
        Fy_xy[8,1] = b[7] / 5040
        Fy_xy[9,1] = b[8] / 40320
        Fy_xy[1,3] = (-b[2] - b[1] * h) / 2
        Fy_xy[2,3] = (-b[3] - b[2] * h + b[1] * h2) / 2
        Fy_xy[3,3] = (-b[4] - b[3] * h + 2 * b[2] * h2 - 2 * b[1] * h3) / 4
        Fy_xy[4,3] = (-b[5] + h * (-b[4] + 3 * h * (b[3] - 2 * b[2] * h + 2 * b[1] * h2))) / 12
        Fy_xy[5,3] = (-b[6] - h * (b[5] + 4 * h * (-b[4] + 3 * h * (b[3] - 2 * b[2] * h + 2 * b[1] * h2)))) / 48
        Fy_xy[6,3] = (-b[7] + h * (-b[6] + 5 * h * (b[5] + 4 * h * (-b[4] + 3 * h * (b[3] - 2 * b[2] * h + 2 * b[1] * h2))))) / 240
        Fy_xy[7,3] = (-b[8] - h * (b[7] + 6 * h * (-b[6] + 5 * h * (b[5] + 4 * h * (-b[4] + 3 * h * (b[3] - 2 * b[2] * h + 2 * b[1] * h2)))))) / 1440
        Fy_xy[8,3] = (h * (-b[8] + 7 * h * (b[7] + 6 * h * (-b[6] + 5 * h * (b[5] + 4 * h * (-b[4] + 3 * h * (b[3] - 2 * b[2] * h + 2 * b[1] * h2))))))) / 10080
        Fy_xy[9,3] = (h2 * (b[8] - 7 * h * (b[7] + 6 * h * (-b[6] + 5 * h * (b[5] + 4 * h * (-b[4] + 3 * h * (b[3] - 2 * b[2] * h + 2 * b[1] * h2))))))) / 10080
    
        Fy_xy[1,5] = (b[4] + h*(2*b[3] + h*(-b[2] + b[1]*h)))/24
        Fy_xy[2,5] = (b[5] + h*(2*b[4] - 3*h*(b[3] - b[2]*h + b[1]*h2)))/24
        Fy_xy[3,5] = (b[6] + h*(2*b[5] + h*(-5*b[4] + 9*b[3]*h - 12*b[2]*h2 + 12*b[1]*h3)))/48
        Fy_xy[4,5] = (b[7] + h*(2*b[6] + h*(-7*b[5] + 19*b[4]*h - 39*b[3]*h2 + 60*b[2]*h3 - 60*b[1]*h4)))/144
        Fy_xy[5,5] = (b[8] + h*(2*b[7] + 3*h*(-3*b[6] + 11*b[5]*h - 32*b[4]*h2 + 72*b[3]*h3 - 120*b[2]*h4 + 120*b[1]*h5)))/576
        Fy_xy[6,5] = (h*(2*b[8] + h*(-11*b[7] - 3*h*(-17*b[6] + 5*h*(13*b[5] + 8*h*(-5*b[4] + 3*h*(4*b[3] - 7*b[2]*h + 7*b[1]*h2)))))))/2880
        Fy_xy[7,5] = (h2*(-13*b[8] + h*(73*b[7] + 12*h*(-29*b[6] + 5*h*(23*b[5] + 2*h*(-37*b[4] + 3*h*(31*b[3] + 56*h*(-b[2] + b[1]*h))))))))/17280
        Fy_xy[1,7] = (-b[6] - 3*h*(b[5] + h*(-b[4] + 2*b[3]*h - 3*b[2]*h2 + 3*b[1]*h3)))/720
        Fy_xy[2,7] = (-b[7] - 3*h*(b[6] + h*(-2*b[5] + h*(4*b[4] - 9*b[3]*h + 15*b[2]*h2 - 15*b[1]*h3))))/720
        Fy_xy[3,7] = (-b[8] - 3*h*(b[7] + h*(-3*b[6] + 8*b[5]*h - 21*b[4]*h2 + 51*b[3]*h3 - 90*b[2]*h4 + 90*b[1]*h5)))/1440
        Fy_xy[4,7] = (h*(-b[8] + h*(4*b[7] + h*(-14*b[6] + 45*b[5]*h - 135*b[4]*h2 + 345*b[3]*h3 - 630*b[2]*h4 + 630*b[1]*h5))))/ 1440
        Fy_xy[5,7] = (h2*(5*b[8] + h*(-22*b[7] + 3*h*(29*b[6] - 5*h*(21*b[5] + 4*h*(-17*b[4] + 45*b[3]*h - 84*b[2]*h2 + 84*b[1]*h3))))))/5760
        Fy_xy[1,9] = (b[8] + h*(4*b[7] + 3*h*(-2*b[6] + 6*b[5]*h - 17*b[4]*h2 + 42*b[3]*h3 - 75*b[2]*h4 + 75*b[1]*h5)))/40320
        Fy_xy[2,9] = (h*(4*b[8] - 5*h*(2*b[7] + 3*h*(-2*b[6] + 7*b[5]*h - 22*b[4]*h2 + 57*b[3]*h3 - 105*b[2]*h4 + 105*b[1]*h5))))/ 40320.;
        Fy_xy[3,9] = (h2*(-14*b[8] + 5*h*(10*b[7] + 3*h*(-13*b[6] + 50*b[5]*h - 167*b[4]*h2 + 447*b[3]*h3 - 840*b[2]*h4 + 840*b[1]*h5))))/80640
        Fy_xy[1,11] = (h2*(2*b[8] + h*(-8*b[7] - 3*h*(-11*b[6] + h*(43*b[5] + 5*h*(-29*b[4] + 3*h*(26*b[3] + 49*h*(-b[2] + b[1]*h))))))))/725760
    end
    # Fx_xy = copy(Fx_xy)
    # Fy_xy = copy(Fy_xy)
    return Fx_xy, Fy_xy, expansionOrder
end

function rotate_coordinates!(coord, angle)
    # Constants for comparing floating point numbers
    epsilon = 1e-12

    if angle == 0 || abs(abs(angle) - 2*pi) < epsilon
        return
    end

    if abs(abs(angle) - pi) < epsilon
        cos_a = -1.0
        sin_a = 0.0
    elseif abs(angle - pi/2) < epsilon
        cos_a = 0.0
        sin_a = 1.0
    elseif abs(angle + pi/2) < epsilon
        cos_a = 0.0
        sin_a = -1.0
    else
        cos_a = cos(angle)
        sin_a = sin(angle)
    end

    x = coord[1]
    xp = coord[2]
    y = coord[3]
    yp = coord[4]
    coord[1] = x * cos_a + y * sin_a
    coord[2] = xp * cos_a + yp * sin_a
    coord[3] = -x * sin_a + y * cos_a
    coord[4] = -xp * sin_a + yp * cos_a
    # new_coord = [x * cos_a + y * sin_a, xp * cos_a + yp * sin_a, -x * sin_a + y * cos_a, -xp * sin_a + yp * cos_a, coord[5], coord[6]]
    # return new_coord
    return nothing
end

function computeEtiltCentroidOffset(rho0, angle, etilt, tilt)
    # compute final offsets due to error-tilt of the magnet
    # see pages 90-93 of notebook 1 about this
        
    dcoord_etilt = zeros(6)
    if etilt == 0
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

    dcoord_etilt[1] = sqrt(q1b^2 + q2b^2)
    dcoord_etilt[2] = tan(atan(tan_alpha) - (pi/2 - angle))
    dcoord_etilt[3] = q3b
    dcoord_etilt[4] = qp3
    dcoord_etilt[5] = dz * sqrt(1 + qp3^2)
    dcoord_etilt[6] = 0
    # dcoord_etilt = [sqrt(q1b^2 + q2b^2), tan(atan(tan_alpha) - (pi/2 - angle)), q3b, qp3, dz * sqrt(1 + qp3^2), 0]

    # rotate by tilt to get into same frame as bend equations.
    rotate_coordinates!(dcoord_etilt, tilt)
    return dcoord_etilt

end

function apply_edge_effects(x::CTPS, xp::CTPS, y::CTPS, yp::CTPS, rho, n, beta, he, psi, which_edge)
    h = 1.0 / rho
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
    T431 = h2 * (2.0 * n + (which_edge == 1 ? sec2_beta : 0.0)) * tan_beta
    T432 = which_edge * h * sec2_beta

    if he != 0
        term = h / 2 * he * sec2_beta * sec_beta
        T211 += term
        T233 -= term
        T431 -= 2 * term
    end

    x0 = CTPS(x)
    xp0 = CTPS(xp)
    y0 = CTPS(y)
    yp0 = CTPS(yp)
    # x0, xp0, y0, yp0 = x, xp, y, yp
    x = x0 + T111 * x0^2 + T133 * y0^2
    xp = xp0 + R21 * x0 + T211 * x0^2 + T221 * x0 * xp0 + T233 * y0^2 + T243 * y0 * yp0
    y = y0 + T331 * x0 * y0
    yp = yp0 + R43 * y0 + T441 * yp0 * x0 + T431 * x0 * y0 + T432 * xp0 * y0
    return x, xp, y, yp
end

function computeCSBENDFields(x::CTPS{T, TPS_Dim, Max_TPS_Degree}, y::CTPS{T, TPS_Dim, Max_TPS_Degree}, 
                            Fx_xy, Fy_xy, expansionOrder1) where {T, TPS_Dim, Max_TPS_Degree}
    sumFX = CTPS(0.0, TPS_Dim, Max_TPS_Degree)
    sumFY = CTPS(0.0, TPS_Dim, Max_TPS_Degree)

    for i in 1:expansionOrder1
        for j in 2:2:(expansionOrder1 - i + 1)
            if Fx_xy[i, j] != 0
                sumFX += Fx_xy[i, j] * x^(i - 1) * y^(j - 1)
            end
        end
    end
    for i in 1:expansionOrder1
        for j in 1:2:(expansionOrder1 - i + 1)
            if Fy_xy[i, j] != 0
                sumFY += Fy_xy[i, j] * x^(i - 1) * y^(j - 1)
            end
        end
    end
    
    
    
    return sumFX, sumFY
end

function convertToDipoleCanonicalCoordinates(x, xp, y, yp, z, dp, rho)
    f = (1.0 + dp) / sqrt((1.0 + x / rho)^2 + xp^2 + yp^2)
    return xp*f, yp*f
end
function convertFromDipoleCanonicalCoordinates(x, xp, y, yp, z, dp, rho)
    f = (1.0 +x/rho)/sqrt((1.0+dp)^2-xp^2-yp^2)
    return xp*f, yp*f
end

function dipoleFringe(x, xp, y, yp, z, dp, h, inFringe, higherOrder)
    # vec = [x, qx, y, qy, s, delta]
    # x = vec[1]
    px = xp
    # y = vec[3]
    py = yp
    delta = dp
    dx = dpx = dy = dpy = ds = 0.0

    a = inFringe * h / (8.0 * (1 + delta))

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
    x += dx
    xp += dpx
    y += dy
    yp += dpy
    z += ds
    # vec[1] += dx
    # vec[2] += dpx
    # vec[3] += dy
    # vec[4] += dpy
    # vec[5] += ds

    # new_vec = [vec[1]+dx, vec[2]+dpx, vec[3]+dy, vec[4]+dpy, vec[5]+ds, vec[6]]

    return x, xp, y, yp, z, dp
end

function addRadiationKick(Qx, Qy, dPoP, sigmaDelta2,
    x, h0, Fx, Fy, ds, radCoef, synch_rad, dsISR, isrCoef,
    distributionBased, includeOpeningAngle, meanPhotonsPerMeter,
    normalizedCriticalEnergy0, Po)

    f = (1.0 + x * h0) / sqrt((1.0 + dPoP)^2 - Qx^2 - Qy^2)
    xp = Qx * f
    yp = Qy * f
    dsFactor = sqrt((1.0 + x * h0)^2 + xp^2 + yp^2)
    F2 = Fx^2 + Fy^2

    if distributionBased == 0 # now must be 0
        deltaFactor = (1.0 + dPoP)^2
        Qx /= (1.0 + dPoP)
        Qy /= (1.0 + dPoP)
        if synch_rad == 1
            dPoP -= radCoef * deltaFactor * F2 * ds * dsFactor
        end
        # if isrCoef > 0
        # # gauss_rn_lim and random_2 are not defined yet
        #     dPoP -= isrCoef * deltaFactor * F2^0.75 * sqrt(dsISR * dsFactor) * gauss_rn_lim(0.0, 1.0, 3.0, random_2())
        # end
        # if !isnothing(sigmaDelta2) # F2^1.5 is not defined yet
        #     sigmaDelta2 += (isrCoef * deltaFactor)^2 * F2^1.5 * dsISR * dsFactor
        # end
        Qx *= (1.0 + dPoP)
        Qy *= (1.0 + dPoP)
    # else
    #     F = sqrt(F2)
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

# function integrate_csbend_ord2(Qi, sigmaDelta2, s, n, rho0, p0, Fx_xy, Fy_xy, rho_actual, 
#                                 rad_coef, isrConstant, synch_rad, expansionOrder)
#     particle_lost = 0
#     s_lost = 0.0

#     if n < 1
#         error("invalid number of steps (integrate_csbend_ord4)")
#     end

#     Qf = [Qi[i] for i in 1:6]

#     dist = 0.0
#     ds = s/n
#     dsh = ds/2
#     for i in 1:n
#         if i == 1
#             # First drift
#             f = (1 + Qf[6])^2 - Qf[4]^2
#             if f <= 0
#                 particle_lost = 1
#                 s_lost = dist
#                 return Qf, particle_lost, s_lost, sigmaDelta2
#             end
#             f = sqrt(f)
#             if abs(Qf[2]/f) > 1
#                 particle_lost = 1
#                 s_lost = dist
#                 return Qf, particle_lost, s_lost, sigmaDelta2
#             end
#             sin_phi = Qf[2] / f
#             phi = asin(sin_phi)
#             sine = sin(dsh / rho0 + phi)
#             cosi = cos(dsh / rho0 + phi)
#             if cosi == 0
#                 particle_lost = 1
#                 s_lost = dist
#                 return Qf, particle_lost, s_lost, sigmaDelta2
#             end
#             tang = sine / cosi
#             cos_phi = cos(phi)
#             Qf[2] = f * sine
#             factor = (rho0+Qf[1])*cos_phi/f*(tang-sin_phi/cos_phi)
#             Qf[3] += Qf[4] * factor
#             dist += factor * (1 + Qf[6])
#             f = cos_phi / cosi
#             Qf[1] = rho0 * (f - 1) + f * Qf[1]
#         end

#         # First kick
#         x = Qf[1]
#         y = Qf[3]
#         Fx, Fy = computeCSBENDFields(x, y, Fx_xy, Fy_xy, expansionOrder+1)

#         Qf[2] += -ds * (1 + Qf[1]/rho0) * Fy / rho_actual
#         Qf[4] += ds * (1 + Qf[1]/rho0) * Fx / rho_actual
#         if synch_rad==1
#             Qx, Qy, dPoP, sigmaDelta2 = addRadiationKick(Qf[2], Qf[4], Qf[6], sigmaDelta2,
#                                         Qf[1], 1/rho0, Fx, Fy, ds, rad_coef, s/3, 0, 0, 0, 0, 0, p0)
#             Qf[2] = Qx
#             Qf[4] = Qy
#             Qf[6] = dPoP
#         end

#         # second drift
#         if i == n
#             f = (1 + Qf[6])^2 - Qf[4]^2
#             if f <= 0
#                 particle_lost = 1
#                 s_lost = dist
#                 return Qf, particle_lost, s_lost, sigmaDelta2
#             end
#             f = sqrt(f)
#             if abs(Qf[2]/f) > 1
#                 particle_lost = 1
#                 s_lost = dist
#                 return Qf, particle_lost, s_lost, sigmaDelta2
#             end
#             sin_phi = Qf[2] / f
#             phi = asin(sin_phi)
#             sine = sin(dsh / rho0 + phi)
#             cosi = cos(dsh / rho0 + phi)
#             if cosi == 0
#                 particle_lost = 1
#                 s_lost = dist
#                 return Qf, particle_lost, s_lost, sigmaDelta2
#             end
#             tang = sine / cosi
#             cos_phi = cos(phi)
#             Qf[2] = f * sine
#             factor = (rho0+Qf[1])*cos_phi/f*(tang-sin_phi/cos_phi)
#             Qf[3] += Qf[4] * factor
#             dist += factor * (1 + Qf[6])
#             f = cos_phi / cosi
#             Qf[1] = rho0 * (f - 1) + f * Qf[1]
#         else
#             f = (1 + Qf[6])^2 - Qf[4]^2
#             if f <= 0
#                 particle_lost = 1
#                 s_lost = dist
#                 return Qf, particle_lost, s_lost, sigmaDelta2
#             end
#             f = sqrt(f)
#             if abs(Qf[2]/f) > 1
#                 particle_lost = 1
#                 s_lost = dist
#                 return Qf, particle_lost, s_lost, sigmaDelta2
#             end
#             sin_phi = Qf[2] / f
#             phi = asin(sin_phi)
#             sine = sin(ds / rho0 + phi)
#             cosi = cos(ds / rho0 + phi)
#             if cosi == 0
#                 particle_lost = 1
#                 s_lost = dist
#                 return Qf, particle_lost, s_lost, sigmaDelta2
#             end
#             tang = sine / cosi
#             cos_phi = cos(phi)
#             Qf[2] = f * sine
#             factor = (rho0+Qf[1])*cos_phi/f*(tang-sin_phi/cos_phi)
#             Qf[3] += Qf[4] * factor
#             dist += factor * (1 + Qf[6])
#             f = cos_phi / cosi
#             Qf[1] = rho0 * (f - 1) + f * Qf[1]
#         end
#     end
#     Qf[5] += dist
#     return Qf, particle_lost, s_lost, sigmaDelta2
# end

function integrate_csbend_ord4(x::CTPS{T, TPS_Dim, Max_TPS_Degree}, xp::CTPS{T, TPS_Dim, Max_TPS_Degree}, 
                                y::CTPS{T, TPS_Dim, Max_TPS_Degree}, yp::CTPS{T, TPS_Dim, Max_TPS_Degree}, 
                                z::CTPS{T, TPS_Dim, Max_TPS_Degree}, dp::CTPS{T, TPS_Dim, Max_TPS_Degree}, 
                                sigmaDelta2::Float64, s::Float64, n::Int, rho0::Float64, p0::Float64, 
                                Fx_xy::Matrix{Float64}, Fy_xy::Matrix{Float64}, rho_actual::Float64, rad_coef::Float64, 
                                isrConstant::Float64, synch_rad::Int, expansionOrder::Int) where {T, TPS_Dim, Max_TPS_Degree}
    BETA = 1.25992104989487316477
    particle_lost = 0
    s_lost = 0.0

    if n < 1
        error("invalid number of steps (integrate_csbend_ord4)")
    end

    # Qf = [Qi[i] for i in 1:6]

    dist = CTPS(0.0, TPS_Dim, Max_TPS_Degree)
    s /= n
    for i in 1:n
        # First drift
        dsh = s / 2.0 / (2.0 - BETA)
        f = (dp + 1.0)^2 - yp^2
        # if f <= 0
        #     particle_lost = 1
        #     s_lost = dist
        #     return Qf, particle_lost, s_lost, sigmaDelta2
        # end
        f = sqrt(f)
        # if abs(xp/f) > 1
        #     particle_lost = 1
        #     s_lost = dist
            
        #     return Qf, particle_lost, s_lost, sigmaDelta2
        # end
        sin_phi = xp / f
        phi = asin(sin_phi)
        sine = sin(dsh / rho0 + phi)
        cosi = cos(dsh / rho0 + phi)
        # if cosi == 0
        #     particle_lost = 1
        #     s_lost = dist
            
        #     return Qf, particle_lost, s_lost, sigmaDelta2
        # end
        tang = sine / cosi
        cos_phi = cos(phi)
        xp = f * sine
        factor = (rho0+x)*cos_phi/f*(tang-sin_phi/cos_phi)
        y += yp * factor
        dist += factor * (1.0 + dp)
        f = cos_phi / cosi
        x = rho0 * (f - 1.0) + f * x

        # First kick
        ds = s / (2.0 - BETA)
        # Assuming computeCSBENDFields and addRadiationKick are defined in Julia
        # x = Qf[1]
        # y = Qf[3]
        Fx, Fy = computeCSBENDFields(x, y, Fx_xy, Fy_xy, expansionOrder+1)

        xp += -ds * (1.0 + x/rho0) * Fy / rho_actual
        yp += ds * (1.0 + x/rho0) * Fx / rho_actual
        if synch_rad == 1
            xp, yp, dp, sigmaDelta2 = addRadiationKick(xp, yp, dp, sigmaDelta2,
                                        x, 1/rho0, Fx, Fy, ds, rad_coef, synch_rad, s/3, 0.0, 0, 0, 0, 0, p0)
        end

        # second drift
        dsh = s*(1 - BETA)/(2 - BETA)/2
        f = (1.0 + dp)^2 - yp^2
        # if f <= 0
        #     particle_lost = 1
        #     s_lost = dist
            
        #     return Qf, particle_lost, s_lost, sigmaDelta2
        # end
        f = sqrt(f)
        # if abs(xp/f) > 1
        #     particle_lost = 1
        #     s_lost = dist
            
        #     return Qf, particle_lost, s_lost, sigmaDelta2
        # end
        sin_phi = xp / f
        phi = asin(sin_phi)
        sine = sin(dsh / rho0 + phi)
        cosi = cos(dsh / rho0 + phi)
        # if cosi == 0
        #     particle_lost = 1
        #     s_lost = dist
            
        #     return Qf, particle_lost, s_lost, sigmaDelta2
        # end
        tang = sine / cosi
        cos_phi = cos(phi)
        xp = f * sine
        factor = (rho0+x)*cos_phi/f*(tang-sin_phi/cos_phi)
        y += yp * factor
        dist += factor * (1.0 + dp)
        f = cos_phi / cosi
        x = rho0 * (f - 1.0) + f * x

        # second kick
        ds = -s*BETA/(2 - BETA)
        # x = Qf[1]
        # y = Qf[3]
        Fx, Fy = computeCSBENDFields(x, y, Fx_xy, Fy_xy, expansionOrder+1)

        xp += -ds * (1.0 + x/rho0) * Fy / rho_actual
        yp += ds * (1.0 + x/rho0) * Fx / rho_actual
        if synch_rad == 1
            xp, yp, dp, sigmaDelta2 = addRadiationKick(xp, yp, dp, sigmaDelta2,
                                        x, 1/rho0, Fx, Fy, ds, rad_coef, synch_rad, s/3, 0.0, 0, 0, 0, 0, p0)
            # xp = Qx
            # yp = Qy
            # dp = dPoP
        end

        # third drift
        dsh = s*(1.0 - BETA)/(2.0 - BETA)/2
        f = (1.0 + dp)^2 - yp^2
        # if f <= 0
        #     particle_lost = 1
        #     s_lost = dist
            
        #     return Qf, particle_lost, s_lost, sigmaDelta2
        # end
        f = sqrt(f)
        # if abs(Qf[2]/f) > 1
        #     particle_lost = 1
        #     s_lost = dist
            
        #     return Qf, particle_lost, s_lost, sigmaDelta2
        # end
        sin_phi = xp / f
        phi = asin(sin_phi)
        sine = sin(dsh / rho0 + phi)
        cosi = cos(dsh / rho0 + phi)
        # if cosi == 0
        #     particle_lost = 1
        #     s_lost = dist
            
        #     return Qf, particle_lost, s_lost, sigmaDelta2
        # end
        tang = sine / cosi
        cos_phi = cos(phi)
        xp = f * sine
        factor = (rho0+x)*cos_phi/f*(tang-sin_phi/cos_phi)
        y += yp * factor 
        dist += factor * (1.0 + dp)
        f = cos_phi / cosi
        x = rho0 * (f - 1.0) + f * x

        # third kick
        ds = s / (2.0 - BETA)
        xp += -ds * (1.0 + x/rho0) * Fy / rho_actual
        yp += ds * (1.0 + x/rho0) * Fx / rho_actual
        if synch_rad == 1 # rad_coef != 0 || isrConstant != 0 not work for Enzyme
            xp, yp, dp, sigmaDelta2 = addRadiationKick(xp, yp, dp, sigmaDelta2,
                                                        x, 1/rho0, Fx, Fy, ds, rad_coef, synch_rad, s/3, 0.0, 0, 0, 0, 0, p0)
            # xp = Qx
            # yp = Qy
            # dp = dPoP
        end

        # fourth drift
        dsh = s / 2 / (2 - BETA)
        f = (1.0 + dp)^2 - yp^2
        # if f <= 0
        #     particle_lost = 1
        #     s_lost = dist
            
        #     return Qf, particle_lost, s_lost, sigmaDelta2
        # end
        f = sqrt(f)
        # if abs(xp/f) > 1
        #     particle_lost = 1
        #     s_lost = dist
            
        #     return Qf, particle_lost, s_lost, sigmaDelta2
        # end
        sin_phi = xp / f
        phi = asin(sin_phi)
        sine = sin(dsh / rho0 + phi)
        cosi = cos(dsh / rho0 + phi)
        # if cosi == 0
        #     particle_lost = 1
        #     s_lost = dist
            
        #     return Qf, particle_lost, s_lost, sigmaDelta2
        # end
        tang = sine / cosi
        cos_phi = cos(phi)
        xp = f * sine
        factor = (rho0+x)*cos_phi/f*(tang-sin_phi/cos_phi)
        y += yp * factor
        dist += factor * (1.0 + dp)
        f = cos_phi / cosi
        x = rho0 * (f - 1.0) + f * x
    end
    z += dist
    
    return x, xp, y, yp, z, dp, particle_lost, s_lost, sigmaDelta2
end

function track_one_part(x0::CTPS, xp0::CTPS, y0::CTPS, yp0::CTPS, z0::CTPS, delta0::CTPS, 
                        n, he1, he2, e1, e2, dxf, dyf, dzf, sin_ttilt, cos_ttilt, dcoord_etilt, 
                        psi1, psi2, csbend, Po, rho0, rho_actual, rad_coef, isrConstant, synch_rad,
                        sigmaDelta2, Fx_xy, Fy_xy, e1_kick_limit, e2_kick_limit, expansionOrder)
    x = x0*cos_ttilt + y0*sin_ttilt
    y = y0*cos_ttilt - x0*sin_ttilt
    xp = xp0*cos_ttilt + yp0*sin_ttilt
    yp = yp0*cos_ttilt - xp0*sin_ttilt
    s = CTPS(z0)
    z = CTPS(z0)
    for i in eachindex(z.map)
        z.map[i] = 0.0
    end
    dp = CTPS(delta0)
    dp0 = CTPS(dp)

    if csbend.edge1_effects != 0
        rho = (1.0+dp)*rho_actual
        if csbend.edge_order < 2 || csbend.edge1_effects > 1
            delta_xp = tan(e1)/rho*x
                if e1_kick_limit > 0 && abs(delta_xp) > e1_kick_limit
                    delta_xp = SIGN(delta_xp)*e1_kick_limit
                end
            xp += delta_xp
            yp -= tan(e1 - psi1/(1.0+dp))/rho*y
        else
            x, xp, y, yp = apply_edge_effects(x, xp, y, yp, rho, n, e1, he1, psi1*(1.0+dp), -1)
        end
    end

    xp *= (1.0+x/rho0)
    yp *= (1.0+x/rho0)

    if csbend.edge1_effects !=0 && e1!=0.0 && synch_rad==1
        Fx, Fy = computeCSBENDFields(x, y, Fx_xy, Fy_xy, expansionOrder+1)
        dp_prime = -rad_coef * (Fx^2 + Fy^2) * (1.0 + dp)^2 * sqrt((1.0+x/rho0)^2 + xp^2 + yp^2)
        dp -= dp_prime*x*tan(e1)
    end

    px, py = convertToDipoleCanonicalCoordinates(x, xp, y, yp, 0.0, dp, rho0)

    if csbend.edge1_effects > 1
        x, px, y, py, z, dp = dipoleFringe(x, px, y, py, 0.0, dp, rho0, -1, csbend.edge1_effects-2)
    end

    particle_lost = 0
    if particle_lost == 0
        if csbend.integration_order == 4
            x, px, y, py, z, dp, particle_lost, s_lost, sigmaDelta2 = integrate_csbend_ord4(x, px, y, py, z, dp, sigmaDelta2, csbend.len, csbend.nSlices, 
                                                rho0, Po, Fx_xy, Fy_xy, rho_actual, rad_coef, isrConstant, synch_rad, expansionOrder)
        # elseif csbend.integration_order == 2
        #     x, xp, y, yp, z, dp, particle_lost, s_lost, sigmaDelta2 = integrate_csbend_ord2(x, xp, y, yp, z, dp, sigmaDelta2, csbend.len, csbend.nSlices, 
        #                                         rho0, Po, Fx_xy, Fy_xy, rho_actual, rad_coef, isrConstant, synch_rad, expansionOrder)
        else
            error("invalid integration order (track_through_csbend)")
        end
    end

    if csbend.edge2_effects > 1
        x, px, y, py, z, dp = dipoleFringe(x, px, y, py, z, dp, rho0, 1, csbend.edge2_effects-2)
    end
    xp, yp = convertFromDipoleCanonicalCoordinates(x, px, y, py, z, dp, rho0)

    ######
    # if particle_lost
    #####

    if csbend.edge2_effects !=0 && e2!=0 && synch_rad==1
        Fx, Fy = computeCSBENDFields(x, y, Fx_xy, Fy_xy, expansionOrder+1)
        dp_prime = -rad_coef * (Fx^2 + Fy^2) * (1.0 + dp)^2 * sqrt((1.0+x/rho0)^2 + xp^2 + yp^2)
        dp -= dp_prime*x*tan(e2)
    end

    # final coordinates
    if synch_rad==1
        p0 = Po*(1.0+dp0)
        beta0 = p0/sqrt(p0^2+1.0)
        p1 = Po*(1.0+dp)
        beta1 = p1/sqrt(p1^2+1.0)
        s = beta1*s/beta0 + z
    else
        s += z
    end
    # x = Qf[1]
    # xp = Qf[2]
    # y = Qf[3]
    # yp = Qf[4]
    # dp = Qf[6]

    xp /= (1.0+x/rho0)
    yp /= (1.0+x/rho0)

    if csbend.edge2_effects != 0
        rho = (1.0+dp)*rho_actual
        if csbend.edge_order < 2 || csbend.edge2_effects > 1
            delta_xp = tan(e2)/rho*x
            if e2_kick_limit > 0 && abs(delta_xp) > e2_kick_limit
                delta_xp = SIGN(delta_xp)*e2_kick_limit
            end
            xp += delta_xp
            yp -= tan(e2 - psi2/(1.0+dp))/rho*y
        else
        x, xp, y, yp = apply_edge_effects(x, xp, y, yp, rho, n, e2, he2, psi2*(1.0+dp), 1)
        end
    end

    x = x*cos_ttilt - y*sin_ttilt + dcoord_etilt[1]
    y = y*cos_ttilt + x*sin_ttilt + dcoord_etilt[3]
    xp = xp*cos_ttilt - yp*sin_ttilt + dcoord_etilt[2]
    yp = yp*cos_ttilt + xp*sin_ttilt + dcoord_etilt[4]

    x += dxf + dzf*xp
    y += dyf + dzf*yp
    s += dzf*sqrt(1.0+xp^2+yp^2)

    # coord[1] = x*cos_ttilt - y*sin_ttilt + dcoord_etilt[1]
    # coord[2] = xp*cos_ttilt - yp*sin_ttilt + dcoord_etilt[2]
    # coord[3] = y*cos_ttilt + x*sin_ttilt + dcoord_etilt[3]
    # coord[4] = yp*cos_ttilt + xp*sin_ttilt + dcoord_etilt[4]
    # coord[5] = s
    # coord[6] = dp

    # coord[1] += dxf + dzf*coord[2]
    # coord[3] += dyf + dzf*coord[4]
    # coord[5] += dzf*sqrt(1.0+coord[2]^2+coord[4]^2)


    return x, xp, y, yp, s, dp, particle_lost, s_lost, sigmaDelta2
end

function pass_TPSA(x::CTPS, xp::CTPS, y::CTPS, yp::CTPS, z::CTPS, delta::CTPS, csbend::CSBEND, Po::Float64, sigmaDelta2::Float64)
    # constants
    particleMass = 9.1093897e-31 # Electron mass (kg)
    particleCharge = 1.60217733e-19 # Electron charge (C)
    c_mks = 299792458.0 # Speed of light (m/s)
    epsilon_o = 8.854187817e-12 # Permittivity of vacuum (F/m)
    me_mev = 0.51099906 # Electron mass (MeV)
    particleRadius = particleCharge^2 / (4 * pi * epsilon_o * particleMass * c_mks^2) # Classical electron radius (m)

    n_part = 1
    largeRhoWarning = 0
    if csbend.angle == 0
        x, xp, y, yp, z, delta = exactDrift(x, xp, y, yp, z, delta, n_part, csbend.len)
        return x, xp, y, yp, z, delta
    end

    rho0 = csbend.len / csbend.angle
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
        rho0 = csbend.len / angle
        for i in 1:2:8
            b[i] = -b[i]
        end
        # b = [i % 2 == 1 ? -b[i] : b[i] for i in 1:8]
    else
        angle = csbend.angle
        e1 = csbend.e1
        e2 = csbend.e2
        etilt = csbend.etilt
        tilt = csbend.tilt
        rho0 = csbend.len / angle
    end

    if rho0 > 1e6
        largeRhoWarning = 1
        println("CSBEND Warning: large bend radius, rho0 = $rho0. Treated as drift.")
        x, xp, y, yp, z, delta = exactDrift(x, xp, y, yp, z, delta, n_part, csbend.len)
        return x, xp, y, yp, z, delta
    end

    fse = csbend.fse
    h = 1.0 / rho0
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
    Kg = 2.0*hgap * fint
    psi1 = Kg/rho_actual/cos(e1)*(1.0+sin(e1)^2)
    psi2 = Kg/rho_actual/cos(e2)*(1.0+sin(e2)^2)
    # rad_coef is d((P-Po)/Po)/ds for the on-axis, on-momentum particle, where po is the momentum of the central particle.
    if csbend.synch_rad == 1
        rad_coef = particleCharge^2*Po^3*(1.0+fse)^2/(6.0*pi*epsilon_o*c_mks^2*particleMass*rho0^2)
    else
        rad_coef = 0.0
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
        cos_ttilt = 1.0
        sin_ttilt = 0.0
    elseif abs(abs(ttilt)-pi)<1e-12
        cos_ttilt = -1.0
        sin_ttilt = 0.0
    elseif abs(abs(ttilt-pi/2))<1e-12
        cos_ttilt = 0.0
        sin_ttilt = 1.0
    elseif abs(abs(ttilt+pi/2))<1e-12
        cos_ttilt = 0.0
        sin_ttilt = -1.0
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
        sigmaDelta2 = 0.0
    end

    if dxi != 0 || dyi != 0 || dzi != 0
        x += dxi + dzi * xp
        y += dyi + dzi * yp
        z += dzi * sqrt(1.0 + xp^2 + yp^2)
    end 


    # for i in 1:n_part
        x, xp, y, yp, z, delta, particle_lost, s_lost, sigmaDelta2 = track_one_part(x, xp, y, yp, z, delta, n, he1, he2, e1, e2, dxf, dyf, dzf, sin_ttilt, cos_ttilt, 
                                            dcoord_etilt, psi1, psi2, csbend, Po, rho0, rho_actual, rad_coef, 
                                            isrConstant, csbend.synch_rad, sigmaDelta2, Fx_xy, Fy_xy, e1_kick_limit, 
                                            e2_kick_limit, expansionOrder)
        # lost_list_Buffer[i] = particle_lost
    # end

    if !isnothing(sigmaDelta2) 
        sigmaDelta2 /= i_top+1
    end
    return x, xp, y, yp, z, delta

end



# function f(xx)
#     x = CTPS(0.0, 1, 6, 3)
#     xp = CTPS(0.0, 2, 6, 3)
#     y = CTPS(0.0, 3, 6, 3)
#     yp = CTPS(0.0, 4, 6, 3)
#     z = CTPS(0.0, 5, 6, 3)
#     delta = CTPS(0.0, 6, 6, 3)
#     L = xx[1]
#     Angle = xx[2]
#     CSB = CSBEND(name="CSB",len=L, angle=Angle, nSlices=4, synch_rad=1, integration_order=4)
#     x,xp,y,yp,z,delta = pass_TPSA(x,xp,y,yp,z,delta, CSB, 19569.50762296901, 0.0)
#     # println(rout)
#     return x.map[2]
# end
# # println(f([0.72, -0.1571]))
# using Enzyme
# using BenchmarkTools
# @btime f([0.72, -0.1571])
# @btime grad = jacobian(Reverse, f, [0.72, -0.1571], Val(7))
# println(grad)

# gradient(Forward, f, [0.72, -0.1571])
