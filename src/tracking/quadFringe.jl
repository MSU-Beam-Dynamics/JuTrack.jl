function quadPartialFringeMatrix(K1::Float64, inFringe::Int, fringeInt::Array{Float64, 1}, part::Int)
    K1sqr = K1^2

    if part == 1
        J1x = inFringe * (K1 * fringeInt[2] - 2 * K1sqr * fringeInt[4] / 3 - K1sqr * fringeInt[1] * fringeInt[3] / 2)
        J2x = inFringe * K1 * fringeInt[3]
        J3x = inFringe * K1sqr * (fringeInt[3] + fringeInt[5] + fringeInt[1] * fringeInt[2])

        K1 = -K1
        J1y = inFringe * (K1 * fringeInt[2] - 2 * K1sqr * fringeInt[4] / 3 - K1sqr * fringeInt[1] * fringeInt[3] / 2)
        J2y = -J2x
        J3y = J3x
    else
        J1x = inFringe * (K1 * fringeInt[2] + K1sqr * fringeInt[1] * fringeInt[3] / 2)
        J2x = inFringe * K1 * fringeInt[3]
        J3x = inFringe * K1sqr * (fringeInt[5] - fringeInt[1] * fringeInt[2])

        K1 = -K1
        J1y = inFringe * (K1 * fringeInt[2] + K1sqr * fringeInt[1] * fringeInt[3] / 2)
        J2y = -J2x
        J3y = J3x
    end

    expJ1x = exp(J1x)
    expJ1y = exp(J1y)

    M = [
        expJ1x        J2x/expJ1x   0             0             0 0;
        expJ1x*J3x    (1+J2x*J3x)/expJ1x 0             0             0 0;
        0             0             expJ1y        J2y/expJ1y   0 0;
        0             0             expJ1y*J3y    (1+J2y*J3y)/expJ1y 0 0;
        0             0             0             0             1 0;
        0             0             0             0             0 1
    ]

    return M
end

function calculate_fringe(vec, K1, inFringe, higherOrder, linearFlag, nonlinearFactor, M1, M2)
    delta = vec[6]

    # Convert from (xp, yp) to (px, py)
    xp = vec[2]
    yp = vec[4]
    denom = sqrt(1 + xp^2 + yp^2)
    # vec[2] = (1 + delta) * xp / denom
    # vec[4] = (1 + delta) * yp / denom
    
    # redefine vec 
    vec = [vec[1], (1 + delta) * xp / denom, vec[3], (1 + delta) * yp / denom, vec[5], vec[6]]

    if linearFlag != 0
        x = M1[1, 1] * vec[1] + M1[1, 2] * vec[2]
        px = M1[2, 1] * vec[1] + M1[2, 2] * vec[2]
        y = M1[3, 3] * vec[3] + M1[3, 4] * vec[4]
        py = M1[4, 3] * vec[3] + M1[4, 4] * vec[4]
    else
        x, px, y, py = vec[1], vec[2], vec[3], vec[4]
    end

    a = -inFringe * K1 / (12 * (1 + delta))
    dx = dpx = dy = dpy = ds = 0

    if abs(higherOrder) > 1
        xpow = [1, 0, 0, 0, 0, 0, 0, 0, 0]
        ypow = [1, 0, 0, 0, 0, 0, 0, 0, 0]
        apow = [1, 0, 0, 0]
        if abs(higherOrder) > 2
            xpow = [x^i for i in 0:8]
            ypow = [y^i for i in 0:8]
            apow = [a^i for i in 0:3]
            dpx = a * (12 * a * (-4 * py * x * y + px * (xpow[3] - 5 * ypow[3])) * (xpow[3] - ypow[3]) 
                    - 24 * (-2 * py * x * y + px * (xpow[3] + ypow[3])) + 4 * apow[3] * (-2 * py * x * y * (3 * xpow[5] 
                    + 50 * xpow[3] * ypow[3] - 21 * ypow[5]) + px * (xpow[7] + 75 * xpow[5] * ypow[3] 
                    - 105 * xpow[3] * ypow[5] - 35 * ypow[7])) + 3 * apow[4] * (-8 * py * x * y * (xpow[7] 
                    - 27 * xpow[5] * ypow[3] + 83 * xpow[3] * ypow[5] + 7 * ypow[7]) + px * (xpow[9] 
                    - 108 * xpow[7] * ypow[3] + 830 * xpow[5] * ypow[5] + 196 * xpow[3] * ypow[7] 
                    + 105 * ypow[9])))/8
            dpy = (a*(-8 * px * x * y * (6 - 6 * a * (xpow[3] - ypow[3]) + apow[3] * (21 * xpow[5] 
                    - 50 * xpow[3] * ypow[3] - 3 * ypow[5]) + 3 * apow[4] * (7 * xpow[7] + 83 * xpow[5] * ypow[3]
                    - 27 * xpow[3] * ypow[5] + ypow[7])) + py * (315 * apow[4] * xpow[9] 
                    + 28 * apow[3] * xpow[7] * (5 + 21 * a * ypow[3]) + 30 * a * xpow[5] * (2 + 14 * a * ypow[3] 
                    + 83 * apow[3] * ypow[5]) + ypow[3] * (24 + 12 * a * ypow[3] - 4 * apow[3] * ypow[5]
                    + 3 * apow[4] * ypow[7]) - 12*xpow[3] * (-2 + 6 * a * ypow[3] + 25*apow[3] * ypow[5] 
                    + 27 * apow[4] * ypow[7]))))/8;
            if higherOrder > 0
                dx = (a * x * (8 * (xpow[3] + 3 * ypow[3]) + 4 * apow[3] * (5 * xpow[7] + 21 * xpow[5] * ypow[3] 
                - 25 * xpow[3] * ypow[5] - ypow[7]) + apow[4] * (35 * xpow[9] + 84 * xpow[7] * ypow[3] 
                + 498 * xpow[5] * ypow[5] - 108 * xpow[3] * ypow[7] + 3 * ypow[9]) + 12 * a * (xpow[3] 
                - ypow[3])^2)) / 8
    
                dy = (a * y * (-8 * (3 * xpow[3] + ypow[3]) + 4 * apow[3] * (xpow[7] + 25 * xpow[5] * ypow[3] 
                - 21 * xpow[3] * ypow[5] - 5 * ypow[7]) + apow[4] * (3 * xpow[9] - 108 * xpow[7] * ypow[3] 
                + 498 * xpow[5] * ypow[5] + 84 * xpow[3] * ypow[7] + 35 * ypow[9]) + 12 * a * (xpow[3] 
                - ypow[3])^2)) / 8
            end
        else
            xpow = [x^i for i in 0:3]
            ypow = [y^i for i in 0:3]
        
            dpx = 3 * a * (-px * xpow[3] + 2 * py * x * y - px * ypow[3])
            dpy = 3 * a * (py * xpow[3] - 2 * px * x * y + py * ypow[3])
            if higherOrder > 0
                dx = a * (xpow[4] + 3 * x * ypow[3])
                dy = a * (-3 * xpow[3] * y - ypow[4])
            end
        end
        ds = (a/(1+delta))*(3*py*y*xpow[3] - px*xpow[4] - 3*px*x*ypow[3] + py*ypow[4])
    end
    
    x = x + nonlinearFactor * dx
    px = px + nonlinearFactor * dpx
    y = y + nonlinearFactor * dy
    py = py + nonlinearFactor * dpy

    if linearFlag != 0
        vec = [M2[1, 1]*x + M2[1, 2]*px, M2[2, 1] * x + M2[2, 2] * px, 
                M2[3, 3]*y + M2[3, 4]*py, M2[4, 3] * y + M2[4, 4] * py, 
                vec[5], vec[6]]
        # vec[1] = M2[1, 1] * x + M2[1, 2] * px
        # vec[2] = M2[2, 1] * x + M2[2, 2] * px
        # vec[3] = M2[3, 3] * y + M2[3, 4] * py
        # vec[4] = M2[4, 3] * y + M2[4, 4] * py
    else
        vec = [x, px, y, py, vec[5], vec[6]]
        # vec[1], vec[2], vec[3], vec[4] = x, px, y, py
    end
    vec = [vec[1], vec[2], vec[3], vec[4], vec[5] - nonlinearFactor * ds, vec[6]]
    # vec[5] = vec[5] - nonlinearFactor * ds

    # convert from (px, py) to (xp, yp)
    px = vec[2]
    py = vec[4]
    denom = (1 + delta)^2 - px^2 - py^2
    if denom > 0
        denom = sqrt(denom)
        vec = [vec[1], px / denom, vec[3], py / denom, vec[5], vec[6]]
        # vec[2] = px / denom
        # vec[4] = py / denom
    else
        DBL_MAX = 1.7976931348623158e+308 # max value 
        vec = [vec[1], DBL_MAX, vec[3], DBL_MAX, vec[5], vec[6]]
        # vec[2] = DBL_MAX
        # vec[4] = DBL_MAX
    end
    return vec
end

function swap_matrix(M1, M2)
    M2_new = [
        M1[2, 2] M1[1, 2] M1[1, 3] M1[1, 4] M1[1, 5] M1[1, 6];
        M1[2, 1] M1[1, 1] M1[2, 3] M1[2, 4] M1[2, 5] M1[2, 6];
        M1[3, 1] M1[3, 2] M1[4, 4] M1[3, 4] M1[3, 5] M1[3, 6];
        M1[4, 1] M1[4, 2] M1[4, 3] M1[3, 3] M1[4, 5] M1[4, 6];
        M1[5, 1] M1[5, 2] M1[5, 3] M1[5, 4] M1[5, 5] M1[5, 6];
        M1[6, 1] M1[6, 2] M1[6, 3] M1[6, 4] M1[6, 5] M1[6, 6]
        ]

    M1_new = [
        M2[2, 2] M2[1, 2] M2[1, 3] M2[1, 4] M2[1, 5] M2[1, 6];
        M2[2, 1] M2[1, 1] M2[2, 3] M2[2, 4] M2[2, 5] M2[2, 6];
        M2[3, 1] M2[3, 2] M2[4, 4] M2[3, 4] M2[3, 5] M2[3, 6];
        M2[4, 1] M2[4, 2] M2[4, 3] M2[3, 3] M2[4, 5] M2[4, 6];
        M2[5, 1] M2[5, 2] M2[5, 3] M2[5, 4] M2[5, 5] M2[5, 6];
        M2[6, 1] M2[6, 2] M2[6, 3] M2[6, 4] M2[6, 5] M2[6, 6]
        ]
    return M1_new, M2_new
end

function quadFringe(coord, np, K1, fringeIntM0, fringeIntP0, backtrack, inFringe, higherOrder, linearFlag, nonlinearFactor)
    if linearFlag != 0
        M1 = quadPartialFringeMatrix(K1, 1, fringeIntM0, 1)
        M2 = quadPartialFringeMatrix(K1, 1, fringeIntP0, 2)

        
        if inFringe * (backtrack ? -1 : 1) == -1
            # M1[1, 1], M1[2, 2] = M1[2, 2], M1[1, 1]
            # M1[3, 3], M1[4, 4] = M1[4, 4], M1[3, 3]
            # M2[1, 1], M2[2, 2] = M2[2, 2], M2[1, 1]
            # M2[3, 3], M2[4, 4] = M2[4, 4], M2[3, 3]
            # M1, M2 = M2, M1
            
            # avoid muatating M1 and M2
            M1_new, M2_new = swap_matrix(M1, M2)
        else
            M1_new, M2_new = M1, M2
        end

        # backtrack is always false. Back track is not in use for now to avoid confusion
        # if backtrack
        #     # Invert the matrices
        #     M1[1, 1], M1[2, 2] = M1[2, 2], M1[1, 1]
        #     M1[3, 3], M1[4, 4] = M1[4, 4], M1[3, 3]
        #     M1[1, 2] *= -1
        #     M1[2, 1] *= -1
        #     M1[3, 4] *= -1
        #     M1[4, 3] *= -1

        #     M2[1, 1], M2[2, 2] = M2[2, 2], M2[1, 1]
        #     M2[3, 3], M2[4, 4] = M2[4, 4], M2[3, 3]
        #     M2[1, 2] *= -1
        #     M2[2, 1] *= -1
        #     M2[3, 4] *= -1
        #     M2[4, 3] *= -1
        # end
    end
    # for ip in 1:np
        # vec = calculate_fringe(coord[ip], K1, inFringe, higherOrder, linearFlag, nonlinearFactor, M1, M2)
        # coord[ip] = vec
    # end
    coord_new = [calculate_fringe(coord[ip], K1, inFringe, higherOrder, linearFlag, nonlinearFactor, M1_new, M2_new) for ip in 1:np]

    return coord_new
end


