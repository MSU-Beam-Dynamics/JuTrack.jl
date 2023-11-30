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

function calculate_fringe(x, xp, y, yp, z, delta, K1, inFringe, higherOrder, linearFlag, nonlinearFactor, M1, M2)
    denom = sqrt(1 + xp^2 + yp^2)
    px = (1 + delta) * xp / denom
    py = (1 + delta) * yp / denom

    if linearFlag != 0
        x = M1[1, 1] * x + M1[1, 2] * px
        px = M1[2, 1] * x + M1[2, 2] * px
        y = M1[3, 3] * y + M1[3, 4] * py
        py = M1[4, 3] * y + M1[4, 4] * py
    # else
    #     x, px, y, py = vec[1], vec[2], vec[3], vec[4]
    end

    a = -inFringe * K1 / (12 * (1 + delta))
    dx = dpx = dy = dpy = ds = 0

    if abs(higherOrder) > 1
        # xpow = [1, 0, 0, 0, 0, 0, 0, 0, 0]
        # ypow = [1, 0, 0, 0, 0, 0, 0, 0, 0]
        # apow = [1, 0, 0, 0]
        if abs(higherOrder) > 2
            # xpow = [x^i for i in 0:8]
            # ypow = [y^i for i in 0:8]
            # apow = [a^i for i in 0:3]
            dpx = a * (12 * a * (-4 * py * x * y + px * (x^2 - 5 * y^2)) * (x^2 - y^2) 
                    - 24 * (-2 * py * x * y + px * (x^2 + y^2)) + 4 * a^2 * (-2 * py * x * y * (3 * x^4 
                    + 50 * x^2 * y^2 - 21 * y^4) + px * (x^6 + 75 * x^4 * y^2 
                    - 105 * x^2 * y^4 - 35 * y^6)) + 3 * a^3 * (-8 * py * x * y * (x^6 
                    - 27 * x^4 * y^2 + 83 * x^2 * y^4 + 7 * y^6) + px * (x^8 
                    - 108 * x^6 * y^2 + 830 * x^4 * y^4 + 196 * x^2 * y^6 
                    + 105 * y^8)))/8
            dpy = (a*(-8 * px * x * y * (6 - 6 * a * (x^2 - y^2) + a^2 * (21 * x^4 
                    - 50 * x^2 * y^2 - 3 * y^4) + 3 * a^3 * (7 * x^6 + 83 * x^4 * y^2
                    - 27 * x^2 * y^4 + y^6)) + py * (315 * a^3 * x^8 
                    + 28 * a^2 * x^6 * (5 + 21 * a * y^2) + 30 * a * x^4 * (2 + 14 * a * y^2 
                    + 83 * a^2 * y^4) + y^2 * (24 + 12 * a * y^2 - 4 * a^2 * y^4
                    + 3 * a^3 * y^6) - 12*x^2 * (-2 + 6 * a * y^2 + 25*a^2 * y^4 
                    + 27 * a^3 * y^6))))/8;
            if higherOrder > 0
                dx = (a * x * (8 * (x^2 + 3 * y^2) + 4 * a^2 * (5 * x^6 + 21 * x^4 * y^2 
                - 25 * x^2 * y^4 - y^6) + a^3 * (35 * x^8 + 84 * x^6 * y^2 
                + 498 * x^4 * y^4 - 108 * x^2 * y^6 + 3 * y^8) + 12 * a * (x^2 
                - y^2)^2)) / 8
    
                dy = (a * y * (-8 * (3 * x^2 + y^2) + 4 * a^2 * (x^6 + 25 * x^4 * y^2 
                - 21 * x^2 * y^4 - 5 * y^6) + a^3 * (3 * x^8 - 108 * x^6 * y^2 
                + 498 * x^4 * y^4 + 84 * x^2 * y^6 + 35 * y^8) + 12 * a * (x^2 
                - y^2)^2)) / 8
            end
        else
            # xpow = [x^i for i in 0:3]
            # ypow = [y^i for i in 0:3]
        
            dpx = 3 * a * (-px * x^2 + 2 * py * x * y - px * y^2)
            dpy = 3 * a * (py * x^2 - 2 * px * x * y + py * y^2)
            if higherOrder > 0
                dx = a * (x^3 + 3 * x * y^2)
                dy = a * (-3 * x^2 * y - y^3)
            end
        end
        ds = (a/(1+delta))*(3*py*y*x^2 - px*x^3 - 3*px*x*y^2 + py*y^3)
    end
    
    x = x + nonlinearFactor * dx
    px = px + nonlinearFactor * dpx
    y = y + nonlinearFactor * dy
    py = py + nonlinearFactor * dpy

    if linearFlag != 0
        x = M2[1, 1] * x + M2[1, 2] * px
        px = M2[2, 1] * x + M2[2, 2] * px
        y = M2[3, 3] * y + M2[3, 4] * py
        py = M2[4, 3] * y + M2[4, 4] * py
    end
    z = z - nonlinearFactor * ds
    # vec[5] = vec[5] - nonlinearFactor * ds

    # convert from (px, py) to (xp, yp)
    denom = (1 + delta)^2 - px^2 - py^2
    denom = sqrt(denom)
    xp = px / denom
    yp = py / denom
    # vec[2] = px / denom
    # vec[4] = py / denom

    return x, xp, y, yp, z, delta
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

function quadFringe(x, xp, y, yp, z, delta, np, K1, fringeIntM0, fringeIntP0, backtrack, 
                    inFringe, higherOrder, linearFlag, nonlinearFactor)
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
    x, xp, y, yp, z, delta = calculate_fringe(x, xp, y, yp, z, delta, K1, inFringe, higherOrder, linearFlag, nonlinearFactor, M1_new, M2_new)
    return x, xp, y, yp, z, delta
end


