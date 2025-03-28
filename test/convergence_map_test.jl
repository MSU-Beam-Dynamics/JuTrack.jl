# Analyze synchro-betatron resonance using 6-D convergence map
using JuTrack
using Serialization
function ESR_map(dim, order, tunes)
    # load lattice of the ESR-EIC
    line = deserialize("src/demo/ESR/esr_main.jls")
    E0 = 18e9
    # turn on all RF cavities
    RFid = findelem(line, RFCA)
    for id in RFid
        line[id].volt = 3.78e6
        line[id].philag = 0.0
    end

    # Courant-Snyder parameters
    M66 = fastfindm66(line, 0.0)
    twiss = twiss_from_6x6(M66)
    alphax, betax, gammax, phix = twiss.horizontal
    alphay, betay, gammay, phiy = twiss.vertical
    alphaz, betaz, gammaz, phiz = twiss.longitudinal

    # construct the 6-D map. The variables of the map are Zx, Zx*, Zy, Zy*, Zz, Zz*.
    map = TPSVar6D(order)
    # reverse the Z variables to phase space variables
    x, px, y, py, z, pz = get_variables(map, betax=betax, alphax=alphax, betay=betay, alphay=alphay, betaz=betaz, alphaz=alphaz)

    # compute transfer map using TPSA.
    rin = [x, px, y, py, z, pz]
    linepass_TPSA!(line, rin, E0=E0)

    x_map = rin[1]
    px_map = rin[2]
    y_map = rin[3]
    py_map = rin[4]
    z_map = rin[5]
    pz_map = rin[6]

    z1 = x_map/sqrt(betax) - 1.0im*alphax*x_map/sqrt(betax) - 1.0im*px_map*sqrt(betax)
    z1c = x_map/sqrt(betax) + 1.0im*alphax*x_map/sqrt(betax) + 1.0im*px_map*sqrt(betax)
    z2 = y_map/sqrt(betay) - 1.0im*alphay*y_map/sqrt(betay) - 1.0im*py_map*sqrt(betay)
    z2c = y_map/sqrt(betay) + 1.0im*alphay*y_map/sqrt(betay) + 1.0im*py_map*sqrt(betay)
    z3 = z_map/sqrt(betaz) - 1.0im*alphaz*z_map/sqrt(betaz) - 1.0im*pz_map*sqrt(betaz)
    z3c = z_map/sqrt(betaz) + 1.0im*alphaz*z_map/sqrt(betaz) + 1.0im*pz_map*sqrt(betaz)

    z1mapfunc = evaluate(z1)
    z1cmapfunc = evaluate(z1c)
    z2mapfunc = evaluate(z2)
    z2cmapfunc = evaluate(z2c)
    z3mapfunc = evaluate(z3)
    z3cmapfunc = evaluate(z3c)

    construct_sqr_matrix(map, [z1, z1c, z2, z2c, z3, z3c])
    return map, function(zs)
        return [
            z1mapfunc(zs),
            z1cmapfunc(zs),
            z2mapfunc(zs),
            z2cmapfunc(zs),
            z3mapfunc(zs),
            z3cmapfunc(zs)
        ]
    end
end

tune_guesses = [0.08, 0.14, 0.052]
CMscan(ESR_map, 3, 3, tune_guesses, 
    -0.0015, 0.0015, 0.0002, 0.0002, -0.015, 0.015,
    100, 1, 100,
    0.0, 0.0, 0.0)
