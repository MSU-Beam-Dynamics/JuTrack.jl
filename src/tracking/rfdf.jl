const c_mks = 2.99792458e8
function track_through_rf_deflector(final, rf_param, initial, np, pc_central)
    omega = 2 * pi * rf_param.freq
    k = omega / c_mks
end
