using JuTrack


function f(x1)
    RLCwake = LongitudinalRLCWake(x1, 5.5e3, 3.0)
    D1 = DRIFT(len=0.1)
    D2 = DRIFT(len=0.1)

    ebeam = Beam(randn(500, 6) / 1e4)
    histogram1DinZ!(ebeam)
    line = [D1,RLCwake]
    plinepass!(line, ebeam)
    get_emittance!(ebeam)
    # scatter(ebeam.r[:, 1], ebeam.r[:, 2], markersize=0.1, label="x-y")
    return ebeam.emittance
end

println(f(180e9))
grad = autodiff(Forward, f, Duplicated,  Duplicated(180e9, 1.0)) 
println(grad)


