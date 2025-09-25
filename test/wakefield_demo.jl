# Using Enzyme to calculate the gradient of a wakefield function
# The implementation of fast TPSA is still under development, so we use the Enzyme AD for now.
using Pkg
Pkg.activate("."); Pkg.instantiate() # change "." to your path of JuTrack.jl
using JuTrack

function f(x1)
    RLCwake = LongitudinalRLCWake(freq=x1, Rshunt=DTPSAD(5.5e3), Q0=DTPSAD(3.0))
    D1 = DRIFT(len=DTPSAD(0.1))
    D2 = DRIFT(len=DTPSAD(0.1))

    ebeam = Beam(DTPSAD.(randn_approx(500, 6)) / 1e4)
    histogram1DinZ!(ebeam)
    line = [D1,RLCwake]
    linepass!(line, ebeam)
    get_emittance!(ebeam)
    # scatter(ebeam.r[:, 1], ebeam.r[:, 2], markersize=0.1, label="x-y")
    return ebeam.emittance[1], ebeam.emittance[2], ebeam.emittance[3]
end
# the grad here may be physically meaningless. it is just a test
grads, results = Jacobian(f, [180e9], true)
println("Results: ", results)
println("Gradient: ", grads)

#############################################################################
# Legacy code for testing the wakefield function with Enzyme AD
# Uncomment the following lines to run the legacy code
#############################################################################
# using JuTrack
# function f(x1)
#     RLCwake = LongitudinalRLCWake(freq=x1, Rshunt=5.5e3, Q0=3.0)
#     D1 = DRIFT(len=0.1)
#     D2 = DRIFT(len=0.1)

#     ebeam = Beam(randn_approx(500, 6) / 1e4)
#     histogram1DinZ!(ebeam)
#     line = [D1,RLCwake]
#     linepass!(line, ebeam)
#     get_emittance!(ebeam)
#     # scatter(ebeam.r[:, 1], ebeam.r[:, 2], markersize=0.1, label="x-y")
#     return ebeam.emittance
# end

# println(f(180e9))
# # the grad here may be physically meaningless. it is just a test
# grad = autodiff(Forward, f,  Duplicated(180e9, 1.0)) 
# println(grad)


