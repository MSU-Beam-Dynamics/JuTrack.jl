# using Enzyme
# Enzyme.API.runtimeActivity!(true)
include("../src/JuTrack.jl")
using .JuTrack
# not work for current version
# function erfcx_AD(z)
#     return [erfcx(z[1])]
# end

# function forward(func::Const{typeof(erfcx_AD)}, ::Type{<:Duplicated}, z::Duplicated)
#     println("Using custom rule for forward mode on erfcx function!")
#     ret = func.val(z.val)  
#     z_derivative = -2.0/sqrt(pi) .+ 2.0 .* z.val .* ret 
#     z.dval .= z_derivative .* z.dval
#     return Duplicated(ret, z.dval)
# end

# x = [5.0 ]
# dx = [1.0 ]
# result = autodiff(Forward, erfcx_AD, Duplicated, Duplicated(x, dx))
# println("Forward mode result: ", result)

# function g(x, y, simga)
#     term1=erfcx_AD([-1.0im*(x+1.0im*y)/simga])
#     return term1[1]
# end
# result = autodiff(Forward, g, Duplicated, Duplicated(0.1, 1.0), Duplicated(0.2, 1.0), Duplicated(0.3, 1.0))
# println("Forward mode result: ", result)

function f(x)
    vbase=3.42*8.5e6
    ϕs=10.0
    vact=vbase/cos(ϕs*π/180.0)
    freq = 591e6
    mainRFe=AccelCavity(freq, vact, 7560.0, π-ϕs*π/180.0)
    tunex, tuney=50.08, 44.14
    αc=3.42/tunex/tunex
    lmap=LongitudinalRFMap(αc, mainRFe)
    opIPp=optics4DUC(0.8, 0.0, 0.072, 0.0)
    opIPe=optics4DUC(0.45,0.0,0.056,0.0)
    
    pstrong=StrongGaussianBeam(1.0, m_p, 1.0, Int(0.688e11), 275e9,  opIPp, [x, 8.5e-6, 0.06], 9)
    initilize_zslice!(pstrong, :gaussian, :evennpar, 7.0)
    crab_ratio=0.33
    overcrab=1.0
    pcrab1st = easyCRABCAVITY(freq=197.0e6, halfthetac=overcrab*12.5e-3*(1+crab_ratio))
    pcrab2nd = easyCRABCAVITY(freq=197.0e6*2.0, halfthetac=-overcrab*12.5e-3*crab_ratio)
    crab_crossing_setup!(pstrong, 12.5e-3, pcrab1st, pcrab2nd)
    
    ebeam=Beam(zeros(5000, 6), np = Int(1.72e11*3), energy = 10e9, emittance=[20e-9, 1.3e-9, 1.36e-4])
    initilize_6DGaussiandist!(ebeam, opIPe, lmap)
    
    rin = matrix_to_array(ebeam.r)
    lumis=pass_lumi_P!(pstrong, rin, ebeam.nmacro, ebeam)
end


vbase=3.42*8.5e6
ϕs=10.0
vact=vbase/cos(ϕs*π/180.0)
mainRFe=AccelCavity(591e6, vact, 7560.0, π-ϕs*π/180.0)
tunex, tuney=50.08, 44.14
αc=3.42/tunex/tunex
lmap=LongitudinalRFMap(αc, mainRFe)
opIPp=optics4DUC(0.8, 0.0, 0.072, 0.0)
opIPe=optics4DUC(0.45,0.0,0.056,0.0)
function f1(x)

    global lmap, opIPe, opIPp
    pstrong=StrongGaussianBeam(1.0, m_p, 1.0, Int(0.688e11), 275e9,  opIPp, [x, 8.5e-6, 0.06], 9)
    initilize_zslice!(pstrong, :gaussian, :evennpar, 7.0)
    crab_ratio=0.33
    overcrab=1.0
    pcrab1st = easyCRABCAVITY(freq=197.0e6, halfthetac=overcrab*12.5e-3*(1+crab_ratio))
    pcrab2nd = easyCRABCAVITY(freq=197.0e6*2.0, halfthetac=-overcrab*12.5e-3*crab_ratio)
    crab_crossing_setup!(pstrong, 12.5e-3, pcrab1st, pcrab2nd)
    
    ebeam=Beam(zeros(5000, 6), np = Int(1.72e11*3), energy = 10e9, emittance=[20e-9, 1.3e-9, 1.36e-4])
    initilize_6DGaussiandist!(ebeam, opIPe, lmap)
    
    rin = matrix_to_array(ebeam.r)
    # pass_P!(pstrong, rin, ebeam.nmacro, ebeam)
    linepass!([pstrong], ebeam)
    get_emittance!(ebeam)
    return ebeam.emittance[1]

end

# println(f(95e-6))
println(f1(95e-6))
@time g = autodiff(Forward, f1, Duplicated, Duplicated(95e-6, 1.0))
print(g)