include("../TPSA_Enzyme/TPSA_fixedmap.jl")
using Enzyme


function ^(ctps::CTPS{T, TPS_Dim, Max_TPS_Degree}, b::Int) where {T, TPS_Dim, Max_TPS_Degree}
    return pow(ctps, b)
end

function f3(xx)
    x = CTPS(xx[1], 1, 6, 3)
    xp = CTPS(xx[2], 2, 6, 3)
    y = CTPS(0.0, 3, 6, 3)
    yp = CTPS(0.0, 4, 6, 3)
    delta = CTPS(0.0, 5, 6, 3)
    z = CTPS(0.0, 6, 6, 3)

    denom1 =  xp^2
    denom2 =  yp^2

    # (1.0 + delta)^2 

    return x.map[2]
end
grad = Enzyme.gradient(Forward, f3, [1.0, 1.0], Val(2))
print(grad)