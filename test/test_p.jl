include("../src/JuTrack.jl")
include("ssrf_ring.jl")
using. JuTrack
using Enzyme
using BenchmarkTools
Enzyme.API.runtimeActivity!(true)


function matrix_to_array(matrix::Matrix{Float64})
    particles = zeros(Float64, size(matrix, 1)*size(matrix, 2))
    for i in 1:size(matrix, 1)
        for j in 1:size(matrix, 2)
            particles[(i-1)*size(matrix, 2)+j] = matrix[i, j]
        end
    end
    return particles
end

# r = zeros(Float64, 100000, 6)
# r[:, 2] .= 0.001
# beam = Beam(r)
# line = ssrf(-1.063770)

# particles6 = matrix_to_array(r)
# noTarray = zeros(6)
# noRmatrix = [1.0 0.0 0.0 0.0 0.0 0.0; 
#              0.0 1.0 0.0 0.0 0.0 0.0; 
#              0.0 0.0 1.0 0.0 0.0 0.0; 
#              0.0 0.0 0.0 1.0 0.0 0.0; 
#              0.0 0.0 0.0 0.0 1.0 0.0; 
#              0.0 0.0 0.0 0.0 0.0 1.0]

# function test_track(xx)
#     r = zeros(Float64, 100000, 6)
#     r[:, 2] .= 0.001
#     beam = Beam(r)
#     D1 = DRIFT(len=0.5)
#     D2 = DRIFT(len=0.5)
#     Q1 = KQUAD(k1=xx, len=0.5)
#     Q2 = KQUAD(k1=0.5, len=0.5)
#     FODO = [Q1, D1, Q2, D2]
#     # line = ssrf(xx)
#     plinepass!(FODO, beam)
#     return beam.r[1,1]
# end
# x = [-1.063770]
# D1 = DRIFT(len=0.5)
# D2 = DRIFT(len=0.5)
# Q1 = KQUAD(k1=-0.5, len=0.5)
# Q2 = KQUAD(k1=0.5, len=0.5)
# FODO = [Q1]
# @btime begin
#     plinepass!(line, beam)
# end

function test_track(xx)
    r = zeros(Float64, 10, 6)
    r[:, 2] .= 0.001
    beam = Beam(r)
    D1 = DRIFT(name="D1", len=0.34)
    D2 = DRIFT(name="D2", len=0.12)
    D3 = DRIFT(name="D3", len=0.475)
    D4 = DRIFT(name="D4", len=0.12)
    D5 = DRIFT(name="D5", len=0.28)
    D6 = DRIFT(name="D6", len=0.59) 
    D7 = DRIFT(name="D7", len=0.36)
    D8 = DRIFT(name="D8", len=0.45)
    D9 = DRIFT(name="D9", len=0.18)
    D33 = DRIFT(name="D33", len=0.58)
    D51 = DRIFT(name="D51", len=0.12)
    D52 = DRIFT(name="D52", len=0.60)
    DL0 = DRIFT(name="DL0", len=0.60)
    DM0 = DRIFT(name="DM0", len=0.325)
    DL = DRIFT(name="DL", len=6.0)
    DM = DRIFT(name="DM", len=3.25)
    
    QL1 = KQUAD(name="QL1", len=0.32, k1=xx )
    QL2H = KQUAD(name="QL2H", len=0.29, k1=1.358860 )
    QL3 = KQUAD(name="QL3", len=0.32, k1=-1.192160 )
    QL4 = KQUAD(name="QL4", len=0.26, k1=-1.077410 )
    QL5 = KQUAD(name="QL5", len=0.32, k1=1.392450 )
    QM1 = KQUAD(name="QM1", len=0.32, k1=-1.562500 )
    QM3 = KQUAD(name="QM3", len=0.32, k1=-1.014220 )
    QM4 = KQUAD(name="QM4", len=0.26, k1=-1.366690 )
    QM5 = KQUAD(name="QM5", len=0.32, k1=1.455000 )
    Q2H = KQUAD(name="Q2H", len=0.29, k1=1.532730 )
    
    # Sextupole strength is 2 times of the AT lattice
    S1 = KSEXT(name="S1", len=0.2, k2=1.555155/0.2 )
    S2 = KSEXT(name="S2", len=0.24, k2=-3.001088/0.24 )
    S3 = KSEXT(name="S3", len=0.2, k2=2.542476/0.2)
    S4 = KSEXT(name="S4", len=0.24, k2=-2.691814/0.24 )
    S5 = KSEXT(name="S5", len=0.2, k2=3.540568/0.2)
    S6 = KSEXT(name="S6", len=0.24, k2=-4.578491/0.24)
    SD = KSEXT(name="SD", len=0.2, k2=-2.424032/0.2 )
    SF = KSEXT(name="SF", len=0.24, k2=3.436611/0.24 )
    
    # For a symmetric rectangular magnet, E1=E2=ANGLE/2
    BendingAngle = pi/20
    BD1 = SBEND(name="BD1", len=0.72, angle=BendingAngle/2, e1=BendingAngle/2, e2=0.0 )
    BD2 = SBEND(name="BD2", len=0.72, angle=BendingAngle/2, e1=0.0, e2=BendingAngle/2 )
    M1 = AbstractElement[DL,  QL1, D1,  S1,   D2,   QL2H, QL2H,  D33, 
    QL3,  D51,  S2,  D52,  BD1, BD2,  D6,   QL4,
    D7,   SD,   D8, QL5,  D9,  SF,
    D9,   QM5,  D8,  SD,   D7,  QM4,  D6,  BD1, BD2,
    D5,   QM3,  D4,  S4,   D3,  Q2H, Q2H,  D2,  S3,  D1, QM1, DM]
    linepass!(M1, beam)
    return beam.r[1,1]
end
println(test_track(-1.063770))
grad = autodiff(Forward, test_track,  Duplicated, Duplicated(-1.063770, 1.0))
println(grad)
@time grad = gradient(Forward, test_track, x, Val(1))



# using SpecialFunctions
# function f(x)
#     return erfcx(x[1])
# end
# grad = gradient(Forward, f, [0.5im], Val(1))
# println(grad)
# println(dawson(1.0im))