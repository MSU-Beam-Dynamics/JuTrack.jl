# Example usage of Frequency Map Analysis (FMA) and plotting the results
using Pkg
Pkg.activate("."); Pkg.instantiate() # change "." to your path of JuTrack.jl
using JuTrack
include("../src/utils/fma.jl")

function testring(rad)
    CENOFSTR01 = MARKER(name="CENOFSTR01")
    CENOFSTR02 = MARKER(name="CENOFSTR02")
    CENOFSTR03 = MARKER(name="CENOFSTR03")
    CENOFSTR04 = MARKER(name="CENOFSTR04")
    CENOFSTR05 = MARKER(name="CENOFSTR05")
    CENOFSTR06 = MARKER(name="CENOFSTR06")
    CENOFSTR07 = MARKER(name="CENOFSTR07")
    CENOFSTR08 = MARKER(name="CENOFSTR08")
    CENOFSTR09 = MARKER(name="CENOFSTR09")
    CENOFSTR10 = MARKER(name="CENOFSTR10")
    CENOFSTR11 = MARKER(name="CENOFSTR11")
    CENOFSTR12 = MARKER(name="CENOFSTR12")
    D11 = DRIFT(len=0.5, name="D11")
    D11A = DRIFT(len=0.535, name="D11A")
    DX = DRIFT(len=0.0375, name="DX")
    D12 = DRIFT(len=0.075, name="D12")
    D15 = DRIFT(len=0.225, name="D15")
    QF1 = KQUAD(len=0.18, k1=0.1384045189E+02, NumIntSteps=10, rad=rad, name="QF1") # 13.80301
    QD1 = KQUAD(len=0.14, k1=-0.1377444068E+02, NumIntSteps=10, rad=rad, name="QD1") # -13.82396
    QF2 = KQUAD(len=0.19, k1=10.19367, NumIntSteps=10, rad=rad, name="QF2")
    QF3 = KQUAD(len=0.115, k1=10.94517, NumIntSteps=10, rad=rad, name="QF3")
    QF4 = KQUAD(len=0.305, k1=15.32064, NumIntSteps=10, rad=rad, name="QF4")
    QF5 = KQUAD(len=0.305, k1=15.8, NumIntSteps=10, rad=rad, name="QF5")
    QF6 = KQUAD(len=0.305, k1=15.68564, NumIntSteps=10, rad=rad, name="QF6")
    SHH = KSEXT(len=0.075, k2=3.514648, NumIntSteps=10, rad=rad, name="SHH")
    SHH2 = KSEXT(len=0.075, k2=-929.5771999999999, NumIntSteps=10, rad=rad, name="SHH2")
    SD = KSEXT(len=0.28, k2=-1367.68410852, NumIntSteps=10, rad=rad, name="SD")
    SF = KSEXT(len=0.28, k2=1610.76982434, NumIntSteps=10, rad=rad, name="SF")
    BEND1 = SBEND(len=0.34, angle=0.05817765, PolynomB=[0.0,-2.827967,0.0,0.0], 
                NumIntSteps=10, rad=rad, MaxOrder=1, name="BEND1")
    BEND2 = SBEND(len=0.5, angle=0.05817765, PolynomB=[0.0,-7.057813,0.0,0.0],
                NumIntSteps=10, rad=rad, MaxOrder=1, name="BEND2")
    BEND3 = SBEND(len=0.5, angle=0.05817765, PolynomB=[0.0,-7.057813,0.0,0.0],
                NumIntSteps=10, rad=rad, MaxOrder=1, name="BEND3")
    STR_A = [D11A, D11, D11, D11, D11]
    STR_B = [D11, D11, D11, D11, D11A]
    ARC = [DX, SHH, D12, QF1, DX, DX, QD1, DX, DX, SHH2, D12, BEND1, 
        DX, DX, SD, D12, QF2, D12, SF, DX, DX, QF3, D15, BEND2, 
        DX, DX, QF4, DX, DX, BEND3, DX, DX, QF5, D12, BEND3, 
        DX, DX, QF6, DX, DX, BEND3, DX, DX, QF6, DX, DX, BEND3, 
        DX, DX, QF5, D12, BEND3, D12, QF4, DX, DX, BEND2, 
        D15, QF3, D12, SF, DX, DX, QF2, D12, SD, DX, DX, BEND1, 
        DX, DX, SHH2, DX, DX, QD1, DX, DX, QF1, D12, SHH, DX]
    RING = [CENOFSTR01, STR_B..., ARC..., STR_A..., CENOFSTR02, STR_B..., ARC..., STR_A..., CENOFSTR03, 
        STR_B..., ARC..., STR_A..., CENOFSTR04, STR_B..., ARC..., STR_A..., CENOFSTR05, STR_B...,
        ARC..., STR_A..., CENOFSTR06, STR_B..., ARC..., STR_A..., CENOFSTR07, STR_B..., ARC..., STR_A..., 
        CENOFSTR08, STR_B..., ARC..., STR_A..., CENOFSTR09, STR_B..., ARC..., STR_A..., CENOFSTR10, 
        STR_B..., ARC..., STR_A..., CENOFSTR11, STR_B..., ARC..., STR_A..., CENOFSTR12, STR_B..., ARC..., STR_A...]
end

# create the ring for FMA test
RING = testring(0)

# get the FMA results
rows = FMA(RING, 1024)

# plot the FMA results
plot_fma(rows)
