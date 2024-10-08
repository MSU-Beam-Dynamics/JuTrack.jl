using JuTrack
using Serialization
function spear3()
    RF = RFCA(name="RF", len=0.0000000000, volt=0.0000000000, h=372.0000000000, freq=4.7630000575e+08, energy=3.0000000000e+09)
    DA1A = DRIFT(name="DA1A", len=3.6792386000)
    BPM = MARKER(name="BPM")
    DA1B = DRIFT(name="DA1B", len=0.1240666500)
    QDX = KQUAD(name="QDX", len=0.3533895000, k1=-1.3864672452)
    DA2A = DRIFT(name="DA2A", len=0.1153052500)
    COR = CORRECTOR(name="COR", len=0.1500000000, xkick=0.0000000000, ykick=0.0000000000)
    DA2B = DRIFT(name="DA2B", len=0.1177344500)
    QFX = KQUAD(name="QFX", len=0.6105311000, k1=1.5731962784)
    DA3A = DRIFT(name="DA3A", len=0.2088992500)
    DA3B = DRIFT(name="DA3B", len=0.0541404500)
    QDY = KQUAD(name="QDY", len=0.3533895000, k1=-0.4606409306)
    DM4 = DRIFT(name="DM4", len=0.2158457200)
    B34 = SBEND(name="B34", len=1.1432900000, angle=0.1385996759, e1=0.0692998379, e2=6.9299837947e-02, PolynomB=[0.0, -3.1537858000e-01, 0.0, 0.0])
    DA5A = DRIFT(name="DA5A", len=0.1139774700)
    DA5B = DRIFT(name="DA5B", len=0.1085630000)
    SDM = KSEXT(name="SDM", len=0.2100000000, k2=-17.0000000000)
    DA6A = DRIFT(name="DA6A", len=0.1266000000)
    DA6B = DRIFT(name="DA6B", len=0.9047682800)
    SFM = KSEXT(name="SFM", len=0.2100000000, k2=15.0000000000)
    DA7A = DRIFT(name="DA7A", len=0.1106966000)
    DA7B = DRIFT(name="DA7B", len=0.0631132500)
    QFY = KQUAD(name="QFY", len=0.5123803000, k1=1.4814937098)
    DM7 = DRIFT(name="DM7", len=0.1738098500)
    DA6C = DRIFT(name="DA6C", len=0.0960000000)
    DA6D = DRIFT(name="DA6D", len=0.9353700000)
    DA5C = DRIFT(name="DA5C", len=0.0518450000)
    DA5D = DRIFT(name="DA5D", len=0.1706954700)
    DA8A = DRIFT(name="DA8A", len=0.3373594700)
    DA8B = DRIFT(name="DA8B", len=0.1284862500)
    QDZ = KQUAD(name="QDZ", len=0.3533895000, k1=-0.8782239377)
    DA9A = DRIFT(name="DA9A", len=0.1093052500)
    DA9B = DRIFT(name="DA9B", len=0.1373052500)
    QFZ = KQUAD(name="QFZ", len=0.3533895000, k1=1.4279020070)
    DA10A = DRIFT(name="DA10A", len=0.1239396500)
    DA10B = DRIFT(name="DA10B", len=3.1459370000)
    DC1A = DRIFT(name="DC1A", len=1.4059340000)
    DC1B = DRIFT(name="DC1B", len=0.1240412500)
    QF = KQUAD(name="QF", len=0.3533895000, k1=1.7686729041)
    DC2A = DRIFT(name="DC2A", len=0.1157652500)
    DC2B = DRIFT(name="DC2B", len=0.1158104500)
    QD = KQUAD(name="QD", len=0.1634591000, k1=-1.5424742304)
    DC3A = DRIFT(name="DC3A", len=0.0532206500)
    DC3B = DRIFT(name="DC3B", len=0.1636824700)
    BND = SBEND(name="BND", len=1.5048000000, angle=0.1847995700, e1=0.0923997850, e2=9.2399785000e-02, PolynomB=[0.0, -3.1537858000e-01, 0.0, 0.0])
    DC4A = DRIFT(name="DC4A", len=0.1592146700)
    DC4B = DRIFT(name="DC4B", len=0.0444180000)
    SD = KSEXT(name="SD", len=0.2500000000, k2=-38.8015300000)
    DC5A = DRIFT(name="DC5A", len=0.0905800000)
    DC5B = DRIFT(name="DC5B", len=0.3613900000)
    SF = KSEXT(name="SF", len=0.2100000000, k2=32.0477093000)
    DC6A = DRIFT(name="DC6A", len=0.1106460000)
    DC6B = DRIFT(name="DC6B", len=0.0631658500)
    QFC = KQUAD(name="QFC", len=0.5123803000, k1=1.7486408311)
    DC5C = DRIFT(name="DC5C", len=0.0958400000)
    DC5D = DRIFT(name="DC5D", len=0.3561300000)
    DC2C = DRIFT(name="DC2C", len=0.1021004500)
    DC2D = DRIFT(name="DC2D", len=0.1294752500)
    DI1 = DRIFT(name="DI1", len=0.9235741000)
    K1 = CORRECTOR(name="K1", len=1.2000000000, xkick=0.0000000000, ykick=0.0000000000)
    DI2 = DRIFT(name="DI2", len=0.6882939000)
    DI3 = DRIFT(name="DI3", len=0.6834939000)
    K2 = CORRECTOR(name="K2", len=0.6000000000, xkick=0.0000000000, ykick=0.0000000000)
    DI4 = DRIFT(name="DI4", len=0.1224401000)
    DI5 = DRIFT(name="DI5", len=1.2401300000)
    DI6 = DRIFT(name="DI6", len=0.1658040000)
    K3 = CORRECTOR(name="K3", len=1.2000000000, xkick=0.0000000000, ykick=0.0000000000)
    DB10B = DRIFT(name="DB10B", len=3.1458354000)
    DB10A = DRIFT(name="DB10A", len=0.1240412500)
    DB9B = DRIFT(name="DB9B", len=0.1233052500)
    DB9A = DRIFT(name="DB9A", len=0.1233052500)
    DB8B = DRIFT(name="DB8B", len=0.1385954500)
    DB8A = DRIFT(name="DB8A", len=0.3272502700)
    DB5D = DRIFT(name="DB5D", len=0.1139772700)
    DB5C = DRIFT(name="DB5C", len=0.1085632000)
    DB6D = DRIFT(name="DB6D", len=0.1270000000)
    DB6C = DRIFT(name="DB6C", len=0.9043700000)
    DB7B = DRIFT(name="DB7B", len=0.1106968000)
    DB7A = DRIFT(name="DB7A", len=0.0631130500)
    DB6B = DRIFT(name="DB6B", len=0.0939985200)
    DB6A = DRIFT(name="DB6A", len=0.9373700000)
    DB5B = DRIFT(name="DB5B", len=0.0518450000)
    DB5A = DRIFT(name="DB5A", len=0.1706954700)
    DB3B = DRIFT(name="DB3B", len=0.1308128500)
    DB3A = DRIFT(name="DB3A", len=0.1322268500)
    DB2B = DRIFT(name="DB2B", len=0.1172344500)
    DB2A = DRIFT(name="DB2A", len=0.1158052500)
    DB1B = DRIFT(name="DB1B", len=0.0562232500)
    DB1A = DRIFT(name="DB1A", len=3.7470820000)
    elements = [RF,DA1A,BPM,DA1B,QDX,DA2A,COR,DA2B,QFX,DA3A,BPM,DA3B,QDY,DM4,
    B34,DA5A,BPM,DA5B,SDM,DA6A,COR,DA6B,SFM,DA7A,BPM,DA7B,QFY,DM7,SFM,DA6C,COR,
    DA6D,SDM,DA5C,BPM,DA5D,B34,DA8A,BPM,DA8B,QDZ,DA9A,COR,DA9B,QFZ,DA10A,BPM,DA10B,
    DC1A,BPM,DC1B,QF,DC2A,COR,DC2B,QD,DC3A,BPM,DC3B,BND,DC4A,BPM,DC4B,SD,DC5A,COR,
    DC5B,SF,DC6A,BPM,DC6B,QFC,DC6B,DC6A,SF,DC5C,COR,DC5D,SD,DC4B,BPM,DC4A,BND,DC3B,
    DC3A,QD,DC2C,COR,DC2D,QF,DC1B,BPM,DI1,K1,DI2,BPM,DC1B,QF,DC2A,COR,DC2B,QD,DC3A,
    BPM,DC3B,BND,DC4A,BPM,DC4B,SD,DC5A,COR,DC5B,SF,DC6A,BPM,DC6B,QFC,DC6B,DC6A,SF,
    DC5C,COR,DC5D,SD,DC4B,BPM,DC4A,BND,DC3B,DC3A,QD,DC2C,COR,DC2D,QF,DC1B,BPM,DI3,
    K2,DI4,DI5,DI6,BPM,DC1B,QF,DC2A,COR,DC2B,QD,DC3A,BPM,DC3B,BND,DC4A,BPM,DC4B,SD,
    DC5A,COR,DC5B,SF,DC6A,BPM,DC6B,QFC,DC6B,DC6A,SF,DC5C,COR,DC5D,SD,DC4B,BPM,DC4A,
    BND,DC3B,DC3A,QD,DC2C,COR,DC2D,QF,DC1B,BPM,DI2,K3,DI1,BPM,DC1B,QF,DC2A,COR,DC2B,QD,
    DC3A,BPM,DC3B,BND,DC4A,BPM,DC4B,SD,DC5A,COR,DC5B,SF,DC6A,BPM,DC6B,QFC,DC6B,DC6A,SF,
    DC5C,COR,DC5D,SD,DC4B,BPM,DC4A,BND,DC3B,DC3A,QD,DC2C,COR,DC2D,QF,DC1B,BPM,DC1A,DC1A,
    BPM,DC1B,QF,DC2A,COR,DC2B,QD,DC3A,BPM,DC3B,BND,DC4A,BPM,DC4B,SD,DC5A,COR,DC5B,SF,DC6A,
    BPM,DC6B,QFC,DC6B,DC6A,SF,DC5C,COR,DC5D,SD,DC4B,BPM,DC4A,BND,DC3B,DC3A,QD,DC2C,COR,DC2D,
    QF,DC1B,BPM,DC1A,DC1A,BPM,DC1B,QF,DC2A,COR,DC2B,QD,DC3A,BPM,DC3B,BND,DC4A,BPM,DC4B,SD,
    DC5A,COR,DC5B,SF,DC6A,BPM,DC6B,QFC,DC6B,DC6A,SF,DC5C,COR,DC5D,SD,DC4B,BPM,DC4A,BND,DC3B,
    DC3A,QD,DC2C,COR,DC2D,QF,DC1B,BPM,DC1A,DC1A,BPM,DC1B,QF,DC2A,COR,DC2B,QD,DC3A,BPM,DC3B,BND,
    DC4A,BPM,DC4B,SD,DC5A,COR,DC5B,SF,DC6A,BPM,DC6B,QFC,DC6B,DC6A,SF,DC5C,COR,DC5D,SD,DC4B,BPM,
    DC4A,BND,DC3B,DC3A,QD,DC2C,COR,DC2D,QF,DC1B,BPM,DC1A,DB10B,BPM,DB10A,QFZ,DB9B,COR,DB9A,QDZ,
    DB8B,BPM,DB8A,B34,DB5D,BPM,DB5C,SDM,DB6D,COR,DB6C,SFM,DB7B,BPM,DB7A,QFY,DM7,SFM,DB6B,COR,DB6A,
    SDM,DB5B,BPM,DB5A,B34,DM4,QDY,DB3B,BPM,DB3A,QFX,DB2B,COR,DB2A,QDX,DB1B,BPM,DB1A,DA1A,BPM,DA1B,
    QDX,DA2A,COR,DA2B,QFX,DA3A,BPM,DA3B,QDY,DM4,B34,DA5A,BPM,DA5B,SDM,DA6A,COR,DA6B,SFM,DA7A,BPM,DA7B,
    QFY,DM7,SFM,DA6C,COR,DA6D,SDM,DA5C,BPM,DA5D,B34,DA8A,BPM,DA8B,QDZ,DA9A,COR,DA9B,QFZ,DA10A,BPM,DA10B,
    DC1A,BPM,DC1B,QF,DC2A,COR,DC2B,QD,DC3A,BPM,DC3B,BND,DC4A,BPM,DC4B,SD,DC5A,COR,DC5B,SF,DC6A,BPM,DC6B,
    QFC,DC6B,DC6A,SF,DC5C,COR,DC5D,SD,DC4B,BPM,DC4A,BND,DC3B,DC3A,QD,DC2C,COR,DC2D,QF,DC1B,BPM,DC1A,DC1A,
    BPM,DC1B,QF,DC2A,COR,DC2B,QD,DC3A,BPM,DC3B,BND,DC4A,BPM,DC4B,SD,DC5A,COR,DC5B,SF,DC6A,BPM,DC6B,QFC,
    DC6B,DC6A,SF,DC5C,COR,DC5D,SD,DC4B,BPM,DC4A,BND,DC3B,DC3A,QD,DC2C,COR,DC2D,QF,DC1B,BPM,DC1A,DC1A,
    BPM,DC1B,QF,DC2A,COR,DC2B,QD,DC3A,BPM,DC3B,BND,DC4A,BPM,DC4B,SD,DC5A,COR,DC5B,SF,DC6A,BPM,DC6B,QFC,
    DC6B,DC6A,SF,DC5C,COR,DC5D,SD,DC4B,BPM,DC4A,BND,DC3B,DC3A,QD,DC2C,COR,DC2D,QF,DC1B,BPM,DC1A,DC1A,
    BPM,DC1B,QF,DC2A,COR,DC2B,QD,DC3A,BPM,DC3B,BND,DC4A,BPM,DC4B,SD,DC5A,COR,DC5B,SF,DC6A,BPM,DC6B,QFC,
    DC6B,DC6A,SF,DC5C,COR,DC5D,SD,DC4B,BPM,DC4A,BND,DC3B,DC3A,QD,DC2C,COR,DC2D,QF,DC1B,BPM,DC1A,DC1A,
    BPM,DC1B,QF,DC2A,COR,DC2B,QD,DC3A,BPM,DC3B,BND,DC4A,BPM,DC4B,SD,DC5A,COR,DC5B,SF,DC6A,BPM,DC6B,
    QFC,DC6B,DC6A,SF,DC5C,COR,DC5D,SD,DC4B,BPM,DC4A,BND,DC3B,DC3A,QD,DC2C,COR,DC2D,QF,DC1B,BPM,DC1A,
    DC1A,BPM,DC1B,QF,DC2A,COR,DC2B,QD,DC3A,BPM,DC3B,BND,DC4A,BPM,DC4B,SD,DC5A,COR,DC5B,SF,DC6A,BPM,
    DC6B,QFC,DC6B,DC6A,SF,DC5C,COR,DC5D,SD,DC4B,BPM,DC4A,BND,DC3B,DC3A,QD,DC2C,COR,DC2D,QF,DC1B,BPM,
    DC1A,DC1A,BPM,DC1B,QF,DC2A,COR,DC2B,QD,DC3A,BPM,DC3B,BND,DC4A,BPM,DC4B,SD,DC5A,COR,DC5B,SF,DC6A,
    BPM,DC6B,QFC,DC6B,DC6A,SF,DC5C,COR,DC5D,SD,DC4B,BPM,DC4A,BND,DC3B,DC3A,QD,DC2C,COR,DC2D,QF,DC1B,
    BPM,DC1A,DB10B,BPM,DB10A,QFZ,DB9B,COR,DB9A,QDZ,DB8B,BPM,DB8A,B34,DB5D,BPM,DB5C,SDM,DB6D,COR,DB6C,
    SFM,DB7B,BPM,DB7A,QFY,DM7,SFM,DB6B,COR,DB6A,SDM,DB5B,BPM,DB5A,B34,DM4,QDY,DB3B,BPM,DB3A,QFX,DB2B,COR,DB2A,QDX,DB1B,BPM,DB1A]
end
# serialize("src/demo/SPEAR3/spear3.jls", elements)