# include("../src/JuTrack.jl")
# using. JuTrack
# K1, L1 = [-1.063770, 0.34]
function ssrf(K1,RAD)
       if RAD == 0
          VOLT = 0.0
       else
          VOLT = 4.0e6
       end
       E0 = 3500.0 # MeV
       gamma0 = E0/0.51099906
       p0 = sqrt(gamma0^2 - 1.0)
   
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
       
       QL1 = KQUAD(name="QL1", len=0.32, k1=K1, rad=RAD)
       QL2H = KQUAD(name="QL2H", len=0.29, k1=1.358860, rad=RAD)
       QL3 = KQUAD(name="QL3", len=0.32, k1=-1.192160, rad=RAD)
       QL4 = KQUAD(name="QL4", len=0.26, k1=-1.077410, rad=RAD)
       QL5 = KQUAD(name="QL5", len=0.32, k1=1.392450, rad=RAD)
       QM1 = KQUAD(name="QM1", len=0.32, k1=-1.562500, rad=RAD)
       QM3 = KQUAD(name="QM3", len=0.32, k1=-1.014220, rad=RAD)
       QM4 = KQUAD(name="QM4", len=0.26, k1=-1.366690, rad=RAD)
       QM5 = KQUAD(name="QM5", len=0.32, k1=1.455000, rad=RAD)
       Q2H = KQUAD(name="Q2H", len=0.29, k1=1.532730, rad=RAD)
       
       # Sextupole strength is 2 times of the AT lattice
       S1 = KSEXT(name="S1", len=0.2, k2=1.555155/0.2*2, rad=RAD)
       S2 = KSEXT(name="S2", len=0.24, k2=-3.001088/0.24*2, rad=RAD)
       S3 = KSEXT(name="S3", len=0.2, k2=2.542476/0.2*2, rad=RAD)
       S4 = KSEXT(name="S4", len=0.24, k2=-2.691814/0.24*2, rad=RAD)
       S5 = KSEXT(name="S5", len=0.2, k2=3.540568/0.2*2, rad=RAD)
       S6 = KSEXT(name="S6", len=0.24, k2=-4.578491/0.24*2, rad=RAD)
       SD = KSEXT(name="SD", len=0.2, k2=-2.424032/0.2*2, rad=RAD)
       SF = KSEXT(name="SF", len=0.24, k2=3.436611/0.24*2, rad=RAD)
       
       # For a symmetric rectangular magnet, E1=E2=ANGLE/2
       BendingAngle = pi/20
       BD1 = SBEND(name="BD1", len=0.72, angle=BendingAngle/2, e1=BendingAngle/2, e2=0.0, rad=RAD)
       BD2 = SBEND(name="BD2", len=0.72, angle=BendingAngle/2, e1=0.0, e2=BendingAngle/2, rad=RAD)
       
       L0 = 432.0+40*1.44*(BendingAngle/2/sin(BendingAngle/2.0)-1)
       C0 =   299792458 	# speed of light [m/s]
       HarmNumber = 720
       # original length 0.0, but cannot be 0.0 in TPSA
       CAV = RFCA(name="CAV1", volt=VOLT, freq=HarmNumber*C0/L0, energy=E0*1e6)
       BPM = MARKER(name="BPM")
       HVC = CORRECTOR(name="HVC")
   
       M1 = [DL,  QL1, BPM, D1,  S1,   D2,   QL2H, QL2H,  D33, 
       QL3,  D51,  S2,  D52,  BD1, BD2,  D6,   QL4,
       D7,   SD,   D8, HVC, QL5,  D9,  SF,
       D9,   QM5,  D8,  SD,  HVC,  D7,  QM4,  D6,  BD1, BD2,
       D5,   QM3,  D4,  S4,   D3,  Q2H, Q2H,  D2,  S3, BPM, D1, QM1, DM]
       
       M2 =[DM,   QM1, BPM, D1,   S3,  D2,   Q2H, Q2H, D3,  
              S4,   D4 ,  QM3, HVC, D5,  BD1, BD2,  D6,   QM4,  
              D7,   SD ,  D8,   QM5, D9,  SF, 
              D9,   QM5,  D8,  SD, HVC, D7,  QM4,  D6,  BD1, BD2, 
              D5,   QM3,  D4,  S6,   D3,  Q2H, Q2H,  D2,  S5,  D1,BPM, QM1, DM]
            
       M3 = [DM,   QM1, BPM, D1,   S5,  D2,   Q2H, Q2H, D3,  
              S6,   D4,   QM3,  D5,  BD1, BD2,  D6,   QM4,  
              D7,   SD,   D8,  HVC,  QM5, D9,  SF, 
              D9,   QM5,  D8,   SD, HVC,  D7,  QM4,  D6,  BD1, BD2, 
              D5,   QM3,  D4,   S6,  D3,  Q2H, Q2H,  D2,  S5, BPM, D1,  QM1, DM]
              
       M4 = [DM,   QM1, BPM,  D1,   S5,  D2,  Q2H, Q2H,  D3,  
              S6,   D4 ,  QM3, HVC, D5,  BD1, BD2,  D6,   QM4,  
              D7,   SD ,  D8,   QM5, D9,  SF, 
              D9,   QM5,  D8,   SD, HVC,  D7,  QM4,  D6,  BD1, BD2, 
              D5,   QM3,  D4,   S4,  D3,  Q2H, Q2H,  D2,  S3, BPM, D1,  QM1, DM]
       
       M5 =[DM,   QM1,  D1, BPM, S3,  D2,  Q2H, Q2H,  D3,  
              S4,   D4,   QM3,  D5, HVC, BD1, BD2,  D6,   QM4,  
              D7,   SD,   D8,   QM5, D9,  SF, 
              D9,   QL5,  D8,   SD,  D7, HVC, QL4,  D6,  BD1, BD2,
              D52,  S2,   D51,  QL3,  D33,  QL2H, QL2H,  D2, BPM, S1,  D1,  QL1,  DL]
   
   
          CELL = [M1..., M2..., M3..., M4..., M5...]
          ELIST = [CAV, CELL..., CELL..., CELL..., CELL...]
   
   
       return ELIST
   end