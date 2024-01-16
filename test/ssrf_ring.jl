include("../lattice/canonical_elements.jl")
function ssrf(K1, L1)
    E0 = 3500.0 # MeV
    gamma0 = E0/0.51099906
    p0 = sqrt(gamma0^2 - 1.0)
    
    D1 = EDRIFT(name="D1", len=L1)
    D2 = EDRIFT(name="D2", len=0.12)
    D3 = EDRIFT(name="D3", len=0.475)
    D4 = EDRIFT(name="D4", len=0.12)
    D5 = EDRIFT(name="D5", len=0.28)
    D6 = EDRIFT(name="D6", len=0.59) 
    D7 = EDRIFT(name="D7", len=0.36)
    D8 = EDRIFT(name="D8", len=0.45)
    D9 = EDRIFT(name="D9", len=0.18)
    D33 = EDRIFT(name="D33", len=0.58)
    D51 = EDRIFT(name="D51", len=0.12)
    D52 = EDRIFT(name="D52", len=0.60)
    DL0 = EDRIFT(name="DL0", len=0.60)
    DM0 = EDRIFT(name="DM0", len=0.325)
    DL = EDRIFT(name="DL", len=6.0)
    DM = EDRIFT(name="DM", len=3.25)
    
    QL1 = KQUAD(name="QL1", len=0.32, k1=K1, nSlices=10, synch_rad=1)
    QL2H = KQUAD(name="QL2H", len=0.29, k1=1.358860, nSlices=10, synch_rad=1)
    QL3 = KQUAD(name="QL3", len=0.32, k1=-1.192160, nSlices=10, synch_rad=1)
    QL4 = KQUAD(name="QL4", len=0.26, k1=-1.077410, nSlices=10, synch_rad=1)
    QL5 = KQUAD(name="QL5", len=0.32, k1=1.392450, nSlices=10, synch_rad=1)
    QM1 = KQUAD(name="QM1", len=0.32, k1=-1.562500, nSlices=10, synch_rad=1)
    QM3 = KQUAD(name="QM3", len=0.32, k1=-1.014220, nSlices=10, synch_rad=1)
    QM4 = KQUAD(name="QM4", len=0.26, k1=-1.366690, nSlices=10, synch_rad=1)
    QM5 = KQUAD(name="QM5", len=0.32, k1=1.455000, nSlices=10, synch_rad=1)
    Q2H = KQUAD(name="Q2H", len=0.29, k1=1.532730, nSlices=10, synch_rad=1)
    
    # Sextupole strength is 2 times of the AT lattice
    S1 = KSEXT(name="S1", len=0.2, k2=1.555155/0.2*2, nSlices=10, synch_rad=1)
    S2 = KSEXT(name="S2", len=0.24, k2=-3.001088/0.24*2, nSlices=10, synch_rad=1)
    S3 = KSEXT(name="S3", len=0.2, k2=2.542476/0.2*2, nSlices=10, synch_rad=1)
    S4 = KSEXT(name="S4", len=0.24, k2=-2.691814/0.24*2, nSlices=10, synch_rad=1)
    S5 = KSEXT(name="S5", len=0.2, k2=3.540568/0.2*2, nSlices=10, synch_rad=1)
    S6 = KSEXT(name="S6", len=0.24, k2=-4.578491/0.24*2, nSlices=10, synch_rad=1)
    SD = KSEXT(name="SD", len=0.2, k2=-2.424032/0.2*2, nSlices=10, synch_rad=1)
    SF = KSEXT(name="SF", len=0.24, k2=3.436611/0.24*2, nSlices=10, synch_rad=1)
    
    # For a symmetric rectangular magnet, E1=E2=ANGLE/2
    BendingAngle = pi/20
    BD1 = CSBEND(name="BD1", len=0.72, angle=BendingAngle/2, e1=BendingAngle/2, e2=0.0, nSlices=10, synch_rad=1)
    BD2 = CSBEND(name="BD2", len=0.72, angle=BendingAngle/2, e1=0.0, e2=BendingAngle/2, nSlices=10, synch_rad=1)
    
    L0 = 432.0+40*1.44*(BendingAngle/2/sin(BendingAngle/2.0)-1)
    C0 =   299792458 	# speed of light [m/s]
    HarmNumber = 720
    # original length 0.0, but cannot be 0.0 in TPSA
    CAV = RFCA(name="CAV1", len=0.1, volt=4.0e6, freq=HarmNumber*C0/L0, phase=0.0, nSlices=10)
    
    # M1 = [DL,  QL1, D1,  S1,   D2,   QL2H, QL2H,  D33, 
    # QL3,  D51,  S2,  D52,  BD1, BD2,  D6,   QL4,
    # D7,   SD,   D8, QL5,  D9,  SF,
    # D9,   QM5,  D8,  SD,   D7,  QM4,  D6,  BD1, BD2,
    # D5,   QM3,  D4,  S4,   D3,  Q2H, Q2H,  D2,  S3,  D1, QM1, DM]
    
    # M2 =[  DM,   QM1, D1,   S3,  D2,   Q2H, Q2H, D3,  
    #        S4,   D4 ,  QM3, D5,  BD1, BD2,  D6,   QM4,  
    #        D7,   SD ,  D8,   QM5, D9,  SF, 
    #        D9,   QM5,  D8,  SD,  D7,  QM4,  D6,  BD1, BD2, 
    #        D5,   QM3,  D4,  S6,   D3,  Q2H, Q2H,  D2,  S5,  D1, QM1, DM]    
         
    # M3 = [ DM,   QM1,  D1,   S5,  D2,   Q2H, Q2H, D3,  
    #        S6,   D4,   QM3,  D5,  BD1, BD2,  D6,   QM4,  
    #        D7,   SD,   D8,    QM5, D9,  SF, 
    #        D9,   QM5,  D8,   SD,  D7,  QM4,  D6,  BD1, BD2, 
    #        D5,   QM3,  D4,   S6,  D3,  Q2H, Q2H,  D2,  S5,  D1,  QM1, DM ] 
           
    # M4 =[  DM,   QM1,   D1,   S5,  D2,  Q2H, Q2H,  D3,  
    #        S6,   D4 ,  QM3,  D5,  BD1, BD2,  D6,   QM4,  
    #        D7,   SD ,  D8,   QM5, D9,  SF, 
    #        D9,   QM5,  D8,   SD,   D7,  QM4,  D6,  BD1, BD2, 
    #        D5,   QM3,  D4,   S4,  D3,  Q2H, Q2H,  D2,  S3,  D1,  QM1, DM ]
    
    # M5 =[  DM,   QM1,  D1,  S3,  D2,  Q2H, Q2H,  D3,  
    #        S4,   D4,   QM3,  D5,  BD1, BD2,  D6,   QM4,  
    #        D7,   SD,   D8,   QM5, D9,  SF, 
    #        D9,   QL5,  D8,   SD,  D7, QL4,  D6,  BD1, BD2,
    #        D52,  S2,   D51,  QL3,  D33,  QL2H, QL2H,  D2, S1,  D1,  QL1,  DL ]   
    # cell = [M1..., M2..., M3..., M4..., M5...]
    # ELIST = [CAV, cell..., cell..., cell..., cell...]

    Mtest = [DL,  QL1, D1]#,  S1,   D2,   QL2H, QL2H,  D33, 
    # QL3,  D51,  S2,  D52,  BD1, BD2,  D6,   QL4,
    # D7,   SD,   D8, QL5,  D9,  SF,
    # D9,   QM5,  D8,  SD,   D7,  QM4,  D6,  BD1, BD2,
    # D5,   QM3,  D4,  S4,   D3,  Q2H, Q2H,  D2,  S3,  D1, QM1, DM]

    return Mtest#ELIST
end