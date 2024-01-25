# include("../src/JuTrack.jl")
using. JuTrack
# K1, L1 = [-1.063770, 0.34]
function ssrf(K1)
    E0 = 3500.0 # MeV
    gamma0 = E0/0.51099906
    p0 = sqrt(gamma0^2 - 1.0)
    
    D1 = EDRIFT(name="D1", len=0.34)
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
    
    QL1 = KQUAD(name="QL1", len=0.32, k1=K1, nSlices=10, synch_rad=0)
    QL2H = KQUAD(name="QL2H", len=0.29, k1=1.358860, nSlices=10, synch_rad=0)
    QL3 = KQUAD(name="QL3", len=0.32, k1=-1.192160, nSlices=10, synch_rad=0)
    QL4 = KQUAD(name="QL4", len=0.26, k1=-1.077410, nSlices=10, synch_rad=0)
    QL5 = KQUAD(name="QL5", len=0.32, k1=1.392450, nSlices=10, synch_rad=0)
    QM1 = KQUAD(name="QM1", len=0.32, k1=-1.562500, nSlices=10, synch_rad=0)
    QM3 = KQUAD(name="QM3", len=0.32, k1=-1.014220, nSlices=10, synch_rad=0)
    QM4 = KQUAD(name="QM4", len=0.26, k1=-1.366690, nSlices=10, synch_rad=0)
    QM5 = KQUAD(name="QM5", len=0.32, k1=1.455000, nSlices=10, synch_rad=0)
    Q2H = KQUAD(name="Q2H", len=0.29, k1=1.532730, nSlices=10, synch_rad=0)
    
    # Sextupole strength is 2 times of the AT lattice
    S1 = KSEXT(name="S1", len=0.2, k2=1.555155/0.2*2, nSlices=10, synch_rad=0)
    S2 = KSEXT(name="S2", len=0.24, k2=-3.001088/0.24*2, nSlices=10, synch_rad=0)
    S3 = KSEXT(name="S3", len=0.2, k2=2.542476/0.2*2, nSlices=10, synch_rad=0)
    S4 = KSEXT(name="S4", len=0.24, k2=-2.691814/0.24*2, nSlices=10, synch_rad=0)
    S5 = KSEXT(name="S5", len=0.2, k2=3.540568/0.2*2, nSlices=10, synch_rad=0)
    S6 = KSEXT(name="S6", len=0.24, k2=-4.578491/0.24*2, nSlices=10, synch_rad=0)
    SD = KSEXT(name="SD", len=0.2, k2=-2.424032/0.2*2, nSlices=10, synch_rad=0)
    SF = KSEXT(name="SF", len=0.24, k2=3.436611/0.24*2, nSlices=10, synch_rad=0)
    
    # For a symmetric rectangular magnet, E1=E2=ANGLE/2
    BendingAngle = pi/20
    BD1 = CSBEND(name="BD1", len=0.72, angle=BendingAngle/2, e1=BendingAngle/2, e2=0.0, nSlices=10, synch_rad=0)
    BD2 = CSBEND(name="BD2", len=0.72, angle=BendingAngle/2, e1=0.0, e2=BendingAngle/2, nSlices=10, synch_rad=0)
    
    L0 = 432.0+40*1.44*(BendingAngle/2/sin(BendingAngle/2.0)-1)
    C0 =   299792458 	# speed of light [m/s]
    HarmNumber = 720
    # original length 0.0, but cannot be 0.0 in TPSA
    CAV = RFCA(name="CAV1", len=0.1, volt=4.0e6, freq=HarmNumber*C0/L0, phase=0.0, nSlices=10)
    
    M1 = (DL,  QL1, D1,  S1,   D2,   QL2H, QL2H,  D33, 
    QL3,  D51,  S2,  D52,  BD1, BD2,  D6,   QL4,
    D7,   SD,   D8, QL5,  D9,  SF,
    D9,   QM5,  D8,  SD,   D7,  QM4,  D6,  BD1, BD2,
    D5,   QM3,  D4,  S4,   D3,  Q2H, Q2H,  D2,  S3,  D1, QM1, DM)
    
    M2 =(DM,   QM1, D1,   S3,  D2,   Q2H, Q2H, D3,  
           S4,   D4 ,  QM3, D5,  BD1, BD2,  D6,   QM4,  
           D7,   SD ,  D8,   QM5, D9,  SF, 
           D9,   QM5,  D8,  SD,  D7,  QM4,  D6,  BD1, BD2, 
           D5,   QM3,  D4,  S6,   D3,  Q2H, Q2H,  D2,  S5,  D1, QM1, DM)   
         
    M3 = (DM,   QM1,  D1,   S5,  D2,   Q2H, Q2H, D3,  
           S6,   D4,   QM3,  D5,  BD1, BD2,  D6,   QM4,  
           D7,   SD,   D8,    QM5, D9,  SF, 
           D9,   QM5,  D8,   SD,  D7,  QM4,  D6,  BD1, BD2, 
           D5,   QM3,  D4,   S6,  D3,  Q2H, Q2H,  D2,  S5,  D1,  QM1, DM) 
           
    M4 = (DM,   QM1,   D1,   S5,  D2,  Q2H, Q2H,  D3,  
           S6,   D4 ,  QM3,  D5,  BD1, BD2,  D6,   QM4,  
           D7,   SD ,  D8,   QM5, D9,  SF, 
           D9,   QM5,  D8,   SD,   D7,  QM4,  D6,  BD1, BD2, 
           D5,   QM3,  D4,   S4,  D3,  Q2H, Q2H,  D2,  S3,  D1,  QM1, DM)
    
    M5 = (DM,   QM1,  D1,  S3,  D2,  Q2H, Q2H,  D3,  
           S4,   D4,   QM3,  D5,  BD1, BD2,  D6,   QM4,  
           D7,   SD,   D8,   QM5, D9,  SF, 
           D9,   QL5,  D8,   SD,  D7, QL4,  D6,  BD1, BD2,
           D52,  S2,   D51,  QL3,  D33,  QL2H, QL2H,  D2, S1,  D1,  QL1,  DL)
# #     ELIST = [CAV, cell..., cell..., cell..., cell...]
       # CELL = (M1..., M2..., M3..., M4..., M5...)
       # CELL = (M1[1], M1[2], M1[3], M1[4], M1[5], M1[6], M1[7], M1[8], M1[9], M1[10], M1[11], M1[12], M1[13], M1[14], M1[15], M1[16], M1[17], M1[18], M1[19], M1[20], M1[21], M1[22], M1[23], M1[24], M1[25], M1[26], M1[27], M1[28], M1[29], M1[30], M1[31], M1[32], M1[33], M1[34], M1[35], M1[36], M1[37], M1[38], M1[39], M1[40], M1[41], M1[42], M1[43],
       #         M2[1], M2[2], M2[3], M2[4], M2[5], M2[6], M2[7], M2[8], M2[9], M2[10], M2[11], M2[12], M2[13], M2[14], M2[15], M2[16], M2[17], M2[18], M2[19], M2[20], M2[21], M2[22], M2[23], M2[24], M2[25], M2[26], M2[27], M2[28], M2[29], M2[30], M2[31], M2[32], M2[33], M2[34], M2[35], M2[36], M2[37], M2[38], M2[39], M2[40], M2[41], M2[42], M2[43],
       #         M3[1], M3[2], M3[3], M3[4], M3[5], M3[6], M3[7], M3[8], M3[9], M3[10], M3[11], M3[12], M3[13], M3[14], M3[15], M3[16], M3[17], M3[18], M3[19], M3[20], M3[21], M3[22], M3[23], M3[24], M3[25], M3[26], M3[27], M3[28], M3[29], M3[30], M3[31], M3[32], M3[33], M3[34], M3[35], M3[36], M3[37], M3[38], M3[39], M3[40], M3[41], M3[42], M3[43],
       #         M4[1], M4[2], M4[3], M4[4], M4[5], M4[6], M4[7], M4[8], M4[9], M4[10], M4[11], M4[12], M4[13], M4[14], M4[15], M4[16], M4[17], M4[18], M4[19], M4[20], M4[21], M4[22], M4[23], M4[24], M4[25], M4[26], M4[27], M4[28], M4[29], M4[30], M4[31], M4[32], M4[33], M4[34], M4[35], M4[36], M4[37], M4[38], M4[39], M4[40], M4[41], M4[42], M4[43],
       #         M5[1], M5[2], M5[3], M5[4], M5[5], M5[6], M5[7], M5[8], M5[9], M5[10], M5[11], M5[12], M5[13], M5[14], M5[15], M5[16], M5[17], M5[18], M5[19], M5[20], M5[21], M5[22], M5[23], M5[24], M5[25], M5[26], M5[27], M5[28], M5[29], M5[30], M5[31], M5[32], M5[33], M5[34], M5[35], M5[36], M5[37], M5[38], M5[39], M5[40], M5[41], M5[42], M5[43])
       CELL = Vector{AbstractElement}(undef, 215)
       for i in eachindex(M1)
              CELL[i] = M1[i]
              CELL[i+43] = M2[i]
              CELL[i+86] = M3[i]
              CELL[i+129] = M4[i]
              CELL[i+172] = M5[i]
       end
       ELIST = Vector{AbstractElement}(undef, 4*length(CELL))
       for i in eachindex(CELL)
              ELIST[i] = CELL[i]
              ELIST[i+length(CELL)] = CELL[i]
              ELIST[i+2*length(CELL)] = CELL[i]
              ELIST[i+3*length(CELL)] = CELL[i]
       end

       # ELIST = (CELL..., CELL..., CELL..., CELL...)


    return ELIST
end