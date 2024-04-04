# ESR lattice of EIC
# translated from MADX, esr-main-18GeV-2IP
# translated by Jinyu Wan, 2024-03-11
using Serialization
include("../src/JuTrack.jl")
using. JuTrack

U0 = 3.75583138918653e7 # Energy loss per turn for reference particle
clight = 299792458.0 # m/s
emass = 0.51099895000e-3 # GeV
A_ELECTRON = 0.0011596521869 
E = 17.846262619763e9
AGAMMA_TOP = 40.5
AGAMMA_MED = 22.5
AGAMMA_BOT = 11.5
GAMMA_TOP = AGAMMA_TOP/A_ELECTRON
GAMMA_MED = AGAMMA_MED/A_ELECTRON
GAMMA_BOT = AGAMMA_BOT/A_ELECTRON
theta2 = (pi/2.0) / AGAMMA_TOP 
theta = (pi/2.0) / AGAMMA_BOT 
theta1 = theta - theta2 
phi1_6 = 0.0
phi2_6 = 1.570796327
phi1_8 = phi1_6
phi2_8 = phi2_6
TH2_DCCN = 2.0

# RF cavity parameters.
HRMN_RF     =   7560.0
VOLT_RF     =      3.78e6
LAG_RF      =      0.5847649913891533 # 0.0
LRF         =      4.01667

# Crab cavity parameters.
crab14 = 2.8713953711409927e6
crab23 = 2.8713762590139e6
# ------- The geometry ------------------------------------
# Number of bends per ARC.
NBENDS_ARC  = 32.0
#
# Number of bends in the 'dispersion suppressor'/'rotator' modules.
NDISP = 6.0
#
# ARC bend magnet parameters.
ROT_ANGLE = theta/12 # Angle of one dipole in non-colliding straight sections
RHO_AVE_CENTER = 380.494
ARC_BEND_ANGLE = (pi/3.0 - 2.0*theta)/NBENDS_ARC

L01_STR = 2.726
L12_STR = 0.2584
L23_STR = 0.8914

TH2_XX = asin(L23_STR*sin(0.5*ARC_BEND_ANGLE)/(2.0*L01_STR + L23_STR))
L01_XX = L01_STR*(0.5*ARC_BEND_ANGLE - TH2_XX)/(sin(0.5*ARC_BEND_ANGLE) - sin(TH2_XX))
L23_XX = L23_STR*TH2_XX/sin(TH2_XX)
L12_XX = L12_STR/cos(TH2_XX)
LBEND_ARC = 2.0*L01_XX + 2.0*L12_XX + L23_XX

# IR.
g_d2er = 0.198*clight/sqrt(18*18-emass*emass)*1e-9
lqir06 = 0.6 
lqir12 = 1.2 

lq0ef_6 = 1.2 
lq1ef_6 = 1.61 
 
lq1er_6 = 1.8 
lq2er_6 = 1.4 

mid_q1er = 6.2
mid_q2er = 8.3
mid_d2er = 12.25
mid_q3er = 38.0

lrd2er   = 5.5
loder_6  = 0.4 
lr_rot_bend = L01_STR
lecrab = 4.0

beg_d2er = mid_d2er - lrd2er/2
end_d2er = mid_d2er + lrd2er/2

beg_q1er = mid_q1er - lq1er_6/2
end_q1er = mid_q1er + lq1er_6/2

beg_q2er = mid_q2er - lq2er_6/2
end_q2er = mid_q2er + lq2er_6/2

beg_q3er = mid_q3er - lqir06/2

#
# Generated from surveyFit-26th.py ...
#   8 mrad crossing ANGLE at IP6 and 11 mrad at IP8.
IP6_ANGLE =  0.016 
IP8_ANGLE =  0.022 
LDF_1 =   2.15422
LDF_2 =   2.08140
LDF_3 =   1.14255
LDF_4 =   1.11651
LDF_9 =   3.96268
LDF_10 =   3.40659
LDF_11 =   3.16505
LDF_12 =  5.23352
LDF20_1 =   6.00508
LDF20_2 =   6.15960
LDF20_3 =  6.77124
LDF20_4 =  6.95336
LDF20_9 =  2.40407
LDF20_10 =   2.87168
LDF20_11 =   2.81883
LDF20_12 =  3.25906
LDF4B_5 =   2.5e-01 
LDF4B_6 =   3.64227e-01
LDF4B_7 =  2.5840e-01 
LDF4B_8 =   2.5e-01
LHORIZ_1 =  1.9718934272101478 
LHORIZ_3 =  1.9701806856414805 
LHORIZ_5 =  -1.1078683276701895 
LHORIZ_7 =  1.4782231393720622 
LHORIZ_9 =  -1.8526683215727076 
LHORIZ_11 =  -0.805701090267803 
LLD09_0 =   1.49354e+01
LLD10_0 =   1.24975e+01
LLD11_0 =   2.35713
LLD12_0 =   2.45587
LLFT0_2 =    4.00002e-01
LLFCN_2 =  6.47220
LD4_3 = 9.90284
LD4_4 =  1.01158e+01
RHO_DISP =  254.98453445947476 
DB1EF_ANG_8 =  -0.0014944764371142855 
DB3ER_ANG_8 =  4.89854e-03 
L3EF_6 = 1.98657 
L4ER_6 =  3.5307579458463034 
L3EF_8 = 3.45318 
L4ER_8 =  9.04032 
LO5A_6f =  1.100 
LO5A_6r =  1.100 
LO5A_8f =  2.58402e-01
LO5A_8r =   5.65960e-01
LDBC_9 =  1.5077463196796046 
LDBC_10 =  1.2017986077902518 
#
# IR dogleg dipole angles.
db1ef_ang_6 = -0.0015
db2ef_ang_6 = -1.23486e-02
db3ef_ang_6 = 0.013 
db4ef_ang_6 = 0.0015
db5ef_ang_6 = (theta2 - db1ef_ang_6 - db2ef_ang_6*2 
- db3ef_ang_6 - db4ef_ang_6*2.0)/4.0
DB2ER_ANG_6 = -2.0*asin(0.5*g_d2er*lrd2er)
DB4ER_ANG_6 =  1.19246e-02 
DB3ER_ANG_6 = (38.79e-3 - DB2ER_ANG_6 - 5*DB4ER_ANG_6)
DB3EF_ANG_8 =  1.65327e-02
DB2EF_ANG_8 = (theta2 - (4.0*DB3EF_ANG_8 + 2.0*DB1EF_ANG_8))/2
DB2ER_ANG_8 = -2.0*asin(0.5*g_d2er*lrd2er)
DB4ER_ANG_8 =  0.01338427324 
DB5ER_ANG_8 = (theta2 - (DB4ER_ANG_8 + DB3ER_ANG_8 + DB2ER_ANG_8))/3.0
#
# Average Arc radius of curvature.
RHO_AVE_1 = RHO_AVE_CENTER + LHORIZ_1
RHO_AVE_3 = RHO_AVE_CENTER + LHORIZ_3
RHO_AVE_5 = RHO_AVE_CENTER + LHORIZ_5
RHO_AVE_7 = RHO_AVE_CENTER + LHORIZ_7
RHO_AVE_9 = RHO_AVE_CENTER + LHORIZ_9
RHO_AVE_11 = RHO_AVE_CENTER + LHORIZ_11
RHO_AVE_ALL = RHO_AVE_CENTER + (LHORIZ_1 + LHORIZ_3 + LHORIZ_5 + LHORIZ_7 +
   LHORIZ_9 + LHORIZ_11)/6.0
#
# Arc FODO lengths
# Sextupole length.
LSX       =     0.24
LSXL      =     0.57
# Corrector length.
LCH       =     0.2
LCV       =     0.2
# Space between magnetic devices.
LDX17     =     0.200
LDF29     =     0.535
# Arc-cell drift lengths
LHQC = 0.064
LDSX = 0.500
LQSX = 0.156
LSSX = 0.104
LQSXL = LQSX + LSX + LSSX/2 - LSXL/2 
LDSXL = LQSX + 2*LSX + LSSX + LDSX - LQSXL - LSXL
LSS1 = 0.1721
# Quad lengths
LQ50 = 0.5
LQ60 = 0.6
LQ80 = 0.8

# Arc-cell lengths
RHO = LBEND_ARC/ARC_BEND_ANGLE
LSTR_ARC = 2*LQ50 + 4*LSX + 2*LSSX + LCH + LCV + 2*LHQC + 2*LDSX + 2*LQSX
LDC_AVE = LBEND_ARC*(RHO_AVE_ALL/RHO - 1.0) - LSTR_ARC/2
LDC_1   = LBEND_ARC*(RHO_AVE_1/RHO - 1.0) - LSTR_ARC/2
LDC_3   = LBEND_ARC*(RHO_AVE_3/RHO - 1.0) - LSTR_ARC/2
LDC_5   = LBEND_ARC*(RHO_AVE_5/RHO - 1.0) - LSTR_ARC/2
LDC_7   = LBEND_ARC*(RHO_AVE_7/RHO - 1.0) - LSTR_ARC/2
LDC_9   = LBEND_ARC*(RHO_AVE_9/RHO - 1.0) - LSTR_ARC/2
LDC_11  = LBEND_ARC*(RHO_AVE_11/RHO - 1.0) - LSTR_ARC/2

#
# Lengths used in the straight sections.
LDBC_1  = LBEND_ARC + LDC_1 - 2*L01_STR - L12_STR 
LDBC_3  = LBEND_ARC + LDC_3 - 2*L01_STR - L12_STR 
LDBC_5  = LBEND_ARC + LDC_5 - 2*L01_STR - L12_STR 
LDBC_11  = LBEND_ARC + LDC_11 - 2*L01_STR - L12_STR 
LLMID_12 =  12.83670696688165 
LDBQ = 2*LSX + LSSX + LDSX + LQSX
#
# IP4 and IP12 dogleg.
ADB12_4  = 0.0327249/3
#
# ======= The dipoles. ====================================
ARC_ANG_M6 = ARC_BEND_ANGLE - IP6_ANGLE/8.0
ARC_ANG_P6 = ARC_BEND_ANGLE + IP6_ANGLE/8.0
ARC_ANG_P8 = ARC_BEND_ANGLE + IP8_ANGLE/8.0
ARC_ANG_M8 = ARC_BEND_ANGLE - IP8_ANGLE/8.0

# Spin rotator solenoid lengths, short and long.
LSOL5  =  2.5
LQSS1  =  0.6480402
LQSS2  =  0.9550568
LQSS3  =  1.634532
LQSS4  =  1.020723
LQSS5  =  0.6861532

#   long rotator
LSOL20 = 6.2
LQLS1 = 9.819319e-01
LQLS2 = 1.8
LQLS3 = 1.8
LQLS4 = 5.187944e-01


LBEND = 0.0

#
# ARC quadrupole strengths and support parameterrs
KD_1               =      -0.3112098116 
KF_1               =       0.3113930481 
KD_3               =      -0.3112112319 
KF_3               =       0.3113944689 
KD_5               =      -0.3137884525 
KF_5               =       0.3139724652 
KD_5a              =      -0.2006935196 
KF_5a              =       0.5299652224 
KD_5b             =       KD_5 
KF_5b             =       KF_5 
KD_5c             =       KD_5 
KF_5c             =       KF_5 
KD_7               =      -0.3116231882 
KF_7               =        0.311806547 
KD_6              =       KD_7 
KF_6              =       KF_7 
KD_6a              =      0.05421254005 
KF_6a              =       0.2691115209 
KD_6b             =       KD_6 
KF_6b             =       KF_6 
KD_6c             =       KD_6 
KF_6c             =       KF_6 
KD_7a              =       0.2730363268 
KF_7a              =     0.003675365321 
KD_7b             =       KD_7 
KF_7b             =       KF_7 
KD_7c             =       KD_7 
KF_7c             =       KF_7 
KD_9               =      -0.3144176426 
KF_9               =        0.314601847 
KD_8              =       KD_9 
KF_8              =       KF_9 
KD_8a              =      0.04587795267 
KF_8a              =       0.1783200006 
KD_8b             =       KD_8 
KF_8b             =       KF_8 
KD_8c             =       KD_8 
KF_8c             =       KF_8 
KD_11              =       -0.313530424 
KF_11              =       0.3137143598 
TUNE_WGHT          =               10.0 
DEL_MUX            =   -5.512180934E-08 
DEL_MUY            =   -3.901123335E-08 
#
# IP6 solenoid and quadrupoles strengths
MUX_IP6            =        6.352125146 
MUY_IP6            =        4.018671283 
KQSS1_5            =      -0.4634519914 
KQSS2_5            =      0.03305713451 
KQSS3_5            =       0.3758491433 
KQSS4_5            =       0.0809412749 
KQSS5_5            =      -0.5222748076 
KQSS1_6            =      -0.3555522202 
KQSS2_6            =    0.0005099694535 
KQSS3_6            =       0.2853875397 
KQSS4_6            =      0.07662062675 
KQSS5_6            =       -0.536478377 
KFF1_5             =       0.3443756116 
KFF2_5             =      -0.2380503976 
KFF3_5             =       0.3775796003 
KFF4_5             =      -0.1949412053 
KFF5_5             =        0.187251159 
KFF6_5             =     -0.08489389511 
KFF1_6             =        0.365384125 
KFF2_6             =      -0.2037488048 
KFF3_6             =       0.4082525179 
KFF4_6             =      -0.2136492148 
KFF5_6             =       0.1945703088 
KFF6_6             =      -0.1017285434 
KQLS1_5            =           0.448239 
KQLS2_5            =           -0.36553 
KQLS3_5            =           0.251782 
KQLS4_5            =           0.397172 
KQLS5_5            =           0.251782 
KQLS6_5            =           -0.36553 
KQLS7_5            =           0.448239 
KQLS1_6            =           0.448239 
KQLS2_6            =           -0.36553 
KQLS3_6            =           0.251782 
KQLS4_6            =           0.397172 
KQLS5_6            =           0.251782 
KQLS6_6            =           -0.36553 
KQLS7_6            =           0.448239 
K0EF_6             =       -0.218192315 
K1EF_6             =       0.1000859995 
K2EF_6             =     -0.07457307348 
K3EF_6             =       0.1367818102 
K4EF_6             =     -0.06677557789 
K5EF_6             =    -0.004101903809 
K6EF_6             =     -0.08410656585 
K7EF_6             =        0.195769614 
K8EF_6             =      -0.1990558344 
K9EF_6             =        0.345229068 
K10EF_6            =      -0.1803467941 
K11EF_6            =       0.2264235945 
K12EF_6            =      -0.1996130663 
K13EF_6            =       0.1615485113 
K14EF_6            =      -0.2208787114 
K1ER_6             =      -0.2278853772 
K2ER_6             =       0.2201156485 
K3ER_6             =       0.1788259727 
K4ER_6             =    -0.007637024497 
K5ER_6             =      -0.1625838028 
K6ER_6             =       0.1774622637 
K7ER_6             =      -0.1354227061 
K8ER_6             =      0.06050609079 
K9ER_6             =       0.1006385785 
K10ER_6            =     -0.03463797301 
K11ER_6            =      -0.1024445065 
K12ER_6            =       0.2282592047 
K13ER_6            =      -0.2580290222 
K14ER_6            =       0.1768808343 
K15ER_6            =      -0.2188168368 
KSOL1_6 = phi1_6/(1.0 + A_ELECTRON)/(2.0*LSOL5)
KSOL2_6 = phi2_6/(1.0 + A_ELECTRON)/(2.0*LSOL20)
#
# IP8 solenoid and quadrupoles strengths
MUX_IP8            =        6.508637025 
MUY_IP8            =        4.587520929 
KQSS1_7            =      -0.2230165281 
KQSS2_7            =      -0.0582811745 
KQSS3_7            =     -0.01193738829 
KQSS4_7            =       0.2137457745 
KQSS5_7            =      -0.1844120558 
KQSS1_8            =      -0.3376839281 
KQSS2_8            =      0.05070009734 
KQSS3_8            =       0.2420319659 
KQSS4_8            =       0.1179833172 
KQSS5_8            =      -0.5909929783 
KFF1_7             =        0.102700806 
KFF2_7             =      -0.1712617411 
KFF3_7             =       0.2272981814 
KFF4_7             =      -0.1936954915 
KFF5_7             =       0.2451513988 
KFF6_7             =     -0.08489389511 
KFF1_8             =       0.4256160053 
KFF2_8             =      -0.2495845973 
KFF3_8             =       0.4057378695 
KFF4_8             =      -0.2618177996 
KFF5_8             =       0.2665763448 
KFF6_8             =      -0.1328195162 
KQLS1_7            =           0.448239 
KQLS2_7            =           -0.36553 
KQLS3_7            =           0.251782 
KQLS4_7            =           0.397172 
KQLS5_7            =           0.251782 
KQLS6_7            =           -0.36553 
KQLS7_7            =           0.448239 
KQLS1_8            =           0.448239 
KQLS2_8            =           -0.36553 
KQLS3_8            =           0.251782 
KQLS4_8            =           0.397172 
KQLS5_8            =           0.251782 
KQLS6_8            =           -0.36553 
KQLS7_8            =           0.448239 
K0EF_8             =      -0.2273064537 
K1EF_8             =       0.1019244209 
K2EF_8             =     -0.08190548605 
K3EF_8             =       0.1463303236 
K4EF_8             =      -0.1397024576 
K5EF_8             =       0.1496742106 
K6EF_8             =      -0.2442159358 
K7EF_8             =      0.04393988393 
K8EF_8             =      -0.2831439398 
K9EF_8             =       0.3083359637 
K10EF_8            =      -0.2708898178 
K11EF_8            =       0.3406620038 
K12EF_8            =      -0.2434214825 
K13EF_8            =        0.298202212 
K14EF_8            =      -0.1542741001 
K15EF_8            =       0.1235640461 
K1ER_8             =      -0.2245628613 
K2ER_8             =       0.2145025935 
K3ER_8             =      0.01485241877 
K4ER_8             =       0.1817930759 
K5ER_8             =      -0.1384591674 
K6ER_8             =      0.05486095398 
K7ER_8             =      0.09864667925 
K8ER_8             =      -0.1628913735 
K9ER_8             =       0.1079272595 
K10ER_8            =       0.1530487875 
K11ER_8            =      -0.2224981329 
K12ER_8            =       0.1781177892 
K13ER_8            =        0.119497545 
K14ER_8            =      -0.1280021585 
K15ER_8            =      0.02964357365 
KSOL1_8 = phi1_8/(1.0 + A_ELECTRON)/(2.0*LSOL5)
KSOL2_8 = phi2_8/(1.0 + A_ELECTRON)/(2.0*LSOL20)
#
# IP10 quadrupoles strengths
MUX_IP10           =        4.112371413 
MUY_IP10           =        4.097089146 
KM22_9             =      -0.1185051943 
KM21_9             =      -0.1323371828 
KM20_9             =       0.2082883818 
KM19_9             =      -0.1797478437 
KM18_9             =      0.05902104803 
KM17_9             =      0.03551808764 
KM16_9             =     -0.09344124638 
KM15_9             =       0.1424311997 
KM14_9             =       -0.145625108 
KM13_9             =       0.2169778696 
KM12_9             =      -0.1221746876 
KDSS_10            =      -0.1566717693 
KFSS_10            =       0.1559175595 
KDSSL_10           =      -0.1170915216 
KFSSL_10           =       0.1360968888 
KM22_10            =      -0.2809701903 
KM21_10            =       0.1651338464 
KM20_10            =      -0.2581690321 
KM19_10            =        0.297798994 
KM18_10            =       -0.262909096 
KM17_10            =      -0.1130247262 
KM16_10            =        0.235462357 
KM15_10            =       -0.335725079 
KM14_10            =       0.2775709987 
KM13_10            =      -0.2444138469 
KM12_10            =       0.1312382732 
#
# IP12 quadrupoles strengths
MUX_IP12           =        3.674982242 
MUY_IP12           =        4.232006539 
KM22_11            =      -0.3359731738 
KM21_11            =       0.1816358361 
KM20_11            =      -0.1485365183 
KM19_11            =       0.2154916461 
KM18_11            =      -0.1500435633 
KM17_11            =     0.003291249802 
KM16_11            =    -0.002162170431 
KM15_11            =     0.004247583095 
KM14_11            =      -0.1327720571 
KM13_11            =       0.2102817626 
KM12_11            =     -0.09862623398 
KDSS_12            =      -0.1483172335 
KFSS_12            =       0.1462599289 
KM22_12            =      -0.3359850684 
KM21_12            =       0.1908516877 
KM20_12            =      -0.2708760985 
KM19_12            =       0.2857040226 
KM18_12            =      -0.3359676021 
KM17_12            =       0.1079293418 
KM16_12            =      -0.3359748537 
KM15_12            =       0.2954917251 
KM14_12            =      -0.1682943092 
#
# IP2 QUADRUPOLES STRENGTHS
MUX_IP2            =        3.062275742 
MUY_IP2            =        2.734790901 
KM22_1             =     -0.02701644447 
KM21_1             =      0.01425207743 
KM20_1             =       0.2030463977 
KM19_1             =      -0.2811627068 
KM18_1             =      0.03017571658 
KM17_1             =       0.1073519579 
KM16_1             =        0.141150677 
KM15_1             =      -0.1756395994 
KM14_1             =       0.2831036567 
KM13_1             =      -0.1648119071 
KM12_1             =        0.140091896 
KDSS_2             =       -0.075634279 
KFSS_2             =      0.06204145252 
KM22_2             =      -0.0628929574 
KM21_2             =        0.335817511 
KM20_2             =      -0.1853975692 
KM19_2             =      -0.2570712752 
KM18_2             =       0.2829300672 
KM17_2             =       0.1274757467 
KM16_2             =       -0.239794824 
KM15_2             =       0.1864067227 
KM14_2             =      0.09504800846 
KM13_2             =      -0.1665854015 
KM12_2             =      0.05881784823 
#
# IP4 QUADRUPOLES STRENGTHS
MUX_IP4            =        2.369616463 
MUY_IP4            =         2.46985865 
KD17_3             =      -0.2877048243 
KF16_3             =       0.1824243838 
KD15_3             =     0.008578802751 
KF14_3             =      -0.1480150886 
KD13_3             =      0.07877982097 
KF12_3             =      0.05174533641 
KD11_3             =      0.02794082844 
KF10_3             =     -0.06070609411 
KD9_3              =     -0.06354608493 
KF8_3              =       0.0794935845 
KD7_3              =      0.05008828685 
KF6_3              =       -0.057826039 
KD5_3              =     -0.05081367491 
KF4_3              =      0.05382316298 
KD3_3              =      0.03358142199 
KF2_3              =      -0.1105862703 
KD1_3              =      0.07357421792 
KF18_4             =                0.0 
KD17_4             =      -0.2013945475 
KF16_4             =      0.09127406468 
KD15_4             =      -0.1150529702 
KF14_4             =       0.1861642754 
KD13_4             =     -0.08356348977 
KF12_4             =     -0.07155545708 
KD11_4             =       0.0933824065 
KF10_4             =      0.06520333211 
KD9_4              =     -0.08791343748 
KF8_4              =   -0.0003306038141 
KD7_4              =      0.07175145111 
KF6_4              =     -0.03914277884 
KD5_4              =     -0.05314628092 
KF4_4              =      0.04838351833 
KD3_4              =      0.04110636443 
KF2_4              =      -0.1041863504 
KD1_4              =      0.04845280893 

# KSEXT STRENGTHS
KSFM1_1            =                  0.0 
KSF00_1            =                  0.0 
KSD00_1            =                  0.0 
KSF01_1            =        3.120884289 
KSD01_1            =       -4.294201621 
KSF02_1            =        3.120884289 
KSD02_1            =       -4.294201621 
KSF03_1            =        3.120884289 
KSD03_1            =       -4.294201621 
KSF04_1            =        3.120884289 
KSD04_1            =       -4.294201621 
KSF05_1            =        3.120884289 
KSD05_1            =       -4.294201621 
KSF06_1            =        3.120884289 
KSD06_1            =       -4.294201621 
KSF07_1            =        3.120884289 
KSD07_1            =       -4.294201621 
KSF08_1            =        3.120884289 
KSD08_1            =       -4.294201621 
KSF09_1            =        3.120884289 
KSD09_1            =       -4.294201621 
KSF10_1            =        3.120884289 
KSD10_1            =       -4.294201621 
KSF11_1            =        3.120884289 
KSD11_1            =       -4.294201621 
KSF12_1            =        3.120884289 
KSD12_1            =       -4.294201621 
KSF13_1            =        3.120884289 
KSD13_1            =       -4.294201621 
KSF14_1            =        3.120884289 
KSD14_1            =       -4.294201621 
KSF15_1            =        3.120884289 
KSD15_1            =       -4.294201621 
KSF16_1            =        3.120884289 
KSD16_1            =       -4.294201621 
KSF17_1            =                  0.0 
KSD17_1            =                  0.0
KSF00_3            =                  0.0
KSD00_3            =                  0.0
KSF01_3            =        3.120884289 
KSD01_3            =       -4.294201621 
KSF02_3            =        3.120884289 
KSD02_3            =       -4.294201621 
KSF03_3            =        3.120884289 
KSD03_3            =       -4.294201621 
KSF04_3            =        3.120884289 
KSD04_3            =       -4.294201621 
KSF05_3            =        3.120884289 
KSD05_3            =       -4.294201621 
KSF06_3            =        3.120884289 
KSD06_3            =       -4.294201621 
KSF07_3            =        3.120884289 
KSD07_3            =       -4.294201621 
KSF08_3            =        3.120884289 
KSD08_3            =       -4.294201621 
KSF09_3            =        3.120884289 
KSD09_3            =       -4.294201621 
KSF10_3            =        3.120884289 
KSD10_3            =       -4.294201621 
KSF11_3            =        3.120884289 
KSD11_3            =       -4.294201621 
KSF12_3            =        3.120884289 
KSD12_3            =       -4.294201621 
KSF13_3            =        3.120884289 
KSD13_3            =       -4.294201621 
KSF14_3            =        3.120884289 
KSD14_3            =       -4.294201621 
KSF15_3            =        3.120884289 
KSD15_3            =       -4.294201621 
KSF16_3            =        3.120884289 
KSD16_3            =       -4.294201621 
KSF17_3            =                  0.0 
KSD17_3            =                  0.0 
KSF18_3            =                  0.0 
KSFM1_5            =                  0.0 
KSF00_5            =                  0.0 
KSD00_5            =                  0.0 
KSF01_5            =        3.120884289 
KSD01_5            =       -4.294201621 
KSF02_5            =        3.120884289 
KSD02_5            =       -4.294201621 
KSF03_5            =        3.120884289 
KSD03_5            =       -4.294201621 
KSF04_5            =        3.120884289 
KSD04_5            =       -4.294201621 
KSF05_5            =        3.120884289 
KSD05_5            =       -4.294201621 
KSF06_5            =        3.120884289 
KSD06_5            =       -4.294201621 
KSF07_5            =        3.120884289 
KSD07_5            =       -4.294201621 
KSF08_5            =        3.120884289 
KSD08_5            =       -4.294201621 
KSF09_5            =        3.120884289 
KSD09_5            =       -4.294201621 
KSF10_5            =        3.120884289 
KSD10_5            =       -4.294201621 
KSF11_5            =        3.120884289 
KSD11_5            =       -4.294201621 
KSF12_5            =        3.120884289 
KSD12_5            =       -4.294201621 
KSF13_5            =        3.120884289 
KSD13_5            =       -4.294201621 
KSF14_5            =        3.120884289 
KSD14_5            =       -4.294201621 
KSF15_5            =        3.120884289 
KSD15_5            =       -4.294201621 
KSF16_5            =        3.120884289 
KSD16_5            =       -4.294201621 
KSF01_7            =        3.120884289 
KSD01_7            =       -4.294201621 
KSF02_7            =        3.120884289 
KSD02_7            =       -4.294201621 
KSF03_7            =        3.120884289 
KSD03_7            =       -4.294201621 
KSF04_7            =        3.120884289 
KSD04_7            =       -4.294201621 
KSF05_7            =        3.120884289 
KSD05_7            =       -4.294201621 
KSF06_7            =        3.120884289 
KSD06_7            =       -4.294201621 
KSF07_7            =        3.120884289 
KSD07_7            =       -4.294201621 
KSF08_7            =        3.120884289 
KSD08_7            =       -4.294201621 
KSF09_7            =        3.120884289 
KSD09_7            =       -4.294201621 
KSF10_7            =        3.120884289 
KSD10_7            =       -4.294201621 
KSF11_7            =        3.120884289 
KSD11_7            =       -4.294201621 
KSF12_7            =        3.120884289 
KSD12_7            =       -4.294201621 
KSF13_7            =        3.120884289 
KSD13_7            =       -4.294201621 
KSF14_7            =        3.120884289 
KSD14_7            =       -4.294201621 
KSF15_7            =        3.120884289 
KSD15_7            =       -4.294201621 
KSF16_7            =        3.120884289 
KSD16_7            =       -4.294201621 
KSF01_9            =        3.120884289 
KSD01_9            =       -4.294201621 
KSF02_9            =        3.120884289 
KSD02_9            =       -4.294201621 
KSF03_9            =        3.120884289 
KSD03_9            =       -4.294201621 
KSF04_9            =        3.120884289 
KSD04_9            =       -4.294201621 
KSF05_9            =        3.120884289 
KSD05_9            =       -4.294201621 
KSF06_9            =        3.120884289 
KSD06_9            =       -4.294201621 
KSF07_9            =        3.120884289 
KSD07_9            =       -4.294201621 
KSF08_9            =        3.120884289 
KSD08_9            =       -4.294201621 
KSF09_9            =        3.120884289 
KSD09_9            =       -4.294201621 
KSF10_9            =        3.120884289 
KSD10_9            =       -4.294201621 
KSF11_9            =        3.120884289 
KSD11_9            =       -4.294201621 
KSF12_9            =        3.120884289 
KSD12_9            =       -4.294201621 
KSF13_9            =        3.120884289 
KSD13_9            =       -4.294201621 
KSF14_9            =        3.120884289 
KSD14_9            =       -4.294201621 
KSF15_9            =        3.120884289 
KSD15_9            =       -4.294201621 
KSF16_9            =        3.120884289 
KSD16_9            =       -4.294201621 
KSF17_9            =                  0.0 
KSD17_9            =                  0.0 
KSF00_11           =                  0.0 
KSD00_11           =                  0.0 
KSF01_11           =        3.120884289 
KSD01_11           =       -4.294201621 
KSF02_11           =        3.120884289 
KSD02_11           =       -4.294201621 
KSF03_11           =        3.120884289 
KSD03_11           =       -4.294201621 
KSF04_11           =        3.120884289 
KSD04_11           =       -4.294201621 
KSF05_11           =        3.120884289 
KSD05_11           =       -4.294201621 
KSF06_11           =        3.120884289 
KSD06_11           =       -4.294201621 
KSF07_11           =        3.120884289 
KSD07_11           =       -4.294201621 
KSF08_11           =        3.120884289 
KSD08_11           =       -4.294201621 
KSF09_11           =        3.120884289 
KSD09_11           =       -4.294201621 
KSF10_11           =        3.120884289 
KSD10_11           =       -4.294201621 
KSF11_11           =        3.120884289 
KSD11_11           =       -4.294201621 
KSF12_11           =        3.120884289 
KSD12_11           =       -4.294201621 
KSF13_11           =        3.120884289 
KSD13_11           =       -4.294201621 
KSF14_11           =        3.120884289 
KSD14_11           =       -4.294201621 
KSF15_11           =        3.120884289 
KSD15_11           =       -4.294201621 
KSF16_11           =        3.120884289 
KSD16_11           =       -4.294201621 
KSF17_11           =                  0.0 
KSD17_11           =                  0.0 
KSF18_11           =                  0.0 
S41_2              =                  0.0 
S42_2              =                  0.0 
S43_2              =                  0.0 
S44_2              =                  0.0 
S45_2              =                  0.0 
S46_2              =                  0.0 
S47_2              =                  0.0 
S48_2              =                  0.0 
S49_2              =                  0.0 
S50_2              =                  0.0 
S51_2              =                  0.0 
S52_2              =                  0.0 
S41_4              =                  0.0 
S42_4              =                  0.0 
S43_4              =                  0.0 
S44_4              =                  0.0 
S45_4              =                  0.0 
S46_4              =                  0.0 
S47_4              =                  0.0 
S48_4              =                  0.0 
S49_4              =                  0.0 
S50_4              =                  0.0 
S51_4              =                  0.0 
S52_4              =                  0.0
 


function CHICANE(ANGLE_ARC, TH2XX, ID) 
   global L23_STR, L01_STR, L12_STR
   if 1.0 < TH2XX 
      TH2_ID = asin(L23_STR*sin(0.5*(ANGLE_ARC))/(2.0*L01_STR + L23_STR))
   else
      TH2_ID = TH2XX
   end
   L01_ID = L01_STR*(0.5*(ANGLE_ARC) - (TH2_ID))/(sin(0.5*(ANGLE_ARC)) - sin(TH2_ID))
   L23_ID = L23_STR*(TH2_ID)/sin(TH2_ID)
   L12_ID = L12_STR/cos(TH2_ID)

   OSB12_ID = DRIFT(name="OSB12_"*ID, len=L12_ID)
   EDGE1_ID = thinMULTIPOLE(name="EDGE1_"*ID, PolynomB=[0.0, -tan(0.5*(ANGLE_ARC))* (0.5*(ANGLE_ARC)-(TH2_ID)) / L01_ID, 0.0, 0.0])
   EDGE2_ID = thinMULTIPOLE(name="EDGE2_"*ID, PolynomB=[0.0, -tan(-TH2_ID)* (0.5*(ANGLE_ARC)-(TH2_ID)) / L01_ID, 0.0, 0.0])
   EDGE3_ID = thinMULTIPOLE(name="EDGE3_"*ID, PolynomB=[0.0, -tan(TH2_ID)* (2*(TH2_ID)) / L23_ID, 0.0, 0.0])

   D01A_ID = SBEND(name="D01A_"*ID, len=L01_ID, angle=0.5*(ANGLE_ARC)-(TH2_ID), rad=1)
   D23_ID = SBEND(name="D23_"*ID, len=L23_ID, angle=2*(TH2_ID), rad=1)
   D01B_ID = SBEND(name="D01B_"*ID, len=L01_ID, angle=0.5*(ANGLE_ARC)-(TH2_ID), rad=1)
   return [EDGE1_ID, D01A_ID, EDGE2_ID, OSB12_ID, EDGE3_ID, D23_ID, 
      EDGE3_ID, OSB12_ID, EDGE2_ID, D01B_ID, EDGE1_ID]
end
DCCN = CHICANE(ARC_BEND_ANGLE, TH2_DCCN, "000")
DMNS_6 = CHICANE(ARC_ANG_M6, TH2_DCCN, "001")
DPLS_6 = CHICANE(ARC_ANG_P6, TH2_DCCN, "002")
DPLS_8 = CHICANE(ARC_ANG_P8, TH2_DCCN, "003")
DMNS_8 = CHICANE(ARC_ANG_M8, TH2_DCCN, "004")
# DCCN = RBEND(name="DCCN", len=LBEND, angle=ARC_BEND_ANGLE)
# DMNS_6 = RBEND(name="DMNS_6", len=LBEND, angle=ARC_BEND_ANGLE - IP6_ANGLE/8.0)
# DPLS_6 = RBEND(name="DPLS_6", len=LBEND, angle=ARC_BEND_ANGLE + IP6_ANGLE/8.0)
# DPLS_8 = RBEND(name="DPLS_8", len=LBEND, angle=ARC_BEND_ANGLE + IP8_ANGLE/8.0)
# DMNS_8 = RBEND(name="DMNS_8", len=LBEND, angle=ARC_BEND_ANGLE - IP8_ANGLE/8.0)

# DB23_1 = RBEND(name="DB23_1", len=L01_STR, angle=ROT_ANGLE, rad=1)
# DB23_2 = RBEND(name="DB23_2", len=L01_STR, angle=ROT_ANGLE, rad=1)
# DB23_3 = RBEND(name="DB23_3", len=L01_STR, angle=ROT_ANGLE, rad=1)
# DB23_4 = RBEND(name="DB23_4", len=L01_STR, angle=ROT_ANGLE, rad=1)
# DB23_5 = RBEND(name="DB23_5", len=3.8, angle=theta1/5, rad=1)
# DB23_6 = RBEND(name="DB23_6", len=3.8, angle=theta1/5, rad=1)
# DB23_7 = RBEND(name="DB23_7", len=3.8, angle=theta1/5, rad=1)
# DB23_8 = RBEND(name="DB23_8", len=3.8, angle=theta1/5, rad=1)
# DB23_9 = RBEND(name="DB23_9", len=L01_STR, angle=ROT_ANGLE, rad=1)
# DB23_10 = RBEND(name="DB23_10", len=L01_STR, angle=ROT_ANGLE, rad=1)
# DB23_11 = RBEND(name="DB23_11", len=L01_STR, angle=ROT_ANGLE, rad=1)
# DB23_12 = RBEND(name="DB23_12", len=L01_STR, angle=ROT_ANGLE, rad=1)
# DB12P = RBEND(name="DB12P", len=L01_STR, angle=ADB12_4, rad=1)
# DB12M = RBEND(name="DB12M", len=L01_STR, angle=-ADB12_4, rad=1)
# DB4P = RBEND(name="DB4P", len=L01_STR, angle=ADB12_4, rad=1)
# DB4M = RBEND(name="DB4M", len=L01_STR, angle=-ADB12_4, rad=1)
# d1ef_6 = RBEND(name="d1ef_6", len=L23_STR, angle=db1ef_ang_6, rad=1)
# d2ef_6 = RBEND(name="d2ef_6", len=L01_STR, angle=db2ef_ang_6, rad=1)
# d3ef_6 = RBEND(name="d3ef_6", len=L01_STR, angle=db3ef_ang_6, rad=1)
# d4ef_6 = RBEND(name="d4ef_6", len=L01_STR, angle=db4ef_ang_6, rad=1)
# d5ef_6 = RBEND(name="d5ef_6", len=L01_STR, angle=db5ef_ang_6, rad=1)
# d2er_6 = RBEND(name="d2er_6", len=lrd2er, angle=DB2ER_ANG_6, rad=1)
# d3er_6 = RBEND(name="d3er_6", len=L01_STR, angle=DB3ER_ANG_6, rad=1)
# d4er_6 = RBEND(name="d4er_6", len=L01_STR, angle=DB4ER_ANG_6, rad=1)
# d1ef_8 = RBEND(name="d1ef_8", len=L23_STR, angle=DB1EF_ANG_8, rad=1)
# d2ef_8 = RBEND(name="d2ef_8", len=L01_STR, angle=DB2EF_ANG_8, rad=1)
# d3ef_8 = RBEND(name="d3ef_8", len=3.8, angle=DB3EF_ANG_8, rad=1)
# d2er_8 = RBEND(name="d2er_8", len=lrd2er, angle=DB2ER_ANG_8, rad=1)
# d3er_8 = RBEND(name="d3er_8", len=L01_STR, angle=DB3ER_ANG_8, rad=1)
# d4er_8 = RBEND(name="d4er_8", len=L01_STR, angle=DB4ER_ANG_8, rad=1)
# d5er_8 = RBEND(name="d5er_8", len=L01_STR, angle=DB5ER_ANG_8, rad=1)
# MADX RBEND length is the chord length, not the arc length
DB23_1 = RBEND(name="DB23_1", len=L01_STR*ROT_ANGLE/2/sin(ROT_ANGLE/2), angle=ROT_ANGLE, rad=1)
DB23_2 = RBEND(name="DB23_2", len=L01_STR*ROT_ANGLE/2/sin(ROT_ANGLE/2), angle=ROT_ANGLE, rad=1)
DB23_3 = RBEND(name="DB23_3", len=L01_STR*ROT_ANGLE/2/sin(ROT_ANGLE/2), angle=ROT_ANGLE, rad=1)
DB23_4 = RBEND(name="DB23_4", len=L01_STR*ROT_ANGLE/2/sin(ROT_ANGLE/2), angle=ROT_ANGLE, rad=1)
DB23_5 = RBEND(name="DB23_5", len=3.8*theta1/5/2/sin(theta1/5/2), angle=theta1/5, rad=1)
DB23_6 = RBEND(name="DB23_6", len=3.8*theta1/5/2/sin(theta1/5/2), angle=theta1/5, rad=1)
DB23_7 = RBEND(name="DB23_7", len=3.8*theta1/5/2/sin(theta1/5/2), angle=theta1/5, rad=1)
DB23_8 = RBEND(name="DB23_8", len=3.8*theta1/5/2/sin(theta1/5/2), angle=theta1/5, rad=1)
DB23_9 = RBEND(name="DB23_9", len=L01_STR*ROT_ANGLE/2/sin(ROT_ANGLE/2), angle=ROT_ANGLE, rad=1)
DB23_10 = RBEND(name="DB23_10", len=L01_STR*ROT_ANGLE/2/sin(ROT_ANGLE/2), angle=ROT_ANGLE, rad=1)
DB23_11 = RBEND(name="DB23_11", len=L01_STR*ROT_ANGLE/2/sin(ROT_ANGLE/2), angle=ROT_ANGLE, rad=1)
DB23_12 = RBEND(name="DB23_12", len=L01_STR*ROT_ANGLE/2/sin(ROT_ANGLE/2), angle=ROT_ANGLE, rad=1)
DB12P = RBEND(name="DB12P", len=L01_STR*ADB12_4/2/sin(ADB12_4/2), angle=ADB12_4, rad=1)
DB12M = RBEND(name="DB12M", len=L01_STR*ADB12_4/2/sin(ADB12_4/2), angle=-ADB12_4, rad=1)
DB4P = RBEND(name="DB4P", len=L01_STR*ADB12_4/2/sin(ADB12_4/2), angle=ADB12_4, rad=1)
DB4M = RBEND(name="DB4M", len=L01_STR*ADB12_4/2/sin(ADB12_4/2), angle=-ADB12_4, rad=1)
d1ef_6 = RBEND(name="d1ef_6", len=L23_STR*db1ef_ang_6/2/sin(db1ef_ang_6/2), angle=db1ef_ang_6, rad=1)
d2ef_6 = RBEND(name="d2ef_6", len=L01_STR*db2ef_ang_6/2/sin(db2ef_ang_6/2), angle=db2ef_ang_6, rad=1)
d3ef_6 = RBEND(name="d3ef_6", len=L01_STR*db3ef_ang_6/2/sin(db3ef_ang_6/2), angle=db3ef_ang_6, rad=1)
d4ef_6 = RBEND(name="d4ef_6", len=L01_STR*db4ef_ang_6/2/sin(db4ef_ang_6/2), angle=db4ef_ang_6, rad=1)
d5ef_6 = RBEND(name="d5ef_6", len=L01_STR*db5ef_ang_6/2/sin(db5ef_ang_6/2), angle=db5ef_ang_6, rad=1)
d2er_6 = RBEND(name="d2er_6", len=lrd2er*DB2ER_ANG_6/2/sin(DB2ER_ANG_6/2), angle=DB2ER_ANG_6, rad=1)
d3er_6 = RBEND(name="d3er_6", len=L01_STR*DB3ER_ANG_6/2/sin(DB3ER_ANG_6/2), angle=DB3ER_ANG_6, rad=1)
d4er_6 = RBEND(name="d4er_6", len=L01_STR*DB4ER_ANG_6/2/sin(DB4ER_ANG_6/2), angle=DB4ER_ANG_6, rad=1)
d1ef_8 = RBEND(name="d1ef_8", len=L23_STR*DB1EF_ANG_8/2/sin(DB1EF_ANG_8/2), angle=DB1EF_ANG_8, rad=1)
d2ef_8 = RBEND(name="d2ef_8", len=L01_STR*DB2EF_ANG_8/2/sin(DB2EF_ANG_8/2), angle=DB2EF_ANG_8, rad=1)
d3ef_8 = RBEND(name="d3ef_8", len=3.8*DB3EF_ANG_8/2/sin(DB3EF_ANG_8/2), angle=DB3EF_ANG_8, rad=1)
d2er_8 = RBEND(name="d2er_8", len=lrd2er*DB2ER_ANG_8/2/sin(DB2ER_ANG_8/2), angle=DB2ER_ANG_8, rad=1)
d3er_8 = RBEND(name="d3er_8", len=L01_STR*DB3ER_ANG_8/2/sin(DB3ER_ANG_8/2), angle=DB3ER_ANG_8, rad=1)
d4er_8 = RBEND(name="d4er_8", len=L01_STR*DB4ER_ANG_8/2/sin(DB4ER_ANG_8/2), angle=DB4ER_ANG_8, rad=1)
d5er_8 = RBEND(name="d5er_8", len=L01_STR*DB5ER_ANG_8/2/sin(DB5ER_ANG_8/2), angle=DB5ER_ANG_8, rad=1)
#
# ======= The quadrupoles =================================
#
# ------- The arc quadrupoles -----------------------------
# HQF_AVE = KQUAD(name="HQF_AVE", len=LQ50/2.0, k1=KF_AVE)
# HQD_AVE = KQUAD(name="HQD_AVE", len=LQ50/2.0, k1=KD_AVE)

HQF_1 = KQUAD(name="HQF_1", len=LQ50/2.0, k1=KF_1, rad=1)
HQD_1 = KQUAD(name="HQD_1", len=LQ50/2.0, k1=KD_1, rad=1)

HQF_3 = KQUAD(name="HQF_3", len=LQ50/2.0, k1=KF_3, rad=1)
HQD_3 = KQUAD(name="HQD_3", len=LQ50/2.0, k1=KD_3, rad=1)

HQF_5 = KQUAD(name="HQF_5", len=LQ50/2.0, k1=KF_5, rad=1)
HQD_5 = KQUAD(name="HQD_5", len=LQ50/2.0, k1=KD_5, rad=1)
HQF_5a = KQUAD(name="HQF_5a", len=LQ50/2.0, k1=KF_5a, rad=1)
HQD_5a = KQUAD(name="HQD_5a", len=LQ50/2.0, k1=KD_5a, rad=1)
HQF_5b = KQUAD(name="HQF_5b", len=LQ50/2.0, k1=KF_5b, rad=1)
HQD_5b = KQUAD(name="HQD_5b", len=LQ50/2.0, k1=KD_5b, rad=1)
HQF_5c = KQUAD(name="HQF_5c", len=LQ50/2.0, k1=KF_5c, rad=1)
HQD_5c = KQUAD(name="HQD_5c", len=LQ50/2.0, k1=KD_5c, rad=1)

HQF_6a = KQUAD(name="HQF_6a", len=LQ50/2.0, k1=KF_6a, rad=1)
HQD_6a = KQUAD(name="HQD_6a", len=LQ50/2.0, k1=KD_6a, rad=1)
HQF_6b = KQUAD(name="HQF_6b", len=LQ50/2.0, k1=KF_6b, rad=1)
HQD_6b = KQUAD(name="HQD_6b", len=LQ50/2.0, k1=KD_6b, rad=1)
HQF_6c = KQUAD(name="HQF_6c", len=LQ50/2.0, k1=KF_6c, rad=1)
HQD_6c = KQUAD(name="HQD_6c", len=LQ50/2.0, k1=KD_6c, rad=1)

HQF_7 = KQUAD(name="HQF_7", len=LQ50/2.0, k1=KF_7, rad=1)
HQD_7 = KQUAD(name="HQD_7", len=LQ50/2.0, k1=KD_7, rad=1)
HQF_7a = KQUAD(name="HQF_7a", len=LQ50/2.0, k1=KF_7a, rad=1)
HQD_7a = KQUAD(name="HQD_7a", len=LQ50/2.0, k1=KD_7a, rad=1)
HQF_7b = KQUAD(name="HQF_7b", len=LQ50/2.0, k1=KF_7b, rad=1)
HQD_7b = KQUAD(name="HQD_7b", len=LQ50/2.0, k1=KD_7b, rad=1)
HQF_7c = KQUAD(name="HQF_7c", len=LQ50/2.0, k1=KF_7c, rad=1)
HQD_7c = KQUAD(name="HQD_7c", len=LQ50/2.0, k1=KD_7c, rad=1)

HQF_8a = KQUAD(name="HQF_8a", len=LQ50/2.0, k1=KF_8a, rad=1)
HQD_8a = KQUAD(name="HQD_8a", len=LQ50/2.0, k1=KD_8a, rad=1)
HQF_8b = KQUAD(name="HQF_8b", len=LQ50/2.0, k1=KF_8b, rad=1)
HQD_8b = KQUAD(name="HQD_8b", len=LQ50/2.0, k1=KD_8b, rad=1)
HQF_8c = KQUAD(name="HQF_8c", len=LQ50/2.0, k1=KF_8c, rad=1)
HQD_8c = KQUAD(name="HQD_8c", len=LQ50/2.0, k1=KD_8c, rad=1)

HQF_9 = KQUAD(name="HQF_9", len=LQ50/2.0, k1=KF_9, rad=1)
HQD_9 = KQUAD(name="HQD_9", len=LQ50/2.0, k1=KD_9, rad=1)

# HQF_9a = KQUAD(name="HQF_9a", len=LQ50/2.0, k1=KF_9a)
# HQD_9a = KQUAD(name="HQD_9a", len=LQ50/2.0, k1=KD_9a)
# HQF_9b = KQUAD(name="HQF_9b", len=LQ50/2.0, k1=KF_9b)
# HQD_9b = KQUAD(name="HQD_9b", len=LQ50/2.0, k1=KD_9b)
# HQF_9c = KQUAD(name="HQF_9c", len=LQ50/2.0, k1=KF_9c)
# HQD_9c = KQUAD(name="HQD_9c", len=LQ50/2.0, k1=KD_9c)

HQF_11 = KQUAD(name="HQF_11", len=LQ50/2.0, k1=KF_11, rad=1)
HQD_11 = KQUAD(name="HQD_11", len=LQ50/2.0, k1=KD_11, rad=1)
#
# IR quads.
qir06 = KQUAD(name="qir06", len=lqir06)
qir12 = KQUAD(name="qir12", len=lqir12)
#
# IP6 quadrupoles.
HQSS1_5 = KQUAD(name="HQSS1_5", len=LQSS1/2.0, k1=KQSS1_5, rad=1)
HQSS2_5 = KQUAD(name="HQSS2_5", len=LQSS2/2.0, k1=KQSS2_5, rad=1)
HQSS3_5 = KQUAD(name="HQSS3_5", len=LQSS3/2.0, k1=KQSS3_5, rad=1)
HQSS4_5 = KQUAD(name="HQSS4_5", len=LQSS4/2.0, k1=KQSS4_5, rad=1)
HQSS5_5 = KQUAD(name="HQSS5_5", len=LQSS5/2.0, k1=KQSS5_5, rad=1)
HQSS1_6 = KQUAD(name="HQSS1_6", len=LQSS1/2.0, k1=KQSS1_6, rad=1)
HQSS2_6 = KQUAD(name="HQSS2_6", len=LQSS2/2.0, k1=KQSS2_6, rad=1)
HQSS3_6 = KQUAD(name="HQSS3_6", len=LQSS3/2.0, k1=KQSS3_6, rad=1)
HQSS4_6 = KQUAD(name="HQSS4_6", len=LQSS4/2.0, k1=KQSS4_6, rad=1)
HQSS5_6 = KQUAD(name="HQSS5_6", len=LQSS5/2.0, k1=KQSS5_6, rad=1)
HQLS1_5 = KQUAD(name="HQLS1_5", len=LQLS1/2.0, k1=KQLS1_5, rad=1)
HQLS2_5 = KQUAD(name="HQLS2_5", len=LQLS2/2.0, k1=KQLS2_5, rad=1)
HQLS3_5 = KQUAD(name="HQLS3_5", len=LQLS3/2.0, k1=KQLS3_5, rad=1)
HQLS4_5 = KQUAD(name="HQLS4_5", len=LQLS4/2.0, k1=KQLS4_5, rad=1)
HQLS5_5 = KQUAD(name="HQLS5_5", len=LQLS3/2.0, k1=KQLS5_5, rad=1)
HQLS6_5 = KQUAD(name="HQLS6_5", len=LQLS2/2.0, k1=KQLS6_5, rad=1)
HQLS7_5 = KQUAD(name="HQLS7_5", len=LQLS1/2.0, k1=KQLS7_5, rad=1)
HQLS1_6 = KQUAD(name="HQLS1_6", len=LQLS1/2.0, k1=KQLS1_6, rad=1)
HQLS2_6 = KQUAD(name="HQLS2_6", len=LQLS2/2.0, k1=KQLS2_6, rad=1)
HQLS3_6 = KQUAD(name="HQLS3_6", len=LQLS3/2.0, k1=KQLS3_6, rad=1)
HQLS4_6 = KQUAD(name="HQLS4_6", len=LQLS4/2.0, k1=KQLS4_6, rad=1)
HQLS5_6 = KQUAD(name="HQLS5_6", len=LQLS3/2.0, k1=KQLS5_6, rad=1)
HQLS6_6 = KQUAD(name="HQLS6_6", len=LQLS2/2.0, k1=KQLS6_6, rad=1)
HQLS7_6 = KQUAD(name="HQLS7_6", len=LQLS1/2.0, k1=KQLS7_6, rad=1)
QFF1_5 = KQUAD(name="QFF1_5", len=LQ50, k1=KFF1_5, rad=1)
QFF2_5 = KQUAD(name="QFF2_5", len=LQ50, k1=KFF2_5, rad=1)
QFF3_5 = KQUAD(name="QFF3_5", len=LQ50, k1=KFF3_5, rad=1)
QFF4_5 = KQUAD(name="QFF4_5", len=LQ50, k1=KFF4_5, rad=1)
QFF5_5 = KQUAD(name="QFF5_5", len=LQ50, k1=KFF5_5, rad=1)
QFF6_5 = KQUAD(name="QFF6_5", len=LQ50, k1=KFF6_5, rad=1)
QFF1_6 = KQUAD(name="QFF1_6", len=LQ50, k1=KFF1_6, rad=1)
QFF2_6 = KQUAD(name="QFF2_6", len=LQ50, k1=KFF2_6, rad=1)
QFF3_6 = KQUAD(name="QFF3_6", len=LQ50, k1=KFF3_6, rad=1)
QFF4_6 = KQUAD(name="QFF4_6", len=LQ50, k1=KFF4_6, rad=1)
QFF5_6 = KQUAD(name="QFF5_6", len=LQ50, k1=KFF5_6, rad=1)
QFF6_6 = KQUAD(name="QFF6_6", len=LQ50, k1=KFF6_6, rad=1)
Q0EF_6 = KQUAD(name="Q0EF_6", len=lq0ef_6, k1=K0EF_6, rad=1)
Q1EF_6 = KQUAD(name="Q1EF_6", len=lq1ef_6, k1=K1EF_6, rad=1)
q2ef_6 = KQUAD(name="q2ef_6", len=0.6, k1=K2EF_6, rad=1)
q3ef_6 = KQUAD(name="q3ef_6", len=0.6, k1=K3EF_6, rad=1)
q4ef_6 = KQUAD(name="q4ef_6", len=1.2, k1=K4EF_6, rad=1)
q5ef_6 = KQUAD(name="q5ef_6", len=1.2, k1=K5EF_6, rad=1)
q6ef_6 = KQUAD(name="q6ef_6", len=1.2, k1=K6EF_6, rad=1)
q7ef_6 = KQUAD(name="q7ef_6", len=1.2, k1=K7EF_6, rad=1)
q8ef_6 = KQUAD(name="q8ef_6", len=1.2, k1=K8EF_6, rad=1)
q9ef_6 = KQUAD(name="q9ef_6", len=1.2, k1=K9EF_6, rad=1)
q10ef_6 = KQUAD(name="q10ef_6", len=1.2, k1=K10EF_6, rad=1)
q11ef_6 = KQUAD(name="q11ef_6", len=1.2, k1=K11EF_6, rad=1)
q12ef_6 = KQUAD(name="q12ef_6", len=1.2, k1=K12EF_6, rad=1)
q13ef_6 = KQUAD(name="q13ef_6", len=1.2, k1=K13EF_6, rad=1)
q14ef_6 = KQUAD(name="q14ef_6", len=1.2, k1=K14EF_6, rad=1)
q1er_6 = KQUAD(name="q1er_6", len=lq1er_6, k1=K1ER_6, rad=1)
q2er_6 = KQUAD(name="q2er_6", len=lq2er_6, k1=K2ER_6, rad=1)
q3er_6 = KQUAD(name="q3er_6", len=0.6, k1=K3ER_6, rad=1)
q4er_6 = KQUAD(name="q4er_6", len=0.6, k1=K4ER_6, rad=1)
q5er_6 = KQUAD(name="q5er_6", len=0.6, k1=K5ER_6, rad=1)
q6er_6 = KQUAD(name="q6er_6", len=1.2, k1=K6ER_6, rad=1)
q7er_6 = KQUAD(name="q7er_6", len=1.2, k1=K7ER_6, rad=1)
q8er_6 = KQUAD(name="q8er_6", len=1.2, k1=K8ER_6, rad=1)
q9er_6 = KQUAD(name="q9er_6", len=1.2, k1=K9ER_6, rad=1)
q10er_6 = KQUAD(name="q10er_6", len=1.2, k1=K10ER_6, rad=1)
q11er_6 = KQUAD(name="q11er_6", len=1.2, k1=K11ER_6, rad=1)
q12er_6 = KQUAD(name="q12er_6", len=1.2, k1=K12ER_6, rad=1)
q13er_6 = KQUAD(name="q13er_6", len=1.2, k1=K13ER_6, rad=1)
q14er_6 = KQUAD(name="q14er_6", len=1.2, k1=K14ER_6, rad=1)
q15er_6 = KQUAD(name="q15er_6", len=1.2, k1=K15ER_6, rad=1)
sq = KQUAD(name="sq", len=0.25, rad=1)
sq2ef_6 = KQUAD(name="sq2ef_6", len=0.25, rad=1)
sq3ef_6 = KQUAD(name="sq3ef_6", len=0.25, rad=1)
sq4ef_6 = KQUAD(name="sq4ef_6", len=0.25, rad=1)
sq5ef_6 = KQUAD(name="sq5ef_6", len=0.25, rad=1)
sq6ef_6 = KQUAD(name="sq6ef_6", len=0.25, rad=1)
sq7ef_6 = KQUAD(name="sq7ef_6", len=0.25, rad=1)
sq8ef_6 = KQUAD(name="sq8ef_6", len=0.25, rad=1)
sq9ef_6 = KQUAD(name="sq9ef_6", len=0.25, rad=1)
sq10ef_6 = KQUAD(name="sq10ef_6", len=0.25, rad=1)
sq11ef_6 = KQUAD(name="sq11ef_6", len=0.25, rad=1)
sq12ef_6 = KQUAD(name="sq12ef_6", len=0.25, rad=1)
sq13ef_6 = KQUAD(name="sq13ef_6", len=0.25, rad=1)
sq14ef_6 = KQUAD(name="sq14ef_6", len=0.25, rad=1)
sq3er_6 = KQUAD(name="sq3er_6", len=0.25, rad=1)
sq4er_6 = KQUAD(name="sq4er_6", len=0.25, rad=1)
sq5er_6 = KQUAD(name="sq5er_6", len=0.25, rad=1)
sq6er_6 = KQUAD(name="sq6er_6", len=0.25, rad=1)
sq7er_6 = KQUAD(name="sq7er_6", len=0.25, rad=1)
sq8er_6 = KQUAD(name="sq8er_6", len=0.25, rad=1)
sq9er_6 = KQUAD(name="sq9er_6", len=0.25, rad=1)
sq10er_6 = KQUAD(name="sq10er_6", len=0.25, rad=1)
#
# IP8 quadrupoles.
HQSS1_7 = KQUAD(name="HQSS1_7", len=LQSS1/2.0, k1=KQSS1_7, rad=1)
HQSS2_7 = KQUAD(name="HQSS2_7", len=LQSS2/2.0, k1=KQSS2_7, rad=1)
HQSS3_7 = KQUAD(name="HQSS3_7", len=LQSS3/2.0, k1=KQSS3_7, rad=1)
HQSS4_7 = KQUAD(name="HQSS4_7", len=LQSS4/2.0, k1=KQSS4_7, rad=1)
HQSS5_7 = KQUAD(name="HQSS5_7", len=LQSS5/2.0, k1=KQSS5_7, rad=1)
HQSS1_8 = KQUAD(name="HQSS1_8", len=LQSS1/2.0, k1=KQSS1_8, rad=1)
HQSS2_8 = KQUAD(name="HQSS2_8", len=LQSS2/2.0, k1=KQSS2_8, rad=1)
HQSS3_8 = KQUAD(name="HQSS3_8", len=LQSS3/2.0, k1=KQSS3_8, rad=1)
HQSS4_8 = KQUAD(name="HQSS4_8", len=LQSS4/2.0, k1=KQSS4_8, rad=1)
HQSS5_8 = KQUAD(name="HQSS5_8", len=LQSS5/2.0, k1=KQSS5_8, rad=1)
HQLS1_7 = KQUAD(name="HQLS1_7", len=LQLS1/2.0, k1=KQLS1_7, rad=1)
HQLS2_7 = KQUAD(name="HQLS2_7", len=LQLS2/2.0, k1=KQLS2_7, rad=1)
HQLS3_7 = KQUAD(name="HQLS3_7", len=LQLS3/2.0, k1=KQLS3_7, rad=1)
HQLS4_7 = KQUAD(name="HQLS4_7", len=LQLS4/2.0, k1=KQLS4_7, rad=1)
HQLS5_7 = KQUAD(name="HQLS5_7", len=LQLS3/2.0, k1=KQLS5_7, rad=1)
HQLS6_7 = KQUAD(name="HQLS6_7", len=LQLS2/2.0, k1=KQLS6_7, rad=1)
HQLS7_7 = KQUAD(name="HQLS7_7", len=LQLS1/2.0, k1=KQLS7_7, rad=1)
HQLS1_8 = KQUAD(name="HQLS1_8", len=LQLS1/2.0, k1=KQLS1_8, rad=1)
HQLS2_8 = KQUAD(name="HQLS2_8", len=LQLS2/2.0, k1=KQLS2_8, rad=1)
HQLS3_8 = KQUAD(name="HQLS3_8", len=LQLS3/2.0, k1=KQLS3_8, rad=1)
HQLS4_8 = KQUAD(name="HQLS4_8", len=LQLS4/2.0, k1=KQLS4_8, rad=1)
HQLS5_8 = KQUAD(name="HQLS5_8", len=LQLS3/2.0, k1=KQLS5_8, rad=1)
HQLS6_8 = KQUAD(name="HQLS6_8", len=LQLS2/2.0, k1=KQLS6_8, rad=1)
HQLS7_8 = KQUAD(name="HQLS7_8", len=LQLS1/2.0, k1=KQLS7_8, rad=1)
QFF1_7 = KQUAD(name="QFF1_7", len=LQ50, k1=KFF1_7, rad=1)
QFF2_7 = KQUAD(name="QFF2_7", len=LQ50, k1=KFF2_7, rad=1)
QFF3_7 = KQUAD(name="QFF3_7", len=LQ50, k1=KFF3_7, rad=1)
QFF4_7 = KQUAD(name="QFF4_7", len=LQ50, k1=KFF4_7, rad=1)
QFF5_7 = KQUAD(name="QFF5_7", len=LQ50, k1=KFF5_7, rad=1)
QFF6_7 = KQUAD(name="QFF6_7", len=LQ50, k1=KFF6_7, rad=1)
QFF1_8 = KQUAD(name="QFF1_8", len=LQ50, k1=KFF1_8, rad=1)
QFF2_8 = KQUAD(name="QFF2_8", len=LQ50, k1=KFF2_8, rad=1)
QFF3_8 = KQUAD(name="QFF3_8", len=LQ50, k1=KFF3_8, rad=1)
QFF4_8 = KQUAD(name="QFF4_8", len=LQ50, k1=KFF4_8, rad=1)
QFF5_8 = KQUAD(name="QFF5_8", len=LQ50, k1=KFF5_8, rad=1)
QFF6_8 = KQUAD(name="QFF6_8", len=LQ50, k1=KFF6_8, rad=1)
q0ef_8 = KQUAD(name="q0ef_8", len=lq0ef_6, k1=K0EF_8, rad=1)
q1ef_8 = KQUAD(name="q1ef_8", len=lq1ef_6, k1=K1EF_8, rad=1)
q2ef_8 = KQUAD(name="q2ef_8", len=0.6, k1=K2EF_8, rad=1)
q3ef_8 = KQUAD(name="q3ef_8", len=0.6, k1=K3EF_8, rad=1)
q4ef_8 = KQUAD(name="q4ef_8", len=1.2, k1=K4EF_8, rad=1)
q5ef_8 = KQUAD(name="q5ef_8", len=1.2, k1=K5EF_8, rad=1)
q6ef_8 = KQUAD(name="q6ef_8", len=1.2, k1=K6EF_8, rad=1)
q7ef_8 = KQUAD(name="q7ef_8", len=1.2, k1=K7EF_8, rad=1)
q8ef_8 = KQUAD(name="q8ef_8", len=1.2, k1=K8EF_8, rad=1)
q9ef_8 = KQUAD(name="q9ef_8", len=1.2, k1=K9EF_8, rad=1)
q10ef_8 = KQUAD(name="q10ef_8", len=1.2, k1=K10EF_8, rad=1)
q11ef_8 = KQUAD(name="q11ef_8", len=1.2, k1=K11EF_8, rad=1)
q12ef_8 = KQUAD(name="q12ef_8", len=1.2, k1=K12EF_8, rad=1)
q13ef_8 = KQUAD(name="q13ef_8", len=1.2, k1=K13EF_8, rad=1)
q14ef_8 = KQUAD(name="q14ef_8", len=1.2, k1=K14EF_8, rad=1)
q15ef_8 = KQUAD(name="q15ef_8", len=1.2, k1=K15EF_8, rad=1)
q1er_8 = KQUAD(name="q1er_8", len=lq1er_6, k1=K1ER_8, rad=1)
q2er_8 = KQUAD(name="q2er_8", len=lq2er_6, k1=K2ER_8, rad=1)
q3er_8 = KQUAD(name="q3er_8", len=0.6, k1=K3ER_8, rad=1)
q4er_8 = KQUAD(name="q4er_8", len=0.6, k1=K4ER_8, rad=1)
q5er_8 = KQUAD(name="q5er_8", len=1.2, k1=K5ER_8, rad=1)
q6er_8 = KQUAD(name="q6er_8", len=1.2, k1=K6ER_8, rad=1)
q7er_8 = KQUAD(name="q7er_8", len=1.2, k1=K7ER_8, rad=1)
q8er_8 = KQUAD(name="q8er_8", len=1.2, k1=K8ER_8, rad=1)
q9er_8 = KQUAD(name="q9er_8", len=1.2, k1=K9ER_8, rad=1)
q10er_8 = KQUAD(name="q10er_8", len=1.2, k1=K10ER_8, rad=1)
q11er_8 = KQUAD(name="q11er_8", len=1.2, k1=K11ER_8, rad=1)
q12er_8 = KQUAD(name="q12er_8", len=1.2, k1=K12ER_8, rad=1)
q13er_8 = KQUAD(name="q13er_8", len=1.2, k1=K13ER_8, rad=1)
q14er_8 = KQUAD(name="q14er_8", len=1.2, k1=K14ER_8, rad=1)
q15er_8 = KQUAD(name="q15er_8", len=1.2, k1=K15ER_8, rad=1)
#
# IP10 quadrupoles.
HQM22_9 = KQUAD(name="HQM22_9", len=LQ80/2.0, k1=KM22_9, rad=1)
HQM21_9 = KQUAD(name="HQM21_9", len=LQ80/2.0, k1=KM21_9, rad=1)
HQM20_9 = KQUAD(name="HQM20_9", len=LQ80/2.0, k1=KM20_9, rad=1)
HQM19_9 = KQUAD(name="HQM19_9", len=LQ80/2.0, k1=KM19_9, rad=1)
HQM18_9 = KQUAD(name="HQM18_9", len=LQ80/2.0, k1=KM18_9, rad=1)
HQM17_9 = KQUAD(name="HQM17_9", len=LQ80/2.0, k1=KM17_9, rad=1)
HQM16_9 = KQUAD(name="HQM16_9", len=LQ80/2.0, k1=KM16_9, rad=1)
HQM15_9 = KQUAD(name="HQM15_9", len=LQ80/2.0, k1=KM15_9, rad=1)
HQM14_9 = KQUAD(name="HQM14_9", len=LQ80/2.0, k1=KM14_9, rad=1)
HQM13_9 = KQUAD(name="HQM13_9", len=LQ80/2.0, k1=KM13_9, rad=1)
HQM12_9 = KQUAD(name="HQM12_9", len=LQ80/2.0, k1=KM12_9, rad=1)
HQDSS_10 = KQUAD(name="HQDSS_10", len=LQ80/2.0, k1=KDSS_10, rad=1)
HQFSS_10 = KQUAD(name="HQFSS_10", len=LQ80/2.0, k1=KFSS_10, rad=1)
HQDLSS_10 = KQUAD(name="HQDLSS_10", len=LQ60, k1=KDSSL_10, rad=1)
HQFLSS_10 = KQUAD(name="HQFLSS_10", len=LQ60, k1=KFSSL_10, rad=1)
HQM22_10 = KQUAD(name="HQM22_10", len=LQ80/2.0, k1=KM22_10, rad=1)
HQM21_10 = KQUAD(name="HQM21_10", len=LQ80/2.0, k1=KM21_10, rad=1)
HQM20_10 = KQUAD(name="HQM20_10", len=LQ80/2.0, k1=KM20_10, rad=1)
HQM19_10 = KQUAD(name="HQM19_10", len=LQ80/2.0, k1=KM19_10, rad=1)
HQM18_10 = KQUAD(name="HQM18_10", len=LQ80/2.0, k1=KM18_10, rad=1)
HQM17_10 = KQUAD(name="HQM17_10", len=LQ80/2.0, k1=KM17_10, rad=1)
HQM16_10 = KQUAD(name="HQM16_10", len=LQ80/2.0, k1=KM16_10, rad=1)
HQM15_10 = KQUAD(name="HQM15_10", len=LQ80/2.0, k1=KM15_10, rad=1)
HQM14_10 = KQUAD(name="HQM14_10", len=LQ80/2.0, k1=KM14_10, rad=1)
HQM13_10 = KQUAD(name="HQM13_10", len=LQ80/2.0, k1=KM13_10, rad=1)
HQM12_10 = KQUAD(name="HQM12_10", len=LQ80/2.0, k1=KM12_10, rad=1)
#
# IP12 quadrupoles.
HQM22_11 = KQUAD(name="HQM22_11", len=LQ80/2.0, k1=KM22_11, rad=1)
HQM21_11 = KQUAD(name="HQM21_11", len=LQ80/2.0, k1=KM21_11, rad=1)
HQM20_11 = KQUAD(name="HQM20_11", len=LQ80/2.0, k1=KM20_11, rad=1)
HQM19_11 = KQUAD(name="HQM19_11", len=LQ80/2.0, k1=KM19_11, rad=1)
HQM18_11 = KQUAD(name="HQM18_11", len=LQ80/2.0, k1=KM18_11, rad=1)
HQM17_11 = KQUAD(name="HQM17_11", len=LQ80/2.0, k1=KM17_11, rad=1)
HQM16_11 = KQUAD(name="HQM16_11", len=LQ80/2.0, k1=KM16_11, rad=1)
HQM15_11 = KQUAD(name="HQM15_11", len=LQ80/2.0, k1=KM15_11, rad=1)
HQM14_11 = KQUAD(name="HQM14_11", len=LQ80/2.0, k1=KM14_11, rad=1)
HQM13_11 = KQUAD(name="HQM13_11", len=LQ80/2.0, k1=KM13_11, rad=1)
HQM12_11 = KQUAD(name="HQM12_11", len=LQ80/2.0, k1=KM12_11, rad=1)
HQDSS_12 = KQUAD(name="HQDSS_12", len=LQ80/2.0, k1=KDSS_12, rad=1)
HQFSS_12 = KQUAD(name="HQFSS_12", len=LQ80/2.0, k1=KFSS_12, rad=1)
HQM22_12 = KQUAD(name="HQM22_12", len=LQ80/2.0, k1=KM22_12, rad=1)
HQM21_12 = KQUAD(name="HQM21_12", len=LQ80/2.0, k1=KM21_12, rad=1)
HQM20_12 = KQUAD(name="HQM20_12", len=LQ80/2.0, k1=KM20_12, rad=1)
HQM19_12 = KQUAD(name="HQM19_12", len=LQ80/2.0, k1=KM19_12, rad=1)
HQM18_12 = KQUAD(name="HQM18_12", len=LQ80/2.0, k1=KM18_12, rad=1)
HQM17_12 = KQUAD(name="HQM17_12", len=LQ80/2.0, k1=KM17_12, rad=1)
HQM16_12 = KQUAD(name="HQM16_12", len=LQ80/2.0, k1=KM16_12, rad=1)
HQM15_12 = KQUAD(name="HQM15_12", len=LQ80/2.0, k1=KM15_12, rad=1)
HQM14_12 = KQUAD(name="HQM14_12", len=LQ80/2.0, k1=KM14_12, rad=1)
#
# IP2 quadrupoles.
HQM22_1 = KQUAD(name="HQM22_1", len=LQ80/2.0, k1=KM22_1)
HQM21_1 = KQUAD(name="HQM21_1", len=LQ80/2.0, k1=KM21_1)
HQM20_1 = KQUAD(name="HQM20_1", len=LQ80/2.0, k1=KM20_1)
HQM19_1 = KQUAD(name="HQM19_1", len=LQ80/2.0, k1=KM19_1)
HQM18_1 = KQUAD(name="HQM18_1", len=LQ80/2.0, k1=KM18_1)
HQM17_1 = KQUAD(name="HQM17_1", len=LQ80/2.0, k1=KM17_1)
HQM16_1 = KQUAD(name="HQM16_1", len=LQ80/2.0, k1=KM16_1)
HQM15_1 = KQUAD(name="HQM15_1", len=LQ80/2.0, k1=KM15_1)
HQM14_1 = KQUAD(name="HQM14_1", len=LQ80/2.0, k1=KM14_1)
HQM13_1 = KQUAD(name="HQM13_1", len=LQ80/2.0, k1=KM13_1)
HQM12_1 = KQUAD(name="HQM12_1", len=LQ80/2.0, k1=KM12_1)
HQDSS_2 = KQUAD(name="HQDSS_2", len=LQ80/2.0, k1=KDSS_2)
HQFSS_2 = KQUAD(name="HQFSS_2", len=LQ80/2.0, k1=KFSS_2)
HQM22_2 = KQUAD(name="HQM22_2", len=LQ80/2.0, k1=KM22_2)
HQM21_2 = KQUAD(name="HQM21_2", len=LQ80/2.0, k1=KM21_2)
HQM20_2 = KQUAD(name="HQM20_2", len=LQ80/2.0, k1=KM20_2)
HQM19_2 = KQUAD(name="HQM19_2", len=LQ80/2.0, k1=KM19_2)
HQM18_2 = KQUAD(name="HQM18_2", len=LQ80/2.0, k1=KM18_2)
HQM17_2 = KQUAD(name="HQM17_2", len=LQ80/2.0, k1=KM17_2)
HQM16_2 = KQUAD(name="HQM16_2", len=LQ80/2.0, k1=KM16_2)
HQM15_2 = KQUAD(name="HQM15_2", len=LQ80/2.0, k1=KM15_2)
HQM14_2 = KQUAD(name="HQM14_2", len=LQ80/2.0, k1=KM14_2)
HQM13_2 = KQUAD(name="HQM13_2", len=LQ80/2.0, k1=KM13_2)
HQM12_2 = KQUAD(name="HQM12_2", len=LQ80/2.0, k1=KM12_2)
#
# IP4 quadrupoles.
HQD17_3 = KQUAD(name="HQD17_3", len=LQ80/2.0, k1=KD17_3, rad=1)
HQF16_3 = KQUAD(name="HQF16_3", len=LQ80/2.0, k1=KF16_3, rad=1)
HQD15_3 = KQUAD(name="HQD15_3", len=LQ80/2.0, k1=KD15_3, rad=1)
HQF14_3 = KQUAD(name="HQF14_3", len=LQ80/2.0, k1=KF14_3, rad=1)
HQD13_3 = KQUAD(name="HQD13_3", len=LQ80/2.0, k1=KD13_3, rad=1)
HQF12_3 = KQUAD(name="HQF12_3", len=LQ80/2.0, k1=KF12_3, rad=1)
HQD11_3 = KQUAD(name="HQD11_3", len=LQ80/2.0, k1=KD11_3, rad=1)
HQF10_3 = KQUAD(name="HQF10_3", len=LQ80/2.0, k1=KF10_3, rad=1)
HQD9_3 = KQUAD(name="HQD9_3", len=LQ80/2.0, k1=KD9_3, rad=1)
HQF8_3 = KQUAD(name="HQF8_3", len=LQ80/2.0, k1=KF8_3, rad=1)
HQD7_3 = KQUAD(name="HQD7_3", len=LQ80/2.0, k1=KD7_3, rad=1)
HQF6_3 = KQUAD(name="HQF6_3", len=LQ80/2.0, k1=KF6_3, rad=1)
HQD5_3 = KQUAD(name="HQD5_3", len=LQ80/2.0, k1=KD5_3, rad=1)
HQF4_3 = KQUAD(name="HQF4_3", len=LQ80/2.0, k1=KF4_3, rad=1)
HQD3_3 = KQUAD(name="HQD3_3", len=LQ80/2.0, k1=KD3_3, rad=1)
HQF2_3 = KQUAD(name="HQF2_3", len=LQ80/2.0, k1=KF2_3, rad=1)
HQD1_3 = KQUAD(name="HQD1_3", len=LQ80/2.0, k1=KD1_3, rad=1)
HQF18_4 = KQUAD(name="HQF18_4", len=LQ80/2.0, k1=KF18_4, rad=1)
HQD17_4 = KQUAD(name="HQD17_4", len=LQ80/2.0, k1=KD17_4, rad=1)
HQF16_4 = KQUAD(name="HQF16_4", len=LQ80/2.0, k1=KF16_4, rad=1)
HQD15_4 = KQUAD(name="HQD15_4", len=LQ80/2.0, k1=KD15_4, rad=1)
HQF14_4 = KQUAD(name="HQF14_4", len=LQ80/2.0, k1=KF14_4, rad=1)
HQD13_4 = KQUAD(name="HQD13_4", len=LQ80/2.0, k1=KD13_4, rad=1)
HQF12_4 = KQUAD(name="HQF12_4", len=LQ80/2.0, k1=KF12_4, rad=1)
HQD11_4 = KQUAD(name="HQD11_4", len=LQ80/2.0, k1=KD11_4, rad=1)
HQF10_4 = KQUAD(name="HQF10_4", len=LQ80/2.0, k1=KF10_4, rad=1)
HQD9_4 = KQUAD(name="HQD9_4", len=LQ80/2.0, k1=KD9_4, rad=1)
HQF8_4 = KQUAD(name="HQF8_4", len=LQ80/2.0, k1=KF8_4, rad=1)
HQD7_4 = KQUAD(name="HQD7_4", len=LQ80/2.0, k1=KD7_4, rad=1)
HQF6_4 = KQUAD(name="HQF6_4", len=LQ80/2.0, k1=KF6_4, rad=1)
HQD5_4 = KQUAD(name="HQD5_4", len=LQ80/2.0, k1=KD5_4, rad=1)
HQF4_4 = KQUAD(name="HQF4_4", len=LQ80/2.0, k1=KF4_4, rad=1)
HQD3_4 = KQUAD(name="HQD3_4", len=LQ80/2.0, k1=KD3_4, rad=1)
HQF2_4 = KQUAD(name="HQF2_4", len=LQ80/2.0, k1=KF2_4, rad=1)
HQD1_4 = KQUAD(name="HQD1_4", len=LQ80/2.0, k1=KD1_4, rad=1)
#
# ======= The sextupoles ==================================
# ARC 1 sextupoles.
SFM1_1 = KSEXT(name="SFM1_1", len=LSX, k2=KSFM1_1, rad=1)
SF00_1 = KSEXT(name="SF00_1", len=LSX, k2=KSF00_1, rad=1)
SD00_1 = KSEXT(name="SD00_1", len=LSX, k2=KSD00_1, rad=1)
SF01_1 = KSEXT(name="SF01_1", len=LSX, k2=KSF01_1, rad=1)
SD01_1 = KSEXT(name="SD01_1", len=LSX, k2=KSD01_1, rad=1)
SF02_1 = KSEXT(name="SF02_1", len=LSX, k2=KSF02_1, rad=1)
SD02_1 = KSEXT(name="SD02_1", len=LSX, k2=KSD02_1, rad=1)
SF03_1 = KSEXT(name="SF03_1", len=LSX, k2=KSF03_1, rad=1)
SD03_1 = KSEXT(name="SD03_1", len=LSX, k2=KSD03_1, rad=1)
SF04_1 = KSEXT(name="SF04_1", len=LSX, k2=KSF04_1, rad=1)
SD04_1 = KSEXT(name="SD04_1", len=LSX, k2=KSD04_1, rad=1)
SF05_1 = KSEXT(name="SF05_1", len=LSX, k2=KSF05_1, rad=1)
SD05_1 = KSEXT(name="SD05_1", len=LSX, k2=KSD05_1, rad=1)
SF06_1 = KSEXT(name="SF06_1", len=LSX, k2=KSF06_1, rad=1)
SD06_1 = KSEXT(name="SD06_1", len=LSX, k2=KSD06_1, rad=1)
SF07_1 = KSEXT(name="SF07_1", len=LSX, k2=KSF07_1, rad=1)
SD07_1 = KSEXT(name="SD07_1", len=LSX, k2=KSD07_1, rad=1)
SF08_1 = KSEXT(name="SF08_1", len=LSX, k2=KSF08_1, rad=1)
SD08_1 = KSEXT(name="SD08_1", len=LSX, k2=KSD08_1, rad=1)
SF09_1 = KSEXT(name="SF09_1", len=LSX, k2=KSF09_1, rad=1)
SD09_1 = KSEXT(name="SD09_1", len=LSX, k2=KSD09_1, rad=1)
SF10_1 = KSEXT(name="SF10_1", len=LSX, k2=KSF10_1, rad=1)
SD10_1 = KSEXT(name="SD10_1", len=LSX, k2=KSD10_1, rad=1)
SF11_1 = KSEXT(name="SF11_1", len=LSX, k2=KSF11_1, rad=1)
SD11_1 = KSEXT(name="SD11_1", len=LSX, k2=KSD11_1, rad=1)
SF12_1 = KSEXT(name="SF12_1", len=LSX, k2=KSF12_1, rad=1)
SD12_1 = KSEXT(name="SD12_1", len=LSX, k2=KSD12_1, rad=1)
SF13_1 = KSEXT(name="SF13_1", len=LSX, k2=KSF13_1, rad=1)
SD13_1 = KSEXT(name="SD13_1", len=LSX, k2=KSD13_1, rad=1)
SF14_1 = KSEXT(name="SF14_1", len=LSX, k2=KSF14_1, rad=1)
SD14_1 = KSEXT(name="SD14_1", len=LSX, k2=KSD14_1, rad=1)
SF15_1 = KSEXT(name="SF15_1", len=LSX, k2=KSF15_1, rad=1)
SD15_1 = KSEXT(name="SD15_1", len=LSX, k2=KSD15_1, rad=1)
SF16_1 = KSEXT(name="SF16_1", len=LSX, k2=KSF16_1, rad=1)
SD16_1 = KSEXT(name="SD16_1", len=LSX, k2=KSD16_1, rad=1)
SF17_1 = KSEXT(name="SF17_1", len=LSX, k2=KSF17_1, rad=1)
#
# ARC 3 sextupoles.
SF00_3 = KSEXT(name="SF00_3", len=LSX, k2=KSF00_3, rad=1)
SD00_3 = KSEXT(name="SD00_3", len=LSX, k2=KSD00_3, rad=1)
SF01_3 = KSEXT(name="SF01_3", len=LSX, k2=KSF01_3, rad=1)
SD01_3 = KSEXT(name="SD01_3", len=LSX, k2=KSD01_3, rad=1)
SF02_3 = KSEXT(name="SF02_3", len=LSX, k2=KSF02_3, rad=1)
SD02_3 = KSEXT(name="SD02_3", len=LSX, k2=KSD02_3, rad=1)
SF03_3 = KSEXT(name="SF03_3", len=LSX, k2=KSF03_3, rad=1)
SD03_3 = KSEXT(name="SD03_3", len=LSX, k2=KSD03_3, rad=1)
SF04_3 = KSEXT(name="SF04_3", len=LSX, k2=KSF04_3, rad=1)
SD04_3 = KSEXT(name="SD04_3", len=LSX, k2=KSD04_3, rad=1)
SF05_3 = KSEXT(name="SF05_3", len=LSX, k2=KSF05_3, rad=1)
SD05_3 = KSEXT(name="SD05_3", len=LSX, k2=KSD05_3, rad=1)
SF06_3 = KSEXT(name="SF06_3", len=LSX, k2=KSF06_3, rad=1)
SD06_3 = KSEXT(name="SD06_3", len=LSX, k2=KSD06_3, rad=1)
SF07_3 = KSEXT(name="SF07_3", len=LSX, k2=KSF07_3, rad=1)
SD07_3 = KSEXT(name="SD07_3", len=LSX, k2=KSD07_3, rad=1)
SF08_3 = KSEXT(name="SF08_3", len=LSX, k2=KSF08_3, rad=1)
SD08_3 = KSEXT(name="SD08_3", len=LSX, k2=KSD08_3, rad=1)
SF09_3 = KSEXT(name="SF09_3", len=LSX, k2=KSF09_3, rad=1)
SD09_3 = KSEXT(name="SD09_3", len=LSX, k2=KSD09_3, rad=1)
SF10_3 = KSEXT(name="SF10_3", len=LSX, k2=KSF10_3, rad=1)
SD10_3 = KSEXT(name="SD10_3", len=LSX, k2=KSD10_3, rad=1)
SF11_3 = KSEXT(name="SF11_3", len=LSX, k2=KSF11_3, rad=1)
SD11_3 = KSEXT(name="SD11_3", len=LSX, k2=KSD11_3, rad=1)
SF12_3 = KSEXT(name="SF12_3", len=LSX, k2=KSF12_3, rad=1)
SD12_3 = KSEXT(name="SD12_3", len=LSX, k2=KSD12_3, rad=1)
SF13_3 = KSEXT(name="SF13_3", len=LSX, k2=KSF13_3, rad=1)
SD13_3 = KSEXT(name="SD13_3", len=LSX, k2=KSD13_3, rad=1)
SF14_3 = KSEXT(name="SF14_3", len=LSX, k2=KSF14_3, rad=1)
SD14_3 = KSEXT(name="SD14_3", len=LSX, k2=KSD14_3, rad=1)
SF15_3 = KSEXT(name="SF15_3", len=LSX, k2=KSF15_3, rad=1)
SD15_3 = KSEXT(name="SD15_3", len=LSX, k2=KSD15_3, rad=1)
SF16_3 = KSEXT(name="SF16_3", len=LSX, k2=KSF16_3, rad=1)
SD16_3 = KSEXT(name="SD16_3", len=LSX, k2=KSD16_3, rad=1)
SF17_3 = KSEXT(name="SF17_3", len=LSX, k2=KSF17_3, rad=1)
SF18_3 = KSEXT(name="SF18_3", len=LSX, k2=KSF18_3, rad=1)
# ARC 5 sextupoles
SFM1_5 = KSEXT(name="SFM1_5", len=LSXL, k2=KSFM1_5, rad=1)
SF00_5 = KSEXT(name="SF00_5", len=LSXL, k2=KSF00_5, rad=1)
SD00_5 = KSEXT(name="SD00_5", len=LSXL, k2=KSD00_5, rad=1)
SF01_5 = KSEXT(name="SF01_5", len=LSX, k2=KSF01_5, rad=1)
SD01_5 = KSEXT(name="SD01_5", len=LSXL, k2=KSD01_5, rad=1)
SF02_5 = KSEXT(name="SF02_5", len=LSX, k2=KSF02_5, rad=1)
SD02_5 = KSEXT(name="SD02_5", len=LSXL, k2=KSD02_5, rad=1)
SF03_5 = KSEXT(name="SF03_5", len=LSX, k2=KSF03_5, rad=1)
SD03_5 = KSEXT(name="SD03_5", len=LSXL, k2=KSD03_5, rad=1)
SF04_5 = KSEXT(name="SF04_5", len=LSX, k2=KSF04_5, rad=1)
SD04_5 = KSEXT(name="SD04_5", len=LSXL, k2=KSD04_5, rad=1)
SF05_5 = KSEXT(name="SF05_5", len=LSX, k2=KSF05_5, rad=1)
SD05_5 = KSEXT(name="SD05_5", len=LSXL, k2=KSD05_5, rad=1)
SF06_5 = KSEXT(name="SF06_5", len=LSX, k2=KSF06_5, rad=1)
SD06_5 = KSEXT(name="SD06_5", len=LSXL, k2=KSD06_5, rad=1)
SF07_5 = KSEXT(name="SF07_5", len=LSX, k2=KSF07_5, rad=1)
SD07_5 = KSEXT(name="SD07_5", len=LSXL, k2=KSD07_5, rad=1)
SF08_5 = KSEXT(name="SF08_5", len=LSX, k2=KSF08_5, rad=1)
SD08_5 = KSEXT(name="SD08_5", len=LSXL, k2=KSD08_5, rad=1)
SF09_5 = KSEXT(name="SF09_5", len=LSX, k2=KSF09_5, rad=1)
SD09_5 = KSEXT(name="SD09_5", len=LSXL, k2=KSD09_5, rad=1)
SF10_5 = KSEXT(name="SF10_5", len=LSX, k2=KSF10_5, rad=1)
SD10_5 = KSEXT(name="SD10_5", len=LSXL, k2=KSD10_5, rad=1)
SF11_5 = KSEXT(name="SF11_5", len=LSX, k2=KSF11_5, rad=1)
SD11_5 = KSEXT(name="SD11_5", len=LSXL, k2=KSD11_5, rad=1)
SF12_5 = KSEXT(name="SF12_5", len=LSX, k2=KSF12_5, rad=1)
SD12_5 = KSEXT(name="SD12_5", len=LSXL, k2=KSD12_5, rad=1)
SF13_5 = KSEXT(name="SF13_5", len=LSX, k2=KSF13_5, rad=1)
SD13_5 = KSEXT(name="SD13_5", len=LSXL, k2=KSD13_5, rad=1)
SF14_5 = KSEXT(name="SF14_5", len=LSX, k2=KSF14_5, rad=1)
SD14_5 = KSEXT(name="SD14_5", len=LSXL, k2=KSD14_5, rad=1)
SF15_5 = KSEXT(name="SF15_5", len=LSX, k2=KSF15_5, rad=1)
SD15_5 = KSEXT(name="SD15_5", len=LSXL, k2=KSD15_5, rad=1)
SF16_5 = KSEXT(name="SF16_5", len=LSX, k2=KSF16_5, rad=1)
SD16_5 = KSEXT(name="SD16_5", len=LSXL, k2=KSD16_5, rad=1)
# ARC 7 sextupoles
SF01_7 = KSEXT(name="SF01_7", len=LSX, k2=KSF01_7, rad=1)
SD01_7 = KSEXT(name="SD01_7", len=LSXL, k2=KSD01_7, rad=1)
SF02_7 = KSEXT(name="SF02_7", len=LSX, k2=KSF02_7, rad=1)
SD02_7 = KSEXT(name="SD02_7", len=LSXL, k2=KSD02_7, rad=1)
SF03_7 = KSEXT(name="SF03_7", len=LSX, k2=KSF03_7, rad=1)
SD03_7 = KSEXT(name="SD03_7", len=LSXL, k2=KSD03_7, rad=1)
SF04_7 = KSEXT(name="SF04_7", len=LSX, k2=KSF04_7, rad=1)
SD04_7 = KSEXT(name="SD04_7", len=LSXL, k2=KSD04_7, rad=1)
SF05_7 = KSEXT(name="SF05_7", len=LSX, k2=KSF05_7, rad=1)
SD05_7 = KSEXT(name="SD05_7", len=LSXL, k2=KSD05_7, rad=1)
SF06_7 = KSEXT(name="SF06_7", len=LSX, k2=KSF06_7, rad=1)
SD06_7 = KSEXT(name="SD06_7", len=LSXL, k2=KSD06_7, rad=1)
SF07_7 = KSEXT(name="SF07_7", len=LSX, k2=KSF07_7, rad=1)
SD07_7 = KSEXT(name="SD07_7", len=LSXL, k2=KSD07_7, rad=1)
SF08_7 = KSEXT(name="SF08_7", len=LSX, k2=KSF08_7, rad=1)
SD08_7 = KSEXT(name="SD08_7", len=LSXL, k2=KSD08_7, rad=1)
SF09_7 = KSEXT(name="SF09_7", len=LSX, k2=KSF09_7, rad=1)
SD09_7 = KSEXT(name="SD09_7", len=LSXL, k2=KSD09_7, rad=1)
SF10_7 = KSEXT(name="SF10_7", len=LSX, k2=KSF10_7, rad=1)
SD10_7 = KSEXT(name="SD10_7", len=LSXL, k2=KSD10_7, rad=1)
SF11_7 = KSEXT(name="SF11_7", len=LSX, k2=KSF11_7, rad=1)
SD11_7 = KSEXT(name="SD11_7", len=LSXL, k2=KSD11_7, rad=1)
SF12_7 = KSEXT(name="SF12_7", len=LSX, k2=KSF12_7, rad=1)
SD12_7 = KSEXT(name="SD12_7", len=LSXL, k2=KSD12_7, rad=1)
SF13_7 = KSEXT(name="SF13_7", len=LSX, k2=KSF13_7, rad=1)
SD13_7 = KSEXT(name="SD13_7", len=LSXL, k2=KSD13_7, rad=1)
SF14_7 = KSEXT(name="SF14_7", len=LSX, k2=KSF14_7, rad=1)
SD14_7 = KSEXT(name="SD14_7", len=LSXL, k2=KSD14_7, rad=1)
SF15_7 = KSEXT(name="SF15_7", len=LSX, k2=KSF15_7, rad=1)
SD15_7 = KSEXT(name="SD15_7", len=LSXL, k2=KSD15_7, rad=1)
SF16_7 = KSEXT(name="SF16_7", len=LSX, k2=KSF16_7, rad=1)
SD16_7 = KSEXT(name="SD16_7", len=LSXL, k2=KSD16_7, rad=1)
# ARC 9 sextupoles
SF01_9 = KSEXT(name="SF01_9", len=LSX, k2=KSF01_9, rad=1)
SD01_9 = KSEXT(name="SD01_9", len=LSXL, k2=KSD01_9, rad=1)
SF02_9 = KSEXT(name="SF02_9", len=LSX, k2=KSF02_9, rad=1)
SD02_9 = KSEXT(name="SD02_9", len=LSXL, k2=KSD02_9, rad=1)
SF03_9 = KSEXT(name="SF03_9", len=LSX, k2=KSF03_9, rad=1)
SD03_9 = KSEXT(name="SD03_9", len=LSXL, k2=KSD03_9, rad=1)
SF04_9 = KSEXT(name="SF04_9", len=LSX, k2=KSF04_9, rad=1)
SD04_9 = KSEXT(name="SD04_9", len=LSXL, k2=KSD04_9, rad=1)
SF05_9 = KSEXT(name="SF05_9", len=LSX, k2=KSF05_9, rad=1)
SD05_9 = KSEXT(name="SD05_9", len=LSXL, k2=KSD05_9, rad=1)
SF06_9 = KSEXT(name="SF06_9", len=LSX, k2=KSF06_9, rad=1)
SD06_9 = KSEXT(name="SD06_9", len=LSXL, k2=KSD06_9, rad=1)
SF07_9 = KSEXT(name="SF07_9", len=LSX, k2=KSF07_9, rad=1)
SD07_9 = KSEXT(name="SD07_9", len=LSXL, k2=KSD07_9, rad=1)
SF08_9 = KSEXT(name="SF08_9", len=LSX, k2=KSF08_9, rad=1)
SD08_9 = KSEXT(name="SD08_9", len=LSXL, k2=KSD08_9, rad=1)
SF09_9 = KSEXT(name="SF09_9", len=LSX, k2=KSF09_9, rad=1)
SD09_9 = KSEXT(name="SD09_9", len=LSXL, k2=KSD09_9, rad=1)
SF10_9 = KSEXT(name="SF10_9", len=LSX, k2=KSF10_9, rad=1)
SD10_9 = KSEXT(name="SD10_9", len=LSXL, k2=KSD10_9, rad=1)
SF11_9 = KSEXT(name="SF11_9", len=LSX, k2=KSF11_9, rad=1)
SD11_9 = KSEXT(name="SD11_9", len=LSXL, k2=KSD11_9, rad=1)
SF12_9 = KSEXT(name="SF12_9", len=LSX, k2=KSF12_9, rad=1)
SD12_9 = KSEXT(name="SD12_9", len=LSXL, k2=KSD12_9, rad=1)
SF13_9 = KSEXT(name="SF13_9", len=LSX, k2=KSF13_9, rad=1)
SD13_9 = KSEXT(name="SD13_9", len=LSXL, k2=KSD13_9, rad=1)
SF14_9 = KSEXT(name="SF14_9", len=LSX, k2=KSF14_9, rad=1)
SD14_9 = KSEXT(name="SD14_9", len=LSXL, k2=KSD14_9, rad=1)
SF15_9 = KSEXT(name="SF15_9", len=LSX, k2=KSF15_9, rad=1)
SD15_9 = KSEXT(name="SD15_9", len=LSXL, k2=KSD15_9, rad=1)
SF16_9 = KSEXT(name="SF16_9", len=LSX, k2=KSF16_9, rad=1)
SD16_9 = KSEXT(name="SD16_9", len=LSXL, k2=KSD16_9, rad=1)
SF17_9 = KSEXT(name="SF17_9", len=LSX, k2=KSF17_9, rad=1)
# ARC 11 sextupoles
SF00_11 = KSEXT(name="SF00_11", len=LSX, k2=KSF00_11, rad=1)
SD00_11 = KSEXT(name="SD00_11", len=LSX, k2=KSD00_11, rad=1)
SF01_11 = KSEXT(name="SF01_11", len=LSX, k2=KSF01_11, rad=1)
SD01_11 = KSEXT(name="SD01_11", len=LSX, k2=KSD01_11, rad=1)
SF02_11 = KSEXT(name="SF02_11", len=LSX, k2=KSF02_11, rad=1)
SD02_11 = KSEXT(name="SD02_11", len=LSX, k2=KSD02_11, rad=1)
SF03_11 = KSEXT(name="SF03_11", len=LSX, k2=KSF03_11, rad=1)
SD03_11 = KSEXT(name="SD03_11", len=LSX, k2=KSD03_11, rad=1)
SF04_11 = KSEXT(name="SF04_11", len=LSX, k2=KSF04_11, rad=1)
SD04_11 = KSEXT(name="SD04_11", len=LSX, k2=KSD04_11, rad=1)
SF05_11 = KSEXT(name="SF05_11", len=LSX, k2=KSF05_11, rad=1)
SD05_11 = KSEXT(name="SD05_11", len=LSX, k2=KSD05_11, rad=1)
SF06_11 = KSEXT(name="SF06_11", len=LSX, k2=KSF06_11, rad=1)
SD06_11 = KSEXT(name="SD06_11", len=LSX, k2=KSD06_11, rad=1)
SF07_11 = KSEXT(name="SF07_11", len=LSX, k2=KSF07_11, rad=1)
SD07_11 = KSEXT(name="SD07_11", len=LSX, k2=KSD07_11, rad=1)
SF08_11 = KSEXT(name="SF08_11", len=LSX, k2=KSF08_11, rad=1)
SD08_11 = KSEXT(name="SD08_11", len=LSX, k2=KSD08_11, rad=1)
SF09_11 = KSEXT(name="SF09_11", len=LSX, k2=KSF09_11, rad=1)
SD09_11 = KSEXT(name="SD09_11", len=LSX, k2=KSD09_11, rad=1)
SF10_11 = KSEXT(name="SF10_11", len=LSX, k2=KSF10_11, rad=1)
SD10_11 = KSEXT(name="SD10_11", len=LSX, k2=KSD10_11, rad=1)
SF11_11 = KSEXT(name="SF11_11", len=LSX, k2=KSF11_11, rad=1)
SD11_11 = KSEXT(name="SD11_11", len=LSX, k2=KSD11_11, rad=1)
SF12_11 = KSEXT(name="SF12_11", len=LSX, k2=KSF12_11, rad=1)
SD12_11 = KSEXT(name="SD12_11", len=LSX, k2=KSD12_11, rad=1)
SF13_11 = KSEXT(name="SF13_11", len=LSX, k2=KSF13_11, rad=1)
SD13_11 = KSEXT(name="SD13_11", len=LSX, k2=KSD13_11, rad=1)
SF14_11 = KSEXT(name="SF14_11", len=LSX, k2=KSF14_11, rad=1)
SD14_11 = KSEXT(name="SD14_11", len=LSX, k2=KSD14_11, rad=1)
SF15_11 = KSEXT(name="SF15_11", len=LSX, k2=KSF15_11, rad=1)
SD15_11 = KSEXT(name="SD15_11", len=LSX, k2=KSD15_11, rad=1)
SF16_11 = KSEXT(name="SF16_11", len=LSX, k2=KSF16_11, rad=1)
SD16_11 = KSEXT(name="SD16_11", len=LSX, k2=KSD16_11, rad=1)
SF17_11 = KSEXT(name="SF17_11", len=LSX, k2=KSF17_11, rad=1)
SF18_11 = KSEXT(name="SF18_11", len=LSX, k2=KSF18_11, rad=1)
# Straight section
SX41_2 = KSEXT(name="SX41_2", len=LSX, k2=S41_2, rad=1)
SX42_2 = KSEXT(name="SX42_2", len=LSX, k2=S42_2, rad=1)
SX43_2 = KSEXT(name="SX43_2", len=LSX, k2=S43_2, rad=1)
SX44_2 = KSEXT(name="SX44_2", len=LSX, k2=S44_2, rad=1)
SX45_2 = KSEXT(name="SX45_2", len=LSX, k2=S45_2, rad=1)
SX46_2 = KSEXT(name="SX46_2", len=LSX, k2=S46_2, rad=1)
SX47_2 = KSEXT(name="SX47_2", len=LSX, k2=S47_2, rad=1)
SX48_2 = KSEXT(name="SX48_2", len=LSX, k2=S48_2, rad=1)
SX49_2 = KSEXT(name="SX49_2", len=LSX, k2=S49_2, rad=1)
SX50_2 = KSEXT(name="SX50_2", len=LSX, k2=S50_2, rad=1)
SX51_2 = KSEXT(name="SX51_2", len=LSX, k2=S51_2, rad=1)
SX52_2 = KSEXT(name="SX52_2", len=LSX, k2=S52_2, rad=1)
#
# ======= Correctors. ======================================
CH00_1 = HKICKER(name="CH00_1", len=LCH)
CV00_1 = VKICKER(name="CV00_1", len=LCV)
CH01_1 = HKICKER(name="CH01_1", len=LCH)
CV01_1 = VKICKER(name="CV01_1", len=LCV)
CH02_1 = HKICKER(name="CH02_1", len=LCH)
CV02_1 = VKICKER(name="CV02_1", len=LCV)
CH03_1 = HKICKER(name="CH03_1", len=LCH)
CV03_1 = VKICKER(name="CV03_1", len=LCV)
CH04_1 = HKICKER(name="CH04_1", len=LCH)
CV04_1 = VKICKER(name="CV04_1", len=LCV)
CH05_1 = HKICKER(name="CH05_1", len=LCH)
CV05_1 = VKICKER(name="CV05_1", len=LCV)
CH06_1 = HKICKER(name="CH06_1", len=LCH)
CV06_1 = VKICKER(name="CV06_1", len=LCV)
CH07_1 = HKICKER(name="CH07_1", len=LCH)
CV07_1 = VKICKER(name="CV07_1", len=LCV)
CH08_1 = HKICKER(name="CH08_1", len=LCH)
CV08_1 = VKICKER(name="CV08_1", len=LCV)
CH09_1 = HKICKER(name="CH09_1", len=LCH)
CV09_1 = VKICKER(name="CV09_1", len=LCV)
CH10_1 = HKICKER(name="CH10_1", len=LCH)
CV10_1 = VKICKER(name="CV10_1", len=LCV)
CH11_1 = HKICKER(name="CH11_1", len=LCH)
CV11_1 = VKICKER(name="CV11_1", len=LCV)
CH12_1 = HKICKER(name="CH12_1", len=LCH)
CV12_1 = VKICKER(name="CV12_1", len=LCV)
CH13_1 = HKICKER(name="CH13_1", len=LCH)
CV13_1 = VKICKER(name="CV13_1", len=LCV)
CH14_1 = HKICKER(name="CH14_1", len=LCH)
CV14_1 = VKICKER(name="CV14_1", len=LCV)
CH15_1 = HKICKER(name="CH15_1", len=LCH)
CV15_1 = VKICKER(name="CV15_1", len=LCV)
CH16_1 = HKICKER(name="CH16_1", len=LCH)
CV16_1 = VKICKER(name="CV16_1", len=LCV)
CH17_1 = HKICKER(name="CH17_1", len=LCH)
CV17_1 = VKICKER(name="CV17_1", len=LCV)
CH00_3 = HKICKER(name="CH00_3", len=LCH)
CV00_3 = HKICKER(name="CV00_3", len=LCV)
CH01_3 = HKICKER(name="CH01_3", len=LCH)
CV01_3 = VKICKER(name="CV01_3", len=LCV)
CH02_3 = HKICKER(name="CH02_3", len=LCH)
CV02_3 = VKICKER(name="CV02_3", len=LCV)
CH03_3 = HKICKER(name="CH03_3", len=LCH)
CV03_3 = VKICKER(name="CV03_3", len=LCV)
CH04_3 = HKICKER(name="CH04_3", len=LCH)
CV04_3 = VKICKER(name="CV04_3", len=LCV)
CH05_3 = HKICKER(name="CH05_3", len=LCH)
CV05_3 = VKICKER(name="CV05_3", len=LCV)
CH06_3 = HKICKER(name="CH06_3", len=LCH)
CV06_3 = VKICKER(name="CV06_3", len=LCV)
CH07_3 = HKICKER(name="CH07_3", len=LCH)
CV07_3 = VKICKER(name="CV07_3", len=LCV)
CH08_3 = HKICKER(name="CH08_3", len=LCH)
CV08_3 = VKICKER(name="CV08_3", len=LCV)
CH09_3 = HKICKER(name="CH09_3", len=LCH)
CV09_3 = VKICKER(name="CV09_3", len=LCV)
CH10_3 = HKICKER(name="CH10_3", len=LCH)
CV10_3 = VKICKER(name="CV10_3", len=LCV)
CH11_3 = HKICKER(name="CH11_3", len=LCH)
CV11_3 = VKICKER(name="CV11_3", len=LCV)
CH12_3 = HKICKER(name="CH12_3", len=LCH)
CV12_3 = VKICKER(name="CV12_3", len=LCV)
CH13_3 = HKICKER(name="CH13_3", len=LCH)
CV13_3 = VKICKER(name="CV13_3", len=LCV)
CH14_3 = HKICKER(name="CH14_3", len=LCH)
CV14_3 = VKICKER(name="CV14_3", len=LCV)
CH15_3 = HKICKER(name="CH15_3", len=LCH)
CV15_3 = VKICKER(name="CV15_3", len=LCV)
CH16_3 = HKICKER(name="CH16_3", len=LCH)
CV16_3 = VKICKER(name="CV16_3", len=LCV)
CH17_3 = HKICKER(name="CH17_3", len=LCH)
CV17_3 = VKICKER(name="CV17_3", len=LCV)
CH00_5 = HKICKER(name="CH00_5", len=LCH)
CV00_5 = VKICKER(name="CV00_5", len=LCV)
CH01_5 = HKICKER(name="CH01_5", len=LCH)
CV01_5 = VKICKER(name="CV01_5", len=LCV)
CH02_5 = HKICKER(name="CH02_5", len=LCH)
CV02_5 = VKICKER(name="CV02_5", len=LCV)
CH03_5 = HKICKER(name="CH03_5", len=LCH)
CV03_5 = VKICKER(name="CV03_5", len=LCV)
CH04_5 = HKICKER(name="CH04_5", len=LCH)
CV04_5 = VKICKER(name="CV04_5", len=LCV)
CH05_5 = HKICKER(name="CH05_5", len=LCH)
CV05_5 = VKICKER(name="CV05_5", len=LCV)
CH06_5 = HKICKER(name="CH06_5", len=LCH)
CV06_5 = VKICKER(name="CV06_5", len=LCV)
CH07_5 = HKICKER(name="CH07_5", len=LCH)
CV07_5 = VKICKER(name="CV07_5", len=LCV)
CH08_5 = HKICKER(name="CH08_5", len=LCH)
CV08_5 = VKICKER(name="CV08_5", len=LCV)
CH09_5 = HKICKER(name="CH09_5", len=LCH)
CV09_5 = VKICKER(name="CV09_5", len=LCV)
CH10_5 = HKICKER(name="CH10_5", len=LCH)
CV10_5 = VKICKER(name="CV10_5", len=LCV)
CH11_5 = HKICKER(name="CH11_5", len=LCH)
CV11_5 = VKICKER(name="CV11_5", len=LCV)
CH12_5 = HKICKER(name="CH12_5", len=LCH)
CV12_5 = VKICKER(name="CV12_5", len=LCV)
CH13_5 = HKICKER(name="CH13_5", len=LCH)
CV13_5 = VKICKER(name="CV13_5", len=LCV)
CH14_5 = HKICKER(name="CH14_5", len=LCH)
CV14_5 = VKICKER(name="CV14_5", len=LCV)
CH15_5 = HKICKER(name="CH15_5", len=LCH)
CV15_5 = VKICKER(name="CV15_5", len=LCV)
CH16_5 = HKICKER(name="CH16_5", len=LCH)
CV16_5 = VKICKER(name="CV16_5", len=LCV)
CH01_7 = HKICKER(name="CH01_7", len=LCH)
CV01_7 = VKICKER(name="CV01_7", len=LCV)
CH02_7 = HKICKER(name="CH02_7", len=LCH)
CV02_7 = VKICKER(name="CV02_7", len=LCV)
CH03_7 = HKICKER(name="CH03_7", len=LCH)
CV03_7 = VKICKER(name="CV03_7", len=LCV)
CH04_7 = HKICKER(name="CH04_7", len=LCH)
CV04_7 = VKICKER(name="CV04_7", len=LCV)
CH05_7 = HKICKER(name="CH05_7", len=LCH)
CV05_7 = VKICKER(name="CV05_7", len=LCV)
CH06_7 = HKICKER(name="CH06_7", len=LCH)
CV06_7 = VKICKER(name="CV06_7", len=LCV)
CH07_7 = HKICKER(name="CH07_7", len=LCH)
CV07_7 = VKICKER(name="CV07_7", len=LCV)
CH08_7 = HKICKER(name="CH08_7", len=LCH)
CV08_7 = VKICKER(name="CV08_7", len=LCV)
CH09_7 = HKICKER(name="CH09_7", len=LCH)
CV09_7 = VKICKER(name="CV09_7", len=LCV)
CH10_7 = HKICKER(name="CH10_7", len=LCH)
CV10_7 = VKICKER(name="CV10_7", len=LCV)
CH11_7 = HKICKER(name="CH11_7", len=LCH)
CV11_7 = VKICKER(name="CV11_7", len=LCV)
CH12_7 = HKICKER(name="CH12_7", len=LCH)
CV12_7 = VKICKER(name="CV12_7", len=LCV)
CH13_7 = HKICKER(name="CH13_7", len=LCH)
CV13_7 = VKICKER(name="CV13_7", len=LCV)
CH14_7 = HKICKER(name="CH14_7", len=LCH)
CV14_7 = VKICKER(name="CV14_7", len=LCV)
CH15_7 = HKICKER(name="CH15_7", len=LCH)
CV15_7 = VKICKER(name="CV15_7", len=LCV)
CH16_7 = HKICKER(name="CH16_7", len=LCH)
CV16_7 = VKICKER(name="CV16_7", len=LCV)
CH01_9 = HKICKER(name="CH01_9", len=LCH)
CV01_9 = VKICKER(name="CV01_9", len=LCV)
CH02_9 = HKICKER(name="CH02_9", len=LCH)
CV02_9 = VKICKER(name="CV02_9", len=LCV)
CH03_9 = HKICKER(name="CH03_9", len=LCH)
CV03_9 = VKICKER(name="CV03_9", len=LCV)
CH04_9 = HKICKER(name="CH04_9", len=LCH)
CV04_9 = VKICKER(name="CV04_9", len=LCV)
CH05_9 = HKICKER(name="CH05_9", len=LCH)
CV05_9 = VKICKER(name="CV05_9", len=LCV)
CH06_9 = HKICKER(name="CH06_9", len=LCH)
CV06_9 = VKICKER(name="CV06_9", len=LCV)
CH07_9 = HKICKER(name="CH07_9", len=LCH)
CV07_9 = VKICKER(name="CV07_9", len=LCV)
CH08_9 = HKICKER(name="CH08_9", len=LCH)
CV08_9 = VKICKER(name="CV08_9", len=LCV)
CH09_9 = HKICKER(name="CH09_9", len=LCH)
CV09_9 = VKICKER(name="CV09_9", len=LCV)
CH10_9 = HKICKER(name="CH10_9", len=LCH)
CV10_9 = VKICKER(name="CV10_9", len=LCV)
CH11_9 = HKICKER(name="CH11_9", len=LCH)
CV11_9 = VKICKER(name="CV11_9", len=LCV)
CH12_9 = HKICKER(name="CH12_9", len=LCH)
CV12_9 = VKICKER(name="CV12_9", len=LCV)
CH13_9 = HKICKER(name="CH13_9", len=LCH)
CV13_9 = VKICKER(name="CV13_9", len=LCV)
CH14_9 = HKICKER(name="CH14_9", len=LCH)
CV14_9 = VKICKER(name="CV14_9", len=LCV)
CH15_9 = HKICKER(name="CH15_9", len=LCH)
CV15_9 = VKICKER(name="CV15_9", len=LCV)
CH16_9 = HKICKER(name="CH16_9", len=LCH)
CV16_9 = VKICKER(name="CV16_9", len=LCV)
CH17_9 = HKICKER(name="CH17_9", len=LCH)
CV17_9 = VKICKER(name="CV17_9", len=LCV)
CH00_11 = HKICKER(name="CH00_11", len=LCH)
CV00_11 = VKICKER(name="CV00_11", len=LCV)
CH01_11 = HKICKER(name="CH01_11", len=LCH)
CV01_11 = VKICKER(name="CV01_11", len=LCV)
CH02_11 = HKICKER(name="CH02_11", len=LCH)
CV02_11 = VKICKER(name="CV02_11", len=LCV)
CH03_11 = HKICKER(name="CH03_11", len=LCH)
CV03_11 = VKICKER(name="CV03_11", len=LCV)
CH04_11 = HKICKER(name="CH04_11", len=LCH)
CV04_11 = VKICKER(name="CV04_11", len=LCV)
CH05_11 = HKICKER(name="CH05_11", len=LCH)
CV05_11 = VKICKER(name="CV05_11", len=LCV)
CH06_11 = HKICKER(name="CH06_11", len=LCH)
CV06_11 = VKICKER(name="CV06_11", len=LCV)
CH07_11 = HKICKER(name="CH07_11", len=LCH)
CV07_11 = VKICKER(name="CV07_11", len=LCV)
CH08_11 = HKICKER(name="CH08_11", len=LCH)
CV08_11 = VKICKER(name="CV08_11", len=LCV)
CH09_11 = HKICKER(name="CH09_11", len=LCH)
CV09_11 = VKICKER(name="CV09_11", len=LCV)
CH10_11 = HKICKER(name="CH10_11", len=LCH)
CV10_11 = VKICKER(name="CV10_11", len=LCV)
CH11_11 = HKICKER(name="CH11_11", len=LCH)
CV11_11 = VKICKER(name="CV11_11", len=LCV)
CH12_11 = HKICKER(name="CH12_11", len=LCH)
CV12_11 = VKICKER(name="CV12_11", len=LCV)
CH13_11 = HKICKER(name="CH13_11", len=LCH)
CV13_11 = VKICKER(name="CV13_11", len=LCV)
CH14_11 = HKICKER(name="CH14_11", len=LCH)
CV14_11 = VKICKER(name="CV14_11", len=LCV)
CH15_11 = HKICKER(name="CH15_11", len=LCH)
CV15_11 = VKICKER(name="CV15_11", len=LCV)
CH16_11 = HKICKER(name="CH16_11", len=LCH)
CV16_11 = VKICKER(name="CV16_11", len=LCV)
CH17_11 = HKICKER(name="CH17_11", len=LCH)
CV17_11 = VKICKER(name="CV17_11", len=LCV)

# ======= Solenoids. ======================================
SOL5_6 = SOLENOID(name="SOL5_6", len=LSOL5, ks=KSOL1_6)
HSOL5_6 = SOLENOID(name="HSOL5_6", len=LSOL5/2, ks=KSOL1_6)
SOL5_8 = SOLENOID(name="SOL5_8", len=LSOL5, ks=KSOL1_8)
HSOL5_8 = SOLENOID(name="HSOL5_8", len=LSOL5/2, ks=KSOL1_8)
SOL20_6 = SOLENOID(name="SOL20_6", len=LSOL20, ks=KSOL2_6)
HSOL20_6 = SOLENOID(name="HSOL20_6", len=LSOL20/2, ks=KSOL2_6)
SOL20_8 = SOLENOID(name="SOL20_8", len=LSOL20, ks=KSOL2_8)
HSOL20_8 = SOLENOID(name="HSOL20_8", len=LSOL20/2, ks=KSOL2_8)
#
# ======= RF Cavities. ====================================
RF0 = RFCA(name="RF0", len=LRF, volt=VOLT_RF, h=HRMN_RF, freq=591.1397738e6, energy=E, philag=LAG_RF)  
# RF_CRAB = CRABCAVITY(name="RF_CRAB", len=lecrab, volt=0.0, freq=394.0e6)
RF_CRAB1 = CRABCAVITY(name="RF_CRAB1", len=lecrab, volt=crab14, freq=394.0e6)
RF_CRAB2 = CRABCAVITY(name="RF_CRAB2", len=lecrab, volt=-crab14, freq=394.0e6)
RF_CRAB3 = CRABCAVITY(name="RF_CRAB3", len=lecrab, volt=crab23, freq=394.0e6)
RF_CRAB4 = CRABCAVITY(name="RF_CRAB4", len=lecrab, volt=-crab23, freq=394.0e6)
#
# ======= Drifts. =========================================
ODB23 = DRIFT(name="ODB23", len=L12_STR)
ODBQ = DRIFT(name="ODBQ", len=LDBQ)
# IR10 Drifts.
ODBC_9 = DRIFT(name="ODBC_9", len=LDBC_9)
ODBC_10 = DRIFT(name="ODBC_10", len=LDBC_10)
ODBC_11 = DRIFT(name="ODBC_11", len=LDBC_11)
OD09_0 = DRIFT(name="OD09_0", len=LLD09_0)
OD10_0 = DRIFT(name="OD10_0", len=LLD10_0)
OD10S_0 = DRIFT(name="OD10S_0", len=8.15)
OD10L_0 = DRIFT(name="OD10L_0", len=9.3)
OD09H_0 = DRIFT(name="OD09H_0", len=6.01312)
OD10H_0 = DRIFT(name="OD10H_0", len=7.52586)
OD10RF_0 = DRIFT(name="OD10RF_0", len=0.3)
OD10IP_0 = DRIFT(name="OD10IP_0", len=0.15)
ODF_9 = DRIFT(name="ODF_9", len=LDF_9)
ODF_10 = DRIFT(name="ODF_10", len=LDF_10)
ODF20_9 = DRIFT(name="ODF20_9", len=LDF20_9)
ODF20_10 = DRIFT(name="ODF20_10", len=LDF20_10)
# IR12 Drifts.
ODBC_1 = DRIFT(name="ODBC_1", len=LDBC_1)
OD11_0 = DRIFT(name="OD11_0", len=LLD11_0)
OD11_M = DRIFT(name="OD11_M", len=0.714288)
OD11_U = DRIFT(name="OD11_U", len=5.21429)
ODMID_12 = DRIFT(name="ODMID_12", len=LLMID_12)
ODMID_12N = DRIFT(name="ODMID_12N", len=19.1 - LLMID_12)
OD12 = DRIFT(name="OD12", len=11.5)
OD12U = DRIFT(name="OD12U", len=7)
OD12_0 = DRIFT(name="OD12_0", len=LLD12_0)
ODF_11 = DRIFT(name="ODF_11", len=LDF_11)
ODF_12 = DRIFT(name="ODF_12", len=LDF_12)
ODF20_11 = DRIFT(name="ODF20_11", len=LDF20_11)
ODF20_12 = DRIFT(name="ODF20_12", len=LDF20_12)
OD12_4 = DRIFT(name="OD12_4", len=0.0975)
# IR2 Drifts.
ODBC_3 = DRIFT(name="ODBC_3", len=LDBC_3)
OLLFT0_2 = DRIFT(name="OLLFT0_2", len=LLFT0_2)
OLLFCA_2 = DRIFT(name="OLLFCA_2", len=LLFCN_2)
OLLFCB_2 = DRIFT(name="OLLFCB_2", len=12.8 - LLFCN_2)
OLLFSX_2 = DRIFT(name="OLLFSX_2", len=12.8 - LSX - LDX17)
ODF_1 = DRIFT(name="ODF_1", len=LDF_1)
ODF_2 = DRIFT(name="ODF_2", len=LDF_2)
ODF20_1 = DRIFT(name="ODF20_1", len=LDF20_1)
ODF20_2 = DRIFT(name="ODF20_2", len=LDF20_2)
# IR4 Drifts.
ODBC_5 = DRIFT(name="ODBC_5", len=LDBC_5)
ODF_3 = DRIFT(name="ODF_3", len=LDF_3)
ODF_4 = DRIFT(name="ODF_4", len=LDF_4)
ODF20_3 = DRIFT(name="ODF20_3", len=LDF20_3)
ODF20_4 = DRIFT(name="ODF20_4", len=LDF20_4)
ODB12 = DRIFT(name="ODB12", len=0.3)
OD4_3 = DRIFT(name="OD4_3", len=LD4_3)
OD3_3 = DRIFT(name="OD3_3", len=0.535)
OD2_3 = DRIFT(name="OD2_3", len=6)
OD1_3 = DRIFT(name="OD1_3", len=23.2385/2)
OD1_4 = DRIFT(name="OD1_4", len=23.2385/2)
OD2_4 = DRIFT(name="OD2_4", len=6)
OD3_4 = DRIFT(name="OD3_4", len=1.2585)
OD4_4 = DRIFT(name="OD4_4", len=LD4_4)
OD5_4 = DRIFT(name="OD5_4", len=0.535)
# Arc Drifts.
ODF29 = DRIFT(name="ODF29", len=LDF29)
ODX17 = DRIFT(name="ODX17", len=LDX17)
OHQC = DRIFT(name="OHQC", len=LHQC)
ODC_1 = DRIFT(name="ODC_1", len=LDC_1)
ODC_3 = DRIFT(name="ODC_3", len=LDC_3)
ODC_5 = DRIFT(name="ODC_5", len=LDC_5)
ODC_7 = DRIFT(name="ODC_7", len=LDC_7)
ODC_9 = DRIFT(name="ODC_9", len=LDC_9)
ODC_11 = DRIFT(name="ODC_11", len=LDC_11)
ODSX = DRIFT(name="ODSX", len=LDSX)
ODSXL = DRIFT(name="ODSXL", len=LDSXL)
OQSX = DRIFT(name="OQSX", len=LQSX)
OQSXL = DRIFT(name="OQSXL", len=LQSXL)
OSSX = DRIFT(name="OSSX", len=LSSX)
OSS1 = DRIFT(name="OSS1", len=LSS1)
# Five-bend drifts.
ODB = DRIFT(name="ODB", len=L12_STR)
ODF4B_5 = DRIFT(name="ODF4B_5", len=LDF4B_5)
ODF4B_6 = DRIFT(name="ODF4B_6", len=LDF4B_6)
ODF4B_7 = DRIFT(name="ODF4B_7", len=LDF4B_7)
ODF4B_8 = DRIFT(name="ODF4B_8", len=LDF4B_8)
# IR6 and IR8 drifts.
oecrab = DRIFT(name="oecrab", len=4.0)
oqc = DRIFT(name="oqc", len=0.3)
oww = DRIFT(name="oww", len=0.4)
oww_ir6f = DRIFT(name="oww_ir6f", len=0.6)
oww_sh = DRIFT(name="oww_sh", len=0.15)
oww_ir6f_sh = DRIFT(name="oww_ir6f_sh", len=oww_ir6f.len - oww_sh.len - sq.len)
ow2c = DRIFT(name="ow2c", len=1.0)
omir = DRIFT(name="omir", len=lr_rot_bend + 2*oww.len)
omir_ir6f = DRIFT(name="omir_ir6f", len=lr_rot_bend + 2*oww_ir6f.len - oww_sh.len - sq.len)
omir_ir6f2 = DRIFT(name="omir_ir6f2", len=0.5)
omir_ir6f1 = DRIFT(name="omir_ir6f1", len=omir_ir6f.len - omir_ir6f2.len)
omir_ir6f4 = DRIFT(name="omir_ir6f4", len=0.25)
omir_ir6f3 = DRIFT(name="omir_ir6f3", len=omir_ir6f.len - omir_ir6f4.len)
omir_ir6r = DRIFT(name="omir_ir6r", len=lr_rot_bend + 2*oww_ir6f.len - oww_sh.len - sq.len)
omir_ir6r2 = DRIFT(name="omir_ir6r2", len=4.01334)
OMIR51 = DRIFT(name="OMIR51", len=lr_rot_bend+1.0)
OMIR52 = DRIFT(name="OMIR52", len=lr_rot_bend-1.0)
olaserip = DRIFT(name="olaserip", len=2 - oww_sh.len - sq.len)
o0ef_6 = DRIFT(name="o0ef_6", len=5.8)
o1ef_6 = DRIFT(name="o1ef_6", len=3.76)
o2ef_6 = DRIFT(name="o2ef_6", len=17.9213)
o3ef_6 = DRIFT(name="o3ef_6", len=L3EF_6)
o0er_6 = DRIFT(name="o0er_6", len=beg_q1er)
o1er_6 = DRIFT(name="o1er_6", len=beg_q2er-end_q1er)
o2er_6 = DRIFT(name="o2er_6", len=beg_d2er-end_q2er)
o3er_6 = DRIFT(name="o3er_6", len=22.1717)
o4er_6 = DRIFT(name="o4er_6", len=L4ER_6)
o0ef_8 = DRIFT(name="o0ef_8", len=5.8)
o1ef_8 = DRIFT(name="o1ef_8", len=3.76)
o2ef_8 = DRIFT(name="o2ef_8", len=2.12718e+01)
o3ef_8 = DRIFT(name="o3ef_8", len=L3EF_8)
o0er_8 = DRIFT(name="o0er_8", len=beg_q1er)
o1er_8 = DRIFT(name="o1er_8", len=beg_q2er-end_q1er)
o2er_8 = DRIFT(name="o2er_8", len=beg_d2er-end_q2er)
o3er_8 = DRIFT(name="o3er_8", len=beg_q3er-end_d2er)
o4er_8 = DRIFT(name="o4er_8", len=L4ER_8)
OSOL85 = DRIFT(name="OSOL85", len=0.85)
OSOL25 = DRIFT(name="OSOL25", len=0.25)
O5A_6f = DRIFT(name="O5A_6f", len=LO5A_6f)
O5A_6r = DRIFT(name="O5A_6r", len=LO5A_6r)
O5A_8f = DRIFT(name="O5A_8f", len=LO5A_8f)
O5A_8r = DRIFT(name="O5A_8r", len=LO5A_8r)
#
# ======= Markers. ========================================
MFF_5 = MARKER(name="MFF_5")
MFF_6 = MARKER(name="MFF_6")
MFF_7 = MARKER(name="MFF_7")
MFF_8 = MARKER(name="MFF_8")
mlrf_6 = MARKER(name="mlrf_6")
mlrr_6 = MARKER(name="mlrr_6")
mlrf_8 = MARKER(name="mlrf_8")
mlrr_8 = MARKER(name="mlrr_8")
MQSS_10 = MARKER(name="MQSS_10")
MQSS_12 = MARKER(name="MQSS_12")
MQSS_2 = MARKER(name="MQSS_2")

MQ1 = MARKER(name="MQ1")
MQ2 = MARKER(name="MQ2")
MQ3 = MARKER(name="MQ3")
MQ4 = MARKER(name="MQ4")
MQ5 = MARKER(name="MQ5")
MQ6 = MARKER(name="MQ6")
MQ7 = MARKER(name="MQ7")
MQ8 = MARKER(name="MQ8")
MQ9 = MARKER(name="MQ9")
MQ10 = MARKER(name="MQ10")
MQ11 = MARKER(name="MQ11")
MQ12 = MARKER(name="MQ12")
MQ13 = MARKER(name="MQ13")
MQ14 = MARKER(name="MQ14")
MQ15 = MARKER(name="MQ15")
MQ16 = MARKER(name="MQ16")
MQ17 = MARKER(name="MQ17")
MQ18 = MARKER(name="MQ18")
MQ19 = MARKER(name="MQ19")
MQ20 = MARKER(name="MQ20")
MQ21 = MARKER(name="MQ21")
MQ22 = MARKER(name="MQ22")
MROT1 = MARKER(name="MROT1")
MROT2 = MARKER(name="MROT2")
MROT3 = MARKER(name="MROT3")
MROT4 = MARKER(name="MROT4")
IP6 = MARKER(name="IP6")
IP8 = MARKER(name="IP8")
IP10 = MARKER(name="IP10")
IP12 = MARKER(name="IP12")
IP4 = MARKER(name="IP4")
IP2 = MARKER(name="IP2")
MBEG = MARKER(name="MBEG")
MMID = MARKER(name="MMID")
MEND = MARKER(name="MEND")
MARC_BEG = MARKER(name="MARC_BEG")
MARC_END = MARKER(name="MARC_END")


MKICK_INJ = MARKER(name="MKICK_INJ")
MLAMB = MARKER(name="MLAMB")
MCOLL_INJ = MARKER(name="MCOLL_INJ")
MCOLL_H1 = MARKER(name="MCOLL_H1")
MCOLL_H2 = MARKER(name="MCOLL_H2")
MCOLL_H3 = MARKER(name="MCOLL_H3")
MCOLL_V1 = MARKER(name="MCOLL_V1")
MCOLL_V2 = MARKER(name="MCOLL_V2")
MCOLL_V3 = MARKER(name="MCOLL_V3")
MCOLL_MASK = MARKER(name="MCOLL_MASK")

M_EDETECT = MARKER(name="M_EDETECT")
M_EDETECT1 = MARKER(name="M_EDETECT1")

IR_6 = [mlrf_6,
   q14ef_6, oww_sh, sq14ef_6, oww_ir6f_sh, d5ef_6, oww_ir6f,
   q13ef_6, oww_sh, sq13ef_6, oww_ir6f_sh, d5ef_6, oww_ir6f,
   q12ef_6, oww_sh, sq12ef_6, oww_ir6f_sh, d5ef_6, oww_ir6f,
   q11ef_6, oww_sh, sq11ef_6, oww_ir6f_sh, d5ef_6, oww_ir6f,
   q10ef_6, oww_sh, sq10ef_6, oww_ir6f_sh, d4ef_6, oww_ir6f,
   q9ef_6, oww_sh, sq9ef_6, olaserip, d4ef_6, oww_ir6f, 
   q8ef_6, oww_sh, sq8ef_6, oww_ir6f_sh, d3ef_6, oww_ir6f,
   q7ef_6, oww_sh, sq7ef_6, omir_ir6f3, M_EDETECT1, omir_ir6f4,
   q6ef_6, oww_sh, sq6ef_6, omir_ir6f1, M_EDETECT, omir_ir6f2,
   q5ef_6, oww_sh, sq5ef_6, omir_ir6f,
   q4ef_6, oww_sh, sq4ef_6, o3ef_6,
   q3ef_6, ow2c,
   RF_CRAB1,
   ow2c, q2ef_6, oww_sh, sq2ef_6,
   oww_sh, d2ef_6, ODB23, d2ef_6, ODB23, d1ef_6, o2ef_6,
   MCOLL_MASK,
   Q1EF_6,
   o1ef_6,
   Q0EF_6,
   o0ef_6]
IR_6N = [o0er_6, q1er_6, o1er_6, 
 q2er_6, o2er_6, d2er_6, o3er_6, 
 sq3er_6, oww_sh, q3er_6, omir_ir6r,
 sq4er_6, oww_sh, q4er_6, omir_ir6r,
 sq5er_6, oww_sh, q5er_6, oww,  
 d3er_6, oww,
 sq6er_6, oww_sh, q6er_6, omir_ir6r2,
 sq7er_6, oww_sh, q7er_6, omir_ir6r2,
 sq8er_6, oww_sh, q8er_6, ow2c, RF_CRAB4, ow2c,
 sq9er_6, oww_sh, q9er_6, omir_ir6r2,
 sq10er_6, oww_sh, q10er_6, oww, d4er_6, oww,
 q11er_6, oww, d4er_6, oww,
 q12er_6, oww, d4er_6, oww,
 q13er_6, oww, d4er_6, oww,
 q14er_6, oww, d4er_6, oww,
 q15er_6, mlrr_6]

IR_8 = [mlrf_8,
   q15ef_8, oww, d3ef_8, oww,
   q14ef_8, oww, d3ef_8, oww,
   q13ef_8, oww, d3ef_8, oww,
   q12ef_8, oww, d3ef_8, oww,
   q11ef_8, oww, d2ef_8, oww,
   q10ef_8, oww, d2ef_8, oww,
   q9ef_8, omir,
   q8ef_8, omir,
   q7ef_8, omir,
   q6ef_8, omir,
   q5ef_8, omir,
   q4ef_8, o3ef_8,
   q3ef_8, oqc,
   RF_CRAB2,
   oqc, q2ef_8,
   oww, d1ef_8,
   oww, d1ef_8,
   o2ef_8,
   q1ef_8,
   o1ef_8,
   q0ef_8,
   o0ef_8]

IR_8N = [o0er_8, q1er_8, o1er_8, q2er_8, o2er_8, d2er_8, o3er_8, q3er_8,
 oww, d3er_8, o4er_8,
 q4er_8, OMIR51,
 q5er_8, OMIR52,
 q6er_8, omir,
 q7er_8, omir,
 q8er_8, omir,
 q9er_8, oqc, RF_CRAB3, oqc,
 q10er_8, omir,
 q11er_8, oww, d4er_8, oww,
 q12er_8, oww, d5er_8, oww,
 q13er_8, oww, d5er_8, oww,
 q14er_8, oww, d5er_8, oww,
 q15er_8, mlrr_8]


# Solenoid sections (NSS5 and NSS20).
ROT5_5 = [OSOL85, HQSS1_5, HQSS1_5,
               OSOL25, HQSS2_5, HQSS2_5,
               OSOL25, HQSS3_5, HQSS3_5,
               OSOL25, HQSS4_5, HQSS4_5,
               OSOL25, HQSS5_5, HQSS5_5, OSOL85]
NSS5_5 = [MROT1, HSOL5_6, HSOL5_6, ROT5_5..., HSOL5_6, HSOL5_6, MROT2]
ROT5_6 = [OSOL85, HQSS1_6, HQSS1_6,
               OSOL25, HQSS2_6, HQSS2_6,
               OSOL25, HQSS3_6, HQSS3_6,
               OSOL25, HQSS4_6, HQSS4_6,
               OSOL25, HQSS5_6, HQSS5_6, OSOL85]
NSS5_6 = [MROT1, HSOL5_6, HSOL5_6, ROT5_6..., HSOL5_6, HSOL5_6, MROT2]
ROT20_5 = [OSOL85, HQLS1_5, HQLS1_5,
               OSOL25, HQLS2_5, HQLS2_5,
               OSOL25, HQLS3_5, HQLS3_5,
               OSOL25, HQLS4_5, HQLS4_5,
               OSOL25, HQLS5_5, HQLS5_5,
               OSOL25, HQLS6_5, HQLS6_5,
               OSOL25, HQLS7_5, HQLS7_5, OSOL85]
NSS20_5 = [MROT3, HSOL20_6, HSOL20_6, ROT20_5..., HSOL20_6, HSOL20_6, MROT4]
ROT20_6 = [OSOL85, HQLS1_6, HQLS1_6,
               OSOL25, HQLS2_6, HQLS2_6,
               OSOL25, HQLS3_6, HQLS3_6,
               OSOL25, HQLS4_6, HQLS4_6,
               OSOL25, HQLS5_6, HQLS5_6,
               OSOL25, HQLS6_6, HQLS6_6,
               OSOL25, HQLS7_6, HQLS7_6, OSOL85]
NSS20_6 = [MROT3, HSOL20_6, HSOL20_6, ROT20_6..., HSOL20_6, HSOL20_6, MROT4]
ROT5_7 = [OSOL85, HQSS1_7, HQSS1_7,
               OSOL25, HQSS2_7, HQSS2_7,
               OSOL25, HQSS3_7, HQSS3_7,
               OSOL25, HQSS4_7, HQSS4_7,
               OSOL25, HQSS5_7, HQSS5_7, OSOL85]
NSS5_7 = [MROT1, HSOL5_8, HSOL5_8, ROT5_7..., HSOL5_8, HSOL5_8, MROT2]
ROT5_8 = [OSOL85, HQSS1_8, HQSS1_8,
               OSOL25, HQSS2_8, HQSS2_8,
               OSOL25, HQSS3_8, HQSS3_8,
               OSOL25, HQSS4_8, HQSS4_8,
               OSOL25, HQSS5_8, HQSS5_8, OSOL85]
NSS5_8 = [MROT1, HSOL5_8, HSOL5_8, ROT5_8..., HSOL5_8, HSOL5_8, MROT2]
ROT20_7 = [OSOL85, HQLS1_7, HQLS1_7,
               OSOL25, HQLS2_7, HQLS2_7,
               OSOL25, HQLS3_7, HQLS3_7,
               OSOL25, HQLS4_7, HQLS4_7,
               OSOL25, HQLS5_7, HQLS5_7,
               OSOL25, HQLS6_7, HQLS6_7,
               OSOL25, HQLS7_7, HQLS7_7, OSOL85]
NSS20_7 = [MROT3, HSOL20_8, HSOL20_8, ROT20_7..., HSOL20_8, HSOL20_8, MROT4]
ROT20_8 = [OSOL85, HQLS1_8, HQLS1_8,
               OSOL25, HQLS2_8, HQLS2_8,
               OSOL25, HQLS3_8, HQLS3_8,
               OSOL25, HQLS4_8, HQLS4_8,
               OSOL25, HQLS5_8, HQLS5_8,
               OSOL25, HQLS6_8, HQLS6_8,
               OSOL25, HQLS7_8, HQLS7_8, OSOL85]
NSS20_8 = [MROT3, HSOL20_8, HSOL20_8, ROT20_8..., HSOL20_8, HSOL20_8, MROT4]

FIVEBND_5 = [QFF1_5, QFF1_5, ODF4B_5, DB23_5, ODF4B_5,
        QFF2_5, QFF2_5, ODF4B_5, DB23_5, ODF4B_5,
        QFF3_5, QFF3_5, ODF4B_5, DB23_5, ODF4B_5,
        QFF4_5, QFF4_5, ODF4B_5, DB23_5, ODF4B_5,
        QFF5_5, QFF5_5, ODF4B_5, DB23_5, ODF4B_5,
        QFF6_5, QFF6_5]
FIVEBND_6 = [QFF1_6, QFF1_6, ODF4B_6, DB23_6, ODF4B_6,
        QFF2_6, QFF2_6, ODF4B_6, DB23_6, ODF4B_6,
        QFF3_6, QFF3_6, ODF4B_6, DB23_6, ODF4B_6,
        QFF4_6, QFF4_6, ODF4B_6, DB23_6, ODF4B_6,
        QFF5_6, QFF5_6, ODF4B_6, DB23_6, ODF4B_6,
        QFF6_6, QFF6_6]
FIVEBND_7 = [QFF1_7, QFF1_7, ODF4B_7, DB23_7, ODF4B_7,
        QFF2_7, QFF2_7, ODF4B_7, DB23_7, ODF4B_7,
        QFF3_7, QFF3_7, ODF4B_7, DB23_7, ODF4B_7,
        QFF4_7, QFF4_7, ODF4B_7, DB23_7, ODF4B_7,
        QFF5_7, QFF5_7, ODF4B_7, DB23_7, ODF4B_7,
        QFF6_7, QFF6_7]
FIVEBND_8 = [QFF1_8, QFF1_8, ODF4B_8, DB23_8, ODF4B_8,
        QFF2_8, QFF2_8, ODF4B_8, DB23_8, ODF4B_8,
        QFF3_8, QFF3_8, ODF4B_8, DB23_8, ODF4B_8,
        QFF4_8, QFF4_8, ODF4B_8, DB23_8, ODF4B_8,
        QFF5_8, QFF5_8, ODF4B_8, DB23_8, ODF4B_8,
        QFF6_8, QFF6_8]

HALF_IR6  = [HQF_5a, O5A_6f, HQD_5a, HQD_5a, OSOL85, NSS5_5..., OSOL85,
   FIVEBND_5..., OSOL85, NSS20_5..., OSOL85]
HALF_IR6N  = [HQF_6a, O5A_6r, HQD_6a, HQD_6a, OSOL85, NSS5_6..., OSOL85,
   FIVEBND_6..., OSOL85, NSS20_6..., OSOL85]
HALF_IR8  =[HQF_7a, O5A_8f, HQD_7a, HQD_7a, OSOL85, NSS5_7..., OSOL85,
   FIVEBND_7..., OSOL85, NSS20_7..., OSOL85]
HALF_IR8N =[HQF_8a, O5A_8r, HQD_8a, HQD_8a, OSOL85, NSS5_8..., OSOL85,
   FIVEBND_8..., OSOL85, NSS20_8..., OSOL85]


function FODOCELL_1(SX1, SX2, HQF1, HQD, HQF2, DCCN, CH, CV, ODC)
    global OHQC, ODSX, OSSX, OQSX, OSS1
    return [HQF1, OHQC, CH, ODC, DCCN..., ODSX, SX2, OSSX, SX2, OQSX, HQD, HQD, OHQC, CV, ODC, DCCN...,
    ODSX, OSS1, SX1, OSS1, OQSX, HQF2]
end
function FODOCELL_3(SX1, SX2, HQF1, HQD, HQF2, DCCN, CH, CV, ODC)
    global OHQC, ODSX, OSSX, OQSX, OSS1
    return [HQF1, OQSX, OSS1, SX1, OSS1, ODSX, DCCN..., ODC, CV, OHQC, HQD, HQD, OQSX, SX2, OSSX, SX2, ODSX, DCCN...,
    ODC, CH, OHQC, HQF2]
end
function FODOCELL_5(SX1, SX2, HQF1, HQD, HQF2, DCCN, CH, CV, ODC)
    global OHQC, ODSX, OSSX, OQSX, OQSXL, ODSXL
    return [HQF1, OHQC, CH, ODC, DCCN..., ODSXL, SX2, OQSXL, HQD, HQD, OHQC, CV, ODC, DCCN...,
    ODSX, SX1, OSSX, SX1, OQSX, HQF2]
end
function FODOCELL_7(SX1, SX2, HQF1, HQD, HQF2, DCCN, CH, CV, ODC)
    global OHQC, ODSX, OSSX, OQSX, OQSXL, ODSXL
    return [HQF1, OQSX, SX1, OSSX, SX1, ODSX, DCCN..., ODC, CV, OHQC, HQD, HQD, OQSXL, SX2, ODSXL, DCCN...,
    ODC, CH, OHQC, HQF2]
end
function FODOCELL_9(SX1, SX2, HQF1, HQD, HQF2, DCCN, CH, CV, ODC)
    global OHQC, ODSX, OSSX, OQSX, OQSXL, ODSXL
    return [HQF1, OHQC, CH, ODC, DCCN..., ODSXL, SX2, OQSXL, HQD, HQD, OHQC, CV, ODC, DCCN...,
    ODSX, SX1, OSSX, SX1, OQSX, HQF2]
end
function FODOCELL_11(SX1, SX2, HQF1, HQD, HQF2, DCCN, CH, CV, ODC)
    global OHQC, ODSX, OSSX, OQSX, OSS1
    return [HQF1, OQSX, OSS1, SX1, OSS1, ODSX, DCCN..., ODC, CV, OHQC, HQD, HQD, OQSX, SX2, OSSX, SX2, ODSX, DCCN...,
    ODC, CH, OHQC, HQF2]
end

ArcCell01_1 = FODOCELL_1(SF01_1, SD01_1,  HQF_1,   HQD_1,   HQF_1,   DCCN,      CH01_1, CV01_1, ODC_1)
ArcCell02_1 = FODOCELL_1(SF02_1, SD02_1,  HQF_1,   HQD_1,   HQF_1,   DCCN,      CH02_1, CV02_1, ODC_1)
ArcCell03_1 = FODOCELL_1(SF03_1, SD03_1,  HQF_1,   HQD_1,   HQF_1,   DCCN,      CH03_1, CV03_1, ODC_1)
ArcCell04_1 = FODOCELL_1(SF04_1, SD04_1,  HQF_1,   HQD_1,   HQF_1,   DCCN,      CH04_1, CV04_1, ODC_1)
ArcCell05_1 = FODOCELL_1(SF05_1, SD05_1,  HQF_1,   HQD_1,   HQF_1,   DCCN,      CH05_1, CV05_1, ODC_1)
ArcCell06_1 = FODOCELL_1(SF06_1, SD06_1,  HQF_1,   HQD_1,   HQF_1,   DCCN,      CH06_1, CV06_1, ODC_1)
ArcCell07_1 = FODOCELL_1(SF07_1, SD07_1,  HQF_1,   HQD_1,   HQF_1,   DCCN,      CH07_1, CV07_1, ODC_1)
ArcCell08_1 = FODOCELL_1(SF08_1, SD08_1,  HQF_1,   HQD_1,   HQF_1,   DCCN,      CH08_1, CV08_1, ODC_1)
ArcCell09_1 = FODOCELL_1(SF09_1, SD09_1,  HQF_1,   HQD_1,   HQF_1,   DCCN,      CH09_1, CV09_1, ODC_1)
ArcCell10_1 = FODOCELL_1(SF10_1, SD10_1,  HQF_1,   HQD_1,   HQF_1,   DCCN,      CH10_1, CV10_1, ODC_1)
ArcCell11_1 = FODOCELL_1(SF11_1, SD11_1,  HQF_1,   HQD_1,   HQF_1,   DCCN,      CH11_1, CV11_1, ODC_1)
ArcCell12_1 = FODOCELL_1(SF12_1, SD12_1,  HQF_1,   HQD_1,   HQF_1,   DCCN,      CH12_1, CV12_1, ODC_1)
ArcCell13_1 = FODOCELL_1(SF13_1, SD13_1,  HQF_1,   HQD_1,   HQF_1,   DCCN,      CH13_1, CV13_1, ODC_1)
ArcCell14_1 = FODOCELL_1(SF14_1, SD14_1,  HQF_1,   HQD_1,   HQF_1,   DCCN,      CH14_1, CV14_1, ODC_1)
ArcCell15_1 = FODOCELL_1(SF15_1, SD15_1,  HQF_1,   HQD_1,   HQF_1,   DCCN,      CH15_1, CV15_1, ODC_1)
ArcCell16_1 = FODOCELL_1(SF16_1, SD16_1,  HQF_1,   HQD_1,   HQF_1,   DCCN,      CH16_1, CV16_1, ODC_1)

ArcCell01_3 = FODOCELL_3(SF01_3, SD01_3,  HQF_3,   HQD_3,   HQF_3,   DCCN,      CH01_3, CV01_3, ODC_3)
ArcCell02_3 = FODOCELL_3(SF02_3, SD02_3,  HQF_3,   HQD_3,   HQF_3,   DCCN,      CH02_3, CV02_3, ODC_3)
ArcCell03_3 = FODOCELL_3(SF03_3, SD03_3,  HQF_3,   HQD_3,   HQF_3,   DCCN,      CH03_3, CV03_3, ODC_3)
ArcCell04_3 = FODOCELL_3(SF04_3, SD04_3,  HQF_3,   HQD_3,   HQF_3,   DCCN,      CH04_3, CV04_3, ODC_3)
ArcCell05_3 = FODOCELL_3(SF05_3, SD05_3,  HQF_3,   HQD_3,   HQF_3,   DCCN,      CH05_3, CV05_3, ODC_3)
ArcCell06_3 = FODOCELL_3(SF06_3, SD06_3,  HQF_3,   HQD_3,   HQF_3,   DCCN,      CH06_3, CV06_3, ODC_3)
ArcCell07_3 = FODOCELL_3(SF07_3, SD07_3,  HQF_3,   HQD_3,   HQF_3,   DCCN,      CH07_3, CV07_3, ODC_3)
ArcCell08_3 = FODOCELL_3(SF08_3, SD08_3,  HQF_3,   HQD_3,   HQF_3,   DCCN,      CH08_3, CV08_3, ODC_3)
ArcCell09_3 = FODOCELL_3(SF09_3, SD09_3,  HQF_3,   HQD_3,   HQF_3,   DCCN,      CH09_3, CV09_3, ODC_3)
ArcCell10_3 = FODOCELL_3(SF10_3, SD10_3,  HQF_3,   HQD_3,   HQF_3,   DCCN,      CH10_3, CV10_3, ODC_3)
ArcCell11_3 = FODOCELL_3(SF11_3, SD11_3,  HQF_3,   HQD_3,   HQF_3,   DCCN,      CH11_3, CV11_3, ODC_3)
ArcCell12_3 = FODOCELL_3(SF12_3, SD12_3,  HQF_3,   HQD_3,   HQF_3,   DCCN,      CH12_3, CV12_3, ODC_3)
ArcCell13_3 = FODOCELL_3(SF13_3, SD13_3,  HQF_3,   HQD_3,   HQF_3,   DCCN,      CH13_3, CV13_3, ODC_3)
ArcCell14_3 = FODOCELL_3(SF14_3, SD14_3,  HQF_3,   HQD_3,   HQF_3,   DCCN,      CH14_3, CV14_3, ODC_3)
ArcCell15_3 = FODOCELL_3(SF15_3, SD15_3,  HQF_3,   HQD_3,   HQF_3,   DCCN,      CH15_3, CV15_3, ODC_3)
ArcCell16_3 = FODOCELL_3(SF16_3, SD16_3,  HQF_3,   HQD_3,   HQF_3,   DCCN,      CH16_3, CV16_3, ODC_3)

ArcCell01_5 = FODOCELL_5(SF01_5, SD01_5,  HQF_5,   HQD_5,   HQF_5,   DCCN,      CH01_5, CV01_5, ODC_5)
ArcCell02_5 = FODOCELL_5(SF02_5, SD02_5,  HQF_5,   HQD_5,   HQF_5,   DCCN,      CH02_5, CV02_5, ODC_5)
ArcCell03_5 = FODOCELL_5(SF03_5, SD03_5,  HQF_5,   HQD_5,   HQF_5,   DCCN,      CH03_5, CV03_5, ODC_5)
ArcCell04_5 = FODOCELL_5(SF04_5, SD04_5,  HQF_5,   HQD_5,   HQF_5,   DCCN,      CH04_5, CV04_5, ODC_5)
ArcCell05_5 = FODOCELL_5(SF05_5, SD05_5,  HQF_5,   HQD_5,   HQF_5,   DCCN,      CH05_5, CV05_5, ODC_5)
ArcCell06_5 = FODOCELL_5(SF06_5, SD06_5,  HQF_5,   HQD_5,   HQF_5,   DCCN,      CH06_5, CV06_5, ODC_5)
ArcCell07_5 = FODOCELL_5(SF07_5, SD07_5,  HQF_5,   HQD_5,   HQF_5,   DCCN,      CH07_5, CV07_5, ODC_5)
ArcCell08_5 = FODOCELL_5(SF08_5, SD08_5,  HQF_5,   HQD_5,   HQF_5,   DCCN,      CH08_5, CV08_5, ODC_5)
ArcCell09_5 = FODOCELL_5(SF09_5, SD09_5,  HQF_5,   HQD_5,   HQF_5,   DCCN,      CH09_5, CV09_5, ODC_5)
ArcCell10_5 = FODOCELL_5(SF10_5, SD10_5,  HQF_5,   HQD_5,   HQF_5,   DCCN,      CH10_5, CV10_5, ODC_5)
ArcCell11_5 = FODOCELL_5(SF11_5, SD11_5,  HQF_5,   HQD_5,   HQF_5,   DCCN,      CH11_5, CV11_5, ODC_5)
ArcCell12_5 = FODOCELL_5(SF12_5, SD12_5,  HQF_5,   HQD_5,   HQF_5,   DCCN,      CH12_5, CV12_5, ODC_5)
ArcCell13_5 = FODOCELL_5(SF13_5, SD13_5,  HQF_5,   HQD_5,   HQF_5,   DCCN,      CH13_5, CV13_5, ODC_5)
ArcCell14_5 = FODOCELL_5(SF14_5, SD14_5,  HQF_5,   HQD_5,   HQF_5c,  DCCN,      CH14_5, CV14_5, ODC_5)
ArcCell15_5 = FODOCELL_5(SF15_5, SD15_5,  HQF_5c,  HQD_5c,  HQF_5b,  DMNS_6,    CH15_5, CV15_5, ODC_5)
ArcCell16_5 = FODOCELL_5(SF16_5, SD16_5,  HQF_5b,  HQD_5b,  HQF_5a,  DMNS_6,    CH16_5, CV16_5, ODC_5)

ArcCell01_7 = FODOCELL_7(SF01_7, SD01_7,  HQF_6a,  HQD_6b,  HQF_6b,  DPLS_6,    CH01_7, CV01_7, ODC_7)
ArcCell02_7 = FODOCELL_7(SF02_7, SD02_7,  HQF_6b,  HQD_6c,  HQF_6c,  DPLS_6,    CH02_7, CV02_7, ODC_7)
ArcCell03_7 = FODOCELL_7(SF03_7, SD03_7,  HQF_6c,  HQD_7,   HQF_7,   DCCN,      CH03_7, CV03_7, ODC_7)
ArcCell04_7 = FODOCELL_7(SF04_7, SD04_7,  HQF_7,   HQD_7,   HQF_7,   DCCN,      CH04_7, CV04_7, ODC_7)
ArcCell05_7 = FODOCELL_7(SF05_7, SD05_7,  HQF_7,   HQD_7,   HQF_7,   DCCN,      CH05_7, CV05_7, ODC_7)
ArcCell06_7 = FODOCELL_7(SF06_7, SD06_7,  HQF_7,   HQD_7,   HQF_7,   DCCN,      CH06_7, CV06_7, ODC_7)
ArcCell07_7 = FODOCELL_7(SF07_7, SD07_7,  HQF_7,   HQD_7,   HQF_7,   DCCN,      CH07_7, CV07_7, ODC_7)
ArcCell08_7 = FODOCELL_7(SF08_7, SD08_7,  HQF_7,   HQD_7,   HQF_7,   DCCN,      CH08_7, CV08_7, ODC_7)
ArcCell09_7 = FODOCELL_7(SF09_7, SD09_7,  HQF_7,   HQD_7,   HQF_7,   DCCN,      CH09_7, CV09_7, ODC_7)
ArcCell10_7 = FODOCELL_7(SF10_7, SD10_7,  HQF_7,   HQD_7,   HQF_7,   DCCN,      CH10_7, CV10_7, ODC_7)
ArcCell11_7 = FODOCELL_7(SF11_7, SD11_7,  HQF_7,   HQD_7,   HQF_7,   DCCN,      CH11_7, CV11_7, ODC_7)
ArcCell12_7 = FODOCELL_7(SF12_7, SD12_7,  HQF_7,   HQD_7,   HQF_7,   DCCN,      CH12_7, CV12_7, ODC_7)
ArcCell13_7 = FODOCELL_7(SF13_7, SD13_7,  HQF_7,   HQD_7,   HQF_7,   DCCN,      CH13_7, CV13_7, ODC_7)
ArcCell14_7 = FODOCELL_7(SF14_7, SD14_7,  HQF_7,   HQD_7,   HQF_7c,  DCCN,      CH14_7, CV14_7, ODC_7)
ArcCell15_7 = FODOCELL_7(SF15_7, SD15_7,  HQF_7c,  HQD_7c,  HQF_7b,  DPLS_8,    CH15_7, CV15_7, ODC_7)
ArcCell16_7 = FODOCELL_7(SF16_7, SD16_7,  HQF_7b,  HQD_7b,  HQF_7a,  DPLS_8,    CH16_7, CV16_7, ODC_7)

ArcCell01_9 = FODOCELL_9(SF01_9, SD01_9,  HQF_8a,  HQD_8b,  HQF_8b,  DMNS_8,    CH01_9, CV01_9, ODC_9)
ArcCell02_9 = FODOCELL_9(SF02_9, SD02_9,  HQF_8b,  HQD_8c,  HQF_8c,  DMNS_8,    CH02_9, CV02_9, ODC_9)
ArcCell03_9 = FODOCELL_9(SF03_9, SD03_9,  HQF_8c,  HQD_9,   HQF_9,   DCCN,      CH03_9, CV03_9, ODC_9)
ArcCell04_9 = FODOCELL_9(SF04_9, SD04_9,  HQF_9,   HQD_9,   HQF_9,   DCCN,      CH04_9, CV04_9, ODC_9)
ArcCell05_9 = FODOCELL_9(SF05_9, SD05_9,  HQF_9,   HQD_9,   HQF_9,   DCCN,      CH05_9, CV05_9, ODC_9)
ArcCell06_9 = FODOCELL_9(SF06_9, SD06_9,  HQF_9,   HQD_9,   HQF_9,   DCCN,      CH06_9, CV06_9, ODC_9)
ArcCell07_9 = FODOCELL_9(SF07_9, SD07_9,  HQF_9,   HQD_9,   HQF_9,   DCCN,      CH07_9, CV07_9, ODC_9)
ArcCell08_9 = FODOCELL_9(SF08_9, SD08_9,  HQF_9,   HQD_9,   HQF_9,   DCCN,      CH08_9, CV08_9, ODC_9)
ArcCell09_9 = FODOCELL_9(SF09_9, SD09_9,  HQF_9,   HQD_9,   HQF_9,   DCCN,      CH09_9, CV09_9, ODC_9)
ArcCell10_9 = FODOCELL_9(SF10_9, SD10_9,  HQF_9,   HQD_9,   HQF_9,   DCCN,      CH10_9, CV10_9, ODC_9)
ArcCell11_9 = FODOCELL_9(SF11_9, SD11_9,  HQF_9,   HQD_9,   HQF_9,   DCCN,      CH11_9, CV11_9, ODC_9)
ArcCell12_9 = FODOCELL_9(SF12_9, SD12_9,  HQF_9,   HQD_9,   HQF_9,   DCCN,      CH12_9, CV12_9, ODC_9)
ArcCell13_9 = FODOCELL_9(SF13_9, SD13_9,  HQF_9,   HQD_9,   HQF_9,   DCCN,      CH13_9, CV13_9, ODC_9)
ArcCell14_9 = FODOCELL_9(SF14_9, SD14_9,  HQF_9,   HQD_9,   HQF_9,   DCCN,      CH14_9, CV14_9, ODC_9)
ArcCell15_9 = FODOCELL_9(SF15_9, SD15_9,  HQF_9,   HQD_9,   HQF_9,   DCCN,      CH15_9, CV15_9, ODC_9)
ArcCell16_9 = FODOCELL_9(SF16_9, SD16_9,  HQF_9,   HQD_9,   HQF_9,   DCCN,      CH16_9, CV16_9, ODC_9)

ArcCell01_11 = FODOCELL_11(SF01_11, SD01_11,  HQF_11,   HQD_11,   HQF_11,   DCCN,      CH01_11, CV01_11, ODC_11)
ArcCell02_11 = FODOCELL_11(SF02_11, SD02_11,  HQF_11,   HQD_11,   HQF_11,   DCCN,      CH02_11, CV02_11, ODC_11)
ArcCell03_11 = FODOCELL_11(SF03_11, SD03_11,  HQF_11,   HQD_11,   HQF_11,   DCCN,      CH03_11, CV03_11, ODC_11)
ArcCell04_11 = FODOCELL_11(SF04_11, SD04_11,  HQF_11,   HQD_11,   HQF_11,   DCCN,      CH04_11, CV04_11, ODC_11)
ArcCell05_11 = FODOCELL_11(SF05_11, SD05_11,  HQF_11,   HQD_11,   HQF_11,   DCCN,      CH05_11, CV05_11, ODC_11)
ArcCell06_11 = FODOCELL_11(SF06_11, SD06_11,  HQF_11,   HQD_11,   HQF_11,   DCCN,      CH06_11, CV06_11, ODC_11)
ArcCell07_11 = FODOCELL_11(SF07_11, SD07_11,  HQF_11,   HQD_11,   HQF_11,   DCCN,      CH07_11, CV07_11, ODC_11)
ArcCell08_11 = FODOCELL_11(SF08_11, SD08_11,  HQF_11,   HQD_11,   HQF_11,   DCCN,      CH08_11, CV08_11, ODC_11)
ArcCell09_11 = FODOCELL_11(SF09_11, SD09_11,  HQF_11,   HQD_11,   HQF_11,   DCCN,      CH09_11, CV09_11, ODC_11)
ArcCell10_11 = FODOCELL_11(SF10_11, SD10_11,  HQF_11,   HQD_11,   HQF_11,   DCCN,      CH10_11, CV10_11, ODC_11)
ArcCell11_11 = FODOCELL_11(SF11_11, SD11_11,  HQF_11,   HQD_11,   HQF_11,   DCCN,      CH11_11, CV11_11, ODC_11)
ArcCell12_11 = FODOCELL_11(SF12_11, SD12_11,  HQF_11,   HQD_11,   HQF_11,   DCCN,      CH12_11, CV12_11, ODC_11)
ArcCell13_11 = FODOCELL_11(SF13_11, SD13_11,  HQF_11,   HQD_11,   HQF_11,   DCCN,      CH13_11, CV13_11, ODC_11)
ArcCell14_11 = FODOCELL_11(SF14_11, SD14_11,  HQF_11,   HQD_11,   HQF_11,   DCCN,      CH14_11, CV14_11, ODC_11)
ArcCell15_11 = FODOCELL_11(SF15_11, SD15_11,  HQF_11,   HQD_11,   HQF_11,   DCCN,      CH15_11, CV15_11, ODC_11)
ArcCell16_11 = FODOCELL_11(SF16_11, SD16_11,  HQF_11,   HQD_11,   HQF_11,   DCCN,      CH16_11, CV16_11, ODC_11)

ARC_1 = [MARC_BEG, ArcCell01_1..., ArcCell02_1..., ArcCell03_1..., ArcCell04_1...,
   ArcCell05_1..., ArcCell06_1..., ArcCell07_1..., ArcCell08_1..., ArcCell09_1..., ArcCell10_1...,
   ArcCell11_1..., ArcCell12_1..., ArcCell13_1..., ArcCell14_1..., ArcCell15_1..., ArcCell16_1...,
   MARC_END]
ARC_3 = [MARC_BEG, ArcCell01_3..., ArcCell02_3..., ArcCell03_3..., ArcCell04_3...,
   ArcCell05_3..., ArcCell06_3..., ArcCell07_3..., ArcCell08_3..., ArcCell09_3..., ArcCell10_3...,
   ArcCell11_3..., ArcCell12_3..., ArcCell13_3..., ArcCell14_3..., ArcCell15_3..., ArcCell16_3...,
   MARC_END]
ARC_5 = [MARC_BEG, ArcCell01_5..., ArcCell02_5..., ArcCell03_5..., ArcCell04_5...,
   ArcCell05_5..., ArcCell06_5..., ArcCell07_5..., ArcCell08_5..., ArcCell09_5..., ArcCell10_5...,
   ArcCell11_5..., ArcCell12_5..., ArcCell13_5..., ArcCell14_5..., ArcCell15_5..., ArcCell16_5...,
   MARC_END]
ARC_7 = [MARC_BEG, ArcCell01_7..., ArcCell02_7..., ArcCell03_7..., ArcCell04_7...,
   ArcCell05_7..., ArcCell06_7..., ArcCell07_7..., ArcCell08_7..., ArcCell09_7..., ArcCell10_7...,
   ArcCell11_7..., ArcCell12_7..., ArcCell13_7..., ArcCell14_7..., ArcCell15_7..., ArcCell16_7...,
   MARC_END]
ARC_9 = [MARC_BEG, ArcCell01_9..., ArcCell02_9..., ArcCell03_9..., ArcCell04_9...,
   ArcCell05_9..., ArcCell06_9..., ArcCell07_9..., ArcCell08_9..., ArcCell09_9..., ArcCell10_9...,
   ArcCell11_9..., ArcCell12_9..., ArcCell13_9..., ArcCell14_9..., ArcCell15_9..., ArcCell16_9...,
   MARC_END]
ARC_11 = [MARC_BEG, ArcCell01_11..., ArcCell02_11..., ArcCell03_11..., ArcCell04_11...,
   ArcCell05_11..., ArcCell06_11..., ArcCell07_11..., ArcCell08_11..., ArcCell09_11...,
   ArcCell10_11..., ArcCell11_11..., ArcCell12_11..., ArcCell13_11..., ArcCell14_11...,
   ArcCell15_11..., ArcCell16_11..., MARC_END]


UTIL_10 = [HQF_9, OHQC, CH17_9, ODBC_9, DB23_9, ODB23, DB23_9, ODBQ, HQD_9, HQD_9,
   OHQC, CV17_9, ODBC_9, DB23_9, ODB23, DB23_9, ODSX, SF17_9, OSSX, SF17_9, OQSX, HQF_9, HQF_9, ODF_9, HQM22_9, MQ22,
   HQM22_9, ODF_9, HQM21_9, MQ21, HQM21_9, ODF29, DB23_9, ODB23, DB23_9, ODF29, HQM20_9, MQ20,
   HQM20_9, ODF29, DB23_9, ODB23, DB23_9, ODF29, HQM19_9, MQ19, HQM19_9, ODF20_9, HQM18_9,
   MQ18, HQM18_9, ODF20_9, HQM17_9, MQ17, HQM17_9, ODF20_9, HQM16_9, MQ16,
   HQM16_9, ODF20_9, HQM15_9, MQ15, HQM15_9, ODF20_9, HQM14_9, MQ14, HQM14_9,
   ODF20_9, HQM13_9, MQ13, HQM13_9, ODF29, DB23_9, ODB23, DB23_9, ODF29, HQM12_9, MQ12,
   HQM12_9, ODF29, DB23_9, ODB23, DB23_9, OD09_0, HQFSS_10, MQSS_10, HQFSS_10, OD10S_0,
   HQDSS_10, MQSS_10, HQDSS_10, OD10S_0, HQFSS_10, MQSS_10, HQFSS_10, OD10S_0,
   HQDSS_10, MQSS_10, HQDSS_10, OD09H_0, HQFLSS_10, MQSS_10, HQFLSS_10,
   OD10RF_0, RF0, OD10RF_0, RF0, OD10RF_0, HQDLSS_10, MQSS_10, HQDLSS_10,
   OD10RF_0, RF0, OD10RF_0, RF0, OD10RF_0, HQFLSS_10, MQSS_10, HQFLSS_10,
   OD10RF_0, RF0, OD10RF_0, RF0, OD10RF_0, HQDLSS_10, MQSS_10, HQDLSS_10,
   OD10RF_0, RF0, OD10RF_0, RF0, OD10RF_0, HQFLSS_10, MQSS_10, HQFLSS_10,
   OD10RF_0, RF0, OD10IP_0]
UTIL_10N = [OD10IP_0, RF0, OD10RF_0, HQDLSS_10, MQSS_10, HQDLSS_10,
   OD10RF_0, RF0, OD10RF_0, RF0, OD10RF_0, HQFLSS_10, MQSS_10, HQFLSS_10,
   OD10RF_0, RF0, OD10RF_0, RF0, OD10RF_0, HQDLSS_10, MQSS_10, HQDLSS_10,
   OD10RF_0, RF0, OD10RF_0, RF0, OD10RF_0, HQFLSS_10, MQSS_10, HQFLSS_10,
   OD10RF_0, RF0, OD10RF_0, RF0, OD10RF_0, HQDLSS_10, MQSS_10, HQDLSS_10, OD10S_0,
   HQFSS_10, MQSS_10, HQFSS_10, OD10S_0, HQDSS_10, MQSS_10, HQDSS_10, OD10S_0,
   HQFSS_10, MQSS_10, HQFSS_10, OD10S_0, HQDSS_10, MQSS_10, HQDSS_10, OD10_0,
   DB23_10, ODB23, DB23_10, ODF29, HQM12_10, MQ12, HQM12_10, ODF29, DB23_10, ODB23, DB23_10, ODF29, HQM13_10,
   MQ13, HQM13_10, ODF20_10, HQM14_10, MQ14, HQM14_10, ODF20_10, HQM15_10, MQ15,
   HQM15_10, ODF20_10, HQM16_10, MQ16, HQM16_10, ODF20_10, HQM17_10, MQ17,
   HQM17_10, ODF20_10, HQM18_10, MQ18, HQM18_10, ODF20_10, HQM19_10, MQ19,
   HQM19_10, ODF20_10, HQM20_10, MQ20, HQM20_10, ODF29, DB23_10, ODB23, DB23_10, ODF29,
   HQM21_10, MQ21, HQM21_10, ODF29, DB23_10, ODB23, DB23_10, ODF29, HQM22_10, MQ22, HQM22_10,
   ODF_10, HQF_11, HQF_11, OQSX, SF00_11, OSSX, SF00_11, ODSX, DB23_10, ODB23, DB23_10, ODBC_10, CV00_11, OHQC, HQD_11,
   HQD_11, ODBQ, DB23_10, ODB23, DB23_10, ODBC_10, CH00_11, OHQC, HQF_11]
   UTIL_12 = [HQF_11, OQSX, SF17_11, OSSX, SF17_11, ODSX, DB23_11, ODB23, DB23_11, ODBC_11, CV17_11, OHQC, HQD_11, HQD_11,
   ODBQ, DB23_11, ODB23, DB23_11, ODBC_11, CH17_11, OHQC, HQF_11, HQF_11, OQSX, SF18_11, OSSX, SF18_11, 
   ODF_11, HQM22_11, MQ22, HQM22_11, ODF_11, HQM21_11, MQ21, HQM21_11, ODF29,
   DB23_11, ODB23, DB23_11, ODF29, HQM20_11, MQ20, HQM20_11, ODF29, DB23_11, ODB23, DB23_11, ODF29, HQM19_11,
   MQ19, HQM19_11, ODF20_11, HQM18_11, MQ18, HQM18_11, ODF20_11, HQM17_11, MQ17,
   HQM17_11, ODF20_11, HQM16_11, MQ16, HQM16_11, ODF20_11, HQM15_11, MQ15,
   HQM15_11, ODF20_11, HQM14_11, MQ14, HQM14_11, ODF20_11, HQM13_11, MQ13,
   HQM13_11, ODF29, DB23_11, ODB23, DB23_11, ODF29, HQM12_11, MQ12, HQM12_11, ODF29, DB23_11, ODB23, DB23_11,
   OD11_0, HQFSS_12, MQSS_12, HQFSS_12, OD12, HQDSS_12, MQSS_12, HQDSS_12, OD12,
   HQFSS_12, MQSS_12, HQFSS_12, OD12, HQDSS_12, MQSS_12, HQDSS_12, OD11_0,
   DB12M, OD12_4, DB12M, OD12_4, DB12M, OD11_U, HQFSS_12, MQSS_12, HQFSS_12, OD12, HQDSS_12, MQSS_12,
   HQDSS_12, OD12, HQFSS_12, MQSS_12, HQFSS_12, ODMID_12]
UTIL_12N = [ODMID_12N, HQDSS_12, MQSS_12, HQDSS_12, OD12, HQFSS_12, MQSS_12,
   HQFSS_12, OD12, HQDSS_12, MQSS_12, HQDSS_12, OD12, HQFSS_12, MQSS_12,
   HQFSS_12, OD11_M, DB12P, OD12_4, DB12P, OD12_4, DB12P, OD12_0, HQDSS_12, MQSS_12, HQDSS_12, MKICK_INJ, OD12,
   HQFSS_12, MQSS_12, HQFSS_12, OD12, HQDSS_12, MQSS_12, HQDSS_12, OD12, MCOLL_INJ,
   HQFSS_12, MQSS_12, HQFSS_12, OD12_0, DB23_12, ODB23, DB23_12, ODF29, HQM14_12, MQ14,
   HQM14_12, ODF29, DB23_12, ODB23, DB23_12, ODF29, HQM15_12, MQ15, HQM15_12, ODF20_12,
   HQM16_12, MQ16, HQM16_12, ODF20_12, HQM17_12, MQ17, HQM17_12, ODF20_12,
   HQM18_12, MQ18, HQM18_12, ODF20_12, HQM19_12, MQ19, HQM19_12, ODF29, DB23_12, ODB23, DB23_12,
   ODF29, HQM20_12, MQ20, HQM20_12, ODF29, DB23_12, ODB23, DB23_12, ODF29, HQM21_12, MQ21,
   HQM21_12, ODF_12, HQM22_12, MQ22, HQM22_12, ODF_12, SFM1_1, OSSX, SFM1_1, ODX17, HQF_1,
   HQF_1, OHQC, CH00_1, ODBC_1, DB23_12, ODB23, DB23_12, ODBQ, HQD_1, HQD_1, OHQC, CV00_1, ODBC_1,
   DB23_12, ODB23, DB23_12, ODSX, SF00_1, OSSX, SF00_1, OQSX, HQF_1]
UTIL_2 = [HQF_1, OHQC, CH17_1, ODBC_1, DB23_1, ODB23, DB23_1, ODBQ, HQD_1, HQD_1,
   OHQC, CV17_1, ODBC_1, DB23_1, ODB23, DB23_1, ODSX, SF17_1, OSSX, SF17_1, OQSX, HQF_1, HQF_1, ODF_1, HQM22_1, MQ22,
   HQM22_1, ODF_1, HQM21_1, MQ21, HQM21_1, ODF_1, HQM20_1, MQ20, HQM20_1, ODF_1,
   HQM19_1, MQ19, HQM19_1, ODF_1, HQM18_1, MQ18, HQM18_1, ODF29, DB23_1, ODB23, DB23_1, ODF29,
   HQM17_1, MQ17, HQM17_1, ODF29, DB23_1, ODB23, DB23_1, ODF29, HQM16_1, MQ16, HQM16_1,
   ODF20_1, HQM15_1, MQ15, HQM15_1, ODF20_1, HQM14_1, MQ14, HQM14_1, ODF20_1,
   HQM13_1, MQ13, HQM13_1, ODF29, DB23_1, ODB23, DB23_1, ODF29, HQM12_1, MQ12, HQM12_1, ODF29,
   DB23_1, ODB23, DB23_1, OLLFT0_2, HQDSS_2, MQSS_2, HQDSS_2, OLLFSX_2, SX41_2, ODX17, HQFSS_2,
   MQSS_2, HQFSS_2, OLLFSX_2, SX42_2, ODX17, HQDSS_2, MQSS_2, HQDSS_2, MCOLL_H1, OLLFSX_2,
   SX43_2, ODX17, HQFSS_2, MQSS_2, HQFSS_2, OLLFSX_2, MCOLL_H2, SX44_2, ODX17, HQDSS_2,
   MQSS_2, HQDSS_2, OLLFSX_2, SX45_2, ODX17, HQFSS_2, MQSS_2, HQFSS_2, OLLFSX_2,
   MCOLL_H3, SX46_2, ODX17, HQDSS_2, MQSS_2, HQDSS_2, OLLFCA_2]
UTIL_2N = [OLLFCB_2, HQFSS_2, MQSS_2, HQFSS_2, ODX17, SX47_2, OLLFSX_2,
   HQDSS_2, MQSS_2, HQDSS_2, ODX17, SX48_2, OLLFSX_2, HQFSS_2, MQSS_2, HQFSS_2,
   ODX17, SX49_2, OLLFSX_2, HQDSS_2, MQSS_2, HQDSS_2, ODX17, SX50_2, MLAMB, OLLFSX_2,
   HQFSS_2, MQSS_2, HQFSS_2, ODX17, SX51_2, OLLFSX_2, HQDSS_2, MQSS_2, HQDSS_2,
   ODX17, SX52_2, OLLFSX_2, HQFSS_2, MQSS_2, HQFSS_2, OLLFT0_2, DB23_2, ODB23, DB23_2, ODF29,
   HQM12_2, MQ12, HQM12_2, ODF29, DB23_2, ODB23, DB23_2, ODF29, HQM13_2, MQ13, HQM13_2,
   ODF20_2, HQM14_2, MQ14, HQM14_2, ODF20_2, HQM15_2, MQ15, HQM15_2, ODF20_2,
   HQM16_2, MQ16, HQM16_2, ODF29, DB23_2, ODB23, DB23_2, ODF29, HQM17_2, MQ17, HQM17_2, ODF29,
   DB23_2, ODB23, DB23_2, ODF29, HQM18_2, MQ18, HQM18_2, ODF_2, HQM19_2, MQ19, HQM19_2, ODF_2,
   HQM20_2, MQ20, HQM20_2, ODF_2, HQM21_2, MQ21, HQM21_2, ODF_2, HQM22_2, MQ22,
   HQM22_2, ODF_2, HQF_3, HQF_3, OQSX, SF00_3, OSSX, SF00_3, ODSX, DB23_2, ODB23, DB23_2, ODBC_3, CV00_3, OHQC, HQD_3,
   HQD_3, ODBQ, DB23_2, ODB23, DB23_2, ODBC_3, CH00_3, OHQC, HQF_3]
UTIL_4 = [HQF_3, OQSX, SF17_3, OSSX, SF17_3, ODSX, DB23_3, ODB23, DB23_3, ODBC_3, CV17_3, OHQC, HQD_3, HQD_3, ODBQ,
   DB23_3, ODB23, DB23_3, ODBC_3, CH17_3, OHQC, HQF_3, HQF_3, OQSX, SF18_3, OSSX, SF18_3, ODF_3, HQD17_3,
   MQ17, HQD17_3, ODF_3, HQF16_3, MQ16, HQF16_3, ODF29, DB23_3, ODB23, DB23_3, ODF29, HQD15_3,
   MQ15, HQD15_3, ODF29, DB23_3, ODB23, DB23_3, ODF29, HQF14_3, MQ14, HQF14_3, ODF20_3,
   HQD13_3, MQ13, HQD13_3, ODF20_3, HQF12_3, MQ12, HQF12_3, ODF20_3, HQD11_3,
   MQ11, HQD11_3, ODF29, DB23_3, ODB23, DB23_3, ODF29, HQF10_3, MQ10, HQF10_3, ODF29, 
   DB23_3, ODB23, DB23_3, OD3_3, HQD9_3,
   MQ9, HQD9_3, OD4_3, HQF8_3, MQ8, HQF8_3, OD4_3, HQD7_3, MQ7, HQD7_3,
   OD4_3, HQF6_3, MQ6, HQF6_3, OD4_3, HQD5_3, MQ5, HQD5_3, OD4_3,
   HQF4_3, MQ4, HQF4_3,
   OD4_3,
   HQD3_3, MQ3, HQD3_3,
   OD3_3,
   DB4P,
   ODB12,
   DB4P,
   ODB12,
   DB4P,
   ODB12,
   DB4P,
   ODB12,
   HQF2_3, MQ2, HQF2_3,
   OD2_3,
   HQD1_3, MQ1, HQD1_3,
   OD1_3 ] 
UTIL_4N = [ OD1_4, HQD1_4, MQ1, HQD1_4,
  OD2_4,
  HQF2_4, MQ2, HQF2_4,
  ODB12,
  DB4M,
  ODB12,
  DB4M,
  ODB12,
  DB4M,
  ODB12,
  DB4M,
  OD3_4,
  HQD3_4, MQ3, HQD3_4,
  OD4_4, HQF4_4, MQ4, HQF4_4, OD4_4, HQD5_4, MQ5, HQD5_4, OD4_4,
  HQF6_4, MQ6, HQF6_4, OD4_4, HQD7_4, MQ7, HQD7_4, OD5_4, DB23_4, ODB23, DB23_4, 
  ODF29, HQF8_4, MQ8, HQF8_4, ODF29, DB23_4, ODB23, DB23_4, ODF29, HQD9_4, MQ9,
  HQD9_4, ODF20_4, HQF10_4, MQ10, HQF10_4, ODF20_4, HQD11_4, MQ11, HQD11_4,
  ODF20_4, HQF12_4, MQ12, HQF12_4, ODF20_4, HQD13_4, MQ13, HQD13_4, ODF20_4,
  HQF14_4, MQ14, HQF14_4, ODF29, DB23_4, ODB23, DB23_4, ODF29, HQD15_4, MQ15, HQD15_4, ODF29,
  DB23_4, ODB23, DB23_4, ODF29, HQF16_4, MQ16, HQF16_4, ODF_4, HQD17_4, MQ17, HQD17_4, ODF_4,
  HQF18_4, MQ18, HQF18_4, ODF_4, SFM1_5, OSSX, SFM1_5, ODX17, HQF_5, HQF_5, OHQC, CH00_5, ODBC_5, DB23_4, ODB23, DB23_4,
  ODBQ, HQD_5, HQD_5, OHQC, CV00_5, ODBC_5,  DB23_4, ODB23, DB23_4, ODSX, SF00_5, OSSX, SF00_5, OQSX,
  HQF_5]

RING  = [
   IP6,   IR_6N...,     reverse(HALF_IR6N)..., ARC_7...,    HALF_IR8..., IR_8...,
   IP8,   IR_8N...,     reverse(HALF_IR8N)..., ARC_9...,    UTIL_10...,
   IP10,  UTIL_10N...,  ARC_11...,             UTIL_12...,
   IP12,  UTIL_12N...,  ARC_1...,              UTIL_2...,
   IP2,   UTIL_2N...,   ARC_3...,              UTIL_4...,
   IP4,   UTIL_4N...,   ARC_5...,              HALF_IR6..., IR_6...,     IP6]
# serialize("esr_main_rad.jls", RING)

