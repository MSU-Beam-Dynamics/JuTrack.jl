include("../JuTrack.jl")
using. JuTrack
using JSON
none = 0
k1er_6 = -0.2278853772
k2er_6 = 0.2201156485
lrd2er = 5.5
g_d2er = 0.0032977170393289
db2er_ang_6 = -2.0*asin(0.5*g_d2er*lrd2er)
k3er_6 = 0.1788259727
k4er_6 = -0.007637024497
k5er_6 = -0.1625838028
l01_str = 2.726
db4er_ang_6 = 0.0119246
db3er_ang_6 = (38.79e-3-db2er_ang_6-5*db4er_ang_6)
k6er_6 = 0.1774622637
k7er_6 = -0.1354227061
k8er_6 = 0.06050609079
k9er_6 = 0.1006385785
k10er_6 = -0.03463797301
k11er_6 = -0.1024445065
k12er_6 = 0.2282592047
k13er_6 = -0.2580290222
k14er_6 = 0.1768808343
k15er_6 = -0.2188168368
lsol20 = 6.2
a_electron = 0.0011596521869
phi2_6 = 1.570796327
ksol2_6 = phi2_6/(1.0+a_electron)/(2.0*lsol20)
lqls1 = 0.9819319
kqls7_6 = 0.448239
lqls2 = 1.8
kqls6_6 = -0.36553
lqls3 = 1.8
kqls5_6 = 0.251782
lqls4 = 0.5187944
kqls4_6 = 0.397172
kqls3_6 = 0.251782
kqls2_6 = -0.36553
kqls1_6 = 0.448239
kff6_6 = -0.1017285434
theta1 = 0.097805890449924
kff5_6 = 0.1945703088
kff4_6 = -0.2136492148
kff3_6 = 0.4082525179
kff2_6 = -0.2037488048
kff1_6 = 0.365384125
lsol5 = 2.5
phi1_6 = 0.0
ksol1_6 = phi1_6/(1.0+a_electron)/(2.0*lsol5)
lqss5 = 0.6861532
kqss5_6 = -0.536478377
lqss4 = 1.020723
kqss4_6 = 0.07662062675
lqss3 = 1.634532
kqss3_6 = 0.2853875397
lqss2 = 0.9550568
kqss2_6 = 0.0005099694535
lqss1 = 0.6480402
kqss1_6 = -0.3555522202
kd_6a = 0.05421254005
kf_6a = 0.2691115209
ksf01_7 = 3.120884289
lcv = 0.2
kd_7 = -0.3116231882
kd_6 = kd_7
kd_6b = kd_6
ksd01_7 = -4.294201621
lch = 0.2
kf_7 = 0.311806547
kf_6 = kf_7
kf_6b = kf_6
ksf02_7 = 3.120884289
kd_6c = kd_6
ksd02_7 = -4.294201621
kf_6c = kf_6
ksf03_7 = 3.120884289
ksd03_7 = -4.294201621
ksf04_7 = 3.120884289
ksd04_7 = -4.294201621
ksf05_7 = 3.120884289
ksd05_7 = -4.294201621
ksf06_7 = 3.120884289
ksd06_7 = -4.294201621
ksf07_7 = 3.120884289
ksd07_7 = -4.294201621
ksf08_7 = 3.120884289
ksd08_7 = -4.294201621
ksf09_7 = 3.120884289
ksd09_7 = -4.294201621
ksf10_7 = 3.120884289
ksd10_7 = -4.294201621
ksf11_7 = 3.120884289
ksd11_7 = -4.294201621
ksf12_7 = 3.120884289
ksd12_7 = -4.294201621
ksf13_7 = 3.120884289
ksd13_7 = -4.294201621
ksf14_7 = 3.120884289
ksd14_7 = -4.294201621
kf_7c = kf_7
ksf15_7 = 3.120884289
kd_7c = kd_7
ksd15_7 = -4.294201621
kf_7b = kf_7
ksf16_7 = 3.120884289
kd_7b = kd_7
ksd16_7 = -4.294201621
kf_7a = 0.003675365321
kd_7a = 0.2730363268
phi1_8 = phi1_6
ksol1_8 = phi1_8/(1.0+a_electron)/(2.0*lsol5)
kqss1_7 = -0.2230165281
kqss2_7 = -0.0582811745
kqss3_7 = -0.01193738829
kqss4_7 = 0.2137457745
kqss5_7 = -0.1844120558
kff1_7 = 0.102700806
kff2_7 = -0.1712617411
kff3_7 = 0.2272981814
kff4_7 = -0.1936954915
kff5_7 = 0.2451513988
kff6_7 = -0.08489389511
phi2_8 = phi2_6
ksol2_8 = phi2_8/(1.0+a_electron)/(2.0*lsol20)
kqls1_7 = 0.448239
kqls2_7 = -0.36553
kqls3_7 = 0.251782
kqls4_7 = 0.397172
kqls5_7 = 0.251782
kqls6_7 = -0.36553
kqls7_7 = 0.448239
k15ef_8 = 0.1235640461
db3ef_ang_8 = 0.0165327
k14ef_8 = -0.1542741001
k13ef_8 = 0.298202212
k12ef_8 = -0.2434214825
k11ef_8 = 0.3406620038
theta2 = 0.038785094488763
db1ef_ang_8 = -0.0014944764371143
db2ef_ang_8 = (theta2-(4.0*db3ef_ang_8+2.0*db1ef_ang_8))/2
k10ef_8 = -0.2708898178
k9ef_8 = 0.3083359637
k8ef_8 = -0.2831439398
k7ef_8 = 0.04393988393
k6ef_8 = -0.2442159358
k5ef_8 = 0.1496742106
k4ef_8 = -0.1397024576
k3ef_8 = 0.1463303236
k2ef_8 = -0.08190548605
l23_str = 0.8914
k1ef_8 = 0.1019244209
k0ef_8 = -0.2273064537
k1er_8 = -0.2245628613
k2er_8 = 0.2145025935
db2er_ang_8 = -2.0*asin(0.5*g_d2er*lrd2er)
k3er_8 = 0.01485241877
db3er_ang_8 = 0.00489854
k4er_8 = 0.1817930759
k5er_8 = -0.1384591674
k6er_8 = 0.05486095398
k7er_8 = 0.09864667925
k8er_8 = -0.1628913735
k9er_8 = 0.1079272595
k10er_8 = 0.1530487875
k11er_8 = -0.2224981329
db4er_ang_8 = 0.01338427324
k12er_8 = 0.1781177892
db5er_ang_8 = (theta2-(db4er_ang_8+db3er_ang_8+db2er_ang_8))/3.0
k13er_8 = 0.119497545
k14er_8 = -0.1280021585
k15er_8 = 0.02964357365
kqls7_8 = 0.448239
kqls6_8 = -0.36553
kqls5_8 = 0.251782
kqls4_8 = 0.397172
kqls3_8 = 0.251782
kqls2_8 = -0.36553
kqls1_8 = 0.448239
kff6_8 = -0.1328195162
kff5_8 = 0.2665763448
kff4_8 = -0.2618177996
kff3_8 = 0.4057378695
kff2_8 = -0.2495845973
kff1_8 = 0.4256160053
kqss5_8 = -0.5909929783
kqss4_8 = 0.1179833172
kqss3_8 = 0.2420319659
kqss2_8 = 0.05070009734
kqss1_8 = -0.3376839281
kd_8a = 0.04587795267
kf_8a = 0.1783200006
ksd01_9 = -4.294201621
kd_9 = -0.3144176426
kd_8 = kd_9
kd_8b = kd_8
ksf01_9 = 3.120884289
kf_9 = 0.314601847
kf_8 = kf_9
kf_8b = kf_8
ksd02_9 = -4.294201621
kd_8c = kd_8
ksf02_9 = 3.120884289
kf_8c = kf_8
ksd03_9 = -4.294201621
ksf03_9 = 3.120884289
ksd04_9 = -4.294201621
ksf04_9 = 3.120884289
ksd05_9 = -4.294201621
ksf05_9 = 3.120884289
ksd06_9 = -4.294201621
ksf06_9 = 3.120884289
ksd07_9 = -4.294201621
ksf07_9 = 3.120884289
ksd08_9 = -4.294201621
ksf08_9 = 3.120884289
ksd09_9 = -4.294201621
ksf09_9 = 3.120884289
ksd10_9 = -4.294201621
ksf10_9 = 3.120884289
ksd11_9 = -4.294201621
ksf11_9 = 3.120884289
ksd12_9 = -4.294201621
ksf12_9 = 3.120884289
ksd13_9 = -4.294201621
ksf13_9 = 3.120884289
ksd14_9 = -4.294201621
ksf14_9 = 3.120884289
ksd15_9 = -4.294201621
ksf15_9 = 3.120884289
ksd16_9 = -4.294201621
ksf16_9 = 3.120884289
rot_angle = 0.011382582078224
ksf17_9 = 0.0
km22_9 = -0.1185051943
km21_9 = -0.1323371828
km20_9 = 0.2082883818
km19_9 = -0.1797478437
km18_9 = 0.05902104803
km17_9 = 0.03551808764
km16_9 = -0.09344124638
km15_9 = 0.1424311997
km14_9 = -0.145625108
km13_9 = 0.2169778696
km12_9 = -0.1221746876
kfss_10 = 0.1559175595
kdss_10 = -0.1566717693
kfssl_10 = 0.1360968888
volt_rf = 0.0
kdssl_10 = -0.1170915216
km12_10 = 0.1312382732
km13_10 = -0.2444138469
km14_10 = 0.2775709987
km15_10 = -0.335725079
km16_10 = 0.235462357
km17_10 = -0.1130247262
km18_10 = -0.262909096
km19_10 = 0.297798994
km20_10 = -0.2581690321
km21_10 = 0.1651338464
km22_10 = -0.2809701903
kf_11 = 0.3137143598
ksf00_11 = 0.0
kd_11 = -0.313530424
ksf01_11 = 3.120884289
ksd01_11 = -4.294201621
ksf02_11 = 3.120884289
ksd02_11 = -4.294201621
ksf03_11 = 3.120884289
ksd03_11 = -4.294201621
ksf04_11 = 3.120884289
ksd04_11 = -4.294201621
ksf05_11 = 3.120884289
ksd05_11 = -4.294201621
ksf06_11 = 3.120884289
ksd06_11 = -4.294201621
ksf07_11 = 3.120884289
ksd07_11 = -4.294201621
ksf08_11 = 3.120884289
ksd08_11 = -4.294201621
ksf09_11 = 3.120884289
ksd09_11 = -4.294201621
ksf10_11 = 3.120884289
ksd10_11 = -4.294201621
ksf11_11 = 3.120884289
ksd11_11 = -4.294201621
ksf12_11 = 3.120884289
ksd12_11 = -4.294201621
ksf13_11 = 3.120884289
ksd13_11 = -4.294201621
ksf14_11 = 3.120884289
ksd14_11 = -4.294201621
ksf15_11 = 3.120884289
ksd15_11 = -4.294201621
ksf16_11 = 3.120884289
ksd16_11 = -4.294201621
ksf17_11 = 0.0
ksf18_11 = 0.0
km22_11 = -0.3359731738
km21_11 = 0.1816358361
km20_11 = -0.1485365183
km19_11 = 0.2154916461
km18_11 = -0.1500435633
km17_11 = 0.003291249802
km16_11 = -0.002162170431
km15_11 = 0.004247583095
km14_11 = -0.1327720571
km13_11 = 0.2102817626
km12_11 = -0.09862623398
kfss_12 = 0.1462599289
kdss_12 = -0.1483172335
km14_12 = -0.1682943092
km15_12 = 0.2954917251
km16_12 = -0.3359748537
km17_12 = 0.1079293418
km18_12 = -0.3359676021
km19_12 = 0.2857040226
km20_12 = -0.2708760985
km21_12 = 0.1908516877
km22_12 = -0.3359850684
ksfm1_1 = 0.0
kf_1 = 0.3113930481
kd_1 = -0.3112098116
ksf00_1 = 0.0
ksd01_1 = -4.294201621
ksf01_1 = 3.120884289
ksd02_1 = -4.294201621
ksf02_1 = 3.120884289
ksd03_1 = -4.294201621
ksf03_1 = 3.120884289
ksd04_1 = -4.294201621
ksf04_1 = 3.120884289
ksd05_1 = -4.294201621
ksf05_1 = 3.120884289
ksd06_1 = -4.294201621
ksf06_1 = 3.120884289
ksd07_1 = -4.294201621
ksf07_1 = 3.120884289
ksd08_1 = -4.294201621
ksf08_1 = 3.120884289
ksd09_1 = -4.294201621
ksf09_1 = 3.120884289
ksd10_1 = -4.294201621
ksf10_1 = 3.120884289
ksd11_1 = -4.294201621
ksf11_1 = 3.120884289
ksd12_1 = -4.294201621
ksf12_1 = 3.120884289
ksd13_1 = -4.294201621
ksf13_1 = 3.120884289
ksd14_1 = -4.294201621
ksf14_1 = 3.120884289
ksd15_1 = -4.294201621
ksf15_1 = 3.120884289
ksd16_1 = -4.294201621
ksf16_1 = 3.120884289
ksf17_1 = 0.0
km22_1 = -0.02701644447
km21_1 = 0.01425207743
km20_1 = 0.2030463977
km19_1 = -0.2811627068
km18_1 = 0.03017571658
km17_1 = 0.1073519579
km16_1 = 0.141150677
km15_1 = -0.1756395994
km14_1 = 0.2831036567
km13_1 = -0.1648119071
km12_1 = 0.140091896
kdss_2 = -0.075634279
s41_2 = 0.0
kfss_2 = 0.06204145252
s42_2 = 0.0
s43_2 = 0.0
s44_2 = 0.0
s45_2 = 0.0
s46_2 = 0.0
s47_2 = 0.0
s48_2 = 0.0
s49_2 = 0.0
s50_2 = 0.0
s51_2 = 0.0
s52_2 = 0.0
km12_2 = 0.05881784823
km13_2 = -0.1665854015
km14_2 = 0.09504800846
km15_2 = 0.1864067227
km16_2 = -0.239794824
km17_2 = 0.1274757467
km18_2 = 0.2829300672
km19_2 = -0.2570712752
km20_2 = -0.1853975692
km21_2 = 0.335817511
km22_2 = -0.0628929574
kf_3 = 0.3113944689
ksf00_3 = 0.0
kd_3 = -0.3112112319
ksf01_3 = 3.120884289
ksd01_3 = -4.294201621
ksf02_3 = 3.120884289
ksd02_3 = -4.294201621
ksf03_3 = 3.120884289
ksd03_3 = -4.294201621
ksf04_3 = 3.120884289
ksd04_3 = -4.294201621
ksf05_3 = 3.120884289
ksd05_3 = -4.294201621
ksf06_3 = 3.120884289
ksd06_3 = -4.294201621
ksf07_3 = 3.120884289
ksd07_3 = -4.294201621
ksf08_3 = 3.120884289
ksd08_3 = -4.294201621
ksf09_3 = 3.120884289
ksd09_3 = -4.294201621
ksf10_3 = 3.120884289
ksd10_3 = -4.294201621
ksf11_3 = 3.120884289
ksd11_3 = -4.294201621
ksf12_3 = 3.120884289
ksd12_3 = -4.294201621
ksf13_3 = 3.120884289
ksd13_3 = -4.294201621
ksf14_3 = 3.120884289
ksd14_3 = -4.294201621
ksf15_3 = 3.120884289
ksd15_3 = -4.294201621
ksf16_3 = 3.120884289
ksd16_3 = -4.294201621
ksf17_3 = 0.0
ksf18_3 = 0.0
kd17_3 = -0.2877048243
kf16_3 = 0.1824243838
kd15_3 = 0.008578802751
kf14_3 = -0.1480150886
kd13_3 = 0.07877982097
kf12_3 = 0.05174533641
kd11_3 = 0.02794082844
kf10_3 = -0.06070609411
kd9_3 = -0.06354608493
kf8_3 = 0.0794935845
kd7_3 = 0.05008828685
kf6_3 = -0.057826039
kd5_3 = -0.05081367491
kf4_3 = 0.05382316298
kd3_3 = 0.03358142199
kf2_3 = -0.1105862703
kd1_3 = 0.07357421792
kd1_4 = 0.04845280893
kf2_4 = -0.1041863504
kd3_4 = 0.04110636443
kf4_4 = 0.04838351833
kd5_4 = -0.05314628092
kf6_4 = -0.03914277884
kd7_4 = 0.07175145111
kf8_4 = -0.0003306038141
kd9_4 = -0.08791343748
kf10_4 = 0.06520333211
kd11_4 = 0.0933824065
kf12_4 = -0.07155545708
kd13_4 = -0.08356348977
kf14_4 = 0.1861642754
kd15_4 = -0.1150529702
kf16_4 = 0.09127406468
kd17_4 = -0.2013945475
kf18_4 = 0.0
ksfm1_5 = 0.0
kf_5 = 0.3139724652
kd_5 = -0.3137884525
ksf00_5 = 0.0
ksd01_5 = -4.294201621
ksf01_5 = 3.120884289
ksd02_5 = -4.294201621
ksf02_5 = 3.120884289
ksd03_5 = -4.294201621
ksf03_5 = 3.120884289
ksd04_5 = -4.294201621
ksf04_5 = 3.120884289
ksd05_5 = -4.294201621
ksf05_5 = 3.120884289
ksd06_5 = -4.294201621
ksf06_5 = 3.120884289
ksd07_5 = -4.294201621
ksf07_5 = 3.120884289
ksd08_5 = -4.294201621
ksf08_5 = 3.120884289
ksd09_5 = -4.294201621
ksf09_5 = 3.120884289
ksd10_5 = -4.294201621
ksf10_5 = 3.120884289
ksd11_5 = -4.294201621
ksf11_5 = 3.120884289
ksd12_5 = -4.294201621
ksf12_5 = 3.120884289
ksd13_5 = -4.294201621
ksf13_5 = 3.120884289
ksd14_5 = -4.294201621
ksf14_5 = 3.120884289
kf_5c = kf_5
ksd15_5 = -4.294201621
kd_5c = kd_5
ksf15_5 = 3.120884289
kf_5b = kf_5
ksd16_5 = -4.294201621
kd_5b = kd_5
ksf16_5 = 3.120884289
kf_5a = 0.5299652224
kd_5a = -0.2006935196
kqss1_5 = -0.4634519914
kqss2_5 = 0.03305713451
kqss3_5 = 0.3758491433
kqss4_5 = 0.0809412749
kqss5_5 = -0.5222748076
kff1_5 = 0.3443756116
kff2_5 = -0.2380503976
kff3_5 = 0.3775796003
kff4_5 = -0.1949412053
kff5_5 = 0.187251159
kff6_5 = -0.08489389511
kqls1_5 = 0.448239
kqls2_5 = -0.36553
kqls3_5 = 0.251782
kqls4_5 = 0.397172
kqls5_5 = 0.251782
kqls6_5 = -0.36553
kqls7_5 = 0.448239
k14ef_6 = -0.2208787114
db1ef_ang_6 = -0.0015
db2ef_ang_6 = -0.0123486
db3ef_ang_6 = 0.013
db4ef_ang_6 = 0.0015
db5ef_ang_6 = (theta2-db1ef_ang_6-db2ef_ang_6*2-db3ef_ang_6-db4ef_ang_6*2)/4
k13ef_6 = 0.1615485113
k12ef_6 = -0.1996130663
k11ef_6 = 0.2264235945
k10ef_6 = -0.1803467941
k9ef_6 = 0.345229068
k8ef_6 = -0.1990558344
k7ef_6 = 0.195769614
k6ef_6 = -0.08410656585
k5ef_6 = -0.004101903809
k4ef_6 = -0.06677557789
k3ef_6 = 0.1367818102
k2ef_6 = -0.07457307348
k1ef_6 = 0.1000859995
k0ef_6 = -0.218192315

# Define elements
ip6 = MARKER(name="ip6")
q1er_6 = KQUAD(name="q1er_6", len=1.8, k1=k1er_6)
q2er_6 = KQUAD(name="q2er_6", len=1.4, k1=k2er_6)
d2er_6 = RBEND(name="d2er_6", len=lrd2er, angle=db2er_ang_6)
sq = KQUAD(name="sq", len=0.25, k1=0.0)
sq3er_6 = KQUAD(name="sq3er_6", len=0.25, k1=0.0)
qir06 = KQUAD(name="qir06", len=0.6, k1=0.0)
q3er_6 = KQUAD(name="q3er_6", len=0.6, k1=k3er_6)
sq4er_6 = KQUAD(name="sq4er_6", len=0.25, k1=0.0)
q4er_6 = KQUAD(name="q4er_6", len=0.6, k1=k4er_6)
sq5er_6 = KQUAD(name="sq5er_6", len=0.25, k1=0.0)
q5er_6 = KQUAD(name="q5er_6", len=0.6, k1=k5er_6)
d3er_6 = RBEND(name="d3er_6", len=l01_str, angle=db3er_ang_6)
sq6er_6 = KQUAD(name="sq6er_6", len=0.25, k1=0.0)
qir12 = KQUAD(name="qir12", len=1.2, k1=0.0)
q6er_6 = KQUAD(name="q6er_6", len=1.2, k1=k6er_6)
sq7er_6 = KQUAD(name="sq7er_6", len=0.25, k1=0.0)
q7er_6 = KQUAD(name="q7er_6", len=1.2, k1=k7er_6)
sq8er_6 = KQUAD(name="sq8er_6", len=0.25, k1=0.0)
q8er_6 = KQUAD(name="q8er_6", len=1.2, k1=k8er_6)
rf_crab = RFCA(name="rf_crab", len=4.0, volt=0.0, freq=394e6)
sq9er_6 = KQUAD(name="sq9er_6", len=0.25, k1=0.0)
q9er_6 = KQUAD(name="q9er_6", len=1.2, k1=k9er_6)
sq10er_6 = KQUAD(name="sq10er_6", len=0.25, k1=0.0)
q10er_6 = KQUAD(name="q10er_6", len=1.2, k1=k10er_6)
d4er_6 = RBEND(name="d4er_6", len=l01_str, angle=db4er_ang_6)
q11er_6 = KQUAD(name="q11er_6", len=1.2, k1=k11er_6)
q12er_6 = KQUAD(name="q12er_6", len=1.2, k1=k12er_6)
q13er_6 = KQUAD(name="q13er_6", len=1.2, k1=k13er_6)
q14er_6 = KQUAD(name="q14er_6", len=1.2, k1=k14er_6)
q15er_6 = KQUAD(name="q15er_6", len=1.2, k1=k15er_6)
mlrr_6 = MARKER(name="mlrr_6")
mrot4 = MARKER(name="mrot4")
hsol20_6 = SOLENOID(name="hsol20_6", len=lsol20 / 2, ks=ksol2_6)
hqls7_6 = KQUAD(name="hqls7_6", len=lqls1 / 2.0, k1=kqls7_6)
hqls6_6 = KQUAD(name="hqls6_6", len=lqls2 / 2.0, k1=kqls6_6)
hqls5_6 = KQUAD(name="hqls5_6", len=lqls3 / 2.0, k1=kqls5_6)
hqls4_6 = KQUAD(name="hqls4_6", len=lqls4 / 2.0, k1=kqls4_6)
hqls3_6 = KQUAD(name="hqls3_6", len=lqls3 / 2.0, k1=kqls3_6)
hqls2_6 = KQUAD(name="hqls2_6", len=lqls2 / 2.0, k1=kqls2_6)
hqls1_6 = KQUAD(name="hqls1_6", len=lqls1 / 2.0, k1=kqls1_6)
mrot3 = MARKER(name="mrot3")
qff6_6 = KQUAD(name="qff6_6", len=0.5, k1=kff6_6)
db23_6 = RBEND(name="db23_6", len=3.8, angle=theta1 / 5)
qff5_6 = KQUAD(name="qff5_6", len=0.5, k1=kff5_6)
qff4_6 = KQUAD(name="qff4_6", len=0.5, k1=kff4_6)
qff3_6 = KQUAD(name="qff3_6", len=0.5, k1=kff3_6)
qff2_6 = KQUAD(name="qff2_6", len=0.5, k1=kff2_6)
qff1_6 = KQUAD(name="qff1_6", len=0.5, k1=kff1_6)
mrot2 = MARKER(name="mrot2")
hsol5_6 = SOLENOID(name="hsol5_6", len=lsol5 / 2, ks=ksol1_6)
hqss5_6 = KQUAD(name="hqss5_6", len=lqss5 / 2.0, k1=kqss5_6)
hqss4_6 = KQUAD(name="hqss4_6", len=lqss4 / 2.0, k1=kqss4_6)
hqss3_6 = KQUAD(name="hqss3_6", len=lqss3 / 2.0, k1=kqss3_6)
hqss2_6 = KQUAD(name="hqss2_6", len=lqss2 / 2.0, k1=kqss2_6)
hqss1_6 = KQUAD(name="hqss1_6", len=lqss1 / 2.0, k1=kqss1_6)
mrot1 = MARKER(name="mrot1")
hqd_6a = KQUAD(name="hqd_6a", len=0.25, k1=kd_6a)
hqf_6a = KQUAD(name="hqf_6a", len=0.25, k1=kf_6a)
marc_beg = MARKER(name="marc_beg")
sf01_7 = KSEXT(name="sf01_7", len=0.24, k2=ksf01_7)
edge1_002 = thinMULTIPOLE(name="edge1_002", PolynomB=[ 0.0,  -5.4058569363803e-05, 0.0, 0.0])
d01a_002 = SBEND(name="d01a_002", len=2.7260903830638, angle=0.011254024541609)
edge2_002 = thinMULTIPOLE(name="edge2_002", PolynomB=[ 0.0,  7.5958885437618e-06, 0.0, 0.0])
edge3_002 = thinMULTIPOLE(name="edge3_002", PolynomB=[ 0.0,  -7.5958885437618e-06, 0.0, 0.0])
d23_002 = SBEND(name="d23_002", len=0.89140050297051, angle=0.003679937833007)
d01b_002 = SBEND(name="d01b_002", len=2.7260903830638, angle=0.011254024541609)
cv01_7 = VKICKER(name="cv01_7", len=lcv)
hqd_6b = KQUAD(name="hqd_6b", len=0.25, k1=kd_6b)
sd01_7 = KSEXT(name="sd01_7", len=0.57, k2=ksd01_7)
ch01_7 = HKICKER(name="ch01_7", len=lch)
hqf_6b = KQUAD(name="hqf_6b", len=0.25, k1=kf_6b)
sf02_7 = KSEXT(name="sf02_7", len=0.24, k2=ksf02_7)
cv02_7 = VKICKER(name="cv02_7", len=lcv)
hqd_6c = KQUAD(name="hqd_6c", len=0.25, k1=kd_6c)
sd02_7 = KSEXT(name="sd02_7", len=0.57, k2=ksd02_7)
ch02_7 = HKICKER(name="ch02_7", len=lch)
hqf_6c = KQUAD(name="hqf_6c", len=0.25, k1=kf_6c)
sf03_7 = KSEXT(name="sf03_7", len=0.24, k2=ksf03_7)
edge1_000 = thinMULTIPOLE(name="edge1_000", PolynomB=[ 0.0,  -4.6116670427111e-05, 0.0, 0.0])
d01a_000 = SBEND(name="d01a_000", len=2.7260771047243, angle=0.010394537621096)
edge2_000 = thinMULTIPOLE(name="edge2_000", PolynomB=[ 0.0,  6.4800348490203e-06, 0.0, 0.0])
edge3_000 = thinMULTIPOLE(name="edge3_000", PolynomB=[ 0.0,  -6.4800348490203e-06, 0.0, 0.0])
d23_000 = SBEND(name="d23_000", len=0.89140042908298, angle=0.0033989116740339)
d01b_000 = SBEND(name="d01b_000", len=2.7260771047243, angle=0.010394537621096)
cv03_7 = VKICKER(name="cv03_7", len=lcv)
hqd_7 = KQUAD(name="hqd_7", len=0.25, k1=kd_7)
sd03_7 = KSEXT(name="sd03_7", len=0.57, k2=ksd03_7)
ch03_7 = HKICKER(name="ch03_7", len=lch)
hqf_7 = KQUAD(name="hqf_7", len=0.25, k1=kf_7)
sf04_7 = KSEXT(name="sf04_7", len=0.24, k2=ksf04_7)
cv04_7 = VKICKER(name="cv04_7", len=lcv)
sd04_7 = KSEXT(name="sd04_7", len=0.57, k2=ksd04_7)
ch04_7 = HKICKER(name="ch04_7", len=lch)
sf05_7 = KSEXT(name="sf05_7", len=0.24, k2=ksf05_7)
cv05_7 = VKICKER(name="cv05_7", len=lcv)
sd05_7 = KSEXT(name="sd05_7", len=0.57, k2=ksd05_7)
ch05_7 = HKICKER(name="ch05_7", len=lch)
sf06_7 = KSEXT(name="sf06_7", len=0.24, k2=ksf06_7)
cv06_7 = VKICKER(name="cv06_7", len=lcv)
sd06_7 = KSEXT(name="sd06_7", len=0.57, k2=ksd06_7)
ch06_7 = HKICKER(name="ch06_7", len=lch)
sf07_7 = KSEXT(name="sf07_7", len=0.24, k2=ksf07_7)
cv07_7 = VKICKER(name="cv07_7", len=lcv)
sd07_7 = KSEXT(name="sd07_7", len=0.57, k2=ksd07_7)
ch07_7 = HKICKER(name="ch07_7", len=lch)
sf08_7 = KSEXT(name="sf08_7", len=0.24, k2=ksf08_7)
cv08_7 = VKICKER(name="cv08_7", len=lcv)
sd08_7 = KSEXT(name="sd08_7", len=0.57, k2=ksd08_7)
ch08_7 = HKICKER(name="ch08_7", len=lch)
sf09_7 = KSEXT(name="sf09_7", len=0.24, k2=ksf09_7)
cv09_7 = VKICKER(name="cv09_7", len=lcv)
sd09_7 = KSEXT(name="sd09_7", len=0.57, k2=ksd09_7)
ch09_7 = HKICKER(name="ch09_7", len=lch)
sf10_7 = KSEXT(name="sf10_7", len=0.24, k2=ksf10_7)
cv10_7 = VKICKER(name="cv10_7", len=lcv)
sd10_7 = KSEXT(name="sd10_7", len=0.57, k2=ksd10_7)
ch10_7 = HKICKER(name="ch10_7", len=lch)
sf11_7 = KSEXT(name="sf11_7", len=0.24, k2=ksf11_7)
cv11_7 = VKICKER(name="cv11_7", len=lcv)
sd11_7 = KSEXT(name="sd11_7", len=0.57, k2=ksd11_7)
ch11_7 = HKICKER(name="ch11_7", len=lch)
sf12_7 = KSEXT(name="sf12_7", len=0.24, k2=ksf12_7)
cv12_7 = VKICKER(name="cv12_7", len=lcv)
sd12_7 = KSEXT(name="sd12_7", len=0.57, k2=ksd12_7)
ch12_7 = HKICKER(name="ch12_7", len=lch)
sf13_7 = KSEXT(name="sf13_7", len=0.24, k2=ksf13_7)
cv13_7 = VKICKER(name="cv13_7", len=lcv)
sd13_7 = KSEXT(name="sd13_7", len=0.57, k2=ksd13_7)
ch13_7 = HKICKER(name="ch13_7", len=lch)
sf14_7 = KSEXT(name="sf14_7", len=0.24, k2=ksf14_7)
cv14_7 = VKICKER(name="cv14_7", len=lcv)
sd14_7 = KSEXT(name="sd14_7", len=0.57, k2=ksd14_7)
ch14_7 = HKICKER(name="ch14_7", len=lch)
hqf_7c = KQUAD(name="hqf_7c", len=0.25, k1=kf_7c)
sf15_7 = KSEXT(name="sf15_7", len=0.24, k2=ksf15_7)
edge1_003 = thinMULTIPOLE(name="edge1_003", PolynomB=[ 0.0,  -5.7199379023574e-05, 0.0, 0.0])
d01a_003 = SBEND(name="d01a_003", len=2.7260956342846, angle=0.011576332594406)
edge2_003 = thinMULTIPOLE(name="edge2_003", PolynomB=[ 0.0,  8.0371713790994e-06, 0.0, 0.0])
edge3_003 = thinMULTIPOLE(name="edge3_003", PolynomB=[ 0.0,  -8.0371713790994e-06, 0.0, 0.0])
d23_003 = SBEND(name="d23_003", len=0.89140053219056, angle=0.0037853217274147)
d01b_003 = SBEND(name="d01b_003", len=2.7260956342846, angle=0.011576332594406)
cv15_7 = VKICKER(name="cv15_7", len=lcv)
hqd_7c = KQUAD(name="hqd_7c", len=0.25, k1=kd_7c)
sd15_7 = KSEXT(name="sd15_7", len=0.57, k2=ksd15_7)
ch15_7 = HKICKER(name="ch15_7", len=lch)
hqf_7b = KQUAD(name="hqf_7b", len=0.25, k1=kf_7b)
sf16_7 = KSEXT(name="sf16_7", len=0.24, k2=ksf16_7)
cv16_7 = VKICKER(name="cv16_7", len=lcv)
hqd_7b = KQUAD(name="hqd_7b", len=0.25, k1=kd_7b)
sd16_7 = KSEXT(name="sd16_7", len=0.57, k2=ksd16_7)
ch16_7 = HKICKER(name="ch16_7", len=lch)
hqf_7a = KQUAD(name="hqf_7a", len=0.25, k1=kf_7a)
marc_end = MARKER(name="marc_end")
hqd_7a = KQUAD(name="hqd_7a", len=0.25, k1=kd_7a)
hsol5_8 = SOLENOID(name="hsol5_8", len=lsol5 / 2, ks=ksol1_8)
hqss1_7 = KQUAD(name="hqss1_7", len=0.3240201, k1=kqss1_7)
hqss2_7 = KQUAD(name="hqss2_7", len=0.4775284, k1=kqss2_7)
hqss3_7 = KQUAD(name="hqss3_7", len=0.817266, k1=kqss3_7)
hqss4_7 = KQUAD(name="hqss4_7", len=0.5103615, k1=kqss4_7)
hqss5_7 = KQUAD(name="hqss5_7", len=0.3430766, k1=kqss5_7)
qff1_7 = KQUAD(name="qff1_7", len=0.5, k1=kff1_7)
db23_7 = RBEND(name="db23_7", len=3.8, angle=theta1 / 5)
qff2_7 = KQUAD(name="qff2_7", len=0.5, k1=kff2_7)
qff3_7 = KQUAD(name="qff3_7", len=0.5, k1=kff3_7)
qff4_7 = KQUAD(name="qff4_7", len=0.5, k1=kff4_7)
qff5_7 = KQUAD(name="qff5_7", len=0.5, k1=kff5_7)
qff6_7 = KQUAD(name="qff6_7", len=0.5, k1=kff6_7)
hsol20_8 = SOLENOID(name="hsol20_8", len=lsol20 / 2, ks=ksol2_8)
hqls1_7 = KQUAD(name="hqls1_7", len=0.49096595, k1=kqls1_7)
hqls2_7 = KQUAD(name="hqls2_7", len=0.9, k1=kqls2_7)
hqls3_7 = KQUAD(name="hqls3_7", len=0.9, k1=kqls3_7)
hqls4_7 = KQUAD(name="hqls4_7", len=0.2593972, k1=kqls4_7)
hqls5_7 = KQUAD(name="hqls5_7", len=0.9, k1=kqls5_7)
hqls6_7 = KQUAD(name="hqls6_7", len=0.9, k1=kqls6_7)
hqls7_7 = KQUAD(name="hqls7_7", len=0.49096595, k1=kqls7_7)
mlrf_8 = MARKER(name="mlrf_8")
q15ef_8 = KQUAD(name="q15ef_8", len=1.2, k1=k15ef_8)
d3ef_8 = RBEND(name="d3ef_8", len=3.8, angle=db3ef_ang_8)
q14ef_8 = KQUAD(name="q14ef_8", len=1.2, k1=k14ef_8)
q13ef_8 = KQUAD(name="q13ef_8", len=1.2, k1=k13ef_8)
q12ef_8 = KQUAD(name="q12ef_8", len=1.2, k1=k12ef_8)
q11ef_8 = KQUAD(name="q11ef_8", len=1.2, k1=k11ef_8)
d2ef_8 = RBEND(name="d2ef_8", len=l01_str, angle=db2ef_ang_8)
q10ef_8 = KQUAD(name="q10ef_8", len=1.2, k1=k10ef_8)
q9ef_8 = KQUAD(name="q9ef_8", len=1.2, k1=k9ef_8)
q8ef_8 = KQUAD(name="q8ef_8", len=1.2, k1=k8ef_8)
q7ef_8 = KQUAD(name="q7ef_8", len=1.2, k1=k7ef_8)
q6ef_8 = KQUAD(name="q6ef_8", len=1.2, k1=k6ef_8)
q5ef_8 = KQUAD(name="q5ef_8", len=1.2, k1=k5ef_8)
q4ef_8 = KQUAD(name="q4ef_8", len=1.2, k1=k4ef_8)
q3ef_8 = KQUAD(name="q3ef_8", len=0.6, k1=k3ef_8)
q2ef_8 = KQUAD(name="q2ef_8", len=0.6, k1=k2ef_8)
d1ef_8 = RBEND(name="d1ef_8", len=l23_str, angle=db1ef_ang_8)
q1ef_8 = KQUAD(name="q1ef_8", len=1.61, k1=k1ef_8)
q0ef_8 = KQUAD(name="q0ef_8", len=1.2, k1=k0ef_8)
ip8 = MARKER(name="ip8")
q1er_8 = KQUAD(name="q1er_8", len=1.8, k1=k1er_8)
q2er_8 = KQUAD(name="q2er_8", len=1.4, k1=k2er_8)
d2er_8 = RBEND(name="d2er_8", len=lrd2er, angle=db2er_ang_8)
q3er_8 = KQUAD(name="q3er_8", len=0.6, k1=k3er_8)
d3er_8 = RBEND(name="d3er_8", len=l01_str, angle=db3er_ang_8)
q4er_8 = KQUAD(name="q4er_8", len=0.6, k1=k4er_8)
q5er_8 = KQUAD(name="q5er_8", len=1.2, k1=k5er_8)
q6er_8 = KQUAD(name="q6er_8", len=1.2, k1=k6er_8)
q7er_8 = KQUAD(name="q7er_8", len=1.2, k1=k7er_8)
q8er_8 = KQUAD(name="q8er_8", len=1.2, k1=k8er_8)
q9er_8 = KQUAD(name="q9er_8", len=1.2, k1=k9er_8)
q10er_8 = KQUAD(name="q10er_8", len=1.2, k1=k10er_8)
q11er_8 = KQUAD(name="q11er_8", len=1.2, k1=k11er_8)
d4er_8 = RBEND(name="d4er_8", len=l01_str, angle=db4er_ang_8)
q12er_8 = KQUAD(name="q12er_8", len=1.2, k1=k12er_8)
d5er_8 = RBEND(name="d5er_8", len=l01_str, angle=db5er_ang_8)
q13er_8 = KQUAD(name="q13er_8", len=1.2, k1=k13er_8)
q14er_8 = KQUAD(name="q14er_8", len=1.2, k1=k14er_8)
q15er_8 = KQUAD(name="q15er_8", len=1.2, k1=k15er_8)
mlrr_8 = MARKER(name="mlrr_8")
hqls7_8 = KQUAD(name="hqls7_8", len=0.49096595, k1=kqls7_8)
hqls6_8 = KQUAD(name="hqls6_8", len=0.9, k1=kqls6_8)
hqls5_8 = KQUAD(name="hqls5_8", len=0.9, k1=kqls5_8)
hqls4_8 = KQUAD(name="hqls4_8", len=0.2593972, k1=kqls4_8)
hqls3_8 = KQUAD(name="hqls3_8", len=0.9, k1=kqls3_8)
hqls2_8 = KQUAD(name="hqls2_8", len=0.9, k1=kqls2_8)
hqls1_8 = KQUAD(name="hqls1_8", len=0.49096595, k1=kqls1_8)
qff6_8 = KQUAD(name="qff6_8", len=0.5, k1=kff6_8)
db23_8 = RBEND(name="db23_8", len=3.8, angle=theta1 / 5)
qff5_8 = KQUAD(name="qff5_8", len=0.5, k1=kff5_8)
qff4_8 = KQUAD(name="qff4_8", len=0.5, k1=kff4_8)
qff3_8 = KQUAD(name="qff3_8", len=0.5, k1=kff3_8)
qff2_8 = KQUAD(name="qff2_8", len=0.5, k1=kff2_8)
qff1_8 = KQUAD(name="qff1_8", len=0.5, k1=kff1_8)
hqss5_8 = KQUAD(name="hqss5_8", len=0.3430766, k1=kqss5_8)
hqss4_8 = KQUAD(name="hqss4_8", len=0.5103615, k1=kqss4_8)
hqss3_8 = KQUAD(name="hqss3_8", len=0.817266, k1=kqss3_8)
hqss2_8 = KQUAD(name="hqss2_8", len=0.4775284, k1=kqss2_8)
hqss1_8 = KQUAD(name="hqss1_8", len=0.3240201, k1=kqss1_8)
hqd_8a = KQUAD(name="hqd_8a", len=0.25, k1=kd_8a)
hqf_8a = KQUAD(name="hqf_8a", len=0.25, k1=kf_8a)
ch01_9 = HKICKER(name="ch01_9", len=lch)
edge1_004 = thinMULTIPOLE(name="edge1_004", PolynomB=[ 0.0,  -3.6226320616364e-05, 0.0, 0.0])
d01a_004 = SBEND(name="d01a_004", len=2.7260605686558, angle=0.0092127457973839)
edge2_004 = thinMULTIPOLE(name="edge2_004", PolynomB=[ 0.0,  5.0903811383167e-06, 0.0, 0.0])
edge3_004 = thinMULTIPOLE(name="edge3_004", PolynomB=[ 0.0,  -5.0903811383167e-06, 0.0, 0.0])
d23_004 = SBEND(name="d23_004", len=0.89140033706547, angle=0.0030124953214579)
d01b_004 = SBEND(name="d01b_004", len=2.7260605686558, angle=0.0092127457973839)
sd01_9 = KSEXT(name="sd01_9", len=0.57, k2=ksd01_9)
hqd_8b = KQUAD(name="hqd_8b", len=0.25, k1=kd_8b)
cv01_9 = VKICKER(name="cv01_9", len=lcv)
sf01_9 = KSEXT(name="sf01_9", len=0.24, k2=ksf01_9)
hqf_8b = KQUAD(name="hqf_8b", len=0.25, k1=kf_8b)
ch02_9 = HKICKER(name="ch02_9", len=lch)
sd02_9 = KSEXT(name="sd02_9", len=0.57, k2=ksd02_9)
hqd_8c = KQUAD(name="hqd_8c", len=0.25, k1=kd_8c)
cv02_9 = VKICKER(name="cv02_9", len=lcv)
sf02_9 = KSEXT(name="sf02_9", len=0.24, k2=ksf02_9)
hqf_8c = KQUAD(name="hqf_8c", len=0.25, k1=kf_8c)
ch03_9 = HKICKER(name="ch03_9", len=lch)
sd03_9 = KSEXT(name="sd03_9", len=0.57, k2=ksd03_9)
hqd_9 = KQUAD(name="hqd_9", len=0.25, k1=kd_9)
cv03_9 = VKICKER(name="cv03_9", len=lcv)
sf03_9 = KSEXT(name="sf03_9", len=0.24, k2=ksf03_9)
hqf_9 = KQUAD(name="hqf_9", len=0.25, k1=kf_9)
ch04_9 = HKICKER(name="ch04_9", len=lch)
sd04_9 = KSEXT(name="sd04_9", len=0.57, k2=ksd04_9)
cv04_9 = VKICKER(name="cv04_9", len=lcv)
sf04_9 = KSEXT(name="sf04_9", len=0.24, k2=ksf04_9)
ch05_9 = HKICKER(name="ch05_9", len=lch)
sd05_9 = KSEXT(name="sd05_9", len=0.57, k2=ksd05_9)
cv05_9 = VKICKER(name="cv05_9", len=lcv)
sf05_9 = KSEXT(name="sf05_9", len=0.24, k2=ksf05_9)
ch06_9 = HKICKER(name="ch06_9", len=lch)
sd06_9 = KSEXT(name="sd06_9", len=0.57, k2=ksd06_9)
cv06_9 = VKICKER(name="cv06_9", len=lcv)
sf06_9 = KSEXT(name="sf06_9", len=0.24, k2=ksf06_9)
ch07_9 = HKICKER(name="ch07_9", len=lch)
sd07_9 = KSEXT(name="sd07_9", len=0.57, k2=ksd07_9)
cv07_9 = VKICKER(name="cv07_9", len=lcv)
sf07_9 = KSEXT(name="sf07_9", len=0.24, k2=ksf07_9)
ch08_9 = HKICKER(name="ch08_9", len=lch)
sd08_9 = KSEXT(name="sd08_9", len=0.57, k2=ksd08_9)
cv08_9 = VKICKER(name="cv08_9", len=lcv)
sf08_9 = KSEXT(name="sf08_9", len=0.24, k2=ksf08_9)
ch09_9 = HKICKER(name="ch09_9", len=lch)
sd09_9 = KSEXT(name="sd09_9", len=0.57, k2=ksd09_9)
cv09_9 = VKICKER(name="cv09_9", len=lcv)
sf09_9 = KSEXT(name="sf09_9", len=0.24, k2=ksf09_9)
ch10_9 = HKICKER(name="ch10_9", len=lch)
sd10_9 = KSEXT(name="sd10_9", len=0.57, k2=ksd10_9)
cv10_9 = VKICKER(name="cv10_9", len=lcv)
sf10_9 = KSEXT(name="sf10_9", len=0.24, k2=ksf10_9)
ch11_9 = HKICKER(name="ch11_9", len=lch)
sd11_9 = KSEXT(name="sd11_9", len=0.57, k2=ksd11_9)
cv11_9 = VKICKER(name="cv11_9", len=lcv)
sf11_9 = KSEXT(name="sf11_9", len=0.24, k2=ksf11_9)
ch12_9 = HKICKER(name="ch12_9", len=lch)
sd12_9 = KSEXT(name="sd12_9", len=0.57, k2=ksd12_9)
cv12_9 = VKICKER(name="cv12_9", len=lcv)
sf12_9 = KSEXT(name="sf12_9", len=0.24, k2=ksf12_9)
ch13_9 = HKICKER(name="ch13_9", len=lch)
sd13_9 = KSEXT(name="sd13_9", len=0.57, k2=ksd13_9)
cv13_9 = VKICKER(name="cv13_9", len=lcv)
sf13_9 = KSEXT(name="sf13_9", len=0.24, k2=ksf13_9)
ch14_9 = HKICKER(name="ch14_9", len=lch)
sd14_9 = KSEXT(name="sd14_9", len=0.57, k2=ksd14_9)
cv14_9 = VKICKER(name="cv14_9", len=lcv)
sf14_9 = KSEXT(name="sf14_9", len=0.24, k2=ksf14_9)
ch15_9 = HKICKER(name="ch15_9", len=lch)
sd15_9 = KSEXT(name="sd15_9", len=0.57, k2=ksd15_9)
cv15_9 = VKICKER(name="cv15_9", len=lcv)
sf15_9 = KSEXT(name="sf15_9", len=0.24, k2=ksf15_9)
ch16_9 = HKICKER(name="ch16_9", len=lch)
sd16_9 = KSEXT(name="sd16_9", len=0.57, k2=ksd16_9)
cv16_9 = VKICKER(name="cv16_9", len=lcv)
sf16_9 = KSEXT(name="sf16_9", len=0.24, k2=ksf16_9)
ch17_9 = HKICKER(name="ch17_9", len=lch)
db23_9 = RBEND(name="db23_9", len=2.726, angle=rot_angle)
cv17_9 = VKICKER(name="cv17_9", len=lcv)
sf17_9 = KSEXT(name="sf17_9", len=0.24, k2=ksf17_9)
hqm22_9 = KQUAD(name="hqm22_9", len=0.4, k1=km22_9)
mq22 = MARKER(name="mq22")
hqm21_9 = KQUAD(name="hqm21_9", len=0.4, k1=km21_9)
mq21 = MARKER(name="mq21")
hqm20_9 = KQUAD(name="hqm20_9", len=0.4, k1=km20_9)
mq20 = MARKER(name="mq20")
hqm19_9 = KQUAD(name="hqm19_9", len=0.4, k1=km19_9)
mq19 = MARKER(name="mq19")
hqm18_9 = KQUAD(name="hqm18_9", len=0.4, k1=km18_9)
mq18 = MARKER(name="mq18")
hqm17_9 = KQUAD(name="hqm17_9", len=0.4, k1=km17_9)
mq17 = MARKER(name="mq17")
hqm16_9 = KQUAD(name="hqm16_9", len=0.4, k1=km16_9)
mq16 = MARKER(name="mq16")
hqm15_9 = KQUAD(name="hqm15_9", len=0.4, k1=km15_9)
mq15 = MARKER(name="mq15")
hqm14_9 = KQUAD(name="hqm14_9", len=0.4, k1=km14_9)
mq14 = MARKER(name="mq14")
hqm13_9 = KQUAD(name="hqm13_9", len=0.4, k1=km13_9)
mq13 = MARKER(name="mq13")
hqm12_9 = KQUAD(name="hqm12_9", len=0.4, k1=km12_9)
mq12 = MARKER(name="mq12")
hqfss_10 = KQUAD(name="hqfss_10", len=0.4, k1=kfss_10)
mqss_10 = MARKER(name="mqss_10")
hqdss_10 = KQUAD(name="hqdss_10", len=0.4, k1=kdss_10)
hqflss_10 = KQUAD(name="hqflss_10", len=0.6, k1=kfssl_10)
rf0 = RFCA(name="rf0", len=4.01667, volt=volt_rf, freq=0.0, lag=0.5, h=7560)
hqdlss_10 = KQUAD(name="hqdlss_10", len=0.6, k1=kdssl_10)
ip10 = MARKER(name="ip10")
db23_10 = RBEND(name="db23_10", len=2.726, angle=rot_angle)
hqm12_10 = KQUAD(name="hqm12_10", len=0.4, k1=km12_10)
hqm13_10 = KQUAD(name="hqm13_10", len=0.4, k1=km13_10)
hqm14_10 = KQUAD(name="hqm14_10", len=0.4, k1=km14_10)
hqm15_10 = KQUAD(name="hqm15_10", len=0.4, k1=km15_10)
hqm16_10 = KQUAD(name="hqm16_10", len=0.4, k1=km16_10)
hqm17_10 = KQUAD(name="hqm17_10", len=0.4, k1=km17_10)
hqm18_10 = KQUAD(name="hqm18_10", len=0.4, k1=km18_10)
hqm19_10 = KQUAD(name="hqm19_10", len=0.4, k1=km19_10)
hqm20_10 = KQUAD(name="hqm20_10", len=0.4, k1=km20_10)
hqm21_10 = KQUAD(name="hqm21_10", len=0.4, k1=km21_10)
hqm22_10 = KQUAD(name="hqm22_10", len=0.4, k1=km22_10)
hqf_11 = KQUAD(name="hqf_11", len=0.25, k1=kf_11)
sf00_11 = KSEXT(name="sf00_11", len=0.24, k2=ksf00_11)
cv00_11 = VKICKER(name="cv00_11", len=lcv)
hqd_11 = KQUAD(name="hqd_11", len=0.25, k1=kd_11)
ch00_11 = HKICKER(name="ch00_11", len=lch)
sf01_11 = KSEXT(name="sf01_11", len=0.24, k2=ksf01_11)
cv01_11 = VKICKER(name="cv01_11", len=lcv)
sd01_11 = KSEXT(name="sd01_11", len=0.24, k2=ksd01_11)
ch01_11 = HKICKER(name="ch01_11", len=lch)
sf02_11 = KSEXT(name="sf02_11", len=0.24, k2=ksf02_11)
cv02_11 = VKICKER(name="cv02_11", len=lcv)
sd02_11 = KSEXT(name="sd02_11", len=0.24, k2=ksd02_11)
ch02_11 = HKICKER(name="ch02_11", len=lch)
sf03_11 = KSEXT(name="sf03_11", len=0.24, k2=ksf03_11)
cv03_11 = VKICKER(name="cv03_11", len=lcv)
sd03_11 = KSEXT(name="sd03_11", len=0.24, k2=ksd03_11)
ch03_11 = HKICKER(name="ch03_11", len=lch)
sf04_11 = KSEXT(name="sf04_11", len=0.24, k2=ksf04_11)
cv04_11 = VKICKER(name="cv04_11", len=lcv)
sd04_11 = KSEXT(name="sd04_11", len=0.24, k2=ksd04_11)
ch04_11 = HKICKER(name="ch04_11", len=lch)
sf05_11 = KSEXT(name="sf05_11", len=0.24, k2=ksf05_11)
cv05_11 = VKICKER(name="cv05_11", len=lcv)
sd05_11 = KSEXT(name="sd05_11", len=0.24, k2=ksd05_11)
ch05_11 = HKICKER(name="ch05_11", len=lch)
sf06_11 = KSEXT(name="sf06_11", len=0.24, k2=ksf06_11)
cv06_11 = VKICKER(name="cv06_11", len=lcv)
sd06_11 = KSEXT(name="sd06_11", len=0.24, k2=ksd06_11)
ch06_11 = HKICKER(name="ch06_11", len=lch)
sf07_11 = KSEXT(name="sf07_11", len=0.24, k2=ksf07_11)
cv07_11 = VKICKER(name="cv07_11", len=lcv)
sd07_11 = KSEXT(name="sd07_11", len=0.24, k2=ksd07_11)
ch07_11 = HKICKER(name="ch07_11", len=lch)
sf08_11 = KSEXT(name="sf08_11", len=0.24, k2=ksf08_11)
cv08_11 = VKICKER(name="cv08_11", len=lcv)
sd08_11 = KSEXT(name="sd08_11", len=0.24, k2=ksd08_11)
ch08_11 = HKICKER(name="ch08_11", len=lch)
sf09_11 = KSEXT(name="sf09_11", len=0.24, k2=ksf09_11)
cv09_11 = VKICKER(name="cv09_11", len=lcv)
sd09_11 = KSEXT(name="sd09_11", len=0.24, k2=ksd09_11)
ch09_11 = HKICKER(name="ch09_11", len=lch)
sf10_11 = KSEXT(name="sf10_11", len=0.24, k2=ksf10_11)
cv10_11 = VKICKER(name="cv10_11", len=lcv)
sd10_11 = KSEXT(name="sd10_11", len=0.24, k2=ksd10_11)
ch10_11 = HKICKER(name="ch10_11", len=lch)
sf11_11 = KSEXT(name="sf11_11", len=0.24, k2=ksf11_11)
cv11_11 = VKICKER(name="cv11_11", len=lcv)
sd11_11 = KSEXT(name="sd11_11", len=0.24, k2=ksd11_11)
ch11_11 = HKICKER(name="ch11_11", len=lch)
sf12_11 = KSEXT(name="sf12_11", len=0.24, k2=ksf12_11)
cv12_11 = VKICKER(name="cv12_11", len=lcv)
sd12_11 = KSEXT(name="sd12_11", len=0.24, k2=ksd12_11)
ch12_11 = HKICKER(name="ch12_11", len=lch)
sf13_11 = KSEXT(name="sf13_11", len=0.24, k2=ksf13_11)
cv13_11 = VKICKER(name="cv13_11", len=lcv)
sd13_11 = KSEXT(name="sd13_11", len=0.24, k2=ksd13_11)
ch13_11 = HKICKER(name="ch13_11", len=lch)
sf14_11 = KSEXT(name="sf14_11", len=0.24, k2=ksf14_11)
cv14_11 = VKICKER(name="cv14_11", len=lcv)
sd14_11 = KSEXT(name="sd14_11", len=0.24, k2=ksd14_11)
ch14_11 = HKICKER(name="ch14_11", len=lch)
sf15_11 = KSEXT(name="sf15_11", len=0.24, k2=ksf15_11)
cv15_11 = VKICKER(name="cv15_11", len=lcv)
sd15_11 = KSEXT(name="sd15_11", len=0.24, k2=ksd15_11)
ch15_11 = HKICKER(name="ch15_11", len=lch)
sf16_11 = KSEXT(name="sf16_11", len=0.24, k2=ksf16_11)
cv16_11 = VKICKER(name="cv16_11", len=lcv)
sd16_11 = KSEXT(name="sd16_11", len=0.24, k2=ksd16_11)
ch16_11 = HKICKER(name="ch16_11", len=lch)
sf17_11 = KSEXT(name="sf17_11", len=0.24, k2=ksf17_11)
db23_11 = RBEND(name="db23_11", len=2.726, angle=rot_angle)
cv17_11 = VKICKER(name="cv17_11", len=lcv)
ch17_11 = HKICKER(name="ch17_11", len=lch)
sf18_11 = KSEXT(name="sf18_11", len=0.24, k2=ksf18_11)
hqm22_11 = KQUAD(name="hqm22_11", len=0.4, k1=km22_11)
hqm21_11 = KQUAD(name="hqm21_11", len=0.4, k1=km21_11)
hqm20_11 = KQUAD(name="hqm20_11", len=0.4, k1=km20_11)
hqm19_11 = KQUAD(name="hqm19_11", len=0.4, k1=km19_11)
hqm18_11 = KQUAD(name="hqm18_11", len=0.4, k1=km18_11)
hqm17_11 = KQUAD(name="hqm17_11", len=0.4, k1=km17_11)
hqm16_11 = KQUAD(name="hqm16_11", len=0.4, k1=km16_11)
hqm15_11 = KQUAD(name="hqm15_11", len=0.4, k1=km15_11)
hqm14_11 = KQUAD(name="hqm14_11", len=0.4, k1=km14_11)
hqm13_11 = KQUAD(name="hqm13_11", len=0.4, k1=km13_11)
hqm12_11 = KQUAD(name="hqm12_11", len=0.4, k1=km12_11)
hqfss_12 = KQUAD(name="hqfss_12", len=0.4, k1=kfss_12)
mqss_12 = MARKER(name="mqss_12")
hqdss_12 = KQUAD(name="hqdss_12", len=0.4, k1=kdss_12)
db12m = RBEND(name="db12m", len=2.726, angle=-0.0109083)
ip12 = MARKER(name="ip12")
db12p = RBEND(name="db12p", len=2.726, angle=0.0109083)
mkick_inj = MARKER(name="mkick_inj")
mcoll_inj = MARKER(name="mcoll_inj")
db23_12 = RBEND(name="db23_12", len=2.726, angle=rot_angle)
hqm14_12 = KQUAD(name="hqm14_12", len=0.4, k1=km14_12)
hqm15_12 = KQUAD(name="hqm15_12", len=0.4, k1=km15_12)
hqm16_12 = KQUAD(name="hqm16_12", len=0.4, k1=km16_12)
hqm17_12 = KQUAD(name="hqm17_12", len=0.4, k1=km17_12)
hqm18_12 = KQUAD(name="hqm18_12", len=0.4, k1=km18_12)
hqm19_12 = KQUAD(name="hqm19_12", len=0.4, k1=km19_12)
hqm20_12 = KQUAD(name="hqm20_12", len=0.4, k1=km20_12)
hqm21_12 = KQUAD(name="hqm21_12", len=0.4, k1=km21_12)
hqm22_12 = KQUAD(name="hqm22_12", len=0.4, k1=km22_12)
sfm1_1 = KSEXT(name="sfm1_1", len=0.24, k2=ksfm1_1)
hqf_1 = KQUAD(name="hqf_1", len=0.25, k1=kf_1)
ch00_1 = HKICKER(name="ch00_1", len=lch)
hqd_1 = KQUAD(name="hqd_1", len=0.25, k1=kd_1)
cv00_1 = VKICKER(name="cv00_1", len=lcv)
sf00_1 = KSEXT(name="sf00_1", len=0.24, k2=ksf00_1)
ch01_1 = HKICKER(name="ch01_1", len=lch)
sd01_1 = KSEXT(name="sd01_1", len=0.24, k2=ksd01_1)
cv01_1 = VKICKER(name="cv01_1", len=lcv)
sf01_1 = KSEXT(name="sf01_1", len=0.24, k2=ksf01_1)
ch02_1 = HKICKER(name="ch02_1", len=lch)
sd02_1 = KSEXT(name="sd02_1", len=0.24, k2=ksd02_1)
cv02_1 = VKICKER(name="cv02_1", len=lcv)
sf02_1 = KSEXT(name="sf02_1", len=0.24, k2=ksf02_1)
ch03_1 = HKICKER(name="ch03_1", len=lch)
sd03_1 = KSEXT(name="sd03_1", len=0.24, k2=ksd03_1)
cv03_1 = VKICKER(name="cv03_1", len=lcv)
sf03_1 = KSEXT(name="sf03_1", len=0.24, k2=ksf03_1)
ch04_1 = HKICKER(name="ch04_1", len=lch)
sd04_1 = KSEXT(name="sd04_1", len=0.24, k2=ksd04_1)
cv04_1 = VKICKER(name="cv04_1", len=lcv)
sf04_1 = KSEXT(name="sf04_1", len=0.24, k2=ksf04_1)
ch05_1 = HKICKER(name="ch05_1", len=lch)
sd05_1 = KSEXT(name="sd05_1", len=0.24, k2=ksd05_1)
cv05_1 = VKICKER(name="cv05_1", len=lcv)
sf05_1 = KSEXT(name="sf05_1", len=0.24, k2=ksf05_1)
ch06_1 = HKICKER(name="ch06_1", len=lch)
sd06_1 = KSEXT(name="sd06_1", len=0.24, k2=ksd06_1)
cv06_1 = VKICKER(name="cv06_1", len=lcv)
sf06_1 = KSEXT(name="sf06_1", len=0.24, k2=ksf06_1)
ch07_1 = HKICKER(name="ch07_1", len=lch)
sd07_1 = KSEXT(name="sd07_1", len=0.24, k2=ksd07_1)
cv07_1 = VKICKER(name="cv07_1", len=lcv)
sf07_1 = KSEXT(name="sf07_1", len=0.24, k2=ksf07_1)
ch08_1 = HKICKER(name="ch08_1", len=lch)
sd08_1 = KSEXT(name="sd08_1", len=0.24, k2=ksd08_1)
cv08_1 = VKICKER(name="cv08_1", len=lcv)
sf08_1 = KSEXT(name="sf08_1", len=0.24, k2=ksf08_1)
ch09_1 = HKICKER(name="ch09_1", len=lch)
sd09_1 = KSEXT(name="sd09_1", len=0.24, k2=ksd09_1)
cv09_1 = VKICKER(name="cv09_1", len=lcv)
sf09_1 = KSEXT(name="sf09_1", len=0.24, k2=ksf09_1)
ch10_1 = HKICKER(name="ch10_1", len=lch)
sd10_1 = KSEXT(name="sd10_1", len=0.24, k2=ksd10_1)
cv10_1 = VKICKER(name="cv10_1", len=lcv)
sf10_1 = KSEXT(name="sf10_1", len=0.24, k2=ksf10_1)
ch11_1 = HKICKER(name="ch11_1", len=lch)
sd11_1 = KSEXT(name="sd11_1", len=0.24, k2=ksd11_1)
cv11_1 = VKICKER(name="cv11_1", len=lcv)
sf11_1 = KSEXT(name="sf11_1", len=0.24, k2=ksf11_1)
ch12_1 = HKICKER(name="ch12_1", len=lch)
sd12_1 = KSEXT(name="sd12_1", len=0.24, k2=ksd12_1)
cv12_1 = VKICKER(name="cv12_1", len=lcv)
sf12_1 = KSEXT(name="sf12_1", len=0.24, k2=ksf12_1)
ch13_1 = HKICKER(name="ch13_1", len=lch)
sd13_1 = KSEXT(name="sd13_1", len=0.24, k2=ksd13_1)
cv13_1 = VKICKER(name="cv13_1", len=lcv)
sf13_1 = KSEXT(name="sf13_1", len=0.24, k2=ksf13_1)
ch14_1 = HKICKER(name="ch14_1", len=lch)
sd14_1 = KSEXT(name="sd14_1", len=0.24, k2=ksd14_1)
cv14_1 = VKICKER(name="cv14_1", len=lcv)
sf14_1 = KSEXT(name="sf14_1", len=0.24, k2=ksf14_1)
ch15_1 = HKICKER(name="ch15_1", len=lch)
sd15_1 = KSEXT(name="sd15_1", len=0.24, k2=ksd15_1)
cv15_1 = VKICKER(name="cv15_1", len=lcv)
sf15_1 = KSEXT(name="sf15_1", len=0.24, k2=ksf15_1)
ch16_1 = HKICKER(name="ch16_1", len=lch)
sd16_1 = KSEXT(name="sd16_1", len=0.24, k2=ksd16_1)
cv16_1 = VKICKER(name="cv16_1", len=lcv)
sf16_1 = KSEXT(name="sf16_1", len=0.24, k2=ksf16_1)
ch17_1 = HKICKER(name="ch17_1", len=lch)
db23_1 = RBEND(name="db23_1", len=2.726, angle=rot_angle)
cv17_1 = VKICKER(name="cv17_1", len=lcv)
sf17_1 = KSEXT(name="sf17_1", len=0.24, k2=ksf17_1)
hqm22_1 = KQUAD(name="hqm22_1", len=0.4, k1=km22_1)
hqm21_1 = KQUAD(name="hqm21_1", len=0.4, k1=km21_1)
hqm20_1 = KQUAD(name="hqm20_1", len=0.4, k1=km20_1)
hqm19_1 = KQUAD(name="hqm19_1", len=0.4, k1=km19_1)
hqm18_1 = KQUAD(name="hqm18_1", len=0.4, k1=km18_1)
hqm17_1 = KQUAD(name="hqm17_1", len=0.4, k1=km17_1)
hqm16_1 = KQUAD(name="hqm16_1", len=0.4, k1=km16_1)
hqm15_1 = KQUAD(name="hqm15_1", len=0.4, k1=km15_1)
hqm14_1 = KQUAD(name="hqm14_1", len=0.4, k1=km14_1)
hqm13_1 = KQUAD(name="hqm13_1", len=0.4, k1=km13_1)
hqm12_1 = KQUAD(name="hqm12_1", len=0.4, k1=km12_1)
hqdss_2 = KQUAD(name="hqdss_2", len=0.4, k1=kdss_2)
mqss_2 = MARKER(name="mqss_2")
sx41_2 = KSEXT(name="sx41_2", len=0.24, k2=s41_2)
hqfss_2 = KQUAD(name="hqfss_2", len=0.4, k1=kfss_2)
sx42_2 = KSEXT(name="sx42_2", len=0.24, k2=s42_2)
mcoll_h1 = MARKER(name="mcoll_h1")
sx43_2 = KSEXT(name="sx43_2", len=0.24, k2=s43_2)
mcoll_h2 = MARKER(name="mcoll_h2")
sx44_2 = KSEXT(name="sx44_2", len=0.24, k2=s44_2)
sx45_2 = KSEXT(name="sx45_2", len=0.24, k2=s45_2)
mcoll_h3 = MARKER(name="mcoll_h3")
sx46_2 = KSEXT(name="sx46_2", len=0.24, k2=s46_2)
ip2 = MARKER(name="ip2")
sx47_2 = KSEXT(name="sx47_2", len=0.24, k2=s47_2)
sx48_2 = KSEXT(name="sx48_2", len=0.24, k2=s48_2)
sx49_2 = KSEXT(name="sx49_2", len=0.24, k2=s49_2)
sx50_2 = KSEXT(name="sx50_2", len=0.24, k2=s50_2)
mlamb = MARKER(name="mlamb")
sx51_2 = KSEXT(name="sx51_2", len=0.24, k2=s51_2)
sx52_2 = KSEXT(name="sx52_2", len=0.24, k2=s52_2)
db23_2 = RBEND(name="db23_2", len=2.726, angle=rot_angle)
hqm12_2 = KQUAD(name="hqm12_2", len=0.4, k1=km12_2)
hqm13_2 = KQUAD(name="hqm13_2", len=0.4, k1=km13_2)
hqm14_2 = KQUAD(name="hqm14_2", len=0.4, k1=km14_2)
hqm15_2 = KQUAD(name="hqm15_2", len=0.4, k1=km15_2)
hqm16_2 = KQUAD(name="hqm16_2", len=0.4, k1=km16_2)
hqm17_2 = KQUAD(name="hqm17_2", len=0.4, k1=km17_2)
hqm18_2 = KQUAD(name="hqm18_2", len=0.4, k1=km18_2)
hqm19_2 = KQUAD(name="hqm19_2", len=0.4, k1=km19_2)
hqm20_2 = KQUAD(name="hqm20_2", len=0.4, k1=km20_2)
hqm21_2 = KQUAD(name="hqm21_2", len=0.4, k1=km21_2)
hqm22_2 = KQUAD(name="hqm22_2", len=0.4, k1=km22_2)
hqf_3 = KQUAD(name="hqf_3", len=0.25, k1=kf_3)
sf00_3 = KSEXT(name="sf00_3", len=0.24, k2=ksf00_3)
cv00_3 = HKICKER(name="cv00_3", len=lcv)
hqd_3 = KQUAD(name="hqd_3", len=0.25, k1=kd_3)
ch00_3 = HKICKER(name="ch00_3", len=lch)
sf01_3 = KSEXT(name="sf01_3", len=0.24, k2=ksf01_3)
cv01_3 = VKICKER(name="cv01_3", len=lcv)
sd01_3 = KSEXT(name="sd01_3", len=0.24, k2=ksd01_3)
ch01_3 = HKICKER(name="ch01_3", len=lch)
sf02_3 = KSEXT(name="sf02_3", len=0.24, k2=ksf02_3)
cv02_3 = VKICKER(name="cv02_3", len=lcv)
sd02_3 = KSEXT(name="sd02_3", len=0.24, k2=ksd02_3)
ch02_3 = HKICKER(name="ch02_3", len=lch)
sf03_3 = KSEXT(name="sf03_3", len=0.24, k2=ksf03_3)
cv03_3 = VKICKER(name="cv03_3", len=lcv)
sd03_3 = KSEXT(name="sd03_3", len=0.24, k2=ksd03_3)
ch03_3 = HKICKER(name="ch03_3", len=lch)
sf04_3 = KSEXT(name="sf04_3", len=0.24, k2=ksf04_3)
cv04_3 = VKICKER(name="cv04_3", len=lcv)
sd04_3 = KSEXT(name="sd04_3", len=0.24, k2=ksd04_3)
ch04_3 = HKICKER(name="ch04_3", len=lch)
sf05_3 = KSEXT(name="sf05_3", len=0.24, k2=ksf05_3)
cv05_3 = VKICKER(name="cv05_3", len=lcv)
sd05_3 = KSEXT(name="sd05_3", len=0.24, k2=ksd05_3)
ch05_3 = HKICKER(name="ch05_3", len=lch)
sf06_3 = KSEXT(name="sf06_3", len=0.24, k2=ksf06_3)
cv06_3 = VKICKER(name="cv06_3", len=lcv)
sd06_3 = KSEXT(name="sd06_3", len=0.24, k2=ksd06_3)
ch06_3 = HKICKER(name="ch06_3", len=lch)
sf07_3 = KSEXT(name="sf07_3", len=0.24, k2=ksf07_3)
cv07_3 = VKICKER(name="cv07_3", len=lcv)
sd07_3 = KSEXT(name="sd07_3", len=0.24, k2=ksd07_3)
ch07_3 = HKICKER(name="ch07_3", len=lch)
sf08_3 = KSEXT(name="sf08_3", len=0.24, k2=ksf08_3)
cv08_3 = VKICKER(name="cv08_3", len=lcv)
sd08_3 = KSEXT(name="sd08_3", len=0.24, k2=ksd08_3)
ch08_3 = HKICKER(name="ch08_3", len=lch)
sf09_3 = KSEXT(name="sf09_3", len=0.24, k2=ksf09_3)
cv09_3 = VKICKER(name="cv09_3", len=lcv)
sd09_3 = KSEXT(name="sd09_3", len=0.24, k2=ksd09_3)
ch09_3 = HKICKER(name="ch09_3", len=lch)
sf10_3 = KSEXT(name="sf10_3", len=0.24, k2=ksf10_3)
cv10_3 = VKICKER(name="cv10_3", len=lcv)
sd10_3 = KSEXT(name="sd10_3", len=0.24, k2=ksd10_3)
ch10_3 = HKICKER(name="ch10_3", len=lch)
sf11_3 = KSEXT(name="sf11_3", len=0.24, k2=ksf11_3)
cv11_3 = VKICKER(name="cv11_3", len=lcv)
sd11_3 = KSEXT(name="sd11_3", len=0.24, k2=ksd11_3)
ch11_3 = HKICKER(name="ch11_3", len=lch)
sf12_3 = KSEXT(name="sf12_3", len=0.24, k2=ksf12_3)
cv12_3 = VKICKER(name="cv12_3", len=lcv)
sd12_3 = KSEXT(name="sd12_3", len=0.24, k2=ksd12_3)
ch12_3 = HKICKER(name="ch12_3", len=lch)
sf13_3 = KSEXT(name="sf13_3", len=0.24, k2=ksf13_3)
cv13_3 = VKICKER(name="cv13_3", len=lcv)
sd13_3 = KSEXT(name="sd13_3", len=0.24, k2=ksd13_3)
ch13_3 = HKICKER(name="ch13_3", len=lch)
sf14_3 = KSEXT(name="sf14_3", len=0.24, k2=ksf14_3)
cv14_3 = VKICKER(name="cv14_3", len=lcv)
sd14_3 = KSEXT(name="sd14_3", len=0.24, k2=ksd14_3)
ch14_3 = HKICKER(name="ch14_3", len=lch)
sf15_3 = KSEXT(name="sf15_3", len=0.24, k2=ksf15_3)
cv15_3 = VKICKER(name="cv15_3", len=lcv)
sd15_3 = KSEXT(name="sd15_3", len=0.24, k2=ksd15_3)
ch15_3 = HKICKER(name="ch15_3", len=lch)
sf16_3 = KSEXT(name="sf16_3", len=0.24, k2=ksf16_3)
cv16_3 = VKICKER(name="cv16_3", len=lcv)
sd16_3 = KSEXT(name="sd16_3", len=0.24, k2=ksd16_3)
ch16_3 = HKICKER(name="ch16_3", len=lch)
sf17_3 = KSEXT(name="sf17_3", len=0.24, k2=ksf17_3)
db23_3 = RBEND(name="db23_3", len=2.726, angle=rot_angle)
cv17_3 = VKICKER(name="cv17_3", len=lcv)
ch17_3 = HKICKER(name="ch17_3", len=lch)
sf18_3 = KSEXT(name="sf18_3", len=0.24, k2=ksf18_3)
hqd17_3 = KQUAD(name="hqd17_3", len=0.4, k1=kd17_3)
hqf16_3 = KQUAD(name="hqf16_3", len=0.4, k1=kf16_3)
hqd15_3 = KQUAD(name="hqd15_3", len=0.4, k1=kd15_3)
hqf14_3 = KQUAD(name="hqf14_3", len=0.4, k1=kf14_3)
hqd13_3 = KQUAD(name="hqd13_3", len=0.4, k1=kd13_3)
hqf12_3 = KQUAD(name="hqf12_3", len=0.4, k1=kf12_3)
hqd11_3 = KQUAD(name="hqd11_3", len=0.4, k1=kd11_3)
mq11 = MARKER(name="mq11")
hqf10_3 = KQUAD(name="hqf10_3", len=0.4, k1=kf10_3)
mq10 = MARKER(name="mq10")
hqd9_3 = KQUAD(name="hqd9_3", len=0.4, k1=kd9_3)
mq9 = MARKER(name="mq9")
hqf8_3 = KQUAD(name="hqf8_3", len=0.4, k1=kf8_3)
mq8 = MARKER(name="mq8")
hqd7_3 = KQUAD(name="hqd7_3", len=0.4, k1=kd7_3)
mq7 = MARKER(name="mq7")
hqf6_3 = KQUAD(name="hqf6_3", len=0.4, k1=kf6_3)
mq6 = MARKER(name="mq6")
hqd5_3 = KQUAD(name="hqd5_3", len=0.4, k1=kd5_3)
mq5 = MARKER(name="mq5")
hqf4_3 = KQUAD(name="hqf4_3", len=0.4, k1=kf4_3)
mq4 = MARKER(name="mq4")
hqd3_3 = KQUAD(name="hqd3_3", len=0.4, k1=kd3_3)
mq3 = MARKER(name="mq3")
db4p = RBEND(name="db4p", len=2.726, angle=0.0109083)
hqf2_3 = KQUAD(name="hqf2_3", len=0.4, k1=kf2_3)
mq2 = MARKER(name="mq2")
hqd1_3 = KQUAD(name="hqd1_3", len=0.4, k1=kd1_3)
mq1 = MARKER(name="mq1")
ip4 = MARKER(name="ip4")
hqd1_4 = KQUAD(name="hqd1_4", len=0.4, k1=kd1_4)
hqf2_4 = KQUAD(name="hqf2_4", len=0.4, k1=kf2_4)
db4m = RBEND(name="db4m", len=2.726, angle=-0.0109083)
hqd3_4 = KQUAD(name="hqd3_4", len=0.4, k1=kd3_4)
hqf4_4 = KQUAD(name="hqf4_4", len=0.4, k1=kf4_4)
hqd5_4 = KQUAD(name="hqd5_4", len=0.4, k1=kd5_4)
hqf6_4 = KQUAD(name="hqf6_4", len=0.4, k1=kf6_4)
hqd7_4 = KQUAD(name="hqd7_4", len=0.4, k1=kd7_4)
db23_4 = RBEND(name="db23_4", len=2.726, angle=rot_angle)
hqf8_4 = KQUAD(name="hqf8_4", len=0.4, k1=kf8_4)
hqd9_4 = KQUAD(name="hqd9_4", len=0.4, k1=kd9_4)
hqf10_4 = KQUAD(name="hqf10_4", len=0.4, k1=kf10_4)
hqd11_4 = KQUAD(name="hqd11_4", len=0.4, k1=kd11_4)
hqf12_4 = KQUAD(name="hqf12_4", len=0.4, k1=kf12_4)
hqd13_4 = KQUAD(name="hqd13_4", len=0.4, k1=kd13_4)
hqf14_4 = KQUAD(name="hqf14_4", len=0.4, k1=kf14_4)
hqd15_4 = KQUAD(name="hqd15_4", len=0.4, k1=kd15_4)
hqf16_4 = KQUAD(name="hqf16_4", len=0.4, k1=kf16_4)
hqd17_4 = KQUAD(name="hqd17_4", len=0.4, k1=kd17_4)
hqf18_4 = KQUAD(name="hqf18_4", len=0.4, k1=kf18_4)
sfm1_5 = KSEXT(name="sfm1_5", len=0.57, k2=ksfm1_5)
hqf_5 = KQUAD(name="hqf_5", len=0.25, k1=kf_5)
ch00_5 = HKICKER(name="ch00_5", len=lch)
hqd_5 = KQUAD(name="hqd_5", len=0.25, k1=kd_5)
cv00_5 = VKICKER(name="cv00_5", len=lcv)
sf00_5 = KSEXT(name="sf00_5", len=0.57, k2=ksf00_5)
ch01_5 = HKICKER(name="ch01_5", len=lch)
sd01_5 = KSEXT(name="sd01_5", len=0.57, k2=ksd01_5)
cv01_5 = VKICKER(name="cv01_5", len=lcv)
sf01_5 = KSEXT(name="sf01_5", len=0.24, k2=ksf01_5)
ch02_5 = HKICKER(name="ch02_5", len=lch)
sd02_5 = KSEXT(name="sd02_5", len=0.57, k2=ksd02_5)
cv02_5 = VKICKER(name="cv02_5", len=lcv)
sf02_5 = KSEXT(name="sf02_5", len=0.24, k2=ksf02_5)
ch03_5 = HKICKER(name="ch03_5", len=lch)
sd03_5 = KSEXT(name="sd03_5", len=0.57, k2=ksd03_5)
cv03_5 = VKICKER(name="cv03_5", len=lcv)
sf03_5 = KSEXT(name="sf03_5", len=0.24, k2=ksf03_5)
ch04_5 = HKICKER(name="ch04_5", len=lch)
sd04_5 = KSEXT(name="sd04_5", len=0.57, k2=ksd04_5)
cv04_5 = VKICKER(name="cv04_5", len=lcv)
sf04_5 = KSEXT(name="sf04_5", len=0.24, k2=ksf04_5)
ch05_5 = HKICKER(name="ch05_5", len=lch)
sd05_5 = KSEXT(name="sd05_5", len=0.57, k2=ksd05_5)
cv05_5 = VKICKER(name="cv05_5", len=lcv)
sf05_5 = KSEXT(name="sf05_5", len=0.24, k2=ksf05_5)
ch06_5 = HKICKER(name="ch06_5", len=lch)
sd06_5 = KSEXT(name="sd06_5", len=0.57, k2=ksd06_5)
cv06_5 = VKICKER(name="cv06_5", len=lcv)
sf06_5 = KSEXT(name="sf06_5", len=0.24, k2=ksf06_5)
ch07_5 = HKICKER(name="ch07_5", len=lch)
sd07_5 = KSEXT(name="sd07_5", len=0.57, k2=ksd07_5)
cv07_5 = VKICKER(name="cv07_5", len=lcv)
sf07_5 = KSEXT(name="sf07_5", len=0.24, k2=ksf07_5)
ch08_5 = HKICKER(name="ch08_5", len=lch)
sd08_5 = KSEXT(name="sd08_5", len=0.57, k2=ksd08_5)
cv08_5 = VKICKER(name="cv08_5", len=lcv)
sf08_5 = KSEXT(name="sf08_5", len=0.24, k2=ksf08_5)
ch09_5 = HKICKER(name="ch09_5", len=lch)
sd09_5 = KSEXT(name="sd09_5", len=0.57, k2=ksd09_5)
cv09_5 = VKICKER(name="cv09_5", len=lcv)
sf09_5 = KSEXT(name="sf09_5", len=0.24, k2=ksf09_5)
ch10_5 = HKICKER(name="ch10_5", len=lch)
sd10_5 = KSEXT(name="sd10_5", len=0.57, k2=ksd10_5)
cv10_5 = VKICKER(name="cv10_5", len=lcv)
sf10_5 = KSEXT(name="sf10_5", len=0.24, k2=ksf10_5)
ch11_5 = HKICKER(name="ch11_5", len=lch)
sd11_5 = KSEXT(name="sd11_5", len=0.57, k2=ksd11_5)
cv11_5 = VKICKER(name="cv11_5", len=lcv)
sf11_5 = KSEXT(name="sf11_5", len=0.24, k2=ksf11_5)
ch12_5 = HKICKER(name="ch12_5", len=lch)
sd12_5 = KSEXT(name="sd12_5", len=0.57, k2=ksd12_5)
cv12_5 = VKICKER(name="cv12_5", len=lcv)
sf12_5 = KSEXT(name="sf12_5", len=0.24, k2=ksf12_5)
ch13_5 = HKICKER(name="ch13_5", len=lch)
sd13_5 = KSEXT(name="sd13_5", len=0.57, k2=ksd13_5)
cv13_5 = VKICKER(name="cv13_5", len=lcv)
sf13_5 = KSEXT(name="sf13_5", len=0.24, k2=ksf13_5)
ch14_5 = HKICKER(name="ch14_5", len=lch)
sd14_5 = KSEXT(name="sd14_5", len=0.57, k2=ksd14_5)
cv14_5 = VKICKER(name="cv14_5", len=lcv)
sf14_5 = KSEXT(name="sf14_5", len=0.24, k2=ksf14_5)
hqf_5c = KQUAD(name="hqf_5c", len=0.25, k1=kf_5c)
ch15_5 = HKICKER(name="ch15_5", len=lch)
edge1_001 = thinMULTIPOLE(name="edge1_001", PolynomB=[ 0.0,  -3.8805440506512e-05, 0.0, 0.0])
d01a_001 = SBEND(name="d01a_001", len=2.7260648807936, angle=0.0095350523664854)
edge2_001 = thinMULTIPOLE(name="edge2_001", PolynomB=[ 0.0,  5.4527671338493e-06, 0.0, 0.0])
edge3_001 = thinMULTIPOLE(name="edge3_001", PolynomB=[ 0.0,  -5.4527671338493e-06, 0.0, 0.0])
d23_001 = SBEND(name="d23_001", len=0.89140036106128, angle=0.0031178821832549)
d01b_001 = SBEND(name="d01b_001", len=2.7260648807936, angle=0.0095350523664854)
sd15_5 = KSEXT(name="sd15_5", len=0.57, k2=ksd15_5)
hqd_5c = KQUAD(name="hqd_5c", len=0.25, k1=kd_5c)
cv15_5 = VKICKER(name="cv15_5", len=lcv)
sf15_5 = KSEXT(name="sf15_5", len=0.24, k2=ksf15_5)
hqf_5b = KQUAD(name="hqf_5b", len=0.25, k1=kf_5b)
ch16_5 = HKICKER(name="ch16_5", len=lch)
sd16_5 = KSEXT(name="sd16_5", len=0.57, k2=ksd16_5)
hqd_5b = KQUAD(name="hqd_5b", len=0.25, k1=kd_5b)
cv16_5 = VKICKER(name="cv16_5", len=lcv)
sf16_5 = KSEXT(name="sf16_5", len=0.24, k2=ksf16_5)
hqf_5a = KQUAD(name="hqf_5a", len=0.25, k1=kf_5a)
hqd_5a = KQUAD(name="hqd_5a", len=0.25, k1=kd_5a)
hqss1_5 = KQUAD(name="hqss1_5", len=lqss1 / 2.0, k1=kqss1_5)
hqss2_5 = KQUAD(name="hqss2_5", len=lqss2 / 2.0, k1=kqss2_5)
hqss3_5 = KQUAD(name="hqss3_5", len=lqss3 / 2.0, k1=kqss3_5)
hqss4_5 = KQUAD(name="hqss4_5", len=lqss4 / 2.0, k1=kqss4_5)
hqss5_5 = KQUAD(name="hqss5_5", len=lqss5 / 2.0, k1=kqss5_5)
qff1_5 = KQUAD(name="qff1_5", len=0.5, k1=kff1_5)
db23_5 = RBEND(name="db23_5", len=3.8, angle=theta1 / 5)
qff2_5 = KQUAD(name="qff2_5", len=0.5, k1=kff2_5)
qff3_5 = KQUAD(name="qff3_5", len=0.5, k1=kff3_5)
qff4_5 = KQUAD(name="qff4_5", len=0.5, k1=kff4_5)
qff5_5 = KQUAD(name="qff5_5", len=0.5, k1=kff5_5)
qff6_5 = KQUAD(name="qff6_5", len=0.5, k1=kff6_5)
hqls1_5 = KQUAD(name="hqls1_5", len=lqls1 / 2.0, k1=kqls1_5)
hqls2_5 = KQUAD(name="hqls2_5", len=lqls2 / 2.0, k1=kqls2_5)
hqls3_5 = KQUAD(name="hqls3_5", len=lqls3 / 2.0, k1=kqls3_5)
hqls4_5 = KQUAD(name="hqls4_5", len=lqls4 / 2.0, k1=kqls4_5)
hqls5_5 = KQUAD(name="hqls5_5", len=lqls3 / 2.0, k1=kqls5_5)
hqls6_5 = KQUAD(name="hqls6_5", len=lqls2 / 2.0, k1=kqls6_5)
hqls7_5 = KQUAD(name="hqls7_5", len=lqls1 / 2.0, k1=kqls7_5)
mlrf_6 = MARKER(name="mlrf_6")
q14ef_6 = KQUAD(name="q14ef_6", len=1.2, k1=k14ef_6)
sq14ef_6 = KQUAD(name="sq14ef_6", len=0.25, k1=0.0)
d5ef_6 = RBEND(name="d5ef_6", len=l01_str, angle=db5ef_ang_6)
q13ef_6 = KQUAD(name="q13ef_6", len=1.2, k1=k13ef_6)
sq13ef_6 = KQUAD(name="sq13ef_6", len=0.25, k1=0.0)
q12ef_6 = KQUAD(name="q12ef_6", len=1.2, k1=k12ef_6)
sq12ef_6 = KQUAD(name="sq12ef_6", len=0.25, k1=0.0)
q11ef_6 = KQUAD(name="q11ef_6", len=1.2, k1=k11ef_6)
sq11ef_6 = KQUAD(name="sq11ef_6", len=0.25, k1=0.0)
q10ef_6 = KQUAD(name="q10ef_6", len=1.2, k1=k10ef_6)
sq10ef_6 = KQUAD(name="sq10ef_6", len=0.25, k1=0.0)
d4ef_6 = RBEND(name="d4ef_6", len=l01_str, angle=db4ef_ang_6)
q9ef_6 = KQUAD(name="q9ef_6", len=1.2, k1=k9ef_6)
sq9ef_6 = KQUAD(name="sq9ef_6", len=0.25, k1=0.0)
q8ef_6 = KQUAD(name="q8ef_6", len=1.2, k1=k8ef_6)
sq8ef_6 = KQUAD(name="sq8ef_6", len=0.25, k1=0.0)
d3ef_6 = RBEND(name="d3ef_6", len=l01_str, angle=db3ef_ang_6)
q7ef_6 = KQUAD(name="q7ef_6", len=1.2, k1=k7ef_6)
sq7ef_6 = KQUAD(name="sq7ef_6", len=0.25, k1=0.0)
m_edetect1 = MARKER(name="m_edetect1")
q6ef_6 = KQUAD(name="q6ef_6", len=1.2, k1=k6ef_6)
sq6ef_6 = KQUAD(name="sq6ef_6", len=0.25, k1=0.0)
m_edetect = MARKER(name="m_edetect")
q5ef_6 = KQUAD(name="q5ef_6", len=1.2, k1=k5ef_6)
sq5ef_6 = KQUAD(name="sq5ef_6", len=0.25, k1=0.0)
q4ef_6 = KQUAD(name="q4ef_6", len=1.2, k1=k4ef_6)
sq4ef_6 = KQUAD(name="sq4ef_6", len=0.25, k1=0.0)
q3ef_6 = KQUAD(name="q3ef_6", len=0.6, k1=k3ef_6)
q2ef_6 = KQUAD(name="q2ef_6", len=0.6, k1=k2ef_6)
sq2ef_6 = KQUAD(name="sq2ef_6", len=0.25, k1=0.0)
d2ef_6 = RBEND(name="d2ef_6", len=l01_str, angle=db2ef_ang_6)
d1ef_6 = RBEND(name="d1ef_6", len=l23_str, angle=db1ef_ang_6)
mcoll_mask = MARKER(name="mcoll_mask")
q1ef_6 = KQUAD(name="q1ef_6", len=1.61, k1=k1ef_6)
q0ef_6 = KQUAD(name="q0ef_6", len=1.2, k1=k0ef_6)


# put all elements in a Dict
vars = Dict()
vars["ip6"] = ip6
vars["q1er_6"] = q1er_6
vars["q2er_6"] = q2er_6
vars["d2er_6"] = d2er_6
vars["sq"] = sq
vars["sq3er_6"] = sq3er_6
vars["qir06"] = qir06
vars["q3er_6"] = q3er_6
vars["sq4er_6"] = sq4er_6
vars["q4er_6"] = q4er_6
vars["sq5er_6"] = sq5er_6
vars["q5er_6"] = q5er_6
vars["d3er_6"] = d3er_6
vars["sq6er_6"] = sq6er_6
vars["qir12"] = qir12
vars["q6er_6"] = q6er_6
vars["sq7er_6"] = sq7er_6
vars["q7er_6"] = q7er_6
vars["sq8er_6"] = sq8er_6
vars["q8er_6"] = q8er_6
vars["rf_crab"] = rf_crab
vars["sq9er_6"] = sq9er_6
vars["q9er_6"] = q9er_6
vars["sq10er_6"] = sq10er_6
vars["q10er_6"] = q10er_6
vars["d4er_6"] = d4er_6
vars["q11er_6"] = q11er_6
vars["q12er_6"] = q12er_6
vars["q13er_6"] = q13er_6
vars["q14er_6"] = q14er_6
vars["q15er_6"] = q15er_6
vars["mlrr_6"] = mlrr_6
vars["mrot4"] = mrot4
vars["hsol20_6"] = hsol20_6
vars["hqls7_6"] = hqls7_6
vars["hqls6_6"] = hqls6_6
vars["hqls5_6"] = hqls5_6
vars["hqls4_6"] = hqls4_6
vars["hqls3_6"] = hqls3_6
vars["hqls2_6"] = hqls2_6
vars["hqls1_6"] = hqls1_6
vars["mrot3"] = mrot3
vars["qff6_6"] = qff6_6
vars["db23_6"] = db23_6
vars["qff5_6"] = qff5_6
vars["qff4_6"] = qff4_6
vars["qff3_6"] = qff3_6
vars["qff2_6"] = qff2_6
vars["qff1_6"] = qff1_6
vars["mrot2"] = mrot2
vars["hsol5_6"] = hsol5_6
vars["hqss5_6"] = hqss5_6
vars["hqss4_6"] = hqss4_6
vars["hqss3_6"] = hqss3_6
vars["hqss2_6"] = hqss2_6
vars["hqss1_6"] = hqss1_6
vars["mrot1"] = mrot1
vars["hqd_6a"] = hqd_6a
vars["hqf_6a"] = hqf_6a
vars["marc_beg"] = marc_beg
vars["sf01_7"] = sf01_7
vars["edge1_002"] = edge1_002
vars["d01a_002"] = d01a_002
vars["edge2_002"] = edge2_002
vars["edge3_002"] = edge3_002
vars["d23_002"] = d23_002
vars["d01b_002"] = d01b_002
vars["cv01_7"] = cv01_7
vars["hqd_6b"] = hqd_6b
vars["sd01_7"] = sd01_7
vars["ch01_7"] = ch01_7
vars["hqf_6b"] = hqf_6b
vars["sf02_7"] = sf02_7
vars["cv02_7"] = cv02_7
vars["hqd_6c"] = hqd_6c
vars["sd02_7"] = sd02_7
vars["ch02_7"] = ch02_7
vars["hqf_6c"] = hqf_6c
vars["sf03_7"] = sf03_7
vars["edge1_000"] = edge1_000
vars["d01a_000"] = d01a_000
vars["edge2_000"] = edge2_000
vars["edge3_000"] = edge3_000
vars["d23_000"] = d23_000
vars["d01b_000"] = d01b_000
vars["cv03_7"] = cv03_7
vars["hqd_7"] = hqd_7
vars["sd03_7"] = sd03_7
vars["ch03_7"] = ch03_7
vars["hqf_7"] = hqf_7
vars["sf04_7"] = sf04_7
vars["cv04_7"] = cv04_7
vars["sd04_7"] = sd04_7
vars["ch04_7"] = ch04_7
vars["sf05_7"] = sf05_7
vars["cv05_7"] = cv05_7
vars["sd05_7"] = sd05_7
vars["ch05_7"] = ch05_7
vars["sf06_7"] = sf06_7
vars["cv06_7"] = cv06_7
vars["sd06_7"] = sd06_7
vars["ch06_7"] = ch06_7
vars["sf07_7"] = sf07_7
vars["cv07_7"] = cv07_7
vars["sd07_7"] = sd07_7
vars["ch07_7"] = ch07_7
vars["sf08_7"] = sf08_7
vars["cv08_7"] = cv08_7
vars["sd08_7"] = sd08_7
vars["ch08_7"] = ch08_7
vars["sf09_7"] = sf09_7
vars["cv09_7"] = cv09_7
vars["sd09_7"] = sd09_7
vars["ch09_7"] = ch09_7
vars["sf10_7"] = sf10_7
vars["cv10_7"] = cv10_7
vars["sd10_7"] = sd10_7
vars["ch10_7"] = ch10_7
vars["sf11_7"] = sf11_7
vars["cv11_7"] = cv11_7
vars["sd11_7"] = sd11_7
vars["ch11_7"] = ch11_7
vars["sf12_7"] = sf12_7
vars["cv12_7"] = cv12_7
vars["sd12_7"] = sd12_7
vars["ch12_7"] = ch12_7
vars["sf13_7"] = sf13_7
vars["cv13_7"] = cv13_7
vars["sd13_7"] = sd13_7
vars["ch13_7"] = ch13_7
vars["sf14_7"] = sf14_7
vars["cv14_7"] = cv14_7
vars["sd14_7"] = sd14_7
vars["ch14_7"] = ch14_7
vars["hqf_7c"] = hqf_7c
vars["sf15_7"] = sf15_7
vars["edge1_003"] = edge1_003
vars["d01a_003"] = d01a_003
vars["edge2_003"] = edge2_003
vars["edge3_003"] = edge3_003
vars["d23_003"] = d23_003
vars["d01b_003"] = d01b_003
vars["cv15_7"] = cv15_7
vars["hqd_7c"] = hqd_7c
vars["sd15_7"] = sd15_7
vars["ch15_7"] = ch15_7
vars["hqf_7b"] = hqf_7b
vars["sf16_7"] = sf16_7
vars["cv16_7"] = cv16_7
vars["hqd_7b"] = hqd_7b
vars["sd16_7"] = sd16_7
vars["ch16_7"] = ch16_7
vars["hqf_7a"] = hqf_7a
vars["marc_end"] = marc_end
vars["hqd_7a"] = hqd_7a
vars["hsol5_8"] = hsol5_8
vars["hqss1_7"] = hqss1_7
vars["hqss2_7"] = hqss2_7
vars["hqss3_7"] = hqss3_7
vars["hqss4_7"] = hqss4_7
vars["hqss5_7"] = hqss5_7
vars["qff1_7"] = qff1_7
vars["db23_7"] = db23_7
vars["qff2_7"] = qff2_7
vars["qff3_7"] = qff3_7
vars["qff4_7"] = qff4_7
vars["qff5_7"] = qff5_7
vars["qff6_7"] = qff6_7
vars["hsol20_8"] = hsol20_8
vars["hqls1_7"] = hqls1_7
vars["hqls2_7"] = hqls2_7
vars["hqls3_7"] = hqls3_7
vars["hqls4_7"] = hqls4_7
vars["hqls5_7"] = hqls5_7
vars["hqls6_7"] = hqls6_7
vars["hqls7_7"] = hqls7_7
vars["mlrf_8"] = mlrf_8
vars["q15ef_8"] = q15ef_8
vars["d3ef_8"] = d3ef_8
vars["q14ef_8"] = q14ef_8
vars["q13ef_8"] = q13ef_8
vars["q12ef_8"] = q12ef_8
vars["q11ef_8"] = q11ef_8
vars["d2ef_8"] = d2ef_8
vars["q10ef_8"] = q10ef_8
vars["q9ef_8"] = q9ef_8
vars["q8ef_8"] = q8ef_8
vars["q7ef_8"] = q7ef_8
vars["q6ef_8"] = q6ef_8
vars["q5ef_8"] = q5ef_8
vars["q4ef_8"] = q4ef_8
vars["q3ef_8"] = q3ef_8
vars["q2ef_8"] = q2ef_8
vars["d1ef_8"] = d1ef_8
vars["q1ef_8"] = q1ef_8
vars["q0ef_8"] = q0ef_8
vars["ip8"] = ip8
vars["q1er_8"] = q1er_8
vars["q2er_8"] = q2er_8
vars["d2er_8"] = d2er_8
vars["q3er_8"] = q3er_8
vars["d3er_8"] = d3er_8
vars["q4er_8"] = q4er_8
vars["q5er_8"] = q5er_8
vars["q6er_8"] = q6er_8
vars["q7er_8"] = q7er_8
vars["q8er_8"] = q8er_8
vars["q9er_8"] = q9er_8
vars["q10er_8"] = q10er_8
vars["q11er_8"] = q11er_8
vars["d4er_8"] = d4er_8
vars["q12er_8"] = q12er_8
vars["d5er_8"] = d5er_8
vars["q13er_8"] = q13er_8
vars["q14er_8"] = q14er_8
vars["q15er_8"] = q15er_8
vars["mlrr_8"] = mlrr_8
vars["hqls7_8"] = hqls7_8
vars["hqls6_8"] = hqls6_8
vars["hqls5_8"] = hqls5_8
vars["hqls4_8"] = hqls4_8
vars["hqls3_8"] = hqls3_8
vars["hqls2_8"] = hqls2_8
vars["hqls1_8"] = hqls1_8
vars["qff6_8"] = qff6_8
vars["db23_8"] = db23_8
vars["qff5_8"] = qff5_8
vars["qff4_8"] = qff4_8
vars["qff3_8"] = qff3_8
vars["qff2_8"] = qff2_8
vars["qff1_8"] = qff1_8
vars["hqss5_8"] = hqss5_8
vars["hqss4_8"] = hqss4_8
vars["hqss3_8"] = hqss3_8
vars["hqss2_8"] = hqss2_8
vars["hqss1_8"] = hqss1_8
vars["hqd_8a"] = hqd_8a
vars["hqf_8a"] = hqf_8a
vars["ch01_9"] = ch01_9
vars["edge1_004"] = edge1_004
vars["d01a_004"] = d01a_004
vars["edge2_004"] = edge2_004
vars["edge3_004"] = edge3_004
vars["d23_004"] = d23_004
vars["d01b_004"] = d01b_004
vars["sd01_9"] = sd01_9
vars["hqd_8b"] = hqd_8b
vars["cv01_9"] = cv01_9
vars["sf01_9"] = sf01_9
vars["hqf_8b"] = hqf_8b
vars["ch02_9"] = ch02_9
vars["sd02_9"] = sd02_9
vars["hqd_8c"] = hqd_8c
vars["cv02_9"] = cv02_9
vars["sf02_9"] = sf02_9
vars["hqf_8c"] = hqf_8c
vars["ch03_9"] = ch03_9
vars["sd03_9"] = sd03_9
vars["hqd_9"] = hqd_9
vars["cv03_9"] = cv03_9
vars["sf03_9"] = sf03_9
vars["hqf_9"] = hqf_9
vars["ch04_9"] = ch04_9
vars["sd04_9"] = sd04_9
vars["cv04_9"] = cv04_9
vars["sf04_9"] = sf04_9
vars["ch05_9"] = ch05_9
vars["sd05_9"] = sd05_9
vars["cv05_9"] = cv05_9
vars["sf05_9"] = sf05_9
vars["ch06_9"] = ch06_9
vars["sd06_9"] = sd06_9
vars["cv06_9"] = cv06_9
vars["sf06_9"] = sf06_9
vars["ch07_9"] = ch07_9
vars["sd07_9"] = sd07_9
vars["cv07_9"] = cv07_9
vars["sf07_9"] = sf07_9
vars["ch08_9"] = ch08_9
vars["sd08_9"] = sd08_9
vars["cv08_9"] = cv08_9
vars["sf08_9"] = sf08_9
vars["ch09_9"] = ch09_9
vars["sd09_9"] = sd09_9
vars["cv09_9"] = cv09_9
vars["sf09_9"] = sf09_9
vars["ch10_9"] = ch10_9
vars["sd10_9"] = sd10_9
vars["cv10_9"] = cv10_9
vars["sf10_9"] = sf10_9
vars["ch11_9"] = ch11_9
vars["sd11_9"] = sd11_9
vars["cv11_9"] = cv11_9
vars["sf11_9"] = sf11_9
vars["ch12_9"] = ch12_9
vars["sd12_9"] = sd12_9
vars["cv12_9"] = cv12_9
vars["sf12_9"] = sf12_9
vars["ch13_9"] = ch13_9
vars["sd13_9"] = sd13_9
vars["cv13_9"] = cv13_9
vars["sf13_9"] = sf13_9
vars["ch14_9"] = ch14_9
vars["sd14_9"] = sd14_9
vars["cv14_9"] = cv14_9
vars["sf14_9"] = sf14_9
vars["ch15_9"] = ch15_9
vars["sd15_9"] = sd15_9
vars["cv15_9"] = cv15_9
vars["sf15_9"] = sf15_9
vars["ch16_9"] = ch16_9
vars["sd16_9"] = sd16_9
vars["cv16_9"] = cv16_9
vars["sf16_9"] = sf16_9
vars["ch17_9"] = ch17_9
vars["db23_9"] = db23_9
vars["cv17_9"] = cv17_9
vars["sf17_9"] = sf17_9
vars["hqm22_9"] = hqm22_9
vars["mq22"] = mq22
vars["hqm21_9"] = hqm21_9
vars["mq21"] = mq21
vars["hqm20_9"] = hqm20_9
vars["mq20"] = mq20
vars["hqm19_9"] = hqm19_9
vars["mq19"] = mq19
vars["hqm18_9"] = hqm18_9
vars["mq18"] = mq18
vars["hqm17_9"] = hqm17_9
vars["mq17"] = mq17
vars["hqm16_9"] = hqm16_9
vars["mq16"] = mq16
vars["hqm15_9"] = hqm15_9
vars["mq15"] = mq15
vars["hqm14_9"] = hqm14_9
vars["mq14"] = mq14
vars["hqm13_9"] = hqm13_9
vars["mq13"] = mq13
vars["hqm12_9"] = hqm12_9
vars["mq12"] = mq12
vars["hqfss_10"] = hqfss_10
vars["mqss_10"] = mqss_10
vars["hqdss_10"] = hqdss_10
vars["hqflss_10"] = hqflss_10
vars["rf0"] = rf0
vars["hqdlss_10"] = hqdlss_10
vars["ip10"] = ip10
vars["db23_10"] = db23_10
vars["hqm12_10"] = hqm12_10
vars["hqm13_10"] = hqm13_10
vars["hqm14_10"] = hqm14_10
vars["hqm15_10"] = hqm15_10
vars["hqm16_10"] = hqm16_10
vars["hqm17_10"] = hqm17_10
vars["hqm18_10"] = hqm18_10
vars["hqm19_10"] = hqm19_10
vars["hqm20_10"] = hqm20_10
vars["hqm21_10"] = hqm21_10
vars["hqm22_10"] = hqm22_10
vars["hqf_11"] = hqf_11
vars["sf00_11"] = sf00_11
vars["cv00_11"] = cv00_11
vars["hqd_11"] = hqd_11
vars["ch00_11"] = ch00_11
vars["sf01_11"] = sf01_11
vars["cv01_11"] = cv01_11
vars["sd01_11"] = sd01_11
vars["ch01_11"] = ch01_11
vars["sf02_11"] = sf02_11
vars["cv02_11"] = cv02_11
vars["sd02_11"] = sd02_11
vars["ch02_11"] = ch02_11
vars["sf03_11"] = sf03_11
vars["cv03_11"] = cv03_11
vars["sd03_11"] = sd03_11
vars["ch03_11"] = ch03_11
vars["sf04_11"] = sf04_11
vars["cv04_11"] = cv04_11
vars["sd04_11"] = sd04_11
vars["ch04_11"] = ch04_11
vars["sf05_11"] = sf05_11
vars["cv05_11"] = cv05_11
vars["sd05_11"] = sd05_11
vars["ch05_11"] = ch05_11
vars["sf06_11"] = sf06_11
vars["cv06_11"] = cv06_11
vars["sd06_11"] = sd06_11
vars["ch06_11"] = ch06_11
vars["sf07_11"] = sf07_11
vars["cv07_11"] = cv07_11
vars["sd07_11"] = sd07_11
vars["ch07_11"] = ch07_11
vars["sf08_11"] = sf08_11
vars["cv08_11"] = cv08_11
vars["sd08_11"] = sd08_11
vars["ch08_11"] = ch08_11
vars["sf09_11"] = sf09_11
vars["cv09_11"] = cv09_11
vars["sd09_11"] = sd09_11
vars["ch09_11"] = ch09_11
vars["sf10_11"] = sf10_11
vars["cv10_11"] = cv10_11
vars["sd10_11"] = sd10_11
vars["ch10_11"] = ch10_11
vars["sf11_11"] = sf11_11
vars["cv11_11"] = cv11_11
vars["sd11_11"] = sd11_11
vars["ch11_11"] = ch11_11
vars["sf12_11"] = sf12_11
vars["cv12_11"] = cv12_11
vars["sd12_11"] = sd12_11
vars["ch12_11"] = ch12_11
vars["sf13_11"] = sf13_11
vars["cv13_11"] = cv13_11
vars["sd13_11"] = sd13_11
vars["ch13_11"] = ch13_11
vars["sf14_11"] = sf14_11
vars["cv14_11"] = cv14_11
vars["sd14_11"] = sd14_11
vars["ch14_11"] = ch14_11
vars["sf15_11"] = sf15_11
vars["cv15_11"] = cv15_11
vars["sd15_11"] = sd15_11
vars["ch15_11"] = ch15_11
vars["sf16_11"] = sf16_11
vars["cv16_11"] = cv16_11
vars["sd16_11"] = sd16_11
vars["ch16_11"] = ch16_11
vars["sf17_11"] = sf17_11
vars["db23_11"] = db23_11
vars["cv17_11"] = cv17_11
vars["ch17_11"] = ch17_11
vars["sf18_11"] = sf18_11
vars["hqm22_11"] = hqm22_11
vars["hqm21_11"] = hqm21_11
vars["hqm20_11"] = hqm20_11
vars["hqm19_11"] = hqm19_11
vars["hqm18_11"] = hqm18_11
vars["hqm17_11"] = hqm17_11
vars["hqm16_11"] = hqm16_11
vars["hqm15_11"] = hqm15_11
vars["hqm14_11"] = hqm14_11
vars["hqm13_11"] = hqm13_11
vars["hqm12_11"] = hqm12_11
vars["hqfss_12"] = hqfss_12
vars["mqss_12"] = mqss_12
vars["hqdss_12"] = hqdss_12
vars["db12m"] = db12m
vars["ip12"] = ip12
vars["db12p"] = db12p
vars["mkick_inj"] = mkick_inj
vars["mcoll_inj"] = mcoll_inj
vars["db23_12"] = db23_12
vars["hqm14_12"] = hqm14_12
vars["hqm15_12"] = hqm15_12
vars["hqm16_12"] = hqm16_12
vars["hqm17_12"] = hqm17_12
vars["hqm18_12"] = hqm18_12
vars["hqm19_12"] = hqm19_12
vars["hqm20_12"] = hqm20_12
vars["hqm21_12"] = hqm21_12
vars["hqm22_12"] = hqm22_12
vars["sfm1_1"] = sfm1_1
vars["hqf_1"] = hqf_1
vars["ch00_1"] = ch00_1
vars["hqd_1"] = hqd_1
vars["cv00_1"] = cv00_1
vars["sf00_1"] = sf00_1
vars["ch01_1"] = ch01_1
vars["sd01_1"] = sd01_1
vars["cv01_1"] = cv01_1
vars["sf01_1"] = sf01_1
vars["ch02_1"] = ch02_1
vars["sd02_1"] = sd02_1
vars["cv02_1"] = cv02_1
vars["sf02_1"] = sf02_1
vars["ch03_1"] = ch03_1
vars["sd03_1"] = sd03_1
vars["cv03_1"] = cv03_1
vars["sf03_1"] = sf03_1
vars["ch04_1"] = ch04_1
vars["sd04_1"] = sd04_1
vars["cv04_1"] = cv04_1
vars["sf04_1"] = sf04_1
vars["ch05_1"] = ch05_1
vars["sd05_1"] = sd05_1
vars["cv05_1"] = cv05_1
vars["sf05_1"] = sf05_1
vars["ch06_1"] = ch06_1
vars["sd06_1"] = sd06_1
vars["cv06_1"] = cv06_1
vars["sf06_1"] = sf06_1
vars["ch07_1"] = ch07_1
vars["sd07_1"] = sd07_1
vars["cv07_1"] = cv07_1
vars["sf07_1"] = sf07_1
vars["ch08_1"] = ch08_1
vars["sd08_1"] = sd08_1
vars["cv08_1"] = cv08_1
vars["sf08_1"] = sf08_1
vars["ch09_1"] = ch09_1
vars["sd09_1"] = sd09_1
vars["cv09_1"] = cv09_1
vars["sf09_1"] = sf09_1
vars["ch10_1"] = ch10_1
vars["sd10_1"] = sd10_1
vars["cv10_1"] = cv10_1
vars["sf10_1"] = sf10_1
vars["ch11_1"] = ch11_1
vars["sd11_1"] = sd11_1
vars["cv11_1"] = cv11_1
vars["sf11_1"] = sf11_1
vars["ch12_1"] = ch12_1
vars["sd12_1"] = sd12_1
vars["cv12_1"] = cv12_1
vars["sf12_1"] = sf12_1
vars["ch13_1"] = ch13_1
vars["sd13_1"] = sd13_1
vars["cv13_1"] = cv13_1
vars["sf13_1"] = sf13_1
vars["ch14_1"] = ch14_1
vars["sd14_1"] = sd14_1
vars["cv14_1"] = cv14_1
vars["sf14_1"] = sf14_1
vars["ch15_1"] = ch15_1
vars["sd15_1"] = sd15_1
vars["cv15_1"] = cv15_1
vars["sf15_1"] = sf15_1
vars["ch16_1"] = ch16_1
vars["sd16_1"] = sd16_1
vars["cv16_1"] = cv16_1
vars["sf16_1"] = sf16_1
vars["ch17_1"] = ch17_1
vars["db23_1"] = db23_1
vars["cv17_1"] = cv17_1
vars["sf17_1"] = sf17_1
vars["hqm22_1"] = hqm22_1
vars["hqm21_1"] = hqm21_1
vars["hqm20_1"] = hqm20_1
vars["hqm19_1"] = hqm19_1
vars["hqm18_1"] = hqm18_1
vars["hqm17_1"] = hqm17_1
vars["hqm16_1"] = hqm16_1
vars["hqm15_1"] = hqm15_1
vars["hqm14_1"] = hqm14_1
vars["hqm13_1"] = hqm13_1
vars["hqm12_1"] = hqm12_1
vars["hqdss_2"] = hqdss_2
vars["mqss_2"] = mqss_2
vars["sx41_2"] = sx41_2
vars["hqfss_2"] = hqfss_2
vars["sx42_2"] = sx42_2
vars["mcoll_h1"] = mcoll_h1
vars["sx43_2"] = sx43_2
vars["mcoll_h2"] = mcoll_h2
vars["sx44_2"] = sx44_2
vars["sx45_2"] = sx45_2
vars["mcoll_h3"] = mcoll_h3
vars["sx46_2"] = sx46_2
vars["ip2"] = ip2
vars["sx47_2"] = sx47_2
vars["sx48_2"] = sx48_2
vars["sx49_2"] = sx49_2
vars["sx50_2"] = sx50_2
vars["mlamb"] = mlamb
vars["sx51_2"] = sx51_2
vars["sx52_2"] = sx52_2
vars["db23_2"] = db23_2
vars["hqm12_2"] = hqm12_2
vars["hqm13_2"] = hqm13_2
vars["hqm14_2"] = hqm14_2
vars["hqm15_2"] = hqm15_2
vars["hqm16_2"] = hqm16_2
vars["hqm17_2"] = hqm17_2
vars["hqm18_2"] = hqm18_2
vars["hqm19_2"] = hqm19_2
vars["hqm20_2"] = hqm20_2
vars["hqm21_2"] = hqm21_2
vars["hqm22_2"] = hqm22_2
vars["hqf_3"] = hqf_3
vars["sf00_3"] = sf00_3
vars["cv00_3"] = cv00_3
vars["hqd_3"] = hqd_3
vars["ch00_3"] = ch00_3
vars["sf01_3"] = sf01_3
vars["cv01_3"] = cv01_3
vars["sd01_3"] = sd01_3
vars["ch01_3"] = ch01_3
vars["sf02_3"] = sf02_3
vars["cv02_3"] = cv02_3
vars["sd02_3"] = sd02_3
vars["ch02_3"] = ch02_3
vars["sf03_3"] = sf03_3
vars["cv03_3"] = cv03_3
vars["sd03_3"] = sd03_3
vars["ch03_3"] = ch03_3
vars["sf04_3"] = sf04_3
vars["cv04_3"] = cv04_3
vars["sd04_3"] = sd04_3
vars["ch04_3"] = ch04_3
vars["sf05_3"] = sf05_3
vars["cv05_3"] = cv05_3
vars["sd05_3"] = sd05_3
vars["ch05_3"] = ch05_3
vars["sf06_3"] = sf06_3
vars["cv06_3"] = cv06_3
vars["sd06_3"] = sd06_3
vars["ch06_3"] = ch06_3
vars["sf07_3"] = sf07_3
vars["cv07_3"] = cv07_3
vars["sd07_3"] = sd07_3
vars["ch07_3"] = ch07_3
vars["sf08_3"] = sf08_3
vars["cv08_3"] = cv08_3
vars["sd08_3"] = sd08_3
vars["ch08_3"] = ch08_3
vars["sf09_3"] = sf09_3
vars["cv09_3"] = cv09_3
vars["sd09_3"] = sd09_3
vars["ch09_3"] = ch09_3
vars["sf10_3"] = sf10_3
vars["cv10_3"] = cv10_3
vars["sd10_3"] = sd10_3
vars["ch10_3"] = ch10_3
vars["sf11_3"] = sf11_3
vars["cv11_3"] = cv11_3
vars["sd11_3"] = sd11_3
vars["ch11_3"] = ch11_3
vars["sf12_3"] = sf12_3
vars["cv12_3"] = cv12_3
vars["sd12_3"] = sd12_3
vars["ch12_3"] = ch12_3
vars["sf13_3"] = sf13_3
vars["cv13_3"] = cv13_3
vars["sd13_3"] = sd13_3
vars["ch13_3"] = ch13_3
vars["sf14_3"] = sf14_3
vars["cv14_3"] = cv14_3
vars["sd14_3"] = sd14_3
vars["ch14_3"] = ch14_3
vars["sf15_3"] = sf15_3
vars["cv15_3"] = cv15_3
vars["sd15_3"] = sd15_3
vars["ch15_3"] = ch15_3
vars["sf16_3"] = sf16_3
vars["cv16_3"] = cv16_3
vars["sd16_3"] = sd16_3
vars["ch16_3"] = ch16_3
vars["sf17_3"] = sf17_3
vars["db23_3"] = db23_3
vars["cv17_3"] = cv17_3
vars["ch17_3"] = ch17_3
vars["sf18_3"] = sf18_3
vars["hqd17_3"] = hqd17_3
vars["hqf16_3"] = hqf16_3
vars["hqd15_3"] = hqd15_3
vars["hqf14_3"] = hqf14_3
vars["hqd13_3"] = hqd13_3
vars["hqf12_3"] = hqf12_3
vars["hqd11_3"] = hqd11_3
vars["mq11"] = mq11
vars["hqf10_3"] = hqf10_3
vars["mq10"] = mq10
vars["hqd9_3"] = hqd9_3
vars["mq9"] = mq9
vars["hqf8_3"] = hqf8_3
vars["mq8"] = mq8
vars["hqd7_3"] = hqd7_3
vars["mq7"] = mq7
vars["hqf6_3"] = hqf6_3
vars["mq6"] = mq6
vars["hqd5_3"] = hqd5_3
vars["mq5"] = mq5
vars["hqf4_3"] = hqf4_3
vars["mq4"] = mq4
vars["hqd3_3"] = hqd3_3
vars["mq3"] = mq3
vars["db4p"] = db4p
vars["hqf2_3"] = hqf2_3
vars["mq2"] = mq2
vars["hqd1_3"] = hqd1_3
vars["mq1"] = mq1
vars["ip4"] = ip4
vars["hqd1_4"] = hqd1_4
vars["hqf2_4"] = hqf2_4
vars["db4m"] = db4m
vars["hqd3_4"] = hqd3_4
vars["hqf4_4"] = hqf4_4
vars["hqd5_4"] = hqd5_4
vars["hqf6_4"] = hqf6_4
vars["hqd7_4"] = hqd7_4
vars["db23_4"] = db23_4
vars["hqf8_4"] = hqf8_4
vars["hqd9_4"] = hqd9_4
vars["hqf10_4"] = hqf10_4
vars["hqd11_4"] = hqd11_4
vars["hqf12_4"] = hqf12_4
vars["hqd13_4"] = hqd13_4
vars["hqf14_4"] = hqf14_4
vars["hqd15_4"] = hqd15_4
vars["hqf16_4"] = hqf16_4
vars["hqd17_4"] = hqd17_4
vars["hqf18_4"] = hqf18_4
vars["sfm1_5"] = sfm1_5
vars["hqf_5"] = hqf_5
vars["ch00_5"] = ch00_5
vars["hqd_5"] = hqd_5
vars["cv00_5"] = cv00_5
vars["sf00_5"] = sf00_5
vars["ch01_5"] = ch01_5
vars["sd01_5"] = sd01_5
vars["cv01_5"] = cv01_5
vars["sf01_5"] = sf01_5
vars["ch02_5"] = ch02_5
vars["sd02_5"] = sd02_5
vars["cv02_5"] = cv02_5
vars["sf02_5"] = sf02_5
vars["ch03_5"] = ch03_5
vars["sd03_5"] = sd03_5
vars["cv03_5"] = cv03_5
vars["sf03_5"] = sf03_5
vars["ch04_5"] = ch04_5
vars["sd04_5"] = sd04_5
vars["cv04_5"] = cv04_5
vars["sf04_5"] = sf04_5
vars["ch05_5"] = ch05_5
vars["sd05_5"] = sd05_5
vars["cv05_5"] = cv05_5
vars["sf05_5"] = sf05_5
vars["ch06_5"] = ch06_5
vars["sd06_5"] = sd06_5
vars["cv06_5"] = cv06_5
vars["sf06_5"] = sf06_5
vars["ch07_5"] = ch07_5
vars["sd07_5"] = sd07_5
vars["cv07_5"] = cv07_5
vars["sf07_5"] = sf07_5
vars["ch08_5"] = ch08_5
vars["sd08_5"] = sd08_5
vars["cv08_5"] = cv08_5
vars["sf08_5"] = sf08_5
vars["ch09_5"] = ch09_5
vars["sd09_5"] = sd09_5
vars["cv09_5"] = cv09_5
vars["sf09_5"] = sf09_5
vars["ch10_5"] = ch10_5
vars["sd10_5"] = sd10_5
vars["cv10_5"] = cv10_5
vars["sf10_5"] = sf10_5
vars["ch11_5"] = ch11_5
vars["sd11_5"] = sd11_5
vars["cv11_5"] = cv11_5
vars["sf11_5"] = sf11_5
vars["ch12_5"] = ch12_5
vars["sd12_5"] = sd12_5
vars["cv12_5"] = cv12_5
vars["sf12_5"] = sf12_5
vars["ch13_5"] = ch13_5
vars["sd13_5"] = sd13_5
vars["cv13_5"] = cv13_5
vars["sf13_5"] = sf13_5
vars["ch14_5"] = ch14_5
vars["sd14_5"] = sd14_5
vars["cv14_5"] = cv14_5
vars["sf14_5"] = sf14_5
vars["hqf_5c"] = hqf_5c
vars["ch15_5"] = ch15_5
vars["edge1_001"] = edge1_001
vars["d01a_001"] = d01a_001
vars["edge2_001"] = edge2_001
vars["edge3_001"] = edge3_001
vars["d23_001"] = d23_001
vars["d01b_001"] = d01b_001
vars["sd15_5"] = sd15_5
vars["hqd_5c"] = hqd_5c
vars["cv15_5"] = cv15_5
vars["sf15_5"] = sf15_5
vars["hqf_5b"] = hqf_5b
vars["ch16_5"] = ch16_5
vars["sd16_5"] = sd16_5
vars["hqd_5b"] = hqd_5b
vars["cv16_5"] = cv16_5
vars["sf16_5"] = sf16_5
vars["hqf_5a"] = hqf_5a
vars["hqd_5a"] = hqd_5a
vars["hqss1_5"] = hqss1_5
vars["hqss2_5"] = hqss2_5
vars["hqss3_5"] = hqss3_5
vars["hqss4_5"] = hqss4_5
vars["hqss5_5"] = hqss5_5
vars["qff1_5"] = qff1_5
vars["db23_5"] = db23_5
vars["qff2_5"] = qff2_5
vars["qff3_5"] = qff3_5
vars["qff4_5"] = qff4_5
vars["qff5_5"] = qff5_5
vars["qff6_5"] = qff6_5
vars["hqls1_5"] = hqls1_5
vars["hqls2_5"] = hqls2_5
vars["hqls3_5"] = hqls3_5
vars["hqls4_5"] = hqls4_5
vars["hqls5_5"] = hqls5_5
vars["hqls6_5"] = hqls6_5
vars["hqls7_5"] = hqls7_5
vars["mlrf_6"] = mlrf_6
vars["q14ef_6"] = q14ef_6
vars["sq14ef_6"] = sq14ef_6
vars["d5ef_6"] = d5ef_6
vars["q13ef_6"] = q13ef_6
vars["sq13ef_6"] = sq13ef_6
vars["q12ef_6"] = q12ef_6
vars["sq12ef_6"] = sq12ef_6
vars["q11ef_6"] = q11ef_6
vars["sq11ef_6"] = sq11ef_6
vars["q10ef_6"] = q10ef_6
vars["sq10ef_6"] = sq10ef_6
vars["d4ef_6"] = d4ef_6
vars["q9ef_6"] = q9ef_6
vars["sq9ef_6"] = sq9ef_6
vars["q8ef_6"] = q8ef_6
vars["sq8ef_6"] = sq8ef_6
vars["d3ef_6"] = d3ef_6
vars["q7ef_6"] = q7ef_6
vars["sq7ef_6"] = sq7ef_6
vars["m_edetect1"] = m_edetect1
vars["q6ef_6"] = q6ef_6
vars["sq6ef_6"] = sq6ef_6
vars["m_edetect"] = m_edetect
vars["q5ef_6"] = q5ef_6
vars["sq5ef_6"] = sq5ef_6
vars["q4ef_6"] = q4ef_6
vars["sq4ef_6"] = sq4ef_6
vars["q3ef_6"] = q3ef_6
vars["q2ef_6"] = q2ef_6
vars["sq2ef_6"] = sq2ef_6
vars["d2ef_6"] = d2ef_6
vars["d1ef_6"] = d1ef_6
vars["mcoll_mask"] = mcoll_mask
vars["q1ef_6"] = q1ef_6
vars["q0ef_6"] = q0ef_6

# build line
filename = "C:/Users/WAN/Desktop/JuTrack.jl/src/unused/line.txt" # list of elements without drifts
pos_file = "C:/Users/WAN/Desktop/JuTrack.jl/src/unused/positions.txt" # list of center positions of elements without drifts
function read_element_names(filename)
    # Open the file and read its content
    file_content = read(filename, String)

    # Normalize the content by replacing commas with newlines, then split by newline
    normalized_content = replace(file_content, ',' => '\n')
    elements = split(normalized_content, '\n', keepempty=false)

    # Trim any leading or trailing whitespace from each element name
    trimmed_elements = strip.(elements)

    return trimmed_elements
end

# Call the function with the filename
element_names = read_element_names(filename)

# Display the extracted element names
# println(element_names)

# read positions
function read_numbers(filename)
    # Open the file and read its content
    file_content = read(filename, String)

    # Split the content by commas to get individual numbers as strings
    number_strings = split(file_content, ",")

    # Convert each string to a number, assuming they are all integers or floats
    numbers = parse.(Float64, number_strings)

    return numbers
end
positions = read_numbers(pos_file)

nele = length(element_names)

DriftDict = Dict()

drift_idx = []
for i in 1:nele-1
    drift_len = positions[i+1] - positions[i] - vars[element_names[i]].len/2 - vars[element_names[i+1]].len/2
    if drift_len < -0.0001
        error("Drift length between $(element_names[i]) and $(element_names[i+1]) is negative: $drift_len")
    end
    if drift_len > 0.0001
        DriftDict["D$(i)"] = DRIFT(name="D$(i)", len=drift_len)
        push!(drift_idx, i)
    end
    
end

# build lattice
lattice = []
for i in 1:nele
    push!(lattice, vars[element_names[i]])
    if i in drift_idx
        push!(lattice, DriftDict["D$(i)"])
    end
end


s = 0
for i in eachindex(lattice)
    global s
    s += lattice[i].len
end

println("Total length of the lattice: $s")

function struct_to_dict(obj)
    return Dict{String, Any}([string(name) => getfield(obj, name) for name in fieldnames(typeof(obj))])
end



myVectorDicts = [struct_to_dict(obj) for obj in lattice]

# Save as JSON
open("ESR_main.json", "w") do f
    JSON.print(f, myVectorDicts)
end

