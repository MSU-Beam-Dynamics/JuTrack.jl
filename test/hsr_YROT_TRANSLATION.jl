# this code was generated from the translation of the MAD-X file in pyAT
# dipole's k0 are treated as field errors
# star_detect_r, star_detect_f, solus and solds are replaced by drift space in pyAT. reconstruct them in here
# YROTATION and TRANSLATION are manually added
# defined CRABCAVITYs
using JuTrack
ip6d = MARKER(name="ip6d")
pbr_sol_r = YROTATION(name="pbr_sol_r",angle= 2.50000000000000014e-02)
star_detect_r = SOLENOID(name="star_detect_r", len=2.0)
pet_sol_r = TRANSLATION(name="pet_sol_r",dx= -5.00104192714922943e-02)
per_sol_r = YROTATION(name="per_sol_r",angle= -2.50000000000000014e-02)
drift_0 = DRIFT(name="drift_0", len=3.2993748371982043)
q1apr = KQUAD(name="q1apr", len=1.8, k1=-0.083018639975902)
drift_1 = DRIFT(name="drift_1", len=0.4999999999999991)
q1bpr = KQUAD(name="q1bpr", len=1.4, k1=-0.083018639975902)
drift_2 = DRIFT(name="drift_2", len=1.5)
q2pr = KQUAD(name="q2pr", len=4.5, k1=0.028769941651938916)
drift_3 = DRIFT(name="drift_3", len=23.000625162801796)
q3pr = KQUAD(name="q3pr", len=1.5, k1=0.0058543264790516325)
drift_4 = DRIFT(name="drift_4", len=0.5)
# o_crab_ip6r = DRIFT(name="o_crab_ip6r", len=15.2)
o_crab_ip6r_d1 = DRIFT(name="o_crab_ip6r_d1", len=(15.2 - 4.0)/2)
o_crab_ip6r = CRABCAVITY(name="o_crab_ip6r", len=4.0, volt=0.0, freq=197e6, energy=275e6)
o_crab_ip6r_d2 = DRIFT(name="o_crab_ip6r_d2", len=(15.2 - 4.0)/2)
drift_5 = DRIFT(name="drift_5", len=0.5)
q4pr = KQUAD(name="q4pr", len=1.5, k1=0.0264808618055792)
drift_6 = DRIFT(name="drift_6", len=16.668678088423803)
b1pr = ESBEND(name="b1pr", len=3.7000445545422247, angle=-0.017, e1=-0.0085, e2=-0.0085)
drift_7 = DRIFT(name="drift_7", len=1.0)
q5pr = KQUAD(name="q5pr", len=1.5, k1=-0.07226658354271716)
drift_8 = DRIFT(name="drift_8", len=7.4733062945146855)
yi6_b4 = MARKER(name="yi6_b4")
drift_9 = DRIFT(name="drift_9", len=0.22623960000001375)
yi6_tq4 = KQUAD(name="yi6_tq4", len=0.75, k1=0.051933404287653466)
drift_10 = DRIFT(name="drift_10", len=0.13007550000000379)
yi6_qd4 = KQUAD(name="yi6_qd4", len=1.811949, k1=0.09162681421434683)
drift_11 = DRIFT(name="drift_11", len=0.5442754999999977)
yi6_tv4 = CORRECTOR(name="yi6_tv4", len=0.0, xkick=0.0, ykick=-0.0)
yi6_oct4 = thinMULTIPOLE(name="yi6_oct4", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yi6_dec4 = thinMULTIPOLE(name="yi6_dec4", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yi6_qs4 = thinMULTIPOLE(name="yi6_qs4", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_12 = DRIFT(name="drift_12", len=0.7592830000000106)
o_rotator_ip6r = DRIFT(name="o_rotator_ip6r", len=12.123674)
drift_13 = DRIFT(name="drift_13", len=0.19738340000000676)
yi6_bh5 = MARKER(name="yi6_bh5")
drift_14 = DRIFT(name="drift_14", len=0.2262275999999872)
yi6_tq5 = KQUAD(name="yi6_tq5", len=0.75, k1=-0.025966702143826733)
drift_15 = DRIFT(name="drift_15", len=0.13104999999998768)
yi6_qf5 = KQUAD(name="yi6_qf5", len=1.11, k1=-0.08803535723678176)
drift_16 = DRIFT(name="drift_16", len=0.54525000000001)
yi6_th5 = CORRECTOR(name="yi6_th5", len=0.0, xkick=-0.0, ykick=0.0)
yi6_oct5 = thinMULTIPOLE(name="yi6_oct5", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yi6_dec5 = thinMULTIPOLE(name="yi6_dec5", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yi6_qgt5 = thinMULTIPOLE(name="yi6_qgt5", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, -0.0, 0.0, 0.0])
drift_17 = DRIFT(name="drift_17", len=1.2592830000000106)
b2pr = ESBEND(name="b2pr", len=3.700002081595341, angle=0.003674538709966448, e1=0.001837269354983224, e2=0.001837269354983224)
drift_18 = DRIFT(name="drift_18", len=0.6973834000000068)
yi6_bv6 = MARKER(name="yi6_bv6")
drift_19 = DRIFT(name="drift_19", len=0.2262275999999872)
yi6_tq6 = KQUAD(name="yi6_tq6", len=0.75, k1=0.025966702143826733)
drift_20 = DRIFT(name="drift_20", len=0.13104999999998768)
yi6_qd6 = KQUAD(name="yi6_qd6", len=1.11, k1=0.08573259111211842)
drift_21 = DRIFT(name="drift_21", len=0.54525000000001)
yi6_tv6 = CORRECTOR(name="yi6_tv6", len=0.0, xkick=0.0, ykick=-0.0)
yi6_oct6 = thinMULTIPOLE(name="yi6_oct6", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yi6_dec6 = thinMULTIPOLE(name="yi6_dec6", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yi6_qs6 = thinMULTIPOLE(name="yi6_qs6", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_22 = DRIFT(name="drift_22", len=2.544731828180005)
yi6_dh5 = ESBEND(name="yi6_dh5", len=6.915999880527922, angle=0.028514815544784158, e1=0.0, e2=0.0)
drift_23 = DRIFT(name="drift_23", len=2.544731828180005)
yi6_th7 = CORRECTOR(name="yi6_th7", len=0.0, xkick=-0.0, ykick=0.0)
yi6_oct7 = thinMULTIPOLE(name="yi6_oct7", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yi6_dec7 = thinMULTIPOLE(name="yi6_dec7", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yi6_qgt7 = thinMULTIPOLE(name="yi6_qgt7", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, -0.0, 0.0, 0.0])
drift_24 = DRIFT(name="drift_24", len=0.5453779999999995)
yi6_qf7 = KQUAD(name="yi6_qf7", len=0.929744, k1=-0.07159153156668992)
drift_25 = DRIFT(name="drift_25", len=0.29623359999999366)
yi6_b7 = MARKER(name="yi6_b7")
drift_26 = DRIFT(name="drift_26", len=0.7146650329019621)
yi6_dh6 = ESBEND(name="yi6_dh6", len=2.94942686017, angle=0.012160550077129212, e1=0.0, e2=0.0)
drift_27 = DRIFT(name="drift_27", len=0.7146650329019621)
yi6_b8 = MARKER(name="yi6_b8")
drift_28 = DRIFT(name="drift_28", len=0.29610559999997577)
yi6_qd8 = KQUAD(name="yi6_qd8", len=1.11, k1=-0.03280775561046414)
drift_29 = DRIFT(name="drift_29", len=0.54525000000001)
yi6_tv8 = CORRECTOR(name="yi6_tv8", len=0.0, xkick=0.0, ykick=-0.0)
yi6_oct8 = thinMULTIPOLE(name="yi6_oct8", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yi6_dec8 = thinMULTIPOLE(name="yi6_dec8", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yi6_qs8 = thinMULTIPOLE(name="yi6_qs8", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_30 = DRIFT(name="drift_30", len=1.276564632901966)
yi6_dh9 = ESBEND(name="yi6_dh9", len=2.94942686017, angle=0.012160550077129212, e1=0.0, e2=0.0)
drift_31 = DRIFT(name="drift_31", len=1.276564632901966)
yi6_th9 = CORRECTOR(name="yi6_th9", len=0.0, xkick=-0.0, ykick=0.0)
yi6_oct9 = thinMULTIPOLE(name="yi6_oct9", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yi6_dec9 = thinMULTIPOLE(name="yi6_dec9", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yi6_qs9 = thinMULTIPOLE(name="yi6_qs9", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_32 = DRIFT(name="drift_32", len=0.54525000000001)
yi6_qf9 = KQUAD(name="yi6_qf9", len=1.11, k1=0.09765754532799617)
drift_33 = DRIFT(name="drift_33", len=1.1072775999999749)
yi6_bh9 = MARKER(name="yi6_bh9")
drift_34 = DRIFT(name="drift_34", len=0.718979989959962)
yi6_dh8 = ESBEND(name="yi6_dh8", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_35 = DRIFT(name="drift_35", len=1.0169694232637312)
yi6_bv10 = MARKER(name="yi6_bv10")
drift_36 = DRIFT(name="drift_36", len=0.22622760000001563)
yi6_sxd10 = KSEXT(name="yi6_sxd10", len=0.75, k2=-0.6290818079014571)
drift_37 = DRIFT(name="drift_37", len=0.1310500000000161)
yi6_qd10 = KQUAD(name="yi6_qd10", len=1.11, k1=-0.08364335555329494)
drift_38 = DRIFT(name="drift_38", len=0.54525000000001)
yi6_tv10 = CORRECTOR(name="yi6_tv10", len=0.0, xkick=0.0, ykick=-0.0)
yi6_oct10 = thinMULTIPOLE(name="yi6_oct10", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yi6_dec10 = thinMULTIPOLE(name="yi6_dec10", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yi6_qs10 = thinMULTIPOLE(name="yi6_qs10", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_39 = DRIFT(name="drift_39", len=1.57618712586617)
yi6_dh10 = ESBEND(name="yi6_dh10", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_40 = DRIFT(name="drift_40", len=1.014159525866205)
yi6_bh11 = MARKER(name="yi6_bh11")
drift_41 = DRIFT(name="drift_41", len=0.22622760000001563)
yi6_sxf11 = KSEXT(name="yi6_sxf11", len=0.75, k2=0.39909333141434805)
drift_42 = DRIFT(name="drift_42", len=0.1310500000000161)
yi6_qf11 = KQUAD(name="yi6_qf11", len=1.11, k1=0.08048071257602774)
drift_43 = DRIFT(name="drift_43", len=0.29525000000001)
yi6_th11 = CORRECTOR(name="yi6_th11", len=0.5, xkick=0.0, ykick=0.0)
drift_44 = DRIFT(name="drift_44", len=1.32618712586617)
yi6_dh11 = ESBEND(name="yi6_dh11", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_45 = DRIFT(name="drift_45", len=1.014159525866205)
yi6_bv12 = MARKER(name="yi6_bv12")
drift_46 = DRIFT(name="drift_46", len=0.22622760000001563)
yi6_sxd12 = KSEXT(name="yi6_sxd12", len=0.75, k2=-0.6290818079014571)
drift_47 = DRIFT(name="drift_47", len=0.1310500000000161)
yi6_qd12 = KQUAD(name="yi6_qd12", len=1.11, k1=-0.08364335555329494)
drift_48 = DRIFT(name="drift_48", len=0.29525000000001)
yi6_tv12 = CORRECTOR(name="yi6_tv12", len=0.5, xkick=0.0, ykick=0.0)
drift_49 = DRIFT(name="drift_49", len=1.32618712586617)
yi6_dh12 = ESBEND(name="yi6_dh12", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_50 = DRIFT(name="drift_50", len=1.014159525866205)
yi6_bh13 = MARKER(name="yi6_bh13")
drift_51 = DRIFT(name="drift_51", len=0.22622760000001563)
yi6_sxf13 = KSEXT(name="yi6_sxf13", len=0.75, k2=0.39909333141434805)
drift_52 = DRIFT(name="drift_52", len=0.1310500000000161)
yi6_qf13 = KQUAD(name="yi6_qf13", len=1.11, k1=0.08048071257602774)
drift_53 = DRIFT(name="drift_53", len=0.29525000000001)
yi6_th13 = CORRECTOR(name="yi6_th13", len=0.5, xkick=0.0, ykick=0.0)
drift_54 = DRIFT(name="drift_54", len=1.32618712586617)
yi6_dh13 = ESBEND(name="yi6_dh13", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_55 = DRIFT(name="drift_55", len=1.014159525866205)
yi6_bv14 = MARKER(name="yi6_bv14")
drift_56 = DRIFT(name="drift_56", len=0.22622760000001563)
yi6_sxd14 = KSEXT(name="yi6_sxd14", len=0.75, k2=-0.6290818079014571)
drift_57 = DRIFT(name="drift_57", len=0.1310500000000161)
yi6_qd14 = KQUAD(name="yi6_qd14", len=1.11, k1=-0.08364335555329494)
drift_58 = DRIFT(name="drift_58", len=0.29525000000001)
yi6_tv14 = CORRECTOR(name="yi6_tv14", len=0.5, xkick=0.0, ykick=0.0)
drift_59 = DRIFT(name="drift_59", len=1.32618712586617)
yi6_dh14 = ESBEND(name="yi6_dh14", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_60 = DRIFT(name="drift_60", len=1.014159525866205)
yi6_bh15 = MARKER(name="yi6_bh15")
drift_61 = DRIFT(name="drift_61", len=0.22622760000001563)
yi6_sxf15 = KSEXT(name="yi6_sxf15", len=0.75, k2=0.39909333141434805)
drift_62 = DRIFT(name="drift_62", len=0.1310500000000161)
yi6_qf15 = KQUAD(name="yi6_qf15", len=1.11, k1=0.08048071257602774)
drift_63 = DRIFT(name="drift_63", len=0.29525000000001)
yi6_th15 = CORRECTOR(name="yi6_th15", len=0.5, xkick=0.0, ykick=0.0)
drift_64 = DRIFT(name="drift_64", len=1.32618712586617)
yi6_dh15 = ESBEND(name="yi6_dh15", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_65 = DRIFT(name="drift_65", len=1.014159525866205)
yi6_bv16 = MARKER(name="yi6_bv16")
drift_66 = DRIFT(name="drift_66", len=0.22622760000001563)
yi6_sxd16 = KSEXT(name="yi6_sxd16", len=0.75, k2=-0.6290818079014571)
drift_67 = DRIFT(name="drift_67", len=0.1310500000000161)
yi6_qd16 = KQUAD(name="yi6_qd16", len=1.11, k1=-0.08364335555329494)
drift_68 = DRIFT(name="drift_68", len=0.29525000000001)
yi6_tv16 = CORRECTOR(name="yi6_tv16", len=0.5, xkick=0.0, ykick=0.0)
drift_69 = DRIFT(name="drift_69", len=1.32618712586617)
yi6_dh16 = ESBEND(name="yi6_dh16", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_70 = DRIFT(name="drift_70", len=1.014159525866262)
yi6_bh17 = MARKER(name="yi6_bh17")
drift_71 = DRIFT(name="drift_71", len=0.22622760000001563)
yi6_sxf17 = KSEXT(name="yi6_sxf17", len=0.75, k2=0.39909333141434805)
drift_72 = DRIFT(name="drift_72", len=0.13104999999995925)
yi6_qf17 = KQUAD(name="yi6_qf17", len=1.11, k1=0.08048071257602774)
drift_73 = DRIFT(name="drift_73", len=0.29525000000001)
yi6_th17 = CORRECTOR(name="yi6_th17", len=0.5, xkick=0.0, ykick=0.0)
drift_74 = DRIFT(name="drift_74", len=1.3261871258662268)
yi6_dh17 = ESBEND(name="yi6_dh17", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_75 = DRIFT(name="drift_75", len=1.014159525866262)
yi6_bv18 = MARKER(name="yi6_bv18")
drift_76 = DRIFT(name="drift_76", len=0.22622760000001563)
yi6_sxd18 = KSEXT(name="yi6_sxd18", len=0.75, k2=-0.6290818079014571)
drift_77 = DRIFT(name="drift_77", len=0.13104999999995925)
yi6_qd18 = KQUAD(name="yi6_qd18", len=1.11, k1=-0.08364335555329494)
drift_78 = DRIFT(name="drift_78", len=0.29525000000001)
yi6_tv18 = CORRECTOR(name="yi6_tv18", len=0.5, xkick=0.0, ykick=0.0)
drift_79 = DRIFT(name="drift_79", len=1.3261871258662268)
yi6_dh18 = ESBEND(name="yi6_dh18", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_80 = DRIFT(name="drift_80", len=1.014159525866262)
yi6_bh19 = MARKER(name="yi6_bh19")
drift_81 = DRIFT(name="drift_81", len=0.22622760000001563)
yi6_sxf19 = KSEXT(name="yi6_sxf19", len=0.75, k2=0.39909333141434805)
drift_82 = DRIFT(name="drift_82", len=0.13104999999995925)
yi6_qf19 = KQUAD(name="yi6_qf19", len=1.11, k1=0.08048071257602774)
drift_83 = DRIFT(name="drift_83", len=0.29525000000001)
yi6_th19 = CORRECTOR(name="yi6_th19", len=0.5, xkick=0.0, ykick=0.0)
drift_84 = DRIFT(name="drift_84", len=1.3261871258662268)
yi6_dh19 = ESBEND(name="yi6_dh19", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_85 = DRIFT(name="drift_85", len=1.014159525866262)
yi6_bv20 = MARKER(name="yi6_bv20")
drift_86 = DRIFT(name="drift_86", len=0.22622760000001563)
yi6_sxd20 = KSEXT(name="yi6_sxd20", len=0.75, k2=-0.6290818079014571)
drift_87 = DRIFT(name="drift_87", len=0.13104999999995925)
yi6_qd20 = KQUAD(name="yi6_qd20", len=1.11, k1=-0.08364335555329494)
drift_88 = DRIFT(name="drift_88", len=0.29525000000001)
yi6_tv20 = CORRECTOR(name="yi6_tv20", len=0.5, xkick=0.0, ykick=0.0)
drift_89 = DRIFT(name="drift_89", len=1.3261871258662268)
yi6_dh20 = ESBEND(name="yi6_dh20", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_90 = DRIFT(name="drift_90", len=1.014159525866262)
yi7_bh21 = MARKER(name="yi7_bh21")
drift_91 = DRIFT(name="drift_91", len=0.22622760000001563)
yi7_sxf21 = KSEXT(name="yi7_sxf21", len=0.75, k2=0.39909333141434805)
drift_92 = DRIFT(name="drift_92", len=0.13104999999995925)
yi7_qf21 = KQUAD(name="yi7_qf21", len=1.11, k1=0.08048071257602774)
drift_93 = DRIFT(name="drift_93", len=0.29525000000001)
yi7_th21 = CORRECTOR(name="yi7_th21", len=0.5, xkick=0.0, ykick=0.0)
drift_94 = DRIFT(name="drift_94", len=1.3261871258662268)
yi7_dh20 = ESBEND(name="yi7_dh20", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_95 = DRIFT(name="drift_95", len=1.014159525866262)
yi7_bv20 = MARKER(name="yi7_bv20")
drift_96 = DRIFT(name="drift_96", len=0.22622760000001563)
yi7_sxd20 = KSEXT(name="yi7_sxd20", len=0.75, k2=-0.6290818079014571)
drift_97 = DRIFT(name="drift_97", len=0.13104999999995925)
yi7_qd20 = KQUAD(name="yi7_qd20", len=1.11, k1=-0.08364335555329494)
drift_98 = DRIFT(name="drift_98", len=0.29525000000001)
yi7_tv20 = CORRECTOR(name="yi7_tv20", len=0.5, xkick=0.0, ykick=0.0)
drift_99 = DRIFT(name="drift_99", len=1.3261871258662268)
yi7_dh19 = ESBEND(name="yi7_dh19", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_100 = DRIFT(name="drift_100", len=1.014159525866262)
yi7_bh19 = MARKER(name="yi7_bh19")
drift_101 = DRIFT(name="drift_101", len=0.22622760000001563)
yi7_sxf19 = KSEXT(name="yi7_sxf19", len=0.75, k2=0.39909333141434805)
drift_102 = DRIFT(name="drift_102", len=0.13104999999995925)
yi7_qf19 = KQUAD(name="yi7_qf19", len=1.11, k1=0.08048071257602774)
drift_103 = DRIFT(name="drift_103", len=0.29525000000001)
yi7_th19 = CORRECTOR(name="yi7_th19", len=0.5, xkick=0.0, ykick=0.0)
drift_104 = DRIFT(name="drift_104", len=1.3261871258662268)
yi7_dh18 = ESBEND(name="yi7_dh18", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_105 = DRIFT(name="drift_105", len=1.014159525866262)
yi7_bv18 = MARKER(name="yi7_bv18")
drift_106 = DRIFT(name="drift_106", len=0.22622760000001563)
yi7_sxd18 = KSEXT(name="yi7_sxd18", len=0.75, k2=-0.6290818079014571)
drift_107 = DRIFT(name="drift_107", len=0.13104999999995925)
yi7_qd18 = KQUAD(name="yi7_qd18", len=1.11, k1=-0.08364335555329494)
drift_108 = DRIFT(name="drift_108", len=0.54525000000001)
yi7_tv18 = CORRECTOR(name="yi7_tv18", len=0.0, xkick=0.0, ykick=-0.0)
yi7_oct18 = thinMULTIPOLE(name="yi7_oct18", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yi7_dec18 = thinMULTIPOLE(name="yi7_dec18", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yi7_qs18 = thinMULTIPOLE(name="yi7_qs18", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_109 = DRIFT(name="drift_109", len=1.5761871258662268)
yi7_dh17 = ESBEND(name="yi7_dh17", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_110 = DRIFT(name="drift_110", len=1.014159525866262)
yi7_bh17 = MARKER(name="yi7_bh17")
drift_111 = DRIFT(name="drift_111", len=0.22622760000001563)
yi7_sxf17 = KSEXT(name="yi7_sxf17", len=0.75, k2=0.39909333141434805)
drift_112 = DRIFT(name="drift_112", len=0.13104999999995925)
yi7_qf17 = KQUAD(name="yi7_qf17", len=1.11, k1=0.08048071257602774)
drift_113 = DRIFT(name="drift_113", len=0.54525000000001)
yi7_th17 = CORRECTOR(name="yi7_th17", len=0.0, xkick=-0.0, ykick=0.0)
yi7_oct17 = thinMULTIPOLE(name="yi7_oct17", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yi7_dec17 = thinMULTIPOLE(name="yi7_dec17", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yi7_qgt17 = thinMULTIPOLE(name="yi7_qgt17", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, -0.0, 0.0, 0.0])
drift_114 = DRIFT(name="drift_114", len=1.5761871258662268)
yi7_dh16 = ESBEND(name="yi7_dh16", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_115 = DRIFT(name="drift_115", len=1.014159525866262)
yi7_bv16 = MARKER(name="yi7_bv16")
drift_116 = DRIFT(name="drift_116", len=0.22622760000001563)
yi7_sxd16 = KSEXT(name="yi7_sxd16", len=0.75, k2=-0.6290818079014571)
drift_117 = DRIFT(name="drift_117", len=0.13104999999995925)
yi7_qd16 = KQUAD(name="yi7_qd16", len=1.11, k1=-0.08364335555329494)
drift_118 = DRIFT(name="drift_118", len=0.54525000000001)
yi7_tv16 = CORRECTOR(name="yi7_tv16", len=0.0, xkick=0.0, ykick=-0.0)
yi7_oct16 = thinMULTIPOLE(name="yi7_oct16", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yi7_dec16 = thinMULTIPOLE(name="yi7_dec16", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yi7_qs16 = thinMULTIPOLE(name="yi7_qs16", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_119 = DRIFT(name="drift_119", len=1.5761871258662268)
yi7_dh15 = ESBEND(name="yi7_dh15", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_120 = DRIFT(name="drift_120", len=1.014159525866262)
yi7_bh15 = MARKER(name="yi7_bh15")
drift_121 = DRIFT(name="drift_121", len=0.22622760000001563)
yi7_sxf15 = KSEXT(name="yi7_sxf15", len=0.75, k2=0.39909333141434805)
drift_122 = DRIFT(name="drift_122", len=0.13104999999995925)
yi7_qf15 = KQUAD(name="yi7_qf15", len=1.11, k1=0.08048071257602774)
drift_123 = DRIFT(name="drift_123", len=0.54525000000001)
yi7_th15 = CORRECTOR(name="yi7_th15", len=0.0, xkick=-0.0, ykick=0.0)
yi7_oct15 = thinMULTIPOLE(name="yi7_oct15", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yi7_dec15 = thinMULTIPOLE(name="yi7_dec15", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yi7_qgt15 = thinMULTIPOLE(name="yi7_qgt15", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, -0.0, 0.0, 0.0])
drift_124 = DRIFT(name="drift_124", len=1.5761871258662268)
yi7_dh14 = ESBEND(name="yi7_dh14", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_125 = DRIFT(name="drift_125", len=1.014159525866262)
yi7_bv14 = MARKER(name="yi7_bv14")
drift_126 = DRIFT(name="drift_126", len=0.22622760000001563)
yi7_sxd14 = KSEXT(name="yi7_sxd14", len=0.75, k2=-0.6290818079014571)
drift_127 = DRIFT(name="drift_127", len=0.13104999999995925)
yi7_qd14 = KQUAD(name="yi7_qd14", len=1.11, k1=-0.08364335555329494)
drift_128 = DRIFT(name="drift_128", len=0.54525000000001)
yi7_tv14 = CORRECTOR(name="yi7_tv14", len=0.0, xkick=0.0, ykick=-0.0)
yi7_oct14 = thinMULTIPOLE(name="yi7_oct14", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yi7_dec14 = thinMULTIPOLE(name="yi7_dec14", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yi7_qs14 = thinMULTIPOLE(name="yi7_qs14", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_129 = DRIFT(name="drift_129", len=1.5761871258662268)
yi7_dh13 = ESBEND(name="yi7_dh13", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_130 = DRIFT(name="drift_130", len=1.014159525866262)
yi7_bh13 = MARKER(name="yi7_bh13")
drift_131 = DRIFT(name="drift_131", len=0.22622760000001563)
yi7_sxf13 = KSEXT(name="yi7_sxf13", len=0.75, k2=0.39909333141434805)
drift_132 = DRIFT(name="drift_132", len=0.13104999999995925)
yi7_qf13 = KQUAD(name="yi7_qf13", len=1.11, k1=0.08048071257602774)
drift_133 = DRIFT(name="drift_133", len=0.54525000000001)
yi7_th13 = CORRECTOR(name="yi7_th13", len=0.0, xkick=-0.0, ykick=0.0)
yi7_oct13 = thinMULTIPOLE(name="yi7_oct13", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yi7_dec13 = thinMULTIPOLE(name="yi7_dec13", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yi7_qgt13 = thinMULTIPOLE(name="yi7_qgt13", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, -0.0, 0.0, 0.0])
drift_134 = DRIFT(name="drift_134", len=1.5761871258662268)
yi7_dh12 = ESBEND(name="yi7_dh12", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_135 = DRIFT(name="drift_135", len=1.014159525866262)
yi7_bv12 = MARKER(name="yi7_bv12")
drift_136 = DRIFT(name="drift_136", len=0.22622760000001563)
yi7_sxd12 = KSEXT(name="yi7_sxd12", len=0.75, k2=-0.6290818079014571)
drift_137 = DRIFT(name="drift_137", len=0.13104999999995925)
yi7_qd12 = KQUAD(name="yi7_qd12", len=1.11, k1=-0.08364335555329494)
drift_138 = DRIFT(name="drift_138", len=0.54525000000001)
yi7_tv12 = CORRECTOR(name="yi7_tv12", len=0.0, xkick=0.0, ykick=-0.0)
yi7_oct12 = thinMULTIPOLE(name="yi7_oct12", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yi7_dec12 = thinMULTIPOLE(name="yi7_dec12", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yi7_qs12 = thinMULTIPOLE(name="yi7_qs12", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_139 = DRIFT(name="drift_139", len=1.5761871258662268)
yi7_dh11 = ESBEND(name="yi7_dh11", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_140 = DRIFT(name="drift_140", len=1.014159525866262)
yi7_bh11 = MARKER(name="yi7_bh11")
drift_141 = DRIFT(name="drift_141", len=0.22622760000001563)
yi7_sxf11 = KSEXT(name="yi7_sxf11", len=0.75, k2=0.39909333141434805)
drift_142 = DRIFT(name="drift_142", len=0.13104999999995925)
yi7_qf11 = KQUAD(name="yi7_qf11", len=1.11, k1=0.08048071257602774)
drift_143 = DRIFT(name="drift_143", len=0.54525000000001)
yi7_th11 = CORRECTOR(name="yi7_th11", len=0.0, xkick=-0.0, ykick=0.0)
yi7_oct11 = thinMULTIPOLE(name="yi7_oct11", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yi7_dec11 = thinMULTIPOLE(name="yi7_dec11", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yi7_qgt11 = thinMULTIPOLE(name="yi7_qgt11", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, -0.0, 0.0, 0.0])
drift_144 = DRIFT(name="drift_144", len=1.5761871258662268)
yi7_dh10 = ESBEND(name="yi7_dh10", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_145 = DRIFT(name="drift_145", len=1.014159525866262)
yi7_bv10 = MARKER(name="yi7_bv10")
drift_146 = DRIFT(name="drift_146", len=0.22622760000001563)
yi7_sxd10 = KSEXT(name="yi7_sxd10", len=0.75, k2=-0.6290818079014571)
drift_147 = DRIFT(name="drift_147", len=0.13104999999995925)
yi7_qd10 = KQUAD(name="yi7_qd10", len=1.11, k1=-0.08364335555329494)
drift_148 = DRIFT(name="drift_148", len=0.54525000000001)
yi7_tv10 = CORRECTOR(name="yi7_tv10", len=0.0, xkick=0.0, ykick=-0.0)
yi7_oct10 = thinMULTIPOLE(name="yi7_oct10", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yi7_dec10 = thinMULTIPOLE(name="yi7_dec10", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yi7_qs10 = thinMULTIPOLE(name="yi7_qs10", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_149 = DRIFT(name="drift_149", len=1.5789970232642645)
yi7_dh9 = ESBEND(name="yi7_dh9", len=2.94942686017, angle=0.012160550077129212, e1=0.0, e2=0.0)
drift_150 = DRIFT(name="drift_150", len=2.9086249048253876)
qds10 = KQUAD(name="qds10", len=1.1, k1=0.031892799274564364)
drift_151 = DRIFT(name="drift_151", len=4.00862490482541)
qds09 = KQUAD(name="qds09", len=1.1, k1=0.06321814098675607)
drift_152 = DRIFT(name="drift_152", len=2.4086249048253876)
bxds9m2 = ESBEND(name="bxds9m2", len=9.440656, angle=0.03745458963, e1=0.0, e2=0.0)
drift_153 = DRIFT(name="drift_153", len=1.813141923499984)
qds08 = KQUAD(name="qds08", len=1.1, k1=-0.07412161569869379)
drift_154 = DRIFT(name="drift_154", len=1.813141923499927)
bxdsds04 = ESBEND(name="bxdsds04", len=2.949, angle=0.011022248537655212, e1=0.005511124268827606, e2=0.005511124268827606)
drift_155 = DRIFT(name="drift_155", len=2.448020151999799)
qds07 = KQUAD(name="qds07", len=1.1, k1=-0.0006750397405645514)
drift_156 = DRIFT(name="drift_156", len=2.448020151999799)
bxds9m1 = ESBEND(name="bxds9m1", len=9.440656, angle=0.03879709293, e1=0.0, e2=0.0)
drift_157 = DRIFT(name="drift_157", len=1.6862876447501094)
qds06 = KQUAD(name="qds06", len=1.1, k1=0.16425817236553306)
drift_158 = DRIFT(name="drift_158", len=2.814000000000078)
yi7_rot_hlx4 = DRIFT(name="yi7_rot_hlx4", len=2.4)
drift_159 = DRIFT(name="drift_159", len=0.21199999999998909)
yi7_rot_hlx3 = DRIFT(name="yi7_rot_hlx3", len=2.4)
drift_160 = DRIFT(name="drift_160", len=0.2239999999999327)
yi7_brot = MARKER(name="yi7_brot")
drift_161 = DRIFT(name="drift_161", len=0.2239999999999327)
yi7_rot_hlx2 = DRIFT(name="yi7_rot_hlx2", len=2.4)
drift_162 = DRIFT(name="drift_162", len=0.21199999999998909)
yi7_rot_hlx1 = DRIFT(name="yi7_rot_hlx1", len=2.4)
drift_163 = DRIFT(name="drift_163", len=2.814000000000078)
qds05 = KQUAD(name="qds05", len=1.1, k1=-0.01639386381088126)
drift_164 = DRIFT(name="drift_164", len=1.6862876447501094)
qds04b = KQUAD(name="qds04b", len=1.5, k1=-0.11265786623022173)
drift_165 = DRIFT(name="drift_165", len=1.1857661100000314)
qds04a = KQUAD(name="qds04a", len=1.5, k1=0.10702121242488563)
# o_crab_ir8w = DRIFT(name="o_crab_ir8w", len=20.0)
o_crab_ir8w_d1 = DRIFT(name="o_crab_ir8w_d1", len=(20.0-4.0)/2.0)
o_crab_ir8w = CRABCAVITY(name="o_crab_ir8w", len=4.0, volt=0.0, freq=197e6, energy=275e6)
o_crab_ir8w_d2 = DRIFT(name="o_crab_ir8w_d2", len=(20.0-4.0)/2.0)
qds03b = KQUAD(name="qds03b", len=1.5, k1=0.10054196528332839)
drift_166 = DRIFT(name="drift_166", len=1.1857661100000314)
qds03a = KQUAD(name="qds03a", len=1.5, k1=-0.09060971382060284)
drift_167 = DRIFT(name="drift_167", len=1.1857661100000314)
bxds02disp1 = ESBEND(name="bxds02disp1", len=3.7, angle=-0.017, e1=-0.0085, e2=-0.0085)
drift_168 = DRIFT(name="drift_168", len=1.1857661100000314)
qds02a = KQUAD(name="qds02a", len=1.5, k1=-0.14573491257755072)
drift_169 = DRIFT(name="drift_169", len=1.1857661100000314)
qds02 = KQUAD(name="qds02", len=1.5, k1=0.09888237063)
drift_170 = DRIFT(name="drift_170", len=1.0)
mpot1o = MARKER(name="mpot1o")
drift_171 = DRIFT(name="drift_171", len=1.0)
mpot1 = MARKER(name="mpot1")
drift_172 = DRIFT(name="drift_172", len=1.0)
mpot1i = MARKER(name="mpot1i")
drift_173 = DRIFT(name="drift_173", len=1.0)
qds01 = KQUAD(name="qds01", len=1.5, k1=0.010086390305574485)
drift_174 = DRIFT(name="drift_174", len=1.0)
bxds01b = ESBEND(name="bxds01b", len=3.0, angle=0.01, e1=0.005, e2=0.005, PolynomB = [3.48259679844040550e-03-0.01/3,0.0,0.0,0.0])
drift_175 = DRIFT(name="drift_175", len=13.001063490084562)
pwt_bxds01a = TRANSLATION(name="pwt_bxds01a",dx= -7.94029646740232881e-02)
pwr_bxds01a = YROTATION(name="pwr_bxds01a",angle= 9.25223789800000071e-03)
bxds01a = ESBEND(name="bxds01a", len=4.8, angle=-0.015, e1=-0.0075, e2=-0.0075, PolynomB = [-3.53948774577505110e-03-(-0.015)/4.8,0.0,0.0,0.0])
pdr_bxds01a = YROTATION(name="pdr_bxds01a",angle= -9.25223789800000071e-03)
pdt_bxds01a = TRANSLATION(name="pdt_bxds01a", dx= 3.50000000000000033e-02)
drift_176 = DRIFT(name="drift_176", len=0.5000145253229675)
pwt_qffds02b = TRANSLATION(name="pwt_qffds02b", dx=-1.95993612047596077e-02)
pwr_qffds02b = YROTATION(name="pwr_qffds02b",angle= 3.47914472500000013e-03)
qffds02b = KQUAD(name="qffds02b", len=2.4, k1=0.03352283848075319)
pdr_qffds02b = YROTATION(name="pdr_qffds02b",angle= -3.47914472500000013e-03)
pdt_qffds02b = TRANSLATION(name="pdt_qffds02b",dx= 1.12494307099999993e-02)
drift_177 = DRIFT(name="drift_177", len=0.5000210444734421)
pwt_qffds02a = TRANSLATION(name="pwt_qffds02a",dx=-9.85196747905143361e-03)
pwr_qffds02a = YROTATION(name="pwr_qffds02a",angle= 4.02344111400000041e-03)
qffds02a = KQUAD(name="qffds02a", len=2.6, k1=0.0372361699)
pdr_qffds02a = YROTATION(name="pdr_qffds02a",angle= -4.02344111400000041e-03)
pdt_qffds02a = TRANSLATION(name="pdt_qffds02a",dx=-6.08951193600000009e-04)
drift_178 = DRIFT(name="drift_178", len=1.0000149480931668)
pwt_qffds01b = TRANSLATION(name="pwt_qffds01b",dx=2.30857655661412253e-03)
pwr_qffds01b = YROTATION(name="pwr_qffds01b",angle= -3.68635198300000001e-03)
qffds01b = KQUAD(name="qffds01b", len=2.2, k1=-0.04855262844)
pdr_qffds01b = YROTATION(name="pdr_qffds01b",angle= 3.68635198300000001e-03)
pdt_qffds01b = TRANSLATION(name="pdt_qffds01b",dx=5.80137943799999972e-03)
drift_179 = DRIFT(name="drift_179", len=0.5000153271050749)
pwt_qffds01a = TRANSLATION(name="pwt_qffds01a",dx=5.38157212818304344e-03)
pwr_qffds01a = YROTATION(name="pwr_qffds01a",angle= -3.91498719000000009e-03)
qffds01a = KQUAD(name="qffds01a", len=2.0, k1=-0.06766538857)
pdr_qffds01a = YROTATION(name="pdr_qffds01a",angle= 3.91498719000000009e-03)
pdt_qffds01a = TRANSLATION(name="pdt_qffds01a", dx=2.44838225000000012e-03)
drift_180 = DRIFT(name="drift_180", len=0.5006143435414288)
pwt_bxsp01 = TRANSLATION(name="pwt_bxsp01",dx=1.91986815824442651e-02)
pwr_bxsp01 = YROTATION(name="pwr_bxsp01",angle= -3.20000000000000007e-02)
bxsp01 = ESBEND(name="bxsp01", len=1.2, angle=-0.006, e1=-0.003, e2=-0.003)
pdr_bxsp01 = YROTATION(name="pdr_bxsp01",angle= 3.20000000000000007e-02)
pdt_bxsp01 = TRANSLATION(name="pdt_bxsp01",dx=1.91967233700000017e-02)
drift_181 = DRIFT(name="drift_181", len=3.7999999999999545)
solds = SOLENOID(name="solds", len=2.0)
ip8 = MARKER(name="ip8")
solus = SOLENOID(name="solus", len=2.0)
drift_182 = DRIFT(name="drift_182", len=3.2999999999999545)
qffus01 = KQUAD(name="qffus01", len=1.8, k1=-0.07764430921850678)
drift_183 = DRIFT(name="drift_183", len=0.5)
qffus02 = KQUAD(name="qffus02", len=1.4, k1=-0.09120782005)
drift_184 = DRIFT(name="drift_184", len=1.5)
qffus03 = KQUAD(name="qffus03", len=4.5, k1=0.02845891824612449)
drift_185 = DRIFT(name="drift_185", len=18.0)
qus01c = KQUAD(name="qus01c", len=1.5, k1=0.007192642992882727)
drift_186 = DRIFT(name="drift_186", len=0.5)
# o_crab_ir8d = DRIFT(name="o_crab_ir8d", len=20.0)
o_crab_ir8d_d1 = DRIFT(name="o_crab_ir8d_d1", len=(20.0-4.0)/2.0)
o_crab_ir8d = CRABCAVITY(name="o_crab_ir8d", len=4.0, volt=0.0, freq=197e6, energy=275e6)
o_crab_ir8d_d2 = DRIFT(name="o_crab_ir8d_d2", len=(20.0-4.0)/2.0)
drift_187 = DRIFT(name="drift_187", len=0.5)
qus02c = KQUAD(name="qus02c", len=1.5, k1=0.09777356951217944)
drift_188 = DRIFT(name="drift_188", len=2.146524238999973)
qus03 = KQUAD(name="qus03", len=1.8, k1=-0.09702430568683512)
drift_189 = DRIFT(name="drift_189", len=2.1685242389999075)
bxus9m1 = ESBEND(name="bxus9m1", len=9.440656, angle=0.03888054621, e1=0.0, e2=0.0)
drift_190 = DRIFT(name="drift_190", len=2.5079789066264766)
qus04 = KQUAD(name="qus04", len=1.8, k1=0.09346958785959181)
drift_191 = DRIFT(name="drift_191", len=3.1999789066267113)
yo8_rot_hlx4 = DRIFT(name="yo8_rot_hlx4", len=2.4)
drift_192 = DRIFT(name="drift_192", len=0.21199999999998909)
yo8_rot_hlx3 = DRIFT(name="yo8_rot_hlx3", len=2.4)
drift_193 = DRIFT(name="drift_193", len=0.2239999999999327)
yo8_brot = MARKER(name="yo8_brot")
drift_194 = DRIFT(name="drift_194", len=0.2239999999999327)
yo8_rot_hlx2 = DRIFT(name="yo8_rot_hlx2", len=2.4)
drift_195 = DRIFT(name="drift_195", len=0.21199999999998909)
yo8_rot_hlx1 = DRIFT(name="yo8_rot_hlx1", len=2.4)
drift_196 = DRIFT(name="drift_196", len=3.1999789066267113)
qus05 = KQUAD(name="qus05", len=1.1, k1=-0.06068119269223554)
drift_197 = DRIFT(name="drift_197", len=2.5079789066264766)
bxus9m2 = ESBEND(name="bxus9m2", len=9.440656, angle=0.03617333923, e1=0.0, e2=0.0)
drift_198 = DRIFT(name="drift_198", len=2.8955731975224808)
qus06 = KQUAD(name="qus06", len=1.1, k1=0.048232583080606704)
drift_199 = DRIFT(name="drift_199", len=3.5875731975227154)
yo8_snk_hlx4 = DRIFT(name="yo8_snk_hlx4", len=2.4)
drift_200 = DRIFT(name="drift_200", len=0.21199999999998909)
yo8_snk_hlx3 = DRIFT(name="yo8_snk_hlx3", len=2.4)
drift_201 = DRIFT(name="drift_201", len=0.2239999999999327)
yo8_bsnk = MARKER(name="yo8_bsnk")
drift_202 = DRIFT(name="drift_202", len=0.2239999999999327)
yo8_snk_hlx2 = DRIFT(name="yo8_snk_hlx2", len=2.4)
drift_203 = DRIFT(name="drift_203", len=0.21199999999998909)
yo8_snk_hlx1 = DRIFT(name="yo8_snk_hlx1", len=2.4)
drift_204 = DRIFT(name="drift_204", len=3.5875731975227154)
qus07 = KQUAD(name="qus07", len=1.1, k1=-0.053096217194990955)
drift_205 = DRIFT(name="drift_205", len=2.8955731975224808)
bxus9m3 = ESBEND(name="bxus9m3", len=9.440656, angle=0.03222004565765377, e1=0.0, e2=0.0)
drift_206 = DRIFT(name="drift_206", len=2.125420175999807)
qus08 = KQUAD(name="qus08", len=1.1, k1=0.08758992046002777)
drift_207 = DRIFT(name="drift_207", len=3.2034201759997813)
qus09 = KQUAD(name="qus09", len=1.1, k1=-0.09702430568683512)
drift_208 = DRIFT(name="drift_208", len=1.6034201759998723)
yo8_dh9 = ESBEND(name="yo8_dh9", len=2.94942686017, angle=0.012160550077129212, e1=0.0, e2=0.0)
drift_209 = DRIFT(name="drift_209", len=1.5844693382359765)
yo8_th10 = CORRECTOR(name="yo8_th10", len=0.0, xkick=-0.0, ykick=0.0)
yo8_oct10 = thinMULTIPOLE(name="yo8_oct10", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yo8_dec10 = thinMULTIPOLE(name="yo8_dec10", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yo8_qs10 = thinMULTIPOLE(name="yo8_qs10", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_210 = DRIFT(name="drift_210", len=0.5452499999998963)
yo8_qf10 = KQUAD(name="yo8_qf10", len=1.11, k1=0.08048071257602774)
drift_211 = DRIFT(name="drift_211", len=0.13104999999995925)
yo8_sxf10 = KSEXT(name="yo8_sxf10", len=0.75, k2=0.39909333141434805)
drift_212 = DRIFT(name="drift_212", len=0.22622759999990194)
yo8_bh10 = MARKER(name="yo8_bh10")
drift_213 = DRIFT(name="drift_213", len=1.0316775497337858)
yo8_dh10 = ESBEND(name="yo8_dh10", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_214 = DRIFT(name="drift_214", len=1.3437051497337507)
yo8_tv11 = CORRECTOR(name="yo8_tv11", len=0.5, xkick=0.0, ykick=0.0)
drift_215 = DRIFT(name="drift_215", len=0.2952499999998963)
yo8_qd11 = KQUAD(name="yo8_qd11", len=1.11, k1=-0.08364335555329494)
drift_216 = DRIFT(name="drift_216", len=0.13104999999995925)
yo8_sxd11 = KSEXT(name="yo8_sxd11", len=0.75, k2=-0.6290818079014571)
drift_217 = DRIFT(name="drift_217", len=0.22622759999990194)
yo8_bv11 = MARKER(name="yo8_bv11")
drift_218 = DRIFT(name="drift_218", len=1.0316775497337858)
yo8_dh11 = ESBEND(name="yo8_dh11", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_219 = DRIFT(name="drift_219", len=1.3437051497337507)
yo8_th12 = CORRECTOR(name="yo8_th12", len=0.5, xkick=0.0, ykick=0.0)
drift_220 = DRIFT(name="drift_220", len=0.2952499999998963)
yo8_qf12 = KQUAD(name="yo8_qf12", len=1.11, k1=0.08048071257602774)
drift_221 = DRIFT(name="drift_221", len=0.13104999999995925)
yo8_sxf12 = KSEXT(name="yo8_sxf12", len=0.75, k2=0.39909333141434805)
drift_222 = DRIFT(name="drift_222", len=0.22622759999990194)
yo8_bh12 = MARKER(name="yo8_bh12")
drift_223 = DRIFT(name="drift_223", len=1.0316775497337858)
yo8_dh12 = ESBEND(name="yo8_dh12", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_224 = DRIFT(name="drift_224", len=1.3437051497337507)
yo8_tv13 = CORRECTOR(name="yo8_tv13", len=0.5, xkick=0.0, ykick=0.0)
drift_225 = DRIFT(name="drift_225", len=0.2952499999998963)
yo8_qd13 = KQUAD(name="yo8_qd13", len=1.11, k1=-0.08364335555329494)
drift_226 = DRIFT(name="drift_226", len=0.13104999999995925)
yo8_sxd13 = KSEXT(name="yo8_sxd13", len=0.75, k2=-0.6290818079014571)
drift_227 = DRIFT(name="drift_227", len=0.22622759999990194)
yo8_bv13 = MARKER(name="yo8_bv13")
drift_228 = DRIFT(name="drift_228", len=1.0316775497337858)
yo8_dh13 = ESBEND(name="yo8_dh13", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_229 = DRIFT(name="drift_229", len=1.3437051497337507)
yo8_th14 = CORRECTOR(name="yo8_th14", len=0.5, xkick=0.0, ykick=0.0)
drift_230 = DRIFT(name="drift_230", len=0.2952499999998963)
yo8_qf14 = KQUAD(name="yo8_qf14", len=1.11, k1=0.08048071257602774)
drift_231 = DRIFT(name="drift_231", len=0.13104999999995925)
yo8_sxf14 = KSEXT(name="yo8_sxf14", len=0.75, k2=0.39909333141434805)
drift_232 = DRIFT(name="drift_232", len=0.22622759999990194)
yo8_bh14 = MARKER(name="yo8_bh14")
drift_233 = DRIFT(name="drift_233", len=1.0316775497337858)
yo8_dh14 = ESBEND(name="yo8_dh14", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_234 = DRIFT(name="drift_234", len=1.3437051497337507)
yo8_tv15 = CORRECTOR(name="yo8_tv15", len=0.5, xkick=0.0, ykick=0.0)
drift_235 = DRIFT(name="drift_235", len=0.2952499999998963)
yo8_qd15 = KQUAD(name="yo8_qd15", len=1.11, k1=-0.08364335555329494)
drift_236 = DRIFT(name="drift_236", len=0.13104999999995925)
yo8_sxd15 = KSEXT(name="yo8_sxd15", len=0.75, k2=-0.6290818079014571)
drift_237 = DRIFT(name="drift_237", len=0.22622759999990194)
yo8_bv15 = MARKER(name="yo8_bv15")
drift_238 = DRIFT(name="drift_238", len=1.0316775497337858)
yo8_dh15 = ESBEND(name="yo8_dh15", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_239 = DRIFT(name="drift_239", len=1.3437051497337507)
yo8_th16 = CORRECTOR(name="yo8_th16", len=0.5, xkick=0.0, ykick=0.0)
drift_240 = DRIFT(name="drift_240", len=0.2952499999998963)
yo8_qf16 = KQUAD(name="yo8_qf16", len=1.11, k1=0.08048071257602774)
drift_241 = DRIFT(name="drift_241", len=0.13104999999995925)
yo8_sxf16 = KSEXT(name="yo8_sxf16", len=0.75, k2=0.39909333141434805)
drift_242 = DRIFT(name="drift_242", len=0.22622759999990194)
yo8_bh16 = MARKER(name="yo8_bh16")
drift_243 = DRIFT(name="drift_243", len=1.0316775497337858)
yo8_dh16 = ESBEND(name="yo8_dh16", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_244 = DRIFT(name="drift_244", len=1.3437051497337507)
yo8_tv17 = CORRECTOR(name="yo8_tv17", len=0.5, xkick=0.0, ykick=0.0)
drift_245 = DRIFT(name="drift_245", len=0.2952499999998963)
yo8_qd17 = KQUAD(name="yo8_qd17", len=1.11, k1=-0.08364335555329494)
drift_246 = DRIFT(name="drift_246", len=0.13104999999995925)
yo8_sxd17 = KSEXT(name="yo8_sxd17", len=0.75, k2=-0.6290818079014571)
drift_247 = DRIFT(name="drift_247", len=0.22622759999990194)
yo8_bv17 = MARKER(name="yo8_bv17")
drift_248 = DRIFT(name="drift_248", len=1.0316775497337858)
yo8_dh17 = ESBEND(name="yo8_dh17", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_249 = DRIFT(name="drift_249", len=1.3437051497337507)
yo8_th18 = CORRECTOR(name="yo8_th18", len=0.5, xkick=0.0, ykick=0.0)
drift_250 = DRIFT(name="drift_250", len=0.2952499999998963)
yo8_qf18 = KQUAD(name="yo8_qf18", len=1.11, k1=0.08048071257602774)
drift_251 = DRIFT(name="drift_251", len=0.13104999999995925)
yo8_sxf18 = KSEXT(name="yo8_sxf18", len=0.75, k2=0.39909333141434805)
drift_252 = DRIFT(name="drift_252", len=0.22622759999990194)
yo8_bh18 = MARKER(name="yo8_bh18")
drift_253 = DRIFT(name="drift_253", len=1.0316775497337858)
yo8_dh18 = ESBEND(name="yo8_dh18", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_254 = DRIFT(name="drift_254", len=1.3437051497337507)
yo8_tv19 = CORRECTOR(name="yo8_tv19", len=0.5, xkick=0.0, ykick=0.0)
drift_255 = DRIFT(name="drift_255", len=0.2952499999998963)
yo8_qd19 = KQUAD(name="yo8_qd19", len=1.11, k1=-0.08364335555329494)
drift_256 = DRIFT(name="drift_256", len=0.13104999999995925)
yo8_sxd19 = KSEXT(name="yo8_sxd19", len=0.75, k2=-0.6290818079014571)
drift_257 = DRIFT(name="drift_257", len=0.22622759999990194)
yo8_bv19 = MARKER(name="yo8_bv19")
drift_258 = DRIFT(name="drift_258", len=1.0316775497337858)
yo8_dh19 = ESBEND(name="yo8_dh19", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_259 = DRIFT(name="drift_259", len=1.3437051497337507)
yo8_th20 = CORRECTOR(name="yo8_th20", len=0.5, xkick=0.0, ykick=0.0)
drift_260 = DRIFT(name="drift_260", len=0.2952499999998963)
yo8_qf20 = KQUAD(name="yo8_qf20", len=1.11, k1=0.08048071257602774)
drift_261 = DRIFT(name="drift_261", len=0.13104999999995925)
yo8_sxf20 = KSEXT(name="yo8_sxf20", len=0.75, k2=0.39909333141434805)
drift_262 = DRIFT(name="drift_262", len=0.22622759999990194)
yo8_bh20 = MARKER(name="yo8_bh20")
drift_263 = DRIFT(name="drift_263", len=1.0316775497337858)
yo8_dh20 = ESBEND(name="yo8_dh20", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_264 = DRIFT(name="drift_264", len=1.3437051497337507)
yo9_tv21 = CORRECTOR(name="yo9_tv21", len=0.5, xkick=0.0, ykick=0.0)
drift_265 = DRIFT(name="drift_265", len=0.2952499999998963)
yo9_qd21 = KQUAD(name="yo9_qd21", len=1.11, k1=-0.08364335555329494)
drift_266 = DRIFT(name="drift_266", len=0.13104999999995925)
yo9_sxd21 = KSEXT(name="yo9_sxd21", len=0.75, k2=-0.6290818079014571)
drift_267 = DRIFT(name="drift_267", len=0.22622759999990194)
yo9_bv21 = MARKER(name="yo9_bv21")
drift_268 = DRIFT(name="drift_268", len=1.0316775497337858)
yo9_dh20 = ESBEND(name="yo9_dh20", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_269 = DRIFT(name="drift_269", len=1.3437051497337507)
yo9_th20 = CORRECTOR(name="yo9_th20", len=0.5, xkick=0.0, ykick=0.0)
drift_270 = DRIFT(name="drift_270", len=0.2952499999998963)
yo9_qf20 = KQUAD(name="yo9_qf20", len=1.11, k1=0.08048071257602774)
drift_271 = DRIFT(name="drift_271", len=0.13104999999995925)
yo9_sxf20 = KSEXT(name="yo9_sxf20", len=0.75, k2=0.39909333141434805)
drift_272 = DRIFT(name="drift_272", len=0.22622759999990194)
yo9_bh20 = MARKER(name="yo9_bh20")
drift_273 = DRIFT(name="drift_273", len=1.0316775497337858)
yo9_dh19 = ESBEND(name="yo9_dh19", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_274 = DRIFT(name="drift_274", len=1.3437051497337507)
yo9_tv19 = CORRECTOR(name="yo9_tv19", len=0.5, xkick=0.0, ykick=0.0)
drift_275 = DRIFT(name="drift_275", len=0.2952499999998963)
yo9_qd19 = KQUAD(name="yo9_qd19", len=1.11, k1=-0.08364335555329494)
drift_276 = DRIFT(name="drift_276", len=0.13104999999995925)
yo9_sxd19 = KSEXT(name="yo9_sxd19", len=0.75, k2=-0.6290818079014571)
drift_277 = DRIFT(name="drift_277", len=0.22622759999990194)
yo9_bv19 = MARKER(name="yo9_bv19")
drift_278 = DRIFT(name="drift_278", len=1.0316775497337858)
yo9_dh18 = ESBEND(name="yo9_dh18", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_279 = DRIFT(name="drift_279", len=1.5937051497337507)
yo9_th18 = CORRECTOR(name="yo9_th18", len=0.0, xkick=-0.0, ykick=0.0)
yo9_oct18 = thinMULTIPOLE(name="yo9_oct18", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yo9_dec18 = thinMULTIPOLE(name="yo9_dec18", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yo9_qgt18 = thinMULTIPOLE(name="yo9_qgt18", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, -0.0, 0.0, 0.0])
drift_280 = DRIFT(name="drift_280", len=0.5452499999998963)
yo9_qf18 = KQUAD(name="yo9_qf18", len=1.11, k1=0.08048071257602774)
drift_281 = DRIFT(name="drift_281", len=0.13104999999995925)
yo9_sxf18 = KSEXT(name="yo9_sxf18", len=0.75, k2=0.39909333141434805)
drift_282 = DRIFT(name="drift_282", len=0.22622759999990194)
yo9_bh18 = MARKER(name="yo9_bh18")
drift_283 = DRIFT(name="drift_283", len=1.0316775497337858)
yo9_dh17 = ESBEND(name="yo9_dh17", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_284 = DRIFT(name="drift_284", len=1.5937051497337507)
yo9_tv17 = CORRECTOR(name="yo9_tv17", len=0.0, xkick=0.0, ykick=-0.0)
yo9_oct17 = thinMULTIPOLE(name="yo9_oct17", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yo9_dec17 = thinMULTIPOLE(name="yo9_dec17", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yo9_qs17 = thinMULTIPOLE(name="yo9_qs17", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_285 = DRIFT(name="drift_285", len=0.5452499999998963)
yo9_qd17 = KQUAD(name="yo9_qd17", len=1.11, k1=-0.08364335555329494)
drift_286 = DRIFT(name="drift_286", len=0.13104999999995925)
yo9_sxd17 = KSEXT(name="yo9_sxd17", len=0.75, k2=-0.6290818079014571)
drift_287 = DRIFT(name="drift_287", len=0.22622759999990194)
yo9_bv17 = MARKER(name="yo9_bv17")
drift_288 = DRIFT(name="drift_288", len=1.0316775497338995)
yo9_dh16 = ESBEND(name="yo9_dh16", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_289 = DRIFT(name="drift_289", len=1.593705149733978)
yo9_th16 = CORRECTOR(name="yo9_th16", len=0.0, xkick=-0.0, ykick=0.0)
yo9_oct16 = thinMULTIPOLE(name="yo9_oct16", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yo9_dec16 = thinMULTIPOLE(name="yo9_dec16", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yo9_qgt16 = thinMULTIPOLE(name="yo9_qgt16", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, -0.0, 0.0, 0.0])
drift_290 = DRIFT(name="drift_290", len=0.5452500000001237)
yo9_qf16 = KQUAD(name="yo9_qf16", len=1.11, k1=0.08048071257602774)
drift_291 = DRIFT(name="drift_291", len=0.13104999999995925)
yo9_sxf16 = KSEXT(name="yo9_sxf16", len=0.75, k2=0.39909333141434805)
drift_292 = DRIFT(name="drift_292", len=0.22622760000012931)
yo9_bh16 = MARKER(name="yo9_bh16")
drift_293 = DRIFT(name="drift_293", len=1.0316775497340132)
yo9_dh15 = ESBEND(name="yo9_dh15", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_294 = DRIFT(name="drift_294", len=1.593705149733978)
yo9_tv15 = CORRECTOR(name="yo9_tv15", len=0.0, xkick=0.0, ykick=-0.0)
yo9_oct15 = thinMULTIPOLE(name="yo9_oct15", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yo9_dec15 = thinMULTIPOLE(name="yo9_dec15", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yo9_qs15 = thinMULTIPOLE(name="yo9_qs15", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_295 = DRIFT(name="drift_295", len=0.5452500000001237)
yo9_qd15 = KQUAD(name="yo9_qd15", len=1.11, k1=-0.08364335555329494)
drift_296 = DRIFT(name="drift_296", len=0.13104999999995925)
yo9_sxd15 = KSEXT(name="yo9_sxd15", len=0.75, k2=-0.6290818079014571)
drift_297 = DRIFT(name="drift_297", len=0.22622760000012931)
yo9_bv15 = MARKER(name="yo9_bv15")
drift_298 = DRIFT(name="drift_298", len=1.0316775497340132)
yo9_dh14 = ESBEND(name="yo9_dh14", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_299 = DRIFT(name="drift_299", len=1.593705149733978)
yo9_th14 = CORRECTOR(name="yo9_th14", len=0.0, xkick=-0.0, ykick=0.0)
yo9_oct14 = thinMULTIPOLE(name="yo9_oct14", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yo9_dec14 = thinMULTIPOLE(name="yo9_dec14", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yo9_qgt14 = thinMULTIPOLE(name="yo9_qgt14", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, -0.0, 0.0, 0.0])
drift_300 = DRIFT(name="drift_300", len=0.5452500000001237)
yo9_qf14 = KQUAD(name="yo9_qf14", len=1.11, k1=0.08048071257602774)
drift_301 = DRIFT(name="drift_301", len=0.13104999999995925)
yo9_sxf14 = KSEXT(name="yo9_sxf14", len=0.75, k2=0.39909333141434805)
drift_302 = DRIFT(name="drift_302", len=0.22622760000012931)
yo9_bh14 = MARKER(name="yo9_bh14")
drift_303 = DRIFT(name="drift_303", len=1.0316775497340132)
yo9_dh13 = ESBEND(name="yo9_dh13", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_304 = DRIFT(name="drift_304", len=1.593705149733978)
yo9_tv13 = CORRECTOR(name="yo9_tv13", len=0.0, xkick=0.0, ykick=-0.0)
yo9_oct13 = thinMULTIPOLE(name="yo9_oct13", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yo9_dec13 = thinMULTIPOLE(name="yo9_dec13", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yo9_qs13 = thinMULTIPOLE(name="yo9_qs13", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_305 = DRIFT(name="drift_305", len=0.5452500000001237)
yo9_qd13 = KQUAD(name="yo9_qd13", len=1.11, k1=-0.08364335555329494)
drift_306 = DRIFT(name="drift_306", len=0.13104999999995925)
yo9_sxd13 = KSEXT(name="yo9_sxd13", len=0.75, k2=-0.6290818079014571)
drift_307 = DRIFT(name="drift_307", len=0.22622760000012931)
yo9_bv13 = MARKER(name="yo9_bv13")
drift_308 = DRIFT(name="drift_308", len=1.0316775497340132)
yo9_dh12 = ESBEND(name="yo9_dh12", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_309 = DRIFT(name="drift_309", len=1.593705149733978)
yo9_th12 = CORRECTOR(name="yo9_th12", len=0.0, xkick=-0.0, ykick=0.0)
yo9_oct12 = thinMULTIPOLE(name="yo9_oct12", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yo9_dec12 = thinMULTIPOLE(name="yo9_dec12", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yo9_qgt12 = thinMULTIPOLE(name="yo9_qgt12", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, -0.0, 0.0, 0.0])
drift_310 = DRIFT(name="drift_310", len=0.5452500000001237)
yo9_qf12 = KQUAD(name="yo9_qf12", len=1.11, k1=0.08048071257602774)
drift_311 = DRIFT(name="drift_311", len=0.13104999999995925)
yo9_sxf12 = KSEXT(name="yo9_sxf12", len=0.75, k2=0.39909333141434805)
drift_312 = DRIFT(name="drift_312", len=0.22622760000012931)
yo9_bh12 = MARKER(name="yo9_bh12")
drift_313 = DRIFT(name="drift_313", len=1.0316775497340132)
yo9_dh11 = ESBEND(name="yo9_dh11", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_314 = DRIFT(name="drift_314", len=1.593705149733978)
yo9_tv11 = CORRECTOR(name="yo9_tv11", len=0.0, xkick=0.0, ykick=-0.0)
yo9_oct11 = thinMULTIPOLE(name="yo9_oct11", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yo9_dec11 = thinMULTIPOLE(name="yo9_dec11", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yo9_qs11 = thinMULTIPOLE(name="yo9_qs11", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_315 = DRIFT(name="drift_315", len=0.5452500000001237)
yo9_qd11 = KQUAD(name="yo9_qd11", len=1.11, k1=-0.08364335555329494)
drift_316 = DRIFT(name="drift_316", len=0.13104999999995925)
yo9_sxd11 = KSEXT(name="yo9_sxd11", len=0.75, k2=-0.6290818079014571)
drift_317 = DRIFT(name="drift_317", len=0.22622760000012931)
yo9_bv11 = MARKER(name="yo9_bv11")
drift_318 = DRIFT(name="drift_318", len=1.0316775497340132)
yo9_dh10 = ESBEND(name="yo9_dh10", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_319 = DRIFT(name="drift_319", len=1.593705149733978)
yo9_th10 = CORRECTOR(name="yo9_th10", len=0.0, xkick=-0.0, ykick=0.0)
yo9_oct10 = thinMULTIPOLE(name="yo9_oct10", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yo9_dec10 = thinMULTIPOLE(name="yo9_dec10", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yo9_qs10 = thinMULTIPOLE(name="yo9_qs10", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_320 = DRIFT(name="drift_320", len=0.5452500000001237)
yo9_qf10 = KQUAD(name="yo9_qf10", len=1.11, k1=0.08048071257602774)
drift_321 = DRIFT(name="drift_321", len=0.13104999999995925)
yo9_sxf10 = KSEXT(name="yo9_sxf10", len=0.75, k2=0.39909333141434805)
drift_322 = DRIFT(name="drift_322", len=0.22622760000012931)
yo9_bh10 = MARKER(name="yo9_bh10")
drift_323 = DRIFT(name="drift_323", len=1.0224417382355568)
yo9_dh9 = ESBEND(name="yo9_dh9", len=2.94942686017, angle=0.012160550077129212, e1=0.0, e2=0.0)
drift_324 = DRIFT(name="drift_324", len=7.5200956482653964)
yo9_bv9 = MARKER(name="yo9_bv9")
drift_325 = DRIFT(name="drift_325", len=1.1072776000000886)
yo9_qd9 = KQUAD(name="yo9_qd9", len=1.11, k1=-0.07911855595309777)
drift_326 = DRIFT(name="drift_326", len=0.5452500000001237)
yo9_tv9 = CORRECTOR(name="yo9_tv9", len=0.0, xkick=0.0, ykick=-0.0)
yo9_oct9 = thinMULTIPOLE(name="yo9_oct9", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yo9_dec9 = thinMULTIPOLE(name="yo9_dec9", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yo9_qs9 = thinMULTIPOLE(name="yo9_qs9", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_327 = DRIFT(name="drift_327", len=1.593705149733978)
yo9_dh8 = ESBEND(name="yo9_dh8", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_328 = DRIFT(name="drift_328", len=1.593705149733978)
yo9_th8 = CORRECTOR(name="yo9_th8", len=0.0, xkick=-0.0, ykick=0.0)
yo9_oct8 = thinMULTIPOLE(name="yo9_oct8", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yo9_dec8 = thinMULTIPOLE(name="yo9_dec8", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yo9_qgt8 = thinMULTIPOLE(name="yo9_qgt8", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, -0.0, 0.0, 0.0])
drift_329 = DRIFT(name="drift_329", len=0.5452500000001237)
yo9_qf8 = KQUAD(name="yo9_qf8", len=1.11, k1=0.07667711731989788)
drift_330 = DRIFT(name="drift_330", len=0.2961055999999189)
yo9_b8 = MARKER(name="yo9_b8")
drift_331 = DRIFT(name="drift_331", len=1.3184179478498663)
yo9_hlx7_4 = DRIFT(name="yo9_hlx7_4", len=2.4)
drift_332 = DRIFT(name="drift_332", len=0.21199999999998909)
yo9_hlx7_3 = DRIFT(name="yo9_hlx7_3", len=2.4)
drift_333 = DRIFT(name="drift_333", len=0.22400000000016007)
yo9_b7_1 = MARKER(name="yo9_b7_1")
drift_334 = DRIFT(name="drift_334", len=0.22400000000016007)
yo9_hlx7_2 = DRIFT(name="yo9_hlx7_2", len=2.4)
drift_335 = DRIFT(name="drift_335", len=0.21199999999998909)
yo9_hlx7_1 = DRIFT(name="yo9_hlx7_1", len=2.4)
drift_336 = DRIFT(name="drift_336", len=1.3184179478498663)
yo9_b7 = MARKER(name="yo9_b7")
drift_337 = DRIFT(name="drift_337", len=0.2962336000000505)
yo9_qd7 = KQUAD(name="yo9_qd7", len=0.929744, k1=-0.08892562201533814)
drift_338 = DRIFT(name="drift_338", len=0.5453779999998005)
yo9_tv7 = CORRECTOR(name="yo9_tv7", len=0.0, xkick=0.0, ykick=-0.0)
yo9_oct7 = thinMULTIPOLE(name="yo9_oct7", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yo9_dec7 = thinMULTIPOLE(name="yo9_dec7", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yo9_qs7 = thinMULTIPOLE(name="yo9_qs7", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_339 = DRIFT(name="drift_339", len=1.5844693382355217)
yo9_dh6 = ESBEND(name="yo9_dh6", len=2.94942686017, angle=0.012160550077129212, e1=0.0, e2=0.0)
drift_340 = DRIFT(name="drift_340", len=8.082123248265361)
yo9_th6 = CORRECTOR(name="yo9_th6", len=0.0, xkick=-0.0, ykick=0.0)
yo9_oct6 = thinMULTIPOLE(name="yo9_oct6", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yo9_dec6 = thinMULTIPOLE(name="yo9_dec6", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yo9_qgt6 = thinMULTIPOLE(name="yo9_qgt6", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, -0.0, 0.0, 0.0])
drift_341 = DRIFT(name="drift_341", len=0.5452500000001237)
yo9_qf6 = KQUAD(name="yo9_qf6", len=1.11, k1=0.09117978411883382)
drift_342 = DRIFT(name="drift_342", len=0.13104999999995925)
yo9_tq6 = KQUAD(name="yo9_tq6", len=0.75, k1=-0.006389820634304611)
drift_343 = DRIFT(name="drift_343", len=0.22622760000012931)
yo9_bh6 = MARKER(name="yo9_bh6")
drift_344 = DRIFT(name="drift_344", len=1.4007982592102053)
yo9_dh5 = ESBEND(name="yo9_dh5", len=8.698449375192075, angle=0.035863892964716426, e1=0.0, e2=0.0)
drift_345 = DRIFT(name="drift_345", len=1.9628259081100623)
yo9_tv5 = CORRECTOR(name="yo9_tv5", len=0.0, xkick=0.0, ykick=-0.0)
yo9_oct5 = thinMULTIPOLE(name="yo9_oct5", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yo9_dec5 = thinMULTIPOLE(name="yo9_dec5", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yo9_qs5 = thinMULTIPOLE(name="yo9_qs5", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_346 = DRIFT(name="drift_346", len=0.5452500000001237)
yo9_qd5 = KQUAD(name="yo9_qd5", len=1.11, k1=-0.09117978411883382)
drift_347 = DRIFT(name="drift_347", len=0.13104999999995925)
yo9_tq5 = KQUAD(name="yo9_tq5", len=0.75, k1=-0.00030909485194469355)
drift_348 = DRIFT(name="drift_348", len=0.22622760000012931)
yo9_bv5 = MARKER(name="yo9_bv5")
drift_349 = DRIFT(name="drift_349", len=4.3674723999997696)
yo9_th4 = CORRECTOR(name="yo9_th4", len=0.0, xkick=-0.0, ykick=0.0)
yo9_oct4 = thinMULTIPOLE(name="yo9_oct4", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yo9_dec4 = thinMULTIPOLE(name="yo9_dec4", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yo9_qs4 = thinMULTIPOLE(name="yo9_qs4", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_350 = DRIFT(name="drift_350", len=0.5442754999999124)
yo9_qf4 = KQUAD(name="yo9_qf4", len=1.811949, k1=0.09150512968949526)
drift_351 = DRIFT(name="drift_351", len=0.13007550000020274)
yo9_tq4 = KQUAD(name="yo9_tq4", len=0.75, k1=-0.012585780072451865)
drift_352 = DRIFT(name="drift_352", len=0.25163960000008956)
yo9_b4 = MARKER(name="yo9_b4")
drift_353 = DRIFT(name="drift_353", len=0.6552910000000338)
yo9_sv4 = MARKER(name="yo9_sv4")
drift_354 = DRIFT(name="drift_354", len=1.0178599999999278)
yo9_dmp3_2 = DRIFT(name="yo9_dmp3_2", len=4.873904)
yo9_dmp3_1 = DRIFT(name="yo9_dmp3_1", len=0.736295)
drift_355 = DRIFT(name="drift_355", len=7.531761000000188)
yo9_sv3_2 = MARKER(name="yo9_sv3_2")
drift_356 = DRIFT(name="drift_356", len=6.612789000000248)
yo9_b3_1 = MARKER(name="yo9_b3_1")
drift_357 = DRIFT(name="drift_357", len=4.816353000000163)
yo9_kfbh3 = CORRECTOR(name="yo9_kfbh3", len=0.0, xkick=0.0, ykick=0.0)
drift_358 = DRIFT(name="drift_358", len=1.1522344999998495)
yo9_ka3_5 = ESBEND(name="yo9_ka3_5", len=1.22, angle=0.0, e1=0.0, e2=0.0)
drift_359 = DRIFT(name="drift_359", len=0.13970000000017535)
yo9_ka3_4 = ESBEND(name="yo9_ka3_4", len=1.22, angle=0.0, e1=0.0, e2=0.0)
drift_360 = DRIFT(name="drift_360", len=0.13970000000017535)
yo9_ka3_3 = ESBEND(name="yo9_ka3_3", len=1.22, angle=0.0, e1=0.0, e2=0.0)
drift_361 = DRIFT(name="drift_361", len=0.13970000000017535)
yo9_ka3_2 = ESBEND(name="yo9_ka3_2", len=1.22, angle=0.0, e1=0.0, e2=0.0)
drift_362 = DRIFT(name="drift_362", len=0.13970000000017535)
yo9_ka3_1 = ESBEND(name="yo9_ka3_1", len=1.22, angle=0.0, e1=0.0, e2=0.0)
drift_363 = DRIFT(name="drift_363", len=0.630580923024354)
yo9_sv3_1 = MARKER(name="yo9_sv3_1")
drift_364 = DRIFT(name="drift_364", len=1.5878609999999753)
yo9_b3 = MARKER(name="yo9_b3")
drift_365 = DRIFT(name="drift_365", len=0.4714282599998114)
yo9_tv3 = CORRECTOR(name="yo9_tv3", len=0.0, xkick=0.0, ykick=-0.0)
yo9_sx3 = thinMULTIPOLE(name="yo9_sx3", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, -0.0, 0.0])
yo9_oct3 = thinMULTIPOLE(name="yo9_oct3", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yo9_dod3 = thinMULTIPOLE(name="yo9_dod3", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_366 = DRIFT(name="drift_366", len=0.5979680000000371)
yo9_qd3 = KQUAD(name="yo9_qd3", len=2.100484, k1=-0.05524805829658597)
drift_367 = DRIFT(name="drift_367", len=0.48977800000011484)
yo9_dods3 = thinMULTIPOLE(name="yo9_dods3", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yo9_octs3 = thinMULTIPOLE(name="yo9_octs3", PolynomA=[0.0, 0.0, 0.0, -0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yo9_sxs3 = thinMULTIPOLE(name="yo9_sxs3", PolynomA=[0.0, 0.0, -0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yo9_qs3 = thinMULTIPOLE(name="yo9_qs3", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_368 = DRIFT(name="drift_368", len=1.3709615000002486)
yo9_qf2 = KQUAD(name="yo9_qf2", len=3.391633, k1=0.056023079966336646)
drift_369 = DRIFT(name="drift_369", len=0.4942034999999123)
yo9_th2 = CORRECTOR(name="yo9_th2", len=0.0, xkick=-0.0, ykick=0.0)
yo9_oct2 = thinMULTIPOLE(name="yo9_oct2", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yo9_dec2 = thinMULTIPOLE(name="yo9_dec2", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yo9_dod2 = thinMULTIPOLE(name="yo9_dod2", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_370 = DRIFT(name="drift_370", len=1.1824280000000726)
yo9_qd1 = KQUAD(name="yo9_qd1", len=1.44, k1=-0.05822175644371885)
drift_371 = DRIFT(name="drift_371", len=0.3369447999998556)
yo9_b1 = MARKER(name="yo9_b1")
drift_372 = DRIFT(name="drift_372", len=0.5715743999999177)
cavity_591mhz = RFCA(name="cavity_591mhz", len=5.211, volt=0.0, freq=591155079.4074855, energy=275000000000.0, h=7560.0, lag=0.0)
drift_373 = DRIFT(name="drift_373", len=2.0)
dsw_ir10h = ESBEND(name="dsw_ir10h", len=6.000003375560012, angle=-0.003674538709966136, e1=-0.001837269354983068, e2=-0.001837269354983068)
drift_374 = DRIFT(name="drift_374", len=22.47907448557453)
dwarm_ir10h = ESBEND(name="dwarm_ir10h", len=6.000003375560012, angle=-0.003674538709966136, e1=-0.001837269354983068, e2=-0.001837269354983068)
drift_375 = DRIFT(name="drift_375", len=7.782574400000158)
bo10_b1 = MARKER(name="bo10_b1")
drift_376 = DRIFT(name="drift_376", len=0.3369447999998556)
bo10_qd1 = KQUAD(name="bo10_qd1", len=1.44, k1=-0.05822175644371885)
drift_377 = DRIFT(name="drift_377", len=1.1824280000000726)
bo10_th2 = CORRECTOR(name="bo10_th2", len=0.0, xkick=-0.0, ykick=0.0)
bo10_oct2 = thinMULTIPOLE(name="bo10_oct2", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
bo10_dec2 = thinMULTIPOLE(name="bo10_dec2", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
bo10_dod2 = thinMULTIPOLE(name="bo10_dod2", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_378 = DRIFT(name="drift_378", len=0.4942034999999123)
bo10_qf2 = KQUAD(name="bo10_qf2", len=3.391633, k1=0.056023079966336646)
drift_379 = DRIFT(name="drift_379", len=1.3709615000002486)
bo10_dods3 = thinMULTIPOLE(name="bo10_dods3", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
bo10_octs3 = thinMULTIPOLE(name="bo10_octs3", PolynomA=[0.0, 0.0, 0.0, -0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
bo10_sxs3 = thinMULTIPOLE(name="bo10_sxs3", PolynomA=[0.0, 0.0, -0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
bo10_qs3 = thinMULTIPOLE(name="bo10_qs3", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_380 = DRIFT(name="drift_380", len=0.48977800000011484)
bo10_qd3 = KQUAD(name="bo10_qd3", len=2.100484, k1=-0.05524805829658597)
drift_381 = DRIFT(name="drift_381", len=0.5979680000000371)
bo10_tv3 = CORRECTOR(name="bo10_tv3", len=0.0, xkick=0.0, ykick=-0.0)
bo10_sx3 = thinMULTIPOLE(name="bo10_sx3", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, -0.0, 0.0])
bo10_oct3 = thinMULTIPOLE(name="bo10_oct3", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
bo10_dod3 = thinMULTIPOLE(name="bo10_dod3", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_382 = DRIFT(name="drift_382", len=0.4714282599998114)
bo10_b3 = MARKER(name="bo10_b3")
drift_383 = DRIFT(name="drift_383", len=1.588039000000208)
bo10_sv3_1 = MARKER(name="bo10_sv3_1")
drift_384 = DRIFT(name="drift_384", len=0.6305815000000621)
bo10_ka3_1 = ESBEND(name="bo10_ka3_1", len=1.22, angle=0.0, e1=0.0, e2=0.0)
drift_385 = DRIFT(name="drift_385", len=0.13970000000017535)
bo10_ka3_2 = ESBEND(name="bo10_ka3_2", len=1.22, angle=0.0, e1=0.0, e2=0.0)
drift_386 = DRIFT(name="drift_386", len=0.13970000000017535)
bo10_ka3_3 = ESBEND(name="bo10_ka3_3", len=1.22, angle=0.0, e1=0.0, e2=0.0)
drift_387 = DRIFT(name="drift_387", len=0.13970000000017535)
bo10_ka3_4 = ESBEND(name="bo10_ka3_4", len=1.22, angle=0.0, e1=0.0, e2=0.0)
drift_388 = DRIFT(name="drift_388", len=0.13970000000017535)
bo10_ka3_5 = ESBEND(name="bo10_ka3_5", len=1.22, angle=0.0, e1=0.0, e2=0.0)
drift_389 = DRIFT(name="drift_389", len=1.2142355000000862)
bo10_kfbh3 = CORRECTOR(name="bo10_kfbh3", len=0.0, xkick=0.0, ykick=0.0)
drift_390 = DRIFT(name="drift_390", len=4.754352000000381)
bo10_b3_1 = MARKER(name="bo10_b3_1")
drift_391 = DRIFT(name="drift_391", len=6.612789000000248)
bo10_sv3_2 = MARKER(name="bo10_sv3_2")
drift_392 = DRIFT(name="drift_392", len=4.520692000000054)
bo10_c3 = DRIFT(name="bo10_c3", len=0.5334)
drift_393 = DRIFT(name="drift_393", len=2.4774910000001)
bo10_dmp3_1 = DRIFT(name="bo10_dmp3_1", len=0.736295)
bo10_dmp3_2 = DRIFT(name="bo10_dmp3_2", len=4.873904)
drift_394 = DRIFT(name="drift_394", len=1.0178604230241035)
bo10_sv4 = MARKER(name="bo10_sv4")
drift_395 = DRIFT(name="drift_395", len=0.65529000000015)
bo10_b4 = MARKER(name="bo10_b4")
drift_396 = DRIFT(name="drift_396", len=0.25163960000008956)
bo10_tq4 = KQUAD(name="bo10_tq4", len=0.75, k1=-0.012585780072451865)
drift_397 = DRIFT(name="drift_397", len=0.13007550000020274)
bo10_qf4 = KQUAD(name="bo10_qf4", len=1.811949, k1=0.09150512968949526)
drift_398 = DRIFT(name="drift_398", len=0.5442754999999124)
bo10_th4 = CORRECTOR(name="bo10_th4", len=0.0, xkick=-0.0, ykick=0.0)
bo10_oct4 = thinMULTIPOLE(name="bo10_oct4", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
bo10_dec4 = thinMULTIPOLE(name="bo10_dec4", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
bo10_qs4 = thinMULTIPOLE(name="bo10_qs4", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_399 = DRIFT(name="drift_399", len=4.3674723999997696)
bo10_bv5 = MARKER(name="bo10_bv5")
drift_400 = DRIFT(name="drift_400", len=0.22622760000012931)
bo10_tq5 = KQUAD(name="bo10_tq5", len=0.75, k1=-0.00030909485194469355)
drift_401 = DRIFT(name="drift_401", len=0.13104999999995925)
bo10_qd5 = KQUAD(name="bo10_qd5", len=1.11, k1=-0.09117978411883382)
drift_402 = DRIFT(name="drift_402", len=0.5452500000001237)
bo10_tv5 = CORRECTOR(name="bo10_tv5", len=0.0, xkick=0.0, ykick=-0.0)
bo10_oct5 = thinMULTIPOLE(name="bo10_oct5", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
bo10_dec5 = thinMULTIPOLE(name="bo10_dec5", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
bo10_qs5 = thinMULTIPOLE(name="bo10_qs5", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_403 = DRIFT(name="drift_403", len=1.9628259081100623)
bo10_dh5 = ESBEND(name="bo10_dh5", len=8.698449375192075, angle=0.035863892964716426, e1=0.0, e2=0.0)
drift_404 = DRIFT(name="drift_404", len=1.4007982592102053)
bo10_bh6 = MARKER(name="bo10_bh6")
drift_405 = DRIFT(name="drift_405", len=0.22622760000012931)
bo10_tq6 = KQUAD(name="bo10_tq6", len=0.75, k1=-0.006389820634304611)
drift_406 = DRIFT(name="drift_406", len=0.13104999999995925)
bo10_qf6 = KQUAD(name="bo10_qf6", len=1.11, k1=0.09117978411883382)
drift_407 = DRIFT(name="drift_407", len=0.5452500000001237)
bo10_th6 = CORRECTOR(name="bo10_th6", len=0.0, xkick=-0.0, ykick=0.0)
bo10_oct6 = thinMULTIPOLE(name="bo10_oct6", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
bo10_dec6 = thinMULTIPOLE(name="bo10_dec6", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
bo10_qgt6 = thinMULTIPOLE(name="bo10_qgt6", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, -0.0, 0.0, 0.0])
drift_408 = DRIFT(name="drift_408", len=8.082123248265361)
bo10_dh6 = ESBEND(name="bo10_dh6", len=2.94942686017, angle=0.012160550077129212, e1=0.0, e2=0.0)
drift_409 = DRIFT(name="drift_409", len=1.5844693382355217)
bo10_tv7 = CORRECTOR(name="bo10_tv7", len=0.0, xkick=0.0, ykick=-0.0)
bo10_oct7 = thinMULTIPOLE(name="bo10_oct7", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
bo10_dec7 = thinMULTIPOLE(name="bo10_dec7", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
bo10_qs7 = thinMULTIPOLE(name="bo10_qs7", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_410 = DRIFT(name="drift_410", len=0.5453779999998005)
bo10_qd7 = KQUAD(name="bo10_qd7", len=0.929744, k1=-0.08892562201533814)
drift_411 = DRIFT(name="drift_411", len=0.2962336000000505)
bo10_b7 = MARKER(name="bo10_b7")
drift_412 = DRIFT(name="drift_412", len=13.10883589569994)
bo10_b8 = MARKER(name="bo10_b8")
drift_413 = DRIFT(name="drift_413", len=0.2961055999999189)
bo10_qf8 = KQUAD(name="bo10_qf8", len=1.11, k1=0.07667711731989788)
drift_414 = DRIFT(name="drift_414", len=0.5452500000001237)
bo10_th8 = CORRECTOR(name="bo10_th8", len=0.0, xkick=-0.0, ykick=0.0)
bo10_oct8 = thinMULTIPOLE(name="bo10_oct8", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
bo10_dec8 = thinMULTIPOLE(name="bo10_dec8", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
bo10_qgt8 = thinMULTIPOLE(name="bo10_qgt8", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, -0.0, 0.0, 0.0])
drift_415 = DRIFT(name="drift_415", len=1.593705149733978)
bo10_dh8 = ESBEND(name="bo10_dh8", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_416 = DRIFT(name="drift_416", len=1.593705149733978)
bo10_tv9 = CORRECTOR(name="bo10_tv9", len=0.0, xkick=0.0, ykick=-0.0)
bo10_oct9 = thinMULTIPOLE(name="bo10_oct9", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
bo10_dec9 = thinMULTIPOLE(name="bo10_dec9", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
bo10_qs9 = thinMULTIPOLE(name="bo10_qs9", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_417 = DRIFT(name="drift_417", len=0.5452500000001237)
bo10_qd9 = KQUAD(name="bo10_qd9", len=1.11, k1=-0.07911855595309777)
drift_418 = DRIFT(name="drift_418", len=1.1072776000000886)
bo10_bv9 = MARKER(name="bo10_bv9")
drift_419 = DRIFT(name="drift_419", len=7.5200956482653964)
bo10_dh9 = ESBEND(name="bo10_dh9", len=2.94942686017, angle=0.012160550077129212, e1=0.0, e2=0.0)
drift_420 = DRIFT(name="drift_420", len=1.0224417382355568)
bo10_bh10 = MARKER(name="bo10_bh10")
drift_421 = DRIFT(name="drift_421", len=0.22622760000012931)
bo10_sxf10 = KSEXT(name="bo10_sxf10", len=0.75, k2=0.39909333141434805)
drift_422 = DRIFT(name="drift_422", len=0.13104999999995925)
bo10_qf10 = KQUAD(name="bo10_qf10", len=1.11, k1=0.08048071257602774)
drift_423 = DRIFT(name="drift_423", len=0.5452500000001237)
bo10_th10 = CORRECTOR(name="bo10_th10", len=0.0, xkick=-0.0, ykick=0.0)
bo10_oct10 = thinMULTIPOLE(name="bo10_oct10", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
bo10_dec10 = thinMULTIPOLE(name="bo10_dec10", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
bo10_qs10 = thinMULTIPOLE(name="bo10_qs10", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_424 = DRIFT(name="drift_424", len=1.593705149733978)
bo10_dh10 = ESBEND(name="bo10_dh10", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_425 = DRIFT(name="drift_425", len=1.0316775497340132)
bo10_bv11 = MARKER(name="bo10_bv11")
drift_426 = DRIFT(name="drift_426", len=0.22622760000012931)
bo10_sxd11 = KSEXT(name="bo10_sxd11", len=0.75, k2=-0.6290818079014571)
drift_427 = DRIFT(name="drift_427", len=0.13104999999995925)
bo10_qd11 = KQUAD(name="bo10_qd11", len=1.11, k1=-0.08364335555329494)
drift_428 = DRIFT(name="drift_428", len=0.5452500000001237)
bo10_tv11 = CORRECTOR(name="bo10_tv11", len=0.0, xkick=0.0, ykick=-0.0)
bo10_oct11 = thinMULTIPOLE(name="bo10_oct11", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
bo10_dec11 = thinMULTIPOLE(name="bo10_dec11", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
bo10_qs11 = thinMULTIPOLE(name="bo10_qs11", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_429 = DRIFT(name="drift_429", len=1.593705149733978)
bo10_dh11 = ESBEND(name="bo10_dh11", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_430 = DRIFT(name="drift_430", len=1.0316775497340132)
bo10_bh12 = MARKER(name="bo10_bh12")
drift_431 = DRIFT(name="drift_431", len=0.22622760000012931)
bo10_sxf12 = KSEXT(name="bo10_sxf12", len=0.75, k2=0.39909333141434805)
drift_432 = DRIFT(name="drift_432", len=0.13104999999995925)
bo10_qf12 = KQUAD(name="bo10_qf12", len=1.11, k1=0.08048071257602774)
drift_433 = DRIFT(name="drift_433", len=0.5452500000001237)
bo10_th12 = CORRECTOR(name="bo10_th12", len=0.0, xkick=-0.0, ykick=0.0)
bo10_oct12 = thinMULTIPOLE(name="bo10_oct12", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
bo10_dec12 = thinMULTIPOLE(name="bo10_dec12", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
bo10_qgt12 = thinMULTIPOLE(name="bo10_qgt12", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, -0.0, 0.0, 0.0])
drift_434 = DRIFT(name="drift_434", len=1.593705149733978)
bo10_dh12 = ESBEND(name="bo10_dh12", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_435 = DRIFT(name="drift_435", len=1.0316775497340132)
bo10_bv13 = MARKER(name="bo10_bv13")
drift_436 = DRIFT(name="drift_436", len=0.22622760000012931)
bo10_sxd13 = KSEXT(name="bo10_sxd13", len=0.75, k2=-0.6290818079014571)
drift_437 = DRIFT(name="drift_437", len=0.13104999999995925)
bo10_qd13 = KQUAD(name="bo10_qd13", len=1.11, k1=-0.08364335555329494)
drift_438 = DRIFT(name="drift_438", len=0.5452500000001237)
bo10_tv13 = CORRECTOR(name="bo10_tv13", len=0.0, xkick=0.0, ykick=-0.0)
bo10_oct13 = thinMULTIPOLE(name="bo10_oct13", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
bo10_dec13 = thinMULTIPOLE(name="bo10_dec13", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
bo10_qs13 = thinMULTIPOLE(name="bo10_qs13", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_439 = DRIFT(name="drift_439", len=1.593705149733978)
bo10_dh13 = ESBEND(name="bo10_dh13", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_440 = DRIFT(name="drift_440", len=1.0316775497340132)
bo10_bh14 = MARKER(name="bo10_bh14")
drift_441 = DRIFT(name="drift_441", len=0.22622760000012931)
bo10_sxf14 = KSEXT(name="bo10_sxf14", len=0.75, k2=0.39909333141434805)
drift_442 = DRIFT(name="drift_442", len=0.13104999999995925)
bo10_qf14 = KQUAD(name="bo10_qf14", len=1.11, k1=0.08048071257602774)
drift_443 = DRIFT(name="drift_443", len=0.5452500000001237)
bo10_th14 = CORRECTOR(name="bo10_th14", len=0.0, xkick=-0.0, ykick=0.0)
bo10_oct14 = thinMULTIPOLE(name="bo10_oct14", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
bo10_dec14 = thinMULTIPOLE(name="bo10_dec14", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
bo10_qgt14 = thinMULTIPOLE(name="bo10_qgt14", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, -0.0, 0.0, 0.0])
drift_444 = DRIFT(name="drift_444", len=1.593705149733978)
bo10_dh14 = ESBEND(name="bo10_dh14", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_445 = DRIFT(name="drift_445", len=1.0316775497340132)
bo10_bv15 = MARKER(name="bo10_bv15")
drift_446 = DRIFT(name="drift_446", len=0.22622760000012931)
bo10_sxd15 = KSEXT(name="bo10_sxd15", len=0.75, k2=-0.6290818079014571)
drift_447 = DRIFT(name="drift_447", len=0.13104999999995925)
bo10_qd15 = KQUAD(name="bo10_qd15", len=1.11, k1=-0.08364335555329494)
drift_448 = DRIFT(name="drift_448", len=0.5452500000001237)
bo10_tv15 = CORRECTOR(name="bo10_tv15", len=0.0, xkick=0.0, ykick=-0.0)
bo10_oct15 = thinMULTIPOLE(name="bo10_oct15", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
bo10_dec15 = thinMULTIPOLE(name="bo10_dec15", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
bo10_qs15 = thinMULTIPOLE(name="bo10_qs15", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_449 = DRIFT(name="drift_449", len=1.593705149733978)
bo10_dh15 = ESBEND(name="bo10_dh15", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_450 = DRIFT(name="drift_450", len=1.0316775497340132)
bo10_bh16 = MARKER(name="bo10_bh16")
drift_451 = DRIFT(name="drift_451", len=0.22622760000012931)
bo10_sxf16 = KSEXT(name="bo10_sxf16", len=0.75, k2=0.39909333141434805)
drift_452 = DRIFT(name="drift_452", len=0.13104999999995925)
bo10_qf16 = KQUAD(name="bo10_qf16", len=1.11, k1=0.08048071257602774)
drift_453 = DRIFT(name="drift_453", len=0.5452500000001237)
bo10_th16 = CORRECTOR(name="bo10_th16", len=0.0, xkick=-0.0, ykick=0.0)
bo10_oct16 = thinMULTIPOLE(name="bo10_oct16", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
bo10_dec16 = thinMULTIPOLE(name="bo10_dec16", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
bo10_qgt16 = thinMULTIPOLE(name="bo10_qgt16", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, -0.0, 0.0, 0.0])
drift_454 = DRIFT(name="drift_454", len=1.593705149733978)
bo10_dh16 = ESBEND(name="bo10_dh16", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_455 = DRIFT(name="drift_455", len=1.0316775497340132)
bo10_bv17 = MARKER(name="bo10_bv17")
drift_456 = DRIFT(name="drift_456", len=0.22622760000012931)
bo10_sxd17 = KSEXT(name="bo10_sxd17", len=0.75, k2=-0.6290818079014571)
drift_457 = DRIFT(name="drift_457", len=0.13104999999995925)
bo10_qd17 = KQUAD(name="bo10_qd17", len=1.11, k1=-0.08364335555329494)
drift_458 = DRIFT(name="drift_458", len=0.5452500000001237)
bo10_tv17 = CORRECTOR(name="bo10_tv17", len=0.0, xkick=0.0, ykick=-0.0)
bo10_oct17 = thinMULTIPOLE(name="bo10_oct17", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
bo10_dec17 = thinMULTIPOLE(name="bo10_dec17", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
bo10_qs17 = thinMULTIPOLE(name="bo10_qs17", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_459 = DRIFT(name="drift_459", len=1.593705149733978)
bo10_dh17 = ESBEND(name="bo10_dh17", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_460 = DRIFT(name="drift_460", len=1.0316775497340132)
bo10_bh18 = MARKER(name="bo10_bh18")
drift_461 = DRIFT(name="drift_461", len=0.22622760000012931)
bo10_sxf18 = KSEXT(name="bo10_sxf18", len=0.75, k2=0.39909333141434805)
drift_462 = DRIFT(name="drift_462", len=0.13104999999995925)
bo10_qf18 = KQUAD(name="bo10_qf18", len=1.11, k1=0.08048071257602774)
drift_463 = DRIFT(name="drift_463", len=0.5452500000001237)
bo10_th18 = CORRECTOR(name="bo10_th18", len=0.0, xkick=-0.0, ykick=0.0)
bo10_oct18 = thinMULTIPOLE(name="bo10_oct18", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
bo10_dec18 = thinMULTIPOLE(name="bo10_dec18", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
bo10_qgt18 = thinMULTIPOLE(name="bo10_qgt18", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, -0.0, 0.0, 0.0])
drift_464 = DRIFT(name="drift_464", len=1.593705149733978)
bo10_dh18 = ESBEND(name="bo10_dh18", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_465 = DRIFT(name="drift_465", len=1.0316775497340132)
bo10_bv19 = MARKER(name="bo10_bv19")
drift_466 = DRIFT(name="drift_466", len=0.22622760000012931)
bo10_sxd19 = KSEXT(name="bo10_sxd19", len=0.75, k2=-0.6290818079014571)
drift_467 = DRIFT(name="drift_467", len=0.13104999999995925)
bo10_qd19 = KQUAD(name="bo10_qd19", len=1.11, k1=-0.08364335555329494)
drift_468 = DRIFT(name="drift_468", len=0.2952500000001237)
bo10_tv19 = CORRECTOR(name="bo10_tv19", len=0.5, xkick=0.0, ykick=-0.0)
drift_469 = DRIFT(name="drift_469", len=1.343705149733978)
bo10_dh19 = ESBEND(name="bo10_dh19", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_470 = DRIFT(name="drift_470", len=1.0316775497340132)
bo10_bh20 = MARKER(name="bo10_bh20")
drift_471 = DRIFT(name="drift_471", len=0.22622760000012931)
bo10_sxf20 = KSEXT(name="bo10_sxf20", len=0.75, k2=0.39909333141434805)
drift_472 = DRIFT(name="drift_472", len=0.13104999999995925)
bo10_qf20 = KQUAD(name="bo10_qf20", len=1.11, k1=0.08048071257602774)
drift_473 = DRIFT(name="drift_473", len=0.2952500000001237)
bo10_th20 = CORRECTOR(name="bo10_th20", len=0.5, xkick=-0.0, ykick=0.0)
drift_474 = DRIFT(name="drift_474", len=1.343705149733978)
bo10_dh20 = ESBEND(name="bo10_dh20", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_475 = DRIFT(name="drift_475", len=1.0316775497340132)
bo11_bv21 = MARKER(name="bo11_bv21")
drift_476 = DRIFT(name="drift_476", len=0.22622760000012931)
bo11_sxd21 = KSEXT(name="bo11_sxd21", len=0.75, k2=-0.6290818079014571)
drift_477 = DRIFT(name="drift_477", len=0.13104999999995925)
bo11_qd21 = KQUAD(name="bo11_qd21", len=1.11, k1=-0.08364335555329494)
drift_478 = DRIFT(name="drift_478", len=0.2952500000001237)
bo11_tv21 = CORRECTOR(name="bo11_tv21", len=0.5, xkick=0.0, ykick=-0.0)
drift_479 = DRIFT(name="drift_479", len=1.343705149733978)
bo11_dh20 = ESBEND(name="bo11_dh20", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_480 = DRIFT(name="drift_480", len=1.0316775497340132)
bo11_bh20 = MARKER(name="bo11_bh20")
drift_481 = DRIFT(name="drift_481", len=0.22622760000012931)
bo11_sxf20 = KSEXT(name="bo11_sxf20", len=0.75, k2=0.39909333141434805)
drift_482 = DRIFT(name="drift_482", len=0.13104999999995925)
bo11_qf20 = KQUAD(name="bo11_qf20", len=1.11, k1=0.08048071257602774)
drift_483 = DRIFT(name="drift_483", len=0.2952500000001237)
bo11_th20 = CORRECTOR(name="bo11_th20", len=0.5, xkick=-0.0, ykick=0.0)
drift_484 = DRIFT(name="drift_484", len=1.343705149733978)
bo11_dh19 = ESBEND(name="bo11_dh19", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_485 = DRIFT(name="drift_485", len=1.0316775497340132)
bo11_bv19 = MARKER(name="bo11_bv19")
drift_486 = DRIFT(name="drift_486", len=0.22622760000012931)
bo11_sxd19 = KSEXT(name="bo11_sxd19", len=0.75, k2=-0.6290818079014571)
drift_487 = DRIFT(name="drift_487", len=0.13104999999995925)
bo11_qd19 = KQUAD(name="bo11_qd19", len=1.11, k1=-0.08364335555329494)
drift_488 = DRIFT(name="drift_488", len=0.2952500000001237)
bo11_tv19 = CORRECTOR(name="bo11_tv19", len=0.5, xkick=0.0, ykick=-0.0)
drift_489 = DRIFT(name="drift_489", len=1.343705149733978)
bo11_dh18 = ESBEND(name="bo11_dh18", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_490 = DRIFT(name="drift_490", len=1.0316775497340132)
bo11_bh18 = MARKER(name="bo11_bh18")
drift_491 = DRIFT(name="drift_491", len=0.22622760000012931)
bo11_sxf18 = KSEXT(name="bo11_sxf18", len=0.75, k2=0.39909333141434805)
drift_492 = DRIFT(name="drift_492", len=0.13104999999995925)
bo11_qf18 = KQUAD(name="bo11_qf18", len=1.11, k1=0.08048071257602774)
drift_493 = DRIFT(name="drift_493", len=0.2952500000001237)
bo11_th18 = CORRECTOR(name="bo11_th18", len=0.5, xkick=-0.0, ykick=0.0)
drift_494 = DRIFT(name="drift_494", len=1.343705149733978)
bo11_dh17 = ESBEND(name="bo11_dh17", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_495 = DRIFT(name="drift_495", len=1.0316775497340132)
bo11_bv17 = MARKER(name="bo11_bv17")
drift_496 = DRIFT(name="drift_496", len=0.22622760000012931)
bo11_sxd17 = KSEXT(name="bo11_sxd17", len=0.75, k2=-0.6290818079014571)
drift_497 = DRIFT(name="drift_497", len=0.13104999999995925)
bo11_qd17 = KQUAD(name="bo11_qd17", len=1.11, k1=-0.08364335555329494)
drift_498 = DRIFT(name="drift_498", len=0.2952500000001237)
bo11_tv17 = CORRECTOR(name="bo11_tv17", len=0.5, xkick=0.0, ykick=-0.0)
drift_499 = DRIFT(name="drift_499", len=1.343705149733978)
bo11_dh16 = ESBEND(name="bo11_dh16", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_500 = DRIFT(name="drift_500", len=1.0316775497340132)
bo11_bh16 = MARKER(name="bo11_bh16")
drift_501 = DRIFT(name="drift_501", len=0.22622760000012931)
bo11_sxf16 = KSEXT(name="bo11_sxf16", len=0.75, k2=0.39909333141434805)
drift_502 = DRIFT(name="drift_502", len=0.13104999999995925)
bo11_qf16 = KQUAD(name="bo11_qf16", len=1.11, k1=0.08048071257602774)
drift_503 = DRIFT(name="drift_503", len=0.2952500000001237)
bo11_th16 = CORRECTOR(name="bo11_th16", len=0.5, xkick=-0.0, ykick=0.0)
drift_504 = DRIFT(name="drift_504", len=1.343705149733978)
bo11_dh15 = ESBEND(name="bo11_dh15", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_505 = DRIFT(name="drift_505", len=1.0316775497340132)
bo11_bv15 = MARKER(name="bo11_bv15")
drift_506 = DRIFT(name="drift_506", len=0.22622760000012931)
bo11_sxd15 = KSEXT(name="bo11_sxd15", len=0.75, k2=-0.6290818079014571)
drift_507 = DRIFT(name="drift_507", len=0.13104999999995925)
bo11_qd15 = KQUAD(name="bo11_qd15", len=1.11, k1=-0.08364335555329494)
drift_508 = DRIFT(name="drift_508", len=0.2952500000001237)
bo11_tv15 = CORRECTOR(name="bo11_tv15", len=0.5, xkick=0.0, ykick=-0.0)
drift_509 = DRIFT(name="drift_509", len=1.343705149733978)
bo11_dh14 = ESBEND(name="bo11_dh14", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_510 = DRIFT(name="drift_510", len=1.0316775497340132)
bo11_bh14 = MARKER(name="bo11_bh14")
drift_511 = DRIFT(name="drift_511", len=0.22622760000012931)
bo11_sxf14 = KSEXT(name="bo11_sxf14", len=0.75, k2=0.39909333141434805)
drift_512 = DRIFT(name="drift_512", len=0.13104999999995925)
bo11_qf14 = KQUAD(name="bo11_qf14", len=1.11, k1=0.08048071257602774)
drift_513 = DRIFT(name="drift_513", len=0.2952500000001237)
bo11_th14 = CORRECTOR(name="bo11_th14", len=0.5, xkick=-0.0, ykick=0.0)
drift_514 = DRIFT(name="drift_514", len=1.343705149733978)
bo11_dh13 = ESBEND(name="bo11_dh13", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_515 = DRIFT(name="drift_515", len=1.0316775497340132)
bo11_bv13 = MARKER(name="bo11_bv13")
drift_516 = DRIFT(name="drift_516", len=0.22622760000012931)
bo11_sxd13 = KSEXT(name="bo11_sxd13", len=0.75, k2=-0.6290818079014571)
drift_517 = DRIFT(name="drift_517", len=0.13104999999995925)
bo11_qd13 = KQUAD(name="bo11_qd13", len=1.11, k1=-0.08364335555329494)
drift_518 = DRIFT(name="drift_518", len=0.2952500000001237)
bo11_tv13 = CORRECTOR(name="bo11_tv13", len=0.5, xkick=0.0, ykick=-0.0)
drift_519 = DRIFT(name="drift_519", len=1.343705149733978)
bo11_dh12 = ESBEND(name="bo11_dh12", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_520 = DRIFT(name="drift_520", len=1.0316775497340132)
bo11_bh12 = MARKER(name="bo11_bh12")
drift_521 = DRIFT(name="drift_521", len=0.22622760000012931)
bo11_sxf12 = KSEXT(name="bo11_sxf12", len=0.75, k2=0.39909333141434805)
drift_522 = DRIFT(name="drift_522", len=0.13104999999995925)
bo11_qf12 = KQUAD(name="bo11_qf12", len=1.11, k1=0.08048071257602774)
drift_523 = DRIFT(name="drift_523", len=0.2952500000001237)
bo11_th12 = CORRECTOR(name="bo11_th12", len=0.5, xkick=-0.0, ykick=0.0)
drift_524 = DRIFT(name="drift_524", len=1.343705149733978)
bo11_dh11 = ESBEND(name="bo11_dh11", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_525 = DRIFT(name="drift_525", len=1.0316775497340132)
bo11_bv11 = MARKER(name="bo11_bv11")
drift_526 = DRIFT(name="drift_526", len=0.22622760000012931)
bo11_sxd11 = KSEXT(name="bo11_sxd11", len=0.75, k2=-0.6290818079014571)
drift_527 = DRIFT(name="drift_527", len=0.13104999999995925)
bo11_qd11 = KQUAD(name="bo11_qd11", len=1.11, k1=-0.08364335555329494)
drift_528 = DRIFT(name="drift_528", len=0.2952500000001237)
bo11_tv11 = CORRECTOR(name="bo11_tv11", len=0.5, xkick=0.0, ykick=-0.0)
drift_529 = DRIFT(name="drift_529", len=1.343705149733978)
bo11_dh10 = ESBEND(name="bo11_dh10", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_530 = DRIFT(name="drift_530", len=1.0316775497340132)
bo11_bh10 = MARKER(name="bo11_bh10")
drift_531 = DRIFT(name="drift_531", len=0.22622760000012931)
bo11_sxf10 = KSEXT(name="bo11_sxf10", len=0.75, k2=0.39909333141434805)
drift_532 = DRIFT(name="drift_532", len=0.13104999999995925)
bo11_qf10 = KQUAD(name="bo11_qf10", len=1.11, k1=0.08048071257602774)
drift_533 = DRIFT(name="drift_533", len=0.5452500000001237)
bo11_th10 = CORRECTOR(name="bo11_th10", len=0.0, xkick=-0.0, ykick=0.0)
bo11_oct10 = thinMULTIPOLE(name="bo11_oct10", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
bo11_dec10 = thinMULTIPOLE(name="bo11_dec10", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
bo11_qs10 = thinMULTIPOLE(name="bo11_qs10", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_534 = DRIFT(name="drift_534", len=1.5844693382355217)
bo11_dh9 = ESBEND(name="bo11_dh9", len=2.94942686017, angle=0.012160550077129212, e1=0.0, e2=0.0)
drift_535 = DRIFT(name="drift_535", len=7.5200956482653964)
bo11_bv9 = MARKER(name="bo11_bv9")
drift_536 = DRIFT(name="drift_536", len=0.22622760000012931)
bo11_sxd9 = KSEXT(name="bo11_sxd9", len=0.75, k2=-0.6290818079014571)
drift_537 = DRIFT(name="drift_537", len=0.13104999999995925)
bo11_qd9 = KQUAD(name="bo11_qd9", len=1.11, k1=-0.08049896111799063)
drift_538 = DRIFT(name="drift_538", len=0.5452500000001237)
bo11_tv9 = CORRECTOR(name="bo11_tv9", len=0.0, xkick=0.0, ykick=-0.0)
bo11_oct9 = thinMULTIPOLE(name="bo11_oct9", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
bo11_dec9 = thinMULTIPOLE(name="bo11_dec9", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
bo11_qs9 = thinMULTIPOLE(name="bo11_qs9", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_539 = DRIFT(name="drift_539", len=1.593705149733978)
bo11_dh8 = ESBEND(name="bo11_dh8", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_540 = DRIFT(name="drift_540", len=1.593705149733978)
bo11_th8 = CORRECTOR(name="bo11_th8", len=0.0, xkick=-0.0, ykick=0.0)
bo11_oct8 = thinMULTIPOLE(name="bo11_oct8", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
bo11_dec8 = thinMULTIPOLE(name="bo11_dec8", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
bo11_qgt8 = thinMULTIPOLE(name="bo11_qgt8", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, -0.0, 0.0, 0.0])
drift_541 = DRIFT(name="drift_541", len=0.5452500000001237)
bo11_qf8 = KQUAD(name="bo11_qf8", len=1.11, k1=0.07952238566471066)
drift_542 = DRIFT(name="drift_542", len=0.2961055999999189)
bo11_b8 = MARKER(name="bo11_b8")
drift_543 = DRIFT(name="drift_543", len=1.3184179478498663)
bo11_snk_hlx4 = DRIFT(name="bo11_snk_hlx4", len=2.4)
drift_544 = DRIFT(name="drift_544", len=0.21199999999998909)
bo11_snk_hlx3 = DRIFT(name="bo11_snk_hlx3", len=2.4)
drift_545 = DRIFT(name="drift_545", len=0.22400000000016007)
bo11_bsnk = MARKER(name="bo11_bsnk")
drift_546 = DRIFT(name="drift_546", len=0.22400000000016007)
bo11_snk_hlx2 = DRIFT(name="bo11_snk_hlx2", len=2.4)
drift_547 = DRIFT(name="drift_547", len=0.21199999999998909)
bo11_snk_hlx1 = DRIFT(name="bo11_snk_hlx1", len=2.4)
drift_548 = DRIFT(name="drift_548", len=1.3184179478498663)
bo11_b7 = MARKER(name="bo11_b7")
drift_549 = DRIFT(name="drift_549", len=0.2962336000000505)
bo11_qd7 = KQUAD(name="bo11_qd7", len=0.929744, k1=-0.08736055902386013)
drift_550 = DRIFT(name="drift_550", len=0.5453779999998005)
bo11_tv7 = CORRECTOR(name="bo11_tv7", len=0.0, xkick=0.0, ykick=-0.0)
bo11_oct7 = thinMULTIPOLE(name="bo11_oct7", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
bo11_dec7 = thinMULTIPOLE(name="bo11_dec7", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
bo11_qs7 = thinMULTIPOLE(name="bo11_qs7", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_551 = DRIFT(name="drift_551", len=1.5844693382355217)
bo11_dh6 = ESBEND(name="bo11_dh6", len=2.94942686017, angle=0.012160550077129212, e1=0.0, e2=0.0)
drift_552 = DRIFT(name="drift_552", len=8.082123248265361)
bo11_th6 = CORRECTOR(name="bo11_th6", len=0.0, xkick=-0.0, ykick=0.0)
bo11_oct6 = thinMULTIPOLE(name="bo11_oct6", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
bo11_dec6 = thinMULTIPOLE(name="bo11_dec6", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
bo11_qgt6 = thinMULTIPOLE(name="bo11_qgt6", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, -0.0, 0.0, 0.0])
drift_553 = DRIFT(name="drift_553", len=0.5452500000001237)
bo11_qf6 = KQUAD(name="bo11_qf6", len=1.11, k1=0.08961372028016185)
drift_554 = DRIFT(name="drift_554", len=0.13104999999995925)
bo11_tq6 = KQUAD(name="bo11_tq6", len=0.75, k1=0.023124629567031638)
drift_555 = DRIFT(name="drift_555", len=0.22622760000012931)
bo11_bh6 = MARKER(name="bo11_bh6")
drift_556 = DRIFT(name="drift_556", len=1.4007982592102053)
bo11_dh5 = ESBEND(name="bo11_dh5", len=8.698449375192075, angle=0.035863892964716426, e1=0.0, e2=0.0)
drift_557 = DRIFT(name="drift_557", len=1.9628259081100623)
bo11_tv5 = CORRECTOR(name="bo11_tv5", len=0.0, xkick=0.0, ykick=-0.0)
bo11_oct5 = thinMULTIPOLE(name="bo11_oct5", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
bo11_dec5 = thinMULTIPOLE(name="bo11_dec5", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
bo11_qs5 = thinMULTIPOLE(name="bo11_qs5", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_558 = DRIFT(name="drift_558", len=0.5452500000001237)
bo11_qd5 = KQUAD(name="bo11_qd5", len=1.11, k1=-0.08961372028016185)
drift_559 = DRIFT(name="drift_559", len=0.13104999999995925)
bo11_tq5 = KQUAD(name="bo11_tq5", len=0.75, k1=-0.016089312473386432)
drift_560 = DRIFT(name="drift_560", len=0.22622760000012931)
bo11_bv5 = MARKER(name="bo11_bv5")
drift_561 = DRIFT(name="drift_561", len=4.3674723999997696)
bo11_th4 = CORRECTOR(name="bo11_th4", len=0.0, xkick=-0.0, ykick=0.0)
bo11_oct4 = thinMULTIPOLE(name="bo11_oct4", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
bo11_dec4 = thinMULTIPOLE(name="bo11_dec4", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
bo11_qs4 = thinMULTIPOLE(name="bo11_qs4", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_562 = DRIFT(name="drift_562", len=0.5442754999999124)
bo11_qf4 = KQUAD(name="bo11_qf4", len=1.811949, k1=0.08993347785850449)
drift_563 = DRIFT(name="drift_563", len=0.13007550000020274)
bo11_tq4 = KQUAD(name="bo11_tq4", len=0.75, k1=-0.009171569727420559)
drift_564 = DRIFT(name="drift_564", len=0.25163960000008956)
bo11_b4 = MARKER(name="bo11_b4")
drift_565 = DRIFT(name="drift_565", len=0.655297000000246)
bo11_sv4 = MARKER(name="bo11_sv4")
drift_566 = DRIFT(name="drift_566", len=6.938098999999966)
bo11_mskh3 = MARKER(name="bo11_mskh3")
drift_567 = DRIFT(name="drift_567", len=26.28698200000008)
scol_h3_2 = MARKER(name="scol_h3_2")
drift_568 = DRIFT(name="drift_568", len=0.5)
bo11_kfbh3 = CORRECTOR(name="bo11_kfbh3", len=0.0, xkick=0.0, ykick=0.0)
drift_569 = DRIFT(name="drift_569", len=0.3055544230242049)
bo11_sv3 = MARKER(name="bo11_sv3")
drift_570 = DRIFT(name="drift_570", len=1.587797000000137)
bo11_b3 = MARKER(name="bo11_b3")
drift_571 = DRIFT(name="drift_571", len=0.4714282599998114)
bo11_tv3 = CORRECTOR(name="bo11_tv3", len=0.0, xkick=0.0, ykick=-0.0)
bo11_sx3 = thinMULTIPOLE(name="bo11_sx3", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, -0.0, 0.0])
bo11_oct3 = thinMULTIPOLE(name="bo11_oct3", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
bo11_dod3 = thinMULTIPOLE(name="bo11_dod3", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_572 = DRIFT(name="drift_572", len=0.5979680000000371)
bo11_qd3 = KQUAD(name="bo11_qd3", len=2.100484, k1=-0.05485242582382659)
drift_573 = DRIFT(name="drift_573", len=0.48977800000011484)
bo11_dods3 = thinMULTIPOLE(name="bo11_dods3", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
bo11_octs3 = thinMULTIPOLE(name="bo11_octs3", PolynomA=[0.0, 0.0, 0.0, -0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
bo11_sxs3 = thinMULTIPOLE(name="bo11_sxs3", PolynomA=[0.0, 0.0, -0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
bo11_qs3 = thinMULTIPOLE(name="bo11_qs3", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_574 = DRIFT(name="drift_574", len=1.3709615000002486)
bo11_qf2 = KQUAD(name="bo11_qf2", len=3.391633, k1=0.05621963315561761)
drift_575 = DRIFT(name="drift_575", len=0.4942034999999123)
bo11_th2 = CORRECTOR(name="bo11_th2", len=0.0, xkick=-0.0, ykick=0.0)
bo11_oct2 = thinMULTIPOLE(name="bo11_oct2", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
bo11_dec2 = thinMULTIPOLE(name="bo11_dec2", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
bo11_dod2 = thinMULTIPOLE(name="bo11_dod2", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_576 = DRIFT(name="drift_576", len=1.1824280000000726)
bo11_qd1 = KQUAD(name="bo11_qd1", len=1.44, k1=-0.05706628739309596)
drift_577 = DRIFT(name="drift_577", len=0.3369447999998556)
bo11_b1 = MARKER(name="bo11_b1")
drift_578 = DRIFT(name="drift_578", len=1.6630552000001444)
dwarm_ir12h = ESBEND(name="dwarm_ir12h", len=6.000008532300472, angle=0.005842017452824684, e1=0.002921008726412342, e2=0.002921008726412342)
drift_579 = DRIFT(name="drift_579", len=34.719843198343824)
dsw_ir12h = ESBEND(name="dsw_ir12h", len=6.000008532300472, angle=-0.005842017452824684, e1=-0.002921008726412342, e2=-0.002921008726412342)
drift_580 = DRIFT(name="drift_580", len=1.6630552000001444)
bi12_b1 = MARKER(name="bi12_b1")
drift_581 = DRIFT(name="drift_581", len=0.3369447999998556)
bi12_qf1 = KQUAD(name="bi12_qf1", len=1.44, k1=0.055778910124163486)
drift_582 = DRIFT(name="drift_582", len=1.1824280000000726)
bi12_tv2 = CORRECTOR(name="bi12_tv2", len=0.0, xkick=0.0, ykick=-0.0)
bi12_oct2 = thinMULTIPOLE(name="bi12_oct2", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
bi12_dec2 = thinMULTIPOLE(name="bi12_dec2", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
bi12_dod2 = thinMULTIPOLE(name="bi12_dod2", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_583 = DRIFT(name="drift_583", len=0.4942034999999123)
bi12_qd2 = KQUAD(name="bi12_qd2", len=3.391633, k1=-0.05493732596646653)
drift_584 = DRIFT(name="drift_584", len=1.3709615000002486)
bi12_dods3 = thinMULTIPOLE(name="bi12_dods3", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
bi12_octs3 = thinMULTIPOLE(name="bi12_octs3", PolynomA=[0.0, 0.0, 0.0, -0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
bi12_sxs3 = thinMULTIPOLE(name="bi12_sxs3", PolynomA=[0.0, 0.0, -0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
bi12_qs3 = thinMULTIPOLE(name="bi12_qs3", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_585 = DRIFT(name="drift_585", len=0.48977800000011484)
bi12_qf3 = KQUAD(name="bi12_qf3", len=2.100484, k1=0.05424921787784084)
drift_586 = DRIFT(name="drift_586", len=0.5979680000000371)
bi12_th3 = CORRECTOR(name="bi12_th3", len=0.0, xkick=-0.0, ykick=0.0)
bi12_sx3 = thinMULTIPOLE(name="bi12_sx3", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, -0.0, 0.0])
bi12_oct3 = thinMULTIPOLE(name="bi12_oct3", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
bi12_dod3 = thinMULTIPOLE(name="bi12_dod3", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_587 = DRIFT(name="drift_587", len=0.4714282599998114)
bi12_b3 = MARKER(name="bi12_b3")
drift_588 = DRIFT(name="drift_588", len=1.5873540000002322)
bi12_sv3_1 = MARKER(name="bi12_sv3_1")
drift_589 = DRIFT(name="drift_589", len=0.3056169999999838)
bi12_kfbh3 = CORRECTOR(name="bi12_kfbh3", len=0.0, xkick=0.0, ykick=0.0)
drift_590 = DRIFT(name="drift_590", len=4.341284000000087)
bi12_ipm3 = MARKER(name="bi12_ipm3")
drift_591 = DRIFT(name="drift_591", len=3.098467000000255)
scol_h3_1 = MARKER(name="scol_h3_1")
drift_592 = DRIFT(name="drift_592", len=3.6375230000003285)
bi12_eld3 = MARKER(name="bi12_eld3")
drift_593 = DRIFT(name="drift_593", len=14.32082500000024)
bi12_ksch3_1 = MARKER(name="bi12_ksch3_1")
drift_594 = DRIFT(name="drift_594", len=1.2616179999999986)
bi12_ksch3_2 = MARKER(name="bi12_ksch3_2")
drift_595 = DRIFT(name="drift_595", len=0.7800339999998869)
pcol_h4 = MARKER(name="pcol_h4")
drift_596 = DRIFT(name="drift_596", len=0.7800339999998869)
bi12_kscv3_1 = MARKER(name="bi12_kscv3_1")
drift_597 = DRIFT(name="drift_597", len=1.1790679999999156)
bi12_kscv3_2 = MARKER(name="bi12_kscv3_2")
drift_598 = DRIFT(name="drift_598", len=1.6361219999998866)
bi12_sv3_2 = MARKER(name="bi12_sv3_2")
drift_599 = DRIFT(name="drift_599", len=0.8097349999998187)
bi12_pol3_1 = MARKER(name="bi12_pol3_1")
drift_600 = DRIFT(name="drift_600", len=0.5079999999998108)
bi12_pol3_2 = MARKER(name="bi12_pol3_2")
drift_601 = DRIFT(name="drift_601", len=1.3727134230243792)
bi12_sv4 = MARKER(name="bi12_sv4")
drift_602 = DRIFT(name="drift_602", len=0.6553349999999227)
bi12_b4 = MARKER(name="bi12_b4")
drift_603 = DRIFT(name="drift_603", len=0.25163960000008956)
bi12_tq4 = KQUAD(name="bi12_tq4", len=0.75, k1=0.03462226952510231)
drift_604 = DRIFT(name="drift_604", len=0.13007550000020274)
bi12_qd4 = KQUAD(name="bi12_qd4", len=1.811949, k1=-0.08993347785850449)
drift_605 = DRIFT(name="drift_605", len=0.5442754999999124)
bi12_tv4 = CORRECTOR(name="bi12_tv4", len=0.0, xkick=0.0, ykick=-0.0)
bi12_oct4 = thinMULTIPOLE(name="bi12_oct4", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
bi12_dec4 = thinMULTIPOLE(name="bi12_dec4", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
bi12_qs4 = thinMULTIPOLE(name="bi12_qs4", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_606 = DRIFT(name="drift_606", len=4.3674723999997696)
bi12_bh5 = MARKER(name="bi12_bh5")
drift_607 = DRIFT(name="drift_607", len=0.22622760000012931)
bi12_tq5 = KQUAD(name="bi12_tq5", len=0.75, k1=0.003959246604291099)
drift_608 = DRIFT(name="drift_608", len=0.13104999999995925)
bi12_qf5 = KQUAD(name="bi12_qf5", len=1.11, k1=0.08961372028016185)
drift_609 = DRIFT(name="drift_609", len=0.5452500000001237)
bi12_th5 = CORRECTOR(name="bi12_th5", len=0.0, xkick=-0.0, ykick=0.0)
bi12_oct5 = thinMULTIPOLE(name="bi12_oct5", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
bi12_dec5 = thinMULTIPOLE(name="bi12_dec5", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
bi12_qgt5 = thinMULTIPOLE(name="bi12_qgt5", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, -0.0, 0.0, 0.0])
drift_610 = DRIFT(name="drift_610", len=2.8397960027023146)
bi12_dh5 = ESBEND(name="bi12_dh5", len=6.915999880527922, angle=0.028514815544784158, e1=0.0, e2=0.0)
drift_611 = DRIFT(name="drift_611", len=2.2777684516026966)
bi12_bv6 = MARKER(name="bi12_bv6")
drift_612 = DRIFT(name="drift_612", len=0.22622760000012931)
bi12_tq6 = KQUAD(name="bi12_tq6", len=0.75, k1=-0.0023796899537632554)
drift_613 = DRIFT(name="drift_613", len=0.13104999999995925)
bi12_qd6 = KQUAD(name="bi12_qd6", len=1.11, k1=-0.08961372028016185)
drift_614 = DRIFT(name="drift_614", len=0.5452500000001237)
bi12_tv6 = CORRECTOR(name="bi12_tv6", len=0.0, xkick=0.0, ykick=-0.0)
bi12_oct6 = thinMULTIPOLE(name="bi12_oct6", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
bi12_dec6 = thinMULTIPOLE(name="bi12_dec6", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
bi12_qs6 = thinMULTIPOLE(name="bi12_qs6", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_615 = DRIFT(name="drift_615", len=8.076650933293877)
bi12_dh6 = ESBEND(name="bi12_dh6", len=2.94942686017, angle=0.012160550077129212, e1=0.0, e2=0.0)
drift_616 = DRIFT(name="drift_616", len=1.5789970232640371)
bi12_th7 = CORRECTOR(name="bi12_th7", len=0.0, xkick=-0.0, ykick=0.0)
bi12_oct7 = thinMULTIPOLE(name="bi12_oct7", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
bi12_dec7 = thinMULTIPOLE(name="bi12_dec7", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
bi12_qgt7 = thinMULTIPOLE(name="bi12_qgt7", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, -0.0, 0.0, 0.0])
drift_617 = DRIFT(name="drift_617", len=0.5453779999998005)
bi12_qf7 = KQUAD(name="bi12_qf7", len=0.929744, k1=0.08736055902386013)
drift_618 = DRIFT(name="drift_618", len=0.2962336000000505)
bi12_b7 = MARKER(name="bi12_b7")
drift_619 = DRIFT(name="drift_619", len=13.10883589569994)
bi12_b8 = MARKER(name="bi12_b8")
drift_620 = DRIFT(name="drift_620", len=0.2961055999999189)
bi12_qd8 = KQUAD(name="bi12_qd8", len=1.11, k1=-0.08049896111799063)
drift_621 = DRIFT(name="drift_621", len=0.5452500000001237)
bi12_tv8 = CORRECTOR(name="bi12_tv8", len=0.0, xkick=0.0, ykick=-0.0)
bi12_oct8 = thinMULTIPOLE(name="bi12_oct8", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
bi12_dec8 = thinMULTIPOLE(name="bi12_dec8", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
bi12_qs8 = thinMULTIPOLE(name="bi12_qs8", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_622 = DRIFT(name="drift_622", len=1.5761871258662268)
bi12_dh8 = ESBEND(name="bi12_dh8", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_623 = DRIFT(name="drift_623", len=1.5761871258664542)
bi12_th9 = CORRECTOR(name="bi12_th9", len=0.0, xkick=-0.0, ykick=0.0)
bi12_oct9 = thinMULTIPOLE(name="bi12_oct9", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
bi12_dec9 = thinMULTIPOLE(name="bi12_dec9", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
bi12_qs9 = thinMULTIPOLE(name="bi12_qs9", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_624 = DRIFT(name="drift_624", len=0.5452500000001237)
bi12_qf9 = KQUAD(name="bi12_qf9", len=1.11, k1=0.07952238566471066)
drift_625 = DRIFT(name="drift_625", len=0.13104999999995925)
bi12_sxf9 = KSEXT(name="bi12_sxf9", len=0.75, k2=0.39909333141434805)
drift_626 = DRIFT(name="drift_626", len=0.22622760000012931)
bi12_bh9 = MARKER(name="bi12_bh9")
drift_627 = DRIFT(name="drift_627", len=7.514623333294367)
bi12_dh9 = ESBEND(name="bi12_dh9", len=2.94942686017, angle=0.012160550077129212, e1=0.0, e2=0.0)
drift_628 = DRIFT(name="drift_628", len=1.5789970232644919)
bi12_tv10 = CORRECTOR(name="bi12_tv10", len=0.0, xkick=0.0, ykick=-0.0)
bi12_oct10 = thinMULTIPOLE(name="bi12_oct10", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
bi12_dec10 = thinMULTIPOLE(name="bi12_dec10", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
bi12_qs10 = thinMULTIPOLE(name="bi12_qs10", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_629 = DRIFT(name="drift_629", len=0.5452500000001237)
bi12_qd10 = KQUAD(name="bi12_qd10", len=1.11, k1=-0.08364335555329494)
drift_630 = DRIFT(name="drift_630", len=0.13104999999995925)
bi12_sxd10 = KSEXT(name="bi12_sxd10", len=0.75, k2=-0.6290818079014571)
drift_631 = DRIFT(name="drift_631", len=0.22622760000012931)
bi12_bv10 = MARKER(name="bi12_bv10")
drift_632 = DRIFT(name="drift_632", len=1.0141595258655798)
bi12_dh10 = ESBEND(name="bi12_dh10", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_633 = DRIFT(name="drift_633", len=1.5761871258664542)
bi12_th11 = CORRECTOR(name="bi12_th11", len=0.0, xkick=-0.0, ykick=0.0)
bi12_oct11 = thinMULTIPOLE(name="bi12_oct11", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
bi12_dec11 = thinMULTIPOLE(name="bi12_dec11", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
bi12_qgt11 = thinMULTIPOLE(name="bi12_qgt11", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, -0.0, 0.0, 0.0])
drift_634 = DRIFT(name="drift_634", len=0.5452500000001237)
bi12_qf11 = KQUAD(name="bi12_qf11", len=1.11, k1=0.08048071257602774)
drift_635 = DRIFT(name="drift_635", len=0.13104999999995925)
bi12_sxf11 = KSEXT(name="bi12_sxf11", len=0.75, k2=0.39909333141434805)
drift_636 = DRIFT(name="drift_636", len=0.22622760000012931)
bi12_bh11 = MARKER(name="bi12_bh11")
drift_637 = DRIFT(name="drift_637", len=1.0141595258655798)
bi12_dh11 = ESBEND(name="bi12_dh11", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_638 = DRIFT(name="drift_638", len=1.5761871258664542)
bi12_tv12 = CORRECTOR(name="bi12_tv12", len=0.0, xkick=0.0, ykick=-0.0)
bi12_oct12 = thinMULTIPOLE(name="bi12_oct12", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
bi12_dec12 = thinMULTIPOLE(name="bi12_dec12", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
bi12_qs12 = thinMULTIPOLE(name="bi12_qs12", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_639 = DRIFT(name="drift_639", len=0.5452500000001237)
bi12_qd12 = KQUAD(name="bi12_qd12", len=1.11, k1=-0.08364335555329494)
drift_640 = DRIFT(name="drift_640", len=0.13104999999995925)
bi12_sxd12 = KSEXT(name="bi12_sxd12", len=0.75, k2=-0.6290818079014571)
drift_641 = DRIFT(name="drift_641", len=0.22622760000012931)
bi12_bv12 = MARKER(name="bi12_bv12")
drift_642 = DRIFT(name="drift_642", len=1.0141595258655798)
bi12_dh12 = ESBEND(name="bi12_dh12", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_643 = DRIFT(name="drift_643", len=1.5761871258664542)
bi12_th13 = CORRECTOR(name="bi12_th13", len=0.0, xkick=-0.0, ykick=0.0)
bi12_oct13 = thinMULTIPOLE(name="bi12_oct13", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
bi12_dec13 = thinMULTIPOLE(name="bi12_dec13", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
bi12_qgt13 = thinMULTIPOLE(name="bi12_qgt13", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, -0.0, 0.0, 0.0])
drift_644 = DRIFT(name="drift_644", len=0.5452500000001237)
bi12_qf13 = KQUAD(name="bi12_qf13", len=1.11, k1=0.08048071257602774)
drift_645 = DRIFT(name="drift_645", len=0.13104999999995925)
bi12_sxf13 = KSEXT(name="bi12_sxf13", len=0.75, k2=0.39909333141434805)
drift_646 = DRIFT(name="drift_646", len=0.22622760000012931)
bi12_bh13 = MARKER(name="bi12_bh13")
drift_647 = DRIFT(name="drift_647", len=1.0141595258655798)
bi12_dh13 = ESBEND(name="bi12_dh13", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_648 = DRIFT(name="drift_648", len=1.5761871258664542)
bi12_tv14 = CORRECTOR(name="bi12_tv14", len=0.0, xkick=0.0, ykick=-0.0)
bi12_oct14 = thinMULTIPOLE(name="bi12_oct14", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
bi12_dec14 = thinMULTIPOLE(name="bi12_dec14", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
bi12_qs14 = thinMULTIPOLE(name="bi12_qs14", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_649 = DRIFT(name="drift_649", len=0.5452500000001237)
bi12_qd14 = KQUAD(name="bi12_qd14", len=1.11, k1=-0.08364335555329494)
drift_650 = DRIFT(name="drift_650", len=0.13104999999995925)
bi12_sxd14 = KSEXT(name="bi12_sxd14", len=0.75, k2=-0.6290818079014571)
drift_651 = DRIFT(name="drift_651", len=0.22622760000012931)
bi12_bv14 = MARKER(name="bi12_bv14")
drift_652 = DRIFT(name="drift_652", len=1.0141595258655798)
bi12_dh14 = ESBEND(name="bi12_dh14", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_653 = DRIFT(name="drift_653", len=1.5761871258664542)
bi12_th15 = CORRECTOR(name="bi12_th15", len=0.0, xkick=-0.0, ykick=0.0)
bi12_oct15 = thinMULTIPOLE(name="bi12_oct15", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
bi12_dec15 = thinMULTIPOLE(name="bi12_dec15", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
bi12_qgt15 = thinMULTIPOLE(name="bi12_qgt15", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, -0.0, 0.0, 0.0])
drift_654 = DRIFT(name="drift_654", len=0.5452500000001237)
bi12_qf15 = KQUAD(name="bi12_qf15", len=1.11, k1=0.08048071257602774)
drift_655 = DRIFT(name="drift_655", len=0.13104999999995925)
bi12_sxf15 = KSEXT(name="bi12_sxf15", len=0.75, k2=0.39909333141434805)
drift_656 = DRIFT(name="drift_656", len=0.22622760000012931)
bi12_bh15 = MARKER(name="bi12_bh15")
drift_657 = DRIFT(name="drift_657", len=1.0141595258655798)
bi12_dh15 = ESBEND(name="bi12_dh15", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_658 = DRIFT(name="drift_658", len=1.5761871258664542)
bi12_tv16 = CORRECTOR(name="bi12_tv16", len=0.0, xkick=0.0, ykick=-0.0)
bi12_oct16 = thinMULTIPOLE(name="bi12_oct16", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
bi12_dec16 = thinMULTIPOLE(name="bi12_dec16", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
bi12_qs16 = thinMULTIPOLE(name="bi12_qs16", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_659 = DRIFT(name="drift_659", len=0.5452500000001237)
bi12_qd16 = KQUAD(name="bi12_qd16", len=1.11, k1=-0.08364335555329494)
drift_660 = DRIFT(name="drift_660", len=0.13104999999995925)
bi12_sxd16 = KSEXT(name="bi12_sxd16", len=0.75, k2=-0.6290818079014571)
drift_661 = DRIFT(name="drift_661", len=0.22622760000012931)
bi12_bv16 = MARKER(name="bi12_bv16")
drift_662 = DRIFT(name="drift_662", len=1.0141595258655798)
bi12_dh16 = ESBEND(name="bi12_dh16", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_663 = DRIFT(name="drift_663", len=1.5761871258664542)
bi12_th17 = CORRECTOR(name="bi12_th17", len=0.0, xkick=-0.0, ykick=0.0)
bi12_oct17 = thinMULTIPOLE(name="bi12_oct17", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
bi12_dec17 = thinMULTIPOLE(name="bi12_dec17", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
bi12_qgt17 = thinMULTIPOLE(name="bi12_qgt17", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, -0.0, 0.0, 0.0])
drift_664 = DRIFT(name="drift_664", len=0.5452500000001237)
bi12_qf17 = KQUAD(name="bi12_qf17", len=1.11, k1=0.08048071257602774)
drift_665 = DRIFT(name="drift_665", len=0.13104999999995925)
bi12_sxf17 = KSEXT(name="bi12_sxf17", len=0.75, k2=0.39909333141434805)
drift_666 = DRIFT(name="drift_666", len=0.22622760000012931)
bi12_bh17 = MARKER(name="bi12_bh17")
drift_667 = DRIFT(name="drift_667", len=1.0141595258655798)
bi12_dh17 = ESBEND(name="bi12_dh17", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_668 = DRIFT(name="drift_668", len=1.5761871258664542)
bi12_tv18 = CORRECTOR(name="bi12_tv18", len=0.0, xkick=0.0, ykick=-0.0)
bi12_oct18 = thinMULTIPOLE(name="bi12_oct18", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
bi12_dec18 = thinMULTIPOLE(name="bi12_dec18", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
bi12_qs18 = thinMULTIPOLE(name="bi12_qs18", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_669 = DRIFT(name="drift_669", len=0.5452500000001237)
bi12_qd18 = KQUAD(name="bi12_qd18", len=1.11, k1=-0.08364335555329494)
drift_670 = DRIFT(name="drift_670", len=0.13104999999995925)
bi12_sxd18 = KSEXT(name="bi12_sxd18", len=0.75, k2=-0.6290818079014571)
drift_671 = DRIFT(name="drift_671", len=0.22622760000012931)
bi12_bv18 = MARKER(name="bi12_bv18")
drift_672 = DRIFT(name="drift_672", len=1.0141595258655798)
bi12_dh18 = ESBEND(name="bi12_dh18", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_673 = DRIFT(name="drift_673", len=1.3261871258664542)
bi12_th19 = CORRECTOR(name="bi12_th19", len=0.5, xkick=-0.0, ykick=0.0)
drift_674 = DRIFT(name="drift_674", len=0.2952500000001237)
bi12_qf19 = KQUAD(name="bi12_qf19", len=1.11, k1=0.08048071257602774)
drift_675 = DRIFT(name="drift_675", len=0.13104999999995925)
bi12_sxf19 = KSEXT(name="bi12_sxf19", len=0.75, k2=0.39909333141434805)
drift_676 = DRIFT(name="drift_676", len=0.22622760000012931)
bi12_bh19 = MARKER(name="bi12_bh19")
drift_677 = DRIFT(name="drift_677", len=1.0141595258655798)
bi12_dh19 = ESBEND(name="bi12_dh19", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_678 = DRIFT(name="drift_678", len=1.3261871258664542)
bi12_tv20 = CORRECTOR(name="bi12_tv20", len=0.5, xkick=0.0, ykick=-0.0)
drift_679 = DRIFT(name="drift_679", len=0.2952500000001237)
bi12_qd20 = KQUAD(name="bi12_qd20", len=1.11, k1=-0.08364335555329494)
drift_680 = DRIFT(name="drift_680", len=0.13104999999995925)
bi12_sxd20 = KSEXT(name="bi12_sxd20", len=0.75, k2=-0.6290818079014571)
drift_681 = DRIFT(name="drift_681", len=0.22622760000012931)
bi12_bv20 = MARKER(name="bi12_bv20")
drift_682 = DRIFT(name="drift_682", len=1.0141595258655798)
bi12_dh20 = ESBEND(name="bi12_dh20", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_683 = DRIFT(name="drift_683", len=1.3261871258664542)
bi1_th21 = CORRECTOR(name="bi1_th21", len=0.5, xkick=-0.0, ykick=0.0)
drift_684 = DRIFT(name="drift_684", len=0.2952500000001237)
bi1_qf21 = KQUAD(name="bi1_qf21", len=1.11, k1=0.08048071257602774)
drift_685 = DRIFT(name="drift_685", len=0.13104999999995925)
bi1_sxf21 = KSEXT(name="bi1_sxf21", len=0.75, k2=0.39909333141434805)
drift_686 = DRIFT(name="drift_686", len=0.22622760000012931)
bi1_bh21 = MARKER(name="bi1_bh21")
drift_687 = DRIFT(name="drift_687", len=1.0141595258655798)
bi1_dh20 = ESBEND(name="bi1_dh20", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_688 = DRIFT(name="drift_688", len=1.3261871258664542)
bi1_tv20 = CORRECTOR(name="bi1_tv20", len=0.5, xkick=0.0, ykick=-0.0)
drift_689 = DRIFT(name="drift_689", len=0.2952500000001237)
bi1_qd20 = KQUAD(name="bi1_qd20", len=1.11, k1=-0.08364335555329494)
drift_690 = DRIFT(name="drift_690", len=0.13104999999995925)
bi1_sxd20 = KSEXT(name="bi1_sxd20", len=0.75, k2=-0.6290818079014571)
drift_691 = DRIFT(name="drift_691", len=0.22622760000012931)
bi1_bv20 = MARKER(name="bi1_bv20")
drift_692 = DRIFT(name="drift_692", len=1.0141595258655798)
bi1_dh19 = ESBEND(name="bi1_dh19", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_693 = DRIFT(name="drift_693", len=1.3261871258664542)
bi1_th19 = CORRECTOR(name="bi1_th19", len=0.5, xkick=-0.0, ykick=0.0)
drift_694 = DRIFT(name="drift_694", len=0.2952500000001237)
bi1_qf19 = KQUAD(name="bi1_qf19", len=1.11, k1=0.08048071257602774)
drift_695 = DRIFT(name="drift_695", len=0.13104999999995925)
bi1_sxf19 = KSEXT(name="bi1_sxf19", len=0.75, k2=0.39909333141434805)
drift_696 = DRIFT(name="drift_696", len=0.22622760000012931)
bi1_bh19 = MARKER(name="bi1_bh19")
drift_697 = DRIFT(name="drift_697", len=1.0141595258655798)
bi1_dh18 = ESBEND(name="bi1_dh18", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_698 = DRIFT(name="drift_698", len=1.3261871258664542)
bi1_tv18 = CORRECTOR(name="bi1_tv18", len=0.5, xkick=0.0, ykick=-0.0)
drift_699 = DRIFT(name="drift_699", len=0.2952500000001237)
bi1_qd18 = KQUAD(name="bi1_qd18", len=1.11, k1=-0.08364335555329494)
drift_700 = DRIFT(name="drift_700", len=0.13104999999995925)
bi1_sxd18 = KSEXT(name="bi1_sxd18", len=0.75, k2=-0.6290818079014571)
drift_701 = DRIFT(name="drift_701", len=0.22622760000012931)
bi1_bv18 = MARKER(name="bi1_bv18")
drift_702 = DRIFT(name="drift_702", len=1.0141595258655798)
bi1_dh17 = ESBEND(name="bi1_dh17", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_703 = DRIFT(name="drift_703", len=1.3261871258664542)
bi1_th17 = CORRECTOR(name="bi1_th17", len=0.5, xkick=-0.0, ykick=0.0)
drift_704 = DRIFT(name="drift_704", len=0.2952500000001237)
bi1_qf17 = KQUAD(name="bi1_qf17", len=1.11, k1=0.08048071257602774)
drift_705 = DRIFT(name="drift_705", len=0.13104999999995925)
bi1_sxf17 = KSEXT(name="bi1_sxf17", len=0.75, k2=0.39909333141434805)
drift_706 = DRIFT(name="drift_706", len=0.22622760000012931)
bi1_bh17 = MARKER(name="bi1_bh17")
drift_707 = DRIFT(name="drift_707", len=1.0141595258655798)
bi1_dh16 = ESBEND(name="bi1_dh16", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_708 = DRIFT(name="drift_708", len=1.3261871258664542)
bi1_tv16 = CORRECTOR(name="bi1_tv16", len=0.5, xkick=0.0, ykick=-0.0)
drift_709 = DRIFT(name="drift_709", len=0.2952500000001237)
bi1_qd16 = KQUAD(name="bi1_qd16", len=1.11, k1=-0.08364335555329494)
drift_710 = DRIFT(name="drift_710", len=0.13104999999995925)
bi1_sxd16 = KSEXT(name="bi1_sxd16", len=0.75, k2=-0.6290818079014571)
drift_711 = DRIFT(name="drift_711", len=0.22622760000012931)
bi1_bv16 = MARKER(name="bi1_bv16")
drift_712 = DRIFT(name="drift_712", len=1.0141595258655798)
bi1_dh15 = ESBEND(name="bi1_dh15", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_713 = DRIFT(name="drift_713", len=1.3261871258664542)
bi1_th15 = CORRECTOR(name="bi1_th15", len=0.5, xkick=-0.0, ykick=0.0)
drift_714 = DRIFT(name="drift_714", len=0.2952500000001237)
bi1_qf15 = KQUAD(name="bi1_qf15", len=1.11, k1=0.08048071257602774)
drift_715 = DRIFT(name="drift_715", len=0.13104999999995925)
bi1_sxf15 = KSEXT(name="bi1_sxf15", len=0.75, k2=0.39909333141434805)
drift_716 = DRIFT(name="drift_716", len=0.22622760000012931)
bi1_bh15 = MARKER(name="bi1_bh15")
drift_717 = DRIFT(name="drift_717", len=1.0141595258655798)
bi1_dh14 = ESBEND(name="bi1_dh14", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_718 = DRIFT(name="drift_718", len=1.3261871258664542)
bi1_tv14 = CORRECTOR(name="bi1_tv14", len=0.5, xkick=0.0, ykick=-0.0)
drift_719 = DRIFT(name="drift_719", len=0.2952500000001237)
bi1_qd14 = KQUAD(name="bi1_qd14", len=1.11, k1=-0.08364335555329494)
drift_720 = DRIFT(name="drift_720", len=0.13104999999995925)
bi1_sxd14 = KSEXT(name="bi1_sxd14", len=0.75, k2=-0.6290818079014571)
drift_721 = DRIFT(name="drift_721", len=0.22622760000012931)
bi1_bv14 = MARKER(name="bi1_bv14")
drift_722 = DRIFT(name="drift_722", len=1.0141595258655798)
bi1_dh13 = ESBEND(name="bi1_dh13", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_723 = DRIFT(name="drift_723", len=1.3261871258664542)
bi1_th13 = CORRECTOR(name="bi1_th13", len=0.5, xkick=-0.0, ykick=0.0)
drift_724 = DRIFT(name="drift_724", len=0.2952500000001237)
bi1_qf13 = KQUAD(name="bi1_qf13", len=1.11, k1=0.08048071257602774)
drift_725 = DRIFT(name="drift_725", len=0.13104999999995925)
bi1_sxf13 = KSEXT(name="bi1_sxf13", len=0.75, k2=0.39909333141434805)
drift_726 = DRIFT(name="drift_726", len=0.22622760000012931)
bi1_bh13 = MARKER(name="bi1_bh13")
drift_727 = DRIFT(name="drift_727", len=1.0141595258655798)
bi1_dh12 = ESBEND(name="bi1_dh12", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_728 = DRIFT(name="drift_728", len=1.3261871258664542)
bi1_tv12 = CORRECTOR(name="bi1_tv12", len=0.5, xkick=0.0, ykick=-0.0)
drift_729 = DRIFT(name="drift_729", len=0.2952500000001237)
bi1_qd12 = KQUAD(name="bi1_qd12", len=1.11, k1=-0.08364335555329494)
drift_730 = DRIFT(name="drift_730", len=0.13104999999995925)
bi1_sxd12 = KSEXT(name="bi1_sxd12", len=0.75, k2=-0.6290818079014571)
drift_731 = DRIFT(name="drift_731", len=0.22622760000012931)
bi1_bv12 = MARKER(name="bi1_bv12")
drift_732 = DRIFT(name="drift_732", len=1.0141595258655798)
bi1_dh11 = ESBEND(name="bi1_dh11", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_733 = DRIFT(name="drift_733", len=1.3261871258664542)
bi1_th11 = CORRECTOR(name="bi1_th11", len=0.5, xkick=-0.0, ykick=0.0)
drift_734 = DRIFT(name="drift_734", len=0.2952500000001237)
bi1_qf11 = KQUAD(name="bi1_qf11", len=1.11, k1=0.08048071257602774)
drift_735 = DRIFT(name="drift_735", len=0.13104999999995925)
bi1_sxf11 = KSEXT(name="bi1_sxf11", len=0.75, k2=0.39909333141434805)
drift_736 = DRIFT(name="drift_736", len=0.22622760000012931)
bi1_bh11 = MARKER(name="bi1_bh11")
drift_737 = DRIFT(name="drift_737", len=1.0141595258655798)
bi1_dh10 = ESBEND(name="bi1_dh10", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_738 = DRIFT(name="drift_738", len=1.5761871258664542)
bi1_tv10 = CORRECTOR(name="bi1_tv10", len=0.0, xkick=0.0, ykick=-0.0)
bi1_oct10 = thinMULTIPOLE(name="bi1_oct10", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
bi1_dec10 = thinMULTIPOLE(name="bi1_dec10", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
bi1_qs10 = thinMULTIPOLE(name="bi1_qs10", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_739 = DRIFT(name="drift_739", len=0.5452500000001237)
bi1_qd10 = KQUAD(name="bi1_qd10", len=1.11, k1=-0.08364335555329494)
drift_740 = DRIFT(name="drift_740", len=0.13104999999995925)
bi1_sxd10 = KSEXT(name="bi1_sxd10", len=0.75, k2=-0.6290818079014571)
drift_741 = DRIFT(name="drift_741", len=0.22622760000012931)
bi1_bv10 = MARKER(name="bi1_bv10")
drift_742 = DRIFT(name="drift_742", len=1.016969423264527)
bi1_dh9 = ESBEND(name="bi1_dh9", len=2.94942686017, angle=0.012160550077129212, e1=0.0, e2=0.0)
drift_743 = DRIFT(name="drift_743", len=7.514623333294367)
bi1_bh9 = MARKER(name="bi1_bh9")
drift_744 = DRIFT(name="drift_744", len=1.1072776000000886)
bi1_qf9 = KQUAD(name="bi1_qf9", len=1.11, k1=0.04618735786934114)
drift_745 = DRIFT(name="drift_745", len=0.5452500000001237)
bi1_th9 = CORRECTOR(name="bi1_th9", len=0.0, xkick=-0.0, ykick=0.0)
bi1_oct9 = thinMULTIPOLE(name="bi1_oct9", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
bi1_dec9 = thinMULTIPOLE(name="bi1_dec9", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
bi1_qs9 = thinMULTIPOLE(name="bi1_qs9", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_746 = DRIFT(name="drift_746", len=1.5761871258664542)
bi1_dh8 = ESBEND(name="bi1_dh8", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_747 = DRIFT(name="drift_747", len=1.5761871258664542)
bi1_tv8 = CORRECTOR(name="bi1_tv8", len=0.0, xkick=0.0, ykick=-0.0)
bi1_oct8 = thinMULTIPOLE(name="bi1_oct8", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
bi1_dec8 = thinMULTIPOLE(name="bi1_dec8", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
bi1_qs8 = thinMULTIPOLE(name="bi1_qs8", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_748 = DRIFT(name="drift_748", len=0.5452500000001237)
bi1_qd8 = KQUAD(name="bi1_qd8", len=1.11, k1=0.014042704606404404)
drift_749 = DRIFT(name="drift_749", len=0.29610560000037367)
bi1_b8 = MARKER(name="bi1_b8")
drift_750 = DRIFT(name="drift_750", len=1.3184179478494116)
bi1_snk_hlx4 = DRIFT(name="bi1_snk_hlx4", len=2.4)
drift_751 = DRIFT(name="drift_751", len=0.21200000000044383)
bi1_snk_hlx3 = DRIFT(name="bi1_snk_hlx3", len=2.4)
drift_752 = DRIFT(name="drift_752", len=0.22400000000016007)
bi1_bsnk = MARKER(name="bi1_bsnk")
drift_753 = DRIFT(name="drift_753", len=0.22400000000016007)
bi1_snk_hlx2 = DRIFT(name="bi1_snk_hlx2", len=2.4)
drift_754 = DRIFT(name="drift_754", len=0.21200000000044383)
bi1_snk_hlx1 = DRIFT(name="bi1_snk_hlx1", len=2.4)
drift_755 = DRIFT(name="drift_755", len=1.8803175478496996)
yo1_th6 = CORRECTOR(name="yo1_th6", len=0.0, xkick=-0.0, ykick=0.0)
yo1_oct6 = thinMULTIPOLE(name="yo1_oct6", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yo1_dec6 = thinMULTIPOLE(name="yo1_dec6", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yo1_qgt6 = thinMULTIPOLE(name="yo1_qgt6", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, -0.0, 0.0, 0.0])
drift_756 = DRIFT(name="drift_756", len=0.5452500000001237)
yo1_qf6 = KQUAD(name="yo1_qf6", len=1.11, k1=-0.09203165083457603)
drift_757 = DRIFT(name="drift_757", len=0.13104999999995925)
yo1_tq6 = KQUAD(name="yo1_tq6", len=0.75, k1=0.0)
drift_758 = DRIFT(name="drift_758", len=0.22622760000012931)
yo1_bh6 = MARKER(name="yo1_bh6")
drift_759 = DRIFT(name="drift_759", len=1.3563004000006913)
bi1_tv6 = CORRECTOR(name="bi1_tv6", len=0.0, xkick=0.0, ykick=-0.0)
bi1_oct6 = thinMULTIPOLE(name="bi1_oct6", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
bi1_dec6 = thinMULTIPOLE(name="bi1_dec6", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
bi1_qs6 = thinMULTIPOLE(name="bi1_qs6", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_760 = DRIFT(name="drift_760", len=0.5452500000001237)
bi1_qd6 = KQUAD(name="bi1_qd6", len=1.11, k1=0.07582674963387082)
drift_761 = DRIFT(name="drift_761", len=0.13104999999995925)
bi1_tq6 = KQUAD(name="bi1_tq6", len=0.75, k1=-0.0)
drift_762 = DRIFT(name="drift_762", len=0.22622760000012931)
bi1_bv6 = MARKER(name="bi1_bv6")
drift_763 = DRIFT(name="drift_763", len=1.0395430279659195)
bi1_dh6 = ESBEND(name="bi1_dh6", len=2.94942686017, angle=0.012160550077129212, e1=0.0, e2=0.0)
drift_764 = DRIFT(name="drift_764", len=0.517281632902268)
cool_kicker_out = DRIFT(name="cool_kicker_out", len=22.736407378961996)
cool_kicker_in = DRIFT(name="cool_kicker_in", len=15.845638500000002)
drift_765 = DRIFT(name="drift_765", len=11.609049000000596)
yo1_th4 = CORRECTOR(name="yo1_th4", len=0.0, xkick=-0.0, ykick=0.0)
yo1_oct4 = thinMULTIPOLE(name="yo1_oct4", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yo1_dec4 = thinMULTIPOLE(name="yo1_dec4", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yo1_qs4 = thinMULTIPOLE(name="yo1_qs4", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_766 = DRIFT(name="drift_766", len=0.5442755000003672)
yo1_qf4 = KQUAD(name="yo1_qf4", len=1.811949, k1=0.04422426805488868)
drift_767 = DRIFT(name="drift_767", len=0.13007550000020274)
yo1_tq4 = KQUAD(name="yo1_tq4", len=0.75, k1=0.0)
drift_768 = DRIFT(name="drift_768", len=0.25163960000008956)
yo1_b4 = MARKER(name="yo1_b4")
drift_769 = DRIFT(name="drift_769", len=1.8819324000005508)
bi1_th5 = CORRECTOR(name="bi1_th5", len=0.0, xkick=-0.0, ykick=0.0)
bi1_oct5 = thinMULTIPOLE(name="bi1_oct5", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
bi1_dec5 = thinMULTIPOLE(name="bi1_dec5", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
bi1_qgt5 = thinMULTIPOLE(name="bi1_qgt5", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, -0.0, 0.0, 0.0])
drift_770 = DRIFT(name="drift_770", len=0.5452500000001237)
bi1_qf5 = KQUAD(name="bi1_qf5", len=1.11, k1=-0.05157130904514161)
drift_771 = DRIFT(name="drift_771", len=0.13104999999995925)
bi1_tq5 = KQUAD(name="bi1_tq5", len=0.75, k1=-0.0)
drift_772 = DRIFT(name="drift_772", len=0.22622760000012931)
bi1_bh5 = MARKER(name="bi1_bh5")
drift_773 = DRIFT(name="drift_773", len=1.7858110860697707)
d5 = ESBEND(name="d5", len=7.807224627859998, angle=0.032189354254750294, e1=0.0, e2=0.0)
drift_774 = DRIFT(name="drift_774", len=1.8760670860701794)
yo1_bv9 = MARKER(name="yo1_bv9")
drift_775 = DRIFT(name="drift_775", len=1.1072776000000886)
yo1_qd9 = KQUAD(name="yo1_qd9", len=1.11, k1=-0.03159359878005887)
drift_776 = DRIFT(name="drift_776", len=0.5452500000001237)
yo1_tv9 = CORRECTOR(name="yo1_tv9", len=0.0, xkick=0.0, ykick=-0.0)
yo1_oct9 = thinMULTIPOLE(name="yo1_oct9", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yo1_dec9 = thinMULTIPOLE(name="yo1_dec9", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yo1_qs9 = thinMULTIPOLE(name="yo1_qs9", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_777 = DRIFT(name="drift_777", len=12.592901873552364)
yo1_tv5 = CORRECTOR(name="yo1_tv5", len=0.0, xkick=0.0, ykick=-0.0)
yo1_oct5 = thinMULTIPOLE(name="yo1_oct5", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yo1_dec5 = thinMULTIPOLE(name="yo1_dec5", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yo1_qs5 = thinMULTIPOLE(name="yo1_qs5", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_778 = DRIFT(name="drift_778", len=0.5452500000001237)
yo1_qd5 = KQUAD(name="yo1_qd5", len=1.11, k1=0.038587509473826616)
drift_779 = DRIFT(name="drift_779", len=0.13104999999995925)
yo1_tq5 = KQUAD(name="yo1_tq5", len=0.75, k1=0.0)
drift_780 = DRIFT(name="drift_780", len=0.22622760000012931)
yo1_bv5 = MARKER(name="yo1_bv5")
drift_781 = DRIFT(name="drift_781", len=12.370147075600471)
bi1_b7 = MARKER(name="bi1_b7")
drift_782 = DRIFT(name="drift_782", len=0.2962336000000505)
bi1_qf7 = KQUAD(name="bi1_qf7", len=0.929744, k1=0.004937753087673467)
drift_783 = DRIFT(name="drift_783", len=0.5453779999998005)
bi1_th7 = CORRECTOR(name="bi1_th7", len=0.0, xkick=-0.0, ykick=0.0)
bi1_oct7 = thinMULTIPOLE(name="bi1_oct7", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
bi1_dec7 = thinMULTIPOLE(name="bi1_dec7", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
bi1_qgt7 = thinMULTIPOLE(name="bi1_qgt7", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, -0.0, 0.0, 0.0])
drift_784 = DRIFT(name="drift_784", len=12.121002675600721)
bo2_bv5 = MARKER(name="bo2_bv5")
drift_785 = DRIFT(name="drift_785", len=0.22622760000012931)
bo2_tq5 = KQUAD(name="bo2_tq5", len=0.75, k1=-0.0)
drift_786 = DRIFT(name="drift_786", len=0.13104999999995925)
bo2_qd5 = KQUAD(name="bo2_qd5", len=1.11, k1=0.038587509473826616)
drift_787 = DRIFT(name="drift_787", len=0.5452500000001237)
bo2_tv5 = CORRECTOR(name="bo2_tv5", len=0.0, xkick=0.0, ykick=-0.0)
bo2_oct5 = thinMULTIPOLE(name="bo2_oct5", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
bo2_dec5 = thinMULTIPOLE(name="bo2_dec5", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
bo2_qs5 = thinMULTIPOLE(name="bo2_qs5", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_788 = DRIFT(name="drift_788", len=12.592901873552364)
bo2_tv9 = CORRECTOR(name="bo2_tv9", len=0.0, xkick=0.0, ykick=-0.0)
bo2_oct9 = thinMULTIPOLE(name="bo2_oct9", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
bo2_dec9 = thinMULTIPOLE(name="bo2_dec9", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
bo2_qs9 = thinMULTIPOLE(name="bo2_qs9", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_789 = DRIFT(name="drift_789", len=0.5452500000001237)
bo2_qd9 = KQUAD(name="bo2_qd9", len=1.11, k1=-0.03159359878005887)
drift_790 = DRIFT(name="drift_790", len=1.1072776000000886)
bo2_bv9 = MARKER(name="bo2_bv9")
drift_791 = DRIFT(name="drift_791", len=1.8760670860701794)
drift_792 = DRIFT(name="drift_792", len=1.7858110860697707)
yi2_bh5 = MARKER(name="yi2_bh5")
drift_793 = DRIFT(name="drift_793", len=0.22622760000012931)
yi2_tq5 = KQUAD(name="yi2_tq5", len=0.75, k1=0.0)
drift_794 = DRIFT(name="drift_794", len=0.13104999999995925)
yi2_qf5 = KQUAD(name="yi2_qf5", len=1.11, k1=-0.05157130904514161)
drift_795 = DRIFT(name="drift_795", len=0.5452500000001237)
yi2_th5 = CORRECTOR(name="yi2_th5", len=0.0, xkick=-0.0, ykick=0.0)
yi2_oct5 = thinMULTIPOLE(name="yi2_oct5", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yi2_dec5 = thinMULTIPOLE(name="yi2_dec5", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yi2_qgt5 = thinMULTIPOLE(name="yi2_qgt5", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, -0.0, 0.0, 0.0])
drift_796 = DRIFT(name="drift_796", len=1.8819324000005508)
bo2_b4 = MARKER(name="bo2_b4")
drift_797 = DRIFT(name="drift_797", len=0.25163960000008956)
bo2_tq4 = KQUAD(name="bo2_tq4", len=0.75, k1=-0.0)
drift_798 = DRIFT(name="drift_798", len=0.13007550000020274)
bo2_qf4 = KQUAD(name="bo2_qf4", len=1.811949, k1=0.04422426805488868)
drift_799 = DRIFT(name="drift_799", len=0.5442755000003672)
bo2_th4 = CORRECTOR(name="bo2_th4", len=0.0, xkick=-0.0, ykick=0.0)
bo2_oct4 = thinMULTIPOLE(name="bo2_oct4", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
bo2_dec4 = thinMULTIPOLE(name="bo2_dec4", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
bo2_qs4 = thinMULTIPOLE(name="bo2_qs4", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_800 = DRIFT(name="drift_800", len=11.609049000000596)
cool_modulator_in = DRIFT(name="cool_modulator_in", len=15.845638500000002)
cool_modulator_out = DRIFT(name="cool_modulator_out", len=22.736407378961996)
drift_801 = DRIFT(name="drift_801", len=0.517281632902268)
yi2_dh6 = ESBEND(name="yi2_dh6", len=2.94942686017, angle=0.012160550077129212, e1=0.0, e2=0.0)
drift_802 = DRIFT(name="drift_802", len=1.0395430279659195)
yi2_bv6 = MARKER(name="yi2_bv6")
drift_803 = DRIFT(name="drift_803", len=0.22622760000012931)
yi2_tq6 = KQUAD(name="yi2_tq6", len=0.75, k1=0.0)
drift_804 = DRIFT(name="drift_804", len=0.13104999999995925)
yi2_qd6 = KQUAD(name="yi2_qd6", len=1.11, k1=0.07582674963387082)
drift_805 = DRIFT(name="drift_805", len=0.5452500000001237)
yi2_tv6 = CORRECTOR(name="yi2_tv6", len=0.0, xkick=0.0, ykick=-0.0)
yi2_oct6 = thinMULTIPOLE(name="yi2_oct6", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yi2_dec6 = thinMULTIPOLE(name="yi2_dec6", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yi2_qs6 = thinMULTIPOLE(name="yi2_qs6", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_806 = DRIFT(name="drift_806", len=1.3563004000006913)
bo2_bh6 = MARKER(name="bo2_bh6")
drift_807 = DRIFT(name="drift_807", len=0.22622760000012931)
bo2_tq6 = KQUAD(name="bo2_tq6", len=0.75, k1=-0.0)
drift_808 = DRIFT(name="drift_808", len=0.13104999999995925)
bo2_qf6 = KQUAD(name="bo2_qf6", len=1.11, k1=-0.09203165083457603)
drift_809 = DRIFT(name="drift_809", len=0.5452500000001237)
bo2_th6 = CORRECTOR(name="bo2_th6", len=0.0, xkick=-0.0, ykick=0.0)
bo2_oct6 = thinMULTIPOLE(name="bo2_oct6", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
bo2_dec6 = thinMULTIPOLE(name="bo2_dec6", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
bo2_qgt6 = thinMULTIPOLE(name="bo2_qgt6", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, -0.0, 0.0, 0.0])
drift_810 = DRIFT(name="drift_810", len=13.670735495699773)
yi2_b8 = MARKER(name="yi2_b8")
drift_811 = DRIFT(name="drift_811", len=0.29610560000037367)
yi2_qd8 = KQUAD(name="yi2_qd8", len=1.11, k1=0.014042704606404404)
drift_812 = DRIFT(name="drift_812", len=0.5452500000001237)
yi2_tv8 = CORRECTOR(name="yi2_tv8", len=0.0, xkick=0.0, ykick=-0.0)
yi2_oct8 = thinMULTIPOLE(name="yi2_oct8", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yi2_dec8 = thinMULTIPOLE(name="yi2_dec8", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yi2_qs8 = thinMULTIPOLE(name="yi2_qs8", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_813 = DRIFT(name="drift_813", len=1.5761871258664542)
yi2_dh8 = ESBEND(name="yi2_dh8", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_814 = DRIFT(name="drift_814", len=1.5761871258664542)
yi2_th9 = CORRECTOR(name="yi2_th9", len=0.0, xkick=-0.0, ykick=0.0)
yi2_oct9 = thinMULTIPOLE(name="yi2_oct9", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yi2_dec9 = thinMULTIPOLE(name="yi2_dec9", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yi2_qs9 = thinMULTIPOLE(name="yi2_qs9", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_815 = DRIFT(name="drift_815", len=0.5452500000001237)
yi2_qf9 = KQUAD(name="yi2_qf9", len=1.11, k1=0.04618735786934114)
drift_816 = DRIFT(name="drift_816", len=1.1072776000000886)
yi2_bh9 = MARKER(name="yi2_bh9")
drift_817 = DRIFT(name="drift_817", len=7.514623333294367)
yi2_dh9 = ESBEND(name="yi2_dh9", len=2.94942686017, angle=0.012160550077129212, e1=0.0, e2=0.0)
drift_818 = DRIFT(name="drift_818", len=1.016969423264527)
yi2_bv10 = MARKER(name="yi2_bv10")
drift_819 = DRIFT(name="drift_819", len=0.22622760000012931)
yi2_sxd10 = KSEXT(name="yi2_sxd10", len=0.75, k2=-0.6290818079014571)
drift_820 = DRIFT(name="drift_820", len=0.13104999999995925)
yi2_qd10 = KQUAD(name="yi2_qd10", len=1.11, k1=-0.08364335555329494)
drift_821 = DRIFT(name="drift_821", len=0.5452500000001237)
yi2_tv10 = CORRECTOR(name="yi2_tv10", len=0.0, xkick=0.0, ykick=-0.0)
yi2_oct10 = thinMULTIPOLE(name="yi2_oct10", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yi2_dec10 = thinMULTIPOLE(name="yi2_dec10", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yi2_qs10 = thinMULTIPOLE(name="yi2_qs10", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_822 = DRIFT(name="drift_822", len=1.5761871258664542)
yi2_dh10 = ESBEND(name="yi2_dh10", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_823 = DRIFT(name="drift_823", len=1.0141595258655798)
yi2_bh11 = MARKER(name="yi2_bh11")
drift_824 = DRIFT(name="drift_824", len=0.22622760000012931)
yi2_sxf11 = KSEXT(name="yi2_sxf11", len=0.75, k2=0.39909333141434805)
drift_825 = DRIFT(name="drift_825", len=0.13104999999995925)
yi2_qf11 = KQUAD(name="yi2_qf11", len=1.11, k1=0.08048071257602774)
drift_826 = DRIFT(name="drift_826", len=0.2952500000001237)
yi2_th11 = CORRECTOR(name="yi2_th11", len=0.5, xkick=0.0, ykick=0.0)
drift_827 = DRIFT(name="drift_827", len=1.3261871258664542)
yi2_dh11 = ESBEND(name="yi2_dh11", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_828 = DRIFT(name="drift_828", len=1.0141595258655798)
yi2_bv12 = MARKER(name="yi2_bv12")
drift_829 = DRIFT(name="drift_829", len=0.22622760000012931)
yi2_sxd12 = KSEXT(name="yi2_sxd12", len=0.75, k2=-0.6290818079014571)
drift_830 = DRIFT(name="drift_830", len=0.13104999999995925)
yi2_qd12 = KQUAD(name="yi2_qd12", len=1.11, k1=-0.08364335555329494)
drift_831 = DRIFT(name="drift_831", len=0.2952500000001237)
yi2_tv12 = CORRECTOR(name="yi2_tv12", len=0.5, xkick=0.0, ykick=0.0)
drift_832 = DRIFT(name="drift_832", len=1.3261871258664542)
yi2_dh12 = ESBEND(name="yi2_dh12", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_833 = DRIFT(name="drift_833", len=1.0141595258655798)
yi2_bh13 = MARKER(name="yi2_bh13")
drift_834 = DRIFT(name="drift_834", len=0.22622760000012931)
yi2_sxf13 = KSEXT(name="yi2_sxf13", len=0.75, k2=0.39909333141434805)
drift_835 = DRIFT(name="drift_835", len=0.13104999999995925)
yi2_qf13 = KQUAD(name="yi2_qf13", len=1.11, k1=0.08048071257602774)
drift_836 = DRIFT(name="drift_836", len=0.2952500000001237)
yi2_th13 = CORRECTOR(name="yi2_th13", len=0.5, xkick=0.0, ykick=0.0)
drift_837 = DRIFT(name="drift_837", len=1.3261871258664542)
yi2_dh13 = ESBEND(name="yi2_dh13", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_838 = DRIFT(name="drift_838", len=1.0141595258655798)
yi2_bv14 = MARKER(name="yi2_bv14")
drift_839 = DRIFT(name="drift_839", len=0.22622760000012931)
yi2_sxd14 = KSEXT(name="yi2_sxd14", len=0.75, k2=-0.6290818079014571)
drift_840 = DRIFT(name="drift_840", len=0.13104999999995925)
yi2_qd14 = KQUAD(name="yi2_qd14", len=1.11, k1=-0.08364335555329494)
drift_841 = DRIFT(name="drift_841", len=0.2952500000001237)
yi2_tv14 = CORRECTOR(name="yi2_tv14", len=0.5, xkick=0.0, ykick=0.0)
drift_842 = DRIFT(name="drift_842", len=1.3261871258664542)
yi2_dh14 = ESBEND(name="yi2_dh14", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_843 = DRIFT(name="drift_843", len=1.0141595258655798)
yi2_bh15 = MARKER(name="yi2_bh15")
drift_844 = DRIFT(name="drift_844", len=0.22622760000012931)
yi2_sxf15 = KSEXT(name="yi2_sxf15", len=0.75, k2=0.39909333141434805)
drift_845 = DRIFT(name="drift_845", len=0.13104999999995925)
yi2_qf15 = KQUAD(name="yi2_qf15", len=1.11, k1=0.08048071257602774)
drift_846 = DRIFT(name="drift_846", len=0.2952500000001237)
yi2_th15 = CORRECTOR(name="yi2_th15", len=0.5, xkick=0.0, ykick=0.0)
drift_847 = DRIFT(name="drift_847", len=1.3261871258664542)
yi2_dh15 = ESBEND(name="yi2_dh15", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_848 = DRIFT(name="drift_848", len=1.0141595258655798)
yi2_bv16 = MARKER(name="yi2_bv16")
drift_849 = DRIFT(name="drift_849", len=0.22622760000012931)
yi2_sxd16 = KSEXT(name="yi2_sxd16", len=0.75, k2=-0.6290818079014571)
drift_850 = DRIFT(name="drift_850", len=0.13104999999995925)
yi2_qd16 = KQUAD(name="yi2_qd16", len=1.11, k1=-0.08364335555329494)
drift_851 = DRIFT(name="drift_851", len=0.2952500000001237)
yi2_tv16 = CORRECTOR(name="yi2_tv16", len=0.5, xkick=0.0, ykick=0.0)
drift_852 = DRIFT(name="drift_852", len=1.3261871258664542)
yi2_dh16 = ESBEND(name="yi2_dh16", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_853 = DRIFT(name="drift_853", len=1.0141595258655798)
yi2_bh17 = MARKER(name="yi2_bh17")
drift_854 = DRIFT(name="drift_854", len=0.22622760000012931)
yi2_sxf17 = KSEXT(name="yi2_sxf17", len=0.75, k2=0.39909333141434805)
drift_855 = DRIFT(name="drift_855", len=0.13104999999995925)
yi2_qf17 = KQUAD(name="yi2_qf17", len=1.11, k1=0.08048071257602774)
drift_856 = DRIFT(name="drift_856", len=0.2952500000001237)
yi2_th17 = CORRECTOR(name="yi2_th17", len=0.5, xkick=0.0, ykick=0.0)
drift_857 = DRIFT(name="drift_857", len=1.3261871258664542)
yi2_dh17 = ESBEND(name="yi2_dh17", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_858 = DRIFT(name="drift_858", len=1.0141595258655798)
yi2_bv18 = MARKER(name="yi2_bv18")
drift_859 = DRIFT(name="drift_859", len=0.22622760000012931)
yi2_sxd18 = KSEXT(name="yi2_sxd18", len=0.75, k2=-0.6290818079014571)
drift_860 = DRIFT(name="drift_860", len=0.13104999999995925)
yi2_qd18 = KQUAD(name="yi2_qd18", len=1.11, k1=-0.08364335555329494)
drift_861 = DRIFT(name="drift_861", len=0.2952500000001237)
yi2_tv18 = CORRECTOR(name="yi2_tv18", len=0.5, xkick=0.0, ykick=0.0)
drift_862 = DRIFT(name="drift_862", len=1.3261871258664542)
yi2_dh18 = ESBEND(name="yi2_dh18", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_863 = DRIFT(name="drift_863", len=1.0141595258655798)
yi2_bh19 = MARKER(name="yi2_bh19")
drift_864 = DRIFT(name="drift_864", len=0.22622760000012931)
yi2_sxf19 = KSEXT(name="yi2_sxf19", len=0.75, k2=0.39909333141434805)
drift_865 = DRIFT(name="drift_865", len=0.13104999999995925)
yi2_qf19 = KQUAD(name="yi2_qf19", len=1.11, k1=0.08048071257602774)
drift_866 = DRIFT(name="drift_866", len=0.2952500000001237)
yi2_th19 = CORRECTOR(name="yi2_th19", len=0.5, xkick=0.0, ykick=0.0)
drift_867 = DRIFT(name="drift_867", len=1.3261871258664542)
yi2_dh19 = ESBEND(name="yi2_dh19", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_868 = DRIFT(name="drift_868", len=1.0141595258655798)
yi2_bv20 = MARKER(name="yi2_bv20")
drift_869 = DRIFT(name="drift_869", len=0.22622760000012931)
yi2_sxd20 = KSEXT(name="yi2_sxd20", len=0.75, k2=-0.6290818079014571)
drift_870 = DRIFT(name="drift_870", len=0.13104999999995925)
yi2_qd20 = KQUAD(name="yi2_qd20", len=1.11, k1=-0.08364335555329494)
drift_871 = DRIFT(name="drift_871", len=0.2952500000001237)
yi2_tv20 = CORRECTOR(name="yi2_tv20", len=0.5, xkick=0.0, ykick=0.0)
drift_872 = DRIFT(name="drift_872", len=1.3261871258664542)
yi2_dh20 = ESBEND(name="yi2_dh20", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_873 = DRIFT(name="drift_873", len=1.0141595258655798)
yi3_bh21 = MARKER(name="yi3_bh21")
drift_874 = DRIFT(name="drift_874", len=0.22622760000012931)
yi3_sxf21 = KSEXT(name="yi3_sxf21", len=0.75, k2=0.39909333141434805)
drift_875 = DRIFT(name="drift_875", len=0.13104999999995925)
yi3_qf21 = KQUAD(name="yi3_qf21", len=1.11, k1=0.08048071257602774)
drift_876 = DRIFT(name="drift_876", len=0.2952500000001237)
yi3_th21 = CORRECTOR(name="yi3_th21", len=0.5, xkick=0.0, ykick=0.0)
drift_877 = DRIFT(name="drift_877", len=1.3261871258664542)
yi3_dh20 = ESBEND(name="yi3_dh20", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_878 = DRIFT(name="drift_878", len=1.0141595258655798)
yi3_bv20 = MARKER(name="yi3_bv20")
drift_879 = DRIFT(name="drift_879", len=0.22622760000012931)
yi3_sxd20 = KSEXT(name="yi3_sxd20", len=0.75, k2=-0.6290818079014571)
drift_880 = DRIFT(name="drift_880", len=0.13104999999995925)
yi3_qd20 = KQUAD(name="yi3_qd20", len=1.11, k1=-0.08364335555329494)
drift_881 = DRIFT(name="drift_881", len=0.2952500000001237)
yi3_tv20 = CORRECTOR(name="yi3_tv20", len=0.5, xkick=0.0, ykick=0.0)
drift_882 = DRIFT(name="drift_882", len=1.3261871258664542)
yi3_dh19 = ESBEND(name="yi3_dh19", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_883 = DRIFT(name="drift_883", len=1.0141595258655798)
yi3_bh19 = MARKER(name="yi3_bh19")
drift_884 = DRIFT(name="drift_884", len=0.22622760000012931)
yi3_sxf19 = KSEXT(name="yi3_sxf19", len=0.75, k2=0.39909333141434805)
drift_885 = DRIFT(name="drift_885", len=0.13104999999995925)
yi3_qf19 = KQUAD(name="yi3_qf19", len=1.11, k1=0.08048071257602774)
drift_886 = DRIFT(name="drift_886", len=0.2952500000001237)
yi3_th19 = CORRECTOR(name="yi3_th19", len=0.5, xkick=0.0, ykick=0.0)
drift_887 = DRIFT(name="drift_887", len=1.3261871258664542)
yi3_dh18 = ESBEND(name="yi3_dh18", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_888 = DRIFT(name="drift_888", len=1.0141595258655798)
yi3_bv18 = MARKER(name="yi3_bv18")
drift_889 = DRIFT(name="drift_889", len=0.22622760000012931)
yi3_sxd18 = KSEXT(name="yi3_sxd18", len=0.75, k2=-0.6290818079014571)
drift_890 = DRIFT(name="drift_890", len=0.13104999999995925)
yi3_qd18 = KQUAD(name="yi3_qd18", len=1.11, k1=-0.08364335555329494)
drift_891 = DRIFT(name="drift_891", len=0.5452500000001237)
yi3_tv18 = CORRECTOR(name="yi3_tv18", len=0.0, xkick=0.0, ykick=-0.0)
yi3_oct18 = thinMULTIPOLE(name="yi3_oct18", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yi3_dec18 = thinMULTIPOLE(name="yi3_dec18", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yi3_qs18 = thinMULTIPOLE(name="yi3_qs18", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_892 = DRIFT(name="drift_892", len=1.5761871258664542)
yi3_dh17 = ESBEND(name="yi3_dh17", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_893 = DRIFT(name="drift_893", len=1.0141595258655798)
yi3_bh17 = MARKER(name="yi3_bh17")
drift_894 = DRIFT(name="drift_894", len=0.22622760000012931)
yi3_sxf17 = KSEXT(name="yi3_sxf17", len=0.75, k2=0.39909333141434805)
drift_895 = DRIFT(name="drift_895", len=0.13104999999995925)
yi3_qf17 = KQUAD(name="yi3_qf17", len=1.11, k1=0.08048071257602774)
drift_896 = DRIFT(name="drift_896", len=0.5452500000001237)
yi3_th17 = CORRECTOR(name="yi3_th17", len=0.0, xkick=-0.0, ykick=0.0)
yi3_oct17 = thinMULTIPOLE(name="yi3_oct17", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yi3_dec17 = thinMULTIPOLE(name="yi3_dec17", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yi3_qgt17 = thinMULTIPOLE(name="yi3_qgt17", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, -0.0, 0.0, 0.0])
drift_897 = DRIFT(name="drift_897", len=1.5761871258664542)
yi3_dh16 = ESBEND(name="yi3_dh16", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_898 = DRIFT(name="drift_898", len=1.0141595258655798)
yi3_bv16 = MARKER(name="yi3_bv16")
drift_899 = DRIFT(name="drift_899", len=0.22622760000012931)
yi3_sxd16 = KSEXT(name="yi3_sxd16", len=0.75, k2=-0.6290818079014571)
drift_900 = DRIFT(name="drift_900", len=0.13104999999995925)
yi3_qd16 = KQUAD(name="yi3_qd16", len=1.11, k1=-0.08364335555329494)
drift_901 = DRIFT(name="drift_901", len=0.5452500000001237)
yi3_tv16 = CORRECTOR(name="yi3_tv16", len=0.0, xkick=0.0, ykick=-0.0)
yi3_oct16 = thinMULTIPOLE(name="yi3_oct16", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yi3_dec16 = thinMULTIPOLE(name="yi3_dec16", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yi3_qs16 = thinMULTIPOLE(name="yi3_qs16", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_902 = DRIFT(name="drift_902", len=1.5761871258664542)
yi3_dh15 = ESBEND(name="yi3_dh15", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_903 = DRIFT(name="drift_903", len=1.0141595258655798)
yi3_bh15 = MARKER(name="yi3_bh15")
drift_904 = DRIFT(name="drift_904", len=0.22622760000012931)
yi3_sxf15 = KSEXT(name="yi3_sxf15", len=0.75, k2=0.39909333141434805)
drift_905 = DRIFT(name="drift_905", len=0.13104999999995925)
yi3_qf15 = KQUAD(name="yi3_qf15", len=1.11, k1=0.08048071257602774)
drift_906 = DRIFT(name="drift_906", len=0.5452500000001237)
yi3_th15 = CORRECTOR(name="yi3_th15", len=0.0, xkick=-0.0, ykick=0.0)
yi3_oct15 = thinMULTIPOLE(name="yi3_oct15", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yi3_dec15 = thinMULTIPOLE(name="yi3_dec15", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yi3_qgt15 = thinMULTIPOLE(name="yi3_qgt15", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, -0.0, 0.0, 0.0])
drift_907 = DRIFT(name="drift_907", len=1.5761871258664542)
yi3_dh14 = ESBEND(name="yi3_dh14", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_908 = DRIFT(name="drift_908", len=1.0141595258655798)
yi3_bv14 = MARKER(name="yi3_bv14")
drift_909 = DRIFT(name="drift_909", len=0.22622760000012931)
yi3_sxd14 = KSEXT(name="yi3_sxd14", len=0.75, k2=-0.6290818079014571)
drift_910 = DRIFT(name="drift_910", len=0.13104999999995925)
yi3_qd14 = KQUAD(name="yi3_qd14", len=1.11, k1=-0.08364335555329494)
drift_911 = DRIFT(name="drift_911", len=0.5452500000001237)
yi3_tv14 = CORRECTOR(name="yi3_tv14", len=0.0, xkick=0.0, ykick=-0.0)
yi3_oct14 = thinMULTIPOLE(name="yi3_oct14", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yi3_dec14 = thinMULTIPOLE(name="yi3_dec14", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yi3_qs14 = thinMULTIPOLE(name="yi3_qs14", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_912 = DRIFT(name="drift_912", len=1.5761871258664542)
yi3_dh13 = ESBEND(name="yi3_dh13", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_913 = DRIFT(name="drift_913", len=1.0141595258655798)
yi3_bh13 = MARKER(name="yi3_bh13")
drift_914 = DRIFT(name="drift_914", len=0.22622760000012931)
yi3_sxf13 = KSEXT(name="yi3_sxf13", len=0.75, k2=0.39909333141434805)
drift_915 = DRIFT(name="drift_915", len=0.13104999999995925)
yi3_qf13 = KQUAD(name="yi3_qf13", len=1.11, k1=0.08048071257602774)
drift_916 = DRIFT(name="drift_916", len=0.5452500000001237)
yi3_th13 = CORRECTOR(name="yi3_th13", len=0.0, xkick=-0.0, ykick=0.0)
yi3_oct13 = thinMULTIPOLE(name="yi3_oct13", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yi3_dec13 = thinMULTIPOLE(name="yi3_dec13", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yi3_qgt13 = thinMULTIPOLE(name="yi3_qgt13", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, -0.0, 0.0, 0.0])
drift_917 = DRIFT(name="drift_917", len=1.5761871258664542)
yi3_dh12 = ESBEND(name="yi3_dh12", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_918 = DRIFT(name="drift_918", len=1.0141595258655798)
yi3_bv12 = MARKER(name="yi3_bv12")
drift_919 = DRIFT(name="drift_919", len=0.22622760000012931)
yi3_sxd12 = KSEXT(name="yi3_sxd12", len=0.75, k2=-0.6290818079014571)
drift_920 = DRIFT(name="drift_920", len=0.13104999999995925)
yi3_qd12 = KQUAD(name="yi3_qd12", len=1.11, k1=-0.08364335555329494)
drift_921 = DRIFT(name="drift_921", len=0.5452500000001237)
yi3_tv12 = CORRECTOR(name="yi3_tv12", len=0.0, xkick=0.0, ykick=-0.0)
yi3_oct12 = thinMULTIPOLE(name="yi3_oct12", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yi3_dec12 = thinMULTIPOLE(name="yi3_dec12", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yi3_qs12 = thinMULTIPOLE(name="yi3_qs12", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_922 = DRIFT(name="drift_922", len=1.5761871258664542)
yi3_dh11 = ESBEND(name="yi3_dh11", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_923 = DRIFT(name="drift_923", len=1.0141595258655798)
yi3_bh11 = MARKER(name="yi3_bh11")
drift_924 = DRIFT(name="drift_924", len=0.22622760000012931)
yi3_sxf11 = KSEXT(name="yi3_sxf11", len=0.75, k2=0.39909333141434805)
drift_925 = DRIFT(name="drift_925", len=0.13104999999995925)
yi3_qf11 = KQUAD(name="yi3_qf11", len=1.11, k1=0.08048071257602774)
drift_926 = DRIFT(name="drift_926", len=0.5452500000001237)
yi3_th11 = CORRECTOR(name="yi3_th11", len=0.0, xkick=-0.0, ykick=0.0)
yi3_oct11 = thinMULTIPOLE(name="yi3_oct11", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yi3_dec11 = thinMULTIPOLE(name="yi3_dec11", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yi3_qgt11 = thinMULTIPOLE(name="yi3_qgt11", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, -0.0, 0.0, 0.0])
drift_927 = DRIFT(name="drift_927", len=1.5761871258664542)
yi3_dh10 = ESBEND(name="yi3_dh10", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_928 = DRIFT(name="drift_928", len=1.0141595258655798)
yi3_bv10 = MARKER(name="yi3_bv10")
drift_929 = DRIFT(name="drift_929", len=0.22622760000012931)
yi3_sxd10 = KSEXT(name="yi3_sxd10", len=0.75, k2=-0.6290818079014571)
drift_930 = DRIFT(name="drift_930", len=0.13104999999995925)
yi3_qd10 = KQUAD(name="yi3_qd10", len=1.11, k1=-0.08364335555329494)
drift_931 = DRIFT(name="drift_931", len=0.5452500000001237)
yi3_tv10 = CORRECTOR(name="yi3_tv10", len=0.0, xkick=0.0, ykick=-0.0)
yi3_oct10 = thinMULTIPOLE(name="yi3_oct10", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yi3_dec10 = thinMULTIPOLE(name="yi3_dec10", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yi3_qs10 = thinMULTIPOLE(name="yi3_qs10", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_932 = DRIFT(name="drift_932", len=1.5789970232644919)
yi3_dh9 = ESBEND(name="yi3_dh9", len=2.94942686017, angle=0.012160550077129212, e1=0.0, e2=0.0)
drift_933 = DRIFT(name="drift_933", len=7.514623333294367)
yi3_bh9 = MARKER(name="yi3_bh9")
drift_934 = DRIFT(name="drift_934", len=0.22622760000012931)
yi3_sxf9 = KSEXT(name="yi3_sxf9", len=0.75, k2=0.39909333141434805)
drift_935 = DRIFT(name="drift_935", len=0.13104999999995925)
yi3_qf9 = KQUAD(name="yi3_qf9", len=1.11, k1=0.07952250035938109)
drift_936 = DRIFT(name="drift_936", len=0.5452500000001237)
yi3_th9 = CORRECTOR(name="yi3_th9", len=0.0, xkick=-0.0, ykick=0.0)
yi3_oct9 = thinMULTIPOLE(name="yi3_oct9", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yi3_dec9 = thinMULTIPOLE(name="yi3_dec9", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yi3_qs9 = thinMULTIPOLE(name="yi3_qs9", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_937 = DRIFT(name="drift_937", len=1.5761871258664542)
yi3_dh8 = ESBEND(name="yi3_dh8", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_938 = DRIFT(name="drift_938", len=1.5761871258664542)
yi3_tv8 = CORRECTOR(name="yi3_tv8", len=0.0, xkick=0.0, ykick=-0.0)
yi3_oct8 = thinMULTIPOLE(name="yi3_oct8", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yi3_dec8 = thinMULTIPOLE(name="yi3_dec8", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yi3_qs8 = thinMULTIPOLE(name="yi3_qs8", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_939 = DRIFT(name="drift_939", len=0.5452500000001237)
yi3_qd8 = KQUAD(name="yi3_qd8", len=1.11, k1=-0.08049907581266105)
drift_940 = DRIFT(name="drift_940", len=0.29610560000037367)
yi3_b8 = MARKER(name="yi3_b8")
drift_941 = DRIFT(name="drift_941", len=1.3184179478494116)
yi3_hlx7_4 = DRIFT(name="yi3_hlx7_4", len=2.4)
drift_942 = DRIFT(name="drift_942", len=0.21200000000044383)
yi3_hlx7_3 = DRIFT(name="yi3_hlx7_3", len=2.4)
drift_943 = DRIFT(name="drift_943", len=0.22400000000016007)
yi3_b7_1 = MARKER(name="yi3_b7_1")
drift_944 = DRIFT(name="drift_944", len=0.22400000000016007)
yi3_hlx7_2 = DRIFT(name="yi3_hlx7_2", len=2.4)
drift_945 = DRIFT(name="drift_945", len=0.21200000000044383)
yi3_hlx7_1 = DRIFT(name="yi3_hlx7_1", len=2.4)
drift_946 = DRIFT(name="drift_946", len=1.3184179478494116)
yi3_b7 = MARKER(name="yi3_b7")
drift_947 = DRIFT(name="drift_947", len=0.2962336000000505)
yi3_qf7 = KQUAD(name="yi3_qf7", len=0.929744, k1=0.08736049017101663)
drift_948 = DRIFT(name="drift_948", len=0.5453779999998005)
yi3_th7 = CORRECTOR(name="yi3_th7", len=0.0, xkick=-0.0, ykick=0.0)
yi3_oct7 = thinMULTIPOLE(name="yi3_oct7", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yi3_dec7 = thinMULTIPOLE(name="yi3_dec7", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yi3_qgt7 = thinMULTIPOLE(name="yi3_qgt7", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, -0.0, 0.0, 0.0])
drift_949 = DRIFT(name="drift_949", len=1.5789970232644919)
yi3_dh6 = ESBEND(name="yi3_dh6", len=2.94942686017, angle=0.012160550077129212, e1=0.0, e2=0.0)
drift_950 = DRIFT(name="drift_950", len=8.076650933295241)
yi3_tv6 = CORRECTOR(name="yi3_tv6", len=0.0, xkick=0.0, ykick=-0.0)
yi3_oct6 = thinMULTIPOLE(name="yi3_oct6", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yi3_dec6 = thinMULTIPOLE(name="yi3_dec6", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yi3_qs6 = thinMULTIPOLE(name="yi3_qs6", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_951 = DRIFT(name="drift_951", len=0.5452500000001237)
yi3_qd6 = KQUAD(name="yi3_qd6", len=1.11, k1=-0.0896136513832874)
drift_952 = DRIFT(name="drift_952", len=0.13104999999995925)
yi3_tq6 = KQUAD(name="yi3_tq6", len=0.75, k1=-0.002379647064218151)
drift_953 = DRIFT(name="drift_953", len=0.22622760000012931)
yi3_bv6 = MARKER(name="yi3_bv6")
drift_954 = DRIFT(name="drift_954", len=2.277768451602242)
yi3_dh5 = ESBEND(name="yi3_dh5", len=6.915999880527922, angle=0.028514815544784158, e1=0.0, e2=0.0)
drift_955 = DRIFT(name="drift_955", len=2.839796002703224)
yi3_th5 = CORRECTOR(name="yi3_th5", len=0.0, xkick=-0.0, ykick=0.0)
yi3_oct5 = thinMULTIPOLE(name="yi3_oct5", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yi3_dec5 = thinMULTIPOLE(name="yi3_dec5", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yi3_qgt5 = thinMULTIPOLE(name="yi3_qgt5", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, -0.0, 0.0, 0.0])
drift_956 = DRIFT(name="drift_956", len=0.5452500000001237)
yi3_qf5 = KQUAD(name="yi3_qf5", len=1.11, k1=0.0896136513832874)
drift_957 = DRIFT(name="drift_957", len=0.13104999999995925)
yi3_tq5 = KQUAD(name="yi3_tq5", len=0.75, k1=0.003958379209394966)
drift_958 = DRIFT(name="drift_958", len=0.22622760000012931)
yi3_bh5 = MARKER(name="yi3_bh5")
drift_959 = DRIFT(name="drift_959", len=4.367472400000224)
yi3_tv4 = CORRECTOR(name="yi3_tv4", len=0.0, xkick=0.0, ykick=-0.0)
yi3_oct4 = thinMULTIPOLE(name="yi3_oct4", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yi3_dec4 = thinMULTIPOLE(name="yi3_dec4", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yi3_qs4 = thinMULTIPOLE(name="yi3_qs4", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_960 = DRIFT(name="drift_960", len=0.5442755000003672)
yi3_qd4 = KQUAD(name="yi3_qd4", len=1.811949, k1=-0.08993340871579382)
drift_961 = DRIFT(name="drift_961", len=0.13007550000020274)
yi3_tq4 = KQUAD(name="yi3_tq4", len=0.75, k1=0.03462226952510231)
drift_962 = DRIFT(name="drift_962", len=0.25163960000008956)
yi3_b4 = MARKER(name="yi3_b4")
drift_963 = DRIFT(name="drift_963", len=0.6552980000005846)
yi3_sv4 = MARKER(name="yi3_sv4")
drift_964 = DRIFT(name="drift_964", len=1.9598239999995712)
yi3_kscv3 = MARKER(name="yi3_kscv3")
drift_965 = DRIFT(name="drift_965", len=2.8562300000003233)
yi3_ksch3_2 = MARKER(name="yi3_ksch3_2")
drift_966 = DRIFT(name="drift_966", len=1.2626340000006167)
yi3_ksch3_1 = MARKER(name="yi3_ksch3_1")
drift_967 = DRIFT(name="drift_967", len=5.0)
cav197_5 = RFCA(name="cav197_5", len=0.874, volt=0.0, freq=197051693.13582847, energy=275000000000.0, h=2520.0, lag=0.0)
drift_968 = DRIFT(name="drift_968", len=0.5)
cav197_4 = RFCA(name="cav197_4", len=0.874, volt=0.0, freq=197051693.13582847, energy=275000000000.0, h=2520.0, lag=0.0)
drift_969 = DRIFT(name="drift_969", len=0.5)
cav197_3 = RFCA(name="cav197_3", len=0.874, volt=0.0, freq=197051693.13582847, energy=275000000000.0, h=2520.0, lag=0.0)
drift_970 = DRIFT(name="drift_970", len=0.5)
cav197_2 = RFCA(name="cav197_2", len=0.874, volt=0.0, freq=197051693.13582847, energy=275000000000.0, h=2520.0, lag=0.0)
drift_971 = DRIFT(name="drift_971", len=0.5)
cav197_1 = RFCA(name="cav197_1", len=0.874, volt=0.0, freq=197051693.13582847, energy=275000000000.0, h=2520.0, lag=0.0)
drift_972 = DRIFT(name="drift_972", len=3.778172999999697)
cav28_2 = RFCA(name="cav28_2", len=3.3948, volt=0.0, freq=24631461.64197856, energy=275000000000.0, h=315.0, lag=0.0)
drift_973 = DRIFT(name="drift_973", len=0.5)
cav28_1 = RFCA(name="cav28_1", len=3.3948, volt=0.0, freq=24631461.64197856, energy=275000000000.0, h=315.0, lag=0.0)
drift_974 = DRIFT(name="drift_974", len=3.651230999999825)
yi3_kfbh3 = CORRECTOR(name="yi3_kfbh3", len=0.0, xkick=0.0, ykick=0.0)
drift_975 = DRIFT(name="drift_975", len=0.7141814230244563)
warm_KKQUAD3 = KQUAD(name="warm_KKQUAD3", len=0.9, k1=0.0)
drift_976 = DRIFT(name="drift_976", len=0.24876899999981106)
yi3_sv3 = MARKER(name="yi3_sv3")
drift_977 = DRIFT(name="drift_977", len=1.587789000000157)
yi3_b3 = MARKER(name="yi3_b3")
drift_978 = DRIFT(name="drift_978", len=0.4714282599998114)
yi3_th3 = CORRECTOR(name="yi3_th3", len=0.0, xkick=-0.0, ykick=0.0)
yi3_sx3 = thinMULTIPOLE(name="yi3_sx3", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, -0.0, 0.0])
yi3_oct3 = thinMULTIPOLE(name="yi3_oct3", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yi3_dod3 = thinMULTIPOLE(name="yi3_dod3", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_979 = DRIFT(name="drift_979", len=0.5979680000000371)
yi3_qf3 = KQUAD(name="yi3_qf3", len=2.100484, k1=0.05424917393519864)
drift_980 = DRIFT(name="drift_980", len=0.48977800000011484)
yi3_dods3 = thinMULTIPOLE(name="yi3_dods3", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yi3_octs3 = thinMULTIPOLE(name="yi3_octs3", PolynomA=[0.0, 0.0, 0.0, -0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yi3_sxs3 = thinMULTIPOLE(name="yi3_sxs3", PolynomA=[0.0, 0.0, -0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yi3_qs3 = thinMULTIPOLE(name="yi3_qs3", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_981 = DRIFT(name="drift_981", len=1.3709614999997939)
yi3_qd2 = KQUAD(name="yi3_qd2", len=3.391633, k1=-0.05493718049433041)
drift_982 = DRIFT(name="drift_982", len=0.49420350000036706)
yi3_tv2 = CORRECTOR(name="yi3_tv2", len=0.0, xkick=0.0, ykick=-0.0)
yi3_oct2 = thinMULTIPOLE(name="yi3_oct2", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yi3_dec2 = thinMULTIPOLE(name="yi3_dec2", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yi3_dod2 = thinMULTIPOLE(name="yi3_dod2", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_983 = DRIFT(name="drift_983", len=1.1824280000000726)
yi3_qf1 = KQUAD(name="yi3_qf1", len=1.44, k1=0.05577876407684902)
drift_984 = DRIFT(name="drift_984", len=0.33694480000031035)
yi3_b1 = MARKER(name="yi3_b1")
drift_985 = DRIFT(name="drift_985", len=2.9277037600004405)
dwarm3 = ESBEND(name="dwarm3", len=5.000007675986592, angle=-0.006069983199868808, e1=-0.003034991599934404, e2=-0.003034991599934404)
drift_986 = DRIFT(name="drift_986", len=2.995289165633949)
sl_kick_mod1b = CORRECTOR(name="sl_kick_mod1b", len=0.9, xkick=0.0, ykick=0.0)
drift_987 = DRIFT(name="drift_987", len=0.1000000000003638)
sl_kick_mod2b = CORRECTOR(name="sl_kick_mod2b", len=0.9, xkick=0.0, ykick=0.0)
drift_988 = DRIFT(name="drift_988", len=0.1000000000003638)
sl_kick_mod3b = CORRECTOR(name="sl_kick_mod3b", len=0.9, xkick=0.0, ykick=0.0)
drift_989 = DRIFT(name="drift_989", len=0.1000000000003638)
sl_kick_mod4b = CORRECTOR(name="sl_kick_mod4b", len=0.9, xkick=0.0, ykick=0.0)
drift_990 = DRIFT(name="drift_990", len=0.1000000000003638)
sl_kick_mod5b = CORRECTOR(name="sl_kick_mod5b", len=0.9, xkick=0.0, ykick=0.0)
drift_991 = DRIFT(name="drift_991", len=0.1000000000003638)
sl_kick_mod6b = CORRECTOR(name="sl_kick_mod6b", len=0.9, xkick=0.0, ykick=0.0)
drift_992 = DRIFT(name="drift_992", len=0.1000000000003638)
sl_kick_mod7b = CORRECTOR(name="sl_kick_mod7b", len=0.9, xkick=0.0, ykick=0.0)
drift_993 = DRIFT(name="drift_993", len=0.1000000000003638)
sl_kick_mod8b = CORRECTOR(name="sl_kick_mod8b", len=0.9, xkick=0.0, ykick=0.0)
drift_994 = DRIFT(name="drift_994", len=0.1000000000003638)
sl_kick_mod9b = CORRECTOR(name="sl_kick_mod9b", len=0.9, xkick=0.0, ykick=0.0)
drift_995 = DRIFT(name="drift_995", len=0.1000000000003638)
sl_kick_mod10b = CORRECTOR(name="sl_kick_mod10b", len=0.9, xkick=0.0, ykick=0.0)
drift_996 = DRIFT(name="drift_996", len=4.200000000000728)
ip4 = MARKER(name="ip4")
drift_997 = DRIFT(name="drift_997", len=4.0)
sl_kick_mod1a = CORRECTOR(name="sl_kick_mod1a", len=0.9, xkick=0.0, ykick=0.0)
drift_998 = DRIFT(name="drift_998", len=0.1000000000003638)
sl_kick_mod2a = CORRECTOR(name="sl_kick_mod2a", len=0.9, xkick=0.0, ykick=0.0)
drift_999 = DRIFT(name="drift_999", len=0.1000000000003638)
sl_kick_mod3a = CORRECTOR(name="sl_kick_mod3a", len=0.9, xkick=0.0, ykick=0.0)
drift_1000 = DRIFT(name="drift_1000", len=0.1000000000003638)
sl_kick_mod4a = CORRECTOR(name="sl_kick_mod4a", len=0.9, xkick=0.0, ykick=0.0)
drift_1001 = DRIFT(name="drift_1001", len=0.1000000000003638)
sl_kick_mod5a = CORRECTOR(name="sl_kick_mod5a", len=0.9, xkick=0.0, ykick=0.0)
drift_1002 = DRIFT(name="drift_1002", len=0.1000000000003638)
sl_kick_mod6a = CORRECTOR(name="sl_kick_mod6a", len=0.9, xkick=0.0, ykick=0.0)
drift_1003 = DRIFT(name="drift_1003", len=0.1000000000003638)
sl_kick_mod7a = CORRECTOR(name="sl_kick_mod7a", len=0.9, xkick=0.0, ykick=0.0)
drift_1004 = DRIFT(name="drift_1004", len=0.1000000000003638)
sl_kick_mod8a = CORRECTOR(name="sl_kick_mod8a", len=0.9, xkick=0.0, ykick=0.0)
drift_1005 = DRIFT(name="drift_1005", len=0.1000000000003638)
sl_kick_mod9a = CORRECTOR(name="sl_kick_mod9a", len=0.9, xkick=0.0, ykick=0.0)
drift_1006 = DRIFT(name="drift_1006", len=0.1000000000003638)
sl_kick_mod10a = CORRECTOR(name="sl_kick_mod10a", len=0.9, xkick=0.0, ykick=0.0)
drift_1007 = DRIFT(name="drift_1007", len=3.1952891656346765)
dwarm4 = ESBEND(name="dwarm4", len=5.000007675986592, angle=0.006069983199868808, e1=0.003034991599934404, e2=0.003034991599934404)
drift_1008 = DRIFT(name="drift_1008", len=0.36306468000020686)
det_pol_hpol_hjet_1 = DRIFT(name="det_pol_hpol_hjet_1", len=1.0)
det_pol_hpol_hjet = MARKER(name="det_pol_hpol_hjet")
drift_1009 = DRIFT(name="drift_1009", len=0.5646390800002337)
yo4_b1 = MARKER(name="yo4_b1")
drift_1010 = DRIFT(name="drift_1010", len=0.33694480000031035)
yo4_qd1 = KQUAD(name="yo4_qd1", len=1.44, k1=-0.05706610893745024)
drift_1011 = DRIFT(name="drift_1011", len=1.1824280000000726)
yo4_th2 = CORRECTOR(name="yo4_th2", len=0.0, xkick=-0.0, ykick=0.0)
yo4_oct2 = thinMULTIPOLE(name="yo4_oct2", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yo4_dec2 = thinMULTIPOLE(name="yo4_dec2", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yo4_dod2 = thinMULTIPOLE(name="yo4_dod2", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_1012 = DRIFT(name="drift_1012", len=0.49420350000036706)
yo4_qf2 = KQUAD(name="yo4_qf2", len=3.391633, k1=0.056219455402784016)
drift_1013 = DRIFT(name="drift_1013", len=1.3709614999997939)
yo4_dods3 = thinMULTIPOLE(name="yo4_dods3", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yo4_octs3 = thinMULTIPOLE(name="yo4_octs3", PolynomA=[0.0, 0.0, 0.0, -0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yo4_sxs3 = thinMULTIPOLE(name="yo4_sxs3", PolynomA=[0.0, 0.0, -0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yo4_qs3 = thinMULTIPOLE(name="yo4_qs3", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_1014 = DRIFT(name="drift_1014", len=0.48977800000011484)
yo4_qd3 = KQUAD(name="yo4_qd3", len=2.100484, k1=-0.05485224779864115)
drift_1015 = DRIFT(name="drift_1015", len=0.5979680000000371)
yo4_tv3 = CORRECTOR(name="yo4_tv3", len=0.0, xkick=0.0, ykick=-0.0)
yo4_sx3 = thinMULTIPOLE(name="yo4_sx3", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, -0.0, 0.0])
yo4_oct3 = thinMULTIPOLE(name="yo4_oct3", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yo4_dod3 = thinMULTIPOLE(name="yo4_dod3", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_1016 = DRIFT(name="drift_1016", len=0.4714282599998114)
yo4_b3 = MARKER(name="yo4_b3")
drift_1017 = DRIFT(name="drift_1017", len=1.587786000000051)
yo4_sv3_1 = MARKER(name="yo4_sv3_1")
drift_1018 = DRIFT(name="drift_1018", len=0.24876599999970495)
warm_KKQUAD4 = KQUAD(name="warm_KKQUAD4", len=0.9, k1=0.0)
drift_1019 = DRIFT(name="drift_1019", len=0.7146169999996346)
yo4_kfbh3 = CORRECTOR(name="yo4_kfbh3", len=0.0, xkick=0.0, ykick=0.0)
drift_1020 = DRIFT(name="drift_1020", len=1.257002499999544)
yo4_wcm3 = MARKER(name="yo4_wcm3")
drift_1021 = DRIFT(name="drift_1021", len=16.0)
septum_ir4 = DRIFT(name="septum_ir4", len=2.0)
drift_1022 = DRIFT(name="drift_1022", len=5.315557923024244)
det_pol_hpol_pc2_1 = DRIFT(name="det_pol_hpol_pc2_1", len=1.25)
det_pol_hpol_pc2 = MARKER(name="det_pol_hpol_pc2")
drift_1023 = DRIFT(name="drift_1023", len=0.5)
det_pol_hpol_pc1_1 = DRIFT(name="det_pol_hpol_pc1_1", len=1.25)
det_pol_hpol_pc1 = MARKER(name="det_pol_hpol_pc1")
drift_1024 = DRIFT(name="drift_1024", len=2.0947100000003047)
yo4_sv4 = MARKER(name="yo4_sv4")
drift_1025 = DRIFT(name="drift_1025", len=0.6552899999996953)
yo4_b4 = MARKER(name="yo4_b4")
drift_1026 = DRIFT(name="drift_1026", len=0.25163960000008956)
yo4_tq4 = KQUAD(name="yo4_tq4", len=0.75, k1=-0.009170217484177364)
drift_1027 = DRIFT(name="drift_1027", len=0.13007550000020274)
yo4_qf4 = KQUAD(name="yo4_qf4", len=1.811949, k1=0.08993340871579382)
drift_1028 = DRIFT(name="drift_1028", len=0.5442755000003672)
yo4_th4 = CORRECTOR(name="yo4_th4", len=0.0, xkick=-0.0, ykick=0.0)
yo4_oct4 = thinMULTIPOLE(name="yo4_oct4", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yo4_dec4 = thinMULTIPOLE(name="yo4_dec4", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yo4_qs4 = thinMULTIPOLE(name="yo4_qs4", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_1029 = DRIFT(name="drift_1029", len=4.367472400000224)
yo4_bv5 = MARKER(name="yo4_bv5")
drift_1030 = DRIFT(name="drift_1030", len=0.22622760000012931)
yo4_tq5 = KQUAD(name="yo4_tq5", len=0.75, k1=-0.016089282251060044)
drift_1031 = DRIFT(name="drift_1031", len=0.13104999999995925)
yo4_qd5 = KQUAD(name="yo4_qd5", len=1.11, k1=-0.0896136513832874)
drift_1032 = DRIFT(name="drift_1032", len=0.5452500000001237)
yo4_tv5 = CORRECTOR(name="yo4_tv5", len=0.0, xkick=0.0, ykick=-0.0)
yo4_oct5 = thinMULTIPOLE(name="yo4_oct5", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yo4_dec5 = thinMULTIPOLE(name="yo4_dec5", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yo4_qs5 = thinMULTIPOLE(name="yo4_qs5", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_1033 = DRIFT(name="drift_1033", len=1.9628259081100623)
yo4_dh5 = ESBEND(name="yo4_dh5", len=8.698449375192075, angle=0.035863892964716426, e1=0.0, e2=0.0)
drift_1034 = DRIFT(name="drift_1034", len=1.4007982592102053)
yo4_bh6 = MARKER(name="yo4_bh6")
drift_1035 = DRIFT(name="drift_1035", len=0.22622760000012931)
yo4_tq6 = KQUAD(name="yo4_tq6", len=0.75, k1=0.02312258002485002)
drift_1036 = DRIFT(name="drift_1036", len=0.13104999999995925)
yo4_qf6 = KQUAD(name="yo4_qf6", len=1.11, k1=0.0896136513832874)
drift_1037 = DRIFT(name="drift_1037", len=0.5452500000001237)
yo4_th6 = CORRECTOR(name="yo4_th6", len=0.0, xkick=-0.0, ykick=0.0)
yo4_oct6 = thinMULTIPOLE(name="yo4_oct6", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yo4_dec6 = thinMULTIPOLE(name="yo4_dec6", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yo4_qgt6 = thinMULTIPOLE(name="yo4_qgt6", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, -0.0, 0.0, 0.0])
drift_1038 = DRIFT(name="drift_1038", len=8.082123248266726)
yo4_dh6 = ESBEND(name="yo4_dh6", len=2.94942686017, angle=0.012160550077129212, e1=0.0, e2=0.0)
drift_1039 = DRIFT(name="drift_1039", len=1.5844693382359765)
yo4_tv7 = CORRECTOR(name="yo4_tv7", len=0.0, xkick=0.0, ykick=-0.0)
yo4_oct7 = thinMULTIPOLE(name="yo4_oct7", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yo4_dec7 = thinMULTIPOLE(name="yo4_dec7", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yo4_qs7 = thinMULTIPOLE(name="yo4_qs7", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_1040 = DRIFT(name="drift_1040", len=0.5453779999998005)
yo4_qd7 = KQUAD(name="yo4_qd7", len=0.929744, k1=-0.08736049017101663)
drift_1041 = DRIFT(name="drift_1041", len=0.2962336000000505)
yo4_b7 = MARKER(name="yo4_b7")
drift_1042 = DRIFT(name="drift_1042", len=13.108835895699485)
yo4_b8 = MARKER(name="yo4_b8")
drift_1043 = DRIFT(name="drift_1043", len=0.29610560000037367)
yo4_qf8 = KQUAD(name="yo4_qf8", len=1.11, k1=0.07952250035938109)
drift_1044 = DRIFT(name="drift_1044", len=0.5452500000001237)
yo4_th8 = CORRECTOR(name="yo4_th8", len=0.0, xkick=-0.0, ykick=0.0)
yo4_oct8 = thinMULTIPOLE(name="yo4_oct8", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yo4_dec8 = thinMULTIPOLE(name="yo4_dec8", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yo4_qgt8 = thinMULTIPOLE(name="yo4_qgt8", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, -0.0, 0.0, 0.0])
drift_1045 = DRIFT(name="drift_1045", len=1.593705149733978)
yo4_dh8 = ESBEND(name="yo4_dh8", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_1046 = DRIFT(name="drift_1046", len=1.593705149733978)
yo4_tv9 = CORRECTOR(name="yo4_tv9", len=0.0, xkick=0.0, ykick=-0.0)
yo4_oct9 = thinMULTIPOLE(name="yo4_oct9", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yo4_dec9 = thinMULTIPOLE(name="yo4_dec9", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yo4_qs9 = thinMULTIPOLE(name="yo4_qs9", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_1047 = DRIFT(name="drift_1047", len=0.5452500000001237)
yo4_qd9 = KQUAD(name="yo4_qd9", len=1.11, k1=-0.08049907581266105)
drift_1048 = DRIFT(name="drift_1048", len=0.13104999999995925)
yo4_sxd9 = KSEXT(name="yo4_sxd9", len=0.75, k2=-0.6290818079014571)
drift_1049 = DRIFT(name="drift_1049", len=0.22622760000012931)
yo4_bv9 = MARKER(name="yo4_bv9")
drift_1050 = DRIFT(name="drift_1050", len=7.520095648265851)
yo4_dh9 = ESBEND(name="yo4_dh9", len=2.94942686017, angle=0.012160550077129212, e1=0.0, e2=0.0)
drift_1051 = DRIFT(name="drift_1051", len=1.5844693382359765)
yo4_th10 = CORRECTOR(name="yo4_th10", len=0.0, xkick=-0.0, ykick=0.0)
yo4_oct10 = thinMULTIPOLE(name="yo4_oct10", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yo4_dec10 = thinMULTIPOLE(name="yo4_dec10", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yo4_qs10 = thinMULTIPOLE(name="yo4_qs10", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_1052 = DRIFT(name="drift_1052", len=0.5452500000001237)
yo4_qf10 = KQUAD(name="yo4_qf10", len=1.11, k1=0.08048071257602774)
drift_1053 = DRIFT(name="drift_1053", len=0.13104999999995925)
yo4_sxf10 = KSEXT(name="yo4_sxf10", len=0.75, k2=0.39909333141434805)
drift_1054 = DRIFT(name="drift_1054", len=0.22622760000012931)
yo4_bh10 = MARKER(name="yo4_bh10")
drift_1055 = DRIFT(name="drift_1055", len=1.0316775497331037)
yo4_dh10 = ESBEND(name="yo4_dh10", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_1056 = DRIFT(name="drift_1056", len=1.343705149733978)
yo4_tv11 = CORRECTOR(name="yo4_tv11", len=0.5, xkick=0.0, ykick=0.0)
drift_1057 = DRIFT(name="drift_1057", len=0.2952500000001237)
yo4_qd11 = KQUAD(name="yo4_qd11", len=1.11, k1=-0.08364335555329494)
drift_1058 = DRIFT(name="drift_1058", len=0.13104999999995925)
yo4_sxd11 = KSEXT(name="yo4_sxd11", len=0.75, k2=-0.6290818079014571)
drift_1059 = DRIFT(name="drift_1059", len=0.22622760000012931)
yo4_bv11 = MARKER(name="yo4_bv11")
drift_1060 = DRIFT(name="drift_1060", len=1.0316775497331037)
yo4_dh11 = ESBEND(name="yo4_dh11", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_1061 = DRIFT(name="drift_1061", len=1.343705149733978)
yo4_th12 = CORRECTOR(name="yo4_th12", len=0.5, xkick=0.0, ykick=0.0)
drift_1062 = DRIFT(name="drift_1062", len=0.2952500000001237)
yo4_qf12 = KQUAD(name="yo4_qf12", len=1.11, k1=0.08048071257602774)
drift_1063 = DRIFT(name="drift_1063", len=0.13104999999995925)
yo4_sxf12 = KSEXT(name="yo4_sxf12", len=0.75, k2=0.39909333141434805)
drift_1064 = DRIFT(name="drift_1064", len=0.22622760000012931)
yo4_bh12 = MARKER(name="yo4_bh12")
drift_1065 = DRIFT(name="drift_1065", len=1.0316775497331037)
yo4_dh12 = ESBEND(name="yo4_dh12", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_1066 = DRIFT(name="drift_1066", len=1.343705149733978)
yo4_tv13 = CORRECTOR(name="yo4_tv13", len=0.5, xkick=0.0, ykick=0.0)
drift_1067 = DRIFT(name="drift_1067", len=0.2952500000001237)
yo4_qd13 = KQUAD(name="yo4_qd13", len=1.11, k1=-0.08364335555329494)
drift_1068 = DRIFT(name="drift_1068", len=0.13104999999995925)
yo4_sxd13 = KSEXT(name="yo4_sxd13", len=0.75, k2=-0.6290818079014571)
drift_1069 = DRIFT(name="drift_1069", len=0.22622760000012931)
yo4_bv13 = MARKER(name="yo4_bv13")
drift_1070 = DRIFT(name="drift_1070", len=1.0316775497331037)
yo4_dh13 = ESBEND(name="yo4_dh13", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_1071 = DRIFT(name="drift_1071", len=1.343705149733978)
yo4_th14 = CORRECTOR(name="yo4_th14", len=0.5, xkick=0.0, ykick=0.0)
drift_1072 = DRIFT(name="drift_1072", len=0.2952500000001237)
yo4_qf14 = KQUAD(name="yo4_qf14", len=1.11, k1=0.08048071257602774)
drift_1073 = DRIFT(name="drift_1073", len=0.13104999999995925)
yo4_sxf14 = KSEXT(name="yo4_sxf14", len=0.75, k2=0.39909333141434805)
drift_1074 = DRIFT(name="drift_1074", len=0.22622760000012931)
yo4_bh14 = MARKER(name="yo4_bh14")
drift_1075 = DRIFT(name="drift_1075", len=1.0316775497331037)
yo4_dh14 = ESBEND(name="yo4_dh14", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_1076 = DRIFT(name="drift_1076", len=1.343705149733978)
yo4_tv15 = CORRECTOR(name="yo4_tv15", len=0.5, xkick=0.0, ykick=0.0)
drift_1077 = DRIFT(name="drift_1077", len=0.2952500000001237)
yo4_qd15 = KQUAD(name="yo4_qd15", len=1.11, k1=-0.08364335555329494)
drift_1078 = DRIFT(name="drift_1078", len=0.13104999999995925)
yo4_sxd15 = KSEXT(name="yo4_sxd15", len=0.75, k2=-0.6290818079014571)
drift_1079 = DRIFT(name="drift_1079", len=0.22622760000012931)
yo4_bv15 = MARKER(name="yo4_bv15")
drift_1080 = DRIFT(name="drift_1080", len=1.0316775497331037)
yo4_dh15 = ESBEND(name="yo4_dh15", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_1081 = DRIFT(name="drift_1081", len=1.343705149733978)
yo4_th16 = CORRECTOR(name="yo4_th16", len=0.5, xkick=0.0, ykick=0.0)
drift_1082 = DRIFT(name="drift_1082", len=0.2952500000001237)
yo4_qf16 = KQUAD(name="yo4_qf16", len=1.11, k1=0.08048071257602774)
drift_1083 = DRIFT(name="drift_1083", len=0.13104999999995925)
yo4_sxf16 = KSEXT(name="yo4_sxf16", len=0.75, k2=0.39909333141434805)
drift_1084 = DRIFT(name="drift_1084", len=0.22622760000012931)
yo4_bh16 = MARKER(name="yo4_bh16")
drift_1085 = DRIFT(name="drift_1085", len=1.0316775497331037)
yo4_dh16 = ESBEND(name="yo4_dh16", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_1086 = DRIFT(name="drift_1086", len=1.343705149733978)
yo4_tv17 = CORRECTOR(name="yo4_tv17", len=0.5, xkick=0.0, ykick=0.0)
drift_1087 = DRIFT(name="drift_1087", len=0.2952500000001237)
yo4_qd17 = KQUAD(name="yo4_qd17", len=1.11, k1=-0.08364335555329494)
drift_1088 = DRIFT(name="drift_1088", len=0.13104999999995925)
yo4_sxd17 = KSEXT(name="yo4_sxd17", len=0.75, k2=-0.6290818079014571)
drift_1089 = DRIFT(name="drift_1089", len=0.22622760000012931)
yo4_bv17 = MARKER(name="yo4_bv17")
drift_1090 = DRIFT(name="drift_1090", len=1.0316775497331037)
yo4_dh17 = ESBEND(name="yo4_dh17", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_1091 = DRIFT(name="drift_1091", len=1.343705149733978)
yo4_th18 = CORRECTOR(name="yo4_th18", len=0.5, xkick=0.0, ykick=0.0)
drift_1092 = DRIFT(name="drift_1092", len=0.2952500000001237)
yo4_qf18 = KQUAD(name="yo4_qf18", len=1.11, k1=0.08048071257602774)
drift_1093 = DRIFT(name="drift_1093", len=0.13104999999995925)
yo4_sxf18 = KSEXT(name="yo4_sxf18", len=0.75, k2=0.39909333141434805)
drift_1094 = DRIFT(name="drift_1094", len=0.22622760000012931)
yo4_bh18 = MARKER(name="yo4_bh18")
drift_1095 = DRIFT(name="drift_1095", len=1.0316775497331037)
yo4_dh18 = ESBEND(name="yo4_dh18", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_1096 = DRIFT(name="drift_1096", len=1.343705149733978)
yo4_tv19 = CORRECTOR(name="yo4_tv19", len=0.5, xkick=0.0, ykick=0.0)
drift_1097 = DRIFT(name="drift_1097", len=0.2952500000001237)
yo4_qd19 = KQUAD(name="yo4_qd19", len=1.11, k1=-0.08364335555329494)
drift_1098 = DRIFT(name="drift_1098", len=0.13104999999995925)
yo4_sxd19 = KSEXT(name="yo4_sxd19", len=0.75, k2=-0.6290818079014571)
drift_1099 = DRIFT(name="drift_1099", len=0.22622760000012931)
yo4_bv19 = MARKER(name="yo4_bv19")
drift_1100 = DRIFT(name="drift_1100", len=1.0316775497331037)
yo4_dh19 = ESBEND(name="yo4_dh19", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_1101 = DRIFT(name="drift_1101", len=1.343705149733978)
yo4_th20 = CORRECTOR(name="yo4_th20", len=0.5, xkick=0.0, ykick=0.0)
drift_1102 = DRIFT(name="drift_1102", len=0.2952500000001237)
yo4_qf20 = KQUAD(name="yo4_qf20", len=1.11, k1=0.08048071257602774)
drift_1103 = DRIFT(name="drift_1103", len=0.13104999999995925)
yo4_sxf20 = KSEXT(name="yo4_sxf20", len=0.75, k2=0.39909333141434805)
drift_1104 = DRIFT(name="drift_1104", len=0.22622760000012931)
yo4_bh20 = MARKER(name="yo4_bh20")
drift_1105 = DRIFT(name="drift_1105", len=1.0316775497331037)
yo4_dh20 = ESBEND(name="yo4_dh20", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_1106 = DRIFT(name="drift_1106", len=1.343705149733978)
yo5_tv21 = CORRECTOR(name="yo5_tv21", len=0.5, xkick=0.0, ykick=0.0)
drift_1107 = DRIFT(name="drift_1107", len=0.2952500000001237)
yo5_qd21 = KQUAD(name="yo5_qd21", len=1.11, k1=-0.08364335555329494)
drift_1108 = DRIFT(name="drift_1108", len=0.13104999999995925)
yo5_sxd21 = KSEXT(name="yo5_sxd21", len=0.75, k2=-0.6290818079014571)
drift_1109 = DRIFT(name="drift_1109", len=0.22622760000012931)
yo5_bv21 = MARKER(name="yo5_bv21")
drift_1110 = DRIFT(name="drift_1110", len=1.0316775497331037)
yo5_dh20 = ESBEND(name="yo5_dh20", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_1111 = DRIFT(name="drift_1111", len=1.343705149733978)
yo5_th20 = CORRECTOR(name="yo5_th20", len=0.5, xkick=0.0, ykick=0.0)
drift_1112 = DRIFT(name="drift_1112", len=0.2952500000001237)
yo5_qf20 = KQUAD(name="yo5_qf20", len=1.11, k1=0.08048071257602774)
drift_1113 = DRIFT(name="drift_1113", len=0.13104999999995925)
yo5_sxf20 = KSEXT(name="yo5_sxf20", len=0.75, k2=0.39909333141434805)
drift_1114 = DRIFT(name="drift_1114", len=0.22622760000012931)
yo5_bh20 = MARKER(name="yo5_bh20")
drift_1115 = DRIFT(name="drift_1115", len=1.0316775497331037)
yo5_dh19 = ESBEND(name="yo5_dh19", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_1116 = DRIFT(name="drift_1116", len=1.343705149733978)
yo5_tv19 = CORRECTOR(name="yo5_tv19", len=0.5, xkick=0.0, ykick=0.0)
drift_1117 = DRIFT(name="drift_1117", len=0.2952500000001237)
yo5_qd19 = KQUAD(name="yo5_qd19", len=1.11, k1=-0.08364335555329494)
drift_1118 = DRIFT(name="drift_1118", len=0.13104999999995925)
yo5_sxd19 = KSEXT(name="yo5_sxd19", len=0.75, k2=-0.6290818079014571)
drift_1119 = DRIFT(name="drift_1119", len=0.22622760000012931)
yo5_bv19 = MARKER(name="yo5_bv19")
drift_1120 = DRIFT(name="drift_1120", len=1.0316775497331037)
yo5_dh18 = ESBEND(name="yo5_dh18", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_1121 = DRIFT(name="drift_1121", len=1.593705149733978)
yo5_th18 = CORRECTOR(name="yo5_th18", len=0.0, xkick=-0.0, ykick=0.0)
yo5_oct18 = thinMULTIPOLE(name="yo5_oct18", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yo5_dec18 = thinMULTIPOLE(name="yo5_dec18", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yo5_qgt18 = thinMULTIPOLE(name="yo5_qgt18", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, -0.0, 0.0, 0.0])
drift_1122 = DRIFT(name="drift_1122", len=0.5452500000001237)
yo5_qf18 = KQUAD(name="yo5_qf18", len=1.11, k1=0.08048071257602774)
drift_1123 = DRIFT(name="drift_1123", len=0.13104999999995925)
yo5_sxf18 = KSEXT(name="yo5_sxf18", len=0.75, k2=0.39909333141434805)
drift_1124 = DRIFT(name="drift_1124", len=0.22622760000012931)
yo5_bh18 = MARKER(name="yo5_bh18")
drift_1125 = DRIFT(name="drift_1125", len=1.0316775497331037)
yo5_dh17 = ESBEND(name="yo5_dh17", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_1126 = DRIFT(name="drift_1126", len=1.593705149733978)
yo5_tv17 = CORRECTOR(name="yo5_tv17", len=0.0, xkick=0.0, ykick=-0.0)
yo5_oct17 = thinMULTIPOLE(name="yo5_oct17", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yo5_dec17 = thinMULTIPOLE(name="yo5_dec17", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yo5_qs17 = thinMULTIPOLE(name="yo5_qs17", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_1127 = DRIFT(name="drift_1127", len=0.5452500000001237)
yo5_qd17 = KQUAD(name="yo5_qd17", len=1.11, k1=-0.08364335555329494)
drift_1128 = DRIFT(name="drift_1128", len=0.13104999999995925)
yo5_sxd17 = KSEXT(name="yo5_sxd17", len=0.75, k2=-0.6290818079014571)
drift_1129 = DRIFT(name="drift_1129", len=0.22622760000012931)
yo5_bv17 = MARKER(name="yo5_bv17")
drift_1130 = DRIFT(name="drift_1130", len=1.0316775497331037)
yo5_dh16 = ESBEND(name="yo5_dh16", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_1131 = DRIFT(name="drift_1131", len=1.593705149733978)
yo5_th16 = CORRECTOR(name="yo5_th16", len=0.0, xkick=-0.0, ykick=0.0)
yo5_oct16 = thinMULTIPOLE(name="yo5_oct16", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yo5_dec16 = thinMULTIPOLE(name="yo5_dec16", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yo5_qgt16 = thinMULTIPOLE(name="yo5_qgt16", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, -0.0, 0.0, 0.0])
drift_1132 = DRIFT(name="drift_1132", len=0.5452500000001237)
yo5_qf16 = KQUAD(name="yo5_qf16", len=1.11, k1=0.08048071257602774)
drift_1133 = DRIFT(name="drift_1133", len=0.13104999999995925)
yo5_sxf16 = KSEXT(name="yo5_sxf16", len=0.75, k2=0.39909333141434805)
drift_1134 = DRIFT(name="drift_1134", len=0.22622760000012931)
yo5_bh16 = MARKER(name="yo5_bh16")
drift_1135 = DRIFT(name="drift_1135", len=1.0316775497331037)
yo5_dh15 = ESBEND(name="yo5_dh15", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_1136 = DRIFT(name="drift_1136", len=1.593705149733978)
yo5_tv15 = CORRECTOR(name="yo5_tv15", len=0.0, xkick=0.0, ykick=-0.0)
yo5_oct15 = thinMULTIPOLE(name="yo5_oct15", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yo5_dec15 = thinMULTIPOLE(name="yo5_dec15", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yo5_qs15 = thinMULTIPOLE(name="yo5_qs15", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_1137 = DRIFT(name="drift_1137", len=0.5452500000001237)
yo5_qd15 = KQUAD(name="yo5_qd15", len=1.11, k1=-0.08364335555329494)
drift_1138 = DRIFT(name="drift_1138", len=0.13104999999995925)
yo5_sxd15 = KSEXT(name="yo5_sxd15", len=0.75, k2=-0.6290818079014571)
drift_1139 = DRIFT(name="drift_1139", len=0.22622760000012931)
yo5_bv15 = MARKER(name="yo5_bv15")
drift_1140 = DRIFT(name="drift_1140", len=1.0316775497331037)
yo5_dh14 = ESBEND(name="yo5_dh14", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_1141 = DRIFT(name="drift_1141", len=1.593705149733978)
yo5_th14 = CORRECTOR(name="yo5_th14", len=0.0, xkick=-0.0, ykick=0.0)
yo5_oct14 = thinMULTIPOLE(name="yo5_oct14", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yo5_dec14 = thinMULTIPOLE(name="yo5_dec14", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yo5_qgt14 = thinMULTIPOLE(name="yo5_qgt14", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, -0.0, 0.0, 0.0])
drift_1142 = DRIFT(name="drift_1142", len=0.5452500000001237)
yo5_qf14 = KQUAD(name="yo5_qf14", len=1.11, k1=0.08048071257602774)
drift_1143 = DRIFT(name="drift_1143", len=0.13104999999995925)
yo5_sxf14 = KSEXT(name="yo5_sxf14", len=0.75, k2=0.39909333141434805)
drift_1144 = DRIFT(name="drift_1144", len=0.22622760000012931)
yo5_bh14 = MARKER(name="yo5_bh14")
drift_1145 = DRIFT(name="drift_1145", len=1.0316775497331037)
yo5_dh13 = ESBEND(name="yo5_dh13", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_1146 = DRIFT(name="drift_1146", len=1.593705149733978)
yo5_tv13 = CORRECTOR(name="yo5_tv13", len=0.0, xkick=0.0, ykick=-0.0)
yo5_oct13 = thinMULTIPOLE(name="yo5_oct13", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yo5_dec13 = thinMULTIPOLE(name="yo5_dec13", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yo5_qs13 = thinMULTIPOLE(name="yo5_qs13", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_1147 = DRIFT(name="drift_1147", len=0.5452500000001237)
yo5_qd13 = KQUAD(name="yo5_qd13", len=1.11, k1=-0.08364335555329494)
drift_1148 = DRIFT(name="drift_1148", len=0.13104999999995925)
yo5_sxd13 = KSEXT(name="yo5_sxd13", len=0.75, k2=-0.6290818079014571)
drift_1149 = DRIFT(name="drift_1149", len=0.22622760000012931)
yo5_bv13 = MARKER(name="yo5_bv13")
drift_1150 = DRIFT(name="drift_1150", len=1.0316775497331037)
yo5_dh12 = ESBEND(name="yo5_dh12", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_1151 = DRIFT(name="drift_1151", len=1.593705149733978)
yo5_th12 = CORRECTOR(name="yo5_th12", len=0.0, xkick=-0.0, ykick=0.0)
yo5_oct12 = thinMULTIPOLE(name="yo5_oct12", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yo5_dec12 = thinMULTIPOLE(name="yo5_dec12", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yo5_qgt12 = thinMULTIPOLE(name="yo5_qgt12", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, -0.0, 0.0, 0.0])
drift_1152 = DRIFT(name="drift_1152", len=0.5452500000001237)
yo5_qf12 = KQUAD(name="yo5_qf12", len=1.11, k1=0.08048071257602774)
drift_1153 = DRIFT(name="drift_1153", len=0.13104999999995925)
yo5_sxf12 = KSEXT(name="yo5_sxf12", len=0.75, k2=0.39909333141434805)
drift_1154 = DRIFT(name="drift_1154", len=0.22622760000012931)
yo5_bh12 = MARKER(name="yo5_bh12")
drift_1155 = DRIFT(name="drift_1155", len=1.0316775497331037)
yo5_dh11 = ESBEND(name="yo5_dh11", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_1156 = DRIFT(name="drift_1156", len=1.593705149733978)
yo5_tv11 = CORRECTOR(name="yo5_tv11", len=0.0, xkick=0.0, ykick=-0.0)
yo5_oct11 = thinMULTIPOLE(name="yo5_oct11", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yo5_dec11 = thinMULTIPOLE(name="yo5_dec11", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yo5_qs11 = thinMULTIPOLE(name="yo5_qs11", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_1157 = DRIFT(name="drift_1157", len=0.5452500000001237)
yo5_qd11 = KQUAD(name="yo5_qd11", len=1.11, k1=-0.08364335555329494)
drift_1158 = DRIFT(name="drift_1158", len=0.13104999999995925)
yo5_sxd11 = KSEXT(name="yo5_sxd11", len=0.75, k2=-0.6290818079014571)
drift_1159 = DRIFT(name="drift_1159", len=0.22622760000012931)
yo5_bv11 = MARKER(name="yo5_bv11")
drift_1160 = DRIFT(name="drift_1160", len=1.0316775497331037)
yo5_dh10 = ESBEND(name="yo5_dh10", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_1161 = DRIFT(name="drift_1161", len=1.593705149733978)
yo5_th10 = CORRECTOR(name="yo5_th10", len=0.0, xkick=-0.0, ykick=0.0)
yo5_oct10 = thinMULTIPOLE(name="yo5_oct10", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yo5_dec10 = thinMULTIPOLE(name="yo5_dec10", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yo5_qs10 = thinMULTIPOLE(name="yo5_qs10", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_1162 = DRIFT(name="drift_1162", len=0.5452500000001237)
yo5_qf10 = KQUAD(name="yo5_qf10", len=1.11, k1=0.08048071257602774)
drift_1163 = DRIFT(name="drift_1163", len=0.13104999999995925)
yo5_sxf10 = KSEXT(name="yo5_sxf10", len=0.75, k2=0.39909333141434805)
drift_1164 = DRIFT(name="drift_1164", len=0.22622760000012931)
yo5_bh10 = MARKER(name="yo5_bh10")
drift_1165 = DRIFT(name="drift_1165", len=0.23988339999959862)
yo5_sv9_2 = MARKER(name="yo5_sv9_2")
drift_1166 = DRIFT(name="drift_1166", len=3.9725006582220885)
h5_tq11 = KQUAD(name="h5_tq11", len=0.75, k1=-0.0)
drift_1167 = DRIFT(name="drift_1167", len=0.13104999999995925)
h5_q11 = KQUAD(name="h5_q11", len=1.11, k1=-0.010124193523445024)
drift_1168 = DRIFT(name="drift_1168", len=0.5452500000001237)
h5_tv11 = CORRECTOR(name="h5_tv11", len=0.0, xkick=0.0, ykick=-0.0)
h5_b311 = thinMULTIPOLE(name="h5_b311", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
h5_b411 = thinMULTIPOLE(name="h5_b411", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
h5_a111 = thinMULTIPOLE(name="h5_a111", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_1169 = DRIFT(name="drift_1169", len=1.207257000000027)
yo5_sv9_1 = MARKER(name="yo5_sv9_1")
drift_1170 = DRIFT(name="drift_1170", len=1.0406167959017694)
yo5_dh9 = ESBEND(name="yo5_dh9", len=2.94942686017, angle=0.012160550077129212, e1=0.0, e2=0.0)
drift_1171 = DRIFT(name="drift_1171", len=1.4627375302361543)
yo5_bv9 = MARKER(name="yo5_bv9")
drift_1172 = DRIFT(name="drift_1172", len=0.29600159999972675)
yo5_qd9 = KQUAD(name="yo5_qd9", len=1.11, k1=-0.021207658757208966)
drift_1173 = DRIFT(name="drift_1173", len=0.5452500000001237)
yo5_tv9 = CORRECTOR(name="yo5_tv9", len=0.0, xkick=0.0, ykick=-0.0)
yo5_oct9 = thinMULTIPOLE(name="yo5_oct9", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yo5_dec9 = thinMULTIPOLE(name="yo5_dec9", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yo5_qs9 = thinMULTIPOLE(name="yo5_qs9", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_1174 = DRIFT(name="drift_1174", len=1.593705149733978)
yo5_dh8 = ESBEND(name="yo5_dh8", len=9.440656, angle=0.038924026765774174, e1=0.0, e2=0.0)
drift_1175 = DRIFT(name="drift_1175", len=1.593705149733978)
yo5_th8 = CORRECTOR(name="yo5_th8", len=0.0, xkick=-0.0, ykick=0.0)
yo5_oct8 = thinMULTIPOLE(name="yo5_oct8", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yo5_dec8 = thinMULTIPOLE(name="yo5_dec8", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yo5_qgt8 = thinMULTIPOLE(name="yo5_qgt8", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, -0.0, 0.0, 0.0])
drift_1176 = DRIFT(name="drift_1176", len=0.5452500000001237)
yo5_qf8 = KQUAD(name="yo5_qf8", len=1.11, k1=-0.04018877120650051)
drift_1177 = DRIFT(name="drift_1177", len=0.29610560000037367)
yo5_b8 = MARKER(name="yo5_b8")
drift_1178 = DRIFT(name="drift_1178", len=1.3184179478494116)
yo5_snk_hlx4 = DRIFT(name="yo5_snk_hlx4", len=2.4)
drift_1179 = DRIFT(name="drift_1179", len=0.21200000000044383)
yo5_snk_hlx3 = DRIFT(name="yo5_snk_hlx3", len=2.4)
drift_1180 = DRIFT(name="drift_1180", len=0.22400000000016007)
yo5_bsnk = MARKER(name="yo5_bsnk")
drift_1181 = DRIFT(name="drift_1181", len=0.22400000000016007)
yo5_snk_hlx2 = DRIFT(name="yo5_snk_hlx2", len=2.4)
drift_1182 = DRIFT(name="drift_1182", len=0.21200000000044383)
yo5_snk_hlx1 = DRIFT(name="yo5_snk_hlx1", len=2.4)
drift_1183 = DRIFT(name="drift_1183", len=1.3184179478494116)
yo5_b7 = MARKER(name="yo5_b7")
drift_1184 = DRIFT(name="drift_1184", len=0.2962336000000505)
yo5_qd7 = KQUAD(name="yo5_qd7", len=0.929744, k1=0.06885631651853527)
drift_1185 = DRIFT(name="drift_1185", len=0.5453779999998005)
yo5_tv7 = CORRECTOR(name="yo5_tv7", len=0.0, xkick=0.0, ykick=-0.0)
yo5_oct7 = thinMULTIPOLE(name="yo5_oct7", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yo5_dec7 = thinMULTIPOLE(name="yo5_dec7", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yo5_qs7 = thinMULTIPOLE(name="yo5_qs7", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_1186 = DRIFT(name="drift_1186", len=1.8803175478496996)
yo5_rot_hlx4 = DRIFT(name="yo5_rot_hlx4", len=2.4)
drift_1187 = DRIFT(name="drift_1187", len=0.21200000000044383)
yo5_rot_hlx3 = DRIFT(name="yo5_rot_hlx3", len=2.4)
drift_1188 = DRIFT(name="drift_1188", len=0.22400000000016007)
yo5_brot = MARKER(name="yo5_brot")
drift_1189 = DRIFT(name="drift_1189", len=0.22400000000016007)
yo5_rot_hlx2 = DRIFT(name="yo5_rot_hlx2", len=2.4)
drift_1190 = DRIFT(name="drift_1190", len=0.21200000000044383)
yo5_rot_hlx1 = DRIFT(name="yo5_rot_hlx1", len=2.4)
drift_1191 = DRIFT(name="drift_1191", len=1.8803175478496996)
bi5_tv8 = CORRECTOR(name="bi5_tv8", len=0.0, xkick=0.0, ykick=-0.0)
bi5_oct8 = thinMULTIPOLE(name="bi5_oct8", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
bi5_dec8 = thinMULTIPOLE(name="bi5_dec8", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
bi5_qs8 = thinMULTIPOLE(name="bi5_qs8", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_1192 = DRIFT(name="drift_1192", len=0.5452500000001237)
bi5_qd8 = KQUAD(name="bi5_qd8", len=1.11, k1=-0.07448144564520073)
drift_1193 = DRIFT(name="drift_1193", len=0.29610560000037367)
bi5_b8 = MARKER(name="bi5_b8")
drift_1194 = DRIFT(name="drift_1194", len=1.014665032902485)
h5_dh4 = ESBEND(name="h5_dh4", len=6.915999880527922, angle=0.028514815544784158, e1=0.0, e2=0.0)
drift_1195 = DRIFT(name="drift_1195", len=1.5765646329027732)
yo5_th4 = CORRECTOR(name="yo5_th4", len=0.0, xkick=-0.0, ykick=0.0)
yo5_oct4 = thinMULTIPOLE(name="yo5_oct4", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yo5_dec4 = thinMULTIPOLE(name="yo5_dec4", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yo5_qs4 = thinMULTIPOLE(name="yo5_qs4", PolynomA=[0.0, -0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_1196 = DRIFT(name="drift_1196", len=0.5442755000003672)
yo5_qf4 = KQUAD(name="yo5_qf4", len=1.811949, k1=0.09200921857209204)
drift_1197 = DRIFT(name="drift_1197", len=0.13007550000020274)
yo5_tq4 = KQUAD(name="yo5_tq4", len=0.75, k1=0.051933404287653466)
drift_1198 = DRIFT(name="drift_1198", len=0.22623959999964427)
yo5_b4 = MARKER(name="yo5_b4")
drift_1199 = DRIFT(name="drift_1199", len=1.406339555202976)
yo5_dh5 = ESBEND(name="yo5_dh5", len=8.698449375192075, angle=0.035863892964716426, e1=0.0, e2=0.0)
drift_1200 = DRIFT(name="drift_1200", len=1.968239155203264)
bo3_th8 = CORRECTOR(name="bo3_th8", len=0.0, xkick=-0.0, ykick=0.0)
bo3_oct8 = thinMULTIPOLE(name="bo3_oct8", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
bo3_dec8 = thinMULTIPOLE(name="bo3_dec8", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
bo3_qgt8 = thinMULTIPOLE(name="bo3_qgt8", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, -0.0, 0.0, 0.0])
drift_1201 = DRIFT(name="drift_1201", len=0.5452500000001237)
bo3_qf8 = KQUAD(name="bo3_qf8", len=1.11, k1=-0.09765754532799617)
drift_1202 = DRIFT(name="drift_1202", len=0.29610560000037367)
bo3_b8 = MARKER(name="bo3_b8")
drift_1203 = DRIFT(name="drift_1203", len=2.1340989004656876)
yi6_tv2 = CORRECTOR(name="yi6_tv2", len=0.0, xkick=0.0, ykick=-0.0)
yi6_oct2 = thinMULTIPOLE(name="yi6_oct2", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, -0.0])
yi6_dec2 = thinMULTIPOLE(name="yi6_dec2", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
yi6_dod2 = thinMULTIPOLE(name="yi6_dod2", PolynomA=[0.0, 0.0, 0.0, 0.0], PolynomB=[0.0, 0.0, 0.0, 0.0])
drift_1204 = DRIFT(name="drift_1204", len=0.49420350000036706)
yi6_qd2 = KQUAD(name="yi6_qd2", len=3.391633, k1=0.03256333415663846)
drift_1205 = DRIFT(name="drift_1205", len=1.1740806999996494)
# o_crab_ip6f = DRIFT(name="o_crab_ip6f", len=15.06)
o_crab_ip6f_d1 = DRIFT(name="o_crab_ip6f_d1", len=(15.06-4.0)/2)
o_crab_ip6f = CRABCAVITY(name="o_crab_ip6f", len=4.0, volt=0.0, freq=197e6, energy=275e6)
o_crab_ip6f_d2 = DRIFT(name="o_crab_ip6f_d2", len=(15.06-4.0)/2)
drift_1206 = DRIFT(name="drift_1206", len=1.0)
b2pf = ESBEND(name="b2pf", len=3.0774350498350653, angle=0.016774494326212253, e1=0.008387247163106126, e2=0.008387247163106126)
drift_1207 = DRIFT(name="drift_1207", len=1.0)
h5_tv3 = CORRECTOR(name="h5_tv3", len=0.25, xkick=0.0, ykick=0.0)
drift_1208 = DRIFT(name="drift_1208", len=0.1499999999996362)
h5_qs3 = KQUAD(name="h5_qs3", len=0.25, k1=0.0)
drift_1209 = DRIFT(name="drift_1209", len=0.1499999999996362)
q3pf = KQUAD(name="q3pf", len=0.75, k1=0.010493191517541332)
drift_1210 = DRIFT(name="drift_1210", len=9.896191084143538)
o_roman_pot_ip6_2_1 = DRIFT(name="o_roman_pot_ip6_2_1", len=0.2)
o_roman_pot_ip6_2 = MARKER(name="o_roman_pot_ip6_2")
drift_1211 = DRIFT(name="drift_1211", len=1.6020090158435778)
o_roman_pot_ip6_1_1 = DRIFT(name="o_roman_pot_ip6_1_1", len=0.2)
o_roman_pot_ip6_1 = MARKER(name="o_roman_pot_ip6_1")
drift_1212 = DRIFT(name="drift_1212", len=1.1015067618827743)
o_off_mom_ip6_2_1 = DRIFT(name="o_off_mom_ip6_2_1", len=0.2)
o_off_mom_ip6_2 = MARKER(name="o_off_mom_ip6_2")
drift_1213 = DRIFT(name="drift_1213", len=1.6020090158435778)
o_off_mom_ip6_1_1 = DRIFT(name="o_off_mom_ip6_1_1", len=0.2)
o_off_mom_ip6_1 = MARKER(name="o_off_mom_ip6_1")
drift_1214 = DRIFT(name="drift_1214", len=0.23092092654042062)
pbr_b1apf = YROTATION(name="pbr_b1apf",angle= 2.38032985038345683e-02)
pbt_b1apf = TRANSLATION(name="pbt_b1apf",dx=-1.53882030013441119e-02)
b1apf = ESBEND(name="b1apf", len=1.5, angle=0.0, e1=0.0, e2=0.0, PolynomB=[-2.94343399274668354e-03-0.0/1.5, 0.0, 0.0, 0.0])
pet_b1apf = TRANSLATION(name="pet_b1apf",dx=-1.70096429725014632e-02)
per_b1apf = YROTATION(name="per_b1apf",angle= -1.93871142166689610e-02)
drift_1215 = DRIFT(name="drift_1215", len=0.5016850837464517)
pbr_b1pf = YROTATION(name="pbr_b1pf",angle= 1.25950476108859111e-02)
pbt_b1pf = TRANSLATION(name="pbt_b1pf",dx=-6.26799876724247647e-03)
b1pf = ESBEND(name="b1pf", len=3.0, angle=0.0, e1=0.0, e2=0.0, PolynomB=[-3.70706499213125252e-03-0.0/3.0, 0.0, 0.0, 0.0])
pet_b1pf = TRANSLATION(name="pet_b1pf",dx=-1.48375342617892581e-02)
per_b1pf = YROTATION(name="per_b1pf",angle= -1.47507561726256458e-03)
drift_1216 = DRIFT(name="drift_1216", len=0.5015444596829184)
pbr_q2pf = YROTATION(name="pbr_q2pf",angle= 1.54668927288380952e-02)
pbt_q2pf = TRANSLATION(name="pbt_q2pf",dx=-1.13897121786721012e-02)
q2pf = KQUAD(name="q2pf", len=3.8, k1=0.037643704139219386)
pet_q2pf = TRANSLATION(name="pet_q2pf",dx=-4.51593211840057920e-02)
per_q2pf = YROTATION(name="per_q2pf",angle= -1.29367560655846096e-02)
drift_1217 = DRIFT(name="drift_1217", len=0.4010213247456704)
pbr_q1bpf = YROTATION(name="pbr_q1bpf",angle= 1.31828003697269328e-02)
pbt_q1bpf = TRANSLATION(name="pbt_q1bpf",dx=4.55987368366134183e-06)
q1bpf = KQUAD(name="q1bpf", len=1.61, k1=-0.05252692321374404)
pet_q1bpf = TRANSLATION(name="pet_q1bpf",dx= -2.17131032452377963e-02)
per_q1bpf = YROTATION(name="per_q1bpf",angle= -1.40902338121864472e-02)
drift_1218 = DRIFT(name="drift_1218", len=0.40059577412102954)
pbr_q1apf = YROTATION(name="pbr_q1apf",angle= 6.76920671694013470e-03)
pbt_q1apf = TRANSLATION(name="pbt_q1apf",dx= 4.05191841125844728e-03)
q1apf = KQUAD(name="q1apf", len=1.46, k1=-0.09073171512043103)
pet_q1apf = TRANSLATION(name="pet_q1apf",dx= -1.46541743061257106e-02)
per_q1apf = YROTATION(name="per_q1apf",angle= -7.98771680129673443e-03)
drift_1219 = DRIFT(name="drift_1219", len=0.4006259721118113)
pbr_b0apf = YROTATION(name="pbr_b0apf",angle= 6.95370413729287911e-03)
pbt_b0apf = TRANSLATION(name="pbt_b0apf",dx= 3.55265096116321247e-03)
b0apf = ESBEND(name="b0apf", len=0.6, angle=0.0, e1=0.0, e2=0.0, PolynomB=[-3.59402012491794346e-03-0.0/0.6, 0.0, 0.0, 0.0])
pet_b0apf = TRANSLATION(name="pet_b0apf",dx= -7.07734721058384929e-03)
per_b0apf = YROTATION(name="per_b0apf",angle= -4.79514821232049273e-03)
drift_1220 = DRIFT(name="drift_1220", len=0.9013554989205659)
pbt_b0pf = TRANSLATION(name="pbt_b0pf",dx= -3.02264219559301335e-02)
pbr_b0pf = YROTATION(name="pbr_b0pf",angle= 2.67043931616361939e-02)
b0pf = ESBEND(name="b0pf", len=1.2, angle=0.0, e1=0.0, e2=0.0, PolynomB=[1.45779264977192957e-01*-0.008854545281644094-0.0/1.2, -0.008854545281644094, 0.0, 0.0])
pet_b0pf = TRANSLATION(name="pet_b0pf",dx= -7.49049089865311268e-04)
per_b0pf = YROTATION(name="per_b0pf",angle= -2.50000000000000014e-02)
drift_1221 = DRIFT(name="drift_1221", len=3.8024379395737924)
pbt_sol_f = TRANSLATION(name="pbt_sol_f",dx= -4.99947918294246646e-02)
pbr_sol_f = YROTATION(name="pbr_sol_f",angle= 2.50000000000000014e-02)
star_detect_f = SOLENOID(name="star_detect_f", len=2.0)
per_sol_f = YROTATION(name="per_sol_f",angle= -2.50000000000000014e-02)
ip6w = MARKER(name="ip6w")
ring = [
    ip6d,
    pbr_sol_r,
    star_detect_r,
    pet_sol_r,
    per_sol_r,
    drift_0,
    q1apr,
    drift_1,
    q1bpr,
    drift_2,
    q2pr,
    drift_3,
    q3pr,
    drift_4,
    o_crab_ip6r_d1,
    o_crab_ip6r,
    o_crab_ip6r_d2,
    drift_5,
    q4pr,
    drift_6,
    b1pr,
    drift_7,
    q5pr,
    drift_8,
    yi6_b4,
    drift_9,
    yi6_tq4,
    drift_10,
    yi6_qd4,
    drift_11,
    yi6_tv4,
    yi6_oct4,
    yi6_dec4,
    yi6_qs4,
    drift_12,
    o_rotator_ip6r,
    drift_13,
    yi6_bh5,
    drift_14,
    yi6_tq5,
    drift_15,
    yi6_qf5,
    drift_16,
    yi6_th5,
    yi6_oct5,
    yi6_dec5,
    yi6_qgt5,
    drift_17,
    b2pr,
    drift_18,
    yi6_bv6,
    drift_19,
    yi6_tq6,
    drift_20,
    yi6_qd6,
    drift_21,
    yi6_tv6,
    yi6_oct6,
    yi6_dec6,
    yi6_qs6,
    drift_22,
    yi6_dh5,
    drift_23,
    yi6_th7,
    yi6_oct7,
    yi6_dec7,
    yi6_qgt7,
    drift_24,
    yi6_qf7,
    drift_25,
    yi6_b7,
    drift_26,
    yi6_dh6,
    drift_27,
    yi6_b8,
    drift_28,
    yi6_qd8,
    drift_29,
    yi6_tv8,
    yi6_oct8,
    yi6_dec8,
    yi6_qs8,
    drift_30,
    yi6_dh9,
    drift_31,
    yi6_th9,
    yi6_oct9,
    yi6_dec9,
    yi6_qs9,
    drift_32,
    yi6_qf9,
    drift_33,
    yi6_bh9,
    drift_34,
    yi6_dh8,
    drift_35,
    yi6_bv10,
    drift_36,
    yi6_sxd10,
    drift_37,
    yi6_qd10,
    drift_38,
    yi6_tv10,
    yi6_oct10,
    yi6_dec10,
    yi6_qs10,
    drift_39,
    yi6_dh10,
    drift_40,
    yi6_bh11,
    drift_41,
    yi6_sxf11,
    drift_42,
    yi6_qf11,
    drift_43,
    yi6_th11,
    drift_44,
    yi6_dh11,
    drift_45,
    yi6_bv12,
    drift_46,
    yi6_sxd12,
    drift_47,
    yi6_qd12,
    drift_48,
    yi6_tv12,
    drift_49,
    yi6_dh12,
    drift_50,
    yi6_bh13,
    drift_51,
    yi6_sxf13,
    drift_52,
    yi6_qf13,
    drift_53,
    yi6_th13,
    drift_54,
    yi6_dh13,
    drift_55,
    yi6_bv14,
    drift_56,
    yi6_sxd14,
    drift_57,
    yi6_qd14,
    drift_58,
    yi6_tv14,
    drift_59,
    yi6_dh14,
    drift_60,
    yi6_bh15,
    drift_61,
    yi6_sxf15,
    drift_62,
    yi6_qf15,
    drift_63,
    yi6_th15,
    drift_64,
    yi6_dh15,
    drift_65,
    yi6_bv16,
    drift_66,
    yi6_sxd16,
    drift_67,
    yi6_qd16,
    drift_68,
    yi6_tv16,
    drift_69,
    yi6_dh16,
    drift_70,
    yi6_bh17,
    drift_71,
    yi6_sxf17,
    drift_72,
    yi6_qf17,
    drift_73,
    yi6_th17,
    drift_74,
    yi6_dh17,
    drift_75,
    yi6_bv18,
    drift_76,
    yi6_sxd18,
    drift_77,
    yi6_qd18,
    drift_78,
    yi6_tv18,
    drift_79,
    yi6_dh18,
    drift_80,
    yi6_bh19,
    drift_81,
    yi6_sxf19,
    drift_82,
    yi6_qf19,
    drift_83,
    yi6_th19,
    drift_84,
    yi6_dh19,
    drift_85,
    yi6_bv20,
    drift_86,
    yi6_sxd20,
    drift_87,
    yi6_qd20,
    drift_88,
    yi6_tv20,
    drift_89,
    yi6_dh20,
    drift_90,
    yi7_bh21,
    drift_91,
    yi7_sxf21,
    drift_92,
    yi7_qf21,
    drift_93,
    yi7_th21,
    drift_94,
    yi7_dh20,
    drift_95,
    yi7_bv20,
    drift_96,
    yi7_sxd20,
    drift_97,
    yi7_qd20,
    drift_98,
    yi7_tv20,
    drift_99,
    yi7_dh19,
    drift_100,
    yi7_bh19,
    drift_101,
    yi7_sxf19,
    drift_102,
    yi7_qf19,
    drift_103,
    yi7_th19,
    drift_104,
    yi7_dh18,
    drift_105,
    yi7_bv18,
    drift_106,
    yi7_sxd18,
    drift_107,
    yi7_qd18,
    drift_108,
    yi7_tv18,
    yi7_oct18,
    yi7_dec18,
    yi7_qs18,
    drift_109,
    yi7_dh17,
    drift_110,
    yi7_bh17,
    drift_111,
    yi7_sxf17,
    drift_112,
    yi7_qf17,
    drift_113,
    yi7_th17,
    yi7_oct17,
    yi7_dec17,
    yi7_qgt17,
    drift_114,
    yi7_dh16,
    drift_115,
    yi7_bv16,
    drift_116,
    yi7_sxd16,
    drift_117,
    yi7_qd16,
    drift_118,
    yi7_tv16,
    yi7_oct16,
    yi7_dec16,
    yi7_qs16,
    drift_119,
    yi7_dh15,
    drift_120,
    yi7_bh15,
    drift_121,
    yi7_sxf15,
    drift_122,
    yi7_qf15,
    drift_123,
    yi7_th15,
    yi7_oct15,
    yi7_dec15,
    yi7_qgt15,
    drift_124,
    yi7_dh14,
    drift_125,
    yi7_bv14,
    drift_126,
    yi7_sxd14,
    drift_127,
    yi7_qd14,
    drift_128,
    yi7_tv14,
    yi7_oct14,
    yi7_dec14,
    yi7_qs14,
    drift_129,
    yi7_dh13,
    drift_130,
    yi7_bh13,
    drift_131,
    yi7_sxf13,
    drift_132,
    yi7_qf13,
    drift_133,
    yi7_th13,
    yi7_oct13,
    yi7_dec13,
    yi7_qgt13,
    drift_134,
    yi7_dh12,
    drift_135,
    yi7_bv12,
    drift_136,
    yi7_sxd12,
    drift_137,
    yi7_qd12,
    drift_138,
    yi7_tv12,
    yi7_oct12,
    yi7_dec12,
    yi7_qs12,
    drift_139,
    yi7_dh11,
    drift_140,
    yi7_bh11,
    drift_141,
    yi7_sxf11,
    drift_142,
    yi7_qf11,
    drift_143,
    yi7_th11,
    yi7_oct11,
    yi7_dec11,
    yi7_qgt11,
    drift_144,
    yi7_dh10,
    drift_145,
    yi7_bv10,
    drift_146,
    yi7_sxd10,
    drift_147,
    yi7_qd10,
    drift_148,
    yi7_tv10,
    yi7_oct10,
    yi7_dec10,
    yi7_qs10,
    drift_149,
    yi7_dh9,
    drift_150,
    qds10,
    drift_151,
    qds09,
    drift_152,
    bxds9m2,
    drift_153,
    qds08,
    drift_154,
    bxdsds04,
    drift_155,
    qds07,
    drift_156,
    bxds9m1,
    drift_157,
    qds06,
    drift_158,
    yi7_rot_hlx4,
    drift_159,
    yi7_rot_hlx3,
    drift_160,
    yi7_brot,
    drift_161,
    yi7_rot_hlx2,
    drift_162,
    yi7_rot_hlx1,
    drift_163,
    qds05,
    drift_164,
    qds04b,
    drift_165,
    qds04a,
    o_crab_ir8w_d1,
    o_crab_ir8w,
    o_crab_ir8w_d2,
    qds03b,
    drift_166,
    qds03a,
    drift_167,
    bxds02disp1,
    drift_168,
    qds02a,
    drift_169,
    qds02,
    drift_170,
    mpot1o,
    drift_171,
    mpot1,
    drift_172,
    mpot1i,
    drift_173,
    qds01,
    drift_174,
    bxds01b,
    drift_175,
    pwt_bxds01a,
    pwr_bxds01a,
    bxds01a,
    pdr_bxds01a,
    pdt_bxds01a,
    drift_176,
    pwt_qffds02b,
    pwr_qffds02b,
    qffds02b,
    pdr_qffds02b,
    pdt_qffds02b,
    drift_177,
    pwt_qffds02a,
    pwr_qffds02a,
    qffds02a,
    pdr_qffds02a,
    pdt_qffds02a,
    drift_178,
    pwt_qffds01b,
    pwr_qffds01b,
    qffds01b,
    pdr_qffds01b,
    pdt_qffds01b,
    drift_179,
    pwt_qffds01a,
    pwr_qffds01a,
    qffds01a,
    pdr_qffds01a,
    pdt_qffds01a,
    drift_180,
    pwt_bxsp01,
    pwr_bxsp01,
    bxsp01,
    pdr_bxsp01,
    pdt_bxsp01,
    drift_181,
    solds,
    ip8,
    solus,
    drift_182,
    qffus01,
    drift_183,
    qffus02,
    drift_184,
    qffus03,
    drift_185,
    qus01c,
    drift_186,
    o_crab_ir8d_d1,
    o_crab_ir8d,
    o_crab_ir8d_d2,
    drift_187,
    qus02c,
    drift_188,
    qus03,
    drift_189,
    bxus9m1,
    drift_190,
    qus04,
    drift_191,
    yo8_rot_hlx4,
    drift_192,
    yo8_rot_hlx3,
    drift_193,
    yo8_brot,
    drift_194,
    yo8_rot_hlx2,
    drift_195,
    yo8_rot_hlx1,
    drift_196,
    qus05,
    drift_197,
    bxus9m2,
    drift_198,
    qus06,
    drift_199,
    yo8_snk_hlx4,
    drift_200,
    yo8_snk_hlx3,
    drift_201,
    yo8_bsnk,
    drift_202,
    yo8_snk_hlx2,
    drift_203,
    yo8_snk_hlx1,
    drift_204,
    qus07,
    drift_205,
    bxus9m3,
    drift_206,
    qus08,
    drift_207,
    qus09,
    drift_208,
    yo8_dh9,
    drift_209,
    yo8_th10,
    yo8_oct10,
    yo8_dec10,
    yo8_qs10,
    drift_210,
    yo8_qf10,
    drift_211,
    yo8_sxf10,
    drift_212,
    yo8_bh10,
    drift_213,
    yo8_dh10,
    drift_214,
    yo8_tv11,
    drift_215,
    yo8_qd11,
    drift_216,
    yo8_sxd11,
    drift_217,
    yo8_bv11,
    drift_218,
    yo8_dh11,
    drift_219,
    yo8_th12,
    drift_220,
    yo8_qf12,
    drift_221,
    yo8_sxf12,
    drift_222,
    yo8_bh12,
    drift_223,
    yo8_dh12,
    drift_224,
    yo8_tv13,
    drift_225,
    yo8_qd13,
    drift_226,
    yo8_sxd13,
    drift_227,
    yo8_bv13,
    drift_228,
    yo8_dh13,
    drift_229,
    yo8_th14,
    drift_230,
    yo8_qf14,
    drift_231,
    yo8_sxf14,
    drift_232,
    yo8_bh14,
    drift_233,
    yo8_dh14,
    drift_234,
    yo8_tv15,
    drift_235,
    yo8_qd15,
    drift_236,
    yo8_sxd15,
    drift_237,
    yo8_bv15,
    drift_238,
    yo8_dh15,
    drift_239,
    yo8_th16,
    drift_240,
    yo8_qf16,
    drift_241,
    yo8_sxf16,
    drift_242,
    yo8_bh16,
    drift_243,
    yo8_dh16,
    drift_244,
    yo8_tv17,
    drift_245,
    yo8_qd17,
    drift_246,
    yo8_sxd17,
    drift_247,
    yo8_bv17,
    drift_248,
    yo8_dh17,
    drift_249,
    yo8_th18,
    drift_250,
    yo8_qf18,
    drift_251,
    yo8_sxf18,
    drift_252,
    yo8_bh18,
    drift_253,
    yo8_dh18,
    drift_254,
    yo8_tv19,
    drift_255,
    yo8_qd19,
    drift_256,
    yo8_sxd19,
    drift_257,
    yo8_bv19,
    drift_258,
    yo8_dh19,
    drift_259,
    yo8_th20,
    drift_260,
    yo8_qf20,
    drift_261,
    yo8_sxf20,
    drift_262,
    yo8_bh20,
    drift_263,
    yo8_dh20,
    drift_264,
    yo9_tv21,
    drift_265,
    yo9_qd21,
    drift_266,
    yo9_sxd21,
    drift_267,
    yo9_bv21,
    drift_268,
    yo9_dh20,
    drift_269,
    yo9_th20,
    drift_270,
    yo9_qf20,
    drift_271,
    yo9_sxf20,
    drift_272,
    yo9_bh20,
    drift_273,
    yo9_dh19,
    drift_274,
    yo9_tv19,
    drift_275,
    yo9_qd19,
    drift_276,
    yo9_sxd19,
    drift_277,
    yo9_bv19,
    drift_278,
    yo9_dh18,
    drift_279,
    yo9_th18,
    yo9_oct18,
    yo9_dec18,
    yo9_qgt18,
    drift_280,
    yo9_qf18,
    drift_281,
    yo9_sxf18,
    drift_282,
    yo9_bh18,
    drift_283,
    yo9_dh17,
    drift_284,
    yo9_tv17,
    yo9_oct17,
    yo9_dec17,
    yo9_qs17,
    drift_285,
    yo9_qd17,
    drift_286,
    yo9_sxd17,
    drift_287,
    yo9_bv17,
    drift_288,
    yo9_dh16,
    drift_289,
    yo9_th16,
    yo9_oct16,
    yo9_dec16,
    yo9_qgt16,
    drift_290,
    yo9_qf16,
    drift_291,
    yo9_sxf16,
    drift_292,
    yo9_bh16,
    drift_293,
    yo9_dh15,
    drift_294,
    yo9_tv15,
    yo9_oct15,
    yo9_dec15,
    yo9_qs15,
    drift_295,
    yo9_qd15,
    drift_296,
    yo9_sxd15,
    drift_297,
    yo9_bv15,
    drift_298,
    yo9_dh14,
    drift_299,
    yo9_th14,
    yo9_oct14,
    yo9_dec14,
    yo9_qgt14,
    drift_300,
    yo9_qf14,
    drift_301,
    yo9_sxf14,
    drift_302,
    yo9_bh14,
    drift_303,
    yo9_dh13,
    drift_304,
    yo9_tv13,
    yo9_oct13,
    yo9_dec13,
    yo9_qs13,
    drift_305,
    yo9_qd13,
    drift_306,
    yo9_sxd13,
    drift_307,
    yo9_bv13,
    drift_308,
    yo9_dh12,
    drift_309,
    yo9_th12,
    yo9_oct12,
    yo9_dec12,
    yo9_qgt12,
    drift_310,
    yo9_qf12,
    drift_311,
    yo9_sxf12,
    drift_312,
    yo9_bh12,
    drift_313,
    yo9_dh11,
    drift_314,
    yo9_tv11,
    yo9_oct11,
    yo9_dec11,
    yo9_qs11,
    drift_315,
    yo9_qd11,
    drift_316,
    yo9_sxd11,
    drift_317,
    yo9_bv11,
    drift_318,
    yo9_dh10,
    drift_319,
    yo9_th10,
    yo9_oct10,
    yo9_dec10,
    yo9_qs10,
    drift_320,
    yo9_qf10,
    drift_321,
    yo9_sxf10,
    drift_322,
    yo9_bh10,
    drift_323,
    yo9_dh9,
    drift_324,
    yo9_bv9,
    drift_325,
    yo9_qd9,
    drift_326,
    yo9_tv9,
    yo9_oct9,
    yo9_dec9,
    yo9_qs9,
    drift_327,
    yo9_dh8,
    drift_328,
    yo9_th8,
    yo9_oct8,
    yo9_dec8,
    yo9_qgt8,
    drift_329,
    yo9_qf8,
    drift_330,
    yo9_b8,
    drift_331,
    yo9_hlx7_4,
    drift_332,
    yo9_hlx7_3,
    drift_333,
    yo9_b7_1,
    drift_334,
    yo9_hlx7_2,
    drift_335,
    yo9_hlx7_1,
    drift_336,
    yo9_b7,
    drift_337,
    yo9_qd7,
    drift_338,
    yo9_tv7,
    yo9_oct7,
    yo9_dec7,
    yo9_qs7,
    drift_339,
    yo9_dh6,
    drift_340,
    yo9_th6,
    yo9_oct6,
    yo9_dec6,
    yo9_qgt6,
    drift_341,
    yo9_qf6,
    drift_342,
    yo9_tq6,
    drift_343,
    yo9_bh6,
    drift_344,
    yo9_dh5,
    drift_345,
    yo9_tv5,
    yo9_oct5,
    yo9_dec5,
    yo9_qs5,
    drift_346,
    yo9_qd5,
    drift_347,
    yo9_tq5,
    drift_348,
    yo9_bv5,
    drift_349,
    yo9_th4,
    yo9_oct4,
    yo9_dec4,
    yo9_qs4,
    drift_350,
    yo9_qf4,
    drift_351,
    yo9_tq4,
    drift_352,
    yo9_b4,
    drift_353,
    yo9_sv4,
    drift_354,
    yo9_dmp3_2,
    yo9_dmp3_1,
    drift_355,
    yo9_sv3_2,
    drift_356,
    yo9_b3_1,
    drift_357,
    yo9_kfbh3,
    drift_358,
    yo9_ka3_5,
    drift_359,
    yo9_ka3_4,
    drift_360,
    yo9_ka3_3,
    drift_361,
    yo9_ka3_2,
    drift_362,
    yo9_ka3_1,
    drift_363,
    yo9_sv3_1,
    drift_364,
    yo9_b3,
    drift_365,
    yo9_tv3,
    yo9_sx3,
    yo9_oct3,
    yo9_dod3,
    drift_366,
    yo9_qd3,
    drift_367,
    yo9_dods3,
    yo9_octs3,
    yo9_sxs3,
    yo9_qs3,
    drift_368,
    yo9_qf2,
    drift_369,
    yo9_th2,
    yo9_oct2,
    yo9_dec2,
    yo9_dod2,
    drift_370,
    yo9_qd1,
    drift_371,
    yo9_b1,
    drift_372,
    cavity_591mhz,
    drift_373,
    dsw_ir10h,
    drift_374,
    dwarm_ir10h,
    drift_375,
    bo10_b1,
    drift_376,
    bo10_qd1,
    drift_377,
    bo10_th2,
    bo10_oct2,
    bo10_dec2,
    bo10_dod2,
    drift_378,
    bo10_qf2,
    drift_379,
    bo10_dods3,
    bo10_octs3,
    bo10_sxs3,
    bo10_qs3,
    drift_380,
    bo10_qd3,
    drift_381,
    bo10_tv3,
    bo10_sx3,
    bo10_oct3,
    bo10_dod3,
    drift_382,
    bo10_b3,
    drift_383,
    bo10_sv3_1,
    drift_384,
    bo10_ka3_1,
    drift_385,
    bo10_ka3_2,
    drift_386,
    bo10_ka3_3,
    drift_387,
    bo10_ka3_4,
    drift_388,
    bo10_ka3_5,
    drift_389,
    bo10_kfbh3,
    drift_390,
    bo10_b3_1,
    drift_391,
    bo10_sv3_2,
    drift_392,
    bo10_c3,
    drift_393,
    bo10_dmp3_1,
    bo10_dmp3_2,
    drift_394,
    bo10_sv4,
    drift_395,
    bo10_b4,
    drift_396,
    bo10_tq4,
    drift_397,
    bo10_qf4,
    drift_398,
    bo10_th4,
    bo10_oct4,
    bo10_dec4,
    bo10_qs4,
    drift_399,
    bo10_bv5,
    drift_400,
    bo10_tq5,
    drift_401,
    bo10_qd5,
    drift_402,
    bo10_tv5,
    bo10_oct5,
    bo10_dec5,
    bo10_qs5,
    drift_403,
    bo10_dh5,
    drift_404,
    bo10_bh6,
    drift_405,
    bo10_tq6,
    drift_406,
    bo10_qf6,
    drift_407,
    bo10_th6,
    bo10_oct6,
    bo10_dec6,
    bo10_qgt6,
    drift_408,
    bo10_dh6,
    drift_409,
    bo10_tv7,
    bo10_oct7,
    bo10_dec7,
    bo10_qs7,
    drift_410,
    bo10_qd7,
    drift_411,
    bo10_b7,
    drift_412,
    bo10_b8,
    drift_413,
    bo10_qf8,
    drift_414,
    bo10_th8,
    bo10_oct8,
    bo10_dec8,
    bo10_qgt8,
    drift_415,
    bo10_dh8,
    drift_416,
    bo10_tv9,
    bo10_oct9,
    bo10_dec9,
    bo10_qs9,
    drift_417,
    bo10_qd9,
    drift_418,
    bo10_bv9,
    drift_419,
    bo10_dh9,
    drift_420,
    bo10_bh10,
    drift_421,
    bo10_sxf10,
    drift_422,
    bo10_qf10,
    drift_423,
    bo10_th10,
    bo10_oct10,
    bo10_dec10,
    bo10_qs10,
    drift_424,
    bo10_dh10,
    drift_425,
    bo10_bv11,
    drift_426,
    bo10_sxd11,
    drift_427,
    bo10_qd11,
    drift_428,
    bo10_tv11,
    bo10_oct11,
    bo10_dec11,
    bo10_qs11,
    drift_429,
    bo10_dh11,
    drift_430,
    bo10_bh12,
    drift_431,
    bo10_sxf12,
    drift_432,
    bo10_qf12,
    drift_433,
    bo10_th12,
    bo10_oct12,
    bo10_dec12,
    bo10_qgt12,
    drift_434,
    bo10_dh12,
    drift_435,
    bo10_bv13,
    drift_436,
    bo10_sxd13,
    drift_437,
    bo10_qd13,
    drift_438,
    bo10_tv13,
    bo10_oct13,
    bo10_dec13,
    bo10_qs13,
    drift_439,
    bo10_dh13,
    drift_440,
    bo10_bh14,
    drift_441,
    bo10_sxf14,
    drift_442,
    bo10_qf14,
    drift_443,
    bo10_th14,
    bo10_oct14,
    bo10_dec14,
    bo10_qgt14,
    drift_444,
    bo10_dh14,
    drift_445,
    bo10_bv15,
    drift_446,
    bo10_sxd15,
    drift_447,
    bo10_qd15,
    drift_448,
    bo10_tv15,
    bo10_oct15,
    bo10_dec15,
    bo10_qs15,
    drift_449,
    bo10_dh15,
    drift_450,
    bo10_bh16,
    drift_451,
    bo10_sxf16,
    drift_452,
    bo10_qf16,
    drift_453,
    bo10_th16,
    bo10_oct16,
    bo10_dec16,
    bo10_qgt16,
    drift_454,
    bo10_dh16,
    drift_455,
    bo10_bv17,
    drift_456,
    bo10_sxd17,
    drift_457,
    bo10_qd17,
    drift_458,
    bo10_tv17,
    bo10_oct17,
    bo10_dec17,
    bo10_qs17,
    drift_459,
    bo10_dh17,
    drift_460,
    bo10_bh18,
    drift_461,
    bo10_sxf18,
    drift_462,
    bo10_qf18,
    drift_463,
    bo10_th18,
    bo10_oct18,
    bo10_dec18,
    bo10_qgt18,
    drift_464,
    bo10_dh18,
    drift_465,
    bo10_bv19,
    drift_466,
    bo10_sxd19,
    drift_467,
    bo10_qd19,
    drift_468,
    bo10_tv19,
    drift_469,
    bo10_dh19,
    drift_470,
    bo10_bh20,
    drift_471,
    bo10_sxf20,
    drift_472,
    bo10_qf20,
    drift_473,
    bo10_th20,
    drift_474,
    bo10_dh20,
    drift_475,
    bo11_bv21,
    drift_476,
    bo11_sxd21,
    drift_477,
    bo11_qd21,
    drift_478,
    bo11_tv21,
    drift_479,
    bo11_dh20,
    drift_480,
    bo11_bh20,
    drift_481,
    bo11_sxf20,
    drift_482,
    bo11_qf20,
    drift_483,
    bo11_th20,
    drift_484,
    bo11_dh19,
    drift_485,
    bo11_bv19,
    drift_486,
    bo11_sxd19,
    drift_487,
    bo11_qd19,
    drift_488,
    bo11_tv19,
    drift_489,
    bo11_dh18,
    drift_490,
    bo11_bh18,
    drift_491,
    bo11_sxf18,
    drift_492,
    bo11_qf18,
    drift_493,
    bo11_th18,
    drift_494,
    bo11_dh17,
    drift_495,
    bo11_bv17,
    drift_496,
    bo11_sxd17,
    drift_497,
    bo11_qd17,
    drift_498,
    bo11_tv17,
    drift_499,
    bo11_dh16,
    drift_500,
    bo11_bh16,
    drift_501,
    bo11_sxf16,
    drift_502,
    bo11_qf16,
    drift_503,
    bo11_th16,
    drift_504,
    bo11_dh15,
    drift_505,
    bo11_bv15,
    drift_506,
    bo11_sxd15,
    drift_507,
    bo11_qd15,
    drift_508,
    bo11_tv15,
    drift_509,
    bo11_dh14,
    drift_510,
    bo11_bh14,
    drift_511,
    bo11_sxf14,
    drift_512,
    bo11_qf14,
    drift_513,
    bo11_th14,
    drift_514,
    bo11_dh13,
    drift_515,
    bo11_bv13,
    drift_516,
    bo11_sxd13,
    drift_517,
    bo11_qd13,
    drift_518,
    bo11_tv13,
    drift_519,
    bo11_dh12,
    drift_520,
    bo11_bh12,
    drift_521,
    bo11_sxf12,
    drift_522,
    bo11_qf12,
    drift_523,
    bo11_th12,
    drift_524,
    bo11_dh11,
    drift_525,
    bo11_bv11,
    drift_526,
    bo11_sxd11,
    drift_527,
    bo11_qd11,
    drift_528,
    bo11_tv11,
    drift_529,
    bo11_dh10,
    drift_530,
    bo11_bh10,
    drift_531,
    bo11_sxf10,
    drift_532,
    bo11_qf10,
    drift_533,
    bo11_th10,
    bo11_oct10,
    bo11_dec10,
    bo11_qs10,
    drift_534,
    bo11_dh9,
    drift_535,
    bo11_bv9,
    drift_536,
    bo11_sxd9,
    drift_537,
    bo11_qd9,
    drift_538,
    bo11_tv9,
    bo11_oct9,
    bo11_dec9,
    bo11_qs9,
    drift_539,
    bo11_dh8,
    drift_540,
    bo11_th8,
    bo11_oct8,
    bo11_dec8,
    bo11_qgt8,
    drift_541,
    bo11_qf8,
    drift_542,
    bo11_b8,
    drift_543,
    bo11_snk_hlx4,
    drift_544,
    bo11_snk_hlx3,
    drift_545,
    bo11_bsnk,
    drift_546,
    bo11_snk_hlx2,
    drift_547,
    bo11_snk_hlx1,
    drift_548,
    bo11_b7,
    drift_549,
    bo11_qd7,
    drift_550,
    bo11_tv7,
    bo11_oct7,
    bo11_dec7,
    bo11_qs7,
    drift_551,
    bo11_dh6,
    drift_552,
    bo11_th6,
    bo11_oct6,
    bo11_dec6,
    bo11_qgt6,
    drift_553,
    bo11_qf6,
    drift_554,
    bo11_tq6,
    drift_555,
    bo11_bh6,
    drift_556,
    bo11_dh5,
    drift_557,
    bo11_tv5,
    bo11_oct5,
    bo11_dec5,
    bo11_qs5,
    drift_558,
    bo11_qd5,
    drift_559,
    bo11_tq5,
    drift_560,
    bo11_bv5,
    drift_561,
    bo11_th4,
    bo11_oct4,
    bo11_dec4,
    bo11_qs4,
    drift_562,
    bo11_qf4,
    drift_563,
    bo11_tq4,
    drift_564,
    bo11_b4,
    drift_565,
    bo11_sv4,
    drift_566,
    bo11_mskh3,
    drift_567,
    scol_h3_2,
    drift_568,
    bo11_kfbh3,
    drift_569,
    bo11_sv3,
    drift_570,
    bo11_b3,
    drift_571,
    bo11_tv3,
    bo11_sx3,
    bo11_oct3,
    bo11_dod3,
    drift_572,
    bo11_qd3,
    drift_573,
    bo11_dods3,
    bo11_octs3,
    bo11_sxs3,
    bo11_qs3,
    drift_574,
    bo11_qf2,
    drift_575,
    bo11_th2,
    bo11_oct2,
    bo11_dec2,
    bo11_dod2,
    drift_576,
    bo11_qd1,
    drift_577,
    bo11_b1,
    drift_578,
    dwarm_ir12h,
    drift_579,
    dsw_ir12h,
    drift_580,
    bi12_b1,
    drift_581,
    bi12_qf1,
    drift_582,
    bi12_tv2,
    bi12_oct2,
    bi12_dec2,
    bi12_dod2,
    drift_583,
    bi12_qd2,
    drift_584,
    bi12_dods3,
    bi12_octs3,
    bi12_sxs3,
    bi12_qs3,
    drift_585,
    bi12_qf3,
    drift_586,
    bi12_th3,
    bi12_sx3,
    bi12_oct3,
    bi12_dod3,
    drift_587,
    bi12_b3,
    drift_588,
    bi12_sv3_1,
    drift_589,
    bi12_kfbh3,
    drift_590,
    bi12_ipm3,
    drift_591,
    scol_h3_1,
    drift_592,
    bi12_eld3,
    drift_593,
    bi12_ksch3_1,
    drift_594,
    bi12_ksch3_2,
    drift_595,
    pcol_h4,
    drift_596,
    bi12_kscv3_1,
    drift_597,
    bi12_kscv3_2,
    drift_598,
    bi12_sv3_2,
    drift_599,
    bi12_pol3_1,
    drift_600,
    bi12_pol3_2,
    drift_601,
    bi12_sv4,
    drift_602,
    bi12_b4,
    drift_603,
    bi12_tq4,
    drift_604,
    bi12_qd4,
    drift_605,
    bi12_tv4,
    bi12_oct4,
    bi12_dec4,
    bi12_qs4,
    drift_606,
    bi12_bh5,
    drift_607,
    bi12_tq5,
    drift_608,
    bi12_qf5,
    drift_609,
    bi12_th5,
    bi12_oct5,
    bi12_dec5,
    bi12_qgt5,
    drift_610,
    bi12_dh5,
    drift_611,
    bi12_bv6,
    drift_612,
    bi12_tq6,
    drift_613,
    bi12_qd6,
    drift_614,
    bi12_tv6,
    bi12_oct6,
    bi12_dec6,
    bi12_qs6,
    drift_615,
    bi12_dh6,
    drift_616,
    bi12_th7,
    bi12_oct7,
    bi12_dec7,
    bi12_qgt7,
    drift_617,
    bi12_qf7,
    drift_618,
    bi12_b7,
    drift_619,
    bi12_b8,
    drift_620,
    bi12_qd8,
    drift_621,
    bi12_tv8,
    bi12_oct8,
    bi12_dec8,
    bi12_qs8,
    drift_622,
    bi12_dh8,
    drift_623,
    bi12_th9,
    bi12_oct9,
    bi12_dec9,
    bi12_qs9,
    drift_624,
    bi12_qf9,
    drift_625,
    bi12_sxf9,
    drift_626,
    bi12_bh9,
    drift_627,
    bi12_dh9,
    drift_628,
    bi12_tv10,
    bi12_oct10,
    bi12_dec10,
    bi12_qs10,
    drift_629,
    bi12_qd10,
    drift_630,
    bi12_sxd10,
    drift_631,
    bi12_bv10,
    drift_632,
    bi12_dh10,
    drift_633,
    bi12_th11,
    bi12_oct11,
    bi12_dec11,
    bi12_qgt11,
    drift_634,
    bi12_qf11,
    drift_635,
    bi12_sxf11,
    drift_636,
    bi12_bh11,
    drift_637,
    bi12_dh11,
    drift_638,
    bi12_tv12,
    bi12_oct12,
    bi12_dec12,
    bi12_qs12,
    drift_639,
    bi12_qd12,
    drift_640,
    bi12_sxd12,
    drift_641,
    bi12_bv12,
    drift_642,
    bi12_dh12,
    drift_643,
    bi12_th13,
    bi12_oct13,
    bi12_dec13,
    bi12_qgt13,
    drift_644,
    bi12_qf13,
    drift_645,
    bi12_sxf13,
    drift_646,
    bi12_bh13,
    drift_647,
    bi12_dh13,
    drift_648,
    bi12_tv14,
    bi12_oct14,
    bi12_dec14,
    bi12_qs14,
    drift_649,
    bi12_qd14,
    drift_650,
    bi12_sxd14,
    drift_651,
    bi12_bv14,
    drift_652,
    bi12_dh14,
    drift_653,
    bi12_th15,
    bi12_oct15,
    bi12_dec15,
    bi12_qgt15,
    drift_654,
    bi12_qf15,
    drift_655,
    bi12_sxf15,
    drift_656,
    bi12_bh15,
    drift_657,
    bi12_dh15,
    drift_658,
    bi12_tv16,
    bi12_oct16,
    bi12_dec16,
    bi12_qs16,
    drift_659,
    bi12_qd16,
    drift_660,
    bi12_sxd16,
    drift_661,
    bi12_bv16,
    drift_662,
    bi12_dh16,
    drift_663,
    bi12_th17,
    bi12_oct17,
    bi12_dec17,
    bi12_qgt17,
    drift_664,
    bi12_qf17,
    drift_665,
    bi12_sxf17,
    drift_666,
    bi12_bh17,
    drift_667,
    bi12_dh17,
    drift_668,
    bi12_tv18,
    bi12_oct18,
    bi12_dec18,
    bi12_qs18,
    drift_669,
    bi12_qd18,
    drift_670,
    bi12_sxd18,
    drift_671,
    bi12_bv18,
    drift_672,
    bi12_dh18,
    drift_673,
    bi12_th19,
    drift_674,
    bi12_qf19,
    drift_675,
    bi12_sxf19,
    drift_676,
    bi12_bh19,
    drift_677,
    bi12_dh19,
    drift_678,
    bi12_tv20,
    drift_679,
    bi12_qd20,
    drift_680,
    bi12_sxd20,
    drift_681,
    bi12_bv20,
    drift_682,
    bi12_dh20,
    drift_683,
    bi1_th21,
    drift_684,
    bi1_qf21,
    drift_685,
    bi1_sxf21,
    drift_686,
    bi1_bh21,
    drift_687,
    bi1_dh20,
    drift_688,
    bi1_tv20,
    drift_689,
    bi1_qd20,
    drift_690,
    bi1_sxd20,
    drift_691,
    bi1_bv20,
    drift_692,
    bi1_dh19,
    drift_693,
    bi1_th19,
    drift_694,
    bi1_qf19,
    drift_695,
    bi1_sxf19,
    drift_696,
    bi1_bh19,
    drift_697,
    bi1_dh18,
    drift_698,
    bi1_tv18,
    drift_699,
    bi1_qd18,
    drift_700,
    bi1_sxd18,
    drift_701,
    bi1_bv18,
    drift_702,
    bi1_dh17,
    drift_703,
    bi1_th17,
    drift_704,
    bi1_qf17,
    drift_705,
    bi1_sxf17,
    drift_706,
    bi1_bh17,
    drift_707,
    bi1_dh16,
    drift_708,
    bi1_tv16,
    drift_709,
    bi1_qd16,
    drift_710,
    bi1_sxd16,
    drift_711,
    bi1_bv16,
    drift_712,
    bi1_dh15,
    drift_713,
    bi1_th15,
    drift_714,
    bi1_qf15,
    drift_715,
    bi1_sxf15,
    drift_716,
    bi1_bh15,
    drift_717,
    bi1_dh14,
    drift_718,
    bi1_tv14,
    drift_719,
    bi1_qd14,
    drift_720,
    bi1_sxd14,
    drift_721,
    bi1_bv14,
    drift_722,
    bi1_dh13,
    drift_723,
    bi1_th13,
    drift_724,
    bi1_qf13,
    drift_725,
    bi1_sxf13,
    drift_726,
    bi1_bh13,
    drift_727,
    bi1_dh12,
    drift_728,
    bi1_tv12,
    drift_729,
    bi1_qd12,
    drift_730,
    bi1_sxd12,
    drift_731,
    bi1_bv12,
    drift_732,
    bi1_dh11,
    drift_733,
    bi1_th11,
    drift_734,
    bi1_qf11,
    drift_735,
    bi1_sxf11,
    drift_736,
    bi1_bh11,
    drift_737,
    bi1_dh10,
    drift_738,
    bi1_tv10,
    bi1_oct10,
    bi1_dec10,
    bi1_qs10,
    drift_739,
    bi1_qd10,
    drift_740,
    bi1_sxd10,
    drift_741,
    bi1_bv10,
    drift_742,
    bi1_dh9,
    drift_743,
    bi1_bh9,
    drift_744,
    bi1_qf9,
    drift_745,
    bi1_th9,
    bi1_oct9,
    bi1_dec9,
    bi1_qs9,
    drift_746,
    bi1_dh8,
    drift_747,
    bi1_tv8,
    bi1_oct8,
    bi1_dec8,
    bi1_qs8,
    drift_748,
    bi1_qd8,
    drift_749,
    bi1_b8,
    drift_750,
    bi1_snk_hlx4,
    drift_751,
    bi1_snk_hlx3,
    drift_752,
    bi1_bsnk,
    drift_753,
    bi1_snk_hlx2,
    drift_754,
    bi1_snk_hlx1,
    drift_755,
    yo1_th6,
    yo1_oct6,
    yo1_dec6,
    yo1_qgt6,
    drift_756,
    yo1_qf6,
    drift_757,
    yo1_tq6,
    drift_758,
    yo1_bh6,
    drift_759,
    bi1_tv6,
    bi1_oct6,
    bi1_dec6,
    bi1_qs6,
    drift_760,
    bi1_qd6,
    drift_761,
    bi1_tq6,
    drift_762,
    bi1_bv6,
    drift_763,
    bi1_dh6,
    drift_764,
    cool_kicker_out,
    cool_kicker_in,
    drift_765,
    yo1_th4,
    yo1_oct4,
    yo1_dec4,
    yo1_qs4,
    drift_766,
    yo1_qf4,
    drift_767,
    yo1_tq4,
    drift_768,
    yo1_b4,
    drift_769,
    bi1_th5,
    bi1_oct5,
    bi1_dec5,
    bi1_qgt5,
    drift_770,
    bi1_qf5,
    drift_771,
    bi1_tq5,
    drift_772,
    bi1_bh5,
    drift_773,
    d5,
    drift_774,
    yo1_bv9,
    drift_775,
    yo1_qd9,
    drift_776,
    yo1_tv9,
    yo1_oct9,
    yo1_dec9,
    yo1_qs9,
    drift_777,
    yo1_tv5,
    yo1_oct5,
    yo1_dec5,
    yo1_qs5,
    drift_778,
    yo1_qd5,
    drift_779,
    yo1_tq5,
    drift_780,
    yo1_bv5,
    drift_781,
    bi1_b7,
    drift_782,
    bi1_qf7,
    drift_783,
    bi1_th7,
    bi1_oct7,
    bi1_dec7,
    bi1_qgt7,
    drift_784,
    bo2_bv5,
    drift_785,
    bo2_tq5,
    drift_786,
    bo2_qd5,
    drift_787,
    bo2_tv5,
    bo2_oct5,
    bo2_dec5,
    bo2_qs5,
    drift_788,
    bo2_tv9,
    bo2_oct9,
    bo2_dec9,
    bo2_qs9,
    drift_789,
    bo2_qd9,
    drift_790,
    bo2_bv9,
    drift_791,
    d5,
    drift_792,
    yi2_bh5,
    drift_793,
    yi2_tq5,
    drift_794,
    yi2_qf5,
    drift_795,
    yi2_th5,
    yi2_oct5,
    yi2_dec5,
    yi2_qgt5,
    drift_796,
    bo2_b4,
    drift_797,
    bo2_tq4,
    drift_798,
    bo2_qf4,
    drift_799,
    bo2_th4,
    bo2_oct4,
    bo2_dec4,
    bo2_qs4,
    drift_800,
    cool_modulator_in,
    cool_modulator_out,
    drift_801,
    yi2_dh6,
    drift_802,
    yi2_bv6,
    drift_803,
    yi2_tq6,
    drift_804,
    yi2_qd6,
    drift_805,
    yi2_tv6,
    yi2_oct6,
    yi2_dec6,
    yi2_qs6,
    drift_806,
    bo2_bh6,
    drift_807,
    bo2_tq6,
    drift_808,
    bo2_qf6,
    drift_809,
    bo2_th6,
    bo2_oct6,
    bo2_dec6,
    bo2_qgt6,
    drift_810,
    yi2_b8,
    drift_811,
    yi2_qd8,
    drift_812,
    yi2_tv8,
    yi2_oct8,
    yi2_dec8,
    yi2_qs8,
    drift_813,
    yi2_dh8,
    drift_814,
    yi2_th9,
    yi2_oct9,
    yi2_dec9,
    yi2_qs9,
    drift_815,
    yi2_qf9,
    drift_816,
    yi2_bh9,
    drift_817,
    yi2_dh9,
    drift_818,
    yi2_bv10,
    drift_819,
    yi2_sxd10,
    drift_820,
    yi2_qd10,
    drift_821,
    yi2_tv10,
    yi2_oct10,
    yi2_dec10,
    yi2_qs10,
    drift_822,
    yi2_dh10,
    drift_823,
    yi2_bh11,
    drift_824,
    yi2_sxf11,
    drift_825,
    yi2_qf11,
    drift_826,
    yi2_th11,
    drift_827,
    yi2_dh11,
    drift_828,
    yi2_bv12,
    drift_829,
    yi2_sxd12,
    drift_830,
    yi2_qd12,
    drift_831,
    yi2_tv12,
    drift_832,
    yi2_dh12,
    drift_833,
    yi2_bh13,
    drift_834,
    yi2_sxf13,
    drift_835,
    yi2_qf13,
    drift_836,
    yi2_th13,
    drift_837,
    yi2_dh13,
    drift_838,
    yi2_bv14,
    drift_839,
    yi2_sxd14,
    drift_840,
    yi2_qd14,
    drift_841,
    yi2_tv14,
    drift_842,
    yi2_dh14,
    drift_843,
    yi2_bh15,
    drift_844,
    yi2_sxf15,
    drift_845,
    yi2_qf15,
    drift_846,
    yi2_th15,
    drift_847,
    yi2_dh15,
    drift_848,
    yi2_bv16,
    drift_849,
    yi2_sxd16,
    drift_850,
    yi2_qd16,
    drift_851,
    yi2_tv16,
    drift_852,
    yi2_dh16,
    drift_853,
    yi2_bh17,
    drift_854,
    yi2_sxf17,
    drift_855,
    yi2_qf17,
    drift_856,
    yi2_th17,
    drift_857,
    yi2_dh17,
    drift_858,
    yi2_bv18,
    drift_859,
    yi2_sxd18,
    drift_860,
    yi2_qd18,
    drift_861,
    yi2_tv18,
    drift_862,
    yi2_dh18,
    drift_863,
    yi2_bh19,
    drift_864,
    yi2_sxf19,
    drift_865,
    yi2_qf19,
    drift_866,
    yi2_th19,
    drift_867,
    yi2_dh19,
    drift_868,
    yi2_bv20,
    drift_869,
    yi2_sxd20,
    drift_870,
    yi2_qd20,
    drift_871,
    yi2_tv20,
    drift_872,
    yi2_dh20,
    drift_873,
    yi3_bh21,
    drift_874,
    yi3_sxf21,
    drift_875,
    yi3_qf21,
    drift_876,
    yi3_th21,
    drift_877,
    yi3_dh20,
    drift_878,
    yi3_bv20,
    drift_879,
    yi3_sxd20,
    drift_880,
    yi3_qd20,
    drift_881,
    yi3_tv20,
    drift_882,
    yi3_dh19,
    drift_883,
    yi3_bh19,
    drift_884,
    yi3_sxf19,
    drift_885,
    yi3_qf19,
    drift_886,
    yi3_th19,
    drift_887,
    yi3_dh18,
    drift_888,
    yi3_bv18,
    drift_889,
    yi3_sxd18,
    drift_890,
    yi3_qd18,
    drift_891,
    yi3_tv18,
    yi3_oct18,
    yi3_dec18,
    yi3_qs18,
    drift_892,
    yi3_dh17,
    drift_893,
    yi3_bh17,
    drift_894,
    yi3_sxf17,
    drift_895,
    yi3_qf17,
    drift_896,
    yi3_th17,
    yi3_oct17,
    yi3_dec17,
    yi3_qgt17,
    drift_897,
    yi3_dh16,
    drift_898,
    yi3_bv16,
    drift_899,
    yi3_sxd16,
    drift_900,
    yi3_qd16,
    drift_901,
    yi3_tv16,
    yi3_oct16,
    yi3_dec16,
    yi3_qs16,
    drift_902,
    yi3_dh15,
    drift_903,
    yi3_bh15,
    drift_904,
    yi3_sxf15,
    drift_905,
    yi3_qf15,
    drift_906,
    yi3_th15,
    yi3_oct15,
    yi3_dec15,
    yi3_qgt15,
    drift_907,
    yi3_dh14,
    drift_908,
    yi3_bv14,
    drift_909,
    yi3_sxd14,
    drift_910,
    yi3_qd14,
    drift_911,
    yi3_tv14,
    yi3_oct14,
    yi3_dec14,
    yi3_qs14,
    drift_912,
    yi3_dh13,
    drift_913,
    yi3_bh13,
    drift_914,
    yi3_sxf13,
    drift_915,
    yi3_qf13,
    drift_916,
    yi3_th13,
    yi3_oct13,
    yi3_dec13,
    yi3_qgt13,
    drift_917,
    yi3_dh12,
    drift_918,
    yi3_bv12,
    drift_919,
    yi3_sxd12,
    drift_920,
    yi3_qd12,
    drift_921,
    yi3_tv12,
    yi3_oct12,
    yi3_dec12,
    yi3_qs12,
    drift_922,
    yi3_dh11,
    drift_923,
    yi3_bh11,
    drift_924,
    yi3_sxf11,
    drift_925,
    yi3_qf11,
    drift_926,
    yi3_th11,
    yi3_oct11,
    yi3_dec11,
    yi3_qgt11,
    drift_927,
    yi3_dh10,
    drift_928,
    yi3_bv10,
    drift_929,
    yi3_sxd10,
    drift_930,
    yi3_qd10,
    drift_931,
    yi3_tv10,
    yi3_oct10,
    yi3_dec10,
    yi3_qs10,
    drift_932,
    yi3_dh9,
    drift_933,
    yi3_bh9,
    drift_934,
    yi3_sxf9,
    drift_935,
    yi3_qf9,
    drift_936,
    yi3_th9,
    yi3_oct9,
    yi3_dec9,
    yi3_qs9,
    drift_937,
    yi3_dh8,
    drift_938,
    yi3_tv8,
    yi3_oct8,
    yi3_dec8,
    yi3_qs8,
    drift_939,
    yi3_qd8,
    drift_940,
    yi3_b8,
    drift_941,
    yi3_hlx7_4,
    drift_942,
    yi3_hlx7_3,
    drift_943,
    yi3_b7_1,
    drift_944,
    yi3_hlx7_2,
    drift_945,
    yi3_hlx7_1,
    drift_946,
    yi3_b7,
    drift_947,
    yi3_qf7,
    drift_948,
    yi3_th7,
    yi3_oct7,
    yi3_dec7,
    yi3_qgt7,
    drift_949,
    yi3_dh6,
    drift_950,
    yi3_tv6,
    yi3_oct6,
    yi3_dec6,
    yi3_qs6,
    drift_951,
    yi3_qd6,
    drift_952,
    yi3_tq6,
    drift_953,
    yi3_bv6,
    drift_954,
    yi3_dh5,
    drift_955,
    yi3_th5,
    yi3_oct5,
    yi3_dec5,
    yi3_qgt5,
    drift_956,
    yi3_qf5,
    drift_957,
    yi3_tq5,
    drift_958,
    yi3_bh5,
    drift_959,
    yi3_tv4,
    yi3_oct4,
    yi3_dec4,
    yi3_qs4,
    drift_960,
    yi3_qd4,
    drift_961,
    yi3_tq4,
    drift_962,
    yi3_b4,
    drift_963,
    yi3_sv4,
    drift_964,
    yi3_kscv3,
    drift_965,
    yi3_ksch3_2,
    drift_966,
    yi3_ksch3_1,
    drift_967,
    cav197_5,
    drift_968,
    cav197_4,
    drift_969,
    cav197_3,
    drift_970,
    cav197_2,
    drift_971,
    cav197_1,
    drift_972,
    cav28_2,
    drift_973,
    cav28_1,
    drift_974,
    yi3_kfbh3,
    drift_975,
    warm_KKQUAD3,
    drift_976,
    yi3_sv3,
    drift_977,
    yi3_b3,
    drift_978,
    yi3_th3,
    yi3_sx3,
    yi3_oct3,
    yi3_dod3,
    drift_979,
    yi3_qf3,
    drift_980,
    yi3_dods3,
    yi3_octs3,
    yi3_sxs3,
    yi3_qs3,
    drift_981,
    yi3_qd2,
    drift_982,
    yi3_tv2,
    yi3_oct2,
    yi3_dec2,
    yi3_dod2,
    drift_983,
    yi3_qf1,
    drift_984,
    yi3_b1,
    drift_985,
    dwarm3,
    drift_986,
    sl_kick_mod1b,
    drift_987,
    sl_kick_mod2b,
    drift_988,
    sl_kick_mod3b,
    drift_989,
    sl_kick_mod4b,
    drift_990,
    sl_kick_mod5b,
    drift_991,
    sl_kick_mod6b,
    drift_992,
    sl_kick_mod7b,
    drift_993,
    sl_kick_mod8b,
    drift_994,
    sl_kick_mod9b,
    drift_995,
    sl_kick_mod10b,
    drift_996,
    ip4,
    drift_997,
    sl_kick_mod1a,
    drift_998,
    sl_kick_mod2a,
    drift_999,
    sl_kick_mod3a,
    drift_1000,
    sl_kick_mod4a,
    drift_1001,
    sl_kick_mod5a,
    drift_1002,
    sl_kick_mod6a,
    drift_1003,
    sl_kick_mod7a,
    drift_1004,
    sl_kick_mod8a,
    drift_1005,
    sl_kick_mod9a,
    drift_1006,
    sl_kick_mod10a,
    drift_1007,
    dwarm4,
    drift_1008,
    det_pol_hpol_hjet_1,
    det_pol_hpol_hjet,
    det_pol_hpol_hjet_1,
    drift_1009,
    yo4_b1,
    drift_1010,
    yo4_qd1,
    drift_1011,
    yo4_th2,
    yo4_oct2,
    yo4_dec2,
    yo4_dod2,
    drift_1012,
    yo4_qf2,
    drift_1013,
    yo4_dods3,
    yo4_octs3,
    yo4_sxs3,
    yo4_qs3,
    drift_1014,
    yo4_qd3,
    drift_1015,
    yo4_tv3,
    yo4_sx3,
    yo4_oct3,
    yo4_dod3,
    drift_1016,
    yo4_b3,
    drift_1017,
    yo4_sv3_1,
    drift_1018,
    warm_KKQUAD4,
    drift_1019,
    yo4_kfbh3,
    drift_1020,
    yo4_wcm3,
    drift_1021,
    septum_ir4,
    drift_1022,
    det_pol_hpol_pc2_1,
    det_pol_hpol_pc2,
    det_pol_hpol_pc2_1,
    drift_1023,
    det_pol_hpol_pc1_1,
    det_pol_hpol_pc1,
    det_pol_hpol_pc1_1,
    drift_1024,
    yo4_sv4,
    drift_1025,
    yo4_b4,
    drift_1026,
    yo4_tq4,
    drift_1027,
    yo4_qf4,
    drift_1028,
    yo4_th4,
    yo4_oct4,
    yo4_dec4,
    yo4_qs4,
    drift_1029,
    yo4_bv5,
    drift_1030,
    yo4_tq5,
    drift_1031,
    yo4_qd5,
    drift_1032,
    yo4_tv5,
    yo4_oct5,
    yo4_dec5,
    yo4_qs5,
    drift_1033,
    yo4_dh5,
    drift_1034,
    yo4_bh6,
    drift_1035,
    yo4_tq6,
    drift_1036,
    yo4_qf6,
    drift_1037,
    yo4_th6,
    yo4_oct6,
    yo4_dec6,
    yo4_qgt6,
    drift_1038,
    yo4_dh6,
    drift_1039,
    yo4_tv7,
    yo4_oct7,
    yo4_dec7,
    yo4_qs7,
    drift_1040,
    yo4_qd7,
    drift_1041,
    yo4_b7,
    drift_1042,
    yo4_b8,
    drift_1043,
    yo4_qf8,
    drift_1044,
    yo4_th8,
    yo4_oct8,
    yo4_dec8,
    yo4_qgt8,
    drift_1045,
    yo4_dh8,
    drift_1046,
    yo4_tv9,
    yo4_oct9,
    yo4_dec9,
    yo4_qs9,
    drift_1047,
    yo4_qd9,
    drift_1048,
    yo4_sxd9,
    drift_1049,
    yo4_bv9,
    drift_1050,
    yo4_dh9,
    drift_1051,
    yo4_th10,
    yo4_oct10,
    yo4_dec10,
    yo4_qs10,
    drift_1052,
    yo4_qf10,
    drift_1053,
    yo4_sxf10,
    drift_1054,
    yo4_bh10,
    drift_1055,
    yo4_dh10,
    drift_1056,
    yo4_tv11,
    drift_1057,
    yo4_qd11,
    drift_1058,
    yo4_sxd11,
    drift_1059,
    yo4_bv11,
    drift_1060,
    yo4_dh11,
    drift_1061,
    yo4_th12,
    drift_1062,
    yo4_qf12,
    drift_1063,
    yo4_sxf12,
    drift_1064,
    yo4_bh12,
    drift_1065,
    yo4_dh12,
    drift_1066,
    yo4_tv13,
    drift_1067,
    yo4_qd13,
    drift_1068,
    yo4_sxd13,
    drift_1069,
    yo4_bv13,
    drift_1070,
    yo4_dh13,
    drift_1071,
    yo4_th14,
    drift_1072,
    yo4_qf14,
    drift_1073,
    yo4_sxf14,
    drift_1074,
    yo4_bh14,
    drift_1075,
    yo4_dh14,
    drift_1076,
    yo4_tv15,
    drift_1077,
    yo4_qd15,
    drift_1078,
    yo4_sxd15,
    drift_1079,
    yo4_bv15,
    drift_1080,
    yo4_dh15,
    drift_1081,
    yo4_th16,
    drift_1082,
    yo4_qf16,
    drift_1083,
    yo4_sxf16,
    drift_1084,
    yo4_bh16,
    drift_1085,
    yo4_dh16,
    drift_1086,
    yo4_tv17,
    drift_1087,
    yo4_qd17,
    drift_1088,
    yo4_sxd17,
    drift_1089,
    yo4_bv17,
    drift_1090,
    yo4_dh17,
    drift_1091,
    yo4_th18,
    drift_1092,
    yo4_qf18,
    drift_1093,
    yo4_sxf18,
    drift_1094,
    yo4_bh18,
    drift_1095,
    yo4_dh18,
    drift_1096,
    yo4_tv19,
    drift_1097,
    yo4_qd19,
    drift_1098,
    yo4_sxd19,
    drift_1099,
    yo4_bv19,
    drift_1100,
    yo4_dh19,
    drift_1101,
    yo4_th20,
    drift_1102,
    yo4_qf20,
    drift_1103,
    yo4_sxf20,
    drift_1104,
    yo4_bh20,
    drift_1105,
    yo4_dh20,
    drift_1106,
    yo5_tv21,
    drift_1107,
    yo5_qd21,
    drift_1108,
    yo5_sxd21,
    drift_1109,
    yo5_bv21,
    drift_1110,
    yo5_dh20,
    drift_1111,
    yo5_th20,
    drift_1112,
    yo5_qf20,
    drift_1113,
    yo5_sxf20,
    drift_1114,
    yo5_bh20,
    drift_1115,
    yo5_dh19,
    drift_1116,
    yo5_tv19,
    drift_1117,
    yo5_qd19,
    drift_1118,
    yo5_sxd19,
    drift_1119,
    yo5_bv19,
    drift_1120,
    yo5_dh18,
    drift_1121,
    yo5_th18,
    yo5_oct18,
    yo5_dec18,
    yo5_qgt18,
    drift_1122,
    yo5_qf18,
    drift_1123,
    yo5_sxf18,
    drift_1124,
    yo5_bh18,
    drift_1125,
    yo5_dh17,
    drift_1126,
    yo5_tv17,
    yo5_oct17,
    yo5_dec17,
    yo5_qs17,
    drift_1127,
    yo5_qd17,
    drift_1128,
    yo5_sxd17,
    drift_1129,
    yo5_bv17,
    drift_1130,
    yo5_dh16,
    drift_1131,
    yo5_th16,
    yo5_oct16,
    yo5_dec16,
    yo5_qgt16,
    drift_1132,
    yo5_qf16,
    drift_1133,
    yo5_sxf16,
    drift_1134,
    yo5_bh16,
    drift_1135,
    yo5_dh15,
    drift_1136,
    yo5_tv15,
    yo5_oct15,
    yo5_dec15,
    yo5_qs15,
    drift_1137,
    yo5_qd15,
    drift_1138,
    yo5_sxd15,
    drift_1139,
    yo5_bv15,
    drift_1140,
    yo5_dh14,
    drift_1141,
    yo5_th14,
    yo5_oct14,
    yo5_dec14,
    yo5_qgt14,
    drift_1142,
    yo5_qf14,
    drift_1143,
    yo5_sxf14,
    drift_1144,
    yo5_bh14,
    drift_1145,
    yo5_dh13,
    drift_1146,
    yo5_tv13,
    yo5_oct13,
    yo5_dec13,
    yo5_qs13,
    drift_1147,
    yo5_qd13,
    drift_1148,
    yo5_sxd13,
    drift_1149,
    yo5_bv13,
    drift_1150,
    yo5_dh12,
    drift_1151,
    yo5_th12,
    yo5_oct12,
    yo5_dec12,
    yo5_qgt12,
    drift_1152,
    yo5_qf12,
    drift_1153,
    yo5_sxf12,
    drift_1154,
    yo5_bh12,
    drift_1155,
    yo5_dh11,
    drift_1156,
    yo5_tv11,
    yo5_oct11,
    yo5_dec11,
    yo5_qs11,
    drift_1157,
    yo5_qd11,
    drift_1158,
    yo5_sxd11,
    drift_1159,
    yo5_bv11,
    drift_1160,
    yo5_dh10,
    drift_1161,
    yo5_th10,
    yo5_oct10,
    yo5_dec10,
    yo5_qs10,
    drift_1162,
    yo5_qf10,
    drift_1163,
    yo5_sxf10,
    drift_1164,
    yo5_bh10,
    drift_1165,
    yo5_sv9_2,
    drift_1166,
    h5_tq11,
    drift_1167,
    h5_q11,
    drift_1168,
    h5_tv11,
    h5_b311,
    h5_b411,
    h5_a111,
    drift_1169,
    yo5_sv9_1,
    drift_1170,
    yo5_dh9,
    drift_1171,
    yo5_bv9,
    drift_1172,
    yo5_qd9,
    drift_1173,
    yo5_tv9,
    yo5_oct9,
    yo5_dec9,
    yo5_qs9,
    drift_1174,
    yo5_dh8,
    drift_1175,
    yo5_th8,
    yo5_oct8,
    yo5_dec8,
    yo5_qgt8,
    drift_1176,
    yo5_qf8,
    drift_1177,
    yo5_b8,
    drift_1178,
    yo5_snk_hlx4,
    drift_1179,
    yo5_snk_hlx3,
    drift_1180,
    yo5_bsnk,
    drift_1181,
    yo5_snk_hlx2,
    drift_1182,
    yo5_snk_hlx1,
    drift_1183,
    yo5_b7,
    drift_1184,
    yo5_qd7,
    drift_1185,
    yo5_tv7,
    yo5_oct7,
    yo5_dec7,
    yo5_qs7,
    drift_1186,
    yo5_rot_hlx4,
    drift_1187,
    yo5_rot_hlx3,
    drift_1188,
    yo5_brot,
    drift_1189,
    yo5_rot_hlx2,
    drift_1190,
    yo5_rot_hlx1,
    drift_1191,
    bi5_tv8,
    bi5_oct8,
    bi5_dec8,
    bi5_qs8,
    drift_1192,
    bi5_qd8,
    drift_1193,
    bi5_b8,
    drift_1194,
    h5_dh4,
    drift_1195,
    yo5_th4,
    yo5_oct4,
    yo5_dec4,
    yo5_qs4,
    drift_1196,
    yo5_qf4,
    drift_1197,
    yo5_tq4,
    drift_1198,
    yo5_b4,
    drift_1199,
    yo5_dh5,
    drift_1200,
    bo3_th8,
    bo3_oct8,
    bo3_dec8,
    bo3_qgt8,
    drift_1201,
    bo3_qf8,
    drift_1202,
    bo3_b8,
    drift_1203,
    yi6_tv2,
    yi6_oct2,
    yi6_dec2,
    yi6_dod2,
    drift_1204,
    yi6_qd2,
    drift_1205,
    o_crab_ip6f_d1,
    o_crab_ip6f,
    o_crab_ip6f_d2,
    drift_1206,
    b2pf,
    drift_1207,
    h5_tv3,
    drift_1208,
    h5_qs3,
    drift_1209,
    q3pf,
    drift_1210,
    o_roman_pot_ip6_2_1,
    o_roman_pot_ip6_2,
    o_roman_pot_ip6_2_1,
    drift_1211,
    o_roman_pot_ip6_1_1,
    o_roman_pot_ip6_1,
    o_roman_pot_ip6_1_1,
    drift_1212,
    o_off_mom_ip6_2_1,
    o_off_mom_ip6_2,
    o_off_mom_ip6_2_1,
    drift_1213,
    o_off_mom_ip6_1_1,
    o_off_mom_ip6_1,
    o_off_mom_ip6_1_1,
    drift_1214,
    pbr_b1apf,
    pbt_b1apf,
    b1apf,
    pet_b1apf,
    per_b1apf,
    drift_1215,
    pbr_b1pf,
    pbt_b1pf,
    b1pf,
    pet_b1pf,
    per_b1pf,
    drift_1216,
    pbr_q2pf,
    pbt_q2pf,
    q2pf,
    pet_q2pf,
    per_q2pf,
    drift_1217,
    pbr_q1bpf,
    pbt_q1bpf,
    q1bpf,
    pet_q1bpf,
    per_q1bpf,
    drift_1218,
    pbr_q1apf,
    pbt_q1apf,
    q1apf,
    pet_q1apf,
    per_q1apf,
    drift_1219,
    pbr_b0apf,
    pbt_b0apf,
    b0apf,
    pet_b0apf,
    per_b0apf,
    drift_1220,
    pbt_b0pf,
    pbr_b0pf,
    b0pf,
    pet_b0pf,
    per_b0pf,
    drift_1221,
    pbt_sol_f,
    pbr_sol_f,
    star_detect_f,
    per_sol_f,
    ip6w,
]
