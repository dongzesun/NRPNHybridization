"""Generated octupole and current quadrupole source jets."""

import math

GRADES = (0, 2, 4, 5, 6)

def i3(potential, radial_speed, speed_squared, nu, log_one, log_infrared):
    common_0 = nu*potential
    common_1 = (13/6)*common_0
    common_2 = 2*nu
    common_3 = common_2 - 1
    common_4 = potential*speed_squared
    common_5 = potential**2
    common_6 = speed_squared**2
    common_7 = common_0*speed_squared
    common_8 = common_5*nu
    common_9 = common_6*nu
    common_10 = radial_speed**2
    common_11 = common_10*potential
    common_12 = common_0*common_10
    common_13 = nu**2
    common_14 = common_13*common_4
    common_15 = common_13*common_5
    common_16 = common_13*common_6
    common_17 = common_11*common_13
    common_18 = common_0*radial_speed
    common_19 = nu*speed_squared
    common_20 = common_19*radial_speed
    common_21 = potential*radial_speed
    common_22 = common_13*common_21
    common_23 = common_10*nu
    common_24 = common_13*speed_squared
    common_25 = common_13*radial_speed
    common_26 = common_8*radial_speed
    common_27 = potential**3
    common_28 = speed_squared**3
    common_29 = common_27*log_infrared
    common_30 = common_27*nu
    common_31 = common_28*nu
    common_32 = radial_speed**4
    common_33 = common_32*potential
    common_34 = common_6*potential
    common_35 = common_5*speed_squared
    common_36 = common_30*log_one
    common_37 = common_0*common_32
    common_38 = common_0*common_6
    common_39 = common_19*common_5
    common_40 = common_10*common_4
    common_41 = common_13*common_27
    common_42 = common_13*common_28
    common_43 = nu**3
    common_44 = common_27*common_43
    common_45 = common_28*common_43
    common_46 = common_10*common_5
    common_47 = common_10*speed_squared
    common_48 = common_0*common_47
    common_49 = common_10*common_8
    common_50 = common_13*potential
    common_51 = common_32*common_50
    common_52 = common_16*potential
    common_53 = common_15*speed_squared
    common_54 = common_33*common_43
    common_55 = common_34*common_43
    common_56 = common_35*common_43
    common_57 = common_10*common_13
    common_58 = common_4*common_57
    common_59 = common_40*common_43
    common_60 = common_10*common_15
    common_61 = common_43*common_46
    common_62 = radial_speed**3
    common_63 = common_62*potential
    common_64 = common_5*radial_speed
    common_65 = radial_speed*speed_squared
    common_66 = common_0*common_65
    common_67 = common_0*common_62
    common_68 = common_9*radial_speed
    common_69 = common_50*common_62
    common_70 = common_43*common_64
    common_71 = common_6*radial_speed
    common_72 = common_43*common_71
    common_73 = (52/7)*common_5*log_infrared
    common_74 = common_10*common_19
    common_75 = common_4*common_43
    common_76 = common_11*common_43
    common_77 = common_43*common_47
    common_78 = common_13*common_65
    common_79 = common_13*common_62
    common_80 = common_4*radial_speed
    common_81 = common_15*radial_speed
    common_82 = -464/15*common_8 + (24/5)*nu*potential*speed_squared
    common_83 = radial_speed**5
    common_84 = (78/7)*common_29
    common_85 = common_0*common_83
    common_86 = common_21*common_6
    common_87 = common_4*common_62
    common_88 = common_5*common_65
    common_89 = common_5*common_62
    common_90 = nu*radial_speed
    common_91 = common_27*common_90
    common_92 = common_91*log_one
    common_93 = common_0*common_71
    common_94 = common_19*common_64
    common_95 = common_25*common_27
    common_96 = common_83*potential
    common_97 = common_43*common_96
    common_98 = common_44*radial_speed
    common_99 = common_4*common_79
    common_100 = common_43*common_6
    common_101 = common_100*common_21
    common_102 = common_15*common_62
    common_103 = common_43*common_5
    common_104 = common_103*common_62
    common_105 = common_25*common_4
    common_106 = common_43*common_80
    common_107 = common_16*radial_speed
    common_108 = common_43*common_63
    common_109 = 3*potential
    common_110 = potential**4
    common_111 = common_110*log_infrared
    common_112 = common_110*nu
    common_113 = radial_speed**6
    common_114 = common_113*potential
    common_115 = common_0*common_28
    common_116 = common_110*common_13
    common_117 = common_10*common_27
    common_118 = common_32*common_7
    common_119 = math.pi**2
    common_120 = common_6*common_8
    common_121 = common_23*common_27
    common_122 = common_113*common_50
    common_123 = common_24*common_27
    common_124 = common_114*common_43
    common_125 = common_45*potential
    common_126 = common_44*speed_squared
    common_127 = common_19*common_27
    common_128 = common_119*common_127
    common_129 = common_12*common_6
    common_130 = common_19*common_46
    common_131 = common_15*common_32
    common_132 = common_27*common_57
    common_133 = common_103*common_32
    common_134 = common_11*common_16
    common_135 = common_10*common_56
    common_136 = common_27*radial_speed
    common_137 = 264*log_one
    common_138 = common_62*common_7
    common_139 = common_62*common_8
    common_140 = common_50*common_83
    common_141 = common_119*common_91
    common_142 = common_16*common_21
    common_143 = common_15*common_65
    common_144 = common_43*common_62
    common_145 = common_144*common_4
    common_146 = common_103*common_65
    common_147 = (520/7)*common_29 + 132*common_36
    common_148 = 12*potential
    common_149 = common_32*common_8
    common_150 = radial_speed**7
    common_151 = common_110*radial_speed
    common_152 = common_0*common_150
    common_153 = common_21*common_28
    common_154 = common_4*common_83
    common_155 = common_27*common_65
    common_156 = common_5*common_83
    common_157 = common_136*common_19
    common_158 = common_29*common_62
    common_159 = common_8*common_83
    common_160 = common_110*common_43
    common_161 = common_160*radial_speed
    common_162 = common_34*common_62
    common_163 = common_35*common_62
    common_164 = common_110*common_90
    common_165 = common_119*common_164
    common_166 = common_38*common_62
    common_167 = common_71*common_8
    common_168 = common_21*common_42
    common_169 = common_14*common_83
    common_170 = common_75*common_83
    common_171 = common_44*common_62
    common_172 = common_30*common_62
    common_173 = common_119*common_172
    common_174 = common_16*common_64
    common_175 = common_53*common_62
    common_176 = common_100*common_64
    common_177 = common_144*common_35
    common_178 = common_28*potential
    common_179 = common_27*speed_squared
    common_180 = common_29*speed_squared
    common_181 = common_0*common_113
    common_182 = common_32*common_4
    common_183 = common_32*common_5
    common_184 = common_28*common_50
    common_185 = common_11*common_6
    common_186 = common_10*common_35
    common_187 = common_16*common_5
    common_188 = common_100*common_5
    common_189 = common_10*common_44
    common_190 = common_15*common_47
    common_191 = common_100*common_11
    common_192 = common_29*radial_speed
    common_193 = 9*common_4
    common_194 = 60*potential
    common_195 = common_19*common_89
    common_196 = potential**5
    common_197 = radial_speed**8
    common_198 = common_197*potential
    common_199 = speed_squared**4
    common_200 = common_196*nu
    common_201 = common_110*common_23
    common_202 = common_150*potential
    return (
        (-1, 0, 0, 0),
        (-common_1 + (19/6)*nu*speed_squared + (5/6)*potential - 5/6*speed_squared, -common_2*radial_speed + radial_speed, common_3, 0),
        ((247/1320)*common_11 - 531/440*common_12 + (17371/1320)*common_14 - 235/66*common_15 - 5501/440*common_16 + (1347/440)*common_17 - 3853/1320*common_4 + (47/33)*common_5 - 257/440*common_6 + (14257/1320)*common_7 + (1591/132)*common_8 + (7319/1320)*common_9, (102/11)*common_13*radial_speed*speed_squared - 8689/660*common_18 - 107/22*common_20 - 1389/220*common_22 + (2461/660)*potential*radial_speed + (13/22)*radial_speed*speed_squared, -62/165*common_0 + (3/11)*common_10*common_13 + (1/11)*common_10 + (483/55)*common_13*potential - 4/11*common_23 - 504/55*common_24 + (519/110)*nu*speed_squared - 1949/330*potential - 61/110*speed_squared, -39/55*common_25 + (52/55)*nu*radial_speed - 13/55*radial_speed),
        ((56/9)*common_26, -232/15*common_8 + (12/5)*nu*potential*speed_squared, 0, 0),
        (-4553429/229320*common_27 - 3235/6864*common_28 + (26/7)*common_29 - 902873/15015*common_30 + (7667/1040)*common_31 - 1/208*common_33 - 11633/2640*common_34 - 40497/20020*common_35 + 22*common_36 + (1379/2288)*common_37 + (93167/2640)*common_38 + (397027/8580)*common_39 + (841/11440)*common_40 + (31673/1716)*common_41 - 10319/286*common_42 - 27085/5148*common_44 + (129707/2288)*common_45 + (60413/120120)*common_46 - 2039/3432*common_48 - 20801/1560*common_49 - 44841/11440*common_51 - 5289/110*common_52 - 120069/2860*common_53 + (63699/11440)*common_54 - 203299/2640*common_55 + (29429/858)*common_56 + (121393/11440)*common_58 - 589981/17160*common_59 + (623351/17160)*common_60 + (70735/3432)*common_61, (1283/330)*common_13*common_5*radial_speed + (34627/1144)*common_13*common_6*radial_speed + (167122/2145)*common_13*potential*radial_speed*speed_squared - 33459/440*common_26 + (299171/17160)*common_43*common_62*potential + (267081/5720)*common_43*potential*radial_speed*speed_squared + (2461/5720)*common_6*radial_speed - 5/572*common_63 - 229/9240*common_64 - 47979/1144*common_66 - 1851/1144*common_67 - 37321/5720*common_68 - 3059/1560*common_69 - 2287/264*common_70 - 127447/2860*common_72 + (80629/17160)*potential*radial_speed*speed_squared, (849/286)*common_10*common_13*speed_squared + (125/104)*common_10*nu*potential + (3037/3432)*common_10*potential + (197/1430)*common_10*speed_squared - 83211/5720*common_14 - 428669/8580*common_15 - 3057/104*common_16 - 62143/17160*common_17 - 91379/17160*common_4 + (62665/2574)*common_43*common_5 + (25071/572)*common_43*common_6 + (75373/8580)*common_5*nu + (34515967/1261260)*common_5 + (35777/5720)*common_6*nu - 465/1144*common_6 - common_73 - 1637/1430*common_74 - 504721/8580*common_75 - 113089/8580*common_76 - 3063/1430*common_77 + (169537/5720)*nu*potential*speed_squared, (1655/156)*common_13*potential*radial_speed - 24367/4290*common_21*common_43 - 26641/8580*common_21 + (448/2145)*common_43*common_62 + (7247/1430)*common_43*radial_speed*speed_squared + (224/715)*common_62*nu - 112/2145*common_62 - 333/1430*common_65 - 1849/286*common_78 - 224/429*common_79 + (8341/2860)*nu*potential*radial_speed + (3181/1430)*nu*radial_speed*speed_squared),
        (0, -3, 0, 0),
        (-common_1*radial_speed - 1/6*common_21, -17/2*common_0 + (15/2)*nu*speed_squared + (7/2)*potential - 3/2*speed_squared, 0, common_3),
        ((10603/1320)*common_13*potential*radial_speed*speed_squared - 1865/132*common_26 + (1153/440)*common_62*nu*potential - 247/440*common_63 - 499/132*common_64 - 5539/1320*common_66 - 1181/440*common_69 - 423/440*common_80 - 455/44*common_81, (4579/1320)*common_11 - 1981/1320*common_12 + (15937/440)*common_14 - 789/44*common_15 - 12423/440*common_16 + (3259/440)*common_17 - 4633/1320*common_4 + (1103/132)*common_5 - 511/440*common_6 + (20677/1320)*common_7 + (3701/66)*common_8 + (5179/440)*common_9, -6023/330*common_18 - 56/55*common_20 - 794/165*common_21 + (887/110)*common_22 + (14/55)*common_65 + (42/55)*common_78, -218/165*common_0 + (84/55)*common_10*nu - 21/55*common_10 + (522/55)*common_13*potential - 543/55*common_24 - 63/55*common_57 + (623/110)*nu*speed_squared - 1871/330*potential - 87/110*speed_squared),
        ((416/45)*common_30 + (172/45)*common_39 - 56/3*common_49, (88/3)*common_26, common_82, 0),
        (-36511/1040*common_101 - 5513/1560*common_102 - 7624/429*common_104 + (60331/1320)*common_13*common_5*radial_speed*speed_squared + (48761/1716)*common_13*common_6*potential*radial_speed + (17671/2288)*common_13*common_83*potential + (25066577/120120)*common_27*nu*radial_speed + (4583147/70070)*common_27*radial_speed + (74213/1320)*common_43*common_5*radial_speed*speed_squared + (8497/440)*common_43*common_62*potential*speed_squared + (280789/8580)*common_5*common_62*nu + (243/40)*common_62*nu*potential*speed_squared + (5/208)*common_83*potential - common_84*radial_speed - 4035/2288*common_85 - 22273/34320*common_86 - 211/880*common_87 - 34343/4620*common_88 - 236297/120120*common_89 - 66*common_92 - 38641/11440*common_93 - 9781/264*common_94 - 43993/1716*common_95 - 7929/2288*common_97 - 113653/3432*common_98 - 3791/176*common_99, -12397261/90090*common_27 - 1023/1040*common_28 + 26*common_29 - 7607753/40040*common_30 + (178369/11440)*common_31 + (27/2288)*common_33 - 261511/34320*common_34 + (352487/60060)*common_35 + 66*common_36 + (3621/2288)*common_37 + (662707/11440)*common_38 + (404803/4290)*common_39 + (4961/3120)*common_40 + (214627/1716)*common_41 - 89201/1144*common_42 - 342977/10296*common_44 + (1435817/11440)*common_45 + (788881/40040)*common_46 - 15831/5720*common_48 - 154659/5720*common_49 - 43319/2288*common_51 - 900499/17160*common_52 - 109001/715*common_53 + (37267/2288)*common_54 - 39991/208*common_55 + (234523/2145)*common_56 + (112795/6864)*common_58 - 457819/5720*common_59 + (2011727/17160)*common_60 + (206407/3432)*common_61, (1578353/17160)*common_105 - 52919/1430*common_106 + (1099/143)*common_107 + (2191/286)*common_108 - 13859/110*common_26 + (1487/3432)*common_63 - 30089/770*common_64 - 11933/520*common_66 - 3533/1144*common_67 - 2023/715*common_68 + (2471/312)*common_69 + (6967/132)*common_70 + (21/65)*common_71 - 4109/715*common_72 + (11233/17160)*common_80 - 30691/660*common_81, -3297/286*common_10*common_24 + (139849/2860)*common_100 + (184087/6435)*common_103 - 2367/440*common_11 + (6773/440)*common_12 + (623/1560)*common_14 - 263903/4290*common_15 - 41023/1144*common_16 + (1757/440)*common_17 - 12419/1560*common_4 - 63/130*common_47 + (338543/11466)*common_5 - 3657/5720*common_6 + (172143/5720)*common_7 - common_73 + (6069/1430)*common_74 - 38579/572*common_75 - 5503/220*common_76 + (12327/1430)*common_77 + (39371/4290)*common_8 + (48501/5720)*common_9),
        (common_109, 0, -6, 0),
        (common_11*common_2 + (1/2)*common_11 + (13/3)*common_4 - 46/3*common_5 - 2/3*common_7 + (14/3)*common_8, -7*common_18 - 13*common_21, -23*common_0 + 15*nu*speed_squared + 10*potential - 3*speed_squared, 4*common_90 - 2*radial_speed),
        ((1180/33)*common_27 - 419/33*common_30 + (247/88)*common_33 + (17/10)*common_34 - 1426/165*common_35 - 329/44*common_37 - 299/30*common_38 - 2876/165*common_39 + (6/5)*common_40 + (229/33)*common_41 + (5671/660)*common_46 + (47/10)*common_48 - 4353/220*common_49 - 38/11*common_51 + (53/30)*common_52 - 263/165*common_53 - 29/10*common_58 + (1821/220)*common_60, (27/22)*common_13*common_62*potential + (181/10)*common_13*potential*radial_speed*speed_squared - 47429/660*common_26 + (235/11)*common_62*nu*potential - 133/11*common_63 - 6829/660*common_64 - 29/5*common_66 - 33/10*common_80 - 8699/220*common_81, (545/44)*common_11 + (361/44)*common_12 + (20099/220)*common_14 - 302/5*common_15 - 2451/44*common_16 + (777/44)*common_17 - 1477/220*common_4 + (398/15)*common_5 - 91/44*common_6 + (23/220)*common_7 + (4571/30)*common_8 + (991/44)*common_9, (607/22)*common_13*potential*radial_speed - 87/22*common_18 - 443/22*common_21 - 23/11*common_65 - 234/11*common_78 + (147/11)*nu*radial_speed*speed_squared),
        ((280/3)*common_5*common_62*nu - 2872/45*common_91 - 316/5*common_94, (224/3)*common_30 + 36*common_39 - 144*common_49, 80*common_26, common_82),
        ((1845/64)*common_10*common_119*common_27*nu + (158423/8580)*common_10*common_13*common_5*speed_squared + (93319/4290)*common_10*common_27*common_43 + (390/7)*common_10*common_27*log_infrared + (24033/2860)*common_10*common_43*common_6*potential + (446821/40040)*common_10*common_5*speed_squared + (3511/2860)*common_10*common_6*potential + (123/16)*common_110*common_119*nu + (12718/1287)*common_110*common_43 + (3942902/315315)*common_110 - 104/7*common_111 - 626281/4290*common_112 + (6615/1144)*common_113*nu*potential - 35/208*common_114 - 7543/520*common_115 - 2689/858*common_116 - 36028931/105105*common_117 - 15879/1144*common_118 - 72037/8580*common_120 - 8674753/17160*common_121 - 24311/1144*common_122 - 238687/4290*common_123 - 4893/572*common_124 - 19771/2860*common_125 - 8707/8580*common_126 - 369/64*common_128 - 7533/2860*common_129 + (584489/17160)*common_13*common_28*potential + (25815/1144)*common_13*common_32*potential*speed_squared + (48899/2145)*common_13*common_5*common_6 - 604613/3432*common_130 - 43529/715*common_131 - 123852/715*common_132 - 100729/2860*common_133 - 3347/5720*common_134 - 155809/17160*common_135 + (5989/260)*common_27*nu*speed_squared + (6097367/105105)*common_27*speed_squared + (3209/2145)*common_28*potential + (13687/572)*common_32*common_43*potential*speed_squared + (185431/2860)*common_32*common_5*nu + (108393/20020)*common_32*common_5 + (1509/1144)*common_32*potential*speed_squared + (57479/8580)*common_43*common_5*common_6 + (19849/30030)*common_5*common_6 - common_84*speed_squared, -214847/2860*common_101 + (507781/4290)*common_102 + (50213/8580)*common_104 + (25403471/38220)*common_136 - common_137*common_91 + (3609/52)*common_138 + (68141/260)*common_139 + (33701/286)*common_140 - 369/32*common_141 + (219339/2860)*common_142 + (210904/2145)*common_143 + (370/143)*common_145 + (178625/1144)*common_146 - 780/7*common_29*radial_speed - 27975/1144*common_85 - 14249/2860*common_86 - 283/52*common_87 - 206645/24024*common_88 - 266153/5005*common_89 + (14692591/10920)*common_91 + (21579/5720)*common_93 - 1836001/17160*common_94 - 1427729/8580*common_95 + (15/1144)*common_96 + (197/52)*common_97 - 204917/1430*common_98 - 118737/572*common_99, common_147 - 121397896/315315*common_27 - 171/104*common_28 - 6048071/15015*common_30 + (32437/1144)*common_31 - 365/286*common_33 - 62187/5720*common_34 + (554821/60060)*common_35 + (19563/1144)*common_37 + (398589/5720)*common_38 + (1057271/8580)*common_39 + (11113/1144)*common_40 + (1293899/2860)*common_41 - 84805/572*common_42 - 1700801/12870*common_44 + (280589/1144)*common_45 + (3891683/60060)*common_46 - 4789/143*common_48 - 200433/2860*common_49 - 8937/104*common_51 + (108263/2860)*common_52 - 353293/780*common_53 + (46309/1144)*common_54 - 2584979/5720*common_55 + (1422029/4290)*common_56 + (65409/1144)*common_58 - 105335/572*common_59 + (329299/780)*common_60 + (121211/780)*common_61, -544763/2860*common_106 - 49815/572*common_107 - 13823/572*common_108 + (4384/143)*common_13*potential*radial_speed*speed_squared - 68974/715*common_26 + (966803/8580)*common_43*common_5*radial_speed + (31257/286)*common_43*common_6*radial_speed + (12937/572)*common_6*nu*radial_speed + (9939/572)*common_62*nu*potential - 733/858*common_63 - 124546/1365*common_64 - 13087/1716*common_69 - 1101/572*common_71 - 7283/260*common_80 - 178207/1716*common_81 + (248433/2860)*nu*potential*radial_speed*speed_squared),
        (-9*common_21, 21*potential, 0, -6),
        (-common_148*common_65 - 43/3*common_26 + (194/3)*common_5*radial_speed - 10*common_62*nu*potential - 5/2*common_63 + 6*nu*potential*radial_speed*speed_squared, -common_109*common_19 + (81/2)*common_11 + 9*common_12 + 18*common_4 - 101*common_5 + 43*common_8, -9*common_18 - 72*common_21, -27*common_0 + 19*nu*speed_squared + 12*potential - 5*speed_squared),
        (-1007/22*common_102 - 23981/132*common_136 - 1175/22*common_138 + (2791/22)*common_139 + (266/11)*common_140 - 111/10*common_142 + (3173/165)*common_143 + (2303/44)*common_85 - 27/10*common_86 + (115/22)*common_87 + (4612/55)*common_88 - 540/11*common_89 + (24523/132)*common_91 + (393/10)*common_93 + (6389/110)*common_94 - 2269/132*common_95 - 1729/88*common_96 + (15/22)*common_99, (1421/44)*common_10*common_13*common_5 + (1269/22)*common_10*common_13*potential*speed_squared + (717/22)*common_10*nu*potential*speed_squared + (3307/44)*common_13*common_27 + (26473/132)*common_27 - 27451/132*common_30 + (551/8)*common_32*potential - 33361/660*common_35 - 427/4*common_37 - 3121/55*common_38 - 37331/660*common_39 - 501/22*common_40 - 5641/132*common_46 - 17729/132*common_49 - 84*common_51 - 351/110*common_52 - 471/220*common_53 + (1313/110)*common_6*potential, -141/44*common_105 - 9478/33*common_26 - 2699/44*common_63 + (1033/33)*common_64 + (1871/44)*common_66 + (2381/44)*common_67 + (2133/44)*common_69 + (151/44)*common_80 - 386/11*common_81, (207/44)*common_11 + (2747/44)*common_12 + (28209/220)*common_14 - 8799/110*common_15 - 3387/44*common_16 + (291/44)*common_17 - 5007/220*common_4 + (12761/330)*common_5 - 183/44*common_6 - 3347/220*common_7 + (27773/165)*common_8 + (1579/44)*common_9),
        ((1120/3)*common_10*common_27*nu + 596*common_10*common_5*nu*speed_squared - 488/45*common_112 - 316/5*common_120 - 1648/45*common_127 - 1960/3*common_149, 1000*common_5*common_62*nu - 8144/15*common_91 - 3216/5*common_94, (3296/15)*common_30 + (784/5)*common_39 - 528*common_49, (304/3)*common_26),
        ((145810/429)*common_110*common_13*radial_speed + (780/7)*common_110*log_infrared*radial_speed + 264*common_110*log_one*nu*radial_speed + (100973/120120)*common_110*nu*radial_speed + (5535/64)*common_119*common_27*nu*radial_speed*speed_squared + (218799/1144)*common_13*common_150*potential + (83937/88)*common_13*common_27*common_62 + (872767/4290)*common_13*common_27*radial_speed*speed_squared + (2435453/5720)*common_13*common_5*common_83 + (106607/1144)*common_13*common_6*common_62*potential + (44037/572)*common_150*common_43*potential + (315/208)*common_150*potential - 5515399/28028*common_151 - 59535/1144*common_152 - 2907/1430*common_153 - 5859/572*common_154 - 177513121/168168*common_155 - 929/20020*common_156 - 15280481/17160*common_157 - 390*common_158 - 3656249/5720*common_159 - 8423/858*common_161 - 493/572*common_162 - 162062/5005*common_163 - 2583/32*common_165 - 24225/572*common_166 - 12389/40*common_167 - 591183/5720*common_168 - 326571/1144*common_169 - 125167/572*common_170 - 175429/1320*common_171 - 12915/64*common_173 - 28537/220*common_174 - 985889/2860*common_175 - 12341/330*common_176 - 20847/286*common_177 + (73559/3432)*common_27*common_43*radial_speed*speed_squared + (2011169/660)*common_27*common_62*nu + (2338333/924)*common_27*common_62 + (1170/7)*common_27*log_infrared*radial_speed*speed_squared + (107379/2860)*common_28*common_43*potential*radial_speed + (218787/5720)*common_28*nu*potential*radial_speed + (280411/1144)*common_43*common_5*common_83 + (30715/572)*common_43*common_6*common_62*potential + (39877/1155)*common_5*common_6*radial_speed + (3819527/2860)*common_5*common_62*nu*speed_squared + (13713/104)*common_83*nu*potential*speed_squared, (5070/7)*common_10*common_29 + (74084729/1261260)*common_110 - 572/7*common_111 + (2091/32)*common_112*common_119 - 16253197/17160*common_112 - 105/176*common_114 - 44667/572*common_115 - 56255/286*common_116 - 7265429/1617*common_117 - 58584/143*common_118 + (16605/64)*common_119*common_121 + (1006963/8580)*common_120 - 8201717/1320*common_121 - 66647/88*common_122 - 2427641/8580*common_123 - 8071/44*common_124 + (1239/52)*common_125 + (604457/17160)*common_126 + (1083211/2860)*common_127 - 3321/64*common_128 + (75261/1144)*common_129 - 22002671/17160*common_130 - 3080631/2860*common_131 - 108083/88*common_132 - 868289/1144*common_133 - 27093/104*common_134 + (626465/1144)*common_135 + (1372555/1144)*common_14*common_32 + (98823/2860)*common_149 + (150022/1287)*common_160 + (901/130)*common_178 + (216711949/280280)*common_179 - 1014/7*common_180 + (1785/11)*common_181 + (8933/286)*common_182 + (6215177/40040)*common_183 + (1090627/5720)*common_184 + (1317/572)*common_185 - 14597039/120120*common_186 - 302681/8580*common_187 - 71013/1144*common_188 + (99221/1320)*common_189 + (3247663/4290)*common_190 - 205221/572*common_191 + (359811/572)*common_32*common_75 + (487423/24024)*common_5*common_6, (127241/5720)*common_101 + (4730279/8580)*common_102 + (119243/286)*common_104 + (114366881/42042)*common_136 + (61470/143)*common_138 + (68999/65)*common_139 + (799753/1144)*common_140 - 1107/16*common_141 + (85583/572)*common_142 + (2872561/8580)*common_143 - 14893/52*common_145 - 50097/1430*common_146 - 3120/7*common_192 - 205245/1144*common_85 + (5061/5720)*common_86 - 51631/1144*common_87 + (6574873/60060)*common_88 - 12199379/60060*common_89 + (273489487/60060)*common_91 - 396*common_92 - 227511/5720*common_93 - 783983/2145*common_94 - 270953/286*common_95 + (3665/572)*common_96 + (98873/1144)*common_97 - 225119/2860*common_98 - 1191075/1144*common_99, common_147 - 22469858/63063*common_27 - 4083/1144*common_28 - 10658419/30030*common_30 + (5301/104)*common_31 + (184/143)*common_33 - 199443/5720*common_34 - 126788/3003*common_35 - 1773/1144*common_37 + (151969/1144)*common_38 + (86936/2145)*common_39 - 17949/5720*common_40 + (1226333/2145)*common_41 - 33655/143*common_42 - 4881631/25740*common_44 + (405617/1144)*common_45 - 421251/20020*common_46 + (144151/1430)*common_48 + (438023/2860)*common_49 - 102631/1144*common_51 + (166999/1430)*common_52 - 2709299/4290*common_53 + (27067/1144)*common_54 - 768617/1144*common_55 + (71688/143)*common_56 - 231665/1144*common_58 - 37527/260*common_59 + (2121763/2860)*common_60 + (180341/2860)*common_61),
        (45*common_10*potential - common_193 - 12*common_5, -90*common_21, common_194, 0),
        (-common_148*common_6 - common_194*common_74 + (253/3)*common_27 - 14/3*common_30 + (35/2)*common_33 + (140/3)*common_35 + 70*common_37 + 6*common_38 - 160/3*common_39 + (105/2)*common_40 - 433*common_46 + 131*common_49, -common_193*radial_speed - 251*common_26 + 709*common_5*radial_speed - 210*common_63 - 75*common_67 + 45*nu*potential*radial_speed*speed_squared, -common_109*speed_squared + 297*common_11 + 18*common_12 - 238*common_5 - 18*common_7 + 140*common_8, 16*common_18 - 146*common_21),
        ((33061/110)*common_10*common_13*common_27 + (633/11)*common_10*common_13*common_6*potential + (269063/165)*common_10*common_27 + (6822/55)*common_10*common_5*nu*speed_squared + (321/11)*common_10*common_6*potential - 8983/33*common_110 - 2165/11*common_112 + (15561/88)*common_113*potential - 7/11*common_116 - 158133/110*common_121 - 2394/11*common_122 - 22816/165*common_123 - 3924/11*common_129 + (5158/11)*common_13*common_32*common_5 + (2555/22)*common_13*common_32*potential*speed_squared + (11966/165)*common_13*common_5*common_6 - 12502/11*common_149 - 27/10*common_178 - 25517/165*common_179 - 20727/44*common_181 - 11865/88*common_182 - 111/10*common_184 - 43203/55*common_186 - 19939/55*common_190 + (60488/165)*common_27*nu*speed_squared + (393/10)*common_28*nu*potential + (4198/11)*common_32*common_5 + (27965/44)*common_32*nu*potential*speed_squared + (1019/55)*common_5*common_6*nu + (3774/55)*common_5*common_6, -7335/22*common_102 - 326342/165*common_136 - 8252/11*common_138 + (25741/22)*common_139 + (7266/11)*common_140 + (1008/11)*common_142 + (6196/55)*common_143 + (9947/11)*common_85 - 984/11*common_86 + (8911/22)*common_87 + (33122/55)*common_88 + (10733/22)*common_89 + (798211/330)*common_91 + (7773/22)*common_93 + (17002/55)*common_94 - 32972/55*common_95 - 23807/44*common_96 - 6846/11*common_99, (5103/11)*common_10*common_13*potential*speed_squared + (2921/11)*common_10*common_5*nu + (13496/55)*common_13*common_27 + (44846/165)*common_27 - 111239/165*common_30 + (4889/11)*common_32*potential - 1592/55*common_35 - 9907/22*common_37 - 15029/110*common_38 - 6502/55*common_39 - 5277/22*common_40 - 10424/11*common_46 - 39/11*common_48 - 5628/11*common_51 - 1191/55*common_52 - 46/55*common_53 + (3013/55)*common_6*potential - 780/11*common_60, -8266/55*common_105 + (9087/55)*common_13*common_5*radial_speed + (1536/11)*common_13*common_62*potential - 33731/55*common_26 + (7154/55)*common_5*radial_speed - 830/11*common_63 - 739/11*common_67 - 2077/55*common_80 + (17933/55)*nu*potential*radial_speed*speed_squared),
        ((3712/5)*common_157 + 5880*common_159 + (4384/45)*common_164 + 1508*common_167 - 7696/3*common_172 - 20356/3*common_195, 3720*common_10*common_27*nu + 8004*common_10*common_5*nu*speed_squared - 448/15*common_112 - 4164/5*common_120 - 1784/5*common_127 - 8960*common_149, 4640*common_5*common_62*nu - 8448/5*common_91 - 14496/5*common_94, (1024/5)*common_30 + (4304/15)*common_39 - 832*common_49),
        ((10701/16)*common_10*common_110*common_119*nu + (1325782/2145)*common_10*common_110*common_43 - 320967/715*common_10*common_110 - 780*common_10*common_111 + (401938/429)*common_10*common_13*common_27*speed_squared + (227751/286)*common_10*common_13*common_28*potential + (94238853/5720)*common_10*common_27*nu*speed_squared + (81429708/5005)*common_10*common_27*speed_squared + (4335/572)*common_10*common_28*potential + (234308/715)*common_10*common_43*common_5*common_6 + (9001071/1430)*common_10*common_5*common_6*nu - 2687273/1144*common_103*common_113 - 7617/286*common_11*common_45 - 369/4*common_110*common_119*common_19 + (1486712/2145)*common_110*common_13*speed_squared - 31666717/8580*common_110*common_57 + (624/7)*common_110*log_infrared*speed_squared + (4656089/8580)*common_110*nu*speed_squared + (2788949/19110)*common_110*speed_squared + (1117683/286)*common_113*common_13*potential*speed_squared - 2292077/572*common_113*common_15 + (717381/286)*common_113*common_43*potential*speed_squared + (783929/104)*common_113*common_5*nu - 3059127/8008*common_113*common_5 - 443583/286*common_113*common_7 + (235179/2288)*common_113*potential*speed_squared - 38745/32*common_117*common_119*common_19 - 123/8*common_119*common_200 + (116235/64)*common_119*common_27*common_32*nu + (5535/64)*common_119*common_27*common_6*nu - 364137/1144*common_12*common_28 - 65711/429*common_13*common_196 + (33171/2860)*common_13*common_28*common_5 + (261387/52)*common_13*common_32*common_5*speed_squared - common_137*common_200 - 10742/165*common_16*common_27 - 22876/11*common_16*common_33 - 68781/143*common_16*common_46 - 1182061/4290*common_160*speed_squared - 21372649/1430*common_183*common_19 + (23653/2574)*common_196*common_43 - 208/7*common_196*log_infrared + (27413343/10010)*common_196*nu + (461973203/630630)*common_196 - 218799/104*common_197*common_50 + (59535/104)*common_197*nu*potential - 44037/52*common_198*common_43 - 3465/208*common_198 + (107379/2860)*common_199*common_43*potential - 591183/5720*common_199*common_50 + (218787/5720)*common_199*nu*potential - 2907/1430*common_199*potential - 396*common_201*log_one - 62745329/20020*common_201 + (19549/10)*common_27*common_32*common_43 + 3510*common_27*common_32*log_infrared - 14714965/616*common_27*common_32 + (44504/165)*common_27*common_43*common_6 + (1170/7)*common_27*common_6*log_infrared - 8948048/8085*common_27*common_6 - 204129/220*common_27*common_9 + (61235/3003)*common_28*common_5 - 473543/1430*common_28*common_8 - 2340*common_29*common_47 - 11111319/440*common_30*common_32 - 497/11*common_32*common_34 - 1921553/220*common_32*common_41 + (1020106/715)*common_32*common_43*common_5*speed_squared + (4060967/10010)*common_32*common_5*speed_squared - 1470*common_32*common_55 + (84105/88)*common_32*common_6*nu*potential - 833801/715*common_44*common_47 - 171511/858*common_45*common_5 - 298409/910*common_46*common_6, (7936989/1144)*common_103*common_83 + (14327591/2860)*common_110*common_25 + (3900/7)*common_111*radial_speed + (16605/16)*common_119*common_157 + (1278747/286)*common_144*common_34 + (226365/26)*common_15*common_83 + (1056762/143)*common_150*common_50 + (63420342/35035)*common_151 - 1849365/1144*common_152 - 31869/1430*common_153 - 289373/1144*common_154 - 1763165242/105105*common_155 + (511419/8008)*common_156 - 52198457/3432*common_157 - 6240*common_158 - 5754711/1144*common_159 + (1825039/286)*common_16*common_63 - 2793383/2145*common_161 + (15850/143)*common_162 + (6469157/10010)*common_163 + 660*common_164*log_one + (90420766/15015)*common_164 - 3321/4*common_165 - 2396343/1144*common_166 - 5855293/1430*common_167 - 401283/286*common_168 - 1973258/143*common_169 - 880929/143*common_170 - 3864989/4290*common_171 + (266946463/5720)*common_172 - 38745/16*common_173 + (138581/260)*common_174 - 1275691/143*common_175 + (382315/286)*common_176 - 4819839/715*common_177 + (2748981/5720)*common_18*common_28 + (9633523/715)*common_195 + (538209/286)*common_202*common_43 + (2835/286)*common_202 - 483627/715*common_21*common_45 + (236732533/6006)*common_27*common_62 + (18720/7)*common_29*common_65 + (43563359/4290)*common_41*common_62 + (2317137/1430)*common_41*common_65 - 114042/715*common_44*common_65 - 389547/10010*common_6*common_64 + (440643/104)*common_7*common_83, (16605/16)*common_10*common_119*common_27*nu + (212514/55)*common_10*common_13*common_5*speed_squared + (25740/7)*common_10*common_27*log_infrared + (262737/55)*common_10*common_43*common_5*speed_squared + (145119/110)*common_10*common_6*nu*potential + 246*common_110*common_119*nu + (163684/495)*common_110*common_43 + (416/7)*common_110*log_infrared - 22039319/24255*common_110 - 505011/110*common_112 + (881475/572)*common_113*nu*potential - 52675/1144*common_114 - 464169/1430*common_115 - 83722/165*common_116 - 472521163/21021*common_117 - 2117877/572*common_118 - 96132347/2860*common_121 - 161497/26*common_122 - 1211447/2145*common_123 - 30387/26*common_124 - 3321/16*common_128 + (991651/1430)*common_13*common_28*potential + (1473185/143)*common_13*common_32*potential*speed_squared - 307829/110*common_130 - 281371/52*common_131 - 63397/4290*common_132 - 1476225/286*common_133 - 33573/11*common_134 - 1341591/572*common_149 - 5148/7*common_180 - 14679/110*common_185 - 364724/385*common_186 - 461707/715*common_187 - 5539/13*common_188 - 3896723/4290*common_189 - 27369/10*common_191 + (60244/195)*common_27*common_43*speed_squared + (32210219/8580)*common_27*nu*speed_squared + (84340297/21021)*common_27*speed_squared + (230249/1430)*common_28*common_43*potential + (4931/130)*common_28*potential + (1148751/286)*common_32*common_43*potential*speed_squared + (1428293/4004)*common_32*common_5 + (366269/1144)*common_32*potential*speed_squared + (330202/715)*common_5*common_6*nu + (268973/1001)*common_5*common_6, -611794/2145*common_102 + (234029/143)*common_13*common_5*radial_speed*speed_squared + (184491/143)*common_13*common_83*potential - 1107/8*common_141 - 48954/55*common_142 - 678981/715*common_145 - 411975/286*common_146 - 4680/7*common_192 + (1820656/2145)*common_27*common_43*radial_speed + (9104153/1430)*common_27*nu*radial_speed + (447288476/105105)*common_27*radial_speed + (958179/715)*common_43*common_5*common_62 + (40122/55)*common_43*common_6*potential*radial_speed + (35661/143)*common_43*common_83*potential + (1467123/715)*common_5*common_62*nu + (42771/110)*common_6*nu*potential*radial_speed + (271551/715)*common_62*nu*potential*speed_squared - 36345/143*common_85 - 2691/55*common_86 - 43717/1430*common_87 - 1381801/10010*common_88 - 1313329/30030*common_89 - 19676/143*common_94 - 5245528/2145*common_95 - 15/572*common_96 - 171338/143*common_99),
    )

def j2(potential, radial_speed, speed_squared, nu, log_one, log_infrared):
    common_0 = nu*potential
    common_1 = (15/7)*common_0
    common_2 = nu*radial_speed
    common_3 = potential*speed_squared
    common_4 = potential**2
    common_5 = speed_squared**2
    common_6 = common_0*speed_squared
    common_7 = common_4*nu
    common_8 = common_5*nu
    common_9 = radial_speed**2
    common_10 = common_9*potential
    common_11 = common_0*common_9
    common_12 = nu**2
    common_13 = common_12*common_3
    common_14 = common_12*common_4
    common_15 = common_12*common_5
    common_16 = common_10*common_12
    common_17 = potential*radial_speed
    common_18 = radial_speed*speed_squared
    common_19 = common_0*radial_speed
    common_20 = common_2*common_4
    common_21 = potential**3
    common_22 = speed_squared**3
    common_23 = common_12*common_22
    common_24 = common_4*speed_squared
    common_25 = nu**3
    common_26 = common_21*common_25
    common_27 = common_21*nu
    common_28 = common_5*potential
    common_29 = radial_speed**4
    common_30 = common_29*potential
    common_31 = common_4*common_9
    common_32 = common_22*nu
    common_33 = common_22*common_25
    common_34 = common_12*common_21
    common_35 = common_15*potential
    common_36 = common_25*common_28
    common_37 = common_14*speed_squared
    common_38 = common_12*common_30
    common_39 = common_0*common_29
    common_40 = common_3*common_9
    common_41 = common_7*common_9
    common_42 = common_25*common_30
    common_43 = common_14*common_9
    common_44 = common_25*common_31
    common_45 = common_24*common_25
    common_46 = common_0*common_5
    common_47 = nu*speed_squared
    common_48 = common_4*common_47
    common_49 = common_25*common_40
    common_50 = common_6*common_9
    common_51 = common_13*common_9
    common_52 = common_21*log_infrared
    common_53 = 22*log_one
    common_54 = math.pi**2
    common_55 = (123/128)*common_54
    common_56 = common_27*common_53 - common_27*common_55 + (214/105)*common_52
    common_57 = common_3*radial_speed
    common_58 = radial_speed**3
    common_59 = common_58*potential
    common_60 = common_4*radial_speed
    common_61 = common_5*radial_speed
    common_62 = common_12*common_59
    common_63 = common_14*radial_speed
    common_64 = common_0*common_58
    common_65 = radial_speed**5
    common_66 = common_65*potential
    common_67 = (214/35)*common_52
    common_68 = common_0*common_65
    common_69 = common_17*common_5
    common_70 = common_18*common_4
    common_71 = common_4*common_58
    common_72 = common_2*common_21
    common_73 = common_47*common_60
    common_74 = common_34*radial_speed
    common_75 = common_14*common_58
    common_76 = common_58*common_7
    common_77 = potential**4
    common_78 = common_77*nu
    common_79 = radial_speed**6
    common_80 = common_0*common_79
    common_81 = common_0*common_22
    common_82 = common_29*common_3
    common_83 = common_21*common_9
    common_84 = common_29*common_6
    common_85 = common_29*common_7
    common_86 = common_4*common_8
    common_87 = common_27*common_9
    common_88 = common_79*potential
    common_89 = common_12*common_88
    common_90 = common_25*common_88
    common_91 = common_33*potential
    common_92 = common_26*speed_squared
    common_93 = common_24*common_9
    common_94 = common_21*common_47
    common_95 = common_54*common_94
    common_96 = common_31*common_47
    common_97 = common_14*common_29
    common_98 = common_15*common_4
    common_99 = common_34*common_9
    common_100 = common_29*common_4
    common_101 = common_100*common_25
    common_102 = common_10*common_15
    common_103 = common_45*common_9
    common_104 = common_12*common_66
    common_105 = common_15*common_17
    common_106 = common_3*common_58
    common_107 = common_21*radial_speed
    common_108 = radial_speed**7
    common_109 = common_77*log_infrared
    common_110 = common_2*common_77
    common_111 = common_25*common_77
    common_112 = common_4*common_5
    common_113 = common_10*common_5
    return (
        (-1, 0),
        (-common_1 + (17/7)*nu*speed_squared - 27/14*potential - 13/28*speed_squared, (5/14)*common_2 - 5/28*radial_speed),
        ((5/252)*common_10 + (241/252)*common_11 + (121/12)*common_13 - 293/84*common_14 - 505/56*common_15 + (335/84)*common_16 - 671/252*common_3 + (43/252)*common_4 - 29/84*common_5 + (1297/126)*common_6 + (1543/126)*common_7 + (11/3)*common_8, -25/14*common_12*common_18 + (173/84)*common_12*potential*radial_speed - 103/63*common_17 - 25/168*common_18 - 337/126*common_19 + (25/24)*nu*radial_speed*speed_squared),
        ((62/7)*common_20, -216/35*common_7 + (4/5)*nu*potential*speed_squared),
        (-23443/7425*common_21 - 709/2464*common_22 - 179615/7392*common_23 - 166931/27720*common_24 - 85031/16632*common_26 - 79067/9240*common_27 - 74083/22176*common_28 + (67/2464)*common_30 + (191/560)*common_31 + (35059/7392)*common_32 + (73517/1848)*common_33 + (120185/5544)*common_34 - 803783/22176*common_35 - 25621/462*common_36 - 9283/264*common_37 - 7813/7392*common_38 - 1667/2464*common_39 - 109/11088*common_40 + (1909/336)*common_41 + (3137/462)*common_42 + (9073/144)*common_43 + (13169/504)*common_44 + (145669/5544)*common_45 + (201835/7392)*common_46 + (273443/5544)*common_48 - 31387/924*common_49 + (2687/1232)*common_50 + (13507/11088)*common_51 + common_56, (17783/1386)*common_12*potential*radial_speed*speed_squared - 49025/7392*common_15*radial_speed - 98485/3696*common_20 + (41215/5544)*common_25*common_4*radial_speed + (32905/3696)*common_25*common_5*radial_speed - 941/66*common_25*common_57 - 967/462*common_25*common_59 + (11815/7392)*common_5*nu*radial_speed - 193/99*common_57 + (3053/1848)*common_58*nu*potential - 1/396*common_59 - 164023/18480*common_60 - 305/2464*common_61 - 14605/5544*common_62 - 194041/11088*common_63 + (6095/1848)*nu*potential*radial_speed*speed_squared),
        (0, -1),
        (-27/28*common_17 - 15/14*common_19, -5/2*common_0 + (39/14)*nu*speed_squared - 7/4*potential - 9/14*speed_squared),
        (-87/56*common_0*common_18 + (99/28)*common_12*potential*radial_speed*speed_squared - 89/9*common_20 - 13/56*common_57 - 5/84*common_59 - 449/126*common_60 - 8/7*common_62 - 379/84*common_63 - 5/168*common_64, -4/21*common_10 + (1115/168)*common_11 + (90/7)*common_13 - 29/6*common_14 - 605/56*common_15 + (109/28)*common_16 - 667/168*common_3 + (275/252)*common_4 - 83/168*common_5 + (1135/168)*common_6 + (2015/126)*common_7 + (113/24)*common_8),
        ((74/35)*common_27 - 186/7*common_41 + (338/35)*common_48, (94/7)*common_20),
        (-795/308*common_0*common_61 + (116401/5544)*common_12*common_4*radial_speed*speed_squared + (7207/462)*common_12*common_5*potential*radial_speed + (565/132)*common_12*common_65*potential - 6511/616*common_13*common_58 + (97973/1260)*common_21*nu*radial_speed + (367141/25200)*common_21*radial_speed + (21299/924)*common_25*common_4*radial_speed*speed_squared + (6595/924)*common_25*common_58*potential*speed_squared - 3215/3696*common_25*common_66 - 55057/3696*common_25*common_69 - 30013/5544*common_25*common_71 - 272/21*common_26*radial_speed + (39185/1848)*common_4*common_58*nu - common_53*common_72 - common_55*common_72 + (1073/616)*common_58*nu*potential*speed_squared + (73/528)*common_58*potential*speed_squared - 335/2464*common_66 - common_67*radial_speed - 205/1232*common_68 - 251/1232*common_69 - 8509/27720*common_70 - 15833/4620*common_71 - 73813/2772*common_73 - 18119/1008*common_74 - 14821/5544*common_75, (653519/831600)*common_21 - 507/1232*common_22 - 7145/231*common_23 - 660211/55440*common_24 - 2503/297*common_26 + (132799/27720)*common_27 - 18545/3696*common_28 + (257/7392)*common_30 + (105547/55440)*common_31 + (23437/3696)*common_32 + (179939/3696)*common_33 + (42257/1232)*common_34 - 18395/924*common_35 - 37927/528*common_36 - 686459/11088*common_37 + (3131/616)*common_38 - 4107/1232*common_39 - 929/3696*common_40 + (131249/3696)*common_41 + (24049/3696)*common_42 + (1060349/11088)*common_43 + (36137/1386)*common_44 + (14851/396)*common_45 + (8915/308)*common_46 + (328931/11088)*common_48 - 32995/924*common_49 + (25007/1848)*common_50 - 673/22*common_51 + common_56),
        (potential, 0),
        ((81/28)*common_10 + (12/7)*common_11 + (19/28)*common_3 - 9/7*common_4 - 6/7*common_6 + (11/7)*common_7, -83/14*common_17 - 1/7*common_19),
        ((155/252)*common_21 + (197/252)*common_24 - 1025/252*common_27 + (19/21)*common_28 + (25/84)*common_30 + (677/126)*common_31 + (185/84)*common_34 + (111/56)*common_35 - 55/21*common_37 + (5/56)*common_38 + (85/42)*common_39 + (29/56)*common_40 + (587/252)*common_41 + (625/84)*common_43 - 173/42*common_46 - 5/252*common_48 - 25/28*common_50 - 27/7*common_51, (67/84)*common_12*common_4*radial_speed + (43/28)*common_12*common_58*potential - 37/28*common_13*radial_speed - 377/18*common_20 - 11/3*common_57 + (43/84)*common_58*potential - 3167/252*common_60 - 1081/84*common_64 + (1607/84)*nu*potential*radial_speed*speed_squared),
        ((302/105)*common_72 - 3042/35*common_73 + (930/7)*common_76, -12/7*common_27 + (184/7)*common_4*nu*speed_squared - 468/7*common_41),
        (-13241/1848*common_101 - 31327/1232*common_102 - 599/33*common_103 + (64177/5544)*common_12*common_21*speed_squared + (661/33)*common_12*common_22*potential + (29255/924)*common_12*common_29*potential*speed_squared + (174569/5544)*common_12*common_4*common_9*speed_squared + (157/231)*common_12*common_77 + (15494/693)*common_21*common_25*common_9 + (1845/128)*common_21*common_54*common_9*nu + (214/7)*common_21*common_9*log_infrared + (14561/528)*common_21*nu*speed_squared + (60169/3300)*common_21*speed_squared + (1297/1848)*common_22*potential + (4625/616)*common_25*common_29*potential*speed_squared + (16883/1848)*common_25*common_4*common_5 + (235/22)*common_25*common_5*common_9*potential + (11183/4158)*common_25*common_77 + (184601/9240)*common_29*common_4 + (24742/3465)*common_4*common_5 + (323/176)*common_5*common_9*nu*potential + (79/77)*common_5*common_9*potential + (287/64)*common_54*common_77*nu - common_67*speed_squared + (428/105)*common_77*log_infrared - 3495893/207900*common_77 - 64271/792*common_78 + (335/352)*common_79*potential - 45/44*common_80 - 9203/1232*common_81 - 10135/7392*common_82 - 567601/5544*common_83 - 515/616*common_84 - 21575/616*common_85 - 3484/693*common_86 - 3885353/11088*common_87 - 10045/528*common_89 - 160/33*common_90 - 4329/616*common_91 - 12589/2772*common_92 - 174199/13860*common_93 - 369/128*common_95 - 74099/2772*common_96 - 9189/308*common_97 - 73421/5544*common_98 - 7342/231*common_99, -16285/1848*common_104 - 68921/924*common_105 - 7109/462*common_106*common_25 + (563561/5544)*common_12*common_4*radial_speed*speed_squared + (39819/616)*common_12*common_58*potential*speed_squared + (1157/154)*common_21*common_25*radial_speed + (1027937/5544)*common_21*nu*radial_speed + (702529/34650)*common_21*radial_speed + (4471/264)*common_25*common_4*common_58 + (17123/1848)*common_25*common_5*potential*radial_speed + (1435/264)*common_25*common_65*potential - 2551/154*common_25*common_70 + (60133/924)*common_4*common_58*nu + (24197/1386)*common_4*nu*radial_speed*speed_squared + (19611/616)*common_5*nu*potential*radial_speed - 428/35*common_52*radial_speed - 369/64*common_54*common_72 - 7813/308*common_58*common_6 + (953/924)*common_58*potential*speed_squared + (590/77)*common_65*nu*potential - 1145/3696*common_66 - 955/308*common_69 - 69809/6930*common_70 - 112943/6930*common_71 - 13411/231*common_74 - 367349/5544*common_75),
        (-3*common_17, potential),
        (-69/7*common_20 + (111/14)*common_4*radial_speed - 405/28*common_59 - 60/7*common_64 + 6*nu*potential*radial_speed*speed_squared + (15/4)*potential*radial_speed*speed_squared, common_1*common_9 + (579/28)*common_10 - 21/4*common_3 + (65/14)*common_4 - common_6 + (12/7)*common_7),
        (-5/8*common_104 - 765/56*common_105 - 235/168*common_106 - 643/126*common_107 + (809/28)*common_12*common_4*radial_speed*speed_squared + (275/14)*common_12*common_58*potential*speed_squared + (8501/252)*common_21*nu*radial_speed + (1139/84)*common_4*radial_speed*speed_squared + (74/7)*common_5*nu*potential*radial_speed + (1055/84)*common_58*nu*potential*speed_squared - 25/12*common_66 - 85/6*common_68 - 47/28*common_69 - 11/42*common_71 - 689/84*common_73 - 1069/42*common_74 - 199/4*common_75 - 289/42*common_76, -1327/126*common_21 - 46/21*common_24 + (1121/252)*common_27 - 58/21*common_28 - 95/42*common_30 - 815/84*common_31 + (47/42)*common_34 + (37/56)*common_35 - 1/14*common_37 - 425/56*common_38 + (5575/84)*common_39 + (731/56)*common_40 + (1337/12)*common_41 + (41/14)*common_43 + (1261/84)*common_46 - 133/6*common_48 - 2713/28*common_50 + (33/7)*common_51),
        ((1300/21)*common_21*nu*speed_squared + (5832/7)*common_4*common_9*nu*speed_squared - 626/105*common_78 - 930*common_85 - 3042/35*common_86 - 1226/7*common_87, (9542/105)*common_72 - 10482/35*common_73 + (3270/7)*common_76),
        (-8740801/11088*common_107*common_47 + (30135/176)*common_108*common_12*potential + (480/11)*common_108*common_25*potential + (405/44)*common_108*nu*potential - 3015/352*common_108*potential - 428/7*common_109*radial_speed - 2583/64*common_110*common_54 - 44*common_110*log_one - 39565/693*common_111*radial_speed - 183971/1848*common_112*common_2 + (1266145/5544)*common_12*common_21*common_58 + (307645/1848)*common_12*common_4*common_5*radial_speed + (258383/1232)*common_12*common_4*common_65 + (937985/3696)*common_12*common_5*common_58*potential + (7900/231)*common_12*common_77*radial_speed - 88645/264*common_13*common_65 - 3/56*common_17*common_22 - 6213/56*common_17*common_23 - 3962209/13860*common_18*common_21 - 129355/924*common_18*common_34 + (2845/33)*common_21*common_25*radial_speed*speed_squared + (5535/128)*common_21*common_54*nu*radial_speed*speed_squared + (8693645/3696)*common_21*common_58*nu + (21323891/27720)*common_21*common_58 + (642/7)*common_21*log_infrared*radial_speed*speed_squared + (2377/56)*common_22*common_25*potential*radial_speed + (2921/112)*common_22*nu*potential*radial_speed - 7185/88*common_25*common_3*common_65 + (60259/616)*common_25*common_4*common_58*speed_squared + (28333/616)*common_25*common_4*common_65 - 50923/616*common_25*common_5*common_60 - 131147/693*common_26*common_58 - 12915/128*common_27*common_54*common_58 - 19615/1848*common_28*common_58 - 1800/77*common_36*common_58 - 161849/462*common_37*common_58 + (17911/264)*common_4*common_58*nu*speed_squared + (301765/1848)*common_4*common_58*speed_squared + (352651/1232)*common_4*common_65*nu - 964247/6160*common_4*common_65 - 15425/1232*common_46*common_58 - 36733/924*common_5*common_60 - 214*common_52*common_58 - 25/88*common_6*common_65 + (16165/1056)*common_65*potential*speed_squared + (24064457/27720)*common_77*nu*radial_speed + (385093/1540)*common_77*radial_speed, (471199/3696)*common_100 - 18035/154*common_101 + (483271/1232)*common_102 + (34913/308)*common_103 + (1712/105)*common_109 - 209161/1232*common_11*common_5 - 13423/4158*common_111 - 4097/308*common_112*common_25 + (1181/308)*common_112 - 38979/616*common_113*common_25 + (2067/154)*common_113 + (64157/2772)*common_12*common_77 - 11075/33*common_13*common_29 + (1610611/34650)*common_21*speed_squared - 403/168*common_22*potential - 4583/84*common_23*potential + (8595/77)*common_25*common_82 - 6161/198*common_26*common_9 - 116957/2772*common_34*speed_squared - 29185/44*common_37*common_9 + (642/7)*common_52*common_9 - 642/35*common_52*speed_squared + (41/4)*common_54*common_78 + (5535/128)*common_54*common_87 - 7069217/207900*common_77 - 151071/616*common_78 - 2405/44*common_80 + (2729/112)*common_81 - 59705/7392*common_82 - 8239229/27720*common_83 + (101215/616)*common_84 - 730725/1232*common_85 - 18019/1848*common_86 - 5538187/3696*common_87 + (3295/1056)*common_88 + (22525/528)*common_89 - 3775/88*common_90 + (47/21)*common_91 + (10247/693)*common_92 - 73361/924*common_93 + (3158255/11088)*common_94 - 1107/128*common_95 + (294647/924)*common_96 + (244781/528)*common_97 + (77285/924)*common_98 + (2380943/5544)*common_99),
    )
