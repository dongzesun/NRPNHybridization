"""Generated higher-rank SH source moments and their derivatives."""

import math

INVENTORY = {'i4': (4, 4, 0, 5, (0, 2, 4, 5)), 'i5': (5, 5, 0, 6, (0, 2, 4)), 'i6': (6, 6, 0, 7, (0, 2)), 'i7': (7, 7, 0, 8, (0, 2)), 'j3': (3, 2, 1, 4, (0, 2, 4, 5)), 'j4': (4, 3, 1, 5, (0, 2, 4)), 'j5': (5, 4, 1, 6, (0, 2)), 'j6': (6, 5, 1, 7, (0, 2))}

TAIL_ORDERS = {'i4': 7, 'j3': 6}

def i4(potential, radial_speed, speed_squared, nu, log_one, log_infrared):
    common_0 = nu*speed_squared
    common_1 = nu**2
    common_2 = common_1*potential
    common_3 = potential**2
    common_4 = speed_squared**2
    common_5 = common_0*potential
    common_6 = common_3*nu
    common_7 = common_4*nu
    common_8 = radial_speed**2
    common_9 = common_8*potential
    common_10 = nu**3
    common_11 = common_10*common_3
    common_12 = common_10*common_4
    common_13 = common_2*common_8
    common_14 = potential*radial_speed
    common_15 = radial_speed*speed_squared
    common_16 = common_2*radial_speed
    common_17 = common_1*speed_squared
    common_18 = common_17*radial_speed
    common_19 = common_10*common_14
    common_20 = nu*potential
    common_21 = (28/13)*common_8
    common_22 = common_10*potential
    common_23 = common_10*speed_squared
    common_24 = common_1*common_8
    common_25 = (2128/715)*radial_speed
    common_26 = common_6*radial_speed
    common_27 = common_1*common_3
    common_28 = common_27*radial_speed
    common_29 = 12*nu
    common_30 = common_20*radial_speed
    common_31 = common_1*radial_speed
    common_32 = (156/11)*common_1 - 156/11*nu + 156/55
    common_33 = radial_speed**3
    common_34 = common_0*common_14
    common_35 = common_20*common_33
    common_36 = common_15*common_2
    common_37 = common_22*common_33
    common_38 = common_11*radial_speed
    common_39 = common_0*radial_speed
    common_40 = common_10*common_15
    common_41 = (1176/715)*common_8
    common_42 = nu*radial_speed
    common_43 = common_10*radial_speed
    common_44 = potential**3
    common_45 = common_44*nu
    common_46 = common_0*common_3
    common_47 = common_1*common_44
    common_48 = common_6*common_8
    common_49 = common_17*common_3
    common_50 = common_27*common_8
    common_51 = 36*nu
    common_52 = potential*speed_squared
    common_53 = common_20*common_8
    common_54 = common_2*speed_squared
    common_55 = radial_speed**4
    common_56 = common_55*potential
    common_57 = common_4*potential
    common_58 = common_9*speed_squared
    common_59 = common_3*common_8
    common_60 = common_2*common_55
    common_61 = common_2*common_4
    common_62 = common_22*common_55
    common_63 = common_11*speed_squared
    common_64 = common_23*common_9
    common_65 = (9408/715)*common_8
    common_66 = common_33*common_6
    common_67 = common_31*common_44
    common_68 = common_3*radial_speed
    common_69 = common_17*common_68
    common_70 = (384/5)*common_1*potential*speed_squared - 22848/55*common_27 + (6912/55)*common_3*nu - 96/5*common_5
    common_71 = 24 - 72*nu
    common_72 = common_2*common_33
    common_73 = radial_speed**5
    common_74 = common_73*potential
    common_75 = common_44*radial_speed
    common_76 = common_20*common_73
    common_77 = common_42*common_44
    common_78 = common_14*common_4
    common_79 = common_33*common_52
    common_80 = common_15*common_3
    common_81 = common_3*common_33
    common_82 = common_30*common_4
    common_83 = common_33*potential
    common_84 = common_0*common_83
    common_85 = common_0*common_68
    common_86 = common_2*common_73
    common_87 = common_22*common_73
    common_88 = common_43*common_44
    common_89 = common_16*common_4
    common_90 = common_33*common_54
    common_91 = common_12*common_14
    common_92 = common_10*common_52
    common_93 = common_33*common_92
    common_94 = common_11*common_15
    common_95 = common_27*common_33
    common_96 = common_11*common_33
    common_97 = common_3*speed_squared
    common_98 = common_20*common_55
    common_99 = common_20*common_4
    common_100 = common_10*common_44
    common_101 = common_0*common_9
    common_102 = common_12*potential
    common_103 = common_13*speed_squared
    common_104 = common_11*common_8
    common_105 = common_52*radial_speed
    common_106 = common_43*common_52
    common_107 = potential**4
    common_108 = common_1*common_107
    common_109 = common_45*common_8
    common_110 = common_17*common_44
    common_111 = common_0*common_59
    common_112 = common_27*common_55
    common_113 = common_27*common_4
    common_114 = common_107*nu
    common_115 = radial_speed**6
    common_116 = common_115*potential
    common_117 = speed_squared**3
    common_118 = common_117*common_20
    common_119 = common_0*common_44
    common_120 = common_10*common_107
    common_121 = common_3*common_55
    common_122 = common_3*common_4
    common_123 = common_44*common_8
    common_124 = common_0*common_56
    common_125 = common_115*common_2
    common_126 = common_115*common_22
    common_127 = common_117*common_22
    common_128 = common_23*common_44
    common_129 = common_24*common_44
    common_130 = common_13*common_4
    common_131 = common_63*common_8
    common_132 = common_107*common_42
    common_133 = common_0*common_75
    common_134 = common_6*common_73
    common_135 = common_33*common_45
    common_136 = common_107*common_31
    common_137 = common_68*common_7
    common_138 = common_0*common_81
    common_139 = common_17*common_75
    common_140 = common_27*common_73
    common_141 = common_33*common_47
    common_142 = common_28*common_4
    common_143 = common_17*common_81
    common_144 = radial_speed**7
    common_145 = potential**5
    return (
        (1 - 3*nu, 0, 0, 0, 0),
        (-147/22*common_0 + (279/22)*common_1*speed_squared - 105/11*common_2 + (61/11)*nu*potential - 10/11*potential + (103/110)*speed_squared, -72/11*common_1*radial_speed + (72/11)*nu*radial_speed - 72/55*radial_speed, (78/11)*common_1 - 78/11*nu + 78/55, 0, 0),
        ((7961/143)*common_1*common_3 + (112357/2860)*common_1*common_4 + (15124/715)*common_1*potential*speed_squared + (9617/715)*common_10*common_8*potential + (46934/715)*common_10*potential*speed_squared - 5829/286*common_11 - 325687/5720*common_12 - 10557/1430*common_13 - 15549/10010*common_3 + (3649/5720)*common_4 - 152489/7150*common_5 - 9457/715*common_6 - 50191/5720*common_7 + (12619/7150)*common_8*nu*potential - 659/3575*common_9 + (11049/3575)*potential*speed_squared, (25796/715)*common_10*radial_speed*speed_squared - 15463/3575*common_14 - 476/715*common_15 - 25606/715*common_16 - 23512/715*common_18 - 18839/715*common_19 + (98374/3575)*nu*potential*radial_speed + (1228/143)*nu*radial_speed*speed_squared, -6913/715*common_0 + common_10*common_21 + (25994/715)*common_17 - 17202/715*common_2 - 72474/3575*common_20 + common_21*nu + (27568/715)*common_22 - 28207/715*common_23 - 56/13*common_24 - 4/13*common_8 + (27818/3575)*potential + (553/715)*speed_squared, (4256/715)*common_1*radial_speed - common_10*common_25 - common_25*nu + (304/715)*radial_speed, (994/715)*common_1 - 497/715*common_10 - 497/715*nu + 71/715),
        (-100/11*common_26 + (300/11)*common_28, (64/5)*common_1*potential*speed_squared - 3808/55*common_27 + (1152/55)*common_3*nu - 16/5*common_5, 0, 0, 0),
        (0, 4 - common_29, 0, 0, 0),
        ((19/55)*common_14 - 102/11*common_16 + (14/11)*common_30, -222/11*common_0 + (486/11)*common_1*speed_squared - 504/11*common_2 + (328/11)*nu*potential - 284/55*potential + (134/55)*speed_squared, -60/11*common_31 + (60/11)*nu*radial_speed - 12/11*radial_speed, common_32, 0),
        ((8523/715)*common_1*common_33*potential + (13085/286)*common_10*potential*radial_speed*speed_squared - 5925/143*common_28 + (42596/3575)*common_3*nu*radial_speed + (96946/25025)*common_3*radial_speed + (1977/3575)*common_33*potential - 229/325*common_34 - 1446/325*common_35 - 21706/715*common_36 - 17337/1430*common_37 - 38208/715*common_38 + (11311/7150)*potential*radial_speed*speed_squared, (46030/143)*common_1*common_3 + (17769/143)*common_1*common_4 + (17334/715)*common_1*potential*speed_squared + (5909/143)*common_10*common_8*potential + (36627/143)*common_10*potential*speed_squared - 83667/715*common_11 - 54819/286*common_12 - 6674/715*common_13 - 283821/25025*common_3 + (2697/1430)*common_4 - 157254/3575*common_5 - 240846/3575*common_6 - 37911/1430*common_7 + (5834/325)*common_8*nu*potential - 17653/3575*common_9 + (20123/3575)*potential*speed_squared, (14099/3575)*common_14 - 762/715*common_15 - 91776/715*common_16 - 24708/715*common_18 + (10489/715)*common_19 + (85218/3575)*common_30 + (7674/715)*common_39 + (24054/715)*common_40, -15954/715*common_0 - common_10*common_41 + (56244/715)*common_17 - 3876/65*common_2 - 124368/3575*common_20 + (59252/715)*common_22 - 5322/65*common_23 + (2352/715)*common_24 - common_41*nu + (168/715)*common_8 + (52696/3575)*potential + (282/143)*speed_squared, (8232/715)*common_1*radial_speed - 4116/715*common_42 - 4116/715*common_43 + (588/715)*radial_speed),
        (-652/55*common_45 - 324/55*common_46 + (2308/55)*common_47 + (300/11)*common_48 + (796/55)*common_49 - 900/11*common_50, -560/11*common_26 + (1680/11)*common_28, (192/5)*common_1*potential*speed_squared - 11424/55*common_27 + (3456/55)*common_3*nu - 48/5*common_5, 0, 0),
        (common_29*potential - 4*potential, 0, 12 - common_51, 0, 0),
        ((108/11)*common_13 + (342/11)*common_27 + (229/11)*common_3 + (236/11)*common_5 - 67/11*common_52 + (24/11)*common_53 - 192/11*common_54 - 782/11*common_6 - 57/55*common_9, -492/11*common_16 - 564/11*common_30 + (1092/55)*potential*radial_speed, -606/11*common_0 + (1398/11)*common_1*speed_squared - 1920/11*common_2 + (1392/11)*nu*potential - 252/11*potential + (342/55)*speed_squared, (192/11)*common_31 - 192/11*common_42 + (192/55)*radial_speed, common_32),
        ((13272/715)*common_1*common_8*potential*speed_squared + (8353/286)*common_10*common_3*common_8 + (8134/143)*common_10*common_4*potential + (663/11)*common_10*common_44 + (159581/7150)*common_3*common_8*nu + (561903/50050)*common_3*speed_squared + (190757/7150)*common_4*nu*potential + (2103/11)*common_44*nu - 3681/77*common_44 - 35336/3575*common_46 - 1954/11*common_47 - 16543/715*common_49 - 60758/715*common_50 + (1917/130)*common_55*nu*potential - 1977/715*common_56 - 9797/3575*common_57 - 22071/7150*common_58 - 150127/25025*common_59 - 2088/143*common_60 - 61671/715*common_61 - 984/143*common_62 - 111091/1430*common_63 - 2298/65*common_64 + (37254/3575)*common_8*nu*potential*speed_squared, (52554/715)*common_1*common_33*potential + (134047/715)*common_10*potential*radial_speed*speed_squared - 246398/715*common_28 + (104362/715)*common_3*nu*radial_speed + (55148/5005)*common_3*radial_speed + (60867/3575)*common_33*potential - 156474/3575*common_34 - 306066/3575*common_35 - 1728/143*common_36 - 14499/715*common_37 - 221236/715*common_38 + (42743/3575)*potential*radial_speed*speed_squared, (963174/715)*common_1*common_3 + (241827/715)*common_1*common_4 + (318/55)*common_1*common_8*potential + (7442/55)*common_10*common_8*potential + (57446/65)*common_10*potential*speed_squared - 386206/715*common_11 - 774177/1430*common_12 - 1314212/25025*common_3 + (597/130)*common_4 - 132564/3575*common_5 - 130758/715*common_54 - 912012/3575*common_6 - 19677/286*common_7 + (26568/275)*common_8*nu*potential - 8366/275*common_9 + (30608/3575)*potential*speed_squared, (175034/3575)*common_14 + (1632/715)*common_15 - 193788/715*common_16 + (67776/715)*common_18 + (31058/143)*common_19 - 369552/3575*common_30 - 18912/715*common_39 - 71328/715*common_40, -4014/143*common_0 - common_10*common_65 + (64476/715)*common_17 - 50868/715*common_2 - 103788/3575*common_20 + (63368/715)*common_22 - 62658/715*common_23 + (18816/715)*common_24 - common_65*nu + (1344/715)*common_8 + (49756/3575)*potential + (1998/715)*speed_squared),
        ((4500/11)*common_1*common_3*common_33 + (5028/55)*common_3*nu*radial_speed*speed_squared + (15152/165)*common_44*nu*radial_speed - 1500/11*common_66 - 16208/55*common_67 - 14556/55*common_69, -7776/55*common_45 - 3392/55*common_46 + (26848/55)*common_47 + (2880/11)*common_48 + (8416/55)*common_49 - 8640/11*common_50, -2160/11*common_26 + (6480/11)*common_28, common_70, 0),
        (-common_14*common_51 + 12*potential*radial_speed, 120*common_20 - 40*potential, 0, common_71, 0),
        (72*common_1*potential*radial_speed*speed_squared - 708/11*common_28 + (3172/11)*common_3*nu*radial_speed + (57/11)*common_33*potential - 60*common_34 - 120/11*common_35 - 4888/55*common_68 - 540/11*common_72 + (81/5)*potential*radial_speed*speed_squared, (720/11)*common_13 + (4116/11)*common_27 + (11288/55)*common_3 + (1592/11)*common_5 - 2252/55*common_52 + (2184/11)*common_53 - 1680/11*common_54 - 7988/11*common_6 - 3504/55*common_9, -1344/11*common_16 - 4992/11*common_30 + (8556/55)*potential*radial_speed, -1404/11*common_0 + (2988/11)*common_1*speed_squared - 4656/11*common_2 + (3600/11)*nu*potential - 3336/55*potential + (876/55)*speed_squared, (504/11)*common_31 - 504/11*common_42 + (504/55)*radial_speed),
        (-5397/55*common_66 + (171746/143)*common_67 + (63716/715)*common_69 + (13839/715)*common_74 + (857828/3575)*common_75 - 13419/130*common_76 - 4045254/3575*common_77 + (1464/715)*common_78 + (1251/286)*common_79 - 2653988/25025*common_80 + (10746/385)*common_81 - 84651/1430*common_82 + (984/143)*common_84 + (656317/3575)*common_85 + (14616/143)*common_86 + (6888/143)*common_87 - 10713/715*common_88 + (211557/715)*common_89 - 21624/143*common_90 - 172566/715*common_91 + (21342/143)*common_93 + (196396/715)*common_94 + (5736/11)*common_95 - 186*common_96, (602868/715)*common_100 + (30606/715)*common_101 + (319764/715)*common_102 - 43224/143*common_103 + (30595/143)*common_104 - 948816/2275*common_44 + (7186924/3575)*common_45 - 569952/3575*common_46 - 381204/143*common_47 - 59482/143*common_48 + (6726/715)*common_49 - 110304/143*common_50 - 13755/143*common_56 - 14748/715*common_57 + (186/65)*common_58 + (623723/5005)*common_59 - 22296/143*common_60 - 461298/715*common_61 - 47352/143*common_62 - 530097/715*common_63 + (6948/143)*common_64 + (2566453/25025)*common_97 + (63213/143)*common_98 + (144017/715)*common_99, (12747/715)*common_105 + (47721/143)*common_106 + (922914/715)*common_26 - 312054/143*common_28 - 29148/143*common_34 - 8064/13*common_35 + (97842/143)*common_36 + (29895/143)*common_37 - 134198/143*common_38 - 524074/5005*common_68 + (30222/143)*common_72 + (20355/143)*common_83, (2429688/715)*common_1*common_3 + (10026/13)*common_1*common_4 + (3840/13)*common_1*common_8*potential + (39642/143)*common_10*common_8*potential + (1581942/715)*common_10*potential*speed_squared - 1075094/715*common_11 - 169101/143*common_12 - 752342/5005*common_3 + (8199/715)*common_4 - 6528/715*common_5 - 674904/715*common_54 - 398808/715*common_6 - 117297/715*common_7 + (83784/715)*common_8*nu*potential - 47214/715*common_9 + (27018/715)*potential*speed_squared, (70722/715)*common_14 + (756/65)*common_15 - 33576/143*common_16 + (46872/143)*common_18 + (56334/143)*common_19 - 214176/715*common_30 - 77868/715*common_39 - 43092/143*common_40),
        ((28056/11)*common_1*common_3*common_8*speed_squared + (20412/11)*common_1*common_44*common_8 + (8176/165)*common_107*nu - 2128/11*common_108 - 6452/11*common_109 - 10068/55*common_110 - 9528/11*common_111 - 31500/11*common_112 - 14556/55*common_113 + (5028/55)*common_3*common_4*nu + (10500/11)*common_3*common_55*nu + (10244/165)*common_44*nu*speed_squared, (61200/11)*common_1*common_3*common_33 + (65424/55)*common_3*nu*radial_speed*speed_squared + (188576/165)*common_44*nu*radial_speed - 20400/11*common_66 - 40672/11*common_67 - 37776/11*common_69, -720*common_45 - 3984/11*common_46 + (12336/5)*common_47 + (15120/11)*common_48 + (51312/55)*common_49 - 45360/11*common_50, -480*common_26 + 1440*common_28, common_70),
        (28*common_3 - common_51*common_52 - 84*common_6 + 180*common_8*nu*potential - 60*common_9 + 12*potential*speed_squared, -504*common_30 + 168*potential*radial_speed, 576*common_20 - 192*potential, 0, common_71),
        (522*common_1*common_3*common_8 + 72*common_1*common_4*potential + (3780/11)*common_1*common_55*potential - 5580/11*common_103 + 597*common_3*common_8 + (2240/11)*common_3*nu*speed_squared + (81/5)*common_4*potential + (7896/11)*common_44*nu - 2512/11*common_44 - 1560/11*common_47 - 1986*common_48 - 2592/11*common_49 + (840/11)*common_55*nu*potential - 399/11*common_56 - 720/11*common_58 + (2940/11)*common_8*nu*potential*speed_squared - 1987/55*common_97 - 60*common_99, (9648/11)*common_1*potential*radial_speed*speed_squared - 17328/11*common_28 + (53232/11)*common_3*nu*radial_speed + (3732/11)*common_33*potential - 3048/11*common_34 - 11400/11*common_35 - 79104/55*common_68 - 5760/11*common_72 + (3312/55)*potential*radial_speed*speed_squared, (2628/11)*common_13 + (22908/11)*common_27 + (51156/55)*common_3 + (3996/11)*common_5 - 4788/55*common_52 + (22716/11)*common_53 - 8220/11*common_54 - 37692/11*common_6 - 7236/11*common_9, -1272/11*common_16 - 19848/11*common_30 + (6504/11)*potential*radial_speed, -1908/11*common_0 + (3492/11)*common_1*speed_squared - 5160/11*common_2 + (4104/11)*nu*potential - 768/11*potential + (276/11)*speed_squared),
        ((197/143)*common_1*common_107 + (211557/715)*common_1*common_117*potential + (1468266/715)*common_1*common_3*common_8*speed_squared + (1826598/715)*common_1*common_44*speed_squared + (224448/143)*common_1*common_55*potential*speed_squared + (462598/715)*common_10*common_3*common_4 + (642963/286)*common_10*common_3*common_55 + (236592/143)*common_10*common_4*common_8*potential + (664514/715)*common_10*common_44*common_8 + (4477296/5005)*common_107 - 15096519/7150*common_111 - 696159/143*common_112 - 70823/715*common_113 - 1819199/715*common_114 + (120771/130)*common_115*nu*potential - 124551/715*common_116 + (1464/715)*common_117*potential - 84651/1430*common_118 - 3430399/3575*common_119 - 29823/143*common_120 - 853089/5005*common_121 - 1569793/25025*common_122 - 7773683/3575*common_123 - 161385/286*common_124 - 131544/143*common_125 - 61992/143*common_126 - 172566/715*common_127 - 577132/715*common_128 - 5573372/715*common_129 - 276429/143*common_130 - 2182632/715*common_131 + (821309/7150)*common_3*common_4*nu + (85563/130)*common_3*common_55*nu + (45347571/50050)*common_3*common_8*speed_squared + (90555/286)*common_4*common_8*nu*potential + (75/26)*common_4*common_8*potential + (32945189/3575)*common_44*common_8*nu + (100473/3575)*common_44*speed_squared - 114954/143*common_55*common_92 + (18921/286)*common_55*potential*speed_squared, (3280116/715)*common_66 + (16125008/715)*common_67 - 147288/715*common_69 + (536781/715)*common_74 + (102812828/25025)*common_75 - 2507673/715*common_76 - 69442712/3575*common_77 + (54192/715)*common_78 - 54564/143*common_79 - 2587204/2275*common_80 - 6537708/5005*common_81 - 540141/715*common_82 + (20562/13)*common_84 + (5914956/3575)*common_85 + (214536/143)*common_86 + (359016/143)*common_87 - 2438288/715*common_88 + (1797882/715)*common_89 + (40440/143)*common_90 - 1580076/715*common_91 - 138780/143*common_93 + (3111536/715)*common_94 + (909768/143)*common_95 - 371736/143*common_96, (3714316/715)*common_100 - 518616/715*common_101 + (1368672/715)*common_102 - 615438/143*common_103 + (68996/143)*common_104 - 36264268/25025*common_44 + (31641312/3575)*common_45 - 1333986/3575*common_46 - 11313864/715*common_47 - 1364904/143*common_48 + (106878/143)*common_49 + (45594/143)*common_50 - 143040/143*common_56 - 90258/715*common_57 + (273222/715)*common_58 + (12707098/5005)*common_59 - 9288/13*common_60 - 1641054/715*common_61 - 465276/143*common_62 - 2785672/715*common_63 + (29772/13)*common_64 + (7239574/25025)*common_97 + (613854/143)*common_98 + (11910/13)*common_99, (3095424/715)*common_1*potential*radial_speed*speed_squared + (189372/143)*common_10*common_33*potential - 273972/715*common_106 - 636888/65*common_28 + (4511112/715)*common_3*nu*radial_speed + (345192/715)*common_33*potential - 65736/65*common_34 - 1296576/715*common_35 - 406216/715*common_38 - 4233016/5005*common_68 - 52236/143*common_72 + (1524/65)*potential*radial_speed*speed_squared, (2663088/715)*common_1*common_3 + (157158/143)*common_1*common_4 + (168804/143)*common_1*common_8*potential + (13800/143)*common_10*common_8*potential + (180072/65)*common_10*potential*speed_squared - 1291244/715*common_11 - 212193/143*common_12 - 212788/1001*common_3 + (3303/143)*common_4 - 129732/715*common_5 - 257196/715*common_53 - 1011624/715*common_54 - 302568/715*common_6 - 39033/143*common_7 - 3528/715*common_9 + (82872/715)*potential*speed_squared),
        (-41376/55*common_132 - 83712/55*common_133 - 94500/11*common_134 + (49072/11)*common_135 + (145248/55)*common_136 - 24084/11*common_137 + (108696/11)*common_138 + (256416/55)*common_139 + (283500/11)*common_140 - 154608/11*common_141 + (70668/11)*common_142 - 322392/11*common_143, 44064*common_1*common_3*common_8*speed_squared + 28032*common_1*common_44*common_8 + (113408/165)*common_107*nu - 142272/55*common_108 - 8768*common_109 - 167936/55*common_110 - 14976*common_111 - 50400*common_112 - 22464/5*common_113 + (7776/5)*common_3*common_4*nu + 16800*common_3*common_55*nu + (32672/33)*common_44*nu*speed_squared, (410400/11)*common_1*common_3*common_33 + (85248/11)*common_3*nu*radial_speed*speed_squared + (70144/11)*common_44*nu*radial_speed - 136800/11*common_66 - 229440/11*common_67 - 246240/11*common_69, -99456/55*common_45 - 68352/55*common_46 + (340608/55)*common_47 + (46080/11)*common_48 + (183936/55)*common_49 - 138240/11*common_50, -6240/11*common_26 + (18720/11)*common_28),
        (-180*common_105 + 720*common_26 + 540*common_34 - 1260*common_35 - 240*common_68 + 420*common_83, 328*common_3 - 648*common_5 - 984*common_6 + 3240*common_8*nu*potential - 1080*common_9 + 216*potential*speed_squared, -3240*common_30 + 1080*potential*radial_speed, 1440*common_20 - 480*potential, 0),
        ((196596/11)*common_66 + (17196/11)*common_67 + (34092/11)*common_69 + (3591/11)*common_74 + (121168/55)*common_75 - 7560/11*common_76 - 77212/11*common_77 - 2331/11*common_78 + (3444/11)*common_79 + (79236/55)*common_80 - 58512/11*common_81 + (9180/11)*common_82 - 17220/11*common_84 - 61020/11*common_85 - 34020/11*common_86 - 15120/11*common_89 + (54180/11)*common_90 - 55620/11*common_95, (152532/11)*common_1*common_3*common_8 + (12816/11)*common_1*common_4*potential + 5040*common_1*common_55*potential - 7200/11*common_101 - 87840/11*common_103 + (137928/11)*common_3*common_8 + (57248/11)*common_3*nu*speed_squared + (6876/55)*common_4*potential + (77496/11)*common_44*nu - 120968/55*common_44 - 20472/11*common_47 - 467220/11*common_48 - 42288/11*common_49 + 7560*common_55*nu*potential - 2520*common_56 + (5004/11)*common_8*potential*speed_squared - 68908/55*common_97 - 5688/11*common_99, (58860/11)*common_1*potential*radial_speed*speed_squared - 9612/11*common_105 - 153960/11*common_28 + (375720/11)*common_3*nu*radial_speed + (47376/11)*common_33*potential - 147780/11*common_35 - 108408/11*common_68 - 30420/11*common_72 + (24300/11)*nu*potential*radial_speed*speed_squared, (4320/11)*common_13 + (61392/11)*common_27 + (106272/55)*common_3 - 384*common_5 + (1104/5)*common_52 + (106560/11)*common_53 - 2016*common_54 - 82512/11*common_6 - 33984/11*common_9, (10944/11)*common_14 + (3240/11)*common_16 - 34920/11*common_30),
        ((15233211/1430)*common_0*common_74 + (8326748/715)*common_100*common_15 - 10326942/715*common_100*common_33 + (2022804/715)*common_107*common_43 - 222378488/25025*common_107*radial_speed - 3728376/143*common_11*common_73 - 639/143*common_117*common_14 - 764415/143*common_117*common_16 + (645750/143)*common_117*common_19 + (265761/286)*common_117*common_30 - 128412/11*common_12*common_68 - 192360/13*common_12*common_83 + (100495402/3575)*common_132 + (103257276/3575)*common_133 - 787659/143*common_134 - 365501172/3575*common_135 - 1264990/143*common_136 - 1529007/275*common_137 + (82391466/3575)*common_138 - 30050802/715*common_139 + (7679832/143)*common_140 + (57190566/715)*common_141 + (338292/55)*common_142 - 30671748/715*common_143 + (131544/13)*common_144*common_2 - 1328481/130*common_144*common_20 + (61992/13)*common_144*common_22 + (124551/65)*common_144*potential - 11406074/2275*common_15*common_44 + (1287522/1001)*common_3*common_73 + (629607378/25025)*common_33*common_44 + (27161766/715)*common_33*common_63 - 223605864/25025*common_33*common_97 - 1279425/286*common_35*common_4 + (4542438/1925)*common_4*common_68 + (2832795/143)*common_4*common_72 + (69909/286)*common_4*common_83 - 2346057/1430*common_52*common_73 - 2809296/143*common_54*common_73 + (662634/143)*common_73*common_92, (191314/143)*common_1*common_107 + (528822/143)*common_1*common_117*potential + (27472944/715)*common_1*common_3*common_8*speed_squared + (27175088/715)*common_1*common_44*speed_squared + (1687392/143)*common_1*common_55*potential*speed_squared + (1396476/143)*common_10*common_3*common_4 + (5367357/143)*common_10*common_3*common_55 + (2110104/143)*common_10*common_4*common_8*potential + (6085364/143)*common_10*common_44*common_8 + (2306724/143)*common_10*common_55*potential*speed_squared + (150148788/25025)*common_107 - 29870964/3575*common_111 - 10767318/143*common_112 - 1476674/715*common_113 - 60058502/3575*common_114 + (1940463/55)*common_115*nu*potential - 409941/55*common_116 + (60048/715)*common_117*potential - 709443/715*common_118 - 83143872/3575*common_119 - 1883364/715*common_120 - 2496996/2275*common_122 - 112978352/2275*common_123 - 401247/13*common_124 - 189000/11*common_125 - 267624/11*common_126 - 454068/143*common_127 - 9651768/715*common_128 - 146340808/715*common_129 - 2782278/143*common_130 - 7833456/143*common_131 + (5639059/3575)*common_3*common_4*nu + (71103213/5005)*common_3*common_55 + (193869006/25025)*common_3*common_8*speed_squared + (1399797/143)*common_4*common_8*nu*potential - 216234/143*common_4*common_9 + (759645628/3575)*common_44*common_8*nu + (97646048/25025)*common_44*speed_squared - 34998924/715*common_55*common_6 + (86961/13)*common_55*potential*speed_squared, (1740240/143)*common_1*common_3*common_33 + (2789838/143)*common_1*common_33*potential*speed_squared + (4162428/715)*common_1*common_4*potential*radial_speed + (125666208/715)*common_1*common_44*radial_speed + (1358784/143)*common_1*common_73*potential + (1698948/55)*common_10*common_3*radial_speed*speed_squared + (4333980/143)*common_10*common_73*potential + (72364896/715)*common_3*common_33*nu + (3652578/143)*common_33*nu*potential*speed_squared + (979794/715)*common_4*potential*radial_speed + (645567656/25025)*common_44*radial_speed - 23400/11*common_69 + (6616743/715)*common_73*potential - 29007909/715*common_76 - 460607244/3575*common_77 - 91734/13*common_79 - 2883096/1925*common_80 - 139582332/5005*common_81 - 84051/13*common_82 - 2332116/275*common_85 - 30987212/715*common_88 - 506484/65*common_91 - 3914904/143*common_93 - 1951200/143*common_96, (1796028/55)*common_1*common_3*common_8 + (211284/143)*common_1*common_55*potential + (3159792/715)*common_10*common_4*potential + (10151264/715)*common_10*common_44 + (9534816/715)*common_10*common_8*potential*speed_squared - 212832/65*common_101 - 18306792/715*common_103 - 164664/55*common_104 + (6405492/385)*common_3*common_8 + (1299516/275)*common_3*nu*speed_squared + (1442544/715)*common_4*nu*potential + (4213584/325)*common_44*nu - 2925008/3575*common_44 - 27938832/715*common_47 - 3670308/55*common_48 - 132972/55*common_49 + (2498544/143)*common_55*nu*potential - 631272/143*common_56 - 301572/715*common_57 - 1857204/715*common_61 - 2109072/143*common_62 - 548172/55*common_63 + (139248/65)*common_8*potential*speed_squared - 1279704/1925*common_97, (255456/715)*common_105 - 2370384/715*common_106 + (8641944/715)*common_26 - 1293696/65*common_28 - 2863656/715*common_34 - 780048/715*common_35 + (8470728/715)*common_36 + (407712/143)*common_37 + (2491488/715)*common_38 - 8763312/5005*common_68 - 488916/143*common_72 + (355776/715)*common_83),
        (-11072/33*common_0*common_107 - 1450764/11*common_0*common_121 + (75776/55)*common_1*common_107*speed_squared + (70668/11)*common_1*common_117*common_3 - 2976/55*common_1*common_145 + (4319028/11)*common_1*common_3*common_55*speed_squared + (156516/55)*common_1*common_4*common_44 + (1327764/11)*common_1*common_44*common_55 - 1245456/55*common_107*common_24 + (353552/55)*common_107*common_8*nu - 283500*common_115*common_27 + 94500*common_115*common_3*nu - 24084/11*common_117*common_6 - 3460152/55*common_123*common_17 + (2144/33)*common_145*nu + (494676/11)*common_3*common_4*common_8*nu - 1461852/11*common_4*common_50 - 10188/11*common_44*common_7 + (1110264/55)*common_44*common_8*nu*speed_squared - 420412/11*common_45*common_55, -1417216/165*common_132 - 1259904/55*common_133 - 2041200/11*common_134 + (831744/11)*common_135 + (1712896/55)*common_136 - 511344/11*common_137 + (2327136/11)*common_138 + (3980352/55)*common_139 + (6123600/11)*common_140 - 2657856/11*common_141 + (1499184/11)*common_142 - 6900096/11*common_143, (3916512/11)*common_1*common_3*common_8*speed_squared + (1977696/11)*common_1*common_44*common_8 + (162432/55)*common_107*nu - 605568/55*common_108 - 608544/11*common_109 - 1072992/55*common_110 - 1330848/11*common_111 - 4536000/11*common_112 - 1972512/55*common_113 + (682848/55)*common_3*common_4*nu + (1512000/11)*common_3*common_55*nu + (326688/55)*common_44*nu*speed_squared, (1512000/11)*common_1*common_3*common_33 + (308736/11)*common_3*nu*radial_speed*speed_squared + (189568/11)*common_44*nu*radial_speed - 504000/11*common_66 - 632064/11*common_67 - 894528/11*common_69, -18720/11*common_45 - 21600/11*common_46 + (64608/11)*common_47 + (64800/11)*common_48 + (60576/11)*common_49 - 194400/11*common_50),
    )

def i4_tail(potential, radial_speed, speed_squared, nu, log_one, log_infrared):
    common_0 = nu*speed_squared
    common_1 = nu**2
    common_2 = common_1*potential
    common_3 = 12*nu
    common_4 = potential*radial_speed
    common_5 = nu*potential
    common_6 = common_5*radial_speed
    common_7 = common_2*radial_speed
    common_8 = (156/11)*common_1 - 156/11*nu + 156/55
    common_9 = 36*nu
    common_10 = potential*speed_squared
    common_11 = potential**2
    common_12 = common_0*potential
    common_13 = common_11*nu
    common_14 = radial_speed**2
    common_15 = common_14*potential
    common_16 = common_14*common_5
    common_17 = common_2*speed_squared
    common_18 = common_1*common_11
    common_19 = common_14*common_2
    common_20 = (192/11)*radial_speed
    common_21 = 24 - 72*nu
    common_22 = radial_speed**3
    common_23 = common_11*radial_speed
    common_24 = common_0*common_4
    common_25 = common_22*common_5
    common_26 = common_2*common_22
    common_27 = common_18*radial_speed
    common_28 = (504/11)*radial_speed
    common_29 = potential**3
    common_30 = radial_speed**4
    common_31 = common_30*potential
    common_32 = speed_squared**2
    common_33 = common_11*speed_squared
    common_34 = common_32*common_5
    common_35 = common_10*common_14
    common_36 = common_1*common_29
    common_37 = common_13*common_14
    common_38 = common_1*speed_squared
    common_39 = common_11*common_38
    common_40 = common_14*common_17
    common_41 = common_4*speed_squared
    common_42 = common_22*potential
    common_43 = common_13*radial_speed
    common_44 = 1440*common_5 - 480*potential
    common_45 = radial_speed**5
    common_46 = common_45*potential
    common_47 = common_29*radial_speed
    common_48 = common_45*common_5
    common_49 = common_29*nu
    common_50 = common_49*radial_speed
    common_51 = common_32*common_4
    common_52 = common_10*common_22
    common_53 = common_23*speed_squared
    common_54 = common_11*common_22
    common_55 = common_32*common_6
    common_56 = common_0*common_42
    common_57 = common_0*common_23
    common_58 = common_13*common_22
    common_59 = common_2*common_45
    common_60 = common_36*radial_speed
    common_61 = common_32*common_7
    common_62 = common_17*common_22
    common_63 = common_23*common_38
    common_64 = common_18*common_22
    common_65 = common_0*common_15
    common_66 = common_32*potential
    common_67 = common_30*common_5
    common_68 = common_0*common_11
    common_69 = common_11*common_14
    common_70 = potential**4
    common_71 = common_70*nu
    common_72 = radial_speed**6
    common_73 = common_72*potential
    common_74 = speed_squared**3
    common_75 = common_74*potential
    common_76 = common_0*common_29
    common_77 = common_1*common_70
    common_78 = common_14*common_29
    common_79 = common_13*common_30
    common_80 = common_13*common_32
    common_81 = common_2*common_74
    common_82 = common_69*speed_squared
    common_83 = common_17*common_30
    common_84 = common_14*common_36
    common_85 = common_38*common_69
    common_86 = radial_speed**7
    return (
        (1 - 3*nu, 0, 0, 0, 0),
        (-147/22*common_0 + (279/22)*common_1*speed_squared - 105/11*common_2 + (61/11)*nu*potential - 10/11*potential + (103/110)*speed_squared, -72/11*common_1*radial_speed + (72/11)*nu*radial_speed - 72/55*radial_speed, (78/11)*common_1 - 78/11*nu + 78/55, 0, 0),
        (0, 4 - common_3, 0, 0, 0),
        ((19/55)*common_4 + (14/11)*common_6 - 102/11*common_7, -222/11*common_0 + (486/11)*common_1*speed_squared - 504/11*common_2 + (328/11)*nu*potential - 284/55*potential + (134/55)*speed_squared, -60/11*common_1*radial_speed + (60/11)*nu*radial_speed - 12/11*radial_speed, common_8, 0),
        (common_3*potential - 4*potential, 0, 12 - common_9, 0, 0),
        (-67/11*common_10 + (229/11)*common_11 + (236/11)*common_12 - 782/11*common_13 - 57/55*common_15 + (24/11)*common_16 - 192/11*common_17 + (342/11)*common_18 + (108/11)*common_19, -564/11*common_6 - 492/11*common_7 + (1092/55)*potential*radial_speed, -606/11*common_0 + (1398/11)*common_1*speed_squared - 1920/11*common_2 + (1392/11)*nu*potential - 252/11*potential + (342/55)*speed_squared, common_1*common_20 - common_20*nu + (192/55)*radial_speed, common_8),
        (-common_4*common_9 + 12*potential*radial_speed, 120*common_5 - 40*potential, 0, common_21, 0),
        (72*common_1*potential*radial_speed*speed_squared + (3172/11)*common_11*nu*radial_speed + (57/11)*common_22*potential - 4888/55*common_23 - 60*common_24 - 120/11*common_25 - 540/11*common_26 - 708/11*common_27 + (81/5)*potential*radial_speed*speed_squared, -2252/55*common_10 + (11288/55)*common_11 + (1592/11)*common_12 - 7988/11*common_13 - 3504/55*common_15 + (2184/11)*common_16 - 1680/11*common_17 + (4116/11)*common_18 + (720/11)*common_19, -4992/11*common_6 - 1344/11*common_7 + (8556/55)*potential*radial_speed, -1404/11*common_0 + (2988/11)*common_1*speed_squared - 4656/11*common_2 + (3600/11)*nu*potential - 3336/55*potential + (876/55)*speed_squared, common_1*common_28 - common_28*nu + (504/55)*radial_speed),
        (-common_10*common_9 + 28*common_11 - 84*common_13 + 180*common_14*nu*potential - 60*common_15 + 12*potential*speed_squared, -504*common_6 + 168*potential*radial_speed, 576*common_5 - 192*potential, 0, common_21),
        (522*common_1*common_11*common_14 + (3780/11)*common_1*common_30*potential + 72*common_1*common_32*potential + 597*common_11*common_14 + (2240/11)*common_11*nu*speed_squared + (2940/11)*common_14*nu*potential*speed_squared + (7896/11)*common_29*nu - 2512/11*common_29 + (840/11)*common_30*nu*potential - 399/11*common_31 + (81/5)*common_32*potential - 1987/55*common_33 - 60*common_34 - 720/11*common_35 - 1560/11*common_36 - 1986*common_37 - 2592/11*common_39 - 5580/11*common_40, (9648/11)*common_1*potential*radial_speed*speed_squared + (53232/11)*common_11*nu*radial_speed + (3732/11)*common_22*potential - 79104/55*common_23 - 3048/11*common_24 - 11400/11*common_25 - 5760/11*common_26 - 17328/11*common_27 + (3312/55)*potential*radial_speed*speed_squared, -4788/55*common_10 + (51156/55)*common_11 + (3996/11)*common_12 - 37692/11*common_13 - 7236/11*common_15 + (22716/11)*common_16 - 8220/11*common_17 + (22908/11)*common_18 + (2628/11)*common_19, -19848/11*common_6 - 1272/11*common_7 + (6504/11)*potential*radial_speed, -1908/11*common_0 + (3492/11)*common_1*speed_squared - 5160/11*common_2 + (4104/11)*nu*potential - 768/11*potential + (276/11)*speed_squared),
        (-240*common_23 + 540*common_24 - 1260*common_25 - 180*common_41 + 420*common_42 + 720*common_43, 328*common_11 - 648*common_12 - 984*common_13 + 3240*common_14*nu*potential - 1080*common_15 + 216*potential*speed_squared, -3240*common_6 + 1080*potential*radial_speed, common_44, 0),
        ((3591/11)*common_46 + (121168/55)*common_47 - 7560/11*common_48 - 77212/11*common_50 - 2331/11*common_51 + (3444/11)*common_52 + (79236/55)*common_53 - 58512/11*common_54 + (9180/11)*common_55 - 17220/11*common_56 - 61020/11*common_57 + (196596/11)*common_58 - 34020/11*common_59 + (17196/11)*common_60 - 15120/11*common_61 + (54180/11)*common_62 + (34092/11)*common_63 - 55620/11*common_64, (152532/11)*common_1*common_11*common_14 + 5040*common_1*common_30*potential + (12816/11)*common_1*common_32*potential + (137928/11)*common_11*common_14 + (57248/11)*common_11*nu*speed_squared + (5004/11)*common_14*potential*speed_squared + (77496/11)*common_29*nu - 120968/55*common_29 + 7560*common_30*nu*potential - 2520*common_31 + (6876/55)*common_32*potential - 68908/55*common_33 - 5688/11*common_34 - 20472/11*common_36 - 467220/11*common_37 - 42288/11*common_39 - 87840/11*common_40 - 7200/11*common_65, (58860/11)*common_1*potential*radial_speed*speed_squared + (375720/11)*common_11*nu*radial_speed + (47376/11)*common_22*potential - 108408/11*common_23 - 147780/11*common_25 - 30420/11*common_26 - 153960/11*common_27 - 9612/11*common_41 + (24300/11)*nu*potential*radial_speed*speed_squared, (1104/5)*common_10 + (106272/55)*common_11 - 384*common_12 - 82512/11*common_13 - 33984/11*common_15 + (106560/11)*common_16 - 2016*common_17 + (61392/11)*common_18 + (4320/11)*common_19, (10944/11)*common_4 - 34920/11*common_6 + (3240/11)*common_7),
        (-88*common_29 - 3780*common_31 - 276*common_33 + 540*common_34 + 2520*common_35 - 6300*common_37 + 264*common_49 - 7560*common_65 - 180*common_66 + 11340*common_67 + 828*common_68 + 2100*common_69, -3360*common_23 + 11880*common_24 - 27720*common_25 - 3960*common_41 + 9240*common_42 + 10080*common_43, 1344*common_11 - 5184*common_12 - 4032*common_13 + 25920*common_14*nu*potential - 8640*common_15 + 1728*potential*speed_squared, -10800*common_6 + 3600*potential*radial_speed, common_44),
        ((654570/11)*common_1*common_11*common_30 + (39960/11)*common_1*common_11*common_32 + (268380/11)*common_1*common_14*common_32*potential + (31728/11)*common_1*common_29*speed_squared + 34020*common_1*common_72*potential - 13041/11*common_10*common_30 + 107622*common_11*common_14*nu*speed_squared + (643437/11)*common_11*common_30 + (16407/11)*common_11*common_32 + (752592/11)*common_14*common_29*nu + (26649/11)*common_14*common_32*potential - 115920/11*common_16*common_32 + (22784/11)*common_29*speed_squared + (117180/11)*common_30*nu*potential*speed_squared + (3832/11)*common_70 - 9964/11*common_71 + 7560*common_72*nu*potential - 3591*common_73 + (9180/11)*common_74*nu*potential - 2331/11*common_75 - 77400/11*common_76 - 2532/11*common_77 - 1153344/55*common_78 - 2173890/11*common_79 - 64512/11*common_80 - 15120/11*common_81 - 150048/5*common_82 - 657720/11*common_83 - 216936/11*common_84 - 46278*common_85, (263844/11)*common_46 + (1341264/55)*common_47 - 778680/11*common_48 - 880848/11*common_50 - 6192/11*common_51 - 12012*common_52 + (2247048/55)*common_53 - 133992*common_54 + (50760/11)*common_55 + 28560*common_56 - 1633248/11*common_57 + 455232*common_58 - 635040/11*common_59 + (292656/11)*common_60 - 300240/11*common_61 + 95760*common_62 + (787392/11)*common_63 - 151200*common_64, (1573596/11)*common_1*common_11*common_14 + (379260/11)*common_1*common_30*potential + (97308/11)*common_1*common_32*potential + (1185840/11)*common_11*common_14 + (535836/11)*common_11*nu*speed_squared + (205200/11)*common_14*potential*speed_squared + (143904/11)*common_29*nu - 43776/11*common_29 + (1283940/11)*common_30*nu*potential - 414792/11*common_31 + (7236/11)*common_32*nu*potential - 717336/55*common_33 - 67872/11*common_36 - 4085820/11*common_37 - 308796/11*common_39 - 649080/11*common_40 - 586440/11*common_65 - 27432/55*common_66, (192888/11)*common_1*potential*radial_speed*speed_squared + (1238256/11)*common_11*nu*radial_speed + (264672/11)*common_22*potential - 1729296/55*common_23 - 828360/11*common_25 - 82440/11*common_26 - 625776/11*common_27 - 472392/55*common_41 + (274392/11)*nu*potential*radial_speed*speed_squared, (58152/11)*common_1*common_11 + (51552/55)*common_11 - 39144/11*common_12 - 47592/11*common_13 + (211320/11)*common_14*nu*potential - 66816/11*common_15 - 18936/11*common_17 - 5400/11*common_19 + (66864/55)*potential*speed_squared),
        (60480*common_11*common_22*nu + 6048*common_11*radial_speed*speed_squared + 113400*common_22*nu*potential*speed_squared + 504*common_29*radial_speed + 6300*common_32*potential*radial_speed + 41580*common_45*potential - 124740*common_48 - 1512*common_50 - 37800*common_52 - 20160*common_54 - 18900*common_55 - 18144*common_57, 32760*common_11*common_14 + 11880*common_11*nu*speed_squared + 65520*common_14*potential*speed_squared + 320*common_29 + 294840*common_30*nu*potential - 98280*common_31 + 14040*common_32*nu*potential - 3960*common_33 - 98280*common_37 - 960*common_49 - 196560*common_65 - 4680*common_66, -15120*common_23 + 113400*common_24 - 264600*common_25 - 37800*common_41 + 88200*common_42 + 45360*common_43, 1008*common_11 - 21168*common_12 - 3024*common_13 + 105840*common_14*nu*potential - 35280*common_15 + 7056*potential*speed_squared, -15120*common_6 + 5040*potential*radial_speed),
        (-71820*common_0*common_46 - 22399560/11*common_0*common_54 + (8777160/11)*common_1*common_11*common_22*speed_squared + (2767716/11)*common_1*common_22*common_29 + 861840*common_1*common_45*potential*speed_squared + (84072/11)*common_1*common_70*radial_speed + (642600/11)*common_1*common_74*potential*radial_speed - 8505*common_10*common_45 + (6405264/11)*common_11*common_22*speed_squared + (3120660/11)*common_11*common_32*nu*radial_speed + 2586060*common_11*common_45*nu - 762480*common_11*common_45 - 827820*common_18*common_45 - 442260*common_2*common_86 + (2245572/11)*common_22*common_29 + (1512000/11)*common_22*common_32*nu*potential - 7507836/11*common_22*common_49 - 4301136/55*common_23*common_32 - 5046300/11*common_26*common_32 - 1495476/11*common_27*common_32 + (1770804/11)*common_29*nu*radial_speed*speed_squared - 292005/11*common_32*common_42 - 680220/11*common_38*common_47 - 532836/11*common_47*speed_squared - 98280*common_5*common_86 - 296100/11*common_6*common_74 + (609344/55)*common_70*radial_speed - 401576/11*common_71*radial_speed + (69615/11)*common_74*potential*radial_speed + 46683*common_86*potential, -6252120/11*common_0*common_31 + (21874860/11)*common_1*common_11*common_30 + (1002960/11)*common_1*common_11*common_32 + (6335280/11)*common_1*common_14*common_32*potential + (216504/11)*common_1*common_29*speed_squared + 771120*common_1*common_72*potential + 3207492*common_11*common_14*nu*speed_squared + (18713736/11)*common_11*common_30 + (2688732/55)*common_11*common_32 + (8100072/11)*common_14*common_29*nu + (123480/11)*common_14*common_32*nu*potential - 246456/11*common_15*common_32 + (910648/55)*common_29*speed_squared + (2456244/11)*common_30*potential*speed_squared + (627104/11)*common_70*nu - 974704/55*common_70 + 808920*common_72*nu*potential - 278208*common_73 + (87480/11)*common_74*nu*potential - 15516/11*common_75 - 573432/11*common_76 - 122688/11*common_77 - 11852568/55*common_78 - 63724140/11*common_79 - 1956528/11*common_80 - 360720/11*common_81 - 4568796/5*common_82 - 15286320/11*common_83 - 3810600/11*common_84 - 1341540*common_85, (4524660/11)*common_46 + (246192/55)*common_47 - 13891500/11*common_48 - 345984/11*common_50 + (419256/11)*common_51 - 3491964/11*common_52 + (24065424/55)*common_53 - 1361304*common_54 - 1056780/11*common_55 + (10183320/11)*common_56 - 17030088/11*common_57 + 4682160*common_58 - 5318460/11*common_59 + (1077792/11)*common_60 - 2685420/11*common_61 + (9220680/11)*common_62 + (7849224/11)*common_63 - 1789200*common_64, (7652232/11)*common_14*common_18 + (1335600/11)*common_2*common_30 + (387504/11)*common_2*common_32 + (1454928/55)*common_29 - 2682288/11*common_31 - 3051432/55*common_33 + (288864/11)*common_34 + (1676808/11)*common_35 + (146544/11)*common_36 - 16700040/11*common_37 - 1194192/11*common_39 - 2509920/11*common_40 - 944880/11*common_49 - 5029920/11*common_65 - 527256/55*common_66 + (8366400/11)*common_67 + (2192112/11)*common_68 + (4742640/11)*common_69, (238896/11)*common_1*potential*radial_speed*speed_squared + (1379952/11)*common_11*nu*radial_speed + 54432*common_22*potential - 1823472/55*common_23 - 171360*common_25 - 5040*common_26 - 936432/11*common_27 - 1341144/55*common_41 + (814464/11)*nu*potential*radial_speed*speed_squared),
    )

def i5(potential, radial_speed, speed_squared, nu, log_one, log_infrared):
    common_0 = nu*speed_squared
    common_1 = nu**2
    common_2 = common_1*potential
    common_3 = common_1*radial_speed
    common_4 = potential**2
    common_5 = speed_squared**2
    common_6 = common_0*potential
    common_7 = common_4*nu
    common_8 = common_5*nu
    common_9 = radial_speed**2
    common_10 = common_9*potential
    common_11 = nu**3
    common_12 = common_11*common_4
    common_13 = common_11*common_5
    common_14 = common_2*common_9
    common_15 = potential*radial_speed
    common_16 = radial_speed*speed_squared
    common_17 = common_2*radial_speed
    common_18 = common_1*speed_squared
    common_19 = common_18*radial_speed
    common_20 = common_11*common_15
    common_21 = common_11*speed_squared
    common_22 = nu*potential
    common_23 = common_11*potential
    common_24 = common_1*common_9
    common_25 = common_9*nu
    common_26 = (80/39)*common_11
    common_27 = -200/39*common_24 + (40/13)*common_25 + common_26*common_9 - 20/39*common_9
    common_28 = (40/13)*radial_speed
    common_29 = (850/273)*common_1 - 340/273*common_11 - 170/91*nu + 85/273
    common_30 = 10*nu
    common_31 = common_22*radial_speed
    common_32 = potential*speed_squared
    common_33 = common_32*radial_speed
    common_34 = radial_speed**3
    common_35 = common_34*potential
    common_36 = common_4*radial_speed
    common_37 = common_0*common_15
    common_38 = common_22*common_34
    common_39 = nu*radial_speed
    common_40 = common_39*common_4
    common_41 = common_16*common_2
    common_42 = common_11*radial_speed
    common_43 = common_32*common_42
    common_44 = common_2*common_34
    common_45 = common_3*common_4
    common_46 = common_23*common_34
    common_47 = common_12*radial_speed
    common_48 = common_2*speed_squared
    common_49 = common_22*common_9
    common_50 = common_1*common_4
    common_51 = (420/13)*common_1 - 560/13*nu + 140/13
    common_52 = potential**3
    common_53 = radial_speed**4
    common_54 = common_53*potential
    common_55 = common_5*potential
    common_56 = common_0*common_4
    common_57 = common_10*speed_squared
    common_58 = common_1*common_52
    common_59 = common_4*common_9
    common_60 = common_7*common_9
    common_61 = common_2*common_53
    common_62 = common_2*common_5
    common_63 = common_23*common_53
    common_64 = common_12*speed_squared
    common_65 = common_10*common_21
    common_66 = common_50*common_9
    common_67 = common_0*radial_speed
    common_68 = common_11*common_16
    common_69 = common_11*common_9
    common_70 = radial_speed**5
    common_71 = common_70*potential
    common_72 = common_52*radial_speed
    common_73 = common_22*common_70
    common_74 = common_39*common_52
    common_75 = common_15*common_5
    common_76 = common_32*common_34
    common_77 = common_16*common_4
    common_78 = common_34*common_4
    common_79 = common_31*common_5
    common_80 = common_0*common_35
    common_81 = common_0*common_36
    common_82 = common_34*common_7
    common_83 = common_2*common_70
    common_84 = common_3*common_52
    common_85 = common_23*common_70
    common_86 = common_42*common_52
    common_87 = common_17*common_5
    common_88 = common_34*common_48
    common_89 = common_18*common_36
    common_90 = common_13*common_15
    common_91 = common_11*common_32
    common_92 = common_34*common_91
    common_93 = common_12*common_16
    common_94 = common_34*common_50
    common_95 = common_12*common_34
    common_96 = common_14*speed_squared
    common_97 = common_1*common_5
    common_98 = common_10*common_11
    common_99 = 120 - 240*nu
    common_100 = common_4*speed_squared
    common_101 = common_22*common_5
    common_102 = common_18*common_4
    common_103 = potential**4
    common_104 = common_103*nu
    common_105 = radial_speed**6
    common_106 = common_105*potential
    common_107 = speed_squared**3
    common_108 = common_52*speed_squared
    common_109 = common_107*common_22
    common_110 = common_0*common_52
    common_111 = common_103*common_11
    common_112 = common_4*common_53
    common_113 = common_4*common_5
    common_114 = common_52*common_9
    common_115 = common_0*common_54
    common_116 = common_53*common_7
    common_117 = common_4*common_8
    common_118 = common_105*common_2
    common_119 = (3600/13)*common_118
    common_120 = common_105*common_23
    common_121 = common_107*common_23
    common_122 = common_21*common_52
    common_123 = common_0*common_59
    common_124 = common_53*common_91
    common_125 = common_50*common_53
    common_126 = common_24*common_52
    common_127 = common_14*common_5
    common_128 = common_64*common_9
    common_129 = common_52*nu
    common_130 = common_22*common_53
    common_131 = common_11*common_52
    common_132 = common_0*common_10
    common_133 = common_13*potential
    common_134 = common_12*common_9
    common_135 = radial_speed**7
    common_136 = common_135*potential
    common_137 = common_103*radial_speed
    common_138 = common_135*common_22
    common_139 = common_103*common_39
    common_140 = common_107*common_15
    common_141 = common_32*common_70
    common_142 = common_16*common_52
    common_143 = common_4*common_70
    common_144 = common_34*common_52
    common_145 = common_107*common_31
    common_146 = common_0*common_71
    common_147 = common_0*common_72
    common_148 = common_7*common_70
    common_149 = common_129*common_34
    common_150 = common_135*common_2
    common_151 = common_103*common_3
    common_152 = common_135*common_23
    common_153 = common_103*common_42
    common_154 = common_35*common_5
    common_155 = common_36*common_5
    common_156 = common_100*common_34
    common_157 = common_38*common_5
    common_158 = common_113*common_39
    common_159 = common_0*common_78
    common_160 = common_107*common_17
    common_161 = common_48*common_70
    common_162 = common_18*common_72
    common_163 = common_107*common_20
    common_164 = common_70*common_91
    common_165 = common_131*common_16
    common_166 = common_50*common_70
    common_167 = common_34*common_58
    common_168 = common_12*common_70
    common_169 = common_131*common_34
    common_170 = common_44*common_5
    common_171 = common_113*common_3
    common_172 = common_18*common_78
    common_173 = common_13*common_35
    common_174 = common_13*common_36
    common_175 = common_34*common_64
    common_176 = common_10*common_5
    common_177 = common_48*common_53
    common_178 = common_5*common_50
    common_179 = common_59*speed_squared
    common_180 = potential**5
    common_181 = radial_speed**8
    common_182 = speed_squared**4
    return (
        (1 - 2*nu, 0, 0, 0, 0, 0),
        (-511/78*common_0 + (137/13)*common_1*speed_squared - 109/13*common_2 + (139/26)*nu*potential - 25/26*potential + (79/78)*speed_squared, -60/13*common_3 + (80/13)*nu*radial_speed - 20/13*radial_speed, (70/13)*common_1 - 280/39*nu + 70/39, 0, 0, 0),
        ((2633/52)*common_1*common_4 + (29195/728)*common_1*common_5 + (14165/2184)*common_1*potential*speed_squared - 1151/6552*common_10 + (38849/3276)*common_11*common_9*potential + (213667/3276)*common_11*potential*speed_squared - 1721/78*common_12 - 19265/364*common_13 - 15053/2184*common_14 - 269/156*common_4 + (1475/2184)*common_5 - 1101/56*common_6 - 1045/52*common_7 - 20687/2184*common_8 + (3499/2184)*common_9*nu*potential + (21053/6552)*potential*speed_squared, (8110/273)*common_11*radial_speed*speed_squared - 15569/3276*common_15 - 190/273*common_16 - 22625/1092*common_17 - 9250/273*common_19 - 36895/1638*common_20 + (27463/1092)*nu*potential*radial_speed + (870/91)*nu*radial_speed*speed_squared, -3365/273*common_0 + (3775/91)*common_18 - 715/21*common_2 - 3225/91*common_21 - 1179/91*common_22 + (2195/63)*common_23 + common_27 + (7841/819)*potential + (275/273)*speed_squared, (200/39)*common_1*radial_speed - common_26*radial_speed - common_28*nu + (20/39)*radial_speed, common_29, 0),
        (0, 5 - common_30, 0, 0, 0, 0),
        ((37/78)*common_15 - 105/13*common_17 + (125/78)*common_31, -2075/78*common_0 + (625/13)*common_1*speed_squared - 625/13*common_2 + (2725/78)*nu*potential - 535/78*potential + (275/78)*speed_squared, -100/13*common_3 + (400/39)*nu*radial_speed - 100/39*radial_speed, (210/13)*common_1 - 280/13*nu + 70/13, 0, 0),
        ((1489/728)*common_33 + (1151/2184)*common_35 + (1907/468)*common_36 + (159/728)*common_37 - 2967/728*common_38 + (3797/156)*common_40 - 23225/728*common_41 + (17819/364)*common_43 + (7465/728)*common_44 - 1445/52*common_45 - 11885/1092*common_46 - 12773/234*common_47, (9105/26)*common_1*common_4 + (363925/2184)*common_1*common_5 - 1877/312*common_10 + (7175/156)*common_11*common_9*potential + (356315/1092)*common_11*potential*speed_squared - 16985/117*common_12 - 256535/1092*common_13 - 995/104*common_14 - 1733/117*common_4 - 41305/2184*common_48 + (5855/2184)*common_5 - 122329/2184*common_6 - 5009/39*common_7 - 82555/2184*common_8 + (1669/104)*common_9*nu*potential + (5783/728)*potential*speed_squared, (2030/39)*common_11*radial_speed*speed_squared - 70/39*common_16 - 36370/273*common_17 - 2450/39*common_19 - 1055/273*common_20 + (13856/273)*nu*potential*radial_speed + (770/39)*nu*radial_speed*speed_squared + (218/91)*potential*radial_speed, -3645/91*common_0 + (35375/273)*common_18 - 10895/91*common_2 - 29585/273*common_21 - 2577/91*common_22 + (30455/273)*common_23 + common_27 + (7361/273)*potential + (965/273)*speed_squared, (6200/273)*common_1*radial_speed - 1240/91*common_39 - 2480/273*common_42 + (620/273)*radial_speed, common_29),
        (common_30*potential - 5*potential, 0, 20 - 40*nu, 0, 0, 0),
        (-37/26*common_10 + (120/13)*common_14 - 314/39*common_32 + (343/13)*common_4 - 340/13*common_48 + (35/13)*common_49 + (470/13)*common_50 + (905/39)*common_6 - 865/13*common_7, -690/13*common_17 - 575/13*common_31 + (355/13)*potential*radial_speed, -1250/13*common_0 + (2400/13)*common_1*speed_squared - 3030/13*common_2 + (7570/39)*nu*potential - 1600/39*potential + (150/13)*speed_squared, common_28 + (120/13)*common_3 - 160/13*common_39, common_51, 0),
        ((12255/364)*common_1*common_4*speed_squared + (2105/182)*common_1*common_9*potential*speed_squared + (1198/63)*common_11*common_4*common_9 + (27203/273)*common_11*common_5*potential + (6821/78)*common_11*common_52 + (47741/3276)*common_4*speed_squared + (3553/91)*common_5*nu*potential + (5687/26)*common_52*nu - 9329/156*common_52 + (4005/364)*common_53*nu*potential - 5755/2184*common_54 - 1136/273*common_55 - 4489/1092*common_56 - 829/182*common_57 - 13335/52*common_58 - 3323/819*common_59 - 200/273*common_60 - 400/91*common_61 - 63775/546*common_62 - 500/273*common_63 - 237869/1638*common_64 - 4313/91*common_65 - 18275/273*common_66 + (2713/182)*common_9*nu*potential*speed_squared, (5545/91)*common_1*common_34*potential + (26870/91)*common_11*potential*radial_speed*speed_squared + (11293/546)*common_34*potential - 13659/182*common_37 - 16041/182*common_38 + (328721/1092)*common_4*nu*radial_speed + (38027/3276)*common_4*radial_speed - 3215/273*common_41 - 366265/1092*common_45 - 11260/273*common_46 - 690815/1638*common_47 + (6416/273)*potential*radial_speed*speed_squared, (177225/91)*common_1*common_4 + (109875/182)*common_1*common_5 + (25/42)*common_1*common_9*potential - 28447/546*common_10 + (16825/91)*common_11*common_9*potential + (396415/273)*common_11*potential*speed_squared - 679730/819*common_12 - 80775/91*common_13 - 5555/63*common_4 - 225215/546*common_48 + (125/14)*common_5 - 10789/182*common_6 - 182855/273*common_7 - 23925/182*common_8 + (22489/182)*common_9*nu*potential + (6427/546)*potential*speed_squared, (7611/91)*common_15 + (60/91)*common_16 - 146345/273*common_17 + (5500/91)*common_19 + (68800/273)*common_20 - 3491/91*common_31 - 1340/91*common_67 - 5140/91*common_68, -8530/91*common_0 + (25650/91)*common_18 - 75820/273*common_2 - 20550/91*common_21 - 3064/91*common_22 + (65090/273)*common_23 + (3200/91)*common_24 - 1920/91*common_25 - 1280/91*common_69 + (320/91)*common_9 + (4559/91)*potential + (850/91)*speed_squared, (3200/91)*common_1*radial_speed - 1920/91*common_39 - 1280/91*common_42 + (320/91)*radial_speed),
        (-30*common_31 + 15*potential*radial_speed, 130*common_22 - 65*potential, 0, 60 - 120*nu, 0, 0),
        ((1260/13)*common_1*potential*radial_speed*speed_squared + (185/26)*common_34*potential - 4442/39*common_36 - 835/13*common_37 - 175/13*common_38 + (10085/39)*common_4*nu*radial_speed - 600/13*common_44 - 750/13*common_45 + (277/13)*potential*radial_speed*speed_squared, -2315/26*common_10 + (1110/13)*common_14 - 2965/39*common_32 + (1040/3)*common_4 - 4070/13*common_48 + (2680/13)*common_49 + 540*common_50 + (8740/39)*common_6 - 2750/3*common_7, -2810/13*common_17 - 21130/39*common_31 + (10840/39)*potential*radial_speed, -3910/13*common_0 + (7320/13)*common_1*speed_squared - 10890/13*common_2 + (9970/13)*nu*potential - 2200/13*potential + (490/13)*speed_squared, (1080/13)*common_3 - 1440/13*common_39 + (360/13)*radial_speed, common_51),
        ((5755/312)*common_71 + (35077/117)*common_72 - 4005/52*common_73 - 69577/52*common_74 + (307/91)*common_75 + (3340/273)*common_76 - 1780/13*common_77 + (3256/273)*common_78 - 7946/91*common_79 - 5555/182*common_80 + (6961/78)*common_81 + (10643/182)*common_82 + (400/13)*common_83 + (105425/78)*common_84 + (500/39)*common_85 + (7378/117)*common_86 + (67985/182)*common_87 - 13725/182*common_88 - 1970/39*common_89 - 35829/91*common_90 + (62695/273)*common_92 + (5008/13)*common_93 + (74545/182)*common_94 - 37685/273*common_95, (4545/7)*common_1*common_4*speed_squared + (4400/91)*common_1*common_53*potential + (59765/273)*common_11*common_4*common_9 + (311515/273)*common_11*common_5*potential + (173210/117)*common_11*common_52 + (123883/546)*common_4*common_9 + (50363/273)*common_4*speed_squared + (70281/182)*common_5*nu*potential + (168929/52)*common_52*nu - 83449/117*common_52 + (153135/364)*common_53*nu*potential - 254635/2184*common_54 - 803/21*common_55 - 7087/42*common_56 - 406/13*common_57 - 26945/6*common_58 - 67848/91*common_60 - 217585/182*common_62 - 69050/273*common_63 - 172090/91*common_64 - 2765/13*common_65 - 146555/182*common_66 + (6497/26)*common_9*nu*potential*speed_squared - 1225/2*common_96, (50587/546)*common_33 + (130513/546)*common_35 - 61820/273*common_36 - 113249/182*common_37 - 11967/14*common_38 + (248870/91)*common_40 + (230095/182)*common_41 + (270575/273)*common_43 + (25135/182)*common_44 - 774370/273*common_45 + (31595/273)*common_46 - 182990/91*common_47, -37509/182*common_10 - 933290/273*common_12 - 247465/91*common_13 + (134065/546)*common_14 + (6329/182)*common_32 - 99356/273*common_4 - 1390135/546*common_48 + (67449/182)*common_49 + (4995/182)*common_5 + (698350/91)*common_50 + (39971/182)*common_6 - 225908/91*common_7 - 74455/182*common_8 + (108205/21)*common_91 + (340625/182)*common_97 + (165365/273)*common_98, (4503/13)*common_15 + (2460/91)*common_16 - 50630/39*common_17 + (68700/91)*common_19 + (40330/39)*common_20 - 7158/13*common_31 - 23580/91*common_67 - 53940/91*common_68, -10450/91*common_0 + (28850/91)*common_18 - 85420/273*common_2 - 21830/91*common_21 - 88/7*common_22 + (68930/273)*common_23 + (9600/91)*common_24 - 5760/91*common_25 - 3840/91*common_69 + (960/91)*common_9 + (4239/91)*potential + (90/7)*speed_squared),
        (-75*common_10 + 50*common_4 - 30*common_6 - 100*common_7 + 150*common_9*nu*potential + 15*potential*speed_squared, -540*common_31 + 270*potential*radial_speed, 880*common_22 - 440*potential, 0, common_99, 0),
        ((5400/13)*common_1*common_4*common_9 + (1260/13)*common_1*common_5*potential + (4200/13)*common_1*common_53*potential - 358/39*common_100 - 835/13*common_101 - 1840/13*common_102 + (9990/13)*common_4*common_9 + (5800/39)*common_4*nu*speed_squared + (277/13)*common_5*potential + (12455/13)*common_52*nu - 5626/13*common_52 + (1225/13)*common_53*nu*potential - 1295/26*common_54 - 2215/26*common_57 - 3670/13*common_58 - 23580/13*common_60 + (3650/13)*common_9*nu*potential*speed_squared - 8100/13*common_96, (20730/13)*common_1*potential*radial_speed*speed_squared + (6250/13)*common_34*potential - 31745/13*common_36 - 7555/13*common_37 - 14275/13*common_38 + (75825/13)*common_4*nu*radial_speed - 8550/13*common_44 - 23670/13*common_45 + (2035/13)*potential*radial_speed*speed_squared, -1190*common_10 + 450*common_14 - 4150/13*common_32 + (91120/39)*common_4 - 27010/13*common_48 + 2720*common_49 + (54200/13)*common_50 + (14000/13)*common_6 - 253720/39*common_7, -7140/13*common_17 - 40920/13*common_31 + (19980/13)*potential*radial_speed, -9260/13*common_0 + (15720/13)*common_1*speed_squared - 1920*common_2 + 1860*nu*potential - 420*potential + (1340/13)*speed_squared, (1920/13)*common_3 - 2560/13*common_39 + (640/13)*radial_speed),
        ((20505/13)*common_1*common_103 + (67985/182)*common_1*common_107*potential + (2185/273)*common_1*common_4*common_5 + (269995/91)*common_1*common_4*common_9*speed_squared + (102895/42)*common_1*common_52*speed_squared + (17725/26)*common_1*common_53*potential*speed_squared + (69986/39)*common_103 - 101137/26*common_104 + (36045/52)*common_105*nu*potential - 17265/104*common_106 + (307/91)*common_107*potential - 163304/819*common_108 - 7946/91*common_109 + (187370/273)*common_11*common_4*common_5 + (51335/26)*common_11*common_4*common_53 + (241840/91)*common_11*common_5*common_9*potential + (347015/819)*common_11*common_52*common_9 - 1194325/1092*common_110 - 22576/39*common_111 - 527/52*common_112 - 4308/91*common_113 - 4542743/1638*common_114 - 8915/52*common_115 - 34215/52*common_116 - 15101/273*common_117 - common_119 - 1500/13*common_120 - 35829/91*common_121 - 417938/819*common_122 - 111415/91*common_123 - 20065/13*common_124 - 214445/52*common_125 - 2163820/273*common_126 - 190550/91*common_127 - 328960/91*common_128 + (101840/91)*common_4*common_9*speed_squared + (62795/182)*common_5*common_9*nu*potential + (1805/91)*common_5*common_9*potential + (11106383/1092)*common_52*common_9*nu + (685/104)*common_53*potential*speed_squared, (47235/52)*common_71 + (5827273/819)*common_72 - 3330*common_73 - 8272601/273*common_74 + (6290/91)*common_75 - 6475/26*common_76 - 14537/7*common_77 - 60305/26*common_78 - 199345/182*common_79 + (3650/13)*common_80 + (94071/91)*common_81 + (206235/26)*common_82 - 2400/13*common_83 + (8951275/273)*common_84 + (23850/13)*common_85 - 2535710/819*common_86 + (29605/7)*common_87 + (37425/13)*common_88 - 1004980/273*common_89 - 529370/91*common_90 + 1200*common_92 + (2350190/273)*common_93 + (174145/26)*common_94 - 39415/13*common_95, (90686/91)*common_100 + (197720/91)*common_101 + (1594070/273)*common_102 + (6445156/273)*common_129 + (73050/13)*common_130 + (10425920/819)*common_131 + (16660/13)*common_132 + (655740/91)*common_133 + (48970/39)*common_134 - 3336328/819*common_52 - 21600/13*common_54 - 1790/7*common_55 - 343858/273*common_56 + (4085/13)*common_57 - 10297100/273*common_58 + (67920/13)*common_59 - 207460/13*common_60 + (15475/13)*common_61 - 574650/91*common_62 - 47000/13*common_63 - 288500/21*common_64 + (15820/13)*common_65 - 35060/39*common_66 - 128465/13*common_96, (129070/13)*common_1*potential*radial_speed*speed_squared + (23370/13)*common_11*common_34*potential + (20210/13)*common_11*potential*radial_speed*speed_squared + (17360/13)*common_34*potential - 102980/39*common_36 - 41970/13*common_37 - 54690/13*common_38 + (229240/13)*common_4*nu*radial_speed - 10520/13*common_44 - 54940/3*common_45 - 71050/13*common_47 + (2350/13)*potential*radial_speed*speed_squared, (722950/39)*common_1*common_4 + (58475/13)*common_1*common_5 + (26415/13)*common_1*common_9*potential - 4050/13*common_10 + 1110*common_11*common_9*potential + (161170/13)*common_11*potential*speed_squared - 112380/13*common_12 - 78410/13*common_13 - 38260/39*common_4 - 103955/13*common_48 - 3165/13*common_49 + (1065/13)*common_5 - 71410/13*common_7 - 14005/13*common_8 + (10445/13)*nu*potential*speed_squared + (3160/13)*potential*speed_squared, (7300/13)*common_15 + (960/13)*common_16 - 17930/13*common_17 + 1600*common_19 + (19860/13)*common_20 - 15290/13*common_31 - 8000/13*common_67 - 15040/13*common_68),
        (-225*common_33 + 525*common_35 - 450*common_36 + 450*common_37 - 1050*common_38 + 900*common_40, -1725*common_10 + 860*common_4 - 690*common_6 - 1720*common_7 + 3450*common_9*nu*potential + 345*potential*speed_squared, -4800*common_31 + 2400*potential*radial_speed, 3600*common_22 - 1800*potential, 0, common_99),
        ((11655/26)*common_71 + (177557/39)*common_72 - 11025/13*common_73 - 394250/39*common_74 - 3600/13*common_75 + (10325/26)*common_76 + (19378/13)*common_77 - 89165/13*common_78 + (11475/13)*common_79 - 20650/13*common_80 - 61015/13*common_81 + (213140/13)*common_82 - 37800/13*common_83 + (34440/13)*common_84 - 22500/13*common_87 + (73500/13)*common_88 + (25530/13)*common_89 - 48900/13*common_94, -54440/39*common_100 - 11730/13*common_101 - 47160/13*common_102 + (609530/39)*common_129 + (106050/13)*common_130 + (13200/13)*common_132 - 266555/39*common_52 - 93975/26*common_54 + (3420/13)*common_55 + (218930/39)*common_56 + (6075/26)*common_57 - 71360/13*common_58 + (276620/13)*common_59 - 667730/13*common_60 + (80850/13)*common_61 + (27030/13)*common_62 + (193380/13)*common_66 - 169800/13*common_96, (175650/13)*common_1*potential*radial_speed*speed_squared - 10350/13*common_33 + (102350/13)*common_34*potential - 944800/39*common_36 - 1500/13*common_37 - 233900/13*common_38 + (2340100/39)*common_4*nu*radial_speed - 63450/13*common_44 - 293500/13*common_45, -106350/13*common_10 + (20250/13)*common_14 - 4070/13*common_32 + (117940/13)*common_4 - 113610/13*common_48 + (238200/13)*common_49 + (244620/13)*common_50 + (31880/13)*common_6 - 346960/13*common_7, -5400/13*common_17 - 137100/13*common_31 + (64500/13)*potential*radial_speed, -11820/13*common_0 + (17640/13)*common_1*speed_squared - 26880/13*common_2 + (26740/13)*nu*potential - 6100/13*potential + (1980/13)*speed_squared),
        ((189915/104)*common_136 - 2406592/117*common_137 - 396495/52*common_138 + (3967721/78)*common_139 + (2075/91)*common_140 - 109755/104*common_141 - 277335/91*common_142 - 32295/52*common_143 + (419564/13)*common_144 + (102525/91)*common_145 + (296505/52)*common_146 + (2574113/84)*common_147 + (445135/52)*common_148 - 5623719/52*common_149 + (39600/13)*common_150 - 1092665/39*common_151 + (16500/13)*common_152 + (626432/117)*common_153 - 225/2*common_154 + (242314/91)*common_155 - 134217/13*common_156 - 80625/26*common_157 - 463809/182*common_158 + (139227/13)*common_159 - 1102125/182*common_160 - 202725/26*common_161 - 22929775/546*common_162 + (662825/91)*common_163 + (171585/13)*common_164 + (636364/91)*common_165 + (2413025/52)*common_166 + (962200/13)*common_167 - 621115/26*common_168 - 104565/13*common_169 + (226000/13)*common_170 + (625295/91)*common_171 - 1277155/26*common_172 - 322100/13*common_173 - 1163646/91*common_174 + (541635/13)*common_175, (1392815/39)*common_1*common_103 + (1109655/182)*common_1*common_107*potential + (7353870/91)*common_1*common_4*common_9*speed_squared + (29544835/546)*common_1*common_52*speed_squared + (2856949/117)*common_103 - 2274838/39*common_104 + (1738665/52)*common_105*nu*potential - 936555/104*common_106 + (7825/91)*common_107*potential - 278805/182*common_109 + (1456740/91)*common_11*common_4*common_5 + (628205/13)*common_11*common_4*common_53 + (4183650/91)*common_11*common_5*common_9*potential + (3750125/91)*common_11*common_52*common_9 - 11807399/364*common_110 - 1430090/117*common_111 - 126416/91*common_113 - 7948292/91*common_114 - 1012575/52*common_115 - 1111960/13*common_116 - 431309/546*common_117 + common_119 - 222150/13*common_120 - 708515/91*common_121 - 3793280/273*common_122 - 90275/13*common_124 - 2188165/26*common_125 - 24484165/91*common_126 - 2091150/91*common_127 - 8555310/91*common_128 - 180825/182*common_176 - 459325/26*common_177 - 954995/273*common_178 + (662385/26)*common_4*common_53 + (23588/7)*common_4*common_9*nu*speed_squared + (1280281/91)*common_4*common_9*speed_squared + (732000/91)*common_5*common_9*nu*potential + (112754333/364)*common_52*common_9*nu + (557783/273)*common_52*speed_squared + (657075/104)*common_53*potential*speed_squared, (351020/13)*common_1*common_34*common_4 + (853925/13)*common_1*common_34*potential*speed_squared + (1464900/91)*common_1*common_5*potential*radial_speed + (33516320/91)*common_1*common_52*radial_speed + (8265350/91)*common_11*common_4*radial_speed*speed_squared + (424400/13)*common_11*common_70*potential + (2267830/13)*common_34*common_4*nu + (223500/13)*common_34*nu*potential*speed_squared + (152160/91)*common_5*potential*radial_speed + (19418068/273)*common_52*radial_speed + (198435/13)*common_70*potential - 684510/13*common_73 - 84302948/273*common_74 - 119775/13*common_76 - 790948/91*common_77 - 749500/13*common_78 - 758610/91*common_79 - 3837146/273*common_81 - 117925/13*common_83 - 18617710/273*common_86 - 12178300/273*common_89 - 3863220/91*common_90 - 204700/13*common_92 - 380400/13*common_95, (269824/273)*common_100 + (707190/91)*common_101 + (6296240/273)*common_102 + (631732/7)*common_129 + (480900/13)*common_130 + (16864730/273)*common_131 + (47940/13)*common_132 + (2634250/91)*common_133 + (36610/13)*common_134 - 2879788/273*common_52 - 151600/13*common_54 - 120700/91*common_55 + (232452/91)*common_56 + (57285/13)*common_57 - 15692280/91*common_58 + (1897940/39)*common_59 - 1964980/13*common_60 + (157525/13)*common_61 - 1683280/91*common_62 - 328050/13*common_63 - 5600460/91*common_64 + (226140/13)*common_65 + (1746880/39)*common_66 - 953325/13*common_96, 41325*common_1*potential*radial_speed*speed_squared + (102450/13)*common_11*common_34*potential + (46870/13)*common_34*potential - 462560/39*common_36 - 182535/13*common_37 - 9165*common_38 + (877270/13)*common_4*nu*radial_speed - 39770/13*common_43 - 109825/13*common_44 - 2921350/39*common_45 - 36160/13*common_47 + 820*potential*radial_speed*speed_squared, (2090/13)*common_10 - 128400/13*common_12 - 93450/13*common_13 + (75705/13)*common_14 + (8860/13)*common_32 - 52480/39*common_4 - 136925/13*common_48 - 40035/13*common_49 + (2025/13)*common_5 + (60340/3)*common_50 + (3795/13)*common_6 - 65080/13*common_7 - 22005/13*common_8 + (190310/13)*common_91 + (79275/13)*common_97 + (4490/13)*common_98),
        (-570*common_100 + 450*common_101 + 820*common_129 + 9450*common_130 - 6300*common_132 - 410*common_52 - 4725*common_54 - 225*common_55 + 1140*common_56 + 3150*common_57 + 4200*common_59 - 8400*common_60, 12600*common_0*common_15 - 6300*common_33 + 14700*common_35 - 9450*common_36 - 29400*common_38 + 18900*common_40, -18900*common_10 + 6440*common_4 - 7560*common_6 - 12880*common_7 + 37800*common_9*nu*potential + 3780*potential*speed_squared, -25200*common_31 + 12600*potential*radial_speed, 8400*common_22 - 4200*potential, 0),
        (1200*common_1*common_103 + 3922*common_103 - 7980*common_104 + (121275/13)*common_105*common_22 - 128205/26*common_106 - 22500/13*common_107*common_2 - 3600/13*common_107*potential + (176593/39)*common_108 + (11475/13)*common_109 - 474205/39*common_110 + (981075/13)*common_112 + (17998/13)*common_113 - 635152/13*common_114 + (130725/13)*common_115 - 2357475/13*common_116 - 62320/13*common_117 + (415800/13)*common_118 + (1228280/13)*common_123 + (555450/13)*common_125 - 415890/13*common_126 + (378000/13)*common_127 + (81375/26)*common_176 - 850500/13*common_177 + (30360/13)*common_178 - 453944/13*common_179 + (81810/13)*common_18*common_52 - 398220/13*common_18*common_59 + (1446385/13)*common_25*common_52 - 17325/13*common_32*common_53 - 142275/13*common_49*common_5, (452025/13)*common_71 + (1181610/13)*common_72 - 1009575/13*common_73 - 2720800/13*common_74 - 29025/13*common_75 - 183400/13*common_76 + (787255/13)*common_77 - 2950710/13*common_78 + (142425/13)*common_79 + (228550/13)*common_80 - 2215945/13*common_81 + (7137270/13)*common_82 - 916650/13*common_83 + (904920/13)*common_84 - 587250/13*common_87 + (1879500/13)*common_88 + (890490/13)*common_89 - 1997100/13*common_94, (2804220/13)*common_1*common_4*common_9 + (283770/13)*common_1*common_5*potential + (767550/13)*common_1*common_53*potential - 325960/13*common_100 - 48420/13*common_101 - 550960/13*common_102 - 641400/13*common_132 + (3403180/13)*common_4*common_9 + (1016800/13)*common_4*nu*speed_squared + (3330/13)*common_5*potential + (3922660/39)*common_52*nu - 1650880/39*common_52 + (2061500/13)*common_53*nu*potential - 904400/13*common_54 - 569800/13*common_58 - 8355520/13*common_60 + (370950/13)*common_9*potential*speed_squared - 1747800/13*common_96, (908280/13)*common_1*potential*radial_speed*speed_squared - 231540/13*common_33 + (838800/13)*common_34*potential - 1734520/13*common_36 - 1892700/13*common_38 + (4438180/13)*common_4*nu*radial_speed - 291600/13*common_44 - 1931460/13*common_45 + (376260/13)*nu*potential*radial_speed*speed_squared, -406200/13*common_10 + (33300/13)*common_14 + (38660/13)*common_32 + (233080/13)*common_4 - 274020/13*common_48 + (899400/13)*common_49 + (597840/13)*common_50 - 22040/13*common_6 - 737320/13*common_7, (97840/13)*common_15 + (17400/13)*common_17 - 218200/13*common_31),
        ((10933263/182)*common_0*common_103 - 1509255/13*common_0*common_106 - 4330205/52*common_0*common_112 - 62077595/91*common_0*common_114 + (127645/13)*common_1*common_180 - 12925/26*common_10*common_107 + (5619555/52)*common_100*common_53 + (225675/4)*common_101*common_53 - 16270165/273*common_103*common_18 + (15633248/819)*common_103*common_21 + (107775610/273)*common_103*common_24 - 171679175/273*common_103*common_25 - 69097832/819*common_103*common_69 + (198502582/819)*common_103*common_9 - 14728744/819*common_103*speed_squared + (4137045/13)*common_105*common_12 + (1268355/52)*common_105*common_32 + (383535/26)*common_105*common_4 + (2784375/26)*common_105*common_48 - 15797475/26*common_105*common_50 - 3184065/26*common_105*common_7 - 1771935/13*common_105*common_91 - 185778/13*common_107*common_12 + (2458125/26)*common_107*common_14 + (33382/13)*common_107*common_4 - 446925/26*common_107*common_49 + (100625/13)*common_107*common_50 - 2442*common_107*common_7 - 1629125/13*common_107*common_98 + (46366/39)*common_11*common_180 + (42607325/52)*common_112*common_18 + (124632325/182)*common_114*common_18 - 10131657/91*common_114*common_21 + (13428654/91)*common_114*speed_squared + (73779799/52)*common_129*common_53 + (606120/91)*common_13*common_52 + (3756825/13)*common_13*common_54 + (22548408/91)*common_13*common_59 + (1609175/13)*common_131*common_53 + (304713/13)*common_180*nu - 650021/39*common_180 - 39600*common_181*common_2 + (396495/4)*common_181*common_22 - 16500*common_181*common_23 - 189915/8*common_181*potential - 1102125/182*common_182*common_2 + (102525/91)*common_182*common_22 + (662825/91)*common_182*common_23 + (2075/91)*common_182*potential - 1220677/273*common_5*common_52 - 443475/104*common_5*common_54 - 5181269/91*common_5*common_59 - 21293215/91*common_5*common_66 - 23176895/52*common_52*common_53 + (35186545/1092)*common_52*common_8 - 599355/14*common_52*common_97 - 11998900/13*common_53*common_58 - 5081625/26*common_53*common_62 - 15209175/26*common_53*common_64 + (5820603/91)*common_59*common_8, (1406460/13)*common_136 - 288297694/819*common_137 - 10553895/26*common_138 + (39266497/42)*common_139 - 209575/91*common_140 - 3020445/26*common_141 - 60019394/273*common_142 - 16603065/52*common_143 + (15723818/13)*common_144 + (5347275/182)*common_145 + (10513845/26)*common_146 + (99452943/91)*common_147 + (4252755/4)*common_148 - 51695908/13*common_149 + (158400/13)*common_150 - 190772980/273*common_151 + (2526150/13)*common_152 + (146485370/819)*common_153 + (823275/26)*common_154 + (4845598/91)*common_155 - 1617783/26*common_156 - 3892275/26*common_157 - 772938/91*common_158 - 8439407/26*common_159 - 9711750/91*common_160 + (1581750/13)*common_161 - 353696255/273*common_162 + (15224000/91)*common_163 + (337500/13)*common_164 + (75938930/273)*common_165 + (60025425/52)*common_166 + (40444750/13)*common_167 - 17790675/26*common_168 - 7117580/13*common_169 + (2302500/13)*common_170 + (22560215/91)*common_171 - 36757945/26*common_172 - 6155250/13*common_173 - 34329240/91*common_174 + (16353435/13)*common_175, (70088560/273)*common_1*common_103 + (1075725/13)*common_1*common_105*potential + (3684210/91)*common_1*common_107*potential + (10981910/13)*common_1*common_4*common_9*speed_squared + (320475/13)*common_1*common_5*common_9*potential + (54159250/91)*common_1*common_52*speed_squared + (79570900/819)*common_103 - 23503930/91*common_104 + 607635*common_105*nu*potential - 4508385/26*common_106 + (183460/91)*common_107*potential - 1316220/91*common_109 + (1208120/7)*common_11*common_4*common_5 + 522690*common_11*common_4*common_53 + (4536000/13)*common_11*common_5*common_9*potential + (32365240/39)*common_11*common_52*common_9 + (3193800/13)*common_11*common_53*potential*speed_squared - 113315407/273*common_110 - 80486600/819*common_111 - 928752/91*common_113 - 43839676/39*common_114 - 5999625/13*common_115 - 27156255/13*common_116 - 2295768/91*common_117 - 4708200/13*common_120 - 6697280/91*common_121 - 15379660/91*common_122 - 7494465/13*common_125 - 154497130/39*common_126 - 14483070/13*common_128 - 39975*common_176 - 7485750/13*common_177 - 5298500/91*common_178 - 912734/13*common_179 + (9069945/13)*common_4*common_53 + 642198*common_4*common_9*nu*speed_squared + (1630650/13)*common_5*common_9*nu*potential + (52584767/13)*common_52*common_9*nu + (20500912/273)*common_52*speed_squared + (332175/2)*common_53*potential*speed_squared, (7958500/13)*common_1*common_34*potential*speed_squared + (87821560/39)*common_1*common_52*radial_speed + (7059150/13)*common_11*common_4*radial_speed*speed_squared + (3569550/13)*common_11*common_70*potential + (23249610/13)*common_34*common_4*nu + (2354400/13)*common_34*nu*potential*speed_squared + 51394*common_4*radial_speed*speed_squared + (231510/13)*common_5*potential*radial_speed + (1036096/3)*common_52*radial_speed + (1656505/13)*common_70*potential - 416910*common_73 - 21226396/13*common_74 - 1252150/13*common_76 - 7867030/13*common_78 - 532320/13*common_79 - 4332242/13*common_81 - 1456450/13*common_83 - 22348480/39*common_86 - 557430/13*common_87 - 2500190/13*common_89 - 2332350/13*common_90 - 3057000/13*common_92 - 3135410/13*common_94 - 1846920/13*common_95, -1244944/91*common_100 + (948810/91)*common_101 + (286010/91)*common_102 + (12737992/91)*common_129 + (1542900/13)*common_130 + (42099220/273)*common_131 + (336000/13)*common_132 + (5971860/91)*common_133 - 219300/13*common_134 + (696944/273)*common_52 - 41350*common_54 - 306955/91*common_55 + (5899214/91)*common_56 + (223200/13)*common_57 - 105919630/273*common_58 + (3015260/13)*common_59 - 9688560/13*common_60 + 72100*common_61 - 1210910/91*common_62 - 1256100/13*common_63 - 13447500/91*common_64 + (1104840/13)*common_65 + (4819720/13)*common_66 - 4053450/13*common_96, (49360/13)*common_33 + (40600/13)*common_35 - 20220*common_36 - 466800/13*common_37 - 28200/13*common_38 + (1460370/13)*common_40 + (1102100/13)*common_41 - 231300/13*common_43 - 325600/13*common_44 - 1796850/13*common_45 + (188700/13)*common_46 + (232800/13)*common_47),
    )

def i6(potential, radial_speed, speed_squared, nu, log_one, log_infrared):
    common_0 = nu**2
    common_1 = nu*potential
    common_2 = nu*speed_squared
    common_3 = common_0*potential
    common_4 = common_0*speed_squared
    common_5 = nu**3
    common_6 = common_5*potential
    common_7 = common_5*speed_squared
    common_8 = 12*radial_speed
    common_9 = common_0*radial_speed
    common_10 = 30*nu
    common_11 = 30*common_0
    common_12 = potential*radial_speed
    common_13 = common_3*radial_speed
    common_14 = common_6*radial_speed
    common_15 = (60/7)*potential
    common_16 = common_5*radial_speed
    common_17 = potential*speed_squared
    common_18 = potential**2
    common_19 = 180*nu
    common_20 = radial_speed**2
    common_21 = common_20*potential
    common_22 = common_3*speed_squared
    common_23 = common_18*common_5
    common_24 = common_20*common_6
    common_25 = common_1*radial_speed
    common_26 = common_12*speed_squared
    common_27 = radial_speed**3
    common_28 = common_18*radial_speed
    common_29 = common_1*speed_squared
    common_30 = common_29*radial_speed
    common_31 = common_1*common_27
    common_32 = common_28*nu
    common_33 = common_13*speed_squared
    common_34 = common_14*speed_squared
    common_35 = common_18*common_9
    common_36 = common_27*common_6
    common_37 = common_16*common_18
    common_38 = common_18*nu
    common_39 = 360*potential
    common_40 = 4800*potential
    common_41 = 360*nu
    common_42 = 360*common_5
    common_43 = 720*common_0 - common_41 - common_42 + 360/7
    common_44 = 450*common_20
    common_45 = common_0*common_18
    common_46 = potential**3
    common_47 = radial_speed**4
    common_48 = common_47*potential
    common_49 = speed_squared**2
    common_50 = common_1*common_49
    common_51 = common_18*common_2
    common_52 = common_17*common_20
    common_53 = common_0*common_46
    common_54 = common_20*common_38
    common_55 = common_18*common_4
    common_56 = common_47*common_6
    common_57 = common_49*common_6
    common_58 = common_20*common_22
    common_59 = common_20*common_23
    common_60 = common_27*potential
    common_61 = common_27*common_3
    common_62 = (5400/7)*speed_squared
    common_63 = common_20*common_3
    common_64 = (1980/7)*speed_squared
    common_65 = nu*radial_speed
    common_66 = 3600*common_0
    common_67 = common_1*common_20
    common_68 = common_66 - 3600*nu + 720
    common_69 = radial_speed**5
    common_70 = common_69*potential
    common_71 = common_46*radial_speed
    common_72 = common_1*common_69
    common_73 = common_46*nu
    common_74 = common_73*radial_speed
    common_75 = common_12*common_49
    common_76 = common_17*common_27
    common_77 = common_28*speed_squared
    common_78 = common_18*common_27
    common_79 = common_25*common_49
    common_80 = common_27*common_29
    common_81 = common_2*common_28
    common_82 = common_27*common_38
    common_83 = common_46*common_9
    common_84 = common_6*common_69
    common_85 = common_16*common_46
    common_86 = common_13*common_49
    common_87 = common_22*common_27
    common_88 = common_28*common_4
    common_89 = common_14*common_49
    common_90 = common_6*speed_squared
    common_91 = common_27*common_90
    common_92 = common_28*common_7
    common_93 = common_27*common_45
    common_94 = common_23*common_27
    common_95 = common_18*speed_squared
    common_96 = common_3*common_47
    common_97 = common_49*potential
    common_98 = common_20*common_29
    common_99 = common_3*common_49
    common_100 = potential**4
    common_101 = common_100*nu
    common_102 = radial_speed**6
    common_103 = common_102*potential
    common_104 = speed_squared**3
    common_105 = common_104*potential
    common_106 = common_2*common_46
    common_107 = common_100*common_5
    common_108 = common_20*common_46
    common_109 = common_38*common_47
    common_110 = common_38*common_49
    common_111 = common_104*common_3
    common_112 = common_102*common_6
    common_113 = common_46*common_7
    common_114 = common_20*common_95
    common_115 = common_20*common_50
    common_116 = common_22*common_47
    common_117 = common_20*common_53
    common_118 = common_23*common_47
    common_119 = common_23*common_49
    common_120 = common_18*common_20
    common_121 = common_120*common_4
    common_122 = common_24*common_49
    common_123 = common_3*common_69
    common_124 = 151200*common_12
    common_125 = 756000*common_13
    common_126 = radial_speed**7
    common_127 = common_46*common_5
    return (
        (5*common_0 - 5*nu + 1, 0, 0, 0, 0, 0, 0),
        (9*common_1 - 21/2*common_2 - 27*common_3 + 33*common_4 + 26*common_6 - 63/2*common_7 - potential + (15/14)*speed_squared, common_5*common_8 + common_8*nu - 24*common_9 - 12/7*radial_speed, 30*common_0 - 15*common_5 - 15*nu + 15/7, 0, 0, 0, 0),
        (0, -common_10 + common_11 + 6, 0, 0, 0, 0, 0),
        ((4/7)*common_12 - 15*common_13 + 25*common_14, 72*common_1 - common_15 - 51*common_2 - 198*common_3 + 174*common_4 + 174*common_6 - 177*common_7 + (33/7)*speed_squared, common_10*radial_speed + 30*common_16 - 60*common_9 - 30/7*radial_speed, 120*common_0 - 60*common_5 - 60*nu + 60/7, 0, 0, 0),
        (-common_11*potential + 30*nu*potential - 6*potential, 0, 150*common_0 - 150*nu + 30, 0, 0, 0, 0),
        (273*common_0*common_18 - 71/7*common_17 - common_18*common_19 + 32*common_18 + 9*common_20*nu*potential - 12/7*common_21 - 129*common_22 - 139*common_23 - 30*common_24 + 112*common_5*potential*speed_squared + 63*nu*potential*speed_squared, (246/7)*common_12 + 60*common_13 + 210*common_14 - 162*common_25, 510*common_1 - 225*common_2 - 1290*common_3 + 810*common_4 + 1020*common_6 - 855*common_7 - 450/7*potential + (135/7)*speed_squared, 0, 360*common_0 - common_19 - 180*common_5 + 180/7, 0, 0),
        (18*common_12 + 90*common_13 - 90*common_25, -480*common_3 + 480*nu*potential - 96*potential, 0, 600*common_0 - 600*nu + 120, 0, 0, 0),
        (common_15*common_27 + 27*common_26 - 976/7*common_28 - 171*common_30 - 45*common_31 + 738*common_32 + 387*common_33 - 396*common_34 - 894*common_35 + 150*common_36 + 182*common_37, 4758*common_0*common_18 - 870/7*common_17 + (3678/7)*common_18 - 630*common_20*common_3 - common_20*common_39*common_5 + 630*common_20*nu*potential - 810/7*common_21 - 1734*common_22 - 2484*common_23 - 3018*common_38 + 1692*common_5*potential*speed_squared + 786*nu*potential*speed_squared, (3090/7)*common_12 + 1770*common_13 + 1140*common_14 - 2190*common_25, 2760*common_1 - 900*common_2 - 6600*common_3 - common_39 + 3240*common_4 + common_40*common_5 - 3420*common_7 + (540/7)*speed_squared, 720*common_0*radial_speed - common_41*radial_speed - common_42*radial_speed + (360/7)*radial_speed, common_43, 0),
        (common_1*common_44 + 18*common_17 + 78*common_18 - 90*common_21 + 90*common_22 - 90*common_29 - common_3*common_44 - 390*common_38 + 390*common_45, 396*common_12 + 1980*common_13 - 1980*common_25, -4200*common_3 + 4200*nu*potential - 840*potential, 0, 1800*common_0 - 1800*nu + 360, 0, 0),
        (6345*common_0*common_18*common_20 + 387*common_0*common_49*potential + (6612/7)*common_18*common_20 + 56*common_18*common_5*speed_squared + (251/7)*common_18*speed_squared + 2430*common_20*common_5*potential*speed_squared + 720*common_20*nu*potential*speed_squared + 1522*common_46*common_5 + 3684*common_46*nu - 698*common_46 + 315*common_47*nu*potential - 60*common_48 + 27*common_49*potential - 171*common_50 - 33*common_51 - 765/7*common_52 - 4644*common_53 - 5094*common_54 - 327*common_55 - 1050*common_56 - 396*common_57 - 1935*common_58 - 1155*common_59, (2124/7)*common_26 - 26076/7*common_28 - 2124*common_30 - 3420*common_31 + 20160*common_32 + 6264*common_33 - 8172*common_34 - 26088*common_35 + 2700*common_36 + 7044*common_37 + 630*common_60 + 3150*common_61, 45420*common_0*common_18 + (32940/7)*common_18 + 10260*common_20*nu*potential - 13320/7*common_21 - 13020*common_22 - 24360*common_23 - 2520*common_24 - 27660*common_38 + 14460*common_5*potential*speed_squared - common_62*potential - 11160*common_63 + 5160*nu*potential*speed_squared, (22440/7)*common_12 + 15120*common_13 + 4440*common_14 - 16200*common_25, 10440*common_1 - 3060*common_2 - 24120*common_3 + 10440*common_4 + 16560*common_6 + common_64 - 10620*common_7 - 9720/7*potential, 2880*common_0*radial_speed - 1440*common_16 - 1440*common_65 + (1440/7)*radial_speed, common_43),
        (3150*common_0*common_27*potential + 3600*common_18*nu*radial_speed - 270*common_26 + 630*common_27*potential - common_28*common_66 - 720*common_28 - 3150*common_31 - 1350*common_33 + 1350*nu*potential*radial_speed*speed_squared, 504*common_17 + 1752*common_18 - 2520*common_21 + 2520*common_22 - 2520*common_29 - 8760*common_38 + 8760*common_45 - 12600*common_63 + 12600*common_67, 4500*common_12 + 22500*common_13 - 22500*common_25, -24000*common_3 - common_40 + 24000*nu*potential, 0, common_68, 0),
        (540*common_70 + (53608/7)*common_71 - 2835*common_72 - 40230*common_74 - 2475/7*common_75 + 525*common_76 + (9612/7)*common_77 - 8448*common_78 + 2295*common_79 - 3780*common_80 - 8910*common_81 + 45810*common_82 + 49080*common_83 + 9450*common_84 - 12980*common_85 - 5805*common_86 + 13545*common_87 + 14310*common_88 + 6840*common_89 - 21210*common_91 - 1890*common_92 - 57330*common_93 + 9570*common_94, 228330*common_0*common_18*common_20 + 8586*common_0*common_49*potential + (227466/7)*common_18*common_20 + 5892*common_18*common_5*speed_squared + 9198*common_18*nu*speed_squared + 63540*common_20*common_5*potential*speed_squared + 4680*common_20*nu*potential*speed_squared - common_21*common_64 + 37968*common_46*common_5 + 80376*common_46*nu - 105072/7*common_46 + 25830*common_47*nu*potential - 4770*common_48 + (3258/7)*common_49*potential - 3150*common_50 - 105456*common_53 - 176886*common_54 - 21114*common_55 - 25200*common_56 - 10548*common_57 - 33480*common_58 - 53460*common_59 - 6906/7*common_95 - 22050*common_96, (180/7)*common_26 - 339060/7*common_28 - 5580*common_30 - 68400*common_31 + 266160*common_32 + 48060*common_33 - 89280*common_34 - 359520*common_35 + 26100*common_36 + 112260*common_37 + (88650/7)*common_60 + 71550*common_61, 277440*common_0*common_18 - 17160/7*common_17 + (188520/7)*common_18 + 91800*common_20*nu*potential - 120600/7*common_21 - 64320*common_22 - 153720*common_23 - 12600*common_24 - 162120*common_38 + 83160*common_5*potential*speed_squared - 100800*common_63 + 19560*nu*potential*speed_squared, (106200/7)*common_12 + 77400*common_13 + 10800*common_14 - 77400*common_25, 24480*common_1 - 7560*common_2 - 55440*common_3 + 23760*common_4 + 36720*common_6 + common_62 - 22680*common_7 - 23040/7*potential, 4320*common_0*radial_speed - 2160*common_16 - 2160*common_65 + (2160/7)*radial_speed),
        (34650*common_0*common_18*common_20 + 18900*common_0*common_20*potential*speed_squared + 6930*common_18*common_20 + 4770*common_18*nu*speed_squared + 3780*common_20*potential*speed_squared + 5160*common_46*nu - 1032*common_46 + 28350*common_47*nu*potential - 5670*common_48 + 1350*common_49*nu*potential - 5160*common_53 - 34650*common_54 - 4770*common_55 - 954*common_95 - 28350*common_96 - 270*common_97 - 18900*common_98 - 1350*common_99, 107100*common_0*common_27*potential + 99000*common_18*nu*radial_speed - 9180*common_26 + 21420*common_27*potential - 19800*common_28 - 107100*common_31 - 45900*common_33 - 99000*common_35 + 45900*nu*potential*radial_speed*speed_squared, 7020*common_17 + 18660*common_18 - 35100*common_21 + 35100*common_22 - 35100*common_29 - 93300*common_38 + 93300*common_45 - 175500*common_63 + 175500*common_67, 32400*common_12 + 162000*common_13 - 162000*common_25, -90000*common_3 + 90000*nu*potential - 18000*potential, 0, common_68),
        (66696*common_0*common_100 + 633735*common_0*common_18*common_47 + 13869*common_0*common_18*common_49 + 81270*common_0*common_20*common_49*potential + 68544*common_0*common_46*speed_squared + 11480*common_100 - 58722*common_101 + 31185*common_102*nu*potential - 5940*common_103 + 6840*common_104*common_5*potential + 2295*common_104*nu*potential - 2475/7*common_105 - 42666*common_106 - 14668*common_107 - 86382*common_108 - 505530*common_109 - 7587*common_110 - 5805*common_111 - 103950*common_112 - 30122*common_113 - 38277*common_114 - 27405*common_115 - 121905*common_116 - 567720*common_117 - 103005*common_118 - 1692*common_119 - 308160*common_121 - 111510*common_122 - 2025*common_17*common_47 + 43965*common_18*common_20*common_5*speed_squared + 222021*common_18*common_20*nu*speed_squared + 93060*common_18*common_47 + (7191/7)*common_18*common_49 + 148590*common_20*common_46*common_5 + 457794*common_20*common_46*nu + 4050*common_20*common_49*potential + (50230/7)*common_46*speed_squared + 238140*common_47*common_5*potential*speed_squared + 19845*common_47*nu*potential*speed_squared, 198450*common_123 + 46170*common_70 + (1482312/7)*common_71 - 249480*common_72 - 1128924*common_74 - 35100/7*common_75 - 13950*common_76 + (552168/7)*common_77 - 346572*common_78 + 38880*common_79 + 47880*common_80 - 471348*common_81 + 1887840*common_82 + 1439520*common_83 + 283500*common_84 - 440880*common_85 - 144720*common_86 + 227430*common_87 + 700740*common_88 + 220860*common_89 - 672840*common_91 - 149580*common_92 - 2422620*common_93 + 519480*common_94, 3769650*common_0*common_18*common_20 + 90990*common_0*common_49*potential + 520740*common_18*common_20 + 130020*common_18*common_5*speed_squared + 239250*common_18*nu*speed_squared + 842400*common_20*common_5*potential*speed_squared + 36450*common_20*potential*speed_squared + 439740*common_46*common_5 + 800280*common_46*nu - 1029060/7*common_46 + 607950*common_47*nu*potential - 112500*common_48 + (16470/7)*common_49*potential - 21330*common_50 - 1099080*common_53 - 2852280*common_54 - 419190*common_55 - 308700*common_56 - 142020*common_57 - 193050*common_58 - 1012050*common_59 - 252990/7*common_95 - 611100*common_96 - 153900*common_98, -27000*common_26 - 380400*common_28 + 102600*common_30 - 732600*common_31 + 2116800*common_32 + 183600*common_33 - 631800*common_34 - 2959200*common_35 + 167400*common_36 + 1023600*common_37 + 136800*common_60 + 790200*common_61, 1068120*common_0*common_18 + 96480*common_18 + 513000*common_20*nu*potential - 97200*common_21 - 198360*common_22 - 619560*common_23 - 43200*common_24 - 596160*common_38 + 319680*common_5*potential*speed_squared - 561600*common_63 + 26280*nu*potential*speed_squared + 360*potential*speed_squared, 44640*common_12 + 244800*common_13 + 7200*common_14 - 230400*common_25, 26640*common_1 - 9720*common_2 - 59760*common_3 + 28080*common_4 + 38880*common_6 - 24840*common_7 - 3600*potential + 1080*speed_squared),
        (311850*common_123 + 62370*common_70 + 17136*common_71 - 311850*common_72 - 85680*common_74 + 9450*common_75 - 56700*common_76 + 24192*common_77 - 75600*common_78 - 47250*common_79 + 283500*common_80 - 120960*common_81 + 378000*common_82 + 85680*common_83 + 47250*common_86 - 283500*common_87 + 120960*common_88 - 378000*common_93, 1121400*common_0*common_18*common_20 + 756000*common_0*common_20*potential*speed_squared + 224280*common_18*common_20 + 151920*common_18*nu*speed_squared + 151200*common_20*potential*speed_squared + 118560*common_46*nu - 23712*common_46 + 1134000*common_47*nu*potential - 226800*common_48 + 54000*common_49*nu*potential - 118560*common_53 - 1121400*common_54 - 151920*common_55 - 30384*common_95 - 1134000*common_96 - 10800*common_97 - 756000*common_98 - 54000*common_99, 1764000*common_0*common_27*potential - common_124*speed_squared - common_125*speed_squared + 1260000*common_18*nu*radial_speed + 352800*common_27*potential - 252000*common_28 - 1764000*common_31 - 1260000*common_35 + 756000*nu*potential*radial_speed*speed_squared, 60480*common_17 + 114240*common_18 - 302400*common_21 + 302400*common_22 - 302400*common_29 - 571200*common_38 + 571200*common_45 - 1512000*common_63 + 1512000*common_67, common_124 + common_125 - 756000*common_25, -201600*common_3 + 201600*nu*potential - 40320*potential, 0),
        (203175*common_0*common_104*potential*radial_speed + 6088500*common_0*common_18*common_27*speed_squared + 7286580*common_0*common_27*common_46 + 1340955*common_0*common_69*potential*speed_squared - 405405*common_1*common_126 + 272144*common_100*common_5*radial_speed - 1186608*common_100*common_9 + 1027152*common_100*nu*radial_speed - 1396768/7*common_100*radial_speed - 270900*common_104*common_14 - 70875*common_104*common_25 + 10575*common_104*potential*radial_speed + 1351350*common_126*common_5*potential + 77220*common_126*potential - 1964340*common_127*common_27 - 13365*common_17*common_69 + 815040*common_18*common_27*speed_squared + 113886*common_18*common_49*common_5*radial_speed - 769662*common_18*common_49*common_9 + 525258*common_18*common_49*nu*radial_speed + 1374030*common_18*common_5*common_69 + 6597450*common_18*common_69*nu - 1213920*common_18*common_69 - 4584060*common_2*common_78 + 1081728*common_27*common_46 + 1956150*common_27*common_49*common_5*potential + 326025*common_27*common_49*nu*potential - 5775588*common_27*common_73 - 44550*common_27*common_97 - 1219050*common_27*common_99 - 611496/7*common_28*common_49 - 31185*common_29*common_69 - 2054952*common_4*common_71 - 8287110*common_45*common_69 + 746496*common_46*common_5*radial_speed*speed_squared + 1432116*common_46*nu*radial_speed*speed_squared - 3243240*common_69*common_90 - 951300*common_7*common_78 - 1780560/7*common_71*speed_squared, 1334016*common_0*common_100 + 30556350*common_0*common_18*common_47 + 795294*common_0*common_18*common_49 + 2182950*common_0*common_20*common_49*potential + 2212224*common_0*common_46*speed_squared + (1548528/7)*common_100 - 1139328*common_101 - 2182950*common_102*common_3 + 2931390*common_102*nu*potential - 543510*common_103 + 261900*common_104*common_5*potential + 52650*common_104*nu*potential - 49950/7*common_105 - 1444512*common_106 - 351408*common_107 - 2839248*common_108 - 23897250*common_109 - 503370*common_110 - 179550*common_111 - 3742200*common_112 - 946272*common_113 - 2034000*common_114 - 292950*common_115 - 1786050*common_116 - 19597320*common_117 - 6306300*common_118 - 169452*common_119 - 15821100*common_121 - 4233600*common_122 + 3237300*common_18*common_20*common_5*speed_squared + 11541600*common_18*common_20*nu*speed_squared + 4388850*common_18*common_47 + (563454/7)*common_18*common_49 + 5908680*common_20*common_46*common_5 + 15225696*common_20*common_46*nu + 17550*common_20*common_49*potential - 1559250*common_29*common_47 + (1750944/7)*common_46*speed_squared + 8901900*common_47*common_5*potential*speed_squared + 344250*common_47*potential*speed_squared, 6492150*common_123 + 1243350*common_70 + (18013080/7)*common_71 - 6718950*common_72 - 13930680*common_74 + (252450/7)*common_75 - 774900*common_76 + (12058920/7)*common_77 - 6493680*common_78 - 6750*common_79 + 3748500*common_80 - 9947160*common_81 + 35481600*common_82 + 18597840*common_83 + 4195800*common_84 - 6552120*common_85 - 1564650*common_86 + 44100*common_87 + 14282280*common_88 + 3499200*common_89 - 10495800*common_91 - 3669840*common_92 - 46330200*common_93 + 11453400*common_94, 7560000*common_1*common_47 + 4916160*common_120 + 2817600*common_127 + 1490760*common_18*common_7 + 36527400*common_20*common_45 + 7030800*common_24*speed_squared - 5263680/7*common_46 - 1407600*common_48 + 17280*common_50 + 2786880*common_51 + 691200*common_52 - 6105600*common_53 - 27079920*common_54 - 4422120*common_55 - 2406600*common_56 - 1199880*common_57 + 680400*common_58 - 10848600*common_59 + 4181760*common_73 - 3218640/7*common_95 - 7975800*common_96 - 123120/7*common_97 - 3326400*common_98 + 547560*common_99, -276480*common_26 - 1844640*common_28 + 1254960*common_30 - 4762800*common_31 + 10392480*common_32 + 22680*common_33 - 2940840*common_34 - 15014160*common_35 + 718200*common_36 + 5680080*common_37 + 896400*common_60 + 5178600*common_61, 2293200*common_0*common_18 + 187200*common_18 + 1723680*common_20*nu*potential - 328320*common_21 - 277200*common_22 - 1436400*common_23 - 75600*common_24 - 110880*common_29 - 1199520*common_38 + 730800*common_5*potential*speed_squared - 1890000*common_63 + 34560*potential*speed_squared, 378000*common_0*potential*radial_speed - 25200*common_14 - 332640*common_25 + 63360*potential*radial_speed),
    )

def i7(potential, radial_speed, speed_squared, nu, log_one, log_infrared):
    common_0 = nu**2
    common_1 = nu*potential
    common_2 = common_0*speed_squared
    common_3 = nu**3
    common_4 = common_3*potential
    common_5 = nu*radial_speed
    common_6 = common_3*radial_speed
    common_7 = 21*common_0
    common_8 = potential*radial_speed
    common_9 = common_1*radial_speed
    common_10 = common_4*radial_speed
    common_11 = (1050/17)*common_0*radial_speed - 630/17*common_5 - 420/17*common_6 + (105/17)*radial_speed
    common_12 = potential*speed_squared
    common_13 = potential**2
    common_14 = common_1*speed_squared
    common_15 = common_13*nu
    common_16 = radial_speed**2
    common_17 = common_16*potential
    common_18 = common_1*common_16
    common_19 = common_2*potential
    common_20 = common_4*speed_squared
    common_21 = common_0*common_13
    common_22 = common_13*common_3
    common_23 = common_0*potential
    common_24 = common_16*common_23
    common_25 = common_16*common_4
    common_26 = common_23*radial_speed
    common_27 = common_8*speed_squared
    common_28 = radial_speed**3
    common_29 = common_28*potential
    common_30 = common_13*common_5
    common_31 = common_2*common_8
    common_32 = common_28*common_4
    common_33 = common_13*common_6
    common_34 = common_0*radial_speed
    common_35 = potential**3
    common_36 = common_35*nu
    common_37 = radial_speed**4
    common_38 = speed_squared**2
    common_39 = common_38*potential
    common_40 = common_13*speed_squared
    common_41 = common_1*common_37
    common_42 = common_3*common_35
    common_43 = common_13*common_16
    common_44 = common_14*common_16
    common_45 = common_23*common_37
    common_46 = common_23*common_38
    common_47 = common_13*common_2
    common_48 = common_16*common_20
    common_49 = common_16*common_21
    common_50 = common_23*common_28
    common_51 = -50400/17*common_0 + (20160/17)*common_3 + (30240/17)*nu - 5040/17
    common_52 = common_13*radial_speed
    common_53 = common_9*speed_squared
    common_54 = common_1*common_28
    common_55 = common_13*common_34
    common_56 = radial_speed**5
    common_57 = common_56*potential
    common_58 = common_35*radial_speed
    common_59 = common_52*speed_squared
    common_60 = common_38*common_9
    common_61 = common_15*common_28
    common_62 = common_34*common_35
    common_63 = common_4*common_56
    common_64 = common_2*common_29
    common_65 = common_2*common_52
    common_66 = common_10*common_38
    common_67 = common_3*speed_squared
    common_68 = common_52*common_67
    common_69 = common_22*common_28
    common_70 = common_37*potential
    common_71 = common_1*common_38
    common_72 = nu*speed_squared
    common_73 = common_13*common_72
    common_74 = common_17*speed_squared
    common_75 = common_0*common_35
    common_76 = common_15*common_16
    common_77 = common_37*common_4
    common_78 = common_38*common_4
    common_79 = common_13*common_67
    common_80 = common_17*common_2
    common_81 = common_16*common_22
    common_82 = -15120*common_0 + 20160*nu - 5040
    common_83 = potential**4
    common_84 = common_83*nu
    common_85 = radial_speed**6
    common_86 = common_85*potential
    common_87 = speed_squared**3
    common_88 = common_87*potential
    common_89 = common_35*speed_squared
    common_90 = common_1*common_85
    common_91 = common_1*common_87
    common_92 = common_35*common_72
    common_93 = common_12*common_37
    common_94 = common_0*common_83
    common_95 = common_3*common_83
    common_96 = common_13*common_37
    common_97 = common_13*common_38
    common_98 = common_16*common_35
    common_99 = common_14*common_37
    common_100 = common_15*common_37
    common_101 = common_15*common_38
    common_102 = common_16*common_36
    common_103 = common_23*common_85
    common_104 = common_23*common_87
    common_105 = common_2*common_35
    common_106 = common_4*common_85
    common_107 = common_4*common_87
    common_108 = common_35*common_67
    common_109 = common_17*common_38
    common_110 = common_43*speed_squared
    common_111 = common_18*common_38
    common_112 = common_43*common_72
    common_113 = common_2*common_70
    common_114 = common_20*common_37
    common_115 = common_21*common_37
    common_116 = common_21*common_38
    common_117 = common_16*common_75
    common_118 = common_22*common_37
    common_119 = common_22*common_38
    common_120 = common_16*common_42
    common_121 = common_24*common_38
    common_122 = common_2*common_43
    common_123 = common_25*common_38
    common_124 = common_43*common_67
    common_125 = common_14*common_28
    common_126 = common_23*common_56
    common_127 = common_38*common_8
    common_128 = common_26*common_38
    common_129 = radial_speed**7
    common_130 = common_129*potential
    common_131 = common_5*common_83
    common_132 = common_8*common_87
    common_133 = common_28*common_35
    common_134 = common_14*common_56
    common_135 = common_58*common_72
    common_136 = common_15*common_56
    common_137 = common_129*common_23
    common_138 = common_129*common_4
    common_139 = common_6*common_83
    common_140 = common_29*common_38
    common_141 = common_28*common_40
    common_142 = common_38*common_54
    common_143 = common_5*common_97
    common_144 = common_13*common_28
    common_145 = common_26*common_87
    common_146 = common_2*common_57
    common_147 = common_58*common_67
    common_148 = common_28*common_75
    common_149 = common_22*common_56
    common_150 = common_144*common_2
    common_151 = common_32*common_38
    common_152 = potential**5
    common_153 = radial_speed**8
    common_154 = speed_squared**4
    return (
        (-3*common_0 + 4*nu - 1, 0, 0, 0, 0, 0, 0, 0),
        ((809/34)*common_0*potential - 295/34*common_1 - 975/34*common_2 + (747/34)*common_3*speed_squared - 633/34*common_4 + (347/34)*nu*speed_squared + (35/34)*potential - 19/17*speed_squared, (315/17)*common_0*radial_speed - 189/17*common_5 - 126/17*common_6 + (63/34)*radial_speed, -420/17*common_0 + (168/17)*common_3 + (252/17)*nu - 42/17, 0, 0, 0, 0, 0),
        (0, -common_7 + 28*nu - 7, 0, 0, 0, 0, 0, 0),
        ((511/34)*common_0*potential*radial_speed - 609/34*common_10 - 11/17*common_8 - 21/34*common_9, (6713/34)*common_0*potential - 2695/34*common_1 - 6195/34*common_2 + (4977/34)*common_3*speed_squared - 4851/34*common_4 + (2051/34)*nu*speed_squared + (175/17)*potential - 203/34*speed_squared, common_11, -2100/17*common_0 + (840/17)*common_3 + (1260/17)*nu - 210/17, 0, 0, 0, 0),
        (common_7*potential - 28*nu*potential + 7*potential, 0, -126*common_0 + 168*nu - 42, 0, 0, 0, 0, 0),
        ((419/34)*common_12 - 640/17*common_13 - 1155/17*common_14 + (3024/17)*common_15 + (33/17)*common_17 - 147/17*common_18 + (2282/17)*common_19 - 1722/17*common_20 - 4088/17*common_21 + (2016/17)*common_22 - 105/34*common_24 + (378/17)*common_25, (147/17)*common_0*potential*radial_speed - 3129/17*common_10 - 735/17*common_8 + (2625/17)*nu*potential*radial_speed, (25389/17)*common_0*potential - 11235/17*common_1 - 17535/17*common_2 + (14511/17)*common_3*speed_squared - 16653/17*common_4 + (5523/17)*nu*speed_squared + (1575/17)*potential - 504/17*speed_squared, common_11, -8400/17*common_0 + (3360/17)*common_3 + (5040/17)*nu - 840/17, 0, 0, 0),
        (-63*common_26 - 21*common_8 + 84*nu*potential*radial_speed, -532*common_1 + 399*common_23 + 133*potential, 0, -630*common_0 + 840*nu - 210, 0, 0, 0, 0),
        ((11746/17)*common_0*common_13*radial_speed + (525/34)*common_0*common_28*potential + (2810/17)*common_13*radial_speed - 1125/34*common_27 + (735/17)*common_28*nu*potential - 165/17*common_29 + (5922/17)*common_3*potential*radial_speed*speed_squared - 12117/17*common_30 - 6951/17*common_31 - 1890/17*common_32 - 2247/17*common_33 + (3171/17)*nu*potential*radial_speed*speed_squared, (6335/34)*common_12 - 12607/17*common_13 - 17934/17*common_14 + (61005/17)*common_15 + (2436/17)*common_17 - 11046/17*common_18 + (38339/17)*common_19 - 31353/17*common_20 - 85253/17*common_21 + (41979/17)*common_22 + (15519/34)*common_24 + (5607/17)*common_25, -21315/17*common_10 - 21147/17*common_26 - 11004/17*common_8 + (43533/17)*nu*potential*radial_speed, (159495/17)*common_0*potential - 75705/17*common_1 - 86625/17*common_2 + (72135/17)*common_3*speed_squared - 96285/17*common_4 + (26985/17)*nu*speed_squared + (11130/17)*potential - 2415/17*speed_squared, -12600/17*common_34 + (7560/17)*common_5 + (5040/17)*common_6 - 1260/17*radial_speed, -25200/17*common_0 + (10080/17)*common_3 + (15120/17)*nu - 2520/17, 0, 0),
        (315*common_0*common_16*potential - 21*common_12 + 448*common_13*nu - 112*common_13 + 105*common_16*potential - 420*common_18 - 63*common_19 - 336*common_21 + 84*nu*potential*speed_squared, -1638*common_26 - 546*common_8 + 2184*nu*potential*radial_speed, -5712*common_1 + 4284*common_23 + 1428*potential, 0, -2520*common_0 + 3360*nu - 840, 0, 0, 0),
        ((71085/34)*common_0*common_16*potential*speed_squared + (81123/17)*common_0*common_35 + (12411/17)*common_13*common_16*common_3 + (84483/17)*common_13*common_16*nu + (6048/17)*common_13*common_3*speed_squared + (4550/17)*common_13*nu*speed_squared + (4635/34)*common_16*potential*speed_squared + (13230/17)*common_3*common_37*potential + (5922/17)*common_3*common_38*potential + (17413/17)*common_35 - 75544/17*common_36 + (1155/17)*common_37*potential + (3171/17)*common_38*nu*potential - 1125/34*common_39 - 1699/17*common_40 - 5145/17*common_41 - 28308/17*common_42 - 19104/17*common_43 - 13650/17*common_44 - 3675/34*common_45 - 6951/17*common_46 - 2506/17*common_47 - 35280/17*common_48 - 83685/17*common_49, (418551/17)*common_0*common_13*radial_speed + (89943/17)*common_13*radial_speed - 504*common_27 + (60375/17)*common_28*nu*potential - 13335/17*common_29 + 8631*common_3*potential*radial_speed*speed_squared - 398643/17*common_30 - 8715*common_31 - 41265/17*common_32 - 103089/17*common_33 - 36960/17*common_50 + 3171*nu*potential*radial_speed*speed_squared, (25956/17)*common_12 - 140868/17*common_13 - 155736/17*common_14 + (699552/17)*common_15 + (47628/17)*common_17 - 212940/17*common_18 + (372372/17)*common_19 - 329448/17*common_20 - 1011696/17*common_21 + (497784/17)*common_22 + (174258/17)*common_24 + (49392/17)*common_25, -110460/17*common_10 - 255780/17*common_26 - 99120/17*common_8 + (401940/17)*nu*potential*radial_speed, (776580/17)*common_0*potential - 385980/17*common_1 - 359100/17*common_2 + (293580/17)*common_3*speed_squared - 440580/17*common_4 + (115500/17)*nu*speed_squared + (58380/17)*potential - 10920/17*speed_squared, -88200/17*common_34 + (52920/17)*common_5 + (35280/17)*common_6 - 8820/17*radial_speed, common_51, 0),
        (315*common_27 - 735*common_29 - 4200*common_30 + 945*common_31 - 2205*common_50 + 1050*common_52 - 1260*common_53 + 2940*common_54 + 3150*common_55, 10395*common_0*common_16*potential - 693*common_12 + 12376*common_13*nu - 3094*common_13 + 3465*common_16*potential - 13860*common_18 - 2079*common_19 - 9282*common_21 + 2772*nu*potential*speed_squared, -22680*common_26 - 7560*common_8 + 30240*nu*potential*radial_speed, -42000*common_1 + 31500*common_23 + 10500*potential, 0, -7560*common_0 + 10080*nu - 2520, 0, 0),
        ((752661/17)*common_0*common_13*common_28 + (105840/17)*common_0*common_38*potential*radial_speed + (33075/34)*common_0*common_56*potential - 1365/2*common_12*common_28 + (171192/17)*common_13*common_28 + (108801/17)*common_13*nu*radial_speed*speed_squared + 17640*common_28*common_3*potential*speed_squared + 4410*common_28*nu*potential*speed_squared + (225771/17)*common_3*common_35*radial_speed + (839349/17)*common_35*nu*radial_speed + (14895/34)*common_38*potential*radial_speed + (46305/17)*common_56*nu*potential - 10395/17*common_57 - 195932/17*common_58 - 18261/17*common_59 - 43155/17*common_60 - 762258/17*common_61 - 839734/17*common_62 - 119070/17*common_63 - 30135/2*common_64 - 119574/17*common_65 - 100170/17*common_66 - 49203/17*common_68 - 90720/17*common_69, 27692*common_35 - 122297*common_36 - 25011/34*common_39 - 4564/17*common_40 - 458640/17*common_41 - 50295*common_42 - 46221*common_43 - 183960/17*common_44 + (491715/34)*common_45 - 196812/17*common_46 + (157850/17)*common_47 - 1104390/17*common_48 - 214431*common_49 + (101430/17)*common_70 + (76104/17)*common_71 - 69958/17*common_73 + (38115/34)*common_74 + 137046*common_75 + 206661*common_76 + (381465/17)*common_77 + (188181/17)*common_78 + (97986/17)*common_79 + (1757385/34)*common_80 + 42651*common_81, (7203210/17)*common_0*common_13*radial_speed + (1447530/17)*common_13*radial_speed - 34020/17*common_27 + (1426950/17)*common_28*nu*potential - 318150/17*common_29 + (1967490/17)*common_3*potential*radial_speed*speed_squared - 6532470/17*common_30 - 1657530/17*common_31 - 494550/17*common_32 - 2009490/17*common_33 - 1093050/17*common_50 + (364770/17)*nu*potential*radial_speed*speed_squared, (131460/17)*common_12 - 1067220/17*common_13 - 895860/17*common_14 + (5439420/17)*common_15 + 31500*common_17 - 138600*common_18 + (2528400/17)*common_19 - 2417940/17*common_20 - 8137500/17*common_21 + (4018980/17)*common_22 + 116550*common_24 + 18900*common_25, -422100/17*common_10 - 1782900/17*common_26 - 617400/17*common_8 + (2526300/17)*nu*potential*radial_speed, 160020*common_0*potential - 81900*common_1 - 1165500/17*common_2 + (916020/17)*common_3*speed_squared - 86940*common_4 + (399420/17)*nu*speed_squared + 12600*potential - 41580/17*speed_squared, -277200/17*common_34 + (166320/17)*common_5 + (110880/17)*common_6 - 27720/17*radial_speed, common_51),
        (2044*common_35 - 8176*common_36 + 315*common_39 + 1428*common_40 - 26460*common_41 - 10290*common_43 + 17640*common_44 + 19845*common_45 + 945*common_46 + 4284*common_47 - 30870*common_49 + 6615*common_70 - 1260*common_71 - 5712*common_73 - 4410*common_74 + 6132*common_75 + 41160*common_76 - 13230*common_80, -29400*common_29 - 141960*common_30 + 37800*common_31 - 88200*common_50 + 35490*common_52 - 50400*common_53 + 117600*common_54 + 106470*common_55 + 12600*potential*radial_speed*speed_squared, 175770*common_0*common_16*potential - 11718*common_12 + 170016*common_13*nu - 42504*common_13 + 58590*common_16*potential - 234360*common_18 - 35154*common_19 - 127512*common_21 + 46872*nu*potential*speed_squared, -207900*common_26 - 69300*common_8 + 277200*nu*potential*radial_speed, -218400*common_1 + 163800*common_23 + 54600*potential, 0, common_82, 0),
        ((8422680/17)*common_100 + (69426/17)*common_101 - 9821448/17*common_102 - 363825/34*common_103 + (105840/17)*common_104 - 1216446/17*common_105 + (1309770/17)*common_106 - 100170/17*common_107 + (451164/17)*common_108 - 86940/17*common_109 + (677754/17)*common_110 + (526995/17)*common_111 - 3292191/17*common_112 + (2388015/17)*common_113 - 3294270/17*common_114 - 8285340/17*common_115 - 86436/17*common_116 + (9878106/17)*common_117 + (890820/17)*common_118 - 79380/17*common_119 - 2457126/17*common_120 - 3018645/34*common_121 + (3329991/17)*common_122 + (1600830/17)*common_123 + (312417/17)*common_124 - 413824/17*common_83 + (1726172/17)*common_84 + (114345/17)*common_86 + (14895/34)*common_88 - 164063/17*common_89 - 509355/17*common_90 - 43155/17*common_91 + (855722/17)*common_92 + (104895/34)*common_93 - 1629040/17*common_94 + (420756/17)*common_95 - 1888110/17*common_96 - 6777/17*common_97 + (2272011/17)*common_98 - 443205/17*common_99, (38387580/17)*common_0*common_13*common_28 + (3482325/17)*common_0*common_38*potential*radial_speed - 22050/17*common_125 - 2096955/17*common_126 + (8383725/17)*common_13*common_28 + (8015175/17)*common_13*nu*radial_speed*speed_squared + (11355750/17)*common_28*common_3*potential*speed_squared + (191100/17)*common_28*potential*speed_squared + (9277212/17)*common_3*common_35*radial_speed + (29833020/17)*common_35*nu*radial_speed + (152775/17)*common_38*potential*radial_speed + (4451895/17)*common_56*nu*potential - 985635/17*common_57 - 6842850/17*common_58 - 1570590/17*common_59 - 1050525/17*common_60 - 37518180/17*common_61 - 31321206/17*common_62 - 4266675/17*common_63 - 6960450/17*common_64 - 8901837/17*common_65 - 3850875/17*common_66 - 283941/17*common_68 - 6556410/17*common_69, (6206634/17)*common_35 - 27957804/17*common_36 - 109053/17*common_39 + (652806/17)*common_40 - 12740490/17*common_41 - 12735660/17*common_42 - 15525972/17*common_43 + (1353240/17)*common_44 + (9126495/17)*common_45 - 2838402/17*common_46 + (5885460/17)*common_47 - 17947440/17*common_48 - 74594142/17*common_49 + (2835630/17)*common_70 + (821394/17)*common_71 - 4222428/17*common_73 - 670005/17*common_74 + (32815902/17)*common_75 + (69883002/17)*common_76 + (5750640/17)*common_77 + (3096576/17)*common_78 + (201936/17)*common_79 + (10280655/17)*common_80 + (17136882/17)*common_81, (77388150/17)*common_0*common_13*radial_speed + (14784210/17)*common_13*radial_speed + (18915750/17)*common_28*nu*potential - 4268250/17*common_29 + (17733870/17)*common_3*potential*radial_speed*speed_squared - 67660110/17*common_30 - 11910150/17*common_31 - 4079250/17*common_32 - 23454690/17*common_33 - 15372000/17*common_50 - 200970/17*common_53 + (506520/17)*potential*radial_speed*speed_squared, (330540/17)*common_12 - 5579280/17*common_13 - 3268440/17*common_14 + (29191680/17)*common_15 + (3994200/17)*common_17 - 17325000/17*common_18 + (12230400/17)*common_19 - 12746160/17*common_20 - 45225600/17*common_21 + (22602720/17)*common_22 + (14559300/17)*common_24 + (1587600/17)*common_25, -1000440/17*common_10 - 8358840/17*common_26 - 2673720/17*common_8 + (11030040/17)*nu*potential*radial_speed, (6070680/17)*common_0*potential - 3162600/17*common_1 - 2608200/17*common_2 + (1942920/17)*common_3*speed_squared - 3207960/17*common_4 + (965160/17)*nu*speed_squared + (491400/17)*potential - 110880/17*speed_squared, -378000/17*common_34 + (226800/17)*common_5 + (151200/17)*common_6 - 37800/17*radial_speed),
        (343980*common_0*common_13*common_28 + 198450*common_0*common_28*potential*speed_squared - 264600*common_125 - 218295*common_126 - 11025*common_127 - 33075*common_128 + 114660*common_13*common_28 + 148176*common_13*nu*radial_speed*speed_squared + 66150*common_28*potential*speed_squared + 144648*common_35*nu*radial_speed + 44100*common_38*nu*potential*radial_speed + 291060*common_56*nu*potential - 72765*common_57 - 36162*common_58 - 37044*common_59 - 458640*common_61 - 108486*common_62 - 111132*common_65, 63826*common_35 - 255304*common_36 + 14805*common_39 + 56322*common_40 - 1243620*common_41 - 410130*common_43 + 829080*common_44 + 932715*common_45 + 44415*common_46 + 168966*common_47 - 1230390*common_49 + 310905*common_70 - 59220*common_71 - 225288*common_73 - 207270*common_74 + 191478*common_75 + 1640520*common_76 - 621810*common_80, 251370*common_27 - 586530*common_29 - 2328480*common_30 + 754110*common_31 - 1759590*common_50 + 582120*common_52 - 1005480*common_53 + 2346120*common_54 + 1746360*common_55, 1918350*common_0*common_16*potential - 127890*common_12 + 1446480*common_13*nu - 361620*common_13 + 639450*common_16*potential - 2557800*common_18 - 383670*common_19 - 1084860*common_21 + 511560*nu*potential*speed_squared, -1323000*common_26 - 441000*common_8 + 1764000*nu*potential*radial_speed, -776160*common_1 + 582120*common_23 + 194040*potential, 0, common_82),
        ((6621615/17)*common_1*common_129 + (3902850/17)*common_10*common_87 + (218295/34)*common_12*common_56 + (24651270/17)*common_13*common_56 - 1486485/17*common_130 - 33391596/17*common_131 - 452025/34*common_132 - 29378034/17*common_133 + (1819125/17)*common_134 - 30824325/17*common_135 - 109978155/17*common_136 + (4729725/34)*common_137 - 17027010/17*common_138 - 7660716/17*common_139 + (992250/17)*common_140 - 15382710/17*common_141 - 6515775/17*common_142 - 7184205/17*common_143 - 158760/17*common_144*common_67 + (71625960/17)*common_144*common_72 - 3759525/17*common_145 - 27359640/17*common_146 - 11493531/17*common_147 - 130174380/17*common_148 - 11291805/17*common_149 - 71306550/17*common_150 - 27584550/17*common_151 + (37323783/17)*common_2*common_58 + (44095590/17)*common_20*common_56 + (108032400/17)*common_21*common_56 + (127997163/17)*common_28*common_36 + (31655421/17)*common_28*common_42 + (31496374/17)*common_34*common_83 + (7452270/17)*common_34*common_97 + (46271925/34)*common_38*common_50 + (1386144/17)*common_38*common_52 + (6497685/17)*common_58*speed_squared + (1610469/17)*common_6*common_97 + (7973282/17)*common_83*radial_speed + (1356075/17)*common_87*common_9, (474773355/17)*common_100 + (7724682/17)*common_101 - 416366559/17*common_102 + (43586235/34)*common_103 + (4223205/17)*common_104 - 3001089*common_105 + (56101815/17)*common_106 - 4552065/17*common_107 + 1105083*common_108 - 1104705/17*common_109 + (45121923/17)*common_110 + (10976490/17)*common_111 - 212937165/17*common_112 + (68875380/17)*common_113 - 146595015/17*common_114 - 481230855/17*common_115 - 8970318/17*common_116 + (435735090/17)*common_117 + (75060405/17)*common_118 - 1523970/17*common_119 - 120233169/17*common_120 - 111645765/34*common_121 + (222608925/17)*common_122 + (72229185/17)*common_123 - 12238947/17*common_124 - 11834410/17*common_83 + (49951076/17)*common_84 + (11642400/17)*common_86 + (409815/34)*common_88 - 448315*common_89 - 52536330/17*common_90 - 1352610/17*common_91 + 2234701*common_92 - 12561885/34*common_93 - 49081102/17*common_94 + (14088564/17)*common_95 - 106176105/17*common_96 - 1368486/17*common_97 + (95042556/17)*common_98 + (19355490/17)*common_99, (907453638/17)*common_0*common_13*common_28 + (55647270/17)*common_0*common_38*potential*radial_speed - 60566940/17*common_125 - 94720185/17*common_126 + (193248846/17)*common_13*common_28 + (23912280/17)*common_13*common_3*radial_speed*speed_squared + (220204656/17)*common_13*nu*radial_speed*speed_squared + (216769140/17)*common_28*common_3*potential*speed_squared + (17179155/17)*common_28*potential*speed_squared + (177526230/17)*common_3*common_35*radial_speed + (499543506/17)*common_35*nu*radial_speed + (121905/17)*common_38*potential*radial_speed + (141375780/17)*common_56*nu*potential - 31434480/17*common_57 - 112652988/17*common_58 - 45382302/17*common_59 - 7703640/17*common_60 - 867290004/17*common_61 - 549050964/17*common_62 - 77355810/17*common_63 - 77221305/17*common_64 - 244791666/17*common_65 - 74483010/17*common_66 - 183941100/17*common_69, (48704880/17)*common_35 - 224381430/17*common_36 - 38745/17*common_39 + (13684860/17)*common_40 - 196112700/17*common_41 - 115426290/17*common_42 - 189369810/17*common_43 + (64518300/17)*common_44 + (153236475/17)*common_45 - 26102160/17*common_46 + (92616300/17)*common_47 - 190644300/17*common_48 - 935056710/17*common_49 + (44055900/17)*common_70 + (3906000/17)*common_71 - 72962820/17*common_73 - 18687375/17*common_74 + (277732560/17)*common_75 + (856955610/17)*common_76 + (57307950/17)*common_77 + (33216750/17)*common_78 - 12006540/17*common_79 + (64838025/17)*common_80 + (236512710/17)*common_81, (552925800/17)*common_0*common_13*radial_speed + (101024280/17)*common_13*radial_speed + (162288000/17)*common_28*nu*potential - 37044000/17*common_29 + (112349160/17)*common_3*potential*radial_speed*speed_squared - 468200880/17*common_30 - 55213200/17*common_31 - 24255000/17*common_32 - 179187120/17*common_33 - 134284500/17*common_50 - 25648560/17*common_53 + (9022860/17)*potential*radial_speed*speed_squared, -502740/17*common_12 - 19068840/17*common_13 - 5080320/17*common_14 + (102717720/17)*common_15 + (20003760/17)*common_17 - 85836240/17*common_18 + (39354840/17)*common_19 - 46269720/17*common_20 - 165798360/17*common_21 + (84971880/17)*common_22 + (71838900/17)*common_24 + (5450760/17)*common_25, -652680/17*common_10 - 25207560/17*common_26 - 7408800/17*common_8 + (30958200/17)*nu*potential*radial_speed, (6448680/17)*common_0*potential - 3389400/17*common_1 - 2986200/17*common_2 + (2094120/17)*common_3*speed_squared - 3359160/17*common_4 + (1191960/17)*nu*speed_squared + (529200/17)*potential - 148680/17*speed_squared),
        (2302020*common_0*common_13*common_16*speed_squared + 1614060*common_0*common_16*common_35 + 893025*common_0*common_16*common_38*potential + 2837835*common_0*common_85*potential - 2152080*common_102 - 33075*common_104 - 166320*common_105 - 1190700*common_111 - 3069360*common_112 - 3274425*common_113 - 4365900*common_115 - 122472*common_116 + 767340*common_13*common_16*speed_squared + 5821200*common_13*common_37*nu + 163296*common_13*common_38*nu + 538020*common_16*common_35 + 297675*common_16*common_38*potential + 221760*common_35*nu*speed_squared + 4365900*common_37*nu*potential*speed_squared + 110656*common_83*nu - 27664*common_83 + 945945*common_85*potential + 44100*common_87*nu*potential - 11025*common_88 - 55440*common_89 - 3783780*common_90 - 1091475*common_93 - 82992*common_94 - 1455300*common_96 - 40824*common_97, 15743700*common_0*common_13*common_28 + 10716300*common_0*common_28*potential*speed_squared - 14288400*common_125 - 11787930*common_126 - 595350*common_127 - 1786050*common_128 + 5247900*common_13*common_28 + 6710256*common_13*nu*radial_speed*speed_squared + 3572100*common_28*potential*speed_squared + 5136768*common_35*nu*radial_speed + 2381400*common_38*nu*potential*radial_speed + 15717240*common_56*nu*potential - 3929310*common_57 - 1284192*common_58 - 1677564*common_59 - 20991600*common_61 - 3852576*common_62 - 5032692*common_65, -4762800*common_17*speed_squared + 885696*common_35 - 3542784*common_36 + 340200*common_39 + 1052352*common_40 - 28576800*common_41 - 7779240*common_43 + 19051200*common_44 + 21432600*common_45 + 1020600*common_46 + 3157056*common_47 - 23337720*common_49 + 7144200*common_70 - 1360800*common_71 - 4209408*common_73 + 2657088*common_75 + 31116960*common_76 - 14288400*common_80, 3175200*common_27 - 7408800*common_29 - 23284800*common_30 + 9525600*common_31 - 22226400*common_50 + 5821200*common_52 - 12700800*common_53 + 29635200*common_54 + 17463600*common_55, 14288400*common_0*common_16*potential - 952560*common_12 + 7902720*common_13*nu - 1975680*common_13 + 4762800*common_16*potential - 19051200*common_18 - 2857680*common_19 - 5927040*common_21 + 3810240*nu*potential*speed_squared, -5715360*common_26 - 1905120*common_8 + 7620480*nu*potential*radial_speed, -1693440*common_1 + 1270080*common_23 + 423360*potential, 0),
        ((1501217550/17)*common_0*common_13*common_37*speed_squared + (7566930/17)*common_0*common_13*common_87 + (19465880/17)*common_0*common_152 + (206487225/34)*common_0*common_16*common_87*potential + (1912604211/17)*common_0*common_35*common_37 + (42048405/17)*common_0*common_35*common_38 + (50734600/17)*common_0*common_83*speed_squared + (744458715/34)*common_0*common_85*potential*speed_squared - 99324225/17*common_1*common_153 - 23648625/34*common_12*common_85 + (303219315/17)*common_13*common_16*common_38*nu + (176164065/17)*common_13*common_3*common_85 + (1681344/17)*common_13*common_3*common_87 + (330756615/17)*common_13*common_37*speed_squared + (1655595585/17)*common_13*common_85*nu - 371313180/17*common_13*common_85 + (1342989/17)*common_13*common_87 - 7123410/17*common_15*common_87 - 3606120/17*common_152*common_3 - 23143512/17*common_152*nu + (5742280/17)*common_152 - 70945875/34*common_153*common_23 + (255405150/17)*common_153*common_3*potential + (22297275/17)*common_153*potential - 3759525/17*common_154*common_23 + (3902850/17)*common_154*common_3*potential + (1356075/17)*common_154*nu*potential - 452025/34*common_154*potential + (256537260/17)*common_16*common_3*common_35*speed_squared + (131683104/17)*common_16*common_3*common_83 + (794564064/17)*common_16*common_35*nu*speed_squared + (539870520/17)*common_16*common_83*nu - 127578408/17*common_16*common_83 + (10021725/34)*common_16*common_87*potential - 172898496/17*common_16*common_89 - 522468996/17*common_16*common_94 - 31752000/17*common_18*common_87 - 904389318/17*common_2*common_98 - 692431740/17*common_20*common_85 - 1627878735/17*common_21*common_85 - 117879300/17*common_25*common_87 + (523908000/17)*common_3*common_37*common_38*potential + (424388097/17)*common_35*common_37 + (6409341/17)*common_35*common_38 - 1859824008/17*common_36*common_37 - 31936716/17*common_36*common_38 + (80769150/17)*common_37*common_38*nu*potential - 20738025/34*common_37*common_39 - 467521740/17*common_37*common_42 - 782587575/34*common_37*common_46 - 14144004/17*common_38*common_42 - 63374400/17*common_38*common_43 - 306695025/17*common_38*common_49 - 15464925/17*common_38*common_81 - 17219832/17*common_67*common_83 - 68028660/17*common_67*common_96 - 44675064/17*common_72*common_83 - 1512999810/17*common_72*common_96 + (9937568/17)*common_83*speed_squared + (22702680/17)*common_85*nu*potential*speed_squared, (591882039/17)*common_0*common_13*common_38*radial_speed + (7015552425/17)*common_0*common_13*common_56 + (939859200/17)*common_0*common_28*common_38*potential + (1809486504/17)*common_0*common_35*radial_speed*speed_squared + (1092705684/17)*common_0*common_83*radial_speed + (729323595/17)*common_129*nu*potential + (455059080/17)*common_13*common_28*common_3*speed_squared + (4844687400/17)*common_13*common_28*nu*speed_squared + (701379/17)*common_13*common_3*common_38*radial_speed + (114158457/17)*common_13*common_38*radial_speed + (1557956295/17)*common_13*common_56 - 161756595/17*common_130 - 1103792928/17*common_131 - 5225850/17*common_132 - 1410208632/17*common_133 - 515394495/17*common_134 - 1514873304/17*common_135 - 6957828045/17*common_136 - 266756490/17*common_137 - 848512665/17*common_138 - 306604200/17*common_139 - 8235675/17*common_140 - 1053042480/17*common_141 - 66976875/17*common_142 - 556536771/17*common_143 - 167524875/17*common_145 - 818387955/17*common_146 - 579310704/17*common_147 - 6527807244/17*common_148 - 1053352755/17*common_149 - 4967950680/17*common_150 - 1429534575/17*common_151 + (1761616836/17)*common_28*common_3*common_35 + (6207383700/17)*common_28*common_35*nu + (2257825185/17)*common_3*common_56*potential*speed_squared + (203642775/17)*common_3*common_87*potential*radial_speed + (323907444/17)*common_35*radial_speed*speed_squared + (139708800/17)*common_56*potential*speed_squared + (260011248/17)*common_83*radial_speed + (40913775/17)*common_87*nu*potential*radial_speed, (12438392400/17)*common_100 + (260289288/17)*common_101 - 8108140068/17*common_102 + 68981220*common_103 + (80986500/17)*common_104 - 974134896/17*common_105 + 69854400*common_106 - 101795400/17*common_107 + (365776656/17)*common_108 + (44055900/17)*common_109 + (1291793832/17)*common_110 - 61916400/17*common_111 - 5998894020/17*common_112 + (634643100/17)*common_113 - 3217271400/17*common_114 - 12858161400/17*common_115 - 296196264/17*common_116 + (8818285644/17)*common_117 + (2420825400/17)*common_118 + (9843120/17)*common_119 - 2673789804/17*common_120 - 956132100/17*common_121 + (6372941652/17)*common_122 + (1605063600/17)*common_123 - 836526852/17*common_124 - 138654432/17*common_83 + (592959360/17)*common_84 + 24449040*common_86 + (1351350/17)*common_88 - 154503216/17*common_89 - 110020680*common_90 - 15819300/17*common_91 + (749107632/17)*common_92 - 349470450/17*common_93 - 612819648/17*common_94 + (202004544/17)*common_95 - 2778318900/17*common_96 - 51349788/17*common_97 + (1826994204/17)*common_98 + (1368114300/17)*common_99, (12974872680/17)*common_0*common_13*common_28 + (538423200/17)*common_0*common_38*potential*radial_speed - 1538913600/17*common_125 - 1852729200/17*common_126 - 36571500/17*common_127 + (2707942860/17)*common_13*common_28 + (642229560/17)*common_13*common_3*radial_speed*speed_squared + (3393721800/17)*common_13*nu*radial_speed*speed_squared + (2647587600/17)*common_28*common_3*potential*speed_squared + (392931000/17)*common_28*potential*speed_squared + (1998437280/17)*common_3*common_35*radial_speed + (4881502080/17)*common_35*nu*radial_speed + (70988400/17)*common_38*nu*potential*radial_speed + (2471893200/17)*common_56*nu*potential - 553675500/17*common_57 - 1081076640/17*common_58 - 717209640/17*common_59 - 12183083640/17*common_61 - 5635437360/17*common_62 - 902550600/17*common_63 - 227026800/17*common_64 - 3784951800/17*common_65 - 919787400/17*common_66 - 2945527200/17*common_69, (225125040/17)*common_35 - 1068863040/17*common_36 + 521640*common_39 + 8250480*common_40 - 1920466800/17*common_41 - 653399040/17*common_42 - 1520497440/17*common_43 + (873180000/17)*common_44 + (1552937400/17)*common_45 - 9389520*common_46 + 50752800*common_47 - 1397088000/17*common_48 - 7707374640/17*common_49 + (435531600/17)*common_70 - 589680*common_71 - 41176800*common_73 - 230995800/17*common_74 + (1422981840/17)*common_75 + (6918866640/17)*common_76 + (399016800/17)*common_77 + 14424480*common_78 - 11128320*common_79 + (132564600/17)*common_80 + (2107380240/17)*common_81, (2580696720/17)*common_0*common_13*radial_speed + (450772560/17)*common_13*radial_speed + (916045200/17)*common_28*nu*potential - 211150800/17*common_29 + (486758160/17)*common_3*potential*radial_speed*speed_squared - 2116270800/17*common_30 - 140026320/17*common_31 - 100018800/17*common_32 - 890855280/17*common_33 - 762048000/17*common_50 - 233377200/17*common_53 + (68584320/17)*potential*radial_speed*speed_squared, -6773760/17*common_12 - 36832320/17*common_13 + (11854080/17)*common_14 + (206599680/17)*common_15 + (62233920/17)*common_17 - 265446720/17*common_18 + (69007680/17)*common_19 - 102453120/17*common_20 - 353928960/17*common_21 + (190512000/17)*common_22 + (222899040/17)*common_24 + (10160640/17)*common_25, (2116800/17)*common_10 - 37679040/17*common_26 - 10039680/17*common_8 + (42759360/17)*common_9),
    )

def j3(potential, radial_speed, speed_squared, nu, log_one, log_infrared):
    common_0 = nu*potential
    common_1 = (16/9)*common_0
    common_2 = nu*speed_squared
    common_3 = nu**2
    common_4 = common_3*potential
    common_5 = (10/9)*radial_speed
    common_6 = potential**2
    common_7 = speed_squared**2
    common_8 = common_2*potential
    common_9 = common_6*nu
    common_10 = common_7*nu
    common_11 = radial_speed**2
    common_12 = common_11*potential
    common_13 = common_0*common_11
    common_14 = nu**3
    common_15 = common_14*common_6
    common_16 = common_14*common_7
    common_17 = common_11*common_4
    common_18 = potential*radial_speed
    common_19 = radial_speed*speed_squared
    common_20 = common_0*radial_speed
    common_21 = common_2*radial_speed
    common_22 = common_4*radial_speed
    common_23 = common_3*speed_squared
    common_24 = common_23*radial_speed
    common_25 = common_14*common_18
    common_26 = common_14*common_19
    common_27 = (98/165)*common_11
    common_28 = common_14*potential
    common_29 = common_14*speed_squared
    common_30 = common_11*common_3
    common_31 = common_9*radial_speed
    common_32 = common_3*common_6
    common_33 = common_32*radial_speed
    common_34 = (32/5)*common_3*potential*speed_squared - 592/15*common_32 + (544/45)*common_6*nu - 8/5*common_8
    common_35 = 6*nu
    common_36 = 2 - common_35
    common_37 = (8/3)*radial_speed
    common_38 = potential*speed_squared
    common_39 = common_38*radial_speed
    common_40 = radial_speed**3
    common_41 = common_40*potential
    common_42 = common_6*radial_speed
    common_43 = common_18*common_2
    common_44 = common_0*common_40
    common_45 = common_19*common_4
    common_46 = common_14*common_38
    common_47 = common_46*radial_speed
    common_48 = common_4*common_40
    common_49 = common_28*common_40
    common_50 = common_15*radial_speed
    common_51 = common_4*speed_squared
    common_52 = potential**3
    common_53 = common_52*nu
    common_54 = common_2*common_6
    common_55 = common_3*common_52
    common_56 = common_11*common_9
    common_57 = common_23*common_6
    common_58 = common_11*common_32
    common_59 = radial_speed**4
    common_60 = common_59*potential
    common_61 = common_7*potential
    common_62 = common_6*speed_squared
    common_63 = common_12*speed_squared
    common_64 = common_11*common_6
    common_65 = common_4*common_7
    common_66 = common_15*speed_squared
    common_67 = common_17*speed_squared
    common_68 = common_12*common_29
    common_69 = common_40*common_9
    common_70 = common_55*radial_speed
    common_71 = common_23*common_42
    common_72 = 24*common_0 - 8*potential
    common_73 = radial_speed**5
    common_74 = common_53*radial_speed
    common_75 = common_40*common_6
    common_76 = common_2*common_41
    common_77 = common_4*common_73
    common_78 = common_14*common_52
    common_79 = common_78*radial_speed
    common_80 = common_16*common_18
    common_81 = common_15*common_40
    common_82 = common_0*common_59
    common_83 = common_28*common_59
    common_84 = potential**4
    common_85 = common_2*common_52
    common_86 = common_3*common_84
    common_87 = common_2*common_64
    common_88 = common_32*common_59
    common_89 = common_32*common_7
    common_90 = common_30*common_52
    common_91 = common_12*common_2
    common_92 = common_84*nu
    common_93 = radial_speed**6
    common_94 = speed_squared**3
    common_95 = common_52*radial_speed
    common_96 = common_2*common_42
    common_97 = common_32*common_40
    return (
        (1 - 3*nu, 0, 0),
        (-common_1 - 77/18*common_2 + (185/18)*common_3*speed_squared - 86/9*common_4 + (14/9)*potential + (41/90)*speed_squared, common_3*common_5 - common_5*nu + (2/9)*radial_speed, (7/9)*common_3 - 7/9*nu + 7/45),
        (-4159/792*common_10 + (32221/1980)*common_11*common_14*potential - 23/396*common_12 - 637/990*common_13 + (21127/396)*common_14*potential*speed_squared - 8099/396*common_15 - 171539/3960*common_16 - 1861/990*common_17 + (9925/198)*common_3*common_6 + (52409/1980)*common_3*common_7 + (6307/198)*common_3*potential*speed_squared - 45/44*common_6 + (1349/3960)*common_7 - 19381/990*common_8 - 988/99*common_9 + (1597/660)*potential*speed_squared, (133/66)*common_18 + (73/495)*common_19 - 81/55*common_20 - 841/495*common_21 - 3914/165*common_22 + (3002/495)*common_24 + (3089/330)*common_25 - 3151/495*common_26, -1792/495*common_0 + (14/165)*common_11 - common_14*common_27 - 259/198*common_2 + (2219/495)*common_23 - common_27*nu + (427/99)*common_28 - 4529/990*common_29 + (196/165)*common_30 - 227/99*common_4 + (751/495)*potential + (119/990)*speed_squared),
        (-112/9*common_31 + (112/3)*common_33, common_34, 0),
        (0, common_36, 0),
        (-23/9*common_20 - 55/9*common_22 + (59/45)*potential*radial_speed, -8/9*common_0 - 29/3*common_2 + (65/3)*common_3*speed_squared - 196/9*common_4 + (116/45)*potential + (17/15)*speed_squared, common_3*common_37 - common_37*nu + (8/15)*radial_speed),
        ((743/990)*common_31 - 15283/330*common_33 + (569/1980)*common_39 + (23/132)*common_41 + (2011/495)*common_42 + (229/330)*common_43 - 287/330*common_44 - 1934/165*common_45 + (51809/1980)*common_47 + (541/330)*common_48 - 4831/660*common_49 - 5488/165*common_50, -8053/660*common_10 + (17047/495)*common_12*common_14 + (7/99)*common_12 - 4502/495*common_13 - 8837/165*common_15 - 61381/660*common_16 + (10621/495)*common_17 + (19471/330)*common_3*common_7 + (66277/495)*common_32 + (89/15)*common_38 + (5561/45)*common_46 + (1363/45)*common_51 - 2456/495*common_6 + (547/660)*common_7 - 1586/45*common_8 - 10319/495*common_9, (449/90)*common_18 + (92/165)*common_19 - 104/9*common_20 - 908/165*common_21 - 1051/45*common_22 + (2872/165)*common_24 + (1741/90)*common_25 - 2756/165*common_26),
        (-40/9*common_53 - 112/9*common_54 + (248/15)*common_55 + (112/3)*common_56 + (536/15)*common_57 - 112*common_58, -304/9*common_31 + (304/3)*common_33, common_34),
        (common_35*potential - 2*potential, 0, common_36),
        (-59/15*common_12 + (32/3)*common_13 + (28/3)*common_17 + (143/9)*common_32 - 82/45*common_38 - 88/9*common_51 + (37/9)*common_6 + (64/9)*common_8 - 149/9*common_9, -320/9*common_20 - 136/9*common_22 + (572/45)*potential*radial_speed, common_1 - 37/3*common_2 + (73/3)*common_23 - 220/9*common_4 + (92/45)*potential + (5/3)*speed_squared),
        ((76907/1980)*common_11*common_14*common_6 + (12602/495)*common_11*common_6*nu + (661/110)*common_11*nu*potential*speed_squared + (3073/99)*common_14*common_52 + (94/33)*common_14*common_59*potential + (14933/495)*common_14*common_7*potential + (472/33)*common_3*common_59*potential + (2435/99)*common_52*nu - 17/11*common_52 - 7477/99*common_55 - 767/330*common_57 - 3721/495*common_58 + (79/132)*common_59*nu*potential + (607/165)*common_6*nu*speed_squared - 115/132*common_60 - 829/495*common_61 - 333/220*common_62 - 56/165*common_63 - 5801/1980*common_64 - 4113/110*common_65 - 9681/220*common_66 - 738/55*common_67 - 5273/165*common_68 + (8687/660)*common_7*nu*potential, (59113/990)*common_14*potential*radial_speed*speed_squared + (4713/55)*common_3*potential*radial_speed*speed_squared - 6157/495*common_31 - 105199/495*common_33 + (138/11)*common_40*nu*potential + (3/22)*common_40*potential - 6734/165*common_43 - 617/11*common_48 - 197/22*common_49 - 40433/495*common_50 + (1211/55)*common_6*radial_speed + (5101/990)*potential*radial_speed*speed_squared, -779/44*common_10 + (187/6)*common_11*common_14*potential + 70*common_11*common_3*potential + (5/2)*common_11*potential - 98/3*common_13 + (150109/990)*common_14*potential*speed_squared - 66893/990*common_15 - 4827/44*common_16 + (8942/55)*common_3*common_6 + (1681/22)*common_3*common_7 - 848/165*common_51 - 7739/990*common_6 + (61/44)*common_7 - 6638/165*common_8 - 1039/55*common_9 + (9733/990)*potential*speed_squared),
        (560*common_3*common_40*common_6 + (896/45)*common_52*nu*radial_speed + (608/5)*common_6*nu*radial_speed*speed_squared - 560/3*common_69 - 208/3*common_70 - 360*common_71, -832/45*common_53 - 928/15*common_54 + (1024/15)*common_55 + 176*common_56 + (896/5)*common_57 - 528*common_58, -128/3*common_31 + 128*common_33),
        (-18*common_20 + 6*potential*radial_speed, common_72, 0),
        (48*common_3*potential*radial_speed*speed_squared - 536/9*common_33 - 12/5*common_39 + (59/3)*common_40*potential - 1154/45*common_42 - 160/3*common_44 - 140/3*common_48 + (848/9)*common_6*nu*radial_speed, -46*common_12 + 134*common_13 + 46*common_17 + (646/9)*common_32 + (26/15)*common_38 - 142/3*common_51 + (334/45)*common_6 + (10/3)*common_8 - 370/9*common_9, -290/3*common_20 - 10/3*common_22 + (94/3)*potential*radial_speed),
        (-553/132*common_0*common_73 + (1883/11)*common_14*common_40*potential*speed_squared + (42224/165)*common_14*common_6*radial_speed*speed_squared - 2044/99*common_19*common_6 - 6043/220*common_20*common_7 - 658/33*common_28*common_73 + (1076/11)*common_3*common_40*common_6 + (4102/33)*common_3*common_40*potential*speed_squared + (261112/495)*common_3*common_52*radial_speed + (9387/110)*common_3*common_7*potential*radial_speed - 59/33*common_38*common_40 + (9821/495)*common_52*radial_speed + (44557/495)*common_6*nu*radial_speed*speed_squared - 952/55*common_69 + (239/55)*common_7*potential*radial_speed - 8924/495*common_71 + (805/132)*common_73*potential - 34957/165*common_74 - 4223/165*common_75 - 1825/66*common_76 - 3304/33*common_77 - 74731/495*common_79 - 8493/55*common_80 - 8896/33*common_81, (21799/99)*common_11*common_14*common_6 + (38011/99)*common_11*common_3*common_6 + (25189/495)*common_11*common_6 + (2146/11)*common_11*nu*potential*speed_squared + (74696/495)*common_14*common_52 + (1600/11)*common_14*common_7*potential + (11684/33)*common_3*common_59*potential + (4292/165)*common_52*nu + (10754/495)*common_52 - 17748/55*common_55 - 94912/495*common_56 - 36259/495*common_57 + (15737/495)*common_6*nu*speed_squared - 80/33*common_60 - 142/33*common_61 - 511/55*common_62 - 173/11*common_63 - 2131/33*common_65 - 108563/495*common_66 - 6239/11*common_67 - 1374/11*common_68 + (3892/165)*common_7*nu*potential - 2278/33*common_82 - 562/33*common_83, (596/495)*common_31 - 238132/495*common_33 + (4088/165)*common_39 - 81/11*common_41 + (19234/495)*common_42 - 10932/55*common_43 + (974/11)*common_44 + (24064/55)*common_45 - 2656/165*common_47 - 2641/11*common_48 + (270/11)*common_49 + (1916/495)*common_50),
        (3480*common_11*common_3*common_6*speed_squared + (64/3)*common_11*common_52*nu + (1528/15)*common_3*common_52*speed_squared + (3920/3)*common_59*common_6*nu + (608/5)*common_6*common_7*nu + (368/45)*common_84*nu - 1648/45*common_85 - 416/15*common_86 - 1168*common_87 - 3920*common_88 - 360*common_89 - 16*common_90, 3760*common_3*common_40*common_6 + (64/9)*common_52*nu*radial_speed + 800*common_6*nu*radial_speed*speed_squared - 3760/3*common_69 - 1184/15*common_70 - 11856/5*common_71, -208/45*common_53 - 1712/15*common_54 + (80/3)*common_55 + 304*common_56 + 336*common_57 - 912*common_58),
        (90*common_11*nu*potential - 30*common_12 - common_35*common_6 + 2*common_6 - 18*common_8 + 6*potential*speed_squared, -108*common_20 + 36*potential*radial_speed, common_72),
        ((1309/3)*common_11*common_3*common_6 + (581/3)*common_11*common_6 + 71*common_11*potential*speed_squared + (980/3)*common_3*common_59*potential + 48*common_3*common_7*potential + (92/9)*common_52 - 298/9*common_53 - 2/9*common_55 - 2131/3*common_56 - 704/9*common_57 + (1120/3)*common_59*nu*potential + (818/9)*common_6*nu*speed_squared - 413/3*common_60 - 12/5*common_61 - 1034/45*common_62 - 380*common_67 - 160*common_91, 330*common_3*potential*radial_speed*speed_squared - 4436/9*common_33 - 102*common_39 + (808/3)*common_40*potential - 1072/9*common_42 - 2330/3*common_44 - 970/3*common_48 + (4460/9)*common_6*nu*radial_speed + 258*nu*potential*radial_speed*speed_squared, -140*common_12 + 424*common_13 + 56*common_17 + (676/9)*common_32 + (496/15)*common_38 - 152/3*common_51 - 1076/45*common_6 - 280/3*common_8 + (500/9)*common_9),
        ((1659/44)*common_0*common_93 - 6043/220*common_0*common_94 + (17741/180)*common_10*common_6 - 26933/99*common_11*common_52 + (1027663/495)*common_11*common_53 - 143666/55*common_11*common_66 + (767293/495)*common_11*common_78 + (14142/11)*common_12*common_16 - 298/11*common_12*common_7 + (2393/44)*common_13*common_7 + (809/33)*common_14*common_84 + (323275/132)*common_15*common_59 + (4303/15)*common_16*common_6 - 1183/22*common_17*common_7 + (7595/44)*common_2*common_60 + (4496/9)*common_23*common_52 + (147931/165)*common_23*common_64 + (1974/11)*common_28*common_93 - 8493/55*common_28*common_94 - 2114/9*common_29*common_52 + (5677/132)*common_38*common_59 + (9912/11)*common_4*common_93 + (9387/110)*common_4*common_94 - 42833/33*common_46*common_59 - 15078/11*common_51*common_59 + (196/3)*common_52*speed_squared + (326651/660)*common_59*common_6 - 16827/22*common_59*common_9 - 901/45*common_6*common_7 + (694/165)*common_64*speed_squared - 1061/11*common_84 - 2998/9*common_85 - 11893/66*common_86 - 162913/330*common_87 - 75247/66*common_88 - 6791/90*common_89 - 2136176/495*common_90 + (23189/66)*common_92 - 2415/44*common_93*potential + (239/55)*common_94*potential, (29660/33)*common_14*common_40*potential*speed_squared + (111778/55)*common_14*common_6*radial_speed*speed_squared + (238/3)*common_14*common_73*potential - 542/55*common_18*common_7 - 42348/55*common_22*common_7 + (148525/33)*common_3*common_40*potential*speed_squared + (1461836/495)*common_3*common_52*radial_speed + (421636/495)*common_3*common_6*radial_speed*speed_squared + (111206/33)*common_40*common_6*nu + (719/11)*common_40*potential*speed_squared + (74428/495)*common_6*radial_speed*speed_squared + (29093/110)*common_7*nu*potential*radial_speed + (2849/6)*common_73*nu*potential + (175/6)*common_73*potential - 44188/495*common_74 - 171044/165*common_75 - 43127/33*common_76 - 8036/3*common_77 - 774038/495*common_79 - 54726/55*common_80 - 23292/11*common_81 - 7984/33*common_95 - 248438/495*common_96 - 75628/33*common_97, -28904/165*common_0*common_7 + (2774/15)*common_11*common_15 + (21344/165)*common_16*potential + (51299/33)*common_4*common_59 + (10712/99)*common_52 - 29620/99*common_53 + (116041/495)*common_54 - 4780/99*common_55 - 1815*common_56 - 345767/495*common_57 + (40927/15)*common_58 + (1135/33)*common_60 + (1126/55)*common_61 - 13139/495*common_62 - 6168/55*common_63 + (1108/3)*common_64 + (61537/165)*common_65 - 31243/165*common_66 - 143002/55*common_67 - 164/55*common_68 + (1544/11)*common_78 - 16888/33*common_82 - 4612/33*common_83 + (58136/55)*common_91),
        (-2944*common_10*common_42 + (40208/3)*common_2*common_75 + (6416/5)*common_2*common_95 - 40040*common_23*common_75 - 18528/5*common_23*common_95 + 35280*common_32*common_73 + 8760*common_33*common_7 - 1808*common_40*common_53 + 5088*common_40*common_55 - 11760*common_73*common_9 + (4928/15)*common_86*radial_speed - 4928/45*common_92*radial_speed, 30096*common_11*common_3*common_6*speed_squared + 1552*common_11*common_52*nu + (8736/5)*common_3*common_52*speed_squared + (34160/3)*common_59*common_6*nu + (5216/5)*common_6*common_7*nu + (4288/45)*common_84*nu - 9184/15*common_85 - 3904/15*common_86 - 10096*common_87 - 34160*common_88 - 15456/5*common_89 - 4272*common_90, -8320/3*common_69 + (4832/5)*common_70 - 26016/5*common_71 - 5312/15*common_74 + (8752/5)*common_96 + 8320*common_97),
    )

def j3_tail(potential, radial_speed, speed_squared, nu, log_one, log_infrared):
    common_0 = nu*potential
    common_1 = (16/9)*common_0
    common_2 = nu*speed_squared
    common_3 = nu**2
    common_4 = common_3*potential
    common_5 = (10/9)*radial_speed
    common_6 = 6*nu
    common_7 = 2 - common_6
    common_8 = common_0*radial_speed
    common_9 = common_4*radial_speed
    common_10 = (8/3)*radial_speed
    common_11 = potential*speed_squared
    common_12 = potential**2
    common_13 = common_2*potential
    common_14 = common_12*nu
    common_15 = radial_speed**2
    common_16 = common_15*potential
    common_17 = common_0*common_15
    common_18 = common_4*speed_squared
    common_19 = common_12*common_3
    common_20 = common_15*common_4
    common_21 = common_3*speed_squared
    common_22 = 24*common_0 - 8*potential
    common_23 = potential*radial_speed*speed_squared
    common_24 = radial_speed**3
    common_25 = common_12*radial_speed
    common_26 = common_0*common_24
    common_27 = common_24*common_4
    common_28 = common_19*radial_speed
    common_29 = potential**3
    common_30 = common_29*nu
    common_31 = radial_speed**4
    common_32 = common_31*potential
    common_33 = speed_squared**2
    common_34 = common_33*potential
    common_35 = common_12*speed_squared
    common_36 = common_29*common_3
    common_37 = common_16*common_2
    common_38 = common_14*common_15
    common_39 = common_12*common_21
    common_40 = common_15*common_18
    common_41 = radial_speed**5
    common_42 = common_29*radial_speed
    common_43 = common_0*common_41
    common_44 = common_11*common_24
    common_45 = common_12*common_24
    common_46 = common_33*common_8
    common_47 = common_2*common_25
    common_48 = common_4*common_41
    common_49 = common_36*radial_speed
    common_50 = common_33*common_9
    common_51 = common_19*common_24
    common_52 = common_0*common_31
    common_53 = common_0*common_33
    common_54 = common_12*common_2
    common_55 = common_11*common_15
    common_56 = common_12*common_15
    common_57 = common_31*common_4
    common_58 = common_33*common_4
    common_59 = common_15*common_19
    common_60 = common_14*radial_speed
    common_61 = potential**4
    common_62 = radial_speed**6
    common_63 = speed_squared**3
    return (
        (1 - 3*nu, 0, 0),
        (-common_1 - 77/18*common_2 + (185/18)*common_3*speed_squared - 86/9*common_4 + (14/9)*potential + (41/90)*speed_squared, common_3*common_5 - common_5*nu + (2/9)*radial_speed, (7/9)*common_3 - 7/9*nu + 7/45),
        (0, common_7, 0),
        (-23/9*common_8 - 55/9*common_9 + (59/45)*potential*radial_speed, -8/9*common_0 - 29/3*common_2 + (65/3)*common_3*speed_squared - 196/9*common_4 + (116/45)*potential + (17/15)*speed_squared, common_10*common_3 - common_10*nu + (8/15)*radial_speed),
        (common_6*potential - 2*potential, 0, common_7),
        (-82/45*common_11 + (37/9)*common_12 + (64/9)*common_13 - 149/9*common_14 - 59/15*common_16 + (32/3)*common_17 - 88/9*common_18 + (143/9)*common_19 + (28/3)*common_20, -320/9*common_8 - 136/9*common_9 + (572/45)*potential*radial_speed, common_1 - 37/3*common_2 + (73/3)*common_21 - 220/9*common_4 + (92/45)*potential + (5/3)*speed_squared),
        (-18*common_8 + 6*potential*radial_speed, common_22, 0),
        ((848/9)*common_12*nu*radial_speed - 12/5*common_23 + (59/3)*common_24*potential - 1154/45*common_25 - 160/3*common_26 - 140/3*common_27 - 536/9*common_28 + 48*common_3*potential*radial_speed*speed_squared, (26/15)*common_11 + (334/45)*common_12 + (10/3)*common_13 - 370/9*common_14 - 46*common_16 + 134*common_17 - 142/3*common_18 + (646/9)*common_19 + 46*common_20, -290/3*common_8 - 10/3*common_9 + (94/3)*potential*radial_speed),
        (-common_12*common_6 + 2*common_12 - 18*common_13 + 90*common_15*nu*potential - 30*common_16 + 6*potential*speed_squared, -108*common_8 + 36*potential*radial_speed, common_22),
        ((1309/3)*common_12*common_15*common_3 + (581/3)*common_12*common_15 + (818/9)*common_12*nu*speed_squared + 71*common_15*potential*speed_squared + (92/9)*common_29 + (980/3)*common_3*common_31*potential + 48*common_3*common_33*potential - 298/9*common_30 + (1120/3)*common_31*nu*potential - 413/3*common_32 - 12/5*common_34 - 1034/45*common_35 - 2/9*common_36 - 160*common_37 - 2131/3*common_38 - 704/9*common_39 - 380*common_40, (4460/9)*common_12*nu*radial_speed - 102*common_23 + (808/3)*common_24*potential - 1072/9*common_25 - 2330/3*common_26 - 970/3*common_27 - 4436/9*common_28 + 330*common_3*potential*radial_speed*speed_squared + 258*nu*potential*radial_speed*speed_squared, (496/15)*common_11 - 1076/45*common_12 - 280/3*common_13 + (500/9)*common_14 - 140*common_16 + 424*common_17 - 152/3*common_18 + (676/9)*common_19 + 56*common_20),
        (-90*common_23 + 210*common_24*potential - 630*common_26 + 270*nu*potential*radial_speed*speed_squared, 48*common_11 - 16*common_12 - 144*common_13 + 48*common_14 - 240*common_16 + 720*common_17, -180*common_8 + 60*potential*radial_speed),
        (6528*common_12*common_24*nu + 1580*common_12*common_3*radial_speed*speed_squared + (2894/5)*common_12*radial_speed*speed_squared + (11900/3)*common_24*common_3*potential*speed_squared + (7840/3)*common_24*nu*potential*speed_squared + (9604/9)*common_29*nu*radial_speed + 154*common_33*potential*radial_speed + 1239*common_41*potential - 15022/45*common_42 - 3360*common_43 - 3143/3*common_44 - 1770*common_45 - 320*common_46 - 2156*common_47 - 2940*common_48 - 1240/9*common_49 - 1000*common_50 - 4200*common_51, (12032/45)*common_29 - 7856/9*common_30 - 6482/3*common_32 - 534/5*common_34 - 746/5*common_35 + (2000/9)*common_36 - 3940*common_37 - 4910*common_38 - 698*common_39 - 3380*common_40 + (18550/3)*common_52 + 258*common_53 + 606*common_54 + 1460*common_55 + 1210*common_56 + (8750/3)*common_57 + 426*common_58 + 4586*common_59, (472/5)*common_12*radial_speed - 2406/5*common_23 + (2908/3)*common_24*potential - 8690/3*common_26 - 1810/3*common_27 - 948*common_28 + 594*common_3*potential*radial_speed*speed_squared - 52*common_60 + 1386*nu*potential*radial_speed*speed_squared),
        (630*common_12*common_15*nu + 42*common_12*speed_squared + 1260*common_15*potential*speed_squared + 16*common_29 - 48*common_30 + 5670*common_31*nu*potential - 1890*common_32 + 270*common_33*nu*potential - 90*common_34 - 3780*common_37 - 126*common_54 - 210*common_56, 360*common_12*radial_speed - 900*common_23 + 2100*common_24*potential - 6300*common_26 - 1080*common_60 + 2700*nu*potential*radial_speed*speed_squared, 108*common_11 - 76*common_12 - 324*common_13 + 228*common_14 - 540*common_16 + 1620*common_17),
        (36960*common_0*common_62 - 320*common_0*common_63 + 15624*common_11*common_31 + 19121*common_12*common_31 + (2868/5)*common_12*common_33 - 71335*common_14*common_31 - 2094*common_14*common_33 + (33934/5)*common_15*common_29 - 22090*common_15*common_30 - 52617/5*common_15*common_35 + 4906*common_15*common_36 - 4221*common_16*common_33 + 10080*common_17*common_33 - 50400*common_18*common_31 + 49945*common_19*common_31 + 1776*common_19*common_33 + (27334/9)*common_2*common_29 - 40320*common_2*common_32 + 39189*common_2*common_56 + 18900*common_20*common_33 - 8206/9*common_21*common_29 - 29841*common_21*common_56 - 41194/45*common_29*speed_squared + (104/9)*common_3*common_61 + 32340*common_4*common_62 - 1000*common_4*common_63 - 308/9*common_61*nu + (22/9)*common_61 - 13629*common_62*potential + 154*common_63*potential, 54320*common_12*common_24*nu + 19164*common_12*common_3*radial_speed*speed_squared + (19872/5)*common_12*radial_speed*speed_squared + 43260*common_24*common_3*potential*speed_squared + 57540*common_24*nu*potential*speed_squared + (180400/9)*common_29*nu*radial_speed + 3762*common_33*potential*radial_speed + 21924*common_41*potential - 271516/45*common_42 - 62370*common_43 - 20958*common_44 - 13276*common_45 - 9810*common_46 - 16380*common_47 - 32130*common_48 - 56872/9*common_49 - 10890*common_50 - 53360*common_51, (18584/45)*common_29 - 12788/9*common_30 - 8946*common_32 - 588*common_34 + (1832/5)*common_35 + (7292/9)*common_36 - 19560*common_37 - 1950*common_38 - 1700*common_39 - 8160*common_40 + 26460*common_52 + 1644*common_53 - 832*common_54 + 6774*common_55 - 342*common_56 + 7140*common_57 + 1020*common_58 + 12066*common_59),
    )

def j4(potential, radial_speed, speed_squared, nu, log_one, log_infrared):
    common_0 = nu*potential
    common_1 = nu*speed_squared
    common_2 = nu**2
    common_3 = common_2*potential
    common_4 = common_2*speed_squared
    common_5 = nu*radial_speed
    common_6 = common_2*radial_speed
    common_7 = (12/11)*common_2 - 16/11*nu + 4/11
    common_8 = potential**2
    common_9 = speed_squared**2
    common_10 = common_1*potential
    common_11 = common_8*nu
    common_12 = common_9*nu
    common_13 = radial_speed**2
    common_14 = common_13*potential
    common_15 = common_0*common_13
    common_16 = nu**3
    common_17 = common_16*common_8
    common_18 = common_16*common_9
    common_19 = common_13*common_3
    common_20 = potential*radial_speed
    common_21 = radial_speed*speed_squared
    common_22 = common_0*radial_speed
    common_23 = common_1*radial_speed
    common_24 = common_3*radial_speed
    common_25 = common_4*radial_speed
    common_26 = common_16*common_20
    common_27 = common_16*common_21
    common_28 = common_13*nu
    common_29 = common_16*potential
    common_30 = common_16*speed_squared
    common_31 = common_13*common_2
    common_32 = common_13*common_16
    common_33 = common_16*radial_speed
    common_34 = 6*nu
    common_35 = potential*speed_squared
    common_36 = common_35*radial_speed
    common_37 = radial_speed**3
    common_38 = common_37*potential
    common_39 = common_8*radial_speed
    common_40 = common_1*common_20
    common_41 = common_0*common_37
    common_42 = common_5*common_8
    common_43 = common_21*common_3
    common_44 = common_33*common_35
    common_45 = common_3*common_37
    common_46 = common_6*common_8
    common_47 = common_29*common_37
    common_48 = common_17*radial_speed
    common_49 = 6 - 12*nu
    common_50 = common_3*speed_squared
    common_51 = common_2*common_8
    common_52 = potential**3
    common_53 = common_52*nu
    common_54 = radial_speed**4
    common_55 = common_54*potential
    common_56 = common_9*potential
    common_57 = common_8*speed_squared
    common_58 = common_0*common_54
    common_59 = common_0*common_9
    common_60 = common_1*common_8
    common_61 = common_14*speed_squared
    common_62 = common_2*common_52
    common_63 = common_16*common_52
    common_64 = common_13*common_8
    common_65 = common_1*common_14
    common_66 = common_11*common_13
    common_67 = common_3*common_54
    common_68 = common_3*common_9
    common_69 = common_4*common_8
    common_70 = common_29*common_54
    common_71 = common_18*potential
    common_72 = common_17*speed_squared
    common_73 = common_19*speed_squared
    common_74 = common_14*common_30
    common_75 = common_13*common_51
    common_76 = common_13*common_17
    common_77 = 21*potential
    common_78 = radial_speed**5
    common_79 = common_0*common_78
    common_80 = common_5*common_52
    common_81 = common_35*common_37
    common_82 = common_21*common_8
    common_83 = common_37*common_8
    common_84 = common_22*common_9
    common_85 = common_1*common_38
    common_86 = common_3*common_78
    common_87 = common_29*common_78
    common_88 = common_33*common_52
    common_89 = common_39*common_4
    common_90 = common_18*common_20
    common_91 = common_17*common_37
    common_92 = 120*common_0 - 60*potential
    common_93 = potential**4
    common_94 = radial_speed**6
    common_95 = common_94*potential
    common_96 = speed_squared**3
    common_97 = common_0*common_94
    common_98 = common_0*common_96
    common_99 = common_1*common_52
    common_100 = common_2*common_93
    common_101 = common_16*common_93
    common_102 = common_8*common_9
    common_103 = common_13*common_52
    common_104 = common_11*common_54
    common_105 = common_12*common_8
    common_106 = common_29*common_96
    common_107 = common_30*common_52
    common_108 = common_14*common_9
    common_109 = common_64*speed_squared
    common_110 = common_15*common_9
    common_111 = common_50*common_54
    common_112 = common_16*common_35*common_54
    common_113 = common_51*common_54
    common_114 = common_51*common_9
    common_115 = common_31*common_52
    common_116 = common_13*common_72
    common_117 = common_1*common_39
    common_118 = common_24*common_9
    common_119 = common_37*common_51
    common_120 = common_78*potential
    common_121 = common_52*radial_speed
    common_122 = common_20*common_9
    common_123 = radial_speed**7
    return (
        (1 - 2*nu, 0, 0, 0),
        ((35/44)*common_0 - 95/22*common_1 - 185/22*common_3 + (195/22)*common_4 + (15/11)*potential + (5/11)*speed_squared, -10/11*common_5 + (15/22)*common_6 + (5/22)*radial_speed, common_7, 0),
        (-55157/2860*common_10 - 1731/130*common_11 - 6311/1144*common_12 + (55/4)*common_13*common_16*potential - 3/40*common_14 - 111/220*common_15 + (31475/572)*common_16*potential*speed_squared - 6353/286*common_17 - 11961/286*common_18 - 159/44*common_19 + (24667/572)*common_2*common_8 + (30885/1144)*common_2*common_9 + (6231/286)*common_2*potential*speed_squared - 7564/5005*common_8 + (391/1144)*common_9 + (13229/5720)*potential*speed_squared, (6077/2860)*common_20 + (71/572)*common_21 + (3739/1430)*common_22 - 881/572*common_23 - 15107/572*common_24 + (2985/572)*common_25 + (2231/286)*common_26 - 2559/572*common_27, -13867/2860*common_0 - 401/143*common_1 + (21/143)*common_13 - 126/143*common_28 + (1028/143)*common_29 - 3689/572*common_3 - 96/13*common_30 + (210/143)*common_31 - 84/143*common_32 + (1275/143)*common_4 + (4801/1430)*potential + (73/286)*speed_squared, (135/143)*common_2*radial_speed - 54/143*common_33 - 81/143*common_5 + (27/286)*radial_speed),
        (0, 3 - common_34, 0, 0),
        (-5/4*common_22 - 6*common_24 + (3/2)*potential*radial_speed, (273/44)*common_0 - 305/22*common_1 - 309/11*common_3 + (300/11)*common_4 + (69/22)*potential + (35/22)*speed_squared, -52/11*common_5 + (39/11)*common_6 + (13/11)*radial_speed, common_7),
        ((2083/5720)*common_36 + (9/40)*common_38 + (88533/20020)*common_39 + (1393/1430)*common_40 - 141/110*common_41 + (8324/715)*common_42 - 7277/572*common_43 + (4551/143)*common_44 + (111/44)*common_45 - 2073/52*common_46 - 84/11*common_47 - 11237/286*common_48, -67409/1430*common_10 - 1017/22*common_11 - 20695/1144*common_12 + (6260/143)*common_13*common_16*potential + (1445/143)*common_13*common_2*potential - 941/5720*common_14 - 14481/1430*common_15 + (26245/143)*common_16*potential*speed_squared - 914/11*common_17 - 74325/572*common_18 + (345/2)*common_2*common_8 + (98625/1144)*common_2*common_9 + (3502/143)*common_2*potential*speed_squared - 2943/308*common_8 + (1315/1144)*common_9 + (42751/5720)*potential*speed_squared, (15563/1430)*common_20 + (301/286)*common_21 - 29131/2860*common_22 - 2989/286*common_23 - 34923/572*common_24 + (8925/286)*common_25 + (4765/143)*common_26 - 7119/286*common_27, -12247/2860*common_0 - 482/143*common_1 + (48/143)*common_13 - 288/143*common_28 + (1082/143)*common_29 - 4229/572*common_3 - 1110/143*common_30 + (480/143)*common_31 - 192/143*common_32 + (1410/143)*common_4 + (2333/715)*potential + (50/143)*speed_squared),
        (common_34*potential - 3*potential, 0, common_49, 0),
        ((423/44)*common_10 - 505/22*common_11 - 9/2*common_14 + (33/4)*common_15 + 9*common_19 - 34/11*common_35 - 168/11*common_50 + (243/11)*common_51 + (81/11)*common_8, -721/22*common_22 - 303/11*common_24 + (218/11)*potential*radial_speed, (43/2)*common_0 - 357/11*common_1 - 63*common_3 + (639/11)*common_4 + 4*potential + (48/11)*speed_squared, -84/11*common_5 + (63/11)*common_6 + (21/11)*radial_speed),
        (-33249/10010*common_52 + (149939/2860)*common_53 - 9/8*common_55 - 1699/715*common_56 - 16903/10010*common_57 + (69/88)*common_58 + (109567/5720)*common_59 + (4412/715)*common_60 - 597/1430*common_61 - 73553/572*common_62 + (15501/286)*common_63 - 18573/20020*common_64 + (26079/2860)*common_65 + (27757/1430)*common_66 + (1365/88)*common_67 - 62479/1144*common_68 + (566/143)*common_69 + (195/44)*common_70 + (32001/572)*common_71 - 24547/286*common_72 - 11475/572*common_73 - 5931/143*common_74 - 1764/143*common_75 + (9583/286)*common_76, (9289/1430)*common_36 + (1671/1430)*common_38 + (125183/4004)*common_39 - 65679/1430*common_40 + (11229/1430)*common_41 + (17943/286)*common_42 + (3613/52)*common_43 + (20298/143)*common_44 - 2211/52*common_45 - 152747/572*common_46 - 3687/143*common_47 - 54143/286*common_48, -233457/2860*common_10 - 72387/715*common_11 - 26673/572*common_12 + (25713/286)*common_13*common_16*potential + (51303/572)*common_13*common_2*potential + (13113/2860)*common_13*potential - 185343/2860*common_15 + (122439/286)*common_16*potential*speed_squared - 29825/143*common_17 - 3702/13*common_18 + (125075/286)*common_2*common_8 + (116475/572)*common_2*common_9 - 29977/572*common_50 - 153611/5005*common_8 + (1917/572)*common_9 + (61367/2860)*potential*speed_squared, (26589/1430)*common_20 + (63/26)*common_21 - 42729/1430*common_22 - 6069/286*common_23 - 9408/143*common_24 + (16485/286)*common_25 + (7203/143)*common_26 - 12327/286*common_27),
        (-18*common_22 + 9*potential*radial_speed, 42*common_0 - common_77, 0, common_49),
        ((702/11)*common_2*potential*radial_speed*speed_squared + (45/2)*common_37*potential - 507/11*common_39 - 543/44*common_40 - 165/4*common_41 - 45*common_45 - 663/11*common_46 + (2615/22)*common_8*nu*radial_speed + (3/11)*potential*radial_speed*speed_squared, (2155/44)*common_10 - 1662/11*common_11 - 1605/22*common_14 + (6207/44)*common_15 + (810/11)*common_19 - 112/11*common_35 - 1293/11*common_50 + (1890/11)*common_51 + (465/11)*common_8, -3943/22*common_22 - 588/11*common_24 + (1025/11)*potential*radial_speed, (641/22)*common_0 - 441/11*common_1 - 756/11*common_3 + (702/11)*common_4 + (23/11)*potential + (69/11)*speed_squared),
        ((32190/143)*common_16*common_37*potential*speed_squared + (4569/13)*common_16*common_8*radial_speed*speed_squared + (39909/286)*common_2*common_37*common_8 + (92865/572)*common_2*common_37*potential*speed_squared + (113765/143)*common_2*common_52*radial_speed + (12867/104)*common_2*common_9*potential*radial_speed + (219/286)*common_37*common_8*nu + (892413/20020)*common_52*radial_speed + (63/8)*common_78*potential - 483/88*common_79 + (85451/1430)*common_8*nu*radial_speed*speed_squared - 1178279/2860*common_80 - 345/143*common_81 - 165086/5005*common_82 - 88623/2002*common_83 - 44877/1144*common_84 - 24285/572*common_85 - 9555/88*common_86 - 1365/44*common_87 - 17933/143*common_88 - 12659/572*common_89 + (900/143)*common_9*potential*radial_speed - 13041/52*common_90 - 69261/286*common_91, (468999/20020)*common_52 + (996903/2860)*common_53 - 10545/1144*common_55 - 2297/143*common_56 - 322727/20020*common_57 - 67965/1144*common_58 + (123009/1144)*common_59 + (23289/715)*common_60 - 4929/286*common_61 - 152310/143*common_62 + (66059/143)*common_63 + (433673/4004)*common_64 + (146313/572)*common_65 - 60993/286*common_66 + (425145/1144)*common_67 - 346715/1144*common_68 + (16571/286)*common_69 + (375/52)*common_70 + (248691/572)*common_71 - 100666/143*common_72 - 385269/572*common_73 - 44235/143*common_74 + (41699/286)*common_75 + (51247/143)*common_76, (35889/143)*common_16*potential*radial_speed*speed_squared + (438023/572)*common_2*potential*radial_speed*speed_squared + (86259/572)*common_37*nu*potential - 6531/572*common_38 - 230475/572*common_40 - 216123/572*common_45 - 328073/286*common_46 - 1893/143*common_47 - 51466/143*common_48 + (92979/286)*common_8*nu*radial_speed + (85951/1001)*common_8*radial_speed + (30229/572)*potential*radial_speed*speed_squared, -50553/572*common_10 - 28053/286*common_11 - 38811/572*common_12 + (11019/143)*common_13*common_16*potential + (125643/572)*common_13*common_2*potential + (1005/52)*common_13*potential - 92343/572*common_15 + (72129/143)*common_16*potential*speed_squared - 35390/143*common_17 - 93771/286*common_18 + (146075/286)*common_2*common_8 + (149445/572)*common_2*common_9 - 90751/572*common_50 - 83381/2002*common_8 + (3303/572)*common_9 + (20431/572)*potential*speed_squared),
        (-18*common_10 - 24*common_11 + 90*common_13*nu*potential - 45*common_14 + 12*common_8 + 9*potential*speed_squared, -180*common_22 + 90*potential*radial_speed, common_92, 0),
        ((4437/11)*common_13*common_2*common_8 + (3888/11)*common_13*common_8 + (1455/22)*common_13*potential*speed_squared + 315*common_2*common_54*potential + (702/11)*common_2*common_9*potential + (2293/22)*common_52*nu - 486/11*common_52 + (1155/4)*common_54*nu*potential - 315/2*common_55 - 266/11*common_57 - 543/44*common_59 - 699/11*common_62 - 1365/22*common_65 - 9942/11*common_66 - 864/11*common_69 - 4995/11*common_73 + (2073/22)*common_8*nu*speed_squared + (3/11)*common_9*potential, (7605/11)*common_2*potential*radial_speed*speed_squared - 1260/11*common_36 + (4755/11)*common_37*potential - 5450/11*common_39 - 9120/11*common_41 - 5535/11*common_45 - 9255/11*common_46 + (29905/22)*common_8*nu*radial_speed + (1080/11)*nu*potential*radial_speed*speed_squared, common_1*common_77 - 3502/11*common_11 - 4680/11*common_14 + (9315/11)*common_15 + (2790/11)*common_19 + 36*common_35 - 372*common_50 + (5844/11)*common_51 + (628/11)*common_8, -4050/11*common_22 - 180/11*common_24 + (1920/11)*potential*radial_speed),
        ((308595/1144)*common_1*common_55 - 313998/715*common_1*common_64 + (1422/11)*common_100 - 1248/11*common_101 - 66596/5005*common_102 - 6706059/10010*common_103 - 105852/143*common_104 + (14693/715)*common_105 - 13041/52*common_106 - 263*common_107 - 5535/143*common_108 + (368736/5005)*common_109 + (78675/1144)*common_110 - 1921185/1144*common_111 - 990045/572*common_112 - 958281/572*common_113 - 493/26*common_114 - 3204741/572*common_115 - 458181/143*common_116 + (1103535/572)*common_14*common_18 + (668517/286)*common_17*common_54 + (60726/143)*common_18*common_8 - 150495/1144*common_19*common_9 + (10057549/2860)*common_28*common_52 + (12285/44)*common_29*common_94 + (85995/88)*common_3*common_94 + (12867/104)*common_3*common_96 + (182470/143)*common_32*common_52 + (64365/1144)*common_35*common_54 + (3133/4)*common_4*common_52 + (377409/286)*common_4*common_64 + (22277/385)*common_52*speed_squared + (1378989/2002)*common_54*common_8 - 1547/55*common_93*nu + (9519/385)*common_93 - 567/8*common_95 + (900/143)*common_96*potential + (4347/88)*common_97 - 44877/1144*common_98 - 90619/220*common_99, -32312/715*common_117 - 9540/143*common_118 - 1407/143*common_119 + (321870/143)*common_16*common_37*potential*speed_squared + (665970/143)*common_16*common_8*radial_speed*speed_squared + (1527615/286)*common_2*common_37*potential*speed_squared + (1328214/143)*common_2*common_52*radial_speed + (1077129/286)*common_37*common_8*nu + (6015/143)*common_37*potential*speed_squared + (815972/5005)*common_52*radial_speed + (20769/52)*common_78*nu*potential + (50421/572)*common_78*potential - 248671/65*common_80 - 39176/5005*common_82 - 1857615/1001*common_83 - 470175/286*common_85 - 837165/286*common_86 - 41055/286*common_87 - 384344/143*common_88 - 329155/572*common_89 + (40797/572)*common_9*nu*potential*radial_speed + (4662/143)*common_9*potential*radial_speed - 765153/286*common_90 - 970575/286*common_91, (1610372/5005)*common_52 + (307759/715)*common_53 + (1005/26)*common_55 - 2205/143*common_56 - 653456/5005*common_57 - 518565/572*common_58 + (14583/572)*common_59 + (375598/715)*common_60 - 32499/143*common_61 - 38409/13*common_62 + (215498/143)*common_63 + (1271916/1001)*common_64 + (652653/286)*common_65 - 1136871/286*common_66 + (1602285/572)*common_67 - 157191/572*common_68 - 75886/143*common_69 - 17430/143*common_70 + (18015/13)*common_71 - 328638/143*common_72 - 1495305/286*common_73 - 129861/143*common_74 + (972357/286)*common_75 + (204366/143)*common_76, (13020/143)*common_16*common_37*potential + (27720/13)*common_2*potential*radial_speed*speed_squared + (78069/143)*common_37*nu*potential - 9924/143*common_38 - 168873/143*common_40 - 4656/143*common_44 - 13170/13*common_45 - 682155/286*common_46 - 4398/143*common_48 + (20247/26)*common_8*nu*radial_speed + (94803/1001)*common_8*radial_speed + (24408/143)*potential*radial_speed*speed_squared),
        (-135*common_36 + 315*common_38 - 90*common_39 + 270*common_40 - 630*common_41 + 180*common_42, -234*common_10 - 132*common_11 + 1170*common_13*nu*potential - 585*common_14 + 66*common_8 + 117*potential*speed_squared, -720*common_22 + 360*potential*radial_speed, common_92),
        ((92463/11)*common_11*common_37 - 28443/11*common_117 - 13500/11*common_118 - 40365/11*common_119 + (2835/2)*common_120 + (2928/11)*common_121 + (1440/11)*common_122 + (48825/11)*common_37*common_50 + (6738/11)*common_52*common_6 - 10395/4*common_79 - 7675/11*common_80 - 24045/22*common_81 + (10155/11)*common_82 - 36324/11*common_83 - 2745/44*common_84 + (34965/22)*common_85 - 2835*common_86 + (16011/11)*common_89, (77526/11)*common_13*common_2*common_8 + (53859/11)*common_13*common_8 + (45495/22)*common_13*potential*speed_squared + (49140/11)*common_2*common_54*potential + (9711/11)*common_2*common_9*potential + (1416/11)*common_52 - 2529/11*common_53 + (293475/44)*common_54*nu*potential - 76965/22*common_55 - 1251/11*common_56 - 5450/11*common_57 - 3210/11*common_62 - 69615/22*common_65 - 143166/11*common_66 - 13248/11*common_69 - 69615/11*common_73 + (16850/11)*common_8*nu*speed_squared + (2691/44)*common_9*nu*potential, 3006*common_2*potential*radial_speed*speed_squared - 1188*common_36 + (32910/11)*common_37*potential - 1684*common_39 - 64815/11*common_41 - 25020/11*common_45 - 4242*common_46 + 5161*common_8*nu*radial_speed + 1827*nu*potential*radial_speed*speed_squared, -3819/11*common_10 + (548/11)*common_11 - 10440/11*common_14 + (21465/11)*common_15 + (3330/11)*common_19 + (2316/11)*common_35 - 4272/11*common_50 + (6024/11)*common_51 - 1292/11*common_8),
        (-4347/8*common_0*common_123 - 2438289/1144*common_1*common_120 - 1494531/715*common_1*common_83 - 3023991/2860*common_102*common_5 - 5330007/286*common_121*common_4 - 12285/4*common_123*common_29 - 85995/8*common_123*common_3 + (6237/8)*common_123*potential + (10841859/286)*common_16*common_37*common_8*speed_squared + (830652/143)*common_16*common_52*radial_speed*speed_squared + (9868635/572)*common_16*common_78*potential*speed_squared + (201482/143)*common_16*common_93*radial_speed + (2924325/572)*common_16*common_96*potential*radial_speed - 7388595/286*common_17*common_78 - 11684925/572*common_18*common_38 - 250641/26*common_18*common_39 + (16696905/286)*common_2*common_37*common_52 + (1052985/52)*common_2*common_78*common_8 + (23998275/1144)*common_2*common_78*potential*speed_squared + (1074393/286)*common_2*common_8*common_9*radial_speed - 15570/143*common_20*common_96 - 1029758/455*common_21*common_52 - 1008675/1144*common_24*common_96 - 1065771/1144*common_35*common_78 + (8434107/910)*common_37*common_52 - 117462057/2860*common_37*common_53 - 196596/13*common_37*common_63 + (34451679/10010)*common_37*common_8*speed_squared + (683655/1144)*common_37*common_9*nu*potential + (141855/286)*common_37*common_9*potential - 6502839/286*common_4*common_83 - 6631275/1144*common_45*common_9 - 3764069/1430*common_5*common_93 + (32324711/2860)*common_52*nu*radial_speed*speed_squared - 172426/143*common_6*common_93 + (4001859/286)*common_78*common_8*nu - 39293181/4004*common_78*common_8 + (1105443/5005)*common_8*common_9*radial_speed + (5770539/5005)*common_93*radial_speed + (381735/1144)*common_96*nu*potential*radial_speed, -281196/143*common_100 - 32066/143*common_101 - 34844/5005*common_102 - 4456427/770*common_103 - 7290732/143*common_104 - 81227/2860*common_105 - 1960659/572*common_106 - 632966/143*common_107 - 21870/143*common_108 - 2929497/385*common_109 - 5814045/1144*common_110 - 65280075/1144*common_111 - 12393045/572*common_112 - 2266947/286*common_113 - 391195/572*common_114 - 2042831/22*common_115 - 566991/11*common_116 + (319784/11)*common_13*common_16*common_52 + (14824575/572)*common_13*common_16*common_9*potential + (296751/22)*common_13*common_2*common_8*speed_squared + (18261495/1144)*common_13*common_2*common_9*potential + (10063677/220)*common_13*common_52*nu + (836931/55)*common_13*common_8*nu*speed_squared + (5459913/143)*common_16*common_54*common_8 + (1782591/286)*common_16*common_8*common_9 + (1218105/572)*common_16*common_94*potential + (3262005/286)*common_2*common_52*speed_squared + (33491745/1144)*common_2*common_94*potential + (348291/1144)*common_2*common_96*potential + (4087373/5005)*common_52*speed_squared + (99922899/4004)*common_54*common_8 + (16375275/1144)*common_54*nu*potential*speed_squared + (360465/1144)*common_54*potential*speed_squared + (1464543/286)*common_93*nu - 1989441/1001*common_93 - 1150821/1144*common_95 + (7362/143)*common_96*potential - 3942729/1144*common_97 - 53037/1144*common_98 - 15343871/2860*common_99, -807316/65*common_117 - 5585967/572*common_118 - 3096876/143*common_119 - 13482/143*common_120 - 1471853/385*common_121 - 49059/143*common_122 + (1223325/143)*common_16*common_37*potential*speed_squared + (80955/143)*common_16*common_78*potential + (3293790/143)*common_16*common_8*radial_speed*speed_squared + (13736325/286)*common_2*common_37*potential*speed_squared + (803233/22)*common_2*common_52*radial_speed + (1719727/286)*common_2*common_8*radial_speed*speed_squared + (7210821/143)*common_37*common_8*nu + (196635/143)*common_37*potential*speed_squared + (4086873/572)*common_78*nu*potential + (20156662/5005)*common_8*radial_speed*speed_squared - 486127/110*common_80 - 19355706/1001*common_83 - 5240745/286*common_85 - 14564655/572*common_86 - 179072/11*common_88 + (2648457/572)*common_9*nu*potential*radial_speed - 1619370/143*common_90 - 2484837/143*common_91, (4630757/5005)*common_52 - 2104867/1430*common_53 + (110295/286)*common_55 + (22203/143)*common_56 - 146717/385*common_57 - 2079945/572*common_58 - 660909/572*common_59 + (256247/110)*common_60 - 135495/143*common_61 - 392163/286*common_62 + (215216/143)*common_63 + (4253796/1001)*common_64 + (2134305/286)*common_65 - 4350981/286*common_66 + (4499685/572)*common_67 + (1062489/572)*common_68 - 86399/22*common_69 - 82530/143*common_70 + (193509/143)*common_71 - 24720/11*common_72 - 4194045/286*common_73 - 76833/143*common_74 + (4151967/286)*common_75 + (212100/143)*common_76),
    )

def j5(potential, radial_speed, speed_squared, nu, log_one, log_infrared):
    common_0 = nu**2
    common_1 = nu*potential
    common_2 = nu*speed_squared
    common_3 = common_0*potential
    common_4 = common_0*speed_squared
    common_5 = nu**3
    common_6 = common_5*potential
    common_7 = common_5*speed_squared
    common_8 = (20/13)*radial_speed
    common_9 = 20*common_0
    common_10 = potential*radial_speed
    common_11 = common_1*radial_speed
    common_12 = common_3*radial_speed
    common_13 = common_6*radial_speed
    common_14 = (168/13)*radial_speed
    common_15 = (216/13)*common_0 - 108/13*common_5 - 108/13*nu + 108/91
    common_16 = 60*nu
    common_17 = 60*common_0
    common_18 = potential*speed_squared
    common_19 = potential**2
    common_20 = common_19*nu
    common_21 = radial_speed**2
    common_22 = common_21*potential
    common_23 = common_3*speed_squared
    common_24 = common_19*common_5
    common_25 = common_21*common_3
    common_26 = common_21*common_6
    common_27 = (552/13)*radial_speed
    common_28 = 120*common_0 - 120*nu + 24
    common_29 = common_10*speed_squared
    common_30 = radial_speed**3
    common_31 = common_30*potential
    common_32 = common_19*radial_speed
    common_33 = common_10*common_2
    common_34 = common_1*common_30
    common_35 = common_20*radial_speed
    common_36 = common_12*speed_squared
    common_37 = common_10*common_7
    common_38 = common_3*common_30
    common_39 = common_0*radial_speed
    common_40 = common_19*common_39
    common_41 = common_30*common_6
    common_42 = common_24*radial_speed
    common_43 = (768/13)*radial_speed
    common_44 = common_1*common_21
    common_45 = common_0*common_19
    common_46 = potential**3
    common_47 = radial_speed**4
    common_48 = common_47*potential
    common_49 = speed_squared**2
    common_50 = common_19*speed_squared
    common_51 = common_1*common_49
    common_52 = common_0*common_46
    common_53 = common_2*common_22
    common_54 = common_20*common_21
    common_55 = common_3*common_47
    common_56 = common_19*common_4
    common_57 = common_47*common_6
    common_58 = common_49*common_6
    common_59 = common_21*common_23
    common_60 = common_21*common_24
    common_61 = common_2*potential
    common_62 = -2400*common_3 + 2400*nu*potential - 480*potential
    common_63 = radial_speed**5
    common_64 = common_1*common_63
    common_65 = common_46*nu
    common_66 = common_65*radial_speed
    common_67 = common_18*common_30
    common_68 = common_19*common_30
    common_69 = common_11*common_49
    common_70 = common_2*common_32
    common_71 = common_46*common_5
    common_72 = common_71*radial_speed
    common_73 = common_12*common_49
    common_74 = common_23*common_30
    common_75 = common_31*common_7
    common_76 = common_32*common_7
    common_77 = common_30*common_45
    common_78 = common_49*potential
    common_79 = common_3*common_49
    common_80 = potential**4
    common_81 = radial_speed**6
    common_82 = speed_squared**3
    common_83 = common_19*common_21
    return (
        (5*common_0 - 5*nu + 1, 0, 0, 0, 0),
        (-138/65*common_1 - 161/26*common_2 - 210/13*common_3 + (317/13)*common_4 + (339/13)*common_6 - 707/26*common_7 + (81/65)*potential + (83/182)*speed_squared, (40/13)*common_0*radial_speed - common_5*common_8 - common_8*nu + (20/91)*radial_speed, (108/13)*common_0 - 54/13*common_5 - 54/13*nu + 54/91, 0, 0),
        (0, common_9 - 20*nu + 4, 0, 0, 0),
        ((738/455)*common_10 - 387/65*common_11 - 74/13*common_12 + (258/13)*common_13, (88/65)*common_1 - 342/13*common_2 - 1096/13*common_3 + (1308/13)*common_4 + (1484/13)*common_6 - 1434/13*common_7 + (1628/455)*potential + (186/91)*speed_squared, (336/13)*common_0*radial_speed - common_14*common_5 - common_14*nu + (24/13)*radial_speed, common_15, 0),
        (-common_9*potential + 20*nu*potential - 4*potential, 0, -common_16 + common_17 + 12, 0, 0),
        (130*common_0*common_19 - 2012/455*common_18 + (54/5)*common_19 - 337/5*common_20 + (1551/65)*common_21*nu*potential - 2214/455*common_22 - 862/13*common_23 - 94*common_24 - 168/13*common_25 - 384/13*common_26 + (912/13)*common_5*potential*speed_squared + (1843/65)*nu*potential*speed_squared, (12344/455)*common_10 - 7976/65*common_11 + (632/13)*common_12 + (1712/13)*common_13, (2724/65)*common_1 - 1194/13*common_2 - 4272/13*common_3 + (4260/13)*common_4 + (4944/13)*common_6 - 4470/13*common_7 + (2424/455)*potential + (726/91)*speed_squared, (1104/13)*common_0*radial_speed - common_27*common_5 - common_27*nu + (552/91)*radial_speed, common_15),
        (-common_10*common_16 + common_10*common_17 + 12*common_10, -200*common_3 + 200*nu*potential - 40*potential, 0, common_28, 0),
        ((1608/455)*common_29 + (2214/91)*common_31 - 4404/65*common_32 - 2427/65*common_33 - 1551/13*common_34 + (24432/65)*common_35 + (2250/13)*common_36 - 3504/13*common_37 + (840/13)*common_38 - 6892/13*common_40 + (1920/13)*common_41 + (2640/13)*common_42, (17792/13)*common_0*common_19 - 1068/35*common_18 + (6592/65)*common_19 - 43076/65*common_20 + (32472/65)*common_21*nu*potential - 45888/455*common_22 - 632*common_23 - 13368/13*common_24 - 4908/13*common_25 - 4332/13*common_26 + 740*common_5*potential*speed_squared + (1112/5)*nu*potential*speed_squared, (12084/65)*common_10 - 57912/65*common_11 + (8376/13)*common_12 + (6108/13)*common_13, (10368/65)*common_1 - 2940/13*common_2 - 10512/13*common_3 + (9624/13)*common_4 + (10872/13)*common_6 - 9492/13*common_7 - 72/455*potential + (2004/91)*speed_squared, (1536/13)*common_0*radial_speed - common_43*common_5 - common_43*nu + (768/91)*radial_speed),
        (-common_16*common_18 + common_17*common_18 + 12*common_18 + 28*common_19 - 140*common_20 - 60*common_22 - 300*common_25 + 300*common_44 + 140*common_45, 168*common_10 - 840*common_11 + 840*common_12, -960*common_3 + 960*nu*potential - 192*potential, 0, common_28),
        ((51330/13)*common_0*common_19*common_21 + (2250/13)*common_0*common_49*potential + (47622/91)*common_19*common_21 + (1984/13)*common_19*common_5*speed_squared + (8763/65)*common_19*nu*speed_squared + (23280/13)*common_21*common_5*potential*speed_squared + (5034/91)*common_21*potential*speed_squared + (7088/13)*common_46*common_5 + (51404/65)*common_46*nu - 9468/65*common_46 + (10857/13)*common_47*nu*potential - 2214/13*common_48 + (1608/455)*common_49*potential - 5812/455*common_50 - 2427/65*common_51 - 14540/13*common_52 - 2226/13*common_53 - 37599/13*common_54 - 5880/13*common_55 - 4566/13*common_56 - 13440/13*common_57 - 3504/13*common_58 - 8730/13*common_59 - 16110/13*common_60, -43692/455*common_29 + (54744/91)*common_31 - 100624/91*common_32 + (11868/65)*common_33 - 38676/13*common_34 + (81280/13)*common_35 + (23832/13)*common_36 - 51540/13*common_37 + (27900/13)*common_38 - 120440/13*common_40 + (29340/13)*common_41 + (51640/13)*common_42, (85896/13)*common_0*common_19 - 19884/455*common_18 + (37020/91)*common_19 - 37332/13*common_20 + (278172/65)*common_21*nu*potential - 391428/455*common_22 - 35784/13*common_23 - 69468/13*common_24 - 46872/13*common_25 - 24300/13*common_26 + (49404/13)*common_5*potential*speed_squared + (38916/65)*nu*potential*speed_squared, (23736/35)*common_10 - 16824/5*common_11 + 3024*common_12 + 840*common_13, (14208/65)*common_1 - 3708/13*common_2 - 12048/13*common_3 + (11160/13)*common_4 + (11640/13)*common_6 - 10260/13*common_7 - 3912/455*potential + (396/13)*speed_squared),
        (2100*common_0*common_30*potential + 1200*common_19*nu*radial_speed - 180*common_29 + 420*common_30*potential - 240*common_32 - 2100*common_34 - 900*common_36 - 1200*common_40 + 900*nu*potential*radial_speed*speed_squared, 216*common_18 + 328*common_19 - 1640*common_20 - 1080*common_22 + 1080*common_23 - 5400*common_25 + 5400*common_44 + 1640*common_45 - 1080*common_61, 1080*common_10 - 5400*common_11 + 5400*common_12, common_62, 0),
        ((8460/91)*common_10*common_49 + (64080/13)*common_13*common_49 + (59010/13)*common_2*common_31 + (353376/13)*common_20*common_30 + (137160/13)*common_24*common_30 + (52920/13)*common_3*common_63 + (138084/13)*common_32*common_4 + (530772/455)*common_32*speed_squared + (148972/13)*common_39*common_46 + (684456/455)*common_46*radial_speed + (120960/13)*common_6*common_63 + (19926/13)*common_63*potential - 97713/13*common_64 - 532044/65*common_66 - 13890/13*common_67 - 64188/13*common_68 - 2025/13*common_69 - 450588/65*common_70 - 67284/13*common_72 - 28710/13*common_73 + (37590/13)*common_74 - 216720/13*common_75 - 43848/13*common_76 - 474780/13*common_77, (1107960/13)*common_0*common_19*common_21 + (32832/13)*common_0*common_49*potential + (141072/13)*common_19*common_21 + (54428/13)*common_19*common_5*speed_squared + (323672/65)*common_19*nu*speed_squared + (438840/13)*common_21*common_5*potential*speed_squared + (32580/13)*common_21*potential*speed_squared + 6736*common_46*common_5 + (32712/5)*common_46*nu - 40408/35*common_46 + (314160/13)*common_47*nu*potential - 63600/13*common_48 + (432/13)*common_49*nu*potential - 344628/455*common_50 - 10584*common_52 - 136800/13*common_53 - 786420/13*common_54 - 218820/13*common_55 - 119048/13*common_56 - 259140/13*common_57 - 65556/13*common_58 - 70380/13*common_59 - 397200/13*common_60 - 7452/91*common_78, -24408/13*common_29 + (79380/13)*common_31 - 97176/13*common_32 + (95040/13)*common_33 - 394200/13*common_34 + (562296/13)*common_35 + (85104/13)*common_36 - 351432/13*common_37 + (318060/13)*common_38 - 882168/13*common_40 + (209520/13)*common_41 + (424704/13)*common_42, (193152/13)*common_0*common_19 + (7200/13)*common_19 - 64752/13*common_20 + (244368/13)*common_21*nu*potential - 48816/13*common_22 - 64416/13*common_23 - 183936/13*common_24 - 221040/13*common_25 - 72000/13*common_26 + (132048/13)*common_5*potential*speed_squared - 10848/13*common_61 + (4848/13)*potential*speed_squared, (14376/13)*common_10 - 73488/13*common_11 + (75840/13)*common_12 + (4200/13)*common_13),
        (10500*common_0*common_19*common_21 + 12600*common_0*common_21*potential*speed_squared + 2100*common_19*common_21 + 1380*common_19*nu*speed_squared + 2520*common_21*potential*speed_squared + 440*common_46*nu - 88*common_46 + 18900*common_47*nu*potential - 3780*common_48 + 900*common_49*nu*potential - 276*common_50 - 440*common_52 - 12600*common_53 - 10500*common_54 - 18900*common_55 - 1380*common_56 - 180*common_78 - 900*common_79, 46200*common_0*common_30*potential + 16800*common_19*nu*radial_speed - 3960*common_29 + 9240*common_30*potential - 3360*common_32 - 46200*common_34 - 19800*common_36 - 16800*common_40 + 19800*nu*potential*radial_speed*speed_squared, 1728*common_18 + 1344*common_19 - 6720*common_20 - 8640*common_22 + 8640*common_23 - 43200*common_25 + 43200*common_44 + 6720*common_45 - 8640*common_61, 3600*common_10 - 18000*common_11 + 18000*common_12, common_62),
        ((5299590/13)*common_0*common_19*common_47 + (138642/13)*common_0*common_19*common_49 + (313740/13)*common_0*common_21*common_49*potential + (146056/13)*common_0*common_46*speed_squared + (60/13)*common_0*common_80 - 2025/13*common_1*common_82 + (849474/13)*common_19*common_21*common_5*speed_squared + (9754404/65)*common_19*common_21*nu*speed_squared + (720990/13)*common_19*common_47 + (509352/455)*common_19*common_49 - 435808/65*common_2*common_46 - 78435*common_2*common_48 - 3961365/13*common_20*common_47 - 437943/65*common_20*common_49 + (726552/13)*common_21*common_46*common_5 + (5226612/65)*common_21*common_46*nu - 945804/65*common_21*common_46 + (191205/13)*common_21*common_49*nu*potential - 1723638/65*common_21*common_50 - 1519716/13*common_21*common_52 - 50130/13*common_22*common_49 - 5670*common_23*common_47 - 1519770/13*common_24*common_47 - 49392/13*common_24*common_49 - 1098720/13*common_26*common_49 - 582120/13*common_3*common_81 - 28710/13*common_3*common_82 - 2764272/13*common_4*common_83 - 90344/13*common_46*common_7 + (523792/455)*common_46*speed_squared + 196560*common_47*common_5*potential*speed_squared + 17280*common_47*potential*speed_squared - 8844/13*common_5*common_80 + (64080/13)*common_5*common_82*potential - 1330560/13*common_6*common_81 + (3828/65)*common_80*nu + (144/65)*common_80 + (1074843/13)*common_81*nu*potential - 219186/13*common_81*potential + (8460/91)*common_82*potential, (3714888/13)*common_0*common_19*radial_speed*speed_squared + (1657544/13)*common_0*common_46*radial_speed + (2181060/13)*common_0*common_63*potential + (4069800/13)*common_19*common_30*common_5 + (8754168/13)*common_19*common_30*nu + (15018816/455)*common_19*radial_speed*speed_squared + (2450280/13)*common_30*nu*potential*speed_squared + (6011648/455)*common_46*radial_speed + (1461780/13)*common_49*common_5*potential*radial_speed + (527220/91)*common_49*potential*radial_speed + (2816100/13)*common_5*common_63*potential + (652104/13)*common_63*potential - 3218292/13*common_64 - 4954592/65*common_66 - 538020/13*common_67 - 1579488/13*common_68 - 283860/13*common_69 - 12378984/65*common_70 - 1075096/13*common_72 - 419760/13*common_73 - 232260/13*common_74 - 4975320/13*common_75 - 1405512/13*common_76 - 12109800/13*common_77, 284760*common_1*common_47 + (3423216/65)*common_19*common_2 + (613476/13)*common_19*common_7 + (9875520/13)*common_21*common_45 + (3702240/13)*common_22*common_7 + (457920/13)*common_22*speed_squared + (414048/455)*common_46 - 57420*common_48 - 3926004/455*common_50 + (96336/13)*common_51 - 156864/13*common_52 - 2068200/13*common_53 - 6716160/13*common_54 - 221760*common_55 - 1177968/13*common_56 - 172620*common_57 - 548100/13*common_58 + (317520/13)*common_59 - 4046220/13*common_60 - 143232/65*common_65 + (343008/13)*common_71 - 193212/91*common_78 + (183600/13)*common_79 + (1188036/13)*common_83, -12384*common_29 + (402840/13)*common_31 - 292560/13*common_32 + 54720*common_33 - 2010240/13*common_34 + (1759632/13)*common_35 - 6048*common_36 - 95616*common_37 + (1741320/13)*common_38 - 3018192/13*common_40 + (779040/13)*common_41 + (1697856/13)*common_42, 9024*common_0*common_19 + (11424/13)*common_0*potential*speed_squared + 672*common_19*nu - 552*common_19 + (464832/13)*common_21*nu*potential - 91944/13*common_22 - 14472*common_24 - 448560/13*common_25 - 84600/13*common_26 + (136248/13)*common_5*potential*speed_squared - 84336/13*common_61 + (19224/13)*potential*speed_squared),
    )

def j6(potential, radial_speed, speed_squared, nu, log_one, log_infrared):
    common_0 = nu**2
    common_1 = common_0*speed_squared
    common_2 = nu**3
    common_3 = common_2*potential
    common_4 = nu*radial_speed
    common_5 = common_0*radial_speed
    common_6 = common_2*radial_speed
    common_7 = 15*common_0
    common_8 = potential*radial_speed
    common_9 = common_3*radial_speed
    common_10 = nu*potential
    common_11 = 15*nu
    common_12 = 25*common_0
    common_13 = 10*common_2
    common_14 = potential*speed_squared
    common_15 = potential**2
    common_16 = nu*speed_squared
    common_17 = common_16*potential
    common_18 = common_15*nu
    common_19 = radial_speed**2
    common_20 = common_19*potential
    common_21 = common_10*common_19
    common_22 = common_1*potential
    common_23 = common_3*speed_squared
    common_24 = common_0*common_15
    common_25 = common_15*common_2
    common_26 = common_0*potential
    common_27 = common_19*common_26
    common_28 = common_19*common_3
    common_29 = -50*common_0 + 20*common_2 + 30*nu - 5
    common_30 = 15*potential
    common_31 = 45*radial_speed
    common_32 = common_8*speed_squared
    common_33 = radial_speed**3
    common_34 = common_33*potential
    common_35 = common_15*common_4
    common_36 = common_1*common_8
    common_37 = common_26*common_33
    common_38 = common_3*common_33
    common_39 = common_15*common_6
    common_40 = common_26*radial_speed
    common_41 = -360*common_0 + 480*nu - 120
    common_42 = potential**3
    common_43 = common_42*nu
    common_44 = radial_speed**4
    common_45 = common_44*potential
    common_46 = speed_squared**2
    common_47 = common_46*potential
    common_48 = common_15*speed_squared
    common_49 = common_10*common_44
    common_50 = common_10*common_46
    common_51 = common_15*common_16
    common_52 = common_20*speed_squared
    common_53 = common_0*common_42
    common_54 = common_2*common_42
    common_55 = common_15*common_19
    common_56 = common_16*common_20
    common_57 = common_18*common_19
    common_58 = common_26*common_44
    common_59 = common_26*common_46
    common_60 = common_1*common_15
    common_61 = common_3*common_44
    common_62 = common_3*common_46
    common_63 = common_2*speed_squared
    common_64 = common_15*common_63
    common_65 = common_1*common_20
    common_66 = common_28*speed_squared
    common_67 = common_19*common_24
    common_68 = common_19*common_25
    common_69 = common_15*radial_speed
    common_70 = common_16*common_8
    common_71 = common_10*common_33
    common_72 = common_15*common_5
    common_73 = radial_speed**5
    common_74 = common_73*potential
    common_75 = common_42*radial_speed
    common_76 = common_46*common_8
    common_77 = common_69*speed_squared
    common_78 = common_16*common_34
    common_79 = common_18*common_33
    common_80 = common_26*common_73
    common_81 = common_42*common_5
    common_82 = common_3*common_73
    common_83 = common_1*common_34
    common_84 = common_1*common_69
    common_85 = common_46*common_9
    common_86 = common_25*common_33
    common_87 = 16800*common_18
    common_88 = -16800*common_10 + 12600*common_26 + 4200*potential
    common_89 = potential**4
    common_90 = common_89*nu
    common_91 = radial_speed**6
    common_92 = common_91*potential
    common_93 = speed_squared**3
    common_94 = common_93*potential
    common_95 = common_42*speed_squared
    common_96 = common_10*common_91
    common_97 = common_10*common_93
    common_98 = common_16*common_42
    common_99 = common_14*common_44
    common_100 = common_0*common_89
    common_101 = common_2*common_89
    common_102 = common_15*common_44
    common_103 = common_15*common_46
    common_104 = common_19*common_42
    common_105 = common_16*common_45
    common_106 = common_18*common_44
    common_107 = common_18*common_46
    common_108 = common_19*common_43
    common_109 = common_26*common_91
    common_110 = common_26*common_93
    common_111 = common_1*common_42
    common_112 = common_3*common_91
    common_113 = common_3*common_93
    common_114 = common_42*common_63
    common_115 = common_20*common_46
    common_116 = common_55*speed_squared
    common_117 = common_21*common_46
    common_118 = common_16*common_55
    common_119 = common_1*common_45
    common_120 = common_23*common_44
    common_121 = common_24*common_44
    common_122 = common_24*common_46
    common_123 = common_19*common_53
    common_124 = common_25*common_44
    common_125 = common_25*common_46
    common_126 = common_19*common_54
    common_127 = common_27*common_46
    common_128 = common_1*common_55
    common_129 = common_28*common_46
    common_130 = common_55*common_63
    common_131 = radial_speed**7
    return (
        (-3*common_0 + 4*nu - 1, 0, 0, 0, 0, 0),
        ((35/2)*common_0*potential - 265/12*common_1 + (58/3)*common_2*speed_squared - 56/3*common_3 + (25/4)*nu*speed_squared - 7/6*potential - 11/24*speed_squared, (5/4)*common_4 - 25/12*common_5 + (5/6)*common_6 - 5/24*radial_speed, -25/3*common_0 + (10/3)*common_2 + 5*nu - 5/6, 0, 0, 0),
        (0, -common_7 + 20*nu - 5, 0, 0, 0, 0),
        ((35/4)*common_0*potential*radial_speed - 41/24*common_8 - 89/6*common_9 + (17/4)*nu*potential*radial_speed, (425/4)*common_0*potential - 225/2*common_1 - 45/4*common_10 + (195/2)*common_2*speed_squared - 605/6*common_3 + (65/2)*nu*speed_squared - 95/24*potential - 5/2*speed_squared, common_11*radial_speed - common_12*radial_speed + common_13*radial_speed - 5/2*radial_speed, common_11 - common_12 + common_13 - 5/2, 0, 0),
        (common_7*potential - 20*nu*potential + 5*potential, 0, -60*common_0 + 80*nu - 20, 0, 0, 0),
        ((139/24)*common_14 - 43/3*common_15 - 133/4*common_17 + 77*common_18 + (41/8)*common_20 - 81/4*common_21 + (305/4)*common_22 - 202/3*common_23 - 135*common_24 + (257/3)*common_25 + (15/4)*common_27 + 22*common_28, (25/2)*common_0*potential*radial_speed - 415/12*common_8 - 385/3*common_9 + (235/2)*nu*potential*radial_speed, 525*common_0*potential - 475*common_1 - 105*common_10 + 400*common_2*speed_squared - 1330/3*common_3 + 145*nu*speed_squared - 35/6*potential - 25/2*speed_squared, 75*common_4 - 125*common_5 + 50*common_6 - 25/2*radial_speed, common_29, 0),
        (-common_26*common_31 - common_30*radial_speed + 60*nu*potential*radial_speed, -260*common_10 + 195*common_26 + 65*potential, 0, -180*common_0 + 240*nu - 60, 0, 0),
        ((935/2)*common_0*common_15*radial_speed + (1081/12)*common_15*radial_speed + 246*common_2*potential*radial_speed*speed_squared - 57/8*common_32 + (405/4)*common_33*nu*potential - 205/8*common_34 - 817/2*common_35 - 885/4*common_36 - 75/4*common_37 - 110*common_38 - 461/3*common_39 + (237/4)*nu*potential*radial_speed*speed_squared, (475/8)*common_14 - 2225/12*common_15 - 1515/4*common_17 + (2075/2)*common_18 + (1035/8)*common_20 - 2055/4*common_21 + (3935/4)*common_22 - 905*common_23 - 3795/2*common_24 + (3610/3)*common_25 + (885/4)*common_27 + 315*common_28, -350*common_40 - 310*common_8 - 660*common_9 + 1140*nu*potential*radial_speed, 1900*common_0*potential - 1550*common_1 - 510*common_10 + 1250*common_2*speed_squared - 1460*common_3 + 510*nu*speed_squared + 15*potential - 50*speed_squared, 210*common_4 - 350*common_5 + 140*common_6 - 35*radial_speed, common_29),
        (225*common_0*common_19*potential + 200*common_15*nu - 50*common_15 + 75*common_19*potential - 300*common_21 - 45*common_22 - 150*common_24 - common_30*speed_squared + 60*nu*potential*speed_squared, -810*common_40 - 270*common_8 + 1080*nu*potential*radial_speed, -1760*common_10 + 1320*common_26 + 440*potential, 0, common_41, 0),
        ((886/3)*common_42 - 1329*common_43 + (1435/8)*common_45 - 57/8*common_47 - 73/6*common_48 - 2835/4*common_49 + (237/4)*common_50 - 39*common_51 - 165/4*common_52 + 1630*common_53 - 2249/3*common_54 - 2795/4*common_55 + (15/2)*common_56 + (6315/2)*common_57 + (525/4)*common_58 - 885/4*common_59 + 155*common_60 + 770*common_61 + 246*common_62 + (166/3)*common_64 + 1050*common_65 - 1560*common_66 - 6855/2*common_67 + 850*common_68, (21635/2)*common_0*common_15*radial_speed + (23455/12)*common_15*radial_speed + 4575*common_2*potential*radial_speed*speed_squared + 3075*common_33*nu*potential - 775*common_34 - 18055/2*common_35 - 3615*common_36 - 1200*common_37 - 2125*common_38 - 11585/3*common_39 + 405*nu*potential*radial_speed*speed_squared + 45*potential*radial_speed*speed_squared, (515/2)*common_14 - 3590/3*common_15 - 2085*common_17 + 7060*common_18 + (2895/2)*common_20 - 5745*common_21 + 6615*common_22 - 6410*common_23 - 13660*common_24 + (26320/3)*common_25 + 3015*common_27 + 2430*common_28, -3250*common_40 - 1665*common_8 - 2140*common_9 + 6390*nu*potential*radial_speed, 4400*common_0*potential - 3450*common_1 - 1380*common_10 + 2640*common_2*speed_squared - 3160*common_3 + 1230*nu*speed_squared + 90*potential - 135*speed_squared, -common_31 + 270*common_4 - 450*common_5 + 180*common_6),
        (225*common_32 - 525*common_34 - 1800*common_35 + 675*common_36 - 1575*common_37 + 450*common_69 - 900*common_70 + 2100*common_71 + 1350*common_72, 5175*common_0*common_19*potential - 345*common_14 + 3440*common_15*nu - 860*common_15 + 1725*common_19*potential - 6900*common_21 - 1035*common_22 - 2580*common_24 + 1380*nu*potential*speed_squared, -7200*common_40 - 2400*common_8 + 9600*nu*potential*radial_speed, -7200*common_10 + 5400*common_26 + 1800*potential, 0, common_41),
        (31395*common_0*common_15*common_33 + (12825/4)*common_0*common_46*potential*radial_speed + 549*common_15*common_2*radial_speed*speed_squared + 6630*common_15*common_33 + 6552*common_15*nu*radial_speed*speed_squared + 14000*common_2*common_33*potential*speed_squared + (20096/3)*common_2*common_42*radial_speed + (4025/4)*common_33*potential*speed_squared - 1125/4*common_4*common_47 + (30187/2)*common_42*nu*radial_speed + (25515/4)*common_73*nu*potential - 12915/8*common_74 - 40801/12*common_75 - 375/8*common_76 - 2577/2*common_77 - 5775/2*common_78 - 29790*common_79 - 4725/4*common_80 - 34565/2*common_81 - 6930*common_82 - 6825*common_83 - 7635*common_84 - 4350*common_85 - 6660*common_86, (52265/12)*common_42 - 40555/2*common_43 + (50575/8)*common_45 + (75/8)*common_47 + (2895/4)*common_48 - 100275/4*common_49 + (2805/4)*common_50 - 9695/2*common_51 - 11025/4*common_52 + (54185/2)*common_53 - 41320/3*common_54 - 76485/4*common_55 + (14475/2)*common_56 + (173805/2)*common_57 + (36225/4)*common_58 - 18885/4*common_59 + (14935/2)*common_60 + 18725*common_61 + 5805*common_62 - 830*common_64 + 19725*common_65 - 37050*common_66 - 194205/2*common_67 + 27390*common_68, 114800*common_0*common_15*radial_speed + 19470*common_15*radial_speed + 42390*common_2*potential*radial_speed*speed_squared + 41025*common_33*nu*potential - 20675/2*common_34 - 91620*common_35 - 28275*common_36 - 19875*common_37 - 20650*common_38 - 44080*common_39 - 3615*common_70 + (4605/2)*potential*radial_speed*speed_squared, (255/2)*common_14 - 4205*common_15 - 5265*common_17 + 27030*common_18 + (18675/2)*common_20 - 37125*common_21 + 26075*common_22 - 27610*common_23 - 57250*common_24 + 38220*common_25 + 21675*common_27 + 11550*common_28, -13750*common_40 - 5325*common_8 - 3700*common_9 + 21150*nu*potential*radial_speed, 4850*common_0*potential - 3900*common_1 - 1650*common_10 + 2820*common_2*speed_squared - 3340*common_3 + 1500*nu*speed_squared + 135*potential - 180*speed_squared),
        (common_19*common_87 + 410*common_42 - 1640*common_43 + 4725*common_45 + 225*common_47 + 570*common_48 - 18900*common_49 - 900*common_50 - 2280*common_51 - 3150*common_52 + 1230*common_53 - 4200*common_55 + 12600*common_56 + 14175*common_58 + 675*common_59 + 1710*common_60 - 9450*common_65 - 12600*common_67, 6300*common_32 - 14700*common_34 - 37800*common_35 + 18900*common_36 - 44100*common_37 + 9450*common_69 - 25200*common_70 + 58800*common_71 + 28350*common_72, 56700*common_0*common_19*potential - 3780*common_14 + 25760*common_15*nu - 6440*common_15 + 18900*common_19*potential - 75600*common_21 - 11340*common_22 - 19320*common_24 + 15120*nu*potential*speed_squared, -37800*common_40 - 12600*common_8 + 50400*nu*potential*radial_speed, common_88, 0),
        (-11450*common_100 + (13844/3)*common_101 - 149925/2*common_102 - 1131*common_103 + 38409*common_104 + (231525/4)*common_105 + 335550*common_106 + 6012*common_107 - 171693*common_108 + (51975/4)*common_109 + (12825/4)*common_110 - 21285*common_111 + 76230*common_112 - 4350*common_113 + (29849/3)*common_114 + (26775/8)*common_115 + 33123*common_116 - 26775/4*common_117 - 155016*common_118 + (222075/4)*common_119 - 160650*common_120 - 349050*common_121 - 7200*common_122 + 196590*common_123 + 69825*common_124 + 174*common_125 - 72441*common_126 - 171675/4*common_127 + 167310*common_128 + 72450*common_129 - 23262*common_130 - 7786/3*common_89 + 10924*common_90 + (142065/8)*common_92 - 375/8*common_94 - 8716/3*common_95 - 280665/4*common_96 - 1125/4*common_97 + 14659*common_98 - 137025/8*common_99, 1048110*common_0*common_15*common_33 + (158175/2)*common_0*common_46*potential*radial_speed + 55245*common_15*common_2*radial_speed*speed_squared + 214620*common_15*common_33 + 248400*common_15*nu*radial_speed*speed_squared + 404250*common_2*common_33*potential*speed_squared + (489080/3)*common_2*common_42*radial_speed + (99225/2)*common_33*potential*speed_squared + 296060*common_42*nu*radial_speed + (19125/2)*common_46*nu*potential*radial_speed + (515025/2)*common_73*nu*potential - 259875/4*common_74 - 194020/3*common_75 - 23175/4*common_76 - 103275/2*common_77 - 165375*common_78 - 967680*common_79 - 174825/2*common_80 - 368950*common_81 - 203175*common_82 - 135975*common_83 - 287160*common_84 - 124875*common_85 - 262500*common_86, (67700/3)*common_42 - 112580*common_43 + 97650*common_45 + 2340*common_47 + 16680*common_48 - 387450*common_49 - 810*common_50 - 88600*common_51 - 53550*common_52 + 177320*common_53 - 323020/3*common_54 - 233310*common_55 + 170100*common_56 + 1070160*common_57 + 175350*common_58 - 47160*common_59 + 121720*common_60 + 219450*common_61 + 65610*common_62 - 33960*common_64 + 160650*common_65 - 422100*common_66 - 1236060*common_67 + 387660*common_68, 676900*common_0*common_15*radial_speed + 106400*common_15*radial_speed + 233100*common_2*potential*radial_speed*speed_squared + 308700*common_33*nu*potential - 77700*common_34 - 512400*common_35 - 119700*common_36 - 168000*common_37 - 119700*common_38 - 281400*common_39 - 69300*common_70 + 25200*potential*radial_speed*speed_squared, -3570*common_14 - 6160*common_15 + 2520*common_17 + 49560*common_18 + 34650*common_20 - 138600*common_21 + 52500*common_22 - 67620*common_23 - 127400*common_24 + 93240*common_25 + 88200*common_27 + 31500*common_28, -25200*common_40 - 7980*common_8 - 1680*common_9 + 32760*nu*potential*radial_speed),
        (132300*common_0*common_15*common_33 + 141750*common_0*common_33*potential*speed_squared + 44100*common_15*common_33 + 55440*common_15*nu*radial_speed*speed_squared + 47250*common_33*potential*speed_squared - 23625*common_40*common_46 + 23520*common_42*nu*radial_speed + 31500*common_46*nu*potential*radial_speed + 207900*common_73*nu*potential - 51975*common_74 - 5880*common_75 - 7875*common_76 - 13860*common_77 - 189000*common_78 - 176400*common_79 - 155925*common_80 - 17640*common_81 - 41580*common_84, 5480*common_42 - 21920*common_43 + 155925*common_45 + 7425*common_47 + 13560*common_48 - 623700*common_49 - 29700*common_50 - 54240*common_51 - 103950*common_52 + 16440*common_53 - 102900*common_55 + 415800*common_56 + 411600*common_57 + 467775*common_58 + 22275*common_59 + 40680*common_60 - 311850*common_65 - 308700*common_67, 81900*common_32 - 191100*common_34 - 336000*common_35 + 245700*common_36 - 573300*common_37 + 84000*common_69 - 327600*common_70 + 764400*common_71 + 252000*common_72, 359100*common_0*common_19*potential - 23940*common_14 + 94080*common_15*nu - 23520*common_15 + 119700*common_19*potential - 478800*common_21 - 71820*common_22 - 70560*common_24 + 95760*nu*potential*speed_squared, -113400*common_40 - 37800*common_8 + 151200*nu*potential*radial_speed, common_88),
        ((820665/2)*common_0*common_15*common_46*radial_speed + (9144225/2)*common_0*common_15*common_73 + (2433375/4)*common_0*common_33*common_46*potential + 623595*common_0*common_42*radial_speed*speed_squared + 149800*common_0*common_89*radial_speed - 3323025*common_1*common_15*common_33 - 2130975/4*common_1*common_74 - 743931/2*common_103*common_4 - 56691*common_103*common_6 - 675675/4*common_131*common_26 - 990990*common_131*common_3 + (3648645/4)*common_131*nu*potential - 1846845/8*common_131*potential + 577920*common_15*common_2*common_33*speed_squared + 3135735*common_15*common_33*nu*speed_squared + (313101/4)*common_15*common_46*radial_speed + (3946905/4)*common_15*common_73 - 4230765/4*common_16*common_74 - 478923*common_16*common_75 - 8814555/2*common_18*common_73 + 899430*common_2*common_33*common_42 + 2224530*common_2*common_73*potential*speed_squared + 175350*common_2*common_93*potential*radial_speed - 922005*common_25*common_73 + 2083620*common_33*common_42*nu - 462315*common_33*common_42 + (1167075/4)*common_33*common_46*nu*potential - 1372935/2*common_33*common_48 - 2422560*common_33*common_53 - 789075/8*common_34*common_46 - 1294650*common_38*common_46 - 129272*common_4*common_89 - 45675/4*common_4*common_94 - 433125/4*common_40*common_93 + (202263/2)*common_42*radial_speed*speed_squared - 220532/3*common_6*common_89 - 268416*common_63*common_75 + (2359665/8)*common_73*potential*speed_squared + (90058/3)*common_89*radial_speed + (56175/8)*common_93*potential*radial_speed, -56660*common_100 + (164440/3)*common_101 - 11138775/4*common_102 - 219675/4*common_103 + 838740*common_104 + (12261375/4)*common_105 + (25010175/2)*common_106 + (538515/2)*common_107 - 3860730*common_108 + (4106025/4)*common_109 + (380475/4)*common_110 - 467565*common_111 + 2616075*common_112 - 146625*common_113 + 248760*common_114 + (1649025/8)*common_115 + (2647425/2)*common_116 - 2386125/4*common_117 - 6110235*common_118 + (4257225/4)*common_119 - 5457375*common_120 - 26751675/2*common_121 - 638535/2*common_122 + 4818390*common_123 + 3211950*common_124 + 61110*common_125 - 2075490*common_126 - 4704525/4*common_127 + 6742695*common_128 + 2449125*common_129 - 1466190*common_130 - 21470/3*common_89 + 31740*common_90 + (6427575/8)*common_92 - 48225/8*common_94 - 125065/2*common_95 - 12733875/4*common_96 + (32625/4)*common_97 + 312365*common_98 - 6855975/8*common_99, 15497580*common_0*common_15*common_33 + 873450*common_0*common_46*potential*radial_speed + 1263960*common_15*common_2*radial_speed*speed_squared + 3052560*common_15*common_33 + 3869100*common_15*nu*radial_speed*speed_squared + 5449500*common_2*common_33*potential*speed_squared + 1613180*common_2*common_42*radial_speed + 963900*common_33*potential*speed_squared + 1999980*common_42*nu*radial_speed + 382500*common_46*nu*potential*radial_speed + 4517100*common_73*nu*potential - 1138725*common_74 - 410120*common_75 - 141975*common_76 - 819150*common_77 - 3402000*common_78 - 13860840*common_79 - 1927800*common_80 - 2935700*common_81 - 2787750*common_82 - 967050*common_83 - 4570680*common_84 - 1671750*common_85 - 4426800*common_86, 2740*common_42 - 82380*common_43 + 836850*common_45 + 32220*common_47 + 141320*common_48 - 3323250*common_49 - 71730*common_50 - 714780*common_51 - 519750*common_52 + 381460*common_53 - 389380*common_54 - 1546230*common_55 + 1782900*common_56 + 7192080*common_57 + 1702050*common_58 - 261180*common_59 + 989560*common_60 + 1496250*common_61 + 429930*common_62 - 383700*common_64 + 576450*common_65 - 2790900*common_66 - 8664180*common_67 + 3015180*common_68, 2128000*common_0*common_15*radial_speed + 299180*common_15*radial_speed + 732060*common_2*potential*radial_speed*speed_squared + 1310400*common_33*nu*potential - 328650*common_34 - 1492680*common_35 - 220500*common_36 - 777000*common_37 - 396900*common_38 - 981120*common_39 - 423360*common_70 + 130410*potential*radial_speed*speed_squared, -11550*common_14 + 1820*common_15 + 35280*common_17 + 58590*common_20 - 236880*common_21 + 27300*common_22 - 69300*common_23 - 102200*common_24 + 94920*common_25 + 163800*common_27 + 36540*common_28 + common_87),
    )
