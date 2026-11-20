"""Generated local SH center-of-mass worldline shift."""
import math

GRADES = (2, 4, 5, 6, 7, 8, 9)

def local_common_shift(potential, radial_speed, speed_squared, nu, log_one, log_two):
    common_0 = (1/2)*potential
    common_1 = -common_0 + (1/2)*speed_squared
    common_2 = potential**2
    common_3 = radial_speed**2
    common_4 = common_3*potential
    common_5 = potential*speed_squared
    common_6 = speed_squared**2
    common_7 = (3/2)*nu
    common_8 = common_6*common_7
    common_9 = common_2*nu
    common_10 = (3/4)*nu
    common_11 = common_10*common_4 + common_5*common_7 + (3/8)*common_6 - common_8 - 1/2*common_9
    common_12 = potential*radial_speed
    common_13 = common_5*radial_speed
    common_14 = radial_speed**3
    common_15 = common_2*radial_speed
    common_16 = common_14*potential
    common_17 = (8/5)*common_2
    common_18 = potential**3
    common_19 = common_2*speed_squared
    common_20 = common_2*common_3
    common_21 = common_18*nu
    common_22 = radial_speed**4
    common_23 = common_22*potential
    common_24 = common_6*potential
    common_25 = common_24*nu
    common_26 = common_9*speed_squared
    common_27 = common_23*nu
    common_28 = common_4*speed_squared
    common_29 = common_3*common_9
    common_30 = common_28*nu
    common_31 = speed_squared**3
    common_32 = nu**2
    common_33 = common_31*common_32
    common_34 = common_31*nu
    common_35 = common_18*common_32
    common_36 = common_18*log_two
    common_37 = common_19*common_32
    common_38 = common_20*common_32
    common_39 = (21/16)*common_23*common_32 - 15/2*common_24*common_32 - 11/2*common_28*common_32 + (5/16)*common_31 + 6*common_33 - 11/4*common_34 - 1/2*common_35 + (22/3)*common_36 + 3*common_37 + 4*common_38
    common_40 = common_16*nu
    common_41 = common_9*radial_speed
    common_42 = radial_speed**5
    common_43 = common_42*potential
    common_44 = 22*common_36
    common_45 = common_21*radial_speed
    common_46 = common_12*common_6
    common_47 = common_15*speed_squared
    common_48 = common_14*common_5
    common_49 = common_48*nu
    common_50 = common_14*common_9
    common_51 = common_32*common_43
    common_52 = common_35*radial_speed
    common_53 = common_32*common_46
    common_54 = common_14*common_2
    common_55 = common_32*common_54
    common_56 = potential**4
    common_57 = common_56*nu
    common_58 = common_18*speed_squared
    common_59 = common_21*speed_squared
    common_60 = common_2*common_22
    common_61 = common_2*common_6
    common_62 = common_18*common_3
    common_63 = common_6*common_9
    common_64 = common_21*common_3
    common_65 = common_19*common_3
    common_66 = common_26*common_3
    common_67 = common_56*log_two
    common_68 = radial_speed**6
    common_69 = common_68*potential
    common_70 = common_32*common_56
    common_71 = common_31*potential
    common_72 = common_57*log_one
    common_73 = common_3*common_36
    common_74 = common_57*log_two
    common_75 = common_22*common_9
    common_76 = common_34*potential
    common_77 = common_32*common_69
    common_78 = common_35*speed_squared
    common_79 = common_4*common_6
    common_80 = math.pi**2
    common_81 = common_57*common_80
    common_82 = common_22*common_5
    common_83 = common_69*nu
    common_84 = common_32*common_61
    common_85 = common_33*potential
    common_86 = common_32*common_60
    common_87 = common_3*common_35
    common_88 = 11*log_one
    common_89 = common_64*common_88
    common_90 = common_64*log_two
    common_91 = common_3*common_37
    common_92 = (123/128)*common_80
    common_93 = common_64*common_92
    common_94 = common_82*nu
    common_95 = common_32*common_79
    common_96 = common_32*common_82
    common_97 = speed_squared**4
    common_98 = common_56*log_one
    common_99 = nu**3
    common_100 = common_56*common_99
    common_101 = -1/2*common_100 + (145/8)*common_32*common_97 + 5*common_58*common_99 + 9*common_60*common_99 - 45/2*common_61*common_99 + (55/4)*common_62*common_99 - 83/2*common_65*common_99 + (55/32)*common_69*common_99 + 42*common_71*common_99 + 36*common_79*common_99 - 213/16*common_82*common_99 - 55/2*common_97*common_99 - 125/32*common_97*nu + (35/128)*common_97 + 16*common_98
    common_102 = common_18*radial_speed
    common_103 = common_36*radial_speed
    common_104 = (44/3)*log_two
    common_105 = potential*radial_speed**7
    common_106 = common_56*radial_speed
    common_107 = common_57*radial_speed
    common_108 = common_12*common_31
    common_109 = common_42*common_5
    common_110 = common_102*speed_squared
    common_111 = common_2*common_42
    common_112 = common_14*common_18
    common_113 = common_45*speed_squared
    common_114 = 55*common_14
    common_115 = common_42*common_9
    common_116 = common_14*common_21
    common_117 = common_70*radial_speed
    common_118 = common_16*common_6
    common_119 = common_15*common_6
    common_120 = common_14*common_19
    common_121 = common_41*common_6
    common_122 = common_111*common_32
    common_123 = common_14*common_35
    common_124 = common_119*common_32
    common_125 = potential**5
    return (
        (common_1, 0, -common_0*radial_speed, common_1),
        (common_11 + (7/4)*common_2 - 1/8*common_4 + (19/8)*common_5, -7/4*common_12, -common_10*common_16 - 9/8*common_13 + (3/8)*common_14*potential - common_15*common_7 - 9/4*common_15 + nu*potential*radial_speed*speed_squared, common_11 + (7/2)*common_2 + (13/8)*common_4 + (5/8)*common_5),
        (0, -common_17 + (4/5)*potential*speed_squared, (8/5)*common_18 - 4/5*common_19, 0),
        (-14351/1260*common_18 + (101/12)*common_19 - 7/3*common_20 + (1/8)*common_21 + (1/16)*common_23 + (53/16)*common_24 - 7*common_25 - 33/8*common_26 - 5/8*common_27 - 5/16*common_28 + (73/8)*common_29 + (21/16)*common_30 + common_39, -15/8*common_13 + (5/12)*common_14*potential - 235/24*common_15 - 19/24*common_40 - 21/4*common_41 + (21/4)*nu*potential*radial_speed*speed_squared, (247/24)*common_14*common_2 + (15/4)*common_14*common_32*potential*speed_squared + (19/16)*common_14*potential*speed_squared + (34547/840)*common_18*radial_speed + (11/2)*common_2*common_32*radial_speed*speed_squared + (11/8)*common_2*nu*radial_speed*speed_squared + (5/4)*common_42*nu*potential - 5/16*common_43 - common_44*radial_speed - 27/8*common_45 - 21/16*common_46 - 12*common_47 - 67/16*common_49 - 5/24*common_50 - 15/16*common_51 - 7/2*common_52 - 7/2*common_53 - 9/2*common_55 + (41/8)*common_6*nu*potential*radial_speed, -21667/2520*common_18 + (9/4)*common_19 + (23/4)*common_20 + (15/8)*common_21 - 19/16*common_23 + (23/16)*common_24 - 7/4*common_25 - 75/8*common_26 + (7/4)*common_27 + (45/16)*common_28 + (127/8)*common_29 - 101/16*common_30 + common_39),
        (-common_14*common_17 + (52/15)*common_18*radial_speed + (68/35)*common_2*nu*radial_speed*speed_squared - 116/35*common_45 - 2/15*common_47, -172/105*common_18 + (6/35)*common_19 + (44/15)*common_20 - 64/35*common_21 + (6/7)*common_24 - 22/7*common_25 + (132/35)*common_26 + (68/35)*common_29, -288/35*common_56 + (68/35)*common_57 + (288/35)*common_58 - 52/7*common_59 + 8*common_60 - 188/105*common_61 - 176/15*common_62 + (38/35)*common_63 + (296/21)*common_64 - 22/5*common_65 + (6/35)*common_66, (8/15)*common_18*radial_speed + (176/35)*common_2*nu*radial_speed*speed_squared + (152/15)*common_2*radial_speed*speed_squared - 164/35*common_45 - 24/7*common_50 - 52/5*common_54),
        (common_101 + common_4*common_8 + common_44*speed_squared + (215279/3600)*common_56 + (22043/720)*common_57 + (2737/1440)*common_58 - common_59*common_88 + common_59*common_92 - 87181/3360*common_59 + (1133/960)*common_60 + (5631/320)*common_61 - 185497/8400*common_62 - 139/4*common_63 - 64347/1120*common_64 - 31/5*common_65 + 26*common_66 - 268/3*common_67 - 5/128*common_69 - 3/2*common_70 + (515/128)*common_71 - 60*common_72 - 11*common_73 + 120*common_74 - 1007/48*common_75 - 749/32*common_76 - 125/64*common_77 - 117/8*common_78 - 53/128*common_79 - 11/2*common_81 + (27/128)*common_82 + (35/64)*common_83 + (71/4)*common_84 + (91/4)*common_85 + (169/24)*common_86 + (495/16)*common_87 + common_89 + 33*common_90 - 541/8*common_91 - common_93 - 115/64*common_94 - 95/8*common_95 + (517/64)*common_96, (14377/280)*common_102 - 110/3*common_103 - common_104*common_45 + (95/4)*common_26*radial_speed + 28*common_32*common_47 + (19/4)*common_32*common_48 + (25/32)*common_43*nu - 13/64*common_43 - 41/64*common_45*common_80 + (22/3)*common_45*log_one + (71509/5040)*common_45 + (199/16)*common_46*nu - 123/64*common_46 - 3733/160*common_47 + (77/96)*common_48 - 187/48*common_49 + (113/24)*common_50 - 17/32*common_51 - 49/4*common_52 - 21*common_53 + (4621/480)*common_54 + (7/12)*common_55, -15/2*common_100*radial_speed - 88*common_103*speed_squared + (175/64)*common_105*common_32 - 35/32*common_105*common_99 - 105/64*common_105*nu + (35/128)*common_105 - 1626731/5040*common_106 + (103183/840)*common_107 + 15*common_108*common_99 - 181/128*common_108 - 875/64*common_109*common_32 + (15/2)*common_109*common_99 + (485/64)*common_109*nu - 165/128*common_109 + (39/2)*common_110*common_99 + (111023/4200)*common_110 - 6*common_111*common_99 - 1941/320*common_111 - 73/4*common_112*common_99 + (94043/1120)*common_112 - common_113*common_88 + common_113*common_92 + 66*common_113*log_two - 89191/672*common_113 + common_114*common_21*log_one + common_114*common_36 + (511/16)*common_115 - 615/128*common_116*common_80 - 165*common_116*log_two + (73421/224)*common_116 - 45/4*common_117 + (445/16)*common_118*common_32 - 75/4*common_118*common_99 + (267/128)*common_118 - 27*common_119*common_99 - 1423/64*common_119 - 103/4*common_12*common_33 + (349/32)*common_12*common_34 + 28*common_120*common_99 + (2639/96)*common_120 + (193/4)*common_121 + (115/8)*common_122 + (815/48)*common_123 + (75/8)*common_124 - 495/8*common_14*common_26 - 211/12*common_14*common_37 - 107/8*common_40*common_6 - 5/4*common_52*speed_squared + 372*common_67*radial_speed + (698/3)*common_72*radial_speed - 1594/3*common_74*radial_speed + (1285/64)*common_81*radial_speed - 64*common_98*radial_speed, common_101 - common_104*common_59 - 44/3*common_36*speed_squared - 377077/25200*common_56 + (70087/5040)*common_57 + (79501/1008)*common_58 + (41/128)*common_59*common_80 - 11/3*common_59*log_one + (10063/2016)*common_59 - 125/4*common_60 - 31/16*common_61 - 4701329/16800*common_62 - 285/16*common_63 - 440149/3360*common_64 + (377/8)*common_65 + 16*common_66 - 158/3*common_67 + (125/128)*common_69 + (1/4)*common_70 + (269/128)*common_71 - 202/3*common_72 + 99*common_73 + (404/3)*common_74 - 1759/48*common_75 - 351/32*common_76 + (45/64)*common_77 - 201/8*common_78 - 181/8*common_79*nu + (501/128)*common_79 - 311/64*common_81 - 411/128*common_82 - 215/64*common_83 + 44*common_84 + (7/4)*common_85 + (409/24)*common_86 + (797/16)*common_87 - common_89 + 77*common_90 - 857/8*common_91 + common_93 + (883/64)*common_94 + (187/8)*common_95 - 565/64*common_96),
        (-62191/945*common_106 - 538/45*common_107 - 223/15*common_112 - 39059/2205*common_113 - 20*common_115 - 116/9*common_117 - 1/42*common_119 - 143/70*common_120 - 5878/441*common_121 - 165716/6615*common_123 - 11798/735*common_124 + (116141/6615)*common_14*common_18*nu + (1098/245)*common_14*common_2*common_32*speed_squared + (9602/245)*common_14*common_2*nu*speed_squared + (83458/2205)*common_18*common_32*radial_speed*speed_squared + (1202/105)*common_18*radial_speed*speed_squared + (4/7)*common_2*common_42, (55477/2835)*common_56 + (2650/189)*common_57 - 323/35*common_58 - 686/45*common_59 + (10/7)*common_60 + (149/30)*common_61 + (49/5)*common_62 + (2309/2205)*common_63 + (34514/2205)*common_64 + (13/42)*common_65 - 30754/2205*common_66 - 304/135*common_70 + (533/630)*common_71 + (2113/735)*common_75 - 3343/490*common_76 + (3184/315)*common_78 - 14629/735*common_84 + (20609/1470)*common_85 + (1121/245)*common_86 + (5116/441)*common_87 - 11974/735*common_91, (4588/945)*common_125*common_32 + (2990/189)*common_125*nu + (192656/2835)*common_125 + (349/5)*common_18*common_22 + (1039/105)*common_18*common_6 + (183/14)*common_19*common_22 - 544/315*common_2*common_31 - 381/70*common_2*common_33 - 4*common_2*common_68 - 212/35*common_20*common_6 - 299884/2205*common_21*common_22 - 76574/2205*common_21*common_6 - 6219/98*common_22*common_26 + (135196/2205)*common_22*common_35 - 579/98*common_22*common_37 + (2477/105)*common_29*common_6 + (324221/945)*common_3*common_56 - 217003/2205*common_3*common_57 + (320869/2205)*common_3*common_59 + (41506/735)*common_3*common_70 - 72964/735*common_3*common_78 + (901/126)*common_31*common_9 + (81091/2205)*common_35*common_6 + (311/35)*common_38*common_6 - 72676/945*common_56*speed_squared + (2130/49)*common_57*speed_squared - 2882/35*common_62*speed_squared + 28*common_68*common_9 - 17926/735*common_70*speed_squared, -112052/945*common_106 - 1328/735*common_107 - 46/7*common_111 - 22*common_112 - 11077/735*common_113 - 5494/147*common_115 - 4076/245*common_117 - 8420/147*common_121 - 680/49*common_122 - 289004/6615*common_123 - 7964/245*common_124 + (222071/6615)*common_14*common_18*nu + 36*common_14*common_2*common_32*speed_squared + (9878/105)*common_14*common_2*nu*speed_squared + (96/35)*common_14*common_2*speed_squared + (40406/735)*common_18*common_32*radial_speed*speed_squared + (3547/105)*common_18*radial_speed*speed_squared + (64/15)*common_2*common_6*radial_speed),
    )
