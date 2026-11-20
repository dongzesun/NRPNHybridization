"""Generated source mass quadrupole jets."""

import math

GRADES = (0, 2, 4, 5, 6, 7)

def derivative_0(potential, radial_speed, speed_squared, nu, log_one, log_infrared):
    common_0 = nu*potential
    common_1 = nu*speed_squared
    common_2 = nu*radial_speed
    common_3 = potential**2
    common_4 = speed_squared**2
    common_5 = common_1*potential
    common_6 = common_3*nu
    common_7 = common_4*nu
    common_8 = radial_speed**2
    common_9 = common_8*potential
    common_10 = nu**2
    common_11 = potential*speed_squared
    common_12 = common_10*common_11
    common_13 = common_10*common_9
    common_14 = potential*radial_speed
    common_15 = radial_speed*speed_squared
    common_16 = (25/63)*common_8
    common_17 = common_10*potential
    common_18 = common_10*speed_squared
    common_19 = common_2*common_3
    common_20 = potential**3
    common_21 = speed_squared**3
    common_22 = (428/105)*log_infrared
    common_23 = radial_speed**4
    common_24 = common_20*nu
    common_25 = common_0*common_4
    common_26 = common_1*common_3
    common_27 = common_10*common_20
    common_28 = nu**3
    common_29 = common_3*common_8
    common_30 = common_10*common_3
    common_31 = common_30*common_8
    common_32 = common_11*radial_speed
    common_33 = radial_speed**3
    common_34 = common_3*radial_speed
    common_35 = common_12*radial_speed
    common_36 = common_30*radial_speed
    common_37 = common_10*common_4
    common_38 = common_8*speed_squared
    return (
        (1, 0, 0),
        ((8/7)*common_0 - 29/14*common_1 - 5/7*potential + (29/42)*speed_squared, (12/7)*common_2 - 4/7*radial_speed, 11/21 - 11/7*nu),
        ((337/252)*common_10*common_3 + (3545/504)*common_10*common_4 - 4883/756*common_12 - 1273/756*common_13 - 355/252*common_3 + (253/504)*common_4 - 5947/756*common_5 - 953/126*common_6 - 1835/504*common_7 + (907/756)*common_8*nu*potential - 131/756*common_9 + (2021/756)*potential*speed_squared, (4057/378)*common_0*radial_speed + (202/63)*common_1*radial_speed + (209/54)*common_10*common_14 - 418/63*common_10*common_15 - 155/54*common_14 - 26/63*common_15, -335/189*common_0 - 337/126*common_1 + common_10*common_16 - common_16*nu - 985/189*common_17 + (733/126)*common_18 + (5/63)*common_8 + (106/27)*potential + (41/126)*speed_squared),
        (-24/7*common_19, (48/7)*common_6, 0),
        (-1745/2772*common_0*common_23 + (117067/5544)*common_10*common_21 + (16319/5544)*common_10*common_23*potential + (434909/16632)*common_10*common_3*speed_squared + (218411/8316)*common_10*common_4*potential - 98311/16632*common_18*common_9 - common_20*common_22 + (13289/8316)*common_20*common_28 + (15502/385)*common_20*nu + (6285233/207900)*common_20 - 328663/11088*common_21*common_28 - 7993/1584*common_21*nu + (4561/11088)*common_21 - 311/99*common_23*common_28*potential + (2/99)*common_23*potential - 44/3*common_24*log_one - 94475/4158*common_25 - 605419/16632*common_26 - 3632/693*common_27 - 54121/5544*common_28*common_29 - 37369/2772*common_28*common_3*speed_squared + (299857/8316)*common_28*common_4*potential + (153407/8316)*common_28*common_8*potential*speed_squared - 8539/20790*common_29 + (52153/4158)*common_3*common_8*nu + (187183/83160)*common_3*speed_squared - 4652/231*common_31 + (307/77)*common_4*potential + (5545/8316)*common_8*nu*potential*speed_squared - 757/5544*common_9*speed_squared, -8611/2772*common_17*common_33 + (4441/1386)*common_28*common_3*radial_speed - 58565/2079*common_28*common_32 - 895/77*common_28*common_33*potential + (40687/1386)*common_28*common_4*radial_speed + (406333/8316)*common_3*nu*radial_speed - 2839/693*common_32 + (3233/2772)*common_33*nu*potential + (305/2772)*common_33*potential - 12587/20790*common_34 - 188063/4158*common_35 - 2713/198*common_36 - 13693/693*common_37*radial_speed + (6103/1386)*common_4*nu*radial_speed - 457/1386*common_4*radial_speed + (237893/8316)*nu*potential*radial_speed*speed_squared, -5807/1386*common_0*common_8 - 1135/1386*common_1*common_8 + (129781/4158)*common_11*common_28 + (587/154)*common_11 + (25660/2079)*common_12 + (515/1386)*common_13 + (1795/693)*common_18*common_8 + common_22*common_3 - 24314/2079*common_28*common_3 - 3445/1386*common_28*common_38 - 139999/5544*common_28*common_4 + (8245/693)*common_28*common_9 - 40716/1925*common_3 + (62576/2079)*common_30 + (45421/2772)*common_37 + (115/1386)*common_38 + (1369/5544)*common_4 - 67933/4158*common_5 - 10762/2079*common_6 - 19351/5544*common_7 + (79/77)*common_9),
        (-1492/63*common_1*common_34 + (27086/315)*common_15*common_30 + (1660/63)*common_2*common_20 + (32/135)*common_27*radial_speed - 4342/63*common_30*common_33 + (120/7)*common_33*common_6, -10064/315*common_24 + (16/105)*common_25 + (3616/315)*common_26 + (13936/945)*common_27 - 272/15*common_30*speed_squared - 104/7*common_31 - 64/105*common_37*potential + (4/7)*common_6*common_8, -8/21*common_1*common_14 - 664/63*common_19 + (32/21)*common_35 + (104/9)*common_36),
    )

def derivative_1(potential, radial_speed, speed_squared, nu, log_one, log_infrared):
    common_0 = potential*radial_speed
    common_1 = common_0*nu
    common_2 = nu*potential
    common_3 = nu*speed_squared
    common_4 = nu*radial_speed
    common_5 = common_0*speed_squared
    common_6 = radial_speed**3
    common_7 = common_6*potential
    common_8 = potential**2
    common_9 = common_8*radial_speed
    common_10 = common_0*common_3
    common_11 = common_2*common_6
    common_12 = common_4*common_8
    common_13 = nu**2
    common_14 = common_13*common_5
    common_15 = common_13*common_7
    common_16 = common_13*common_9
    common_17 = speed_squared**2
    common_18 = common_8*nu
    common_19 = common_17*nu
    common_20 = radial_speed**2
    common_21 = common_20*potential
    common_22 = common_13*common_21
    common_23 = radial_speed*speed_squared
    common_24 = potential**3
    common_25 = common_24*nu
    common_26 = common_3*common_8
    common_27 = radial_speed**5
    common_28 = common_27*potential
    common_29 = common_24*radial_speed
    common_30 = common_24*common_4
    common_31 = common_0*common_17
    common_32 = common_7*speed_squared
    common_33 = common_9*speed_squared
    common_34 = common_6*common_8
    common_35 = common_3*common_9
    common_36 = common_13*common_29
    common_37 = nu**3
    common_38 = common_13*common_17
    common_39 = common_13*common_8
    common_40 = common_39*common_6
    common_41 = speed_squared**3
    common_42 = radial_speed**4
    common_43 = common_42*potential
    common_44 = common_8*speed_squared
    common_45 = common_17*common_2
    common_46 = common_13*common_24
    common_47 = common_20*common_8
    common_48 = common_20*common_39
    common_49 = common_17*radial_speed
    common_50 = potential**4
    return (
        (0, 2, 0),
        (-2/21*common_0 + (9/7)*common_1, (26/7)*common_2 - 17/7*common_3 - 40/21*potential + (17/21)*speed_squared, -10/7*common_4 + (10/21)*radial_speed),
        ((577/252)*common_10 - 601/252*common_11 + (3641/378)*common_12 - 715/252*common_14 + (355/252)*common_15 + (923/189)*common_16 + (25/252)*common_5 + (131/252)*common_7 + (106/27)*common_9, (1873/252)*common_13*common_17 + (172/27)*common_13*common_8 - 205/21*common_13*potential*speed_squared + (149/252)*common_17 - 10235/378*common_18 - 1027/252*common_19 + (58/21)*common_20*nu*potential - 29/21*common_21 - 79/21*common_22 - 20/7*common_3*potential - 1114/189*common_8 + (37/21)*potential*speed_squared, -943/126*common_0*common_13 + (559/126)*common_0 + (55/126)*common_1 + (365/63)*common_13*common_23 + (25/63)*common_23 - 185/63*common_3*radial_speed),
        ((72/7)*common_20*common_8*nu - 24/7*common_25 - 24/7*common_26, -96/7*common_12, (48/7)*common_18),
        ((11965/5544)*common_0*common_19 - 7919/924*common_0*common_38 - 31105/5544*common_13*common_28 + (2953/264)*common_15*speed_squared - 31774/2079*common_16*speed_squared - 44255/2079*common_18*common_6 + (9035/5544)*common_2*common_27 + (11345/5544)*common_28*common_37 - 10/99*common_28 + (28393/2079)*common_29*common_37 + (428/35)*common_29*log_infrared - 1673842/17325*common_29 - 155/33*common_3*common_7 + 44*common_30*log_one - 822689/5940*common_30 + (19513/1848)*common_31*common_37 + (211/5544)*common_31 - 11999/1386*common_32*common_37 + (2719/5544)*common_32 - 149753/8316*common_33*common_37 + (41401/6930)*common_33 + (79819/8316)*common_34*common_37 + (20353/10395)*common_34 + (50768/2079)*common_35 + (31940/2079)*common_36 + (20311/8316)*common_40, (2641/1386)*common_13*common_17*potential + (20765/924)*common_13*common_41 + (2141/126)*common_13*common_42*potential + (137507/2772)*common_13*common_8*speed_squared + (123023/2772)*common_17*common_37*potential + (3211/924)*common_17*potential - 95/63*common_2*common_42 + (21607/693)*common_20*common_37*potential*speed_squared + (15805/594)*common_20*common_8*nu + (3578/693)*common_20*nu*potential*speed_squared - 97/154*common_21*speed_squared - 601/22*common_22*speed_squared + (3052/297)*common_24*common_37 - 1712/105*common_24*log_infrared + (3460937/41580)*common_24*nu + (883664/7425)*common_24 - 88/3*common_25*log_one - 14477/308*common_26 - 18435/616*common_37*common_41 - 2129/252*common_37*common_43 - 15145/594*common_37*common_44 - 52708/2079*common_37*common_47 - 10513/1848*common_41*nu + (911/1848)*common_41 - 73/252*common_43 - 152161/41580*common_44 - 9881/693*common_45 - 78467/2079*common_46 - 673213/41580*common_47 - 112615/2079*common_48, -13247/924*common_10 - 1075/308*common_11 + (182041/8316)*common_12 - 5359/1386*common_14 + (7349/2772)*common_15 + (99515/4158)*common_16 - 3895/924*common_17*common_4 - 8045/308*common_37*common_49 + (4228/99)*common_37*common_5 + (8797/1386)*common_37*common_7 - 102871/4158*common_37*common_9 + (8405/462)*common_38*radial_speed + (305/924)*common_49 + (496/99)*common_5 + (521/2772)*common_7 + (547727/20790)*common_9),
        ((29366/315)*common_13*common_17*common_8 + (7310/189)*common_13*common_20*common_24 + (21710/63)*common_13*common_42*common_8 - 512/63*common_13*common_50 - 1172/45*common_17*common_18 - 600/7*common_18*common_42 - 36668/315*common_20*common_25 + (13268/105)*common_20*common_8*nu*speed_squared + (2992/105)*common_24*nu*speed_squared - 54826/945*common_46*speed_squared - 10004/21*common_48*speed_squared + (676/35)*common_50*nu, (46444/315)*common_13*common_8*radial_speed*speed_squared + (47504/315)*common_24*nu*radial_speed - 16832/315*common_35 - 3952/135*common_36 - 5876/63*common_40 + (228/7)*common_6*common_8*nu, (32/35)*common_13*common_17*potential + (17272/945)*common_13*common_24 + (748/63)*common_20*common_8*nu - 2776/105*common_25 - 16/45*common_26 - 968/315*common_39*speed_squared - 8/35*common_45 - 1856/63*common_48),
    )

def derivative_2(potential, radial_speed, speed_squared, nu, log_one, log_infrared):
    common_0 = potential*speed_squared
    common_1 = (61/21)*common_0
    common_2 = potential**2
    common_3 = common_2*nu
    common_4 = common_0*nu
    common_5 = radial_speed**2
    common_6 = common_5*potential
    common_7 = common_6*nu
    common_8 = potential*radial_speed
    common_9 = nu*potential
    common_10 = nu*speed_squared
    common_11 = potential**3
    common_12 = common_11*nu
    common_13 = radial_speed**4
    common_14 = common_13*potential
    common_15 = speed_squared**2
    common_16 = common_15*potential
    common_17 = nu**2
    common_18 = common_11*common_17
    common_19 = common_2*common_5
    common_20 = common_4*common_5
    common_21 = common_2*speed_squared
    common_22 = common_17*common_21
    common_23 = common_17*common_19
    common_24 = common_0*radial_speed
    common_25 = radial_speed**3
    common_26 = common_25*potential
    common_27 = common_2*radial_speed
    common_28 = common_3*radial_speed
    common_29 = common_15*nu
    common_30 = common_0*common_17
    common_31 = common_12*radial_speed
    common_32 = common_3*speed_squared
    common_33 = common_32*radial_speed
    common_34 = common_25*common_3
    common_35 = potential**4
    common_36 = radial_speed**6
    common_37 = speed_squared**3
    common_38 = common_37*potential
    common_39 = common_11*speed_squared
    common_40 = common_35*nu
    common_41 = common_10*common_11
    common_42 = common_0*common_13
    common_43 = common_17*common_35
    common_44 = nu**3
    common_45 = math.pi**2
    common_46 = common_11*common_5
    common_47 = common_13*common_3
    common_48 = common_15*common_3
    common_49 = common_19*speed_squared
    common_50 = common_15*common_17
    common_51 = common_2*common_50
    common_52 = common_12*common_5
    common_53 = radial_speed**5
    common_54 = common_53*potential
    common_55 = common_11*radial_speed
    common_56 = common_15*common_8
    common_57 = common_0*common_25
    common_58 = common_21*radial_speed
    common_59 = common_2*common_25
    common_60 = common_18*radial_speed
    common_61 = common_17*common_2
    common_62 = common_25*common_61
    common_63 = common_3*common_5
    common_64 = common_18*speed_squared
    return (
        (-2*potential, 0, 2),
        (-common_1 + 10*common_2 - common_3 - 16/7*common_4 + (2/7)*common_5*potential - 6/7*common_7, (18/7)*common_8*nu + (50/7)*common_8, -27/7*common_10 + (36/7)*common_9 - 50/21*potential + (9/7)*speed_squared),
        (-common_1*common_17*common_5 - 1514/63*common_11 - 181/18*common_12 + (265/63)*common_13*common_17*potential + (515/63)*common_13*nu*potential - 655/252*common_14 + (316/63)*common_15*common_17*potential + (23/63)*common_15*nu*potential - 82/63*common_16 - 79/63*common_18 - 9421/756*common_19 + (10295/756)*common_2*common_5*nu + (12175/756)*common_2*nu*speed_squared + (5521/756)*common_2*speed_squared - 80/21*common_20 - 4819/756*common_22 - 4079/756*common_23 + (53/42)*common_5*potential*speed_squared, -383/126*common_17*common_24 - 437/126*common_17*common_26 + (221/27)*common_17*common_27 - 43/126*common_24 - 2293/126*common_25*common_9 + (653/126)*common_26 + (1639/189)*common_27 + (4666/189)*common_28 + (2033/126)*common_4*radial_speed, (1111/84)*common_15*common_17 + (83/84)*common_15 + (4157/378)*common_17*common_2 - 181/126*common_17*common_6 - 455/54*common_2 - 589/84*common_29 - 6100/189*common_3 - 2363/126*common_30 - 113/126*common_6 - 1397/126*common_7 + (65/126)*nu*potential*speed_squared + (671/126)*potential*speed_squared),
        ((3344/105)*common_31 + (1416/35)*common_33 - 360/7*common_34, -576/35*common_12 + (432/7)*common_2*common_5*nu - 832/35*common_32, -144/7*common_28),
        ((96674/2079)*common_11*common_17*common_5 + (83035/1386)*common_11*common_17*speed_squared + (123/32)*common_11*common_45*nu*speed_squared + (75247/231)*common_11*common_5*nu + (3000119/5940)*common_11*common_5 + (428/35)*common_11*log_infrared*speed_squared + (29817/616)*common_13*common_17*common_2 + (143821/5544)*common_13*common_2*common_44 - 19736/3465*common_13*common_2 - 730/99*common_13*common_30 + (10475/693)*common_13*nu*potential*speed_squared + (408689/16632)*common_15*common_2*common_44 + (82189/27720)*common_15*common_2 + (183/14)*common_15*common_44*common_5*potential + (19/14)*common_15*common_5*potential + (11173/2079)*common_17*common_2*common_5*speed_squared + (3445/198)*common_17*common_36*potential + (12161/1386)*common_17*common_37*potential + (1909007/16632)*common_2*common_5*nu*speed_squared - 29/28*common_29*common_6 - 2963/2079*common_35*common_44 + (428/105)*common_35*log_infrared + (1375951/6930)*common_35*nu + (1025428/51975)*common_35 + (1495/198)*common_36*common_44*potential - 2785/396*common_36*common_9 + (70/99)*common_36*potential + (3077/1386)*common_37*nu*potential - 2699/154*common_38*common_44 - 725/693*common_38 - 16393/1188*common_39*common_44 - 20125129/207900*common_39 - 41/8*common_40*common_45 - 44/3*common_40*log_one - 49355/1386*common_41 - 33115/1386*common_42*common_44 - 5465/1848*common_42 - 1615/308*common_43 - 167663/8316*common_44*common_46 - 85601/16632*common_44*common_49 - 615/32*common_45*common_52 - 428/7*common_46*log_infrared - 50377/924*common_47 - 163883/16632*common_48 - 459707/27720*common_49 - 919/42*common_50*common_6 - 34621/1188*common_51, -269585/2772*common_17*common_54 + (83009/2079)*common_22*radial_speed + (470933/2772)*common_25*common_30 - 34327/693*common_25*common_4 + (3101/132)*common_29*common_8 + (123/16)*common_31*common_45 + 88*common_31*log_one - 89953/165*common_31 + (66427/594)*common_33 - 10986/77*common_34 - 1805/396*common_44*common_54 + (41654/2079)*common_44*common_55 + (42109/2772)*common_44*common_56 + (1459/154)*common_44*common_57 - 90511/4158*common_44*common_58 - 1823/77*common_44*common_59 - 5587/66*common_50*common_8 + (47755/2772)*common_53*common_9 + (3455/2772)*common_54 + (2568/35)*common_55*log_infrared - 29933986/51975*common_55 + (1081/2772)*common_56 + (4745/2772)*common_57 - 57136/3465*common_58 + (136967/2310)*common_59 + (140146/2079)*common_60 - 53353/462*common_62, (27787/1386)*common_11*common_44 - 1712/105*common_11*log_infrared + (516037/6930)*common_11*nu + (3677287/34650)*common_11 - 88/3*common_12*log_one + (57065/2772)*common_13*common_17*potential + (3925/1386)*common_13*nu*potential - 3664/693*common_14*common_44 - 169/198*common_14 + (20816/231)*common_15*common_44*potential - 16675/693*common_15*common_9 + (10753/1386)*common_15*potential - 15173/1386*common_16*common_17 + (753793/8316)*common_17*common_2*speed_squared + (12525/308)*common_17*common_37 + (18503/924)*common_17*common_5*potential*speed_squared - 41350/693*common_18 - 25469/2772*common_19*common_44 - 172961/6930*common_19 + (205501/13860)*common_2*speed_squared - 437/22*common_20 - 233045/4158*common_21*common_44 - 77647/693*common_23 - 269267/8316*common_32 - 34525/616*common_37*common_44 - 6101/616*common_37*nu + (507/616)*common_37 + (19394/693)*common_44*common_5*potential*speed_squared + (2069/2772)*common_5*potential*speed_squared - 12163/462*common_63),
        (-32192/63*common_10*common_55 + (25224/35)*common_12*common_25 + (10984/35)*common_15*common_28 - 238136/315*common_18*common_25 + (238412/63)*common_22*common_25 - 2240/3*common_25*common_32 - 26630/21*common_27*common_50 + 376*common_3*common_53 - 14188/45*common_40*radial_speed + (70072/945)*common_43*radial_speed - 21710/9*common_53*common_61 + (268384/315)*common_64*radial_speed, (10400/9)*common_13*common_61 + (105928/315)*common_18*common_5 - 180256/105*common_23*speed_squared + (3200/63)*common_40 + (1664/7)*common_41 + (32/135)*common_43 - 660/7*common_47 - 7684/105*common_48 + (109784/315)*common_51 - 92648/105*common_52 + (8816/35)*common_63*speed_squared - 99784/315*common_64, (73384/315)*common_11*nu*radial_speed + (5540/63)*common_17*common_2*radial_speed*speed_squared - 2984/105*common_33 - 64/21*common_34 - 144/7*common_60 - 44/9*common_62),
    )

def derivative_3(potential, radial_speed, speed_squared, nu, log_one, log_infrared):
    common_0 = potential*radial_speed
    common_1 = common_0*speed_squared
    common_2 = radial_speed**3
    common_3 = common_2*potential
    common_4 = potential**2
    common_5 = common_4*radial_speed
    common_6 = common_1*nu
    common_7 = common_3*nu
    common_8 = common_5*nu
    common_9 = potential*speed_squared
    common_10 = common_9*nu
    common_11 = common_4*nu
    common_12 = radial_speed**2
    common_13 = common_12*potential
    common_14 = common_13*nu
    common_15 = common_0*nu
    common_16 = radial_speed**5
    common_17 = common_16*potential
    common_18 = potential**3
    common_19 = common_18*radial_speed
    common_20 = common_17*nu
    common_21 = common_19*nu
    common_22 = speed_squared**2
    common_23 = common_0*common_22
    common_24 = common_3*speed_squared
    common_25 = common_5*speed_squared
    common_26 = common_2*common_4
    common_27 = common_15*common_22
    common_28 = common_7*speed_squared
    common_29 = common_8*speed_squared
    common_30 = common_11*common_2
    common_31 = nu**2
    common_32 = common_17*common_31
    common_33 = common_19*common_31
    common_34 = common_23*common_31
    common_35 = common_24*common_31
    common_36 = common_25*common_31
    common_37 = common_26*common_31
    common_38 = common_18*nu
    common_39 = radial_speed**4
    common_40 = common_18*common_31
    common_41 = common_12*common_4
    common_42 = common_31*common_4*speed_squared
    common_43 = common_31*common_41
    common_44 = common_3*common_31
    common_45 = common_31*common_5
    common_46 = potential**4
    common_47 = common_46*nu
    common_48 = common_12*common_38
    common_49 = common_11*common_12
    common_50 = common_49*speed_squared
    common_51 = potential*radial_speed**7
    common_52 = common_46*radial_speed
    common_53 = common_47*radial_speed
    common_54 = speed_squared**3
    common_55 = common_0*common_54
    common_56 = common_16*common_9
    common_57 = common_19*speed_squared
    common_58 = common_16*common_4
    common_59 = common_18*common_2
    common_60 = common_21*speed_squared
    common_61 = common_11*common_16
    common_62 = common_2*common_38
    common_63 = common_31*common_52
    common_64 = nu**3
    common_65 = common_22*common_3
    common_66 = common_22*common_5
    common_67 = common_26*speed_squared
    common_68 = math.pi**2
    common_69 = common_22*common_8
    common_70 = common_30*speed_squared
    common_71 = common_33*speed_squared
    common_72 = common_31*common_58
    common_73 = common_2*common_40
    common_74 = common_22*common_45
    common_75 = common_37*speed_squared
    common_76 = radial_speed**6
    common_77 = common_76*potential
    common_78 = common_54*potential
    common_79 = common_18*speed_squared
    common_80 = common_38*speed_squared
    common_81 = common_39*common_9
    common_82 = common_31*common_46
    common_83 = common_39*common_4
    common_84 = common_22*common_4
    common_85 = common_12*common_18
    common_86 = common_11*common_22
    common_87 = common_43*speed_squared
    common_88 = potential**5
    common_89 = common_12*common_40
    return (
        (6*common_0, -8*potential, 0),
        ((65/7)*common_1 - 10/7*common_3 - 880/21*common_5 + (36/7)*common_6 + (30/7)*common_7 + (54/7)*common_8, -44/7*common_10 - 48/7*common_11 - 146/7*common_13 - 24/7*common_14 + (706/21)*common_4 - 110/21*common_9, -20/7*common_15 + (482/21)*potential*radial_speed),
        ((655/36)*common_17 + (24461/189)*common_19 - 515/9*common_20 - 16973/378*common_21 + (45/7)*common_23 - 2105/126*common_24 - 4804/63*common_25 + (5024/63)*common_26 - 61/7*common_27 + (3260/63)*common_28 - 3517/63*common_29 - 5905/63*common_30 - 265/9*common_32 + (2683/189)*common_33 - 146/7*common_34 + (1975/63)*common_35 + (214/9)*common_36 + (1429/63)*common_37, (617/9)*common_12*common_4*nu + (401/21)*common_12*potential*speed_squared - 128/3*common_13*common_31*speed_squared - 1864/21*common_14*speed_squared - 12542/189*common_18 + (1240/63)*common_22*common_31*potential + (170/9)*common_22*nu*potential - 472/63*common_22*potential + (3040/63)*common_31*common_39*potential - 2090/189*common_38 + (6290/63)*common_39*nu*potential - 280/9*common_39*potential + (1261/63)*common_4*nu*speed_squared + (1315/63)*common_4*speed_squared - 1318/189*common_40 - 475/63*common_41 - 1639/63*common_42 - 452/63*common_43, (1382/63)*common_1*common_31 + (415/63)*common_1 + (496/63)*common_3 - 1432/63*common_44 - 178/7*common_45 - 2/21*common_5 - 1292/63*common_6 + (58/63)*common_7 + (160/3)*common_8),
        ((1592/105)*common_18*nu*speed_squared + (1416/35)*common_22*common_4*nu + 360*common_39*common_4*nu - 1616/105*common_47 - 1256/7*common_48 - 2496/7*common_50, (2752/21)*common_21 + (2064/7)*common_29 - 2880/7*common_30, -1776/35*common_11*speed_squared + (864/7)*common_12*common_4*nu - 528/35*common_38),
        (-29305/198*common_10*common_16 - 2017/231*common_15*common_54 + (110755/1386)*common_22*common_44 + (181955/2772)*common_22*common_7 - 3445/22*common_31*common_51 - 10793/154*common_31*common_55 + (15445/99)*common_31*common_56 - 1495/22*common_51*common_64 + (2785/44)*common_51*nu - 70/11*common_51 + (283/11)*common_52*common_64 - 1311956/3465*common_52 + (861/16)*common_53*common_68 - 90667/99*common_53 + (12123/154)*common_55*common_64 + (1352/231)*common_55 + (42085/198)*common_56*common_64 + (2195/88)*common_56 + (16265/198)*common_57*common_64 - 1284/7*common_57*log_infrared + (11114123/6930)*common_57 - 94475/396*common_58*common_64 + (15263/3465)*common_58 + (175481/2079)*common_59*common_64 + 428*common_59*log_infrared - 38542237/10395*common_59 - 1845/32*common_60*common_68 + (2002681/2772)*common_60 + (44068/77)*common_61 + (4305/32)*common_62*common_68 - 16121773/8316*common_62 - 42368/693*common_63 - 223045/1386*common_64*common_65 - 98774/693*common_64*common_66 + (50461/231)*common_64*common_67 - 4300/231*common_65 - 5906/99*common_66 + (24091/231)*common_67 + (46268/231)*common_69 - 73068/77*common_70 - 156757/462*common_71 - 952799/2772*common_72 - 236447/594*common_73 + (294017/1386)*common_74 + (35767/462)*common_75, (301993/2079)*common_12*common_18*common_31 + (11341331/4158)*common_12*common_18*nu + (5883619/1485)*common_12*common_18 + (5832/11)*common_12*common_22*common_31*potential + (42725/231)*common_12*common_22*common_64*potential + (1543/231)*common_12*common_22*potential + (263981/693)*common_12*common_4*nu*speed_squared + (131041/693)*common_12*common_4*speed_squared - 40798/231*common_14*common_22 + (362080/2079)*common_18*common_31*speed_squared + (369/16)*common_18*common_68*nu*speed_squared + (3424/35)*common_18*log_infrared*speed_squared + (26795/693)*common_22*common_31*common_4 + (52375/462)*common_22*common_4*common_64 + (387941/396)*common_31*common_39*common_4 + (7390/11)*common_31*common_76*potential - 34918/693*common_31*common_78 - 822265/693*common_31*common_81 + (266081/693)*common_39*common_4*common_64 + (52205/924)*common_39*common_4*nu + (65185/198)*common_39*nu*potential*speed_squared - 194125/693*common_41*common_64*speed_squared - 779/297*common_46*common_64 - 3424/105*common_46*log_infrared + (2049043/2310)*common_46*nu + (3066433/7425)*common_46 - 451/16*common_47*common_68 - 176/3*common_47*log_one - 1845/16*common_48*common_68 + (51451/1386)*common_54*nu*potential + (2995/33)*common_64*common_76*potential - 48497/693*common_64*common_78 - 156026/2079*common_64*common_79 - 28595/99*common_64*common_81 - 1259/2079*common_64*common_85 - 1385/11*common_77*nu - 965/132*common_77 - 410/77*common_78 - 39315118/51975*common_79 - 732404/2079*common_80 - 7615/924*common_81 - 38260/693*common_82 - 809917/2772*common_83 - 40015/1386*common_84 - 3424/7*common_85*log_infrared - 33896/693*common_86 - 195571/231*common_87, -1205/33*common_17*common_64 + (5095/924)*common_17 - 223561/2079*common_19*common_64 + (856/7)*common_19*log_infrared - 1366126/1485*common_19 + (3195/154)*common_20 + (369/16)*common_21*common_68 - 236047/297*common_21 - 59467/693*common_23*common_64 + (689/77)*common_23 + (86138/693)*common_24*common_64 - 607/154*common_24 + (117968/693)*common_25*common_64 - 32372/693*common_25 - 1702/9*common_26*common_64 + (9566/63)*common_26 - 46213/1386*common_27 - 22355/693*common_28 + (57017/693)*common_29 - 6614/21*common_30 - 53420/231*common_32 + (371024/2079)*common_33 + (709/99)*common_34 + (183607/693)*common_35 + (730/77)*common_36 - 5885/63*common_37),
        (-3384*common_11*common_76 + (96992/15)*common_12*common_18*nu*speed_squared + (123854/7)*common_12*common_22*common_31*common_4 + (644044/315)*common_12*common_46*nu - 359824/945*common_12*common_82 + (77798/45)*common_18*common_22*common_31 + (203194/21)*common_18*common_31*common_39 - 34156/45*common_22*common_38 - 26664/7*common_22*common_49 - 26630/21*common_31*common_4*common_54 + 21710*common_31*common_4*common_76 - 13736/315*common_31*common_88 - 700876/105*common_38*common_39 + (21320/3)*common_39*common_4*nu*speed_squared - 115654/3*common_39*common_42 + (10984/35)*common_4*common_54*nu + (3908/105)*common_46*nu*speed_squared - 427088/945*common_82*speed_squared + (11412/35)*common_88*nu - 84992/7*common_89*speed_squared, -153448/315*common_53 - 1347872/315*common_60 + 964*common_61 + (149528/21)*common_62 + (185536/945)*common_63 + (46252/35)*common_69 - 56152/21*common_70 + (1802488/315)*common_71 - 38740/3*common_72 - 111880/21*common_73 - 248508/35*common_74 + (437272/21)*common_75, (8420/21)*common_11*common_39 + 1180*common_31*common_83 + (16476/35)*common_31*common_84 - 133676/315*common_40*speed_squared - 4168/105*common_47 - 586672/315*common_48 - 7072/35*common_50 + (29696/63)*common_80 + (12752/945)*common_82 - 268/7*common_86 - 71464/35*common_87 + (159428/315)*common_89),
    )

def derivative_4(potential, radial_speed, speed_squared, nu, log_one, log_infrared):
    common_0 = potential*speed_squared
    common_1 = potential**2
    common_2 = radial_speed**2
    common_3 = common_2*potential
    common_4 = potential*radial_speed
    common_5 = potential**3
    common_6 = common_5*nu
    common_7 = radial_speed**4
    common_8 = speed_squared**2
    common_9 = common_1*speed_squared
    common_10 = common_7*potential
    common_11 = common_10*nu
    common_12 = common_3*speed_squared
    common_13 = common_12*nu
    common_14 = common_1*common_2
    common_15 = common_14*nu
    common_16 = common_0*radial_speed
    common_17 = radial_speed**3
    common_18 = common_17*potential
    common_19 = common_1*radial_speed
    common_20 = common_19*nu
    common_21 = common_0*nu
    common_22 = common_1*nu
    common_23 = potential**4
    common_24 = radial_speed**6
    common_25 = speed_squared**3
    common_26 = common_25*potential
    common_27 = nu**2
    common_28 = common_1*common_7
    common_29 = common_2*common_5
    common_30 = radial_speed**5
    common_31 = common_30*potential
    common_32 = common_5*radial_speed
    common_33 = common_6*radial_speed
    common_34 = common_4*common_8
    common_35 = common_0*common_17
    common_36 = common_9*radial_speed
    common_37 = common_1*common_17
    common_38 = common_9*nu
    common_39 = common_38*radial_speed
    common_40 = common_17*common_22
    common_41 = common_8*potential
    common_42 = common_23*nu
    common_43 = common_6*speed_squared
    return (
        (6*common_0 + 2*common_1 - 30*common_3, 36*common_4, -8*potential),
        (282*common_1*common_2 + (104/7)*common_1*nu*speed_squared - 30*common_11 - 355/7*common_12 - 90/7*common_13 - 81*common_15 + (2/7)*common_5 - 34/7*common_6 + 10*common_7*potential + (36/7)*common_8*nu*potential + (65/7)*common_8*potential - 923/21*common_9, (156/7)*common_16*nu - 52/7*common_16 + (180/7)*common_18*nu + (710/7)*common_18 - 244*common_19 + 84*common_20, (32/3)*common_1 + (36/7)*common_2*nu*potential - 64/7*common_21 - 4*common_22 - 628/7*common_3 + (124/7)*potential*speed_squared),
        (-1100/3*common_0*common_27*common_7 + (983/14)*common_1*common_2*nu*speed_squared + (11481/14)*common_1*common_2*speed_squared + (428/21)*common_1*common_27*common_8 + (72199/84)*common_1*common_7*nu - 1664/21*common_1*common_8 - 772/7*common_14*common_27*speed_squared + (595/3)*common_2*common_27*common_8*potential + (239483/378)*common_2*common_5*nu + (4175/21)*common_2*common_8*nu*potential - 1945/3*common_21*common_7 - 71*common_22*common_8 - 347/63*common_23*common_27 + (5438/63)*common_23*nu - 4927/63*common_23 + 265*common_24*common_27*potential + 515*common_24*nu*potential - 655/4*common_24*potential + (45/7)*common_25*potential - 146/7*common_26*common_27 - 61/7*common_26*nu - 2047/12*common_27*common_28 - 39787/189*common_27*common_29 + (2944/189)*common_27*common_5*speed_squared - 57251/84*common_28 - 225866/189*common_29 - 3455/42*common_3*common_8 + (2053/378)*common_5*nu*speed_squared + (39554/189)*common_5*speed_squared + (2495/12)*common_7*potential*speed_squared, (2840/3)*common_17*common_21 - 1190/3*common_27*common_31 + (4750/27)*common_27*common_32 - 3908/21*common_27*common_34 + (9850/21)*common_27*common_35 + (1380/7)*common_27*common_36 + (694/21)*common_27*common_37 - 2440/3*common_31*nu + (1525/6)*common_31 + (130316/189)*common_32 - 68443/189*common_33 - 5284/21*common_34*nu + (1544/21)*common_34 - 760/3*common_35 - 6910/21*common_36 + (524/21)*common_37 + (596/21)*common_39 - 14950/21*common_40, (3400/21)*common_10*common_27 - 1480/21*common_10 + (2000/21)*common_11 - 530/3*common_12*common_27 + (482/21)*common_12 - 514/21*common_13 + (922/7)*common_14*common_27 + (4861/21)*common_14 - 5126/21*common_15 + (874/21)*common_27*common_41 + (344/27)*common_27*common_5 - 1361/21*common_27*common_9 + (195/7)*common_38 - 34/21*common_41*nu - 19/21*common_41 + (4828/189)*common_5 - 5654/189*common_6 - 184/21*common_9),
        (3936*common_17*common_38 + 1120*common_17*common_6 - 6408/7*common_20*common_8 - 3240*common_22*common_30 + (1808/7)*common_42*radial_speed - 1488/7*common_43*radial_speed, 3600*common_1*common_7*nu + (13152/35)*common_1*common_8*nu - 23952/7*common_15*speed_squared - 5232/7*common_2*common_6 - 3264/35*common_42 - 96/5*common_43, (4832/105)*common_33 + (24288/35)*common_39 - 7200/7*common_40),
    )

def derivative_5(potential, radial_speed, speed_squared, nu, log_one, log_infrared):
    common_0 = potential*speed_squared
    common_1 = common_0*radial_speed
    common_2 = radial_speed**3
    common_3 = common_2*potential
    common_4 = potential**2
    common_5 = radial_speed**2
    common_6 = 240*potential
    common_7 = potential*radial_speed
    common_8 = radial_speed**5
    common_9 = common_8*potential
    common_10 = potential**3
    common_11 = common_10*radial_speed
    common_12 = speed_squared**2
    common_13 = common_12*common_7
    common_14 = common_2*common_4
    common_15 = common_13*nu
    common_16 = common_3*speed_squared
    common_17 = common_16*nu
    common_18 = common_4*speed_squared
    common_19 = common_18*radial_speed
    common_20 = common_19*nu
    common_21 = common_10*nu
    common_22 = radial_speed**4
    common_23 = common_22*nu
    common_24 = 60*nu
    common_25 = common_4*common_5
    common_26 = common_25*nu
    common_27 = common_4*radial_speed
    common_28 = common_27*nu
    common_29 = potential*radial_speed**7
    common_30 = potential**4
    common_31 = common_30*radial_speed
    common_32 = common_31*nu
    common_33 = speed_squared**3
    common_34 = common_33*common_7
    common_35 = common_0*common_8
    common_36 = common_11*speed_squared
    common_37 = common_4*common_8
    common_38 = common_10*common_2
    common_39 = common_11*nu
    common_40 = common_39*speed_squared
    common_41 = common_37*nu
    common_42 = common_2*common_21
    common_43 = nu**2
    common_44 = common_12*common_3
    common_45 = common_12*common_27
    common_46 = common_14*speed_squared
    common_47 = common_12*common_28
    common_48 = common_14*nu
    common_49 = common_48*speed_squared
    common_50 = radial_speed**6
    common_51 = common_33*potential
    common_52 = common_12*common_4
    common_53 = common_10*common_5
    common_54 = common_26*speed_squared
    common_55 = common_30*nu
    return (
        (-90*common_1 + 210*common_3, -16*common_4 - common_5*common_6 + 48*potential*speed_squared, 60*common_7),
        ((408/7)*common_10*nu*radial_speed - 5924/21*common_11 - 1035/7*common_13 - 2536*common_14 - 360/7*common_15 - 30*common_17 + 870*common_2*common_4*nu + 395*common_2*potential*speed_squared - 1794/7*common_20 + (6716/7)*common_4*radial_speed*speed_squared + 270*common_8*nu*potential - 90*common_9, -common_0*common_24*common_5 + (6368/21)*common_10 + (228/7)*common_12*nu*potential + (78/7)*common_12*potential - 7978/21*common_18 - 320/7*common_21 - 690*common_22*potential - common_23*common_6 - 912*common_26 + 2206*common_4*common_5 + (348/7)*common_4*nu*speed_squared + 240*common_5*potential*speed_squared, common_1*common_24 - 620/3*common_27 + 140*common_28 + 550*common_3 - common_6*radial_speed*speed_squared),
        (-2915*common_29*common_43 - 5665*common_29*nu + (7205/4)*common_29 + (30512/189)*common_31*common_43 + (556432/189)*common_31 - 203048/189*common_32 + (10520/21)*common_34*common_43 + (9265/21)*common_34*nu - 590/3*common_34 + 4890*common_35*common_43 + 8925*common_35*nu - 11415/4*common_35 - 4432/9*common_36*common_43 - 361859/63*common_36 + (4255/3)*common_37*common_43 + (153050/21)*common_37 + (62833/21)*common_38*common_43 + (42986/3)*common_38 + (74603/126)*common_40 - 205435/21*common_41 - 397589/42*common_42 - 2855*common_43*common_44 - 3814/7*common_43*common_45 + (16985/21)*common_43*common_46 - 3985*common_44*nu + (2815/2)*common_44 + (16718/7)*common_45 - 218138/21*common_46 + (5633/7)*common_47 + (68506/21)*common_49, -6000*common_0*common_22*common_43 - 11990*common_0*common_23 - 1084/9*common_10*common_43*speed_squared + (469676/63)*common_10*common_5*nu + (22594/63)*common_10*nu*speed_squared + (112726/63)*common_10*speed_squared + (8584/21)*common_12*common_4*common_43 + (19140/7)*common_12*common_43*common_5*potential + (31470/7)*common_12*common_5*nu*potential - 9045/7*common_12*common_5*potential + (14/3)*common_22*common_4*common_43 + (184813/21)*common_22*common_4*nu - 11825/21*common_22*common_4 + 3460*common_22*potential*speed_squared - 2158*common_25*common_43*speed_squared - 11420/189*common_30*common_43 + (3452/189)*common_30*nu - 371890/189*common_30 + (1814/21)*common_33*potential + (27163/7)*common_4*common_5*speed_squared + 4100*common_43*common_50*potential - 4784/21*common_43*common_51 - 164768/63*common_43*common_53 + 8350*common_50*nu*potential - 2615*common_50*potential - 5650/21*common_51*nu - 142/21*common_52*nu - 12332/21*common_52 - 548878/63*common_53 - 23269/7*common_54, (2672/63)*common_11*common_43 - 21566/63*common_11 - 4650/7*common_13*common_43 + (855/7)*common_13 - 3170/3*common_14*common_43 - 7690/3*common_14 - 2070/7*common_15 + 2000*common_16*common_43 - 650*common_16 + 1450*common_17 + (20408/21)*common_19*common_43 + (11636/21)*common_19 - 10160/21*common_20 + (9655/63)*common_39 - 1530*common_43*common_9 + (3230/3)*common_48 - 1480*common_9*nu + (1495/2)*common_9),
        ((10776/35)*common_10*common_12*nu + (1872/5)*common_10*common_5*nu*speed_squared + 18216*common_12*common_4*common_5*nu - 51624*common_18*common_23 - 6472*common_21*common_22 + (15136/35)*common_30*nu*speed_squared - 6408/7*common_33*common_4*nu + 35640*common_4*common_50*nu - 21664/5*common_5*common_55 - 5776/35*nu*potential**5, (256864/105)*common_32 + (74496/35)*common_40 - 38880*common_41 + 2720*common_42 - 73872/7*common_47 + 46224*common_49, (4720/7)*common_10*common_5*nu + (7488/7)*common_12*common_4*nu - 13472/21*common_21*speed_squared + 10800*common_22*common_4*nu - 69840/7*common_54 - 1312/21*common_55),
    )
