"""Generated SH gauge moments and leading multipole jets."""

import math

INVENTORY = {'w': (1, 0, 0, 7, 2), 'w1': (2, 1, 0, 7, 0), 'w2': (3, 2, 0, 4, 0), 'x': (2, 0, 0, 3, 0), 'y1': (1, 1, 0, 6, 0), 'y2': (2, 2, 0, 6, 0), 'z1': (2, 0, 1, 6, 0), 'i2': (2, 2, 0, 8, 2), 'i3': (3, 3, 0, 6, 2), 'j2': (2, 1, 1, 6, 2)}

def w(potential, radial_speed, speed_squared, nu, log_one, log_infrared):
    common_0 = (1/3)*radial_speed
    common_1 = potential*radial_speed
    common_2 = radial_speed*speed_squared
    common_3 = common_1*nu
    common_4 = (9/10)*nu
    common_5 = (1/3)*potential
    common_6 = potential**2
    common_7 = speed_squared**2
    common_8 = common_6*nu
    common_9 = radial_speed**2
    common_10 = common_9*potential
    common_11 = radial_speed**3
    common_12 = common_11*potential
    common_13 = common_8*radial_speed
    common_14 = potential**3
    common_15 = common_14*nu
    common_16 = radial_speed**4
    common_17 = common_16*potential
    common_18 = common_7*potential
    common_19 = common_6*speed_squared
    common_20 = common_17*nu
    common_21 = common_8*speed_squared
    common_22 = common_10*speed_squared
    common_23 = common_6*common_9
    common_24 = common_8*common_9
    common_25 = radial_speed**5
    common_26 = common_25*potential
    common_27 = potential*speed_squared
    common_28 = potential**4
    common_29 = common_28*nu
    common_30 = radial_speed**6
    common_31 = common_30*potential
    common_32 = speed_squared**3
    common_33 = common_32*potential
    common_34 = common_14*speed_squared
    common_35 = common_16*common_27
    common_36 = common_16*common_6
    common_37 = common_6*common_7
    common_38 = common_14*common_9
    common_39 = common_27*nu
    common_40 = common_15*common_9
    common_41 = common_10*common_7
    common_42 = common_23*speed_squared
    common_43 = radial_speed**7
    common_44 = potential**5
    common_45 = potential*radial_speed**8
    common_46 = potential*speed_squared**4
    common_47 = common_10*common_32
    return (
        (common_0,),
        ((7/5)*common_1 - common_2*common_4 + (3/10)*common_2 + (17/15)*common_3,),
        (-common_5 + (1/3)*speed_squared,),
        (-2/3*common_10 - common_4*common_7 - 1/15*common_6 + (3/10)*common_7 - 7/15*common_8 + (1/2)*common_9*nu*potential + (31/30)*nu*potential*speed_squared + (23/30)*potential*speed_squared,),
        (-common_0*potential,),
        (2*common_11*potential - 1/2*common_12*nu - 4/5*common_13 - 13/10*common_2*potential + (13/5)*common_6*radial_speed + (7/30)*nu*potential*radial_speed*speed_squared,),
        (common_10 - common_5*speed_squared + (1/3)*common_6,),
        (-59/15*common_14 + (2/15)*common_15 - 10*common_17 + (7/30)*common_18*nu - 13/10*common_18 + (127/30)*common_19 + (5/2)*common_20 - 1/30*common_21 - 11/5*common_22*nu + (99/10)*common_22 - 227/15*common_23 + (22/5)*common_24,),
        (-5*common_12 - 8/3*common_6*radial_speed + 3*potential*radial_speed*speed_squared,),
        (-179/2*common_11*common_27 + 119*common_11*common_6 - 34*common_11*common_8 + 21*common_11*nu*potential*speed_squared + (234/5)*common_14*radial_speed - 101/15*common_15*radial_speed - 329/5*common_2*common_6 + 70*common_25*potential - 35/2*common_26*nu - 51/10*common_3*common_7 + (146/15)*common_6*nu*radial_speed*speed_squared + (237/10)*common_7*potential*radial_speed,),
        ((8/3)*common_14 + 35*common_17 + 3*common_18 - 17/3*common_19 - 30*common_22 + 25*common_23,),
        (-37/15*common_15*speed_squared - 469/2*common_16*common_39 + 334*common_16*common_8 - 379/2*common_21*common_9 - 862/15*common_28 + (7/5)*common_29 + (315/2)*common_31*nu - 630*common_31 - 51/10*common_33*nu + (237/10)*common_33 + (1909/15)*common_34 + (1953/2)*common_35 - 1183*common_36 - 185/2*common_37 - 1799/3*common_38 + 113*common_40 + (177/2)*common_41*nu - 387*common_41 + (1941/2)*common_42 + (35/6)*common_7*common_8,),
        (-75*common_1*common_7 - 280*common_11*common_6 + 350*common_11*potential*speed_squared - 172/3*common_14*radial_speed - 315*common_26 + 132*common_6*radial_speed*speed_squared,),
        (-1785/2*common_1*common_32 + 8708*common_11*common_14 - 1856*common_11*common_15 - 15414*common_11*common_19 + 3364*common_11*common_6*nu*speed_squared + 6615*common_11*common_7*potential - 3115/2*common_12*common_7*nu - 2202/5*common_13*common_7 - 19686/5*common_14*common_2 + (2048/5)*common_14*nu*radial_speed*speed_squared - 25137/2*common_25*common_27 + 14217*common_25*common_6 - 3976*common_25*common_8 + (6111/2)*common_25*nu*potential*speed_squared + (7796/5)*common_28*radial_speed - 2324/15*common_29*radial_speed + (405/2)*common_32*nu*potential*radial_speed - 3465/2*common_43*nu*potential + 6930*common_43*potential + (16119/5)*common_6*common_7*radial_speed,),
        ((172/3)*common_28 + 3465*common_31 - 75*common_33 - 568/3*common_34 - 4725*common_35 + 3675*common_36 + 207*common_37 + 1092*common_38 + 1575*common_41 - 2646*common_42,),
        (-142905*common_14*common_16 - 7593*common_14*common_7 + 32634*common_15*common_16 + 304*common_15*common_7 + (538965/2)*common_16*common_19 - 124635/2*common_16*common_21 - 244755/2*common_17*common_7 + 29295*common_20*common_7 - 443862/5*common_23*common_7 + (77301/5)*common_24*common_7 + (373527/2)*common_27*common_30 - 545624/15*common_28*common_9 + (91226/15)*common_28*speed_squared + (25874/5)*common_29*common_9 - 1928/15*common_29*speed_squared - 45738*common_30*common_39 - 200277*common_30*common_6 + 55566*common_30*common_8 + (41913/10)*common_32*common_6 - 4179/10*common_32*common_8 + (500328/5)*common_38*speed_squared - 74754/5*common_40*speed_squared + (604/15)*common_44*nu - 26828/15*common_44 + (45045/2)*common_45*nu - 90090*common_45 + (405/2)*common_46*nu - 1785/2*common_46 - 6090*common_47*nu + (52185/2)*common_47,),
    )

def w1(potential, radial_speed, speed_squared, nu, log_one, log_infrared):
    common_0 = potential*radial_speed
    common_1 = potential**2
    common_2 = radial_speed**2
    common_3 = common_2*potential
    common_4 = common_0*speed_squared
    common_5 = radial_speed**3
    common_6 = common_5*potential
    common_7 = common_1*radial_speed
    common_8 = potential**3
    common_9 = radial_speed**4
    common_10 = common_9*potential
    common_11 = speed_squared**2
    common_12 = common_11*potential
    common_13 = common_1*common_2
    common_14 = potential*radial_speed**5
    common_15 = common_8*radial_speed
    common_16 = common_0*common_11
    common_17 = common_5*potential*speed_squared
    common_18 = common_7*speed_squared
    common_19 = common_1*common_5
    return (
        (-3/10*radial_speed, 1/10,),
        ((1/5)*potential - 3/10*speed_squared, -1/10*radial_speed,),
        ((1/2)*common_0, (3/10)*potential - 2/5*speed_squared,),
        (-4/5*common_1 - 3/2*common_3 + (9/10)*potential*speed_squared, common_0,),
        (-57/10*common_4 + (15/2)*common_6 + (17/5)*common_7, -9/5*common_1 - 9/2*common_3 + (19/10)*potential*speed_squared,),
        ((36/5)*common_1*speed_squared - 105/2*common_10 - 57/10*common_12 - 27*common_13 + 51*common_2*potential*speed_squared - 8/5*common_8, -102/5*common_4 + 30*common_6 + (79/5)*common_7,),
        ((945/2)*common_14 + 35*common_15 + (261/2)*common_16 - 567*common_17 - 156*common_18 + 294*common_19, (217/5)*common_1*speed_squared - 525/2*common_10 - 261/10*common_12 - 171*common_13 + 243*common_2*potential*speed_squared - 87/5*common_8,),
        (-1302/5*common_1*common_11 + 3066*common_1*common_2*speed_squared - 3906*common_1*common_9 - 5229/2*common_11*common_3 - 714*common_2*common_8 + (738/5)*common_8*speed_squared + (14931/2)*common_9*potential*speed_squared - 88/5*potential**4 - 10395/2*potential*radial_speed**6 + (261/2)*potential*speed_squared**3, 2835*common_14 + 412*common_15 + 747*common_16 - 3318*common_17 - 1140*common_18 + 2226*common_19,),
    )

def w2(potential, radial_speed, speed_squared, nu, log_one, log_infrared):
    common_0 = (5/7)*nu
    common_1 = (2/7)*nu
    common_2 = common_1 - 2/21
    common_3 = nu*potential
    common_4 = nu*radial_speed
    common_5 = potential*radial_speed
    common_6 = common_3*radial_speed
    common_7 = nu*speed_squared
    common_8 = potential*speed_squared
    common_9 = potential**2
    common_10 = common_9*nu
    common_11 = radial_speed**2
    common_12 = common_11*common_3
    common_13 = radial_speed**3
    return (
        (-common_0*radial_speed + (5/21)*radial_speed, common_2, 0,),
        (-common_0*speed_squared + (3/7)*common_3 - 1/7*potential + (5/21)*speed_squared, -6/7*common_4 + (2/7)*radial_speed, common_2,),
        (-13/21*common_5 + (13/7)*common_6, (8/7)*common_3 - 16/7*common_7 - 8/21*potential + (16/21)*speed_squared, -common_1*radial_speed + (2/21)*radial_speed,),
        (-3*common_10 + (13/7)*common_11*potential - 39/7*common_12 - 29/21*common_8 + common_9 + (29/7)*nu*potential*speed_squared, -18/7*common_5 + (54/7)*common_6, (10/7)*common_3 - 18/7*common_7 - 10/21*potential + (6/7)*speed_squared,),
        ((195/7)*common_13*common_3 - 65/7*common_13*potential + (50/7)*common_4*common_9 + (55/7)*common_5*speed_squared - 165/7*common_6*speed_squared - 50/21*common_9*radial_speed, -116/7*common_10 + (80/7)*common_11*potential - 240/7*common_12 - 148/21*common_8 + (116/21)*common_9 + (148/7)*nu*potential*speed_squared, -80/21*common_5 + (80/7)*common_6,),
    )

def x(potential, radial_speed, speed_squared, nu, log_one, log_infrared):
    common_0 = (1/15)*potential
    common_1 = nu*potential
    common_2 = nu*speed_squared
    common_3 = radial_speed**2
    common_4 = common_3*nu
    common_5 = radial_speed*speed_squared
    common_6 = potential**2
    common_7 = speed_squared**2
    common_8 = common_6*nu
    return (
        (-common_0 + (1/30)*common_1 + (1/20)*common_2 + (1/20)*common_3 - 3/20*common_4 - 1/60*speed_squared,),
        ((7/30)*common_1*radial_speed - 1/5*common_2*radial_speed + (1/15)*common_5 - 2/15*potential*radial_speed,),
        ((13/30)*common_2*potential + (1/6)*common_4*potential + (2/15)*common_6 - 1/5*common_7*nu + (1/15)*common_7 - 7/30*common_8 - 1/5*potential*speed_squared,),
        (-common_0*common_5 - 1/2*common_1*radial_speed**3 + (2/15)*common_6*radial_speed - 11/15*common_8*radial_speed + (7/10)*nu*potential*radial_speed*speed_squared,),
    )

def y1(potential, radial_speed, speed_squared, nu, log_one, log_infrared):
    common_0 = potential**2
    common_1 = -2/5*common_0 + (1/5)*potential*speed_squared
    common_2 = potential*speed_squared
    common_3 = common_2*radial_speed
    common_4 = common_0*radial_speed
    common_5 = potential**3
    common_6 = speed_squared**2
    common_7 = common_6*potential
    common_8 = radial_speed**2
    common_9 = common_0*common_8
    common_10 = common_5*radial_speed
    common_11 = common_7*radial_speed
    common_12 = radial_speed**3
    common_13 = common_12*common_2
    common_14 = common_4*speed_squared
    common_15 = common_0*common_12
    common_16 = radial_speed**4
    return (
        ((1/10)*potential + (1/10)*speed_squared, -3/10*radial_speed,),
        (0, (2/5)*potential - 1/5*speed_squared,),
        (common_1, 0,),
        (-3/5*common_3 + (6/5)*common_4, common_1,),
        ((8/5)*common_0*speed_squared - 4/5*common_5 - 3/5*common_7 + 3*common_8*potential*speed_squared - 6*common_9, -6/5*common_3 + (12/5)*common_4,),
        (12*common_10 + 9*common_11 - 21*common_13 - 24*common_14 + 42*common_15, (26/5)*common_0*speed_squared - 16/5*common_5 - 9/5*common_7 + 9*common_8*potential*speed_squared - 18*common_9,),
        (-378*common_0*common_16 - 156/5*common_0*common_6 + 336*common_0*common_8*speed_squared + 189*common_16*potential*speed_squared - 168*common_5*common_8 + (154/5)*common_5*speed_squared - 126*common_7*common_8 - 44/5*potential**4 + 9*potential*speed_squared**3, 60*common_10 + 36*common_11 - 84*common_13 - 102*common_14 + 168*common_15,),
    )

def y2(potential, radial_speed, speed_squared, nu, log_one, log_infrared):
    common_0 = (2/7)*speed_squared
    common_1 = nu*potential
    common_2 = (6/7)*nu
    common_3 = common_2*speed_squared
    common_4 = nu*radial_speed
    common_5 = potential*radial_speed
    common_6 = common_1*radial_speed
    common_7 = potential*speed_squared
    common_8 = potential**2
    common_9 = common_8*nu
    common_10 = radial_speed**2
    common_11 = common_1*common_10
    common_12 = nu*speed_squared
    common_13 = common_5*speed_squared
    common_14 = radial_speed**3
    common_15 = common_14*potential
    common_16 = common_8*radial_speed
    common_17 = common_6*speed_squared
    common_18 = common_1*common_14
    common_19 = common_4*common_8
    common_20 = potential**3
    common_21 = common_20*nu
    common_22 = radial_speed**4
    common_23 = speed_squared**2
    common_24 = common_8*speed_squared
    common_25 = common_1*common_22
    common_26 = common_1*common_23
    common_27 = common_10*common_7
    common_28 = common_10*common_9
    common_29 = radial_speed**5
    common_30 = common_29*potential
    common_31 = common_20*common_4
    common_32 = common_23*common_5
    common_33 = common_14*common_8
    common_34 = common_18*speed_squared
    common_35 = common_24*common_4
    common_36 = potential**4
    common_37 = radial_speed**6
    common_38 = speed_squared**3
    return (
        (-common_0 + (2/7)*common_1 + common_3 + (1/14)*potential, -18/7*common_4 + (6/7)*radial_speed, common_2 - 2/7,),
        (-5/14*common_5 + (4/7)*common_6, common_0 + (10/7)*common_1 - common_3 - 1/7*potential, -6/7*nu*radial_speed + (2/7)*radial_speed,),
        ((15/14)*common_10*potential - 12/7*common_11 - 9/14*common_7 + (1/2)*common_8 - 2*common_9 + (10/7)*nu*potential*speed_squared, -12/7*common_5 + (22/7)*common_6, (16/7)*common_1 - 12/7*common_12 - 3/7*potential + (4/7)*speed_squared,),
        ((57/14)*common_13 - 75/14*common_15 - 8/7*common_16 - 54/7*common_17 + (60/7)*common_18 + (38/7)*common_19, (51/7)*common_10*potential - 90/7*common_11 - 29/7*common_7 + (25/7)*common_8 - 82/7*common_9 + (66/7)*nu*potential*speed_squared, -17/7*common_5 + (30/7)*common_6,),
        ((15/2)*common_10*common_8 + (450/7)*common_11*speed_squared - 17/7*common_20 + (44/7)*common_21 + (75/2)*common_22*potential + (57/14)*common_23*potential - 15/14*common_24 - 60*common_25 - 54/7*common_26 - 255/7*common_27 - 30*common_28 + (26/7)*common_9*speed_squared, (246/7)*common_13 - 330/7*common_15 - 18*common_16 - 486/7*common_17 + (570/7)*common_18 + 56*common_19, (102/7)*common_10*potential - 180/7*common_11 - 46/7*common_7 + 6*common_8 - 16*common_9 + (96/7)*nu*potential*speed_squared,),
        (270*common_14*common_8*nu + 405*common_14*potential*speed_squared + (155/7)*common_20*radial_speed + (1170/7)*common_23*nu*potential*radial_speed + 540*common_29*nu*potential - 675/2*common_30 - 332/7*common_31 - 1305/14*common_32 - 90*common_33 - 690*common_34 - 774/7*common_35 + (300/7)*common_8*radial_speed*speed_squared, 165*common_10*common_8 + 720*common_10*nu*potential*speed_squared + (8/7)*common_20 - 80/7*common_21 + 405*common_22*potential + (303/7)*common_23*potential - 295/7*common_24 - 690*common_25 - 594/7*common_26 - 390*common_27 - 450*common_28 + (738/7)*common_8*nu*speed_squared, 84*common_13 - 120*common_15 - 58*common_16 - 162*common_17 + 210*common_18 + 144*common_19,),
        (8910*common_1*common_22*speed_squared - 5940*common_1*common_37 + (1170/7)*common_1*common_38 - 180*common_10*common_20 + 288*common_10*common_21 + (3735/2)*common_10*common_23*potential - 1065*common_10*common_24 - 3240*common_11*common_23 - 296/7*common_12*common_20 + (150/7)*common_20*speed_squared - 10665/2*common_22*common_7 + (2745/2)*common_22*common_8 - 3330*common_22*common_9 + (1299/14)*common_23*common_8 - 1350/7*common_23*common_9 + 2376*common_28*speed_squared + (412/7)*common_36*nu - 163/7*common_36 + (7425/2)*common_37*potential - 1305/14*common_38*potential, 7290*common_1*common_29 + 5160*common_14*common_7 + 5040*common_14*common_9 + (7752/7)*common_16*speed_squared - 654/7*common_20*radial_speed + (15390/7)*common_23*common_6 - 4320*common_30 + (2704/7)*common_31 - 8280/7*common_32 - 2100*common_33 - 9180*common_34 - 17712/7*common_35, 705*common_10*common_8 + 2160*common_10*nu*potential*speed_squared + (414/7)*common_20 - 1088/7*common_21 + 1245*common_22*potential + (891/7)*common_23*potential - 1289/7*common_24 - 2160*common_25 - 1728/7*common_26 - 1170*common_27 - 1620*common_28 + (2880/7)*common_8*nu*speed_squared,),
    )

def z1(potential, radial_speed, speed_squared, nu, log_one, log_infrared):
    common_0 = (3/5)*nu
    common_1 = common_0*radial_speed
    common_2 = (1/5)*potential
    common_3 = potential**2
    common_4 = radial_speed**2
    common_5 = common_4*potential
    common_6 = radial_speed*speed_squared
    common_7 = radial_speed**3
    common_8 = common_7*potential
    common_9 = common_3*radial_speed
    common_10 = (27/5)*nu
    common_11 = common_8*nu
    common_12 = (24/5)*nu
    common_13 = potential**3
    common_14 = radial_speed**4
    common_15 = speed_squared**2
    common_16 = common_15*potential
    common_17 = potential*radial_speed**5
    common_18 = common_13*radial_speed
    common_19 = common_16*radial_speed
    common_20 = common_9*speed_squared
    common_21 = common_3*common_7
    return (
        (-common_1 + (1/5)*radial_speed,),
        (common_0*potential - common_0*speed_squared - common_2 + (1/5)*speed_squared,),
        (common_1*potential - 1/5*potential*radial_speed,),
        (-common_0*common_3 - common_2*speed_squared + (1/5)*common_3 + (3/5)*common_4*potential - 9/5*common_5*nu + (3/5)*nu*potential*speed_squared,),
        (-common_10*common_6*potential + 9*common_11 + common_12*common_9 + (9/5)*common_6*potential - 3*common_8 - 8/5*common_9,),
        (-common_10*common_16 - common_12*common_13 + (8/5)*common_13 - 63*common_14*nu*potential + 21*common_14*potential + (9/5)*common_15*potential - 45*common_3*common_4*nu + 15*common_3*common_4 + (51/5)*common_3*nu*speed_squared - 17/5*common_3*speed_squared + 54*common_4*nu*potential*speed_squared - 18*common_5*speed_squared,),
        (-630*common_11*speed_squared + 567*common_17*nu - 189*common_17 + (516/5)*common_18*nu - 172/5*common_18 + 135*common_19*nu - 45*common_19 - 1188/5*common_20*nu + (396/5)*common_20 + 504*common_21*nu - 168*common_21 + 210*common_8*speed_squared,),
    )

def i2(potential, radial_speed, speed_squared, nu, log_one, log_infrared):
    common_0 = nu*potential
    common_1 = nu*speed_squared
    common_2 = nu*radial_speed
    common_3 = potential*radial_speed
    common_4 = common_0*radial_speed
    common_5 = potential*speed_squared
    common_6 = potential**2
    common_7 = common_6*nu
    common_8 = common_1*potential
    common_9 = radial_speed**2
    common_10 = common_0*common_9
    common_11 = -8*potential
    common_12 = common_3*speed_squared
    common_13 = radial_speed**3
    common_14 = common_13*potential
    common_15 = common_6*radial_speed
    common_16 = common_1*common_3
    common_17 = common_0*common_13
    common_18 = common_7*radial_speed
    common_19 = common_9*potential
    common_20 = potential**3
    common_21 = common_20*nu
    common_22 = radial_speed**4
    common_23 = speed_squared**2
    common_24 = common_6*speed_squared
    common_25 = common_0*common_22
    common_26 = common_5*common_9
    common_27 = common_1*common_19
    common_28 = common_7*common_9
    common_29 = 60*common_3
    common_30 = radial_speed**5
    common_31 = common_30*potential
    common_32 = common_20*radial_speed
    common_33 = common_23*common_3
    common_34 = common_13*common_6
    common_35 = common_23*common_4
    common_36 = common_1*common_14
    common_37 = common_18*speed_squared
    common_38 = common_22*potential
    common_39 = 240*common_25
    common_40 = common_23*potential
    common_41 = common_6*common_9
    common_42 = potential**4
    common_43 = common_42*nu
    common_44 = radial_speed**6
    common_45 = speed_squared**3
    common_46 = common_45*potential
    common_47 = common_20*speed_squared
    common_48 = common_0*common_44
    common_49 = common_0*common_45
    common_50 = common_22*common_5
    common_51 = common_22*common_7
    common_52 = common_23*common_7
    common_53 = common_21*common_9
    common_54 = common_41*speed_squared
    common_55 = common_0*common_30
    common_56 = common_2*common_20
    common_57 = common_13*common_5
    common_58 = common_15*speed_squared
    common_59 = common_13*common_7
    common_60 = common_7*speed_squared
    common_61 = radial_speed**7
    common_62 = common_61*potential
    common_63 = common_42*radial_speed
    common_64 = common_0*common_61
    common_65 = common_2*common_42
    common_66 = common_3*common_45
    common_67 = common_30*common_5
    common_68 = common_32*speed_squared
    common_69 = common_30*common_6
    common_70 = common_13*common_20
    common_71 = common_4*common_45
    common_72 = common_1*common_31
    common_73 = common_1*common_32
    common_74 = common_30*common_7
    common_75 = common_13*common_21
    common_76 = common_14*common_23
    common_77 = common_15*common_23
    common_78 = common_13*common_24
    common_79 = common_17*common_23
    common_80 = common_18*common_23
    common_81 = common_13*common_60
    common_82 = common_44*potential
    common_83 = common_19*common_23
    common_84 = common_22*common_6
    common_85 = common_23*common_6
    common_86 = common_20*common_9
    common_87 = potential**5
    common_88 = radial_speed**8
    common_89 = speed_squared**4
    return (
        (1, 0, 0,),
        ((8/7)*common_0 - 29/14*common_1 - 5/7*potential + (29/42)*speed_squared, (12/7)*common_2 - 4/7*radial_speed, 11/21 - 11/7*nu,),
        (0, 2, 0,),
        (-2/21*common_3 + (9/7)*common_4, (26/7)*common_0 - 17/7*common_1 - 40/21*potential + (17/21)*speed_squared, -10/7*common_2 + (10/21)*radial_speed,),
        (-2*potential, 0, 2,),
        (-6/7*common_10 - 61/21*common_5 + 10*common_6 - common_7 - 16/7*common_8 + (2/7)*common_9*potential, (50/7)*common_3 + (18/7)*common_4, (36/7)*common_0 - 27/7*common_1 - 50/21*potential + (9/7)*speed_squared,),
        (6*common_3, common_11, 0,),
        ((65/7)*common_12 - 10/7*common_14 - 880/21*common_15 + (36/7)*common_16 + (30/7)*common_17 + (54/7)*common_18, -24/7*common_10 - 146/7*common_19 - 110/21*common_5 + (706/21)*common_6 - 48/7*common_7 - 44/7*common_8, -20/7*common_4 + (482/21)*potential*radial_speed,),
        (-30*common_19 + 6*common_5 + 2*common_6, 36*common_3, common_11,),
        ((2/7)*common_20 - 34/7*common_21 + 10*common_22*potential + (36/7)*common_23*nu*potential + (65/7)*common_23*potential - 923/21*common_24 - 30*common_25 - 355/7*common_26 - 90/7*common_27 - 81*common_28 + 282*common_6*common_9 + (104/7)*common_6*nu*speed_squared, -52/7*common_12 + (710/7)*common_14 - 244*common_15 + (156/7)*common_16 + (180/7)*common_17 + 84*common_18, -628/7*common_19 + (32/3)*common_6 - 4*common_7 - 64/7*common_8 + (36/7)*common_9*nu*potential + (124/7)*potential*speed_squared,),
        (-90*common_12 + 210*common_14, -240*common_19 - 16*common_6 + 48*potential*speed_squared, common_29,),
        (870*common_13*common_6*nu + 395*common_13*potential*speed_squared + (408/7)*common_20*nu*radial_speed + 270*common_30*nu*potential - 90*common_31 - 5924/21*common_32 - 1035/7*common_33 - 2536*common_34 - 360/7*common_35 - 30*common_36 - 1794/7*common_37 + (6716/7)*common_6*radial_speed*speed_squared, (6368/21)*common_20 - 320/7*common_21 + (228/7)*common_23*nu*potential + (78/7)*common_23*potential - 7978/21*common_24 - 60*common_27 - 912*common_28 - 690*common_38 - common_39 + 2206*common_6*common_9 + (348/7)*common_6*nu*speed_squared + 240*common_9*potential*speed_squared, common_1*common_29 - 240*common_12 + 550*common_14 - 620/3*common_15 + 140*common_18,),
        (16*common_20 + 42*common_24 + 1260*common_26 - 1890*common_38 - 90*common_40 - 210*common_41, -900*common_12 + 2100*common_14 + 360*common_15, -540*common_19 - 76*common_6 + 108*potential*speed_squared,),
        (6862*common_20*common_9 + (1602/7)*common_20*nu*speed_squared + 28230*common_22*common_6 + 1620*common_22*nu*potential*speed_squared + (7967/7)*common_23*common_6 + 270*common_23*common_9*nu*potential + 2220*common_23*common_9*potential - 596/7*common_42 - 312/7*common_43 + 990*common_44*potential - 1035/7*common_46 - 21286/21*common_47 - 2970*common_48 - 360/7*common_49 - 4005*common_50 - 10695*common_51 - 780/7*common_52 - 1314*common_53 - 17407*common_54 + 4863*common_6*common_9*nu*speed_squared, 6030*common_31 - 49320/7*common_32 + (900/7)*common_33 - 24420*common_34 - 2700/7*common_35 - 600*common_36 - 15828/7*common_37 + 2700*common_55 + (9696/7)*common_56 - 3650*common_57 + (65840/7)*common_58 + 11160*common_59, (15748/21)*common_20 - 460/7*common_21 + (648/7)*common_23*nu*potential - 2846/7*common_24 - 360*common_27 - 2142*common_28 - 4540*common_38 - common_39 - 1602/7*common_40 + 2996*common_6*common_9 - 352/7*common_60 + 3090*common_9*potential*speed_squared,),
        (5040*common_13*common_6 + 3150*common_23*potential*radial_speed + 20790*common_30*potential - 168*common_32 - 18900*common_57 - 2016*common_58, -176*common_20 - 22680*common_38 - 1080*common_40 - 6720*common_41 + 1128*common_6*speed_squared + 15120*common_9*potential*speed_squared, -2520*common_12 + 5880*common_14 + 1680*common_15,),
        (-12870*common_62 - 77776/21*common_63 + 38610*common_64 + (8952/7)*common_65 + 5475*common_66 + 49995*common_67 + (343380/7)*common_68 - 372960*common_69 - 142608*common_70 + 900*common_71 - 35640*common_72 - 65988/7*common_73 + 150930*common_74 + 31344*common_75 - 36000*common_76 - 352104/7*common_77 + 320880*common_78 + 4050*common_79 + (53982/7)*common_80 - 97740*common_81, 142152*common_20*common_9 + (15488/7)*common_20*nu*speed_squared + 320490*common_22*common_6 + 22140*common_22*nu*potential*speed_squared + 54990*common_22*potential*speed_squared + (88866/7)*common_23*common_6 + 1440*common_23*common_9*nu*potential + (130360/21)*common_42 - 6488/7*common_43 - 1170/7*common_46 - 432800/21*common_47 - 35640*common_48 - 3420/7*common_49 - 156480*common_51 - 1620/7*common_52 - 34056*common_53 - 196884*common_54 + 61272*common_6*common_9*nu*speed_squared - 64350*common_82 - 7410*common_83, 46890*common_31 - 452296/21*common_32 + (52170/7)*common_33 - 45048*common_34 - 10980/7*common_35 + 960*common_36 - 29160/7*common_37 + 4860*common_55 + (33640/7)*common_56 - 43440*common_57 + (106152/7)*common_58 + 33000*common_59,),
        (344*common_42 + 3150*common_46 + 720*common_47 + 311850*common_50 + 64260*common_54 - 270270*common_82 - 85050*common_83 - 103950*common_84 - 4086*common_85 - 2520*common_86, 291060*common_31 + 9072*common_32 + 44100*common_33 + 126000*common_34 - 264600*common_57 - 47376*common_58, -1856*common_20 - 75600*common_38 - 3600*common_40 - 32760*common_41 + 5328*common_6*speed_squared + 50400*common_9*potential*speed_squared,),
        (-579150*common_0*common_88 - 157275*common_19*common_45 + 2990814*common_20*common_22 + (670794/7)*common_20*common_23 + (2369700/7)*common_20*common_9*nu*speed_squared - 746910*common_21*common_22 - 70722/7*common_21*common_23 + 645975*common_22*common_23*potential - 6297795*common_22*common_24 + 2082510*common_22*common_6*nu*speed_squared - 222750*common_23*common_25 - 2315475/7*common_23*common_28 + (11481300/7)*common_23*common_6*common_9 + (1741648/7)*common_42*common_9 + (54244/7)*common_42*nu*speed_squared - 742484/21*common_42*speed_squared - 335388/7*common_43*common_9 - 740025*common_44*common_5 + 5691690*common_44*common_6 - 2419065*common_44*common_7 + 733590*common_44*nu*potential*speed_squared + (7632/7)*common_45*common_6*nu - 403749/7*common_45*common_6 + 4050*common_45*common_9*nu*potential - 11735184/7*common_86*speed_squared - 1040*common_87*nu - 3880*common_87 + 193050*common_88*potential + 900*common_89*nu*potential + 5475*common_89*potential, 810810*common_62 - 2110448/7*common_63 + 540540*common_64 + (354240/7)*common_65 - 2700*common_66 - 891000*common_67 + (6898056/7)*common_68 - 4863060*common_69 - 2922504*common_70 + 8100*common_71 - 528660*common_72 - 1131912/7*common_73 + 2493180*common_74 + 805824*common_75 + 214650*common_76 - 4462452/7*common_77 + 4140240*common_78 + 83700*common_79 + (638172/7)*common_80 - 1504080*common_81, 1104*common_1*common_20 + 37800*common_1*common_38 + 15300*common_10*common_23 + 180612*common_28*speed_squared + (723776/21)*common_42 - 16608/7*common_43 + (51000/7)*common_46 - 69072*common_47 - 89100*common_48 - 14400/7*common_49 + 680400*common_50 - 552600*common_51 + (33120/7)*common_52 - 150336*common_53 - 415836*common_54 - 580140*common_82 - 189900*common_83 + 741000*common_84 + (160488/7)*common_85 + 511368*common_86,),
    )

def i3(potential, radial_speed, speed_squared, nu, log_one, log_infrared):
    common_0 = nu*potential
    common_1 = (13/6)*common_0
    common_2 = 2*nu
    common_3 = common_2 - 1
    common_4 = potential*radial_speed
    common_5 = 3*potential
    common_6 = potential*speed_squared
    common_7 = potential**2
    common_8 = common_0*speed_squared
    common_9 = common_7*nu
    common_10 = radial_speed**2
    common_11 = common_10*potential
    common_12 = common_0*radial_speed
    common_13 = nu*radial_speed
    common_14 = 9*common_4
    common_15 = 12*potential
    common_16 = radial_speed**3
    common_17 = common_16*potential
    common_18 = common_9*radial_speed
    common_19 = nu*speed_squared
    common_20 = common_0*common_10
    common_21 = 60*potential
    common_22 = potential**3
    common_23 = common_22*nu
    common_24 = radial_speed**4
    common_25 = common_24*potential
    common_26 = speed_squared**2
    common_27 = common_7*speed_squared
    common_28 = common_0*common_24
    common_29 = common_0*common_26
    common_30 = common_19*common_7
    common_31 = common_11*speed_squared
    common_32 = common_10*common_7
    common_33 = common_10*common_9
    common_34 = common_0*common_16
    common_35 = common_4*speed_squared
    common_36 = common_7*radial_speed
    common_37 = radial_speed**5
    common_38 = common_37*potential
    common_39 = common_0*common_37
    common_40 = common_36*speed_squared
    common_41 = common_12*common_26
    common_42 = common_16*common_9
    common_43 = common_26*potential
    common_44 = common_10*common_8
    common_45 = potential**4
    common_46 = radial_speed**6
    common_47 = speed_squared**3
    return (
        (-1, 0, 0, 0,),
        (-common_1 + (19/6)*nu*speed_squared + (5/6)*potential - 5/6*speed_squared, -common_2*radial_speed + radial_speed, common_3, 0,),
        (0, -3, 0, 0,),
        (-common_1*radial_speed - 1/6*common_4, -17/2*common_0 + (15/2)*nu*speed_squared + (7/2)*potential - 3/2*speed_squared, 0, common_3,),
        (common_5, 0, -6, 0,),
        (common_11*common_2 + (1/2)*common_11 + (13/3)*common_6 - 46/3*common_7 - 2/3*common_8 + (14/3)*common_9, -7*common_12 - 13*common_4, -23*common_0 + 15*nu*speed_squared + 10*potential - 3*speed_squared, 4*common_13 - 2*radial_speed,),
        (-common_14, 21*potential, 0, -6,),
        (-common_15*radial_speed*speed_squared - 10*common_16*nu*potential - 5/2*common_17 - 43/3*common_18 + (194/3)*common_7*radial_speed + 6*nu*potential*radial_speed*speed_squared, (81/2)*common_11 - common_19*common_5 + 9*common_20 + 18*common_6 - 101*common_7 + 43*common_9, -9*common_12 - 72*common_4, -27*common_0 + 19*nu*speed_squared + 12*potential - 5*speed_squared,),
        (45*common_10*potential - 9*common_6 - 12*common_7, -90*common_4, common_21, 0,),
        (-common_10*common_19*common_21 - common_15*common_26 + (253/3)*common_22 - 14/3*common_23 + (35/2)*common_25 + (140/3)*common_27 + 70*common_28 + 6*common_29 - 160/3*common_30 + (105/2)*common_31 - 433*common_32 + 131*common_33, -common_14*speed_squared - 210*common_17 - 251*common_18 - 75*common_34 + 709*common_7*radial_speed + 45*nu*potential*radial_speed*speed_squared, 297*common_11 + 18*common_20 - common_5*speed_squared - 238*common_7 - 18*common_8 + 140*common_9, 16*common_12 - 146*common_4,),
        (-315*common_17 + 135*common_35 + 90*common_36, 585*common_10*potential - 117*common_6 - 66*common_7, -360*common_4, common_21,),
        (-595/2*common_16*common_6 + 3859*common_16*common_7 + 700*common_16*nu*potential*speed_squared + (277/3)*common_22*nu*radial_speed - 1796/3*common_22*radial_speed + 165*common_26*potential*radial_speed - 315/2*common_38 - 630*common_39 - 1248*common_40 - 150*common_41 - 1340*common_42 + 723*common_7*nu*radial_speed*speed_squared, 140*common_22 + 17*common_23 + (3045/2)*common_25 + 834*common_27 + 735*common_28 + 63*common_29 - 510*common_30 - 855/2*common_31 - 6219*common_32 + 2403*common_33 - 45*common_43 - 630*common_44, -1905*common_17 - 1350*common_18 - 240*common_34 + 2700*common_7*radial_speed + 180*nu*potential*radial_speed*speed_squared + 585*potential*radial_speed*speed_squared, 735*common_10*potential - common_2*common_6 - 30*common_20 - 149*common_6 + 124*common_7*nu - 92*common_7,),
        (-24*common_22 + 2835*common_24*potential + 135*common_26*potential - 1890*common_31 - 630*common_32 + 72*common_7*speed_squared, -5040*common_17 + 2160*common_35 + 450*common_36, 2970*common_11 - 594*common_6 + 48*common_7, -540*common_4,),
        (-150*common_0*common_47 + 2526*common_10*common_22 - 984*common_10*common_23 + 3150*common_10*common_26*nu*potential + 23931*common_10*common_7*speed_squared - 4095/2*common_11*common_26 - 470/3*common_19*common_22 - 11187*common_19*common_32 - 410/3*common_22*speed_squared + 16170*common_24*common_7*nu - 42510*common_24*common_7 - 9450*common_24*common_8 + 1890*common_24*potential*speed_squared + 756*common_26*common_7*nu - 1386*common_26*common_7 - 184/3*common_45*nu + (1664/3)*common_45 + 6930*common_46*nu*potential + (3465/2)*common_46*potential + 165*common_47*potential, -364*common_13*common_22 + 66924*common_16*common_7 + 9450*common_16*nu*potential*speed_squared + 8190*common_16*potential*speed_squared + 3194*common_22*radial_speed - 135*common_26*common_4 - 14175*common_38 - 8505*common_39 - 22941*common_40 - 2025*common_41 - 27630*common_42 + 10737*common_7*nu*radial_speed*speed_squared, -2864*common_22 + 652*common_23 + 16380*common_25 + 4410*common_27 + 3150*common_28 + 306*common_29 - 2004*common_30 - 9495*common_31 - 30618*common_32 + 15246*common_33 + 495*common_43 - 2880*common_44, -5580*common_17 - 2142*common_18 - 90*common_34 + 2616*common_7*radial_speed + 126*nu*potential*radial_speed*speed_squared + 2502*potential*radial_speed*speed_squared,),
    )

def j2(potential, radial_speed, speed_squared, nu, log_one, log_infrared):
    common_0 = nu*potential
    common_1 = (15/7)*common_0
    common_2 = nu*radial_speed
    common_3 = potential*radial_speed
    common_4 = common_0*radial_speed
    common_5 = potential*speed_squared
    common_6 = potential**2
    common_7 = common_0*speed_squared
    common_8 = common_6*nu
    common_9 = radial_speed**2
    common_10 = common_9*potential
    common_11 = common_0*common_9
    common_12 = radial_speed**3
    common_13 = common_12*potential
    common_14 = common_0*common_12
    common_15 = common_2*common_6
    common_16 = potential**3
    common_17 = common_16*nu
    common_18 = radial_speed**4
    common_19 = common_18*potential
    common_20 = speed_squared**2
    common_21 = common_20*potential
    common_22 = common_6*speed_squared
    common_23 = common_0*common_18
    common_24 = common_0*common_20
    common_25 = nu*speed_squared
    common_26 = common_25*common_6
    common_27 = common_5*common_9
    common_28 = common_6*common_9
    common_29 = common_7*common_9
    common_30 = common_8*common_9
    common_31 = common_6*radial_speed
    common_32 = radial_speed**5
    common_33 = common_32*potential
    common_34 = common_0*common_32
    common_35 = common_16*common_2
    common_36 = common_20*common_3
    common_37 = common_20*common_4
    common_38 = common_12*common_8
    common_39 = potential**4
    common_40 = radial_speed**6
    common_41 = speed_squared**3
    return (
        (-1, 0,),
        (-common_1 + (17/7)*nu*speed_squared - 27/14*potential - 13/28*speed_squared, (5/14)*common_2 - 5/28*radial_speed,),
        (0, -1,),
        (-27/28*common_3 - 15/14*common_4, -5/2*common_0 + (39/14)*nu*speed_squared - 7/4*potential - 9/14*speed_squared,),
        (potential, 0,),
        ((81/28)*common_10 + (12/7)*common_11 + (19/28)*common_5 - 9/7*common_6 - 6/7*common_7 + (11/7)*common_8, -83/14*common_3 - 1/7*common_4,),
        (-3*common_3, potential,),
        (-405/28*common_13 - 60/7*common_14 - 69/7*common_15 + (111/14)*common_6*radial_speed + 6*nu*potential*radial_speed*speed_squared + (15/4)*potential*radial_speed*speed_squared, common_1*common_9 + (579/28)*common_10 - 21/4*common_5 + (65/14)*common_6 - common_7 + (12/7)*common_8,),
        (15*common_10 - 3*common_5 + 2*common_6, -6*common_3,),
        (-144/7*common_16 + (29/7)*common_17 + (405/4)*common_19 + (15/4)*common_21 + (80/7)*common_22 + 60*common_23 + 6*common_24 - 62/7*common_26 - 435/7*common_27 - 789/14*common_28 - 390/7*common_29 + (558/7)*common_30, -825/7*common_13 - 135/7*common_14 - 23*common_15 - 67/2*common_31 + (93/7)*nu*potential*radial_speed*speed_squared + (426/7)*potential*radial_speed*speed_squared,),
        (-105*common_13 - 30*common_31 + 45*potential*radial_speed*speed_squared, 45*common_10 - 9*common_5 + 8*common_6,),
        (468*common_12*common_6 + 630*common_12*nu*potential*speed_squared + 840*common_12*potential*speed_squared + (4863/14)*common_16*radial_speed - 1308/7*common_31*speed_squared - 3645/4*common_33 - 540*common_34 - 810/7*common_35 - 4005/28*common_36 - 990/7*common_37 - 810*common_38 + (1755/7)*common_6*nu*radial_speed*speed_squared, -155/14*common_16 + (106/7)*common_17 + (3705/4)*common_19 + (1809/28)*common_21 - 1077/14*common_22 + 195*common_23 + (135/7)*common_24 - 190/7*common_26 - 720*common_27 + (609/2)*common_28 - 180*common_29 + 276*common_30,),
        (22*common_16 + 945*common_19 + 45*common_21 - 66*common_22 - 630*common_27 + 420*common_28, -420*common_13 - 150*common_31 + 180*potential*radial_speed*speed_squared,),
        (5940*common_0*common_40 - 990/7*common_0*common_41 + (14085/4)*common_10*common_20 + 2880*common_11*common_20 - 1409/7*common_16*common_25 - 5421*common_16*common_9 + (5440/7)*common_16*speed_squared + 2421*common_17*common_9 - 48465/4*common_18*common_5 - 4410*common_18*common_6 - 8370*common_18*common_7 + 9915*common_18*common_8 - 1011/7*common_20*common_6 + (1854/7)*common_20*common_8 - 5346*common_25*common_28 + 2571*common_28*speed_squared + (396/7)*common_39*nu - 2970/7*common_39 + (40095/4)*common_40*potential - 4005/28*common_41*potential, -3513*common_12*common_6 + 2670*common_12*nu*potential*speed_squared + 9585*common_12*potential*speed_squared + (2252/7)*common_16*radial_speed - 18495/2*common_33 - 2295*common_34 - 4252/7*common_35 - 26685/14*common_36 - 4185/7*common_37 - 3690*common_38 + (7731/7)*common_6*nu*radial_speed*speed_squared + (12945/7)*common_6*radial_speed*speed_squared,),
    )
