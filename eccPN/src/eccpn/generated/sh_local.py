"""Generated local SH force polynomials. Run generate/sh_local.py."""

import math


def a1(potential, radial_speed, speed_squared, nu):
    return (
        + potential * (-2*(nu + 2))
        + radial_speed**2 * (-3/2*nu)
        + speed_squared * (3*nu + 1)
    )

def a2(potential, radial_speed, speed_squared, nu):
    return (
        + potential**2 * (3*((29/4)*nu + 3))
        + potential*radial_speed**2 * (-(2*nu**2 + 25*nu + 2))
        + potential*speed_squared * (nu*(2*nu - 13/2))
        + radial_speed**4 * ((15/8)*nu*(1 - 3*nu))
        + radial_speed**2*speed_squared * (3*nu*(2*nu - 3/2))
        + speed_squared**2 * (nu*(3 - 4*nu))
    )

def a3(potential, radial_speed, speed_squared, nu, log_one):
    return (
        + potential**3 * (-71/2*nu**2 - 437/4*nu + (41/16)*math.pi**2*nu - 16)
        + potential**2*radial_speed**2 * (-110*log_one*nu - 7*nu**3 + (11/8)*nu**2 + (615/64)*math.pi**2*nu + (32573/168)*nu + 1)
        + potential**2*speed_squared * (nu*(22*log_one + nu**2 - 26987/840 - 123/64*math.pi**2))
        + potential*radial_speed**4 * (nu*(-30*nu**2 - 69/2*nu + 79))
        + potential*radial_speed**2*speed_squared * (nu*(20*nu**2 + 16*nu - 121))
        + potential*speed_squared**2 * (nu*(-10*nu**2 + 8*nu + 75/4))
        + radial_speed**6 * ((35/16)*nu*(-5*nu**2 + 5*nu - 1))
        + radial_speed**4*speed_squared * ((15/8)*nu*(17*nu**2 - 18*nu + 4))
        + radial_speed**2*speed_squared**2 * ((3/8)*nu*(-60*nu**2 + 79*nu - 20))
        + speed_squared**3 * (nu*(13*nu**2 - 49/4*nu + 11/4))
    )

def a4_0(radial_speed, speed_squared, nu):
    return (
        + radial_speed**8 * ((315/128)*nu*(-7*nu**3 + 14*nu**2 - 7*nu + 1))
        + radial_speed**6*speed_squared * ((35/16)*nu*(42*nu**3 - 69*nu**2 + 34*nu - 5))
        + radial_speed**4*speed_squared**2 * ((15/16)*nu*(-178*nu**3 + 269*nu**2 - 125*nu + 18))
        + radial_speed**2*speed_squared**3 * (3*nu*(33*nu**3 - 109/2*nu**2 + (397/16)*nu - 7/2))
        + speed_squared**4 * (nu*(-54*nu**3 + 61*nu**2 - 175/8*nu + 21/8))
    )

def a4_1(potential, radial_speed, speed_squared, nu):
    return (
        + potential*radial_speed**6 * (nu*(-86*nu**3 + (181/2)*nu**2 + 407*nu + 2973/40))
        + potential*radial_speed**4*speed_squared * (nu*(228*nu**3 - 81*nu**2 - 1627/2*nu + 1497/32))
        + potential*radial_speed**2*speed_squared**2 * (nu*(-104*nu**3 + 47*nu**2 + (1009/2)*nu - 2583/16))
        + potential*speed_squared**3 * (nu*(58*nu**3 - 44*nu**2 - 58*nu + 1067/32))
    )

def a4_2(potential, radial_speed, speed_squared, nu, log_one, log_two):
    return (
        + potential**2*radial_speed**4 * (nu*(385*log_one*nu - 1540*log_two*nu + 385*log_two - 126*nu**3 - 1959/32*nu**2 - 4305/128*math.pi**2*nu + (326101/96)*nu + (45255/1024)*math.pi**2 + 2094751/960))
        + potential**2*radial_speed**2*speed_squared * (nu*(-275*log_one*nu - 330*log_one + 1100*log_two*nu - 275*log_two + 63*nu**3 - 309/4*nu**2 - 255461/112*nu + (3075/128)*math.pi**2*nu - 1636681/1120 - 12585/512*math.pi**2))
        + potential**2*speed_squared**2 * (nu*(22*log_one*nu + 66*log_one - 88*log_two*nu + 22*log_two - 16*nu**3 + (145/2)*nu**2 - 123/64*math.pi**2*nu + (7263/70)*nu + (1155/1024)*math.pi**2 + 1096941/11200))
    )

def a4_3(potential, radial_speed, speed_squared, nu, log_one, log_two):
    return (
        + potential**3*radial_speed**2 * (3100*log_one*nu**2 - 920*log_one*nu + 384*log_one - 6464*log_two*nu**2 + 3152*log_two*nu - 384*log_two - 24*nu**4 - 17*nu**3 + (1255151/840)*nu**2 + (7095/32)*math.pi**2*nu**2 - 2969/16*math.pi**2*nu + (1297943/8400)*nu - 2)
        + potential**3*speed_squared * (-1616/3*log_one*nu**2 + 80*log_one*nu - 64*log_one + (3232/3)*log_two*nu**2 - 1576/3*log_two*nu + 64*log_two + (19/2)*nu**3 - 229/8*math.pi**2*nu**2 - 693947/2520*nu**2 + (1237279/25200)*nu + (3835/96)*math.pi**2*nu)
    )

def a4_4(potential, nu, log_one, log_two):
    return (
        + potential**4 * (98*log_one*nu**2 - 20*log_one*nu + 16*log_one - 808/3*log_two*nu**2 + (394/3)*log_two*nu - 16*log_two + (3/4)*math.pi**2*nu**2 + (477763/720)*nu**2 - 4543/96*math.pi**2*nu + (6625537/12600)*nu + 25)
    )

def b1(radial_speed, nu):
    return (
        + radial_speed * (2*(nu - 2))
    )

def b2(potential, radial_speed, speed_squared, nu):
    return (
        + potential*radial_speed * (4*nu**2 + (41/2)*nu + 2)
        + radial_speed**3 * (3*nu*(nu + 3/2))
        + radial_speed*speed_squared * (-nu*(2*nu + 15/2))
    )

def b3(potential, radial_speed, speed_squared, nu, log_one):
    return (
        + potential**2*radial_speed * (44*log_one*nu + 8*nu**3 + 25*nu**2 - 123/32*math.pi**2*nu - 18169/840*nu - 4)
        + potential*radial_speed**3 * (nu*(18*nu**2 + (59/2)*nu + 329/6))
        + potential*radial_speed*speed_squared * (-nu*(10*nu**2 + 27*nu + 15))
        + radial_speed**5 * (15*nu*((1/4)*nu**2 + nu - 3/8))
        + radial_speed**3*speed_squared * (3*nu*(-4*nu**2 - 37/4*nu + 4))
        + radial_speed*speed_squared**2 * (nu*(6*nu**2 + 19*nu - 65/8))
    )

def b4_0(radial_speed, speed_squared, nu):
    return (
        + radial_speed**7 * ((35/16)*nu*(2*nu**3 + 11*nu**2 - 14*nu + 3))
        + radial_speed**5*speed_squared * ((15/16)*nu*(-28*nu**3 - 115*nu**2 + 111*nu - 22))
        + radial_speed**3*speed_squared**2 * (3*nu*(18*nu**3 + 43*nu**2 - 623/16*nu + 15/2))
        + radial_speed*speed_squared**3 * (nu*(-24*nu**3 - 69*nu**2 + 54*nu - 157/16))
    )

def b4_1(potential, radial_speed, speed_squared, nu):
    return (
        + potential*radial_speed**5 * (nu*(30*nu**3 + 60*nu**2 - 901/8*nu - 54319/160))
        + potential*radial_speed**3*speed_squared * (nu*(-98*nu**3 - 349/2*nu**2 + (1199/12)*nu + 25943/48))
        + potential*radial_speed*speed_squared**2 * (nu*(44*nu**3 + 118*nu**2 - 389/8*nu - 5725/32))
    )

def b4_2(potential, radial_speed, speed_squared, nu, log_one, log_two):
    return (
        + potential**2*radial_speed**3 * (nu*(-330*log_one*nu + 660*log_one + 880*log_two*nu - 220*log_two + 74*nu**3 + (209/2)*nu**2 - 184613/112*nu + (1845/64)*math.pi**2*nu - 9130111/3360 - 4695/256*math.pi**2))
        + potential**2*radial_speed*speed_squared * (nu*(176*log_one*nu - 264*log_one - 440*log_two*nu + 110*log_two - 34*nu**3 - 70*nu**2 - 123/8*math.pi**2*nu + (58557/70)*nu + (1455/256)*math.pi**2 + 8692601/5600))
    )

def b4_3(potential, radial_speed, nu, log_one, log_two):
    return (
        + potential**3*radial_speed * (-3188/3*log_one*nu**2 + 600*log_one*nu - 128*log_one + (6992/3)*log_two*nu**2 - 3284/3*log_two*nu + 128*log_two + 16*nu**4 + 66*nu**3 - 2201/32*math.pi**2*nu**2 - 28406/45*nu**2 - 619267/525*nu + (791/16)*math.pi**2*nu + 2)
    )

def rr35n(potential, radial_speed, speed_squared, nu):
    return (
        + potential**2*radial_speed * (-92/35*(14*nu + 43))
        + potential*radial_speed**3 * (-2/5*(188*nu + 147))
        + potential*radial_speed*speed_squared * ((4/105)*(1267*nu - 519))
        + radial_speed**5 * (-112)
        + radial_speed**3*speed_squared * (6*(2*nu + 19))
        + radial_speed*speed_squared**2 * (-6*(2*nu + 61/35))
    )

def rr35v(potential, radial_speed, speed_squared, nu):
    return (
        + potential**2 * ((4/105)*(546*nu + 1325))
        + potential*radial_speed**2 * ((2/15)*(424*nu + 205))
        + potential*speed_squared * (-4/105*(777*nu + 205))
        + radial_speed**4 * (120)
        + radial_speed**2*speed_squared * (-6/5*(2*nu + 113))
        + speed_squared**2 * ((2/35)*(42*nu + 313))
    )

def rr45mattern(potential, radial_speed, speed_squared, nu):
    return (
        + potential**3*radial_speed * ((2/945)*(-16344*nu**2 + 278694*nu + 168461))
        + potential**2*radial_speed**3 * ((1/945)*(-270684*nu**2 + 290668*nu - 83177))
        + potential**2*radial_speed*speed_squared * ((1/315)*(13908*nu**2 - 58468*nu - 129769))
        + potential*radial_speed**5 * ((1/105)*(-37296*nu**2 + 111636*nu + 261883))
        + potential*radial_speed**3*speed_squared * ((1/210)*(96168*nu**2 - 285858*nu - 774655))
        + potential*radial_speed*speed_squared**2 * ((1/210)*(-26590*nu**2 + 92674*nu + 214605))
        + radial_speed**7 * (-36*(14*nu + 5))
        + radial_speed**5*speed_squared * (7*(3*nu**2 + 158*nu + 47/2))
        + radial_speed**3*speed_squared**2 * (2*(-33*nu**2 - 2181/7*nu + 44/7))
        + radial_speed*speed_squared**3 * (45*nu**2 + (1248/35)*nu - 1643/210)
    )

def rr45matterv(potential, radial_speed, speed_squared, nu):
    return (
        + potential**3 * ((2/2835)*(11016*nu**2 - 390096*nu - 249643))
        + potential**2*radial_speed**2 * ((1/315)*(69060*nu**2 - 97228*nu + 85991))
        + potential**2*speed_squared * ((1/315)*(-16164*nu**2 + 66632*nu + 47459))
        + potential*radial_speed**4 * ((1/105)*(17880*nu**2 - 31938*nu - 86323))
        + potential*radial_speed**2*speed_squared * ((1/210)*(-39424*nu**2 + 43062*nu + 224351))
        + potential*speed_squared**2 * ((1/210)*(10614*nu**2 - 7546*nu - 32241))
        + radial_speed**6 * (70*(6*nu + 5))
        + radial_speed**4*speed_squared * (-(3*nu**2 + 984*nu + 1007/2))
        + radial_speed**2*speed_squared**2 * ((6/35)*(63*nu**2 + 3625*nu + 963))
        + speed_squared**3 * (-1/630*(4914*nu**2 + 37548*nu + 4873))
    )

def incoming_conservative(speed, nu):
    return (
        (1, 1, 0),
        ((3/2)*nu*speed**2 + speed**2, (7/2)*nu*speed**2 - 3*speed**2, -2*nu - 4),
        (-29/8*nu**2*speed**4 + (3/8)*nu*speed**4, -21/8*nu**2*speed**4 - 21/8*nu*speed**4, 4*nu**2*speed**2 - 11*nu*speed**2),
        ((183/16)*nu**3*speed**6 - 87/16*nu**2*speed**6 + (9/16)*nu*speed**6, (147/16)*nu**3*speed**6 + (13/16)*nu**2*speed**6 - 19/16*nu*speed**6, -12*nu**3*speed**4 - 8*nu**2*speed**4 + (199/12)*nu*speed**4),
        (-6045/128*nu**4*speed**8 + (2125/64)*nu**3*speed**8 - 957/128*nu**2*speed**8 + (67/128)*nu*speed**8, -5005/128*nu**4*speed**8 + (605/64)*nu**3*speed**8 + (403/128)*nu**2*speed**8 - 109/128*nu*speed**8, 72*nu**4*speed**6 + 16*nu**3*speed**6 - 64/3*nu**2*speed**6 + (1811/120)*nu*speed**6),
    )
