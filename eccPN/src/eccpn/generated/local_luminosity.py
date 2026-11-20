"""Generated local canonical luminosity through relative 4PN."""

import math


def flux_0(potential, radial_speed, speed_squared, nu):
    return (
        + potential**3*radial_speed**2 * (-88/15*nu**2)
        + potential**3*speed_squared * ((32/5)*nu**2)
    )

def flux_2(potential, radial_speed, speed_squared, nu):
    return (
        + potential**5 * ((32/105)*nu**2*(1 - 4*nu))
        + potential**4*radial_speed**2 * ((16/105)*nu**2*(367 - 15*nu))
        + potential**4*speed_squared * ((64/21)*nu**2*(nu - 17))
        + potential**3*radial_speed**4 * ((2/35)*nu**2*(687 - 620*nu))
        + potential**3*radial_speed**2*speed_squared * ((4/105)*nu**2*(1392*nu - 1487))
        + potential**3*speed_squared**2 * ((2/105)*nu**2*(785 - 852*nu))
    )

def flux_4(potential, radial_speed, speed_squared, nu):
    return (
        + potential**6 * ((16/945)*nu**2*(-56*nu**2 + 1026*nu - 253))
        + potential**5*radial_speed**2 * (-2/945*nu**2*(5376*nu**2 + 9798*nu + 106319))
        + potential**5*speed_squared * ((2/2835)*nu**2*(4368*nu**2 + 81828*nu + 281473))
        + potential**4*radial_speed**4 * ((8/945)*nu**2*(-14290*nu**2 + 60971*nu - 33510))
        + potential**4*radial_speed**2*speed_squared * ((8/105)*nu**2*(2165*nu**2 - 8513*nu + 4987))
        + potential**4*speed_squared**2 * ((8/315)*nu**2*(-1393*nu**2 + 5237*nu - 4446))
        + potential**3*radial_speed**6 * ((4/315)*nu**2*(-8404*nu**2 + 20234*nu - 2501))
        + potential**3*radial_speed**4*speed_squared * ((4/105)*nu**2*(7572*nu**2 - 15207*nu + 2018))
        + potential**3*radial_speed**2*speed_squared**2 * ((4/105)*nu**2*(-6292*nu**2 + 10278*nu - 1719))
        + potential**3*speed_squared**3 * ((4/315)*nu**2*(4430*nu**2 - 5497*nu + 1692))
    )

def flux_5(potential, radial_speed, speed_squared, nu):
    return (
        + potential**6*radial_speed * ((192/175)*nu**3)
        + potential**5*radial_speed**3 * ((182816/1575)*nu**3)
        + potential**5*radial_speed*speed_squared * (-4832/45*nu**3)
        + potential**4*radial_speed**5 * (-17312/35*nu**3)
        + potential**4*radial_speed**3*speed_squared * ((459584/525)*nu**3)
        + potential**4*radial_speed*speed_squared**2 * (-595232/1575*nu**3)
    )

def flux_6(potential, radial_speed, speed_squared, nu, log_one, log_two):
    return (
        + potential**7 * ((8/31185)*nu**2*(-3464*nu**3 - 126798*nu**2 - 418936*nu + 112713))
        + potential**6*radial_speed**2 * ((1/1559250)*nu**2*(-97574400*log_one*nu - 54236160*log_two - 68828800*nu**3 + 310915800*nu**2 - 91631925*math.pi**2*nu + 2736887560*nu + 1383523632))
        + potential**6*speed_squared * ((1/1559250)*nu**2*(146361600*log_one*nu + 81354240*log_two - 3700400*nu**3 + 80460000*nu**2 - 2796768040*nu + 70322175*math.pi**2*nu - 1574943488))
        + potential**5*radial_speed**4 * ((1/51975)*nu**2*(-23728320*log_two - 16506620*nu**3 + 76954850*nu**2 - 7245315*math.pi**2*nu + 17921820*nu + 231255634))
        + potential**5*radial_speed**2*speed_squared * ((1/259875)*nu**2*(162708480*log_two + 99216200*nu**3 - 407922300*nu**2 - 203497650*nu + 44750475*math.pi**2*nu - 1567676396))
        + potential**5*speed_squared**2 * ((1/1559250)*nu**2*(-244062720*log_two - 64713200*nu**3 + 390272100*nu**2 - 57536325*math.pi**2*nu + 191500800*nu + 2394459764))
        + potential**4*radial_speed**6 * ((2/17325)*nu**2*(-4302385*nu**3 + 26025095*nu**2 - 30229425*nu + 5476951))
        + potential**4*radial_speed**4*speed_squared * ((2/10395)*nu**2*(6968200*nu**3 - 33049287*nu**2 + 35158037*nu - 6436413))
        + potential**4*radial_speed**2*speed_squared**2 * ((2/3465)*nu**2*(-1862198*nu**3 + 6699393*nu**2 - 6688869*nu + 1230255))
        + potential**4*speed_squared**3 * ((2/3465)*nu**2*(378040*nu**3 - 955013*nu**2 + 1048531*nu - 185007))
        + potential**3*radial_speed**8 * ((1/1386)*nu**2*(-209632*nu**3 + 1099688*nu**2 - 1433696*nu + 301585))
        + potential**3*radial_speed**6*speed_squared * ((2/3465)*nu**2*(1447280*nu**3 - 5288564*nu**2 + 5179198*nu - 1005979))
        + potential**3*radial_speed**4*speed_squared**2 * ((1/1155)*nu**2*(-1751152*nu**3 + 4709506*nu**2 - 3522149*nu + 613047))
        + potential**3*radial_speed**2*speed_squared**3 * ((2/3465)*nu**2*(1858800*nu**3 - 3754884*nu**2 + 2239826*nu - 346489))
        + potential**3*speed_squared**4 * ((1/6930)*nu**2*(-1630560*nu**3 + 2416340*nu**2 - 1388854*nu + 240945))
    )

def flux_7(potential, radial_speed, speed_squared, nu):
    return (
        + potential**7*radial_speed * ((32/99225)*nu**3*(64817*nu - 144374))
        + potential**6*radial_speed**3 * ((64/99225)*nu**3*(222394*nu + 40349))
        + potential**6*radial_speed*speed_squared * (-32/33075*nu**3*(22553*nu + 229339))
        + potential**5*radial_speed**5 * ((8/3675)*nu**3*(5013565 - 2300058*nu))
        + potential**5*radial_speed**3*speed_squared * ((32/11025)*nu**3*(3222115*nu - 6366253))
        + potential**5*radial_speed*speed_squared**2 * ((8/33075)*nu**3*(30925927 - 17550874*nu))
        + potential**4*radial_speed**7 * ((8/6615)*nu**3*(4702853 - 11650712*nu))
        + potential**4*radial_speed**5*speed_squared * ((8/11025)*nu**3*(48107006*nu - 18052091))
        + potential**4*radial_speed**3*speed_squared**2 * ((8/11025)*nu**3*(13139981 - 38494004*nu))
        + potential**4*radial_speed*speed_squared**3 * ((8/33075)*nu**3*(29413042*nu - 8821111))
    )

def flux_8(potential, radial_speed, speed_squared, nu, log_one, log_two):
    return (
        + potential**8 * ((2/19864845)*nu**2*(-532756224*log_one*nu**2 + 133189056*log_one*nu - 59963904*log_two*nu + 14990976*log_two + 1329664*nu**4 - 335360312*nu**3 - 31026996*math.pi**2*nu**2 + 8077760950*nu**2 + 7756749*math.pi**2*nu + 3367686190*nu - 1341546796))
        + potential**7*radial_speed**2 * ((1/4767562800)*nu**2*(-3270154567680*log_one*nu**2 + 11885307033600*log_one*nu - 848203259904*log_two*nu + 5999056487424*log_two - 652641533040*nu**4 + 1897645061760*nu**3 - 10498253957688*nu**2 + 650340971280*math.pi**2*nu**2 - 131583553270384*nu + 6729255027495*math.pi**2*nu - 50513396971248))
        + potential**7*speed_squared * ((1/3972969000)*nu**2*(5622999782400*log_one*nu**2 - 10936032307200*log_one*nu - 100831610880*log_two*nu - 4792269081600*log_two - 87380857200*nu**4 + 727334244000*nu**3 - 826226200800*math.pi**2*nu**2 + 7328137276880*nu**2 - 3946969701675*math.pi**2*nu + 90486027937768*nu + 41861442180352))
        + potential**6*radial_speed**4 * ((1/529729200)*nu**2*(19534394880*log_one*nu**2 + 1914370698240*log_one*nu + 337061716992*log_two*nu + 4139515860480*log_two - 392921850720*nu**4 + 1717923833920*nu**3 - 227272745700*math.pi**2*nu**2 + 5386905201888*nu**2 - 22667861813232*nu + 2075053015035*math.pi**2*nu - 31186278486496))
        + potential**6*radial_speed**2*speed_squared * ((1/993242250)*nu**2*(-25527902400*log_one*nu**2 - 4657177324800*log_one*nu - 1382773392000*log_two*nu - 9692052289920*log_two + 693658102200*nu**4 - 1706431095600*nu**3 - 13393907703680*nu**2 + 382213806975*math.pi**2*nu**2 - 4415338602675*math.pi**2*nu + 46472922663852*nu + 72692615180448))
        + potential**6*speed_squared**2 * ((1/7945938000)*nu**2*(1491717427200*log_one*nu**2 + 5966869708800*log_one*nu + 5268105722880*log_two*nu + 15144183774720*log_two + 111059480000*nu**4 + 61301979200*nu**3 - 430499569500*math.pi**2*nu**2 + 11372142483520*nu**2 - 64215919478896*nu + 6820164125775*math.pi**2*nu - 114443962622208))
        + potential**5*radial_speed**6 * (nu**2*(1152*log_two*nu + (192240/49)*log_two - 8735219/8190*nu**4 + (370344998/27027)*nu**3 - 6109/7*math.pi**2*nu**2 - 511326653/90090*nu**2 + (238723/448)*math.pi**2*nu + (571186117207/8108100)*nu - 702363704057/22072050))
        + potential**5*radial_speed**4*speed_squared * ((1/15891876000)*nu**2*(-36933152256000*log_two*nu - 111993176655360*log_two + 50804224296000*nu**4 - 448125461968800*nu**3 + 150832890348000*nu**2 + 25038785772000*math.pi**2*nu**2 - 2142689486704432*nu - 13031920076175*math.pi**2*nu + 888461722972032))
        + potential**5*radial_speed**2*speed_squared**2 * ((1/10594584000)*nu**2*(12527658823680*log_two*nu + 35757075594240*log_two - 27754335078400*nu**4 + 176366172987200*nu**3 - 8917158650400*math.pi**2*nu**2 - 54454040208960*nu**2 + 3582326823075*math.pi**2*nu + 776175807288464*nu - 275276020253312))
        + potential**5*speed_squared**3 * ((1/7945938000)*nu**2*(-18819440640*log_two*nu - 2474860227840*log_two + 3518807982200*nu**4 - 19810339830400*nu**3 + 8949761421720*nu**2 + 1019236818600*math.pi**2*nu**2 - 71109054213336*nu - 267111219375*math.pi**2*nu + 18704109794584))
        + potential**4*radial_speed**8 * ((1/3153150)*nu**2*(523115565*nu**4 + 29768629520*nu**3 - 87641025745*nu**2 + 30142228881*nu - 1470546975))
        + potential**4*radial_speed**6*speed_squared * ((1/450450)*nu**2*(1081096705*nu**4 - 16322753840*nu**3 + 35852312875*nu**2 - 11944898911*nu + 625848879))
        + potential**4*radial_speed**4*speed_squared**2 * ((1/135135)*nu**2*(-983746018*nu**4 + 6319709261*nu**3 - 10543843493*nu**2 + 3401965282*nu - 206844747))
        + potential**4*radial_speed**2*speed_squared**3 * ((1/630630)*nu**2*(3858809241*nu**4 - 14407126888*nu**3 + 18635659273*nu**2 - 5797819389*nu + 435737025))
        + potential**4*speed_squared**4 * ((1/3153150)*nu**2*(-4240272785*nu**4 + 9360852620*nu**3 - 10604586055*nu**2 + 3128250653*nu - 198144961))
        + potential**3*radial_speed**10 * ((1/15015)*nu**2*(1994360*nu**4 + 25639360*nu**3 - 117441020*nu**2 + 43404065*nu - 4044729))
        + potential**3*radial_speed**8*speed_squared * ((1/9009)*nu**2*(5148040*nu**4 - 95166598*nu**3 + 280637620*nu**2 - 103349531*nu + 10059068))
        + potential**3*radial_speed**6*speed_squared**2 * ((1/45045)*nu**2*(-187807952*nu**4 + 1100106680*nu**3 - 2134586230*nu**2 + 763721359*nu - 75939146))
        + potential**3*radial_speed**4*speed_squared**3 * ((1/30030)*nu**2*(224160096*nu**4 - 779619264*nu**3 + 1019816990*nu**2 - 342794749*nu + 33890680))
        + potential**3*radial_speed**2*speed_squared**4 * ((1/45045)*nu**2*(-230878512*nu**4 + 553715176*nu**3 - 509690226*nu**2 + 156129980*nu - 14988355))
        + potential**3*speed_squared**5 * ((1/90090)*nu**2*(99345120*nu**4 - 170292092*nu**3 + 127910138*nu**2 - 39107757*nu + 4149616))
    )
