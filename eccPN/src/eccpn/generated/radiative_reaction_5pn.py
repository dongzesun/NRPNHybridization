"""Exact generated SH nonlinear radiative sectors at grade ten."""


def failed_tail_n(potential, radial_speed, speed_squared, nu):
    return (
        + potential**5*radial_speed**2 * (-1072/45*nu**2)
        + potential**5*speed_squared * (-2216/45*nu**2)
        + potential**4*radial_speed**4 * ((3448/3)*nu**2)
        + potential**4*radial_speed**2*speed_squared * (-22408/15*nu**2)
        + potential**4*speed_squared**2 * ((1792/9)*nu**2)
        + potential**3*radial_speed**6 * (3570*nu**2)
        + potential**3*radial_speed**4*speed_squared * (-6076*nu**2)
        + potential**3*radial_speed**2*speed_squared**2 * (2658*nu**2)
        + potential**3*speed_squared**3 * (-904/5*nu**2)
    )

def failed_tail_v(potential, radial_speed, speed_squared, nu):
    return (
        + potential**5*radial_speed * ((1096/15)*nu**2)
        + potential**4*radial_speed**3 * (-1244/3*nu**2)
        + potential**4*radial_speed*speed_squared * ((25204/45)*nu**2)
        + potential**3*radial_speed**5 * (-4970/3*nu**2)
        + potential**3*radial_speed**3*speed_squared * ((7316/3)*nu**2)
        + potential**3*radial_speed*speed_squared**2 * (-3766/5*nu**2)
    )

def memory_n(potential, radial_speed, speed_squared, nu):
    return (
        + potential**6 * (-1004/21*nu**2)
        + potential**5*radial_speed**2 * ((9032/105)*nu**2)
        + potential**5*speed_squared * (-26126/105*nu**2)
        + potential**4*radial_speed**4 * ((726412/105)*nu**2)
        + potential**4*radial_speed**2*speed_squared * (-272788/35*nu**2)
        + potential**4*speed_squared**2 * ((106744/105)*nu**2)
        + potential**3*radial_speed**6 * (15822*nu**2)
        + potential**3*radial_speed**4*speed_squared * (-26438*nu**2)
        + potential**3*radial_speed**2*speed_squared**2 * ((79738/7)*nu**2)
        + potential**3*speed_squared**3 * (-26794/35*nu**2)
    )

def memory_v(potential, radial_speed, speed_squared, nu):
    return (
        + potential**5*radial_speed * (-19280/63*nu**2)
        + potential**4*radial_speed**3 * (-341576/105*nu**2)
        + potential**4*radial_speed*speed_squared * ((97472/35)*nu**2)
        + potential**3*radial_speed**5 * (-6866*nu**2)
        + potential**3*radial_speed**3*speed_squared * (9760*nu**2)
        + potential**3*radial_speed*speed_squared**2 * (-2926*nu**2)
    )

def reaction_square_n(potential, radial_speed, speed_squared, nu):
    return (
        + potential**6 * ((14576/225)*nu**2)
        + potential**5*radial_speed**2 * (-4184/225*nu**2)
        + potential**5*speed_squared * ((23464/225)*nu**2)
        + potential**4*radial_speed**4 * (-21264/5*nu**2)
        + potential**4*radial_speed**2*speed_squared * ((347584/75)*nu**2)
        + potential**4*speed_squared**2 * (-47632/75*nu**2)
        + potential**3*radial_speed**6 * (-11088*nu**2)
        + potential**3*radial_speed**4*speed_squared * ((92624/5)*nu**2)
        + potential**3*radial_speed**2*speed_squared**2 * (-39936/5*nu**2)
        + potential**3*speed_squared**3 * ((2688/5)*nu**2)
    )

def reaction_square_v(potential, radial_speed, speed_squared, nu):
    return (
        + potential**5*radial_speed * ((11136/25)*nu**2)
        + potential**4*radial_speed**3 * ((8336/3)*nu**2)
        + potential**4*radial_speed*speed_squared * (-167824/75*nu**2)
        + potential**3*radial_speed**5 * ((17528/3)*nu**2)
        + potential**3*radial_speed**3*speed_squared * (-41984/5*nu**2)
        + potential**3*radial_speed*speed_squared**2 * ((12728/5)*nu**2)
    )
