"""Generated balance expressions from PNEvolution.py; do not edit."""
from numba import njit

@njit(cache=True)
def balance6(Cons, Vars):
    Flux = Vars.Fcal_coeff*(Cons.Fcal_0 + Vars.v**2*(Cons.Fcal_2 + Vars.v*(Cons.Fcal_3 + Vars.Fcal_SO_3 + Vars.v*(Cons.Fcal_4 + Vars.Fcal_SQ_4 + Vars.v*(Cons.Fcal_5 + Vars.Fcal_SO_5 + Vars.v*(Cons.Fcal_6 + Vars.Fcal_SO_6 + Vars.Fcal_SQ_6 + Cons.Fcal_lnv_6*Vars.logv + Vars.v*(Cons.Fcal_7 + Vars.Fcal_SO_7 + Vars.Fcal_SQ_7 + Vars.v*(Cons.Fcal_8 + Vars.Fcal_SO_8 + Cons.Fcal_lnv_8*Vars.logv + Vars.v*(Cons.Fcal_9 + Cons.Fcal_lnv_9*Vars.logv + Vars.v*(Cons.Fcal_10 + Cons.Fcal_lnv_10*Vars.logv + Vars.v*(Cons.Fcal_11 + Cons.Fcal_lnv_11*Vars.logv + Vars.v*(Cons.Fcal_12 + Cons.Fcal_lnv2_12*Vars.logv**2 + Cons.Fcal_lnv_12*Vars.logv))))))))))))
    dEdV = -0.5*Cons.M*Cons.nu*Vars.v*(2.0*Cons.E_0 + Vars.v**2*(4.0*Cons.E_2 + Vars.v*(5.0*Vars.E_SO_3 + Vars.v*(6.0*Cons.E_4 + 6.0*Vars.E_SQ_4 + Vars.v*(7.0*Vars.E_SO_5 + Vars.v*(8.0*Cons.E_6 + 8.0*Vars.E_SQ_6 + Vars.v*(9.0*Vars.E_SO_7 + 9.0*Vars.E_SQ_7 + Vars.v*(10.0*Cons.E_8 + Cons.E_lnv_8*(10.0*Vars.logv + 1.0) + Vars.v**2*(12.0*Cons.E_10 + Cons.E_lnv_10*(12.0*Vars.logv + 1.0) + Vars.v*(13.0*Cons.E_11 + Vars.v*(14.0*Cons.E_12 + Cons.E_lnv_12*(14.0*Vars.logv + 1.0))))))))))))
    Absorption = Vars.Fcal_coeff*Vars.v**5*(Vars.MDot_Alvi_5 + Cons.MDot_Alvi_8*Vars.v**3)
    return Flux[0], dEdV[0], Absorption[0]
