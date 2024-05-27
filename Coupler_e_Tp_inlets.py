
from cpstate import State, Fluid
import numpy as np
import CoolProp as CP
from ioutils import *
import ioutils as io
import math
import matplotlib.pyplot as plt
from time import process_time




def coupler(LCESin, PTESin, e, p_loss, c, npt = 100):
    if c == 'c':
        output = io.Dict({'desc':'Output from coupler (charging)'})
    elif c == 'd':
        output = io.Dict({'desc':'Output from coupler (discharging)'})
    if type(LCESin) == Dict:
        LCESin = LCESin.unpack('state2')[0]
    if type(PTESin) == Dict:
        PTESin = PTESin.unpack('state2')[0]

    State.T0 = PTESin.T
    #print("Coupler_new", State.T0)
    #PTESin = State('pT', PTESin.p, PTESin.T, fluid = PTESin.fluid)
    #LCESin = State('pT', LCESin.p, LCESin.T, fluid = LCESin.fluid)
    
    #Finding a valid TQ diagram
    PTESout_T = LCESin.T
    LCESout_T = PTESin.T

    if LCESin.fluid.trl.T > LCESout_T:
        #LCESmin_T = State('pC', LCESin.p, 'f', fluid=LCESin.fluid)
        #print(LCESmin_T.t)

        LCESout_T = LCESin.fluid.trl.T

    hx_valid = False
    iterated = False
    dh_factor = 1
    i_time = 0
    while hx_valid == False:
        h1 = LCESin.h
        p1 = LCESin.p
        LCESout = State('pT', LCESin.p*(1-p_loss), LCESout_T, fluid = LCESin.fluid)
        LCES_h = np.zeros(npt); LCES_T = np.zeros(npt); LCES_s = np.zeros(npt)
        LCES_cp = np.zeros(npt); LCES_dT = np.zeros(npt)
        dh = dh_factor * (LCESout.h - h1) / (npt - 1)
        dp = (LCESin.p * p_loss) / (npt - 1)

        LCESout.update('pT', LCESin.p, LCESin.T)
        for n in range(npt):
            T1 = LCESout.T
            try:
                LCESout.update('hp', h1, p1)
            except ValueError:
                LCESout.update('pT', p1, T1+LCES_dT[n-1])
            LCES_h[n] = LCESout.h; LCES_T[n] = LCESout.T
            LCES_cp[n] = LCESout.fluid.cpmass()
            LCES_dT[n] = LCESout.T - T1
            p1, h1 = p1 - dp, h1 + dh
    
            LCES_h_norm = ((LCES_h - LCES_h[0])/(LCES_h[-1]-LCES_h[0]))*100

        h1 = PTESin.h
        p1 = PTESin.p
        PTESout = State('pT', PTESin.p*(1-p_loss), PTESout_T, fluid = PTESin.fluid)
        PTES_h = np.zeros(npt); PTES_T = np.zeros(npt); PTES_s = np.zeros(npt)
        PTES_cp = np.zeros(npt); PTES_dT = np.zeros(npt)
        dh = dh_factor * (PTESout.h - h1) / (npt - 1)
        dp = (PTESin.p * p_loss) / (npt - 1)

        PTESout.update('pT', PTESin.p, PTESin.T)
        for n in range(npt):
            T1 = PTESout.T
            try:
                PTESout.update('hp', h1, p1)
            except ValueError:
                PTESout.update('pT', p1, T1+PTES_dT[n-1])
            PTES_h[n] = PTESout.h; PTES_T[n] = PTESout.T
            PTES_cp[n] = PTESout.fluid.cpmass()
            PTES_dT[n] = PTESout.T - T1
            p1, h1 = p1 - dp, h1 + dh

        PTES_h = np.flip(PTES_h); PTES_T = np.flip(PTES_T); PTES_s = np.flip(PTES_s)
        PTES_h_norm = ((PTES_h - PTES_h[0])/(PTES_h[-1]-PTES_h[0]))*100

        hx_valid = True
        if c == 'c':
            T_diff = np.round(LCES_T - PTES_T, 3)
        elif c == 'd':
            T_diff = np.round(PTES_T - LCES_T, 3)
        for n in range(len(T_diff)):
            if T_diff[n] < -0.001:
                hx_valid = False
                LCES_top = False
                dh_factor -= 0.05
                iterated = True
                i_time += 1
                break

    output.iteration = i_time
    mr = (LCESout.h - LCESin.h) / (PTESin.h - PTESout.h)
    
    #If iteration is not required to find a valid TQ diagram, pinch point must be at either ends
    if iterated == False:
        #LCES stream Q_max
        qa = 0
        for n in range(len(LCES_dT)):
            qa += LCES_cp[n] * LCES_dT[n]

        #PTES stream Q_max
        qb = 0
        for n in range(len(PTES_dT)):
            qb += PTES_cp[n] * PTES_dT[n]
        qb = mr*qb

        #Q_max being the smallest
        if abs(qa) < abs(qb):
            qmax = abs(qa)
            pinch_index = len(LCES_dT) - 1
        else:
            qmax = abs(qb)
            pinch_index = 0
    
    #If iteration is required to find a valid TQ diagram, pinch point is definitely somewhere in between
    else:
        pinch_index = (T_diff).argmin()
        qmax = 0
        for n in range(len(LCES_dT)):
            if n > pinch_index:
                qmax += mr*abs(PTES_dT[n]*PTES_cp[n])
            else:
                qmax += abs(LCES_dT[n]*LCES_cp[n])

    #Finding TQ diagram using Q_actual = e * Q_max
    qact = e * qmax

    #State.T0 = 298.15

    h1 = LCESin.h
    p1 = LCESin.p
    if c == 'c':
        LCESout = State('hp', h1 - qact, LCESin.p*(1-p_loss), fluid = LCESin.fluid)
    elif c == 'd':
        LCESout = State('hp', h1 + qact, LCESin.p*(1-p_loss), fluid = LCESin.fluid)
    LCES_h = np.zeros(npt); LCES_T = np.zeros(npt); LCES_s = np.zeros(npt); 
    LCES_cp = np.zeros(npt); LCES_dT = np.zeros(npt)
    LCES_dQT = 0
    dh = (LCESout.h - h1) / (npt - 1)
    dp = (LCESin.p * p_loss) / (npt - 1)

    LCESout.update('pT', LCESin.p, LCESin.T)
    for n in range(npt):
        T1 = LCESout.T
        LCES_dQT += dh/T1
        LCESout.update('hp', h1, p1)
        LCES_h1 = h1; LCES_p1 = p1
        LCES_h[n] = LCESout.h; LCES_T[n] = LCESout.T; LCES_s[n] = LCESout.s
        LCES_cp[n] = LCESout.fluid.cpmass()
        LCES_dT[n] = LCESout.T - T1
        p1, h1 = p1 - dp, h1 + dh
    
    LCES_h_norm = ((LCES_h - LCES_h[0])/(LCES_h[-1]-LCES_h[0]))*100

    h1 = PTESin.h
    p1 = PTESin.p
    if c == 'c':
        PTESout = State('hp', h1 + qact/mr, PTESin.p*(1-p_loss), fluid = PTESin.fluid)
    elif c == 'd':
        PTESout = State('hp', h1 - qact/mr, PTESin.p*(1-p_loss), fluid = PTESin.fluid)
    PTES_h = np.zeros(npt); PTES_T = np.zeros(npt); PTES_s = np.zeros(npt)
    PTES_cp = np.zeros(npt); PTES_dT = np.zeros(npt)
    PTES_dQT = 0
    dh = (PTESout.h - h1) / (npt - 1)
    dp = (PTESin.p * p_loss) / (npt - 1)

    PTESout.update('pT', PTESin.p, PTESin.T)
    for n in range(npt):
        T1 = PTESout.T
        PTES_dQT += dh/T1
        PTESout.update('hp', h1, p1)
        PTES_h1 = h1; PTES_p1 = p1
        PTES_h[n] = PTESout.h; PTES_T[n] = PTESout.T; PTES_s[n] = PTESout.s
        PTES_cp[n] = PTESout.fluid.cpmass()
        PTES_dT[n] = PTESout.T - T1
        p1, h1 = p1 - dp, h1 + dh

    PTES_h = np.flip(PTES_h); PTES_T = np.flip(PTES_T); PTES_s = np.flip(PTES_s)
    PTES_h_norm = ((PTES_h - PTES_h[0])/(PTES_h[-1]-PTES_h[0]))*100

    mr = (LCESout.h - LCESin.h) / (PTESin.h - PTESout.h)
    T_diff = np.round(LCES_T - PTES_T, 3)
    pinch_index = (T_diff).argmin()

    liquid_T = LCESout.fluid.T_critical()

    output.mass_r = mr
    
    eta = mr*(PTESout.b - PTESin.b)/(LCESout.b - LCESin.b)
    output.eta = abs(eta)

    if LCESin.fluid.trl.T > PTESin.T:
        PTESin_T = LCESin.fluid.trl.T
    else:
        PTESin_T = PTESin.T
    LCES_min = State('pT', LCESin.p, PTESin_T, fluid = LCESin.fluid)
    zeta = mr*(PTESout.b - PTESin.b)/(LCESin.b - LCES_min.b)
    output.zeta = abs(zeta)

    State.T0 = 298.15
    LCESout.update('hp', LCES_h1, LCES_p1)
    PTESout.update('hp', PTES_h1, PTES_p1)

    output.LCESout = LCESout
    output.LCESin = LCESin
    output.LCESmin = LCES_min
    output.PTESout = PTESout
    output.PTESin = PTESin
    output.LCESdQT = LCES_dQT
    output.PTESdQT = PTES_dQT

    wlirr_LCES = LCESout.T0 * (LCESout.s - LCESin.s)# - LCES_dQT)
    wlirr_PTES = PTESout.T0 * (PTESout.s - PTESin.s)# - PTES_dQT)
    wlirr = (wlirr_LCES + mr * wlirr_PTES)

    output.wlirr_per_unit_LCES_flow = wlirr
    output.wlirr_LCES = wlirr_LCES
    output.wlirr_PTES_per_unit_PTES_flow = wlirr_PTES
    output.b_delta_LCES = LCESout.b - LCESin.b
    output.b_delta_PTES = PTESout.b - PTESin.b
    output.b_delta_PTES_per_unit_LCES_flow = mr * (PTESout.b - PTESin.b)
    output.h_delta_PTES_per_unit_LCES_flow = mr * (PTESout.h - PTESin.h)
    output.b_delta = (LCESout.b - LCESin.b) + mr * (PTESout.b - PTESin.b)
    output.h_delta = (LCESout.h - LCESin.h) + mr * (PTESout.h - PTESin.h)

    h_bal = io.Dict({'desc'                             : 'Coupler Enthalpy Balance',
    'Availability decrease in LCES side'                 : LCESout.b - LCESin.b,
    'Work loss due to irreversibility in LCES side'      : LCESout.T0 * (LCESout.s - LCESin.s),
    'LCES side sum'                                      : LCESout.b - LCESin.b + LCESout.T0 * (LCESout.s - LCESin.s),
    'Enthalpy decrease in LCES side'                     : LCESout.h - LCESin.h,
    'Availability increase in PTES side'                : mr*(PTESout.b - PTESin.b),
    'Work loss due to irreversibility in PTES side'     : mr*(PTESout.T0 * (PTESout.s - PTESin.s)),
    'PTES side sum'                                     : mr*(PTESout.b - PTESin.b + PTESout.T0 * (PTESout.s - PTESin.s)),
    'Enthalpy increase in PTES side'                    : mr*(PTESout.h - PTESin.h)})
    output.h_bal = h_bal

    
    output.LCES_T = LCES_T; output.LCES_s = LCES_s; output.LCES_h = LCES_h
    output.PTES_T = PTES_T; output.PTES_s = PTES_s; output.PTES_h = PTES_h
    
    font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 16}

    plt.rc('font', **font)
    plt.rcParams["figure.figsize"]  = [8,8]
    #plt.axhline(liquid_T, label = "Critical temperature of CO2")
    plt.plot(PTES_h_norm, PTES_T, c = "b", label = 'Nitrogen')
    plt.plot(LCES_h_norm, LCES_T, c = "r", label = 'CO$_{2}$')
    plt.scatter([pinch_index + pinch_index*(1/(1+npt)), pinch_index + pinch_index*(1/(1+npt))], [LCES_T[pinch_index], PTES_T[pinch_index]], label = "Pinch Point")
    plt.legend()
    #for n in range(len(PTES_T)):
        #plt.plot([LCES_h_norm[n],PTES_h_norm[n]], [LCES_T[n], PTES_T[n]])
    #plt.title("Pinch point temp diff:" + str(min(T_diff)))
    #plt.legend()
    plt.savefig("test.png", bbox_inches='tight')
    plt.ylabel('Temperature, K')
    plt.xlabel('q$_{0}$, %')
    plt.xlim(1, 100)
    plt.show()

    #print("LCES mCp:", LCES_cp)
    #print("LCES mCp:", mr*PTES_cp)
    #print(qa)
    #print(qb)

    return output