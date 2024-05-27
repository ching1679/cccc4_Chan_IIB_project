
from cpstate import State, Fluid
import numpy as np
import CoolProp as CP
from ioutils import *
import ioutils as io
import math
import matplotlib.pyplot as plt
from time import process_time
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 14}

#plt.rc('font', **font)
def HX_dchg(hotin, coldin, e, mr_max, q_max, p_loss, npt = 100):
    #mr defined as mass flow of storage fluid over mass flow of working fluid
    output = io.Dict({'desc':'Output from HX_dchg'})
    if type(hotin) == Dict:
        hotin = hotin.unpack('state2')[0]
    if type(coldin) == Dict:
        coldin = coldin.unpack('state2')[0]
    State.T0 = coldin.T

    q_act = q_max * e

    hot_h = np.zeros(npt); hot_T = np.zeros(npt); hot_s = np.zeros(npt)
    hot_cp = np.zeros(npt); hot_dT = np.zeros(npt)
    
    dh_flow = q_act / npt   #for working fluid
    dh_store = q_act / (npt * mr_max)



def HX(wflin, stflin, e, p_loss, c, npt = 100):
    if c == 'c':
        output = io.Dict({'desc':'Output from HX (charging)'})
    if c == 'd':
        output = io.Dict({'desc':'Output from HX (discharging)'})

    if type(wflin) == Dict:
        wflin = wflin.unpack('state2')[0]
    if type(stflin) == Dict:
        stflin = stflin.unpack('state2')[0]

    #set T0 to stfl in temperature to evaluate Eta and Zeta correctly
    State.T0 = stflin.T
    #print("HX_new", State.T0)
    #stflin = State('pT', stflin.p, stflin.T, fluid = stflin.fluid)
    #wflin = State('pT', wflin.p, wflin.T, fluid = wflin.fluid)
    #Finding a valid TQ diagram
    stflout_T = wflin.T
    wflout_T = stflin.T

    hx_valid = False
    iterated = False
    dh_factor = 1
    i_time = 0
    while hx_valid == False:
        h1 = wflin.h
        p1 = wflin.p
        wflout = State('pT', wflin.p*(1-p_loss), wflout_T, fluid = wflin.fluid)
        wfl_h = np.zeros(npt); wfl_T = np.zeros(npt); wfl_s = np.zeros(npt)
        wfl_cp = np.zeros(npt); wfl_dT = np.zeros(npt)
        dh = dh_factor * (wflout.h - h1) / (npt - 1)
        dp = (wflin.p * p_loss) / (npt - 1)

        wflout.update('pT', wflin.p, wflin.T)
        for n in range(npt):
            T1 = wflout.T
            try:
                wflout.update('hp', h1, p1)
            except ValueError:
                wflout.update('pT', p1, T1+wfl_dT[n-1])
            wfl_h[n] = wflout.h; wfl_T[n] = wflout.T
            wfl_cp[n] = wflout.fluid.cpmass()
            wfl_dT[n] = wflout.T - T1
            p1, h1 = p1 - dp, h1 + dh

        h1 = stflin.h
        p1 = stflin.p
        stflout = State('pT', stflin.p, stflout_T, fluid = stflin.fluid)
        stfl_h = np.zeros(npt); stfl_T = np.zeros(npt); stfl_s = np.zeros(npt)
        stfl_cp = np.zeros(npt); stfl_dT = np.zeros(npt)
        dh = dh_factor * (stflout.h - h1) / (npt - 1)
        dp = 0 # (stflin.p * p_loss) / (npt - 1)

        stflout.update('pT', stflin.p, stflin.T)
        for n in range(npt):
            T1 = stflout.T
            try:
                stflout.update('hp', h1, p1)
            except ValueError:
                stflout.update('pT', p1, T1+stfl_dT[n-1])
            stfl_h[n] = stflout.h; stfl_T[n] = stflout.T
            stfl_cp[n] = stflout.fluid.cpmass()
            stfl_dT[n] = stflout.T - T1
            p1, h1 = p1 - dp, h1 + dh

        stfl_h = np.flip(stfl_h); stfl_T = np.flip(stfl_T); stfl_s = np.flip(stfl_s)

        hx_valid = True
        if c == 'c':
            T_diff = np.round(wfl_T - stfl_T, 3)
        elif c == 'd':
            T_diff = np.round(stfl_T - wfl_T, 3)
        for n in range(len(T_diff)):
            if T_diff[n] < -0.001:
                hx_valid = False
                wfl_top = False
                dh_factor -= 0.025
                iterated = True
                i_time += 1
                break

    output.iteration = i_time
    mr = (wflout.h - wflin.h) / (stflin.h - stflout.h)
    
    #If iteration is not required to find a valid TQ diagram, pinch point must be at either ends
    if iterated == False:
        #wfl stream Q_max
        qa = 0
        for n in range(len(wfl_dT)):
            qa += wfl_cp[n] * wfl_dT[n]

        #stfl stream Q_max
        qb = 0
        for n in range(len(stfl_dT)):
            qb += stfl_cp[n] * stfl_dT[n]
        qb = mr*qb

        #Q_max being the smallest
        if abs(qa) < abs(qb):
            qmax = abs(qa)
            pinch_index = len(wfl_dT) - 1
        else:
            qmax = abs(qb)
            pinch_index = 0
    
    #If iteration is required to find a valid TQ diagram, pinch point is definitely somewhere in between
    elif iterated == True:
        pinch_index = (T_diff).argmin()
        qmax = 0
        for n in range(len(wfl_dT)):
            if n > pinch_index:
                qmax += mr*abs(stfl_dT[n]*stfl_cp[n])
            else:
                qmax += abs(wfl_dT[n]*wfl_cp[n])

    #Finding TQ diagram using Q_actual = e * Q_max
    qact = e * qmax

    #TQ diagram of wfl stream
    h1 = wflin.h
    p1 = wflin.p
    if c == 'c':
        wflout = State('hp', h1 - qact, wflin.p*(1-p_loss), fluid = wflin.fluid)
    elif c == 'd':
        wflout = State('hp', h1 + qact, wflin.p*(1-p_loss), fluid = wflin.fluid)
    wfl_h = np.zeros(npt); wfl_T = np.zeros(npt); wfl_s = np.zeros(npt); 
    wfl_cp = np.zeros(npt); wfl_dT = np.zeros(npt)
    wfl_dQT = 0
    dh = (wflout.h - h1) / (npt - 1)
    dp = (wflin.p * p_loss) / (npt - 1)

    wflout.update('pT', wflin.p, wflin.T)
    for n in range(npt):
        T1 = wflout.T
        wfl_dQT += dh/T1
        wflout.update('hp', h1, p1)
        wfl_h1 = h1; wfl_p1 = p1
        wfl_h[n] = wflout.h; wfl_T[n] = wflout.T; wfl_s[n] = wflout.s
        wfl_cp[n] = wflout.fluid.cpmass()
        wfl_dT[n] = wflout.T - T1
        p1, h1 = p1 - dp, h1 + dh
    
    wfl_h_norm = ((wfl_h - wfl_h[0])/(wfl_h[-1]-wfl_h[0]))*100

    #TQ diagram of stfl stream
    h1 = stflin.h
    p1 = stflin.p
    if c == 'c':
        stflout = State('hp', h1 + qact/mr, stflin.p, fluid = stflin.fluid)
    elif c == 'd':
        stflout = State('hp', h1 - qact/mr, stflin.p, fluid = stflin.fluid)
    stfl_h = np.zeros(npt); stfl_T = np.zeros(npt); stfl_s = np.zeros(npt)
    stfl_cp = np.zeros(npt); stfl_dT = np.zeros(npt)
    stfl_dQT = 0
    dh = (stflout.h - h1) / (npt - 1)
    dp = 0 #(stflin.p * p_loss) / (npt - 1)

    stflout.update('pT', stflin.p, stflin.T)
    for n in range(npt):
        T1 = stflout.T
        stfl_dQT += dh/T1
        stflout.update('hp', h1, p1)
        stfl_h1 = h1; stfl_p1 = p1
        stfl_h[n] = stflout.h; stfl_T[n] = stflout.T; stfl_s[n] = stflout.s
        stfl_cp[n] = stflout.fluid.cpmass()
        stfl_dT[n] = stflout.T - T1
        p1, h1 = p1 - dp, h1 + dh

    stfl_h = np.flip(stfl_h); stfl_T = np.flip(stfl_T); stfl_s = np.flip(stfl_s)
    stfl_h_norm = ((stfl_h - stfl_h[0])/(stfl_h[-1]-stfl_h[0]))*100

    mr = (wflout.h - wflin.h) / (stflin.h - stflout.h)
    T_diff = abs(np.round(wfl_T - stfl_T, 3))
    pinch_index = (T_diff).argmin()

    liquid_T = wflout.fluid.T_critical()

    output.mass_r = mr

    #uncomment folow section for single (non-loop mode)
    '''
    h_bal = io.Dict({'desc'                             : 'Heat Exchanger Enthalpy Balance',
    'Availability decrease in wfl side'                 : wflout.b - wflin.b,
    'Work loss due to irreversibility in wfl side'      : wflout.T0 * (wflout.s - wflin.s),
    'wfl side sum'                                      : wflout.b - wflin.b + wflout.T0 * (wflout.s - wflin.s),
    'Enthalpy decrease in wfl side'                     : wflout.h - wflin.h,
    'Availability increase in stfl side'                : mr*(stflout.b - stflin.b),
    'Work loss due to irreversibility in stfl side'     : mr*(stflout.T0 * (stflout.s - stflin.s)),
    'stfl side sum'                                     : mr*(stflout.b - stflin.b + stflout.T0 * (stflout.s - stflin.s)),
    'Enthalpy increase in stfl side'                    : mr*(stflout.h - stflin.h)})

    output.h_bal = h_bal'''
    wfl_min = State('pT', wflin.p, stflin.T, fluid = wflin.fluid)
    eta = mr*(stflout.b - stflin.b)/(wflout.b - wflin.b)
    zeta = mr*(stflout.b - stflin.b)/(wflin.b - wfl_min.b)


    output.eta = abs(eta)
    output.zeta = abs(zeta)

    State.T0 = 298.15
    stflout.update('hp', stfl_h1, stfl_p1)
    wflout.update('hp', wfl_h1, wfl_p1)

    output.wflout = wflout
    output.wflin = wflin
    output.wflmin = wfl_min
    output.stflout = stflout
    output.stflin = stflin
    output.wfldQT = wfl_dQT
    output.stfldQT = stfl_dQT

    wlirr_wfl = wflout.T0 * ((wflout.s - wflin.s))# - wfl_dQT)
    wlirr_stfl = stflout.T0 * ((stflout.s - stflin.s))# - stfl_dQT)
    wlirr = (wlirr_wfl + mr * wlirr_stfl)

    output.wlirr_per_unit_wfl_flow = wlirr
    output.wlirr_wfl = wlirr_wfl
    output.wlirr_stfl_per_unit_stfl_flow = wlirr_stfl
    output.b_delta_wfl = wflout.b - wflin.b
    output.b_delta_stfl = stflout.b - stflin.b
    output.b_delta_stfl_per_unit_wfl_flow = mr * (stflout.b - stflin.b)
    output.h_delta_stfl_per_unit_wfl_flow = mr * (stflout.h - stflin.h)
    output.b_delta = (wflout.b - wflin.b) + mr * (stflout.b - stflin.b)
    output.h_delta = (wflout.h - wflin.h) + mr * (stflout.h - stflin.h)
    
    output.wfl_T = wfl_T; output.wfl_s = wfl_s; output.wfl_h = wfl_h
    output.stfl_T = stfl_T; output.stfl_s = stfl_s; output.stfl_h = stfl_h

    '''if c == 'c':
        #plt.axhline(liquid_T, label = "Critical temperature of %s" %wflin.fluid.name)
        plt.plot(stfl_h_norm, stfl_T, c = "b", label = 'Cold fluid')
        plt.plot(wfl_h_norm, wfl_T, c = "r", label = 'Hot fluid')
        #plt.scatter([pinch_index + pinch_index*(1/(1+npt)), pinch_index + pinch_index*(1/(1+npt))], [wfl_T[pinch_index], stfl_T[pinch_index]], label = "Pinch Point")
        #for n in range(len(cold_T)):
            #plt.plot([wfl_h_norm[n],cold_h_norm[n]], [wfl_T[n], cold_T[n]])
        #plt.title("Pinch point temp diff:" + str(min(T_diff)))
        #plt.legend()
        plt.xlim(0, 100)
        plt.xlabel('Normalised heat transferred, %')
        plt.ylabel('Temperature, K')
        plt.legend()

        plt.show()'''

    #print("wfl mCp:", wfl_cp)
    #print("wfl mCp:", mr*cold_cp)
    #print(qa)
    #print(qb)

    return output