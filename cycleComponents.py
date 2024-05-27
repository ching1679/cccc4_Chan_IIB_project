import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import math
import numpy as np
from cpstate import State as ST
import CoolProp as CP

class state:
    def __init__(self, n, t, p, h, s, gam, cv, cp, R):
        self.n = n
        self.t = t #temperature in Kelvin
        self.p = p #pressure in Pa
        self.h = h #enthalpy
        self.s = s #entropy
        self.gam = gam
        self.cv = cv
        self.cp = cp
        self.R = R

    def __str__(self):
        return f"{self.n}({self.t})"

def poly_comp(pre_st, pr, poly_eff, mode):
    step = 100
    if mode == 'Ratio':
        dP = (pr*pre_st.p - pre_st.p)/step
    elif mode == 'Absolute':
        dP = (pr - pre_st.p)/step

    inter_st = pre_st

    for i in range(step):
        inter_h = inter_st.h
        inter_d = inter_st.deffl.rhomass()
        pred_h = inter_h + (dP/(inter_d*poly_eff))
        pred_st = ST('hp', pred_h, inter_st.p + dP)
        av_d = 0.5 * (inter_d + pred_st.deffl.rhomass())
        inter_h = inter_h + (dP/(av_d*poly_eff))
        inter_st = ST('hp', inter_h, inter_st.p + dP)

    dh = inter_st.h - pre_st.h
    dq = ((inter_st.T + pre_st.T) * (inter_st.s - pre_st.s))/2
    new_state = {
        "st"    : inter_st,
        "dh"    : dh,
        "dq"    : dq,
        "type"  : 'w in'
    }
    return new_state

def compressor(pre_st, pr, eff):
    cp = pre_st.fluid.cpmass()*1000   #Constant Pressure Heat Capacity in kJ/kg/K
    cv = pre_st.fluid.cvmass()*1000   #Cosntant Volume Heat Capacity in kJ/kg/K
    gam = cp/cv
    new_T = pre_st.T*(((pr**((gam-1)/gam)-1)/eff)+1)
    new_p = pre_st.p * pr
    label = pre_st.desc[:-1] + str(int(pre_st.desc[-1])+1)
    if eff == 1:
        new_st = ST('ps', new_p, pre_st.s, pre_st.fluid.name())
    else:
        new_st = ST('pT', new_p, new_T, pre_st.fluid.name())
    dh = (new_st.h - pre_st.h)      #Change in enthalpy in kJ/kg
    dq = ((new_st.T + pre_st.T) * (new_st.s - pre_st.s))/2      #Heat stored in fluid
    process = "w in"
    st_process = {
        "st"    : new_st,
        "dh"    : dh,
        "dq"    : dq,
        "type"  : process
    }
    return st_process

def turbine(pre_st, pr, eff):
    global w_out
    cp = pre_st.fluid.cpmass()*1000   #Constant Pressure Heat Capacity in kJ/kg/K
    cv = pre_st.fluid.cvmass()*1000   #Cosntant Volume Heat Capacity in kJ/kg/K
    gam = cp/cv
    new_T = pre_st.T*(1-eff*(1-(pr**((1-gam)/gam))))
    new_p = pre_st.p * (1/pr)
    label = pre_st.desc[:-1] + str(int(pre_st.desc[-1])+1)
    if eff == 1:
        new_st = ST('ps', new_p, pre_st.s, pre_st.fluid.name())
    else:
        new_st = ST('pT', new_p, new_T, pre_st.fluid.name())
    dh = (new_st.h - pre_st.h)      #Change in enthalpy in kJ/kg
    dq = ((new_st.T + pre_st.T) * (new_st.s - pre_st.s))/2      #Heat removed from fluid
    process = "w out"
    st_process = {
        "st"    : new_st,
        "dh"    : dh,
        "dq"    :dq,
        "type"  : process
    }
    return st_process

def heat_addition(pre_st, out_t, p_loss=None):
    global q_in
    if out_t < pre_st.t:
        print('Error in heat addition to state', pre_st.n)
        print('Outlet temperature is lower than inlet temperature')
        return False
    else:
        new_T = out_t + 273.15
        if p_loss == None:
            new_p = pre_st.p
        else:
            new_p = pre_st.p*(1-p_loss)
        new_st = ST('pT', new_p, new_T, pre_st.fluid.name())
        dh = (new_st.h - pre_st.h)      #Change in enthalpy in kJ/kg
        dq = dh      #Heat stored in fluid
        process = "q in"
        st_process = {
            "st"    : new_st,
            "dh"    : dh,
            "dq"    : dq,
            "type"  : process
        }
    return st_process

def heat_rejection(pre_st, out_t, p_loss=None):
    global q_out
    if out_t > pre_st.t:
        print('Error in heat addition to state', pre_st.n)
        print('Outlet temperature is lower than inlet temperature')
        return False
    else:
        new_T = out_t + 273.15
        if p_loss == None:
            new_p = pre_st.p
        else:
            new_p = pre_st.p*(1-p_loss)
        new_st = ST('pT', new_p, new_T, pre_st.fluid.name())
        dh = (new_st.h - pre_st.h)     #Change in enthalpy in kJ/kg
        dq = dh      #Heat removed from fluid
        process = "q out"
        st_process = {
            "st"    : new_st,
            "dh"    : dh,
            "dq"    :dq,
            "type"  : process
        }
    return st_process

def chg_cyc(list):
    w_in = 0
    w_out = 0
    q_in = 0
    q_out = 0

    clean_list = []

    for i in list:
        if type(i) == dict:
            if i['type'] == 'w in':
                w_in += i['dq']
            elif i['type'] == 'w out':
                w_out += i['dq']
            elif i['type'] == 'q in':
                q_in += i['dq']
            elif i['type'] == 'q out':
                q_out += i['dq']

    net = w_in + w_out + q_in + q_out

    energies = {
        'Tot W in'  : w_in,
        'Tot W out' : w_out,
        'Tot Q in'  : q_in,
        'Tot Q out' : q_out,
        'Net': net
    }

    return energies
