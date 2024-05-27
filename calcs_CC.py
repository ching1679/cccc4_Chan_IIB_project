import CoolProp as CP

from cpstate import State as ST

import numpy as np
from cptest import *
from cycleComponents import *
from propertiesCalc import *
from plot import *
from ioutils import *
from satlineAS import *
import matplotlib.pyplot as plt
from Coupler_e_Tp_inlets import *

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 16}

plt.rc('font', **font)
plt.rcParams["figure.figsize"]  = [16,8]

def calcs_CC(LCES, PTES, mode):
    output = io.Dict({'desc':'Cycle Calculations'})
    initial_b = 0
    storage_per_co2 = 0
    ch_h_LCES = np.array([])
    ch_T_LCES = np.array([])
    ch_s_LCES = np.array([])
    dch_h_LCES = np.array([])
    dch_T_LCES = np.array([])
    dch_s_LCES = np.array([])
    wx_hin_LCES = np.array([])
    wx_bin_LCES = np.array([])
    wx_hout_LCES = np.array([])
    wx_bout_LCES = np.array([])
    wx_ch_wlirr_LCES = np.array([])
    wx_dch_wlirr_LCES = np.array([])
    hx_hin_LCES = np.array([])
    hx_bin_LCES = np.array([])
    hx_hout_LCES = np.array([])
    hx_bout_LCES = np.array([])
    hx_vol_LCES = np.array([])
    hx_ch_wlirr_LCES = np.array([])
    hx_dch_wlirr_LCES = np.array([])
    coup_b_LCES = 0
    LCES_ch_st_T = np.array([])
    LCES_ch_st_s = np.array([])
    LCES_dch_st_T = np.array([])
    LCES_dch_st_s = np.array([])
    LCES_wx_ct = 0
    LCES_hx_ct = 0
    LCES_ehx_ct = 0
    h_ch_rej_LCES = np.array([])
    b_ch_rej_LCES = np.array([])
    h_dch_rej_LCES = np.array([])
    b_dch_rej_LCES = np.array([])
    LCES_h_bal = io.Dict({'desc':'Enthalpy Balance for LCES'})
    LCES_b_bal = io.Dict({'desc':'Exergy Balance for LCES'})

    coupler_mr = 0
    heat_storage_LCES = np.array([])
    for i in LCES:
        if type(i) == State:
            initial_b = i.b
            initial_co2 = i
            ch_h_LCES = np.append(ch_h_LCES, i.h)
            ch_T_LCES = np.append(ch_T_LCES, i.T)
            ch_s_LCES = np.append(ch_s_LCES, i.s)
            LCES_ch_st_T = np.append(LCES_ch_st_T, i.T)
            LCES_ch_st_s = np.append(LCES_ch_st_s, i.s)
        elif type(i) == Dict:
            if i.desc == "Output from polytrope (charging)":
                LCES_wx_ct += 1
                for n in range(len(i.h)):
                    if i.T[n] != 0 and i.s[n] != 0:
                        ch_h_LCES = np.append(ch_h_LCES, i.h[n])    
                        ch_T_LCES = np.append(ch_T_LCES, i.T[n])   
                        ch_s_LCES = np.append(ch_s_LCES, i.s[n])

                wx_hin_LCES = np.append(wx_hin_LCES, i.wx)
                wx_bin_LCES = np.append(wx_bin_LCES, i.b_delta)
                LCES_ch_st_T = np.append(LCES_ch_st_T, i.state2.T)
                LCES_ch_st_s = np.append(LCES_ch_st_s, i.state2.s)
                wx_ch_wlirr_LCES = np.append(wx_ch_wlirr_LCES, i.wlirr_mass)

                k1 = 'Enthalpy input by compressor %i (charging)' %LCES_wx_ct
                k2 = 'Exergy input by compressor %i (charging)' %LCES_wx_ct
                k3 = 'Lost power by compressor %i (charging)' %LCES_wx_ct

                LCES_h_bal.add_kv(k1, i.wx)
                LCES_b_bal.add_kv(k1, i.wx)
                LCES_b_bal.add_kv(k2, i.b_delta)
                LCES_b_bal.add_kv(k3, i.wlirr_mass)
                #if i.wlirr_mass > 0:
                    #print(i)
                
            if i.desc == "Output from polytrope (discharging)":
                for n in range(len(i.h)):
                    if i.T[n] != 0 and i.s[n] != 0:
                        dch_h_LCES = np.append(dch_h_LCES, i.h[n])    
                        dch_T_LCES = np.append(dch_T_LCES, i.T[n])   
                        dch_s_LCES = np.append(dch_s_LCES, i.s[n])
                
                wx_hout_LCES = np.append(wx_hout_LCES, -i.wx)
                wx_bout_LCES = np.append(wx_bout_LCES, -i.b_delta)
                LCES_dch_st_T = np.append(LCES_dch_st_T, i.state2.T)
                LCES_dch_st_s = np.append(LCES_dch_st_s, i.state2.s)
                wx_dch_wlirr_LCES = np.append(wx_dch_wlirr_LCES, i.wlirr_mass)

                k1 = 'Enthalpy input by compressor %i (discharging)' %LCES_wx_ct
                k2 = 'Exergy input by compressor %i (discharging)' %LCES_wx_ct
                k3 = 'Lost power by compressor %i (discharging)' %LCES_wx_ct

                LCES_h_bal.add_kv(k1, i.wx)
                LCES_b_bal.add_kv(k1, i.wx)
                LCES_b_bal.add_kv(k2, i.b_delta)
                LCES_b_bal.add_kv(k3, i.wlirr_mass)
                
                LCES_wx_ct -= 1
                #if i.wlirr_mass > 0:
                    #print(i)

            if i.desc == "Output from cryo turbine":
                ch_h_LCES = np.append(ch_h_LCES, i.h)
                ch_T_LCES = np.append(ch_T_LCES, i.T)
                ch_s_LCES = np.append(ch_s_LCES, i.s)
                wx_hin_LCES = np.append(wx_hin_LCES, i.wx)
                wx_bin_LCES = np.append(wx_bin_LCES, i.b_delta)        
                LCES_ch_st_T = np.append(LCES_ch_st_T, i.state2.T)
                LCES_ch_st_s = np.append(LCES_ch_st_s, i.state2.s)
                wx_ch_wlirr_LCES = np.append(wx_ch_wlirr_LCES, i.wlirr_mass)
                #if i.wlirr_mass > 0:
                    #print(i)
                storage_co2 = i.state2
                storage_per_co2 += storage_co2.b - initial_co2.b
                k1 = 'Enthalpy output by cryo turbine'
                k2 = 'Exergy output by cryo turbine'
                k3 = 'Lost power by cryo turbine'
                k4 = 'Storage-Ambient enthalpy difference'
                k5 = 'Storage-Ambient exergy difference'
                LCES_h_bal.add_kv(k1, i.wx)
                LCES_b_bal.add_kv(k1, i.wx)
                LCES_b_bal.add_kv(k2, i.b_delta)
                LCES_b_bal.add_kv(k3, i.wlirr_mass)
                LCES_h_bal.add_kv(k4, storage_co2.h - initial_co2.h)
                LCES_b_bal.add_kv(k5, storage_co2.b - initial_co2.b)

                LCES_wx_ct += 1
                vden_co2 = 1/storage_co2.d

            if i.desc == 'Output from HX (charging)':
                LCES_hx_ct += 1
                ch_h_LCES = np.append(ch_h_LCES, i.wfl_h)
                ch_T_LCES = np.append(ch_T_LCES, i.wfl_T)
                ch_s_LCES = np.append(ch_s_LCES, i.wfl_s)

                hx_hin_LCES = np.append(hx_hin_LCES, i.h_delta_stfl_per_unit_wfl_flow)
                hx_bin_LCES = np.append(hx_bin_LCES, i.b_delta_stfl_per_unit_wfl_flow)
                #hx_vol_LCES = np.append(hx_vol_LCES, i.mass_r/i.stflin.d + i.mass_r/i.stflout.d)   #account for both hot and cold tanks
                LCES_ch_st_T = np.append(LCES_ch_st_T, i.wflout.T)
                LCES_ch_st_s = np.append(LCES_ch_st_s, i.wflout.s)
                hx_ch_wlirr_LCES = np.append(hx_ch_wlirr_LCES, i.wlirr_per_unit_wfl_flow)
                storage_per_co2 += (i.stflout.b - i.stflin.b) * i.mass_r
                heat_storage_LCES = np.append(heat_storage_LCES, (i.stflout.h - i.stflin.h) * i.mass_r)

                k1 = 'Enthalpy rejected to heat storage %i (charging)' %LCES_hx_ct
                k2 = 'Exergy rejected by working fluid side of heat storage %i (charging)' %LCES_hx_ct
                k3 = 'Exergy increase in storage fluid side of heat storage %i (charging)' %LCES_hx_ct
                k4 = 'Lost power in heat exchanger %i (charging)' %LCES_hx_ct

                LCES_h_bal.add_kv(k1, i.wflout.h - i.wflin.h)
                LCES_b_bal.add_kv(k2, abs(i.b_delta_wfl))
                LCES_b_bal.add_kv(k3, i.b_delta_stfl_per_unit_wfl_flow)
                #LCES_b_bal.add_kv(k1, i.h_delta)
                LCES_b_bal.add_kv(k4, i.wlirr_per_unit_wfl_flow)


            if i.desc == 'Output from HX (discharging)':
                dch_h_LCES = np.append(dch_h_LCES, i.wfl_h)
                dch_T_LCES = np.append(dch_T_LCES, i.wfl_T)
                dch_s_LCES = np.append(dch_s_LCES, i.wfl_s)

                hx_hout_LCES = np.append(hx_hout_LCES, i.h_delta_stfl_per_unit_wfl_flow)
                hx_bout_LCES = np.append(hx_bout_LCES, i.b_delta_stfl_per_unit_wfl_flow)
                hx_vol_LCES = np.append(hx_vol_LCES, i.mass_r/i.stflin.d + i.mass_r/i.stflout.d)   #account for both hot and cold tanks
                LCES_dch_st_T = np.append(LCES_dch_st_T, i.wflout.T)
                LCES_dch_st_s = np.append(LCES_dch_st_s, i.wflout.s)
                hx_dch_wlirr_LCES = np.append(hx_dch_wlirr_LCES, i.wlirr_per_unit_wfl_flow)
                storage_per_co2 += (i.stflout.b - i.stflin.b) * i.mass_r
                heat_storage_LCES = np.append(heat_storage_LCES, (i.stflout.h - i.stflin.h) * i.mass_r)

                k1 = 'Enthalpy rejected to working fluid %i (discharging)' %LCES_hx_ct
                k2 = 'Exergy increase in working fluid side of heat storage %i (discharging)' %LCES_hx_ct
                k3 = 'Exergy rejected by storage fluid side of heat storage %i (discharging)' %LCES_hx_ct
                k4 = 'Lost power in heat exchanger %i (discharging)' %LCES_hx_ct

                LCES_h_bal.add_kv(k1, i.wflout.h - i.wflin.h)
                LCES_b_bal.add_kv(k2, abs(i.b_delta_wfl))
                LCES_b_bal.add_kv(k3, abs(i.b_delta_stfl_per_unit_wfl_flow))
                #LCES_b_bal.add_kv(k1, i.h_delta)
                LCES_b_bal.add_kv(k4, i.wlirr_per_unit_wfl_flow)

                LCES_hx_ct -= 1

            if i.desc == 'Output from set temp HX (charging)':
                LCES_ehx_ct += 1
                ch_h_LCES = np.append(ch_h_LCES, i.h)
                ch_T_LCES = np.append(ch_T_LCES, i.T)
                ch_s_LCES = np.append(ch_s_LCES, i.s)
                
                h_ch_rej_LCES = np.append(h_ch_rej_LCES, abs(i.state2.h - i.state1.h))
                b_ch_rej_LCES = np.append(b_ch_rej_LCES, abs(i.state2.b - i.state1.b))
                LCES_ch_st_T = np.append(LCES_ch_st_T, i.state2.T)
                LCES_ch_st_s = np.append(LCES_ch_st_s, i.state2.s)

                k1 = 'Enthalpy gained/rejected by LCES in external heat exchanger %i' %LCES_ehx_ct
                k2 = 'Enthalpy gained/rejected by LCES in external heat exchanger %i' %LCES_ehx_ct

                LCES_h_bal.add_kv(k1, (i.state2.h - i.state1.h))
                LCES_b_bal.add_kv(k1, (i.state2.b - i.state1.b))
            
            if i.desc == 'Output from set temp HX (discharging)':
                LCES_ehx_ct += 1
                dch_h_LCES = np.append(dch_h_LCES, i.h)
                dch_T_LCES = np.append(dch_T_LCES, i.T)
                dch_s_LCES = np.append(dch_s_LCES, i.s)
                
                h_dch_rej_LCES = np.append(h_dch_rej_LCES, abs(i.state2.h - i.state1.h))
                b_dch_rej_LCES = np.append(b_dch_rej_LCES, abs(i.state2.b - i.state1.b))
                LCES_dch_st_T = np.append(LCES_dch_st_T, i.state2.T)
                LCES_dch_st_s = np.append(LCES_dch_st_s, i.state2.s)

                k1 = 'Enthalpy gained/rejected by LCES in external heat exchanger %i' %LCES_ehx_ct
                k2 = 'Enthalpy gained/rejected by LCES in external heat exchanger %i' %LCES_ehx_ct

                LCES_h_bal.add_kv(k1, (i.state2.h - i.state1.h))
                LCES_b_bal.add_kv(k1, (i.state2.b - i.state1.b))

            if i.desc == 'Output from coupler (charging)':
                ch_h_LCES = np.append(ch_h_LCES, i.LCES_h)
                ch_T_LCES = np.append(ch_T_LCES, i.LCES_T)
                ch_s_LCES = np.append(ch_s_LCES, i.LCES_s)
                ch_coupler_mr = i.mass_r
                coup_b_LCES = i.LCESout.b - i.LCESin.b
                coup_h_LCES = i.LCESout.h - i.LCESin.h
                LCES_ch_st_T = np.append(LCES_ch_st_T, i.LCESout.T)
                LCES_ch_st_s = np.append(LCES_ch_st_s, i.LCESout.s)
                coup_ch_wlirr_LCES = i.wlirr_per_unit_LCES_flow

                k1 = 'Enthalpy rejected by PTES in coupler (charging)'
                k2 = 'Exergy increase in LCES side in coupler (charging)'
                k3 = 'Exergy change by PTES side in coupler (charging)'
                k4 = 'Lost power in coupler (charging)'

                LCES_h_bal.add_kv(k1, i.LCESout.h - i.LCESin.h)
                LCES_b_bal.add_kv(k3, abs(i.b_delta_PTES_per_unit_LCES_flow))
                LCES_b_bal.add_kv(k2, i.b_delta_LCES)
                #LCES_b_bal.add_kv(k1, (i.LCESout.h - i.LCESin.h))
                LCES_b_bal.add_kv(k4, i.wlirr_per_unit_LCES_flow)

            if i.desc == 'Output from coupler (discharging)':
                dch_h_LCES = np.append(dch_h_LCES, i.LCES_h)
                dch_T_LCES = np.append(dch_T_LCES, i.LCES_T)
                dch_s_LCES = np.append(dch_s_LCES, i.LCES_s)
                dch_coupler_mr = i.mass_r
                coup_b_LCES = i.LCESout.b - i.LCESin.b
                coup_h_LCES = i.LCESout.h - i.LCESin.h
                LCES_dch_st_T = np.append(LCES_dch_st_T, i.LCESout.T)
                LCES_dch_st_s = np.append(LCES_dch_st_s, i.LCESout.s)
                coup_dch_wlirr_LCES = i.wlirr_per_unit_LCES_flow

                k1 = 'Enthalpy rejected by LCES in coupler (discharging)'
                k2 = 'Exergy increase in PTES side in coupler (discharging)'
                k3 = 'Exergy change by LCES side in coupler (discharging)'
                k4 = 'Lost power in coupler (discharging)'

                LCES_h_bal.add_kv(k1, i.LCESout.h - i.LCESin.h)
                LCES_b_bal.add_kv(k3, abs(i.b_delta_PTES_per_unit_LCES_flow))
                LCES_b_bal.add_kv(k2, i.b_delta_LCES)
                #LCES_b_bal.add_kv(k1, (i.LCESout.h - i.LCESin.h))
                
                LCES_b_bal.add_kv(k4, i.wlirr_per_unit_LCES_flow)

    ch_h_PTES = np.array([])
    ch_T_PTES = np.array([])
    ch_s_PTES = np.array([])
    dch_h_PTES = np.array([])
    dch_T_PTES = np.array([])
    dch_s_PTES = np.array([])
    wx_hin_PTES = np.array([])
    wx_bin_PTES = np.array([])
    wx_hout_PTES = np.array([])
    wx_bout_PTES = np.array([])
    wx_ch_wlirr_PTES = np.array([])
    wx_dch_wlirr_PTES = np.array([])
    hx_hin_PTES = np.array([])
    hx_bin_PTES = np.array([])
    hx_hout_PTES = np.array([])
    hx_bout_PTES = np.array([])
    hx_vol_PTES = np.array([])
    hx_ch_wlirr_PTES = np.array([])
    hx_dch_wlirr_PTES = np.array([])
    h_ch_rej_PTES = np.array([])
    b_ch_rej_PTES = np.array([])
    h_dch_rej_PTES = np.array([])
    b_dch_rej_PTES = np.array([])
    heat_storage_PTES = np.array([])
    coup_b_PTES = 0
    PTES_ch_st_T = np.array([])
    PTES_ch_st_s = np.array([])
    PTES_dch_st_T = np.array([])
    PTES_dch_st_s = np.array([])
    PTES_wx_ct = 0
    PTES_hx_ct = 0
    PTES_ehx_ct = 0
    PTES_h_bal = io.Dict({'desc':'Enthalpy Balance for PTES'})
    PTES_b_bal = io.Dict({'desc':'Exergy Balance for PTES'})

    for i in PTES:
        if type(i) == State:
            ch_h_PTES = np.append(ch_h_PTES, i.h)
            ch_T_PTES = np.append(ch_T_PTES, i.T)
            ch_s_PTES = np.append(ch_s_PTES, i.s)
            PTES_ch_st_T = np.append(PTES_ch_st_T, i.T)
            PTES_ch_st_s = np.append(PTES_ch_st_s, i.s)
        elif type(i) == Dict:
            if i.desc == "Output from polytrope (charging)":
                PTES_wx_ct += 1
                for n in range(len(i.h)):
                    if i.T[n] != 0 and i.s[n] != 0:
                        ch_h_PTES = np.append(ch_h_PTES, i.h[n])    
                        ch_T_PTES = np.append(ch_T_PTES, i.T[n])   
                        ch_s_PTES = np.append(ch_s_PTES, i.s[n])
                wx_hin_PTES = np.append(wx_hin_PTES, ch_coupler_mr * i.wx)
                wx_bin_PTES = np.append(wx_bin_PTES, ch_coupler_mr * i.b_delta)
                PTES_ch_st_T = np.append(PTES_ch_st_T, i.state2.T)
                PTES_ch_st_s = np.append(PTES_ch_st_s, i.state2.s)
                wx_ch_wlirr_PTES = np.append(wx_ch_wlirr_PTES, ch_coupler_mr * i.wlirr_mass)

                k1 = 'Enthalpy input by compressor %i (charging)' %PTES_wx_ct
                k2 = 'Exergy input by compressor %i (charging)' %PTES_wx_ct
                k3 = 'Lost power by compressor %i (charging)' %PTES_wx_ct
                PTES_h_bal.add_kv(k1, ch_coupler_mr * i.wx)
                PTES_b_bal.add_kv(k2, ch_coupler_mr * i.b_delta)
                PTES_b_bal.add_kv(k1, ch_coupler_mr * i.wx)
                PTES_b_bal.add_kv(k3, ch_coupler_mr * i.wlirr_mass)

            if i.desc == "Output from polytrope (discharging)":
                for n in range(len(i.h)):
                    if i.T[n] != 0 and i.s[n] != 0:
                        dch_h_PTES = np.append(dch_h_PTES, i.h[n])    
                        dch_T_PTES = np.append(dch_T_PTES, i.T[n])   
                        dch_s_PTES = np.append(dch_s_PTES, i.s[n])

                wx_hout_PTES = np.append(wx_hout_PTES, -dch_coupler_mr * i.wx)
                wx_bout_PTES = np.append(wx_bout_PTES, -dch_coupler_mr * i.b_delta)
                PTES_dch_st_T = np.append(PTES_dch_st_T, i.state2.T)
                PTES_dch_st_s = np.append(PTES_dch_st_s, i.state2.s)
                wx_dch_wlirr_PTES = np.append(wx_dch_wlirr_PTES, dch_coupler_mr * i.wlirr_mass)

                k1 = 'Enthalpy input by compressor %i (discharging)' %PTES_wx_ct
                k2 = 'Exergy input by compressor %i (discharging)' %PTES_wx_ct
                k3 = 'Lost power by compressor %i (discharging)' %PTES_wx_ct
                PTES_h_bal.add_kv(k1, dch_coupler_mr * i.wx)
                PTES_b_bal.add_kv(k2, dch_coupler_mr * i.b_delta)
                PTES_b_bal.add_kv(k1, dch_coupler_mr * i.wx)
                PTES_b_bal.add_kv(k3, dch_coupler_mr * i.wlirr_mass)
                PTES_wx_ct -= 1

            if i.desc == 'Output from HX (charging)':
                PTES_hx_ct += 1
                ch_h_PTES = np.append(ch_h_PTES, i.wfl_h)
                ch_T_PTES = np.append(ch_T_PTES, i.wfl_T)
                ch_s_PTES = np.append(ch_s_PTES, i.wfl_s)

                hx_hin_PTES = np.append(hx_hin_PTES, ch_coupler_mr * i.h_delta_stfl_per_unit_wfl_flow)
                hx_bin_PTES = np.append(hx_bin_PTES, ch_coupler_mr * i.b_delta_stfl_per_unit_wfl_flow)
                #hx_vol_PTES = np.append(hx_vol_PTES, coupler_mr*i.mass_r/i.stflin.d + coupler_mr*i.mass_r/i.stflout.d)   #account for both hot and cold tanks
                PTES_ch_st_T = np.append(PTES_ch_st_T, i.wflout.T)
                PTES_ch_st_s = np.append(PTES_ch_st_s, i.wflout.s)
                storage_per_co2 += ch_coupler_mr * (i.stflout.b - i.stflin.b) * i.mass_r
                heat_storage_PTES = np.append(heat_storage_PTES, ch_coupler_mr * (i.stflout.h - i.stflin.h) * i.mass_r)
                hx_ch_wlirr_PTES = np.append(hx_ch_wlirr_PTES, ch_coupler_mr * i.wlirr_per_unit_wfl_flow)

                k1 = 'Enthalpy rejected to heat storage %i' %PTES_hx_ct
                k2 = 'Exergy rejected by working fluid side of heat storage %i' %PTES_hx_ct
                k3 = 'Exergy increase in storage fluid side of heat storage %i' %PTES_hx_ct
                k4 = 'Lost power in heat exchanger %i' %PTES_hx_ct

                PTES_h_bal.add_kv(k1, ch_coupler_mr * (i.wflout.h - i.wflin.h))
                PTES_b_bal.add_kv(k2, ch_coupler_mr * abs(i.b_delta_wfl))
                PTES_b_bal.add_kv(k3, ch_coupler_mr * i.b_delta_stfl_per_unit_wfl_flow)
                PTES_b_bal.add_kv(k4, ch_coupler_mr * i.wlirr_per_unit_wfl_flow)

            if i.desc == 'Output from HX (discharging)':
                dch_h_PTES = np.append(dch_h_PTES, i.wfl_h)
                dch_T_PTES = np.append(dch_T_PTES, i.wfl_T)
                dch_s_PTES = np.append(dch_s_PTES, i.wfl_s)

                hx_hout_PTES = np.append(hx_hout_PTES, dch_coupler_mr * i.h_delta_stfl_per_unit_wfl_flow)
                hx_bout_PTES = np.append(hx_bout_PTES, dch_coupler_mr * i.b_delta_stfl_per_unit_wfl_flow)
                hx_vol_PTES = np.append(hx_vol_PTES, dch_coupler_mr*i.mass_r/i.stflin.d + dch_coupler_mr*i.mass_r/i.stflout.d)   #account for both hot and cold tanks
                PTES_dch_st_T = np.append(PTES_dch_st_T, i.wflout.T)
                PTES_dch_st_s = np.append(PTES_dch_st_s, i.wflout.s)
                storage_per_co2 += dch_coupler_mr * (i.stflout.b - i.stflin.b) * i.mass_r
                heat_storage_PTES = np.append(heat_storage_PTES, dch_coupler_mr * (i.stflout.h - i.stflin.h) * i.mass_r)
                hx_dch_wlirr_PTES = np.append(hx_dch_wlirr_PTES, dch_coupler_mr * i.wlirr_per_unit_wfl_flow)

                k1 = 'Enthalpy rejected to heat storage %i' %PTES_hx_ct
                k2 = 'Exergy rejected by working fluid side of heat storage %i' %PTES_hx_ct
                k3 = 'Exergy increase in storage fluid side of heat storage %i' %PTES_hx_ct
                k4 = 'Lost power in heat exchanger %i' %PTES_hx_ct

                PTES_h_bal.add_kv(k1, dch_coupler_mr * (i.wflout.h - i.wflin.h))
                PTES_b_bal.add_kv(k2, dch_coupler_mr * abs(i.b_delta_wfl))
                PTES_b_bal.add_kv(k3, dch_coupler_mr * i.b_delta_stfl_per_unit_wfl_flow)
                PTES_b_bal.add_kv(k4, dch_coupler_mr * i.wlirr_per_unit_wfl_flow)

                PTES_hx_ct -= 1

            if i.desc == 'Output from set temp HX (charging)':
                PTES_ehx_ct += 1
                ch_h_PTES = np.append(ch_h_PTES, i.h)
                ch_T_PTES = np.append(ch_T_PTES, i.T)
                ch_s_PTES = np.append(ch_s_PTES, i.s)
                
                h_ch_rej_PTES = np.append(h_ch_rej_PTES, abs(ch_coupler_mr * (i.state2.h - i.state1.h)))
                b_ch_rej_PTES = np.append(b_ch_rej_PTES, ch_coupler_mr * abs(i.state2.b - i.state1.b))
                PTES_ch_st_T = np.append(PTES_ch_st_T, i.state2.T)
                PTES_ch_st_s = np.append(PTES_ch_st_s, i.state2.s)

                k1 = 'Enthalpy gained/rejected by PTES in external heat exchanger %i' %PTES_ehx_ct
                k2 = 'Enthalpy gained/rejected by PTES in external heat exchanger %i' %PTES_ehx_ct

                PTES_h_bal.add_kv(k1, ch_coupler_mr * (i.state2.h - i.state1.h))
                PTES_b_bal.add_kv(k1, ch_coupler_mr * (i.state2.b - i.state1.b))

            if i.desc == 'Output from set temp HX (discharging)':
                PTES_ehx_ct += 1
                dch_h_PTES = np.append(dch_h_PTES, i.h)
                dch_T_PTES = np.append(dch_T_PTES, i.T)
                dch_s_PTES = np.append(dch_s_PTES, i.s)
                
                h_dch_rej_PTES = np.append(h_dch_rej_PTES, abs(dch_coupler_mr * (i.state2.h - i.state1.h)))
                b_dch_rej_PTES = np.append(b_dch_rej_PTES, dch_coupler_mr * abs(i.state2.b - i.state1.b))
                PTES_dch_st_T = np.append(PTES_dch_st_T, i.state2.T)
                PTES_dch_st_s = np.append(PTES_dch_st_s, i.state2.s)

                k1 = 'Enthalpy gained/rejected by PTES in external heat exchanger %i' %PTES_ehx_ct
                k2 = 'Enthalpy gained/rejected by PTES in external heat exchanger %i' %PTES_ehx_ct

                PTES_h_bal.add_kv(k1, dch_coupler_mr * (i.state2.h - i.state1.h))
                PTES_b_bal.add_kv(k1, dch_coupler_mr * (i.state2.b - i.state1.b))

            if i.desc == 'Output from coupler (charging)':
                PTES_h = np.flip(i.PTES_h)
                PTES_T = np.flip(i.PTES_T)
                PTES_s = np.flip(i.PTES_s)
                ch_h_PTES = np.append(ch_h_PTES, PTES_h)
                ch_T_PTES = np.append(ch_T_PTES, PTES_T)
                ch_s_PTES = np.append(ch_s_PTES, PTES_s)
                ch_coupler_mr = i.mass_r
                coup_b_PTES = i.b_delta_PTES_per_unit_LCES_flow
                coup_h_PTES = ch_coupler_mr*(i.PTESout.h - i.PTESin.h)
                PTES_ch_st_T = np.append(PTES_ch_st_T, i.PTESout.T)
                PTES_ch_st_s = np.append(PTES_ch_st_s, i.PTESout.s)
                coup_wlirr_PTES = ch_coupler_mr*(i.wlirr_PTES_per_unit_PTES_flow)

                k1 = 'Enthalpy gained by PTES in coupler (charging)'
                k2 = 'Exergy increase in LCES side in coupler (charging)'
                k3 = 'Exergy rejected by PTES side in coupler (charging)'
                k4 = 'Lost power in coupler (charging)'

                PTES_h_bal.add_kv(k1, i.h_delta_PTES_per_unit_LCES_flow)
                PTES_b_bal.add_kv(k3, abs(i.b_delta_PTES_per_unit_LCES_flow))
                PTES_b_bal.add_kv(k2, i.b_delta_LCES)
                #PTES_b_bal.add_kv(k1, (i.LCESout.h - i.LCESin.h))
                PTES_b_bal.add_kv(k4, i.wlirr_per_unit_LCES_flow)
            
            if i.desc == 'Output from coupler (discharging)':
                PTES_h = np.flip(i.PTES_h)
                PTES_T = np.flip(i.PTES_T)
                PTES_s = np.flip(i.PTES_s)
                dch_h_PTES = np.append(dch_h_PTES, PTES_h)
                dch_T_PTES = np.append(dch_T_PTES, PTES_T)
                dch_s_PTES = np.append(dch_s_PTES, PTES_s)
                dch_coupler_mr = i.mass_r
                coup_b_PTES = i.b_delta_PTES_per_unit_LCES_flow
                coup_h_PTES = dch_coupler_mr*(i.PTESout.h - i.PTESin.h)
                PTES_dch_st_T = np.append(PTES_dch_st_T, i.PTESout.T)
                PTES_dch_st_s = np.append(PTES_dch_st_s, i.PTESout.s)
                coup_wlirr_PTES = dch_coupler_mr*(i.wlirr_PTES_per_unit_PTES_flow)

                k1 = 'Enthalpy gained by LCES in coupler (discharging)'
                k2 = 'Exergy increase in PTES side in coupler (discharging)'
                k3 = 'Exergy rejected by LCES side in coupler (discharging)'
                k4 = 'Lost power in coupler (discharging)'

                PTES_h_bal.add_kv(k1, i.h_delta_PTES_per_unit_LCES_flow)
                PTES_b_bal.add_kv(k3, abs(i.b_delta_PTES_per_unit_LCES_flow))
                PTES_b_bal.add_kv(k2, i.b_delta_LCES)
                #PTES_b_bal.add_kv(k1, (i.LCESout.h - i.LCESin.h))
                PTES_b_bal.add_kv(k4, i.wlirr_per_unit_LCES_flow)
    
    PTES_fl = PTES[0].state2.fluid
    PTES_sat = PTES_fl.satline()
    #PTES_crit = PTES_fl.crit_isobar()
    PTES_s_sat = PTES_sat.unpack('s')[0]
    PTES_T_sat = PTES_sat.unpack('T')[0]
    #PTES_s_crit = PTES_crit.unpack('s')[0]

    co2 = FL('CO2')
    co2_sat = co2.satline()
    #co2_crit = co2.crit_isobar()
    co2_s_sat = co2_sat.unpack('s')[0]
    co2_T_sat = co2_sat.unpack('T')[0]
    #co2_s_crit = co2_crit.unpack('s')[0]
    #co2_T_crit = co2_crit.unpack('T')[0]

    losses = io.Dict({'desc': 'Exergetic losses analysis'})
    losses.add_kv("PTES Turbomachinery Losses (charging) / J/kgCO2: ", sum(wx_ch_wlirr_PTES))
    losses.add_kv("PTES Turbomachinery Losses (charging) / % ", sum(wx_ch_wlirr_PTES)/(sum(wx_bin_LCES)+sum(wx_bin_PTES)))
    losses.add_kv("PTES Heat Exchanger Losses (charging) / J/kgCO2: ", sum(hx_ch_wlirr_PTES))
    losses.add_kv("PTES Heat Exchanger Losses (charging) / % ", sum(hx_ch_wlirr_PTES)/(sum(wx_bin_LCES)+sum(wx_bin_PTES)))
    losses.add_kv("PTES Heat Rejection (charging) / J/kgCO2: ", sum(h_ch_rej_PTES))
    losses.add_kv("PTES Heat Rejection (charging) / % ", sum(h_ch_rej_PTES)/(sum(wx_bin_LCES)+sum(wx_bin_PTES)))
    losses.add_kv("LCES Turbomachinery Losses (charging) / J/kgCO2: ", sum(wx_ch_wlirr_LCES))
    losses.add_kv("LCES Turbomachinery Losses (charging) / % ", sum(wx_ch_wlirr_LCES)/(sum(wx_bin_LCES)+sum(wx_bin_PTES)))
    losses.add_kv("LCES Heat Exchanger Losses (charging) / J/kgCO2: ", sum(hx_ch_wlirr_LCES))
    losses.add_kv("LCES Heat Exchanger Losses (charging) / % ", sum(hx_ch_wlirr_LCES)/(sum(wx_bin_LCES)+sum(wx_bin_PTES)))
    losses.add_kv("LCES Heat Rejection (charging) / J/kgCO2: ", sum(h_ch_rej_LCES))
    losses.add_kv("LCES Heat Rejection (charging) / % ", sum(h_ch_rej_LCES)/(sum(wx_bin_LCES)+sum(wx_bin_PTES)))
    losses.add_kv("Coupler Losses (charging) / J/kgCO2: ", abs(coup_ch_wlirr_LCES))
    losses.add_kv("Coupler Losses (charging) / % ", abs(coup_ch_wlirr_LCES)/(sum(wx_bin_LCES)+sum(wx_bin_PTES)))
    
    losses.add_kv("PTES Turbomachinery Losses (discharging) / J/kgCO2: ", sum(wx_dch_wlirr_PTES))
    losses.add_kv("PTES Turbomachinery Losses (discharging) / % ", sum(wx_dch_wlirr_PTES)/(sum(wx_bin_LCES)+sum(wx_bin_PTES)))
    losses.add_kv("PTES Heat Exchanger Losses (discharging) / J/kgCO2: ", sum(hx_dch_wlirr_PTES))
    losses.add_kv("PTES Heat Exchanger Losses (discharging) / % ", sum(hx_dch_wlirr_PTES)/(sum(wx_bin_LCES)+sum(wx_bin_PTES)))
    losses.add_kv("PTES Heat Rejection (discharging) / J/kgCO2: ", sum(h_dch_rej_PTES))
    losses.add_kv("PTES Heat Rejection (discharging) / % ", sum(h_dch_rej_PTES)/(sum(wx_bin_LCES)+sum(wx_bin_PTES)))
    losses.add_kv("LCES Turbomachinery Losses (discharging) / J/kgCO2: ", sum(wx_dch_wlirr_LCES))
    losses.add_kv("LCES Turbomachinery Losses (discharging) / % ", sum(wx_dch_wlirr_LCES)/(sum(wx_bin_LCES)+sum(wx_bin_PTES)))
    losses.add_kv("LCES Heat Exchanger Losses (discharging) / J/kgCO2: ", sum(hx_dch_wlirr_LCES))
    losses.add_kv("LCES Heat Exchanger Losses (discharging) / % ", sum(hx_dch_wlirr_LCES)/(sum(wx_bin_LCES)+sum(wx_bin_PTES)))
    losses.add_kv("LCES Heat Rejection (discharging) / J/kgCO2: ", sum(h_dch_rej_LCES))
    losses.add_kv("LCES Heat Rejection (discharging) / % ", sum(h_dch_rej_LCES)/(sum(wx_bin_LCES)+sum(wx_bin_PTES)))
    losses.add_kv("Coupler Losses (discharging) / J/kgCO2: ", abs(coup_dch_wlirr_LCES))
    losses.add_kv("Coupler Losses (discharging) / % ", abs(coup_dch_wlirr_LCES)/(sum(wx_bin_LCES)+sum(wx_bin_PTES)))

    if mode == 's':
        plt.plot(ch_s_PTES - 5500, ch_T_PTES, label = 'Charging', c = 'indianred')
        plt.plot(dch_s_PTES - 5500, dch_T_PTES, label = 'Discharging', c = 'tomato', linestyle = 'dotted')
        plt.scatter(PTES_ch_st_s - 5500, PTES_ch_st_T, c = 'darkred', marker = '^')
        plt.scatter(PTES_dch_st_s - 5500, PTES_dch_st_T, c = 'red', marker = 'v')
        #plt.plot(PTES_s_sat - 5500, PTES_T_sat - 273.15, label = 'Saturation line for PTES fluid, %s' %PTES_fl.name())

    
        plt.plot(ch_s_LCES, ch_T_LCES, label = 'Charging', c = 'blue')
        plt.plot(dch_s_LCES, dch_T_LCES, label = 'Discharging', c = 'cornflowerblue', linestyle = 'dotted')
        plt.plot(co2_s_sat, co2_T_sat, c = 'black', linewidth = 0.5, label = 'Saturation line')
        plt.scatter(LCES_ch_st_s, LCES_ch_st_T, c = 'midnightblue', marker = '^')
        plt.scatter(LCES_dch_st_s, LCES_dch_st_T, c = 'royalblue', marker = 'v')
        plt.legend()
        plt.ylabel('Temperature, K')
        plt.xlabel('Specific Entropy, J/Kg.K')
        plt.ylim(200, 600)
        plt.show()
    
    #all values below are already referred to per unit flow of CO2 in LCES

    b_in_LCES = np.sum(wx_bin_LCES)
    b_in_PTES = np.sum(wx_bin_PTES)
    h_in_LCES = np.sum(wx_hin_LCES)
    h_in_PTES = np.sum(wx_hin_PTES)

    b_out_LCES = np.sum(wx_bout_LCES)
    b_out_PTES = np.sum(wx_bout_PTES)
    h_out_LCES = np.sum(wx_hout_LCES)
    h_out_PTES = np.sum(wx_hout_PTES)

    h_vden_co2 = (h_out_LCES + h_out_PTES)/(1/storage_co2.d)
    h_vden_LCES = h_out_LCES/(1/storage_co2.d + np.sum(hx_vol_LCES))
    h_vden_PTES = h_out_PTES/np.sum(hx_vol_PTES)
    h_vden_CC = (h_out_LCES + h_out_PTES)/(1/storage_co2.d + np.sum(hx_vol_LCES) + np.sum(hx_vol_PTES))

    b_vden_co2 = abs(storage_co2.b - initial_co2.b)/(1/storage_co2.d)
    b_vden_LCES = b_out_PTES/(1/storage_co2.d + np.sum(hx_vol_LCES))
    b_vden_PTES = b_out_PTES/np.sum(hx_vol_PTES)
    b_vden_CC = (b_out_LCES + b_out_PTES)/(1/storage_co2.d + np.sum(hx_vol_LCES) + np.sum(hx_vol_PTES))

    b_eff_LCES = b_out_LCES/b_in_LCES
    b_eff_PTES = b_out_PTES/b_in_PTES
    b_eff_CC = (b_out_LCES + b_out_PTES)/(b_in_PTES + b_in_LCES)
    #b_eff_test = storage_per_co2/(b_in_LCES + b_in_PTES)

    h_eff_LCES = h_out_LCES/h_in_LCES
    h_eff_PTES = h_out_PTES/h_in_PTES
    h_eff_CC = (h_out_LCES + h_out_PTES)/(h_in_LCES + h_in_PTES)

    output.wx_hin_LCES = wx_hin_LCES
    output.wx_bin_LCES = wx_bin_LCES
    output.wx_hout_LCES = wx_hout_LCES
    output.wx_bout_LCES = wx_bout_LCES
    output.wx_ch_wlirr_LCES = wx_ch_wlirr_LCES
    output.wx_dch_wlirr_LCES = wx_dch_wlirr_LCES
    output.hx_hin_LCES = hx_hin_LCES
    output.hx_bin_LCES = hx_bin_LCES
    output.hx_hout_LCES = hx_hout_LCES
    output.hx_bout_LCES = hx_bout_LCES
    output.hx_vol_LCES = hx_vol_LCES
    output.hx_ch_wlirr_LCES = hx_ch_wlirr_LCES
    output.hx_dch_wlirr_LCES = hx_dch_wlirr_LCES
    output.LCES_enthalpy_diff = storage_co2.h - initial_co2.h
    output.LCES_exergy_diff = storage_co2.b - initial_co2.b
    output.heat_storage_LCES = heat_storage_LCES

    output.coup_ch_wlirr_LCES = coup_ch_wlirr_LCES
    output.co2_storage = storage_co2
    output.co2_ambient = initial_co2

    output.LCES_h_bal = LCES_h_bal
    output.LCES_b_bal = LCES_b_bal
    output.PTES_h_bal = PTES_h_bal
    output.PTES_b_bal = PTES_b_bal
    
    output.PTES_ch_wx_loss = sum(wx_ch_wlirr_PTES)/(sum(wx_bin_LCES)+sum(wx_bin_PTES))
    output.PTES_ch_hx_loss = sum(hx_ch_wlirr_PTES)/(sum(wx_bin_LCES)+sum(wx_bin_PTES))
    output.PTES_ch_hrej_loss = sum(b_ch_rej_PTES)/(sum(wx_bin_LCES)+sum(wx_bin_PTES))
    output.LCES_ch_wx_loss = sum(wx_ch_wlirr_LCES)/(sum(wx_bin_LCES)+sum(wx_bin_PTES))
    output.LCES_ch_hx_loss = sum(hx_ch_wlirr_LCES)/(sum(wx_bin_LCES)+sum(wx_bin_PTES))
    output.LCES_ch_hrej_loss = sum(b_ch_rej_LCES)/(sum(wx_bin_LCES)+sum(wx_bin_PTES))
    output.coup_ch_loss = abs(coup_ch_wlirr_LCES)/(sum(wx_bin_LCES)+sum(wx_bin_PTES))
    output.PTES_dch_wx_loss = sum(wx_dch_wlirr_PTES)/(sum(wx_bin_LCES)+sum(wx_bin_PTES))
    output.PTES_dch_hx_loss = sum(hx_dch_wlirr_PTES)/(sum(wx_bin_LCES)+sum(wx_bin_PTES))
    output.PTES_dch_hrej_loss = sum(b_dch_rej_PTES)/(sum(wx_bin_LCES)+sum(wx_bin_PTES))
    output.LCES_dch_wx_loss = sum(wx_dch_wlirr_LCES)/(sum(wx_bin_LCES)+sum(wx_bin_PTES))
    output.LCES_dch_hx_loss = sum(hx_dch_wlirr_LCES)/(sum(wx_bin_LCES)+sum(wx_bin_PTES))
    output.LCES_dch_hrej_loss = sum(b_dch_rej_LCES)/(sum(wx_bin_LCES)+sum(wx_bin_PTES))
    output.coup_dch_loss = abs(coup_dch_wlirr_LCES)/(sum(wx_bin_LCES)+sum(wx_bin_PTES))
    
    output.losses = losses



    output.wx_hin_PTES = wx_hin_PTES
    output.wx_bin_PTES = wx_bin_PTES
    output.wx_hout_PTES = wx_hout_PTES
    output.wx_bout_PTES = wx_bout_PTES
    output.wx_ch_wlirr_PTES = wx_ch_wlirr_PTES
    output.wx_dch_wlirr_PTES = wx_dch_wlirr_PTES
    output.hx_hin_PTES = hx_hin_PTES
    output.hx_bin_PTES = hx_bin_PTES
    output.hx_hout_PTES = hx_hout_PTES
    output.hx_bout_PTES = hx_bout_PTES
    output.hx_vol_PTES = hx_vol_PTES
    output.hx_ch_wlirr_PTES = hx_ch_wlirr_PTES
    output.hx_dch_wlirr_PTES = hx_dch_wlirr_PTES
    output.heat_storage_PTES = heat_storage_PTES
    output.h_ch_rej_PTES = h_ch_rej_PTES
    output.b_ch_rej_PTES = b_ch_rej_PTES
    output.h_ch_rej_LCES = h_ch_rej_LCES
    output.b_ch_rej_LCES = b_ch_rej_LCES
    output.h_dch_rej_PTES = h_dch_rej_PTES
    output.b_dch_rej_PTES = b_dch_rej_PTES
    output.h_dch_rej_LCES = h_dch_rej_LCES
    output.b_dch_rej_LCES = b_dch_rej_LCES
    output.coup_b_LCES = coup_b_LCES
    output.coup_h_LCES = coup_h_LCES
    output.coup_b_PTES = coup_b_PTES
    output.coup_h_PTES = coup_h_PTES
    output.coup_wlirr_PTES = coup_wlirr_PTES
    output.coup_mr = dch_coupler_mr

    output.h_eff_LCES = h_eff_LCES
    output.h_eff_PTES = h_eff_PTES
    output.h_eff_CC = h_eff_CC
    output.b_eff_LCES = b_eff_LCES
    output.b_eff_PTES = b_eff_PTES
    output.b_eff_CC = b_eff_CC
    output.h_vden_co2 = h_vden_co2
    output.h_vden_LCES = h_vden_LCES
    output.h_vden_PTES = h_vden_PTES
    output.h_vden_CC = h_vden_CC
    output.b_vden_co2 = b_vden_co2
    output.b_vden_LCES = b_vden_LCES
    output.b_vden_PTES = b_vden_PTES
    output.b_vden_CC = b_vden_CC

    return output

    print("storage_per_co2", storage_per_co2)
    print("h_LCES", h_LCES)
    print("T_LCES", T_LCES)
    print("s_LCES", s_LCES)

    print("h_PTES", h_PTES)
    print("T_PTES", T_PTES)
    print("s_PTES", s_PTES)



