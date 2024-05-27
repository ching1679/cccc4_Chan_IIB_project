import CoolProp as CP

from cpstate import State as ST
from cpstate import Fluid as FL
from cpstate import IncompFluid as IFL

import numpy as np
import pandas as pd
from cptest import *
from cycleComponents import *
from propertiesCalc import *
from plot import *
from ioutils import *
from satlineAS import *
import matplotlib.pyplot as plt
from Coupler_e_Tp_inlets import *
from HX import *
from calcs_CC import *

P0 = 100000
T0 = 298.15

co2 = FL('CO2')
co = FL('CO')
NE = FL('Neon')
C12 = FL('C12')
HE = FL('Helium')
N = FL('Nitrogen')
Air_fl = FL('Air')
Na = IFL('LiqNa')
Oil = IFL('PHR')
Salt = IFL('NaK')
Water_FL = IFL('Water')

Neon = ST('pT', P0, 220, 'Air_1', NE)
Helium = ST('pT', P0, 220, 'Air_1', HE)
Dodecane = ST('pT', P0, T0, 'Air_1', C12)
Nitrogen = ST('crit', P0, 330, 'Air_1', N)
Air = ST('pT', P0, T0, 'Air_1', Air_fl)
Water = ST('pT', 200*P0, T0, 'Air_1', Water_FL)

Sodium = ST('pT', P0, 400, 'Liquid Sodium', Na)
ThermOil = ST('pT', 103000, T0, 'Thermal Oil', Oil)
HT_ThermOil = ST('pT', P0, T0 + 50, 'Thermal Oil', Oil)
ThermSalt = ST('pT', P0, T0 + 275, 'Thermal Salt', Salt)

loop_prm = io.Dict('looping inputs.dat')
print(loop_prm)
loop_output = io.Dict({'Desc': 'Looping Output'})

cyc_input = io.Dict({'desc': 'Inputs'})
cyc_input.LCES_ov_pr = 150
cyc_input.PTES_ov_pr = 4
cyc_input.PTES_CD_PR = 2
cyc_input.hx_eff = 0.95
cyc_input.wx_eff = 0.9
cyc_input.hx_dp = 0.01
cyc_input.PTES_fl = FL(loop_prm.PTES_g)

pt_n = loop_prm.p2_ln * loop_prm.p1_pt
pt_done = 0

LCES_ov_pr = np.array([])
PTES_ov_pr = np.array([])
PTES_CD_PR = np.array([])
HX_eff = np.array([])
HX_dp = np.array([])
WX_eff = np.array([])
LCES_h_eff = np.array([])
LCES_b_eff = np.array([])
LCES_h_vden = np.array([])
PTES_h_eff = np.array([])
PTES_b_eff = np.array([])
PTES_h_vden = np.array([])
CC_h_eff = np.array([])
CC_b_eff = np.array([])
CC_h_vden = np.array([])
PTES_ch_wx_loss = np.array([])
PTES_ch_hx_loss = np.array([])
PTES_ch_hrej_loss = np.array([])
LCES_ch_wx_loss = np.array([])
LCES_ch_hx_loss = np.array([])
LCES_ch_hrej_loss = np.array([])
coup_ch_loss = np.array([])
PTES_dch_wx_loss = np.array([])
PTES_dch_hx_loss = np.array([])
PTES_dch_hrej_loss = np.array([])
LCES_dch_wx_loss = np.array([])
LCES_dch_hx_loss = np.array([])
LCES_dch_hrej_loss = np.array([])
coup_dch_loss = np.array([])


for i in range(loop_prm.p2_ln):
    if loop_prm.p2 == 'PTES_P':         #PTES Charging Pressure Ratio
        cyc_input.PTES_ov_pr = loop_prm.p2_st + i * (loop_prm.p2_en - loop_prm.p2_st)/(loop_prm.p2_ln - 1)
    
    elif loop_prm.p2 == 'PTES_CD_PR':   #PTES Charging Pressure Ratio to PTES Discharging Pressure Ratio
        cyc_input.PTES_CD_PR = loop_prm.p2_st + i * (loop_prm.p2_en - loop_prm.p2_st)/(loop_prm.p2_ln - 1)

    elif loop_prm.p2 == 'LCES_P':       #LCES Overall Pressure Ratio
        cyc_input.LCES_ov_pr = loop_prm.p2_st + i * (loop_prm.p2_en - loop_prm.p2_st)/(loop_prm.p2_ln - 1)

    elif loop_prm.p2 == 'hx_eff':       #Heat Exchange Effectiveness
        cyc_input.hx_eff = loop_prm.p2_st + i * (loop_prm.p2_en - loop_prm.p2_st)/(loop_prm.p2_ln - 1)

    elif loop_prm.p2 == 'wx_eff':       #Turbomachinery Polytropic Efficiency
        cyc_input.wx_eff = loop_prm.p2_st + i * (loop_prm.p2_en - loop_prm.p2_st)/(loop_prm.p2_ln - 1)

    elif loop_prm.p2 == 'hx_dp':        #Pressure drop across heat exchange
        cyc_input.hx_dp = loop_prm.p2_st + i * (loop_prm.p2_en - loop_prm.p2_st)/(loop_prm.p2_ln - 1)

    for j in range(loop_prm.p1_pt):
            if loop_prm.p1 == 'PTES_P':
                cyc_input.PTES_ov_pr = loop_prm.p1_st + j * (loop_prm.p1_en - loop_prm.p1_st)/(loop_prm.p1_pt - 1)
    
            elif loop_prm.p1 == 'PTES_CD_PR':
                cyc_input.PTES_CD_PR = loop_prm.p1_st + j * (loop_prm.p1_en - loop_prm.p1_st)/(loop_prm.p1_pt - 1)            
                
            elif loop_prm.p1 == 'LCES_P':
                cyc_input.LCES_ov_pr = loop_prm.p1_st + j * (loop_prm.p1_en - loop_prm.p1_st)/(loop_prm.p1_pt - 1)

            elif loop_prm.p1 == 'hx_eff':
                cyc_input.hx_eff = loop_prm.p1_st + j * (loop_prm.p1_en - loop_prm.p1_st)/(loop_prm.p1_pt - 1)

            elif loop_prm.p1 == 'wx_eff':
                cyc_input.wx_eff = loop_prm.p1_st + j * (loop_prm.p1_en - loop_prm.p1_st)/(loop_prm.p1_pt - 1)

            elif loop_prm.p1 == 'hx_dp':
                cyc_input.hx_dp = loop_prm.p1_st + j * (loop_prm.p1_en - loop_prm.p1_st)/(loop_prm.p1_pt - 1)

            #print("PTES_ov_PR", cyc_input.PTES_ov_pr)
            #print("LCES_ov_PR", cyc_input.LCES_ov_pr)
            '''N_1c = ST('pT', P0, T0, 'Initial PTES', cyc_input.PTES_fl)
            N_2c = polytrope(N_1c, cyc_input.PTES_ov_pr*N_1c.p, cyc_input.wx_eff, 'c')
            N_3c = HX(N_2c.state2, ThermOil, cyc_input.hx_eff, cyc_input.hx_dp, 'c')
            N_4c = polytrope(N_3c.wflout, N_1c.p*((1 + cyc_input.hx_dp)**3), cyc_input.wx_eff, 'c')
            N_5c = ex_hx(N_4c.state2, 220, cyc_input.hx_dp, 'c')'''

            #Charging
            co2_1c = ST('pT', P0, T0, 'Initial LCES', co2)
            co2_2c = polytrope(co2_1c, (cyc_input.LCES_ov_pr**0.5)*co2_1c.p, cyc_input.wx_eff, 'c')
            co2_3c = HX(co2_2c.state2, ThermOil, cyc_input.hx_eff, cyc_input.hx_dp, 'c')
            co2_3cc = ex_hx(co2_3c.wflout, T0, cyc_input.hx_dp, 'c')
            co2_4c = polytrope(co2_3cc.state2, cyc_input.LCES_ov_pr*co2_1c.p, cyc_input.wx_eff, 'c')
            co2_5ca = HX(co2_4c.state2, HT_ThermOil, cyc_input.hx_eff, cyc_input.hx_dp, 'c')
            
            N_5c = ST('pT', (1 + 0.01)*P0, 220, 'Initial PTES', N)
            co2_6c = coupler(co2_5ca.wflout, N_5c, cyc_input.hx_eff, cyc_input.hx_dp, 'c')
            co2_7c = cryotb_test(co2_6c.LCESout, cyc_input.wx_eff, 'c')

            N_1c = coupler(co2_5ca.wflout, N_5c, cyc_input.hx_eff, cyc_input.hx_dp, 'c')
            N_2c = polytrope(N_1c.PTESout, cyc_input.PTES_ov_pr*N_1c.PTESout.p, cyc_input.wx_eff, 'c')
            N_3c = HX(N_2c, ThermOil, cyc_input.hx_eff, cyc_input.hx_dp, 'c')
            N_4c = polytrope(N_3c.wflout, N_1c.PTESout.p*(1 + cyc_input.hx_dp)**2, cyc_input.wx_eff, 'c')
            if N_4c.state2.T < 220:
                N_5c = ex_hx(N_4c.state2, 220, 0, 'c')
            else:
                N_5c = ex_hx(N_4c.state2, N_4c.state2.T, 0, 'c')

            co2_6c = coupler(co2_5ca.wflout, N_5c, cyc_input.hx_eff, cyc_input.hx_dp, 'c')
            co2_7c = cryotb_test(co2_6c.LCESout, cyc_input.wx_eff, 'c')

            #Discharging
            N_3di = ST('pT', (cyc_input.PTES_ov_pr/cyc_input.PTES_CD_PR)*P0, T0, 'Initial PTES', cyc_input.PTES_fl)
            N_2d = HX(N_3di, N_3c.stflout, cyc_input.hx_eff, cyc_input.hx_dp, 'd')
            N_1d = polytrope(N_2d.wflout, P0, cyc_input.wx_eff, 'd') 
            co2_6d = polytrope(co2_7c.state2, co2_6c.LCESout.p, cyc_input.wx_eff, 'd')
            N_5d = coupler(co2_6d.state2, N_1d, cyc_input.hx_eff, cyc_input.hx_dp, 'd')
            N_3d = polytrope(N_5d.PTESout, N_3di.p*(1 + cyc_input.hx_dp), cyc_input.wx_eff, 'd')
            N_3dd = ex_hx(N_3d.state2, T0, cyc_input.hx_dp, 'd')

            co2_5d = coupler(co2_6d.state2, N_1d.state2, cyc_input.hx_eff, cyc_input.hx_dp, 'd')
            co2_4da = HX(co2_5d.LCESout, co2_5ca.stflout, cyc_input.hx_eff, cyc_input.hx_dp, 'd')
            co2_3d = polytrope(co2_4da.wflout, co2_2c.state2.p, cyc_input.wx_eff, 'd')  
            co2_2d = HX(co2_3d.state2, co2_3c.stflout, cyc_input.hx_eff, cyc_input.hx_dp, 'd')
            co2_1d = polytrope(co2_2d.wflout, co2_1c.p*(1 + cyc_input.hx_dp), cyc_input.wx_eff, 'd')
            co2_1 = ex_hx(co2_1d.state2, T0, cyc_input.hx_dp, 'd')

            co2_list = [co2_1c, co2_2c, co2_3c, co2_3cc, co2_4c, co2_5ca, co2_6c, co2_7c, co2_6d, co2_5d, co2_4da, co2_3d, co2_2d, co2_1d, co2_1]
            N_list = [N_5c, N_1c, N_2c, N_3c, N_4c, N_5c, N_2d, N_1d, N_5d, N_3d, N_3dd]
            cycle = calcs_CC(co2_list, N_list, 'l')            





            LCES_ov_pr = np.append(LCES_ov_pr, cyc_input.LCES_ov_pr)
            PTES_ov_pr = np.append(PTES_ov_pr, cyc_input.PTES_ov_pr)
            PTES_CD_PR = np.append(PTES_CD_PR, cyc_input.PTES_CD_PR)
            HX_eff = np.append(HX_eff, cyc_input.hx_eff)
            HX_dp = np.append(HX_dp, cyc_input.hx_dp)
            WX_eff = np.append(WX_eff, cyc_input.wx_eff)
            LCES_h_eff = np.append(LCES_h_eff, cycle.h_eff_LCES)
            LCES_b_eff = np.append(LCES_b_eff, cycle.b_eff_LCES)
            LCES_h_vden = np.append(LCES_h_vden, cycle.h_vden_LCES)
            PTES_h_eff = np.append(PTES_h_eff, cycle.h_eff_PTES)
            PTES_b_eff = np.append(PTES_b_eff, cycle.b_eff_PTES)
            PTES_h_vden = np.append(PTES_h_vden, cycle.h_vden_PTES)
            CC_h_eff = np.append(CC_h_eff, cycle.h_eff_CC)
            CC_b_eff = np.append(CC_b_eff, cycle.b_eff_CC)
            CC_h_vden = np.append(CC_h_vden, cycle.h_vden_CC)
            PTES_ch_wx_loss = np.append(PTES_ch_wx_loss, cycle.PTES_ch_wx_loss)
            PTES_ch_hx_loss = np.append(PTES_ch_hx_loss, cycle.PTES_ch_hx_loss)
            PTES_ch_hrej_loss = np.append(PTES_ch_hrej_loss, cycle.PTES_ch_hrej_loss)            
            LCES_ch_wx_loss = np.append(LCES_ch_wx_loss, cycle.LCES_ch_wx_loss)
            LCES_ch_hx_loss = np.append(LCES_ch_hx_loss, cycle.LCES_ch_hx_loss)
            LCES_ch_hrej_loss = np.append(LCES_ch_hrej_loss, cycle.LCES_ch_hrej_loss)
            coup_ch_loss = np.append(coup_ch_loss, cycle.coup_ch_loss)
            PTES_dch_wx_loss = np.append(PTES_dch_wx_loss, cycle.PTES_dch_wx_loss)
            PTES_dch_hx_loss = np.append(PTES_dch_hx_loss, cycle.PTES_dch_hx_loss)
            PTES_dch_hrej_loss = np.append(PTES_dch_hrej_loss, cycle.PTES_dch_hrej_loss)            
            LCES_dch_wx_loss = np.append(LCES_dch_wx_loss, cycle.LCES_dch_wx_loss)
            LCES_dch_hx_loss = np.append(LCES_dch_hx_loss, cycle.LCES_dch_hx_loss)
            LCES_dch_hrej_loss = np.append(LCES_dch_hrej_loss, cycle.LCES_dch_hrej_loss)
            coup_dch_loss = np.append(coup_dch_loss, cycle.coup_dch_loss)

            pt_done += 1
            print(pt_done)
            #print(t2-t1,t3-t2,t4-t3,t5-t4,t6-t5,t7-t6,t8-t7,t9-t8)

results = {'LCES Overall Pressure Ratio': LCES_ov_pr,
           'PTES Overall Pressure Ratio': PTES_ov_pr,
           'PTES Charging-Discharging Pressure Ratio': PTES_CD_PR,
           'Heat Exchanger Efficiency': HX_eff,
           'Turbomachinery Efficiency': WX_eff,
           'Pressure drop across heat exchanger': HX_dp,
           'LCES Energy Efficiency': LCES_h_eff,
           'LCES Exergy Efficiency': LCES_b_eff,
           'LCES Energy Density': LCES_h_vden,
           'PTES Energy Efficiency': PTES_h_eff,
           'PTES Exergy Efficiency': PTES_b_eff,
           'PTES Energy Density': PTES_h_vden,
           'CC Energy Efficiency': CC_h_eff,
           'CC Exergy Efficiency': CC_b_eff,
           'CC Energy Density': CC_h_vden,
           "PTES Turbomachinery Losses (charging) / % ": PTES_ch_wx_loss,
           "PTES Heat Exchanger Losses (charging) / % ": PTES_ch_hx_loss,
           "PTES Heat Rejection (charging) / % ": PTES_ch_hrej_loss,
           "LCES Turbomachinery Losses (charging) / % ": LCES_ch_wx_loss,
           "LCES Heat Exchanger Losses (charging) / % ": LCES_ch_hx_loss,
           "LCES Heat Rejection (charging) / % ": LCES_ch_hrej_loss,
           "Coupler Losses (charging) / % ": coup_ch_loss,
           "PTES Turbomachinery Losses (discharging) / % ": PTES_dch_wx_loss,
           "PTES Heat Exchanger Losses (discharging) / % ": PTES_dch_hx_loss,
           "PTES Heat Rejection (discharging) / % ": PTES_dch_hrej_loss,
           "LCES Turbomachinery Losses (discharging) / % ": LCES_dch_wx_loss,
           "LCES Heat Exchanger Losses (discharging) / % ": LCES_dch_hx_loss,
           "LCES Heat Rejection (discharging) / % ": LCES_dch_hrej_loss,
           "Coupler Losses (discharging) / % ": coup_dch_loss,

           }
results_df = pd.DataFrame(data=results)
results_df.to_csv('Loop Dodecane.csv')
print(results_df)
