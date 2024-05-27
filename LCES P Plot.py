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
from time import process_time
import pprint

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 16}

plt.rc('font', **font)
plt.rcParams["figure.figsize"]  = [8,8]

df = pd.read_csv('config_2_LCES_100_250_PTES_3_5_PHR.csv')
df = df.transpose()

LCES_P = df.loc["LCES Overall Pressure Ratio"]
LCES_P = np.array(LCES_P[0:151])

cc_eff = df.loc["CC Energy Efficiency"]
cc_eff_PTES_3b = np.array(cc_eff[0:151])
cc_eff_PTES_4b = np.array(cc_eff[151:302])
cc_eff_PTES_5b = np.array(cc_eff[302:453])

PTES_eff = df.loc["PTES Energy Efficiency"]
PTES_eff_PTES_3b = np.array(PTES_eff[0:151])
PTES_eff_PTES_4b = np.array(PTES_eff[151:302])
PTES_eff_PTES_5b = np.array(PTES_eff[302:453])

LCES_eff = df.loc["LCES Energy Efficiency"]
LCES_eff_PTES_3b = np.array(LCES_eff[0:151])
LCES_eff_PTES_4b = np.array(LCES_eff[151:302])
LCES_eff_PTES_5b = np.array(LCES_eff[302:453])

hden = df.loc["CC Energy Density"]
hden_PTES_3b = np.array(hden[0:151])
hden_PTES_4b = np.array(hden[151:302])
hden_PTES_5b = np.array(hden[302:453])

coup_loss_dch = np.array(df.loc['Coupler Losses (discharging) / % '])
coup_loss_dch_PTES_4b = coup_loss_dch[151:302]

coup_loss_ch = np.array(df.loc['Coupler Losses (charging) / % '])
coup_loss_ch_PTES_4b = coup_loss_ch[151:302]

f1 = plt.figure(1)
plt.plot(LCES_P, cc_eff_PTES_3b*100, label = '3')
plt.plot(LCES_P, cc_eff_PTES_4b*100, label = '4')
plt.plot(LCES_P, cc_eff_PTES_5b*100, label = '5')
plt.legend(title = 'PTES Charging\nPressure Ratio')
plt.ylabel('Round Trip Combined Cycle \n Energy Efficiency, %')
plt.xlabel('LCES Overall Pressure Ratio')
plt.show()

f2 = plt.figure(2)
plt.plot(LCES_P, hden_PTES_3b/3.6e6, label = '3')
plt.plot(LCES_P, hden_PTES_4b/3.6e6, label = '4')
plt.plot(LCES_P, hden_PTES_5b/3.6e6, label = '5')
plt.legend(title = 'PTES Charging\nPressure Ratio')
plt.ylabel('Combined Cycle Energy Density, kWh/m$^{3}$')
plt.xlabel('LCES Overall Pressure Ratio')
plt.show()


f3 = plt.figure(3)
plt.plot(LCES_P, LCES_eff_PTES_3b*100, label = '3')
plt.plot(LCES_P, LCES_eff_PTES_4b*100, label = '4')
plt.plot(LCES_P, LCES_eff_PTES_5b*100, label = '5')
plt.legend(title = 'PTES Charging\nPressure Ratio')
plt.ylabel('Round Trip LCES Energy Efficiency, %')
plt.xlabel('LCES Overall Pressure Ratio')
plt.show()

f4 = plt.figure(4)
plt.plot(LCES_P, PTES_eff_PTES_3b*100, label = '3')
plt.plot(LCES_P, PTES_eff_PTES_4b*100, label = '4')
plt.plot(LCES_P, PTES_eff_PTES_5b*100, label = '5')
plt.legend(title = 'PTES Charging\nPressure Ratio')
plt.ylabel('Round Trip PTES Energy Efficiency, %')
plt.xlabel('LCES Overall Pressure Ratio')
plt.show()

plt.rcParams["figure.figsize"]  = [8,8]
f5 = plt.figure(5)
plt.plot(LCES_P, coup_loss_ch_PTES_4b*100, label = 'Charging')
plt.plot(LCES_P, coup_loss_dch_PTES_4b*100, label = 'Discharging')
plt.legend()
plt.ylabel('Lost Work in Coupler, %')
plt.xlabel('LCES Overall Pressure Ratio')
plt.show()

