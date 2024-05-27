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
from matplotlib.ticker import FormatStrFormatter
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

df = pd.read_csv('config_2_hx_eff_85_95_LCES_100_200_PHR.csv')
df = df.transpose()

hx_eff = df.loc["Heat Exchanger Efficiency"]
hx_eff = np.array(hx_eff[50:100])

cc_eff = df.loc["CC Energy Efficiency"]
cc_eff_LCES_120 = np.array(cc_eff[50:100])
cc_eff_LCES_160 = np.array(cc_eff[150:200])
cc_eff_LCES_200 = np.array(cc_eff[250:300])

PTES_eff = df.loc["PTES Energy Efficiency"]
PTES_eff_LCES_120 = np.array(PTES_eff[50:100])
PTES_eff_LCES_160 = np.array(PTES_eff[150:200])
PTES_eff_LCES_200 = np.array(PTES_eff[250:300])

LCES_eff = df.loc["LCES Energy Efficiency"]
LCES_eff_LCES_120 = np.array(LCES_eff[50:100])
LCES_eff_LCES_160 = np.array(LCES_eff[150:200])
LCES_eff_LCES_200 = np.array(LCES_eff[250:300])

hden = df.loc["CC Energy Density"]
hden_LCES_120 = np.array(hden[50:100])
hden_LCES_160 = np.array(hden[150:200])
hden_LCES_200 = np.array(hden[250:300])

PTES_hx_loss_dch = np.array(df.loc['PTES Heat Exchanger Losses (discharging) / % '])    #losses taken at LCES 160
PTES_hx_loss_ch = np.array(df.loc['PTES Heat Exchanger Losses (charging) / % '])
PTES_hx_loss = PTES_hx_loss_dch + PTES_hx_loss_ch
PTES_hx_loss = PTES_hx_loss[150:200]

LCES_hx_loss_dch = np.array(df.loc['LCES Heat Exchanger Losses (discharging) / % '])
LCES_hx_loss_ch = np.array(df.loc['LCES Heat Exchanger Losses (charging) / % '])
LCES_hx_loss = LCES_hx_loss_dch + LCES_hx_loss_ch
LCES_hx_loss = LCES_hx_loss[150:200]

coup_loss_dch = np.array(df.loc['Coupler Losses (discharging) / % '])
coup_loss_dch = coup_loss_dch[150:200]

coup_loss_ch = np.array(df.loc['Coupler Losses (charging) / % '])
coup_loss_ch = coup_loss_ch[150:200]
print(coup_loss_dch)

f1 = plt.figure(1)
plt.plot(hx_eff, cc_eff_LCES_120*100, label = '120')
plt.plot(hx_eff, cc_eff_LCES_160*100, label = '160')
plt.plot(hx_eff, cc_eff_LCES_200*100, label = '200')
plt.legend(title = 'LCES Charging\nPressure Ratio')
plt.ylabel('Round Trip Combined Cycle \n Energy Efficiency, %')
plt.xlabel('Heat Exchanger Effectiveness')
plt.show()

f2 = plt.figure(2)
plt.plot(hx_eff, hden_LCES_120/3.6e6, label = '120')
plt.plot(hx_eff, hden_LCES_160/3.6e6, label = '160')
plt.plot(hx_eff, hden_LCES_200/3.6e6, label = '200')
plt.legend(title = 'LCES Charging\nPressure Ratio')
plt.ylabel('Combined Cycle Energy Density, kWh/m$^{3}$')
plt.xlabel('Heat Exchanger Effectiveness')
plt.show()


f3 = plt.figure(3)
plt.plot(hx_eff, LCES_eff_LCES_120*100, label = '120')
plt.plot(hx_eff, LCES_eff_LCES_160*100, label = '160')
plt.plot(hx_eff, LCES_eff_LCES_200*100, label = '200')
plt.legend(title = 'LCES Charging\nPressure Ratio')
plt.ylabel('Round Trip LCES \n Energy Efficiency, %')
plt.xlabel('Heat Exchanger Effectiveness')
plt.show()

f4 = plt.figure(4)
plt.plot(hx_eff, PTES_eff_LCES_120*100, label = '120')
plt.plot(hx_eff, PTES_eff_LCES_160*100, label = '160')
plt.plot(hx_eff, PTES_eff_LCES_200*100, label = '200')
plt.legend(title = 'LCES Charging\nPressure Ratio')
plt.ylabel('Round Trip PTES \n Energy Efficiency, %')
plt.xlabel('Heat Exchanger Effectiveness')
plt.show()


f5 = plt.figure(5)
plt.plot(hx_eff, PTES_hx_loss*100, label = 'PTES')
plt.plot(hx_eff, LCES_hx_loss*100, label = 'LCES')
plt.legend()
plt.ylabel('Round Trip Lost Work in Heat Exchanger\n with Thermal Storage Medium, %')
plt.xlabel('Heat Exchanger Effectiveness')
plt.show()

f6 = plt.figure(6)
plt.plot(hx_eff, coup_loss_ch*100, label = 'Charging')
plt.plot(hx_eff, coup_loss_dch*100, label = 'Discharging')
plt.legend()
plt.ylabel('Lost Work in Coupler, %')
plt.xlabel('Heat Exchanger Effectiveness')
plt.show()
