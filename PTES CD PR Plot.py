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
plt.rcParams["figure.figsize"]  = [9,8]

df = pd.read_csv('config_2_PTES_CD_PR_3_1-2_PTES_3_5_PHR.csv')
df = df.transpose()

PTES_CD_PR = df.loc['PTES Charging-Discharging Pressure Ratio']
PTES_CD_PR = np.array(PTES_CD_PR[100:200])

cc_eff = df.loc['CC Energy Efficiency']
cc_eff = np.array(cc_eff[100:200])

cc_hden = df.loc['CC Energy Density']
cc_hden = np.array(cc_hden[100:200])

coup_loss_dch = df.loc['Coupler Losses (discharging) / % ']
coup_loss_dch = np.array(coup_loss_dch[100:200])

hrej_loss_dch = df.loc['PTES Heat Rejection (discharging) / % ']
hrej_loss_dch = np.array(hrej_loss_dch[100:200])


f1 = plt.figure(1)
plt.plot(PTES_CD_PR, cc_eff*100)
plt.ylabel('Round Trip Combined Cycle \n Energy Efficiency, %')
plt.xlabel('Ratio of PTES Charging-Discharging Pressure Ratio')
plt.show()

f2 = plt.figure(2)
plt.plot(PTES_CD_PR, cc_hden/3.6e6)
plt.ylabel('Combined Cycle Energy Density, kWh/m$^{3}$')
plt.xlabel('Ratio of PTES Charging-Discharging Pressure Ratio')
plt.show()

f3 = plt.figure(3)
plt.plot(PTES_CD_PR, coup_loss_dch*100)
plt.ylabel('Lost Work in Coupler (Discharging), %')
plt.xlabel('Ratio of PTES Charging-Discharging Pressure Ratio')
plt.show()

f3 = plt.figure(3)
plt.plot(PTES_CD_PR, hrej_loss_dch*100)
plt.ylabel('Lost Work in PTES Heat Rejection (Discharging), %')
plt.xlabel('Ratio of PTES Charging-Discharging Pressure Ratio')
plt.show()

