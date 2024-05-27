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

PHR = pd.read_csv('config_2_LCES_100_250_PTES_3_5_PHR.csv')
PHR = PHR.transpose()

LCES_P = PHR.loc['LCES Overall Pressure Ratio']
LCES_P = np.array(LCES_P)[151:302]

PHR_cc_eff = PHR.loc['CC Energy Efficiency']
PHR_cc_eff = np.array(PHR_cc_eff[151:302])

PHR_cc_hden = PHR.loc['CC Energy Density']
PHR_cc_hden = np.array(PHR_cc_hden[151:302])

S800 = pd.read_csv('config_2_LCES_100_250_PTES_3_5_S800.csv')
S800 = S800.transpose()

S800_cc_eff = S800.loc['CC Energy Efficiency']
S800_cc_eff = np.array(S800_cc_eff[151:302])

S800_cc_hden = S800.loc['CC Energy Density']
S800_cc_hden = np.array(S800_cc_hden[151:302])

T66 = pd.read_csv('config_2_LCES_100_250_PTES_3_5_T66.csv')
T66 = T66.transpose()

T66_cc_eff = T66.loc['CC Energy Efficiency']
T66_cc_eff = np.array(T66_cc_eff[151:302])

T66_cc_hden = T66.loc['CC Energy Density']
T66_cc_hden = np.array(T66_cc_hden[151:302])

f1 = plt.figure(1)
plt.plot(LCES_P, PHR_cc_eff*100, label = 'Paratherm HR')
plt.plot(LCES_P, S800_cc_eff*100, label = 'Syltherm 800')
plt.plot(LCES_P, T66_cc_eff*100, label = 'Therminol66')
plt.legend(title = 'Thermo Storage Medium')
plt.ylabel('Round Trip Combined Cycle \n Energy Efficiency, %')
plt.xlabel('LCES Overall Pressure Ratio')
plt.show()

f1 = plt.figure(1)
plt.plot(LCES_P, PHR_cc_hden/3.6e6, label = 'Paratherm HR')
plt.plot(LCES_P, S800_cc_hden/3.6e6, label = 'Syltherm 800')
plt.plot(LCES_P, T66_cc_hden/3.6e6, label = 'Therminol66')
plt.legend(title = 'Thermo Storage Medium')
plt.ylabel('Combined Cycle Energy Density, kWh/m$^{3}$')
plt.xlabel('LCES Overall Pressure Ratio')
plt.show()
