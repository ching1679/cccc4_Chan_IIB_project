import CoolProp as CP

from cpstate import State as ST
from cpstate import Fluid as FL
from cpstate import IncompFluid as IFL

import numpy as np
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
plt.rcParams["figure.figsize"]  = [16,8]

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
Water_FL = FL('Water')

Neon = ST('pT', P0, 220, 'Air_1', NE)
Helium = ST('pT', P0, 220, 'Air_1', HE)
Dodecane = ST('pT', P0, T0, 'Air_1', C12)
Nitrogen = ST('crit', P0, 330, 'Air_1', N)
Air = ST('pT', P0, 1400, 'Air_1', Air_fl)
Water = ST('pT', 10*P0, T0, 'Air_1', Water_FL)

Sodium = ST('pT', P0, 400, 'Liquid Sodium', Na)
ThermOil = ST('pT', P0, T0, 'Thermal Oil', Oil)
HT_ThermOil = ST('pT', P0, T0 + 100, 'Thermal Oil', Oil)
ThermSalt = ST('pT', P0, T0 + 275, 'Thermal Salt', Salt)

cp_80 = np.array([])
cp_120 = np.array([])
cp_160 = np.array([])


co2_80 = ST('pT', 80*P0, 500, 'co2', co2)
co2_120 = polytrope(co2_80, 120*P0, 1, 'c')
co2_120 = ex_hx(co2_120.state2, 500, 0, 'c').state2
co2_160 = polytrope(co2_80, 160*P0, 1, 'c')
co2_160 = ex_hx(co2_160.state2, 500, 0, 'c').state2

co2_crit = polytrope(co2_80, co2.p_critical(), 1, 'c').state2
co2_crit.update('pT', co2.p_critical(), 500)

N_1 = ST('pT', P0, 220, 'N', N)

co2_crit_hx = coupler(co2_crit, N_1, 0.95, 0, 'c')
co2_160_hx = coupler(co2_160, N_1, 0.95, 0, 'c')

print(co2_crit_hx)
print(co2_160_hx)

final_T = 220
for i in range(280):
    co2_80.update('pT', 80*P0, co2_80.T - 1)
    cp_80 = np.append(cp_80, co2_80.fluid.cpmass())
    co2_120.update('pT', 120*P0, co2_120.T - 1)
    cp_120 = np.append(cp_120, co2_120.fluid.cpmass())
    co2_160.update('pT', 160*P0, co2_160.T - 1)
    cp_160 = np.append(cp_160, co2_160.fluid.cpmass())

T = np.linspace(500, 220, 280)
plt.plot(T, cp_80, label = '80 bar')
plt.plot(T, cp_120, label = '120 bar')
plt.plot(T, cp_160, label = '160 bar')
plt.ylim(0, 10000)
plt.xlim(200, 500)
plt.axvline(co2.T_critical(), c = 'grey', linestyle = 'dotted')
plt.savefig("test.png", bbox_inches='tight')
plt.annotate('Critical Temperature of CO$_{2}$', (co2.T_critical() + 2, 100))
plt.legend()
plt.ylabel('Specific heat capacity, J/Kg.K')
plt.xlabel('Temeprature, K')
plt.show()



