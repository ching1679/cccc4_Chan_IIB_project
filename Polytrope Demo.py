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

P0 = 100000

#T0 = 298.15


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
ThermOil = ST('pT', P0, T0, 'Thermal Oil', Oil)
HT_ThermOil = ST('pT', P0, T0 + 100, 'Thermal Oil', Oil)
ThermSalt = ST('pT', P0, T0 + 275, 'Thermal Salt', Salt)

co2 = ST('pT', P0, T0, 'co2', co2)


pT = 'pT'

co2_1 = polytrope(co2, 10*co2.p, 0.78, 'c')
co2_2 = polytrope(co2, 10*co2.p, 0.82, 'c')
co2_3 = polytrope(co2, 10*co2.p, 0.86, 'c')
co2_4 = polytrope(co2, 10*co2.p, 0.9, 'c')
co2_5 = polytrope(co2, 10*co2.p, 0.94, 'c')

co2_lp = ex_hx(co2, 450, 0, 'c')
co2_lp = ex_hx(co2_lp.state2, 280, 0, 'c')

co2_hp = ex_hx(co2_5.state2, 550, 0, 'c')
co2_hp = ex_hx(co2_hp.state2, 450, 0, 'c')

np.insert(co2_1.s, 0, co2_1.state1.s)

plt.plot(np.insert(co2_1.s, 0, co2_1.state1.s), np.insert(co2_1.T, 0, co2_1.state1.T), label = '$\eta$ = 0.78')
plt.plot(np.insert(co2_2.s, 0, co2_2.state1.s), np.insert(co2_2.T, 0, co2_2.state1.T), label = '$\eta$ = 0.82')
plt.plot(np.insert(co2_3.s, 0, co2_3.state1.s), np.insert(co2_3.T, 0, co2_3.state1.T), label = '$\eta$ = 0.86')
plt.plot(np.insert(co2_4.s, 0, co2_4.state1.s), np.insert(co2_4.T, 0, co2_4.state1.T), label = '$\eta$ = 0.90')
plt.plot(np.insert(co2_5.s, 0, co2_5.state1.s), np.insert(co2_5.T, 0, co2_5.state1.T), label = '$\eta$ = 0.94')
plt.plot(co2_lp.s, co2_lp.T, c = 'black', linestyle = 'solid', label = '1 bar isobar')
plt.plot(co2_hp.s, co2_hp.T, c = 'black', linestyle = 'dashed', label = '10 bar isobar')
plt.xlabel('Specific Entropy, J/kg.K')
plt.ylabel('Temperature, K')
plt.legend()
plt.show()



print (co2_1)

