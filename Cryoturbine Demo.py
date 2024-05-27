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
        'size'   : 24}

plt.rc('font', **font)
plt.rcParams["figure.figsize"]  = [16,16]

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

co2_0 = ST('pT', P0, T0, 'co2', co2)


pT = 'pT'
co2_sat_T = co2.satline().unpack('T')[0]
co2_sat_s = co2.satline().unpack('s')[0]
co2_1 = polytrope(co2_0, 150*co2_0.p, 1, 'c')
co2_2 = ex_hx(co2_1.state2, 220, 0, 'c')
co2_3 = ex_hx(co2_2.state2, 240, 0, 'c')
co2_4 = cryotb_test(co2_3, 0.78)
co2_5 = cryotb_test(co2_3, 0.82)
co2_6 = cryotb_test(co2_3, 0.86)
co2_7 = cryotb_test(co2_3, 0.9)
co2_8 = cryotb_test(co2_3, 0.94)

plt.plot(co2_sat_s, co2_sat_T, c = 'black', linestyle = 'dotted', label = 'Saturation Line')
plt.plot(co2_2.s, co2_2.T, c = 'black', linestyle = 'dashed', label = "Supercritical Isobar")
plt.plot(co2_4.s, co2_4.T, linestyle = 'solid', label = '$\eta$ = 0.78')
plt.plot(co2_5.s, co2_5.T, linestyle = 'solid', label = '$\eta$ = 0.82')
plt.plot(co2_6.s, co2_6.T, linestyle = 'solid', label = '$\eta$ = 0.86')
plt.plot(co2_7.s, co2_7.T, linestyle = 'solid', label = '$\eta$ = 0.90')
plt.plot(co2_8.s, co2_8.T, linestyle = 'solid', label = '$\eta$ = 0.94')
plt.xlabel('Specific Entropy, J/kg.K')
plt.ylabel('Temperature, K')
plt.ylim(200, 550)
plt.xlim(500, 2300)
plt.legend()
plt.show()
