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

'''
e_vol_HP_list = []
t_list = []

d_state = ST('pT', P0, T0, 'dead state', co2)

X13_co2_LP = ST('pT', P0, T0, 'X13_LP', co2)
LP_rhomass = X13_co2_LP.deffl.rhomass()
e_mass_LP = (X13_co2_LP.b - d_state.b) / 3600000    #kWh/kg
e_vol_LP = e_mass_LP * LP_rhomass   #kWh/m3

X13_co2_HP = ST('pT', 10000000, 313.15, 'X13_HP', co2)
HP_rhomass = X13_co2_HP.deffl.rhomass()
e_mass_HP = (X13_co2_HP.b - d_state.b) / 3600000  #kWh/kg
e_vol_HP = e_mass_HP * HP_rhomass   #kWh/m3


for i in range(200):
    t = 298.15 + i * 2
    t_list.append(t)
    X13_co2_HP = ST('pT', 10000000, t, 'X13_HP', co2)
    HP_rhomass = X13_co2_HP.deffl.rhomass()
    e_mass_HP = (X13_co2_HP.b - d_state.b) / 3600000  #kWh/kg
    e_vol_HP = e_mass_HP * HP_rhomass   #kWh/m3
    e_vol_HP_list.append(e_vol_HP)


V_tot = (1 / HP_rhomass) + 0.01
e_vol_eff = (e_mass_HP - e_mass_LP) / V_tot

print('Exergy density by mass:', e_mass_LP, e_mass_HP)
print('Exergy density by volume:', e_vol_LP, e_vol_HP)
print(e_vol_eff)
t_x13 = X13_co2_LP.table()

print(t_x13)
'''
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
ThermOil = ST('pT', 101000, T0, 'Thermal Oil', Oil)
HT_ThermOil = ST('pT', P0, T0 + 50, 'Thermal Oil', Oil)
ThermSalt = ST('pT', P0, T0 + 275, 'Thermal Salt', Salt)


#Charging cycle
'''N_1c = ST('pT', P0, T0, 'Initial PTES', N)
N_2c = polytrope(N_1c, 4*N_1c.p, 0.9, 'c')
N_3c = HX(N_2c, ThermOil, 0.95, 0.01, 'c')
N_4c = polytrope(N_3c.wflout, N_1c.p, 0.9, 'c')
N_5c = ex_hx(N_4c.state2, 220, 0.01, 'c')'''

co2_1c = ST('pT', P0, T0, 'Initial LCES', co2)
co2_2c = polytrope(co2_1c, 9*co2_1c.p, 0.9, 'c')
co2_3c = HX(co2_2c, ThermOil, 0.95, 0.01, 'c')
co2_3cc = ex_hx(co2_3c.wflout, T0, 0.01, 'c')
co2_4c = polytrope(co2_3cc.state2, 9*co2_3c.wflout.p, 0.9, 'c')
co2_5ca = HX(co2_4c, HT_ThermOil, 0.95, 0.01, 'c')

N_5c = ST('pT', P0, 220, 'Initial PTES', N)
co2_6c = coupler(co2_5ca.wflout, N_5c, 0.95, 0.01, 'c')
co2_7c = cryotb_test(co2_6c.LCESout, 0.9, 'c')

N_1c = coupler(co2_5ca.wflout, N_5c, 0.95, 0.01, 'c')
N_2c = polytrope(N_1c.PTESout, 4*N_1c.PTESout.p, 0.9, 'c')
N_3c = HX(N_2c, ThermOil, 0.95, 0.01, 'c')
N_4c = polytrope(N_3c.wflout, N_1c.PTESout.p*(1 + 0.01), 0.9, 'c')
if N_4c.state2.T < 220:
    N_5c = ex_hx(N_4c.state2, 220, 0, 'c')
else:
    N_5c = ex_hx(N_4c.state2, N_4c.state2.T, 0, 'c')

co2_6c = coupler(co2_5ca.wflout, N_5c, 0.95, 0.01, 'c')
co2_7c = cryotb_test(co2_6c.LCESout, 0.9, 'c')

N_1c = coupler(co2_5ca.wflout, N_5c, 0.95, 0.01, 'c')
N_2c = polytrope(N_1c.PTESout, 4*N_1c.PTESout.p, 0.9, 'c')
N_3c = HX(N_2c, ThermOil, 0.95, 0.01, 'c')
N_4c = polytrope(N_3c.wflout, N_1c.PTESout.p*(1 + 0.01), 0.9, 'c')



#Discharging cycle
N_3di = ST('pT', 2*P0, T0, 'Initial PTES', N)
N_2d = HX(N_3di, N_3c.stflout, 0.95, 0.01, 'd')
N_1d = polytrope(N_2d.wflout, P0, 0.9, 'd')
co2_6d = polytrope(co2_7c.state2, co2_6c.LCESout.p, 0.9, 'd')
N_5d = coupler(co2_6d.state2, N_1d, 0.95, 0.01, 'd')
N_3d = polytrope(N_5d.PTESout, N_3di.p, 0.9, 'd')
N_3dd = ex_hx(N_3d.state2, T0, 0, 'd')

co2_5d = coupler(co2_6d.state2, N_1d.state2, 0.95, 0.01, 'd')
co2_4da = HX(co2_5d.LCESout, co2_5ca.stflout, 0.95, 0.01, 'd')
co2_3d = polytrope(co2_4da.wflout, co2_2c.state2.p, 0.9, 'd')
co2_3dd = ex_hx(co2_3d.state2, T0, 0, 'd')
co2_2d = HX(co2_3dd.state2, co2_3c.stflout, 0.95, 0.01, 'd')
co2_1d = polytrope(co2_2d.wflout, co2_1c.p, 0.9, 'd')
co2_1 = ex_hx(co2_1d.state2, T0, 0, 'd')


looping = io.Dict({'desc' : 'Parameters for looping mode'})

co2_list = [co2_1c, co2_2c, co2_3c, co2_3cc, co2_4c, co2_5ca, co2_6c, co2_7c, co2_6d, co2_5d, co2_4da, co2_3d, co2_3dd, co2_2d, co2_1d, co2_1]
N_list = [N_5c, N_1c, N_2c, N_3c, N_4c, N_5c, N_2d, N_1d, N_5d, N_3d, N_3dd]
cycle = calcs_CC(co2_list, N_list, 's')
print(cycle)

'''#print(co2_8.LCESout)
#print(co2_9.state2)
#print(co2_10.state2)
#print(t)


air_c1 = ST('pT', P0, T0, 'Air_1', Air_fl)
air_c2 = polytrope(air_c1, 5.3*air_c1.p, 0.9, 'c')
air_c3 = ex_hx(air_c2, T0, 0.01, 'c')
air_c4 = polytrope(air_c3, 5.3*air_c3.state2.p, 0.9, 'c')
air_c5 = ex_hx(air_c4, T0, 0.01, 'c')
air_c6 = polytrope(air_c5, 5.3*air_c5.state2.p, 0.9, 'c')
air_c7 = ex_hx(air_c6, T0, 0.01, 'c')
air_c8 = ex_hx(air_c7, 120, 0.01, 'c')
#air_c9 = cryotb_test(air_c8, 0.9)
air_c_list = [air_c1, air_c2, air_c3, air_c4, air_c5, air_c6, air_c7, air_c8, "Charging for air"]

co2_c = FL("co2")
co2_c0 = ST('pT', P0, T0, 'co2_1', co2_c)
co2_c1 = polytrope(co2_c0, 5.3*co2_c0.p, 0.9, 'c')
print(co2_c1)
co2_c2 = ex_hx(co2_c1, T0, 0.01, 'c')
co2_c3 = polytrope(co2_c2, 5.3*co2_c2.state2.p, 0.9, 'c')
co2_c4 = ex_hx(co2_c3, T0, 0.01, 'c')
co2_c5 = polytrope(co2_c4, 5.3*co2_c4.state2.p, 0.9, 'c')
co2_c6 = ex_hx(co2_c5, 220, 0.01, 'c')

co2_list = [co2_c0, co2_c1, co2_c2, co2_c3, co2_c4, co2_c5, co2_c6, "Charging for CO2"]

#print(co2_3.unpack('state2')[0].C)
cycplot([air_c_list, co2_list], 's','t')
#cycplot([co2_list], 's','t', co2)
#print(air_c0.fluid.name())
#cycplot([state_list], 's', 't', co2)

#cycplot([co2_list], 's','t')
#hxplot([co2_3])
#hxplot([co2_7_d, co2_8_d])


air = FL("nitrogen")
co2 = FL("co2")
co2_s_sat = co2.satline().unpack('s')[0]
co2_T_sat = co2.satline().unpack('T')[0]
co2_s_crit = co2.crit_isobar().unpack('s')[0]
co2_T_crit = co2.crit_isobar().unpack('T')[0]

plt.plot(co2_s_sat, co2_T_sat, c = 'r', linestyle = 'solid', label = "Saturation Line for Carbon Dioxide")
plt.plot(co2_s_crit, co2_T_crit, c = 'r', linestyle = 'dashed', label = "Critical Isobar for Carbon Dioxide")
plt.ylabel('Temperature, K')
plt.xlabel('Specific entropy, Kj/K.kg')
plt.title('Saturation line and critical isobar for carbon dioxide on T-s diagram')
plt.legend()
plt.show()'''