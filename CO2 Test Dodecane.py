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
from Code.HX import *
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
Oil = IFL('TX22')
Salt = IFL('NaK')

Neon = ST('pT', P0, 220, 'Air_1', NE)
Helium = ST('pT', P0, 220, 'Air_1', HE)
Dodecane = ST('pT', P0, T0, 'Air_1', C12)
Nitrogen = ST('crit', P0, 220, 'Air_1', N)
Air = ST('pT', P0, T0, 'Air_1', Air_fl)
co2_1 = ST('pT', P0, T0, 'Initial LCES', co2)
Sodium = ST('pT', P0, 400, 'Liquid Sodium', Na)
ThermOil = ST('pT', P0, T0, 'Thermal Oil', Oil)
ThermSalt = ST('pT', P0, T0 + 300, 'Thermal Salt', Salt)

print(Sodium.T0)
print(ThermOil.T0)
print(ThermSalt.T0)
print(Nitrogen)
#MoltenSalt = ST('pT', P0, T0 + 300, 'Nitrate Salt', NiSalt)
#co2_1 = hx(co2_0, T0, 0.01)


N_1 = ST('pT', P0, T0, 'Initial PTES',N)
N_2 = polytrope(N_1, 2*N_1.p, 0.9)
N_3 = HX_chg(N_2, Dodecane, 0.95, 10, 0)
N_4 = polytrope(N_3.hotout, 2*N_3.hotout.p, 0.9)
N_5 = HX_chg(N_4, Dodecane, 0.95, 10, 0)
N_6 = polytrope(N_5.hotout, 2*N_5.hotout.p, 0.9)
N_7 = HX_chg(N_6, Dodecane, 0.95, 10, 0)
N_8 = polytrope(N_7.hotout, N_1.p, 0.9)
N_9 = hx(N_8, 220, 0)


co2_2 = polytrope(co2_1, 30*co2_1.p, 0.9)
co2_3 = HX_chg(co2_2, Dodecane, 0.95, 10, 0)
co2_4 = polytrope(co2_3.hotout, 5*co2_3.hotout.p, 0.9)
co2_5 = HX_chg(co2_4, Dodecane, 0.95, 10, 0)
#co2_6 = polytrope(co2_5.hotout, 5.75*co2_5.hotout.p, 0.9)
#co2_7 = HX_chg(co2_6, Dodecane, 0.95, 10, 0)
co2_8 = coupler_new(co2_5.hotout, N_9, 0.95, 10, 0)
co2_9 = cryotb_test(co2_8.hotout, 0.9, co2_8.hotout.T0)


#print(MoltenSalt.fluid.name())

N_10 = coupler_new(co2_5.hotout, N_9, 0.95, 10, 0)
N_11 = hx(N_10.coldout, T0, 0)
#co2_8 = coupler(co2_7.hotout, might b230, Neon, 10, 0)

#coupler_new(Dodecane, Nitrogen, 0.9, 10, 0)

#coupler_new(co2_2, Dodecane, 0.95, 10, 0)
#coupler_new(co2_4, Dodecane, 0.95, 10, 0)
#print(coupler_new(co2_6, Nitrogen, 0.95, 10, 0))
#print(coupler_new(co2_6, Neon, 0.95, 10, 0))
#print(coupler_new(co2_6, Helium, 0.95, 10, 0))
#print(coupler_new(co2_6, Nitrogen, 0.95, 10, 0))
#print('NE_11', NE_11)

looping = io.Dict({'desc' : 'Parameters for looping mode'})

'''
couplerplot(co2_3)
couplerplot(co2_5)
couplerplot(co2_7)
couplerplot(co2_8)

couplerplot(NE_9)

#co2_p_2 = polytrope(co2_p_1, 12*co2_p_1.p, 0.9)
#co2_p_3 = hx(co2_p_2, T0, 0.01)
#co2_p_4 = polytrope(co2_p_3, P0, 0.9)

#co2_7_d = coupler(co2_7, co2_6.state2.T, 123, 0.01)

#co2_8 = cryotb_test(co2_7, 0.9, T0)
'''
co2_list = [co2_1, co2_2, co2_3, co2_4, co2_5, co2_8, co2_9]
N_list = [N_1, N_2, N_3, N_4, N_5, N_6, N_7, N_8, N_9, N_10, N_11]
cycle = calcs_CC(co2_list, N_list)
print(cycle)

#print(t)

'''
air = FL('Air')

air_c1 = ST('pT', P0, T0, 'Air_1', air)
air_c2 = polytrope(air_c1, 5.3*air_c1.p, 0.9)
air_c3 = hx(air_c2, T0, 0.01)
air_c4 = polytrope(air_c3, 5.3*air_c3.state2.p, 0.9)
air_c5 = hx(air_c4, T0, 0.01)
air_c6 = polytrope(air_c5, 5.3*air_c5.state2.p, 0.9)
air_c7 = hx(air_c6, T0, 0.01)
air_c8 = hx(air_c7, 120, 0.01)
air_c9 = cryotb_test(air_c8, 0.9)
air_c_list = [air_c1, air_c2, air_c3, air_c4, air_c5, air_c6, air_c7, air_c8, air_c9, "Charging for air"]

air_c = FL("Air")
air_c0 = ST('pT', 500000, T0, 'Air_1', air_c)
air_c1 = polytrope(air_c0, 1000000, 0.9)
air_c2 = hx(air_c1, T0, 0.01)
air_c3 = polytrope(air_c2, 500000, 0.9)
air_c4 = hx(air_c3, T0, 0.01)

co2_c = FL("co2")
co2_c0 = ST('pT', 500000, T0, 'co2_1', co2_c)
co2_c1 = polytrope(co2_c0, 1000000, 0.9)
co2_c2 = hx(co2_c1, T0, 0.01)
co2_c3 = polytrope(co2_c2, 500000, 0.9)
co2_c4 = hx(co2_c3, T0, 0.01)


air_d1 = hx(air_c0, air_c3.unpack('state2')[0].T+20, 0.01)
air_d2 = polytrope(air_d1, 2000000, 0.9)
air_d3 = hx(air_d2, air_c1.unpack('state2')[0].T-20, 0.01)
air_d4 = polytrope(air_d3, P0, 0.9)
air_d5 = hx(air_d4, T0, 0.01)


#print(co2_3.unpack('state2')[0].C)
#cycplot([air_c_list, co2_list], 's','t')
#cycplot([co2_list], 's','t', co2)
#print(air_c0.fluid.name())
#cycplot([state_list], 's', 't', co2)
'''
#cycplot([co2_list], 's','t')
#hxplot([co2_3])
#hxplot([co2_7_d, co2_8_d])
'''

air = FL("air")
co2 = FL("co2")
air_s_sat = air.satline().unpack('s')[0]
air_T_sat = air.satline().unpack('T')[0]
co2_s_sat = co2.satline().unpack('s')[0]
co2_T_sat = co2.satline().unpack('T')[0]

plt.plot(air_s_sat, air_T_sat-273.15, label = "Air")
plt.plot(co2_s_sat, co2_T_sat-273.15, label = "Carbon Dioxide")
plt.legend()
plt.show()'''