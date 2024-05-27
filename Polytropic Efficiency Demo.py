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
plt.rcParams["figure.figsize"]  = [8,8]

def comp_isen(rp):
    comp_ieff = (rp**(31/131)-1)/(rp**(31/(131*0.9))-1)
    return comp_ieff

def tb_isen(rp):
    tb_ieff = (1-((1/rp)**(0.9*(31/131))))/(1-((1/rp)**(31/131)))
    return tb_ieff

rp = np.arange(1.00001, 21, 0.1)
comp_ieff = np.array([])
tb_ieff = np.array([])

for i in rp:
    comp_ieff = np.append(comp_ieff, comp_isen(i))
    tb_ieff = np.append(tb_ieff, tb_isen(i))


rt = comp_ieff * tb_ieff

plt.legend()
plt.plot(rp, comp_ieff, label = 'Compressor')
plt.plot(rp, tb_ieff, label = 'Turbine')
plt.plot(rp, rt)
plt.ylabel('Isentropic Efficiency')
plt.xlabel('Pressure Ratio')
plt.legend()
plt.show()



