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
plt.rcParams["figure.figsize"]  = [12,8]

species = ['LCES$_{ch}$', 'LCES$_{dch}$', 'PTES$_{ch}$', 'PTES$_{dch}$', 'Coupler$_{ch}$', 'Coupler$_{dch}$']
turbo   = np.array([5.5, 5.1, 10.9, 5.0, 0, 0])
hx      = np.array([2.4, 1.3, 3.8, 3.0, 4.6, 4.5])
hrej    = np.array([2.9, 3.3, 0, 13, 0, 0])

plt.bar(species, turbo, color = 'darkred', label = 'Turbomachinery')
plt.bar(species, hx, bottom = turbo, color = 'darkolivegreen', label = 'Heat Exchange')
plt.bar(species, hrej, bottom = turbo + hx, color = 'navy', label = 'Heat Rejection')
plt.ylabel('Lost Work, %')
plt.legend()
plt.show()
