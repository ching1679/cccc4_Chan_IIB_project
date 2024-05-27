import math
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from scipy.optimize import curve_fit
import numpy as np
import CoolProp as CP
from cpstate import State
from cpstate import Fluid as FL
from cycleComponents import *
from satlineAS import *
from ioutils import *
import random

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 16}

plt.rc('font', **font)
plt.rcParams["figure.figsize"]  = [18,9]

def cycplot(list, x, y):

    plt.xlabel("Specific entropy, kJ/K.kg")
    plt.ylabel("Temperature, C")
    plt.axhline(y = 20, linewidth = 1, linestyle = '--', label = "Ambient Temperature")

    #ax.set_xticks([])
    #ax.set_yticks([])
    
    for l in list:
        x_list = np.array([])
        y_list = np.array([])
        legend = None

        dead_st = l[0]
        if type(l[-1]) == str:
            storage_st = l[-2]
            legend = l[-1]
        else:
            storage_st = l[-1]

        if type(dead_st) == Dict:
            dead_st = dead_st.unpack('state2')[0]
        if type(storage_st) == Dict:
            storage_st = storage_st.unpack('state2')[0]

        col = (np.random.random(), np.random.random(), np.random.random())
        fluid_FL = dead_st.fluid
        sat = fluid_FL.satline()
        crit = fluid_FL.crit_isobar()

        s_sat = sat.unpack('s')[0]
        T_sat = sat.unpack('T')[0]
        s_crit = crit.unpack('s')[0]
        T_crit = crit.unpack('T')[0]

        plt.plot(s_sat, T_sat-273.15, c = col, linewidth = 2, linestyle = '--', label = 'Saturation Line for ' + fluid_FL.name())
        #plt.plot(s_crit, T_crit-273.15, linewidth = 2, label = 'Critical Line ' + fluid_FL.name())
        #exergy = storage_st.b - dead_st.b
        #exergy_density = exergy * storage_st.d

        if legend != None:
            for i in l[:-1]:
                if type(i) == State:
                    x_list = np.append(x_list, i.s)
                    print(i.s)
                    y_list = np.append(y_list, i.T)
                elif type(i) == Dict:
                    for j in i.unpack('s')[0]:
                        x_list = np.append(x_list, j)
                    for k in i.unpack('T')[0]:
                        y_list = np.append(y_list, k)
        else:
            for i in l:
                if type(i) == State:
                    x_list = np.append(x_list, i.s)
                    print(i.s)
                    y_list = np.append(y_list, i.T)
                elif type(i) == Dict:
                    for j in i.unpack('s')[0]:
                        x_list = np.append(x_list, j)
                    for k in i.unpack('T')[0]:
                        y_list = np.append(y_list, k)

        greyness = np.random.random() 
        plt.plot(x_list, y_list-273.15, linewidth = 2, c = col, label = legend)        



    '''
    if x == 'h':
        for i in clean_list:
            x_list.append(i.h/1000)
        plt.xlabel('Enthalpy (kJ/kg/K)', fontsize = 15)
    elif x == 's':
        for i in clean_list:
            x_list.append(i.s/1000)
        plt.xlabel('Entropy (kJ/kg/K)', fontsize = 15)
    elif x == 't':
        for i in clean_list:
            x_list.append(i.t)
        plt.xlabel('Temperature'+ u'\u2103', fontsize = 15)

    if y == 'h':
        for i in clean_list:
            y_list.append(i.h/1000)
        plt.ylabel('Enthalpy (kJ/kg/K)', fontsize = 15)
    elif y == 's':
        for i in clean_list:
            y_list.append(i.s/1000)
        plt.ylabel('Entropy (kJ/kg/K)', fontsize = 15)
    elif y == 't':
        for i in clean_list:
            y_list.append(i.t)
        plt.ylabel('Temperature'+ u'\u2103', fontsize = 15)'''
    


    '''
    ax.xaxis.set_major_locator(MultipleLocator(5))
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.yaxis.set_major_locator(MultipleLocator(10))
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax.grid(which = 'major', color='dimgrey', linestyle='-', linewidth=0.25, alpha=0.5)
    ax.grid(which = 'minor', color='lightgrey', linestyle='-', linewidth=0.25, alpha=0.5)
    plt.rc('legend', fontsize=12) 
    plt.rc('xtick', labelsize=12) 
    plt.rc('ytick', labelsize=12)
    plt.rc('axes', labelsize=15)
    plt.rc('figure', titlesize=15)'''
    #plt.plot(s,Ti, label = 'Inversion')
 
    #textstr = '\n'.join((
    #'Working Fluid: %s' % (fluid.name()),
    #'Exergy Density: %s kWh/m3' % int(exergy_density/3600000)))

    # these are matplotlib.patch.Patch properties
    #props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

    # place a text box in upper left in axes coords
    #ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=10,
    #    verticalalignment='top', bbox=props)

    plt.title('T-s diagram for typical gas liquefaction')
    plt.legend(loc = 'upper left')
    plt.show()

def hxplot(list):
    fig, ax = plt.subplots()
    x_list = np.array([])
    y_list = np.array([])


    for i in list:
        if type(i) == State:
            x_list = np.append(x_list, i.h)
            y_list = np.append(y_list, i.T)
        elif type(i) == Dict:
            for j in i.unpack('h')[0]:
                x_list = np.append(x_list, j)
            for k in i.unpack('T')[0]:
                y_list = np.append(y_list, k)    
    min = x_list[0]
    max = x_list[-1]
    x_list = ((x_list - min)/(max-min))*100
    plt.xlabel("Percentage of heat transferred")
    plt.ylabel("Temperature, C")
    plt.title("T-Q diagram for heat exchanger 3")

    plt.plot(x_list, y_list-273.15, linewidth = 1, c = 'green')
    plt.show()

def couplerplot(dict):
    fig, ax = plt.subplots()
    hot_x_list = np.array(dict.hot_h)
    cold_x_list = np.array(dict.cold_h)
    hot_y_list = np.array(dict.hot_T)
    cold_y_list = np.array(dict.cold_T)

    min = hot_x_list[0]
    max = hot_x_list[-1]
    hot_x_list = ((hot_x_list - min)/(max-min))*100
    plt.plot(hot_x_list, hot_y_list-273.15, linewidth = 1, c = 'red', label = dict.hotout.fluid.name())

    min = cold_x_list[0]
    max = cold_x_list[-1]
    cold_x_list = ((cold_x_list - min)/(max-min))*100
    #cold_x_list = np.flip(cold_x_list)
    plt.plot(cold_x_list, cold_y_list-273.15, linewidth = 1, c = 'blue', label = dict.coldout.fluid.name())
    pinch = str((hot_y_list - cold_y_list).min())
    plt.legend()
    plt.xlabel("Percentage of heat transferred")
    plt.ylabel("Temperature, C")
    plt.title(pinch)


    plt.show()