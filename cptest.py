# Simple State Object
from cpstate import State, Fluid
import numpy as np
from ioutils import *
import ioutils as io
import math
import matplotlib.pyplot as plt
from time import process_time

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 4}

plt.rc('font', **font)

def polytrope (state1, p2, eta, c, line=True,npt=100,props='hTs'):
    if c == 'c':
        output = io.Dict({'desc':'Output from polytrope (charging)'})
    elif c == 'd':
        output = io.Dict({'desc':'Output from polytrope (discharging)'})
    if type(state1) == Dict:
        state1 = state1.unpack('state2')[0]
    p1 = state1.p
    
    #print("Polytrope", State.T0)
    if p2 > p1 : eta = 1 / eta
    state2 = State('hp',state1.h,state1.p,fluid=state1.fluid)
    prps = np.zeros((len(props),npt))
    dp = (p2 - p1) / npt
    ro1, h1 = state1.d, state1.h

    for n in range(npt):
        p1 = state2.p; h1 = state2.h; d1 = state2.d
        p2 = p1 + dp
        for it in range(2):
            dav = 0.5* (d1 + state2.d)
            h2 = h1 + eta * dp / dav
            #print(p2, h2)
            try:
                state2.update('hp',h2,p2)
            except ValueError:
                h2 += eta * dp / dav
                p2 += dp
                state2.update('hp',h2,p2)
                n += 1

        for i, prp in enumerate(props):
            prps[i,n] = getattr(state2,prp)
        p1, h1, ro1 = p2, h2, state2.d

        if n > 99:
            print(n)

    h_bal = io.Dict({'desc':    'Polytropic Process Enthalpy Balance',
    'Availablity increase in fluid      ': state2.b - state1.b,
    'Work loss due to irreversibility   ': state1.T0 * (state2.s - state1.s),
    'Sum                                ': state2.b - state1.b + state1.T0 * (state2.s - state1.s),
    
    'Work input from motor              ': state2.h - state1.h})
    output.h_bal = h_bal
    output.wx = state2.h - state1.h
    output.wlirr_mass = state1.T0 * (state2.s - state1.s)
    output.h_delta= state2.h - state1.h
    output.b_delta = state2.b - state1.b
    output.b_analysis = (state2.h - state1.h) - (state1.T0 * (state2.s - state1.s))

    output.state1 = state1
    output.state2 = state2

    for i, prp in enumerate(props):
        output.pack(prp,prps[i])
    return output

def cryotb_test (state1, eta, line=True,npt =1000, props = 'hTs'):
    output = io.Dict({'desc':'Output from cryo turbine'})
    if type(state1) == Dict:
        state1 = state1.unpack('state2')[0]
    p1 = state1.p
    
    state2 = State('hp',state1.h,state1.p,fluid=state1.fluid)
    prps = np.empty((len(props), npt))
    dp = - 50000
    ro1, h1 = state1.d, state1.h
    n = 0
    #while state2.C != 'Liquid':
        
    while n < npt:
        p2 = p1 + dp
        h2 = h1 + eta * dp / ro1
        state2.update('hp',h2,p2)
        if state2.C == "Two phase":
            p2 = p2 - dp
            h2 = h2 - eta * dp / ro1
            state2.update('hp',h2,p2)
            break
        roav = 0.5 * (ro1 + state2.d)            
        h2 = h1 + eta * dp / roav
        state2.update('hp',h2,p2)
        for i, prp in enumerate(props):
            prps[i,n] = getattr(state2,prp)
        p1, h1, ro1 = p2, h2, state2.d
        n += 1
    
    output.wx = state2.h - state1.h
    output.wlirr_mass = state1.T0 * (state2.s - state1.s)
    output.state1 = state1
    output.state2 = state2
    output.h_delta = state2.h - state1.h
    #output.b_delta = (state2.h - state2.T0*state2.s) - (state1.h - state1.T0*state1.s)
    output.b_delta = state2.b - state1.b

    prps_new = np.array([])
    prps_new = prps[:,:n-1]

    h_bal = io.Dict({'desc':    'Polytropic Process Enthalpy Balance',
    'Availablity increase in fluid      ': state2.b - state1.b,
    'Work loss due to irreversibility   ': state1.T0 * (state2.s - state1.s),
    'Sum                                ': state2.b - state1.b + state1.T0 * (state2.s - state1.s),
    
    'Work input from motor              ': state2.h - state1.h})
    
    output.h_bal = h_bal
    for i, prp in enumerate(props):
        output.pack(prp,prps_new[i])
    return output

def hx_shape (prps, npt = 100, props = 'hTs'):
    output = io.Dict({'desc':'Output from heat exchanger pinch point'})
    h = prps[0]
    T = prps[1]

    #fig,axes = plt.subplots(1,2,figsize=(12,8))
    #ax1 = axes[0]
    #ax2 = axes[1]

    end_pinch = True
    mid_pinch = True
    mid_out_pinch = False
    in_mid_pinch = False

    convex = False
    concave = False
    convex_concave = False
    concave_convex = False

    h_norm = ((h - h[0])/(h[-1]-h[0]))*100
    T_1d = np.zeros(len(T))
    T_2d = np.zeros(len(T))
    T_d_in = np.zeros(len(T))
    T_d_out = np.zeros(len(T))
    
    for n in range(len(T)):     #Calculate gradient of TQ diagram at each point
        if n == 0:
            T_1d[n] = (T[n+1] - T[n]) / (h_norm[n+1] - h_norm[n])            
        elif n == len(T) - 1:
            T_1d[n] = (T[n] - T[n-1]) / (h_norm[n] - h_norm[n-1])
        else:           
            T_1d[n] = (T[n+1] - T[n-1]) / (h_norm[n+1] - h_norm[n-1])

    for n in range(len(T_1d)):  #See whether TQ curve is getting flatter or steeper
        if n == 0:
            T_2d[n] = (abs(T_1d[n+1]) - abs(T_1d[n])) / (h_norm[n+1] - h_norm[n])            
        elif n == len(T) - 1:
            T_2d[n] = (abs(T_1d[n]) - abs(T_1d[n-1])) / (h_norm[n] - h_norm[n-1])
        else:           
            T_2d[n] = (abs(T_1d[n+1]) - abs(T_1d[n-1])) / (h_norm[n+1] - h_norm[n-1])

    for n in range(1, len(T_2d)-1):     #determining whether pinch point(s) are at the inlet, outlet or in the middle of TQ diagram
        if T_2d[n] > 0:
            convex = True
            if T_2d[n+1] < 0: 
                convex = False
                convex_concave = True
                break
        else:
            concave = True
            if T_2d[n+1] > 0:
                concave = False
                concave_convex = True
                break
    
    if convex == True:
        TQ_shape = "convex"
    elif concave == True:
        TQ_shape = "concave"
    elif convex_concave == True:
        TQ_shape = "convex_concave"
    elif concave_convex == True:
        TQ_shape = "concave_convex"

    '''
    if in_mid_pinch == True:        #determining the overall gradient to the inlet for each point
        for n in range(len(T)):
            if n != 0:
                T_d_in[n] = (T[n] - T[0]) / (h_norm[n] - h_norm[0])
            else:
                T_d_in[n] = (T[1] - T[0]) / (h_norm[1] - h_norm[0])

        pinch_type = "Inlet-Mid"
        T_grad_diff = abs(T_1d - T_d_in)        #find the difference between gradient on the TQ line and to the inlet for each point

        for n in range(1, len(T_grad_diff)-1):
            if T_grad_diff[n-1] > T_grad_diff[n] and T_grad_diff[n] < T_grad_diff[n+1]:
                pinch = n        #the index with minimum difference and the inlet is the pinch point
                TQ_Grad = (T[n] - T[0]) / (h_norm[n] - h_norm[0])   #gradient to the inlet is the overall TQ gradient for the other side of the HX

    if mid_out_pinch == True:     #determining the overall gradient to the outlet for each point
        for n in range(len(T)):
            if n != len(T):
                T_d_out[n] = (T[-1] - T[n]) / (h_norm[-1] - h_norm[n])
            else:
                T_d_out[n] = (T[-1] - T[-2]) / (h_norm[-1] - h_norm[-2])

        pinch_type = "Mid-Outlet"
        T_grad_diff = abs(T_1d - T_d_out)       #find the difference between gradient on the TQ line and to the outlet for each point

        for n in range(1, len(T_grad_diff)-1):  #find the point on TQ that has 
            if T_grad_diff[n-1] > T_grad_diff[n] and T_grad_diff[n] < T_grad_diff[n+1]:
                pinch = n   #the index with minimum difference and the outlet is the pinch point
                TQ_Grad = (T[-1] - T[n]) / (h_norm[-1] - h_norm[n])     #gradient to the outlet is the overall TQ gradient for the other side of the HX

    else:
        TQ_Grad = (T[-1] - T[0]) / (h_norm[-1] - h_norm[0])
        T_grad_diff = abs(T_1d - TQ_Grad)
        if end_pinch == True:
            pinch = len(T) - 1
            pinch_type = "Inlet-Outlet"
        elif mid_pinch == True:
            for n in range(1, len(T_grad_diff)-1):
                if T_grad_diff[n-1] > T_grad_diff[n] and T_grad_diff[n] < T_grad_diff[n+1]:
                    pinch = n 
            pinch_type = "Mid"'''
                       
                       
    output.TQ_shape = TQ_shape

    #for n in (pinch):
        #T_grad_diff_pinch.append(T_grad_diff[n])
        #h_norm_pinch.append(h_norm[n])

    
    #ax1.plot(h_norm, T, label = "T")
    #ax2.scatter(h_norm, T_grad_diff, label = "np gradient T")
    #ax2.scatter(h_norm_pinch, T_grad_diff_pinch, c = 'r', label = "Minimum gradient different to outlet")
    #ax2.plot(h_norm, T_1d, label = "T_1d")
    #ax2.plot(h_norm, T_2d, label = "T_2d")
    #if T_d_in[50] == 0:
        #ax2.plot(h_norm, T_d_out, label = "T_d_out")
    #else:
        #ax2.plot(h_norm, T_d_in, label = "T_d_in")
    #ax2.legend()
    #plt.show()
    return output

def ex_hx (state1, T2, p_loss, c, line = True, npt = 100, props = 'hTs'):
    t1 = process_time()
    if c == 'c':
        output = io.Dict({'desc':'Output from set temp HX (charging)'})
    elif c == 'd':
        output = io.Dict({'desc':'Output from set temp HX (discharging)'})
    if type(state1) == Dict:
        state1 = state1.unpack('state2')[0]

    output.state1 = state1
    T1 = state1.T
    p1 = state1.p
    state2 = State('pT',state1.p,state1.T,fluid=state1.fluid)
    prps = np.zeros((len(props),npt))
    dT = (T2 - T1) / npt
    dp = (state1.p * p_loss) / npt

    for n in range(npt):
        T2 = T1 + dT
        p2 = p1 - dp
        state2.update('pT', p2, T2)
        for i, prp in enumerate(props):
            prps[i,n] = getattr(state2,prp)
        p1, T1 = p2, T2
        #print(state1.fluid.name())
        #print("T2: %f" %T2)
    
    h_bal = io.Dict({'desc':    'Heat Exchanger Enthalpy Balance',
    'Availablity increase in fluid'     : state2.b - state1.b,
    'Work loss due to irreversibility'  : state1.T0 * (state2.s - state1.s),
    'Sum'                               : state2.b - state1.b + state1.T0 * (state2.s - state1.s),
    
    'Enthalpy change in fluid'          : state2.h - state1.h})
    
    output.h_bal = h_bal
    output.state2 = state2
    for i, prp in enumerate(props):
        output.pack(prp,prps[i])
    t2 = process_time()
    output.time_elapsed = t2 - t1
    return output

def coupler_pinch(hot_prps, cold_fl, TP):
    output = io.Dict({'desc':'Output from coupler pinch point'})
    hot_h = hot_prps[0]
    hot_T = hot_prps[1]
    hot_s = hot_prps[2]

    cold_h = np.zeros(len(hot_h))
    cold_T = np.zeros(len(hot_T))
    cold_s = np.zeros(len(hot_T))

    dh_pinch = np.zeros(len(hot_h))
    s_pinch = np.full(len(hot_h), 100000000)
    mr_pinch = np.zeros(len(hot_h))

    d_T = np.zeros(len(hot_T))
    
    hot_h_norm = ((hot_h - hot_h[0])/(hot_h[-1]-hot_h[0]))*100

    coldin = State('pT', 100000, hot_T[0], fluid = cold_fl)             #Find pinch point based on same inlet-outlet temperature
    coldout = State('pT', 100000, hot_T[-1], fluid = cold_fl)
    dh = (coldout.h - coldin.h) / (len(hot_T) - 1)

    h1 = coldin.h
    for n in range(len(hot_T)):        
        coldout.update('hp', h1, 100000)
        cold_h[n] = coldout.h; cold_T[n] = coldout.T; cold_s[n] = coldout.s
        h1 = h1 + dh

    hot_top = True
    order = hot_T - cold_T
    for i in order:
        if i < 0:
            hot_top = False
            break

    coldshape = hx_shape([cold_h, cold_T])
    hotshape = hx_shape([hot_h, hot_T])

    print(coldshape)
    print(hotshape)
    print("hot_top:", hot_top)

    def_pinch_index = (abs(hot_T - cold_T)).argmax()
    if hot_top == True:
        def_pinch_index = 0

    cold_h_norm = ((cold_h - cold_h[0])/(cold_h[-1] - cold_h[0]))*100
    #plt.scatter(cold_h_norm, cold_T)


    #reset cold TQ to correct temperature
    coldpinch = State('pT', 100000, hot_T[def_pinch_index] - TP, fluid = cold_fl)
    cold_h[def_pinch_index] = coldpinch.h; cold_T[def_pinch_index] = coldpinch.T; cold_s[def_pinch_index] = coldpinch.s

    coldin.update('pT', 100000, hot_T[def_pinch_index] - TP)
    coldout.update('pT', 100000, hot_T[def_pinch_index] - TP)

    n = def_pinch_index - 1
    while n > -1:
        print(n, coldout.T)
        h = coldout.h - dh
        coldout.update('hp', h, 100000)
        cold_h[n] = coldout.h; cold_T[n] = coldout.T; cold_s[n] = coldout.s
        n -= 1
    
    n = def_pinch_index + 1
    while n < len(hot_T):
        h = coldin.h + dh
        coldin.update('hp', h, 100000)
        cold_h[n] = coldin.h; cold_T[n] = coldin.T; cold_s[n] = coldin.s
        n += 1

    #cold_h = np.flip(cold_h); cold_T = np.flip(cold_T); cold_s = np.flip(cold_s)
    cold_h_norm = ((cold_h - cold_h[0])/(cold_h[-1] - cold_h[0]))*100
    mr_pinch[def_pinch_index] = (hot_h[0]-hot_h[-1])/(cold_h[0]-cold_h[-1])
    s_pinch[def_pinch_index] = mr_pinch[def_pinch_index]*(cold_s[0]-cold_s[-1]) + (hot_s[-1]-hot_s[0])
    dh_pinch[def_pinch_index] = dh

    def_hot_h_norm = hot_h_norm;  def_hot_T = hot_T;  def_cold_h_norm = cold_h_norm;    def_cold_T = cold_T
    def_dh = dh

    if hotshape.TQ_shape == "convex_concave" or hotshape.TQ_shape == "concave_convex":
        if coldshape.TQ_shape == "concave":
            #test pitch point closer to cold outlet
            for i in range(1,11):
                new_pinch = def_pinch_index - i
                pinch_factor = abs(def_hot_T[new_pinch] - def_cold_T[new_pinch])/TP
                dh = dh_pinch[new_pinch + 1] * pinch_factor    #step change in enthalpy is now bigger, making cold side TQ steeper
                coldpinch = State('pT', 100000, def_hot_T[new_pinch] - TP, fluid = cold_fl)
                coldout = State('hp', coldpinch.h - ((new_pinch - 1)*dh), 100000, fluid = cold_fl)
                coldin = State('hp', coldpinch.h + ((len(hot_h) - new_pinch - 1)*dh), 100000, fluid = cold_fl)
                if coldout.T > def_hot_T[0] - TP:
                    break

                mr_pinch[new_pinch] = (hot_h[0]-hot_h[-1])/(coldout.h-coldin.h)
                s_pinch[new_pinch] = mr_pinch[new_pinch]*(coldout.s - coldin.s) + (hot_s[-1]-hot_s[0])
                dh_pinch[new_pinch] = dh
            
            #test pitch point closer to cold inlet
            for i in range(1,11):
                new_pinch = def_pinch_index + i
                pinch_factor = abs(def_hot_T[new_pinch] - def_cold_T[new_pinch])/TP
                dh = dh_pinch[new_pinch - 1] / pinch_factor    #step change in enthalpy is now smaller, making cold side TQ less steep
                coldpinch = State('pT', 100000, def_hot_T[new_pinch] - TP, fluid = cold_fl)
                coldout = State('hp', coldpinch.h - (new_pinch*dh), 100000, fluid = cold_fl)
                coldin = State('hp', coldpinch.h + ((len(hot_h) - new_pinch)*dh), 100000, fluid = cold_fl)
                if coldin.T > def_hot_T[-1] - TP:
                    break

                mr_pinch[new_pinch] = (hot_h[0]-hot_h[-1])/(coldout.h-coldin.h)
                s_pinch[new_pinch] = mr_pinch[new_pinch]*(coldout.s - coldin.s) + (hot_s[-1]-hot_s[0])
                dh_pinch[new_pinch] = dh

    #print("hot_h:", hot_h)
    #print("cold_h:", cold_h)
    #test pitch point closer to cold outlet
    print("s_pinch:", s_pinch)
    print("dh_pinch:", dh_pinch)
    print("mr_pinch:", mr_pinch)

    min_s_pinch = (s_pinch).argmin()
    output.pinch_pt = min_s_pinch
    output.mass_r = mr_pinch[min_s_pinch]
    output.dh = dh_pinch[min_s_pinch]
    


    pinch_temp = str(min(abs(hot_T - cold_T)))


    #for n in range(len(hot_h_norm)):
        #plt.plot([cold_h_norm[n], hot_h_norm[n]],[(cold_T + (hot_T-cold_T)[pinch_index] - 10)[n], hot_T[n]], c = "black")
    #plt.scatter(cold_h_norm, cold_T)
    #plt.scatter(def_hot_h_norm, def_hot_T, label = "hot")
    #plt.scatter(def_cold_h_norm, def_cold_T, label = "cold")
    #plt.title(str(min(abs(def_hot_T - def_cold_T))))
    #plt.legend()
    #plt.scatter([hot_h_norm[pinch_index], cold_h_norm[pinch_index]],[hot_T[pinch_index], (cold_T + (hot_T-cold_T)[pinch_index] - 10)[pinch_index]])
    #plt.show()

    return output

def coupler (hotin, T2, coldin, TP, p_loss, line = True, npt = 100, props = 'hTs'):        #currently only works for charging
    output = io.Dict({'desc':'Output from coupler'})
    if type(hotin) == Dict:
        hotin = hotin.unpack('state2')[0]
    if type(coldin) == Dict:
        coldin = coldin.unpack('state2')[0]
    h1 = hotin.h
    p1 = hotin.p
    hotout = State('pT', hotin.p, T2, fluid = hotin.fluid)
    hot_h = np.zeros(npt); hot_T = np.zeros(npt); hot_s = np.zeros(npt)
    dh = (hotout.h - h1) / (npt - 1)
    dp = (hotin.p * p_loss) / (npt - 1)

    hotout = State('pT', hotin.p, hotin.T, fluid = hotin.fluid)
    for n in range(npt):
        hotout.update('hp', h1, p1)
        hot_h[n] = hotout.h; hot_T[n] = hotout.T; hot_s[n] = hotout.s
        p1, h1 = p1 - dp, h1 + dh
    
    output.hotout = hotout
    output.hot_h = hot_h
    output.hot_T = hot_T
    output.hot_s = hot_s
    
    prps = np.array([hot_h, hot_T, hot_s])
    pinch_setting = coupler_pinch(prps, coldin.fluid, TP)
    print(pinch_setting)
    

    cold_h = np.zeros(npt); cold_T = np.zeros(npt); cold_s = np.zeros(npt)
    
    '''
    if hx_setting.pinch_type == "Mid-Outlet" or hx_setting.pinch_type == "Inlet-Outlet":
        coldin = State('pT', coldin.p, hotout.T-TP, fluid = coldin.fluid)
        T_coldout = coldin.T - (hx_setting.TQ_gradient * 100)
    elif hx_setting.pinch_type == "Mid" or hx_setting.pinch_type == "Inlet-Mid":
        n = hx_setting.pinch_index
        T_pinch = hot_T[n] - TP
        T_coldin = T_pinch + hx_setting.TQ_gradient * (100 - n)
        coldin = State('pT', coldin.p, T_coldin, fluid = coldin.fluid)
        T_coldout = coldin.T - (hx_setting.TQ_gradient * 100)'''
    pinch_pt = pinch_setting.pinch_pt
    dh = pinch_setting.dh

    coldpinch = State('pT', 100000, hot_T[pinch_pt] - TP, fluid = coldin.fluid)
    cold_h[pinch_pt] = coldpinch.h; cold_T[pinch_pt] = coldpinch.T; cold_s[pinch_pt] = coldpinch.s

    coldin = State('pT', 100000, hot_T[pinch_pt] - TP, fluid = coldin.fluid)
    coldout = State('pT', 100000, hot_T[pinch_pt] - TP, fluid = coldin.fluid)
    n = pinch_pt - 1
    while n > -1:
        h = coldout.h - dh
        coldout.update('hp', h, 100000)
        cold_h[n] = coldout.h; cold_T[n] = coldout.T; cold_s[n] = coldout.s
        n -= 1
    
    n = pinch_pt + 1
    while n < len(hot_T):
        h = coldin.h + dh
        coldin.update('hp', h, 100000)
        cold_h[n] = coldin.h; cold_T[n] = coldin.T; cold_s[n] = coldin.s
        n += 1

    
    output.coldout = coldout
    output.coldin = coldin
    output.cold_h = cold_h; output.cold_T = cold_T; output.cold_s = cold_s

    print("cold_h:", cold_h)
    print(dh)
    mass_r = (hotin.h - hotout.h) / (coldout.h - coldin.h)
    print("mass_r", mass_r)
    print("mr_pinch", pinch_setting.mass_r)
    print(pinch_pt)
    output.mass_r = mass_r
    
    wlirr_hot = (hotout.s - hotin.s)
    wlirr_cold = (coldout.s - coldin.s)
    wlirr = wlirr_hot + mass_r * wlirr_cold
    output.wlirr_mass = wlirr
    output.wlirr_hot = wlirr_hot
    output.wlirr_cold_per_cold_flow = wlirr_cold

    output.b_delta_hot = hotout.b - hotin.b
    output.b_delta_cold = coldout.b - coldin.b
    output.b_delta = (hotout.b - hotin.b) + mass_r * (coldout.b - coldin.b)
    output.b_analysis_hot = (hotout.h - hotin.h) - (hotin.T0 * (hotout.s - hotin.s))
    output.b_analysis_cold = (coldout.h - coldin.h) - (hotin.T0 * (coldout.s - coldin.s))
    output.b_analysis = (hotout.h - hotin.h) - (hotin.T0 * (hotout.s - hotin.s)) + mass_r * ((coldout.h - coldin.h) - (hotin.T0 * (coldout.s - coldin.s)))
    return output


def valve (state1,p2,line=True,npt=100,props='hTs'):
    output = io.Dict({'desc':'Output from valve'})
    if type(state1) == Dict:
        state1 = state1.unpack('state2')[0]
    p1 = state1.p
    h1 = state1.h

    dp = (p2 - p1) / npt
    state2 = State('hp',state1.h,state1.p,fluid=state1.fluid)
    prps = np.zeros((len(props),npt))
    for n in range(npt):
        p2 = p1 + dp
        state2.update('hp',h1,p2)
        for i, prp in enumerate(props):
            prps[i,n] = getattr(state2,prp)
        p1 = p2

        print (state2.s)
    output.state2 = state2
    for i, prp in enumerate(props):
        output.pack(prp,prps[i])
    return output

    
        
    



