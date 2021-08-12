# -*- coding: utf-8 -*-
"""
Created on Thu May 27 11:16:39 2021

@author: Ramon Velazquez
"""
"""
------------------------------------------------------------------------------
reset spyder console
------------------------------------------------------------------------------
"""
from IPython import get_ipython
get_ipython().magic('reset -sf')


"""
------------------------------------------------------------------------------
importing zone
------------------------------------------------------------------------------
"""
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from sim_db import nodeDof2idx, sfti_time, search_time
from scipy.fft import fft
plt.rcParams.update({'font.size': 15})

"""
------------------------------------------------------------------------------
classes
------------------------------------------------------------------------------
"""

class MyClass(): # just an example
    pass

"""
------------------------------------------------------------------------------
LOW LEVEL PLOT functions
------------------------------------------------------------------------------
"""

#Plot u data

def plt_ut(struCase, dofDict, ax, **kwargs):
    """
    Plot DOFs as f(t), u_avr as default.
    
    Inputs: struCase Stru Class Obj
            dof_Dict dofDict dict {'NODE': [DOFs]}
            ax matplotlib.pyplot Axes obj
            kwargs: 'u_type': raw or avr (default)
                    'vel': False (default) or True (in order to calculate and plot velocities)
                    'env': False (default) or True (in order to plot the envelope)
                    
    """
    if 'u_type' in kwargs:
        u_type = kwargs.get('u_type')
    else:
        u_type = 'avr'
    
    if 'vel' in kwargs:
        vel = kwargs.get('vel')
    else:
        vel = False
    
    if 'env' in kwargs:
        env = kwargs.get('env')
    else:
        env = False

    t=struCase.t[struCase.plot_timeInds[0]:struCase.plot_timeInds[1]]
    desired_inds = nodeDof2idx(struCase, dofDict)
    original_inds = flatten_values_list(dofDict.values())
    node_labels = label_asoc(dofDict) #OJO. No contempla posibles errores (q pida algo que no tengo) y esto daría problemas. Parece no importar.
    for i in range(len(desired_inds)):
        if u_type=='avr':
            u = struCase.u_avr[desired_inds[i],struCase.plot_timeInds[0]:struCase.plot_timeInds[1]]
        elif u_type=='raw':
            u = struCase.u_raw[desired_inds[i],struCase.plot_timeInds[0]:struCase.plot_timeInds[1]]
        else:
            print('Warning: Bad u_type def')
            
        if vel:
            u=np.gradient(u,t) #NOTA: Agregar al plot que es una velocidad
        ax.plot(t,u, label= node_labels[i] +' - DOF: ' + str(original_inds[i])) #NOTA: Creo que no es necesario el transpose, lo detecta sólo.
        if env:
            high_idx, low_idx = hl_envelopes_idx(u)
            ax.plot(t[high_idx], u[high_idx])
            ax.plot(t[low_idx], u[low_idx])
    ax.legend(title='Node(s): ' + keys_to_str(dofDict)) #NOTA: ¿Quiero esta info como titulo?
    return(ax) #NOTA: ¿Necesito hacer el return? Quizá para actualizar
    
#Plot q data

def plt_qt(struCase, modal_inds, ax, **kwargs):
    """
    Plot modal coordinates as f(t).
    
    Inputs: struCase is a Stru Class Obj
            modal_inds is a list (indexes)
            ax is a matplotlib.pyplot Axes obj
            kwargs: 'vel': False (default) or True (in order to calculate and plot modal velocities)
                    'env': False (default) or True (in order to plot the envelope)
    """
    if 'vel' in kwargs:
        vel = kwargs.get('vel')
    else:
        vel = False
    if 'env' in kwargs:
        env = kwargs.get('env')
    else:
        env = False

    t=struCase.t[struCase.plot_timeInds[0]:struCase.plot_timeInds[1]]
    for loc_ind in modal_inds:
        u = struCase.q[loc_ind-1,struCase.plot_timeInds[0]:struCase.plot_timeInds[1]]
        if vel:
            u=np.gradient(u,t)
        ax.plot(t,u, label=str(loc_ind)) #NOTA: Creo que no es necesario el transpose, lo detecta sólo.
        if env:
            high_idx, low_idx = hl_envelopes_idx(u)
            ax.plot(t[high_idx], u[high_idx])
            ax.plot(t[low_idx], u[low_idx])
    ax.legend(title='plt_q')
    return(ax) #NOTA: ¿Necesito hacer el return? Quizá para actualizar

#Plot one dof for all nodes for a single t val

def plt_us(struCase, tdof_dict,ax,**kwargs):
    
    """
    Plot all DOFs in particular instants of time, u_avr as default.
    
    Inputs: struCase is a Stru Class Obj
            tdof_dict is a dict ['DOF':[t_vals]]
            ax is a matplotlib.pyplot Axes obj
            kwargs: 'u_type': raw or avr (default)
                    'vel': False (default) or True (in order to calculate and plot velocities)
                    
    """
    if 'u_type' in kwargs:
        u_type = kwargs.get('u_type')
    else:
        u_type = 'avr'
    
    if 'vel' in kwargs:
        vel = kwargs.get('vel')
    else:
        vel = False

    desired_t = flatten_values_list(tdof_dict.values())
    inds_t = []
    for des_t in desired_t:
        inds_t.append(search_time(struCase.t,[des_t,0])[0])
    dofDict = tdof2dofDict(struCase,tdof_dict)
    stru_inds = np.array(nodeDof2idx(struCase, dofDict))
    node_labels = list(dofDict.keys())
    
    for i in range(len(inds_t)):
        if u_type=='avr':
            u = struCase.u_avr[stru_inds,inds_t[i]]
        elif u_type=='raw':
            u = struCase.u_raw[stru_inds,inds_t[i]]
        else:
            print('Warning: Bad u_type def')
            
        ax.plot(node_labels,u, label= '{0:.2f}, {1:.3f}'.format(desired_t[i],struCase.t[inds_t[i]])) #NOTA: Creo que no es necesario el transpose, lo detecta sólo.

    ax.legend(title='Time instants (des vs act):') #NOTA: ¿Quiero esta info como titulo?
    return(ax) #NOTA: ¿Necesito hacer el return? Quizá para actualizar

#Plot all modes for a single t val

def plt_qs(struCase, tmode_dict,ax,**kwargs):
    
    """
    Plot all modal DOFs in particular instants of time, u_avr as default.
    
    Inputs: struCase is a Stru Class Obj
            tmode_dict is a dict, ['MODE':[t_vals]]
            ax is a matplotlib.pyplot Axes obj
            kwargs: 'u_type': raw or avr (default)
                    'vel': False (default) or True (in order to calculate and plot velocities)
                    
    """
    if 'vel' in kwargs:
        vel = kwargs.get('vel')
    else:
        vel = False
    
    modal_inds = np.linspace(1,len(struCase.q),len(struCase.q))
    inds_t = []
    t_lst = flatten_values_list(tmode_dict.values())
    for des_t in t_lst:
        inds_t.append(search_time(struCase.t,[des_t,0])[0])
    for i in range(len(inds_t)):
        u = struCase.q[:,inds_t[i]]
        ax.plot(modal_inds,u, label= '{0:.2f}, {1:.3f}'.format(t_lst[i],struCase.t[inds_t[i]])) #NOTA: Creo que no es necesario el transpose, lo detecta sólo.

    ax.legend(title='Time instants (des vs act):') #NOTA: ¿Quiero esta info como titulo?
    return(ax) #NOTA: ¿Necesito hacer el return? Quizá para actualizar
    
#Fourier 

def plt_uFFT(struCase, dofDict, ax, **kwargs):
    """
    Plot the FFT of a signal
    inputs: struCase stru class obj
            dofDict dict {'NODE': [DOFs]}
            ax matplotlib.pyplot Axes obj
    kwargs (may contain):
            graphs_pack, standard dict for plot customization
        
    returns:
            ax obj
    """
    if 'u_type' in kwargs:
        u_type = kwargs.get('u_type')
    else:
        u_type = 'avr'
        
    if 'vel' in kwargs:
        vel = kwargs.get('vel')
    else:
        vel = False
    if 'graphs_pack' in kwargs:
        graphs_pack = kwargs.get('graphs_pack')
    else:
        graphs_pack = handle_graph_info(**kwargs)
    
    t = struCase.t[struCase.plot_timeInds[0]:struCase.plot_timeInds[1]]
    desired_inds = nodeDof2idx(struCase,dofDict)
    original_inds = flatten_values_list(dofDict.values())
    node_labels = label_asoc(dofDict) #OJO. No contempla posibles errores
    fDef = 1/(t[-1]-t[0])
    for i in range(len(desired_inds)):
        if u_type=='avr':
            u = struCase.u_avr[desired_inds[i],struCase.plot_timeInds[0]:struCase.plot_timeInds[1]]
        elif u_type=='raw':
            u = struCase.u_raw[desired_inds[i],struCase.plot_timeInds[0]:struCase.plot_timeInds[1]]
        else:
            print('Warning: Bad u_type def')
            
        if vel:
            u=np.gradient(u,t) #NOTA: Agregar al plot que es una velocidad
        y_f = abs(fft(u))
        loc_m = max(y_f)
        if not loc_m== 0:
            y_f = y_f/loc_m
        else:
            print('Warning: 1/0 found')
        x_f = np.arange(0,fDef*(len(t)-1),fDef)
        x_f = x_f*2*np.pi
        ax.plot(x_f[:(len(t)-1)//2],y_f[:(len(t)-1)//2], label= node_labels[i] +' - DOF: ' + str(original_inds[i])) #NOTA: Creo que no es necesario el transpose, lo detecta sólo.

    ax.legend(title='Node(s): ' + keys_to_str(dofDict)) #NOTA: ¿Quiero esta info como titulo?
    ax.set_xlabel(graphs_pack['x_label'])
    ax.set_ylabel(graphs_pack['y_label'])
    ax.legend(title=graphs_pack['legend_title'])
    return(ax)
    
def plt_qFFT(struCase, modal_inds, ax, **kwargs):
    """
    Plot the FFT of (a) q-signal(s)
    inputs: struCase stru class obj
            modal_inds list of modal indexes
            ax matplotlib.pyplot Axes obj
    kwargs (may contain):
            vel, for in order to calculate and plot the FFT of the modal velocities
            graphs_pack, standard dict for plot customization
        
    returns:
            ax obj
    """
        
    if 'vel' in kwargs:
        vel = kwargs.get('vel')
    else:
        vel = False
    if 'graphs_pack' in kwargs:
        graphs_pack = kwargs.get('graphs_pack')
    else:
        graphs_pack = handle_graph_info(**kwargs)
    
    t = struCase.t[struCase.plot_timeInds[0]:struCase.plot_timeInds[1]]
    fDef = 1/(t[-1]-t[0])
    for i in modal_inds:
        q = struCase.q[i-1,struCase.plot_timeInds[0]:struCase.plot_timeInds[1]]
        if vel:
            q=np.gradient(q,t) #NOTA: Agregar al plot que es una velocidad
        y_f = abs(fft(q))
        loc_m = max(y_f)
        if not loc_m== 0:
            y_f = y_f/loc_m
        else:
            print('Warning: 1/0 found')
        
        x_f = np.arange(0,fDef*(len(t)-1),fDef)
        x_f = x_f*2*np.pi
        ax.plot(x_f[:(len(t)-1)//2],y_f[:(len(t)-1)//2], label=' - MODO: ' + str(i)) #NOTA: Creo que no es necesario el transpose, lo detecta sólo.

    ax.legend(title='plt_qFFT') #NOTA: ¿Quiero esta info como titulo?
    ax.set_xlabel(graphs_pack['x_label'])
    ax.set_ylabel(graphs_pack['y_label'])
    ax.legend(title=graphs_pack['legend_title'])
    return(ax)

#Retratos de fase

def plt_uPP(struCase, dofDict,ax,**kwargs):
    """
    Plots phase-plane portraits, du/dt vs u
    Inputs: struCase is a Stru Class Obj
            dof_lst is a dict (o lista, ver cual dejar), {node:[DOFs]}
            ax is a matplotlib.pyplot Axes obj
            kwargs: 'u_type': raw or avr (default)
    """
    if 'u_type' in kwargs:
        u_type = kwargs.get('u_type')
    else:
        u_type = 'avr'

    desired_inds = nodeDof2idx(struCase, dofDict)
    for loc_ind in desired_inds:
        #NOTA: Esto se puede mejorar tomando u = todos y luego plot(u[desired])
        if u_type=='avr':
            u = struCase.u_avr[loc_ind,struCase.plot_timeInds[0]:struCase.plot_timeInds[1]]
        elif u_type=='raw':
            u = struCase.u_raw[loc_ind,struCase.plot_timeInds[0]:struCase.plot_timeInds[1]]
        else:
            print('Warning: Bad u_type def')
        
        du = np.gradient(u,struCase.t[struCase.plot_timeInds[0]:struCase.plot_timeInds[1]])
        ax.plot(u,du, label=str(loc_ind)) #NOTA: Creo que no es necesario el transpose, lo detecta sólo.
        
    return(ax) #NOTA: ¿Necesito hacer el return? Quizá para actualizar

def plt_qPP(struCase, modal_inds,ax,**kwargs):
    """
    Plots phase-plane portraits, u vs du/dt
    Inputs: struCase is a Stru Class Obj
            modal_inds list of modal indexes
            ax is a matplotlib.pyplot Axes obj
            kwargs: 'u_type': raw or avr (default)
    """
    for loc_ind in modal_inds:
        #NOTA: Esto se puede mejorar tomando u = todos y luego plot(u[desired])
        dq = np.gradient(struCase.q[loc_ind-1,struCase.plot_timeInds[0]:struCase.plot_timeInds[1]],struCase.t[struCase.plot_timeInds[0]:struCase.plot_timeInds[1]])
        ax.plot(struCase.q[loc_ind-1,struCase.plot_timeInds[0]:struCase.plot_timeInds[1]],dq, label=str(loc_ind)) #NOTA: Creo que no es necesario el transpose, lo detecta sólo.
    ax.legend(title='plt_qPP')
    return(ax) #NOTA: ¿Necesito hacer el return? Quizá para actualizar

#Spectrogram

def plt_uspectr(struCase, dofDict, fig, ax, **kwargs):
    """
    Plots spectrogram
    Inputs: struCase is a Stru Class Obj
            dof_lst is a dict (o lista, ver cual dejar), {node:[DOFs]}
            ax is a matplotlib.pyplot Axes obj
            kwargs: 'u_type': raw or avr (default)
    """
    if 'u_type' in kwargs:
        u_type = kwargs.get('u_type')
    else:
        u_type = 'avr'
    if 'vel' in kwargs:
        vel = kwargs.get('vel')
    else:
        vel = False
    if 'graphs_pack' in kwargs:
        graphs_pack = kwargs.get('graphs_pack')
    else:
        graphs_pack = handle_graph_info(**kwargs)
    t = struCase.t[struCase.plot_timeInds[0]:struCase.plot_timeInds[1]]
    D_t = t[-1] - t[0] #NOTA: Esto asume un único t para el struCase
    if 'SP_Winsize' in kwargs:
        WinSize = int(len(t)/kwargs.get('SP_Winsize'))
    else:
        WinSize = int(len(t)/20) #NOTA: Esto asume un único winsize para el struCase. Ver si agregar info en un dict aparte para más customización
    if 'SP_OvrLapFactor' in kwargs:
        OverLapFactor = kwargs.get('SP_OvrLapFactor')
    else:
        OverLapFactor = 80 #NOTA: Default 80%
        
    if 'SP_WinType' in kwargs:
        WinType = kwargs.get('SP_WinType')
    else:
        WinType = 'hann' #NOTA: Default Hann's window. ¿Sirve?
        
    if 'SP_normalize' in kwargs:
        b_norm = kwargs.get('SP_normalize')
    else:
        b_norm = True #NOTA: Default normalizar el spectrogram
        
    OverLap = np.round(WinSize*OverLapFactor/100)
    fDef = len(t)/D_t
    desired_inds = nodeDof2idx(struCase,dofDict)
    original_inds = flatten_values_list(dofDict.values())
    node_labels = label_asoc(dofDict) #OJO. No contempla posibles errores
    for i in range(len(desired_inds)):
        if u_type=='avr':
            u = struCase.u_avr[desired_inds[i],struCase.plot_timeInds[0]:struCase.plot_timeInds[1]]
        elif u_type=='raw':
            u = struCase.u_raw[desired_inds[i],struCase.plot_timeInds[0]:struCase.plot_timeInds[1]]
        else:
            print('Warning: Bad u_type def')
        if vel:
            u=np.gradient(u,t) #NOTA: Agregar al plot que es una velocidad
        F, T, S = signal.spectrogram(u, fDef,window=WinType, noverlap=OverLap,nperseg=WinSize)
        if b_norm:
            loc_m = 0
            for j in range(len(S)):
                if max(S[j,:]) > loc_m:
                    loc_m = max(S[j,:])
            if not loc_m == 0:
                S = abs(S)/loc_m
            else:
                print('Warning: 1/0 found')
        c = ax.pcolormesh(T,F*2*np.pi,S,shading = 'auto', cmap='gray_r')
        fig.colorbar(c, ax = ax, label= 'DOF: ' + str(original_inds[i]))
        ax.set_title('Node(s): ' + keys_to_str(dofDict))
    ax.set_xlabel(graphs_pack['x_label'])
    ax.set_ylabel(graphs_pack['y_label'])
    return(ax)

def plt_qspectr(struCase, modal_inds, fig, ax, **kwargs):
    """
    Plots spectrogram
    Inputs: struCase is a Stru Class Obj
            dof_lst is a dict (o lista, ver cual dejar), {node:[DOFs]}
            ax is a matplotlib.pyplot Axes obj
            kwargs may contain:
                'vel': bool, default False - In order to calculate and plot the modal vel's spectrogram
                'SP_Winsize': str, default Dt/20 
                'SP_OvrLapFactor': int, detault 80 (%)
                'SP_WinType': str, default 'Hann' - Check supported FFT-Windows in scipy.signal
                'SP_Normalize': bool, default True - In order to normalize the spectrogram
    """
    if 'vel' in kwargs:
        vel = kwargs.get('vel')
    else:
        vel = False
    if 'graphs_pack' in kwargs:
        graphs_pack = kwargs.get('graphs_pack')
    else:
        graphs_pack = handle_graph_info(**kwargs)
    t = struCase.t[struCase.plot_timeInds[0]:struCase.plot_timeInds[1]]
    D_t = t[-1] - t[0] #NOTA: Esto asume un único t para el struCase
    if 'SP_Winsize' in kwargs:
        WinSize = int(len(t)/kwargs.get('SP_Winsize'))
    else:
        WinSize = int(len(t)/20) #NOTA: Esto asume un único winsize para el struCase. Ver si agregar info en un dict aparte para más customización
    if 'SP_OvrLapFactor' in kwargs:
        OverLapFactor = kwargs.get('SP_OvrLapFactor')
    else:
        OverLapFactor = 80 #NOTA: Default 80%
        
    if 'SP_WinType' in kwargs:
        WinType = kwargs.get('SP_WinType')
    else:
        WinType = 'hann' #NOTA: Default Hann's window. ¿Sirve?
        
    if 'SP_Normalize' in kwargs:
        b_norm = kwargs.get('SP_Normalize')
    else:
        b_norm = True #NOTA: Default normalizar el spectrogram
        
    OverLap = np.round(WinSize*OverLapFactor/100)
    fDef = len(t)/D_t
    

    for i in modal_inds:
        q = struCase.q[i-1,struCase.plot_timeInds[0]:struCase.plot_timeInds[1]]
        if vel:
            q=np.gradient(q,t) #NOTA: Agregar al plot que es una velocidad
        F, T, S = signal.spectrogram(q, fDef,window=WinType, noverlap=OverLap,nperseg=WinSize)
        if b_norm:
            loc_m = 0
            for j in range(len(S)):
                if max(S[j,:]) > loc_m:
                    loc_m = max(S[j,:])
            if not loc_m == 0:
                S = abs(S)/loc_m
            else:
                print('Warning: 1/0 found')
        c = ax.pcolormesh(T,F*2*np.pi,S,shading = 'auto', cmap='gray_r')
        fig.colorbar(c, ax = ax, label= 'MODE: ' + str(i))
    ax.set_xlabel(graphs_pack['x_label'])
    ax.set_ylabel(graphs_pack['y_label'])
    return(ax)

#Vs DOFs

def plt_uxuy(struCase, vsDict, ax, **kwargs):
    '''
    Plots two DOFs for all struCase nodes for a single t val (one per plot)
    inputs:
        struCase, stru class obj - 
        vsDict, dict - Contains the two DOFs info and the time vals ['DOFs':[indexes],'t_vals':[t]]
        ax, ax obj
    kwargs:
    returns:
        ax obj
    '''
    
    if 'u_type' in kwargs:
        u_type = kwargs.get('u_type')
    else:
        u_type = 'avr'
    if u_type=='avr':
        u = struCase.u_avr
    elif u_type=='raw':
        u = struCase.u_raw
    else:
        print('Warning: Bad u_type def')
        
    desired_t = vsDict['t_vals']
    inds_t = []
    for des_t in desired_t:
        inds_t.append(search_time(struCase.t,[des_t,0])[0])
    ux_Dict = dof2dofDict(struCase,vsDict['DOFs'][0])
    ux_inds = np.array(nodeDof2idx(struCase,ux_Dict))
    uy_Dict = dof2dofDict(struCase,vsDict['DOFs'][1])
    uy_inds = np.array(nodeDof2idx(struCase,uy_Dict))
    
    for i in range(len(inds_t)):
        ax.plot(u[ux_inds,inds_t[i]],u[uy_inds,inds_t[i]],label='{0:.2f}, {1:.3f}'.format(desired_t[i],struCase.t[inds_t[i]]))
        
    ax.legend(title='Time (des vs act)')
    return(ax)
"""
------------------------------------------------------------------------------
HIGH LEVEL PLOT functions
------------------------------------------------------------------------------
"""

def fig_uxuy(struCase,vsLIST, **kwargs):
    '''
    Arranges plots of DOF vs DOF @t fixed
    
    inputs:
        struCase, stru class obj
        vsLIST, list of vsDicts or a single vsDict {'DOFs':[ux,uy],'t_vals':[t]}
    kwargs:
    returns:
        fig, fig obj
    '''
    if 'sharex' in kwargs:
        sharex = kwargs.get('sharex')
    else:
        sharex = "col"
        
    if 'p_prow' in kwargs:
        p_prow = kwargs.get('p_prows')
    else:
        p_prow = 1
    graphs_pack = handle_graph_info(**kwargs)
    
    if type(vsLIST) == dict: #Para comodidad end-user
        vsLIST = [vsLIST]

    n = len(vsLIST)
    
    fig, axs = plt.subplots(n, p_prow, sharex = sharex)
    if n == 1: #Esto falla si ax no es un iterable (cuando n = 1 es sólo ax, no ax[:])
        axs = plt_uxuy(struCase, vsLIST[0], axs, **kwargs)
        axs.set_xlabel(graphs_pack['x_label'])
        axs.set_ylabel(graphs_pack['y_label'])
        axs.grid()
    else:
        for ax, vs_dict in zip(axs, vsLIST):
            ax = plt_uxuy(struCase, vs_dict, ax)
            ax.set_xlabel(graphs_pack['x_label'])
            ax.set_ylabel(graphs_pack['y_label'])
            ax.grid()
    fig.suptitle(graphs_pack['fig_title'])
    return(fig)

def fig_us(struCase, tdofLIST, **kwargs):
    '''
    Arranges plots of DOF(nodes) @t fixed
    
    struCase:   stru class object
    tdofLIST:    list of tdofDicts or tdofDict {DOF: [t_instants]} 
    kwargs: may contain
        #General:
        sharex: matplotlib.pyplot.subplots() argument - default 'col'
        p_prow: plots per row for the global figure - default 1
        limit_tvals or limit_tinds: list or ndarray - time limits for plotting, values or indexes
        #Plot customization:
            fig_title
            x_label
            y_label
            legend_title
    '''
    
    if 'sharex' in kwargs:
        sharex = kwargs.get('sharex')
    else:
        sharex = "col"
        
    if 'p_prow' in kwargs:
        p_prow = kwargs.get('p_prows')
    else:
        p_prow = 1
    graphs_pack = handle_graph_info(**kwargs)
    
    if type(tdofLIST) == dict: #Para comodidad end-user
        tdofLIST = [tdofLIST]

    n = len(tdofLIST)
    
    fig, axs = plt.subplots(n, p_prow, sharex = sharex)
    if n == 1: #Esto falla si ax no es un iterable (cuando n = 1 es sólo ax, no ax[:])
        axs = plt_us(struCase, tdofLIST[0], axs, **kwargs)
        axs.set_xlabel(graphs_pack['x_label'])
        axs.set_ylabel(graphs_pack['y_label'])
        axs.grid()
    else:
        for ax, tdof_dict in zip(axs, tdofLIST):
            ax = plt_us(struCase, tdof_dict, ax)
            ax.set_xlabel(graphs_pack['x_label'])
            ax.set_ylabel(graphs_pack['y_label'])
            ax.grid()
    fig.suptitle(graphs_pack['fig_title'])
    return(fig)

def fig_qs(struCase, tmodeLIST, **kwargs):
    '''
    Arranges plots of mode(nodes) @t fixed
    
    struCase:   stru class object
    tdofLIST:    list of tmodeDicts or a sigle tmodeDict {MODE: [t_instants]} 
    kwargs: may contain
        #General:
        sharex: matplotlib.pyplot.subplots() argument - default 'col'
        p_prow: plots per row for the global figure - default 1
        limit_tvals or limit_tinds: list or ndarray - time limits for plotting, values or indexes
        #Plot customization:
            fig_title
            x_label
            y_label
            legend_title
    '''
    
    if 'sharex' in kwargs:
        sharex = kwargs.get('sharex')
    else:
        sharex = "col"
        
    if 'p_prow' in kwargs:
        p_prow = kwargs.get('p_prows')
    else:
        p_prow = 1
    graphs_pack = handle_graph_info(**kwargs)

    n = len(tmodeLIST)
    
    fig, axs = plt.subplots(n, p_prow, sharex = sharex)
    if n == 1: #Esto falla si ax no es un iterable (cuando n = 1 es sólo ax, no ax[:])
        axs = plt_qs(struCase, tmodeLIST, axs, **kwargs)
        axs.set_xlabel(graphs_pack['x_label'])
        axs.set_ylabel(graphs_pack['y_label'])
        axs.grid()
    else:
        for ax, tmode_dict in zip(axs, tmodeLIST):
            ax = plt_qs(struCase, tmode_dict, ax)
            ax.set_xlabel(graphs_pack['x_label'])
            ax.set_ylabel(graphs_pack['y_label'])
            ax.grid()
    fig.suptitle(graphs_pack['fig_title'])
    return(fig)

def fig_ut(struCase, dofLIST, **kwargs):
    '''
    Arranges plots of u(t)
    
    struCase:   stru class object
    dofLIST:    list of dofDicts or dofDict {NODE: [DOFs]} 
    kwargs: may contain
        #General:
        sharex: matplotlib.pyplot.subplots() argument - default 'col'
        p_prow: plots per row for the global figure - default 1
        limit_tvals or limit_tinds: list or ndarray - time limits for plotting, values or indexes
        #Plot customization:
            fig_title
            x_label
            y_label
            legend_title
    '''
    
    if 'sharex' in kwargs:
        sharex = kwargs.get('sharex')
    else:
        sharex = "col"
        
    if 'p_prow' in kwargs:
        p_prow = kwargs.get('p_prows')
    else:
        p_prow = 1
    if 'limit_tvals' in kwargs:
        struCase = sfti_time(struCase,time_vals=kwargs.get('limit_tvals'))
    elif 'limit_tinds' in kwargs:
        struCase = sfti_time(struCase,indexes = kwargs.get('limit_tinds'))
    else:
        pass
        #Nada, quedan los indexes que ya vienen con el objeto
    if 'aspect_ratio' in kwargs:
        aspect_ratio = kwargs.get('aspect_ratio')
    else:
        aspect_ratio = 'auto'
    graphs_pack = handle_graph_info(**kwargs)
    
    if type(dofLIST) == dict: #Para comodidad end-user
        dofLIST = [dofLIST]

    n = len(dofLIST)
    
    fig, axs = plt.subplots(n, p_prow, sharex = sharex)
    if n == 1: #Esto falla si ax no es un iterable (cuando n = 1 es sólo ax, no ax[:])
        axs = plt_ut(struCase, dofLIST[0], axs, **kwargs)
        axs.set_xlabel(graphs_pack['x_label'])
        axs.set_ylabel(graphs_pack['y_label'])
        axs.set_aspect(aspect_ratio)
        axs.grid()
    else:
        for ax, dof_dict in zip(axs, dofLIST):
            ax = plt_ut(struCase, dof_dict, ax)
            ax.set_xlabel(graphs_pack['x_label'])
            ax.set_ylabel(graphs_pack['y_label'])
            ax.set_aspect(aspect_ratio)
            ax.grid()
    fig.suptitle(graphs_pack['fig_title'])
    return(fig)

def fig_qt(struCase, modeLIST, **kwargs):
    '''
    Arranges plots of q(t)
    
    struCase:   stru class object
    dofLIST:    list of modal_inds or modal_inds list of modal indexes
    kwargs: may contain
        #General:
        sharex: matplotlib.pyplot.subplots() argument - default 'col'
        p_prow: plots per row for the global figure - default 1
        limit_tvals or limit_tinds: list or ndarray - time limits for plotting, values or indexes
        #Plot customization:
            fig_title
            x_label
            y_label
            legend_title
    '''
    
    if 'sharex' in kwargs:
        sharex = kwargs.get('sharex')
    else:
        sharex = "col"
        
    if 'p_prow' in kwargs:
        p_prow = kwargs.get('p_prows')
    else:
        p_prow = 1
    if 'limit_tvals' in kwargs:
        struCase = sfti_time(struCase,time_vals=kwargs.get('limit_tvals'))
    elif 'limit_tinds' in kwargs:
        struCase = sfti_time(struCase,indexes = kwargs.get('limit_tinds'))
    else:
        pass
        #Nada, quedan los indexes que ya vienen con el objeto
    graphs_pack = handle_graph_info(**kwargs)


    n = len(modeLIST)
    
    fig, axs = plt.subplots(n, p_prow, sharex = sharex)
    
    if n == 1: #Esto falla si ax no es un iterable (cuando n = 1 es sólo ax, no ax[:])
        if type(modeLIST[0]) == int:
            modeLIST = [modeLIST]
        axs = plt_qt(struCase, modeLIST[0], axs, **kwargs)
        axs.set_xlabel(graphs_pack['x_label'])
        axs.set_ylabel(graphs_pack['y_label'])
        axs.grid()
    else:
        for ax, dof_dict in zip(axs, modeLIST):
            ax = plt_qt(struCase, dof_dict, ax)
            ax.set_xlabel(graphs_pack['x_label'])
            ax.set_ylabel(graphs_pack['y_label'])
            ax.grid()
    fig.suptitle(graphs_pack['fig_title'])
    return(fig)


def fig_u_FFT(struCase, dofLIST, **kwargs):
    '''
    Arranges plots of FFT(u(t))
    
    struCase: stru class obj
    dofLIST: list of dofDicts or a single dofDict: {'NODE':[DOFs]}
    kwargs: may contain
        #General:
        sharex: matplotlib.pyplot.subplots() argument - default 'col'
        p_prow: plots per row for the global figure - default 1
        limit_tvals or limit_tinds: list or ndarray - time limits for plotting, values or indexes
        #Plot customization:
            fig_title
            x_label
            y_label
            legend_title
    '''
    if 'sharex' in kwargs:
        sharex = kwargs.get('sharex')
    else:
        sharex = "col"
        
    if 'p_prow' in kwargs:
        p_prow = kwargs.get('p_prows')
    else:
        p_prow = 1
    if 'limit_tvals' in kwargs:
        struCase = sfti_time(struCase,time_vals=kwargs.get('limit_tvals'))
    elif 'limit_tinds' in kwargs:
        struCase = sfti_time(struCase,indexes = kwargs.get('limit_tinds'))
    else:
        pass
        #Nada, quedan los indexes que ya vienen con el objeto
    graphs_pack = handle_graph_info(**kwargs)
    if type(dofLIST) == dict: #Para comodidad end-user
        dofLIST = [dofLIST]
    n = len(dofLIST)
    fig, axs = plt.subplots(n,p_prow, sharex=sharex)
    if n == 1: #Esto falla si ax no es un iterable (cuando n = 1 es sólo ax, no ax[:])
        axs = plt_uFFT(struCase, dofLIST[0], axs, **kwargs)
        axs.set_xlabel(graphs_pack['x_label'])
        axs.set_ylabel(graphs_pack['y_label'])
        axs.grid()
    else:
        for ax, dof_dict in zip(axs, dofLIST):
            ax = plt_uFFT(struCase, dof_dict, ax)
            ax.set_xlabel(graphs_pack['x_label'])
            ax.set_ylabel(graphs_pack['y_label'])
            ax.grid()
    fig.suptitle(graphs_pack['fig_title'])
    return(fig)

def fig_q_FFT(struCase, modeLIST, **kwargs):
    '''
    Arranges plots of FFT(q(t))
    
    struCase: stru class obj
    dofLIST: list of modal_inds or modal_inds list of modal indexes
    kwargs: may contain
        #General:
        sharex: matplotlib.pyplot.subplots() argument - default 'col'
        p_prow: plots per row for the global figure - default 1
        limit_tvals or limit_tinds: list or ndarray - time limits for plotting, values or indexes
        #Plot customization:
            fig_title
            x_label
            y_label
            legend_title
    '''
    if 'sharex' in kwargs:
        sharex = kwargs.get('sharex')
    else:
        sharex = "col"
        
    if 'p_prow' in kwargs:
        p_prow = kwargs.get('p_prows')
    else:
        p_prow = 1
    if 'limit_tvals' in kwargs:
        struCase = sfti_time(struCase,time_vals=kwargs.get('limit_tvals'))
    elif 'limit_tinds' in kwargs:
        struCase = sfti_time(struCase,indexes = kwargs.get('limit_tinds'))
    else:
        pass
        #Nada, quedan los indexes que ya vienen con el objeto
        
    graphs_pack = handle_graph_info(**kwargs)
    if type(modeLIST) == dict: #Para comodidad end-user
        modeLIST = [modeLIST]
    n = len(modeLIST)
    fig, axs = plt.subplots(n,p_prow, sharex=sharex)
    if n == 1: #Esto falla si ax no es un iterable (cuando n = 1 es sólo ax, no ax[:])
        if type(modeLIST[0]) == int:
            modeLIST = [modeLIST]
        axs = plt_qFFT(struCase, modeLIST[0], axs, **kwargs)
        axs.set_xlabel(graphs_pack['x_label'])
        axs.set_ylabel(graphs_pack['y_label'])
        axs.grid()
    else:
        for ax, mode_dict in zip(axs, modeLIST):
            ax = plt_qFFT(struCase, mode_dict, ax)
            ax.set_xlabel(graphs_pack['x_label'])
            ax.set_ylabel(graphs_pack['y_label'])
            ax.grid()
    fig.suptitle(graphs_pack['fig_title'])
    return(fig)


def fig_u_spect(struCase, dofLIST, **kwargs):
    '''
    Arranges plots of Spectrogram(u(t))
    
    struCase: stru class obj
    dofLIST: list of dofDicts or a single dofDict: {'NODE':[DOFs]}
    kwargs: may contain
        #General:
        sharex: matplotlib.pyplot.subplots() argument - default 'col'
        p_prow: plots per row for the global figure - default 1
        limit_tvals or limit_tinds: list or ndarray - time limits for plotting, values or indexes
        #Plot customization:
            fig_title
            x_label
            y_label
            legend_title
    '''
    if 'sharex' in kwargs:
        sharex = kwargs.get('sharex')
    else:
        sharex = "col"
        
    if 'p_prow' in kwargs:
        p_prow = kwargs.get('p_prows')
    else:
        p_prow = 1

    if 'limit_tvals' in kwargs:
        struCase = sfti_time(struCase,time_vals=kwargs.get('limit_tvals'))
    elif 'limit_tinds' in kwargs:
        struCase = sfti_time(struCase,indexes = kwargs.get('limit_tinds'))
    else:
        pass
        #Nada, quedan los indexes que ya vienen con el objeto
        
    graphs_pack = handle_graph_info(**kwargs)
    if type(dofLIST) == dict: #Para comodidad end-user
        dofLIST = [dofLIST]
    n = len(dofLIST)
    fig, axs = plt.subplots(n,p_prow, sharex=sharex)
    if n == 1: #Esto falla si ax no es un iterable (cuando n = 1 es sólo ax, no ax[:])
        axs = plt_uspectr(struCase, dofLIST[0], fig, axs, **kwargs)
        axs.set_xlabel(graphs_pack['x_label'])
        axs.set_ylabel(graphs_pack['y_label'])
        axs.grid()
    else:
        for ax, dof_dict in zip(axs, dofLIST):
            ax = plt_uspectr(struCase, dof_dict, fig, ax)
            ax.set_xlabel(graphs_pack['x_label'])
            ax.set_ylabel(graphs_pack['y_label'])
            ax.grid()
    fig.suptitle(graphs_pack['fig_title'])
    return(fig)

def fig_q_spect(struCase, modeLIST, **kwargs):
    '''
    Arranges plots of Spectrogram(q(t))
    
    struCase: stru class obj
    dofLIST: list of modal_indexes or modal_inds list of modal indexes
    kwargs: may contain
        #General:
        sharex: matplotlib.pyplot.subplots() argument - default 'col'
        p_prow: plots per row for the global figure - default 1
        limit_tvals or limit_tinds: list or ndarray - time limits for plotting, values or indexes
        #Plot customization:
            fig_title
            x_label
            y_label
            legend_title
    '''
    if 'sharex' in kwargs:
        sharex = kwargs.get('sharex')
    else:
        sharex = "col"
        
    if 'p_prow' in kwargs:
        p_prow = kwargs.get('p_prows')
    else:
        p_prow = 1
        
    if 'limit_tvals' in kwargs:
        struCase = sfti_time(struCase,time_vals=kwargs.get('limit_tvals'))
    elif 'limit_tinds' in kwargs:
        struCase = sfti_time(struCase,indexes = kwargs.get('limit_tinds'))
    else:
        pass
        #Nada, quedan los indexes que ya vienen con el objeto
        
    graphs_pack = handle_graph_info(**kwargs)
    if type(modeLIST) == dict: #Para comodidad end-user
        modeLIST = [modeLIST]
    n = len(modeLIST)
    fig, axs = plt.subplots(n,p_prow, sharex=sharex)
    if n == 1: #Esto falla si ax no es un iterable (cuando n = 1 es sólo ax, no ax[:])
        if type(modeLIST[0]) == int:
            modeLIST = [modeLIST]
        axs = plt_qspectr(struCase, modeLIST[0], fig, axs, **kwargs)
        axs.set_xlabel(graphs_pack['x_label'])
        axs.set_ylabel(graphs_pack['y_label'])
        axs.grid()
    else:
        for ax, mode_dict in zip(axs, modeLIST):
            ax = plt_qspectr(struCase, mode_dict, fig, ax)
            ax.set_xlabel(graphs_pack['x_label'])
            ax.set_ylabel(graphs_pack['y_label'])
            ax.grid()
    fig.suptitle(graphs_pack['fig_title'])
    return(fig)

def fig_ut_vt_pp(struCase, dofDict, **kwargs): #NOTA: No sé si esto refleja lo solicitado en el repo.
    """
    Arranges plots of u(t), v(t) and PP
    
    struCase: stru class obj
    dofLIST: a single dofDict: {'NODE':[DOFs]}
    kwargs: may contain
        #General:
        sharex: matplotlib.pyplot.subplots() argument - default 'col'
        p_prow: plots per row for the global figure - default 1
        limit_tvals or limit_tinds: list or ndarray - time limits for plotting, values or indexes
        #Plot customization:
            fig_title
            x_label
            y_label
            legend_title
    """  
    
    if 'sharex' in kwargs:
        sharex = kwargs.get('sharex')
    else:
        sharex = "col"
        
    if 'p_prow' in kwargs:
        p_prow = kwargs.get('p_prows')
    else:
        p_prow = 1
        
    if 'limit_tvals' in kwargs:
        struCase = sfti_time(struCase,time_vals=kwargs.get('limit_tvals'))
    elif 'limit_tinds' in kwargs:
        struCase = sfti_time(struCase,indexes = kwargs.get('limit_tinds'))
    else:
        pass
        #Nada, quedan los indexes que ya vienen con el objeto
        
    graphs_pack = handle_graph_info(**kwargs)
    gs = gridspec.GridSpec(2, 3)
    fig = plt.figure()
    ax1 = fig.add_subplot(gs[0,0:2])
    ax2 = fig.add_subplot(gs[1,0:2])
    ax3 = fig.add_subplot(gs[:,2])
    axs = [ax1,ax2,ax3]
    plt_ut(struCase, dofDict, ax1,**kwargs)
    kwargs['vel'] = True
    plt_ut(struCase,dofDict, ax2, **kwargs)
    kwargs['vel'] = False
    plt_uPP(struCase, dofDict, ax3,**kwargs)
    for ax in axs:
        ax.grid()
    fig.suptitle(graphs_pack['fig_title'])
    return(fig)

def fig_qt_vt_pp(struCase, modal_inds, **kwargs): #NOTA: No sé si esto refleja lo solicitado en el repo.
    """
    Arranges plots of q(t), q_dot(t) and PP
    
    struCase: stru class obj
    dofLIST: a single modal_inds list of modal indexes
    kwargs: may contain
        #General:
        sharex: matplotlib.pyplot.subplots() argument - default 'col'
        p_prow: plots per row for the global figure - default 1
        limit_tvals or limit_tinds: list or ndarray - time limits for plotting, values or indexes
        #Plot customization:
            fig_title
            x_label
            y_label
            legend_title
    """  
    
    if 'sharex' in kwargs:
        sharex = kwargs.get('sharex')
    else:
        sharex = "col"
        
    if 'p_prow' in kwargs:
        p_prow = kwargs.get('p_prows')
    else:
        p_prow = 1
        
    if 'limit_tvals' in kwargs:
        struCase = sfti_time(struCase,time_vals=kwargs.get('limit_tvals'))
    elif 'limit_tinds' in kwargs:
        struCase = sfti_time(struCase,indexes = kwargs.get('limit_tinds'))
    else:
        pass
        #Nada, quedan los indexes que ya vienen con el objeto
    
    if type(modal_inds[0]) == list:
        print('qt_vt_pp plot > Warning: 1D list only!')
        return()
    graphs_pack = handle_graph_info(**kwargs)
    gs = gridspec.GridSpec(2, 3)
    fig = plt.figure()
    ax1 = fig.add_subplot(gs[0,0:2])
    ax2 = fig.add_subplot(gs[1,0:2])
    ax3 = fig.add_subplot(gs[:,2])
    axs = [ax1,ax2,ax3]
    
    plt_qt(struCase, modal_inds, ax1,**kwargs)
    kwargs['vel'] = True
    plt_qt(struCase,modal_inds, ax2, **kwargs)
    kwargs['vel'] = False
    plt_qPP(struCase, modal_inds, ax3,**kwargs)
    for ax in axs:
        ax.grid()
    fig.suptitle(graphs_pack['fig_title'])
    return(fig)


"""
------------------------------------------------------------------------------
GENERAL TOOLS functions
------------------------------------------------------------------------------
"""

#tDofDict to dofDict
def dof2dofDict(struCase, dof):
    '''
    Generates a dofDict for a desired DOF using all the available nodes in struCase
    inputs: 
        struCase, stru class obj -
        dof, int - desired dof to look for
    returns:
        dofdict, dict
    '''
    dofdict = {}
    for a in struCase.nodes:
        dofdict[str(a)]=[dof]
    return(dofdict)
    
def tdof2dofDict(struCase, tdof_dict):
    '''Generates a dofDict for a desired DOF using all the available nodes in struCase
    inputs:
        struCase, stru class obj -
        tdofDict, dict - contains the desired dof as key and some time instants {DOF:[t_vals]}
    returns:
        dofDict, dict 
    '''
    dofdict = {}
    loc_lst = []
    for j in list(tdof_dict.keys()):
            loc_lst.append(int(j))
    for a in struCase.nodes:
        dofdict[str(a)]=loc_lst
    return(dofdict)

#Search closest t index

def tldxs(struCase, desired_t):
    """
    Generates a list of indexes for the closest struCase.t time values in t arg (list, floats)
    input: struCase, a sim.stru obj and desired_t, a desired time [s]
    output: closest struCase.t index
    """
    des_ind = []
    for loc_t in desired_t:
        des_ind.append((np.abs(struCase.t-loc_t)).argmin())

    return(des_ind)

def keys_to_str(dct_keys):
    """
    Simple tool for str generation from dict_keys data
    input: dict.keys
    output: str
    """
    word = ''
    for loc_nam in list(dct_keys):
        word+=loc_nam
        word+='-'
    return(word[:-1])

def flatten_values_list(dct_values):
    """
    Simple tool for fusion nested lists (coming from dicts_values)
    input: dict.values
    output: list (dim 1)
    """
    flat_values = []
    if type(dct_values) == list:
        for sub in dct_values:
            for item in sub:
                flat_values.append(item)
    else:
        for sub in list(dct_values):
            for item in sub:
                flat_values.append(item)
    return(flat_values)

def label_asoc(dct):
    """
    Simple tool for generate a list of labels from dict keys
    input: dict {key:list}
    output: list of str
    """
    label_list = []
    for i in dct:
        label_list.append(len(dct[i])*[str(i)])
    label_list = flatten_values_list(label_list)
    return(label_list)
    
def handle_graph_info(**kwargs):
    """
    Prepares labels, titles and general things for plots
    input: just the kwargs
    output: graphs_info, dict
    """
    #General basic pack
    if 'graphs_pack' in kwargs:
        graphs_pack = kwargs.get('graphs_pack')
    else:
        graphs_pack = {'fig_title':'Default title', 'x_label': 'Eje x', 'y_label': 'Eje y', 'legend_title':'default'}
        
    #Also, add custom things:
    if 'fig_title' in kwargs:
        graphs_pack['fig_title'] = kwargs.get('fig_title')
    
    if 'x_label' in kwargs:
        graphs_pack['x_label'] = kwargs.get('x_label')
        
    if 'y_label' in kwargs:
        graphs_pack['y_label'] = kwargs.get('y_label')
        
    if 'legend_title' in kwargs:
        graphs_pack['legend_title'] = kwargs.get('legend_title')
        
    return(graphs_pack)

def hl_envelopes_idx(y, dmin=1, dmax=1, split=False):
    """
    Input : array-1D
    y: 1d-array, data signal from which to extract high and low envelopes
    dmin, dmax: int, optional, size of chunks, use this if the size of the input signal is too big
    split: bool, optional, if True, split the signal in half along its mean, might help to generate the envelope in some cases
    Output :
    lmin,lmax : high/low envelope idx of input signal s
    """

    # locals min      
    lmin = (np.diff(np.sign(np.diff(y))) > 0).nonzero()[0] + 1 
    # locals max
    lmax = (np.diff(np.sign(np.diff(y))) < 0).nonzero()[0] + 1 
    

    if split:
        # s_mid is zero if s centered around x-axis or more generally mean of signal
        s_mid = np.mean(y) 
        # pre-sorting of locals min based on relative position with respect to s_mid 
        lmin = lmin[y[lmin]<s_mid]
        # pre-sorting of local max based on relative position with respect to s_mid 
        lmax = lmax[y[lmax]>s_mid]


    # global max of dmax-chunks of locals max 
    lmin = lmin[[i+np.argmin(y[lmin[i:i+dmin]]) for i in range(0,len(lmin),dmin)]]
    # global min of dmin-chunks of locals min 
    lmax = lmax[[i+np.argmax(y[lmax[i:i+dmax]]) for i in range(0,len(lmax),dmax)]]
    
    return lmin,lmax

def dte_ut(struCase, dofDict, **kwargs):
    """
    Dof temporal evolution - Search & handle DOFs as f(t)
    
    Inputs: struCase Stru Class Obj
            dof_Dict dofDict dict {'NODE': [DOFs]}
            kwargs: 'u_type': str, raw or avr (default)
                    'vel': bool, False (default) or True (in order to calculate and plot velocities)
                    't_pref': list, ['type',t_in, t_f] (type:str, 'index' or 'vals')
    returns:
            ndarray
    """
    if 'u_type' in kwargs:
        u_type = kwargs.get('u_type')
    else:
        u_type = 'avr'
    
    if 'vel' in kwargs:
        vel = kwargs.get('vel')
    else:
        vel = False
    
    if 't_pref' in kwargs:
        t_pref = kwargs.get('t_pref')
        if t_pref[0] == 'index':
            plt_t_inds = t_pref[1:]
            if plt_t_inds[1] == -1:
                plt_t_inds[1] = None
            else:
                plt_t_inds[1] += 1
            t = struCase.t[plt_t_inds[0]:plt_t_inds[1]]   
        elif t_pref[0] == 'vals':
            plt_t_inds = search_time(struCase.t,t_pref[1:])
            if plt_t_inds[1] == -1:
                plt_t_inds[1] = None
            else:
                plt_t_inds[1] += 1
            t = struCase.t[plt_t_inds[0]:plt_t_inds[1]]     
    else:
        t = struCase.t
        plt_t_inds = [0,None]

    desired_inds = nodeDof2idx(struCase, dofDict)
    desired_u = []
    desired_u.append(t)
    for i in range(len(desired_inds)):
        if u_type=='avr':
            u = struCase.u_avr[desired_inds[i],plt_t_inds[0]:plt_t_inds[1]]
        elif u_type=='raw':
            u = struCase.u_raw[desired_inds[i],plt_t_inds[0]:plt_t_inds[1]]
        else:
            print('Warning: Bad u_type def')
            
        if vel:
            u=np.gradient(u,t) #NOTA: Agregar al plot que es una velocidad
        desired_u.append(u)
    desired_u=np.array(desired_u)
    return(desired_u)

def dte_qt(struCase, modal_inds, **kwargs):
    """
    Dof temporal evolution - Search & handle MODALs as f(t)
    
    Inputs: struCase Stru Class Obj
            modal_inds, list: Modal indexes (WARNING: Use normal indexes, not Python´s)
            kwargs: 'vel': bool, False (default) or True (in order to calculate and plot modal velocities)
                    't_pref': list, ['type',t_in, t_f] (type:str, 'index' or 'vals')
    returns:
            ndarray
    """
    if 'vel' in kwargs:
        vel = kwargs.get('vel')
    else:
        vel = False
    
    if 't_pref' in kwargs:
        t_pref = kwargs.get('t_pref')
        if t_pref[0] == 'index':
            plt_t_inds = t_pref[1:]
            if plt_t_inds[1] == -1:
                plt_t_inds[1] = None
            else:
                plt_t_inds[1] += 1
            t = struCase.t[plt_t_inds[0]:plt_t_inds[1]]  
        elif t_pref[0] == 'vals':
            plt_t_inds = search_time(struCase.t,t_pref[1:])
            if plt_t_inds[1] == -1:
                plt_t_inds[1] = None
            else:
                plt_t_inds[1] += 1
            t = struCase.t[plt_t_inds[0]:plt_t_inds[1]]
    else:
        t = struCase.t
        plt_t_inds = [0,None]

    desired_q = []
    desired_q.append(t)
    for loc_ind in modal_inds:
        q = struCase.q[loc_ind-1,plt_t_inds[0]:plt_t_inds[1]]
        if vel:
            q=np.gradient(q,t) #NOTA: Agregar al plot que es una velocidad
        desired_q.append(q)
    desired_q=np.array(desired_q)
    return(desired_q)

"""
------------------------------------------------------------------------------
runing things
------------------------------------------------------------------------------
"""

if __name__ == '__main__':
    #Examples
    pass    