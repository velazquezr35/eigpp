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
from sim_db import nodeDof2idx, sfti_time, search_time, modalDof2idx
from scipy.fft import fft
import os
plt.rcParams.update({'font.size': 15})

"""
------------------------------------------------------------------------------
general opts
------------------------------------------------------------------------------
"""
if True:
    plt.close('all')

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

#Plot u / FC data
def plt_ut(struCase, dofDict, ax, **kwargs): #NOTA: ¿Cambiar nombre? (algo mixto)
    """
    Plot DOFs (u_mdr as default) or FCs as f(t)

    inputs: struCase Stru Class Obj
            dof_Dict dofDict dict {'NODE': [DOFs]} (also for FCs)
            ax matplotlib.pyplot Axes obj
    kwargs: 
            'data_type', str: 'FCS' or 'UDS' (loads or u-dofs, u-dofs as default)
            if 'UDS', may contain:
                'u_type': raw or mdr (default)
                'vel': False (default) or True (in order to calculate and plot velocities)
                'deg',bool - False (default) or True (in order to plot rotational dofs in degs) this does not change the values stored
            elif 'FCS', may contain:
            global, may contain:
                'env': False (default) or True (in order to plot the envelope)
                    
    """
    if 'data_type' in kwargs:
        data_type = kwargs.get('data_type')
    else:
        data_type = 'UDS'
    if 'u_type' in kwargs:
        u_type = kwargs.get('u_type')
    else:
        u_type = 'mdr'
    if 'vel' in kwargs:
        vel = kwargs.get('vel')
    else:
        vel = False
    if 'env' in kwargs:
        env = kwargs.get('env')
    else:
        env = False
    if 'deg' in kwargs:
        deg = kwargs.get('deg')
    else:
        deg = False

    t=struCase.t[struCase.plot_timeInds[0]:struCase.plot_timeInds[1]]
    desired_inds = nodeDof2idx(struCase, dofDict)
    original_inds = flatten_values_list(dofDict.values())
    node_labels = label_asoc(dofDict) #OJO. No contempla posibles errores (q pida algo que no tengo) y esto daría problemas. Parece no importar.
    for i in range(len(desired_inds)):
        if data_type == 'UDS':
            if u_type=='mdr':
                y = struCase.u_mdr[desired_inds[i],struCase.plot_timeInds[0]:struCase.plot_timeInds[1]]
            elif u_type=='raw':
                y = struCase.u_raw[desired_inds[i],struCase.plot_timeInds[0]:struCase.plot_timeInds[1]]
            else:
                print('Warning: Bad u_type def')
            if deg:
                if original_inds[i] in struCase.rot_inds:
                    print('rot dof to deg', original_inds[i])
                    y = np.rad2deg(y)
            if vel:
                y=np.gradient(y,t) #NOTA: Agregar al plot que es una velocidad
        elif data_type == 'FCS':
            
            y = struCase.aLoad[desired_inds[i],struCase.plot_timeInds[0]:struCase.plot_timeInds[1]]
    
        ax.plot(t,y, label= node_labels[i] +' - DOF: ' + str(original_inds[i])) #NOTA: Creo que no es necesario el transpose, lo detecta sólo.
        if env:
            high_idx, low_idx = hl_envelopes_idx(y)
            ax.plot(t[high_idx], y[high_idx])
            ax.plot(t[low_idx], y[low_idx])
    ax.legend()
    return(ax)

#Plot q data
def plt_qt(struCase, modal_inds, ax, **kwargs):
    """
    Plot modal coordinates or modal external loads as f(t).
    
    Inputs: struCase is a Stru Class Obj
            modal_inds is a list (indexes)
            ax is a matplotlib.pyplot Axes obj
    kwargs (may contain): 
        'data_type': attr name (q, Q, etc.)
        'vel': False (default) or True (in order to calculate and plot modal velocities)
        'env': False (default) or True (in order to plot the envelope)
        'modal_inds_type': 'relative' (in order to use MOI´s inds), 'absolute' (phi´s inds, default)
    """
    if 'data_type' in kwargs:
        data_type = kwargs.get('data_type')
    else:
        data_type = 'q'
    if 'vel' in kwargs:
        vel = kwargs.get('vel')
    else:
        vel = False
    if 'env' in kwargs:
        env = kwargs.get('env')
    else:
        env = False
        
    if 'modal_inds_type' in kwargs:
        modal_inds_type = kwargs.get('modal_inds_type')
    else:
        modal_inds_type = 'absolute'

    if type(modal_inds) == int:
        modal_inds = [modal_inds]
    modal_inds = handle_modal_inds(struCase, modal_inds, **kwargs)
    t=struCase.t[struCase.plot_timeInds[0]:struCase.plot_timeInds[1]]
    for loc_ind in modal_inds:
        y = getattr(struCase,data_type)[loc_ind-1,struCase.plot_timeInds[0]:struCase.plot_timeInds[1]]
        if vel:
            y=np.gradient(y,t)
        ax.plot(t,y, label=str(loc_ind))
        if env:
            high_idx, low_idx = hl_envelopes_idx(y)
            ax.plot(t[high_idx], y[high_idx])
            ax.plot(t[low_idx], y[low_idx])
    ax.legend()
    return(ax)
    
#Plot one dof for all nodes for a single t val
def plt_us(struCase, tdof_dict,ax,**kwargs):
    
    """
    Plot a (or some) DOF (u_mdr as default) or FC for all nodes in particular instants of time.
    
    Inputs: struCase is a Stru Class Obj
            tdof_dict is a dict ['DOF':[t_vals]]
            ax is a matplotlib.pyplot Axes obj
    kwargs: 
            'data_type', str: 'FCS' or 'UDS' (loads or u-dofs, u-dofs as default)
            if 'UDS', may contain:
                'u_type': raw or mdr (default)
                'vel': False (default) or True (in order to calculate and plot velocities)
                'deg',bool - False (default) or True (in order to plot rotational dofs in degs) this does not change the values stored
            elif 'FCS', may contain:
            global, may contain:           
    """
    if 'data_type' in kwargs:
        data_type = kwargs.get('data_type')
    else:
        data_type = 'UDS'
    if 'u_type' in kwargs:
        u_type = kwargs.get('u_type')
    else:
        u_type = 'mdr'
    
    if 'vel' in kwargs:
        vel = kwargs.get('vel')
    else:
        vel = False
        
    if 'deg' in kwargs:
        deg = kwargs.get('deg')
    else:
        deg = False

    desired_t = flatten_values_list(tdof_dict.values())
    inds_t = []
    for des_t in desired_t:
        inds_t.append(search_time(struCase.t,[des_t,0])[0])
    dofDict = tdof2dofDict(struCase,tdof_dict)
    stru_inds = np.array(nodeDof2idx(struCase, dofDict))
    raw_node_labels = list(dofDict.keys())
    node_labels = []
    for loc_node_label in raw_node_labels:
        pre_node_labels = [loc_node_label]*len(dofDict[loc_node_label])
        node_labels.append(pre_node_labels)
    node_labels = flatten_values_list(node_labels)
    original_inds = flatten_values_list(list(dofDict.values()))
    for i in range(len(inds_t)):
        if data_type == 'UDS':
            if u_type == 'mdr':
                y = struCase.u_mdr[stru_inds,inds_t[i]]
            elif u_type =='raw':
                y = struCase.u_raw[stru_inds,inds_t[i]]
            else:
                print('Warning: Bad u_type def')
            if deg:
                for j in range(len(original_inds)):
                    if original_inds[j] in struCase.rot_inds:
                        y[j] = np.rad2deg(y[j])
        elif data_type == 'FCS':
            y = struCase.aLoad[stru_inds, inds_t[i]]
        ax.plot(node_labels,y, label= '{0:.2f}, {1:.3f}'.format(desired_t[i],struCase.t[inds_t[i]]))
    ax.legend(title='Time instants (des vs act):')
    return(ax)

#Plot modal shapes

def plt_phi(struCase, modedofDict, ax, **kwargs):
    '''
    Plot modal shapes
    inputs:
        struCase, stru class obj
        modedofDict, dict {'mode_ind':[DOFs]}
        ax, matplotlib.pyplot Axes obj
    kwargs:
        graphs_pack, dict - Standard graphs pack
    returns:
        ax, matplotlib.pyplot Axes obj
    '''
    if 'graphs_pack' in kwargs:
        graphs_pack = kwargs.get('graphs_pack')
    else:
        graphs_pack = handle_graph_info(**kwargs)
    x_labels = nodes2labels(struCase, **kwargs)
    
    mode_inds = list(modedofDict.keys())
    
    for loc_mode in mode_inds:
        # if modedofDict[loc_mode] == int:
        for loc_ind in modedofDict[loc_mode]:
            des_inds = modalDof2idx(struCase, loc_ind, **kwargs)
            ax.plot(x_labels, struCase.phi[des_inds,int(loc_mode)-1], label = 'M: '+loc_mode + ',dof: '+str(loc_ind))
    ax.legend(title=graphs_pack['legend_title'])
    return(ax)

#Plot all modes for a single t val

def plt_qs(struCase, tmode_dict,ax,**kwargs):
    
    """
    Plot all modal coords or modal external loads in particular instants of time.
    
    Inputs: struCase is a Stru Class Obj
            tmode_dict is a dict, ['MODE':[t_vals]]
            ax is a matplotlib.pyplot Axes obj

    kwargs (may contain): 
        'data_type': 'mod_desp' or 'mod_aLoad' (mod_desp as default)
        'vel': False (default) or True (in order to calculate and plot velocities)
                    
    """
    if 'data_type' in kwargs:
        data_type = kwargs.get('data_type')
    else:
        data_type = 'mod_desp'
    if 'vel' in kwargs:
        vel = kwargs.get('vel')
    else:
        vel = False
    if data_type == 'modal_desp':
        modal_inds = np.linspace(1,len(struCase.q),len(struCase.q))
    elif data_type == 'modal_aLoad':
        modal_inds = np.linspace(1,len(struCase.Q),len(struCase.Q))
    inds_t = []
    t_lst = flatten_values_list(tmode_dict.values())
    for des_t in t_lst:
        inds_t.append(search_time(struCase.t,[des_t,0])[0])
    for i in range(len(inds_t)):
        if data_type == 'modal_desp':
            y = struCase.q[:,inds_t[i]]
        elif data_type == 'modal_aLoad':
           y = struCase.Q[:,inds_t[i]]
        ax.plot(modal_inds,y, label= '{0:.2f}, {1:.3f}'.format(t_lst[i],struCase.t[inds_t[i]]))
    ax.legend(title='Time instants (des vs act):')
    return(ax)
    
#Fourier 

def plt_uFFT(struCase, dofDict, ax, **kwargs):
    """
    Plot the FFT of a signal, u-DOF (u_mdr as default) or FC
    inputs: struCase stru class obj
            dofDict dict {'NODE': [DOFs]}
            ax matplotlib.pyplot Axes obj
    kwargs (may contain):
        'data_type', str: 'FCS' or 'UDS' (loads or u-dofs, u-dofs as default)
        if 'UDS' (may contain):
            u_type: 'mdr' or 'raw'
            vel, for in order to calculate and plot the FFT of DOF velocities
        elif 'FCS' (may contain):
        global (may contain):
            x_units, str: 'rad/s' or 'Hz' (default) - x freqs units
            graphs_pack, standard dict for plot customization     
            'x_lims', 'y_lims', list: Plotting lims
    returns:
            ax obj
    """
    if 'data_type' in kwargs:
        data_type = kwargs.get('data_type')
    else:
        data_type = 'UDS'
    if 'u_type' in kwargs:
        u_type = kwargs.get('u_type')
    else:
        u_type = 'mdr'
        
    if 'vel' in kwargs:
        vel = kwargs.get('vel')
    else:
        vel = False
    if 'x_units' in kwargs:
        x_units = kwargs.get('x_units')
    else:
        x_units = 'Hz'
    if 'graphs_pack' in kwargs:
        graphs_pack = kwargs.get('graphs_pack')
    else:
        graphs_pack = handle_graph_info(**kwargs)
    if 'y_lims' in kwargs:
        y_lims = kwargs.get('y_lims')
    else:
        y_lims = False
    if 'x_lims' in kwargs:
        x_lims = kwargs.get('x_lims')
    else:
        x_lims = False
    t = struCase.t[struCase.plot_timeInds[0]:struCase.plot_timeInds[1]]
    desired_inds = nodeDof2idx(struCase,dofDict)
    original_inds = flatten_values_list(dofDict.values())
    node_labels = label_asoc(dofDict) #OJO. No contempla posibles errores
    fDef = 1/(t[-1]-t[0])
    for i in range(len(desired_inds)):
        if data_type == 'UDS':
            if u_type=='mdr':
                y = struCase.u_mdr[desired_inds[i],struCase.plot_timeInds[0]:struCase.plot_timeInds[1]]
            elif u_type=='raw':
                y = struCase.u_raw[desired_inds[i],struCase.plot_timeInds[0]:struCase.plot_timeInds[1]]
            else:
                print('Warning: Bad u_type def')
                
            if vel:
                y=np.gradient(y,t) #NOTA: Agregar al plot que es una velocidad
        elif data_type == 'FCS':
            y = struCase.aLoad[desired_inds[i],struCase.plot_timeInds[0]:struCase.plot_timeInds[1]]
        y_f = abs(fft(y))
        loc_m = max(y_f)
        if not loc_m== 0:
            y_f = y_f/loc_m
        else:
            print('Norm Warning: 1/0 found')
        x_f = np.arange(0,fDef*(len(t)-1),fDef)
        if x_units == 'rad/s':
            x_f = x_f*2*np.pi
        ax.plot(x_f[:(len(t)-1)//2],y_f[:(len(t)-1)//2], label= node_labels[i] +' - DOF: ' + str(original_inds[i]))
    ax.set_ylabel(graphs_pack['y_label'])
    if y_lims:
        ax.set_ylim(y_lims)
    if x_lims:
        ax.set_xlim(x_lims)
    ax.legend(title=graphs_pack['legend_title'])
    return(ax)
    
def plt_qFFT(struCase, modal_inds, ax, **kwargs):
    """
    Plot the FFT of (a) q-signal(s) (modal coords or external modal loads)
    inputs: struCase stru class obj
            modal_inds list of modal indexes
            ax matplotlib.pyplot Axes obj
    kwargs (may contain): 
        'data_type': attr name (q, Q, etc...)
        'vel': False (default) or True (in order to calculate and plot velocities)
        'x_units', str: 'rad/s' or 'Hz' (default) - x freqs units
        'vel', for in order to calculate and plot the FFT of the modal velocities
        'graphs_pack', standard dict for plot customization
        'modal_inds_type': 'relative' (in order to use MOI´s inds), 'absolute' (phi´s inds, default)
        'x_lims', 'y_lims', list: Plotting lims
    returns:
            ax obj
    """
    if 'data_type' in kwargs:
        data_type = kwargs.get('data_type')
    else:
        data_type = 'q'   
    if 'vel' in kwargs:
        vel = kwargs.get('vel')
    else:
        vel = False
    if 'x_units' in kwargs:
        x_units = kwargs.get('x_units')
    else:
        x_units = 'Hz'
    if 'graphs_pack' in kwargs:
        graphs_pack = kwargs.get('graphs_pack')
    else:
        graphs_pack = handle_graph_info(**kwargs)
    if 'modal_inds_type' in kwargs:
        modal_inds_type = kwargs.get('modal_inds_type')
    else:
        modal_inds_type = 'relative'
    if 'y_lims' in kwargs:
        y_lims = kwargs.get('y_lims')
    else:
        y_lims = False
    if 'x_lims' in kwargs:
        x_lims = kwargs.get('x_lims')
    else:
        x_lims = False
    if type(modal_inds) == int:
        modal_inds = [modal_inds]
    modal_inds = handle_modal_inds(struCase, modal_inds, **kwargs)
    
    t = struCase.t[struCase.plot_timeInds[0]:struCase.plot_timeInds[1]]
    fDef = 1/(t[-1]-t[0])
    for i in modal_inds:
        y = getattr(struCase, data_type)
        if len(y.shape) == 1:
            y = y[struCase.plot_timeInds[0]:struCase.plot_timeInds[1]]
        else:
            y = y[i-1,struCase.plot_timeInds[0]:struCase.plot_timeInds[1]]
        if vel:
            y=np.gradient(y,t) #NOTA: Agregar al plot que es una velocidad
        y_f = abs(fft(y))
        loc_m = max(y_f)
        if not loc_m== 0:
            y_f = y_f/loc_m
        else:
            print('Warning: 1/0 found')
        x_f = np.arange(0,fDef*(len(t)-1),fDef)
        if x_units == 'rad/s':
            x_f = x_f*2*np.pi
        x_values = x_f[:(len(t)-1)//2]
        y_values = y_f[:(len(t)-1)//2]
        ax.plot(x_values,y_values, label=data_type + '-' + str(i)) #NOTA: Creo que no es necesario el transpose, lo detecta sólo.
    ax.set_ylabel(graphs_pack['y_label'])
    if y_lims:
        if type(y_lims) == str:
            if y_lims == 'individual':
                if not 'y_lims_FFT' in kwargs:
                    print('y_lims not set - FFT')
                else:
                    y_lims = kwargs.get('y_lims_FFT')
                    ax.set_ylim(y_lims)
        else:
            ax.set_ylim(y_lims)
    if x_lims:
        if type(x_lims) == str:
            if x_lims == 'individual':
                if not 'x_lims_FFT' in kwargs:
                    print('x_lims not set - FFT')
                else:
                    x_lims = kwargs.get('x_lims_FFT')
                    ax.set_xlim(x_lims)
        else:
            ax.set_xlim(x_lims)
    return(ax)

#Retratos de fase

def plt_uPP(struCase, dofDict,ax,**kwargs):
    """
    Plots phase-plane portraits, du/dt vs u or dFCS/dt vs FCS
    Inputs: 
            struCase is a Stru Class Obj
            dof_lst is a dict {node:[DOFs]}
            ax is a matplotlib.pyplot Axes obj
    kwargs (may contain):
        'data_type', str: 'FCS' or 'UDS' (loads or u-dofs, u-dofs as default)
        if 'UDS' (may contain):
            u_type: 'mdr' (default) or 'raw'
            'deg', bool - False (defaul) or True (in order to plot rot dofs in degs)
        elif 'FCS' (may contain):
        global (may contain):
            graphs_pack, standard dict for plot customization        
    returns:
            ax obj
    """
    if 'data_type' in kwargs:
        data_type = kwargs.get('data_type')
    else:
        data_type = 'UDS'
    if 'u_type' in kwargs:
        u_type = kwargs.get('u_type')
    else:
        u_type = 'mdr'
    if 'deg' in kwargs:
        deg = kwargs.get('deg')
    else:
        deg = False
    if 'graphs_pack' in kwargs:
        graphs_pack = kwargs.get('graphs_pack')
    else:
        graphs_pack = handle_graph_info(**kwargs)

    desired_inds = nodeDof2idx(struCase, dofDict)
    original_inds = flatten_values_list(dofDict.values())
    for i in range(len(desired_inds)):
        #NOTA: Esto se puede mejorar tomando u = todos y luego plot(u[desired])
        if data_type == 'UDS':
            if u_type=='mdr':
                y = struCase.u_mdr[desired_inds[i],struCase.plot_timeInds[0]:struCase.plot_timeInds[1]]
            elif u_type=='raw':
                y = struCase.u_raw[desired_inds[i],struCase.plot_timeInds[0]:struCase.plot_timeInds[1]]
            else:
                print('Warning: Bad u_type def')
            
            if deg:
                if original_inds[i] in struCase.rot_inds:
                    print('rot rad 2 deg,', original_inds[i])
                    y = np.rad2deg(y)
        elif data_type=='FCS':
            y = struCase.aLoad[desired_inds[i],struCase.plot_timeInds[0]:struCase.plot_timeInds[1]]
        
        dy = np.gradient(y,struCase.t[struCase.plot_timeInds[0]:struCase.plot_timeInds[1]])
        ax.plot(y,dy, label=str(original_inds[i])) #NOTA: Creo que no es necesario el transpose, lo detecta sólo.
    return(ax) #NOTA: ¿Necesito hacer el return? Quizá para actualizar

def plt_qPP(struCase, modal_inds,ax,**kwargs):
    """
    Plots phase-plane portraits, q vs dq/dt or external modal loads vs its d/dt
    Inputs: struCase is a Stru Class Obj
            modal_inds list of modal indexes
            ax is a matplotlib.pyplot Axes obj
    kwargs (may contain): 
        'data_type': 'mod_desp' or 'mod_aLoad' (mod_desp as default)
        graphs_pack, standard dict for plot customization   
        'modal_inds_type': 'relative' (in order to use MOI´s inds), 'absolute' (phi´s inds, default)
    returns:
            ax obj
                
    """
    if 'data_type' in kwargs:
        data_type = kwargs.get('data_type')
    else:
        data_type = 'mod_desp'  
    if 'graphs_pack' in kwargs:
        graphs_pack = kwargs.get('graphs_pack')
    else:
        graphs_pack = handle_graph_info(**kwargs)
    if 'modal_inds_type' in kwargs:
        modal_inds_type = kwargs.get('modal_inds_type')
    else:
        modal_inds_type = 'absolute'

    if type(modal_inds) == int:
        modal_inds = [modal_inds]
    modal_inds = handle_modal_inds(struCase, modal_inds, **kwargs)
    for loc_ind in modal_inds:
        #NOTA: Esto se puede mejorar tomando u = todos y luego plot(u[desired])
        if data_type == 'mod_desp':
            dy = np.gradient(struCase.q[loc_ind-1,struCase.plot_timeInds[0]:struCase.plot_timeInds[1]],struCase.t[struCase.plot_timeInds[0]:struCase.plot_timeInds[1]])
            y = struCase.q[loc_ind-1,struCase.plot_timeInds[0]:struCase.plot_timeInds[1]]
        elif data_type == 'mod_aLoad':
            dy = np.gradient(struCase.Q[loc_ind-1,struCase.plot_timeInds[0]:struCase.plot_timeInds[1]],struCase.t[struCase.plot_timeInds[0]:struCase.plot_timeInds[1]])
            y = struCase.Q[loc_ind-1,struCase.plot_timeInds[0]:struCase.plot_timeInds[1]]
        ax.plot(y,dy, label=str(loc_ind))
    ax.legend(graphs_pack['legend_title'])
    return(ax)

#Spectrogram

def plt_uspectr(struCase, dofDict, fig, ax, **kwargs):
    """
    Plots spectrogram of an u-signal or FCS-signal
    Inputs: struCase is a Stru Class Obj
            dof_lst is a dict (o lista, ver cual dejar), {node:[DOFs]}
            ax is a matplotlib.pyplot Axes obj
    kwargs (may contain):
        'data_type', str: 'FCS' or 'UDS' (loads or u-dofs, u-dofs as default)
        if 'UDS' (may contain):
            u_type: 'mdr' (default) or 'raw'
            'vel': bool, default False - In order to calculate and plot the modal vel's spectrogram
        elif 'FCS' (may contain):
        global (may contain):
            'SP_Winsize': str, default Dt/20 
            'SP_OvrLapFactor': int, detault 80 (%)
            'SP_WinType': str, default 'Hann' - Check supported FFT-Windows in scipy.signal
            'SP_Normalize': bool, default True - In order to normalize the spectrogram
            'y_units': str, 'Hz' or 'rad/s' - freq units
            graphs_pack, standard dict for plot customization 
            'x_lims', 'y_lims', list: Plotting lims
    returns:
            ax obj
    """
    if 'data_type' in kwargs:
        data_type = kwargs.get('data_type')
    else:
        data_type = 'UDS'
    if 'u_type' in kwargs:
        u_type = kwargs.get('u_type')
    else:
        u_type = 'mdr'
    if 'vel' in kwargs:
        vel = kwargs.get('vel')
    else:
        vel = False
    if 'graphs_pack' in kwargs:
        graphs_pack = kwargs.get('graphs_pack')
    else:
        graphs_pack = handle_graph_info(**kwargs)
    t = struCase.t[struCase.plot_timeInds[0]:struCase.plot_timeInds[1]]
    D_t = t[-1] - t[0] #NOTA: Esto asume un único dt para el struCase
    if 'SP_Winsize' in kwargs:
        WinSize = int(len(t)/kwargs.get('SP_Winsize'))
    else:
        WinSize = int(len(t)/5) #NOTA: Esto asume un único winsize para el struCase. Ver si agregar info en un dict aparte para más customización
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
    if 'y_lims' in kwargs:
        y_lims = kwargs.get('y_lims')
    else:
        y_lims = False
    if 'x_lims' in kwargs:
        x_lims = kwargs.get('x_lims')
    else:
        x_lims = False
    if 'y_units' in kwargs:
        y_units = kwargs.get('y_units')
    else:
        print('Missing f units. Set to default (Hz)')
        
    OverLap = np.round(WinSize*OverLapFactor/100)
    fDef = len(t)/D_t
    desired_inds = nodeDof2idx(struCase,dofDict)
    original_inds = flatten_values_list(dofDict.values())
    node_labels = label_asoc(dofDict) #OJO. No contempla posibles errores
    for i in range(len(desired_inds)):
        if data_type =='UDS':
            if u_type=='mdr':
                y = struCase.u_mdr[desired_inds[i],struCase.plot_timeInds[0]:struCase.plot_timeInds[1]]
            elif u_type=='raw':
                y = struCase.u_raw[desired_inds[i],struCase.plot_timeInds[0]:struCase.plot_timeInds[1]]
            else:
                print('Warning: Bad u_type def')
            if vel:
                y=np.gradient(y,t) #NOTA: Agregar al plot que es una velocidad
        elif data_type =='FCS':
            y = struCase.aLoad[desired_inds[i],struCase.plot_timeInds[0]:struCase.plot_timeInds[1]]
        F, T, S = signal.spectrogram(y, fDef,window=WinType, noverlap=OverLap,nperseg=WinSize)
        if b_norm:
            loc_m = 0
            for j in range(len(S)):
                if max(S[j,:]) > loc_m:
                    loc_m = max(S[j,:])
            if not loc_m == 0:
                S = abs(S)/loc_m
            else:
                print('Warning: 1/0 found')
        if y_units == 'rad/s':
            print('f in rad/s!')
            c = ax.pcolormesh(T,F*2*np.pi,S,shading = 'auto', cmap='gray_r')
        elif y_units == 'Hz':
            print('f in hz!')
            c = ax.pcolormesh(T,F,S,shading = 'auto', cmap='gray_r')
        fig.colorbar(c, ax = ax)
        ax.set_title('Node(s): ' + keys_to_str(dofDict) + ', DOF(s):' + str(original_inds[i]))
    if y_lims:
        ax.set_ylim(y_lims)
    if x_lims:
        ax.set_xlim(x_lims)
    ax.set_ylabel(graphs_pack['y_label'])
    return(ax)

def plt_q_spectr(struCase, modal_inds, fig, ax, **kwargs):
    """
    Plots spectrogram of modal coords or modal external loads
    Inputs: struCase is a Stru Class Obj
            dof_lst is a dict (o lista, ver cual dejar), {node:[DOFs]}
            ax is a matplotlib.pyplot Axes obj
    kwargs (may contain):
        'data_type': attr name (q, Q, etc...)
        'vel': bool, default False - In order to calculate and plot the modal vel's spectrogram
        'SP_Winsize': str, default Dt/20 
        'SP_OvrLapFactor': int, detault 80 (%)
        'SP_WinType': str, default 'Hann' - Check supported FFT-Windows in scipy.signal
        'SP_Normalize': bool, default True - In order to normalize the spectrogram
        'y_units': str, 'Hz' or 'rad/s' - freq units
        'modal_inds_type': 'relative' (in order to use MOI´s inds), 'absolute' (phi´s inds, default)
        'x_lims', 'y_lims', list: Plotting lims
    """
    if 'data_type' in kwargs:
        data_type = kwargs.get('data_type')
    else:
        data_type = 'q'  
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
        
    if 'y_lims' in kwargs:
        y_lims = kwargs.get('y_lims')
    else:
        y_lims = False
    if 'y_units' in kwargs:
        y_units = kwargs.get('y_units')
    else:
        print('Missing f units. Set to default (Hz)')
        y_units = 'Hz'
    if 'modal_inds_type' in kwargs:
        modal_inds_type = kwargs.get('modal_inds_type')
    else:
        modal_inds_type = 'absolute'
    OverLap = np.round(WinSize*OverLapFactor/100)
    fDef = len(t)/D_t
    
    if type(modal_inds) == int:
        modal_inds = [modal_inds]
    modal_inds = handle_modal_inds(struCase, modal_inds, **kwargs)
    for loc_ind in modal_inds:
        y = getattr(struCase, data_type)[loc_ind-1,struCase.plot_timeInds[0]:struCase.plot_timeInds[1]]
        if vel:
            y=np.gradient(y,t) #NOTA: Agregar al plot que es una velocidad
        F, T, S = signal.spectrogram(y, fDef,window=WinType, noverlap=OverLap,nperseg=WinSize)
        if b_norm:
            loc_m = np.max(S)
            if not loc_m == 0:
                S = abs(S)/loc_m
            else:
                print('Warning: 1/0 found')
        if y_units == 'rad/s':
            print('f in rad/s!')
            y_values = F*2*np.pi
        elif y_units == 'Hz':
            print('f in hz!')
            y_values = F
        c = ax.pcolormesh(T,y_values,S,shading = 'auto', cmap='gray_r')
        fig.colorbar(c, ax = ax)
    ax.set_title('Modo: ' + lst2str(modal_inds))
    ax.set_ylabel(graphs_pack['y_label'])
    if y_lims:
        if type(y_lims) == str:
            if y_lims == 'individual':
                if not 'y_lims_SPEC' in kwargs:
                    print('y_lims not set - SPECT')
                else:
                    y_lims = kwargs.get('y_lims_SPEC')
                ax.set_ylim(y_lims)
        else:
            ax.set_ylim(y_lims)
    return(ax)

#Vs DOFs

def plt_uxuy(struCase, vsDict, ax, **kwargs):
    '''
    Plots two DOFs or FCS for all struCase nodes for a single t val (one per curve)
    inputs:
        struCase, stru class obj - 
        vsDict, dict - Contains the two DOFs info and the time vals ['DOFs':[indexes],'t_vals':[t]]
        ax, ax obj
    kwargs:
        'data_type', str: 'FCS' or 'UDS' (loads or u-dofs, u-dofs as default)
        if 'UDS' (may contain):
            'u_type', str - 'raw' or 'mdr'
            'deg', bool - False (default) or True (in order to plot rot dofs in deg)
        elif 'FCS' (may contain):
        global (may contain):
            
    returns:
        ax obj
    '''
    if 'data_type' in kwargs:
        data_type = kwargs.get('data_type')
    else:
        data_type = 'UDS'
    if 'u_type' in kwargs:
        u_type = kwargs.get('u_type')
    else:
        u_type = 'mdr'
    if 'deg' in kwargs:
        deg = kwargs.get('deg')
    else:
        deg = False
    
    if data_type == 'UDS':
        if u_type=='mdr':
            y = struCase.u_mdr
        elif u_type=='raw':
            y = struCase.u_raw
        else:
            print('Warning: Bad u_type def')
    elif data_type == 'FCS':
        y = struCase.aLoad
    desired_t = vsDict['t_vals']
    inds_t = []
    for des_t in desired_t:
        inds_t.append(search_time(struCase.t,[des_t,0])[0])
    ux_Dict = dof2dofDict(struCase,vsDict['DOFs'][0])
    ux_inds = np.array(nodeDof2idx(struCase,ux_Dict))
    uy_Dict = dof2dofDict(struCase,vsDict['DOFs'][1])
    uy_inds = np.array(nodeDof2idx(struCase,uy_Dict))
    if deg and data_type == 'UDS':
        if vsDict['DOFs'][0] in struCase.rot_inds:
            y[ux_inds] = np.rad2deg(y[ux_inds])
        if vsDict['DOFs'][1] in struCase.rot_inds:
            y[uy_inds] = np.rad2deg(y[uy_inds])
    
    for i in range(len(inds_t)):
        ax.plot(y[ux_inds,inds_t[i]],y[uy_inds,inds_t[i]],label='{0:.2f}, {1:.3f}'.format(desired_t[i],struCase.t[inds_t[i]]))
    ax.legend()
    return(ax)

def plt_general(struCase, inds, ax, **kwargs):
    '''
    Plots general props as f(t) using some inds (not Python´s)
    inputs:
        struCase, stru class obj
        inds, list: list of inds to plot
        ax, matplotlib.pyplot ax obj
    kwargs:
        (must contain)
        special_mode, str or bool - Flag for special plots (default False), 'max_vs_q'
        loc_data_type, str: 'E, QW, LW', etc - Data type
        (may contain)
        vel, bool - False(default) or True (in order to compute and plot dy/dt)
        (if special_mode = 'max_vs_q')
        x_vals, list or numpy.array - x values to plot
    returns:
        ax, matplotlib.pyplot ax obj
    '''
    if 'marker' in kwargs:
        marker = kwargs.get('marker')
    else:
        marker = None
    if 'loc_data_type' in kwargs:
        loc_data_type = kwargs.get('loc_data_type')
    else:
        raise NameError('Data type str error')
    if 'special_mode' in kwargs:
        special_mode = kwargs.get('special_mode')
    else:
        special_mode = False
    graphs_pack = handle_graph_info(**kwargs)
    if not special_mode:
        if 'vel' in kwargs:
            vel = kwargs.get('vel')
        else:
            vel = False
        t=struCase.t[struCase.plot_timeInds[0]:struCase.plot_timeInds[1]]
        for i in range(len(inds)):
            try:
                raw_y = getattr(struCase,loc_data_type[i])
                if len(raw_y.shape) == 1:
                    y = raw_y[struCase.plot_timeInds[0]:struCase.plot_timeInds[1]]
                else:
                    y = raw_y[inds[i]-1,struCase.plot_timeInds[0]:struCase.plot_timeInds[1]]
            except:
                raise NameError('Wrong data type code - pls check')
            if vel:
                y = np.gradient(y, t)
            ax.plot(t, y, label = loc_data_type[i] + ' - ' + str(inds[i]), marker = marker)
    elif special_mode == 'max_vs_q':
        if 'x_vals' in kwargs:
            x_vals = kwargs.get('x_vals')
        else:
            if struCase.moi == []:
                x_vals = np.arange(1,struCase.nm+1,1)
            else:
                x_vals = struCase.moi
        if 'extra_inds' in kwargs:
            extra_inds = kwargs.get('extra_inds')
        else:
            extra_inds = [0,None]
        for i in range(len(inds)):
            try:
                y_vals = getattr(struCase,loc_data_type[i])
            except:
                raise NameError('Wrong data type code - pls check')
            if len(y_vals.shape) == 2:
                ax.plot(x_vals, y_vals[:,0], label = 'envMAX - ' + str(inds[i]), marker = marker)
                ax.plot(x_vals, y_vals[:,1], label = 'envMIN - ' + str(inds[i]), marker = marker)
            else:
                ax.plot(x_vals, y_vals, label = loc_data_type[i] + ' - ' + str(inds[i]), marker = marker)
    ax.legend(title=graphs_pack['legend_title'])
    return(ax)
"""
------------------------------------------------------------------------------
HIGH LEVEL PLOT functions
------------------------------------------------------------------------------
"""

def fig_special_comb(struCase, indsLIST, **kwargs):
    '''
    Plots special axs schemes using presets
    inputs:
        struCase, stru class obj
        indsLIST, nested list or list, inds (not Python´s)
    kwargs: 
        (may contain)
        #General:
        fig_save, bool - For saving purp.
            fig_save_opts, dict - Folder, filecode, etc
        sharex: matplotlib.pyplot.subplots() argument - default 'col'
        p_prow: plots per row for the global figure - default 1
        limit_tvals or limit_tinds: list or ndarray - time limits for plotting, values or indexes
        #Plot customization:
            fig_title
            x_label
            y_label
            y_label_setting
            legend_title
    returns:
        fig, fig obj
        
    returns:
    '''
    if 'preset' in kwargs:
        preset = kwargs.get('preset')
    else:
        preset = 'default'        
    graphs_pack = handle_graph_info(**kwargs)
    
    if len(np.shape(indsLIST)) == 1:
        indsLIST = [indsLIST]
        
    n = np.shape(indsLIST)[0]
    
    if preset == 'q_spectr_FFT':
        fig, axs = plt.subplots(3,n)
        if n == 1:
            axs[0] = plt_qt(struCase, indsLIST[0], axs[0], **kwargs)
            axs[1] = plt_q_spectr(struCase, indsLIST[0],fig, axs[1], **kwargs)
            axs[2] = plt_qFFT(struCase, indsLIST[0], axs[2], **kwargs)
        else:
            for ax, inds in zip(axs, indsLIST):
                ax = plt_qt(struCase, inds, ax, **kwargs)
                ax = plt_q_spectr(struCase, inds,ax, **kwargs)
                ax = plt_qFFT(struCase, inds, ax, **kwargs)
        fig.suptitle(graphs_pack['fig_title'])
        for ax in axs:
            ax.grid()
    else:
        print('Other presets under construction')
    return(fig)
def fig_general(struCase, indLIST, **kwargs):
    '''
    Arranges plots of general props as f(t) using some inds (not Python´s)
    inputs:
        struCase, stru class obj
        indLIST, nested list, [[IND]]
    kwargs: 
        (must contain)
        special_mode, str or bool - Flag for special plots (default False), 'max_vs_q'
        data_type, nested list of str: 'E, W_load, W_modal' - Data type (for each 0-list)
        (may contain)
        #General:
        fig_save, bool - For saving purp.
            fig_save_opts, dict - Folder, filecode, etc
        sharex: matplotlib.pyplot.subplots() argument - default 'col'
        p_prow: plots per row for the global figure - default 1
        limit_tvals or limit_tinds: list or ndarray - time limits for plotting, values or indexes
        #Plot customization:
            fig_title
            x_label
            y_label
            y_label_setting
            legend_title
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
    if 'fig_save' in kwargs:
        fig_save = kwargs.get('fig_save')
        if 'fig_save_opts' in kwargs:
            fig_save_opts = kwargs.get('fig_save_opts')
        else:
            fig_save_opts = {}
    else:
        fig_save = False
    graphs_pack = handle_graph_info(**kwargs)
    
    if not 'data_type' in kwargs:
        raise NameError('Wrong data type')
    else:
        data_type = kwargs.get('data_type')
    
    if type(data_type) == str:
        kwargs['data_type'] = list(kwargs['data_type'])
        
    
    if type(indLIST)==int:
        indLIST=list(indLIST)
    n = len(indLIST)
    
    fig, axs = plt.subplots(n, p_prow, sharex = sharex)
    if n == 1: #Esto falla si ax no es un iterable (cuando n = 1 es sólo ax, no ax[:])
        kwargs['loc_data_type'] = data_type[0]
        axs = plt_general(struCase, indLIST[0], axs, **kwargs)
        axs.set_xlabel(graphs_pack['x_label'])
        axs.set_ylabel(graphs_pack['y_label'])
        axs.grid()
    else:
        for ax, inds in zip(axs, indLIST):
            kwargs['loc_data_type'] = data_type[np.where(axs==ax)[0][0]]
            ax = plt_general(struCase, inds, ax,**kwargs)
            ax.set_ylabel(graphs_pack['y_label'])
            ax.grid()
        axs[-1].set_xlabel(graphs_pack['x_label'])
    fig.suptitle(graphs_pack['fig_title'])
    if fig_save:
        save_figure(fig,fig_save_opts,**kwargs)
    return(fig, axs)


def fig_uxuy(struCase,vsLIST, **kwargs):
    '''
    Arranges plots of DOF vs DOF or FC vs FC @t fixed
    inputs:
        struCase, stru class obj
        vsLIST, list of dicts or dict - ['DOFs':[indexes],'t_vals':[t]]
    kwargs: may contain
        #General:
        fig_save, bool - For saving purp.
            fig_save_opts, dict - Folder, filecode, etc
        sharex: matplotlib.pyplot.subplots() argument - default 'col'
        p_prow: plots per row for the global figure - default 1
        limit_tvals or limit_tinds: list or ndarray - time limits for plotting, values or indexes
        #Plot customization:
            fig_title
            x_label
            y_label
            legend_title
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
    if 'fig_save' in kwargs:
        fig_save = kwargs.get('fig_save')
        if 'fig_save_opts' in kwargs:
            fig_save_opts = kwargs.get('fig_save_opts')
        else:
            fig_save_opts = {}
    else:
        fig_save = False
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
            ax = plt_uxuy(struCase, vs_dict, ax,**kwargs)
            ax.set_ylabel(graphs_pack['y_label'])
            ax.grid()
        axs[-1].set_xlabel(graphs_pack['x_label'])
    fig.suptitle(graphs_pack['fig_title'])
    if fig_save:
        save_figure(fig,fig_save_opts,**kwargs)
    return(fig)

def fig_us(struCase, tdofLIST, **kwargs):
    '''
    Arranges plots of DOF(nodes) or FC(nodes) @t fixed
    
    struCase:   stru class object
    tdofLIST:    list of tdofDicts or tdofDict {DOF: [t_instants]} 
    kwargs: may contain
        #General:
        fig_save, bool - For saving purp.
            fig_save_opts, dict - Folder, filecode, etc
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
        
    if 'fig_save' in kwargs:
        fig_save = kwargs.get('fig_save')
        if 'fig_save_opts' in kwargs:
            fig_save_opts = kwargs.get('fig_save_opts')
        else:
            fig_save_opts = {}
    else:
        fig_save = False
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
            ax = plt_us(struCase, tdof_dict, ax,**kwargs)
            ax.set_ylabel(graphs_pack['y_label'])
            ax.grid()
        axs[-1].set_xlabel(graphs_pack['x_label'])
    fig.suptitle(graphs_pack['fig_title'])
    if fig_save:
        save_figure(fig,fig_save_opts,**kwargs)
    return(fig)

def fig_qs(struCase, tmodeLIST, **kwargs):
    '''
    Arranges plots of mode(nodes) or modal external loads (nodes) @t fixed
    
    struCase:   stru class object
    tdofLIST:    list of tmodeDicts or a sigle tmodeDict {MODE: [t_instants]} 
    kwargs: may contain
        #General:
        fig_save, bool - For saving purp.
            fig_save_opts, dict - Folder, filecode, etc
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
    if 'fig_save' in kwargs:
        fig_save = kwargs.get('fig_save')
        if 'fig_save_opts' in kwargs:
            fig_save_opts = kwargs.get('fig_save_opts')
        else:
            fig_save_opts = {}
    else:
        fig_save = False
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
            ax = plt_qs(struCase, tmode_dict, ax,**kwargs)
            ax.set_ylabel(graphs_pack['y_label'])
            ax.grid()
        axs[-1].set_xlabel(graphs_pack['x_label'])
    fig.suptitle(graphs_pack['fig_title'])
    if fig_save:
        save_figure(fig,fig_save_opts,**kwargs)
    return(fig)

def fig_ut(struCase, dofLIST, **kwargs):
    '''
    Arranges plots of u(t) or load(t)
    
    struCase:   stru class object
    dofLIST:    list of dofDicts or dofDict {NODE: [DOFs]} 
    kwargs: may contain
        #General:
        fig_save, bool - For saving purp.
            fig_save_opts, dict - Folder, filecode, etc
        sharex: matplotlib.pyplot.subplots() argument - default 'col'
        p_prow: plots per row for the global figure - default 1
        limit_tvals or limit_tinds: list or ndarray - time limits for plotting, values or indexes
        #Plot customization:
            fig_title
            x_label
            y_label
            legend_title
            plot_exp, dict - fig saving purposes
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
    if 'fig_save' in kwargs:
        fig_save = kwargs.get('fig_save')
        if 'fig_save_opts' in kwargs:
            fig_save_opts = kwargs.get('fig_save_opts')
        else:
            fig_save_opts = {}
    else:
        fig_save = False
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
            ax = plt_ut(struCase, dof_dict, ax,**kwargs)
            ax.set_ylabel(graphs_pack['y_label'])
            ax.set_aspect(aspect_ratio)
            ax.grid()
        axs[-1].set_xlabel(graphs_pack['x_label'])
    fig.suptitle(graphs_pack['fig_title'])
    if fig_save:
        save_figure(fig,fig_save_opts,**kwargs)
    return(fig)

def fig_qt(struCase, modeLIST, **kwargs):
    '''
    Arranges plots of q(t) or Q(t)
    
    struCase:   stru class object
    modeLIST:    list of modal_inds or modal_inds list of modal indexes
    kwargs: may contain
        #General:
        fig_save, bool - For saving purp.
            fig_save_opts, dict - Folder, filecode, etc
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
    if 'fig_save' in kwargs:
        fig_save = kwargs.get('fig_save')
        if 'fig_save_opts' in kwargs:
            fig_save_opts = kwargs.get('fig_save_opts')
        else:
            fig_save_opts = {}
    else:
        fig_save = False
        
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
            ax = plt_qt(struCase, dof_dict, ax,**kwargs)
            ax.set_ylabel(graphs_pack['y_label'])
            ax.grid()
        axs[-1].set_xlabel(graphs_pack['x_label'])
    fig.suptitle(graphs_pack['fig_title'])
    if fig_save:
        save_figure(fig,fig_save_opts,**kwargs)
    return(fig)

def fig_u_FFT(struCase, dofLIST, **kwargs):
    '''
    Arranges plots of FFT(u(t)) or FFT(load(t))
    
    struCase: stru class obj
    dofLIST: list of dofDicts or a single dofDict: {'NODE':[DOFs]}
    kwargs: may contain
        x_units, str: 'Hz' (Default) or 'rad/s' - x axis units
        #General:
        fig_save, bool - For saving purp.
            fig_save_opts, dict - Folder, filecode, etc
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
    if not 'x_units' in kwargs:
        kwargs['x_units'] = 'Hz'
    
    if not 'x_label' in kwargs:
        kwargs['x_label'] = '$\Omega$' + ' [' + kwargs['x_units'] + ']'
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
    if 'fig_save' in kwargs:
        fig_save = kwargs.get('fig_save')
        if 'fig_save_opts' in kwargs:
            fig_save_opts = kwargs.get('fig_save_opts')
        else:
            fig_save_opts = {}
    else:
        fig_save = False
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
            ax = plt_uFFT(struCase, dof_dict, ax,**kwargs)
            ax.set_ylabel(graphs_pack['y_label'])
            ax.grid()
        axs[-1].set_xlabel(graphs_pack['x_label'])
    fig.suptitle(graphs_pack['fig_title'])
    if fig_save:
        save_figure(fig,fig_save_opts,**kwargs)
    return(fig)

def fig_q_FFT(struCase, modeLIST, **kwargs):
    '''
    Arranges plots of FFT(q(t)) or FFT(Q(t))
    
    struCase: stru class obj
    dofLIST: list of modal_inds or modal_inds list of modal indexes
    kwargs: may contain
        x_units, str: 'Hz' (Default) or 'rad/s' - x axis units
        #General:
        fig_save, bool - For saving purp.
            fig_save_opts, dict - Folder, filecode, etc
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
    if not 'x_units' in kwargs:
        kwargs['x_units'] = 'Hz'
    if not 'x_label' in kwargs:
        kwargs['x_label'] = '$\Omega$' + ' [' + kwargs['x_units'] + ']'
    if 'fig_save' in kwargs:
        fig_save = kwargs.get('fig_save')
        if 'fig_save_opts' in kwargs:
            fig_save_opts = kwargs.get('fig_save_opts')
        else:
            fig_save_opts = {}
    else:
        fig_save = False
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
        axs.legend()
    else:
        for ax, mode_dict in zip(axs, modeLIST):
            ax = plt_qFFT(struCase, mode_dict, ax,**kwargs)
            ax.set_ylabel(graphs_pack['y_label'])
            ax.grid()
            ax.legend()
        axs[-1].set_xlabel(graphs_pack['x_label'])
    fig.suptitle(graphs_pack['fig_title'])
    if fig_save:
        save_figure(fig,fig_save_opts,**kwargs)
    return(fig)


def fig_u_spect(struCase, dofLIST, **kwargs):
    '''
    Arranges plots of Spectrogram(u(t)) or Spectrogram(load(t))
    
    struCase: stru class obj
    dofLIST: list of dofDicts or a single dofDict: {'NODE':[DOFs]}
    kwargs: may contain
        y_units, str: 'Hz' or 'rad/s' - spectrogram f units
        #General:
        fig_save, bool - For saving purp.
            fig_save_opts, dict - Folder, filecode, etc
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
        pass #Nada, quedan los indexes que ya vienen con el objeto
    if not 'y_units' in kwargs:
       kwargs['y_units'] = 'Hz'
    if not 'y_label' in kwargs:
        kwargs['y_label'] = '$\Omega$' + ' [' + kwargs['y_units'] + ']'
    if 'fig_save' in kwargs:
        fig_save = kwargs.get('fig_save')
        if 'fig_save_opts' in kwargs:
            fig_save_opts = kwargs.get('fig_save_opts')
        else:
            fig_save_opts = {}
    else:
        fig_save = False
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
            ax = plt_uspectr(struCase, dof_dict, fig, ax,**kwargs)
            ax.set_ylabel(graphs_pack['y_label'])
            ax.grid()
        axs[-1].set_xlabel(graphs_pack['x_label'])
    # fig.suptitle(graphs_pack['fig_title'])
    if fig_save:
        save_figure(fig,fig_save_opts,**kwargs)
    return(fig)

def fig_q_spect(struCase, modeLIST, **kwargs):
    '''
    Arranges plots of Spectrogram(q(t)) or Spectrogram(Q(t))
    
    struCase: stru class obj
    dofLIST: list of modal_indexes or modal_inds list of modal indexes
    kwargs: may contain
        y_units, str: 'Hz' or 'rad' - Spectrogram f units
        #General:
        fig_save, bool - For saving purp.
            fig_save_opts, dict - Folder, filecode, etc
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
    if not 'y_units' in kwargs:
       kwargs['y_units'] = 'Hz'
    if not 'y_label' in kwargs:
        kwargs['y_label'] = '$\Omega$' + ' [' + kwargs['y_units'] + ']'
    if 'fig_save' in kwargs:
        fig_save = kwargs.get('fig_save')
        if 'fig_save_opts' in kwargs:
            fig_save_opts = kwargs.get('fig_save_opts')
        else:
            fig_save_opts = {}
    else:
        fig_save = False 
    graphs_pack = handle_graph_info(**kwargs)
    if type(modeLIST) == dict: #Para comodidad end-user
        modeLIST = [modeLIST]
    n = len(modeLIST)
    fig, axs = plt.subplots(n,p_prow, sharex=sharex)
    if n == 1: #Esto falla si ax no es un iterable (cuando n = 1 es sólo ax, no ax[:])
        if type(modeLIST[0]) == int:
            modeLIST = [modeLIST]
        axs = plt_q_spectr(struCase, modeLIST[0], fig, axs, **kwargs)
        axs.set_xlabel(graphs_pack['x_label'])
        axs.set_ylabel(graphs_pack['y_label'])
        axs.grid()
    else:
        for ax, mode_dict in zip(axs, modeLIST):
            ax = plt_q_spectr(struCase, mode_dict, fig, ax,**kwargs)
            ax.set_ylabel(graphs_pack['y_label'])
            ax.grid()
        axs[-1].set_xlabel(graphs_pack['x_label'])
    # fig.suptitle(graphs_pack['fig_title'])
    if fig_save:
        save_figure(fig,fig_save_opts,**kwargs)
    return(fig)

def fig_u_PP(struCase, dofLIST, **kwargs):
    '''
    Arranges plots of PP(u(t)) or PP(load(t))
    
    struCase: stru class obj
    dofLIST: a single dofDict: {'NODE':[DOFs]}
    kwargs: may contain
        #General:
        fig_save, bool - For saving purp.
            fig_save_opts, dict - Folder, filecode, etc
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
    if 'fig_save' in kwargs:
        fig_save = kwargs.get('fig_save')
        if 'fig_save_opts' in kwargs:
            fig_save_opts = kwargs.get('fig_save_opts')
        else:
            fig_save_opts = {}
    else:
        fig_save = False 
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
        axs = plt_uPP(struCase, dofLIST[0], axs, **kwargs)
        axs.set_xlabel(graphs_pack['x_label'])
        axs.set_ylabel(graphs_pack['y_label'])
        axs.set_aspect(aspect_ratio)
        axs.grid()
    else:
        for ax, dof_dict in zip(axs, dofLIST):
            ax = plt_uPP(struCase, dof_dict, ax,**kwargs)
            ax.set_ylabel(graphs_pack['y_label'])
            ax.set_aspect(aspect_ratio)
            ax.grid()
        axs[-1].set_xlabel(graphs_pack['x_label'])
    fig.suptitle(graphs_pack['fig_title'])
    if fig_save:
        save_figure(fig,fig_save_opts,**kwargs)
    return(fig)

def fig_ut_vt_pp(struCase, dofDict, **kwargs): #NOTA: No sé si esto refleja lo solicitado en el repo.
    """
    Arranges plots of u(t), v(t) and PP or load(t), vload(t) and PP
    
    struCase: stru class obj
    dofLIST: a single dofDict: {'NODE':[DOFs]}
    kwargs: may contain
        #General:
        fig_save, bool - For saving purp.
            fig_save_opts, dict - Folder, filecode, etc
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
    if 'fig_save' in kwargs:
        fig_save = kwargs.get('fig_save')
        if 'fig_save_opts' in kwargs:
            fig_save_opts = kwargs.get('fig_save_opts')
        else:
            fig_save_opts = {}
    else:
        fig_save = False 
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
    axs[-1].set_xlabel(graphs_pack['x_label'])
    fig.suptitle(graphs_pack['fig_title'])
    if fig_save:
        save_figure(fig,fig_save_opts,**kwargs)
    return(fig)

def fig_qt_vt_pp(struCase, modal_inds, **kwargs): #NOTA: No sé si esto refleja lo solicitado en el repo.
    """
    Arranges plots of q(t), q_dot(t) and PP or Q(t), Q_dot(t) and PP
    
    struCase: stru class obj
    dofLIST: a single modal_inds list of modal indexes
    kwargs: may contain
        fig_save, bool - For saving purp.
            fig_save_opts, dict - Folder, filecode, etc
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
    if 'fig_save' in kwargs:
        fig_save = kwargs.get('fig_save')
        if 'fig_save_opts' in kwargs:
            fig_save_opts = kwargs.get('fig_save_opts')
        else:
            fig_save_opts = {}
    else:
        fig_save = False
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
    axs[-1].set_xlabel(graphs_pack['x_label'])
    fig.suptitle(graphs_pack['fig_title'])
    if fig_save:
        save_figure(fig,fig_save_opts,**kwargs)
    return(fig)

#Plot modal shapes

def fig_phi(struCase, modedofLIST, **kwargs):
    '''
    Arranges plots of modal shapes
    
    struCase:   stru class object
    mdofLIST: list or dict or a single mdofLIST [{'mode_ind':[DOFs]}]
    kwargs: may contain
        #General:
        fig_save, bool - For saving purp.
        fig_save_opts, dict - Folder, filecode, etc
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
    if 'fig_save' in kwargs:
        fig_save = kwargs.get('fig_save')
        if 'fig_save_opts' in kwargs:
            fig_save_opts = kwargs.get('fig_save_opts')
        else:
            fig_save_opts = {}
    else:
        fig_save = False
    graphs_pack = handle_graph_info(**kwargs)

    if type(modedofLIST)==dict:
        modedofLIST = [modedofLIST]
    n = len(modedofLIST)
    fig, axs = plt.subplots(n, p_prow, sharex = sharex)
    if n == 1: #Esto falla si ax no es un iterable (cuando n = 1 es sólo ax, no ax[:])
        axs = plt_phi(struCase, modedofLIST[0], axs, **kwargs)
        axs.set_xlabel(graphs_pack['x_label'])
        axs.set_ylabel(graphs_pack['y_label'])
        axs.grid()
    else:
        for ax, shape_dict in zip(axs, modedofLIST):
            ax = plt_phi(struCase, shape_dict, ax,**kwargs)
            ax.set_ylabel(graphs_pack['y_label'])
            ax.grid()
        axs[-1].set_xlabel(graphs_pack['x_label'])
    fig.suptitle(graphs_pack['fig_title'])
    if fig_save:
        save_figure(fig,fig_save_opts,**kwargs)
    return(fig)
"""
------------------------------------------------------------------------------
GENERAL TOOLS functions
------------------------------------------------------------------------------
"""

def lst2str(lst, **kwargs):
    '''
    Generates a str from a lst
    inputs:
        lst, list - items
    kwargs (may contain):
        spacer, str - change prnt style
    returns:
        str
    '''
    if 'spacer' in kwargs:
        spacer = kwargs.get('spacer') #O se dice separator? no me suena
    else:
        spacer = ','
    
    if type(lst) == list:
        ret = ''
        for loc_itm in lst:
            ret += str(loc_itm)
            if not loc_itm == lst[-1]:
                ret += ' '+spacer
    else:
        ret = str(lst)
    return(ret)
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

#Special tool
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

#Nodes to x-labels (or y)
def nodes2labels(struCase, **kwargs):
    '''
    Generates a list of strings from struCase.nodes
    inputs:
        struCase, stru class obj
    kwargs:
    returns:
        labels
    '''
    loc_lst = []
    for loc_node in struCase.nodes:
        loc_lst.append(str(loc_node))
    return(loc_lst)
    

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
        graphs_pack = {'fig_title':' ', 'x_label': 'Eje x', 'y_label': 'Eje y', 'legend_title':''}
        
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
            kwargs: 'u_type': str, raw or mdr (default)
                    'vel': bool, False (default) or True (in order to calculate and plot velocities)
                    't_pref': list, ['type',t_in, t_f] (type:str, 'index' or 'vals')
    returns:
            ndarray
    """
    if 'u_type' in kwargs:
        u_type = kwargs.get('u_type')
    else:
        u_type = 'mdr'
    
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
        if u_type=='mdr':
            u = struCase.u_mdr[desired_inds[i],plt_t_inds[0]:plt_t_inds[1]]
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

def save_figure(fig, opts,**kwargs):
    '''
    Save plots in desired dir \n
    inputs:
        fig, matplotlib fig obj - Figure
        opts, dict - Options: {'folder','filecode','filename','dpi','fig_format', 'bbox_inches','fig_dpi'}
    kwargs:
        close, bool - Close fig
    returns:
        Done
    '''
    if 'close' in kwargs:
        fig_close = kwargs.get('close')
    else:
        fig_close = False
    if 'folder' in opts:
        if not opts['folder'][-1] == '/':
            opts['folder']+='/'
        s_name = opts['folder']
    else:
        s_name = ''
    if 'filecode' in opts:
        s_name = s_name + opts['filecode']
    if 'filename' in opts:
        s_name = s_name + '_'+opts['filename']
    elif 'fig_name' in kwargs:
        s_name = s_name + '_'+kwargs.get('fig_name')
    else:
        s_name = s_name + '_name_please'
    
    if 'bbox_inches' in opts:
        bbox_inches = opts['bbox_inches']
    else:
        bbox_inches = 'tight'
    
    if 'fig_format' in opts:
        fig_format = opts['fig_format']
    else:
        fig_format = None #png

    if 'fig_dpi' in opts:
        fig_dpi = opts['fig_dpi']
    else:
        fig_dpi = None #default
    try:
        if not opts['folder'][:-1] in os.listdir():
            os.mkdir(opts['folder'])
        fig.savefig(s_name, bbox_inches = bbox_inches, format=fig_format, dpi = fig_dpi)
    except:
        raise NameError('Check fig exp')
    if fig_close:
        plt.close(fig)
    return('Done!')

#handle moi modal inds
def handle_modal_inds(struCase, modal_inds, **kwargs):
    '''
    Looks for the real modal index in struCase.moi, for labeling purposes
    '''
    if 'modal_inds_type' in kwargs:
        modal_inds_type = kwargs.get('modal_inds_type')
    else:
        modal_inds_type = 'relative'
    
    updated_modal_inds = []
    
    if modal_inds_type == 'absolute':
        for loc_ind in modal_inds:
            if loc_ind in struCase.moi:
              updated_modal_inds.append(struCase.moi.index(loc_ind)+1)
            else:
                 print('Abs modal ind not found in MOI: ', loc_ind)
        if len(updated_modal_inds)==0:
            raise ValueError('Empty modal inds! - ¿MOI? or try modal_inds_type = relative')
        else:
            return(updated_modal_inds)
    else:
        return(modal_inds)
    
def do_grid_and_labels(ax, **kwargs):
    '''
    Auto grid and legend for a set of ax
    inputs:
        ax, a set of 
    '''
    try:
        for loc_ax in ax:
            loc_ax.grid()
            loc_ax.legend()
    except:
        try:
            ax.grid()
            ax.legend()
        except:
            raise NameError('Warning: Error handling ax')
    return(ax)

"""
------------------------------------------------------------------------------
runing things
------------------------------------------------------------------------------
"""

if __name__ == '__main__':
    #Examples
    pass    