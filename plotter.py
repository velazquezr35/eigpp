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
import sim_db
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
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

def plt_temp_ev(struCase, data_indic, ax, **kwargs):
    """
    Plot DOFs (u_mdr as default) or FCs (or modal coordinates or modal ext. loads) as f(t)

    Inputs: 
        struCase 'stru' class obj
        data_indic:
            dof_Dict dofDict dict {'NODE': [DOFs]} (also for FCs)
            modal_inds is a list (indexes)
        ax: matplotlib.pyplot Axes obj
    **kwargs (may contain): 
            'data_type', str: 'FCS' or 'UDS' (loads or u-dofs, u-dofs as default) or 'q' or 'Q' (modal coords. or m. ext. loads)
            if 'UDS', may contain:
                'u_type': raw or mdr (default)
                'deg',bool - False (default) or True (in order to plot rotational dofs in degs) this does not change the values stored
            elif 'q' or 'Q':
                'modal_inds_type': 'relative' (in order to use MOI´s inds), 'absolute' (phi´s inds, default)
            'env': False (default) or True (in order to plot the envelope)
            'vel': False (default) or True (in order to calculate and plot velocities)
    Returns:
        ax: matplotlib.pyplot Axes obj
    """
    
    if 'vel' in kwargs:
        vel = kwargs.get('vel')
    else:
        vel = False
    if 'env' in kwargs:
        env = kwargs.get('env')
    else:
        env = False
        
    if isinstance(data_indic, dict):
        plt_type = 'DOF_FCS'
    elif isinstance(data_indic, list):
        plt_type = 'MODAL_Q'
    else:
        raise TypeError('Wrong data_indic')
        
    t=struCase.t[struCase.plot_timeInds[0]:struCase.plot_timeInds[1]]
    
    if plt_type == 'DOF_FCS':
        data_type = kwargs.get('data_type')
        if 'data_type' == 'UDS':
            if 'u_type' in kwargs:
                u_type = kwargs.get('u_type')
            else:
                u_type = 'mdr'

        if 'deg' in kwargs:
            deg = kwargs.get('deg')
        else:
            deg = False
        desired_inds = sim_db.nodeDof2idx(struCase, data_indic)
        original_inds = flatten_values_list(data_indic.values())
        try:
            node_labels = label_asoc(data_indic) #OJO. No contempla posibles errores (q pida algo que no tengo) y esto daría problemas. Parece no importar.
        except:
            raise NameError('Warning - Check dofDict')
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
                if struCase.aLoad.shape == (0,):
                    raise TypeError('Empty aLoad prop!')
                y = struCase.aLoad[desired_inds[i],struCase.plot_timeInds[0]:struCase.plot_timeInds[1]]
        
            ax.plot(t,y, label= node_labels[i] +' - DOF: ' + str(original_inds[i])) #NOTA: Creo que no es necesario el transpose, lo detecta sólo.
            if env:
                high_idx, low_idx = sim_db.hl_envelopes_idx(y)
                ax.plot(t[high_idx], y[high_idx])
                ax.plot(t[low_idx], y[low_idx])
                
    elif plt_type == 'MODAL_Q':
        
        if 'data_type' in kwargs:
            data_type = kwargs.get('data_type')
        else:
            data_type = 'q'
            
        if 'modal_inds_type' in kwargs:
            modal_inds_type = kwargs.get('modal_inds_type')
        else:
            modal_inds_type = 'absolute'
            
        if type(data_indic) == int:
            data_indic = [data_indic]
        data_indic= handle_modal_inds(struCase, data_indic, **kwargs)
        t=struCase.t[struCase.plot_timeInds[0]:struCase.plot_timeInds[1]]
        for loc_ind in data_indic:
            y = getattr(struCase,data_type)[loc_ind-1,struCase.plot_timeInds[0]:struCase.plot_timeInds[1]]
            if vel:
                y=np.gradient(y,t)
            ax.plot(t,y, label=str(loc_ind))
            if env:
                high_idx, low_idx = sim_db.hl_envelopes_idx(y)
                ax.plot(t[high_idx], y[high_idx])
                ax.plot(t[low_idx], y[low_idx])
    ax.legend()
    return ax
    
def plt_tfixed(struCase, data_indic, ax, **kwargs):
    """
    Plot all dof coords (or ext. loads) or modal coords (or modal ext. loads) in particular instants of time.
    
    Inputs: 
        struCase: 'stru' class obj
        data_indic: dict, {'MODE/DOF':[t_vals]}
        ax: matplotlib.pyplot Axes obj

    **kwargs (must contain): 
        'data_type': str, 'FCS' or 'UDS' or 'q' or 'Q'
        (may contain):
            'vel': False (default) or True (in order to calculate and plot velocities)
            if 'UDS', may contain:
                'u_type': raw or mdr (default)
                'deg',bool - False (default) or True (in order to plot rotational dofs in degs) this does not change the values stored
            elif 'FCS', may contain:
            global, may contain:
                general plot customization kwargs
    Returns:
        ax: matplotlib.pyplot ax obj
    """    
    data_type = kwargs.get('data_type')
    
    if data_type == 'q' or data_type == 'Q':
        if data_type == 'q':
            data_sourc = struCase.q
            modal_inds = np.linspace(1,len(struCase.q),len(struCase.q))
        elif data_type == 'Q':
            data_sourc = struCase.Q
            modal_inds = np.linspace(1,len(struCase.Q),len(struCase.Q))
        inds_t = []
        t_lst = flatten_values_list(data_indic.values())
        for des_t in t_lst:
            inds_t.append(sim_db.search_time(struCase.t,[des_t,0])[0])
        for i in range(len(inds_t)):
            y = data_sourc[:, inds_t[i]]
            ax.plot(modal_inds,y, label= '{0:.2f}, {1:.3f}'.format(t_lst[i],struCase.t[inds_t[i]]))
        ax.legend(title='Time instants (des vs act):')
        
    elif data_type == 'UDS' or data_type == 'FCS':
        if 'u_type' in kwargs:
            u_type = kwargs.get('u_type')
        else:
            u_type = 'mdr'
        if 'deg' in kwargs:
            deg = kwargs.get('deg')
        else:
            deg = False
    
        desired_t = flatten_values_list(data_indic.values())
        inds_t = []
        for des_t in desired_t:
            inds_t.append(sim_db.search_time(struCase.t,[des_t,0])[0])
        dofDict = tdof2dofDict(struCase,data_indic)
        stru_inds = np.array(sim_db.nodeDof2idx(struCase, dofDict))
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
    else:
        raise NameError('Wrong data_type')    
    return ax



def plt_phi(struCase, modedofDict, ax, **kwargs):
    '''
    Plots modal shapes
    Inputs:
        struCase: 'stru' stru class obj
        modedofDict: dict, {'mode_ind':[DOFs]}
            mode_ind: Python´s ind
            DOFs: Real ind
        ax: matplotlib.pyplot Axes obj
    **kwargs (may contain):
        mode_ind_type, str, 'real' or 'py' (default)
        graphs_pack: dict - Standard graphs pack
    Returns:
        ax: matplotlib.pyplot Axes obj
    '''
    
    if 'graphs_pack' in kwargs:
        graphs_pack = kwargs.get('graphs_pack')
    else:
        graphs_pack = handle_graph_info(**kwargs)
    if 'color_scheme' in kwargs:
        color_scheme = kwargs.get('color_scheme')
    else:
        color_scheme = None
    if 'mode_ind_type' in kwargs:
        mode_ind_type = kwargs.get('mode_ind_type')
    else:
        mode_ind_type = 'py'
    
    x_labels = nodes2labels(struCase, **kwargs)
    
    mode_inds = list(modedofDict.keys())
    if mode_ind_type == 'py':
        m_offset = 0
    else:
        m_offset = 1
    
    for loc_mode in mode_inds:
        for loc_ind in modedofDict[loc_mode]:
            des_inds = sim_db.modalDof2idx(struCase, loc_ind)
            if not color_scheme == None:
                ax.plot(x_labels, struCase.phi[des_inds,int(loc_mode)-m_offset], label = 'M: '+str(int(loc_mode)-m_offset) + ',dof: '+str(loc_ind), color = color_scheme[loc_ind-1])
            else:
                ax.plot(x_labels, struCase.phi[des_inds,int(loc_mode)-m_offset], label = 'M: '+str(int(loc_mode)-m_offset) + ',dof: '+str(loc_ind))
    ax.legend(title=graphs_pack['legend_title'])
    return(ax)

def plt_FFT(struCase, data_indic, ax, **kwargs):
    """
    Plots the FFT of a signal: u-DOF (u_mdr as default) (or FC) or q (or Q)

    Inputs:
        struCase, 'stru' class obj
        data_indic:
            (DOF an): dict {'NODE':[DOFs]}
            (modal an): list [modal_inds]
        ax: matplotlib.pyplot ax obj
    **kwargs (must contain):
        'data_type', str: 'FCS' or 'UDS' or 'q' or 'Q'
        (may contain):
            if modal:
                'modal_inds_type': 'relative' (in order to use MOI´s inds), 'absolute' (phi´s inds, default)
            if 'UDS' (may contain):
                u_type: 'mdr' or 'raw'
                vel, for in order to calculate and plot the FFT of DOF velocities
            'vel': False (default) or True (in order to calculate and plot velocities)
            'x_units', str: 'rad/s' or 'Hz' (default) - x freqs units
            'vel', for in order to calculate and plot the FFT of the modal velocities
            'graphs_pack', standard dict for plot customization
    Returns:
            ax: matplotlib.pyplot ax obj
    """
    data_type = kwargs.get('data_type')
        
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
    fDef = 1/(t[-1]-t[0])
    y_data = []
    if data_type == 'q' or data_type == 'Q':
        if 'modal_inds_type' in kwargs:
            modal_inds_type = kwargs.get('modal_inds_type')
        else:
            modal_inds_type = 'relative'
        if type(data_indic) == int:
            data_indic = [data_indic]
        data_indic = handle_modal_inds(struCase, data_indic, **kwargs)
        for i in data_indic:
            y = getattr(struCase, data_type)
            if len(y.shape) == 1:
                y = y[struCase.plot_timeInds[0]:struCase.plot_timeInds[1]]
            else:
                y = y[i-1,struCase.plot_timeInds[0]:struCase.plot_timeInds[1]]
            if vel:
                y=np.gradient(y,t) #NOTA: Agregar al plot que es una velocidad
            y_data.append(y)
        original_inds = None
        node_labels = None

    elif data_type == 'UDS' or data_type == 'FCS':
        if 'u_type' in kwargs:
            u_type = kwargs.get('u_type')
        else:
            u_type = 'mdr'
        desired_inds = sim_db.nodeDof2idx(struCase,data_indic)
        original_inds = flatten_values_list(data_indic.values())
        node_labels = label_asoc(data_indic) #OJO. No contempla posibles errores
        
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
            y_data.append(y)
    for i in range(len(y_data)):
        y_f = abs(fft(y_data[i]))
        loc_m = max(y_f)
        if not loc_m== 0:
            y_f = y_f/loc_m
        else:
            print('Warning: 1/0 found - Data_ind: ' + str(i))
        x_f = np.arange(0,fDef*(len(t)-1),fDef)
        if x_units == 'rad/s':
            x_f = x_f*2*np.pi
        x_values = x_f[:(len(t)-1)//2]
        y_values = y_f[:(len(t)-1)//2]
        ax.plot(x_values,y_values, label=FFT_labl_gen(data_type, i, original_inds, node_labels))
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

def plt_PP(struCase, data_indic, ax, **kwargs):
    """
    Plots phase-plane portraits, du/dt vs u or dFCS/dt vs FCS (DOF an) or q vs dq/dt or Q vs dQ/dt (modal an)
    
    Inputs:
        struCase: 'stru' class obj
        data_indic: dict {node:[DOFs]} (DOF an) or list [modal_inds] (modal an)
        ax: matplotlib.pyplot ax obj
    **kwargs (must contain):
        data_type: str, 'FCS' or 'UDS' or 'q' or 'Q'
        (may contain):
            if modal an:
                'modal_inds_type': 'relative' (in order to use MOI´s inds), 'absolute' (phi´s inds, default)
            if 'UDS' (may contain):
                u_type: 'mdr' (default) or 'raw'
                'deg', bool - False (defaul) or True (in order to plot rot dofs in degs)
            graphs_pack, standard dict for plot customization        
    Returns:
            ax obj
    """
    
    data_type = kwargs.get('data_type')
    if 'deg' in kwargs:
        deg = kwargs.get('deg')
    else:
        deg = False
    if 'graphs_pack' in kwargs:
        graphs_pack = kwargs.get('graphs_pack')
    else:
        graphs_pack = handle_graph_info(**kwargs)
    
    y_data = []
    if data_type == 'UDS' or data_type == 'FCS':
        plot_type = 'DOF'
        if 'u_type' in kwargs:
            u_type = kwargs.get('u_type')
        else:
            u_type = 'mdr'
        desired_inds = sim_db.nodeDof2idx(struCase, data_indic)
        original_inds = flatten_values_list(data_indic.values())
        for i in range(len(desired_inds)):
            if data_type == 'UDS':
                if u_type=='mdr':
                    y = struCase.u_mdr[desired_inds[i],struCase.plot_timeInds[0]:struCase.plot_timeInds[1]]
                elif u_type=='raw':
                    y = struCase.u_raw[desired_inds[i],struCase.plot_timeInds[0]:struCase.plot_timeInds[1]]
                else:
                    raise NameError('Warning: Bad u_type def')
                
                if deg:
                    if original_inds[i] in struCase.rot_inds:
                        print('rot rad 2 deg,', original_inds[i])
                        y = np.rad2deg(y)
            elif data_type=='FCS':
                y = struCase.aLoad[desired_inds[i],struCase.plot_timeInds[0]:struCase.plot_timeInds[1]]
            y_data.append(y)           
    elif data_type == 'q' or data_type == 'Q':
        plot_type = 'modal'
        if 'modal_inds_type' in kwargs:
            modal_inds_type = kwargs.get('modal_inds_type')
        else:
            modal_inds_type = 'absolute'
        if type(data_indic) == int:
            data_indic = [data_indic]
        data_indic = handle_modal_inds(struCase, data_indic, **kwargs)
        for loc_ind in data_indic:
            if data_type == 'q':
                y = struCase.q[loc_ind-1,struCase.plot_timeInds[0]:struCase.plot_timeInds[1]]
            elif data_type == 'Q':
                y = struCase.Q[loc_ind-1,struCase.plot_timeInds[0]:struCase.plot_timeInds[1]]
            y_data.append(y)
    counter = 0
    for y in y_data:
        dy = np.gradient(y,struCase.t[struCase.plot_timeInds[0]:struCase.plot_timeInds[1]])
        if plot_type == 'modal':
            lab = str(data_indic[counter])
        elif plot_type == 'DOF':
            lab = str(original_inds[counter])
        ax.plot(y, dy, label = lab)
        counter +=1
    ax.legend(graphs_pack['legend_title'])
    return(ax)

#Spectrogram

def plt_spectr(struCase, data_indic, fig, ax, **kwargs):
    """
    Plots spectrograms of u(t) (mdr def) or F(t) or q(t) or Q(t)
    Inputs:
        struCase: 'stru' class obj
        data_indic: dict (DOF an) {node:[DOFs]} or list (modal an) [inds]
        ax: matplotlib.pyplot ax obj
    **kwargs (must contain):
        data_type: str, data indicator
        (may contain):
            if 'UDS' (may contain):
                u_type: 'mdr' (default) or 'raw'
            if modal an:
                'modal_inds_type': 'relative' (in order to use MOI´s inds), 'absolute' (phi´s inds, default)
            'vel': bool, default False - In order to use the corresponding data vel
            'SP_Winsize': str, default Dt/20 
            'SP_OvrLapFactor': int, detault 80 (%)
            'SP_WinType': str, default 'hann' - Check supported FFT-Windows in scipy.signal
            'SP_Normalize': bool, default True - In order to normalize the spectrogram
            'y_units': str, 'Hz' or 'rad/s' - freq units
            'x_lims', 'y_lims', list: Plotting lims
            graphs_pack, standard dict for plot customization 
    Returns:
        ax: matplotlib.pyplot ax obj
    """
    data_type = kwargs.get('data_type')
    if 'vel' in kwargs:
        vel = kwargs.get('vel')
    else:
        vel = False
    if 'graphs_pack' in kwargs:
        graphs_pack = kwargs.get('graphs_pack')
    else:
        graphs_pack = handle_graph_info(**kwargs)
    t = struCase.t[struCase.plot_timeInds[0]:struCase.plot_timeInds[1]]
    D_t = t[-1] - t[0]
    if 'SP_Winsize' in kwargs:
        WinSize = int(len(t)/kwargs.get('SP_Winsize'))
    else:
        WinSize = int(len(t)/20) #NOTA: Esto asume un único winsize para el struCase. Ver si agregar info en un dict aparte para más customización
    if 'SP_OvrLapFactor' in kwargs:
        OverLapFactor = kwargs.get('SP_OvrLapFactor')
    else:
        OverLapFactor = 80
        
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

    OverLap = np.round(WinSize*OverLapFactor/100)
    fDef = len(t)/D_t
    y_data = []
    if data_type == 'q' or data_type == 'Q':
        plot_type = 'modal'
        if 'modal_inds_type' in kwargs:
            modal_inds_type = kwargs.get('modal_inds_type')
        else:
            modal_inds_type = 'absolute' #Ampliar esto a mois o relativos
        if type(data_indic) == int:
            data_indic = [data_indic]
        data_indic = handle_modal_inds(struCase, data_indic, **kwargs)
        for loc_ind in data_indic:
            y_data.append(getattr(struCase, data_type)[loc_ind-1,struCase.plot_timeInds[0]:struCase.plot_timeInds[1]])
    elif data_type == 'UDS' or data_type == 'FCS':
        plot_type = 'DOF'
        if 'u_type' in kwargs:
            u_type = kwargs.get('u_type')
        else:
            u_type = 'mdr'
        desired_inds = sim_db.nodeDof2idx(struCase,data_indic)
        original_inds = flatten_values_list(data_indic.values())
        node_labels = label_asoc(data_indic) #OJO. No contempla posibles errores
        for i in range(len(desired_inds)):
            if data_type =='UDS':
                if u_type=='mdr':
                    y_data.append(struCase.u_mdr[desired_inds[i],struCase.plot_timeInds[0]:struCase.plot_timeInds[1]])
                elif u_type=='raw':
                    y_data.append(struCase.u_raw[desired_inds[i],struCase.plot_timeInds[0]:struCase.plot_timeInds[1]])
                else:
                    raise NameError('Wrong u_type def')
            elif data_type =='FCS':
                y_data.append(struCase.aLoad[desired_inds[i],struCase.plot_timeInds[0]:struCase.plot_timeInds[1]])
    else:
        raise NameError('Bad data_type')
    counter = 0
    for y in y_data:
        if vel:
            y=np.gradient(y,t) #NOTA: Agregar al plot que es una velocidad
        F, T, S = signal.spectrogram(y, fDef,window=WinType, noverlap=OverLap,nperseg=WinSize)
        if b_norm:
            loc_m = np.max(S)
            if not loc_m == 0:
                S = abs(S)/loc_m
            else:
                print('Warning: 1/0 found, data_indic: '+str(i))
        if y_units == 'rad/s':
            print('f in rad/s!')
            y_values = F*2*np.pi
        elif y_units == 'Hz':
            print('f in hz!')
            y_values = F
        c = ax.pcolormesh(T,y_values,S,shading = 'auto', cmap='gray_r')
        fig.colorbar(c, ax = ax)
        if plot_type == 'modal':
            ax.set_title('Modo: ' + lst2str(data_indic))
            ax.set_ylabel(graphs_pack['y_label'])
        elif plot_type == 'DOF':
            ax.set_title('Node(s): ' + keys_to_str(data_indic) + ', DOF(s):' + str(original_inds[counter]))
            counter +=1
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

def plt_uxuy(struCase, vsDict, ax, **kwargs):
    '''
    Plots two DOFs or FCS for all struCase nodes for a single t val (one per curve)
    
    Inputs:
        struCase: 'stru' class obj
        vsDict: dict, contains the two DOFs info and the desired time valsn ['DOFs':[indexes],'t_vals':[t]]
        ax: matplotlib.pyplot ax obj
    **kwargs (may contain):
        'data_type', str: 'FCS' or 'UDS' (loads or u-dofs, u-dofs as default)
        if 'UDS' (may contain):
            'u_type', str - 'raw' or 'mdr'
            'deg', bool - False (default) or True (in order to plot rot dofs in deg)            
    Returns:
        ax: matplotlib.pyplot ax obj
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
        inds_t.append(sim_db.search_time(struCase.t,[des_t,0])[0])
    ux_Dict = dof2dofDict(struCase,vsDict['DOFs'][0])
    ux_inds = np.array(sim_db.nodeDof2idx(struCase,ux_Dict))
    uy_Dict = dof2dofDict(struCase,vsDict['DOFs'][1])
    uy_inds = np.array(sim_db.nodeDof2idx(struCase,uy_Dict))
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
    Inputs:
        struCase: 'stru' class obj
        inds: list, list of inds to plot
        ax: matplotlib.pyplot ax obj
    **kwargs:
        (must contain)
        special_mode: str or bool - Flag for special plots (default False), 'max_vs_q'
        loc_data_type:, str, 'E, QW, LW', etc - Data type
        (may contain)
        vel: bool, False(default) or True (in order to compute and plot dy/dt)
        (if special_mode = 'max_vs_q')
        x_vals: list or numpy.array - x values to plot
    Returns:
        ax: matplotlib.pyplot ax obj
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
    if graphs_pack['legend_bool']:
        ax.legend(title=graphs_pack['legend_title'])
    return(ax)

"""
------------------------------------------------------------------------------
HIGH LEVEL PLOT functions
------------------------------------------------------------------------------
"""

def fig_general(struCase, indLIST, **kwargs):
    '''
    Arranges plots of general props as f(t) using some inds (not Python´s)
    Inputs:
        struCase: 'stru' class obj
        indLIST: list, nested [[INDs]]
    **kwargs: 
        (must contain)
        special_mode: str or bool, flag for special plots (default False), 'max_vs_q'
        data_type: nested list of str, 'E, W_load, W_modal' - Data type (for each 0-list)
        (may contain):
            #General:
            fig_save: bool - For saving purp.
                fig_save_opts: dict - Folder, filecode, etc
            sharex: matplotlib.pyplot.subplots() argument - default 'col'
            p_prow: plots per row for the global figure - default 1
            limit_tvals or limit_tinds: list or ndarray - time limits for plotting, values or indexes
            #Plot customization:
                fig_title
                x_label
                y_label
                y_label_setting
                legend_title
    Returns:
        fig: matplotlib.pyplot fig obj
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

def fig_uqs(struCase, DATA_indc, **kwargs):
    '''
    Arranges plots of DOF(nodes) or FC(nodes) (or modes (nodes) or modal ext. loads (nodes)) @t fixed
    
    Inputs:
        struCase: 'stru' class object
        DATA_indc: list of tdataDicts or a single tdataDict {DOF/MODE: [t_instants]}
    **kwargs (must contain):
        data_type: str, data type flag
        (may contain):
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
    if not 'data_type' in kwargs:
        print('Warning - No data_type in kwargs, using DEFAULT') #Acomodar
        kwargs['data_type'] = 'modal_desp' #ACOMODAR Y SET UN DEF
        
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
    
    if kwargs['data_type'] == 'q' or kwargs['data_type'] == 'Q':
        fig_type = 'modal'
    elif kwargs['data_type'] == 'UDS' or kwargs['data_type'] == 'FCS':
        fig_type = 'DOF'
    if type(DATA_indc) == dict: #Para comodidad end-user
        DATA_indc = [DATA_indc]
    n = len(DATA_indc)
        
    fig, axs = plt.subplots(n, p_prow, sharex = sharex)
    if n == 1: #Esto falla si ax no es un iterable (cuando n = 1 es sólo ax, no ax[:])
        axs = plt_tfixed(struCase, DATA_indc[0], axs, **kwargs)
        axs.set_xlabel(graphs_pack['x_label'])
        axs.set_ylabel(graphs_pack['y_label'])
        axs.grid()
    else:
        for ax, exp_DATA_indc in zip(axs, DATA_indc):
            ax = plt_tfixed(struCase,  exp_DATA_indc, ax,**kwargs)
            ax.set_ylabel(graphs_pack['y_label'])
            ax.grid()
        axs[-1].set_xlabel(graphs_pack['x_label'])
    fig.suptitle(graphs_pack['fig_title'])
    if fig_save:
        save_figure(fig,fig_save_opts,**kwargs)
    return(fig)


def fig_uqt(struCase, DATA_indc, **kwargs):
    '''
    Arranges plots of u(t) or load(t) (or q(t) or Q(t))
    
    Inputs:
        struCase: 'stru' class object
        DATA_indc: lst of dofDicts or dofDict {NODE: [DOFs]} (DOF an) or lst of modal_inds or [[fig_1 inds], ...,[fig_n inds]]
    **kwargs (must contain):
            data_type: str, data indicator (UDS, FCS, q Q)
        (may contain):
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
    Returns:
        fig: matplotlib.pyplot figure obj
    '''
    if not 'data_type' in kwargs:
        raise NameError('No data_type in kwargs!')
    
    if 'sharex' in kwargs:
        sharex = kwargs.get('sharex')
    else:
        sharex = "col"
        
    if 'p_prow' in kwargs:
        p_prow = kwargs.get('p_prows')
    else:
        p_prow = 1
    if 'limit_tvals' in kwargs:
        struCase = sim_db.sfti_time(struCase,time_vals=kwargs.get('limit_tvals'))
    elif 'limit_tinds' in kwargs:
        struCase = sim_db.sfti_time(struCase,indexes = kwargs.get('limit_tinds'))
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
    
    if kwargs['data_type'] == 'FCS' or kwargs['data_type'] == 'UDS':
        fig_type = 'DOF'
        if type(DATA_indc) == dict: #Para comodidad end-user
            DATA_indc = [DATA_indc]
    elif kwargs['data_type'] == 'q' or kwargs['data_type'] == 'Q':
        fig_type = 'modal'
    n = len(DATA_indc)
    
    fig, axs = plt.subplots(n, p_prow, sharex = sharex)
    if n == 1: #Esto falla si ax no es un iterable (cuando n = 1 es sólo ax, no ax[:])
        if fig_type == 'modal':
            if type(DATA_indc[0]) == int:
                DATA_indc = [DATA_indc]
        axs = plt_temp_ev(struCase, DATA_indc[0], axs, **kwargs)
        axs.set_xlabel(graphs_pack['x_label'])
        axs.set_ylabel(graphs_pack['y_label'])
        axs.set_aspect(aspect_ratio)
        axs.grid()
    else:
        for ax, exp_DATA_indc in zip(axs, DATA_indc):
            ax = plt_temp_ev(struCase, exp_DATA_indc, ax,**kwargs)
            ax.set_ylabel(graphs_pack['y_label'])
            ax.set_aspect(aspect_ratio)
            ax.grid()
        axs[-1].set_xlabel(graphs_pack['x_label'])
    fig.suptitle(graphs_pack['fig_title'])
    if fig_save:
        save_figure(fig,fig_save_opts,**kwargs)
    return(fig)

def fig_FFT(struCase, DATA_indic, **kwargs):
    '''
    Arranges plots of FFT(u(t)) or FFT(load(t)) or FFT(q(t)) or FFT(Q(t))
    
    Inputs:
        struCase: 'stru' class obj
        DATA_indic: lisf of dofDicts or a single one (DOF an) {'NODE':[DOFs]} or list of modal inds list [[modal_ind]] (modal an)
    **kwargs (must contain):
        data_type: str, data indicator (UDS, FCS, q, Q)
    (may contain):
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
    Returns:
        fig: matplotlib.pyplot figure obj
    '''
    if not 'data_type' in kwargs:
        raise NameError('Warning: No data_type in kwargs')
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
        struCase = sim_db.sfti_time(struCase,time_vals=kwargs.get('limit_tvals'))
    elif 'limit_tinds' in kwargs:
        struCase = sim_db.sfti_time(struCase,indexes=kwargs.get('limit_tinds'))
    else:
        pass
    if 'fig_save' in kwargs:
        fig_save = kwargs.get('fig_save')
        if 'fig_save_opts' in kwargs:
            fig_save_opts = kwargs.get('fig_save_opts')
        else:
            fig_save_opts = {}
    else:
        fig_save = False
    graphs_pack = handle_graph_info(**kwargs)

    if type(DATA_indic) == dict: #Para comodidad end-user
        DATA_indic = [DATA_indic]            
    n = len(DATA_indic)
    
    fig, axs = plt.subplots(n,p_prow, sharex=sharex)
    if n == 1: #Esto falla si ax no es un iterable (cuando n = 1 es sólo ax, no ax[:])
        axs = plt_FFT(struCase, DATA_indic[0], axs, **kwargs)
        axs.set_xlabel(graphs_pack['x_label'])
        axs.set_ylabel(graphs_pack['y_label'])
        axs.grid()
    else:
        for ax, exp_DATA_indic in zip(axs, DATA_indic):
            ax = plt_FFT(struCase, exp_DATA_indic, ax,**kwargs)
            ax.set_ylabel(graphs_pack['y_label'])
            ax.grid()
        axs[-1].set_xlabel(graphs_pack['x_label'])
    fig.suptitle(graphs_pack['fig_title'])
    if fig_save:
        save_figure(fig,fig_save_opts,**kwargs)
    return(fig)

def fig_spect(struCase, DATA_indic, **kwargs):
    '''
    Arranges plots of Spectrogram(u(t)) or Spectrogram(load(t)) or Spectrogram(q(t)) or Spectrogram(Q(t))
    
    Inputs:
        struCase: 'stru' class obj
        DATA_indic: list of dofDics {'NODE':[DOFs]} (or a single one) or list of modal inds lists [[modal inds]]
    **kwargs (must contain):
        data_type: str, data type indicator
    (may contain):
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
        general kwargs for downstream funcs
    Returns:
        fig: matplotlib.pyplot figure obj
    '''
    
    if 'data_type' in kwargs:
        data_type = kwargs.get('data_type')
    else:
        raise NameError('Warning: No data_type in kwargs')
        
    if 'sharex' in kwargs:
        sharex = kwargs.get('sharex')
    else:
        sharex = "col"
        
    if 'p_prow' in kwargs:
        p_prow = kwargs.get('p_prows')
    else:
        p_prow = 1

    if 'limit_tvals' in kwargs:
        struCase = sim_db.sfti_time(struCase,time_vals=kwargs.get('limit_tvals'))
    elif 'limit_tinds' in kwargs:
        struCase = sim_db.sfti_time(struCase,indexes=kwargs.get('limit_tinds'))
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
    
    if type(DATA_indic) == dict: #Para comodidad end-user
        DATA_indic = [DATA_indic]
    n = len(DATA_indic)
    fig, axs = plt.subplots(n,p_prow, sharex=sharex)
    if n == 1: #Esto falla si ax no es un iterable (cuando n = 1 es sólo ax, no ax[:])
        axs = plt_spectr(struCase, DATA_indic[0], fig, axs, **kwargs)
        axs.set_xlabel(graphs_pack['x_label'])
        axs.set_ylabel(graphs_pack['y_label'])
        axs.grid()
    else:
        for ax, exp_DATA_indic in zip(axs, DATA_indic):
            ax = plt_spectr(struCase, exp_DATA_indic, fig, ax,**kwargs)
            ax.set_ylabel(graphs_pack['y_label'])
            ax.grid()
        axs[-1].set_xlabel(graphs_pack['x_label'])
    if fig_save:
        save_figure(fig,fig_save_opts,**kwargs)
    return(fig)

def fig_PP(struCase, DATA_indic, **kwargs):
    '''
    Arranges plots of PP(u(t)) or PP(load(t)) or PP(q(t)) or PP(Q(t))
    
    Inputs:
        struCase: 'stru' class obj
        DATA_indic: list of dofDicts {'NODE':[DOFs]} (or a single one, DOF an) or a list of modal inds lists [[inds]] (modal_an)
    **kwargs (must contain):
        data_type: str, data type indicator (UDS, FCS, q, Q)
        (may contain):
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
            general kwargs for downstream funcs
    '''
    
    if 'data_type' in kwargs:
        data_type = kwargs.get('data_type')
    else:
        raise NameError('Warning: No data_type in kwargs')
    if 'sharex' in kwargs:
        sharex = kwargs.get('sharex')
    else:
        sharex = "col"
        
    if 'p_prow' in kwargs:
        p_prow = kwargs.get('p_prows')
    else:
        p_prow = 1
        
    if 'limit_tvals' in kwargs:
        struCase = sim_db.sfti_time(struCase,time_vals=kwargs.get('limit_tvals'))
    elif 'limit_tinds' in kwargs:
        struCase = sim_db.sfti_time(struCase,indexes=kwargs.get('limit_tinds'))
    else:
        pass
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
    
    if type(DATA_indic) == dict: #Para comodidad end-user
        DATA_indic = [DATA_indic]

    n = len(DATA_indic)
    
    fig, axs = plt.subplots(n, p_prow, sharex = sharex)
    if n == 1: #Esto falla si ax no es un iterable (cuando n = 1 es sólo ax, no ax[:])
        axs = plt_PP(struCase, DATA_indic[0], axs, **kwargs)
        axs.set_xlabel(graphs_pack['x_label'])
        axs.set_ylabel(graphs_pack['y_label'])
        axs.set_aspect(aspect_ratio)
        axs.grid()
    else:
        for ax, exp_DATA_indic in zip(axs, DATA_indic):
            ax = plt_PP(struCase, exp_DATA_indic, ax,**kwargs)
            ax.set_ylabel(graphs_pack['y_label'])
            ax.set_aspect(aspect_ratio)
            ax.grid()
        axs[-1].set_xlabel(graphs_pack['x_label'])
    fig.suptitle(graphs_pack['fig_title'])
    if fig_save:
        save_figure(fig,fig_save_opts,**kwargs)
    return(fig)

def fig_phi(struCase, modedofLIST, **kwargs):
    '''
    Arranges plots of modal shapes
    
    Inputs:
        struCase: 'stru' class object
        mdofLIST: list or dict or a single mdofLIST [{'mode_ind':[DOFs]}]
    **kwargs (may contain):
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
    Returns:
        fig: matplotlib.pyplot figure obj
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
    Inputs:
        lst: list - items
    **kwargs (may contain):
        spacer: str - change prnt style
    Returns:
        ret: str
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

def dof2dofDict(struCase, dof):
    '''
    Generates a dofDict for a desired DOF using all the available nodes in struCase
    Inputs: 
        struCase: 'stru' class obj
        dof: int, desired dof to look for
    Returns:
        dofdict: dict
    '''
    dofdict = {}
    for a in struCase.nodes:
        dofdict[str(a)]=[dof]
    return(dofdict)

def tdof2dofDict(struCase, tdof_dict):
    '''
    Generates a dofDict for a desired DOF using all the available nodes in struCase
    Inputs:
        struCase: 'stru' class obj
        tdofDict: dict, contains the desired dof as key and some time instants {DOF:[t_vals]}
    Returns:
        dofDict: dict 
    '''
    dofdict = {}
    loc_lst = []
    for j in list(tdof_dict.keys()):
            loc_lst.append(int(j))
    for a in struCase.nodes:
        dofdict[str(a)]=loc_lst
    return(dofdict)

def nodes2labels(struCase, **kwargs):
    '''
    Generates a list of strings from struCase.nodes
    Inputs:
        struCase: 'stru' class obj
    Returns:
        labels, list of str
    '''
    loc_lst = []
    for loc_node in struCase.nodes:
        loc_lst.append(str(loc_node))
    return(loc_lst)
    
def tldxs(struCase, desired_t):
    """
    Generates a list of indexes for the closest struCase.t time values in t arg (list, floats)
    Inputs:
        struCase: 'stru' class obj
        desired_t: float, a desired time [s]
    Returns:
        closest struCase.t index: int
    """
    des_ind = []
    for loc_t in desired_t:
        des_ind.append((np.abs(struCase.t-loc_t)).argmin())

    return(des_ind)

def keys_to_str(dct_keys):
    """
    Simple tool for str generation from dict_keys data
    
    Inputs: 
        dct_keys: dict.keys
    Returns: 
        word: str
    """
    word = ''
    for loc_nam in list(dct_keys):
        word+=loc_nam
        word+='-'
    return(word[:-1])

def flatten_values_list(dct_values):
    """
    Simple tool for fusion nested lists (coming from dicts_values)
    
    Inputs:
        dct_values: dict.values
    Returns:
        flat_values, list (1D)
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
    
    Inputs:
        dct: dict {key:list}
    Returns:
        label_list: list of str
    """
    label_list = []
    for i in dct:
        label_list.append(len(dct[i])*[str(i)])
    label_list = flatten_values_list(label_list)
    return(label_list)
    
def handle_graph_info(**kwargs):
    """
    Prepares labels, titles and general things for plots
    
    Inputs:
    **kwargs(may contain):
        general graph info kwargs
    Returns:
        graphs_info: dict
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
    
    if 'legend_bool' in kwargs:
        graphs_pack['legend_bool'] = kwargs.get('legend_bool')
    else:
        graphs_pack['legend_bool'] = False
    return(graphs_pack)

def save_figure(fig, opts,**kwargs):
    '''
    Save plots in desired dir
    
    Inputs:
        fig: matplotlib fig obj - Figure
        opts: dict
            Options: {'folder','filecode','filename','dpi','fig_format', 'bbox_inches','fig_dpi'}
    **kwargs (may contain):
        close, bool - Close fig
    Returns:
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

def handle_modal_inds(struCase, modal_inds, **kwargs):
    '''
    Looks for the real modal index in struCase.moi, for labeling purposes
    
    Inputs:
        struCase: 'stru' class obj
        modal_inds: list, modal inds
    **kwargs (may contain):
        modal_inds_type: str, relative or absolute (MOIs or full)
    Returns:
        modal_inds: list, updated modal inds
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

def hide_axes(ax, **kwargs):
    '''
    Auto ax hide
    Inputs:
        ax, matplotlib.pyplot ax obj (or a set)
    **kwargs (may contain):
        hide_x: bool, default False
        hide_y: bool, default False
        inds: list, desired ax inds to hide, default False
    Returns:
        ax, matplotlib.pyplot ax obj (or a set)  
    '''
    if 'hide_x' in kwargs:
        hide_x = kwargs.get('hide_x')
    else:
        hide_x = False
    if 'hide_y' in kwargs:
        hide_y = kwargs.get('hide_y')
    else:
        hide_y = False
    if 'inds' in kwargs:
        inds = kwargs.get('inds')
    else:
        inds = False
    
    if not inds:
        try:
            for loc_ax in ax[:-1]:
                if hide_x:
                    loc_ax.axes.xaxis.set_visible(False)
                if hide_y:
                    loc_ax.axes.yaxis.set_visible(False)
        except:
            try:
                if hide_x:
                    ax.axes.xaxis.set_visible(False)
                if hide_y:
                    ax.axes.yaxis.set_visible(False)
            except:
                raise NameError('Warning: Error handling ax')
    return(ax)
                        
def do_grid_and_labels(ax, **kwargs):
    '''
    Auto grid and legend
    Inputs:
        ax, matplotlib.pyplot ax obj (or a set)
    **kwargs (may contain):
        grid: bool, default False
        legend: bool, default False
        y_label: str, default False
    Returns:
        ax, matplotlib.pyplot ax obj (or a set)
    '''
    if 'grid' in kwargs:
        grid = kwargs.get('grid')
    else:
        grid = True
    if 'legend' in kwargs:
        legend = kwargs.get('legend')
    else:
        legend = True
    if 'y_label' in kwargs:
        y_label_bl = True
    else:
        y_label_bl = False
        
    try:
        for loc_ax in ax:
            if grid:
                loc_ax.grid()
            if legend:
                loc_ax.legend()
            if y_label_bl:
                loc_ax.set_ylabel(kwargs.get('y_label'))
    except:
        try:
            if grid:
                ax.grid()
            if legend:
                ax.legend()
            if y_label_bl:
                ax.set_ylabel(kwargs.get('y_label'))
        except:
            raise NameError('Warning: Error handling ax')
    return(ax)

def FFT_labl_gen(data_type, i, original_inds, node_labels):
    '''
    Gen a custom label for FFT plots
    
    Inputs:
        data_type: str, UDS or FCS or q or Q
        i: int, counter
        original_inds: list, only for DOF plots
        node_labels: list, only for DOF plots
    Returns:
        labl: str, final label

    '''
    if data_type == 'UDS' or data_type == 'FCS':
        labl = node_labels[i] +' - DOF: ' + str(original_inds[i])
    else:
        labl = data_type + '-' + str(i)
    return labl

def gen_aloneinds_prop(inds):
    '''
    Generates a nested list for indsLIST kwarg to use in fig_general, using a single list of inds.
    Inputs:
        inds: list, desired inds
    Returns:
        indsLIST: list, nested inds list
    '''
    indsLIST = []
    for ind in inds:
        indsLIST.append([ind])
    return(indsLIST)

def general_envs(struCase,**kwargs):
    '''
    Gens envelopes - Plots (max(f[i])) vs some inds
    
    Inputs:
        struCase: 'stru' class obj
    **kwargs (may contain):
        data_type: str - desired data
        t_window: list - desired inds window [t_i, t_f]
        new_name: str - new attr name
    Returns:
        struCase: 'stru' class obj, updated
    '''
    if 'data_type' in kwargs:
        data_type = kwargs.get('data_type')
    else:
        raise NameError('Data type str error')
        
    if 't_window' in kwargs:
        t_window = kwargs.get('t_window')
    else:
        t_window = [0,None]
    if 'new_name' in kwargs:
        new_name = kwargs.get('new_name')
    else:
        new_name = data_type + '_envs'
    prop = getattr(struCase,data_type)[:,t_window[0]:t_window[1]]
    
    new_prop = []
    for i in range(prop.shape[0]):
        lmin,lmax = sim_db.hl_envelopes_idx(prop[i])
        raw_data_max = prop[i,lmax]
        raw_data_min = prop[i,lmin]
        new_prop.append([max(raw_data_max),max(abs(raw_data_min))])
    
    new_prop = np.array(new_prop)
    setattr(struCase, new_name, new_prop)
    return struCase

def add_ticks(ax,**kwargs):
    '''
    Add ticks to an ax obj
    
    Inputs:
        ax, matplotlib.pyplot ax obj
    **kwargs (may contain):
        x: array or list, data ticks steps
            xt_type: str, 'none' or 'double' (number of ticks)
            xt_major_step: int, ticks steps
    Returns:
        ax, matplotlib.pyplot ax obj
    '''
    if 'x' in kwargs:
        x = kwargs.get('x')
        if 'xt_type' in kwargs:
            x_ticks = kwargs.get('xt_type')
        else:
            x_ticks = 'double'
        if 'xt_major_step' in kwargs:
            xt_major_step = kwargs.get('xt_major_step')
        else:
            xt_major_step = 2
    
        if not x_ticks == 'none':
            major_ticks  = np.arange(x[0],x[-1],abs(x[xt_major_step]-x[0]))
            ax.set_xticks(major_ticks)
            if x_ticks == 'double':
                minor_ticks = np.arange(x[0],x[-1], abs(x[1]-x[0])/2)
                ax.set_xticks(minor_ticks, minor=True)
    return(ax)

def gen_single_prop(desired_shape, prop):
    '''
    Generates a list for data_type kwarg to use in fig_general, using a single struCase.attr.
    
    Inputs:
        desired_shape: list, inds list or nested inds list
        prop: str, desired attr name
    Returns:
        propt_lst: list, fig_general data_type kwarg
    '''
    
    prop_lst = []
    
    for frst_lev in desired_shape:
        if type(frst_lev) == list:
            loc_prop_lst = []
            for secnd_lev in frst_lev:
                loc_prop_lst.append(prop)
            prop_lst.append(loc_prop_lst)
        else:
            prop_lst.append(prop)
    return(prop_lst)


"""
------------------------------------------------------------------------------
runing things
------------------------------------------------------------------------------
"""

if __name__ == '__main__':
    #Examples
    pass    