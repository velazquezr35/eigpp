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
from sim_db import nodeDof2idx
from scipy.fft import fft, fftfreq
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
    
    Inputs: struCase is a Stru Class Obj
            dof_Dict is a dict (o lista, ver cual dejar), {node:[DOFs]}
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

    t=struCase.t
    desired_inds = nodeDof2idx(struCase, dofDict)
    original_inds = flatten_values_list(dofDict.values())
    node_labels = label_asoc(dofDict) #OJO. No contempla posibles errores (q pida algo que no tengo) y esto daría problemas. Parece no importar.
    for i in range(len(desired_inds)):
        if u_type=='avr':
            u = struCase.u_avr[desired_inds[i],:]
        elif u_type=='raw':
            u = struCase.u_raw[desired_inds[i],:]
        else:
            print('Warning: Bad u_type def')
            
        if vel:
            u=np.gradient(u,struCase.t) #NOTA: Agregar al plot que es una velocidad
        ax.plot(t,u, label= node_labels[i] +' - DOF: ' + str(original_inds[i])) #NOTA: Creo que no es necesario el transpose, lo detecta sólo.
    ax.legend(title='Node(s): ' + keys_to_str(dofDict)) #NOTA: ¿Quiero esta info como titulo?
    return(ax) #NOTA: ¿Necesito hacer el return? Quizá para actualizar
    
#Plot q data

def plt_qt(struCase, modeList, ax, **kwargs):
    """
    Plot modal coordinates as f(t).
    
    Inputs: struCase is a Stru Class Obj
            dof_lst is a dict (o lista, ver cual dejar), {node:[DOFs]}
            ax is a matplotlib.pyplot Axes obj
            kwargs: 'vel': False (default) or True (in order to calculate and plot modal velocities)
                    
    """
    if 'vel' in kwargs:
        vel = kwargs.get('vel')
    else:
        vel = False

    t=struCase.t
    for loc_ind in modeList:
        u = struCase.q[loc_ind,:]
        if vel:
            u=np.gradient(u,struCase.t)
        ax.plot(t,u, label=str(loc_ind)) #NOTA: Creo que no es necesario el transpose, lo detecta sólo.

    ax.legend(title='COMPLETAME')
        
    return(ax) #NOTA: ¿Necesito hacer el return? Quizá para actualizar
    

def plt_us(struCase, t_lst,ax,**kwargs):
    
    """
    Plot all DOFs in particular instants of time, u_avr as default.
    
    Inputs: struCase is a Stru Class Obj
            t_lst is a Stru.t index list
            ax is a matplotlib.pyplot Axes obj
            kwargs: 'u_type': raw or avr (default)
                    'vel': False (default) or True (in order to calculate and plot velocities)
                    
    """
    
    print('Coming soon...')
    
    return()

def plt_qs(struCase, t_lst,ax,**kwargs):
    
    """
    Plot all DOFs in particular instants of time, u_avr as default.
    
    Inputs: struCase is a Stru Class Obj
            t_lst is a Stru.t index list
            ax is a matplotlib.pyplot Axes obj
            kwargs: 'u_type': raw or avr (default)
                    'vel': False (default) or True (in order to calculate and plot velocities)
                    
    """
    print('Coming soon...')
    return()



#Fourier 

def plt_uFFT(struCase, dofDict, ax, **kwargs):
    """
    Plots FFT for u signal
    Inputs: struCase is a Stru Class Obj
            dof_lst is a dict (o lista, ver cual dejar), {node:[DOFs]}
            ax is a matplotlib.pyplot Axes obj
            kwargs: 'u_type': raw or avr (default)
    """
    if 'u_type' in kwargs:
        u_type = kwargs.get('u_type')
    else:
        u_type = 'avr'
    
    #NOTA: Nuevamente, ciclar sobre los DOFs de interés a sacar la FFT
    
    if u_type == 'avr':
        u = struCase.u_avr#[CUAL]
    elif u_type == 'raw':
        u = struCase.u_raw#[CUAL]
    else:
        print('Warning: Bad u_type def')
    
    loc_tf = max(struCase.t)
    loc_ti = min(struCase.t)
    loc_Deltat = loc_tf-loc_ti
    Dfreq = 1/loc_Deltat
    
    #REVISAR CLASE GRABADA Y ACOMODAR FACTORES
    '''
    N = 1000
    factor= 1000
    T = 1.0 /1000
    yf = fft(y)
    xf = fftfreq(N, T)[:N//2]
    fig4,ax4 = plt.subplots()
    ax4.plot(xf/factor, np.abs(yf[0:N//2])/max(abs(yf)))
    #multiplicar xf por len(t)-1
    Probar distintos delta t y valores inicio / final
    input t
    ti=min(t)
    tf=max(t)
    Dfreq=1/(tf-ti)
    freq=0:Dfreq:Dfreq*(len(t)-1)'''
    return()
    

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
            u = struCase.u_avr[loc_ind,:]
        elif u_type=='raw':
            u = struCase.u_raw[loc_ind,:]
        else:
            print('Warning: Bad u_type def')
        
        du = np.gradient(u,struCase.t)
        ax.plot(u,du, label=str(loc_ind)) #NOTA: Creo que no es necesario el transpose, lo detecta sólo.
        
    return(ax) #NOTA: ¿Necesito hacer el return? Quizá para actualizar

def plt_qPP(struCase, modeList,ax,**kwargs):
    """
    Plots phase-plane portraits, u vs du/dt
    Inputs: struCase is a Stru Class Obj
            modeList
            ax is a matplotlib.pyplot Axes obj
            kwargs: 'u_type': raw or avr (default)
    """

    for loc_ind in modeList:
        #NOTA: Esto se puede mejorar tomando u = todos y luego plot(u[desired])
        dq = np.gradient(struCase.q[loc_ind,:],struCase.t)
        ax.plot(struCase.q[loc_ind,:],dq, label=str(loc_ind)) #NOTA: Creo que no es necesario el transpose, lo detecta sólo.
        
    return(ax) #NOTA: ¿Necesito hacer el return? Quizá para actualizar

#Spectrogram

def plt_spectr(struCase, dofDict, ax,**kwargs):
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
        
    """ #COMPLETAR CON LO DE LA REUNIÓN
        WinSize = 1000
    #Tomar def. igual delta t /20
    OvrLap = 10
    dy = np.diff(y)
    
    fs = len(f_time)/(f_time[-1]-f_time[0])
    
    F, T, S = signal.spectrogram(y,fs,window='hann',nfft=WinSize, noverlap=np.round((WinSize*OvrLap)/100))
    loc_m = 0
    for i in range(len(S)):
        if max(S[i,:]) > loc_m:
            loc_m = max(S[i,:])
    
    S = abs(S)/loc_m
    fig, ax = plt.subplots()
    
    c = plt.pcolormesh(T,F*2*np.pi,S,shading = 'auto', cmap='gray_r')
    ax.grid(linewidth=0.5)
    # ax.set_ylim(2,4)
    ax.set_xlabel('Tiempo [s]')
    ax.set_ylabel('w [rad/s]')
    fig.colorbar(c, ax=ax)
    """
    
        
    return(ax)
"""
------------------------------------------------------------------------------
HIGH LEVEL PLOT functions
------------------------------------------------------------------------------
"""

def fig_ut(struCase, dofLIST, **kwargs):
    '''
    Arranges plots of u(t)
    
    struCase:   stru class object
    dofLIST:    list of dofDicts or dofDict {NODE: [DOFs]} 
    kwargs: may contain
        #General:
        sharex: matplotlib.pyplot.subplots() argument - default 'col'
        p_prow: plots per row for the global figure - default 1
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
    
    if type(dofLIST) == dict: #Para comodidad end-user
        dofLIST = [dofLIST]

    n = len(dofLIST)
    
    fig, axs = plt.subplots(n, p_prow, sharex = sharex)
    
    if n == 1: #Esto falla si ax no es un iterable (cuando n = 1 es sólo ax, no ax[:])
        axs = plt_ut(struCase, dofLIST[0], axs, **kwargs)
        axs.set_xlabel(graphs_pack['x_label'])
        axs.set_ylabel(graphs_pack['y_label'])
        axs.grid()
    else:
        for ax, dof_dict in zip(axs, dofLIST):
            ax = plt_ut(struCase, dof_dict, ax)
            ax.set_xlabel(graphs_pack['x_label'])
            ax.set_ylabel(graphs_pack['y_label'])
            ax.grid()
    fig.suptitle(graphs_pack['fig_title'])
    return(fig)

#TEST VERSION. Creo que no será necesaria, y también creo que no hace exactamente lo solicitado en el Issue con su pseudocódigo
#IDEM PARA fig_one
def fig_ut_col(struCase, dofDict, **kwargs):
    """
    Arranges plots of u(t)
    
    struCase:   stru class object
    dofDict:    dict (o lista, ver cual dejar), {node:[DOFs]}
    kwargs: may contain
        sharex: matplotlib.pyplot.subplots() argument - default 'col'
        p_prow: plots per row for the global figure - default 1
    """
    if 'sharex' in kwargs:
        sharex = kwargs.get('sharex')
    else:
        sharex = "col"
        
    if 'p_prow' in kwargs:
        p_prow = kwargs.get('p_prows')
    else:
        p_prow = 1
        
    graphs_pack = handle_graph_info(**kwargs)

    #n subplots for n lists in dofLIST
    if len(dofDict) > 1:
        raise ValueError('Warning: More than 1 node(?) in dofDict. Use fig_ut instead')
    else:
        
        n = len(flatten_values_list(dofDict.values()))
        key = list(dofDict.keys())[0]
        fig, axs = plt.subplots(n, p_prow, sharex = sharex)
        
        if n == 1: #El zip falla si axs no es iterable, nuevamente.
            axs = plt_ut(struCase, dofDict, axs,**kwargs)
            axs.set_xlabel(graphs_pack['x_label'])
            axs.set_ylabel(graphs_pack['y_label'])
            axs.grid()
        else:
            for ax, i in zip(axs, range(n)):
                ax = plt_ut(struCase,{key:[flatten_values_list(dofDict.values())[i]]},ax,**kwargs)
                ax.set_xlabel(graphs_pack['x_label'])
                ax.set_ylabel(graphs_pack['y_label'])
                ax.grid()
        fig.suptitle(graphs_pack['fig_title'])
        return(fig)

def fig_qt():
    print('Coming soon...')
    pass

def fig_u_FFT():
    print('Coming soon...')
    pass

def fig_ut_spect():
    print('Coming soon...')
    pass

def fig_ut_vt_pp():
    print('Coming soon...')
    pass


"""
------------------------------------------------------------------------------
GENERAL TOOLS functions
------------------------------------------------------------------------------
"""
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
    
    
    
    
###################################################################################################################
###################################################################################################################

"""
------------------------------------------------------------------------------
OLD functions PRE REUNIÓN
------------------------------------------------------------------------------
"""

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

def spect_plot(t, y, **kwargs):
    
    WinSize = 1000
    #Tomar def. igual delta t /20
    OvrLap = 10
    dy = np.diff(y)
    
    fs = len(f_time)/(f_time[-1]-f_time[0])
    
    F, T, S = signal.spectrogram(y,fs,window='hann',nfft=WinSize, noverlap=np.round((WinSize*OvrLap)/100))
    loc_m = 0
    for i in range(len(S)):
        if max(S[i,:]) > loc_m:
            loc_m = max(S[i,:])
    
    S = abs(S)/loc_m
    fig, ax = plt.subplots()
    
    c = plt.pcolormesh(T,F*2*np.pi,S,shading = 'auto', cmap='gray_r')
    ax.grid(linewidth=0.5)
    # ax.set_ylim(2,4)
    ax.set_xlabel('Tiempo [s]')
    ax.set_ylabel('w [rad/s]')
    fig.colorbar(c, ax=ax)
    return('None')

def env_plot(t,y,**kwargs):
    fig2,ax2 = plt.subplots()
    ax2.plot(f_time, y, label='signal')
    high_idx, low_idx = hl_envelopes_idx(y)
    plt.plot(f_time[high_idx], y[high_idx], 'r', label='low')
    plt.plot(f_time[low_idx], y[low_idx], 'g', label='high')
    ax2.legend()
    ax2.set_xlabel('Tiempo [s]')
    ax2.set_ylabel('f(t)')
    ax2.grid()
    
def FFT_plot(t,y,**kwargs):

    N = 1000
    factor= 1000
    T = 1.0 /1000
    yf = fft(y)
    xf = fftfreq(N, T)[:N//2]
    fig4,ax4 = plt.subplots()
    ax4.plot(xf/factor, np.abs(yf[0:N//2])/max(abs(yf)))
    #multiplicar xf por len(t)-1
    # ax4.set_xlim(-10,100)
    ax4.set_xlabel('Eje x')
    ax4.set_ylabel('Eje y')
    ax4.set_title('Intervalo t1 to t2')
    ax4.grid()
    
    '''
    Probar distintos delta t y valores inicio / final
    input t
    ti=min(t)
    tf=max(t)
    Dfreq=1/(tf-ti)
    freq=0:Dfreq:Dfreq*(len(t)-1)'''
    
"""
------------------------------------------------------------------------------
runing things
------------------------------------------------------------------------------
"""

if __name__ == '__main__':
    #Examples
    f_time = np.linspace(0,1000,2000)
    y = np.sin(0.3*f_time) + np.sin(0.7*f_time)
    
    spect_plot(f_time,y)
    env_plot(f_time,y)
    FFT_plot(f_time,y)