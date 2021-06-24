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
            dof_lst is a dict (o lista, ver cual dejar), {node:[DOFs]}
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
    for loc_ind in desired_inds:
        if u_type=='avr':
            u = struCase.u_avr[loc_ind,:]
        elif u_type=='raw':
            u = struCase.u_raw[loc_ind,:]
        else:
            print('Warning: Bad u_type def')
            
        if vel:
            u=np.gradient(u,struCase.t)
        ax.plot(t,u, label=str(loc_ind)) #NOTA: Creo que no es necesario el transpose, lo detecta sólo.
        
    return(ax) #NOTA: ¿Necesito hacer el return? Quizá para actualizar
    #NOTA: Ver si directamente plotear las leyendas acá o devolverlas y agregarlas luego manualmente (no recomiendo). Creo que lo mejor es plotearlas acá directo.
    
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
        
    return(ax) #NOTA: ¿Necesito hacer el return? Quizá para actualizar
    #NOTA: Ver si directamente plotear las leyendas acá o devolverlas y agregarlas luego manualmente (no recomiendo). Creo que lo mejor es plotearlas acá directo.
    



def plt_us(struCase, t_lst,ax,**kwargs):
    
    """
    Plot all DOFs in particular instants of time, u_avr as default.
    
    Inputs: struCase is a Stru Class Obj
            t_lst is a Stru.t index list
            ax is a matplotlib.pyplot Axes obj
            kwargs: 'u_type': raw or avr (default)
                    'vel': False (default) or True (in order to calculate and plot velocities)
                    
    """
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
    dofLIST:    list (level 3) of lists (level 2 - called dofList) of lists (level 1 - pairs [node,DoF])
                e.g.:   dofLIST = [ dofList1, dofList2, dofList3 ]
                            dofList1 = [ [node1,3] ]
                            dofList2 = [ [node1,3], [node1,1] ]
                            dofList3 = [ [node2,1], [node3,1] ]
                            node1, node2 and node3 being integers (external labels of nodes included in struCase.nodes)
                        the figure will have 3 axes arranged on a column
                            top axes will have one curve (time evolution of the 3rd DoF of node node1)
                            mid axes will have two curves (time evolution of the 3rd and the 1st DoF of node node1)
                            bot axes will have two curves (time evolution of the 1st DoF of node node2 and node3)
    kwargs: may contain
        sharex: matplotlib.pyplot.subplots() argument - default 'col'
    '''
    
    if 'sharex' in kwargs:
        sharex = kwargs.get('sharex')
    else:
        sharex = "col"
    
    
    fig, axs = plt.subplots( len(dofLIST), 1, sharex=sharex ) # esta llamada podría generalizarse con más especificaciones que vengan en kwargs, pero creo que esta es la única interesante para esta función - tal vez para otras puede ser interesante usar otros argumentos - lo veremos cuando lo estemos usando
    for ax,lst in zip(axs,dofLIST):
        ax = plt_ut( struCase, lst, ax, **kwargs)
    
    return(fig)

def fig_qt():
    pass

def fig_u_FFT():
    pass

def fig_ut_spect():
    pass

def fig_ut_vt_pp():
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
    """
    des_ind = []
    for loc_t in desired_t:
        des_ind.append((np.abs(struCase.t-loc_t)).argmin())

    return(des_ind)



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