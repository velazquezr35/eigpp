"""
Created on 2021 05 29
@author: Mauro Maza
Departamento de Estructuras
Facultad de Ciencias Exactas, Físicas y Naturales
Universidad Nacioal de Córdoba
Córdoba, Argentina

Exampli Gratia
    script showing how to use eigpp to analize data
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
from sim_db import sim
from eigpp  import epp
import plotter

"""
------------------------------------------------------------------------------
functions
------------------------------------------------------------------------------
"""
def basic_general_plots(case1, **kwargs):
    '''
    General plots...
    
    inputs:
        case1, sim class obj
    
    '''
    #Valores cualquiera
    desired = {'200000':[2]}
    modes = [0]
    modes_2 = [1]
    

    plotter.fig_us(case1.stru, {'2':[100]})
    plotter.fig_qs(case1.stru,{'1':[150]})
    plotter.fig_uxuy(case1.stru,{'DOFs':[2,3],'t_vals':[150,200]}, fig_title = 'Simple versus plot')
    plotter.fig_q_FFT(case1.stru,modes,fig_title = 'Simple FFT(q) plot')
    plotter.fig_qt_vt_pp(case1.stru,modes,fig_title='Simple comb. qt, vqt, phase p.')
    plotter.fig_ut(case1.stru, [desired,desired], fig_title='Simple u(t) plot',x_label ='t',y_label='u')
    plotter.fig_u_FFT(case1.stru, desired, fig_title='Simple FFT(u) plot', limit_tinds=[0,-1])
    plotter.fig_ut_vt_pp(case1.stru,desired, fig_title='Simple comb. ut,vt, phase p.')
    plotter.fig_u_spect(case1.stru,desired, fig_title = 'Simple spectrogram', x_label = 't', y_label='f')
    plotter.fig_qt(case1.stru, [modes, modes_2],fig_title='Simple q(t) plot')
    plotter.fig_q_spect(case1.stru,modes, fig_title='Simple q spectrogram')
    return('Done!')

"""
------------------------------------------------------------------------------
runing things
------------------------------------------------------------------------------
"""

# create simulation data object and populate minimum structural info
case1=sim() # this will contain data corresponding to a particular simulation (a "study case") under analysis
case1.fName='eg' # output binary file name (without extension)
case1.stru.name='str ex raw data'
case1.stru.nodes=[200001] # nodes of interest - only information relative to these nodes is going to be read
case1.stru.struRdOpt = 'raw' # set flag for reading ASCII file - default 'bin'
case1.stru.p11FN='pcolgante.@1' # binary *.p11 file name (without extension - Simpact output) from wich extract generalized displacements
case1.stru.rsnDe='pcolgante' # ASCII *.rsn file name (without extension) - Delta output
case1.stru.loadsFN = 'AeroFcsOnStruc' #ASCII *.dat file name - Loads
case1.stru.loadRdOpt = 'non'
case1.stru.struEigOpt = False
case1.stru.loadEigOpt = False
#data dir: use / instead of \
# case1 = epp(case1, **{'subDir_P11':'C:/Users/Ramon/Documents/ayudantia/code/28-08-21/eg_advance/020 aeroelástica vigas rígido/frecuencias y modos/', 'subDir_RSN':'C:/Users/Ramon/Documents/ayudantia/code/28-08-21/eg_advance/020 aeroelástica vigas rígido/frecuencias y modos/', 'subDir_FCS':'C:/Users/Ramon/Documents/ayudantia/code/28-08-21/eg_advance/020 aeroelástica vigas rígido/frecuencias y modos/', 'glob_print_output': True, 'BN_mode':'preserve'})
case1 = epp(case1, **{'subDir_P11':'C:/Users/Ramon/Documents/ayudantia/code/28-08-21/eg/subF/', 'subDir_RSN':'C:/Users/Ramon/Documents/ayudantia/code/28-08-21/eg/subF/', 'subDir_FCS':'C:/Users/Ramon/Documents/ayudantia/code/28-08-21/eg/subF/', 'BN_mode':'preserve'})


## Plot

#plot opts
save_opts = {'folder':'figures', 'filecode':'170', 'fig_format':None}
giro_1 = {str(case1.stru.nodes[0]):[6]}
giro_2 = {str(case1.stru.nodes[0]):[3,4,5]}
fig_geom_phit = plotter.fig_ut(case1.stru, [giro_2,giro_2], deg = True, fig_title = r'MPC - $/phi$', x_label = 't', y_label = r'$/phi$ [rad]',fig_save = True, fig_save_opts = save_opts, fig_name='geom_phit')
plotter.fig_ut_vt_pp(case1.stru,giro_1, deg = True)
plotter.fig_us(case1.stru,{'6':[100,150],'1':[100,120]},deg=True)
plotter.fig_uxuy(case1.stru,{'DOFs':[1,6],'t_vals':[100,150]},deg=True)
plotter.fig_u_spect(case1.stru,[giro_1,giro_1], y_units = 'rad/s')
plotter.fig_q_spect(case1.stru,[1,2],y_units='rad/s')
plotter.fig_u_FFT(case1.stru,giro_1,x_units='rad/s')
plotter.fig_q_FFT(case1.stru,[1],x_units='rad/s')

# tambien es posible guardar manualmente:
# plotter.save_figure(fig_phi_spect,save_opts,close=True)
