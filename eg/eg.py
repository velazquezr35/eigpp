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
from sim_db import sim, time_slice, upd_phiR
from eigpp  import epp
import plotter

"""
------------------------------------------------------------------------------
functions
------------------------------------------------------------------------------
"""

def my_fun():
    pass

"""
------------------------------------------------------------------------------
runing things
------------------------------------------------------------------------------
"""

# create simulation data object and populate minimum structural info
case1=sim() # this will contain data corresponding to a particular simulation (a "study case") under analysis
case1.fName='eg' # output binary file name (without extension)
case1.stru.name='str ex raw data'
case1.stru.nodes=[200001, 200002] # nodes of interest - only information relative to these nodes is going to be read
case1.stru.struRdOpt = 'bin' # set flag for reading ASCII file - default 'bin'
case1.stru.p11FN='pcolgante.@1' # binary *.p11 file name (without extension - Simpact output) from wich extract generalized displacements
case1.stru.rsnDe='pcolgante' # ASCII *.rsn file name (without extension) - Delta output
case1.stru.loadsFN = 'AeroFcsOnStruc' #ASCII *.dat file name - Loads
case1.stru.loadRdOpt = 'bin'
case1.stru.struEigOpt = True
case1.stru.loadEigOpt = True
case1.stru.intLabOffset = 24
#data dir: use / instead of \
# case1 = epp(case1, **{'subDir_P11':'C:/Users/Ramon/Documents/ayudantia/code/28-08-21/eg_advance/020 aeroelástica vigas rígido/frecuencias y modos/', 'subDir_RSN':'C:/Users/Ramon/Documents/ayudantia/code/28-08-21/eg_advance/020 aeroelástica vigas rígido/frecuencias y modos/', 'subDir_FCS':'C:/Users/Ramon/Documents/ayudantia/code/28-08-21/eg_advance/020 aeroelástica vigas rígido/frecuencias y modos/', 'glob_print_output': True, 'BN_mode':'preserve'})
case1 = epp(case1, **{'subDir_P11':'C:/Users/Ramon/Documents/ayudantia/code/28-08-21/eg/subF/', 'subDir_RSN':'C:/Users/Ramon/Documents/ayudantia/code/28-08-21/eg/subF/', 'subDir_FCS':'C:/Users/Ramon/Documents/ayudantia/code/28-08-21/eg/subF/', 'BN_mode':'preserve'})

case1.stru = upd_phiR(case1.stru, **{'MOI':[1],'subDir_P11':'C:/Users/Ramon/Documents/ayudantia/code/28-08-21/eg/subF/', 'subDir_RSN':'C:/Users/Ramon/Documents/ayudantia/code/28-08-21/eg/subF/', 'subDir_FCS':'C:/Users/Ramon/Documents/ayudantia/code/28-08-21/eg/subF/', 'BN_mode':'preserve'})

plotter.fig_qt(case1.stru, [[1,2]])
# case1.stru = time_slice(case1.stru)
# case1.stru = modal_w(case1.stru)

## Plot
#plot opts
# save_opts = {'folder':'figures', 'filecode':'170', 'fig_format':None}
# giro_1 = {str(case1.stru.nodes[0]):[6]}
# giro_2 = {str(case1.stru.nodes[0]):[3,4,5]}
# fig_geom_phit = plotter.fig_ut(case1.stru, [giro_2,giro_2], deg = True, fig_title = r'MPC - $/phi$', x_label = 't', y_label = r'$/phi$ [rad]',fig_save = True, fig_save_opts = save_opts, fig_name='geom_phit')
# plotter.fig_ut_vt_pp(case1.stru,giro_1, deg = True)
# plotter.fig_us(case1.stru,{'6':[100,150],'1':[100,120]},deg=True)
# plotter.fig_uxuy(case1.stru,{'DOFs':[1,6],'t_vals':[100,150]},deg=True)
# plotter.fig_u_spect(case1.stru,[giro_1,giro_1], y_units = 'rad/s')
# plotter.fig_q_spect(case1.stru,[1,2],y_units='rad/s')
# plotter.fig_u_FFT(case1.stru,giro_1,x_units='rad/s')
# plotter.fig_q_FFT(case1.stru,[1],x_units='rad/s')

