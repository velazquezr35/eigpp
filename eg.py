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
runing things
------------------------------------------------------------------------------
"""

# create simulation data object and populate minimum structural info
case1=sim() # this will contain data corresponding to a particular simulation (a "study case") under analysis
case1.fName='eg' # output binary file name (without extension)
case1.stru.name='str ex raw data'
case1.stru.nodes=[200000, 200010] # nodes of interest - only information relative to these nodes is going to be read
case1.stru.struRdOpt='bin' # set flag for reading ASCII file - default 'bin'
case1.stru.p11FN='pcolgante.@1' # binary *.p11 file name (without extension - Simpact output) from wich extract generalized displacements
case1.stru.rsnDe='pcolgante' # ASCII *.rsn file name (without extension) - Delta output
case1.stru.loadsFN = 'AeroFcsOnStruc' #ASCII *.dat file name - Loads
case1.stru.loadRdOpt  = 'bin'

# eigen modes postprocess
#   call the function
#   asking for
#       work over data in a subfolder called "subF"
#       print status messages
case1 = epp(case1, **{'data_folder': 'subF/', 'glob_print_output': True, 'BN_mode':'preserve'})

###############################################################################################################

# Tests

###############################################################################################################

#Usar un dict resultó más limpio para el end-user, aunque internamente se hagan 2 pasos extras.
#Por ejemplo:
desired = {'200000':[1]}

#Tests
import numpy as np
case1.stru.t = np.linspace(0,100,5000)
case1.stru.u_avr = np.zeros((2,5000))
case1.stru.q = np.copy(case1.stru.u_avr)
case1.stru.u_avr[1] = np.sin(8*case1.stru.t)
case1.stru.q = np.copy(case1.stru.u_avr)

# # plotter.fig_ut_vt_pp(case1.stru,desired,fig_title='gg')
# # plotter.fig_u_spect(case1.stru, desired, fig_title = 'FFT')

# # plotter.fig_qt(case1.stru, [{'1':[0,1]},{'2':[0]}])
# # plotter.fig_qPP(case1.stru, [{'1':[0,1]},{'2':[0]}])
# plotter.fig_qt_vt_pp(case1.stru,{'a':[1]})