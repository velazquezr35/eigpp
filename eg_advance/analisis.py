"""
Created on 2021 08 12
@author: Mauro Maza
Departamento de Estructuras
Facultad de Ciencias Exactas, Físicas y Naturales
Universidad Nacioal de Córdoba
Córdoba, Argentina

analisis
    analize data from a set of simulations of one aeroelastic model at many 
    free stream velocity values
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
from sim_db import sim, rd_u, rd_rsn_De, svBin, rdBin
from eigpp  import modalDecomp
import plotter

"""
------------------------------------------------------------------------------
functions
------------------------------------------------------------------------------
"""
# def basic_general_plots(case1, **kwargs):
    # '''
    # General plots...
    
    # inputs:
        # case1, sim class obj
    
    # '''
    # #Valores cualquiera
    # desired = {'200000':[2]}
    # modes = [0]
    # modes_2 = [1]
    

    # plotter.fig_ut(case1.stru, [desired,desired], fig_title='Simple u(t) plot',x_label ='t',y_label='u')
    # plotter.fig_u_FFT(case1.stru, desired, fig_title='Simple FFT(u) plot', limit_tinds=[0,-1])
    # plotter.fig_u_spect(case1.stru,desired, fig_title = 'Simple spectrogram', x_label = 't', y_label='f')
    # plotter.fig_qt(case1.stru, [modes, modes_2],fig_title='Simple q(t) plot')
    
    
    # plotter.fig_us(case1.stru, {'2':[100]})
    # plotter.fig_qs(case1.stru,{'1':[150]})
    # plotter.fig_uxuy(case1.stru,{'DOFs':[2,3],'t_vals':[150,200]}, fig_title = 'Simple versus plot')
    # plotter.fig_q_FFT(case1.stru,modes,fig_title = 'Simple FFT(q) plot')
    # plotter.fig_qt_vt_pp(case1.stru,modes,fig_title='Simple comb. qt, vqt, phase p.')
    # plotter.fig_ut_vt_pp(case1.stru,desired, fig_title='Simple comb. ut,vt, phase p.')
    # plotter.fig_q_spect(case1.stru,modes, fig_title='Simple q spectrogram')
    # return('Done!')

"""
------------------------------------------------------------------------------
runing things
------------------------------------------------------------------------------
"""


# ANALIZE ONE CASE (ONE SIMULATION)
# case: v_inf = 170 ft/s
fName='R1700'
if False: # creat object, read raw data and save bin
    # 1. create object with basic data
    R1700=sim() # create a simulation object
    R1700.name='170.0'
    R1700.descr='v_\infty = 170.0 \, ft/s'
    R1700.fName=fName
    R1700.stru.nodes=[200000, 200001]
    # 2. read raw data
    #   2.1 DoFs data
    R1700.stru.p11FN='pcolgante.@1' # binary *.p11 file name (without extension - Simpact output) from wich extract generalized displacements
    R1700.stru=rd_u(R1700.stru, **{'subDir_P11': '020 aeroelástica vigas rígido/020 10 08x40/R100840_1700/'})
    #   2.2 mass and eigen modes
    R1700.stru.rsnDe='pcolgante' # ASCII *.rsn file name (without extension) - Delta output
    R1700.stru=rd_rsn_De(R1700.stru, **{'subDir_RSN': '020 aeroelástica vigas rígido/frecuencias y modos/'})
    # 3. modal decomposition
    R1700 = modalDecomp(R1700, **{'subDir_BIN': '020 aeroelástica vigas rígido/020 10 08x40/R100840_1700/'})
    # 4. save bin
    svBin(R1700, **{'subDir_BIN': '020 aeroelástica vigas rígido/020 10 08x40/R100840_1700/'})
else:
    # 1. read bin file
    R1700 = rdBin(fName, **{'subDir_BIN': '020 aeroelástica vigas rígido/020 10 08x40/R100840_1700/'})


# print some figures
# f1=plotter.fig_ut(R1700.stru, [{200000: [3]}, {200000: [4]}], **{'deg':[4]})
# f2=plotter.fig_u_FFT(R1700.stru, [{200000: [3]}, {200000: [4]}])
# f3=plotter.fig_u_spect(R1700.stru,[{200000: [3]}, {200000: [4]}])
# f4=plotter.fig_qt(R1700.stru, [1,2])
