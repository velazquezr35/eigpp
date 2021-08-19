# -*- coding: utf-8 -*-
"""
Created on Thu Feb 25 17:52:06 2021

@author:    Rodrigo Velazquez
            Mauro Maza

Ayudantía 2021 - Departamento de EStructuras FCEFyN
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
from sim_db import svBin, rdBin, sim, check_BN_files, rd_rawRespData, ae_Ftable #NOTA RV: Clases con mayuscula
import sim_db


"""
------------------------------------------------------------------------------
classes
------------------------------------------------------------------------------
"""

class MyClass(): # just an example
    pass


"""
------------------------------------------------------------------------------
functions
------------------------------------------------------------------------------
"""
def rd_data(case, **kwargs):
    """
    reads data
    
    case: "sim" class object
    """
    
    if 'subDir_P11' in kwargs:
        subDir_P11 = kwargs.get('subDir_P11')
    else:
        subDir_P11=''
    
    if 'glob_print_output' in kwargs:
        glob_print_output = kwargs.get('glob_print_output')
    else:
        glob_print_output = False
    
    # read response data
    #   - generalized displacements over time
    #   - eigen modes and eigen frequencies
    #   - mass matrix
    
    #If raw
    if case.stru.struRdOpt == 'raw':
        case = check_BN_files(case, **kwargs) #NOTA: Ver si usar un nombre loc o cambiarlo en case

        #Check if available:
        ##NOTA: Ahora desde una función, para ahorrar líneas
        sim_db.check_case_attris(case) #Esto debería imprimir alertas si falta algo

        case.stru = rd_rawRespData(case.stru, **kwargs)
        svBin(case, **kwargs)
        if glob_print_output:
            print("Raw response data read")
            
    #elif from bin file
    elif case.stru.struRdOpt == 'bin':
        if case.fName == '':
            print('Warning: fName empty!') ##NOTA: Esto podría hacerse con un 'kit' de opciones a verificar desde check_attris
        temp_case = rdBin(case.fName, **kwargs) # expects "sim" class object
        case = sim_db.update_BN_objs(case, temp_case, **kwargs) #Se compara y actualiza para no perder info
        if type(case) != sim: ##NOTA: Idem
            print('Warning: Not a sim class obj')
        if glob_print_output:
            print("Response data read from binary preprocessed file")
    
    # read external load data
    if case.stru.loadRdOpt == 'raw':
        case = check_BN_files(case, **kwargs) #NOTA: Ver si usar un nombre loc o cambiarlo en case
        sim_db.check_case_attris(case) #Esto debería imprimir alertas si falta algo
        case.stru = ae_Ftable(case.stru, **kwargs)
        svBin(case, **kwargs)
        if glob_print_output:
            print("Raw external load data used")
    elif case.stru.loadRdOpt == 'bin':
        if case.fName == '':
            print('Warning, empty Filename')
        temp_case = rdBin(case.fName,**kwargs) # expects "sim" class object
        case = sim_db.update_BN_objs(case,temp_case,**kwargs) #Nuevamente
        if type(case) != sim:
            print('Warning: Not a sim class obj')
        if glob_print_output:
            print("External load data read from binary preprocessed file")
    elif case.stru.loadRdOpt == 'non':
        if glob_print_output:
            print("no external load data read")
    return case


def modalDecomp(case,**kwargs):
    """
    Applies modal decomposition

    
    case: "sim" class object
    """
    if 'subDir_P11' in kwargs:
        subDir_P11 = kwargs.get('subDir_P11')
    else:
        subDir_P11=''
    if len(case.stru.auxMD) == 0: ##NOTA: ¿Por qué si 0??
        case.stru.auxMD = np.zeros((len(case.stru.phi[0,:]),len(case.stru.mass)))
        #NOTA: Else, la recuperamos de la clase?
        for i in range(len(case.stru.auxMD)):
            case.stru.auxMD[i] = np.multiply(case.stru.mass, case.stru.phi[:,i])
    
    if case.stru.struEigOpt:
        case.stru.q = np.matmul( case.stru.auxMD, case.stru.u_mdr )
    
    if case.stru.loadEigOpt:
        case.stru.Q = np.matmul(case.stru.auxMD, case.stru.eLoad)
    
    case = check_BN_files(case, **kwargs)
    svBin(case,**kwargs) #Se exporta todo, finalmente.
    return case


"""
------------------------------------------------------------------------------
end-user functs
------------------------------------------------------------------------------
"""
def epp(case, **kwargs):
    """
    global eigen analysis
    
    case: "sim" class object
    kwargs:
        -subDir_RSN, subDir_P11, subDir_FCS:  str
                        path to data from current folder
                        e.g.: 'subFolder1/subFolder2/folderWhereDataIs/'
                        default: ''
        - glob_print_output:    bool
                                info printing flag
                                default: False
    """

    case = rd_data(case, **kwargs)
    if case.stru.struEigOpt or case.stru.loadEigOpt:
        case = modalDecomp(case,**kwargs)
        

    # NOTA: acá va todo lo que sigue, si es que ponemos algo más, como graficar cosas o imprimir un informe de algún tipo
    # si es que lo ponemos acá o si hacemos otra "end-user funct" para hacer un reporte de algo
    # hay que pensarlo una vez que tengamos algo de salida (gráficos o cosas por el estilo) andando, porqeu recién ahí nos vamos a dar cuenta qué es lo mejor
    return(case)


"""
------------------------------------------------------------------------------
runing things
------------------------------------------------------------------------------
"""
if __name__ == '__main__':
    
    pass

