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
from sim_db import ae_Ftable, rd_mass, rd_eig, rd_u, save_bin, read_bin, sim, check_BN_files #NOTA RV: Clases con mayuscula
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

def rd_rawRespData(struCase, **kwargs):
    """
    reads response data from Simpact and/or Delta output
        - generalized displacements over time
        - eigen modes and eigen frequencies
        - mass matrix
    
    struCase: "stru" class object
    """
    
    struCase = rd_mass(struCase, **kwargs)
    struCase = rd_eig(struCase, **kwargs)
    struCase = rd_u(struCase, **kwargs)
    
    return struCase


def rd_rawLoadData(struCase, **kwargs):
    """
    reads external load data
    
    struCase: "stru" class object
    """
    struCase = ae_Ftable(struCase, **kwargs)
    
    return struCase


def rd_data(case, **kwargs):
    """
    reads data
    
    case: "sim" class object
    """
    
    if 'data_folder' in kwargs:
        data_folder = kwargs.get('data_folder')
    else:
        data_folder=''
    
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
        save_bin(data_folder+case.fName, case, glob_print_output)
        if glob_print_output:
            print("Raw response data read")
            
    #elif from bin file
    elif case.stru.struRdOpt == 'bin':
        if case.fName == '':
            print('Warning: fName empty!') ##NOTA: Esto podría hacerse con un 'kit' de opciones a verificar desde check_attris
        temp_case = read_bin(data_folder+case.fName, glob_print_output) # expects "sim" class object
        case = sim_db.update_BN_objs(case, temp_case, **kwargs) #Se compara y actualiza para no perder info
        if type(case) != sim: ##NOTA: Idem
            print('Warning: Not a sim class obj')
        if glob_print_output:
            print("Response data read from binary preprocessed file")
    
    # read external load data
    if case.stru.loadRdOpt == 'raw':
        case = check_BN_files(case, **kwargs) #NOTA: Ver si usar un nombre loc o cambiarlo en case
        sim_db.check_case_attris(case) #Esto debería imprimir alertas si falta algo
        case.stru = rd_rawLoadData(case.stru, **kwargs)
        save_bin(data_folder+case.fName, case, glob_print_output)
        if glob_print_output:
            print("Raw external load data used")
    elif case.stru.loadRdOpt == 'bin':
        if case.fName == '':
            print('Warning, empty Filename')
        temp_case = read_bin(data_folder+case.fName, glob_print_output) # expects "sim" class object
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
    if 'data_folder' in kwargs:
        data_folder = kwargs.get('data_folder')
    else:
        data_folder=''
    if len(case.stru.auxMD) == 0: ##NOTA: ¿Por qué si 0??
        case.stru.auxMD = np.zeros((len(case.stru.phi[0,:]),len(case.stru.mass)))
        #NOTA: Else, la recuperamos de la clase?
        for i in range(len(case.stru.auxMD)):
            case.stru.auxMD[i] = np.multiply(case.stru.mass, case.stru.phi[:,i])
    
    if case.stru.struEigOpt:
        case.stru.q = np.matmul( case.stru.auxMD, case.stru.u_avr )
    
    if case.stru.loadEigOpt:
        case.stru.Q = np.matmul(case.stru.auxMD, case.stru.eLoad)
    
    case = check_BN_files(case, **kwargs)
    save_bin(data_folder+case.fName, case) #Se exporta todo, finalmente.
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
        - data_folder:  str
                        path to data from current folder
                        e.g.: 'subFolder1/subFolder2/folderWhereDataIs/'
                        default: ''
        - glob_print_output:    bool
                                info printing flag
                                default: False
    """

    # read data
    case = rd_data(case, **kwargs)
    # apply modal decomposition
    #pero si ya la leí de antes...?
    # if case.stru.struRdOpt == 'raw' or case.stru.loadRdOpt == 'raw':
    if case.stru.struEigOpt or case.stru.loadEigOpt:
        case = modalDecomp(case,**kwargs)
        

    # NOTA: acá va todo lo que sigue, si es que ponemos algo más, como graficar cosas o imprimir un informe de algún tipo
    # si es qeu lo ponemos acá o si hacemos otra "end-user funct" para hacer un reporte de algo
    # hay que pensarlo una vez que tengamos algo de salida (gráficos o cosas por el estilo) andando, porqeu recién ahí nos vamos a dar cuenta qué es lo mejor
    return(case)


"""
------------------------------------------------------------------------------
runing things
------------------------------------------------------------------------------
"""
if __name__ == '__main__':
    
    pass

