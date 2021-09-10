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
        case = sim_db.check_BN_files(case, **kwargs) #NOTA: Ver si usar un nombre loc o cambiarlo en case

        #Check if available:
        ##NOTA: Ahora desde una función, para ahorrar líneas
        sim_db.check_case_attris(case) #Esto debería imprimir alertas si falta algo

        case.stru = sim_db.rd_rawRespData(case.stru, **kwargs)
        sim_db.svBin(case, **kwargs)
        if glob_print_output:
            print("Raw response data read")
            
    #elif from bin file
    elif case.stru.struRdOpt == 'bin':
        if case.fName == '':
            print('Warning: fName empty!') ##NOTA: Esto podría hacerse con un 'kit' de opciones a verificar desde check_attris
        temp_case = sim_db.rdBin(case.fName, **kwargs) # expects "sim" class object
        case = sim_db.update_BN_objs(case, temp_case, **kwargs) #Se compara y actualiza para no perder info
        if type(case) != sim_db.sim: ##NOTA: Idem
            print('Warning: Not a sim class obj')
        if glob_print_output:
            print("Response data read from binary preprocessed file")
    
    # read external load data
    if case.stru.loadRdOpt == 'raw':
        case = sim_db.check_BN_files(case, **kwargs) #NOTA: Ver si usar un nombre loc o cambiarlo en case
        sim_db.check_case_attris(case) #Esto debería imprimir alertas si falta algo
        case.stru = sim_db.ae_Ftable(case.stru, **kwargs)
        sim_db.svBin(case, **kwargs)
        if glob_print_output:
            print("Raw external load data used")
    elif case.stru.loadRdOpt == 'bin':
        if case.fName == '':
            print('Warning, empty Filename')
        temp_case = sim_db.rdBin(case.fName,**kwargs) # expects "sim" class object
        case = sim_db.update_BN_objs(case,temp_case,**kwargs) #Nuevamente
        if type(case) != sim_db.sim:
            print('Warning: Not a sim class obj')
        if glob_print_output:
            print("External load data read from binary preprocessed file")
    elif case.stru.loadRdOpt == 'non':
        if glob_print_output:
            print("no external load data read")
    return case


def modalDecomp(struCase,**kwargs):
    """
    Applies modal decomposition
    input:
        case: 'sim' class obj
    kwargs: (may contain)
    returns:
        case, 'sim' class obj
    """
    
    if 'glob_print_output' in kwargs:
        glob_print_output = kwargs.get('glob_print_output')
    else:
        glob_print_output = False
    
    
    if len(struCase.auxMD) == 0:
        if (len(struCase.mass)!=0) and (struCase.phiR.shape[0]!=0):
            struCase.auxMD = np.zeros(struCase.phiR.T.shape)
            for i in range(struCase.auxMD.shape[0]):
               struCase.auxMD[i] = np.multiply(struCase.mass, struCase.phiR[:,i])
        elif glob_print_output:
            print("no mass and modal data for modal decomposition")
    
    if struCase.struEigOpt:
        struCase.q = np.matmul(struCase.auxMD, struCase.u_mdr)
    
    if struCase.loadEigOpt:
        if (struCase.eLoad.shape[0]!=0):
            struCase.Q = np.matmul(struCase.auxMD, struCase.eLoad)
        
        else:
            if glob_print_output:
                print("no external load data for modal decomposition")
    
    struCase = sim_db.modal_updnorm(struCase, 'norm_modalphiR', **kwargs)        
    return struCase

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
        case.stru = modalDecomp(case.stru,**kwargs)
    
    if case.stru.EigWorkOpt:
        case.stru = sim_db.modal_w(case.stru, **kwargs)
        #NOTA: Acá podría ir otro svBIN
        
    case = sim_db.check_BN_files(case, **kwargs) #NOTA: ¿Por qué esto acá? Si leo de binario no sirve re-guardar, o sí?
    sim_db.svBin(case,**kwargs) #Se exporta todo, finalmente.
    
    return(case)

"""
------------------------------------------------------------------------------
runing things
------------------------------------------------------------------------------
"""
if __name__ == '__main__':
    
    pass

