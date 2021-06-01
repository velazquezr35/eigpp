# -*- coding: utf-8 -*-
"""
Created on Thu Feb 25 17:52:06 2021

@author:    Rodrigo Velazquez
            Mauro Maza

Ayudantía 2021 - Departamento de Estructuras FCEFyN
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
from sim_db import rd_mass, rd_eig, rd_u, save_bin, read_bin


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
    # NOTA: hay que hacer lo mismo que se hizo con el resto de las cosas
    pass


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
    if case.stru.struRdOpt == 'raw':
        # NOTA: acá hay que corregir: debe revisar primero si ya existe un archivo con el nombre que va a escribir
        #       y después tomar decisiones:
        #           - pedir permiso para sobreescribir,
        #           - guardar el nvo con otro nombre,
        #           - acutalizar solo la info que va a leer, <<-- creo que esta es la mejor, combinada con la anterior (gaurdar un nvo. archivo con un agregado, onda nombreViejo_fechaYHora)
        #           - una combinación de ellas,
        #       lo que sea (tampoco matarse programando)
        #       además, debería avisar si no tiene algunos datos que son necesarios, a saber:
        #           - case.stru.nodes
        #           - case.stru.p11FN
        #           - case.stru.rsnDe
        #       (creo que no falta ninguno)
        case.stru = rd_rawRespData(case.stru, **kwargs)
        save_bin( data_folder+case.fName, case, glob_print_output)
        if glob_print_output:
            print("Raw response data read")
    elif case.stru.struRdOpt == 'bin':
        # NOTA: esto debería revisar que exista el campo case.fName e imprmir un warning si no
        case = read_bin( data_folder+case.fName, glob_print_output) # expects "sim" class object
        # NOTA: esto debería imprimir un warning si "case" no es "sim class"
        if glob_print_output:
            print("Response data read from binary preprocessed file")
    
    # read external load data
    if case.stru.loadRdOpt == 'raw':
        # NOTA: acá hay que corregir lo mismo que arriba
        #       los datos necesarios en este caso son:
        #           - ?
        #       (creo que no falta ninguno)
        case.stru = rd_rawLoadData(case.stru, **kwargs)
        save_bin( data_folder+case.fName, case, glob_print_output)
        if glob_print_output:
            print("Raw external load data read")
    elif case.stru.loadRdOpt == 'bin':
        # NOTA: esto debería revisar que exista el campo case.fName e imprmir un warning si no
        case = read_bin( data_folder+case.fName, glob_print_output) # expects "sim" class object
        # NOTA: esto debería imprimir un warning si "case" no es "sim class"
        #       y otro si no se encuentra el archivo
        #       en teoría siempre va a encontrar el archivo porque lo lee de lo que se guardó antes (para la "respuesta")
        if glob_print_output:
            print("External load data read from binary preprocessed file")
    elif case.stru.loadRdOpt == 'non':
        if glob_print_output:
            print("no external load data read")
    
    return case


def modalDecomp(case):
    """
    applies modal decomposition
    
    case: "sim" class object
    """
    
    # preliminary 
    if case.stru.auxMD.size == 0:
        case.stru.auxMD = np.zeros((len(case.stru.stru.phi[0,:]),len(case.stru.stru.mass)))
        for i in range(len(case.stru.auxMD)):
            case.stru.auxMD[i] = np.multiply( case.stru.stru.mass, case.stru.stru.phi[:,i] )
    
    if case.stru.struEigOpt:
        case.stru.stru.q = np.matmul( case.stru.auxMD, case.stru.stru.u_avr )
    
    if case.stru.loadEigOpt:
        case.stru.stru.Q = np.matmul( case.stru.auxMD, case.stru.stru.eLoad )
    
    # NOTA: acá hay que guardar el binario más completo tmb, con los mismos cuidados que mencioné en la parte de lectura de datos
    
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
    # NOTA: no lo tengo claro, pero creo que data_folder y bin_folder deberían ser lo mismo
    
    # read data
    case = rd_data(case, **kwargs)
    
    # apply modal decomposition
    if case.stru.struEigOpt or case.stru.loadEigOpt:
        case = modalDecomp(case)
    
    # NOTA: acá va todo lo que sigue, si es que ponemos algo más, como graficar cosas o imprimir un informe de algún tipo
    # si es qeu lo ponemos acá o si hacemos otra "end-user funct" para hacer un reporte de algo
    # hay que pensarlo una vez que tengamos algo de salida (gráficos o cosas por el estilo) andando, porqeu recién ahí nos vamos a dar cuenta qué es lo mejor


"""
------------------------------------------------------------------------------
runing things
------------------------------------------------------------------------------
"""
if __name__ == '__main__':
    
    pass
