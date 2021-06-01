# -*- coding: utf-8 -*-
"""
Created on 2021 05 29
@author:    Rodrigo Velazquez
            Mauro Maza
Departamento de Estructuras
Facultad de Ciencias Exactas, Físicas y Naturales
Universidad Nacioal de Córdoba
Córdoba, Argentina

Defines classes and some related functions to store and manipulate
aeroelastic simulation data
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
import pickle
import os
import re


"""
------------------------------------------------------------------------------
classes
------------------------------------------------------------------------------
"""
class stru:
    """
    structural model info
    """
    
    # instance attributes
    def __init__(self):
        self.name   = ''                            # short description
        self.descr  = ''                            # description
                    
        self.nnode  = 0                             # 1         - number of nodes considered
        self.nodes  = []                            # nnode     - label of nodes considered (external)
        self.iLabl  = np.array([], dtype=int)       # nnode     - label of nodes considered (internal)
        self.ndof   = 0                             # 1         - number of DoF considered (including those where BCs are applied)
        self.nt     = 0                             # 1         - number of time steps
        self.t      = np.array([], dtype=float)     # nt        - temporal grid
        self.u_raw  = np.array([], dtype=float)     # ndof x nt - generalized displacements as functions of time - as read from "curvas"
                                                    #           - rotational DoFs (if exist) are expresed as 3-1-3 Euler angles rotations [rad]
        self.u_avr  = np.array([], dtype=float)     # ndof x nt - generalized displacements as functions of time - avr: axial vector rotations
                                                    #           - rotational DoFs (if exist) are expresed as axial vector rotations
        self.eLoad  = np.array([], dtype=float)     # ndof x nt - external loads as functions of time
                                                    #           - NOTA: qué onda con los momentos leídos acá? necesitan un tratamiento especial?
                    
        self.mass   = np.array([], dtype=float)     # ndof      - lumped mass matrix
        self.nm     = 0                             # 1         - number of modes read
        self.om     = np.array([], dtype=float)     # nm        - ordered natural frequencies
        self.phi    = np.array([], dtype=float)     # ndof x nm - modal matrix
        self.auxMD  = np.array([], dtype=float)     # nm x ndof - auxiliary matix PHI^T * M, used for modal decomposition
        self.q      = np.array([], dtype=float)     # nm x nt   - modal coordinates as functions of time
        self.Q      = np.array([], dtype=float)     # nm x nt   - modal external loads as functions of time
                    
        self.p11FN  = ''                            # binary *.p11 file name (without extension - Simpact output) from wich extract generalized displacements (and/or other data)
        self.rsnSi  = ''                            # ASCII *.rsn file name (without extension) - Simpact output
        self.rsnDe  = ''                            # ASCII *.rsn file name (without extension) - Delta output
        
        self.eqInfo = np.array([], dtype=float)     # Information Relative to the Equations Numbers
        
        self.rdof       = True   # True if rotational DoFs exist
        self.struRdOpt  = 'bin'  # reading data flag for structural response: 
                                 #   'raw': from ASCII data files
                                 #   'bin': from binary file with preprocessed data
        self.loadRdOpt  = 'non'  # reading data flag for external loading: 
                                 #   'raw': from ASCII data files
                                 #   'bin': from binary file with preprocessed data
                                 #   'non': no external load data available
        self.struEigOpt = False  # True if modal decomposition should be done over generalized displacements
        self.loadEigOpt = False  # True if modal decomposition should be done over external loads


class aero:
    """
    aerodynamic model info
    """
    
    pass
    
    # instance attributes
    # def __init__(self):
        # self.name  = ''                             # short description
        # self.descr = ''                             # description


class sim:
    """
    full model info
    container for simulation subparts
    """
    
    # instance attributes
    def __init__(self):
        self.name  = ''                             # short description
        self.descr = ''                             # description
        
        self.fName = ''                             # binary file name (without extension)
        
        self.stru       = stru() # 1 "stru" class objects - active
        self.aero       = aero() # 1 "aero" class objects - active
        
        self.struList   = []    # list of "stru" class objects
        self.aeroList   = []    # list of "aero" class objects


"""
------------------------------------------------------------------------------
functions
------------------------------------------------------------------------------
"""

# read Simpact and Delta output files ----------------------------------------

def callbat(file_data, nodo, dof, name_otp):
    """
    extracts generalized displacements data from *.p11 bin file
    """
    cmd = "(echo " + file_data + " 1 0 && echo d && echo " + nodo + " " + dof + " " + name_otp + " && echo s) | curvas"  
    os.system(cmd)


def line_spliter(line):
    """
    splits a line... how??? why???
    """
    ret = []
    local_start = 0
    for i in range(0,len(line)):
        if line[i] == "-" or line[i] == " " or line[i] == "\n":
            if not i == 0:
                if not line[i-1] == "E":
                    a = line[local_start:i]
                    ret.append(a)
                    local_start = i
    filt = []
    for i in range(len(ret)):
        try:
            filt.append(float(ret[i]))
        except:
            pass
    filt = np.array(filt)
    return(filt)


def search_string(x_dat, frase):
    # NOTA: agregar descripción y traducir comentarios
    lin_count = 0  # contador de lineas
    locs = []  # ubicaciones donde encuentra (line number python = -1 realidad)
    for txt in x_dat:
        lin_count = lin_count + 1
        x = re.search(frase, txt)
        if x != None:
            locs.append(lin_count)
            #print(lin_count)
    # devolvemos la lista con las lineas de x_dat que tienen la coincidencia
    
    return(locs)


def rd_SimpactTable(x_dat, start_line, **kwargs):
    """
    locates and reads a table from *.rsn ASCII Simpact and Alpha files
    """
    
    if 'glob_print_output' in kwargs:
        glob_print_output = kwargs.get('glob_print_output')
    else:
        glob_print_output = False
        
    stop_flag = True
    counter = 0
    while stop_flag:
        loc_arr = line_spliter(x_dat[start_line+counter])
        b_cond = np.isnan(loc_arr)
        # print(b_cond)
        if b_cond.any():
            stop_flag = False
            if glob_print_output:
                print("Nan encontrado")
                # NOTA: en general, todos los resultados desde que aparece un NaN en adelante no sirver,
                #       así que el tratamiento de esos casos debería ser borrar todos los datos que existan con t>=t_NaN
                #       habría que avisar tmb que se encontró un NaN y en qué tiempo y en qué variable y decir que todo se va a borrar desde ahí en adelante
            # break
        else:
            # print(loc_arr)
            if counter == 0:  # primer ciclo
                table_gen = np.array([loc_arr])
                # print(loc_arr)
            else:
                # print(loc_arr)
                table_gen = np.append(table_gen, [loc_arr], axis=0)
                
            counter = counter+1
            # print(len(x_dat))
            # print(loc_arr)
            try: 
                if x_dat[start_line + counter] == '\n':
                    stop_flag = False
            except:
                stop_flag=False
    
    return(table_gen)


def euler2axial(cols):
    """
    converts rotations expresed as 3-1-3 Euler Angles
    to axial vector form using SciPy class Rotation
    """
    
    from scipy.spatial.transform import Rotation as Rot
    
    for i in range(len(cols[:,0])):
        R = Rot.from_euler('ZXZ',cols[i,:],degrees=False)
        cols[i,:] = R.as_rotvec()
    return(cols)


def rd_u(struCase, **kwargs):
    """
    extracts generalized displacements data from *.p11 bin file
    and imports data to "stru" class object
    also creates "u_avr" field if necessary
    
    struCase: "stru" class object
    """
    if 'data_folder' in kwargs:
        data_folder = kwargs.get('data_folder')
    else:
        data_folder=''
    
    glob_u_raw = []
    glob_u_avr = []
    for node_name in struCase.nodes:
        callbat( data_folder+struCase.p11FN, str(node_name), "0", data_folder+"temp_file")
        loc_lines = open(data_folder+"temp_file",'r')
        loc_x_dat = loc_lines.readlines()
        loc_table_raw = rd_SimpactTable(loc_x_dat,0)
        loc_table_avr = np.copy(loc_table_raw)
        if struCase.rdof: # if rotational DoFs, create field u_avr
            loc_table_avr[:,4:]= euler2axial(loc_table_avr[:,4:])
        if struCase.nodes.index(node_name) == 0:
            #Acomodar posibilidad de que el tiempo de 1 nodo sea menor (NaN antes q resto)???
            glob_time = loc_table_raw[:,0]
            total_ntime = len(glob_time)
            glob_u_raw = np.transpose(loc_table_raw)[1:,:]
            glob_u_avr = np.transpose(loc_table_avr)[1:,:]
        else:
            glob_u_raw = np.append(glob_u_raw,np.transpose(loc_table_raw)[1:,:],axis=0)
            glob_u_avr = np.append(glob_u_avr,np.transpose(loc_table_avr)[1:,:],axis=0)
    
    loc_lines.close()
    os.remove(data_folder+"temp_file")
    
    struCase.nt = total_ntime
    struCase.t  = glob_time
    struCase.u_raw = glob_u_raw
    struCase.u_avr = glob_u_avr
    
    return struCase


def rd_eqInfo(struCase, **kwargs):
    """
    extracts Information Relative to the Equations Numbers from ASCII *.rsn file
    (Simpact or Delta output, default Delta)
    
    struCase: "stru" class object
    """
    
    if 'data_folder' in kwargs:
        data_folder = kwargs.get('data_folder')
    else:
        data_folder=''
        
    if 'c_info_eq' in kwargs:
        c_info_eq = kwargs.get('c_info_eq')
    else:
        c_info_eq=1
    
    if not struCase.rsnDe:
        struCase.rsnDe = struCase.rsnSi
    
    
    # open Delta *.rsn output file
    y = open(data_folder+struCase.rsnDe+'.rsn', 'r')
    x_dat = y.readlines()

    # find reference line and read table
    locs = search_string(x_dat, "Information Relative to the Equations Numbers")
    struCase.eqInfo = rd_SimpactTable(x_dat, locs[0]+c_info_eq)
    
    # close file
    y.close()
    
    return struCase


def rd_mass(struCase, **kwargs):
    """
    extracts lumped mass matrix from ASCII *.rsn file
    (Simpact or Delta output, default Delta)
    
    struCase: "stru" class object
    """
    
    if 'data_folder' in kwargs:
        data_folder = kwargs.get('data_folder')
    else:
        data_folder=''
    
    if 'glob_print_output' in kwargs:
        glob_print_output = kwargs.get('glob_print_output')
    else:
        glob_print_output = False
        
    if 'c_lumped_m' in kwargs:
        c_lumped_m = kwargs.get('c_lumped_m')
    else:
        c_lumped_m=1
    
    if not struCase.rsnDe:
        struCase.rsnDe = struCase.rsnSi
    
    struCase.nnode = len(struCase.nodes)
    
    if struCase.eqInfo.size == 0:
        struCase = rd_eqInfo(struCase, **kwargs)
    
    # open Delta *.rsn output file
    y = open(data_folder+struCase.rsnDe+'.rsn', 'r')
    x_dat = y.readlines()

    # find reference line and read table
    locs_lumped = search_string(x_dat, "  Lumped Mass Matrix")
    raw_lumped_matrix = rd_SimpactTable(x_dat, locs_lumped[0]+c_lumped_m)
    struCase.mass = []
    struCase.iLabl = []
    for j in range(0, len(struCase.eqInfo[:, 0])): #desplazamos fila a fila
        for a in range(0, struCase.nnode): #desplazamos nodo a nodo (de interés)
            if struCase.eqInfo[j, 0] == struCase.nodes[a]:
                # save data
                struCase.iLabl.append([struCase.nodes[a], struCase.eqInfo[j,-1],struCase.eqInfo[j,-2]]) # NOTA: hay que modificar esto para que lea solo la última columna, así queda funcionando también apra leer los *.rsn de Simpact
                struCase.mass = np.append(struCase.mass, raw_lumped_matrix[j])
                if glob_print_output:
                    print("Lumped Mass Matrix")
                    print(struCase.mass)
    
    # close file
    y.close()
    
    struCase.iLabl = np.array(struCase.iLabl)
    
    return struCase


def rd_eig(struCase, **kwargs):
    """
    extracts eigen modes and eigen frequencies from ASCII *.rsn file
    (Delta output)
    
    struCase: "stru" class object
    """
    
    if 'data_folder' in kwargs:
        data_folder = kwargs.get('data_folder')
    else:
        data_folder=''
    
    if 'glob_print_output' in kwargs:
        glob_print_output = kwargs.get('glob_print_output')
    else:
        glob_print_output = False
        
    if 'c_modes' in kwargs:
        c_modes = kwargs.get('c_modes')
    else:
        c_modes=6
    
    struCase.nnode = len(struCase.nodes)
    
    
    # open Delta *.rsn output file
    y = open(data_folder+struCase.rsnDe+'.rsn', 'r')
    x_dat = y.readlines()

    # find reference line and read table
    locs_modes = search_string(x_dat, "linear dynamic eigen-mode analysis")
    
    
    struCase.nm = len(locs_modes)
    struCase.om = np.zeros((struCase.nm, 3))
    struCase.phi = np.zeros((struCase.nm, struCase.nnode*6))
    for i in range(struCase.nm):
        struCase.om[i] = line_spliter(x_dat[locs_modes[i]])
        raw_mode_table = rd_SimpactTable(x_dat, locs_modes[i]+c_modes)
        if glob_print_output:
            print("tabla modo full")
            print(raw_mode_table)
            print(" modos------------------------------")
        local_mode_row = []
        for j in range(0, len(raw_mode_table[:, 0])): # number of rows in modes' table
            for a in range(0, struCase.nnode):
                if raw_mode_table[j, 0] == struCase.nodes[a]:
                    local_mode_row = np.append(local_mode_row, raw_mode_table[j,1:])
                    if glob_print_output:
                        print("aporte local ---------------")
                        print(local_mode_row)
        struCase.phi[i] = local_mode_row
    
    # close file
    y.close()
    
    struCase.phi = np.transpose(struCase.phi)
    
    return struCase


"""
def ae_Ftable(fname, noderefs):
    # NOTA: esta sólo la copié, pero no la revisé ni actualicé - hay que hacerlo - por eso está todo comentado
    x = open(fname, 'r')
    x_dat = x.readlines()
    init_line = 0
    init_cond = True
    while init_cond:
        try:
            loc_header = x_dat[init_line].split()
            float(loc_header[0]) + float(loc_header[1])
            init_cond = False
        except:
            init_line = init_line+1
    
    n_nods = int(x_dat[init_line+1])
    if glob_print_output:
        print("N nodes",n_nods)
    loc_ittab = []
    loc_ftab = []
    
    tab_cond = True
    act_line = init_line
    while tab_cond:
        try:
            loc_ftab.append(rd_SimpactTable(x_dat, act_line+2))
            loc_ittab.append([float(x_dat[act_line].split()[0]),float(x_dat[act_line].split()[1])])
            act_line = act_line + n_nods + 4
        except:
            if glob_print_output:
                print("Fin de tabla fzas, line: ", act_line)
            tab_cond = False
    #Ahora se guarda solamente lo que interesa de los nodos de refs:
    #Agrego manualmente para test. Pruebas con F y RSN que no se cruzan.
    noderefs[0,1] = 27
    noderefs[1,1] = 35

    counter = 0
    for a in loc_ftab: #por cada instante tiempo
        loc_frow_filt = np.array([])
        for i in range(len(a)): #por cada fila de tabla fzas (nodos)
           for j in range(len(noderefs[:,0])): #por cada nodo de interés
               if a[i,0] == noderefs[j,1]:
                   loc_frow_filt = np.append(loc_frow_filt, a[i,1:])
        if counter == 0:
           loc_ftab_filt = np.array([loc_frow_filt])
        else:
            loc_ftab_filt = np.append(loc_ftab_filt, [loc_frow_filt],axis=0)
        counter = counter + 1
    loc_ftab_filt = np.transpose(loc_ftab_filt)
    loc_ittab = np.array(loc_ittab)
    #Devolvemos la info
    #step, instante, tabla fuerza
    return(loc_ittab[:,0], loc_ittab[:,1], loc_ftab_filt)
"""


# handle postprocessed data files --------------------------------------------

def save_bin(file, data, msg:bool=False):
    
    """
    file: str - file name without extension
    data: variable to be saved
    msg: bool - True if printing message
    """
    
    f = open(file+'.sim', 'wb')
    pickle.dump(data, f)
    f.close()
    if msg:
        print ('bin data file saved (save_bin funct)')


def read_bin(file, msg:bool=False):
    
    """
    file: str - file name without extension
    data: variable read
    msg: bool - True if printing message
    """
    
    f = open(file+'.sim','rb')
    data = pickle.load(f)
    f.close()
    if msg:
        print ('bin data file read (read_bin funct)')
    
    return(data)


"""
------------------------------------------------------------------------------
Unit tests
------------------------------------------------------------------------------
"""
# NOTA: hay que hacer más unit tests
def uTest1():
    """
    determine if str class is fine to work with "pickle"
    """
    
    # create str object and populate
    a=str()
    a.nnode=2
    a.nodes=[200001,200003]
    a.ndof=6
    
    
    # str.t as list 
    a.t=[0.0, 0.1, 0.3]
    print('type str.t should be "list" at this point')
    print( type(a.t) )
    
    file1Name='file1.bin'
    print('save')
    save_bin(file1Name, a)    
    print('read')
    b=read_bin(file1Name, True)
    print('type str.t should be "list" at this point')
    print( type(b.t) )
    os.remove(file1Name)
    
    
    # str.t as numpy array 
    a.t=np.array([0.0, 0.1, 0.3], dtype=float)
    print('type str.t should be "numpy.ndarray" at this point')
    print( type(a.t) )
    
    file1Name='file1.bin'
    print('save')
    save_bin(file1Name, a)    
    print('read')
    b=read_bin(file1Name, True)
    print('type str.t should be "numpy.ndarray" at this point')
    print( type(b.t) )
    os.remove(file1Name)
    

"""
------------------------------------------------------------------------------
runing things
------------------------------------------------------------------------------
"""
if __name__ == '__main__':
    
    pass
    