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
from scipy.spatial.transform import Rotation as Rot

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
        self.name   = ''                                    # short description
        self.descr  = ''                                    # description
                    
        self.nnode  = 0                                     # 1          - number of nodes considered
        self.nodes  = []                                    # nnode      - label of nodes considered (external)
        self.iLabl  = np.array([], dtype=int)               # nnode      - label of nodes considered (internal)
        self.ndof   = 0                                     # 1          - Total number of DoF considered (including those where BCs are applied)
        self.nt     = 0                                     # 1          - number of time steps
        self.t      = np.array([], dtype=float)             # nt         - Response temporal grid
        self.u_raw  = np.array([], dtype=float)             # ndof x nt  - generalized displacements as functions of time - as read from "curvas"
                                                            #            - rotational DoFs (if exist) are expresed as 3-1-3 Euler angles rotations [rad]
        self.u_mdr  = np.array([], dtype=np.longdouble)     # ndof x nt  - generalized displacements as functions of time - mdr: modal decomposition ready
                                                            #            - rotational DoFs (if exist) are expresed as axial vector rotations relative to the initial orientation of each node's local system
        self.aLoad  = np.array([], dtype=float)             # ndof x nt  - aerodynamic loads as functions of time
        self.LW     = np.array([], dtype=np.longdouble)     # ndof x nt  - work from external loads (over GDoFs), as functions of time
        self.LWtot  = np.array([], dtype=np.longdouble)     # 1 x nt     - sum of LW over the GDoFs, as a function of time
        
        self.mass   = np.array([], dtype=float)             # ndof       - lumped mass matrix
        self.nm     = 0                                     # 1          - number of modes read
        self.om     = np.array([], dtype=float)             # nm         - ordered natural frequencies
        self.phi    = np.array([], dtype=np.longdouble)     # ndof x nm  - modal matrix - as read from *.RSN - must be normalized w.r.t. the mass
        self.moi    = []                                    # 1          - modes of interest indices
        self.nmoi   = 0                                     # 1          - number of modes of interest indices
        self.mnorm  = 'mass'                                # str        - moi normalization criteria
        self.omR    = np.array([], dtype=float)             # nmoi       - reduced ordered natural frequencies (using moi) - if moi=[] => omR=om
        self.phiR   = np.array([], dtype=np.longdouble)     # ndof x nmoi- reduced modal matix (using moi) - if moi=[] => phiR=phi
        self.mmass  = np.array([], dtype=np.longdouble)     # nmoi       - modal mass matrix (diagonal) - reduced (for moi only)
        self.mstif  = np.array([], dtype=np.longdouble)     # nmoi       - modal stiffness matrix (diagonal) - reduced (for moi only)
        self.auxMD  = np.array([], dtype=np.longdouble)     # nmoi x ndof- auxiliary matix PHI^T * M, used for modal decomposition
        self.q      = np.array([], dtype=np.longdouble)     # nmoi x nt  - modal coordinates as functions of time
        self.Q      = np.array([], dtype=np.longdouble)     # nmoi x nt  - modal external loads as functions of time
        self.mKE    = np.array([], dtype=np.longdouble)     # nmoi x nt  - modal Kinetic Energy as functions of time
        self.mPE    = np.array([], dtype=np.longdouble)     # nmoi x nt  - modal Potential Energy as functions of time
        self.mME    = np.array([], dtype=np.longdouble)     # nmoi x nt  - modal Mechanical Energy as functions of time
        self.mKEtot = np.array([], dtype=np.longdouble)     # 1 x nt     - sum of mKE over the MDoFs, as a function of time
        self.mPEtot = np.array([], dtype=np.longdouble)     # 1 x nt     - sum of mPE over the MDoFs, as a function of time
        self.mMEtot = np.array([], dtype=np.longdouble)     # 1 x nt     - sum of mME over the MDoFs, as a function of time
        self.QW     = np.array([], dtype=np.longdouble)     # nmoi x nt  - modal work from external loads, as function of time
        self.QWtot  = np.array([], dtype=np.longdouble)     # 1 x nt     - sum of QW over the MDoFs, as a function of time
        
        self.p11FN  = ''                            # binary *.p11 file name (without extension - Simpact output) from wich extract generalized displacements (and/or other data)
        self.rsnSi  = ''                            # ASCII *.rsn file name (without extension) - Simpact output
        self.rsnDe  = ''                            # ASCII *.rsn file name (without extension) - Delta output
        self.loadsFN = ''                           # ASCII *. ??? file name (without extension) - Loads on stru
        self.t_Nan = np.inf                         # inf       - NaN minimum time
        
        
        self.eqInfo = np.array([], dtype=float)             # Information Relative to the Equations Numbers
        
        self.rdof       = True                              # True if rotational DoFs exist
        self.struRdOpt  = 'raw'                             # reading data flag for structural response: 
                                                            #   'raw': from ASCII data files
                                                            #   'bin': from binary file with preprocessed data
        self.loadRdOpt  = 'raw'                             #reading data flag for external loading: 
                                                            #   'raw': from ASCII data files
                                                            #   'bin': from binary file with preprocessed data
                                                            #   'non': no external load data available
        self.struEigOpt = True                              # True if modal decomposition should be done over generalized displacements
        self.loadEigOpt = True                              # True if modal decomposition should be done over external loads
        self.EigWorkOpt = True                              # True if modal work from external loads should be computed
        self.plot_timeInds = np.array([0,None])               # desired plot indexes
        self.plot_timeVals = np.array([np.inf,np.inf])      # desired plot time values
        self.intLabOffset = 0                               # offset node labels
        self.rot_inds = [4,5,6]                             # rotational DOFs inds (not Python´s)
    
    #Methods
    #Coming soon...

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

def upd_phiR(struCase):
    '''
    Updates struCase.phiR and struCase.omR.
    
    Inputs:
        struCase, stru class obj
    Returns:
        struCase, updated
    '''

    if len(struCase.moi) == 0:
        struCase.phiR = np.copy(struCase.phi)
        struCase.omR = np.copy(struCase.om)
    else:
        moi_inds = []
        for i in range(len(struCase.moi)):
            moi_inds.append(struCase.moi[i]-1)
        struCase.phiR = struCase.phi[:,moi_inds]
        struCase.omR = struCase.om[moi_inds]
        
    return(struCase)

def upd_modalMK(struCase):
    '''
    Determines struCase.mmass and struCase.mstif from struCase.phiR and struCase.mass.
    Inputs:
        struCase, stru class obj
    Returns:
        struCase, updated
    '''
    
    struCase.mmass = np.diag( np.matmul( struCase.phiR.T, np.matmul( np.diag(struCase.mass), struCase.phiR ) ) )
    struCase.mstif = np.multiply( struCase.mmass, np.power( struCase.omR[:,1], 2 ) )
        
    return(struCase)

def upd_mnorm(struCase):
    '''
    Updates struCase.phiR, struCase.q, struCase.Q,
    struCase.mmass and struCase.mstif from struCase.mnorm.
    
    struCase.phiR, struCase.q, struCase.Q before this
    function is called must have been determined using
    mass-normalized modes.
    
    Inputs:
        struCase, stru class obj
    Returns:
        struCase, updated
    '''

    if struCase.mnorm == 'mass':
        alpha = 1

    elif struCase.mnorm == 'stiff':
        alpha = 1/struCase.om[:,2]
        
    elif struCase.mnorm == 'norm':
        alpha = 1/np.linalg.norm(struCase.phiR,axis=0)
    
    elif struCase.mnorm == 'max':
        alpha = 1/np.max(struCase.phiR, axis=0)
        
    struCase.phiR = struCase.phiR*alpha
    
    if struCase.struEigOpt:
        struCase.q = np.multiply(struCase.q,np.transpose([1/alpha]))
    
    if struCase.loadEigOpt:
        struCase.Q = np.multiply(struCase.Q,np.transpose([1/alpha]))
        
    struCase = upd_modalMK(struCase)
    
    return(struCase)


def modalDecomp(struCase,**kwargs):
    """
    Applies modal decomposition over stru´s DOFs and LOADS (if available).
    Modal decomposition can be by mass (def), stiff, max or norm.
    
    Inputs:
        struCase: 'stru' class obj
        **kwargs: (may contain)
            general kwargs for downstream funcs
    Returns:
        struCase, 'stru' class obj, updated
    """
    
    if 'glob_print_output' in kwargs:
        glob_print_output = kwargs.get('glob_print_output')
    else:
        glob_print_output = False
        
    struCase = upd_phiR(struCase)
        
    if (len(struCase.mass)!=0) and (struCase.phiR.shape[0]!=0):
        struCase.auxMD = np.zeros(struCase.phiR.T.shape)
        for i in range(struCase.auxMD.shape[0]):
           struCase.auxMD[i] = np.multiply(struCase.mass, struCase.phiR[:,i])
    elif glob_print_output:
        print("no mass and modal data for modal decomposition")
    
    if struCase.struEigOpt:
        struCase.q = np.matmul(struCase.auxMD, struCase.u_mdr)
    
    if struCase.loadEigOpt:
        if (struCase.aLoad.shape[0]!=0):
            struCase.Q = np.matmul(struCase.phiR.T, struCase.aLoad)
        else:
            if glob_print_output:
                print("no external load data for modal decomposition")
    
    if struCase.mnorm == 'mass':
        struCase = upd_modalMK(struCase)
    else:
        struCase = upd_mnorm(struCase)
        
    return struCase


def modal_mechEnergy(struCase):
    '''
    Determines mechanical energy related to each mode.
    Inputs:
        struCase: 'stru' class obj
    Returns:
        struCase, updated
    '''
    
    dq = np.gradient( struCase.q, struCase.t, axis=1 )
    struCase.mKE = 0.5*np.multiply( struCase.mmass, np.power(dq,2).T ).T
    struCase.mPE = 0.5*np.multiply( struCase.mstif, np.power(struCase.q,2).T ).T
    struCase.mME = struCase.mKE + struCase.mPE
    
    struCase.mKEtot = np.sum(struCase.mKE, axis=0 )
    struCase.mPEtot = np.sum(struCase.mPE, axis=0 )
    struCase.mMEtot = np.sum(struCase.mME, axis=0 )
        
    return struCase


def search_time(t_array, t_values, **kwargs):
    '''
    Searchs for indexes (equal or max) in a t_array corresponding to some t_values
    
    Inputs:
        t_array, numpy ndarray 1D
        t_values, list [start, end]
    
    Returns:
        [ind_start, ind_end], 2-len list
    '''
    t_abs = abs(t_array-t_values[0])
    ind_start = list(t_abs).index(min(t_abs))
    t_abs = abs(t_array-t_values[1])
    ind_end = list(t_abs).index(min(t_abs))
    return([ind_start,ind_end])

def sfti_time(struCase, *reset, **kwargs):
    '''
    Searchs for time indexes - Updates the desired indexes for plotting purposes
    Inputs:
        struCase stru class obj
    **kwargs:
        must contain:
                indexes, list, time indexes [start, end] for struCase.t
            or
                time_vals, list, time values [start, end] for struCase.t
    
    Returns:
        struCase, 'stru' class obj, updated
    '''
    if not reset:
        if 'indexes' in kwargs:
            time_inds = kwargs.get('indexes')
            struCase.plot_timeInds[:] = time_inds[:]
            return(struCase)
        elif 'time_vals' in kwargs:
            time_vals = kwargs.get('time_vals')
            time_inds = search_time(struCase.t,time_vals)
            struCase.plot_timeInds[:] = time_inds[:]
            struCase.plot_timeVals[0] = struCase.t[time_inds[0]]
            struCase.plot_timeVals[1] = struCase.t[time_inds[1]]
            return(struCase)
        else:
            print('Warning: Time interval not set')
    else:
        struCase.plot_timeInds = [0,None]
        struCase.plot_timeVals = [struCase.t[0],struCase.t[-1]]
    return(struCase)
    
    
def nodeDof2idx(struCase, nodeDOFs):
    """
    Returns indexes for the nodes and nodes DOFs of interest
    Inputs:
        struCase: sim.stru obj
        nodeDOFs: dict, keys: nodes, values: list of DOFs per node
    Returns:
        list of ints (indexes), sorted
    """
    loc_indexes = []
    nodes = list(nodeDOFs.keys())
    dofs = list(nodeDOFs.values())
    
    for i in range(len(nodes)):
        for j in range(len(dofs[i])):
            try:
                if dofs[i][j] > 6 or dofs[i][j] < 1:
                    print('DOF local mal definido')
                    raise ValueError()
                loc_indexes.append(6*struCase.nodes.index(int(nodes[i]))+(dofs[i][j]-1))
            except:
                print('Error - Revisar datos: nodo:', nodes[i], ', DOF: ', (dofs[i][j]-1))
                break
    return loc_indexes


def modalDof2idx(struCase, modalDOF):
    '''
    Returns modal indexes for a particular shape-DOF
    Inputs:
        struCase, stru class obj
        modalDOF, int - desired DOF
    Returns:
        loc_indexes, numpy ndarray
    '''
    loc_indexes = []
    if modalDOF > 6 or modalDOF < 1:
        raise ValueError('Warning: Bad modal DOF', modalDOF)
    for i in range(len(struCase.nodes)):
        loc_indexes.append(6*i+modalDOF-1)
    return np.array(loc_indexes)
    


def callbat(file_folder, file_name, nodo, dof, output_name):
    """
    Extracts generalized displacements data from *.p11 bin file. Sends windows cmd commands.
    
    Inputs:
        file_folder: str, file location
        file_name: str, file name
        nodo: str, desired node
        dof: str, desired dof
        output_name: str, final name
    Returns:
        
    """
    cmd = "cd " + file_folder + " && (echo " + file_name + " 1 0 && echo d && echo " + nodo + " " + dof + " " + output_name + " && echo s) | curvas"  
    os.system(cmd)


def line_spliter(line):
    """
    Splits a string line. For general use with SimpactTable. 
    Inputs:
        line: str, desired string
    Returns:
        filt: numpy ndarray, filtered indexes
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
        if ret[i][0] == '*':
            filt.append(-1)
        else:
            try:
                filt.append(float(ret[i]))
            except:
                pass
    filt = np.array(filt)
    return(filt)


def search_string(x_dat, frase):
    '''
    Searchs for a keyword or string in a list of strings
    Inputs:
        x_dat: list, string container
        frase: str, desired string
    Returns:
        locs:list, desired match
    '''
    lin_count = 0  #Line count
    locs = []  #Index. In Python, first = [0]
    for txt in x_dat:
        lin_count = lin_count + 1
        x = re.search(frase, txt)
        if x != None:
            locs.append(lin_count) #Save the local match
    return(locs)


def rd_SimpactTable(x_dat, start_line, **kwargs):
    """
    Locates and reads a table from *.rsn ASCII Simpact and Alpha files
    Inputs:
        x_dat: list, string container
        start_line: int, starting index
    **kwargs (may contain):
        'STable_mode', 'normal' or 'with_counter'
            'with_counter' - also returns the Nan counter
    Returns:
        table_gen: numpy ndarray, data read
    """
    
    if 'glob_print_output' in kwargs:
        glob_print_output = kwargs.get('glob_print_output')
    else:
        glob_print_output = False
        
    if 'STable_mode' in kwargs:
        S_mode = kwargs.get('STable_mode')
    else:
        S_mode = 'normal'
        
    stop_flag = True
    counter = 0
    loc_Nan_counter = 0
    while stop_flag:
        loc_arr = line_spliter(x_dat[start_line+counter])
        b_cond = np.isnan(loc_arr)
        if b_cond.any():
            stop_flag = False
            print('NaN found after ', counter, ' timesteps')
            loc_Nan_counter = counter
        else:
            if counter == 0:  #First cycle
                table_gen = np.array([loc_arr])
            else:
                table_gen = np.append(table_gen, [loc_arr], axis=0)
                
            counter = counter+1

            try: 
                if x_dat[start_line + counter] == '\n' or x_dat[start_line + counter][:4] == '  iw':
                    stop_flag = False
            except:
                stop_flag=False
    if glob_print_output:
        print(' ')
        print('Table:')
        print(table_gen)
        print(' ')
    if S_mode == 'normal':
        return table_gen
    elif S_mode == 'with_counter':
        return table_gen, loc_Nan_counter

def rd_rsn_De(struCase, **kwargs):
    """
    Reads data from Delta *.rsn output (eigen modes and eigen frequencies, mass matrix) and performs modal decomposition.
    
    Inputs:
        struCase, 'stru' class object
    **kwargs:
        General kwargs for downstream funcs
    Returns:
        struCase, 'stru' class object, updated
    """
    
    struCase = rd_mass(struCase, **kwargs)
    struCase = rd_eig(struCase, **kwargs)
    struCase = modalDecomp(struCase,**kwargs)
    
    return struCase

def euler2axial(cols):
    """
    Converts rotations expresed as 3-1-3 Euler Angles to axial vector form using SciPy class Rotation
    
    Inputs:
        cols: numpy ndarray, data
    Returns:
        cols: numpy ndarray, data (rot)
    """
    
    from scipy.spatial.transform import Rotation as Rot
    
    for i in range(len(cols[:,0])):
        R = Rot.from_euler('ZXZ',cols[i,:],degrees=False)
        cols[i,:] = R.as_rotvec()
    return(cols)

def rotModalDec(cols):
    """
    Prepares rotational data for modal decomposition. Represent rotations as incremental rotation vectors expressed in initial local system.
    
    Inputs:
        cols: numpy ndarray, data
    Returns:
        cols: numpy ndarray, data (rot)
    
    Performs:
        o_0 = M_0 * g
        o_i = M_i * g = M_r * M_0 * g
        M_i * M_0^T = M_r
    """

    R_0 = Rot.from_euler('ZXZ',cols[0,:],degrees=False) # rotation to initial orientation
    cols[0,:] = Rot.as_rotvec(Rot.from_matrix(np.diag([1,1,1])))
    for i in range(1,len(cols[:,0])):
        R_i = Rot.from_euler('ZXZ',cols[i,:],degrees=False) # global rotation to instantaneous orientation
        R_r = R_i*R_0.inv() # rotation relative to initial orientation
        cols[i,:] = R_0.inv().apply( R_r.as_rotvec() ) # expressed as vector in local initial orientation
    return(cols)

def rd_u(struCase, **kwargs):
    """
    Extracts generalized displacements data from *.p11 bin file and imports data to "stru" class object also creates "u_mdr" field if necessary
    
    Inputs:
        struCase: 'stru' class object
    **kwargs (may contain):
        subDir_P11: str, P11 file dir
    
    Returns:
        struCase: 'stru' class object, updated
    """
    if 'glob_print_output' in kwargs:
        glob_print_output = kwargs.get('glob_print_output')
    else:
        glob_print_output = False
    if 'subDir_P11' in kwargs:
        subDir_P11 = kwargs.get('subDir_P11')
    else:
        subDir_P11=''
    
    glob_u_raw = []
    glob_u_mdr = []
    loc_unf_raw = []
    
    glob_step_Nan = 0
    for node_name in struCase.nodes:
        callbat(subDir_P11,struCase.p11FN, str(node_name), "0", "temp_file")
        loc_lines = open(subDir_P11+"temp_file",'r')
        loc_x_dat = loc_lines.readlines()
        loc_table_raw, loc_step_Nan = rd_SimpactTable(loc_x_dat,0, STable_mode='with_counter')
        if loc_step_Nan > glob_step_Nan:
            glob_step_Nan = loc_step_Nan
            if glob_print_output:
                print('t_Nan update:', loc_step_Nan,', for node: ',node_name)
        #Save all before NaN filter
        loc_unf_raw.append(loc_table_raw)
        loc_lines.close()
    
    loc_unf_raw = NaN_filter(loc_unf_raw, glob_step_Nan, **kwargs) #Delete all t>t_Nan data
    for a in range(len(struCase.nodes)): #Continue. Esto debería reproducir nuevamente el bucle interrumpido por el filtrado
        loc_table_raw = loc_unf_raw[a]
        loc_table_mdr = np.copy(loc_table_raw)
        if struCase.rdof: # if rotational DoFs, create field u_mdr
            loc_table_mdr[:,4:]= rotModalDec(loc_table_mdr[:,4:])
        if a == 0:
            glob_time = loc_table_raw[:,0]
            total_ntime = len(glob_time)
            glob_u_raw = np.transpose(loc_table_raw)[1:,:]
            glob_u_mdr = np.transpose(loc_table_mdr)[1:,:]
        else:
            glob_u_raw = np.append(glob_u_raw,np.transpose(loc_table_raw)[1:,:],axis=0)
            glob_u_mdr = np.append(glob_u_mdr,np.transpose(loc_table_mdr)[1:,:],axis=0)
    
    
    os.remove(subDir_P11+"temp_file")
    
    struCase.nt = total_ntime
    struCase.t  = glob_time
    struCase.u_raw = glob_u_raw
    struCase.u_mdr = glob_u_mdr
    
    return struCase

def NaN_filter(full_data, Nan_step, **kwargs):
    """
    Deletes all data with index > earliest NaN
    
    Inputs:
        bulk_data: list of numpy arrays, each one is a stru disp. case
        Nan_step: int, Nan´s first index
    Returns:
        full_data, list of numpy arrays, delimited
    """
    if 'glob_print_output' in kwargs:
        glob_print_output = kwargs.get('glob_print_output')
    else:
        glob_print_output = False
    
    if Nan_step == 0:
        return(full_data)
    else:
        
        for i in range(len(full_data)):
            if Nan_step < len(full_data[i]):
                full_data[i] = np.delete(full_data[i],range(Nan_step,len(full_data[i])))
                if glob_print_output:
                    print(len(full_data[i])-Nan_step, ' steps deleted')
            else:
                if glob_print_output:
                    print('No data was deleted')
        return(full_data)
    
def rd_eqInfo(struCase, **kwargs):
    """
    Extracts Information Relative to the Equations Numbers from ASCII *.rsn file (Simpact or Delta output, default Delta)
    
    Inputs:
        struCase: 'stru' class object
    **kwargs (may contain):
        c_info_eq: int, line indicator for info eqs (in .rsn file), default 1
        
    Returns:
        struCase: 'stru' class object, updated
    """
    
    if 'subDir_RSN' in kwargs:
        subDir_RSN = kwargs.get('subDir_RSN')
    else:
        subDir_RSN=''
        
    if 'c_info_eq' in kwargs:
        c_info_eq = kwargs.get('c_info_eq')
    else:
        c_info_eq=1
    
    if not struCase.rsnDe:
        struCase.rsnDe = struCase.rsnSi
    y = open(subDir_RSN+struCase.rsnDe+'.rsn', 'r')
    x_dat = y.readlines()
    locs = search_string(x_dat, "Information Relative to the Equations Numbers")
    struCase.eqInfo = rd_SimpactTable(x_dat, locs[0]+c_info_eq)
    y.close()
    
    return struCase


def rd_mass(struCase, **kwargs):
    """
    Extracts lumped mass matrix from ASCII *.rsn file (Simpact or Delta output, default Delta)
    
    Inputs:
        struCase: 'stru' class object
    **kwargs (may contain):
        subDir_RSN: str, .rsn file dir
        c_lumped_m: int, lumped matrix line indic. (.rsn file), default 1
    Returns:
        struCase: 'stru' class object, updated
    """
    
    if 'subDir_RSN' in kwargs:
        subDir_RSN = kwargs.get('subDir_RSN')
    else:
        subDir_RSN=''
    
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
    y = open(subDir_RSN+struCase.rsnDe+'.rsn', 'r')
    x_dat = y.readlines()

    # find reference line and read table
    locs_lumped = search_string(x_dat, "  Lumped Mass Matrix")
    raw_lumped_matrix = rd_SimpactTable(x_dat, locs_lumped[0]+c_lumped_m)
    struCase.mass = []
    struCase.iLabl = []
    for j in range(0, len(struCase.eqInfo[:, 0])): #File to file
        for a in range(0, struCase.nnode): #Node to node (of interest)
            if struCase.eqInfo[j, 0] == struCase.nodes[a]:
                # save data
                struCase.iLabl.append([struCase.nodes[a], struCase.eqInfo[j,-1]]) # NOTA: hay que modificar esto para que lea solo la última columna, así queda funcionando también apra leer los *.rsn de Simpact
                if len(raw_lumped_matrix[j]) == 6:
                    struCase.mass = np.append(struCase.mass, raw_lumped_matrix[j])
                else:
                    struCase.mass = np.append(struCase.mass, raw_lumped_matrix[j,1:])
                if glob_print_output:
                    print("Lumped Mass Matrix")
                    print(struCase.mass)
    
    # close file
    y.close()
    
    struCase.iLabl = np.array(struCase.iLabl) + struCase.intLabOffset
    
    return struCase


def rd_eig(struCase, **kwargs):
    """
    Extracts eigen modes and eigen frequencies from ASCII *.rsn file (Delta output)
    Inputs:
        struCase: 'stru' class object
    **kwargs (may contain):
        subDir_RSN: str, .rsn file dir
        c_modes: int, line indicator, default 6
    Returns:
        struCase: 'stru' class object, updated
    """
    
    if 'subDir_RSN' in kwargs:
        subDir_RSN = kwargs.get('subDir_RSN')
    else:
        subDir_RSN=''
    
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
    y = open(subDir_RSN+struCase.rsnDe+'.rsn', 'r')
    x_dat = y.readlines()

    # find reference line and read table
    locs_modes = search_string(x_dat, "linear dynamic eigen-mode analysis")
    struCase = clean_eqInfo(struCase)    
    struCase.nm = len(locs_modes)
    struCase.om = np.zeros((struCase.nm, 3))
    struCase.phi = np.zeros((struCase.nm, struCase.nnode*6))
    for i in range(struCase.nm):
        struCase.om[i] = line_spliter(x_dat[locs_modes[i]])
        raw_mode_table = rd_SimpactTable(x_dat, locs_modes[i]+c_modes)
        if glob_print_output:
            print("Full mode table")
            print(raw_mode_table)
        local_mode_row = []
        for j in range(0, len(raw_mode_table[:, 0])): # number of rows in modes' table
            for a in range(0, struCase.nnode):
                if struCase.eqInfo[j, 0] == struCase.nodes[a]:
                    if len(raw_mode_table[j]) == 6:
                        local_mode_row = np.append(local_mode_row, raw_mode_table[j,:])
                    else:
                        local_mode_row = np.append(local_mode_row, raw_mode_table[j,1:])
                    if glob_print_output:
                        print("Local mode row ---------------")
                        print(local_mode_row)
        struCase.phi[i] = local_mode_row
    
    # close file
    y.close()
    
    struCase.phi = np.transpose(struCase.phi)
    
    return struCase

def ae_Ftable(struCase, **kwargs): ##NOTA: Si el nombre no gusta, lo cambio
    '''
    Extracts loads from .DAT files
    
    Inputs:
        struCase: 'stru' class obj
    **kwargs (may contain):
        subDir_FCS: str, .dat file loc
    Returns:
        struCase: 'stru' class obj, updated
    '''   
    if 'subDir_FCS' in kwargs:
        subDir_FCS = kwargs.get('subDir_FCS')
    else:
        subDir_FCS=''
    if 'glob_print_output' in kwargs:
        glob_print_output = kwargs.get('glob_print_output')
    else:
        glob_print_output = False
    
    x = open(subDir_FCS+struCase.loadsFN+'.dat','r') #Read file
    x_dat = x.readlines() #Read lines
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
                print("End of table, line: ", act_line)
            tab_cond = False
    counter = 0
    for a in loc_ftab: #For each \Delta t
        loc_frow_filt = np.array([])
        for i in range(len(a)): #For each row
           for j in range(len(struCase.iLabl[:,0])): #For each label (node of interest)
               if a[i,0] == struCase.iLabl[j,1]:
                   loc_frow_filt = np.append(loc_frow_filt, a[i,1:])
        if counter == 0:
           loc_ftab_filt = np.array([loc_frow_filt])
        else:
            loc_ftab_filt = np.append(loc_ftab_filt, [loc_frow_filt],axis=0)
        counter = counter + 1
    loc_ftab_filt = np.transpose(loc_ftab_filt)
    loc_ittab = np.array(loc_ittab)
    struCase = FTable_fit(struCase, loc_ftab_filt, loc_ittab[:,1], **kwargs)
    return struCase



def FTable_fit(struCase, y_loads, t_loads):
    '''
    Fits the aLoads to the stru-shaped time arr
    Inputs:
        struCase: 'stru' class obj
        y_loads: numpy array or table, loads data
        t_loads: numpy array, time data (loads)
    Returns:
        struCase: 'stru' class obj, updated
    '''
    final_y = []
    t_stru = struCase.t
    pre_count = 0
    len_t_stru = len(t_stru)
    for i in range(1,len(t_loads)):
        corr_count = 0
        if not pre_count + corr_count > len_t_stru:
            while t_stru[pre_count+corr_count] < t_loads[i]:
                final_y.append(y_loads[:,i-1])
                corr_count +=1
        pre_count = pre_count+corr_count
    final_y.append(y_loads[:,-1]) #NOTA: Esto supone que terminan en igual tf
    struCase.aLoad = np.transpose(final_y)
    return(struCase)



def rdBin(file, **kwargs):
    """
    Reads bin file
    
    Inputs:
        file: str, file name without extension
    **kwargs (may contain):
        subDir_BIN or subDir_P11: str, dir
        glob_print_output: bool, for msg output
    Returns:
        data: general type
    """
    
    if 'subDir_BIN' in kwargs:
        subDir = kwargs.get('subDir_BIN')
    elif 'subDir_P11' in kwargs:
        subDir = kwargs.get('subDir_P11')
    else:
        subDir=''
    if 'glob_print_output' in kwargs:
        print_output = kwargs.get('glob_print_output')
    else:
        print_output = False
    try:
        if file[-4:] == '.sim':
            f = open(subDir+file,'rb')
        else:
            f = open(subDir+file+'.sim','rb')
    except:
        raise NameError('File not found: '+file)
    data = pickle.load(f)
    f.close()
    if print_output:
        print ('bin data file read (read_bin funct)')
    
    return(data)  

def svBin(data, **kwargs):
    '''
    Saves data to BIN file
    Inputs:
        data: 'sim' class object, variable to be saved
    kwargs (may contain):
        glob_print_output: bool - print msg
        subDir_BIN or subDir_P11: str - dir
    Returns:
        None
    '''
    
    if 'subDir_BIN' in kwargs:
        subDir = kwargs.get('subDir_BIN')
    elif 'subDir_P11' in kwargs:
        subDir = kwargs.get('subDir_P11')
    else:
        subDir=''
    if 'glob_print_output' in kwargs:
        print_output = kwargs.get('glob_print_output')
    else:
        print_output = False
    if isinstance(data, sim):
        f = open(subDir+data.fName+'.sim', 'wb')
        pickle.dump(data, f)
        f.close()
        if print_output:
            print ('bin data file saved (save_bin funct)')

def clean_eqInfo(struCase):
    '''
    Deletes non-useful nodes from eqInfo Table
    Inputs:
        struCase, 'stru' class obj
    Returns:
        strCase, 'stru' class obj, updated
    '''
    rem_inds = []
    for i in range(len(struCase.eqInfo)):
        if sum(struCase.eqInfo[i,1:-2])==0:
            rem_inds.append(i)
    struCase.eqInfo = np.delete(struCase.eqInfo, rem_inds,axis=0)
    return(struCase)


"""
Funcs from análisis.py
"""

def mult_int(struCase, **kwargs):
    '''
    Searchs for integer multiples in data, computing a[i] % a[:]. Also, saves the original a[i] value in the main diag
    Inputs:
        struCase, 'stru' class obj
    **kwargs (may contain):
        data_type: str, Type of attri, default 'om'
        new_name: str, New name for the attr, default data_type + '_mults_matrix'
        extra_inds: lst or int, Extra inds for tables, if req.
    Returns:
        struCase, 'stru' class obj, updated
    '''
    if 'data_type' in kwargs:
        data_type = kwargs.get('data_type')
    else:
        data_type = 'om'
    if 'new_name' in kwargs:
        new_name = kwargs.get('new_name')
    else:
        new_name = data_type +'_mults_matrix'
    
    if data_type == 'om':
        data = getattr(struCase, data_type)[:,1]
    if 's_inds_0' in kwargs:
        s_inds_0 = kwargs.get('s_inds_0')
        data = data[s_inds_0]
    mults_matrix = []
    if not type(data) == np.ndarray:
        data = np.array(data)
    for i in range(len(data)):
        mults_matrix.append(data/data[i])
    setattr(struCase, new_name, np.array(mults_matrix)) 
    
    return(struCase)

def geomDOF_comp(struCase, dofDict, **kwargs):
    '''
    Computes the modal composition u = \Phi \cdot q but element-wise
    Inputs:
        struCase, 'stru' class obj
        dofDict, dict: {'node':[DOF]} - Desired geom. DOF
    **kwargs (may contain):
        norm, bool: Default False (norm using u_mdr)
        new_name, str: New attr name, default: 'qcomp_DOF'
        otp, str: 'default' or 'comp' - For internal use ('default' updates struCase)
    Returns:
        struCase, 'stru' class obj, updated
    '''
    if 'norm' in kwargs:
        norm = kwargs.get('norm')
    else:
        norm = False
    if 'new_name' in kwargs:
        new_name = kwargs.get('new_name')
    else:
        new_name = None    
    
    if 'otp' in kwargs:
        otp = kwargs.get('otp')
    else:
        otp = 'default'
    
    desired_ind = nodeDof2idx(struCase, dofDict)
    loc_phi = struCase.phiR[desired_ind]
    dec_u = np.multiply(np.transpose(loc_phi), struCase.q)
    if norm:
        dec_u *= 1/struCase.u_mdr[desired_ind]
    if new_name == None:
        new_name = 'qcomp_'+ str(list(dofDict.keys())[0]) + '_'+str(list(dofDict.values())[0][0]) + '_g' + str(desired_ind[0]+1)
    
    if otp == 'default':
        setattr(struCase, new_name, dec_u)
        return struCase
    elif otp == 'comp':
        return(dec_u)
    
def act_mINDS(struCase, dofDict, **kwargs):
    '''
    Determines active modal inds
    Inputs:
        struCase: 'stru' class obj
        dofDict: dict, {'node':[DOF]}
    **kwargs (may contain):
        tol: float, desired tolerance (default 1e-2)
        des_name: str, desired attr name (default: 'node_dof_amINDS')
    Returns:
        struCase, 'stru' class obj, updated
    '''
    if 'tol' in kwargs:
        tol = kwargs.get('tol')
    else:
        tol = 1e-3
    if 'des_name' in kwargs:
        des_name = kwargs.get('des_name')
    else:
        des_name = None
    
    dec_u = geomDOF_comp(struCase, dofDict, otp = 'comp')
    inds = []
    for i in range(dec_u.shape[0]):
        if np.max(dec_u[i]) >= tol:
            inds.append(i+1)
    
    if des_name == None:
        dof_name = str(list(dofDict.keys())[0])
        des_name = 'amINDS_' + dof_name + '_' +str(list(dofDict.values())[0][0]) + '_g' + str(nodeDof2idx(struCase, dofDict)[0]+1)
    if not inds == []:
        setattr(struCase, des_name, inds)
        struCase.part_active_modes.append({des_name:inds})
    else:
        print('Not active modes found - ', dofDict)
    return struCase

def sum_data(struCase, **kwargs):
    '''
    Sums some stru.data[inds,:] and saves it as a new attr
    Inputs:
        struCase: 'stru' class obj
    **kwargs (may contain):
        data_type: str, data to sum
        inds: list, absolute-inds pos (not Python´s)
        sum_name: str, stru.sum_name = sum(QW[inds,:])
    Returns:
        struCase, 'stru' class obj, updated
    '''
    if 'data_type' in kwargs:
        data_type = kwargs.get('data_type')
    else:
        raise NameError('Data type str error')
    if 'inds' in kwargs:
        inds = kwargs.get('inds')
    else:
        raise NameError('inds error')
    if 'sum_name' in kwargs:
        sum_name = kwargs.get('sum_name')
    else:
        raise NameError('sum_name erorr')
        
    for i in range(len(inds)):
        try:
            if i == 0:
                y = getattr(struCase,data_type)[inds[i]-1,:]
            else:
                y = y + getattr(struCase,data_type)[inds[i]-1,:]
        except:
            raise NameError('Wron data_type, pls check')
    setattr(struCase, sum_name, y)
    return(struCase)

def case_tag(**kwargs):
    '''
    Creates fnames and tags
    Inputs:
        None
    **kwargs (may contain):
        'case_type': str, 'R' (rigid), '' ()
        'case_vel': str, 'xxxx' (vel)
        'vel_otp': str, '.' or '' (default)
    Returns:
        tag: str
    '''
    tag = ''
    if 'case_type' in kwargs:
        case_type = kwargs.get('case_type')
    else:
        case_type = ''
    if 'case_vel' in kwargs:
        case_vel = kwargs.get('case_vel')
    else:
        case_vel = ''
    if 'vel_otp' in kwargs:
        vel_otp = kwargs.get('vel_otp')
    else:
        vel_otp = ''
    tag += case_type
    if not case_type == '':
        tag += '_'
    if vel_otp == '.':
        try:
            case_vel = float(case_vel)
            case_vel = case_vel // 10
            case_vel = str(case_vel)
        except:
            print('Bad vel to float')
    tag += case_vel
    if tag =='':
        tag = 'no_name'
    return(tag)

def amp_search(struCase, **kwargs):
    '''
    Determines the amplitude of a signal
    Inputs:
        struCase, 'stru' class obj or single data: ndarray
    **kwargs (may contain):
        mode: str, 'normal' or 'object' (default 'object')
            if 'object':
                data_type, str - signal attr name
                pos_ind, int - row selector ind, real, not Python´s
        start_type: int, start index (default -1)
        reverse: int, search dir (default -1, end > start)
        otp: str, returns just 'amp' (amplitude) or 'full' (amplitude, inds), 'inds' (inds), 'obj_attr' (attr, default)
    Returns:
        amplitude or inds or mid_value (or all of them) or a dict in obj.attr {'amplitude', 'mid_value', 'inds'}
    '''
    if 'mode' in kwargs:
        mode = kwargs.get('mode')
    else:
        mode = 'object'
    
    if mode == 'object':
        if 'data_type' in kwargs:
            data_type = kwargs.get('data_type')
        else:
            raise NameError('No data_type (attr)')
        if 'pos_ind' in kwargs:
            pos_ind = kwargs.get('pos_ind')
        else:
            pos_ind = None
        
        try:
            y = getattr (struCase,data_type)
            if not pos_ind == None:
                y = y[pos_ind-1]
        except:
            raise NameError('Wrong data type or pos ind')

    if 'start_type' in kwargs:
        start_type = kwargs.get('start_type')
    else:
        start_type= -1
    if 'reverse' in kwargs:
        reverse = kwargs.get('reverse')
    else:
        reverse = -1
    
    if 'otp' in kwargs:
        otp = kwargs.get('otp')
    else:
        otp = 'obj_attr'
    if 'use' in kwargs:
        use = kwargs.get('use')
    else:
        use = 'env'
    
    if start_type == -1:
        start_ind = len(y)-2
    elif start_type == 0:
        start_ind = 1
    else:
        start_ind = start_type
    
    inds = []
    
    if use == 'loc_extr':
        #max local    
        while y[start_ind-1] >= y[start_ind] or y[start_ind+1] >= y[start_ind]:
            start_ind += reverse
        inds.append(start_ind)
        #Min local
        while y[start_ind-1] <= y[start_ind] or y[start_ind+1] <= y[start_ind]:
            start_ind += reverse
        inds.append(start_ind)
    elif use == 'env':
        lmin,lmax = hl_envelopes_idx(y)
        inds.append(lmax[-2])
        inds.append(lmin[-2])
        
    amplitude = abs(y[inds[1]]-y[inds[0]])
    mid_value = y[inds[1]]+amplitude/2
    
    if otp == 'full':
        return(amplitude, mid_value, inds)
    elif otp == 'amp':
        return(amplitude)
    elif otp == 'inds':
        return(inds)
    elif otp == 'mid_value':
        return(mid_value)
    elif otp == 'obj_attr':
        setattr(struCase, data_type + '_amp_data', {'amplitude':amplitude, 'mid_value':mid_value, 'inds':inds})
        return(struCase)
    
def hl_envelopes_idx(y, dmin=1, dmax=1, split=False):
    """
    Extracts envelopes, both high and low
    
    Inputs:
        y: 1D numpy array, data
        dmin, dmax: int, optional, size of chunks, use this if the size of the input signal is too big
        split: bool, optional, if True, split the signal in half along its mean, might help to generate the envelope in some cases
    Returns_
        lmin,lmax : high/low envelope idx of input signal s
    """
 
    lmin = (np.diff(np.sign(np.diff(y))) > 0).nonzero()[0] + 1 
    lmax = (np.diff(np.sign(np.diff(y))) < 0).nonzero()[0] + 1 
    if split:
        # s_mid is zero if s centered around x-axis or more generally mean of signal
        s_mid = np.mean(y) 
        # pre-sorting of locals min based on relative position with respect to s_mid 
        lmin = lmin[y[lmin]<s_mid]
        # pre-sorting of local max based on relative position with respect to s_mid 
        lmax = lmax[y[lmax]>s_mid]
    lmin = lmin[[i+np.argmin(y[lmin[i:i+dmin]]) for i in range(0,len(lmin),dmin)]]
    lmax = lmax[[i+np.argmax(y[lmax[i:i+dmax]]) for i in range(0,len(lmax),dmax)]]
    
    return lmin,lmax
    
def handle_act_modes(struCase, **kwargs):
    '''
    Creates a single list of inds (real, not Python´s) for active and pasive modes.
    Inputs:
        struCase, 'stru' class obj
    **kwargs (may contain):
        top, int - Max pasiv modal (real, not Python´s), default -1 (last)
    Returns:
        struCase, 'stru' class obj, updated (active and pasive _minds)
    '''
    if 'top' in kwargs:
        top = kwargs.get('top')
    else:
        top = -1
    activ_inds = []
    pasiv_inds = []
    if not struCase.part_active_modes == []:
        for loc_dct in struCase.part_active_modes:
            for loc_ind in list(loc_dct.values())[0]:
                if not loc_ind in activ_inds:
                    activ_inds.append(loc_ind)
    
    if struCase.moi == []:
        for loc_ind in np.arange(1,struCase.nm+1,1):
            if top > loc_ind:
                if not loc_ind in activ_inds:
                    pasiv_inds.append(loc_ind)
    else:
        for loc_ind in struCase.moi:
            if top > loc_ind:
                if not loc_ind in activ_inds:
                    pasiv_inds.append(loc_ind)

    setattr(struCase, 'active_minds', activ_inds)
    setattr(struCase, 'pasive_minds', pasiv_inds)
    return struCase


"""
General tools from análisis.py
"""

def lst_av_dirs(path):
    '''
    Lists available cases
    
    Inputs:
        path: str, abs path to data dir
    Returns:
        lst: list, subdir names
    '''
    if path == '':
        return(os.listdir())
    else:
        try:
            return(os.listdir(path))
        except:
            raise NameError('bad path')

def delete_av_bins(path, **kwargs):
    '''
    Deletes all avaiable bin files
    Inputs:
        path: str, global or rel path
    **kwargs (may contain):
        filecode: str, 'R' or 'D' (only those, default: all)
    Returns
        None. Files are gone
    '''
    if 'filecode' in kwargs:
        filecode = kwargs.get('filecode')
    else:
        filecode = ''
    
    av_files = lst_av_dirs(path, **kwargs)
    for loc_file in av_files:
        if not filecode == '':
            if loc_file[0] == filecode:
                os.remove(path+loc_file)
        else:
            os.remove(path+loc_file) 
    return None

"""
Ev cases from eigpp and análisis
"""

def rec_cases(av_cases, dirs, **kwargs):
    """
    Simple recursive function, calls eigen_an over a set of cases
    
    Inputs:
        av_cases, list of str - p11 dir
        dirs, dict - files location
        **kwargs
    Returns:
        None

    """
    for loc_case in av_cases:
        eigen_an(loc_case, dirs, **kwargs)
    return None

def eigen_an(loc_case, dirs, **kwargs):
    """
    Simplified general analysis. Creates case and .sim file. Reads and processes displacement and load files. Performs modal decomposition. Exports to .BIN
    
    Inputs:
        loc_case: str, available p11 file name
        dirs: dict - assoc. files names and location
            keys: 'r_p11_subdir', 'r_RSN_subdir', 'r_FCS_fname', 'bin_path'
    **kwargs (may contain):
        general kwargs for downstream funcs
        (must contain):
        case_type: str, R, D or S.
        case_extra_info: dict, extra info for loc_sim gen
            nodes
            intLabOffset
            sti
    Returns:
        None
    """
    if 'case_type' in kwargs:
        case_type = kwargs.get('case_type')
    else:
        raise NameError('Warning: No case type - cmd will fail')
        case_type = 'NC'
    if not 'case_extra_info' in kwargs:
        print('Warning: No case extra info')
        case_extra_info = None
    else:
        case_extra_info = kwargs.get('case_extra_info')
    
    loc_sim = sim()
    loc_sim.name = case_tag(case_type=case_type,case_vel = loc_case[-4:], vel_otp ='.')
    loc_sim.fName = case_tag(case_type=case_type,case_vel = loc_case[-4:])
    loc_sim.descr = 'v_\infty =' + case_tag(case_vel = loc_case[-4:], vel_otp ='.') +'\, ft/s'
    
    if case_extra_info:
        for loc_key in kwargs['case_extra_info']:
            setattr(loc_sim.stru, loc_key, kwargs['case_extra_info'][loc_key])
    
    loc_sim.stru.p11FN = 'pcolgante.@1' #binary .p11 fname
    loc_sim.stru = rd_u(loc_sim.stru, **{'subDir_P11':dirs['r_path']+dirs['r_p11_subdir']+loc_case+'/'})
    loc_sim.stru.rsnDe= 'pcolgante' #Temporal
    loc_sim.stru = rd_rsn_De(loc_sim.stru, **{'subDir_RSN':dirs['r_path']+dirs['r_RSN_subdir']})

    if dirs['r_FCS_fname'] in lst_av_dirs(dirs['r_path']+dirs['r_p11_subdir']+loc_case+'/'):
        loc_sim.stru.loadsFN = dirs['r_FCS_fname'][:-4]
        loc_sim.stru = ae_Ftable(loc_sim.stru,**{'subDir_FCS':dirs['r_path'] + dirs['r_p11_subdir']+loc_case+'/'})
        print('FCS file analized', loc_case)
    else:
        print('FCS file not found, case:', loc_case)
    loc_sim.stru = modalDecomp(loc_sim.stru, **{'subDir_BIN':dirs['bin_path']})
    svBin(loc_sim, **{'subDir_BIN': dirs['bin_path']})
        
    return None

"""
------------------------------------------------------------------------------
Working examples
------------------------------------------------------------------------------

#Read and gen R-case

#File dirs data
direc = {'r_path':'D:/Archivos/UNC/ayudantia/data/020 aeroelástica vigas rígido/'} #global or rel path to p11 and RSN (rigid)
direc['r_p11_subdir'] = '020 10 08x40/' #P11(s) subdir
direc['r_RSN_subdir'] = 'frecuencias y modos/' #RSN subdir
direc['r_FCS_fname'] =  'AeroFcsOnStruc.dat' #Searchs for this file, use .dat!
direc['bin_path'] = 'D:/Archivos/UNC/ayudantia/data/BINS/'

#Extra data
case_extra = {}
case_extra['nodes'] = [200000, 200001]
case_extra['intLabOffset'] = 6
case_extra['sti'] = np.array([127300,127300,127300,221448514,221448514,221448514,0,0,0,0,0,0])

#Available files
av_cases = sim_db.lst_av_dirs(direc['r_path']+direc['r_p11_subdir'])

#Run
sim_db.eigen_an(av_cases[3], direc, case_extra_info = case_extra, case_type = 'R')
"""

"""
------------------------------------------------------------------------------
Unit tests
------------------------------------------------------------------------------
"""
def uTest1():
    pass
    

"""
------------------------------------------------------------------------------
runing things
------------------------------------------------------------------------------
"""
if __name__ == '__main__':
    
    pass
    
