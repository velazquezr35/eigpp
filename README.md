# eigpp tool

## sim_db.py

### classes

`class sim: `     full model info, container for simulation subparts

`
    def __init__(self):
    
        self.name  = ''                             # short description
        self.descr = ''                             # description
        
        self.fName = ''                             # binary file name (without extension)
        
        self.stru       = stru() # 1 "stru" class objects - active
        self.aero       = aero() # 1 "aero" class objects - active
        
        self.struList   = []    # list of "stru" class objects
        self.aeroList   = []    # list of "aero" class objects`

`class stru: ` structural model info

`  
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
        self.rot_inds = [4,5,6]                             # rotational DOFs inds (not Python´s)`

`class aero: ` aerodynamic model info

### funcs

`upd_phiR(struCase): `    Updates struCase.phiR and struCase.omR.

`upd_modalMK(struCase): `    Determines struCase.mmass and struCase.mstif from struCase.phiR and struCase.mass.

`upd_mnorm(struCase): `    Updates struCase.phiR, struCase.q, struCase.Q,

`modalDecomp(struCase,**kwargs): `    Applies modal decomposition over stru´s DOFs and LOADS (if available).

`modal_mechEnergy(struCase): `    Determines mechanical energy related to each mode.

`search_time(t_array, t_values, **kwargs): `    Searchs for indexes (equal or max) in a t_array corresponding to some t_values

`sfti_time(struCase, *reset, **kwargs): `    Searchs for time indexes - Updates the desired indexes for plotting purposes

`nodeDof2idx(struCase, nodeDOFs): `    Returns indexes for the nodes and nodes DOFs of interest

`modalDof2idx(struCase, modalDOF): `    Returns modal indexes for a particular shape-DOF

`callbat(file_folder, file_name, nodo, dof, output_name): `    Extracts generalized displacements data from *.p11 bin file. Sends windows cmd commands.

`line_spliter(line): `    Splits a string line. For general use with SimpactTable. 

`search_string(x_dat, frase): `    Searchs for a keyword or string in a list of strings

`rd_SimpactTable(x_dat, start_line, **kwargs): `    Locates and reads a table from *.rsn ASCII Simpact and Alpha files

`rd_rsn_De(struCase, **kwargs): `    Reads data from Delta *.rsn output (eigen modes and eigen frequencies, mass matrix) and performs modal decomposition.

`euler2axial(cols): `    Converts rotations expresed as 3-1-3 Euler Angles to axial vector form using SciPy class Rotation

`rotModalDec(cols): `    Prepares rotational data for modal decomposition. Represent rotations as incremental rotation vectors expressed in initial local system.

`rd_u(struCase, **kwargs): `    Extracts generalized displacements data from *.p11 bin file and imports data to "stru" class object also creates "u_mdr" field if necessary

`NaN_filter(full_data, Nan_step, **kwargs): `    Deletes all data with index > earliest NaN

`rd_eqInfo(struCase, **kwargs): `    Extracts Information Relative to the Equations Numbers from ASCII *.rsn file (Simpact or Delta output, default Delta)

`rd_mass(struCase, **kwargs): `    Extracts lumped mass matrix from ASCII *.rsn file (Simpact or Delta output, default Delta)

`rd_eig(struCase, **kwargs): `    Extracts eigen modes and eigen frequencies from ASCII *.rsn file (Delta output)

`ae_Ftable(struCase, **kwargs): ##NOTA: Si el nombre no gusta, lo cambi: `    Extracts loads from .DAT files

`FTable_fit(struCase, y_loads, t_loads): `    Fits the aLoads to the stru-shaped time arr

`rdBin(file, **kwargs): `    Reads bin file

`svBin(data, **kwargs): `    Saves data to BIN file

`clean_eqInfo(struCase): `    Deletes non-useful nodes from eqInfo Table

`mult_int(struCase, **kwargs): `    Searchs for integer multiples in data, computing a[i] % a[:]. Also, saves the original a[i] value in the main diag

`geomDOF_comp(struCase, dofDict, **kwargs): `    Computes the modal composition u = \Phi \cdot q but element-wise

`act_mINDS(struCase, dofDict, **kwargs): `    Determines active modal inds

`sum_data(struCase, **kwargs): `    Sums some stru.data[inds,:] and saves it as a new attr

`case_tag(**kwargs): `    Creates fnames and tags

`amp_search(struCase, **kwargs): `    Determines the amplitude of a signal

`hl_envelopes_idx(y, dmin=1, dmax=1, split=False): `    Extracts envelopes, both high and low

`handle_act_modes(struCase, **kwargs): `    Creates a single list of inds (real, not Python´s) for active and pasive modes.

`lst_av_dirs(path): `    Lists available cases

`delete_av_bins(path, **kwargs): `    Deletes all avaiable bin files

`rec_cases(av_cases, dirs, **kwargs): `    Simple recursive function, calls eigen_an over a set of cases

`eigen_an(loc_case, dirs, **kwargs): `    Simplified general analysis. Creates case and .sim file. Reads and processes displacement and load files. Performs modal decomposition. Exports to .BIN

(Created using autodoc.py)

## plotter.py

### funcs

`plt_temp_ev(struCase, data_indic, ax, **kwargs): `    Plot DOFs (u_mdr as default) or FCs (or modal coordinates or modal ext. loads) as f(t)

`plt_tfixed(struCase, data_indic, ax, **kwargs): `    Plot all dof coords (or ext. loads) or modal coords (or modal ext. loads) in particular instants of time.

`plt_phi(struCase, modedofDict, ax, **kwargs): `    Plots modal shapes

`plt_FFT(struCase, data_indic, ax, **kwargs): `    Plots the FFT of a signal: u-DOF (u_mdr as default) (or FC) or q (or Q)

`plt_PP(struCase, data_indic, ax, **kwargs): `    Plots phase-plane portraits, du/dt vs u or dFCS/dt vs FCS (DOF an) or q vs dq/dt or Q vs dQ/dt (modal an)

`plt_spectr(struCase, data_indic, fig, ax, **kwargs): `    Plots spectrograms of u(t) (mdr def) or F(t) or q(t) or Q(t)

`plt_uxuy(struCase, vsDict, ax, **kwargs): `    Plots two DOFs or FCS for all struCase nodes for a single t val (one per curve)

`plt_general(struCase, inds, ax, **kwargs): `    Plots general props as f(t) using some inds (not Python´s)

`fig_general(struCase, indLIST, **kwargs): `    Arranges plots of general props as f(t) using some inds (not Python´s)

`fig_uxuy(struCase,vsLIST, **kwargs): `    Arranges plots of DOF vs DOF or FC vs FC @t fixed

`fig_uqs(struCase, DATA_indc, **kwargs): `    Arranges plots of DOF(nodes) or FC(nodes) (or modes (nodes) or modal ext. loads (nodes)) @t fixed

`fig_uqt(struCase, DATA_indc, **kwargs): `    Arranges plots of u(t) or load(t) (or q(t) or Q(t))

`fig_FFT(struCase, DATA_indic, **kwargs): `    Arranges plots of FFT(u(t)) or FFT(load(t)) or FFT(q(t)) or FFT(Q(t))

`fig_spect(struCase, DATA_indic, **kwargs): `    Arranges plots of Spectrogram(u(t)) or Spectrogram(load(t)) or Spectrogram(q(t)) or Spectrogram(Q(t))

`fig_PP(struCase, DATA_indic, **kwargs): `    Arranges plots of PP(u(t)) or PP(load(t)) or PP(q(t)) or PP(Q(t))

`fig_phi(struCase, modedofLIST, **kwargs): `    Arranges plots of modal shapes

`lst2str(lst, **kwargs): `    Generates a str from a lst

`dof2dofDict(struCase, dof): `    Generates a dofDict for a desired DOF using all the available nodes in struCase

`tdof2dofDict(struCase, tdof_dict): `    Generates a dofDict for a desired DOF using all the available nodes in struCase

`nodes2labels(struCase, **kwargs): `    Generates a list of strings from struCase.nodes

`tldxs(struCase, desired_t): `    Generates a list of indexes for the closest struCase.t time values in t arg (list, floats)

`keys_to_str(dct_keys): `    Simple tool for str generation from dict_keys data

`flatten_values_list(dct_values): `    Simple tool for fusion nested lists (coming from dicts_values)

`label_asoc(dct): `    Simple tool for generate a list of labels from dict keys

`handle_graph_info(**kwargs): `    Prepares labels, titles and general things for plots

`save_figure(fig, opts,**kwargs): `    Save plots in desired dir

`handle_modal_inds(struCase, modal_inds, **kwargs): `    Looks for the real modal index in struCase.moi, for labeling purposes

`hide_axes(ax, **kwargs): `    Auto ax hide

`do_grid_and_labels(ax, **kwargs): `    Auto grid and legend

`FFT_labl_gen(data_type, i, original_inds, node_labels): `    Gen a custom label for FFT plots

`gen_aloneinds_prop(inds): `    Generates a nested list for indsLIST kwarg to use in fig_general, using a single list of inds.

`general_envs(struCase,**kwargs): `    Gens envelopes - Plots (max(f[i])) vs some inds

`add_ticks(ax,**kwargs): `    Add ticks to an ax obj

`gen_single_prop(desired_shape, prop): `    Generates a list for data_type kwarg to use in fig_general, using a single struCase.attr.

Created using autodoc.py ;)
