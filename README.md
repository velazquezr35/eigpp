# sim_db.py

## Docs

### Classes

`sim: ` structural model info

`stru: ` full model info. container for simulation subparts

`aero: `  aerodynamic model info

### Funcs

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

`ae_Ftable(struCase, **kwargs): `    Extracts loads from .DAT files

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

-- Doc extracted using auto_doc.py --
