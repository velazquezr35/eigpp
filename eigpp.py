# -*- coding: utf-8 -*-
"""
Created on Thu Feb 25 17:52:06 2021

@author: Rodrigo Velazquez

Ayudantía 2021 - Departamento de Estructuras FCEFyN
"""

#Zona para importar libs y módulos
import pickle
import numpy as np
import re
import os
from scipy.spatial.transform import Rotation as Rot
import matplotlib.pyplot as plt

data_folder = 'data/'
bin_folder = 'ie_data/'

glob_print_output = True

# Class offsets - Modificar para posibles parts. a leer

class offsets:
    c_info_eq = 1
    c_lumped_m = 1
    c_modes = 6

#------------------------------------------------------------------------------
# Definición de funciones

##Funciones originales para el análisis estructural

def search_string(x_dat, frase):
    lin_count = 0  # contador de lineas
    locs = []  # ubicaciones dónde encuentra (line number python = -1 realidad)
    for txt in x_dat:
        lin_count = lin_count + 1
        x = re.search(frase, txt)
        if x != None:
            locs.append(lin_count)
            #print(lin_count)
    # devolvemos la lista con las lineas de x_dat que tienen la coincidencia
    return(locs)

def rd_SimpactTable(x_dat, start_line, **kwargs):
    #función que devuelve una tabla asociada genérica. Pensada para los modos

    # Directamente trabajamos con start_line dándo este valor como offset al llamarla
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


def rd_rawData(fname, nodos_interes):
    if glob_print_output:
        print("Raw data análisis completo")
    return *rd_eig(fname,nodos_interes),*rd_u(loc_sim.stru)

def rd_eig(fname, nodos_interes):
    
    # leemos archivo así queda cargado y evitamos re-leer
    y = open(fname, 'r')
    x_dat = y.readlines()

    # Buscamos la linea que tiene la info de los nodos
    locs = search_string(x_dat, "Information Relative to the Equations Numbers")
    raw_info_matrix = rd_SimpactTable(x_dat, locs[0]+offsets.c_info_eq)
    # buscamos la matriz de masa
    locs_lumped = search_string(x_dat, "  Lumped Mass Matrix")
    raw_lumped_matrix = rd_SimpactTable(x_dat, locs_lumped[0]+offsets.c_lumped_m)
    n_rows = len(nodos_interes)
    M_i_matrix = []
    info_matrix = []
    for j in range(0, len(raw_info_matrix[:, 0])): #desplazamos fila a fila
        for a in range(0, n_rows): #desplazamos nodo a nodo (de interés)
            if raw_info_matrix[j, 0] == nodos_interes[a]:
                #guardamos la info que nos interesa
                info_matrix.append([nodos_interes[a], raw_info_matrix[j,-1],raw_info_matrix[j,-2]])
                M_i_matrix = np.append(M_i_matrix, raw_lumped_matrix[j])
                if glob_print_output:
                    print("Matriz de masa")
                    print(M_i_matrix)

    # por último, loopear para leer todos los modos
    info_matrix = np.array(info_matrix)
    locs_modes = search_string(x_dat, "linear dynamic eigen-mode analysis")

    # guardamos orden,frecuencias y períodos
    m_O_f_p = np.zeros((len(locs_modes), 3)) #sim.str.omsim.str.om
    Modes_reshaped_rows = np.zeros((len(locs_modes), len(nodos_interes)*6))

    for i in range(len(locs_modes)):
        # guardamos O,f,P
        m_O_f_p[i] = line_spliter(x_dat[locs_modes[i]])
        # guardamos y adecuamos tablas tablas
        raw_mode_table = rd_SimpactTable(x_dat, locs_modes[i]+offsets.c_modes)
        if glob_print_output:
            print("tabla modo full")
            print(raw_mode_table)
            print(" modos------------------------------")
        local_mode_row = []
        for j in range(0, len(raw_mode_table[:, 0])): #cantidad de filas de 1 tabla modo
            for a in range(0, n_rows): #cantidad de nodos de interes
                if raw_mode_table[j, 0] == nodos_interes[a]: #si el nodo es el que quiero
                    local_mode_row = np.append(local_mode_row, raw_mode_table[j,1:]) #apilar la columna correspondiente del modo
                    if glob_print_output:
                        print("aporte local ---------------")
                        print(local_mode_row)
        Modes_reshaped_rows[i] = local_mode_row
        
    y.close()
    if glob_print_output:
        print("MATRIX INFO")
        print(info_matrix)
    return M_i_matrix, np.transpose(Modes_reshaped_rows), m_O_f_p[:,1], info_matrix

def rd_u(stru_clase, **kwargs):
    '''input: class stru \n
        kwargs: 
        output: '''
    glob_u_raw = []
    glob_u_avr = []
    for node_name in stru_clase.nodes:
        callbat(data_folder+"pcolgante.@1", str(node_name), "0", data_folder+"temp_file")
        loc_lines = open(data_folder+"temp_file",'r')
        loc_x_dat = loc_lines.readlines()
        loc_table_raw = rd_SimpactTable(loc_x_dat,0)
        loc_table_avr = np.copy(loc_table_raw)
        if stru_clase.rdof: #si hay grados de libertad, rotar
            loc_table_avr[:,4:]= euler2axial(loc_table_avr[:,4:])
        if stru_clase.nodes.index(node_name) == 0:
            #Acomodar posibilidad de que el tiempo de 1 nodo sea menor (NaN antes q resto)
            glob_time = loc_table_raw[:,0]
            total_ntime = len(glob_time)
            glob_u_raw = np.transpose(loc_table_raw)[1:,:]
            glob_u_avr = np.transpose(loc_table_avr)[1:,:]
        else:
            glob_u_raw = np.append(glob_u_raw,np.transpose(loc_table_raw)[1:,:],axis=0)
            glob_u_avr = np.append(glob_u_avr,np.transpose(loc_table_avr)[1:,:],axis=0)
        
    return glob_u_raw, glob_u_avr,glob_time,total_ntime

def euler2axial(cols):
    for i in range(len(cols[:,0])):
        R = Rot.from_euler('ZXZ',cols[i,:],degrees=False)
        cols[i,:] = R.as_rotvec()
    return(cols)
    
def callbat(arch_data, nodo, dof, nombre_otp):
    cmd = "(echo " + arch_data + " 1 0 && echo d && echo " + nodo + " " + dof + " " + nombre_otp + " && echo s) | curvas"  
    os.system(cmd)

def line_spliter(line):
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

def save_bin(archivo, data):
    # global data
    f = open(archivo, 'wb')
    pickle.dump(data, f)
    f.close()
    print ('Data saved - BIN file')
    
## Funciones agregadas para el análisis de cargas

def ae_Ftable(fname, noderefs):
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
------------------------------------------------------------------------------
Zona para definir opciones y clases
------------------------------------------------------------------------------
"""

#Definición modelo propuesta por consigna

# sim.str.mass:   reales, ngl         - matriz de masa diagonal del modelo
# sim.str.om:     reales, nm          - vector de frecuencias naturales ordeandas de menor a mayor
# sim.str.phi:    reales, ngl x nm    - matriz modal (cada columna es un modo)
# sim.str.q:      reales, nm=ngl x nt    - matriz con coordenadas modales a lo largo del tiempo, cada columna tiene el vector q de un instante
# sim.str.t:      reales, nt          - vector con instantes donde se conoce la solución (malla o grilla temporal)
# sim.str.nt:     entero, 1           - cantidad de instantes de tiempo
# sim.str.u_raw:  reales, ngl x nt    - matriz con desplazamientos a lo largo del tiempo, el vector u para un instante se guarda en una columna - como se obtienen de Simpact/curvas (si existen GL rotacionales, se expresan con ángulos de Euler y no se puede usar para hacer descomposición modal)
# sim.str.u_avr:  reales, ngl x nt    - como sim.str.u_raw pero transformando los GL rotacionales de ángulos de Euler a vector axial (avr: axial vector rotations)


class sim(): #Más fácil definir una clase que contenga todo lo de interés
    def __init__(self):
        self.name = "SIM_TESTS" #Nombre de la simulación (se usa luego para guardar archivos y demás)
    #subclase para análisis estructural
    class stru:
        
        def __init__(self):
            self.nodes = [200001,200003]  #etiquetas de nodos de interés - lista de enteros
            self.nnode = 2       #cantidad de nodos considerados
            self.ndof = 0        #cantidad de GL considerados (incluyendo los que son restringidos por condiciones de borde)
            self.strurdOpt = 'raw'  #bandera para leer datos: 
                    #'raw' desde archivos provenientes de Simpact y delta1
                    #'bin' a partir de archivo binario guardado previamente
            self.rdof= True     #bandera que indique si hay GL rotacinales
            self.strueigOpt = True   #bandera para hacer (o no) descomposición modal
            self.loadsrdOpt = 'raw' #bandera para leer datos de cargas
            self.loadseigOpt = True #bandera para descomponer cargas

        def resultados_loads(clase_aero, steps, t, fzas):
            clase_aero.steps = steps
            clase_aero.loads_t = t
            clase_aero.fzas = fzas
            clase_aero.loads_nt = int(steps[-1])
            clase_aero.loads_q = []

        def resultados(clase,mass, phi, om, refer, u_raw,u_avr,t,nt):
            clase.mass = mass
            clase.phi = phi
            clase.om = om
            clase.refer = refer
            clase.u_raw = u_raw
            clase.u_avr = u_avr
            clase.t = t
            clase.nt = nt
            clase.q = []
            
            
        def plot(self): #Pruebas para plotear
            fig, ax = plt.subplots()
            ax.plot(self.t,self.u_raw[:,0])

"""
------------------------------------------------------------------------------
Código Principal
------------------------------------------------------------------------------
"""

if __name__ == '__main__':
    #Vemos si calcular y exportar o directamente cargar análisis previo
    #generamos una clase para trabajar. Podríamos generar varias y analizar en serie.
    loc_sim = sim()
    loc_sim.stru = sim.stru()    
    #hacemos los cálculos, si corresponde:
    if loc_sim.stru.strurdOpt == 'raw':
        loc_sim.stru.resultados(*rd_rawData(data_folder+"red_pcolgante_2_modos.txt", loc_sim.stru.nodes))
        save_bin(bin_folder+"SIM_BN", loc_sim.stru) #Exportamos la info
    if loc_sim.stru.strurdOpt=='bin':
        f = open(bin_folder+"SIM_BN",'rb')
        loc_sim.stru = pickle.load(f)
        loc_sim.stru.strueigOpt = True
        print("Data loaded - BIN file")
        f.close()   
        
    #nuevamente, si corresponde:
    if loc_sim.stru.loadsrdOpt == 'raw':
        loc_sim.stru.resultados_loads(*ae_Ftable(data_folder+"AeroFcsOnStruc.dat",loc_sim.stru.refer))
        save_bin(bin_folder+"SIM_BN", loc_sim.stru) #Exportamos la info
    
    elif loc_sim.stru.loadsrdOpt =='bin':
        f = open(bin_folder+"SIM_BN",'rb')
        loc_sim.stru = pickle.load(f)
        loc_sim.stru.eigOpt = True
        print("Data loaded - BIN file")
        f.close()  
    
    
    #Ya tenemos la info cargada y lista. Ahora ver si corresponde hacer la descomposición, en sendos casos:
    if loc_sim.stru.strueigOpt:
        #Versión más eficiente (sirve luego p/descomponer cargas)
        aux = np.zeros((len(loc_sim.stru.phi[0,:]),len(loc_sim.stru.mass)))
        for i in range(len(aux)):
            aux[i] = np.multiply(loc_sim.stru.mass,loc_sim.stru.phi[:,i])
            
        #Descomponemos tiempo a tiempo
        for i in range(loc_sim.stru.nt):
            loc_sim.stru.q.append(np.matmul(aux, loc_sim.stru.u_avr[:,i]))
        
        loc_sim.stru.q = np.transpose(np.array(loc_sim.stru.q))
        
    #Idem para la descomposición de fzas
    if loc_sim.stru.loadseigOpt:
        if not aux.any(): #Repetir cálculo, si no reutilizamos
            aux = np.zeros((len(loc_sim.stru.phi[0,:]),len(loc_sim.stru.mass)))
            for i in range(len(aux)):
                aux[i] = np.multiply(loc_sim.stru.mass,loc_sim.stru.phi[:,i])
            #Descomponemos tiempo a tiempo
        for i in range(loc_sim.stru.loads_nt):
            loc_sim.stru.loads_q.append(np.matmul(aux, loc_sim.stru.fzas[:,i]))
        
        loc_sim.stru.loads_q = np.transpose(np.array(loc_sim.stru.loads_q))
    
def test_aux_vs_norm():
    '''Función para testear diferencias entre los 2 métodos propuestos para la descomposición modal'''
    listado = np.linspace(100, 10e3,15)
    tiempos_full = []
    tiempos_prev = []
    loc_sim.stru.nt = 1
    for j in range(len(listado)):
        size = int(listado[j])
        loc_sim.stru.mass = np.ones(size)
        loc_sim.stru.u_avr = np.ones((loc_sim.stru.nt,size))
        loc_sim.stru.phi = np.ones((size,6))
        loc_sim.stru.q = []
        
        if loc_sim.stru.eigOpt:
            import time
            t_1 = time.time()
            test = []
            
            for i in range(loc_sim.stru.nt):
                loc_u = loc_sim.stru.u_avr[i]
                loc_prod = np.matmul(np.diag(loc_sim.stru.mass),loc_u)
                loc_prod = np.matmul(np.transpose(loc_sim.stru.phi),loc_prod)
                test.append(loc_prod)
            
            dt_full = time.time()-t_1
            test = np.array(test)
            
            t_2 = time.time()
            #Versión más eficiente (sirve luego p/descomponer cargas)
            aux = np.zeros((size,len(loc_sim.stru.phi[0])))
            #Calculamos producto auxiliar
            for i in range(len(loc_sim.stru.phi)):
                aux[i] = loc_sim.stru.phi[i] * loc_sim.stru.mass[i]
                
            #Descomponemos tiempo a tiempo
            for i in range(loc_sim.stru.nt):
                loc_sim.stru.q.append(np.matmul(loc_sim.stru.u_avr[i], aux))
            
            loc_sim.stru.q = np.array(loc_sim.stru.q)
            dt_par = time.time()-t_2
        tiempos_full.append(dt_full)
        tiempos_prev.append(dt_par)
        print(dt_full, dt_par)
        
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    ax.plot(listado, tiempos_full, marker='o', label = 'Sin desc.')
    ax.plot(listado, tiempos_prev, marker='x', label = 'Con pre-prod')
    ax.set_xlabel('size')
    ax.set_ylabel('t [s]')
    ax.set_yscale('log')
    ax.legend()
    ax.grid()
    
    #Tests: 
        #Cantidad de modos = cantidad de DOFs
        #Tiempos de u_raw iguales?