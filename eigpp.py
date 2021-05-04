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

# Class offsets - Modificar para posibles parts. a leer

class offsets:

    c_info_eq = 1
    c_lumped_m = 1
    c_modes = 6

#------------------------------------------------------------------------------
# Definición de funciones

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
                    # print("Fin de tabla")
            except:
                stop_flag=False
            
    return(table_gen)


def rd_rawData(fname, nodos_interes):
    print("Raw data análisis completo")
    return *rd_eig(fname,nodos_interes),*rd_u(loc_sim.stru)
    # rd_eig(fname,nodos_interes)
    # rd_u(sim.stru)
    
    # return



def rd_eig(fname, nodos_interes):
    
    # leemos archivo así queda cargado y evitamos re-leer
    y = open(fname, 'r')
    x_dat = y.readlines()

    # Buscamos la linea que tiene la info de los nodos
    locs = search_string(x_dat, "Information Relative to the Equations Numbers")
    # print(x_dat[locs[0]+offsets.info_eq])
    raw_info_matrix = rd_SimpactTable(x_dat, locs[0]+offsets.c_info_eq)
    # print("paso 1: tabla completa de info eqs")
    # print(raw_info_matrix)
    # print("fin paso 1")
    # buscamos la matriz de masa
    locs_lumped = search_string(x_dat, "  Lumped Mass Matrix")
    raw_lumped_matrix = rd_SimpactTable(x_dat, locs_lumped[0]+offsets.c_lumped_m)
    # print("paso 2: matriz masa completa")
    # print(raw_lumped_matrix)
    # print("fin paso 2")
    
    # Ahora nos quedamos con lo que nos interesa
    # nodos de interés alimentados, sin últimas dos columnas
    # la matriz asociada tendrá len(nodos_interés) filas x 6 columnas
    
    n_rows = len(nodos_interes)
    # print(len(raw_info_matrix[:,0]))
    M_i_matrix = []
    for j in range(0, len(raw_info_matrix[:, 0])): #desplazamos fila a fila
        for a in range(0, n_rows): #desplazamos nodo a nodo (de interés)
            if raw_info_matrix[j, 0] == nodos_interes[a]:
                #guardamos la info que nos interesa
                M_i_matrix = np.append(M_i_matrix, raw_lumped_matrix[j])
                print("Matriz de masa leída")
                print(M_i_matrix)
                # print(' ')
                # M_i_matrix[a,:] = raw_lumped_matrix[j,1:]
    
    # por último, loopear para leer todos los modos
    locs_modes = search_string(x_dat, "linear dynamic eigen-mode analysis")

    # guardamos orden,frecuencias y períodos
    m_O_f_p = np.zeros((len(locs_modes), 3)) #sim.str.omsim.str.om
    Modes_reshaped_rows = np.zeros((len(locs_modes), len(nodos_interes)*6))

    for i in range(len(locs_modes)):
        # guardamos O,f,P
        m_O_f_p[i] = line_spliter(x_dat[locs_modes[i]])
        # guardamos y adecuamos tablas tablas
        raw_mode_table = rd_SimpactTable(x_dat, locs_modes[i]+offsets.c_modes)
        print("tabla modo full")
        print(raw_mode_table)
        print(" modos------------------------------")
        local_mode_row = []
        for j in range(0, len(raw_mode_table[:, 0])): #cantidad de filas de 1 tabla modo
            for a in range(0, n_rows): #cantidad de nodos de interes
                if raw_mode_table[j, 0] == nodos_interes[a]: #si el nodo es el que quiero
                    local_mode_row = np.append(local_mode_row, raw_mode_table[j,1:]) #apilar la columna correspondiente del modo
                    print("aporte local ---------------")
                    print(local_mode_row)
        Modes_reshaped_rows[i] = local_mode_row
        
    print("matriz modal phi finall --------------------------")
    print("ENDDDDDDDDD ---------")
    print(" ")
    # sim.stru.phi = np.transpose(Modes_reshaped_rows)
    y.close()
    # sim.stru.mass = M_i_matrix #matriz de masa adecuada
    # sim.stru.om = m_O_f_p[:,1] #frecuencias
    return M_i_matrix, np.transpose(Modes_reshaped_rows), m_O_f_p[:,1]

def rd_u(stru_clase,**kwargs):
    #espera input tipo sim.stru
    glob_u_raw = []
    glob_u_avr = []
    for a in range(0,len(stru_clase.nodes)):
        local_name = str(stru_clase.nodes[a])
        # callbat("pcolgante.@1", local_name, "0", "temp_file.txt")
        local_lines = open("temp_file.txt",'r')
        local_x_dat = local_lines.readlines()
        local_table_raw = rd_SimpactTable(local_x_dat,0)
        local_table_avr = np.copy(local_table_raw)
        # print(local_table_raw)
        if a == 0: #guardamos el tiempo si es
            glob_time = local_table_raw[:,0]
            total_ntime = len(glob_time) #cuantos instantes tengo

        if stru_clase.rdof: #si hay grados de libertad, rotar
            # print(local_table_avr[:,4:])
            local_table_avr[:,4:]= euler2axial(local_table_avr[:,4:]) #rotamo

        
        glob_u_raw.append(local_table_raw[:,1:]) #guardamos todo menos el t
        glob_u_avr.append(local_table_avr[:,1:]) #guardamos todo menos el t

    return glob_u_raw[0], glob_u_avr[0],glob_time,total_ntime


def euler2axial(cols):
    for i in range(len(cols[:,0])):
        # print(cols[i,:])
        R = Rot.from_euler('ZXZ',cols[i,:],degrees=False)
        cols[i,:] = R.as_rotvec()
    return(cols)
    
def callbat(arch_data, nodo, dof, nombre_otp):
    cmd = "(echo " + arch_data + " 1 0 && echo d && echo " + nodo + " " + dof + " " + nombre_otp + " && echo s) | curvas"  
    # "(echo pcolgante.@1 1 2 && echo d && echo 200000 0 d4.t && echo s) | curvas"
    # cmd = "(echo pcolgante.@1 1 2 && echo d && echo 200000 0 d4.t && echo s) | curvas"
    os.system(cmd)

def line_spliter(line):
    ret = []
    local_start = 0
    for i in range(0,len(line)):
        # local = []
        if line[i] == "-" or line[i] == " " or line[i] == "\n":
            if not i == 0:
                if not line[i-1] == "E":
                    # print(line[local_start:i])
                    a = line[local_start:i]
                    ret.append(a)
                    local_start = i
    # print(ret)
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
    print ('data saved')

"""
------------------------------------------------------------------------------
Zona para definir opciones y clases
------------------------------------------------------------------------------
"""

#Definición modelo propuesta por consigna

# sim.str.mass:   reales, ngl         - matriz de masa diagonal del modelo
# sim.str.om:     reales, nm          - vector de frecuencias naturales ordeandas de menor a mayor
# sim.str.phi:    reales, ngl x nm    - matriz modal (cada columna es un modo)
# sim.str.q:      reales, ngl x nt    - matriz con coordenadas modales a lo largo del tiempo, cada columna tiene el vector q de un instante
# sim.str.t:      reales, nt          - vector con instantes donde se conoce la solución (malla o grilla temporal)
# sim.str.nt:     entero, 1           - cantidad de instantes de tiempo
# sim.str.u_raw:  reales, ngl x nt    - matriz con desplazamientos a lo largo del tiempo, el vector u para un instante se guarda en una columna - como se obtienen de Simpact/curvas (si existen GL rotacionales, se expresan con ángulos de Euler y no se puede usar para hacer descomposición modal)
# sim.str.u_avr:  reales, ngl x nt    - como sim.str.u_raw pero transformando los GL rotacionales de ángulos de Euler a vector axial (avr: axial vector rotations)


class sim(): #Más fácil definir una clase que contenga todo lo de interés
    def __init__(self):
        self.name = "nombre" #Nombre de la simulación (se usa luego para guardar archivos y demás)
    #subclase para análisis estructural
    class stru:
        
        def __init__(self):
            self.nodes = [200001]  #etiquetas de nodos de interés - lista de enteros
            self.nnode = 1       #cantidad de nodos considerados
            self.ndof = 0        #cantidad de GL considerados (incluyendo los que son restringidos por condiciones de borde)
            self.rdOpt = 'bin'  #bandera para leer datos: 
                    #'raw' desde archivos provenientes de Simpact y delta1
                    #'bin' a partir de archivo binario guardado previamente
            self.rdof= True     #bandera que indique si hay GL rotacinales
            self.eigOpt = False   #bandera para hacer (o no) descomposición modal

        def resultados(clase,mass, phi, om, u_raw,u_avr,t,nt):
            clase.mass = mass
            clase.phi = phi
            clase.om = om
            clase.u_raw = u_raw
            clase.u_avr = u_avr
            clase.t = t
            clase.nt = nt
            clase.q = []

    # class aero:
    #Agregar para el caso aero

"""
------------------------------------------------------------------------------
Código Principal
------------------------------------------------------------------------------
"""

#Vemos si calcular y exportar o directamente cargar análisis previo
#generamos una clase para trabajar. Podríamos generar varias y analizar en serie.

loc_sim = sim()
loc_sim.stru = sim.stru()

#hacemos los cálculos, si corresponde:
    
if loc_sim.stru.rdOpt == 'raw':
    loc_sim.stru.resultados(*rd_rawData("test_leer_v1.txt", loc_sim.stru.nodes))
    #Hacemos las cuentas
    save_bin("hh", loc_sim.stru) #Exportamos la info

#Y si no, importamos
    
if loc_sim.stru.rdOpt=='bin':
    f = open("hh",'rb')
    loc_sim.stru = pickle.load(f)
    print("data loaded")
    f.close()   
        
#Ya tenemos la info cargada y lista. Ahora ver si corresponde hacer la descomposición.
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
aux = np.zeros((6,len(loc_sim.stru.phi[0])))
#Calculamos producto auxiliar
for i in range(len(loc_sim.stru.phi)):
    aux[i] = loc_sim.stru.phi[i] * loc_sim.stru.mass[i]
    
#Descomponemos tiempo a tiempo
for i in range(loc_sim.stru.nt):
    loc_sim.stru.q.append(np.matmul(loc_sim.stru.u_avr[i], aux))

loc_sim.stru.q = np.array(loc_sim.stru.q)
dt_par = time.time()-t_2