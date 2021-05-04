# -*- coding: utf-8 -*-
"""
Created on Sat Feb 27 10:21:04 2021

@author: Ramón Velazquez
"""

# Archivo de funciones - módulo

# imports necesarios p/ hacer las cuentas

import numpy as np
import re

# Clase offset


class offsets:

    # re encuentra numero de linea del ARCHIVO
    # en python es -1 dicho valor

    class info_eq:
        corr = 1
        est = [0, 8, 16, 24, 32, 40, 48, 56, 64, 72]

    class lumped_m:
        corr = 1

    class modes:
        corr = 6
# Main


def search_string(x_dat, frase):
    lin_count = 0  # contador de lineas
    locs = []  # ubicaciones dónde encuentra (line number python = -1 realidad)
    for txt in x_dat:
        lin_count = lin_count + 1
        x = re.search(frase, txt)
        if x != None:
            locs.append(lin_count)
            print(lin_count)
    # devolvemos la lista con las lineas de x_dat que tienen la coincidencia
    return(locs)


# función que devuelve una tabla asociada (pensada p modos norm)
def rd_SimpactTable(x_dat, start_line, **kwargs):

    # importante
    # Directamente trabajamos con start_line dándo este valor como offset al llamarla

    modo_lect = 1  # idea para automatizar acá mismo filtro de tablas con cosas que no nos sirvan
    est_flag = False

    if 'modo' in kwargs:
        modo_lect = kwargs.get('modo')

    if 'est_particular' in kwargs:
        est_flag = True
        est_particular = kwargs.get('est_particular')

    stop_flag = True
    # 2 en eigen modes
    # distinto en demás? Testear
    counter = 0
    while stop_flag:
        loc_arr = []

        if est_flag:  # debemos usar estructura lugar a lugar
            for a in range(0, len(est_particular)-1):
                try:
                    loc_arr.append(
                        float(x_dat[start_line+counter][est_particular[a]:est_particular[a+1]]))
                except:
                    loc_arr.append(x_dat[a])
            if counter == 0:  # primer ciclo
                table_gen = np.array([loc_arr])
            else:
                table_gen = np.append(table_gen, [loc_arr], axis=0)
            counter = counter+1

        else:  # lectura normal spliteable
            for a in x_dat[start_line+counter].split():
                try:
                    # Necesitamos un try por si no hay números (ejemplo lumped M)
                    loc_arr.append(float(a))
                except:
                    loc_arr.append(a)
            if counter == 0:  # primer ciclo
                table_gen = np.array([loc_arr])
            else:
                table_gen = np.append(table_gen, [loc_arr], axis=0)
            counter = counter+1

        print(loc_arr)
        if x_dat[start_line + counter] == '\n':
            stop_flag = False
            print("Fin de tabla")
    return(table_gen)


def rd_rawData(nodos_interes):
    rd_eig("pcolgante.rsn", nodos_interes)


def rd_eig(fname, nodos_interes):
    # leemos archivo así queda cargado y evitamos re-leer
    y = open(fname, 'r')
    x_dat = y.readlines()

    # Buscamos la linea que tiene la info de los nodos
    locs = search_string(
        x_dat, "Information Relative to the Equations Numbers")
    # print(x_dat[locs[0]+offsets.info_eq])
    raw_info_matrix = rd_SimpactTable(
        x_dat, locs[0]+offsets.info_eq.corr, est_particular=offsets.info_eq.est)

    # buscamos la matriz de masa
    locs_lumped = search_string(x_dat, "  Lumped Mass Matrix")
    raw_lumped_matrix = rd_SimpactTable(
        x_dat, locs_lumped[0]+offsets.lumped_m.corr)

    # Ahora nos quedamos con lo que nos interesa
    # nodos de interés alimentados, sin últimas dos columnas
    # la matriz asociada tendrá len(nodos_interés) filas x 6 columnas
    n_rows = len(nodos_interes)

    M_i_matrix = []
    for j in range(0, len(raw_info_matrix[:, 0])):
        for a in range(0, n_rows):
            if raw_info_matrix[j, 0] == nodos_interes[a]:
                M_i_matrix = np.append(M_i_matrix, raw_lumped_matrix[j, 1:])
                # print(M_i_matrix)
                # print(' ')
                # M_i_matrix[a,:] = raw_lumped_matrix[j,1:]

    # por último, loopear para leer todos los modos
    locs_modes = search_string(x_dat, "linear dynamic eigen-mode analysis")

    # guardamos orden,frecuencias y períodos
    m_O_f_p = np.zeros((len(locs_modes), 3))
    Modes_reshaped_rows = np.zeros((len(locs_modes), len(nodos_interes)*6))

    for i in range(len(locs_modes)):
        # guardamos O,f,P

        # guardamos y adecuamos tablas tablas
        raw_mode_table = rd_SimpactTable(
            x_dat, locs_modes[i]+offsets.modes.corr)
        local_mode_row = []
        for j in range(0, len(raw_mode_table[:, 0])):
            for a in range(0, n_rows):
                if raw_mode_table[j, 0] == nodos_interes[a]:
                    local_mode_row = np.append(
                        local_mode_row, raw_mode_table[j, 1:])
                    # print(local_mode_row)
        Modes_reshaped_rows[i] = local_mode_row
    print("modeeee")
    print(Modes_reshaped_rows)
    mass = 'asignar multiples salidas a sim.str.mass y phi'
    return(M_i_matrix)

def rd_u(nodos_interes, file, folder,**kwargs):
    header_tipo = [0,1,2] #relativo posición y cantidad de lineas a saltear segun header
    #llamar generador de .bat y crear archivo
    #consultar dudas
    
    #for i in (nodos_interes):
    #conservar nombres
    loc_fname = "test_leer_v1.txt" #nodos_interes[i]
    op_file = open(loc_fname,'r')
    dat = op_file.readlines()

    fname = "uz_200001.dat"
    # fname = "test_leer_v1.txt"
    
    y = open(fname,'r')
    x_dat = y.readlines()
    
    #ver si tipo 0 o tipo tecplot
    if 'tipo' in kwargs:
        tipo = kwargs.get('tipo')
    else:
        tipo = 0 #por default
    
    n_cols = len(dat[header_tipo[tipo]].split())

    data_arr = []
    
    for a in dat[header_tipo[tipo]:]:
        valores = [0,0]
        for j in range(0,n_cols):
            valores[j] = float(a.split()[j])
        data_arr.append(valores)
    data_arr = np.array(data_arr)

    return(data_arr)

# Test
nodes = [200002, 200003]
gg = rd_rawData(nodes)
# y = open("pcolgante.rsn",'r')
# x_dat = y.readlines()
