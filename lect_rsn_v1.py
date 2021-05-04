# -*- coding: utf-8 -*-
"""
Created on Wed Feb 24 20:47:46 2021

@author: Rodrigo Velazquez
"""

#pruebas tentativas de análisis del archivo pcolgante.rsn

#importar
import re
import numpy as np

#definir
fname = "pcolgante.rsn"

#trabajar
y = open(fname,'r')
x_dat = y.readlines()

def search_string(raw_data, frase):
    lin_count = 0 #contador de lineas
    locs = [] #ubicaciones dónde encuentra (line number python = -1 realidad)
    for txt in raw_data:
        lin_count = lin_count + 1
        x = re.search(frase, txt)
        if x != None:
            locs.append(lin_count)
            res_1 = x.span()
            line = lin_count
            # print(line)
    return(locs) #devolvemos la lista con las lineas de x_dat que tienen la coincidencia

def read_table(start_loc,**kwargs): #función que devuelve una tabla asociada (pensada p modos norm)
    modo_lect = 1 #idea para automatizar acá mismo filtro de tablas con cosas que no nos sirvan 
    if 'modo' in kwargs:
        modo_lect = kwargs.get('modo')
        
    stop_flag = True
    start_line = start_loc + 1 #corrimiento de lineas para arrancar tabla
    #2 en eigen modes
    # distinto en demás? Testear
    counter = 0
    while stop_flag:
        loc_arr = []
        for a in x_dat[start_line+counter].split():
          try:  
              loc_arr.append(float(a)) #Necesitamos un try por si no hay números (full lumped)
          except:
              loc_arr.append(a)
        if counter == 0: #primer ciclo
            table_gen = np.array([loc_arr])
        else:
            table_gen = np.append(table_gen,[loc_arr],axis=0)
        counter = counter+1
        print(loc_arr)
        if x_dat[start_line + counter] == '\n' and x_dat[start_line+counter+1]=='\n':
            stop_flag = False
            print("Fin de tabla")
    return(table_gen)
        
