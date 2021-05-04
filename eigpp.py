# -*- coding: utf-8 -*-
"""
Created on Thu Feb 25 17:52:06 2021

@author: Ramón Velazquez 



"""

#Zona para importar libs y módulos

import numpy as np
import mods_v1 as md


#Zona para definir opciones

class sim: #Más fácil definir una clase que contenga todo lo de interés
    
    name = "nombre simulación" #Nombre de la simulación (se usa luego para guardar archivos y demás)
    #subclase para estr
    class struct:
        nodes = [200002,200003]  #etiquetas de nodos de interés - lista de enteros
        nnode = 1       #cantidad de nodos considerados
        ndof = 0        #cantidad de GL considerados (incluyendo los que son restringidos por condiciones de borde)
        rdOpt = 'raw'  #bandera para leer datos: 
                #'raw' desde archivos provenientes de Simpact y delta1
                #'bin' a partir de archivo binario guardado previamente
        rdof= True      #bandera que indique si hay GL rotacinales
        eigOpt = True   #bandera para hacer (o no) descomposición modal
        # def __init__(self, name):
        #     self.name = name
        ## idea para automatizar nombres de clases / archivos a analizar
    
    # class aero:
    #Agregar para el caso aero

#Main

if sim.struct.rdOpt == 'raw':
    md.rd_rawData(sim.struct.nodes)