# -*- coding: utf-8 -*-
"""
Created on Mon Mar  1 18:30:41 2021

@author: Ramon Velazquez
"""

#Test fácil para escribir .bat ejecutable y ejecutarlo

import os
from datetime import datetime
import subprocess

# datetime object containing current date and time

# now = datetime.now()
 
# print("now =", now)

# # dd/mm/YY H:M:S
# dt_string = now.strftime("%d/%m/%Y %H:%M:%S")
# print("date and time =", dt_string)	

# try:
#     os.makedirs('my_folder')

# except:
#     print("gg")


myBat = open(r'my_folder\ejecutame.bat','w+')

#Estructura comandos

org_name = "pcolgante.@1"
num_1 = 1
num_2 = 0
disp = "d"
node = str(200005)
dofs = str(0)
filename = "my_folder\outp"


start = "@ECHO off"
cont_1 = "("
cont_2 = "echo " + str(org_name) + " " + str(num_1) + " " + str(num_2) 
cont_3 = "echo " + str(disp)
cont_4 = "echo " + node + " " + dofs + " " + filename
cont_5 = "echo s"
cont_6 = ") | curvas"

comandos = [start,cont_1,cont_2, cont_3, cont_4, cont_5, cont_6]
for a in comandos:
    myBat.write(a)
    myBat.write('\n')
myBat.close()

subprocess.call(['my_folder\\ejecutame.bat'])
