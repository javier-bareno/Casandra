#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 21 23:21:46 2019

@author: MIMAT_JB

Quick test of pusputt.py
"""

from pysputt import *
import matplotlib.pyplot as plt
import numpy as np

input_data = CalcInput()

input_data.Gas_symbol = 'Ar'
#retrieve gas atom properties
input_data.MG = 39.948 #atomic mass
input_data.ZG = 18 #atomic number
input_data.RG = 0.98 # atomic radius
input_data.UG = 5 #
input_data.QG = 1
input_data.DENG = 5
input_data.KAPA = 0.0002

input_data.Target_symbol= 'Si'
input_data.MA = 28 #atomic mass
input_data.ZA = 14
input_data.RA = 2
input_data.U0 = 4.63
input_data.QZ = 0.75
input_data.DEN = 2.33

input_data.disch_volt = 100.
input_data.disch_cur = 1.1
input_data.pres = 1.2
input_data.Eavg = 10.
input_data.T0 = 300.

input_data.Tcorrect = True

input_data.in_rad = 20
input_data.out_rad = 30
input_data.dist = 150
input_data.L5 = 2
input_data.N1 = 10
input_data.N2 = 20
input_data.grid_R2 = 2
        
calc_sput = Sputt_plan_circ(input_data)
#plt.plot(calc_sput.pR, calc_sput.DEP1, '.', calc_sput.pR, calc_sput.DEP22, '.')

# load cpp data
fname = '/Users/MIMAT_JB/Desktop/Github/pyCasandra/Casw2019/Test_03_noTC_out copy.txt'
with open(fname,'r')as ifile:
    dt = ifile.readlines()

l0 =0
for i in range(len(dt)):
    if dt[i][0:4] == "R(mm":
        l0 = i
        break

cpp_R, cpp_total, cpp_balistic, cpp_diff, cpp_resp = [],[],[],[],[]
for l in dt[l0+1::]:
    l=l.split()
    cpp_R.append(float(l[0]))
    cpp_total.append(float(l[1]))
    cpp_balistic.append(float(l[2]))
    cpp_diff.append(float(l[3]))
    cpp_resp.append(float(l[4]))

cpp_R = np.array(cpp_R)
cpp_total = np.array(cpp_total)
cpp_balistic = np.array(cpp_balistic)
cpp_diff = np.array(cpp_diff)
cpp_resp = np.array(cpp_resp)

R= calc_sput.pR
total = calc_sput.pTot
balistic = calc_sput.DEP1
diff = calc_sput.DEP21
resp = calc_sput.DEP22

d_tot = cpp_total-total
d_bal = cpp_balistic - balistic  
d_dif = cpp_diff -diff
d_res = cpp_resp -resp

dr_tot = cpp_total / total
dr_bal = cpp_balistic / balistic  
dr_dif = cpp_diff  /diff
dr_res = cpp_resp /resp
#print(dr_res)

fname = '/Users/MIMAT_JB/Desktop/Github/pyCasandra/Casw2019/'
fname += 'compare.csv'
with open(fname, 'w') as ofile:
    ofile.write('Comparisson of new and old Casandra\n\n')
    ofile.write('R(mm)\t')
    for l in ['Total', 'Ballistic', 'Diffusive', 'Resputtered']:
        for ll in ['(old)', '(new)', '(ratio)']:
            ofile.write(l + ll + '\t')
    ofile.write('\n')
    for i in range(len(R)):
        ofile.write(str( R[i]) + '\t')
        ofile.write(str( cpp_total[i]) + '\t')
        ofile.write(str( total[i] ) + '\t')
        ofile.write(str( dr_tot[i] ) + '\t')
        ofile.write(str( cpp_balistic[i] ) + '\t')
        ofile.write(str( balistic[i] ) + '\t')
        ofile.write(str( dr_bal[i] ) + '\t')
        ofile.write(str( cpp_diff[i] )+ '\t')
        ofile.write(str( diff[i] ) + '\t')
        ofile.write(str( dr_dif[i] ) + '\t')
        ofile.write(str( cpp_resp[i] ) + '\t')
        ofile.write(str( resp[i] ) + '\t')
        ofile.write(str( dr_res[i] ) + '\n')

maxdr=np.max([dr_tot, dr_bal, dr_dif, dr_res ])
mindr=np.min([dr_tot, dr_bal, dr_dif, dr_res ])
print(maxdr-1, 1-mindr)
        
    
            
                
                

