#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 20 15:55:00 2019

@author: MIMAT_JB

Goal is to test translation to python of Casandra
In this test, compare J calc from python vectorized version 
to cpp version in file
"""

import numpy as np
import pysputt as ps
import py2c

# set up calculation with same parameters as cpp
inputdata = ps.CalcInput()

inputdata.Gas_symbol = 'Ar'
inputdata.MG = 39.948
inputdata.ZG = 18
inputdata.RG = 0.98
inputdata.UG = 5
inputdata.QG = 1
inputdata.DENG = 5
inputdata.KAPA = 0.0002

inputdata.Target_symbol = 'Si'
inputdata.MA = 28
inputdata.ZA = 14
inputdata.RA = 2
inputdata.UO = 4.63
inputdata.QZ = 0.75
inputdata.DENA = 2.33

inputdata.disch_volt = 100
inputdata.Eavg = 10
inputdata.disch_cur = 1.1
inputdata.pres = 1.2
inputdata.T0 = 300

inputdata.TC = 0
inputdata.Tcorrect = False
inputdata.in_rad =20
inputdata.out_rad =30
inputdata.N1 = 10
inputdata.N2 = 20
inputdata.L5 = 2
inputdata.grid_R2 = 2.
inputdata.dist = 150

sput_calc = ps.Sputt_plan_circ(inputdata)
fname = "/Users/MIMAT_JB/Desktop/Github/VSc++_test/files/cpp_J.txt"
J = py2c.arr_from_AsCII(fname)
r =J /sput_calc.J
print('Comparing calculated J:\n')
print('Min J_cpp / Jpy :', np.min(r))
print('Max J_cpp / Jpy :', np.max(r))
    