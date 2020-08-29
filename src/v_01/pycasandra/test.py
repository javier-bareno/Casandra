#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 17 11:54:11 2019

@author: MIMAT_JBAr
"""

import pickle
from atomCOE import *
import pysputt

def test1():
    f_coe_name = '/Users/MIMAT_JB/Desktop/Github/pyCasandra/COE files/atoms.coe'
    with open(f_coe_name, 'rb') as ifile:
        atCOE = pickle.load(ifile)
        
    input_data = pysputt.CalcInput()
    input_data.iter_input(atCOE)
    with open('test_input_calc.dat', 'w') as ofile:
        input_data.toASCII(ofile)

def test2():
    input_data = pysputt.CalcInput()
    with open('test_input_calc.dat', 'r') as ifile:
        input_data.fromASCII(ifile)
    print(input_data)
    
if __name__ == "__main__":
    test2()