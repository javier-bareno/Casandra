#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 20 10:13:39 2019

@author: MIMAT_JB
"""

import numpy as np

def arr_from_AsCII( fname):
    ''' Dump a numpy array in ordered fashion to ASCII file'''
    with open(fname, 'r') as ifile:
        lines = ifile.readlines()
        
    # first line has shape of array:
    sh=[]
    print(lines[0].split())
    for Ni in lines[0].split():
        sh.append(int(Ni))
    sh=tuple(sh)
    ndims = len(sh)
    
    ret_arr = np.zeros(sh)
        
    #each following line starts with indices, followed by content
    for line in lines[1:]:
        sh=[]
        lsp = line.split();
        for Ni in lsp[0:ndims]:
            sh.append(int(Ni))
        sh=tuple(sh)
        ret_arr[sh] = float(lsp[-1]) 
            
    return(ret_arr)
    
if __name__ == "__main__":
    ''' Test '''
    x=np.arange(16).reshape(4,4)
    print('2D array')
    print(x)
    fname = 'test_arr.dat'
    with open(fname, 'w') as ofile:
        ofile.write('4\t4\n')
        for i in range(4):
            for j in range(4):
                outstr = str(i) + '\t' + str(j) + '\t' + str(x[i,j]) + '\n' 
                ofile.write(outstr)
    print('Load array')
    x2 = arr_from_AsCII( fname)
    print(x2)
    
    x=np.arange(16).reshape(2,4,2)
    print('3D array')
    print(x)
    fname = 'test_arr.dat'
    with open(fname, 'w') as ofile:
        ofile.write('2\t4\t2\n')
        for i in range(2):
            for j in range(4):
                for k in range(2):
                    outstr = str(i) + '\t' + str(j) + '\t' 
                    outstr += str(k) + '\t' +  str(x[i,j,k]) + '\n' 
                    ofile.write(outstr)
    print('Load array')
    x2 = arr_from_AsCII( fname)
    print(x2)
    