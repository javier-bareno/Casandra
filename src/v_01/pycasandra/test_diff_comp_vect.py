#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 19 13:24:52 2019

@author: MIMAT_JB
"""

import numpy as np
import timeit
from test_diff_source import test as getJ
from py2c import arr_from_AsCII

class Id():
    def __init__(self, L5=2, R1=20, R2=30, Z=150,
                 LAMBDA=103.194, N1=10, N2=20 ):
        self.L5 = L5
        self.R1 = R1
        self.R2 = R2
        self.grid_R2 = 2
        self.dist = Z
        self.mfp_lambda =[LAMBDA]
        self.N1 = N1
        self.N2 = N2
        

class test():
    def __init__(self, inputData):
        self.inputData = inputData
        self.J = getJ(inputData).vec_difusion_source()

    def for_difusion_comp(self):
        
        L5 = self.inputData.L5
        R1 = self.inputData.R1
        R2 = self.inputData.R2
        grid_R2 = self.inputData.grid_R2
        Z = self.inputData.dist
        LAMBDA = self.inputData.mfp_lambda[-1]
        N1 = self.inputData.N1
        N2 = self.inputData.N2
        J = self.J
        
        nr = int(np.ceil(L5*R2/grid_R2)) #num points in R axis of substrate
        drv = L5*R2/nr #= grid_R2, step between points
        nz = int(np.ceil(Z/2)) #number of integration steps in distance between
        #target and substrate. Step is assumed 2 mm
        # ensure dz < lambda
        dz = Z/nz
        if (dz > LAMBDA):
            nz = int(np.ceil(Z/LAMBDA))
            dz = Z/nz
        dfi = np.pi/N2
            
        log_dump =[]
        
        #code of interest
        
        #DEP21 = np.zeros(nr+1)
        r = grid_R2 * np.arange(nr+1) # r[iir]
        zc = dz * (np.arange(nz)+ 0.5) #zc[{iz}]
        rc = (np.arange(nr)+0.5)*drv # rc[{ir}]
        fi = (np.arange(N2)+0.5)*dfi #fi[{k2}]
        
        rop0 = rc.reshape(1,1,nr,1) * np.cos(fi).reshape((1,1,1,N2))
        rop1 = r.reshape(nr+1,1,1,1) * rop0
        rop0 = (r**2).reshape(nr+1, 1,1,1) + (rc**2).reshape(1,1,nr,1) -2 *rop1
        #rop = np.sqrt(rop0)
        #C= rop /(2*Z)
        C2 = rop0 /(4*Z**2) # sh = (nr+1, 1, nr, N2)
        #rop0 = rop1 = False free some memory
                        
        A = (Z-zc.reshape((1,nz,1,1)))/(2*Z)
        B = (Z+zc.reshape((1,nz,1,1)))/(2*Z)
        
        kk1 = A/np.power(A**2 + C2, 1.5) # sh = (nr+1, nz, nr, N2)
        kk2 = B/np.power(B**2 + C2, 1.5)
        s0 = kk1 - kk2
        s1 = (1+A) /np.power((1+A)**2 + C2, 1.5) 
        s1 += (1-B) /np.power((1-B)**2 +C2, 1.5)
        s1 -= (1-A) /np.power((1-A)**2 +C2, 1.5) 
        s1 -= (1+B) /np.power((1+B)**2 +C2, 1.5)
        s0 += s1
        s2 = J.T.reshape(1, nz, 1,nr) * rc.reshape(1,1,1,nr) 
        sz = np.matmul(s2, s0)
        DEP21 = 2* drv *dz *dfi / (16 * np.pi * Z**2 ) *np.sum(sz, axis=(-1,-2,-3))
        
          # now DEP22[]
        #DEP22 = np.zeros(nr+1)
        #for ir in range(nr+1):
         #   r = grid_R2*ir
          #  sr = 0 #accumulator on r substrate
           # for iz in range(nz):
            #    zc = (iz+0.5)*dz
             #   sz=0 # accumulator on z
              #  for ir in range(nr):
               #     rc = (ir+0.5)*drv
                #    for k2 in range(N2):
                 #       fi = (k2+0.5)*dfi
                       # rop = np.sqrt(r*r +rc*rc - 2*r*rc*np.cos(fi))
                    
        A2= (0.0-zc.reshape((1,nz,1,1)))/(2*Z)
        B2= (0.0+zc.reshape((1,nz,1,1)))/(2*Z)
        s0= A2/np.power(A2**2 + C2, 1.5) 
        s0 -= B2/np.power(B2**2 +C2, 1.5)
        s1= (1+A2) /np.power((1+A2)**2 + C2, 1.5) 
        s1 += (1-B2) /np.power((1-B2)**2 +C2, 1.5)
        s1 -= (1-A2) /np.power((1-A2)**2 + C2, 1.5)  
        s1 -= (1+B2) /np.power((1+B)**2 + C2, 1.5)
        s0 += s1
        s2 = J.T.reshape(1, nz, 1,nr) * rc.reshape(1,1,1,nr) 
        sz = np.matmul(s2, s0)
        DEP22 = 2* drv *dz *dfi / (16 * np.pi * Z**2 ) *np.sum(sz, axis=(-1,-2,-3))
        
        return(DEP21, DEP22)
        

if __name__ == "__main__":
    tst = test(Id())
    #r, DEP21=tst.for_difusion_comp()
    
    fname = '/Users/MIMAT_JB/Desktop/Github/VSc++_test/files/cpp_DEP21_calc.txt'
    #cppcalc = arr_from_AsCII(fname)
    
    DEP21, DEP22 = tst.for_difusion_comp()
    print('timeit:', timeit.timeit(tst.for_difusion_comp, number=10))