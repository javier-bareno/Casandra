#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 19 13:24:52 2019

@author: MIMAT_JB
"""

import numpy as np
import time
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
        #calc DEP21
        DEP21 = np.zeros(nr+1)
        for iir in range(nr +1):
            r = grid_R2*iir
            sr = 0. #accumulator on r substrate
            ld4 =[]
            for iz in range(nz):
                zc = (iz+0.5)*dz
                sz=0. # accumulator on z
                ld3=[]
                for ir in range(nr):
                    rc = (ir+0.5)*drv
                    ld2=[]
                    for k2 in range(N2):
                        fi = (k2+0.5)*dfi
                        rop = np.sqrt(r*r +rc*rc - 2*r*rc*np.cos(fi))
                        A = (Z-zc)/(2*Z)
                        B = (Z+zc)/(2*Z)
                        C= rop /(2*Z)
                        kk1 = A/np.power(A*A+C*C, 1.5)
                        kk2 = B/np.power(B*B+C*C, 1.5)
                        s0 = kk1 - kk2
                        s1 = (1+A) /np.power((1+A)*(1+A)+C*C, 1.5) 
                        s1 += (1-B) /np.power((1-B)*(1-B)+C*C, 1.5)
                        s1 -= (1-A) /np.power((1-A)*(1-A)+C*C, 1.5) 
                        s1 -= (1+B) /np.power((1+B)*(1+B)+C*C, 1.5)
                        s0 += s1
                        s2= J[ir][iz]*drv*dz*dfi*rc*s0/(16*Z*Z*np.pi)
                        ld1=[rop, A, B, C, kk1, kk2, s0, s1, s2]
                        ld2.append(ld1)
                        sz+=s2
                    ld3.append(ld2)
                ld4.append(ld3)
                sr = sr + sz
            log_dump.append(ld4)
            DEP21[iir] = 2*sr
                        
        '''
        # now DEP22[]
        DEP22 = np.zeros(nr+1)
        for ir in range(nr+1):
            r = grid_R2*ir
            sr = 0 #accumulator on r substrate
            for iz in range(nz):
                zc = (iz+0.5)*dz
                sz=0 # accumulator on z
                for ir in range(nr):
                    rc = (ir+0.5)*drv
                    for k2 in range(N2):
                        fi = (k2+0.5)*dfi
                        rop = np.sqrt(r*r +rc*rc - 2*r*rc*np.cos(fi))
                        A= (0.0-zc)/(2*Z)
                        B= (0.0+zc)/(2*Z)
                        C= rop /(2*Z)
                        s0= A/np.power(A*A+C*C, 1.5) 
                        s0 -= B/np.power(B*B+C*C, 1.5)
                        s1= (1+A) /np.power((1+A)*(1+A)+C*C, 1.5) 
                        s1 += (1-B) /np.power((1-B)*(1-B)+C*C, 1.5)
                        s1 -= (1-A) /np.power((1-A)*(1-A)+C*C, 1.5)  
                        s1 -= (1+B) /np.power((1+B)*(1+B)+C*C, 1.5)
                        s0 += s1
                        s2= J[ir, iz]*drv*dz*dfi*rc*s0/(16*Z*Z*np.pi)
                        sz+=s2
                sr+=sz
            DEP22[ir]=2*sr
            
        return(DEP21, DEP22)
        '''
        return(np.array(log_dump), DEP21)

if __name__ == "__main__":
    tst = test(Id())
    r, DEP21=tst.for_difusion_comp()
    
    fname = '/Users/MIMAT_JB/Desktop/Github/VSc++_test/files/cpp_DEP21_calc.txt'
    cppcalc = arr_from_AsCII(fname)
    print('Compare python and cpp intermediate results:\n')
    print('Var\tcpp\tpython\tratio\n')
    print('rop', r[0,0,0,0,0], cppcalc[0,0,0,0,0] )
    print('A', r[0,0,0,0,1], cppcalc[0,0,0,0,1])
    print('B', r[0,0,0,0,2], cppcalc[0,0,0,0,2])
    print('C', r[0,0,0,0,3], cppcalc[0,0,0,0,3])
    print('kk1', r[0,0,0,0,4], cppcalc[0,0,0,0,4])
    print('kk2', r[0,0,0,0,5], cppcalc[0,0,0,0,5])
    print('s0', r[0,0,0,0,6], cppcalc[0,0,0,0,6])
    print('s1', r[0,0,0,0,7], cppcalc[0,0,0,0,7])
    print('s2', r[0,0,0,0,8], cppcalc[0,0,0,0,8])