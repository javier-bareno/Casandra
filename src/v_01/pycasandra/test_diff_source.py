#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 18 18:56:21 2019

@author: MIMAT_JB
"""

import numpy as np
import time

class Id():
    def __init__(self, L5, R1, R2, Z, LAMBDA, N1, N2 ):
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

    def for_difusion_source(self):
        '''Calculate currfent J(r,z)'''
        L5 = self.inputData.L5
        R1 = self.inputData.R1
        R2 = self.inputData.R2
        grid_R2 = self.inputData.grid_R2
        Z = self.inputData.dist
        LAMBDA = self.inputData.mfp_lambda[-1]
        N1 = self.inputData.N1
        N2 = self.inputData.N2
        
        
        nr = int(np.ceil(L5*R2/grid_R2)) #num points in R axis of substrate
        drv = L5*R2/nr #= grid_R2, step between points
        
        nz = int(np.ceil(Z/2)) #number of integration steps in distance between
        #target and substrate. Step is assumed 2 mm
        # ensure dz < lambda
        dz = Z/nz
        if (dz > LAMBDA):
            nz = int(np.ceil(Z/LAMBDA))
            dz = Z/nz
        
        J=np.zeros((nr, nz))
        for ir in range(nr):
            r = (0.5+ir)*drv
            for iz in range(nz):
                z = (iz+0.5) *dz
                s=0
                for k1 in range(N1):
                    rc = R1 + (k1+0.5) *(R2-R1)/N1
                    for k2 in range(N2):
                        fi = (k2 +0.5) * np.pi/N2
                        r0 = np.sqrt(
                                z*z + r*r + rc*rc  - 2*r*rc*np.cos(fi))
                        s += rc* z * np.exp(-r0/LAMBDA) /(
                                LAMBDA * np.power(r0, 3))
                J[ir, iz] = ((R2-R1)/N1) * (np.pi/N2) * s * 2 /np.pi
        return(J)
        
    def vec_difusion_source(self):
        '''Calculate currfent J(r,z)'''
        L5 = self.inputData.L5
        R1 = self.inputData.R1
        R2 = self.inputData.R2
        grid_R2 = self.inputData.grid_R2
        Z = self.inputData.dist
        LAMBDA = self.inputData.mfp_lambda[-1]
        N1 = self.inputData.N1
        N2 = self.inputData.N2
        
        nr = int(np.ceil(L5*R2/grid_R2)) #num points in R axis of substrate
        drv = L5*R2/nr #= grid_R2, step between points
        
        nz = int(np.ceil(Z/2)) #number of integration steps in distance between
        #target and substrate. Step is assumed 2 mm
        # ensure dz < lambda
        dz = Z/nz
        if (dz > LAMBDA):
            nz = int(np.ceil(Z/LAMBDA))
            dz = Z/nz
        
        #J=np.zeros(nr, nz)
        
        Mr =  drv * (np.arange(nr)+0.5)
        Mz = dz * (np.arange(nz)+0.5)
        Mrc =  R1 + ((R2-R1)/N1) * (np.arange(N1)+0.5) 
        Mfi = (np.pi/N2) * (np.arange(N2)+0.5)
        
        Mr00 = Mr.reshape((nr,1,1,1)) * Mrc.reshape((1,1,N1,1))
        Mr01 = Mr00 * np.cos(Mfi).reshape((1,1,1,N2))
        Mr00 = (Mz**2).reshape(1,nz,1,1) - 2 * Mr01
        Mr01 = Mr00 + (Mr**2).reshape((nr,1,1,1)) 
        Mr00 = Mr01 + (Mrc**2).reshape((1,1,N1,1))
        Mr0 = np.sqrt(Mr00)
        
        S0 = np.exp(-Mr0/LAMBDA
                   ) / (LAMBDA * np.power(Mr0, 3))
        S1 = np.sum(S0, -1)
        S0 = np.matmul(S1, Mrc)
        S1 = S0 * Mz.reshape(1,nz)
        
        J= (((R2-R1)/N1) * (np.pi/N2) * 2 /np.pi) * S1
        return(J)
        
                
'''       #for ir in range(nr):
         #   r = (0.5+ir)*drv
            #for iz in range(nz):
                #z = (iz+0.5) *dz
                s=0
                
                
                
                #for k1 in range(N1):
                    #rc = R1 + (k1+0.5) *(R2-R1)/N1
                    #for k2 in range(N2):
                        #fi = (k2 +0.5) * np.Pi/N2
                        r0 = np.sqrt(
                                z*z + r*r + rc*rc  - 2*r*rc*np.cos(fi))
                        s += rc* z * np.exp(-r0/LAMBDA) /(
                                LAMBDA * np.pow(r0, 3))
                J[ir, iz] = ((R2-R1)/N1) * (np.pi/N2) * s * 2 /np.pi'''
        
                        
if __name__ == "__main__":
    print('Lets go')
    tst = test(Id(2, 10, 20, 100, 5, 10, 12))
    t0 = time.time()
    tv= tst.vec_difusion_source()
    t1=time.time()
    tf = tst.for_difusion_source()
    t2 = time.time()
    print(np.max(np.abs(tv/tf)))
    print(t1-t0, t2-t0)
    print(np.max(tf/tv))
    print(np.min(tf/tv))