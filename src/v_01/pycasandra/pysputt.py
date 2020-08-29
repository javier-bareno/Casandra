#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 17 10:43:20 2019

@author: Javier Bareno and Ivan Petrov

This is a port of Ivan's sputtering code (Casandra in cpp version)

This module defines basic functionality. Other modules will build on
it to make user programs with GUI etc.
"""
import numpy as np

class CalcInput():
    ''' Class to handle input data for calculation. Saving and loading
    input data files, user input, etc.
    '''
    
    def __init__(self):
        pass
        
    def iter_input(self, atomCOE):
        ''' INTERACTIVE INPUT OF DATA
        Also used to doc parameters
        '''
        
        self.atCOE =atomCOE 
         #object with table of atomic coeficients 
        
        print('DEPOSITION RATE CALCULATION FOR CIRCULAR, PLANAR MAGNETRON')
        print('Let\'s talk about the plasma:')
        self.Gas_symbol = input('ENTER CHEMICAL SYMBOL OF SPUTTERING GAS: ')
        #retrieve gas atom properties
        self.MG = self.atCOE.get_G(self.Gas_symbol, 'M') #atomic mass
        self.ZG = self.atCOE.get_G(self.Gas_symbol, 'Z') #atomic number
        self.RG = self.atCOE.get_G(self.Gas_symbol, 'R') # atomic radius
        self.UG = self.atCOE.get_G(self.Gas_symbol, 'U') #
        self.QG = self.atCOE.get_G(self.Gas_symbol, 'Q')
        self.DENG = self.atCOE.get_G(self.Gas_symbol, 'DEN') 
        self.KAPA = self.atCOE.get_G(self.Gas_symbol, 'KAPA')
        
        self.Target_symbol=input('ENTER CHEMICAL SYMBOL OF TARGET ELEMENT: ')
        self.MA = self.atCOE.get_A(self.Target_symbol, 'M') #atomic mass
        self.ZA = self.atCOE.get_A(self.Target_symbol, 'Z')
        self.RA = self.atCOE.get_A(self.Target_symbol, 'R')
        self.U0 = self.atCOE.get_A(self.Target_symbol, 'U')
        self.QZ = self.atCOE.get_A(self.Target_symbol, 'Q')
        self.DEN = self.atCOE.get_A(self.Target_symbol, 'DEN') #superseeds gas
        #self.RA = self.robinson_R(self.ZA)
        
        self.disch_volt = float(input('ENTER DISCHARGE VOLTAGE  IN VOLTS: '))
        self.disch_cur = float(input('ENTER DISCHARGE CURRENT IN AMPS: '))
        self.pres = float(input('ENTER GAS PRESSURE IN PASCALS: '))
        self.Eavg = float(input('ENTER AVERAGE ENERGY OF SPUTTERED ATOM IN eV: '))
        self.T0 = float(input('ENTER GAS TEMPERATURE IN KELVIN: '))

        self.Tcorrect = input(
                'DO YOU WANT TO PERFORM TEMPERATURE CORRECTION - (Y)ES?') == 'Y'

        print('\n Let\'s talk about the calculation geomery.')
        self.in_rad = float(input('ENTER INNER RADIUS OF EROSION ZONE IN mm: '))
        self.out_rad = float(input('ENTER OUTER RADIUS IN mm: '))
        self.dist = float(input('TARGET SUBSTRATE DISTANCE Z IN mm: '))
        hstr = 'HOW FAR AWAY FROM CENTER OF SUBSTRATE DO YOU WANT TO CALCULATE\n'
        hstr += 'DEP RATE? (NUMBER OF OUTER EROSION ZONE RADII): '
        self.L5 = float(input(hstr))
        self.N1 = float(input('NUMBER OF STEPS FOR RADIAL INTEGRATION AT TARGET: '))
        self.N2 = float(input('NUMBER OF STEPS FOR POLAR INTEGRATION AT TARGET: '))
        self.grid_R2 = float(input('HOW MANY mm BETWEEN CALCULATED DATA POINTS? '))
        
        return(self)
        
    def __str__(self):
        ret_str = 'Sputtering calculation input parameters:\n'
        ret_str += 'Atom\tM\tZ\tR\tU\tQ\tDEN\tKAPA\n'
        ret_str += self.Gas_symbol +'\t'
        ret_str += str(self.MG) +'\t'
        ret_str += str(self.ZG) +'\t'
        ret_str += str(self.RG) +'\t'
        ret_str += str(self.UG) +'\t'
        ret_str += str(self.QG) +'\t'
        ret_str += str(self.DENG) +'\t'
        ret_str += str(self.KAPA) +'\n'
        ret_str += str(self.Target_symbol) +'\t'
        ret_str += str(self.MA) +'\t'
        ret_str += str(self.ZA) +'\t'
        ret_str += str(self.RA) +'\t'
        ret_str += str(self.U0) +'\t'
        ret_str += str(self.DEN) +'\t'
        ret_str += str(self.QZ) + '\n'
        
        ret_str += 'Discharge voltage (V): '+ str(self.disch_volt)
        ret_str += '\nDischarge current (A): '+ str(self.disch_cur)
        ret_str += '\nGas pressure (Pa): '+ str(self.pres)
        ret_str += '\nAverage energy of sputtered atom (eV): '+ str(self.Eavg)
        ret_str += '\nGas temperature (K): '+ str(self.T0)
        ret_str += '\nPerform T correction?: '
        ret_str += 'Yes' if self.Tcorrect else 'No'
        
        ret_str += '\nErosion disk inner radius (mm): '+ str(self.in_rad)
        ret_str += '\nErosion disk outer radius (mm): '+ str(self.out_rad)
        ret_str += '\nTarget to substrate distance(mm): '+ str(self.dist)
        ret_str += '\nCalculated radius / outer erosion radius: '+ str(self.L5)
        ret_str += '\nNumber of steps for radial integration at target: '
        ret_str +=  str(self.N1)
        ret_str += '\nNumber of steps for polar integration at target: '
        ret_str +=  str(self.N2)
        ret_str += '\nmm between calculated data points: '
        ret_str +=  str(self.grid_R2)
    
        return(ret_str)
        
    def toASCII(self, ofile):
        '''Write calc parameters into open ASCII file'''
        ofile.write(self.__str__())
        return()
    
    def fromASCII(self, ifile):
        '''Load data from open ASCII file'''
        data_lines = ifile.readlines()
        ln =data_lines[2].split() #Gas atomic params
        self.Gas_symbol = ln[0]
        self.MG  = float(ln[1])
        self.ZG = float(ln[2])
        self.RG = float(ln[3])
        self.UG = float(ln[4])
        self.QG = float(ln[5])
        self.DENG = float(ln[6])
        self.KAPA = float(ln[7])
        ln =data_lines[3].split() #Target atomic params
        self.Target_symbol = ln[0]
        self.MA = float(ln[1])
        self.ZA = float(ln[2])
        self.RA = float(ln[3])
        self.U0 = float(ln[4])
        self.DEN = float(ln[5])
        self.QZ = float(ln[6])
        
        self.disch_volt = float(data_lines[4].split(':')[-1]) 
        self.disch_cur = float(data_lines[5].split(':')[-1])
        self.pres = float(data_lines[6].split(':')[-1])
        self.Eavg = float(data_lines[7].split(':')[-1])
        self.T0 = float(data_lines[8].split(':')[-1])
        cT = data_lines[9].split(':')[-1]
        self.Tcorrect = cT == 'Yes'
        
        self.in_rad = float(data_lines[10].split(':')[-1])
        self.out_rad = float(data_lines[11].split(':')[-1])
        self.dist = float(data_lines[12].split(':')[-1])
        self.L5 = float(data_lines[13].split(':')[-1])
        self.N1 = float(data_lines[14].split(':')[-1])
        self.N2 = float(data_lines[15].split(':')[-1])
        self.grid_R2 = float(data_lines[16].split(':')[-1])
        
        return(self)
        
        

class Sputt_plan_circ():
    '''Class to perfrom the calculation.'''
    
    def __init__(self, inputData):
        '''A calculation can only be created by supplying data and
        it will be performed apon creation'''
        self.inputData = inputData
        self.mfp_lambda = []
        self.eta = [] 
        self.TC=[]
        self.LALA=[]
        self.R99 =[]
        self.Y=[]
        
        self.calculate() #perfroms the calculation
        
    def atom_radius(self, z):
        '''Interpolation of Robinson's data'''
        if (z < 39):
             r = np.sqrt(0.77*0.77 + ((0.93*0.93 - 0.77*0.77)/18) * (z - 18))
        else:
            r = np.sqrt(0.93*0.93 + ((1.04*1.04 - 0.93*0.93)/18) * (z - 36))
        return(r)

    def mean_free_path(self, T):    
        '''Mean free path of sputtered atoma'''
        m = self.inputData.MA/self.inputData.MG
        P = self.inputData.pres
        RG = self.inputData.RG
        RA = self.inputData.RA
        ret_lambda = 1/(2.276516*(P/T)*(RG+RA)*(RG+RA)*np.sqrt(1+m))
        self.mfp_lambda.append(ret_lambda)
        return(ret_lambda)

    def avg_collision_number_to_thermalize(self, T):
        '''calculates eta, avg energy of sputt atoms is 10 eV'''
        m = self.inputData.MG / self.inputData.MA
        E9 = self.inputData.Eavg
        #double vper;
        #double eta;
        vper = np.log(np.sqrt(1+m)+np.sqrt(m))/(4 * m**1.5 * np.sqrt(1+m))
        vper += (2* m**4 + 5* m**3 + 3*m*m - m -1)/(4*m* (1+m)**3)
        vper = (1-m)/(1+m) + 2 * m * vper/(1+m)
        eta = np.log(np.sqrt(8.6174E-5*T/(E9*m)))
        eta /= np.log(vper)
        self.eta.append(eta)
        return(eta)

    def collisionless_transport(self):
        '''Vectorized version'''
        
        self.inputData.R1 = self.inputData.in_rad
        self.inputData.R2 = self.inputData.out_rad

        dr=(self.inputData.R2 - self.inputData.R1) / self.inputData.N1
        dfi= np.pi/self.inputData.N2
        grid_R2 = self.inputData.grid_R2
        maxk1 = int(self.inputData.N1 -1)
        maxk2 = int(self.inputData.N2 -1)
        R1 = self.inputData.R1 
        Z = self.inputData.dist
        LAMBDA = self.mfp_lambda[-1]
       
        maxi=self.inputData.L5 * self.inputData.R2
        maxi /= self.inputData.grid_R2
        maxi=int(np.ceil(maxi))
          
        R = grid_R2 *np.arange(maxi+1).reshape(maxi+1,1)
        lR = len(R)
        R=R.reshape(lR,1)
        RC = R1 + dr * (0.5 + np.arange(maxk1+1))
        lRC=len(RC)
        RC = RC.reshape(lRC,1)
        fi = (np.arange(maxk2+1) +0.5) * dfi
        lfi=len(fi)
        fi=fi.reshape(lfi,1)
        
        R0 = np.dot(np.cos(fi), RC.T)
        #print(R0.shape)
        R0=R0.reshape(lfi, lRC, 1)
        #print(R0.shape)
        R0 = -2 *np.dot(R0,R.T)
        #print('R0 shape:', R0.shape)
        R2 = np.dot(np.ones((lfi, lRC, 1)), (R**2).T)
        RC2 =np.dot(np.ones((lfi, 1)), (RC**2).T)
        RC2 =np.dot(RC2.reshape(lfi, lRC, 1), np.ones((1,lR)))
        #print('R0', R0.shape)
        #print('R2', R2.shape)
        #print('RC2', RC2.shape)
        R0 += RC2 + R2 + Z**2
        
        eR0 = np.exp(-np.sqrt(R0)/LAMBDA)
        eR0 /= R0**2
        #print('e0 shape:', eR0.shape)
        S=np.dot(np.swapaxes(eR0,1,2), RC)
        #print('S: ', S.shape)
        S=np.sum(S,0).flatten()
        DEP1 = (Z*Z*dr*dfi*2/np.pi) * S
        self.DEP1 = DEP1
        return(DEP1)


    def collisionless_dep_rate(self):
        maxj=self.inputData.L5 * self.inputData.R2
        maxj=int(np.ceil(maxj / self.inputData.grid_R2))
        A = 2 * np.pi * self.inputData.grid_R2**2 *np.arange(maxj+1)
        A[0] = np.pi * self.inputData.grid_R2**2 / 4
        RTOT = np.dot(A, self.DEP1)
        R99 = RTOT / np.pi 
        R99 /= self.inputData.R2**2 - self.inputData.R1**2
        self.R99.append( R99)
        return(R99)


    def sputtering_yield(self):
        MA = self.inputData.MA
        ZA = self.inputData.ZA
        MG = self.inputData.MG
        ZG = self.inputData.ZG
        QA = self.inputData.QZ
        UA = self.inputData.U0
        EI = self.inputData.disch_volt
        
        H = 0.834 if MA/MG > 1 else 0.18
        TH = 1.5 * UA * np.power((1+ 1.38 * np.power(MG/MA,H)),2)
        TH /= 4*MA*MG/np.power((MA+MG),2)
        ALPHA = 0.1 + 0.155*np.power(MA/MG,0.73)
        PT = (1+MG/MA)*ZA*ZG
        PT *= np.sqrt(np.power(ZA,2.0/3)+np.power(ZG,2.0/3)) /0.0325
        CP = 3.56 * MG * ZG * ZA * ALPHA / ((MA+MG) * UA)
        CP /=  np.sqrt(np.power(ZA,2.0/3)+np.power(ZG,2.0/3))
        EPS = EI/PT
        SN = 3.441 * np.sqrt(EPS) * np.log(EPS+2.718)
        SN /= (1 + 6.355*np.sqrt(EPS)+EPS*(6.882*np.sqrt(EPS)-1.708))
        Y = QA * CP *SN * np.power(1-np.sqrt(TH/EI),2)
        self.Y.append(Y)
        return(Y)
        

    def correct_T(self, T):
        '''Recursive T correction. T is current guess'''
        #Upkeep needed even if no T correction wanted:
        T0 = self.inputData.T0
        LAMBDA = self.mean_free_path(T) #updates self.mfp_lambda
        ETA = self.avg_collision_number_to_thermalize(T) #updates self.eta
        LAMBDA *= ETA
        self.mfp_lambda[-1]=LAMBDA #used by collisionless_transport()
        LALA = LAMBDA /(10*T)
        self.LALA.append(LALA)
        self.collisionless_transport() # writes self.DEP1[]
        R99 = self.collisionless_dep_rate()
        Y = self.sputtering_yield()
        
        if not self.inputData.Tcorrect:  
            TC=self.inputData.T0
            self.TC.append(TC)
            return(TC)
        else:    
            CUR = self.inputData.disch_cur
            E9 = self.inputData.Eavg
            KAPPA = self.inputData.KAPA
            AA = CUR * E9 * Y * (1-R99) / (4*np.pi*KAPPA*LALA)
            TC = (T0 + np.sqrt(T0*T0 + 4*AA))/2
            self.TC.append(TC)
            if (np.abs(T-TC)>0.1): # Need to correct T?
                return(self.correct_T(TC))
            else:
                return(TC)   
    
    def difusion_source(self):
        '''Calculate currfent J(r,z)'''
        L5 = self.inputData.L5
        R1 = self.inputData.R1
        R2 = self.inputData.R2
        grid_R2 = self.inputData.grid_R2
        Z = self.inputData.dist
        LAMBDA = self.mfp_lambda[-1]
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
        self.J = J
        return(J)
        
    def diffusion_component(self):
        L5 = self.inputData.L5
        R1 = self.inputData.R1
        R2 = self.inputData.R2
        grid_R2 = self.inputData.grid_R2
        Z = self.inputData.dist
        LAMBDA = self.mfp_lambda[-1]
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
        
        #DEP21 = np.zeros(nr+1)
        r = grid_R2 * np.arange(nr+1) # r[iir]
        zc = dz * (np.arange(nz)+ 0.5) #zc[{iz}]
        rc = (np.arange(nr)+0.5)*drv # rc[{ir}]
        fi = (np.arange(N2)+0.5)*dfi #fi[{k2}]
        
        rop0 = rc.reshape(1,1,nr,1) * np.cos(fi).reshape((1,1,1,N2))
        rop1 = r.reshape(nr+1,1,1,1) * rop0
        rop0 = (r**2).reshape(nr+1, 1,1,1) + (rc**2).reshape(1,1,nr,1) -2 *rop1
        rop = np.sqrt(rop0)
        #rop0 = rop1 = False free some memory
                        
        A = (Z-zc.reshape((1,nz,1,1)))/(2*Z)
        B = (Z+zc.reshape((1,nz,1,1)))/(2*Z)
        C= rop /(2*Z)
        C2 = C**2
        kk1 = A/np.power(A**2 + C2, 1.5)
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
        
        self.DEP21 =DEP21
        
        A2= (0.0-zc.reshape((1,nz,1,1)))/(2*Z)
        B2= (0.0+zc.reshape((1,nz,1,1)))/(2*Z)
        s0= A2/np.power(A2**2 + C2, 1.5) 
        s0 -= B2/np.power(B2**2 +C2, 1.5)
        s1= (1+A2) /np.power((1+A2)**2 + C2, 1.5) 
        s1 += (1-B2) /np.power((1-B2)**2 +C2, 1.5)
        s1 -= (1-A2) /np.power((1-A2)**2 + C2, 1.5)  
        s1 -= (1+B2) /np.power((1+B2)**2 + C2, 1.5)
        s0 += s1
        s2 = J.T.reshape(1, nz, 1,nr) * rc.reshape(1,1,1,nr) 
        sz = np.matmul(s2, s0)
        DEP22 = 2* drv *dz *dfi / (16 * np.pi * Z**2 ) *np.sum(sz, axis=(-1,-2,-3))
        self.DEP22 = DEP22
        return(DEP21, DEP22)
               
        
    def calculate(self):
        self.inputData.RG = self.atom_radius(self.inputData.ZG)
        self.inputData.RA = self.atom_radius(self.inputData.ZA)
        self.TG = self.correct_T(self.inputData.T0) #Calculates lambda using eta
        self.difusion_source()
        self.diffusion_component()
        
        R1 = self.inputData.R1
        R2 = self.inputData.R2
        grid_R2 = self.inputData.grid_R2
        ir=int(np.ceil(self.inputData.L5 *R2/grid_R2)+1)
        CUR = self.inputData.disch_cur
        Y = self.Y[-1]
        
        self.pR = np.arange(ir) * self.inputData.grid_R2
        
        # renormalize dep rates form fraction of sputter current at target
            # to nm per sec
        c = CUR * Y / (np.pi * ((R2*R2 - R1*R1)/100) * 1.6022E-19)
        c *= (self.inputData.MA / self.inputData.DEN) * (1.0E7 / 6.022E23)
        self.DEP1 *= c
        self.DEP21 *= c
        self.pTot = self.DEP1+ self.DEP21
        self.DEP22 *= -1
            
        return(self)
    

    



        