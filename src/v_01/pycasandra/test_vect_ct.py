## test of vectorization of collisionless transport
import numpy as np
import time

def loops(maxi, grid_R2, maxk1, maxk2, R1, dr, dfi, Z, LAMBDA ):
    DEP1 = []
    for i in range(maxi+1):
        R = grid_R2*i
        S =0
        for k1 in range(maxk1 +1):
            RC = R1 + (k1 + 0.5) * dr
            for k2 in range(maxk2+1):
                fi = (k2+0.5)*dfi
                R0 = Z*Z + R*R + RC*RC - 2*RC*R*np.cos(fi)
                S += RC * np.exp(-np.sqrt(R0)/LAMBDA)/(R0*R0)
        DEP1.append(Z*Z*dr*dfi*S*2/np.pi)
    return(DEP1)
    
def vect(maxi, grid_R2, maxk1, maxk2, R1, dr, dfi, Z, LAMBDA ):
    #DEP1 = []
    
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
    print(R0.shape)
    R0=R0.reshape(lfi, lRC, 1)
    print(R0.shape)
    R0 = -2 *np.dot(R0,R.T)
    print('R0 shape:', R0.shape)
    R2 = np.dot(np.ones((lfi, lRC, 1)), (R**2).T)
    RC2 =np.dot(np.ones((lfi, 1)), (RC**2).T)
    RC2 =np.dot(RC2.reshape(lfi, lRC, 1), np.ones((1,lR)))
    print('R0', R0.shape)
    print('R2', R2.shape)
    print('RC2', RC2.shape)
    R0 += RC2 + R2 + Z**2
    
    eR0 = np.exp(-np.sqrt(R0)/LAMBDA)
    eR0 /= R0**2
    print('e0 shape:', eR0.shape)
    S=np.dot(np.swapaxes(eR0,1,2), RC)
    print('S: ', S.shape)
    S=np.sum(S,0).flatten()
    DEP1 = (Z*Z*dr*dfi*2/np.pi) * S
    return(DEP1)

if __name__ == "__main__":
     t0 = time.time()
     l=loops(100,2,50, 30, 20, 2, 3, 100, 30)
     t1=time.time()
     #print(l)
     m=vect(100,2,50, 30, 20, 2, 3, 100, 30)
     t2=time.time()
     print(m.shape)
     #print(np.array(l)/m)
     print(t1-t0)
     print(t2-t1)