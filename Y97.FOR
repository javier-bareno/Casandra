
C
C	DEPOSITION RATE CALCULATIONS FOR PLANAR MAGNETRONS
C
C	DECLARATIONS
	INTEGER NZ
	REAL MA,MG,RG,RA,P,TG,R1,R2,DR,DFI,Z,M,S,ETA,LAMBDA,DRV
	REAL A,B,C,DEP1(30),J(100,100),DEP2(2,30),DEP(30),LALA,KAPA
	REAL QZ
	CHARACTER*2 IO,EL
	CHARACTER*6 FIL,COEFIL
	CHARACTER*1 W,WW
	PARAMETER(PI=3.1415927)
C
C	INTERACTIVE INPUT OF DATA
C
 1	CONTINUE
	TYPE*,'DEPOSITION RATE CALCULATIONS FOR PLANAR MAGNETRONS'
	TYPE*,'FOR INERT GASES AND ELEMENTAL TARGETS'
	TYPE*,'  '
	TYPE*,'ENTER CHEMICAL SYMBOL OF SPUTTERING GAS'
        TYPE*,'note: a coefficient file for the element has to be'
	TYPE*,'created first by running coe.exe'
	READ(*,2) IO
	TYPE*,'CHEMICAL SYMBOL OF TARGET ELEMENT'
	READ(*,2) EL
 2	FORMAT(A)
C	RETRIEVAL OF COEFFICIENT FILES
	COEFIL=IO//'.COE'
	OPEN(UNIT=1,FILE=COEFIL,TYPE='OLD')
	READ(1,*) MG,ZG,RG,UG,Q,DEN,KAPA
	CLOSE(UNIT=1)
	KAPA=KAPA*1E-6
C	CALCULATION OF GAS ATOM RADIUS FROM ROBINSON DATA AT 5 EV
C	IN THIS VERSION THE VALUE FOR GAS ATOM RADIUS IN THE COEFFICIENT
C	FILE IS NOT USED
	RG=SQRT(0.77*0.77 +(0.93*0.93-0.77*0.77)/18*(ZG-18))
	IF (ZG.LT.39) GOTO 4
	RG=SQRT(0.93*0.93+(1.04*1.04-0.93*0.93)/18*(ZG-36))
4	TYPE*,'ION=',IO,MG,ZG,RG,KAPA
	COEFIL=EL//'.COE'
	OPEN(UNIT=1,FILE=COEFIL,TYPE='OLD')
	READ(1,*) MA,ZA,RA,U0,QZ,DEN
	CLOSE(UNIT=1)	
 	OPEN (UNIT=4,NAME='DEPRATE.OUT',TYPE='NEW')
C	CALCULATION OF SPUTTERED ATOM RADIUS BY INTERPOLATION OF ROBINSON'S DATA
C	AT 5 eV 
	RA=SQRT(0.77*0.77+(0.93*0.93-0.77*0.77)/18*(ZA-18))
	IF (ZA.LT.39) GOTO 6
	RA=SQRT(0.93*0.93+(1.04*1.04-0.93*0.93)/18*(ZA-36))
6	TYPE*,'TARGET = ',EL,MA,ZA,RA,U0,QZ,DEN
	TYPE*,'DISCHARGE VOLTAGE  IN VOLTS  '
	ACCEPT*,EI
	TYPE*,'ENTER SPUTTERED ATOM AVERAGE ENERGY'
        ACCEPT*,E9
	TYPE*,'DISCHARGE CURRENT IN AMPS'
	ACCEPT*,CUR
	TYPE*,'GAS PRESSURE IN PASCALS'
	ACCEPT*,P
	TYPE*,'GAS TEMPERATURE IN KELVIN'
	ACCEPT*,T0
	TYPE*,'DO YOU WANT TEMPERATURE CORRECTION - (Y)ES?'
	READ(*,2) WW
	TYPE*,'INPUT ACTIVE TARGET AREA IN MILIMETERS'
	TYPE*,'INNER RADIUS OF EROSION ZONE IN mm R1'
	ACCEPT*,R1
	TYPE*,'OUTER RADIUS R2'
	ACCEPT*,R2
	TYPE*,'INPUT N1, (number of segments for intergration over target)'
	ACCEPT*,N1
	TYPE*,'INPUT RADIUS OF CALCULATED VOLUME - ?*R2'
	ACCEPT*,L5
	NR=L5*R2/10
	TYPE*,'dFI=PI/N2   ENTER N2'
	ACCEPT*,N2
	TYPE*,'TARGET SUBSTRATE DISTANCE Z IN mm'
	ACCEPT*,Z
	NZ=Z/2
C
C	CALCULATION OF THE MEAN FREE PATH IN mm
C
	TG=T0
	M=MG/MA
 3	LAMBDA=1/(2.276516*P/TG*(RG+RA)**2*SQRT(1+MA/MG))
	TYPE*,' LAMBDA = ',LAMBDA
C	
C	CALCULATION OF AVERAGE NUMBER OF COLLISIONS BEFORE THERMALIZATION'
C       BASED ON WESTWOOD'S PAPER
C
	VPER=LOG(SQRT(1+M)+SQRT(M))/4/M**1.5/SQRT(1+M)
	VPER=VPER+(2*M**4+5*M**3+3*M*M-M-1)/4/M/(1+M)**3
	VPER=(1-M)/(1+M)+2*M/(1+M)*VPER
	ETA=LOG(SQRT(8.6174E-5*TG/E9/M))
	TYPE*,' LOG(VG/V0)=',ETA,'   V1/V=',VPER
	ETA=ETA/LOG(VPER)
C 	IN THE EXPRESSION FOR ETA WE HAVE FIXED THE AVERAGE ENERGY OF SPUTTERED
C	ATOMS TO BE E9 eV
	TYPE*,' ETA=',ETA,'   M=',M
	LAMBDA=LAMBDA*ETA
C	THIS FACTOR IS USED IN THE TEMPERATURE CORRECTION IN cm/K
	LALA=LAMBDA/TG/10
	TYPE*,'LAMBDA*ETA=',LAMBDA

C
C	CALCULATION OF COLLISIONLESS TRANSPORT  	
C	
	DR=(R2-R1)/N1
	DFI=PI/N2	
	DO 40 R=0,2*R2,10
	S=0
	   	DO 30 K1=0,N1-1,1
	   	RC=R1+(K1+.5)*DR
			DO 20 K2=0,N2-1,1
			FI=(K2+.5)*DFI
			RO=Z*Z+R*R+RC*RC-2*RC*R*COS(FI)
			S=S+RC*EXP(-SQRT(RO)/LAMBDA)/(RO*RO)
 20	CONTINUE
 30	CONTINUE
	DEP1(R/10+1)=Z*Z*DR*DFI*S*2/PI
 	TYPE*,R,DEP1(R/10+1),DEP1(R/10+1)/DEP1(1)
 40	CONTINUE
C	
C	CALCULATION OF THE TOTAL COLISIONLESS DEPOSITION RATE
C
	RTOT=DEP1(1)
	DO 45 I=2,2*R2/10+1,1
	RTOT=RTOT+DEP1(I)*(I-1)*8
 45 	CONTINUE
	RTOT=RTOT*PI*25
	R99=RTOT/PI/(R2**2-R1**2)
	TYPE*,'RDTOT/RSTOT',R99
C
C	SPUTTERING YIELD FROM YAMAMURA'S FORMULA
C
	IF (MA/MG.GT.1) THEN
	H=.834
	ELSE 
	H=.18
 	END IF
	TH=1.5*U0*(1+1.38*(MG/MA)**H)**2/(4*MA*MG/(MA+MG)**2)
	ALFA=.1+.155*(MA/MG)**.73
	PT=(1+MG/MA)*ZA*ZG*SQRT(ZA**.66666+ZG**.666666)/.0325
	CP=3.56*MG*ZG*ZA*ALFA/(MA+MG)/SQRT(ZA**.666666+ZG**.66666)/U0
	EPS=EI/PT
	SN=3.441*SQRT(EPS)*LOG(EPS+2.718)
	SN=SN/(1+6.355*SQRT(EPS)+EPS*(6.882*SQRT(EPS)-1.708))
	Y=QZ*CP*SN*(1-SQRT(TH/EI))**2
	TYPE*,'E=',EI,'   Y=',Y,' TH=',TH,' ALFA=',ALFA
	TYPE*,'EPS=',EPS,' SN=',SN
C
C	TEMPERATURE CORRECTION (MODIFIED VERSION OF ROSSNAGEL'S APPROACH)
C
	IF (WW.NE.'Y') GOTO 8
	TYPE*,'KAPA=',KAPA,' LALA=',LALA
	AA=CUR*E9*Y*(1-R99)/4/PI/KAPA/LALA
	T1=(T0+SQRT(T0**2+4*AA))/2
	IF (ABS(TG-T1).LT.1) GOTO 8
	TG=T1
	TYPE*,'TG=',TG,TG,TG,TG
	GOTO 3
8	TYPE*,'FINAL TEMPERATURE IS ',TG,' K'	
C
C	ENSURING THAT DZ<LAMBDA
C
	DZ=Z/NZ
	IF(DZ.GT.LAMBDA) GOTO 9
	GOTO 10
 9	NZ=Z/LAMBDA
	TYPE*,'NEW VALUE OF NZ--->',NZ
	DZ=Z/NZ
 10	CONTINUE
	TYPE*,'NZ=',NZ
C
C	DIFUSION SOURCE CALCULATION
C
	DRV=L5*R2/NR
	DO 80 IR=0,NR-1,1
	R=(IR+.5)*DRV
	  DO 70 IZ=0,NZ-1,1
	  ZC=(IZ+.5)*DZ
	  S=0
	    DO 65 K1=0,N1-1,1
            RC=R1+(K1+.5)*DR
	      DO 60 K2=0,N2-1,1
	      FI=(K2+.5)*DFI
	      RO=SQRT(ZC**2+R**2+RC**2-2*RC*R*COS(FI))
	      S=S+RC*ZC*EXP(-RO/LAMBDA)/LAMBDA/RO**3
 60	CONTINUE	 	
 65	CONTINUE
	J(IR+1,IZ+1)=DR*DFI*S*2/PI
 70 	CONTINUE
 80	CONTINUE
C	
C	DIFFUSION COMPONENT OF THE DEPOSITION RATE
C
	Z9=Z
	I9=1
 100 	CONTINUE
	DO 200 R=0,2*R2,10
	SR=0
	 DO 190 IZ=0,NZ-1,1
	 ZC=(IZ+.5)*DZ
	 SZ=0
	  DO 180 IR=0,NR-1,1
	  RC=(IR+.5)*DRV
	   DO 170 K2=0,N2-1,1
  	   FI=(K2+.5)*DFI
	   ROP=SQRT(R*R+RC*RC-2*RC*R*COS(FI))
	   A=(Z9-ZC)/(2*Z)
	   B=(Z9+ZC)/(2*Z)
	   C=ROP/2/Z
	   S0=A/(A*A+C*C)**1.5-B/(B*B+C*C)**1.5
	    DO 160 I=1,1,1
	    S1=(I+A)/((I+A)*(I+A)+C*C)**1.5+(I-B)/((I-B)*(I-B)+C*C)**1.5
	    S1=S1-(I-A)/((I-A)*(I-A)+C*C)**1.5-(I+B)/((I+B)*(I+B)+C*C)**1.5
	    S0=S0+S1
 160	    CONTINUE
	   S2=J(IR+1,IZ+1)*DRV*DZ*DFI*RC*S0/16/Z/Z/PI
	   SZ=SZ+S2
 170       CONTINUE
 180      CONTINUE
	 SR=SR+SZ
 190	 CONTINUE
	DEP2(I9,R/10+1)=2*SR
	TYPE*,'R=',R,'  DEP2=',DEP2(I9,R/10+1)
 200	CONTINUE
	IF (I9.EQ.2) GOTO 210
	I9=2
	Z9=0
	GO TO 100	
 210    CONTINUE
	
C
C	RESULTS
C
	DO 300 R=0,2*R2,10
	DEP(R/10+1)=DEP1(R/10+1)+DEP2(1,R/10+1)
	TYPE*,'R=',R,'mm  DSUB=',DEP(R/10+1),' DTAR=',DEP2(2,R/10+1)
	TYPE*,'DIRECT FLUX=',DEP1(R/10+1),'  DIFUSED=',DEP2(1,R/10+1)
 300	CONTINUE
C
C	TOTAL DEPOSITION AND REDEPOSITION
C	DEPOSITION
	DTOT=DEP(1)
	DO 245 I=2,2*R2/10+1,1
	DTOT=DTOT+DEP(I)*(I-1)*8
 245 	CONTINUE
	DTOT=DTOT*PI*25
	TYPE*,'DTOT/SPTOT',DTOT/PI/(R2**2-R1**2)
C	REDEPOSITION
	PTOT=DEP2(2,1)
	DO 255 I=2,2*R2/10+1,1
	PTOT=PTOT+DEP2(2,I)*(I-1)*8
 255 	CONTINUE
	PTOT=PTOT*PI*25
	TYPE*,'BACKTOT/SPTOT',PTOT/PI/(R2**2-R1**2)
C	REDEPOSITION ON MAGNETRON ACTIVE AREA
	ACTOT=0
	DO 257 I=R1/10+1,R2/10,1
	ACTOT=ACTOT+DEP2(2,I)*(I-1)*8
 257 	CONTINUE
	ACTOT=ACTOT*PI*25
	TYPE*,'BACKTOT ON ACTIVE AREA/SPTOT',ACTOT/PI/(R2**2-R1**2)
	TYPE*,'LAMBDA=',LAMBDA
C
C	OUTPUT OF RESULTS
C
260	CONTINUE
	TYPE*,'DO YOU WANT THESE RESULTS IN A FILE? (Y)'
	TYPE*,'IF YOU WANT TO END -- INPUT  E, TO START ANEW - RETURN'
	READ(*,2) W
	IF(W.EQ.'E') GOTO 315
	IF(W.NE.'Y') GOTO 1
	WRITE(4,312)'   '
	WRITE(4,304)'MG=',MG,' ZG=',ZG,' RG=',RG,' THERM.COND.=',KAPA
	WRITE(4,304)'MA=',MA,' ZA=',ZA,' RA=',RA,' U0=',U0
	WRITE(4,302)'LAMBDA=',LAMBDA/ETA,' ,10eV,TERM.DIST.=',LAMBDA
 	WRITE(4,302)'ION CURRENT=',CUR,'A VOLTAGE=',EI
	WRITE(4,302)'PRESSURE=',P,'PA Y=',Y
	WRITE(4,304)'TG=',TG,'K R2=',R2,' R1=',R1,' Z=',Z
	WRITE(4,307)'N1=',N1,' N2=',N2,' NR=',NR,' NZ=',NZ
	DO 350 R=0,2*R2,10
	WRITE(4,303)'R=',R,'mm  DSUB=',DEP(R/10+1),' DTAR=',DEP2(2,R/10+1)
	WRITE(4,302)'DIRECT FLUX=',DEP1(R/10+1),' (1)DIFUSED=',DEP2(1,R/10+1)
 350	CONTINUE
	WRITE(4,312)'   '
	WRITE(4,301)'TOTAL BACKDIFFUSED/TOTAL SPUTTERED',PTOT/PI/(R2**2-R1**2)
	WRITE(4,301)'TOTAL DEPOSITED/TOTAL SPUTTERED',DTOT/PI/(R2**2-R1**2)
	PP=(PTOT-ACTOT)/PI/(R2**2-R1**2)
	WRITE(4,301)'BACKTOT ON INACTIVE AREA/SPTOT',PP
	TYPE*,CUR*MA*Y/(R2**2-R1**2)/DEN*32989.4
	TYPE*,CUR*Y/(R2**2-R1**2)*2E20

	WRITE(4,301)'SRATE IN ANSTROMS/SEC=',CUR*MA*Y/(R2**2-R1**2)/DEN*32989.4
	WRITE(4,311)'SRATE IN AT/CM**2/SEC=',CUR*Y/(R2**2-R1**2)*2E20
	WRITE(4,312)'   '
 301 	FORMAT(' ',A,F14.7)
 302	FORMAT(' ',A,F14.7,A,F14.7)
 303 	FORMAT(' ',A,F14.7,A,F14.7,A,F14.7)
 304    FORMAT(' ',A,F14.7,A,F14.7,A,F14.7,A,F14.7)	
 307	FORMAT(' ',A,I3,A,I3,A,I3,A,I3)
 311 	FORMAT(' ',A,E14.7)
 312	FORMAT(' ',A)
	GOTO 260
 315	CONTINUE
	CLOSE(UNIT=4)
	STOP
	END
