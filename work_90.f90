{\rtf1\ansi\ansicpg1252\cocoartf2511
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\margl1440\margr1440\vieww10800\viewh8400\viewkind0
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs24 \cf0 !*==AA0001.spg  processed by SPAG 6.72Dc at 01:14 on 17 Dec 2019\
      IMPLICIT NONE\
!*--AA00013\
!*** Start of declarations inserted by SPAG\
      REAL aa , alfa , cp , cur , den , dz , ei , eps , fi , h , PI ,   &\
         & pt , q , r , r99 , rc , ro , rop , rtot , s0\
      REAL s1 , s2 , sn , sr , sz , t0 , t1 , th , u0 , ug , vper , y , &\
         & z9 , za , zc , zg\
      INTEGER i , ir , iz , jr , jz , k1 , k2 , n1 , n2 , nr\
!*** End of declarations inserted by SPAG\
!	GAS PHASE DENSITY OF SPUTTERED ATOMS\
!	DEPOSITION RATE CALCULATIONS FOR PLANAR MAGNETRONS\
!\
!	DECLARATIONS\
      INTEGER nz\
      REAL ma , mg , rg , ra , p , tg , r1 , r2 , dr , dfi , z , m , s ,&\
         & eta , lambda , drv\
      REAL a , b , c , dep1(30) , j(100,100) , dep2(2,30) , dep(30) ,   &\
         & lala , kapa\
      REAL qz , conc(100,100) , dif\
      CHARACTER*2 io , el\
      CHARACTER*6 fil , coefil\
      CHARACTER*1 w , ww\
      PARAMETER (PI=3.1415927)\
!\
!	INTERACTIVE INPUT OF DATA\
!\
 100  PRINT * , 'DEPOSITION RATE CALCULATIONS FOR PLANAR MAGNETRONS'\
      PRINT * , 'ENTER CHEMICAL SYMBOL OF SPUTTERING GAS'\
      READ (*,99006) io\
      PRINT * , 'CHEMICAL SYMBOL OF TARGET ELEMENT'\
      READ (*,99006) el\
!	RETRIEVAL OF COEFFICIENT FILES\
      coefil = io//'.COE'\
      OPEN (UNIT=1,FILE=coefil,TYPE='OLD')\
      READ (1,*) mg , zg , rg , ug , q , den , kapa\
      CLOSE (UNIT=1)\
      kapa = kapa*1E-6\
!	CALCULATION OF GAS ATOM RADIUS FROM ROBINSON DATA AT 10/2 EV\
      rg = SQRT(0.77*0.77+(0.93*0.93-0.77*0.77)/18*(zg-18))\
      IF ( zg.GE.39 ) rg = SQRT(0.93*0.93+(1.04*1.04-0.93*0.93)         &\
                         & /18*(zg-36))\
      PRINT * , 'ION=' , io , mg , zg , rg , kapa\
      coefil = el//'.COE'\
      OPEN (UNIT=1,FILE=coefil,TYPE='OLD')\
      READ (1,*) ma , za , ra , u0 , qz , den\
      CLOSE (UNIT=1)\
      OPEN (UNIT=4,NAME='DEPRATE.OUT',TYPE='NEW')\
!	CALCULATION OF SPUTTERED ATOM RADIUS BY INTERPOLATION OF ROBINSON'S DATA\
!	AT 10/2 eV\
      ra = SQRT(0.77*0.77+(0.93*0.93-0.77*0.77)/18*(za-18))\
      IF ( za.GE.39 ) ra = SQRT(0.93*0.93+(1.04*1.04-0.93*0.93)         &\
                         & /18*(za-36))\
      PRINT * , 'TARGET = ' , el , ma , za , ra , u0 , qz , den\
      PRINT * , 'DISCHARGE VOLTAGE  IN VOLTS  '\
      READ * , ei\
      PRINT * , 'DISCHARGE CURRENT IN AMPS'\
      READ * , cur\
      PRINT * , 'GAS PRESSURE IN PASCALS'\
      READ * , p\
      PRINT * , 'GAS TEMPERATURE IN KELVIN'\
      READ * , t0\
      PRINT * , 'DO YOU WANT TEMPERATURE CORRECTION - (Y)ES?'\
      READ (*,99006) ww\
      PRINT * , 'INPUT ACTIVE TARGET AREA IN MILIMETERS'\
      PRINT * , 'INNER RADIUS OF EROSION ZONE IN mm R1'\
      READ * , r1\
      PRINT * , 'OUTER RADIUS R2'\
      READ * , r2\
      PRINT * , 'INPUT N1'\
      READ * , n1\
!	TYPE*,'INPUT RADIUS OF CALCULATED VOLUME - ?*R2'\
!	ACCEPT*,L5\
      nr = 3*r2/2\
      PRINT * , 'dFI=PI/N2   ENTER N2'\
      READ * , n2\
      PRINT * , 'TARGET SUBSTRATE DISTANCE Z IN mm'\
      READ * , z\
      nz = z/2\
!\
!	CALCULATION OF THE MEAN FREE PATH IN mm\
!\
      tg = t0\
      m = mg/ma\
 200  lambda = 1/(2.276516*p/tg*(rg+ra)**2*SQRT(1+ma/mg))\
      PRINT * , ' LAMBDA = ' , lambda\
!\
!	CALCULATION OF AVERAGE NUMBER OF COLLISIONS BEFORE THERMALIZATION'\
!\
      vper = LOG(SQRT(1+m)+SQRT(m))/4/m**1.5/SQRT(1+m)\
      vper = vper + (2*m**4+5*m**3+3*m*m-m-1)/4/m/(1+m)**3\
      vper = (1-m)/(1+m) + 2*m/(1+m)*vper\
      eta = LOG(SQRT(8.6174E-5*tg/10/m))\
      PRINT * , ' LOG(VG/V0)=' , eta , '   V1/V=' , vper\
      eta = eta/LOG(vper)\
! 	IN THE EXPRESSION FOR ETA WE HAVE FIXED THE AVERAGE ENERGY OF SPUTTERED\
!	ATOMS TO BE 10 eV\
      PRINT * , ' ETA=' , eta , '   M=' , m\
      lambda = lambda*eta\
      lala = lambda/tg/10\
      PRINT * , 'LAMBDA*ETA=' , lambda\
 \
!\
!	CALCULATION OF COLLISIONLESS TRANSPORT\
!\
      dr = (r2-r1)/n1\
      dfi = PI/n2\
      DO r = 0 , 2*r2 , 10\
         s = 0\
         DO k1 = 0 , n1 - 1 , 1\
            rc = r1 + (k1+.5)*dr\
            DO k2 = 0 , n2 - 1 , 1\
               fi = (k2+.5)*dfi\
               ro = z*z + r*r + rc*rc - 2*rc*r*COS(fi)\
               s = s + rc*EXP(-SQRT(ro)/lambda)/(ro*ro)\
            ENDDO\
         ENDDO\
         dep1(r/10+1) = z*z*dr*dfi*s*2/PI\
         PRINT * , r , dep1(r/10+1) , dep1(r/10+1)/dep1(1)\
      ENDDO\
!\
!	CALCULATION OF THE TOTAL COLISIONLESS DEPOSITION RATE\
!\
      rtot = dep1(1)\
      DO i = 2 , 2*r2/10 + 1 , 1\
         rtot = rtot + dep1(i)*(i-1)*8\
      ENDDO\
      rtot = rtot*PI*25\
      r99 = rtot/PI/(r2**2-r1**2)\
      PRINT * , 'RDTOT/RSTOT' , r99\
!\
!	SPUTTERING YIELD\
!\
      IF ( ma/mg.GT.1 ) THEN\
         h = .834\
      ELSE\
         h = .18\
      ENDIF\
      th = 1.5*u0*(1+1.38*(mg/ma)**h)**2/(4*ma*mg/(ma+mg)**2)\
      alfa = .1 + .155*(ma/mg)**.73\
      pt = (1+mg/ma)*za*zg*SQRT(za**.66666+zg**.666666)/.0325\
      cp = 3.56*mg*zg*za*alfa/(ma+mg)/SQRT(za**.666666+zg**.66666)/u0\
      eps = ei/pt\
      sn = 3.441*SQRT(eps)*LOG(eps+2.718)\
      sn = sn/(1+6.355*SQRT(eps)+eps*(6.882*SQRT(eps)-1.708))\
      y = qz*cp*sn*(1-SQRT(th/ei))**2\
      PRINT * , 'E=' , ei , '   Y=' , y , ' TH=' , th , ' ALFA=' , alfa\
      PRINT * , 'EPS=' , eps , ' SN=' , sn\
!\
!	TEMPERATURE CORRECTION\
!\
      IF ( ww.EQ.'Y' ) THEN\
         PRINT * , 'KAPA=' , kapa , ' LALA=' , lala\
         aa = cur*10*y*(1-r99)/4/PI/kapa/lala\
         t1 = (t0+SQRT(t0**2+4*aa))/2\
         IF ( ABS(tg-t1).GE.1 ) THEN\
            tg = t1\
            PRINT * , 'TG=' , tg , tg , tg , tg\
            GOTO 200\
         ENDIF\
      ENDIF\
      PRINT * , 'FINAL TEMPERATURE IS ' , tg , ' K'\
!\
!	ENSURING THAT DZ<LAMBDA\
!\
      dz = z/nz\
      IF ( dz.GT.lambda ) THEN\
         nz = z/lambda\
         PRINT * , 'NEW VALUE OF NZ--->' , nz\
         dz = z/nz\
      ENDIF\
      PRINT * , 'NZ=' , nz\
!\
!	DIFUSION SOURCE CALCULATION\
!\
      drv = 3*r2/nr\
      DO ir = 0 , nr - 1 , 1\
         r = (ir+.5)*drv\
         DO iz = 0 , nz - 1 , 1\
            zc = (iz+.5)*dz\
            s = 0\
            DO k1 = 0 , n1 - 1 , 1\
               rc = r1 + (k1+.5)*dr\
               DO k2 = 0 , n2 - 1 , 1\
                  fi = (k2+.5)*dfi\
                  ro = SQRT(zc**2+r**2+rc**2-2*rc*r*COS(fi))\
                  s = s + rc*zc*EXP(-ro/lambda)/lambda/ro**3\
               ENDDO\
            ENDDO\
            j(ir+1,iz+1) = dr*dfi*s*2/PI\
         ENDDO\
      ENDDO\
!\
!	CONCENTRATION CALCULATIONS\
 \
!\
      DO jz = 0 , nz , 1\
         z9 = jz*dz\
         DO jr = 0 , r2 , 1\
            r = jr*2\
            sr = 0\
            DO iz = 0 , nz - 1 , 1\
               zc = (iz+.5)*dz\
               sz = 0\
               DO ir = 0 , nr - 1 , 1\
                  rc = (ir+.5)*drv\
                  DO k2 = 0 , n2 - 1 , 1\
                     fi = (k2+.5)*dfi\
                     rop = SQRT(r*r+rc*rc-2*rc*r*COS(fi))\
                     a = (z9-zc)/(2*z)\
                     b = (z9+zc)/(2*z)\
                     c = rop/2/z\
                     s0 = 1./(a*a+c*c)**0.5 - 1./(b*b+c*c)**0.5\
                     DO i = 1 , 1 , 1\
                        s1 = 1./((i+a)*(i+a)+c*c)                       &\
                           & **0.5 - 1./((i-b)*(i-b)+c*c)**0.5\
                        s1 = s1 + 1./((i-a)*(i-a)+c*c)                  &\
                           & **0.5 - 1./((i+b)*(i+b)+c*c)**0.5\
                        s0 = s0 + s1\
                     ENDDO\
                     s2 = j(ir+1,iz+1)*drv*dz*dfi*rc*s0/8/z/PI\
                     sz = sz + s2\
                  ENDDO\
               ENDDO\
               sr = sr + sz\
            ENDDO\
            conc(jz+1,jr+1) = 2*sr\
            PRINT * , 'Z=' , z9 , 'R=' , r , '  CONC=' , conc(jz+1,jr+1)\
         ENDDO\
      ENDDO\
!\
!\
!\
!	OUTPUT OF RESULTS\
!\
 300  PRINT * , 'DO YOU WANT THESE RESULTS IN A FILE? (Y)'\
      PRINT * , 'IF YOU WANT TO END -- INPUT  E, TO START ANEW - RETURN'\
      READ (*,99006) w\
      IF ( w.EQ.'E' ) THEN\
         CLOSE (UNIT=4)\
      ELSE\
         IF ( w.NE.'Y' ) GOTO 100\
         WRITE (4,99009) '   '\
         WRITE (4,99008) 'MG=' , mg , ' ZG=' , zg , ' RG=' , rg ,       &\
                        &' THERM.COND.=' , kapa\
         WRITE (4,99008) 'MA=' , ma , ' ZA=' , za , ' RA=' , ra ,       &\
                       & ' U0=' , u0\
         WRITE (4,99007) 'LAMBDA=' , lambda/eta , ' ,10eV,TERM.DIST.=' ,&\
                       & lambda\
         WRITE (4,99007) 'ION CURRENT=' , cur , 'A VOLTAGE=' , ei\
         WRITE (4,99007) 'PRESSURE=' , p , 'PA Y=' , y\
         WRITE (4,99008) 'TG=' , tg , 'K R2=' , r2 , ' R1=' , r1 ,      &\
                       & ' Z=' , z\
         WRITE (4,99001) 'N1=' , n1 , ' N2=' , n2 , ' NR=' , nr ,       &\
                       & ' NZ=' , nz\
99001    FORMAT (' ',A,I3,A,I3,A,I3,A,I3)\
         DO r = 0 , 2*r2 , 10\
            WRITE (4,99002) 'R=' , r , 'mm  DSUB=' , dep(r/10+1) ,      &\
                           &' DTAR=' , dep2(2,r/10+1)\
99002       FORMAT (' ',A,F14.7,A,F14.7,A,F14.7)\
            WRITE (4,99007) 'DIRECT FLUX=' , dep1(r/10+1) ,             &\
                           &' (1)DIFUSED=' , dep2(1,r/10+1)\
         ENDDO\
         WRITE (4,99009) '   '\
 \
         DO ir = 0 , r2 , 1\
            r = ir*2\
            DO iz = 0 , nz , 1\
               zc = iz*dz\
               WRITE (4,99003) r , zc , conc(iz+1,ir+1)\
99003          FORMAT (' 'F14.7,F14.7,F14.7,F14.7)\
               PRINT * , conc(iz+1,ir+1)\
            ENDDO\
         ENDDO\
99004    FORMAT (' ',A,F14.7)\
99005    FORMAT (' ',A,E14.7)\
 \
         GOTO 300\
      ENDIF\
99006 FORMAT (A)\
99007 FORMAT (' ',A,F14.7,A,F14.7)\
99008 FORMAT (' ',A,F14.7,A,F14.7,A,F14.7,A,F14.7)\
99009 FORMAT (' ',A)\
      END\
}