	program gerSigQCD_CT14lo
c
c Gera Sig_QCD
c
c Obs1: Eq.9 do PRD92,074039(2015)
c
c Obs2: Calcula a também a parte imaginária
c
c Obs3: Usa a RDD 
c
c-----------------------------------------------------------------------
c Importante: Foi minimizado com uma parametrização real
c-----------------------------------------------------------------------
c
	implicit none
	double precision b1,b2,b3,b4,b5,b6,b7,b8,b9,b10
        double precision ssroot0,ssroot,deltas,ss
        double precision pi,mbfactor
        double precision ReSigQCDGeV2p100_Der,ReSigQCDGeV2m100_Der
        double precision ImSigQCDGeV2p100_Der,ImSigQCDGeV2m100_Der
        double precision d1,d2,d3
        double precision Xp
        double precision ww,w1,w2,w3,w4
        integer npts,iesse
c  
	open(unit=1,file='gerReSigQCDder_CT14lo.dat')
	open(unit=2,file='gerImSigQCDder_CT14lo.dat')
c
c--------------------------------------------
c Definição das constantes:
c      
	pi=4.d0*datan(1.d0)
        mbfactor=0.389379d0
c
ccc 10mil Calls
c
c FCN=   10.86528     FROM MIGRAD    STATUS=CONVERGED    128 CALLS     8901 TOTAL
c                     EDM=  0.28E-04    STRATEGY= 1      ERR MATRIX NOT POS-DEF
c
c  EXT PARAMETER                APPROXIMATE        STEP         FIRST   
c  NO.   NAME        VALUE          ERROR          SIZE      DERIVATIVE 
c   1     A1         99.155        1.7246       0.65540E-04  -0.25358E-02
c   2     A2        0.47865E-03   0.16683E-03   0.11669E-07   -1311.1    
c   3     A3         1.7757       0.41126E-01   0.30192E-05   -4.9699    
c   4     A4         1.8653       0.19808E-01   0.15147E-05   -9.7134    
c   5     A5       -0.40858E-10   0.94613E-11   0.73328E-15  -0.72349E+10
c   6     A6         22.715       0.18991       0.14108E-04   0.39347    
c   7     A7        0.24000       0.86567E-02   0.61397E-06    10.778    
c   8     A8        0.11370       0.82822E-02   0.56991E-06   -18.889    
c   9     A9         2.0616       0.15468E-01   0.11112E-05   -10.362    
c  10     A10        1.3202       0.63173E-02   0.62954E-06   -25.527    
c                               ERR DEF=  11.5  
c
	b1=99.155d0
	b2=0.47865d-3
	b3=1.7757d0
	b4=1.8653d0
	b5=-0.40858d-10
	b6=22.715d0
	b7=0.24000d0
	b8=0.11370d0
	b9=2.0616d0
	b10=1.3202d0
c
c--------------------------------------------	
c       
	npts=10000!2000
c	deltas=100.d0
c	ssroot0=1.d0
c
	do iesse=0,npts
         if (iesse.le.399) then
          ssroot0=3.d0
          deltas=0.25d0
         else if (400.le.iesse .and. iesse.le.799) then
          ssroot0=-900.d0
          deltas=2.5d0
         else if (800.le.iesse .and. iesse.le.1199) then
          ssroot0=-18900.d0
          deltas=25.d0
         else if (1200.le.iesse .and. iesse.le.1400) then
          ssroot0=-288900.d0
          deltas=250.d0
         end if
        ssroot=ssroot0+iesse*deltas
        ss=ssroot*ssroot
c        
c-----------------------------------------------------------------------
c Cálculo de ReSigQCD(s)
c-----------------------------------------------------------------------
c
        Xp=dlog(dlog(ss))
c 
        w1=b1
        w2=b2*dexp(b3*(Xp**(1.01d0*b4)))
        w3=b5*dexp(b6*(Xp**(1.05d0*b7)))
        w4=b8*dexp(b9*(Xp**(1.09d0*b10)))
        ww=w1+w2+w3+w4
c
	ReSigQCDGeV2m100_Der=ww-100.d0
c        
c-----------------------------------------------------------------------
c Cálculo de ImSigQCD(s) via RDD
c-----------------------------------------------------------------------
c
c 
	d1=b2*dexp(b3*Xp**(1.01d0*b4))*b3*1.01d0*b4*
     & Xp**(1.01d0*b4-1.d0)*1.d0/(ss*dlog(ss))
c     
	d2=b5*dexp(b6*Xp**(1.05d0*b7))*b6*1.05d0*b7*
     & Xp**(1.05d0*b7-1.d0)*1.d0/(ss*dlog(ss))
c     
	d3=b8*dexp(b9*Xp**(1.09d0*b10))*b9*1.09d0*b10*
     & Xp**(1.09d0*b10-1.d0)*1.d0/(ss*dlog(ss))     
c
	ImSigQCDGeV2m100_Der=-(pi*ss/2.d0)*(d1+d2+d3)
c
	write(*,*)ssroot,ReSigQCDGeV2m100_Der,ImSigQCDGeV2m100_Der
c
	write(1,*)ssroot,ReSigQCDGeV2m100_Der
	write(2,*)ssroot,-ImSigQCDGeV2m100_Der
c
       end do !ss
       stop
       end
