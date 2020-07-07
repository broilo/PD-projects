	program gerSigQCD_MMHTlo
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
	open(unit=1,file='gerReSigQCDder_MMHTlo.dat')
	open(unit=2,file='gerImSigQCDder_MMHTlo.dat')
c
c--------------------------------------------
c Definição das constantes:
c      
	pi=4.d0*datan(1.d0)
        mbfactor=0.389379d0
c
c FCN=  0.1021195     FROM MIGRAD    STATUS=CONVERGED    128 CALLS     1032 TOTAL
c                     EDM=  0.54E-04    STRATEGY= 1      ERR MATRIX NOT POS-DEF
c
c  EXT PARAMETER                APPROXIMATE        STEP         FIRST   
c  NO.   NAME        VALUE          ERROR          SIZE      DERIVATIVE 
c   1     A1         99.934        7.4205       0.47150E-04  -0.28936    
c   2     A2         4.0982       0.22117       0.14849E-04  -0.23073    
c   3     A3        0.18627       0.18949E-02   0.10239E-06   -17.663    
c   4     A4         3.2762       0.83259E-02   0.15622E-05   -3.5016    
c   5     A5       -0.14961E-07   0.26586E-07   0.12697E-11  -0.85203E+07
c   6     A6         18.019        1.5692       0.74765E-04   0.14322    
c   7     A7        0.15125       0.84351E-01   0.47439E-05    2.0335    
c   8     A8        -1.2886        2.7312       0.13654E-03  -0.10045    
c   9     A9        0.90585        1.1746       0.55213E-04   0.22753    
c  10     A10       0.82852        1.4612       0.71107E-04   0.14314    
c                               ERR DEF=  11.5    
c
	b1=99.934d0
	b2=4.0982d0
	b3=0.18627d0
	b4=3.2762d0
	b5=-0.14961d-7
	b6=18.019d0
	b7=0.15125d0
	b8=-1.2886d0
	b9=0.90585d0
	b10=0.82852d0
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
