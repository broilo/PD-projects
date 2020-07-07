	program gerSigQCD_CTEQ6Llo
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
	open(unit=1,file='gerReSigQCDder_CTEQ6Llo.dat')
	open(unit=2,file='gerImSigQCDder_CTEQ6Llo.dat')
c
c--------------------------------------------
c Definição das constantes:
c      
	pi=4.d0*datan(1.d0)
        mbfactor=0.389379d0
c
c FCN=  0.4646056     FROM MIGRAD    STATUS=CONVERGED    128 CALLS     6364 TOTAL
c                     EDM=  0.14E-03    STRATEGY= 1      ERR MATRIX NOT POS-DEF
c
c  EXT PARAMETER                APPROXIMATE        STEP         FIRST   
c  NO.   NAME        VALUE          ERROR          SIZE      DERIVATIVE 
c   1     A1         97.819        2.0494       0.51227E-04   0.16270E-02
c   2     A2        0.26229E-01   0.27068E-02   0.13845E-06    106.92    
c   3     A3         1.5642       0.13473E-01   0.74585E-06    19.284    
c   4     A4         1.8017       0.73426E-02   0.85912E-06    32.176    
c   5     A5       -0.30293E-08   0.59900E-09   0.33822E-13   0.18721E+09
c   6     A6         19.854       0.16671       0.94672E-05  -0.77968    
c   7     A7        0.22379       0.92146E-02   0.47361E-06   -21.135    
c   8     A8        0.62907       0.81892E-01   0.39964E-05    2.3335    
c   9     A9         1.5934       0.29808E-01   0.15517E-05    7.7118    
c  10     A10        1.3127       0.15895E-01   0.88864E-06    15.137    
c                               ERR DEF=  11.5 
c
	b1=97.819d0
	b2=0.26229d-1
	b3=1.5642d0
	b4=1.8017d0
	b5=-0.30293d-8
	b6= 19.854d0
	b7=0.22379d0
	b8=0.62907d0
	b9=1.5934d0
	b10=1.3127d0
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
