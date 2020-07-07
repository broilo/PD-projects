	program gerSigQCD_MMHTlo
c
c Gera Sig_QCD
c
c Obs1: Eq.9 do PRD92,074039(2015)
c
c Obs2: Calcula a também a parte imaginária
c
c Obs3: Usa a prescrição s -> -is
c
c-----------------------------------------------------------------------
c Importante: Foi minimizado com uma parametrização complexa cujo a
c parte real é exatamente igual aos pontos gerados
c-----------------------------------------------------------------------
c
	implicit none
	double precision b1,b2,b3,b4,b5,b6,b7,b8,b9,b10
        double precision ssroot0,ssroot,deltas,ss
        double precision pi,mbfactor
        double precision ReSigQCDGeV2p100_Presc,ReSigQCDGeV2m100_Presc
        double precision ImSigQCDGeV2p100_Presc,ImSigQCDGeV2m100_Presc
        double complex iss,imagi,ii,iXp
        double complex iww,iw1,iw2,iw3,iw4
        integer npts,iesse
c  
	open(unit=1,file='gerReSigQCDpresc_MMHTlo.dat')
	open(unit=2,file='gerImSigQCDpresc_MMHTlo.dat')
c
c--------------------------------------------
c Definição das constantes:
c      
	pi=4.d0*datan(1.d0)
        mbfactor=0.389379d0
        imagi=dcmplx(0.d0,1.d0)
        ii=imagi
c
c FCN=  0.7422312     FROM MIGRAD    STATUS=CONVERGED     42 CALLS     6734 TOTAL
c                     EDM=  0.33E-06  STRATEGY=1  ERROR MATRIX UNCERTAINTY=  1.3%
c
c  EXT PARAMETER                                   STEP         FIRST   
c  NO.   NAME        VALUE          ERROR          SIZE      DERIVATIVE 
c   1     A1         95.284        6.5128       0.74013E-05   0.21794E-01
c   2     A2        0.37224       0.14810      -0.12498E-06   0.54876    
c   3     A3        0.59990       0.24360E-02  -0.10180E-06    3.7301    
c   4     A4         2.4959       0.28733E-01   0.72729E-07    1.0398    
c   5     A5       -0.25543E-05   0.32482E-06   0.48345E-11   0.10403E+06
c   6     A6         14.281       0.11520       0.16623E-05  -0.34380    
c   7     A7        0.28071       0.61193E-01   0.10582E-05   -5.3795    
c   8     A8        0.90856E-01   0.22149E-01   0.84415E-07    5.7174    
c   9     A9         4.2904       0.16081       0.23726E-05   0.80930    
c  10     A10       0.67342       0.59826E-01   0.15390E-05    2.3202    
c                               ERR DEF=  11.5    
c
	b1=95.284d0
	b2=0.37224d0
	b3=0.59990d0
	b4=2.4959d0
	b5=-0.25543d-5
	b6=14.281d0
	b7=0.28071d0
	b8=0.90856d-1
	b9=4.2904d0
	b10=0.67342d0
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
        iss=-ii*ssroot*ssroot
        ss=ssroot*ssroot
c        
c-----------------------------------------------------------------------
c Cálculo de SigQCD(s) via prescrição s -> -is
c-----------------------------------------------------------------------
c
        iXp=cdlog(cdlog(iss))
c 
        iw1=b1
        iw2=b2*cdexp(b3*(iXp**(1.01d0*b4)))
        iw3=b5*cdexp(b6*(iXp**(1.05d0*b7)))
        iw4=b8*cdexp(b9*(iXp**(1.09d0*b10)))
        iww=iw1+iw2+iw3+iw4
c
	ReSigQCDGeV2p100_Presc=dreal(iww)
	ReSigQCDGeV2m100_Presc=dreal(iww-100.d0)
	ImSigQCDGeV2p100_Presc=dimag(iww)
	ImSigQCDGeV2m100_Presc=dimag(iww-100.d0)   
c        
	write(*,*)ssroot,ReSigQCDGeV2m100_Presc,ImSigQCDGeV2m100_Presc
c
	write(1,*)ssroot,ReSigQCDGeV2m100_Presc
	write(2,*)ssroot,-ImSigQCDGeV2m100_Presc
c
c
       end do !ss
       stop
       end

