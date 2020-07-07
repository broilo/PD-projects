	program gerSigQCD_CTEQ6Llo
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
	open(unit=1,file='gerReSigQCDpresc_CTEQ6Llo.dat')
	open(unit=2,file='gerImSigQCDpresc_CTEQ6Llo.dat')
c
c--------------------------------------------
c Definição das constantes:
c      
	pi=4.d0*datan(1.d0)
        mbfactor=0.389379d0
        imagi=dcmplx(0.d0,1.d0)
        ii=imagi
c
c FCN=  0.3727160     FROM MIGRAD    STATUS=CONVERGED    128 CALLS     3267 TOTAL
c                     EDM=  0.99E-05    STRATEGY= 1      ERR MATRIX NOT POS-DEF
c
c  EXT PARAMETER                APPROXIMATE        STEP         FIRST   
c  NO.   NAME        VALUE          ERROR          SIZE      DERIVATIVE 
c   1     A1         97.005        1.8767       0.52957E-04   0.32854E-01
c   2     A2        0.28009E-01   0.26942E-02   0.14199E-06    11.206    
c   3     A3         1.6986       0.14280E-01   0.80994E-06    1.8379    
c   4     A4         1.7362       0.72495E-02   0.82788E-06    3.0838    
c   5     A5       -0.14947E-05   0.13699E-06   0.76678E-11  -0.12058E+06
c   6     A6         14.140       0.72785E-01   0.67422E-05   0.24483    
c   7     A7        0.31972       0.55084E-02   0.27450E-06    3.4023    
c   8     A8        0.83551E-01   0.66240E-02   0.31975E-06   -2.8321    
c   9     A9         3.8127       0.33703E-01   0.18180E-05  -0.48168    
c  10     A10       0.80990       0.71546E-02   0.41359E-06   -1.5070    
c                               ERR DEF=  11.5 
c
	b1=97.005d0
	b2=0.28009d-1
	b3=1.6986d0
	b4=1.7362d0
	b5=-0.14947d-5
	b6=14.140d0
	b7=0.31972d0
	b8=0.83551d-1
	b9=3.8127d0
	b10=0.80990d0
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

