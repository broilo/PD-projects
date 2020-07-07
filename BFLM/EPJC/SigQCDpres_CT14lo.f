	program gerSigQCD_CT14lo
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
	open(unit=1,file='gerReSigQCDpresc_CT14lo.dat')
	open(unit=2,file='gerImSigQCDpresc_CT14lo.dat')
c
c--------------------------------------------
c Definição das constantes:
c      
	pi=4.d0*datan(1.d0)
        mbfactor=0.389379d0
        imagi=dcmplx(0.d0,1.d0)
        ii=imagi
c
ccc 10mil Calls
c
c FCN=   10.79836     FROM MIGRAD    STATUS=CONVERGED    128 CALLS     3149 TOTAL
c                     EDM=  0.12E-04    STRATEGY= 1      ERR MATRIX NOT POS-DEF
c
c  EXT PARAMETER                APPROXIMATE        STEP         FIRST   
c  NO.   NAME        VALUE          ERROR          SIZE      DERIVATIVE 
c   1     A1         100.22        1.4924       0.63580E-04  -0.62420E-01
c   2     A2        0.43425E-01   0.54975E-02   0.39069E-06   -5.6496    
c   3     A3         1.2736       0.14536E-01   0.11133E-05   -2.4184    
c   4     A4         1.9189       0.98502E-02   0.91498E-06   -3.8546    
c   5     A5        0.12198E-07   0.21820E-08   0.16978E-12  -0.31305E+07
c   6     A6         14.050       0.11146       0.80548E-05  -0.66261E-01
c   7     A7        0.50348       0.75773E-02   0.52459E-06   -1.0744    
c   8     A8         3699.4        4543.5       0.40065      -0.25060E-06
c   9     A9        -80.280        23.950       0.20371E-02   0.53892E-03
c  10     A10       -2.6318       0.30144       0.23009E-04  -0.65531E-01
c                               ERR DEF=  11.5 
c
	b1=100.22d0
	b2=0.43425d-1
	b3=1.2736d0
	b4=1.9189d0
	b5=0.12198d-7
	b6=14.050d0
	b7=0.50348d0
	b8=3699.4d0
	b9=-80.280d0
	b10=-2.6318d0
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

