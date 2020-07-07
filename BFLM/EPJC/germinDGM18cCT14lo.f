	program germinDGM18cCT14lo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c O FATOR DE FORMA NÃO ESTÁ SENDO CONSIDERADO NA RDD!!!
c
c ImXsh é obtida via s->-is
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Gera St e Rh
c
c Obs1: Os fatores de forma K3(x) estão sendo obtidos via integração
c por Regra de Simpson modificada
c
c Obs2: Utiliza os novos conjuntos de dados do LHC
c
c Obs3: TOTEM + ATLAS
c
c Obs4: Usa uma parametrização para SigQCD com 10 parâmetros livres
c
c Obs5: PDF: CT14lo
c
c Obs6: Usa a estrutura do Durand & Pi
c
c Obs7: SigQCD(s) foi parametrizada por uma função real
c
c Obs8: A parte soft da eiconal é composta por termos Xqq e Xqg, conforme
c o DGM antigo, só que com um log2(s/s0).
c
c Obs9: mu^{-}=0.5 está fixo.
c
c--------------------------------------------
c Importante: s_{0}=25 GeV^{2} !4*m^{2}_{p} GeV^{2}
c             mg=0.4d0 !GeV -> mg=500+/-200 MeV for lmb=300 MeV
c (MSTW) lbd=0.318920 GeV lambda equivale a alfa_LO(Mz^2)=0.13939~0.139 fazendo o match; m_b=4.18 GeV
c (MMHT) lbd=0.267471 GeV lambda equivale a alfa_LO(Mz^2)=0.13499~0.135 fazendo o match; m_b=4.18 GeV
c (CTEQ6L) lbd=0.3260 GeV lambda equivale a alfa_NLO(Mz^2)=0.118 para 4 sabores
c (CT14)   lbd=0.3260 GeV lambda equivale a alfa_NLO(Mz^2)=0.118 para 4 sabores
c--------------------------------------------
c
c------------------------------------
c chi2/gl=FCN/(Dados-Livres)=195.7392/(173-7)=1.17915
c------------------------------------
c
	implicit none
	double precision dgauss,bmin,bmax,epsabs
	double precision mbfactor,pi,pi2,pi4,aa0,mp,mp2
	double precision nf,lbd,lbd2,Bz,alfa0,alfa02,SIGMA
	double precision y1,y2,y3,y4,y5,y6,y7,y8,y9
	double precision b1,b2,b3,b4,b5,b6,b7,b8,b9,b10
	double precision argImFpp,ImFpp,argImFpa,ImFpa
	double precision argReFpp,ReFpp,argReFpa,ReFpa
	double precision Stpp,Stpa,Rhpp,Rhpa
	double precision ss,ss0,ssroot,ssroot0,deltas,passo,s
	double complex ii,iss,is
	integer iesse,npts
	external argImFpp,argImFpa,argReFpp,argReFpa
	common/cte/pi,pi2,pi4,aa0
	common/imag/ii
	common/par/y1,y2,y3,y4,y5,y6,y7,y8,y9
	common/sig/b1,b2,b3,b4,b5,b6,b7,b8,b9,b10
	common/ener/ss,iss,ss0,s,is
c
	open(unit=1,file='Stpp_minDGM18cCT14lo.dat')
	open(unit=2,file='Stpa_minDGM18cCT14lo.dat')
	open(unit=3,file='Rhpp_minDGM18cCT14lo.dat')
	open(unit=4,file='Rhpa_minDGM18cCT14lo.dat')
c
c--------------------------------------------
c Definição das constantes comuns:
c
	mbfactor=0.389379d0
	pi=4.d0*datan(1.d0)
	pi2=pi/2.d0
	pi4=pi/4.d0
	ii=dcmplx(0.d0,1.d0)
cc	mp=0.938272046d0!massa do próton em GeV
cc	mp2=mp*mp
	ss0=25.d0!4.d0*mp2
	aa0=1.d0
c
c--------------------------------------------	
c Parâmetros ajustados:
c
c minDGM18cCT14lo.min
c
ccc 1 sigma:
c
c MIGRAD MINIMIZATION HAS CONVERGED.
c
c FCN=   195.7392     FROM MIGRAD    STATUS=CONVERGED     21 CALLS      469 TOTAL
c                     EDM=  0.14E-06  STRATEGY=1  ERROR MATRIX UNCERTAINTY=  1.0%
c
c  EXT PARAMETER                                   STEP         FIRST   
c  NO.   NAME        VALUE          ERROR          SIZE      DERIVATIVE 
c   1      A1        2.3859       0.55094        0.0000      -0.18735E-02
c   2      A2       0.58388E-01   0.35970E-01    0.0000       0.22637E-01
c   3      A3        123.49        2.5460        0.0000       0.74028E-03
c   4      A4        42.372        8.0563        0.0000       0.13065E-03
c   5      A5       0.73205       0.14260        0.0000       0.61728E-02
c   6      A6       0.89990       0.17868        0.0000      -0.84634E-03
c   7      A7        24.186        1.4300        0.0000      -0.80169E-04
c                               ERR DEF=  8.18    
c
c EXTERNAL ERROR MATRIX.    NDIM=  50    NPAR=  7    ERR DEF=  8.18    
c  0.304E+00
c  0.198E-01 0.129E-02
c -0.690E+00-0.413E-01 0.648E+01
c  0.227E+01 0.138E+00-0.154E+02 0.649E+02
c  0.568E-01 0.354E-02-0.234E+00 0.997E+00 0.203E-01
c  0.429E-03 0.364E-04 0.208E-01 0.311E-01 0.804E-03 0.333E-03
c -0.238E-01-0.139E-02 0.599E+00-0.949E+00-0.872E-02 0.149E-02 0.204E+01
c
	y1=2.3859d0   !nu1
	y2=0.58388d-1 !nu2
	y3=123.49d0   !A
	y4=42.372d0   !B
	y5=0.73205d0  !C
	y6=0.89990d0  !mu^{+}_{soft}
	y7=24.186d0   !D
	y8=1.d0       !K_{f}  ### fixado em 1.0 ###
	y9=0.5d0      !mu^{-} ### Pode ser fixado em 0.5 ###	
c	
ccc 2 sigmas:
c

c

c
c Parâmetros ajustados na parametrização complexa de SigQCD:
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
	epsabs=1.d-8
	bmin=0.d0
	bmax=30.d0
	passo=5.d0
	npts=149!0
c
	do iesse=0,npts
	ssroot0=8.d0!1960.d0!14000.d0!8.d0
	deltas=iesse*iesse*passo/100.d0
	ssroot=ssroot0+iesse*deltas
	ss=ssroot*ssroot
	iss=-ii*ss
	s=ss/ss0
	is=-ii*s
c
	ImFpp=dgauss(argImFpp,bmin,bmax,epsabs)
	Stpp=4.d0*pi*ImFpp*mbfactor
c
	ImFpa=dgauss(argImFpa,bmin,bmax,epsabs)
	Stpa=4.d0*pi*ImFpa*mbfactor
c
	ReFpp=dgauss(argReFpp,bmin,bmax,epsabs)
	Rhpp=-ReFpp/ImFpp
c
	ReFpa=dgauss(argReFpa,bmin,bmax,epsabs)
	Rhpa=-ReFpa/ImFpa
c
	write(*,*)ssroot,Stpp,Stpa,Rhpp,Rhpa
	write(1,*)ssroot,Stpp
	write(2,*)ssroot,Stpa
	write(3,*)ssroot,Rhpp
	write(4,*)ssroot,Rhpa	
c
	end do! ss
	stop
	end
c
c
c
	function argImFpp(bb)
	implicit none
	double precision argImFpp,bb
	double precision Wdip
	double complex SigQCD
	double precision nu
	double complex chie,chio,chipp
	double complex Xsoft,Xsh
	double precision pi,pi2,pi4,aa0
	double precision y1,y2,y3,y4,y5,y6,y7,y8,y9
	double precision ss,ss0,s
	double complex ii,iss,is
	external Wdip,SigQCD
	common/cte/pi,pi2,pi4,aa0
	common/imag/ii
	common/par/y1,y2,y3,y4,y5,y6,y7,y8,y9
	common/ener/ss,iss,ss0,s,is
c
	nu=y1-y2*dlog(s)
c
	Xsoft=0.5d0*(y3+y4*cdexp(ii*pi4)/dsqrt(s)
     & +y5*(dlog(s)-ii*pi2)**2.d0)*Wdip(bb,y6)
c
	Xsh=0.5d0*(y8*Wdip(bb,nu)*SigQCD(aa0))
c
	chie=Xsoft+Xsh
c
	chio=0.5d0*y7*cdexp(-ii*pi4)/dsqrt(s)*Wdip(bb,y9)
c
	chipp=(chie-chio) !pp
c
	argImFpp=bb*(1.d0-dexp(-dreal(chipp))*dcos(dimag(chipp)))
c
	return
	end
c
c
c
	function argImFpa(bb)
	implicit none
	double precision argImFpa,bb
	double precision Wdip
	double complex SigQCD
	double precision nu
	double complex chie,chio,chipa
	double complex Xsoft,Xsh
	double precision pi,pi2,pi4,aa0
	double precision y1,y2,y3,y4,y5,y6,y7,y8,y9
	double precision ss,ss0,s
	double complex ii,iss,is
	external Wdip,SigQCD
	common/cte/pi,pi2,pi4,aa0
	common/imag/ii
	common/par/y1,y2,y3,y4,y5,y6,y7,y8,y9
	common/ener/ss,iss,ss0,s,is
c
	nu=y1-y2*dlog(s)
c
	Xsoft=0.5d0*(y3+y4*cdexp(ii*pi4)/dsqrt(s)
     & +y5*(dlog(s)-ii*pi2)**2.d0)*Wdip(bb,y6)
c
	Xsh=0.5d0*(y8*Wdip(bb,nu)*SigQCD(aa0))
c
	chie=Xsoft+Xsh
c
	chio=0.5d0*y7*cdexp(-ii*pi4)/dsqrt(s)*Wdip(bb,y9)
c
	chipa=(chie+chio) !pa
c
	argImFpa=bb*(1.d0-dexp(-dreal(chipa))*dcos(dimag(chipa)))
c
	return
	end
c
c
c
	function argReFpp(bb)
	implicit none
	double precision argReFpp,bb
	double precision Wdip
	double complex SigQCD
	double precision nu
	double complex chie,chio,chipp
	double complex Xsoft,Xsh
	double precision pi,pi2,pi4,aa0
	double precision y1,y2,y3,y4,y5,y6,y7,y8,y9
	double precision ss,ss0,s
	double complex ii,iss,is
	external Wdip,SigQCD
	common/cte/pi,pi2,pi4,aa0
	common/imag/ii
	common/par/y1,y2,y3,y4,y5,y6,y7,y8,y9
	common/ener/ss,iss,ss0,s,is
c
	nu=y1-y2*dlog(s)
c
	Xsoft=0.5d0*(y3+y4*cdexp(ii*pi4)/dsqrt(s)
     & +y5*(dlog(s)-ii*pi2)**2.d0)*Wdip(bb,y6)
c
	Xsh=0.5d0*(y8*Wdip(bb,nu)*SigQCD(aa0))
c
	chie=Xsoft+Xsh
c
	chio=0.5d0*y7*cdexp(-ii*pi4)/dsqrt(s)*Wdip(bb,y9)
c
	chipp=(chie-chio) !pp
c
	argReFpp=bb*(dexp(-dreal(chipp))*dsin(dimag(chipp)))
c
	return
	end
c
c
c
	function argReFpa(bb)
	implicit none
	double precision argReFpa,bb
	double precision Wdip
	double complex SigQCD
	double precision nu
	double complex chie,chio,chipa
	double complex Xsoft,Xsh
	double precision pi,pi2,pi4,aa0
	double precision y1,y2,y3,y4,y5,y6,y7,y8,y9
	double precision ss,ss0,s
	double complex ii,iss,is
	external Wdip,SigQCD
	common/cte/pi,pi2,pi4,aa0
	common/imag/ii
	common/par/y1,y2,y3,y4,y5,y6,y7,y8,y9
	common/ener/ss,iss,ss0,s,is
c
	nu=y1-y2*dlog(s)
c
	Xsoft=0.5d0*(y3+y4*cdexp(ii*pi4)/dsqrt(s)
     & +y5*(dlog(s)-ii*pi2)**2.d0)*Wdip(bb,y6)
c
	Xsh=0.5d0*(y8*Wdip(bb,nu)*SigQCD(aa0))
c
	chie=Xsoft+Xsh
c
	chio=0.5d0*y7*cdexp(-ii*pi4)/dsqrt(s)*Wdip(bb,y9)
c
	chipa=(chie+chio) !pa
c
	argReFpa=bb*(dexp(-dreal(chipa))*dsin(dimag(chipa)))
c
	return
	end	
c
c
c
c----------------------------------------------------------------------
c Calcula a seção de choque de QCD com a PDF: CT14lo ij -> ij
c----------------------------------------------------------------------
c
c
c
	function SigQCD(aa)
	implicit none
	double complex SigQCD
	double precision aa
	double precision b1,b2,b3,b4,b5,b6,b7,b8,b9,b10
	double complex ww,w0,w1,w2,w3,w4
	double precision ss,ss0,ssroot,s
	double complex iss,is,Xp
	common/sig/b1,b2,b3,b4,b5,b6,b7,b8,b9,b10
	common/ener/ss,iss,ss0,s,is
c
	Xp=cdlog(cdlog(iss))
c 
	w1=b1
	w2=b2*cdexp(b3*(Xp**(1.01d0*b4)))
	w3=b5*cdexp(b6*(Xp**(1.05d0*b7)))
	w4=b8*cdexp(b9*(Xp**(1.09d0*b10)))
c
	w0=w1+w2+w3+w4
	ww=(w0-100.d0)
	SigQCD=aa*ww
c
	return
	end
c
c
c
c-----------------------------------------------------------------------
c Fator de forma do tipo dipolo
c-----------------------------------------------------------------------
c
c
c
	function Wdip(b,mu)
	implicit none
	double precision Wdip,b,mu
	double precision zz,zz3,xx1,xx2,xx3,nn3
	double precision ff,H,U
	double precision BesselK3w1,BesselK3w2,BesselK3w3
	double precision funcK3,IntK3,BesselK3
	double precision pi,pi2,pi4,aa0
	integer N,Ndiv
	common/cte/pi,pi2,pi4,aa0
c
	zz=mu*b
	zz3=zz**3.d0
c
	nn3=3.0
	ff=60.d0 !limite superior de integração
	Ndiv=200
c
	do N=1,Ndiv
	 H=ff/(6.d0*Ndiv)
	 U=ff/Ndiv
	 xx1=(N-1.d0)*U !xx's são as variáveis de integração
	 xx2=((N-1.d0)*U+N*U)/2.d0
	 xx3=N*U
c
	 BesselK3w1=dexp(-zz*dcosh(xx1))*dcosh(nn3*xx1)
	 BesselK3w2=dexp(-zz*dcosh(xx2))*dcosh(nn3*xx2)
	 BesselK3w3=dexp(-zz*dcosh(xx3))*dcosh(nn3*xx3)
c
	 funcK3=H*(BesselK3w1+4.d0*BesselK3w2+BesselK3w3)
	 IntK3=IntK3+funcK3
	end do !N
c
	BesselK3=IntK3
	Wdip=(mu**2.d0)*zz3*BesselK3/(96.d0*pi)
	IntK3=0.d0
c
	return
	end

