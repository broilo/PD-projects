/* PROGRAM: BSG20 based on Lipari-Lusignoli original model

   ---> Ref. Original model: PHYSICAL REVIEW D 80, 074014 (2009)
        SciHub link: https://sci-hub.tw/10.1103/PhysRevD.80.074014
        
   Modified version: SigEik=SigSoft+SigPQCD

------------------------------------------------------------------------
   *** Main characteristics ***

	> SigEik = SigSoft + SigPQCD
	> SigSoft = A1*s^{-del1} +/- A2*s^{-del2} + Sig0
	> SigPQCD = B*SigPQCD
	
	> Calculates Total, Elastic, Diff and Eik xsec pp and pbp

   *** Observations ***

 Obs1: TOTEM + ATLAS dataset
 Obs2: Real parametrization for SigQCD with 10 free-parameters
 Obs3: PDF: MMHT_member_39
 Obs4: Free: w, Sig0, A1, A2, del1, del2
 Obs5: B is fixed in 1.0

   *** Important *** 
   (See: https://arxiv.org/abs/1908.01040 Chapter 8, sections 8.4 and 8.8)
   
   mg=0.4d0 !GeV -> mg=500+/-200 MeV for lmb=300 MeV
 
 --> (MSTW) lbd=0.318920 GeV lambda equivalent to alpha_LO(Mz^2)=0.13939~0.139 macthing-prescrition squeme; m_b=4.18 GeV
 
 --> (MMHT) lbd=0.267471 GeV lambda equivalent to alpha_LO(Mz^2)=0.13499~0.135 macthing-prescrition squeme; m_b=4.18 GeV
 
 --> (CTEQ6L) lbd=0.3260 GeV lambda equivalent to alpha_NLO(Mz^2)=0.118 for 4 flavors
 
 --> (CT14)   lbd=0.3260 GeV lambda equivalent to alpha_NLO(Mz^2)=0.118 for 4 flavors

   *** How to run the code ***
 
 root -l -q /usr/lib/x86_64-linux-gnu/libgsl.so /usr/lib/x86_64-linux-gnu/libgslcblas.so BSG20orig_w002b_MMHT_member_39.C | cat > BSG20orig_w002b_MMHT_member_39.min
 
  MMHT_member_39 -> See: https://lhapdf.hepforge.org/index.html
   ---------------------------------------------------------------------
   Codes by: M.Broilo*, P.V.R.G.Silva, V.P.B.Gonçalves
   mateus.broilo@ufrgs.br | mateus.broilo90@gmail.com
   High and Medium Energy Group
   Grupo de Altas e Medias Energias (GAME)
   Universidade Federal de Pelotas, Pelotas - RS (Brasil)
   ---------------------------------------------------------------------
   
   Creation: 18/dez/2019 (Pelotas, RS, Brazil)
   Last update: 20/mar/2020 (Pelotas, RS, Brazil)
   
*/
#include <chrono> 
#include "gsl/gsl_sf.h"
#include "gsl/gsl_sf_bessel.h"
#include "gsl/gsl_sf_hyperg.h"
using namespace chrono;
using namespace ROOT;
using namespace TMath;
/*
   *** Some convertions ***
   GeV^{-1} = 0.1975 fm ~0.19 fm
   GeV^{-2} = 0.389379 mb
   mb -> 0.1 fm2
   fm2 = 10 mb
*/
const double mbfactor=0.389379; //mbGeV^{2}
const double mb_to_fm2=0.1;
const double fm2_to_mb=10.;
//----------------------------------------------------------------------
//constant parameters throughout the model
const double s0=25.;
const double r0=0.19;//fm
//----------------------------------------------------------------------
//number of fit parameters and data points 
const int numpar=6;
const int npSigTotPBP=59;
const int npSigTotPP=118;
const int npapfit=npSigTotPBP;
const int nppfit=npSigTotPP;
const int npMin=npapfit+nppfit;

const int npSigElPBP=141;
const int npSigElPP=155;

const int npSigDiffpp=5.;
const int npSigDiffpbp=16.;

// Plots
const int npSigTotPPplot=104;
const int npSigTotPPTOTEM=12;
const int npSigTotPPATLAS=2;

const int npSigElPPplot=147;
const int npSigElPPTOTEM=6;
const int npSigElPPATLAS=2;

const int npSigDiff1=1;
const int npSigDiff2=1;
const int npSigDiff3=3;
const int npSigDiff4=8;
const int npSigDiff5=3;
const int npSigDiff6=2;
const int npSigDiff7=2;
const int npSigDiff8=1;
//----------------------------------------------------------------------
//parameters to control energy and b range (for plots and integrations)
const double Wmin=5.;
const double Wmax=1.e6;
const double bmin=0.;
const double bmax=30.;
//----------------------------------------------------------------------
//QCD cross section - analytical parametrization for MMHT_member_39
//----------------------------------------------------------------------
const double c1=1.01;
const double c2=1.05;
const double c3=1.09;
const double bkg=100.;
// SigQCD(GG e QG) parameters
const double b1=1.03566e+02;
const double b2=2.18981e-01;
const double b3=3.63494e-01;
const double b4=2.87487e+00;
const double b5=5.49698e-02;
const double b6=1.34048e+00;
const double b7=1.75436e+00;
const double b8=-2.43193e+00;
const double b9=4.92022e-01;
const double b10=8.42211e-01; 
/*
 NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
   1  b1           1.03566e+02   2.31937e+00   2.29081e-04  -8.87488e-03
   2  b2           2.18981e-01   2.73959e-03   3.11883e-06  -3.57109e+00
   3  b3           3.63494e-01   7.91437e-04   7.56400e-07  -7.63553e+00
   4  b4           2.87487e+00   2.33079e-03   2.16897e-06  -1.25987e+00
   5  b5           5.49698e-02   7.71680e-03   6.26027e-06  -6.68544e+00
   6  b6           1.34048e+00   7.20983e-03  -6.77725e-06  -2.30007e+00
   7  b7           1.75436e+00   1.22936e-02  -1.09265e-05  -4.62734e+00
   8  b8          -2.43193e+00   8.59676e-01  -8.63908e-04   1.25987e-02
   9  b9           4.92022e-01   6.99392e-02   6.88148e-05   2.28629e-02
  10  b10          8.42211e-01   5.36634e-02   4.36627e-05   1.53402e-01
*/
//----------------------------------------------------------------------
/*
                              Analysis 
*/
//----------------------------------------------------------------------
// Sigma pQCD
//----------------------------------------------------------------------
double SigPQCD(double W)
{
   double s=W*W;
        
    double Y=Log(s);
    double X=Log(Y);

    double sig=(b1-bkg)+b2*Exp(b3*Power(X,c1*b4))+b5*Exp(b6*Power(X,c2*b7))+b8*Exp(b9*Power(X,c3*b10));
          
    return sig;
}//----------------------------------------------------------------------
//Integration kernels
//----------------------------------------------------------------------
//Forward observables - total, elastic and diff xsection
double KerForward(double *x,double *par)
{
    double b=x[0];      
    double ww=par[0]; 
    double Sig0=par[1];
    double A1=par[2];
    double del1=par[3];
    double A2=par[4];
    double del2=par[5];        
    double B=par[6];
    double reac=par[7];
    double obs=par[8];
    double W=par[9];        
        
    double s=W*W;
    double winv=1./ww;
    /*******************************************************/
    //Form factors
    /*******************************************************/
    double K3=gsl_sf_bessel_Kn(3,b/r0);
    double Wdip=Power(b,3)*K3/(96.*Pi()*Power(r0,5));
    /*******************************************************/
    //SigPQCD
    /*******************************************************/
    double SigMinJet=B*SigPQCD(W)*mbfactor;
    /*******************************************************/
    //SigSoft
    /*******************************************************/
    double SigR1=A1*Power(s/s0,-del1);
    double SigR2=A2*Power(s/s0,-del2);
    //pp
    double SigSoftPP=(Sig0+SigR1-SigR2);
    //pbp
    double SigSoftPBP=(Sig0+SigR1+SigR2);
    /*******************************************************/
    //SigEik
    /*******************************************************/
    //pp
    double SigEikPP=SigSoftPP+SigMinJet;
    /*******************************************************/
    //pbp
    double SigEikPBP=SigSoftPBP+SigMinJet;
    /*******************************************************/
    double SigEik;
    double KerSig;
    double nn,xx,XSecEik;
    
    	if(reac==1) //pp
    	{
    	/*******************************************************/
    	//Average number of interactions
    	/*******************************************************/
	nn=(SigEikPP*mb_to_fm2)*Wdip;
	xx=(nn*ww)/2.;
    	}
    	else if(reac==2) //pbp
    	{
	nn=(SigEikPBP*mb_to_fm2)*Wdip;
	xx=(nn*ww)/2.;
    	}
    	else
    	{
    	nn=0.;
    	xx=0.;
    	}
    	if(obs==1) //St
    	{
            KerSig=2.*Pi()*b*(2.-2.*Power((1.+xx),(-winv)))*fm2_to_mb; //in mb 
            return KerSig;
        }
        else if(obs==2) //SigEl
        {
            KerSig=2.*Pi()*b*Power((1.-Power((1.+xx),(-winv))),2.)*fm2_to_mb; //in mb 
            return KerSig;
        }
        else if(obs==3) //Diff
        {
            KerSig=2.*Pi()*b*(Power((1.+nn*ww),(-winv))-Power((1.+xx),-2.*winv))*fm2_to_mb; //in mb 
            return KerSig;
        }
        else
        {
            return 0.;
        }
}
// eikonal xsection
double KerThings(double *x,double *par)
{
    double W=x[0];      
    double ww=par[0]; 
    double Sig0=par[1];
    double A1=par[2];
    double del1=par[3];
    double A2=par[4];
    double del2=par[5];        
    double B=par[6];
    double reac=par[7];
    double obs=par[8];
        
    double s=W*W;
    double winv=1./ww;
    /*******************************************************/
    //SigPQCD
    /*******************************************************/
    double SigMinJet=B*SigPQCD(W)*mbfactor;
    /*******************************************************/
    //SigSoft
    /*******************************************************/
    double SigR1=A1*Power(s/s0,-del1);
    double SigR2=A2*Power(s/s0,-del2);
    //pp
    double SigSoftPP=(Sig0+SigR1-SigR2);
    //pbp
    double SigSoftPBP=(Sig0+SigR1+SigR2);
    /*******************************************************/
    //SigEik
    /*******************************************************/
    //pp
    double SigEikPP=SigSoftPP+SigMinJet;
    /*******************************************************/
    //pbp
    double SigEikPBP=SigSoftPBP+SigMinJet;
    /*******************************************************/
    double SigEik;
    double KerSig;
    double nn,xx;
    double XSecSoft,XSecHard,XSecEik;
    
    	if(reac==1) //pp
    	{
	XSecEik=SigEikPP;
	XSecSoft=SigSoftPP;
    	}
    	else if(reac==2) //pbp
    	{
	XSecEik=SigEikPBP;
	XSecSoft=SigSoftPBP;
    	}
    	else
    	{
    	    return 0.;
    	}
    	
    	XSecHard=SigMinJet;
    	
    	if(obs==1) //XSecEik
    	{
    	KerSig=XSecEik; //in mb
    	}
    	if(obs==2) //XSecSoft
    	{
    	KerSig=XSecSoft; //in mb
    	}
    	if(obs==3) //XSecHard
    	{
    	KerSig=XSecHard; //in mb
    	}
            return KerSig;
}
//----------------------------------------------------------------------
// Plotando a seção de choque total pp e pbp
//----------------------------------------------------------------------
void PlotTotXSec(double *Wcm,double *XSecTOTpp,double *XSecTOTpbp,int npfit)
{
    FILE *aq1,*aq2,*aq3,*aq4;
// Plotting SigTot data     
    aq1 = fopen("paw_Stpa.dat","r");
    aq2 = fopen("paw_Stpp.dat","r");
    aq3 = fopen("paw_StTOTEM.dat","r");
    aq4 = fopen("paw_StATLAS.dat","r");

    // pbp
    const int npStpa = npSigTotPBP;
    double W21[npStpa],Sig21Exp[npStpa],uW21[npStpa],uSig21Exp[npStpa];
   
    for(int i=0;i<npStpa;i++)     
        fscanf(aq1,"%lg %lg %lg %lg",&W21[i],&Sig21Exp[i],&uSig21Exp[i],&uW21[i]);
        
     // pp     
     const int npStpp = npSigTotPPplot;
     double W11[npStpp],Sig11Exp[npStpp],uW11[npStpp],uSig11Exp[npStpp];
     
     for(int i=0;i<npStpp;i++)   
        fscanf(aq2,"%lg %lg %lg %lg",&W11[i],&Sig11Exp[i],&uSig11Exp[i],&uW11[i]);
       
       //---> TOTEM data <---
       
     const int npStTOTEM = npSigTotPPTOTEM;
     double W11totem[npStTOTEM],Sig11ExpTOTEM[npStTOTEM],uW11totem[npStTOTEM],uSig11ExpTOTEM[npStTOTEM];  
       
     for(int i=0;i<npStTOTEM;i++)   
        fscanf(aq3,"%lg %lg %lg %lg",&W11totem[i],&Sig11ExpTOTEM[i],&uSig11ExpTOTEM[i],&uW11totem[i]);
        
       //---> ATLAS data <---
       
     const int npStATLAS = npSigTotPPATLAS;
     double W11atlas[npStATLAS],Sig11ExpATLAS[npStATLAS],uW11atlas[npStATLAS],uSig11ExpATLAS[npStATLAS];  
       
     for(int i=0;i<npStATLAS;i++)   
        fscanf(aq4,"%lg %lg %lg %lg",&W11atlas[i],&Sig11ExpATLAS[i],&uSig11ExpATLAS[i],&uW11atlas[i]);
      
     fclose(aq1);
     fclose(aq2);
     fclose(aq3);
     fclose(aq4);
     
      TCanvas *canv1 = new TCanvas("c1","Total Cross Section",200,10,900,600);
      canv1->SetFillColor(0);
      canv1->SetLogx();
       
       double x_min(5e00),x_max(3.e04);
       double y_min(3.0e01),y_max(1.4e02);

      TGraphErrors *gr0 = new TGraphErrors(npStpa,W21,Sig21Exp,uW21,uSig21Exp);
      gr0->SetMarkerStyle(24);
      gr0->SetMarkerColor(1);
      gr0->SetLineColor(1);
      gr0->SetMarkerSize(1.0);
      gr0->SetTitle();
      gr0->Draw("apz");    
           
      TGraphErrors *gr1 = new TGraphErrors(npStpp,W11,Sig11Exp,uW11,uSig11Exp);
      gr1->SetMarkerStyle(20);
      gr1->SetMarkerColor(1);
      gr1->SetLineColor(1);
      gr1->SetMarkerSize(1.0);
      gr1->Draw("pz");         
      
      TGraphErrors *gr2 = new TGraphErrors(npStTOTEM,W11totem,Sig11ExpTOTEM,uW11totem,uSig11ExpTOTEM);
      gr2->SetMarkerStyle(33);
      gr2->SetMarkerColor(2);
      gr2->SetLineColor(2);
      gr2->SetMarkerSize(1.3);
      gr2->Draw("pz");  
      
      TGraphErrors *gr3 = new TGraphErrors(npStATLAS,W11atlas,Sig11ExpATLAS,uW11atlas,uSig11ExpATLAS);
      gr3->SetMarkerStyle(21);
      gr3->SetMarkerColor(4);
      gr3->SetLineColor(4);
      gr3->SetMarkerSize(0.8);
      gr3->Draw("pz");
      
      TGraph *grSigTotPBARP = new TGraph(npfit,Wcm,XSecTOTpbp);
      grSigTotPBARP->SetLineColor(4);
      grSigTotPBARP->SetLineWidth(1);
      grSigTotPBARP->SetLineStyle(1);
      grSigTotPBARP->Draw("c");
      
      TGraph *grSigTotPP = new TGraph(npfit,Wcm,XSecTOTpp);
      grSigTotPP->SetLineColor(2);
      grSigTotPP->SetLineWidth(1);
      grSigTotPP->SetLineStyle(2);
      grSigTotPP->Draw("c"); 

      gr0->GetYaxis()->SetTitle("#bf{#sigma_{tot} [mb]}");
      gr0->GetYaxis()->SetRangeUser(y_min,y_max);
      gr0->GetYaxis()->SetTitleOffset(1.2);
      gr0->GetYaxis()->SetTitleFont(42);
      gr0->GetXaxis()->SetTitle("#bf{#sqrt{s} [GeV]}");
      gr0->GetXaxis()->SetTitleOffset(1.3);
      gr0->GetXaxis()->SetLimits(x_min,x_max);   
      gr0->GetXaxis()->SetTitleFont(42);
           
      TLegend *leg1 = new TLegend(0.65,0.15,0.89,0.42);
      leg1->AddEntry(gr0,"#bar{p}p data ","p");
      leg1->AddEntry(gr1,"pp data below LHC","p");
      leg1->AddEntry(gr2,"TOTEM","p");
      leg1->AddEntry(gr3,"ATLAS","p");
      leg1->AddEntry(grSigTotPBARP,"#bar{p}p","l");
      leg1->AddEntry(grSigTotPP,"pp","l");
      leg1->SetFillColor(0);
      leg1->SetShadowColor(0);
      leg1->SetBorderSize(0);
      leg1->SetTextSize(0.035);
      leg1->Draw();    

/* ---> desenhando um miniframe em um pad dentro do canvas canv1 <--- */

    TPad *pad1 = new TPad("pad1","This is pad1",0.15,0.51,0.55,0.89);
    pad1->SetFillColor(0);
    pad1->SetLogx();
    pad1->Draw();     
    pad1->cd();
    
    double xi = 5e3;
    double xf = 5e4;
    double yi = 90.;
    double yf = 120.;

    gr2->Draw("apz");
    gr3->Draw("pz");
    gr2->SetTitle("");
    gr2->GetYaxis()->SetRangeUser(yi,yf);
    gr2->GetYaxis()->SetLabelSize(0.06);
    gr2->GetXaxis()->SetLabelSize(0.06);
    gr2->GetYaxis()->SetTitleFont(42);
    gr2->GetXaxis()->SetRangeUser(xi,xf); 
    grSigTotPBARP->Draw("csame");
    grSigTotPP->Draw("csame");
      
      canv1->SaveAs("XSECtot_BSG20orig_w002b_MMHT_member_39.eps");    
}
//----------------------------------------------------------------------
// Plotando a seção de choque elástica pp e pbp
//----------------------------------------------------------------------
void PlotElXSec(double *Wcm,double *XSecElpp,double *XSecElpbp,int npfit)
{
    FILE *aq1,*aq2,*aq3,*aq4;
// Plotting SigEl data     
    aq1 = fopen("SigElPBP.dat","r");
    aq2 = fopen("SigElPP.dat","r");
    aq3 = fopen("SigElPP_TOTEM.dat","r");
    aq4 = fopen("SigElPP_ATLAS.dat","r");

    // pbp
    const int npSigElpa = npSigElPBP;
    double W22[npSigElpa],Sig22Exp[npSigElpa],uW22[npSigElpa],uSig22Exp[npSigElpa];
   
    for(int i=0;i<npSigElpa;i++)     
        fscanf(aq1,"%lg %lg %lg %lg",&W22[i],&Sig22Exp[i],&uSig22Exp[i],&uW22[i]);
        
     // pp     
     const int npSigElPP = npSigElPPplot;
     double W12[npSigElPP],Sig12Exp[npSigElPP],uW12[npSigElPP],uSig12Exp[npSigElPP];
     
     for(int i=0;i<npSigElPP;i++)   
        fscanf(aq2,"%lg %lg %lg %lg",&W12[i],&Sig12Exp[i],&uSig12Exp[i],&uW12[i]);
       
       //---> TOTEM data <---
       
     const int npSigElTOTEM = npSigElPPTOTEM;
     double W12totem[npSigElTOTEM],Sig12ExpTOTEM[npSigElTOTEM],uW12totem[npSigElTOTEM],uSig12ExpTOTEM[npSigElTOTEM];  
       
     for(int i=0;i<npSigElTOTEM;i++)   
        fscanf(aq3,"%lg %lg %lg %lg",&W12totem[i],&Sig12ExpTOTEM[i],&uSig12ExpTOTEM[i],&uW12totem[i]);
        
       //---> ATLAS data <---
       
     const int npSigElATLAS = npSigElPPATLAS;
     double W12atlas[npSigElATLAS],Sig12ExpATLAS[npSigElATLAS],uW12atlas[npSigElATLAS],uSig12ExpATLAS[npSigElATLAS];  
       
     for(int i=0;i<npSigElATLAS;i++)   
        fscanf(aq4,"%lg %lg %lg %lg",&W12atlas[i],&Sig12ExpATLAS[i],&uSig12ExpATLAS[i],&uW12atlas[i]);
      
     fclose(aq1);
     fclose(aq2);
     fclose(aq3);
     fclose(aq4);
     
      TCanvas *canv1 = new TCanvas("c1","Elastic Cross Section",200,10,900,600);
      canv1->SetFillColor(0);
      canv1->SetLogx();

       double x_min(5e00),x_max(3.e04);
       double y_min(0.0),y_max(45.);
              
      TGraphErrors *gr0 = new TGraphErrors(npSigElpa,W22,Sig22Exp,uW22,uSig22Exp);
      gr0->SetMarkerStyle(24);
      gr0->SetMarkerColor(1);
      gr0->SetLineColor(1);
      gr0->SetMarkerSize(1.0);
      gr0->SetTitle();
      gr0->Draw("apz");
          
      TGraphErrors *gr1 = new TGraphErrors(npSigElPP,W12,Sig12Exp,uW12,uSig12Exp);
      gr1->SetMarkerStyle(20);
      gr1->SetMarkerColor(1);
      gr1->SetLineColor(1);
      gr1->SetMarkerSize(1.0);
      gr1->Draw("pz");         
      
      TGraphErrors *gr2 = new TGraphErrors(npSigElTOTEM,W12totem,Sig12ExpTOTEM,uW12totem,uSig12ExpTOTEM);
      gr2->SetMarkerStyle(33);
      gr2->SetMarkerColor(2);
      gr2->SetLineColor(2);
      gr2->SetMarkerSize(1.3);
      gr2->Draw("pz");  
      
      TGraphErrors *gr3 = new TGraphErrors(npSigElATLAS,W12atlas,Sig12ExpATLAS,uW12atlas,uSig12ExpATLAS);
      gr3->SetMarkerStyle(21);
      gr3->SetMarkerColor(4);
      gr3->SetLineColor(4);
      gr3->SetMarkerSize(0.8);
      gr3->Draw("pz");
      
      TGraph *grSigElPBARP = new TGraph(npfit,Wcm,XSecElpbp);
      grSigElPBARP->SetLineColor(4);
      grSigElPBARP->SetLineWidth(1);
      grSigElPBARP->SetLineStyle(1);
      grSigElPBARP->Draw("c");
      
      TGraph *grSigElPP = new TGraph(npfit,Wcm,XSecElpp);
      grSigElPP->SetLineColor(2);
      grSigElPP->SetLineWidth(1);
      grSigElPP->SetLineStyle(2);
      grSigElPP->Draw("c"); 

      gr0->GetYaxis()->SetTitle("#bf{#sigma_{el} [mb]}");
      gr0->GetYaxis()->SetRangeUser(y_min,y_max);
      gr0->GetYaxis()->SetTitleOffset(1.2);
      gr0->GetYaxis()->SetTitleFont(42);
      gr0->GetXaxis()->SetTitle("#bf{#sqrt{s} [GeV]}");
      gr0->GetXaxis()->SetTitleOffset(1.3);
      gr0->GetXaxis()->SetLimits(x_min,x_max);   
      gr0->GetXaxis()->SetTitleFont(42);
           
      TLegend *leg1 = new TLegend(0.65,0.15,0.89,0.42);
      leg1->AddEntry(gr0,"#bar{p}p data ","p");
      leg1->AddEntry(gr1,"pp data below LHC","p");
      leg1->AddEntry(gr2,"TOTEM","p");
      leg1->AddEntry(gr3,"ATLAS","p");
      leg1->AddEntry(grSigElPBARP,"#bar{p}p","l");
      leg1->AddEntry(grSigElPP,"pp","l");
      leg1->SetFillColor(0);
      leg1->SetShadowColor(0);
      leg1->SetBorderSize(0);
      leg1->SetTextSize(0.035);
      leg1->Draw();    

/* ---> desenhando um miniframe em um pad dentro do canvas canv1 <--- */

    TPad *pad1 = new TPad("pad1","This is pad1",0.15,0.51,0.55,0.89);
    pad1->SetFillColor(0);
    pad1->SetLogx();
    pad1->Draw();     
    pad1->cd();
    
    double xi = 5e3;
    double xf = 5e4;
    double yi = 20.0;
    double yf = 40.;

    gr2->Draw("apz");
    gr3->Draw("pz");
    gr2->SetTitle("");
    gr2->GetYaxis()->SetRangeUser(yi,yf);
    gr2->GetYaxis()->SetLabelSize(0.06);
    gr2->GetXaxis()->SetLabelSize(0.06);
    gr2->GetYaxis()->SetTitleFont(42);
    gr2->GetXaxis()->SetRangeUser(xi,xf); 
    grSigElPBARP->Draw("csame");
    grSigElPP->Draw("csame");
      
      canv1->SaveAs("XSECel_BSG20orig_w002b_MMHT_member_39.eps");    
}
//----------------------------------------------------------------------
// Plotando a seção de choque difrativa pp e pbp
//----------------------------------------------------------------------
void PlotDiffXSec(double *Wcm,double *XSecDiffpp,double *XSecDiffpbp,int npfit)
{
    FILE *aq1,*aq2,*aq3,*aq4,*aq5,*aq6,*aq7,*aq8;
// Plotting SigDiff data    
    //pp 
    aq1 = fopen("sigma_SD_pp_TOTEM_Cartiglia2013.dat","r");
    aq2 = fopen("sigma_SD_pp_CMS_Khachatryan2015.dat","r");
    aq3 = fopen("sigma_SD_pp_ALICE_Abelev2013.dat","r");
    //pbp
    aq4 = fopen("sigma_SD_pp_ISR_Armitage1982.dat","r");
    aq5 = fopen("sigma_SD_ppbar_UA5_Ansorge1986_Alner1987.dat","r");
    aq6 = fopen("sigma_SD_ppbar_CDF_Abe1994.dat","r");
    aq7 = fopen("sigma_SD_ppbar_E710_Amos1990_Amos1993.dat","r");
    aq8 = fopen("sigma_SD_ppbar_UA4_Bernard1987.dat","r");

    // pbp
    const int NPaq4 = npSigDiff4;
    double W23aq4[NPaq4],Sig23Exp4[NPaq4],uW23aq4[NPaq4],uSig23Exp4[NPaq4];
   
    for(int i=0;i<NPaq4;i++)     
        fscanf(aq4,"%lg %lg %lg %lg",&W23aq4[i],&Sig23Exp4[i],&uSig23Exp4[i],&uW23aq4[i]);
        
    const int NPaq5 = npSigDiff5;
    double W23aq5[NPaq5],Sig23Exp5[NPaq5],uW23aq5[NPaq5],uSig23Exp5[NPaq5];
   
    for(int i=0;i<NPaq5;i++)     
        fscanf(aq5,"%lg %lg %lg %lg",&W23aq5[i],&Sig23Exp5[i],&uSig23Exp5[i],&uW23aq5[i]);  
        
    const int NPaq6 = npSigDiff6;
    double W23aq6[NPaq6],Sig23Exp6[NPaq6],uW23aq6[NPaq6],uSig23Exp6[NPaq6];
   
    for(int i=0;i<NPaq6;i++)     
        fscanf(aq6,"%lg %lg %lg %lg",&W23aq6[i],&Sig23Exp6[i],&uSig23Exp6[i],&uW23aq6[i]); 
        
    const int NPaq7 = npSigDiff7;
    double W23aq7[NPaq7],Sig23Exp7[NPaq7],uW23aq7[NPaq7],uSig23Exp7[NPaq7];
   
    for(int i=0;i<NPaq7;i++)     
        fscanf(aq7,"%lg %lg %lg %lg",&W23aq7[i],&Sig23Exp7[i],&uSig23Exp7[i],&uW23aq7[i]);                  

    const int NPaq8 = npSigDiff8;
    double W23aq8[NPaq8],Sig23Exp8[NPaq8],uW23aq8[NPaq8],uSig23Exp8[NPaq8];
   
    for(int i=0;i<NPaq8;i++)     
        fscanf(aq8,"%lg %lg %lg %lg",&W23aq8[i],&Sig23Exp8[i],&uSig23Exp8[i],&uW23aq8[i]); 
        
    // pp     
    const int NPaq1 = npSigDiff1;
    double W13aq1[NPaq1],Sig13Exp1[NPaq1],uW13aq1[NPaq1],uSig13Exp1[NPaq1];
   
    for(int i=0;i<NPaq1;i++)     
        fscanf(aq1,"%lg %lg %lg %lg",&W13aq1[i],&Sig13Exp1[i],&uSig13Exp1[i],&uW13aq1[i]);
        
    const int NPaq2 = npSigDiff2;
    double W13aq2[NPaq2],Sig13Exp2[NPaq2],uW13aq2[NPaq2],uSig13Exp2[NPaq2];
   
    for(int i=0;i<NPaq2;i++)     
        fscanf(aq2,"%lg %lg %lg %lg",&W13aq2[i],&Sig13Exp2[i],&uSig13Exp2[i],&uW13aq2[i]);
        
    const int NPaq3 = npSigDiff3;
    double W13aq3[NPaq3],Sig13Exp3[NPaq3],uW13aq3[NPaq3],uSig13Exp3[NPaq3];
   
    for(int i=0;i<NPaq3;i++)     
        fscanf(aq3,"%lg %lg %lg %lg",&W13aq3[i],&Sig13Exp3[i],&uSig13Exp3[i],&uW13aq3[i]);              
       
     fclose(aq1);
     fclose(aq2);
     fclose(aq3);
     fclose(aq4);
     fclose(aq5);
     fclose(aq6);
     fclose(aq7);
     fclose(aq8);
     
      TCanvas *canv1 = new TCanvas("c1","Diffractive Cross Section",200,10,900,600);
      canv1->SetFillColor(0);
      canv1->SetLogx();
       
       double x_min(5e00),x_max(3.e04);
       double y_min(0.0),y_max(2.5e01);

      TGraphErrors *graq4 = new TGraphErrors(NPaq4,W23aq4,Sig23Exp4,uW23aq4,uSig23Exp4);
      graq4->SetMarkerStyle(20);
      graq4->SetMarkerColor(1);
      graq4->SetLineColor(1);
      graq4->SetMarkerSize(1.0);
      graq4->SetTitle();
      graq4->Draw("apz");
      
      TGraphErrors *graq5 = new TGraphErrors(NPaq5,W23aq5,Sig23Exp5,uW23aq5,uSig23Exp5);
      graq5->SetMarkerStyle(26);
      graq5->SetMarkerColor(4);
      graq5->SetLineColor(4);
      graq5->SetMarkerSize(1.0);
      graq5->SetTitle();
      graq5->Draw("pz");
      
      TGraphErrors *graq6 = new TGraphErrors(NPaq6,W23aq6,Sig23Exp6,uW23aq6,uSig23Exp6);
      graq6->SetMarkerStyle(28);
      graq6->SetMarkerColor(4);
      graq6->SetLineColor(4);
      graq6->SetMarkerSize(1.0);
      graq6->SetTitle();
      graq6->Draw("pz");
      
      TGraphErrors *graq7 = new TGraphErrors(NPaq7,W23aq7,Sig23Exp7,uW23aq7,uSig23Exp7);
      graq7->SetMarkerStyle(24);
      graq7->SetMarkerColor(4);
      graq7->SetLineColor(4);
      graq7->SetMarkerSize(1.0);
      graq7->SetTitle();
      graq7->Draw("pz");            
      
      TGraphErrors *graq8 = new TGraphErrors(NPaq8,W23aq8,Sig23Exp8,uW23aq8,uSig23Exp8);
      graq8->SetMarkerStyle(25);
      graq8->SetMarkerColor(4);
      graq8->SetLineColor(4);
      graq8->SetMarkerSize(1.0);
      graq8->SetTitle();
      graq8->Draw("pz");   
      
      TGraphErrors *graq1 = new TGraphErrors(NPaq1,W13aq1,Sig13Exp1,uW13aq1,uSig13Exp1);
      graq1->SetMarkerStyle(29);
      graq1->SetMarkerColor(2);
      graq1->SetLineColor(2);
      graq1->SetMarkerSize(1.0);
      graq1->SetTitle();
      graq1->Draw("pz");   
      
      TGraphErrors *graq2 = new TGraphErrors(NPaq2,W13aq2,Sig13Exp2,uW13aq2,uSig13Exp2);
      graq2->SetMarkerStyle(21);
      graq2->SetMarkerColor(2);
      graq2->SetLineColor(2);
      graq2->SetMarkerSize(1.0);
      graq2->SetTitle();
      graq2->Draw("pz");            
      
      TGraphErrors *graq3 = new TGraphErrors(NPaq3,W13aq3,Sig13Exp3,uW13aq3,uSig13Exp3);
      graq3->SetMarkerStyle(22);
      graq3->SetMarkerColor(2);
      graq3->SetLineColor(2);
      graq3->SetMarkerSize(1.0);
      graq3->SetTitle();
      graq3->Draw("pz");      
      
      TGraph *grSigDiffPBARP = new TGraph(npfit,Wcm,XSecDiffpbp);
      grSigDiffPBARP->SetLineColor(4);
      grSigDiffPBARP->SetLineWidth(1);
      grSigDiffPBARP->SetLineStyle(1);
      grSigDiffPBARP->Draw("c");
      
      TGraph *grSigDiffPP = new TGraph(npfit,Wcm,XSecDiffpp);
      grSigDiffPP->SetLineColor(2);
      grSigDiffPP->SetLineWidth(1);
      grSigDiffPP->SetLineStyle(2);
      grSigDiffPP->Draw("c");
      
      graq4->GetYaxis()->SetTitle("#bf{#sigma_{Diff} [mb]}");
      graq4->GetYaxis()->SetRangeUser(y_min,y_max);
      graq4->GetYaxis()->SetTitleOffset(1.2);
      graq4->GetYaxis()->SetTitleFont(42);
      graq4->GetXaxis()->SetTitle("#bf{#sqrt{s} [GeV]}");
      graq4->GetXaxis()->SetTitleOffset(1.3);
      graq4->GetXaxis()->SetLimits(x_min,x_max);   
      graq4->GetXaxis()->SetTitleFont(42);
      
      TLegend *leg1 = new TLegend(0.125,0.4,0.4,0.875);
      leg1->AddEntry(graq1,"TOTEM for 3.4 GeV < M_{SD} < 1100 GeV","p");
      leg1->AddEntry(graq2,"CMS for M_{X}^{2}/s < 0.05 or M_{Y}^{2}/s < 0.05","p");
      leg1->AddEntry(graq3,"ALICE for M_{X} < 200 GeV","p");      
      leg1->AddEntry(graq4,"ISR for M_{X}^{2}/s < 0.05","p");
      leg1->AddEntry(graq5,"UA5 for M_{X}^{2}/s < 0.05","p");
      leg1->AddEntry(graq6,"CDF for M_{X}^{2}/s < 0.2","p");
      leg1->AddEntry(graq7,"E710 for M_{X}^{2}/s < 0.05","p");
      leg1->AddEntry(graq8,"UA4 for M_{X}^{2}/s < 0.05","p");
      leg1->AddEntry(grSigDiffPBARP,"#bar{p}p","l");
      leg1->AddEntry(grSigDiffPP,"pp","l");
      leg1->SetFillColor(0);
      leg1->SetShadowColor(0);
      leg1->SetBorderSize(0);
      leg1->SetTextSize(0.035);
      leg1->Draw();    
    
      canv1->SaveAs("XSECdiff_BSG20orig_w002b_MMHT_member_39.eps");
}
//----------------------------------------------------------------------
// Plotando a seção de choque eik pp e pbp
//----------------------------------------------------------------------
void PlotEikXSec(double *Wcm,double *XSecEikpp,double *XSecEikpbp,int npfit)
{   
      TCanvas *canv1 = new TCanvas("c1","Eikonal Cross Section",200,10,900,600);
      canv1->SetFillColor(0);
      canv1->SetLogx();
       
       double x_min(5e00),x_max(3.e04);
       double y_min(5.0e01),y_max(8.0e02);
  
      TGraph *grSigEikPBARP = new TGraph(npfit,Wcm,XSecEikpbp);
      grSigEikPBARP->SetTitle();
      grSigEikPBARP->Draw("apz");
      grSigEikPBARP->SetLineColor(4);
      grSigEikPBARP->SetLineWidth(1);
      grSigEikPBARP->SetLineStyle(1);
      grSigEikPBARP->Draw("c");
      
      TGraph *grSigEikPP = new TGraph(npfit,Wcm,XSecEikpp);
      grSigEikPP->SetLineColor(2);
      grSigEikPP->SetLineWidth(1);
      grSigEikPP->SetLineStyle(2);
      grSigEikPP->Draw("c"); 

      grSigEikPBARP->GetYaxis()->SetTitle("#bf{#sigma_{eik} [mb]}");
      grSigEikPBARP->GetYaxis()->SetRangeUser(y_min,y_max);
      grSigEikPBARP->GetYaxis()->SetTitleOffset(1.2);
      grSigEikPBARP->GetYaxis()->SetTitleFont(42);
      grSigEikPBARP->GetXaxis()->SetTitle("#bf{#sqrt{s} [GeV]}");
      grSigEikPBARP->GetXaxis()->SetTitleOffset(1.3);
      grSigEikPBARP->GetXaxis()->SetLimits(x_min,x_max);   
      grSigEikPBARP->GetXaxis()->SetTitleFont(42);
           
      TLegend *leg1 = new TLegend(0.15,0.72,0.4,0.87);
      leg1->AddEntry(grSigEikPBARP,"#bar{p}p","l");
      leg1->AddEntry(grSigEikPP,"pp","l");
      leg1->SetFillColor(0);
      leg1->SetShadowColor(0);
      leg1->SetBorderSize(0);
      leg1->SetTextSize(0.035);
      leg1->Draw();    
     
      canv1->SaveAs("XSECeik_BSG20orig_w002b_MMHT_member_39.eps");    
}


//----------------------------------------------------------------------
// Plotando a seção de choque soft pp e pbp e hard
//----------------------------------------------------------------------
void PlotSoftHardXSec(double *Wcm,double *XSecSoftpp,double *XSecSoftpbp,double *XSecHard,int npfit)
{   
      TCanvas *canv1 = new TCanvas("c1","Soft and Hard Cross Section",200,10,900,600);
      canv1->SetFillColor(0);
      canv1->SetLogx();
       
       double x_min(5e00),x_max(3.e04);
       double y_min(0.0),y_max(1.5e02);
  
      TGraph *grSigSoftPBARP = new TGraph(npfit,Wcm,XSecSoftpbp);
      grSigSoftPBARP->SetTitle();
      grSigSoftPBARP->Draw("apz");
      grSigSoftPBARP->SetLineColor(4);
      grSigSoftPBARP->SetLineWidth(1);
      grSigSoftPBARP->SetLineStyle(1);
      grSigSoftPBARP->Draw("c");
      
      TGraph *grSigSoftPP = new TGraph(npfit,Wcm,XSecSoftpp);
      grSigSoftPP->SetLineColor(2);
      grSigSoftPP->SetLineWidth(1);
      grSigSoftPP->SetLineStyle(2);
      grSigSoftPP->Draw("c"); 

      TGraph *grSigHard = new TGraph(npfit,Wcm,XSecHard);
      grSigHard->SetLineColor(1);
      grSigHard->SetLineWidth(1);
      grSigHard->SetLineStyle(3);
      grSigHard->Draw("c"); 

      grSigSoftPBARP->GetYaxis()->SetTitle("#bf{#sigma [mb]}");
      grSigSoftPBARP->GetYaxis()->SetRangeUser(y_min,y_max);
      grSigSoftPBARP->GetYaxis()->SetTitleOffset(1.2);
      grSigSoftPBARP->GetYaxis()->SetTitleFont(42);
      grSigSoftPBARP->GetXaxis()->SetTitle("#bf{#sqrt{s} [GeV]}");
      grSigSoftPBARP->GetXaxis()->SetTitleOffset(1.3);
      grSigSoftPBARP->GetXaxis()->SetLimits(x_min,x_max);   
      grSigSoftPBARP->GetXaxis()->SetTitleFont(42);
           
      //TLegend *leg1 = new TLegend(0.15,0.72,0.4,0.87);
      TLegend *leg1 = new TLegend(0.72,0.15,0.88,0.30);
      leg1->AddEntry(grSigSoftPBARP,"#sigma^{#bar{p}p}_{soft}","l");
      leg1->AddEntry(grSigSoftPP,"#sigma^{pp}_{soft}","l");
      leg1->AddEntry(grSigHard,"#sigma_{pQCD}","l");
      leg1->SetFillColor(0);
      leg1->SetShadowColor(0);
      leg1->SetBorderSize(0);
      leg1->SetTextSize(0.035);
      leg1->Draw();    
     
      canv1->SaveAs("XSECsofthard_BSG20orig_w002b_MMHT_member_39.eps");    
}

//----------------------------------------------------------------------
void fcn(int &npar, double *gin, double &f, double *par, int iflag)
{    
    /*******************************************************/ 
    //fit parameters
    double ww=par[0]; 
    double Sig0=par[1];
    double A1=par[2];
    double del1=par[3];
    double A2=par[4];  
    double del2=par[5];
    double B=1.0;
    /*******************************************************/
    //reading data, then fitting

    //
    // St pp e pbp
    // 
    FILE *aq1,*aq2;
 
    // pbp   
    aq1 = fopen("paw_Stpa.dat","r");
    aq2 = fopen("paw_Stpp.dat","r");    

    const int npStpa = npSigTotPBP;
    double W21[npStpa],Sig21Exp[npStpa],uW21[npStpa],uSig21Exp[npStpa];
   
    for(int i=0;i<npStpa;i++)     
        fscanf(aq1,"%lg %lg %lg %lg",&W21[i],&Sig21Exp[i],&uSig21Exp[i],&uW21[i]);
        
     // pp   
     const int npStpp = npSigTotPP;
     double W11[npStpp],Sig11Exp[npStpp],uW11[npStpp],uSig11Exp[npStpp];
     
     for(int i=0;i<npStpp;i++)   
        fscanf(aq2,"%lg %lg %lg %lg",&W11[i],&Sig11Exp[i],&uSig11Exp[i],&uW11[i]);
      
     fclose(aq1);
     fclose(aq2);

    /*******************************************************/
    //Total XSection - pp
       
    TF1 SigTotPP("SigTotPP",KerForward,bmin,bmax,10);
    SigTotPP.SetParameter(0,ww);
    SigTotPP.SetParameter(1,Sig0);
    SigTotPP.SetParameter(2,A1);
    SigTotPP.SetParameter(3,del1);
    SigTotPP.SetParameter(4,A2);
    SigTotPP.SetParameter(5,del2);
    SigTotPP.SetParameter(6,B);
    SigTotPP.SetParameter(7,1);
    SigTotPP.SetParameter(8,1);    
    
    /*******************************************************/
    //Total XSection - pbp
       
    TF1 SigTotPBARP("SigTotPBARP",KerForward,bmin,bmax,10);
    SigTotPBARP.SetParameter(0,ww);
    SigTotPBARP.SetParameter(1,Sig0);
    SigTotPBARP.SetParameter(2,A1);
    SigTotPBARP.SetParameter(3,del1);
    SigTotPBARP.SetParameter(4,A2);
    SigTotPBARP.SetParameter(5,del2);
    SigTotPBARP.SetParameter(6,B);
    SigTotPBARP.SetParameter(7,2);
    SigTotPBARP.SetParameter(8,1);
        
    /*******************************************************/
    
    //Chi-squared calculation
    int np = 50;
    double *x=new double[np];
    double *w=new double[np];
    double Stpp,Stpbarp;    
    
    double delta;
    double chisq=0.;
    
    for(int i=0;i<npStpp;i++){
        SigTotPP.SetParameter(9,W11[i]);
        SigTotPP.CalcGaussLegendreSamplingPoints(np,x,w,1e-15);
        Stpp=SigTotPP.IntegralFast(np,x,w,bmin,bmax);
        delta=(Sig11Exp[i]-Stpp)/uSig11Exp[i];
        chisq+=Power(delta,2);          
   }
       
    delta=0.;

    for(int i=0;i<npStpa;i++){
        SigTotPBARP.SetParameter(9,W21[i]);
        SigTotPBARP.CalcGaussLegendreSamplingPoints(np,x,w,1e-15);
        Stpbarp=SigTotPBARP.IntegralFast(np,x,w,bmin,bmax);
        delta=(Sig21Exp[i]-Stpbarp)/uSig21Exp[i];
        chisq+=Power(delta,2);           
   }
   
    delta=0.; 
       
   f = chisq;        
}
//----------------------------------------------------------------------
void BSG20orig_w002b_MMHT_member_39()
{
    /******************************************/
    //time control - start
    auto start = high_resolution_clock::now();    
    /******************************************/ 
         
    TMinuit *gMinuit = new TMinuit(numpar);  //initialize TMinuit with a maximum of 4 params
    gMinuit->SetFCN(fcn);
    
    double arglist[numpar];
    int ierflg = 0;
  
    arglist[0] = 7.04; // --> número qualaquer //fitting with CL = 1\sigma
    gMinuit->mnexcm("SET ERR",arglist,1,ierflg); 
      
    double vstart[numpar] = 
{
    1.42511e+00,
    8.53035e+01,
    2.53558e+01,
    1.60817e+00,
    2.27573e+01,
    5.75454e-01}; //start values */
    double step[numpar] = {1.e-3,1.e-3,1.e-3,1.e-3,1.e-3,1.e-3}; //steps
    gMinuit->mnparm(0, "w",    vstart[0],   step[0], 0,0,ierflg);
    gMinuit->mnparm(1, "Sig0", vstart[1],   step[1], 0,0,ierflg);
    gMinuit->mnparm(2, "A1",   vstart[2],   step[2], 0,0,ierflg);
    gMinuit->mnparm(3, "del1", vstart[3],   step[3], 0,0,ierflg);
    gMinuit->mnparm(4, "A2",   vstart[4],   step[4], 0,0,ierflg);
    gMinuit->mnparm(5, "del2", vstart[5],   step[5], 0,0,ierflg);
//    gMinuit->mnparm(3, "B",    vstart[3],   step[3], 0,0,ierflg);

/*
 FCN=1029.54 FROM MIGRAD    STATUS=CONVERGED      25 CALLS         332 TOTAL
                     EDM=8.05383e-11    STRATEGY= 1  ERROR MATRIX UNCERTAINTY   1.8 per cent
  EXT PARAMETER                                   STEP         FIRST   
  NO.   NAME      VALUE            ERROR          SIZE      DERIVATIVE 
   1  w            1.42511e+00   9.90409e-02  -5.90193e-07   7.71800e-04
   2  Sig0         8.53035e+01   2.02440e+00   8.97189e-07  -5.26809e-05
   3  A1           2.53558e+01   2.42270e+00  -6.54695e-05  -1.51312e-05
   4  del1         1.60817e+00   1.25338e-01   2.65169e-06   4.55755e-04
   5  A2           2.27573e+01   1.50058e+00  -5.89936e-05   7.31768e-05
   6  del2         5.75454e-01   4.24671e-02  -1.85834e-06  -1.63778e-03
                               ERR DEF= 7.04
*/

        gMinuit->FixParameter(0);
        gMinuit->FixParameter(1);
        gMinuit->FixParameter(2);
        gMinuit->FixParameter(3);
        gMinuit->FixParameter(4);
        gMinuit->FixParameter(5);

    //start minimizing data
    arglist[0] = 1;
    gMinuit->mnexcm("SET STR", arglist ,1,ierflg);
    arglist[0] = 500;    
    gMinuit->mnexcm("MIGRAD", arglist ,1,ierflg);    
     gMinuit->mnexcm("MIGRAD", arglist ,1,ierflg);
     gMinuit->mnexcm("MIGRAD", arglist ,1,ierflg);
//    gMinuit->mnexcm("MIGRAD", arglist ,1,ierflg);   
   
      //status of minimization 
   double amin,edm,errdef;
   int nvpar,nparx,icstat;  
   
   gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
  //void mnstat(double &fmin, double &fedm, double &errdef, int &npari, int &nparx, int &istat) 
  //*-*-*-*-*Returns concerning the current status of the minimization*-*-*-*-*
  //*-*      =========================================================
  //*-*       User-called
  //*-*          Namely, it returns:
  //*-*        FMIN: the best function value found so far
  //*-*        FEDM: the estimated vertical distance remaining to minimum
  //*-*        ERRDEF: the value of UP defining parameter uncertainties
  //*-*        NPARI: the number of currently variable parameters
  //*-*        NPARX: the highest (external) parameter number defined by user
  //*-*        ISTAT: a status integer indicating how good is the covariance
  //*-*           matrix:  0= not calculated at all
  //*-*                    1= approximation only, not accurate
  //*-*                    2= full matrix, but forced positive-definite
  //*-*                    3= full accurate covariance matrix
  //*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
   
  cout << "\n";
  
  cout << " **Statistical information**\n\n";
  cout << " Minimum chi square = " << amin << "\n";
  cout << " Degrees of freedom = " << (npMin-numpar) << "\n";
  cout << " Chi²/DOF = " << amin/(npMin-numpar) << "\n";
  cout << " Integrated probability = "<< Prob(amin,npMin-numpar) << "\n";
  cout << " Estimated vert. distance to min. = " << edm << "\n";
  cout << " Number of variable parameters = " << nvpar << "\n";
  cout << " Highest number of parameters defined by user = " << nparx << "\n";
  cout << " Status of covariance matrix = " << icstat << "\n\n";  
     
   double outpar[numpar], err[numpar];
  
   for (int i=0; i<numpar; i++)
    gMinuit->GetParameter(i,outpar[i],err[i]);     
    
    /*******************************************************/ 
    //fit parameters
    double ww=outpar[0]; 
    double Sig0=outpar[1];
    double A1=outpar[2];
    double del1=outpar[3];
    double A2=outpar[4];  
    double del2=outpar[5];
    double B=1.0;
    /*******************************************************/
    //Total XSection - pp
       
    TF1 SigTotPP("SigTotPP",KerForward,bmin,bmax,10);
    SigTotPP.SetParameter(0,ww);
    SigTotPP.SetParameter(1,Sig0);
    SigTotPP.SetParameter(2,A1);
    SigTotPP.SetParameter(3,del1);
    SigTotPP.SetParameter(4,A2);
    SigTotPP.SetParameter(5,del2);
    SigTotPP.SetParameter(6,B);
    SigTotPP.SetParameter(7,1);  
    SigTotPP.SetParameter(8,1);  
    
    /*******************************************************/
    //Total XSection - pbp
       
    TF1 SigTotPBARP("SigTotPBARP",KerForward,bmin,bmax,10);
    SigTotPBARP.SetParameter(0,ww);
    SigTotPBARP.SetParameter(1,Sig0);
    SigTotPBARP.SetParameter(2,A1);
    SigTotPBARP.SetParameter(3,del1);
    SigTotPBARP.SetParameter(4,A2);
    SigTotPBARP.SetParameter(5,del2);
    SigTotPBARP.SetParameter(6,B);
    SigTotPBARP.SetParameter(7,2);
    SigTotPBARP.SetParameter(8,1);
    
    /*******************************************************/    
    //Elastic XSection - pp
       
    TF1 SigElPP("SigElPP",KerForward,bmin,bmax,10);
    SigElPP.SetParameter(0,ww);
    SigElPP.SetParameter(1,Sig0);
    SigElPP.SetParameter(2,A1);
    SigElPP.SetParameter(3,del1);
    SigElPP.SetParameter(4,A2);
    SigElPP.SetParameter(5,del2);
    SigElPP.SetParameter(6,B);
    SigElPP.SetParameter(7,1);  
    SigElPP.SetParameter(8,2);  
    
    /*******************************************************/
    //Elastic XSection - pbp
       
    TF1 SigElPBARP("SigElPBARP",KerForward,bmin,bmax,10);
    SigElPBARP.SetParameter(0,ww);
    SigElPBARP.SetParameter(1,Sig0);
    SigElPBARP.SetParameter(2,A1);
    SigElPBARP.SetParameter(3,del1);
    SigElPBARP.SetParameter(4,A2);
    SigElPBARP.SetParameter(5,del2);
    SigElPBARP.SetParameter(6,B);
    SigElPBARP.SetParameter(7,2);
    SigElPBARP.SetParameter(8,2);
    
    /*******************************************************/    
    //Diff XSection - pp
       
    TF1 SigDiffPP("SigDiffPP",KerForward,bmin,bmax,10);
    SigDiffPP.SetParameter(0,ww);
    SigDiffPP.SetParameter(1,Sig0);
    SigDiffPP.SetParameter(2,A1);
    SigDiffPP.SetParameter(3,del1);
    SigDiffPP.SetParameter(4,A2);
    SigDiffPP.SetParameter(5,del2);
    SigDiffPP.SetParameter(6,B);
    SigDiffPP.SetParameter(7,1);  
    SigDiffPP.SetParameter(8,3);  
    
    /*******************************************************/
    //Diff XSection - pbp
       
    TF1 SigDiffPBARP("SigDiffPBARP",KerForward,bmin,bmax,10);
    SigDiffPBARP.SetParameter(0,ww);
    SigDiffPBARP.SetParameter(1,Sig0);
    SigDiffPBARP.SetParameter(2,A1);
    SigDiffPBARP.SetParameter(3,del1);
    SigDiffPBARP.SetParameter(4,A2);
    SigDiffPBARP.SetParameter(5,del2);
    SigDiffPBARP.SetParameter(6,B);
    SigDiffPBARP.SetParameter(7,2);
    SigDiffPBARP.SetParameter(8,3); 
    
    /*******************************************************/    
    //Eik XSection - pp
       
    TF1 SigEikPP("SigEikPP",KerThings,Wmin,Wmax,9);
    SigEikPP.SetParameter(0,ww);
    SigEikPP.SetParameter(1,Sig0);
    SigEikPP.SetParameter(2,A1);
    SigEikPP.SetParameter(3,del1);
    SigEikPP.SetParameter(4,A2);
    SigEikPP.SetParameter(5,del2);
    SigEikPP.SetParameter(6,B);
    SigEikPP.SetParameter(7,1);
    SigEikPP.SetParameter(8,1);
    
    /*******************************************************/
    //Eik XSection - pbp
       
    TF1 SigEikPBARP("SigEikPBARP",KerThings,Wmin,Wmax,9);
    SigEikPBARP.SetParameter(0,ww);
    SigEikPBARP.SetParameter(1,Sig0);
    SigEikPBARP.SetParameter(2,A1);
    SigEikPBARP.SetParameter(3,del1);
    SigEikPBARP.SetParameter(4,A2);
    SigEikPBARP.SetParameter(5,del2);
    SigEikPBARP.SetParameter(6,B);
    SigEikPBARP.SetParameter(7,2);     
    SigEikPBARP.SetParameter(8,1);    
        
    /*******************************************************/  
    
    //Soft XSection - pp
       
    TF1 SigSoftPP("SigSoftPP",KerThings,Wmin,Wmax,9);
    SigSoftPP.SetParameter(0,ww);
    SigSoftPP.SetParameter(1,Sig0);
    SigSoftPP.SetParameter(2,A1);
    SigSoftPP.SetParameter(3,del1);
    SigSoftPP.SetParameter(4,A2);
    SigSoftPP.SetParameter(5,del2);
    SigSoftPP.SetParameter(6,B);
    SigSoftPP.SetParameter(7,1);
    SigSoftPP.SetParameter(8,2);
    
    /*******************************************************/
    //Soft XSection - pbp
       
    TF1 SigSoftPBARP("SigSoftPBARP",KerThings,Wmin,Wmax,9);
    SigSoftPBARP.SetParameter(0,ww);
    SigSoftPBARP.SetParameter(1,Sig0);
    SigSoftPBARP.SetParameter(2,A1);
    SigSoftPBARP.SetParameter(3,del1);
    SigSoftPBARP.SetParameter(4,A2);
    SigSoftPBARP.SetParameter(5,del2);
    SigSoftPBARP.SetParameter(6,B);
    SigSoftPBARP.SetParameter(7,2);     
    SigSoftPBARP.SetParameter(8,2);    
    
    /*******************************************************/  
    
    //Hard XSection
       
    TF1 SigHard("SigHard",KerThings,Wmin,Wmax,9);
    SigHard.SetParameter(0,ww);
    SigHard.SetParameter(1,Sig0);
    SigHard.SetParameter(2,A1);
    SigHard.SetParameter(3,del1);
    SigHard.SetParameter(4,A2);
    SigHard.SetParameter(5,del2);
    SigHard.SetParameter(6,B);
    SigHard.SetParameter(7,2);
    SigHard.SetParameter(8,3);    
        
    /*******************************************************/     
//----------------------------------------------------------------------
   //calculates total, elastic, diff and eik xsection, creating files to make the plots
    FILE *aq1,*aq2,*aq3,*aq4,*aq5,*aq6,*aq7,*aq8,*aq9,*aq10,*aq11;
 
    aq1=fopen("XSecTotPP_BSG20orig_w002b_MMHT_member_39.dat","w"); 
    aq2=fopen("XSecTotPBP_BSG20orig_w002b_MMHT_member_39.dat","w");
    aq3=fopen("XSecElPP_BSG20orig_w002b_MMHT_member_39.dat","w"); 
    aq4=fopen("XSecElPBP_BSG20orig_w002b_MMHT_member_39.dat","w");
    aq5=fopen("XSecDiffPP_BSG20orig_w002b_MMHT_member_39.dat","w"); 
    aq6=fopen("XSecDiffPBP_BSG20orig_w002b_MMHT_member_39.dat","w");  
    aq7=fopen("XSecEikPP_BSG20orig_w002b_MMHT_member_39.dat","w"); 
    aq8=fopen("XSecEikPBP_BSG20orig_w002b_MMHT_member_39.dat","w"); 
    aq9=fopen("XSecSoftPP_BSG20orig_w002b_MMHT_member_39.dat","w"); 
    aq10=fopen("XSecSoftPBP_BSG20orig_w002b_MMHT_member_39.dat","w");             
    aq11=fopen("XSecHard_BSG20orig_w002b_MMHT_member_39.dat","w"); 
           
     //number of directive computing for GaussLegendreIntegration
    int np = 50;
    double *x=new double[np];
    double *w=new double[np];
    
    double i=0.;
    const double di=0.05;
    double Ecm;
    double Stpp,Stpbarp,Selpp,Selpbarp,Sdiffpp,Sdiffpbarp,Seikpp,Seikpbarp,Ssoftpp,Ssoftpbarp,Shard;
    
     do{
            Ecm=5.*pow(10.,i);
            /*************************************************/
            //pp
            SigTotPP.SetParameter(9,Ecm);
            SigTotPP.CalcGaussLegendreSamplingPoints(np,x,w,1e-15);
            Stpp = SigTotPP.IntegralFast(np,x,w,bmin,bmax);   
            
            SigElPP.SetParameter(9,Ecm);
            SigElPP.CalcGaussLegendreSamplingPoints(np,x,w,1e-15);
            Selpp = SigElPP.IntegralFast(np,x,w,bmin,bmax);
            
            SigDiffPP.SetParameter(9,Ecm);
            SigDiffPP.CalcGaussLegendreSamplingPoints(np,x,w,1e-15);
            Sdiffpp = SigDiffPP.IntegralFast(np,x,w,bmin,bmax);
            
            Seikpp = SigEikPP.Eval(Ecm);
            
            Ssoftpp = SigSoftPP.Eval(Ecm);
                 
            fprintf(aq1,"%.2e   %.4lf\n",Ecm,Stpp); 
            fprintf(aq3,"%.2e   %.4lf\n",Ecm,Selpp);
            fprintf(aq5,"%.2e   %.4lf\n",Ecm,Sdiffpp);
            fprintf(aq7,"%.2e   %.4lf\n",Ecm,Seikpp);
            fprintf(aq9,"%.2e   %.4lf\n",Ecm,Ssoftpp);
            /*************************************************/
            //pbarp
            SigTotPBARP.SetParameter(9,Ecm);
            SigTotPBARP.CalcGaussLegendreSamplingPoints(np,x,w,1e-15);
            Stpbarp = SigTotPBARP.IntegralFast(np,x,w,bmin,bmax);
            
            SigElPBARP.SetParameter(9,Ecm);
            SigElPBARP.CalcGaussLegendreSamplingPoints(np,x,w,1e-15);
            Selpbarp = SigElPBARP.IntegralFast(np,x,w,bmin,bmax);
            
            SigDiffPBARP.SetParameter(9,Ecm);
            SigDiffPBARP.CalcGaussLegendreSamplingPoints(np,x,w,1e-15);
            Sdiffpbarp = SigDiffPBARP.IntegralFast(np,x,w,bmin,bmax);
            
            Seikpbarp = SigEikPBARP.Eval(Ecm);
            
            Ssoftpbarp = SigSoftPBARP.Eval(Ecm);
            
            fprintf(aq2,"%.2e   %.4lf\n",Ecm,Stpbarp);  
            fprintf(aq4,"%.2e   %.4lf\n",Ecm,Selpbarp);
            fprintf(aq6,"%.2e   %.4lf\n",Ecm,Sdiffpbarp); 
            fprintf(aq8,"%.2e   %.4lf\n",Ecm,Seikpbarp);
            fprintf(aq10,"%.2e   %.4lf\n",Ecm,Ssoftpbarp);    
            /*************************************************/      
            
            Shard = SigHard.Eval(Ecm);
            
            fprintf(aq11,"%.2e   %.4lf\n",Ecm,Shard);
            /*************************************************/
                                   
            i+=di;
                       
       //cout << Ecm << "\t" << Seikpp << "\t" << Seikpbarp << endl;
            
       }while(Ecm<Wmax); 
       
      fclose(aq1);
      fclose(aq2); 
      fclose(aq3);
      fclose(aq4);
      fclose(aq5);
      fclose(aq6);   
      fclose(aq7);
      fclose(aq8); 
      fclose(aq9);   
      fclose(aq10);
      fclose(aq11);           
       
      aq1 = fopen("XSecTotPP_BSG20orig_w002b_MMHT_member_39.dat","r");
      aq2 = fopen("XSecTotPBP_BSG20orig_w002b_MMHT_member_39.dat","r");
      aq3 = fopen("XSecElPP_BSG20orig_w002b_MMHT_member_39.dat","r");
      aq4 = fopen("XSecElPBP_BSG20orig_w002b_MMHT_member_39.dat","r");
      aq5 = fopen("XSecDiffPP_BSG20orig_w002b_MMHT_member_39.dat","r"); 
      aq6 = fopen("XSecDiffPBP_BSG20orig_w002b_MMHT_member_39.dat","r");  
      aq7 = fopen("XSecEikPP_BSG20orig_w002b_MMHT_member_39.dat","r"); 
      aq8 = fopen("XSecEikPBP_BSG20orig_w002b_MMHT_member_39.dat","r"); 
      aq9 = fopen("XSecSoftPP_BSG20orig_w002b_MMHT_member_39.dat","r"); 
      aq10 = fopen("XSecSoftPBP_BSG20orig_w002b_MMHT_member_39.dat","r");             
      aq11 = fopen("XSecHard_BSG20orig_w002b_MMHT_member_39.dat","r");                   
      
      const int npfit=108;
      double Wcm[npfit],XSecTOTpp[npfit],XSecTOTpbp[npfit],XSecElpp[npfit],XSecElpbp[npfit],XSecDiffpp[npfit],XSecDiffpbp[npfit],XSecEikpp[npfit],XSecEikpbp[npfit],XSecSoftpp[npfit],XSecSoftpbp[npfit],XSecHard[npfit];       
            
      for(int j=0;j<npfit;j++){
          fscanf(aq1,"%lg %lg",&Wcm[j],&XSecTOTpp[j]);
          fscanf(aq2,"%lg %lg",&Wcm[j],&XSecTOTpbp[j]);
          fscanf(aq3,"%lg %lg",&Wcm[j],&XSecElpp[j]);
          fscanf(aq4,"%lg %lg",&Wcm[j],&XSecElpbp[j]);
          fscanf(aq5,"%lg %lg",&Wcm[j],&XSecDiffpp[j]);
          fscanf(aq6,"%lg %lg",&Wcm[j],&XSecDiffpbp[j]);
          fscanf(aq7,"%lg %lg",&Wcm[j],&XSecEikpp[j]);
          fscanf(aq8,"%lg %lg",&Wcm[j],&XSecEikpbp[j]);          
          fscanf(aq9,"%lg %lg",&Wcm[j],&XSecSoftpp[j]);
          fscanf(aq10,"%lg %lg",&Wcm[j],&XSecSoftpbp[j]);                              
          fscanf(aq11,"%lg %lg",&Wcm[j],&XSecHard[j]);          
      }
      
      fclose(aq1);
      fclose(aq2);
      fclose(aq3);
      fclose(aq4);    
      fclose(aq5);
      fclose(aq6);
      fclose(aq7);
      fclose(aq8);
      fclose(aq9);   
      fclose(aq10);
      fclose(aq11);
             
//----------------------------------------------------------------------
    //Plotting xsection data         
    PlotTotXSec(Wcm,XSecTOTpp,XSecTOTpbp,npfit);      
    PlotElXSec(Wcm,XSecElpp,XSecElpbp,npfit);
    PlotDiffXSec(Wcm,XSecDiffpp,XSecDiffpbp,npfit);
    PlotEikXSec(Wcm,XSecEikpp,XSecEikpbp,npfit);
    PlotSoftHardXSec(Wcm,XSecSoftpp,XSecSoftpbp,XSecHard,npfit);
    /******************************************/
    auto stop = high_resolution_clock::now(); 
    auto duration = duration_cast<seconds>(stop - start); 
    cout << "Time(in min): "<<duration.count()/6.e1 << endl;
    /******************************************/         
}
