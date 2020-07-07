/* PROGRAM: BSG20 based on Lipari-Lusignoli original model

   ---> Ref. Original model: PHYSICAL REVIEW D 80, 074014 (2009)
        SciHub link: https://sci-hub.tw/10.1103/PhysRevD.80.074014
        
   Modified version: SigEik=SigSoft+SigPQCD

------------------------------------------------------------------------
   *** Main characteristics ***

	> SigEik = SigSoft + SigPQCD
	> SigSoft = A1*s^{-del1} +/- A2*s^{-del2} + Sig0
	> SigPQCD = B*SigPQCD
	
	> Calculates Total and Elastic xsec pp and pbp

   *** Observations ***

 Obs1: TOTEM + ATLAS dataset
 Obs2: Real parametrization for SigQCD with 10 free-parameters
 Obs3: PDF: CT14
 Obs4: Free: w, Sig0, A1, A2, del1, del2
 Obs5: B is fixed in 2.0

   *** Important *** 
   (See: https://arxiv.org/abs/1908.01040 Chapter 8, sections 8.4 and 8.8)
   
   mg=0.4d0 !GeV -> mg=500+/-200 MeV for lmb=300 MeV
 
 --> (MSTW) lbd=0.318920 GeV lambda equivalent to alpha_LO(Mz^2)=0.13939~0.139 macthing-prescrition squeme; m_b=4.18 GeV
 
 --> (MMHT) lbd=0.267471 GeV lambda equivalent to alpha_LO(Mz^2)=0.13499~0.135 macthing-prescrition squeme; m_b=4.18 GeV
 
 --> (CTEQ6L) lbd=0.3260 GeV lambda equivalent to alpha_NLO(Mz^2)=0.118 for 4 flavors
 
 --> (CT14)   lbd=0.3260 GeV lambda equivalent to alpha_NLO(Mz^2)=0.118 for 4 flavors

   *** How to run the code ***
 
 root -l -q /usr/lib/x86_64-linux-gnu/libgsl.so /usr/lib/x86_64-linux-gnu/libgslcblas.so BSG20orig_w002teste.C | cat > BSG20orig_w002teste.min
 
   ---------------------------------------------------------------------
   Codes by: M.Broilo*, P.V.R.G.Silva, V.P.B.Gonçalves
   mateus.broilo@ufrgs.br | mateus.broilo90@gmail.com
   High and Medium Energy Group
   Grupo de Altas e Medias Energias (GAME)
   Universidade Federal de Pelotas, Pelotas - RS (Brasil)
   ---------------------------------------------------------------------
   
   Creation: 18/dez/2019 (Pelotas, RS, Brazil)
   Last update: 11/mar/2020 (Pelotas, RS, Brazil)
   
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

// Plots
const int npSigTotPPplot=104;
const int npSigTotPPTOTEM=12;
const int npSigTotPPATLAS=2;
const int npSigElPPplot=147;
const int npSigElPPTOTEM=6;
const int npSigElPPATLAS=2;
//----------------------------------------------------------------------
//parameters to control energy and b range (for plots and integrations)
const double Wmin=5.;
const double Wmax=1.e6;
const double bmin=0.;
const double bmax=30.;
//----------------------------------------------------------------------
//QCD cross section - analytical parametrization for CT14
//----------------------------------------------------------------------
const double c1=1.01;
const double c2=1.05;
const double c3=1.09;
const double bkg=100.;
// SigQCD(GG e QG) parameters
const double b1=100.22;
const double b2=0.43425e-1;
const double b3=1.2736;
const double b4=1.9189;
const double b5=0.12198e-7;
const double b6=14.050;
const double b7=0.50348;
const double b8=3699.4;
const double b9=-80.280;
const double b10=-2.6318;  
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
//Forward observables - total xsection
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
    double nn,xx;
    
    	if(reac==1) //pp
    	{
    	/*******************************************************/
    	//Average number of interactions
    	/*******************************************************/
	nn=(SigEikPP*mb_to_fm2)*Wdip;
	xx=(nn*ww*ww)/2.;
    	}
    	else if(reac==2) //pbp
    	{
	nn=(SigEikPBP*mb_to_fm2)*Wdip;
	xx=(nn*ww*ww)/2.;
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
        else
        {
            return 0.;
        }
}
//----------------------------------------------------------------------
// Plotando a seção de choque total pp e pbp
//----------------------------------------------------------------------
void PlotTotXSec(double *Wcm,double *XSecTOTpp,double *XSecTOTpbp,int npfit)
{
    FILE *aq1,*aq2,*aq3,*aq4;
// Plotting St data     
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
       
       double x_min(5e00),x_max(2.e05);
       double y_min(3.0e01),y_max(1.8e02);

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
      
      TGraph *grSigPBARP_Eik = new TGraph(npfit,Wcm,XSecTOTpbp);
      grSigPBARP_Eik->SetLineColor(4);
      grSigPBARP_Eik->SetLineWidth(1);
      grSigPBARP_Eik->SetLineStyle(1);
      grSigPBARP_Eik->Draw("c");
      
      TGraph *grSigPP_Eik = new TGraph(npfit,Wcm,XSecTOTpp);
      grSigPP_Eik->SetLineColor(2);
      grSigPP_Eik->SetLineWidth(1);
      grSigPP_Eik->SetLineStyle(2);
      grSigPP_Eik->Draw("c"); 

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
      leg1->AddEntry(grSigPBARP_Eik,"#bar{p}p","l");
      leg1->AddEntry(grSigPP_Eik,"pp","l");
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
    grSigPBARP_Eik->Draw("csame");
    grSigPP_Eik->Draw("csame");
      
      canv1->SaveAs("XSECtot_BSG20orig_w002teste.eps");    
}
//----------------------------------------------------------------------
// Plotando a seção de choque elástica pp e pbp
//----------------------------------------------------------------------
void PlotElXSec(double *Wcm,double *XSecElpp,double *XSecElpbp,int npfit)
{
    FILE *aq1,*aq2,*aq3,*aq4;
// Plotting St data     
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
       
     const int npStTOTEM = npSigElPPTOTEM;
     double W12totem[npStTOTEM],Sig12ExpTOTEM[npStTOTEM],uW12totem[npStTOTEM],uSig12ExpTOTEM[npStTOTEM];  
       
     for(int i=0;i<npStTOTEM;i++)   
        fscanf(aq3,"%lg %lg %lg %lg",&W12totem[i],&Sig12ExpTOTEM[i],&uSig12ExpTOTEM[i],&uW12totem[i]);
        
       //---> ATLAS data <---
       
     const int npStATLAS = npSigElPPATLAS;
     double W12atlas[npStATLAS],Sig12ExpATLAS[npStATLAS],uW12atlas[npStATLAS],uSig12ExpATLAS[npStATLAS];  
       
     for(int i=0;i<npStATLAS;i++)   
        fscanf(aq4,"%lg %lg %lg %lg",&W12atlas[i],&Sig12ExpATLAS[i],&uSig12ExpATLAS[i],&uW12atlas[i]);
      
     fclose(aq1);
     fclose(aq2);
     fclose(aq3);
     fclose(aq4);
     
      TCanvas *canv1 = new TCanvas("c1","Total Cross Section",200,10,900,600);
      canv1->SetFillColor(0);
      canv1->SetLogx();

       double x_min(5e00),x_max(2.e05);
       double y_min(0.0),y_max(60.);
              
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
      
      TGraphErrors *gr2 = new TGraphErrors(npStTOTEM,W12totem,Sig12ExpTOTEM,uW12totem,uSig12ExpTOTEM);
      gr2->SetMarkerStyle(33);
      gr2->SetMarkerColor(2);
      gr2->SetLineColor(2);
      gr2->SetMarkerSize(1.3);
      gr2->Draw("pz");  
      
      TGraphErrors *gr3 = new TGraphErrors(npStATLAS,W12atlas,Sig12ExpATLAS,uW12atlas,uSig12ExpATLAS);
      gr3->SetMarkerStyle(21);
      gr3->SetMarkerColor(4);
      gr3->SetLineColor(4);
      gr3->SetMarkerSize(0.8);
      gr3->Draw("pz");
      
      TGraph *grSigPBARP_Eik = new TGraph(npfit,Wcm,XSecElpbp);
      grSigPBARP_Eik->SetLineColor(4);
      grSigPBARP_Eik->SetLineWidth(1);
      grSigPBARP_Eik->SetLineStyle(1);
      grSigPBARP_Eik->Draw("c");
      
      TGraph *grSigPP_Eik = new TGraph(npfit,Wcm,XSecElpp);
      grSigPP_Eik->SetLineColor(2);
      grSigPP_Eik->SetLineWidth(1);
      grSigPP_Eik->SetLineStyle(2);
      grSigPP_Eik->Draw("c"); 

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
      leg1->AddEntry(grSigPBARP_Eik,"#bar{p}p","l");
      leg1->AddEntry(grSigPP_Eik,"pp","l");
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
    grSigPBARP_Eik->Draw("csame");
    grSigPP_Eik->Draw("csame");
      
      canv1->SaveAs("XSECel_BSG20orig_w002teste.eps");    
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
    double B=2.0;//par[3];
    /*******************************************************/
    //reading data, then fitting

    //
    // St pp e pbp
    // 
    FILE *aq1,*aq2;//,*aq3,*aq4;
 
    // pbp   
    aq1 = fopen("paw_Stpa.dat","r");
    aq2 = fopen("paw_Stpp.dat","r");
    //aq3 = fopen("SigElPBP.dat","r");
    //aq4 = fopen("SigElPP.dat","r");
    

    const int npStpa = npSigTotPBP;
    double W21[npStpa],Sig21Exp[npStpa],uW21[npStpa],uSig21Exp[npStpa];
   
    for(int i=0;i<npStpa;i++)     
        fscanf(aq1,"%lg %lg %lg %lg",&W21[i],&Sig21Exp[i],&uSig21Exp[i],&uW21[i]);
        
     // pp   
     const int npStpp = npSigTotPP;
     double W11[npStpp],Sig11Exp[npStpp],uW11[npStpp],uSig11Exp[npStpp];
     
     for(int i=0;i<npStpp;i++)   
        fscanf(aq2,"%lg %lg %lg %lg",&W11[i],&Sig11Exp[i],&uSig11Exp[i],&uW11[i]);
/*  
    // pbp 
    const int npSigElpa = npSigElPBP;
    double W22[npSigElpa],Sig22Exp[npSigElpa],uW22[npSigElpa],uSig22Exp[npSigElpa];
   
    for(int i=0;i<npSigElpa;i++)     
        fscanf(aq1,"%lg %lg %lg %lg",&W22[i],&Sig22Exp[i],&uSig22Exp[i],&uW22[i]);
        
     // pp   
     const int npSigElPP = npSigElPP;
     double W12[npSigElPP],Sig12Exp[npSigElPP],uW12[npSigElPP],uSig12Exp[npSigElPP];
     
     for(int i=0;i<npSigElPP;i++)   
        fscanf(aq2,"%lg %lg %lg %lg",&W12[i],&Sig12Exp[i],&uSig12Exp[i],&uW12[i]); 
*/       
     fclose(aq1);
     fclose(aq2);
     //fclose(aq3);
     //fclose(aq4);   

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

void BSG20orig_w002teste()
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
    1.329,
    119.6,
    42.3,
    1.365,
    36.55,
    0.585}; //start values */
    double step[numpar] = {1.e-3,1.e-3,1.e-3,1.e-3,1.e-3,1.e-3}; //steps
    gMinuit->mnparm(0, "w",    vstart[0],   step[0], 0,0,ierflg);
    gMinuit->mnparm(1, "Sig0", vstart[1],   step[1], 0,0,ierflg);
    gMinuit->mnparm(2, "A1",   vstart[2],   step[2], 0,0,ierflg);
    gMinuit->mnparm(3, "del1", vstart[3],   step[3], 0,0,ierflg);
    gMinuit->mnparm(4, "A2",   vstart[4],   step[4], 0,0,ierflg);
    gMinuit->mnparm(5, "del2", vstart[5],   step[5], 0,0,ierflg);
//    gMinuit->mnparm(3, "B",    vstart[3],   step[3], 0,0,ierflg);

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
    double B=2.0;//outpar[3];
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
//----------------------------------------------------------------------
   //calculates total xsection and elastic xsection, creating files to make the plots
    FILE *aq1,*aq2,*aq3,*aq4;
   
    aq1=fopen("Stpp_BSG20orig_w002teste.dat","w"); 
    aq2=fopen("Stpbp_BSG20orig_w002teste.dat","w");
    aq3=fopen("SigElpp_BSG20orig_SigEl_V002.dat","w"); 
    aq4=fopen("SigElpbp_BSG20orig_SigEl_V002.dat","w");
           
     //number of directive computing for GaussLegendreIntegration
    int np = 50;
    double *x=new double[np];
    double *w=new double[np];
    
    double i=0.;
    const double di=0.05;
    double Ecm;
    double Stpp,Stpbarp,Selpp,Selpbarp;
    
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
                     
            fprintf(aq1,"%.2e   %.4lf\n",Ecm,Stpp); 
            fprintf(aq3,"%.2e   %.4lf\n",Ecm,Selpp);
            /*************************************************/
            //pbarp
            SigTotPBARP.SetParameter(9,Ecm);
            SigTotPBARP.CalcGaussLegendreSamplingPoints(np,x,w,1e-15);
            Stpbarp = SigTotPBARP.IntegralFast(np,x,w,bmin,bmax);
            
            SigElPBARP.SetParameter(9,Ecm);
            SigElPBARP.CalcGaussLegendreSamplingPoints(np,x,w,1e-15);
            Selpbarp = SigElPBARP.IntegralFast(np,x,w,bmin,bmax);
            
            fprintf(aq2,"%.2e   %.4lf\n",Ecm,Stpbarp);  
            fprintf(aq4,"%.2e   %.4lf\n",Ecm,Selpbarp);          
            /*************************************************/                             
            i+=di;
            
                   //cout << Ecm << "\t" << Stpp << "\t" << Stpbarp << endl;
            
       }while(Ecm<Wmax); 
       
      fclose(aq1);
      fclose(aq2); 
      fclose(aq3);
      fclose(aq4);
       
      aq1 = fopen("Stpp_BSG20orig_w002teste.dat","r");
      aq2 = fopen("Stpbp_BSG20orig_w002teste.dat","r");
      aq3 = fopen("SigElpp_BSG20orig_SigEl_V002.dat","r");
      aq4 = fopen("SigElpbp_BSG20orig_SigEl_V002.dat","r");
      
      const int npfit=108;
      double Wcm[npfit],XSecTOTpp[npfit],XSecTOTpbp[npfit],XSecElpp[npfit],XSecElpbp[npfit];       
            
      for(int j=0;j<npfit;j++){
          fscanf(aq1,"%lg %lg",&Wcm[j],&XSecTOTpp[j]);
          fscanf(aq2,"%lg %lg",&Wcm[j],&XSecTOTpbp[j]);
          fscanf(aq3,"%lg %lg",&Wcm[j],&XSecElpp[j]);
          fscanf(aq4,"%lg %lg",&Wcm[j],&XSecElpbp[j]);
      }
      
      fclose(aq1);
      fclose(aq2);
      fclose(aq3);
      fclose(aq4);    
       
//----------------------------------------------------------------------
    //Plotting total xsection data         
    PlotTotXSec(Wcm,XSecTOTpp,XSecTOTpbp,npfit);      
    PlotElXSec(Wcm,XSecElpp,XSecElpbp,npfit);
    /******************************************/
    auto stop = high_resolution_clock::now(); 
    auto duration = duration_cast<seconds>(stop - start); 
    cout << "Time(in min): "<<duration.count()/6.e1 << endl;
    /******************************************/ 
        
}
