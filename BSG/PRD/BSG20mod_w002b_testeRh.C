/* PROGRAM: Two-channel QCD minijet-model

   ---> Ref. Original model: PHYSICAL REVIEW D 80, 074014 (2009)
        SciHub link: https://sci-hub.tw/10.1103/PhysRevD.80.074014
        
   Modified version: Good-Walker with two indexes

------------------------------------------------------------------------
   *** Main characteristics ***

	> SigEik = SigSoft + SigPQCD
	> SigSoft = A1*s^{-del1} +/- A2*s^{-del2} + Sig0
	> SigPQCD = B*SigPQCD
	
	> Calculates St pp e pbp

   *** Observations ***

 Obs1: TOTEM + ATLAS dataset
 Obs2: Real parametrization for SigQCD with 10 free-parameters
 Obs3: PDF: CT14
 Obs4: Free: w, Sig0, A1, A2
 Obs5: B is fixed in 2.0
 Obs6: Parâmetros fixos no fit do Paulo (curva central). 

 *** Bandas de incerteza ***
 
 Estratégia 2 -> Essa é a curva central.

   *** Important *** 
   (See: https://arxiv.org/abs/1908.01040 Chapter 8, sections 8.4 and 8.8)
   
   mg=0.4d0 !GeV -> mg=500+/-200 MeV for lmb=300 MeV
 
 --> (MSTW) lbd=0.318920 GeV lambda equivalent to alpha_LO(Mz^2)=0.13939~0.139 macthing-prescrition squeme; m_b=4.18 GeV
 
 --> (MMHT) lbd=0.267471 GeV lambda equivalent to alpha_LO(Mz^2)=0.13499~0.135 macthing-prescrition squeme; m_b=4.18 GeV
 
 --> (CTEQ6L) lbd=0.3260 GeV lambda equivalent to alpha_NLO(Mz^2)=0.118 for 4 flavors
 
 --> (CT14)   lbd=0.3260 GeV lambda equivalent to alpha_NLO(Mz^2)=0.118 for 4 flavors

   *** How to run the code ***
 
 root -l -q /usr/lib/x86_64-linux-gnu/libgsl.so /usr/lib/x86_64-linux-gnu/libgslcblas.so BSG20mod_w002b_testeRh.C | cat > BSG20mod_w002b_testeRh.min
 
   ---------------------------------------------------------------------
   Codes by: M.Broilo*, P.V.R.G.Silva, V.P.B.Gonçalves
   mateus.broilo@ufrgs.br | mateus.broilo90@gmail.com
   High and Medium Energy Group
   Grupo de Altas e Medias Energias (GAME)
   Universidade Federal de Pelotas, Pelotas - RS (Brasil)
   ---------------------------------------------------------------------
   
   Creation: 18/dez/2019 (Pelotas, RS, Brazil)
   Last update: 4/mar/2020 (Pelotas, RS, Brazil)
   
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
const int npRap=12;
const int npRpp=67;
const int npSap=59;
const int npSpp=118;
const int npapfit=npRap+npSap;
const int nppfit=npRpp+npSpp;
const int npMin=npapfit+nppfit;
// Plots
const int npSppPlot=104;
const int npSppTOTEM=12;
const int npSppATLAS=2;
const int npRppPlot=64;
const int npRppTOTEM=4;
//-------------------------------------------------
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
}

double DerSigPQCD(double W)
{    
   double s=W*W;
        
    double Y=Log(s);
    double X=Log(Y);

    double dersig=c1*b2*b3*b4*Exp(b3*Power(X,c1*b4))*Power(X,(c1*b4-1.))/(s*Y)+c2*b5*b6*b7*Exp(b6*Power(X,c2*b7))*Power(X,(c2*b7-1.))/(s*Y)+c3*b8*b9*b10*Exp(b9*Power(X,c3*b10))*Power(X,(c3*b10-1.))/(s*Y);
          
    return dersig;
}
//----------------------------------------------------------------------
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
    double obs=par[7];
    double W=par[8];        
        
    double s=W*W;
    double winv=1./ww;
    double K = 0.;
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
    double nnPP,nnPBP,xxPP,CHFpp,xxPBP,CHFpbp;
    /*******************************************************/
    //Average number of interactions
    /*******************************************************/
	nnPP=(SigEikPP*mb_to_fm2)*Wdip;
	nnPBP=(SigEikPBP*mb_to_fm2)*Wdip;
    /*******************************************************/
    //Confluent Hypergeometric function of second kind
    /*******************************************************/  
	xxPP=2./(nnPP*ww*ww);  
	CHFpp=gsl_sf_hyperg_U(winv,1.0,xxPP);
	xxPBP=2./(nnPBP*ww*ww);  
	CHFpbp=gsl_sf_hyperg_U(winv,1.0,xxPBP);
    /*******************************************************/
    // Expressões para o cálculo de Rh
    	
	double CHFpp2=gsl_sf_hyperg_U(winv+1.,2.,xxPP);
	double CHFpbp2=gsl_sf_hyperg_U(winv+1.,2.,xxPBP);
	double GammaNum = gsl_sf_gamma(winv+1.);
	double GammaDenum = gsl_sf_gamma(winv);
	
	double DerSoftPP = -del1*A1*Power(s/s0,-del1-1.)+del2*A2*Power(s/s0,-del2-1.);
	double DerSoftPBP = -del1*A1*Power(s/s0,-del1-1.)-del2*A2*Power(s/s0,-del2-1.);
		
	double Term1PP = Power(xxPP,winv)*(GammaNum/GammaDenum)*CHFpp2;
	double Term1PBP = Power(xxPBP,winv)*(GammaNum/GammaDenum)*CHFpbp2;
	
	double Term2PP = -CHFpp*winv*Power(xxPP,winv-1.)*(2./(ww*ww))*1./Power(nnPP,2.)*Wdip*(DerSigPQCD(W)+DerSoftPP);
	double Term2PBP = -CHFpbp*winv*Power(xxPBP,winv-1.)*(2./(ww*ww))*1./Power(nnPBP,2.)*Wdip*(DerSigPQCD(W)+DerSoftPBP);
	
	double DerRePerfilPP = b*(Term1PP + Term2PP);
	double DerRePerfilPBP = b*(Term1PBP + Term2PBP);
	
	double SigPP=b*(1.-Power(xxPP,winv)*CHFpp);
	double SigPBP=b*(1.-Power(xxPBP,winv)*CHFpbp);
    
        double ReFeven = K + (Pi()/4.)*(-(SigPP+SigPBP)+s*(DerRePerfilPP+DerRePerfilPBP));
        double ReFodd = (s*Pi()/4.)*(DerRePerfilPP-DerRePerfilPBP);    
    
    /*******************************************************/	
    
    	if(obs==1) //Stpp
    	{
            KerSig=2.*Pi()*b*(2.*(1.-Power(xxPP,winv)*CHFpp))*fm2_to_mb; //in mb 
            return KerSig;
    	}
    	else if(obs==2) //Stpbp
    	{
            KerSig=2.*Pi()*b*(2.*(1.-Power(xxPBP,winv)*CHFpbp))*fm2_to_mb; //in mb 
            return KerSig;
    	}
    	else if(obs==3) //Rhpp
    	{
    	    KerSig=4.*Pi()*(ReFeven+ReFodd);
    	    return KerSig;
    	}
    	else if(obs==4) //Rhpp
    	{
    	    KerSig=4.*Pi()*(ReFeven-ReFodd);
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
void PlotTotXSec(double *Wcm,double *sigtotPPEIK,double *sigtotPBARPEIK,int npfit)
{
    FILE *aq1,*aq2,*aq3,*aq4;
// Plotting St data     
    aq1 = fopen("paw_Stpa.dat","r");
    aq2 = fopen("paw_Stpp.dat","r");
    aq3 = fopen("paw_StTOTEM.dat","r");
    aq4 = fopen("paw_StATLAS.dat","r");

    // pbp
    const int npdpa = npSap;
    double Wpap[npdpa],SigExpPaP[npdpa],uWpap[npdpa],uSigExpPaP[npdpa];
   
    for(int i=0;i<npdpa;i++)     
        fscanf(aq1,"%lg %lg %lg %lg",&Wpap[i],&SigExpPaP[i],&uSigExpPaP[i],&uWpap[i]);
        
     // pp     
     const int npdpp = npSppPlot;
     double Wpp[npdpp],SigExpPP[npdpp],uWpp[npdpp],uSigExpPP[npdpp];
     
     for(int i=0;i<npdpp;i++)   
        fscanf(aq2,"%lg %lg %lg %lg",&Wpp[i],&SigExpPP[i],&uSigExpPP[i],&uWpp[i]);
       
       //---> TOTEM data <---
       
     const int npTOTEM = npSppTOTEM;
     double Wpp2[npTOTEM],SigExpPP2[npTOTEM],uWpp2[npTOTEM],uSigExpPP2[npTOTEM];  
       
     for(int i=0;i<npTOTEM;i++)   
        fscanf(aq3,"%lg %lg %lg %lg",&Wpp2[i],&SigExpPP2[i],&uSigExpPP2[i],&uWpp2[i]);
        
       //---> ATLAS data <---
       
     const int npATLAS = npSppATLAS;
     double Wpp3[npATLAS],SigExpPP3[npATLAS],uWpp3[npATLAS],uSigExpPP3[npATLAS];  
       
     for(int i=0;i<npATLAS;i++)   
        fscanf(aq4,"%lg %lg %lg %lg",&Wpp3[i],&SigExpPP3[i],&uSigExpPP3[i],&uWpp3[i]);
      
     fclose(aq1);
     fclose(aq2);
     fclose(aq3);
     fclose(aq4);
     
      TCanvas *canv1 = new TCanvas("c1","Total Cross Section",200,10,900,600);
      canv1->SetFillColor(0);
      canv1->SetLogx();
       
       double x_min(5e00),x_max(2.e05);
       double y_min(3.0e01),y_max(1.8e02);
      
      TGraphErrors *gr0 = new TGraphErrors(npdpa,Wpap,SigExpPaP,uWpap,uSigExpPaP);
      gr0->SetMarkerStyle(24);
      gr0->SetMarkerColor(1);
      gr0->SetLineColor(1);
      gr0->SetMarkerSize(1.0);
      gr0->SetTitle();
      gr0->Draw("apz");    
           
      TGraphErrors *gr1 = new TGraphErrors(npdpp,Wpp,SigExpPP,uWpp,uSigExpPP);
      gr1->SetMarkerStyle(20);
      gr1->SetMarkerColor(1);
      gr1->SetLineColor(1);
      gr1->SetMarkerSize(1.0);
      gr1->Draw("pz");         
      
      TGraphErrors *gr2 = new TGraphErrors(npTOTEM,Wpp2,SigExpPP2,uWpp2,uSigExpPP2);
      gr2->SetMarkerStyle(33);
      gr2->SetMarkerColor(2);
      gr2->SetLineColor(2);
      gr2->SetMarkerSize(1.3);
      gr2->Draw("pz");  
      
      TGraphErrors *gr3 = new TGraphErrors(npATLAS,Wpp3,SigExpPP3,uWpp3,uSigExpPP3);
      gr3->SetMarkerStyle(21);
      gr3->SetMarkerColor(4);
      gr3->SetLineColor(4);
      gr3->SetMarkerSize(0.8);
      gr3->Draw("pz");
      
      TGraph *grSigPBARP_Eik = new TGraph(npfit,Wcm,sigtotPBARPEIK);
      grSigPBARP_Eik->SetLineColor(4);
      grSigPBARP_Eik->SetLineWidth(1);
      grSigPBARP_Eik->SetLineStyle(1);
      grSigPBARP_Eik->Draw("c");
      
      TGraph *grSigPP_Eik = new TGraph(npfit,Wcm,sigtotPPEIK);
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
      
      canv1->SaveAs("St_BSG20mod_w002b_testeRh.eps");    
}

/*******************************************************************************/
// Plotando o parâmetro Rho pp e pbp
/*******************************************************************************/
void PlotRho(double *Wcm,double *rhoPPEIK,double *rhoPBARPEIK,int npfit)
{
    FILE *aq1,*aq2,*aq3;   
// Plotting rho data        
    aq1 = fopen("paw_Rhpa.dat","r");
    aq2 = fopen("paw_Rhpp.dat","r");
    aq3 = fopen("paw_RhTOTEM.dat","r");
    
    // pbp
    const int nr1 = npRap;
    double Er1[nr1],uEr1[nr1],rho1[nr1],urho1[nr1];
   
    for(int i=0;i<nr1;i++)
       fscanf(aq1,"%lg %lg %lg %lg",&Er1[i],&rho1[i],&urho1[i],&uEr1[i]);
      
    // pp 
    const int nr2 = npRppPlot;
    double Er2[nr2],uEr2[nr2],rho2[nr2],urho2[nr2];
   
    for(int i=0;i<nr2;i++)
       fscanf(aq2,"%lg %lg %lg %lg",&Er2[i],&rho2[i],&urho2[i],&uEr2[i]);
     
       //---> TOTEM data <---
       
     const int npTOTEM = npRppTOTEM;
     double Er3[npTOTEM],uEr3[npTOTEM],rho3[npTOTEM],urho3[npTOTEM];  
       
     for(int i=0;i<npTOTEM;i++)   
        fscanf(aq3,"%lg %lg %lg %lg",&Er3[i],&rho3[i],&urho3[i],&uEr3[i]);    
     
     fclose(aq1);
     fclose(aq2);
     fclose(aq3);
    
     TCanvas *canv2 = new TCanvas("c1","Rho parameter",200,10,900,600);
      canv2->SetFillColor(0);
      canv2->SetLogx();
//       c1->SetGrid();
       
      double x_min(5e0),x_max(2.e4);
      double y_min(-0.2),y_max(0.2);    
       
      TGraphErrors *gr00 = new TGraphErrors(nr1,Er1,rho1,uEr1,urho1);
      gr00->SetMarkerStyle(24);
      gr00->SetMarkerColor(1);
      gr00->SetLineColor(1);
      gr00->SetMarkerSize(1.0);
      gr00->SetTitle();
      gr00->Draw("apz");     
     
      TGraphErrors *gr11= new TGraphErrors(nr2,Er2,rho2,uEr2,urho2);
      gr11->SetMarkerStyle(20);
      gr11->SetMarkerColor(1);
      gr11->SetLineColor(1);
      gr11->SetMarkerSize(1.0);
      gr11->SetTitle();
      gr11->Draw("pz"); 
      
      TGraphErrors *gr22 = new TGraphErrors(npTOTEM,Er3,rho3,uEr3,urho3);
      gr22->SetMarkerStyle(33);
      gr22->SetMarkerColor(2);
      gr22->SetLineColor(2);
      gr22->SetMarkerSize(1.3);
      gr22->Draw("pz");

      TGraph *grRhoPBARP_EIK = new TGraph(npfit,Wcm,rhoPBARPEIK);
      grRhoPBARP_EIK->SetLineColor(4);
      grRhoPBARP_EIK->SetLineWidth(1);
      grRhoPBARP_EIK->SetLineStyle(1);
      grRhoPBARP_EIK->Draw("c");
      
      TGraph *grRhoPP_EIK = new TGraph(npfit,Wcm,rhoPPEIK);
      grRhoPP_EIK->SetLineColor(2);
      grRhoPP_EIK->SetLineWidth(1);
      grRhoPP_EIK->SetLineStyle(2);
      grRhoPP_EIK->Draw("c");  
      
      gr00->GetYaxis()->SetTitle("#bf{#rho}");
      gr00->GetYaxis()->SetRangeUser(y_min,y_max);
      gr00->GetYaxis()->SetTitleOffset(1.2);
      gr00->GetYaxis()->SetTitleOffset(1.2);
      gr00->GetYaxis()->SetTitleFont(42);
      gr00->GetXaxis()->SetTitle("#bf{#sqrt{s} [GeV]}");
      gr00->GetXaxis()->SetTitleOffset(1.3);
      gr00->GetXaxis()->SetLimits(x_min,x_max);   
      gr00->GetXaxis()->SetTitleFont(42);           
     
//      TLegend *leg2 = new TLegend(0.15,0.65,0.42,0.89);
      TLegend *leg2 = new TLegend(0.65,0.15,0.89,0.42);  
      leg2->AddEntry(gr00,"#bar{p}p data ","p");
      leg2->AddEntry(gr11,"pp data below LHC","p");
      leg2->AddEntry(gr22,"TOTEM","p");
      leg2->AddEntry(grRhoPBARP_EIK,"#bar{p}p","l");
      leg2->AddEntry(grRhoPP_EIK,"pp","l");
      leg2->SetShadowColor(0);
      leg2->SetBorderSize(0);
      leg2->SetTextSize(0.035);
      leg2->Draw();     

/* ---> desenhando um miniframe em um pad dentro do canvas canv2 <--- */
/*
    TPad *pad1 = new TPad("pad1","This is pad1",0.51,0.15,0.89,0.55);
    pad1->SetFillColor(0);
    pad1->SetLogx();
    pad1->Draw();     
    pad1->cd();
    
    double xi = 5e3;
    double xf = 15e3;
    double yi = 0.07;
    double yf = 0.16;

    gr22->Draw("apz");
    gr22->SetTitle("");
    gr22->GetYaxis()->SetRangeUser(yi,yf);
    gr22->GetYaxis()->SetLabelSize(0.06);
    gr22->GetXaxis()->SetLabelSize(0.06);
    gr22->GetYaxis()->SetTitleFont(42);
    gr22->GetXaxis()->SetRangeUser(xi,xf); 
    grRhoPBARP_EIK->Draw("csame");
    grRhoPP_EIK->Draw("csame");
*/           
      canv2->SaveAs("Rh_BSG20mod_w002b_testeRh.eps");      
      canv2->Clear();
}

//----------------------------------------------------------------------
void fcn(int &npar, double *gin, double &f, double *par, int iflag)
{    
   /*******************************************************/ 
    //fit parameters
    double ww=1.72;//par[0]; 
    double Sig0=146.6;//par[1];
    double A1=54.6;//par[2];
    double del1=1.39;//par[3];
    double A2=47.5;//par[4];  
    double del2=0.588;//par[5];
    double B=2.0;//par[3];
    /*******************************************************/
    //reading data, then fitting

    //
    // St pp e pbp
    // 
    FILE *aq1,*aq2;
 
    // pbp   
    aq1 = fopen("paw_Stpa.dat","r");
    aq2 = fopen("paw_Stpp.dat","r");

    const int npdpa = npSap;
    double Wpap[npdpa],SigExpPaP[npdpa],uWpap[npdpa],uSigExpPaP[npdpa];
   
    for(int i=0;i<npdpa;i++)     
        fscanf(aq1,"%lg %lg %lg %lg",&Wpap[i],&SigExpPaP[i],&uSigExpPaP[i],&uWpap[i]);
        
     // pp   
     const int npdpp = npSpp;
     double Wpp[npdpp],SigExpPP[npdpp],uWpp[npdpp],uSigExpPP[npdpp];
     
     for(int i=0;i<npdpp;i++)   
        fscanf(aq2,"%lg %lg %lg %lg",&Wpp[i],&SigExpPP[i],&uSigExpPP[i],&uWpp[i]);
       
     fclose(aq1);
     fclose(aq2);
     
    //
    // Rh pp e pbp
    //             
    aq1 = fopen("paw_Rhpa.dat","r");
    aq2 = fopen("paw_Rhpp.dat","r");
    
    //pbp
    const int nr1 = npRap;
    double Er1[nr1],uEr1[nr1],rho1[nr1],urho1[nr1];
   
    for(int i=0;i<nr1;i++)
       fscanf(aq1,"%lg %lg %lg %lg",&Er1[i],&rho1[i],&urho1[i],&uEr1[i]);
      
    //pp 
    const int nr2 = npRpp;
    double Er2[nr2],uEr2[nr2],rho2[nr2],urho2[nr2];
   
    for(int i=0;i<nr2;i++)
       fscanf(aq2,"%lg %lg %lg %lg",&Er2[i],&rho2[i],&urho2[i],&uEr2[i]);
     
     fclose(aq1);
     fclose(aq2);   

    /*******************************************************/
    //Total XSection - pp
       
    TF1 SigTotPP("SigTotPP",KerForward,bmin,bmax,9);
    SigTotPP.SetParameter(0,ww);
    SigTotPP.SetParameter(1,Sig0);
    SigTotPP.SetParameter(2,A1);
    SigTotPP.SetParameter(3,del1);
    SigTotPP.SetParameter(4,A2);
    SigTotPP.SetParameter(5,del2);
    SigTotPP.SetParameter(6,B);
    SigTotPP.SetParameter(7,1);    
    
    /*******************************************************/
    //Total XSection - pbp
       
    TF1 SigTotPBARP("SigTotPBARP",KerForward,bmin,bmax,9);
    SigTotPBARP.SetParameter(0,ww);
    SigTotPBARP.SetParameter(1,Sig0);
    SigTotPBARP.SetParameter(2,A1);
    SigTotPBARP.SetParameter(3,del1);
    SigTotPBARP.SetParameter(4,A2);
    SigTotPBARP.SetParameter(5,del2);
    SigTotPBARP.SetParameter(6,B);
    SigTotPBARP.SetParameter(7,2);

    /*******************************************************/
    //Rho-parameter - pp
       
    TF1 RhoPP("RhoPP",KerForward,bmin,bmax,9);
    RhoPP.SetParameter(0,ww);
    RhoPP.SetParameter(1,Sig0);
    RhoPP.SetParameter(2,A1);
    RhoPP.SetParameter(3,del1);
    RhoPP.SetParameter(4,A2);
    RhoPP.SetParameter(5,del2);
    RhoPP.SetParameter(6,B);
    RhoPP.SetParameter(7,3);    
    
    /*******************************************************/
    //Rho-parameter - pbp
       
    TF1 RhoPBARP("RhoPBARP",KerForward,bmin,bmax,9);
    RhoPBARP.SetParameter(0,ww);
    RhoPBARP.SetParameter(1,Sig0);
    RhoPBARP.SetParameter(2,A1);
    RhoPBARP.SetParameter(3,del1);
    RhoPBARP.SetParameter(4,A2);
    RhoPBARP.SetParameter(5,del2);
    RhoPBARP.SetParameter(6,B);
    RhoPBARP.SetParameter(7,4);
        
    /*******************************************************/
    
    //Chi-squared calculation
    int np = 50;
    double *x=new double[np];
    double *w=new double[np];
    double sigpp,sigpbarp;    
    
    double delta;
    double chisq=0.;
    
    for(int i=0;i<npdpp;i++){
        SigTotPP.SetParameter(8,Wpp[i]);
        SigTotPP.CalcGaussLegendreSamplingPoints(np,x,w,1e-15);
        sigpp=SigTotPP.IntegralFast(np,x,w,bmin,bmax);
        delta=(SigExpPP[i]-sigpp)/uSigExpPP[i];
        chisq+=Power(delta,2);          
   }
       
    delta=0.;

    for(int i=0;i<npdpa;i++){
        SigTotPBARP.SetParameter(8,Wpap[i]);
        SigTotPBARP.CalcGaussLegendreSamplingPoints(np,x,w,1e-15);
        sigpbarp=SigTotPBARP.IntegralFast(np,x,w,bmin,bmax);
        delta=(SigExpPaP[i]-sigpbarp)/uSigExpPaP[i];
        chisq+=Power(delta,2);           
   }
   
    delta=0.; 
       
   f = chisq;        
}
//----------------------------------------------------------------------

void BSG20mod_w002b_testeRh()
{
    /******************************************/
    //time control - start
    auto start = high_resolution_clock::now();    
    /******************************************/ 
         
    TMinuit *gMinuit = new TMinuit(numpar);  //initialize TMinuit with a maximum of 4 params
    gMinuit->SetFCN(fcn);
    
    double arglist[numpar];
    int ierflg = 0;
  
    arglist[0] = 7.04;
    gMinuit->mnexcm("SET ERR",arglist,1,ierflg); 
      
    double vstart[numpar] = {
    1.72,
    146.6,
    54.6,
    1.39,
    47.5,
    0.588
    }; //start values
    double step[numpar] = {1.e-3,1.e-3,1.e-3,1.e-3,1.e-3,1.e-3}; //steps
//    gMinuit->mnparm(0, "w",    vstart[0],   step[0], 0,0,ierflg);
//    gMinuit->mnparm(1, "Sig0", vstart[1],   step[1], 0,0,ierflg);
//    gMinuit->mnparm(2, "A1",   vstart[2],   step[2], 0,0,ierflg);
//    gMinuit->mnparm(3, "del1", vstart[3],   step[3], 0,0,ierflg);
//    gMinuit->mnparm(4, "A2",   vstart[4],   step[4], 0,0,ierflg);
//    gMinuit->mnparm(5, "del2", vstart[5],   step[5], 0,0,ierflg);
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
    double ww=1.72;//outpar[0]; 
    double Sig0=146.6;//outpar[1];
    double A1=54.6;//outpar[2];
    double del1=1.39;//outpar[3];
    double A2=47.5;//outpar[4];  
    double del2=0.588;//outpar[5];
    double B=2.0;//outpar[3];     
    /*******************************************************/
    //Total XSection - pp
       
    TF1 SigTotPP("SigTotPP",KerForward,bmin,bmax,9);
    SigTotPP.SetParameter(0,ww);
    SigTotPP.SetParameter(1,Sig0);
    SigTotPP.SetParameter(2,A1);
    SigTotPP.SetParameter(3,del1);
    SigTotPP.SetParameter(4,A2);
    SigTotPP.SetParameter(5,del2);
    SigTotPP.SetParameter(6,B);
    SigTotPP.SetParameter(7,1);    
    
    /*******************************************************/
    //Total XSection - pbp
       
    TF1 SigTotPBARP("SigTotPBARP",KerForward,bmin,bmax,9);
    SigTotPBARP.SetParameter(0,ww);
    SigTotPBARP.SetParameter(1,Sig0);
    SigTotPBARP.SetParameter(2,A1);
    SigTotPBARP.SetParameter(3,del1);
    SigTotPBARP.SetParameter(4,A2);
    SigTotPBARP.SetParameter(5,del2);
    SigTotPBARP.SetParameter(6,B);
    SigTotPBARP.SetParameter(7,2);

    /*******************************************************/
    //Rho-parameter - pp
       
    TF1 RhoPP("RhoPP",KerForward,bmin,bmax,9);
    RhoPP.SetParameter(0,ww);
    RhoPP.SetParameter(1,Sig0);
    RhoPP.SetParameter(2,A1);
    RhoPP.SetParameter(3,del1);
    RhoPP.SetParameter(4,A2);
    RhoPP.SetParameter(5,del2);
    RhoPP.SetParameter(6,B);
    RhoPP.SetParameter(7,3);    
    
    /*******************************************************/
    //Rho-parameter - pbp
       
    TF1 RhoPBARP("RhoPBARP",KerForward,bmin,bmax,9);
    RhoPBARP.SetParameter(0,ww);
    RhoPBARP.SetParameter(1,Sig0);
    RhoPBARP.SetParameter(2,A1);
    RhoPBARP.SetParameter(3,del1);
    RhoPBARP.SetParameter(4,A2);
    RhoPBARP.SetParameter(5,del2);
    RhoPBARP.SetParameter(6,B);
    RhoPBARP.SetParameter(7,4);
        
    /*******************************************************/  
//----------------------------------------------------------------------
   //calculates total xsection and rho, creating files to make the plots
    FILE *aq1,*aq2,*aq3,*aq4;

    aq1=fopen("Stpp_BSG20mod_w002b_testeRh.dat","w"); 
    aq2=fopen("Rhpp_BSG20mod_w002b_testeRh.dat","w"); 
    aq3=fopen("Stpbp_BSG20mod_w002b_testeRh.dat","w");
    aq4=fopen("Rhpbp_BSG20mod_w002b_testeRh.dat","w"); 
           
     //number of directive computing for GaussLegendreIntegration
    int np = 50;
    double *x=new double[np];
    double *w=new double[np];
    
    double i=0.;
    const double di=0.05;
    double Ecm;
    double sigpp,rhopp,sigpbarp,rhopbarp;
    
     do{
            Ecm=5.*pow(10.,i);
            /*************************************************/
            //pp
            SigTotPP.SetParameter(8,Ecm);
            RhoPP.SetParameter(8,Ecm);
            SigTotPP.CalcGaussLegendreSamplingPoints(np,x,w,1e-15);
            RhoPP.CalcGaussLegendreSamplingPoints(np,x,w,1e-15);
            sigpp = SigTotPP.IntegralFast(np,x,w,bmin,bmax);
            rhopp = RhoPP.IntegralFast(np,x,w,bmin,bmax);
            rhopp=rhopp/sigpp;
            fprintf(aq1,"%.2e   %.4lf\n",Ecm,sigpp);             
            fprintf(aq2,"%.2e   %.4lf\n",Ecm,rhopp);
            /*************************************************/
            //pbarp
            SigTotPBARP.SetParameter(8,Ecm);
            RhoPBARP.SetParameter(8,Ecm);
            SigTotPBARP.CalcGaussLegendreSamplingPoints(np,x,w,1e-15);
            RhoPBARP.CalcGaussLegendreSamplingPoints(np,x,w,1e-15);
            sigpbarp = SigTotPBARP.IntegralFast(np,x,w,bmin,bmax);
            rhopbarp = RhoPBARP.IntegralFast(np,x,w,bmin,bmax);
            rhopbarp=rhopbarp/sigpbarp;
            fprintf(aq3,"%.2e   %.4lf\n",Ecm,sigpbarp);             
            fprintf(aq4,"%.2e   %.4lf\n",Ecm,rhopbarp);            
            /*************************************************/                            
            i+=di;
            
                   cout << Ecm << "\t" << rhopp << "\t" << rhopbarp << endl;
            
       }while(Ecm<Wmax); 

      fclose(aq1);
      fclose(aq2);
      fclose(aq3);
      fclose(aq4); 
       
    aq1=fopen("Stpp_BSG20mod_w002b_testeRh.dat","r"); 
    aq2=fopen("Rhpp_BSG20mod_w002b_testeRh.dat","r"); 
    aq3=fopen("Stpbp_BSG20mod_w002b_testeRh.dat","r");
    aq4=fopen("Rhpbp_BSG20mod_w002b_testeRh.dat","r");
      
      const int npfit=108;
      double Wcm[npfit],sigtotPPEIK[npfit],sigtotPBARPEIK[npfit],rhoPPEIK[npfit],rhoPBARPEIK[npfit];     
            
      for(int j=0;j<npfit;j++){
          fscanf(aq1,"%lg %lg",&Wcm[j],&sigtotPPEIK[j]);
          fscanf(aq2,"%lg %lg",&Wcm[j],&rhoPPEIK[j]);
          fscanf(aq3,"%lg %lg",&Wcm[j],&sigtotPBARPEIK[j]);
          fscanf(aq4,"%lg %lg",&Wcm[j],&rhoPBARPEIK[j]);
      }
      
      fclose(aq1);
      fclose(aq2);
      fclose(aq3);
      fclose(aq4);
       
//----------------------------------------------------------------------
    //Plotting total xsection data         
    PlotTotXSec(Wcm,sigtotPPEIK,sigtotPBARPEIK,npfit);   
    
    //Plotting rho-parameter data         
    PlotRho(Wcm,rhoPPEIK,rhoPBARPEIK,npfit);   
   
    /******************************************/
    auto stop = high_resolution_clock::now(); 
    auto duration = duration_cast<seconds>(stop - start); 
    cout << "Time(in min): "<<duration.count()/6.e1 << endl;
    /******************************************/ 
        
}
