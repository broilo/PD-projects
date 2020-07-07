/* PROGRAM: Two-channel QCD minijet-model

   ---> Ref. Original model: PHYSICAL REVIEW D 80, 074014 (2009)
        SciHub link: https://sci-hub.tw/10.1103/PhysRevD.80.074014
        
   Modified version: Good-Walker with two indexes

------------------------------------------------------------------------
   *** Main characteristics ***

	> SigEik = SigSoft + SigPQCD
	> SigSoft = A1*s^{-del1} +/- A2*s^{-del2} + Sig0
	> SigPQCD = B*SigPQCD
	
	> Calculates SigSD e SigDD pp e pbp

   *** Observations ***

 Obs1: TOTEM + ATLAS dataset
 Obs2: Real parametrization for SigQCD with 10 free-parameters
 Obs3: PDF: CT14
 Obs4: Free: w, Sig0, A1, A2
 Obs5: B is fixed in 2.0
 Obs6: Parâmetros fixos no fit do Paulo (curva central). 

 *** Bandas de incerteza ***
 
 Estratégia 1 -> Manter os parâmetros nos valores ajustados e estudar a incerteza em w.
 
 Intervalo w=[1.3,(1.9)]

   *** Important *** 
   (See: https://arxiv.org/abs/1908.01040 Chapter 8, sections 8.4 and 8.8)
   
   mg=0.4d0 !GeV -> mg=500+/-200 MeV for lmb=300 MeV
 
 --> (MSTW) lbd=0.318920 GeV lambda equivalent to alpha_LO(Mz^2)=0.13939~0.139 macthing-prescrition squeme; m_b=4.18 GeV
 
 --> (MMHT) lbd=0.267471 GeV lambda equivalent to alpha_LO(Mz^2)=0.13499~0.135 macthing-prescrition squeme; m_b=4.18 GeV
 
 --> (CTEQ6L) lbd=0.3260 GeV lambda equivalent to alpha_NLO(Mz^2)=0.118 for 4 flavors
 
 --> (CT14)   lbd=0.3260 GeV lambda equivalent to alpha_NLO(Mz^2)=0.118 for 4 flavors

   *** How to run the code ***
 
 root -l -q /usr/lib/x86_64-linux-gnu/libgsl.so /usr/lib/x86_64-linux-gnu/libgslcblas.so BSG20mod_Diffs_V002zc.C | cat > BSG20mod_Diffs_V002zc.min
 
   ---------------------------------------------------------------------
   Codes by: M.Broilo*, P.V.R.G.Silva, V.P.B.Gonçalves
   mateus.broilo@ufrgs.br | mateus.broilo90@gmail.com
   High and Medium Energy Group
   Grupo de Altas e Medias Energias (GAME)
   Universidade Federal de Pelotas, Pelotas - RS (Brasil)
   ---------------------------------------------------------------------
   
   Creation: 18/fev/2020 (Pelotas, RS, Brazil)
   Last update: 18/fev/2020 (Pelotas, RS, Brazil)
   
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
const int numpar=4;
const int npSigSDpp=5.;
const int npSigSDpbp=16.;
const int npapfit=npSigSDpbp;
const int nppfit=npSigSDpp;
const int npMin=npapfit+nppfit;
// Plots
const int npSigSD1=1.;
const int npSigSD2=1.;
const int npSigSD3=3.;
const int npSigSD4=8.;
const int npSigSD5=3.;
const int npSigSD6=2.;
const int npSigSD7=2.;
const int npSigSD8=1.;

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
    double nn,xx,CHF;
    
    	if(reac==1) //pp
    	{
    	/*******************************************************/
    	//Average number of interactions
    	/*******************************************************/
	nn=(SigEikPP*mb_to_fm2)*Wdip;
    	}
    	else if(reac==2) //pbp
    	{
	nn=(SigEikPBP*mb_to_fm2)*Wdip;
    	}
    	else
    	{
    	nn=0.;
    	}
    	/*******************************************************/
	//Confluent Hypergeometric function of second kind
	/*******************************************************/  
	xx=2./(nn*ww*ww);  	
	double CHFa=gsl_sf_hyperg_U(winv,1.0,xx);
	double CHFb=gsl_sf_hyperg_U(winv,1.0-winv,xx);
	double CHFc=gsl_sf_hyperg_U(winv,1.0,xx/2.);
    	/*******************************************************/
	//Mean values over Target and Projectile
	/*******************************************************/  	
	double avg_braTket2_PT = 1.-2.*Power(xx,winv)*CHFa + Power(xx,winv)*CHFb;	
	double avg_braTket_PT = 1. - Power(xx,winv)*CHFa;	
	double avg_braT2ket_PT = 1. - 2.*Power(xx,winv)*CHFa + Power(xx/2.,winv)*CHFc;
	/*******************************************************/
	if(obs==1) //SigSD
	{
            KerSig=2.*Pi()*b*(2.*avg_braTket2_PT-2.*Power(avg_braTket_PT,2.))*fm2_to_mb; //in mb 
        }
        else if(obs==2) //SigDD
        {
            KerSig=2.*Pi()*b*(avg_braT2ket_PT-2.*avg_braTket2_PT+Power(avg_braTket_PT,2.))*fm2_to_mb; //in mb
        }        
        else
        {
            KerSig=0;
        }
            return KerSig;    
}
/*
	SigSD = <<t>^{2}_{T}>_{P} + <<t>^{2}_{P}>_{T} - 2*<t>^{2}_{PT}
        SigDD = <t^{2}>_{PT} - <<t>^{2}_{T}>_{P} - <<t>^{2}_{P}>_{T} + <t>^{2}_{PT}

 <<t>^{2}_{T}>_{P} = <<t>^{2}_{P}>_{T} = 1 - 2*x**(1/w)*U(1/w,1,x) + x**(1/w)*U(1/w,1-(1/w),x)

 <t>_{PT} = 1 - x**(1/w)*U(1/w,1,x)
 
 <t^{2}>_{PT} = 1 - 2*x**(1/w)*U(1/w,1,x) + (x/2)**(1/w)*U(1/w,1,x/2)

*/

//----------------------------------------------------------------------
// Plotando a seção de choque total elástica pp e pbp
//----------------------------------------------------------------------
void PlotTotXSec(double *Wcm,double *sigtotsdPP,double *sigtotsdPBARP,double *sigtotddPP,double *sigtotddPBARP,int npfit)
{
    FILE *aq1,*aq2,*aq3,*aq4,*aq5,*aq6,*aq7,*aq8;
// Plotting St data    
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
    const int npdpa4 = npSigSD4;
    double Wpap4[npdpa4],SigExpPaP4[npdpa4],uWpap4[npdpa4],uSigExpPaP4[npdpa4];
   
    for(int i=0;i<npdpa4;i++)     
        fscanf(aq4,"%lg %lg %lg %lg",&Wpap4[i],&SigExpPaP4[i],&uSigExpPaP4[i],&uWpap4[i]);
        
    const int npdpa5 = npSigSD5;
    double Wpap5[npdpa5],SigExpPaP5[npdpa5],uWpap5[npdpa5],uSigExpPaP5[npdpa5];
   
    for(int i=0;i<npdpa5;i++)     
        fscanf(aq5,"%lg %lg %lg %lg",&Wpap5[i],&SigExpPaP5[i],&uSigExpPaP5[i],&uWpap5[i]);
        
    const int npdpa6 = npSigSD6;
    double Wpap6[npdpa6],SigExpPaP6[npdpa6],uWpap6[npdpa6],uSigExpPaP6[npdpa6];
   
    for(int i=0;i<npdpa6;i++)     
        fscanf(aq6,"%lg %lg %lg %lg",&Wpap6[i],&SigExpPaP6[i],&uSigExpPaP6[i],&uWpap6[i]);        

    const int npdpa7 = npSigSD7;
    double Wpap7[npdpa7],SigExpPaP7[npdpa7],uWpap7[npdpa7],uSigExpPaP7[npdpa7];
   
    for(int i=0;i<npdpa7;i++)     
        fscanf(aq7,"%lg %lg %lg %lg",&Wpap7[i],&SigExpPaP7[i],&uSigExpPaP7[i],&uWpap7[i]);
       
    const int npdpa8 = npSigSD8;
    double Wpap8[npdpa8],SigExpPaP8[npdpa8],uWpap8[npdpa8],uSigExpPaP8[npdpa8];
   
    for(int i=0;i<npdpa8;i++)     
        fscanf(aq8,"%lg %lg %lg %lg",&Wpap8[i],&SigExpPaP8[i],&uSigExpPaP8[i],&uWpap8[i]);       
       
     // pp     
     const int npdpp1 = npSigSD1;
     double Wpp1[npdpp1],SigExpPP1[npdpp1],uWpp1[npdpp1],uSigExpPP1[npdpp1];
     
     for(int i=0;i<npdpp1;i++)   
        fscanf(aq1,"%lg %lg %lg %lg",&Wpp1[i],&SigExpPP1[i],&uSigExpPP1[i],&uWpp1[i]);

     const int npdpp2 = npSigSD2;
     double Wpp2[npdpp2],SigExpPP2[npdpp2],uWpp2[npdpp2],uSigExpPP2[npdpp2];
     
     for(int i=0;i<npdpp2;i++)   
        fscanf(aq2,"%lg %lg %lg %lg",&Wpp2[i],&SigExpPP2[i],&uSigExpPP2[i],&uWpp2[i]);

     const int npdpp3 = npSigSD3;
     double Wpp3[npdpp3],SigExpPP3[npdpp3],uWpp3[npdpp3],uSigExpPP3[npdpp3];
     
     for(int i=0;i<npdpp3;i++)   
        fscanf(aq3,"%lg %lg %lg %lg",&Wpp3[i],&SigExpPP3[i],&uSigExpPP3[i],&uWpp3[i]);
      
     fclose(aq1);
     fclose(aq2);
     fclose(aq3);
     fclose(aq4);
     fclose(aq5);
     fclose(aq6);
     fclose(aq7);
     fclose(aq8);
     
      TCanvas *canv1 = new TCanvas("c1","Total SD and DD Cross Section",200,10,900,600);
      canv1->SetFillColor(0);
      canv1->SetLogx();
       
       double x_min(5e00),x_max(2.e05);
       double y_min(0.0),y_max(4.0e01);
              
      TGraphErrors *gr0 = new TGraphErrors(npdpa4,Wpap4,SigExpPaP4,uWpap4,uSigExpPaP4);
      gr0->SetMarkerStyle(20);
      gr0->SetMarkerColor(1);
      gr0->SetLineColor(1);
      gr0->SetMarkerSize(1.0);
      gr0->SetTitle();
      gr0->Draw("apz");
      
      TGraphErrors *gr1 = new TGraphErrors(npdpa5,Wpap5,SigExpPaP5,uWpap5,uSigExpPaP5);
      gr1->SetMarkerStyle(26);
      gr1->SetMarkerColor(4);
      gr1->SetLineColor(1);
      gr1->SetMarkerSize(1.0);
      gr1->SetTitle();
      gr1->Draw("apz");
      
      TGraphErrors *gr6 = new TGraphErrors(npdpa6,Wpap6,SigExpPaP6,uWpap6,uSigExpPaP6);
      gr6->SetMarkerStyle(28);
      gr6->SetMarkerColor(4);
      gr6->SetLineColor(1);
      gr6->SetMarkerSize(1.0);
      gr6->SetTitle();
      gr6->Draw("apz");
      
      TGraphErrors *gr7 = new TGraphErrors(npdpa7,Wpap7,SigExpPaP7,uWpap7,uSigExpPaP7);
      gr7->SetMarkerStyle(24);
      gr7->SetMarkerColor(4);
      gr7->SetLineColor(1);
      gr7->SetMarkerSize(1.0);
      gr7->SetTitle();
      gr7->Draw("apz");
      
      TGraphErrors *gr8 = new TGraphErrors(npdpa8,Wpap8,SigExpPaP8,uWpap8,uSigExpPaP8);
      gr8->SetMarkerStyle(25);
      gr8->SetMarkerColor(4);
      gr8->SetLineColor(1);
      gr8->SetMarkerSize(1.0);
      gr8->SetTitle();
      gr8->Draw("apz");                                
           
      TGraphErrors *gr1 = new TGraphErrors(npdpp1,Wpp1,SigExpPP1,uWpp1,uSigExpPP1);
      gr1->SetMarkerStyle(29);
      gr1->SetMarkerColor(2);
      gr1->SetLineColor(1);
      gr1->SetMarkerSize(1.0);
      gr1->Draw("pz");         
      
      TGraphErrors *gr2 = new TGraphErrors(npdpp2,Wpp2,SigExpPP2,uWpp2,uSigExpPP2);
      gr2->SetMarkerStyle(21);
      gr2->SetMarkerColor(2);
      gr2->SetLineColor(1);
      gr2->SetMarkerSize(1.0);
      gr2->Draw("pz"); 
      
      TGraphErrors *gr3 = new TGraphErrors(npdpp3,Wpp3,SigExpPP3,uWpp3,uSigExpPP3);
      gr3->SetMarkerStyle(22);
      gr3->SetMarkerColor(2);
      gr3->SetLineColor(1);
      gr3->SetMarkerSize(1.0);
      gr3->Draw("pz");       
      
      TGraph *grSigSDPBARP = new TGraph(npfit,Wcm,sigtotsdPBARP);
      grSigSDPBARP->SetLineColor(1);
      grSigSDPBARP->SetLineWidth(1);
      grSigSDPBARP->SetLineStyle(2);
      grSigSDPBARP->Draw("c");
      
      TGraph *grSigSDPP = new TGraph(npfit,Wcm,sigtotsdPP);
      grSigSDPP->SetLineColor(1);
      grSigSDPP->SetLineWidth(1);
      grSigSDPP->SetLineStyle(1);
      grSigSDPP->Draw("c");
      
      TGraph *grSigDDPBARP = new TGraph(npfit,Wcm,sigtotddPBARP);
      grSigDDPBARP->SetLineColor(1);
      grSigDDPBARP->SetLineWidth(1);
      grSigDDPBARP->SetLineStyle(2);
      grSigDDPBARP->Draw("c");
      
      TGraph *grSigDDPP = new TGraph(npfit,Wcm,sigtotddPP);
      grSigDDPP->SetLineColor(1);
      grSigDDPP->SetLineWidth(1);
      grSigDDPP->SetLineStyle(1);
      grSigDDPP->Draw("c");  
           
      gr0->GetYaxis()->SetTitle("#bf{#sigma_{SD},#sigma_{DD} [mb]}");
      gr0->GetYaxis()->SetRangeUser(y_min,y_max);
      gr0->GetYaxis()->SetTitleOffset(1.2);
      gr0->GetYaxis()->SetTitleFont(42);
      gr0->GetXaxis()->SetTitle("#bf{#sqrt{s} [GeV]}");
      gr0->GetXaxis()->SetTitleOffset(1.3);
      gr0->GetXaxis()->SetLimits(x_min,x_max);   
      gr0->GetXaxis()->SetTitleFont(42);
           
      TLegend *leg1 = new TLegend(0.65,0.15,0.89,0.42);
      leg1->AddEntry(gr1,"TOTEM for 3.4 GeV < M_{SD} < 1100 GeV","p");
      leg1->AddEntry(gr2,"CMS for M_{X}^{2}/s < 0.05 or M_{Y}^{2}/s < 0.05","p");
      leg1->AddEntry(gr3,"ALICE for M_{X} < 200 GeV","p");
      leg1->AddEntry(gr0,"ISR for M_{X}^{2}/s < 0.05","p");
      leg1->AddEntry(gr1,"UA5 for M_{X}^{2}/s < 0.05","p");
      leg1->AddEntry(gr6,"CDF for M_{X}^{2}/s < 0.2","p");
      leg1->AddEntry(gr7,"E710 for M_{X}^{2}/s < 0.05","p");
      leg1->AddEntry(gr8,"UA4 for M_{X}^{2}/s < 0.05","p");                      
      leg1->AddEntry(grSigSDPP,"SD pp","l");
      leg1->AddEntry(grSigSDPBARP,"SD #bar{p}p","l");
      leg1->AddEntry(grSigDDPP,"SD pp","l");
      leg1->AddEntry(grSigDDPBARP,"SD #bar{p}p","l");
      leg1->SetFillColor(0);
      leg1->SetShadowColor(0);
      leg1->SetBorderSize(0);
      leg1->SetTextSize(0.035);
      leg1->Draw();    
    
      canv1->SaveAs("SigDiffs_BSG20mod_V002zc.eps");    
}
//----------------------------------------------------------------------
void fcn(int &npar, double *gin, double &f, double *par, int iflag)
{    
   /*******************************************************/ 
    //fit parameters
    double ww=1.9;//par[0]; 
    double Sig0=135.9;//par[1];
    double A1=32.3;//par[2];
    double del1=0.5;//par[3];
    double A2=39.4;//par[3];  
    double del2=0.5;//par[5];
    double B=2.0;//par[3];
    /*******************************************************/
    //reading data, then fitting

    //
    // St pp e pbp
    // 
    FILE *aq1,*aq2;
 
    // pbp   
    aq1 = fopen("SigSDpbp.dat","r");
    aq2 = fopen("SigSDpp.dat","r");

    const int npdpa = npSigSDpbp;
    double Wpap[npdpa],SigExpPaP[npdpa],uWpap[npdpa],uSigExpPaP[npdpa];
   
    for(int i=0;i<npdpa;i++)     
        fscanf(aq1,"%lg %lg %lg %lg",&Wpap[i],&SigExpPaP[i],&uSigExpPaP[i],&uWpap[i]);
        
     // pp   
     const int npdpp = npSigSDpp;
     double Wpp[npdpp],SigExpPP[npdpp],uWpp[npdpp],uSigExpPP[npdpp];
     
     for(int i=0;i<npdpp;i++)   
        fscanf(aq2,"%lg %lg %lg %lg",&Wpp[i],&SigExpPP[i],&uSigExpPP[i],&uWpp[i]);
       
     fclose(aq1);
     fclose(aq2);   

    /*******************************************************/
    //Total SD XSection - pp
       
    TF1 SigTotSDPP("SigTotSDPP",KerForward,bmin,bmax,10);
    SigTotSDPP.SetParameter(0,ww);
    SigTotSDPP.SetParameter(1,Sig0);
    SigTotSDPP.SetParameter(2,A1);
    SigTotSDPP.SetParameter(3,del1);
    SigTotSDPP.SetParameter(4,A2);
    SigTotSDPP.SetParameter(5,del2);
    SigTotSDPP.SetParameter(6,B);
    SigTotSDPP.SetParameter(7,1);
    SigTotSDPP.SetParameter(8,1);    
    
    /*******************************************************/
    //Total SD XSection - pbp
       
    TF1 SigTotSDPBARP("SigTotSDPBARP",KerForward,bmin,bmax,10);
    SigTotSDPBARP.SetParameter(0,ww);
    SigTotSDPBARP.SetParameter(1,Sig0);
    SigTotSDPBARP.SetParameter(2,A1);
    SigTotSDPBARP.SetParameter(3,del1);
    SigTotSDPBARP.SetParameter(4,A2);
    SigTotSDPBARP.SetParameter(5,del2);
    SigTotSDPBARP.SetParameter(6,B);
    SigTotSDPBARP.SetParameter(7,2);
    SigTotSDPBARP.SetParameter(8,1);

    /*******************************************************/
    //Total DD XSection - pp
       
    TF1 SigTotDDPP("SigTotDDPP",KerForward,bmin,bmax,10);
    SigTotDDPP.SetParameter(0,ww);
    SigTotDDPP.SetParameter(1,Sig0);
    SigTotDDPP.SetParameter(2,A1);
    SigTotDDPP.SetParameter(3,del1);
    SigTotDDPP.SetParameter(4,A2);
    SigTotDDPP.SetParameter(5,del2);
    SigTotDDPP.SetParameter(6,B);
    SigTotDDPP.SetParameter(7,1);
    SigTotDDPP.SetParameter(8,2);

    /*******************************************************/
    //Total DD XSection - pbp
       
    TF1 SigTotDDPBARP("SigTotDDPBARP",KerForward,bmin,bmax,10);
    SigTotDDPBARP.SetParameter(0,ww);
    SigTotDDPBARP.SetParameter(1,Sig0);
    SigTotDDPBARP.SetParameter(2,A1);
    SigTotDDPBARP.SetParameter(3,del1);
    SigTotDDPBARP.SetParameter(4,A2);
    SigTotDDPBARP.SetParameter(5,del2);
    SigTotDDPBARP.SetParameter(6,B);
    SigTotDDPBARP.SetParameter(7,2);
    SigTotDDPBARP.SetParameter(8,2);
        
    /*******************************************************/
    
    //Chi-squared calculation
    int np = 50;
    double *x=new double[np];
    double *w=new double[np];
    double sigsdpp,sigsdpbarp,sigddpp,sigddpbarp;    
    
    double delta;
    double chisq=0.;
    
    for(int i=0;i<npdpp;i++){
        SigTotSDPP.SetParameter(9,Wpp[i]);
        SigTotSDPP.CalcGaussLegendreSamplingPoints(np,x,w,1e-15);
        sigsdpp=SigTotSDPP.IntegralFast(np,x,w,bmin,bmax);
        delta=(SigExpPP[i]-sigsdpp)/uSigExpPP[i];
        chisq+=Power(delta,2);          
   }
       
    delta=0.;

    for(int i=0;i<npdpa;i++){
        SigTotSDPBARP.SetParameter(9,Wpap[i]);
        SigTotSDPBARP.CalcGaussLegendreSamplingPoints(np,x,w,1e-15);
        sigsdpbarp=SigTotSDPBARP.IntegralFast(np,x,w,bmin,bmax);
        delta=(SigExpPaP[i]-sigsdpbarp)/uSigExpPaP[i];
        chisq+=Power(delta,2);           
   }
   
    delta=0.; 
/*
    for(int i=0;i<npdpp;i++){
        SigTotDDPP.SetParameter(9,Wpp[i]);
        SigTotDDPP.CalcGaussLegendreSamplingPoints(np,x,w,1e-15);
        sigddpp=SigTotDDPP.IntegralFast(np,x,w,bmin,bmax);
        delta=(SigExpPP[i]-sigddpp)/uSigExpPP[i];
        chisq+=Power(delta,2);          
   }
       
    delta=0.;

    for(int i=0;i<npdpa;i++){
        SigTotDDPBARP.SetParameter(9,Wpap[i]);
        SigTotDDPBARP.CalcGaussLegendreSamplingPoints(np,x,w,1e-15);
        sigddpbarp=SigTotDDPBARP.IntegralFast(np,x,w,bmin,bmax);
        delta=(SigExpPaP[i]-sigddpbarp)/uSigExpPaP[i];
        chisq+=Power(delta,2);           
   }
   
    delta=0.;
*/       
   f = chisq;        
}
//----------------------------------------------------------------------

void BSG20mod_Diffs_V002zc()
{
    /******************************************/
    //time control - start
    auto start = high_resolution_clock::now();    
    /******************************************/ 
         
    TMinuit *gMinuit = new TMinuit(numpar);  //initialize TMinuit with a maximum of 4 params
    gMinuit->SetFCN(fcn);
    
    double arglist[numpar];
    int ierflg = 0;
  
    arglist[0] = 4.72; // --> número qualaquer //fitting with CL = 1\sigma
    gMinuit->mnexcm("SET ERR",arglist,1,ierflg); 
/*      
    double ww=1.4;//par[0]; 
    double Sig0=135.7;//par[1];
    double A1=32.3;//par[2];
    double del1=0.5;//par[3];
    double A2=39.3;//par[3];  
    double del2=0.5;//par[5];
    double B=2.0;//par[3];     
*/      
    double vstart[numpar] = {
    1.9,
    135.9,
    32.3,
    39.4,
    }; //start values
    double step[numpar] = {1.e-3,1.e-3,1.e-3,1.e-3}; //steps
//    gMinuit->mnparm(0, "w",    vstart[0],   step[0], 0,0,ierflg);
//    gMinuit->mnparm(1, "Sig0", vstart[1],   step[1], 0,0,ierflg);
//    gMinuit->mnparm(2, "A1",   vstart[2],   step[2], 0,0,ierflg);
//    gMinuit->mnparm(3, "del1", vstart[3],   step[3], 0,0,ierflg);
//    gMinuit->mnparm(3, "A2",   vstart[3],   step[3], 0,0,ierflg);
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
    double ww=1.9;//outpar[0]; 
    double Sig0=135.9;//outpar[1];
    double A1=32.3;//outpar[2];
    double del1=0.5;//outpar[3];
    double A2=39.4;//outpar[3];  
    double del2=0.5;//outpar[5];
    double B=2.0;//outpar[3];     
/*    
\begin{table}[ht!]
 \centering
 \caption{\label{tab:res_pp_ppbar_Bfixed}Results of fits to $\sigmatot$ data from $pp$ and $\bar{p}p$ scattering with $B$ fixed and also reggeon intercepts $\delta_1 = \delta_2 = 0.5$ fixed. We considered all TOTEM and ATLAS data.}
 \begin{tabular}{c|c|c|c}\hline\hline
 Parameters     & \multicolumn{3}{c}{$B$ fixed}\\\hline
 $A_1$ (mb)     & 24.9(1.4)  & 32.3(1.8)  & 39.7(2.1)  \\
 $\delta_1$     & 0.5 (fixed)& 0.5 (fixed)& 0.5 (fixed)\\
 $A_2$ (mb)     & 31.2(1.1)  & 39.4(1.4)  & 47.2(1.7)  \\
 $\delta_2$     & 0.5 (fixed)& 0.5 (fixed)& 0.5 (fixed)\\
 $B$            & 1.5 (fixed)& 2.0 (fixed)& 2.5 (fixed)\\
 $w$            & 1.309(40)  & 1.660(43)  & 1.940(44)  \\
 $\sigma_0$ (mb)& 114.3(2.6) & 135.9(3.3) & 155.7(3.9) \\\hline
 $\chi^2/$dof   & 7.16       & 6.78       & 6.49       \\
 dof            & 173        & 173        & 173        \\\hline\hline
 \end{tabular}
\end{table}
*/ 
    /*******************************************************/
    //Total SD XSection - pp
       
    TF1 SigTotSDPP("SigTotSDPP",KerForward,bmin,bmax,10);
    SigTotSDPP.SetParameter(0,ww);
    SigTotSDPP.SetParameter(1,Sig0);
    SigTotSDPP.SetParameter(2,A1);
    SigTotSDPP.SetParameter(3,del1);
    SigTotSDPP.SetParameter(4,A2);
    SigTotSDPP.SetParameter(5,del2);
    SigTotSDPP.SetParameter(6,B);
    SigTotSDPP.SetParameter(7,1);
    SigTotSDPP.SetParameter(8,1);    
    
    /*******************************************************/
    //Total SD XSection - pbp
       
    TF1 SigTotSDPBARP("SigTotSDPBARP",KerForward,bmin,bmax,10);
    SigTotSDPBARP.SetParameter(0,ww);
    SigTotSDPBARP.SetParameter(1,Sig0);
    SigTotSDPBARP.SetParameter(2,A1);
    SigTotSDPBARP.SetParameter(3,del1);
    SigTotSDPBARP.SetParameter(4,A2);
    SigTotSDPBARP.SetParameter(5,del2);
    SigTotSDPBARP.SetParameter(6,B);
    SigTotSDPBARP.SetParameter(7,2);
    SigTotSDPBARP.SetParameter(8,1);

    /*******************************************************/
    //Total DD XSection - pp
       
    TF1 SigTotDDPP("SigTotDDPP",KerForward,bmin,bmax,10);
    SigTotDDPP.SetParameter(0,ww);
    SigTotDDPP.SetParameter(1,Sig0);
    SigTotDDPP.SetParameter(2,A1);
    SigTotDDPP.SetParameter(3,del1);
    SigTotDDPP.SetParameter(4,A2);
    SigTotDDPP.SetParameter(5,del2);
    SigTotDDPP.SetParameter(6,B);
    SigTotDDPP.SetParameter(7,1);
    SigTotDDPP.SetParameter(8,2);

    /*******************************************************/
    //Total DD XSection - pbp
       
    TF1 SigTotDDPBARP("SigTotDDPBARP",KerForward,bmin,bmax,10);
    SigTotDDPBARP.SetParameter(0,ww);
    SigTotDDPBARP.SetParameter(1,Sig0);
    SigTotDDPBARP.SetParameter(2,A1);
    SigTotDDPBARP.SetParameter(3,del1);
    SigTotDDPBARP.SetParameter(4,A2);
    SigTotDDPBARP.SetParameter(5,del2);
    SigTotDDPBARP.SetParameter(6,B);
    SigTotDDPBARP.SetParameter(7,2);
    SigTotDDPBARP.SetParameter(8,2);
        
    /*******************************************************/  
//----------------------------------------------------------------------
   //calculates total xsection and rho, creating files to make the plots
    FILE *aq1,*aq2,*aq3,*aq4;
   
    aq1=fopen("SigSDpp_BSG20mod_Diffs_V002zc.dat","w"); 
    aq2=fopen("SigSDpbp_BSG20mod_Diffs_V002zc.dat","w"); 
    aq3=fopen("SigDDpp_BSG20mod_Diffs_V002zc.dat","w"); 
    aq4=fopen("SigDDpbp_BSG20mod_Diffs_V002zc.dat","w"); 
           
     //number of directive computing for GaussLegendreIntegration
    int np = 50;
    double *x=new double[np];
    double *w=new double[np];
    
    double i=0.;
    const double di=0.05;
    double Ecm;
    double sigsdpp,sigsdpbarp,sigddpp,sigddpbarp;
    
     do{
            Ecm=5.*pow(10.,i);
            /*************************************************/
            //pp
            SigTotSDPP.SetParameter(9,Ecm);
            SigTotSDPP.CalcGaussLegendreSamplingPoints(np,x,w,1e-15);
            sigsdpp = SigTotSDPP.IntegralFast(np,x,w,bmin,bmax);            
            fprintf(aq1,"%.2e   %.4lf\n",Ecm,sigsdpp); 
            SigTotDDPP.SetParameter(9,Ecm);
            SigTotDDPP.CalcGaussLegendreSamplingPoints(np,x,w,1e-15);
            sigddpp = SigTotDDPP.IntegralFast(np,x,w,bmin,bmax);            
            fprintf(aq3,"%.2e   %.4lf\n",Ecm,sigddpp);            
            /*************************************************/
            //pbarp
            SigTotSDPBARP.SetParameter(9,Ecm);
            SigTotSDPBARP.CalcGaussLegendreSamplingPoints(np,x,w,1e-15);
            sigsdpbarp = SigTotSDPBARP.IntegralFast(np,x,w,bmin,bmax);
            fprintf(aq2,"%.2e   %.4lf\n",Ecm,sigsdpbarp);
            SigTotDDPBARP.SetParameter(9,Ecm);
            SigTotDDPBARP.CalcGaussLegendreSamplingPoints(np,x,w,1e-15);
            sigddpbarp = SigTotDDPBARP.IntegralFast(np,x,w,bmin,bmax);
            fprintf(aq4,"%.2e   %.4lf\n",Ecm,sigddpbarp);            
            /*************************************************/                             
            i+=di;
            
                   //cout << Ecm << "\t" <<sigpp << "\t" <<sigpbarp<< endl;
            
       }while(Ecm<Wmax); 
       
      fclose(aq1);
      fclose(aq2); 
      fclose(aq3);
      fclose(aq4);
       
    aq1=fopen("SigSDpp_BSG20mod_Diffs_V002zc.dat","r"); 
    aq2=fopen("SigSDpbp_BSG20mod_Diffs_V002zc.dat","r"); 
    aq3=fopen("SigDDpp_BSG20mod_Diffs_V002zc.dat","r"); 
    aq4=fopen("SigDDpbp_BSG20mod_Diffs_V002zc.dat","r");
      
      const int npfit=108;
      double Wcm[npfit],sigtotsdPP[npfit],sigtotsdPBARP[npfit],sigtotddPP[npfit],sigtotddPBARP[npfit];       
            
      for(int j=0;j<npfit;j++){
          fscanf(aq1,"%lg %lg",&Wcm[j],&sigtotsdPP[j]);
          fscanf(aq2,"%lg %lg",&Wcm[j],&sigtotsdPBARP[j]);
          fscanf(aq3,"%lg %lg",&Wcm[j],&sigtotddPP[j]);
          fscanf(aq4,"%lg %lg",&Wcm[j],&sigtotddPBARP[j]);          
      }
      
      fclose(aq1);
      fclose(aq2); 
      fclose(aq3);
      fclose(aq4);    
       
//----------------------------------------------------------------------
    //Plotting total elastic xsection data         
    PlotTotXSec(Wcm,sigtotsdPP,sigtotsdPBARP,sigtotddPP,sigtotddPBARP,npfit);      
   
    /******************************************/
    auto stop = high_resolution_clock::now(); 
    auto duration = duration_cast<seconds>(stop - start); 
    cout << "Time(in min): "<<duration.count()/6.e1 << endl;
    /******************************************/ 
        
}
