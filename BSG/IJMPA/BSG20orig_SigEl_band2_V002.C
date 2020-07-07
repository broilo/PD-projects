/* PROGRAM: BSG20 based on Lipari-Lusignoli original model

   ---> Ref. Original model: PHYSICAL REVIEW D 80, 074014 (2009)
        SciHub link: https://sci-hub.tw/10.1103/PhysRevD.80.074014
        
   Modified version: SigEik=SigSoft+SigPQCD

------------------------------------------------------------------------
   *** Main characteristics ***

	> SigEik = SigSoft + SigPQCD
	> SigSoft = A1*s^{-del1} +/- A2*s^{-del2} + Sig0
	> SigPQCD = B*SigPQCD
	
	> Calculates SigEl

   *** Observations ***

 Obs1: TOTEM + ATLAS dataset
 Obs2: Real parametrization for SigQCD with 10 free-parameters
 Obs3: PDF: CT14
 Obs4: Free: w, Sig0, A1, A2 -> Fixados no ajuste de BSGorig_V002.C
 Obs5: B is fixed in 0.76

 *** Bandas de incerteza ***

  Estratégia 1 -> Manter os parâmetros nos valores ajustados e estudar a incerteza em w.
 
 Intervalo w=[1.0,(1.5)]

   *** Important *** 
   (See: https://arxiv.org/abs/1908.01040 Chapter 8, sections 8.4 and 8.8)
   
   mg=0.4d0 !GeV -> mg=500+/-200 MeV for lmb=300 MeV
 
 --> (MSTW) lbd=0.318920 GeV lambda equivalent to alpha_LO(Mz^2)=0.13939~0.139 macthing-prescrition squeme; m_b=4.18 GeV
 
 --> (MMHT) lbd=0.267471 GeV lambda equivalent to alpha_LO(Mz^2)=0.13499~0.135 macthing-prescrition squeme; m_b=4.18 GeV
 
 --> (CTEQ6L) lbd=0.3260 GeV lambda equivalent to alpha_NLO(Mz^2)=0.118 for 4 flavors
 
 --> (CT14)   lbd=0.3260 GeV lambda equivalent to alpha_NLO(Mz^2)=0.118 for 4 flavors

   *** How to run the code ***
 
 root -l -q /usr/lib/x86_64-linux-gnu/libgsl.so /usr/lib/x86_64-linux-gnu/libgslcblas.so BSG20orig_SigEl_band2_V002.C | cat > BSG20orig_SigEl_band2_V002.min
 
   ---------------------------------------------------------------------
   Codes by: M.Broilo*, P.V.R.G.Silva, V.P.B.Gonçalves
   mateus.broilo@ufrgs.br | mateus.broilo90@gmail.com
   High and Medium Energy Group
   Grupo de Altas e Medias Energias (GAME)
   Universidade Federal de Pelotas, Pelotas - RS (Brasil)
   ---------------------------------------------------------------------
   
   Creation: 29/jan/2020 (Pelotas, RS, Brazil)
   Last update: 29/jan/2020 (Pelotas, RS, Brazil)
   
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
const double s0=1.;
const double r0=0.19;//fm
//----------------------------------------------------------------------
//number of fit parameters and data points 
const int numpar=4;
const int npSap=141;
const int npSpp=155;
const int npapfit=npSap;
const int nppfit=npSpp;
const int npMin=npapfit+nppfit;
// Plots
const int npSppPlot=147;
const int npSppTOTEM=6;
const int npSppATLAS=2;
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
    double W=par[8];        
        
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
    double SigR1=A1*Power(s,-del1);
    double SigR2=A2*Power(s,-del2);
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
	xx=(nn*ww*ww)/2.;
	/*******************************************************/
            KerSig=2.*Pi()*b*Power((1.-Power((1.+xx),(-winv))),2.)*fm2_to_mb; //in mb 
            return KerSig;    
}
//----------------------------------------------------------------------
// Plotando a seção de choque total pp e pbp
//----------------------------------------------------------------------
void PlotElXSec(double *Wcm,double *sigelPPEIK,double *sigelPBARPEIK,int npfit)
{
    FILE *aq1,*aq2,*aq3,*aq4;
// Plotting St data     
    aq1 = fopen("SigElPBP.dat","r");
    aq2 = fopen("SigElPP.dat","r");
    aq3 = fopen("SigElPP_TOTEM.dat","r");
    aq4 = fopen("SigElPP_ATLAS.dat","r");

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
       
     for(int i=0;i<npTOTEM;i++)   
        fscanf(aq4,"%lg %lg %lg %lg",&Wpp3[i],&SigExpPP3[i],&uSigExpPP3[i],&uWpp3[i]);
      
     fclose(aq1);
     fclose(aq2);
     fclose(aq3);
     fclose(aq4);
     
      TCanvas *canv1 = new TCanvas("c1","Total Cross Section",200,10,900,600);
      canv1->SetFillColor(0);
      canv1->SetLogx();
       
      TGraphErrors *gr0 = new TGraphErrors(npdpa,Wpap,SigExpPaP,uWpap,uSigExpPaP);
      gr0->SetMarkerStyle(24);
      gr0->SetMarkerColor(1);
      gr0->SetLineColor(1);
      gr0->SetMarkerSize(1.0);
      gr0->SetTitle();
      gr0->Draw("apz");

       double x_min(5e00),x_max(2.e05);
       double y_min(0.0),y_max(60.);
       
      gr0->GetYaxis()->SetTitle("#bf{#sigma_{el} [mb]}");
      gr0->GetYaxis()->SetRangeUser(y_min,y_max);
      gr0->GetYaxis()->SetTitleOffset(1.2);
      gr0->GetYaxis()->SetTitleFont(42);
      gr0->GetXaxis()->SetTitle("#bf{#sqrt{s} [GeV]}");
      gr0->GetXaxis()->SetTitleOffset(1.3);
      gr0->GetXaxis()->SetLimits(x_min,x_max);   
      gr0->GetXaxis()->SetTitleFont(42);     
           
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
      
      TGraph *grSigPBARP_Eik = new TGraph(npfit,Wcm,sigelPBARPEIK);
      grSigPBARP_Eik->SetLineColor(4);
      grSigPBARP_Eik->SetLineWidth(1);
      grSigPBARP_Eik->SetLineStyle(1);
      grSigPBARP_Eik->Draw("c");
      
      TGraph *grSigPP_Eik = new TGraph(npfit,Wcm,sigelPPEIK);
      grSigPP_Eik->SetLineColor(2);
      grSigPP_Eik->SetLineWidth(1);
      grSigPP_Eik->SetLineStyle(2);
      grSigPP_Eik->Draw("c"); 
           
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
      
      canv1->SaveAs("SigEl_BSG20orig_SigEl_band2_V002.eps");    
}
//----------------------------------------------------------------------
void fcn(int &npar, double *gin, double &f, double *par, int iflag)
{    
   /*******************************************************/ 
    //fit parameters
    double ww=1.0;//par[0]; 
    double Sig0=7.35489e+01;//par[0];
    double A1=6.51125e+01;//par[1];
    double del1=0.5;//par[3];
    double A2=8.60594e+01;//par[2];  
    double del2=0.5;//par[5];
    double B=7.55541e-01;//par[3];
    /*******************************************************/
    //reading data, then fitting

    //
    // St pp e pbp
    // 
    FILE *aq1,*aq2;
 
    // pbp   
    aq1 = fopen("SigElPBP.dat","r");
    aq2 = fopen("SigElPP.dat","r");

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

    /*******************************************************/
    //Total XSection - pp
       
    TF1 SigElPP("SigElPP",KerForward,bmin,bmax,9);
    SigElPP.SetParameter(0,ww);
    SigElPP.SetParameter(1,Sig0);
    SigElPP.SetParameter(2,A1);
    SigElPP.SetParameter(3,del1);
    SigElPP.SetParameter(4,A2);
    SigElPP.SetParameter(5,del2);
    SigElPP.SetParameter(6,B);
    SigElPP.SetParameter(7,1);    
    
    /*******************************************************/
    //Total XSection - pbp
       
    TF1 SigElPBARP("SigElPBARP",KerForward,bmin,bmax,9);
    SigElPBARP.SetParameter(0,ww);
    SigElPBARP.SetParameter(1,Sig0);
    SigElPBARP.SetParameter(2,A1);
    SigElPBARP.SetParameter(3,del1);
    SigElPBARP.SetParameter(4,A2);
    SigElPBARP.SetParameter(5,del2);
    SigElPBARP.SetParameter(6,B);
    SigElPBARP.SetParameter(7,2);
        
    /*******************************************************/
    
    //Chi-squared calculation
    int np = 50;
    double *x=new double[np];
    double *w=new double[np];
    double sigpp,sigpbarp;    
    
    double delta;
    double chisq=0.;
    
    for(int i=0;i<npdpp;i++){
        SigElPP.SetParameter(8,Wpp[i]);
        SigElPP.CalcGaussLegendreSamplingPoints(np,x,w,1e-15);
        sigpp=SigElPP.IntegralFast(np,x,w,bmin,bmax);
        delta=(SigExpPP[i]-sigpp)/uSigExpPP[i];
        chisq+=Power(delta,2);          
   }
       
    delta=0.;

    for(int i=0;i<npdpa;i++){
        SigElPBARP.SetParameter(8,Wpap[i]);
        SigElPBARP.CalcGaussLegendreSamplingPoints(np,x,w,1e-15);
        sigpbarp=SigElPBARP.IntegralFast(np,x,w,bmin,bmax);
        delta=(SigExpPaP[i]-sigpbarp)/uSigExpPaP[i];
        chisq+=Power(delta,2);           
   }
   
    delta=0.; 
       
   f = chisq;        
}
//----------------------------------------------------------------------

void BSG20orig_SigEl_band2_V002()
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
      
    double vstart[numpar] = {77.0,35.0,70.0,0.68};  
/*    
    double ww=1.0;//outpar[0]; 
    double Sig0=77.;//outpar[1];
    double A1=35.;//outpar[2];
    double del1=0.5;//outpar[3];
    double A2=70.;//outpar[4];  
    double del2=0.5;//outpar[5];
    double B=0.68;//outpar[6]
{
    1.329,
    119.6,
    42.3,
    1.365,
    36.55,
    0.585,
    1.5}; //start values */
    double step[numpar] = {1.e-3,1.e-3,1.e-3,1.e-3}; //steps
//    gMinuit->mnparm(0, "w",    vstart[0],   step[0], 0,0,ierflg);
//    gMinuit->mnparm(0, "Sig0", vstart[0],   step[0], 0,0,ierflg);
//    gMinuit->mnparm(1, "A1",   vstart[1],   step[1], 0,0,ierflg);
//    gMinuit->mnparm(3, "del1", vstart[3],   step[3], 0,0,ierflg);
//    gMinuit->mnparm(2, "A2",   vstart[2],   step[2], 0,0,ierflg);
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
    double ww=1.5;//outpar[0]; 
    double Sig0=5.98107e+01;//outpar[0];
    double A1=5.69369e+01;//outpar[1];
    double del1=0.5;//outpar[3];
    double A2=7.28361e+01;//outpar[2];  
    double del2=0.5;//outpar[5];
    double B=7.55541e-01;//outpar[3];
        /*******************************************************/
    //Total XSection - pp
       
    TF1 SigElPP("SigElPP",KerForward,bmin,bmax,9);
    SigElPP.SetParameter(0,ww);
    SigElPP.SetParameter(1,Sig0);
    SigElPP.SetParameter(2,A1);
    SigElPP.SetParameter(3,del1);
    SigElPP.SetParameter(4,A2);
    SigElPP.SetParameter(5,del2);
    SigElPP.SetParameter(6,B);
    SigElPP.SetParameter(7,1);    
    
    /*******************************************************/
    //Total XSection - pbp
       
    TF1 SigElPBARP("SigElPBARP",KerForward,bmin,bmax,9);
    SigElPBARP.SetParameter(0,ww);
    SigElPBARP.SetParameter(1,Sig0);
    SigElPBARP.SetParameter(2,A1);
    SigElPBARP.SetParameter(3,del1);
    SigElPBARP.SetParameter(4,A2);
    SigElPBARP.SetParameter(5,del2);
    SigElPBARP.SetParameter(6,B);
    SigElPBARP.SetParameter(7,2);
        
    /*******************************************************/  
//----------------------------------------------------------------------
   //calculates total xsection and rho, creating files to make the plots
    FILE *aq1,*aq2;
   
    aq1=fopen("SigElpp_BSG20orig_SigEl_band2_V002.dat","w"); 
    aq2=fopen("SigElpbp_BSG20orig_SigEl_band2_V002.dat","w"); 
           
     //number of directive computing for GaussLegendreIntegration
    int np = 50;
    double *x=new double[np];
    double *w=new double[np];
    
    double i=0.;
    const double di=0.05;
    double Ecm;
    double sigpp,sigpbarp;
    
     do{
            Ecm=5.*pow(10.,i);
            /*************************************************/
            //pp
            SigElPP.SetParameter(8,Ecm);
            SigElPP.CalcGaussLegendreSamplingPoints(np,x,w,1e-15);
            sigpp = SigElPP.IntegralFast(np,x,w,bmin,bmax);            
            fprintf(aq1,"%.2e   %.4lf\n",Ecm,sigpp); 
            /*************************************************/
            //pbarp
            SigElPBARP.SetParameter(8,Ecm);
            SigElPBARP.CalcGaussLegendreSamplingPoints(np,x,w,1e-15);
            sigpbarp = SigElPBARP.IntegralFast(np,x,w,bmin,bmax);
            fprintf(aq2,"%.2e   %.4lf\n",Ecm,sigpbarp);           
            /*************************************************/                             
            i+=di;
            
                   //cout << Ecm << "\t" <<sigpp << "\t" <<sigpbarp<< endl;
            
       }while(Ecm<Wmax); 
       
      fclose(aq1);
      fclose(aq2); 
       
      aq1 = fopen("SigElpp_BSG20orig_SigEl_band2_V002.dat","r");
      aq2 = fopen("SigElpbp_BSG20orig_SigEl_band2_V002.dat","r");
      
      const int npfit=108;
      double Wcm[npfit],sigelPPEIK[npfit],sigelPBARPEIK[npfit];       
            
      for(int j=0;j<npfit;j++){
          fscanf(aq1,"%lg %lg",&Wcm[j],&sigelPPEIK[j]);
          fscanf(aq2,"%lg %lg",&Wcm[j],&sigelPBARPEIK[j]);
      }
      
      fclose(aq1);
      fclose(aq2);    
       
//----------------------------------------------------------------------
    //Plotting total xsection data         
    PlotElXSec(Wcm,sigelPPEIK,sigelPBARPEIK,npfit);      
   
    /******************************************/
    auto stop = high_resolution_clock::now(); 
    auto duration = duration_cast<seconds>(stop - start); 
    cout << "Time(in min): "<<duration.count()/6.e1 << endl;
    /******************************************/ 
        
}
