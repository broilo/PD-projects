/* DGM19cCT14cut7_FFpresc -> DGM19cCT14cut7_FFpresc_mexido
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c O FATOR DE FORMA ESTÁ SENDO CONSIDERADO NA RDD!!!
c
c ImXsh é obtida via s->-is
c
c Foi retirado do conjunto de dados os pontos de StRh(8e13TeV)
c
c OBS: Eu mudei a incerteza de Rh(7TeV) para ver se Rh(S) continua crescendo
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
c Obs7: SigQCD(s) foi parametrizada por uma função complexa
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
*/
#include <chrono> 
#include "Math/WrappedTF1.h"
using namespace chrono;
using namespace TMath;
//
const double mbfactor=0.389379; //mbGeV^{2}
//
/********************************************************************************/
//constant parameters throughout the model
const double s0=25.; 
const double muoddsoft=0.5;
const double Kf=1.; 
/********************************************************************************/
//number of fit parameters and data points 
const int numpar = 7;
const int npRap=12;
const int npRpp=49;//52;
const int npSap=30;
const int npSpp=72;//80;
const int npapfit=npRap+npSap;
const int nppfit=npRpp+npSpp;
const int npMin=npapfit+nppfit;
// Plots
const int npSppPlot=80;
const int npRppPlot=52;
/********************************************************************************/
//parameters to control energy and b range (for plots and integrations)
const double Wmin=5.;
const double Wmax=1.e6;
const double bmin=0.;
const double bmax=30.;   
/*******************************************************************************/
//QCD cross section - Mateus' analytical parametrization for CT14LO PDF
/*******************************************************************************/
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
const double c1=1.01;
const double c2=1.05;
const double c3=1.09;  
const double bkg=100.;
/*******************************************************************************/
// Modified Bessel function of first kind and third order
/*******************************************************************************/
double BesselK3(double z)
{
   if(z>0.)
       return BesselK(3,z);
   else
       return 0.;     
}
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc/
/*
                             Forward Analysis 
*/
//cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc/
// Aqui começa a ser escrito parte do código que lida com os observáveis frontais
/*******************************************************************************/
// Sigma QCD - Complex
/*******************************************************************************/
//real part
double ReSigQCD(double W)
{
   TComplex s=-TComplex::I()*W*W;
        
    TComplex Y=TComplex::Log(s);
    TComplex X=TComplex::Log(Y);
    
    TComplex sig=(b1-bkg)+b2*TComplex::Exp(b3*TComplex::Power(c1*X,b4))+b5*TComplex::Exp(b6*TComplex::Power(c2*X,b7))+b8*TComplex::Exp(b9*TComplex::Power(c3*X,b10));
      
    return sig.Re();
}
//imaginary part
double ImSigQCD(double W)
{    
    TComplex s=-TComplex::I()*W*W;
        
    TComplex Y=TComplex::Log(s);
    TComplex X=TComplex::Log(Y);

    TComplex sig=(b1-bkg)+b2*TComplex::Exp(b3*TComplex::Power(c1*X,b4))+b5*TComplex::Exp(b6*TComplex::Power(c2*X,b7))+b8*TComplex::Exp(b9*TComplex::Power(c3*X,b10));
          
    return sig.Im();
}
/*******************************************************************************/
//Form factors
/*******************************************************************************/
/*
Complex K3(x)
*/
//
// -> Real K3(x)
//
double ReK3(double *x,double *par)
{
     double b=x[0];
     double W=par[0];
     double nu1=par[1];
     double nu2=par[2];

     TComplex s=-TComplex::I()*W*W;
     TComplex Y=TComplex::Log(s/s0);
     TComplex mush=nu1-nu2*Y;
     
     double A=(mush.Re())*b;
     double B=(mush.Im())*b;
     
     const double tmin=0.;
     const double tmax=10.;
     
     TF1 WF_RE("WF_RE","TMath::Exp(-[0]*TMath::CosH(x))*TMath::Cos([1]*TMath::CosH(x))*TMath::CosH(3.*x)",tmin,tmax);
     WF_RE.SetParameter(0,A);
     WF_RE.SetParameter(1,B);   
     
     ROOT::Math::WrappedTF1 wf1(WF_RE);
   // Create the Integrator
     ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE);// GSL Adaptive 
     ig.SetFunction(wf1);
     ig.SetRelTolerance(1e-5);    
     
     return ig.Integral(tmin,tmax);
}
//
// -> Imaginary K3(x)
//
double ImK3(double *x,double *par)
{
     double b=x[0];
     double W=par[0];
     double nu1=par[1];
     double nu2=par[2];

     TComplex s=-TComplex::I()*W*W;
     TComplex Y=TComplex::Log(s/s0);
     TComplex mush=nu1-nu2*Y;
     
     double A=(mush.Re())*b;
     double B=(mush.Im())*b;
     
     const double tmin=0.;
     const double tmax=10.;
     
     TF1 WF_RE("WF_RE","-TMath::Exp(-[0]*TMath::CosH(x))*TMath::Sin([1]*TMath::CosH(x))*TMath::CosH(3.*x)",tmin,tmax);
     WF_RE.SetParameter(0,A);
     WF_RE.SetParameter(1,B);   
     
     ROOT::Math::WrappedTF1 wf1(WF_RE);
   // Create the Integrator
     ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE);// GSL Adaptive 
     ig.SetFunction(wf1);
     ig.SetRelTolerance(1e-5);    
     
     return ig.Integral(tmin,tmax); 
}
//semi-hard
double ReWSH(double b,double W,double nu1,double nu2)
{    
    TComplex s=-TComplex::I()*W*W;
    TComplex Y=TComplex::Log(s/s0);
    TComplex mush=nu1-nu2*Y;
    
    TF1 RealK3("RealK3",ReK3,bmin,bmax,3);
    RealK3.SetParameter(0,W);
    RealK3.SetParameter(1,nu1);
    RealK3.SetParameter(2,nu2);
    
    TF1 ImagK3("ImagK3",ImK3,bmin,bmax,3);
    ImagK3.SetParameter(0,W);
    ImagK3.SetParameter(1,nu1);
    ImagK3.SetParameter(2,nu2);    
    
    TComplex K3Complex = RealK3.Eval(b)+TComplex::I()*ImagK3.Eval(b);    
    
    TComplex ws=TComplex::Power(mush,5)*Power(b,3)*K3Complex/(96.*Pi());  
    
    return ws.Re();    
}
double ImWSH(double b,double W,double nu1,double nu2)
{    
    TComplex s=-TComplex::I()*W*W;
    TComplex Y=TComplex::Log(s/s0);
    TComplex mush=nu1-nu2*Y;
    
    TF1 RealK3("RealK3",ReK3,bmin,bmax,3);
    RealK3.SetParameter(0,W);
    RealK3.SetParameter(1,nu1);
    RealK3.SetParameter(2,nu2);
    
    TF1 ImagK3("ImagK3",ImK3,bmin,bmax,3);
    ImagK3.SetParameter(0,W);
    ImagK3.SetParameter(1,nu1);
    ImagK3.SetParameter(2,nu2);    
    
    TComplex K3Complex = RealK3.Eval(b)+TComplex::I()*ImagK3.Eval(b);    
    
    TComplex ws=TComplex::Power(mush,5)*Power(b,3)*K3Complex/(96.*Pi());  
    
    return ws.Im();   
}
//soft
double WSoft(double b,double musoft)
{
    double ws=Power(musoft,5)*Power(b,3)*BesselK3(musoft*b)/(96.*Pi());  
    
    return ws;    
}
/*******************************************************************************/
//Integration kernels
/*******************************************************************************/
//Forward observables - total xsection and rho
double KerForward(double *x,double *par)
{
    double b=x[0];      
    double muevsoft=par[0];
    double A=par[1];
    double B=par[2];
    double C=par[3];
    double D=par[4];
    double nu1=par[5];
    double nu2=par[6];        
    int reac=par[7]; 
    int obs=par[8];         
    double W=par[9]; 
        
    double s=W*W;
    double Y=Log(s/s0);
    double Y2=Y*Y;      
    /*******************************************************/
    //Form factors
    /*******************************************************/
    double ReWSh=ReWSH(b,W,nu1,nu2);
    double ImWSh=ImWSH(b,W,nu1,nu2);
    double WEvS=WSoft(b,muevsoft); 
    double WOddS=WSoft(b,muoddsoft); 
    /*******************************************************/
    //SigQCD
    /*******************************************************/
    double resig=ReSigQCD(W);
    double imsig=ImSigQCD(W);
    /*******************************************************/
    //Eikonal - real part
    /*******************************************************/
    double ReChiEvenSH=0.5*Kf*(ReWSh*resig-ImWSh*imsig);
    /*******************************************************/
    double ReChiEvenSoft=0.5*WEvS*(A+(B*Cos(PiOver4())/Sqrt(s/s0))+C*(Y2-Pi()*PiOver4()));     
    /*******************************************************/
    double ReChiOddSoft=0.5*WOddS*D*Cos(PiOver4())/Sqrt(s/s0);     
    /*******************************************************/
    double ReChiEven=ReChiEvenSH+ReChiEvenSoft;
    double ReChiOdd=ReChiOddSoft;         
    /*******************************************************/
    //Eikonal - imaginary part
    /*******************************************************/    
    double ImChiEvenSH=0.5*Kf*(ReWSh*imsig+ImWSh*resig);
    /*******************************************************/
    double ImChiEvenSoft=0.5*WEvS*((B*Sin(PiOver4())/Sqrt(s/s0))-C*Pi()*Y);
    /*******************************************************/
    double ImChiOddSoft=-0.5*WOddS*D*Sin(PiOver4())/Sqrt(s/s0);
    /*******************************************************/    
    double ImChiEven=ImChiEvenSH+ImChiEvenSoft;    
    double ImChiOdd=ImChiOddSoft;    
    /*******************************************************/  
    double ReChi;
    double ImChi;
    
    double KerSig;
    double KerRho;
    
        if(reac==1) //pp
        {
            ReChi=ReChiEven-ReChiOdd;
            ImChi=ImChiEven-ImChiOdd;
        }else if(reac==2) //pbarp
        {
            ReChi=ReChiEven+ReChiOdd;
            ImChi=ImChiEven+ImChiOdd;
        }else{ 
            ReChi=0.;
            ImChi=0.;
        }
      
        if(obs==1) //total xsection
        {      
            KerSig=4.*Pi()*mbfactor*b*(1.-Cos(ImChi)*Exp(-ReChi)); //in mb      
            return KerSig;
        
        }else if(obs==2)//rho
        {
           KerRho=-4.*Pi()*mbfactor*b*Sin(ImChi)*Exp(-ReChi); //in mb    
           return KerRho;
        }else{
           return 0.; 
        }      
    
}
/*******************************************************************************/
// Plotando a seção de choque total pp e pbp
/*******************************************************************************/
void PlotTotXSec(double *Wcm,double *sigtotPPCT14,double *sigtotPBARPCT14,int npfit)
{
    FILE *aq1,*aq2;
// Plotting St data     
    aq1 = fopen("paw_Stpac10GeV.dat","r");
    aq2 = fopen("paw_StppTAc10GeV.dat","r");

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
       
     fclose(aq1);
     fclose(aq2);
     
      TCanvas *canv1 = new TCanvas("c1","Total Cross Section",200,10,900,600);
      canv1->SetFillColor(0);
      canv1->SetLogx();
       
      TGraphErrors *gr0 = new TGraphErrors(npdpa,Wpap,SigExpPaP,uWpap,uSigExpPaP);
      gr0->SetMarkerStyle(24);
      gr0->SetMarkerColor(1);
      gr0->SetLineColor(1);
      gr0->SetMarkerSize(1.2);
      gr0->SetTitle();
      gr0->Draw("apz");

       double x_min(8e00),x_max(2.e05);
       double y_min(3.0e01),y_max(1.8e02);
       
      gr0->GetYaxis()->SetTitle("#bf{#sigma_{tot} [mb]}");
      gr0->GetYaxis()->SetRangeUser(y_min,y_max);
      gr0->GetYaxis()->SetTitleOffset(1.2);
      gr0->GetYaxis()->SetTitleFont(42);
      gr0->GetXaxis()->SetTitle("#bf{#sqrt{s} [GeV]}");
      gr0->GetXaxis()->SetTitleOffset(1.3);
      gr0->GetXaxis()->SetLimits(x_min,x_max);   
      gr0->GetXaxis()->SetTitleFont(42);     
           
      //
      // ppCR
      //
              
      const int npt_cosmics = 2;
      double x_cosmics[npt_cosmics] = {57e03,95e03};
      double y_cosmics[npt_cosmics] = {133.00,170.0};
      double erx_lowucosmics[npt_cosmics] = {0.0,0.0};
      double erx_uppucosmics[npt_cosmics] = {0.0,0.0};
      double ery_lowucosmics[npt_cosmics] = {28.7,51.0};
      double ery_uppucosmics[npt_cosmics] = {26.7,51.0};
           
      TGraphAsymmErrors *gr_cosmics = new TGraphAsymmErrors(npt_cosmics,x_cosmics,y_cosmics,erx_lowucosmics,erx_uppucosmics,ery_lowucosmics,ery_uppucosmics);
      gr_cosmics->SetMarkerStyle(22);
      gr_cosmics->SetMarkerSize(1.2);
      gr_cosmics->SetMarkerColorAlpha(kGray,0.8);
      gr_cosmics->SetLineColorAlpha(kGray,0.8);
      gr_cosmics->Draw("pz"); 
      
       TGraphErrors *gr1 = new TGraphErrors(npdpp,Wpp,SigExpPP,uWpp,uSigExpPP);
      gr1->SetMarkerStyle(20);
      gr1->SetMarkerColor(1);
      gr1->SetLineColor(1);
      gr1->SetMarkerSize(1.2);
      gr1->Draw("pz");         
      
      TGraph *grSigPP_CT14 = new TGraph(npfit,Wcm,sigtotPPCT14);
      grSigPP_CT14->SetLineColor(2);
      grSigPP_CT14->SetLineWidth(1);
      grSigPP_CT14->SetLineStyle(1);
      grSigPP_CT14->Draw("c"); 
      
      TGraph *grSigPBARP_CT14 = new TGraph(npfit,Wcm,sigtotPBARPCT14);
      grSigPBARP_CT14->SetLineColor(2);
      grSigPBARP_CT14->SetLineWidth(1);
      grSigPBARP_CT14->SetLineStyle(2);
      grSigPBARP_CT14->Draw("c"); 
           
      TLegend *leg1 = new TLegend(0.63,0.15,0.88,0.32);
      leg1->AddEntry(gr0,"#bar{p}p accelerator data ","p");
      leg1->AddEntry(gr1,"pp accelerator data","p");
      leg1->AddEntry(gr_cosmics,"pp cosmic rays","p");
      leg1->SetFillColor(0);
      leg1->SetShadowColor(0);
      leg1->SetBorderSize(0);
      leg1->SetTextSize(0.035);
      leg1->Draw();    

    //desenhando um miniframa em um pad dentro do canvas canv1
    TPad *pad1 = new TPad("pad1","This is pad1",0.15,0.51,0.55,0.89);
    pad1->SetFillColor(0);
    pad1->SetLogx();
    pad1->Draw();     
    pad1->cd();
    
    double xi = 5e3;
    double xf = 5e4;
    double yi = 90;
    double yf = 120;

    gr1->Draw("apz");
    gr1->SetTitle("");
    gr1->GetYaxis()->SetRangeUser(yi,yf);
    gr1->GetYaxis()->SetLabelSize(0.06);
    gr1->GetXaxis()->SetLabelSize(0.06);
    gr1->GetYaxis()->SetTitleFont(42);
    gr1->GetXaxis()->SetRangeUser(xi,xf); 
    gr_cosmics->Draw("pz");
    grSigPP_CT14->Draw("csame"); 
    grSigPBARP_CT14->Draw("csame");
      
      canv1->SaveAs("St_DGM19cCT14cut7_FFpresc_mexido.eps");    
//       canv1->Clear();
}
/*******************************************************************************/
// Plotando o parâmetro Rho pp e pbp
/*******************************************************************************/
void PlotRho(double *Wcm,double *rhoPPCT14,double *rhoPBARPCT14,int npfit)
{
    FILE *aq1,*aq2;   
// Plotting rho data        
    aq1 = fopen("paw_Rhpa.dat","r");
    aq2 = fopen("paw_Rhppc10GeV.dat","r");
    
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
     
     fclose(aq1);
     fclose(aq2);
    
     TCanvas *canv2 = new TCanvas("c1","Rho parameter",200,10,900,600);
      canv2->SetFillColor(0);
      canv2->SetLogx();
//       c1->SetGrid();
       
      TGraphErrors *gr00 = new TGraphErrors(nr1,Er1,rho1,uEr1,urho1);
      gr00->SetMarkerStyle(24);
      gr00->SetMarkerColor(1);
      gr00->SetLineColor(1);
      gr00->SetMarkerSize(1.2);
      gr00->SetTitle();
      gr00->Draw("apz");     
    
      double x_min(8e0),x_max(2.e4);
      double y_min(-0.3),y_max(0.3);
         
      gr00->GetYaxis()->SetTitle("#bf{#rho}");
      gr00->GetYaxis()->SetRangeUser(y_min,y_max);
      gr00->GetYaxis()->SetTitleOffset(1.2);
      gr00->GetYaxis()->SetTitleOffset(1.2);
      gr00->GetYaxis()->SetTitleFont(42);
      gr00->GetXaxis()->SetTitle("#bf{#sqrt{s} [GeV]}");
      gr00->GetXaxis()->SetTitleOffset(1.3);
      gr00->GetXaxis()->SetLimits(x_min,x_max);   
      gr00->GetXaxis()->SetTitleFont(42);     
     
      TGraphErrors *gr11= new TGraphErrors(nr2,Er2,rho2,uEr2,urho2);
      gr11->SetMarkerStyle(20);
      gr11->SetMarkerColor(1);
      gr11->SetLineColor(1);
      gr11->SetMarkerSize(1.2);
      gr11->SetTitle();
      gr11->Draw("pz"); 
      
      TGraph *grRhoPP_CT14 = new TGraph(npfit,Wcm,rhoPPCT14);
      grRhoPP_CT14->SetLineColor(2);
      grRhoPP_CT14->SetLineWidth(1);
      grRhoPP_CT14->SetLineStyle(1);
      grRhoPP_CT14->Draw("c"); 
      
      TGraph *grRhoPBARP_CT14 = new TGraph(npfit,Wcm,rhoPBARPCT14);
      grRhoPBARP_CT14->SetLineColor(2);
      grRhoPBARP_CT14->SetLineWidth(1);
      grRhoPBARP_CT14->SetLineStyle(2);
      grRhoPBARP_CT14->Draw("c");       
     
      TLegend *leg2 = new TLegend(0.12,0.75,0.45,0.88);
      leg2->AddEntry(gr00,"#bar{p}p accelerator data ","p");
      leg2->AddEntry(gr11,"pp accelerator data","p");
      leg2->SetShadowColor(0);
      leg2->SetBorderSize(0);
      leg2->SetTextSize(0.035);
      leg2->Draw();     

    //desenhando um miniframa em um pad dentro do canvas canv1
    TPad *pad1 = new TPad("pad1","This is pad1",0.51,0.15,0.89,0.55);
    pad1->SetFillColor(0);
    pad1->SetLogx();
    pad1->Draw();     
    pad1->cd();
    
    double xi = 5e3;
    double xf = 15e3;
    double yi = 0.07;
    double yf = 0.16;

    gr11->Draw("apz");
    gr11->SetTitle("");
    gr11->GetYaxis()->SetRangeUser(yi,yf);
    gr11->GetYaxis()->SetLabelSize(0.06);
    gr11->GetXaxis()->SetLabelSize(0.06);
    gr11->GetYaxis()->SetTitleFont(42);
    gr11->GetXaxis()->SetRangeUser(xi,xf); 
    grRhoPP_CT14->Draw("csame");
    grRhoPBARP_CT14->Draw("csame");
           
      canv2->SaveAs("Rh_DGM19cCT14cut7_FFpresc_mexido.eps");      
      canv2->Clear();
}
/***********************************************************************************************************/
void fcn(int &npar, double *gin, double &f, double *par, int iflag)
{    
   /*******************************************************/ 
    //fit parameters
    double muevsoft=par[0];
    double A=par[1];
    double B=par[2];
    double C=par[3];
    double D=par[4];
    double nu1=par[5];
    double nu2=par[6];         
    /*******************************************************/
    //reading data, then fitting

    //
    // St pp e pbp
    // 
    FILE *aq1,*aq2;
 
    // pbp   
    aq1 = fopen("paw_Stpac10GeV.dat","r");
    aq2 = fopen("paw_StppTAc10GeV.dat","r");

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
    aq2 = fopen("paw_Rhppc10GeV_mexido.dat","r");
    
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
       
    TF1 SigTotPP("SigTotPP",KerForward,bmin,bmax,10); //este 10 significa o número máximo de parâmetros 
    SigTotPP.SetParameter(0,muevsoft);
    SigTotPP.SetParameter(1,A);
    SigTotPP.SetParameter(2,B);
    SigTotPP.SetParameter(3,C);
    SigTotPP.SetParameter(4,D);
    SigTotPP.SetParameter(5,nu1);
    SigTotPP.SetParameter(6,nu2);
    SigTotPP.SetParameter(7,1);
    SigTotPP.SetParameter(8,1); 
    /*******************************************************/
    //Total XSection - pbarp
       
    TF1 SigTotPBARP("SigTotPBARP",KerForward,bmin,bmax,10);
    SigTotPBARP.SetParameter(0,muevsoft);
    SigTotPBARP.SetParameter(1,A);
    SigTotPBARP.SetParameter(2,B);
    SigTotPBARP.SetParameter(3,C);
    SigTotPBARP.SetParameter(4,D);
    SigTotPBARP.SetParameter(5,nu1);
    SigTotPBARP.SetParameter(6,nu2);
    SigTotPBARP.SetParameter(7,2);
    SigTotPBARP.SetParameter(8,1);      
    
    /*******************************************************/
    //Rho - pp  
   
    TF1 RhoPP("RhoPP",KerForward,bmin,bmax,10);
    RhoPP.SetParameter(0,muevsoft);
    RhoPP.SetParameter(1,A);
    RhoPP.SetParameter(2,B);
    RhoPP.SetParameter(3,C);
    RhoPP.SetParameter(4,D);
    RhoPP.SetParameter(5,nu1);
    RhoPP.SetParameter(6,nu2);
    RhoPP.SetParameter(7,1);
    RhoPP.SetParameter(8,2);   
     /*******************************************************/
    //Rho - pbarp  
   
    TF1 RhoPBARP("RhoPBARP",KerForward,bmin,bmax,10);
    RhoPBARP.SetParameter(0,muevsoft);
    RhoPBARP.SetParameter(1,A);
    RhoPBARP.SetParameter(2,B);
    RhoPBARP.SetParameter(3,C);
    RhoPBARP.SetParameter(4,D);
    RhoPBARP.SetParameter(5,nu1);
    RhoPBARP.SetParameter(6,nu2);
    RhoPBARP.SetParameter(7,2);
    RhoPBARP.SetParameter(8,2);   
    /*******************************************************/
    //Chi-squared calculation
    int np = 50;
    double *x=new double[np];
    double *w=new double[np];
    double sigpp,rhopp,sigpbarp,rhopbarp;     
    
    double delta;
    double chisq=0.;
    
    for(int i=0;i<npdpa;i++){
        SigTotPBARP.SetParameter(9,Wpap[i]);
        SigTotPBARP.CalcGaussLegendreSamplingPoints(np,x,w,1e-15);
        sigpbarp=SigTotPBARP.IntegralFast(np,x,w,bmin,bmax);
        delta=(SigExpPaP[i]-sigpbarp)/uSigExpPaP[i];
        chisq+=Power(delta,2);           
   }
   
    delta=0.;
    
    for(int i=0;i<npdpp;i++){
        SigTotPP.SetParameter(9,Wpp[i]);
        SigTotPP.CalcGaussLegendreSamplingPoints(np,x,w,1e-15);
        sigpp=SigTotPP.IntegralFast(np,x,w,bmin,bmax);
        delta=(SigExpPP[i]-sigpp)/uSigExpPP[i];
        chisq+=Power(delta,2);          
   }
       
    delta=0.;
       
    for(int i=0;i<nr1;i++){
        SigTotPBARP.SetParameter(9,Er1[i]);
        SigTotPBARP.CalcGaussLegendreSamplingPoints(np,x,w,1e-15);
        sigpbarp=SigTotPBARP.IntegralFast(np,x,w,bmin,bmax);
        RhoPBARP.SetParameter(9,Er1[i]);
        RhoPBARP.CalcGaussLegendreSamplingPoints(np,x,w,1e-15);
        rhopbarp=RhoPBARP.IntegralFast(np,x,w,bmin,bmax);
        rhopbarp=rhopbarp/sigpbarp;
        delta=(rho1[i]-rhopbarp)/urho1[i];
        chisq+= Power(delta,2);     
   }
   
    delta=0.;
       
    for(int i=0;i<nr2;i++){
        SigTotPP.SetParameter(9,Er2[i]);
        SigTotPP.CalcGaussLegendreSamplingPoints(np,x,w,1e-15);
        sigpp=SigTotPP.IntegralFast(np,x,w,bmin,bmax);
        RhoPP.SetParameter(9,Er2[i]);
        RhoPP.CalcGaussLegendreSamplingPoints(np,x,w,1e-15);
        rhopp=RhoPP.IntegralFast(np,x,w,bmin,bmax);
        rhopp=rhopp/sigpp;
        delta=(rho2[i]-rhopp)/urho2[i];
        chisq+= Power(delta,2);     
   } 
   
   f = chisq;        
}
/***********************************************************************************************************/

void DGM19cCT14cut7_FFpresc_mexido()
{
    /******************************************/
    //time control - start
    auto start = high_resolution_clock::now();    
    /******************************************/ 
         
    TMinuit *gMinuit = new TMinuit(numpar);  //initialize TMinuit with a maximum of 4 params
    gMinuit->SetFCN(fcn);
    
    double arglist[numpar];
    int ierflg = 0;
  
    arglist[0] = 8.18; //fitting with CL = 1\sigma
    gMinuit->mnexcm("SET ERR",arglist,1,ierflg); 
      
    double vstart[numpar] = {
    0.897727,
    123.525,
    41.6346,
    0.62000,
    24.1728,
    2.29757000,
    0.05414270}; //start values
    double step[numpar] = {1.e-3,1.e-3,1.e-3,1.e-3,1.e-3,1.e-3,1.e-3}; //steps
    gMinuit->mnparm(0, "MuSoft", vstart[0], step[0], 0.7,0.9,ierflg);
    gMinuit->mnparm(1, "A", vstart[1],   step[1], 0,0,ierflg);
    gMinuit->mnparm(2, "B", vstart[2],   step[2], 0,0,ierflg);
    gMinuit->mnparm(3, "C", vstart[3],   step[3], 0,0,ierflg);
    gMinuit->mnparm(4, "D", vstart[4],   step[4], 0,0,ierflg);
    gMinuit->mnparm(5, "Nu1", vstart[5], step[5], 0,0,ierflg);
    gMinuit->mnparm(6, "Nu2", vstart[6], step[6], 0,0,ierflg);
    
    
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
    double muevsoft=outpar[0];
    double A=outpar[1];
    double B=outpar[2];
    double C=outpar[3];
    double D=outpar[4];
    double nu1=outpar[5];
    double nu2=outpar[6];         
    /*******************************************************/
    //Total XSection - pp
       
    TF1 SigTotPP("SigTotPP",KerForward,bmin,bmax,10);
    SigTotPP.SetParameter(0,muevsoft);
    SigTotPP.SetParameter(1,A);
    SigTotPP.SetParameter(2,B);
    SigTotPP.SetParameter(3,C);
    SigTotPP.SetParameter(4,D);
    SigTotPP.SetParameter(5,nu1);
    SigTotPP.SetParameter(6,nu2);
    SigTotPP.SetParameter(7,1);
    SigTotPP.SetParameter(8,1); 
    /*******************************************************/
    //Total XSection - pbarp
       
    TF1 SigTotPBARP("SigTotPBARP",KerForward,bmin,bmax,10);
    SigTotPBARP.SetParameter(0,muevsoft);
    SigTotPBARP.SetParameter(1,A);
    SigTotPBARP.SetParameter(2,B);
    SigTotPBARP.SetParameter(3,C);
    SigTotPBARP.SetParameter(4,D);
    SigTotPBARP.SetParameter(5,nu1);
    SigTotPBARP.SetParameter(6,nu2);
    SigTotPBARP.SetParameter(7,2);
    SigTotPBARP.SetParameter(8,1);    
     /*******************************************************/
    //Rho - pp  
   
    TF1 RhoPP("RhoPP",KerForward,bmin,bmax,10);
    RhoPP.SetParameter(0,muevsoft);
    RhoPP.SetParameter(1,A);
    RhoPP.SetParameter(2,B);
    RhoPP.SetParameter(3,C);
    RhoPP.SetParameter(4,D);
    RhoPP.SetParameter(5,nu1);
    RhoPP.SetParameter(6,nu2);
    RhoPP.SetParameter(7,1);
    RhoPP.SetParameter(8,2);   
     /*******************************************************/
    //Rho - pbarp  
   
    TF1 RhoPBARP("RhoPBARP",KerForward,bmin,bmax,10);
    RhoPBARP.SetParameter(0,muevsoft);
    RhoPBARP.SetParameter(1,A);
    RhoPBARP.SetParameter(2,B);
    RhoPBARP.SetParameter(3,C);
    RhoPBARP.SetParameter(4,D);
    RhoPBARP.SetParameter(5,nu1);
    RhoPBARP.SetParameter(6,nu2);
    RhoPBARP.SetParameter(7,2);
    RhoPBARP.SetParameter(8,2);      
/*************************************************************************************************/
   //calculates total xsection and rho, creating files to make the plots
    FILE *aq1,*aq2,*aq3,*aq4;
   
    aq1=fopen("Stpp_DGM19cCT14cut7_FFpresc_mexido.dat","w"); 
    aq2=fopen("Rhpp_DGM19cCT14cut7_FFpresc_mexido.dat","w"); 
    aq3=fopen("Stpa_DGM19cCT14cut7_FFpresc_mexido.dat","w"); 
    aq4=fopen("Rhpa_DGM19cCT14cut7_FFpresc_mexido.dat","w");  
       
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
            SigTotPP.SetParameter(9,Ecm);
            RhoPP.SetParameter(9,Ecm);
            SigTotPP.CalcGaussLegendreSamplingPoints(np,x,w,1e-15);
            RhoPP.CalcGaussLegendreSamplingPoints(np,x,w,1e-15);
            sigpp = SigTotPP.IntegralFast(np,x,w,bmin,bmax);
            rhopp = RhoPP.IntegralFast(np,x,w,bmin,bmax);
            rhopp=rhopp/sigpp;
            fprintf(aq1,"%.2e   %.4lf\n",Ecm,sigpp);             
            fprintf(aq2,"%.2e   %.4lf\n",Ecm,rhopp);
            /*************************************************/
            //pbarp
            SigTotPBARP.SetParameter(9,Ecm);
            RhoPBARP.SetParameter(9,Ecm);
            SigTotPBARP.CalcGaussLegendreSamplingPoints(np,x,w,1e-15);
            RhoPBARP.CalcGaussLegendreSamplingPoints(np,x,w,1e-15);
            sigpbarp = SigTotPBARP.IntegralFast(np,x,w,bmin,bmax);
            rhopbarp = RhoPBARP.IntegralFast(np,x,w,bmin,bmax);
            rhopbarp=rhopbarp/sigpbarp;
            fprintf(aq3,"%.2e   %.4lf\n",Ecm,sigpbarp);             
            fprintf(aq4,"%.2e   %.4lf\n",Ecm,rhopbarp); 
            /*************************************************/                             
            i+=di;
       }while(Ecm<Wmax); 
       
      fclose(aq1);
      fclose(aq2);
      fclose(aq3);
      fclose(aq4);  
       
      aq1 = fopen("Stpp_DGM19cCT14cut7_FFpresc_mexido.dat","r");
      aq2 = fopen("Rhpp_DGM19cCT14cut7_FFpresc_mexido.dat","r");
      aq3 = fopen("Stpa_DGM19cCT14cut7_FFpresc_mexido.dat","r");
      aq4 = fopen("Rhpa_DGM19cCT14cut7_FFpresc_mexido.dat","r");
      
      const int npfit=108;
      double Wcm[npfit],sigtotPPCT14[npfit],rhoPPCT14[npfit],sigtotPBARPCT14[npfit],rhoPBARPCT14[npfit];       
            
      for(int j=0;j<npfit;j++){
          fscanf(aq1,"%lg %lg",&Wcm[j],&sigtotPPCT14[j]);
          fscanf(aq2,"%lg %lg",&Wcm[j],&rhoPPCT14[j]);
          fscanf(aq3,"%lg %lg",&Wcm[j],&sigtotPBARPCT14[j]);
          fscanf(aq4,"%lg %lg",&Wcm[j],&rhoPBARPCT14[j]);
      }
      
      fclose(aq1);
      fclose(aq2);
      fclose(aq3);
      fclose(aq4);      
       
    /************************************************************************************************/
    //Plotting total xsection data         
    PlotTotXSec(Wcm,sigtotPPCT14,sigtotPBARPCT14,npfit);      
    /************************************************************************************************/
    //Plotting rho data     
     PlotRho(Wcm,rhoPPCT14,rhoPBARPCT14,npfit);     
      
    /******************************************/
    auto stop = high_resolution_clock::now(); 
    auto duration = duration_cast<seconds>(stop - start); 
    cout << "Time(in min): "<<duration.count()/6.e1 << endl;
    /******************************************/ 
        
}   
    
    
    
    
    
    
    
    
    
    
    
