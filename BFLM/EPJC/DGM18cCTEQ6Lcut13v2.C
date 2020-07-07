/*
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
c Obs10: Foi adicionado St(13TeV)= 110.3 +/- 3.5.
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

using namespace TMath;

const double hc2=0.389379; //mbGeV^{2}

/********************************************************************************/
//constant parameters throughout the model
const double s0=25.; 
const double muoddsoft=0.5;
const double Kf=1.; 
/********************************************************************************/
//number of fit parameters and data points 
const int numpar = 7;
const int npRap=12;
const int npRpp=52;
const int npSap=30;
const int npSpp=80;
const int npapfit=npRap+npSap;
const int nppfit=npRpp+npSpp;
const int npMin=npapfit+nppfit;
/********************************************************************************/
//parameters to control energy and b range (for plots and integrations)
const double Wmin=5.;
const double Wmax=1.e6;
const double bmin=0.;
const double bmax=30.;   
/*******************************************************************************/
//QCD cross section - Mateus' analytical parametrization for CTEQ6LLO PDF
/*******************************************************************************/
const double b1=97.005;
const double b2=0.28009e-1;
const double b3=1.6986;
const double b4=1.7362;
const double b5=-0.14947e-5;
const double b6=14.140;
const double b7=0.31972;
const double b8=0.83551e-1;
const double b9=3.8127;
const double b10=0.80990;   
const double c1=1.01;
const double c2=1.05;
const double c3=1.09;  
const double bkg=100.;

double BesselK3(double z)
{
   if(z>0.)
       return BesselK(3,z);
   else
       return 0.; 
    
}

/*******************************************************************************/
//Sigma QCD - Complex
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
//semi-hard
double WSH(double b,double W,double nu1,double nu2)
{    
    double s=W*W;
    double Y=Log(s/s0);
    double mush=nu1-nu2*Y;
    
    double ws=Power(mush,5)*Power(b,3)*BesselK3(mush*b)/(96.*Pi());  
    
    return ws;    
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
    
    double W=par[0];
    double muevsoft=par[1];
    double A=par[2];
    double B=par[3];
    double C=par[4];
    double D=par[5];
    double nu1=par[6];
    double nu2=par[7];    
    int reac=par[8]; 
    int obs=par[9]; 
        
    double s=W*W;
    double Y=Log(s/s0);
    double Y2=Y*Y;      
    /*******************************************************/
    //Form factors
    /*******************************************************/
    double WSh=WSH(b,W,nu1,nu2);
    double WEvS=WSoft(b,muevsoft); 
    double WOddS=WSoft(b,muoddsoft); 
    /*******************************************************/
    //Eikonal pp - real part
    /*******************************************************/
    double resig=ReSigQCD(W);
    double ReChiEvenSH=0.5*Kf*WSh*resig;   
    /*******************************************************/
    double ReChiEvenSoft=0.5*WEvS*(A+(B*Cos(PiOver4())/Sqrt(s/s0))+C*(Y2-Pi()*PiOver4()));     
    /*******************************************************/
    double ReChiOddSoft=0.5*WOddS*D*Cos(PiOver4())/Sqrt(s/s0);     
    /*******************************************************/
    double ReChiEven=ReChiEvenSH+ReChiEvenSoft;
    double ReChiOdd=ReChiOddSoft;         
    /*******************************************************/
    //Eikonal pp - imaginary part
    /*******************************************************/    
    double imsig=ImSigQCD(W);
    double ImChiEvenSH=0.5*Kf*WSh*imsig;
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
      
        if(obs==1)//total xsection
        {      
            KerSig=4.*Pi()*hc2*b*(1.-Cos(ImChi)*Exp(-ReChi)); //in mb      
            return KerSig;
        
        }else if(obs==2)//rho
        {
           KerRho=-4.*Pi()*hc2*b*Sin(ImChi)*Exp(-ReChi); //in mb    
           return KerRho;
        }else{
           return 0.; 
        }      
    
}

void PlotTotXSec(double *Wcm,double *sigtotPPCTEQ6L,double *sigtotPBARPCTEQ6L,int npfit)
{
    FILE *aq1,*aq2;
    
    aq1 = fopen("paw_Stpa.dat","r");
    aq2 = fopen("paw_Stpp.dat","r");

    const int npdpa = npSap;
    double Wpap[npdpa],SigExpPaP[npdpa],uWpap[npdpa],uSigExpPaP[npdpa];
   
    for(int i=0;i<npdpa;i++)     
        fscanf(aq1,"%lg %lg %lg %lg",&Wpap[i],&SigExpPaP[i],&uSigExpPaP[i],&uWpap[i]);
     
     const int npdpp = npSpp;
     double Wpp[npdpp],SigExpPP[npdpp],uWpp[npdpp],uSigExpPP[npdpp];
     
     for(int i=0;i<npdpp;i++)   
        fscanf(aq2,"%lg %lg %lg %lg",&Wpp[i],&SigExpPP[i],&uSigExpPP[i],&uWpp[i]);
       
     fclose(aq1);
     fclose(aq2);
     
      TCanvas *canv1 = new TCanvas("c1","Total Cross Section",200,10,900,600);
      canv1->SetFillColor(0);
      canv1->SetLogx();
       
      TGraphErrors *gr0 = new TGraphErrors(npdpa,Wpap,SigExpPaP,uWpap,uSigExpPaP);
      gr0->SetMarkerStyle(20);
      gr0->SetMarkerColor(2);
      gr0->SetLineColor(1);
      gr0->SetMarkerSize(1.1);
      gr0->SetTitle();
      gr0->Draw("apz");
      
      double x_min(1e1),x_max(6.e4);
      double y_min(3.e1),y_max(1.5e2);
       
      gr0->GetYaxis()->SetTitle("#bf{#sigma_{tot} [mb]}");
      gr0->GetYaxis()->SetRangeUser(y_min,y_max);
      gr0->GetYaxis()->SetTitleOffset(1.2);
      gr0->GetYaxis()->SetTitleFont(42);
      gr0->GetXaxis()->SetTitle("#bf{#sqrt{s} [GeV]}");
      gr0->GetXaxis()->SetTitleOffset(1.3);
      gr0->GetXaxis()->SetLimits(x_min,x_max);   
      gr0->GetXaxis()->SetTitleFont(42);    
       
      const int npt_cosmics = 13;
      double x_cosmics[npt_cosmics] = {86.453220,137.000283,193.510517,273.332639,386.084553,6164.111067,8125.766425,10711.911104,14120.962078,18615.117942,24539.416782,30000.154666,57e03};
      double y_cosmics[npt_cosmics] = {43.000,48.000,54.000,56.000,54.000,93.000,101.00,117.00,104.00,100.00,124.00,120.00,133.00};
      double erx_lowucosmics[npt_cosmics] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
      double erx_uppucosmics[npt_cosmics] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
      double ery_lowucosmics[npt_cosmics] = {5.968718,6.923698,7.161257,8.306528,9.033825,14.000000,16.000000,18.000000,26.000000,27.000000,34.000000,15.000000,28.7};
      double ery_uppucosmics[npt_cosmics] = {5.968718,6.923698,7.161257,8.306528,9.033825,14.000000,16.000000,18.000000,26.000000,27.000000,34.000000,15.000000,26.7};
      
      TGraphAsymmErrors *gr_cosmics = new TGraphAsymmErrors(npt_cosmics,x_cosmics,y_cosmics,erx_lowucosmics,erx_uppucosmics,ery_lowucosmics,ery_uppucosmics);
      gr_cosmics->SetMarkerStyle(22);
      gr_cosmics->SetMarkerSize(1);
      gr_cosmics->SetMarkerColorAlpha(kGray,0.8);
      gr_cosmics->SetLineColorAlpha(kGray,0.8);
      gr_cosmics->Draw("pz"); 
      
      TGraphErrors *gr1 = new TGraphErrors(npdpp,Wpp,SigExpPP,uWpp,uSigExpPP);
      gr1->SetMarkerStyle(20);
      gr1->SetMarkerColor(4);
      gr1->SetLineColor(1);
      gr1->SetMarkerSize(1.1);
      gr1->Draw("pz");     

      
      TGraph *grSigPP_CTEQ6L = new TGraph(npfit,Wcm,sigtotPPCTEQ6L);
      grSigPP_CTEQ6L->SetLineColor(1);
      grSigPP_CTEQ6L->SetLineWidth(1);
      grSigPP_CTEQ6L->SetLineStyle(1);
      grSigPP_CTEQ6L->Draw("c"); 
      
      TGraph *grSigPBARP_CTEQ6L = new TGraph(npfit,Wcm,sigtotPBARPCTEQ6L);
      grSigPBARP_CTEQ6L->SetLineColor(1);
      grSigPBARP_CTEQ6L->SetLineWidth(1);
      grSigPBARP_CTEQ6L->SetLineStyle(2);
      grSigPBARP_CTEQ6L->Draw("c"); 
      
      
      TLegend *leg1 = new TLegend(0.15,0.58,0.48,0.88);
      leg1->AddEntry(gr0,"#bar{p}p accelerator data ","p");
      leg1->AddEntry(gr1,"pp accelerator data","p");
      leg1->AddEntry(gr_cosmics,"pp cosmic rays","p");
      leg1->AddEntry(grSigPP_CTEQ6L,"CTEQ6L - pp","l");
      leg1->AddEntry(grSigPBARP_CTEQ6L,"CTEQ6L - #bar{p}p","l");
      leg1->SetFillColor(0);
      leg1->SetShadowColor(0);
      leg1->SetBorderSize(0);
      leg1->SetTextSize(0.028);
      leg1->Draw();    
      
      canv1->SaveAs("St_DGM18cCTEQ6Lcut13v2.eps");    
//       canv1->Clear();
}

void PlotRho(double *Wcm,double *rhoPPCTEQ6L,double *rhoPBARPCTEQ6L,int npfit)
{
    FILE *aq1,*aq2;
    
  //Plotting rho data        
    aq1 = fopen("paw_Rhpa.dat","r");
    aq2 = fopen("paw_Rhpp.dat","r");
    
    //pbarp
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
    
     TCanvas *canv2 = new TCanvas("c1","Rho parameter",200,10,900,600);
      canv2->SetFillColor(0);
      canv2->SetLogx();
//       c1->SetGrid();
       
      TGraphErrors *gr00 = new TGraphErrors(nr1,Er1,rho1,uEr1,urho1);
      gr00->SetMarkerStyle(20);
      gr00->SetMarkerColor(2);
      gr00->SetLineColor(1);
      gr00->SetMarkerSize(1.1);
      gr00->SetTitle();
      gr00->Draw("apz");      
    
      double x_min(1e1),x_max(6.e4);
      double y_min(-0.2),y_max(0.25);
         
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
      gr11->SetMarkerColor(4);
      gr11->SetLineColor(1);
      gr11->SetMarkerSize(1.1);
      gr11->SetTitle();
      gr11->Draw("pz"); 
      
      TGraph *grRhoPP_CTEQ6L = new TGraph(npfit,Wcm,rhoPPCTEQ6L);
      grRhoPP_CTEQ6L->SetLineColor(1);
      grRhoPP_CTEQ6L->SetLineWidth(1);
      grRhoPP_CTEQ6L->SetLineStyle(1);
      grRhoPP_CTEQ6L->Draw("c"); 
      
      TGraph *grRhoPBARP_CTEQ6L = new TGraph(npfit,Wcm,rhoPBARPCTEQ6L);
      grRhoPBARP_CTEQ6L->SetLineColor(1);
      grRhoPBARP_CTEQ6L->SetLineWidth(1);
      grRhoPBARP_CTEQ6L->SetLineStyle(2);
      grRhoPBARP_CTEQ6L->Draw("c");       
     
      TLegend *leg2 = new TLegend(0.58,0.15,0.88,0.48);
      leg2->AddEntry(gr00,"#bar{p}p accelerator data ","p");
      leg2->AddEntry(gr11,"pp accelerator data","p");
      leg2->AddEntry(grRhoPP_CTEQ6L,"CTEQ6L - pp","l");
      leg2->AddEntry(grRhoPBARP_CTEQ6L,"CTEQ6L - #bar{p}p","l");
      leg2->SetFillColor(0);
      leg2->SetShadowColor(0);
      leg2->SetBorderSize(0);
      leg2->SetTextSize(0.028);
      leg2->Draw();        
      
      canv2->SaveAs("Rh_DGM18cCTEQ6Lcut13v2.eps");      
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
    FILE *aq1,*aq2;
    
    aq1 = fopen("paw_Stpa.dat","r");
    aq2 = fopen("paw_Stpp.dat","r");

    const int npdpa = npSap;
    double Wpap[npdpa],SigExpPaP[npdpa],uWpap[npdpa],uSigExpPaP[npdpa];
   
    for(int i=0;i<npdpa;i++)     
        fscanf(aq1,"%lg %lg %lg %lg",&Wpap[i],&SigExpPaP[i],&uSigExpPaP[i],&uWpap[i]);
     
     const int npdpp = npSpp;
     double Wpp[npdpp],SigExpPP[npdpp],uWpp[npdpp],uSigExpPP[npdpp];
     
     for(int i=0;i<npdpp;i++)   
        fscanf(aq2,"%lg %lg %lg %lg",&Wpp[i],&SigExpPP[i],&uSigExpPP[i],&uWpp[i]);
       
     fclose(aq1);
     fclose(aq2);   
           
    aq1 = fopen("paw_Rhpa.dat","r");
    aq2 = fopen("paw_Rhpp.dat","r");
    
    //pbarp
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
       
    TF1 SigTotPP("SigTotPP",KerForward,bmin,bmax,10);
    SigTotPP.SetParameter(1,muevsoft);
    SigTotPP.SetParameter(2,A);
    SigTotPP.SetParameter(3,B);
    SigTotPP.SetParameter(4,C);
    SigTotPP.SetParameter(5,D);
    SigTotPP.SetParameter(6,nu1);
    SigTotPP.SetParameter(7,nu2);
    SigTotPP.SetParameter(8,1);
    SigTotPP.SetParameter(9,1); 
    /*******************************************************/
    //Total XSection - pbarp
       
    TF1 SigTotPBARP("SigTotPBARP",KerForward,bmin,bmax,10);
    SigTotPBARP.SetParameter(1,muevsoft);
    SigTotPBARP.SetParameter(2,A);
    SigTotPBARP.SetParameter(3,B);
    SigTotPBARP.SetParameter(4,C);
    SigTotPBARP.SetParameter(5,D);
    SigTotPBARP.SetParameter(6,nu1);
    SigTotPBARP.SetParameter(7,nu2);
    SigTotPBARP.SetParameter(8,2);
    SigTotPBARP.SetParameter(9,1);      
    
    /*******************************************************/
    //Rho - pp  
   
    TF1 RhoPP("RhoPP",KerForward,bmin,bmax,10);
    RhoPP.SetParameter(1,muevsoft);
    RhoPP.SetParameter(2,A);
    RhoPP.SetParameter(3,B);
    RhoPP.SetParameter(4,C);
    RhoPP.SetParameter(5,D);
    RhoPP.SetParameter(6,nu1);
    RhoPP.SetParameter(7,nu2);
    RhoPP.SetParameter(8,1);
    RhoPP.SetParameter(9,2);   
     /*******************************************************/
    //Rho - pbarp  
   
    TF1 RhoPBARP("RhoPBARP",KerForward,bmin,bmax,10);
    RhoPBARP.SetParameter(1,muevsoft);
    RhoPBARP.SetParameter(2,A);
    RhoPBARP.SetParameter(3,B);
    RhoPBARP.SetParameter(4,C);
    RhoPBARP.SetParameter(5,D);
    RhoPBARP.SetParameter(6,nu1);
    RhoPBARP.SetParameter(7,nu2);
    RhoPBARP.SetParameter(8,2);
    RhoPBARP.SetParameter(9,2);   
    /*******************************************************/
    //Chi-squared calculation
    int np = 50;
    double *x=new double[np];
    double *w=new double[np];
    double sigpp,rhopp,sigpbarp,rhopbarp;     
    
    double delta;
    double chisq=0.;
    
    for(int i=0;i<npdpa;i++){
        SigTotPBARP.SetParameter(0,Wpap[i]);
        SigTotPBARP.CalcGaussLegendreSamplingPoints(np,x,w,1e-15);
        sigpbarp=SigTotPBARP.IntegralFast(np,x,w,bmin,bmax);
        delta=(SigExpPaP[i]-sigpbarp)/uSigExpPaP[i];
        chisq+=Power(delta,2);       
    
   }
   
    delta=0.;
    
    for(int i=0;i<npdpp;i++){
        SigTotPP.SetParameter(0,Wpp[i]);
        SigTotPP.CalcGaussLegendreSamplingPoints(np,x,w,1e-15);
        sigpp=SigTotPP.IntegralFast(np,x,w,bmin,bmax);
        delta=(SigExpPP[i]-sigpp)/uSigExpPP[i];
        chisq+=Power(delta,2);       
    
   }
       
    delta=0.;
       
    for(int i=0;i<nr1;i++){
        SigTotPBARP.SetParameter(0,Er1[i]);
        SigTotPBARP.CalcGaussLegendreSamplingPoints(np,x,w,1e-15);
        sigpbarp=SigTotPBARP.IntegralFast(np,x,w,bmin,bmax);
        RhoPBARP.SetParameter(0,Er1[i]);
        RhoPBARP.CalcGaussLegendreSamplingPoints(np,x,w,1e-15);
        rhopbarp=RhoPBARP.IntegralFast(np,x,w,bmin,bmax);
        rhopbarp=rhopbarp/sigpbarp;
        delta=(rho1[i]-rhopbarp)/urho1[i];
        chisq+= Power(delta,2);     
   }
   
    delta=0.;
       
    for(int i=0;i<nr2;i++){
        SigTotPP.SetParameter(0,Er2[i]);
        SigTotPP.CalcGaussLegendreSamplingPoints(np,x,w,1e-15);
        sigpp=SigTotPP.IntegralFast(np,x,w,bmin,bmax);
        RhoPP.SetParameter(0,Er2[i]);
        RhoPP.CalcGaussLegendreSamplingPoints(np,x,w,1e-15);
        rhopp=RhoPP.IntegralFast(np,x,w,bmin,bmax);
        rhopp=rhopp/sigpp;
        delta=(rho2[i]-rhopp)/urho2[i];
        chisq+= Power(delta,2);     
   } 
   
   f = chisq;   
     
}
/***********************************************************************************************************/

void DGM18cCTEQ6Lcut13v2()
{
    time_t t0,tf;  //variáveis para controle de tempo de execucao
    double texec;
    time(&t0); //tempo de iní­cio do programa
         
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
    SigTotPP.SetParameter(1,muevsoft);
    SigTotPP.SetParameter(2,A);
    SigTotPP.SetParameter(3,B);
    SigTotPP.SetParameter(4,C);
    SigTotPP.SetParameter(5,D);
    SigTotPP.SetParameter(6,nu1);
    SigTotPP.SetParameter(7,nu2);
    SigTotPP.SetParameter(8,1);
    SigTotPP.SetParameter(9,1); 
    /*******************************************************/
    //Total XSection - pbarp
       
    TF1 SigTotPBARP("SigTotPBARP",KerForward,bmin,bmax,10);
    SigTotPBARP.SetParameter(1,muevsoft);
    SigTotPBARP.SetParameter(2,A);
    SigTotPBARP.SetParameter(3,B);
    SigTotPBARP.SetParameter(4,C);
    SigTotPBARP.SetParameter(5,D);
    SigTotPBARP.SetParameter(6,nu1);
    SigTotPBARP.SetParameter(7,nu2);
    SigTotPBARP.SetParameter(8,2);
    SigTotPBARP.SetParameter(9,1);  
   
     /*******************************************************/
    //Rho - pp  
   
    TF1 RhoPP("RhoPP",KerForward,bmin,bmax,10);
    RhoPP.SetParameter(1,muevsoft);
    RhoPP.SetParameter(2,A);
    RhoPP.SetParameter(3,B);
    RhoPP.SetParameter(4,C);
    RhoPP.SetParameter(5,D);
    RhoPP.SetParameter(6,nu1);
    RhoPP.SetParameter(7,nu2);
    RhoPP.SetParameter(8,1);
    RhoPP.SetParameter(9,2);   
     /*******************************************************/
    //Rho - pbarp  
   
    TF1 RhoPBARP("RhoPBARP",KerForward,bmin,bmax,10);
    RhoPBARP.SetParameter(1,muevsoft);
    RhoPBARP.SetParameter(2,A);
    RhoPBARP.SetParameter(3,B);
    RhoPBARP.SetParameter(4,C);
    RhoPBARP.SetParameter(5,D);
    RhoPBARP.SetParameter(6,nu1);
    RhoPBARP.SetParameter(7,nu2);
    RhoPBARP.SetParameter(8,2);
    RhoPBARP.SetParameter(9,2);   
    
   /*************************************************************************************************/
   //calculates total xsection and rho, creating files to make the plots
    FILE *aq1,*aq2,*aq3,*aq4;
   
    aq1=fopen("Stpp_DGM18cCTEQ6Lcut13v2.dat","w"); 
    aq2=fopen("Rhpp_DGM18cCTEQ6Lcut13v2.dat","w"); 
    aq3=fopen("Stpa_DGM18cCTEQ6Lcut13v2.dat","w"); 
    aq4=fopen("Rhpa_DGM18cCTEQ6Lcut13v2.dat","w");  
       
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
            SigTotPP.SetParameter(0,Ecm);
            RhoPP.SetParameter(0,Ecm);
            SigTotPP.CalcGaussLegendreSamplingPoints(np,x,w,1e-15);
            RhoPP.CalcGaussLegendreSamplingPoints(np,x,w,1e-15);
            sigpp = SigTotPP.IntegralFast(np,x,w,bmin,bmax);
            rhopp = RhoPP.IntegralFast(np,x,w,bmin,bmax);
            rhopp=rhopp/sigpp;
            fprintf(aq1,"%.2e   %.4lf\n",Ecm,sigpp);             
            fprintf(aq2,"%.2e   %.4lf\n",Ecm,rhopp);
            /*************************************************/
            //pbarp
            SigTotPBARP.SetParameter(0,Ecm);
            RhoPBARP.SetParameter(0,Ecm);
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
       
      aq1 = fopen("Stpp_DGM18cCTEQ6Lcut13v2.dat","r");
      aq2 = fopen("Rhpp_DGM18cCTEQ6Lcut13v2.dat","r");
      aq3 = fopen("Stpa_DGM18cCTEQ6Lcut13v2.dat","r");
      aq4 = fopen("Rhpa_DGM18cCTEQ6Lcut13v2.dat","r");
      
      const int npfit=108;
      double Wcm[npfit],sigtotPPCTEQ6L[npfit],rhoPPCTEQ6L[npfit],sigtotPBARPCTEQ6L[npfit],rhoPBARPCTEQ6L[npfit];       
            
      for(int j=0;j<npfit;j++){
          fscanf(aq1,"%lg %lg",&Wcm[j],&sigtotPPCTEQ6L[j]);
          fscanf(aq2,"%lg %lg",&Wcm[j],&rhoPPCTEQ6L[j]);
          fscanf(aq3,"%lg %lg",&Wcm[j],&sigtotPBARPCTEQ6L[j]);
          fscanf(aq4,"%lg %lg",&Wcm[j],&rhoPBARPCTEQ6L[j]);
      }
      
      fclose(aq1);
      fclose(aq2);
      fclose(aq3);
      fclose(aq4);      
       
    /************************************************************************************************/
    //Plotting total xsection data         
    PlotTotXSec(Wcm,sigtotPPCTEQ6L,sigtotPBARPCTEQ6L,npfit);      
    /************************************************************************************************/
    //Plotting rho data     
     PlotRho(Wcm,rhoPPCTEQ6L,rhoPBARPCTEQ6L,npfit);     
        
    
   time(&tf); //tempo de término do programa
   texec=difftime (tf,t0);
   printf("Tempo de execucao do programa (em segundos): \n\%.2lf\n",texec); 
        
}   
    
    
    
    
    
    
    
    
    
    
    
