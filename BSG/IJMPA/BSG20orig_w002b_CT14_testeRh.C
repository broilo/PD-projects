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
    Obs3: PDF: CT14
    Obs4: Free: w, Sig0, A1, A2, del1, del2
    Obs5: B is fixed in 1.0

    *** Important *** 
    (See: https://arxiv.org/abs/1908.01040 Chapter 8, sections 8.4 and 8.8)

    mg=0.4d0 !GeV -> mg=500+/-200 MeV for lmb=300 MeV

    --> (MSTW) lbd=0.318920 GeV lambda equivalent to alpha_LO(Mz^2)=0.13939~0.139 macthing-prescrition squeme; m_b=4.18 GeV

    --> (MMHT) lbd=0.267471 GeV lambda equivalent to alpha_LO(Mz^2)=0.13499~0.135 macthing-prescrition squeme; m_b=4.18 GeV

    --> (CTEQ6L) lbd=0.3260 GeV lambda equivalent to alpha_NLO(Mz^2)=0.118 for 4 flavors

    --> (CT14)   lbd=0.3260 GeV lambda equivalent to alpha_NLO(Mz^2)=0.118 for 4 flavors

    ------------------------------------------------------------------------

    *** Notice that ***
    This is a temp code to fit XSecTot from BSG20orig_w002b.C (CT14 analysis) by
    by means of a Regge-Gribov parametrization (Phys.Lett.B 781 (2018) 616-620, 
    SciHub link: https://sci-hub.tw/https://doi.org/10.1016/j.physletb.2018.04.045).

    ------------------------------------------------------------------------


    *** How to run the code *** 

    root -l -q /usr/lib/x86_64-linux-gnu/libgsl.so /usr/lib/x86_64-linux-gnu/libgslcblas.so BSG20orig_w002b_CT14_testeRh.C | cat > BSG20orig_w002b_CT14_testeRh.min

    ---------------------------------------------------------------------
    Codes by: M.Broilo*, P.V.R.G.Silva, V.P.B.Gonçalves
    mateus.broilo@ufrgs.br | mateus.broilo90@gmail.com
    High and Medium Energy Group
    Grupo de Altas e Medias Energias (GAME)
    Universidade Federal de Pelotas, Pelotas - RS (Brasil)
    ---------------------------------------------------------------------

    Creation: 18/dez/2019 (Pelotas, RS, Brazil)
    Last update: 27/abr/2020 (Pelotas, RS, Brazil)
   
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

//----------------------------------------------------------------------
//constant parameters throughout the model
const double s0=25.;

//----------------------------------------------------------------------
//number of fit parameters and data points
const int numpar=6;
const int npSigTotPBP=59;
const int npSigTotPP=118;

const int npSigTotPPcurve = 108;
const int npSigTotPBPcurve = 108;
const int nppfit = npSigTotPPcurve;
const int npapfit = npSigTotPBPcurve;
const int npMin=npapfit+nppfit;

// Plots
const int npSigTotPPplot=104;
const int npSigTotPPTOTEM=12;
const int npSigTotPPATLAS=2;

const int npRppPlot=64;
const int npRppTOTEM=4;

//----------------------------------------------------------------------
//parameters to control energy range (for plots and integrations)
const double Wmin=5.;
const double Wmax=1.e6;

//----------------------------------------------------------------------
/*
                              Analysis 
*/
//----------------------------------------------------------------------
//Forward observables - total xsection
double KerForward(double *x,double *par)
{
    double W=x[0];      
    double a1=par[0]; 
    double b1=par[1];
    double a2=par[2];
    double b2=par[3];
    double A=par[4];
    double B=par[5]; 
    double reac=par[6];     
        
    double s=W*W;
    double Y=Log(s/s0);
    double Y2=Power(Y,2.);
    double tau;
    /*******************************************************/
    // C-even and C-Odd
    /*******************************************************/
    double SigEven = a1*Power((s/s0),-b1) + A*Y + B*Y2;
    double SigOdd = a2*Power((s/s0),-b2);
    double SigTot;
    /*******************************************************/    
    double KerSig;
    
    	if(reac==1) //pp
    	{
        SigTot = SigEven - SigOdd;//tau = -1.;
    	}
    	else if(reac==2) //pbp
    	{
        SigTot = SigEven + SigOdd;//tau = 1.;
    	}
    	else
    	{
    	tau=0.;
    	}
            KerSig=SigTot*mbfactor; //in mb 
            return KerSig;
}

//----------------------------------------------------------------------
// Plotando a seção de choque total pp e pbp
//----------------------------------------------------------------------
void PlotTotXSec(double *Wcm,double *XSecTOTpp,double *XSecTOTpbp,int npfit)
{
    FILE *aq1,*aq2,*aq3,*aq4,*aq5,*aq6;

    // --------- Plotting SigTot data ---------
        aq1 = fopen("paw_Stpa.dat","r");
        aq2 = fopen("paw_Stpp.dat","r");
        aq3 = fopen("paw_StTOTEM.dat","r");
        aq4 = fopen("paw_StATLAS.dat","r");

        //---> pbp <---
        const int npStpa = npSigTotPBP;
        double W21[npStpa],Sig21Exp[npStpa],uW21[npStpa],uSig21Exp[npStpa];
    
        for(int i=0;i<npStpa;i++)     
        fscanf(aq1,"%lg %lg %lg %lg",&W21[i],&Sig21Exp[i],&uSig21Exp[i],&uW21[i]);
            
        //---> pp <---     
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
    
    // --------- Plotting the real-result curve ---------

        aq5 = fopen("paramet_Stpp_BSG20orig_w002b.dat","r");//XSecTotPP_BSG20orig_w002b.dat
        aq6 = fopen("paramet_Stpbp_BSG20orig_w002b.dat","r");//XSecTotPBP_BSG20orig_w002b.dat    

    //---> The real-result <--- 
        
        //pp 
        const int npResPP = npSigTotPPcurve;
        double W11res[npResPP],Sig11res[npResPP],uW11res[npResPP],uSig11res[npResPP];  
        
        for(int i=0;i<npResPP;i++)   
        fscanf(aq5,"%lg %lg %lg %lg",&W11res[i],&Sig11res[i],&uSig11res[i],&uW11res[i]);
            
        //pbp 
        const int npResPBP = npSigTotPBPcurve;
        double W21res[npResPBP],Sig21res[npResPBP],uW21res[npResPBP],uSig21res[npResPBP];  
        
        for(int i=0;i<npResPBP;i++)   
        fscanf(aq6,"%lg %lg %lg %lg",&W21res[i],&Sig21res[i],&uSig21res[i],&uW21res[i]);

        fclose(aq1);
        fclose(aq2);
        fclose(aq3);
        fclose(aq4);
        fclose(aq5);
        fclose(aq6);

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

        TGraphErrors *gr4 = new TGraphErrors(npResPP,W11res,Sig11res,uW11res,uSig11res);
            gr4->SetMarkerStyle(24);
            gr4->SetMarkerColor(13);
            gr4->SetLineColor(1);
            gr4->SetMarkerSize(1.0);
            gr4->SetTitle();
            gr4->Draw("pz");
        
        TGraphErrors *gr5 = new TGraphErrors(npResPBP,W21res,Sig21res,uW21res,uSig21res);
            gr5->SetMarkerStyle(24);
            gr5->SetMarkerColor(13);
            gr5->SetLineColor(1);
            gr5->SetMarkerSize(1.0);
            gr5->Draw("pz");

        TGraph *grSigTotPBARP = new TGraph(npfit,Wcm,XSecTOTpbp);
            grSigTotPBARP->SetLineColor(4);
            grSigTotPBARP->SetLineWidth(1);
            grSigTotPBARP->SetLineStyle(1);
            grSigTotPBARP->Draw("c");
        
        TGraph *grSigTotPP = new TGraph(npfit,Wcm,XSecTOTpp);
            grSigTotPP->SetLineColor(2);
            grSigTotPP->SetLineWidth(1);
            grSigTotPP->SetLineStyle(1);
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

            canv1->SaveAs("XSECtot_BSG20orig_w002b_CT14_testeRh.eps");    
}

//----------------------------------------------------------------------
void fcn(int &npar, double *gin, double &f, double *par, int iflag)
{    
    /*******************************************************/ 
    //fit parameters
    double a1=par[0]; 
    double b1=par[1];
    double a2=par[2];
    double b2=par[3];
    double A=par[4];  
    double B=par[5];
    /*******************************************************/
    //reading data, then fitting

    //
    // St pp e pbp
    // 
    FILE *aq1,*aq2;
 
    // pbp   
    aq1 = fopen("paramet_Stpbp_BSG20orig_w002b.dat","r");
    aq2 = fopen("paramet_Stpp_BSG20orig_w002b.dat","r");    

    const int npdpa = npSigTotPBPcurve;
    double W21res[npdpa],Sig21res[npdpa],uW21res[npdpa],uSig21res[npdpa];
   
    for(int i=0;i<npdpa;i++)     
        fscanf(aq1,"%lg %lg %lg %lg",&W21res[i],&Sig21res[i],&uSig21res[i],&uW21res[i]);
        
     // pp   
     const int npdpp = npSigTotPPcurve;
     double W11res[npdpp],Sig11res[npdpp],uW11res[npdpp],uSig11res[npdpp];
     
     for(int i=0;i<npdpp;i++)   
        fscanf(aq2,"%lg %lg %lg %lg",&W11res[i],&Sig11res[i],&uSig11res[i],&uW11res[i]);
      
     fclose(aq1);
     fclose(aq2);

    /*******************************************************/
    //Total XSection - pp
       
    TF1 SigTotPP("SigTotPP",KerForward,Wmin,Wmax,8);
    SigTotPP.SetParameter(0,a1);
    SigTotPP.SetParameter(1,b1);
    SigTotPP.SetParameter(2,a2);
    SigTotPP.SetParameter(3,b2);
    SigTotPP.SetParameter(4,A);
    SigTotPP.SetParameter(5,B);
    SigTotPP.SetParameter(6,1);   
    
    /*******************************************************/
    //Total XSection - pbp
       
    TF1 SigTotPBARP("SigTotPBARP",KerForward,Wmin,Wmax,8);
    SigTotPBARP.SetParameter(0,a1);
    SigTotPBARP.SetParameter(1,b1);
    SigTotPBARP.SetParameter(2,a2);
    SigTotPBARP.SetParameter(3,b2);
    SigTotPBARP.SetParameter(4,A);
    SigTotPBARP.SetParameter(5,B);
    SigTotPBARP.SetParameter(6,2); 
        
    /*******************************************************/

    //Chi-squared calculation
    int np = 50;
    double *x=new double[np];
    double *w=new double[np];
    double Stpp,Stpbarp;    
    
    double delta;
    double chisq=0.;

    for(int i=0;i<npdpp;i++){
        Stpp=SigTotPP.Eval(W11res[i]);
        delta=(Sig11res[i]-Stpp)/uSig11res[i];
        chisq+=Power(delta,2);          
    }
       
    delta=0.;

    for(int i=0;i<npdpa;i++){
        Stpbarp=SigTotPBARP.Eval(W21res[i]);
        delta=(Sig21res[i]-Stpbarp)/uSig21res[i];
        chisq+=Power(delta,2);           
    }
   
    delta=0.; 
       
    f = chisq;        
}

//----------------------------------------------------------------------
void BSG20orig_w002b_CT14_testeRh()
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
        58.8,
        0.231,
        17.1,
        0.548,
        3.76,
        0.122
    }; //start values */
    double step[numpar] = {1.e-3,1.e-3,1.e-3,1.e-3,1.e-3,1.e-3}; //steps
    gMinuit->mnparm(0, "a1",    vstart[0],   step[0], 0,0,ierflg);
    gMinuit->mnparm(1, "b1",    vstart[1],   step[1], 0,0,ierflg);
    gMinuit->mnparm(2, "a2",    vstart[2],   step[2], 0,0,ierflg);
    gMinuit->mnparm(3, "b2",    vstart[3],   step[3], 0,0,ierflg);
    gMinuit->mnparm(4, "A",     vstart[4],   step[4], 0,0,ierflg);
    gMinuit->mnparm(5, "B",     vstart[5],   step[5], 0,0,ierflg);

    //start minimizing data
    arglist[0] = 1;
    gMinuit->mnexcm("SET STR", arglist ,1,ierflg);
    arglist[0] = 500;    
    gMinuit->mnexcm("MIGRAD", arglist ,1,ierflg);    
    gMinuit->mnexcm("MIGRAD", arglist ,1,ierflg);
    gMinuit->mnexcm("MIGRAD", arglist ,1,ierflg);
    //gMinuit->mnexcm("MIGRAD", arglist ,1,ierflg);   
   
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
    double a1=outpar[0]; 
    double b1=outpar[1];
    double a2=outpar[2];
    double b2=outpar[3];
    double A=outpar[4];  
    double B=outpar[5];
    /*******************************************************/
    //Total XSection - pp
       
    TF1 SigTotPP("SigTotPP",KerForward,Wmin,Wmax,8);
    SigTotPP.SetParameter(0,a1);
    SigTotPP.SetParameter(1,b1);
    SigTotPP.SetParameter(2,a2);
    SigTotPP.SetParameter(3,b2);
    SigTotPP.SetParameter(4,A);
    SigTotPP.SetParameter(5,B);
    SigTotPP.SetParameter(6,1);   
    
    /*******************************************************/
    //Total XSection - pbp
       
    TF1 SigTotPBARP("SigTotPBARP",KerForward,Wmin,Wmax,8);
    SigTotPBARP.SetParameter(0,a1);
    SigTotPBARP.SetParameter(1,b1);
    SigTotPBARP.SetParameter(2,a2);
    SigTotPBARP.SetParameter(3,b2);
    SigTotPBARP.SetParameter(4,A);
    SigTotPBARP.SetParameter(5,B);
    SigTotPBARP.SetParameter(6,2); 
        
    /*******************************************************/   
    //----------------------------------------------------------------------
    //calculates total creating files to make the plots
    FILE *aq1,*aq2;
 
    aq1=fopen("XSecTotPP_BSG20orig_w002b_CT14_testeRh.dat","w"); 
    aq2=fopen("XSecTotPBP_BSG20orig_w002b_CT14_testeRh.dat","w");
           
     //number of directive computing for GaussLegendreIntegration
    int np = 50;
    double *x=new double[np];
    double *w=new double[np];
    
    double i=0.;
    const double di=0.05;
    double Ecm;
    double Stpp,Stpbarp;
    
     do{
            Ecm=5.*pow(10.,i);
            /*************************************************/
            //pp
            Stpp = SigTotPP.Eval(Ecm);   

            fprintf(aq1,"%.2e   %.4lf\n",Ecm,Stpp); 
            /*************************************************/
            //pbarp
            Stpbarp = SigTotPBARP.Eval(Ecm);
     
            fprintf(aq2,"%.2e   %.4lf\n",Ecm,Stpbarp);   
            /*************************************************/      
                                   
            i+=di;
                       
       //cout << Ecm << "\t" << Seikpp << "\t" << Seikpbarp << endl;
            
        }while(Ecm<Wmax); 
       
    fclose(aq1);
    fclose(aq2);      
    
    aq1 = fopen("XSecTotPP_BSG20orig_w002b_CT14_testeRh.dat","r");
    aq2 = fopen("XSecTotPBP_BSG20orig_w002b_CT14_testeRh.dat","r");            
    
    const int npfit=108;
    double Wcm[npfit],XSecTOTpp[npfit],XSecTOTpbp[npfit];       
        
    for(int j=0;j<npfit;j++)
        {
          fscanf(aq1,"%lg %lg",&Wcm[j],&XSecTOTpp[j]);
          fscanf(aq2,"%lg %lg",&Wcm[j],&XSecTOTpbp[j]);  
        }
      
    fclose(aq1);
    fclose(aq2);
             
    //----------------------------------------------------------------------
    //Plotting xsection data         
    PlotTotXSec(Wcm,XSecTOTpp,XSecTOTpbp,npfit);   
    /******************************************/
    auto stop = high_resolution_clock::now(); 
    auto duration = duration_cast<seconds>(stop - start); 
    cout << "Time(in min): " << duration.count()/6.e1 << endl;
    /******************************************/         
}