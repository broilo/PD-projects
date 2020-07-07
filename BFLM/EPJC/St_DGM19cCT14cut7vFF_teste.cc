/********************************************************************************/
//número de dados pp e pbarp
const int npSap=30;
const int npSpp=80;
/********************************************************************************/
void St_DGM19cCT14cut7vFF_teste()
{
   
    FILE *aq1,*aq2,*aq3,*aq4,*aq5,*aq6;
    /********************************************************************************/
    //lê e armazena os dados (de aceleradores) de SigTot pp e pbarp acima de 10 GeV

    aq1 = fopen("paw_Stpac10GeV.dat","r");
    aq2 = fopen("paw_StppTAc10GeV.dat","r");

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
    /********************************************************************************/

    /********************************************************************************/
    //cria ponteiro do tipo TCanvas para desenhar o gráfico 
     
      TCanvas *canv1 = new TCanvas("c1","Total Cross Section",200,10,900,600);
      canv1->SetFillColor(0);
      canv1->SetLogx(); //escala log em x
    /********************************************************************************/
    //cria ponteiro do tipo TGraphError para desenhar os dados de SigTot pbarp
   
      TGraphErrors *gr0 = new TGraphErrors(npdpa,Wpap,SigExpPaP,uWpap,uSigExpPaP);
      gr0->SetMarkerStyle(24);
      gr0->SetMarkerColor(1);
      gr0->SetLineColor(1);
      gr0->SetMarkerSize(1.2);
      gr0->SetTitle();
      gr0->Draw("apz");
    /********************************************************************************/    
    //configurações dos eixos x e y do gráfico (títulos, fontes, etc)
      
       double x_min(8e00),x_max(2.e05);
       double y_min(3.0e01),y_max(1.8e02);
      
     // double x_min(1e1),x_max(2.e4);
     // double y_min(3.e1),y_max(1.2e2);
       
      gr0->GetYaxis()->SetTitle("#bf{#sigma_{tot} [mb]}");
      gr0->GetYaxis()->SetRangeUser(y_min,y_max);
      gr0->GetYaxis()->SetTitleOffset(1.2);
      gr0->GetYaxis()->SetTitleFont(42);
      gr0->GetXaxis()->SetTitle("#bf{#sqrt{s} [GeV]}");
      gr0->GetXaxis()->SetTitleOffset(1.3);
      gr0->GetXaxis()->SetLimits(x_min,x_max);   
      gr0->GetXaxis()->SetTitleFont(42);    
    /********************************************************************************/   
      //dados de raios cósmicos (pontos com erros assimétricos, que não entram no fit)  
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
    /********************************************************************************/   
   
    /********************************************************************************/   
    //cria ponteiro do tipo TGraphError para desenhar os dados de SigTot pp
      TGraphErrors *gr1 = new TGraphErrors(npdpp,Wpp,SigExpPP,uWpp,uSigExpPP);
      gr1->SetMarkerStyle(20);
      gr1->SetMarkerColor(1);
      gr1->SetLineColor(1);
      gr1->SetMarkerSize(1.2);
      gr1->Draw("pz");     
    /********************************************************************************/   
    
    /********************************************************************************/
    //lê e armazena curvas de SigTot pp e pbarp calculadas com PDF CT14
    
    //Mateus: aqui vc deve substituir meus arquivos pelos seus resultados, lembrando de criar mais ponteiros do tipo "FILE" para ler os resultados da CTEQ e MMHT 
      aq1 = fopen("Stpp_DGM19cCT14cut7vFF_teste.dat","r"); 
      aq2 = fopen("Stpa_DGM19cCT14cut7vFF_teste.dat","r");
          
      const int npfit=160; //esse é o número de pontos que eu calculei para a minha curva, caso vc tenha um número diferente, deves alterar.
      
      double Wcm[npfit],sigtotPPCT14[npfit],sigtotPBARPCT14[npfit],sigtotPPCTEQ6L[npfit],sigtotPBARPCTEQ6L[npfit],sigtotPPMMHT[npfit],sigtotPBARPMMHT[npfit];       
            
      for(int j=0;j<npfit;j++){
          fscanf(aq1,"%lg %lg",&Wcm[j],&sigtotPPCT14[j]);
          fscanf(aq2,"%lg %lg",&Wcm[j],&sigtotPBARPCT14[j]);       
      }
      
      fclose(aq1);
      fclose(aq2);
    /********************************************************************************/
    
    /********************************************************************************/
    //cria ponteiro do tipo TGraph para desenhar curvas de SigTot pp

      TGraph *grSigPP_CT14 = new TGraph(npfit,Wcm,sigtotPPCT14);
      grSigPP_CT14->SetLineColor(2);
      grSigPP_CT14->SetLineWidth(1);
      grSigPP_CT14->SetLineStyle(1);
      grSigPP_CT14->Draw("c"); 

    /********************************************************************************/
  
    /********************************************************************************/
    //cria ponteiro do tipo TGraph para desenhar curvas de SigTot pbarp  
   
      TGraph *grSigPBARP_CT14 = new TGraph(npfit,Wcm,sigtotPBARPCT14);
      grSigPBARP_CT14->SetLineColor(2);
      grSigPBARP_CT14->SetLineWidth(1);
      grSigPBARP_CT14->SetLineStyle(2);
      grSigPBARP_CT14->Draw("c");

    /********************************************************************************/
    //cria legenda do gráfico    
      
      TLegend *leg1 = new TLegend(0.63,0.15,0.88,0.32);
      leg1->AddEntry(gr0,"#bar{p}p accelerator data ","p");
      leg1->AddEntry(gr1,"pp accelerator data","p");
      leg1->AddEntry(gr_cosmics,"pp cosmic rays","p");
//      leg1->AddEntry(grSigPP_CT14,"CT14 ","l");
//      leg1->AddEntry(grSigPBARP_CT14,"CT14LO - #bar{p}p","l");
      leg1->SetFillColor(0);
      leg1->SetShadowColor(0);
      leg1->SetBorderSize(0);
      leg1->SetTextSize(0.035);
      leg1->Draw();
      
    /********************************************************************************/
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

    /********************************************************************************/
    //cria arquivo .eps com o gráfico de saída
      
      canv1->SaveAs("St_DGM19cCT14cut7vFF_teste.eps"); 
    /********************************************************************************/

}
