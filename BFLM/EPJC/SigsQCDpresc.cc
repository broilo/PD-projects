/********************************************************************************/
//número de dados pp e pbarp
const int np=30;
/********************************************************************************/

void SigsQCDpresc()
{
   
    FILE *aq1,*aq2,*aq3,*aq4,*aq5,*aq6,*aq7,*aq8,*aq9;
    /********************************************************************************/
    //lê e armazena os dados de SigQCD

    aq1 = fopen("sig_sumij_real_CTEQ6Llo_m400q13V001.dat","r");
    aq2 = fopen("sig_sumij_real_CT14lo_m400q13V001.dat","r");
    aq3 = fopen("sig_sumij_real_MMHTlo_m400q13V001.dat","r");

    const int npCTEQ6L = np;
    double WCTEQ6L[npCTEQ6L],SigCTEQ6L[npCTEQ6L],uSigCTEQ6L[npCTEQ6L];
   
    for(int i=0;i<npCTEQ6L;i++)     
        fscanf(aq1,"%lg %lg %lg",&WCTEQ6L[i],&SigCTEQ6L[i],&uSigCTEQ6L[i]);
        
    const int npCT14 = np;
    double WCT14[npCT14],SigCT14[npCT14],uSigCT14[npCT14];
   
    for(int i=0;i<npCT14;i++)     
        fscanf(aq2,"%lg %lg %lg",&WCT14[i],&SigCT14[i],&uSigCT14[i]);
        
    const int npMMHT = np;
    double WMMHT[npMMHT],SigMMHT[npMMHT],uSigMMHT[npMMHT];
   
    for(int i=0;i<npMMHT;i++)     
        fscanf(aq3,"%lg %lg %lg",&WMMHT[i],&SigMMHT[i],&uSigMMHT[i]);    
       
     fclose(aq1);
     fclose(aq2);
     fclose(aq3);
     
    /********************************************************************************/
    //cria ponteiro do tipo TCanvas para desenhar o gráfico 
     
      TCanvas *canv1 = new TCanvas("c1","QCD Cross Section",200,10,900,600);
      canv1->SetFillColor(0);
      canv1->SetLogx(); //escala log em x
    /********************************************************************************/
    //cria ponteiro do tipo TGraphError para desenhar os dados de SigQCD
   
      TGraphErrors *gr0 = new TGraphErrors(npCTEQ6L,WCTEQ6L,SigCTEQ6L,uSigCTEQ6L);
      gr0->SetMarkerStyle(20);
      gr0->SetMarkerColor(1);
      gr0->SetLineColor(1);
      gr0->SetMarkerSize(0.8);
      gr0->SetTitle();
      gr0->Draw("apz");
      
      /********************************************************************************/    
      //configurações dos eixos x e y do gráfico (títulos, fontes, etc)
      
      double x_min(8e00),x_max(1.e05);
      double y_min(0.),y_max(15e3);
       
      gr0->GetYaxis()->SetTitle("#bf{#sigma_{QCD}(s) [GeV^{-2}]}");
      gr0->GetYaxis()->SetRangeUser(y_min,y_max);
      gr0->GetYaxis()->SetTitleOffset(1.5);
      gr0->GetYaxis()->SetTitleFont(42);
      gr0->GetXaxis()->SetTitle("#bf{#sqrt{s} [GeV]}");
      gr0->GetXaxis()->SetTitleOffset(1.3);
      gr0->GetXaxis()->SetLimits(x_min,x_max);   
      gr0->GetXaxis()->SetTitleFont(42); 
      
      /********************************************************************************/  
      
      TGraphErrors *gr1 = new TGraphErrors(npCT14,WCT14,SigCT14,uSigCT14);
      gr1->SetMarkerStyle(20);
      gr1->SetMarkerColor(2);
      gr1->SetLineColor(1);
      gr1->SetMarkerSize(0.8);
      gr1->SetTitle();
      gr1->Draw("pz");
      
      TGraphErrors *gr2 = new TGraphErrors(npMMHT,WMMHT,SigMMHT,uSigMMHT);
      gr2->SetMarkerStyle(20);
      gr2->SetMarkerColor(4);
      gr2->SetLineColor(1);
      gr2->SetMarkerSize(0.8);
      gr2->SetTitle();
      gr2->Draw("pz");        

      /********************************************************************************/   
      //lê e armazena curvas de SigQCD
 
      aq4 = fopen("gerReSigQCDpresc_CTEQ6Llo.dat","r"); 
      aq5 = fopen("gerImSigQCDpresc_CTEQ6Llo.dat","r");
      aq6 = fopen("gerReSigQCDpresc_CT14lo.dat","r"); 
      aq7 = fopen("gerImSigQCDpresc_CT14lo.dat","r");
      aq8 = fopen("gerReSigQCDpresc_MMHTlo.dat","r"); 
      aq9 = fopen("gerImSigQCDpresc_MMHTlo.dat","r");
                 
      const int npfit=10001; //esse é o número de pontos que eu calculei para a minha curva, caso vc tenha um número diferente, deves alterar.
      
      double wReCTEQ6L[npfit],ResigCTEQ6L[npfit],
             wImCTEQ6L[npfit],ImsigCTEQ6L[npfit],
             wReCT14[npfit],ResigCT14[npfit],
             wImCT14[npfit],ImsigCT14[npfit],
             wReMMHT[npfit],ResigMMHT[npfit],
             wImMMHT[npfit],ImsigMMHT[npfit];       
            
      for(int j=0;j<npfit;j++){
          fscanf(aq4,"%lg %lg",&wReCTEQ6L[j],&ResigCTEQ6L[j]);
          fscanf(aq5,"%lg %lg",&wImCTEQ6L[j],&ImsigCTEQ6L[j]);
          fscanf(aq6,"%lg %lg",&wReCT14[j],&ResigCT14[j]);
          fscanf(aq7,"%lg %lg",&wImCT14[j],&ImsigCT14[j]);
          fscanf(aq8,"%lg %lg",&wReMMHT[j],&ResigMMHT[j]);
          fscanf(aq9,"%lg %lg",&wImMMHT[j],&ImsigMMHT[j]);
                 
      }
      
      fclose(aq4);
      fclose(aq5);
      fclose(aq6);
      fclose(aq7);
      fclose(aq8);
      fclose(aq9);

     /********************************************************************************/
     //cria ponteiro do tipo TGraph para desenhar curvas de SigQCD
           
      TGraph *grReSigCTEQ6L = new TGraph(npfit,wReCTEQ6L,ResigCTEQ6L);
      grReSigCTEQ6L->SetLineColor(1);
      grReSigCTEQ6L->SetLineWidth(1);
      grReSigCTEQ6L->SetLineStyle(1);
      grReSigCTEQ6L->Draw("c"); 
      
      TGraph *grImSigCTEQ6L = new TGraph(npfit,wImCTEQ6L,ImsigCTEQ6L);
      grImSigCTEQ6L->SetLineColor(1);
      grImSigCTEQ6L->SetLineWidth(1);
      grImSigCTEQ6L->SetLineStyle(2);
      grImSigCTEQ6L->Draw("c"); 
      
      TGraph *grReSigCT14 = new TGraph(npfit,wReCT14,ResigCT14);
      grReSigCT14->SetLineColor(2);
      grReSigCT14->SetLineWidth(1);
      grReSigCT14->SetLineStyle(1);
      grReSigCT14->Draw("c"); 
      
      TGraph *grImSigCT14 = new TGraph(npfit,wImCT14,ImsigCT14);
      grImSigCT14->SetLineColor(2);
      grImSigCT14->SetLineWidth(1);
      grImSigCT14->SetLineStyle(2);
      grImSigCT14->Draw("c");
      
      TGraph *grReSigMMHT = new TGraph(npfit,wReMMHT,ResigMMHT);
      grReSigMMHT->SetLineColor(4);
      grReSigMMHT->SetLineWidth(1);
      grReSigMMHT->SetLineStyle(1);
      grReSigMMHT->Draw("c"); 
      
      TGraph *grImSigMMHT = new TGraph(npfit,wImMMHT,ImsigMMHT);
      grImSigMMHT->SetLineColor(4);
      grImSigMMHT->SetLineWidth(1);
      grImSigMMHT->SetLineStyle(2);
      grImSigMMHT->Draw("c");

      /********************************************************************************/
      //Add texto no plot
      
      TLatex *latex_1 = new TLatex();   
      latex_1->SetText(20,8e3,"#bf{Complex Parametrization}");
      latex_1->SetTextFont(42);
      latex_1->SetTextSize(0.030);
      latex_1->Draw();
      
      TLatex *latex_2 = new TLatex();   
      latex_2->SetText(20,7e3,"#bf{Q_{min} = 1.3 GeV}");
      latex_2->SetTextFont(42);
      latex_2->SetTextSize(0.030);
      latex_2->Draw();
      
      TLatex *latex_3 = new TLatex();   
      latex_3->SetText(20,6e3,"#bf{m_{g} = 400 MeV}");
      latex_3->SetTextFont(42);
      latex_3->SetTextSize(0.030);
      latex_3->Draw();

      /********************************************************************************/
      //cria legenda do gráfico    
      
      TLegend *leg1 = new TLegend(0.15,0.58,0.36,0.88);
/*
      leg1->AddEntry(grReSigCTEQ6L,"(CTEQ6L) Re #sigma_{QCD}(s) Parametrized ","l");
      leg1->AddEntry(grImSigCTEQ6L,"(CTEQ6L) -Im #sigma_{QCD}(s) s #rightarrow -is ","l");
      leg1->AddEntry(grReSigCT14,"(CT14) Re #sigma_{QCD}(s) Parametrized ","l");
      leg1->AddEntry(grImSigCT14,"(CT14) -Im #sigma_{QCD}(s) s #rightarrow -is ","l");
      leg1->AddEntry(grReSigMMHT,"(MMHT) Re #sigma_{QCD}(s) Parametrized ","l");
      leg1->AddEntry(grImSigMMHT,"(MMHT) -Im #sigma_{QCD}(s) s #rightarrow -is ","l");
      leg1->SetFillColor(0);
      leg1->SetShadowColor(0);
      leg1->SetBorderSize(0);
      leg1->SetTextSize(0.028);
      leg1->Draw();
*/
      leg1->AddEntry(grReSigCTEQ6L,"(CTEQ6L) Re #sigma_{QCD}(s)","l");
      leg1->AddEntry(grImSigCTEQ6L,"(CTEQ6L) -Im #sigma_{QCD}(s)","l");
      leg1->AddEntry(grReSigCT14,"(CT14) Re #sigma_{QCD}(s)","l");
      leg1->AddEntry(grImSigCT14,"(CT14) -Im #sigma_{QCD}(s)","l");
      leg1->AddEntry(grReSigMMHT,"(MMHT) Re #sigma_{QCD}(s)","l");
      leg1->AddEntry(grImSigMMHT,"(MMHT) -Im #sigma_{QCD}(s)","l");
      leg1->SetFillColor(0);
      leg1->SetShadowColor(0);
      leg1->SetBorderSize(0);
      leg1->SetTextSize(0.028);
      leg1->Draw();
 
    
      
    /********************************************************************************/
    //cria arquivo .eps com o gráfico de saída
      
      canv1->SaveAs("SigsQCDprescV4.eps"); 
    /********************************************************************************/

}
