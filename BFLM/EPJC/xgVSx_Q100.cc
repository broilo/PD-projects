/********************************************************************************/

void xgVSx_Q100()
{
   
   //cria ponteiro do tipo TCanvas para desenhar o gráfico 
     
      TCanvas *canv1 = new TCanvas("c1","xPDF",200,10,900,600);
      canv1->SetFillColor(0);
      canv1->SetLogx(); //escala log em x
   
    FILE *aq1,*aq2,*aq3,*aq4,*aq5,*aq6;
    /********************************************************************************/
    //configurações dos eixos x e y do gráfico (títulos, fontes, etc)
      
       double x_min(0.1),x_max(1.0);
       double y_min(-0.05),y_max(0.5);
      
     // double x_min(1e1),x_max(2.e4);
     // double y_min(3.e1),y_max(1.2e2);
 
    /********************************************************************************/
    //lê e armazena curvas de SigTot pp e pbarp calculadas com PDF CT14
    
    //Mateus: aqui vc deve substituir meus arquivos pelos seus resultados, lembrando de criar mais ponteiros do tipo "FILE" para ler os resultados da CTEQ e MMHT 
    aq1 = fopen("xgVSx_CTEQ6L_Q100.dat","r");
    aq2 = fopen("xgVSx_CT14_Q100.dat","r");
    aq3 = fopen("xgVSx_MMHT_Q100.dat","r");
           
      const int npfit=77; //esse é o número de pontos que eu calculei para a minha curva, caso vc tenha um número diferente, deves alterar.
      
      double xx[npfit],CTEQ6L_Q100[npfit],CT14_Q100[npfit],MMHT_Q100[npfit];       
            
      for(int j=0;j<npfit;j++){
          fscanf(aq1,"%lg %lg",&xx[j],&CTEQ6L_Q100[j]);
          fscanf(aq2,"%lg %lg",&xx[j],&CT14_Q100[j]); 
          fscanf(aq3,"%lg %lg",&xx[j],&MMHT_Q100[j]);
      }
      
      fclose(aq1);
      fclose(aq2);
      fclose(aq3);
     
    /********************************************************************************/
    
    /********************************************************************************/
    //cria ponteiro do tipo TGraph para desenhar curvas 

      
      TGraph *grCTEQ6L_Q100 = new TGraph(npfit,xx,CTEQ6L_Q100);
      grCTEQ6L_Q100->SetTitle();
      grCTEQ6L_Q100->Draw("apz");
      grCTEQ6L_Q100->SetLineColor(1);
      grCTEQ6L_Q100->SetLineWidth(1);
      grCTEQ6L_Q100->SetLineStyle(1);
      grCTEQ6L_Q100->Draw("c"); 
      
      TGraph *grCT14_Q100 = new TGraph(npfit,xx,CT14_Q100);
      grCT14_Q100->SetLineColor(4);
      grCT14_Q100->SetLineWidth(1);
      grCT14_Q100->SetLineStyle(3);
      grCT14_Q100->Draw("c");
      
      TGraph *grMMHT_Q100 = new TGraph(npfit,xx,MMHT_Q100);
      grMMHT_Q100->SetLineColor(2);
      grMMHT_Q100->SetLineWidth(1);
      grMMHT_Q100->SetLineStyle(2);
      grMMHT_Q100->Draw("c");

      grCTEQ6L_Q100->GetYaxis()->SetTitle("#bf{xg(x,Q^{2})}");
      grCTEQ6L_Q100->GetYaxis()->SetRangeUser(y_min,y_max);
      grCTEQ6L_Q100->GetYaxis()->SetTitleOffset(1.2);
      grCTEQ6L_Q100->GetYaxis()->SetTitleFont(42);
      grCTEQ6L_Q100->GetXaxis()->SetTitle("#bf{x}");
      grCTEQ6L_Q100->GetXaxis()->SetTitleOffset(1.3);
      grCTEQ6L_Q100->GetXaxis()->SetLimits(x_min,x_max);   
      grCTEQ6L_Q100->GetXaxis()->SetTitleFont(42);  

   TLatex *latex_1 = new TLatex();
    
   latex_1->SetText(0.5,0.25,"#bf{Q^{2}= 10^{4} GeV^{2}}");
   
   latex_1->SetTextFont(42);
   latex_1->SetTextSize(0.035);
   latex_1->Draw();

    /********************************************************************************/
    //cria legenda do gráfico    
      
      TLegend *leg1 = new TLegend(0.68,0.60,0.88,0.82);
      leg1->AddEntry(grCTEQ6L_Q100,"CTEQ6L ","l");
      leg1->AddEntry(grCT14_Q100,"CT14 ","l");
      leg1->AddEntry(grMMHT_Q100,"MMHT ","l");
//      leg1->AddEntry(grSigPBARP_CT14,"CT14LO - #bar{p}p","l");
      leg1->SetFillColor(0);
      leg1->SetShadowColor(0);
      leg1->SetBorderSize(0);
      leg1->SetTextSize(0.028);
      leg1->Draw();
      
      
    /********************************************************************************/
    //cria arquivo .eps com o gráfico de saída
      
      canv1->SaveAs("xgVSxQ100v3.eps"); 
    /********************************************************************************/

}
