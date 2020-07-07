/********************************************************************************/
const double Norm=1.e-1;
const double NormA=0.035;

void xgVSx_Q1_3()
{
   
   //cria ponteiro do tipo TCanvas para desenhar o gráfico 
     
      TCanvas *canv1 = new TCanvas("c1","xPDF",200,10,900,600);
      canv1->SetFillColor(0);
      canv1->SetLogx(); //escala log em x
   
    FILE *aq1,*aq2,*aq3,*aq4,*aq5,*aq6;
    /********************************************************************************/
    //configurações dos eixos x e y do gráfico (títulos, fontes, etc)
      
       double x_min(0.00000001),x_max(1.0);
       double y_min(0.0),y_max(10.5);
      
     // double x_min(1e1),x_max(2.e4);
     // double y_min(3.e1),y_max(1.2e2);
 
    /********************************************************************************/
    //lê e armazena curvas de SigTot pp e pbarp calculadas com PDF CT14
    
    //Mateus: aqui vc deve substituir meus arquivos pelos seus resultados, lembrando de criar mais ponteiros do tipo "FILE" para ler os resultados da CTEQ e MMHT 
    aq1 = fopen("xgVSx_CTEQ6L_Q1_3v2.dat","r");
    aq2 = fopen("xgVSx_CT14_Q1_3v2.dat","r");
    aq3 = fopen("xgVSx_MMHT_Q1_3v2.dat","r");
           
      const int npfit=50000; //10000 //esse é o número de pontos que eu calculei para a minha curva, caso vc tenha um número diferente, deves alterar.
      
      double xx[npfit],CTEQ6L_Q1_3[npfit],CT14_Q1_3[npfit],MMHT_Q1_3[npfit];       
            
      for(int j=0;j<npfit;j++){
          fscanf(aq1,"%lg %lg",&xx[j],&CTEQ6L_Q1_3[j]);
          fscanf(aq2,"%lg %lg",&xx[j],&CT14_Q1_3[j]); 
          fscanf(aq3,"%lg %lg",&xx[j],&MMHT_Q1_3[j]);
          
          CTEQ6L_Q1_3[j]*=Norm;
          CT14_Q1_3[j]*=Norm;
          MMHT_Q1_3[j]*=Norm;
      }
      
      fclose(aq1);
      fclose(aq2);
      fclose(aq3);
     
    /********************************************************************************/
    
    /********************************************************************************/
    //cria ponteiro do tipo TGraph para desenhar curvas 

      
      TGraph *grCTEQ6L_Q1_3 = new TGraph(npfit,xx,CTEQ6L_Q1_3);
      grCTEQ6L_Q1_3->SetTitle();
      grCTEQ6L_Q1_3->Draw("apz");
      grCTEQ6L_Q1_3->SetLineColor(2);
      grCTEQ6L_Q1_3->SetLineWidth(2);
      grCTEQ6L_Q1_3->SetLineStyle(1);
      grCTEQ6L_Q1_3->Draw("c"); 
      
      TGraph *grCT14_Q1_3 = new TGraph(npfit,xx,CT14_Q1_3);
      grCT14_Q1_3->SetLineColor(1);
      grCT14_Q1_3->SetLineWidth(2);
      grCT14_Q1_3->SetLineStyle(1);
      grCT14_Q1_3->Draw("c");
      
      TGraph *grMMHT_Q1_3 = new TGraph(npfit,xx,MMHT_Q1_3);
      grMMHT_Q1_3->SetLineColor(4);
      grMMHT_Q1_3->SetLineWidth(2);
      grMMHT_Q1_3->SetLineStyle(1);
      grMMHT_Q1_3->Draw("c");

      grCTEQ6L_Q1_3->GetYaxis()->SetTitle("#bf{xg(x,Q^{2}) [#times 10^{-1}]}");
      grCTEQ6L_Q1_3->GetYaxis()->SetRangeUser(y_min,y_max);
      grCTEQ6L_Q1_3->GetYaxis()->SetTitleOffset(1.2);
      grCTEQ6L_Q1_3->GetYaxis()->SetTitleFont(42);
      grCTEQ6L_Q1_3->GetXaxis()->SetTitle("#bf{x}");
      grCTEQ6L_Q1_3->GetXaxis()->SetTitleOffset(1.3);
      grCTEQ6L_Q1_3->GetXaxis()->SetLimits(x_min,x_max);   
      grCTEQ6L_Q1_3->GetXaxis()->SetTitleFont(42);  

   TLatex *latex_1 = new TLatex();
    
   latex_1->SetText(0.008,5.5,"#bf{Q= 1.3 GeV}");
   
   latex_1->SetTextFont(42);
   latex_1->SetTextSize(0.035);
   latex_1->Draw();

    /********************************************************************************/
    //cria legenda do gráfico    
      
      TLegend *leg1 = new TLegend(0.68,0.60,0.88,0.82);
      leg1->AddEntry(grCT14_Q1_3,"CT14 ","l");
      leg1->AddEntry(grCTEQ6L_Q1_3,"CTEQ6L ","l");
      leg1->AddEntry(grMMHT_Q1_3,"MMHT ","l");
//      leg1->AddEntry(grSigPBARP_CT14,"CT14LO - #bar{p}p","l");
      leg1->SetFillColor(0);
      leg1->SetShadowColor(0);
      leg1->SetBorderSize(0);
      leg1->SetTextSize(0.028);
      leg1->Draw();
      
      
    /********************************************************************************/
    //cria arquivo .eps com o gráfico de saída
      
      canv1->SaveAs("xgVSxQ1_3.eps"); 
    /********************************************************************************/

}
