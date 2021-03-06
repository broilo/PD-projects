/********************************************************************************/
//número de dados pp e pbarp
const int npRap=12;
const int npRpp=52;
/********************************************************************************/
void Rh_DGM19cCT14cut7vFF_teste()
{
    FILE *aq1,*aq2,*aq3,*aq4,*aq5,*aq6;
    /********************************************************************************/
    //lê e armazena os dados (de aceleradores) de SigTot pp e pbarp acima de 10 GeV
         
    aq1 = fopen("paw_Rhpa.dat","r");
    aq2 = fopen("paw_Rhppc10GeV.dat","r");
    
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
    /********************************************************************************/

    /********************************************************************************/
    //cria ponteiro do tipo TCanvas para desenhar o gráfico 
    
     TCanvas *canv2 = new TCanvas("c1","Rho parameter",200,10,900,600);
      canv2->SetFillColor(0);
      canv2->SetLogx();
//       c1->SetGrid();
    /********************************************************************************/
    //cria ponteiro do tipo TGraphError para desenhar os dados de Rho pbarp
       
      TGraphErrors *gr00 = new TGraphErrors(nr1,Er1,rho1,uEr1,urho1);
      gr00->SetMarkerStyle(24);
      gr00->SetMarkerColor(1);
      gr00->SetLineColor(1);
      gr00->SetMarkerSize(1.2);
      gr00->SetTitle();
      gr00->Draw("apz");
    /********************************************************************************/    
    //configurações dos eixos x e y do gráfico (títulos, fontes, etc)      
    
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
     /********************************************************************************/   
   
     /********************************************************************************/   
     //cria ponteiro do tipo TGraphError para desenhar os dados de Rho pp
      TGraphErrors *gr11= new TGraphErrors(nr2,Er2,rho2,uEr2,urho2);
      gr11->SetMarkerStyle(20);
      gr11->SetMarkerColor(1);
      gr11->SetLineColor(1);
      gr11->SetMarkerSize(1.2);
      gr11->SetTitle();
      gr11->Draw("pz"); 
      /********************************************************************************/   
    
     /********************************************************************************/
     //lê e armazena curvas de Rho pp e pbarp calculadas com PDF CT14
    
     //Mateus: aqui vc deve substituir meus arquivos pelos seus resultados, lembrando de criar mais ponteiros do tipo "FILE" para ler os resultados da CTEQ e MMHT 
      aq1 = fopen("Rhpp_DGM19cCT14cut7vFF_teste.dat","r"); 
      aq2 = fopen("Rhpa_DGM19cCT14cut7vFF_teste.dat","r");
            
      const int npfit=90;//esse é o número de pontos que eu calculei para a minha curva, caso vc tenha um número diferente, deves alterar.
      double Wcm[npfit],rhoPPCT14[npfit],rhoPBARPCT14[npfit],rhoPPCTEQ6L[npfit],rhoPBARPCTEQ6L[npfit],rhoPPMMHT[npfit],rhoPBARPMMHT[npfit];       
            
      for(int j=0;j<npfit;j++){
          fscanf(aq1,"%lg %lg",&Wcm[j],&rhoPPCT14[j]);
          fscanf(aq2,"%lg %lg",&Wcm[j],&rhoPBARPCT14[j]);
      }
      
      fclose(aq1);
      fclose(aq2);
     /********************************************************************************/
    
     /********************************************************************************/
     //cria ponteiro do tipo TGraph para desenhar curvas de Rho pp      

      TGraph *grRhoPP_CT14 = new TGraph(npfit,Wcm,rhoPPCT14);
      grRhoPP_CT14->SetLineColor(2);
      grRhoPP_CT14->SetLineWidth(1);
      grRhoPP_CT14->SetLineStyle(1);
      grRhoPP_CT14->Draw("c");

      /********************************************************************************/
  
     /********************************************************************************/
     //cria ponteiro do tipo TGraph para desenhar curvas de Rho pbarp  
        
      TGraph *grRhoPBARP_CT14 = new TGraph(npfit,Wcm,rhoPBARPCT14);
      grRhoPBARP_CT14->SetLineColor(2);
      grRhoPBARP_CT14->SetLineWidth(1);
      grRhoPBARP_CT14->SetLineStyle(2);
      grRhoPBARP_CT14->Draw("c");

     /********************************************************************************/
     //cria legenda do gráfico
     
      TLegend *leg2 = new TLegend(0.12,0.75,0.45,0.88);
      leg2->AddEntry(gr00,"#bar{p}p accelerator data ","p");
      leg2->AddEntry(gr11,"pp accelerator data","p");
//      leg2->AddEntry(grRhoPP_CT14,"CT14 ","l");
//      leg2->SetFillColor(0);
      leg2->SetShadowColor(0);
      leg2->SetBorderSize(0);
      leg2->SetTextSize(0.035);
      leg2->Draw();
     
        /********************************************************************************/
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
     
      /********************************************************************************/
      //cria arquivo .eps com o gráfico de saída
      
      canv2->SaveAs("Rh_DGM19cCT14cut7vFF_teste.eps");    
      /********************************************************************************/
   
}
