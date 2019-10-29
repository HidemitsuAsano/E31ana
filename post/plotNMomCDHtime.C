void plotNMomCDHtime(const char* filename="evanaIMpisigma_v163.root")
{
  TFile *file = new TFile(filename,"READ");


  for(int iseg=0;iseg<36;iseg++){
    TH2D *his2 = (TH2D*)file->Get(Form("NMomCDHtime%d",iseg+1));
    TH1D *his = his2->ProjectionY();
    his->SetLineColor(iseg%10);
    his->SetLineStyle(iseg/10);
    his->Draw("same");
  }






}
