void plotCDHdE(const char* filename="evanaIMpisigma_v156.root")
{
  TFile *file = new TFile(filename,"READ");

  TCanvas *cdEall = new TCanvas("cdEall","cdEall");
  cdEall->cd()->SetLogy();
  TH2D* CDHdE = (TH2D*)file->Get("CDHdE");
  CDHdE->ProjectionY()->Draw("HE");
  TH2D* CDHdE_wt = (TH2D*)file->Get("CDHdE_wt");
  TH1D *CDHdE_wt_py = (TH1D*)CDHdE_wt->ProjectionY();
  CDHdE_wt_py->SetLineColor(2);
  CDHdE_wt_py->Draw("HEsame");

  TCanvas *cdEeach = new TCanvas("cdEeach","cdEeach",1400,1000);
  cdEeach->cd();
  //cdEeach->Divide(6,6);
  double max = 0;
  for(int iseg=0;iseg<36;iseg++){
    //cdEeach->cd(iseg+1)->SetLogy();
    //cdEeach->cd(iseg+1);
    TH1D *his = (TH1D*)CDHdE->ProjectionY(Form("h%d",iseg+1),iseg+1,iseg+1);
    his->GetXaxis()->SetRangeUser(0,5);
    //his->Draw("HE");
    TH1D *his2 = (TH1D*)CDHdE_wt->ProjectionY(Form("h_wt%d",iseg+1),iseg+1,iseg+1);
    //if(iseg==0) max = his2->GetMaximum()+10000;
    his2->SetMaximum(200000);
    his2->SetLineColor(2);
    his2->GetXaxis()->SetRangeUser(0,10);
    //his2->Draw("HEsame");
    his2->SetLineColor((iseg+1)%10);
    his2->SetLineStyle((iseg+1)/10);
    if(iseg==0)his2->Draw("HE");
    else his2->Draw("HEsame");

    his2->GetXaxis()->SetRangeUser(0,10);
    int maxbin = his2->GetMaximumBin();
    std::cout << iseg+1 << "   "  << his2->GetBinCenter(maxbin) << std::endl;
  }
 
  
  TCanvas *cCDHtime = new TCanvas("cCDHtime","cCDHtime",1400,1000);
  cCDHtime->cd();
  TH2D* CDHtime = (TH2D*)file->Get("CDHtime");
  for(int iseg=0;iseg<36;iseg++){
    TH1D *his = (TH1D*)CDHtime->ProjectionY(Form("h%d",iseg+1),iseg+1,iseg+1);
    his->GetXaxis()->SetRangeUser(0,100);
    his->SetLineColor(iseg%10);
    his->SetLineStyle(iseg/10);
    if(iseg==0)his->Draw("HE");
    else his->Draw("HEsame");
  }

}
