void plot_MML(){
  gStyle->SetOptStat("i");

  TFile *_file0 = TFile::Open("evanaIMLambda_v1.root");
  TCanvas *c1 = new TCanvas("c1","c1",1000,1000);
  c1->Divide(2,2,0,0);
  c1->cd(3);
  MMass_IMppim->SetXTitle("IM(p#pi^{-}) [GeV/c^{2}]");
  MMass_IMppim->GetXaxis()->CenterTitle();
  MMass_IMppim->SetYTitle("MM(d(K^{-},p#pi^{-})\"X\") [GeV/c^{2}]");
  MMass_IMppim->GetYaxis()->CenterTitle();
  MMass_IMppim->GetYaxis()->SetTitleOffset(1.2);
  MMass_IMppim->Draw("colz");
  MMass_IMppim->SetTitle("");

  c1->SetLogz();
  TH1D* IMppim = (TH1D*)MMass_IMppim->ProjectionX("IMppim");
  TH1D* MMppim = (TH1D*)MMass_IMppim->ProjectionY("MMppim");

  c1->cd(1);
  IMppim->GetXaxis()->SetLabelSize(0);
  IMppim->SetMarkerStyle(20);
  IMppim->Draw("E");
  
  c1->cd(4);
  TGraphErrors *gr_MMppim = new TGraphErrors();  
  for(int ibin=0;ibin<MMppim->GetNbinsX();ibin++){
    double cont = MMppim->GetBinContent(ibin);
    double err = MMppim->GetBinError(ibin);
    double bincenter = MMppim->GetBinCenter(ibin);
    gr_MMppim->SetPoint(ibin,cont, bincenter);
    gr_MMppim->SetPointError(ibin,err,0);
  }
  gr_MMppim->SetMarkerStyle(20);
  gr_MMppim->GetYaxis()->SetRangeUser(0,2);
  gr_MMppim->GetYaxis()->SetLabelSize(0);
  gr_MMppim->Draw("AP");

  TCanvas *c2 = new TCanvas("c2","c2");
  c2->cd();
  c2->SetLogy();
  q_MMass->RebinX(5);
  q_MMass->SetXTitle("MM(d(K^{-},p#pi^{-})\"X\") [GeV/c^{2}]");
  q_MMass->GetXaxis()->CenterTitle();
  q_MMass->ProjectionX()->Draw("HE");

}
