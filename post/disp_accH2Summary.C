void disp_accH2Summary()
{
  gStyle->SetOptStat(0);
  
  TCanvas *c1 = new TCanvas("c1","c1");
  TFile *_file0 = TFile::Open("accH2dE2.root");
  TH1F* accCosSp2 = (TH1F*)_file0->Get("accCosSp2");
  accCosSp2->SetLineColor(1);
  accCosSp2->SetTitle("#Sigma^{+}#pi^{-} acc. #times eff.");
  accCosSp2->Draw("");

  TFile *_file1 = TFile::Open("accH2dE4.root");
  TH1F* accCosSp2dE4 = (TH1F*)_file1->Get("accCosSp2");
  accCosSp2dE4->SetLineColor(2);
  accCosSp2dE4->Draw("same");
  
  TFile *_file2 = TFile::Open("accH2dE6.root");
  TH1F* accCosSp2dE6 = (TH1F*)_file2->Get("accCosSp2");
  accCosSp2dE6->SetLineColor(3);
  accCosSp2dE6->Draw("same");

  TLegend *leg = new TLegend(0.3, 0.6, 0.6, 0.9);
  leg->AddEntry(accCosSp2, "dE < 2 MeVee", "LP");
  leg->AddEntry(accCosSp2dE4, "dE < 4 MeVee", "LP");
  leg->AddEntry(accCosSp2dE6, "dE < 6 MeVee", "LP");
  leg->Draw();

  TCanvas *c2 = new TCanvas("c2","c2");
  TH1F* accCosSm2 = (TH1F*)_file0->Get("accCosSm2");
  accCosSm2->SetTitle("#Sigma^{-}#pi^{+} acc. #times eff.");
  accCosSm2->SetLineColor(1);
  accCosSm2->Draw("");

  TH1F* accCosSm2dE4 = (TH1F*)_file1->Get("accCosSm2");
  accCosSm2dE4->SetLineColor(2);
  accCosSm2dE4->Draw("same");
  
  TH1F* accCosSm2dE6 = (TH1F*)_file2->Get("accCosSm2");
  accCosSm2dE6->SetLineColor(3);
  accCosSm2dE6->Draw("same");

  TLegend *leg2 = new TLegend(0.3, 0.6, 0.6, 0.9);
  leg2->AddEntry(accCosSm2, "dE < 2 MeVee", "LP");
  leg2->AddEntry(accCosSm2dE4, "dE < 4 MeVee", "LP");
  leg2->AddEntry(accCosSm2dE6, "dE < 6 MeVee", "LP");
  leg2->Draw();

}
