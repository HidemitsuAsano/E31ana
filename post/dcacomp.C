void dcacomp(){
  
  TCanvas *c1 = new TCanvas("c1","c1");
  TFile *_file0 = TFile::Open("evanaIMpisigma_v43.root");
  DCA_pippim->Rebin(5);
  DCA_pippim->Draw();
  TH1F* real = (TH1F*)DCA_pippim->Clone();
  TFile *_file1 = TFile::Open("../simpost/simIMpisigma_nSppim_DoraAir_v15.root");
  DCA_pippim->SetLineColor(2);
  DCA_pippim->Rebin(5);
  DCA_pippim->Draw("same");
  TH1F* sp = (TH1F*)DCA_pippim->Clone();
 
  TFile *_file2 = TFile::Open("../simpost/simIMpisigma_nSmpip_DoraAir_v15.root");
  DCA_pippim->SetLineColor(3);
  DCA_pippim->Rebin(5);
  DCA_pippim->Draw("same");
  TH1F* sm = (TH1F*)DCA_pippim->Clone();
  real->GetXaxis()->SetTitle("DCA (#pi^{+}-#pi^{-}) [cm]");
  real->GetXaxis()->CenterTitle();
  real->Draw();
  sp->Scale(real->GetMaximum()/sp->GetMaximum());
  sp->Draw("same");
  sm->Scale(real->GetMaximum()/sm->GetMaximum());
  sm->Draw("same");
 
  TCanvas *c2 = new TCanvas("c2","c2");
  _file0->cd();
  TH1F* real2 = (TH1F*)DCA_pippim_SigmaPM->Clone();
  real2->Rebin(5);
 
  real2->GetXaxis()->SetTitle("DCA (#pi^{+}-#pi^{-}) [cm]");
  real2->GetXaxis()->CenterTitle();
  real2->Draw();
  _file1->cd();
  DCA_pippim_SigmaPM->SetLineColor(2);
  TH1F* sp2 = (TH1F*)DCA_pippim_SigmaPM->Clone();
  sp2->Rebin(5);
  sp2->Scale(real2->GetMaximum()/sp2->GetMaximum());
  sp2->Draw("same");
  _file2->cd();
  DCA_pippim_SigmaPM->SetLineColor(3);
  TH1F* sm2 = (TH1F*)DCA_pippim_SigmaPM->Clone();
  sm2->Rebin(5);
  sm2->Scale(real2->GetMaximum()/sm2->GetMaximum());
  sm2->Draw("same");

  TCanvas *c3 = new TCanvas("c3","c3");
  _file0->cd();
  TH1F* rpip = (TH1F*)DCA_pip->Clone();
  TH1F* rpim = (TH1F*)DCA_pim->Clone();
  rpip->Draw();
  rpim->SetLineColor(4);
  rpim->Draw("same");
  _file1->cd();
  //TH1F* rpipSp 

}
