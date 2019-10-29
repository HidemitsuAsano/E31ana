void plot_dEnmom()
{
  TFile *_file0 = TFile::Open("evanaIMpisigma_npippim_v176_outncutK015.root");
  TCanvas *c1 = new TCanvas("c1","c1");
  TH2F *dE_nmom_fid_beta_wK0 = (TH2F*)_file0->Get("dE_nmom_fid_beta_wK0");
  dE_nmom_fid_beta_wK0->SetTitle("dE_nmom_fid_beta_wK0 real data");
  dE_nmom_fid_beta_wK0->Draw("colz");


  TCanvas *c2 = new TCanvas("c2","c2");
  double nmom100MeVbin = dE_nmom_fid_beta_wK0->GetXaxis()->FindBin(0.1);
  double nmom200MeVbin = dE_nmom_fid_beta_wK0->GetXaxis()->FindBin(0.2);
  double nmom300MeVbin = dE_nmom_fid_beta_wK0->GetXaxis()->FindBin(0.3);
  double nmom400MeVbin = dE_nmom_fid_beta_wK0->GetXaxis()->FindBin(0.4);
  double nmom500MeVbin = dE_nmom_fid_beta_wK0->GetXaxis()->FindBin(0.5);
  TH1D *dEnmom0MeV = (TH1D*)dE_nmom_fid_beta_wK0->ProjectionY("dE0MeV",2,nmom100MeVbin);
  TH1D *dEnmom100MeV = (TH1D*)dE_nmom_fid_beta_wK0->ProjectionY("dE100MeV",nmom100MeVbin,nmom200MeVbin);
  dEnmom100MeV->SetLineColor(2);
  dEnmom100MeV->Draw("HE");
  dEnmom0MeV->Draw("HEsame");
  TH1D *dEnmom200MeV = (TH1D*)dE_nmom_fid_beta_wK0->ProjectionY("dE200MeV",nmom200MeVbin,nmom300MeVbin);
  dEnmom200MeV->SetLineColor(3);
  dEnmom200MeV->Draw("HEsame");

  TH1D *dEnmom300MeV = (TH1D*)dE_nmom_fid_beta_wK0->ProjectionY("dE300MeV",nmom300MeVbin,nmom400MeVbin);
  dEnmom300MeV->SetLineColor(4);
  dEnmom300MeV->Draw("HEsame");

  TH1D *dEnmom400MeV = (TH1D*)dE_nmom_fid_beta_wK0->ProjectionY("dE400MeV",nmom400MeVbin,nmom500MeVbin);
  dEnmom400MeV->SetLineColor(5);
  dEnmom400MeV->Draw("HEsame");
  


  TFile *_file1 = TFile::Open("../simpost/simIMpisigma_K0n_ns_pippimn_v16_outncutK015.root");
  _file1->cd();
  
  TCanvas *c3 = new TCanvas("c3","c3");
  TH2F* dE_nmom_fid_beta_wK0_mc = (TH2F*)_file1->Get("dE_nmom_fid_beta_wK0");
  dE_nmom_fid_beta_wK0_mc->Draw("colz");


  TCanvas *c4 = new TCanvas("c4","c4");
  double nmom100MeVbin = dE_nmom_fid_beta_wK0_mc->GetXaxis()->FindBin(0.1);
  double nmom200MeVbin = dE_nmom_fid_beta_wK0_mc->GetXaxis()->FindBin(0.2);
  double nmom300MeVbin = dE_nmom_fid_beta_wK0_mc->GetXaxis()->FindBin(0.3);
  double nmom400MeVbin = dE_nmom_fid_beta_wK0_mc->GetXaxis()->FindBin(0.4);
  double nmom500MeVbin = dE_nmom_fid_beta_wK0_mc->GetXaxis()->FindBin(0.5);
  TH1D *dEnmom0MeV_mc = (TH1D*)dE_nmom_fid_beta_wK0_mc->ProjectionY("dE0MeV_mc",2,nmom100MeVbin);
  dEnmom0MeV_mc->Draw("HE");
  TH1D *dEnmom100MeV_mc = (TH1D*)dE_nmom_fid_beta_wK0_mc->ProjectionY("dE100MeV_mc",nmom100MeVbin,nmom200MeVbin);
  dEnmom100MeV_mc->SetLineColor(2);
  dEnmom100MeV_mc->Draw("HEsame");
  TH1D *dEnmom200MeV_mc = (TH1D*)dE_nmom_fid_beta_wK0_mc->ProjectionY("dE200MeV_mc",nmom200MeVbin,nmom300MeVbin);
  dEnmom200MeV_mc->SetLineColor(3);
  dEnmom200MeV_mc->Draw("HEsame");
  TH1D *dEnmom300MeV_mc = (TH1D*)dE_nmom_fid_beta_wK0_mc->ProjectionY("dE300MeV_mc",nmom300MeVbin,nmom400MeVbin);
  dEnmom300MeV_mc->SetLineColor(4);
  dEnmom300MeV_mc->Draw("HEsame");
  TH1D *dEnmom400MeV_mc = (TH1D*)dE_nmom_fid_beta_wK0_mc->ProjectionY("dE400MeV_mc",nmom400MeVbin,nmom500MeVbin);
  dEnmom400MeV_mc->SetLineColor(5);
  dEnmom400MeV_mc->Draw("HEsame");

  TCanvas *c5 = new TCanvas("c5","c5");
  TH1D* dEnmom100MeV_clone = (TH1D*)dEnmom100MeV->Clone("dEnmom100MeV_clone");
  dEnmom100MeV->SetTitle("1 GeV/c < nmon < 2 GeV/c") ;
  dEnmom100MeV->SetLineColor(1);
  dEnmom100MeV->Draw("HE");
  TH1D* dEnmom100MeV_mc_clone = (TH1D*)dEnmom100MeV_mc->Clone("dEnmom100MeV_mc_clone");
  dEnmom100MeV_mc_clone->SetLineColor(2);
  dEnmom100MeV_mc_clone->Draw("hesame");

  TCanvas *c6 = new TCanvas("c6","c6");
  TH1D* dEnmom200MeV_clone = (TH1D*)dEnmom200MeV->Clone("dEnmom200MeV_clone");
  dEnmom200MeV_clone->SetTitle("2 GeV/c < nmon < 3 GeV/c") ;
  dEnmom200MeV_clone->SetLineColor(1);
  dEnmom200MeV_clone->Draw("HE");
  TH1D* dEnmom200MeV_mc_clone = (TH1D*)dEnmom200MeV_mc->Clone("dEnmom200MeV_mc_clone");
  dEnmom200MeV_mc_clone->SetLineColor(2);
  dEnmom200MeV_mc_clone->Scale(1.9);
  dEnmom200MeV_mc_clone->Draw("hesame");

  TCanvas *c7 = new TCanvas("c7","c7");
  TH1D* dEnmom300MeV_clone = (TH1D*)dEnmom300MeV->Clone("dEnmom300MeV_clone");
  dEnmom300MeV_clone->SetTitle("3 GeV/c < nmon < 4 GeV/c") ;
  dEnmom300MeV_clone->SetLineColor(1);
  dEnmom300MeV_clone->Draw("HE");
  TH1D* dEnmom300MeV_mc_clone = (TH1D*)dEnmom300MeV_mc->Clone("dEnmom300MeV_mc_clone");
  dEnmom300MeV_mc_clone->SetLineColor(2);
  dEnmom300MeV_mc_clone->Scale(1.9);
  dEnmom300MeV_mc_clone->Draw("hesame");

  TCanvas *c8 = new TCanvas("c8","c8");
  TH1D* dEnmom400MeV_clone = (TH1D*)dEnmom400MeV->Clone("dEnmom400MeV_clone");
  dEnmom400MeV_clone->SetTitle("4 GeV/c < nmon < 5 GeV/c") ;
  dEnmom400MeV_clone->SetLineColor(1);
  dEnmom400MeV_clone->Draw("HE");
  TH1D* dEnmom400MeV_mc_clone = (TH1D*)dEnmom400MeV_mc->Clone("dEnmom400MeV_mc_clone");
  dEnmom400MeV_mc_clone->SetLineColor(2);
  dEnmom400MeV_mc_clone->Scale(1.9);
  dEnmom400MeV_mc_clone->Draw("hesame");


}
