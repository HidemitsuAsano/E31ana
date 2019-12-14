void plot_qvsIMnpipi()
{
  TFile *file0 = TFile::Open("evanaIMpisigma_npippim_v195_out.root","READ");
  file0->cd();
 
  TH2F* q_IMnpipi_wSid_n_Sp = (TH2F*)file0->Get("q_IMnpipi_wSid_n_Sp");
  TH2F* q_IMnpipi_wSid_n_Sm = (TH2F*)file0->Get("q_IMnpipi_wSid_n_Sm");
  TH2F* q_IMnpipi_wK0_wSid_n_Sp = (TH2F*)file0->Get("q_IMnpipi_wK0_wSid_n_Sp");
  TH2F* q_IMnpipi_wK0_wSid_n_Sm = (TH2F*)file0->Get("q_IMnpipi_wK0_wSid_n_Sm");


  TCanvas *c1 = new TCanvas("c1","c1");
  c1->cd();
  q_IMnpipi_wSid_n_Sp->Draw("colz");

  TCanvas *c2 = new TCanvas("c2","c2");
  c2->cd();
  q_IMnpipi_wSid_n_Sm->Draw("colz");

  TCanvas *c3 = new TCanvas("c3","c3");
  c3->cd();
  q_IMnpipi_wK0_wSid_n_Sp->Draw("colz");

  TCanvas *c4 = new TCanvas("c4","c4");
  c4->cd();
  q_IMnpipi_wK0_wSid_n_Sm->Draw("colz");

  TCanvas *c5 = new TCanvas("c5","c5");
  c5->cd();
  TH1D* IMnpipi_Sp = q_IMnpipi_wSid_n_Sp->ProjectionX("IMnpipi_Sp");
  IMnpipi_Sp->SetLineColor(2);
  IMnpipi_Sp->Draw("HE");
  TH1D* IMnpipi_Sm = q_IMnpipi_wSid_n_Sm->ProjectionX("IMnpipi_Sm");
  IMnpipi_Sm->SetLineColor(3);
  IMnpipi_Sm->Draw("HEsame");





}
