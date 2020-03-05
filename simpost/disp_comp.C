void disp_comp(const char *filename="comp_fakedata_out.root")
{
  TFile *f = TFile::Open(filename,"READ");
  f->cd();

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  TH1D* q_woK0_woSid_won_rdata = (TH1D*)f->Get("q_woK0_woSid_won_rdata");
  TH1D* q_woK0_woSid_won_mc = (TH1D*)f->Get("q_woK0_woSid_won_mc");
  TH1D* q_wK0_woSid_won_rdata = (TH1D*)f->Get("q_wK0_woSid_won_rdata");
  TH1D* q_wK0_woSid_won_mc = (TH1D*)f->Get("q_wK0_woSid_won_mc");
  TH1D* q_woK0_woSid_won_ratio = f->Get("q_woK0_woSid_won_ratio");
  TH1D* q_wK0_woSid_won_ratio = f->Get("q_wK0_woSid_won_ratio");
  
  TCanvas *cq = new TCanvas("cq","cq",1000,1000);
  cq->Divide(2,2);
  cq->cd(1);
  q_woK0_woSid_won_rdata->Draw("HE");
  q_woK0_woSid_won_mc->Draw("HEsame");
  cq->cd(2);
  q_wK0_woSid_won_rdata->Draw("HE");
  q_wK0_woSid_won_mc->Draw("HEsame");
  cq->cd(3);
  q_woK0_woSid_won_ratio->Draw("HE");
  cq->cd(4);
  q_wK0_woSid_won_ratio->Draw("HE");
  cq->Update();
  cq->Modified();
  
  TH1D* MMnmiss_woK0_woSid_won_rdata = (TH1D*)f->Get("MMnmiss_woK0_woSid_won_rdata");
  TH1D* MMnmiss_woK0_woSid_won_mc = (TH1D*)f->Get("MMnmiss_woK0_woSid_won_mc");
  TH1D* MMnmiss_wK0_woSid_won_rdata = (TH1D*)f->Get("MMnmiss_wK0_woSid_won_rdata");
  TH1D* MMnmiss_wK0_woSid_won_mc = (TH1D*)f->Get("MMnmiss_wK0_woSid_won_mc");
  TH1D* MMnmiss_woK0_woSid_won_ratio = f->Get("MMnmiss_woK0_woSid_won_ratio");
  TH1D* MMnmiss_wK0_woSid_won_ratio = f->Get("MMnmiss_wK0_woSid_won_ratio");
  
  TCanvas *cMMnmiss = new TCanvas("cMMnmiss","cMMnmiss",1000,1000);
  cMMnmiss->Divide(2,2);
  cMMnmiss->cd(1);
  MMnmiss_woK0_woSid_won_rdata->Draw("HE");
  MMnmiss_woK0_woSid_won_mc->Draw("HEsame");
  cMMnmiss->cd(2);
  MMnmiss_wK0_woSid_won_rdata->Draw("HE");
  MMnmiss_wK0_woSid_won_mc->Draw("HEsame");
  cMMnmiss->cd(3);
  MMnmiss_woK0_woSid_won_ratio->Draw("HE");
  cMMnmiss->cd(4);
  MMnmiss_wK0_woSid_won_ratio->Draw("HE");

  TH1D* IMnpip_woK0_woSid_won_rdata = (TH1D*)f->Get("IMnpip_woK0_woSid_won_rdata");
  TH1D* IMnpip_woK0_woSid_won_mc = (TH1D*)f->Get("IMnpip_woK0_woSid_won_mc");
  TH1D* IMnpip_wK0_woSid_won_rdata = (TH1D*)f->Get("IMnpip_wK0_woSid_won_rdata");
  TH1D* IMnpip_wK0_woSid_won_mc = (TH1D*)f->Get("IMnpip_wK0_woSid_won_mc");
  TH1D* IMnpip_woK0_woSid_won_ratio = f->Get("IMnpip_woK0_woSid_won_ratio");
  TH1D* IMnpip_wK0_woSid_won_ratio = f->Get("IMnpip_wK0_woSid_won_ratio");
  
  TCanvas *cIMnpip = new TCanvas("cIMnpip","cIMnpip",1000,1000);
  cIMnpip->Divide(2,2);
  cIMnpip->cd(1);
  IMnpip_woK0_woSid_won_rdata->Draw("HE");
  IMnpip_woK0_woSid_won_mc->Draw("HEsame");
  cIMnpip->cd(2);
  IMnpip_wK0_woSid_won_rdata->Draw("HE");
  IMnpip_wK0_woSid_won_mc->Draw("HEsame");
  cIMnpip->cd(3);
  IMnpip_woK0_woSid_won_ratio->Draw("HE");
  cIMnpip->cd(4);
  IMnpip_wK0_woSid_won_ratio->Draw("HE");

  TH1D* IMnpim_woK0_woSid_won_rdata = (TH1D*)f->Get("IMnpim_woK0_woSid_won_rdata");
  TH1D* IMnpim_woK0_woSid_won_mc = (TH1D*)f->Get("IMnpim_woK0_woSid_won_mc");
  TH1D* IMnpim_wK0_woSid_won_rdata = (TH1D*)f->Get("IMnpim_wK0_woSid_won_rdata");
  TH1D* IMnpim_wK0_woSid_won_mc = (TH1D*)f->Get("IMnpim_wK0_woSid_won_mc");
  TH1D* IMnpim_woK0_woSid_won_ratio = f->Get("IMnpim_woK0_woSid_won_ratio");
  TH1D* IMnpim_wK0_woSid_won_ratio = f->Get("IMnpim_wK0_woSid_won_ratio");
  
  TCanvas *cIMnpim = new TCanvas("cIMnpim","cIMnpim",1000,1000);
  cIMnpim->Divide(2,2);
  cIMnpim->cd(1);
  IMnpim_woK0_woSid_won_rdata->Draw("HE");
  IMnpim_woK0_woSid_won_mc->Draw("HEsame");
  cIMnpim->cd(2);
  IMnpim_wK0_woSid_won_rdata->Draw("HE");
  IMnpim_wK0_woSid_won_mc->Draw("HEsame");
  cIMnpim->cd(3);
  IMnpim_woK0_woSid_won_ratio->Draw("HE");
  cIMnpim->cd(4);
  IMnpim_wK0_woSid_won_ratio->Draw("HE");





















}
