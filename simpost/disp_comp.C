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
  TH1D* q_woK0_woSid_won_ratio = (TH1D*)f->Get("q_woK0_woSid_won_ratio");
  TH1D* q_wK0_woSid_won_ratio = (TH1D*)f->Get("q_wK0_woSid_won_ratio");
  
  TCanvas *cq = new TCanvas("cq","cq",1000,1000);
  cq->Divide(2,2);
  cq->cd(1);
  q_woK0_woSid_won_rdata->Draw("HE");
  q_woK0_woSid_won_mc->Draw("HEsame");
  cq->cd(2);
  q_wK0_woSid_won_rdata->Draw("HE");
  q_wK0_woSid_won_mc->Draw("HEsame");
  cq->cd(3);
  q_woK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,2);
  q_woK0_woSid_won_ratio->Draw("HE");
  cq->cd(4);
  q_wK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,2);
  q_wK0_woSid_won_ratio->Draw("HE");
  cq->Update();
  cq->Modified();
  
  TH1D* MMnmiss_woK0_woSid_won_rdata = (TH1D*)f->Get("MMnmiss_woK0_woSid_won_rdata");
  TH1D* MMnmiss_woK0_woSid_won_mc = (TH1D*)f->Get("MMnmiss_woK0_woSid_won_mc");
  TH1D* MMnmiss_wK0_woSid_won_rdata = (TH1D*)f->Get("MMnmiss_wK0_woSid_won_rdata");
  TH1D* MMnmiss_wK0_woSid_won_mc = (TH1D*)f->Get("MMnmiss_wK0_woSid_won_mc");
  TH1D* MMnmiss_woK0_woSid_won_ratio = (TH1D*)f->Get("MMnmiss_woK0_woSid_won_ratio");
  TH1D* MMnmiss_wK0_woSid_won_ratio = (TH1D*)f->Get("MMnmiss_wK0_woSid_won_ratio");
  
  TCanvas *cMMnmiss = new TCanvas("cMMnmiss","cMMnmiss",1000,1000);
  cMMnmiss->Divide(2,2);
  cMMnmiss->cd(1);
  MMnmiss_woK0_woSid_won_rdata->Draw("HE");
  MMnmiss_woK0_woSid_won_mc->Draw("HEsame");
  cMMnmiss->cd(2);
  MMnmiss_wK0_woSid_won_rdata->Draw("HE");
  MMnmiss_wK0_woSid_won_mc->Draw("HEsame");
  cMMnmiss->cd(3);
  MMnmiss_woK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,2);
  MMnmiss_woK0_woSid_won_ratio->Draw("HE");
  cMMnmiss->cd(4);
  MMnmiss_wK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,2);
  MMnmiss_wK0_woSid_won_ratio->Draw("HE");

  TH1D* IMnpip_woK0_woSid_won_rdata = (TH1D*)f->Get("IMnpip_woK0_woSid_won_rdata");
  TH1D* IMnpip_woK0_woSid_won_mc = (TH1D*)f->Get("IMnpip_woK0_woSid_won_mc");
  TH1D* IMnpip_wK0_woSid_won_rdata = (TH1D*)f->Get("IMnpip_wK0_woSid_won_rdata");
  TH1D* IMnpip_wK0_woSid_won_mc = (TH1D*)f->Get("IMnpip_wK0_woSid_won_mc");
  TH1D* IMnpip_woK0_woSid_won_ratio = (TH1D*)f->Get("IMnpip_woK0_woSid_won_ratio");
  TH1D* IMnpip_wK0_woSid_won_ratio = (TH1D*)f->Get("IMnpip_wK0_woSid_won_ratio");
  
  TCanvas *cIMnpip = new TCanvas("cIMnpip","cIMnpip",1000,1000);
  cIMnpip->Divide(2,2);
  cIMnpip->cd(1);
  IMnpip_woK0_woSid_won_rdata->Draw("HE");
  IMnpip_woK0_woSid_won_mc->Draw("HEsame");
  cIMnpip->cd(2);
  IMnpip_wK0_woSid_won_rdata->Draw("HE");
  IMnpip_wK0_woSid_won_mc->Draw("HEsame");
  cIMnpip->cd(3);
  IMnpip_woK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,2);
  IMnpip_woK0_woSid_won_ratio->Draw("HE");
  cIMnpip->cd(4);
  IMnpip_wK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,2);
  IMnpip_wK0_woSid_won_ratio->Draw("HE");

  TH1D* IMnpim_woK0_woSid_won_rdata = (TH1D*)f->Get("IMnpim_woK0_woSid_won_rdata");
  TH1D* IMnpim_woK0_woSid_won_mc = (TH1D*)f->Get("IMnpim_woK0_woSid_won_mc");
  TH1D* IMnpim_wK0_woSid_won_rdata = (TH1D*)f->Get("IMnpim_wK0_woSid_won_rdata");
  TH1D* IMnpim_wK0_woSid_won_mc = (TH1D*)f->Get("IMnpim_wK0_woSid_won_mc");
  TH1D* IMnpim_woK0_woSid_won_ratio = (TH1D*)f->Get("IMnpim_woK0_woSid_won_ratio");
  TH1D* IMnpim_wK0_woSid_won_ratio = (TH1D*)f->Get("IMnpim_wK0_woSid_won_ratio");
  
  TCanvas *cIMnpim = new TCanvas("cIMnpim","cIMnpim",1000,1000);
  cIMnpim->Divide(2,2);
  cIMnpim->cd(1);
  IMnpim_woK0_woSid_won_rdata->Draw("HE");
  IMnpim_woK0_woSid_won_mc->Draw("HEsame");
  cIMnpim->cd(2);
  IMnpim_wK0_woSid_won_rdata->Draw("HE");
  IMnpim_wK0_woSid_won_mc->Draw("HEsame");
  cIMnpim->cd(3);
  IMnpim_woK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,2);
  IMnpim_woK0_woSid_won_ratio->Draw("HE");
  cIMnpim->cd(4);
  IMnpim_wK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,2);
  IMnpim_wK0_woSid_won_ratio->Draw("HE");


  TH1D* Momnpip_woK0_woSid_won_rdata = (TH1D*)f->Get("Momnpip_woK0_woSid_won_rdata");
  TH1D* Momnpip_woK0_woSid_won_mc = (TH1D*)f->Get("Momnpip_woK0_woSid_won_mc");
  TH1D* Momnpip_wK0_woSid_won_rdata = (TH1D*)f->Get("Momnpip_wK0_woSid_won_rdata");
  TH1D* Momnpip_wK0_woSid_won_mc = (TH1D*)f->Get("Momnpip_wK0_woSid_won_mc");
  TH1D* Momnpip_woK0_woSid_won_ratio = (TH1D*)f->Get("Momnpip_woK0_woSid_won_ratio");
  TH1D* Momnpip_wK0_woSid_won_ratio = (TH1D*)f->Get("Momnpip_wK0_woSid_won_ratio");
  
  TCanvas *cMomnpip = new TCanvas("cMomnpip","cMomnpip",1000,1000);
  cMomnpip->Divide(2,2);
  cMomnpip->cd(1);
  Momnpip_woK0_woSid_won_rdata->Draw("HE");
  Momnpip_woK0_woSid_won_mc->Draw("HEsame");
  cMomnpip->cd(2);
  Momnpip_wK0_woSid_won_rdata->Draw("HE");
  Momnpip_wK0_woSid_won_mc->Draw("HEsame");
  cMomnpip->cd(3);
  Momnpip_woK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,2);
  Momnpip_woK0_woSid_won_ratio->Draw("HE");
  cMomnpip->cd(4);
  Momnpip_wK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,2);
  Momnpip_wK0_woSid_won_ratio->Draw("HE");

  TH1D* Momnpim_woK0_woSid_won_rdata = (TH1D*)f->Get("Momnpim_woK0_woSid_won_rdata");
  TH1D* Momnpim_woK0_woSid_won_mc = (TH1D*)f->Get("Momnpim_woK0_woSid_won_mc");
  TH1D* Momnpim_wK0_woSid_won_rdata = (TH1D*)f->Get("Momnpim_wK0_woSid_won_rdata");
  TH1D* Momnpim_wK0_woSid_won_mc = (TH1D*)f->Get("Momnpim_wK0_woSid_won_mc");
  TH1D* Momnpim_woK0_woSid_won_ratio = (TH1D*)f->Get("Momnpim_woK0_woSid_won_ratio");
  TH1D* Momnpim_wK0_woSid_won_ratio = (TH1D*)f->Get("Momnpim_wK0_woSid_won_ratio");
  
  TCanvas *cMomnpim = new TCanvas("cMomnpim","cMomnpim",1000,1000);
  cMomnpim->Divide(2,2);
  cMomnpim->cd(1);
  Momnpim_woK0_woSid_won_rdata->Draw("HE");
  Momnpim_woK0_woSid_won_mc->Draw("HEsame");
  cMomnpim->cd(2);
  Momnpim_wK0_woSid_won_rdata->Draw("HE");
  Momnpim_wK0_woSid_won_mc->Draw("HEsame");
  cMomnpim->cd(3);
  Momnpim_woK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,2);
  Momnpim_woK0_woSid_won_ratio->Draw("HE");
  cMomnpim->cd(4);
  Momnpim_wK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,2);
  Momnpim_wK0_woSid_won_ratio->Draw("HE");


  TH1D* IMnpipi_woK0_woSid_won_rdata = (TH1D*)f->Get("IMnpipi_woK0_woSid_won_rdata");
  TH1D* IMnpipi_woK0_woSid_won_mc = (TH1D*)f->Get("IMnpipi_woK0_woSid_won_mc");
  TH1D* IMnpipi_wK0_woSid_won_rdata = (TH1D*)f->Get("IMnpipi_wK0_woSid_won_rdata");
  TH1D* IMnpipi_wK0_woSid_won_mc = (TH1D*)f->Get("IMnpipi_wK0_woSid_won_mc");
  TH1D* IMnpipi_woK0_woSid_won_ratio = (TH1D*)f->Get("IMnpipi_woK0_woSid_won_ratio");
  TH1D* IMnpipi_wK0_woSid_won_ratio = (TH1D*)f->Get("IMnpipi_wK0_woSid_won_ratio");

  TCanvas *cIMnpipi = new TCanvas("cIMnpipi","cIMnpipi",1000,1000);
  cIMnpipi->Divide(2,2);
  cIMnpipi->cd(1);
  IMnpipi_woK0_woSid_won_rdata->Draw("HE");
  IMnpipi_woK0_woSid_won_mc->Draw("HEsame");
  cIMnpipi->cd(2);
  IMnpipi_wK0_woSid_won_rdata->Draw("HE");
  IMnpipi_wK0_woSid_won_mc->Draw("HEsame");
  cIMnpipi->cd(3);
  IMnpipi_woK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,2);
  IMnpipi_woK0_woSid_won_ratio->Draw("HE");
  cIMnpipi->cd(4);
  IMnpipi_wK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,2);
  IMnpipi_wK0_woSid_won_ratio->Draw("HE");

  /*
  TH1D* Momnpipi_woK0_woSid_won_rdata = (TH1D*)f->Get("Momnpipi_woK0_woSid_won_rdata");
  TH1D* Momnpipi_woK0_woSid_won_mc = (TH1D*)f->Get("Momnpipi_woK0_woSid_won_mc");
  TH1D* Momnpipi_wK0_woSid_won_rdata = (TH1D*)f->Get("Momnpipi_wK0_woSid_won_rdata");
  TH1D* Momnpipi_wK0_woSid_won_mc = (TH1D*)f->Get("Momnpipi_wK0_woSid_won_mc");
  TH1D* Momnpipi_woK0_woSid_won_ratio = (TH1D*)f->Get("Momnpipi_woK0_woSid_won_ratio");
  TH1D* Momnpipi_wK0_woSid_won_ratio = (TH1D*)f->Get("Momnpipi_wK0_woSid_won_ratio");

  TCanvas *cMomnpipi = new TCanvas("cMomnpipi","cMomnpipi",1000,1000);
  cMomnpipi->Divide(2,2);
  cMomnpipi->cd(1);
  Momnpipi_woK0_woSid_won_rdata->Draw("HE");
  Momnpipi_woK0_woSid_won_mc->Draw("HEsame");
  cMomnpipi->cd(2);
  Momnpipi_wK0_woSid_won_rdata->Draw("HE");
  Momnpipi_wK0_woSid_won_mc->Draw("HEsame");
  cMomnpipi->cd(3);
  Momnpipi_woK0_woSid_won_ratio->Draw("HE");
  cMomnpipi->cd(4);
  Momnpipi_wK0_woSid_won_ratio->Draw("HE");
  */

  TH1D* IMpippim_woK0_woSid_won_rdata = (TH1D*)f->Get("IMpippim_woK0_woSid_won_rdata");
  TH1D* IMpippim_woK0_woSid_won_mc = (TH1D*)f->Get("IMpippim_woK0_woSid_won_mc");
  TH1D* IMpippim_wK0_woSid_won_rdata = (TH1D*)f->Get("IMpippim_wK0_woSid_won_rdata");
  TH1D* IMpippim_wK0_woSid_won_mc = (TH1D*)f->Get("IMpippim_wK0_woSid_won_mc");
  TH1D* IMpippim_woK0_woSid_won_ratio = (TH1D*)f->Get("IMpippim_woK0_woSid_won_ratio");
  TH1D* IMpippim_wK0_woSid_won_ratio = (TH1D*)f->Get("IMpippim_wK0_woSid_won_ratio");

  TCanvas *cIMpippim = new TCanvas("cIMpippim","cIMpippim",1000,1000);
  cIMpippim->Divide(2,2);
  cIMpippim->cd(1);
  IMpippim_woK0_woSid_won_rdata->Draw("HE");
  IMpippim_woK0_woSid_won_mc->Draw("HEsame");
  cIMpippim->cd(2);
  IMpippim_wK0_woSid_won_rdata->Draw("HE");
  IMpippim_wK0_woSid_won_mc->Draw("HEsame");
  cIMpippim->cd(3);
  IMpippim_woK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,2);
  IMpippim_woK0_woSid_won_ratio->Draw("HE");
  cIMpippim->cd(4);
  IMpippim_wK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,2);
  IMpippim_wK0_woSid_won_ratio->Draw("HE");


  TH1D* Mompippim_woK0_woSid_won_rdata = (TH1D*)f->Get("Mompippim_woK0_woSid_won_rdata");
  TH1D* Mompippim_woK0_woSid_won_mc = (TH1D*)f->Get("Mompippim_woK0_woSid_won_mc");
  TH1D* Mompippim_wK0_woSid_won_rdata = (TH1D*)f->Get("Mompippim_wK0_woSid_won_rdata");
  TH1D* Mompippim_wK0_woSid_won_mc = (TH1D*)f->Get("Mompippim_wK0_woSid_won_mc");
  TH1D* Mompippim_woK0_woSid_won_ratio = (TH1D*)f->Get("Mompippim_woK0_woSid_won_ratio");
  TH1D* Mompippim_wK0_woSid_won_ratio = (TH1D*)f->Get("Mompippim_wK0_woSid_won_ratio");

  TCanvas *cMompippim = new TCanvas("cMompippim","cMompippim",1000,1000);
  cMompippim->Divide(2,2);
  cMompippim->cd(1);
  Mompippim_woK0_woSid_won_rdata->Draw("HE");
  Mompippim_woK0_woSid_won_mc->Draw("HEsame");
  cMompippim->cd(2);
  Mompippim_wK0_woSid_won_rdata->Draw("HE");
  Mompippim_wK0_woSid_won_mc->Draw("HEsame");
  cMompippim->cd(3);
  Mompippim_woK0_woSid_won_ratio->Draw("HE");
  cMompippim->cd(4);
  Mompippim_wK0_woSid_won_ratio->Draw("HE");

  TH1D* nmom_woK0_woSid_won_rdata = (TH1D*)f->Get("nmom_woK0_woSid_won_rdata");
  TH1D* nmom_woK0_woSid_won_mc = (TH1D*)f->Get("nmom_woK0_woSid_won_mc");
  TH1D* nmom_wK0_woSid_won_rdata = (TH1D*)f->Get("nmom_wK0_woSid_won_rdata");
  TH1D* nmom_wK0_woSid_won_mc = (TH1D*)f->Get("nmom_wK0_woSid_won_mc");
  TH1D* nmom_woK0_woSid_won_ratio = (TH1D*)f->Get("nmom_woK0_woSid_won_ratio");
  TH1D* nmom_wK0_woSid_won_ratio = (TH1D*)f->Get("nmom_wK0_woSid_won_ratio");

  TCanvas *cnmom = new TCanvas("cnmom","cnmom",1000,1000);
  cnmom->Divide(2,2);
  cnmom->cd(1);
  nmom_woK0_woSid_won_rdata->Draw("HE");
  nmom_woK0_woSid_won_mc->Draw("HEsame");
  cnmom->cd(2);
  nmom_wK0_woSid_won_rdata->Draw("HE");
  nmom_wK0_woSid_won_mc->Draw("HEsame");
  cnmom->cd(3);
  nmom_woK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,2);
  nmom_woK0_woSid_won_ratio->Draw("HE");
  cnmom->cd(4);
  nmom_wK0_woSid_won_ratio->GetYaxis()->SetRangeUser(0,2);
  nmom_wK0_woSid_won_ratio->Draw("HE");


}
