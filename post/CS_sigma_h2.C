void CS_sigma_h2(const char *filename="evanaIMsigma_npi_h2_v4_out_iso_nostop.root")
{
  TFile *f = TFile::Open(filename);

  TCanvas *c1 = new TCanvas("c1","c1",1000,800);
  TH2F* MMnpi_IMnpip = (TH2F*)f->Get("MMnpi_IMnpip");
  MMnpi_IMnpip->GetXaxis()->SetRangeUser(1,1.5);
  MMnpi_IMnpip->GetYaxis()->SetRangeUser(-0.6,0.6);
  MMnpi_IMnpip->Draw("colz");
  
  TCanvas *c11 = new TCanvas("c11","c11",1000,800);
  TH2F* MM2npi_IMnpip = (TH2F*)f->Get("MM2npi_IMnpip");
  MM2npi_IMnpip->GetXaxis()->SetRangeUser(1,1.5);
  MM2npi_IMnpip->GetYaxis()->SetRangeUser(-0.6,0.4);
  MM2npi_IMnpip->Draw("colz");

  TCanvas *c2 = new TCanvas("c2","c2",1000,800);
  TH2F* MMpi_IMnpip_pi = (TH2F*)f->Get("MMpi_IMnpip_pi");
  MMpi_IMnpip_pi->GetXaxis()->SetRangeUser(1,1.5);
  MMpi_IMnpip_pi->GetYaxis()->SetRangeUser(1,1.7);
  MMpi_IMnpip_pi->Draw("colz");
  
  TCanvas *c2_2 = new TCanvas("c2_2","c2_2",1000,800);
  const int sigma_low = MMpi_IMnpip_pi->GetXaxis()->FindBin(1.17);
  const int sigma_hi = MMpi_IMnpip_pi->GetXaxis()->FindBin(1.20);
  TH1D* MMpi_sp = MMpi_IMnpip_pi->ProjectionY("MMpi_sp",sigma_low,sigma_hi);
  MMpi_sp->Draw("E");


  TCanvas *c3 = new TCanvas("c3","c3",1000,800);
  TH2F* MMn_IMnpip_pi = (TH2F*)f->Get("MMn_IMnpip_pi");
  MMn_IMnpip_pi->GetXaxis()->SetRangeUser(1,1.5);
  MMn_IMnpip_pi->GetYaxis()->SetRangeUser(0,0.7);
  MMn_IMnpip_pi->Draw("colz");

  TCanvas *c4 = new TCanvas("c4","c4",1000,800);
  TH2F* Cospicm_IMnpip_pi = (TH2F*)f->Get("Cospicm_IMnpip_pi");
  Cospicm_IMnpip_pi->GetXaxis()->SetRangeUser(1,1.5);
  //Cospicm_IMnpip->GetYaxis()->SetRangeUser(0,0.7);
  Cospicm_IMnpip_pi->Draw("colz");


  TCanvas *c5 = new TCanvas("c5","c5",1000,800);
  TH2F* MMnpi_IMnpim = (TH2F*)f->Get("MMnpi_IMnpim");
  MMnpi_IMnpim->GetXaxis()->SetRangeUser(1,1.5);
  MMnpi_IMnpim->GetYaxis()->SetRangeUser(-0.6,0.6);
  MMnpi_IMnpim->Draw("colz");
  
  TCanvas *c55 = new TCanvas("c55","c55",1000,800);
  TH2F* MM2npi_IMnpim = (TH2F*)f->Get("MM2npi_IMnpim");
  MM2npi_IMnpim->GetXaxis()->SetRangeUser(1,1.5);
  MM2npi_IMnpim->GetYaxis()->SetRangeUser(-0.6,0.4);
  MM2npi_IMnpim->Draw("colz");

  TCanvas *c6 = new TCanvas("c6","c6",1000,800);
  TH2F* MMpi_IMnpim = (TH2F*)f->Get("MMpi_IMnpim");
  MMpi_IMnpim->GetXaxis()->SetRangeUser(1,1.5);
  MMpi_IMnpim->GetYaxis()->SetRangeUser(1,1.7);
  MMpi_IMnpim->Draw("colz");

  TCanvas *c7 = new TCanvas("c7","c7",1000,800);
  TH2F* MMn_IMnpim = (TH2F*)f->Get("MMn_IMnpim");
  MMn_IMnpim->GetXaxis()->SetRangeUser(1,1.5);
  MMn_IMnpim->GetYaxis()->SetRangeUser(0,0.7);
  MMn_IMnpim->Draw("colz");

  TCanvas *c8 = new TCanvas("c8","c8",1000,800);
  TH2F* Cospicm_IMnpim_pi = (TH2F*)f->Get("Cospicm_IMnpim_pi");
  Cospicm_IMnpim_pi->Draw("colz");




}
