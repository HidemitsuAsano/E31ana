//plot woSid 3 plots
void plot_woSid(const char *filename="evanaIMpisigma_npippim_v202_out_iso.root")
{
  gStyle->SetOptStat(0);
  TFile *f = new TFile(filename,"READ");
  
  TCanvas *c1 = new TCanvas("c1","c1",800,800);
  TH2D* MMnmiss_IMnpip_dE_woSid = (TH2D*)f->Get("MMnmiss_IMnpip_dE_woSid");
  //MMnmiss_IMnpip_dE_woSid->RebinY(2);
  MMnmiss_IMnpip_dE_woSid->Draw("colz");
    
  TCanvas *c2 = new TCanvas("c2","c2",800,800);
  TH2D* MMnmiss_IMnpim_dE_woSid = (TH2D*)f->Get("MMnmiss_IMnpim_dE_woSid");
  //MMnmiss_IMnpim_dE_woSid->RebinY(2);
  MMnmiss_IMnpim_dE_woSid->Draw("colz");
    
  TCanvas *c3 = new TCanvas("c3","c3",800,800);
  TH2D* IMnpim_IMnpip_dE_woSid = (TH2D*)f->Get("IMnpim_IMnpip_dE_woSid");
  IMnpim_IMnpip_dE_woSid->Draw("colz");



}
