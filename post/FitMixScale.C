
TH1D* MM2npi_fmix_pip;
TH1D* IMnpip_fmix;
TH1D* MM2npi_fmix_pim;
TH1D* IMnpim_fmix;


void FitMixScale()
{
  const int dEcut=2;
  const int Version=9;

  TFile *f = TFile::Open(Form("evanaIMsigma_npi_h2_v%d_out_dE%d_iso_nostop.root",Version,dEcut));
  TFile *fmix = TFile::Open(Form("evanaIMsigma_npi_h2_v%d_MIX_out_dE%d_iso_nostop.root",Version,dEcut));

  TH2F* MM2npi_IMnpip_vici_f = (TH2F*)f->Get("MM2npi_IMnpip_vici");
  TH2F* MM2npi_IMnpim_vici_f = (TH2F*)f->Get("MM2npi_IMnpim_vici"); 

  TH2F* MM2npi_IMnpip_vici_fmix = (TH2F*)fmix->Get("MM2npi_IMnpip_vici");
  TH2F* MM2npi_IMnpim_vici_fmix = (TH2F*)fmix->Get("MM2npi_IMnpim_vici"); 

  TCanvas *c1 = new TCanvas("c1","c1");
  MM2npi_IMnpip_vici_f->Draw("colz");

  TCanvas *c2 = new TCanvas("c2","c2");
  MM2npi_IMnpim_vici_f->Draw("colz");

  TCanvas *c3 = new TCanvas("c3","c3");
  MM2npi_IMnpip_vici_fmix->Draw("colz");

  TCanvas *c4 = new TCanvas("c4","c4");
  MM2npi_IMnpim_vici_fmix->Draw("colz");

  TH1D* MM2npi_f_pip = (TH1D*)MM2npi_IMnpip_vici_f->ProjectionY("MM2npi_f_pip");
  TH1D* IMnpip_f = (TH1D*)MM2npi_IMnpip_vici_f->ProjectionX("IMnpip_f");
  TH1D* MM2npi_f_pim = (TH1D*)MM2npi_IMnpim_vici_f->ProjectionY("MM2npi_f_pim");
  TH1D* IMnpim_f = (TH1D*)MM2npi_IMnpim_vici_f->ProjectionX("IMnpim_f");
  MM2npi_fmix_pip = (TH1D*)MM2npi_IMnpip_vici_fmix->ProjectionY("MM2npi_fmix_pip");
  IMnpip_fmix = (TH1D*)MM2npi_IMnpip_vici_fmix->ProjectionX("IMnpip_fmix");
  MM2npi_fmix_pim = (TH1D*)MM2npi_IMnpim_vici_fmix->ProjectionY("MM2npi_fmix_pim");
  IMnpim_fmix = (TH1D*)MM2npi_IMnpim_vici_fmix->ProjectionX("IMnpim_fmix");
  
  TCanvas *c5 = new TCanvas("c5","c5");
  MM2npi_f_pip->Draw("HE");
  MM2npi_fmix_pip->SetLineColor(2);
  MM2npi_fmix_pip->Draw("HEsame");

}

double fit()
{
  static bool isInit = false;
  if(!isInit){
     = TFile::Open(Form("evanaIMsigma_npi_h2_v%d_MIX_out_dE%d_iso_nostop.root",Version,dEcut));




}

