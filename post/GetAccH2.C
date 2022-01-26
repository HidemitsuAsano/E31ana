void GetAccH2()
{
  TFile *fSp = TFile::Open("../simpost/simIMsigma_H2_Sppim_npi_v1_out_iso_rej_nostop.root");
  TFile *fSm = TFile::Open("../simpost/simIMsigma_H2_Smpip_npi_v1_out_iso_rej_nostop.root");


  TH2F* Cospicm_IMnpip_pi_Sp = (TH2F*)fSp->Get("Cospicm_IMnpip_pi");
  TH2F* Cospicm_IMnpip_pi_Sm = (TH2F*)fSm->Get("Cospicm_IMnpim_pi");



}
