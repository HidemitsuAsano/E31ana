void Fit2DK0()
{
  TFile *fr = TFile::Open("evanaIMpisigma_npippim_v202_out_iso_qhi_sub.root","READ");

  TH2F* IMnpim_IMnpip_dE_wK0_woSid_n = (TH2F*)fr->Get("IMnpim_IMnpip_dE_wK0_woSid_n");
  TCanvas *cIMnpim_IMnpip_dE_wK0_woSid_n = new TCanvas("cIMnpim_IMnpip_dE_wK0_woSid_n","cIMnpim_IMnpip_dE_wK0_woSid_n");
  IMnpim_IMnpip_dE_wK0_woSid_n->RebinX(5);
  IMnpim_IMnpip_dE_wK0_woSid_n->RebinY(5);
  IMnpim_IMnpip_dE_wK0_woSid_n->Draw("colz");














}
