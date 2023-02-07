
const int Version = 245;

void plotOverlap()
{
  TFile *fr = TFile::Open(Form("evanaIMpisigma_npippim_v%d_out_dE2_iso_nostop_sys0_sub.root",Version),"READ");
  TFile *flo = TFile::Open(Form("evanaIMpisigma_npippim_v%d_out_dE2_iso_qlo_nostop_sys0_sub.root",Version),"READ");
  TFile *fhi = TFile::Open(Form("evanaIMpisigma_npippim_v%d_out_dE2_iso_qhi_nostop_sys0_sub.root",Version),"READ");
 
  TH2F* IMnpim_IMnpip_dE_wK0orwSid_n = (TH2F*)fr->Get("IMnpim_IMnpip_dE_wK0orwSid_n");
  gStyle->SetOptStat(0);   
  IMnpim_IMnpip_dE_wK0orwSid_n->Rebin2D(4,4);
  const bool gridon=false;
  gStyle->SetPadGridX(gridon);
  gStyle->SetPadGridY(gridon);
  gStyle->SetTitleYOffset(1.6);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetPadRightMargin(0.15);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.17);
  gStyle->SetPadLeftMargin(0.17);
  //gStyle->SetNdivisions(505,"x");
  gStyle->SetNdivisions(505,"y");
  
  TCanvas *c1 = new TCanvas("c1","",1000,800);
  IMnpim_IMnpip_dE_wK0orwSid_n->SetMinimum(0);
  IMnpim_IMnpip_dE_wK0orwSid_n->SetTitle("");
  //IMnpim_IMnpip_dE_wK0orwSid_n->GetXaxis()->SetRangeUser(0.9,1.7);
  //IMnpim_IMnpip_dE_wK0orwSid_n->GetYaxis()->SetRangeUser(0.9,1.7);
  IMnpim_IMnpip_dE_wK0orwSid_n->GetXaxis()->SetTitleSize(0.06);
  IMnpim_IMnpip_dE_wK0orwSid_n->GetYaxis()->SetTitleSize(0.06);
  IMnpim_IMnpip_dE_wK0orwSid_n->GetXaxis()->SetLabelSize(0.05);
  IMnpim_IMnpip_dE_wK0orwSid_n->GetYaxis()->SetLabelSize(0.05);
  IMnpim_IMnpip_dE_wK0orwSid_n->GetXaxis()->SetTitleOffset(1.31);
  IMnpim_IMnpip_dE_wK0orwSid_n->GetYaxis()->SetTitleOffset(1.42);
  IMnpim_IMnpip_dE_wK0orwSid_n->Draw("colz");

  TH2F* IMnpim_IMnpip_dE_wK0orwSid_n_qlo = (TH2F*)flo->Get("IMnpim_IMnpip_dE_wK0orwSid_n");
  IMnpim_IMnpip_dE_wK0orwSid_n_qlo->Rebin2D(2,2);
  TCanvas *c2 = new TCanvas("c2","",1000,800);
  IMnpim_IMnpip_dE_wK0orwSid_n_qlo->SetMinimum(0);
  IMnpim_IMnpip_dE_wK0orwSid_n_qlo->SetTitle("");
  IMnpim_IMnpip_dE_wK0orwSid_n_qlo->GetXaxis()->SetTitleSize(0.06);
  IMnpim_IMnpip_dE_wK0orwSid_n_qlo->GetYaxis()->SetTitleSize(0.06);
  IMnpim_IMnpip_dE_wK0orwSid_n_qlo->GetXaxis()->SetLabelSize(0.05);
  IMnpim_IMnpip_dE_wK0orwSid_n_qlo->GetYaxis()->SetLabelSize(0.05);
  IMnpim_IMnpip_dE_wK0orwSid_n_qlo->GetXaxis()->SetTitleOffset(1.31);
  IMnpim_IMnpip_dE_wK0orwSid_n_qlo->GetYaxis()->SetTitleOffset(1.42);
  IMnpim_IMnpip_dE_wK0orwSid_n_qlo->Draw("colz");

  TH2F* IMnpim_IMnpip_dE_wK0orwSid_n_qhi = (TH2F*)fhi->Get("IMnpim_IMnpip_dE_wK0orwSid_n");
  IMnpim_IMnpip_dE_wK0orwSid_n_qhi->Rebin2D(2,2);
  TCanvas *c3 = new TCanvas("c3","",1000,800);
  IMnpim_IMnpip_dE_wK0orwSid_n_qhi->SetMinimum(0);
  IMnpim_IMnpip_dE_wK0orwSid_n_qhi->SetTitle("");
  IMnpim_IMnpip_dE_wK0orwSid_n_qhi->GetXaxis()->SetTitleSize(0.06);
  IMnpim_IMnpip_dE_wK0orwSid_n_qhi->GetYaxis()->SetTitleSize(0.06);
  IMnpim_IMnpip_dE_wK0orwSid_n_qhi->GetXaxis()->SetLabelSize(0.05);
  IMnpim_IMnpip_dE_wK0orwSid_n_qhi->GetYaxis()->SetLabelSize(0.05);
  IMnpim_IMnpip_dE_wK0orwSid_n_qhi->GetXaxis()->SetTitleOffset(1.31);
  IMnpim_IMnpip_dE_wK0orwSid_n_qhi->GetYaxis()->SetTitleOffset(1.42);
  IMnpim_IMnpip_dE_wK0orwSid_n_qhi->Draw("colz");

}
