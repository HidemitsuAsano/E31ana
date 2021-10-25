#include <../post/anacuts.h>

void plot_S0Pim()
{
  gStyle->SetOptStat("te");
  TFile *_file0 = TFile::Open("simIMLpim_pS0pim_v2.root");
  TCanvas *cReact_q_IMS0Pim = new TCanvas("cReact_q_IMS0Pim","cReact_q_IMS0Pim",1000,800);
  TH2D* React_q_IMS0Pim = _file0->Get("React_q_IMS0Pim");
  React_q_IMS0Pim->SetXTitle("true IM(#Sigma^{0}#pi^{-}) [GeV/c^{2}]");
  React_q_IMS0Pim->SetYTitle("true Mom. Transfer [GeV/c]");
  React_q_IMS0Pim->GetXaxis()->CenterTitle();
  React_q_IMS0Pim->GetYaxis()->CenterTitle();
  React_q_IMS0Pim->RebinX(15);
  React_q_IMS0Pim->RebinY(3);
  React_q_IMS0Pim->Draw("colz");

  TCanvas *cMMass = new TCanvas("cMMass","cMMass",1000,800);
  TFile *freco = TFile::Open("simIMLpim_ppimpim_pS0pim_v2_out.root");
  TH1D* MMass_wL_or = (TH1D*)freco->Get("MMass_wL_or");
  
  MMass_wL_or->Draw("HE");
  double mmax = MMass_wL_or->GetMaximum();
  TLine *phigh_narrow = new TLine(anacuts::Proton_MAX_narrow,0,anacuts::Proton_MAX_narrow,mmax);
  phigh_narrow->SetLineColor(5);
  phigh_narrow->SetLineWidth(2.0);
  phigh_narrow->SetLineStyle(10);
  phigh_narrow->Draw();
  TLine *plow_narrow = new TLine(anacuts::Proton_MIN_narrow,0,anacuts::Proton_MIN_narrow,mmax);
  plow_narrow->SetLineColor(5);
  plow_narrow->SetLineWidth(2.0);
  plow_narrow->SetLineStyle(10);
  plow_narrow->Draw();

  TCanvas *cq_IMppipi_p_wL_sum  = new TCanvas("cq_IMppipi_p_wL_sum","cq_IMppipi_p_wL_sum",1000,800);
  TH2F* q_IMppipi_p_wL_sum = (TH2F*)freco->Get("q_IMppipi_p_wL_sum");
  q_IMppipi_p_wL_sum->Draw("colz");

  TFile *facc = TFile::Open("../simpost/accmapLpimv21.root","READ");
  TH2F* q_IMppipi_p_wL_acc = (TH2F*)facc->Get("q_IMppipi_p_wL_acc");
   
  TH2F* CS_q_IMppipi_p_wL_sum = (TH2F*)q_IMppipi_p_wL_sum->Clone("CS_q_IMppipi_p_wL_sum");
  CS_q_IMppipi_p_wL_sum->Divide(q_IMppipi_p_wL_acc);
  TCanvas *cCS = new TCanvas("cCS","cCS",1000,800);
//  CS_q_IMppipi_p_wL_sum->SetMaximum(0.02);
  CS_q_IMppipi_p_wL_sum->SetXTitle("IM(#pi^{-}#Lambda) [GeV/c^{2}]");
  CS_q_IMppipi_p_wL_sum->SetMaximum(30000);
  CS_q_IMppipi_p_wL_sum->Draw("colz");
  

  TCanvas *cCS_q = new TCanvas("cCS_q","cCS_q",1000,800);
  const int bin1360 = CS_q_IMppipi_p_wL_sum->GetXaxis()->FindBin(1.360);
  const int bin1410 = CS_q_IMppipi_p_wL_sum->GetXaxis()->FindBin(1.410);
  std::cout << CS_q_IMppipi_p_wL_sum->GetXaxis()->GetBinLowEdge(bin1360) << std::endl;
  std::cout << CS_q_IMppipi_p_wL_sum->GetXaxis()->GetBinLowEdge(bin1410+1) << std::endl;
  TH1D* CS_q = (TH1D*)CS_q_IMppipi_p_wL_sum->ProjectionY("CS_q",bin1360,bin1410);
  TH1D* CS_q_true = (TH1D*)React_q_IMS0Pim->ProjectionY("CS_q_true",bin1360,bin1410);
  CS_q->SetMarkerStyle(20);
  CS_q->SetYTitle("d#rho/dM [#mu b (MeV/c^{2})]");
  CS_q->Draw("E");
  CS_q_true->SetLineColor(2);
  CS_q_true->Draw("Esame");

}
