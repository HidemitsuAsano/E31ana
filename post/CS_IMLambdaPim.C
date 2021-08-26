void CS_IMLambdaPim()
{

  TFile *file = TFile::Open("evanaIMLambdaPim_ppimpim_v10_out.root","READ");
  TFile *facc = TFile::Open("../simpost/accmapLpim.root");
  TFile *flumi = TFile::Open("InteLumi.root");
  TParameter<double>*IntegLumi = flumi->Get("IntegLumi");
  TParameter<double>*Err = flumi->Get("Err");
  double lumi = IntegLumi->GetVal();
  double lumierr = Err->GetVal();
  double trigScale = 2.0;
  std::cout << "Lumi:  " << lumi << std::endl;
  std::cout << "Err:   " << lumierr << std::endl;
  TCanvas *cq_IMppipi_p_wL_sum  = new TCanvas("cq_IMppipi_p_wL_sum","cq_IMppipi_p_wL_sum",800,800);
  TH2F* q_IMppipi_p_wL_sum = (TH2F*)file->Get("q_IMppipi_p_wL_sum");
  q_IMppipi_p_wL_sum->Draw("colz");

  TCanvas *cacc = new TCanvas("cacc","cacc",800,800);
  TH2F* q_IMppipi_p_wL_acc = (TH2F*)facc->Get("q_IMppipi_p_wL_acc");
  q_IMppipi_p_wL_acc->Draw("colz");  


  TH2F* CS_q_IMppipi_p_wL_sum = (TH2F*)q_IMppipi_p_wL_sum->Clone("CS_q_IMppipi_p_wL_sum");
  CS_q_IMppipi_p_wL_sum->Divide(q_IMppipi_p_wL_acc);
  TCanvas *cCS = new TCanvas("cCS","cCS",800,800);
  double binwidth = CS_q_IMppipi_p_wL_sum->ProjectionX()->GetBinWidth(1)*1000.0;
  CS_q_IMppipi_p_wL_sum->Scale(1.0/binwidth/trigScale/lumi);
  CS_q_IMppipi_p_wL_sum->SetMaximum(0.02);
  CS_q_IMppipi_p_wL_sum->Draw("colz");

  TCanvas *cCS_px = new TCanvas("cCS_px","cCS_px",800,800);
  TH1D* CS_IMppipi_p_wL_sum = (TH1D*)CS_q_IMppipi_p_wL_sum->ProjectionX("CS_IMppipi_p_wL_sum");
  CS_IMppipi_p_wL_sum->Draw("HE");
  
  const int bin350 = CS_q_IMppipi_p_wL_sum->GetYaxis()->FindBin(0.35);
  TCanvas *cCS_px_0 = new TCanvas("cCS_px_0","cCS_px_0",800,800);
  TH1D* CS_IMppipi_p_wL_sum_0 = (TH1D*)CS_q_IMppipi_p_wL_sum->ProjectionX("CS_IMppipi_p_wL_sum_0",1,bin350-1);
  CS_IMppipi_p_wL_sum_0->GetXaxis()->SetRangeUser(1.2,1.6);
  CS_IMppipi_p_wL_sum_0->Draw("HE");

  TCanvas *cCS_px_350 = new TCanvas("cCS_px_350","cCS_px_350",800,800);
  TH1D* CS_IMppipi_p_wL_sum_350 = (TH1D*)CS_q_IMppipi_p_wL_sum->ProjectionX("CS_IMppipi_p_wL_sum_350",bin350,600);
  CS_IMppipi_p_wL_sum_350->GetXaxis()->SetRangeUser(1.2,1.6);
  CS_IMppipi_p_wL_sum_350->Draw("HE");

}
