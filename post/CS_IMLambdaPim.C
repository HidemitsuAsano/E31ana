void CS_IMLambdaPim()
{

  TFile *file = TFile::Open("evanaIMLambdaPim_ppimpim_v10_out.root","READ");
  TFile *facc = TFile::Open("../simpost/accmapLpim.root");

  
  TCanvas *cq_IMppipi_p_wL_sum  = new TCanvas("cq_IMppipi_p_wL_sum","cq_IMppipi_p_wL_sum",800,800);
  TH2F* q_IMppipi_p_wL_sum = (TH2F*)file->Get("q_IMppipi_p_wL_sum");
  q_IMppipi_p_wL_sum->Draw("colz");

  TCanvas *cacc = new TCanvas("cacc","cacc",800,800);
  TH2F* q_IMppipi_p_wL_acc = (TH2F*)facc->Get("q_IMppipi_p_wL_acc");
  q_IMppipi_p_wL_acc->Draw("colz");  


  TH2F* CS_q_IMppipi_p_wL_sum = (TH2F*)q_IMppipi_p_wL_sum->Clone("CS_q_IMppipi_p_wL_sum");
  CS_q_IMppipi_p_wL_sum->Divide(q_IMppipi_p_wL_acc);
  TCanvas *cCS = new TCanvas("cCS","cCS",800,800);
  CS_q_IMppipi_p_wL_sum->SetMaximum(5000);
  CS_q_IMppipi_p_wL_sum->Draw("colz");


}
