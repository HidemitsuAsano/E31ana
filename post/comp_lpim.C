void comp_lpim()
{

  TFile *fdata = TFile::Open("evanaIMLambdaPim_ppimpim_v13_out.root","READ");
  TFile *fsim = TFile::Open("../simpost/simIMLpim_ppimpim_v15_out.root","READ");


  TH2F* IMppim1_IMppim2_data = (TH2F*)fdata->Get("IMppim1_IMppim2");
  TH2F* IMppim1_IMppim2_sim = (TH2F*)fsim->Get("IMppim1_IMppim2");
  TH1D* IMppim1_data = (TH1D*)IMppim1_IMppim2_data->ProjectionY("IMppim1_data");
  TH1D* IMppim1_sim = (TH1D*)IMppim1_IMppim2_sim->ProjectionY("IMppim1_sim");
  TH1D* IMppim2_data = (TH1D*)IMppim1_IMppim2_data->ProjectionX("IMppim2_data");
  TH1D* IMppim2_sim = (TH1D*)IMppim1_IMppim2_sim->ProjectionX("IMppim2_sim");

  TCanvas *clpim = new TCanvas("clpim","clpim");
  IMppim1_data->Draw("HE");
  IMppim1_sim->SetLineColor(2);
  double maxdata = IMppim1_data->GetMaximum();
  double maxsim = IMppim1_sim->GetMaximum();
  IMppim1_sim->Scale(maxdata/maxsim);
  IMppim1_sim->Draw("HEsame");

  TCanvas *clpim2 = new TCanvas("clpim2","clpim2");
  IMppim2_data->Draw("HE");
  IMppim2_sim->SetLineColor(2);
  double maxdata2 = IMppim2_data->GetMaximum();
  double maxsim2 = IMppim2_sim->GetMaximum();
  IMppim2_sim->Scale(maxdata2/maxsim2);
  IMppim2_sim->Draw("HEsame");


  TCanvas *cpcos = new TCanvas("cpcos","cpcos");
  TH1D* pcos_data = (TH1D*)fdata->Get("pcos");
  TH1D* pcos_sim = (TH1D*)fsim->Get("pcos");
  pcos_data->Draw("HE");
  double maxdata3 = pcos_data->GetMaximum();
  double maxsim3 = pcos_sim->GetMaximum();
  pcos_sim->Scale(maxdata3/maxsim3);
  pcos_sim->SetLineColor(2);
  pcos_sim->Draw("HEsame");

  TCanvas *cpmisscos = new TCanvas("cpmisscos","cpmisscos");
  TH1D* pmisscos_data = (TH1D*)fdata->Get("pmisscos");
  TH1D* pmisscos_sim = (TH1D*)fsim->Get("pmisscos");
  pmisscos_data->Draw("HE");
  double maxdata4 = pmisscos_data->GetMaximum();
  double maxsim4 = pmisscos_sim->GetMaximum();
  pmisscos_sim->Scale(maxdata4/maxsim4);
  pmisscos_sim->SetLineColor(2);
  pmisscos_sim->Draw("HEsame");

  TCanvas *cpmom = new TCanvas("cpmom","cpmom");
  TH2D* q_PMom_data = (TH2D*)fdata->Get("q_PMom");
  TH2D* q_PMom_sim = (TH2D*)fsim->Get("q_PMom");
  TH1D* PMom_data = (TH1D*)q_PMom_data->ProjectionX("PMom_data");
  TH1D* PMom_sim = (TH1D*)q_PMom_sim->ProjectionX("PMom_sim");
  PMom_data->Draw("HE");
  double maxdata5 = PMom_data->GetMaximum();
  double maxsim5 = PMom_sim->GetMaximum();
  PMom_sim->Scale(maxdata5/maxsim5);
  PMom_sim->SetLineColor(2);
  PMom_sim->Draw("HEsame");

  
  TCanvas *cpphi


}
