void plot_miss()
{

  TFile *_file0 = TFile::Open("simIMpisigma_npipiL_pippimn_v11_out.root");
  TH2D* MMnmiss_IMnpipi_woK0_wSid_lmiss = (TH2D*) _file0->Get("MMnmiss_IMnpipi_woK0_wSid");
  TH1D* lmiss =  (TH1D*) MMnmiss_IMnpipi_woK0_wSid_lmiss->ProjectionY("lmiss");
  lmiss->SetLineColor(2);
  double nlmiss = lmiss->Integral();

  TFile *_file1 = TFile::Open("simIMpisigma_nSppim_pippimn_v52_out.root");
  TH2D* MMnmiss_IMnpipi_woK0_wSid_Sp = (TH2D*) _file1->Get("MMnmiss_IMnpipi_woK0_wSid");
  TH1D* Spmiss =  (TH1D*) MMnmiss_IMnpipi_woK0_wSid_Sp->ProjectionY("Spmiss");
  Spmiss->SetLineColor(3);
  //Spmiss->Draw("HEsame");
  double nSpmiss = Spmiss->Integral();
  
  TFile *_file2 = TFile::Open("simIMpisigma_nSmpip_pippimn_v52_out.root");
  TH2D* MMnmiss_IMnpipi_woK0_wSid_Sm = (TH2D*) _file2->Get("MMnmiss_IMnpipi_woK0_wSid");
  TH1D* Smmiss =  (TH1D*) MMnmiss_IMnpipi_woK0_wSid_Sm->ProjectionY("Smmiss");
  Smmiss->SetLineColor(4);
  //Smmiss->Draw("HEsame");
  double nSmmiss = Smmiss->Integral();


  TFile *_file3 = TFile::Open("simIMpisigma_nS0pippim_pippimn_v1_out.root");
  TH2D* MMnmiss_IMnpipi_woK0_wSid_S0 = (TH2D*) _file3->Get("MMnmiss_IMnpipi_woK0_wSid");
  TH1D* S0miss =  (TH1D*) MMnmiss_IMnpipi_woK0_wSid_S0->ProjectionY("S0miss");
  S0miss->SetLineColor(6);
  //Smmiss->Draw("HEsame");
  double nS0miss = S0miss->Integral();
  
  
  Spmiss->Scale(nlmiss/nSpmiss/4.*1875./2200.);

  Smmiss->Scale(nlmiss/nSmmiss/4.*1875./2200.);

  TH1D* mcsum = (TH1D*)lmiss->Clone("mcsum");
  mcsum->Add(Spmiss);
  mcsum->Add(Smmiss);
  mcsum->Add(S0miss);
  mcsum->SetLineColor(1);
  mcsum->Draw("HE");

  Spmiss->Draw("HEsame");
  Smmiss->Draw("HEsame");
  lmiss->Draw("HEsame");
  S0miss->Draw("HEsame");
}
