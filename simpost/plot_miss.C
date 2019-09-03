void plot_miss()
{

  TFile *_file0 = TFile::Open("simIMpisigma_npipiL_pippimn_v16_out.root");
  TH2D* MMnmiss_IMnpipi_woK0_wSid_lmiss = (TH2D*) _file0->Get("MMnmiss_IMnpipi_woK0_wSid");
  TH1D* lmiss =  (TH1D*) MMnmiss_IMnpipi_woK0_wSid_lmiss->ProjectionY("lmiss");
  lmiss->SetLineColor(4);
  double nlmiss = lmiss->Integral();

  TFile *_file1 = TFile::Open("simIMpisigma_nSppim_pippimn_v72_out.root");
  TH2D* MMnmiss_IMnpipi_woK0_wSid_Sp = (TH2D*) _file1->Get("MMnmiss_IMnpipi_woK0_wSid");
  TH1D* Spmiss =  (TH1D*) MMnmiss_IMnpipi_woK0_wSid_Sp->ProjectionY("Spmiss");
  Spmiss->SetLineColor(2);
  //Spmiss->Draw("HEsame");
  double nSpmiss = Spmiss->Integral();
  
  TFile *_file2 = TFile::Open("simIMpisigma_nSmpip_pippimn_v72_out.root");
  TH2D* MMnmiss_IMnpipi_woK0_wSid_Sm = (TH2D*) _file2->Get("MMnmiss_IMnpipi_woK0_wSid");
  TH1D* Smmiss =  (TH1D*) MMnmiss_IMnpipi_woK0_wSid_Sm->ProjectionY("Smmiss");
  Smmiss->SetLineColor(3);
  //Smmiss->Draw("HEsame");
  double nSmmiss = Smmiss->Integral();


  TFile *_file3 = TFile::Open("simIMpisigma_nS0pippim_pippimn_v5_out.root");
  TH2D* MMnmiss_IMnpipi_woK0_wSid_S0 = (TH2D*) _file3->Get("MMnmiss_IMnpipi_woK0_wSid");
  TH1D* S0miss =  (TH1D*) MMnmiss_IMnpipi_woK0_wSid_S0->ProjectionY("S0miss");
  S0miss->SetLineColor(5);
  //Smmiss->Draw("HEsame");
  double nS0miss = S0miss->Integral();
  S0miss->Scale(0.9/2.0*0.8*0.49);
  
  Spmiss->Scale(nlmiss/nSpmiss/4.*1875./2200./2.0*0.92);

  Smmiss->Scale(nlmiss/nSmmiss/4.*1875./2200./2.0*0.92);
  lmiss->Scale(1.0/2.0*0.60);

  TFile *_file5 = TFile::Open("simIMpisigma_nSmpippi0_pippimn_v2_out.root");
  TH2D* MMnmiss_IMnpipi_woK0_wSid_Smpi0 = (TH2D*) _file5->Get("MMnmiss_IMnpipi_woK0_wSid");
  TH1D* Smpi0 =  (TH1D*) MMnmiss_IMnpipi_woK0_wSid_Smpi0->ProjectionY("Smpi0");
  Smpi0->SetLineColor(7);
  Smpi0->Scale(0.04);

  TFile *_file6 = TFile::Open("simIMpisigma_nSppimpi0_pippimn_v2_out.root");
  TH2D* MMnmiss_IMnpipi_woK0_wSid_Sppi0 = (TH2D*) _file6->Get("MMnmiss_IMnpipi_woK0_wSid");
  TH1D* Sppi0 =  (TH1D*) MMnmiss_IMnpipi_woK0_wSid_Sppi0->ProjectionY("Sppi0");
  Sppi0->SetLineColor(9);
  Sppi0->Scale(0.05);

  TH1D* mcsum = (TH1D*)lmiss->Clone("mcsum");
  mcsum->Add(Spmiss);
  mcsum->Add(Smmiss);
  mcsum->Add(S0miss);
  mcsum->Add(Smpi0);
  mcsum->Add(Sppi0);
  mcsum->SetLineColor(6);
  //TFile *_file4 = TFile::Open("../post/evanaIMpisigma_npippim_v156_out.root");
  TFile *_file4 = TFile::Open("../post/evanaIMpisigma_npippim_v162_out.root");
  TH2D* MMnmiss_IMnpipi_woK0_wSid_rdata = (TH2D*) _file4->Get("MMnmiss_IMnpipi_woK0_wSid");
  TH1D* rdata = (TH1D*) MMnmiss_IMnpipi_woK0_wSid_rdata->ProjectionY("real");
  rdata->SetLineColor(1);
  //rdata->Scale(3000./1600.);
  //rdata->Scale(3000./1600.);
  //rdata->Scale(1000./1600.);
  //rdata->Scale(1.2);
  //rdata->Scale(0.98);
  //rdata->Scale(0.99);
  //rdata->Scale(0.99);
  rdata->RebinX(2);
  mcsum->RebinX(2);
  Spmiss->RebinX(2);
  Smmiss->RebinX(2);
  lmiss->RebinX(2);
  S0miss->RebinX(2);
  Sppi0->RebinX(2);
  Smpi0->RebinX(2);
  rdata->Draw("HE");
  mcsum->Draw("HEsame");

  Spmiss->Draw("HEsame");
  Smmiss->Draw("HEsame");
  lmiss->Draw("HEsame");
  S0miss->Draw("HEsame");
  Sppi0->Draw("HEsame");
  Smpi0->Draw("HEsame");
}
