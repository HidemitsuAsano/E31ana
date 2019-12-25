void plot_miss2D()
{
  TH1::SetDefaultSumw2();
  //Lambda
  TFile *_file0 = TFile::Open("simIMpisigma_npipiL_pippimn_v23_out.root");
  //Sp
  TH2D* MMnmiss_IMnpipi_woK0_wSid_lmiss_Sp = (TH2D*) _file0->Get("MMnmiss_IMnpipi_woK0_wSid_Sp");
  TH1D* lmiss_Sp =  (TH1D*) MMnmiss_IMnpipi_woK0_wSid_lmiss_Sp->ProjectionY("lmiss_Sp");
  lmiss_Sp->SetLineColor(4);
  double nlmiss_Sp = lmiss_Sp->Integral();
  //Sm
  TH2D* MMnmiss_IMnpipi_woK0_wSid_lmiss_Sm = (TH2D*) _file0->Get("MMnmiss_IMnpipi_woK0_wSid_Sm");
  TH1D* lmiss_Sm =  (TH1D*) MMnmiss_IMnpipi_woK0_wSid_lmiss_Sm->ProjectionY("lmiss_Sm");
  lmiss_Sm->SetLineColor(4);
  double nlmiss_Sm = lmiss_Sm->Integral();
  

  //Sigma+
  TFile *_file1 = TFile::Open("simIMpisigma_nSppim_pippimn_v107_out.root");
  //Sp
  TH2D* MMnmiss_IMnpipi_woK0_wSid_Sigmap_Sp = (TH2D*) _file1->Get("MMnmiss_IMnpipi_woK0_wSid_Sp");
  TH1D* Sigmap_Sp =  (TH1D*) MMnmiss_IMnpipi_woK0_wSid_Sigmap_Sp->ProjectionY("Sigmap_Sp");
  Sigmap_Sp->SetLineColor(2);
  double nSigmap_Sp = Sigmap_Sp->Integral();
  //Sm
  TH2D* MMnmiss_IMnpipi_woK0_wSid_Sigmap_Sm = (TH2D*) _file1->Get("MMnmiss_IMnpipi_woK0_wSid_Sm");
  TH1D* Sigmap_Sm =  (TH1D*) MMnmiss_IMnpipi_woK0_wSid_Sigmap_Sm->ProjectionY("Sigmap_Sm");
  Sigmap_Sm->SetLineColor(2);
  double nSigmap_Sm = Sigmap_Sm->Integral();


  //Sigma-
  TFile *_file2 = TFile::Open("simIMpisigma_nSmpip_pippimn_v107_out.root");
  //Sp
  TH2D* MMnmiss_IMnpipi_woK0_wSid_Sigmam_Sp = (TH2D*) _file2->Get("MMnmiss_IMnpipi_woK0_wSid_Sp");
  TH1D* Sigmam_Sp =  (TH1D*) MMnmiss_IMnpipi_woK0_wSid_Sigmam_Sp->ProjectionY("Sigmam_Sp");
  Sigmam_Sp->SetLineColor(3);
  double nSigmam_Sp = Sigmam_Sp->Integral();
  //Sm
  TH2D* MMnmiss_IMnpipi_woK0_wSid_Sigmam_Sm = (TH2D*) _file2->Get("MMnmiss_IMnpipi_woK0_wSid_Sm");
  TH1D* Sigmam_Sm =  (TH1D*) MMnmiss_IMnpipi_woK0_wSid_Sigmam_Sm->ProjectionY("Sigmam_Sm");
  Sigmam_Sm->SetLineColor(3);
  double nSigmam_Sm = Sigmam_Sm->Integral();
  

  //nS0pi+pi-
  TFile *_file3 = TFile::Open("simIMpisigma_nS0pippim_pippimn_v5_out.root");
  //Sp
  TH2D* MMnmiss_IMnpipi_woK0_wSid_S0_Sp = (TH2D*) _file3->Get("MMnmiss_IMnpipi_woK0_wSid_Sp");
  TH1D* S0miss_Sp =  (TH1D*) MMnmiss_IMnpipi_woK0_wSid_S0_Sp->ProjectionY("S0miss_Sp");
  S0miss_Sp->SetLineColor(5);
  double nS0miss_Sp = S0miss_Sp->Integral();
  //Sm
  TH2D* MMnmiss_IMnpipi_woK0_wSid_S0_Sm = (TH2D*) _file3->Get("MMnmiss_IMnpipi_woK0_wSid_Sm");
  TH1D* S0miss_Sm =  (TH1D*) MMnmiss_IMnpipi_woK0_wSid_S0_Sm->ProjectionY("S0miss_Sm");
  S0miss_Sm->SetLineColor(5);
  double nS0miss_Sm = S0miss_Sm->Integral();

  
  //nSppi-pi0
  TFile *_file5 = TFile::Open("simIMpisigma_nSppimpi0_pippimn_v2_out.root");
  //Sp
  TH2D* MMnmiss_IMnpipi_woK0_wSid_Sppi0_Sp = (TH2D*) _file5->Get("MMnmiss_IMnpipi_woK0_wSid_Sp");
  TH1D* Sppi0_Sp =  (TH1D*) MMnmiss_IMnpipi_woK0_wSid_Sppi0_Sp->ProjectionY("Sppi0_Sp");
  Sppi0_Sp->SetLineColor(7);
  double nSppi0_Sp= Sppi0_Sp->Integral();
  //Sm
  TH2D* MMnmiss_IMnpipi_woK0_wSid_Sppi0_Sm = (TH2D*) _file5->Get("MMnmiss_IMnpipi_woK0_wSid_Sm");
  TH1D* Sppi0_Sm =  (TH1D*) MMnmiss_IMnpipi_woK0_wSid_Sppi0_Sm->ProjectionY("Sppi0_Sm");
  Sppi0_Sm->SetLineColor(7);
  double nSppi0_Sm = Sppi0_Sm->Integral();
  
  //nSmpi+pi0
  TFile *_file6 = TFile::Open("simIMpisigma_nSmpippi0_pippimn_v2_out.root");
  //Sp
  TH2D* MMnmiss_IMnpipi_woK0_wSid_Smpi0_Sp = (TH2D*) _file6->Get("MMnmiss_IMnpipi_woK0_wSid_Sp");
  TH1D* Smpi0_Sp =  (TH1D*) MMnmiss_IMnpipi_woK0_wSid_Smpi0_Sp->ProjectionY("Smpi0_Sp");
  Smpi0_Sp->SetLineColor(9);
  double nSmpi0_Sp = Smpi0_Sp->Integral();
  //Sm
  TH2D* MMnmiss_IMnpipi_woK0_wSid_Smpi0_Sm = (TH2D*) _file6->Get("MMnmiss_IMnpipi_woK0_wSid_Sm");
  TH1D* Smpi0_Sm =  (TH1D*) MMnmiss_IMnpipi_woK0_wSid_Smpi0_Sm->ProjectionY("Smpi0_Sm");
  Smpi0_Sm->SetLineColor(9);
  double nSmpi0_Sm = Smpi0_Sm->Integral();

  TFile *_file4 = TFile::Open("../post/evanaIMpisigma_npippim_v196_out.root");
  //Sp
  TH2D* MMnmiss_IMnpipi_woK0_wSid_rdata_Sp = (TH2D*) _file4->Get("MMnmiss_IMnpipi_woK0_wSid_Sp");
  TH1D* rdata_Sp = (TH1D*) MMnmiss_IMnpipi_woK0_wSid_rdata_Sp->ProjectionY("real_Sp");
  double nrdata_Sp = rdata_Sp->Integral();
  //Sm
  TH2D* MMnmiss_IMnpipi_woK0_wSid_rdata_Sm = (TH2D*) _file4->Get("MMnmiss_IMnpipi_woK0_wSid_Sm");
  TH1D* rdata_Sm = (TH1D*) MMnmiss_IMnpipi_woK0_wSid_rdata_Sm->ProjectionY("real_Sm");
  double nrdata_Sm = rdata_Sm->Integral();
  
  TCanvas *miss_Sp = new TCanvas("miss_Sp","miss_Sp",800,800);
  miss_Sp->cd();
  lmiss_Sp->Scale(0.40*nrdata_Sp/nlmiss_Sp);
  Sigmap_Sp->Scale(0.18*nrdata_Sp/nSigmap_Sp);
  Sigmam_Sp->Scale(0.01*nrdata_Sp/nSigmam_Sp);
  S0miss_Sp->Scale(0.24*nrdata_Sp/nS0miss_Sp);
  Sppi0_Sp->Scale(0.03*nrdata_Sp/nSppi0_Sp);
  Smpi0_Sp->Scale(0.03*nrdata_Sp/nSmpi0_Sp);
  TH1D* mcsum_Sp = (TH1D*)lmiss_Sp->Clone("mcsum_Sp");
  mcsum_Sp->Add(Sigmap_Sp);
  mcsum_Sp->Add(Sigmam_Sp);
  mcsum_Sp->Add(S0miss_Sp);
  mcsum_Sp->Add(Sppi0_Sp);
  mcsum_Sp->Add(Smpi0_Sp);
  mcsum_Sp->SetLineColor(6);
  //TFile *_file4 = TFile::Open("../post/evanaIMpisigma_npippim_v156_out.root");
  rdata_Sp->SetLineColor(1);
  rdata_Sp->Draw("HE");
  mcsum_Sp->Draw("HEsame");
  Sigmap_Sp->Draw("HEsame");
  Sigmam_Sp->Draw("HEsame");
  lmiss_Sp->Draw("HEsame");
  S0miss_Sp->Draw("HEsame");
  Sppi0_Sp->Draw("HEsame");
  Smpi0_Sp->Draw("HEsame");
  
  TCanvas *miss_Sm = new TCanvas("miss_Sm","miss_Sm",800,0,800,800);
  miss_Sm->cd();
  lmiss_Sm->Scale(0.40*nrdata_Sm/nlmiss_Sm);
  Sigmap_Sm->Scale(0.01*nrdata_Sm/nSigmap_Sm);
  Sigmam_Sm->Scale(0.21*nrdata_Sm/nSigmam_Sm);
  S0miss_Sm->Scale(0.2*nrdata_Sm/nS0miss_Sm);
  Sppi0_Sm->Scale(0.03*nrdata_Sm/nSppi0_Sm);
  Smpi0_Sm->Scale(0.03*nrdata_Sm/nSmpi0_Sm);
  TH1D* mcsum_Sm = (TH1D*)lmiss_Sm->Clone("mcsum_Sm");
  mcsum_Sm->Add(Sigmap_Sm);
  mcsum_Sm->Add(Sigmam_Sm);
  mcsum_Sm->Add(S0miss_Sm);
  mcsum_Sm->Add(Sppi0_Sm);
  mcsum_Sm->Add(Smpi0_Sm);
  mcsum_Sm->SetLineColor(6);
  //TFile *_file4 = TFile::Open("../post/evanaIMpisigma_npippim_v156_out.root");
  rdata_Sm->SetLineColor(1);
  rdata_Sm->Draw("HE");
  mcsum_Sm->Draw("HEsame");
  Sigmap_Sm->Draw("HEsame");
  Sigmam_Sm->Draw("HEsame");
  lmiss_Sm->Draw("HEsame");
  S0miss_Sm->Draw("HEsame");
  Sppi0_Sm->Draw("HEsame");
  Smpi0_Sm->Draw("HEsame");

  //superimpose 2d hists of miss mass vs IM(npi+) or IM(npi-)
  //Lambda
  TH2D *MMnmiss_IMnpip_dE_woK0_lmiss = (TH2D*)_file0->Get("MMnmiss_IMnpip_dE_woK0");
  TH2D *MMnmiss_IMnpim_dE_woK0_lmiss = (TH2D*)_file0->Get("MMnmiss_IMnpim_dE_woK0");
  TH2D *MMnmiss_IMnpip_dE_woK0_woSm_lmiss = (TH2D*)_file0->Get("MMnmiss_IMnpip_dE_woK0_woSm");
  TH2D *MMnmiss_IMnpim_dE_woK0_woSp_lmiss = (TH2D*)_file0->Get("MMnmiss_IMnpim_dE_woK0_woSp");
  TH2D *MMnmiss_IMnpip_dE_woK0_woSidn_lmiss = (TH2D*)_file0->Get("MMnmiss_IMnpip_dE_woK0_woSidn");
  TH2D *MMnmiss_IMnpim_dE_woK0_woSidn_lmiss = (TH2D*)_file0->Get("MMnmiss_IMnpim_dE_woK0_woSidn");
  TH2D *MMnmiss_IMpippim_dE_lmiss = (TH2D*)_file0->Get("MMnmiss_IMpippim_dE");
  TH2D *MMnmiss_IMpippim_dE_wSid_lmiss = (TH2D*)_file0->Get("MMnmiss_IMpippim_wSid_dE");
  TH2D *MMnmiss_IMpippim_dE_woK0_woSidn_lmiss = (TH2D*)_file0->Get("MMnmiss_IMpippim_dE_woK0_woSidn");
  //TH2D *IMnpim_IMnpip_dE_lmiss = (TH2D*)_file0->Get("IMnpim_IMnpip_dE");
  double Nnpip_lmiss = MMnmiss_IMnpip_dE_woK0_lmiss->Integral();
  double Nnpim_lmiss = MMnmiss_IMnpim_dE_woK0_lmiss->Integral();
  //nSppi-
  TH2D *MMnmiss_IMnpip_dE_woK0_Sigmap = (TH2D*)_file1->Get("MMnmiss_IMnpip_dE_woK0");
  TH2D *MMnmiss_IMnpim_dE_woK0_Sigmap = (TH2D*)_file1->Get("MMnmiss_IMnpim_dE_woK0");
  TH2D *MMnmiss_IMnpip_dE_woK0_woSm_Sigmap = (TH2D*)_file1->Get("MMnmiss_IMnpip_dE_woK0_woSm");
  TH2D *MMnmiss_IMnpim_dE_woK0_woSp_Sigmap = (TH2D*)_file1->Get("MMnmiss_IMnpim_dE_woK0_woSp");
  TH2D *MMnmiss_IMnpip_dE_woK0_woSidn_Sigmap = (TH2D*)_file1->Get("MMnmiss_IMnpip_dE_woK0_woSidn");
  TH2D *MMnmiss_IMnpim_dE_woK0_woSidn_Sigmap = (TH2D*)_file1->Get("MMnmiss_IMnpim_dE_woK0_woSidn");
  double Nnpip_Sigmap = MMnmiss_IMnpip_dE_woK0_Sigmap->Integral();
  double Nnpim_Sigmap = MMnmiss_IMnpim_dE_woK0_Sigmap->Integral();
  //nSmpi+
  TH2D *MMnmiss_IMnpip_dE_woK0_Sigmam = (TH2D*)_file2->Get("MMnmiss_IMnpip_dE_woK0");
  TH2D *MMnmiss_IMnpim_dE_woK0_Sigmam = (TH2D*)_file2->Get("MMnmiss_IMnpim_dE_woK0");
  TH2D *MMnmiss_IMnpip_dE_woK0_woSm_Sigmam = (TH2D*)_file2->Get("MMnmiss_IMnpip_dE_woK0_woSm");
  TH2D *MMnmiss_IMnpim_dE_woK0_woSp_Sigmam = (TH2D*)_file2->Get("MMnmiss_IMnpim_dE_woK0_woSp");
  TH2D *MMnmiss_IMnpip_dE_woK0_woSidn_Sigmam = (TH2D*)_file2->Get("MMnmiss_IMnpip_dE_woK0_woSidn");
  TH2D *MMnmiss_IMnpim_dE_woK0_woSidn_Sigmam = (TH2D*)_file2->Get("MMnmiss_IMnpim_dE_woK0_woSidn");
  double Nnpip_Sigmam = MMnmiss_IMnpip_dE_woK0_Sigmam->Integral();
  double Nnpim_Sigmam = MMnmiss_IMnpim_dE_woK0_Sigmam->Integral();
  //Sigma0
  TH2D *MMnmiss_IMnpip_dE_woK0_S0 = (TH2D*)_file3->Get("MMnmiss_IMnpip_dE_woK0");
  TH2D *MMnmiss_IMnpim_dE_woK0_S0 = (TH2D*)_file3->Get("MMnmiss_IMnpim_dE_woK0");
  TH2D *MMnmiss_IMnpip_dE_woK0_woSm_S0 = (TH2D*)_file3->Get("MMnmiss_IMnpip_dE_woK0_woSm");
  TH2D *MMnmiss_IMnpim_dE_woK0_woSp_S0 = (TH2D*)_file3->Get("MMnmiss_IMnpim_dE_woK0_woSp");
  TH2D *MMnmiss_IMnpip_dE_woK0_woSidn_S0 = (TH2D*)_file3->Get("MMnmiss_IMnpip_dE_woK0_woSidn");
  TH2D *MMnmiss_IMnpim_dE_woK0_woSidn_S0 = (TH2D*)_file3->Get("MMnmiss_IMnpim_dE_woK0_woSidn");
  double Nnpip_S0 = MMnmiss_IMnpip_dE_woK0_S0->Integral();
  double Nnpim_S0 = MMnmiss_IMnpim_dE_woK0_S0->Integral();
  //Sppi0
  TH2D *MMnmiss_IMnpip_dE_woK0_Sppi0 = (TH2D*)_file5->Get("MMnmiss_IMnpip_dE_woK0");
  TH2D *MMnmiss_IMnpim_dE_woK0_Sppi0 = (TH2D*)_file5->Get("MMnmiss_IMnpim_dE_woK0");
  double Nnpip_Sppi0 = MMnmiss_IMnpip_dE_woK0_Sppi0->Integral();
  double Nnpim_Sppi0 = MMnmiss_IMnpim_dE_woK0_Sppi0->Integral();
  //Smpi0
  TH2D *MMnmiss_IMnpip_dE_woK0_Smpi0 = (TH2D*)_file6->Get("MMnmiss_IMnpip_dE_woK0");
  TH2D *MMnmiss_IMnpim_dE_woK0_Smpi0 = (TH2D*)_file6->Get("MMnmiss_IMnpim_dE_woK0");
  double Nnpip_Smpi0 = MMnmiss_IMnpip_dE_woK0_Smpi0->Integral();
  double Nnpim_Smpi0 = MMnmiss_IMnpim_dE_woK0_Smpi0->Integral();
  //real data
  TH2D *MMnmiss_IMnpip_dE_woK0_rdata = (TH2D*)_file4->Get("MMnmiss_IMnpip_dE_woK0");
  TH2D *MMnmiss_IMnpim_dE_woK0_rdata = (TH2D*)_file4->Get("MMnmiss_IMnpim_dE_woK0");
  TH2D *MMnmiss_IMnpip_dE_woK0_woSm_rdata = (TH2D*)_file4->Get("MMnmiss_IMnpip_dE_woK0_woSm");
  TH2D *MMnmiss_IMnpim_dE_woK0_woSp_rdata = (TH2D*)_file4->Get("MMnmiss_IMnpim_dE_woK0_woSp");
  TH2D *MMnmiss_IMnpip_dE_woK0_woSidn_rdata = (TH2D*)_file4->Get("MMnmiss_IMnpip_dE_woK0_woSidn");
  TH2D *MMnmiss_IMnpim_dE_woK0_woSidn_rdata = (TH2D*)_file4->Get("MMnmiss_IMnpim_dE_woK0_woSidn");
  double Nnpip_rdata = MMnmiss_IMnpip_dE_woK0_rdata->Integral();
  double Nnpim_rdata = MMnmiss_IMnpim_dE_woK0_rdata->Integral();
  
  //normalization 
  MMnmiss_IMnpip_dE_woK0_lmiss->Scale(0.40*nrdata_Sp/nlmiss_Sp);
  MMnmiss_IMnpip_dE_woK0_Sigmap->Scale(0.18*nrdata_Sp/nSigmap_Sp);
  MMnmiss_IMnpip_dE_woK0_Sigmam->Scale(0.01*nrdata_Sp/nSigmam_Sp);
  MMnmiss_IMnpip_dE_woK0_S0->Scale(0.24*nrdata_Sp/nS0miss_Sp);
  MMnmiss_IMnpip_dE_woK0_Sppi0->Scale(0.03*nrdata_Sp/nSppi0_Sp);
  MMnmiss_IMnpip_dE_woK0_Smpi0->Scale(0.03*nrdata_Sp/nSmpi0_Sp);
  
  MMnmiss_IMnpip_dE_woK0_woSm_lmiss->Scale(0.40*nrdata_Sp/nlmiss_Sp);
  MMnmiss_IMnpip_dE_woK0_woSm_Sigmap->Scale(0.18*nrdata_Sp/nSigmap_Sp);
  MMnmiss_IMnpip_dE_woK0_woSm_Sigmam->Scale(0.01*nrdata_Sp/nSigmam_Sp);
  MMnmiss_IMnpip_dE_woK0_woSm_S0->Scale(0.24*nrdata_Sp/nS0miss_Sp);
  //MMnmiss_IMnpip_dE_woK0_woSm_Sppi0->Scale(0.03*nrdata_Sp/nSppi0_Sp);
  //MMnmiss_IMnpip_dE_woK0_woSm_Smpi0->Scale(0.03*nrdata_Sp/nSmpi0_Sp);
  
  MMnmiss_IMnpip_dE_woK0_woSidn_lmiss->Scale(0.40*nrdata_Sp/nlmiss_Sp);
  MMnmiss_IMnpip_dE_woK0_woSidn_Sigmap->Scale(0.18*nrdata_Sp/nSigmap_Sp);
  MMnmiss_IMnpip_dE_woK0_woSidn_Sigmam->Scale(0.01*nrdata_Sp/nSigmam_Sp);
  MMnmiss_IMnpip_dE_woK0_woSidn_S0->Scale(0.24*nrdata_Sp/nS0miss_Sp);
  //MMnmiss_IMnpip_dE_woK0_wSidn_Sppi0->Scale(0.03*nrdata_Sp/nSppi0_Sp);
  //MMnmiss_IMnpip_dE_woK0_wSidn_Smpi0->Scale(0.03*nrdata_Sp/nSmpi0_Sp);
  

  TH2D *MMnmiss_IMnpip_mc = (TH2D*)MMnmiss_IMnpip_dE_woK0_lmiss->Clone("MMnmiss_IMnpip_mc");
  MMnmiss_IMnpip_mc->Add(MMnmiss_IMnpip_dE_woK0_Sigmap);
  MMnmiss_IMnpip_mc->Add(MMnmiss_IMnpip_dE_woK0_Sigmam);
  MMnmiss_IMnpip_mc->Add(MMnmiss_IMnpip_dE_woK0_S0);
  MMnmiss_IMnpip_mc->Add(MMnmiss_IMnpip_dE_woK0_Sppi0);
  MMnmiss_IMnpip_mc->Add(MMnmiss_IMnpip_dE_woK0_Smpi0);

  TCanvas *cMMnmiss_IMnpip = new TCanvas("cMMnmiss_IMnpip","cMMnmiss_IMnpip",1200,800);
  cMMnmiss_IMnpip->Divide(2,1);
  cMMnmiss_IMnpip->cd(1);
  MMnmiss_IMnpip_dE_woK0_rdata->SetTitle("real data");
  MMnmiss_IMnpip_dE_woK0_rdata->Draw("colz");
  cMMnmiss_IMnpip->cd(2);
  MMnmiss_IMnpip_mc->SetTitle("MC sum");
  MMnmiss_IMnpip_mc->Draw("colz");
  
  //
  TH2D *MMnmiss_IMnpip_woSm_mc = (TH2D*)MMnmiss_IMnpip_dE_woK0_woSm_lmiss->Clone("MMnmiss_IMnpip_woSm_mc");
  MMnmiss_IMnpip_woSm_mc->Add(MMnmiss_IMnpip_dE_woK0_woSm_Sigmap);
  MMnmiss_IMnpip_woSm_mc->Add(MMnmiss_IMnpip_dE_woK0_woSm_Sigmam);
  MMnmiss_IMnpip_woSm_mc->Add(MMnmiss_IMnpip_dE_woK0_woSm_S0);
  //MMnmiss_IMnpip_mc->Add(MMnmiss_IMnpip_dE_woK0_Sppi0);
  //MMnmiss_IMnpip_mc->Add(MMnmiss_IMnpip_dE_woK0_Smpi0);

  TCanvas *cMMnmiss_IMnpip_woSm = new TCanvas("cMMnmiss_IMnpip_woSm","cMMnmiss_IMnpip_woSm",1200,800);
  cMMnmiss_IMnpip_woSm->Divide(2,1);
  cMMnmiss_IMnpip_woSm->cd(1);
  MMnmiss_IMnpip_dE_woK0_woSm_rdata->SetTitle("real data");
  MMnmiss_IMnpip_dE_woK0_woSm_rdata->Draw("colz");
  cMMnmiss_IMnpip_woSm->cd(2);
  MMnmiss_IMnpip_woSm_mc->SetTitle("MC sum");
  MMnmiss_IMnpip_woSm_mc->Draw("colz");
  
  //
  TH2D *MMnmiss_IMnpip_woSidn_mc = (TH2D*)MMnmiss_IMnpip_dE_woK0_woSidn_lmiss->Clone("MMnmiss_IMnpip_woSidn_mc");
  MMnmiss_IMnpip_woSidn_mc->Add(MMnmiss_IMnpip_dE_woK0_woSidn_Sigmap);
  MMnmiss_IMnpip_woSidn_mc->Add(MMnmiss_IMnpip_dE_woK0_woSidn_Sigmam);
  MMnmiss_IMnpip_woSidn_mc->Add(MMnmiss_IMnpip_dE_woK0_woSidn_S0);
  //MMnmiss_IMnpip_mc->Add(MMnmiss_IMnpip_dE_woK0_Sppi0);
  //MMnmiss_IMnpip_mc->Add(MMnmiss_IMnpip_dE_woK0_Smpi0);

  TCanvas *cMMnmiss_IMnpip_woSidn = new TCanvas("cMMnmiss_IMnpip_woSidn","cMMnmiss_IMnpip_woSidn",1200,800);
  cMMnmiss_IMnpip_woSidn->Divide(2,1);
  cMMnmiss_IMnpip_woSidn->cd(1);
  MMnmiss_IMnpip_dE_woK0_woSidn_rdata->SetTitle("real data");
  MMnmiss_IMnpip_dE_woK0_woSidn_rdata->Draw("colz");
  cMMnmiss_IMnpip_woSidn->cd(2);
  MMnmiss_IMnpip_woSidn_mc->SetTitle("MC sum");
  MMnmiss_IMnpip_woSidn_mc->Draw("colz");
  
  //projection to Missing mass (miss n & Sigma+/-)
  TCanvas *cMMnmiss_woSidn = new TCanvas("cMMnmiss_woSidn","cMMnmiss_woSidn");
  cMMnmiss_woSidn->cd();
  MMnmiss_IMnpip_dE_woK0_woSidn_rdata->ProjectionY("MMnmiss_IMnpip_dE_woK0_woSidn_rdata_py")->Draw("HE");
  TH1D* MMnmiss_woSidn_mc = (TH1D*)MMnmiss_IMnpip_woSidn_mc->ProjectionY("MMnmiss_IMnpip_woSidn_mc_py");
  MMnmiss_woSidn_mc->SetLineColor(6);
  MMnmiss_woSidn_mc->Draw("HEsame");

  TH1D* MMnmiss_woSidn_lmiss = (TH1D*)MMnmiss_IMnpip_dE_woK0_woSidn_lmiss->ProjectionY("MMnmiss_IMnpip_dE_woK0_woSidn_lmiss_py");
  MMnmiss_woSidn_lmiss->SetLineColor(4); 
  MMnmiss_woSidn_lmiss->Draw("HEsame"); 
  TH1D* MMnmiss_woSidn_S0 = (TH1D*)MMnmiss_IMnpip_dE_woK0_woSidn_S0->ProjectionY("MMnmiss_IMnpip_dE_woK0_woSidn_S0_py");
  MMnmiss_woSidn_S0->SetLineColor(5); 
  MMnmiss_woSidn_S0->Draw("HEsame");; 
  TH1D* MMnmiss_woSidn_Sigmap = (TH1D*)MMnmiss_IMnpip_dE_woK0_woSidn_Sigmap->ProjectionY("MMnmiss_IMnpip_dE_woK0_woSidn_Sigmap_py");
  MMnmiss_woSidn_Sigmap->SetLineColor(2); 
  MMnmiss_woSidn_Sigmap->Draw("HEsame"); 
  TH1D* MMnmiss_woSidn_Sigmam = (TH1D*)MMnmiss_IMnpip_dE_woK0_woSidn_Sigmam->ProjectionY("MMnmiss_IMnpip_dE_woK0_woSidn_Sigmam_py");
  MMnmiss_woSidn_Sigmam->SetLineColor(3); 
  MMnmiss_woSidn_Sigmam->Draw("HEsame"); 

  TCanvas *cMMnmiss_woSidn_ratio = new TCanvas("cMMnmiss_woSidn_ratio","cMMnmiss_woSid_ratio");
  TH1D* MMnmiss_woSidn_ratio = MMnmiss_woSidn_mc->Clone("MMnmiss_woSidn_ratio");
  MMnmiss_woSidn_ratio->Divide(MMnmiss_IMnpip_dE_woK0_woSidn_rdata_py);
  MMnmiss_woSidn_ratio->Draw("HE");

  
  //projection to IMnpip (miss n & Sigma+/-)
  TCanvas *cIMnpip_woSidn = new TCanvas("cIMnpip_woSidn","cIMnpip_woSidn");
  cIMnpip_woSidn->cd();
  MMnmiss_IMnpip_dE_woK0_woSidn_rdata->ProjectionX("MMnmiss_IMnpip_dE_woK0_woSidn_rdata_px")->Draw("HE");
  TH1D* IMnpip_woSidn_mc = (TH1D*)MMnmiss_IMnpip_woSidn_mc->ProjectionX("MMnmiss_IMnpip_woSidn_mc_px");
  IMnpip_woSidn_mc->SetLineColor(6);
  IMnpip_woSidn_mc->Draw("HEsame");
  TH1D* IMnpip_woSidn_lmiss = (TH1D*)MMnmiss_IMnpip_dE_woK0_woSidn_lmiss->ProjectionX("MMnmiss_IMnpip_dE_woK0_woSidn_lmiss_px");
  IMnpip_woSidn_lmiss->SetLineColor(4); 
  IMnpip_woSidn_lmiss->Draw("HEsame"); 
  TH1D* IMnpip_woSidn_S0 = (TH1D*)MMnmiss_IMnpip_dE_woK0_woSidn_S0->ProjectionX("MMnmiss_IMnpip_dE_woK0_woSidn_S0_px");
  IMnpip_woSidn_S0->SetLineColor(5); 
  IMnpip_woSidn_S0->Draw("HEsame");; 
  TH1D* IMnpip_woSidn_Sigmap = (TH1D*)MMnmiss_IMnpip_dE_woK0_woSidn_Sigmap->ProjectionX("MMnmiss_IMnpip_dE_woK0_woSidn_Sigmap_px");
  IMnpip_woSidn_Sigmap->SetLineColor(2); 
  IMnpip_woSidn_Sigmap->Draw("HEsame"); 
  TH1D* IMnpip_woSidn_Sigmam = (TH1D*)MMnmiss_IMnpip_dE_woK0_woSidn_Sigmam->ProjectionX("MMnmiss_IMnpip_dE_woK0_woSidn_Sigmam_px");
  IMnpip_woSidn_Sigmam->SetLineColor(3); 
  IMnpip_woSidn_Sigmam->Draw("HEsame"); 



  //Sigma-
  MMnmiss_IMnpim_dE_woK0_lmiss->Scale(0.40*nrdata_Sm/nlmiss_Sm);
  MMnmiss_IMnpim_dE_woK0_Sigmap->Scale(0.01*nrdata_Sm/nSigmap_Sm);
  MMnmiss_IMnpim_dE_woK0_Sigmam->Scale(0.21*nrdata_Sm/nSigmam_Sm);
  MMnmiss_IMnpim_dE_woK0_S0->Scale(0.2*nrdata_Sm/nS0miss_Sm);
  MMnmiss_IMnpim_dE_woK0_Sppi0->Scale(0.03*nrdata_Sm/nSppi0_Sm);
  MMnmiss_IMnpim_dE_woK0_Smpi0->Scale(0.03*nrdata_Sm/nSmpi0_Sm);
  
  MMnmiss_IMnpim_dE_woK0_woSp_lmiss->Scale(0.40*nrdata_Sm/nlmiss_Sm);
  MMnmiss_IMnpim_dE_woK0_woSp_Sigmap->Scale(0.01*nrdata_Sm/nSigmap_Sm);
  MMnmiss_IMnpim_dE_woK0_woSp_Sigmam->Scale(0.21*nrdata_Sm/nSigmam_Sm);
  MMnmiss_IMnpim_dE_woK0_woSp_S0->Scale(0.2*nrdata_Sm/nS0miss_Sm);
  //MMnmiss_IMnpim_dE_woK0_woSp_Sppi0->Scale(0.03*nrdata_Sm/nSppi0_Sm);
  //MMnmiss_IMnpim_dE_woK0_woSp_Smpi0->Scale(0.03*nrdata_Sm/nSmpi0_Sm);
  
  MMnmiss_IMnpim_dE_woK0_woSidn_lmiss->Scale(0.40*nrdata_Sm/nlmiss_Sm);
  MMnmiss_IMnpim_dE_woK0_woSidn_Sigmap->Scale(0.01*nrdata_Sm/nSigmap_Sm);
  MMnmiss_IMnpim_dE_woK0_woSidn_Sigmam->Scale(0.21*nrdata_Sm/nSigmam_Sm);
  MMnmiss_IMnpim_dE_woK0_woSidn_S0->Scale(0.2*nrdata_Sm/nS0miss_Sm);
  //MMnmiss_IMnpim_dE_woK0_wSidn_Sppi0->Scale(0.03*nrdata_Sm/nSppi0_Sm);
  //MMnmiss_IMnpim_dE_woK0_wSidn_Smpi0->Scale(0.03*nrdata_Sm/nSmpi0_Sm);

  TH2D *MMnmiss_IMnpim_mc = (TH2D*)MMnmiss_IMnpim_dE_woK0_lmiss->Clone("MMnmiss_IMnpim_mc");
  MMnmiss_IMnpim_mc->Add(MMnmiss_IMnpim_dE_woK0_Sigmap);
  MMnmiss_IMnpim_mc->Add(MMnmiss_IMnpim_dE_woK0_Sigmam);
  MMnmiss_IMnpim_mc->Add(MMnmiss_IMnpim_dE_woK0_S0);
  MMnmiss_IMnpim_mc->Add(MMnmiss_IMnpim_dE_woK0_Sppi0);
  MMnmiss_IMnpim_mc->Add(MMnmiss_IMnpim_dE_woK0_Smpi0);

  TCanvas *cMMnmiss_IMnpim = new TCanvas("cMMnmiss_IMnpim","cMMnmiss_IMnpim",1200,800);
  cMMnmiss_IMnpim->Divide(2,1);
  cMMnmiss_IMnpim->cd(1);
  MMnmiss_IMnpim_dE_woK0_rdata->SetTitle("real data");
  MMnmiss_IMnpim_dE_woK0_rdata->Draw("colz");
  cMMnmiss_IMnpim->cd(2);
  MMnmiss_IMnpim_mc->SetTitle("MC sum");
  MMnmiss_IMnpim_mc->Draw("colz");

  TH2D *MMnmiss_IMnpim_woSidn_mc = (TH2D*)MMnmiss_IMnpim_dE_woK0_woSidn_lmiss->Clone("MMnmiss_IMnpim_woSidn_mc");
  MMnmiss_IMnpim_woSidn_mc->Add(MMnmiss_IMnpim_dE_woK0_woSidn_Sigmap);
  MMnmiss_IMnpim_woSidn_mc->Add(MMnmiss_IMnpim_dE_woK0_woSidn_Sigmam);
  MMnmiss_IMnpim_woSidn_mc->Add(MMnmiss_IMnpim_dE_woK0_woSidn_S0);
  //MMnmiss_IMnpim_woSidn_mc->Add(MMnmiss_IMnpim_dE_woK0_woSidn_Sppi0);
  //MMnmiss_IMnpim_woSidn_mc->Add(MMnmiss_IMnpim_dE_woK0_woSidn_Smpi0);

  TCanvas *cMMnmiss_IMnpim_woSidn = new TCanvas("cMMnmiss_IMnpim_woSidn","cMMnmiss_IMnpim_woSidn",1200,800);
  cMMnmiss_IMnpim_woSidn->Divide(2,1);
  cMMnmiss_IMnpim_woSidn->cd(1);
  MMnmiss_IMnpim_dE_woK0_woSidn_rdata->SetTitle("real data");
  MMnmiss_IMnpim_dE_woK0_woSidn_rdata->Draw("colz");
  cMMnmiss_IMnpim_woSidn->cd(2);
  MMnmiss_IMnpim_woSidn_mc->SetTitle("MC sum");
  MMnmiss_IMnpim_woSidn_mc->Draw("colz");
  
  TCanvas *cIMnpim_woSidn = new TCanvas("IMnpim_woSidn","IMnpim_woSidn");
  cIMnpim_woSidn->cd();
  MMnmiss_IMnpim_dE_woK0_woSidn_rdata->ProjectionX("MMnmiss_IMnpim_dE_woK0_woSidn_rdata_px")->Draw("HE");
  TH1D* IMnpim_woSidn_mc = (TH1D*)MMnmiss_IMnpim_woSidn_mc->ProjectionX("MMnmiss_IMnpim_woSidn_mc_px");
  IMnpim_woSidn_mc->SetLineColor(6);
  IMnpim_woSidn_mc->Draw("HEsame");
  TH1D* IMnpim_woSidn_lmiss = (TH1D*)MMnmiss_IMnpim_dE_woK0_woSidn_lmiss->ProjectionX("MMnmiss_IMnpim_dE_woK0_woSidn_lmiss_px");
  IMnpim_woSidn_lmiss->SetLineColor(4); 
  IMnpim_woSidn_lmiss->Draw("HEsame"); 
  TH1D* IMnpim_woSidn_S0 = (TH1D*)MMnmiss_IMnpim_dE_woK0_woSidn_S0->ProjectionX("MMnmiss_IMnpim_dE_woK0_woSidn_S0_px");
  IMnpim_woSidn_S0->SetLineColor(5); 
  IMnpim_woSidn_S0->Draw("HEsame");; 
  TH1D* IMnpim_woSidn_Sigmap = (TH1D*)MMnmiss_IMnpim_dE_woK0_woSidn_Sigmap->ProjectionX("MMnmiss_IMnpim_dE_woK0_woSidn_Sigmap_px");
  IMnpim_woSidn_Sigmap->SetLineColor(2); 
  IMnpim_woSidn_Sigmap->Draw("HEsame"); 
  TH1D* IMnpim_woSidn_Sigmam = (TH1D*)MMnmiss_IMnpip_dE_woK0_woSidn_Sigmam->ProjectionX("MMnmiss_IMnpim_dE_woK0_woSidn_Sigmam_px");
  IMnpim_woSidn_Sigmam->SetLineColor(3); 
  IMnpim_woSidn_Sigmam->Draw("HEsame"); 


  //weighting and decomposition of BG
  int nbinx_lmiss = MMnmiss_IMnpip_dE_woK0_woSidn_lmiss->GetNbinsX();
  int nbiny_lmiss = MMnmiss_IMnpip_dE_woK0_woSidn_lmiss->GetNbinsY();
  int nbinx_S0 = MMnmiss_IMnpip_dE_woK0_woSidn_S0->GetNbinsX();
  int nbiny_S0 = MMnmiss_IMnpip_dE_woK0_woSidn_S0->GetNbinsY();
  TH2D* MMnmiss_IMnpip_dE_woK0_woSidn_lmiss_mod = (TH1D*)MMnmiss_IMnpip_dE_woK0_woSidn_lmiss->Clone("MMnmiss_IMnpip_dE_woK0_woSidn_lmiss_mod");
  TH2D* MMnmiss_IMnpim_dE_woK0_woSidn_lmiss_mod = (TH1D*)MMnmiss_IMnpim_dE_woK0_woSidn_lmiss->Clone("MMnmiss_IMnpim_dE_woK0_woSidn_lmiss_mod");
  TH2D* MMnmiss_IMnpip_dE_woK0_woSidn_S0_mod = (TH1D*)MMnmiss_IMnpip_dE_woK0_woSidn_S0->Clone("MMnmiss_IMnpip_dE_woK0_woSidn_S0_mod");
  TH2D* MMnmiss_IMnpim_dE_woK0_woSidn_S0_mod = (TH1D*)MMnmiss_IMnpim_dE_woK0_woSidn_S0->Clone("MMnmiss_IMnpim_dE_woK0_woSidn_S0_mod");
  
  for(int ix=0;ix<nbinx_lmiss;ix++){
    for(int iy=0;iy<nbiny_lmiss;iy++){
      double cont = MMnmiss_IMnpip_dE_woK0_woSidn_lmiss_mod->GetBinContent(ix,iy);
      double weight = 
      
    }
  }



}
