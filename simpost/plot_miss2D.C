void plot_miss2D()
{
  TH1::SetDefaultSumw2();
                                //rdata,nSp,nSm,Lambda,S0,nSppi0,nSmpi0,mcsum
  const unsigned int colordef[8]={1,2,3,4,5,7,9,6};

  //real data
  TFile *filerdata = TFile::Open("../post/evanaIMpisigma_npippim_v196_out.root");
  //Sp
  TH2D* MMnmiss_IMnpipi_woK0_wSid_Sp_rdata = (TH2D*) filerdata->Get("MMnmiss_IMnpipi_woK0_wSid_Sp");
  TH1D* MMnmiss_woK0_wSid_Sp_rdata = (TH1D*)MMnmiss_IMnpipi_woK0_wSid_Sp_rdata->ProjectionY("MMnmiss_woK0_wSid_Sp_rdata");
  double  nrdata_Sp = MMnmiss_IMnpipi_woK0_wSid_Sp_rdata->Integral();
  TH1D* IMnpipi_woK0_wSid_Sp_rdata = (TH1D*)MMnmiss_IMnpipi_woK0_wSid_Sp_rdata->ProjectionX("IMnpipi_woK0_wSid_Sp_rdata");
  IMnpipi_woK0_wSid_Sp_rdata->SetLineColor(1);
  //Sm
  TH2D* MMnmiss_IMnpipi_woK0_wSid_Sm_rdata = (TH2D*) filerdata->Get("MMnmiss_IMnpipi_woK0_wSid_Sm");
  TH1D* MMnmiss_woK0_wSid_Sm_rdata = (TH1D*) MMnmiss_IMnpipi_woK0_wSid_Sm_rdata->ProjectionY("MMnmiss_woK0_wSid_Sm_rdata");
  double nrdata_Sm = MMnmiss_IMnpipi_woK0_wSid_Sm_rdata->Integral();
  TH1D* IMnpipi_woK0_wSid_Sm_rdata = (TH1D*)MMnmiss_IMnpipi_woK0_wSid_Sm_rdata->ProjectionX("IMnpipi_woK0_wSid_Sm_rdata");
  IMnpipi_woK0_wSid_Sm_rdata->SetLineColor(1);
  
  
  //Lambda pi+ pi- n sim.
  TFile *fileLambda = TFile::Open("simIMpisigma_npipiL_pippimn_v23_out.root");
  //Sp
  TH2D* MMnmiss_IMnpipi_woK0_wSid_Sp_Lambda = (TH2D*) fileLambda->Get("MMnmiss_IMnpipi_woK0_wSid_Sp");
  double nLambda_Sp = MMnmiss_IMnpipi_woK0_wSid_Sp_Lambda->Integral();
  MMnmiss_IMnpipi_woK0_wSid_Sp_Lambda->Scale(0.40*nrdata_Sp/nLambda_Sp);
  TH1D* MMnmiss_woK0_wSid_Sp_Lambda = (TH1D*)MMnmiss_IMnpipi_woK0_wSid_Sp_Lambda->ProjectionY("MMnmiss_woK0_wSid_Sp_Lambda");
  MMnmiss_woK0_wSid_Sp_Lambda->SetLineColor(4);
  TH1D* IMnpipi_woK0_wSid_Sp_Lambda = (TH1D*)MMnmiss_IMnpipi_woK0_wSid_Sp_Lambda->ProjectionX("IMnpipi_woK0_wSid_Sp_Lambda");
  IMnpipi_woK0_wSid_Sp_Lambda->SetLineColor(4);
  //Sm
  TH2D* MMnmiss_IMnpipi_woK0_wSid_Sm_Lambda = (TH2D*) fileLambda->Get("MMnmiss_IMnpipi_woK0_wSid_Sm");
  double nLambda_Sm = MMnmiss_IMnpipi_woK0_wSid_Sm_Lambda->Integral();
  MMnmiss_IMnpipi_woK0_wSid_Sm_Lambda->Scale(0.40*nrdata_Sm/nLambda_Sm);
  TH1D* MMnmiss_woK0_wSid_Sm_Lambda =  (TH1D*)MMnmiss_IMnpipi_woK0_wSid_Sm_Lambda->ProjectionY("MMnmiss_woK0_wSid_Sm_Lambda");
  MMnmiss_woK0_wSid_Sm_Lambda->SetLineColor(4);
  TH1D* IMnpipi_woK0_wSid_Sm_Lambda = (TH1D*)MMnmiss_IMnpipi_woK0_wSid_Sm_Lambda->ProjectionX("IMnpipi_woK0_wSid_Sm_Lambda");
  IMnpipi_woK0_wSid_Sm_Lambda->SetLineColor(4);
  

  //Sigma+
  TFile *filenSp = TFile::Open("simIMpisigma_nSppim_pippimn_v108_out.root");
  //Sp
  TH2D* MMnmiss_IMnpipi_woK0_wSid_Sp_Sigmap = (TH2D*) filenSp->Get("MMnmiss_IMnpipi_woK0_wSid_Sp");
  double nSigmap_Sp = MMnmiss_IMnpipi_woK0_wSid_Sp_Sigmap->Integral();
  MMnmiss_IMnpipi_woK0_wSid_Sp_Sigmap->Scale(0.18*nrdata_Sp/nSigmap_Sp);
  TH1D* MMnmiss_woK0_wSid_Sp_Sigmap = (TH1D*)MMnmiss_IMnpipi_woK0_wSid_Sp_Sigmap->ProjectionY("MMnmiss_woK0_wSid_Sp_Sigmap");
  MMnmiss_woK0_wSid_Sp_Sigmap->SetLineColor(2);
  TH1D* IMnpipi_woK0_wSid_Sp_Sigmap = (TH1D*)MMnmiss_IMnpipi_woK0_wSid_Sp_Sigmap->ProjectionX("IMnpipi_woK0_wSid_Sp_Sigmap");
  IMnpipi_woK0_wSid_Sp_Sigmap->SetLineColor(2);
  //Sm
  TH2D* MMnmiss_IMnpipi_woK0_wSid_Sm_Sigmap = (TH2D*) filenSp->Get("MMnmiss_IMnpipi_woK0_wSid_Sm");
  double nSigmap_Sm = MMnmiss_IMnpipi_woK0_wSid_Sm_Sigmap->Integral();
  MMnmiss_IMnpipi_woK0_wSid_Sm_Sigmap->Scale(0.01*nrdata_Sm/nSigmap_Sm);
  TH1D* MMnmiss_woK0_wSid_Sm_Sigmap = (TH1D*)MMnmiss_IMnpipi_woK0_wSid_Sm_Sigmap->ProjectionY("MMnmiss_woK0_wSid_Sm_Sigmap");
  MMnmiss_woK0_wSid_Sm_Sigmap->SetLineColor(2);
  TH1D* IMnpipi_woK0_wSid_Sm_Sigmap = (TH1D*)MMnmiss_IMnpipi_woK0_wSid_Sm_Sigmap->ProjectionX("IMnpipi_woK0_wSid_Sm_Sigmap");
  IMnpipi_woK0_wSid_Sm_Sigmap->SetLineColor(2);


  //Sigma-
  TFile *filenSm = TFile::Open("simIMpisigma_nSmpip_pippimn_v108_out.root");
  //Sp
  TH2D* MMnmiss_IMnpipi_woK0_wSid_Sp_Sigmam = (TH2D*) filenSm->Get("MMnmiss_IMnpipi_woK0_wSid_Sp");
  double nSigmam_Sp = MMnmiss_IMnpipi_woK0_wSid_Sp_Sigmam->Integral();
  MMnmiss_IMnpipi_woK0_wSid_Sp_Sigmam->Scale(0.01*nrdata_Sp/nSigmam_Sp); 
  TH1D* MMnmiss_woK0_wSid_Sp_Sigmam =  (TH1D*)MMnmiss_IMnpipi_woK0_wSid_Sp_Sigmam->ProjectionY("MMnmiss_woK0_wSid_Sp_Sigmam");
  MMnmiss_woK0_wSid_Sp_Sigmam->SetLineColor(3);
  TH1D* IMnpipi_woK0_wSid_Sp_Sigmam = (TH1D*)MMnmiss_IMnpipi_woK0_wSid_Sp_Sigmam->ProjectionX("IMnpipi_woK0_wSid_Sp_Sigmam");
  IMnpipi_woK0_wSid_Sp_Sigmam->SetLineColor(3);
  //Sm
  TH2D* MMnmiss_IMnpipi_woK0_wSid_Sm_Sigmam = (TH2D*) filenSm->Get("MMnmiss_IMnpipi_woK0_wSid_Sm");
  double nSigmam_Sm = MMnmiss_IMnpipi_woK0_wSid_Sm_Sigmam->Integral();
  MMnmiss_IMnpipi_woK0_wSid_Sm_Sigmam->Scale(0.21*nrdata_Sm/nSigmam_Sm);
  TH1D* MMnmiss_woK0_wSid_Sm_Sigmam =  (TH1D*) MMnmiss_IMnpipi_woK0_wSid_Sm_Sigmam->ProjectionY("MMnmiss_woK0_wSid_Sm_Sigmam");
  MMnmiss_woK0_wSid_Sm_Sigmam->SetLineColor(3);
  TH1D* IMnpipi_woK0_wSid_Sm_Sigmam = (TH1D*)MMnmiss_IMnpipi_woK0_wSid_Sm_Sigmam->ProjectionX("IMnpipi_woK0_wSid_Sm_Sigmam");
  IMnpipi_woK0_wSid_Sm_Sigmam->SetLineColor(3);
  

  //nS0pi+pi-
  TFile *fileS0pipi = TFile::Open("simIMpisigma_nS0pippim_pippimn_v5_out.root");
  //Sp
  TH2D* MMnmiss_IMnpipi_woK0_wSid_Sp_S0 = (TH2D*) fileS0pipi->Get("MMnmiss_IMnpipi_woK0_wSid_Sp");
  double nS0miss_Sp = MMnmiss_IMnpipi_woK0_wSid_Sp_S0->Integral();
  MMnmiss_IMnpipi_woK0_wSid_Sp_S0->Scale(0.24*nrdata_Sp/nS0miss_Sp);
  TH1D* MMnmiss_woK0_wSid_Sp_S0 =  (TH1D*) MMnmiss_IMnpipi_woK0_wSid_Sp_S0->ProjectionY("MMnmiss_woK0_wSid_Sp_S0");
  MMnmiss_woK0_wSid_Sp_S0->SetLineColor(5);
  TH1D* IMnpipi_woK0_wSid_Sp_S0 = (TH1D*)MMnmiss_IMnpipi_woK0_wSid_Sp_S0->ProjectionX("IMnpipi_woK0_wSid_Sp_S0");
  IMnpipi_woK0_wSid_Sp_S0->SetLineColor(5);
  //Sm
  TH2D* MMnmiss_IMnpipi_woK0_wSid_Sm_S0 = (TH2D*) fileS0pipi->Get("MMnmiss_IMnpipi_woK0_wSid_Sm");
  double nS0miss_Sm = MMnmiss_IMnpipi_woK0_wSid_Sm_S0->Integral();
  MMnmiss_IMnpipi_woK0_wSid_Sm_S0->Scale(0.2*nrdata_Sm/nS0miss_Sm);
  TH1D* MMnmiss_woK0_wSid_Sm_S0 = (TH1D*)MMnmiss_IMnpipi_woK0_wSid_Sm_S0->ProjectionY("MMnmiss_woK0_wSid_Sm_S0");
  MMnmiss_woK0_wSid_Sm_S0->SetLineColor(5);
  TH1D* IMnpipi_woK0_wSid_Sm_S0 = (TH1D*)MMnmiss_IMnpipi_woK0_wSid_Sm_S0->ProjectionX("IMnpipi_woK0_wSid_Sm_S0");
  IMnpipi_woK0_wSid_Sm_S0->SetLineColor(5);

  
  //nSppi-pi0
  TFile *fileSppi0 = TFile::Open("simIMpisigma_nSppimpi0_pippimn_v4_out.root");
  //Sp
  TH2D* MMnmiss_IMnpipi_woK0_wSid_Sp_Sppi0 = (TH2D*) fileSppi0->Get("MMnmiss_IMnpipi_woK0_wSid_Sp");
  double nSppi0_Sp= MMnmiss_IMnpipi_woK0_wSid_Sp_Sppi0->Integral();
  MMnmiss_IMnpipi_woK0_wSid_Sp_Sppi0->Scale(0.03*nrdata_Sp/nSppi0_Sp);
  TH1D* MMnmiss_woK0_wSid_Sp_Sppi0 = (TH1D*)MMnmiss_IMnpipi_woK0_wSid_Sp_Sppi0->ProjectionY("MMnmiss_woK0_wSid_Sp_Sppi0");
  MMnmiss_woK0_wSid_Sp_Sppi0->SetLineColor(7);
  TH1D* IMnpipi_woK0_wSid_Sp_Sppi0 = (TH1D*)MMnmiss_IMnpipi_woK0_wSid_Sp_Sppi0->ProjectionX("IMnpipi_woK0_wSid_Sp_Sppi0");
  IMnpipi_woK0_wSid_Sp_Sppi0->SetLineColor(7);
  //Sm
  TH2D* MMnmiss_IMnpipi_woK0_wSid_Sm_Sppi0 = (TH2D*) fileSppi0->Get("MMnmiss_IMnpipi_woK0_wSid_Sm");
  double nSppi0_Sm = MMnmiss_IMnpipi_woK0_wSid_Sm_Sppi0->Integral();
  MMnmiss_IMnpipi_woK0_wSid_Sm_Sppi0->Scale(0.03*nrdata_Sm/nSppi0_Sm);
  TH1D* MMnmiss_woK0_wSid_Sm_Sppi0 = (TH1D*)MMnmiss_IMnpipi_woK0_wSid_Sm_Sppi0->ProjectionY("MMnmiss_woK0_wSid_Sm_Sppi0");
  MMnmiss_woK0_wSid_Sm_Sppi0->SetLineColor(7);
  TH1D* IMnpipi_woK0_wSid_Sm_Sppi0 = (TH1D*)MMnmiss_IMnpipi_woK0_wSid_Sm_Sppi0->ProjectionX("IMnpipi_woK0_wSid_Sm_Sppi0");
  IMnpipi_woK0_wSid_Sm_Sppi0->SetLineColor(7);
  
  //nSmpi+pi0
  TFile *fileSmpi0 = TFile::Open("simIMpisigma_nSmpippi0_pippimn_v4_out.root");
  //Sp
  TH2D* MMnmiss_IMnpipi_woK0_wSid_Sp_Smpi0 = (TH2D*) fileSmpi0->Get("MMnmiss_IMnpipi_woK0_wSid_Sp");
  double nSmpi0_Sp = MMnmiss_IMnpipi_woK0_wSid_Sp_Smpi0->Integral();
  MMnmiss_IMnpipi_woK0_wSid_Sp_Smpi0->Scale(0.03*nrdata_Sp/nSmpi0_Sp);
  TH1D* MMnmiss_woK0_wSid_Sp_Smpi0 = (TH1D*)MMnmiss_IMnpipi_woK0_wSid_Sp_Smpi0->ProjectionY("MMnmiss_woK0_wSid_Sp_Smpi0");
  MMnmiss_woK0_wSid_Sp_Smpi0->SetLineColor(9);
  TH1D* IMnpipi_woK0_wSid_Sp_Smpi0 = (TH1D*)MMnmiss_IMnpipi_woK0_wSid_Sp_Smpi0->ProjectionX("IMnpipi_woK0_wSid_Sp_Smpi0");
  IMnpipi_woK0_wSid_Sp_Smpi0->SetLineColor(9);
  //Sm
  TH2D* MMnmiss_IMnpipi_woK0_wSid_Sm_Smpi0 = (TH2D*) fileSmpi0->Get("MMnmiss_IMnpipi_woK0_wSid_Sm");
  double nSmpi0_Sm = MMnmiss_IMnpipi_woK0_wSid_Sm_Smpi0->Integral();
  MMnmiss_IMnpipi_woK0_wSid_Sm_Smpi0->Scale(0.03*nrdata_Sm/nSmpi0_Sm);
  TH1D* MMnmiss_woK0_wSid_Sm_Smpi0 = (TH1D*)MMnmiss_IMnpipi_woK0_wSid_Sm_Smpi0->ProjectionY("MMnmiss_woK0_wSid_Sm_Smpi0");
  MMnmiss_woK0_wSid_Sm_Smpi0->SetLineColor(9);
  TH1D* IMnpipi_woK0_wSid_Sm_Smpi0 = (TH1D*)MMnmiss_IMnpipi_woK0_wSid_Sm_Smpi0->ProjectionX("IMnpipi_woK0_wSid_Sm_Smpi0");
  IMnpipi_woK0_wSid_Sm_Smpi0->SetLineColor(9);

  
  //missing mass and IMnpip/npim w/ Sp or Sm selection w/o K0
  TCanvas *miss_Sp = new TCanvas("miss_Sp","miss_Sp",800,800);
  miss_Sp->cd();
  
  TH1D* MMnmiss_woK0_wSid_Sp_mcsum = (TH1D*)MMnmiss_woK0_wSid_Sp_Lambda->Clone("MMnmiss_woK0_wSid_Sp_mcsum");
  MMnmiss_woK0_wSid_Sp_mcsum->Add(MMnmiss_woK0_wSid_Sp_Sigmap);
  MMnmiss_woK0_wSid_Sp_mcsum->Add(MMnmiss_woK0_wSid_Sp_Sigmam);
  MMnmiss_woK0_wSid_Sp_mcsum->Add(MMnmiss_woK0_wSid_Sp_S0);
  MMnmiss_woK0_wSid_Sp_mcsum->Add(MMnmiss_woK0_wSid_Sp_Sppi0);
  MMnmiss_woK0_wSid_Sp_mcsum->Add(MMnmiss_woK0_wSid_Sp_Smpi0);
  MMnmiss_woK0_wSid_Sp_mcsum->SetLineColor(6);
  //TFile *filerdata = TFile::Open("../post/evanaIMpisigma_npippim_v156_out.root");
  MMnmiss_woK0_wSid_Sp_rdata->SetLineColor(1);
  MMnmiss_woK0_wSid_Sp_rdata->Draw("HE");
  MMnmiss_woK0_wSid_Sp_mcsum->Draw("HEsame");
  MMnmiss_woK0_wSid_Sp_Sigmap->Draw("HEsame");
  MMnmiss_woK0_wSid_Sp_Sigmam->Draw("HEsame");
  MMnmiss_woK0_wSid_Sp_Lambda->Draw("HEsame");
  MMnmiss_woK0_wSid_Sp_S0->Draw("HEsame");
  MMnmiss_woK0_wSid_Sp_Sppi0->Draw("HEsame");
  MMnmiss_woK0_wSid_Sp_Smpi0->Draw("HEsame");
  
  TCanvas *miss_Sm = new TCanvas("miss_Sm","miss_Sm",800,0,800,800);
  miss_Sm->cd();
  TH1D* MMnmiss_woK0_wSid_Sm_mcsum = (TH1D*)MMnmiss_woK0_wSid_Sm_Lambda->Clone("MMnmiss_woK0_wSid_Sm_mcsum");
  MMnmiss_woK0_wSid_Sm_mcsum->Add(MMnmiss_woK0_wSid_Sm_Sigmap);
  MMnmiss_woK0_wSid_Sm_mcsum->Add(MMnmiss_woK0_wSid_Sm_Sigmam);
  MMnmiss_woK0_wSid_Sm_mcsum->Add(MMnmiss_woK0_wSid_Sm_S0);
  MMnmiss_woK0_wSid_Sm_mcsum->Add(MMnmiss_woK0_wSid_Sm_Sppi0);
  MMnmiss_woK0_wSid_Sm_mcsum->Add(MMnmiss_woK0_wSid_Sm_Smpi0);
  MMnmiss_woK0_wSid_Sm_mcsum->SetLineColor(6);
  //TFile *filerdata = TFile::Open("../post/evanaIMpisigma_npippim_v156_out.root");
  MMnmiss_woK0_wSid_Sm_rdata->SetLineColor(1);
  MMnmiss_woK0_wSid_Sm_rdata->Draw("HE");
  MMnmiss_woK0_wSid_Sm_mcsum->Draw("HEsame");
  MMnmiss_woK0_wSid_Sm_Sigmap->Draw("HEsame");
  MMnmiss_woK0_wSid_Sm_Sigmam->Draw("HEsame");
  MMnmiss_woK0_wSid_Sm_Lambda->Draw("HEsame");
  MMnmiss_woK0_wSid_Sm_S0->Draw("HEsame");
  MMnmiss_woK0_wSid_Sm_Sppi0->Draw("HEsame");
  MMnmiss_woK0_wSid_Sm_Smpi0->Draw("HEsame");
  
  TCanvas *cIMnpipi_woK0_wSid_Sp  = new TCanvas("cIMnpipi_woK0_wSid_Sp","cIMnpipi_woK0_wSid_Sp");
  cIMnpipi_woK0_wSid_Sp->cd();
  IMnpipi_woK0_wSid_Sp_rdata->Draw("HE");
  TH1D* IMnpipi_woK0_wSid_Sp_mcsum = IMnpipi_woK0_wSid_Sp_Lambda->Clone("IMnpipi_woK0_wSid_Sp_mcsum");
  IMnpipi_woK0_wSid_Sp_mcsum->Add(IMnpipi_woK0_wSid_Sp_S0);
  IMnpipi_woK0_wSid_Sp_mcsum->Add(IMnpipi_woK0_wSid_Sp_Sigmap);
  IMnpipi_woK0_wSid_Sp_mcsum->Add(IMnpipi_woK0_wSid_Sp_Sigmam);
  IMnpipi_woK0_wSid_Sp_mcsum->Add(IMnpipi_woK0_wSid_Sp_Sppi0);
  IMnpipi_woK0_wSid_Sp_mcsum->Add(IMnpipi_woK0_wSid_Sp_Smpi0);
  IMnpipi_woK0_wSid_Sp_mcsum->SetLineColor(6);
  IMnpipi_woK0_wSid_Sp_mcsum->Draw("HEsame");
  IMnpipi_woK0_wSid_Sp_Lambda->Draw("HEsame");
  IMnpipi_woK0_wSid_Sp_S0->Draw("HEsame");
  IMnpipi_woK0_wSid_Sp_Sigmap->Draw("HEsame");
  IMnpipi_woK0_wSid_Sp_Sigmam->Draw("HEsame");
  IMnpipi_woK0_wSid_Sp_Sppi0->Draw("HEsame");
  IMnpipi_woK0_wSid_Sp_Smpi0->Draw("HEsame");

  TCanvas *cIMnpipi_woK0_wSid_Sm  = new TCanvas("cIMnpipi_woK0_wSid_Sm","cIMnpipi_woK0_wSid_Sm");
  cIMnpipi_woK0_wSid_Sm->cd();
  IMnpipi_woK0_wSid_Sm_rdata->Draw("HE");
  TH1D* IMnpipi_woK0_wSid_Sm_mcsum = IMnpipi_woK0_wSid_Sm_Lambda->Clone("IMnpipi_woK0_wSid_Sm_mcsum");
  IMnpipi_woK0_wSid_Sm_mcsum->Add(IMnpipi_woK0_wSid_Sm_S0);
  IMnpipi_woK0_wSid_Sm_mcsum->Add(IMnpipi_woK0_wSid_Sm_Sigmap);
  IMnpipi_woK0_wSid_Sm_mcsum->Add(IMnpipi_woK0_wSid_Sm_Sigmam);
  IMnpipi_woK0_wSid_Sm_mcsum->Add(IMnpipi_woK0_wSid_Sm_Sppi0);
  IMnpipi_woK0_wSid_Sm_mcsum->Add(IMnpipi_woK0_wSid_Sm_Smpi0);
  IMnpipi_woK0_wSid_Sm_mcsum->SetLineColor(6);
  IMnpipi_woK0_wSid_Sm_mcsum->Draw("HEsame");
  IMnpipi_woK0_wSid_Sm_Lambda->Draw("HEsame");
  IMnpipi_woK0_wSid_Sm_S0->Draw("HEsame");
  IMnpipi_woK0_wSid_Sm_Sigmap->Draw("HEsame");
  IMnpipi_woK0_wSid_Sm_Sigmam->Draw("HEsame");
  IMnpipi_woK0_wSid_Sm_Sppi0->Draw("HEsame");
  IMnpipi_woK0_wSid_Sm_Smpi0->Draw("HEsame");

  //superimpose 2d hists of miss mass vs IM(npi+) or IM(npi-)
  //Lambda
  TH2D *MMnmiss_IMnpip_dE_woK0_Lambda = (TH2D*)fileLambda->Get("MMnmiss_IMnpip_dE_woK0");
  TH2D *MMnmiss_IMnpim_dE_woK0_Lambda = (TH2D*)fileLambda->Get("MMnmiss_IMnpim_dE_woK0");
  TH2D *MMnmiss_IMnpip_dE_woK0_woSm_Lambda = (TH2D*)fileLambda->Get("MMnmiss_IMnpip_dE_woK0_woSm");
  TH2D *MMnmiss_IMnpim_dE_woK0_woSp_Lambda = (TH2D*)fileLambda->Get("MMnmiss_IMnpim_dE_woK0_woSp");
  TH2D *MMnmiss_IMnpip_dE_woK0_woSm_n_Lambda = (TH2D*)fileLambda->Get("MMnmiss_IMnpip_dE_woK0_woSm_n");
  TH2D *MMnmiss_IMnpim_dE_woK0_woSp_n_Lambda = (TH2D*)fileLambda->Get("MMnmiss_IMnpim_dE_woK0_woSp_n");
  TH2D *MMnmiss_IMnpip_dE_woK0_woSidn_Lambda = (TH2D*)fileLambda->Get("MMnmiss_IMnpip_dE_woK0_woSidn");
  TH2D *MMnmiss_IMnpim_dE_woK0_woSidn_Lambda = (TH2D*)fileLambda->Get("MMnmiss_IMnpim_dE_woK0_woSidn");
  TH2D *MMnmiss_IMpippim_dE_Lambda = (TH2D*)fileLambda->Get("MMnmiss_IMpippim_dE");
  TH2D *MMnmiss_IMpippim_dE_wSid_Lambda = (TH2D*)fileLambda->Get("MMnmiss_IMpippim_wSid_dE");
  TH2D *MMnmiss_IMpippim_dE_woK0_woSidn_Lambda = (TH2D*)fileLambda->Get("MMnmiss_IMpippim_dE_woK0_woSidn");
  TH2D *IMnpim_IMnpip_dE_Lambda = (TH2D*)fileLambda->Get("IMnpim_IMnpip_dE");
  TH2D *IMnpim_IMnpip_dE_woK0_woSidn_Lambda = (TH2D*)fileLambda->Get("IMnpim_IMnpip_dE_woK0_woSidn");
  double Nnpip_Lambda = MMnmiss_IMnpip_dE_woK0_Lambda->Integral();
  double Nnpim_Lambda = MMnmiss_IMnpim_dE_woK0_Lambda->Integral();
  //nSppi-
  TH2D *MMnmiss_IMnpip_dE_woK0_Sigmap = (TH2D*)filenSp->Get("MMnmiss_IMnpip_dE_woK0");
  TH2D *MMnmiss_IMnpim_dE_woK0_Sigmap = (TH2D*)filenSp->Get("MMnmiss_IMnpim_dE_woK0");
  TH2D *MMnmiss_IMnpip_dE_woK0_woSm_Sigmap = (TH2D*)filenSp->Get("MMnmiss_IMnpip_dE_woK0_woSm");
  TH2D *MMnmiss_IMnpim_dE_woK0_woSp_Sigmap = (TH2D*)filenSp->Get("MMnmiss_IMnpim_dE_woK0_woSp");
  TH2D *MMnmiss_IMnpip_dE_woK0_woSm_n_Sigmap = (TH2D*)filenSp->Get("MMnmiss_IMnpip_dE_woK0_woSm_n");
  TH2D *MMnmiss_IMnpim_dE_woK0_woSp_n_Sigmap = (TH2D*)filenSp->Get("MMnmiss_IMnpim_dE_woK0_woSp_n");
  TH2D *MMnmiss_IMnpip_dE_woK0_woSidn_Sigmap = (TH2D*)filenSp->Get("MMnmiss_IMnpip_dE_woK0_woSidn");
  TH2D *MMnmiss_IMnpim_dE_woK0_woSidn_Sigmap = (TH2D*)filenSp->Get("MMnmiss_IMnpim_dE_woK0_woSidn");
  TH2D *MMnmiss_IMpippim_dE_Sigmap = (TH2D*)filenSp->Get("MMnmiss_IMpippim_dE");
  TH2D *MMnmiss_IMpippim_dE_wSid_Sigmap = (TH2D*)filenSp->Get("MMnmiss_IMpippim_wSid_dE");
  TH2D *MMnmiss_IMpippim_dE_woK0_woSidn_Sigmap = (TH2D*)filenSp->Get("MMnmiss_IMpippim_dE_woK0_woSidn");
  TH2D *IMnpim_IMnpip_dE_Sigmap = (TH2D*)filenSp->Get("IMnpim_IMnpip_dE");
  TH2D *IMnpim_IMnpip_dE_woK0_woSidn_Sigmap = (TH2D*)filenSp->Get("IMnpim_IMnpip_dE_woK0_woSidn");
  double Nnpip_Sigmap = MMnmiss_IMnpip_dE_woK0_Sigmap->Integral();
  double Nnpim_Sigmap = MMnmiss_IMnpim_dE_woK0_Sigmap->Integral();
  //nSmpi+
  TH2D *MMnmiss_IMnpip_dE_woK0_Sigmam = (TH2D*)filenSm->Get("MMnmiss_IMnpip_dE_woK0");
  TH2D *MMnmiss_IMnpim_dE_woK0_Sigmam = (TH2D*)filenSm->Get("MMnmiss_IMnpim_dE_woK0");
  TH2D *MMnmiss_IMnpip_dE_woK0_woSm_Sigmam = (TH2D*)filenSm->Get("MMnmiss_IMnpip_dE_woK0_woSm");
  TH2D *MMnmiss_IMnpim_dE_woK0_woSp_Sigmam = (TH2D*)filenSm->Get("MMnmiss_IMnpim_dE_woK0_woSp");
  TH2D *MMnmiss_IMnpip_dE_woK0_woSm_n_Sigmam = (TH2D*)filenSm->Get("MMnmiss_IMnpip_dE_woK0_woSm_n");
  TH2D *MMnmiss_IMnpim_dE_woK0_woSp_n_Sigmam = (TH2D*)filenSm->Get("MMnmiss_IMnpim_dE_woK0_woSp_n");
  TH2D *MMnmiss_IMnpip_dE_woK0_woSidn_Sigmam = (TH2D*)filenSm->Get("MMnmiss_IMnpip_dE_woK0_woSidn");
  TH2D *MMnmiss_IMnpim_dE_woK0_woSidn_Sigmam = (TH2D*)filenSm->Get("MMnmiss_IMnpim_dE_woK0_woSidn");
  TH2D *MMnmiss_IMpippim_dE_Sigmam = (TH2D*)filenSm->Get("MMnmiss_IMpippim_dE");
  TH2D *MMnmiss_IMpippim_dE_wSid_Sigmam = (TH2D*)filenSm->Get("MMnmiss_IMpippim_wSid_dE");
  TH2D *MMnmiss_IMpippim_dE_woK0_woSidn_Sigmam = (TH2D*)filenSm->Get("MMnmiss_IMpippim_dE_woK0_woSidn");
  TH2D *IMnpim_IMnpip_dE_Sigmam = (TH2D*)filenSm->Get("IMnpim_IMnpip_dE");
  TH2D *IMnpim_IMnpip_dE_woK0_woSidn_Sigmam = (TH2D*)filenSm->Get("IMnpim_IMnpip_dE_woK0_woSidn");
  double Nnpip_Sigmam = MMnmiss_IMnpip_dE_woK0_Sigmam->Integral();
  double Nnpim_Sigmam = MMnmiss_IMnpim_dE_woK0_Sigmam->Integral();
  //Sigma0
  TH2D *MMnmiss_IMnpip_dE_woK0_S0 = (TH2D*)fileS0pipi->Get("MMnmiss_IMnpip_dE_woK0");
  TH2D *MMnmiss_IMnpim_dE_woK0_S0 = (TH2D*)fileS0pipi->Get("MMnmiss_IMnpim_dE_woK0");
  TH2D *MMnmiss_IMnpip_dE_woK0_woSm_S0 = (TH2D*)fileS0pipi->Get("MMnmiss_IMnpip_dE_woK0_woSm");
  TH2D *MMnmiss_IMnpim_dE_woK0_woSp_S0 = (TH2D*)fileS0pipi->Get("MMnmiss_IMnpim_dE_woK0_woSp");
  TH2D *MMnmiss_IMnpip_dE_woK0_woSm_n_S0 = (TH2D*)fileS0pipi->Get("MMnmiss_IMnpip_dE_woK0_woSm_n");
  TH2D *MMnmiss_IMnpim_dE_woK0_woSp_n_S0 = (TH2D*)fileS0pipi->Get("MMnmiss_IMnpim_dE_woK0_woSp_n");
  TH2D *MMnmiss_IMnpip_dE_woK0_woSidn_S0 = (TH2D*)fileS0pipi->Get("MMnmiss_IMnpip_dE_woK0_woSidn");
  TH2D *MMnmiss_IMnpim_dE_woK0_woSidn_S0 = (TH2D*)fileS0pipi->Get("MMnmiss_IMnpim_dE_woK0_woSidn");
  TH2D *MMnmiss_IMpippim_dE_S0 = (TH2D*)fileS0pipi->Get("MMnmiss_IMpippim_dE");
  TH2D *MMnmiss_IMpippim_dE_wSid_S0 = (TH2D*)fileS0pipi->Get("MMnmiss_IMpippim_wSid_dE");
  TH2D *MMnmiss_IMpippim_dE_woK0_woSidn_S0 = (TH2D*)fileS0pipi->Get("MMnmiss_IMpippim_dE_woK0_woSidn");
  TH2D *IMnpim_IMnpip_dE_S0 = (TH2D*)fileS0pipi->Get("IMnpim_IMnpip_dE");
  TH2D *IMnpim_IMnpip_dE_woK0_woSidn_S0 = (TH2D*)fileS0pipi->Get("IMnpim_IMnpip_dE_woK0_woSidn");
  double Nnpip_S0 = MMnmiss_IMnpip_dE_woK0_S0->Integral();
  double Nnpim_S0 = MMnmiss_IMnpim_dE_woK0_S0->Integral();
  //Sppi0
  TH2D *MMnmiss_IMnpip_dE_woK0_Sppi0 = (TH2D*)fileSppi0->Get("MMnmiss_IMnpip_dE_woK0");
  TH2D *MMnmiss_IMnpim_dE_woK0_Sppi0 = (TH2D*)fileSppi0->Get("MMnmiss_IMnpim_dE_woK0");
  TH2D *MMnmiss_IMnpip_dE_woK0_woSm_Sppi0 = (TH2D*)fileSppi0->Get("MMnmiss_IMnpip_dE_woK0_woSm");
  TH2D *MMnmiss_IMnpim_dE_woK0_woSp_Sppi0 = (TH2D*)fileSppi0->Get("MMnmiss_IMnpim_dE_woK0_woSp");
  TH2D *MMnmiss_IMnpip_dE_woK0_woSm_n_Sppi0 = (TH2D*)fileSppi0->Get("MMnmiss_IMnpip_dE_woK0_woSm_n");
  TH2D *MMnmiss_IMnpim_dE_woK0_woSp_n_Sppi0 = (TH2D*)fileSppi0->Get("MMnmiss_IMnpim_dE_woK0_woSp_n");
  TH2D *MMnmiss_IMnpip_dE_woK0_woSidn_Sppi0 = (TH2D*)fileSppi0->Get("MMnmiss_IMnpip_dE_woK0_woSidn");
  TH2D *MMnmiss_IMnpim_dE_woK0_woSidn_Sppi0 = (TH2D*)fileSppi0->Get("MMnmiss_IMnpim_dE_woK0_woSidn");
  TH2D *MMnmiss_IMpippim_dE_Sppi0 = (TH2D*)fileSppi0->Get("MMnmiss_IMpippim_dE");
  TH2D *MMnmiss_IMpippim_dE_wSid_Sppi0 = (TH2D*)fileSppi0->Get("MMnmiss_IMpippim_wSid_dE");
  TH2D *MMnmiss_IMpippim_dE_woK0_woSidn_Sppi0 = (TH2D*)fileSppi0->Get("MMnmiss_IMpippim_dE_woK0_woSidn");
  TH2D *IMnpim_IMnpip_dE_Sppi0 = (TH2D*)fileSppi0->Get("IMnpim_IMnpip_dE");
  TH2D *IMnpim_IMnpip_dE_woK0_woSidn_Sppi0 = (TH2D*)fileSppi0->Get("IMnpim_IMnpip_dE_woK0_woSidn");
  double Nnpip_Sppi0 = MMnmiss_IMnpip_dE_woK0_Sppi0->Integral();
  double Nnpim_Sppi0 = MMnmiss_IMnpim_dE_woK0_Sppi0->Integral();
  //Smpi0
  TH2D *MMnmiss_IMnpip_dE_woK0_Smpi0 = (TH2D*)fileSmpi0->Get("MMnmiss_IMnpip_dE_woK0");
  TH2D *MMnmiss_IMnpim_dE_woK0_Smpi0 = (TH2D*)fileSmpi0->Get("MMnmiss_IMnpim_dE_woK0");
  TH2D *MMnmiss_IMnpip_dE_woK0_woSm_Smpi0 = (TH2D*)fileSmpi0->Get("MMnmiss_IMnpip_dE_woK0_woSm");
  TH2D *MMnmiss_IMnpim_dE_woK0_woSp_Smpi0 = (TH2D*)fileSmpi0->Get("MMnmiss_IMnpim_dE_woK0_woSp");
  TH2D *MMnmiss_IMnpip_dE_woK0_woSm_n_Smpi0 = (TH2D*)fileSmpi0->Get("MMnmiss_IMnpip_dE_woK0_woSm_n");
  TH2D *MMnmiss_IMnpim_dE_woK0_woSp_n_Smpi0 = (TH2D*)fileSmpi0->Get("MMnmiss_IMnpim_dE_woK0_woSp_n");
  TH2D *MMnmiss_IMnpip_dE_woK0_woSidn_Smpi0 = (TH2D*)fileSmpi0->Get("MMnmiss_IMnpip_dE_woK0_woSidn");
  TH2D *MMnmiss_IMnpim_dE_woK0_woSidn_Smpi0 = (TH2D*)fileSmpi0->Get("MMnmiss_IMnpim_dE_woK0_woSidn");
  TH2D *MMnmiss_IMpippim_dE_Smpi0 = (TH2D*)fileSmpi0->Get("MMnmiss_IMpippim_dE");
  TH2D *MMnmiss_IMpippim_dE_wSid_Smpi0 = (TH2D*)fileSmpi0->Get("MMnmiss_IMpippim_wSid_dE");
  TH2D *MMnmiss_IMpippim_dE_woK0_woSidn_Smpi0 = (TH2D*)fileSmpi0->Get("MMnmiss_IMpippim_dE_woK0_woSidn");
  TH2D *IMnpim_IMnpip_dE_Smpi0 = (TH2D*)fileSmpi0->Get("IMnpim_IMnpip_dE");
  TH2D *IMnpim_IMnpip_dE_woK0_woSidn_Smpi0 = (TH2D*)fileSmpi0->Get("IMnpim_IMnpip_dE_woK0_woSidn");
  double Nnpip_Smpi0 = MMnmiss_IMnpip_dE_woK0_Smpi0->Integral();
  double Nnpim_Smpi0 = MMnmiss_IMnpim_dE_woK0_Smpi0->Integral();
  //real data
  TH2D *MMnmiss_IMnpip_dE_woK0_rdata = (TH2D*)filerdata->Get("MMnmiss_IMnpip_dE_woK0");
  TH2D *MMnmiss_IMnpim_dE_woK0_rdata = (TH2D*)filerdata->Get("MMnmiss_IMnpim_dE_woK0");
  TH2D *MMnmiss_IMnpip_dE_woK0_woSm_rdata = (TH2D*)filerdata->Get("MMnmiss_IMnpip_dE_woK0_woSm");
  TH2D *MMnmiss_IMnpim_dE_woK0_woSp_rdata = (TH2D*)filerdata->Get("MMnmiss_IMnpim_dE_woK0_woSp");
  TH2D *MMnmiss_IMnpip_dE_woK0_woSm_n_rdata = (TH2D*)filerdata->Get("MMnmiss_IMnpip_dE_woK0_woSm_n");
  TH2D *MMnmiss_IMnpim_dE_woK0_woSp_n_rdata = (TH2D*)filerdata->Get("MMnmiss_IMnpim_dE_woK0_woSp_n");
  TH2D *MMnmiss_IMnpip_dE_woK0_woSidn_rdata = (TH2D*)filerdata->Get("MMnmiss_IMnpip_dE_woK0_woSidn");
  TH2D *MMnmiss_IMnpim_dE_woK0_woSidn_rdata = (TH2D*)filerdata->Get("MMnmiss_IMnpim_dE_woK0_woSidn");
  TH2D *MMnmiss_IMpippim_dE_rdata = (TH2D*)filerdata->Get("MMnmiss_IMpippim_dE");
  TH2D *MMnmiss_IMpippim_dE_wSid_rdata = (TH2D*)filerdata->Get("MMnmiss_IMpippim_wSid_dE");
  TH2D *MMnmiss_IMpippim_dE_woK0_woSidn_rdata = (TH2D*)filerdata->Get("MMnmiss_IMpippim_dE_woK0_woSidn");
  TH2D *IMnpim_IMnpip_dE_rdata = (TH2D*)filerdata->Get("IMnpim_IMnpip_dE");
  TH2D *IMnpim_IMnpip_dE_woK0_woSidn_rdata = (TH2D*)filerdata->Get("IMnpim_IMnpip_dE_woK0_woSidn");
  double Nnpip_rdata = MMnmiss_IMnpip_dE_woK0_rdata->Integral();
  double Nnpim_rdata = MMnmiss_IMnpim_dE_woK0_rdata->Integral();
  
  const double scaleFactor_Sp[7]={1.0,
                                  0.18*nrdata_Sp/nSigmap_Sp,
                                  0.01*nrdata_Sp/nSigmam_Sp,
                                  0.40*nrdata_Sp/nLambda_Sp,
                                  0.24*nrdata_Sp/nS0miss_Sp,
                                  0.03*nrdata_Sp/nSppi0_Sp,
                                  0.03*nrdata_Sp/nSmpi0_Sp};

  //array 
  TH2D* MMnmiss_IMnpip_dE_woK0[7]={
    MMnmiss_IMnpip_dE_woK0_rdata,
    MMnmiss_IMnpip_dE_woK0_Sigmap,
    MMnmiss_IMnpip_dE_woK0_Sigmam,
    MMnmiss_IMnpip_dE_woK0_Lambda,
    MMnmiss_IMnpip_dE_woK0_S0,
    MMnmiss_IMnpip_dE_woK0_Sppi0,
    MMnmiss_IMnpip_dE_woK0_Smpi0};

  TH2D* MMnmiss_IMnpip_dE_woK0_woSm[7]={
    MMnmiss_IMnpip_dE_woK0_woSm_rdata,
    MMnmiss_IMnpip_dE_woK0_woSm_Sigmap,
    MMnmiss_IMnpip_dE_woK0_woSm_Sigmam,
    MMnmiss_IMnpip_dE_woK0_woSm_Lambda,
    MMnmiss_IMnpip_dE_woK0_woSm_S0,
    MMnmiss_IMnpip_dE_woK0_woSm_Sppi0,
    MMnmiss_IMnpip_dE_woK0_woSm_Smpi0};
  
  TH2D* MMnmiss_IMnpip_dE_woK0_woSm_n[7]={
    MMnmiss_IMnpip_dE_woK0_woSm_n_rdata,
    MMnmiss_IMnpip_dE_woK0_woSm_n_Sigmap,
    MMnmiss_IMnpip_dE_woK0_woSm_n_Sigmam,
    MMnmiss_IMnpip_dE_woK0_woSm_n_Lambda,
    MMnmiss_IMnpip_dE_woK0_woSm_n_S0,
    MMnmiss_IMnpip_dE_woK0_woSm_n_Sppi0,
    MMnmiss_IMnpip_dE_woK0_woSm_n_Smpi0};

  TH2D* MMnmiss_IMnpip_dE_woK0_woSidn[7]={
    MMnmiss_IMnpip_dE_woK0_woSidn_rdata,
    MMnmiss_IMnpip_dE_woK0_woSidn_Sigmap,
    MMnmiss_IMnpip_dE_woK0_woSidn_Sigmam,
    MMnmiss_IMnpip_dE_woK0_woSidn_Lambda,
    MMnmiss_IMnpip_dE_woK0_woSidn_S0,
    MMnmiss_IMnpip_dE_woK0_woSidn_Sppi0,
    MMnmiss_IMnpip_dE_woK0_woSidn_Smpi0};

  
  for(int i=0;i<7;i++){
    MMnmiss_IMnpip_dE_woK0[i]->Scale(scaleFactor_Sp[i]);
    MMnmiss_IMnpip_dE_woK0_woSm[i]->Scale(scaleFactor_Sp[i]);
    MMnmiss_IMnpip_dE_woK0_woSm_n[i]->Scale(scaleFactor_Sp[i]);
    MMnmiss_IMnpip_dE_woK0_woSidn[i]->Scale(scaleFactor_Sp[i]);
  }

  TH2D *MMnmiss_IMnpip_mc = MMnmiss_IMnpip_dE_woK0[1]->Clone("MMnmiss_IMnpip_mc");
  for(int i=2;i<7;i++)  MMnmiss_IMnpip_mc->Add(MMnmiss_IMnpip_dE_woK0[i]);

  TCanvas *cMMnmiss_IMnpip = new TCanvas("cMMnmiss_IMnpip","cMMnmiss_IMnpip",1200,800);
  cMMnmiss_IMnpip->Divide(2,1);
  cMMnmiss_IMnpip->cd(1);
  MMnmiss_IMnpip_dE_woK0_rdata->SetTitle("real data");
  MMnmiss_IMnpip_dE_woK0_rdata->Draw("colz");
  cMMnmiss_IMnpip->cd(2);
  MMnmiss_IMnpip_mc->SetTitle("MC sum");
  MMnmiss_IMnpip_mc->Draw("colz");
  
  //
  TH2D *MMnmiss_IMnpip_woK0_woSm_mc = (TH2D*)MMnmiss_IMnpip_dE_woK0_woSm[2]->Clone("MMnmiss_IMnpip_woK0_woSm_mc");
  for(int i=3;i<7;i++)MMnmiss_IMnpip_woK0_woSm_mc->Add(MMnmiss_IMnpip_dE_woK0_woSm[i]);
  
  //
  TH2D *MMnmiss_IMnpip_woK0_woSm_n_mc = (TH2D*)MMnmiss_IMnpip_dE_woK0_woSm_n[2]->Clone("MMnmiss_IMnpip_woK0_woSm_n_mc");
  for(int i=3;i<7;i++)MMnmiss_IMnpip_woK0_woSm_n_mc->Add(MMnmiss_IMnpip_dE_woK0_woSm_n[i]);

  TCanvas *cMMnmiss_IMnpip_woK0_woSm = new TCanvas("cMMnmiss_IMnpip_woK0_woSm","cMMnmiss_IMnpip_woK0_woSm",1200,800);
  cMMnmiss_IMnpip_woK0_woSm->Divide(2,1);
  cMMnmiss_IMnpip_woK0_woSm->cd(1);
  MMnmiss_IMnpip_dE_woK0_woSm_rdata->SetTitle("real data");
  MMnmiss_IMnpip_dE_woK0_woSm_rdata->Draw("colz");
  cMMnmiss_IMnpip_woK0_woSm->cd(2);
  MMnmiss_IMnpip_woK0_woSm_mc->SetTitle("MC sum");
  MMnmiss_IMnpip_woK0_woSm_mc->Draw("colz");
  
  const char name[][10]={"rdata","Sigmap","Sigman","Lambda","S0","Sppi0","Smpi0"};
  TCanvas *cIMnpip_woK0_woSm_n = new TCanvas("cIMnpip_woK0_woSm_n","cIMnpip_woK0_woSm_n");
  cIMnpip_woK0_woSm_n->cd();
  TH1D *IMnpip_woK0_woSm_n[7];
  for(int i=0;i<7;i++){
    IMnpip_woK0_woSm_n[i] = (TH1D*)MMnmiss_IMnpip_dE_woK0_woSm_n[i]->ProjectionX(Form("IMnpip_woK0_woSm_n_%s",name[i]));
    IMnpip_woK0_woSm_n[i]->SetLineColor(colordef[i]);
  }
  TH1D* IMnpip_woK0_woSm_n_mc = IMnpip_woK0_woSm_n[2]->Clone("IMnpip_woK0_woSm_n_mc");
  for(int i=3;i<7;i++)IMnpip_woK0_woSm_n_mc->Add(IMnpip_woK0_woSm_n[i]);
  IMnpip_woK0_woSm_n[0]->Draw("HE");
  for(int i=2;i<7;i++)IMnpip_woK0_woSm_n[i]->Draw("HEsame");
  IMnpip_woK0_woSm_n_mc->SetLineColor(6);
  IMnpip_woK0_woSm_n_mc->Draw("HEsame");

  /*
  TH1D* IMnpip_woK0_woSm_n_rdata = MMnmiss_IMnpip_dE_woK0_woSm_n_rdata->ProjectionX("IMnpip_woK0_woSm_n_rdata");
  IMnpip_woK0_woSm_n_rdata->Draw("HE");
  TH1D* IMnpip_woK0_woSm_n_mc = MMnmiss_IMnpip_woK0_woSm_n_mc->ProjectionX("IMnpip_woK0_woSm_n_mc");
  IMnpip_woK0_woSm_n_mc->SetLineColor(6);
  IMnpip_woK0_woSm_n_mc->Draw("HEsame");
  TH1D* IMnpip_woK0_woSm_n[7]
  IMnpip_woK0_woSm_n_Sigmam->SetLineColor(3);
  IMnpip_woK0_woSm_n_S0->SetLineColor(4);
  IMnpip_woK0_woSm_n_Sppi0->SetLineColor(4);
  */

  //
  TH2D *MMnmiss_IMnpip_woSidn_mc = (TH2D*)MMnmiss_IMnpip_dE_woK0_woSidn_Lambda->Clone("MMnmiss_IMnpip_woSidn_mc");
  MMnmiss_IMnpip_woSidn_mc->Add(MMnmiss_IMnpip_dE_woK0_woSidn_Sigmap);
  MMnmiss_IMnpip_woSidn_mc->Add(MMnmiss_IMnpip_dE_woK0_woSidn_Sigmam);
  MMnmiss_IMnpip_woSidn_mc->Add(MMnmiss_IMnpip_dE_woK0_woSidn_S0);
  MMnmiss_IMnpip_woSidn_mc->Add(MMnmiss_IMnpip_dE_woK0_woSidn_Sppi0);
  MMnmiss_IMnpip_woSidn_mc->Add(MMnmiss_IMnpip_dE_woK0_woSidn_Smpi0);

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
  TH1D* MMnmiss_dE_woK0_woSidn_rdata = MMnmiss_IMnpip_dE_woK0_woSidn_rdata->ProjectionY("MMnmiss_dE_woK0_woSidn_rdata");
  MMnmiss_dE_woK0_woSidn_rdata->Draw("HE");
  TH1D* MMnmiss_woSidn_mc = (TH1D*)MMnmiss_IMnpip_woSidn_mc->ProjectionY("MMnmiss_IMnpip_woSidn_mc_py");
  MMnmiss_woSidn_mc->SetLineColor(6);
  MMnmiss_woSidn_mc->Draw("HEsame");

  TH1D* MMnmiss_woSidn_Lambda = (TH1D*)MMnmiss_IMnpip_dE_woK0_woSidn_Lambda->ProjectionY("MMnmiss_IMnpip_dE_woK0_woSidn_Lambda_py");
  MMnmiss_woSidn_Lambda->SetLineColor(4); 
  MMnmiss_woSidn_Lambda->Draw("HEsame"); 
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
  MMnmiss_woSidn_ratio->Divide(MMnmiss_dE_woK0_woSidn_rdata);
  MMnmiss_woSidn_ratio->Draw("HE");

  
  //projection to IMnpip (miss n & Sigma+/-)
  TCanvas *cIMnpip_woSidn = new TCanvas("cIMnpip_woSidn","cIMnpip_woSidn");
  cIMnpip_woSidn->cd();
  TH1D* IMnpip_woSidn_rdata = (TH1D*) MMnmiss_IMnpip_dE_woK0_woSidn_rdata->ProjectionX("IMnpip_dE_woK0_woSidn_rdata");
  IMnpip_woSidn_rdata->Draw("HE");
  TH1D* IMnpip_woSidn_mc = (TH1D*)MMnmiss_IMnpip_woSidn_mc->ProjectionX("MMnmiss_IMnpip_woSidn_mc_px");
  IMnpip_woSidn_mc->SetLineColor(6);
  IMnpip_woSidn_mc->Draw("HEsame");
  TH1D* IMnpip_woSidn_Lambda = (TH1D*)MMnmiss_IMnpip_dE_woK0_woSidn_Lambda->ProjectionX("MMnmiss_IMnpip_dE_woK0_woSidn_Lambda_px");
  IMnpip_woSidn_Lambda->SetLineColor(4); 
  IMnpip_woSidn_Lambda->Draw("HEsame"); 
  TH1D* IMnpip_woSidn_S0 = (TH1D*)MMnmiss_IMnpip_dE_woK0_woSidn_S0->ProjectionX("MMnmiss_IMnpip_dE_woK0_woSidn_S0_px");
  IMnpip_woSidn_S0->SetLineColor(5); 
  IMnpip_woSidn_S0->Draw("HEsame");; 
  TH1D* IMnpip_woSidn_Sigmap = (TH1D*)MMnmiss_IMnpip_dE_woK0_woSidn_Sigmap->ProjectionX("MMnmiss_IMnpip_dE_woK0_woSidn_Sigmap_px");
  IMnpip_woSidn_Sigmap->SetLineColor(2); 
  IMnpip_woSidn_Sigmap->Draw("HEsame"); 
  TH1D* IMnpip_woSidn_Sigmam = (TH1D*)MMnmiss_IMnpip_dE_woK0_woSidn_Sigmam->ProjectionX("MMnmiss_IMnpip_dE_woK0_woSidn_Sigmam_px");
  IMnpip_woSidn_Sigmam->SetLineColor(3); 
  IMnpip_woSidn_Sigmam->Draw("HEsame"); 

   
  TCanvas *cIMnpip_woSidn_ratio = new TCanvas("cIMnpip_woSidn_ratio","cIMnpip_woSidn_ratio");
  cIMnpip_woSidn_ratio->cd();
  TH1D* IMnpip_woSidn_ratio = IMnpip_woSidn_mc->Clone("IMnpip_woSidn_ratio");
  IMnpip_woSidn_ratio->Divide(IMnpip_woSidn_rdata);
  IMnpip_woSidn_ratio->Draw("HE");
  TF1* IMnpip_mod_gaus = new TF1("IMnpip_mod_gaus","gaus",1.1,1.7);
  IMnpip_mod_gaus->SetParameters(0,0.8);
  IMnpip_mod_gaus->SetParameters(1,1.2);
  IMnpip_mod_gaus->SetParameters(2,0.2);
  IMnpip_mod_gaus->SetLineColor(3);
  IMnpip_woSidn_ratio->Fit("IMnpip_mod_gaus","","",1.1,1.7);


  //Sigma-
  MMnmiss_IMnpim_dE_woK0_Lambda->Scale(0.40*nrdata_Sm/nLambda_Sm);
  MMnmiss_IMnpim_dE_woK0_Sigmap->Scale(0.01*nrdata_Sm/nSigmap_Sm);
  MMnmiss_IMnpim_dE_woK0_Sigmam->Scale(0.21*nrdata_Sm/nSigmam_Sm);
  MMnmiss_IMnpim_dE_woK0_S0->Scale(0.2*nrdata_Sm/nS0miss_Sm);
  MMnmiss_IMnpim_dE_woK0_Sppi0->Scale(0.03*nrdata_Sm/nSppi0_Sm);
  MMnmiss_IMnpim_dE_woK0_Smpi0->Scale(0.03*nrdata_Sm/nSmpi0_Sm);
  //
  MMnmiss_IMnpim_dE_woK0_woSp_Lambda->Scale(0.40*nrdata_Sm/nLambda_Sm);
  MMnmiss_IMnpim_dE_woK0_woSp_Sigmap->Scale(0.01*nrdata_Sm/nSigmap_Sm);
  MMnmiss_IMnpim_dE_woK0_woSp_Sigmam->Scale(0.21*nrdata_Sm/nSigmam_Sm);
  MMnmiss_IMnpim_dE_woK0_woSp_S0->Scale(0.2*nrdata_Sm/nS0miss_Sm);
  MMnmiss_IMnpim_dE_woK0_woSp_Sppi0->Scale(0.03*nrdata_Sm/nSppi0_Sm);
  MMnmiss_IMnpim_dE_woK0_woSp_Smpi0->Scale(0.03*nrdata_Sm/nSmpi0_Sm);
  //
  MMnmiss_IMnpim_dE_woK0_woSp_n_Lambda->Scale(0.40*nrdata_Sm/nLambda_Sm);
  MMnmiss_IMnpim_dE_woK0_woSp_n_Sigmap->Scale(0.01*nrdata_Sm/nSigmap_Sm);
  MMnmiss_IMnpim_dE_woK0_woSp_n_Sigmam->Scale(0.21*nrdata_Sm/nSigmam_Sm);
  MMnmiss_IMnpim_dE_woK0_woSp_n_S0->Scale(0.2*nrdata_Sm/nS0miss_Sm);
  MMnmiss_IMnpim_dE_woK0_woSp_n_Sppi0->Scale(0.03*nrdata_Sm/nSppi0_Sm);
  MMnmiss_IMnpim_dE_woK0_woSp_n_Smpi0->Scale(0.03*nrdata_Sm/nSmpi0_Sm);
  //
  MMnmiss_IMnpim_dE_woK0_woSidn_Lambda->Scale(0.40*nrdata_Sm/nLambda_Sm);
  MMnmiss_IMnpim_dE_woK0_woSidn_Sigmap->Scale(0.01*nrdata_Sm/nSigmap_Sm);
  MMnmiss_IMnpim_dE_woK0_woSidn_Sigmam->Scale(0.21*nrdata_Sm/nSigmam_Sm);
  MMnmiss_IMnpim_dE_woK0_woSidn_S0->Scale(0.2*nrdata_Sm/nS0miss_Sm);
  MMnmiss_IMnpim_dE_woK0_woSidn_Sppi0->Scale(0.03*nrdata_Sm/nSppi0_Sm);
  MMnmiss_IMnpim_dE_woK0_woSidn_Smpi0->Scale(0.03*nrdata_Sm/nSmpi0_Sm);

  TH2D *MMnmiss_IMnpim_mc = (TH2D*)MMnmiss_IMnpim_dE_woK0_Lambda->Clone("MMnmiss_IMnpim_mc");
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

  TH2D *MMnmiss_IMnpim_woSidn_mc = (TH2D*)MMnmiss_IMnpim_dE_woK0_woSidn_Lambda->Clone("MMnmiss_IMnpim_woSidn_mc");
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
  TH1D* IMnpim_woSidn_rdata = (TH1D*) MMnmiss_IMnpim_dE_woK0_woSidn_rdata->ProjectionX("IMnpim_dE_woK0_woSidn_rdata");
  IMnpim_woSidn_rdata->Draw("HE");
  TH1D* IMnpim_woSidn_mc = (TH1D*)MMnmiss_IMnpim_woSidn_mc->ProjectionX("MMnmiss_IMnpim_woSidn_mc_px");
  IMnpim_woSidn_mc->SetLineColor(6);
  IMnpim_woSidn_mc->Draw("HEsame");
  TH1D* IMnpim_woSidn_Lambda = (TH1D*)MMnmiss_IMnpim_dE_woK0_woSidn_Lambda->ProjectionX("MMnmiss_IMnpim_dE_woK0_woSidn_Lambda_px");
  IMnpim_woSidn_Lambda->SetLineColor(4); 
  IMnpim_woSidn_Lambda->Draw("HEsame"); 
  TH1D* IMnpim_woSidn_S0 = (TH1D*)MMnmiss_IMnpim_dE_woK0_woSidn_S0->ProjectionX("MMnmiss_IMnpim_dE_woK0_woSidn_S0_px");
  IMnpim_woSidn_S0->SetLineColor(5); 
  IMnpim_woSidn_S0->Draw("HEsame");; 
  TH1D* IMnpim_woSidn_Sigmap = (TH1D*)MMnmiss_IMnpim_dE_woK0_woSidn_Sigmap->ProjectionX("MMnmiss_IMnpim_dE_woK0_woSidn_Sigmap_px");
  IMnpim_woSidn_Sigmap->SetLineColor(2); 
  IMnpim_woSidn_Sigmap->Draw("HEsame"); 
  TH1D* IMnpim_woSidn_Sigmam = (TH1D*)MMnmiss_IMnpip_dE_woK0_woSidn_Sigmam->ProjectionX("MMnmiss_IMnpim_dE_woK0_woSidn_Sigmam_px");
  IMnpim_woSidn_Sigmam->SetLineColor(3); 
  IMnpim_woSidn_Sigmam->Draw("HEsame"); 

  TCanvas *cIMnpim_woSidn_ratio = new TCanvas("cIMnpim_woSidn_ratio","cIMnpim_woSidn_ratio");
  cIMnpim_woSidn_ratio->cd();
  TH1D* IMnpim_woSidn_ratio = (TH1D*)IMnpim_woSidn_mc->Clone("IMnpim_woSidn_ratio");
  IMnpim_woSidn_ratio->Divide(IMnpim_woSidn_rdata);
  IMnpim_woSidn_ratio->Draw("HE");
  



  //weighting and decomposition of BG
  int nbinx_Lambda = MMnmiss_IMnpip_dE_woK0_woSidn_Lambda->GetNbinsX();
  int nbiny_Lambda = MMnmiss_IMnpip_dE_woK0_woSidn_Lambda->GetNbinsY();
  int nbinx_S0 = MMnmiss_IMnpip_dE_woK0_woSidn_S0->GetNbinsX();
  int nbiny_S0 = MMnmiss_IMnpip_dE_woK0_woSidn_S0->GetNbinsY();
  TH2D* MMnmiss_IMnpip_dE_woK0_woSidn_Lambda_mod = (TH2D*)MMnmiss_IMnpip_dE_woK0_woSidn_Lambda->Clone("MMnmiss_IMnpip_dE_woK0_woSidn_Lambda_mod");
  TH2D* MMnmiss_IMnpim_dE_woK0_woSidn_Lambda_mod = (TH2D*)MMnmiss_IMnpim_dE_woK0_woSidn_Lambda->Clone("MMnmiss_IMnpim_dE_woK0_woSidn_Lambda_mod");
  TH2D* MMnmiss_IMnpip_dE_woK0_woSidn_S0_mod = (TH2D*)MMnmiss_IMnpip_dE_woK0_woSidn_S0->Clone("MMnmiss_IMnpip_dE_woK0_woSidn_S0_mod");
  TH2D* MMnmiss_IMnpim_dE_woK0_woSidn_S0_mod = (TH2D*)MMnmiss_IMnpim_dE_woK0_woSidn_S0->Clone("MMnmiss_IMnpim_dE_woK0_woSidn_S0_mod");
  MMnmiss_IMnpip_dE_woK0_woSidn_Lambda_mod->Clear();
  MMnmiss_IMnpim_dE_woK0_woSidn_Lambda_mod->Clear();
  for(int ix=0;ix<nbinx_Lambda;ix++){
    for(int iy=0;iy<nbiny_Lambda;iy++){
      double cont1 = MMnmiss_IMnpip_dE_woK0_woSidn_Lambda->GetBinContent(ix,iy);
      double cont2 = MMnmiss_IMnpim_dE_woK0_woSidn_Lambda->GetBinContent(ix,iy);
      double xval = MMnmiss_IMnpip_dE_woK0_woSidn_Lambda->GetXaxis()->GetBinCenter(ix);
      double yval = MMnmiss_IMnpip_dE_woK0_woSidn_Lambda->GetYaxis()->GetBinCenter(iy);
      double weight = 1.0;
      //if(yval<1.16) weight = TMath::Exp(yval);
      //else weight = 1;
      MMnmiss_IMnpip_dE_woK0_woSidn_Lambda_mod->Fill(xval,yval,cont1*weight);
      MMnmiss_IMnpim_dE_woK0_woSidn_Lambda_mod->Fill(xval,yval,cont2*weight);
    }
  }
  
  TCanvas *cMMnmiss_dE_woK0_woSidn_Lambda_mod = new TCanvas("MMnmiss_dE_woK0_woSidn_Lambda_mod","MMnmiss_dE_woK0_woSidn_Lambda_mod");
  cMMnmiss_dE_woK0_woSidn_Lambda_mod->cd(); 
  MMnmiss_dE_woK0_woSidn_rdata->Draw("HE");
  double rdatamax = MMnmiss_dE_woK0_woSidn_rdata->GetMaximum();
  TH1D *MMnmiss_dE_woK0_woSidn_Lambda_mod = (TH1D*) MMnmiss_IMnpip_dE_woK0_woSidn_Lambda_mod->ProjectionY("MMnmiss_dE_woK0_woSidn_Lambda_mod");
  MMnmiss_dE_woK0_woSidn_Lambda_mod->SetLineColor(4);
  double Lambdamax = MMnmiss_dE_woK0_woSidn_Lambda_mod->GetMaximum();
  MMnmiss_dE_woK0_woSidn_Lambda_mod->Scale(rdatamax/Lambdamax);
  MMnmiss_dE_woK0_woSidn_Lambda_mod->Draw("HEsame");
  TH1D *MMnmiss_dE_woK0_woSidn_S0_mod = (TH1D*) MMnmiss_IMnpip_dE_woK0_woSidn_S0_mod->ProjectionY("MMnmiss_dE_woK0_woSidn_S0_mod");
  MMnmiss_dE_woK0_woSidn_S0_mod->SetLineColor(5);
  MMnmiss_dE_woK0_woSidn_S0_mod->Draw("HEsame");
  MMnmiss_woSidn_Sigmap->Draw("HEsame");
  MMnmiss_woSidn_Sigmam->Draw("HEsame");
  TH1D* MMnmiss_dE_woK0_woSidn_mc_mod = (TH1D*)MMnmiss_dE_woK0_woSidn_Lambda_mod->Clone();
  MMnmiss_dE_woK0_woSidn_mc_mod->Add(MMnmiss_dE_woK0_woSidn_S0_mod);
  MMnmiss_dE_woK0_woSidn_mc_mod->Add(MMnmiss_woSidn_Sigmap);
  MMnmiss_dE_woK0_woSidn_mc_mod->Add(MMnmiss_woSidn_Sigmam);
  MMnmiss_dE_woK0_woSidn_mc_mod->SetLineColor(6);
  MMnmiss_dE_woK0_woSidn_mc_mod->Draw("HEsame");
  
  TCanvas *cMMnmiss_woK0_woSidn_mod_ratio = new TCanvas("MMnmiss_woK0_woSidn_mod_ratio","MMnmiss_woK0_woSidn_mod_ratio");
  cMMnmiss_woK0_woSidn_mod_ratio->cd(); 
  TH1D* MMnmiss_dE_woK0_woSidn_mc_mod_ratio = (TH1D*)MMnmiss_dE_woK0_woSidn_mc_mod->Clone("MMnmiss_woK0_woSidn_mod_ratio");
  MMnmiss_dE_woK0_woSidn_mc_mod_ratio->Divide(MMnmiss_dE_woK0_woSidn_rdata);
  MMnmiss_dE_woK0_woSidn_mc_mod_ratio->Draw("HE");

}
