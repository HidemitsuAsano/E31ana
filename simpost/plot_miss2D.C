//plot_miss2D
//purpose : plot Missing mass, IM(npip) and so on. BG fit by modifying MC data.

void plot_miss2D()
{
  TH1::SetDefaultSumw2();
                                //rdata,nSp,nSm,Lambda,S0pi+pi-,nSppi0,nSmpi0,mcsum,K0nn
  const unsigned int colordef[9]={1,2,3,4,5,7,9,6,28};

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
  

  TFile *fileK0nn = TFile::Open("simIMpisigma_K0nn_pippimn_v8_out.root");

  
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
                                  0.38*nrdata_Sp/nLambda_Sp,//  0.40*nrdata_Sp/nLambda_Sp,
                                  0.18*nrdata_Sp/nS0miss_Sp,//  0.24*nrdata_Sp/nS0miss_Sp,
                                  0.02*nrdata_Sp/nSppi0_Sp,
                                  0.02*nrdata_Sp/nSmpi0_Sp};
  
  const double scaleFactor_Sm[7]={1.0,
                                  0.01*nrdata_Sm/nSigmap_Sm,
                                  0.21*nrdata_Sm/nSigmam_Sm,
                                  0.40*nrdata_Sm/nLambda_Sm,
                                  0.2*nrdata_Sm/nS0miss_Sm,
                                  0.03*nrdata_Sm/nSppi0_Sm,
                                  0.03*nrdata_Sm/nSmpi0_Sm};

  //array Sp mode
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


  TH2D* MMnmiss_IMnpim_dE_woK0[7]={
    MMnmiss_IMnpim_dE_woK0_rdata,
    MMnmiss_IMnpim_dE_woK0_Sigmap,
    MMnmiss_IMnpim_dE_woK0_Sigmam,
    MMnmiss_IMnpim_dE_woK0_Lambda,
    MMnmiss_IMnpim_dE_woK0_S0,
    MMnmiss_IMnpim_dE_woK0_Sppi0,
    MMnmiss_IMnpim_dE_woK0_Smpi0};

  TH2D* MMnmiss_IMnpim_dE_woK0_woSm[7]={
    MMnmiss_IMnpim_dE_woK0_woSp_rdata,
    MMnmiss_IMnpim_dE_woK0_woSp_Sigmap,
    MMnmiss_IMnpim_dE_woK0_woSp_Sigmam,
    MMnmiss_IMnpim_dE_woK0_woSp_Lambda,
    MMnmiss_IMnpim_dE_woK0_woSp_S0,
    MMnmiss_IMnpim_dE_woK0_woSp_Sppi0,
    MMnmiss_IMnpim_dE_woK0_woSp_Smpi0};
  
  TH2D* MMnmiss_IMnpim_dE_woK0_woSp[7]={
    MMnmiss_IMnpim_dE_woK0_woSp_rdata,
    MMnmiss_IMnpim_dE_woK0_woSp_Sigmap,
    MMnmiss_IMnpim_dE_woK0_woSp_Sigmam,
    MMnmiss_IMnpim_dE_woK0_woSp_Lambda,
    MMnmiss_IMnpim_dE_woK0_woSp_S0,
    MMnmiss_IMnpim_dE_woK0_woSp_Sppi0,
    MMnmiss_IMnpim_dE_woK0_woSp_Smpi0};
  
  TH2D* MMnmiss_IMnpim_dE_woK0_woSp_n[7]={
    MMnmiss_IMnpim_dE_woK0_woSp_n_rdata,
    MMnmiss_IMnpim_dE_woK0_woSp_n_Sigmap,
    MMnmiss_IMnpim_dE_woK0_woSp_n_Sigmam,
    MMnmiss_IMnpim_dE_woK0_woSp_n_Lambda,
    MMnmiss_IMnpim_dE_woK0_woSp_n_S0,
    MMnmiss_IMnpim_dE_woK0_woSp_n_Sppi0,
    MMnmiss_IMnpim_dE_woK0_woSp_n_Smpi0};

  TH2D* MMnmiss_IMnpim_dE_woK0_woSidn[7]={
    MMnmiss_IMnpim_dE_woK0_woSidn_rdata,
    MMnmiss_IMnpim_dE_woK0_woSidn_Sigmap,
    MMnmiss_IMnpim_dE_woK0_woSidn_Sigmam,
    MMnmiss_IMnpim_dE_woK0_woSidn_Lambda,
    MMnmiss_IMnpim_dE_woK0_woSidn_S0,
    MMnmiss_IMnpim_dE_woK0_woSidn_Sppi0,
    MMnmiss_IMnpim_dE_woK0_woSidn_Smpi0};

  //Sp
  for(int i=0;i<7;i++){
    MMnmiss_IMnpip_dE_woK0[i]->Scale(scaleFactor_Sp[i]);
    MMnmiss_IMnpip_dE_woK0_woSm[i]->Scale(scaleFactor_Sp[i]);
    MMnmiss_IMnpip_dE_woK0_woSm_n[i]->Scale(scaleFactor_Sp[i]);
    MMnmiss_IMnpip_dE_woK0_woSidn[i]->Scale(scaleFactor_Sp[i]);
  }
  
  //Sm
  for(int i=0;i<7;i++){
    MMnmiss_IMnpim_dE_woK0[i]->Scale(scaleFactor_Sm[i]);
    MMnmiss_IMnpim_dE_woK0_woSp[i]->Scale(scaleFactor_Sm[i]);
    MMnmiss_IMnpim_dE_woK0_woSp_n[i]->Scale(scaleFactor_Sm[i]);
    MMnmiss_IMnpim_dE_woK0_woSidn[i]->Scale(scaleFactor_Sm[i]);
  }
  
  //
  //ploting Sp mode
  //
  
  //w/o K0, including Sp/Sm mode
  TH2D *MMnmiss_IMnpip_woK0_mc = MMnmiss_IMnpip_dE_woK0[1]->Clone("MMnmiss_IMnpip_woK0_mc");
  //adding all MC data
  for(int i=2;i<7;i++)  MMnmiss_IMnpip_woK0_mc->Add(MMnmiss_IMnpip_dE_woK0[i]);
  
  TCanvas *cMMnmiss_IMnpip_woK0 = new TCanvas("cMMnmiss_IMnpip_woK0","cMMnmiss_IMnpip_woK0",1200,800);
  cMMnmiss_IMnpip_woK0->Divide(2,1);
  cMMnmiss_IMnpip_woK0->cd(1);
  MMnmiss_IMnpip_dE_woK0_rdata->SetTitle("real data");
  MMnmiss_IMnpip_dE_woK0_rdata->Draw("colz");
  cMMnmiss_IMnpip_woK0->cd(2);
  MMnmiss_IMnpip_woK0_mc->SetTitle("MC sum");
  MMnmiss_IMnpip_woK0_mc->Draw("colz");
  
  //w/o K0, w/o Sm mode
  TH2D *MMnmiss_IMnpip_woK0_woSm_mc = (TH2D*)MMnmiss_IMnpip_dE_woK0_woSm[1]->Clone("MMnmiss_IMnpip_woK0_woSm_mc");
  //adding all MC data
  for(int i=2;i<7;i++)MMnmiss_IMnpip_woK0_woSm_mc->Add(MMnmiss_IMnpip_dE_woK0_woSm[i]);
  
  //w/o K0, w/o Sm mode, selecting missing neutron
  TH2D *MMnmiss_IMnpip_woK0_woSm_n_mc = (TH2D*)MMnmiss_IMnpip_dE_woK0_woSm_n[1]->Clone("MMnmiss_IMnpip_woK0_woSm_n_mc");
  for(int i=2;i<7;i++)MMnmiss_IMnpip_woK0_woSm_n_mc->Add(MMnmiss_IMnpip_dE_woK0_woSm_n[i]);

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
  TH1D* IMnpip_woK0_woSm_n_mc = (TH1D*) MMnmiss_IMnpip_woK0_woSm_n_mc->ProjectionX("IMnpip_woK0_woSm_n_mc");
  IMnpip_woK0_woSm_n[0]->Draw("HE");
  for(int i=1;i<7;i++)IMnpip_woK0_woSm_n[i]->Draw("HEsame");
  IMnpip_woK0_woSm_n_mc->SetLineColor(6);
  IMnpip_woK0_woSm_n_mc->Draw("HEsame");

  //w/o K0, w/o ((Sp or Sm) & missing neutron)
  TH2D *MMnmiss_IMnpip_woK0_woSidn_mc = (TH2D*)MMnmiss_IMnpip_dE_woK0_woSidn[1]->Clone("MMnmiss_IMnpip_woK0_woSidn_mc");
  for(int i=2;i<7;i++) MMnmiss_IMnpip_woK0_woSidn_mc->Add(MMnmiss_IMnpip_dE_woK0_woSidn[i]);

  TCanvas *cMMnmiss_IMnpip_woSidn = new TCanvas("cMMnmiss_IMnpip_woSidn","cMMnmiss_IMnpip_woSidn",1200,800);
  cMMnmiss_IMnpip_woSidn->Divide(2,1);
  cMMnmiss_IMnpip_woSidn->cd(1);
  MMnmiss_IMnpip_dE_woK0_woSidn_rdata->SetTitle("real data");
  MMnmiss_IMnpip_dE_woK0_woSidn_rdata->Draw("colz");
  cMMnmiss_IMnpip_woSidn->cd(2);
  MMnmiss_IMnpip_woK0_woSidn_mc->SetTitle("MC sum");
  MMnmiss_IMnpip_woK0_woSidn_mc->Draw("colz");
  
  //projection to Missing mass (miss n & Sigma+/-)
  TCanvas *cMMnmiss_woK0_woSidn = new TCanvas("cMMnmiss_woK0_woSidn","cMMnmiss_woK0_woSidn");
  cMMnmiss_woK0_woSidn->cd();
  TH1D* MMnmiss_woK0_woSidn[7];
  for(int i=0;i<7;i++) MMnmiss_woK0_woSidn[i] = (TH1D*)MMnmiss_IMnpip_dE_woK0_woSidn[i]->ProjectionY(Form("MMnmiss_dE_woK0_woSidn_%s",name[i]));
  MMnmiss_woK0_woSidn[0]->Draw("HE");
  TH1D* MMnmiss_woK0_woSidn_mc = (TH1D*)MMnmiss_IMnpip_woK0_woSidn_mc->ProjectionY("MMnmiss_woK0_woSidn_mc");
  MMnmiss_woK0_woSidn_mc->SetLineColor(6);
  MMnmiss_woK0_woSidn_mc->Draw("HEsame");
  
  for(int i=1;i<7;i++){
    MMnmiss_woK0_woSidn[i]->SetLineColor(colordef[i]);
    MMnmiss_woK0_woSidn[i]->Draw("HEsame");
  }
  
  //Data/MC before modifying MC data
  TCanvas *cMMnmiss_woK0_woSidn_ratio = new TCanvas("cMMnmiss_woK0_woSidn_ratio","cMMnmiss_woK0_woSid_ratio");
  TH1D* MMnmiss_woK0_woSidn_ratio = MMnmiss_woK0_woSidn[0]->Clone("MMnmiss_woK0_woSidn_ratio");
  MMnmiss_woK0_woSidn_ratio->Divide(MMnmiss_woK0_woSidn_mc);
  MMnmiss_woK0_woSidn_ratio->SetTitle("Data/MC");
  MMnmiss_woK0_woSidn_ratio->Draw("HE");


  //projection to IMnpip (miss n & Sigma+/-)
  TCanvas *cIMnpip_woK0_woSidn = new TCanvas("cIMnpip_woK0_woSidn","cIMnpip_woK0_woSidn");
  cIMnpip_woK0_woSidn->cd();
  TH1D* IMnpip_woK0_woSidn[7];
  for(int i=0;i<7;i++)IMnpip_woK0_woSidn[i] = (TH1D*)MMnmiss_IMnpip_dE_woK0_woSidn[i]->ProjectionX(Form("IMnpip_dE_woK0_woSidn_%s",name[i]));
  IMnpip_woK0_woSidn[0]->Draw("HE");//rdata
  TH1D* IMnpip_woSidn_mc = (TH1D*)MMnmiss_IMnpip_woK0_woSidn_mc->ProjectionX("IMnpip_woK0_woSidn_mc");
  IMnpip_woSidn_mc->SetLineColor(6);
  IMnpip_woSidn_mc->Draw("HEsame");
  for(int i=1;i<7;i++){
    IMnpip_woK0_woSidn[i]->SetLineColor(colordef[i]);
    IMnpip_woK0_woSidn[i]->Draw("HEsame");
  }
  
  TCanvas *cIMnpip_woK0_woSidn_ratio = new TCanvas("cIMnpip_woK0_woSidn_ratio","cIMnpip_woK0_woSidn_ratio");
  cIMnpip_woK0_woSidn_ratio->cd();
  TH1D* IMnpip_woK0_woSidn_ratio = IMnpip_woK0_woSidn[0]->Clone("IMnpip_woK0_woSidn_ratio");
  IMnpip_woK0_woSidn_ratio->Divide(IMnpip_woK0_woSidn_mc);
  IMnpip_woK0_woSidn_ratio->SetTitle("IMnpip_woK0_woSidn Data/MC");
  IMnpip_woK0_woSidn_ratio->GetYaxis()->SetRangeUser(0,3);
  IMnpip_woK0_woSidn_ratio->Draw("HE");
  TF1* IMnpip_mod_gaus = new TF1("IMnpip_mod_gaus","gaus",1.1,1.7);
  IMnpip_mod_gaus->SetParameters(0,0.8);
  IMnpip_mod_gaus->SetParameters(1,1.2);
  IMnpip_mod_gaus->SetParameters(2,0.2);
  IMnpip_mod_gaus->SetLineColor(3);
  IMnpip_woK0_woSidn_ratio->Fit("IMnpip_mod_gaus","","",1.1,1.7);

  //
  //Sigma-
  TH2D *MMnmiss_IMnpim_woK0_mc = (TH2D*)MMnmiss_IMnpim_dE_woK0[1]->Clone("MMnmiss_IMnpim_woK0_mc");
  //adding all MC data
  for(int i=2;i<7;i++)MMnmiss_IMnpim_woK0_mc->Add(MMnmiss_IMnpim_dE_woK0[i]);

  TCanvas *cMMnmiss_IMnpim_woK0 = new TCanvas("cMMnmiss_IMnpim_woK0","cMMnmiss_IMnpim_woK0",1200,800);
  cMMnmiss_IMnpim_woK0->Divide(2,1);
  cMMnmiss_IMnpim_woK0->cd(1);
  MMnmiss_IMnpim_dE_woK0_rdata->SetTitle("real data");
  MMnmiss_IMnpim_dE_woK0_rdata->Draw("colz");
  cMMnmiss_IMnpim_woK0->cd(2);
  MMnmiss_IMnpim_woK0_mc->SetTitle("MC sum");
  MMnmiss_IMnpim_woK0_mc->Draw("colz");
  
  //w/o K0, w/o Sp mode
  TH2D *MMnmiss_IMnpim_woK0_woSp_mc = (TH2D*)MMnmiss_IMnpim_dE_woK0_woSp[1]->Clone("MMnmiss_IMnpim_woK0_woSp_mc");
  //adding all MC data
  for(int i=2;i<7;i++){
    MMnmiss_IMnpim_dE_woK0_woSp[i]->Print("base");
    MMnmiss_IMnpim_woK0_woSp_mc->Add(MMnmiss_IMnpim_dE_woK0_woSp[i]);
  }
  //w/o K0, w/o Sp mode, selecting missing neutron
  TH2D *MMnmiss_IMnpim_woK0_woSp_n_mc = (TH2D*)MMnmiss_IMnpim_dE_woK0_woSp_n[1]->Clone("MMnmiss_IMnpim_woK0_woSp_n_mc");
  for(int i=2;i<7;i++)MMnmiss_IMnpim_woK0_woSp_n_mc->Add(MMnmiss_IMnpim_dE_woK0_woSp_n[i]);

  TCanvas *cMMnmiss_IMnpim_woK0_woSp = new TCanvas("cMMnmiss_IMnpim_woK0_woSp","cMMnmiss_IMnpim_woK0_woSp",1200,800);
  cMMnmiss_IMnpim_woK0_woSp->Divide(2,1);
  cMMnmiss_IMnpim_woK0_woSp->cd(1);
  MMnmiss_IMnpim_dE_woK0_woSp_rdata->SetTitle("real data");
  MMnmiss_IMnpim_dE_woK0_woSp_rdata->Draw("colz");
  cMMnmiss_IMnpim_woK0_woSp->cd(2);
  MMnmiss_IMnpim_woK0_woSp_mc->SetTitle("MC sum");
  MMnmiss_IMnpim_woK0_woSp_mc->Draw("colz");

  TCanvas *cIMnpim_woK0_woSp_n = new TCanvas("cIMnpim_woK0_woSp_n","cIMnpim_woK0_woSp_n");
  cIMnpim_woK0_woSp_n->cd();
  TH1D *IMnpim_woK0_woSp_n[7];
  for(int i=0;i<7;i++){
    IMnpim_woK0_woSp_n[i] = (TH1D*)MMnmiss_IMnpim_dE_woK0_woSp_n[i]->ProjectionX(Form("IMnpim_woK0_woSp_n_%s",name[i]));
    IMnpim_woK0_woSp_n[i]->SetLineColor(colordef[i]);
  }
  TH1D* IMnpim_woK0_woSp_n_mc = (TH1D*) MMnmiss_IMnpim_woK0_woSp_n_mc->ProjectionX("IMnpim_woK0_woSp_n_mc");
  IMnpim_woK0_woSp_n[0]->Draw("HE");
  for(int i=1;i<7;i++)IMnpim_woK0_woSp_n[i]->Draw("HEsame");
  IMnpim_woK0_woSp_n_mc->SetLineColor(6);
  IMnpim_woK0_woSp_n_mc->Draw("HEsame");


  //w/o K0, w/o ((Sp or Sm) & missing neutron)
  TH2D *MMnmiss_IMnpim_woK0_woSidn_mc = (TH2D*)MMnmiss_IMnpim_dE_woK0_woSidn[1]->Clone("MMnmiss_IMnpim_woK0_woSidn_mc");
  for(int i=2;i<7;i++)MMnmiss_IMnpim_woK0_woSidn_mc->Add(MMnmiss_IMnpim_dE_woK0_woSidn[i]);

  TCanvas *cMMnmiss_IMnpim_woK0_woSidn = new TCanvas("cMMnmiss_IMnpim_woK0_woSidn","cMMnmiss_IMnpim_woK0_woSidn",1200,800);
  cMMnmiss_IMnpim_woK0_woSidn->Divide(2,1);
  cMMnmiss_IMnpim_woK0_woSidn->cd(1);
  MMnmiss_IMnpim_dE_woK0_woSidn_rdata->SetTitle("real data");
  MMnmiss_IMnpim_dE_woK0_woSidn_rdata->Draw("colz");
  cMMnmiss_IMnpim_woK0_woSidn->cd(2);
  MMnmiss_IMnpim_woK0_woSidn_mc->SetTitle("MC sum");
  MMnmiss_IMnpim_woK0_woSidn_mc->Draw("colz");
  
  TCanvas *cIMnpim_woK0_woSidn = new TCanvas("IMnpim_woK0_woSidn","IMnpim_woK0_woSidn");
  cIMnpim_woK0_woSidn->cd();
  TH1D* IMnpim_woK0_woSidn[7];
  for(int i=0;i<7;i++)IMnpim_woK0_woSidn[i] = (TH1D*)MMnmiss_IMnpim_dE_woK0_woSidn[i]->ProjectionX(Form("IMnpim_dE_woK0_woSidn_%s",name[i]));
  TH1D* IMnpim_woK0_woSidn_mc = MMnmiss_IMnpim_woK0_woSidn_mc->ProjectionX("IMnpim_woK0_woSidn_mc");
  IMnpim_woK0_woSidn[0]->Draw("HE");
  IMnpim_woK0_woSidn_mc->SetLineColor(6);
  IMnpim_woK0_woSidn_mc->Draw("HEsame");
  for(int i=1;i<7;i++){
    IMnpim_woK0_woSidn[i]->SetLineColor(colordef[i]);
    IMnpim_woK0_woSidn[i]->Draw("HEsame");
  }
  
  TCanvas *cIMnpim_woK0_woSidn_ratio = new TCanvas("cIMnpim_woK0_woSidn_ratio","cIMnpim_woK0_woSidn_ratio");
  cIMnpim_woK0_woSidn_ratio->cd();
  TH1D* IMnpim_woK0_woSidn_ratio = (TH1D*)IMnpim_woK0_woSidn[0]->Clone("IMnpim_woK0_woSidn_ratio");
  IMnpim_woK0_woSidn_ratio->Divide(IMnpim_woK0_woSidn_mc);
  IMnpim_woK0_woSidn_ratio->SetTitle("IMnpim_woK0_woSidn Data/MC");
  IMnpim_woK0_woSidn_ratio->GetYaxis()->SetRangeUser(0,3);
  IMnpim_woK0_woSidn_ratio->Draw("HE");
  



  //weighting and decomposition of BG
  int nbinx_MM_npip = MMnmiss_IMnpip_dE_woK0_woSidn[0]->GetNbinsX();
  int nbiny_MM_npip = MMnmiss_IMnpip_dE_woK0_woSidn[0]->GetNbinsY();
  TH2D* MMnmiss_IMnpip_dE_woK0_woSidn_mod[7];
  TH2D* MMnmiss_IMnpim_dE_woK0_woSidn_mod[7];
  for(int i=0;i<7;i++){
    MMnmiss_IMnpip_dE_woK0_woSidn_mod[i] = (TH2D*)MMnmiss_IMnpip_dE_woK0_woSidn[i]->Clone(Form("MMnmiss_IMnpip_dE_woK0_woSidn_%s_mod",name[i]));
    MMnmiss_IMnpim_dE_woK0_woSidn_mod[i] = (TH2D*)MMnmiss_IMnpim_dE_woK0_woSidn[i]->Clone(Form("MMnmiss_IMnpim_dE_woK0_woSidn_%s_mod",name[i]));
  }  
  //for(int i=3;i<7;i++){
  for(int i=1;i<4;i++){
    MMnmiss_IMnpip_dE_woK0_woSidn_mod[i]->Reset("ICESM");
    MMnmiss_IMnpim_dE_woK0_woSidn_mod[i]->Reset("ICESM");
  }

  for(int ix=0;ix<nbinx_MM_npip;ix++){
    for(int iy=0;iy<nbiny_MM_npip;iy++){
      double cont1 = MMnmiss_IMnpip_dE_woK0_woSidn[3]->GetBinContent(ix,iy);
      double cont2 = MMnmiss_IMnpim_dE_woK0_woSidn[3]->GetBinContent(ix,iy);
      double xval = MMnmiss_IMnpip_dE_woK0_woSidn[3]->GetXaxis()->GetBinCenter(ix);
      double yval = MMnmiss_IMnpip_dE_woK0_woSidn[3]->GetYaxis()->GetBinCenter(iy);
      double weight = 1.0;
      if(yval<1.16){
        weight = TMath::Exp(yval/2.0);
      }else{
        weight = 1.0;
      }
      MMnmiss_IMnpip_dE_woK0_woSidn_mod[3]->Fill(xval,yval,cont1*weight);
      MMnmiss_IMnpim_dE_woK0_woSidn_mod[3]->Fill(xval,yval,cont2*weight);
    }
  }
  
  TCanvas *cMMnmiss_dE_woK0_woSidn_mod = new TCanvas("MMnmiss_dE_woK0_woSidn_mod","MMnmiss_dE_woK0_woSidn_mod");
  cMMnmiss_dE_woK0_woSidn_mod->cd(); 
  MMnmiss_woK0_woSidn[0]->Draw("HE");
  double rdatamax = MMnmiss_woK0_woSidn[0]->GetMaximum();
  TH1D* MMnmiss_dE_woK0_woSidn_mod[7];
  for(int i=0;i<7;i++){
    MMnmiss_dE_woK0_woSidn_mod[i] = (TH1D*)MMnmiss_IMnpip_dE_woK0_woSidn_mod[i]->ProjectionY(Form("MMnmiss_dE_woK0_woSidn_mod_%s",name[i]));
    MMnmiss_dE_woK0_woSidn_mod[i]->SetLineColor(colordef[i]);
  }
  //double Lambdamax = MMnmiss_dE_woK0_woSidn_mod[3]->GetMaximum();
  //MMnmiss_dE_woK0_woSidn_mod[3]->Scale(rdatamax/Lambdamax);
  //MMnmiss_dE_woK0_woSidn_mod[3]->Draw("HEsame");
  for(int i=3;i<7;i++)MMnmiss_dE_woK0_woSidn_mod[i]->Draw("HEsame");
  TH1D* MMnmiss_dE_woK0_woSidn_mc_mod = (TH1D*)MMnmiss_dE_woK0_woSidn_mod[1]->Clone();
  for(int i=3;i<7;i++)MMnmiss_dE_woK0_woSidn_mc_mod->Add(MMnmiss_dE_woK0_woSidn_mod[i]);
  MMnmiss_dE_woK0_woSidn_mc_mod->SetLineColor(6);
  MMnmiss_dE_woK0_woSidn_mc_mod->Draw("HEsame");
  
  TCanvas *cMMnmiss_woK0_woSidn_mod_ratio = new TCanvas("MMnmiss_woK0_woSidn_mod_ratio","MMnmiss_woK0_woSidn_mod_ratio");
  cMMnmiss_woK0_woSidn_mod_ratio->cd(); 
  TH1D* MMnmiss_dE_woK0_woSidn_mc_mod_ratio = (TH1D*)MMnmiss_woK0_woSidn[0]->Clone("MMnmiss_woK0_woSidn_mod_ratio");
  MMnmiss_dE_woK0_woSidn_mc_mod_ratio->Divide(MMnmiss_dE_woK0_woSidn_mc_mod);
  MMnmiss_dE_woK0_woSidn_mc_mod_ratio->SetTitle("MMnmiss_dE_woK0_woSidn Data/MC");
  MMnmiss_dE_woK0_woSidn_mc_mod_ratio->GetYaxis()->SetRangeUser(-0.1,3);
  MMnmiss_dE_woK0_woSidn_mc_mod_ratio->Draw("HE");

  
  TCanvas *c = NULL;
  TSeqCollection *SCol = gROOT->GetListOfCanvases();
  int size = SCol->GetSize();
  TIter next(SCol);
  TString pdfname="plot_miss2D_out.pdf";
  for(int i=0;i<size;i++){
    c= (TCanvas*)next();
    c->Draw();
    c->cd();
    if(i==0) c->Print(pdfname+"(",Form("Title:%s",c->GetTitle()));
    else if(i==size-1)c->Print(pdfname+")",Form("Title:%s",c->GetTitle())); 
    else c->Print(pdfname,Form("Title:%s",c->GetTitle())); 
  }


}
