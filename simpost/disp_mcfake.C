void disp_mcfake()
{

  TFile *fGSp = TFile::Open("simIMpisigma_nSppim_pippimn_v114_out.root","READ");

  TCanvas *c1 = new TCanvas();
  TH2F* IMnpim_IMnpip_dE_n_Sp = (TH2F*)fGSp->Get("IMnpim_IMnpip_dE_n");
  TH2F* IMnpim_IMnpip_dE_n_fake_Sp = (TH2F*)fGSp->Get("IMnpim_IMnpip_dE_n_fake");
  c1->cd();
  TH1D* IMnpip_true_Sp = (TH1D*)IMnpim_IMnpip_dE_n_Sp->ProjectionX("IMnpip_true_Sp");
  TH1D* IMnpip_fake_Sp = (TH1D*)IMnpim_IMnpip_dE_n_fake_Sp->ProjectionX("IMnpip_fake_Sp");
  IMnpip_fake_Sp->Print("base");
  IMnpip_true_Sp->Print("base");
  IMnpip_true_Sp->Add(IMnpip_fake_Sp,-1);
  IMnpip_true_Sp->GetXaxis()->SetRangeUser(0,1.7);
  IMnpip_true_Sp->SetMinimum(0);
  IMnpip_true_Sp->Draw("HE");
  IMnpip_fake_Sp->SetLineColor(2);
  IMnpip_fake_Sp->Draw("HEsame");

  
  TCanvas *c3 = new TCanvas();
  TH2F* MMnmiss_IMpippim_dE_wSid_Sp = (TH2F*)fGSp->Get("MMnmiss_IMpippim_dE_wSid");
  TH2F* MMnmiss_IMpippim_dE_wSid_fake_Sp = (TH2F*)fGSp->Get("MMnmiss_IMpippim_dE_wSid_fake");
  TH1D* MMnmiss_Sp_true = (TH1D*)MMnmiss_IMpippim_dE_wSid_Sp->ProjectionY("MMnmis_Sp_true");
  TH1D* MMnmiss_Sp_fake = (TH1D*)MMnmiss_IMpippim_dE_wSid_fake_Sp->ProjectionY("MMnmis_Sp_fake");
  MMnmiss_Sp_true->Add(MMnmiss_Sp_fake,-1);
  MMnmiss_Sp_true->Draw("HE");
  MMnmiss_Sp_fake->SetLineColor(2);
  MMnmiss_Sp_fake->Draw("HEsame");
  
  TCanvas *c5 = new TCanvas();
  TH2F* vtxr_diffMass_npip_ncan_wSid_n_mc_Sp = (TH2F*)fGSp->Get("vtxr_diffMass_npip_ncan_wSid_n_mc");
  const int bin20cm = vtxr_diffMass_npip_ncan_wSid_n_mc_Sp->GetYaxis()->FindBin(20);
  const int bin58cm = vtxr_diffMass_npip_ncan_wSid_n_mc_Sp->GetYaxis()->FindBin(58);
  const int bin120cm = vtxr_diffMass_npip_ncan_wSid_n_mc_Sp->GetYaxis()->FindBin(120);
  

  TH1D* diffMass_Sp_inCDH = (TH1D*)vtxr_diffMass_npip_ncan_wSid_n_mc_Sp->ProjectionX("diffMass_Sp_inCDH",0,bin58cm);
  TH1D* diffMass_Sp_outCDH = (TH1D*)vtxr_diffMass_npip_ncan_wSid_n_mc_Sp->ProjectionX("diffMass_Sp_outCDH",bin58cm+1,9999);
  diffMass_Sp_inCDH->Draw("HE");
  diffMass_Sp_outCDH->SetLineColor(2);
  diffMass_Sp_outCDH->Draw("HEsame");


  TFile *fGSm = TFile::Open("simIMpisigma_nSmpip_pippimn_v114_out.root","READ");

  TCanvas *c2 = new TCanvas();
  TH2F* IMnpim_IMnpip_dE_n_Sm = (TH2F*)fGSm->Get("IMnpim_IMnpip_dE_n");
  TH2F* IMnpim_IMnpip_dE_n_fake_Sm = (TH2F*)fGSm->Get("IMnpim_IMnpip_dE_n_fake");
  c2->cd();
  TH1D* IMnpip_true_Sm = (TH1D*)IMnpim_IMnpip_dE_n_Sm->ProjectionY("IMnpip_true_Sm");
  TH1D* IMnpip_fake_Sm = (TH1D*)IMnpim_IMnpip_dE_n_fake_Sm->ProjectionY("IMnpip_fake_Sm");
  IMnpip_fake_Sm->Print("base");
  IMnpip_true_Sm->Print("base");
  IMnpip_true_Sm->Add(IMnpip_fake_Sm,-1);
  IMnpip_true_Sm->GetXaxis()->SetRangeUser(0,1.7);
  IMnpip_true_Sm->SetMinimum(0);
  IMnpip_true_Sm->Draw("HE");
  IMnpip_fake_Sm->SetLineColor(2);
  IMnpip_fake_Sm->Draw("HEsame");


  TCanvas *c4 = new TCanvas();
  TH2F* MMnmiss_IMpippim_dE_wSid_Sm = (TH2F*)fGSm->Get("MMnmiss_IMpippim_dE_wSid");
  TH2F* MMnmiss_IMpippim_dE_wSid_fake_Sm = (TH2F*)fGSm->Get("MMnmiss_IMpippim_dE_wSid_fake");
  TH1D* MMnmiss_Sm_true = (TH1D*)MMnmiss_IMpippim_dE_wSid_Sm->ProjectionY("MMnmis_Sm_true");
  TH1D* MMnmiss_Sm_fake = (TH1D*)MMnmiss_IMpippim_dE_wSid_fake_Sm->ProjectionY("MMnmis_Sm_fake");
  MMnmiss_Sm_true->Add(MMnmiss_Sm_fake,-1);
  MMnmiss_Sm_true->Draw("HE");
  MMnmiss_Sm_fake->SetLineColor(2);
  MMnmiss_Sm_fake->Draw("HEsame");

  TCanvas *c6 = new TCanvas();
  TH2F* vtxr_diffMass_npim_ncan_wSid_n_mc_Sm = (TH2F*)fGSm->Get("vtxr_diffMass_npim_ncan_wSid_n_mc");
  const int bin20cm = vtxr_diffMass_npim_ncan_wSid_n_mc_Sm->GetYaxis()->FindBin(20);
  const int bin58cm = vtxr_diffMass_npim_ncan_wSid_n_mc_Sm->GetYaxis()->FindBin(58);
  const int bin120cm = vtxr_diffMass_npim_ncan_wSid_n_mc_Sm->GetYaxis()->FindBin(120);
  

  TH1D* diffMass_Sm_inCDH = (TH1D*)vtxr_diffMass_npim_ncan_wSid_n_mc_Sm->ProjectionX("diffMass_Sm_inCDH",0,bin58cm);
  TH1D* diffMass_Sm_outCDH = (TH1D*)vtxr_diffMass_npim_ncan_wSid_n_mc_Sm->ProjectionX("diffMass_Sm_outCDH",bin58cm+1,9999);
  diffMass_Sm_inCDH->Draw("HE");
  diffMass_Sm_outCDH->SetLineColor(2);
  diffMass_Sm_outCDH->Draw("HEsame");


};
