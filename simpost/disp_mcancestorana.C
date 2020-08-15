void disp_mcancestorana()
{

  TFile *fGSp = TFile::Open("simIMpisigma_nSppim_v111.root","READ");
  fGSp->cd();

  TH2F* generation_diffmom_npip_ncan_select_sigma_Sp 
  = (TH2F*)fGSp->Get("generation_diffmom_npip_ncan_select_sigma");
  
  TCanvas *c1 = new TCanvas();
  generation_diffmom_npip_ncan_select_sigma_Sp->SetXTitle("diff. MC data - react. IM(n#pi^{+}) [GeV/c^{2}]");
  generation_diffmom_npip_ncan_select_sigma_Sp->SetYTitle("generation of neutron candidate");
  generation_diffmom_npip_ncan_select_sigma_Sp->Draw("colz");
  c1->SetLogz();

  TCanvas *c3 = new TCanvas();
  int bin2nd = generation_diffmom_npip_ncan_select_sigma_Sp->GetYaxis()->FindBin(2);
  TH1D* diffmom_npip_ncan_2g =  (TH1D*)generation_diffmom_npip_ncan_select_sigma_Sp->ProjectionX("Sp2g",bin2nd,bin2nd);
  TH1D* diffmom_npip_ncan_3g =  (TH1D*)generation_diffmom_npip_ncan_select_sigma_Sp->ProjectionX("Sp3g",bin2nd+1,bin2nd+1);
  TH1D* diffmom_npip_ncan_4g =  (TH1D*)generation_diffmom_npip_ncan_select_sigma_Sp->ProjectionX("Sp4g",bin2nd+2,bin2nd+2);
  TH1D* diffmom_npip_ncan_5g =  (TH1D*)generation_diffmom_npip_ncan_select_sigma_Sp->ProjectionX("Sp5g",bin2nd+3,bin2nd+3);
  TH1D* diffmom_npip_ncan_6g =  (TH1D*)generation_diffmom_npip_ncan_select_sigma_Sp->ProjectionX("Sp6g",bin2nd+4,bin2nd+4);
  diffmom_npip_ncan_3g->Draw("HE");
  diffmom_npip_ncan_2g->SetLineColor(2);
  diffmom_npip_ncan_2g->Draw("HEsame");
  diffmom_npip_ncan_4g->SetLineColor(3);
  diffmom_npip_ncan_4g->Draw("HEsame");
  diffmom_npip_ncan_5g->SetLineColor(4);
  diffmom_npip_ncan_5g->Draw("HEsame");
  diffmom_npip_ncan_6g->SetLineColor(5);
  diffmom_npip_ncan_6g->Draw("HEsame");
  c3->SetLogy();
  
  TH2F* vtxr_diffmom_npip_ncan_select_sigma_Sp
  = (TH2F*)fGSp->Get("vtxr_diffmom_npip_ncan_select_sigma");
  TCanvas *c5 = new TCanvas();
  vtxr_diffmom_npip_ncan_select_sigma_Sp->Draw("colz");
  c5->SetLogz();




  TFile *fGSm = TFile::Open("simIMpisigma_nSmpip_v111.root","READ");
  fGSm->cd();

  TH2F* generation_diffmom_npim_ncan_select_sigma_Sm 
  = (TH2F*)fGSm->Get("generation_diffmom_npim_ncan_select_sigma");
  
  TCanvas *c2 = new TCanvas();
  generation_diffmom_npim_ncan_select_sigma_Sm->SetXTitle("diff. MC data - react. IM(n#pi^{-}) [GeV/c^{2}]");
  generation_diffmom_npim_ncan_select_sigma_Sm->SetYTitle("generation of neutron candidate");
  generation_diffmom_npim_ncan_select_sigma_Sm->Draw("colz");
  c2->SetLogz();


  TCanvas *c4 = new TCanvas();
  TH1D* diffmom_npim_ncan_2g =  (TH1D*)generation_diffmom_npim_ncan_select_sigma_Sm->ProjectionX("Sm2g",bin2nd,bin2nd);
  TH1D* diffmom_npim_ncan_3g =  (TH1D*)generation_diffmom_npim_ncan_select_sigma_Sm->ProjectionX("Sm3g",bin2nd+1,bin2nd+1);
  TH1D* diffmom_npim_ncan_4g =  (TH1D*)generation_diffmom_npim_ncan_select_sigma_Sm->ProjectionX("Sm4g",bin2nd+2,bin2nd+2);
  TH1D* diffmom_npim_ncan_5g =  (TH1D*)generation_diffmom_npim_ncan_select_sigma_Sm->ProjectionX("Sm5g",bin2nd+3,bin2nd+3);
  TH1D* diffmom_npim_ncan_6g =  (TH1D*)generation_diffmom_npim_ncan_select_sigma_Sm->ProjectionX("Sm6g",bin2nd+4,bin2nd+4);
  diffmom_npim_ncan_3g->Draw("HE");
  diffmom_npim_ncan_2g->SetLineColor(2);
  diffmom_npim_ncan_2g->Draw("HEsame");
  diffmom_npim_ncan_4g->SetLineColor(3);
  diffmom_npim_ncan_4g->Draw("HEsame");
  diffmom_npim_ncan_5g->SetLineColor(4);
  diffmom_npim_ncan_5g->Draw("HEsame");
  diffmom_npim_ncan_6g->SetLineColor(5);
  diffmom_npim_ncan_6g->Draw("HEsame");
  c4->SetLogy();

  TH2F* vtxr_diffmom_npim_ncan_select_sigma_Sm
  = (TH2F*)fGSm->Get("vtxr_diffmom_npim_ncan_select_sigma");
  TCanvas *c6 = new TCanvas();
  vtxr_diffmom_npim_ncan_select_sigma_Sm->Draw("colz");
  c6->SetLogz();



};
