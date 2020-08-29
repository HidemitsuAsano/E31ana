void disp_mcancestorana()
{

  //TFile *fGSp = TFile::Open("simIMpisigma_nSppim_v113.root","READ");
  TFile *fGSp = TFile::Open("simIMpisigma_nSppim_v114.root","READ");
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
  TH1D* diffmom_npip_ncan_except3g = (TH1D*)generation_diffmom_npip_ncan_select_sigma_Sp->ProjectionX("Spex",bin2nd+2,100);
  diffmom_npip_ncan_except3g->Add(diffmom_npip_ncan_2g);
  diffmom_npip_ncan_3g->Draw("HE");
  /*
  diffmom_npip_ncan_2g->SetLineColor(2);
  diffmom_npip_ncan_2g->Draw("HEsame");
  diffmom_npip_ncan_4g->SetLineColor(3);
  diffmom_npip_ncan_4g->Draw("HEsame");
  diffmom_npip_ncan_5g->SetLineColor(4);
  diffmom_npip_ncan_5g->Draw("HEsame");
  diffmom_npip_ncan_6g->SetLineColor(5);
  diffmom_npip_ncan_6g->Draw("HEsame");
  */
  diffmom_npip_ncan_except3g->SetLineColor(2);
  diffmom_npip_ncan_except3g->Draw("HEsame");
  c3->SetLogy();
  TLegend *leg = new TLegend(0.6,0.7,0.9,0.9);
  leg->AddEntry(diffmom_npip_ncan_3g,"3rd generation");
  leg->AddEntry(diffmom_npip_ncan_except3g,"others");
  leg->Draw();

  TH2F* vtxr_diffmom_npip_ncan_select_sigma_Sp
  = (TH2F*)fGSp->Get("vtxr_diffmom_npip_ncan_select_sigma");
  
  TCanvas *c5 = new TCanvas();
  vtxr_diffmom_npip_ncan_select_sigma_Sp->Draw("colz");
  c5->SetLogz();
  
  TCanvas *c7 = new TCanvas();
  int bin20cm = vtxr_diffmom_npip_ncan_select_sigma_Sp->GetYaxis()->FindBin(20);
  int bin58cm = vtxr_diffmom_npip_ncan_select_sigma_Sp->GetYaxis()->FindBin(58);
  int bin120cm = vtxr_diffmom_npip_ncan_select_sigma_Sp->GetYaxis()->FindBin(120);

  TH1D* pxfidinSp = (TH1D*)vtxr_diffmom_npip_ncan_select_sigma_Sp->ProjectionX("pxfidinSp",0,bin20cm);
  TH1D* pxfidoutSp = (TH1D*)vtxr_diffmom_npip_ncan_select_sigma_Sp->ProjectionX("pxfidoutSp",bin20cm,bin120cm);
  pxfidinSp->Draw("HE");
  pxfidoutSp->SetLineColor(2);
  pxfidoutSp->Draw("HEsame");
  c7->SetLogy();
  
  TLegend *leg7 = new TLegend(0.6,0.7,0.9,0.9);
  leg7->AddEntry(pxfidinSp,"n_{CDS} origin R<20 cm");
  leg7->AddEntry(pxfidoutSp,"n_{CDS} origin R>20 cm");
  leg7->Draw();
  

  //TFile *fGSm = TFile::Open("simIMpisigma_nSmpip_v113.root","READ");
  TFile *fGSm = TFile::Open("simIMpisigma_nSmpip_v114.root","READ");
  fGSm->cd();

  TH2F* generation_diffmom_npim_ncan_select_sigma_Sm 
  = (TH2F*)fGSm->Get("generation_diffmom_npim_ncan_select_sigma");
  
  TH2F* generation_diffMass_npim_ncan_select_sigma_Sm 
  = (TH2F*)fGSm->Get("generation_diffMass_npim_ncan_select_sigma");


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
  TH1D* diffmom_npim_ncan_except3g = (TH1D*)generation_diffmom_npim_ncan_select_sigma_Sm->ProjectionX("Smex",bin2nd+2,100);
  diffmom_npim_ncan_except3g->Add(diffmom_npim_ncan_2g);
  diffmom_npim_ncan_3g->Draw("HE");
  diffmom_npim_ncan_except3g->SetLineColor(2);
  diffmom_npim_ncan_except3g->Draw("HEsame");
  /*
  diffmom_npim_ncan_2g->SetLineColor(2);
  diffmom_npim_ncan_2g->Draw("HEsame");
  diffmom_npim_ncan_4g->SetLineColor(3);
  diffmom_npim_ncan_4g->Draw("HEsame");
  diffmom_npim_ncan_5g->SetLineColor(4);
  diffmom_npim_ncan_5g->Draw("HEsame");
  diffmom_npim_ncan_6g->SetLineColor(5);
  diffmom_npim_ncan_6g->Draw("HEsame");
  */


  c4->SetLogy();
  
  TLegend *leg2 = new TLegend(0.6,0.7,0.9,0.9);
  leg2->AddEntry(diffmom_npim_ncan_3g,"3rd generation");
  leg2->AddEntry(diffmom_npim_ncan_except3g,"others");
  leg2->Draw();

  TH2F* vtxr_diffmom_npim_ncan_select_sigma_Sm
  = (TH2F*)fGSm->Get("vtxr_diffmom_npim_ncan_select_sigma");
  TCanvas *c6 = new TCanvas();
  vtxr_diffmom_npim_ncan_select_sigma_Sm->Draw("colz");
  c6->SetLogz();


  TCanvas *c8 = new TCanvas();
  int bin20cm = vtxr_diffmom_npim_ncan_select_sigma_Sm->GetYaxis()->FindBin(20);
  int bin120cm = vtxr_diffmom_npim_ncan_select_sigma_Sm->GetYaxis()->FindBin(120);

  TH1D* pxfidinSm = (TH1D*)vtxr_diffmom_npim_ncan_select_sigma_Sm->ProjectionX("pxfidinSm",0,bin20cm);
  TH1D* pxfidoutSm = (TH1D*)vtxr_diffmom_npim_ncan_select_sigma_Sm->ProjectionX("pxfidoutSm",bin20cm,bin120cm);
  pxfidinSm->Draw("HE");
  pxfidoutSm->SetLineColor(2);
  pxfidoutSm->Draw("HEsame");
  c8->SetLogy();

  TLegend *leg8 = new TLegend(0.5,0.7,0.9,0.9);
  leg8->AddEntry(pxfidinSm,"n_{CDS} origin R<20 cm");
  leg8->AddEntry(pxfidoutSm,"n_{CDS} origin R>20 cm");
  leg8->Draw();
};
