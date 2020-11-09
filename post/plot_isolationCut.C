void plot_isolationCut(const char *filename="evanaIMpisigma_npippim_v201_out_iso.root"){
  
  TFile *f = new TFile(filename,"READ");
  f->cd();

  TCanvas *c1 = new  TCanvas("c1","c1");
  TH2D* diff2d_CDC_CDH_pim = (TH2D*)f->Get("diff2d_CDC_CDH_pim");
  diff2d_CDC_CDH_pim->Draw("colz");

  TCanvas *c2 = new  TCanvas("c2","c2");
  TH2D* diff2d_CDC_CDH_pip = (TH2D*)f->Get("diff2d_CDC_CDH_pip");
  diff2d_CDC_CDH_pip->Draw("colz");
  
  TCanvas *c3 = new TCanvas("c3","c3");
  TH1D *pimpx = (TH1D*) diff2d_CDC_CDH_pim->ProjectionX();;
  TH1D *pippx = (TH1D*) diff2d_CDC_CDH_pip->ProjectionX();
  pimpx->Draw("HE");
  pippx->SetLineColor(2);
  pippx->Draw("HEsame");
  
  TCanvas *c4 = new TCanvas("c4","c4");
  TH1D* pimpy = (TH1D*) diff2d_CDC_CDH_pim->ProjectionY();
  TH1D* pippy = (TH1D*) diff2d_CDC_CDH_pip->ProjectionY();
  pimpy->Draw("HE");
  pippy->SetLineColor(2);
  pippy->Draw("HEsame");

  TCanvas *c5 = new TCanvas("c5","c5");
  TH2D* diff2d_CDC_CDH_pim_phi_tof = (TH2D*)f->Get("diff2d_CDC_CDH_pim_phi_tof");
  diff2d_CDC_CDH_pim_phi_tof->RebinY(2);
  diff2d_CDC_CDH_pim_phi_tof->Draw("colz");

  TCanvas *c6 = new TCanvas("c6","c6");
  TH2D* diff2d_CDC_CDH_pim_z_tof = (TH2D*)f->Get("diff2d_CDC_CDH_pim_z_tof");
  diff2d_CDC_CDH_pim_z_tof->RebinY(2);
  diff2d_CDC_CDH_pim_z_tof->Draw("colz");

  TCanvas *c7 = new TCanvas("c7","c7");
  TH2D* diff2d_CDC_CDH_pip_phi_tof = (TH2D*)f->Get("diff2d_CDC_CDH_pip_phi_tof");
  diff2d_CDC_CDH_pip_phi_tof->RebinY(2);
  diff2d_CDC_CDH_pip_phi_tof->Draw("colz");

  TCanvas *c8 = new TCanvas("c8","c8");
  TH2D* diff2d_CDC_CDH_pip_z_tof = (TH2D*)f->Get("diff2d_CDC_CDH_pip_z_tof");
  diff2d_CDC_CDH_pip_z_tof->RebinY(2);
  diff2d_CDC_CDH_pip_z_tof->Draw("colz");

  TCanvas *c1_1 = new  TCanvas("c1_1","c1_1");
  TH2D* diff2d_CDC_CDH_pim_woSid_won = (TH2D*)f->Get("diff2d_CDC_CDH_pim_woSid_won");
  diff2d_CDC_CDH_pim_woSid_won->Draw("colz");

  TCanvas *c2_2 = new  TCanvas("c2_2","c2_2");
  TH2D* diff2d_CDC_CDH_pip_woSid_won = (TH2D*)f->Get("diff2d_CDC_CDH_pip_woSid_won");
  diff2d_CDC_CDH_pip_woSid_won->Draw("colz");
  
  TCanvas *c3_3 = new TCanvas("c3_3","c3_3");
  TH1D *pimpx_woSid_won = (TH1D*) diff2d_CDC_CDH_pim_woSid_won->ProjectionX();;
  TH1D *pippx_woSid_won = (TH1D*) diff2d_CDC_CDH_pip_woSid_won->ProjectionX();
  pimpx_woSid_won->Draw("HE");
  pippx_woSid_won->SetLineColor(2);
  pippx_woSid_won->Draw("HEsame");
  
  TCanvas *c4_4 = new TCanvas("c4_4","c4_4");
  TH1D* pimpy_woSid_won = (TH1D*) diff2d_CDC_CDH_pim_woSid_won->ProjectionY();
  TH1D* pippy_woSid_won = (TH1D*) diff2d_CDC_CDH_pip_woSid_won->ProjectionY();
  pimpy_woSid_won->Draw("HE");
  pippy_woSid_won->SetLineColor(2);
  pippy_woSid_won->Draw("HEsame");
  

  TCanvas *c5_5 = new TCanvas("c5_5","c5_5");
  TH1D* npimtof_woSid_won = (TH1D*) diff2d_CDC_CDH_pim_phi_tof->ProjectionY();
  TH1D* npiptof_woSid_won = (TH1D*) diff2d_CDC_CDH_pip_phi_tof->ProjectionY();
  npimtof_woSid_won->Draw("HE");
  npiptof_woSid_won->SetLineColor(2);
  npiptof_woSid_won->Draw("HEsame");

  /*
  TCanvas *c5 = new TCanvas("c5","c5");
  MMnmiss_diffphi_CDC_CDH_pim->Draw("colz");

  TCanvas *c6 = new TCanvas("c6","c6");
  MMnmiss_diffphi_CDC_CDH_pip->Draw("colz");
  
  TCanvas *c7 = new TCanvas("c7","c7");
  MMnmiss_diffphi_CDC_CDH_pim_woK0_wSid->Draw("colz");

  TCanvas *c8 = new TCanvas("c8","c8");
  MMnmiss_diffphi_CDC_CDH_pip_woK0_wSid->Draw("colz");
  */

  
  TCanvas *c9 = new TCanvas("c9","c9");
  pimmom_diffphi_CDC_CDH_pim->Draw("colz");
  
  TCanvas *c10 = new TCanvas("c10","c10");
  pipmom_diffphi_CDC_CDH_pip->Draw("colz");

  TCanvas *c11 = new TCanvas("c11","c11");
  nmom_diffphi_CDC_CDH_pim->Draw("colz");

  TCanvas *c12 = new TCanvas("c12","c12");
  nmom_diffphi_CDC_CDH_pip->Draw("colz");
  
  TCanvas *c13 = new TCanvas("c13","c13");
  nmom_diffphi_CDC_CDH_pim_woK0_wSid->Draw("colz");

  TCanvas *c14 = new TCanvas("c14","c14");
  nmom_diffphi_CDC_CDH_pip_woK0_wSid->Draw("colz");

  TCanvas *c15 = new TCanvas("c15","c15");
  pimmom_diffz_CDC_CDH_pim->Draw("colz");

  TCanvas *c16 = new TCanvas("c16","c16");
  pipmom_diffz_CDC_CDH_pip->Draw("colz");

  TCanvas *c17 = new TCanvas("c17","c17");
  nmom_diffz_CDC_CDH_pim->Draw("colz");
  
  TCanvas *c18 = new TCanvas("c18","c18");
  nmom_diffz_CDC_CDH_pip->Draw("colz");
  

  /*
  TCanvas *c = NULL;
  TSeqCollection *SCol = gROOT->GetListOfCanvases();
  int size = SCol->GetSize();
  //TIter next(gROOT->GetListOfCanvases());
  TIter next(SCol);
  TString pdfname = std::string("isolation.pdf");
  for(int i=0; i<size; i++) {
    c= (TCanvas*)next();
    c->Draw();
    c->cd();
    c->Modified();
    c->Update();
    if(i==0) c->Print(pdfname+"(",Form("Title:%s",c->GetTitle()));
    else if(i==size-1)c->Print(pdfname+")",Form("Title:%s",c->GetTitle()));
    else c->Print(pdfname,Form("Title:%s",c->GetTitle()));
  }
 */

}

