void plot_isolationCut(){
  TCanvas *c1 = new  TCanvas("c1","c1");
  diff2d_CDC_CDH_pim->Draw("colz");

  TCanvas *c2 = new  TCanvas("c2","c2");
  diff2d_CDC_CDH_pip->Draw("colz");
  
  TCanvas *c3 = new TCanvas("c3","c3");
  diff2d_CDC_CDH_pim->ProjectionX()->Draw("HE");
  diff2d_CDC_CDH_pip->ProjectionX("pippx");
  pippx->SetLineColor(2);
  pippx->Draw("HEsame");
  
  TCanvas *c4 = new TCanvas("c4","c4");
  diff2d_CDC_CDH_pim->ProjectionY()->Draw("HE");
  diff2d_CDC_CDH_pip->ProjectionY("pippy");
  pippy->SetLineColor(2);
  pippy->Draw("HEsame");

  TCanvas *c5 = new TCanvas("c5","c5");
  MMnmiss_diffphi_CDC_CDH_pim->Draw("colz");

  TCanvas *c6 = new TCanvas("c6","c6");
  MMnmiss_diffphi_CDC_CDH_pip->Draw("colz");
  
  TCanvas *c7 = new TCanvas("c7","c7");
  MMnmiss_diffphi_CDC_CDH_pim_woK0_wSid->Draw("colz");

  TCanvas *c8 = new TCanvas("c8","c8");
  MMnmiss_diffphi_CDC_CDH_pip_woK0_wSid->Draw("colz");
  
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
}
