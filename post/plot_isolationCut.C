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
}
