void comp_mat(){

  //bobbin + Iron solenoid
  //TFile *_file0 = TFile::Open("simIMpisigma_nSmpip_v53.root");
  //TH1D* mul_CDH0 = _file0->Get("mul_CDH");
  //TH2D* MMom_MMass_fid_beta_dE_woK0_wSid_0 = (TH2D*) _file0->Get("MMom_MMass_fid_beta_dE_woK0_wSid");
  //TH1D* px0 =  (TH1D*) MMom_MMass_fid_beta_dE_woK0_wSid_0->ProjectionX("px0");

  //bobbin + Coil + Iron
  TFile *_file1 = TFile::Open("simIMpisigma_nSmpip_v56.root");
  TH1D* mul_CDH1 = _file1->Get("mul_CDH");
  TH2D* MMom_MMass_fid_beta_dE_woK0_wSid_1 = (TH2D*) _file1->Get("MMom_MMass_fid_beta_dE_woK0_wSid");
  TH1D* px1 =  (TH1D*) MMom_MMass_fid_beta_dE_woK0_wSid_1->ProjectionX("px1");

  //All Iron
  TFile *_file2 = TFile::Open("simIMpisigma_nSmpip_v55.root");
  TH1D* mul_CDH2 = _file2->Get("mul_CDH");
  TH2D* MMom_MMass_fid_beta_dE_woK0_wSid_2 = (TH2D*) _file2->Get("MMom_MMass_fid_beta_dE_woK0_wSid");
  TH1D* px2 =  (TH1D*) MMom_MMass_fid_beta_dE_woK0_wSid_2->ProjectionX("px2");

  //DoraAir 
  TFile *_file3 = TFile::Open("simIMpisigma_nSmpip_DoraAir_v70.root");
  TH1D* mul_CDH3 = _file3->Get("mul_CDH");
  TH2D* MMom_MMass_fid_beta_dE_woK0_wSid_3 = (TH2D*) _file3->Get("MMom_MMass_fid_beta_dE_woK0_wSid");
  TH1D* px3 =  (TH1D*) MMom_MMass_fid_beta_dE_woK0_wSid_3->ProjectionX("px3");

  TCanvas *cmul = new TCanvas("cmul","cmul");
  cmul->cd();
  mul_CDH3->SetXTitle("CDH nhits");
  mul_CDH3->Draw("HE");
  mul_CDH2->SetLineColor(4);
  mul_CDH2->SetMarkerStyle(22);
  mul_CDH2->SetMarkerColor(4);
  mul_CDH2->Draw("HEsame");
  mul_CDH1->SetLineColor(2);
  mul_CDH1->SetMarkerStyle(22);
  mul_CDH1->SetMarkerColor(2);
  mul_CDH1->Draw("HEsame");
  //mul_CDH0->SetLineColor(3);
  //mul_CDH0->SetMarkerStyle(23);
  //mul_CDH0->SetMarkerColor(3);
  //mul_CDH0->Draw("HEsame");
  std::cout << "mean 3 " << mul_CDH3->GetMean() << std::endl;
  std::cout << "mean 2 " << mul_CDH2->GetMean() << std::endl;
  std::cout << "mean 1 " << mul_CDH1->GetMean() << std::endl;
  //std::cout << "mean 0 " << mul_CDH0->GetMean() << std::endl;

  
  TCanvas *cmiss = new TCanvas("cmiss","cmiss");
  cmiss->cd();
  px3->SetXTitle("Missing Mass [GeV/c^{2}]");
  px3->Draw("HE");
  double max3 = px3->GetMaximum();
  double max2 = px2->GetMaximum();
  double max1 = px1->GetMaximum();
  //double max0 = px0->GetMaximum();
  px2->Scale(max3/max2);
  px2->SetLineColor(2);
  px2->Draw("HEsame");
  px1->SetLineColor(4);
  px1->Scale(max3/max1);
  px1->Draw("HEsame");
  //px0->SetLineColor(3);
  //px0->Scale(max3/max0);
  //px0->Draw("HEsame");

  

}
