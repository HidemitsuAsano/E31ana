void comp_isolationcut(){

  TFile *file0 = TFile::Open("../simpost/simIMpisigma_npipiL_pippimn_v22_out.root");
  file0->cd();
  TH2D* MMnmiss_IMnpipi_wSid_0 = file0->Get("MMnmiss_IMnpipi_wSid");
  TH1D* MMnmiss_0 = (TH1D*)MMnmiss_IMnpipi_wSid_0->ProjectionY();

  TFile *file1 = TFile::Open("../simpost/simIMpisigma_npipiL_pippimn_v22_out_iso.root");
  TH2D* MMnmiss_IMnpipi_wSid_1 = file1->Get("MMnmiss_IMnpipi_wSid");
  TH1D* MMnmiss_1 = (TH1D*)MMnmiss_IMnpipi_wSid_1->ProjectionY();
  
  TFile *file2 = TFile::Open("../simpost/simIMpisigma_npipiL_pippimn_v22_out_iso2.root");
  TH2D* MMnmiss_IMnpipi_wSid_2 = file2->Get("MMnmiss_IMnpipi_wSid");
  TH1D* MMnmiss_2 = (TH1D*)MMnmiss_IMnpipi_wSid_2->ProjectionY();

  TFile *file3 = TFile::Open("../simpost/simIMpisigma_npipiL_pippimn_noconv_v9_out.root");
  TH2D* MMnmiss_IMnpipi_wSid_3 = file3->Get("MMnmiss_IMnpipi_wSid");
  TH1D* MMnmiss_3 = (TH1D*)MMnmiss_IMnpipi_wSid_3->ProjectionY();
  
  TFile *file4 = TFile::Open("../simpost/simIMpisigma_npipiL_pippimn_noconv_v9_out_iso.root");
  TH2D* MMnmiss_IMnpipi_wSid_4 = file4->Get("MMnmiss_IMnpipi_wSid");
  TH1D* MMnmiss_4 = (TH1D*)MMnmiss_IMnpipi_wSid_4->ProjectionY();

  TFile *file5 = TFile::Open("../simpost/simIMpisigma_npipiL_pippimn_noconv_v9_out_iso2.root");
  TH2D* MMnmiss_IMnpipi_wSid_5 = file5->Get("MMnmiss_IMnpipi_wSid");
  TH1D* MMnmiss_5 = (TH1D*)MMnmiss_IMnpipi_wSid_5->ProjectionY();

  const double max0 = MMnmiss_0->GetMaximum();
  const double max1 = MMnmiss_1->GetMaximum();
  const double max2 = MMnmiss_2->GetMaximum();
  const double max3 = MMnmiss_3->GetMaximum();
  const double max4 = MMnmiss_4->GetMaximum();
  const double max5 = MMnmiss_5->GetMaximum();
  
  const double scale1 = max0/max1;
  std::cout << "scale1 " << scale1 << std::endl;
  MMnmiss_1->Scale(scale1);
  const double scale2 = max0/max2;
  std::cout << "scale2 " << scale2 << std::endl;
  MMnmiss_2->Scale(scale2);
  const double scale3 = max3/max3;
  std::cout << "scale3 " << scale3 << std::endl;
  MMnmiss_3->Scale(scale3);
  const double scale4 = max3/max4;
  std::cout << "scale4 " << scale4 << std::endl;
  MMnmiss_4->Scale(scale4);
  const double scale5 = max3/max5;
  std::cout << "scale5 " << scale5 << std::endl;
  MMnmiss_5->Scale(scale5);

  TCanvas *c1 = new TCanvas("c1","c1");
  c1->cd();
  MMnmiss_0->Draw("HE");
  MMnmiss_1->SetLineColor(2);
  MMnmiss_1->Draw("HEsame");
  MMnmiss_2->SetLineColor(3);
  MMnmiss_2->Draw("HEsame");
  
  TCanvas *c2 = new TCanvas("c2","c2");
  MMnmiss_3->SetLineColor(4);
  MMnmiss_3->Draw("HE");
  MMnmiss_4->SetLineColor(5);
  //MMnmiss_4->Scale(1.18);
  MMnmiss_4->Draw("HEsame");
  MMnmiss_5->SetLineColor(6);
  //MMnmiss_5->Scale(1.5);
  MMnmiss_5->Draw("HEsame");



}
