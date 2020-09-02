void plot_ndEmom_mc(){


  //TFile *_file0 = TFile::Open("simIMpisigma_nSppim_v101.root");
  TFile *_file0 = TFile::Open("simIMpisigma_K0n_ns_v22.root");
  std::cout << _file0->Print() << std::endl;

  TCanvas *c1 = new TCanvas("c1","c1");
  CDHdE_mom_ncan_select->Draw("colz");
 
  TCanvas *c11 = new TCanvas("c11","c11");
  CDHdE_mom_ncan_select->ProjectionY("py",6,200)->Draw();
  

  double dEmom100MeVcbin = CDHdE_mom_ncan_select->GetXaxis()->FindBin(0.1);
  double dEmom200MeVcbin = CDHdE_mom_ncan_select->GetXaxis()->FindBin(0.2);
  double dEmom300MeVcbin = CDHdE_mom_ncan_select->GetXaxis()->FindBin(0.3);
  double dEmom400MeVcbin = CDHdE_mom_ncan_select->GetXaxis()->FindBin(0.4);
  double dEmom500MeVcbin = CDHdE_mom_ncan_select->GetXaxis()->FindBin(0.5);
  std::cout << "100 MeV bin " << dEmom100MeVcbin << std::endl;
  std::cout << "200 MeV bin " << dEmom200MeVcbin << std::endl;
  std::cout << "300 MeV bin " << dEmom300MeVcbin << std::endl;
  std::cout << "400 MeV bin " << dEmom400MeVcbin << std::endl;
  std::cout << "500 MeV bin " << dEmom500MeVcbin << std::endl;
  
  TH1D* dEnocut  = (TH1D*)CDHdE_mom_ncan_select->ProjectionY("dEnocut",2,100);
  TH1D* dE0MeV = (TH1D*)CDHdE_mom_ncan_select->ProjectionY("dE0MeV",2,dEmom100MeVcbin);
  TH1D* dE100MeV = (TH1D*)CDHdE_mom_ncan_select->ProjectionY("dE100MeV",dEmom100MeVcbin,dEmom200MeVcbin);
  TH1D* dE200MeV = (TH1D*)CDHdE_mom_ncan_select->ProjectionY("dE200MeV",dEmom200MeVcbin,dEmom300MeVcbin);
  TH1D* dE300MeV = (TH1D*)CDHdE_mom_ncan_select->ProjectionY("dE300MeV",dEmom300MeVcbin,dEmom400MeVcbin);
  TH1D* dE400MeV = (TH1D*)CDHdE_mom_ncan_select->ProjectionY("dE400MeV",dEmom400MeVcbin,dEmom500MeVcbin);
  
  TCanvas *c2 = new TCanvas("c2","c2");
  //dEnocut->SetLineColor(6);
  //dEnocut->SetMinimum(10);
  //dEnocut->GetYaxis()->SetRangeUser(10,12000);
  dEnocut->Draw("EH");
  dE0MeV->SetLineColor(2);
  dE0MeV->Draw("HEsame");
  dE100MeV->SetLineColor(3);
  dE200MeV->SetLineColor(4);
  dE300MeV->SetLineColor(5);
  dE400MeV->SetLineColor(6);
  dE100MeV->Draw("HEsame");
  dE200MeV->Draw("HEsame");
  dE300MeV->Draw("HEsame");
  dE400MeV->Draw("HEsame");
  
  TCanvas *cmom = new TCanvas("cmom","cmom");
  TH2D* ncan_mom_parentvtxr_select = (TH2D*)_file0->Get("ncan_mom_parentvtxr_select");
  int nbinfid = ncan_mom_parentvtxr_select->GetYaxis()->FindBin(10);
  TH1D* hmom_fid = (TH1D*)ncan_mom_parentvtxr_select->ProjectionX("hmom_fid",0,nbinfid);
  hmom_fid->Draw();
  TH1D* hmom_out = (TH1D*)ncan_mom_parentvtxr_select->ProjectionX("hmom_out",nbinfid+1,480);
  hmom_out->SetLineColor(2);
  hmom_out->Draw("same");

}
