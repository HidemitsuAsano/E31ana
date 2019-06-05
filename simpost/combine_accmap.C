void combine_accmap(){

  TFile *file1 = new TFile("simIMpisigma_nSppim_DoraAir_v45_v46_acc.root","READ");
  TFile *file2 = new TFile("simIMpisigma_nSppim_DoraAir_v47_v48_acc.root","READ");
  
  TH1::SetDefaultSumw2();
  
  TH2F* q_IMnpipi_wSid_n_acc_nocut = (TH2F*)file1->Get("q_IMnpipi_wSid_n_acc");
  TCanvas *cacc_nocut = new TCanvas("cacc_nocut","acc_nocut");
  q_IMnpipi_wSid_n_acc_nocut->Draw("colz");
  
  TH2F* q_IMnpipi_woK0_wSid_n_acc_nocut = (TH2F*)file1->Get("q_IMnpipi_woK0_wSid_n_acc");
  TCanvas *cacc_woK0_nocut = new TCanvas("cacc_woK0_nocut","acc_woK0_nocut");
  q_IMnpipi_woK0_wSid_n_acc_nocut->Draw("colz");

  TH2F* q_IMnpipi_woK0_wSid_n_Sp_acc_nocut = (TH2F*)file1->Get("q_IMnpipi_woK0_wSid_n_Sp_acc");
  TCanvas *cacc_woK0_Sp_nocut = new TCanvas("cacc_woK0_Sp_nocut","acc_woK0_Sp_nocut");
  q_IMnpipi_woK0_wSid_n_Sp_acc_nocut->Draw("colz");

  TH2F* q_IMnpipi_wSid_n_acc_cut1400 = (TH2F*)file1->Get("q_IMnpipi_wSid_n_acc");
  TCanvas *cacc_cut1400 = new TCanvas("cacc_cut1400","acc_cut1400");
  q_IMnpipi_wSid_n_acc_cut1400->Draw("colz");
  
  TH2F* q_IMnpipi_woK0_wSid_n_acc_cut1400 = (TH2F*)file1->Get("q_IMnpipi_woK0_wSid_n_acc");
  TCanvas *cacc_woK0_cut1400 = new TCanvas("cacc_woK0_cut1400","acc_woK0_cut1400");
  q_IMnpipi_woK0_wSid_n_acc_cut1400->Draw("colz");

  TH2F* q_IMnpipi_woK0_wSid_n_Sp_acc_cut1400 = (TH2F*)file2->Get("q_IMnpipi_woK0_wSid_n_Sp_acc");
  TCanvas *cacc_woK0_Sp_cut1400 = new TCanvas("cacc_woK0_Sp_cut1400","acc_woK0_Sp_cut1400");
  q_IMnpipi_woK0_wSid_n_Sp_acc_cut1400->Draw("colz");
   
  //ADD
  //including K0
  TH2F* q_IMnpipi_wSid_n_acc_add = (TH2F*)q_IMnpipi_wSid_n_acc_nocut->Clone();
  q_IMnpipi_wSid_n_acc_add->Add(q_IMnpipi_wSid_n_acc_cut1400);
  TCanvas *cacc_add = new TCanvas("cacc_add","acc_add");
  q_IMnpipi_wSid_n_acc_add->Draw("colz");
  
  //w/o K0
  TH2F* q_IMnpipi_woK0_wSid_n_acc_add = (TH2F*)q_IMnpipi_woK0_wSid_n_acc_nocut->Clone();
  q_IMnpipi_woK0_wSid_n_acc_add->Add(q_IMnpipi_woK0_wSid_n_acc_cut1400);
  TCanvas *cacc_woK0_add = new TCanvas("cacc_woK0_add","acc_woK0_add");
  q_IMnpipi_wSid_n_acc_add->Draw("colz");

  //w/o K0 S+/S- selection
  TH2F* q_IMnpipi_woK0_wSid_n_Sp_acc_add = (TH2F*)q_IMnpipi_woK0_wSid_n_Sp_acc_nocut->Clone();
  q_IMnpipi_woK0_wSid_n_Sp_acc_add->Add(q_IMnpipi_woK0_wSid_n_Sp_acc_cut1400);
  TCanvas *cacc_woK0_Sp_add = new TCanvas("cacc_woK0_Sp_add","acc_woK0_Sp_add");
  q_IMnpipi_woK0_wSid_n_Sp_acc_add->Draw("colz");
  

  TH2F* React_q_IMPiSigma_Sp_nocut = (TH2F*)file1->Get("React_q_IMPiSigma_Sp");
  TCanvas *cgen_Sp_nocut = new TCanvas("cgen_Sp_nocut","gen_Sp_nocut");
  React_q_IMPiSigma_Sp_nocut->Draw("colz");

  TH2F* React_q_IMPiSigma_Sp_cut1400 = (TH2F*)file2->Get("React_q_IMPiSigma_Sp");
  TCanvas *cgen_Sp_cut1400 = new TCanvas("cgen_Sp_cut1400","gen_Sp_cut1400");
  React_q_IMPiSigma_Sp_cut1400->Draw("colz");

  TH2F* React_q_IMPiSigma_Sp_add = (TH2F*)React_q_IMPiSigma_Sp_nocut->Clone();
  React_q_IMPiSigma_Sp_add->Add(React_q_IMPiSigma_Sp_cut1400);
  TCanvas *cgen_Sp_add = new TCanvas("cgen_Sp_add","gen_Sp_add");
  React_q_IMPiSigma_Sp_add->Draw("colz");
  
  /*
  TEfficiency *pEff;
  pEff = new TEfficiency(*q_IMnpipi_woK0_wSid_n_Sp_acc_add,*React_q_IMPiSigma_Sp_add);
  pEff->SetTitle("eff Sp test");
  pEff->SetName("heff_Sp_add");
  TCanvas *ceff = new TCanvas("ceff","ceff");
  ceff->cd();
  pEff->Draw("EY");
  ceff->Update();
  //TCanvas *ceff_hist = new TCanvas("ceff_hist","ceff_hist");
  auto heff = pEff->GetPaintedHistogram();
  heff->SetMinimum(0.0);
  heff->SetMaximum(0.008);
  ceff->Update();
  heff->Draw("colz");
  heff->Print("base");
  ceff->Modified();
  */

  TH2F* heff_add = new TH2F("heff_add","heff_add",500,1,2,25,0,1.5);
  TCanvas *ceff_hist = new TCanvas("ceff_hist","ceff_hist");
  heff_add->RebinX(4);
  q_IMnpipi_woK0_wSid_n_Sp_acc_add->RebinX(4);
  React_q_IMPiSigma_Sp_add->RebinX(4);

  heff_add->Divide(q_IMnpipi_woK0_wSid_n_Sp_acc_add,React_q_IMPiSigma_Sp_add,1.0,1.0,"b");
  heff_add->Draw("colz");

  TH2F* heff_err_add = new TH2F("heff_err_add","heff_err_add",500,1,2,25,0,1.5);
  heff_err_add->RebinX(4);
  for(int ix=0;ix<heff_add->GetNbinsX();ix++){
    for(int iy=0;iy<heff_add->GetNbinsY();iy++){
      double err = heff_add->GetBinErrorUp(ix,iy);
      double cont = heff_add->GetBinContent(ix,iy);
      if(cont>0.00)heff_err_add->SetBinContent(ix,iy,err/cont);
    }
  }
  heff_err_add->SetMaximum(1.0);
  TCanvas *cacc_err = new TCanvas("cacc_err","acc_err");
  cacc_err->cd();
  heff_err_add->Draw("colz");

  



}
