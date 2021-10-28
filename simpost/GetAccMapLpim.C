const bool RemoveNotEnough = true;
const double UncertCut = 0.20;

void GetAccMapLpim()
{
   
  gStyle->SetPalette(1);
  gStyle->SetOptStat("ei");
  gStyle->SetStatX(0.9);
  gStyle->SetStatY(0.9);
  TH1::SetDefaultSumw2(true);
  TFile *fLpim=NULL;
  TFile *fgen=NULL;

  fLpim = TFile::Open("simIMLpim_ppimpim_v21_out.root");
  fgen = TFile::Open("simIMLpim_v21.root");
  TFile *fkin = TFile::Open("NumericalRootFinderLPim.root","READ");
//  fLpim = TFile::Open("simIMLpim_ppimpim_pS0pim_v1_out.root");
//  fgen = TFile::Open("simIMLpim_pS0pim_v1.root");
  
  TH2F* q_IMLPim_gen=NULL;
  q_IMLPim_gen = (TH2F*) fgen->Get("React_q_IMLPim");
//  q_IMLPim_gen = (TH2F*) fgen->Get("React_q_IMS0Pim");
  q_IMLPim_gen->SetTitle("generated evt. ");
  q_IMLPim_gen->SetXTitle("true IM(#pi^{-}#Lambda) [GeV/c^{2}]");
  q_IMLPim_gen->SetYTitle("true Mom. Transfer [GeV/c]");
  q_IMLPim_gen->GetXaxis()->CenterTitle();
  q_IMLPim_gen->GetYaxis()->CenterTitle();
  q_IMLPim_gen->Print("base");
  q_IMLPim_gen->RebinX(15);
  q_IMLPim_gen->RebinY(6);
  q_IMLPim_gen->Print("base");
   

   
  TH2F* q_IMppipi_p_wL_sum;
  q_IMppipi_p_wL_sum= (TH2F*)fLpim->Get("q_IMppipi_p_wL_sum");
  q_IMppipi_p_wL_sum->GetYaxis()->SetRangeUser(0,1.5);
  q_IMppipi_p_wL_sum->SetTitle("reco. evt.");
  TH1F* BLAnaPassed = (TH1F*)fgen->Get("BLAnaPassed");
  double SimBeamSurvivalOK = BLAnaPassed->GetBinContent(2);//passed
  double SimBeamSurvivalFail = BLAnaPassed->GetBinContent(1);//not passed
  double SimBeamSurvivalRate = SimBeamSurvivalOK / (SimBeamSurvivalOK+SimBeamSurvivalFail);
  std::cout << "L " << __LINE__ << " SimBeamSurvivalRate " << SimBeamSurvivalRate << std::endl;
  q_IMppipi_p_wL_sum->Scale(1./SimBeamSurvivalRate);
  //q_IMppipi_p_wL_sum->RebinX(3);
  //q_IMppipi_p_wL_sum->RebinY(3);

  TH2F* q_IMppipi_p_wL_acc;
  q_IMppipi_p_wL_acc = (TH2F*)q_IMppipi_p_wL_sum->Clone(Form("q_IMppipi_p_wL_acc"));
  q_IMppipi_p_wL_acc->SetTitle(Form("q_IMppipi_p_wL acc. "));
  q_IMppipi_p_wL_acc->Print("base");
  q_IMppipi_p_wL_acc->Divide(q_IMppipi_p_wL_acc,q_IMLPim_gen,1.0,1.0,"b");

  TH2F* q_IMppipi_p_wL_accerr=NULL;
  q_IMppipi_p_wL_accerr = (TH2F*)q_IMppipi_p_wL_acc->Clone(Form("q_IMppipi_p_wL_accerr"));
  q_IMppipi_p_wL_accerr->Reset();
  q_IMppipi_p_wL_accerr->SetTitle(Form("q_IMppipi_p_wL precision"));
  for(int ix=0;ix<q_IMppipi_p_wL_accerr->GetNbinsX();ix++){
    for(int iy=0;iy<q_IMppipi_p_wL_accerr->GetNbinsY();iy++){
      double cont = q_IMppipi_p_wL_acc->GetBinContent(ix,iy);
      double err = q_IMppipi_p_wL_acc->GetBinError(ix,iy);
      if(cont!=0){
        q_IMppipi_p_wL_accerr->SetBinContent(ix,iy,err/cont);  
      }
    }
  }
  
  
  TCanvas *cLpim;
  cLpim = new TCanvas(Form("cLpim"),Form("cLpim"),2000,1200);
  cLpim->Divide(2,2);
  cLpim->cd(1);
  q_IMLPim_gen->Draw("colz");
 
  TGraph *gth = (TGraph*)fkin->Get("th");
  TGraph *gr_0 = (TGraph*)fkin->Get("gr_0");
  TGraph *gr_100 = (TGraph*)fkin->Get("gr_100");
  TGraph *gr_65 = (TGraph*)fkin->Get("gr_65");
 
  //gth->Draw("pc");
  //gr_0->Draw("pc");
  //gr_65->Draw("pc");
  //gr_100->Draw("pc");
  cLpim->cd(2);
  q_IMppipi_p_wL_sum->Draw("colz");
  //gth->Draw("pc");
  //gr_0->Draw("pc");
  //gr_100->Draw("pc");
  //gr_65->Draw("pc");
  cLpim->cd(3);
  q_IMppipi_p_wL_acc->SetMaximum(0.055);
  q_IMppipi_p_wL_acc->Draw("colz");
  //gth->Draw("pc");
  //gr_0->Draw("pc");
  //gr_100->Draw("pc");
  //gr_65->Draw("pc");
  cLpim->cd(4);
  q_IMppipi_p_wL_accerr->SetMaximum(0.5);
  q_IMppipi_p_wL_accerr->Draw("colz");
   
  //test for reco.events/acc map 
  TCanvas *cLpim_gentest = new TCanvas("cLpim_gentest","cLpim_gentest",2000,800);
  cLpim_gentest->Divide(2,1);
  cLpim_gentest->cd(1);
  TH2F* q_IMppipi_p_wL_gentest = (TH2F*)q_IMppipi_p_wL_sum->Clone("q_IMppipi_p_wL_gentest");
  q_IMppipi_p_wL_gentest->Divide(q_IMppipi_p_wL_acc);
  q_IMppipi_p_wL_gentest->SetTitle("TEST reco/acc map");
  q_IMppipi_p_wL_gentest->Draw("colz");
  cLpim_gentest->cd(2);
  q_IMLPim_gen->Draw("colz");
   
  TH1F* htest = new TH1F("htest","htest",8000,-6000,2000);
  for(int ix=0;ix<q_IMppipi_p_wL_gentest->GetNbinsX();ix++){
    for(int iy=0;iy<q_IMppipi_p_wL_gentest->GetNbinsY();iy++){
      double cont = q_IMppipi_p_wL_gentest->GetBinContent(ix,iy);
      double gen = q_IMLPim_gen->GetBinContent(ix,iy);
      htest->Fill(cont-gen);
    }
  }
  TCanvas *ctest  = new TCanvas("ctest","ctest",1000,800);
  htest->Draw("HE");

  TCanvas *ccleanup = new TCanvas("ccleanup","ccleanup",1000,800);
  TH2F* q_IMppipi_p_wL_gentest_clean = (TH2F*)q_IMppipi_p_wL_gentest->Clone("q_IMppipi_p_wL_gentest_clean");
  TH2F* q_IMppipi_p_wL_acc_clean = (TH2F*)q_IMppipi_p_wL_acc->Clone("q_IMppipi_p_wL_acc_clean");
  q_IMppipi_p_wL_acc_clean->SetName("q_IMppipi_p_wL_acc_clean");
  q_IMppipi_p_wL_acc_clean->SetTitle("q_IMppipi_p_wL_acc_clean");
  for(int ix=0;ix<q_IMppipi_p_wL_gentest->GetNbinsX();ix++){
    for(int iy=0;iy<q_IMppipi_p_wL_gentest->GetNbinsY();iy++){
      double cont = q_IMppipi_p_wL_gentest->GetBinContent(ix,iy);
      double gen = q_IMLPim_gen->GetBinContent(ix,iy);
      if(cont -gen < -100 ){
        q_IMppipi_p_wL_gentest_clean->SetBinContent(ix,iy,0);
        q_IMppipi_p_wL_gentest_clean->SetBinError(ix,iy,0);
        q_IMppipi_p_wL_acc_clean->SetBinContent(ix,iy,0);
        q_IMppipi_p_wL_acc_clean->SetBinError(ix,iy,0);
      }
    }
  }
  
  
  for(int ix=0;ix<q_IMppipi_p_wL_accerr->GetNbinsX();ix++){
    for(int iy=0;iy<q_IMppipi_p_wL_accerr->GetNbinsY();iy++){
      double err = q_IMppipi_p_wL_accerr->GetBinContent(ix,iy);
      if( RemoveNotEnough && (err>UncertCut) 
        //  || q_IMppipi_p_wL_acc->GetBinContent(ix,iy)>0.053
        //  || q_IMppipi_p_wL_acc->GetXaxis()->GetBinCenter(ix)>1.895 
          || q_IMppipi_p_wL_acc->GetXaxis()->GetBinCenter(ix)<1.260) {
        q_IMppipi_p_wL_acc_clean->SetBinContent(ix,iy,0);
        q_IMppipi_p_wL_acc_clean->SetBinError(ix,iy,0);
      }
    }
  }
  q_IMppipi_p_wL_gentest_clean->SetMaximum(q_IMppipi_p_wL_gentest->GetMaximum());
  q_IMppipi_p_wL_gentest_clean->Draw("colz");

  TCanvas *cacccleanup = new TCanvas("cacccleanup","cacccleanup",1000,800);
  q_IMppipi_p_wL_acc_clean->SetMaximum(q_IMppipi_p_wL_acc->GetMaximum());
  q_IMppipi_p_wL_acc_clean->Draw("colz");
  gr_0->Draw("pc");
  gr_100->Draw("pc");
  gr_65->Draw("pc");
  gth->Draw("pc");
   
  //////////////////////////////////////////////////////////
  //acceptance map for costheta_missing p v.s. IM(pi-Lambda) 
  //////////////////////////////////////////////////////////
  TH2F* costhetap_IMLpim_gen=NULL;
  costhetap_IMLpim_gen = (TH2F*)fgen->Get("React_costhetap_IMLPim");
  costhetap_IMLpim_gen->SetTitle("generated evt. ");
  costhetap_IMLpim_gen->SetXTitle("true IM(#pi^{-}#Lambda) [GeV/c^{2}]");
  costhetap_IMLpim_gen->SetYTitle("true proton CosTheta");
  costhetap_IMLpim_gen->GetXaxis()->CenterTitle();
  costhetap_IMLpim_gen->GetYaxis()->CenterTitle();
  costhetap_IMLpim_gen->Print("base");
  costhetap_IMLpim_gen->RebinX(10);
  costhetap_IMLpim_gen->Print("base");
 
  TH2F* CosTheta_IMppipi_p_wL_sum=NULL;
  CosTheta_IMppipi_p_wL_sum = (TH2F*)fLpim->Get("CosTheta_IMppipi_p_wL_sum");
  CosTheta_IMppipi_p_wL_sum ->SetTitle("reco. evt.");
  CosTheta_IMppipi_p_wL_sum->Scale(1./SimBeamSurvivalRate);

  TH2F* CosTheta_IMppipi_p_wL_acc=NULL;
  CosTheta_IMppipi_p_wL_acc = (TH2F*)CosTheta_IMppipi_p_wL_sum->Clone("CosTheta_IMppipi_p_wL_acc");
  CosTheta_IMppipi_p_wL_acc->SetTitle("CosTheta_IMppipi_p_wL_acc");
  CosTheta_IMppipi_p_wL_acc->Print("base");
  CosTheta_IMppipi_p_wL_acc->Divide(CosTheta_IMppipi_p_wL_acc,costhetap_IMLpim_gen,1.0,1.0,"b");


  TH2F* CosTheta_IMppipi_p_wL_accerr=NULL;
  CosTheta_IMppipi_p_wL_accerr = (TH2F*)CosTheta_IMppipi_p_wL_acc->Clone("CosTheta_IMppipi_p_wL_accerr");
  CosTheta_IMppipi_p_wL_accerr->Reset();
  CosTheta_IMppipi_p_wL_accerr->SetTitle("CosTheta_IMppipi_p_wL_acc precision");
  
  for(int ix=0;ix<CosTheta_IMppipi_p_wL_accerr->GetNbinsX();ix++){
    for(int iy=0;iy<CosTheta_IMppipi_p_wL_accerr->GetNbinsY();iy++){
      double cont = CosTheta_IMppipi_p_wL_acc->GetBinContent(ix,iy);
      double err = CosTheta_IMppipi_p_wL_acc->GetBinError(ix,iy);
      if(cont!=0){
        CosTheta_IMppipi_p_wL_accerr->SetBinContent(ix,iy,err/cont);  
      }
    }
  }
  


  TCanvas *cLpimCos;
  cLpimCos = new TCanvas(Form("cLpimCos"),Form("cLpimCos"),2000,1200);
  cLpimCos->Divide(2,2);
  cLpimCos->cd(1);
  costhetap_IMLpim_gen->GetYaxis()->SetRangeUser(0,1);
  costhetap_IMLpim_gen->Draw("colz");
  gPad->SetLogz();
  cLpimCos->cd(2);
  CosTheta_IMppipi_p_wL_sum->GetYaxis()->SetRangeUser(0,1);
  CosTheta_IMppipi_p_wL_sum->Draw("colz");
  gPad->SetLogz();
  cLpimCos->cd(3);
  //CosThetaIMppipi_p_wL_acc->SetMaximum(0.05);
  CosTheta_IMppipi_p_wL_acc->GetYaxis()->SetRangeUser(0,1);
  CosTheta_IMppipi_p_wL_acc->SetMaximum(0.05);
  CosTheta_IMppipi_p_wL_acc->Draw("colz");
 // gPad->SetLogz();
  cLpimCos->cd(4);
  //CosThetaIMppipi_p_wL_accerr->SetMaximum(0.5);
  CosTheta_IMppipi_p_wL_accerr->GetYaxis()->SetRangeUser(0,1);
  CosTheta_IMppipi_p_wL_accerr->Draw("colz");
  
  //TCanvas *cLpimCosTest = new TCanvas("cLpimCosTest","cLpimCosTest",1000,800);
  //TH2F* CosTheta_IMppipi_p_wL_test = (TH2F*)CosTheta_IMppipi_p_wL_sum->Clone("CosTheta_IMppipi_p_wL_test");
  //CosTheta_IMppipi_p_wL_test->Divide(CosTheta_IMppipi_p_wL_acc);
  //CosTheta_IMppipi_p_wL_test->Draw("colz");
  
  /*
  TCanvas *cLpimCosTestDiff = new TCanvas("cLpimCosTestDiff","cLpimCosTestDiff",1000,800);
  TH1D* hLpimCosTestDiff = new TH1D("hLpimCosTestDiff","hLpimCosTestDiff",10000,-0.0001,0.0001);
  for(int ix=0;ix<CosTheta_IMppipi_p_wL_test->GetNbinsX();ix++){
    for(int iy=0;iy<CosTheta_IMppipi_p_wL_test->GetNbinsY();iy++){
      double conttest = CosTheta_IMppipi_p_wL_test->GetBinContent(ix,iy);
      double contgen = costhetap_IMLpim_gen->GetBinContent(ix,iy);
      if(contgen!=0){
        hLpimCosTestDiff->Fill((conttest-contgen)/contgen*100.0);
        if(fabs((conttest-contgen)/contgen*100.0)>0.02e-03){
          CosTheta_IMppipi_p_wL_acc_clean->SetBinContent(ix,iy,0);
          CosTheta_IMppipi_p_wL_acc_clean->SetBinError(ix,iy,0);
          //std::cout << "ix:iy" << ix << " " << iy << std::endl;
        }
      }
    }
  }
  hLpimCosTestDiff->Draw("HE");
  */
  

  TH2F* CosTheta_IMppipi_p_wL_acc_clean = (TH2F*)CosTheta_IMppipi_p_wL_acc->Clone("CosTheta_IMppipi_p_wL_acc_clean");
  for(int ix=0;ix<CosTheta_IMppipi_p_wL_acc->GetNbinsX();ix++){
    for(int iy=0;iy<CosTheta_IMppipi_p_wL_acc->GetNbinsY();iy++){
      double error = CosTheta_IMppipi_p_wL_accerr->GetBinContent(ix,iy);
      if(error>0.30){
        CosTheta_IMppipi_p_wL_acc_clean->SetBinContent(ix,iy,0);
        CosTheta_IMppipi_p_wL_acc_clean->SetBinError(ix,iy,0);
      }
    }
  }

  TCanvas *cLpimCosClean = new TCanvas("cLpimCosClean","cLpimCosClean",1000,800);
  CosTheta_IMppipi_p_wL_acc_clean->Draw("colz");
        



  //  TFile *fout = new TFile("accmapLpim_pS0pim.root","RECREATE");
  TFile *fout = new TFile("accmapLpimv21.root","RECREATE");
  q_IMppipi_p_wL_acc->Write();
  q_IMppipi_p_wL_accerr->Write();
  q_IMppipi_p_wL_acc_clean->Write();
  CosTheta_IMppipi_p_wL_acc->Write();
  CosTheta_IMppipi_p_wL_acc_clean->Write();
  fout->Close();
}
