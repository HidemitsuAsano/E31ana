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

  fLpim = TFile::Open("simIMLpim_ppimpim_v17_out.root");
  fgen = TFile::Open("simIMLpim_v17.root");
  TFile *fkin = TFile::Open("NumericalRootFinderLPim.root","READ");
//  fLpim = TFile::Open("simIMLpim_ppimpim_pS0pim_v1_out.root");
//  fgen = TFile::Open("simIMLpim_pS0pim_v1.root");
  
  TH2F* q_IMLPim_gen;
  q_IMLPim_gen = (TH2F*) fgen->Get("React_q_IMLPim");
//  q_IMLPim_gen = (TH2F*) fgen->Get("React_q_IMS0Pim");
  q_IMLPim_gen->SetTitle("generated evt. ");
  q_IMLPim_gen->SetXTitle("true IM(#pi^{-}#Lambda) [GeV/c^{2}]");
  q_IMLPim_gen->SetYTitle("true Mom. Transfer [GeV/c]");
  q_IMLPim_gen->GetXaxis()->CenterTitle();
  q_IMLPim_gen->GetYaxis()->CenterTitle();
  q_IMLPim_gen->Print("base");
  q_IMLPim_gen->RebinX(15);
  q_IMLPim_gen->RebinY(3);
  q_IMLPim_gen->Print("base");
  
   
  TH2F* q_IMppipi_p_wL_sum;
  q_IMppipi_p_wL_sum= (TH2F*)fLpim->Get("q_IMppipi_p_wL_sum");
  q_IMppipi_p_wL_sum->GetYaxis()->SetRangeUser(0,1.5);
  q_IMppipi_p_wL_sum->SetTitle("reco. evt.");
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
  
  for(int ix=0;ix<q_IMppipi_p_wL_accerr->GetNbinsX();ix++){
    for(int iy=0;iy<q_IMppipi_p_wL_accerr->GetNbinsY();iy++){
      double err = q_IMppipi_p_wL_accerr->GetBinContent(ix,iy);
      if( RemoveNotEnough && (err>UncertCut)){
        q_IMppipi_p_wL_acc->SetBinContent(ix,iy,0);
        q_IMppipi_p_wL_acc->SetBinError(ix,iy,0);
      }
    }
  }
  
  TCanvas *cLpim;
  cLpim = new TCanvas(Form("cLpim"),Form("cLpim"),1300,1000);
  cLpim->Divide(2,2);
  cLpim->cd(1);
  q_IMLPim_gen->Draw("colz");
  TGraph *gth = (TGraph*)fkin->Get("th");
  gth->Draw("pc");
  TGraph *gr_0 = (TGraph*)fkin->Get("gr_0");
  gr_0->Draw("pc");
  TGraph *gr_100 = (TGraph*)fkin->Get("gr_100");
  gr_100->Draw("pc");
  TGraph *gr_65 = (TGraph*)fkin->Get("gr_65");
  gr_65->Draw("pc");
  cLpim->cd(2);
  q_IMppipi_p_wL_sum->Draw("colz");
  gth->Draw("pc");
  gr_0->Draw("pc");
  gr_100->Draw("pc");
  gr_65->Draw("pc");
  cLpim->cd(3);
  q_IMppipi_p_wL_acc->SetMaximum(0.05);
  q_IMppipi_p_wL_acc->Draw("colz");
  gth->Draw("pc");
  gr_0->Draw("pc");
  gr_100->Draw("pc");
  gr_65->Draw("pc");
  cLpim->cd(4);
  q_IMppipi_p_wL_accerr->SetMaximum(0.5);
  q_IMppipi_p_wL_accerr->Draw("colz");
    
 
//  TFile *fout = new TFile("accmapLpim_pS0pim.root","RECREATE");
  TFile *fout = new TFile("accmapLpimv17.root","RECREATE");
  q_IMppipi_p_wL_acc->Write();
  q_IMppipi_p_wL_accerr->Write();
  fout->Close();
}
