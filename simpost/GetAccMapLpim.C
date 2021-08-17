const bool RemoveNotEnough = false;
const double UncertCut = 0.25;

void GetAccMapLpim()
{
   
  gStyle->SetPalette(1);
  gStyle->SetOptStat("ei");
  gStyle->SetStatX(0.9);
  gStyle->SetStatY(0.9);
  TH1::SetDefaultSumw2(true);
  TFile *fLpim=NULL;
  TFile *fgen=NULL;

  fLpim = TFile::Open("simIMLpim_ppimpim_v7_out.root");
  fgen = TFile::Open("simIMLpim_ppimpim_v7.root");
  
  TH2F* q_IMLPim_gen;
  q_IMLPim_gen = (TH2F*) fgen->Get("React_q_IMLPim");
  q_IMLPim_gen->SetTitle("generated evt. ");
  q_IMLPim_gen->SetXTitle("true IM(#pi^{-}#Lambda) [GeV/c^{2}]");
  q_IMLPim_gen->SetYTitle("true Mom. Transfer [GeV/c]");
  q_IMLPim_gen->GetXaxis()->CenterTitle();
  q_IMLPim_gen->GetYaxis()->CenterTitle();
  q_IMLPim_gen->print("base");
  q_IMLPim_gen->RebinX(15);
  q_IMLPim_gen->RebinY(3);
  q_IMLPim_gen->print("base");
  
   
  TH2F* q_IMppipi_p_wL_sum;
  q_IMppipi_p_wL_sum= (TH2F*)fLpim->Get("q_IMppipi_p_wL_sum");
  q_IMppipi_p_wL_sum->GetYaxis()->SetRangeUser(0,1.5);
  q_IMppipi_p_wL_sum->SetTitle("reco. evt.");
  q_IMppipi_p_wL_sum->RebinX(3);
  //q_IMppipi_p_wL_sum->RebinY(3);

  TH2F* q_IMppipi_p_wL_acc;
  q_IMppipi_p_wL_acc = (TH2F*)q_IMppipi_p_wL_sum->Clone(Form("q_IMppipi_p_wL_acc"));
  q_IMppipi_p_wL_acc->SetTitle(Form("q_IMppipi_p_wL acc. "));
  q_IMppipi_p_wL_acc->Print("base");
  q_IMppipi_p_wL_acc->Divide(q_IMppipi_p_wL_acc,q_IMLPim_gen,1.0,1.0,"b");

  TH2F* q_IMppipi_p_wL_accerr;
  q_IMnpipi_K0_accerr[iq] = (TH2F*)q_IMnpipi_K0_acc[0]->Clone(Form("q_IMnpipi_K0_accerr_%d",iq));
  q_IMnpipi_K0_accerr[iq]->Reset();
  q_IMnpipi_K0_accerr[iq]->SetTitle(Form("q_IMnpipi_K0 precision",iq));
  for(int ix=0;ix<q_IMnpipi_Sp_acc[0]->GetNbinsX();ix++){
    for(int iy=0;iy<q_IMnpipi_Sp_acc[0]->GetNbinsY();iy++){
      double contSp = q_IMnpipi_Sp_acc[iq]->GetBinContent(ix,iy);
      double errSp = q_IMnpipi_Sp_acc[iq]->GetBinError(ix,iy);
      if(contSp!=0){
        q_IMnpipi_Sp_accerr[iq]->SetBinContent(ix,iy,errSp/contSp);  
      }
      double contSm = q_IMnpipi_Sm_acc[iq]->GetBinContent(ix,iy);
      double errSm = q_IMnpipi_Sm_acc[iq]->GetBinError(ix,iy);
      if(contSm!=0){
        q_IMnpipi_Sm_accerr[iq]->SetBinContent(ix,iy,errSm/contSm);  
      }
      double contK0 = q_IMnpipi_K0_acc[iq]->GetBinContent(ix,iy);
      double errK0 = q_IMnpipi_K0_acc[iq]->GetBinError(ix,iy);
      if(contK0!=0){
        q_IMnpipi_K0_accerr[iq]->SetBinContent(ix,iy,errK0/contK0);  
      }
    }
  }
  
  for(int ix=0;ix<q_IMnpipi_Sp_accerr[0]->GetNbinsX();ix++){
    for(int iy=0;iy<q_IMnpipi_Sp_accerr[0]->GetNbinsY();iy++){
      double err = q_IMnpipi_Sp_accerr[0]->GetBinContent(ix,iy);
      if( RemoveNotEnough && (err>UncertCut)){
        q_IMnpipi_Sp_acc[0]->SetBinContent(ix,iy,0);
        q_IMnpipi_Sp_acc[0]->SetBinError(ix,iy,0);
      }
    }
  }
  
  for(int ix=0;ix<q_IMnpipi_Sm_accerr[0]->GetNbinsX();ix++){
    for(int iy=0;iy<q_IMnpipi_Sm_accerr[0]->GetNbinsY();iy++){
      double err = q_IMnpipi_Sm_accerr[0]->GetBinContent(ix,iy);
      if(RemoveNotEnough && err>UncertCut){
        q_IMnpipi_Sm_acc[0]->SetBinContent(ix,iy,0);
        q_IMnpipi_Sm_acc[0]->SetBinError(ix,iy,0);
      }
    }
  }

  TCanvas *cSp[nqcut];
  TCanvas *cSm[nqcut];
  TCanvas *cK0[nqcut];
  for(int iq=0;iq<1;iq++){
    cSp[iq] = new TCanvas(Form("cSp%d",iq),Form("cSp%d",iq),1300,1000);
    cSp[iq]->Divide(2,2);
    cSp[iq]->cd(1);
    q_IMnpipi_gen_Sp[iq]->Draw("colz");
    cSp[iq]->cd(2);
    q_IMnpipi_wSid_n_Sp_reco[iq]->Draw("colz");
    cSp[iq]->cd(3);
    q_IMnpipi_Sp_acc[iq]->SetMaximum(0.005);
    q_IMnpipi_Sp_acc[iq]->Draw("colz");
    cSp[iq]->cd(4);
    q_IMnpipi_Sp_accerr[iq]->SetMaximum(0.5);
    q_IMnpipi_Sp_accerr[iq]->Draw("colz");
    
    cSm[iq] = new TCanvas(Form("cSm%d",iq),Form("cSm%d",iq),1300,1000);
    cSm[iq]->Divide(2,2);
    cSm[iq]->cd(1);
    q_IMnpipi_gen_Sm[iq]->Draw("colz");
    cSm[iq]->cd(2);
    q_IMnpipi_wSid_n_Sm_reco[iq]->Draw("colz");
    cSm[iq]->cd(3);
    q_IMnpipi_Sm_acc[iq]->SetMaximum(0.010);
    q_IMnpipi_Sm_acc[iq]->Draw("colz");
    cSm[iq]->cd(4);
    q_IMnpipi_Sm_accerr[iq]->SetMaximum(0.5);
    q_IMnpipi_Sm_accerr[iq]->Draw("colz");
    
    cK0[iq] = new TCanvas(Form("cK0%d",iq),Form("cK0%d",iq),1300,1000);
    cK0[iq]->Divide(2,2);
    cK0[iq]->cd(1);
    q_IMnpipi_gen_K0[iq]->Draw("colz");
    cK0[iq]->cd(2);
    q_IMnpipi_wK0_n_K0_reco[iq]->Draw("colz");
    cK0[iq]->cd(3);
    q_IMnpipi_K0_acc[iq]->SetMaximum(0.005);
    q_IMnpipi_K0_acc[iq]->Draw("colz");
    cK0[iq]->cd(4);
    q_IMnpipi_K0_accerr[iq]->SetMaximum(0.5);
    q_IMnpipi_K0_accerr[iq]->Draw("colz");
  }
 
  TFile *fout = new TFile("accmap.root","RECREATE");
  for(int iq=0;iq<nqcut;iq++){
    q_IMnpipi_Sp_acc[iq]->Write();
    q_IMnpipi_Sp_accerr[iq]->Write();
    q_IMnpipi_Sm_acc[iq]->Write();
    q_IMnpipi_Sm_accerr[iq]->Write();
    q_IMnpipi_K0_acc[iq]->Write();
    q_IMnpipi_K0_accerr[iq]->Write();
  }
  fout->Close();
}
