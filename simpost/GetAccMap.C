const bool RemoveNotEnough = true;
const double UncertCut = 0.25;

void GetAccMap()
{
   
  gStyle->SetPalette(1);
  gStyle->SetOptStat("ei");
  gStyle->SetStatX(0.9);
  gStyle->SetStatY(0.9);
  TH1::SetDefaultSumw2(true);
  TFile *fSp[4]={NULL};
  TFile *fSm[4]={NULL};
  TFile *fK0[4]={NULL};
  TFile *fSpgen=NULL;
  TFile *fSmgen=NULL;
  TFile *fK0gen=NULL;
  const int versionSigma = 152;
  const int versionK0 = 27;
  const int dEcut=2;

  fSp[0] = TFile::Open(Form("simIMpisigma_nSppim_pippimn_v%d_out_dE%d_iso_rej_nostop.root",versionSigma,dEcut));
  fSp[1] = TFile::Open(Form("simIMpisigma_nSppim_pippimn_v%d_out_dE%d_iso_qlo_rej_nostop.root",versionSigma,dEcut));
  fSp[2] = TFile::Open(Form("simIMpisigma_nSppim_pippimn_v%d_out_dE%d_iso_qhi_rej_nostop.root",versionSigma,dEcut));
  fSp[3] = TFile::Open(Form("simIMpisigma_nSppim_pippimn_v%d_out_dE%d_iso_theta15_rej_nostop.root",versionSigma,dEcut));

  fSm[0] = TFile::Open(Form("simIMpisigma_nSmpip_pippimn_v%d_out_dE%d_iso_rej_nostop.root",versionSigma,dEcut));
  fSm[1] = TFile::Open(Form("simIMpisigma_nSmpip_pippimn_v%d_out_dE%d_iso_qlo_rej_nostop.root",versionSigma,dEcut));
  fSm[2] = TFile::Open(Form("simIMpisigma_nSmpip_pippimn_v%d_out_dE%d_iso_qhi_rej_nostop.root",versionSigma,dEcut));
  fSm[3] = TFile::Open(Form("simIMpisigma_nSmpip_pippimn_v%d_out_dE%d_iso_theta15_rej_nostop.root",versionSigma,dEcut));

  fK0[0] = TFile::Open(Form("simIMpisigma_K0nn_pippimn_v%d_out_dE%d_iso_rej_nostop.root",versionK0,dEcut));
  fK0[1] = TFile::Open(Form("simIMpisigma_K0nn_pippimn_v%d_out_dE%d_iso_qlo_rej_nostop.root",versionK0,dEcut));
  fK0[2] = TFile::Open(Form("simIMpisigma_K0nn_pippimn_v%d_out_dE%d_iso_qhi_rej_nostop.root",versionK0,dEcut));
  fK0[3] = TFile::Open(Form("simIMpisigma_K0nn_pippimn_v%d_out_dE%d_iso_theta15_rej_nostop.root",versionK0,dEcut));
  
  fSpgen = TFile::Open(Form("simIMpisigma_nSppim_v%d.root",versionSigma));
  fSmgen = TFile::Open(Form("simIMpisigma_nSmpip_v%d.root",versionSigma));
  fK0gen = TFile::Open(Form("simIMpisigma_K0nn_v%d.root",versionK0));
  
  const int nqcut=1;
  TH2F* q_IMnpipi_gen_Sp[nqcut];
  TH2F* q_IMnpipi_gen_Sm[nqcut];
  TH2F* q_IMnpipi_gen_K0[nqcut];
  for(int iq=0;iq<nqcut;iq++){
    q_IMnpipi_gen_Sp[iq] = (TH2F*) fSpgen->Get("React_q_IMPiSigma");
    q_IMnpipi_gen_Sp[iq]->SetTitle("generated evt. Sp");
    q_IMnpipi_gen_Sp[iq]->SetXTitle("true IM(#pi^{-}#Sigma^{+}) [GeV/c^{2}]");
    q_IMnpipi_gen_Sp[iq]->SetYTitle("true Mom. Transfer [GeV/c]");
    q_IMnpipi_gen_Sp[iq]->GetXaxis()->CenterTitle();
    q_IMnpipi_gen_Sp[iq]->GetYaxis()->CenterTitle();

    q_IMnpipi_gen_Sm[iq] = (TH2F*) fSmgen->Get("React_q_IMPiSigma");
    q_IMnpipi_gen_Sm[iq]->SetTitle("generated evt. Sm");
    q_IMnpipi_gen_Sm[iq]->SetXTitle("true IM(#pi^{+}#Sigma^{-}) [GeV/c^{2}]");
    q_IMnpipi_gen_Sm[iq]->SetYTitle("true Mom. Transfer [GeV/c]");
    q_IMnpipi_gen_Sm[iq]->GetXaxis()->CenterTitle();
    q_IMnpipi_gen_Sm[iq]->GetYaxis()->CenterTitle();
    
    q_IMnpipi_gen_K0[iq] = (TH2F*) fK0gen->Get("React_q_IMnpipi");
    q_IMnpipi_gen_K0[iq]->SetTitle("generated evt. K0");
    q_IMnpipi_gen_K0[iq]->SetXTitle("true IM(nK^{0}) [GeV/c^{2}]");
    q_IMnpipi_gen_K0[iq]->SetYTitle("true Mom. Transfer [GeV/c]");
    q_IMnpipi_gen_K0[iq]->GetXaxis()->CenterTitle();
    q_IMnpipi_gen_K0[iq]->GetYaxis()->CenterTitle();
  }
  //for some reason, Rebin iq = 0-4 histograms
  q_IMnpipi_gen_Sp[0]->RebinX(15);
  q_IMnpipi_gen_Sp[0]->RebinY(10);
  //q_IMnpipi_gen_Sp[0]->GetXaxis()->SetRangeUser(1.2,2.0);
  
  q_IMnpipi_gen_Sm[0]->RebinX(15);
  q_IMnpipi_gen_Sm[0]->RebinY(10);
  //q_IMnpipi_gen_Sm[0]->GetXaxis()->SetRangeUser(1.2,2.0);
  q_IMnpipi_gen_K0[0]->RebinX(15);
  q_IMnpipi_gen_K0[0]->RebinY(10);
  //q_IMnpipi_gen_K0[0]->GetXaxis()->SetRangeUser(1.2,2.0);
   
  TH2F* q_IMnpipi_wSid_n_Sp_reco[nqcut];
  TH2F* q_IMnpipi_wSid_n_Sm_reco[nqcut];
  TH2F* q_IMnpipi_wK0_n_K0_reco[nqcut];
  TH1F* BLAnaPassedSp = (TH1F*)fSpgen->Get("BLAnaPassed");
  double SimBeamSurvivalOKSp = BLAnaPassedSp->GetBinContent(2);//passed
  double SimBeamSurvivalFailSp = BLAnaPassedSp->GetBinContent(1);//not passed
  double SimBeamSurvivalRateSp = SimBeamSurvivalOKSp / (SimBeamSurvivalOKSp+SimBeamSurvivalFailSp);
  TH1F* BLAnaPassedSm = (TH1F*)fSmgen->Get("BLAnaPassed");
  double SimBeamSurvivalOKSm = BLAnaPassedSm->GetBinContent(2);//passed
  double SimBeamSurvivalFailSm = BLAnaPassedSm->GetBinContent(1);//not passed
  double SimBeamSurvivalRateSm = SimBeamSurvivalOKSm / (SimBeamSurvivalOKSm+SimBeamSurvivalFailSm);
  TH1F* BLAnaPassedK0 = (TH1F*)fK0gen->Get("BLAnaPassed");
  double SimBeamSurvivalOKK0 = BLAnaPassedK0->GetBinContent(2);//passed
  double SimBeamSurvivalFailK0 = BLAnaPassedK0->GetBinContent(1);//not passed
  double SimBeamSurvivalRateK0 = SimBeamSurvivalOKK0 / (SimBeamSurvivalOKK0+SimBeamSurvivalFailK0);
  for(int iq=0;iq<nqcut;iq++){
    q_IMnpipi_wSid_n_Sp_reco[iq] = (TH2F*)fSp[iq]->Get("q_IMnpipi_wSid_n_Sp");
    q_IMnpipi_wSid_n_Sp_reco[iq]->SetTitle("reco. evt. Sp");
    q_IMnpipi_wSid_n_Sm_reco[iq] = (TH2F*)fSm[iq]->Get("q_IMnpipi_wSid_n_Sm");
    q_IMnpipi_wSid_n_Sm_reco[iq]->SetTitle("reco. evt. Sm");
    q_IMnpipi_wK0_n_K0_reco[iq] = (TH2F*)fK0[iq]->Get("q_IMnpipi_wK0_n");
    q_IMnpipi_wK0_n_K0_reco[iq]->GetYaxis()->SetRangeUser(0,1.5);
    q_IMnpipi_wK0_n_K0_reco[iq]->SetTitle("reco. evt. K0");
    q_IMnpipi_wSid_n_Sp_reco[iq]->RebinX(3);
    q_IMnpipi_wSid_n_Sp_reco[iq]->Scale(1./SimBeamSurvivalRateSp);
    q_IMnpipi_wSid_n_Sm_reco[iq]->RebinX(3);
    q_IMnpipi_wSid_n_Sm_reco[iq]->Scale(1./SimBeamSurvivalRateSm);
    q_IMnpipi_wK0_n_K0_reco[iq]->RebinX(3);
    q_IMnpipi_wK0_n_K0_reco[iq]->Scale(1./SimBeamSurvivalRateK0);
  }
  

  TH2F* q_IMnpipi_Sp_acc[nqcut];
  TH2F* q_IMnpipi_Sm_acc[nqcut];
  TH2F* q_IMnpipi_K0_acc[nqcut];
  for(int iq=0;iq<nqcut;iq++){
    q_IMnpipi_Sp_acc[iq] = (TH2F*)q_IMnpipi_wSid_n_Sp_reco[iq]->Clone(Form("q_IMnpipi_Sp_accp_%d",iq));
    q_IMnpipi_Sm_acc[iq] = (TH2F*)q_IMnpipi_wSid_n_Sm_reco[iq]->Clone(Form("q_IMnpipi_Sm_accp_%d",iq));
    q_IMnpipi_K0_acc[iq] = (TH2F*)q_IMnpipi_wK0_n_K0_reco[iq]->Clone(Form("q_IMnpipi_K0_accp_%d",iq));
    q_IMnpipi_Sp_acc[iq]->SetTitle(Form("q_IMnpipi Sp acc. ",iq));
    q_IMnpipi_Sm_acc[iq]->SetTitle(Form("q_IMnpipi Sm acc. ",iq));
    q_IMnpipi_K0_acc[iq]->SetTitle(Form("q_IMnpipi K0 acc. ",iq));
    
    
    q_IMnpipi_Sp_acc[iq]->Print("base");
    q_IMnpipi_gen_Sp[iq]->Print("base");
    q_IMnpipi_Sp_acc[iq]->Divide(q_IMnpipi_Sp_acc[iq],q_IMnpipi_gen_Sp[iq],1.0,1.0,"b");
    q_IMnpipi_Sm_acc[iq]->Print("base");
    q_IMnpipi_gen_Sm[iq]->Print("base");
    q_IMnpipi_Sm_acc[iq]->Divide(q_IMnpipi_Sm_acc[iq],q_IMnpipi_gen_Sm[iq],1.0,1.0,"b");
    q_IMnpipi_K0_acc[iq]->Print("base");
    q_IMnpipi_gen_K0[iq]->Print("base");
    q_IMnpipi_K0_acc[iq]->Divide(q_IMnpipi_K0_acc[iq],q_IMnpipi_gen_K0[iq],1.0,1.0,"b");
  }

  TH2F* q_IMnpipi_Sp_accerr[nqcut];
  TH2F* q_IMnpipi_Sm_accerr[nqcut];
  TH2F* q_IMnpipi_K0_accerr[nqcut];
  for(int iq=0;iq<nqcut;iq++){
    q_IMnpipi_Sp_accerr[iq] = (TH2F*)q_IMnpipi_Sp_acc[0]->Clone(Form("q_IMnpipi_Sp_accerr_%d",iq));
    q_IMnpipi_Sm_accerr[iq] = (TH2F*)q_IMnpipi_Sm_acc[0]->Clone(Form("q_IMnpipi_Sm_accerr_%d",iq));
    q_IMnpipi_K0_accerr[iq] = (TH2F*)q_IMnpipi_K0_acc[0]->Clone(Form("q_IMnpipi_K0_accerr_%d",iq));
    q_IMnpipi_Sp_accerr[iq]->Reset();
    q_IMnpipi_Sm_accerr[iq]->Reset();
    q_IMnpipi_K0_accerr[iq]->Reset();
    q_IMnpipi_Sp_accerr[iq]->SetTitle(Form("q_IMnpipi_Sp precision",iq));
    q_IMnpipi_Sm_accerr[iq]->SetTitle(Form("q_IMnpipi_Sm precision",iq));
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
  
  for(int ix=0;ix<q_IMnpipi_K0_accerr[0]->GetNbinsX();ix++){
    for(int iy=0;iy<q_IMnpipi_K0_accerr[0]->GetNbinsY();iy++){
      double err = q_IMnpipi_K0_accerr[0]->GetBinContent(ix,iy);
      if(RemoveNotEnough && err>UncertCut){
        q_IMnpipi_K0_acc[0]->SetBinContent(ix,iy,0);
        q_IMnpipi_K0_acc[0]->SetBinError(ix,iy,0);
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
 
  TFile *fout = new TFile(Form("accmapv%d_%d.root",versionSigma,versionK0), "RECREATE");
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
