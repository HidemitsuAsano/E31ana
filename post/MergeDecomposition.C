void MergeDecomposition()
{
  TFile *fr[4] = {NULL};
  TFile *fmix[4] = {NULL};

  fr[0] = TFile::Open("evanaIMpisigma_npippim_v202_out_iso.root","READ");
  fmix[0] = TFile::Open("evanaIMpisigma_npippim_v202_MIX_cut4_out_iso.root","READ");
  fr[1] = TFile::Open("evanaIMpisigma_npippim_v202_out_iso_qlo.root","READ");
  fmix[1] = TFile::Open("evanaIMpisigma_npippim_v202_MIX_cut4_out_iso_qlo.root","READ");
  fr[2] = TFile::Open("evanaIMpisigma_npippim_v202_out_iso_qhi.root","READ");
  fmix[2] = TFile::Open("evanaIMpisigma_npippim_v202_MIX_cut4_out_iso_qhi.root","READ");
  fr[3] = TFile::Open("evanaIMpisigma_npippim_v202_out_iso_theta15.root","READ");
  fmix[3] = TFile::Open("evanaIMpisigma_npippim_v202_MIX_cut4_out_iso_theta15.root","READ");
  fr[0]->Print();
  fmix[0]->Print();
  fr[1]->Print();
  fmix[1]->Print();
  fr[2]->Print();
  fmix[2]->Print();
  fr[3]->Print();
  fmix[3]->Print();

  TFile *fdecompos = TFile::Open("K0SigmaTemp.root","READ");
  fdecompos->Print();

  const int nqcut=4;
  const char cqcut[][10]= {"all","qlo","qhi","theta"};
  TH2D* q_IMnpipi_woK0_wSid_n_woSm_data[nqcut];
  TH2D* q_IMnpipi_woK0_wSid_n_woSm_mix[nqcut];
  TH2D* q_IMnpipi_woK0_wSid_n_woSm_sub[nqcut];
  TH1D* IMnpipi_woK0_wSid_n_woSm_sub[nqcut];
  TH1D* IMnpipi_woK0_wSid_n_woSm_sub_merge[nqcut];
  TH2D* q_IMnpipi_woK0_wSid_n_woSp_data[nqcut];
  TH2D* q_IMnpipi_woK0_wSid_n_woSp_mix[nqcut];
  TH2D* q_IMnpipi_woK0_wSid_n_woSp_sub[nqcut];
  TH1D* IMnpipi_woK0_wSid_n_woSp_sub[nqcut];
  TH1D* IMnpipi_woK0_wSid_n_woSp_sub_merge[nqcut];
  TH2D* q_IMnpipi_wK0_woSid_n_data[nqcut];
  TH2D* q_IMnpipi_wK0_woSid_n_mix[nqcut];
  TH2D* q_IMnpipi_wK0_woSid_n_sub[nqcut];
  TH1D* IMnpipi_wK0_woSid_n_sub[nqcut];
  TH1D* IMnpipi_wK0_woSid_n_sub_merge[nqcut];
  TH1D* IMnpipi_wK0_wSid_n_SpSm[nqcut];//for double-counting subtraction
  //checking histogram
  TH2D* q_IMnpipi_wSid_n_Sp_data[nqcut];
  TH2D* q_IMnpipi_wSid_n_Sp_mix[nqcut];
  TH2D* q_IMnpipi_wSid_n_Sp_sub[nqcut];
  TH2D* q_IMnpipi_wK0_wSid_n_Sp_data[nqcut];
  TH2D* q_IMnpipi_wK0_wSid_n_Sp_mix[nqcut];
  TH2D* q_IMnpipi_wK0_wSid_n_Sp_sub[nqcut];
  TH2D* q_IMnpipi_wSid_n_Sm_data[nqcut];
  TH2D* q_IMnpipi_wSid_n_Sm_mix[nqcut];
  TH2D* q_IMnpipi_wSid_n_Sm_sub[nqcut];
  TH2D* q_IMnpipi_wK0_wSid_n_Sm_data[nqcut];
  TH2D* q_IMnpipi_wK0_wSid_n_Sm_mix[nqcut];
  TH2D* q_IMnpipi_wK0_wSid_n_Sm_sub[nqcut];
  TH2D* q_IMnpipi_wSid_n_SpSm_data[nqcut];
  TH2D* q_IMnpipi_wSid_n_SpSm_mix[nqcut];
  TH2D* q_IMnpipi_wSid_n_SpSm_sub[nqcut];
  TH2D* q_IMnpipi_wK0_n_data[nqcut];
  TH2D* q_IMnpipi_wK0_n_mix[nqcut];
  TH2D* q_IMnpipi_wK0_n_sub[nqcut];



  for(int iqcut=0;iqcut<nqcut;iqcut++){
    q_IMnpipi_woK0_wSid_n_woSm_data[iqcut] = (TH2D*)fr[iqcut]->Get("q_IMnpipi_woK0_wSid_n_woSm");
    q_IMnpipi_woK0_wSid_n_woSm_mix[iqcut] = (TH2D*)fmix[iqcut]->Get("q_IMnpipi_woK0_wSid_n_woSm");
    q_IMnpipi_woK0_wSid_n_woSm_sub[iqcut] = (TH2D*)q_IMnpipi_woK0_wSid_n_woSm_data[iqcut]->Clone(Form("q_IMnpipi_woK0_wSid_n_woSm_sub_%s",cqcut[iqcut]));
    q_IMnpipi_woK0_wSid_n_woSm_sub[iqcut]->Add(q_IMnpipi_woK0_wSid_n_woSm_mix[iqcut],-1.0);

    q_IMnpipi_woK0_wSid_n_woSp_data[iqcut] = (TH2D*)fr[iqcut]->Get("q_IMnpipi_woK0_wSid_n_woSp");
    q_IMnpipi_woK0_wSid_n_woSp_mix[iqcut] = (TH2D*)fmix[iqcut]->Get("q_IMnpipi_woK0_wSid_n_woSp");
    q_IMnpipi_woK0_wSid_n_woSp_sub[iqcut] = (TH2D*)q_IMnpipi_woK0_wSid_n_woSp_data[iqcut]->Clone(Form("q_IMnpipi_woK0_wSid_n_woSp_sub_%s",cqcut[iqcut]));
    q_IMnpipi_woK0_wSid_n_woSp_sub[iqcut]->Add(q_IMnpipi_woK0_wSid_n_woSp_mix[iqcut],-1.0);

    q_IMnpipi_wK0_woSid_n_data[iqcut] = (TH2D*)fr[iqcut]->Get("q_IMnpipi_wK0_woSid_n");
    q_IMnpipi_wK0_woSid_n_mix[iqcut] = (TH2D*)fmix[iqcut]->Get("q_IMnpipi_wK0_woSid_n");
    q_IMnpipi_wK0_woSid_n_sub[iqcut] = (TH2D*)q_IMnpipi_wK0_woSid_n_data[iqcut]->Clone(Form("q_IMnpipi_wK0_woSid_n_sub_%s",cqcut[iqcut]));
    q_IMnpipi_wK0_woSid_n_sub[iqcut]->Add(q_IMnpipi_wK0_woSid_n_mix[iqcut],-1.0);
  
    q_IMnpipi_wSid_n_Sp_data[iqcut] = (TH2D*)fr[iqcut]->Get("q_IMnpipi_wSid_n_Sp");
    q_IMnpipi_wSid_n_Sp_mix[iqcut] = (TH2D*)fmix[iqcut]->Get("q_IMnpipi_wSid_n_Sp");
    q_IMnpipi_wSid_n_Sp_sub[iqcut] = (TH2D*)q_IMnpipi_wSid_n_Sp_data[iqcut]->Clone(Form("q_IMnpipi_wSid_n_Sp_sub_%s",cqcut[iqcut]));
    q_IMnpipi_wSid_n_Sp_sub[iqcut]->Add(q_IMnpipi_wSid_n_Sp_mix[iqcut],-1.0);
  
    q_IMnpipi_wK0_wSid_n_Sp_data[iqcut] = (TH2D*)fr[iqcut]->Get("q_IMnpipi_wK0_wSid_n_Sp");
    q_IMnpipi_wK0_wSid_n_Sp_mix[iqcut] = (TH2D*)fmix[iqcut]->Get("q_IMnpipi_wK0_wSid_n_Sp");
    q_IMnpipi_wK0_wSid_n_Sp_sub[iqcut] = (TH2D*)q_IMnpipi_wK0_wSid_n_Sp_data[iqcut]->Clone(Form("q_IMnpipi_wK0_wSid_n_Sp_sub_%s",cqcut[iqcut]));
    q_IMnpipi_wK0_wSid_n_Sp_sub[iqcut]->Add(q_IMnpipi_wK0_wSid_n_Sp_mix[iqcut],-1.0);
  
    q_IMnpipi_wSid_n_Sm_data[iqcut] = (TH2D*)fr[iqcut]->Get("q_IMnpipi_wSid_n_Sm");
    q_IMnpipi_wSid_n_Sm_mix[iqcut] = (TH2D*)fmix[iqcut]->Get("q_IMnpipi_wSid_n_Sm");
    q_IMnpipi_wSid_n_Sm_sub[iqcut] = (TH2D*)q_IMnpipi_wSid_n_Sm_data[iqcut]->Clone(Form("q_IMnpipi_wSid_n_Sm_sub_%s",cqcut[iqcut]));
    q_IMnpipi_wSid_n_Sm_sub[iqcut]->Add(q_IMnpipi_wSid_n_Sm_mix[iqcut],-1.0);
  
    q_IMnpipi_wK0_wSid_n_Sm_data[iqcut] = (TH2D*)fr[iqcut]->Get("q_IMnpipi_wK0_wSid_n_Sm");
    q_IMnpipi_wK0_wSid_n_Sm_mix[iqcut] = (TH2D*)fmix[iqcut]->Get("q_IMnpipi_wK0_wSid_n_Sm");
    q_IMnpipi_wK0_wSid_n_Sm_sub[iqcut] = (TH2D*)q_IMnpipi_wK0_wSid_n_Sm_data[iqcut]->Clone(Form("q_IMnpipi_wK0_wSid_n_Sm_sub_%s",cqcut[iqcut]));
    q_IMnpipi_wK0_wSid_n_Sm_sub[iqcut]->Add(q_IMnpipi_wK0_wSid_n_Sm_mix[iqcut],-1.0);
  
    q_IMnpipi_wSid_n_SpSm_data[iqcut] = (TH2D*)fr[iqcut]->Get("q_IMnpipi_wSid_n_SpSm");
    q_IMnpipi_wSid_n_SpSm_mix[iqcut] = (TH2D*)fmix[iqcut]->Get("q_IMnpipi_wSid_n_SpSm");
    q_IMnpipi_wSid_n_SpSm_sub[iqcut] = (TH2D*)q_IMnpipi_wSid_n_SpSm_data[iqcut]->Clone(Form("q_IMnpipi_wSid_n_SpSm_sub_%s",cqcut[iqcut]));
    q_IMnpipi_wSid_n_SpSm_sub[iqcut]->Add(q_IMnpipi_wSid_n_SpSm_mix[iqcut],-1.0);
  
    q_IMnpipi_wK0_n_data[iqcut] = (TH2D*)fr[iqcut]->Get("q_IMnpipi_wK0_n");
    q_IMnpipi_wK0_n_mix[iqcut] = (TH2D*)fmix[iqcut]->Get("q_IMnpipi_wK0_n");
    q_IMnpipi_wK0_n_sub[iqcut] = (TH2D*)q_IMnpipi_wK0_n_data[iqcut]->Clone(Form("q_IMnpipi_wK0_n_sub_%s",cqcut[iqcut]));
    q_IMnpipi_wK0_n_sub[iqcut]->Add(q_IMnpipi_wK0_n_mix[iqcut],-1.0);
  }
  
  TCanvas *cOverlapCheck[3][nqcut];
  TH1D* IMnpipi_check_Sp[nqcut];
  TH1D* IMnpipi_check_SpSm[nqcut];
  TH1D* IMnpipi_check_wK0_Sp[nqcut];
  TH1D* IMnpipi_check_Sm[nqcut];
  TH1D* IMnpipi_check_wK0_Sm[nqcut];
  TH1D* IMnpipi_check_K0[nqcut];
  for(int iq=0;iq<nqcut;iq++){
    cOverlapCheck[0][iq] = new TCanvas(Form("cOverlapCheck_Sp_%d",iq),Form("cOverlapCheck_Sp_%d",iq));
    cOverlapCheck[0][iq]->cd();
    IMnpipi_check_Sp[iq] = (TH1D*)q_IMnpipi_wSid_n_Sp_sub[iq]->ProjectionX(Form("IMnpipi_check_Sp_%d",iq));
    IMnpipi_check_Sp[iq]->Draw("HE");
    IMnpipi_check_SpSm[iq] = (TH1D*)q_IMnpipi_wSid_n_SpSm_sub[iq]->ProjectionX(Form("IMnpipi_check_SpSm_%d",iq));
    IMnpipi_check_SpSm[iq]->SetLineColor(2);
    IMnpipi_check_SpSm[iq]->Draw("HEsame");
    IMnpipi_check_wK0_Sp[iq] = (TH1D*)q_IMnpipi_wK0_wSid_n_Sp_sub[iq]->ProjectionX(Form("IMnpipi_check_wK0_Sp_%d",iq));
    IMnpipi_check_wK0_Sp[iq]->SetLineColor(3);
    IMnpipi_check_wK0_Sp[iq]->Draw("HEsame");
    cOverlapCheck[1][iq] = new TCanvas(Form("cOverlapCheck_Sm_%d",iq),Form("cOverlapCheck_Sm_%d",iq));
    cOverlapCheck[1][iq]->cd();
    IMnpipi_check_Sm[iq] = (TH1D*)q_IMnpipi_wSid_n_Sm_sub[iq]->ProjectionX(Form("IMnpipi_check_Sm_%d",iq));
    IMnpipi_check_Sm[iq]->Draw("HE");
    IMnpipi_check_SpSm[iq]->Draw("HEsame");
    IMnpipi_check_wK0_Sm[iq] = (TH1D*)q_IMnpipi_wK0_wSid_n_Sm_sub[iq]->ProjectionX(Form("IMnpipi_check_wK0_Sm_%d",iq));
    IMnpipi_check_wK0_Sm[iq]->SetLineColor(3);
    IMnpipi_check_wK0_Sm[iq]->Draw("HEsame");
    //cOverlapCheck[2][iq] = new TCanvas(Form("cOverlapCheck_K0_%d",iq),Form("cOverlapCheck_K0_%d",iq));
    //cOverlapCheck[2][iq]->cd();
    //IMnpipi_check_K0[iq] = (TH1D*)q_IMnpipi_wK0_n_sub[iq]->ProjectionX(Form("IMnpipi_check_K0_%d",iq));
    //IMnpipi_check_K0[iq]->Draw("HE");
    //IMnpipi_check_wK0_Sp[iq]->Draw("HEsame");
    //IMnpipi_check_wK0_Sm[iq]->Draw("HEsame");
  }

  
  return;
  const int ngroup = 3;
  TH1D* IMnpipi_overlapdeco_K0[ngroup][nqcut];
  TH1D* IMnpipi_overlapdeco_Sp[ngroup][nqcut];
  TH1D* IMnpipi_overlapdeco_Sm[ngroup][nqcut];
  
  for(int ig=0;ig<ngroup;ig++){
    for(int iq=0;iq<nqcut;iq++){
      std::cout << ig << " " << iq << std::endl;
      if(ig!=2)IMnpipi_overlapdeco_K0[ig][iq] = (TH1D*)fdecompos->Get(Form("IMnpipi_overlapdeco_K0_g%d_%d",ig,iq));
      if(ig!=1)IMnpipi_overlapdeco_Sp[ig][iq] = (TH1D*)fdecompos->Get(Form("IMnpipi_overlapdeco_Sp_g%d_%d",ig,iq));
      if(ig!=0)IMnpipi_overlapdeco_Sm[ig][iq] = (TH1D*)fdecompos->Get(Form("IMnpipi_overlapdeco_Sm_g%d_%d",ig,iq));
    }
  }

  for(int iq=0;iq<nqcut;iq++){
    IMnpipi_wK0_wSid_n_SpSm[iq] = (TH1D*)fdecompos->Get(Form("IMnpipi_wK0_wSid_n_SpSm_%d",iq));
  }


  
  TH1D* IMnpipi_overlapToK0[nqcut];
  TH1D* IMnpipi_overlapToSp[nqcut];
  TH1D* IMnpipi_overlapToSm[nqcut];

  for(int iqcut=0;iqcut<nqcut;iqcut++){
    IMnpipi_woK0_wSid_n_woSm_sub[iqcut] = (TH1D*)q_IMnpipi_woK0_wSid_n_woSm_sub[iqcut]->ProjectionX(Form("IMnpipi_woK0_wSid_n_woSm_sub_%d",iqcut));
    IMnpipi_woK0_wSid_n_woSm_sub_merge[iqcut] = (TH1D*)IMnpipi_woK0_wSid_n_woSm_sub[iqcut]->Clone(Form("IMnpipi_woK0_wSid_n_woSm_sub_merge_%d",iqcut));
    IMnpipi_overlapToSp[iqcut] = (TH1D*)IMnpipi_overlapdeco_Sp[0][iqcut]->Clone(Form("IMnpipi_overlapToSp_%d",iqcut));
    IMnpipi_overlapToSp[iqcut]->Add(IMnpipi_overlapdeco_Sp[2][iqcut]);
    //subtract double-counting 
    IMnpipi_overlapToSp[iqcut]->Add(IMnpipi_wK0_wSid_n_SpSm[iqcut],-1);
    IMnpipi_woK0_wSid_n_woSm_sub_merge[iqcut]->Add(IMnpipi_overlapToSp[iqcut]);

    IMnpipi_woK0_wSid_n_woSp_sub[iqcut] = (TH1D*)q_IMnpipi_woK0_wSid_n_woSp_sub[iqcut]->ProjectionX(Form("IMnpipi_woK0_wSid_n_woSp_sub_%d",iqcut));
    IMnpipi_woK0_wSid_n_woSp_sub_merge[iqcut] = (TH1D*)IMnpipi_woK0_wSid_n_woSp_sub[iqcut]->Clone(Form("IMnpipi_woK0_wSid_n_woSp_sub_merge_%d",iqcut));
    IMnpipi_overlapToSm[iqcut] = (TH1D*)IMnpipi_overlapdeco_Sm[1][iqcut]->Clone(Form("IMnpipi_overlapToSm_%d",iqcut));
    IMnpipi_overlapToSm[iqcut]->Add(IMnpipi_overlapdeco_Sm[2][iqcut]); 
    //subtract double-counting 
    IMnpipi_overlapToSm[iqcut]->Add(IMnpipi_wK0_wSid_n_SpSm[iqcut],-1);
    IMnpipi_woK0_wSid_n_woSp_sub_merge[iqcut]->Add(IMnpipi_overlapToSm[iqcut]);
    IMnpipi_wK0_woSid_n_sub[iqcut] = (TH1D*)q_IMnpipi_wK0_woSid_n_sub[iqcut]->ProjectionX(Form("IMnpipi_wK0_woSid_n_sub_%d",iqcut));
    IMnpipi_wK0_woSid_n_sub_merge[iqcut] = (TH1D*)IMnpipi_wK0_woSid_n_sub[iqcut]->Clone(Form("IMnpipi_wK0_woSid_n_sub_merge_%d",iqcut));
    IMnpipi_overlapToK0[iqcut] = (TH1D*)IMnpipi_overlapdeco_K0[0][iqcut]->Clone(Form("IMnpipi_overlapToK0_%d",iqcut));
    IMnpipi_overlapToK0[iqcut]->Add(IMnpipi_overlapdeco_K0[1][iqcut]); 
    //subtract double-counting 
    IMnpipi_overlapToK0[iqcut]->Add(IMnpipi_wK0_wSid_n_SpSm[iqcut],-1);
    IMnpipi_wK0_woSid_n_sub_merge[iqcut]->Add(IMnpipi_overlapToK0[iqcut]);
  }
  
  TCanvas *cmerge_Sp[nqcut];
  for(int iqcut=0;iqcut<nqcut;iqcut++){
    cmerge_Sp[iqcut] = new TCanvas(Form("cmerge_Sp%d",iqcut),Form("cmerge_Sp%d",iqcut));
    cmerge_Sp[iqcut]->cd();
    IMnpipi_woK0_wSid_n_woSm_sub_merge[iqcut]->SetLineColor(3);
    IMnpipi_woK0_wSid_n_woSm_sub_merge[iqcut]->SetMarkerColor(3);
    IMnpipi_woK0_wSid_n_woSm_sub_merge[iqcut]->SetMarkerStyle(20);
    IMnpipi_woK0_wSid_n_woSm_sub_merge[iqcut]->Draw("E");
    IMnpipi_woK0_wSid_n_woSm_sub[iqcut]->SetMarkerStyle(20);
    IMnpipi_woK0_wSid_n_woSm_sub[iqcut]->Draw("Esame");
    IMnpipi_overlapToSp[iqcut]->SetLineColor(2);
    IMnpipi_overlapToSp[iqcut]->SetMarkerStyle(20);
    IMnpipi_overlapToSp[iqcut]->SetMarkerColor(2);
    IMnpipi_overlapToSp[iqcut]->Draw("Esame");
  }
  
  TCanvas *cmerge_Sm[nqcut];
  for(int iqcut=0;iqcut<nqcut;iqcut++){
    cmerge_Sm[iqcut] = new TCanvas(Form("cmerge_Sm%d",iqcut),Form("cmerge_Sm%d",iqcut));
    cmerge_Sm[iqcut]->cd();
    IMnpipi_woK0_wSid_n_woSp_sub_merge[iqcut]->SetLineColor(3);
    IMnpipi_woK0_wSid_n_woSp_sub_merge[iqcut]->SetMarkerColor(3);
    IMnpipi_woK0_wSid_n_woSp_sub_merge[iqcut]->SetMarkerStyle(20);
    IMnpipi_woK0_wSid_n_woSp_sub_merge[iqcut]->Draw("E");
    IMnpipi_woK0_wSid_n_woSp_sub[iqcut]->SetMarkerStyle(20);
    IMnpipi_woK0_wSid_n_woSp_sub[iqcut]->Draw("Esame");
    IMnpipi_overlapToSm[iqcut]->SetLineColor(2);
    IMnpipi_overlapToSm[iqcut]->SetMarkerColor(2);
    IMnpipi_overlapToSm[iqcut]->SetMarkerStyle(20);
    IMnpipi_overlapToSm[iqcut]->Draw("Esame");
  }

  TCanvas *cmerge_K0[nqcut];
  for(int iqcut=0;iqcut<nqcut;iqcut++){
    cmerge_K0[iqcut] = new TCanvas(Form("cmerge_K0%d",iqcut),Form("cmerge_K0%d",iqcut));
    cmerge_K0[iqcut]->cd();
    IMnpipi_wK0_woSid_n_sub_merge[iqcut]->SetLineColor(3);
    IMnpipi_wK0_woSid_n_sub_merge[iqcut]->SetMarkerColor(3);
    IMnpipi_wK0_woSid_n_sub_merge[iqcut]->SetMarkerStyle(20);
    IMnpipi_wK0_woSid_n_sub_merge[iqcut]->Draw("E");
    IMnpipi_wK0_woSid_n_sub[iqcut]->SetMarkerStyle(20);
    IMnpipi_wK0_woSid_n_sub[iqcut]->Draw("Esame");
    IMnpipi_overlapToK0[iqcut]->SetLineColor(2);
    IMnpipi_overlapToK0[iqcut]->SetMarkerColor(2);
    IMnpipi_overlapToK0[iqcut]->SetMarkerStyle(20);
    IMnpipi_overlapToK0[iqcut]->Draw("Esame");
  }
  

}
