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

  TFile *fdecompos = TFile::Open("K0SigmaTemp_Rebin.root","READ");
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
  
  //overlap events
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
  
  //no K0 Sigma+/- events check
  TH2D* q_IMnpipi_woK0_woSid_n_data[nqcut];
  TH2D* q_IMnpipi_woK0_woSid_n_mix[nqcut];
  TH2D* q_IMnpipi_woK0_woSid_n_sub[nqcut];
  TH2D* IMnpim_IMnpip_woK0_woSid_n_data[nqcut];
  TH2D* IMnpim_IMnpip_woK0_woSid_n_mix[nqcut];
  TH2D* IMnpim_IMnpip_woK0_woSid_n_sub[nqcut];


  for(int iq=0;iq<nqcut;iq++){
    q_IMnpipi_woK0_wSid_n_woSm_data[iq] = (TH2D*)fr[iq]->Get("q_IMnpipi_woK0_wSid_n_woSm");
    q_IMnpipi_woK0_wSid_n_woSm_mix[iq] = (TH2D*)fmix[iq]->Get("q_IMnpipi_woK0_wSid_n_woSm");
    q_IMnpipi_woK0_wSid_n_woSm_sub[iq] = (TH2D*)q_IMnpipi_woK0_wSid_n_woSm_data[iq]->Clone(Form("q_IMnpipi_woK0_wSid_n_woSm_sub_%s",cqcut[iq]));
    q_IMnpipi_woK0_wSid_n_woSm_sub[iq]->Add(q_IMnpipi_woK0_wSid_n_woSm_mix[iq],-1.0);

    q_IMnpipi_woK0_wSid_n_woSp_data[iq] = (TH2D*)fr[iq]->Get("q_IMnpipi_woK0_wSid_n_woSp");
    q_IMnpipi_woK0_wSid_n_woSp_mix[iq] = (TH2D*)fmix[iq]->Get("q_IMnpipi_woK0_wSid_n_woSp");
    q_IMnpipi_woK0_wSid_n_woSp_sub[iq] = (TH2D*)q_IMnpipi_woK0_wSid_n_woSp_data[iq]->Clone(Form("q_IMnpipi_woK0_wSid_n_woSp_sub_%s",cqcut[iq]));
    q_IMnpipi_woK0_wSid_n_woSp_sub[iq]->Add(q_IMnpipi_woK0_wSid_n_woSp_mix[iq],-1.0);

    q_IMnpipi_wK0_woSid_n_data[iq] = (TH2D*)fr[iq]->Get("q_IMnpipi_wK0_woSid_n");
    q_IMnpipi_wK0_woSid_n_mix[iq] = (TH2D*)fmix[iq]->Get("q_IMnpipi_wK0_woSid_n");
    q_IMnpipi_wK0_woSid_n_sub[iq] = (TH2D*)q_IMnpipi_wK0_woSid_n_data[iq]->Clone(Form("q_IMnpipi_wK0_woSid_n_sub_%s",cqcut[iq]));
    q_IMnpipi_wK0_woSid_n_sub[iq]->Add(q_IMnpipi_wK0_woSid_n_mix[iq],-1.0);
  
    q_IMnpipi_wSid_n_Sp_data[iq] = (TH2D*)fr[iq]->Get("q_IMnpipi_wSid_n_Sp");
    q_IMnpipi_wSid_n_Sp_mix[iq] = (TH2D*)fmix[iq]->Get("q_IMnpipi_wSid_n_Sp");
    q_IMnpipi_wSid_n_Sp_sub[iq] = (TH2D*)q_IMnpipi_wSid_n_Sp_data[iq]->Clone(Form("q_IMnpipi_wSid_n_Sp_sub_%s",cqcut[iq]));
    q_IMnpipi_wSid_n_Sp_sub[iq]->Add(q_IMnpipi_wSid_n_Sp_mix[iq],-1.0);
  
    q_IMnpipi_wK0_wSid_n_Sp_data[iq] = (TH2D*)fr[iq]->Get("q_IMnpipi_wK0_wSid_n_Sp");
    q_IMnpipi_wK0_wSid_n_Sp_mix[iq] = (TH2D*)fmix[iq]->Get("q_IMnpipi_wK0_wSid_n_Sp");
    q_IMnpipi_wK0_wSid_n_Sp_sub[iq] = (TH2D*)q_IMnpipi_wK0_wSid_n_Sp_data[iq]->Clone(Form("q_IMnpipi_wK0_wSid_n_Sp_sub_%s",cqcut[iq]));
    q_IMnpipi_wK0_wSid_n_Sp_sub[iq]->Add(q_IMnpipi_wK0_wSid_n_Sp_mix[iq],-1.0);
  
    q_IMnpipi_wSid_n_Sm_data[iq] = (TH2D*)fr[iq]->Get("q_IMnpipi_wSid_n_Sm");
    q_IMnpipi_wSid_n_Sm_mix[iq] = (TH2D*)fmix[iq]->Get("q_IMnpipi_wSid_n_Sm");
    q_IMnpipi_wSid_n_Sm_sub[iq] = (TH2D*)q_IMnpipi_wSid_n_Sm_data[iq]->Clone(Form("q_IMnpipi_wSid_n_Sm_sub_%s",cqcut[iq]));
    q_IMnpipi_wSid_n_Sm_sub[iq]->Add(q_IMnpipi_wSid_n_Sm_mix[iq],-1.0);
  
    q_IMnpipi_wK0_wSid_n_Sm_data[iq] = (TH2D*)fr[iq]->Get("q_IMnpipi_wK0_wSid_n_Sm");
    q_IMnpipi_wK0_wSid_n_Sm_mix[iq] = (TH2D*)fmix[iq]->Get("q_IMnpipi_wK0_wSid_n_Sm");
    q_IMnpipi_wK0_wSid_n_Sm_sub[iq] = (TH2D*)q_IMnpipi_wK0_wSid_n_Sm_data[iq]->Clone(Form("q_IMnpipi_wK0_wSid_n_Sm_sub_%s",cqcut[iq]));
    q_IMnpipi_wK0_wSid_n_Sm_sub[iq]->Add(q_IMnpipi_wK0_wSid_n_Sm_mix[iq],-1.0);
  
    q_IMnpipi_wSid_n_SpSm_data[iq] = (TH2D*)fr[iq]->Get("q_IMnpipi_wSid_n_SpSm");
    q_IMnpipi_wSid_n_SpSm_mix[iq] = (TH2D*)fmix[iq]->Get("q_IMnpipi_wSid_n_SpSm");
    q_IMnpipi_wSid_n_SpSm_sub[iq] = (TH2D*)q_IMnpipi_wSid_n_SpSm_data[iq]->Clone(Form("q_IMnpipi_wSid_n_SpSm_sub_%s",cqcut[iq]));
    q_IMnpipi_wSid_n_SpSm_sub[iq]->Add(q_IMnpipi_wSid_n_SpSm_mix[iq],-1.0);
  
    q_IMnpipi_wK0_n_data[iq] = (TH2D*)fr[iq]->Get("q_IMnpipi_wK0_n");
    q_IMnpipi_wK0_n_mix[iq] = (TH2D*)fmix[iq]->Get("q_IMnpipi_wK0_n");
    q_IMnpipi_wK0_n_sub[iq] = (TH2D*)q_IMnpipi_wK0_n_data[iq]->Clone(Form("q_IMnpipi_wK0_n_sub_%s",cqcut[iq]));
    q_IMnpipi_wK0_n_sub[iq]->Add(q_IMnpipi_wK0_n_mix[iq],-1.0);
  
    q_IMnpipi_woK0_woSid_n_data[iq] = (TH2D*)fr[iq]->Get("q_IMnpipi_woK0_woSid_n");
    q_IMnpipi_woK0_woSid_n_mix[iq]  = (TH2D*)fmix[iq]->Get("q_IMnpipi_woK0_woSid_n");
    q_IMnpipi_woK0_woSid_n_sub[iq] = (TH2D*)q_IMnpipi_woK0_woSid_n_data[iq]->Clone(Form("q_IMnpipi_woK0_woSid_n_sub_%s",cqcut[iq]));
    q_IMnpipi_woK0_woSid_n_sub[iq]->Add(q_IMnpipi_woK0_woSid_n_mix[iq],-1.0);
    
    IMnpim_IMnpip_woK0_woSid_n_data[iq] = (TH2D*)fr[iq]->Get("IMnpim_IMnpip_dE_woK0_woSid_n");
    IMnpim_IMnpip_woK0_woSid_n_mix[iq]  = (TH2D*)fmix[iq]->Get("IMnpim_IMnpip_dE_woK0_woSid_n");
    IMnpim_IMnpip_woK0_woSid_n_sub[iq] = (TH2D*)IMnpim_IMnpip_woK0_woSid_n_data[iq]->Clone(Form("IMnpim_IMnpip_woK0_woSid_n_sub_%s",cqcut[iq]));
    IMnpim_IMnpip_woK0_woSid_n_sub[iq]->Add(IMnpim_IMnpip_woK0_woSid_n_mix[iq],-1.0);
  }
  
  TCanvas *cOverlapCheck[3][nqcut];
  TH1D* IMnpipi_overlapevt_Sp[nqcut];
  TH1D* IMnpipi_overlapevt_SpSm[nqcut];
  TH1D* IMnpipi_overlapevt_SpK0[nqcut];
  TH1D* IMnpipi_overlapevt_Sm[nqcut];
  TH1D* IMnpipi_overlapevt_SmK0[nqcut];
  TH1D* IMnpipi_overlapevt_K0[nqcut];
  for(int iq=0;iq<nqcut;iq++){
    cOverlapCheck[0][iq] = new TCanvas(Form("cOverlapCheck_Sp_%d",iq),Form("cOverlapCheck_Sp_%d",iq));
    cOverlapCheck[0][iq]->cd();
    IMnpipi_overlapevt_Sp[iq] = (TH1D*)q_IMnpipi_wSid_n_Sp_sub[iq]->ProjectionX(Form("IMnpipi_overlapevt_Sp_%d",iq));
    IMnpipi_overlapevt_Sp[iq]->Draw("HE");
    IMnpipi_overlapevt_SpSm[iq] = (TH1D*)q_IMnpipi_wSid_n_SpSm_sub[iq]->ProjectionX(Form("IMnpipi_overlapevt_SpSm_%d",iq));
    IMnpipi_overlapevt_SpSm[iq]->SetLineColor(2);
    IMnpipi_overlapevt_SpSm[iq]->Draw("HEsame");
    IMnpipi_overlapevt_SpK0[iq] = (TH1D*)q_IMnpipi_wK0_wSid_n_Sp_sub[iq]->ProjectionX(Form("IMnpipi_overlapevt_SpK0_%d",iq));
    IMnpipi_overlapevt_SpK0[iq]->SetLineColor(3);
    IMnpipi_overlapevt_SpK0[iq]->Draw("HEsame");
    cOverlapCheck[1][iq] = new TCanvas(Form("cOverlapCheck_Sm_%d",iq),Form("cOverlapCheck_Sm_%d",iq));
    cOverlapCheck[1][iq]->cd();
    IMnpipi_overlapevt_Sm[iq] = (TH1D*)q_IMnpipi_wSid_n_Sm_sub[iq]->ProjectionX(Form("IMnpipi_overlapevt_Sm_%d",iq));
    IMnpipi_overlapevt_Sm[iq]->Draw("HE");
    IMnpipi_overlapevt_SpSm[iq]->Draw("HEsame");
    IMnpipi_overlapevt_SmK0[iq] = (TH1D*)q_IMnpipi_wK0_wSid_n_Sm_sub[iq]->ProjectionX(Form("IMnpipi_overlapevt_SmK0_%d",iq));
    IMnpipi_overlapevt_SmK0[iq]->SetLineColor(3);
    IMnpipi_overlapevt_SmK0[iq]->Draw("HEsame");
    cOverlapCheck[2][iq] = new TCanvas(Form("cOverlapCheck_K0_%d",iq),Form("cOverlapCheck_K0_%d",iq));
    cOverlapCheck[2][iq]->cd();
    IMnpipi_overlapevt_K0[iq] = (TH1D*)q_IMnpipi_wK0_n_sub[iq]->ProjectionX(Form("IMnpipi_overlapevt_K0_%d",iq));
    IMnpipi_overlapevt_K0[iq]->Draw("HE");
    IMnpipi_overlapevt_SpK0[iq]->Draw("HEsame");
    IMnpipi_overlapevt_SmK0[iq]->Draw("HEsame");
  }
  
  TCanvas *cwoK0_woSid_ncheck[nqcut];
  TCanvas *cwoK0_woSid_ncheck2[nqcut];
  TH1D* IMnpipi_woK0_woSid_n[nqcut];
  for(iq=0;iq<nqcut;iq++){ 
   cwoK0_woSid_ncheck[iq] = new TCanvas(Form("cwoK0_woSid_ncheck_%d",iq),Form("cwoK0_woSid_ncheck_%d",iq));
   IMnpipi_woK0_woSid_n[iq] = (TH1D*)q_IMnpipi_woK0_woSid_n_sub[iq]->ProjectionX(Form("IMnpipi_woK0_woSid_n_%d",iq));
   IMnpipi_woK0_woSid_n[iq]->Draw("HE");
   cwoK0_woSid_ncheck2[iq] = new TCanvas(Form("cwoK0_woSid_ncheck2_%d",iq),Form("cwoK0_woSid_ncheck2_%d",iq));
   IMnpim_IMnpip_woK0_woSid_n_sub[iq]->Draw("colz");
  }

  
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

  for(int iq=0;iq<nqcut;iq++){
    IMnpipi_woK0_wSid_n_woSm_sub[iq] = (TH1D*)q_IMnpipi_woK0_wSid_n_woSm_sub[iq]->ProjectionX(Form("IMnpipi_woK0_wSid_n_woSm_sub_%d",iq));
    IMnpipi_woK0_wSid_n_woSm_sub_merge[iq] = (TH1D*)IMnpipi_woK0_wSid_n_woSm_sub[iq]->Clone(Form("IMnpipi_woK0_wSid_n_woSm_sub_merge_%d",iq));
    IMnpipi_overlapToSp[iq] = (TH1D*)IMnpipi_overlapdeco_Sp[0][iq]->Clone(Form("IMnpipi_overlapToSp_%d",iq));
    IMnpipi_overlapToSp[iq]->Add(IMnpipi_overlapdeco_Sp[2][iq]);
    //subtract double-counting 
    IMnpipi_overlapToSp[iq]->Add(IMnpipi_wK0_wSid_n_SpSm[iq],-1);
    IMnpipi_woK0_wSid_n_woSm_sub_merge[iq]->Add(IMnpipi_overlapToSp[iq]);

    IMnpipi_woK0_wSid_n_woSp_sub[iq] = (TH1D*)q_IMnpipi_woK0_wSid_n_woSp_sub[iq]->ProjectionX(Form("IMnpipi_woK0_wSid_n_woSp_sub_%d",iq));
    IMnpipi_woK0_wSid_n_woSp_sub_merge[iq] = (TH1D*)IMnpipi_woK0_wSid_n_woSp_sub[iq]->Clone(Form("IMnpipi_woK0_wSid_n_woSp_sub_merge_%d",iq));
    IMnpipi_overlapToSm[iq] = (TH1D*)IMnpipi_overlapdeco_Sm[1][iq]->Clone(Form("IMnpipi_overlapToSm_%d",iq));
    IMnpipi_overlapToSm[iq]->Add(IMnpipi_overlapdeco_Sm[2][iq]); 
    //subtract double-counting 
    IMnpipi_overlapToSm[iq]->Add(IMnpipi_wK0_wSid_n_SpSm[iq],-1);
    IMnpipi_woK0_wSid_n_woSp_sub_merge[iq]->Add(IMnpipi_overlapToSm[iq]);
    IMnpipi_wK0_woSid_n_sub[iq] = (TH1D*)q_IMnpipi_wK0_woSid_n_sub[iq]->ProjectionX(Form("IMnpipi_wK0_woSid_n_sub_%d",iq));
    IMnpipi_wK0_woSid_n_sub_merge[iq] = (TH1D*)IMnpipi_wK0_woSid_n_sub[iq]->Clone(Form("IMnpipi_wK0_woSid_n_sub_merge_%d",iq));
    IMnpipi_overlapToK0[iq] = (TH1D*)IMnpipi_overlapdeco_K0[0][iq]->Clone(Form("IMnpipi_overlapToK0_%d",iq));
    IMnpipi_overlapToK0[iq]->Add(IMnpipi_overlapdeco_K0[1][iq]); 
    //subtract double-counting 
    IMnpipi_overlapToK0[iq]->Add(IMnpipi_wK0_wSid_n_SpSm[iq],-1);
    IMnpipi_wK0_woSid_n_sub_merge[iq]->Add(IMnpipi_overlapToK0[iq]);
  }
  
  /*
  TCanvas *cmerge_Sp[nqcut];
  for(int iq=0;iq<nqcut;iq++){
    cmerge_Sp[iq] = new TCanvas(Form("cmerge_Sp%d",iq),Form("cmerge_Sp%d",iq));
    cmerge_Sp[iq]->cd();
    IMnpipi_woK0_wSid_n_woSm_sub_merge[iq]->SetLineColor(3);
    IMnpipi_woK0_wSid_n_woSm_sub_merge[iq]->SetMarkerColor(3);
    IMnpipi_woK0_wSid_n_woSm_sub_merge[iq]->SetMarkerStyle(20);
    IMnpipi_woK0_wSid_n_woSm_sub_merge[iq]->Draw("E");
    IMnpipi_woK0_wSid_n_woSm_sub[iq]->SetMarkerStyle(20);
    IMnpipi_woK0_wSid_n_woSm_sub[iq]->Draw("Esame");
    IMnpipi_overlapToSp[iq]->SetLineColor(2);
    IMnpipi_overlapToSp[iq]->SetMarkerStyle(20);
    IMnpipi_overlapToSp[iq]->SetMarkerColor(2);
    IMnpipi_overlapToSp[iq]->Draw("Esame");
  }
  
  TCanvas *cmerge_Sm[nqcut];
  for(int iq=0;iq<nqcut;iq++){
    cmerge_Sm[iq] = new TCanvas(Form("cmerge_Sm%d",iq),Form("cmerge_Sm%d",iq));
    cmerge_Sm[iq]->cd();
    IMnpipi_woK0_wSid_n_woSp_sub_merge[iq]->SetLineColor(3);
    IMnpipi_woK0_wSid_n_woSp_sub_merge[iq]->SetMarkerColor(3);
    IMnpipi_woK0_wSid_n_woSp_sub_merge[iq]->SetMarkerStyle(20);
    IMnpipi_woK0_wSid_n_woSp_sub_merge[iq]->Draw("E");
    IMnpipi_woK0_wSid_n_woSp_sub[iq]->SetMarkerStyle(20);
    IMnpipi_woK0_wSid_n_woSp_sub[iq]->Draw("Esame");
    IMnpipi_overlapToSm[iq]->SetLineColor(2);
    IMnpipi_overlapToSm[iq]->SetMarkerColor(2);
    IMnpipi_overlapToSm[iq]->SetMarkerStyle(20);
    IMnpipi_overlapToSm[iq]->Draw("Esame");
  }

  TCanvas *cmerge_K0[nqcut];
  for(int iq=0;iq<nqcut;iq++){
    cmerge_K0[iq] = new TCanvas(Form("cmerge_K0%d",iq),Form("cmerge_K0%d",iq));
    cmerge_K0[iq]->cd();
    IMnpipi_wK0_woSid_n_sub_merge[iq]->SetLineColor(3);
    IMnpipi_wK0_woSid_n_sub_merge[iq]->SetMarkerColor(3);
    IMnpipi_wK0_woSid_n_sub_merge[iq]->SetMarkerStyle(20);
    IMnpipi_wK0_woSid_n_sub_merge[iq]->Draw("E");
    IMnpipi_wK0_woSid_n_sub[iq]->SetMarkerStyle(20);
    IMnpipi_wK0_woSid_n_sub[iq]->Draw("Esame");
    IMnpipi_overlapToK0[iq]->SetLineColor(2);
    IMnpipi_overlapToK0[iq]->SetMarkerColor(2);
    IMnpipi_overlapToK0[iq]->SetMarkerStyle(20);
    IMnpipi_overlapToK0[iq]->Draw("Esame");
  }
 */
  const int initbin[10]={40,52};
  const int startbin[10]={41,53};
  const int endbin[10]={52,100};
  const float startval[10]={1.40,1.52};
  const float endval[10]={1.52,2.00};
  TH1D* IMnpipi_Sp_ratio_wK0_merge[nqcut];
  TH1D* IMnpipi_K0_ratio_wSp_merge[nqcut];
  TH1D* IMnpipi_Sm_ratio_wK0_merge[nqcut];
  TH1D* IMnpipi_K0_ratio_wSm_merge[nqcut];
  TH1D* IMnpipi_Sp_ratio_wSm_merge[nqcut];
  TH1D* IMnpipi_Sm_ratio_wSp_merge[nqcut];
  
  for(int iq=0;iq<nqcut;iq++){
    IMnpipi_Sp_ratio_wK0_merge[iq] = (TH1D*)fdecompos->Get(Form("IMnpipi_Sp_ratio_wK0_merge_%d",iq));
    IMnpipi_K0_ratio_wSp_merge[iq] = (TH1D*)fdecompos->Get(Form("IMnpipi_K0_ratio_wSp_merge_%d",iq));
    IMnpipi_Sm_ratio_wK0_merge[iq] = (TH1D*)fdecompos->Get(Form("IMnpipi_Sm_ratio_wK0_merge_%d",iq));
    IMnpipi_K0_ratio_wSm_merge[iq] = (TH1D*)fdecompos->Get(Form("IMnpipi_K0_ratio_wSm_merge_%d",iq));
    IMnpipi_Sp_ratio_wSm_merge[iq] = (TH1D*)fdecompos->Get(Form("IMnpipi_Sp_ratio_wSm_merge_%d",iq));
    IMnpipi_Sm_ratio_wSp_merge[iq] = (TH1D*)fdecompos->Get(Form("IMnpipi_Sm_ratio_wSp_merge_%d",iq));
  }
  
  TCanvas *cratio[nqcut][6];
  for(int iq=0;iq<nqcut;iq++){
    cratio[iq][0] = new TCanvas(Form("cratio_q%d_0",iq),Form("cratio_q%d_0",iq));
    IMnpipi_Sp_ratio_wK0_merge[iq]->Draw("HE");

    cratio[iq][1] = new TCanvas(Form("cratio_q%d_1",iq),Form("cratio_q%d_1",iq));
    IMnpipi_K0_ratio_wSp_merge[iq]->Draw("HE");

    cratio[iq][2] = new TCanvas(Form("cratio_q%d_2",iq),Form("cratio_q%d_2",iq));
    IMnpipi_Sm_ratio_wK0_merge[iq]->Draw("HE");

    cratio[iq][3] = new TCanvas(Form("cratio_q%d_3",iq),Form("cratio_q%d_3",iq));
    IMnpipi_K0_ratio_wSm_merge[iq]->Draw("HE");

    cratio[iq][4] = new TCanvas(Form("cratio_q%d_4",iq),Form("cratio_q%d_4",iq));
    IMnpipi_Sp_ratio_wSm_merge[iq]->Draw("HE");

    cratio[iq][5] = new TCanvas(Form("cratio_q%d_5",iq),Form("cratio_q%d_5",iq));
    IMnpipi_Sm_ratio_wSp_merge[iq]->Draw("HE");
  }


  //decomposition using larger bin
  TH1D* IMnpipi_overlapToK02[nqcut];
  TH1D* IMnpipi_overlapToSp2[nqcut];
  TH1D* IMnpipi_overlapToSm2[nqcut];
  TH1D* IMnpipi_ToSp_wK0[nqcut];
  TH1D* IMnpipi_ToK0_wSp[nqcut];
  TH1D* IMnpipi_ToSm_wK0[nqcut];
  TH1D* IMnpipi_ToK0_wSm[nqcut];
  TH1D* IMnpipi_ToSp_wSm[nqcut];
  TH1D* IMnpipi_ToSm_wSp[nqcut];
  TH1D* IMnpipi_Sp_merge2[nqcut];
  TH1D* IMnpipi_Sm_merge2[nqcut];
  TH1D* IMnpipi_K0_merge2[nqcut];
  for(int iq=0;iq<nqcut;iq++){
    //Sp
    IMnpipi_Sp_merge2[iq] = (TH1D*)IMnpipi_woK0_wSid_n_woSm_sub[iq]->Clone(Form("IMnpipi_Sp_merge2_%d",iq));
    IMnpipi_ToSp_wK0[iq] = (TH1D*)IMnpipi_overlapevt_SpK0[iq]->Clone(Form("IMnpipi_ToSp_wK0_%d",iq));
    for(int imerge=0;imerge<2;imerge++){
      for(int ibin=initbin[imerge];ibin<endbin[imerge];ibin++){
        double cont = IMnpipi_ToSp_wK0[iq]->GetBinContent(ibin);
        double ratio = IMnpipi_Sp_ratio_wK0_merge[iq]->GetBinContent(ibin);
        IMnpipi_ToSp_wK0[iq]->SetBinContent(ibin,cont*ratio);
      }
    }
    
    IMnpipi_ToSp_wSm[iq] = (TH1D*)IMnpipi_overlapevt_SpSm[iq]->Clone(Form("IMnpipi_ToSp_wSm_%d",iq));
    for(int imerge=0;imerge<2;imerge++){
      for(int ibin=initbin[imerge];ibin<endbin[imerge];ibin++){
        double cont = IMnpipi_ToSp_wSm[iq]->GetBinContent(ibin);
        double ratio = IMnpipi_Sp_ratio_wSm_merge[iq]->GetBinContent(ibin);
        IMnpipi_ToSp_wSm[iq]->SetBinContent(ibin,cont*ratio);
      }
    }
    IMnpipi_Sp_merge2[iq]->Add(IMnpipi_ToSp_wK0[iq],1.0);
    IMnpipi_Sp_merge2[iq]->Add(IMnpipi_ToSp_wSm[iq],1.0);
    //subtract double-counting 
    IMnpipi_Sp_merge2[iq]->Add(IMnpipi_wK0_wSid_n_SpSm[iq],-1.0);
    
    //Sm
    IMnpipi_Sm_merge2[iq] = (TH1D*)IMnpipi_woK0_wSid_n_woSp_sub[iq]->Clone(Form("IMnpipi_Sm_merge2_%d",iq));
    IMnpipi_ToSm_wK0[iq] = (TH1D*)IMnpipi_overlapevt_SmK0[iq]->Clone(Form("IMnpipi_ToSm_wK0_%d",iq));
    //IMnpipi_ToSm_wK0[iq]->Multiply(IMnpipi_Sm_ratio_wK0_merge[iq]);
    for(int imerge=0;imerge<2;imerge++){
      for(int ibin=initbin[imerge];ibin<endbin[imerge];ibin++){
        double cont = IMnpipi_ToSm_wK0[iq]->GetBinContent(ibin);
        double ratio = IMnpipi_Sm_ratio_wK0_merge[iq]->GetBinContent(ibin);
        IMnpipi_ToSm_wK0[iq]->SetBinContent(ibin,cont*ratio);
      }
    }
    IMnpipi_ToSm_wSp[iq] = (TH1D*)IMnpipi_overlapevt_SpSm[iq]->Clone(Form("IMnpipi_ToSm_wSp_%d",iq));
    for(int imerge=0;imerge<2;imerge++){
      for(int ibin=initbin[imerge];ibin<endbin[imerge];ibin++){
        double cont = IMnpipi_ToSm_wSp[iq]->GetBinContent(ibin);
        double ratio = IMnpipi_Sm_ratio_wSp_merge[iq]->GetBinContent(ibin);
        IMnpipi_ToSm_wSp[iq]->SetBinContent(ibin,cont*ratio);
      }
    }
    IMnpipi_Sm_merge2[iq]->Add(IMnpipi_ToSm_wK0[iq],1.0);
    IMnpipi_Sm_merge2[iq]->Add(IMnpipi_ToSm_wSp[iq],1.0);
    //subtract double-counting 
    IMnpipi_Sm_merge2[iq]->Add(IMnpipi_wK0_wSid_n_SpSm[iq],-1.0);

  
    //K0
    IMnpipi_K0_merge2[iq] = (TH1D*)IMnpipi_wK0_woSid_n_sub[iq]->Clone(Form("IMnpipi_K0_merge2_%d",iq));
    IMnpipi_ToK0_wSp[iq] = (TH1D*)IMnpipi_overlapevt_SpK0[iq]->Clone(Form("IMnpipi_ToK0_wSp_%d",iq));
    for(int imerge=0;imerge<2;imerge++){
      for(int ibin=initbin[imerge];ibin<endbin[imerge];ibin++){
        double cont = IMnpipi_ToK0_wSp[iq]->GetBinContent(ibin);
        double ratio = IMnpipi_K0_ratio_wSp_merge[iq]->GetBinContent(ibin);
        IMnpipi_ToK0_wSp[iq]->SetBinContent(ibin,cont*ratio);
      }
    }
    IMnpipi_ToK0_wSm[iq] = (TH1D*)IMnpipi_overlapevt_SmK0[iq]->Clone(Form("IMnpipi_ToK0_wSm_%d",iq));
    for(int imerge=0;imerge<2;imerge++){
      for(int ibin=initbin[imerge];ibin<endbin[imerge];ibin++){
        double cont = IMnpipi_ToK0_wSm[iq]->GetBinContent(ibin);
        double ratio = IMnpipi_K0_ratio_wSm_merge[iq]->GetBinContent(ibin);
        IMnpipi_ToK0_wSm[iq]->SetBinContent(ibin,cont*ratio);
      }
    }
    IMnpipi_K0_merge2[iq]->Add(IMnpipi_ToK0_wSp[iq],1.0);
    IMnpipi_K0_merge2[iq]->Add(IMnpipi_ToK0_wSm[iq],1.0);
    //subtract double-counting 
    IMnpipi_K0_merge2[iq]->Add(IMnpipi_wK0_wSid_n_SpSm[iq],-1.0);
  }

  
  TCanvas *cmerge2_Sp[nqcut];
  TLegend *legSp[nqcut]; 
  for(int iq=0;iq<nqcut;iq++){
    cmerge2_Sp[iq] = new TCanvas(Form("cmerge2_Sp%d",iq),Form("cmerge2_Sp%d",iq));
    cmerge2_Sp[iq]->cd();
    //IMnpipi_overlapevt_Sp[iq]->SetMarkerStyle(24);
    //IMnpipi_overlapevt_Sp[iq]->Draw("E");
    IMnpipi_Sp_merge2[iq]->SetLineColor(2);
    IMnpipi_Sp_merge2[iq]->SetMarkerColor(2);
    IMnpipi_Sp_merge2[iq]->SetMarkerStyle(20);
    IMnpipi_Sp_merge2[iq]->Draw("E");
    IMnpipi_woK0_wSid_n_woSm_sub[iq]->SetMarkerStyle(24);
    IMnpipi_woK0_wSid_n_woSm_sub[iq]->Draw("Esame");    
    IMnpipi_ToSp_wK0[iq]->SetMarkerStyle(24);
    IMnpipi_ToSp_wK0[iq]->SetMarkerColor(3);
    IMnpipi_ToSp_wK0[iq]->SetLineColor(3);
    IMnpipi_ToSp_wK0[iq]->Draw("Esame");
    IMnpipi_ToSp_wSm[iq]->SetMarkerStyle(24);
    IMnpipi_ToSp_wSm[iq]->SetMarkerColor(4);
    IMnpipi_ToSp_wSm[iq]->SetLineColor(4);
    IMnpipi_ToSp_wSm[iq]->Draw("Esame");
    legSp[iq] = new TLegend(0.6,0.7,0.9,0.9);   
    legSp[iq]->AddEntry(IMnpipi_Sp_merge2[iq],"#Sigma^{+} Total","L");
    //legSp[iq]->AddEntry(IMnpipi_overlapevt_Sp[iq],"#Sigma^{+} like events","L");
    legSp[iq]->AddEntry(IMnpipi_woK0_wSid_n_woSm_sub[iq],"#Sigma^{+} No overlap","L");
    legSp[iq]->AddEntry(IMnpipi_ToSp_wK0[iq],"Overlap with K^{0}","L");
    legSp[iq]->AddEntry(IMnpipi_ToSp_wSm[iq],"Overlap with #Sigma^{-}","L");
    legSp[iq]->Draw();
    
  }
  
  TCanvas *cmerge2_Sm[nqcut];
  TLegend *legSm[nqcut]; 
  for(int iq=0;iq<nqcut;iq++){
    cmerge2_Sm[iq] = new TCanvas(Form("cmerge2_Sm%d",iq),Form("cmerge2_Sm%d",iq));
    cmerge2_Sm[iq]->cd();
    //IMnpipi_overlapevt_Sm[iq]->SetMarkerStyle(24);
    //IMnpipi_overlapevt_Sm[iq]->Draw("E");
    IMnpipi_Sm_merge2[iq]->SetLineColor(2);
    IMnpipi_Sm_merge2[iq]->SetMarkerColor(2);
    IMnpipi_Sm_merge2[iq]->SetMarkerStyle(20);
    IMnpipi_Sm_merge2[iq]->Draw("E");
    
    IMnpipi_woK0_wSid_n_woSp_sub[iq]->SetMarkerStyle(24);
    IMnpipi_woK0_wSid_n_woSp_sub[iq]->Draw("Esame");
    IMnpipi_ToSm_wK0[iq]->SetMarkerStyle(24);
    IMnpipi_ToSm_wK0[iq]->SetMarkerColor(3);
    IMnpipi_ToSm_wK0[iq]->SetLineColor(3);
    IMnpipi_ToSm_wK0[iq]->Draw("Esame");
    IMnpipi_ToSm_wSp[iq]->SetMarkerStyle(24);
    IMnpipi_ToSm_wSp[iq]->SetMarkerColor(4);
    IMnpipi_ToSm_wSp[iq]->SetLineColor(4);
    IMnpipi_ToSm_wSp[iq]->Draw("Esame");
    legSm[iq] = new TLegend(0.6,0.7,0.9,0.9);   
    legSm[iq]->AddEntry(IMnpipi_Sm_merge2[iq],"#Sigma^{-} Total","L");
    //legSm[iq]->AddEntry(IMnpipi_overlapevt_Sm[iq],"#Sigma^{-} like events","L");
    legSm[iq]->AddEntry(IMnpipi_woK0_wSid_n_woSp_sub[iq],"#Sigma^{-} No Overlap","L");
    legSm[iq]->AddEntry(IMnpipi_ToSm_wK0[iq],"Overlap with K^{0}","L");
    legSm[iq]->AddEntry(IMnpipi_ToSm_wSp[iq],"Overlap with #Sigma^{+}","L");
    legSm[iq]->Draw();
  }
  
  TCanvas *cmerge2_K0[nqcut];
  TLegend *legK0[nqcut]; 
  for(int iq=0;iq<nqcut;iq++){
    cmerge2_K0[iq] = new TCanvas(Form("cmerge2_K0%d",iq),Form("cmerge2_K0%d",iq));
    cmerge2_K0[iq]->cd();
    //IMnpipi_overlapevt_K0[iq]->SetMarkerStyle(24);
    //IMnpipi_overlapevt_K0[iq]->Draw("E");
    IMnpipi_K0_merge2[iq]->SetLineColor(2);
    IMnpipi_K0_merge2[iq]->SetMarkerColor(2);
    IMnpipi_K0_merge2[iq]->SetMarkerStyle(20);
    IMnpipi_K0_merge2[iq]->Draw("E");
    IMnpipi_wK0_woSid_n_sub[iq]->SetMarkerStyle(24);
    IMnpipi_wK0_woSid_n_sub[iq]->Draw("Esame");
    
    IMnpipi_ToK0_wSp[iq]->SetMarkerStyle(24);
    IMnpipi_ToK0_wSp[iq]->SetMarkerColor(3);
    IMnpipi_ToK0_wSp[iq]->SetLineColor(3);
    IMnpipi_ToK0_wSp[iq]->Draw("Esame");
    IMnpipi_ToK0_wSm[iq]->SetMarkerStyle(24);
    IMnpipi_ToK0_wSm[iq]->SetMarkerColor(4);
    IMnpipi_ToK0_wSm[iq]->SetLineColor(4);
    IMnpipi_ToK0_wSm[iq]->Draw("Esame");
    legK0[iq] = new TLegend(0.6,0.7,0.9,0.9);   
    legK0[iq]->AddEntry(IMnpipi_K0_merge2[iq],"K^{0} Total","L");
    legK0[iq]->AddEntry(IMnpipi_wK0_woSid_n_sub[iq],"K^{0} No Overlap","L");
    //legK0[iq]->AddEntry(IMnpipi_overlapevt_K0[iq],"K^{0} like events","L");
    legK0[iq]->AddEntry(IMnpipi_ToK0_wSp[iq],"Overlap with #Sigma^{+}","L");
    legK0[iq]->AddEntry(IMnpipi_ToK0_wSm[iq],"Overlap with #Sigma^{-}","L");
    legK0[iq]->Draw();
    }
  
  TCanvas *c = NULL;
  TSeqCollection *SCol = gROOT->GetListOfCanvases();
  int size = SCol->GetSize();
  TIter next(SCol);
  TString pdfname = "Merge.pdf";
  for(int i=0;i<size;i++){
    //pdf->NewPage();
    c= (TCanvas*)next();
    c->Draw();
    c->cd();
    //inside the canvas
    //TPaveText *pt = new TPaveText(.74,.81,0.9,0.90,"NDC");
    c->Modified();
    c->Update();
    //std::cout << c->GetName() << std::endl;
    //make 1 pdf file
    if(i==0) c->Print(pdfname+"(",Form("pdf Title:%s",c->GetTitle()));
    else if(i==size-1)c->Print(pdfname+")",Form("pdf Title:%s",c->GetTitle())); 
    else c->Print(pdfname,Form("pdf Title:%s",c->GetTitle())); 
    
    //make separated pdf files
    //c->Print(Form("pdf/%s.pdf",c->GetTitle()));
  }

  TIter nexthist2(gDirectory->GetList());
  TFile *fout = new TFile("Merge.root","RECREATE");
  fout->Print();
  fout->cd();
  TObject *obj = NULL;
  while( (obj = (TObject*)nexthist2())!=NULL) {
    obj->Write();
  }
  fout->cd();
  fout->Close();

}
