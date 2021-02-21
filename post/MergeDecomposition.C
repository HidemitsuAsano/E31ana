void MergeDecomposition()
{
  TFile *fr[3] = {NULL};
  TFile *fmix[3] = {NULL};

  fr[0] = TFile::Open("evanaIMpisigma_npippim_v202_out_iso.root","READ");
  fmix[0] = TFile::Open("evanaIMpisigma_npippim_v202_MIX_cut4_out_iso.root","READ");
  fr[1] = TFile::Open("evanaIMpisigma_npippim_v202_out_iso_qlo.root","READ");
  fmix[1] = TFile::Open("evanaIMpisigma_npippim_v202_MIX_cut4_out_iso_qlo.root","READ");
  fr[2] = TFile::Open("evanaIMpisigma_npippim_v202_out_iso_qhi.root","READ");
  fmix[2] = TFile::Open("evanaIMpisigma_npippim_v202_MIX_cut4_out_iso_qhi.root","READ");
  fr[0]->Print();
  fmix[0]->Print();
  fr[1]->Print();
  fmix[1]->Print();
  fr[2]->Print();
  fmix[2]->Print();

  TFile *fdecompos = TFile::Open("K0SigmaTemp.root","READ");
  fdecompos->Print();

  const int nqcut=3;
  const char cqcut[][6]= {"all","qlo","qhi"};
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
 
  for(int iqcut=0;iqcut<nqcut;iqcut++){
    IMnpipi_woK0_wSid_n_woSm_sub[iqcut] = (TH1D*)q_IMnpipi_woK0_wSid_n_woSm_sub[iqcut]->ProjectionX(Form("IMnpipi_woK0_wSid_n_woSm_sub_%d",iqcut));
    IMnpipi_woK0_wSid_n_woSm_sub_merge[iqcut] = (TH1D*)IMnpipi_woK0_wSid_n_woSm_sub[iqcut]->Clone(Form("IMnpipi_woK0_wSid_n_woSm_sub_merge_%d",iqcut));
    IMnpipi_woK0_wSid_n_woSm_sub_merge[iqcut]->Add(IMnpipi_overlapdeco_Sp[0][iqcut]);
    IMnpipi_woK0_wSid_n_woSm_sub_merge[iqcut]->Add(IMnpipi_overlapdeco_Sp[2][iqcut]);
    IMnpipi_woK0_wSid_n_woSp_sub[iqcut] = (TH1D*)q_IMnpipi_woK0_wSid_n_woSp_sub[iqcut]->ProjectionX(Form("IMnpipi_woK0_wSid_n_woSp_sub_%d",iqcut));
    IMnpipi_woK0_wSid_n_woSp_sub_merge[iqcut] = (TH1D*)IMnpipi_woK0_wSid_n_woSp_sub[iqcut]->Clone(Form("IMnpipi_woK0_wSid_n_woSp_sub_merge_%d",iqcut));
    IMnpipi_woK0_wSid_n_woSp_sub_merge[iqcut]->Add(IMnpipi_overlapdeco_Sm[1][iqcut]);
    IMnpipi_woK0_wSid_n_woSp_sub_merge[iqcut]->Add(IMnpipi_overlapdeco_Sm[2][iqcut]);
    IMnpipi_wK0_woSid_n_sub[iqcut] = (TH1D*)q_IMnpipi_wK0_woSid_n_sub[iqcut]->ProjectionX(Form("IMnpipi_wK0_woSid_n_sub_%d",iqcut));
    IMnpipi_wK0_woSid_n_sub_merge[iqcut] = (TH1D*)IMnpipi_wK0_woSid_n_sub[iqcut]->Clone(Form("IMnpipi_wK0_woSid_n_sub_merge_%d",iqcut));
    IMnpipi_wK0_woSid_n_sub_merge[iqcut]->Add(IMnpipi_overlapdeco_K0[0][iqcut]);
    IMnpipi_wK0_woSid_n_sub_merge[iqcut]->Add(IMnpipi_overlapdeco_K0[1][iqcut]);
  }
  
  TCanvas *cmerge_Sp[3];
  for(int iqcut=0;iqcut<nqcut;iqcut++){
    cmerge_Sp[iqcut] = new TCanvas(Form("cmerge_Sp%d",iqcut),Form("cmerge_Sp%d",iqcut));
    cmerge_Sp[iqcut]->cd();
    IMnpipi_woK0_wSid_n_woSm_sub_merge[iqcut]->Draw("HE");
  }
  
  TCanvas *cmerge_Sm[3];
  for(int iqcut=0;iqcut<nqcut;iqcut++){
    cmerge_Sm[iqcut] = new TCanvas(Form("cmerge_Sm%d",iqcut),Form("cmerge_Sm%d",iqcut));
    cmerge_Sm[iqcut]->cd();
    IMnpipi_woK0_wSid_n_woSp_sub_merge[iqcut]->Draw("HE");
  }

  TCanvas *cmerge_K0[3];
  for(int iqcut=0;iqcut<nqcut;iqcut++){
    cmerge_K0[iqcut] = new TCanvas(Form("cmerge_K0%d",iqcut),Form("cmerge_K0%d",iqcut));
    cmerge_K0[iqcut]->cd();
    IMnpipi_wK0_woSid_n_sub_merge[iqcut]->Draw("HE");
  }
  

}
