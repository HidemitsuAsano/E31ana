void GetAccMap()
{
  
  TFile *fSp[4]={NULL};
  TFile *fSm[4]={NULL};
  TFile *fK0[4]={NULL};
  TFile *fK0gen=NULL;

  fSp[0] = TFile::Open("simIMpisigma_nSppim_pippimn_v132_out_iso_rej.root");
  fSp[1] = TFile::Open("simIMpisigma_nSppim_pippimn_v132_out_iso_qlo_rej.root");
  fSp[2] = TFile::Open("simIMpisigma_nSppim_pippimn_v132_out_iso_qhi_rej.root");
  fSp[3] = TFile::Open("simIMpisigma_nSppim_pippimn_v132_out_iso_theta15_rej.root");

  fSm[0] = TFile::Open("simIMpisigma_nSmpip_pippimn_v132_out_iso_rej.root");
  fSm[1] = TFile::Open("simIMpisigma_nSmpip_pippimn_v132_out_iso_qlo_rej.root");
  fSm[2] = TFile::Open("simIMpisigma_nSmpip_pippimn_v132_out_iso_qhi_rej.root");
  fSm[3] = TFile::Open("simIMpisigma_nSmpip_pippimn_v132_out_iso_theta15_rej.root");

  fK0[0] = TFile::Open("simIMpisigma_K0nn_pippimn_v11_out_iso_rej.root");
  fK0[1] = TFile::Open("simIMpisigma_K0nn_pippimn_v11_out_iso_qlo_rej.root");
  fK0[2] = TFile::Open("simIMpisigma_K0nn_pippimn_v11_out_iso_qhi_rej.root");
  fK0[3] = TFile::Open("simIMpisigma_K0nn_pippimn_v11_out_iso_theta15_rej.root");

  fK0gen = TFile::Open("simIMpisigma_K0nn_v11.root");
  
  const int nqcut=4;
  TH2D* q_IMnpipi_gen_Sp[nqcut];
  TH2D* q_IMnpipi_gen_Sm[nqcut];
  TH2D* q_IMnpipi_gen_K0[nqcut];
  for(int iq=0;iq<nqcut;iq++){
    q_IMnpipi_gen_Sp[iq] = (TH2D*) fSp[iq]->Get("q_IMpiSigma_gen");
    q_IMnpipi_gen_Sm[iq] = (TH2D*) fSm[iq]->Get("q_IMpiSigma_gen");
    q_IMnpipi_gen_K0[iq] = (TH2D*) fK0gen->Get("React_q_IMnpipi");
    q_IMnpipi_gen_Sp[iq]->RebinY(
  }
  
  
  
  TH2D* q_IMnpipi_wSid_n_Sp_reco[nqcut];
  TH2D* q_IMnpipi_wSid_n_Sm_reco[nqcut];
  TH2D* q_IMnpipi_wK0_n_K0_reco[nqcut];
  for(int iq=0;iq<nqcut;iq++){
    q_IMnpipi_wSid_n_Sp_reco[iq] = (TH2D*)fSp[iq]->Get("q_IMnpipi_wSid_n_Sp");
    q_IMnpipi_wSid_n_Sm_reco[iq] = (TH2D*)fSm[iq]->Get("q_IMnpipi_wSid_n_Sm");
    q_IMnpipi_wK0_n_K0_reco[iq] = (TH2D*)fK0[iq]->Get("q_IMnpipi_wK0_n");
  }
  
  TH2D* q_IMnpipi_Sp_acc[nqcut];
  TH2D* q_IMnpipi_Sm_acc[nqcut];
  TH2D* q_IMnpipi_K0_acc[nqcut];
  TH2D* q_IMnpipi_Sp_accerr[nqcut];
  TH2D* q_IMnpipi_Sm_accerr[nqcut];
  TH2D* q_IMnpipi_K0_accerr[nqcut];
  
  for(int iq=0;iq<nqcut;iq++){
    q_IMnpipi_Sp_acc[iq] = (TH2D*)q_IMnpipi_wSid_n_Sp_reco[iq]->Clone(Form("q_IMnpipi_Sp_acc_%d",iq));
    q_IMnpipi_Sm_acc[iq] = (TH2D*)q_IMnpipi_wSid_n_Sm_reco[iq]->Clone(Form("q_IMnpipi_Sm_acc_%d",iq));
    q_IMnpipi_K0_acc[iq] = (TH2D*)q_IMnpipi_wK0_n_K0_reco[iq]->Clone(Form("q_IMnpipi_K0_acc_%d",iq));
    q_IMnpipi_Sp_acc[iq]->Divide(q_IMnpipi_Sp_acc[iq],q_IMnpipi_gen_Sp[iq],1.0,1.0,"b");
    q_IMnpipi_Sm_acc[iq]->Divide(q_IMnpipi_Sp_acc[iq],q_IMnpipi_gen_Sm[iq],1.0,1.0,"b");
    q_IMnpipi_K0_acc[iq]->Divide(q_IMnpipi_K0_acc[iq],q_IMnpipi_gen_K0[iq],1.0,1.0,"b");
  }


  TCanvas *cSp[nqcut];
  TCanvas *cSm[nqcut];
  TCanvas *cK0[nqcut];

 



}
