 #include "HistManwMC.h"

void HistManwMC::initHistCom()
{
  rtFile-> cd();
  initHistCDH();

  //*** for Forward Charge ***//
  new TH1F("FC_mass2",     "FC mass^{2}",           5000, -0.5, 4.5);
  new TH2F("FC_mass2_mom", "FC mass^{2} vs mom",    5000, -0.5, 4.5, 2000, 0.0, 2.0);
  new TH2F("FC_mom_USWK_TOF", "FC mom USWK vs TOF", 2000,  0.0, 0.0, 2000, 0.0, 2.0);
  new TH1F("FC_angle",     "FC angle", 1000, -1.0, 1.0);

  //*** for CDS PID study n & CDS +- charged ***//
  new TH2F("CDS_mmDCA_wn_pm", "CDS_DCA_npipi", 500, 0.0, 10.0, 500, 0.0, 10.0);
  new TH2F("CDS_mmDCA_npipi", "CDS_DCA_npipi", 500, 0.0, 10.0, 500, 0.0, 10.0);
  new TH2F("CDS_mmDCA_tailp", "CDS_DCA_npipi", 500, 0.0, 10.0, 500, 0.0, 10.0);
  new TH2F("CDS_mmDCA_tailm", "CDS_DCA_npipi", 500, 0.0, 10.0, 500, 0.0, 10.0);

  new TH2F("CDS_DCA_wn_pm", "CDS_DCA_npipi", 500, 0.0, 10.0, 500, 0.0, 10.0);
  new TH2F("CDS_DCA_npipi", "CDS_DCA_npipi", 500, 0.0, 10.0, 500, 0.0, 10.0);
  new TH2F("CDS_DCA_tailp", "CDS_DCA_npipi", 500, 0.0, 10.0, 500, 0.0, 10.0);
  new TH2F("CDS_DCA_tailm", "CDS_DCA_npipi", 500, 0.0, 10.0, 500, 0.0, 10.0);

  new TH1F("CDH_clus_size", "CDH cluster size",         5, -0.5, 4.5);
  new TH1F("CDH_time_diff", "CDH cluster time diff", 1000, -10.0, 10.0);

  new TH2F("hitpatNC", "NC hit pattern", 16, 0.5, 16.5, 7, 0.5, 7.5);
  new TH2F("overbeta_NCseg", "1/#beta vs NC seg",   10000, 0.5, 1.5, 16, 0.5, 16.5);
  new TH2F("overbeta_NClay", "1/#beta vs NC layer", 10000, 0.5, 1.5, 7,  0.5, 7.5);

  new TH2F("KN_MM_pipi_NClay", "d(K^{-}, n)\"X\" vs NC layer", 400, 0.0, 2.0,  7, 0.5, 7.5);
  new TH2F("KN_MM_pipi_NCseg", "d(K^{-}, n)\"X\" vs NC seg",   400, 0.0, 2.0, 16, 0.5, 16.5);

  new TH2F("KN_MM_pipi_wN_NClay", "d(K^{-}, n)\"X\" vs NC layer", 400, 0.0, 2.0,  7, 0.5, 7.5);
  new TH2F("KN_MM_pipi_wN_NCseg", "d(K^{-}, n)\"X\" vs NC seg",   400, 0.0, 2.0, 16, 0.5, 16.5);

  new TH2F("KN_MM_NClay", "d(K^{-}, n)\"X\" vs NC layer", 400, 0.0, 2.0,  7, 0.5, 7.5);
  new TH2F("KN_MM_NCseg", "d(K^{-}, n)\"X\" vs NC seg",   400, 0.0, 2.0, 16, 0.5, 16.5);

  new TH2F("nCDC",       "CDC multiplicity",      10, -0.5, 9.5, NumOfCDCLayers, 0.5, NumOfCDCLayers+0.5);
  new TH2F("nCDC_tail",  "CDC multiplicity",      10, -0.5, 9.5, NumOfCDCLayers, 0.5, NumOfCDCLayers+0.5);
  new TH2F("nCDC_wn_pm", "CDC multiplicity",      10, -0.5, 9.5, NumOfCDCLayers, 0.5, NumOfCDCLayers+0.5);
  new TH2F("nCDC_wn_pm_tail", "CDC multiplicity", 10, -0.5, 9.5, NumOfCDCLayers, 0.5, NumOfCDCLayers+0.5);

  for( int i=0; i<20; i++ ){
    new TH1F(Form("trig_%d", i), Form("trigger pattern %d", i), 4096, -0.5, 4095.5);
  }
  new TH1F("BeamMom_diff_MC", "Beam Mom  MC val - Ana val", 200, -0.1, 0.1);
  new TH2F("Vtx_XY_diff_MC", "Vertex XY  MC - Ana", 2000, -2.0, 2.0, 2000, -2.0, 2.0);
  new TH1F("Vtx_Z_diff_MC", "Vertex Z  MC - Ana", 2000, -10, 10);

  new TH1F("KN_MM_diff_MC", "d(K^{-}, n)\"X\" MC val - Ana -val", 200, -0.1, 0.1);
  new TH1F("KN_MM_pipi_diff_MC", "d(K^{-}, n)\"X\" MC val - Ana -val", 200, -0.1, 0.1);
  new TH1F("KN_MM_pipi_wK0_diff_MC", "d(K^{-}, n)\"X\" MC val - Ana -val", 200, -0.1, 0.1);
  new TH1F("KN_MM_MC", "d(K{-}, n)\"X\" by MC", 3000, 0.0, 3.0);

  new TH1F("T0NC_FL_diff_MC", "T0-NC FL diff", 1000, -25, 25);

  for( int lay=1; lay<=NumOfCDCLayers; lay++ ){
    new TH1F(Form("CDC_res_%d", lay), Form("CDC residual layer%d", lay), 1000, -1, 1);
    new TH1F(Form("CDC_res_%d_pip", lay), Form("CDC residual layer%d", lay), 1000, -1, 1);
    new TH1F(Form("CDC_res_%d_pim", lay), Form("CDC residual layer%d", lay), 1000, -1, 1);
    new TH1F(Form("CDC_res_%d_km", lay),  Form("CDC residual layer%d", lay), 1000, -1, 1);
    new TH1F(Form("CDC_res_%d_p", lay),   Form("CDC residual layer%d", lay), 1000, -1, 1);
  }
  new TH1F("CDC_nGoodTrack_wo", "CDC n track", 10, -0.5, 9.5);
  new TH1F("CDC_nGoodTrack",    "CDC n track", 10, -0.5, 9.5);

  new TH1F("CDC_trig", "CDC trigger",   36, 0.5, 36.5);
  new TH1F("CDC_eff",  "CDC effective", 36, 0.5, 36.5);
  new TH1F("CDC_eff2", "CDC effective", 36, 0.5, 36.5);

  new TH1F("CDC_IH_trig", "CDC trigger",   24, 0.5, 24.5);
  new TH1F("CDC_IH_eff",  "CDC effective", 24, 0.5, 24.5);
  new TH1F("CDC_IH_eff2", "CDC effective", 24, 0.5, 24.5);

  new TH1F("T0CDH_TOF", "T0-CDH TOF", 1100, -10, 100);
  new TH2F("CDS_beta_overmom",      "CDS #beta vs 1/mom", 1200, 0.0, 1.2, 2000, -10.0, 10.0);
  new TH2F("CDS_beta_overmom_pi",   "CDS #beta vs 1/mom", 1200, 0.0, 1.2, 2000, -10.0, 10.0);

  new TH2F("CDS_overbeta_mom",      "CDS 1/#beta vs mom", 2000, 0.0, 9.5, 2000, -1.0, 1.0);
  new TH2F("CDS_overbeta_mom_pi",   "CDS 1/#beta vs mom", 2000, 0.0, 9.5, 2000, -1.0, 1.0);
  new TH2F("CDS_overbeta_mom_e",    "CDS 1/#beta vs mom", 2000, 0.0, 9.5, 2000, -1.0, 1.0);
  
  new TH2F("CDS_mass2_mom",     "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
  new TH2F("CDS_mass2_mom_pi",  "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
  new TH2F("CDS_mass2_mom_e",   "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
  new TH2F("CDS_mass2_mom_pim", "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
  new TH2F("CDS_mass2_mom_pip", "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
  new TH2F("CDS_mass2_mom_K",   "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
  new TH2F("CDS_mass2_mom_p",   "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
  new TH2F("CDS_mass2_mom_d",   "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);

  new TH2F("CDS_mass2_mom_CDHshare", "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
  new TH2F("CDS_mass2_mom_CDH2",     "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
  for( int seg=1; seg<=36; seg++ ){
    new TH2F(Form("CDS_mass2_mom_CDH_seg%d", seg), Form("CDS mass2 vs mom seg%d", seg), 500, -0.5, 2.3, 250, -1.0, 1.0);
  }

  new TH2F("CDS_mass2_mom_1p2m", "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);

  if( simReader ){
    new TH2F("CDS_mass2_mom_MC_e",     "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_MC_ep",    "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_MC_pim",   "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_MC_pip",   "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_MC_km",    "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_MC_kp",    "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_MC_p",     "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_MC_mum",   "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_MC_mup",   "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_MC_other", "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);

    new TH2F("CDS_mass2_mom_MC_e_2hit",     "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_MC_ep_2hit",    "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_MC_pim_2hit",   "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_MC_pip_2hit",   "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_MC_km_2hit",    "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_MC_kp_2hit",    "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_MC_p_2hit",     "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_MC_mum_2hit",   "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_MC_mup_2hit",   "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);

    new TH2F("CDS_mass2_mom_MC_e_2tra",     "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_MC_ep_2tra",    "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_MC_pim_2tra",   "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_MC_pip_2tra",   "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_MC_km_2tra",    "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_MC_kp_2tra",    "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_MC_p_2tra",     "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_MC_mum_2tra",   "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_MC_mup_2tra",   "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);

    new TH2F("CDS_mass2_mom_MC_e_2",     "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_MC_ep_3",    "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_MC_mum_10",  "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_MC_mum_11",  "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_MC_mup_15",  "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_MC_mup_16",  "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_MC_pim_0",   "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_MC_pim_20",  "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_MC_pim_21",  "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_MC_pim_22",  "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_MC_pip_0",   "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_MC_pip_30",  "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_MC_pip_31",  "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_MC_pip_32",  "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_MC_km_0",    "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_MC_p_0",     "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_MC_p_50",    "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_MC_p_51",    "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_MC_p_52",    "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_MC_p_53",    "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);

    new TH2F("CDS_mass2_mom_wn_pm_MC_e",     "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_wn_pm_MC_ep",    "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_wn_pm_MC_pim",   "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_wn_pm_MC_pip",   "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_wn_pm_MC_km",    "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_wn_pm_MC_kp",    "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_wn_pm_MC_p",     "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_wn_pm_MC_mum",   "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_wn_pm_MC_mup",   "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_wn_pm_MC_other", "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);

    new TH2F("CDS_mass2_mom_wn_pm_MC_e_2hit",     "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_wn_pm_MC_ep_2hit",    "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_wn_pm_MC_pim_2hit",   "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_wn_pm_MC_pip_2hit",   "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_wn_pm_MC_km_2hit",    "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_wn_pm_MC_kp_2hit",    "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_wn_pm_MC_p_2hit",     "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_wn_pm_MC_mum_2hit",   "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_wn_pm_MC_mup_2hit",   "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);

    new TH2F("CDS_mass2_mom_wn_pm_MC_e_2tra",     "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_wn_pm_MC_ep_2tra",    "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_wn_pm_MC_pim_2tra",   "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_wn_pm_MC_pip_2tra",   "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_wn_pm_MC_km_2tra",    "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_wn_pm_MC_kp_2tra",    "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_wn_pm_MC_p_2tra",     "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_wn_pm_MC_mum_2tra",   "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_wn_pm_MC_mup_2tra",   "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);

    new TH2F("CDS_mass2_mom_wn_pm_MC_e_2",     "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_wn_pm_MC_ep_3",    "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_wn_pm_MC_mum_10",  "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_wn_pm_MC_mum_11",  "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_wn_pm_MC_mup_15",  "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_wn_pm_MC_mup_16",  "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_wn_pm_MC_pim_0",   "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_wn_pm_MC_pim_20",   "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_wn_pm_MC_pim_21",   "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_wn_pm_MC_pim_22",   "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_wn_pm_MC_pip_0",    "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_wn_pm_MC_pip_30",   "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_wn_pm_MC_pip_31",   "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_wn_pm_MC_pip_32",   "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_wn_pm_MC_km_0",    "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_wn_pm_MC_p_0",      "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_wn_pm_MC_p_50",     "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_wn_pm_MC_p_51",     "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_wn_pm_MC_p_52",     "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
    new TH2F("CDS_mass2_mom_wn_pm_MC_p_53",     "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
  }

  new TH1F("CDS_mass2_plus",     "CDS mass2 mom>0", 2000, -0.5, 4.5);
  new TH1F("CDS_mass2_minus",    "CDS mass2 mom<0", 2000, -0.5, 4.5);
  new TH1F("CDS_mass2_plus_pi",  "CDS mass2 mom>0", 2000, -0.5, 4.5);
  new TH1F("CDS_mass2_minus_pi", "CDS mass2 mom<0", 2000, -0.5, 4.5);
  new TH1F("CDS_mass2_pim", "CDS mass2 #pi^{-}", 2000, -0.5, 4.5);
  new TH1F("CDS_mass2_K",   "CDS mass2 K^{-}",   2000, -0.5, 4.5);
  new TH1F("CDS_mass2_pip", "CDS mass2 #pi^{+}", 2000, -0.5, 4.5);
  new TH1F("CDS_mass2_p",   "CDS mass2 p", 2000, -0.5, 4.5);
  new TH1F("CDS_mass2_d",   "CDS mass2 d", 2000, -0.5, 4.5);
  new TH2F("CDS_mass2_mom_chi30",     "CDS mass2 vs mom", 250, -0.5, 4.5, 250, -1.0, 1.0);

  new TH1F("CDS_chi2",     "CDS chi2", 1000, 0.0, 100.);

  new TH2F("CDS_mass2_mom_wn_pm",     "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
  new TH2F("CDS_mass2_mom_wn_pm_3",   "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
  new TH2F("CDS_mass2_mom_wn_pm_pip", "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 0.0);
  new TH2F("CDS_mass2_mom_wn_pm_km",  "CDS mass2 vs mom", 500, -0.5, 2.5, 250, 0.0, 1.0);
  new TH2F("CDS_mass2_mom_wn_pm_p",   "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 0.0);
  new TH2F("CDS_mass2_mom_wn_pm_CDHshare", "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);
  new TH2F("CDS_mass2_mom_wn_pm_CDH2",     "CDS mass2 vs mom", 500, -0.5, 2.5, 250, -1.0, 1.0);

  new TH1F("KN_MM_km_2track", "d(K^{-}, n) w/ K^{-}", 2000, 0.0, 2.0);

  new TH2F("CDS_mom_diff_pip_MC", "CDS mom vs mom diff #pi^{+}", 1000, 0.0, 1.0, 1000, -0.1, 0.1);
  new TH2F("CDS_mom_diff_p_MC", "  CDS mom vs mom diff p",       1000, 0.0, 1.0, 1000, -0.1, 0.1);
  new TH2F("CDS_mom_diff_pim_MC", "CDS mom vs mom diff #pi^{-}", 1000, 0.0, 1.0, 1000, -0.1, 0.1);
  new TH2F("CDS_mom_diff_km_MC",  "CDS mom vs mom diff K^{-}",   1000, 0.0, 1.0, 1000, -0.1, 0.1);

  new TH1F("Vtx_Z_r2", "Vtx Z R<2", 1000, -50, 50);
  new TH2F("Vtx_XY_z0", "Vtx XY z=", 500, -25, 25, 500, -25, 25);

  new TH2F("Vtx_XY", "Vertex XY Plane", 1000, -30, 30, 1000,-30, 30);
  new TH2F("Vtx_ZX", "Vertex ZX Plane", 1000, -50, 50, 1000,-30, 30);
  new TH2F("Vtx_ZY", "Vertex ZY Plane", 1000, -50, 50, 1000,-30, 30);

  new TH2F("Vtx_XY_wtar", "Vertex XY Plane", 1000, -30, 30, 1000,-30, 30);
  new TH2F("Vtx_ZX_wtar", "Vertex ZX Plane", 1000, -50, 50, 1000,-30, 30);
  new TH2F("Vtx_ZY_wtar", "Vertex ZY Plane", 1000, -50, 50, 1000,-30, 30);
  new TH2F("Vtx_XY_wtar1", "Vertex XY Plane", 1000, -30, 30, 1000,-30, 30);
  new TH2F("Vtx_ZX_wtar1", "Vertex ZX Plane", 1000, -50, 50, 1000,-30, 30);
  new TH2F("Vtx_ZY_wtar1", "Vertex ZY Plane", 1000, -50, 50, 1000,-30, 30);
  new TH2F("Vtx_XY_wtar2", "Vertex XY Plane", 1000, -30, 30, 1000,-30, 30);
  new TH2F("Vtx_ZX_wtar2", "Vertex ZX Plane", 1000, -50, 50, 1000,-30, 30);
  new TH2F("Vtx_ZY_wtar2", "Vertex ZY Plane", 1000, -50, 50, 1000,-30, 30);
  new TH2F("Vtx_XY_wtar3", "Vertex XY Plane", 1000, -30, 30, 1000,-30, 30);
  new TH2F("Vtx_ZX_wtar3", "Vertex ZX Plane", 1000, -50, 50, 1000,-30, 30);
  new TH2F("Vtx_ZY_wtar3", "Vertex ZY Plane", 1000, -50, 50, 1000,-30, 30);
  new TH2F("Vtx_XY_wtar4", "Vertex XY Plane", 1000, -30, 30, 1000,-30, 30);
  new TH2F("Vtx_ZX_wtar4", "Vertex ZX Plane", 1000, -50, 50, 1000,-30, 30);
  new TH2F("Vtx_ZY_wtar4", "Vertex ZY Plane", 1000, -50, 50, 1000,-30, 30);
  new TH2F("Vtx_XY_wtar5", "Vertex XY Plane", 1000, -30, 30, 1000,-30, 30);
  new TH2F("Vtx_ZX_wtar5", "Vertex ZX Plane", 1000, -50, 50, 1000,-30, 30);
  new TH2F("Vtx_ZY_wtar5", "Vertex ZY Plane", 1000, -50, 50, 1000,-30, 30);

  new TH1F("Vtx_Z_min",       "#pi^{+} #pi^{-} vtx Z",  1000, -50, 50);
  new TH2F("Vtx_XY_min_z0",   "#pi^{+} #pi^{-} vtx XY", 1000, -20, 20, 1000, -20, 20);

  new TH2F("Vtx_XY_pipi",   "#pi^{+} #pi^{-} vtx XY", 1000, -20, 20, 1000, -20, 20);
  new TH1F("Vtx_Z_pipi",       "#pi^{+} #pi^{-} vtx Z",  1000, -50, 50);

  new TH1F("DCA_pipi_beam",     "DCA #pi^{+} #pi^{-} & beam", 1000, 0.0, 50.0);
  new TH1F("DCA_pipi_beam_woF", "DCA #pi^{+} #pi^{-} & beam", 1000, 0.0, 50.0);

  new TH1F("nNC", "NC multiplicity", 113, -0.5, 112.5);
  new TH1F("nlayNC", "NC hit layer", 7, -0.5, 7.5);

  new TH1F("nNCCluster", "Num of NC Cluster", 10, -0.5, 9.5);
  new TH1F("NCClusterSize", "NC Cluster Size", 20, 0.5, 20.5);
  new TH2F("nNCCluster_Size", "NC n Cluster vs Size", 10, -0.5, 9.5, 20, 0.5, 20.5);

  new TH1F("NC_hitpos", "NC hitpos", 2000, -100, 100);
  new TH1F("NC_overbeta", "NC 1/#beta dE>8MeVee", 1000, 0.0, 10.0);
  new TH1F("NC_overbeta_wtar", "NC 1/#beta dE>8MeVee", 1000, 0.0, 10.0);
  new TH1F("NC_overbeta_8MeVee", "NC 1/#beta dE>8MeVee", 10000, 0.0, 10.0);
  new TH1F("NC_overbeta_8MeVee_wtar", "NC 1/#beta dE>8MeVee", 1000, 0.0, 10.0);

  new TH1F("NC_overbeta_8MeVee_gamma", "NC 1/#beta dE>8MeVee",      1000, 0.9, 1.1);
  new TH1F("NC_overbeta_8MeVee_gamma_wtar", "NC 1/#beta dE>8MeVee", 1000, 0.9, 1.1);
  for( int seg=1; seg<=112; seg++ ){
    new TH1F(Form("NC_overbeta_8MeVee_gamma_seg%d", seg), Form("NC 1/#beta dE>8MeVee NC seg%d", seg),      1000, 0.9, 1.1);
    new TH1F(Form("NC_overbeta_8MeVee_gamma_wtar_seg%d", seg), Form("NC 1/#beta dE>8MeVee NC seg%d", seg), 1000, 0.9, 1.1);
  }

  //*******************//
  //*** NC S/N test ***//
  //*******************//
  for( int i=0; i<30; i++ ){
    new TH1F(Form("NC_overbeta_%d", 2*i+1), "NC 1/#beta dE>8MeVee", 1000, 0.0, 10.0);
    new TH1F(Form("NC_overbeta_%d_wtar", 2*i+1),"NC 1/#beta dE>8MeVee", 1000, 0.0, 10.0);
  }

  new TH2F("NC_overbeta_dE", "NC 1/#beta",            2000, 0.0, 5.0, 50, 0.0, 200);
  new TH2F("NC_overbeta_dE_wtar", "NC 1/#beta vs dE", 2000, 0.0, 5.0, 50, 0.0, 200);

  //*************************//
  //*** for BPC pos check ***//
  //*************************//
  new TH1F("DCA_pkm", "p K^{-} DCA", 1000, 0, 10);

  new TH2F("CDS2_BPC_dx_dy", "CDS2 & BPC dx vs dy", 1000, -10, 10, 1000, -10, 10);
  new TH2F("CDS2_BPC_z_dx",  "CDS2 & BPC z vs dx",  1000, -30, 30, 1000, -10, 10);
  new TH2F("CDS2_BPC_z_dy",  "CDS2 & BPC z vs dy",  1000, -30, 30, 1000, -10, 10);
  new TH2F("CDS2_BPC_x_dx",  "CDS2 & BPC z vs dx",  1000, -30, 30, 1000, -10, 10);
  new TH2F("CDS2_BPC_y_dy",  "CDS2 & BPC z vs dy",  1000, -30, 30, 1000, -10, 10);

  new TH1F("CDS_IM_pipi", "CDS #pi^{+} #pi^{-} IM",       10000, 0.0, 1.0);
  new TH2F("pim_mom_CDS_IM_pipi", "#pi^{-} mom vs #pi^{+} #pi^{-} IM", 1000, 0.0, 1.0, 1000, 0.0, 1.0);
  new TH2F("pip_mom_CDS_IM_pipi", "#pi^{+} mom vs #pi^{+} #pi^{-} IM", 1000, 0.0, 1.0, 1000, 0.0, 1.0);

  new TH1F("CDS_IM_ppim", "CDS p #pi^{-} IM",             10000, 1.0, 2.0);
  new TH2F("p_mom_CDS_IM_ppim",   "p mom vs p #pi^{-} IM",       1000, 0.0, 1.0, 1000, 1.0, 2.0);
  new TH2F("pim_mom_CDS_IM_ppim", "#pi^{-} mom vs p #pi^{-} IM", 1000, 0.0, 1.0, 1000, 1.0, 2.0);

  new TH1F("CDS_IM_pkm",     "CDS p K^{-} IM",               2000, 1.0, 3.0);
  new TH1F("KCDSppim_MM",    "d(K^{-}, p #pi^{-}) MM",       2000, 0.0, 2.0);
  new TH1F("KCDSppim_MM2",   "d(K^{-}, p #pi^{-}) MM",       2000, 0.0, 2.0);
  new TH1F("KCDSppim_MM_wL", "d(K^{-}, p #pi^{-}) MM",       2000, 0.0, 2.0);
  new TH1F("KCDSppim_MM2_wL", "d(K^{-}, p #pi^{-}) MM",      2000, -1.0, 1.0);
  new TH1F("KCDSpkm_MM",     "d(K^{-}, p K^{-}) MM",         2000, 0.0, 2.0);

  new TH1F("KN_MM_wo",      "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pim_wo",  "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pip_wo",  "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_km_wo",   "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_p_wo",    "d(K^{-}, n)", 3000, 0.0, 3.0);

  new TH1F("KN_MM",      "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pim",  "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pip",  "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_km",   "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_p",    "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_wSp",  "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_wSm",  "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_wS",   "d(K^{-}, n)", 3000, 0.0, 3.0);

  new TH1F("Npip_IM", "n #pi^{+} IM", 1000, 1.0, 2.0);
  new TH1F("Npim_IM", "n #pi^{-} IM", 1000, 1.0, 2.0);
  new TH1F("Nkm_IM",  "n K^{-} IM",   1000, 1.0, 2.0);
  new TH1F("KNkm_MM", "d(K^{-}, n K^{-})\"X\"", 2000, 0.0, 2.0);

  initHistCDC();
  initHistT0NC();
  initHistIH();
}
