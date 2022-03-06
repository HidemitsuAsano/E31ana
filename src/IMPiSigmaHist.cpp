#include "IMPiSigmaHist.h"
#include "Tools.h"
#include <TMath.h>

void InitBasicHist(const bool MCFlag)
{
  // geneneral informantion **//
  Tools::newTH1F( Form("Time"), 3000, -0.5, 2999.5 );
  Tools::newTH1F( Form("nTrack"), 11, -0.5, 10.5 );
  Tools::newTH1F( Form("nTrack_If2GoodTracks"), 11, -0.5, 10.5 );
  Tools::newTH1F( Form("nGoodTrack"), 11, -0.5, 10.5 );
  Tools::newTH1F( Form("Scaler"), 41, -0.5, 40.5 );
  Tools::newTH1F( Form("Trigger"), 10,-0.5,9.5);
  Tools::newTH1F( Form("Trigmode"), 20,-0.5,19.5);
  Tools::newTH1F( Form("Trigmode_Kf"), 20,-0.5,19.5);

  // CDC and CDH information from CDC-trackig file **//
  Tools::newTH1F( Form("mul_CDH"),Form("CDH multiplicity"), 11, -0.5, 10.5 );
  Tools::newTH1F( Form("mul_CDH_iso"),Form("CDH multiplicity after isolation cut"), 11, -0.5, 10.5 );
  Tools::newTH1F( Form("mul_CDH_assoc"),Form("n CDH associated with CDCtrack" ),11, -0.5, 10.5 );
  Tools::newTH1F( Form("npimangle"),628, 0, 2*3.14);
  Tools::newTH1F( Form("npipangle"),628, 0, 2*3.14);
  Tools::newTH1F( Form("CDCInner3Mul"),20,0,20);
  Tools::newTH1F( Form("NCDCOutHit"),30,0,30);

  if(MCFlag){
    //vertex reso//
    Tools::newTH1F( Form("vertex_diff_X"), 100, -50, 50 );
    Tools::newTH1F( Form("vertex_diff_Y"), 100, -50, 50 );
    Tools::newTH1F( Form("vertex_diff_Z"), 100, -50, 50 );
  }
  //===== CDH hit pos study =====//
  Tools::newTH2F( Form("CDH_mom_diffpos_pi_phi"), 100, -10, 10, 100, -1.0, 1.0 );
  Tools::newTH2F( Form("CDH_mom_diffpos_pi_z"), 1000, -50, 50, 100, -1.0, 1.0 );
  
  Tools::newTH2F( Form("CDH_mom_TOF_pi"),   100, -2, 2, 200, -1.0, 1.0 );
  Tools::newTH2F( Form("CDH_mom_diffpos_p_phi"), 100, -10, 10, 100, -1.0, 1.0 );
  Tools::newTH2F( Form("CDH_mom_diffpos_p_z"), 1000, -50, 50, 100, -1.0, 1.0 );
  
  Tools::newTH2F( Form("CDH_mom_TOF_p"),   100, -2, 2, 200, -1.0, 1.0 );
  for(int iseg=0;iseg<36;iseg++){
    Tools::newTH2F( Form("CDH%d_mom_TOF_pi",iseg+1),       100, -2, 2, 200, -1.0, 1.0 );
    Tools::newTH2F( Form("CDH%d_mom_TOF_p",iseg+1),       100, -2, 2, 200, -1.0, 1.0 );
  }
  Tools::newTH2F( Form("CDH_diffpos_z_p_z"), 1000,-50,50,1000,-50,50);
  
  //** beam line **//
  Tools::newTH1F( Form("mul_BHD"), 12, -0.5, 11.5 );
  Tools::newTH1F( Form("mul_T0"),   6, -0.5, 5.5 );
  Tools::newTH1F( Form("T0time"),   1000, -50, 50 );
  Tools::newTH1F( Form("tof_T0BHD"), 2000, 20, 40 );
  Tools::newTH1F( Form("tracktime_BPC"),  1200, -200, 400 );
  Tools::newTH1F( Form("trackchi2_BPC"),  500, 0, 50 );
  Tools::newTH1F( Form("ntrack_BPC"),  6, -0.5, 5.5 );
  Tools::newTH1F( Form("tracktime_BLC1"), 1200, -200, 400 );
  Tools::newTH1F( Form("tracktime_BLC2"), 1200, -200, 400 );
  Tools::newTH1F( Form("trackchi2_BLC1"), 500, 0, 50 );
  Tools::newTH1F( Form("trackchi2_BLC2"), 500, 0, 50 );
  Tools::newTH1F( Form("ntrack_BLC1"), 6, -0.5, 5.5 );
  Tools::newTH1F( Form("ntrack_BLC2"), 6, -0.5, 5.5 );
  Tools::newTH2F( Form("dydx_BLC2BPC"),     500, -5, 5, 500, -5, 5 );
  Tools::newTH2F( Form("dydzdxdz_BLC2BPC"), 175, -0.035, 0.035, 175, -0.035, 0.035 );
  Tools::newTH1F( Form("trackchi2_beam"), 500, 0, 50 );
  Tools::newTH2F( Form("D5chi2_mom"),180,0.92,1.10, 500, 0, 50 );
  Tools::newTH1F( Form("momentum_beam"), 180, 0.92, 1.10 );
  Tools::newTH1F( Form("PID_beam"),5,-1.5,3.5);
  //** CDS **//
  Tools::newTH1F( Form("trackchi2_CDC"), 1000, 0, 50 );
  Tools::newTH2F( Form("PID_CDS_beta"), 2000, 0, 10., 1000, -1.2, 1.2 );
  Tools::newTH2F( Form("PID_CDS_beta_select"), 2000, 0, 10., 1000, -1.2, 1.2 );
  Tools::newTH2F( Form("PID_CDS_beta_select2"), 2000, 0, 10., 1000, -1.2, 1.2 );
  Tools::newTH2F( Form("PID_CDS"), 1000, -0.6, 5, 1000, -1.2, 1.2 );
  Tools::newTH2F( Form("PID_CDS_select"), 1000, -0.6, 5, 1000, -1.2, 1.2 );
  Tools::newTH2F( Form("PID_CDS_select2"), 1000, -0.6, 5, 1000, -1.2, 1.2 );
  Tools::newTH2F( Form("PID_CDS_PIM_beta"), 400, 0, 10, 200, -1.2, 1.2 );
  Tools::newTH2F( Form("PID_CDS_PIM_beta_select"), 400, 0, 10, 200, -1.2, 1.2 );
  Tools::newTH2F( Form("PID_CDS_PIM_beta_select2"), 400, 0, 10, 200, -1.2, 1.2 );
  Tools::newTH2F( Form("PID_CDS_PIM"), 200, -0.6, 5, 200, -1.2, 1.2 );
  Tools::newTH2F( Form("PID_CDS_PIM_select"), 200, -0.6, 5, 200, -1.2, 1.2 );
  Tools::newTH2F( Form("PID_CDS_PIM_select2"), 200, -0.6, 5, 200, -1.2, 1.2 );
  Tools::newTH2F( Form("PID_CDS_PIP_beta"), 400, 0, 10, 200, -1.2, 1.2 );
  Tools::newTH2F( Form("PID_CDS_PIP_beta_select"), 400, 0, 10, 200, -1.2, 1.2 );
  Tools::newTH2F( Form("PID_CDS_PIP_beta_select2"), 400, 0, 10, 200, -1.2, 1.2 );
  Tools::newTH2F( Form("PID_CDS_PIP"), 200, -0.6, 5, 200, -1.2, 1.2 );
  Tools::newTH2F( Form("PID_CDS_PIP_select"), 200, -0.6, 5, 200, -1.2, 1.2 );
  Tools::newTH2F( Form("PID_CDS_PIP_select2"), 200, -0.6, 5, 200, -1.2, 1.2 );
  Tools::newTH2F( Form("PID_CDS_Proton"), 200, -0.6, 5, 200, -1.2, 1.2 );
  Tools::newTH2F( Form("PID_CDS_Proton_select"), 200, -0.6, 5, 200, -1.2, 1.2 );
  Tools::newTH2F( Form("PID_CDS_Proton_select2"), 200, -0.6, 5, 200, -1.2, 1.2 );
  Tools::newTH2F( Form("PID_CDS_Kaon"), 200, -0.6, 5, 200, -1.2, 1.2 );
  Tools::newTH2F( Form("PID_CDS_Kaon_select"), 200, -0.6, 5, 200, -1.2, 1.2 );
  Tools::newTH2F( Form("PID_CDS_Kaon_select2"), 200, -0.6, 5, 200, -1.2, 1.2 );
  Tools::newTH1F( Form("ntrack_CDS"), 6, -0.5, 5.5 );
  Tools::newTH1F( Form("ntrack_pi_plus"), 6, -0.5, 5.5 );
  Tools::newTH1F( Form("ntrack_proton"), 6, -0.5, 5.5 );
  //Tools::newTH1F( Form("ntrack_deuteron"), 6, -0.5, 5.5 );
  Tools::newTH1F( Form("ntrack_pi_minus"), 6, -0.5, 5.5 );
  Tools::newTH1F( Form("ntrack_K_minus"), 6, -0.5, 5.5 );
  Tools::newTH1F( Form("CDHNeutralSeg"),36, 0.5, 36.5);
  Tools::newTH2F( Form("CDHseg_MMass_fid_beta_dE_woK0"),36,0.5,36.5,140,0.4,1.8);
  Tools::newTH2F( Form("CDHz_MMass_fid_beta_dE_woK0"),100,-50,50,140,0.4,1.8);
  //Tools::newTH2F( Form("zVTX_MMass_fid_beta_dE_woK0"),100,-12,2,140,0.4,1.8);
  Tools::newTH2F( Form("CDHz_IMnpip_fid_beta_dE_woK0_n"),100,-50,50,140,1,1.7);
  Tools::newTH2F( Form("CDHz_IMnpim_fid_beta_dE_woK0_n"),100,-50,50,140,1,1.7);
  Tools::newTH2F( Form("CDH_diffpos_z_n"),1000,-50,50,1000,-50,50);

  //** forward counters **//
  Tools::newTH1F( Form("mul_BVC"), 9, -0.5, 8.5 );
  Tools::newTH1F( Form("mul_CVC"), 11, -0.5, 10.5 );
  Tools::newTH1F( Form("mul_PC"), 11, -0.5, 10.5 );
  Tools::newTH1F( Form("mul_NC"), 11, -0.5, 10.5 );
  //target fiducial 
  Tools::newTH2F( Form("Vtx_ZX"),1000,-25,25,500,-12.5,12.5);
  Tools::newTH2F( Form("Vtx_ZY"),1000,-25,25,500,-12.5,12.5);
  Tools::newTH2F( Form("Vtx_XY"),500,-12.5,12.5,500,-12.5,12.5);
  Tools::newTH2F( Form("Vtx_ZX_primfid"),1000,-25,25,500,-12.5,12.5);
  Tools::newTH2F( Form("Vtx_ZY_primfid"),1000,-25,25,500,-12.5,12.5);
  Tools::newTH2F( Form("Vtx_XY_primfid"),500,-12.5,12.5,500,-12.5,12.5);
  Tools::newTH2F( Form("Vtx_ZX_fid"),1000,-25,25,500,-12.5,12.5);
  Tools::newTH2F( Form("Vtx_ZY_fid"),1000,-25,25,500,-12.5,12.5);
  Tools::newTH2F( Form("Vtx_XY_fid"),500,-12.5,12.5,500,-12.5,12.5);
  Tools::newTH2F( Form("Vtx_ZX_nofid"),1000,-25,25,500,-12.5,12.5);
  Tools::newTH2F( Form("Vtx_ZY_nofid"),1000,-25,25,500,-12.5,12.5);
  Tools::newTH2F( Form("Vtx_XY_nofid"),500,-12.5,12.5,500,-12.5,12.5);
  Tools::newTH2F( Form("CDHtime"),36,0.5,36.5,400,0,200);
  Tools::newTH2F( Form("CDHdE"),36,0.5,36.5,100,0,10);
  Tools::newTH2F( Form("CDHdE_wt"),36,0.5,36.5,100,0,10);
  Tools::newTH2F( Form("CDHdE_pippim"),36,0.5,36.5,100,0,10);
  Tools::newTH2F( Form("CDHdE_woK0_wSid"),36,0.5,36.5,100,0,10);
  Tools::newTH2F( Form("CDHdE_woK0_wSmid"),36,0.5,36.5,100,0,10);
  Tools::newTH2F( Form("CDHdE_woK0_wSid_n"),36,0.5,36.5,100,0,10);
  Tools::newTH2F( Form("dE_CDHtime"),150,0.,150,100,0,50);
  Tools::newTH2F( Form("dE_CDHtime_2track"),150,0.,150,100,0,50);
}

void InitIMPiSigmaHist()
{
  
  //CDH
  Tools::newTH2F( Form("dE_CDHtime_pippimn"),150,0.,150,100,0,50);
  
  for(int it0=0;it0<5;it0++){
    for(int iseg=0;iseg<36;iseg++){
      Tools::newTH1F( Form("CDH%d_T0%d_TOF_Neutral",iseg+1,it0+1), 400, 0, 100);
    }
  }
  Tools::newTH2F( Form("NeutraltimeEnergy"),100,0,100,150,0,150);
  Tools::newTH2F( Form("CDHzNeutraltime"),100,-50,50,150,0,150);
  Tools::newTH2F( Form("NMomCDHtime"),300,0,150,150,0,1.5);
  Tools::newTH2F( Form("NMomCDHtime_wK0"),300,0,150,150,0,1.5);
  for(int iseg=0;iseg<36;iseg++){
    Tools::newTH2F( Form("CDH%dzNeutraltime",iseg+1),100,-50,50,150,0,150);
    Tools::newTH2F( Form("NMomCDHtime%d",iseg+1),300,0,150,150,0,1.5);
    Tools::newTH2F( Form("NMomCDHtime%d_wK0",iseg+1),300,0,150,150,0,1.5);
  }
  
  Tools::newTH2F( Form("ntof_nlen"),1000,0,100,1000,0,100);

  Tools::newTH1F( Form("diff_CDH"), 73, -36.5, 36.5 );
  Tools::newTH1F( Form("diff_CDH_pippim"), 73, -36.5, 36.5 );
  Tools::newTH2F( Form("diff2D_CDH"), 73, -36.5, 36.5,500,-50,50 );
  Tools::newTH2F( Form("diff2D_CDH_pippim"), 73, -36.5, 36.5,500,-50,50 );
  Tools::newTH1F( Form("diff_CDH_pip"), 73, -36.5, 36.5 );
  Tools::newTH1F( Form("diff_CDH_pim"), 73, -36.5, 36.5 );
  Tools::newTH1F( Form("diff_CDH_CDC"), 181, 0, 181 );
  //Tools::newTH2F( Form("diff2D_CDH_CDC"), 181, 0, 181, 500,-50,50);
  Tools::newTH1F( Form("diff_CDH_CDC_pip"), 181, 0, 181 );
  Tools::newTH1F( Form("diff_CDH_CDC_pim"), 181, 0, 181 );
  
  Tools::newTH2F("diff2d_CDC_CDH_pim",100,-1.*TMath::Pi(),TMath::Pi(),100,-100,100);
  Tools::newTH2F("diff2d_CDC_CDH_pip",100,-1.*TMath::Pi(),TMath::Pi(),100,-100,100);

  //pi+ pi- X event Neutron ID
  Tools::newTH2F( Form("dE_betainv"), 500, 0, 50, 200, 0, 50);
  Tools::newTH2F( Form("dE_betainv_fid"), 500, 0, 50, 200, 0, 50);
  Tools::newTH2F( Form("dE_betainv_fid_beta"), 500, 0, 50, 200, 0, 50);
  Tools::newTH2F( Form("dE_betainv_fid_beta_dE"), 500, 0, 50, 200, 0, 50);
  Tools::newTH2F( Form("dE_betainv_fid_beta_dE_woK0"), 500, 0, 50, 200, 0, 50);
  Tools::newTH2F( Form("dE_betainv_fid_beta_dE_wK0"), 500, 0, 50, 200, 0, 50);
  Tools::newTH2F( Form("dE_betainv_fid_beta_dE_woK0_n"), 500, 0, 50, 200, 0, 50);
  Tools::newTH2F( Form("dE_MMom_fid_beta_woK0"), 100, 0, 1.5, 200, 0, 50);
  Tools::newTH2F( Form("dE_MMass_fid_beta_woK0"), 140, 0.4, 1.8, 200, 0, 50);
  
  //Miss mom. vs Miss mass
  Tools::newTH2F( Form("MMom_MMass"), 140, 0.4, 1.8, 100, 0, 1.5 );
  Tools::newTH2F( Form("MMom_MMass_fid"), 140, 0.4, 1.8, 100, 0, 1.5 );
  Tools::newTH2F( Form("MMom_MMass_fid_beta"), 140, 0.4, 1.8, 100, 0, 1.5 );
  Tools::newTH2F( Form("MMom_MMass_fid_beta_dE"), 140, 0.4, 1.8, 100, 0, 1.5 );
  Tools::newTH2F( Form("MMom_MMass_fid_beta_dE_woK0"), 140, 0.4, 1.8, 100, 0, 1.5 );
  Tools::newTH2F( Form("MMom_MMass_fid_beta_dE_woK0_wSid"), 140, 0.4, 1.8, 100, 0, 1.5 );
  Tools::newTH2F( Form("MMom_MMass_fid_beta_dE_wK0"), 140, 0.4, 1.8, 100, 0, 1.5 );
  Tools::newTH2F( Form("MMom_MMass_fid_beta_dE_woK0_n"), 140, 0.4, 1.8, 100, 0, 1.5 );
  //
  Tools::newTH1F( Form("IMpipi_dE"), 200, 0.4, 0.6 );
  Tools::newTH2F( Form("IMpipi_NMom_dE"),100,0, 1.5, 200, 0.4,0.6);
  Tools::newTH2F( Form("NMom_NMom_fid_beta_dE_woK0_n"), 100, 0, 1.5, 100, 0, 1.5 );
  Tools::newTH2F( Form("IMnpim_IMnpip_dE_woK0"), 140, 1, 1.7, 140, 1, 1.7 );
  Tools::newTH2F( Form("IMnpim_IMnpip_dE_woK0_n"), 140, 1, 1.7, 140, 1, 1.7 );
  Tools::newTH2F( Form("IMmnpim_IMmnpip_woK0_wSid_n"), 70, 1, 1.7, 70, 1, 1.7 );
  Tools::newTH2F( Form("MMnpip_MMnpim_woK0_wSid_n"), 70, 1, 1.7, 70, 1, 1.7 );
  Tools::newTH2F( Form("Cosn_IMnpipi_woK0_wSid_n"), 100, 1, 2, 50, -1, 1 );
  Tools::newTH1F( Form("IMnpipi_n"), 100, 1, 2 );
  Tools::newTH1F( Form("IMnpipi_wSid_n"), 100, 1, 2 );
  Tools::newTH1F( Form("IMnpipi_woK0_wSid_n"), 100, 1, 2 );
  Tools::newTH2F( Form("dE_IMnpipi_woK0_wSid_n"), 100, 1, 2, 200, 0, 50);
  Tools::newTH2F( Form("MMnmiss_IMnpipi_n"),100,1,2,140,0.4,1.8);
  Tools::newTH2F( Form("MMnmiss_IMnpipi_wSid_n"),100,1,2,140,0.4,1.8);
  Tools::newTH2F( Form("MMnmiss_IMnpipi_woK0_wSid_n"),100,1,2,140,0.4,1.8);
  Tools::newTH2F( Form("nmom_IMnpipi_woK0_wSid_n"),80,1.2,2,100,0,1.0);
  Tools::newTH2F( Form("q_IMnpipi_woK0_wSid_n"),80,1.2,2,300,0,1.5);
  Tools::newTH1F( Form("DCA_pip"), 3000, 0, 30 );
  Tools::newTH1F( Form("DCA_pim"), 3000, 0, 30 );
  Tools::newTH1F( Form("DCA_pip_SigmaP"), 3000, 0, 30 );
  Tools::newTH1F( Form("DCA_pim_SigmaP"), 3000, 0, 30 );
  Tools::newTH1F( Form("DCA_pippim_SigmaP"), 3000, 0, 30 );
  Tools::newTH1F( Form("DCA_pip_SigmaM"), 3000, 0, 30 );
  Tools::newTH1F( Form("DCA_pim_SigmaM"), 3000, 0, 30 );
  Tools::newTH1F( Form("DCA_pippim_SigmaM"), 3000, 0, 30 );
  Tools::newTH1F( Form("DCA_pippim"), 3000, 0, 30);
  Tools::newTH1F( Form("DCA_pip_SigmaPM"), 3000, 0, 30);
  Tools::newTH1F( Form("DCA_pim_SigmaPM"), 3000, 0, 30);
  Tools::newTH1F( Form("DCA_pippim_SigmaPM"), 3000, 0, 30);

  Tools::newTH2F( Form("KFchi2_vs"),100,0,100,100,0,100);
  Tools::SetXTitleH2(Form("KFchi2_vs"),"chi2/NDF S+");
  Tools::SetYTitleH2(Form("KFchi2_vs"),"chi2/NDF S-");
  Tools::newTH1F( Form("KF_decision"), 2, -0.5, 1.5 );
}




void InitIMLambdaPimHist() {
  Tools::newTH2F( Form("MMom_MMass"), 140, 0.4, 1.8, 100, 0, 1.5 );
  Tools::newTH2F( Form("MMom_MMass_fid"),140, 0.4, 1.8, 100, 0, 1.5);
  Tools::newTH2F( Form("MMom_MMass_fid_p"),140, 0.4, 1.8, 100, 0, 1.5);
  Tools::newTH2F( Form("MMom_MMass_fid_wL"),140, 0.4, 1.8, 100, 0, 1.5);
  Tools::newTH2F( Form("MMom_PMom_fid"),100, 0, 1.5, 100, 0, 1.5);
  Tools::newTH2F( Form("MMom_PMom_fid_p"),100, 0, 1.5, 100, 0, 1.5);
  Tools::newTH2F( Form("IMppim1_IMppim2"),600, 1, 2.5, 600, 1, 2.5);
  Tools::newTH2F( Form("IMppim1_IMppim2_p"),600, 1, 2.5, 600, 1, 2.5);
  Tools::newTH1F( Form("DCA_pim1"),3000,0,30);
  Tools::newTH1F( Form("DCA_pim2"),3000,0,30);
  Tools::newTH1F( Form("DCA_pim1p"),3000,0,30);
  Tools::newTH1F( Form("DCA_pim2p"),3000,0,30);
  Tools::newTH1F( Form("DCA_pim1pim2"),3000,0,30);
  Tools::newTH2F( Form("q_IMppipi"),100,1,2,300,0,1.5);
  Tools::newTH2F( Form("q_IMppipi_p"),100,1,2,300,0,1.5);
  Tools::newTH2F( Form("q_IMppipi_wL_p"),100,1,2,300,0,1.5);
}
