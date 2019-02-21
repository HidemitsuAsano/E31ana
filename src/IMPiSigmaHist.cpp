#include "IMPiSigmaHist.h"
#include "Tools.h"

void InitBasicHist(const bool MCFlag)
{
  // geneneral informantion **//
  Tools::newTH1F( Form("Time"), 3000, -0.5, 2999.5 );
  Tools::newTH1F( Form("nGoodTrack"), 11, -0.5, 10.5 );
  Tools::newTH1F( Form("Scaler"), 41, -0.5, 40.5 );

  // CDC and CDH information from CDC-trackig file **//
  Tools::newTH1F( Form("mul_CDH"),Form("CDH multiplicity"), 11, -0.5, 10.5 );
  Tools::newTH1F( Form("mul_CDH_assoc"), 11, -0.5, 10.5 );
  Tools::newTH1F( Form("npimangle"),628, 0, 2*3.14);
  Tools::newTH1F( Form("npipangle"),628, 0, 2*3.14);

  if(MCFlag){
    //vertex reso//
    Tools::newTH1F( Form("vertex_diff_X"), 100, -50, 50 );
    Tools::newTH1F( Form("vertex_diff_Y"), 100, -50, 50 );
    Tools::newTH1F( Form("vertex_diff_Z"), 100, -50, 50 );

    //===== CDH hit pos study =====//
    Tools::newTH2F( Form("CDH_mom_diffpos_pi_phi"), 100, -10, 10, 100, 0, 1.0 );
    Tools::newTH2F( Form("CDH_mom_diffpos_pi_z"), 100, -10, 10, 100, 0, 1.0 );
    Tools::newTH2F( Form("CDH_mom_TOF_pi"),       100, -2, 2, 100, 0, 1.0 );
  }
  //** beam line **//
  Tools::newTH1F( Form("mul_BHD"), 12, -0.5, 11.5 );
  Tools::newTH1F( Form("mul_T0"),   6, -0.5, 5.5 );
  Tools::newTH1F( Form("tof_T0BHD"), 2000, 20, 40 );
  Tools::newTH1F( Form("tracktime_BPC"),  1200, -200, 400 );
  Tools::newTH1F( Form("trackchi2_BPC"),  200, 0, 20 );
  Tools::newTH1F( Form("ntrack_BPC"),  6, -0.5, 5.5 );
  Tools::newTH1F( Form("tracktime_BLC1"), 1200, -200, 400 );
  Tools::newTH1F( Form("tracktime_BLC2"), 1200, -200, 400 );
  Tools::newTH1F( Form("trackchi2_BLC1"), 200, 0, 20 );
  Tools::newTH1F( Form("trackchi2_BLC2"), 200, 0, 20 );
  Tools::newTH1F( Form("ntrack_BLC1"), 6, -0.5, 5.5 );
  Tools::newTH1F( Form("ntrack_BLC2"), 6, -0.5, 5.5 );
  Tools::newTH2F( Form("dydx_BLC2BPC"),     130, -1.3, 1.3, 130, -1.3, 1.3 );
  Tools::newTH2F( Form("dydzdxdz_BLC2BPC"), 175, -0.035, 0.035, 175, -0.035, 0.035 );
  Tools::newTH1F( Form("trackchi2_beam"), 400, 0, 40 );
  Tools::newTH1F( Form("momentum_beam"), 180, 0.92, 1.10 );
  Tools::newTH1F( Form("PID_beam"),5,-1.5,3.5);
  //** CDS **//
  Tools::newTH1F( Form("trackchi2_CDC"), 1000, 0, 50 );
  Tools::newTH2F( Form("PID_CDS_beta"), 2000, 0, 10., 1000, -1.2, 1.2 );
  Tools::newTH2F( Form("PID_CDS"), 1000, -0.6, 5, 1000, -1.2, 1.2 );
  Tools::newTH1F( Form("ntrack_CDS"), 6, -0.5, 5.5 );
  Tools::newTH1F( Form("ntrack_pi_plus"), 6, -0.5, 5.5 );
  Tools::newTH1F( Form("ntrack_proton"), 6, -0.5, 5.5 );
  //Tools::newTH1F( Form("ntrack_deuteron"), 6, -0.5, 5.5 );
  Tools::newTH1F( Form("ntrack_pi_minus"), 6, -0.5, 5.5 );
  Tools::newTH1F( Form("ntrack_K_minus"), 6, -0.5, 5.5 );

  //** forward counters **//
  Tools::newTH1F( Form("mul_BVC"), 9, -0.5, 8.5 );
  Tools::newTH1F( Form("mul_CVC"), 11, -0.5, 10.5 );
  Tools::newTH1F( Form("mul_PC"), 11, -0.5, 10.5 );
  //target fiducial 
  Tools::newTH2F( Form("Vtx_ZX"),1000,-25,25,500,-12.5,12.5);
  Tools::newTH2F( Form("Vtx_ZY"),1000,-25,25,500,-12.5,12.5);
  Tools::newTH2F( Form("Vtx_XY"),500,-12.5,12.5,500,-12.5,12.5);
  Tools::newTH2F( Form("Vtx_ZX_fid"),1000,-25,25,500,-12.5,12.5);
  Tools::newTH2F( Form("Vtx_ZY_fid"),1000,-25,25,500,-12.5,12.5);
  Tools::newTH2F( Form("Vtx_XY_fid"),500,-12.5,12.5,500,-12.5,12.5);
  Tools::newTH1F( Form("CDHNeutralSeg"),36, 0.5, 36.5);
  Tools::newTH2F( Form("Vtx_ZX_nofid"),1000,-25,25,500,-12.5,12.5);
  Tools::newTH2F( Form("Vtx_ZY_nofid"),1000,-25,25,500,-12.5,12.5);
  Tools::newTH2F( Form("Vtx_XY_nofid"),500,-12.5,12.5,500,-12.5,12.5);
}

void InitIMPiSigmaHist()
{
  
  //CDH
  Tools::newTH2F( Form("CDHtime"),36,0.5,36.5,400,0,200);
  Tools::newTH2F( Form("dE_CDHtime"),100,0.,100,100,0,50);
  Tools::newTH2F( Form("dE_CDHtime_2track"),100,0.,100,100,0,50);
  Tools::newTH2F( Form("dE_CDHtime_pippimn"),100,0.,100,100,0,50);

  Tools::newTH2F( Form("NeutraltimeEnergy"),100,0,100,200,0,50);
  Tools::newTH1F( Form("diff_CDH"), 73, -36.5, 36.5 );
  Tools::newTH1F( Form("diff_CDH_CDC"), 181, 0, 181 );
  
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
  Tools::newTH2F( Form("MMnmiss_IMnpipi_n"),100,1,2,100,0,1.5);
  Tools::newTH2F( Form("MMnmiss_IMnpipi_wSid_n"),100,1,2,100,0,1.5);
  Tools::newTH2F( Form("MMnmiss_IMnpipi_woK0_wSid_n"),100,1,2,100,0,1.5);
  Tools::newTH2F( Form("nmom_IMnpipi_woK0_wSid_n"),100,1,2,100,0,1.0);
  Tools::newTH2F( Form("q_IMnpipi_woK0_wSid_n"),100,1,2,300,0,1.5);
  Tools::newTH1F( Form("DCA_pip"), 3000, 0, 30 );
  Tools::newTH1F( Form("DCA_pim"), 3000, 0, 30 );
  Tools::newTH1F( Form("DCA_pip_SigmaP"), 500, 0, 5 );
  Tools::newTH1F( Form("DCA_pim_SigmaP"), 500, 0, 5 );
  Tools::newTH1F( Form("DCA_pip_SigmaM"), 500, 0, 5 );
  Tools::newTH1F( Form("DCA_pim_SigmaM"), 500, 0, 5 );
  Tools::newTH1F( Form("DCA_pippim"), 3000, 0, 30);

  Tools::newTH2F( Form("KFchi2_vs"),100,0,100,100,0,100);
  Tools::SetXTitleH2(Form("KFchi2_vs"),"chi2/NDF S+");
  Tools::SetYTitleH2(Form("KFchi2_vs"),"chi2/NDF S-");
  Tools::newTH1F( Form("KF_decision"), 2, -0.5, 1.5 );//TODO implement
}




