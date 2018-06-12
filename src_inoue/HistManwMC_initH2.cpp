#include "HistManwMC.h"

void HistManwMC::initHist()
{
  rtFile-> cd();
  std::cout<<"===== HistManwMC::initHist START  ====="<<std::endl;
  initHistCom();
  new TH1F("LP_num", "#Lambda p number", 10, -0.5, 9.5);

  new TH1F("KNpim_MM",     "p(K^{-}, n #pi^{-})\"X\"", 1000,  0.0, 1.0);
  new TH1F("KNpim_MM2",    "p(K^{-}, n #pi^{-})\"X\"", 2000, -1.0, 1.0);
  new TH1F("KNpim_MM_wSm", "p(K^{-}, n #pi^{-})\"X\"", 1000,  0.0, 1.0);
  new TH2F("KNpim_MM_Npim_IM", "p(K^{-}, n #pi^{-})\"X\" vs n #pi^{-}", 1000, 0.0, 1.0, 1000, 1.0, 2.0);
  new TH2F("Kpim_MM_Npim_IM", "p(K^{-}, #pi^{-})\"X\" vs n #pi^{-} IM", 2000, 0.0, 2.0, 1000, 1.0, 2.0);
  new TH2F("Kpim_MM_Npim_IM_mm_pip", "p(K^{-}, #pi^{-})\"X\" vs n #pi^{-} IM", 2000, 0.0, 2.0, 1000, 1.0, 2.0);
  new TH2F("KN_MM_Npim_IM_mm_pip",  "p(K^{-}, n)\"X\" vs n #pi^{-} w/ #pi^{+}", 1000, 0.0, 1.0, 1000, 1.0, 2.0);
  new TH2F("KN_MM_Npim_IM_mm_pip_wSm_woK0",  "p(K^{-}, n)\"X\" vs n #pi^{-} w/ #pi^{+}", 1000, 0.0, 1.0, 1000, 1.0, 2.0);
  new TH1F("Sm_ang_Npim", "#Sigma^{-} cos#theta", 1000, -1.0, 1.0);

  new TH1F("KNpip_MM",     "p(K^{-}, n #pi^{+})\"X\"", 1000,  0.0, 1.0);
  new TH1F("KNpip_MM2",    "p(K^{-}, n #pi^{+})\"X\"", 2000, -1.0, 1.0);
  new TH1F("KNpip_MM_wSp", "p(K^{-}, n #pi^{+})\"X\"", 1000,  0.0, 1.0);
  new TH2F("KNpip_MM_Npip_IM", "p(K^{+}, n #pi^{-})\"X\" vs n #pi^{+}", 1000, 0.0, 1.0, 1000, 1.0, 2.0);
  new TH2F("Kpip_MM_Npip_IM", "p(K^{+}, #pi^{+})\"X\" vs n #pi^{+} IM", 2000, 0.0, 2.0, 1000, 1.0, 2.0);
  new TH2F("Kpip_MM_Npip_IM_mm_pim", "p(K^{+}, #pi^{+})\"X\" vs n #pi^{+} IM", 2000, 0.0, 2.0, 1000, 1.0, 2.0);
  new TH2F("KN_MM_Npip_IM_mm_pim",  "p(K^{-}, n)\"X\" vs n #pi^{+} w/ #pi^{-}", 1000, 0.0, 1.0, 1000, 1.0, 2.0);
  new TH2F("KN_MM_Npip_IM_mm_pim_wSp_woK0",  "p(K^{-}, n)\"X\" vs n #pi^{+} w/ #pi^{-}", 1000, 0.0, 1.0, 1000, 1.0, 2.0);
  new TH1F("Sp_ang_Npip", "#Sigma^{+} cos#theta", 1000, -1.0, 1.0);

  new TH1F("KPpipi_MM_noCut",      "p(K^{-}, #pi^{+} #pi^{-})\"X\"", 2000, 0.0, 2.0);
  new TH1F("KPpipi_MM_noCut_wp",   "p(K^{-}, #pi^{+} #pi^{-})\"X\"", 2000, 0.0, 2.0);
  new TH1F("KPpipi_MM_noCut_woK0", "p(K^{-}, #pi^{+} #pi^{-})\"X\"", 2000, 0.0, 2.0);
  new TH1F("CDS_IM_pipi_wL",        "CDS #pi^{+} #pi^{-} IM",        1000, 0.0, 1.0);
  new TH2F("KPpip_KPpim_MM_mmL",    "p(K^{-}, #pi^{+})\"X\" vs p(K^{-}, #pi^{-})\"X\"", 1000, 1.0, 2.0, 1000, 1.0, 2.0);
  new TH2F("KPpip_KPpim_MM_mmL_woK0", "p(K^{-}, #pi^{+})\"X\" vs p(K^{-}, #pi^{-})\"X\"", 1000, 1.0, 2.0, 1000, 1.0, 2.0);
  new TH2F("KPpipi_mom_mmL", "#Lambda mom", 1000, -1.0, 1.0, 2000, 0.0, 2.0);
  new TH2F("KPpipi_mom_wp_mmL", "#Lambda mom", 1000, -1.0, 1.0, 2000, 0.0, 2.0);
  new TH1F("KPpipip_MM2", "p(K^{-}, p #pi^{+} #pi^{-} p)MM^{2}", 2000, -1.0, 1.0);
  new TH2F("KPpipip_MM2_KPpipi_PL_MM_mmL", "p(K^{-}, p #pi^{+} #pi^{-} p)MM^{2} vs p(\"#Lambda\", p)\"X\"", 2000, -1.0, 1.0, 2000, 0.0, 2.0);

  new TH1F("KPpipi_PL_MM_vtx", "p(#Lambda, p)\"X\"", 2000, 0.0, 2.0);
  new TH2F("LP_vtx_XY", "#Lambda p vertex XY Plane", 1000, -25, 25, 1000, -25, 25);
  new TH2F("LP_vtx_ZX", "#Lambda p vertex ZX Plane", 1000, -25, 25, 1000, -25, 25);
  new TH2F("LP_vtx_ZY", "#Lambda p vertex ZY Plane", 1000, -25, 25, 1000, -25, 25);
  new TH1F("DCA_LP",    "#Lambda p DCA",             1000, 0.,   10);
  new TH1F("FL_L",      "#Lambda flight length",     1000, -25., 25);

  new TH1F("KN_MM_pipi_wN_NC0", "d(K^{-}, n)\"X\"", 500, 0.0, 2.0);
  new TH1F("KN_MM_pipi_wN_NC1", "d(K^{-}, n)\"X\"", 500, 0.0, 2.0);
  new TH1F("KN_MM_pipi_wN_NC2", "d(K^{-}, n)\"X\"", 500, 0.0, 2.0);
  new TH1F("KN_MM_pipi_wN_NC3", "d(K^{-}, n)\"X\"", 500, 0.0, 2.0);
  new TH1F("KN_MM_pipi_wN_NC4", "d(K^{-}, n)\"X\"", 500, 0.0, 2.0);
  new TH1F("KN_MM_pipi_wN_NC5", "d(K^{-}, n)\"X\"", 500, 0.0, 2.0);
  new TH1F("KN_MM_pipi_wN_NC6", "d(K^{-}, n)\"X\"", 500, 0.0, 2.0);

  new TH1F("KN_MM_pipi_wK0_wN_NC0", "d(K^{-}, n)\"X\"", 500, 0.0, 2.0);
  new TH1F("KN_MM_pipi_wK0_wN_NC1", "d(K^{-}, n)\"X\"", 500, 0.0, 2.0);
  new TH1F("KN_MM_pipi_wK0_wN_NC2", "d(K^{-}, n)\"X\"", 500, 0.0, 2.0);
  new TH1F("KN_MM_pipi_wK0_wN_NC3", "d(K^{-}, n)\"X\"", 500, 0.0, 2.0);
  new TH1F("KN_MM_pipi_wK0_wN_NC4", "d(K^{-}, n)\"X\"", 500, 0.0, 2.0);
  new TH1F("KN_MM_pipi_wK0_wN_NC5", "d(K^{-}, n)\"X\"", 500, 0.0, 2.0);
  new TH1F("KN_MM_pipi_wK0_wN_NC6", "d(K^{-}, n)\"X\"", 500, 0.0, 2.0);

  new TH2F("KN_MM_pipi_wK0_KPpipi_MM", "p(K^{-}, #pi^{-} #pi^{+})\"X\"", 1000, 0.0, 1.0, 2000, 0.0, 2.0); 

  new TH1F("KL_MM", "d(K^{-}, #Lambda)\"X\"",   1000, 0.0, 1.0);
  new TH1F("KL_MM2", "d(K^{-}, #Lambda)MM^{2}", 2000, -1.0, 1.0);
  new TH2F("KL_pi0_kin", "#pi^{0} cos#theta vs mom", 1000, -1.0, 1.0, 2000, 0.0, 2.0);
  new TH2F("Kpi0_L_kin", "#Lambda cos#theta vs mom", 1000, -1.0, 1.0, 2000, 0.0, 2.0);

  new TH1F("CDS_IM_pipi_wN",     "CDS #pi^{+} #pi^{-} IM", 1000, 0.0, 1.0);
  new TH1F("CDS_IM_pipi_wN_woS", "CDS #pi^{+} #pi^{-} IM", 1000, 0.0, 1.0);

  new TH1F("KPpipi_MM",         "p(K^{-}, #pi^{+} #pi^{-})\"X\"",                   2000, 0.0, 2.0);
  new TH1F("KPpipi_MM_wK0",     "p(K^{-}, #pi^{+} #pi^{-})\"X\"",                   2000, 0.0, 2.0);

  new TH2F("KPpip_KPpim_MM",           "p(K^{-}, #pi^{+})\"X\" vs p(K^{-}, #pi^{-})\"X\"", 1000, 1.0, 2.0, 1000, 1.0, 2.0);
  new TH2F("KPpip_KPpim_MM_wN",        "p(K^{-}, #pi^{+})\"X\" vs p(K^{-}, #pi^{-})\"X\"", 1000, 1.0, 2.0, 1000, 1.0, 2.0);
  new TH2F("KPpip_KPpim_MM_wN_woK0",   "p(K^{-}, #pi^{+})\"X\" vs p(K^{-}, #pi^{-})\"X\"", 1000, 1.0, 2.0, 1000, 1.0, 2.0);
  new TH2F("KPpip_KPpim_MM_wN_woK0_3", "p(K^{-}, #pi^{+})\"X\" vs p(K^{-}, #pi^{-})\"X\"", 1000, 1.0, 2.0, 1000, 1.0, 2.0);
  new TH2F("KPpip_KPpim_MM_wN_woK0_5", "p(K^{-}, #pi^{+})\"X\" vs p(K^{-}, #pi^{-})\"X\"", 1000, 1.0, 2.0, 1000, 1.0, 2.0);

  new TH2F("n_hitpos",      "Neutron Hit Position", 1000, -500, 500, 1000, -500, 500);
  new TH2F("n_hitpos_MC",   "Neutron Hit Position", 1000, -500, 500, 1000, -500, 500);
  new TH2F("n_hitpos_whit", "Neutron Hit Position", 1000, -500, 500, 1000, -500, 500);
  new TH2F("n_hitpos_wMChit", "Neutron Hit Position", 1000, -500, 500, 1000, -500, 500);
  new TH2F("n_hitpos_wBVC", "Neutron Hit Position", 1000, -500, 500, 1000, -500, 500);
  new TH2F("n_hitpos_wCVC", "Neutron Hit Position", 1000, -500, 500, 1000, -500, 500);

  new TH1F("K0_ang", "K0 cos#theta", 1000, -1.0, 1.0);
  new TH1F("K0_ang_nhit", "K0 cos#theta", 1000, -1.0, 1.0);
  new TH1F("Sm_ang", "#Sigma^{-} cos#theta", 1000, -1.0, 1.0);
  new TH1F("Sp_ang", "#Sigma^{+} cos#theta", 1000, -1.0, 1.0);

  new TH1F("K0_ang_CM", "K0 cos#theta", 1000, -1.0, 1.0);
  new TH1F("K0_ang_CM_nhit", "K0 cos#theta", 1000, -1.0, 1.0);
  new TH1F("Sm_ang_CM", "#Sigma^{-} cos#theta", 1000, -1.0, 1.0);
  new TH1F("Sp_ang_CM", "#Sigma^{+} cos#theta", 1000, -1.0, 1.0);

  new TH1F("Npip_IM_Sp", "n #pi^{+} IM", 1000, 1.0, 2.0);
  new TH1F("Npim_IM_Sm", "n #pi^{-} IM", 1000, 1.0, 2.0);

  new TH1F("NC_trig",  "NC efficiency trigger", 10, 0, 10);
  new TH1F("NC_eff",   "NC efficiency effective", 10, 0, 10);
  for( int i=0; i<10; i++ ){
    new TH1F(Form("NC_eff_ev%d", i), Form("NC eff event selection condition%d", i), 10, 0, 10);
  }
  new TH1F("NC_eff0",  "NC efficiency effective", 10, 0, 10);

  new TH2F("NC_eff2", "NC efficiency dE check", 10, -0.5, 9.5, 50, -0.5, 49.5);
  new TH2F("n_mom_diff", "neutron mom diff", 1000, 0.5, 1.5, 500, -0.1, 0.1); 

  for( int i=0; i<50; i++ ){
    new TH1F(Form("NC_overbeta_dE%d",i), Form("NC overbeta dE>%d", 2*i), 5000, 0.0, 5.0);
  }
  new TH1F("KN_MM_pipi_wK0_wN", "d(K^{-}, n)", 1000, 0.0, 1.0);

  if( simReader ){
    new TH2F("diff_n_pos", "n diff", 1000, -50, 50, 1000, -50, 50);
  }

  //*** for TTree ***//
  // fNpipiData = new NpipiData();
  // new TTree("NpipiTree", "NpipiTree");
  // TTree *tree = (TTree*)rtFile->Get("NpipiTree");
  // tree-> Branch("NpipiData", &fNpipiData);

#if CALIB
  if( !simReader ) initCalib();
#endif
  if( simReader ) simReader->initHist(rtFile);

  std::cout<<"===== HistManwMC::initHist FINISH ====="<<std::endl;
}
