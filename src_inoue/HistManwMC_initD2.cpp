#include "HistManwMC.h"

void HistManwMC::initHist()
{
  rtFile-> cd();
  std::cout<<"===== HistManwMC::initHist START  ====="<<std::endl;
  initHistCom();
  new TH1F("CDS_IM_Lpim_C_mm_p",        "CDS #Lambda #pi^{-} IM", 1000, 1.0, 2.0);
  new TH1F("CDS_IM_Lpim_C_mm_p_wFDC1",  "CDS #Lambda #pi^{-} IM", 1000, 1.0, 2.0);
  new TH1F("CDS_IM_Lpim_C_mm_p_FC",     "CDS #Lambda #pi^{-} IM", 1000, 1.0, 2.0);
  new TH1F("CDS_IM_Lpim_C_mm_p_fp",     "CDS #Lambda #pi^{-} IM", 1000, 1.0, 2.0);
  new TH1F("CDS_IM_Lpim_C_mm_p_wCVCPC", "CDS #Lambda #pi^{-} IM", 1000, 1.0, 2.0);
  new TH1F("Lpim_C_mm_p_mass2",         "FC mass^{2}", 2000, 0.0, 2.0);

  new TH2F("CDS_Lpim_p_ang_CDS_IM_Lpim", "\"p\" cos#theta vs #Lambda #pi^{-} IM", 500, -1, 1, 1000, 1.0, 2.0);
  new TH2F("CDS_Lpim_p_ang_CDS_IM_Lpim_CDH2", "\"p\" cos#theta vs #Lambda #pi^{-} IM", 500, -1, 1, 1000, 1.0, 2.0);
  new TH2F("CDS_Lpim_p_ang_CDS_IM_Lpim_C", "\"p\" cos#theta vs #Lambda #pi^{-} IM", 500, -1, 1, 1000, 1.0, 2.0);

  new TH2F("CDS_ppimpim_vtxL_XY",  "CDS L vtx in p #pi^{-} #pi^{-}", 500, -25, 25, 500, -25, 25);
  new TH2F("CDS_ppimpim_vtxL_ZX",  "CDS L vtx in p #pi^{-} #pi^{-}", 500, -50, 50, 500, -25, 25);
  new TH2F("CDS_ppimpim_vtxL_ZY",  "CDS L vtx in p #pi^{-} #pi^{-}", 500, -50, 50, 500, -25, 25);

  new TH2F("CDS_IM_ppim_ppim",     "CDS p #pi^{-} vs p #pi^{-} IM", 1000, 1.0, 2.0, 1000, 1.0, 2.0);
  new TH1F("CDS_IM_Lpim",          "CDS #Lambda #pi^{-} IM", 1000, 1.0, 2.0);
  new TH1F("KLpim_MM",             "d(K^{-}, #Lambda #pi^{-})\"X\"", 2000, 0.0, 2.0);
  new TH2F("CDS_IM_Lpim_KLpim_MM", "CDS #Lambda #pi^{-} vs d(K^{-}, #Lambda #pi^{-})\"X\"", 1000, 1.0, 2.0, 2000, 0.0, 2.0);
  new TH1F("CDS_Lpim_p_ang",       "\"p\"cos#theta", 1000, -1.0, 1.0);
  new TH1F("CDS_Lpim_p_ang_lab",   "\"p\"cos#theta", 1000, -1.0, 1.0);

  new TH2F("CDS_IM_ppim_ppim_CDH2",     "CDS p #pi^{-} vs p #pi^{-} IM", 1000, 1.0, 2.0, 1000, 1.0, 2.0);
  new TH1F("CDS_IM_Lpim_CDH2",          "CDS #Lambda #pi^{-} IM", 1000, 1.0, 2.0);
  new TH1F("KLpim_MM_CDH2",             "d(K^{-}, #Lambda #pi^{-})\"X\"", 2000, 0.0, 2.0);
  new TH2F("CDS_IM_Lpim_KLpim_MM_CDH2", "CDS #Lambda #pi^{-} vs d(K^{-}, #Lambda #pi^{-})\"X\"", 1000, 1.0, 2.0, 2000, 0.0, 2.0);
  new TH1F("CDS_Lpim_p_ang_CDH2",       "\"p\"cos#theta", 1000, -1.0, 1.0);
  new TH1F("CDS_Lpim_p_ang_lab_CDH2",   "\"p\"cos#theta", 1000, -1.0, 1.0);

  new TH2F("CDS_IM_ppim_ppim_C",     "CDS p #pi^{-} vs p #pi^{-} IM", 1000, 1.0, 2.0, 1000, 1.0, 2.0);
  new TH1F("CDS_IM_Lpim_C",          "CDS #Lambda #pi^{-} IM", 1000, 1.0, 2.0);
  new TH1F("KLpim_MM_C",             "d(K^{-}, #Lambda #pi^{-})\"X\"", 2000, 0.0, 2.0);
  new TH2F("CDS_IM_Lpim_KLpim_MM_C", "CDS #Lambda #pi^{-} vs d(K^{-}, #Lambda #pi^{-})\"X\"", 1000, 1.0, 2.0, 2000, 0.0, 2.0);
  new TH1F("CDS_Lpim_p_ang_C",       "\"p\"cos#theta", 1000, -1.0, 1.0);
  new TH1F("CDS_Lpim_p_ang_lab_C",   "\"p\"cos#theta", 1000, -1.0, 1.0);

  new TH1F("CDS_IM_Lpip", "CDS #Lambda #pi^{+}", 1000, 1.0, 2.0);
  new TH1F("KLpip_MM", "d(K^{-}, #Lambda #pi^{+})\"X\"", 1000, 1.0, 2.0);

  new TH1F("CDS_ppim_IM_wN", "CDS p #pi^{-} IM", 1000, 1.0, 2.0);
  new TH1F("KN_MM_wL", "d(K^{-}, n) w/ #Lambda", 2000, 0.0, 2.0);
  new TH2F("KN_MM_KNppim_MM",  "d(K^{-}, n)\"X\" vs d(K^{-}, n #Lambda)\"X\"", 2000, 0.0, 2.0, 1000,  0.0, 1.0);
  new TH2F("KN_MM_KNppim_MM2", "d(K^{-}, n)\"X\" vs d(K^{-}, n #Lambda)\"X\"", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("KN_MM_KNppim_MM_wL",  "d(K^{-}, n)\"X\" vs d(K^{-}, n #Lambda)\"X\"", 2000, 0.0, 2.0, 1000,  0.0, 1.0);
  new TH2F("KN_MM_KNppim_MM2_wL", "d(K^{-}, n)\"X\" vs d(K^{-}, n #Lambda)\"X\"", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("KN_MM_KNL_MM",  "d(K^{-}, n)\"X\" vs d(K^{-}, n #Lambda)\"X\"", 2000, 0.0, 2.0, 1000,  0.0, 1.0);
  new TH2F("KN_MM_KNL_MM2", "d(K^{-}, n)\"X\" vs d(K^{-}, n #Lambda)\"X\"", 2000, 0.0, 2.0, 2000, -1.0, 1.0);

  // 2015/12/16 Add **************************************************************//
  new TH2F("CDS_chi2_npipin", "CDS chi2", 250, 0.0, 50, 250, 0.0, 50);
  //  new TNtuple("tupNpipi", "tupNpipi", "mm_dkn:mm_dknpipi:mm_dknpim:mm_dknpip:im_pipi:im_npim:im_npip:n_mom:pim_mom:pip_mmom:pipi_mom:n_theta");
  new TH2F("KNpipi_MM_pipi_mom", "d(K^{-}, n)\"X\" vs #pi^{+} #pi^{-} mom", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  //===== for d(K-, n) slice =====//
  // new TH2F("KN_MM_pipi_wN_Npim_IM",         "d(K^{-}, n)\"X\" vs n #pi^{-} IM",                     2000, 0.0, 2.0, 1000, 1.0, 2.0);
  // new TH2F("KN_MM_pipi_wN_Npip_IM",         "d(K^{-}, n)\"X\" vs n #pi^{+} IM",                     2000, 0.0, 2.0, 1000, 1.0, 2.0);
  // new TH2F("KN_MM_pipi_wN_CDS_IM_pipi",     "d(K^{-}, n)\"X\" vs #pi^{+} #pi^{-} IM",               2000, 0.0, 2.0, 1000, 0.0, 1.0);
  // new TH2F("KN_MM_pipi_wN_KNpim_MM", "d(K^{-}, n)\"X\" vs d(K^{-}, n #pi^{-})\"X\"", 2000, 0.0, 2.0, 1000, 1.0, 2.0);
  // new TH2F("KN_MM_pipi_wN_KNpip_MM", "d(K^{-}, n)\"X\" vs d(K^{-}, n #pi^{+})\"X\"", 2000, 0.0, 2.0, 1000, 1.0, 2.0);
  // new TH2F("KN_MM_pipi_wN_n_theta",         "d(K^{-}, n)\"X vs \"n\" cos#theta", 2000, 0.0, 2.0, 1000, -1.0, 1.0);
  // new TH2F("KN_MM_pipi_wN_mmN_mom",         "d(K^{-}, n)\"X vs \"n\" mom", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  // new TH2F("KN_MM_pipi_wN_pim_mom",   "d(K^{-}, n)\"X\" vs #pi^{-} mom", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  // new TH2F("KN_MM_pipi_wN_pip_mom",   "d(K^{-}, n)\"X\" vs #pi^{+} mom", 2000, 0.0, 2.0, 1000, 0.0, 1.0);

  new TH2F("KN_MM_pipi_wK0_K0_theta",  "d(K^{-}, n)\"X\" vs K^{0} cos#theta",   2000, 0.0, 2.0, 1000, -1.0, 1.0);
  new TH2F("KN_MM_pipi_wSm_pip_theta", "d(K^{-}, n)\"X\" vs #pi^{+} cos#theta", 2000, 0.0, 2.0, 1000, -1.0, 1.0);
  new TH2F("KN_MM_pipi_wSp_pim_theta", "d(K^{-}, n)\"X\" vs #pi^{-} cos#theta", 2000, 0.0, 2.0, 1000, -1.0, 1.0);

  //===== for d(K-, n) w/pi+ pi- w/o K0 S_F =====//
  new TH2F("KN_MM_pipi_woAll_wN_KNpim_MM",  "d(K^{-}, n)\"X\" vs d(K^{-}, n #pi^{-})\"X\"",         2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("KN_MM_pipi_woAll_wN_KNpip_MM",  "d(K^{-}, n)\"X\" vs d(K^{-}, n #pi^{+})\"X\"",         2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("KN_MM_pipi_woAll_wN_n_theta",   "d(K^{-}, n)\"X\" vs \"n\" mom", 2000, 0.0, 2.0, 1000, -1.0, 1.0);
  new TH2F("KN_MM_pipi_woAll_wN_mmN_mom",   "d(K^{-}, n)\"X\" vs \"n\" mom", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("KN_MM_pipi_woAll_wN_pim_mom",   "d(K^{-}, n)\"X\" vs #pi^{-} mom", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("KN_MM_pipi_woAll_wN_pip_mom",   "d(K^{-}, n)\"X\" vs #pi^{+} mom", 2000, 0.0, 2.0, 1000, 0.0, 1.0);

  // 2015/12/12 Add **************************************************************//
  new TH1F("KN_MM_pipi_wN_woK0",       "d(K^{-}, n)\"X\"", 2000, 0.0, 2.0);
  new TH1F("KN_MM_pipi_wN_woK0_25",    "d(K^{-}, n)\"X\"", 2000, 0.0, 2.0);
  new TH1F("KN_MM_pipi_wN_woK0_3",     "d(K^{-}, n)\"X\"", 2000, 0.0, 2.0);
  new TH1F("KN_MM_pipi_wN_woK0_4",     "d(K^{-}, n)\"X\"", 2000, 0.0, 2.0);

  new TH1F("KN_MM_pipi_wN_woSm",       "d(K^{-}, n)\"X\"", 2000, 0.0, 2.0);
  new TH1F("KN_MM_pipi_wN_woSm_25",    "d(K^{-}, n)\"X\"", 2000, 0.0, 2.0);
  new TH1F("KN_MM_pipi_wN_woSm_3",     "d(K^{-}, n)\"X\"", 2000, 0.0, 2.0);
  new TH1F("KN_MM_pipi_wN_woSm_4",     "d(K^{-}, n)\"X\"", 2000, 0.0, 2.0);

  new TH1F("KN_MM_pipi_wN_woSp",       "d(K^{-}, n)\"X\"", 2000, 0.0, 2.0);
  new TH1F("KN_MM_pipi_wN_woSp_25",    "d(K^{-}, n)\"X\"", 2000, 0.0, 2.0);
  new TH1F("KN_MM_pipi_wN_woSp_3",     "d(K^{-}, n)\"X\"", 2000, 0.0, 2.0);
  new TH1F("KN_MM_pipi_wN_woSp_4",     "d(K^{-}, n)\"X\"", 2000, 0.0, 2.0);

  new TH1F("KN_MM_pipi_wN_woAll25", "d(K^{-}, n)\"X\"", 2000, 0.0, 2.0);
  new TH1F("KN_MM_pipi_wN_woAll3",  "d(K^{-}, n)\"X\"", 2000, 0.0, 2.0);
  new TH1F("KN_MM_pipi_wN_woAll4",  "d(K^{-}, n)\"X\"", 2000, 0.0, 2.0);

  // new TH2F("KN_MM_KNpipi_mom_woAll25", "d(K^{-}, n)\"X\"", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  // new TH2F("KN_MM_KNpipi_mom_woAll3",  "d(K^{-}, n)\"X\"", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  // new TH2F("KN_MM_KNpipi_mom_woAll4",  "d(K^{-}, n)\"X\"", 2000, 0.0, 2.0, 1000, 0.0, 1.0);

  // 2015/11/11 Add **************************************************************//
  new TH1F("KNpipi_MM_woAll_woF",     "d(K^{-}, n #pi^{-} #pi^{+})\"X\"",            2000, 0.0, 2.0);
  new TH1F("KNpipi_MM_woAll_vtx0",    "d(K^{-}, n #pi^{-} #pi^{+})\"X\"",            2000, 0.0, 2.0);
  new TH1F("KNpipi_MM_woAll_vtx1",    "d(K^{-}, n #pi^{-} #pi^{+})\"X\"",            2000, 0.0, 2.0);
  new TH1F("KNpipi_MM_woAll_vtx2",    "d(K^{-}, n #pi^{-} #pi^{+})\"X\"",            2000, 0.0, 2.0);
  new TH1F("KNpipi_MM_woAll_vtx3",    "d(K^{-}, n #pi^{-} #pi^{+})\"X\"",            2000, 0.0, 2.0);
  new TH1F("KNpipi_MM_woAll_vtx4",    "d(K^{-}, n #pi^{-} #pi^{+})\"X\"",            2000, 0.0, 2.0);
  new TH1F("KNpipi_MM_woAll_vtx5",    "d(K^{-}, n #pi^{-} #pi^{+})\"X\"",            2000, 0.0, 2.0);

  new TH1F("KN_MM_pipi_woAll_woF",  "d(K^{-}, n)\"X\"", 2000, 0.0, 2.0);
  new TH1F("KN_MM_pipi_woAll_vtx0", "d(K^{-}, n)\"X\"", 2000, 0.0, 2.0);
  new TH1F("KN_MM_pipi_woAll_vtx1", "d(K^{-}, n)\"X\"", 2000, 0.0, 2.0);
  new TH1F("KN_MM_pipi_woAll_vtx2", "d(K^{-}, n)\"X\"", 2000, 0.0, 2.0);
  new TH1F("KN_MM_pipi_woAll_vtx3", "d(K^{-}, n)\"X\"", 2000, 0.0, 2.0);
  new TH1F("KN_MM_pipi_woAll_vtx4", "d(K^{-}, n)\"X\"", 2000, 0.0, 2.0);
  new TH1F("KN_MM_pipi_woAll_vtx5", "d(K^{-}, n)\"X\"", 2000, 0.0, 2.0);

  new TH1F("DCA_pipi_beam_npipi",     "DCA #pi^{+} #pi^{-} & beam", 1000, 0.0, 50.0);
  new TH1F("DCA_pipi_beam_npipi_woF", "DCA #pi^{+} #pi^{-} & beam", 1000, 0.0, 50.0);

  new TH2F("KN_MM_KNpipi_mom",       "KN_MM_KNpipi_mom", 1000, 1.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("KN_MM_KNpipi_mom_wK0",   "KN_MM_KNpipi_mom", 1000, 1.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("KN_MM_KNpipi_mom_wSm",   "KN_MM_KNpipi_mom", 1000, 1.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("KN_MM_KNpipi_mom_wSp",   "KN_MM_KNpipi_mom", 1000, 1.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("KN_MM_KNpipi_mom_woAll", "KN_MM_KNpipi_mom", 1000, 1.0, 2.0, 1000, 0.0, 1.0);

  // 2015/11/04 Add **************************************************************//
  new TH1F("Vtx_Z_pipi_n",       "#pi^{+} #pi^{-} vtx Z",  1000, -50, 50);
  new TH1F("Vtx_Z_pipi_n_wN",    "#pi^{+} #pi^{-} vtx Z",  1000, -50, 50);
  new TH1F("Vtx_Z_pipi_n_tail",  "#pi^{+} #pi^{-} vtx Z",  1000, -50, 50);
  new TH2F("Vtx_XY_pipi_n",      "#pi^{+} #pi^{-} vtx XY", 1000, -20, 20, 1000, -20, 20);
  new TH2F("Vtx_XY_pipi_n_wN",   "#pi^{+} #pi^{-} vtx XY", 1000, -20, 20, 1000, -20, 20);
  new TH2F("Vtx_XY_pipi_n_tail", "#pi^{+} #pi^{-} vtx XY", 1000, -20, 20, 1000, -20, 20);

  new TH1F("Vtx_Z_pipi_n_woF",       "#pi^{+} #pi^{-} vtx Z",  1000, -50, 50);
  new TH1F("Vtx_Z_pipi_n_wN_woF",    "#pi^{+} #pi^{-} vtx Z",  1000, -50, 50);
  new TH1F("Vtx_Z_pipi_n_tail_woF",  "#pi^{+} #pi^{-} vtx Z",  1000, -50, 50);
  new TH2F("Vtx_XY_pipi_n_woF",      "#pi^{+} #pi^{-} vtx XY", 1000, -20, 20, 1000, -20, 20);
  new TH2F("Vtx_XY_pipi_n_wN_woF",   "#pi^{+} #pi^{-} vtx XY", 1000, -20, 20, 1000, -20, 20);
  new TH2F("Vtx_XY_pipi_n_tail_woF", "#pi^{+} #pi^{-} vtx XY", 1000, -20, 20, 1000, -20, 20);
  //******************************************************************************//
  new TH2F("Vtx_XY_npipi", "Vertex XY Plane", 1000, -30, 30, 1000,-30, 30);
  new TH2F("Vtx_ZX_npipi", "Vertex ZX Plane", 1000, -50, 50, 1000,-30, 30);
  new TH2F("Vtx_ZY_npipi", "Vertex ZY Plane", 1000, -50, 50, 1000,-30, 30);

  new TH1F("KCDSpim_MM", "d(K^{-}, #pi^{-})\"X\"", 2000,  0.0, 2.0);
  new TH1F("KCDSpip_MM", "d(K^{-}, #pi^{+})\"X\"", 2000,  0.0, 2.0);
  new TH1F("KCDSkm_MM",  "d(K^{-}, K^{-})\"X\"",   2000,  0.0, 2.0);
  new TH1F("KCDSp_MM",   "d(K^{-}, p)\"X\"",       2000,  0.0, 2.0);
  new TH1F("KCDSd_MM",   "d(K^{-}, d)\"X\"",       2000, -1.0, 1.0);

  new TH2F("KN_MM_woAll_wN_MM_Npim_IM", "d(K^{-}, n)\"X\" vs \"n\" #pi^{-} IM", 600, 0.0, 3.0, 200, 1.0, 2.0);
  new TH2F("KN_MM_woAll_wN_MM_Npip_IM", "d(K^{-}, n)\"X\" vs \"n\" #pi^{+} IM", 600, 0.0, 3.0, 200, 1.0, 2.0);

  //*************************//
  //*** for p p pi- event ***//
  //*************************//
  new TH1F("CDS_IM_ppim_2p", "CDS p #pi^{-} IM", 1000, 1.0, 2.0);
  new TH1F("CDS_IM_Lp", "CDS #Lambda p IM", 2000, 2.0, 4.0);
  new TH2F("CDS_IM_Lp_KCDSLp_MM", "CDS #Lambda p IM vs d(K~{-}, #Lambda p)\"X\"", 2000, 2.0, 4.0, 2000, -1.0, 1.0);
  new TH1F("Missing_mom_Lp", "Missing mom of d(K^{-}, #Lambda p)", 2000, 0.0, 2.0);
  new TH1F("Missing_ang_Lp", "Angle of d(K^{-}, #Lambda p)", 1000, -1.0, 1.0);
  new TH2F("Missing_ang_mom_Lp", "d(K^{-}, #Lambda p)\"X\" cos#theta vs mom", 1000, -1.0, 1.0, 2000, 0.0, 2.0);

  //**********************************//
  //*** for n pi- pi+ tagged event ***//
  //**********************************//
  //  new TNtuple("tupNpipi", "n pi+ pi- run check", "evnum");
  new TH1F("CDS_GoodTrack_npipi", "CDS Good Track n #pi^{+} #pi^{-} event", 10, -0.5, 9.5);

  new TH1F("KNpipi_status",      "n_{Forward} #pi^{+} #pi^{-} event status", 20, -0.5, 19.5);
  new TH1F("KNpipi_event",       "n_{Forward} #pi^{+} #pi^{-} event",        20, -0.5, 19.5);
  new TH1F("KNpipi_event25",     "n_{Forward} #pi^{+} #pi^{-} event",        20, -0.5, 19.5);
  new TH1F("KNpipi_event3",      "n_{Forward} #pi^{+} #pi^{-} event",        20, -0.5, 19.5);
  new TH1F("KNpipi_event4",      "n_{Forward} #pi^{+} #pi^{-} event",        20, -0.5, 19.5);
  new TH1F("KNpipi_woAll_event", "n_{Forward} #pi^{+} #pi^{-} event",        20, -0.5, 19.5);

  new TH1F("KNpipi_n_event",    "d(K^{-}, n #pi^{+} #pi^{-})\"n\" event", 20, -0.5, 19.5);
  new TH1F("KNpipi_n_event25",  "d(K^{-}, n #pi^{+} #pi^{-})\"n\" event", 20, -0.5, 19.5);
  new TH1F("KNpipi_n_event3",   "d(K^{-}, n #pi^{+} #pi^{-})\"n\" event", 20, -0.5, 19.5);
  new TH1F("KNpipi_n_event4",   "d(K^{-}, n #pi^{+} #pi^{-})\"n\" event", 20, -0.5, 19.5);

  new TH1F("CDS_IM_pipi_wN",        "CDS #pi^{+} #pi^{-} IM",       1000, 0.0, 1.0);
  new TH1F("CDS_IM_pipi_wN_wN",     "CDS #pi^{+} #pi^{-} IM",       1000, 0.0, 1.0);
  new TH1F("CDS_IM_pipi_wN_wN_wSp", "CDS #pi^{+} #pi^{-} IM",       1000, 0.0, 1.0);
  new TH1F("CDS_IM_pipi_wN_wN_wSm", "CDS #pi^{+} #pi^{-} IM",       1000, 0.0, 1.0);
  new TH1F("CDS_IM_pipi_wN_wN_wS",  "CDS #pi^{+} #pi^{-} IM",       1000, 0.0, 1.0);
  new TH1F("CDS_IM_pipi_wN_wSm",    "CDS #pi^{+} #pi^{-} IM",       1000, 0.0, 1.0);
  new TH1F("CDS_IM_pipi_wN_wSp",    "CDS #pi^{+} #pi^{-} IM",       1000, 0.0, 1.0);
  new TH1F("CDS_IM_pipi_wN_wSf",    "CDS #pi^{+} #pi^{-} IM",       1000, 0.0, 1.0);

  new TH1F("Npim_IM_pip",           "n #pi^{-} w/ #pi^{+}", 1000, 1.0, 2.0);
  new TH1F("Npim_IM_pip_wSp",       "n #pi^{-} w/ #pi^{+}", 1000, 1.0, 2.0);
  new TH1F("Npim_IM_pip_wK0",       "n #pi^{-} w/ #pi^{+}", 1000, 1.0, 2.0);
  new TH1F("Npim_IM_pip_wSp_K0",    "n #pi^{-} w/ #pi^{+}", 1000, 1.0, 2.0);
  new TH1F("Npim_IM_pip_wN",        "n #pi^{-} w/ #pi^{+}", 1000, 1.0, 2.0);
  new TH1F("Npim_IM_pip_wN_wK0",    "n #pi^{-} w/ #pi^{+}", 1000, 1.0, 2.0);
  new TH1F("Npim_IM_pip_wN_wSp",    "n #pi^{-} w/ #pi^{+}", 1000, 1.0, 2.0);
  new TH1F("Npim_IM_pip_wN_wK0_Sp", "n #pi^{-} w/ #pi^{+}", 1000, 1.0, 2.0);

  new TH1F("Npip_IM_pim",           "n #pi^{-} w/ #pi^{+}", 1000, 1.0, 2.0);
  new TH1F("Npip_IM_pim_wSm",       "n #pi^{-} w/ #pi^{+}", 1000, 1.0, 2.0);
  new TH1F("Npip_IM_pim_wK0",       "n #pi^{-} w/ #pi^{+}", 1000, 1.0, 2.0);
  new TH1F("Npip_IM_pim_wSm_K0",    "n #pi^{-} w/ #pi^{+}", 1000, 1.0, 2.0);
  new TH1F("Npip_IM_pim_wN",        "n #pi^{-} w/ #pi^{+}", 1000, 1.0, 2.0);
  new TH1F("Npip_IM_pim_wN_wK0",    "n #pi^{-} w/ #pi^{+}", 1000, 1.0, 2.0);
  new TH1F("Npip_IM_pim_wN_wSm",    "n #pi^{-} w/ #pi^{+}", 1000, 1.0, 2.0);
  new TH1F("Npip_IM_pim_wN_wK0_Sm", "n #pi^{-} w/ #pi^{+}", 1000, 1.0, 2.0);

  new TH2F("Npim_Npip_IM",        "n #pi^{-} IM vs #pi^{+} IM", 1000, 1.0, 2.0, 1000, 1.0, 2.0);
  new TH2F("Npim_Npip_IM_wN",     "n #pi^{-} IM vs #pi^{+} IM", 1000, 1.0, 2.0, 1000, 1.0, 2.0);
  new TH2F("Npim_Npip_IM_wN_wK0", "n #pi^{-} IM vs #pi^{+} IM", 1000, 1.0, 2.0, 1000, 1.0, 2.0);
  new TH2F("Npim_Npip_IM_wN_woK0","n #pi^{-} IM vs #pi^{+} IM", 1000, 1.0, 2.0, 1000, 1.0, 2.0);
  new TH2F("Npim_Npip_IM_wK0",    "n #pi^{-} IM vs #pi^{+} IM", 1000, 1.0, 2.0, 1000, 1.0, 2.0);
  new TH2F("Npim_Npip_IM_woK0",   "n #pi^{-} IM vs #pi^{+} IM", 1000, 1.0, 2.0, 1000, 1.0, 2.0);

  // for Vertex study 
  new TH2F("Vtx_bpim_bpip",             "beam #pi^{-} Vtx vs beam #pi^{+} Vtx", 1000, 0.0, 10, 1000, 0.0, 10);
  new TH2F("Vtx_bpim_bpip_woAll",       "beam #pi^{-} Vtx vs beam #pi^{+} Vtx", 1000, 0.0, 10, 1000, 0.0, 10);
  new TH2F("Vtx_bpim_bpip_woAll_wN",    "beam #pi^{-} Vtx vs beam #pi^{+} Vtx", 1000, 0.0, 10, 1000, 0.0, 10);
  new TH2F("Vtx_bpim_bpip_woAll_wN_wS", "beam #pi^{-} Vtx vs beam #pi^{+} Vtx", 1000, 0.0, 10, 1000, 0.0, 10);

  new TH1F("Vtxz_pim_pip",             "Vertex Z #pi^{-}-#pi^{+}", 2000, -10, 10);
  new TH1F("Vtxz_pim_pip_woAll",       "Vertex Z #pi^{-}-#pi^{+}", 2000, -10, 10);
  new TH1F("Vtxz_pim_pip_woAll_wN",    "Vertex Z #pi^{-}-#pi^{+}", 2000, -10, 10);
  new TH1F("Vtxz_pim_pip_woAll_wN_wS", "Vertex Z #pi^{-}-#pi^{+}", 2000, -10, 10);

  new TH1F("MinDCA_npipi", "Min DCA n pi pi event", 1000, 0.0, 10);

  //*** d(K-, n) Analysis *********************************************************************************************************//
  new TH2F("CDS_mass2_mom_pim_npipi", "CDS mass2 vs mom", 2000, -0.5, 9.5, 2000, -1.0, 1.0);
  new TH2F("CDS_mass2_mom_pip_npipi", "CDS mass2 vs mom", 2000, -0.5, 9.5, 2000, -1.0, 1.0);

  new TH2F("KNpim_mom", "d(K^{-}, n#pi^{-}) cos#theta vs mom", 1000, -1, 1, 1000, 0.0, 1.0);
  new TH2F("KNpip_mom", "d(K^{-}, n#pi^{+}) cos#theta vs mom", 1000, -1, 1, 1000, 0.0, 1.0);

  new TH2F("KNpim_mom_CM", "d(K^{-}, n#pi^{-}) cos#theta vs mom", 1000, -1, 1, 1000, 0.0, 1.0);
  new TH2F("KNpip_mom_CM", "d(K^{-}, n#pi^{+}) cos#theta vs mom", 1000, -1, 1, 1000, 0.0, 1.0);

  new TH2F("KN_mom", "d(K^{-}, n#pi^{+}) cos#theta vs mom", 1000, -1, 1, 1000, 0.0, 1.0);

  new TH2F("KN_mom_CM", "d(K^{-}, n#pi^{+}) cos#theta vs mom", 1000, -1, 1, 1000, 0.0, 1.0);


  new TH1F("MinMom_MMSA_npipi", "minimum momentum by MMSA", 1000, -1.0, 1.0);

  new TH1F("KN_MM_pipi",               "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_wK0",           "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_wK0_SB0",       "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_wK0_SB1",       "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_wK0_SB2",       "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_wK0_SB3",       "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_wK0_SB4",       "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_wK0_SB5",       "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_wK0_25",        "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_wK0_3",         "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_wK0_4",         "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_wS",            "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_wSm",           "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_wSp",           "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_rcSm",          "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_rcSp",          "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_wS_25",         "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_wS_3",          "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_wS_4",          "d(K^{-}, n)", 3000, 0.0, 3.0);

  new TH2F("KN_MM_pipi_wN_dE",        "d(K^{-}, n) vs dE", 1000, 1.0, 2.0, 500, 0.0, 200);
  new TH2F("KN_MM_pipi_wK0_wN_dE",    "d(K^{-}, n) vs dE", 1000, 1.0, 2.0, 500, 0.0, 200);
  for( int seg=1; seg<=5; seg++ ){
    new TH1F(Form("KN_MM_pipi_wN_T0%d", seg), "d(K^{-}, n)", 3000, 0.0, 3.0);
  }

  new TH1F("KN_MM_pipi_wN",           "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_wK0_wN",       "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_wK0_wN_SB0",   "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_wK0_wN_SB1",   "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_wK0_wN_SB2",   "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_wK0_wN_SB3",   "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_wK0_wN_SB4",   "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_wK0_wN_SB5",   "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_wS_wN",        "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_wSp_wN",       "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_wSm_wN",       "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_wSp_woSm_wN",  "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_wSm_woSp_wN",  "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_wAll_wN",      "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_woAll",        "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_woAll_wN",     "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_woAll_wN_wS",  "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_woAll2",       "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_woAll_woS2",   "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_woAll2_wN",    "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_woAll_woS2_wN", "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_rcSp_wN",       "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_rcSm_wN",       "d(K^{-}, n)", 3000, 0.0, 3.0);

  new TH1F("KNn_MM_npipi",           "\"n\"(K^{-}, n)", 2000, -1.0, 1.0);
  new TH1F("KNn_MM_npipi_cut",       "\"n\"(K^{-}, n)", 2000, -1.0, 1.0);
  new TH1F("KNn_MM_npipi_cut_wK0",   "\"n\"(K^{-}, n)", 2000, -1.0, 1.0);
  new TH1F("KNn_MM_npipi_cut_wSm",   "\"n\"(K^{-}, n)", 2000, -1.0, 1.0);
  new TH1F("KNn_MM_npipi_cut_wSp",   "\"n\"(K^{-}, n)", 2000, -1.0, 1.0);
  new TH1F("KNn_MM_npipi_cut_woAll", "\"n\"(K^{-}, n)", 2000, -1.0, 1.0);

  new TH1F("pim_ang_npimSp",    "#pi^{-} ang", 1000, -1.0, 1.0);
  new TH1F("pim_ang_npimSp_Kp", "#pi^{-} ang", 1000, -1.0, 1.0);
  new TH1F("pip_ang_npipSm",    "#pi^{-} ang", 1000, -1.0, 1.0);
  new TH1F("pip_ang_npipSm_Kp", "#pi^{-} ang", 1000, -1.0, 1.0);

  new TH1F("K0_ang_npipin",    "K^{0} cos#theta in Lab frame", 1000, -1.0, 1.0);
  new TH1F("K0_ang_npipin_Kp", "K^{0} cos#theta in K^{-}p CM", 1000, -1.0, 1.0);
  new TH1F("K0_ang_npipin_SB0",    "K^{0} cos#theta in Lab frame", 1000, -1.0, 1.0);
  new TH1F("K0_ang_npipin_Kp_SB0", "K^{0} cos#theta in K^{-}p CM", 1000, -1.0, 1.0);
  new TH1F("K0_ang_npipin_SB1",    "K^{0} cos#theta in Lab frame", 1000, -1.0, 1.0);
  new TH1F("K0_ang_npipin_Kp_SB1", "K^{0} cos#theta in K^{-}p CM", 1000, -1.0, 1.0);
  new TH1F("K0_ang_npipin_SB2",    "K^{0} cos#theta in Lab frame", 1000, -1.0, 1.0);
  new TH1F("K0_ang_npipin_Kp_SB2", "K^{0} cos#theta in K^{-}p CM", 1000, -1.0, 1.0);
  new TH1F("K0_ang_npipin_SB3",    "K^{0} cos#theta in Lab frame", 1000, -1.0, 1.0);
  new TH1F("K0_ang_npipin_Kp_SB3", "K^{0} cos#theta in K^{-}p CM", 1000, -1.0, 1.0);
  new TH1F("K0_ang_npipin_SB4",    "K^{0} cos#theta in Lab frame", 1000, -1.0, 1.0);
  new TH1F("K0_ang_npipin_Kp_SB4", "K^{0} cos#theta in K^{-}p CM", 1000, -1.0, 1.0);
  new TH1F("K0_ang_npipin_SB5",    "K^{0} cos#theta in Lab frame", 1000, -1.0, 1.0);
  new TH1F("K0_ang_npipin_Kp_SB5", "K^{0} cos#theta in K^{-}p CM", 1000, -1.0, 1.0);

  new TH1F("KN_MM_pipi_woAll_wN_wSp",  "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_woAll_wN_wSm",  "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_woAll_wN_wSp_woSm",  "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_woAll_wN_wSm_woSp",  "d(K^{-}, n)", 3000, 0.0, 3.0);

  new TH1F("KN_MM_pipi_woAll_wN0",    "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_woAll_wN0_wS", "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_woAll_wN1",    "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_woAll_wN1_wS", "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_woAll_wN2",    "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_woAll_wN2_wS", "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_woAll_wN0_wA", "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_woAll_wN0_wB", "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_woAll_wN0_wC", "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_woAll_wN0_wD", "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_woAll_wN_wA", "d(K^{-}, n)",  3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_woAll_wN_wB", "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_woAll_wN_wC", "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_woAll_wN_wD", "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_woAll_wN1_wA", "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_woAll_wN1_wB", "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_woAll_wN1_wC", "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_woAll_wN1_wD", "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_woAll_wN2_wA", "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_woAll_wN2_wB", "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_woAll_wN2_wC", "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_woAll_wN2_wD", "d(K^{-}, n)", 3000, 0.0, 3.0);

  new TH1F("KN_MM_pipi_woAll_wN0_wSm_up",   "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_woAll_wN0_wSm_down", "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_woAll_wN_wSm_up",    "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_woAll_wN_wSm_down",  "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_woAll_wN1_wSm_up",   "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_woAll_wN1_wSm_down", "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_woAll_wN2_wSm_up",   "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_woAll_wN2_wSm_down", "d(K^{-}, n)", 3000, 0.0, 3.0);

  new TH1F("KN_MM_pipi_woAll_wN0_wSp_up",   "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_woAll_wN0_wSp_down", "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_woAll_wN_wSp_up",    "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_woAll_wN_wSp_down",  "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_woAll_wN1_wSp_up",   "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_woAll_wN1_wSp_down", "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_woAll_wN2_wSp_up",   "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_woAll_wN2_wSp_down", "d(K^{-}, n)", 3000, 0.0, 3.0);

  //*** K- p threshold study ***//
  new TH1F("KNpipi_mom",       "d(K^{-}, n #pi^{+} #pi^{-})\"n\"", 1000, 0.0, 1.0);
  new TH1F("KNpipi_mom_wK0",    "d(K^{-}, n #pi^{+} #pi^{-})\"n\"", 1000, 0.0, 1.0);
  new TH1F("KNpipi_mom_wK0_SB0",    "d(K^{-}, n #pi^{+} #pi^{-})\"n\"", 1000, 0.0, 1.0);
  new TH1F("KNpipi_mom_wK0_SB1",    "d(K^{-}, n #pi^{+} #pi^{-})\"n\"", 1000, 0.0, 1.0);
  new TH1F("KNpipi_mom_wK0_SB2",    "d(K^{-}, n #pi^{+} #pi^{-})\"n\"", 1000, 0.0, 1.0);
  new TH1F("KNpipi_mom_wK0_SB3",    "d(K^{-}, n #pi^{+} #pi^{-})\"n\"", 1000, 0.0, 1.0);
  new TH1F("KNpipi_mom_wK0_SB4",    "d(K^{-}, n #pi^{+} #pi^{-})\"n\"", 1000, 0.0, 1.0);
  new TH1F("KNpipi_mom_wK0_SB5",    "d(K^{-}, n #pi^{+} #pi^{-})\"n\"", 1000, 0.0, 1.0);
  new TH1F("KNpipi_mom_wSm",    "d(K^{-}, n #pi^{+} #pi^{-})\"n\"", 1000, 0.0, 1.0);
  new TH1F("KNpipi_mom_wSp",    "d(K^{-}, n #pi^{+} #pi^{-})\"n\"", 1000, 0.0, 1.0);
  new TH1F("KNpipi_mom_woAll", "d(K^{-}, n #pi^{+} #pi^{-})\"n\"", 1000, 0.0, 1.0);

  new TH2F("KN_MM_pipi_KNpipi_MM",          "d(K^{-}, n)\"X\" vs d(K^{-}, n #pi^{+} #pi^{-})\"X\"", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("KN_MM_pipi_woAll_KNpipi_MM",    "d(K^{-}, n)\"X\" vs d(K^{-}, n #pi^{+} #pi^{-})\"X\"", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("KN_MM_pipi_wK0_KNpipi_MM",      "d(K^{-}, n)\"X\" vs d(K^{-}, n #pi^{+} #pi^{-})\"X\"", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("KN_MM_pipi_wSm_KNpipi_MM",      "d(K^{-}, n)\"X\" vs d(K^{-}, n #pi^{+} #pi^{-})\"X\"", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("KN_MM_pipi_wSp_KNpipi_MM",      "d(K^{-}, n)\"X\" vs d(K^{-}, n #pi^{+} #pi^{-})\"X\"", 2000, 0.0, 2.0, 2000, 0.0, 2.0);

  new TH2F("KN_MM_KNpim_MM_wN",           "d(K^{-}, n)\"X\" vs d(K^{-}, n #pi^{-})\"X\"", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("KN_MM_KNpip_MM_wN",           "d(K^{-}, n)\"X\" vs d(K^{-}, n #pi^{+})\"X\"", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("KN_MM_KNpim_MM_woAll_wN",     "d(K^{-}, n)\"X\" vs d(K^{-}, n #pi^{-})\"X\"", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("KN_MM_KNpip_MM_woAll_wN",     "d(K^{-}, n)\"X\" vs d(K^{-}, n #pi^{+})\"X\"", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("KN_MM_KNpim_MM_woAll_wN_200", "d(K^{-}, n)\"X\" vs d(K^{-}, n #pi^{-})\"X\"", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("KN_MM_KNpip_MM_woAll_wN_200", "d(K^{-}, n)\"X\" vs d(K^{-}, n #pi^{+})\"X\"", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("KN_MM_KNpim_mom_woAll",  "d(K^{-}, n)\"X\" vs d(K^{-}, n #pi^{-}) mom",  400, 0.0, 2.0, 400, 0.0, 1.0);
  new TH2F("KN_MM_KNpip_mom_woAll",  "d(K^{-}, n)\"X\" vs d(K^{-}, n #pi^{+}) mom",  400, 0.0, 2.0, 400, 0.0, 1.0);
  new TH2F("KN_MM_beam_mom_woAll",   "d(K^{-}, n)\"X\" vs beam mom",                 400, 0.0, 2.0, 400, 0.0, 2.0);
  new TH2F("KN_MM_n_mom_woAll",      "d(K^{-}, n)\"X\" vs n mom",                    400, 0.0, 2.0, 400, 0.0, 2.0);
  new TH2F("KN_MM_n_overbeta_woAll", "d(K^{-}, n)\"X\" vs n mom",                    400, 0.0, 2.0, 400, 0.0, 2.0);
  new TH2F("KN_MM_pim_mom_woAll",    "d(K^{-}, n)\"X\" vs #pi^{-} mom",              400, 0.0, 2.0, 400, 0.0, 1.0);
  new TH2F("KN_MM_pip_mom_woAll",    "d(K^{-}, n)\"X\" vs #pi^{+} mom",              400, 0.0, 2.0, 400, 0.0, 1.0);

  new TH1F("CDC_pim_chi2_Npipi",      "CDC #pi^{-} chi2", 1000, 0.0, 100.0);
  new TH1F("CDC_pip_chi2_Npipi",      "CDC #pi^{+} chi2", 1000, 0.0, 100.0);
  new TH1F("KNpipi_MM",               "d(K^{-}, n #pi^{-} #pi^{+})\"X\"",            3000, -1.0, 2.0);
  new TH1F("KNpipi_MM_wK0",           "d(K^{-}, n #pi^{-} #pi^{+})\"X\"",            3000, -1.0, 2.0);
  new TH1F("KNpipi_MM_wSp",           "d(K^{-}, n #pi^{-} #pi^{+})\"X\"",            3000, -1.0, 2.0);
  new TH1F("KNpipi_MM_wSp_rcSp",      "d(K^{-}, n #pi^{-} #pi^{+})\"X\"",            3000, -1.0, 2.0);
  new TH1F("KNpipi_MM_wSm",           "d(K^{-}, n #pi^{-} #pi^{+})\"X\"",            3000, -1.0, 2.0);
  new TH1F("KNpipi_MM_wSm_rcSm",      "d(K^{-}, n #pi^{-} #pi^{+})\"X\"",            3000, -1.0, 2.0);
  new TH1F("KNpipi_MM_woK0",          "d(K^{-}, n #pi^{-} #pi^{+})\"X\"",            3000, -1.0, 2.0);
  new TH1F("KNpipi_MM_woAll",         "d(K^{-}, n #pi^{-} #pi^{+})\"X\"",            3000, -1.0, 2.0);
  new TH1F("KNpipi_MM_woAll_wPi",     "d(K^{-}, n #pi^{-} #pi^{+})\"X\"",            3000, -1.0, 2.0);
  new TH1F("KNpipi_MM_woAll_BHD1hit", "d(K^{-}, n #pi^{-} #pi^{+})\"X\"",            3000, -1.0, 2.0);
  new TH1F("KNpipi_MM_woAll_chi2_50", "d(K^{-}, n #pi^{-} #pi^{+})\"X\"",            3000, -1.0, 2.0);
  new TH1F("KNpipi_MM_woAll_chi2_30", "d(K^{-}, n #pi^{-} #pi^{+})\"X\"",            3000, -1.0, 2.0);
  new TH1F("KNpipi_MM_woAll_DCA_3",   "d(K^{-}, n #pi^{-} #pi^{+})\"X\"",            3000, -1.0, 2.0); 
  new TH1F("KNpipi_MM_woAll_DCA_2",   "d(K^{-}, n #pi^{-} #pi^{+})\"X\"",            3000, -1.0, 2.0);
  new TH1F("KNpipi_MM_woAll_DCA_1",   "d(K^{-}, n #pi^{-} #pi^{+})\"X\"",            3000, -1.0, 2.0);
  new TH1F("KNpipi_MM_woAll2",        "d(K^{-}, n #pi^{-} #pi^{+})\"X\"",            3000, -1.0, 2.0);
  new TH1F("KNpipi_MM_woAll_woS2",    "d(K^{-}, n #pi^{-} #pi^{+})\"X\"",            3000, -1.0, 2.0);
  new TH1F("KNpipi_MM_woAll_rcS",     "d(K^{-}, n #pi^{-} #pi^{+})\"X\"",            3000, -1.0, 2.0);
  new TH1F("KNpipi_MM_woAll_rcS2",    "d(K^{-}, n #pi^{-} #pi^{+})\"X\"",            3000, -1.0, 2.0);
  new TH1F("KNpipi_MM_woAll_rcSm",    "d(K^{-}, n #pi^{-} #pi^{+})\"X\"",            3000, -1.0, 2.0);
  new TH1F("KNpipi_MM_woAll_rcSp",    "d(K^{-}, n #pi^{-} #pi^{+})\"X\"",            3000, -1.0, 2.0);

  new TH2F("KNpipi_MM_woAll_DCA",     "d(K^{-}, n #pi^{-} #pi^{+})\"X\" vs min DCA", 3000, -1.0, 2.0, 1000, 0.0, 1.0);

  //"p"(K-, pi+ pi-)n Reconstruction START
  new TH2F("KNpipi_MM_KP_pipi_MM", "KNpipi_MM_KP_pipi_MM", 2000, 0.0, 2.0, 2000, 0.0, 2.0);

  new TH1F("KN_MM_pipi_woAll_rcSm",     "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_woAll_rcSp",     "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_woAll_wN_rcSm",  "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_woAll_wN_rcSp",  "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH2F("KN_MM_woAll_MM_Npim_IM", "d(K^{-}, n)\"X\" vs \"n\" #pi^{-} IM", 600, 0.0, 3.0, 200, 1.0, 2.0);
  new TH2F("KN_MM_woAll_MM_Npip_IM", "d(K^{-}, n)\"X\" vs \"n\" #pi^{+} IM", 600, 0.0, 3.0, 200, 1.0, 2.0);

  new TH1F("MM_Npim_IM",             "n #pi^{-} IM", 1000, 1.0, 2.0);
  new TH1F("MM_Npip_IM",             "n #pi^{+} IM", 1000, 1.0, 2.0);
  new TH1F("MM_Npim_IM_wSm",         "n #pi^{-} IM", 1000, 1.0, 2.0);
  new TH1F("MM_Npip_IM_wSp",         "n #pi^{+} IM", 1000, 1.0, 2.0);
  new TH1F("MM_Npim_IM_no_cut",      "n #pi^{-} IM", 1000, 1.0, 2.0);
  new TH1F("MM_Npip_IM_no_cut",      "n #pi^{+} IM", 1000, 1.0, 2.0);
  new TH1F("MM_Npim_IM_woSm",        "n #pi^{-} IM", 1000, 1.0, 2.0);
  new TH1F("MM_Npip_IM_woSp",        "n #pi^{+} IM", 1000, 1.0, 2.0);

  new TH1F("Npim_IM_pip_MMSA", "n #pi^{-} IM", 1000, 1.0, 2.0);
  new TH1F("Npip_IM_pim_MMSA", "n #pi^{+} IM", 1000, 1.0, 2.0);

  new TH2F("KNpipi_MM_Npim_IM",           "d(K^{-}, n)\"X\" vs \"n\" #pi^{-} IM",                  3000, 0.0, 3.0, 1000, 1.0, 2.0);
  new TH2F("KNpipi_MM_Npip_IM",           "d(K^{-}, n)\"X\" vs \"n\" #pi^{+} IM",                  3000, 0.0, 3.0, 1000, 1.0, 2.0);
  new TH2F("MM_Npim_MM_Npip_MM",          "\"n\" #pi^{-} IM vs \"n\" #pi^{+} IM",                  1000, 1.0, 2.0, 1000, 1.0, 2.0);
  new TH2F("Npim_IM_MM_Npim_IM",          "n #pi^{-} vs \"n\" #pi^{-} IM",                         1000, 1.0, 2.0, 1000, 1.0, 2.0);
  new TH2F("Npip_IM_MM_Npip_IM",          "n #pi^{+} vs \"n\" #pi^{+} IM",                         1000, 1.0, 2.0, 1000, 1.0, 2.0);
  new TH2F("KNpipi_MM_MM_Npim_IM",        "d(K^{-}, n #pi^{+} #pi^{-})\"X\" vs \"n\" #pi^{-} IM",  2000, 0.0, 2.0, 1000, 1.0, 2.0);
  new TH2F("KNpipi_MM_MM_Npip_IM",        "d(K^{-}, n #pi^{+} #pi^{-})\"X\" vs \"n\" #pi^{+} IM",  2000, 0.0, 2.0, 1000, 1.0, 2.0);
  new TH2F("Npim_IM_MM_Npim_IM_no_cut",   "n #pi^{-} vs \"n\" #pi^{-} IM",                         1000, 1.0, 2.0, 1000, 1.0, 2.0);
  new TH2F("Npip_IM_MM_Npip_IM_no_cut",   "n #pi^{+} vs \"n\" #pi^{+} IM",                         1000, 1.0, 2.0, 1000, 1.0, 2.0);
  new TH2F("KNpipi_MM_MM_Npim_IM_no_cut", "d(K^{-}, n #pi^{+} #pi^{-})\"X\" vs \"n\" #pi^{-} IM",  2000, 0.0, 2.0, 1000, 1.0, 2.0);
  new TH2F("KNpipi_MM_MM_Npip_IM_no_cut", "d(K^{-}, n #pi^{+} #pi^{-})\"X\" vs \"n\" #pi^{+} IM",  2000, 0.0, 2.0, 1000, 1.0, 2.0);

  new TH1F("KNpipi_MM_rcSm",      "d(K^{-}, n #pi^{+} #pi^{-})\"X\"", 2000, 0.0, 2.0);
  new TH1F("KNpipi_MM_rcSp",      "d(K^{-}, n #pi^{+} #pi^{-})\"X\"", 2000, 0.0, 2.0);
  new TH1F("KNpipi_MM_rcSm_woSm", "d(K^{-}, n #pi^{+} #pi^{-})\"X\"", 2000, 0.0, 2.0);
  new TH1F("KNpipi_MM_rcSp_woSp", "d(K^{-}, n #pi^{+} #pi^{-})\"X\"", 2000, 0.0, 2.0);
  new TH1F("KNpipi_MM_rcS",       "d(K^{-}, n #pi^{+} #pi^{-})\"X\"", 2000, 0.0, 2.0);
  new TH1F("KNpipi_MM_wo_rcS",       "d(K^{-}, n #pi^{+} #pi^{-})\"X\"", 2000, 0.0, 2.0);
  new TH1F("KNpipi_MM_wo_rcS2",       "d(K^{-}, n #pi^{+} #pi^{-})\"X\"", 2000, 0.0, 2.0);
  new TH1F("KNpipi_MM_wo_rcS3",       "d(K^{-}, n #pi^{+} #pi^{-})\"X\"", 2000, 0.0, 2.0);
  new TH1F("KNpipi_MM_wo_rcS4",       "d(K^{-}, n #pi^{+} #pi^{-})\"X\"", 2000, 0.0, 2.0);
  new TH1F("KNpipi_MM_wo_rcS5",       "d(K^{-}, n #pi^{+} #pi^{-})\"X\"", 2000, 0.0, 2.0);
  new TH1F("KNpipi_MM_rcS_woS",   "d(K^{-}, n #pi^{+} #pi^{-})\"X\"", 2000, 0.0, 2.0);

  new TH1F("Npim_IM_rcSm",    "n #pi^{-} IM", 1000, 1.0, 2.0);
  new TH1F("Npim_IM_wN_rcSm", "n #pi^{-} IM", 1000, 1.0, 2.0);

  new TH1F("Npip_IM_rcSp",    "n #pi^{+} IM", 1000, 1.0, 2.0);
  new TH1F("Npip_IM_wN_rcSp", "n #pi^{+} IM", 1000, 1.0, 2.0);

  //"p"(K-, pi+ pi-)n Reconstruction FINISH

  new TH1F("KNpipi_MM_woAll_down", "d(K^{-}, n #pi^{+} #pi^{-})\"X\"", 3000, -1.0, 2.0);
  new TH1F("KNpipi_MM_woAll_up",   "d(K^{-}, n #pi^{+} #pi^{-})\"X\"", 3000, -1.0, 2.0);
  for( int mm=130; mm<161; mm++ ){ 
    new TH1F(Form("KNpipi_MM_woAll_%d", mm), "d(K^{-}, n #pi^{-} #pi^{+})\"X\"", 3000, -1.0, 2.0);
  }

  new TH2F("KNpipi_MM_missing_ene",       "d(K^{-}, n #pi^{-} #pi^{+})\"X\" vs E_{missing}", 3000, -1.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("KNpipi_MM_missing_ene_woAll", "d(K^{-}, n #pi^{-} #pi^{+})\"X\" vs E_{missing}", 3000, -1.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("KNpipi_MM_missing_ene_wK0",   "d(K^{-}, n #pi^{-} #pi^{+})\"X\" vs E_{missing}", 3000, -1.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("KNpipi_MM_missing_ene_wSp",   "d(K^{-}, n #pi^{-} #pi^{+})\"X\" vs E_{missing}", 3000, -1.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("KNpipi_MM_missing_ene_wSm",   "d(K^{-}, n #pi^{-} #pi^{+})\"X\" vs E_{missing}", 3000, -1.0, 2.0, 2000, -1.0, 1.0);

  new TH2F("KNpipi_MM_mom",       "d(K^{-}, n #pi^{-} #pi^{+})\"X\" vs mom", 3000, -1.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("KNpipi_MM_mom_wK0",   "d(K^{-}, n #pi^{-} #pi^{+})\"X\" vs mom", 3000, -1.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("KNpipi_MM_mom_woAll", "d(K^{-}, n #pi^{-} #pi^{+})\"X\" vs mom", 3000, -1.0, 2.0, 1000, 0.0, 1.0);

  new TH2F("KNpim_KNpip_MM",           "d(K^{-}, n #pi^{-}) MM vs d(K^{-}, n #pi^{+}) MM", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("KNpim_KNpip_MM_wN",        "d(K^{-}, n #pi^{-}) MM vs d(K^{-}, n #pi^{+}) MM", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("KNpim_KNpip_MM_woAll",     "d(K^{-}, n #pi^{-}) MM vs d(K^{-}, n #pi^{+}) MM", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("KNpim_KNpip_MM_woAll2",    "d(K^{-}, n #pi^{-}) MM vs d(K^{-}, n #pi^{+}) MM", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("KNpim_KNpip_MM_woAll_wN",  "d(K^{-}, n #pi^{-}) MM vs d(K^{-}, n #pi^{+}) MM", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("KNpim_KNpip_MM_woAll_wN0", "d(K^{-}, n #pi^{-}) MM vs d(K^{-}, n #pi^{+}) MM", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("KNpim_KNpip_MM_woAll_wN1", "d(K^{-}, n #pi^{-}) MM vs d(K^{-}, n #pi^{+}) MM", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("KNpim_KNpip_MM_woAll_wN2", "d(K^{-}, n #pi^{-}) MM vs d(K^{-}, n #pi^{+}) MM", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("KNpim_KNpip_MM_woAll2_wN",  "d(K^{-}, n #pi^{-}) MM vs d(K^{-}, n #pi^{+}) MM", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("KNpim_KNpip_MM_woAll2_wN0", "d(K^{-}, n #pi^{-}) MM vs d(K^{-}, n #pi^{+}) MM", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("KNpim_KNpip_MM_woAll2_wN1", "d(K^{-}, n #pi^{-}) MM vs d(K^{-}, n #pi^{+}) MM", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("KNpim_KNpip_MM_woAll2_wN2", "d(K^{-}, n #pi^{-}) MM vs d(K^{-}, n #pi^{+}) MM", 2000, 0.0, 2.0, 2000, 0.0, 2.0);

  new TH2F("KNpim_KNpip_MM_woAll_wN_up", "d(K^{-}, n #pi^{-}) MM vs d(K^{-}, n #pi^{+}) MM", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("KNpim_KNpip_MM_woAll_wN_down", "d(K^{-}, n #pi^{-}) MM vs d(K^{-}, n #pi^{+}) MM", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  for( int mm=1300; mm<1600; mm+=5 ){
    new TH2F(Form("KNpim_KNpip_MM_woAll_wN_%d", mm), "d(K^{-}, n #pi^{-}) MM vs d(K^{-}, n #pi^{+}) MM", 1000, 1.0, 2.0, 1000, 1.0, 2.0);
  }

  new TH2F("KNpim_KNpip_MM_woAll_wN_lt144", "d(K^{-}, n #pi^{-}) MM vs d(K^{-}, n #pi^{+}) MM", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("KNpim_KNpip_MM_woAll_wN_lt143", "d(K^{-}, n #pi^{-}) MM vs d(K^{-}, n #pi^{+}) MM", 2000, 0.0, 2.0, 2000, 0.0, 2.0);

  new TH1F("KNpim_ang_CM_npipi",           "#Sigma^{+} angle", 1000, -1.0, 1.0);
  new TH1F("KNpip_ang_CM_npipi",           "#Sigma^{-} angle", 1000, -1.0, 1.0);
  new TH1F("KNpim_ang_CM_npipi_mmSp",      "#Sigma^{+} angle", 1000, -1.0, 1.0);
  new TH1F("KNpip_ang_CM_npipi_mmSm",      "#Sigma^{+} angle", 1000, -1.0, 1.0);
  new TH1F("KNpim_ang_CM_npipi_mmSp_woSm", "#Sigma^{+} angle", 1000, -1.0, 1.0);
  new TH1F("KNpip_ang_CM_npipi_mmSm_woSp", "#Sigma^{+} angle", 1000, -1.0, 1.0);

  new TH1F("KNpim_ang_piS_npipi",           "#Sigma^{+} angle", 1000, -1.0, 1.0);
  new TH1F("KNpip_ang_piS_npipi",           "#Sigma^{-} angle", 1000, -1.0, 1.0);
  new TH1F("KNpim_ang_piS_npipi_mmSp",      "#Sigma^{+} angle", 1000, -1.0, 1.0);
  new TH1F("KNpip_ang_piS_npipi_mmSm",      "#Sigma^{+} angle", 1000, -1.0, 1.0);
  new TH1F("KNpim_ang_piS_npipi_mmSp_woSm", "#Sigma^{+} angle", 1000, -1.0, 1.0);
  new TH1F("KNpip_ang_piS_npipi_mmSm_woSp", "#Sigma^{+} angle", 1000, -1.0, 1.0);

  new TH1F("pim_ang_piS_npipi",           "#Sigma^{+} angle", 1000, -1.0, 1.0);
  new TH1F("pip_ang_piS_npipi",           "#Sigma^{-} angle", 1000, -1.0, 1.0);
  new TH1F("pim_ang_piS_npipi_mmSp",      "#Sigma^{+} angle", 1000, -1.0, 1.0);
  new TH1F("pip_ang_piS_npipi_mmSm",      "#Sigma^{+} angle", 1000, -1.0, 1.0);
  new TH1F("pim_ang_piS_npipi_mmSp_woSm", "#Sigma^{+} angle", 1000, -1.0, 1.0);
  new TH1F("pip_ang_piS_npipi_mmSm_woSp", "#Sigma^{+} angle", 1000, -1.0, 1.0);

  new TH1F("KNpim_ang_piS_npipi_mmSp_144",      "#Sigma^{+} angle", 1000, -1.0, 1.0);
  new TH1F("KNpip_ang_piS_npipi_mmSm_144",      "#Sigma^{+} angle", 1000, -1.0, 1.0);
  new TH1F("KNpim_ang_piS_npipi_mmSp_woSm_144", "#Sigma^{+} angle", 1000, -1.0, 1.0);
  new TH1F("KNpip_ang_piS_npipi_mmSm_woSp_144", "#Sigma^{+} angle", 1000, -1.0, 1.0);

  new TH1F("pim_ang_piS_npipi_mmSp_144",      "#Sigma^{+} angle", 1000, -1.0, 1.0);
  new TH1F("pip_ang_piS_npipi_mmSm_144",      "#Sigma^{+} angle", 1000, -1.0, 1.0);
  new TH1F("pim_ang_piS_npipi_mmSp_woSm_144", "#Sigma^{+} angle", 1000, -1.0, 1.0);
  new TH1F("pip_ang_piS_npipi_mmSm_woSp_144", "#Sigma^{+} angle", 1000, -1.0, 1.0);


  // for K0 Sf cut 2.5sigma
  new TH1F("KN_MM_pipi_woAll25",        "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_woAll25_wN",     "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_woAll25_wN_wS",  "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_woAll25_wN0",    "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_woAll25_wN0_wS", "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_woAll25_wN1",    "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_woAll25_wN1_wS", "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_woAll25_wN2",    "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_woAll25_wN2_wS", "d(K^{-}, n)", 3000, 0.0, 3.0);

  new TH1F("KNpipi_MM_woAll25",         "d(K^{-}, n #pi^{-} #pi^{+})\"X\"", 3000, 0.0, 3.0);
  new TH2F("KNpipi_MM_mom_woAll25", "d(K^{-}, n #pi^{-} #pi^{+})\"X\" vs mom", 3000, 0.0, 3.0, 1000, 0.0, 1.0);

  new TH2F("KNpim_KNpip_MM_woAll25",     "d(K^{-}, n #pi^{-}) MM vs d(K^{-}, n #pi^{+}) MM", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("KNpim_KNpip_MM_woAll25_wN",  "d(K^{-}, n #pi^{-}) MM vs d(K^{-}, n #pi^{+}) MM", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("KNpim_KNpip_MM_woAll25_wN0", "d(K^{-}, n #pi^{-}) MM vs d(K^{-}, n #pi^{+}) MM", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("KNpim_KNpip_MM_woAll25_wN1", "d(K^{-}, n #pi^{-}) MM vs d(K^{-}, n #pi^{+}) MM", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("KNpim_KNpip_MM_woAll25_wN2", "d(K^{-}, n #pi^{-}) MM vs d(K^{-}, n #pi^{+}) MM", 2000, 0.0, 2.0, 2000, 0.0, 2.0);

  // for K0 Sf cut 3sigma
  new TH1F("KN_MM_pipi_woAll3",        "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_woAll3_wN",     "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_woAll3_wN_wS",  "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_woAll3_wN0",    "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_woAll3_wN0_wS", "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_woAll3_wN1",    "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_woAll3_wN1_wS", "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_woAll3_wN2",    "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_woAll3_wN2_wS", "d(K^{-}, n)", 3000, 0.0, 3.0);

  new TH1F("KNpipi_MM_woAll3",         "d(K^{-}, n #pi^{-} #pi^{+})\"X\"", 3000, 0.0, 3.0);
  new TH2F("KNpipi_MM_mom_woAll3", "d(K^{-}, n #pi^{-} #pi^{+})\"X\" vs mom", 2000, 0.0, 2.0, 1000, 0.0, 1.0);

  new TH2F("KNpim_KNpip_MM_woAll3",     "d(K^{-}, n #pi^{-}) MM vs d(K^{-}, n #pi^{+}) MM", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("KNpim_KNpip_MM_woAll3_wN",  "d(K^{-}, n #pi^{-}) MM vs d(K^{-}, n #pi^{+}) MM", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("KNpim_KNpip_MM_woAll3_wN0", "d(K^{-}, n #pi^{-}) MM vs d(K^{-}, n #pi^{+}) MM", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("KNpim_KNpip_MM_woAll3_wN1", "d(K^{-}, n #pi^{-}) MM vs d(K^{-}, n #pi^{+}) MM", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("KNpim_KNpip_MM_woAll3_wN2", "d(K^{-}, n #pi^{-}) MM vs d(K^{-}, n #pi^{+}) MM", 2000, 0.0, 2.0, 2000, 0.0, 2.0);

  // for K0 Sf cut 4sigma
  new TH1F("KN_MM_pipi_woAll4",        "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_woAll4_wN",     "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_woAll4_wN_wS",  "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_woAll4_wN0",    "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_woAll4_wN0_wS", "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_woAll4_wN1",    "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_woAll4_wN1_wS", "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_woAll4_wN2",    "d(K^{-}, n)", 3000, 0.0, 3.0);
  new TH1F("KN_MM_pipi_woAll4_wN2_wS", "d(K^{-}, n)", 3000, 0.0, 3.0);

  new TH1F("KNpipi_MM_woAll4",         "d(K^{-}, n #pi^{-} #pi^{+})\"X\"", 2000, 0.0, 2.0);
  new TH2F("KNpipi_MM_mom_woAll4", "d(K^{-}, n #pi^{-} #pi^{+})\"X\" vs mom", 2000, 0.0, 2.0, 1000, 0.0, 1.0);

  new TH2F("KNpim_KNpip_MM_woAll4",     "d(K^{-}, n #pi^{-}) MM vs d(K^{-}, n #pi^{+}) MM", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("KNpim_KNpip_MM_woAll4_wN",  "d(K^{-}, n #pi^{-}) MM vs d(K^{-}, n #pi^{+}) MM", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("KNpim_KNpip_MM_woAll4_wN0", "d(K^{-}, n #pi^{-}) MM vs d(K^{-}, n #pi^{+}) MM", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("KNpim_KNpip_MM_woAll4_wN1", "d(K^{-}, n #pi^{-}) MM vs d(K^{-}, n #pi^{+}) MM", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("KNpim_KNpip_MM_woAll4_wN2", "d(K^{-}, n #pi^{-}) MM vs d(K^{-}, n #pi^{+}) MM", 2000, 0.0, 2.0, 2000, 0.0, 2.0);

  //*** for d(K-, n pi+ pi-) n tail test *****************************************************************************************//
  new TH2F("KNpipi_KN_MM",    "d(K^{-}, n #pi^{+} #pi^{-})\"X\" vs d(K^{-}, n)\"X\"",              2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("KNpipi_KNpim_MM", "d(K^{-}, n #pi^{+} #pi^{-})\"X\" vs d(K^{-}, n #pi^{-})\"X\"",      2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("KNpipi_KNpip_MM", "d(K^{-}, n #pi^{+} #pi^{-})\"X\" vs d(K^{-}, n #pi^{+})\"X\"",      2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("KNpipi_MM_beam_mom",    "d(K^{-}, n #pi^{+} #pi^{-})\"X\" vs beam mom",                2000, 0.0, 2.0, 1000, 0.5, 1.5);
  new TH2F("KNpipi_MM_nc_beta",     "d(K^{-}, n #pi^{+} #pi^{-})\"X\" vs n #beta",                 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("KNpipi_MM_nc_overbeta", "d(K^{-}, n #pi^{+} #pi^{-})\"X\" vs n 1/#beta",               2000, 0.0, 2.0, 5000, 0.0, 5.0);
  new TH2F("KNpipi_MM_pim_mom",     "d(K^{-}, n #pi^{-} #pi^{+})\"X\" vs #pi^{-} mom",             2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("KNpipi_MM_pip_mom",     "d(K^{-}, n #pi^{-} #pi^{+})\"X\" vs #pi^{+} mom",             2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("KNpipi_MM_KNpim_mom",   "d(K^{-}, n #pi^{+} #pi^{-})\"X\" vs d(K^{-}, n #pi^{-}) mom", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("KNpipi_MM_KNpip_mom",   "d(K^{-}, n #pi^{+} #pi^{-})\"X\" vs d(K^{-}, n #pi^{+}) mom", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("KNpipi_MM_n_mom",   "d(K^{-}, n #pi^{-} #pi^{+})\"X\" vs n mom",                       2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("KNpipi_MM_NC_seg",  "d(K^{-}, n #pi^{-} #pi^{+})\"X\" vs n mom",                       2000, 0.0, 2.0, 112,  0.5, 112.5);
  new TH2F("KNpipi_MM_NC_dE",  "d(K^{-}, n #pi^{-} #pi^{+})\"X\" vs NC dE",                        2000, 0.0, 2.0, 2000, 0.0, 200.0);

  new TH2F("KNpipi_KN_MM_wK0", "d(K^{-}, n #pi^{+} #pi^{-})\"X\" vs d(K^{-}, n)\"X\"",                   2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("KNpipi_KN_MM_wSm", "d(K^{-}, n #pi^{+} #pi^{-})\"X\" vs d(K^{-}, n)\"X\"",                   2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("KNpipi_KN_MM_wSp", "d(K^{-}, n #pi^{+} #pi^{-})\"X\" vs d(K^{-}, n)\"X\"",                   2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("KNpipi_KN_MM_woAll", "d(K^{-}, n #pi^{+} #pi^{-})\"X\" vs d(K^{-}, n)\"X\"",                 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("KNpipi_KNpim_MM_woAll", "d(K^{-}, n #pi^{+} #pi^{-})\"X\" vs d(K^{-}, n #pi^{-})\"X\"",      2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("KNpipi_KNpip_MM_woAll", "d(K^{-}, n #pi^{+} #pi^{-})\"X\" vs d(K^{-}, n #pi^{+})\"X\"",      2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("KNpipi_MM_beam_mom_woAll",    "d(K^{-}, n #pi^{+} #pi^{-})\"X\" vs beam mom",                2000, 0.0, 2.0, 1000, 0.5, 1.5);
  new TH2F("KNpipi_MM_nc_beta_woAll",     "d(K^{-}, n #pi^{+} #pi^{-})\"X\" vs n #beta",                 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("KNpipi_MM_nc_overbeta_woAll", "d(K^{-}, n #pi^{+} #pi^{-})\"X\" vs n 1/#beta",               2000, 0.0, 2.0, 5000, 0.0, 5.0);
  new TH2F("KNpipi_MM_pim_mom_woAll", "d(K^{-}, n #pi^{-} #pi^{+})\"X\" vs #pi^{-} mom",                 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("KNpipi_MM_pip_mom_woAll", "d(K^{-}, n #pi^{-} #pi^{+})\"X\" vs #pi^{+} mom",                 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("KNpipi_MM_KNpim_mom_woAll",   "d(K^{-}, n #pi^{+} #pi^{-})\"X\" vs d(K^{-}, n #pi^{-}) mom", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("KNpipi_MM_KNpip_mom_woAll",   "d(K^{-}, n #pi^{+} #pi^{-})\"X\" vs d(K^{-}, n #pi^{+}) mom", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("KNpipi_MM_n_mom_woAll",   "d(K^{-}, n #pi^{-} #pi^{+})\"X\" vs n mom",                       2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("KNpipi_MM_NC_seg_woAll",  "d(K^{-}, n #pi^{-} #pi^{+})\"X\" vs n mom",                       2000, 0.0, 2.0, 112,  0.5, 112.5);
  new TH2F("KNpipi_MM_NC_dE_woAll",  "d(K^{-}, n #pi^{-} #pi^{+})\"X\" vs NC dE",                        2000, 0.0, 2.0, 2000, 0.0, 200.0);

  new TH1F("CDS_IM_pipi_test", "CDS #pi^{+} #pi^{-} IM", 1000, 0.0, 1.0);
  new TH1F("Npim_IM_test", "n #pi^{-} IM", 1000, 1.0, 2.0);
  new TH1F("Npip_IM_test", "n #pi^{+} IM", 1000, 1.0, 2.0);
  new TH2F("KNpim_KNpip_MM_test",     "d(K^{-}, n #pi^{-}) MM vs d(K^{-}, n #pi^{+}) MM", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("CDS_mass2_mom_pim_test", "CDS mass2 vs mom", 2000, -0.5, 9.5, 2000, -1.0, 1.0);
  new TH2F("CDS_mass2_mom_pip_test", "CDS mass2 vs mom", 2000, -0.5, 9.5, 2000, -1.0, 1.0);
  new TH1F("pim_dmom_test", "#pi^{-} dmom", 1000, -0.1, 0.4);
  new TH1F("pip_dmom_test", "#pi^{+} dmom", 1000, -0.1, 0.4);

  //*** n pi+ pi- IM Analysis ****************************************************************************************************//
  new TH2F("KN_MM_Npipi_IM",       "d(K^{-}, n) MM vs n #pi^{+} #pi^{-} IM", 3000, 0.0, 3.0, 1000, 1.0, 2.0);
  new TH2F("KN_MM_Npipi_IM_wN",    "d(K^{-}, n) MM vs n #pi^{+} #pi^{-} IM", 3000, 0.0, 3.0, 1000, 1.0, 2.0);
  new TH2F("KN_MM_Npipi_IM_wN_K0", "d(K^{-}, n) MM vs n #pi^{+} #pi^{-} IM", 3000, 0.0, 3.0, 1000, 1.0, 2.0);
  new TH2F("KN_MM_Npipi_IM_wN_Sp", "d(K^{-}, n) MM vs n #pi^{+} #pi^{-} IM", 3000, 0.0, 3.0, 1000, 1.0, 2.0);
  new TH2F("KN_MM_Npipi_IM_wN_Sm", "d(K^{-}, n) MM vs n #pi^{+} #pi^{-} IM", 3000, 0.0, 3.0, 1000, 1.0, 2.0);

  new TH1F("Npipi_IM", "n #pi^{+} #pi^{-} IM", 1000, 1.0, 2.0);
  new TH1F("Npipi_IM_wK0", "n #pi^{+} #pi^{-} IM", 1000, 1.0, 2.0);
  new TH1F("Npipi_IM_wSf", "n #pi^{+} #pi^{-} IM", 1000, 1.0, 2.0);
  new TH1F("Npipi_IM_wSp", "n #pi^{+} #pi^{-} IM", 1000, 1.0, 2.0);
  new TH1F("Npipi_IM_wSm", "n #pi^{+} #pi^{-} IM", 1000, 1.0, 2.0);
  new TH1F("Npipi_IM_wSp_woSm", "n #pi^{+} #pi^{-} IM", 1000, 1.0, 2.0);
  new TH1F("Npipi_IM_wSm_woSp", "n #pi^{+} #pi^{-} IM", 1000, 1.0, 2.0);

  new TH1F("Npipi_IM_wN", "n #pi^{+} #pi^{-} IM", 1000, 1.0, 2.0);
  new TH1F("Npipi_IM_wN_wK0", "n #pi^{+} #pi^{-} IM", 1000, 1.0, 2.0);
  new TH1F("Npipi_IM_wN_wK0_SB0", "n #pi^{+} #pi^{-} IM", 1000, 1.0, 2.0);
  new TH1F("Npipi_IM_wN_wK0_SB1", "n #pi^{+} #pi^{-} IM", 1000, 1.0, 2.0);
  new TH1F("Npipi_IM_wN_wK0_SB2", "n #pi^{+} #pi^{-} IM", 1000, 1.0, 2.0);
  new TH1F("Npipi_IM_wN_wK0_SB3", "n #pi^{+} #pi^{-} IM", 1000, 1.0, 2.0);
  new TH1F("Npipi_IM_wN_wK0_SB4", "n #pi^{+} #pi^{-} IM", 1000, 1.0, 2.0);
  new TH1F("Npipi_IM_wN_wK0_SB5", "n #pi^{+} #pi^{-} IM", 1000, 1.0, 2.0);
  new TH1F("Npipi_IM_wN_wSf", "n #pi^{+} #pi^{-} IM", 1000, 1.0, 2.0);
  new TH1F("Npipi_IM_wN_wSp", "n #pi^{+} #pi^{-} IM", 1000, 1.0, 2.0);
  new TH1F("Npipi_IM_wN_wSm", "n #pi^{+} #pi^{-} IM", 1000, 1.0, 2.0);
  new TH1F("Npipi_IM_wN_wSp_woSm", "n #pi^{+} #pi^{-} IM", 1000, 1.0, 2.0);
  new TH1F("Npipi_IM_wN_wSm_woSp", "n #pi^{+} #pi^{-} IM", 1000, 1.0, 2.0);

  new TH2F("Npipi_IM_KNpipi_MM", "n #pi^{+} #pi^{-} IM vs d(K^{-}, n #pi^{-} #pi^{+})\"X\"",          1000, 1.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("Npipi_IM_KNpipi_MM_wK0", "n #pi^{+} #pi^{-} IM vs d(K^{-}, n #pi^{-} #pi^{+})\"X\"",      1000, 1.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("Npipi_IM_KNpipi_MM_wSf", "n #pi^{+} #pi^{-} IM vs d(K^{-}, n #pi^{-} #pi^{+})\"X\"",      1000, 1.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("Npipi_IM_KNpipi_MM_wSp", "n #pi^{+} #pi^{-} IM vs d(K^{-}, n #pi^{-} #pi^{+})\"X\"",      1000, 1.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("Npipi_IM_KNpipi_MM_wSm", "n #pi^{+} #pi^{-} IM vs d(K^{-}, n #pi^{-} #pi^{+})\"X\"",      1000, 1.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("Npipi_IM_KNpipi_MM_wSp_woSm", "n #pi^{+} #pi^{-} IM vs d(K^{-}, n #pi^{-} #pi^{+})\"X\"", 1000, 1.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("Npipi_IM_KNpipi_MM_wSm_woSp", "n #pi^{+} #pi^{-} IM vs d(K^{-}, n #pi^{-} #pi^{+})\"X\"", 1000, 1.0, 2.0, 2000, 0.0, 2.0);

  new TH2F("Npipi_IM_ME", "n #pi^{+} #pi^{-} IM vs d(K^{-}, n #pi^{-} #pi^{+})\"X\"",          1000, 1.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("Npipi_IM_ME_wK0", "n #pi^{+} #pi^{-} IM vs d(K^{-}, n #pi^{-} #pi^{+})\"X\"",      1000, 1.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("Npipi_IM_ME_wSf", "n #pi^{+} #pi^{-} IM vs d(K^{-}, n #pi^{-} #pi^{+})\"X\"",      1000, 1.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("Npipi_IM_ME_wSp", "n #pi^{+} #pi^{-} IM vs d(K^{-}, n #pi^{-} #pi^{+})\"X\"",      1000, 1.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("Npipi_IM_ME_wSm", "n #pi^{+} #pi^{-} IM vs d(K^{-}, n #pi^{-} #pi^{+})\"X\"",      1000, 1.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("Npipi_IM_ME_wSp_woSm", "n #pi^{+} #pi^{-} IM vs d(K^{-}, n #pi^{-} #pi^{+})\"X\"", 1000, 1.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("Npipi_IM_ME_wSm_woSp", "n #pi^{+} #pi^{-} IM vs d(K^{-}, n #pi^{-} #pi^{+})\"X\"", 1000, 1.0, 2.0, 2000, -1.0, 1.0);

  new TH1F("KP_Npipi_MM", "\"p\"(K^{-} n #pi^{+} #pi^{-})\"X\"", 2000, -1.0, 1.0);
  new TH2F("Npipi_IM_KP_Npipi_MM", "n #pi^{+} #pi^{-} IM vs \"p\"(K^{-}, n #pi^{-} #pi^{+})\"X\"", 1000, 1.0, 2.0, 2000, -1.0, 1.0);

  //**********************************//
  //*** d(K-, n)p pi- tagged event ***//
  //**********************************//
  new TH1F("CDS_IM_ppim_wN", "CDS p #pi^{-} IM", 1000, 1.0, 2.0);
  new TH1F("KN_MM_ppim",  "d(K^{-}, n)\"X\"", 3000, 0.0, 3.0);
  new TH1F("KN_MM_L",  "d(K^{-}, n)\"X\"", 3000, 0.0, 3.0);
  new TH1F("KNppim_MM",  "d(K^{-}, n)\"X\"", 3000, 0.0, 3.0);
  new TH1F("KNL_MM",     "d(K^{-}, n)\"X\"", 3000, 0.0, 3.0);
  new TH2F("KN_MM_KNL_MM",  "d(K^{-}, n)\"n\" vs d(K^{-}, n)\"X\"", 3000, 0.0, 3.0, 3000, 0.0, 3.0);

  //*** for K- d-> L(1405) n Sim *************************************************************************************************//
  new TH2F("L1405_KN_MM_res", "d(K^{-}, n) diff", 1000, 1.0, 2.0, 1000, -100, 100);
  new TH1F("KN_MM_res", "d(K^{-}, n) diff", 1000, -100, 100);

  new TH1F("L1405_eff",         "#Lambda(1405) effective event", 2000, 0.0, 2.0);
  new TH1F("L1405_pipi_hit",    "#Lambda(1405) effective event", 2000, 0.0, 2.0);
  new TH1F("L1405_pipi_hit_K0", "#Lambda(1405) effective event", 2000, 0.0, 2.0);
  new TH1F("L1405_pipi_hit_Sf", "#Lambda(1405) effective event", 2000, 0.0, 2.0);
  new TH1F("L1405_pipi_hit_woAll", "#Lambda(1405) effective event", 2000, 0.0, 2.0);
  new TH1F("L1405_pipi_hit_woAll_wN", "#Lambda(1405) effective event", 2000, 0.0, 2.0);
  new TH1F("L1405_pipi_hit_woAll_wN_wS", "#Lambda(1405) effective event", 2000, 0.0, 2.0);

  new TH1F("L1405_eff_Sm",                  "#Lambda(1405) effective event", 2000, 0.0, 2.0);
  new TH1F("L1405_pipi_hit_Sm",             "#Lambda(1405) effective event", 2000, 0.0, 2.0);
  new TH1F("L1405_pipi_hit_K0_Sm",          "#Lambda(1405) effective event", 2000, 0.0, 2.0);
  new TH1F("L1405_pipi_hit_Sf_Sm",          "#Lambda(1405) effective event", 2000, 0.0, 2.0);
  new TH1F("L1405_pipi_hit_woAll_Sm",       "#Lambda(1405) effective event", 2000, 0.0, 2.0);
  new TH1F("L1405_pipi_hit_woAll_wN_Sm",    "#Lambda(1405) effective event", 2000, 0.0, 2.0);
  new TH1F("L1405_pipi_hit_woAll_wN_wS_Sm", "#Lambda(1405) effective event", 2000, 0.0, 2.0);
  new TH1F("L1405_pipi_hit_woAll_wN_wSm_Sm", "#Lambda(1405) effective event", 2000, 0.0, 2.0);
  new TH1F("L1405_pipi_hit_woAll_wN_wSp_Sm", "#Lambda(1405) effective event", 2000, 0.0, 2.0);
  new TH1F("L1405_pipi_hit_woAll_wN_wSm_woSp_Sm", "#Lambda(1405) effective event", 2000, 0.0, 2.0);
  new TH1F("L1405_pipi_hit_woAll_wN_wSp_woSm_Sm", "#Lambda(1405) effective event", 2000, 0.0, 2.0);

  new TH1F("L1405_eff_Sp",                  "#Lambda(1405) effective event", 2000, 0.0, 2.0);
  new TH1F("L1405_pipi_hit_Sp",             "#Lambda(1405) effective event", 2000, 0.0, 2.0);
  new TH1F("L1405_pipi_hit_K0_Sp",          "#Lambda(1405) effective event", 2000, 0.0, 2.0);
  new TH1F("L1405_pipi_hit_Sf_Sp",          "#Lambda(1405) effective event", 2000, 0.0, 2.0);
  new TH1F("L1405_pipi_hit_woAll_Sp",       "#Lambda(1405) effective event", 2000, 0.0, 2.0);
  new TH1F("L1405_pipi_hit_woAll_wN_Sp",    "#Lambda(1405) effective event", 2000, 0.0, 2.0);
  new TH1F("L1405_pipi_hit_woAll_wN_wS_Sp", "#Lambda(1405) effective event", 2000, 0.0, 2.0);
  new TH1F("L1405_pipi_hit_woAll_wN_wSm_Sp", "#Lambda(1405) effective event", 2000, 0.0, 2.0);
  new TH1F("L1405_pipi_hit_woAll_wN_wSp_Sp", "#Lambda(1405) effective event", 2000, 0.0, 2.0);
  new TH1F("L1405_pipi_hit_woAll_wN_wSm_woSp_Sp", "#Lambda(1405) effective event", 2000, 0.0, 2.0);
  new TH1F("L1405_pipi_hit_woAll_wN_wSp_woSm_Sp", "#Lambda(1405) effective event", 2000, 0.0, 2.0);

  new TH1F("n_mom_res",       "n mom diff",          1000, -100, 100);

  //*** for TTree ***//
  // fNpipiData = new NpipiData();
  // new TTree("NpipiTree", "NpipiTree");
  // TTree *tree = (TTree*)rtFile->Get("NpipiTree");
  // tree-> Branch("NpipiData", &fNpipiData);
  
  if( simReader ) simReader->initHist(rtFile);

  std::cout<<"===== HistManwMC::initHist FINISH ====="<<std::endl;
}
