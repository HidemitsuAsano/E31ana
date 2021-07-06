
  //= = = = pipipnn final-sample tree = = = =//
  TLorentzVector *LVec_beam=NULL;   // 4-momentum(beam)
  TLorentzVector *LVec_beam_Sp=NULL;   // 4-momentum(beam),Sp mode assumption
  TLorentzVector *LVec_beam_Sm=NULL;   // 4-momentum(beam),Sm mode assumption
  TLorentzVector *LVec_beam_K0=NULL;   // 4-momentum(beam),K0 mode assumption
  TLorentzVector *LVec_target=NULL; // 4-momentum(target)
  TLorentzVector *LVec_pip=NULL;    // 4-momentum(pi+)
  TLorentzVector *LVec_pim=NULL;    // 4-momentum(pi-)
  TLorentzVector *LVec_n=NULL;      // 4-momentum(neutron)
  TLorentzVector *LVec_n_beam=NULL;      // 4-momentum(neutron)
  TLorentzVector *LVec_n_Sp=NULL;      // 4-momentum(neutron),Sp mode assumption
  TLorentzVector *LVec_n_Sm=NULL;      // 4-momentum(neutron),Sm mode assumption
  TLorentzVector *LVec_n_K0=NULL;      // 4-momentum(neutron),K0 mode assumption
  TLorentzVector *mcmom_beam=NULL;   // generated 4-momentum(beam)
  TLorentzVector *mcmom_pip=NULL;    // generated 4-momentum(pi+)
  TLorentzVector *mcmom_pim=NULL;    // generated 4-momentum(pi-)
  TLorentzVector *mcmom_ncds=NULL;      // generated 4-momentum(neutron)
  TLorentzVector *mcmom_nmiss=NULL;      // generated 4-momentum(neutron)
  TLorentzVector *react_nmiss=NULL;      // generated 4-momentum(neutron)
  TLorentzVector *react_Sigma=NULL;      // generated 4-momentum(neutron)
  TLorentzVector *react_pi=NULL;      // generated 4-momentum(neutron)
  double NeutralBetaCDH; // velocity of neutral particle on CDH
  double NeutralBetaCDH_beam; // velocity of neutral particle on CDH
  double NeutralBetaCDH_vtx[2]; // velocity of neutral particle on CDH,0: Spmode 1:Smmode
  double tofpim;
  double tofpip;
  double tofn;
  double dE;   // energy deposit on CDH
  int neutralseg;
  int nhitOutCDC;
  int ForwardCharge;
  double mcncanvtxr;
  double mcncanvtxz;
  int mcncdsgen;
  int mcpattern;
  TVector3 *vtx_reaction = NULL; // vertex(reaction) 
  TVector3 *vtx_displaced = NULL; // vertex(displaced) 
  TVector3 *vtx_pip_beam = NULL; //C.A.P of pip-beam beam side
  TVector3 *vtx_pim_beam = NULL; //C.A.P of pim-beam beam side
  TVector3 *vtx_pip_cdc = NULL;//C.A.P of pip-beam pip side
  TVector3 *vtx_pim_cdc = NULL;//C.A.P of pim-beam pim side
  TVector3 *CA_pip = NULL;//C.A.P of pip-pim pip side
  TVector3 *CA_pim = NULL;//C.A.P of pip-pim pim side
  TVector3 *CDH_Pos = NULL;
  TVector3 *CDH_Pos_pim = NULL;
  TVector3 *CDH_Pos_pip = NULL;
  TVector3 *mc_vtx = NULL;
  TVector3 *mc_disvtx = NULL;
  //int run_num;   // run number
  //int event_num; // event number
  //int block_num; // block number
  TLorentzVector *kfSpmode_mom_beam=NULL;   // 4-momentum(beam) after kinematical refit for pi- Sigma+
  TLorentzVector *kfSpmode_mom_pip=NULL;    // 4-momentum(pi+) after kinematical refit for pi- Sigma+
  TLorentzVector *kfSpmode_mom_pim=NULL;    // 4-momentum(pi-) after kinematical refit for pi- Sigma+
  TLorentzVector *kfSpmode_mom_n=NULL;      // 4-momentum(neutron) after kinematical refit for pi- Sigma+
  double kfSpmode_chi2;   // chi2 of kinematical refit
  double kfSpmode_NDF;    // NDF of kinematical refit
  double kfSpmode_status; // status of kinematical refit -> details can be found in this code
  double kfSpmode_pvalue; // p-value of kinematical refit
  TLorentzVector *kfSmmode_mom_beam=NULL;   // 4-momentum(beam) after kinematical refit for pi+ Sigma-
  TLorentzVector *kfSmmode_mom_pip=NULL;    // 4-momentum(pi+) after kinematical refit for pi+ Sigma-
  TLorentzVector *kfSmmode_mom_pim=NULL;    // 4-momentum(pi-) after kinematical refit for pi+ Sigma-
  TLorentzVector *kfSmmode_mom_n=NULL;      // 4-momentum(neutron) after kinematical refit for pi+ Sigma-
  double kfSmmode_chi2;   // chi2 of kinematical refit
  double kfSmmode_NDF;    // NDF of kinematical refit
  double kfSmmode_status; // status of kinematical refit -> details can be found in this code
  double kfSmmode_pvalue; // p-value of kinematical refit
  int kf_flag; // flag of correct pair reconstruction, etc


  

