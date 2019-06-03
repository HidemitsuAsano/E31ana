
  //= = = = pipipnn final-sample tree = = = =//
  TLorentzVector *LVec_beam=nullptr;   // 4-momentum(beam)
  TLorentzVector *LVec_beam_Sp=nullptr;   // 4-momentum(beam),Sp mode assumption
  TLorentzVector *LVec_beam_Sm=nullptr;   // 4-momentum(beam),Sm mode assumption
  TLorentzVector *LVec_target=nullptr; // 4-momentum(target)
  TLorentzVector *LVec_pip=nullptr;    // 4-momentum(pi+)
  TLorentzVector *LVec_pim=nullptr;    // 4-momentum(pi-)
  TLorentzVector *LVec_n=nullptr;      // 4-momentum(neutron)
  TLorentzVector *LVec_n_Sp=nullptr;      // 4-momentum(neutron),Sp mode assumption
  TLorentzVector *LVec_n_Sm=nullptr;      // 4-momentum(neutron),Sm mode assumption
  TLorentzVector *mcmom_beam=nullptr;   // generated 4-momentum(beam)
  TLorentzVector *mcmom_pip=nullptr;    // generated 4-momentum(pi+)
  TLorentzVector *mcmom_pim=nullptr;    // generated 4-momentum(pi-)
  TLorentzVector *mcmom_ncds=nullptr;      // generated 4-momentum(neutron)
  TLorentzVector *mcmom_nmiss=nullptr;      // generated 4-momentum(neutron)
  double NeutralBetaCDH; // velocity of neutral particle on CDH
  double NeutralBetaCDH_vtx[2]; // velocity of neutral particle on CDH,0: Spmode 1:Smmode
  double dE;   // energy deposit on CDH
  TVector3 *vtx_reaction = nullptr; // vertex(reaction) 
  TVector3 *vtx_pip_beam = nullptr; //C.A.P of pip-beam beam side
  TVector3 *vtx_pim_beam = nullptr; //C.A.P of pim-beam beam side
  TVector3 *vtx_pip_cdc = nullptr;//C.A.P of pip-beam pip side
  TVector3 *vtx_pim_cdc = nullptr;//C.A.P of pim-beam pim side
  TVector3 *CA_pip = nullptr;//C.A.P of pip-pim pip side
  TVector3 *CA_pim = nullptr;//C.A.P of pip-pim pim side
  //int run_num;   // run number
  //int event_num; // event number
  //int block_num; // block number
  TLorentzVector *kfSpmode_mom_beam=nullptr;   // 4-momentum(beam) after kinematical refit for pi- Sigma+
  TLorentzVector *kfSpmode_mom_pip=nullptr;    // 4-momentum(pi+) after kinematical refit for pi- Sigma+
  TLorentzVector *kfSpmode_mom_pim=nullptr;    // 4-momentum(pi-) after kinematical refit for pi- Sigma+
  TLorentzVector *kfSpmode_mom_n=nullptr;      // 4-momentum(neutron) after kinematical refit for pi- Sigma+
  double kfSpmode_chi2;   // chi2 of kinematical refit
  double kfSpmode_NDF;    // NDF of kinematical refit
  double kfSpmode_status; // status of kinematical refit -> details can be found in this code
  double kfSpmode_pvalue; // p-value of kinematical refit
  TLorentzVector *kfSmmode_mom_beam=nullptr;   // 4-momentum(beam) after kinematical refit for pi+ Sigma-
  TLorentzVector *kfSmmode_mom_pip=nullptr;    // 4-momentum(pi+) after kinematical refit for pi+ Sigma-
  TLorentzVector *kfSmmode_mom_pim=nullptr;    // 4-momentum(pi-) after kinematical refit for pi+ Sigma-
  TLorentzVector *kfSmmode_mom_n=nullptr;      // 4-momentum(neutron) after kinematical refit for pi+ Sigma-
  double kfSmmode_chi2;   // chi2 of kinematical refit
  double kfSmmode_NDF;    // NDF of kinematical refit
  double kfSmmode_status; // status of kinematical refit -> details can be found in this code
  double kfSmmode_pvalue; // p-value of kinematical refit
  int kf_flag; // flag of correct pair reconstruction, etc

