// MyAnalysisHeKpippfwd.cpp

#include "MyAnalysisHeKpippfwd.h"

MyAnalysisHeKpippfwd::MyAnalysisHeKpippfwd(TFile* rt, ConfMan* conf)
{
  Initialize(conf);
  CutCondition();
  Clear();
}

MyAnalysisHeKpippfwd::~MyAnalysisHeKpippfwd()
{
  Clear();
  rtFile->cd();
  rtFile->Write();
  rtFile->Close();
}

void MyAnalysisHeKpippfwd::Clear()
{
}

bool MyAnalysisHeKpippfwd::DoAnalysis(ConfMan* conf, EventHeader* header, BeamLineHitMan* blMan, BeamLineTrackMan* bltrackMan, CDSHitMan* cdsMan, CDSTrackingMan* cdstrackMan, Particle* particle)
{
  rtFile->cd();

  FillHist("EventNumber",0); /* All Events */
  DetectorList *dlist=DetectorList::GetInstance();
  if(particle->nBeam()!=1) return false;
  pBeam* beam = particle->beam(0);
  double beamtof = beam->bhdt0tof();

  if(cdstrackMan->nGoodTrack()!=2) return true;
  FillHist("EventNumber",1); /* 2 Good Tracks in CDC */

  /* Event selection */
  int MulCDH=0;
  for(int i=0; i<cdsMan->nCDH(); i++){
    HodoscopeLikeHit* hit = cdsMan->CDH(i);
    if(hit->CheckRange()) MulCDH++;
  }
  if(MulCDH!=2) return true;

  if(particle->nCDS()!=2) return true;
  FillHist("EventNumber",2); /* CDH 2 hits and No FWD hit */

  if(particle->nProton()!=1) return true;
  if(particle->nPiminus()!=1) return true;
  FillHist("EventNumber",3); /* proton and pion in CDS */

  if(particle->nPC()!=1) return true;
  FillHist("EventNumber",4); /* proton in FWD */

  pCDS* pim = particle->pim(0);
  pCDS* p1 = particle->proton(0);
  pPC*  p2 = particle->pc(0);

  // ############## //
  // CDC Chi-square //
  // ############## //
  if(pim->chi()>30) return true;
  if(p1->chi()>30) return true;
  FillHist("EventNumber",5); /* CDC Chi-square < 30 */

  // ######################### //
  // 2 Helix calculation       //
  // ######################### //
  pCDS* pip1 = 0;
  {
    for(int it=0; it<particle->nProduct(); it++){
      pCDS* product = particle->product(it);
      int comb = product->comb();
      if(comb==pow(2,CDS_Proton)+pow(2,CDS_PiMinus)){
        if(product->daughter1()==p1->id()||product->daughter2()==p1->id()){
          pip1 = product;
        }
      }
    }
  }
  if(pip1==0){ return true; }
  FillHist("EventNumber",5); /* 2-track Combination */

  bool pip1flag = false;
  pCDS* lam = 0;
  if(lamll<pip1->mass()&&pip1->mass()<lamul){
    pip1flag = true;
    lam = pip1;
  }
  if(!pip1flag){ return true; }
  if(lam==0){ return true; }
  FillHist("EventNumber",6); /* Lambda Reconstruction */

  // ################ //
  // Vertex decision1 //
  // ################ //
  TVector3 vertex_pim;
  vertex = pim->vbeam();

  // ################ //
  // Vertex decision2 //
  // ################ //
  TVector3 vertex_lam = lam->vbeam();
  if(GeomTools::GetID(vertex_lam)!=CID_Fiducial){ return true; }
  FillHist("EventNumber",7); /* Vertex Decision 2 */

  TVector3 vertex = vertex_lam

  // ##################### //
  // FWD Charged  Analysis //
  // ##################### //
  for(int ipc=0; ipc<particle->nPC(); ipc++){
    pPC* pc = particle->pc(ipc);
    int pid = -1;
    pc->CalcMom(beam,vertex,pid,false);
  }

  // ######################### //
  // Checking for CDS particle //
  // ######################### //
  int ncds[6] = {0,0,0,0,0,0};
  ncds[0] = particle->nPiplus();
  ncds[1] = particle->nPiminus();
  ncds[2] = particle->nKaon();
  ncds[3] = particle->nProton();
  ncds[4] = particle->nDeuteron();
  ncds[5] = particle->nTriton();
  ncds[5] += particle->nHelium3();
  ncds[5] += particle->nOther();
  TString cdsstr[6] = {"#pi^{+}","#pi^{-}","K^{-}","p","d","Other"};
  FillHist("CDS_NumOfParticle",particle->nCDS());
  for(int i=0; i<6; i++){
    FillHist("CDS_Particle",cdsstr[i],Form("%d track",particle->nCDS()),ncds[i]);
  }
  TVector3 vtxb;
  /* Pi Minus */
  vtxb = pim->vbeam();
  FillHist("pim_Vertex_XY",vtxb.X(),vtxb.Y());
  FillHist("pim_Vertex_ZX",vtxb.Z(),vtxb.X());
  FillHist("pim_Vertex_ZY",vtxb.Z(),vtxb.Y());
  FillHist("pim_VBDIS",pim->vbdis());
  FillHist("pim_OverbetavsMomentum",1.0/pim->beta(),pim->mom());
  FillHist("pim_Mass2vsMomentum",pim->mass2(),pim->mom());
  /* Proton1 */
  vtxb = p1->vbeam();
  FillHist("p1_Vertex_XY",vtxb.X(),vtxb.Y());
  FillHist("p1_Vertex_ZX",vtxb.Z(),vtxb.X());
  FillHist("p1_Vertex_ZY",vtxb.Z(),vtxb.Y());
  FillHist("p1_VBDIS",p1->vbdis());
  FillHist("p1_OverbetavsMomentum",1.0/p1->beta(),p1->mom());
  FillHist("p1_Mass2vsMomentum",p1->mass2(),p1->mom());
  /* Pi+ and Proton1 */
  vtxb = pip1->vbeam();
  FillHist("pip1_Vertex_XY",vtxb.X(),vtxb.Y());
  FillHist("pip1_Vertex_ZX",vtxb.Z(),vtxb.X());
  FillHist("pip1_Vertex_ZY",vtxb.Z(),vtxb.Y());
  FillHist("pip1_Mass",pip1->vdis());
  FillHist("pip1_OA",pip1->oa());
  FillHist("pip1_VDIS",pip1->vdis());
  FillHist("pip1_VBDIS",pip1->vbdis());
  FillHist("pip1_PBDCA",pip1->pbdca());

  // ############ //
  // DCA Decision //
  // ############ //
  if(lam->pbdca()>1.0){ return true; }
  FillHist("EventNumber",8); /* DCA Decision */

  // Trigger Pattern
  FillHist("TriggerPattern",0);
  for( int i=0; i<20; i++ ){
    int val = header->pattern(i);
    if( 0<val ){
      FillHist("TriggerPattern",i);
    }
  }

  TString frame[2] = {"Lab","CM"};
  TString pname[8] = {"pip","p","d","t","he","pim","k","o"};

  // Target
  TVector3 ZeroV;
  TLorentzVector Ltgt; Ltgt.SetVectM(ZeroV,ThreeHeMass);
  TLorentzVector Lbeam = beam->GetLorentzVector(vertex);
  TVector3 boost = (Ltgt+Lbeam).BoostVector();
  TLorentzVector cLtgt = Ltgt; cLtgt.Boost(-boost); TLorentzVector cLbeam = Lbeam; cLbeam.Boost(-boost);

  // Pi+, Pi-, and Pi+ Pi- pair //
  TLorentzVector Lpim = pim->GetLorentzVector();
  TLorentzVector cLpim = pim->GetLorentzVector(); cLpim.Boost(-boost);
  TLorentzVector Lp1 = p1->GetLorentzVector();
  TLorentzVector cLp1 = p1->GetLorentzVector(); cLp1.Boost(-boost);
  TLorentzVector Lp2 = p2->GetLorentzVector();
  TLorentzVector cLp2 = p2->GetLorentzVector(); cLp2.Boost(-boost);
  TLorentzVector Lpip1 = pip1->GetLorentzVector();
  TLorentzVector cLpip1 = pip1->GetLorentzVector(); cLpip1.Boost(-boost);
  TLorentzVector Llam = lam->GetLorentzVector();
  TLorentzVector cLlam = lam->GetLorentzVector(); cLlam.Boost(-boost);
  TLorentzVector Lpipp, cLpipp;
  TLorentzVector Lp, cLp;
  if(pip1flag){
    Lpipp = Llam + Lp2;
    cLpipp = cLlam + cLp2;
    Lp = Lp2;
    cLp = cLp2;
  }

  double pim_mass[2]  = { Lpim.M()                           , cLpim.M()};
  double pim_mom[2]   = { Lpim.Vect().Mag()                  , cLpim.Vect().Mag()};
  double pim_cost[2]  = { Lpim.Vect().CosTheta()             , cLpim.Vect().CosTheta()};
  double pim_phi[2]   = { Lpim.Vect().Phi()                  , cLpim.Vect().Phi()};
  double pim_mmass[2] = { (Ltgt+Lbeam-Lpim).M()              , (cLtgt+cLbeam-cLpim).M()};
  double pim_mmom[2]  = { (Ltgt+Lbeam-Lpim).Vect().Mag()     , (cLtgt+cLbeam-cLpim).Vect().Mag()};
  double pim_mcost[2] = { (Ltgt+Lbeam-Lpim).Vect().CosTheta(), (cLtgt+cLbeam-cLpim).Vect().CosTheta()};

  double p_mass[2]  = { Lp.M()                           , cLp.M()};
  double p_mom[2]   = { Lp.Vect().Mag()                  , cLp.Vect().Mag()};
  double p_cost[2]  = { Lp.Vect().CosTheta()             , cLp.Vect().CosTheta()};
  double p_phi[2]   = { Lp.Vect().Phi()                  , cLp.Vect().Phi()};
  double p_mmass[2] = { (Ltgt+Lbeam-Lp).M()              , (cLtgt+cLbeam-cLp).M()};
  double p_mmom[2]  = { (Ltgt+Lbeam-Lp).Vect().Mag()     , (cLtgt+cLbeam-cLp).Vect().Mag()};
  double p_mcost[2] = { (Ltgt+Lbeam-Lp).Vect().CosTheta(), (cLtgt+cLbeam-cLp).Vect().CosTheta()};

  double p1_mass[2]  = { Lp1.M()                           , cLp1.M()};
  double p1_mom[2]   = { Lp1.Vect().Mag()                  , cLp1.Vect().Mag()};
  double p1_cost[2]  = { Lp1.Vect().CosTheta()             , cLp1.Vect().CosTheta()};
  double p1_phi[2]   = { Lp1.Vect().Phi()                  , cLp1.Vect().Phi()};
  double p1_mmass[2] = { (Ltgt+Lbeam-Lp1).M()              , (cLtgt+cLbeam-cLp1).M()};
  double p1_mmom[2]  = { (Ltgt+Lbeam-Lp1).Vect().Mag()     , (cLtgt+cLbeam-cLp1).Vect().Mag()};
  double p1_mcost[2] = { (Ltgt+Lbeam-Lp1).Vect().CosTheta(), (cLtgt+cLbeam-cLp1).Vect().CosTheta()};

  double pip1_mass[2]  = { Lpip1.M()                           , cLpip1.M()};
  double pip1_mom[2]   = { Lpip1.Vect().Mag()                  , cLpip1.Vect().Mag()};
  double pip1_cost[2]  = { Lpip1.Vect().CosTheta()             , cLpip1.Vect().CosTheta()};
  double pip1_phi[2]   = { Lpip1.Vect().Phi()                  , cLpip1.Vect().Phi()};
  double pip1_mmass[2] = { (Ltgt+Lbeam-Lpip1).M()              , (cLtgt+cLbeam-cLpip1).M()};
  double pip1_mmom[2]  = { (Ltgt+Lbeam-Lpip1).Vect().Mag()     , (cLtgt+cLbeam-cLpip1).Vect().Mag()};
  double pip1_mcost[2] = { (Ltgt+Lbeam-Lpip1).Vect().CosTheta(), (cLtgt+cLbeam-cLpip1).Vect().CosTheta()};

  double lam_mass[2]  = { Llam.M()                           , cLlam.M()};
  double lam_mom[2]   = { Llam.Vect().Mag()                  , cLlam.Vect().Mag()};
  double lam_cost[2]  = { Llam.Vect().CosTheta()             , cLlam.Vect().CosTheta()};
  double lam_phi[2]   = { Llam.Vect().Phi()                  , cLlam.Vect().Phi()};
  double lam_mmass[2] = { (Ltgt+Lbeam-Llam).M()              , (cLtgt+cLbeam-cLlam).M()};
  double lam_mmom[2]  = { (Ltgt+Lbeam-Llam).Vect().Mag()     , (cLtgt+cLbeam-cLlam).Vect().Mag()};
  double lam_mcost[2] = { (Ltgt+Lbeam-Llam).Vect().CosTheta(), (cLtgt+cLbeam-cLlam).Vect().CosTheta()};

  /* PC */
  TLorentzVector Lpc  = p2 ->GetLorentzVector();
  TLorentzVector cLpc  = p2 ->GetLorentzVector(); cLpc.Boost(-boost);
  double nc_mom[2]     = { Lpc.Vect().Mag()                  , cLpc.Vect().Mag()};
  double nc_cost[2]    = { Lpc.Vect().CosTheta()             , cLpc.Vect().CosTheta()};
  double nc_phi[2]     = { Lpc.Vect().Phi()                  , cLpc.Vect().Phi()};
  double nc_mmass[2]   = { (Ltgt+Lbeam-Lpc).M()              , (cLtgt+cLbeam-cLpc).M()};
  double nc_mmom[2]    = { (Ltgt+Lbeam-Lpc).Vect().Mag()     , (cLtgt+cLbeam-cLpc).Vect().Mag()};
  double nc_mcost[2]   = { (Ltgt+Lbeam-Lpc).Vect().CosTheta(), (cLtgt+cLbeam-cLpc).Vect().CosTheta()};
  double pnc_mass[2]    = { (Lp+Lpc).M()                           , (cLp+cLpc).M()};
  double pnc_mom[2]     = { (Lp+Lpc).Vect().Mag()                  , (cLp+cLpc).Vect().Mag()};
  double pnc_cost[2]    = { (Lp+Lpc).Vect().CosTheta()             , (cLp+cLpc).Vect().CosTheta()};
  double pnc_phi[2]     = { (Lp+Lpc).Vect().Phi()                  , (cLp+cLpc).Vect().Phi()};
  double pnc_mmass[2]   = { (Ltgt+Lbeam-Lp-Lpc).M()              , (cLtgt+cLbeam-cLp-cLpc).M()};
  double pnc_mmass2[2]   = { (Ltgt+Lbeam-Lp-Lpc).M2()              , (cLtgt+cLbeam-cLp-cLpc).M2()};
  double pnc_mmom[2]    = { (Ltgt+Lbeam-Lp-Lpc).Vect().Mag()     , (cLtgt+cLbeam-cLp-cLpc).Vect().Mag()};
  double pnc_mcost[2]   = { (Ltgt+Lbeam-Lp-Lpc).Vect().CosTheta(), (cLtgt+cLbeam-cLp-cLpc).Vect().CosTheta()};
  double pipppc_mass[2]    = { (Lpipp+Lpc).M()                           , (cLpipp+cLpc).M()};
  double pipppc_mom[2]     = { (Lpipp+Lpc).Vect().Mag()                  , (cLpipp+cLpc).Vect().Mag()};
  double pipppc_cost[2]    = { (Lpipp+Lpc).Vect().CosTheta()             , (cLpipp+cLpc).Vect().CosTheta()};
  double pipppc_phi[2]     = { (Lpipp+Lpc).Vect().Phi()                  , (cLpipp+cLpc).Vect().Phi()};
  double pipppc_mmass[2]   = { (Ltgt+Lbeam-Lpipp-Lpc).M()              , (cLtgt+cLbeam-cLpipp-cLpc).M()};
  double pipppc_mmass2[2]   = { (Ltgt+Lbeam-Lpipp-Lpc).M2()              , (cLtgt+cLbeam-cLpipp-cLpc).M2()};
  double pipppc_mmom[2]    = { (Ltgt+Lbeam-Lpipp-Lpc).Vect().Mag()     , (cLtgt+cLbeam-cLpipp-cLpc).Vect().Mag()};
  double pipppc_mcost[2]   = { (Ltgt+Lbeam-Lpipp-Lpc).Vect().CosTheta(), (cLtgt+cLbeam-cLpipp-cLpc).Vect().CosTheta()};

  /* flag */
  bool mnflag = false; /* mm flag */
  bool sblmnflag = false; /* lower sideband of mm flag */
  bool sbumnflag = false; /* upper sideband of mm flag */
  {
    double pro_mmass = pipp_mmass[1];
    if(mnll<pro_mmass&&pro_mmass<mnul) mnflag = true;
    if(sblmnll<pro_mmass&&pro_mmass<sblmnul) sblmnflag = true;
    if(sbumnll<pro_mmass&&pro_mmass<sbumnul) sbumnflag = true;
  }

  //if(!lamflag1&&!lamflag2){ return true; }
  //if(lamflag1){ return true; }
  //if(!lamflag2){ return true; }
  if(!mnflag){ return true; }
  //if(mnflag){ return true; }
  FillHist("EventNumber",10); /* Missing Neutron */

  /*-----------------------------------------*/
  /* All                                     */ 
  /*-----------------------------------------*/
  // ######################### //
  // Invariant/Missing Mass    //
  // ######################### //
  // pi-
  for(int f=1; f<2; f++){
    FillHist(Form("pim_Momentum_%s",frame[f].Data()),pim_mom[f]);
    FillHist(Form("pim_CosT_%s",frame[f].Data()),pim_cost[f]);
    FillHist(Form("pim_Phi_%s",frame[f].Data()),pim_phi[f]);
    FillHist(Form("pim_MMass_%s",frame[f].Data()),pim_mmass[f]);
    FillHist(Form("pim_MMomentum_%s",frame[f].Data()),pim_mmom[f]);
    FillHist(Form("pim_MCosT_%s",frame[f].Data()),pim_mcost[f]);
  }
  // proton
  for(int f=1; f<2; f++){
    FillHist(Form("p_Momentum_%s",frame[f].Data()),p_mom[f]);
    FillHist(Form("p_CosT_%s",frame[f].Data()),p_cost[f]);
    FillHist(Form("p_Phi_%s",frame[f].Data()),p_phi[f]);
    FillHist(Form("p_MMass_%s",frame[f].Data()),p_mmass[f]);
    FillHist(Form("p_MMomentum_%s",frame[f].Data()),p_mmom[f]);
    FillHist(Form("p_MCosT_%s",frame[f].Data()),p_mcost[f]);
  }
  // proton1
  for(int f=1; f<2; f++){
    FillHist(Form("p1_Momentum_%s",frame[f].Data()),p1_mom[f]);
    FillHist(Form("p1_CosT_%s",frame[f].Data()),p1_cost[f]);
    FillHist(Form("p1_Phi_%s",frame[f].Data()),p1_phi[f]);
    FillHist(Form("p1_MMass_%s",frame[f].Data()),p1_mmass[f]);
    FillHist(Form("p1_MMomentum_%s",frame[f].Data()),p1_mmom[f]);
    FillHist(Form("p1_MCosT_%s",frame[f].Data()),p1_mcost[f]);
  }
  // proton2
  for(int f=1; f<2; f++){
    FillHist(Form("p2_Momentum_%s",frame[f].Data()),p2_mom[f]);
    FillHist(Form("p2_CosT_%s",frame[f].Data()),p2_cost[f]);
    FillHist(Form("p2_Phi_%s",frame[f].Data()),p2_phi[f]);
    FillHist(Form("p2_MMass_%s",frame[f].Data()),p2_mmass[f]);
    FillHist(Form("p2_MMomentum_%s",frame[f].Data()),p2_mmom[f]);
    FillHist(Form("p2_MCosT_%s",frame[f].Data()),p2_mcost[f]);
  }
  // pi- proton1
  for(int f=1; f<2; f++){
    FillHist(Form("pip1_Mass_%s",frame[f].Data()),pip1_mass[f]);
    FillHist(Form("pip1_Momentum_%s",frame[f].Data()),pip1_mom[f]);
    FillHist(Form("pip1_CosT_%s",frame[f].Data()),pip1_cost[f]);
    FillHist(Form("pip1_Phi_%s",frame[f].Data()),pip1_phi[f]);
    FillHist(Form("pip1_MMass_%s",frame[f].Data()),pip1_mmass[f]);
    FillHist(Form("pip1_MMomentum_%s",frame[f].Data()),pip1_mmom[f]);
    FillHist(Form("pip1_MCosT_%s",frame[f].Data()),pip1_mcost[f]);
  }
  // pi- proton2
  for(int f=1; f<2; f++){
    FillHist(Form("pip2_Mass_%s",frame[f].Data()),pip2_mass[f]);
    FillHist(Form("pip2_Momentum_%s",frame[f].Data()),pip2_mom[f]);
    FillHist(Form("pip2_CosT_%s",frame[f].Data()),pip2_cost[f]);
    FillHist(Form("pip2_Phi_%s",frame[f].Data()),pip2_phi[f]);
    FillHist(Form("pip2_MMass_%s",frame[f].Data()),pip2_mmass[f]);
    FillHist(Form("pip2_MMomentum_%s",frame[f].Data()),pip2_mmom[f]);
    FillHist(Form("pip2_MCosT_%s",frame[f].Data()),pip2_mcost[f]);
  }
  // lambda
  for(int f=1; f<2; f++){
    FillHist(Form("lam_Mass_%s",frame[f].Data()),lam_mass[f]);
    FillHist(Form("lam_Momentum_%s",frame[f].Data()),lam_mom[f]);
    FillHist(Form("lam_CosT_%s",frame[f].Data()),lam_cost[f]);
    FillHist(Form("lam_Phi_%s",frame[f].Data()),lam_phi[f]);
    FillHist(Form("lam_MMass_%s",frame[f].Data()),lam_mmass[f]);
    FillHist(Form("lam_MMomentum_%s",frame[f].Data()),lam_mmom[f]);
    FillHist(Form("lam_MCosT_%s",frame[f].Data()),lam_mcost[f]);
  }
  // pi- proton1 proton2
  for(int f=1; f<2; f++){
    FillHist(Form("pipp_Mass_%s",frame[f].Data()),pipp_mass[f]);
    FillHist(Form("pipp_Momentum_%s",frame[f].Data()),pipp_mom[f]);
    FillHist(Form("pipp_CosT_%s",frame[f].Data()),pipp_cost[f]);
    FillHist(Form("pipp_Phi_%s",frame[f].Data()),pipp_phi[f]);
    FillHist(Form("pipp_MMass_%s",frame[f].Data()),pipp_mmass[f]);
    FillHist(Form("pipp_MMomentum_%s",frame[f].Data()),pipp_mmom[f]);
    FillHist(Form("pipp_MCosT_%s",frame[f].Data()),pipp_mcost[f]);
    FillHist(Form("pipp_MTheta_%s",frame[f].Data()),TMath::ACos(pipp_mcost[f])*180.0/TMath::Pi());
  }
  // Momentum transfer
  for(int f=1; f<2; f++){
    FillHist(Form("pipp_MomentumTransfer_%s",frame[f].Data()),q[f]);
  }

  /*-----------------------------------------*/
  /* 2D Plots                                */ 
  /*-----------------------------------------*/
  // IM(pip1) vs. IM(pip2)
  for(int f=1; f<2; f++){
    FillHist(Form("pip1_IMvspip2_IM_%s",frame[f].Data()),pip1_mass[f],pip2_mass[f]);
  }
  // IM(pip1) vs. IM(pip2)
  for(int f=1; f<2; f++){
    FillHist(Form("pip1_IM2vspip2_IM2_%s",frame[f].Data()),pip1_mass[f]*pip1_mass[f],pip2_mass[f]*pip2_mass[f]);
  }
  // IM(pipp) vs. MM(pipp)
  for(int f=1; f<2; f++){
    FillHist(Form("pipp_IMvsHeKpipp_MM_%s",frame[f].Data()),pipp_mass[f],pipp_mmass[f]);
  }
  // IM(pipp) vs. MCosT(pipp)
  for(int f=1; f<2; f++){
    FillHist(Form("pipp_IMvsHeKpipp_MCosT_%s",frame[f].Data()),pipp_mass[f],pipp_mcost[f]);
  }
  // IM(pipp) vs. MTheta(pipp)
  for(int f=1; f<2; f++){
    FillHist(Form("pipp_IMvsHeKpipp_MTheta_%s",frame[f].Data()),pipp_mass[f],TMath::ACos(pipp_mcost[f])*180.0/TMath::Pi());
  }
  // IM(pipp) vs. MomentumTransfer
  for(int f=1; f<2; f++){
    FillHist(Form("pipp_IMvsMomentumTransfer_%s",frame[f].Data()),pipp_mass[f],q[f]);
  }
  // Dalitz Plot
  for(int f=1; f<2; f++){
    FillHist(Form("pipp_DalitzPlot_%s",frame[f].Data()),(tp-tl)/sqrt(3)/qvalue,tn/qvalue);
  }

  if(nc==0) return true;
  FillHist("EventNumber",11); /* NC Hit */

  // n
  for(int f=1; f<2; f++){
    FillHist(Form("n_Momentum_%s",frame[f].Data()),nc_mom[f]);
    FillHist(Form("n_CosT_%s",frame[f].Data()),nc_cost[f]);
    FillHist(Form("n_Phi_%s",frame[f].Data()),nc_phi[f]);
    FillHist(Form("n_MMass_%s",frame[f].Data()),nc_mmass[f]);
    FillHist(Form("n_MMomentum_%s",frame[f].Data()),nc_mmom[f]);
    FillHist(Form("n_MCosT_%s",frame[f].Data()),nc_mcost[f]);
  }
  // n + p
  for(int f=1; f<2; f++){
    FillHist(Form("pn_Mass_%s",frame[f].Data()),pnc_mass[f]);
    FillHist(Form("pn_Momentum_%s",frame[f].Data()),pnc_mom[f]);
    FillHist(Form("pn_CosT_%s",frame[f].Data()),pnc_cost[f]);
    FillHist(Form("pn_Phi_%s",frame[f].Data()),pnc_phi[f]);
    FillHist(Form("pn_MMass_%s",frame[f].Data()),pnc_mmass[f]);
    FillHist(Form("pn_MMass2_%s",frame[f].Data()),pnc_mmass2[f]);
    FillHist(Form("pn_MMomentum_%s",frame[f].Data()),pnc_mmom[f]);
    FillHist(Form("pn_MCosT_%s",frame[f].Data()),pnc_mcost[f]);
  }
  // n + pipp
  for(int f=1; f<2; f++){
    FillHist(Form("pippn_Mass_%s",frame[f].Data()),pippnc_mass[f]);
    FillHist(Form("pippn_Momentum_%s",frame[f].Data()),pippnc_mom[f]);
    FillHist(Form("pippn_CosT_%s",frame[f].Data()),pippnc_cost[f]);
    FillHist(Form("pippn_Phi_%s",frame[f].Data()),pippnc_phi[f]);
    FillHist(Form("pippn_MMass_%s",frame[f].Data()),pippnc_mmass[f]);
    FillHist(Form("pippn_MMass2_%s",frame[f].Data()),pippnc_mmass2[f]);
    FillHist(Form("pippn_MMomentum_%s",frame[f].Data()),pippnc_mmom[f]);
    FillHist(Form("pippn_MCosT_%s",frame[f].Data()),pippnc_mcost[f]);
  }
  // pipp (FWD)
  for(int f=1; f<2; f++){
    FillHist(Form("pipp_Mass_n_%s",frame[f].Data()),pipp_mass[f]);
    FillHist(Form("pipp_Momentum_n_%s",frame[f].Data()),pipp_mom[f]);
    FillHist(Form("pipp_CosT_n_%s",frame[f].Data()),pipp_cost[f]);
    FillHist(Form("pipp_Phi_n_%s",frame[f].Data()),pipp_phi[f]);
    FillHist(Form("pipp_MMass_n_%s",frame[f].Data()),pipp_mmass[f]);
    FillHist(Form("pipp_MMomentum_n_%s",frame[f].Data()),pipp_mmom[f]);
    FillHist(Form("pipp_MCosT_n_%s",frame[f].Data()),pipp_mcost[f]);
    FillHist(Form("pipp_MTheta_n_%s",frame[f].Data()),TMath::ACos(pipp_mcost[f])*180.0/TMath::Pi());
  }
  // Momentum transfer (FWD)
  for(int f=1; f<2; f++){
    FillHist(Form("pipp_MomentumTransfer_n_%s",frame[f].Data()),q_2);
  }
  /*-----------------------------------------*/
  /* 2D Plots                                */ 
  /*-----------------------------------------*/
  // IM(pipp) vs. MM2(pippn)
  for(int f=1; f<2; f++){
    FillHist(Form("pipp_IMvsHeKpippn_MM2_%s",frame[f].Data()),pipp_mass[f],pippnc_mmass2[f]);
  }
  // IM(pipp) vs. MM(n)
  for(int f=1; f<2; f++){
    FillHist(Form("pipp_IMvsHeKn_MM_%s",frame[f].Data()),pipp_mass[f],nc_mmass[f]);
  }
  // MMom(pipp) vs. Mom(n)
  for(int f=1; f<2; f++){
    FillHist(Form("HeKpipp_MMomvsn_Mom_%s",frame[f].Data()),pipp_mom[f],nc_mom[f]);
  }
  // MM(pipp) vs. MM2(pippn)
  for(int f=1; f<2; f++){
    FillHist(Form("HeKpipp_MMvsHeKpippn_MM2_%s",frame[f].Data()),pipp_mmass[f],pippnc_mmass2[f]);
  }
  // MM(pn) vs. MM2(pippn)
  for(int f=1; f<2; f++){
    FillHist(Form("HeKpn_MMvsHeKpippn_MM2_%s",frame[f].Data()),pnc_mmass[f],pippnc_mmass2[f]);
  }
  // IM(pipp) vs. MM(pn) 
  for(int f=1; f<2; f++){
    FillHist(Form("pipp_IMvsHeKpn_MM_%s",frame[f].Data()),pipp_mass[f],pnc_mmass[f]);
  }
  // MM(n) vs. MM(pn) 
  for(int f=1; f<2; f++){
    FillHist(Form("HeKn_MMvsHeKpn_MM_%s",frame[f].Data()),nc_mmass[f],pnc_mmass[f]);
  }
  // MM2(pn) vs. MM2(pippn)
  for(int f=1; f<2; f++){
    FillHist(Form("HeKpn_MM2vsHeKpippn_MM2_%s",frame[f].Data()),pnc_mmass2[f],pippnc_mmass2[f]);
  }
  // IM(pipp) vs. MM(pipp) (FWD)
  for(int f=1; f<2; f++){
    FillHist(Form("pipp_IMvsHeKpipp_MM_n_%s",frame[f].Data()),pipp_mass[f],pipp_mmass[f]);
  }
  // IM(pipp) vs. MCosT(pipp) (FWD)
  for(int f=1; f<2; f++){
    FillHist(Form("pipp_IMvsHeKpipp_MCosT_n_%s",frame[f].Data()),pipp_mass[f],pipp_mcost[f]);
  }
  // IM(pipp) vs. MomentumTransfer (FWD)
  for(int f=1; f<2; f++){
    FillHist(Form("pipp_IMvsMomentumTransfer_n_%s",frame[f].Data()),pipp_mass[f],q_2);
  }
  // Dalitz Plot (FWD)
  for(int f=1; f<2; f++){
    FillHist(Form("pipp_DalitzPlot_n_%s",frame[f].Data()),(tp-tl)/sqrt(3)/qvalue,tn/qvalue);
  }

  return true;

}

bool MyAnalysisHeKpippfwd::FillHist(TString name, double val1, int weight)
{
  TH1F* h1 = (TH1F*)gFile -> Get(name);
  if(h1&&weight>0){ 
    for(int i=0; i<weight; i++){
      h1 -> Fill(val1,1);
    }
    return true;
  }
  else {
    return false;
  }
}

bool MyAnalysisHeKpippfwd::FillHist(TString name, TString val1, int weight)
{
  TH1F* h1 = (TH1F*)gFile -> Get(name);
  if(h1&&weight>0){
    for(int i=0; i<weight; i++){
      h1 -> Fill(val1,1);
    }
    return true;
  }
  else {
    return false;
  }
}

bool MyAnalysisHeKpippfwd::FillHist(TString name, double val1, double val2, int weight)
{
  TH2F* h2 = (TH2F*)gFile -> Get(name);
  if(h2&&weight>0){
    for(int i=0; i<weight; i++){
      h2 -> Fill(val1,val2,1);
    }
    return true;
  }
  else {
    return false;
  }
}

bool MyAnalysisHeKpippfwd::FillHist(TString name, TString val1, TString val2, int weight)
{
  TH2F* h2 = (TH2F*)gFile -> Get(name);
  if(h2&&weight>0){
    for(int i=0; i<weight; i++){
      h2 -> Fill(val1,val2,1);
    }
    return true;
  }
  else {
    return false;
  }
}

void MyAnalysisHeKpippfwd::CutCondition()
{
  /* K0 mass cut */ 
  double mean  = 0.4976;
  double sigma = 0.0069;
  k0ll = mean-2*sigma; k0ul = mean+2*sigma;
  sblk0ll = mean-5*sigma; sblk0ul = mean-3*sigma;
  sbuk0ll = mean+3*sigma; sbuk0ul = mean+5*sigma;
  /* Missing neutron mass cut */ 
  mean  = 0.9400;
  sigma = 0.0450;
  mnll = 0.85; mnul = 1.03;
  //mnll = mean-2*sigma; mnul = mean+2*sigma;
  sblmnll = mean-5*sigma; sblmnul = mean-3*sigma;
  sbumnll = mean+3*sigma; sbumnul = mean+5*sigma;
  /* SigmaPlus mass cut */
  mean  = 1.1888;
  sigma = 0.0044;
  spll = mean-2*sigma; spul = mean+2*sigma;
  sblspll = mean-5*sigma; sblspul = mean-3*sigma;
  sbuspll = mean+3*sigma; sbuspul = mean+5*sigma;
  /* SigmaMinus mass cut */
  mean  = 1.1970;
  sigma = 0.0054;
  smll = mean-2*sigma; smul = mean+2*sigma;
  sblsmll = mean-5*sigma; sblsmul = mean-3*sigma;
  sbusmll = mean+3*sigma; sbusmul = mean+5*sigma;
  /* Lambda mass cut */ 
  mean  = 1.1155;
  sigma = 0.0020;
  lamll = mean-2*sigma; lamul = mean+2*sigma;
  sbllamll = mean-5*sigma; sbllamul = mean-3*sigma;
  sbulamll = mean+3*sigma; sbulamul = mean+5*sigma;
}

bool MyAnalysisHeKpippfwd::Initialize(ConfMan* confMan)
{
  std::cout << "### MyAnalysisHeKpippfwd::Initialize ###" << std::endl;

  std::string ofname = confMan->GetOutFileName();
  ofname.insert(ofname.find(".root"),"_anaHeKpipp");

  rtFile =  new TFile( Form("%s",ofname.c_str()), "RECREATE");
  rtFile -> cd();

  TString pname[8]  = {"pip","p","d","t","he","pim","k","o"};
  TString pname2[8] = {"#pi^{+}","p","d","t","he","#pi^{-}","K^{-}","o"};
  TString proname[6]  = {"pipi","k0","sp","sm","s","l"};
  TString proname2[6] = {"#pi^{+}#pi^{-}","K^{0}","#Sigma^{+}","#Sigma^{-}","#Sigma^{#pm}","#Lambda"};

  new TH1F( Form("CDC_Chi2"), Form("#chi^{2} CDC;#chi^{2};Counts"), 5000, 0.0, 100.0 );
  new TH1F( "EventNumber", "Number of Events", 20, 0, 20);
  TH1F* h1 = new TH1F( "TriggerPattern", "Trigger Pattern", 18, -0.5, 17.5);
  h1->GetXaxis()->SetBinLabel(1,"All");
  h1->GetXaxis()->SetBinLabel(2,"Beam");
  h1->GetXaxis()->SetBinLabel(3,"Kaon");
  h1->GetXaxis()->SetBinLabel(4,"KCDH1f");
  h1->GetXaxis()->SetBinLabel(5,"Pion");
  h1->GetXaxis()->SetBinLabel(6,"Proton");
  h1->GetXaxis()->SetBinLabel(7,"KCDH1");
  h1->GetXaxis()->SetBinLabel(8,"KCDH2");
  h1->GetXaxis()->SetBinLabel(9,"PivBVC");
  h1->GetXaxis()->SetBinLabel(10,"PiCDH1");
  h1->GetXaxis()->SetBinLabel(11,"PiCDH2");
  h1->GetXaxis()->SetBinLabel(12,"Kf");
  h1->GetXaxis()->SetBinLabel(13,"1stMix");
  h1->GetXaxis()->SetBinLabel(14,"Charged");
  h1->GetXaxis()->SetBinLabel(15,"Neutral");
  h1->GetXaxis()->SetBinLabel(16,"Cosmic");
  h1->GetXaxis()->SetBinLabel(17,"Reject");
  h1->GetXaxis()->SetBinLabel(18,"SIM");

  // CDS Particles
  std::cout << "Define Histograms for CDS particles" << std::endl;
  new TH1F( "CDS_NumOfParticle", "Number of CDS tracks", 10, 0, 10);
  TH2F* h2;
  h2 = new TH2F("CDS_Particle","Detected particles by CDS;Particle;Number of all tracks",6,0,6,6,0,6);
  h2->GetYaxis()->SetBinLabel(1,"1 track");
  h2->GetYaxis()->SetBinLabel(2,"2 track");
  h2->GetYaxis()->SetBinLabel(3,"3 track");
  h2->GetYaxis()->SetBinLabel(4,"4 track");
  h2->GetYaxis()->SetBinLabel(5,"5 track");
  h2->GetYaxis()->SetBinLabel(6,"6 track");
  h2->GetXaxis()->SetBinLabel(1,"#pi^{+}");
  h2->GetXaxis()->SetBinLabel(2,"#pi^{-}");
  h2->GetXaxis()->SetBinLabel(3,"K^{-}");
  h2->GetXaxis()->SetBinLabel(4,"p");
  h2->GetXaxis()->SetBinLabel(5,"d");
  h2->GetXaxis()->SetBinLabel(6,"Other");
  // DCA
  std::cout << "Define Histograms for DCA" << std::endl;
  new TH1F( "pim_VBDIS", "VBDIS of #pi^{-};VBDIS (cm);Coutns", 100, 0, 10);
  new TH1F( "p1_VBDIS", "VBDIS of p_{1};VBDIS (cm);Coutns", 100, 0, 10);
  new TH1F( "p2_VBDIS", "VBDIS of p_{2};VBDIS (cm);Coutns", 100, 0, 10);
  new TH1F( "p_VBDIS", "VBDIS of p;VBDIS (cm);Coutns", 100, 0, 10);
  new TH2F( "p1_VBDISvsp2_VBDIS", "VBDIS of p_{1} vs p_{2};VBDIS (cm);VBDIS (cm)", 100, 0, 10, 100, 0, 10);
  new TH1F( "pip1_VDIS", "VDIS of #pi^{-}p_{1};VDIS (cm);Coutns", 100, 0, 10);
  new TH1F( "pip2_VDIS", "VDIS of #pi^{-}p_{2};VDIS (cm);Coutns", 100, 0, 10);
  new TH1F( "lam_VDIS", "VDIS of #Lambda;VDIS (cm);Coutns", 100, 0, 10);
  new TH2F( "pip1_VDISvspip2_VDIS", "VDIS of #pi^{-}p_{1} vs #pi^{-}p_{2};VDIS (cm);VDIS (cm)", 100, 0, 10, 10, 0, 10);
  new TH1F( "pip1_VBDIS", "VBDIS of #pi^{-}p_{1};VBDIS (cm);Coutns", 100, 0, 10);
  new TH1F( "pip2_VBDIS", "VBDIS of #pi^{-}p_{2};VBDIS (cm);Coutns", 100, 0, 10);
  new TH1F( "lam_VBDIS", "VBDIS of #pi^{-}p_{2};VBDIS (cm);Coutns", 100, 0, 10);
  new TH2F( "pip1_VBDISvspip2_VBDIS", "VBDIS of #pi^{-}p_{1} vs #pi^{-}p_{2};VBDIS (cm);VBDIS (cm)", 100, 0, 10, 100, 0, 10);
  new TH1F( "pip1_PBDCA", "PBDCA of #pi^{-}p_{1};PBDCA (cm);Coutns", 100, 0, 10);
  new TH1F( "pip2_PBDCA", "PBDCA of #pi^{-}p_{2};PBDCA (cm);Coutns", 100, 0, 10);
  new TH1F( "lam_PBDCA", "PBDCA of #Lambda;PBDCA (cm);Coutns", 100, 0, 10);
  new TH2F( "pip1_PBDCAvspip2_PBDCA", "PBDCA of #pi^{-}p_{1} vs #pi^{-}p_{2};PBDCA (cm);PBDCA (cm)", 100, 0, 10, 100, 0, 10);
  new TH2F( "pim_Vertex_XY", "Vertex XY plane #pi^{-};X (cm);Y (cm)", 300, -15, 15, 300, -15, 15);
  new TH2F( "pim_Vertex_ZX", "Vertex ZX plane #pi^{-};Z (cm);X (cm)", 600, -30, 30, 300, -15, 15);
  new TH2F( "pim_Vertex_ZY", "Vertex ZY plane #pi^{-};Z (cm);Y (cm)", 600, -30, 30, 300, -15, 15);
  new TH2F( "p1_Vertex_XY", "Vertex XY plane p_{1};X (cm);Y (cm)", 300, -15, 15, 300, -15, 15);
  new TH2F( "p1_Vertex_ZX", "Vertex ZX plane p_{1};Z (cm);X (cm)", 600, -30, 30, 300, -15, 15);
  new TH2F( "p1_Vertex_ZY", "Vertex ZY plane p_{1};Z (cm);Y (cm)", 600, -30, 30, 300, -15, 15);
  new TH2F( "p2_Vertex_XY", "Vertex XY plane p_{2};X (cm);Y (cm)", 300, -15, 15, 300, -15, 15);
  new TH2F( "p2_Vertex_ZX", "Vertex ZX plane p_{2};Z (cm);X (cm)", 600, -30, 30, 300, -15, 15);
  new TH2F( "p2_Vertex_ZY", "Vertex ZY plane p_{2};Z (cm);Y (cm)", 600, -30, 30, 300, -15, 15);
  new TH2F( "pip1_Vertex_XY", "Vertex XY plane #pi^{-}p_{1};X (cm);Y (cm)", 300, -15, 15, 300, -15, 15);
  new TH2F( "pip1_Vertex_ZX", "Vertex ZX plane #pi^{-}p_{1};Z (cm);X (cm)", 600, -30, 30, 300, -15, 15);
  new TH2F( "pip1_Vertex_ZY", "Vertex ZY plane #pi^{-}p_{1};Z (cm);Y (cm)", 600, -30, 30, 300, -15, 15);
  new TH2F( "pip2_Vertex_XY", "Vertex XY plane #pi^{-}p_{2};X (cm);Y (cm)", 300, -15, 15, 300, -15, 15);
  new TH2F( "pip2_Vertex_ZX", "Vertex ZX plane #pi^{-}p_{2};Z (cm);X (cm)", 600, -30, 30, 300, -15, 15);
  new TH2F( "pip2_Vertex_ZY", "Vertex ZY plane #pi^{-}p_{2};Z (cm);Y (cm)", 600, -30, 30, 300, -15, 15);
  new TH2F( "lam_Vertex_XY", "Vertex XY plane #Lambda;X (cm);Y (cm)", 300, -15, 15, 300, -15, 15);
  new TH2F( "lam_Vertex_ZX", "Vertex ZX plane #Lambda;Z (cm);X (cm)", 600, -30, 30, 300, -15, 15);
  new TH2F( "lam_Vertex_ZY", "Vertex ZY plane #Lambda;Z (cm);Y (cm)", 600, -30, 30, 300, -15, 15);
  new TH2F( "lam_PVertex_XY", "Vertex XY plane #Lambda;X (cm);Y (cm)", 300, -15, 15, 300, -15, 15);
  new TH2F( "lam_PVertex_ZX", "Vertex ZX plane #Lambda;Z (cm);X (cm)", 600, -30, 30, 300, -15, 15);
  new TH2F( "lam_PVertex_ZY", "Vertex ZY plane #Lambda;Z (cm);Y (cm)", 600, -30, 30, 300, -15, 15);

  //// FWD neutral Particles
  //std::cout << "Define Histograms for FWD neutral particles" << std::endl;
  //new TH2F( "FWD_OverbetavsMomentum", "1/#beta vs. Momentum;1/#beta;Momentum (GeV/c)", 1000, -0, 10, 500, -5, 5 );
  //new TH2F( "FWD_OverbetavsEnergy", "1/#beta vs. Energy deposit;1/#beta;Energy deposit (MeVee)", 1000, -0, 10, 5000, 0, 100);
  //new TH2F( "FWD_TOFvsMomentum", "TOF vs. Momentum;TOF (ns);Momentum (GeV/c)", 4000, 30, 70, 500, -5, 5 );
  //new TH2F( "FWD_HitPosition", "Hit position at NC (X vs. Y);X position (cm);Y position (cm)",40,-400,400,25,-250,250);
  //new TH2F( "FWD_HitSegment", "NC segment;Segment in layer;Layer",16,0.5,16.5,7,0.5,7.5);


  /*-----------------------------------------*/
  /* All                                     */ 
  /*-----------------------------------------*/
  std::cout << "Define Histograms for IM/MM" << std::endl;
  /* 1D Plots */
  // pi-
  new TH1F("pim_Momentum_CM","Momentum of #pi^{-};Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pim_CosT_CM","cos#theta of #pi^{-};cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("pim_Phi_CM","#phi of #pi^{-};#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("pim_MMass_CM", "^{3}He(K^{-},#pi^{-})X missing mass;MM(#pi^{-}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pim_MMomentum_CM", "^{3}He(K^{-},#pi^{-})X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pim_MCosT_CM", "^{3}He(K^{-},#pi^{-})X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // proton1
  new TH1F("p1_Momentum_CM","Momentum of p_{1};Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("p1_CosT_CM","cos#theta of p_{1};cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("p1_Phi_CM","#phi of p_{1};#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("p1_MMass_CM", "^{3}He(K^{-},p_{1})X missing mass;MM(p_{1}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("p1_MMomentum_CM", "^{3}He(K^{-},p_{1})X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("p1_MCosT_CM", "^{3}He(K^{-},p_{1})X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // proton2
  new TH1F("p2_Momentum_CM","Momentum of p_{2};Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("p2_CosT_CM","cos#theta of p_{2};cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("p2_Phi_CM","#phi of p_{2};#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("p2_MMass_CM", "^{3}He(K^{-},p_{2})X missing mass;MM(p_{2}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("p2_MMomentum_CM", "^{3}He(K^{-},p_{2})X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("p2_MCosT_CM", "^{3}He(K^{-},p_{2})X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // pi-,proton1
  new TH1F("pip1_Mass_CM","Invariant mass of #pi^{-}p_{1};IM(#pi^{-}p_{1}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pip1_Momentum_CM","Momentum of #pi^{-}p_{1};Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pip1_CosT_CM","cos#theta of #pi^{-}p_{1};cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("pip1_Phi_CM","#phi of #pi^{-}p_{1};#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("pip1_MMass_CM", "^{3}He(K^{-},#pi^{-}p_{1})X missing mass;MM(#pi^{-}p_{1}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pip1_MMomentum_CM", "^{3}He(K^{-},#pi^{-}p_{1})X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pip1_MCosT_CM", "^{3}He(K^{-},#pi^{-}p_{1})X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // pi-,proton2
  new TH1F("pip2_Mass_CM","Invariant mass of #pi^{-}p_{2};IM(#pi^{-}p_{2}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pip2_Momentum_CM","Momentum of #pi^{-}p_{2};Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pip2_CosT_CM","cos#theta of #pi^{-}p_{2};cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("pip2_Phi_CM","#phi of #pi^{-}p_{2};#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("pip2_MMass_CM", "^{3}He(K^{-},#pi^{-}p_{2})X missing mass;MM(#pi^{-}p_{2}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pip2_MMomentum_CM", "^{3}He(K^{-},#pi^{-}p_{2})X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pip2_MCosT_CM", "^{3}He(K^{-},#pi^{-}p_{2})X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // proton
  new TH1F("p_Momentum_CM","Momentum of p;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("p_CosT_CM","cos#theta of p;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("p_Phi_CM","#phi of p;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("p_MMass_CM", "^{3}He(K^{-},p)X missing mass;MM(p) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("p_MMomentum_CM", "^{3}He(K^{-},p)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("p_MCosT_CM", "^{3}He(K^{-},p)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // Lambda
  new TH1F("lam_Mass_CM","Invariant mass of #Lambda;IM(#Lambda) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("lam_Momentum_CM","Momentum of #Lambda;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("lam_CosT_CM","cos#theta of #Lambda;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("lam_Phi_CM","#phi of #Lambda;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("lam_MMass_CM", "^{3}He(K^{-},#Lambda)X missing mass;MM(#Lambda) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("lam_MMomentum_CM", "^{3}He(K^{-},#Lambda)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("lam_MCosT_CM", "^{3}He(K^{-},#Lambda)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // pi-,proton2
  new TH1F("lam_Mass_CM","Invariant mass of #pi^{-}p_{2};IM(#pi^{-}p_{2}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("lam_Momentum_CM","Momentum of #pi^{-}p_{2};Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("lam_CosT_CM","cos#theta of #pi^{-}p_{2};cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("lam_Phi_CM","#phi of #pi^{-}p_{2};#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("lam_MMass_CM", "^{3}He(K^{-},#pi^{-}p_{2})X missing mass;MM(#pi^{-}p_{2}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("lam_MMomentum_CM", "^{3}He(K^{-},#pi^{-}p_{2})X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("lam_MCosT_CM", "^{3}He(K^{-},#pi^{-}p_{2})X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // n
  new TH1F("n_Mass_CM","Invariant mass of n;IM(n) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("n_Momentum_CM","Momentum of n;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("n_CosT_CM","cos#theta of n;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("n_Phi_CM","#phi of n;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("n_MMass_CM", "^{3}He(K^{-},n)X missing mass;MM(n) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("n_MMass2_CM", "^{3}He(K^{-},n)X missing mass^{2};MM^{2}(n) (GeV/c^{2})^{2};Counts", 6000, -1.0, 5.0);
  new TH1F("n_MMomentum_CM", "^{3}He(K^{-},n)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("n_MCosT_CM", "^{3}He(K^{-},n)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("n_MomentumTransfer_CM","Momentum transfer to n;Momentum transfer (GeV/c);Counts", 2000, 0.0, 2.0);
  // pi-,proton1,proton2
  new TH1F("pipp_Mass_CM","Invariant mass of #pi^{-}pp;IM(#pi^{-}pp) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pipp_Momentum_CM","Momentum of #pi^{-}pp;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pipp_CosT_CM","cos#theta of #pi^{-}pp;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("pipp_Phi_CM","#phi of #pi^{-}pp;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("pipp_MMass_CM", "^{3}He(K^{-},#pi^{-}pp)X missing mass;MM(#pi^{-}pp) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pipp_MMomentum_CM", "^{3}He(K^{-},#pi^{-}pp)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pipp_MCosT_CM", "^{3}He(K^{-},#pi^{-}pp)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("pipp_MTheta_CM", "^{3}He(K^{-},#pi^{-}pp)X #theta;#theta;Counts", 200, -10.0, 190.0);
  new TH1F("pipp_MomentumTransfer_CM","Momentum transfer to #pi^{-}pp;Momentum transfer (GeV/c);Counts", 2000, 0.0, 2.0);
  // proton + n
  new TH1F("pn_Mass_CM","Invariant mass of pn;IM(pn) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pn_Momentum_CM","Momentum of pn;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pn_CosT_CM","cos#theta of pn;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("pn_Phi_CM","#phi of pn;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("pn_MMass_CM", "^{3}He(K^{-},pn)X missing mass;MM(pn) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pn_MMass2_CM", "^{3}He(K^{-},pn)X missing mass^{2};MM^{2}(pn) (GeV/c^{2})^{2};Counts", 6000, -1.0, 2.0);
  new TH1F("pn_MMomentum_CM", "^{3}He(K^{-},pn)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pn_MCosT_CM", "^{3}He(K^{-},pn)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("pn_MomentumTransfer_CM","Momentum transfer to pn;Momentum transfer (GeV/c);Counts", 2000, 0.0, 2.0);
  // pi-,proton1,proton2 + n
  new TH1F("pippn_Mass_CM","Invariant mass of #pi^{-}ppn;IM(#pi^{-}ppn) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pippn_Momentum_CM","Momentum of #pi^{-}ppn;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pippn_CosT_CM","cos#theta of #pi^{-}ppn;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("pippn_Phi_CM","#phi of #pi^{-}ppn;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("pippn_MMass_CM", "^{3}He(K^{-},#pi^{-}ppn)X missing mass;MM(#pi^{-}ppn) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pippn_MMass2_CM", "^{3}He(K^{-},#pi^{-}ppn)X missing mass^{2};MM^{2}(#pi^{-}ppn) (GeV/c^{2})^{2};Counts", 6000, -1.0, 2.0);
  new TH1F("pippn_MMomentum_CM", "^{3}He(K^{-},#pi^{-}ppn)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pippn_MCosT_CM", "^{3}He(K^{-},#pi^{-}ppn)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("pippn_MomentumTransfer_CM","Momentum transfer to #pi^{-}ppn;Momentum transfer (GeV/c);Counts", 2000, 0.0, 2.0);
  // pi-,proton1,proton2 (FWD)
  new TH1F("pipp_Mass_n_CM","Invariant mass of #pi^{-}pp;IM(#pi^{-}pp) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pipp_Momentum_n_CM","Momentum of #pi^{-}pp;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pipp_CosT_n_CM","cos#theta of #pi^{-}pp;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("pipp_Phi_n_CM","#phi of #pi^{-}pp;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("pipp_MMass_n_CM", "^{3}He(K^{-},#pi^{-}pp)X missing mass;MM(#pi^{-}pp) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pipp_MMomentum_n_CM", "^{3}He(K^{-},#pi^{-}pp)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pipp_MCosT_n_CM", "^{3}He(K^{-},#pi^{-}pp)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("pipp_MTheta_n_CM", "^{3}He(K^{-},#pi^{-}pp)X #theta;#theta;Counts", 200, -10.0, 190.0);
  new TH1F("pipp_MomentumTransfer_n_CM","Momentum transfer to #pi^{-}pp;Momentum transfer (GeV/c);Counts", 2000, 0.0, 2.0);
  /* 2D Plots */
  new TH2F("pip1_IMvspip2_IM_All", "(#pi^{-}p_{1}) invariant mass vs. (#pi^{-}p_{2}) invariant mass;IM(#pi^{-}p_{1}) (GeV/c^{2});IM(#pi^{-}p_{2}) (GeV/c^{2})", 2000, 1.0, 3.0, 2000, 1.0, 3.0);
  new TH2F("pip1_IMvspip2_IM_CM", "(#pi^{-}p_{1}) invariant mass vs. (#pi^{-}p_{2}) invariant mass;IM(#pi^{-}p_{1}) (GeV/c^{2});IM(#pi^{-}p_{2}) (GeV/c^{2})", 2000, 1.0, 3.0, 2000, 1.0, 3.0);
  new TH2F("pip1_IM2vspip2_IM2_All", "(#pi^{-}p_{1}) invariant mass vs. (#pi^{-}p_{2}) invariant mass;IM(#pi^{-}p_{1}) (GeV/c^{2});IM(#pi^{-}p_{2}) (GeV/c^{2})", 2000, 1.0, 5.0, 2000, 1.0, 5.0);
  new TH2F("pip1_IM2vspip2_IM2_CM", "(#pi^{-}p_{1}) invariant mass vs. (#pi^{-}p_{2}) invariant mass;IM(#pi^{-}p_{1}) (GeV/c^{2});IM(#pi^{-}p_{2}) (GeV/c^{2})", 2000, 1.0, 5.0, 2000, 1.0, 5.0);
  new TH2F("pipp_IMvsHeKpipp_MM_CM", "(#pi^{-}p_{1}p_{2}) invariant mass vs. ^{3}He(K^{-},#pi^{-}p_{1}p_{2})X missing mass;IM(#pi^{-}p_{1}p_{2}) (GeV/c^{2});MM(#pi^{-}p_{1}p_{2}) (GeV/c^{2})", 2000, 2.0, 4.0, 2000, 0.0, 2.0);
  new TH2F("pipp_IMvsHeKpipp_MCosT_CM", "(#pi^{-}p_{1}p_{2}) invariant mass vs. cos#theta_{n};IM(#pi^{-}p_{1}p_{2}) (GeV/c^{2});cos#theta_{n}", 2000, 2.0, 4.0, 3000, -1.5, 1.5);
  new TH2F("pipp_IMvsHeKpipp_MTheta_CM", "(#pi^{-}p_{1}p_{2}) invariant mass vs. #theta_{n};IM(#pi^{-}p_{1}p_{2}) (GeV/c^{2});#theta_{n}", 2000, 2.0, 4.0, 200, -10.5, 190.0);
  new TH2F("pipp_IMvsMomentumTransfer_CM", "(#pi^{-}p_{1}p_{2}) invariant mass vs. Momentum transfer;IM(#pi^{-}p_{1}p_{2}) (GeV/c^{2});Momentum Transfer (GeV/c)", 2000, 2.0, 4.0, 2000, 0.0, 2.0);
  new TH2F("pipp_DalitzPlot_CM", "Dalitz Plot;(T_{p}-T_{#Lambda})/#sqrt{3}Q;T_{n}/Q", 1000, -1.0, 1.0, 1000, -1.0, 1.0);
  /* 2D Plots  (FWD)*/
  new TH2F("pipp_IMvsHeKpippn_MM2_CM", "(#pi^{-}p_{1}p_{2}) invariant mass vs. ^{3}He(K^{-},#pi^{-}p_{1}p_{2}n)X missing mass^{2};IM(#pi^{-}p_{1}p_{2}) (GeV/c^{2});MM^{2}(#pi^{-}p_{1}p_{2}n) (GeV/c^{2})^{2}", 2000, 2.0, 4.0, 2000, -1.0, 1.0);
  new TH2F("pipp_IMvsHeKn_MM_CM", "(#pi^{-}p_{1}p_{2}) invariant mass vs. ^{3}He(K^{-},n)X missing mass^{2};IM(#pi^{-}p_{1}p_{2}) (GeV/c^{2});MM(n) (GeV/c^{2})", 2000, 2.0, 4.0, 2000, 2.0, 4.0);
  new TH2F("HeKpipp_MMomvsn_Mom_CM", "(#pi^{-}p_{1}p_{2}) momentum vs. n momentum;Missing Momentum (#pi^{-}p_{1}p_{2}) (GeV/c);Momentum (n) (GeV/c)", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("HeKpipp_MMvsHeKpippn_MM2_CM", "^{3}He(K^{-},(#pi^{-}p_{1}p_{2})X missing mass vs. ^{3}He(K^{-},#pi^{-}p_{1}p_{2}n)X missing mass^{2};MM(#pi^{-}p_{1}p_{2}) (GeV/c^{2});MM^{2}(#pi^{-}p_{1}p_{2}n) (GeV/c^{2})^{2}", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("HeKpn_MMvsHeKpippn_MM2_CM", "^{3}He(K^{-},(pn)X missing mass vs. ^{3}He(K^{-},#pi^{-}p_{1}p_{2}n)X missing mass^{2};MM(pn) (GeV/c^{2});MM^{2}(#pi^{-}p_{1}p_{2}n) (GeV/c^{2})^{2}", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("HeKpn_MM2vsHeKpippn_MM2_CM", "^{3}He(K^{-},(pn)X missing mass^{2} vs. ^{3}He(K^{-},#pi^{-}p_{1}p_{2}n)X missing mass^{2};MM^{2}(pn) (GeV/c^{2});MM^{2}(#pi^{-}p_{1}p_{2}n) (GeV/c^{2})^{2}", 4000, 0.0, 4.0, 2000, -1.0, 1.0);
  new TH2F("pipp_IMvsHeKpn_MM_CM", "(#pi^{-}p_{1}p_{2}) invariant mass vs. ^{3}He(K^{-},pn)X missing mass^{2};IM(#pi^{-}p_{1}p_{2}) (GeV/c^{2});MM^{2}(pn) (GeV/c^{2})^{2}", 2000, 2.0, 4.0, 1000, 1.0, 2.0);
  new TH2F("HeKn_MMvsHeKpn_MM_CM", "^{3}He(K^{-},n)X missing mass vs. ^{3}He(K^{-},pn)X missing mass^{2};MM(n) (GeV/c^{2});MM^{2}(pn) (GeV/c^{2})^{2}", 2000, 2.0, 4.0, 1000, 1.0, 2.0);
  new TH2F("pipp_IMvsHeKpipp_MM_n_CM", "(#pi^{-}p_{1}p_{2}) invariant mass vs. ^{3}He(K^{-},#pi^{-}p_{1}p_{2})X missing mass;IM(#pi^{-}p_{1}p_{2}) (GeV/c^{2});MM(#pi^{-}p_{1}p_{2}) (GeV/c^{2})", 2000, 2.0, 4.0, 2000, 0.0, 2.0);
  new TH2F("pipp_IMvsHeKpipp_MCosT_n_CM", "(#pi^{-}p_{1}p_{2}) invariant mass vs. cos#theta_{n};IM(#pi^{-}p_{1}p_{2}) (GeV/c^{2});cos#theta_{n}", 2000, 2.0, 4.0, 3000, -1.5, 1.5);
  new TH2F("pipp_IMvsMomentumTransfer_n_CM", "(#pi^{-}p_{1}p_{2}) invariant mass vs. Momentum transfer;IM(#pi^{-}p_{1}p_{2}) (GeV/c^{2});Momentum Transfer (GeV/c)", 2000, 2.0, 4.0, 2000, 0.0, 2.0);
  new TH2F("pipp_DalitzPlot_n_CM", "Dalitz Plot;(T_{p}-T_{#Lambda})/#sqrt{3}Q;T_{n}/Q", 1000, -1.0, 1.0, 1000, -1.0, 1.0);

  std::cout << "== Finish Histogram Initialization ==" << std::endl;
  return true;

}
