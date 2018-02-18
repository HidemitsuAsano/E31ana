// MyAnalysisHeKppin.cpp

#include "MyAnalysisHeKppin.h"

MyAnalysisHeKppin::MyAnalysisHeKppin(TFile* rt, ConfMan* conf)
{
  Initialize(conf);
  CutCondition();
  Clear();
}

MyAnalysisHeKppin::~MyAnalysisHeKppin()
{
  Clear();
  rtFile->cd();
  rtFile->Write();
  rtFile->Close();
}

void MyAnalysisHeKppin::Clear()
{
}

bool MyAnalysisHeKppin::DoAnalysis(ConfMan* conf, EventHeader* header, BeamLineHitMan* blMan, BeamLineTrackMan* bltrackMan, CDSHitMan* cdsMan, CDSTrackingMan* cdstrackMan, Particle* particle)
{
  rtFile->cd();

  FillHist("EventNumber",0); /* All Events */
  DetectorList *dlist=DetectorList::GetInstance();
  if(particle->nBeam()!=1) return false;
  pBeam* beam = particle->beam(0);
  double beamtof = beam->bhdt0tof();

  if(cdstrackMan->nGoodTrack()!=2) return true;
  FillHist("EventNumber",1); /* 2 Good Tracks in CDC */

  //bool chiflag = true;
  //for(int i=0; i<cdstrackMan->nGoodTrack(); i++){
  //  FillHist("CDC_Chi2",cdstrackMan->GoodTrack(i)->Chi());
  //  if(cdstrackMan->GoodTrack(i)->Chi()>30){
  //    chiflag = false;
  //  }
  //}
  //if(chiflag) return true;

  /* Event selection */
  int MulBVC=0;
  for(int i=0; i<blMan->nBVC(); i++){
    HodoscopeLikeHit* hit = blMan->BVC(i);
    if(hit->CheckRange()) MulBVC++;
  }
  int MulCVC=0;
  for(int i=0; i<blMan->nCVC(); i++){
    HodoscopeLikeHit* hit = blMan->CVC(i);
    if(hit->CheckRange()) MulCVC++;
  }
  int MulPC=0;
  for(int i=0; i<blMan->nPC(); i++){
    HodoscopeLikeHit* hit = blMan->PC(i);
    if(hit->CheckRange()) MulPC++;
  }
  int MulCDH=0;
  for(int i=0; i<cdsMan->nCDH(); i++){
    HodoscopeLikeHit* hit = cdsMan->CDH(i);
    if(hit->CheckRange()) MulCDH++;
  }
  if(MulBVC!=0||MulCVC!=0||MulPC!=0||MulCDH!=2) return true;

  if(particle->nCDS()!=2) return true;
  FillHist("EventNumber",2); /* CDH 2 hits and No FWD hit */

  if(particle->nProton()!=1) return true;
  if(particle->nPiminus()!=1&&particle->nPiplus()!=1) return true;
  FillHist("EventNumber",3); /* 1 protons and pion in CDS */

  int picharge = 0;
  if(particle->nPiminus()==1){ picharge = -1; }
  else if(particle->nPiplus()==1){ picharge = 1; }
  pCDS* pi = 0;
  if(picharge>0){
    pi = particle->pip(0);
  }
  else if(picharge<0){
    pi = particle->pim(0);
  }
  pCDS* p = particle->proton(0);

  // ########################################### //
  // Proton/Pion in FIDUCIAL and Vertex decision //
  // ########################################### //
  /* For Proton */
  TVector3 vertex_p = p->vbeam();
  double vdis_p = p->vbdis();
  if(GeomTools::GetID(vertex_p)!=CID_Fiducial){ return true; }
  /* For Pion */
  TVector3 vertex_pi = pi->vbeam();
  double vdis_pi = pi->vbdis();
  if(GeomTools::GetID(vertex_pi)!=CID_Fiducial){ return true; }
  TVector3 vertex;
  if(vdis_p>vdis_pi){ vertex = pi->vbeam(); }
  else{ vertex = p->vbeam(); }
  FillHist("EventNumber",4); /* Proton and Pion in FIDUCIAL */

  // ############## //
  // CDC Chi-square //
  // ############## //
  if(pi->chi()>30) return true;
  if(p->chi()>30) return true;
  FillHist("EventNumber",5); /* CDC Chi-square < 30 */

  // ######################### //
  // 2 Helix calculation       //
  // ######################### //
  pCDS* pip = 0;
  {
    for(int it=0; it<particle->nProduct(); it++){
      pCDS* product = particle->product(it);
      int comb = product->comb();
      if((comb==pow(2,CDS_Proton)+pow(2,CDS_PiMinus))||(comb==pow(2,CDS_Proton)+pow(2,CDS_PiPlus))){
        pip = product;
      }
    }
  }

  // ############### //
  // NC hit decision //
  // ############### //
  double time = 9999;
  int fnc = -1;
  for(int it=0; it<particle->nNC(); it++){
    pNC* nc = particle->nc(it);
    nc->CalcMom(beam,vertex);
    FillHist("FWD_OverbetavsMomentum",1.0/nc->beta(),nc->mom().Mag());
    FillHist("FWD_OverbetavsEnergy",1.0/nc->beta(),nc->energy());
    FillHist("FWD_TOFvsMomentum",nc->tof(),nc->mom().Mag());
    FillHist("FWD_HitPosition",nc->hitpos().X(),nc->hitpos().Y());
    double seg = (nc->seg()-1)%16+1;
    double lay = (nc->seg()-1)/16+1;
    FillHist("FWD_HitSegment",seg,lay);
    if(nc->pid()==F_Neutron) {
      if(nc->energy()>8.0 && (time>9998||time>nc->time())){
        time = nc->time();
        fnc = it;
      }
    }
  }
  pNC*  nc =0;
  if(fnc!=-1)  nc =particle->nc(fnc);
  if(nc==0) return true;
  FillHist("EventNumber",6); /* Neutron detection */

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
  if(picharge<0){
    vtxb = pi->vbeam();
    FillHist("pim_Vertex_XY",vtxb.X(),vtxb.Y());
    FillHist("pim_Vertex_ZX",vtxb.Z(),vtxb.X());
    FillHist("pim_Vertex_ZY",vtxb.Z(),vtxb.Y());
    FillHist("pim_VBDIS",pi->vbdis());
    FillHist("pim_OverbetavsMomentum",1.0/pi->beta(),pi->mom());
    FillHist("pim_Mass2vsMomentum",pi->mass2(),pi->mom());
  }
  /* Pi Plus */
  if(picharge>0){
    vtxb = pi->vbeam();
    FillHist("pip_Vertex_XY",vtxb.X(),vtxb.Y());
    FillHist("pip_Vertex_ZX",vtxb.Z(),vtxb.X());
    FillHist("pip_Vertex_ZY",vtxb.Z(),vtxb.Y());
    FillHist("pip_VBDIS",pi->vbdis());
    FillHist("pip_OverbetavsMomentum",1.0/pi->beta(),pi->mom());
    FillHist("pip_Mass2vsMomentum",pi->mass2(),pi->mom());
  }
  /* Proton */
  vtxb = p->vbeam();
  FillHist("p_Vertex_XY",vtxb.X(),vtxb.Y());
  FillHist("p_Vertex_ZX",vtxb.Z(),vtxb.X());
  FillHist("p_Vertex_ZY",vtxb.Z(),vtxb.Y());
  FillHist("p_VBDIS",p->vbdis());
  FillHist("p_OverbetavsMomentum",1.0/p->beta(),p->mom());
  FillHist("p_Mass2vsMomentum",p->mass2(),p->mom());
  /* Pi- and Proton */
  if(picharge<0){
    vtxb = pip->vbeam();
    FillHist("ppim_Vertex_XY",vtxb.X(),vtxb.Y());
    FillHist("ppim_Vertex_ZX",vtxb.Z(),vtxb.X());
    FillHist("ppim_Vertex_ZY",vtxb.Z(),vtxb.Y());
    FillHist("ppim_Mass",pip->vdis());
    FillHist("ppim_OA",pip->oa());
    FillHist("ppim_VDIS",pip->vdis());
    FillHist("ppim_VBDIS",pip->vbdis());
    FillHist("ppim_PBDCA",pip->pbdca());
  }
  /* Pi+ and Proton */
  if(picharge>0){
    vtxb = pip->vbeam();
    FillHist("ppip_Vertex_XY",vtxb.X(),vtxb.Y());
    FillHist("ppip_Vertex_ZX",vtxb.Z(),vtxb.X());
    FillHist("ppip_Vertex_ZY",vtxb.Z(),vtxb.Y());
    FillHist("ppip_Mass",pip->vdis());
    FillHist("ppip_OA",pip->oa());
    FillHist("ppip_VDIS",pip->vdis());
    FillHist("ppip_VBDIS",pip->vbdis());
    FillHist("ppip_PBDCA",pip->pbdca());
  }

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
  TLorentzVector Lpi = pi->GetLorentzVector();
  TLorentzVector cLpi = pi->GetLorentzVector(); cLpi.Boost(-boost);
  TLorentzVector Lp = p->GetLorentzVector();
  TLorentzVector cLp = p->GetLorentzVector(); cLp.Boost(-boost);
  //TLorentzVector Lpip = pip->GetLorentzVector();
  //TLorentzVector cLpip = pip->GetLorentzVector(); cLpip.Boost(-boost);
  TLorentzVector Lpip = Lpi + Lp;
  TLorentzVector cLpip = cLpi + cLp;
  TLorentzVector Llam = pip->GetLorentzVector();
  TLorentzVector cLlam = pip->GetLorentzVector(); cLpip.Boost(-boost);
  /* NC */
  TLorentzVector Ln  = nc ->GetLorentzVector();
  TLorentzVector cLn  = nc ->GetLorentzVector(); cLn.Boost(-boost);

  double pi_mass[2]  = { Lpi.M()                           , cLpi.M()};
  double pi_mom[2]   = { Lpi.Vect().Mag()                  , cLpi.Vect().Mag()};
  double pi_cost[2]  = { Lpi.Vect().CosTheta()             , cLpi.Vect().CosTheta()};
  double pi_phi[2]   = { Lpi.Vect().Phi()                  , cLpi.Vect().Phi()};
  double pi_mmass[2] = { (Ltgt+Lbeam-Lpi).M()              , (cLtgt+cLbeam-cLpi).M()};
  double pi_mmom[2]  = { (Ltgt+Lbeam-Lpi).Vect().Mag()     , (cLtgt+cLbeam-cLpi).Vect().Mag()};
  double pi_mcost[2] = { (Ltgt+Lbeam-Lpi).Vect().CosTheta(), (cLtgt+cLbeam-cLpi).Vect().CosTheta()};

  double p_mass[2]  = { Lp.M()                           , cLp.M()};
  double p_mom[2]   = { Lp.Vect().Mag()                  , cLp.Vect().Mag()};
  double p_cost[2]  = { Lp.Vect().CosTheta()             , cLp.Vect().CosTheta()};
  double p_phi[2]   = { Lp.Vect().Phi()                  , cLp.Vect().Phi()};
  double p_mmass[2] = { (Ltgt+Lbeam-Lp).M()              , (cLtgt+cLbeam-cLp).M()};
  double p_mmom[2]  = { (Ltgt+Lbeam-Lp).Vect().Mag()     , (cLtgt+cLbeam-cLp).Vect().Mag()};
  double p_mcost[2] = { (Ltgt+Lbeam-Lp).Vect().CosTheta(), (cLtgt+cLbeam-cLp).Vect().CosTheta()};

  double pip_mass[2]  = { Lpip.M()                           , cLpip.M()};
  double pip_mom[2]   = { Lpip.Vect().Mag()                  , cLpip.Vect().Mag()};
  double pip_cost[2]  = { Lpip.Vect().CosTheta()             , cLpip.Vect().CosTheta()};
  double pip_phi[2]   = { Lpip.Vect().Phi()                  , cLpip.Vect().Phi()};
  double pip_mmass[2] = { (Ltgt+Lbeam-Lpip).M()              , (cLtgt+cLbeam-cLpip).M()};
  double pip_mmom[2]  = { (Ltgt+Lbeam-Lpip).Vect().Mag()     , (cLtgt+cLbeam-cLpip).Vect().Mag()};
  double pip_mcost[2] = { (Ltgt+Lbeam-Lpip).Vect().CosTheta(), (cLtgt+cLbeam-cLpip).Vect().CosTheta()};

  double lam_mass[2]  = { Llam.M()                           , cLlam.M()};
  double lam_mom[2]   = { Llam.Vect().Mag()                  , cLlam.Vect().Mag()};
  double lam_cost[2]  = { Llam.Vect().CosTheta()             , cLlam.Vect().CosTheta()};
  double lam_phi[2]   = { Llam.Vect().Phi()                  , cLlam.Vect().Phi()};
  double lam_mmass[2] = { (Ltgt+Lbeam-Llam).M()              , (cLtgt+cLbeam-cLlam).M()};
  double lam_mmom[2]  = { (Ltgt+Lbeam-Llam).Vect().Mag()     , (cLtgt+cLbeam-cLlam).Vect().Mag()};
  double lam_mcost[2] = { (Ltgt+Lbeam-Llam).Vect().CosTheta(), (cLtgt+cLbeam-cLlam).Vect().CosTheta()};

  double nc_mom[2]     = { Ln.Vect().Mag()                  , cLn.Vect().Mag()};
  double nc_cost[2]    = { Ln.Vect().CosTheta()             , cLn.Vect().CosTheta()};
  double nc_phi[2]     = { Ln.Vect().Phi()                  , cLn.Vect().Phi()};
  double nc_mmass[2]   = { (Ltgt+Lbeam-Ln).M()              , (cLtgt+cLbeam-cLn).M()};
  double nc_mmom[2]    = { (Ltgt+Lbeam-Ln).Vect().Mag()     , (cLtgt+cLbeam-cLn).Vect().Mag()};
  double nc_mcost[2]   = { (Ltgt+Lbeam-Ln).Vect().CosTheta(), (cLtgt+cLbeam-cLn).Vect().CosTheta()};

  double pnc_mass[2]    = { (Lp+Ln).M()                           , (cLp+cLn).M()};
  double pnc_mom[2]     = { (Lp+Ln).Vect().Mag()                  , (cLp+cLn).Vect().Mag()};
  double pnc_cost[2]    = { (Lp+Ln).Vect().CosTheta()             , (cLp+cLn).Vect().CosTheta()};
  double pnc_phi[2]     = { (Lp+Ln).Vect().Phi()                  , (cLp+cLn).Vect().Phi()};
  double pnc_mmass[2]   = { (Ltgt+Lbeam-Lp-Ln).M()              , (cLtgt+cLbeam-cLp-cLn).M()};
  double pnc_mmass2[2]  = { (Ltgt+Lbeam-Lp-Ln).M2()              , (cLtgt+cLbeam-cLp-cLn).M2()};
  double pnc_mmom[2]    = { (Ltgt+Lbeam-Lp-Ln).Vect().Mag()     , (cLtgt+cLbeam-cLp-cLn).Vect().Mag()};
  double pnc_mcost[2]   = { (Ltgt+Lbeam-Lp-Ln).Vect().CosTheta(), (cLtgt+cLbeam-cLp-cLn).Vect().CosTheta()};

  double pinc_mass[2]    = { (Lpi+Ln).M()                           , (cLpi+cLn).M()};
  double pinc_mom[2]     = { (Lpi+Ln).Vect().Mag()                  , (cLpi+cLn).Vect().Mag()};
  double pinc_cost[2]    = { (Lpi+Ln).Vect().CosTheta()             , (cLpi+cLn).Vect().CosTheta()};
  double pinc_phi[2]     = { (Lpi+Ln).Vect().Phi()                  , (cLpi+cLn).Vect().Phi()};
  double pinc_mmass[2]   = { (Ltgt+Lbeam-Lpi-Ln).M()              , (cLtgt+cLbeam-cLpi-cLn).M()};
  double pinc_mmass2[2]  = { (Ltgt+Lbeam-Lpi-Ln).M2()              , (cLtgt+cLbeam-cLpi-cLn).M2()};
  double pinc_mmom[2]    = { (Ltgt+Lbeam-Lpi-Ln).Vect().Mag()     , (cLtgt+cLbeam-cLpi-cLn).Vect().Mag()};
  double pinc_mcost[2]   = { (Ltgt+Lbeam-Lpi-Ln).Vect().CosTheta(), (cLtgt+cLbeam-cLpi-cLn).Vect().CosTheta()};

  double pipnc_mass[2]    = { (Lpip+Ln).M()                           , (cLpip+cLn).M()};
  double pipnc_mom[2]     = { (Lpip+Ln).Vect().Mag()                  , (cLpip+cLn).Vect().Mag()};
  double pipnc_cost[2]    = { (Lpip+Ln).Vect().CosTheta()             , (cLpip+cLn).Vect().CosTheta()};
  double pipnc_phi[2]     = { (Lpip+Ln).Vect().Phi()                  , (cLpip+cLn).Vect().Phi()};
  double pipnc_mmass[2]   = { (Ltgt+Lbeam-Lpip-Ln).M()              , (cLtgt+cLbeam-cLpip-cLn).M()};
  double pipnc_mmass2[2]  = { (Ltgt+Lbeam-Lpip-Ln).M2()              , (cLtgt+cLbeam-cLpip-cLn).M2()};
  double pipnc_mmom[2]    = { (Ltgt+Lbeam-Lpip-Ln).Vect().Mag()     , (cLtgt+cLbeam-cLpip-cLn).Vect().Mag()};
  double pipnc_mcost[2]   = { (Ltgt+Lbeam-Lpip-Ln).Vect().CosTheta(), (cLtgt+cLbeam-cLpip-cLn).Vect().CosTheta()};

  /* flag */
  bool lamflag = false; /* lam flag */
  bool sbllamflag = false; /* lower sideband of lam flag */
  bool sbulamflag = false; /* upper sideband of lam flag */
  if(picharge<0){
    double pro_mass = lam_mass[1];
    if(lamll<pro_mass&&pro_mass<lamul) lamflag = true;
  }
  bool smflag = false; /* sm flag */
  bool sblsmflag = false; /* lower sideband of sm flag */
  bool sbusmflag = false; /* upper sideband of sm flag */
  if(picharge<0){
    double pro_mass = pinc_mass[1];
    if(smll<pro_mass&&pro_mass<smul) smflag = true;
  }
  bool spflag = false; /* sp flag */
  bool sblspflag = false; /* lower sideband of sp flag */
  bool sbuspflag = false; /* upper sideband of sp flag */
  if(picharge>0){
    double pro_mass = pinc_mass[1];
    if(spll<pro_mass&&pro_mass<spul) spflag = true;
  }

  if(spflag){ return true; }
  if(smflag){ return true; }
  if(lamflag){ return true; }
  //if(!mnflag){ return true; }
  //FillHist("EventNumber",10); /* Missing Neutron */

  /*-----------------------------------------*/
  /* All                                     */ 
  /*-----------------------------------------*/
  // ######################### //
  // Invariant/Missing Mass    //
  // ######################### //
  // pi-
  if(picharge<0){
    for(int f=1; f<2; f++){
      FillHist(Form("pim_Momentum_%s",frame[f].Data()),pi_mom[f]);
      FillHist(Form("pim_CosT_%s",frame[f].Data()),pi_cost[f]);
      FillHist(Form("pim_Phi_%s",frame[f].Data()),pi_phi[f]);
      FillHist(Form("pim_MMass_%s",frame[f].Data()),pi_mmass[f]);
      FillHist(Form("pim_MMomentum_%s",frame[f].Data()),pi_mmom[f]);
      FillHist(Form("pim_MCosT_%s",frame[f].Data()),pi_mcost[f]);
    }
  }
  // pi+
  if(picharge>0){
    for(int f=1; f<2; f++){
      FillHist(Form("pip_Momentum_%s",frame[f].Data()),pi_mom[f]);
      FillHist(Form("pip_CosT_%s",frame[f].Data()),pi_cost[f]);
      FillHist(Form("pip_Phi_%s",frame[f].Data()),pi_phi[f]);
      FillHist(Form("pip_MMass_%s",frame[f].Data()),pi_mmass[f]);
      FillHist(Form("pip_MMomentum_%s",frame[f].Data()),pi_mmom[f]);
      FillHist(Form("pip_MCosT_%s",frame[f].Data()),pi_mcost[f]);
    }
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
  // pi- proton
  if(picharge<0){
    for(int f=1; f<2; f++){
      FillHist(Form("ppim_Mass_%s",frame[f].Data()),pip_mass[f]);
      FillHist(Form("ppim_Momentum_%s",frame[f].Data()),pip_mom[f]);
      FillHist(Form("ppim_CosT_%s",frame[f].Data()),pip_cost[f]);
      FillHist(Form("ppim_Phi_%s",frame[f].Data()),pip_phi[f]);
      FillHist(Form("ppim_MMass_%s",frame[f].Data()),pip_mmass[f]);
      FillHist(Form("ppim_MMomentum_%s",frame[f].Data()),pip_mmom[f]);
      FillHist(Form("ppim_MCosT_%s",frame[f].Data()),pip_mcost[f]);
    }
  }
  // pi+ proton
  if(picharge>0){
    for(int f=1; f<2; f++){
      FillHist(Form("ppip_Mass_%s",frame[f].Data()),pip_mass[f]);
      FillHist(Form("ppip_Momentum_%s",frame[f].Data()),pip_mom[f]);
      FillHist(Form("ppip_CosT_%s",frame[f].Data()),pip_cost[f]);
      FillHist(Form("ppip_Phi_%s",frame[f].Data()),pip_phi[f]);
      FillHist(Form("ppip_MMass_%s",frame[f].Data()),pip_mmass[f]);
      FillHist(Form("ppip_MMomentum_%s",frame[f].Data()),pip_mmom[f]);
      FillHist(Form("ppip_MCosT_%s",frame[f].Data()),pip_mcost[f]);
    }
  }

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
  // n + pi-
  if(picharge<0){
    for(int f=1; f<2; f++){
      FillHist(Form("pimn_Mass_%s",frame[f].Data()),pinc_mass[f]);
      FillHist(Form("pimn_Momentum_%s",frame[f].Data()),pinc_mom[f]);
      FillHist(Form("pimn_CosT_%s",frame[f].Data()),pinc_cost[f]);
      FillHist(Form("pimn_Phi_%s",frame[f].Data()),pinc_phi[f]);
      FillHist(Form("pimn_MMass_%s",frame[f].Data()),pinc_mmass[f]);
      FillHist(Form("pimn_MMass2_%s",frame[f].Data()),pinc_mmass2[f]);
      FillHist(Form("pimn_MMomentum_%s",frame[f].Data()),pinc_mmom[f]);
      FillHist(Form("pimn_MCosT_%s",frame[f].Data()),pinc_mcost[f]);
    }
  }
  // n + pi+
  if(picharge>0){
    for(int f=1; f<2; f++){
      FillHist(Form("pipn_Mass_%s",frame[f].Data()),pinc_mass[f]);
      FillHist(Form("pipn_Momentum_%s",frame[f].Data()),pinc_mom[f]);
      FillHist(Form("pipn_CosT_%s",frame[f].Data()),pinc_cost[f]);
      FillHist(Form("pipn_Phi_%s",frame[f].Data()),pinc_phi[f]);
      FillHist(Form("pipn_MMass_%s",frame[f].Data()),pinc_mmass[f]);
      FillHist(Form("pipn_MMass2_%s",frame[f].Data()),pinc_mmass2[f]);
      FillHist(Form("pipn_MMomentum_%s",frame[f].Data()),pinc_mmom[f]);
      FillHist(Form("pipn_MCosT_%s",frame[f].Data()),pinc_mcost[f]);
    }
  }
  // n + pi- + p
  if(picharge<0){
    for(int f=1; f<2; f++){
      FillHist(Form("pimpn_Mass_%s",frame[f].Data()),pipnc_mass[f]);
      FillHist(Form("pimpn_Momentum_%s",frame[f].Data()),pipnc_mom[f]);
      FillHist(Form("pimpn_CosT_%s",frame[f].Data()),pipnc_cost[f]);
      FillHist(Form("pimpn_Phi_%s",frame[f].Data()),pipnc_phi[f]);
      FillHist(Form("pimpn_MMass_%s",frame[f].Data()),pipnc_mmass[f]);
      FillHist(Form("pimpn_MMass2_%s",frame[f].Data()),pipnc_mmass2[f]);
      FillHist(Form("pimpn_MMomentum_%s",frame[f].Data()),pipnc_mmom[f]);
      FillHist(Form("pimpn_MCosT_%s",frame[f].Data()),pipnc_mcost[f]);
    }
  }
  // n + pi+ + p
  if(picharge>0){
    for(int f=1; f<2; f++){
      FillHist(Form("pippn_Mass_%s",frame[f].Data()),pipnc_mass[f]);
      FillHist(Form("pippn_Momentum_%s",frame[f].Data()),pipnc_mom[f]);
      FillHist(Form("pippn_CosT_%s",frame[f].Data()),pipnc_cost[f]);
      FillHist(Form("pippn_Phi_%s",frame[f].Data()),pipnc_phi[f]);
      FillHist(Form("pippn_MMass_%s",frame[f].Data()),pipnc_mmass[f]);
      FillHist(Form("pippn_MMass2_%s",frame[f].Data()),pipnc_mmass2[f]);
      FillHist(Form("pippn_MMomentum_%s",frame[f].Data()),pipnc_mmom[f]);
      FillHist(Form("pippn_MCosT_%s",frame[f].Data()),pipnc_mcost[f]);
    }
  }

  return true;

}

bool MyAnalysisHeKppin::FillHist(TString name, double val1, int weight)
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

bool MyAnalysisHeKppin::FillHist(TString name, TString val1, int weight)
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

bool MyAnalysisHeKppin::FillHist(TString name, double val1, double val2, int weight)
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

bool MyAnalysisHeKppin::FillHist(TString name, TString val1, TString val2, int weight)
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

void MyAnalysisHeKppin::CutCondition()
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

bool MyAnalysisHeKppin::Initialize(ConfMan* confMan)
{
  std::cout << "### MyAnalysisHeKppin::Initialize ###" << std::endl;

  std::string ofname = confMan->GetOutFileName();
  ofname.insert(ofname.find(".root"),"_anaHeKppin");

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
  new TH1F( "pip_VBDIS", "VBDIS of #pi^{+};VBDIS (cm);Coutns", 100, 0, 10);
  new TH1F( "p_VBDIS", "VBDIS of p;VBDIS (cm);Coutns", 100, 0, 10);
  new TH1F( "ppim_VDIS", "VDIS of #pi^{-}p};VDIS (cm);Coutns", 100, 0, 10);
  new TH1F( "ppip_VDIS", "VDIS of #pi^{+}p;VDIS (cm);Coutns", 100, 0, 10);
  new TH1F( "ppim_VBDIS", "VBDIS of #pi^{-}p};VBDIS (cm);Coutns", 100, 0, 10);
  new TH1F( "ppip_VBDIS", "VBDIS of #pi^{+}p;VBDIS (cm);Coutns", 100, 0, 10);
  new TH1F( "ppim_PBDCA", "PBDCA of #pi^{-}p;PBDCA (cm);Coutns", 100, 0, 10);
  new TH1F( "ppip_PBDCA", "PBDCA of #pi^{+}p;PBDCA (cm);Coutns", 100, 0, 10);
  new TH2F( "pim_Vertex_XY", "Vertex XY plane #pi^{-};X (cm);Y (cm)", 300, -15, 15, 300, -15, 15);
  new TH2F( "pim_Vertex_ZX", "Vertex ZX plane #pi^{-};Z (cm);X (cm)", 600, -30, 30, 300, -15, 15);
  new TH2F( "pim_Vertex_ZY", "Vertex ZY plane #pi^{-};Z (cm);Y (cm)", 600, -30, 30, 300, -15, 15);
  new TH2F( "pip_Vertex_XY", "Vertex XY plane #pi^{+};X (cm);Y (cm)", 300, -15, 15, 300, -15, 15);
  new TH2F( "pip_Vertex_ZX", "Vertex ZX plane #pi^{+};Z (cm);X (cm)", 600, -30, 30, 300, -15, 15);
  new TH2F( "pip_Vertex_ZY", "Vertex ZY plane #pi^{+};Z (cm);Y (cm)", 600, -30, 30, 300, -15, 15);
  new TH2F( "p_Vertex_XY", "Vertex XY plane p;X (cm);Y (cm)", 300, -15, 15, 300, -15, 15);
  new TH2F( "p_Vertex_ZX", "Vertex ZX plane p;Z (cm);X (cm)", 600, -30, 30, 300, -15, 15);
  new TH2F( "p_Vertex_ZY", "Vertex ZY plane p;Z (cm);Y (cm)", 600, -30, 30, 300, -15, 15);
  new TH2F( "ppim_Vertex_XY", "Vertex XY plane #pi^{-}p;X (cm);Y (cm)", 300, -15, 15, 300, -15, 15);
  new TH2F( "ppim_Vertex_ZX", "Vertex ZX plane #pi^{-}p;Z (cm);X (cm)", 600, -30, 30, 300, -15, 15);
  new TH2F( "ppim_Vertex_ZY", "Vertex ZY plane #pi^{-}p;Z (cm);Y (cm)", 600, -30, 30, 300, -15, 15);
  new TH2F( "ppip_Vertex_XY", "Vertex XY plane #pi^{+}p;X (cm);Y (cm)", 300, -15, 15, 300, -15, 15);
  new TH2F( "ppip_Vertex_ZX", "Vertex ZX plane #pi^{+}p;Z (cm);X (cm)", 600, -30, 30, 300, -15, 15);
  new TH2F( "ppip_Vertex_ZY", "Vertex ZY plane #pi^{+}p;Z (cm);Y (cm)", 600, -30, 30, 300, -15, 15);

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
  // pi+
  new TH1F("pip_Momentum_CM","Momentum of #pi^{+};Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pip_CosT_CM","cos#theta of #pi^{+};cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("pip_Phi_CM","#phi of #pi^{+};#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("pip_MMass_CM", "^{3}He(K^{-},#pi^{+})X missing mass;MM(#pi^{+}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pip_MMomentum_CM", "^{3}He(K^{-},#pi^{+})X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pip_MCosT_CM", "^{3}He(K^{-},#pi^{+})X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // proton
  new TH1F("p_Momentum_CM","Momentum of p;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("p_CosT_CM","cos#theta of p;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("p_Phi_CM","#phi of p;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("p_MMass_CM", "^{3}He(K^{-},p)X missing mass;MM(p) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("p_MMomentum_CM", "^{3}He(K^{-},p)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("p_MCosT_CM", "^{3}He(K^{-},p)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // pi-,proton
  new TH1F("ppim_Mass_CM","Invariant mass of #pi^{-}p;IM(#pi^{-}p) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("ppim_Momentum_CM","Momentum of #pi^{-}p;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("ppim_CosT_CM","cos#theta of #pi^{-}p;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("ppim_Phi_CM","#phi of #pi^{-}p;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("ppim_MMass_CM", "^{3}He(K^{-},#pi^{-}p)X missing mass;MM(#pi^{-}p) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("ppim_MMomentum_CM", "^{3}He(K^{-},#pi^{-}p)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("ppim_MCosT_CM", "^{3}He(K^{-},#pi^{-}p)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // pi+,proton
  new TH1F("ppip_Mass_CM","Invariant mass of #pi^{+}p;IM(#pi^{+}p) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("ppip_Momentum_CM","Momentum of #pi^{+}p;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("ppip_CosT_CM","cos#theta of #pi^{+}p;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("ppip_Phi_CM","#phi of #pi^{+}p;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("ppip_MMass_CM", "^{3}He(K^{-},#pi^{+}p)X missing mass;MM(#pi^{+}p) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("ppip_MMomentum_CM", "^{3}He(K^{-},#pi^{+}p)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("ppip_MCosT_CM", "^{3}He(K^{-},#pi^{+}p)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
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
  // pi- + n
  new TH1F("pimn_Mass_CM","Invariant mass of #pi^{-}n;IM(#pi^{-}n) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pimn_Momentum_CM","Momentum of #pi^{-}n;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pimn_CosT_CM","cos#theta of #pi^{-}n;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("pimn_Phi_CM","#phi of #pi^{-}n;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("pimn_MMass_CM", "^{3}He(K^{-},#pi^{-}n)X missing mass;MM(#pi^{-}n) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pimn_MMass2_CM", "^{3}He(K^{-},#pi^{-}n)X missing mass^{2};MM^{2}(#pi^{-}n) (GeV/c^{2})^{2};Counts", 6000, -1.0, 2.0);
  new TH1F("pimn_MMomentum_CM", "^{3}He(K^{-},#pi^{-}n)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pimn_MCosT_CM", "^{3}He(K^{-},#pi^{-}n)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("pimn_MomentumTransfer_CM","Momentum transfer to #pi^{-}n;Momentum transfer (GeV/c);Counts", 2000, 0.0, 2.0);
  // pi+ + n
  new TH1F("pipn_Mass_CM","Invariant mass of #pi^{+}n;IM(#pi^{+}n) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pipn_Momentum_CM","Momentum of #pi^{+}n;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pipn_CosT_CM","cos#theta of #pi^{+}n;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("pipn_Phi_CM","#phi of #pi^{+}n;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("pipn_MMass_CM", "^{3}He(K^{-},#pi^{+}n)X missing mass;MM(#pi^{+}n) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pipn_MMass2_CM", "^{3}He(K^{-},#pi^{+}n)X missing mass^{2};MM^{2}(#pi^{+}n) (GeV/c^{2})^{2};Counts", 6000, -1.0, 2.0);
  new TH1F("pipn_MMomentum_CM", "^{3}He(K^{-},#pi^{+}n)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pipn_MCosT_CM", "^{3}He(K^{-},#pi^{+}n)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("pipn_MomentumTransfer_CM","Momentum transfer to #pi^{+}n;Momentum transfer (GeV/c);Counts", 2000, 0.0, 2.0);
  // pi- + p + n
  new TH1F("pimpn_Mass_CM","Invariant mass of #pi^{-}pn;IM(#pi^{-}pn) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pimpn_Momentum_CM","Momentum of #pi^{-}pn;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pimpn_CosT_CM","cos#theta of #pi^{-}pn;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("pimpn_Phi_CM","#phi of #pi^{-}pn;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("pimpn_MMass_CM", "^{3}He(K^{-},#pi^{-}pn)X missing mass;MM(#pi^{-}pn) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pimpn_MMass2_CM", "^{3}He(K^{-},#pi^{-}pn)X missing mass^{2};MM^{2}(#pi^{-}pn) (GeV/c^{2})^{2};Counts", 6000, -1.0, 2.0);
  new TH1F("pimpn_MMomentum_CM", "^{3}He(K^{-},#pi^{-}pn)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pimpn_MCosT_CM", "^{3}He(K^{-},#pi^{-}pn)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("pimpn_MomentumTransfer_CM","Momentum transfer to #pi^{-}pn;Momentum transfer (GeV/c);Counts", 2000, 0.0, 2.0);
  // pi+ + p + n
  new TH1F("pippn_Mass_CM","Invariant mass of #pi^{+}pn;IM(#pi^{+}pn) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pippn_Momentum_CM","Momentum of #pi^{+}pn;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pippn_CosT_CM","cos#theta of #pi^{+}pn;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("pippn_Phi_CM","#phi of #pi^{+}pn;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("pippn_MMass_CM", "^{3}He(K^{-},#pi^{+}pn)X missing mass;MM(#pi^{+}pn) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pippn_MMass2_CM", "^{3}He(K^{-},#pi^{+}pn)X missing mass^{2};MM^{2}(#pi^{+}pn) (GeV/c^{2})^{2};Counts", 6000, -1.0, 2.0);
  new TH1F("pippn_MMomentum_CM", "^{3}He(K^{-},#pi^{+}pn)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pippn_MCosT_CM", "^{3}He(K^{-},#pi^{+}pn)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("pippn_MomentumTransfer_CM","Momentum transfer to #pi^{+}pn;Momentum transfer (GeV/c);Counts", 2000, 0.0, 2.0);

  std::cout << "== Finish Histogram Initialization ==" << std::endl;
  return true;

}
