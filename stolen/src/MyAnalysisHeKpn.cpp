// MyAnalysisHeKpn.cpp

#include "MyAnalysisHeKpn.h"

MyAnalysisHeKpn::MyAnalysisHeKpn(TFile* rt, ConfMan* conf)
{
  Initialize(conf);
  CutCondition();
  Clear();
}

MyAnalysisHeKpn::~MyAnalysisHeKpn()
{
  Clear();
  rtFile->cd();
  rtFile->Write();
  rtFile->Close();
}

void MyAnalysisHeKpn::Clear()
{
}

bool MyAnalysisHeKpn::DoAnalysis(ConfMan* conf, EventHeader* header, BeamLineHitMan* blMan, BeamLineTrackMan* bltrackMan, CDSHitMan* cdsMan, CDSTrackingMan* cdstrackMan, Particle* particle)
{
  DetectorList *dlist=DetectorList::GetInstance();
  int ievent=0;
  rtFile->cd();

  FillHist("EventNumber",ievent); ievent++; /* All Events */

  if(particle->nBeam()!=1) return false;
  FillHist("EventNumber",ievent); ievent++; /* Beam Analysis */
  pBeam* beam = particle->beam(0);
  double beamtof = beam->bhdt0tof();

  if(!header->IsTrig(Trig_Neutral)){ return true; }
  FillHist("EventNumber",ievent); ievent++; /* Neutral Trigger */

  int MulCDH=0;
  for(int i=0; i<cdsMan->nCDH(); i++){
    HodoscopeLikeHit* hit = cdsMan->CDH(i);
    if(hit->CheckRange()) MulCDH++;
  }
  FillHist("CDH_Multiplicity",MulCDH);
  if(MulCDH!=3){ return true; }
  FillHist("EventNumber",ievent); ievent++; /* 3 hits in CDH */

  if(cdstrackMan->nGoodTrack()==0) return true;
  FillHist("EventNumber",ievent); ievent++; /* Least 1 Good Tracks in CDC */

  /* Event selection */
  int MulBVC=0;
  for(int i=0; i<blMan->nBVC(); i++){
    HodoscopeLikeHit* hit = blMan->BVC(i);
    if(hit->CheckRange()) MulBVC++;
  }
  FillHist("BVC_Multiplicity",MulBVC);
  int MulCVC=0;
  for(int i=0; i<blMan->nCVC(); i++){
    HodoscopeLikeHit* hit = blMan->CVC(i);
    if(hit->CheckRange()) MulCVC++;
  }
  FillHist("CVC_Multiplicity",MulCVC);
  int MulPC=0;
  for(int i=0; i<blMan->nPC(); i++){
    HodoscopeLikeHit* hit = blMan->PC(i);
    if(hit->CheckRange()) MulPC++;
  }
  FillHist("PC_Multiplicity",MulPC);
  if(MulBVC!=0||MulCVC!=0||MulPC!=0) return true;

  int MulIH=0;
  for(int i=0; i<cdsMan->nIH(); i++){
    HodoscopeLikeHit* hit = cdsMan->IH(i);
    if(hit->CheckRange()) MulIH++;
  }
  FillHist("IH_Multiplicity",MulIH);
  //if(MulIH!=3){ return true; }
  FillHist("EventNumber",ievent); ievent++; /* 3 hits in IH */

  if(particle->nCDS()==0) return true;
  FillHist("EventNumber",ievent); ievent++; /* No FWD charged hit */

  if(particle->nProton()==0) return true;
  if(particle->nProton()>2) return true;
  //if(particle->nPiminus()!=1) return true;
  FillHist("EventNumber",ievent); ievent++; /* 2 protons and pion in CDS */

  // ############## //
  // p1/p2 decision //
  // ############## //
  TVector3 vertex;
  double vdis = 9999;
  bool FIDUCIAL = false;
  int fcds = -1;
  for(int it=0; it<particle->nProton(); it++){
    pCDS* cds = particle->proton(it);
    if(vdis>9998||vdis>cds->vdis()){
      vdis = cds->vdis();
      vertex = cds->vbeam();
      fcds=it;
    }
  }
  if(fcds<0){ return true; }

  pCDS* pim = particle->pim(0);
  pCDS* p = particle->proton(fcds);

  // ############## //
  // CDC Chi-square //
  // ############## //
  if(p->chi()>30) return true;
  FillHist("EventNumber",ievent); ievent++; /* CDC Chi-square < 30 */

  // ################ //
  // Vertex decision1 //
  // ################ //
  vertex = p->vbeam();
  vdis = p->vbdis();
  if(GeomTools::GetID(vertex)!=CID_Fiducial){ return true; }
  FillHist("EventNumber",ievent); ievent++; /* Vertex Decision 1 */

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
  /* Proton */
  vtxb = p->vbeam();
  FillHist("p_Vertex_XY",vtxb.X(),vtxb.Y());
  FillHist("p_Vertex_ZX",vtxb.Z(),vtxb.X());
  FillHist("p_Vertex_ZY",vtxb.Z(),vtxb.Y());
  FillHist("p_VBDIS",p->vbdis());

  // ############ //
  // DCA Decision //
  // ############ //
  //if(proton->pbdca()>1.0){ return true; }
  //if(lam->pbdca()>1.0){ return true; }
  //FillHist("EventNumber",ievent); ievent++; /* DCA Decision */

  // Trigger Pattern
  FillHist("TriggerPattern",0);
  for( int i=0; i<20; i++ ){
    int val = header->pattern(i);
    if( 0<val ){
      FillHist("TriggerPattern",i);
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

  TString frame[2] = {"Lab","CM"};
  TString pname[8] = {"pip","p","d","t","he","pim","k","o"};

  // Target
  TVector3 ZeroV;
  TLorentzVector Ltgt; Ltgt.SetVectM(ZeroV,ThreeHeMass);
  TLorentzVector Lbeam = beam->GetLorentzVector(vertex);
  TVector3 boost = (Ltgt+Lbeam).BoostVector();
  TLorentzVector cLtgt = Ltgt; cLtgt.Boost(-boost); TLorentzVector cLbeam = Lbeam; cLbeam.Boost(-boost);

  // proton, neutron //
  TLorentzVector Lp = p->GetLorentzVector();
  TLorentzVector cLp = p->GetLorentzVector(); cLp.Boost(-boost);
  
  TLorentzVector Ln  = nc ->GetLorentzVector();
  TLorentzVector cLn  = nc ->GetLorentzVector(); cLn.Boost(-boost);

  double p_mass[2]  = { Lp.M()                           , cLp.M()};
  double p_mom[2]   = { Lp.Vect().Mag()                  , cLp.Vect().Mag()};
  double p_cost[2]  = { Lp.Vect().CosTheta()             , cLp.Vect().CosTheta()};
  double p_phi[2]   = { Lp.Vect().Phi()                  , cLp.Vect().Phi()};
  double p_mmass[2] = { (Ltgt+Lbeam-Lp).M()              , (cLtgt+cLbeam-cLp).M()};
  double p_mmom[2]  = { (Ltgt+Lbeam-Lp).Vect().Mag()     , (cLtgt+cLbeam-cLp).Vect().Mag()};
  double p_mcost[2] = { (Ltgt+Lbeam-Lp).Vect().CosTheta(), (cLtgt+cLbeam-cLp).Vect().CosTheta()};

  double n_mom[2]     = { Ln.Vect().Mag()                  , cLn.Vect().Mag()};
  double n_cost[2]    = { Ln.Vect().CosTheta()             , cLn.Vect().CosTheta()};
  double n_phi[2]     = { Ln.Vect().Phi()                  , cLn.Vect().Phi()};
  double n_mmass[2]   = { (Ltgt+Lbeam-Ln).M()              , (cLtgt+cLbeam-cLn).M()};
  double n_mmom[2]    = { (Ltgt+Lbeam-Ln).Vect().Mag()     , (cLtgt+cLbeam-cLn).Vect().Mag()};
  double n_mcost[2]   = { (Ltgt+Lbeam-Ln).Vect().CosTheta(), (cLtgt+cLbeam-cLn).Vect().CosTheta()};

  double pn_mom[2]     = { (Lp+Ln).Vect().Mag()                  , (cLp+cLn).Vect().Mag()};
  double pn_cost[2]    = { (Lp+Ln).Vect().CosTheta()             , (cLp+cLn).Vect().CosTheta()};
  double pn_phi[2]     = { (Lp+Ln).Vect().Phi()                  , (cLp+cLn).Vect().Phi()};
  double pn_mmass[2]   = { (Ltgt+Lbeam-Lp-Ln).M()              , (cLtgt+cLbeam-cLp-cLn).M()};
  double pn_mmom[2]    = { (Ltgt+Lbeam-Lp-Ln).Vect().Mag()     , (cLtgt+cLbeam-cLp-cLn).Vect().Mag()};
  double pn_mcost[2]   = { (Ltgt+Lbeam-Lp-Ln).Vect().CosTheta(), (cLtgt+cLbeam-cLp-cLn).Vect().CosTheta()};

  /*-----------------------------------------*/
  /* All                                     */ 
  /*-----------------------------------------*/
  // ######################### //
  // Invariant/Missing Mass    //
  // ######################### //
  // proton
  for(int f=1; f<2; f++){
    FillHist(Form("p_Momentum_%s",frame[f].Data()),p_mom[f]);
    FillHist(Form("p_CosT_%s",frame[f].Data()),p_cost[f]);
    FillHist(Form("p_Phi_%s",frame[f].Data()),p_phi[f]);
    FillHist(Form("p_MMass_%s",frame[f].Data()),p_mmass[f]);
    FillHist(Form("p_MMomentum_%s",frame[f].Data()),p_mmom[f]);
    FillHist(Form("p_MCosT_%s",frame[f].Data()),p_mcost[f]);
  }
  // neutron
  for(int f=1; f<2; f++){
    FillHist(Form("n_Momentum_%s",frame[f].Data()),n_mom[f]);
    FillHist(Form("n_CosT_%s",frame[f].Data()),n_cost[f]);
    FillHist(Form("n_Phi_%s",frame[f].Data()),n_phi[f]);
    FillHist(Form("n_MMass_%s",frame[f].Data()),n_mmass[f]);
    FillHist(Form("n_MMomentum_%s",frame[f].Data()),n_mmom[f]);
    FillHist(Form("n_MCosT_%s",frame[f].Data()),n_mcost[f]);
  }
  // proton + neutron
  for(int f=1; f<2; f++){
    FillHist(Form("pn_Momentum_%s",frame[f].Data()),pn_mom[f]);
    FillHist(Form("pn_CosT_%s",frame[f].Data()),pn_cost[f]);
    FillHist(Form("pn_Phi_%s",frame[f].Data()),pn_phi[f]);
    FillHist(Form("pn_MMass_%s",frame[f].Data()),pn_mmass[f]);
    FillHist(Form("pn_MMomentum_%s",frame[f].Data()),pn_mmom[f]);
    FillHist(Form("pn_MCosT_%s",frame[f].Data()),pn_mcost[f]);
  }
  /*-----------------------------------------*/
  /* 2D Plots                                */ 
  /*-----------------------------------------*/
  // MM(n) vs. MM(pn)
  for(int f=1; f<2; f++){
    FillHist(Form("n_MMvspn_MM_%s",frame[f].Data()),n_mmass[f],pn_mmass[f]);
  }


  /* flag */
  bool pflag = false; /* p flag */
  {
    for(int it=0; it<particle->nProton(); it++){
      if(it==fcds) continue;

      pCDS* cds = particle->proton(it);
      if(cds->chi()<30) pflag=true;
    }
  }
  bool pimflag = false; /* pi- flag */
  {
    for(int it=0; it<particle->nPiminus(); it++){
      pCDS* cds = particle->pim(it);
      if(cds->chi()<30) pimflag=true;
    }
  }
  bool pipflag = false; /* pi+ flag */
  {
    for(int it=0; it<particle->nPiplus(); it++){
      pCDS* cds = particle->pip(it);
      if(cds->chi()<30) pipflag=true;
    }
  }

  if(pflag){
    // proton
    for(int f=1; f<2; f++){
      FillHist(Form("p_Momentum_withp_%s",frame[f].Data()),p_mom[f]);
      FillHist(Form("p_CosT_withp_%s",frame[f].Data()),p_cost[f]);
      FillHist(Form("p_Phi_withp_%s",frame[f].Data()),p_phi[f]);
      FillHist(Form("p_MMass_withp_%s",frame[f].Data()),p_mmass[f]);
      FillHist(Form("p_MMomentum_withp_%s",frame[f].Data()),p_mmom[f]);
      FillHist(Form("p_MCosT_withp_%s",frame[f].Data()),p_mcost[f]);
    }
    // neutron
    for(int f=1; f<2; f++){
      FillHist(Form("n_Momentum_withp_%s",frame[f].Data()),n_mom[f]);
      FillHist(Form("n_CosT_withp_%s",frame[f].Data()),n_cost[f]);
      FillHist(Form("n_Phi_withp_%s",frame[f].Data()),n_phi[f]);
      FillHist(Form("n_MMass_withp_%s",frame[f].Data()),n_mmass[f]);
      FillHist(Form("n_MMomentum_withp_%s",frame[f].Data()),n_mmom[f]);
      FillHist(Form("n_MCosT_withp_%s",frame[f].Data()),n_mcost[f]);
    }
    // proton + neutron
    for(int f=1; f<2; f++){
      FillHist(Form("pn_Momentum_withp_%s",frame[f].Data()),pn_mom[f]);
      FillHist(Form("pn_CosT_withp_%s",frame[f].Data()),pn_cost[f]);
      FillHist(Form("pn_Phi_withp_%s",frame[f].Data()),pn_phi[f]);
      FillHist(Form("pn_MMass_withp_%s",frame[f].Data()),pn_mmass[f]);
      FillHist(Form("pn_MMomentum_withp_%s",frame[f].Data()),pn_mmom[f]);
      FillHist(Form("pn_MCosT_withp_%s",frame[f].Data()),pn_mcost[f]);
    }
    /*-----------------------------------------*/
    /* 2D Plots                                */ 
    /*-----------------------------------------*/
    // MM(n) vs. MM(pn)
    for(int f=1; f<2; f++){
      FillHist(Form("n_MMvspn_MM_withp_%s",frame[f].Data()),n_mmass[f],pn_mmass[f]);
    }
  }
  if(pimflag){
    // proton
    for(int f=1; f<2; f++){
      FillHist(Form("p_Momentum_withpim_%s",frame[f].Data()),p_mom[f]);
      FillHist(Form("p_CosT_withpim_%s",frame[f].Data()),p_cost[f]);
      FillHist(Form("p_Phi_withpim_%s",frame[f].Data()),p_phi[f]);
      FillHist(Form("p_MMass_withpim_%s",frame[f].Data()),p_mmass[f]);
      FillHist(Form("p_MMomentum_withpim_%s",frame[f].Data()),p_mmom[f]);
      FillHist(Form("p_MCosT_withpim_%s",frame[f].Data()),p_mcost[f]);
    }
    // neutron
    for(int f=1; f<2; f++){
      FillHist(Form("n_Momentum_withpim_%s",frame[f].Data()),n_mom[f]);
      FillHist(Form("n_CosT_withpim_%s",frame[f].Data()),n_cost[f]);
      FillHist(Form("n_Phi_withpim_%s",frame[f].Data()),n_phi[f]);
      FillHist(Form("n_MMass_withpim_%s",frame[f].Data()),n_mmass[f]);
      FillHist(Form("n_MMomentum_withpim_%s",frame[f].Data()),n_mmom[f]);
      FillHist(Form("n_MCosT_withpim_%s",frame[f].Data()),n_mcost[f]);
    }
    // proton + neutron
    for(int f=1; f<2; f++){
      FillHist(Form("pn_Momentum_withpim_%s",frame[f].Data()),pn_mom[f]);
      FillHist(Form("pn_CosT_withpim_%s",frame[f].Data()),pn_cost[f]);
      FillHist(Form("pn_Phi_withpim_%s",frame[f].Data()),pn_phi[f]);
      FillHist(Form("pn_MMass_withpim_%s",frame[f].Data()),pn_mmass[f]);
      FillHist(Form("pn_MMomentum_withpim_%s",frame[f].Data()),pn_mmom[f]);
      FillHist(Form("pn_MCosT_withpim_%s",frame[f].Data()),pn_mcost[f]);
    }
    /*-----------------------------------------*/
    /* 2D Plots                                */ 
    /*-----------------------------------------*/
    // MM(n) vs. MM(pn)
    for(int f=1; f<2; f++){
      FillHist(Form("n_MMvspn_MM_withpim_%s",frame[f].Data()),n_mmass[f],pn_mmass[f]);
    }
  }
  if(pipflag){
    // proton
    for(int f=1; f<2; f++){
      FillHist(Form("p_Momentum_withpip_%s",frame[f].Data()),p_mom[f]);
      FillHist(Form("p_CosT_withpip_%s",frame[f].Data()),p_cost[f]);
      FillHist(Form("p_Phi_withpip_%s",frame[f].Data()),p_phi[f]);
      FillHist(Form("p_MMass_withpip_%s",frame[f].Data()),p_mmass[f]);
      FillHist(Form("p_MMomentum_withpip_%s",frame[f].Data()),p_mmom[f]);
      FillHist(Form("p_MCosT_withpip_%s",frame[f].Data()),p_mcost[f]);
    }
    // neutron
    for(int f=1; f<2; f++){
      FillHist(Form("n_Momentum_withpip_%s",frame[f].Data()),n_mom[f]);
      FillHist(Form("n_CosT_withpip_%s",frame[f].Data()),n_cost[f]);
      FillHist(Form("n_Phi_withpip_%s",frame[f].Data()),n_phi[f]);
      FillHist(Form("n_MMass_withpip_%s",frame[f].Data()),n_mmass[f]);
      FillHist(Form("n_MMomentum_withpip_%s",frame[f].Data()),n_mmom[f]);
      FillHist(Form("n_MCosT_withpip_%s",frame[f].Data()),n_mcost[f]);
    }
    // proton + neutron
    for(int f=1; f<2; f++){
      FillHist(Form("pn_Momentum_withpip_%s",frame[f].Data()),pn_mom[f]);
      FillHist(Form("pn_CosT_withpip_%s",frame[f].Data()),pn_cost[f]);
      FillHist(Form("pn_Phi_withpip_%s",frame[f].Data()),pn_phi[f]);
      FillHist(Form("pn_MMass_withpip_%s",frame[f].Data()),pn_mmass[f]);
      FillHist(Form("pn_MMomentum_withpip_%s",frame[f].Data()),pn_mmom[f]);
      FillHist(Form("pn_MCosT_withpip_%s",frame[f].Data()),pn_mcost[f]);
    }
    /*-----------------------------------------*/
    /* 2D Plots                                */ 
    /*-----------------------------------------*/
    // MM(n) vs. MM(pn)
    for(int f=1; f<2; f++){
      FillHist(Form("n_MMvspn_MM_withpip_%s",frame[f].Data()),n_mmass[f],pn_mmass[f]);
    }
  }
  if(pflag&&pimflag){
    // proton
    for(int f=1; f<2; f++){
      FillHist(Form("p_Momentum_withppim_%s",frame[f].Data()),p_mom[f]);
      FillHist(Form("p_CosT_withppim_%s",frame[f].Data()),p_cost[f]);
      FillHist(Form("p_Phi_withppim_%s",frame[f].Data()),p_phi[f]);
      FillHist(Form("p_MMass_withppim_%s",frame[f].Data()),p_mmass[f]);
      FillHist(Form("p_MMomentum_withppim_%s",frame[f].Data()),p_mmom[f]);
      FillHist(Form("p_MCosT_withppim_%s",frame[f].Data()),p_mcost[f]);
    }
    // neutron
    for(int f=1; f<2; f++){
      FillHist(Form("n_Momentum_withppim_%s",frame[f].Data()),n_mom[f]);
      FillHist(Form("n_CosT_withppim_%s",frame[f].Data()),n_cost[f]);
      FillHist(Form("n_Phi_withppim_%s",frame[f].Data()),n_phi[f]);
      FillHist(Form("n_MMass_withppim_%s",frame[f].Data()),n_mmass[f]);
      FillHist(Form("n_MMomentum_withppim_%s",frame[f].Data()),n_mmom[f]);
      FillHist(Form("n_MCosT_withppim_%s",frame[f].Data()),n_mcost[f]);
    }
    // proton + neutron
    for(int f=1; f<2; f++){
      FillHist(Form("pn_Momentum_withppim_%s",frame[f].Data()),pn_mom[f]);
      FillHist(Form("pn_CosT_withppim_%s",frame[f].Data()),pn_cost[f]);
      FillHist(Form("pn_Phi_withppim_%s",frame[f].Data()),pn_phi[f]);
      FillHist(Form("pn_MMass_withppim_%s",frame[f].Data()),pn_mmass[f]);
      FillHist(Form("pn_MMomentum_withppim_%s",frame[f].Data()),pn_mmom[f]);
      FillHist(Form("pn_MCosT_withppim_%s",frame[f].Data()),pn_mcost[f]);
    }
    /*-----------------------------------------*/
    /* 2D Plots                                */ 
    /*-----------------------------------------*/
    // MM(n) vs. MM(pn)
    for(int f=1; f<2; f++){
      FillHist(Form("n_MMvspn_MM_withppim_%s",frame[f].Data()),n_mmass[f],pn_mmass[f]);
    }
  }


  return true;

}

bool MyAnalysisHeKpn::FillHist(TString name, double val1, int weight)
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

bool MyAnalysisHeKpn::FillHist(TString name, TString val1, int weight)
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

bool MyAnalysisHeKpn::FillHist(TString name, double val1, double val2, int weight)
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

bool MyAnalysisHeKpn::FillHist(TString name, TString val1, TString val2, int weight)
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

void MyAnalysisHeKpn::CutCondition()
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

bool MyAnalysisHeKpn::Initialize(ConfMan* confMan)
{
  std::cout << "### MyAnalysisHeKpn::Initialize ###" << std::endl;

  std::string ofname = confMan->GetOutFileName();
  ofname.insert(ofname.find(".root"),"_anaHeKpn");

  rtFile =  new TFile( Form("%s",ofname.c_str()), "RECREATE");
  rtFile -> cd();

  TString pname[8]  = {"pip","p","d","t","he","pim","k","o"};
  TString pname2[8] = {"#pi^{+}","p","d","t","he","#pi^{-}","K^{-}","o"};
  TString proname[6]  = {"pipi","k0","sp","sm","s","l"};
  TString proname2[6] = {"#pi^{+}#pi^{-}","K^{0}","#Sigma^{+}","#Sigma^{-}","#Sigma^{#pm}","#Lambda"};

  new TH1F( Form("CDH_Multiplicity"), Form("Multiplicity CDH;Multiplicity;Counts"), 36+1, -0.5, 36+0.5 );
  new TH1F( Form("BVC_Multiplicity"), Form("Multiplicity BVC;Multiplicity;Counts"), 8+1, -0.5, 8+0.5 );
  new TH1F( Form("CVC_Multiplicity"), Form("Multiplicity CVC;Multiplicity;Counts"), 34+1, -0.5, 34+0.5 );
  new TH1F( Form("PC_Multiplicity"), Form("Multiplicity PC;Multiplicity;Counts"), 27+1, -0.5, 27+0.5 );
  new TH1F( Form("IH_Multiplicity"), Form("Multiplicity IH;Multiplicity;Counts"), 24+1, -0.5, 24+0.5 );

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

  // proton
  new TH1F("p_Momentum_CM","Momentum of p;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("p_CosT_CM","cos#theta of p;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("p_Phi_CM","#phi of p;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("p_MMass_CM", "^{3}He(K^{-},p)X missing mass;MM(p) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("p_MMomentum_CM", "^{3}He(K^{-},p)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("p_MCosT_CM", "^{3}He(K^{-},p)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // neutron
  new TH1F("n_Momentum_CM","Momentum of n;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("n_CosT_CM","cos#theta of n;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("n_Phi_CM","#phi of n;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("n_MMass_CM", "^{3}He(K^{-},n)X missing mass;MM(p) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("n_MMomentum_CM", "^{3}He(K^{-},n)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("n_MCosT_CM", "^{3}He(K^{-},n)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // proton + neutron
  new TH1F("pn_Momentum_CM","Momentum of pn;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pn_CosT_CM","cos#theta of pn;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("pn_Phi_CM","#phi of pn;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("pn_MMass_CM", "^{3}He(K^{-},pn)X missing mass;MM(p) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pn_MMomentum_CM", "^{3}He(K^{-},pn)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pn_MCosT_CM", "^{3}He(K^{-},pn)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  /* 2D Plots */
  new TH2F("n_MMvspn_MM_CM", "^{3}He(K^{-},n)X missing mass vs. ^{3}He(K^{-},pn)X;MM(n) (GeV/c^{2});MM(pn) (GeV/c^{2})", 200, 2.0, 4.0, 200, 1.0, 3.0);

  // with proton
  // proton
  new TH1F("p_Momentum_withp_CM","Momentum of p;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("p_CosT_withp_CM","cos#theta of p;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("p_Phi_withp_CM","#phi of p;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("p_MMass_withp_CM", "^{3}He(K^{-},p)X missing mass;MM(p) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("p_MMomentum_withp_CM", "^{3}He(K^{-},p)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("p_MCosT_withp_CM", "^{3}He(K^{-},p)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // neutron
  new TH1F("n_Momentum_withp_CM","Momentum of n;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("n_CosT_withp_CM","cos#theta of n;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("n_Phi_withp_CM","#phi of n;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("n_MMass_withp_CM", "^{3}He(K^{-},n)X missing mass;MM(p) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("n_MMomentum_withp_CM", "^{3}He(K^{-},n)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("n_MCosT_withp_CM", "^{3}He(K^{-},n)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // proton + neutron
  new TH1F("pn_Momentum_withp_CM","Momentum of pn;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pn_CosT_withp_CM","cos#theta of pn;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("pn_Phi_withp_CM","#phi of pn;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("pn_MMass_withp_CM", "^{3}He(K^{-},pn)X missing mass;MM(p) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pn_MMomentum_withp_CM", "^{3}He(K^{-},pn)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pn_MCosT_withp_CM", "^{3}He(K^{-},pn)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  /* 2D Plots */
  new TH2F("n_MMvspn_MM_withp_CM", "^{3}He(K^{-},n)X missing mass vs. ^{3}He(K^{-},pn)X;MM(n) (GeV/c^{2});MM(pn) (GeV/c^{2})", 200, 2.0, 4.0, 200, 1.0, 3.0);

  // with pi-
  // proton
  new TH1F("p_Momentum_withpim_CM","Momentum of p;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("p_CosT_withpim_CM","cos#theta of p;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("p_Phi_withpim_CM","#phi of p;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("p_MMass_withpim_CM", "^{3}He(K^{-},p)X missing mass;MM(p) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("p_MMomentum_withpim_CM", "^{3}He(K^{-},p)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("p_MCosT_withpim_CM", "^{3}He(K^{-},p)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // neutron
  new TH1F("n_Momentum_withpim_CM","Momentum of n;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("n_CosT_withpim_CM","cos#theta of n;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("n_Phi_withpim_CM","#phi of n;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("n_MMass_withpim_CM", "^{3}He(K^{-},n)X missing mass;MM(p) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("n_MMomentum_withpim_CM", "^{3}He(K^{-},n)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("n_MCosT_withpim_CM", "^{3}He(K^{-},n)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // proton + neutron
  new TH1F("pn_Momentum_withpim_CM","Momentum of pn;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pn_CosT_withpim_CM","cos#theta of pn;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("pn_Phi_withpim_CM","#phi of pn;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("pn_MMass_withpim_CM", "^{3}He(K^{-},pn)X missing mass;MM(p) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pn_MMomentum_withpim_CM", "^{3}He(K^{-},pn)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pn_MCosT_withpim_CM", "^{3}He(K^{-},pn)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  /* 2D Plots */
  new TH2F("n_MMvspn_MM_withpim_CM", "^{3}He(K^{-},n)X missing mass vs. ^{3}He(K^{-},pn)X;MM(n) (GeV/c^{2});MM(pn) (GeV/c^{2})", 200, 2.0, 4.0, 200, 1.0, 3.0);

  // with pi+
  // proton
  new TH1F("p_Momentum_withpip_CM","Momentum of p;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("p_CosT_withpip_CM","cos#theta of p;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("p_Phi_withpip_CM","#phi of p;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("p_MMass_withpip_CM", "^{3}He(K^{-},p)X missing mass;MM(p) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("p_MMomentum_withpip_CM", "^{3}He(K^{-},p)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("p_MCosT_withpip_CM", "^{3}He(K^{-},p)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // neutron
  new TH1F("n_Momentum_withpip_CM","Momentum of n;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("n_CosT_withpip_CM","cos#theta of n;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("n_Phi_withpip_CM","#phi of n;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("n_MMass_withpip_CM", "^{3}He(K^{-},n)X missing mass;MM(p) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("n_MMomentum_withpip_CM", "^{3}He(K^{-},n)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("n_MCosT_withpip_CM", "^{3}He(K^{-},n)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // proton + neutron
  new TH1F("pn_Momentum_withpip_CM","Momentum of pn;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pn_CosT_withpip_CM","cos#theta of pn;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("pn_Phi_withpip_CM","#phi of pn;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("pn_MMass_withpip_CM", "^{3}He(K^{-},pn)X missing mass;MM(p) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pn_MMomentum_withpip_CM", "^{3}He(K^{-},pn)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pn_MCosT_withpip_CM", "^{3}He(K^{-},pn)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  /* 2D Plots */
  new TH2F("n_MMvspn_MM_withpip_CM", "^{3}He(K^{-},n)X missing mass vs. ^{3}He(K^{-},pn)X;MM(n) (GeV/c^{2});MM(pn) (GeV/c^{2})", 200, 2.0, 4.0, 200, 1.0, 3.0);

  // with proton & pi-
  // proton
  new TH1F("p_Momentum_withppim_CM","Momentum of p;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("p_CosT_withppim_CM","cos#theta of p;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("p_Phi_withppim_CM","#phi of p;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("p_MMass_withppim_CM", "^{3}He(K^{-},p)X missing mass;MM(p) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("p_MMomentum_withppim_CM", "^{3}He(K^{-},p)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("p_MCosT_withppim_CM", "^{3}He(K^{-},p)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // neutron
  new TH1F("n_Momentum_withppim_CM","Momentum of n;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("n_CosT_withppim_CM","cos#theta of n;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("n_Phi_withppim_CM","#phi of n;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("n_MMass_withppim_CM", "^{3}He(K^{-},n)X missing mass;MM(p) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("n_MMomentum_withppim_CM", "^{3}He(K^{-},n)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("n_MCosT_withppim_CM", "^{3}He(K^{-},n)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // proton + neutron
  new TH1F("pn_Momentum_withppim_CM","Momentum of pn;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pn_CosT_withppim_CM","cos#theta of pn;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("pn_Phi_withppim_CM","#phi of pn;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("pn_MMass_withppim_CM", "^{3}He(K^{-},pn)X missing mass;MM(p) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pn_MMomentum_withppim_CM", "^{3}He(K^{-},pn)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pn_MCosT_withppim_CM", "^{3}He(K^{-},pn)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  /* 2D Plots */
  new TH2F("n_MMvspn_MM_withppim_CM", "^{3}He(K^{-},n)X missing mass vs. ^{3}He(K^{-},pn)X;MM(n) (GeV/c^{2});MM(pn) (GeV/c^{2})", 200, 2.0, 4.0, 200, 1.0, 3.0);

  std::cout << "== Finish Histogram Initialization ==" << std::endl;
  return true;

}
