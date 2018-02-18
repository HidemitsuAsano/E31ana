// MyAnalysisHeKpipp.cpp

#include "MyAnalysisHeKpipp.h"

MyAnalysisHeKpipp::MyAnalysisHeKpipp(TFile* rt, ConfMan* conf)
{
  Initialize(conf);
  CutCondition();
  Clear();
}

MyAnalysisHeKpipp::~MyAnalysisHeKpipp()
{
  Clear();
  rtFile->cd();
  rtFile->Write();
  rtFile->Close();
}

void MyAnalysisHeKpipp::Clear()
{
}

bool MyAnalysisHeKpipp::DoAnalysis(ConfMan* conf, EventHeader* header, BeamLineHitMan* blMan, BeamLineTrackMan* bltrackMan, CDSHitMan* cdsMan, CDSTrackingMan* cdstrackMan, Particle* particle)
{
  rtFile->cd();

  FillHist("EventNumber",0);
  DetectorList *dlist=DetectorList::GetInstance();
  if(particle->nBeam()!=1) return false;
  pBeam* beam = particle->beam(0);
  FillHist("EventNumber",1);

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
  if(MulBVC!=0||MulCVC!=0||MulPC!=0||MulCDH!=3) return true;

  if(particle->nCDS()!=3) return true;
  FillHist("EventNumber",2);
  if(particle->nProton()!=2) return true;
  if(particle->nPiminus()!=1) return true;
  FillHist("EventNumber",3);

  // ############### //
  // Vertex decision //
  // ############### //
  TVector3 vertex;
  double vdis = 9999;
  bool FIDUCIAL = false;
  int fcds = -1;
  for(int it=0; it<particle->nProton(); it++){
    pCDS* cds = particle->proton(it);
    if(vdis>9998||vdis>cds->vdis()){
      vdis = cds->vdis();
      vertex = cds->vbeam();
      if(GeomTools::GetID(vertex)==CID_Fiducial){
        FIDUCIAL = true;
        fcds=it;
      }
    }
  }
  if(!FIDUCIAL){ return true; }
  if(fcds<0){ return true; }
  FillHist("EventNumber",4);

  pCDS* pim = particle->pim(0);
  pCDS* p1 = particle->proton(fcds);
  pCDS* p2 = particle->proton(1-fcds);

  /* CDC Track Chi-Square Cut */
  if(pim->chi()>30) { return true; }
  if(p1->chi()>30) { return true; }
  if(p2->chi()>30) { return true; }
  /* CDC Track Chi-Square Cut */
  FillHist("EventNumber",5);

  // ######################### //
  // 2 Helix calculation       //
  // ######################### //
  //pCDS* pip1 = 0;
  //pip1 = TrackTools::Calc2HelixAll(beam,pim,p1,cdstrackMan,true);
  //pCDS* pip2 = 0;
  //pip2 = TrackTools::Calc2HelixAll(beam,pim,p2,cdstrackMan,true);
  //if(pip1==0||pip2==0){ return true; }
  //FillHist("EventNumber",6);

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
  //TLorentzVector Lpip1 = pip1->GetLorentzVector();
  //TLorentzVector cLpip1 = pip1->GetLorentzVector(); cLpip1.Boost(-boost);
  //TLorentzVector Lpip2 = pip2->GetLorentzVector();
  //TLorentzVector cLpip2 = pip2->GetLorentzVector(); cLpip2.Boost(-boost);
  TLorentzVector Lpip1 = Lpim + Lp1;
  TLorentzVector cLpip1 = cLpim + cLp1;
  TLorentzVector Lpip2 = Lpim + Lp2;
  TLorentzVector cLpip2 = cLpim + cLp2;
  TLorentzVector Lpipp = Lpim + Lp1 + Lp2;
  TLorentzVector cLpipp = cLpim + cLp1 + cLp2;

  double pim_mass[2]  = { Lpim.M()                           , cLpim.M()};
  double pim_mom[2]   = { Lpim.Vect().Mag()                  , cLpim.Vect().Mag()};
  double pim_cost[2]  = { Lpim.Vect().CosTheta()             , cLpim.Vect().CosTheta()};
  double pim_phi[2]   = { Lpim.Vect().Phi()                  , cLpim.Vect().Phi()};
  double pim_mmass[2] = { (Ltgt+Lbeam-Lpim).M()              , (cLtgt+cLbeam-cLpim).M()};
  double pim_mmom[2]  = { (Ltgt+Lbeam-Lpim).Vect().Mag()     , (cLtgt+cLbeam-cLpim).Vect().Mag()};
  double pim_mcost[2] = { (Ltgt+Lbeam-Lpim).Vect().CosTheta(), (cLtgt+cLbeam-cLpim).Vect().CosTheta()};

  double p1_mass[2]  = { Lp1.M()                           , cLp1.M()};
  double p1_mom[2]   = { Lp1.Vect().Mag()                  , cLp1.Vect().Mag()};
  double p1_cost[2]  = { Lp1.Vect().CosTheta()             , cLp1.Vect().CosTheta()};
  double p1_phi[2]   = { Lp1.Vect().Phi()                  , cLp1.Vect().Phi()};
  double p1_mmass[2] = { (Ltgt+Lbeam-Lp1).M()              , (cLtgt+cLbeam-cLp1).M()};
  double p1_mmom[2]  = { (Ltgt+Lbeam-Lp1).Vect().Mag()     , (cLtgt+cLbeam-cLp1).Vect().Mag()};
  double p1_mcost[2] = { (Ltgt+Lbeam-Lp1).Vect().CosTheta(), (cLtgt+cLbeam-cLp1).Vect().CosTheta()};

  double p2_mass[2]  = { Lp2.M()                           , cLp2.M()};
  double p2_mom[2]   = { Lp2.Vect().Mag()                  , cLp2.Vect().Mag()};
  double p2_cost[2]  = { Lp2.Vect().CosTheta()             , cLp2.Vect().CosTheta()};
  double p2_phi[2]   = { Lp2.Vect().Phi()                  , cLp2.Vect().Phi()};
  double p2_mmass[2] = { (Ltgt+Lbeam-Lp2).M()              , (cLtgt+cLbeam-cLp2).M()};
  double p2_mmom[2]  = { (Ltgt+Lbeam-Lp2).Vect().Mag()     , (cLtgt+cLbeam-cLp2).Vect().Mag()};
  double p2_mcost[2] = { (Ltgt+Lbeam-Lp2).Vect().CosTheta(), (cLtgt+cLbeam-cLp2).Vect().CosTheta()};

  double pip1_mass[2]  = { Lpip1.M()                           , cLpip1.M()};
  double pip1_mom[2]   = { Lpip1.Vect().Mag()                  , cLpip1.Vect().Mag()};
  double pip1_cost[2]  = { Lpip1.Vect().CosTheta()             , cLpip1.Vect().CosTheta()};
  double pip1_phi[2]   = { Lpip1.Vect().Phi()                  , cLpip1.Vect().Phi()};
  double pip1_mmass[2] = { (Ltgt+Lbeam-Lpip1).M()              , (cLtgt+cLbeam-cLpip1).M()};
  double pip1_mmom[2]  = { (Ltgt+Lbeam-Lpip1).Vect().Mag()     , (cLtgt+cLbeam-cLpip1).Vect().Mag()};
  double pip1_mcost[2] = { (Ltgt+Lbeam-Lpip1).Vect().CosTheta(), (cLtgt+cLbeam-cLpip1).Vect().CosTheta()};

  double pip2_mass[2]  = { Lpip2.M()                           , cLpip2.M()};
  double pip2_mom[2]   = { Lpip2.Vect().Mag()                  , cLpip2.Vect().Mag()};
  double pip2_cost[2]  = { Lpip2.Vect().CosTheta()             , cLpip2.Vect().CosTheta()};
  double pip2_phi[2]   = { Lpip2.Vect().Phi()                  , cLpip2.Vect().Phi()};
  double pip2_mmass[2] = { (Ltgt+Lbeam-Lpip2).M()              , (cLtgt+cLbeam-cLpip2).M()};
  double pip2_mmom[2]  = { (Ltgt+Lbeam-Lpip2).Vect().Mag()     , (cLtgt+cLbeam-cLpip2).Vect().Mag()};
  double pip2_mcost[2] = { (Ltgt+Lbeam-Lpip2).Vect().CosTheta(), (cLtgt+cLbeam-cLpip2).Vect().CosTheta()};

  double pipp_mass[2]  = { Lpipp.M()                           , cLpipp.M()};
  double pipp_mom[2]   = { Lpipp.Vect().Mag()                  , cLpipp.Vect().Mag()};
  double pipp_cost[2]  = { Lpipp.Vect().CosTheta()             , cLpipp.Vect().CosTheta()};
  double pipp_phi[2]   = { Lpipp.Vect().Phi()                  , cLpipp.Vect().Phi()};
  double pipp_mmass[2] = { (Ltgt+Lbeam-Lpipp).M()              , (cLtgt+cLbeam-cLpipp).M()};
  double pipp_mmom[2]  = { (Ltgt+Lbeam-Lpipp).Vect().Mag()     , (cLtgt+cLbeam-cLpipp).Vect().Mag()};
  double pipp_mcost[2] = { (Ltgt+Lbeam-Lpipp).Vect().CosTheta(), (cLtgt+cLbeam-cLpipp).Vect().CosTheta()};

  /* flag */
  bool lamflag1 = false; /* lam flag */
  bool sbllamflag1 = false; /* lower sideband of lam flag */
  bool sbulamflag1 = false; /* upper sideband of lam flag */
  {
    double pro_mass = pip1_mass[1];
    if(lamll<pro_mass&&pro_mass<lamul) lamflag1 = true;
    if(sbllamll<pro_mass&&pro_mass<sbllamul) sbllamflag1 = true;
    if(sbulamll<pro_mass&&pro_mass<sbulamul) sbulamflag1 = true;
  }
  bool lamflag2 = false; /* lam flag */
  bool sbllamflag2 = false; /* lower sideband of lam flag */
  bool sbulamflag2 = false; /* upper sideband of lam flag */
  {
    double pro_mass = pip2_mass[1];
    if(lamll<pro_mass&&pro_mass<lamul) lamflag2 = true;
    if(sbllamll<pro_mass&&pro_mass<sbllamul) sbllamflag2 = true;
    if(sbulamll<pro_mass&&pro_mass<sbulamul) sbulamflag2 = true;
  }
  bool mnflag = false; /* mm flag */
  bool sblmnflag = false; /* lower sideband of mm flag */
  bool sbumnflag = false; /* upper sideband of mm flag */
  {
    double pro_mmass = pipp_mmass[1];
    if(mnll<pro_mmass&&pro_mmass<mnul) mnflag = true;
    if(sblmnll<pro_mmass&&pro_mmass<sblmnul) sblmnflag = true;
    if(sbumnll<pro_mmass&&pro_mmass<sbumnul) sbumnflag = true;
  }

  //if(lamflag1){ return true; }
  //if(lamflag2){ return true; }
  //if(k0flag){ return true; }
  //if(!k0flag){ return true; }
  //if(k0flag){ return true; }
  //if(!smflag){ return true; }
  //if(smflag){ return true; }
  //if(!spflag){ return true; }
  //if(spflag){ return true; }
  //if(!lamflag1&&!lamflag2){ return true; }
  //if(!spflag&&!smflag){ return true; }
  //if(spflag||smflag){ return true; }
  //if(!mnflag){ return true; }
  //if(mnflag){ return true; }
  // ############### //
  // NC hit decision //
  // ############### //
  //double time = 9999;
  //int fnc = -1;
  //for(int it=0; it<particle->nNC(); it++){
  //  pNC* nc = particle->nc(it);
  //  nc->CalcMom(beam,vertex);
  //  FillHist("FWD_OverbetavsMomentum",1.0/nc->beta(),nc->mom().Mag());
  //  FillHist("FWD_OverbetavsEnergy",1.0/nc->beta(),nc->energy());
  //  FillHist("FWD_TOFvsMomentum",nc->tof(),nc->mom().Mag());
  //  FillHist("FWD_HitPosition",nc->hitpos().X(),nc->hitpos().Y());
  //  double seg = (nc->seg()-1)%16+1;
  //  double lay = (nc->seg()-1)/16+1;
  //  FillHist("FWD_HitSegment",seg,lay);
  //  if(nc->pid()==F_Neutron) {
  //    if(nc->energy()>8.0 && (time>9998||time>nc->time())){
  //      time = nc->time();
  //      fnc = it;
  //    }
  //  }
  //}
  //pNC*  nc =0;
  //if(fnc!=-1)  nc =particle->nc(fnc);
  //if(nc==0) return true;

  /*-----------------------------------------*/
  /* All                                     */ 
  /*-----------------------------------------*/
  // ######################### //
  // DCA                       //
  // ######################### //
  /* pim DCA */
  FillHist("pim_VDIS",pim->vdis());
  /* p1 DCA */
  FillHist("p1_VDIS",p1->vdis());
  /* p2 DCA */
  FillHist("p2_VDIS",p1->vdis());
  /* p1 vs p2 DCA */
  FillHist("p1_VDISvsp2_VDIS",p1->vdis(),p2->vdis());
  /* Vertex */
  FillHist("p1_Vertex_XY",p1->vertex().X(),p1->vertex().Y());
  FillHist("p1_Vertex_ZX",p1->vertex().Z(),p1->vertex().X());
  FillHist("p1_Vertex_ZY",p1->vertex().Z(),p1->vertex().Y());
  FillHist("p2_Vertex_XY",p2->vertex().X(),p2->vertex().Y());
  FillHist("p2_Vertex_ZX",p2->vertex().Z(),p2->vertex().X());
  FillHist("p2_Vertex_ZY",p2->vertex().Z(),p2->vertex().Y());
  ///* Vertex */
  //FillHist("pip1_Vertex_XY",pip1->vertex().X(),pip1->vertex().Y());
  //FillHist("pip1_Vertex_ZX",pip1->vertex().Z(),pip1->vertex().X());
  //FillHist("pip1_Vertex_ZY",pip1->vertex().Z(),pip1->vertex().Y());
  //FillHist("pip2_Vertex_XY",pip2->vertex().X(),pip2->vertex().Y());
  //FillHist("pip2_Vertex_ZX",pip2->vertex().Z(),pip2->vertex().X());
  //FillHist("pip2_Vertex_ZY",pip2->vertex().Z(),pip2->vertex().Y());
  ///* VDIS */
  //FillHist("pip1_VDIS",pip1->vdis());
  //FillHist("pip2_VDIS",pip2->vdis());
  //FillHist("pip1_VDISvspip2_VDIS",pip1->vdis(),pip2->vdis());
  ///* VBDIS */
  //FillHist("pip1_VBDIS",pip1->vbdis());
  //FillHist("pip2_VBDIS",pip2->vbdis());
  //FillHist("pip1_VBDISvspip2_VBDIS",pip1->vbdis(),pip2->vbdis());
  ///* PBDCA */
  //FillHist("pip1_PBDCA",pip1->pbdca());
  //FillHist("pip2_PBDCA",pip2->pbdca());
  //FillHist("pip1_PBDCAvspip2_PBDCA",pip1->pbdca(),pip2->pbdca());

  // ######################### //
  // Invariant/Missing Mass    //
  // ######################### //
  // pi-
  for(int f=1; f<2; f++){
    FillHist(Form("pim_CosTvsMomentum_%s",frame[f].Data()),pim_cost[f],pim_mom[f]);
    FillHist(Form("pim_CosTvsPhi_%s",frame[f].Data()),pim_cost[f],pim_phi[f]);
    FillHist(Form("HeKpim_MMvsMomentum_%s",frame[f].Data()),pim_mmass[f],pim_mmom[f]);
    FillHist(Form("HeKpim_MMvsCosT_%s",frame[f].Data()),pim_mmass[f],pim_mcost[f]);
    FillHist(Form("pim_IMvsMomentum_%s",frame[f].Data()),pim_mass[f],pim_mom[f]);
    FillHist(Form("pim_IMvsCosT_%s",frame[f].Data()),pim_mass[f],pim_cost[f]);
  }
  // proton1
  for(int f=1; f<2; f++){
    FillHist(Form("p1_CosTvsMomentum_%s",frame[f].Data()),p1_cost[f],p1_mom[f]);
    FillHist(Form("p1_CosTvsPhi_%s",frame[f].Data()),p1_cost[f],p1_phi[f]);
    FillHist(Form("H3Kp1_MMvsMomentum_%s",frame[f].Data()),p1_mmass[f],p1_mmom[f]);
    FillHist(Form("H3Kp1_MMvsCosT_%s",frame[f].Data()),p1_mmass[f],p1_mcost[f]);
    FillHist(Form("p1_IMvsMomentum_%s",frame[f].Data()),p1_mass[f],p1_mom[f]);
    FillHist(Form("p1_IMvsCosT_%s",frame[f].Data()),p1_mass[f],p1_cost[f]);
  }
  // proton2
  for(int f=1; f<2; f++){
    FillHist(Form("p2_CosTvsMomentum_%s",frame[f].Data()),p2_cost[f],p2_mom[f]);
    FillHist(Form("p2_CosTvsPhi_%s",frame[f].Data()),p2_cost[f],p2_phi[f]);
    FillHist(Form("H3Kp2_MMvsMomentum_%s",frame[f].Data()),p2_mmass[f],p2_mmom[f]);
    FillHist(Form("H3Kp2_MMvsCosT_%s",frame[f].Data()),p2_mmass[f],p2_mcost[f]);
    FillHist(Form("p2_IMvsMomentum_%s",frame[f].Data()),p2_mass[f],p2_mom[f]);
    FillHist(Form("p2_IMvsCosT_%s",frame[f].Data()),p2_mass[f],p2_cost[f]);
  }
  // pi- proton1
  for(int f=1; f<2; f++){
    FillHist(Form("pip1_CosTvsMomentum_%s",frame[f].Data()),pip1_cost[f],pip1_mom[f]);
    FillHist(Form("pip1_CosTvsPhi_%s",frame[f].Data()),pip1_cost[f],pip1_phi[f]);
    FillHist(Form("HeKpip1_MMvsMomentum_%s",frame[f].Data()),pip1_mmass[f],pip1_mmom[f]);
    FillHist(Form("HeKpip1_MMvsCosT_%s",frame[f].Data()),pip1_mmass[f],pip1_mcost[f]);
    FillHist(Form("pip1_IMvsMomentum_%s",frame[f].Data()),pip1_mass[f],pip1_mom[f]);
    FillHist(Form("pip1_IMvsCosT_%s",frame[f].Data()),pip1_mass[f],pip1_cost[f]);
  }
  // pi- proton2
  for(int f=1; f<2; f++){
    FillHist(Form("pip2_CosTvsMomentum_%s",frame[f].Data()),pip2_cost[f],pip2_mom[f]);
    FillHist(Form("pip2_CosTvsPhi_%s",frame[f].Data()),pip2_cost[f],pip2_phi[f]);
    FillHist(Form("HeKpip2_MMvsMomentum_%s",frame[f].Data()),pip2_mmass[f],pip2_mmom[f]);
    FillHist(Form("HeKpip2_MMvsCosT_%s",frame[f].Data()),pip2_mmass[f],pip2_mcost[f]);
    FillHist(Form("pip2_IMvsMomentum_%s",frame[f].Data()),pip2_mass[f],pip2_mom[f]);
    FillHist(Form("pip2_IMvsCosT_%s",frame[f].Data()),pip2_mass[f],pip2_cost[f]);
  }
  // pi- proton1 proton2
  for(int f=1; f<2; f++){
    FillHist(Form("pipp_CosTvsMomentum_%s",frame[f].Data()),pipp_cost[f],pipp_mom[f]);
    FillHist(Form("pipp_CosTvsPhi_%s",frame[f].Data()),pipp_cost[f],pipp_phi[f]);
    FillHist(Form("HeKpipp_MMvsMomentum_%s",frame[f].Data()),pipp_mmass[f],pipp_mmom[f]);
    FillHist(Form("HeKpipp_MMvsCosT_%s",frame[f].Data()),pipp_mmass[f],pipp_mcost[f]);
    FillHist(Form("pipp_IMvsMomentum_%s",frame[f].Data()),pipp_mass[f],pipp_mom[f]);
    FillHist(Form("pipp_IMvsCosT_%s",frame[f].Data()),pipp_mass[f],pipp_cost[f]);
  }

  /*-----------------------------------------*/
  /* Advanced                                */ 
  /*-----------------------------------------*/
  // IM(pip1) vs. IM(pip2)
  for(int f=1; f<2; f++){
    FillHist(Form("pip1_IMvspip2_IM_%s",frame[f].Data()),pip1_mass[f],pip2_mass[f]);
  }
  // IM(pipp) vs. MM(pipp)
  for(int f=1; f<2; f++){
    FillHist(Form("pipp_IMvsHeKpipp_MM_%s",frame[f].Data()),pipp_mass[f],pipp_mmass[f]);
  }

  return true;

}

bool MyAnalysisHeKpipp::FillHist(TString name, double val1, int weight)
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

bool MyAnalysisHeKpipp::FillHist(TString name, TString val1, int weight)
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

bool MyAnalysisHeKpipp::FillHist(TString name, double val1, double val2, int weight)
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

bool MyAnalysisHeKpipp::FillHist(TString name, TString val1, TString val2, int weight)
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

void MyAnalysisHeKpipp::CutCondition()
{
  /* K0 mass cut */ 
  double mean  = 0.4976;
  double sigma = 0.0069;
  k0ll = mean-2*sigma; k0ul = mean+2*sigma;
  sblk0ll = mean-5*sigma; sblk0ul = mean-3*sigma;
  sbuk0ll = mean+3*sigma; sbuk0ul = mean+5*sigma;
  /* Missing neutron mass cut */ 
  mean  = 0.9379;
  sigma = 0.0230;
  mnll = mean-2*sigma; mnul = mean+2*sigma;
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
  mean  = 1.1160;
  sigma = 0.0029;
  lamll = mean-2*sigma; lamul = mean+2*sigma;
  sbllamll = mean-5*sigma; sbllamul = mean-3*sigma;
  sbulamll = mean+3*sigma; sbulamul = mean+5*sigma;
}

bool MyAnalysisHeKpipp::Initialize(ConfMan* confMan)
{
  std::cout << "### MyAnalysisHeKpipp::Initialize ###" << std::endl;

  std::string ofname = confMan->GetOutFileName();
  ofname.insert(ofname.find(".root"),"_anaHeKpipp");

  rtFile =  new TFile( Form("%s",ofname.c_str()), "RECREATE");
  rtFile -> cd();

  TString pname[8]  = {"pip","p","d","t","he","pim","k","o"};
  TString pname2[8] = {"#pi^{+}","p","d","t","he","#pi^{-}","K^{-}","o"};
  TString proname[6]  = {"pipi","k0","sp","sm","s","l"};
  TString proname2[6] = {"#pi^{+}#pi^{-}","K^{0}","#Sigma^{+}","#Sigma^{-}","#Sigma^{#pm}","#Lambda"};

  new TH1F( "EventNumber", "Number of Events", 10, 0, 10);
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
  new TH1F( "pim_VDIS", "VDIS of #pi^{-};VDIS (cm);Coutns", 100, 0, 10);
  new TH1F( "p1_VDIS", "VDIS of p_{1};VDIS (cm);Coutns", 100, 0, 10);
  new TH1F( "p2_VDIS", "VDIS of p_{2};VDIS (cm);Coutns", 100, 0, 10);
  new TH2F( "p1_VDISvsp2_VDIS", "VDIS of p_{1} vs p_{2};VDIS (cm);VDIS (cm)", 100, 0, 10, 100, 0, 10);
  new TH1F( "pip1_VDIS", "VDIS of #pi^{-}p_{1};VDIS (cm);Coutns", 100, 0, 10);
  new TH1F( "pip2_VDIS", "VDIS of #pi^{-}p_{2};VDIS (cm);Coutns", 100, 0, 10);
  new TH2F( "pip1_VDISvspip2_VDIS", "VDIS of #pi^{-}p_{1} vs #pi^{-}p_{2};VDIS (cm);VDIS (cm)", 100, 0, 10, 10, 0, 10);
  new TH1F( "pip1_VBDIS", "VBDIS of #pi^{-}p_{1};VBDIS (cm);Coutns", 100, 0, 10);
  new TH1F( "pip2_VBDIS", "VBDIS of #pi^{-}p_{2};VBDIS (cm);Coutns", 100, 0, 10);
  new TH2F( "pip1_VBDISvspip2_VBDIS", "VBDIS of #pi^{-}p_{1} vs #pi^{-}p_{2};VBDIS (cm);VBDIS (cm)", 100, 0, 10, 100, 0, 10);
  new TH1F( "pip1_PBDCA", "PBDCA of #pi^{-}p_{1};PBDCA (cm);Coutns", 100, 0, 10);
  new TH1F( "pip2_PBDCA", "PBDCA of #pi^{-}p_{2};PBDCA (cm);Coutns", 100, 0, 10);
  new TH2F( "pip1_PBDCAvspip2_PBDCA", "PBDCA of #pi^{-}p_{1} vs #pi^{-}p_{2};PBDCA (cm);PBDCA (cm)", 100, 0, 10, 100, 0, 10);
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
  // pi-
  new TH2F("pim_CosTvsMomentum_CM","cos#theta vs. momentum of #pi^{-};cos#theta;Momentum (GeV/c)", 2000, -1.0, 1.0, 2000, 0.0, 2.0);
  new TH2F("pim_CosTvsPhi_CM","cos#theta vs. #phi of #pi^{-};cos#theta;#phi", 2000, -1.0, 1.0, 1600, 0.0, 3.2);
  new TH2F("HeKpim_MMvsMomentum_CM", "^{3}He(K^{-},#pi^{-})X missing mass vs. missing momentum;MM(#pi^{-}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 2.0, 4.0, 2000, 0.0, 2.0);
  new TH2F("HeKpim_MMvsCosT_CM", "^{3}He(K^{-},#pi^{-})X missing mass vs. cos#theta;MM(#pi^{-}) (GeV/c^{2});cos#theta", 2000, 2.0, 4.0, 2000, -1.0, 1.0);
  new TH2F("pim_IMvsMomentum_CM", "(#pi^{-}) invariant mass vs. momentum;IM(#pi^{-}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("pim_IMvsCosT_CM", "(#pi^{-}) invariant mass vs. cos#theta;IM(#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  // proton1
  new TH2F("p1_CosTvsMomentum_CM","cos#theta vs. momentum of p_{1};cos#theta;Momentum (GeV/c)", 2000, -1.0, 1.0, 2000, 0.0, 2.0);
  new TH2F("p1_CosTvsPhi_CM","cos#theta vs. #phi of p_{1};cos#theta;#phi", 2000, -1.0, 1.0, 1600, 0.0, 3.2);
  new TH2F("HeKp1_MMvsMomentum_CM", "^{3}He(K^{-},p_{1})X missing mass vs. missing momentum;MM(p_{1}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 2.0, 4.0, 2000, 0.0, 2.0);
  new TH2F("HeKp1_MMvsCosT_CM", "^{3}He(K^{-},p_{1})X missing mass vs. cos#theta;MM(p_{1}) (GeV/c^{2});cos#theta", 2000, 2.0, 4.0, 2000, -1.0, 1.0);
  new TH2F("p1_IMvsMomentum_CM", "(p_{1}) invariant mass vs. momentum;IM(p_{1}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("p1_IMvsCosT_CM", "(p_{1}) invariant mass vs. cos#theta;IM(p_{1}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  // proton2
  new TH2F("p2_CosTvsMomentum_CM","cos#theta vs. momentum of p_{2};cos#theta;Momentum (GeV/c)", 2000, -1.0, 1.0, 2000, 0.0, 2.0);
  new TH2F("p2_CosTvsPhi_CM","cos#theta vs. #phi of p_{2};cos#theta;#phi", 2000, -1.0, 1.0, 1600, 0.0, 3.2);
  new TH2F("HeKp2_MMvsMomentum_CM", "^{3}He(K^{-},p_{2})X missing mass vs. missing momentum;MM(p_{2}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 2.0, 4.0, 2000, 0.0, 2.0);
  new TH2F("HeKp2_MMvsCosT_CM", "^{3}He(K^{-},p_{2})X missing mass vs. cos#theta;MM(p_{2}) (GeV/c^{2});cos#theta", 2000, 2.0, 4.0, 2000, -1.0, 1.0);
  new TH2F("p2_IMvsMomentum_CM", "(p_{2}) invariant mass vs. momentum;IM(p_{2}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("p2_IMvsCosT_CM", "(p_{2}) invariant mass vs. cos#theta;IM(p_{2}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  // pi- proton1
  new TH2F("pip1_CosTvsMomentum_CM","cos#theta vs. momentum of #pi^{-}p_{1};cos#theta;Momentum (GeV/c)", 2000, -1.0, 1.0, 2000, 0.0, 2.0);
  new TH2F("pip1_CosTvsPhi_CM","cos#theta vs. #phi of #pi^{-}p_{1};cos#theta;#phi", 2000, -1.0, 1.0, 1600, 0.0, 3.2);
  new TH2F("HeKpip1_MMvsMomentum_CM", "^{3}He(K^{-},#pi^{-}p_{1})X missing mass vs. missing momentum;MM(#pi^{-}p_{1}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 1.0, 3.0, 2000, 0.0, 2.0);
  new TH2F("HeKpip1_MMvsCosT_CM", "^{3}He(K^{-},#pi^{-}p_{1})X missing mass vs. cos#theta;MM(#pi^{-}p_{1}) (GeV/c^{2});cos#theta", 2000, 1.0, 3.0, 2000, -1.0, 1.0);
  new TH2F("pip1_IMvsMomentum_CM", "(#pi^{-}p_{1}) invariant mass vs. momentum;IM(#pi^{-}p_{1}) (GeV/c^{2});Momentum (GeV/c)", 2000, 1.0, 3.0, 2000, 0.0, 2.0);
  new TH2F("pip1_IMvsCosT_CM", "(#pi^{-}p_{1}) invariant mass vs. cos#theta;IM(#pi^{-}p_{1}) (GeV/c^{2});cos#theta", 2000, 1.0, 3.0, 2000, -1.0, 1.0);
  // pi- proton2
  new TH2F("pip2_CosTvsMomentum_CM","cos#theta vs. momentum of #pi^{-}p_{2};cos#theta;Momentum (GeV/c)", 2000, -1.0, 1.0, 2000, 0.0, 2.0);
  new TH2F("pip2_CosTvsPhi_CM","cos#theta vs. #phi of #pi^{-}p_{2};cos#theta;#phi", 2000, -1.0, 1.0, 1600, 0.0, 3.2);
  new TH2F("HeKpip2_MMvsMomentum_CM", "^{3}He(K^{-},#pi^{-}p_{2})X missing mass vs. missing momentum;MM(#pi^{-}p_{2}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 1.0, 3.0, 2000, 0.0, 2.0);
  new TH2F("HeKpip2_MMvsCosT_CM", "^{3}He(K^{-},#pi^{-}p_{2})X missing mass vs. cos#theta;MM(#pi^{-}p_{2}) (GeV/c^{2});cos#theta", 2000, 1.0, 3.0, 2000, -1.0, 1.0);
  new TH2F("pip2_IMvsMomentum_CM", "(#pi^{-}p_{2}) invariant mass vs. momentum;IM(#pi^{-}p_{2}) (GeV/c^{2});Momentum (GeV/c)", 2000, 1.0, 4.0, 2000, 0.0, 2.0);
  new TH2F("pip2_IMvsCosT_CM", "(#pi^{-}p_{2}) invariant mass vs. cos#theta;IM(#pi^{-}p_{2}) (GeV/c^{2});cos#theta", 2000, 1.0, 3.0, 2000, -1.0, 1.0);
  // pi- proton1 proton2
  new TH2F("pipp_CosTvsMomentum_CM","cos#theta vs. momentum of #pi^{-}p_{1}p_{2};cos#theta;Momentum (GeV/c)", 2000, -1.0, 1.0, 2000, 0.0, 2.0);
  new TH2F("pipp_CosTvsPhi_CM","cos#theta vs. #phi of #pi^{-}p_{1}p_{2};cos#theta;#phi", 2000, -1.0, 1.0, 1600, 0.0, 3.2);
  new TH2F("HeKpipp_MMvsMomentum_CM", "^{3}He(K^{-},#pi^{-}p_{1}p_{2})X missing mass vs. missing momentum;MM(#pi^{-}p_{1}p_{2}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("HeKpipp_MMvsCosT_CM", "^{3}He(K^{-},#pi^{-}p_{1}p_{2})X missing mass vs. cos#theta;MM(#pi^{-}p_{1}p_{2}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pipp_IMvsMomentum_CM", "(#pi^{-}p_{1}p_{2}) invariant mass vs. momentum;IM(#pi^{-}p_{1}p_{2}) (GeV/c^{2});Momentum (GeV/c)", 2000, 2.0, 4.0, 2000, 0.0, 2.0);
  new TH2F("pipp_IMvsCosT_CM", "(#pi^{-}p_{1}p_{2}) invariant mass vs. cos#theta;IM(#pi^{-}p_{1}p_{2}) (GeV/c^{2});cos#theta", 2000, 2.0, 4.0, 2000, -1.0, 1.0);

  new TH2F("pip1_IMvspip2_IM_CM", "(#pi^{-}p_{1}) invariant mass vs. (#pi^{-}p_{2}) invariant mass;IM(#pi^{-}p_{1}) (GeV/c^{2});IM(#pi^{-}p_{2}) (GeV/c^{2})", 2000, 1.0, 3.0, 2000, 1.0, 3.0);
  new TH2F("pipp_IMvsHeKpipp_MM_CM", "(#pi^{-}p_{1}p_{2}) invariant mass vs. ^{3}He(K^{-},#pi^{-}p_{1}p_{2})X missing mass;IM(#pi^{-}p_{1}p_{2}) (GeV/c^{2});MM(#pi^{-}p_{1}p_{2}) (GeV/c^{2})", 2000, 2.0, 4.0, 2000, 0.0, 2.0);

  std::cout << "== Finish Histogram Initialization ==" << std::endl;
  return true;

}
