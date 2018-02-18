// MyAnalysisHeKpipin.cpp

#include "MyAnalysisHeKpipin.h"

MyAnalysisHeKpipin::MyAnalysisHeKpipin(TFile* rt, ConfMan* conf)
{
  Initialize(conf);
  CutCondition();
  Clear();
}

MyAnalysisHeKpipin::~MyAnalysisHeKpipin()
{
  Clear();
  rtFile->cd();
  rtFile->Write();
  rtFile->Close();
}

void MyAnalysisHeKpipin::Clear()
{
}

bool MyAnalysisHeKpipin::DoAnalysis(ConfMan* conf, EventHeader* header, BeamLineHitMan* blMan, BeamLineTrackMan* bltrackMan, CDSHitMan* cdsMan, CDSTrackingMan* cdstrackMan, Particle* particle)
{
  rtFile->cd();

  int ievent=0;
  FillHist("EventNumber",ievent); ievent++; /* All Events */
  DetectorList *dlist=DetectorList::GetInstance();
  if(particle->nBeam()!=1) return false;
  pBeam* beam = particle->beam(0);
  double beamtof = beam->bhdt0tof();

  /* Event selection */
  if(particle->nCDS()==0) return true;
  FillHist("EventNumber",ievent); ievent++; /* Least 1 hit in CDS*/
  bool neutral=false;
  if(particle->nNC()>0) neutral = true;
  if(!neutral) return true;
  FillHist("EventNumber",ievent); ievent++; /* Neutral Trigger */
  if(particle->nPiminus()!=1) return true;
  if(particle->nPiplus()!=1) return true;
  FillHist("EventNumber",ievent); ievent++; /* pi+ and pi- in CDS */
  pCDS* pip=particle->pip(0);
  pCDS* pim=particle->pim(0);

  // ############### //
  // Vertex decision //
  // ############### //
  TVector3 vertex;
  double vdis = 9999;
  int fcds = -1;
  for(int it=0; it<particle->nCDS(); it++){
    pCDS* cds = particle->cdsi(it);
    if(cds->chi()>30) continue;
    if(vdis>9998||vdis>cds->vdis()){
      vdis = cds->vdis();
      vertex = cds->vbeam();
      if(GeomTools::GetID(vertex)==CID_Fiducial) fcds=it;
    }
  }
  if(fcds==-1) return true;
  FillHist("EventNumber",ievent); ievent++; /* Vertex decision */

  if(particle->nProton()==1){
    DoOneProtonAna(particle,beam,cds1,vertex);
  }
  if(particle->nProton()==2){
    DoTwoProtonAna(particle,beam,cds1,cds2,vertex);
  }

  return true;

}

bool MyAnalysisHeKpipin::DoOneProtonAna(Particle* particle, pBeam* beam, pCDS* cds, TVector3 vertex)
{
  TString mode = "1p";

  // ############### //
  // NC hit decision //
  // ############### //
  double time = 9999;
  int fnc = -1;
  for(int it=0; it<particle->nNC(); it++){
    pNC* nc = particle->nc(it);
    nc->CalcMom(beam,vertex);
    double seg = (nc->seg()-1)%16+1;
    double lay = (nc->seg()-1)/16+1;
    FillHist("FWDN_Overbeta",1.0/nc->beta());
    FillHist("FWDN_OverbetavsEnergy",1.0/nc->beta(),nc->energy());
    FillHist("NCHitPosition_XY",nc->hitpos().X(),nc->hitpos().Y());
    FillHist("NCHitSegment",seg,lay);
    if(nc->energy()>8.0 && (time>9998||time>nc->time())){
      time = nc->time();
      fnc = it;
      FillHist("FWDN_Overbeta_withdECut",1.0/nc->beta());
      FillHist("FWDN_OverbetavsEnergy_withdECut",1.0/nc->beta(),nc->energy());
    }
  }
  if(fnc==-1) return true;
  pNC*  nc =0;
  if(fnc!=-1) nc = particle->nc(fnc);
  FillHist("FWDN_Overbeta_Selected",1.0/nc->beta());
  if(nc->pid()!=F_Neutron){ return true; }

  TString frame[2] = {"Lab","CM"};

  // Target
  TVector3 ZeroV;
  TLorentzVector Ltgt; Ltgt.SetVectM(ZeroV,ThreeHeMass);
  TLorentzVector Lbeam = beam->GetLorentzVector(vertex);
  TVector3 boost = (Ltgt+Lbeam).BoostVector();
  TLorentzVector cLtgt = Ltgt; cLtgt.Boost(-boost); TLorentzVector cLbeam = Lbeam; cLbeam.Boost(-boost);
  /* CDS */
  TLorentzVector Lcds = cds->GetLorentzVector();
  TLorentzVector cLcds = cds->GetLorentzVector(); cLcds.Boost(-boost);
  /* NC */
  TLorentzVector Ln  = nc ->GetLorentzVector();
  TLorentzVector cLn  = nc ->GetLorentzVector(); cLn.Boost(-boost);
  double nc_mass[2]    = { Ln.M()                           , cLn.M()};
  double nc_mom[2]     = { Ln.Vect().Mag()                  , cLn.Vect().Mag()};
  double nc_cost[2]    = { Ln.Vect().CosTheta()             , cLn.Vect().CosTheta()};
  double nc_phi[2]     = { Ln.Vect().Phi()                  , cLn.Vect().Phi()};
  double nc_mmass[2]   = { (Ltgt+Lbeam-Ln).M()              , (cLtgt+cLbeam-cLn).M()};
  double nc_mmom[2]    = { (Ltgt+Lbeam-Ln).Vect().Mag()     , (cLtgt+cLbeam-cLn).Vect().Mag()};
  double nc_mcost[2]   = { (Ltgt+Lbeam-Ln).Vect().CosTheta(), (cLtgt+cLbeam-cLn).Vect().CosTheta()};
  /* CDS */
  double cds_mass[2]  = { (Lcds).M()                           , (cLcds).M()};
  double cds_mom[2]   = { (Lcds).Vect().Mag()                  , (cLcds).Vect().Mag()};
  double cds_cost[2]  = { (Lcds).Vect().CosTheta()             , (cLcds).Vect().CosTheta()};
  double cds_phi[2]   = { (Lcds).Vect().Phi()                  , (cLcds).Vect().Phi()};
  double cds_mmass[2] = { (Ltgt+Lbeam-(Lcds)).M()              , (cLtgt+cLbeam-(cLcds)).M()};
  double cds_mmom[2]  = { (Ltgt+Lbeam-(Lcds)).Vect().Mag()     , (cLtgt+cLbeam-(cLcds)).Vect().Mag()};
  double cds_mcost[2] = { (Ltgt+Lbeam-(Lcds)).Vect().CosTheta(), (cLtgt+cLbeam-(cLcds)).Vect().CosTheta()};
  /* CDS + NC */
  double ncds_mass[2]  = { (Lcds+Ln).M()                           , (cLcds+cLn).M()};
  double ncds_mom[2]   = { (Lcds+Ln).Vect().Mag()                  , (cLcds+cLn).Vect().Mag()};
  double ncds_cost[2]  = { (Lcds+Ln).Vect().CosTheta()             , (cLcds+cLn).Vect().CosTheta()};
  double ncds_phi[2]   = { (Lcds+Ln).Vect().Phi()                  , (cLcds+cLn).Vect().Phi()};
  double ncds_mmass[2] = { (Ltgt+Lbeam-(Lcds+Ln)).M()              , (cLtgt+cLbeam-(cLcds+cLn)).M()};
  double ncds_mmom[2]  = { (Ltgt+Lbeam-(Lcds+Ln)).Vect().Mag()     , (cLtgt+cLbeam-(cLcds+cLn)).Vect().Mag()};
  double ncds_mcost[2] = { (Ltgt+Lbeam-(Lcds+Ln)).Vect().CosTheta(), (cLtgt+cLbeam-(cLcds+cLn)).Vect().CosTheta()};

  // n
  for(int f=1; f<2; f++){
    FillHist(Form("n_Mass_%s_%s",frame[f].Data(),mode.Data()),nc_mass[f]);
    FillHist(Form("n_Momentum_%s_%s",frame[f].Data(),mode.Data()),nc_mom[f]);
    FillHist(Form("n_CosT_%s_%s",frame[f].Data(),mode.Data()),nc_cost[f]);
    FillHist(Form("n_Phi_%s_%s",frame[f].Data(),mode.Data()),nc_phi[f]);
    FillHist(Form("n_MMass_%s_%s",frame[f].Data(),mode.Data()),nc_mmass[f]);
    FillHist(Form("n_MMomentum_%s_%s",frame[f].Data(),mode.Data()),nc_mmom[f]);
    FillHist(Form("n_MCosT_%s_%s",frame[f].Data(),mode.Data()),nc_mcost[f]);
  }

  // proton
  for(int f=1; f<2; f++){
    FillHist(Form("p_Mass_%s_%s",frame[f].Data(),mode.Data()),cds_mass[f]);
    FillHist(Form("p_Momentum_%s_%s",frame[f].Data(),mode.Data()),cds_mom[f]);
    FillHist(Form("p_CosT_%s_%s",frame[f].Data(),mode.Data()),cds_cost[f]);
    FillHist(Form("p_Phi_%s_%s",frame[f].Data(),mode.Data()),cds_phi[f]);
    FillHist(Form("p_MMass_%s_%s",frame[f].Data(),mode.Data()),cds_mmass[f]);
    FillHist(Form("p_MMomentum_%s_%s",frame[f].Data(),mode.Data()),cds_mmom[f]);
    FillHist(Form("p_MCosT_%s_%s",frame[f].Data(),mode.Data()),cds_mcost[f]);
  }
  // n + proton
  for(int f=1; f<2; f++){
    FillHist(Form("np_Mass_%s_%s",frame[f].Data(),mode.Data()),ncds_mass[f]);
    FillHist(Form("np_Momentum_%s_%s",frame[f].Data(),mode.Data()),ncds_mom[f]);
    FillHist(Form("np_CosT_%s_%s",frame[f].Data(),mode.Data()),ncds_cost[f]);
    FillHist(Form("np_Phi_%s_%s",frame[f].Data(),mode.Data()),ncds_phi[f]);
    FillHist(Form("np_MMass_%s_%s",frame[f].Data(),mode.Data()),ncds_mmass[f]);
    FillHist(Form("np_MMomentum_%s_%s",frame[f].Data(),mode.Data()),ncds_mmom[f]);
    FillHist(Form("np_MCosT_%s_%s",frame[f].Data(),mode.Data()),ncds_mcost[f]);
  }

  /* 2D plot */
  // MM(np) vs. MM(n)
  for(int f=1; f<2; f++){
    FillHist(Form("n_MMass_vs_np_MMass_%s_%s",frame[f].Data(),mode.Data()),nc_mmass[f],ncds_mmass[f]);
  }

  return true;

}

bool MyAnalysisHeKpipin::DoTwoProtonAna(Particle* particle, pBeam* beam, pCDS* cds1, pCDS* cds2, TVector3 vertex)
{
  TString mode = "2p";

  // ############### //
  // NC hit decision //
  // ############### //
  double time = 9999;
  int fnc = -1;
  for(int it=0; it<particle->nNC(); it++){
    pNC* nc = particle->nc(it);
    nc->CalcMom(beam,vertex);
    double seg = (nc->seg()-1)%16+1;
    double lay = (nc->seg()-1)/16+1;
    if(nc->energy()>8.0 && (time>9998||time>nc->time())){
      time = nc->time();
      fnc = it;
    }
  }
  if(fnc==-1) return true;
  pNC*  nc =0;
  if(fnc!=-1) nc = particle->nc(fnc);
  if(nc->pid()!=F_Neutron){ return true; }

  TString frame[2] = {"Lab","CM"};

  // Target
  TVector3 ZeroV;
  TLorentzVector Ltgt; Ltgt.SetVectM(ZeroV,ThreeHeMass);
  TLorentzVector Lbeam = beam->GetLorentzVector(vertex);
  TVector3 boost = (Ltgt+Lbeam).BoostVector();
  TLorentzVector cLtgt = Ltgt; cLtgt.Boost(-boost); TLorentzVector cLbeam = Lbeam; cLbeam.Boost(-boost);
  /* CDS */
  TLorentzVector Lcds1 = cds1->GetLorentzVector();
  TLorentzVector cLcds1 = cds1->GetLorentzVector(); cLcds1.Boost(-boost);
  TLorentzVector Lcds2 = cds2->GetLorentzVector();
  TLorentzVector cLcds2 = cds2->GetLorentzVector(); cLcds2.Boost(-boost);
  /* NC */
  TLorentzVector Ln  = nc ->GetLorentzVector();
  TLorentzVector cLn  = nc ->GetLorentzVector(); cLn.Boost(-boost);
  double nc_mass[2]    = { Ln.M()                           , cLn.M()};
  double nc_mom[2]     = { Ln.Vect().Mag()                  , cLn.Vect().Mag()};
  double nc_cost[2]    = { Ln.Vect().CosTheta()             , cLn.Vect().CosTheta()};
  double nc_phi[2]     = { Ln.Vect().Phi()                  , cLn.Vect().Phi()};
  double nc_mmass[2]   = { (Ltgt+Lbeam-Ln).M()              , (cLtgt+cLbeam-cLn).M()};
  double nc_mmom[2]    = { (Ltgt+Lbeam-Ln).Vect().Mag()     , (cLtgt+cLbeam-cLn).Vect().Mag()};
  double nc_mcost[2]   = { (Ltgt+Lbeam-Ln).Vect().CosTheta(), (cLtgt+cLbeam-cLn).Vect().CosTheta()};
  /* CDS1 */
  double cds1_mass[2]  = { (Lcds1).M()                           , (cLcds1).M()};
  double cds1_mom[2]   = { (Lcds1).Vect().Mag()                  , (cLcds1).Vect().Mag()};
  double cds1_cost[2]  = { (Lcds1).Vect().CosTheta()             , (cLcds1).Vect().CosTheta()};
  double cds1_phi[2]   = { (Lcds1).Vect().Phi()                  , (cLcds1).Vect().Phi()};
  double cds1_mmass[2] = { (Ltgt+Lbeam-(Lcds1)).M()              , (cLtgt+cLbeam-(cLcds1)).M()};
  double cds1_mmom[2]  = { (Ltgt+Lbeam-(Lcds1)).Vect().Mag()     , (cLtgt+cLbeam-(cLcds1)).Vect().Mag()};
  double cds1_mcost[2] = { (Ltgt+Lbeam-(Lcds1)).Vect().CosTheta(), (cLtgt+cLbeam-(cLcds1)).Vect().CosTheta()};
  /* CDS1 + NC */
  double ncds1_mass[2]  = { (Lcds1+Ln).M()                           , (cLcds1+cLn).M()};
  double ncds1_mom[2]   = { (Lcds1+Ln).Vect().Mag()                  , (cLcds1+cLn).Vect().Mag()};
  double ncds1_cost[2]  = { (Lcds1+Ln).Vect().CosTheta()             , (cLcds1+cLn).Vect().CosTheta()};
  double ncds1_phi[2]   = { (Lcds1+Ln).Vect().Phi()                  , (cLcds1+cLn).Vect().Phi()};
  double ncds1_mmass[2] = { (Ltgt+Lbeam-(Lcds1+Ln)).M()              , (cLtgt+cLbeam-(cLcds1+cLn)).M()};
  double ncds1_mmom[2]  = { (Ltgt+Lbeam-(Lcds1+Ln)).Vect().Mag()     , (cLtgt+cLbeam-(cLcds1+cLn)).Vect().Mag()};
  double ncds1_mcost[2] = { (Ltgt+Lbeam-(Lcds1+Ln)).Vect().CosTheta(), (cLtgt+cLbeam-(cLcds1+cLn)).Vect().CosTheta()};
  /* CDS2 */
  double cds2_mass[2]  = { (Lcds2).M()                           , (cLcds2).M()};
  double cds2_mom[2]   = { (Lcds2).Vect().Mag()                  , (cLcds2).Vect().Mag()};
  double cds2_cost[2]  = { (Lcds2).Vect().CosTheta()             , (cLcds2).Vect().CosTheta()};
  double cds2_phi[2]   = { (Lcds2).Vect().Phi()                  , (cLcds2).Vect().Phi()};
  double cds2_mmass[2] = { (Ltgt+Lbeam-(Lcds2)).M()              , (cLtgt+cLbeam-(cLcds2)).M()};
  double cds2_mmom[2]  = { (Ltgt+Lbeam-(Lcds2)).Vect().Mag()     , (cLtgt+cLbeam-(cLcds2)).Vect().Mag()};
  double cds2_mcost[2] = { (Ltgt+Lbeam-(Lcds2)).Vect().CosTheta(), (cLtgt+cLbeam-(cLcds2)).Vect().CosTheta()};
  /* CDS2 + NC */
  double ncds2_mass[2]  = { (Lcds2+Ln).M()                           , (cLcds2+cLn).M()};
  double ncds2_mom[2]   = { (Lcds2+Ln).Vect().Mag()                  , (cLcds2+cLn).Vect().Mag()};
  double ncds2_cost[2]  = { (Lcds2+Ln).Vect().CosTheta()             , (cLcds2+cLn).Vect().CosTheta()};
  double ncds2_phi[2]   = { (Lcds2+Ln).Vect().Phi()                  , (cLcds2+cLn).Vect().Phi()};
  double ncds2_mmass[2] = { (Ltgt+Lbeam-(Lcds2+Ln)).M()              , (cLtgt+cLbeam-(cLcds2+cLn)).M()};
  double ncds2_mmom[2]  = { (Ltgt+Lbeam-(Lcds2+Ln)).Vect().Mag()     , (cLtgt+cLbeam-(cLcds2+cLn)).Vect().Mag()};
  double ncds2_mcost[2] = { (Ltgt+Lbeam-(Lcds2+Ln)).Vect().CosTheta(), (cLtgt+cLbeam-(cLcds2+cLn)).Vect().CosTheta()};
  /* CDS */
  double cds_mass[2]  = { (Lcds1+Lcds2).M()                           , (cLcds1+cLcds2).M()};
  double cds_mom[2]   = { (Lcds1+Lcds2).Vect().Mag()                  , (cLcds1+cLcds2).Vect().Mag()};
  double cds_cost[2]  = { (Lcds1+Lcds2).Vect().CosTheta()             , (cLcds1+cLcds2).Vect().CosTheta()};
  double cds_phi[2]   = { (Lcds1+Lcds2).Vect().Phi()                  , (cLcds1+cLcds2).Vect().Phi()};
  double cds_mmass[2] = { (Ltgt+Lbeam-(Lcds1+Lcds2)).M()              , (cLtgt+cLbeam-(cLcds1+cLcds2)).M()};
  double cds_mmom[2]  = { (Ltgt+Lbeam-(Lcds1+Lcds2)).Vect().Mag()     , (cLtgt+cLbeam-(cLcds1+cLcds2)).Vect().Mag()};
  double cds_mcost[2] = { (Ltgt+Lbeam-(Lcds1+Lcds2)).Vect().CosTheta(), (cLtgt+cLbeam-(cLcds1+cLcds2)).Vect().CosTheta()};
  /* CDS + NC */
  double ncds_mass[2]  = { (Lcds1+Lcds2+Ln).M()                           , (cLcds1+cLcds2+cLn).M()};
  double ncds_mom[2]   = { (Lcds1+Lcds2+Ln).Vect().Mag()                  , (cLcds1+cLcds2+cLn).Vect().Mag()};
  double ncds_cost[2]  = { (Lcds1+Lcds2+Ln).Vect().CosTheta()             , (cLcds1+cLcds2+cLn).Vect().CosTheta()};
  double ncds_phi[2]   = { (Lcds1+Lcds2+Ln).Vect().Phi()                  , (cLcds1+cLcds2+cLn).Vect().Phi()};
  double ncds_mmass[2] = { (Ltgt+Lbeam-(Lcds1+Lcds2+Ln)).M()              , (cLtgt+cLbeam-(cLcds1+cLcds2+cLn)).M()};
  double ncds_mmom[2]  = { (Ltgt+Lbeam-(Lcds1+Lcds2+Ln)).Vect().Mag()     , (cLtgt+cLbeam-(cLcds1+cLcds2+cLn)).Vect().Mag()};
  double ncds_mcost[2] = { (Ltgt+Lbeam-(Lcds1+Lcds2+Ln)).Vect().CosTheta(), (cLtgt+cLbeam-(cLcds1+cLcds2+cLn)).Vect().CosTheta()};

  // n
  for(int f=1; f<2; f++){
    FillHist(Form("n_Mass_%s_%s",frame[f].Data(),mode.Data()),nc_mass[f]);
    FillHist(Form("n_Momentum_%s_%s",frame[f].Data(),mode.Data()),nc_mom[f]);
    FillHist(Form("n_CosT_%s_%s",frame[f].Data(),mode.Data()),nc_cost[f]);
    FillHist(Form("n_Phi_%s_%s",frame[f].Data(),mode.Data()),nc_phi[f]);
    FillHist(Form("n_MMass_%s_%s",frame[f].Data(),mode.Data()),nc_mmass[f]);
    FillHist(Form("n_MMomentum_%s_%s",frame[f].Data(),mode.Data()),nc_mmom[f]);
    FillHist(Form("n_MCosT_%s_%s",frame[f].Data(),mode.Data()),nc_mcost[f]);
  }

  // proton1
  for(int f=1; f<2; f++){
    FillHist(Form("p1_Mass_%s_%s",frame[f].Data(),mode.Data()),cds1_mass[f]);
    FillHist(Form("p1_Momentum_%s_%s",frame[f].Data(),mode.Data()),cds1_mom[f]);
    FillHist(Form("p1_CosT_%s_%s",frame[f].Data(),mode.Data()),cds1_cost[f]);
    FillHist(Form("p1_Phi_%s_%s",frame[f].Data(),mode.Data()),cds1_phi[f]);
    FillHist(Form("p1_MMass_%s_%s",frame[f].Data(),mode.Data()),cds1_mmass[f]);
    FillHist(Form("p1_MMomentum_%s_%s",frame[f].Data(),mode.Data()),cds1_mmom[f]);
    FillHist(Form("p1_MCosT_%s_%s",frame[f].Data(),mode.Data()),cds1_mcost[f]);
  }
  // n + proton1
  for(int f=1; f<2; f++){
    FillHist(Form("np1_Mass_%s_%s",frame[f].Data(),mode.Data()),ncds1_mass[f]);
    FillHist(Form("np1_Momentum_%s_%s",frame[f].Data(),mode.Data()),ncds1_mom[f]);
    FillHist(Form("np1_CosT_%s_%s",frame[f].Data(),mode.Data()),ncds1_cost[f]);
    FillHist(Form("np1_Phi_%s_%s",frame[f].Data(),mode.Data()),ncds1_phi[f]);
    FillHist(Form("np1_MMass_%s_%s",frame[f].Data(),mode.Data()),ncds1_mmass[f]);
    FillHist(Form("np1_MMomentum_%s_%s",frame[f].Data(),mode.Data()),ncds1_mmom[f]);
    FillHist(Form("np1_MCosT_%s_%s",frame[f].Data(),mode.Data()),ncds1_mcost[f]);
  }

  // proton2
  for(int f=1; f<2; f++){
    FillHist(Form("p2_Mass_%s_%s",frame[f].Data(),mode.Data()),cds2_mass[f]);
    FillHist(Form("p2_Momentum_%s_%s",frame[f].Data(),mode.Data()),cds2_mom[f]);
    FillHist(Form("p2_CosT_%s_%s",frame[f].Data(),mode.Data()),cds2_cost[f]);
    FillHist(Form("p2_Phi_%s_%s",frame[f].Data(),mode.Data()),cds2_phi[f]);
    FillHist(Form("p2_MMass_%s_%s",frame[f].Data(),mode.Data()),cds2_mmass[f]);
    FillHist(Form("p2_MMomentum_%s_%s",frame[f].Data(),mode.Data()),cds2_mmom[f]);
    FillHist(Form("p2_MCosT_%s_%s",frame[f].Data(),mode.Data()),cds2_mcost[f]);
  }
  // n + proton2
  for(int f=1; f<2; f++){
    FillHist(Form("np2_Mass_%s_%s",frame[f].Data(),mode.Data()),ncds2_mass[f]);
    FillHist(Form("np2_Momentum_%s_%s",frame[f].Data(),mode.Data()),ncds2_mom[f]);
    FillHist(Form("np2_CosT_%s_%s",frame[f].Data(),mode.Data()),ncds2_cost[f]);
    FillHist(Form("np2_Phi_%s_%s",frame[f].Data(),mode.Data()),ncds2_phi[f]);
    FillHist(Form("np2_MMass_%s_%s",frame[f].Data(),mode.Data()),ncds2_mmass[f]);
    FillHist(Form("np2_MMomentum_%s_%s",frame[f].Data(),mode.Data()),ncds2_mmom[f]);
    FillHist(Form("np2_MCosT_%s_%s",frame[f].Data(),mode.Data()),ncds2_mcost[f]);
  }

  // proton1 + proton2
  for(int f=1; f<2; f++){
    FillHist(Form("pp_Mass_%s_%s",frame[f].Data(),mode.Data()),cds_mass[f]);
    FillHist(Form("pp_Momentum_%s_%s",frame[f].Data(),mode.Data()),cds_mom[f]);
    FillHist(Form("pp_CosT_%s_%s",frame[f].Data(),mode.Data()),cds_cost[f]);
    FillHist(Form("pp_Phi_%s_%s",frame[f].Data(),mode.Data()),cds_phi[f]);
    FillHist(Form("pp_MMass_%s_%s",frame[f].Data(),mode.Data()),cds_mmass[f]);
    FillHist(Form("pp_MMomentum_%s_%s",frame[f].Data(),mode.Data()),cds_mmom[f]);
    FillHist(Form("pp_MCosT_%s_%s",frame[f].Data(),mode.Data()),cds_mcost[f]);
  }
  // n + proton1 + proton2
  for(int f=1; f<2; f++){
    FillHist(Form("npp_Mass_%s_%s",frame[f].Data(),mode.Data()),ncds_mass[f]);
    FillHist(Form("npp_Momentum_%s_%s",frame[f].Data(),mode.Data()),ncds_mom[f]);
    FillHist(Form("npp_CosT_%s_%s",frame[f].Data(),mode.Data()),ncds_cost[f]);
    FillHist(Form("npp_Phi_%s_%s",frame[f].Data(),mode.Data()),ncds_phi[f]);
    FillHist(Form("npp_MMass_%s_%s",frame[f].Data(),mode.Data()),ncds_mmass[f]);
    FillHist(Form("npp_MMomentum_%s_%s",frame[f].Data(),mode.Data()),ncds_mmom[f]);
    FillHist(Form("npp_MCosT_%s_%s",frame[f].Data(),mode.Data()),ncds_mcost[f]);
  }


  /* 2D plot */
  // MM(npp) vs. MM(n)
  for(int f=1; f<2; f++){
    FillHist(Form("n_MMass_vs_npp_MMass_%s_%s",frame[f].Data(),mode.Data()),nc_mmass[f],ncds_mmass[f]);
  }
  // MM(npp) vs. MM(np1)
  for(int f=1; f<2; f++){
    FillHist(Form("np1_MMass_vs_npp_MMass_%s_%s",frame[f].Data(),mode.Data()),ncds1_mmass[f],ncds_mmass[f]);
  }
  // MM(npp) vs. MM(np2)
  for(int f=1; f<2; f++){
    FillHist(Form("np2_MMass_vs_npp_MMass_%s_%s",frame[f].Data(),mode.Data()),ncds2_mmass[f],ncds_mmass[f]);
  }
  // MM(np1) vs. MM(n)
  for(int f=1; f<2; f++){
    FillHist(Form("n_MMass_vs_np1_MMass_%s_%s",frame[f].Data(),mode.Data()),nc_mmass[f],ncds1_mmass[f]);
  }
  // MM(np2) vs. MM(n)
  for(int f=1; f<2; f++){
    FillHist(Form("n_MMass_vs_np2_MMass_%s_%s",frame[f].Data(),mode.Data()),nc_mmass[f],ncds2_mmass[f]);
  }
  
  return true;

}


bool MyAnalysisHeKpipin::FillHist(TString name, double val1, int weight)
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

bool MyAnalysisHeKpipin::FillHist(TString name, TString val1, int weight)
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

bool MyAnalysisHeKpipin::FillHist(TString name, double val1, double val2, int weight)
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

bool MyAnalysisHeKpipin::FillHist(TString name, TString val1, TString val2, int weight)
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

void MyAnalysisHeKpipin::CutCondition()
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

bool MyAnalysisHeKpipin::Initialize(ConfMan* confMan)
{
  std::cout << "### MyAnalysisHeKpipin::Initialize ###" << std::endl;

  std::string ofname = confMan->GetOutFileName();
  ofname.insert(ofname.find(".root"),"_anaHeKpipin");

  rtFile =  new TFile( Form("%s",ofname.c_str()), "RECREATE");
  rtFile -> cd();

  new TH1F( "EventNumber", "Number of Events", 20, 0, 20);
  // FWD neutral Particles
  std::cout << "Define Histograms for FWD neutral particles" << std::endl;
  new TH1F( "FWDN_Overbeta", "1/#beta;1/#beta;Counts", 5000, 0, 5);
  new TH1F( "FWDN_Overbeta_withdECut", "1/#beta;1/#beta;Counts", 5000, 0, 5);
  new TH1F( "FWDN_Overbeta_Selected", "1/#beta;1/#beta;Counts", 5000, 0, 5);
  new TH2F( "FWDN_OverbetavsEnergy", "1/#beta vs. Energy deposit;1/#beta;Energy deposit (MeVee)", 1000, -0, 10, 5000, 0, 100);
  new TH2F( "FWDN_OverbetavsEnergy_withdECut", "1/#beta vs. Energy deposit;1/#beta;Energy deposit (MeVee)", 1000, -0, 10, 5000, 0, 100);
  new TH2F( "NCHitPosition_XY", "Hit position at NC (X vs. Y);X position (cm);Y position (cm)",16,-160,160,15,-75,75);
  new TH2F( "NCHitSegment", "NC segment;Segment in layer;Layer",16,0.5,16.5,7,0.5,7.5);

  std::cout << "Define Histograms for 1p coin." << std::endl;
  // n
  new TH1F("n_Mass_CM_1p","Invariant mass of n;IM(n) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("n_Momentum_CM_1p","Momentum of n;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("n_CosT_CM_1p","cos#theta of n;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("n_Phi_CM_1p","#phi of n;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("n_MMass_CM_1p", "^{3}He(K^{-},n)X missing mass;MM(n) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("n_MMomentum_CM_1p", "^{3}He(K^{-},n)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("n_MCosT_CM_1p", "^{3}He(K^{-},n)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // proton
  new TH1F("p_Mass_CM_1p","Invariant mass of p;IM(p) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("p_Momentum_CM_1p","Momentum of p;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("p_CosT_CM_1p","cos#theta of p;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("p_Phi_CM_1p","#phi of p;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("p_MMass_CM_1p", "^{3}He(K^{-},p)X missing mass;MM(p) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("p_MMomentum_CM_1p", "^{3}He(K^{-},p)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("p_MCosT_CM_1p", "^{3}He(K^{-},p)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // n + proton
  new TH1F("np_Mass_CM_1p","Invariant mass of np;IM(np) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("np_Momentum_CM_1p","Momentum of np;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("np_CosT_CM_1p","cos#theta of np;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("np_Phi_CM_1p","#phi of np;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("np_MMass_CM_1p", "^{3}He(K^{-},np)X missing mass;MM(np) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("np_MMomentum_CM_1p", "^{3}He(K^{-},np)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("np_MCosT_CM_1p", "^{3}He(K^{-},np)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);

  new TH2F("n_MMass_vs_np_MMass_CM_1p", "^{3}He(K^{-},n)X missing mass vs. ^{3}He(K^{-},np)X missing mass;MM(n) (GeV/c^{2});MM(np) (GeV/c^{2})", 1000, 2.0, 3.0, 1000, 1.0, 2.0);

  std::cout << "Define Histograms for 2p coin." << std::endl;
  // n
  new TH1F("n_Mass_CM_2p","Invariant mass of n;IM(n) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("n_Momentum_CM_2p","Momentum of n;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("n_CosT_CM_2p","cos#theta of n;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("n_Phi_CM_2p","#phi of n;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("n_MMass_CM_2p", "^{3}He(K^{-},n)X missing mass;MM(n) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("n_MMomentum_CM_2p", "^{3}He(K^{-},n)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("n_MCosT_CM_2p", "^{3}He(K^{-},n)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // proton1
  new TH1F("p1_Mass_CM_2p","Invariant mass of p1;IM(p1) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("p1_Momentum_CM_2p","Momentum of p1;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("p1_CosT_CM_2p","cos#theta of p1;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("p1_Phi_CM_2p","#phi of p1;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("p1_MMass_CM_2p", "^{3}He(K^{-},p1)X missing mass;MM(p1) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("p1_MMomentum_CM_2p", "^{3}He(K^{-},p1)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("p1_MCosT_CM_2p", "^{3}He(K^{-},p1)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // proton2
  new TH1F("p2_Mass_CM_2p","Invariant mass of p2;IM(p2) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("p2_Momentum_CM_2p","Momentum of p2;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("p2_CosT_CM_2p","cos#theta of p2;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("p2_Phi_CM_2p","#phi of p2;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("p2_MMass_CM_2p", "^{3}He(K^{-},p2)X missing mass;MM(p2) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("p2_MMomentum_CM_2p", "^{3}He(K^{-},p2)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("p2_MCosT_CM_2p", "^{3}He(K^{-},p2)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // proton1 + proton2
  new TH1F("pp_Mass_CM_2p","Invariant mass of pp;IM(pp) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pp_Momentum_CM_2p","Momentum of pp;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pp_CosT_CM_2p","cos#theta of pp;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("pp_Phi_CM_2p","#phi of pp;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("pp_MMass_CM_2p", "^{3}He(K^{-},pp)X missing mass;MM(pp) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pp_MMomentum_CM_2p", "^{3}He(K^{-},pp)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pp_MCosT_CM_2p", "^{3}He(K^{-},pp)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // n + proton1
  new TH1F("np1_Mass_CM_2p","Invariant mass of np1;IM(np1) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("np1_Momentum_CM_2p","Momentum of np1;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("np1_CosT_CM_2p","cos#theta of np1;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("np1_Phi_CM_2p","#phi of np1;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("np1_MMass_CM_2p", "^{3}He(K^{-},np1)X missing mass;MM(np1) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("np1_MMomentum_CM_2p", "^{3}He(K^{-},np1)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("np1_MCosT_CM_2p", "^{3}He(K^{-},np1)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // n + proton2
  new TH1F("np2_Mass_CM_2p","Invariant mass of np2;IM(np2) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("np2_Momentum_CM_2p","Momentum of np2;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("np2_CosT_CM_2p","cos#theta of np2;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("np2_Phi_CM_2p","#phi of np2;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("np2_MMass_CM_2p", "^{3}He(K^{-},np2)X missing mass;MM(np2) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("np2_MMomentum_CM_2p", "^{3}He(K^{-},np2)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("np2_MCosT_CM_2p", "^{3}He(K^{-},np2)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // n + proton1 + proton2
  new TH1F("npp_Mass_CM_2p","Invariant mass of npp;IM(npp) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("npp_Momentum_CM_2p","Momentum of npp;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("npp_CosT_CM_2p","cos#theta of npp;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("npp_Phi_CM_2p","#phi of npp;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("npp_MMass_CM_2p", "^{3}He(K^{-},npp)X missing mass;MM(npp) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("npp_MMomentum_CM_2p", "^{3}He(K^{-},npp)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("npp_MCosT_CM_2p", "^{3}He(K^{-},npp)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);

  new TH2F("n_MMass_vs_npp_MMass_CM_2p", "^{3}He(K^{-},n)X missing mass vs. ^{3}He(K^{-},npp)X missing mass;MM(n) (GeV/c^{2});MM(npp) (GeV/c^{2})", 1000, 2.0, 3.0, 1000, 0.0, 1.0);
  new TH2F("np1_MMass_vs_npp_MMass_CM_2p", "^{3}He(K^{-},np1)X missing mass vs. ^{3}He(K^{-},npp)X missing mass;MM(np1) (GeV/c^{2});MM(npp) (GeV/c^{2})", 1000, 1.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("np2_MMass_vs_npp_MMass_CM_2p", "^{3}He(K^{-},np2)X missing mass vs. ^{3}He(K^{-},npp)X missing mass;MM(np2) (GeV/c^{2});MM(npp) (GeV/c^{2})", 1000, 1.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("n_MMass_vs_np1_MMass_CM_2p", "^{3}He(K^{-},n)X missing mass vs. ^{3}He(K^{-},np1)X missing mass;MM(n) (GeV/c^{2});MM(np1) (GeV/c^{2})", 1000, 2.0, 3.0, 1000, 1.0, 2.0);
  new TH2F("n_MMass_vs_np2_MMass_CM_2p", "^{3}He(K^{-},n)X missing mass vs. ^{3}He(K^{-},np2)X missing mass;MM(n) (GeV/c^{2});MM(np2) (GeV/c^{2})", 1000, 2.0, 3.0, 1000, 1.0, 2.0);

  std::cout << "=== End of [pipin::Initialize] === " << std::endl;

  return true;
}
