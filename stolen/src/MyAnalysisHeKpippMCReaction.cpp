// MyAnalysisHeKpippMCReaction.cpp

#include "MyAnalysisHeKpippMCReaction.h"

MyAnalysisHeKpippMCReaction::MyAnalysisHeKpippMCReaction(TFile* rt, ConfMan* conf)
{
  Initialize(conf);
  CutCondition();
  Clear();
}

MyAnalysisHeKpippMCReaction::~MyAnalysisHeKpippMCReaction()
{
  Clear();
  rtFile->cd();
  rtFile->Write();
  rtFile->Close();
}

void MyAnalysisHeKpippMCReaction::Clear()
{
}

bool MyAnalysisHeKpippMCReaction::DoAnalysis(ConfMan* conf, EventHeader* header, BeamLineHitMan* blMan, BeamLineTrackMan* bltrackMan, CDSHitMan* cdsMan, CDSTrackingMan* cdstrackMan, Particle* particle, MCData* mcdata, ReactionData* reaction, DetectorData* detector)
{
  if(reaction==0) return true;

  DetectorList *dlist=DetectorList::GetInstance();
  rtFile->cd();

  int lamID = -1;
  for(int i=0; i<mcdata->trackSize(); i++){
    Track* track = mcdata->track(i);
    if(track->pdgID()==3122&&track->parentTrackID()==0){
      lamID = track->trackID();
    }
  }
  int lamseg=-1;
  int notlamseg=-1;
  for(int i=0; i<detector->detectorHitSize(); i++){
    DetectorHit* hit = detector->detectorHit(i);
    int pdg = hit->pdg();
    int parentID = hit->parentID();
    int detectorID = hit->detectorID();
    int layerID = hit->layerID();
    int channelID = hit->channelID();
    int seg = channelID+1;
    if(detectorID!=CID_CDH) continue;
    if(parentID==lamID&&pdg==2212){
      lamseg = seg;
    }
    if(parentID==0&&pdg==2212){
      notlamseg = seg;
    }
  }

  int size_cm = reaction->CMparticleSize();
  int size_lab = reaction->ParticleSize();
  int size_init = reaction->InitParticleSize();

  TLorentzVector Ltgt_mc;
  TLorentzVector Lbeam_mc;
  TLorentzVector Llambda_mc;
  TLorentzVector Lproton_mc;
  TLorentzVector Lneutron_mc;

  for(int i=0; i<reaction->InitParticleSize(); i++){
    if(reaction->InitPDG(i)==1000020030){
      Ltgt_mc = reaction->GetInitParticle(i);
      Ltgt_mc.SetPxPyPzE(Ltgt_mc.Px()/1000.0,Ltgt_mc.Py()/1000.0,Ltgt_mc.Pz()/1000.0,Ltgt_mc.E()/1000.0);
    }
    if(reaction->InitPDG(i)==-321){
      Lbeam_mc = reaction->GetInitParticle(i);
      Lbeam_mc.SetPxPyPzE(Lbeam_mc.Px()/1000.0,Lbeam_mc.Py()/1000.0,Lbeam_mc.Pz()/1000.0,Lbeam_mc.E()/1000.0);
    }
  }
  for(int i=0; i<reaction->ParticleSize(); i++){
    if(reaction->PDG(i)==3122){
      Llambda_mc = reaction->GetParticle(i);
      Llambda_mc.SetPxPyPzE(Llambda_mc.Px()/1000.0,Llambda_mc.Py()/1000.0,Llambda_mc.Pz()/1000.0,Llambda_mc.E()/1000.0);
    }
    if(reaction->PDG(i)==2212){
      Lproton_mc = reaction->GetParticle(i);
      Lproton_mc.SetPxPyPzE(Lproton_mc.Px()/1000.0,Lproton_mc.Py()/1000.0,Lproton_mc.Pz()/1000.0,Lproton_mc.E()/1000.0);
    }
    if(reaction->PDG(i)==2112){
      Lneutron_mc = reaction->GetParticle(i);
      Lneutron_mc.SetPxPyPzE(Lneutron_mc.Px()/1000.0,Lneutron_mc.Py()/1000.0,Lneutron_mc.Pz()/1000.0,Lneutron_mc.E()/1000.0);
    }
  }

  TVector3 boost_mc = (Ltgt_mc+Lbeam_mc).BoostVector();
  TLorentzVector cLtgt_mc = Ltgt_mc; cLtgt_mc.Boost(-boost_mc);
  TLorentzVector cLbeam_mc = Lbeam_mc; cLbeam_mc.Boost(-boost_mc);
  TLorentzVector cLlambda_mc = Llambda_mc; cLlambda_mc.Boost(-boost_mc);
  TLorentzVector cLproton_mc = Lproton_mc; cLproton_mc.Boost(-boost_mc);
  TLorentzVector cLneutron_mc = Lneutron_mc; cLneutron_mc.Boost(-boost_mc);

  TString frame[2] = {"Lab","CM"};

  double lam_mc_mass[2]  = { Llambda_mc.M()                           , cLlambda_mc.M()};
  double lam_mc_mom[2]   = { Llambda_mc.Vect().Mag()                  , cLlambda_mc.Vect().Mag()};
  double lam_mc_cost[2]  = { Llambda_mc.Vect().CosTheta()             , cLlambda_mc.Vect().CosTheta()};
  double lam_mc_phi[2]   = { Llambda_mc.Vect().Phi()                  , cLlambda_mc.Vect().Phi()};
  double lam_mc_mmass[2] = { (Ltgt_mc+Lbeam_mc-Llambda_mc).M()              , (cLtgt_mc+cLbeam_mc-cLlambda_mc).M()};
  double lam_mc_mmom[2]  = { (Ltgt_mc+Lbeam_mc-Llambda_mc).Vect().Mag()     , (cLtgt_mc+cLbeam_mc-cLlambda_mc).Vect().Mag()};
  double lam_mc_mcost[2] = { (Ltgt_mc+Lbeam_mc-Llambda_mc).Vect().CosTheta(), (cLtgt_mc+cLbeam_mc-cLlambda_mc).Vect().CosTheta()};

  double p_mc_mass[2]  = { Lproton_mc.M()                           , cLproton_mc.M()};
  double p_mc_mom[2]   = { Lproton_mc.Vect().Mag()                  , cLproton_mc.Vect().Mag()};
  double p_mc_cost[2]  = { Lproton_mc.Vect().CosTheta()             , cLproton_mc.Vect().CosTheta()};
  double p_mc_phi[2]   = { Lproton_mc.Vect().Phi()                  , cLproton_mc.Vect().Phi()};
  double p_mc_mmass[2] = { (Ltgt_mc+Lbeam_mc-Lproton_mc).M()              , (cLtgt_mc+cLbeam_mc-cLproton_mc).M()};
  double p_mc_mmom[2]  = { (Ltgt_mc+Lbeam_mc-Lproton_mc).Vect().Mag()     , (cLtgt_mc+cLbeam_mc-cLproton_mc).Vect().Mag()};
  double p_mc_mcost[2] = { (Ltgt_mc+Lbeam_mc-Lproton_mc).Vect().CosTheta(), (cLtgt_mc+cLbeam_mc-cLproton_mc).Vect().CosTheta()};

  double n_mc_mass[2]  = { Lneutron_mc.M()                           , cLneutron_mc.M()};
  double n_mc_mom[2]   = { Lneutron_mc.Vect().Mag()                  , cLneutron_mc.Vect().Mag()};
  double n_mc_cost[2]  = { Lneutron_mc.Vect().CosTheta()             , cLneutron_mc.Vect().CosTheta()};
  double n_mc_phi[2]   = { Lneutron_mc.Vect().Phi()                  , cLneutron_mc.Vect().Phi()};
  double n_mc_mmass[2] = { (Ltgt_mc+Lbeam_mc-Lneutron_mc).M()              , (cLtgt_mc+cLbeam_mc-cLneutron_mc).M()};
  double n_mc_mmom[2]  = { (Ltgt_mc+Lbeam_mc-Lneutron_mc).Vect().Mag()     , (cLtgt_mc+cLbeam_mc-cLneutron_mc).Vect().Mag()};
  double n_mc_mcost[2] = { (Ltgt_mc+Lbeam_mc-Lneutron_mc).Vect().CosTheta(), (cLtgt_mc+cLbeam_mc-cLneutron_mc).Vect().CosTheta()};

  TLorentzVector Lplam_mc = Llambda_mc + Lproton_mc; TLorentzVector cLplam_mc = cLlambda_mc + cLproton_mc;

  double plam_mc_mass[2]  = { Lplam_mc.M()                           , cLplam_mc.M()};
  double plam_mc_mom[2]   = { Lplam_mc.Vect().Mag()                  , cLplam_mc.Vect().Mag()};
  double plam_mc_cost[2]  = { Lplam_mc.Vect().CosTheta()             , cLplam_mc.Vect().CosTheta()};
  double plam_mc_phi[2]   = { Lplam_mc.Vect().Phi()                  , cLplam_mc.Vect().Phi()};
  double plam_mc_mmass[2] = { (Ltgt_mc+Lbeam_mc-Lplam_mc).M()              , (cLtgt_mc+cLbeam_mc-cLplam_mc).M()};
  double plam_mc_mmom[2]  = { (Ltgt_mc+Lbeam_mc-Lplam_mc).Vect().Mag()     , (cLtgt_mc+cLbeam_mc-cLplam_mc).Vect().Mag()};
  double plam_mc_mcost[2] = { (Ltgt_mc+Lbeam_mc-Lplam_mc).Vect().CosTheta(), (cLtgt_mc+cLbeam_mc-cLplam_mc).Vect().CosTheta()};

  TVector3 boost_plam_mc = (Lplam_mc).BoostVector();
  TLorentzVector cLproton_plam_mc = Lproton_mc; cLproton_plam_mc.Boost(-boost_plam_mc);
  TLorentzVector cLlambda_plam_mc = Llambda_mc; cLlambda_plam_mc.Boost(-boost_plam_mc);

  double lam_plam_mc_mass[2]  = { Llambda_mc.M()                           , cLlambda_plam_mc.M()};
  double lam_plam_mc_mom[2]   = { Llambda_mc.Vect().Mag()                  , cLlambda_plam_mc.Vect().Mag()};
  double lam_plam_mc_cost[2]  = { Llambda_mc.Vect().CosTheta()             , cLlambda_plam_mc.Vect().CosTheta()};
  double lam_plam_mc_phi[2]   = { Llambda_mc.Vect().Phi()                  , cLlambda_plam_mc.Vect().Phi()};

  double p_plam_mc_mass[2]  = { Lproton_mc.M()                           , cLproton_plam_mc.M()};
  double p_plam_mc_mom[2]   = { Lproton_mc.Vect().Mag()                  , cLproton_plam_mc.Vect().Mag()};
  double p_plam_mc_cost[2]  = { Lproton_mc.Vect().CosTheta()             , cLproton_plam_mc.Vect().CosTheta()};
  double p_plam_mc_phi[2]   = { Lproton_mc.Vect().Phi()                  , cLproton_plam_mc.Vect().Phi()};
  double q_mc[2] = {(Lbeam_mc.Vect()-Lneutron_mc.Vect()).Mag(), (Lbeam_mc.Vect()-Lneutron_mc.Vect()).Mag()};

  // lambda
  for(int f=1; f<2; f++){
    FillHist(Form("lam_Mass_%s",frame[f].Data()),lam_mc_mass[f]);
    FillHist(Form("lam_Momentum_%s",frame[f].Data()),lam_mc_mom[f]);
    FillHist(Form("lam_CosT_%s",frame[f].Data()),lam_mc_cost[f]);
    FillHist(Form("lam_Phi_%s",frame[f].Data()),lam_mc_phi[f]);
    FillHist(Form("lam_MMass_%s",frame[f].Data()),lam_mc_mmass[f]);
    FillHist(Form("lam_MMomentum_%s",frame[f].Data()),lam_mc_mmom[f]);
    FillHist(Form("lam_MCosT_%s",frame[f].Data()),lam_mc_mcost[f]);
  }
  // proton
  for(int f=1; f<2; f++){
    FillHist(Form("p_Mass_%s",frame[f].Data()),p_mc_mass[f]);
    FillHist(Form("p_Momentum_%s",frame[f].Data()),p_mc_mom[f]);
    FillHist(Form("p_CosT_%s",frame[f].Data()),p_mc_cost[f]);
    FillHist(Form("p_Phi_%s",frame[f].Data()),p_mc_phi[f]);
    FillHist(Form("p_MMass_%s",frame[f].Data()),p_mc_mmass[f]);
    FillHist(Form("p_MMomentum_%s",frame[f].Data()),p_mc_mmom[f]);
    FillHist(Form("p_MCosT_%s",frame[f].Data()),p_mc_mcost[f]);
  }
  // lambda in pLam frame
  for(int f=1; f<2; f++){
    FillHist(Form("lam_plam_Mass_%s",frame[f].Data()),lam_plam_mc_mass[f]);
    FillHist(Form("lam_plam_Momentum_%s",frame[f].Data()),lam_plam_mc_mom[f]);
    FillHist(Form("lam_plam_CosT_%s",frame[f].Data()),lam_plam_mc_cost[f]);
    FillHist(Form("lam_plam_Phi_%s",frame[f].Data()),lam_plam_mc_phi[f]);
  }
  // proton in pLam frame
  for(int f=1; f<2; f++){
    FillHist(Form("p_plam_Mass_%s",frame[f].Data()),p_plam_mc_mass[f]);
    FillHist(Form("p_plam_Momentum_%s",frame[f].Data()),p_plam_mc_mom[f]);
    FillHist(Form("p_plam_CosT_%s",frame[f].Data()),p_plam_mc_cost[f]);
    FillHist(Form("p_plam_Phi_%s",frame[f].Data()),p_plam_mc_phi[f]);
  }
  // neutron
  for(int f=1; f<2; f++){
    FillHist(Form("n_Mass_%s",frame[f].Data()),n_mc_mass[f]);
    FillHist(Form("n_Momentum_%s",frame[f].Data()),n_mc_mom[f]);
    FillHist(Form("n_CosT_%s",frame[f].Data()),n_mc_cost[f]);
    FillHist(Form("n_Phi_%s",frame[f].Data()),n_mc_phi[f]);
    FillHist(Form("n_MMass_%s",frame[f].Data()),n_mc_mmass[f]);
    FillHist(Form("n_MMomentum_%s",frame[f].Data()),n_mc_mmom[f]);
    FillHist(Form("n_MCosT_%s",frame[f].Data()),n_mc_mcost[f]);
  }
  // p + lambda
  for(int f=1; f<2; f++){
    FillHist(Form("plam_Mass_%s",frame[f].Data()),plam_mc_mass[f]);
    FillHist(Form("plam_Momentum_%s",frame[f].Data()),plam_mc_mom[f]);
    FillHist(Form("plam_CosT_%s",frame[f].Data()),plam_mc_cost[f]);
    FillHist(Form("plam_Phi_%s",frame[f].Data()),plam_mc_phi[f]);
    FillHist(Form("plam_MMass_%s",frame[f].Data()),plam_mc_mmass[f]);
    FillHist(Form("plam_MMomentum_%s",frame[f].Data()),plam_mc_mmom[f]);
    FillHist(Form("plam_MCosT_%s",frame[f].Data()),plam_mc_mcost[f]);
  }
	// IM(pipp) vs. MomentumTransfer
	for(int f=1; f<2; f++){
		FillHist(Form("pipp_IMvsMomentumTransfer_%s",frame[f].Data()),plam_mc_mass[f],q_mc[f]);
		FillHist(Form("pipp_IMvsHeKpipp_MCosT_%s",frame[f].Data()),plam_mc_mass[f],plam_mc_mcost[f]);
	}

  int reduction = 1;
  if(particle->nBeam()!=1) return false;
  // p + lambda
  for(int f=1; f<2; f++){
    FillHist(Form("plam_Mass_%s_%02d",frame[f].Data(),reduction),plam_mc_mass[f]);
    FillHist(Form("plam_Momentum_%s_%02d",frame[f].Data(),reduction),plam_mc_mom[f]);
    FillHist(Form("plam_CosT_%s_%02d",frame[f].Data(),reduction),plam_mc_cost[f]);
    FillHist(Form("plam_Phi_%s_%02d",frame[f].Data(),reduction),plam_mc_phi[f]);
    FillHist(Form("plam_MMass_%s_%02d",frame[f].Data(),reduction),plam_mc_mmass[f]);
    FillHist(Form("plam_MMomentum_%s_%02d",frame[f].Data(),reduction),plam_mc_mmom[f]);
    FillHist(Form("plam_MCosT_%s_%02d",frame[f].Data(),reduction),plam_mc_mcost[f]);
    FillHist(Form("n_Mass_%s_%02d",frame[f].Data(),reduction),n_mc_mass[f]);
    FillHist(Form("n_Momentum_%s_%02d",frame[f].Data(),reduction),n_mc_mom[f]);
    FillHist(Form("n_CosT_%s_%02d",frame[f].Data(),reduction),n_mc_cost[f]);
    FillHist(Form("n_Phi_%s_%02d",frame[f].Data(),reduction),n_mc_phi[f]);
    FillHist(Form("n_MMass_%s_%02d",frame[f].Data(),reduction),n_mc_mmass[f]);
    FillHist(Form("n_MMomentum_%s_%02d",frame[f].Data(),reduction),n_mc_mmom[f]);
    FillHist(Form("n_MCosT_%s_%02d",frame[f].Data(),reduction),n_mc_mcost[f]);
    FillHist(Form("lam_plam_Mass_%s_%02d",frame[f].Data(),reduction),lam_plam_mc_mass[f]);
    FillHist(Form("lam_plam_Momentum_%s_%02d",frame[f].Data(),reduction),lam_plam_mc_mom[f]);
    FillHist(Form("lam_plam_CosT_%s_%02d",frame[f].Data(),reduction),lam_plam_mc_cost[f]);
    FillHist(Form("lam_plam_Phi_%s_%02d",frame[f].Data(),reduction),lam_plam_mc_phi[f]);
  }
  reduction++;
  pBeam* beam = particle->beam(0);
  double beamtof = beam->bhdt0tof();

  //if(!header->IsTrig(Trig_KCDH3)){ return true; }
  //for(int f=1; f<2; f++){
  //  FillHist(Form("plam_Mass_%s_%02d",frame[f].Data(),reduction),plam_mc_mass[f]);
  //  FillHist(Form("plam_Momentum_%s_%02d",frame[f].Data(),reduction),plam_mc_mom[f]);
  //  FillHist(Form("plam_CosT_%s_%02d",frame[f].Data(),reduction),plam_mc_cost[f]);
  //  FillHist(Form("plam_Phi_%s_%02d",frame[f].Data(),reduction),plam_mc_phi[f]);
  //  FillHist(Form("plam_MMass_%s_%02d",frame[f].Data(),reduction),plam_mc_mmass[f]);
  //  FillHist(Form("plam_MMomentum_%s_%02d",frame[f].Data(),reduction),plam_mc_mmom[f]);
  //  FillHist(Form("plam_MCosT_%s_%02d",frame[f].Data(),reduction),plam_mc_mcost[f]);
  //}
  //reduction++;

  int MulCDH=0;
  for(int i=0; i<cdsMan->nCDH(); i++){
    HodoscopeLikeHit* hit = cdsMan->CDH(i);
    if(hit->CheckRange()) MulCDH++;
  }
  if(MulCDH!=3){ return true; }
  for(int f=1; f<2; f++){
    FillHist(Form("plam_Mass_%s_%02d",frame[f].Data(),reduction),plam_mc_mass[f]);
    FillHist(Form("plam_Momentum_%s_%02d",frame[f].Data(),reduction),plam_mc_mom[f]);
    FillHist(Form("plam_CosT_%s_%02d",frame[f].Data(),reduction),plam_mc_cost[f]);
    FillHist(Form("plam_Phi_%s_%02d",frame[f].Data(),reduction),plam_mc_phi[f]);
    FillHist(Form("plam_MMass_%s_%02d",frame[f].Data(),reduction),plam_mc_mmass[f]);
    FillHist(Form("plam_MMomentum_%s_%02d",frame[f].Data(),reduction),plam_mc_mmom[f]);
    FillHist(Form("plam_MCosT_%s_%02d",frame[f].Data(),reduction),plam_mc_mcost[f]);
    FillHist(Form("n_Mass_%s_%02d",frame[f].Data(),reduction),n_mc_mass[f]);
    FillHist(Form("n_Momentum_%s_%02d",frame[f].Data(),reduction),n_mc_mom[f]);
    FillHist(Form("n_CosT_%s_%02d",frame[f].Data(),reduction),n_mc_cost[f]);
    FillHist(Form("n_Phi_%s_%02d",frame[f].Data(),reduction),n_mc_phi[f]);
    FillHist(Form("n_MMass_%s_%02d",frame[f].Data(),reduction),n_mc_mmass[f]);
    FillHist(Form("n_MMomentum_%s_%02d",frame[f].Data(),reduction),n_mc_mmom[f]);
    FillHist(Form("n_MCosT_%s_%02d",frame[f].Data(),reduction),n_mc_mcost[f]);
    FillHist(Form("lam_plam_Mass_%s_%02d",frame[f].Data(),reduction),lam_plam_mc_mass[f]);
    FillHist(Form("lam_plam_Momentum_%s_%02d",frame[f].Data(),reduction),lam_plam_mc_mom[f]);
    FillHist(Form("lam_plam_CosT_%s_%02d",frame[f].Data(),reduction),lam_plam_mc_cost[f]);
    FillHist(Form("lam_plam_Phi_%s_%02d",frame[f].Data(),reduction),lam_plam_mc_phi[f]);
  }
  reduction++;

  if(cdstrackMan->nGoodTrack()!=3) return true;
  for(int f=1; f<2; f++){
    FillHist(Form("plam_Mass_%s_%02d",frame[f].Data(),reduction),plam_mc_mass[f]);
    FillHist(Form("plam_Momentum_%s_%02d",frame[f].Data(),reduction),plam_mc_mom[f]);
    FillHist(Form("plam_CosT_%s_%02d",frame[f].Data(),reduction),plam_mc_cost[f]);
    FillHist(Form("plam_Phi_%s_%02d",frame[f].Data(),reduction),plam_mc_phi[f]);
    FillHist(Form("plam_MMass_%s_%02d",frame[f].Data(),reduction),plam_mc_mmass[f]);
    FillHist(Form("plam_MMomentum_%s_%02d",frame[f].Data(),reduction),plam_mc_mmom[f]);
    FillHist(Form("plam_MCosT_%s_%02d",frame[f].Data(),reduction),plam_mc_mcost[f]);
    FillHist(Form("n_Mass_%s_%02d",frame[f].Data(),reduction),n_mc_mass[f]);
    FillHist(Form("n_Momentum_%s_%02d",frame[f].Data(),reduction),n_mc_mom[f]);
    FillHist(Form("n_CosT_%s_%02d",frame[f].Data(),reduction),n_mc_cost[f]);
    FillHist(Form("n_Phi_%s_%02d",frame[f].Data(),reduction),n_mc_phi[f]);
    FillHist(Form("n_MMass_%s_%02d",frame[f].Data(),reduction),n_mc_mmass[f]);
    FillHist(Form("n_MMomentum_%s_%02d",frame[f].Data(),reduction),n_mc_mmom[f]);
    FillHist(Form("n_MCosT_%s_%02d",frame[f].Data(),reduction),n_mc_mcost[f]);
    FillHist(Form("lam_plam_Mass_%s_%02d",frame[f].Data(),reduction),lam_plam_mc_mass[f]);
    FillHist(Form("lam_plam_Momentum_%s_%02d",frame[f].Data(),reduction),lam_plam_mc_mom[f]);
    FillHist(Form("lam_plam_CosT_%s_%02d",frame[f].Data(),reduction),lam_plam_mc_cost[f]);
    FillHist(Form("lam_plam_Phi_%s_%02d",frame[f].Data(),reduction),lam_plam_mc_phi[f]);
  }
  reduction++;

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
  if(MulBVC!=0||MulCVC!=0||MulPC!=0) return true;
  for(int f=1; f<2; f++){
    FillHist(Form("plam_Mass_%s_%02d",frame[f].Data(),reduction),plam_mc_mass[f]);
    FillHist(Form("plam_Momentum_%s_%02d",frame[f].Data(),reduction),plam_mc_mom[f]);
    FillHist(Form("plam_CosT_%s_%02d",frame[f].Data(),reduction),plam_mc_cost[f]);
    FillHist(Form("plam_Phi_%s_%02d",frame[f].Data(),reduction),plam_mc_phi[f]);
    FillHist(Form("plam_MMass_%s_%02d",frame[f].Data(),reduction),plam_mc_mmass[f]);
    FillHist(Form("plam_MMomentum_%s_%02d",frame[f].Data(),reduction),plam_mc_mmom[f]);
    FillHist(Form("plam_MCosT_%s_%02d",frame[f].Data(),reduction),plam_mc_mcost[f]);
    FillHist(Form("n_Mass_%s_%02d",frame[f].Data(),reduction),n_mc_mass[f]);
    FillHist(Form("n_Momentum_%s_%02d",frame[f].Data(),reduction),n_mc_mom[f]);
    FillHist(Form("n_CosT_%s_%02d",frame[f].Data(),reduction),n_mc_cost[f]);
    FillHist(Form("n_Phi_%s_%02d",frame[f].Data(),reduction),n_mc_phi[f]);
    FillHist(Form("n_MMass_%s_%02d",frame[f].Data(),reduction),n_mc_mmass[f]);
    FillHist(Form("n_MMomentum_%s_%02d",frame[f].Data(),reduction),n_mc_mmom[f]);
    FillHist(Form("n_MCosT_%s_%02d",frame[f].Data(),reduction),n_mc_mcost[f]);
    FillHist(Form("lam_plam_Mass_%s_%02d",frame[f].Data(),reduction),lam_plam_mc_mass[f]);
    FillHist(Form("lam_plam_Momentum_%s_%02d",frame[f].Data(),reduction),lam_plam_mc_mom[f]);
    FillHist(Form("lam_plam_CosT_%s_%02d",frame[f].Data(),reduction),lam_plam_mc_cost[f]);
    FillHist(Form("lam_plam_Phi_%s_%02d",frame[f].Data(),reduction),lam_plam_mc_phi[f]);
  }
  reduction++;

  int MulIH=0;
  for(int i=0; i<cdsMan->nIH(); i++){
    HodoscopeLikeHit* hit = cdsMan->IH(i);
    if(hit->CheckRange()) MulIH++;
  }
  //if(MulIH<=3){ return true; }

  if(particle->nCDS()!=3) return true;
  for(int f=1; f<2; f++){
    FillHist(Form("plam_Mass_%s_%02d",frame[f].Data(),reduction),plam_mc_mass[f]);
    FillHist(Form("plam_Momentum_%s_%02d",frame[f].Data(),reduction),plam_mc_mom[f]);
    FillHist(Form("plam_CosT_%s_%02d",frame[f].Data(),reduction),plam_mc_cost[f]);
    FillHist(Form("plam_Phi_%s_%02d",frame[f].Data(),reduction),plam_mc_phi[f]);
    FillHist(Form("plam_MMass_%s_%02d",frame[f].Data(),reduction),plam_mc_mmass[f]);
    FillHist(Form("plam_MMomentum_%s_%02d",frame[f].Data(),reduction),plam_mc_mmom[f]);
    FillHist(Form("plam_MCosT_%s_%02d",frame[f].Data(),reduction),plam_mc_mcost[f]);
    FillHist(Form("n_Mass_%s_%02d",frame[f].Data(),reduction),n_mc_mass[f]);
    FillHist(Form("n_Momentum_%s_%02d",frame[f].Data(),reduction),n_mc_mom[f]);
    FillHist(Form("n_CosT_%s_%02d",frame[f].Data(),reduction),n_mc_cost[f]);
    FillHist(Form("n_Phi_%s_%02d",frame[f].Data(),reduction),n_mc_phi[f]);
    FillHist(Form("n_MMass_%s_%02d",frame[f].Data(),reduction),n_mc_mmass[f]);
    FillHist(Form("n_MMomentum_%s_%02d",frame[f].Data(),reduction),n_mc_mmom[f]);
    FillHist(Form("n_MCosT_%s_%02d",frame[f].Data(),reduction),n_mc_mcost[f]);
    FillHist(Form("lam_plam_Mass_%s_%02d",frame[f].Data(),reduction),lam_plam_mc_mass[f]);
    FillHist(Form("lam_plam_Momentum_%s_%02d",frame[f].Data(),reduction),lam_plam_mc_mom[f]);
    FillHist(Form("lam_plam_CosT_%s_%02d",frame[f].Data(),reduction),lam_plam_mc_cost[f]);
    FillHist(Form("lam_plam_Phi_%s_%02d",frame[f].Data(),reduction),lam_plam_mc_phi[f]);
  }
  reduction++;

  if(particle->nProton()!=2) return true;
  for(int f=1; f<2; f++){
    FillHist(Form("plam_Mass_%s_%02d",frame[f].Data(),reduction),plam_mc_mass[f]);
    FillHist(Form("plam_Momentum_%s_%02d",frame[f].Data(),reduction),plam_mc_mom[f]);
    FillHist(Form("plam_CosT_%s_%02d",frame[f].Data(),reduction),plam_mc_cost[f]);
    FillHist(Form("plam_Phi_%s_%02d",frame[f].Data(),reduction),plam_mc_phi[f]);
    FillHist(Form("plam_MMass_%s_%02d",frame[f].Data(),reduction),plam_mc_mmass[f]);
    FillHist(Form("plam_MMomentum_%s_%02d",frame[f].Data(),reduction),plam_mc_mmom[f]);
    FillHist(Form("plam_MCosT_%s_%02d",frame[f].Data(),reduction),plam_mc_mcost[f]);
    FillHist(Form("n_Mass_%s_%02d",frame[f].Data(),reduction),n_mc_mass[f]);
    FillHist(Form("n_Momentum_%s_%02d",frame[f].Data(),reduction),n_mc_mom[f]);
    FillHist(Form("n_CosT_%s_%02d",frame[f].Data(),reduction),n_mc_cost[f]);
    FillHist(Form("n_Phi_%s_%02d",frame[f].Data(),reduction),n_mc_phi[f]);
    FillHist(Form("n_MMass_%s_%02d",frame[f].Data(),reduction),n_mc_mmass[f]);
    FillHist(Form("n_MMomentum_%s_%02d",frame[f].Data(),reduction),n_mc_mmom[f]);
    FillHist(Form("n_MCosT_%s_%02d",frame[f].Data(),reduction),n_mc_mcost[f]);
    FillHist(Form("lam_plam_Mass_%s_%02d",frame[f].Data(),reduction),lam_plam_mc_mass[f]);
    FillHist(Form("lam_plam_Momentum_%s_%02d",frame[f].Data(),reduction),lam_plam_mc_mom[f]);
    FillHist(Form("lam_plam_CosT_%s_%02d",frame[f].Data(),reduction),lam_plam_mc_cost[f]);
    FillHist(Form("lam_plam_Phi_%s_%02d",frame[f].Data(),reduction),lam_plam_mc_phi[f]);
  }
  reduction++;
  if(particle->nPiminus()!=1) return true;
  for(int f=1; f<2; f++){
    FillHist(Form("plam_Mass_%s_%02d",frame[f].Data(),reduction),plam_mc_mass[f]);
    FillHist(Form("plam_Momentum_%s_%02d",frame[f].Data(),reduction),plam_mc_mom[f]);
    FillHist(Form("plam_CosT_%s_%02d",frame[f].Data(),reduction),plam_mc_cost[f]);
    FillHist(Form("plam_Phi_%s_%02d",frame[f].Data(),reduction),plam_mc_phi[f]);
    FillHist(Form("plam_MMass_%s_%02d",frame[f].Data(),reduction),plam_mc_mmass[f]);
    FillHist(Form("plam_MMomentum_%s_%02d",frame[f].Data(),reduction),plam_mc_mmom[f]);
    FillHist(Form("plam_MCosT_%s_%02d",frame[f].Data(),reduction),plam_mc_mcost[f]);
    FillHist(Form("n_Mass_%s_%02d",frame[f].Data(),reduction),n_mc_mass[f]);
    FillHist(Form("n_Momentum_%s_%02d",frame[f].Data(),reduction),n_mc_mom[f]);
    FillHist(Form("n_CosT_%s_%02d",frame[f].Data(),reduction),n_mc_cost[f]);
    FillHist(Form("n_Phi_%s_%02d",frame[f].Data(),reduction),n_mc_phi[f]);
    FillHist(Form("n_MMass_%s_%02d",frame[f].Data(),reduction),n_mc_mmass[f]);
    FillHist(Form("n_MMomentum_%s_%02d",frame[f].Data(),reduction),n_mc_mmom[f]);
    FillHist(Form("n_MCosT_%s_%02d",frame[f].Data(),reduction),n_mc_mcost[f]);
    FillHist(Form("lam_plam_Mass_%s_%02d",frame[f].Data(),reduction),lam_plam_mc_mass[f]);
    FillHist(Form("lam_plam_Momentum_%s_%02d",frame[f].Data(),reduction),lam_plam_mc_mom[f]);
    FillHist(Form("lam_plam_CosT_%s_%02d",frame[f].Data(),reduction),lam_plam_mc_cost[f]);
    FillHist(Form("lam_plam_Phi_%s_%02d",frame[f].Data(),reduction),lam_plam_mc_phi[f]);
  }
  reduction++;

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
  for(int f=1; f<2; f++){
    FillHist(Form("plam_Mass_%s_%02d",frame[f].Data(),reduction),plam_mc_mass[f]);
    FillHist(Form("plam_Momentum_%s_%02d",frame[f].Data(),reduction),plam_mc_mom[f]);
    FillHist(Form("plam_CosT_%s_%02d",frame[f].Data(),reduction),plam_mc_cost[f]);
    FillHist(Form("plam_Phi_%s_%02d",frame[f].Data(),reduction),plam_mc_phi[f]);
    FillHist(Form("plam_MMass_%s_%02d",frame[f].Data(),reduction),plam_mc_mmass[f]);
    FillHist(Form("plam_MMomentum_%s_%02d",frame[f].Data(),reduction),plam_mc_mmom[f]);
    FillHist(Form("plam_MCosT_%s_%02d",frame[f].Data(),reduction),plam_mc_mcost[f]);
    FillHist(Form("n_Mass_%s_%02d",frame[f].Data(),reduction),n_mc_mass[f]);
    FillHist(Form("n_Momentum_%s_%02d",frame[f].Data(),reduction),n_mc_mom[f]);
    FillHist(Form("n_CosT_%s_%02d",frame[f].Data(),reduction),n_mc_cost[f]);
    FillHist(Form("n_Phi_%s_%02d",frame[f].Data(),reduction),n_mc_phi[f]);
    FillHist(Form("n_MMass_%s_%02d",frame[f].Data(),reduction),n_mc_mmass[f]);
    FillHist(Form("n_MMomentum_%s_%02d",frame[f].Data(),reduction),n_mc_mmom[f]);
    FillHist(Form("n_MCosT_%s_%02d",frame[f].Data(),reduction),n_mc_mcost[f]);
    FillHist(Form("lam_plam_Mass_%s_%02d",frame[f].Data(),reduction),lam_plam_mc_mass[f]);
    FillHist(Form("lam_plam_Momentum_%s_%02d",frame[f].Data(),reduction),lam_plam_mc_mom[f]);
    FillHist(Form("lam_plam_CosT_%s_%02d",frame[f].Data(),reduction),lam_plam_mc_cost[f]);
    FillHist(Form("lam_plam_Phi_%s_%02d",frame[f].Data(),reduction),lam_plam_mc_phi[f]);
  }
  reduction++;

  pCDS* pim = particle->pim(0);
  //pCDS* p1 = particle->proton(fcds);
  //pCDS* p2 = particle->proton(1-fcds);
  // #################### //
  // p1/p2 decision for MC//
  // #################### //
  int lamnum = -1;
  if(particle->proton(0)->cdhseg(0)==lamseg){
    lamnum = 0;
  }
  else{
    lamnum = 1;
  }
  if(lamnum<0) return true;
  pCDS* p1 = particle->proton(lamnum);
  pCDS* p2 = particle->proton(1-lamnum);

  // ############## //
  // CDC Chi-square //
  // ############## //
  if(pim->chi()>30) return true;
  for(int f=1; f<2; f++){
    FillHist(Form("plam_Mass_%s_%02d",frame[f].Data(),reduction),plam_mc_mass[f]);
    FillHist(Form("plam_Momentum_%s_%02d",frame[f].Data(),reduction),plam_mc_mom[f]);
    FillHist(Form("plam_CosT_%s_%02d",frame[f].Data(),reduction),plam_mc_cost[f]);
    FillHist(Form("plam_Phi_%s_%02d",frame[f].Data(),reduction),plam_mc_phi[f]);
    FillHist(Form("plam_MMass_%s_%02d",frame[f].Data(),reduction),plam_mc_mmass[f]);
    FillHist(Form("plam_MMomentum_%s_%02d",frame[f].Data(),reduction),plam_mc_mmom[f]);
    FillHist(Form("plam_MCosT_%s_%02d",frame[f].Data(),reduction),plam_mc_mcost[f]);
    FillHist(Form("n_Mass_%s_%02d",frame[f].Data(),reduction),n_mc_mass[f]);
    FillHist(Form("n_Momentum_%s_%02d",frame[f].Data(),reduction),n_mc_mom[f]);
    FillHist(Form("n_CosT_%s_%02d",frame[f].Data(),reduction),n_mc_cost[f]);
    FillHist(Form("n_Phi_%s_%02d",frame[f].Data(),reduction),n_mc_phi[f]);
    FillHist(Form("n_MMass_%s_%02d",frame[f].Data(),reduction),n_mc_mmass[f]);
    FillHist(Form("n_MMomentum_%s_%02d",frame[f].Data(),reduction),n_mc_mmom[f]);
    FillHist(Form("n_MCosT_%s_%02d",frame[f].Data(),reduction),n_mc_mcost[f]);
    FillHist(Form("lam_plam_Mass_%s_%02d",frame[f].Data(),reduction),lam_plam_mc_mass[f]);
    FillHist(Form("lam_plam_Momentum_%s_%02d",frame[f].Data(),reduction),lam_plam_mc_mom[f]);
    FillHist(Form("lam_plam_CosT_%s_%02d",frame[f].Data(),reduction),lam_plam_mc_cost[f]);
    FillHist(Form("lam_plam_Phi_%s_%02d",frame[f].Data(),reduction),lam_plam_mc_phi[f]);
  }
  reduction++;
  if(p1->chi()>30) return true;
  for(int f=1; f<2; f++){
    FillHist(Form("plam_Mass_%s_%02d",frame[f].Data(),reduction),plam_mc_mass[f]);
    FillHist(Form("plam_Momentum_%s_%02d",frame[f].Data(),reduction),plam_mc_mom[f]);
    FillHist(Form("plam_CosT_%s_%02d",frame[f].Data(),reduction),plam_mc_cost[f]);
    FillHist(Form("plam_Phi_%s_%02d",frame[f].Data(),reduction),plam_mc_phi[f]);
    FillHist(Form("plam_MMass_%s_%02d",frame[f].Data(),reduction),plam_mc_mmass[f]);
    FillHist(Form("plam_MMomentum_%s_%02d",frame[f].Data(),reduction),plam_mc_mmom[f]);
    FillHist(Form("plam_MCosT_%s_%02d",frame[f].Data(),reduction),plam_mc_mcost[f]);
    FillHist(Form("n_Mass_%s_%02d",frame[f].Data(),reduction),n_mc_mass[f]);
    FillHist(Form("n_Momentum_%s_%02d",frame[f].Data(),reduction),n_mc_mom[f]);
    FillHist(Form("n_CosT_%s_%02d",frame[f].Data(),reduction),n_mc_cost[f]);
    FillHist(Form("n_Phi_%s_%02d",frame[f].Data(),reduction),n_mc_phi[f]);
    FillHist(Form("n_MMass_%s_%02d",frame[f].Data(),reduction),n_mc_mmass[f]);
    FillHist(Form("n_MMomentum_%s_%02d",frame[f].Data(),reduction),n_mc_mmom[f]);
    FillHist(Form("n_MCosT_%s_%02d",frame[f].Data(),reduction),n_mc_mcost[f]);
    FillHist(Form("lam_plam_Mass_%s_%02d",frame[f].Data(),reduction),lam_plam_mc_mass[f]);
    FillHist(Form("lam_plam_Momentum_%s_%02d",frame[f].Data(),reduction),lam_plam_mc_mom[f]);
    FillHist(Form("lam_plam_CosT_%s_%02d",frame[f].Data(),reduction),lam_plam_mc_cost[f]);
    FillHist(Form("lam_plam_Phi_%s_%02d",frame[f].Data(),reduction),lam_plam_mc_phi[f]);
  }
  reduction++;
  if(p2->chi()>30) return true;
  for(int f=1; f<2; f++){
    FillHist(Form("plam_Mass_%s_%02d",frame[f].Data(),reduction),plam_mc_mass[f]);
    FillHist(Form("plam_Momentum_%s_%02d",frame[f].Data(),reduction),plam_mc_mom[f]);
    FillHist(Form("plam_CosT_%s_%02d",frame[f].Data(),reduction),plam_mc_cost[f]);
    FillHist(Form("plam_Phi_%s_%02d",frame[f].Data(),reduction),plam_mc_phi[f]);
    FillHist(Form("plam_MMass_%s_%02d",frame[f].Data(),reduction),plam_mc_mmass[f]);
    FillHist(Form("plam_MMomentum_%s_%02d",frame[f].Data(),reduction),plam_mc_mmom[f]);
    FillHist(Form("plam_MCosT_%s_%02d",frame[f].Data(),reduction),plam_mc_mcost[f]);
    FillHist(Form("n_Mass_%s_%02d",frame[f].Data(),reduction),n_mc_mass[f]);
    FillHist(Form("n_Momentum_%s_%02d",frame[f].Data(),reduction),n_mc_mom[f]);
    FillHist(Form("n_CosT_%s_%02d",frame[f].Data(),reduction),n_mc_cost[f]);
    FillHist(Form("n_Phi_%s_%02d",frame[f].Data(),reduction),n_mc_phi[f]);
    FillHist(Form("n_MMass_%s_%02d",frame[f].Data(),reduction),n_mc_mmass[f]);
    FillHist(Form("n_MMomentum_%s_%02d",frame[f].Data(),reduction),n_mc_mmom[f]);
    FillHist(Form("n_MCosT_%s_%02d",frame[f].Data(),reduction),n_mc_mcost[f]);
    FillHist(Form("lam_plam_Mass_%s_%02d",frame[f].Data(),reduction),lam_plam_mc_mass[f]);
    FillHist(Form("lam_plam_Momentum_%s_%02d",frame[f].Data(),reduction),lam_plam_mc_mom[f]);
    FillHist(Form("lam_plam_CosT_%s_%02d",frame[f].Data(),reduction),lam_plam_mc_cost[f]);
    FillHist(Form("lam_plam_Phi_%s_%02d",frame[f].Data(),reduction),lam_plam_mc_phi[f]);
  }
  reduction++;

  // ######################### //
  // 2 Helix calculation       //
  // ######################### //
  pCDS* pip1 = 0;
  pCDS* pip2 = 0;
  {
    for(int it=0; it<particle->nProduct(); it++){
      pCDS* product = particle->product(it);
      int comb = product->comb();
      if(comb==TMath::Power(2,CDS_Proton)+TMath::Power(2,CDS_PiMinus)){
        if(product->daughter1()==p1->id()||product->daughter2()==p1->id()){
          pip1 = product;
        }
        else if(product->daughter1()==p2->id()||product->daughter2()==p2->id()){
          pip2 = product;
        }
      }
    }
  }
  if(pip1==0||pip2==0){ return true; }
  for(int f=1; f<2; f++){
    FillHist(Form("plam_Mass_%s_%02d",frame[f].Data(),reduction),plam_mc_mass[f]);
    FillHist(Form("plam_Momentum_%s_%02d",frame[f].Data(),reduction),plam_mc_mom[f]);
    FillHist(Form("plam_CosT_%s_%02d",frame[f].Data(),reduction),plam_mc_cost[f]);
    FillHist(Form("plam_Phi_%s_%02d",frame[f].Data(),reduction),plam_mc_phi[f]);
    FillHist(Form("plam_MMass_%s_%02d",frame[f].Data(),reduction),plam_mc_mmass[f]);
    FillHist(Form("plam_MMomentum_%s_%02d",frame[f].Data(),reduction),plam_mc_mmom[f]);
    FillHist(Form("plam_MCosT_%s_%02d",frame[f].Data(),reduction),plam_mc_mcost[f]);
    FillHist(Form("n_Mass_%s_%02d",frame[f].Data(),reduction),n_mc_mass[f]);
    FillHist(Form("n_Momentum_%s_%02d",frame[f].Data(),reduction),n_mc_mom[f]);
    FillHist(Form("n_CosT_%s_%02d",frame[f].Data(),reduction),n_mc_cost[f]);
    FillHist(Form("n_Phi_%s_%02d",frame[f].Data(),reduction),n_mc_phi[f]);
    FillHist(Form("n_MMass_%s_%02d",frame[f].Data(),reduction),n_mc_mmass[f]);
    FillHist(Form("n_MMomentum_%s_%02d",frame[f].Data(),reduction),n_mc_mmom[f]);
    FillHist(Form("n_MCosT_%s_%02d",frame[f].Data(),reduction),n_mc_mcost[f]);
    FillHist(Form("lam_plam_Mass_%s_%02d",frame[f].Data(),reduction),lam_plam_mc_mass[f]);
    FillHist(Form("lam_plam_Momentum_%s_%02d",frame[f].Data(),reduction),lam_plam_mc_mom[f]);
    FillHist(Form("lam_plam_CosT_%s_%02d",frame[f].Data(),reduction),lam_plam_mc_cost[f]);
    FillHist(Form("lam_plam_Phi_%s_%02d",frame[f].Data(),reduction),lam_plam_mc_phi[f]);
  }
  reduction++;


  // Target
  TVector3 tmpZeroV;
  TLorentzVector tmpLtgt; tmpLtgt.SetVectM(tmpZeroV,ThreeHeMass);
  TLorentzVector tmpLbeam = beam->GetLorentzVector(vertex);

  TLorentzVector tmpLpim = pim->GetLorentzVector();
  TLorentzVector tmpLp1 = p1->GetLorentzVector();
  TLorentzVector tmpLp2 = p2->GetLorentzVector();
  double tmp_im = (tmpLpim+tmpLp1+tmpLp2).M();
  double tmp_mm = (tmpLtgt+tmpLbeam-tmpLpim-tmpLp1-tmpLp2).M();

  bool pip1flag = false, pip2flag = false;
  pCDS* lam = 0;
  pCDS* nlam = 0;
  pCDS* proton = 0;
  pCDS* nproton = 0;
  /* For pip1 */
	TVector3 v_pip1 = pip1->vbeam();
	TVector3 v_p2 = p2->vbeam();
	double pip1p2dis = (v_pip1-v_p2).Mag();
  double tmpper1[5] = {pip1->mass(),pip1->vdis(),pip1->pbdca(),p2->vbdis(),pip1p2dis};
  double tmppdf1[5] = {0,0,0,0,0};
  /* Input : mass, dca_pip, dca_Lk, dca_pk, dcaLp */
  /*Output : mass, dca_pip, dca_Lk, dca_pk, dcaLp */
  TrackTools::PDFLambda(tmpper1,tmppdf1);
  //tmppdf1[4] = 1.0;
	/* For pip2 */
	TVector3 v_pip2 = pip2->vbeam();
	TVector3 v_p1 = p1->vbeam();
	double pip2p1dis = (v_pip2-v_p1).Mag();
  double tmpper2[5] = {pip2->mass(),pip2->vdis(),pip2->pbdca(),p1->vbdis(),pip2p1dis};
  double tmppdf2[5] = {0,0,0,0,0};
  /* Input : mass, dca_pip, dca_Lk, dca_pk, dcaLp */
  /*Output : mass, dca_pip, dca_Lk, dca_pk, dcaLp */
  TrackTools::PDFLambda(tmpper2,tmppdf2);
  //tmppdf2[4] = 1.0;
  /* All */
  double pdf1 = -TMath::Log(tmppdf1[0]*tmppdf1[1]*tmppdf1[2]*tmppdf1[3]*tmppdf1[4]);
  double pdf2 = -TMath::Log(tmppdf2[0]*tmppdf2[1]*tmppdf2[2]*tmppdf2[3]*tmppdf2[4]);
  //double pdf1 = -TMath::Log(tmppdf1[0]*tmppdf1[1]*tmppdf1[2]*tmppdf1[3]);
  //double pdf2 = -TMath::Log(tmppdf2[0]*tmppdf2[1]*tmppdf2[2]*tmppdf2[3]);
  /* Except DCA_pip */
  //double pdf1 = -TMath::Log(tmppdf1[0]*tmppdf1[2]*tmppdf1[3]);
  //double pdf2 = -TMath::Log(tmppdf2[0]*tmppdf2[2]*tmppdf2[3]);
  ///* Except DCA_Lk */
  //double pdf1 = -TMath::Log(tmppdf1[0]*tmppdf1[1]*tmppdf1[3]);
  //double pdf2 = -TMath::Log(tmppdf2[0]*tmppdf2[1]*tmppdf2[3]);
  /* Except DCA_pk */
  //double pdf1 = -TMath::Log(tmppdf1[0]*tmppdf1[1]*tmppdf1[2]);
  //double pdf2 = -TMath::Log(tmppdf2[0]*tmppdf2[1]*tmppdf2[2]);
  /* Only Mass */
  //double pdf1 = -TMath::Log(tmppdf1[0]);
  //double pdf2 = -TMath::Log(tmppdf2[0]);
  /* Only DCA_pip */
  //double pdf1 = -TMath::Log(tmppdf1[1]);
  //double pdf2 = -TMath::Log(tmppdf2[1]);
  /* Mass*DCA_pip */
  //double pdf1 = -TMath::Log(tmppdf1[0]*tmppdf1[1]);
  //double pdf2 = -TMath::Log(tmppdf2[0]*tmppdf2[1]);
  /* ALL DCA */
  //double pdf1 = -TMath::Log(tmppdf1[1]*tmppdf1[2]*tmppdf1[3]);
  //double pdf2 = -TMath::Log(tmppdf2[1]*tmppdf2[2]*tmppdf2[3]);

  double pdf=-1;
  double npdf=-1;
  //std::cout << "pdf1 : " << pdf1 <<std::endl;
  //std::cout << "pdf2 : " << pdf2 <<std::endl;
  FillHist("pip1_PDF",pdf1);
  FillHist("pip2_PDF",pdf2);
  FillHist("pip1_PDFvspip2_PDF",pdf1,pdf2);

  double pdfll =  0.0;
  //double pdful = 12.0;
  double pdful = 15.75;
  //double pdfll = 12.0;
  //double pdful = 40.0;
  //double pdfll =  0.0;
  //double pdful = 40.0;
  //double pdfll = 40.0;
  //double pdful = 999999999;

  if(pdf1<pdf2){
    FillHist("pip1_PDF_SemiSelected",pdf1);
    FillHist("pip1_PDFvspip2_SemiSelected",pdf1,pdf2);
    FillHist("lam_PDF",pdf1);
    FillHist("nlam_PDF",pdf2);
    if(pdfll<pdf1&&pdf1<pdful){
    FillHist("pip1_PDF_Selected",pdf1);
    pip1flag = true;
    lam = pip1;
    nlam = pip2;
    proton = p2;
    nproton = p1;
    pdf=pdf1;
    npdf=pdf2;
    }
  }
  else{
    FillHist("pip2_PDF_SemiSelected",pdf2);
    FillHist("pip1_PDFvspip2_SemiSelected",pdf1,pdf2);
    FillHist("lam_PDF",pdf2);
    FillHist("nlam_PDF",pdf1);
    if(pdfll<pdf2&&pdf2<pdful){
    FillHist("pip2_PDF_Selected",pdf2);
    pip2flag = true;
    pip1flag = false;
    lam = pip2;
    nlam = pip1;
    proton = p1;
    nproton = p2;
    pdf=pdf2;
    npdf=pdf1;
    }
  }
  //if(lamll<pip1->mass()&&pip1->mass()<lamul){
  //  pip1flag = true;
  //  lam = pip1;
  //  proton = p2;
  //}
  //if(lamll<pip2->mass()&&pip2->mass()<lamul){
  //  pip2flag = true;
  //  pip1flag = false;
  //  lam = pip2;
  //  proton = p1;
  //}
  if(!pip1flag&&!pip2flag){ return true; }
  for(int f=1; f<2; f++){
    FillHist(Form("plam_Mass_%s_%02d",frame[f].Data(),reduction),plam_mc_mass[f]);
    FillHist(Form("plam_Momentum_%s_%02d",frame[f].Data(),reduction),plam_mc_mom[f]);
    FillHist(Form("plam_CosT_%s_%02d",frame[f].Data(),reduction),plam_mc_cost[f]);
    FillHist(Form("plam_Phi_%s_%02d",frame[f].Data(),reduction),plam_mc_phi[f]);
    FillHist(Form("plam_MMass_%s_%02d",frame[f].Data(),reduction),plam_mc_mmass[f]);
    FillHist(Form("plam_MMomentum_%s_%02d",frame[f].Data(),reduction),plam_mc_mmom[f]);
    FillHist(Form("plam_MCosT_%s_%02d",frame[f].Data(),reduction),plam_mc_mcost[f]);
    FillHist(Form("n_Mass_%s_%02d",frame[f].Data(),reduction),n_mc_mass[f]);
    FillHist(Form("n_Momentum_%s_%02d",frame[f].Data(),reduction),n_mc_mom[f]);
    FillHist(Form("n_CosT_%s_%02d",frame[f].Data(),reduction),n_mc_cost[f]);
    FillHist(Form("n_Phi_%s_%02d",frame[f].Data(),reduction),n_mc_phi[f]);
    FillHist(Form("n_MMass_%s_%02d",frame[f].Data(),reduction),n_mc_mmass[f]);
    FillHist(Form("n_MMomentum_%s_%02d",frame[f].Data(),reduction),n_mc_mmom[f]);
    FillHist(Form("n_MCosT_%s_%02d",frame[f].Data(),reduction),n_mc_mcost[f]);
    FillHist(Form("lam_plam_Mass_%s_%02d",frame[f].Data(),reduction),lam_plam_mc_mass[f]);
    FillHist(Form("lam_plam_Momentum_%s_%02d",frame[f].Data(),reduction),lam_plam_mc_mom[f]);
    FillHist(Form("lam_plam_CosT_%s_%02d",frame[f].Data(),reduction),lam_plam_mc_cost[f]);
    FillHist(Form("lam_plam_Phi_%s_%02d",frame[f].Data(),reduction),lam_plam_mc_phi[f]);
  }
  reduction++;
  if(lam==0){ return true; }
  for(int f=1; f<2; f++){
    FillHist(Form("plam_Mass_%s_%02d",frame[f].Data(),reduction),plam_mc_mass[f]);
    FillHist(Form("plam_Momentum_%s_%02d",frame[f].Data(),reduction),plam_mc_mom[f]);
    FillHist(Form("plam_CosT_%s_%02d",frame[f].Data(),reduction),plam_mc_cost[f]);
    FillHist(Form("plam_Phi_%s_%02d",frame[f].Data(),reduction),plam_mc_phi[f]);
    FillHist(Form("plam_MMass_%s_%02d",frame[f].Data(),reduction),plam_mc_mmass[f]);
    FillHist(Form("plam_MMomentum_%s_%02d",frame[f].Data(),reduction),plam_mc_mmom[f]);
    FillHist(Form("plam_MCosT_%s_%02d",frame[f].Data(),reduction),plam_mc_mcost[f]);
    FillHist(Form("n_Mass_%s_%02d",frame[f].Data(),reduction),n_mc_mass[f]);
    FillHist(Form("n_Momentum_%s_%02d",frame[f].Data(),reduction),n_mc_mom[f]);
    FillHist(Form("n_CosT_%s_%02d",frame[f].Data(),reduction),n_mc_cost[f]);
    FillHist(Form("n_Phi_%s_%02d",frame[f].Data(),reduction),n_mc_phi[f]);
    FillHist(Form("n_MMass_%s_%02d",frame[f].Data(),reduction),n_mc_mmass[f]);
    FillHist(Form("n_MMomentum_%s_%02d",frame[f].Data(),reduction),n_mc_mmom[f]);
    FillHist(Form("n_MCosT_%s_%02d",frame[f].Data(),reduction),n_mc_mcost[f]);
    FillHist(Form("lam_plam_Mass_%s_%02d",frame[f].Data(),reduction),lam_plam_mc_mass[f]);
    FillHist(Form("lam_plam_Momentum_%s_%02d",frame[f].Data(),reduction),lam_plam_mc_mom[f]);
    FillHist(Form("lam_plam_CosT_%s_%02d",frame[f].Data(),reduction),lam_plam_mc_cost[f]);
    FillHist(Form("lam_plam_Phi_%s_%02d",frame[f].Data(),reduction),lam_plam_mc_phi[f]);
  }
  reduction++;
  FillHist("pip1_PDFvspip2_PDF_Selected",pdf1,pdf2);
  FillHist("lam_PDF_Selected",pdf);
  FillHist("nlam_PDF_Selected",npdf);

  // ################ //
  // Vertex decision1 //
  // ################ //
  const int fiducial = CID_Fiducial;
  //const int fiducial = CID_CellTube;
  //const int fiducial = CID_CellAlBe;
  //const int fiducial = CID_DEF;
  if(pip1flag){
    vertex = p2->vbeam();
    vdis = p2->vbdis();
    if(GeomTools::GetID(vertex)!=fiducial){ return true; }
    for(int f=1; f<2; f++){
      FillHist(Form("plam_Mass_%s_%02d",frame[f].Data(),reduction),plam_mc_mass[f]);
      FillHist(Form("plam_Momentum_%s_%02d",frame[f].Data(),reduction),plam_mc_mom[f]);
      FillHist(Form("plam_CosT_%s_%02d",frame[f].Data(),reduction),plam_mc_cost[f]);
      FillHist(Form("plam_Phi_%s_%02d",frame[f].Data(),reduction),plam_mc_phi[f]);
      FillHist(Form("plam_MMass_%s_%02d",frame[f].Data(),reduction),plam_mc_mmass[f]);
      FillHist(Form("plam_MMomentum_%s_%02d",frame[f].Data(),reduction),plam_mc_mmom[f]);
      FillHist(Form("plam_MCosT_%s_%02d",frame[f].Data(),reduction),plam_mc_mcost[f]);
      FillHist(Form("n_Mass_%s_%02d",frame[f].Data(),reduction),n_mc_mass[f]);
      FillHist(Form("n_Momentum_%s_%02d",frame[f].Data(),reduction),n_mc_mom[f]);
      FillHist(Form("n_CosT_%s_%02d",frame[f].Data(),reduction),n_mc_cost[f]);
      FillHist(Form("n_Phi_%s_%02d",frame[f].Data(),reduction),n_mc_phi[f]);
      FillHist(Form("n_MMass_%s_%02d",frame[f].Data(),reduction),n_mc_mmass[f]);
      FillHist(Form("n_MMomentum_%s_%02d",frame[f].Data(),reduction),n_mc_mmom[f]);
      FillHist(Form("n_MCosT_%s_%02d",frame[f].Data(),reduction),n_mc_mcost[f]);
      FillHist(Form("lam_plam_Mass_%s_%02d",frame[f].Data(),reduction),lam_plam_mc_mass[f]);
      FillHist(Form("lam_plam_Momentum_%s_%02d",frame[f].Data(),reduction),lam_plam_mc_mom[f]);
      FillHist(Form("lam_plam_CosT_%s_%02d",frame[f].Data(),reduction),lam_plam_mc_cost[f]);
      FillHist(Form("lam_plam_Phi_%s_%02d",frame[f].Data(),reduction),lam_plam_mc_phi[f]);
    }
    reduction++;
  }
  if(pip2flag){
    vertex = p1->vbeam();
    vdis = p1->vbdis();
    if(GeomTools::GetID(vertex)!=fiducial){ return true; }
    for(int f=1; f<2; f++){
      FillHist(Form("plam_Mass_%s_%02d",frame[f].Data(),reduction),plam_mc_mass[f]);
      FillHist(Form("plam_Momentum_%s_%02d",frame[f].Data(),reduction),plam_mc_mom[f]);
      FillHist(Form("plam_CosT_%s_%02d",frame[f].Data(),reduction),plam_mc_cost[f]);
      FillHist(Form("plam_Phi_%s_%02d",frame[f].Data(),reduction),plam_mc_phi[f]);
      FillHist(Form("plam_MMass_%s_%02d",frame[f].Data(),reduction),plam_mc_mmass[f]);
      FillHist(Form("plam_MMomentum_%s_%02d",frame[f].Data(),reduction),plam_mc_mmom[f]);
      FillHist(Form("plam_MCosT_%s_%02d",frame[f].Data(),reduction),plam_mc_mcost[f]);
      FillHist(Form("n_Mass_%s_%02d",frame[f].Data(),reduction),n_mc_mass[f]);
      FillHist(Form("n_Momentum_%s_%02d",frame[f].Data(),reduction),n_mc_mom[f]);
      FillHist(Form("n_CosT_%s_%02d",frame[f].Data(),reduction),n_mc_cost[f]);
      FillHist(Form("n_Phi_%s_%02d",frame[f].Data(),reduction),n_mc_phi[f]);
      FillHist(Form("n_MMass_%s_%02d",frame[f].Data(),reduction),n_mc_mmass[f]);
      FillHist(Form("n_MMomentum_%s_%02d",frame[f].Data(),reduction),n_mc_mmom[f]);
      FillHist(Form("n_MCosT_%s_%02d",frame[f].Data(),reduction),n_mc_mcost[f]);
      FillHist(Form("lam_plam_Mass_%s_%02d",frame[f].Data(),reduction),lam_plam_mc_mass[f]);
      FillHist(Form("lam_plam_Momentum_%s_%02d",frame[f].Data(),reduction),lam_plam_mc_mom[f]);
      FillHist(Form("lam_plam_CosT_%s_%02d",frame[f].Data(),reduction),lam_plam_mc_cost[f]);
      FillHist(Form("lam_plam_Phi_%s_%02d",frame[f].Data(),reduction),lam_plam_mc_phi[f]);
    }
    reduction++;
  }

  // ################ //
  // Vertex decision2 //
  // ################ //
  TVector3 vertex_lam = lam->vbeam();
  if(GeomTools::GetID(vertex_lam)!=fiducial){ return true; }
  for(int f=1; f<2; f++){
    FillHist(Form("plam_Mass_%s_%02d",frame[f].Data(),reduction),plam_mc_mass[f]);
    FillHist(Form("plam_Momentum_%s_%02d",frame[f].Data(),reduction),plam_mc_mom[f]);
    FillHist(Form("plam_CosT_%s_%02d",frame[f].Data(),reduction),plam_mc_cost[f]);
    FillHist(Form("plam_Phi_%s_%02d",frame[f].Data(),reduction),plam_mc_phi[f]);
    FillHist(Form("plam_MMass_%s_%02d",frame[f].Data(),reduction),plam_mc_mmass[f]);
    FillHist(Form("plam_MMomentum_%s_%02d",frame[f].Data(),reduction),plam_mc_mmom[f]);
    FillHist(Form("plam_MCosT_%s_%02d",frame[f].Data(),reduction),plam_mc_mcost[f]);
    FillHist(Form("n_Mass_%s_%02d",frame[f].Data(),reduction),n_mc_mass[f]);
    FillHist(Form("n_Momentum_%s_%02d",frame[f].Data(),reduction),n_mc_mom[f]);
    FillHist(Form("n_CosT_%s_%02d",frame[f].Data(),reduction),n_mc_cost[f]);
    FillHist(Form("n_Phi_%s_%02d",frame[f].Data(),reduction),n_mc_phi[f]);
    FillHist(Form("n_MMass_%s_%02d",frame[f].Data(),reduction),n_mc_mmass[f]);
    FillHist(Form("n_MMomentum_%s_%02d",frame[f].Data(),reduction),n_mc_mmom[f]);
    FillHist(Form("n_MCosT_%s_%02d",frame[f].Data(),reduction),n_mc_mcost[f]);
    FillHist(Form("lam_plam_Mass_%s_%02d",frame[f].Data(),reduction),lam_plam_mc_mass[f]);
    FillHist(Form("lam_plam_Momentum_%s_%02d",frame[f].Data(),reduction),lam_plam_mc_mom[f]);
    FillHist(Form("lam_plam_CosT_%s_%02d",frame[f].Data(),reduction),lam_plam_mc_cost[f]);
    FillHist(Form("lam_plam_Phi_%s_%02d",frame[f].Data(),reduction),lam_plam_mc_phi[f]);
  }
  reduction++;

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
  /* Proton1 */
  vtxb = p1->vbeam();
  FillHist("p1_Vertex_XY",vtxb.X(),vtxb.Y());
  FillHist("p1_Vertex_ZX",vtxb.Z(),vtxb.X());
  FillHist("p1_Vertex_ZY",vtxb.Z(),vtxb.Y());
  FillHist("p1_VBDIS",p1->vbdis());
  /* Proton2 */
  vtxb = p2->vbeam();
  FillHist("p2_Vertex_XY",vtxb.X(),vtxb.Y());
  FillHist("p2_Vertex_ZX",vtxb.Z(),vtxb.X());
  FillHist("p2_Vertex_ZY",vtxb.Z(),vtxb.Y());
  FillHist("p2_VBDIS",p2->vbdis());
  /* Pi+ and Proton1 */
  vtxb = pip1->vbeam();
  FillHist("pip1_Vertex_XY",vtxb.X(),vtxb.Y());
  FillHist("pip1_Vertex_ZX",vtxb.Z(),vtxb.X());
  FillHist("pip1_Vertex_ZY",vtxb.Z(),vtxb.Y());
  FillHist("pip1_Mass",pip1->vdis());
  FillHist("pip1_VDIS",pip1->vdis());
  FillHist("pip1_VBDIS",pip1->vbdis());
  FillHist("pip1_PBDCA",pip1->pbdca());
  /* Pi+ and Proton2 */
  vtxb = pip2->vbeam();
  FillHist("pip2_Vertex_XY",vtxb.X(),vtxb.Y());
  FillHist("pip2_Vertex_ZX",vtxb.Z(),vtxb.X());
  FillHist("pip2_Vertex_ZY",vtxb.Z(),vtxb.Y());
  FillHist("pip2_Mass",pip2->vdis());
  FillHist("pip2_VDIS",pip2->vdis());
  FillHist("pip2_VBDIS",pip2->vbdis());
  FillHist("pip2_PBDCA",pip2->pbdca());
  /* 2D Plots */
  FillHist("p1_VBDISvsp2_VBDIS",p1->vbdis(),p2->vbdis());
  FillHist("pip1_VDISvspip2_VDIS",pip1->vdis(),pip2->vdis());
  FillHist("pip1_VBDISvspip2_VBDIS",pip1->vbdis(),pip2->vbdis());
  FillHist("pip1_PBDCAvspip2_PBDCA",pip1->pbdca(),pip2->pbdca());
  /* Proton */
  vtxb = proton->vbeam();
  FillHist("p_Vertex_XY",vtxb.X(),vtxb.Y());
  FillHist("p_Vertex_ZX",vtxb.Z(),vtxb.X());
  FillHist("p_Vertex_ZY",vtxb.Z(),vtxb.Y());
  FillHist("p_VBDIS",proton->vbdis());
  /* Lambda */
  vtxb = lam->vbeam();
  FillHist("lam_Vertex_XY",vtxb.X(),vtxb.Y());
  FillHist("lam_Vertex_ZX",vtxb.Z(),vtxb.X());
  FillHist("lam_Vertex_ZY",vtxb.Z(),vtxb.Y());
  vtxb = lam->vertex();
  FillHist("lam_PVertex_XY",vtxb.X(),vtxb.Y());
  FillHist("lam_PVertex_ZX",vtxb.Z(),vtxb.X());
  FillHist("lam_PVertex_ZY",vtxb.Z(),vtxb.Y());
  FillHist("lam_Mass",lam->mass());
  FillHist("lam_VDIS",lam->vdis());
  FillHist("lam_VBDIS",lam->vbdis());
  FillHist("lam_PBDCA",lam->pbdca());
  /* Not Lambda */
  vtxb = nlam->vbeam();
  FillHist("nlam_Vertex_XY",vtxb.X(),vtxb.Y());
  FillHist("nlam_Vertex_ZX",vtxb.Z(),vtxb.X());
  FillHist("nlam_Vertex_ZY",vtxb.Z(),vtxb.Y());
  vtxb = nlam->vertex();
  FillHist("nlam_PVertex_XY",vtxb.X(),vtxb.Y());
  FillHist("nlam_PVertex_ZX",vtxb.Z(),vtxb.X());
  FillHist("nlam_PVertex_ZY",vtxb.Z(),vtxb.Y());
  FillHist("nlam_Mass",nlam->mass());
  FillHist("nlam_VDIS",nlam->vdis());
  FillHist("nlam_VBDIS",nlam->vbdis());
  FillHist("nlam_PBDCA",nlam->pbdca());

  FillHist("lam_PDFvslam_Mass",pdf,lam->mass());
  FillHist("lam_PDFvsnlam_Mass",pdf,nlam->mass());
  FillHist("lam_Massvsnlam_Mass",lam->mass(),nlam->mass());

  FillHist("pip1p2_VDIS",pip1p2dis);
  FillHist("pip2p1_VDIS",pip2p1dis);
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
    if(nc->pid()==F_Neutron) {
      if(nc->energy()>8.0 && (time>9998||time>nc->time())){
        time = nc->time();
        fnc = it;
      }
    }
  }
  pNC*  nc =0;
  if(fnc!=-1)  nc =particle->nc(fnc);
  //if(nc==0) return true;

  TString pname[8] = {"pip","p","d","t","he","pim","k","o"};

  //-----------------------------------------//
  //--- covariance matrices for KinFitter ---//
  //-----------------------------------------//
  double covVal[6][16] = {
    // ### obtained from (p_meas[j]-p_gene[j])*(p_meas[k]-p_gene[k])
    // ###  using G4-data with TH1F(Form("cov_%d_%d_%d", i, j, k), 100, -cov_MAX, cov_MAX);
    // {beam kaon}, {lambda}, {neutron},
    // {proton}, {proton from lambda}, {pi- from lambda}
    { 5.89783e-05, 3.466e-05, 7.87064e-05, 7.0336e-05,
      3.466e-05, 5.65764e-05, 7.58042e-05, 6.77393e-05,
      7.87064e-05, 7.58042e-05, 5.28858e-05, 4.73357e-05,
      7.0336e-05, 6.77393e-05, 4.73357e-05, 4.23708e-05 },
    { 0.000207257, 0.000179973, 0.000148674, 0.000199125,
      0.000179973, 0.000213038, 0.000155228, 0.000208214,
      0.000148674, 0.000155228, 0.000188014, 0.000144255,
      0.000199125, 0.000208214, 0.000144255, 0.000195619 },
    { 0.000310575, 0.000310411, 0.000313857, 0.000302319,
      0.000310411, 0.000309396, 0.000326896, 0.000312685,
      0.000313857, 0.000326896, 0.000302884, 0.000307791,
      0.000302319, 0.000312685, 0.000307791, 0.000326694 },
    { 0.000249203, 0.000218144, 0.000188326, 0.000237145,
      0.000218144, 0.000258923, 0.000189825, 0.000245481,
      0.000188326, 0.000189825, 0.000213973, 0.000186421,
      0.000237145, 0.000245481, 0.000186421, 0.000232521 },
    { 0.000211854, 0.000185869, 0.000139497, 0.000190273,
      0.000185869, 0.000205783, 0.000155463, 0.000208168,
      0.000139497, 0.000155463, 0.000167076, 0.000140599,
      0.000190273, 0.000208168, 0.000140599, 0.000192585 },
    { 2.38346e-05, 2.27798e-05, 1.05223e-05, 1.19754e-05,
      2.27798e-05, 1.52574e-05, 1.0513e-05, 1.84231e-05,
      1.05223e-05, 1.0513e-05, 1.79304e-05, 8.66373e-06,
      1.19754e-05, 1.84231e-05, 8.66373e-06, 1.05972e-05 }
  };
  TMatrixD *covZero = new TMatrixD(4, 4);
  covZero->Zero();
  covZero->ResizeTo(3, 3); // resize from 4x4 to 3x3
  TMatrixD *covParticle[6];
  for( int i=0; i<6; i++ ){
    covParticle[i] = new TMatrixD(4, 4);
    int n = 0;
    for( int j=0; j<4; j++ ){
      for( int k=0; k<4; k++ ){
	if( j==k ) (*covParticle[i])[j][k] = covVal[i][n]; // only diagonal elements
	else       (*covParticle[i])[j][k] = 0;
	n++;
      }
    }
    covParticle[i]->ResizeTo(3, 3); // resize from 4x4 to 3x3
    covParticle[i]->Print(); // Print all
  }
  //-----------------------------------------//
  //--- covariance matrices for KinFitter ---//
  //-----------------------------------------//

  // Target
  TVector3 ZeroV;
  TLorentzVector Ltgt; Ltgt.SetVectM(ZeroV,ThreeHeMass);
  TLorentzVector Lbeam = beam->GetLorentzVector(vertex);

  // Pi+, Pi-, and Pi+ Pi- pair //
  TLorentzVector Lpim = pim->GetLorentzVector();
  TLorentzVector Llam = lam->GetLorentzVector();
  TLorentzVector Lp = proton->GetLorentzVector();
  TLorentzVector Lnp = nproton->GetLorentzVector();
  TLorentzVector Lpipp = Llam + Lp;
  TLorentzVector Lmn = Ltgt+Lbeam-Lpipp;

  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
  // %%% Kinematical Fit using KinFitter %%% //
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
  //--- set TLorentzVector ---//
  // beam_K(K+), L, n, p, p from L, pi- from L
  //  = TLorentzVector L3_beam, L_Llab, L_nlab, L_plab, L_pL, L_piL
  TLorentzVector TL_meas[6]; // measured
  TLorentzVector TL_kfit[6]; // kinematical fitted
  TL_meas[0] = Lbeam;
  TL_meas[1] = Llam; //L_Llab;
  TL_meas[2] = Lmn; //L_nlab;
  TL_meas[3] = Lp;
  TL_meas[4] = Lnp;
  TL_meas[5] = Lpim;
  // L3_target is defined as (0, 0, 0, M_3He)
  TVector3 TV_target = Ltgt.Vect();
  TVector3 TV_meas[6];
  for( int i=0; i<6; i++ ){
    TV_meas[i] = TL_meas[i].Vect();
  }

  TDatabasePDG *pdg = new TDatabasePDG();
  pdg->ReadPDGTable("pdg_table.txt");
  int PDG[6] = {321, 3122, 2112, 2212, 2212, -211};

  //--- KinFitter :: initialization ---//
  //*** definition of fit particles in cartesian coordinates ***//
  TString str_particle[6] = {"Lbeam", "Llam", "Lmn", "Lp", "Lnp", "Lpim"};
  TFitParticlePxPyPz ParticleTgt = TFitParticlePxPyPz("target", "target", &TV_target,
      pdg->GetParticle("He3")->Mass(), covZero);
  TFitParticlePxPyPz Particle[6];
  for( int i=0; i<6; i++ ){
    Particle[i] = TFitParticlePxPyPz(str_particle[i], str_particle[i], &TV_meas[i],
        pdg->GetParticle(PDG[i])->Mass(), covParticle[i]);
  }
  //*** definition of constraints ***//
  // constraint :: mass of Lambda
  TFitConstraintM ConstML = TFitConstraintM("M_L", "M_L", 0, 0, pdg->GetParticle(PDG[1])->Mass());
  ConstML.addParticles1(&Particle[4], &Particle[5]);
  // constraint :: 4-momentum conservation
  TFitConstraintEp ConstEp[4];
  TString str_constEp[4]  = {"Px", "Py", "Pz", "E"};
  for( int i=0; i<4; i++ ){
    ConstEp[i] = TFitConstraintEp(str_constEp[i], str_constEp[i], 0, TFitConstraintEp::component(i), 0);
    ConstEp[i].addParticles1(&ParticleTgt, &Particle[0]);
    ConstEp[i].addParticles2(&Particle[2], &Particle[3], &Particle[4], &Particle[5]);
  }

  //--- KinFitter :: execution ---//
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // Naively,
  //  measured values are momenta (3-vectors) of K-, p, p, pi- -> 3*4=12 
  //  constraints for kinematical fit are masses of lambda and missing-neutron -> 1*2=2
  //    where energy and momentum of missing-neutron is obtained from 4-momentum conservation of K- 3He -> L p n
  //   => number of parameters is 12-2=10
  //   => DOF is 12-10=2
  //
  // In the kinematical fit routine, KinFitter,
  //  fitting values are momenta (3-vectors) of K-, p, p, pi-, n -> 3*5=15
  //    where mass of neutron is fixed to the PDG value
  //  constraints for kinematical fit are mass of lambda and 4-momentum conservation -> 1+4=5
  //   => number of parameters is 15-5=10
  //  measured values are momenta (3-vectors) of K-, p, p, pi- -> 3*4=12 
  //   => DOF is 12-10=2
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //*** definition of the fitter ***//
  TKinFitter kinfitter;
  // add measured particles
  kinfitter.addMeasParticles(&Particle[0], &Particle[3], &Particle[4], &Particle[5]); // K, p, p, pi-
  kinfitter.addUnmeasParticles(&Particle[2]); // n
  // add constraints
  kinfitter.addConstraint(&ConstML); // mass of Lambda
  for( int i=0; i<4; i++ ){
    kinfitter.addConstraint(&ConstEp[i]); // 4-momentum conservation
  }
  //*** perform the fit ***//
  kinfitter.setMaxNbIter(50);       // max number of iterations
  kinfitter.setMaxDeltaS(5e-5);     // max delta chi2
  kinfitter.setMaxF(1e-4);          // max sum of constraints
  //kinfitter.setVerbosity(KFDEBUG);  // verbosity level
  kinfitter.fit();
  //*** copy fit results ***//
  for( int i=0; i<6; i++ ){
    TL_kfit[i] = (*Particle[i].getCurr4Vec());
  }
  TL_kfit[1] = TL_kfit[4]+TL_kfit[5];
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
  // %%% Kinematical Fit using KinFitter %%% //
  // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
  //Lbeam = TL_kfit[0];
  //Llam = TL_kfit[1];
  //Lmn = TL_kfit[2];
  //Lp = TL_kfit[3];
  //Lnp = TL_kfit[4];
  //Lpim = TL_kfit[5];

  // beam_K(K+), L, n, p, p from L, pi- from L
  //  = TLorentzVector L3_beam, L_Llab, L_nlab, L_plab, L_pL, L_piL
  // Target
  TVector3 boost = (Ltgt+Lbeam).BoostVector();
  TLorentzVector cLtgt = Ltgt; cLtgt.Boost(-boost); TLorentzVector cLbeam = Lbeam; cLbeam.Boost(-boost);

  // Pi+, Pi-, and Pi+ Pi- pair //
  TLorentzVector cLpim = pim->GetLorentzVector(); cLpim.Boost(-boost);
  TLorentzVector cLlam = lam->GetLorentzVector(); cLlam.Boost(-boost);
  TLorentzVector cLp = proton->GetLorentzVector(); cLp.Boost(-boost);
  TLorentzVector cLnp = nproton->GetLorentzVector(); cLnp.Boost(-boost);
  TLorentzVector cLpipp = cLlam + cLp;
  TLorentzVector cLmn = Lmn; cLmn.Boost(-boost);
  TLorentzVector Lpipn = Llam + Lmn, cLpipn = cLlam + cLmn;


  TLorentzVector Lp1 = Lp, cLp1 = cLp;
  TLorentzVector Lp2 = Lnp, cLp2 = cLnp;
  TLorentzVector Lpip1 = Lp1 + Lpim, Lpip2 = Lp2 + Lpim;
  TLorentzVector cLpip1 = cLp1 + cLpim, cLpip2 = cLp2 + cLpim;
  TLorentzVector Lnlam = Lp1 + cLpim, cLnlam  = cLp1 + cLpim;
  TLorentzVector Lnpipp = Lnlam + Lp2, cLnpipp = cLnlam + cLp2;


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

  double lam_mass[2]  = { Llam.M()                           , cLlam.M()};
  double lam_mom[2]   = { Llam.Vect().Mag()                  , cLlam.Vect().Mag()};
  double lam_cost[2]  = { Llam.Vect().CosTheta()             , cLlam.Vect().CosTheta()};
  double lam_phi[2]   = { Llam.Vect().Phi()                  , cLlam.Vect().Phi()};
  double lam_mmass[2] = { (Ltgt+Lbeam-Llam).M()              , (cLtgt+cLbeam-cLlam).M()};
  double lam_mmom[2]  = { (Ltgt+Lbeam-Llam).Vect().Mag()     , (cLtgt+cLbeam-cLlam).Vect().Mag()};
  double lam_mcost[2] = { (Ltgt+Lbeam-Llam).Vect().CosTheta(), (cLtgt+cLbeam-cLlam).Vect().CosTheta()};

  double nlam_mass[2]  = { Lnlam.M()                           , cLnlam.M()};
  double nlam_mom[2]   = { Lnlam.Vect().Mag()                  , cLnlam.Vect().Mag()};
  double nlam_cost[2]  = { Lnlam.Vect().CosTheta()             , cLnlam.Vect().CosTheta()};
  double nlam_phi[2]   = { Lnlam.Vect().Phi()                  , cLnlam.Vect().Phi()};
  double nlam_mmass[2] = { (Ltgt+Lbeam-Lnlam).M()              , (cLtgt+cLbeam-cLnlam).M()};
  double nlam_mmom[2]  = { (Ltgt+Lbeam-Lnlam).Vect().Mag()     , (cLtgt+cLbeam-cLnlam).Vect().Mag()};
  double nlam_mcost[2] = { (Ltgt+Lbeam-Lnlam).Vect().CosTheta(), (cLtgt+cLbeam-cLnlam).Vect().CosTheta()};

  double pipp_mass[2]  = { Lpipp.M()                           , cLpipp.M()};
  double pipp_mom[2]   = { Lpipp.Vect().Mag()                  , cLpipp.Vect().Mag()};
  double pipp_cost[2]  = { Lpipp.Vect().CosTheta()             , cLpipp.Vect().CosTheta()};
  double pipp_phi[2]   = { Lpipp.Vect().Phi()                  , cLpipp.Vect().Phi()};
  double pipp_mmass[2] = { (Ltgt+Lbeam-Lpipp).M()              , (cLtgt+cLbeam-cLpipp).M()};
  double pipp_mmom[2]  = { (Ltgt+Lbeam-Lpipp).Vect().Mag()     , (cLtgt+cLbeam-cLpipp).Vect().Mag()};
  double pipp_mcost[2] = { (Ltgt+Lbeam-Lpipp).Vect().CosTheta(), (cLtgt+cLbeam-cLpipp).Vect().CosTheta()};

  double npipp_mass[2]  = { Lnpipp.M()                           , cLnpipp.M()};
  double npipp_mom[2]   = { Lnpipp.Vect().Mag()                  , cLnpipp.Vect().Mag()};
  double npipp_cost[2]  = { Lnpipp.Vect().CosTheta()             , cLnpipp.Vect().CosTheta()};
  double npipp_phi[2]   = { Lnpipp.Vect().Phi()                  , cLnpipp.Vect().Phi()};
  double npipp_mmass[2] = { (Ltgt+Lbeam-Lnpipp).M()              , (cLtgt+cLbeam-cLnpipp).M()};
  double npipp_mmom[2]  = { (Ltgt+Lbeam-Lnpipp).Vect().Mag()     , (cLtgt+cLbeam-cLnpipp).Vect().Mag()};
  double npipp_mcost[2] = { (Ltgt+Lbeam-Lnpipp).Vect().CosTheta(), (cLtgt+cLbeam-cLnpipp).Vect().CosTheta()};

  double pipn_mass[2]  = { Lpipn.M()                           , cLpipn.M()};
  double pipn_mom[2]   = { Lpipn.Vect().Mag()                  , cLpipn.Vect().Mag()};
  double pipn_cost[2]  = { Lpipn.Vect().CosTheta()             , cLpipn.Vect().CosTheta()};
  double pipn_phi[2]   = { Lpipn.Vect().Phi()                  , cLpipn.Vect().Phi()};
  double pipn_mmass[2] = { (Ltgt+Lbeam-Lpipn).M()              , (cLtgt+cLbeam-cLpipn).M()};
  double pipn_mmom[2]  = { (Ltgt+Lbeam-Lpipn).Vect().Mag()     , (cLtgt+cLbeam-cLpipn).Vect().Mag()};
  double pipn_mcost[2] = { (Ltgt+Lbeam-Lpipn).Vect().CosTheta(), (cLtgt+cLbeam-cLpipn).Vect().CosTheta()};

  //double q[2] = {(Lbeam.Vect()-(Ltgt+Lbeam-Lpipp).Vect()).Mag(), (cLbeam.Vect()-(cLtgt+cLbeam-cLpipp).Vect()).Mag()};
  double q[2] = {(Lbeam.Vect()-Lmn.Vect()).Mag(), (cLbeam.Vect()-cLmn.Vect()).Mag()};
  //double tn = (cLtgt+cLbeam-cLpipp).E() - 0.939565;
  //double tn = (cLmn).E() - 0.939565;
  //double tp = cLp.E() - 0.938272;
  //double tl = (cLlam).E() - 1.115683;
  double tn = cLmn.E() - cLmn.M();
  double tp = cLp.E() - cLp.M();
  double tl = cLlam.E() - cLlam.M();
  double qvalue = tn+tp+tl;
  //double qvalue = (Lbeam+Ltgt).E() - 1.115683 - 0.938272 - 0.939565;
  double cos_Lp = (cLp.Vect()*cLlam.Vect())/(cLp.Vect().Mag())/(cLlam.Vect().Mag());

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
  for(int f=1; f<2; f++){
    FillHist(Form("plam_Mass_%s_%02d",frame[f].Data(),reduction),plam_mc_mass[f]);
    FillHist(Form("plam_Momentum_%s_%02d",frame[f].Data(),reduction),plam_mc_mom[f]);
    FillHist(Form("plam_CosT_%s_%02d",frame[f].Data(),reduction),plam_mc_cost[f]);
    FillHist(Form("plam_Phi_%s_%02d",frame[f].Data(),reduction),plam_mc_phi[f]);
    FillHist(Form("plam_MMass_%s_%02d",frame[f].Data(),reduction),plam_mc_mmass[f]);
    FillHist(Form("plam_MMomentum_%s_%02d",frame[f].Data(),reduction),plam_mc_mmom[f]);
    FillHist(Form("plam_MCosT_%s_%02d",frame[f].Data(),reduction),plam_mc_mcost[f]);
    FillHist(Form("n_Mass_%s_%02d",frame[f].Data(),reduction),n_mc_mass[f]);
    FillHist(Form("n_Momentum_%s_%02d",frame[f].Data(),reduction),n_mc_mom[f]);
    FillHist(Form("n_CosT_%s_%02d",frame[f].Data(),reduction),n_mc_cost[f]);
    FillHist(Form("n_Phi_%s_%02d",frame[f].Data(),reduction),n_mc_phi[f]);
    FillHist(Form("n_MMass_%s_%02d",frame[f].Data(),reduction),n_mc_mmass[f]);
    FillHist(Form("n_MMomentum_%s_%02d",frame[f].Data(),reduction),n_mc_mmom[f]);
    FillHist(Form("n_MCosT_%s_%02d",frame[f].Data(),reduction),n_mc_mcost[f]);
    FillHist(Form("lam_plam_Mass_%s_%02d",frame[f].Data(),reduction),lam_plam_mc_mass[f]);
    FillHist(Form("lam_plam_Momentum_%s_%02d",frame[f].Data(),reduction),lam_plam_mc_mom[f]);
    FillHist(Form("lam_plam_CosT_%s_%02d",frame[f].Data(),reduction),lam_plam_mc_cost[f]);
    FillHist(Form("lam_plam_Phi_%s_%02d",frame[f].Data(),reduction),lam_plam_mc_phi[f]);
		FillHist(Form("pipp_IMvsMomentumTransfer_%s_%02d",frame[f].Data(),reduction),plam_mc_mass[f],q_mc[f]);
		FillHist(Form("pipp_IMvsHeKpipp_MCosT_%s_%02d",frame[f].Data(),reduction),plam_mc_mass[f],plam_mc_mcost[f]);
  }
  reduction++;
  //if(mnflag){ return true; }


  return true;

}

bool MyAnalysisHeKpippMCReaction::FillHist(TString name, double val1, int weight)
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

bool MyAnalysisHeKpippMCReaction::FillHist(TString name, TString val1, int weight)
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

bool MyAnalysisHeKpippMCReaction::FillHist(TString name, double val1, double val2, int weight)
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

bool MyAnalysisHeKpippMCReaction::FillHist(TString name, TString val1, TString val2, int weight)
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

void MyAnalysisHeKpippMCReaction::CutCondition()
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

bool MyAnalysisHeKpippMCReaction::Initialize(ConfMan* confMan)
{
  std::cout << "### MyAnalysisHeKpippMCReaction::Initialize ###" << std::endl;

  std::string ofname = confMan->GetOutFileName();
  ofname.insert(ofname.find(".root"),"_anaHeKpippMCReaction");

  rtFile =  new TFile( Form("%s",ofname.c_str()), "RECREATE");
  rtFile -> cd();

  // PDf
  new TH1F( "lam_PDF", "PDF of #Lambda;PDF;Coutns", 400, 0, 40);
  new TH1F( "nlam_PDF", "PDF of Not-#Lambda;PDF;Coutns", 400, 0, 40);
  new TH1F( "pip1_PDF", "PDF of #pi^{-}p_{1};PDF;Coutns", 400, 0, 40);
  new TH1F( "pip2_PDF", "PDF of #pi^{-}p_{2};PDF;Coutns", 400, 0, 40);
  new TH2F( "pip1_PDFvspip2_PDF", "PDF of #pi^{-}p_{1} vs #pi^{-}p_{2};PDF;PDF", 400, 0, 40, 400, 0, 40);
  new TH1F( "pip1_PDF_SemiSelected", "PDF of #pi^{-}p_{1};PDF;Coutns", 400, 0, 40);
  new TH1F( "pip2_PDF_SemiSelected", "PDF of #pi^{-}p_{2};PDF;Coutns", 400, 0, 40);
  new TH2F( "pip1_PDFvspip2_PDF_SemiSelected", "PDF of #pi^{-}p_{1} vs #pi^{-}p_{2};PDF;PDF", 400, 0, 40, 400, 0, 40);
  new TH1F( "lam_PDF_Selected", "PDF of #Lambda;PDF;Coutns", 400, 0, 40);
  new TH1F( "nlam_PDF_Selected", "PDF of Not-#Lambda;PDF;Coutns", 400, 0, 40);
  new TH1F( "pip1_PDF_Selected", "PDF of #pi^{-}p_{1};PDF;Coutns", 400, 0, 40);
  new TH1F( "pip2_PDF_Selected", "PDF of #pi^{-}p_{2};PDF;Coutns", 400, 0, 40);
  new TH2F( "pip1_PDFvspip2_PDF_Selected", "PDF of #pi^{-}p_{1} vs #pi^{-}p_{2};PDF;PDF", 400, 0, 40, 400, 0, 40);
  new TH2F("lam_PDFvslam_Mass","PDF vs Invariant mass of #Lambda;PDF;IM(#Lambda) (GeV/c^{2})", 400,0,40,1000, 1.0, 2.0);
  new TH2F("lam_PDFvsnlam_Mass","PDF vs Invariant mass of Not #Lambda;PDF;IM(#Lambda) (GeV/c^{2})", 400,0,40,1000, 1.0, 2.0);
  new TH2F("lam_Massvsnlam_Mass","Invariant mass of #Lambda vs. Not #Lambda;IM(#Lambda) (GeV/c^{2});IM(Not #Lambda) (GeV/c^{2})", 1000, 1.0, 2.0, 1000, 1.0, 2.0);
  // DCA
  std::cout << "Define Histograms for DCA" << std::endl;
  new TH1F("lam_Mass","Invariant mass of #Lambda;IM(#Lambda) (GeV/c^{2});Counts", 1000, 1.0, 2.0);
  new TH1F("nlam_Mass","Invariant mass of Not #Lambda;IM(#Lambda) (GeV/c^{2});Counts", 1000, 1.0, 2.0);
  new TH1F( "pim_VBDIS", "VBDIS of #pi^{-};VBDIS (cm);Coutns", 1000, 0, 10);
  new TH1F( "p1_VBDIS", "VBDIS of p_{1};VBDIS (cm);Coutns", 1000, 0, 10);
  new TH1F( "p2_VBDIS", "VBDIS of p_{2};VBDIS (cm);Coutns", 1000, 0, 10);
  new TH1F( "p_VBDIS", "VBDIS of p;VBDIS (cm);Coutns", 1000, 0, 10);
  new TH2F( "p1_VBDISvsp2_VBDIS", "VBDIS of p_{1} vs p_{2};VBDIS (cm);VBDIS (cm)", 100, 0, 10, 100, 0, 10);
  new TH1F( "pip1_VDIS", "VDIS of #pi^{-}p_{1};VDIS (cm);Coutns", 1000, 0, 10);
  new TH1F( "pip2_VDIS", "VDIS of #pi^{-}p_{2};VDIS (cm);Coutns", 1000, 0, 10);
  new TH1F( "lam_VDIS", "VDIS of #Lambda;VDIS (cm);Coutns", 1000, 0, 10);
  new TH1F( "nlam_VDIS", "VDIS of Not #Lambda;VDIS (cm);Coutns", 1000, 0, 10);
  new TH2F( "pip1_VDISvspip2_VDIS", "VDIS of #pi^{-}p_{1} vs #pi^{-}p_{2};VDIS (cm);VDIS (cm)", 100, 0, 10, 10, 0, 10);
  new TH1F( "pip1_VBDIS", "VBDIS of #pi^{-}p_{1};VBDIS (cm);Coutns", 1000, 0, 10);
  new TH1F( "pip2_VBDIS", "VBDIS of #pi^{-}p_{2};VBDIS (cm);Coutns", 1000, 0, 10);
  new TH1F( "lam_VBDIS", "VBDIS of #Lambda;VBDIS (cm);Coutns", 1000, 0, 10);
  new TH1F( "nlam_VBDIS", "VBDIS of Not #Lambda;VBDIS (cm);Coutns", 1000, 0, 10);
  new TH2F( "pip1_VBDISvspip2_VBDIS", "VBDIS of #pi^{-}p_{1} vs #pi^{-}p_{2};VBDIS (cm);VBDIS (cm)", 100, 0, 10, 100, 0, 10);
  new TH1F( "pip1_PBDCA", "PBDCA of #pi^{-}p_{1};PBDCA (cm);Coutns", 1000, 0, 10);
  new TH1F( "pip2_PBDCA", "PBDCA of #pi^{-}p_{2};PBDCA (cm);Coutns", 1000, 0, 10);
  new TH1F( "lam_PBDCA", "PBDCA of #Lambda;PBDCA (cm);Coutns", 1000, 0, 10);
  new TH1F( "nlam_PBDCA", "PBDCA of Not #Lambda;PBDCA (cm);Coutns", 1000, 0, 10);
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
  new TH2F( "nlam_Vertex_XY", "Vertex XY plane Not #Lambda;X (cm);Y (cm)", 300, -15, 15, 300, -15, 15);
  new TH2F( "nlam_Vertex_ZX", "Vertex ZX plane Not #Lambda;Z (cm);X (cm)", 600, -30, 30, 300, -15, 15);
  new TH2F( "nlam_Vertex_ZY", "Vertex ZY plane Not #Lambda;Z (cm);Y (cm)", 600, -30, 30, 300, -15, 15);
  new TH2F( "nlam_PVertex_XY", "Vertex XY plane Not #Lambda;X (cm);Y (cm)", 300, -15, 15, 300, -15, 15);
  new TH2F( "nlam_PVertex_ZX", "Vertex ZX plane Not #Lambda;Z (cm);X (cm)", 600, -30, 30, 300, -15, 15);
  new TH2F( "nlam_PVertex_ZY", "Vertex ZY plane Not #Lambda;Z (cm);Y (cm)", 600, -30, 30, 300, -15, 15);
	new TH1F( "pip1p2_VDIS", "VDIS of pip1-p2;VDIS (cm);Coutns", 1000, 0, 10);
	new TH1F( "pip2p1_VDIS", "VDIS of pip2-p1;VDIS (cm);Coutns", 1000, 0, 10);


  // Lambda
  new TH1F("lam_Mass_CM","Invariant mass of #Lambda;IM(#Lambda) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("lam_Momentum_CM","Momentum of #Lambda;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("lam_CosT_CM","cos#theta of #Lambda;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("lam_Phi_CM","#phi of #Lambda;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("lam_MMass_CM", "^{3}He(K^{-},#Lambda)X missing mass;MM(#Lambda) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("lam_MMomentum_CM", "^{3}He(K^{-},#Lambda)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("lam_MCosT_CM", "^{3}He(K^{-},#Lambda)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // proton
  new TH1F("p_Momentum_CM","Momentum of p;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("p_CosT_CM","cos#theta of p;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("p_Phi_CM","#phi of p;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("p_MMass_CM", "^{3}He(K^{-},p)X missing mass;MM(p) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("p_MMomentum_CM", "^{3}He(K^{-},p)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("p_MCosT_CM", "^{3}He(K^{-},p)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // Lambda
  new TH1F("lam_plam_Mass_CM","Invariant mass of #Lambda;IM(#Lambda) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("lam_plam_Momentum_CM","Momentum of #Lambda;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("lam_plam_CosT_CM","cos#theta of #Lambda;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("lam_plam_Phi_CM","#phi of #Lambda;#phi;Counts", 1600, 0.0, 3.2);
  // proton
  new TH1F("p_plam_Momentum_CM","Momentum of p;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("p_plam_CosT_CM","cos#theta of p;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("p_plam_Phi_CM","#phi of p;#phi;Counts", 1600, 0.0, 3.2);
  // neutron
  new TH1F("n_Momentum_CM","Momentum of n;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("n_CosT_CM","cos#theta of n;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("n_Phi_CM","#phi of n;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("n_MMass_CM", "^{3}He(K^{-},n)X missing mass;MM(n) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("n_MMomentum_CM", "^{3}He(K^{-},n)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("n_MCosT_CM", "^{3}He(K^{-},n)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // p + Lambda
  new TH1F("plam_Mass_CM","Invariant mass of #Lambdap;IM(#Lambdap) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("plam_Momentum_CM","Momentum of #Lambdap;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("plam_CosT_CM","cos#theta of #Lambdap;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("plam_Phi_CM","#phi of #Lambdap;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("plam_MMass_CM", "^{3}He(K^{-},#Lambdap)X missing mass;MM(#Lambdap) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("plam_MMomentum_CM", "^{3}He(K^{-},#Lambdap)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("plam_MCosT_CM", "^{3}He(K^{-},#Lambdap)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);

  // Checking reduction fuctor
  for(int i=1; i<30; i++){
    // Lambda in pLam frame
    new TH1F(Form("lam_plam_Mass_CM_%02d",i),"Invariant mass of #Lambda;IM(#Lambda) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
    new TH1F(Form("lam_plam_Momentum_CM_%02d",i),"Momentum of #Lambda;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
    new TH1F(Form("lam_plam_CosT_CM_%02d",i),"cos#theta of #Lambda;cos#theta;Counts", 3000, -1.5, 1.5);
    new TH1F(Form("lam_plam_Phi_CM_%02d",i),"#phi of #Lambda;#phi;Counts", 1600, 0.0, 3.2);
    // p + Lambda
    new TH1F(Form("plam_Mass_CM_%02d",i),"Invariant mass of #Lambdap;IM(#Lambdap) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
    new TH1F(Form("plam_Momentum_CM_%02d",i),"Momentum of #Lambdap;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
    new TH1F(Form("plam_CosT_CM_%02d",i),"cos#theta of #Lambdap;cos#theta;Counts", 3000, -1.5, 1.5);
    new TH1F(Form("plam_Phi_CM_%02d",i),"#phi of #Lambdap;#phi;Counts", 1600, 0.0, 3.2);
    new TH1F(Form("plam_MMass_CM_%02d",i), "^{3}He(K^{-},#Lambdap)X missing mass;MM(#Lambdap) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
    new TH1F(Form("plam_MMomentum_CM_%02d",i), "^{3}He(K^{-},#Lambdap)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
    new TH1F(Form("plam_MCosT_CM_%02d",i), "^{3}He(K^{-},#Lambdap)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
    // neutron
    new TH1F(Form("n_Mass_CM_%02d",i),"Invariant mass of n;IM(n) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
    new TH1F(Form("n_Momentum_CM_%02d",i),"Momentum of n;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
    new TH1F(Form("n_CosT_CM_%02d",i),"cos#theta of n;cos#theta;Counts", 3000, -1.5, 1.5);
    new TH1F(Form("n_Phi_CM_%02d",i),"#phi of n;#phi;Counts", 1600, 0.0, 3.2);
    new TH1F(Form("n_MMass_CM_%02d",i), "^{3}He(K^{-},n)X missing mass;MM(n) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
    new TH1F(Form("n_MMomentum_CM_%02d",i), "^{3}He(K^{-},n)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
		new TH1F(Form("n_MCosT_CM_%02d",i), "^{3}He(K^{-},n)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
		new TH2F(Form("pipp_IMvsMomentumTransfer_CM_%02d",i), "(#pi^{-}p_{1}p_{2}) invariant mass vs. Momentum transfer;IM(#pi^{-}p_{1}p_{2}) (GeV/c^{2});Momentum Transfer (GeV/c)", 200, 2.0, 4.0, 200, 0.0, 2.0);
		new TH2F(Form("pipp_IMvsHeKpipp_MCosT_CM_%02d",i), "(#pi^{-}p_{1}p_{2}) invariant mass vs. cos#theta_{n};IM(#pi^{-}p_{1}p_{2}) (GeV/c^{2});cos#theta_{n}", 200, 2.0, 4.0, 300, -1.5, 1.5);
	}
	new TH2F("pipp_IMvsMomentumTransfer_CM", "(#pi^{-}p_{1}p_{2}) invariant mass vs. Momentum transfer;IM(#pi^{-}p_{1}p_{2}) (GeV/c^{2});Momentum Transfer (GeV/c)", 200, 2.0, 4.0, 200, 0.0, 2.0);
	new TH2F("pipp_IMvsHeKpipp_MCosT_CM", "(#pi^{-}p_{1}p_{2}) invariant mass vs. cos#theta_{n};IM(#pi^{-}p_{1}p_{2}) (GeV/c^{2});cos#theta_{n}", 200, 2.0, 4.0, 300, -1.5, 1.5);

	std::cout << "== Finish Histogram Initialization ==" << std::endl;
	return true;

}
