// MyAnalysisdKnppim.cpp

#include "MyAnalysisdKnppim.h"

MyAnalysisdKnppim::MyAnalysisdKnppim(TFile* rt, ConfMan* conf)
{
  Initialize(conf);
  CutCondition();
  Clear();
}

MyAnalysisdKnppim::~MyAnalysisdKnppim()
{
  Clear();
  rtFile->cd();
  rtFile->Write();
  rtFile->Close();
}

void MyAnalysisdKnppim::Clear()
{
}

bool MyAnalysisdKnppim::DoAnalysis(ConfMan* conf, EventHeader* header, BeamLineHitMan* blMan, BeamLineTrackMan* bltrackMan, CDSHitMan* cdsMan, CDSTrackingMan* cdstrackMan, Particle* particle)
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
  if(MulBVC!=0||MulCVC!=0||MulPC!=0) return true;

  if(particle->nCDS()!=2) return true;
  FillHist("EventNumber",2);
  if(particle->nProton()!=1) return true;
  if(particle->nPiminus()!=1) return true;
  FillHist("EventNumber",3);

  pCDS* proton = particle->proton(0);
  pCDS* pim = particle->pim(0);
  pCDS* pro = particle->product(0);

  // ############### //
  // Vertex decision //
  // ############### //
  TVector3 vertex;
  double vdis = 9999;
  bool FIDUCIAL = false;
  vdis = pro->vdis();
  vertex = pro->vbeam();
  if(GeomTools::GetID(vertex)==CID_Fiducial){ FIDUCIAL = true; }

  if(!FIDUCIAL){ return true; }
  FillHist("EventNumber",4);

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
  FillHist("EventNumber",5);

  // Target
  TString frame[2] = {"Lab","CM"};
  TString pname[8] = {"pip","p","d","t","he","pim","k","o"};
  TVector3 ZeroV;
  TLorentzVector Ltgt; Ltgt.SetVectM(ZeroV,dMass);
  TLorentzVector Lbeam = beam->GetLorentzVector(vertex);
  TVector3 boost = (Ltgt+Lbeam).BoostVector();
  TLorentzVector cLtgt = Ltgt; cLtgt.Boost(-boost); TLorentzVector cLbeam = Lbeam; cLbeam.Boost(-boost);

  // Pi+, Pi-, and Pi+ Pi- pair //
  TLorentzVector Lproton = proton->GetLorentzVector();
  TLorentzVector cLproton = proton->GetLorentzVector(); cLproton.Boost(-boost);
  TLorentzVector Lpim = pim->GetLorentzVector();
  TLorentzVector cLpim = pim->GetLorentzVector(); cLpim.Boost(-boost);
  TLorentzVector Lppim = Lproton + Lpim;
  TLorentzVector cLppim = cLproton + cLpim;
  double proton_mass[2]  = { Lproton.M()                           , cLproton.M()};
  double proton_mom[2]   = { Lproton.Vect().Mag()                  , cLproton.Vect().Mag()};
  double proton_cost[2]  = { Lproton.Vect().CosTheta()             , cLproton.Vect().CosTheta()};
  double proton_phi[2]   = { Lproton.Vect().Phi()                  , cLproton.Vect().Phi()};
  double proton_mmass[2] = { (Ltgt+Lbeam-Lproton).M()              , (cLtgt+cLbeam-cLproton).M()};
  double proton_mmom[2]  = { (Ltgt+Lbeam-Lproton).Vect().Mag()     , (cLtgt+cLbeam-cLproton).Vect().Mag()};
  double proton_mcost[2] = { (Ltgt+Lbeam-Lproton).Vect().CosTheta(), (cLtgt+cLbeam-cLproton).Vect().CosTheta()};

  double pim_mass[2]  = { Lpim.M()                           , cLpim.M()};
  double pim_mom[2]   = { Lpim.Vect().Mag()                  , cLpim.Vect().Mag()};
  double pim_cost[2]  = { Lpim.Vect().CosTheta()             , cLpim.Vect().CosTheta()};
  double pim_phi[2]   = { Lpim.Vect().Phi()                  , cLpim.Vect().Phi()};
  double pim_mmass[2] = { (Ltgt+Lbeam-Lpim).M()              , (cLtgt+cLbeam-cLpim).M()};
  double pim_mmom[2]  = { (Ltgt+Lbeam-Lpim).Vect().Mag()     , (cLtgt+cLbeam-cLpim).Vect().Mag()};
  double pim_mcost[2] = { (Ltgt+Lbeam-Lpim).Vect().CosTheta(), (cLtgt+cLbeam-cLpim).Vect().CosTheta()};

  double ppim_mass[2]  = { Lppim.M()                           , cLppim.M()};
  double ppim_mom[2]   = { Lppim.Vect().Mag()                  , cLppim.Vect().Mag()};
  double ppim_cost[2]  = { Lppim.Vect().CosTheta()             , cLppim.Vect().CosTheta()};
  double ppim_phi[2]   = { Lppim.Vect().Phi()                  , cLppim.Vect().Phi()};
  double ppim_mmass[2] = { (Ltgt+Lbeam-Lppim).M()              , (cLtgt+cLbeam-cLppim).M()};
  double ppim_mmom[2]  = { (Ltgt+Lbeam-Lppim).Vect().Mag()     , (cLtgt+cLbeam-cLppim).Vect().Mag()};
  double ppim_mcost[2] = { (Ltgt+Lbeam-Lppim).Vect().CosTheta(), (cLtgt+cLbeam-cLppim).Vect().CosTheta()};

  double proton_vdis = proton->vdis();
  double pim_vdis = pim->vdis();
  double pro_pbdca = pro->pbdca();
  double pro_vbdis = pro->vbdis();
  double pro_vdis  = pro->vdis();

  TLorentzVector Ln  = nc ->GetLorentzVector();
  TLorentzVector cLn  = nc ->GetLorentzVector(); cLn.Boost(-boost);
  double nc_mass[2]    = { Ln.M()                           , cLn.M()};
  double nc_mom[2]     = { Ln.Vect().Mag()                  , cLn.Vect().Mag()};
  double nc_cost[2]    = { Ln.Vect().CosTheta()             , cLn.Vect().CosTheta()};
  double nc_phi[2]     = { Ln.Vect().Phi()                  , cLn.Vect().Phi()};
  double nc_mmass[2]   = { (Ltgt+Lbeam-Ln).M()              , (cLtgt+cLbeam-cLn).M()};
  double nc_mmom[2]    = { (Ltgt+Lbeam-Ln).Vect().Mag()     , (cLtgt+cLbeam-cLn).Vect().Mag()};
  double nc_mcost[2]   = { (Ltgt+Lbeam-Ln).Vect().CosTheta(), (cLtgt+cLbeam-cLn).Vect().CosTheta()};

  double nproton_mass[2]  = { (Lproton+Ln).M()                           , (cLproton+cLn).M()};
  double nproton_mom[2]   = { (Lproton+Ln).Vect().Mag()                  , (cLproton+cLn).Vect().Mag()};
  double nproton_cost[2]  = { (Lproton+Ln).Vect().CosTheta()             , (cLproton+cLn).Vect().CosTheta()};
  double nproton_phi[2]   = { (Lproton+Ln).Vect().Phi()                  , (cLproton+cLn).Vect().Phi()};
  double nproton_mmass[2] = { (Ltgt+Lbeam-(Lproton+Ln)).M()              , (cLtgt+cLbeam-(cLproton+cLn)).M()};
  double nproton_mmom[2]  = { (Ltgt+Lbeam-(Lproton+Ln)).Vect().Mag()     , (cLtgt+cLbeam-(cLproton+cLn)).Vect().Mag()};
  double nproton_mcost[2] = { (Ltgt+Lbeam-(Lproton+Ln)).Vect().CosTheta(), (cLtgt+cLbeam-(cLproton+cLn)).Vect().CosTheta()};

  double npim_mass[2]  = { (Lpim+Ln).M()                           , (cLpim+cLn).M()};
  double npim_mom[2]   = { (Lpim+Ln).Vect().Mag()                  , (cLpim+cLn).Vect().Mag()};
  double npim_cost[2]  = { (Lpim+Ln).Vect().CosTheta()             , (cLpim+cLn).Vect().CosTheta()};
  double npim_phi[2]   = { (Lpim+Ln).Vect().Phi()                  , (cLpim+cLn).Vect().Phi()};
  double npim_mmass[2] = { (Ltgt+Lbeam-(Lpim+Ln)).M()              , (cLtgt+cLbeam-(cLpim+cLn)).M()};
  double npim_mmom[2]  = { (Ltgt+Lbeam-(Lpim+Ln)).Vect().Mag()     , (cLtgt+cLbeam-(cLpim+cLn)).Vect().Mag()};
  double npim_mcost[2] = { (Ltgt+Lbeam-(Lpim+Ln)).Vect().CosTheta(), (cLtgt+cLbeam-(cLpim+cLn)).Vect().CosTheta()};

  double nppim_mass[2] = { (Lppim+Ln).M()                           , (Lppim+cLn).M()};
  double nppim_mom[2]  = { (Lppim+Ln).Vect().Mag()                  , (Lppim+cLn).Vect().Mag()};
  double nppim_cost[2] = { (Lppim+Ln).Vect().CosTheta()             , (Lppim+cLn).Vect().CosTheta()};
  double nppim_phi[2]  = { (Lppim+Ln).Vect().Phi()                  , (Lppim+cLn).Vect().Phi()};
  double nppim_mmass[2] = { (Ltgt+Lbeam-(Lppim+Ln)).M()              , (cLtgt+cLbeam-(cLppim+cLn)).M()};
  double nppim_mmom[2]  = { (Ltgt+Lbeam-(Lppim+Ln)).Vect().Mag()     , (cLtgt+cLbeam-(cLppim+cLn)).Vect().Mag()};
  double nppim_mcost[2] = { (Ltgt+Lbeam-(Lppim+Ln)).Vect().CosTheta(), (cLtgt+cLbeam-(cLppim+cLn)).Vect().CosTheta()};

  double seg = (nc->seg()-1)%16+1;
  double lay = (nc->seg()-1)/16+1;

  /* flag */
  bool lamflag = false; /* lam flag */
  bool sbllamflag = false; /* lower sideband of lam flag */
  bool sbulamflag = false; /* upper sideband of lam flag */
  {
    double pro_mass = ppim_mass[1];
    if(lamll<pro_mass&&pro_mass<lamul) lamflag = true;
    if(sbllamll<pro_mass&&pro_mass<sbllamul) sbllamflag = true;
    if(sbulamll<pro_mass&&pro_mass<sbulamul) sbulamflag = true;
  }
  bool smflag = false; /* sm flag */
  bool sblsmflag = false; /* lower sideband of sm flag */
  bool sbusmflag = false; /* upper sideband of sm flag */
  {
    double pro_mass = npim_mass[1];
    if(smll<pro_mass&&pro_mass<smul) smflag = true;
    if(sblsmll<pro_mass&&pro_mass<sblsmul) sblsmflag = true;
    if(sbusmll<pro_mass&&pro_mass<sbusmul) sbusmflag = true;
  }

  //if(!smflag){ return true; }
  //if(smflag){ return true; }
  //if(!lamflag){ return true; }
  //if(lamflag){ return true; }

  FillHist("FWD_n_OverbetavsMomentum",1.0/nc->beta(),nc->mom().Mag());
  FillHist("FWD_n_OverbetavsEnergy",1.0/nc->beta(),nc->energy());
  FillHist("FWD_n_TOFvsMomentum",nc->tof(),nc->mom().Mag());
  FillHist("FWD_n_HitPosition",nc->hitpos().X(),nc->hitpos().Y());
  FillHist("FWD_n_HitSegment",seg,lay);

  /*-----------------------------------------*/
  /* All                                     */ 
  /*-----------------------------------------*/
  // proton
  for(int f=1; f<2; f++){
    FillHist(Form("p_CosTvsMomentum_%s",frame[f].Data()),proton_cost[f],proton_mom[f]);
    FillHist(Form("p_CosTvsPhi_%s",frame[f].Data()),proton_cost[f],proton_phi[f]);
    FillHist(Form("dKp_MMvsMomentum_%s",frame[f].Data()),proton_mmass[f],proton_mmom[f]);
    FillHist(Form("dKp_MMvsCosT_%s",frame[f].Data()),proton_mmass[f],proton_mcost[f]);
    FillHist(Form("p_IMvsMomentum_%s",frame[f].Data()),proton_mass[f],proton_mom[f]);
    FillHist(Form("p_IMvsCosT_%s",frame[f].Data()),proton_mass[f],proton_cost[f]);
  }
  // pi-
  for(int f=1; f<2; f++){
    FillHist(Form("pim_CosTvsMomentum_%s",frame[f].Data()),pim_cost[f],pim_mom[f]);
    FillHist(Form("pim_CosTvsPhi_%s",frame[f].Data()),pim_cost[f],pim_phi[f]);
    FillHist(Form("dKpim_MMvsMomentum_%s",frame[f].Data()),pim_mmass[f],pim_mmom[f]);
    FillHist(Form("dKpim_MMvsCosT_%s",frame[f].Data()),pim_mmass[f],pim_mcost[f]);
    FillHist(Form("pim_IMvsMomentum_%s",frame[f].Data()),pim_mass[f],pim_mom[f]);
    FillHist(Form("pim_IMvsCosT_%s",frame[f].Data()),pim_mass[f],pim_cost[f]);
  }
  // neutron
  for(int f=1; f<2; f++){
    FillHist(Form("n_CosTvsMomentum_%s",frame[f].Data()),nc_cost[f],nc_mom[f]);
    FillHist(Form("n_CosTvsPhi_%s",frame[f].Data()),nc_cost[f],nc_phi[f]);
    FillHist(Form("dKn_MMvsMomentum_%s",frame[f].Data()),nc_mmass[f],nc_mmom[f]);
    FillHist(Form("dKn_MMvsCosT_%s",frame[f].Data()),nc_mmass[f],nc_mcost[f]);
    FillHist(Form("n_IMvsMomentum_%s",frame[f].Data()),nc_mass[f],nc_mom[f]);
    FillHist(Form("n_IMvsCosT_%s",frame[f].Data()),nc_mass[f],nc_cost[f]);
  }
  // proton pi-
  for(int f=1; f<2; f++){
    FillHist(Form("ppim_CosTvsMomentum_%s",frame[f].Data()),ppim_cost[f],ppim_mom[f]);
    FillHist(Form("ppim_CosTvsPhi_%s",frame[f].Data()),ppim_cost[f],ppim_phi[f]);
    FillHist(Form("dKppim_MMvsMomentum_%s",frame[f].Data()),ppim_mmass[f],ppim_mmom[f]);
    FillHist(Form("dKppim_MMvsCosT_%s",frame[f].Data()),ppim_mmass[f],ppim_mcost[f]);
    FillHist(Form("ppim_IMvsMomentum_%s",frame[f].Data()),ppim_mass[f],ppim_mom[f]);
    FillHist(Form("ppim_IMvsCosT_%s",frame[f].Data()),ppim_mass[f],ppim_cost[f]);
  }
  // neutron proton
  for(int f=1; f<2; f++){
    FillHist(Form("np_CosTvsMomentum_%s",frame[f].Data()),nproton_cost[f],nproton_mom[f]);
    FillHist(Form("np_CosTvsPhi_%s",frame[f].Data()),nproton_cost[f],nproton_phi[f]);
    FillHist(Form("dKnp_MMvsMomentum_%s",frame[f].Data()),nproton_mmass[f],nproton_mmom[f]);
    FillHist(Form("dKnp_MMvsCosT_%s",frame[f].Data()),nproton_mmass[f],nproton_mcost[f]);
    FillHist(Form("np_IMvsMomentum_%s",frame[f].Data()),nproton_mass[f],nproton_mom[f]);
    FillHist(Form("np_IMvsCosT_%s",frame[f].Data()),nproton_mass[f],nproton_cost[f]);
  }
  // neutron pi-
  for(int f=1; f<2; f++){
    FillHist(Form("npim_CosTvsMomentum_%s",frame[f].Data()),npim_cost[f],npim_mom[f]);
    FillHist(Form("npim_CosTvsPhi_%s",frame[f].Data()),npim_cost[f],npim_phi[f]);
    FillHist(Form("dKnpim_MMvsMomentum_%s",frame[f].Data()),npim_mmass[f],npim_mmom[f]);
    FillHist(Form("dKnpim_MMvsCosT_%s",frame[f].Data()),npim_mmass[f],npim_mcost[f]);
    FillHist(Form("npim_IMvsMomentum_%s",frame[f].Data()),npim_mass[f],npim_mom[f]);
    FillHist(Form("npim_IMvsCosT_%s",frame[f].Data()),npim_mass[f],npim_cost[f]);
  }
  // neutron proton pi-
  for(int f=1; f<2; f++){
    FillHist(Form("nppim_CosTvsMomentum_%s",frame[f].Data()),nppim_cost[f],nppim_mom[f]);
    FillHist(Form("nppim_CosTvsPhi_%s",frame[f].Data()),nppim_cost[f],nppim_phi[f]);
    FillHist(Form("dKnppim_MMvsMomentum_%s",frame[f].Data()),nppim_mmass[f],nppim_mmom[f]);
    FillHist(Form("dKnppim_MMvsCosT_%s",frame[f].Data()),nppim_mmass[f],nppim_mcost[f]);
    FillHist(Form("nppim_IMvsMomentum_%s",frame[f].Data()),nppim_mass[f],nppim_mom[f]);
    FillHist(Form("nppim_IMvsCosT_%s",frame[f].Data()),nppim_mass[f],nppim_cost[f]);
  }

  /*-----------------------------------------*/
  /* Advanced                                */ 
  /*-----------------------------------------*/
  // MM(n) vs. MM(nppim)
  for(int f=1; f<2; f++){
    FillHist(Form("dKn_MMvsdKnppim_MM_%s",frame[f].Data()),nc_mmass[f],nppim_mmass[f]);
  }

  //FillHist(Form("%s_DCA",pname[cds->pid()].Data()),cds->vdis());
  //FillHist(Form("all_DCA"),cds->vdis());

  return true;

}

bool MyAnalysisdKnppim::FillHist(TString name, double val1, int weight)
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

bool MyAnalysisdKnppim::FillHist(TString name, TString val1, int weight)
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

bool MyAnalysisdKnppim::FillHist(TString name, double val1, double val2, int weight)
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

bool MyAnalysisdKnppim::FillHist(TString name, TString val1, TString val2, int weight)
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

void MyAnalysisdKnppim::CutCondition()
{
  /* K0 mass cut */ 
  double mean  = 0.4976;
  double sigma = 0.0069;
  k0ll = mean-2*sigma; k0ul = mean+2*sigma;
  sblk0ll = mean-5*sigma; sblk0ul = mean-3*sigma;
  sbuk0ll = mean+3*sigma; sbuk0ul = mean+5*sigma;
  /* Lambda mass cut */ 
  mean  = 1.1160;
  sigma = 0.0029;
  lamll = mean-2*sigma; lamul = mean+2*sigma;
  sbllamll = mean-5*sigma; sbllamul = mean-3*sigma;
  sbulamll = mean+3*sigma; sbulamul = mean+5*sigma;
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
}

bool MyAnalysisdKnppim::Initialize(ConfMan* confMan)
{
  std::cout << "### MyAnalysisdKnppim::Initialize ###" << std::endl;

  std::string ofname = confMan->GetOutFileName();
  ofname.insert(ofname.find(".root"),"_anadKnppim");

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
  new TH1F( "all_DCA", "DCA;DCA (cm);Coutns", 100, 0, 10);
  new TH1F( "pip_DCA", "DCA of #pi^{+};DCA (cm);Coutns", 100, 0, 10);
  new TH1F( "pim_DCA", "DCA of #pi^{-};DCA (cm);Coutns", 100, 0, 10);
  new TH1F( "k_DCA", "DCA of K^{-};DCA (cm);Coutns", 100, 0, 10);
  new TH1F( "p_DCA", "DCA of p;DCA (cm);Coutns", 100, 0, 10);
  // FWD neutral Particles
  std::cout << "Define Histograms for FWD neutral particles" << std::endl;
  new TH2F( "FWD_OverbetavsMomentum", "1/#beta vs. Momentum;1/#beta;Momentum (GeV/c)", 1000, -0, 10, 500, -5, 5 );
  new TH2F( "FWD_OverbetavsEnergy", "1/#beta vs. Energy deposit;1/#beta;Energy deposit (MeVee)", 1000, -0, 10, 5000, 0, 100);
  new TH2F( "FWD_TOFvsMomentum", "TOF vs. Momentum;TOF (ns);Momentum (GeV/c)", 4000, 30, 70, 500, -5, 5 );
  new TH2F( "FWD_HitPosition", "Hit position at NC (X vs. Y);X position (cm);Y position (cm)",40,-400,400,25,-250,250);
  new TH2F( "FWD_HitSegment", "NC segment;Segment in layer;Layer",16,0.5,16.5,7,0.5,7.5);

  /*-----------------------------------------*/
  /* All                                     */ 
  /*-----------------------------------------*/
  // proton
  new TH2F("p_CosTvsMomentum_CM","cos#theta vs. momentum of p;cos#theta;Momentum (GeV/c)", 2000, -1.0, 1.0, 2000, 0.0, 2.0);
  new TH2F("p_CosTvsPhi_CM","cos#theta vs. #phi of p;cos#theta;#phi", 2000, -1.0, 1.0, 1600, 0.0, 3.2);
  new TH2F("dKp_MMvsMomentum_CM", "d(K^{-},p)X missing mass vs. missing momentum;MM_{d}(p) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("dKp_MMvsCosT_CM", "d(K^{-},p)X missing mass vs. cos#theta;MM_{d}(p) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("p_IMvsMomentum_CM", "(p) invariant mass vs. momentum;IM(p) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("p_IMvsCosT_CM", "(p) invariant mass vs. cos#theta;IM(p) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  // pi-
  new TH2F("pim_CosTvsMomentum_CM","cos#theta vs. momentum of #pi^{-};cos#theta;Momentum (GeV/c)", 2000, -1.0, 1.0, 2000, 0.0, 2.0);
  new TH2F("pim_CosTvsPhi_CM","cos#theta vs. #phi of #pi^{-};cos#theta;#phi", 2000, -1.0, 1.0, 1600, 0.0, 3.2);
  new TH2F("dKpim_MMvsMomentum_CM", "d(K^{-},#pi^{-})X missing mass vs. missing momentum;MM_{d}(#pi^{-}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 1.0, 3.0, 2000, 0.0, 2.0);
  new TH2F("dKpim_MMvsCosT_CM", "d(K^{-},#pi^{-})X missing mass vs. cos#theta;MM_{d}(#pi^{-}) (GeV/c^{2});cos#theta", 2000, 1.0, 3.0, 2000, -1.0, 1.0);
  new TH2F("pim_IMvsMomentum_CM", "(#pi^{-}) invariant mass vs. momentum;IM(#pi^{-}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("pim_IMvsCosT_CM", "(#pi^{-}) invariant mass vs. cos#theta;IM(#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  // neutron
  new TH2F("n_CosTvsMomentum_CM","cos#theta vs. momentum of n;cos#theta;Momentum (GeV/c)", 2000, -1.0, 1.0, 2000, 0.0, 2.0);
  new TH2F("n_CosTvsPhi_CM","cos#theta vs. #phi of n;cos#theta;#phi", 2000, -1.0, 1.0, 1600, 0.0, 3.2);
  new TH2F("dKn_MMvsMomentum_CM", "d(K^{-},n)X missing mass vs. missing momentum;MM_{d}(n) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("dKn_MMvsCosT_CM", "d(K^{-},n)X missing mass vs. cos#theta;MM_{d}(n) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("n_IMvsMomentum_CM", "(n) invariant mass vs. momentum;IM(n) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("n_IMvsCosT_CM", "(n) invariant mass vs. cos#theta;IM(n) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  // proton pi-
  new TH2F("ppim_CosTvsMomentum_CM","cos#theta vs. momentum of p#pi^{-};cos#theta;Momentum (GeV/c)", 2000, -1.0, 1.0, 2000, 0.0, 2.0);
  new TH2F("ppim_CosTvsPhi_CM","cos#theta vs. #phi of p#pi^{-};cos#theta;#phi", 2000, -1.0, 1.0, 1600, 0.0, 3.2);
  new TH2F("dKppim_MMvsMomentum_CM", "d(K^{-},p#pi^{-})X missing mass vs. missing momentum;MM_{d}(p#pi^{-}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("dKppim_MMvsCosT_CM", "d(K^{-},p#pi^{-})X missing mass vs. cos#theta;MM_{d}(p#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("ppim_IMvsMomentum_CM", "(p#pi^{-}) invariant mass vs. momentum;IM(p#pi^{-}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("ppim_IMvsCosT_CM", "(p#pi^{-}) invariant mass vs. cos#theta;IM(p#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  // neutron proton
  new TH2F("np_CosTvsMomentum_CM","cos#theta vs. momentum of np;cos#theta;Momentum (GeV/c)", 2000, -1.0, 1.0, 2000, 0.0, 2.0);
  new TH2F("np_CosTvsPhi_CM","cos#theta vs. #phi of np;cos#theta;#phi", 2000, -1.0, 1.0, 1600, 0.0, 3.2);
  new TH2F("dKnp_MMvsMomentum_CM", "d(K^{-},np)X missing mass vs. missing momentum;MM_{d}(np) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("dKnp_MMvsCosT_CM", "d(K^{-},np)X missing mass vs. cos#theta;MM_{d}(np) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("np_IMvsMomentum_CM", "(np) invariant mass vs. momentum;IM(np) (GeV/c^{2});Momentum (GeV/c)", 2000, 1.0, 3.0, 2000, 0.0, 2.0);
  new TH2F("np_IMvsCosT_CM", "(np) invariant mass vs. cos#theta;IM(np) (GeV/c^{2});cos#theta", 2000, 1.0, 3.0, 2000, -1.0, 1.0);
  // neutron pi-
  new TH2F("npim_CosTvsMomentum_CM","cos#theta vs. momentum of n#pi^{-};cos#theta;Momentum (GeV/c)", 2000, -1.0, 1.0, 2000, 0.0, 2.0);
  new TH2F("npim_CosTvsPhi_CM","cos#theta vs. #phi of n#pi^{-};cos#theta;#phi", 2000, -1.0, 1.0, 1600, 0.0, 3.2);
  new TH2F("dKnpim_MMvsMomentum_CM", "d(K^{-},n#pi^{-})X missing mass vs. missing momentum;MM_{d}(n#pi^{-}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("dKnpim_MMvsCosT_CM", "d(K^{-},n#pi^{-})X missing mass vs. cos#theta;MM_{d}(n#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("npim_IMvsMomentum_CM", "(n#pi^{-}) invariant mass vs. momentum;IM(n#pi^{-}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("npim_IMvsCosT_CM", "(n#pi^{-}) invariant mass vs. cos#theta;IM(n#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  // neutron proton pi-
  new TH2F("nppim_CosTvsMomentum_CM","cos#theta vs. momentum of np#pi^{-};cos#theta;Momentum (GeV/c)", 2000, -1.0, 1.0, 2000, 0.0, 2.0);
  new TH2F("nppim_CosTvsPhi_CM","cos#theta vs. #phi of np#pi^{-};cos#theta;#phi", 2000, -1.0, 1.0, 1600, 0.0, 3.2);
  new TH2F("dKnppim_MMvsMomentum_CM", "d(K^{-},np#pi^{-})X missing mass vs. missing momentum;MM_{d}(np#pi^{-}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("dKnppim_MMvsCosT_CM", "d(K^{-},nppi^{-})X missing mass vs. cos#theta;MM_{d}(np#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("nppim_IMvsMomentum_CM", "(np#pi^{-}) invariant mass vs. momentum;IM(np#pi^{-}) (GeV/c^{2});Momentum (GeV/c)", 2000, 1.0, 3.0, 2000, 0.0, 2.0);
  new TH2F("nppim_IMvsCosT_CM", "(np#pi^{-}) invariant mass vs. cos#theta;IM(np#pi^{-}) (GeV/c^{2});cos#theta", 2000, 1.0, 3.0, 2000, -1.0, 1.0);

  new TH2F("dKn_MMvsdKnppim_MM_CM", "d(K^{-},n)X missing mass vs. d(K^{-},np#pi^{-})X missing mass;MM_{d}(n) (GeV/c^{2});MM_{d}(np#pi^{-}) (GeV/c^{2});", 2000, 0.0, 2.0, 2000, 0.0, 2.0);

  return true;

}
