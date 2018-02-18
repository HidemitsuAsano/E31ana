// MyAnalysisdKnpipi.cpp

#include "MyAnalysisdKnpipi.h"

MyAnalysisdKnpipi::MyAnalysisdKnpipi(TFile* rt, ConfMan* conf)
{
  Initialize(conf);
  CutCondition();
  Clear();
}

MyAnalysisdKnpipi::~MyAnalysisdKnpipi()
{
  Clear();
  rtFile->cd();
  rtFile->Write();
  rtFile->Close();
}

void MyAnalysisdKnpipi::Clear()
{
}

bool MyAnalysisdKnpipi::DoAnalysis(ConfMan* conf, EventHeader* header, BeamLineHitMan* blMan, BeamLineTrackMan* bltrackMan, CDSHitMan* cdsMan, CDSTrackingMan* cdstrackMan, Particle* particle)
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
  if(particle->nPiplus()!=1) return true;
  if(particle->nPiminus()!=1) return true;
  FillHist("EventNumber",3);

  pCDS* pip = particle->pip(0);
  pCDS* pim = particle->pim(0);
  pCDS* pro = particle->product(0);

  // ############### //
  // Vertex decision //
  // ############### //
  TVector3 vertex;
  double vdis = 9999;
  bool FIDUCIAL = false;
  int ChargeOfPi =  0;
  if(pim->vdis()>pip->vdis()){
    vdis = pip->vdis();
    vertex = pip->vbeam();
    ChargeOfPi = 1;
    if(GeomTools::GetID(vertex)==CID_Fiducial) FIDUCIAL = true;
  }
  else {
    vdis = pim->vdis();
    vertex = pim->vbeam();
    ChargeOfPi = -1;
    if(GeomTools::GetID(vertex)==CID_Fiducial) FIDUCIAL = true;
  }

  if(!FIDUCIAL) return true;
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
  TLorentzVector Lpip = pip->GetLorentzVector();
  TLorentzVector cLpip = pip->GetLorentzVector(); cLpip.Boost(-boost);
  TLorentzVector Lpim = pim->GetLorentzVector();
  TLorentzVector cLpim = pim->GetLorentzVector(); cLpim.Boost(-boost);
  TLorentzVector Lpipi = Lpip + Lpim;
  TLorentzVector cLpipi = cLpip + cLpim;
  double pip_mass[2]  = { Lpip.M()                           , cLpip.M()};
  double pip_mom[2]   = { Lpip.Vect().Mag()                  , cLpip.Vect().Mag()};
  double pip_cost[2]  = { Lpip.Vect().CosTheta()             , cLpip.Vect().CosTheta()};
  double pip_phi[2]   = { Lpip.Vect().Phi()                  , cLpip.Vect().Phi()};
  double pip_mmass[2] = { (Ltgt+Lbeam-Lpip).M()              , (cLtgt+cLbeam-cLpip).M()};
  double pip_mmom[2]  = { (Ltgt+Lbeam-Lpip).Vect().Mag()     , (cLtgt+cLbeam-cLpip).Vect().Mag()};
  double pip_mcost[2] = { (Ltgt+Lbeam-Lpip).Vect().CosTheta(), (cLtgt+cLbeam-cLpip).Vect().CosTheta()};

  double pim_mass[2]  = { Lpim.M()                           , cLpim.M()};
  double pim_mom[2]   = { Lpim.Vect().Mag()                  , cLpim.Vect().Mag()};
  double pim_cost[2]  = { Lpim.Vect().CosTheta()             , cLpim.Vect().CosTheta()};
  double pim_phi[2]   = { Lpim.Vect().Phi()                  , cLpim.Vect().Phi()};
  double pim_mmass[2] = { (Ltgt+Lbeam-Lpim).M()              , (cLtgt+cLbeam-cLpim).M()};
  double pim_mmom[2]  = { (Ltgt+Lbeam-Lpim).Vect().Mag()     , (cLtgt+cLbeam-cLpim).Vect().Mag()};
  double pim_mcost[2] = { (Ltgt+Lbeam-Lpim).Vect().CosTheta(), (cLtgt+cLbeam-cLpim).Vect().CosTheta()};

  double pipi_mass[2]  = { Lpipi.M()                           , cLpipi.M()};
  double pipi_mom[2]   = { Lpipi.Vect().Mag()                  , cLpipi.Vect().Mag()};
  double pipi_cost[2]  = { Lpipi.Vect().CosTheta()             , cLpipi.Vect().CosTheta()};
  double pipi_phi[2]   = { Lpipi.Vect().Phi()                  , cLpipi.Vect().Phi()};
  double pipi_mmass[2] = { (Ltgt+Lbeam-Lpipi).M()              , (cLtgt+cLbeam-cLpipi).M()};
  double pipi_mmom[2]  = { (Ltgt+Lbeam-Lpipi).Vect().Mag()     , (cLtgt+cLbeam-cLpipi).Vect().Mag()};
  double pipi_mcost[2] = { (Ltgt+Lbeam-Lpipi).Vect().CosTheta(), (cLtgt+cLbeam-cLpipi).Vect().CosTheta()};

  double pip_vdis = pip->vdis();
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

  double npip_mass[2]  = { (Lpip+Ln).M()                           , (cLpip+cLn).M()};
  double npip_mom[2]   = { (Lpip+Ln).Vect().Mag()                  , (cLpip+cLn).Vect().Mag()};
  double npip_cost[2]  = { (Lpip+Ln).Vect().CosTheta()             , (cLpip+cLn).Vect().CosTheta()};
  double npip_phi[2]   = { (Lpip+Ln).Vect().Phi()                  , (cLpip+cLn).Vect().Phi()};
  double npip_mmass[2] = { (Ltgt+Lbeam-(Lpip+Ln)).M()              , (cLtgt+cLbeam-(cLpip+cLn)).M()};
  double npip_mmom[2]  = { (Ltgt+Lbeam-(Lpip+Ln)).Vect().Mag()     , (cLtgt+cLbeam-(cLpip+cLn)).Vect().Mag()};
  double npip_mcost[2] = { (Ltgt+Lbeam-(Lpip+Ln)).Vect().CosTheta(), (cLtgt+cLbeam-(cLpip+cLn)).Vect().CosTheta()};

  double npim_mass[2]  = { (Lpim+Ln).M()                           , (cLpim+cLn).M()};
  double npim_mom[2]   = { (Lpim+Ln).Vect().Mag()                  , (cLpim+cLn).Vect().Mag()};
  double npim_cost[2]  = { (Lpim+Ln).Vect().CosTheta()             , (cLpim+cLn).Vect().CosTheta()};
  double npim_phi[2]   = { (Lpim+Ln).Vect().Phi()                  , (cLpim+cLn).Vect().Phi()};
  double npim_mmass[2] = { (Ltgt+Lbeam-(Lpim+Ln)).M()              , (cLtgt+cLbeam-(cLpim+cLn)).M()};
  double npim_mmom[2]  = { (Ltgt+Lbeam-(Lpim+Ln)).Vect().Mag()     , (cLtgt+cLbeam-(cLpim+cLn)).Vect().Mag()};
  double npim_mcost[2] = { (Ltgt+Lbeam-(Lpim+Ln)).Vect().CosTheta(), (cLtgt+cLbeam-(cLpim+cLn)).Vect().CosTheta()};

  double npipi_mass[2] = { (Lpipi+Ln).M()                           , (Lpipi+cLn).M()};
  double npipi_mom[2]  = { (Lpipi+Ln).Vect().Mag()                  , (Lpipi+cLn).Vect().Mag()};
  double npipi_cost[2] = { (Lpipi+Ln).Vect().CosTheta()             , (Lpipi+cLn).Vect().CosTheta()};
  double npipi_phi[2]  = { (Lpipi+Ln).Vect().Phi()                  , (Lpipi+cLn).Vect().Phi()};
  double npipi_mmass[2] = { (Ltgt+Lbeam-(Lpipi+Ln)).M()              , (cLtgt+cLbeam-(cLpipi+cLn)).M()};
  double npipi_mmom[2]  = { (Ltgt+Lbeam-(Lpipi+Ln)).Vect().Mag()     , (cLtgt+cLbeam-(cLpipi+cLn)).Vect().Mag()};
  double npipi_mcost[2] = { (Ltgt+Lbeam-(Lpipi+Ln)).Vect().CosTheta(), (cLtgt+cLbeam-(cLpipi+cLn)).Vect().CosTheta()};

  double seg = (nc->seg()-1)%16+1;
  double lay = (nc->seg()-1)/16+1;

  /* flag */
  bool k0flag = false; /* k0 flag */
  bool sblk0flag = false; /* lower sideband of k0 flag */
  bool sbuk0flag = false; /* upper sideband of k0 flag */
  {
    double pro_mass = pipi_mass[1];
    if(k0ll<pro_mass&&pro_mass<k0ul) k0flag = true;
    if(sblk0ll<pro_mass&&pro_mass<sblk0ul) sblk0flag = true;
    if(sbuk0ll<pro_mass&&pro_mass<sbuk0ul) sbuk0flag = true;
  }
  bool spflag = false; /* sp flag */
  bool sblspflag = false; /* lower sideband of sp flag */
  bool sbuspflag = false; /* upper sideband of sp flag */
  {
    double pro_mass = npip_mass[1];
    if(spll<pro_mass&&pro_mass<spul) spflag = true;
    if(sblspll<pro_mass&&pro_mass<sblspul) sblspflag = true;
    if(sbuspll<pro_mass&&pro_mass<sbuspul) sbuspflag = true;
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
  bool mnflag = false; /* mm flag */
  bool sblmnflag = false; /* lower sideband of mm flag */
  bool sbumnflag = false; /* upper sideband of mm flag */
  {
    double pro_mmass = npipi_mmass[1];
    if(mnll<pro_mmass&&pro_mmass<mnul) mnflag = true;
    if(sblmnll<pro_mmass&&pro_mmass<sblmnul) sblmnflag = true;
    if(sbumnll<pro_mmass&&pro_mmass<sbumnul) sbumnflag = true;
  }

  //if(!k0flag){ return true; }
  //if(k0flag){ return true; }
  //if(!smflag){ return true; }
  //if(smflag){ return true; }
  //if(!spflag){ return true; }
  //if(spflag){ return true; }
  //if(!spflag&&!smflag){ return true; }
  //if(spflag||smflag){ return true; }
  //if(!mnflag){ return true; }
  //if(mnflag){ return true; }

  FillHist("FWD_n_OverbetavsMomentum",1.0/nc->beta(),nc->mom().Mag());
  FillHist("FWD_n_OverbetavsEnergy",1.0/nc->beta(),nc->energy());
  FillHist("FWD_n_TOFvsMomentum",nc->tof(),nc->mom().Mag());
  FillHist("FWD_n_HitPosition",nc->hitpos().X(),nc->hitpos().Y());
  FillHist("FWD_n_HitSegment",seg,lay);

  /*-----------------------------------------*/
  /* All                                     */ 
  /*-----------------------------------------*/
  // pi+
  for(int f=1; f<2; f++){
    FillHist(Form("pip_CosTvsMomentum_%s",frame[f].Data()),pip_cost[f],pip_mom[f]);
    FillHist(Form("pip_CosTvsPhi_%s",frame[f].Data()),pip_cost[f],pip_phi[f]);
    FillHist(Form("dKpip_MMvsMomentum_%s",frame[f].Data()),pip_mmass[f],pip_mmom[f]);
    FillHist(Form("dKpip_MMvsCosT_%s",frame[f].Data()),pip_mmass[f],pip_mcost[f]);
    FillHist(Form("pip_IMvsMomentum_%s",frame[f].Data()),pip_mass[f],pip_mom[f]);
    FillHist(Form("pip_IMvsCosT_%s",frame[f].Data()),pip_mass[f],pip_cost[f]);
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
  // pi+ pi-
  for(int f=1; f<2; f++){
    FillHist(Form("pipi_CosTvsMomentum_%s",frame[f].Data()),pipi_cost[f],pipi_mom[f]);
    FillHist(Form("pipi_CosTvsPhi_%s",frame[f].Data()),pipi_cost[f],pipi_phi[f]);
    FillHist(Form("dKpipi_MMvsMomentum_%s",frame[f].Data()),pipi_mmass[f],pipi_mmom[f]);
    FillHist(Form("dKpipi_MMvsCosT_%s",frame[f].Data()),pipi_mmass[f],pipi_mcost[f]);
    FillHist(Form("pipi_IMvsMomentum_%s",frame[f].Data()),pipi_mass[f],pipi_mom[f]);
    FillHist(Form("pipi_IMvsCosT_%s",frame[f].Data()),pipi_mass[f],pipi_cost[f]);
  }
  // neutron pi+
  for(int f=1; f<2; f++){
    FillHist(Form("npip_CosTvsMomentum_%s",frame[f].Data()),npip_cost[f],npip_mom[f]);
    FillHist(Form("npip_CosTvsPhi_%s",frame[f].Data()),npip_cost[f],npip_phi[f]);
    FillHist(Form("dKnpip_MMvsMomentum_%s",frame[f].Data()),npip_mmass[f],npip_mmom[f]);
    FillHist(Form("dKnpip_MMvsCosT_%s",frame[f].Data()),npip_mmass[f],npip_mcost[f]);
    FillHist(Form("npip_IMvsMomentum_%s",frame[f].Data()),npip_mass[f],npip_mom[f]);
    FillHist(Form("npip_IMvsCosT_%s",frame[f].Data()),npip_mass[f],npip_cost[f]);
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
  // neutron pi+ pi-
  for(int f=1; f<2; f++){
    FillHist(Form("npipi_CosTvsMomentum_%s",frame[f].Data()),npipi_cost[f],npipi_mom[f]);
    FillHist(Form("npipi_CosTvsPhi_%s",frame[f].Data()),npipi_cost[f],npipi_phi[f]);
    FillHist(Form("dKnpipi_MMvsMomentum_%s",frame[f].Data()),npipi_mmass[f],npipi_mmom[f]);
    FillHist(Form("dKnpipi_MMvsCosT_%s",frame[f].Data()),npipi_mmass[f],npipi_mcost[f]);
    FillHist(Form("npipi_IMvsMomentum_%s",frame[f].Data()),npipi_mass[f],npipi_mom[f]);
    FillHist(Form("npipi_IMvsCosT_%s",frame[f].Data()),npipi_mass[f],npipi_cost[f]);
  }

  //FillHist(Form("%s_DCA",pname[cds->pid()].Data()),cds->vdis());
  //FillHist(Form("all_DCA"),cds->vdis());

  return true;

}

bool MyAnalysisdKnpipi::FillHist(TString name, double val1, int weight)
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

bool MyAnalysisdKnpipi::FillHist(TString name, TString val1, int weight)
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

bool MyAnalysisdKnpipi::FillHist(TString name, double val1, double val2, int weight)
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

bool MyAnalysisdKnpipi::FillHist(TString name, TString val1, TString val2, int weight)
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

void MyAnalysisdKnpipi::CutCondition()
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
}

bool MyAnalysisdKnpipi::Initialize(ConfMan* confMan)
{
  std::cout << "### MyAnalysisdKnpipi::Initialize ###" << std::endl;

  std::string ofname = confMan->GetOutFileName();
  ofname.insert(ofname.find(".root"),"_anadKnpipi");

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
  // pi+
  new TH2F("pip_CosTvsMomentum_CM","cos#theta vs. momentum of #pi^{+};cos#theta;Momentum (GeV/c)", 2000, -1.0, 1.0, 2000, 0.0, 2.0);
  new TH2F("pip_CosTvsPhi_CM","cos#theta vs. #phi of #pi^{+};cos#theta;#phi", 2000, -1.0, 1.0, 1600, 0.0, 3.2);
  new TH2F("dKpip_MMvsMomentum_CM", "d(K^{-},#pi^{+})X missing mass vs. missing momentum;MM_{d}(#pi^{+}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 1.0, 3.0, 2000, 0.0, 2.0);
  new TH2F("dKpip_MMvsCosT_CM", "d(K^{-},#pi^{+})X missing mass vs. cos#theta;MM_{d}(#pi^{+}) (GeV/c^{2});cos#theta", 2000, 1.0, 3.0, 2000, -1.0, 1.0);
  new TH2F("pip_IMvsMomentum_CM", "(#pi^{+}) invariant mass vs. momentum;IM(#pi^{+}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("pip_IMvsCosT_CM", "(#pi^{+}) invariant mass vs. cos#theta;IM(#pi^{+}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
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
  // pi+ pi-
  new TH2F("pipi_CosTvsMomentum_CM","cos#theta vs. momentum of #pi^{+}#pi^{-};cos#theta;Momentum (GeV/c)", 2000, -1.0, 1.0, 2000, 0.0, 2.0);
  new TH2F("pipi_CosTvsPhi_CM","cos#theta vs. #phi of #pi^{+}#pi^{-};cos#theta;#phi", 2000, -1.0, 1.0, 1600, 0.0, 3.2);
  new TH2F("dKpipi_MMvsMomentum_CM", "d(K^{-},#pi^{+}#pi^{-})X missing mass vs. missing momentum;MM_{d}(#pi^{+}#pi^{-}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 1.0, 3.0, 2000, 0.0, 2.0);
  new TH2F("dKpipi_MMvsCosT_CM", "d(K^{-},#pi^{+}#pi^{-})X missing mass vs. cos#theta;MM_{d}(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 1.0, 3.0, 2000, -1.0, 1.0);
  new TH2F("pipi_IMvsMomentum_CM", "(#pi^{+}#pi^{-}) invariant mass vs. momentum;IM(#pi^{+}#pi^{-}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("pipi_IMvsCosT_CM", "(#pi^{+}#pi^{-}) invariant mass vs. cos#theta;IM(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  // neutron pi+
  new TH2F("npip_CosTvsMomentum_CM","cos#theta vs. momentum of n#pi^{+};cos#theta;Momentum (GeV/c)", 2000, -1.0, 1.0, 2000, 0.0, 2.0);
  new TH2F("npip_CosTvsPhi_CM","cos#theta vs. #phi of n#pi^{+};cos#theta;#phi", 2000, -1.0, 1.0, 1600, 0.0, 3.2);
  new TH2F("dKnpip_MMvsMomentum_CM", "d(K^{-},n#pi^{+})X missing mass vs. missing momentum;MM_{d}(n#pi^{+}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("dKnpip_MMvsCosT_CM", "d(K^{-},n#pi^{+})X missing mass vs. cos#theta;MM_{d}(n#pi^{+}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("npip_IMvsMomentum_CM", "(n#pi^{+}) invariant mass vs. momentum;IM(n#pi^{+}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("npip_IMvsCosT_CM", "(n#pi^{+}) invariant mass vs. cos#theta;IM(n#pi^{+}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  // neutron pi-
  new TH2F("npim_CosTvsMomentum_CM","cos#theta vs. momentum of n#pi^{-};cos#theta;Momentum (GeV/c)", 2000, -1.0, 1.0, 2000, 0.0, 2.0);
  new TH2F("npim_CosTvsPhi_CM","cos#theta vs. #phi of n#pi^{-};cos#theta;#phi", 2000, -1.0, 1.0, 1600, 0.0, 3.2);
  new TH2F("dKnpim_MMvsMomentum_CM", "d(K^{-},n#pi^{-})X missing mass vs. missing momentum;MM_{d}(n#pi^{-}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("dKnpim_MMvsCosT_CM", "d(K^{-},n#pi^{-})X missing mass vs. cos#theta;MM_{d}(n#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("npim_IMvsMomentum_CM", "(n#pi^{-}) invariant mass vs. momentum;IM(n#pi^{-}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("npim_IMvsCosT_CM", "(n#pi^{-}) invariant mass vs. cos#theta;IM(n#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  // neutron pi+ pi-
  new TH2F("npipi_CosTvsMomentum_CM","cos#theta vs. momentum of n#pi^{+}#pi^{-};cos#theta;Momentum (GeV/c)", 2000, -1.0, 1.0, 2000, 0.0, 2.0);
  new TH2F("npipi_CosTvsPhi_CM","cos#theta vs. #phi of n#pi^{+}#pi^{-};cos#theta;#phi", 2000, -1.0, 1.0, 1600, 0.0, 3.2);
  new TH2F("dKnpipi_MMvsMomentum_CM", "d(K^{-},n#pi^{+}#pi^{-})X missing mass vs. missing momentum;MM_{d}(n#pi^{+}#pi^{-}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("dKnpipi_MMvsCosT_CM", "d(K^{-},n#pi^{+}#pi^{-})X missing mass vs. cos#theta;MM_{d}(n#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("npipi_IMvsMomentum_CM", "(n#pi^{+}#pi^{-}) invariant mass vs. momentum;IM(n#pi^{+}#pi^{-}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("npipi_IMvsCosT_CM", "(n#pi^{+}#pi^{-}) invariant mass vs. cos#theta;IM(n#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);


  return true;

}
