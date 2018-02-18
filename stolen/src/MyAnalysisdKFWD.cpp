// MyAnalysisdKFWD.cpp

#include "MyAnalysisdKFWD.h"

MyAnalysisdKFWD::MyAnalysisdKFWD(TFile* rt, ConfMan* conf)
{
  Initialize(conf);
  CutCondition();
  Clear();
}

MyAnalysisdKFWD::~MyAnalysisdKFWD()
{
  Clear();
  rtFile->cd();
  rtFile->Write();
  rtFile->Close();
}

void MyAnalysisdKFWD::Clear()
{
}

bool MyAnalysisdKFWD::DoAnalysis(ConfMan* conf, EventHeader* header, BeamLineHitMan* blMan, BeamLineTrackMan* bltrackMan, CDSHitMan* cdsMan, CDSTrackingMan* cdstrackMan, Particle* particle)
{
  rtFile->cd();

  FillHist("EventNumber",0); /* All Events */
  DetectorList *dlist=DetectorList::GetInstance();
  if(particle->nBeam()!=1) return false;
  pBeam* beam = particle->beam(0);
  double beamtof = beam->bhdt0tof();

  /* Event selection */
  if(particle->nCDS()==0) return true;
  FillHist("EventNumber",1); /* Least 1 hit in CDS*/
  bool neutral=false; bool charged=false;
  if(header->IsTrig(Trig_Neutral)&&particle->nNC()>0) neutral = true;
  if(header->IsTrig(Trig_Charged)&&particle->nPC()>0) charged = true;
  if(neutral&&charged) return true;
  FillHist("EventNumber",2); /* Both Neutral and Charged */

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
  if(fcds==-1) return true; /* Vertex decision */
  pCDS* cds=0;
  if(fcds!=-1) cds = particle->cdsi(fcds);

  if(neutral) DoNeutralAna(particle,beam,cds,vertex);
  if(charged) DoChargedAna(particle,beam,cds,vertex);

  return true;

}

bool MyAnalysisdKFWD::DoNeutralAna(Particle* particle, pBeam* beam, pCDS* cds, TVector3 vertex)
{
  DetectorList *dlist=DetectorList::GetInstance();

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
  TLorentzVector Ltgt;
  if(dlist->GetMaterial(CID_Target)=="LHydrogen"){
    Ltgt.SetVectM(ZeroV,pMass);
  }
  if(dlist->GetMaterial(CID_Target)=="LDeuterium"){
    Ltgt.SetVectM(ZeroV,dMass);
  }
  if(dlist->GetMaterial(CID_Target)=="LHelium-3"){
    Ltgt.SetVectM(ZeroV,ThreeHeMass);
  }
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
    FillHist(Form("n_Mass_%s",frame[f].Data()),nc_mass[f]);
    FillHist(Form("n_Momentum_%s",frame[f].Data()),nc_mom[f]);
    FillHist(Form("n_CosT_%s",frame[f].Data()),nc_cost[f]);
    FillHist(Form("n_Phi_%s",frame[f].Data()),nc_phi[f]);
    FillHist(Form("n_MMass_%s",frame[f].Data()),nc_mmass[f]);
    FillHist(Form("n_MMomentum_%s",frame[f].Data()),nc_mmom[f]);
    FillHist(Form("n_MCosT_%s",frame[f].Data()),nc_mcost[f]);
  }

  // n + pi+
  if(cds->pid()==CDS_PiPlus){
    for(int f=1; f<2; f++){
      FillHist(Form("npip_Mass_%s",frame[f].Data()),ncds_mass[f]);
      FillHist(Form("npip_Momentum_%s",frame[f].Data()),ncds_mom[f]);
      FillHist(Form("npip_CosT_%s",frame[f].Data()),ncds_cost[f]);
      FillHist(Form("npip_Phi_%s",frame[f].Data()),ncds_phi[f]);
      FillHist(Form("npip_MMass_%s",frame[f].Data()),ncds_mmass[f]);
      FillHist(Form("npip_MMomentum_%s",frame[f].Data()),ncds_mmom[f]);
      FillHist(Form("npip_MCosT_%s",frame[f].Data()),ncds_mcost[f]);
    }
  }
  // n + pi-
  if(cds->pid()==CDS_PiMinus){
    for(int f=1; f<2; f++){
      FillHist(Form("npim_Mass_%s",frame[f].Data()),ncds_mass[f]);
      FillHist(Form("npim_Momentum_%s",frame[f].Data()),ncds_mom[f]);
      FillHist(Form("npim_CosT_%s",frame[f].Data()),ncds_cost[f]);
      FillHist(Form("npim_Phi_%s",frame[f].Data()),ncds_phi[f]);
      FillHist(Form("npim_MMass_%s",frame[f].Data()),ncds_mmass[f]);
      FillHist(Form("npim_MMomentum_%s",frame[f].Data()),ncds_mmom[f]);
      FillHist(Form("npim_MCosT_%s",frame[f].Data()),ncds_mcost[f]);
    }
  }
  // n + k-
  if(cds->pid()==CDS_Kaon){
    for(int f=1; f<2; f++){
      FillHist(Form("nk_Mass_%s",frame[f].Data()),ncds_mass[f]);
      FillHist(Form("nk_Momentum_%s",frame[f].Data()),ncds_mom[f]);
      FillHist(Form("nk_CosT_%s",frame[f].Data()),ncds_cost[f]);
      FillHist(Form("nk_Phi_%s",frame[f].Data()),ncds_phi[f]);
      FillHist(Form("nk_MMass_%s",frame[f].Data()),ncds_mmass[f]);
      FillHist(Form("nk_MMomentum_%s",frame[f].Data()),ncds_mmom[f]);
      FillHist(Form("nk_MCosT_%s",frame[f].Data()),ncds_mcost[f]);
    }
  }

  return true;

}

bool MyAnalysisdKFWD::DoChargedAna(Particle* particle, pBeam* beam, pCDS* cds, TVector3 vertex)
{
  DetectorList *dlist=DetectorList::GetInstance();

  // ##################### //
  // FWD Charged  Analysis //
  // ##################### //
  for(int ipc=0; ipc<particle->nPC(); ipc++){
    pPC* pc = particle->pc(ipc);
    int pid = -1;
    pc->CalcMom(beam,vertex,pid,false);
    double r = pc->r();
    double mass2 = pc->mass2();
    double momr = pc->momr();
    double angle = pc->angle();
    double mom = pc->mom().Mag();
    double tof = pc->tof();
    double fl = pc->fl();
    double energy = pc->energy();
    double beta = pc->beta();
    
    FillHist("FWDC_Beta",beta);
    FillHist("FWDC_Overbeta",1.0/beta);
    FillHist("FWDC_OverbetavsEnergy",1.0/beta,energy);
    FillHist("FWDC_MomentumR",momr);
    FillHist("FWDC_Momentum",mom);
    FillHist("FWDC_MomentumvsMomentumR",mom,momr);
    FillHist("FWDC_Mass2",mass2);
    FillHist("FWDC_Mass2vsMomentum",mass2,momr);
  }

  if(particle->nPC()!=1) return true;
  pPC* pc = 0;
  pc = particle->pc(0);
  if(pc->pid()!=F_Proton) return true;
  
  TString frame[2] = {"Lab","CM"};

  // Target
  TVector3 ZeroV;
  TLorentzVector Ltgt;
  if(dlist->GetMaterial(CID_Target)=="LHydrogen"){
    Ltgt.SetVectM(ZeroV,pMass);
  }
  if(dlist->GetMaterial(CID_Target)=="LDeuterium"){
    Ltgt.SetVectM(ZeroV,dMass);
  }
  if(dlist->GetMaterial(CID_Target)=="LHelium-3"){
    Ltgt.SetVectM(ZeroV,ThreeHeMass);
  }
  TLorentzVector Lbeam = beam->GetLorentzVector(vertex);
  TVector3 boost = (Ltgt+Lbeam).BoostVector();
  TLorentzVector cLtgt = Ltgt; cLtgt.Boost(-boost); TLorentzVector cLbeam = Lbeam; cLbeam.Boost(-boost);
  /* CDS */
  TLorentzVector Lcds = cds->GetLorentzVector();
  TLorentzVector cLcds = cds->GetLorentzVector(); cLcds.Boost(-boost);
  /* PC */
  TLorentzVector Lp  = pc ->GetLorentzVector();
  TLorentzVector cLp  = pc ->GetLorentzVector(); cLp.Boost(-boost);
  double pc_mass[2]    = { Lp.M()                           , cLp.M()};
  double pc_mom[2]     = { Lp.Vect().Mag()                  , cLp.Vect().Mag()};
  double pc_cost[2]    = { Lp.Vect().CosTheta()             , cLp.Vect().CosTheta()};
  double pc_phi[2]     = { Lp.Vect().Phi()                  , cLp.Vect().Phi()};
  double pc_mmass[2]   = { (Ltgt+Lbeam-Lp).M()              , (cLtgt+cLbeam-cLp).M()};
  double pc_mmom[2]    = { (Ltgt+Lbeam-Lp).Vect().Mag()     , (cLtgt+cLbeam-cLp).Vect().Mag()};
  double pc_mcost[2]   = { (Ltgt+Lbeam-Lp).Vect().CosTheta(), (cLtgt+cLbeam-cLp).Vect().CosTheta()};
  /* CDS + PC */
  double pcds_mass[2]  = { (Lcds+Lp).M()                           , (cLcds+cLp).M()};
  double pcds_mom[2]   = { (Lcds+Lp).Vect().Mag()                  , (cLcds+cLp).Vect().Mag()};
  double pcds_cost[2]  = { (Lcds+Lp).Vect().CosTheta()             , (cLcds+cLp).Vect().CosTheta()};
  double pcds_phi[2]   = { (Lcds+Lp).Vect().Phi()                  , (cLcds+cLp).Vect().Phi()};
  double pcds_mmass[2] = { (Ltgt+Lbeam-(Lcds+Lp)).M()              , (cLtgt+cLbeam-(cLcds+cLp)).M()};
  double pcds_mmom[2]  = { (Ltgt+Lbeam-(Lcds+Lp)).Vect().Mag()     , (cLtgt+cLbeam-(cLcds+cLp)).Vect().Mag()};
  double pcds_mcost[2] = { (Ltgt+Lbeam-(Lcds+Lp)).Vect().CosTheta(), (cLtgt+cLbeam-(cLcds+cLp)).Vect().CosTheta()};

  // p
  for(int f=1; f<2; f++){
    FillHist(Form("pfwd_Mass_%s",frame[f].Data()),pc_mass[f]);
    FillHist(Form("pfwd_Momentum_%s",frame[f].Data()),pc_mom[f]);
    FillHist(Form("pfwd_CosT_%s",frame[f].Data()),pc_cost[f]);
    FillHist(Form("pfwd_Phi_%s",frame[f].Data()),pc_phi[f]);
    FillHist(Form("pfwd_MMass_%s",frame[f].Data()),pc_mmass[f]);
    FillHist(Form("pfwd_MMomentum_%s",frame[f].Data()),pc_mmom[f]);
    FillHist(Form("pfwd_MCosT_%s",frame[f].Data()),pc_mcost[f]);
  }

  // n + pi-
  if(cds->pid()==CDS_PiMinus){
    for(int f=1; f<2; f++){
      FillHist(Form("pfwdpim_Mass_%s",frame[f].Data()),pcds_mass[f]);
      FillHist(Form("pfwdpim_Momentum_%s",frame[f].Data()),pcds_mom[f]);
      FillHist(Form("pfwdpim_CosT_%s",frame[f].Data()),pcds_cost[f]);
      FillHist(Form("pfwdpim_Phi_%s",frame[f].Data()),pcds_phi[f]);
      FillHist(Form("pfwdpim_MMass_%s",frame[f].Data()),pcds_mmass[f]);
      FillHist(Form("pfwdpim_MMomentum_%s",frame[f].Data()),pcds_mmom[f]);
      FillHist(Form("pfwdpim_MCosT_%s",frame[f].Data()),pcds_mcost[f]);
    }
  }
  // n + k-
  if(cds->pid()==CDS_Kaon){
    for(int f=1; f<2; f++){
      FillHist(Form("pfwdk_Mass_%s",frame[f].Data()),pcds_mass[f]);
      FillHist(Form("pfwdk_Momentum_%s",frame[f].Data()),pcds_mom[f]);
      FillHist(Form("pfwdk_CosT_%s",frame[f].Data()),pcds_cost[f]);
      FillHist(Form("pfwdk_Phi_%s",frame[f].Data()),pcds_phi[f]);
      FillHist(Form("pfwdk_MMass_%s",frame[f].Data()),pcds_mmass[f]);
      FillHist(Form("pfwdk_MMomentum_%s",frame[f].Data()),pcds_mmom[f]);
      FillHist(Form("pfwdk_MCosT_%s",frame[f].Data()),pcds_mcost[f]);
    }
  }
  return true;

}

bool MyAnalysisdKFWD::FillHist(TString name, double val1, int weight)
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

bool MyAnalysisdKFWD::FillHist(TString name, TString val1, int weight)
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

bool MyAnalysisdKFWD::FillHist(TString name, double val1, double val2, int weight)
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

bool MyAnalysisdKFWD::FillHist(TString name, TString val1, TString val2, int weight)
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

void MyAnalysisdKFWD::CutCondition()
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

bool MyAnalysisdKFWD::Initialize(ConfMan* confMan)
{
  std::cout << "### MyAnalysisdKFWD::Initialize ###" << std::endl;

  std::string ofname = confMan->GetOutFileName();
  ofname.insert(ofname.find(".root"),"_anadKFWD");

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
  // Missing mass p(K-,n)X
  std::cout << "Define Histograms for (K-,n)X M.M." << std::endl;
  new TH1F("n_Mass_CM","Invariant mass of n;IM(n) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("n_Momentum_CM","Momentum of n;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("n_CosT_CM","cos#theta of n;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("n_Phi_CM","#phi of n;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("n_MMass_CM", "^{3}He(K^{-},n)X missing mass;MM(n) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("n_MMomentum_CM", "^{3}He(K^{-},n)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("n_MCosT_CM", "^{3}He(K^{-},n)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // Invariant mass with FWD
  std::cout << "Define Histograms for product I.M. with FWD" << std::endl;
  // n + pi+
  new TH1F("npip_Mass_CM","Invariant mass of n#pi^{+};IM(n#pi^{+}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("npip_Momentum_CM","Momentum of n#pi^{+};Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("npip_CosT_CM","cos#theta of n#pi^{+};cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("npip_Phi_CM","#phi of n#pi^{+};#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("npip_MMass_CM", "^{3}He(K^{-},n#pi^{+})X missing mass;MM(n#pi^{+}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("npip_MMomentum_CM", "^{3}He(K^{-},n#pi^{+})X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("npip_MCosT_CM", "^{3}He(K^{-},n#pi^{+})X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // n + pi-
  new TH1F("npim_Mass_CM","Invariant mass of n#pi^{-};IM(n#pi^{-}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("npim_Momentum_CM","Momentum of n#pi^{-};Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("npim_CosT_CM","cos#theta of n#pi^{-};cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("npim_Phi_CM","#phi of n#pi^{-};#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("npim_MMass_CM", "^{3}He(K^{-},n#pi^{-})X missing mass;MM(n#pi^{-}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("npim_MMomentum_CM", "^{3}He(K^{-},n#pi^{-})X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("npim_MCosT_CM", "^{3}He(K^{-},n#pi^{-})X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // n + k-
  new TH1F("nk_Mass_CM","Invariant mass of nK^{-};IM(nK^{-}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("nk_Momentum_CM","Momentum of nK^{-};Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("nk_CosT_CM","cos#theta of nK^{-};cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("nk_Phi_CM","#phi of nK^{-};#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("nk_MMass_CM", "^{3}He(K^{-},nK^{-})X missing mass;MM(nK^{-}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("nk_MMomentum_CM", "^{3}He(K^{-},nK^{-})X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("nk_MCosT_CM", "^{3}He(K^{-},nK^{-})X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);

  std::cout << "=== End of [FWDNeutral::Initialize] === " << std::endl;


  std::cout << "Define Histograms for FWD Charged particles" << std::endl;
  // FDC Tracking
  new TH1F( "FDC1_nTrack", "Number of Tracks in FDC1;Num. of Tracks;Counts", 20, 0, 20 );
  new TH1F( "FDC1_Chi2", "Chi-squared of FDC1;Chi-square;Counts", 500, 0, 50);
  new TH2F( "FDC1_Position", "Position at FDC1;X (cm);Y (cm)", 1000,-50,50, 1000,-50,50);
  // FWD Charged
  new TH1F( "FWDC_Beta", "#beta of FWD Charged;#beta;Counts", 1000,0.0,1.0);
  new TH1F( "FWDC_Overbeta", "1/#beta of FWD Charged;1/#beta;Counts", 5000,0.0,5.0);
  new TH2F( "FWDC_OverbetavsEnergy", "1/#beta vs. Energy deposit;1/#beta;Energy deposit (MeVee)", 1000, -0, 5, 500, 0, 50);
  new TH1F( "FWDC_Momentum", "Momentum of FWD Charged;Momentum (GeV/c);Counts", 2000,0.0,2.0);
  new TH1F( "FWDC_MomentumR", "Momentum of FWD Charged;Momentum (GeV/c);Counts", 2000,0.0,2.0);
  new TH2F( "FWDC_MomentumvsMomentumR", "Momentum vs Mometuntum_{R} of FWD Charged;Momentum (GeV/c);Momentum_{R} (GeV/c)", 2000,0.0,2.0, 2000,0.0,2.0);
  new TH1F( "FWDC_Mass2", "Mass2 of FWD Charged;Mass^{2} (GeV^{2}/c^{4});Counts", 600,-1.0,5.0);
  new TH2F( "FWDC_Mass2vsMomentum", "Mass2 vs Momentum of FWD Charged;Mass^{2} (GeV^{2}/c^{4});Momentum", 600,-1.0,5.0, 200, 0.0, 2.0);

  // p
  new TH1F("pfwd_Mass_CM","Invariant mass of n;IM(p) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pfwd_Momentum_CM","Momentum of p;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pfwd_CosT_CM","cos#theta of p;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("pfwd_Phi_CM","#phi of p;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("pfwd_MMass_CM", "^{3}He(K^{-},p)X missing mass;MM(p) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pfwd_MMomentum_CM", "^{3}He(K^{-},p)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pfwd_MCosT_CM", "^{3}He(K^{-},p)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // p + pi-
  new TH1F("pfwdpim_Mass_CM","Invariant mass of p#pi^{-};IM(p#pi^{-}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pfwdpim_Momentum_CM","Momentum of p#pi^{-};Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pfwdpim_CosT_CM","cos#theta of p#pi^{-};cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("pfwdpim_Phi_CM","#phi of p#pi^{-};#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("pfwdpim_MMass_CM", "^{3}He(K^{-},p#pi^{-})X missing mass;MM(p#pi^{-}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pfwdpim_MMomentum_CM", "^{3}He(K^{-},p#pi^{-})X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pfwdpim_MCosT_CM", "^{3}He(K^{-},p#pi^{-})X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // p + k-
  std::cout << "Define Histograms for product I.M. with FWD" << std::endl;
  new TH1F("pfwdk_Mass_CM","Invariant mass of pK^{-};IM(pK^{-}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pfwdk_Momentum_CM","Momentum of pK^{-};Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pfwdk_CosT_CM","cos#theta of pK^{-};cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("pfwdk_Phi_CM","#phi of pK^{-};#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("pfwdk_MMass_CM", "^{3}He(K^{-},pK^{-})X missing mass;MM(pK^{-}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pfwdk_MMomentum_CM", "^{3}He(K^{-},pK^{-})X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pfwdk_MCosT_CM", "^{3}He(K^{-},pK^{-})X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);

  std::cout << "=== End of [FWDCharged::Initialize] === " << std::endl;

  return true;
}
