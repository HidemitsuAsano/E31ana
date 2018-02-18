// MyAnalysisFWDCharged.cpp

#include "MyAnalysisFWDCharged.h"

MyAnalysisFWDCharged::MyAnalysisFWDCharged(TFile* rtFile, ConfMan* conf, bool simflag)
{
  Initialize(rtFile, conf);
  SIMULATION = simflag;
  Clear();
}

MyAnalysisFWDCharged::~MyAnalysisFWDCharged()
{
}

void MyAnalysisFWDCharged::Clear()
{
}

bool MyAnalysisFWDCharged::DoAnalysis(ConfMan* conf, EventHeader* header, BeamLineHitMan* blMan, BeamLineTrackMan* bltrackMan, CDSHitMan* cdsMan, CDSTrackingMan* cdstrackMan, Particle* particle)
{
  DetectorList *dlist=DetectorList::GetInstance();

  if(particle->nBeam()!=1) return false;
  pBeam* beam = particle->beam(0);

  // ################# //
  // Charged selection //
  // ################# //
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
  if(MulCVC==0&&MulPC==0) return true;

  // ############# //
  // FDC Selection //
  // ############# //
  FillHist("FDC1_nTrack",bltrackMan->ntrackBLDC(CID_FDC1));
  if(bltrackMan->ntrackBLDC(CID_FDC1)!=1) return true;

  LocalTrack* fdc = bltrackMan->trackBLDC(CID_FDC1,0);
  FillHist("FDC1_Chi2",fdc->chi2all());
  TVector3 fdcpos = fdc->GetPosatZ(169.6);
  FillHist("FDC1_Position",fdcpos.X(),fdcpos.Y());

  // ############ //
  // CVC analysis //
  // ############ //
  for(int i=0; i<blMan->nCVC(); i++){
    HodoscopeLikeHit* hit = blMan->CVC(i);
    if(hit->CheckRange()){
      pPC* pc = new pPC();
      pc->SetSegment(hit->seg());
      pc->SetHitPosition(hit->pos().X(),hit->hitpos(),hit->pos().Z());
      pc->SetTime(hit->ctmean());
      pc->SetEnergy(hit->emean());
      pc->SetFDC1Pos(fdcpos);
      particle->AddPC(*pc);
      delete pc;
    }
  }

  // ########### //
  // PC analysis //
  // ########### //
  for(int i=0; i<blMan->nPC(); i++){
    HodoscopeLikeHit* hit = blMan->PC(i);
    if(hit->CheckRange()){
      pPC* pc = new pPC();
      pc->SetSegment(hit->seg());
      pc->SetHitPosition(hit->pos().X(),hit->hitpos(),hit->pos().Z());
      pc->SetTime(hit->ctmean());
      pc->SetEnergy(hit->emean());
      pc->SetFDC1Pos(fdcpos);
      particle->AddPC(*pc);
      delete pc;
    }
  }

  // ############### //
  // Vertex decision //
  // ############### //
  TVector3 vertex;
  double vdis = 9999;
  int fcds = -1;
  for(int it=0; it<particle->nCDS(); it++){
    pCDS* cds = particle->cdsi(it);
    if(vdis>9998||vdis>cds->vdis()){
      vdis = cds->vdis();
      vertex = cds->vbeam();
      if(GeomTools::GetID(vertex)==CID_Fiducial) fcds=it;
    }
  }
  if(fcds==-1) return true; /* there is no vertex in fiducial volume */
  pCDS* cds=0;
  if(fcds!=-1) cds = particle->cdsi(fcds);

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
  if(dlist->GetMaterial(CID_Target)=="Vacuum"){
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

  // p + pi-
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
  // p + k-
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

bool MyAnalysisFWDCharged::FillHist(TString name, double val1)
{
  TH1F* h1 = (TH1F*)gFile -> Get(name);
  if(h1){
    h1 -> Fill(val1);
    return true;
  }
  else {
    return false;
  }
}

bool MyAnalysisFWDCharged::FillHist(TString name, double val1, double val2)
{
  TH2F* h2 = (TH2F*)gFile -> Get(name);
  if(h2){ h2 -> Fill(val1,val2);
    return true;
  }
  else {
    return false;
  }
}

bool MyAnalysisFWDCharged::Initialize(TFile* rtFile, ConfMan* confMan)
{
  rtFile -> cd();

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
  new TH2F( "FWDC_MomentumvsMomentumR", "Momentum vs Mometuntum_{R} of FWD Charged;Momentum (GeV/c);Momentum_{R} (GeV/c)", 200,0.0,2.0, 200,0.0,2.0);
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

  std::cout << "=== End of [MyAnalysisFWDCharged::Initialize] === " << std::endl;

  return true;
}
