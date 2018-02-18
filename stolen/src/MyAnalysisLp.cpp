// MyAnalysisLp.cpp

#include "MyAnalysisLp.h"

MyAnalysisLp::MyAnalysisLp(TFile* rt, ConfMan* conf)
{
  Initialize(conf);
  CutCondition();
  Clear();
}

MyAnalysisLp::~MyAnalysisLp()
{
  Clear();
  rtFile->cd();
  rtFile->Write();
  rtFile->Close();
}

void MyAnalysisLp::Clear()
{
}

bool MyAnalysisLp::DoAnalysis(ConfMan* conf, EventHeader* header, BeamLineHitMan* blMan, BeamLineTrackMan* bltrackMan, CDSHitMan* cdsMan, CDSTrackingMan* cdstrackMan, Particle* particle)
{
  rtFile->cd();

  FillHist("EventNumber",0);
  DetectorList *dlist=DetectorList::GetInstance();
  if(particle->nBeam()!=1) return false;
  pBeam* beam = particle->beam(0);
  FillHist("EventNumber",1);

  /* Event selection */
  if(particle->nCDS()!=3) return true;
  if(particle->nProton()!=2) return true;
  if(particle->nPiMinus()!=1) return true;

  // ############### //
  // Vertex decision //
  // ############### //
  TVector3 vertex;
  double vdis = 9999;
  int fproton1 = -1;
  int fproton2 = -1;
  for(int it=0; it<particle->nProton(); it++){
    pCDS* cds = particle->proton(it);
    if(vdis>9998||vdis>cds->vdis()){
      vdis = cds->vdis();
      vertex = cds->vbeam();
      if(GeomTools::GetID(vertex)==CID_Fiducial) fproton1=it;
    }
  }
  if(fproton1==-1) return true; /* there is no vertex in fiducial volume */
  fproton2=1-fproton1;
  FillHist("EventNumber",2);

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

  pCDS* proton1 = particle->proton(fproton1);
  pCDS* proton2 = particle->proton(fproton2);
  pCDS* pim     = particle->pim(0);
  pCDS* ppim1 = 0;
  pCDS* ppim2 = 0;
  for(int it=0; it<particle->nProduct(); it++){
    pCDS* pro = particle->product(it);
    if(pro->comb()!=pow(2,CDS_PiMinus)+pow(2,CDS_Proton)) continue;
    if(ppim1==0) ppim1=pro;
    elseif(ppim2==0) ppim2=pro;
  }
  if(proton1==0||proton2==0||pim==0||ppim1==0||ppim2==0) return true;

  /* Set Lorentz vector */
  /* 3He Target */
  TVector3 ZeroV;
  TLorentzVector Ltgt; Ltgt.SetVectM(ZeroV,ThreeHeMass);
  TLorentzVector Lbeam = beam->GetLorentzVector(vertex);
  TVector3 boost = (Ltgt+Lbeam).BoostVector();
  TLorentzVector cLtgt = Ltgt; cLtgt.Boost(-boost); TLorentzVector cLbeam = Lbeam; cLbeam.Boost(-boost);

  /* PiMinus */
  TLorentzVector Lpim = pim->GetLorentzVector();
  TLorentzVector cLpim = pim->GetLorentzVector(); cLpim.Boost(-boost);
  double pim_mass[2]  = { Lpim.M()                           , cLpim.M()};
  double pim_mom[2]   = { Lpim.Vect().Mag()                  , cLpim.Vect().Mag()};
  double pim_cost[2]  = { Lpim.Vect().CosTheta()             , cLpim.Vect().CosTheta()};
  double pim_phi[2]   = { Lpim.Vect().Phi()                  , cLpim.Vect().Phi()};
  double pim_mmass[2] = { (Ltgt+Lbeam-Lpim).M()              , (cLtgt+cLbeam-cLpim).M()};
  double pim_mmom[2]  = { (Ltgt+Lbeam-Lpim).Vect().Mag()     , (cLtgt+cLbeam-cLpim).Vect().Mag()};
  double pim_mcost[2] = { (Ltgt+Lbeam-Lpim).Vect().CosTheta(), (cLtgt+cLbeam-cLpim).Vect().CosTheta()};
  double pim_pbdca = pim->pbdca();

  /* Proton1 */
  TLorentzVector Lproton1 = proton1->GetLorentzVector();
  TLorentzVector cLproton1 = proton1->GetLorentzVector(); cLproton1.Boost(-boost);
  double proton1_mass[2]  = { Lproton1.M()                           , cLproton1.M()};
  double proton1_mom[2]   = { Lproton1.Vect().Mag()                  , cLproton1.Vect().Mag()};
  double proton1_cost[2]  = { Lproton1.Vect().CosTheta()             , cLproton1.Vect().CosTheta()};
  double proton1_phi[2]   = { Lproton1.Vect().Phi()                  , cLproton1.Vect().Phi()};
  double proton1_mmass[2] = { (Ltgt+Lbeam-Lproton1).M()              , (cLtgt+cLbeam-cLproton1).M()};
  double proton1_mmom[2]  = { (Ltgt+Lbeam-Lproton1).Vect().Mag()     , (cLtgt+cLbeam-cLproton1).Vect().Mag()};
  double proton1_mcost[2] = { (Ltgt+Lbeam-Lproton1).Vect().CosTheta(), (cLtgt+cLbeam-cLproton1).Vect().CosTheta()};
  double proton1_pbdca = proton1->pbdca();

  /* Proton2 */
  TLorentzVector Lproton2 = proton2->GetLorentzVector();
  TLorentzVector cLproton2 = proton2->GetLorentzVector(); cLproton2.Boost(-boost);
  double proton2_mass[2]  = { Lproton2.M()                           , cLproton2.M()};
  double proton2_mom[2]   = { Lproton2.Vect().Mag()                  , cLproton2.Vect().Mag()};
  double proton2_cost[2]  = { Lproton2.Vect().CosTheta()             , cLproton2.Vect().CosTheta()};
  double proton2_phi[2]   = { Lproton2.Vect().Phi()                  , cLproton2.Vect().Phi()};
  double proton2_mmass[2] = { (Ltgt+Lbeam-Lproton2).M()              , (cLtgt+cLbeam-cLproton2).M()};
  double proton2_mmom[2]  = { (Ltgt+Lbeam-Lproton2).Vect().Mag()     , (cLtgt+cLbeam-cLproton2).Vect().Mag()};
  double proton2_mcost[2] = { (Ltgt+Lbeam-Lproton2).Vect().CosTheta(), (cLtgt+cLbeam-cLproton2).Vect().CosTheta()};
  double proton2_pbdca = proton2->pbdca();

  /* PiMinus+Proton1 */
  TLorentzVector Lppim1 = ppim1->GetLorentzVector();
  TLorentzVector cLppim1 = ppim1->GetLorentzVector(); cLppim1.Boost(-boost);
  double ppim1_mass[2]  = { Lppim1.M()                           , cLppim1.M()};
  double ppim1_mom[2]   = { Lppim1.Vect().Mag()                  , cLppim1.Vect().Mag()};
  double ppim1_cost[2]  = { Lppim1.Vect().CosTheta()             , cLppim1.Vect().CosTheta()};
  double ppim1_phi[2]   = { Lppim1.Vect().Phi()                  , cLppim1.Vect().Phi()};
  double ppim1_mmass[2] = { (Ltgt+Lbeam-Lppim1).M()              , (cLtgt+cLbeam-cLppim1).M()};
  double ppim1_mmom[2]  = { (Ltgt+Lbeam-Lppim1).Vect().Mag()     , (cLtgt+cLbeam-cLppim1).Vect().Mag()};
  double ppim1_mcost[2] = { (Ltgt+Lbeam-Lppim1).Vect().CosTheta(), (cLtgt+cLbeam-cLppim1).Vect().CosTheta()};
  double ppim1_pbdca = ppim1->pbdca();
  double ppim1_vbdis = ppim1->vbdis();
  double ppim1_vdis  = ppim1->vdis();


  /* PiMinus+Proton2 */
  TLorentzVector Lppim2 = ppim2->GetLorentzVector();
  TLorentzVector cLppim2 = ppim2->GetLorentzVector(); cLppim2.Boost(-boost);
  double ppim2_mass[2]  = { Lppim2.M()                           , cLppim2.M()};
  double ppim2_mom[2]   = { Lppim2.Vect().Mag()                  , cLppim2.Vect().Mag()};
  double ppim2_cost[2]  = { Lppim2.Vect().CosTheta()             , cLppim2.Vect().CosTheta()};
  double ppim2_phi[2]   = { Lppim2.Vect().Phi()                  , cLppim2.Vect().Phi()};
  double ppim2_mmass[2] = { (Ltgt+Lbeam-Lppim2).M()              , (cLtgt+cLbeam-cLppim2).M()};
  double ppim2_mmom[2]  = { (Ltgt+Lbeam-Lppim2).Vect().Mag()     , (cLtgt+cLbeam-cLppim2).Vect().Mag()};
  double ppim2_mcost[2] = { (Ltgt+Lbeam-Lppim2).Vect().CosTheta(), (cLtgt+cLbeam-cLppim2).Vect().CosTheta()};
  double ppim2_pbdca = ppim2->pbdca();
  double ppim2_vbdis = ppim2->vbdis();
  double ppim2_vdis  = ppim2->vdis();

  /* Proton1+(PiMinus+Proton2) */
  pCDS* pppim = Calcpppim(pBeam,proton1,ppim2);
  TLorentzVector Lpppim = pppim->GetLorentzVector();
  TLorentzVector cLpppim = pppim->GetLorentzVector(); cLpppim.Boost(-boost);
  double pppim_mass[2]  = { Lpppim.M()                           , cLpppim.M()};
  double pppim_mom[2]   = { Lpppim.Vect().Mag()                  , cLpppim.Vect().Mag()};
  double pppim_cost[2]  = { Lpppim.Vect().CosTheta()             , cLpppim.Vect().CosTheta()};
  double pppim_phi[2]   = { Lpppim.Vect().Phi()                  , cLpppim.Vect().Phi()};
  double pppim_mmass[2] = { (Ltgt+Lbeam-Lpppim).M()              , (cLtgt+cLbeam-cLpppim).M()};
  double pppim_mmom[2]  = { (Ltgt+Lbeam-Lpppim).Vect().Mag()     , (cLtgt+cLbeam-cLpppim).Vect().Mag()};
  double pppim_mcost[2] = { (Ltgt+Lbeam-Lpppim).Vect().CosTheta(), (cLtgt+cLbeam-cLpppim).Vect().CosTheta()};
  double pppim_pbdca = pppim->pbdca();
  double pppim_vbdis = pppim->vbdis();
  double pppim_vdis  = pppim->vdis();

  return true;
}

pCDS* MyAnalysisLp::Calcpppim(pBeam* beam,pCDS* proton,pCDS* lambda,CDSTrackMan* cdstrackMan)
{
  TLorentzVector L1 = proton->GetLorentzVector();
  TLorentzVector L2 = lambda->GetLorentzVector();
  double im = (L1+L2).M();
  TVector3 Pp = (L1+L2).Vect();

  CDSTrack* cds = cdstrackMan->GoodTrack(proton->id());
  TVector3 lpos = lambda->vertex();
  TVector3 ldir = lambda->momdir();
  TVector3 vtx1,vtx2;
  double dis;
  if(!TrackTools::CalcLinehilixVertex(lpos,ldir,cds,vtx1,vtx2,dis)) return 0;

  double dist,dltmp=0;
  TVector3 xest,nest;
  MathTools::LineToLine( vtx,Pp.Unit(),beam->bpcpos(), beam->bpcdir(),dltmp,dist,xest,nest );
  pCDS* pro=new pCDS();
    
  //pro->SetDaughterID1(cds1->id());
  //pro->SetDaughterID2(cds2->id());
  //pro->SetCombID((int)(pow(2,cds1->pid())+pow(2,cds2->pid())));
  pro->SetMomentum(Pp.Mag());
  pro->SetMass(im);
  pro->SetVertex(vtx);
  pro->SetVertexBeam(xest);
  pro->SetMomDir(Pp.Unit());
  pro->SetVDis(dist1);
  pro->SetVBDis((xest-vtx).Mag());
  pro->SetPBDCA(dist);
  pro->SetOA(Pp1.Angle(Pp2));
  pro->SetAngleLab(beam->bpcdir().Angle(Pp));
  
  return pro;
}

bool MyAnalysisLp::FillHist(TString name, double val1, int weight)
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

bool MyAnalysisLp::FillHist(TString name, TString val1, int weight)
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

bool MyAnalysisLp::FillHist(TString name, double val1, double val2, int weight)
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

bool MyAnalysisLp::FillHist(TString name, TString val1, TString val2, int weight)
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

void MyAnalysisLp::CutCondition()
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

bool MyAnalysisLp::Initialize(ConfMan* confMan)
{
  std::cout << "### MyAnalysisLp::Initialize ###" << std::endl;

  std::string ofname = confMan->GetOutFileName();
  ofname.insert(ofname.find(".root"),"_anaLp");

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
  new TH2F( "FWD_n_OverbetavsMomentum", "1/#beta vs. Momentum;1/#beta;Momentum (GeV/c)", 1000, -0, 10, 500, -5, 5 );
  new TH2F( "FWD_n_OverbetavsEnergy", "1/#beta vs. Energy deposit;1/#beta;Energy deposit (MeVee)", 1000, -0, 10, 5000, 0, 100);
  new TH2F( "FWD_n_TOFvsMomentum", "TOF vs. Momentum;TOF (ns);Momentum (GeV/c)", 4000, 30, 70, 500, -5, 5 );
  new TH2F( "FWD_n_HitPosition", "Hit position at NC (X vs. Y);X position (cm);Y position (cm)",40,-400,400,25,-250,250);
  new TH2F( "FWD_n_HitSegment", "NC segment;Segment in layer;Layer",16,0.5,16.5,7,0.5,7.5);
  new TH2F( "FWD_k0_OverbetavsMomentum", "1/#beta vs. Momentum;1/#beta;Momentum (GeV/c)", 1000, -0, 10, 500, -5, 5 );
  new TH2F( "FWD_k0_OverbetavsEnergy", "1/#beta vs. Energy deposit;1/#beta;Energy deposit (MeVee)", 1000, -0, 10, 5000, 0, 100);
  new TH2F( "FWD_k0_TOFvsMomentum", "TOF vs. Momentum;TOF (ns);Momentum (GeV/c)", 4000, 30, 70, 500, -5, 5 );
  new TH2F( "FWD_k0_HitPosition", "Hit position at NC (X vs. Y);X position (cm);Y position (cm)",40,-400,400,25,-250,250);
  new TH2F( "FWD_k0_HitSegment", "NC segment;Segment in layer;Layer",16,0.5,16.5,7,0.5,7.5);
  // neutron
  new TH2F("n_CosTvsMomentum_CM","cos#theta vs. momentum of n;cos#theta;Momentum (GeV/c)", 2000, -1.0, 1.0, 1000, 0.0, 1.0);
  new TH2F("n_CosTvsPhi_CM","cos#theta vs. #phi of n;cos#theta;#phi", 2000, -1.0, 1.0, 1600, 0.0, 3.2);
  for(int ip=0; ip<8; ip++){
    new TH2F(Form("n_%s_CosTvsMomentum_CM",pname[ip].Data()),"cos#theta vs. momentum of n;cos#theta;Momentum (GeV/c)", 2000, -1.0, 1.0, 1000, 0.0, 1.0);
    new TH2F(Form("n_%s_CosTvsPhi_CM",pname[ip].Data()),"cos#theta vs. #phi of n;cos#theta;#phi", 2000, -1.0, 1.0, 1600, 0.0, 3.2);
  }
  for(int ip=0; ip<6; ip++){
    new TH2F(Form("n_%s_CosTvsMomentum_CM",proname[ip].Data()),"cos#theta vs. momentum of n;cos#theta;Momentum (GeV/c)", 2000, -1.0, 1.0, 1000, 0.0, 1.0);
    new TH2F(Form("n_%s_CosTvsPhi_CM",proname[ip].Data()),"cos#theta vs. #phi of n;cos#theta;#phi", 2000, -1.0, 1.0, 1600, 0.0, 3.2);
  }
  // Invariant mass
  std::cout << "Define Histograms for IM(nX) M.M." << std::endl;
  for(int ip=0; ip<8; ip++){
    new TH2F(Form("n%s_IMvsMomentum_CM",pname[ip].Data()), Form("Invariant mass (n%s) vs Momentum;IM(n%s) (GeV/c^{2});Momentum (GeV/c)",pname2[ip].Data(),pname2[ip].Data()), 2000, 0.0, 2.0, 1000, 0.0, 1.0 );
    new TH2F(Form("n%s_IMvsCosT_CM",pname[ip].Data()), Form("Invariant mass (n%s) vs cos#theta;IM(n%s) (GeV/c^{2});cos#theta",pname2[ip].Data(),pname2[ip].Data()), 2000, 0.0, 2.0, 2000, -1.0, 1.0 );
  }
  // Missing mass
  std::cout << "Define Histograms for p(K-,n)X M.M." << std::endl;
  new TH2F("pKn_MMvsMomentum_CM", "p(K^{-},n)X missing mass vs. missing momentum;MM_{p}(n) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKn_MMvsCosT_CM", "p(K^{-},n)X missing mass vs. cos#theta;MM_{p}(n) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  for(int ip=0; ip<8; ip++){
    new TH2F(Form("pKn_%s_MMvsMomentum_CM",pname[ip].Data()), "p(K^{-},n)X missing mass vs. missing momentum;MM_{p}(n) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
    new TH2F(Form("pKn_%s_MMvsCosT_CM",pname[ip].Data()), "p(K^{-},n)X missing mass vs. cos#theta;MM_{p}(n) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
    new TH2F(Form("pKn%s_MMvsMomentum_CM",pname[ip].Data()), Form("p(K^{-},n%s)X missing mass vs. missing momentum;MM_{p}(n%s) (GeV/c^{2});Missing momentum (GeV/c)",pname2[ip].Data(),pname2[ip].Data()), 2000, 0.0, 2.0, 1000, 0.0, 1.0);
    new TH2F(Form("pKn%s_MMvsCosT_CM",pname[ip].Data()), Form("p(K^{-},n%s)X missing mass vs. cos#theta;MM_{p}(n%s) (GeV/c^{2});cos#theta",pname2[ip].Data(),pname2[ip].Data()), 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  }
  for(int ip=0; ip<6; ip++){
    new TH2F(Form("pKn_%s_MMvsMomentum_CM",proname[ip].Data()), "p(K^{-},n)X missing mass vs. missing momentum;MM_{p}(n) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
    new TH2F(Form("pKn_%s_MMvsCosT_CM",proname[ip].Data()), "p(K^{-},n)X missing mass vs. cos#theta;MM_{p}(n) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
    new TH2F(Form("pKn%s_MMvsMomentum_CM",proname[ip].Data()), Form("p(K^{-},n%s)X missing mass vs. missing momentum;MM_{p}(n%s) (GeV/c^{2});Missing momentum (GeV/c)",proname2[ip].Data(),proname2[ip].Data()), 2000, 0.0, 2.0, 1000, 0.0, 1.0);
    new TH2F(Form("pKn%s_MMvsCosT_CM",proname[ip].Data()), Form("p(K^{-},n%s)X missing mass vs. cos#theta;MM_{p}(n%s) (GeV/c^{2});cos#theta",proname2[ip].Data(),proname2[ip].Data()), 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  }

  return true;
}
