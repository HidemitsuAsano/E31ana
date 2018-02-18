// MyAnalysisOverVeto.cpp

#include "MyAnalysisOverVeto.h"

MyAnalysisOverVeto::MyAnalysisOverVeto(TFile* rt, ConfMan* conf)
{
  Initialize(conf);
  CutCondition();
  Clear();
}

MyAnalysisOverVeto::~MyAnalysisOverVeto()
{
  Clear();
  rtFile->cd();
  rtFile->Write();
  rtFile->Close();
}

void MyAnalysisOverVeto::Clear()
{
}

bool MyAnalysisOverVeto::DoAnalysis(ConfMan* conf, EventHeader* header, BeamLineHitMan* blMan, BeamLineTrackMan* bltrackMan, CDSHitMan* cdsMan, CDSTrackingMan* cdstrackMan, Particle* particle)
{
  rtFile->cd();

  FillHist("EventNumber",0);
  DetectorList *dlist=DetectorList::GetInstance();
  if(particle->nBeam()!=1) return false;
  pBeam* beam = particle->beam(0);
  FillHist("EventNumber",1);

  /* Event selection */
  //if(particle->nCDS()!=2) return true;

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

  // ####### //
  // Only NC //
  // ####### //
  /* charge veto with NC layer 1 */
  bool Nflag_NC = true;
  for(int i=0; i<blMan->nNC(); i++){
    HodoscopeLikeHit* hit = blMan->NC(i);
    int seg = (hit->seg()-1)%16+1;
    int lay = (hit->seg()-1)/16+1;
    if(lay==1&&hit->CheckRange()){
      FillHist(Form("NCL1_Time"),hit->ctmean());
      FillHist(Form("NCL1S%02d_Time",seg),hit->ctmean());
      if(60.0<hit->ctmean()&&hit->ctmean()<80.0){
        Nflag_NC = false;
        break;
      }
    }
  }
  if(!Nflag_NC) return true;
  /* NC hit refill */
  Particle* particleNC = new Particle();
  for(int i=0; i<blMan->nNC(); i++){
    HodoscopeLikeHit* hit = blMan->NC(i);
    int lay = (hit->seg()-1)/16+1;
    if(lay==1) continue;
    if(hit->CheckRange()){
      pNC* nc = new pNC();
      nc->SetSegment(hit->seg());
      nc->SetHitPosition(hit->pos().X(),hit->hitpos(),hit->pos().Z());
      nc->SetTime(hit->ctmean());
      nc->SetEnergy(hit->emean());
      particleNC->AddNC(*nc);
      delete nc;
    }
  }

  // ############### //
  // NC hit decision //
  // ############### //
  double time = 9999;
  int fnc = -1;
  for(int it=0; it<particleNC->nNC(); it++){
    pNC* nc = particleNC->nc(it);
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

  pCDS* cds=0;
  if(fcds!=-1) {
    cds=particle->cdsi(fcds);
  }
  if(cds==0) return false;

  TVector3 ZeroV;
  TLorentzVector Ltgt; Ltgt.SetVectM(ZeroV,pMass);
  TLorentzVector Lbeam = beam->GetLorentzVector(vertex);
  TVector3 boost = (Ltgt+Lbeam).BoostVector();
  TLorentzVector cLtgt = Ltgt; cLtgt.Boost(-boost); TLorentzVector cLbeam = Lbeam; cLbeam.Boost(-boost);
  TLorentzVector Lcds = cds->GetLorentzVector();
  TLorentzVector cLcds = cds->GetLorentzVector(); cLcds.Boost(-boost);
  TString frame[2] = {"Lab","CM"};
  TString pname[8] = {"pip","p","d","t","he","pim","k","o"};

  pNC*  nc =0;
  if(fnc!=-1)  nc =particleNC->nc(fnc);
  if(nc==0) return true;
  FillHist("EventNumber",5);

  TLorentzVector Ln  = nc ->GetLorentzVector();
  TLorentzVector cLn  = nc ->GetLorentzVector(); cLn.Boost(-boost);
  double seg = (nc->seg()-1)%16+1;
  double lay = (nc->seg()-1)/16+1;

  FillHist("FWD_n_OverbetavsMomentum",1.0/nc->beta(),nc->mom().Mag());
  FillHist("FWD_n_OverbetavsEnergy",1.0/nc->beta(),nc->energy());
  FillHist("FWD_n_TOFvsMomentum",nc->tof(),nc->mom().Mag());
  FillHist("FWD_n_HitPosition",nc->hitpos().X(),nc->hitpos().Y());
  FillHist("FWD_n_HitSegment",seg,lay);


  /* BVC Veto */
  int MulBVC=0;
  for(int i=0; i<blMan->nBVC(); i++){
    HodoscopeLikeHit* hit = blMan->BVC(i);
    if(hit->CheckRange()){
      MulBVC++;
      FillHist("BVC_Time",hit->ctmean());
    }
  }
  if(MulBVC==0){
    FillHist("FWD_BVC_OverbetavsMomentum",1.0/nc->beta(),nc->mom().Mag());
    FillHist("FWD_BVC_OverbetavsEnergy",1.0/nc->beta(),nc->energy());
    FillHist("FWD_BVC_TOFvsMomentum",nc->tof(),nc->mom().Mag());
    FillHist("FWD_BVC_HitPosition",nc->hitpos().X(),nc->hitpos().Y());
    FillHist("FWD_BVC_HitSegment",seg,lay);
  }
  /* CVC Veto */
  int MulCVC=0;
  for(int i=0; i<blMan->nCVC(); i++){
    HodoscopeLikeHit* hit = blMan->CVC(i);
    if(hit->CheckRange()){
      MulCVC++;
      FillHist("CVC_Time",hit->ctmean());
    }
  }
  if(MulCVC==0){
    FillHist("FWD_CVC_OverbetavsMomentum",1.0/nc->beta(),nc->mom().Mag());
    FillHist("FWD_CVC_OverbetavsEnergy",1.0/nc->beta(),nc->energy());
    FillHist("FWD_CVC_TOFvsMomentum",nc->tof(),nc->mom().Mag());
    FillHist("FWD_CVC_HitPosition",nc->hitpos().X(),nc->hitpos().Y());
    FillHist("FWD_CVC_HitSegment",seg,lay);
  }
  /* BVC and CVC Veto */
  if(MulBVC==0&&MulCVC==0){
    FillHist("FWD_BVCCVC_OverbetavsMomentum",1.0/nc->beta(),nc->mom().Mag());
    FillHist("FWD_BVCCVC_OverbetavsEnergy",1.0/nc->beta(),nc->energy());
    FillHist("FWD_BVCCVC_TOFvsMomentum",nc->tof(),nc->mom().Mag());
    FillHist("FWD_BVCCVC_HitPosition",nc->hitpos().X(),nc->hitpos().Y());
    FillHist("FWD_BVCCVC_HitSegment",seg,lay);
  }



  return true;

}

bool MyAnalysisOverVeto::FillHist(TString name, double val1, int weight)
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

bool MyAnalysisOverVeto::FillHist(TString name, TString val1, int weight)
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

bool MyAnalysisOverVeto::FillHist(TString name, double val1, double val2, int weight)
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

bool MyAnalysisOverVeto::FillHist(TString name, TString val1, TString val2, int weight)
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

void MyAnalysisOverVeto::CutCondition()
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
  sigma = 0.0040;
  spll = mean-2*sigma; spul = mean+2*sigma;
  sblspll = mean-5*sigma; sblspul = mean-3*sigma;
  sbuspll = mean+3*sigma; sbuspul = mean+5*sigma;
  /* SigmaMinus mass cut */
  mean  = 1.1972;
  sigma = 0.0071;
  smll = mean-2*sigma; smul = mean+2*sigma;
  sblsmll = mean-5*sigma; sblsmul = mean-3*sigma;
  sbusmll = mean+3*sigma; sbusmul = mean+5*sigma;
}

bool MyAnalysisOverVeto::Initialize(ConfMan* confMan)
{
  std::cout << "### MyAnalysisOverVeto::Initialize ###" << std::endl;

  std::string ofname = confMan->GetOutFileName();
  ofname.insert(ofname.find(".root"),"_anaOverVeto");

  rtFile =  new TFile( Form("%s",ofname.c_str()), "RECREATE");
  rtFile -> cd();

  TString pname[8]  = {"pip","p","d","t","he","pim","k","o"};
  TString pname2[8] = {"#pi^{+}","p","d","t","he","#pi^{-}","K^{-}","o"};

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
  // NC hit 
  new TH1F(Form("NCL1_Time"),Form("NC Layer-1 Time;Time (ns);Counts"),3000,0,300);
  for(int seg=1; seg<=16; seg++){
    new TH1F(Form("NCL1S%02d_Time",seg),Form("NC Layer-1 Segment-%d Time;Time (ns);Counts",seg),3000,0,300);
  }
  // BVC,CVC hit 
  new TH1F("BVC_Time","BVC Time;Time (ns);Counts",3000,0,300);
  new TH1F("CVC_Time","CVC Time;Time (ns);Counts",3000,0,300);
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
  new TH2F( "FWD_BVC_OverbetavsMomentum", "1/#beta vs. Momentum;1/#beta;Momentum (GeV/c)", 1000, -0, 10, 500, -5, 5 );
  new TH2F( "FWD_BVC_OverbetavsEnergy", "1/#beta vs. Energy deposit;1/#beta;Energy deposit (MeVee)", 1000, -0, 10, 5000, 0, 100);
  new TH2F( "FWD_BVC_TOFvsMomentum", "TOF vs. Momentum;TOF (ns);Momentum (GeV/c)", 4000, 30, 70, 500, -5, 5 );
  new TH2F( "FWD_BVC_HitPosition", "Hit position at NC (X vs. Y);X position (cm);Y position (cm)",40,-400,400,25,-250,250);
  new TH2F( "FWD_BVC_HitSegment", "NC segment;Segment in layer;Layer",16,0.5,16.5,7,0.5,7.5);
  new TH2F( "FWD_CVC_OverbetavsMomentum", "1/#beta vs. Momentum;1/#beta;Momentum (GeV/c)", 1000, -0, 10, 500, -5, 5 );
  new TH2F( "FWD_CVC_OverbetavsEnergy", "1/#beta vs. Energy deposit;1/#beta;Energy deposit (MeVee)", 1000, -0, 10, 5000, 0, 100);
  new TH2F( "FWD_CVC_TOFvsMomentum", "TOF vs. Momentum;TOF (ns);Momentum (GeV/c)", 4000, 30, 70, 500, -5, 5 );
  new TH2F( "FWD_CVC_HitPosition", "Hit position at NC (X vs. Y);X position (cm);Y position (cm)",40,-400,400,25,-250,250);
  new TH2F( "FWD_CVC_HitSegment", "NC segment;Segment in layer;Layer",16,0.5,16.5,7,0.5,7.5);
  new TH2F( "FWD_BVCCVC_OverbetavsMomentum", "1/#beta vs. Momentum;1/#beta;Momentum (GeV/c)", 1000, -0, 10, 500, -5, 5 );
  new TH2F( "FWD_BVCCVC_OverbetavsEnergy", "1/#beta vs. Energy deposit;1/#beta;Energy deposit (MeVee)", 1000, -0, 10, 5000, 0, 100);
  new TH2F( "FWD_BVCCVC_TOFvsMomentum", "TOF vs. Momentum;TOF (ns);Momentum (GeV/c)", 4000, 30, 70, 500, -5, 5 );
  new TH2F( "FWD_BVCCVC_HitPosition", "Hit position at NC (X vs. Y);X position (cm);Y position (cm)",40,-400,400,25,-250,250);
  new TH2F( "FWD_BVCCVC_HitSegment", "NC segment;Segment in layer;Layer",16,0.5,16.5,7,0.5,7.5);

  return true;
}
