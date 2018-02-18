// MyAnalysisNCCalib.cpp

#include "MyAnalysisNCCalib.h"

MyAnalysisNCCalib::MyAnalysisNCCalib(TFile* rt, ConfMan* conf)
{
  Initialize(conf);
  CutCondition();
  Clear();
}

MyAnalysisNCCalib::~MyAnalysisNCCalib()
{
  Clear();
  rtFile->cd();
  rtFile->Write();
  rtFile->Close();
}

void MyAnalysisNCCalib::Clear()
{
}

bool MyAnalysisNCCalib::DoAnalysis(ConfMan* conf, EventHeader* header, BeamLineHitMan* blMan, BeamLineTrackMan* bltrackMan, CDSHitMan* cdsMan, CDSTrackingMan* cdstrackMan, Particle* particle)
{
  rtFile->cd();

  DetectorList *dlist=DetectorList::GetInstance();
  if(particle->nBeam()!=1) return false;
  pBeam* beam = particle->beam(0);

  int T0Segment = beam->t0seg();
  double T0Timing = beam->t0time();
  double T0dEu=0,T0dEd=0,T0dEm=0;
  int MulT0=0;
  for(int i=0; i<blMan->nT0(); i++){
    HodoscopeLikeHit* hit = blMan->T0(i);
    if(hit->CheckRange()){
      MulT0++;
      T0dEu = hit->ene(0);
      T0dEd = hit->ene(1);
      T0dEm = hit->emean();
      T0Timing = hit->ctmean();
    }
  }
  if(MulT0!=1) return true;
  //if(T0Segment!=3) return true;
  beam->SetT0Time(T0Timing);

  // ############### //
  // Vertex decision //
  // ############### //
  TVector3 vertex;
  double vdis = 9999;
  int fpro = -1;
  for(int it=0; it<particle->nCDS(); it++){
    pCDS* pro = particle->cdsi(it);
    if(vdis>9998||vdis>pro->vdis()){
      vdis = pro->vdis();
      vertex = pro->vbeam();
      if(GeomTools::GetID(vertex)==CID_Fiducial) fpro=it;
    }
  }
  if(fpro==-1) return true; /* there is no vertex in fiducial volume */
  pCDS* cds = particle->cdsi(fpro);

  // ################# //
  // Neutral selection //
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
  if(MulBVC!=0||MulCVC!=0||MulPC!=0) return true;

  // ########### //
  // NC analysis //
  // ########### //
  HodoscopeLikeHit* nchit = 0;
  int nhit = 0;
  double tmp_t = 99999.9;
  const int cid = CID_NC;
  const char* name = dlist -> GetName(cid).c_str();
  const int nsegs = dlist -> GetNsegs(cid);
  double time_beam = beam->CalcVertexTime(vertex);
  for(int i=0; i<blMan->nNC(); i++){
    HodoscopeLikeHit* hit = blMan->NC(i);
    if(!hit->CheckRange()) continue;
    TVector3 hitpos(hit->pos().X(),hit->hitpos(),hit->pos().Z());
    double fl = (vertex-hitpos).Mag();
    double tof = hit->ctmean() - time_beam;
    double beta = fl/tof/(Const*100.0);
    FillHist("FWD_TOF_All",tof,1);
    FillHist("FWD_Overbeta_All",1.0/beta,1);
    FillHist("FWD_TOFvsdE_All",tof,hit->emean(),1);
    FillHist("FWD_OverbetavsdE_All",1.0/beta,hit->emean(),1);
    if(tmp_t>hit->ctmean()){
      tmp_t = hit->ctmean();
      nchit = hit;
    }

    if(1.0/beta<0.9||1.1<1.0/beta) continue;
    double tof_c = fl/(Const*100.0);
    double offset = tof-tof_c;
    FillHist(Form("%s_TOF",name),tof);
    FillHist(Form("%s_dEvsTOF",name),hit->emean(),tof);
    FillHist(Form("%s_dEvsTimeOffset",name),hit->emean(),offset);

    FillHist(Form("%s%d_TOF",name,hit->seg()),tof);
    FillHist(Form("%s%d_dEvsTOF",name,hit->seg()),hit->emean(),tof);
    FillHist(Form("%s%d_dEvsTimeOffset",name,hit->seg()),hit->emean(),offset);
    FillHist(Form("%s%d_CTimeSub",name,hit->seg()),hit->ctsub());

    FillHist(Form("%su%d_TOF",name,hit->seg()),tof);
    FillHist(Form("%su%d_dEvsTOF",name,hit->seg()),hit->ene(0),tof);
    FillHist(Form("%sd%d_TOF",name,hit->seg()),tof);
    FillHist(Form("%sd%d_dEvsTOF",name,hit->seg()),hit->ene(1),tof);
    FillHist(Form("%su%d_TimeOffset",name,hit->seg()),offset);
    FillHist(Form("%su%d_dEvsTimeOffset",name,hit->seg()),hit->ene(0),offset);
    FillHist(Form("%sd%d_dEvsTimeOffset",name,hit->seg()),hit->ene(1),offset);
    if(hit->emean()>8.0){
      FillHist(Form("%s_TimeOffset",name),offset);
      FillHist(Form("%s%d_TimeOffset",name,hit->seg()),offset);
      FillHist(Form("%sd%d_TimeOffset",name,hit->seg()),offset);
    }

    // For T0
    FillHist(Form("T0%d_TOF",T0Segment),tof);
    FillHist(Form("T0%d_dEvsTOF",T0Segment),T0dEm,tof);
    FillHist(Form("T0%d_TimeOffset",T0Segment),-offset);
    FillHist(Form("T0%d_dEvsTimeOffset",T0Segment),T0dEm,-offset);

    FillHist(Form("T0u%d_TOF",T0Segment),tof);
    FillHist(Form("T0u%d_dEvsTOF",T0Segment),T0dEu,tof);
    FillHist(Form("T0d%d_TOF",T0Segment),tof);
    FillHist(Form("T0d%d_dEvsTOF",T0Segment),T0dEd,tof);
    FillHist(Form("T0u%d_TimeOffset",T0Segment),-offset);
    FillHist(Form("T0u%d_dEvsTimeOffset",T0Segment),T0dEu,-offset);
    FillHist(Form("T0d%d_TimeOffset",T0Segment),-offset);
    FillHist(Form("T0d%d_dEvsTimeOffset",T0Segment),T0dEd,-offset);


  }
  if(nchit==0) return true;

  pNC* nc = new pNC();
  nc->SetSegment(nchit->seg());
  nc->SetHitPosition(nchit->pos().X(),nchit->hitpos(),nchit->pos().Z());
  nc->SetTime(nchit->ctmean());
  nc->SetEnergy(nchit->emean());
  nc->CalcMom(beam,vertex);

  FillHist("FWD_TOF",nc->tof(),1);
  FillHist("FWD_Overbeta",1.0/nc->beta(),1);
  FillHist("FWD_TOFvsdE",nc->tof(),nc->energy(),1);
  FillHist("FWD_OverbetavsdE",1.0/nc->beta(),nc->energy(),1);
  FillHist("FWD_HitPosition",nc->hitpos().X(),nc->hitpos().Y(),1);
  double seg = (nc->seg()-1)%16+1;
  double lay = (nc->seg()-1)/16+1;
  FillHist("FWD_HitSegment",seg,lay,1);

  if(nc->energy()<8.0) return true;
  FillHist(Form("NC%d_Overbeta",nc->seg()),1.0/nc->beta(),1);
  FillHist(Form("T0%d_Overbeta",T0Segment),1.0/nc->beta(),1);
  if(nc->pid()!=F_Neutron) return true;

  /* Target */
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
  /* Beam */
  TLorentzVector Lbeam = beam->GetLorentzVector(vertex);
  TVector3 boost = (Ltgt+Lbeam).BoostVector();
  TLorentzVector cLtgt = Ltgt; cLtgt.Boost(-boost); TLorentzVector cLbeam = Lbeam; cLbeam.Boost(-boost);
  /* CDS */
  TLorentzVector Lcds = cds->GetLorentzVector();
  TLorentzVector cLcds = cds->GetLorentzVector(); cLcds.Boost(-boost);
  /* NC */
  TLorentzVector Lnc = nc->GetLorentzVector();
  TLorentzVector cLnc = nc->GetLorentzVector(); cLnc.Boost(-boost);

  TString frame[2] = {"Lab","CM"};
  TString pname[8] = {"pip","p","d","t","he","pim","k","o"};
  double cds_mass[2]  = { (Lcds).M()                    , (cLcds).M()};
  double cds_mom[2]   = { (Lcds).Vect().Mag()           , (cLcds).Vect().Mag()};
  double cds_cost[2]  = { (Lcds).Vect().CosTheta()      , (cLcds).Vect().CosTheta()};
  double cds_phi[2]   = { (Lcds).Vect().Phi()           , (cLcds).Vect().Phi()};
  double cds_mmass[2] = { (Ltgt+Lbeam-Lcds).M()              , (cLtgt+cLbeam-cLcds).M()};
  double cds_mmom[2]  = { (Ltgt+Lbeam-Lcds).Vect().Mag()     , (cLtgt+cLbeam-cLcds).Vect().Mag()};
  double cds_mcost[2] = { (Ltgt+Lbeam-Lcds).Vect().CosTheta(), (cLtgt+cLbeam-cLcds).Vect().CosTheta()};
  double nc_mass[2]    = { Lnc.M()                  , cLnc.M()};
  double nc_mom[2]    = { Lnc.Vect().Mag()                  , cLnc.Vect().Mag()};
  double nc_cost[2]   = { Lnc.Vect().CosTheta()             , cLnc.Vect().CosTheta()};
  double nc_phi[2]    = { Lnc.Vect().Phi()                  , cLnc.Vect().Phi()};
  double nc_mmass[2]  = { (Ltgt+Lbeam-Lnc).M()              , (cLtgt+cLbeam-cLnc).M()};
  double nc_mmom[2]   = { (Ltgt+Lbeam-Lnc).Vect().Mag()     , (cLtgt+cLbeam-cLnc).Vect().Mag()};
  double nc_mcost[2]  = { (Ltgt+Lbeam-Lnc).Vect().CosTheta(), (cLtgt+cLbeam-cLnc).Vect().CosTheta()};
  double ncds_mass[2] = { (Lcds+Lnc).M()                           , (cLcds+cLnc).M()};
  double ncds_mom[2]  = { (Lcds+Lnc).Vect().Mag()                  , (cLcds+cLnc).Vect().Mag()};
  double ncds_cost[2] = { (Lcds+Lnc).Vect().CosTheta()             , (cLcds+cLnc).Vect().CosTheta()};
  double ncds_phi[2]  = { (Lcds+Lnc).Vect().Phi()                  , (cLcds+cLnc).Vect().Phi()};
  double ncds_mmass[2]= { (Ltgt+Lbeam-(Lcds+Lnc)).M()              , (cLtgt+cLbeam-(cLcds+cLnc)).M()};
  double ncds_mmom[2] = { (Ltgt+Lbeam-(Lcds+Lnc)).Vect().Mag()     , (cLtgt+cLbeam-(cLcds+cLnc)).Vect().Mag()};
  double ncds_mcost[2]= { (Ltgt+Lbeam-(Lcds+Lnc)).Vect().CosTheta(), (cLtgt+cLbeam-(cLcds+cLnc)).Vect().CosTheta()};

  for(int f=1; f<2; f++){
    FillHist(Form("n_Mass_%s",pname[cds->pid()].Data(),frame[f].Data()),nc_mass[f]);
    FillHist(Form("n_Momentum_%s",pname[cds->pid()].Data(),frame[f].Data()),nc_mom[f]);
    FillHist(Form("n_CosT_%s",pname[cds->pid()].Data(),frame[f].Data()),nc_cost[f]);
    FillHist(Form("n_MMass_%s",pname[cds->pid()].Data(),frame[f].Data()),nc_mmass[f]);
    FillHist(Form("n_MMomentum_%s",pname[cds->pid()].Data(),frame[f].Data()),nc_mmom[f]);
    FillHist(Form("n_MCosT_%s",pname[cds->pid()].Data(),frame[f].Data()),nc_mcost[f]);
  }
  for(int f=1; f<2; f++){
    FillHist(Form("%s_Mass_%s",pname[cds->pid()].Data(),frame[f].Data()),ncds_mass[f]);
    FillHist(Form("%s_Momentum_%s",pname[cds->pid()].Data(),frame[f].Data()),ncds_mom[f]);
    FillHist(Form("%s_CosT_%s",pname[cds->pid()].Data(),frame[f].Data()),ncds_cost[f]);
    FillHist(Form("%s_MMass_%s",pname[cds->pid()].Data(),frame[f].Data()),ncds_mmass[f]);
    FillHist(Form("%s_MMomentum_%s",pname[cds->pid()].Data(),frame[f].Data()),ncds_mmom[f]);
    FillHist(Form("%s_MCosT_%s",pname[cds->pid()].Data(),frame[f].Data()),ncds_mcost[f]);
  }


  return true;
}

bool MyAnalysisNCCalib::FillHist(TString name, double val1, int weight)
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

bool MyAnalysisNCCalib::FillHist(TString name, TString val1, int weight)
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

bool MyAnalysisNCCalib::FillHist(TString name, double val1, double val2, int weight)
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

bool MyAnalysisNCCalib::FillHist(TString name, TString val1, TString val2, int weight)
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

void MyAnalysisNCCalib::CutCondition()
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
  /* Missing k0 mass cut */ 
  mean  = 0.4970;
  sigma = 0.0130;
  mk0ll = mean-2*sigma; mk0ul = mean+2*sigma;
  sblmk0ll = mean-5*sigma; sblmk0ul = mean-3*sigma;
  sbumk0ll = mean+3*sigma; sbumk0ul = mean+5*sigma;
}

bool MyAnalysisNCCalib::Initialize(ConfMan* confMan)
{
  std::cout << "### MyAnalysisNCCalib::Initialize ###" << std::endl;

  DetectorList *dlist=DetectorList::GetInstance();
  std::string ofname = confMan->GetOutFileName();
  ofname.insert(ofname.find(".root"),"_anaNCCalib");

  rtFile =  new TFile( Form("%s",ofname.c_str()), "RECREATE");
  rtFile -> cd();

  // FWD neutral Particles
  std::cout << "Define Histograms for FWD neutral particles" << std::endl;
  //new TH2F( "FWD_OverbetavsdE_All", "1/#beta vs. Energy deposit;1/#beta;Energy deposit (MeVee)", 5000, 0, 10, 1500, 0, 150);
  //new TH2F( "FWD_TOFvsdE_All", "TOF vs. Energy deposit;TOF (ns);Energy deposit (MeVee)", 2000, 30, 130, 1500, 0, 150 );
  //new TH1F( "FWD_Overbeta_All", "1/#beta;1/#beta;Counts", 5000, 0, 10);
  //new TH1F( "FWD_TOF_All", "TOF;TOF (ns);Counts", 2000, 30, 130);
  new TH2F( "FWD_OverbetavsdE", "1/#beta vs. Energy deposit;1/#beta;Energy deposit (MeVee)", 3000, 0, 5, 1500, 0, 150);
  new TH2F( "FWD_TOFvsdE", "TOF vs. Energy deposit;TOF (ns);Energy deposit (MeVee)", 1000, 40, 100, 1500, 0, 150 );
  new TH1F( "FWD_Overbeta", "1/#beta;1/#beta;Counts", 3000, 0, 5);
  new TH1F( "FWD_TOF", "TOF;TOF (ns);Counts", 1000, 40, 100);
  new TH2F( "FWD_HitPosition", "Hit position at NC (X vs. Y);X position (cm);Y position (cm)",40,-400,400,25,-250,250);
  new TH2F( "FWD_HitSegment", "NC segment;Segment in layer;Layer",18,-0.5,17.5,9,-0.5,8.5);

  // CDS 1 particle
  TString pname[9]  = {"pip","p","d","t","he","pim","k","o","n"};
  TString pname2[9] = {"#pi^{+}","p","d","t","he","#pi^{-}","K^{-}","o","n"};
  for(int ip=8; ip<9; ip++){
    new TH1F( Form("%s_Mass_CM",pname[ip].Data()), Form("Invariant Mass (%s);IM(%s) (GeV/c^{2});Counts",pname2[ip].Data(),pname2[ip].Data()), 5000, 0.0, 5.0);
    new TH1F( Form("%s_Momentum_CM",pname[ip].Data()), Form("Momentum (%s);Momentum (%s) (GeV/c);Counts",pname2[ip].Data(),pname2[ip].Data()), 5000, 0.0, 5.0);
    new TH1F( Form("%s_CosT_CM",pname[ip].Data()), Form("cos#theta (%s);cos#theta(%s);Counts",pname2[ip].Data(),pname2[ip].Data()), 3000, -1.5, 1.5 );
    new TH1F( Form("%s_MMass_CM",pname[ip].Data()), Form("Missing mass (%s);Missing Mass(%s) (GeV/c^{2});Counts",pname2[ip].Data(),pname2[ip].Data()), 5000, 0.0, 5.0);
    new TH1F( Form("%s_MMomentum_CM",pname[ip].Data()), Form("Missing Momentum (%s);Missing Momentum(%s) (GeV/c);Counts",pname2[ip].Data(),pname2[ip].Data()), 5000, 0.0, 5.0);
    new TH1F( Form("%s_MCosT_CM",pname[ip].Data()), Form("Missing cos#theta (%s);Missing cos#theta(%s);Counts",pname2[ip].Data(),pname2[ip].Data()), 3000, -1.5, 1.5);
  }
  for(int ip=0; ip<7; ip++){
    new TH1F( Form("%sn_Mass_CM",pname[ip].Data()), Form("Invariant Mass (%s);IM(%s) (GeV/c^{2});Counts",pname2[ip].Data(),pname2[ip].Data()), 5000, 0.0, 5.0);
    new TH1F( Form("%sn_Momentum_CM",pname[ip].Data()), Form("Momentum (%s);Momentum (%s) (GeV/c);Counts",pname2[ip].Data(),pname2[ip].Data()), 5000, 0.0, 5.0);
    new TH1F( Form("%sn_CosT_CM",pname[ip].Data()), Form("cos#theta (%s);cos#theta(%s);Counts",pname2[ip].Data(),pname2[ip].Data()), 3000, -1.5, 1.5 );
    new TH1F( Form("%sn_MMass_CM",pname[ip].Data()), Form("Missing mass (%s);Missing Mass(%s) (GeV/c^{2});Counts",pname2[ip].Data(),pname2[ip].Data()), 5000, 0.0, 5.0);
    new TH1F( Form("%sn_MMomentum_CM",pname[ip].Data()), Form("Missing Momentum (%s);Missing Momentum(%s) (GeV/c);Counts",pname2[ip].Data(),pname2[ip].Data()), 5000, 0.0, 5.0);
    new TH1F( Form("%sn_MCosT_CM",pname[ip].Data()), Form("Missing cos#theta (%s);Missing cos#theta(%s);Counts",pname2[ip].Data(),pname2[ip].Data()), 3000, -1.5, 1.5);
  }

  /* Hodoscope */
  std::cout << "Define Histograms for Hodoscopes" << std::endl;
  const int nhodo = 2;
  int hodoid[nhodo]={
    CID_NC,CID_T0
  };
  for(int ihodo=0; ihodo<nhodo; ihodo++){
    const int cid = hodoid[ihodo];
    const char* name = dlist -> GetName(cid).c_str();
    const int nsegs = dlist -> GetNsegs(cid);
    for( int seg=1; seg<=nsegs; seg++ ){
      const char* ud[2] = {"u","d"};
      const char* udname[2] = {"up","down"};
      for(int iud=0; iud<2; iud++){
        new TH1F( Form("%s%s%d_TOF",name,ud[iud],seg),   Form("TOF %s-%d(%s);Time (ns);Counts",name,seg,udname[iud]),    1000,    40, 100 );
        new TH1F( Form("%s%s%d_TimeOffset",name,ud[iud],seg),   Form("Time Offset %s-%d(%s);Time (ns);Counts",name,seg,udname[iud]),    1000,    -20, 20 );
        //new TH2F( Form("%s%s%d_dEvsTOF",name,ud[iud],seg),   Form("dE TOF corr. %s-%d(%s);dE (MeV);TOF (ns)",name,seg,udname[iud]),     500,    0, 100,  1000,    -20, 20 );
        //new TH2F( Form("%s%s%d_dEvsTimeOffset",name,ud[iud],seg),   Form("dE Time Offset corr. %s-%d(%s);dE (MeV);Time (ns)",name,seg,udname[iud]),   500,    0, 100,  1000,    -20, 20 );
      }
      //new TH1F( Form("%s%d_TOF",name,seg),   Form("TOF %s-%d;Time (ns);Counts",name,seg),    1000,    40, 100 );
      new TH1F( Form("%s%d_TimeOffset",name,seg),   Form("Time Offset %s-%d;Time (ns);Counts",name,seg),    1000,    -20, 20 );
      //new TH2F( Form("%s%d_dEvsTOF",name,seg),   Form("dEMean TOF corr. %s-%d;dE (MeV);TOF (ns)",name,seg),     500,    0, 100,  1000,    -20, 20 );
      //new TH2F( Form("%s%d_dEvsTimeOffset",name,seg),   Form("dEMean Time Offset corr. %s-%d;dE (MeV);Time (ns)",name,seg),     500,    0, 100,  1000,    -20, 20 );
      new TH1F( Form("%s%d_CTimeSub",name,seg),   Form("Subtract CTime %s-%d;Subtract Time (ns);Counts",name,seg),    1000,    -10, 10 );
      new TH1F( Form("%s%d_Overbeta",name,seg), "1/#beta;1/#beta;Counts", 3000, 0, 5);
    }
    new TH1F( Form("%s_TOF",name),   Form("TOF %s;Time (ns);Counts",name),    1000,    40, 100 );
    new TH1F( Form("%s_TimeOffset",name),   Form("Time Offset %s;Time (ns);Counts",name),    1000,    -20, 20 );
    new TH2F( Form("%s_dEvsTOF",name),   Form("dEMean TOF corr. %s;dE (MeV);TOF (ns)",name),     500,    0, 100,  1000,    -20, 20 );
    new TH2F( Form("%s_dEvsTimeOffset",name),   Form("dEMean Time Offset corr. %s;dE (MeV);Time (ns)",name),     500,    0, 100,  1000,    -20, 20 );
  }

  return true;
}
