// MyAnalysisCDSCalib.cpp

#include "MyAnalysisCDSCalib.h"

MyAnalysisCDSCalib::MyAnalysisCDSCalib(TFile* rt, ConfMan* conf)
{
  Initialize(conf);
  CutCondition();
  Clear();
}

MyAnalysisCDSCalib::~MyAnalysisCDSCalib()
{
  Clear();
  rtFile->cd();
  rtFile->Write();
  rtFile->Close();
}

void MyAnalysisCDSCalib::Clear()
{
}

bool MyAnalysisCDSCalib::DoAnalysis(ConfMan* conf, EventHeader* header, BeamLineHitMan* blMan, BeamLineTrackMan* bltrackMan, CDSHitMan* cdsMan, CDSTrackingMan* cdstrackMan, Particle* particle)
{
  rtFile->cd();

  DetectorList *dlist=DetectorList::GetInstance();
  if(particle->nBeam()!=1) return false;
  pBeam* beam = particle->beam(0);

  /* Trigger selection */
  //if(!header->trigmode(Mode_KCDH2)) return true;
  //if(header->IsTrig(Trig_KCDH2)&&header->IsTrig(Trig_1stMix)){
  //if(header->IsTrig(Trig_KCDH2)&&header->IsTrig(Trig_1stMix)){
  //}
  //else {
  //  return true;
  //}
  if(header->IsTrig(Trig_Cosmic)) return true;
  
  /* Event selection */
  int MulCDH=0;
  for(int i=0; i<cdsMan->nCDH(); i++){
    HodoscopeLikeHit* hit = cdsMan->CDH(i);
    if(hit->CheckRange()) MulCDH++;
  }
  if(MulCDH==0) return true;
  if(cdstrackMan->nGoodTrack()==0) return true;
  
  //if(particle->nCDS()==0) return true;
  //if(particle->nPiplus()==0&&particle->nPiminus()==0) return true;

  // Trigger Pattern
  FillHist("TriggerPattern",0);
  for( int i=0; i<20; i++ ){
    int val = header->pattern(i);
    if( 0<val ){
      FillHist("TriggerPattern",i);
      FillHist(Form("Trigger%02d_TDC",i),val);
    }
  }
  
  for(int id=0; id<cdstrackMan->nGoodTrack(); id++){
    CDSTrack* cdc=cdstrackMan->GoodTrack(id);
    double cdhtime=0,dis=0;
    int cdhseg=0;
    if(!cdc->GetCDHHit(cdsMan,cdhseg,cdhtime)) continue;
    if(cdc->Chi()>30) continue;

    TVector3 vtx1,vtx2;
    double par[5];
    cdc->GetParameters(CID_CDC,par,vtx1);
    double mom=cdc->Momentum();
    TrackTools::CalcLineHelixVertex(beam,cdc,vtx1,vtx2,dis);

    //if(GeomTools::GetID(vtx1)!=CID_Fiducial) continue;

    double time_beam=beam->CalcVertexTime(vtx1)-beam->t0time();
    TVector3 cdhvtx=cdc->CDHVertex();
    double cdc_dis=MathTools::CalcHelixArc(par,cdhvtx,vtx2);
    double beta=cdc_dis/(cdhtime-beam->t0time()-time_beam)/(Const*100);
    double mass2=mom*mom*(1/(beta*beta)-1);
    FillHist("CDS_Mass2vsMomentum",mass2,mom);
    double beta_calc,tofvtxcdc;
    double tof=cdhtime-beam->t0time();
    TrackTools::FindMass2(cdc,beam,tof,beta_calc,mass2,tofvtxcdc);
    FillHist("CDS_Mass2vsMomentum_AfterFindMass",mass2,mom);
    int pid=TrackTools::PID(mom,mass2);

    cdc->SetPID(pid);

    double tmpl;
    cdc->CalcVertexTimeLength(beam->t0pos(),beam->bpcdir(),cdsMass[pid],vtx2,vtx1,tofvtxcdc,tmpl,true);

    TVector3 Pd;
    if(!cdc->GetMomentum(vtx1,Pd,true,true)) continue;
    mom = mom/TMath::Abs(mom)*Pd.Mag();
    FillHist("CDS_Mass2vsMomentum_AfterELossCorr",mass2,mom);

    //if(pid!=CDS_PiPlus&&pid!=CDS_PiMinus) continue;
    if(pid==CDS_Other) continue;
    //if(TMath::Abs(mom)<0.2) continue;
    //if(cdhvtx.Z()<-5.0||5.0<cdhvtx.Z()) continue;

    /* Time Offset and Slewing */
    HodoscopeLikeHit* cdh = cdc->CDHHit(cdsMan,0);
    TVector3 posl;
    cdc->GetParameters(CID_CDCCFRP,par,posl);
    cdc->GetParameters(CID_CDC,par,vtx1);
    cdc_dis=MathTools::CalcHelixArc(par,cdhvtx,posl);
    double beta_c = TMath::Abs(mom)/TMath::Sqrt(mom*mom+cdsMass[pid]*cdsMass[pid]);
    double tof_c = cdc_dis/beta_c/(Const*100);
    double offset = (cdh->ctmean()-beam->t0time()) - (tof_c+tofvtxcdc+time_beam);

    FillHist("CDS_OverbetavsMomentum_Calc",1./beta_c,mom);
    FillHist("CDS_TOFvsMomentum_Calc",tof_c,mom);
    FillHist("CDS_OverbetavsMomentum",1./beta,mom);
    FillHist("CDS_TOFvsMomentum",(cdh->ctmean()-beam->t0time())-tofvtxcdc-time_beam,mom);
    FillHist("CDS_TimeBeam",time_beam);
    FillHist("CDS_TOFVTXCDC",tofvtxcdc);

    FillHist(Form("CDH_TimeOffset"),offset,1);
    FillHist(Form("CDH_dEMeanvsTimeOffset"),cdh->emean(),offset,1);
    FillHist(Form("CDH%d_TimeOffset",cdh->seg()),offset,1);
    FillHist(Form("CDH%d_dEMeanvsTimeOffset",cdh->seg()),cdh->emean(),offset,1);
    FillHist(Form("CDHu%d_TimeOffset",cdh->seg()),offset,1);
    FillHist(Form("CDHu%d_dEvsTimeOffset",cdh->seg()),cdh->eu(),offset,1);
    FillHist(Form("CDHu%d_dEMeanvsTimeOffset",cdh->seg()),cdh->emean(),offset,1);
    FillHist(Form("CDHd%d_TimeOffset",cdh->seg()),offset,1);
    FillHist(Form("CDHd%d_dEvsTimeOffset",cdh->seg()),cdh->ed(),offset,1);
    FillHist(Form("CDHd%d_dEMeanvsTimeOffset",cdh->seg()),cdh->emean(),offset,1);

  }

  return true;
}

bool MyAnalysisCDSCalib::FillHist(TString name, double val1, int weight)
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

bool MyAnalysisCDSCalib::FillHist(TString name, TString val1, int weight)
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

bool MyAnalysisCDSCalib::FillHist(TString name, double val1, double val2, int weight)
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

bool MyAnalysisCDSCalib::FillHist(TString name, TString val1, TString val2, int weight)
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

void MyAnalysisCDSCalib::CutCondition()
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

bool MyAnalysisCDSCalib::Initialize(ConfMan* confMan)
{
  std::cout << "### MyAnalysisCDSCalib::Initialize ###" << std::endl;

  DetectorList *dlist=DetectorList::GetInstance();
  std::string ofname = confMan->GetOutFileName();
  ofname.insert(ofname.find(".root"),"_anaCDSCalib");

  rtFile =  new TFile( Form("%s",ofname.c_str()), "RECREATE");
  rtFile -> cd();

  // CDS Particles
  std::cout << "Define Histograms for CDS particles" << std::endl;
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
  for(int i=1; i<=18; i++){
    new TH1F( Form("Trigger%02d_TDC",i), Form("Trigger%02d TDC;TDC ch.;Counts",i), 4096, 0, 4096);
  }

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

  new TH2F( "CDS_Mass2vsMomentum", "Mass^{2} vs. Momentum;Mass^{2} (GeV/c^{2})^{2};Momentum (GeV/c)", 3000, -5, 10, 800, -2, 2 );
  new TH2F( "CDS_Mass2vsMomentum_AfterFindMass", "Mass^{2} vs. Momentum;Mass^{2} (GeV/c^{2})^{2};Momentum (GeV/c)", 3000, -5, 10, 800, -2, 2 );
  new TH2F( "CDS_Mass2vsMomentum_AfterELossCorr", "Mass^{2} vs. Momentum;Mass^{2} (GeV/c^{2})^{2};Momentum (GeV/c)", 3000, -5, 10, 800, -2, 2 );
  new TH2F( "CDS_OverbetavsMomentum", "1/#beta vs. Momentum;1/#beta;Momentum (GeV/c)", 1000, 0, 10, 800, -2, 2 );
  new TH2F( "CDS_TOFvsMomentum", "TOF vs. Momentum;TOF (ns);Momentum (GeV/c)", 1000, 0, 40, 800, -2, 2 );
  new TH2F( "CDS_OverbetavsMomentum_Calc", "1/#beta vs. Momentum;1/#beta;Momentum (GeV/c)", 1000, 0, 10, 800, -2, 2 );
  new TH2F( "CDS_TOFvsMomentum_Calc", "TOF vs. Momentum;TOF (ns);Momentum (GeV/c)", 1000, 0, 40, 800, -2, 2 );

  new TH1F( "CDS_TimeBeam", "Time Beam;Time Beam (ns);Counts", 1000, 0, 10);
  new TH1F( "CDS_TOFVTXCDC", "TOFVTXCDC;TOFVTXCDC (ns);Counts", 1000, 0, 10);

  /* Hodoscope */
  std::cout << "Define Histograms for Hodoscopes" << std::endl;
  const int nhodo = 2;
  int hodoid[nhodo]={
    CID_CDH,CID_IH
  };
  for(int ihodo=0; ihodo<1; ihodo++){
    const int cid = hodoid[ihodo];
    const char* name = dlist -> GetName(cid).c_str();
    const int nsegs = dlist -> GetNsegs(cid);
    for( int seg=1; seg<=nsegs; seg++ ){
      const char* ud[2] = {"u","d"};
      const char* udname[2] = {"up","down"};
      for(int iud=0; iud<2; iud++){
        new TH1F( Form("%s%s%d_TimeOffset",name,ud[iud],seg),   Form("Time Offset %s-%d(%s);Time (ns);Counts",name,seg,udname[iud]),    400,    -10, 10 );
        new TH2F( Form("%s%s%d_dEvsTimeOffset",name,ud[iud],seg),   Form("dE Time Offset corr. %s-%d(%s);dE (MeV);Time (ns)",name,seg,udname[iud]),     1000,    0, 100,  400,    -10, 10 );
        new TH2F( Form("%s%s%d_dEMeanvsTimeOffset",name,ud[iud],seg),   Form("dEMean Time Offset corr. %s-%d(%s);dE (MeV);Time (ns)",name,seg,udname[iud]),     1000,    0, 100,  400,    -10, 10 );
      }
      new TH1F( Form("%s%d_TimeOffset",name,seg),   Form("Time Offset %s-%d;Time (ns);Counts",name,seg),    400,    -10, 10 );
      new TH2F( Form("%s%d_dEMeanvsTimeOffset",name,seg),   Form("dEMean Time Offset corr. %s-%d;dE (MeV);Time (ns)",name,seg),     1000,    0, 100,  400,    -10, 10 );
    }

  }



  return true;
}
