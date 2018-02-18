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
  //}
  //else {
  //  return true;
  //}
  
  /* Event selection */
  //int MulCDH=0;
  //for(int i=0; i<cdsMan->nCDH(); i++){
  //  HodoscopeLikeHit* hit = cdsMan->CDH(i);
  //  if(hit->CheckRange()) MulCDH++;
  //}
  //if(MulCDH!=3) return true;
  //if(cdstrackMan->nGoodTrack()!=3) return true;
  
  if(particle->nCDS()==0) return true;
  if(particle->nPiplus()==0&&particle->nPiminus()==0) return true;

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
    double cdhtime,dis;
    int cdhseg;
    if(!cdc->GetCDHHit(cdsMan,cdhseg,cdhtime)) return 0;
    TVector3 vtx1,vtx2;
    double par[5];
    cdc->GetParameters(CID_CDC,par,vtx1);
    double mom=cdc->Momentum();
    TrackTools::CalcLineHelixVertex(beam,cdc,vtx1,vtx2,dis);

    pCDS* cdstrack=new pCDS();
    cdstrack->SetParameters(par);
    cdstrack->SetTOF(cdhtime);
    cdstrack->SetCDHSeg(cdhseg);

    double time_beam=(beam->t0pos()-cdstrack->vbeam()).Mag()/beam->beta()/(Const*100);
    TVector3 cdhvtx=cdc->CDHVertex();
    double cdc_dis=MathTools::CalcHelixArc(par,cdhvtx,vtx2);
    double beta=cdc_dis/(cdhtime-beam->t0time()-time_beam)/(Const*100);
    double mass2=mom*mom*(1/(beta*beta)-1);

    double beta_calc,tofvtxcdc;
    double tof=cdhtime-beam->t0time();
    TrackTools::FindMass2(cdc,beam,tof,beta_calc,mass2,tofvtxcdc);

    FillHist("CDS_Mass2vsMomentum",mass2,mom);
  }



  // ############### //
  // Vertex decision //
  // ############### //
  for(int it=0; it<particle->nCDS(); it++){
    TVector3 vertex;
    double vdis = 9999;
    pCDS* pro = particle->cdsi(it);
    if(GeomTools::GetID(vertex)!=CID_Fiducial) continue;

    if(pro->pid()!=CDS_PiPlus&&pro->pid()!=CDS_PiMinus) continue;
    //if(pro->pid()!=CDS_Proton) continue;
    
    if(0.25<TMath::Abs(pro->mom())) continue;

    CDSTrack* track = cdstrackMan->GoodTrack(pro->id());
    if(track->nCDHHit()!=1) continue;
    HodoscopeLikeHit* cdh = track->CDHHit(cdsMan,0);
    int seg = cdh->seg();

    /* Fill data */
    int T0Segment = beam->t0seg();
    double T0Timing = beam->t0time();
    const int nhodo = 2;
    int hodoid[nhodo]={
      CID_CDH,CID_IH
    };
    for(int ihodo=0; ihodo<nhodo; ihodo++){
      const int cid = hodoid[ihodo];
      const char* name = dlist -> GetName(cid).c_str();
      const int nsegs = dlist -> GetNsegs(cid);
      int nhit=0;
      for( int i=0; i<cdsMan->nHodo(cid); i++ ){
        HodoscopeLikeHit *hit = cdsMan->Hodoi(cid,i);
        int seg = hit->seg();
        int adc[2]; adc[0]=hit->adcu(); adc[1]=hit->adcd();
        int tdc[2]; tdc[0]=hit->tdcu(); tdc[1]=hit->tdcd();
        double time[2]; time[0]=hit->time(0), time[1] = hit->time(1);
        double ctime[2]; ctime[0]=hit->ctime(0); ctime[1]=hit->ctime(1);
        double ene[2]; ene[0]=hit->ene(0); ene[1]=hit->ene(1);
        double hitpos = hit->hitpos();
        double emean = hit->emean();
        double tmean = hit->tmean();
        double ctmean = hit->ctmean();
        double tsub = hit->tsub();
        double ctsub = hit->ctsub();
        if(cid==CID_IH){
          emean = ene[0];
          tmean = time[0];
          ctmean = ctime[0];
          tsub = time[0];
          ctsub = time[0];
        }
        double tof = ctmean - T0Timing;

        const char* ud[2] = {"u","d"};
        for(int iud=0; iud<2; iud++){
          FillHist(Form("%s%s%d_ADC",name,ud[iud],seg),adc[iud],1);
          if(tdc[iud]>0){
            FillHist(Form("%s%s%d_TDC",name,ud[iud],seg),tdc[iud],1);
            FillHist(Form("%s%s%d_Time",name,ud[iud],seg),time[iud],1);
            FillHist(Form("%s%s%d_CTime",name,ud[iud],seg),ctime[iud],1);
          }
        }
        if( hit->CheckRange() ){
          for(int iud=0; iud<2; iud++){
            FillHist(Form("%s%s%d_ADCwithTDC",name,ud[iud],seg),adc[iud],1);
            FillHist(Form("%s%s%d_dE",name,ud[iud],seg),ene[iud],1);
            FillHist(Form("%s%s%d_ADCvsTDC",name,ud[iud],seg),adc[iud],tdc[iud],1);
            FillHist(Form("%s%s%d_dEvsTime",name,ud[iud],seg),ene[iud],time[iud],1);
            FillHist(Form("%s%s%d_dEvsCTime",name,ud[iud],seg),ene[iud],ctime[iud],1);
            FillHist(Form("%s%s%d_dEMeanvsCTime",name,ud[iud],seg),emean,ctime[iud],1);
          }
          FillHist(Form("%s%d_dEMean",name,seg),emean,1);
          FillHist(Form("%s%d_TimeMean",name,seg),tmean,1);
          FillHist(Form("%s%d_TimeSub",name,seg),tsub,1);
          FillHist(Form("%s%d_CTimeMean",name,seg),ctmean,1);
          FillHist(Form("%s_CTimeMean",name),ctmean,1);
          FillHist(Form("%s%d_EventvsCTimeMean",name,seg),header->ev(),ctmean,1);
          FillHist(Form("%s_EventvsCTimeMean",name),header->ev(),ctmean,1);
          FillHist(Form("%s%d_CTimeSub",name,seg),ctsub,1);
          FillHist(Form("%s%d_TOF",name,seg),tof,1);
          FillHist(Form("%s_TOF",name),tof,1);
          FillHist(Form("%s_TOFvsT0Segment",name),tof,T0Segment,1);
          FillHist(Form("%s_TOFvsCTimeMean",name),tof,ctmean,1);
          FillHist(Form("%s_TOFvsT0CTimeMean",name),tof,T0Timing,1);
          FillHist(Form("%s%d_Position",name,seg),hitpos,1);
          FillHist(Form("%s%d_dEMeanvsCTimeMean",name,seg),emean,ctmean,1);
          FillHist(Form("%s_HitPattern",name),seg,1);
          nhit++;
          for( int j=0; j<20; j++ ){
            int val = header->pattern(j);
            if( 0<val ){
              FillHist(Form("%s_TriggerPatternvsTOF",name),j,tof,1);
            }
          }
        }
      }
      FillHist(Form("%s_Multiplicity",name),nhit,1);
    }

    /* Time Sub Offset */
    TVector3 vcdhtrack = track->CDHVertex();
    double lv = 9999.9;
    if(!conf->GetGeomMapManager()->GetLightVelocity(CID_CDH,seg,lv)) continue;
    double calcsubtime = (vcdhtrack.Z()-79.0/2.0)/lv;
    double suboffs = cdh->ctsub() - calcsubtime;
    FillHist(Form("CDH%d_CTimeSubOffset",seg),suboffs,1);
    /* Time Offset and Slewing */
    if(vcdhtrack.Z()<-5.0||vcdhtrack.Z()<5.0) continue;
    double dt = pro->dt();
    double tof = cdh->ctmean() - beam->t0time();
    double utof = cdh->ctime(0) - beam->t0time();
    double dtof = cdh->ctime(1) - beam->t0time();
    double offset = tof - dt;
    double uoffset = utof - dt;
    double doffset = dtof - dt;
    FillHist(Form("CDH%d_TimeOffset",seg),offset,1);
    FillHist(Form("CDH%d_dEMeanvsTimeOffset",seg),cdh->emean(),offset,1);
    FillHist(Form("CDHu%d_TimeOffset",seg),uoffset,1);
    FillHist(Form("CDHu%d_dEvsTimeOffset",seg),cdh->eu(),offset,1);
    FillHist(Form("CDHu%d_dEMeanvsTimeOffset",seg),cdh->emean(),offset,1);
    FillHist(Form("CDHd%d_TimeOffset",seg),doffset,1);
    FillHist(Form("CDHd%d_dEvsTimeOffset",seg),cdh->ed(),offset,1);
    FillHist(Form("CDHd%d_dEMeanvsTimeOffset",seg),cdh->emean(),offset,1);
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

  new TH2F( "CDS_Mass2vsMomentum", "Mass^{2} vs. Momentum;Mass^{2} (GeV/c^{2})^{2};Momentum (GeV/c)", 2000, -5, 5, 4000, -2, 18 );


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
        new TH1F( Form("%s%s%d_ADC",name,ud[iud],seg),   Form("ADC %s-%d(%s);ADC ch.;Counts",name,seg,udname[iud]),    2000,    0, 4000 );
        new TH1F( Form("%s%s%d_ADCwithTDC",name,ud[iud],seg), Form("ADC with TDC %s-%d(%s);ADC ch.;Counts",name,seg,udname[iud]),  2000,    0, 4000 );
        new TH1F( Form("%s%s%d_dE",name,ud[iud],seg),   Form("dE %s-%d(%s);dE (MeV);Counts",name,seg,udname[iud]),    600,    0, 60 );
        new TH1F( Form("%s%s%d_TDC",name,ud[iud],seg),   Form("TDC %s-%d(%s);TDC ch.;Counts",name,seg,udname[iud]),    2000,    0, 4000 );
        new TH1F( Form("%s%s%d_Time",name,ud[iud],seg),   Form("Time %s-%d(%s);Time (ns);Counts",name,seg,udname[iud]),    2000,    -100, 100 );
        new TH1F( Form("%s%s%d_CTime",name,ud[iud],seg),   Form("CTime %s-%d(%s);Time (ns);Counts",name,seg,udname[iud]),    2000,    -100, 100 );
        new TH2F( Form("%s%s%d_ADCvsTDC",name,ud[iud],seg),   Form("ADC TDC corr. %s-%d(%s);ADC ch.;TDC ch.",name,seg,udname[iud]),     400,    0, 4000,  400,    0, 4000 );
        new TH2F( Form("%s%s%d_dEvsTime",name,ud[iud],seg),   Form("dE Time corr. %s-%d(%s);dE (MeV);Time (ns)",name,seg,udname[iud]),     400,    0, 60,  2000,    -100, 100 );
        new TH2F( Form("%s%s%d_dEvsCTime",name,ud[iud],seg),   Form("dE CTime corr. %s-%d(%s);dE (MeV);Time (ns)",name,seg,udname[iud]),     200,    0, 20,  2000,    -100, 100 );
        new TH2F( Form("%s%s%d_dEMeanvsCTime",name,ud[iud],seg),   Form("dEMean CTime corr. %s-%d(%s);dE (MeV);Time (ns)",name,seg,udname[iud]),     400,    0, 60,  2000,    -100, 100 );
        new TH1F( Form("%s%s%d_TimeOffset",name,ud[iud],seg),   Form("Time Offset %s-%d(%s);Time (ns);Counts",name,seg,udname[iud]),    2000,    -100, 100 );
        new TH2F( Form("%s%s%d_dEvsTimeOffset",name,ud[iud],seg),   Form("dE Time Offset corr. %s-%d(%s);dE (MeV);Time (ns)",name,seg,udname[iud]),     400,    0, 60,  2000,    -100, 100 );
        new TH2F( Form("%s%s%d_dEMeanvsTimeOffset",name,ud[iud],seg),   Form("dEMean Time Offset corr. %s-%d(%s);dE (MeV);Time (ns)",name,seg,udname[iud]),     400,    0, 60,  2000,    -100, 100 );
      }
      new TH1F( Form("%s%d_dEMean",name,seg),   Form("Mean dE %s-%d;dE (MeV);Counts",name,seg),    600,    0, 60 );
      new TH1F( Form("%s%d_TimeMean",name,seg),   Form("Mean Time %s-%d;Time (ns);Counts",name,seg),    2000,    -100, 100 );
      new TH1F( Form("%s%d_TimeSub",name,seg),   Form("Subtract Time %s-%d;Subtract Time (ns);Counts",name,seg),    2000,    -100, 100 );
      new TH1F( Form("%s%d_CTimeMean",name,seg),   Form("Mean CTime %s-%d;Time (ns);Counts",name,seg),    2000,    -100, 100 );
      new TH1F( Form("%s%d_CTimeSub",name,seg),   Form("Subtract CTime %s-%d;Subtract Time (ns);Counts",name,seg),    2000,    -100, 100 );
      new TH1F( Form("%s%d_TOF",name,seg),   Form("TOF between T0 and %s-%d;Time of flight T0-%s%d (ns);Counts",name,seg,name,seg),    3000,    -50, 100 );
      new TH1F( Form("%s%d_Position",name,seg),   Form("Hit position %s-%d;Position (cm);Counts",name,seg),    1000,    -50, 50 );
      new TH2F( Form("%s%d_dEMeanvsCTimeMean",name,seg),   Form("dEMean CTimeMean corr. %s-%d;dE (MeV);Time (ns)",name,seg),     400,    0, 60,  2000,    -100, 100 );
      new TH2F( Form("%s%d_EventvsCTimeMean",name,seg),   Form("Event number vs. Mean CTime %s-%d;Event number;Time (ns);",name,seg), 800, 0, 800000,   2000,    -100, 100 );
      new TH1F( Form("%s%d_TimeOffset",name,seg),   Form("Time Offset %s-%d;Time (ns);Counts",name,seg),    2000,    -100, 100 );
      new TH2F( Form("%s%d_dEMeanvsTimeOffset",name,seg),   Form("dEMean Time Offset corr. %s-%d;dE (MeV);Time (ns)",name,seg),     400,    0, 60,  2000,    -100, 100 );
      new TH1F( Form("%s%d_CTimeSubOffset",name,seg),   Form("Offset of Subtract CTime %s-%d;Subtract Time (ns);Counts",name,seg),    2000,    -100, 100 );
    }
    new TH1F( Form("%s_HitPattern",name), Form("Hit pattern %s;%s Segment;Counts",name,name), nsegs, 0.5, nsegs+0.5);
    new TH1F( Form("%s_Multiplicity",name), Form("Multiplicity %s;Multiplicity;Counts",name), nsegs+1, -0.5, nsegs+0.5 );
    new TH1F( Form("%s_CTimeMean",name),   Form("Mean CTime %s;Time (ns);Counts",name),    2000,    -100, 100 );
    new TH2F( Form("%s_EventvsCTimeMean",name),   Form("Event number vs. Mean CTime %s;Event number;Time (ns);",name), 800, 0, 800000,   2000,    -100, 100 );
    new TH1F( Form("%s_TOF",name),   Form("TOF between T0 and %s;Time of flight T0-%s (ns);Counts",name,name),    3000,    -50, 100 );
    new TH2F( Form("%s_TOFvsT0Segment",name),   Form("TOF between T0 and %s vs. T0 segment;Time of flight T0-%s (ns);T0 segment",name,name),    3000,    -50, 100, 5, 0.5, 5.5 );
    new TH2F( Form("%s_TOFvsCTimeMean",name),   Form("TOF between T0 and %s vs. CTimeMean;Time of flight T0-%s (ns);Time (ns)",name,name),    3000,    -50, 100, 2000, -100, 100 );
    new TH2F( Form("%s_TOFvsT0CTimeMean",name),   Form("TOF between T0 and %s vs. T0 CTimeMean;Time of flight T0-%s (ns);T0 Time (ns)",name,name),    3000,    -50, 100, 200, -10, 10 );
    TH2F* h2 = new TH2F( Form("%s_TriggerPatternvsTOF",name),   Form("TriggerPattern vs. TOF between T0 and %s;;Time of flight T0-%s (ns)",name,name), 18, -0.5, 17.5 ,    3000,    -50, 100);
    h2->GetXaxis()->SetBinLabel(1,"All");
    h2->GetXaxis()->SetBinLabel(2,"Beam");
    h2->GetXaxis()->SetBinLabel(3,"Kaon");
    h2->GetXaxis()->SetBinLabel(4,"KCDH1f");
    h2->GetXaxis()->SetBinLabel(5,"Pion");
    h2->GetXaxis()->SetBinLabel(6,"Proton");
    h2->GetXaxis()->SetBinLabel(7,"KCDH1");
    h2->GetXaxis()->SetBinLabel(8,"KCDH2");
    h2->GetXaxis()->SetBinLabel(9,"PivBVC");
    h2->GetXaxis()->SetBinLabel(10,"PiCDH1");
    h2->GetXaxis()->SetBinLabel(11,"PiCDH2");
    h2->GetXaxis()->SetBinLabel(12,"Kf");
    h2->GetXaxis()->SetBinLabel(13,"1stMix");
    h2->GetXaxis()->SetBinLabel(14,"Charged");
    h2->GetXaxis()->SetBinLabel(15,"Neutral");
    h2->GetXaxis()->SetBinLabel(16,"Cosmic");
    h2->GetXaxis()->SetBinLabel(17,"Reject");
    h2->GetXaxis()->SetBinLabel(18,"SIM");
  }



  return true;
}
