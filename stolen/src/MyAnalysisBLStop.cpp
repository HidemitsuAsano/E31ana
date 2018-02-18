// MyAnalysisBL.cpp

#include "MyAnalysisBLStop.h"

MyAnalysisBL::MyAnalysisBL(TFile* rtFile, ConfMan* conf, bool simflag)
{
  D5 = new BeamSpectrometer(conf);
  Initialize(rtFile, conf);
  SIMULATION = simflag;
  Clear();
}

MyAnalysisBL::~MyAnalysisBL()
{
  if(D5) delete D5;
}

bool MyAnalysisBL::Clear()
{
  T0Timing = -9999.0;
  T0Segment = -1;
  BHDSegment = -1;
  BeamPID = Beam_Other;
  BeamMom = -9999.0; 
  BeamChi2 = -9999.0;
  BeamTOF = -9999.0;
  trackblc1 = 0;
  trackblc2 = 0;
  trackbpc  = 0;
  D5 -> Clear();
  FIDUCIAL = false;

  return true;
}

bool MyAnalysisBL::DoAnalysis(ConfMan* conf, EventHeader* header, BeamLineHitMan* blMan, BeamLineTrackMan* bltrackMan, Particle* particle, double _mom)
{
  int icount = 0;
  DetectorList *dlist=DetectorList::GetInstance();
  double E0Timing = -9999.0;
  double E0Energy = -9999.0;
  double E0Segment = -1;
  int MulE0 = 0;
  double DEFTiming = -9999.0;
  double DEFEnergy = -9999.0;
  double DEFSegment = -1;
  int MulDEF = 0;

  // ############# //
  // Trigger check //
  // ############# //
  /* ## Select K/f trigger ## */
  //if(!header->IsTrig(Trig_Kaon)) return false;
  //if(!header->IsTrig(Trig_KCDH1)&&!header->IsTrig(Trig_PiCDH1)) return false;
  //if(!header->IsTrig(Trig_KCDH1)) return false;
  //if(!header->IsTrig(Trig_Kf)) return false;
  //if(!header->IsTrig(Trig_Kf)||!header->IsTrig(Trig_1stMix)) return false;

  /* ## All events ## */
  FillHist("Beam_Count",icount); icount++;

  // ############# //
  // T0 timing set //
  // ############# //
  {
    T0Timing = -9999.; 
    int MulT0 = 0;
    for(int i=0; i<blMan->nT0(); i++){
      if(blMan->T0(i)->CheckRange()){
        MulT0++;
        T0Timing = blMan->T0(i)->ctmean();
        T0Segment = blMan->T0(i)->seg();
      }
    }
    FillHist("T0_Mult",MulT0);
    /* ## T0 1hit request ## */
    if(MulT0!=1) return false;
    FillHist("Beam_Count",icount); icount++;
    FillHist("T0_MultwCut",MulT0);
  }

  // ###################### //
  // TOF between BHD and T0 //
  // ###################### //
  {
    int MulBHD = 0;
    int seg[20]; for(int i=0; i<20; i++) {seg[i] = -1;}
    double tof[20]; for(int i=0; i<20; i++) {tof[i] = -9999.0;}
    for(int i=0; i<blMan->nBHD(); i++){
      if(blMan->BHD(i)->CheckRange()){
        if(blMan->BHD(i)->seg()<6||15<blMan->BHD(i)->seg()) continue;
        MulBHD++;
        seg[MulBHD-1] = blMan->BHD(i)->seg();
        tof[MulBHD-1] = T0Timing - blMan->BHD(i)->ctmean();;
        FillHist("T0BHD_TOF",seg[MulBHD-1],tof[MulBHD-1]);
        FillHist(Form("Beam_TOF"), tof[MulBHD-1]);
      }
    }
    FillHist("BHD_Mult",MulBHD);
    /* ## TOF kaon request ## */
    double tofk[2] = {3.0000,4.5000};
    //double tofk[2] = {27.9637,30.4648};
    double tofpi[2] = {25.0, 27.0};
    bool TOFflag=false;
    for(int i=0; i<MulBHD; i++){
      if(tofk[0]<tof[i]&&tof[i]<tofk[1]) {
        TOFflag=true;
        BeamTOF = tof[i];
        BHDSegment = seg[i];
        TVector3 gpos;
        conf->GetGeomMapManager()->GetGPos(CID_BHD,BHDSegment,gpos);
        BHDX = gpos.X();
        break;
      }
    }
    if(!TOFflag) return false;
    FillHist("Beam_Count",icount); icount++;
    FillHist("BHD_MultwCut",MulBHD);
    FillHist(Form("Beam_TOFwCut"), BeamTOF);
  }
  // #################### //
  // DC tracking          //
  // #################### //
  {
    bltrackMan->DoTracking(blMan,conf,true,true);
    const int ndc = 3;
    int dccid[ndc] = {CID_BLC1,CID_BLC2,CID_SDC};
    int nhitdc[ndc] = {0,0,0};
    int trackid[ndc] = {-1,-1,-1};
    double ll[ndc] = {-30.0,-30.0,-30.0};
    double ul[ndc] = {100.0,100.0,100.0};
    double trigll[ndc] = {-10.0,-10.0,-10.0};
    double trigul[ndc] = { 10.0, 10.0, 10.0};
    double chiul[ndc] = {10.0, 10.0, 20.0};
    int idc;
    // ==== BLC1 and BLC2 analysis ==== // 
    for(idc=0;idc<2;idc++){
      const char* dcname = dlist->GetName(dccid[idc]).data();
      FillHist(Form("%s_nTrack",dcname),bltrackMan->ntrackBLDC(dccid[idc]));
      for(int i=0; i<bltrackMan->ntrackBLDC(dccid[idc]); i++){
        if(i>=1) break;
        LocalTrack* track = bltrackMan->trackBLDC(dccid[idc],i);
        double tracktime = track -> GetTrackTime();
        double chi2 = track -> chi2all();
        FillHist(Form("%s_TrackTime",dcname),tracktime);
        if(ll[idc]<tracktime&&tracktime<ul[idc]){
          FillHist(Form("%s_Chi2",dcname),chi2);
          if(chi2<chiul[idc]){
            nhitdc[idc]++;
            trackid[idc] = i;
            FillHist(Form("%s_TrackTimewCut",dcname),tracktime);
            FillHist(Form("%s_Chi2wCut",dcname),chi2);
          }
        }
      }
      FillHist(Form("%s_nGoodTrack",dcname),nhitdc[idc]);
    }
    bool BLCGoodflag = false;
    if(nhitdc[0]==1&&nhitdc[1]==1){
      LocalTrack* blc1 = bltrackMan->trackBLC1(trackid[0]);
      LocalTrack* blc2 = bltrackMan->trackBLC2(trackid[1]);
      double blc1time = blc1->GetTrackTime();
      double blc2time = blc2->GetTrackTime();
      if(trigll[0]<blc1time&&blc1time<trigul[0]
          && trigll[1]<blc2time&&blc2time<trigul[1])
        BLCGoodflag = true;
    }
    /* ## BLC1/2 1 good track request ## */
    if(!BLCGoodflag) return false;
    FillHist("Beam_Count",icount); icount++;
    trackblc1 = bltrackMan->trackBLC1(trackid[0]);
    trackblc2 = bltrackMan->trackBLC2(trackid[1]);
    // ==== Beam momentum analysis ==== // 
    D5->TMinuitFit(trackblc1,trackblc2,conf);
    BeamChi2 = D5->chisquare();
    BeamMom = D5->mom();
    FillHist(Form("Beam_Momentum"), BeamMom);
    FillHist(Form("Beam_MomentumvsTOF"), BeamMom,BeamTOF);
    FillHist(Form("Beam_Chi2"), BeamChi2);
    FillHist(Form("BHD_SegmentvsBeamMomentum"),BHDSegment,BeamMom); 
    /* ## Beam chisquare < 20 request ## */
    if(!(BeamChi2<20)) return false;
    FillHist("Beam_Count",icount); icount++;
    FillHist(Form("Beam_MomentumwCut"), BeamMom);
    FillHist(Form("Beam_MomentumvsTOFwCut"), BeamMom,BeamTOF);
    FillHist(Form("Beam_Chi2wCut"), BeamChi2);
    FillHist(Form("BHD_SegmentvsBeamMomentumwCut"),BHDSegment,BeamMom); 

    // ==== BVC analysis ==== // 
    double BVCTiming = -9999.0;
    double BVCEnergy = -9999.0;
    double BVCSegment = -1;
    int MulBVC = 0;
    for(int i=0; i<blMan->nBVC(); i++){
      if(blMan->BVC(i)->CheckRange()){
        MulBVC++;
        BVCTiming = blMan->BVC(i)->ctmean();
        BVCEnergy = blMan->BVC(i)->emean();
        BVCSegment = blMan->BVC(i)->seg();
      }
    }
    double tof_bvc = BVCTiming - T0Timing;
    FillHist("BVC_Mult",MulBVC);
    FillHist("BVC_TOF",tof_bvc);
    FillHist("BVC_dE",BVCEnergy);
    FillHist("BVC_dEvsTOF",tof_bvc,BVCEnergy);
    /* ## BVC veto request ## */
    if(MulBVC!=0) return false;
    FillHist("Beam_Count",icount); icount++;

    // ==== E0 analysis ==== // 
    for(int i=0; i<blMan->nE0(); i++){
      if(blMan->E0(i)->CheckRange()){
        MulE0++;
        E0Timing = blMan->E0(i)->ctmean();
        E0Energy = blMan->E0(i)->emean();
        E0Segment = blMan->E0(i)->seg();
      }
    }
    FillHist("E0_Mult",MulE0);
    /* ## E0 1hit request ## */
    if(MulE0!=1) return false;
    FillHist("Beam_Count",icount); icount++;
    double tof_e0 = E0Timing - T0Timing;
    FillHist("E0_TOF",tof_e0);
    FillHist("E0_dE",E0Energy);
    FillHist("E0_dEvsTOF",tof_e0,E0Energy);

    // ==== DEF analysis ==== // 
    for(int i=0; i<blMan->nDEF(); i++){
      if(blMan->DEF(i)->CheckRange()){
        MulDEF++;
        DEFTiming = blMan->DEF(i)->ctmean();
        DEFEnergy = blMan->DEF(i)->emean();
        DEFSegment = blMan->DEF(i)->seg();
      }
    }
    FillHist("DEF_Mult",MulDEF);
    /* ## DEF 1hit request ## */
    if(MulDEF!=1) return false;
    FillHist("Beam_Count",icount); icount++;
    double tof_def = DEFTiming - T0Timing;
    FillHist("DEF_TOF",tof_def);
    FillHist("DEF_dE",DEFEnergy);
    FillHist("DEF_dEvsTOF",tof_def,DEFEnergy);

    FillHist("E0DEF_TOF",DEFTiming-E0Timing);

    // ==== SDC analysis ==== // 
    idc=2;
    const char* dcname = dlist->GetName(dccid[idc]).data();
    FillHist(Form("%s_nTrack",dcname),bltrackMan->ntrackBLDC(dccid[idc]));
    for(int i=0; i<bltrackMan->ntrackBLDC(dccid[idc]); i++){
      if(i>=1) break;
      LocalTrack* track = bltrackMan->trackBLDC(dccid[idc],i);
      double tracktime = track -> GetTrackTime();
      double chi2 = track -> chi2all();
      FillHist(Form("%s_TrackTime",dcname),tracktime);
      if(ll[idc]<tracktime&&tracktime<ul[idc]){
        FillHist(Form("%s_Chi2",dcname),chi2);
        if(chi2<chiul[idc]){
          nhitdc[idc]++;
          trackid[idc] = i;
          FillHist(Form("%s_TrackTimewCut",dcname),tracktime);
          FillHist(Form("%s_Chi2wCut",dcname),chi2);
        }
      }
    }
    FillHist(Form("%s_nGoodTrack",dcname),nhitdc[idc]);
    bool SDCGoodflag = false;
    // SDC 1 track
    if(nhitdc[idc]==1){
      LocalTrack* track = bltrackMan->trackBLDC(dccid[idc],trackid[idc]);
      double tracktime = track -> GetTrackTime();
      if(trigll[idc]<tracktime&&tracktime<trigul[idc])
        SDCGoodflag = true;
    }
    /* ## SDC 1 good track request ## */
    if(!SDCGoodflag) return false;
    FillHist("Beam_Count",icount); icount++;
    trackbpc  = bltrackMan->trackBPC(trackid[2]);
    TVector3 gpos;
    conf->GetGeomMapManager()->GetGPos(CID_T0,T0Segment,gpos);
    T0Position = trackbpc->GetPosatZ(gpos.Z()-0.5);

    BeamP = trackbpc->GetMomDir();
    BeamP.SetMag(BeamMom);
    BeamPID = Beam_Kaon;
    BeamMass = kpMass;
    BeamL.SetVectM( BeamP, BeamMass );
  }
  // ==== BLC2 and SDC maching ==== // 
  {
    double xll  = -0.75; double xul  = 0.75;
    double yll  = -0.75; double yul  = 0.75;
    double dxll = -0.02; double dxul = 0.02;
    double dyll = -0.02; double dyul = 0.02;
    TVector3 posblc2 = trackblc2->GetPosatZ(-75.0);
    TVector3 posbpc  = trackbpc->GetPosatZ(-75.0);
    TVector3 dirblc2 = trackblc2->GetMomDir();
    TVector3 dirbpc  = trackbpc->GetMomDir();
    TVector3 dist = posblc2-posbpc;
    double dxdz = dirblc2.X()/dirblc2.Z() - dirbpc.X()/dirbpc.Z();
    double dydz = dirblc2.Y()/dirblc2.Z() - dirbpc.Y()/dirbpc.Z();
    FillHist("BLC2_75Position",posblc2.X(),posblc2.Y());
    FillHist("SDC_75Position",posbpc.X(),posbpc.Y());
    FillHist("BLC2_Direction",dirblc2.X()/dirblc2.Z(),dirblc2.Y()/dirblc2.Z());
    FillHist("SDC_Direction",dirbpc.X()/dirbpc.Z(),dirbpc.Y()/dirbpc.Z());
    FillHist("BLC2SDC_Distance",dist.X(),dist.Y());
    FillHist("BLC2SDC_Distance",dist.X(),dist.Y());
    FillHist("BLC2SDC_Angle",dxdz,dydz);
    TVector3 posffblc2 = trackblc2->GetPosatZ(0.0);
    TVector3 posffbpc  = trackbpc->GetPosatZ(0.0);
    FillHist("BLC2_FFPosition",posffblc2.X(),posffblc2.Y());
    FillHist("SDC_FFPosition",posffbpc.X(),posffbpc.Y());
    /* ## BLC2 and SDC track maching ## */
    {
      //if( !((xll<dist.X()&&dist.X()<xul) && (yll<dist.Y()&&dist.Y()<yul)
      //      && (dxll<dxdz&&dxdz<dxul) && (dyll<dydz&&dydz<dyul)) )
      //  return false;
      //FillHist("Beam_Count",icount); icount++;
      FillHist("BLC2SDC_DistancewCut",dist.X(),dist.Y());
      FillHist("BLC2SDC_AnglewCut",dxdz,dydz);
    }
  }
  // ==== Beam position at FF ==== // 
  {
    TVector3 posff  = trackbpc->GetPosatZ(0.0);
    TVector3 posffblc2 = trackblc2->GetPosatZ(0.0);
    double dxdz = BeamP.X()/BeamP.Z();
    double dydz = BeamP.Y()/BeamP.Z();
    FillHist("Beam_FFPosition",posff.X(),posff.Y());
    FillHist("Beam_Angle",dxdz,dydz);
    /* ## On target selection ## */
    if(GeomTools::GetID(posff)==CID_Fiducial){
      FIDUCIAL = true;
      //FillHist("Beam_Count",icount); icount++;
      FillHist("Beam_FFPositionwCut",posff.X(),posff.Y());
      FillHist("Beam_AnglewCut",dxdz,dydz);
    }
    double par_cellring1[10]; GeomTools::GetParam(CID_CellRing,1,conf,par_cellring1);
    double par_cellring2[10]; GeomTools::GetParam(CID_CellRing,2,conf,par_cellring2);
    TVector3 poscr1 = trackbpc->GetPosatZ(par_cellring1[2]);
    FillHist("Beam_CellRing1Position",poscr1.X(),poscr1.Y());
    TVector3 poscr2 = trackbpc->GetPosatZ(par_cellring2[2]);
    FillHist("Beam_CellRing2Position",poscr2.X(),poscr2.Y());
  }

  pBeam* beam = new pBeam();
  beam->SetPID(Beam_Kaon);
  beam->SetT0seg(T0Segment);
  beam->SetT0Time(T0Timing);
  beam->SetT0Pos(T0Position);
  beam->SetBHDseg(BHDSegment);
  beam->SetBHDT0TOF(BeamTOF);
  beam->SetBHDX(BHDX);
  beam->SetBPCPos(trackbpc->GetPosatZ(0.0));
  beam->SetBPCDir(trackbpc->GetMomDir());
  TVector3 dir=trackblc2->GetMomDir();
  beam->SetBLCPos(trackblc2->GetPosatZ(0.0));
  beam->SetBLCDir(dir);
  beam->SetMomentum(BeamMom);

  particle->AddBeam(*beam);

  // ==== Start analysis ==== // 
  double StartTiming = -9999.0;
  double StartEnergy = -9999.0;
  double StartSegment = -1;
  int MulStart = 0;
  for(int i=0; i<blMan->nPA(); i++){
    if(blMan->PA(i)->CheckRange()){
      MulStart++;
      StartTiming = blMan->PA(i)->ctmean();
      StartEnergy = blMan->PA(i)->emean();
      StartSegment = blMan->PA(i)->seg();
    }
  }
  FillHist("Start_Mult",MulStart);
  /* ## Start 1hit request ## */
  if(MulStart!=1) return false;
  FillHist("Beam_Count",icount); icount++;
  double tof_start = StartTiming - T0Timing;
  FillHist("Start_TOF",tof_start);
  FillHist("Start_dE",StartEnergy);
  FillHist("Start_dEvsTOF",tof_start,StartEnergy);
  double tof_defstart = DEFTiming - StartTiming;
  FillHist("DEFStart_TOF",tof_defstart);

  // ==== Stop analysis ==== // 
  double StopTiming = -9999.0;
  double StopEnergy = -9999.0;
  double StopSegment = -1;
  int MulStop = 0;
  for(int i=0; i<blMan->nPC(); i++){
    if(blMan->PC(i)->CheckRange()){
      MulStop++;
      StopTiming = blMan->PC(i)->ctmean();
      StopEnergy = blMan->PC(i)->emean();
      StopSegment = blMan->PC(i)->seg();
    }
  }
  FillHist("Stop_Mult",MulStop);
  /* ## Stop 1hit request ## */
  if(MulStop!=1) return false;
  FillHist("Beam_Count",icount); icount++;
  double tof_stop = StopTiming - T0Timing;
  FillHist("Stop_TOF",tof_stop);
  FillHist("Stop_dE",StopEnergy);
  FillHist("Stop_dEvsTOF",tof_stop,StopEnergy);
  double tof_defstop = DEFTiming - StopTiming;
  FillHist("DEFStop_TOF",tof_defstop);
  double tof_startstop = StopTiming - StartTiming;
  FillHist("StartStop_TOF",tof_startstop);

  // ==== FDC1 analysis ==== // 
  // == RETIMING == //
  {
    const int cid= CID_FDC1;
    const char* name= dlist->GetName(cid).data();
    const int nlays= dlist->GetNlayers(cid);
    for(int layer=1;layer<=nlays;layer++){
      for(int i=0;i<blMan->nBLDC(cid,layer);i++){
        ChamberLikeHit *hit = blMan->BLDC(cid,layer,i);
        hit->SetTimeOffset(DEFTiming);
      }
    }
  }
  bltrackMan->LocalTracking(blMan,conf,CID_FDC1,"");

  const int ndc = 1;
  int dccid[ndc] = {CID_FDC1};
  int nhitdc[ndc] = {0};
  int trackid[ndc] = {-1};
  double ll[ndc] = {-30.0};
  double ul[ndc] = {100.0};
  double trigll[ndc] = {-20.0};
  double trigul[ndc] = { 20.0};
  double chiul[ndc] = {20.0};
  int idc=0;
  const char* dcname = dlist->GetName(dccid[idc]).data();
  FillHist(Form("%s_nTrack",dcname),bltrackMan->ntrackBLDC(dccid[idc]));
  for(int i=0; i<bltrackMan->ntrackBLDC(dccid[idc]); i++){
    if(i>=1) break;
    LocalTrack* track = bltrackMan->trackBLDC(dccid[idc],i);
    for(int ihit=0; ihit<track->nhit(); ihit++){
      ChamberLikeHit* hit = track->hit(ihit);
    }

    double tracktime = track -> GetTrackTime();
    double chi2 = track -> chi2all();
    FillHist(Form("%s_TrackTime",dcname),tracktime);
    if(ll[idc]<tracktime&&tracktime<ul[idc]){
      FillHist(Form("%s_Chi2",dcname),chi2);
      if(chi2<chiul[idc]){
        nhitdc[idc]++;
        trackid[idc] = i;
        FillHist(Form("%s_TrackTimewCut",dcname),tracktime);
        FillHist(Form("%s_Chi2wCut",dcname),chi2);
      }
    }
  }
  FillHist(Form("%s_nGoodTrack",dcname),nhitdc[idc]);
  bool FDC1Goodflag = false;
  // FDC1 1 track
  if(nhitdc[idc]==1){
    LocalTrack* track = bltrackMan->trackBLDC(dccid[idc],trackid[idc]);
    double tracktime = track -> GetTrackTime();
    if(trigll[idc]<tracktime&&tracktime<trigul[idc])
      FDC1Goodflag = true;
  }
  /* ## FDC1 1 good track request ## */
  if(!FDC1Goodflag) return false;
  FillHist("Beam_Count",icount); icount++;
  LocalTrack* trackfdc1  = bltrackMan->trackBPC(trackid[idc]);

  // ==== Vertex analysis ==== // 
  TVector3 sdcpos=trackbpc->GetPosatZ(0.0);
  TVector3 sdcdir=trackbpc->GetMomDir();
  TVector3 fdcpos=trackfdc1->GetPosatZ(0.0);
  TVector3 fdcdir=trackfdc1->GetMomDir();
  TVector3 xest=DEFVECT, nest=DEFVECT;
  double dl=0,dist=0; 
  MathTools::LineToLine(sdcpos,sdcdir,fdcpos,fdcdir,dl,dist,xest,nest);
  TVector3 vertex = xest;
  FillHist("Vertex_XY", vertex.X(), vertex.Y());
  FillHist("Vertex_ZX", vertex.Z(), vertex.X());
  FillHist("Vertex_ZY", vertex.Z(), vertex.Y());

  std::cout << "dX/dY = " << sdcdir.X()/sdcdir.Y() << ", dZ/dY = " << sdcdir.Z()/sdcdir.Y() << std::endl;


  return true;
}

bool MyAnalysisBL::FillHist(TString name, double val1)
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

bool MyAnalysisBL::FillHist(TString name, double val1, double val2)
{
  TH2F* h2 = (TH2F*)gFile -> Get(name);
  if(h2){
    h2 -> Fill(val1,val2);
    return true;
  }
  else {
    return false;
  }
}

bool MyAnalysisBL::Initialize(TFile* rtFile, ConfMan* confMan)
{
  rtFile -> cd();

  // Scaler
  new TH1F( "Scaler", "Scaler", 50, 0, 50 );

  // Hodoscope
  const int nhodo = 7;
  int NumOfSegments[nhodo] = {20,5,3,5,1,8,4};
  TString CounterName[nhodo] = {"BHD","T0","E0","DEF","BVC","Start","Stop"};
  std::cout << "Define Histograms for Hodoscope" << std::endl;
  for(int ihodo=0; ihodo<nhodo; ihodo++){
    new TH1F( Form("%s_HitPat",CounterName[ihodo].Data()), Form("Hit Pattern %s;%s segment",CounterName[ihodo].Data(),CounterName[ihodo].Data()), NumOfSegments[ihodo]+1, 0, NumOfSegments[ihodo]+1 );
    new TH1F( Form("%s_Mult",CounterName[ihodo].Data()), Form("Multiplicity %s;%s Multiplicity;Counts",CounterName[ihodo].Data(),CounterName[ihodo].Data()), NumOfSegments[ihodo]+1, 0, NumOfSegments[ihodo]+1 );
    new TH1F( Form("%s_MultwCut",CounterName[ihodo].Data()), Form("Multiplicity %s;%s Multiplicity;Counts",CounterName[ihodo].Data(),CounterName[ihodo].Data()), NumOfSegments[ihodo]+1, 0, NumOfSegments[ihodo]+1 );
  }
  // Chamber
  const int ndc=4;
  int dccid[ndc]={CID_BLC1,CID_BLC2,CID_BPC,CID_FDC1};
  TString dcname[ndc]={"BLC1","BLC2","SmallDC","FDC1"};
  std::cout << "Define Histgram for DC" << std::endl;
  for(int idc=0;idc<ndc;idc++){
    new TH1F( Form("%s_nTrack",dcname[idc].Data()), Form("Num. of tracks %s;Number of Tracks;Counts",dcname[idc].Data()), 10, 0.0, 10.0 );
    new TH1F( Form("%s_nGoodTrack",dcname[idc].Data()), Form("Num. of tracks %s;Number of Good Tracks;Counts",dcname[idc].Data()), 10, 0.0, 10.0 );
    new TH1F( Form("%s_nBeamTrack",dcname[idc].Data()), Form("Num. of tracks %s;Number of Beam Tracks;Counts",dcname[idc].Data()), 10, 0.0, 10.0 );
    new TH1F( Form("%s_TrackTime",dcname[idc].Data()), Form("Track time %s;Time (ns);Counts",dcname[idc].Data()), 1200, -200.0, 400.0 );
    new TH1F( Form("%s_TrackTimewCut",dcname[idc].Data()), Form("Track time %s;Time (ns);Counts",dcname[idc].Data()), 1200, -200.0, 400.0 );
    new TH1F( Form("%s_Chi2",dcname[idc].Data()), Form("#chi^{2} %s;#chi^{2}/ndf;Counts",dcname[idc].Data()), 1000, 0.0, 100.0 );
    new TH1F( Form("%s_Chi2wCut",dcname[idc].Data()), Form("#chi^{2} %s;#chi^{2}/ndf;Counts",dcname[idc].Data()), 1000, 0.0, 100.0 );
  }
  // TOF
  std::cout << "Define Histograms for TOF" << std::endl;
  new TH2F( Form("T0%s_TOF",CounterName[0].Data()), Form("TOF T0-%s;%s segment;TOF (ns)",CounterName[0].Data(),CounterName[0].Data()), NumOfSegments[0], 1, NumOfSegments[0]+1, 4000, -50, 50 );
  new TH1F( "E0_TOF", "E0 TOF;TOF between E0 and T0 (ns);Counts",1000, -5.0, 20.0 );
  new TH1F( "E0_dE", "E0 Energy Deposit;Energy Deposit (MeVee);Counts",1000, 0.0, 30.0 );
  new TH2F( "E0_dEvsTOF", "E0 Energy Deposit vs. TOF;TOF between E0 and T0 (ns);Energy Deposit (MeVee)",1000, -5, 20, 1000, 0.0, 30.0 );
  new TH1F( "DEF_TOF", "DEF TOF;TOF between DEF and T0 (ns);Counts",1000, -5.0, 20.0 );
  new TH1F( "DEF_dE", "DEF Energy Deposit;Energy Deposit (MeVee);Counts",1000, 0.0, 5.0 );
  new TH2F( "DEF_dEvsTOF", "DEF Energy Deposit vs. TOF;TOF between DEF and T0 (ns);Energy Deposit (MeVee)",1000, -5, 20, 1000, 0.0, 5.0 );
  new TH1F( "E0DEF_TOF", "E0 to DEF TOF;TOF between DEF and E0 (ns);Counts",1000, -5.0, 20.0 );
  new TH1F( "BVC_TOF", "BVC TOF;TOF between BVC and T0 (ns);Counts",1000, -5.0, 20.0 );
  new TH1F( "BVC_dE", "BVC Energy Deposit;Energy Deposit (MeVee);Counts",1000, 0.0, 20.0 );
  new TH2F( "BVC_dEvsTOF", "BVC Energy Deposit vs. TOF;TOF between BVC and T0 (ns);Energy Deposit (MeVee)",1000, -5, 20, 1000, 0.0, 20.0 );
  new TH1F( "Start_TOF", "Start TOF;TOF between Start and T0 (ns);Counts",1000, -5.0, 20.0 );
  new TH1F( "Start_dE", "Start Energy Deposit;Energy Deposit (MeVee);Counts",1000, 0.0, 20.0 );
  new TH2F( "Start_dEvsTOF", "Start Energy Deposit vs. TOF;TOF between Start and T0 (ns);Energy Deposit (MeVee)",1000, -5, 20, 1000, 0.0, 20.0 );
  new TH1F( "DEFStart_TOF", "DEF Start TOF;TOF between Start and DEF (ns);Counts",1000, -5.0, 20.0 );
  new TH1F( "Stop_TOF", "Stop TOF;TOF between Stop and T0 (ns);Counts",1000, -5.0, 20.0 );
  new TH1F( "Stop_dE", "Stop Energy Deposit;Energy Deposit (MeVee);Counts",1000, 0.0, 25.0 );
  new TH2F( "Stop_dEvsTOF", "Stop Energy Deposit vs. TOF;TOF between Stop and T0 (ns);Energy Deposit (MeVee)",1000, -5, 20, 1000, 0.0, 25.0 );
  new TH1F( "DEFStop_TOF", "DEF Stop TOF;TOF between Stop and DEF (ns);Counts",1000, -5.0, 20.0 );
  new TH1F( "StartStop_TOF", "Start to Stop TOF;TOF between Stop and Start (ns);Counts",1000, -5.0, 20.0 );
  // Beam analysis
  std::cout << "Define Histograms for Beam analysis" << std::endl;
  new TH1F( "Beam_Count", "Beam counts;Selection;Counts", 15, 0, 15 );
  new TH1F( "Beam_Momentum", "Beam momentum;Beam momentum (GeV/c);Counts",1000, 0.5, 1.5 );
  new TH1F( "Beam_MomentumwCut", "Beam momentum;Beam momentum (GeV/c);Counts",1000, 0.5, 1.5 );
  new TH1F( "Beam_TOF", "Beam TOF;TOF between BHD and T0 (ns);Counts",1000, -10.0, 10.0 );
  new TH1F( "Beam_TOFwCut", "Beam TOF;TOF between BHD and T0 (ns);Counts",1000, -10.0, 10.0 );
  new TH2F( "Beam_MomentumvsTOF", "Beam Momentum vs. TOF;Beam Momentum (GeV/c);TOF between BHD and T0 (ns)",1000,0.5,1.5,1000, 20.0, 40.0 );
  new TH2F( "Beam_MomentumvsTOFwCut", "Beam Momentum vs. TOF;Beam Momentum (GeV/c);TOF between BHD and T0 (ns)",1000,0.5,1.5,1000, 20.0, 40.0 );
  new TH1F( "Beam_Chi2", "Beam mom #chi^{2};#chi^{2}/ndf;Counts",1000, 0.0, 100.0 );
  new TH1F( "Beam_Chi2wCut", "Beam mom #chi^{2};#chi^{2}/ndf;Counts",1000, 0.0, 100.0 );
  new TH2F( "Beam_FFPosition", "Position at FF;X position (cm);Y position (cm)",600,-15,15,600,-15,15);
  new TH2F( "Beam_FFPositionwCut", "Position at FF;X position (cm);Y position (cm)",600,-15,15,600,-15,15);
  new TH2F( "Beam_CellRing1Position", "Position at CellRing1;X position (cm);Y position (cm)",600,-15,15,600,-15,15);
  new TH2F( "Beam_CellRing2Position", "Position at CellRing2;X position (cm);Y position (cm)",600,-15,15,600,-15,15);
  new TH2F( "Beam_Angle", "Beam Angle;dx/dz;dy/dz",1000,-0.05,0.05,1000,-0.05,0.05);
  new TH2F( "Beam_AnglewCut", "Beam Angle;dx/dz;dy/dz",1000,-0.05,0.05,1000,-0.05,0.05);
  new TH2F( "BHD_SegmentvsBeamMomentum", "BHD sement vs. Beam momentum;BHD segment;Beam momentum (GeV/c)", 20, 1, 21, 1000, 0.5, 1.5 );
  new TH2F( "BHD_SegmentvsBeamMomentumwCut", "BHD sement vs. Beam momentum;BHD segment;Beam momentum (GeV/c)", 20, 1, 21, 1000, 0.5, 1.5 );

  new TH2F( "BLC2_FFPosition", "BLC2 Track Position at FF;X position (cm);Y position (cm)",600,-15,15,600,-15,15);
  new TH2F( "BLC2_75Position", "BLC2 Track Position at z=-75.0 cm;X position (cm);Y position (cm)",600,-15,15,600,-15,15);
  new TH2F( "SDC_FFPosition", "SDC Track Position at FF;X position (cm);Y position (cm)",600,-15,15,600,-15,15);
  new TH2F( "SDC_75Position", "SDC Track Position at z=-75.0 cm;X position (cm);Y position (cm)",600,-15,15,600,-15,15);
  new TH2F( "BLC2_Direction", "BLC2 Track Direction;dX/dZ;dY/dZ",1000,-0.05,0.05,1000,-0.05,0.05);
  new TH2F( "SDC_Direction", "SDC Track Direction;dX/dZ;dY/dZ",1000,-0.05,0.05,1000,-0.05,0.05);
  new TH2F( "BLC2SDC_Distance", "Distance between BLC2 and SDC;X distance (cm);Y distance (cm)",1000,-10.0,10.0,1000,-10.0,10.0 );
  new TH2F( "BLC2SDC_DistancewCut", "Distance between BLC2 and SDC;X distance (cm);Y distance (cm)",1000,-10.0,10.0,1000,-10.0,10.0 );
  new TH2F( "BLC2SDC_Angle", "Angle between BLC2 and SDC;dx/dz;dy/dz",1000,-0.20,0.20,1000,-0.20,0.20);
  new TH2F( "BLC2SDC_AnglewCut", "Angle between BLC2 and SDC;dx/dz;dy/dz",1000,-0.20,0.20,1000,-0.20,0.20);
  // Vertex
  new TH2F( "Vertex_XY", "Vertex XY plane;x (cm);y (cm)", 600, -15, 15, 600, -15, 15 );
  new TH2F( "Vertex_ZX", "Vertex ZX plane;z (cm);x (cm)", 1200, -30, 30, 600, -15, 15 );
  new TH2F( "Vertex_ZY", "Vertex ZY plane;z (cm);y (cm)", 1200, -30, 30, 600, -15, 15 );

  std::cout << "=== End of [MyAnalysisBL::Initialize] === " << std::endl;

  return true;
}

