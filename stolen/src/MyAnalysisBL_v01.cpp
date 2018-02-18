// MyAnalysisBL.cpp

#include "MyAnalysisBL.h"

MyAnalysisBL::MyAnalysisBL(TFile* rtFile, ConfMan* conf)
{
  D5 = new BeamSpectrometer(conf);
  Initialize(rtFile, conf);
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

  return true;
}

bool MyAnalysisBL::DoAnalysis(ConfMan* conf, EventHeader* header, BeamLineHitMan* blMan, BeamLineTrackMan* bltrackMan)
{
  DetectorList *dlist=DetectorList::GetInstance();

  // ############# //
  // Trigger check //
  // ############# //
  /* ## Select K/f trigger ## */
  if(!header->IsTrig(Trig_Kaon)) return false;
  //if(!header->IsTrig(Trig_Kf)) return false;

  /* ## All events ## */
  FillHist("Beam_Count",0);

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
    FillHist("Beam_Count",1);
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
        MulBHD++;
        seg[MulBHD-1] = blMan->BHD(i)->seg();
        tof[MulBHD-1] = T0Timing - blMan->BHD(i)->ctmean();;
        FillHist("T0BHD_TOF",seg[MulBHD-1],tof[MulBHD-1]);
      }
    }
    FillHist("BHD_Mult",MulBHD);
    /* ## TOF kaon request ## */
    double tofk[2] = {27.0, 31.0};
    bool TOFflag=false;
    for(int i=0; i<MulBHD; i++){
      if(tofk[0]<tof[i]&&tof[i]<tofk[1]) {
        TOFflag=true;
        BeamTOF = tof[i];
        BHDSegment = seg[i];
        break;
      }
    }
    FillHist(Form("Beam_TOF"), BeamTOF);
    if(!TOFflag) return false;
    FillHist("Beam_Count",2);
    FillHist("BHD_MultwCut",MulBHD);
    FillHist(Form("Beam_TOFwCut"), BeamTOF);
  }
  // #################### //
  // DC tracking          //
  // #################### //
  {
    bltrackMan->DoTracking(blMan,conf,true,true);
    const int ndc = 3;
    int dccid[ndc] = {CID_BLC1,CID_BLC2,CID_BPC};
    int nhitdc[ndc] = {0,0,0};
    int trackid[ndc] = {-1,-1,-1};
    double ll[ndc] = {-30.0,-30.0,-30.0};
    double ul[ndc] = {100.0,100.0,100.0};
    double trigll[ndc] = {-20.0,-20.0,-20.0};
    double trigul[ndc] = { 20.0, 20.0, 20.0};
    double chiul[ndc] = {30.0, 30.0, 30.0};
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
    FillHist("Beam_Count",3);
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
    /* ## Beam chisquare < 30 request ## */
    if(!(BeamChi2<30)) return false;
    FillHist("Beam_Count",4);
    FillHist(Form("Beam_MomentumwCut"), BeamMom);
    FillHist(Form("Beam_MomentumvsTOFwCut"), BeamMom,BeamTOF);
    FillHist(Form("Beam_Chi2wCut"), BeamChi2);
    FillHist(Form("BHD_SegmentvsBeamMomentumwCut"),BHDSegment,BeamMom); 

    // ==== BPC analysis ==== // 
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
    bool BPCGoodflag = false;
    // BPC 1 track
    if(nhitdc[idc]==1){
      LocalTrack* track = bltrackMan->trackBLDC(dccid[idc],trackid[idc]);
      double tracktime = track -> GetTrackTime();
      if(trigll[idc]<tracktime&&tracktime<trigul[idc])
        BPCGoodflag = true;
    }
    /* ## BPC 1 good track request ## */
    if(!BPCGoodflag) return false;
    FillHist("Beam_Count",5);
    trackbpc  = bltrackMan->trackBPC(trackid[2]);
  }
  BeamP = trackbpc->GetMomDir();
  BeamP.SetMag(BeamMom);
  BeamPID = Beam_Kaon;
  BeamMass = kpMass;
  BeamL.SetVectM( BeamP, BeamMass );
  // ==== BLC2 and BPC maching ==== // 
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
    FillHist("BLC2BPC_Distance",dist.X(),dist.Y());
    FillHist("BLC2BPC_Angle",dxdz,dydz);
    /* ## BLC2 and BPC track maching ## */
    if( !((xll<dist.X()&&dist.X()<xul) && (yll<dist.Y()&&dist.Y()<yul)
          && (dxll<dxdz&&dxdz<dxul) && (dyll<dydz&&dydz<dyul)) )
      return false;
    FillHist("Beam_Count",6);
    FillHist("BLC2BPC_DistancewCut",dist.X(),dist.Y());
    FillHist("BLC2BPC_AnglewCut",dxdz,dydz);
  }
  // ==== Beam position at FF ==== // 
  {
    TVector3 posff  = trackbpc->GetPosatZ(0.0);
    double dxdz = BeamP.X()/BeamP.Z();
    double dydz = BeamP.Y()/BeamP.Z();
    FillHist("Beam_FFPosition",posff.X(),posff.Y());
    FillHist("Beam_Angle",dxdz,dydz);
    /* ## BLC2 and BPC track maching ## */
    //if(GeomTools::GetID(posff)!=CID_Fiducial) return false;
    FillHist("Beam_Count",7);
    FillHist("Beam_FFPositionwCut",posff.X(),posff.Y());
    FillHist("Beam_AnglewCut",dxdz,dydz);
  }

  //  // ################ //
  //  // Hodoscope        //
  //  // ################ //
  //  const int nhodo=4;
  //  int hodoid[nhodo]={CID_BHD,CID_T0,CID_BPD,CID_DEF};
  //  int hodonHit[nhodo]; for(int i=0; i<nhodo; i++) hodonHit[i]=0;
  //  int hodoseg[nhodo][100];
  //  double hodoemean[nhodo][100];
  //  double hodoctmean[nhodo][100];
  //  double hodohitpos[nhodo][100];
  //  double hodotof[nhodo][100];
  //  TVector3 hodopos[nhodo][100];
  //  for(int ihodo=0; ihodo<nhodo; ihodo++){
  //    const int cid = hodoid[ihodo];
  //    const char* name = dlist->GetName(cid).data();
  //    const int nsegs = dlist->GetNsegs(cid);
  //    int ihit = 0;
  //    for(int i=0; i<blMan->nHodo(cid); i++){
  //      HodoscopeLikeHit* hit = blMan->Hodoi(cid,i);
  //      if(hit->CheckRange()){
  //        hodoseg[ihodo][ihit] = hit->seg();
  //        hodoemean[ihodo][ihit] = hit->emean();
  //        hodoctmean[ihodo][ihit] = hit->ctmean();
  //        hodohitpos[ihodo][ihit] = hit->hitpos();
  //        FillHist(Form("%s_HitPat",name),hodoseg[ihodo][ihit]);
  //        FillHist(Form("%s_SegmentvsPosition",name),hodoseg[ihodo][ihit],hodohitpos[ihodo][ihit]);
  //        // === TOF calculation === //
  //        if(T0Timing>-900){
  //          hodotof[ihodo][i] = hodoctmean[ihodo][ihit]-T0Timing;
  //          if(ihodo==0){ hodotof[ihodo][i]=hodotof[ihodo][i]*-1.; }
  //          FillHist(Form("T0%s_TOF",name),hodoseg[ihodo][ihit],hodotof[ihodo][i]);
  //        }
  //        // === TOF calculation === //
  //        hodonHit[ihodo]++;
  //        ihit++;
  //      }
  //    }
  //    FillHist(Form("%s_Mult",name),hodonHit[ihodo]);
  //  }
  //
  //  // ################ //
  //  // Chamber          //
  //  // ################ //
  //  bltrackMan->DoTracking(blMan,conf,true,false);
  //  const int ndc = 3;
  //  int dccid[ndc] = {CID_BLC1,CID_BLC2,CID_BPC};
  //  double ll[ndc] = {-20.0,-20.0,-20.0};
  //  double ul[ndc] = { 20.0, 20.0, 20.0};
  //  double chiul[ndc] = {30.0, 30.0, 30.0};
  //  int nhitdc[ndc] = { 0, 0, 0};
  //  int trackid[ndc] = { -1, -1, -1};
  //  for(int idc=0;idc<ndc;idc++){
  //    const int cid = dccid[idc];
  //    const char* name = dlist->GetName(cid).data();
  //    for(int i=0; i<bltrackMan->ntrackBLDC(cid); i++){
  //      LocalTrack* track = bltrackMan->trackBLDC(cid,i);
  //      double tracktime = track -> GetTrackTime();
  //      double chi2 = track -> chi2all();
  //      FillHist(Form("%s_TrackTime",name),tracktime);
  //      FillHist(Form("%s_Chi2",name),chi2);
  //      if(ll[idc]<=tracktime&&tracktime<=ul[idc] && chi2<=chiul[idc]){
  //        nhitdc[idc]++;        /* Num. of hit in cut region */
  //        trackid[idc] = i;     /* Selected track ID */
  //        FillHist(Form("%s_TrackTimewCut",name),tracktime);
  //        FillHist(Form("%s_Chi2wCut",name),chi2);
  //      }
  //    }
  //    FillHist(Form("%s_nTrack",name),nhitdc[idc]);
  //  }
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

  // Hodoscope
  const int nhodo = 4;
  int NumOfSegments[nhodo] = {20,5,70,8};
  TString CounterName[nhodo] = {"BHD","T0","BPD","DEF"};
  std::cout << "Define Histograms for Hodoscope" << std::endl;
  for(int ihodo=0; ihodo<2; ihodo++){
    new TH1F( Form("%s_HitPat",CounterName[ihodo].Data()), Form("Hit Pattern %s;%s segment",CounterName[ihodo].Data(),CounterName[ihodo].Data()), NumOfSegments[ihodo]+1, 0, NumOfSegments[ihodo]+1 );
    new TH1F( Form("%s_Mult",CounterName[ihodo].Data()), Form("Multiplicity %s;Multiplicity;Counts",CounterName[ihodo].Data()), NumOfSegments[ihodo]+1, 0, NumOfSegments[ihodo]+1 );
    new TH1F( Form("%s_MultwCut",CounterName[ihodo].Data()), Form("Multiplicity %s;Multiplicity;Counts",CounterName[ihodo].Data()), NumOfSegments[ihodo]+1, 0, NumOfSegments[ihodo]+1 );
  }
  // Chamber
  const int ndc=3;
  int dccid[ndc]={CID_BLC1,CID_BLC2,CID_BPC};
  TString dcname[ndc]={"BLC1","BLC2","BPC"};
  std::cout << "Define Histgram for DC" << std::endl;
  for(int idc=0;idc<ndc;idc++){
    new TH1F( Form("%s_nTrack",dcname[idc].Data()), Form("Num. of tracks %s;Time (ns);Counts",dcname[idc].Data()), 10, 0.0, 10.0 );
    new TH1F( Form("%s_nGoodTrack",dcname[idc].Data()), Form("Num. of tracks %s;Time (ns);Counts",dcname[idc].Data()), 10, 0.0, 10.0 );
    new TH1F( Form("%s_nBeamTrack",dcname[idc].Data()), Form("Num. of tracks %s;Time (ns);Counts",dcname[idc].Data()), 10, 0.0, 10.0 );
    new TH1F( Form("%s_TrackTime",dcname[idc].Data()), Form("Track time %s;Time (ns);Counts",dcname[idc].Data()), 1200, -200.0, 400.0 );
    new TH1F( Form("%s_TrackTimewCut",dcname[idc].Data()), Form("Track time %s;Time (ns);Counts",dcname[idc].Data()), 1200, -200.0, 400.0 );
    new TH1F( Form("%s_Chi2",dcname[idc].Data()), Form("#chi^{2} %s;#chi^{2}/ndf;Counts",dcname[idc].Data()), 1000, 0.0, 100.0 );
    new TH1F( Form("%s_Chi2wCut",dcname[idc].Data()), Form("#chi^{2} %s;#chi^{2}/ndf;Counts",dcname[idc].Data()), 1000, 0.0, 100.0 );
  }
  // TOF
  std::cout << "Define Histograms for TOF" << std::endl;
  new TH2F( Form("T0%s_TOF",CounterName[0].Data()), Form("TOF T0-%s;%s segment;TOF (ns)",CounterName[0].Data(),CounterName[0].Data()), NumOfSegments[0], 1, NumOfSegments[0]+1, 4000, -50, 50 );
  // Beam analysis
  std::cout << "Define Histograms for Beam analysis" << std::endl;
  new TH1F( "Beam_Count", "Beam counts;Selection;Counts", 10, 0, 10 );
  new TH1F( "Beam_Momentum", "Beam momentum;Beam momentum (GeV/c);Counts",1000, 0.5, 1.5 );
  new TH1F( "Beam_MomentumwCut", "Beam momentum;Beam momentum (GeV/c);Counts",1000, 0.5, 1.5 );
  new TH1F( "Beam_TOF", "Beam TOF;TOF between BHD and T0 (ns);Counts",1000, 20.0, 40.0 );
  new TH1F( "Beam_TOFwCut", "Beam TOF;TOF between BHD and T0 (ns);Counts",1000, 20.0, 40.0 );
  new TH2F( "Beam_MomentumvsTOF", "Beam Momentum vs. TOF;Beam Momentum (GeV/c);TOF between BHD and T0 (ns)",1000,0.5,1.5,1000, 20.0, 40.0 );
  new TH2F( "Beam_MomentumvsTOFwCut", "Beam Momentum vs. TOF;Beam Momentum (GeV/c);TOF between BHD and T0 (ns)",1000,0.5,1.5,1000, 20.0, 40.0 );
  new TH1F( "Beam_Chi2", "Beam mom #chi^{2};#chi^{2}/ndf;Counts",1000, 0.0, 100.0 );
  new TH1F( "Beam_Chi2wCut", "Beam mom #chi^{2};#chi^{2}/ndf;Counts",1000, 0.0, 100.0 );
  new TH2F( "Beam_FFPosition", "Position at FF;X position (cm);Y position (cm)",600,-15,15,600,-15,15);
  new TH2F( "Beam_FFPositionwCut", "Position at FF;X position (cm);Y position (cm)",600,-15,15,600,-15,15);
  new TH2F( "Beam_Angle", "Beam Angle;dx/dz;dy/dz",1000,-0.05,0.05,1000,-0.05,0.05);
  new TH2F( "Beam_AnglewCut", "Beam Angle;dx/dz;dy/dz",1000,-0.05,0.05,1000,-0.05,0.05);
  new TH2F( "BHD_SegmentvsBeamMomentum", "BHD sement vs. Beam momentum;BHD segment;Beam momentum (GeV/c)", 20, 1, 21, 1000, 0.5, 1.5 );
  new TH2F( "BHD_SegmentvsBeamMomentumwCut", "BHD sement vs. Beam momentum;BHD segment;Beam momentum (GeV/c)", 20, 1, 21, 1000, 0.5, 1.5 );

  new TH2F( "BLC2BPC_Distance", "Distance between BLC2 and BPC;X distance (cm);Y distance (cm)",1000,-5.0,5.0,1000,-5.0,5.0 );
  new TH2F( "BLC2BPC_DistancewCut", "Distance between BLC2 and BPC;X distance (cm);Y distance (cm)",1000,-5.0,5.0,1000,-5.0,5.0 );
  new TH2F( "BLC2BPC_Angle", "Angle between BLC2 and BPC;dx/dz;dy/dz",1000,-0.05,0.05,1000,-0.05,0.05);
  new TH2F( "BLC2BPC_AnglewCut", "Angle between BLC2 and BPC;dx/dz;dy/dz",1000,-0.05,0.05,1000,-0.05,0.05);

  return true;
}

