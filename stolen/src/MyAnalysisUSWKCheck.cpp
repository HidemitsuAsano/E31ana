// MyAnalysisUSWKCheck.cpp

#include "MyAnalysisUSWKCheck.h"

MyAnalysisUSWKCheck::MyAnalysisUSWKCheck(TFile* rt, ConfMan* conf)
{
  Initialize(conf);
  Clear();
}

MyAnalysisUSWKCheck::~MyAnalysisUSWKCheck()
{
  Clear();
  rtFile->cd();
  rtFile->Write();
  rtFile->Close();
}

void MyAnalysisUSWKCheck::Clear()
{
  return;
}

bool MyAnalysisUSWKCheck::DoAnalysis(ConfMan* conf, EventHeader* header, BeamLineHitMan* blMan, BeamLineTrackMan* bltrackMan)
{
  rtFile->cd();
  DetectorList *dlist=DetectorList::GetInstance();

  /* ======================================= */
  /* ========== Trigger selection ========== */
  if(!header->IsTrig(Trig_Pion)){
    return true;
  }
  //if(!header->IsTrig(Trig_Kaon)){
  //  return true;
  //}
  /* ========== Trigger selection ========== */
  /* ======================================= */

  /* ===================================== */
  /* ========== Event selection ========== */
  double T0Timing = -9999.; 
  int T0Segment = -1;
  int MulT0 = 0;
  for(int i=0; i<blMan->nT0(); i++){
    if(blMan->T0(i)->CheckRange()){
      MulT0++;
      T0Timing = blMan->T0(i)->ctmean();
      T0Segment = blMan->T0(i)->seg();
    }
  }
  if(MulT0!=1){
    return true;
  }
  //if(T0Segment!=3){
  //  return true;
  //}
  int MulBHD = 0;
  double BHDTiming = -9999.;
  for(int i=0; i<blMan->nBHD(); i++){
    if(blMan->BHD(i)->CheckRange()){
      if(6<=blMan->BHD(i)->seg()&&blMan->BHD(i)->seg()<=15){
        BHDTiming = blMan->BHD(i)->ctmean();
        MulBHD++;
      }
    }
  }
  if(MulBHD!=1){
    return true;
  }
  double BeamTOF = T0Timing-BHDTiming;
  FillHist(Form("Beam_TOF"), BeamTOF);
  if((T0Timing-BHDTiming)<25||27<(T0Timing-BHDTiming)){ /* Pion */
    return true;
  }
  //if((T0Timing-BHDTiming)<34||36<(T0Timing-BHDTiming)){ /* Proton */
  //  return true;
  //}
  FillHist(Form("Beam_TOFwithCut"), BeamTOF);
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
  int MulNC=0;
  for(int i=0; i<blMan->nNC(); i++){
    HodoscopeLikeHit* hit = blMan->NC(i);
    if(hit->CheckRange()) MulNC++;
  }
  if(MulCVC==0&&MulPC==0) return true;

  // #################### //
  // DC tracking          //
  // #################### //
  bltrackMan->DoTracking(blMan,conf,true,true);
  const int ndc = 4;
  int dccid[ndc] = {CID_BLC1,CID_BLC2,CID_BPC,CID_FDC1};
  int nhitdc[ndc] = {0,0,0,0};
  int trackid[ndc] = {-1,-1,-1,-1};
  double ll[ndc] = {-30.0,-30.0,-30.0,-30.0};
  double ul[ndc] = {100.0,100.0,100.0,100.0};
  double trigll[ndc] = {-10.0,-10.0,-10.0,-10.0};
  double trigul[ndc] = { 10.0, 10.0, 10.0, 10.0};
  double chiul[ndc] = {10.0, 10.0, 10.0, 10.0};
  int idc;
  const char* dcname;
  // ==== BLC1 and BLC2 analysis ==== // 
  for(idc=0;idc<2;idc++){
    dcname = dlist->GetName(dccid[idc]).data();
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
  LocalTrack* trackblc1 = bltrackMan->trackBLC1(trackid[0]);
  LocalTrack* trackblc2 = bltrackMan->trackBLC2(trackid[1]);
  // ==== Beam momentum analysis ==== // 
  BeamSpectrometer* D5 = new BeamSpectrometer(conf);
  D5->TMinuitFit(trackblc1,trackblc2,conf);
  double BeamChi2 = D5->chisquare();
  double BeamMom = D5->mom();
  FillHist(Form("Beam_Momentum"), BeamMom);
  FillHist(Form("Beam_MomentumvsTOF"), BeamMom,BeamTOF);
  FillHist(Form("Beam_Chi2"), BeamChi2);
  delete D5;
  /* ## Beam chisquare < 20 request ## */
  if(!(BeamChi2<20)) return false;
  FillHist(Form("Beam_MomentumwCut"), BeamMom);
  FillHist(Form("Beam_MomentumvsTOFwCut"), BeamMom,BeamTOF);
  FillHist(Form("Beam_Chi2wCut"), BeamChi2);

  // ==== BPC analysis ==== // 
  idc=2;
  dcname = dlist->GetName(dccid[idc]).data();
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
  LocalTrack* trackbpc  = bltrackMan->trackBPC(trackid[2]);
  TVector3 gpos;
  conf->GetGeomMapManager()->GetGPos(CID_T0,T0Segment,gpos);
  TVector3 T0Position = trackbpc->GetPosatZ(gpos.Z()-0.5);
  TVector3 VertexPosition = trackbpc->GetPosatZ(0.0);
  TVector3 BeamP = trackbpc->GetMomDir();
  BeamP.SetMag(BeamMom);
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
    FillHist("BLC2_75Position",posblc2.X(),posblc2.Y());
    FillHist("BPC_75Position",posbpc.X(),posbpc.Y());
    FillHist("BLC2_Direction",dirblc2.X()/dirblc2.Z(),dirblc2.Y()/dirblc2.Z());
    FillHist("BPC_Direction",dirbpc.X()/dirbpc.Z(),dirbpc.Y()/dirbpc.Z());
    FillHist("BLC2BPC_Distance",dist.X(),dist.Y());
    FillHist("BLC2BPC_Distance",dist.X(),dist.Y());
    FillHist("BLC2BPC_Angle",dxdz,dydz);
    TVector3 posffblc2 = trackblc2->GetPosatZ(0.0);
    TVector3 posffbpc  = trackbpc->GetPosatZ(0.0);
    FillHist("BLC2_FFPosition",posffblc2.X(),posffblc2.Y());
    FillHist("BPC_FFPosition",posffbpc.X(),posffbpc.Y());
    /* ## BLC2 and BPC track maching ## */
    if( !((xll<dist.X()&&dist.X()<xul) && (yll<dist.Y()&&dist.Y()<yul)
          && (dxll<dxdz&&dxdz<dxul) && (dyll<dydz&&dydz<dyul)) )
      return false;
    FillHist("Beam_Count",6);
    FillHist("BLC2BPC_DistancewCut",dist.X(),dist.Y());
    FillHist("BLC2BPC_AnglewCut",dxdz,dydz);
  }
  // ==== Beam position at FF ==== // 
  bool FIDUCIAL = false;
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
      FillHist("Beam_Count",7);
      FillHist("Beam_FFPositionwCut",posff.X(),posff.Y());
      FillHist("Beam_AnglewCut",dxdz,dydz);
    }
  }
  if(!FIDUCIAL){ return false; }
  TVector3 vertex = trackbpc->GetPosatZ(0.0);

  // ==== FDC analysis ==== // 
  idc=3;
  dcname = dlist->GetName(dccid[idc]).data();
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
  LocalTrack* trackfdc1  = bltrackMan->trackFDC1(trackid[3]);
  TVector3 USWKPosition = trackfdc1->GetPosatZ(169.6);
  // ==== BPC and FDC1 maching ==== // 
  {
    double xll  = -2.00; double xul  = 2.00;
    double yll  = -2.00; double yul  = 2.00;
    double dxll = -9999.9; double dxul = 9999.9;
    double dyll = -9999.9; double dyul = 9999.9;
    TVector3 posfdc1 = trackfdc1->GetPosatZ(169.6);
    TVector3 posbpc  = trackbpc->GetPosatZ(169.6);
    TVector3 dirfdc1 = trackfdc1->GetMomDir();
    TVector3 dirbpc  = trackbpc->GetMomDir();
    TVector3 dist = posfdc1-posbpc;
    double dxdz = dirfdc1.X()/dirfdc1.Z() - dirbpc.X()/dirbpc.Z();
    double dydz = dirfdc1.Y()/dirfdc1.Z() - dirbpc.Y()/dirbpc.Z();
    FillHist("FDC1_USWKPosition",posfdc1.X(),posfdc1.Y());
    FillHist("BPC_USWKPosition",posbpc.X(),posbpc.Y());
    FillHist("FDC1_Direction",dirfdc1.X()/dirfdc1.Z(),dirfdc1.Y()/dirfdc1.Z());
    FillHist("FDC1BPC_Distance",dist.X(),dist.Y());
    FillHist("FDC1BPC_Distance",dist.X(),dist.Y());
    FillHist("FDC1BPC_Angle",dxdz,dydz);
    /* ## BPC and FDC1 track maching ## */
    if( !((xll<dist.X()&&dist.X()<xul) && (yll<dist.Y()&&dist.Y()<yul)
          && (dxll<dxdz&&dxdz<dxul) && (dyll<dydz&&dydz<dyul)) )
      return false;
    FillHist("FDC1BPC_DistancewCut",dist.X(),dist.Y());
    FillHist("FDC1BPC_AnglewCut",dxdz,dydz);
  }
  // ==== Forward particle position and direction ==== // 
  {
    TVector3 fdcpos = trackfdc1->GetPosatZ(169.6);
    TVector3 fwddir  = (fdcpos-vertex).Unit();
    TVector3 posuswk = fdcpos+(250-169.6)/fwddir.Z()*fwddir;
    double dxdz = fwddir.X()/fwddir.Z();
    double dydz = fwddir.Y()/fwddir.Z();
    FillHist("FWD_USWKPosition",posuswk.X(),posuswk.Y());
    FillHist("FWD_Angle",dxdz,dydz);
  }

  TVector3 fdcpos = trackfdc1->GetPosatZ(169.6);
  TVector3 fwddir  = (fdcpos-vertex).Unit();
  TVector3 posuswk = fdcpos+(250-169.6)/fwddir.Z()*fwddir;

  const int nhodo = 3;
  int hodoid[nhodo]={
    CID_CVC,CID_PC,CID_NC
  };
  if(MulNC>=1) MulNC = 1;
  int MulFWD[nhodo] = {MulCVC,MulPC,MulNC};
  // ########################## //
  // FWD Hodo Analysis          //
  // ########################## //
  double leff = 100.0; /* cm */
  for(int ihodo=0; ihodo<nhodo; ihodo++){
    const int cid = hodoid[ihodo];
    const char* name = dlist -> GetName(cid).c_str();
    const int nsegs = dlist -> GetNsegs(cid);
    for( int i=0; i<blMan->nHodo(cid); i++ ){
      if(MulFWD[ihodo]!=1){ break; }
      HodoscopeLikeHit *hit = blMan->Hodoi(cid,i);
      if(hit->CheckRange()){
        int seg = hit->seg();
        double ene[2]; ene[0]=hit->ene(0); ene[1]=hit->ene(1);
        double emean = hit->emean();
        double ctmean = hit->ctmean();
        double tof = ctmean - T0Timing;
        TVector3 hitpos(hit->pos().X(),hit->hitpos(),hit->pos().Z());

        double angle = (posuswk-vertex).Angle((hitpos-posuswk));
        leff = leff/TMath::Cos(fwddir.Y()/fwddir.Z());
        TVector3 tmpfwddir = fwddir;
        tmpfwddir.RotateY(fwddir.Y()/fwddir.Z());
        double alpha = tmpfwddir.X()/tmpfwddir.Z();
        double radius = leff/2.0/TMath::Cos(alpha)/TMath::Tan(angle/2.0);
        double larc = radius*angle;
        double line = leff/TMath::Cos(alpha);
        double diff = line - larc;
        double flength = (posuswk-vertex).Mag()+(hitpos-posuswk).Mag() - diff;

        FillHist("FWD_BendingAngle",angle);
        FillHist("FWD_IncidentAngle",alpha);
        FillHist("FWD_Radius",radius);
        FillHist("FWD_ArcLength",larc);
        FillHist("FWD_LineLength",line);
        FillHist("FWD_LengthDifference",diff);
        FillHist("FWD_FlightLength",flength);

        double tmpmom1=0.0;
        double tmptof1=0.0;
        double mass = piMass;
        //mass = pMass;
        
        ELossTools::CalcElossBeamTGeo(T0Position,vertex,BeamMom,mass,tmpmom1,tmptof1);
        double tmpmom2=0.0;
        double tmptof2=0.0;
        ELossTools::CalcElossForwardTGeo(vertex,fdcpos,flength,tmpmom1,mass,tmpmom2,tmptof2);

        double tof_calc = tmptof1+tmptof2;
        double tofdiff = tof_calc - tof;
        FillHist("FWD_Momentum",tmpmom2);
        FillHist("FWD_TOFCalc",tof_calc);
        FillHist("FWD_TOF",tof);
        FillHist("FWD_TOFDifference",tofdiff);

        const char* ud[2] = {"u","d"};
        for(int iud=0; iud<2; iud++){
          FillHist(Form("%s%s%d_dEvsTimeOffset",name,ud[iud],seg),ene[iud],tofdiff,1);
          FillHist(Form("%s%s%d_dEvsTOF",name,ud[iud],seg),ene[iud],tof,1);
        }
        FillHist(Form("%s%d_TOF",name,seg),tof,1);
        FillHist(Form("%s_TOF",name),tof,1);
        FillHist(Form("%s_TOFvsT0Segment",name),tof,T0Segment,1);
        FillHist(Form("%s_TOFvsT0CTimeMean",name),tof,T0Timing,1);
        FillHist(Form("%s%d_dEMeanvsTOF",name,seg),emean,tof,1);
        FillHist(Form("%s%d_TimeOffset",name,seg),tofdiff,1);
        FillHist(Form("%s_TimeOffset",name),tofdiff,1);
        FillHist(Form("%s_TimeOffsetvsT0Segment",name),tofdiff,T0Segment,1);
        FillHist(Form("%s_TimeOffsetvsT0CTimeMean",name),tofdiff,T0Timing,1);
        FillHist(Form("%s%d_dEMeanvsTimeOffset",name,seg),emean,tofdiff,1);
        FillHist(Form("%s_HitPattern",name),seg,1);
      }
    }
  }
}

bool MyAnalysisUSWKCheck::FillHist(TString name, double val1, int weight)
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

bool MyAnalysisUSWKCheck::FillHist(TString name, TString val1, int weight)
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

bool MyAnalysisUSWKCheck::FillHist(TString name, double val1, double val2, int weight)
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

bool MyAnalysisUSWKCheck::FillHist(TString name, TString val1, TString val2, int weight)
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

bool MyAnalysisUSWKCheck::Initialize(ConfMan* confMan)
{
  std::cout << "### MyAnalysisUSWKCheck::Initialize ###" << std::endl;
  DetectorList *dlist=DetectorList::GetInstance();
  confFile = confMan;
  std::string ofname = confMan->GetOutFileName();
  ofname.insert(ofname.find(".root"),"_anaUSWKCheck");
  rtFile =  new TFile( Form("%s",ofname.c_str()), "RECREATE");
  rtFile -> cd();

  // Chamber
  const int ndc=4;
  int dccid[ndc]={CID_BLC1,CID_BLC2,CID_BPC,CID_FDC1};
  TString dcname[ndc]={"BLC1","BLC2","BPC","FDC1"};
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
  new TH2F( "FWD_USWKPosition", "Position at USWK;X position (cm);Y position (cm)",600,-25,25,600,-25,25);
  new TH2F( "FWD_Angle", "FWD Angle;dx/dz;dy/dz",1000,-0.05,0.05,1000,-0.05,0.05);

  new TH2F( "BLC2_FFPosition", "BLC2 Track Position at FF;X position (cm);Y position (cm)",600,-15,15,600,-15,15);
  new TH2F( "BLC2_75Position", "BLC2 Track Position at z=-75.0 cm;X position (cm);Y position (cm)",600,-15,15,600,-15,15);
  new TH2F( "BPC_FFPosition", "BPC Track Position at FF;X position (cm);Y position (cm)",600,-15,15,600,-15,15);
  new TH2F( "BPC_75Position", "BPC Track Position at z=-75.0 cm;X position (cm);Y position (cm)",600,-15,15,600,-15,15);
  new TH2F( "BPC_USWKPosition", "BPC Track Position at USWK;X position (cm);Y position (cm)",600,-30,30,600,-30,30);
  new TH2F( "FDC1_USWKPosition", "FDC1 Track Position at USWK;X position (cm);Y position (cm)",600,-30,30,600,-30,30);
  new TH2F( "BLC2_Direction", "BLC2 Track Direction;dX/dZ;dY/dZ",1000,-0.05,0.05,1000,-0.05,0.05);
  new TH2F( "BPC_Direction", "BPC Track Direction;dX/dZ;dY/dZ",1000,-0.05,0.05,1000,-0.05,0.05);
  new TH2F( "FDC1_Direction", "FDC1 Track Direction;dX/dZ;dY/dZ",1000,-0.05,0.05,1000,-0.05,0.05);
  new TH2F( "BLC2BPC_Distance", "Distance between BLC2 and BPC;X distance (cm);Y distance (cm)",1000,-5.0,5.0,1000,-5.0,5.0 );
  new TH2F( "BLC2BPC_DistancewCut", "Distance between BLC2 and BPC;X distance (cm);Y distance (cm)",1000,-5.0,5.0,1000,-5.0,5.0 );
  new TH2F( "BLC2BPC_Angle", "Angle between BLC2 and BPC;dx/dz;dy/dz",1000,-0.05,0.05,1000,-0.05,0.05);
  new TH2F( "BLC2BPC_AnglewCut", "Angle between BLC2 and BPC;dx/dz;dy/dz",1000,-0.05,0.05,1000,-0.05,0.05);
  new TH2F( "FDC1BPC_Distance", "Distance between FDC1 and BPC;X distance (cm);Y distance (cm)",1000,-5.0,5.0,1000,-5.0,5.0 );
  new TH2F( "FDC1BPC_DistancewCut", "Distance between FDC1 and BPC;X distance (cm);Y distance (cm)",1000,-5.0,5.0,1000,-5.0,5.0 );
  new TH2F( "FDC1BPC_Angle", "Angle between FDC1 and BPC;dx/dz;dy/dz",1000,-0.05,0.05,1000,-0.05,0.05);
  new TH2F( "FDC1BPC_AnglewCut", "Angle between FDC1 and BPC;dx/dz;dy/dz",1000,-0.05,0.05,1000,-0.05,0.05);

  // FWD analysis
  std::cout << "Define Histograms for FWD analysis" << std::endl;
  new TH1F( "FWD_BendingAngle", "FWD Bending angle;Bending angle;Counts",1500, 0.0, 1.5 );
  new TH1F( "FWD_IncidentAngle", "FWD Incident angle;Bending angle;Counts",1500, 0.0, 1.5 );
  new TH1F( "FWD_Radius", "FWD Radius;Radius (cm);Counts",2000, 0.0, 2000 );
  new TH1F( "FWD_ArcLength", "FWD Arc length;Arc length (cm);Counts",2000, 0.0, 400 );
  new TH1F( "FWD_LineLength", "FWD Line length;Line length (cm);Counts",2000, 0.0, 400 );
  new TH1F( "FWD_LengthDifference", "FWD Length difference between arc and line;Length difference (cm);Counts",2000, 0.0, 100 );
  new TH1F( "FWD_FlightLength", "FWD Flight length;Flight length (cm);Counts",2000, 0.0, 2000 );
  new TH1F( "FWD_Momentum", "FWD Momentum;Momentum (GeV/c);Counts",1000, 0.5, 1.5 );
  new TH1F( "FWD_TOFCalc", "FWD calclated TOF;TOF_{Calc} (ns);Counts",1000, 20.0, 70.0 );
  new TH1F( "FWD_TOF", "FWD measured TOF;TOF_{Meas} (ns);Counts",1000, 20.0, 70.0 );
  new TH1F( "FWD_TOFDifference", "FWD TOF Difference;TOF Diddfrence (ns);Counts",1000, -20.0, 20.0 );

  std::cout << "Define Histograms for Hodoscopes" << std::endl;
  const int nhodo = 3;
  int hodoid[nhodo]={
    CID_CVC,CID_PC,CID_NC
  };
  for(int ihodo=0; ihodo<nhodo; ihodo++){
    const int cid = hodoid[ihodo];
    const char* name = dlist -> GetName(cid).c_str();
    const int nsegs = dlist -> GetNsegs(cid);
    for( int seg=1; seg<=nsegs; seg++ ){
      const char* ud[2] = {"u","d"};
      const char* udname[2] = {"up","down"};
      for(int iud=0; iud<2; iud++){
        new TH2F( Form("%s%s%d_dEvsTimeOffset",name,ud[iud],seg),   Form("dE and TimeOffset corr. %s-%d(%s);dE (MeV);Time (ns)",name,seg,udname[iud]),     200,    0, 100,  2000,    -20, 20 );
        new TH2F( Form("%s%s%d_dEvsTOF",name,ud[iud],seg),   Form("dE and TOF corr. %s-%d(%s);dE (MeV);TOF (ns)",name,seg,udname[iud]),     200,    0, 100,  2000,    20, 70 );
      }
      new TH1F( Form("%s%d_TimeOffset",name,seg),   Form("TimeOffset %s-%d;Time (ns);Counts",name,seg),    2000,    -20, 20 );
      new TH1F( Form("%s%d_TOF",name,seg),   Form("TOF between T0 and %s-%d;Time of flight T0-%s%d (ns);Counts",name,seg,name,seg),    3000,    -50, 100 );
      new TH2F( Form("%s%d_dEMeanvsTOF",name,seg),   Form("dEMean TOF corr. %s-%d;dE (MeV);TOF (ns)",name,seg),     200,    0, 100,  2000,    -100, 100 );
      new TH2F( Form("%s%d_dEMeanvsTimeOffset",name,seg),   Form("dEMean TimeOffset corr. %s-%d;dE (MeV);TimeOffset (ns)",name,seg),     200,    0, 100,  2000,    -20, 20 );
    }
    new TH1F( Form("%s_HitPattern",name), Form("Hit pattern %s;%s Segment;Counts",name,name), nsegs, 0.5, nsegs+0.5);
    new TH1F( Form("%s_TOF",name),   Form("TOF between T0 and %s;Time of flight T0-%s (ns);Counts",name,name),    3000,    -50, 100 );
    new TH2F( Form("%s_TOFvsT0Segment",name),   Form("TOF between T0 and %s vs. T0 segment;Time of flight T0-%s (ns);T0 segment",name,name),    3000,    -50, 100, 5, 0.5, 5.5 );
    new TH2F( Form("%s_TOFvsSegment",name),   Form("TOF between T0 and %s vs. segment;Time of flight T0-%s (ns);T0 segment",name,name),    3000,    -50, 100, 5, 0.5, 5.5 );
    new TH1F( Form("%s_TimeOffset",name),   Form("TimeOffset;TimeOffset (ns);Counts"),    3000,    -20, 20 );
    new TH2F( Form("%s_TimeOffsetvsT0Segment",name),   Form("TimeOffset vs. T0 segment;TimeOffset (ns);T0 segment"),    3000,    -20, 20, 5, 0.5, 5.5 );
    new TH2F( Form("%s_TimeOffsetvsSegment",name),   Form("TimeOffset vs. segment;TimeOffset (ns);Segment"),    3000,    -20, 20, nsegs, 0.5, nsegs+0.5 );
  }

  std::cout << "=== End of [MyAnalysisUSWKCheck::Initialize] === " << std::endl;

  return true;
}
