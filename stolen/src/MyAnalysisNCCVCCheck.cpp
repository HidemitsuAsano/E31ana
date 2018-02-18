// MyAnalysisNCCVCCheck.cpp

#include "MyAnalysisNCCVCCheck.h"

MyAnalysisNCCVCCheck::MyAnalysisNCCVCCheck(TFile* rt, ConfMan* conf)
{
  Initialize(conf);
  Clear();
}

MyAnalysisNCCVCCheck::~MyAnalysisNCCVCCheck()
{
  Clear();
  rtFile->cd();
  rtFile->Write();
  rtFile->Close();
}

void MyAnalysisNCCVCCheck::Clear()
{
  return;
}

bool MyAnalysisNCCVCCheck::DoAnalysis(ConfMan* conf, EventHeader* header, BeamLineHitMan* blMan, BeamLineTrackMan* bltrackMan, CDSHitMan* cdsMan, CDSTrackingMan* cdstrackMan)
{
  rtFile->cd();
  DetectorList *dlist=DetectorList::GetInstance();

  /* ======================================= */
  /* ========== Trigger selection ========== */
  if(!header->IsTrig(Trig_Kaon)){
    return true;
  }
  if(!header->IsTrig(Trig_KCDH1)){
    return true;
  }
  /* ========== Trigger selection ========== */
  /* ======================================= */

  /* ===================================== */
  /* ========== Event selection ========== */
  double T0Timing = -9999.; 
  double T0dEMean = -9999.;
  double T0dE[2] = {-9999.,-9999.};
  int T0Segment = -1;
  int MulT0 = 0;
  for(int i=0; i<blMan->nT0(); i++){
    if(blMan->T0(i)->CheckRange()){
      MulT0++;
      T0Timing = blMan->T0(i)->ctmean();
      T0Segment = blMan->T0(i)->seg();
      T0dEMean = blMan->T0(i)->ene(0);
      T0dE[0] = blMan->T0(i)->ene(0);
      T0dE[1] = blMan->T0(i)->ene(1);
    }
  }
  if(MulT0!=1){
    return true;
  }
  if(T0Segment!=3){
    return true;
  }
  int MulBHD = 0;
  int BHDSegment = -1;
  double BHDTiming = -9999.;
  double BHDX = -9999.;
  for(int i=0; i<blMan->nBHD(); i++){
    if(blMan->BHD(i)->CheckRange()){
      if(6<=blMan->BHD(i)->seg()&&blMan->BHD(i)->seg()<=15){
        MulBHD++;
        BHDTiming = blMan->BHD(i)->ctmean();
        BHDSegment = blMan->BHD(i)->seg();
        TVector3 tmppos; conf->GetGeomMapManager()->GetGPos(CID_BHD,BHDSegment,tmppos);
        BHDX = tmppos.X();
      }
    }
  }
  if(MulBHD!=1){
    return true;
  }
  double BeamTOF = T0Timing-BHDTiming;
  FillHist(Form("Beam_TOF"), BeamTOF);
  //if((T0Timing-BHDTiming)<27.8||29.5<(T0Timing-BHDTiming)){ /* Kaon */
  if((BeamTOF<27.8)||(29.5<BeamTOF)){ /* Kaon */
    return true;
  }
  FillHist(Form("Beam_TOFwithCut"), BeamTOF);

  int MulCDH=0;
  for(int i=0; i<cdsMan->nCDH(); i++){
    HodoscopeLikeHit* hit = cdsMan->CDH(i);
    if(hit->CheckRange()) MulCDH++;
  }
  int MulIH=0;
  for(int i=0; i<cdsMan->nIH(); i++){
    HodoscopeLikeHit* hit = cdsMan->IH(i);
    if(hit->CheckRange()) MulIH++;
  }
  if(MulCDH==0||MulIH==0) return true;
  if(cdstrackMan->nGoodTrack()==0) return true;
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
  if(MulBVC!=0||MulCVC!=0) return true;
  int MulNC=0;
  for(int i=0; i<blMan->nNC(); i++){
    HodoscopeLikeHit* hit = blMan->NC(i);
    if(hit->CheckRange()) MulNC++;
  }
  if(MulNC==0) return true;

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
    for(int i=0; i<bltrackMan->ntrackBLDC(dccid[idc]); i++){
      if(i>=1) break;
      LocalTrack* track = bltrackMan->trackBLDC(dccid[idc],i);
      double tracktime = track -> GetTrackTime();
      double chi2 = track -> chi2all();
      if(ll[idc]<tracktime&&tracktime<ul[idc]){
        if(chi2<chiul[idc]){
          nhitdc[idc]++;
          trackid[idc] = i;
        }
      }
    }
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
  delete D5;
  /* ## Beam chisquare < 20 request ## */
  if(!(BeamChi2<20)) return false;

  // ==== BPC analysis ==== // 
  idc=2;
  dcname = dlist->GetName(dccid[idc]).data();
  for(int i=0; i<bltrackMan->ntrackBLDC(dccid[idc]); i++){
    if(i>=1) break;
    LocalTrack* track = bltrackMan->trackBLDC(dccid[idc],i);
    double tracktime = track -> GetTrackTime();
    double chi2 = track -> chi2all();
    if(ll[idc]<tracktime&&tracktime<ul[idc]){
      if(chi2<chiul[idc]){
        nhitdc[idc]++;
        trackid[idc] = i;
      }
    }
  }
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
    TVector3 posffblc2 = trackblc2->GetPosatZ(0.0);
    TVector3 posffbpc  = trackbpc->GetPosatZ(0.0);
    /* ## BLC2 and BPC track maching ## */
    if( !((xll<dist.X()&&dist.X()<xul) && (yll<dist.Y()&&dist.Y()<yul)
          && (dxll<dxdz&&dxdz<dxul) && (dyll<dydz&&dydz<dyul)) )
      return false;
  }
  // ==== Beam position at FF ==== // 
  bool FIDUCIAL = false;
  {
    TVector3 posff  = trackbpc->GetPosatZ(0.0);
    TVector3 posffblc2 = trackblc2->GetPosatZ(0.0);
    double dxdz = BeamP.X()/BeamP.Z();
    double dydz = BeamP.Y()/BeamP.Z();
    /* ## On target selection ## */
    if(GeomTools::GetID(posff)==CID_Fiducial){
      FIDUCIAL = true;
    }
  }
  if(!FIDUCIAL){ return false; }

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
  beam->SetBLCPos(trackblc2->GetPosatZ(0.0));
  TVector3 tmpdir=trackblc2->GetMomDir();
  beam->SetBLCDir(tmpdir);
  beam->SetMomentum(BeamMom);

  double tmpdis = 9999.;
  TVector3 vertex;
  for(int id=0; id<cdstrackMan->nGoodTrack(); id++){
    CDSTrack* cdc=cdstrackMan->GoodTrack(id);
    if(cdc->nCDHHit()!=1) continue;
    if(cdc->nIHHit()!=1) continue;
    double cdhtime=0,dis=0;
    int cdhseg=0;
    if(!cdc->GetCDHHit(cdsMan,cdhseg,cdhtime)) continue;
    if(cdc->Chi()>30) continue;

    TVector3 cdsvtx,beamvtx;
    double par[5];
    cdc->GetParameters(CID_CDC,par,cdsvtx);
    double mom=cdc->Momentum();
    if(!cdc->GetVertex(beam->t0pos(),beam->bpcdir(),beamvtx,cdsvtx)) continue;
    dis = (beamvtx-cdsvtx).Mag();

    double beam_out=0, beam_tof=0;
    ELossTools::CalcElossBeamTGeo(beam->t0pos(),beamvtx,beam->mom(),kpMass,beam_out,beam_tof);
    TVector3 cdhvtx=cdc->CDHVertex();
    double cdc_dis=MathTools::CalcHelixArc(par,cdhvtx,cdsvtx);
    double beta=cdc_dis/(cdhtime-beam->t0time()-beam_tof)/(Const*100);
    double mass2=mom*mom*(1/(beta*beta)-1);

    double calc_beta,tofvtxcdc;
    if(!TrackTools::FindMass2(cdc,beam,cdhtime-beam->t0time(),calc_beta,mass2,tofvtxcdc)) continue;
    beta = calc_beta;

    int pid=TrackTools::PID(mom,mass2);
    //pid=TrackTools::PIDcorr3(mom,mass2);
    cdc->SetPID(pid);

    TVector3 cdsvtx2,beamvtx2;
    double tof=0,tmpl=0;
    if(!cdc->CalcVertexTimeLength(beam->t0pos(),beam->bpcdir(),cdsMass[pid],cdsvtx2,beamvtx2,tof,tmpl,true)) continue;
    FillHist("CDS_OverbetavsMomentum_AfterELossCorr",1./beta,mom);
    FillHist("CDS_Mass2vsMomentum_AfterELossCorr",mass2,mom);
    TString pname[8] = {"pip","p","d","t","he","pim","k","o"};
    TVector3 vtxb = beamvtx2;
    double vbdis = (beamvtx2-cdsvtx2).Mag();
    if(tmpdis>vbdis){
      if(GeomTools::GetID(vtxb)==CID_Fiducial){
        tmpdis = vbdis;
        vertex = vtxb;
      }
    }
  }
  FillHist("BVDIS",tmpdis);
  if(tmpdis>10){ return true; }

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
      FillHist(Form("%s_Overbeta",name),1.0/beta);
      FillHist(Form("%s%d_Overbeta",name,hit->seg()),1.0/beta);
    }

    // For T0
    FillHist(Form("T0%d_TOF",T0Segment),tof);
    FillHist(Form("T0%d_dEvsTOF",T0Segment),T0dEMean,tof);
    FillHist(Form("T0%d_TimeOffset",T0Segment),-offset);
    FillHist(Form("T0%d_dEvsTimeOffset",T0Segment),T0dEMean,-offset);

    FillHist(Form("T0u%d_TOF",T0Segment),tof);
    FillHist(Form("T0u%d_dEvsTOF",T0Segment),T0dE[0],tof);
    FillHist(Form("T0d%d_TOF",T0Segment),tof);
    FillHist(Form("T0d%d_dEvsTOF",T0Segment),T0dE[1],tof);
    FillHist(Form("T0u%d_TimeOffset",T0Segment),-offset);
    FillHist(Form("T0u%d_dEvsTimeOffset",T0Segment),T0dE[0],-offset);
    FillHist(Form("T0d%d_TimeOffset",T0Segment),-offset);
    FillHist(Form("T0d%d_dEvsTimeOffset",T0Segment),T0dE[1],-offset);


  }

  delete beam;
  return true;
}

bool MyAnalysisNCCVCCheck::FillHist(TString name, double val1, int weight)
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

bool MyAnalysisNCCVCCheck::FillHist(TString name, TString val1, int weight)
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

bool MyAnalysisNCCVCCheck::FillHist(TString name, double val1, double val2, int weight)
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

bool MyAnalysisNCCVCCheck::FillHist(TString name, TString val1, TString val2, int weight)
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

bool MyAnalysisNCCVCCheck::Initialize(ConfMan* confMan)
{
  std::cout << "### MyAnalysisNCCVCCheck::Initialize ###" << std::endl;
  DetectorList *dlist=DetectorList::GetInstance();
  confFile = confMan;
  std::string ofname = confMan->GetOutFileName();
  ofname.insert(ofname.find(".root"),"_anaNCCVCCheck");
  rtFile =  new TFile( Form("%s",ofname.c_str()), "RECREATE");
  rtFile -> cd();

  // CDS Particles
  std::cout << "Define Histograms for CDS particles" << std::endl;
  new TH1F( Form("VBDIS"), Form("VBDIS ;VBDIS (cm);Counts"), 3000, 0.0, 30.0);

  // FWD neutral Particles
  std::cout << "Define Histograms for FWD neutral particles" << std::endl;
  new TH2F( "FWD_OverbetavsdE", "1/#beta vs. Energy deposit;1/#beta;Energy deposit (MeVee)", 3000, 0, 5, 1500, 0, 150);
  new TH2F( "FWD_TOFvsdE", "TOF vs. Energy deposit;TOF (ns);Energy deposit (MeVee)", 1000, 40, 100, 1500, 0, 150 );
  new TH1F( "FWD_Overbeta", "1/#beta;1/#beta;Counts", 3000, 0, 5);
  new TH1F( "FWD_TOF", "TOF;TOF (ns);Counts", 1000, 40, 100);
  new TH2F( "FWD_HitPosition", "Hit position at NC (X vs. Y);X position (cm);Y position (cm)",40,-400,400,25,-250,250);
  new TH2F( "FWD_HitSegment", "NC segment;Segment in layer;Layer",18,-0.5,17.5,9,-0.5,8.5);

  /* Hodoscope */
  std::cout << "Define Histograms for Hodoscopes" << std::endl;
  const int nhodo = 3;
  int hodoid[nhodo]={
    CID_NC,CID_CVC,CID_T0
  };
  for(int ihodo=0; ihodo<nhodo; ihodo++){
    const int cid = hodoid[ihodo];
    const char* name = dlist -> GetName(cid).c_str();
    const int nsegs = dlist -> GetNsegs(cid);
    if(cid==CID_CVC){ continue; }
    for( int seg=1; seg<=nsegs; seg++ ){
      const char* ud[2] = {"u","d"};
      const char* udname[2] = {"up","down"};
      for(int iud=0; iud<2; iud++){
        new TH1F( Form("%s%s%d_TOF",name,ud[iud],seg),   Form("TOF %s-%d(%s);Time (ns);Counts",name,seg,udname[iud]),    1000,    40, 100 );
        new TH1F( Form("%s%s%d_TimeOffset",name,ud[iud],seg),   Form("Time Offset %s-%d(%s);Time (ns);Counts",name,seg,udname[iud]),    1000,    -20, 20 );
        if(cid==CID_T0){
          new TH2F( Form("%s%s%d_dEvsTimeOffset",name,ud[iud],seg),   Form("dE Time Offset corr. %s-%d(%s);dE (MeV);Time (ns)",name,seg,udname[iud]),   200,    0, 10,  100,    -1, 1 );
        }
        else{
          new TH2F( Form("%s%s%d_dEvsTimeOffset",name,ud[iud],seg),   Form("dE Time Offset corr. %s-%d(%s);dE (MeV);Time (ns)",name,seg,udname[iud]),   200,    0, 100,  100,   -1, 1 );
        }
      }
      new TH1F( Form("%s%d_TimeOffset",name,seg),   Form("Time Offset %s-%d;Time (ns);Counts",name,seg),    1000,    -20, 20 );
      if(cid==CID_T0){
        new TH2F( Form("%s%d_dEvsTimeOffset",name,seg),   Form("dEMean Time Offset corr. %s-%d;dE (MeV);Time (ns)",name,seg),     500,    0, 10,  100,    -1, 1 );
      }
      else{
        new TH2F( Form("%s%d_dEvsTimeOffset",name,seg),   Form("dEMean Time Offset corr. %s-%d;dE (MeV);Time (ns)",name,seg),     500,    0, 100,  100,   -1, 1 );
      }
      new TH1F( Form("%s%d_CTimeSub",name,seg),   Form("Subtract CTime %s-%d;Subtract Time (ns);Counts",name,seg),    1000,    -10, 10 );
      new TH1F( Form("%s%d_Overbeta",name,seg), "1/#beta;1/#beta;Counts", 3000, 0, 5);
    }
    new TH1F( Form("%s_Overbeta",name), "1/#beta;1/#beta;Counts", 3000, 0, 5);
    new TH1F( Form("%s_TOF",name),   Form("TOF %s;Time (ns);Counts",name),    1000,    40, 100 );
    new TH1F( Form("%s_TimeOffset",name),   Form("Time Offset %s;Time (ns);Counts",name),    1000,    -20, 20 );
    new TH2F( Form("%s_dEvsTimeOffset",name),   Form("dEMean Time Offset corr. %s;dE (MeV);Time (ns)",name),     500,    0, 100,  100,    -1, 1 );
  }

  std::cout << "=== End of [MyAnalysisNCCVCCheck::Initialize] === " << std::endl;

  return true;
}
