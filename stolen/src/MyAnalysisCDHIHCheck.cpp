// MyAnalysisCDHIHCheck.cpp

#include "MyAnalysisCDHIHCheck.h"

MyAnalysisCDHIHCheck::MyAnalysisCDHIHCheck(TFile* rt, ConfMan* conf)
{
  Initialize(conf);
  Clear();
}

MyAnalysisCDHIHCheck::~MyAnalysisCDHIHCheck()
{
  Clear();
  rtFile->cd();
  rtFile->Write();
  rtFile->Close();
}

void MyAnalysisCDHIHCheck::Clear()
{
  return;
}

bool MyAnalysisCDHIHCheck::DoAnalysis(ConfMan* conf, EventHeader* header, BeamLineHitMan* blMan, BeamLineTrackMan* bltrackMan, CDSHitMan* cdsMan, CDSTrackingMan* cdstrackMan)
{
  rtFile->cd();
  DetectorList *dlist=DetectorList::GetInstance();

  /* ======================================= */
  /* ========== Trigger selection ========== */
  if(header->IsTrig(Trig_Cosmic)){
    return true;
  }
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
    TVector3 posff  = trackbpc->GetPosatZ(-5.0);
    TVector3 posffblc2 = trackblc2->GetPosatZ(-5.0);
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
  //if(!FIDUCIAL){ return false; }
  TVector3 vertex = trackbpc->GetPosatZ(0.0);

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
    FillHist("CDS_OverbetavsMomentum",1./beta,mom);
    FillHist("CDS_Mass2vsMomentum",mass2,mom);

    double calc_beta,tofvtxcdc;
    if(!TrackTools::FindMass2(cdc,beam,cdhtime-beam->t0time(),calc_beta,mass2,tofvtxcdc)) continue;
    beta = calc_beta;
    FillHist("CDS_OverbetavsMomentum_AfterFindMass",1./beta,mom);
    FillHist("CDS_Mass2vsMomentum_AfterFindMass",mass2,mom);

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

    TVector3 cfrpvtx;
    double param[5];  
    cdc->GetParameters(CID_CDCCFRP,param,cfrpvtx);
    cdc->GetParameters(param);
    TVector3 vtxcdh = cdc->CDHVertex();
    cdc->GetParameters(param);
    double fl = MathTools::CalcHelixArc(param,cfrpvtx,vtxcdh);
    double mom2 = TMath::Abs(cdc->Momentum(CID_CDCCFRP));
    double mass = cdsMass[pid];
    double momout=0,tmptof=0;
    ELossTools::CalcdE(mom2,mass,fl,"CDCGas",momout,-1,tmptof);
    double offset = cdhtime - beam->t0time() - beam_tof - tof - tmptof;
    if(pid==CDS_Other){ continue; }
    //if(pid!=CDS_PiPlus&&pid!=CDS_PiMinus){ continue; }

    HodoscopeLikeHit* cdh = cdc->CDHHit(cdsMan,0);
    FillHist("CDS_OverbetavsEnergy",1./beta,cdh->emean());
    if((pid==CDS_PiPlus||pid==CDS_PiMinus)&&TMath::Abs(mom)>0.2&&(-5.0<vtxcdh.Z()&&vtxcdh.Z()<5.0)){
      FillHist(Form("CDHu%d_dE",cdh->seg()),cdh->eu(),1);
      FillHist(Form("CDHd%d_dE",cdh->seg()),cdh->ed(),1);
      FillHist(Form("CDH%d_dEMean",cdh->seg()),cdh->emean(),1);
    }
    FillHist(Form("CDH_TimeOffset"),offset,1);
    FillHist(Form("CDH_TimeOffsetvsMomentum"),offset,TMath::Abs(mom),1);
    FillHist(Form("CDH_dEMeanvsTimeOffset"),cdh->emean(),offset,1);
    FillHist(Form("CDH%d_TimeOffset",cdh->seg()),offset,1);
    FillHist(Form("CDH%d_TimeOffsetvsMomentum",cdh->seg()),offset,TMath::Abs(mom),1);
    FillHist(Form("CDH%d_dEMeanvsTimeOffset",cdh->seg()),cdh->emean(),offset,1);
    FillHist(Form("CDH%d_CTimeSub",cdh->seg()),cdh->ctsub(),1);
    FillHist(Form("CDHu%d_TimeOffset",cdh->seg()),offset,1);
    FillHist(Form("CDHu%d_dEvsTimeOffset",cdh->seg()),cdh->eu(),offset,1);
    FillHist(Form("CDHu%d_dEMeanvsTimeOffset",cdh->seg()),cdh->emean(),offset,1);
    FillHist(Form("CDHd%d_TimeOffset",cdh->seg()),offset,1);
    FillHist(Form("CDHd%d_dEvsTimeOffset",cdh->seg()),cdh->ed(),offset,1);
    FillHist(Form("CDHd%d_dEMeanvsTimeOffset",cdh->seg()),cdh->emean(),offset,1);
    FillHist(Form("CDH%d_TimeOffset_%s",cdh->seg(),pname[pid].Data()),offset,1);
    FillHist(Form("CDH%d_TimeOffsetvsMomentum_%s",cdh->seg(),pname[pid].Data()),offset,TMath::Abs(mom),1);
    FillHist(Form("CDS_Mass2vsdEMean",pname[pid].Data()),mass2,cdh->emean());

    FillHist(Form("T0%d_TimeOffset",T0Segment),-offset,1);
    FillHist(Form("T0_dEMeanvsTimeOffset"),T0dEMean,-offset,1);
    FillHist(Form("T0u%d_TimeOffset",T0Segment),-offset,1);
    FillHist(Form("T0u%d_dEvsTimeOffset",T0Segment),T0dE[0],-offset,1);
    FillHist(Form("T0u%d_dEMeanvsTimeOffset",T0Segment),T0dEMean,-offset,1);
    FillHist(Form("T0d%d_TimeOffset",T0Segment),-offset,1);
    FillHist(Form("T0d%d_dEvsTimeOffset",T0Segment),T0dE[1],-offset,1);
    FillHist(Form("T0d%d_dEMeanvsTimeOffset",T0Segment),T0dEMean,-offset,1);

    HodoscopeLikeHit* ih = cdc->IHHit(cdsMan,0);
    if((pid==CDS_PiPlus||pid==CDS_PiMinus)&&TMath::Abs(mom)>0.2&&(-5.0<vtxcdh.Z()&&vtxcdh.Z()<5.0)){
      FillHist(Form("IHu%d_dE",ih->seg()),ih->eu(),1);
      FillHist(Form("IHd%d_dE",ih->seg()),ih->ed(),1);
      FillHist(Form("IH%d_dEMean",ih->seg()),ih->emean(),1);
    }
  }

  delete beam;
  return true;
}

bool MyAnalysisCDHIHCheck::FillHist(TString name, double val1, int weight)
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

bool MyAnalysisCDHIHCheck::FillHist(TString name, TString val1, int weight)
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

bool MyAnalysisCDHIHCheck::FillHist(TString name, double val1, double val2, int weight)
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

bool MyAnalysisCDHIHCheck::FillHist(TString name, TString val1, TString val2, int weight)
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

bool MyAnalysisCDHIHCheck::Initialize(ConfMan* confMan)
{
  std::cout << "### MyAnalysisCDHIHCheck::Initialize ###" << std::endl;
  DetectorList *dlist=DetectorList::GetInstance();
  confFile = confMan;
  std::string ofname = confMan->GetOutFileName();
  ofname.insert(ofname.find(".root"),"_anaCDHIHCheck");
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
  new TH1F( "Beam_TOFwithCut", "Beam TOF;TOF between BHD and T0 (ns);Counts",1000, 20.0, 40.0 );
  new TH2F( "Beam_MomentumvsTOF", "Beam Momentum vs. TOF;Beam Momentum (GeV/c);TOF between BHD and T0 (ns)",1000,0.5,1.5,1000, 20.0, 40.0 );
  new TH2F( "Beam_MomentumvsTOFwCut", "Beam Momentum vs. TOF;Beam Momentum (GeV/c);TOF between BHD and T0 (ns)",1000,0.5,1.5,1000, 20.0, 40.0 );
  new TH1F( "Beam_Chi2", "Beam mom #chi^{2};#chi^{2}/ndf;Counts",1000, 0.0, 100.0 );
  new TH1F( "Beam_Chi2wCut", "Beam mom #chi^{2};#chi^{2}/ndf;Counts",1000, 0.0, 100.0 );
  new TH2F( "Beam_FFPosition", "Position at FF;X position (cm);Y position (cm)",600,-15,15,600,-15,15);
  new TH2F( "Beam_FFPositionwCut", "Position at FF;X position (cm);Y position (cm)",600,-15,15,600,-15,15);
  new TH2F( "Beam_Angle", "Beam Angle;dx/dz;dy/dz",1000,-0.05,0.05,1000,-0.05,0.05);
  new TH2F( "Beam_AnglewCut", "Beam Angle;dx/dz;dy/dz",1000,-0.05,0.05,1000,-0.05,0.05);
  new TH2F( "BHD_SegmentvsBeamMomentum", "BHD segment vs. Beam momentum;BHD segment;Beam momentum (GeV/c)", 20, 1, 21, 1000, 0.5, 1.5 );
  new TH2F( "BHD_SegmentvsBeamMomentumwCut", "BHD segment vs. Beam momentum;BHD segment;Beam momentum (GeV/c)", 20, 1, 21, 1000, 0.5, 1.5 );

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

  /* CDS Analysis */
  std::cout << "Define Histograms for CDS Analysis" << std::endl;
  new TH2F( "CDS_Mass2vsMomentum", "Mass^{2} vs. Momentum;Mass^{2} (GeV/c^{2})^{2};Momentum (GeV/c)", 3000, -5, 10, 800, -2, 2 );
  new TH2F( "CDS_Mass2vsMomentum_AfterFindMass", "Mass^{2} vs. Momentum;Mass^{2} (GeV/c^{2})^{2};Momentum (GeV/c)", 3000, -5, 10, 800, -2, 2 );
  new TH2F( "CDS_Mass2vsMomentum_AfterELossCorr", "Mass^{2} vs. Momentum;Mass^{2} (GeV/c^{2})^{2};Momentum (GeV/c)", 3000, -5, 10, 800, -2, 2 );
  new TH2F( "CDS_Mass2vsdEMean", "Mass^{2} vs. dEMean;Mass^{2} (GeV/c^{2})^{2};dE Mean (MeVee)", 3000, -5, 10, 1000, 0, 100 );
  new TH2F( "CDS_OverbetavsMomentum", "1/#beta vs. Momentum;1/#beta;Momentum (GeV/c)", 1000, 0, 10, 800, -2, 2 );
  new TH2F( "CDS_OverbetavsMomentum_AfterFindMass", "1/#beta vs. Momentum;1/#beta;Momentum (GeV/c)", 1000, 0, 10, 800, -2, 2 );
  new TH2F( "CDS_TOFvsMomentum", "TOF vs. Momentum;TOF (ns);Momentum (GeV/c)", 1000, 0, 40, 800, -2, 2 );
  new TH2F( "CDS_OverbetavsMomentum_Calc", "1/#beta vs. Momentum;1/#beta;Momentum (GeV/c)", 1000, 0, 10, 800, -2, 2 );
  new TH2F( "CDS_TOFvsMomentum_Calc", "TOF vs. Momentum;TOF (ns);Momentum (GeV/c)", 1000, 0, 40, 800, -2, 2 );
  new TH2F( "CDS_OverbetavsEnergy", "1/#beta vs. Momentum;1/#beta;dE (MeV)", 1000, 0, 10, 800, 80, 0 );


  /* Hodoscope */
  std::cout << "Define Histograms for Hodoscopes" << std::endl;
  TString pname[8] = {"pip","p","d","t","he","pim","k","o"};
  const int nhodo = 3;
  int hodoid[nhodo]={
    CID_T0,CID_CDH,CID_IH
  };
  for(int ihodo=0; ihodo<2; ihodo++){
    const int cid = hodoid[ihodo];
    const char* name = dlist -> GetName(cid).c_str();
    const int nsegs = dlist -> GetNsegs(cid);
    for( int seg=1; seg<=nsegs; seg++ ){
      const char* ud[2] = {"u","d"};
      const char* udname[2] = {"up","down"};
      for(int iud=0; iud<2; iud++){
        if(cid==CID_CDH){
          new TH1F( Form("%s%s%d_dE",name,ud[iud],seg),   Form("dE %s-%d(%s);dE (MeV);Counts",name,seg,udname[iud]),    800,    0, 40 );
        }
        else if(cid==CID_IH){
          new TH1F( Form("%s%s%d_dE",name,ud[iud],seg),   Form("dE %s-%d(%s);dE (MeV);Counts",name,seg,udname[iud]),    800,    0, 4 );
        }
        new TH1F( Form("%s%s%d_TimeOffset",name,ud[iud],seg),   Form("Time Offset %s-%d(%s);Time (ns);Counts",name,seg,udname[iud]),    400,    -10, 10 );
        if(cid!=CID_T0){
          new TH2F( Form("%s%s%d_dEvsTimeOffset",name,ud[iud],seg),   Form("dE Time Offset corr. %s-%d(%s);dE (MeV);Time (ns)",name,seg,udname[iud]),     200,    0, 100,  200,    -1, 1 );
          new TH2F( Form("%s%s%d_dEMeanvsTimeOffset",name,ud[iud],seg),   Form("dEMean Time Offset corr. %s-%d(%s);dE (MeV);Time (ns)",name,seg,udname[iud]),     100,    0, 100,  400,    -10, 10 );
        }
        else{
          new TH2F( Form("%s%s%d_dEvsTimeOffset",name,ud[iud],seg),   Form("dE Time Offset corr. %s-%d(%s);dE (MeV);Time (ns)",name,seg,udname[iud]),     1000,    0, 10,  1000,    -5, 5 );
          new TH2F( Form("%s%s%d_dEMeanvsTimeOffset",name,ud[iud],seg),   Form("dEMean Time Offset corr. %s-%d(%s);dE (MeV);Time (ns)",name,seg,udname[iud]),     1000,    0, 10,  400,    -10, 10 );
        }
      }
      if(cid==CID_CDH){
        new TH1F( Form("%s%d_dEMean",name,seg),   Form("Mean dE %s-%d;dE (MeV);Counts",name,seg),    800,    0, 40 );
      }
      else if(cid==CID_IH){
        new TH1F( Form("%s%d_dEMean",name,seg),   Form("Mean dE %s-%d;dE (MeV);Counts",name,seg),    800,    0, 4 );
      }
      new TH1F( Form("%s%d_TimeOffset",name,seg),   Form("Time Offset %s-%d;Time (ns);Counts",name,seg),    1000,    -10, 10 );
      new TH2F( Form("%s%d_TimeOffsetvsMomentum",name,seg),   Form("Time Offset %s-%d vs. Momentum;Time (ns);Momentum (GeV/c)",name,seg),    1000,    -10, 10, 800, -2, 2 );
      new TH2F( Form("%s%d_dEMeanvsTimeOffset",name,seg),   Form("dEMean Time Offset corr. %s-%d;dE (MeV);Time (ns)",name,seg),     1000,    0, 100,  400,    -10, 10 );
      new TH1F( Form("%s%d_CTimeSub",name,seg),   Form("CTime Subtract %s-%d;Time (ns);Counts",name,seg),    1000,    -10, 10 );
      new TH1F( Form("%s%d_Multiplicity",name,seg),   Form("Multiplicity %s-%d;Multiplicity;Counts",name,seg),    11,    0, 10 );
      for(int ip=0; ip<8; ip++){
        new TH1F( Form("%s%d_TimeOffset_%s",name,seg,pname[ip].Data()),   Form("Time Offset %s-%d;Time (ns);Counts",name,seg),    5000,    -5, 5 );
        new TH2F( Form("%s%d_TimeOffsetvsMomentum_%s",name,seg,pname[ip].Data()),   Form("Time Offset %s-%d vs. Momentum;Time (ns);Momentum (GeV/c)",name,seg),    1000,    -5, 5, 800, -2, 2 );
      }
    }
    new TH1F( Form("%s_TimeOffset",name),   Form("Time Offset %s;Time (ns);Counts",name),    5000,    -5, 5 );
  }

  std::cout << "=== End of [MyAnalysisCDHIHCheck::Initialize] === " << std::endl;

  return true;
}
