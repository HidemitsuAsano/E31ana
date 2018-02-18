// MyAnalysisCDCBPCCheck.cpp

#include "MyAnalysisCDCBPCCheck.h"

MyAnalysisCDCBPCCheck::MyAnalysisCDCBPCCheck(TFile* rt, ConfMan* conf)
{
  Initialize(conf);
  Clear();
}

MyAnalysisCDCBPCCheck::~MyAnalysisCDCBPCCheck()
{
  Clear();
  rtFile->cd();
  rtFile->Write();
  rtFile->Close();
}

void MyAnalysisCDCBPCCheck::Clear()
{
  return;
}

bool MyAnalysisCDCBPCCheck::DoAnalysis(ConfMan* conf, EventHeader* header, BeamLineHitMan* blMan, BeamLineTrackMan* bltrackMan, CDSHitMan* cdsMan, CDSTrackingMan* cdstrackMan)
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
  if(!header->IsTrig(Trig_KCDH2)){
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
  //if(T0Segment!=3){
  //  return true;
  //}
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
    //if( !((xll<dist.X()&&dist.X()<xul) && (yll<dist.Y()&&dist.Y()<yul)
    //      && (dxll<dxdz&&dxdz<dxul) && (dyll<dydz&&dydz<dyul)) )
    //  return false;
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

  Particle* particle = new Particle();

  // ###################### //
  // CDS all track analysis //
  // ###################### //
  for( int id=0; id<cdstrackMan->nGoodTrack(); id++ ){
    CDSTrack* cdc=cdstrackMan->GoodTrack(id);

    if(cdc->Chi()>30) continue;
    if(cdc->nCDHHit()!=1) continue;
    double cdhtime=0,dis=0;
    int cdhseg=0;
    if(!cdc->GetCDHHit(cdsMan,cdhseg,cdhtime)) continue;
    HodoscopeLikeHit* cdh=0;
    for(int i=0; i<cdsMan->nCDH(); i++){
      HodoscopeLikeHit* hit = cdsMan->CDH(i);
      if(hit->CheckRange()){
        if(hit->seg()==cdhseg) cdh=hit;
      }
    }
    double de = cdh->emean();

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
    FillHist("CDS_dEvsMomentum",de,mom);

    double calc_beta,tofvtxcdc;
    if(!TrackTools::FindMass2(cdc,beam,cdhtime-beam->t0time(),calc_beta,mass2,tofvtxcdc)) continue;
    beta = calc_beta;
    FillHist("CDS_OverbetavsMomentum_AfterFindMass",1./beta,mom);
    FillHist("CDS_Mass2vsMomentum_AfterFindMass",mass2,mom);
    FillHist("CDS_dEvsMomentum_AfterFindMass",de,mom);

    int pid=TrackTools::PID(mom,mass2);
    pid=TrackTools::PIDcorr3(mom,mass2);
    cdc->SetPID(pid);

    TVector3 cdsvtx2,beamvtx2;
    double tof=0,tmpl=0;
    if(!cdc->CalcVertexTimeLength(beam->t0pos(),beam->bpcdir(),cdsMass[pid],beamvtx2,cdsvtx2,tof,tmpl,true)) continue;
    TVector3 Pp;
    if(!cdc->GetMomentum(cdsvtx2,Pp,true,true)) continue;
    mom = mom/TMath::Abs(mom)*Pp.Mag();
    FillHist("CDS_OverbetavsMomentum_AfterELossCorr",1./beta,mom);
    FillHist("CDS_Mass2vsMomentum_AfterELossCorr",mass2,mom);

    /* Set pCDS */
    pCDS* cdstrack = new pCDS();
    cdstrack->SetTrackID(id);
    cdstrack->SetPID(pid);
    cdstrack->SetMomentum(mom);
    cdstrack->SetMass(TMath::Sqrt(mass2));
    cdstrack->SetMass2(mass2);
    cdstrack->SetBeta(beta);
    cdstrack->SetTOF(cdhtime-beam->t0time());
    cdstrack->SetFL(tmpl);
    cdstrack->SetChi2(cdc->Chi());
    cdstrack->SetVertexCDC(cdsvtx2);
    cdstrack->SetVertexBeam(beamvtx2);
    cdstrack->SetMomDir(Pp);
    cdstrack->SetVBDis((beamvtx2-cdsvtx2).Mag());
    for(int icdh=0; icdh<cdc->nCDHHit(); icdh++){
      cdstrack->SetCDHSeg(cdc->CDHHit(cdsMan,icdh)->seg());
    }
    for(int iih=0; iih<cdc->nIHHit(); iih++){
      cdstrack->SetIHSeg(cdc->IHHit(cdsMan,iih)->seg());
    }

    particle->AddCDS(*cdstrack);
  }

  for(int it=0; it<particle->nCDS(); it++){
    pCDS* cdstrack = particle->cdsi(it);
    double mass2 = cdstrack->mass2();
    double beta = cdstrack->beta();
    double mom = cdstrack->mom();
    TVector3 vtxb = cdstrack->vbeam();
    TVector3 vtxh = cdstrack->vcdc();
    FillHist("Vertex_XY", vtxb.X(), vtxb.Y());
    FillHist("Vertex_ZX", vtxb.Z(), vtxb.X());
    FillHist("Vertex_ZY", vtxb.Z(), vtxb.Y());
    FillHist("Vertex_Z", vtxb.Z());
    if(-0.5<vtxb.Z()&&vtxb.Z()<0.5){
      FillHist("Vertex_XY_ZCut", vtxb.X(), vtxb.Y());
    }
    if(-0.5<vtxb.X()&&vtxb.X()<0.5){
      FillHist("Vertex_ZY_XCut", vtxb.Z(), vtxb.Y());
    }
    if(-0.5<vtxb.Y()&&vtxb.Y()<0.5){
      FillHist("Vertex_ZX_YCut", vtxb.Z(), vtxb.X());
    }
    if(-0.5<vtxb.X()&&vtxb.X()<0.5&&-0.5<vtxb.Y()&&vtxb.Y()<0.5){
      FillHist("Vertex_Z_XYCut", vtxb.Z());
    }

    if(GeomTools::GetID(vtxb)==CID_Fiducial){
      FillHist("Vertex_XY_Fiducial", vtxb.X(), vtxb.Y());
      FillHist("Vertex_ZX_Fiducial", vtxb.Z(), vtxb.X());
      FillHist("Vertex_ZY_Fiducial", vtxb.Z(), vtxb.Y());
      FillHist("Vertex_Z_Fiducial", vtxb.Z());
      if(-0.5<vtxb.Z()&&vtxb.Z()<0.5){
        FillHist("Vertex_XY_FiducialZCut", vtxb.X(), vtxb.Y());
      }
      if(-0.5<vtxb.X()&&vtxb.X()<0.5){
        FillHist("Vertex_ZY_FiducialXCut", vtxb.Z(), vtxb.Y());
      }
      if(-0.5<vtxb.Y()&&vtxb.Y()<0.5){
        FillHist("Vertex_ZX_FiducialYCut", vtxb.Z(), vtxb.X());
      }
      if(-0.5<vtxb.X()&&vtxb.X()<0.5&&-0.5<vtxb.Y()&&vtxb.Y()<0.5){
        FillHist("Vertex_Z_FiducialXYCut", vtxb.Z());
      }
      FillHist("CDS_OverbetavsMomentum_Fiducial",1./beta,mom);
      FillHist("CDS_Mass2vsMomentum_Fiducial",mass2,mom);
    }
    if(GeomTools::GetID(vtxb)==CID_TarCell){
      FillHist("Vertex_XY_TargetCell", vtxb.X(), vtxb.Y());
      FillHist("Vertex_ZX_TargetCell", vtxb.Z(), vtxb.X());
      FillHist("Vertex_ZY_TargetCell", vtxb.Z(), vtxb.Y());
    }
    if(GeomTools::GetID(vtxb)==CID_CellRing){
      FillHist("Vertex_XY_CellRing", vtxb.X(), vtxb.Y());
      FillHist("Vertex_ZX_CellRing", vtxb.Z(), vtxb.X());
      FillHist("Vertex_ZY_CellRing", vtxb.Z(), vtxb.Y());
      if(-0.5<vtxb.Z()&&vtxb.Z()<0.5){
        FillHist("Vertex_XY_CellRingZCut", vtxb.X(), vtxb.Y());
      }
      if(-0.5<vtxb.X()&&vtxb.X()<0.5){
        FillHist("Vertex_ZY_CellRingXCut", vtxb.Z(), vtxb.Y());
      }
      if(-0.5<vtxb.Y()&&vtxb.Y()<0.5){
        FillHist("Vertex_ZX_CellRingYCut", vtxb.Z(), vtxb.X());
      }
    }
    if(GeomTools::GetID(vtxb)==CID_CellTube){
      FillHist("Vertex_XY_CellTube", vtxb.X(), vtxb.Y());
      FillHist("Vertex_ZX_CellTube", vtxb.Z(), vtxb.X());
      FillHist("Vertex_ZY_CellTube", vtxb.Z(), vtxb.Y());
      if(-0.5<vtxb.Z()&&vtxb.Z()<0.5){
        FillHist("Vertex_XY_CellTubeZCut", vtxb.X(), vtxb.Y());
      }
      if(-0.5<vtxb.X()&&vtxb.X()<0.5){
        FillHist("Vertex_ZY_CellTubeXCut", vtxb.Z(), vtxb.Y());
      }
      if(-0.5<vtxb.Y()&&vtxb.Y()<0.5){
        FillHist("Vertex_ZX_CellTubeYCut", vtxb.Z(), vtxb.X());
      }
    }
    if(GeomTools::GetID(vtxb)==CID_CellAlBe){
      FillHist("Vertex_XY_CellAlBe", vtxb.X(), vtxb.Y());
      FillHist("Vertex_ZX_CellAlBe", vtxb.Z(), vtxb.X());
      FillHist("Vertex_ZY_CellAlBe", vtxb.Z(), vtxb.Y());
      if(-0.5<vtxb.Z()&&vtxb.Z()<0.5){
        FillHist("Vertex_XY_CellAlBeZCut", vtxb.X(), vtxb.Y());
      }
      if(-0.5<vtxb.X()&&vtxb.X()<0.5){
        FillHist("Vertex_ZY_CellAlBeXCut", vtxb.Z(), vtxb.Y());
      }
      if(-0.5<vtxb.Y()&&vtxb.Y()<0.5){
        FillHist("Vertex_ZX_CellAlBeYCut", vtxb.Z(), vtxb.X());
      }
    }
  }

  // ###################### //
  // CDS 2 track analysis //
  // ###################### //
  for( int id1=0; id1<particle->nCDS(); id1++ ){
    for( int id2=id1+1; id2<particle->nCDS(); id2++ ){
      pCDS* cds1 = particle->cdsi(id1);
      pCDS* cds2 = particle->cdsi(id2);
      CDSTrack* cdc1=cdstrackMan->GoodTrack(cds1->id());
      CDSTrack* cdc2=cdstrackMan->GoodTrack(cds2->id());
      TVector3 vtx=DEFVECT,vtx1=DEFVECT,vtx2=DEFVECT;
      if(!TrackTools::Calc2HelixVertex(cdc1,cdc2,vtx1,vtx2)) continue;
      vtx = (vtx1+vtx2)*0.5;

      TVector3 Pp1=DEFVECT,Pp2=DEFVECT;
      if(!cdc1->GetMomentum(vtx1,Pp1,true,true)) continue;
      if(!cdc2->GetMomentum(vtx2,Pp2,true,true)) continue;
      TVector3 Pp = Pp1+Pp2;
      TLorentzVector L1; L1.SetVectM(Pp1,cds1->pdgmass());
      TLorentzVector L2; L2.SetVectM(Pp2,cds2->pdgmass());
      double im = (L1+L2).M();
      double dist=0, dltmp=0;
      TVector3 xest=DEFVECT, nest=DEFVECT;
      MathTools::LineToLine(beam->bpcpos(),beam->bpcdir(),vtx,Pp.Unit(),dltmp,dist,xest,nest);

      pCDS* pro = new pCDS();
      pro->SetDaughterID1(cds1->id());
      pro->SetDaughterID2(cds2->id());
      pro->SetCombID((int)(pow(2,cds1->pid())+pow(2,cds2->pid())));
      pro->SetMomentum(Pp.Mag());
      pro->SetMass(im);
      pro->SetVertex(vtx);
      pro->SetVertexBeam(xest);
      pro->SetMomDir(Pp.Unit());
      pro->SetVDis((vtx1-vtx2).Mag());
      pro->SetVBDis((xest-vtx).Mag());
      pro->SetPBDCA(dist);
      pro->SetOA(Pp1.Angle(Pp2));
      particle->AddProduct(*pro);
    }
  }

  for(int it=0; it<particle->nProduct(); it++){
    pCDS* cds = particle->product(it);
    TVector3 vtxb = cds->vbeam();
    TVector3 vertex_cdc = cds->vertex();
    double bpcz = vertex_cdc.Z();
    double bpcx=0,bpcy=0;
    beam->bpcpos(bpcz,bpcx,bpcy);
    TVector3 vertex_bpc(bpcx,bpcy,bpcz);
    double diffx = vertex_bpc.X()-vertex_cdc.X();
    double diffy = vertex_bpc.Y()-vertex_cdc.Y();
    TVector3 vertex_diff(diffx,diffy,bpcz);
  
    FillHist("VertexCDC_XY", vertex_cdc.X(), vertex_cdc.Y());
    FillHist("VertexCDC_ZX", vertex_cdc.Z(), vertex_cdc.X());
    FillHist("VertexCDC_ZY", vertex_cdc.Z(), vertex_cdc.Y());
    if(-0.5<vertex.Z()&&vertex.Z()<0.5){
      FillHist("VertexCDC_XY_ZCut", vertex_cdc.X(), vertex_cdc.Y());
    }
    if(-0.5<vertex.Y()&&vertex.Y()<0.5){
      FillHist("VertexCDC_ZX_XCut", vertex_cdc.Z(), vertex_cdc.X());
    }
    if(-0.5<vertex.X()&&vertex.X()<0.5){
      FillHist("VertexCDC_ZY_XCut", vertex_cdc.Z(), vertex_cdc.Y());
    }

    FillHist("VertexBPC_XY", vertex_bpc.X(), vertex_bpc.Y());
    FillHist("VertexBPC_ZX", vertex_bpc.Z(), vertex_bpc.X());
    FillHist("VertexBPC_ZY", vertex_bpc.Z(), vertex_bpc.Y());
    if(-0.5<vertex.Z()&&vertex.Z()<0.5){
      FillHist("VertexBPC_XY_ZCut", vertex_bpc.X(), vertex_bpc.Y());
    }
    if(-0.5<vertex.Y()&&vertex.Y()<0.5){
      FillHist("VertexBPC_ZX_XCut", vertex_bpc.Z(), vertex_bpc.X());
    }
    if(-0.5<vertex.X()&&vertex.X()<0.5){
      FillHist("VertexBPC_ZY_XCut", vertex_bpc.Z(), vertex_bpc.Y());
    }

    FillHist("VertexDiffBPC_XY", vertex_diff.X(), vertex_diff.Y());
    FillHist("VertexDiffBPC_ZX", vertex_diff.Z(), vertex_diff.X());
    FillHist("VertexDiffBPC_ZY", vertex_diff.Z(), vertex_diff.Y());
    if(-0.5<vertex.Z()&&vertex.Z()<0.5){
      FillHist("VertexDiffBPC_XY_ZCut", vertex_diff.X(), vertex_diff.Y());
    }
    if(-0.5<vertex.Y()&&vertex.Y()<0.5){
      FillHist("VertexDiffBPC_ZX_XCut", vertex_diff.Z(), vertex_diff.X());
    }
    if(-0.5<vertex.X()&&vertex.X()<0.5){
      FillHist("VertexDiffBPC_ZY_XCut", vertex_diff.Z(), vertex_diff.Y());
    }

    FillHist("VertexCDC_XvsVertexDiffBPC_X", vertex_cdc.X(), vertex_diff.X());
    FillHist("VertexCDC_XvsVertexDiffBPC_Y", vertex_cdc.X(), vertex_diff.Y());
    FillHist("VertexCDC_YvsVertexDiffBPC_X", vertex_cdc.Y(), vertex_diff.X());
    FillHist("VertexCDC_YvsVertexDiffBPC_Y", vertex_cdc.Y(), vertex_diff.Y());

  }

  for(int it=0; it<particle->nProduct(); it++){
    pCDS* cds = particle->product(it);
    TVector3 vtxb = cds->vbeam();
    TVector3 vertex_cdc = cds->vertex();
    double bpcz = vertex_cdc.Z();
    double bpcx=0,bpcy=0;
    beam->blcpos(bpcz,bpcx,bpcy);
    TVector3 vertex_bpc(bpcx,bpcy,bpcz);
    double diffx = vertex_bpc.X()-vertex_cdc.X();
    double diffy = vertex_bpc.Y()-vertex_cdc.Y();
    TVector3 vertex_diff(diffx,diffy,bpcz);
  
    FillHist("VertexCDC_XY", vertex_cdc.X(), vertex_cdc.Y());
    FillHist("VertexCDC_ZX", vertex_cdc.Z(), vertex_cdc.X());
    FillHist("VertexCDC_ZY", vertex_cdc.Z(), vertex_cdc.Y());
    if(-0.5<vertex.Z()&&vertex.Z()<0.5){
      FillHist("VertexCDC_XY_ZCut", vertex_cdc.X(), vertex_cdc.Y());
    }
    if(-0.5<vertex.Y()&&vertex.Y()<0.5){
      FillHist("VertexCDC_ZX_XCut", vertex_cdc.Z(), vertex_cdc.X());
    }
    if(-0.5<vertex.X()&&vertex.X()<0.5){
      FillHist("VertexCDC_ZY_XCut", vertex_cdc.Z(), vertex_cdc.Y());
    }

    FillHist("VertexBLC_XY", vertex_bpc.X(), vertex_bpc.Y());
    FillHist("VertexBLC_ZX", vertex_bpc.Z(), vertex_bpc.X());
    FillHist("VertexBLC_ZY", vertex_bpc.Z(), vertex_bpc.Y());
    if(-0.5<vertex.Z()&&vertex.Z()<0.5){
      FillHist("VertexBLC_XY_ZCut", vertex_bpc.X(), vertex_bpc.Y());
    }
    if(-0.5<vertex.Y()&&vertex.Y()<0.5){
      FillHist("VertexBLC_ZX_XCut", vertex_bpc.Z(), vertex_bpc.X());
    }
    if(-0.5<vertex.X()&&vertex.X()<0.5){
      FillHist("VertexBLC_ZY_XCut", vertex_bpc.Z(), vertex_bpc.Y());
    }

    FillHist("VertexDiffBLC_XY", vertex_diff.X(), vertex_diff.Y());
    FillHist("VertexDiffBLC_ZX", vertex_diff.Z(), vertex_diff.X());
    FillHist("VertexDiffBLC_ZY", vertex_diff.Z(), vertex_diff.Y());
    if(-0.5<vertex.Z()&&vertex.Z()<0.5){
      FillHist("VertexDiffBLC_XY_ZCut", vertex_diff.X(), vertex_diff.Y());
    }
    if(-0.5<vertex.Y()&&vertex.Y()<0.5){
      FillHist("VertexDiffBLC_ZX_XCut", vertex_diff.Z(), vertex_diff.X());
    }
    if(-0.5<vertex.X()&&vertex.X()<0.5){
      FillHist("VertexDiffBLC_ZY_XCut", vertex_diff.Z(), vertex_diff.Y());
    }

    FillHist("VertexCDC_XvsVertexDiffBLC_X", vertex_cdc.X(), vertex_diff.X());
    FillHist("VertexCDC_XvsVertexDiffBLC_Y", vertex_cdc.X(), vertex_diff.Y());
    FillHist("VertexCDC_YvsVertexDiffBLC_X", vertex_cdc.Y(), vertex_diff.X());
    FillHist("VertexCDC_YvsVertexDiffBLC_Y", vertex_cdc.Y(), vertex_diff.Y());

  }

  delete beam;
  delete particle;
  return true;
}

bool MyAnalysisCDCBPCCheck::FillHist(TString name, double val1, int weight)
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

bool MyAnalysisCDCBPCCheck::FillHist(TString name, TString val1, int weight)
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

bool MyAnalysisCDCBPCCheck::FillHist(TString name, double val1, double val2, int weight)
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

bool MyAnalysisCDCBPCCheck::FillHist(TString name, TString val1, TString val2, int weight)
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

bool MyAnalysisCDCBPCCheck::Initialize(ConfMan* confMan)
{
  std::cout << "### MyAnalysisCDCBPCCheck::Initialize ###" << std::endl;
  DetectorList *dlist=DetectorList::GetInstance();
  confFile = confMan;
  std::string ofname = confMan->GetOutFileName();
  ofname.insert(ofname.find(".root"),"_anaCDCBPCCheck");
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
  new TH2F( "CDS_Mass2vsMomentum_Fiducial", "Mass^{2} vs. Momentum;Mass^{2} (GeV/c^{2})^{2};Momentum (GeV/c)", 3000, -5, 10, 800, -2, 2 );
  new TH2F( "CDS_OverbetavsMomentum", "1/#beta vs. Momentum;1/#beta;Momentum (GeV/c)", 1000, 0, 10, 800, -2, 2 );
  //new TH2F( "CDS_OverbetavsMomentum_AfterFindMass", "1/#beta vs. Momentum;1/#beta;Momentum (GeV/c)", 1000, 0, 10, 800, -2, 2 );
  //new TH2F( "CDS_OverbetavsMomentum_AfterELossCorr", "1/#beta vs. Momentum;1/#beta;Momentum (GeV/c)", 1000, 0, 10, 800, -2, 2 );
  new TH2F( "CDS_OverbetavsMomentum_Fiducial", "1/#beta vs. Momentum;1/#beta;Momentum (GeV/c)", 1000, 0, 10, 800, -2, 2 );

  // Vertex
  std::cout << "Define Histograms for Vertex" << std::endl;
  new TH2F( "Vertex_XY", "Vertex XY plane;x (cm);y (cm)", 600, -15, 15, 600, -15, 15 );
  new TH2F( "Vertex_ZX", "Vertex ZX plane;z (cm);x (cm)", 1200, -30, 30, 600, -15, 15 );
  new TH2F( "Vertex_ZY", "Vertex ZY plane;z (cm);y (cm)", 1200, -30, 30, 600, -15, 15 );
  new TH1F( "Vertex_Z", "Vertex Z axis;z (cm);Counts", 1200, -30, 30 );
  new TH2F( "Vertex_XY_ZCut", "Vertex XY plane;x (cm);y (cm)", 600, -15, 15, 600, -15, 15 );
  new TH2F( "Vertex_ZX_YCut", "Vertex ZX plane;z (cm);x (cm)", 1200, -30, 30, 600, -15, 15 );
  new TH2F( "Vertex_ZY_XCut", "Vertex ZY plane;z (cm);y (cm)", 1200, -30, 30, 600, -15, 15 );
  new TH1F( "Vertex_Z_XYCut", "Vertex Z axis;z (cm);Counts", 1200, -30, 30 );
  new TH2F( "Vertex_XY_Fiducial", "Vertex XY plane;x (cm);y (cm)", 600, -15, 15, 600, -15, 15 );
  new TH2F( "Vertex_ZX_Fiducial", "Vertex ZX plane;z (cm);x (cm)", 1200, -30, 30, 600, -15, 15 );
  new TH2F( "Vertex_ZY_Fiducial", "Vertex ZY plane;z (cm);y (cm)", 1200, -30, 30, 600, -15, 15 );
  new TH1F( "Vertex_Z_Fiducial", "Vertex Z axis;z (cm);Counts", 1200, -30, 30 );
  new TH2F( "Vertex_XY_FiducialZCut", "Vertex XY plane;x (cm);y (cm)", 600, -15, 15, 600, -15, 15 );
  new TH2F( "Vertex_ZX_FiducialYCut", "Vertex ZX plane;z (cm);x (cm)", 1200, -30, 30, 600, -15, 15 );
  new TH2F( "Vertex_ZY_FiducialXCut", "Vertex ZY plane;z (cm);y (cm)", 1200, -30, 30, 600, -15, 15 );
  new TH1F( "Vertex_Z_FiducialXYCut", "Vertex Z axis;z (cm);Counts", 1200, -30, 30 );
  new TH2F( "Vertex_XY_TargetCell", "Vertex XY plane;x (cm);y (cm)", 600, -15, 15, 600, -15, 15 );
  new TH2F( "Vertex_ZX_TargetCell", "Vertex ZX plane;z (cm);x (cm)", 1200, -30, 30, 600, -15, 15 );
  new TH2F( "Vertex_ZY_TargetCell", "Vertex ZY plane;z (cm);y (cm)", 1200, -30, 30, 600, -15, 15 );
  new TH2F( "Vertex_XY_CellRing", "Vertex XY plane;x (cm);y (cm)", 600, -15, 15, 600, -15, 15 );
  new TH2F( "Vertex_ZX_CellRing", "Vertex ZX plane;z (cm);x (cm)", 1200, -30, 30, 600, -15, 15 );
  new TH2F( "Vertex_ZY_CellRing", "Vertex ZY plane;z (cm);y (cm)", 1200, -30, 30, 600, -15, 15 );
  new TH2F( "Vertex_XY_CellRingZCut", "Vertex XY plane;x (cm);y (cm)", 600, -15, 15, 600, -15, 15 );
  new TH2F( "Vertex_ZX_CellRingYCut", "Vertex ZX plane;z (cm);x (cm)", 1200, -30, 30, 600, -15, 15 );
  new TH2F( "Vertex_ZY_CellRingXCut", "Vertex ZY plane;z (cm);y (cm)", 1200, -30, 30, 600, -15, 15 );
  new TH2F( "Vertex_XY_CellTube", "Vertex XY plane;x (cm);y (cm)", 600, -15, 15, 600, -15, 15 );
  new TH2F( "Vertex_ZX_CellTube", "Vertex ZX plane;z (cm);x (cm)", 1200, -30, 30, 600, -15, 15 );
  new TH2F( "Vertex_ZY_CellTube", "Vertex ZY plane;z (cm);y (cm)", 1200, -30, 30, 600, -15, 15 );
  new TH2F( "Vertex_XY_CellTubeZCut", "Vertex XY plane;x (cm);y (cm)", 600, -15, 15, 600, -15, 15 );
  new TH2F( "Vertex_ZX_CellTubeYCut", "Vertex ZX plane;z (cm);x (cm)", 1200, -30, 30, 600, -15, 15 );
  new TH2F( "Vertex_ZY_CellTubeXCut", "Vertex ZY plane;z (cm);y (cm)", 1200, -30, 30, 600, -15, 15 );
  new TH2F( "Vertex_XY_CellAlBe", "Vertex XY plane;x (cm);y (cm)", 600, -15, 15, 600, -15, 15 );
  new TH2F( "Vertex_ZX_CellAlBe", "Vertex ZX plane;z (cm);x (cm)", 1200, -30, 30, 600, -15, 15 );
  new TH2F( "Vertex_ZY_CellAlBe", "Vertex ZY plane;z (cm);y (cm)", 1200, -30, 30, 600, -15, 15 );
  new TH2F( "Vertex_XY_CellAlBeZCut", "Vertex XY plane;x (cm);y (cm)", 600, -15, 15, 600, -15, 15 );
  new TH2F( "Vertex_ZX_CellAlBeYCut", "Vertex ZX plane;z (cm);x (cm)", 1200, -30, 30, 600, -15, 15 );
  new TH2F( "Vertex_ZY_CellAlBeXCut", "Vertex ZY plane;z (cm);y (cm)", 1200, -30, 30, 600, -15, 15 );

  new TH2F( "VertexCDC_XY", "Vertex XY plane;x (cm);y (cm)", 600, -15, 15, 600, -15, 15 );
  new TH2F( "VertexCDC_ZX", "Vertex ZX plane;z (cm);x (cm)", 1200, -30, 30, 600, -15, 15 );
  new TH2F( "VertexCDC_ZY", "Vertex ZY plane;z (cm);y (cm)", 1200, -30, 30, 600, -15, 15 );
  new TH2F( "VertexCDC_XY_ZCut", "Vertex XY plane;x (cm);y (cm)", 600, -15, 15, 600, -15, 15 );
  new TH2F( "VertexCDC_ZX_YCut", "Vertex ZX plane;z (cm);x (cm)", 1200, -30, 30, 600, -15, 15 );
  new TH2F( "VertexCDC_ZY_XCut", "Vertex ZY plane;z (cm);y (cm)", 1200, -30, 30, 600, -15, 15 );

  new TH2F( "VertexBPC_XY", "Vertex XY plane;x (cm);y (cm)", 600, -15, 15, 600, -15, 15 );
  new TH2F( "VertexBPC_ZX", "Vertex ZX plane;z (cm);x (cm)", 1200, -30, 30, 600, -15, 15 );
  new TH2F( "VertexBPC_ZY", "Vertex ZY plane;z (cm);y (cm)", 1200, -30, 30, 600, -15, 15 );
  new TH2F( "VertexBPC_XY_ZCut", "Vertex XY plane;x (cm);y (cm)", 600, -15, 15, 600, -15, 15 );
  new TH2F( "VertexBPC_ZX_YCut", "Vertex ZX plane;z (cm);x (cm)", 1200, -30, 30, 600, -15, 15 );
  new TH2F( "VertexBPC_ZY_XCut", "Vertex ZY plane;z (cm);y (cm)", 1200, -30, 30, 600, -15, 15 );

  new TH2F( "VertexBLC_XY", "Vertex XY plane;x (cm);y (cm)", 600, -15, 15, 600, -15, 15 );
  new TH2F( "VertexBLC_ZX", "Vertex ZX plane;z (cm);x (cm)", 1200, -30, 30, 600, -15, 15 );
  new TH2F( "VertexBLC_ZY", "Vertex ZY plane;z (cm);y (cm)", 1200, -30, 30, 600, -15, 15 );
  new TH2F( "VertexBLC_XY_ZCut", "Vertex XY plane;x (cm);y (cm)", 600, -15, 15, 600, -15, 15 );
  new TH2F( "VertexBLC_ZX_YCut", "Vertex ZX plane;z (cm);x (cm)", 1200, -30, 30, 600, -15, 15 );
  new TH2F( "VertexBLC_ZY_XCut", "Vertex ZY plane;z (cm);y (cm)", 1200, -30, 30, 600, -15, 15 );

  new TH2F( "VertexDiffBPC_XY", "Vertex XY plane;x (cm);y (cm)", 600, -1, 1, 600, -1, 1 );
  new TH2F( "VertexDiffBPC_ZX", "Vertex ZX plane;z (cm);x (cm)", 1200, -30, 30, 600, -1, 1 );
  new TH2F( "VertexDiffBPC_ZY", "Vertex ZY plane;z (cm);y (cm)", 1200, -30, 30, 600, -1, 1 );
  new TH2F( "VertexDiffBPC_XY_ZCut", "Vertex XY plane;x (cm);y (cm)", 600, -1, 1, 600, -1, 1 );
  new TH2F( "VertexDiffBPC_ZX_YCut", "Vertex ZX plane;z (cm);x (cm)", 1200, -30, 30, 600, -1, 1 );
  new TH2F( "VertexDiffBPC_ZY_XCut", "Vertex ZY plane;z (cm);y (cm)", 1200, -30, 30, 600, -1, 1 );

  new TH2F( "VertexCDC_XvsVertexDiffBPC_X", "VertexCDC X vs. VertexDiffBPC X plane;x (cm);y (cm)", 600, -15, 15, 600, -1, 1 );
  new TH2F( "VertexCDC_XvsVertexDiffBPC_Y", "VertexCDC X vs. VertexDiffBPC Y plane;x (cm);y (cm)", 600, -15, 15, 600, -1, 1 );
  new TH2F( "VertexCDC_YvsVertexDiffBPC_X", "VertexCDC Y vs. VertexDiffBPC X plane;x (cm);y (cm)", 600, -15, 15, 600, -1, 1 );
  new TH2F( "VertexCDC_YvsVertexDiffBPC_Y", "VertexCDC Y vs. VertexDiffBPC Y plane;x (cm);y (cm)", 600, -15, 15, 600, -1, 1 );

  new TH2F( "VertexDiffBLC_XY", "Vertex XY plane;x (cm);y (cm)", 600, -1, 1, 600, -1, 1 );
  new TH2F( "VertexDiffBLC_ZX", "Vertex ZX plane;z (cm);x (cm)", 1200, -30, 30, 600, -1, 1 );
  new TH2F( "VertexDiffBLC_ZY", "Vertex ZY plane;z (cm);y (cm)", 1200, -30, 30, 600, -1, 1 );
  new TH2F( "VertexDiffBLC_XY_ZCut", "Vertex XY plane;x (cm);y (cm)", 600, -1, 1, 600, -1, 1 );
  new TH2F( "VertexDiffBLC_ZX_YCut", "Vertex ZX plane;z (cm);x (cm)", 1200, -30, 30, 600, -1, 1 );
  new TH2F( "VertexDiffBLC_ZY_XCut", "Vertex ZY plane;z (cm);y (cm)", 1200, -30, 30, 600, -1, 1 );

  new TH2F( "VertexCDC_XvsVertexDiffBLC_X", "VertexCDC X vs. VertexDiffBLC X plane;x (cm);y (cm)", 600, -15, 15, 600, -1, 1 );
  new TH2F( "VertexCDC_XvsVertexDiffBLC_Y", "VertexCDC X vs. VertexDiffBLC Y plane;x (cm);y (cm)", 600, -15, 15, 600, -1, 1 );
  new TH2F( "VertexCDC_YvsVertexDiffBLC_X", "VertexCDC Y vs. VertexDiffBLC X plane;x (cm);y (cm)", 600, -15, 15, 600, -1, 1 );
  new TH2F( "VertexCDC_YvsVertexDiffBLC_Y", "VertexCDC Y vs. VertexDiffBLC Y plane;x (cm);y (cm)", 600, -15, 15, 600, -1, 1 );

  new TH2F( "PVertex_XY", "Vertex XY plane;x (cm);y (cm)", 600, -15, 15, 600, -15, 15 );
  new TH2F( "PVertex_ZX", "Vertex ZX plane;z (cm);x (cm)", 1200, -30, 30, 600, -15, 15 );
  new TH2F( "PVertex_ZY", "Vertex ZY plane;z (cm);y (cm)", 1200, -30, 30, 600, -15, 15 );
  new TH1F( "PVertex_Z", "Vertex Z axis;z (cm);Counts", 1200, -30, 30 );
  new TH2F( "PVertex_XY_ZCut", "Vertex XY plane;x (cm);y (cm)", 600, -15, 15, 600, -15, 15 );
  new TH1F( "PVertex_Z_XYCut", "Vertex Z axis;z (cm);Counts", 1200, -30, 30 );
  new TH2F( "PVertex_XY_Fiducial", "Vertex XY plane;x (cm);y (cm)", 600, -15, 15, 600, -15, 15 );
  new TH2F( "PVertex_ZX_Fiducial", "Vertex ZX plane;z (cm);x (cm)", 1200, -30, 30, 600, -15, 15 );
  new TH2F( "PVertex_ZY_Fiducial", "Vertex ZY plane;z (cm);y (cm)", 1200, -30, 30, 600, -15, 15 );
  new TH1F( "PVertex_Z_Fiducial", "Vertex Z axis;z (cm);Counts", 1200, -30, 30 );
  new TH2F( "PVertex_XY_FiducialZCut", "Vertex XY plane;x (cm);y (cm)", 600, -15, 15, 600, -15, 15 );
  new TH1F( "PVertex_Z_FiducialXYCut", "Vertex Z axis;z (cm);Counts", 1200, -30, 30 );

  std::cout << "=== End of [MyAnalysisCDCBPCCheck::Initialize] === " << std::endl;

  return true;
}
