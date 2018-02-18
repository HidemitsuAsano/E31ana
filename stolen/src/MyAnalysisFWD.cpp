// MyAnalysisFWD.cpp

#include "MyAnalysisFWD.h"

MyAnalysisFWD::MyAnalysisFWD(TFile* rtFile, ConfMan* conf)
{
  Initialize(rtFile, conf);
  Clear();
}

MyAnalysisFWD::~MyAnalysisFWD()
{
}

bool MyAnalysisFWD::Clear()
{
  BeamPID = Beam_Other;
  BeamMom = -9999.0;
  T0Segment = -1;
  T0Timing = -9999.0;
  trackblc1 = 0;
  trackblc2 = 0;
  trackbpc  = 0;
}

bool MyAnalysisFWD::DoAnalysis(ConfMan* conf, EventHeader* header, BeamLineHitMan* blMan, BeamLineTrackMan* bltrackMan, CDSHitMan* cdsMan, CDSTrackingMan* cdstrackMan)
{
  DetectorList *dlist=DetectorList::GetInstance();

  // T0 timing
  for(int i=0; i<blMan->nT0(); i++){
    HodoscopeLikeHit* hit = blMan->T0(i);
    int seg = hit->seg();
    T0Segment = seg;
    if(hit->CheckRange()){
      T0Timing = hit->ctmean();
    }
  }
  // ################ //
  // Hodoscope        //
  // ################ //
  const int nhodo=2;
  int hodoid[nhodo]={CID_IH,CID_CDH};
  int hodonHit[nhodo];
  int hodoseg[nhodo][100];
  double hodoemean[nhodo][100];
  double hodoctmean[nhodo][100];
  double hodohitpos[nhodo][100];
  double hodotof[nhodo][100];
  TVector3 hodopos[nhodo][100];
  for(int ihodo=0; ihodo<nhodo; ihodo++){
    const int cid = hodoid[ihodo];
    const char* name = dlist->GetName(cid).data();
    const int nsegs = dlist->GetNsegs(cid);
    int ihit = 0;
    for(int i=0; i<cdsMan->nHodo(cid); i++){
      HodoscopeLikeHit* hit = cdsMan->Hodoi(cid,i);
      if(hit->CheckRange()){
        hodoseg[ihodo][ihit] = hit->seg();
        hodoemean[ihodo][ihit] = hit->emean();
        hodoctmean[ihodo][ihit] = hit->ctmean();
        hodohitpos[ihodo][ihit] = hit->hitpos();
        FillHist(Form("%s_HitPat",name),hodoseg[ihodo][ihit]);
        FillHist(Form("%s_SegmentvsPosition",name),hodoseg[ihodo][ihit],hodohitpos[ihodo][ihit]);

        // === TOF calculation === //
        if(T0Timing>-900){
          hodotof[ihodo][i] = hodoctmean[ihodo][ihit]-T0Timing;
          if(ihodo==0){ hodotof[ihodo][i]=hodotof[ihodo][i]*-1.; }
          FillHist(Form("T0%s_TOF",name),hodoseg[ihodo][ihit],hodotof[ihodo][i]);
        }
        // === TOF calculation === //

        hodonHit[ihodo]++;
        ihit++;
      }
    }
    FillHist(Form("%s_Mult",name),hodonHit[ihodo]);
  }

  // ################ //
  // CDS track        //
  // ################ //
  if(trackblc1==0&&trackblc2==0&&trackbpc==0) return false;
  TVector3 gpos;
  conf -> GetGeomMapManager() -> GetGPos(CID_T0, T0Segment, gpos);
  TVector3 T0Position = trackblc2 -> GetPosatZ(gpos.Z()-0.5);
  FillHist("CDC_nGoodTrack",cdstrackMan->nGoodTrack());

  for( int it=0; it<cdstrackMan->nGoodTrack(); it++ ){
    CDSTrack *track = cdstrackMan->Track( cdstrackMan->GoodTrackID(it) );
    double param[5];
    track->GetParameters(param);
    double chi2 = track->Chi();
    double mom = track->Momentum();
    TVector3 vtxb1, vtxb2;
    double tmpdis;
    TrackTools::CalcLineHelixVertex(trackbpc,track,vtxb1,vtxb2,tmpdis);
    track->SetPID(-1);
    if(!track->CDHFlag()) continue;

    double beta=-1.;
    double tof=999.;
    double mass2=-999.;
    double beamtof = 0.0;
    HodoscopeLikeHit* CDH = 0;

    for(int icdh=0;icdh<track->nCDHHit();icdh++){
      HodoscopeLikeHit *cdhhit=track->CDHHit(cdsMan,icdh);
      double tmptof    = cdhhit->ctmean()-T0Timing;      
      if(tmptof<tof||tof==999.){
        tof=tmptof;
        CDH = cdhhit;
      }
    }
    TrackTools::CalcBetaMass(vtxb1,trackbpc,track,conf,BeamPID,tof,beta,mass2);
    //if(TrackTools::FindMass2C(track, trackbpc, tof, BeamMom, BeamPID, beta, mass2, beamtof)) continue;
    double mass=sqrt(fabs(mass2));
    FillHist("OverbetavsMomentum",1.0/beta,mom);
    FillHist("TOFvsMomentum",tof,mom);
    FillHist("Mass2vsMomentum",mass2,mom);
    FillHist("CDC_Chi2",chi2);
    int pid=TrackTools::PID(mom,mass2);      
    track->SetPID(pid);

    TVector3 cvtxb1,cvtxb2,cPd;
    double tofvtxcdc,tmpl;
    track->GetVertex(trackbpc->GetPosatZ(0),trackbpc->GetMomDir(),cvtxb2,cvtxb1);
    track->CalcVertexTimeLength(trackbpc->GetPosatZ(-110.5),
        trackbpc->GetMomDir(),
        cdsMass[pid],
        cvtxb2,cvtxb1,tofvtxcdc,tmpl,true);
    FillHist("Vertex_Z",cvtxb2.Z());
    FillHist("Vertex_ZX",cvtxb2.Z(),cvtxb2.X());
    FillHist("Vertex_ZY",cvtxb2.Z(),cvtxb2.Y());
    FillHist("Vertex_XY",cvtxb2.X(),cvtxb2.Y());

    // == Hodoscope caliblation == //
    if(pid == CDS_PiPlus || pid == CDS_PiMinus){
      HodoscopeLikeHit* hit = CDH;
      int seg = hit->seg();
      int au = hit->adcu(), ad = hit->adcd();
      int tu = hit->tdcu(), td = hit->tdcd();
      double hitpos = hit->hitpos();
      double timeu = hit->time(0), timed = hit->time(1);
      double ctimeu = hit->ctime(0), ctimed = hit->ctime(1);
      double eneu = hit->ene(0), ened = hit->ene(1);
      double emean = hit->emean();
      double tmean = hit->tmean();
      double ctmean = hit->ctmean();

      double cdc_dis = MathTools::CalcHelixArc(param, vtxb2, track->CDHVertex());
      double v = 100.*Const*fabs(mom)/piMass/sqrt(1.0+mom*mom/(piMass*piMass));
      double cdc_calc_tof = cdc_dis/v;
      double beam_tof, out_mom;
      ELossTools::CalcElossBeamTGeo(trackbpc->GetPosatZ(-110.5), vtxb1, BeamMom, parMass[BeamPID], out_mom, beam_tof);
      double offset = hit->ctmean()-T0Timing-cdc_calc_tof-beam_tof;

      FillHist( Form("CDHu%d_ADC",seg), au );
      FillHist( Form("CDHd%d_ADC",seg), ad );
      FillHist( Form("CDHu%d_Time",seg), timeu );
      FillHist( Form("CDHu%d_CTime",seg), ctimeu );
      FillHist( Form("CDHd%d_Time",seg), timed );
      FillHist( Form("CDHd%d_CTime",seg), ctimed );
      FillHist( Form("CDHu%d_dE",seg), eneu );
      FillHist( Form("CDHd%d_dE",seg), ened );
      FillHist( Form("CDH%d_dEMean",seg), emean );
      FillHist( Form("CDH%d_TimeMean",seg), tmean );
      FillHist( Form("CDH%d_CTimeMean",seg), ctmean );
      FillHist( Form("CDHu%d_dEvTime",seg), eneu, timeu );
      FillHist( Form("CDHd%d_dEvTime",seg), ened, timed );
      FillHist( Form("CDHu%d_dEvCTime",seg), eneu, ctimeu );
      FillHist( Form("CDHd%d_dEvCTime",seg), ened, ctimed );
      FillHist( Form("CDHu%d_dEvCTimeMean",seg), eneu, ctmean );
      FillHist( Form("CDHd%d_dEvCTimeMean",seg), ened, ctmean );
      FillHist( Form("CDHu%d_dEvOffset",seg), eneu, offset );
      FillHist( Form("CDHd%d_dEvOffset",seg), ened, offset );
      FillHist( Form("CDH%d_dEMeanvCTimeMean",seg), emean, ctmean );
      FillHist( Form("CDH%d_Offset",seg), offset );
      FillHist( Form("CDH%d_dEMeanvOffset",seg), emean, offset );
    }
    // == Hodoscope caliblation == //
  }

  // CDS 2 tracks
  for( int it1=0; it1<cdstrackMan->nGoodTrack(); it1++ ){
    for( int it2=it1+1; it2<cdstrackMan->nGoodTrack(); it2++ ){
      CDSTrack *track1 = cdstrackMan->Track( cdstrackMan->GoodTrackID(it1) );
      CDSTrack *track2 = cdstrackMan->Track( cdstrackMan->GoodTrackID(it2) );
      int pid1=track1->PID();
      int pid2=track2->PID();
      if(pid1<0||pid2<0) continue;
      double mass1=cdsMass[pid1];
      double mass2=cdsMass[pid2];

      // Vertex selection
      TVector3 cvtx1,cvtx2;
      if( !TrackTools::Calc2HelixVertex(track1,track2,cvtx1,cvtx2) )continue;
      TVector3 cvtx=(cvtx1+cvtx2)*0.5;
      //if( GeomTools::GetID(cvtx)!=CID_Fiducial ) continue;

      // IM calculation
      TVector3 cPp1,cPp2;
      if( !track1->GetMomentum(cvtx1,cPp1,true,true) ) continue;
      if( !track2->GetMomentum(cvtx2,cPp2,true,true) ) continue;
      TVector3 cPp= cPp1+cPp2;
      TLorentzVector cL1; cL1.SetVectM( cPp1, mass1 );
      TLorentzVector cL2; cL2.SetVectM( cPp2, mass2 );
      double cim = (cL1+cL2).M();	
      int combid=pow(2,pid1)+pow(2,pid2);
      if(combid==pow(2,CDS_PiPlus)+pow(2,CDS_PiMinus)){
        FillHist(Form("IM_pipi"), cim );
      }
      if(combid==pow(2,CDS_Proton)+pow(2,CDS_PiMinus)){
        FillHist(Form("IM_pip"), cim );
      }
    }
  }	

  return true;

}

bool MyAnalysisFWD::FillHist(TString name, double val1)
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

bool MyAnalysisFWD::FillHist(TString name, double val1, double val2)
{
  TH2F* h2 = (TH2F*)gFile -> Get(name);
  if(h2){ h2 -> Fill(val1,val2);
    return true;
  }
  else {
    return false;
  }
}

bool MyAnalysisFWD::Initialize(TFile* rtFile, ConfMan* confMan)
{
  rtFile -> cd();

  // FWD counter
  const int nhodo = 2;
  int NumOfSegments[nhodo] = {34,114};
  TString CounterName[nhodo] = {"CVC","NC"};
  for(int ihodo=0; ihodo<nhodo; ihodo++){
    std::cout << "Define Histograms for " << CounterName[ihodo].Data() << std::endl;
    new TH1F( Form("%s_HitPat",CounterName[ihodo].Data()), Form("Hit Pattern %s;%s segment",CounterName[ihodo].Data(),CounterName[ihodo].Data()), NumOfSegments[ihodo]+1, 0, NumOfSegments[ihodo]+1 );
    new TH2F( Form("%s_SegmentvPosition",CounterName[ihodo].Data()), Form("Hit Position %s;%s segment;Position (cm)",CounterName[ihodo].Data(),CounterName[ihodo].Data()), NumOfSegments[ihodo]+1, 0, NumOfSegments[ihodo]+1, 201, -50.25, 50.25 );
    new TH1F( Form("%s_Mult",CounterName[ihodo].Data()), Form("Multiplicity %s;Multiplicity;Counts",CounterName[ihodo].Data()), NumOfSegments[ihodo]+1, 0, NumOfSegments[ihodo]+1 );
  }
  // TOF
  std::cout << "Define Histograms for TOF" << std::endl;
  for(int ihodo=0; ihodo<nhodo; ihodo++){
    new TH2F( Form("T0%s_TOF",CounterName[ihodo].Data()), Form("TOF T0-%s;%s segment;TOF (ns)",CounterName[ihodo].Data(),CounterName[ihodo].Data()), NumOfSegments[ihodo], 1, NumOfSegments[ihodo]+1, 4000, -50, 50 );
  }
  // FWD Particles
  std::cout << "Define Histograms for CDS particles" << std::endl;
  new TH2F( "Mass2vsMomentum", "Mass^{2} vs. Momentum;Mass^{2} (GeV/c^{2})^{2};Momentum (GeV/c)", 1500, -5, 10, 500, -5, 5 );
  new TH2F( "OverbetavsdE", "1/#beta vs. dE;1/#beta;Momentum (GeV/c)", 1000, -0, 10, 500, 0, 100 );

  // Vertex
  std::cout << "Define Histograms for Vertex" << std::endl;
  new TH2F( "Vertex_XY", "Vertex XY plane;x (cm);y (cm)", 600, -15, 15, 600, -15, 15 );
  new TH2F( "Vertex_ZX", "Vertex ZX plane;z (cm);x (cm)", 1200, -30, 30, 600, -15, 15 );
  new TH2F( "Vertex_ZY", "Vertex ZY plane;z (cm);y (cm)", 1200, -30, 30, 600, -15, 15 );
  new TH1F( "Vertex_Z", "Vertex Z axis;z (cm);Counts", 1200, -30, 30 );
  // Invariant mass
  std::cout << "Define Histograms for Invariant mass" << std::endl;
  new TH1F( "IM_pipi", "Invariant mass #pi^{-}#pi^{+};IM(#pi^{-}#pi^{+}) (GeV/c^{2});Counts", 2000, 0.0, 2.0);
  new TH1F( "IM_piK", "Invariant mass #pi^{-}K;IM(#pi^{-}K^{+}) (GeV/c^{2});Counts", 2000, 0.0, 2.0);
  new TH1F( "IM_pip", "Invariant mass #pi^{-}p;IM(#pi^{-}p) (GeV/c^{2});Counts", 2000, 0.0, 2.0);
  new TH1F( "IM_Kpi", "Invariant mass K^{-}#pi^{+};IM(K^{-}#pi^{+}) (GeV/c^{2});Counts", 2000, 0.0, 2.0);
  new TH1F( "IM_KK", "Invariant mass K^{-}K^{+};IM(K^{-}K^{+}) (GeV/c^{2});Counts", 2000, 0.0, 2.0);
  new TH1F( "IM_Kp", "Invariant mass K^{-}p;IM(K^{-}p) (GeV/c^{2});Counts", 2000, 0.0, 2.0);
}


