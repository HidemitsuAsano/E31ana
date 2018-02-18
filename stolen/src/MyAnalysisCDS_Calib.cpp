// MyAnalysisCDS.cpp

#include "MyAnalysisCDS.h"

MyAnalysisCDS::MyAnalysisCDS(TFile* rtFile, ConfMan* conf)
{
  Initialize(rtFile, conf);
  Clear();
}

MyAnalysisCDS::~MyAnalysisCDS()
{
}

bool MyAnalysisCDS::Clear()
{
  BeamPID = Beam_Other;
  BeamMom = -9999.0;
  T0Segment = -1;
  T0Timing = -9999.0;
  trackblc1 = 0;
  trackblc2 = 0;
  trackbpc  = 0;
}

bool MyAnalysisCDS::DoAnalysis(ConfMan* conf, EventHeader* header, BeamLineHitMan* blMan, BeamLineTrackMan* bltrackMan, CDSHitMan* cdsMan, CDSTrackingMan* cdstrackMan)
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
    if( GeomTools::GetID(cvtxb2)==CID_Fiducial ){
      FillHist("VertexwTRG_Z",cvtxb2.Z());
      FillHist("VertexwTRG_ZX",cvtxb2.Z(),cvtxb2.X());
      FillHist("VertexwTRG_ZY",cvtxb2.Z(),cvtxb2.Y());
      FillHist("VertexwTRG_XY",cvtxb2.X(),cvtxb2.Y());
    }

    // MM calculation
    if( GeomTools::GetID(cvtxb2)==CID_Fiducial ){
      TVector3 cPp;
      if( track->GetMomentum(cvtxb1,cPp,true,true) ){
        TLorentzVector cL; cL.SetVectM( cPp, cdsMass[pid] );
        int combid=pow(2,pid);

        TVector3 ZVector;
        TLorentzVector dL; dL.SetVectM( ZVector, dMass );
        TLorentzVector pL; pL.SetVectM( ZVector, pMass );
        double dmm = ( dL+BeamL - cL).M();
        double pmm = ( pL+BeamL - cL).M();
        if(combid==pow(2,CDS_PiPlus)){
          FillHist(Form("MM_dpiplus"), dmm );
          FillHist(Form("MM_ppiplus"), pmm );
        }
        if(combid==pow(2,CDS_PiMinus)){
          FillHist(Form("MM_dpiminus"), dmm );
          FillHist(Form("MM_ppiminus"), pmm );
        }
        if(combid==pow(2,CDS_Kaon)){
          FillHist(Form("MM_dK"), dmm );
          FillHist(Form("MM_pK"), pmm );
        }
        if(combid==pow(2,CDS_Proton)){
          FillHist(Form("MM_dp"), dmm );
          FillHist(Form("MM_pp"), pmm );
        }
        if(combid==pow(2,CDS_Deuteron)){
          FillHist(Form("MM_dd"), dmm );
          FillHist(Form("MM_pd"), pmm );
        }
      }
    }

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
    if(cdstrackMan->nGoodTrack()!=2) continue;
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
      if( GeomTools::GetID(cvtx)!=CID_Fiducial ) continue;

      // IM calculation
      double MK0[2] = {0.470,0.525};
      double ML[2] = {1.109,1.121};
      double Mn[2] = {0.850,1.050};
      bool K0flag = false;
      bool Lflag = false;
      bool MissNflag = false;
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
        if(MK0[0]<=cim&&cim<=MK0[1]) K0flag = true;
      }
      if(combid==pow(2,CDS_Proton)+pow(2,CDS_PiMinus)){
        FillHist(Form("IM_pip"), cim );
        if(ML[0]<=cim&&cim<=ML[1]) Lflag = true;
      }
      // MM calculation
      TVector3 ZVector;
      TLorentzVector dL; dL.SetVectM( ZVector, dMass );
      TLorentzVector pL; pL.SetVectM( ZVector, pMass );
      double dmm = ( dL+BeamL - (cL1+cL2)).M();
      double pmm = ( pL+BeamL - (cL1+cL2)).M();
      if(combid==pow(2,CDS_PiPlus)+pow(2,CDS_PiMinus)){
        FillHist(Form("MM_dpipi"), dmm );
        FillHist(Form("MM_ppipi"), pmm );
        if(K0flag){
          FillHist(Form("MM_dK0"), dmm );
          FillHist(Form("MM_pK0"), pmm );
          if(Mn[0]<=pmm&&pmm<=Mn[1]) MissNflag = true;
        }
      }
      if(combid==pow(2,CDS_Proton)+pow(2,CDS_PiMinus)){
        FillHist(Form("MM_dpip"), dmm );
        FillHist(Form("MM_ppip"), pmm );
        if(Lflag){
          FillHist(Form("MM_dL"), dmm );
          FillHist(Form("MM_pL"), pmm );
        }
      }

      // K0
      if(K0flag){
        TVector3 np = (cL1+cL2).Vect();
        double cost = np.CosTheta();
        double dxdz = np.X()/np.Z();
        double dydz = np.Y()/np.Z();
        FillHist("K0_cosT",cost);
        FillHist("K0_dxdz",dxdz);
        FillHist("K0_dydz",dydz);
      }
      // Lambda
      if(Lflag){
        TVector3 np = (cL1+cL2).Vect();
        double cost = np.CosTheta();
        double dxdz = np.X()/np.Z();
        double dydz = np.Y()/np.Z();
        FillHist("Lambda_cosT",cost);
        FillHist("Lambda_dxdz",dxdz);
        FillHist("Lambda_dydz",dydz);
      }
      // Missing nuetron
      if(MissNflag){
        TVector3 np = (pL+BeamL - (cL1+cL2)).Vect();
        double cost = np.CosTheta();
        double dxdz = np.X()/np.Z();
        double dydz = np.Y()/np.Z();
        FillHist("MissN_cosT",cost);
        FillHist("MissN_dxdz",dxdz);
        FillHist("MissN_dydz",dydz);
      }

    }
  }	

  return true;

}

bool MyAnalysisCDS::FillHist(TString name, double val1)
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

bool MyAnalysisCDS::FillHist(TString name, double val1, double val2)
{
  TH2F* h2 = (TH2F*)gFile -> Get(name);
  if(h2){ h2 -> Fill(val1,val2);
    return true;
  }
  else {
    return false;
  }
}

bool MyAnalysisCDS::Initialize(TFile* rtFile, ConfMan* confMan)
{
  rtFile -> cd();

  // IH and CDH
  const int nhodo = 2;
  int NumOfSegments[nhodo] = {24,36};
  TString CounterName[nhodo] = {"IH","CDH"};
  for(int ihodo=0; ihodo<nhodo; ihodo++){
    std::cout << "Define Histograms for " << CounterName[ihodo].Data() << std::endl;
    new TH1F( Form("%s_HitPat",CounterName[ihodo].Data()), Form("Hit Pattern %s;%s segment",CounterName[ihodo].Data(),CounterName[ihodo].Data()), NumOfSegments[ihodo]+1, 0, NumOfSegments[ihodo]+1 );
    new TH2F( Form("%s_SegmentvPosition",CounterName[ihodo].Data()), Form("Hit Position %s;%s segment;Position (cm)",CounterName[ihodo].Data(),CounterName[ihodo].Data()), NumOfSegments[ihodo]+1, 0, NumOfSegments[ihodo]+1, 201, -50.25, 50.25 );
    new TH1F( Form("%s_Mult",CounterName[ihodo].Data()), Form("Multiplicity %s;Multiplicity;Counts",CounterName[ihodo].Data()), NumOfSegments[ihodo]+1, 0, NumOfSegments[ihodo]+1 );
  }
  for(int seg=1; seg<=NumOfSegments[1]; seg++){
    char* ud[2] = {"u","d"};
    for(int iud=0; iud<2; iud++){
      new TH1F( Form("CDH%s%d_ADC",ud[iud],seg),   Form("ADC CDHU%s%d;ADC ch.;Counts",ud[iud],seg),    4000,    0, 4000 );
      new TH1F( Form("CDH%s%d_dE",ud[iud],seg),   Form("dE CDHU%s%d;dE (MeV);Counts",ud[iud],seg),    400,    0, 40 );
      new TH1F( Form("CDH%s%d_Time",ud[iud],seg),   Form("Time CDH%s%d;Time (ns);Counts",ud[iud],seg),    4000,    -100, 100 );
      new TH1F( Form("CDH%s%d_CTime",ud[iud],seg),   Form("CTime CDHU%s%d;Time (ns);Counts",ud[iud],seg),    4000,    -100, 100 );
      new TH2F( Form("CDH%s%d_dEvCTime",ud[iud],seg),   Form("dE CTime corr. CDH%s%d;dE (MeV);Time (ns)",ud[iud],seg),     400,    0, 40,  4000,    -100, 100 );
      new TH2F( Form("CDH%s%d_dEvCTimeMean",ud[iud],seg),   Form("dE CTimeMean corr. CDH%s%d;dE (MeV);Time (ns)",ud[iud],seg),     400,    0, 40,  4000,    -100, 100 );
      new TH2F( Form("CDH%s%d_dEvOffset",ud[iud],seg),   Form("dE CTimeMean corr. CDH%s%d;dE (MeV);Offset (ns)",ud[iud],seg),     400,    0, 40,  4000,    -100, 100 );
    }
    new TH1F( Form("CDH%d_dEMean",seg),   Form("Mean dE CDH%d;dE (MeV);Counts",seg),    400,    0, 40 );
    new TH1F( Form("CDH%d_TimeMean",seg),   Form("Mean Time CDH%d;Time (ns);Counts",seg),    4000,    -100, 100 );
    new TH1F( Form("CDH%d_CTimeMean",seg),   Form("Mean CTime CDH%d;Time (ns);Counts",seg),    4000,    -100, 100 );
    new TH2F( Form("CDH%d_dEMeanvCTimeMean",seg),   Form("dEMean CTimeMean corr. CDH%d;dE (MeV);Time (ns)",seg),     400,    0, 40,  4000,    -100, 100 );
    new TH2F( Form("CDH%d_dEMeanvOffset",seg),   Form("dEMean Offset corr. CDH%d;dE (MeV);Offset (ns)",seg),     400,    0, 40,  4000,    -100, 100 );
    new TH1F( Form("CDH%d_Offset",seg),   Form("Offset CDH%d;Offset (ns);Counts",seg),    4000,    -100, 100 );
  }
  // CDC
  const int ndc=1;
  int dccid[ndc]={CID_CDC};
  TString dcname[ndc]={"CDC"};
  int NumOfLayers[ndc]={NumOfCDCLayers};

  for(int idc=0;idc<ndc;idc++){
    std::cout << "Define Histgram for " << dcname[idc] << std::endl;
    new TH1F( Form("%s_nGoodTrack",dcname[idc].Data()), Form("N good track %s;# of good tracks;Counts",dcname[idc].Data()), 10, 0, 10 );
    new TH1F( Form("%s_Chi2",dcname[idc].Data()), Form("#chi^{2} %s;#chi^{2};Counts",dcname[idc].Data()), 101, 0.0, 101.0 );
  }
  // TOF
  std::cout << "Define Histograms for TOF" << std::endl;
  for(int ihodo=0; ihodo<nhodo; ihodo++){
    new TH2F( Form("T0%s_TOF",CounterName[ihodo].Data()), Form("TOF T0-%s;%s segment;TOF (ns)",CounterName[ihodo].Data(),CounterName[ihodo].Data()), NumOfSegments[ihodo], 1, NumOfSegments[ihodo]+1, 4000, -50, 50 );
  }
  // CDS Particles
  std::cout << "Define Histograms for CDS particles" << std::endl;
  new TH2F( "Mass2vsMomentum", "Mass^{2} vs. Momentum;Mass^{2} (GeV/c^{2})^{2};Momentum (GeV/c)", 1500, -5, 10, 500, -5, 5 );
  new TH2F( "OverbetavsMomentum", "1/#beta vs. Momentum;1/#beta;Momentum (GeV/c)", 1000, -0, 10, 500, -5, 5 );
  new TH2F( "TOFvsMomentum", "TOF vs. Momentum;TOF (ns);Momentum (GeV/c)", 1000, 0, 20, 500, -5, 5 );

  // Vertex
  std::cout << "Define Histograms for Vertex" << std::endl;
  new TH2F( "Vertex_XY", "Vertex XY plane;x (cm);y (cm)", 600, -15, 15, 600, -15, 15 );
  new TH2F( "Vertex_ZX", "Vertex ZX plane;z (cm);x (cm)", 1200, -30, 30, 600, -15, 15 );
  new TH2F( "Vertex_ZY", "Vertex ZY plane;z (cm);y (cm)", 1200, -30, 30, 600, -15, 15 );
  new TH1F( "Vertex_Z", "Vertex Z axis;z (cm);Counts", 1200, -30, 30 );
  new TH2F( "VertexwTRG_XY", "Vertex XY plane;x (cm);y (cm)", 600, -15, 15, 600, -15, 15 );
  new TH2F( "VertexwTRG_ZX", "Vertex ZX plane;z (cm);x (cm)", 1200, -30, 30, 600, -15, 15 );
  new TH2F( "VertexwTRG_ZY", "Vertex ZY plane;z (cm);y (cm)", 1200, -30, 30, 600, -15, 15 );
  new TH1F( "VertexwTRG_Z", "Vertex Z axis;z (cm);Counts", 1200, -30, 30 );
  // Invariant mass
  std::cout << "Define Histograms for Invariant mass" << std::endl;
  new TH1F( "IM_pipi", "Invariant mass #pi^{-}#pi^{+};IM(#pi^{-}#pi^{+}) (GeV/c^{2});Counts", 2000, 0.0, 2.0);
  new TH1F( "IM_piK", "Invariant mass #pi^{-}K;IM(#pi^{-}K^{+}) (GeV/c^{2});Counts", 2000, 0.0, 2.0);
  new TH1F( "IM_pip", "Invariant mass #pi^{-}p;IM(#pi^{-}p) (GeV/c^{2});Counts", 2000, 0.0, 2.0);
  new TH1F( "IM_Kpi", "Invariant mass K^{-}#pi^{+};IM(K^{-}#pi^{+}) (GeV/c^{2});Counts", 2000, 0.0, 2.0);
  new TH1F( "IM_KK", "Invariant mass K^{-}K^{+};IM(K^{-}K^{+}) (GeV/c^{2});Counts", 2000, 0.0, 2.0);
  new TH1F( "IM_Kp", "Invariant mass K^{-}p;IM(K^{-}p) (GeV/c^{2});Counts", 2000, 0.0, 2.0);
  // Missing mass
  std::cout << "Define Histograms for Missing mass" << std::endl;
  // 1track
  new TH1F( "MM_dpiplus", "Missing mass d(K^{-},#pi^{+})X;MM_{d}(#pi^{+}) (GeV/c^{2});Counts", 2000, 1.0, 3.0);
  new TH1F( "MM_dpiminus", "Missing mass d(K^{-},#pi^{-})X;MM_{d}(#pi^{-}) (GeV/c^{2});Counts", 2000, 1.0, 3.0);
  new TH1F( "MM_dK", "Missing mass d(K^{-},K^{-})X;MM_{d}(#K^{-}) (GeV/c^{2});Counts", 2000, 1.0, 3.0);
  new TH1F( "MM_dp", "Missing mass d(K^{-},p)X;MM_{d}(p) (GeV/c^{2});Counts", 2000, 0.0, 2.0);
  new TH1F( "MM_dd", "Missing mass d(K^{-},d)X;MM_{d}(d) (GeV/c^{2});Counts", 2000, 0.0, 2.0);
  new TH1F( "MM_ppiplus", "Missing mass p(K^{-},#pi^{+})X;MM_{p}(#pi^{+}) (GeV/c^{2});Counts", 2000, 0.0, 2.0);
  new TH1F( "MM_ppiminus", "Missing mass p(K^{-},#pi^{-})X;MM_{p}(#pi^{-}) (GeV/c^{2});Counts", 2000, 0.0, 2.0);
  new TH1F( "MM_pK", "Missing mass p(K^{-},K^{-})X;MM_{p}(#K^{-}) (GeV/c^{2});Counts", 2000, 0.0, 2.0);
  new TH1F( "MM_pp", "Missing mass p(K^{-},p)X;MM_{p}(p) (GeV/c^{2});Counts", 2000, 0.0, 2.0);
  new TH1F( "MM_pd", "Missing mass p(K^{-},d)X;MM_{p}(d) (GeV/c^{2});Counts", 2000, 0.0, 2.0);
  // 2track
  new TH1F( "MM_dpipi", "Missing mass d(K^{-},#pi^{-}#pi^{+})X;MM_{d}(#pi^{-}#pi^{+}) (GeV/c^{2});Counts", 2000, 1.0, 3.0);
  new TH1F( "MM_dpip", "Missing mass d(K^{-},#pi^{-}p)X;MM_{d}(#pi^{-}p) (GeV/c^{2});Counts", 2000, 0.0, 2.0);
  new TH1F( "MM_dK0", "Missing mass d(K^{-},K^{0})X;MM_{d}(K^{0}) (GeV/c^{2});Counts", 2000, 1.0, 3.0);
  new TH1F( "MM_dL", "Missing mass d(K^{-},#Lambda)X;MM_{d}(#Lambda) (GeV/c^{2});Counts", 2000, 0.0, 2.0);
  new TH1F( "MM_ppipi", "Missing mass p(K^{-},#pi^{-}#pi^{+})X;MM_{p}(#pi^{-}#pi^{+}) (GeV/c^{2});Counts", 2000, 0.0, 2.0);
  new TH1F( "MM_ppip", "Missing mass p(K^{-},#pi^{-}p)X;MM_{p}(#pi^{-}p) (GeV/c^{2});Counts", 2000, 0.0, 2.0);
  new TH1F( "MM_pK0", "Missing mass p(K^{-},K^{0})X;MM_{p}(K^{0}) (GeV/c^{2});Counts", 2000, 0.0, 2.0);
  new TH1F( "MM_pL", "Missing mass p(K^{-},#Lambda)X;MM_{p}(#Lambda) (GeV/c^{2});Counts", 2000, 0.0, 2.0);
  // Angular distribution
  new TH1F( "K0_cosT", "K^{0} cos#theta;cos#theta;Counts", 200, -1.0, 1.0);
  new TH1F( "K0_dxdz", "K^{0} dx/dz;dx/dz;Counts", 200, -1.0, 1.0);
  new TH1F( "K0_dydz", "K^{0} dy/dz;dx/dz;Counts", 200, -1.0, 1.0);
  new TH1F( "Lambda_cosT", "#Lambda cos#theta;cos#theta;Counts", 200, -1.0, 1.0);
  new TH1F( "Lambda_dxdz", "#Lambda dx/dz;dx/dz;Counts", 200, -1.0, 1.0);
  new TH1F( "Lambda_dydz", "#Lambda dy/dz;dx/dz;Counts", 200, -1.0, 1.0);
  new TH1F( "MissN_cosT", "Missing n cos#theta;cos#theta;Counts", 200, -1.0, 1.0);
  new TH1F( "MissN_dxdz", "Missing n dx/dz;dx/dz;Counts", 200, -1.0, 1.0);
  new TH1F( "MissN_dydz", "Missing n dy/dz;dx/dz;Counts", 200, -1.0, 1.0);
  
}

