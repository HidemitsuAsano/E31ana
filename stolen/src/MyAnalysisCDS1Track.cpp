// MyAnalysisCDS1Track.cpp

#include "MyAnalysisCDS1Track.h"

MyAnalysisCDS1Track::MyAnalysisCDS1Track(TFile* rtFile, ConfMan* conf)
{
  Initialize(rtFile, conf);
  Clear();
}

MyAnalysisCDS1Track::~MyAnalysisCDS1Track()
{
}

bool MyAnalysisCDS1Track::Clear()
{
  BeamPID = Beam_Other;
  BeamMom = -9999.0;
  T0Segment = -1;
  T0Timing = -9999.0;
  trackblc1 = 0;
  trackblc2 = 0;
  trackbpc  = 0;
}

bool MyAnalysisCDS1Track::DoAnalysis(ConfMan* conf, EventHeader* header, BeamLineHitMan* blMan, BeamLineTrackMan* bltrackMan, CDSHitMan* cdsMan, CDSTrackingMan* cdstrackMan)
{
  DetectorList *dlist=DetectorList::GetInstance();
  if(trackblc1==0&&trackblc2==0&&trackbpc==0) return false;

  // ######################## //
  // CDS 1 track  check       //
  // ######################## //
  if(cdstrackMan->nGoodTrack()!=1) return false;

  // ################ //
  // CDS track        //
  // ################ //
  TVector3 gpos;
  conf -> GetGeomMapManager() -> GetGPos(CID_T0, T0Segment, gpos);
  TVector3 T0Position = trackblc2 -> GetPosatZ(gpos.Z()-0.5);

  CDSTrack* trackcds = cdstrackMan->Track( cdstrackMan->GoodTrackID(0) );
  double param[5];
  trackcds->GetParameters(param);
  double chi2 = trackcds->Chi();
  double mom = trackcds->Momentum();
  TVector3 vtxb1, vtxb2;
  double tmpdis;
  TrackTools::CalcLineHelixVertex(trackbpc,trackcds,vtxb1,vtxb2,tmpdis);
  trackcds->SetPID(-1);
  if(!trackcds->CDHFlag()) continue;

  double beta=-1.;
  double tof=999.;
  double mass2=-999.;
  double beamtof = 0.0;
  HodoscopeLikeHit* CDH = 0;

  for(int icdh=0;icdh<trackcds->nCDHHit();icdh++){
    HodoscopeLikeHit *cdhhit=trackcds->CDHHit(cdsMan,icdh);
    double tmptof    = cdhhit->ctmean()-T0Timing;      
    if(tmptof<tof||tof==999.){
      tof=tmptof;
      CDH = cdhhit;
    }
  }
  TrackTools::CalcBetaMass(vtxb1,trackbpc,trackcds,conf,BeamPID,tof,beta,mass2);
  if(!TrackTools::FindMass2(trackcds, trackbpc, tof, BeamMom, BeamPID, beta, mass2, beamtof)) continue;
  double mass=sqrt(fabs(mass2));
  FillHist("OverbetavsMomentum",1.0/beta,mom);
  FillHist("TOFvsMomentum",tof,mom);
  FillHist("Mass2vsMomentum",mass2,mom);
  FillHist("CDC_Chi2",chi2);
  int pid=TrackTools::PID(mom,mass2);      
  trackcds->SetPID(pid);

  TVector3 cvtxb1,cvtxb2,cPd;
  double tofvtxcdc,tmpl;
  trackcds->GetVertex(trackbpc->GetPosatZ(0),trackbpc->GetMomDir(),cvtxb2,cvtxb1);
  trackcds->CalcVertexTimeLength(trackbpc->GetPosatZ(-110.5),
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

  return true;
}

bool MyAnalysisCDS1Track::FillHist(TString name, double val1)
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

bool MyAnalysisCDS1Track::FillHist(TString name, double val1, double val2)
{
  TH2F* h2 = (TH2F*)gFile -> Get(name);
  if(h2){ h2 -> Fill(val1,val2);
    return true;
  }
  else {
    return false;
  }
}

bool MyAnalysisCDS1Track::Initialize(TFile* rtFile, ConfMan* confMan)
{
  rtFile -> cd();

  //  // IH and CDH
  const int nhodo = 2;
  int NumOfSegments[nhodo] = {24,36};
  TString CounterName[nhodo] = {"IH","CDH"};
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
  //// TOF
  std::cout << "Define Histograms for TOF" << std::endl;
  new TH2F( Form("T0%s_TOF",CounterName[1].Data()), Form("TOF T0-%s;%s segment;TOF (ns)",CounterName[1].Data(),CounterName[1].Data()), NumOfSegments[1], 1, NumOfSegments[1]+1, 4000, -50, 50 );
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
  new TH2F( "IMvsMom_pipi", "Invariant mass vs. Momentum #pi^{-}#pi^{+};IM(#pi^{-}#pi^{+}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F( "IMvsMom_pip", "Invariant mass vs. Momentum #pi^{-}p;IM(#pi^{-}p) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F( "IMvsMomCM_pipi", "Invariant mass vs. Momentum #pi^{-}#pi^{+};IM(#pi^{-}#pi^{+}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F( "IMvsMomCM_pip", "Invariant mass vs. Momentum #pi^{-}p;IM(#pi^{-}p) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  // Missing mass
  std::cout << "Define Histograms for Missing mass" << std::endl;
  // 2track
  // Mass vs. Momentum
  new TH2F( "MMvsMom_ppipi", "Missing mass  vs. Momentum p(K^{-},#pi^{-}#pi^{+})X;MM_{p}(#pi^{-}#pi^{+}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F( "MMvsMom_ppip", "Missing mass  vs. Momentum p(K^{-},#pi^{-}p)X;MM_{p}(#pi^{-}p) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F( "MMvsMom_pK0", "Missing mass  vs. Momentum p(K^{-},K^{0})X;MM_{p}(K^{0}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F( "MMvsMom_pL", "Missing mass  vs. Momentum p(K^{-},#Lambda)X;MM_{p}(#Lambda) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F( "MMvsMom_pK0_LowSide", "Missing mass  vs. Momentum p(K^{-},K^{0})X;MM_{p}(K^{0}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F( "MMvsMom_pL_LowSide", "Missing mass  vs. Momentum p(K^{-},#Lambda)X;MM_{p}(#Lambda) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F( "MMvsMom_pK0_UpSide", "Missing mass  vs. Momentum p(K^{-},K^{0})X;MM_{p}(K^{0}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F( "MMvsMom_pL_UpSide", "Missing mass  vs. Momentum p(K^{-},#Lambda)X;MM_{p}(#Lambda) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F( "MMvsMomCM_ppipi", "Missing mass  vs. Momentum p(K^{-},#pi^{-}#pi^{+})X;MM_{p}(#pi^{-}#pi^{+}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F( "MMvsMomCM_ppip", "Missing mass  vs. Momentum p(K^{-},#pi^{-}p)X;MM_{p}(#pi^{-}p) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F( "MMvsMomCM_pK0", "Missing mass  vs. Momentum p(K^{-},K^{0})X;MM_{p}(K^{0}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F( "MMvsMomCM_pL", "Missing mass  vs. Momentum p(K^{-},#Lambda)X;MM_{p}(#Lambda) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F( "MMvsMomCM_pK0_LowSide", "Missing mass  vs. Momentum p(K^{-},K^{0})X;MM_{p}(K^{0}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F( "MMvsMomCM_pL_LowSide", "Missing mass  vs. Momentum p(K^{-},#Lambda)X;MM_{p}(#Lambda) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F( "MMvsMomCM_pK0_UpSide", "Missing mass  vs. Momentum p(K^{-},K^{0})X;MM_{p}(K^{0}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F( "MMvsMomCM_pL_UpSide", "Missing mass  vs. Momentum p(K^{-},#Lambda)X;MM_{p}(#Lambda) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
}
