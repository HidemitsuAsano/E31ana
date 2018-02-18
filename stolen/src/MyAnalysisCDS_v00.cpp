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
        //FillHist(Form("%s_HitPat",name),hodoseg[ihodo][ihit]);
        //FillHist(Form("%s_SegmentvsPosition",name),hodoseg[ihodo][ihit],hodohitpos[ihodo][ihit]);

        // === TOF calculation === //
        if(T0Timing>-900){
          hodotof[ihodo][i] = hodoctmean[ihodo][ihit]-T0Timing;
          if(ihodo==0){ hodotof[ihodo][i]=hodotof[ihodo][i]*-1.; }
          //FillHist(Form("T0%s_TOF",name),hodoseg[ihodo][ihit],hodotof[ihodo][i]);
        }
        // === TOF calculation === //

        hodonHit[ihodo]++;
        ihit++;
      }
    }
    //FillHist(Form("%s_Mult",name),hodonHit[ihodo]);
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
    if(!TrackTools::FindMass2C(track, trackbpc, tof, BeamMom, BeamPID, beta, mass2, beamtof)) continue;
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

    //    // == Hodoscope caliblation == //
    //    if(pid == CDS_PiPlus || pid == CDS_PiMinus){
    //      HodoscopeLikeHit* hit = CDH;
    //      int seg = hit->seg();
    //      double emean = hit->emean();
    //      double ctmean = hit->ctmean();
    //
    //      double cdc_dis = MathTools::CalcHelixArc(param, vtxb2, track->CDHVertex());
    //      double v = 100.*Const*fabs(mom)/piMass/sqrt(1.0+mom*mom/(piMass*piMass));
    //      double cdc_calc_tof = cdc_dis/v;
    //      double beam_tof, out_mom;
    //      ELossTools::CalcElossBeamTGeo(trackbpc->GetPosatZ(-110.5), vtxb1, BeamMom, parMass[BeamPID], out_mom, beam_tof);
    //      double offset = hit->ctmean()-T0Timing-cdc_calc_tof-beam_tof;
    //
    //      FillHist( Form("CDH%d_Offset",seg), offset );
    //      FillHist( Form("CDH%d_dEMeanvOffset",seg), emean, offset );
    //    }
    //    // == Hodoscope caliblation == //
  }

  // CDS 2 tracks
  if(cdstrackMan->nGoodTrack()!=2) return false;
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
      if( GeomTools::GetID(cvtx)!=CID_Fiducial ) continue;

      // IM calculation
      double MK0[2] = {0.470,0.525};
      double MLK0[2] = {0.415,0.470};
      double MUK0[2] = {0.525,0.580};
      double ML[2] = {1.110,1.120};
      double MLL[2] = {1.100,1.110};
      double MUL[2] = {1.120,1.130};
      double Mn[2] = {0.850,1.050};
      bool K0flag = false;
      bool K0Lflag = false;
      bool K0Uflag = false;
      bool Lflag = false;
      bool LLflag = false;
      bool LUflag = false;
      bool MissNflag = false;


      TVector3 cPp1,cPp2;
      if( !track1->GetMomentum(cvtx1,cPp1,true,true) ) continue;
      if( !track2->GetMomentum(cvtx2,cPp2,true,true) ) continue;
      TVector3 cPp= cPp1+cPp2;
      TLorentzVector cL1; cL1.SetVectM( cPp1, mass1 );
      TLorentzVector cL2; cL2.SetVectM( cPp2, mass2 );

      TVector3 ZVector;
      TLorentzVector dL; dL.SetVectM( ZVector, dMass );
      TLorentzVector pL; pL.SetVectM( ZVector, pMass );
      TVector3 pbv = (pL+BeamL).BoostVector();

      double cim = (cL1+cL2).M();	
      TVector3 cp = (cL1+cL2).Vect();
      double pmom = cp.Mag();
      double cost = cp.CosTheta();
      double dxdz = cp.X()/cp.Z();
      double dydz = cp.Y()/cp.Z();

      TLorentzVector pbL = cL1+cL2; pbL.Boost(-pbv);
      TVector3 pbp = pbL.Vect();
      double bpmom = pbp.Mag();
      double bcost = pbp.CosTheta();
      double bdxdz = pbp.X()/pbp.Z();
      double bdydz = pbp.Y()/pbp.Z();

      int combid=pow(2,pid1)+pow(2,pid2);
      if(combid==pow(2,CDS_PiPlus)+pow(2,CDS_PiMinus)){
        FillHist(Form("IMvsMom_pipi"), cim, pmom );
        FillHist(Form("IMvsMomCM_pipi"), cim, bpmom );
        if(MK0[0]<cim&&cim<MK0[1]) { K0flag = true; }
        if(MLK0[0]<cim&&cim<MLK0[1]) { K0Lflag = true; }
        if(MUK0[0]<cim&&cim<MUK0[1]) { K0Uflag = true; }
      }
      if(combid==pow(2,CDS_Proton)+pow(2,CDS_PiMinus)){
        FillHist(Form("IMvsMom_pip"), cim, pmom );
        FillHist(Form("IMvsMomCM_pip"), cim, bpmom );
        if(ML[0]<cim&&cim<ML[1]) { Lflag = true; }
        if(MLL[0]<cim&&cim<MLL[1]) { LLflag = true; }
        if(MUL[0]<cim&&cim<MUL[1]) { LUflag = true; }
      }
      // MM calculation
      double pmm = ( pL+BeamL - (cL1+cL2)).M();
      TVector3 mp = (pL+BeamL - (cL1+cL2)).Vect();
      double mpmom = mp.Mag();
      double mcost = mp.CosTheta();
      double mdxdz = mp.X()/mp.Z();
      double mdydz = mp.Y()/mp.Z();
      TLorentzVector mpbL = (pL+BeamL - (cL1+cL2)); mpbL.Boost(-pbv);
      TVector3 mpbp = mpbL.Vect();
      double mbpmom = mpbp.Mag();
      double mbcost = mpbp.CosTheta();
      double mbdxdz = mpbp.X()/mpbp.Z();
      double mbdydz = mpbp.Y()/mpbp.Z();
      if(combid==pow(2,CDS_PiPlus)+pow(2,CDS_PiMinus)){
        FillHist(Form("MMvsMom_ppipi"), pmm, mpmom );
        FillHist(Form("MMvsMomCM_ppipi"), pmm, mbpmom );
        if(K0flag){
          FillHist(Form("MMvsMom_pK0"), pmm, mpmom );
          FillHist(Form("MMvsMomCM_pK0"), pmm, mbpmom );
          if(Mn[0]<=pmm&&pmm<=Mn[1]) MissNflag = true;
        }
        if(K0Lflag){
          FillHist(Form("MMvsMom_pK0_LowSide"), pmm, mpmom );
          FillHist(Form("MMvsMomCM_pK0_LowSide"), pmm, mbpmom );
        }
        if(K0Uflag){
          FillHist(Form("MMvsMom_pK0_UpSide"), pmm, mpmom );
          FillHist(Form("MMvsMomCM_pK0_UpSide"), pmm, mbpmom );
        }
      }
      if(combid==pow(2,CDS_Proton)+pow(2,CDS_PiMinus)){
        FillHist(Form("MMvsMom_ppip"), pmm, mpmom );
        FillHist(Form("MMvsMomCM_ppip"), pmm, mbpmom );
        if(Lflag){
          FillHist(Form("MMvsMom_pL"), pmm, mpmom );
          FillHist(Form("MMvsMomCM_pL"), pmm, mbpmom );
        }
        if(LLflag){
          FillHist(Form("MMvsMom_pL_LowSide"), pmm, mpmom );
          FillHist(Form("MMvsMomCM_pL_LowSide"), pmm, mbpmom );
        }
        if(LUflag){
          FillHist(Form("MMvsMom_pL_UpSide"), pmm, mpmom );
          FillHist(Form("MMvsMomCM_pL_UpSide"), pmm, mbpmom );
        }
      }

      // K0
      if(K0flag){
        FillHist("K0_cosTvsMomentum",cost, pmom);
        FillHist("K0_dxdzvdydz",dxdz,dydz);
        FillHist("K0_cosTvsMomentum_CM",bcost,bpmom);
        FillHist("K0_dxdzvdydz_CM",bdxdz,bdydz);
      }
      if(K0Lflag){
        FillHist("K0_LowSide_cosTvsMomentum",cost, pmom);
        FillHist("K0_LowSide_dxdzvdydz",dxdz,dydz);
        FillHist("K0_LowSide_cosTvsMomentum_CM",bcost,bpmom);
        FillHist("K0_LowSide_dxdzvdydz_CM",bdxdz,bdydz);
      }
      if(K0Uflag){
        FillHist("K0_UpSide_cosTvsMomentum",cost, pmom);
        FillHist("K0_UpSide_dxdzvdydz",dxdz,dydz);
        FillHist("K0_UpSide_cosTvsMomentum_CM",bcost,bpmom);
        FillHist("K0_UpSide_dxdzvdydz_CM",bdxdz,bdydz);
      }
      // Lambda
      if(Lflag){
        FillHist("Lambda_cosTvsMomentum",cost,pmom);
        FillHist("Lambda_dxdzvdydz",dxdz,dydz);
        FillHist("Lambda_cosTvsMomentum_CM",bcost,bpmom);
        FillHist("Lambda_dxdzvdydz_CM",dxdz,dydz);
      }
      if(LLflag){
        FillHist("Lambda_LowSide_cosTvsMomentum",cost,pmom);
        FillHist("Lambda_LowSide_dxdzvdydz",dxdz,dydz);
        FillHist("Lambda_LowSide_cosTvsMomentum_CM",bcost,bpmom);
        FillHist("Lambda_LowSide_dxdzvdydz_CM",dxdz,dydz);
      }
      if(LUflag){
        FillHist("Lambda_UpSide_cosTvsMomentum",cost,pmom);
        FillHist("Lambda_UpSide_dxdzvdydz",dxdz,dydz);
        FillHist("Lambda_UpSide_cosTvsMomentum_CM",bcost,bpmom);
        FillHist("Lambda_UpSide_dxdzvdydz_CM",dxdz,dydz);
      }
      // Missing nuetron
      if(MissNflag){
        FillHist("MissN_cosTvsMomentum",mcost,mpmom);
        FillHist("MissN_dxdzvdydz",mdxdz,mdydz);
        FillHist("MissN_cosTvsMomentum_CM",mbcost,mbpmom);
        FillHist("MissN_dxdzvdydz_CM",mbdxdz,mbdydz);
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

  //  // IH and CDH
  //  const int nhodo = 2;
  //  int NumOfSegments[nhodo] = {24,36};
  //  TString CounterName[nhodo] = {"IH","CDH"};
  //  for(int ihodo=0; ihodo<nhodo; ihodo++){
  //    std::cout << "Define Histograms for " << CounterName[ihodo].Data() << std::endl;
  //    new TH1F( Form("%s_HitPat",CounterName[ihodo].Data()), Form("Hit Pattern %s;%s segment",CounterName[ihodo].Data(),CounterName[ihodo].Data()), NumOfSegments[ihodo]+1, 0, NumOfSegments[ihodo]+1 );
  //    new TH2F( Form("%s_SegmentvPosition",CounterName[ihodo].Data()), Form("Hit Position %s;%s segment;Position (cm)",CounterName[ihodo].Data(),CounterName[ihodo].Data()), NumOfSegments[ihodo]+1, 0, NumOfSegments[ihodo]+1, 201, -50.25, 50.25 );
  //    new TH1F( Form("%s_Mult",CounterName[ihodo].Data()), Form("Multiplicity %s;Multiplicity;Counts",CounterName[ihodo].Data()), NumOfSegments[ihodo]+1, 0, NumOfSegments[ihodo]+1 );
  //  }
  //  for(int seg=1; seg<=NumOfSegments[1]; seg++){
  //    new TH2F( Form("CDH%d_dEMeanvOffset",seg),   Form("dEMean Offset corr. CDH%d;dE (MeV);Offset (ns)",seg),     400,    0, 40,  4000,    -100, 100 );
  //    new TH1F( Form("CDH%d_Offset",seg),   Form("Offset CDH%d;Offset (ns);Counts",seg),    4000,    -100, 100 );
  //  }
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
  //std::cout << "Define Histograms for TOF" << std::endl;
  //for(int ihodo=0; ihodo<nhodo; ihodo++){
  //  new TH2F( Form("T0%s_TOF",CounterName[ihodo].Data()), Form("TOF T0-%s;%s segment;TOF (ns)",CounterName[ihodo].Data(),CounterName[ihodo].Data()), NumOfSegments[ihodo], 1, NumOfSegments[ihodo]+1, 4000, -50, 50 );
  //}
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
  // Angular distribution
  new TH2F( "K0_cosTvsMomentum", "K^{0} cos#theta vs Momentum;cos#theta;Momentum (GeV/c)", 200, -1.0, 1.0, 2000, 0.0, 2.0);
  new TH2F( "K0_dxdzvdydz", "K^{0} dx/dz vs. dy/dz;dx/dz;dydz", 200, -1.0, 1.0, 200, -1.0, 1.0);
  new TH2F( "K0_LowSide_cosTvsMomentum", "K^{0} cos#theta vs Momentum;cos#theta;Momentum (GeV/c)", 200, -1.0, 1.0, 2000, 0.0, 2.0);
  new TH2F( "K0_LowSide_dxdzvdydz", "K^{0} dx/dz vs. dy/dz;dx/dz;dydz", 200, -1.0, 1.0, 200, -1.0, 1.0);
  new TH2F( "K0_UpSide_cosTvsMomentum", "K^{0} cos#theta vs Momentum;cos#theta;Momentum (GeV/c)", 200, -1.0, 1.0, 2000, 0.0, 2.0);
  new TH2F( "K0_UpSide_dxdzvdydz", "K^{0} dx/dz vs. dy/dz;dx/dz;dydz", 200, -1.0, 1.0, 200, -1.0, 1.0);
  new TH2F( "Lambda_cosTvsMomentum", "#Lambda cos#theta vs Momentum;cos#theta;Momentum (GeV/c)", 200, -1.0, 1.0, 2000, 0.0, 2.0);
  new TH2F( "Lambda_dxdzvdydz", "#Lambda dx/dz vs. dy/dz;dx/dz;dydz", 200, -1.0, 1.0, 200, -1.0, 1.0);
  new TH2F( "Lambda_LowSide_cosTvsMomentum", "#Lambda cos#theta vs Momentum;cos#theta;Momentum (GeV/c)", 200, -1.0, 1.0, 2000, 0.0, 2.0);
  new TH2F( "Lambda_LowSide_dxdzvdydz", "#Lambda dx/dz vs. dy/dz;dx/dz;dydz", 200, -1.0, 1.0, 200, -1.0, 1.0);
  new TH2F( "Lambda_UpSide_cosTvsMomentum", "#Lambda cos#theta vs Momentum;cos#theta;Momentum (GeV/c)", 200, -1.0, 1.0, 2000, 0.0, 2.0);
  new TH2F( "Lambda_UpSide_dxdzvdydz", "#Lambda dx/dz vs. dy/dz;dx/dz;dydz", 200, -1.0, 1.0, 200, -1.0, 1.0);
  new TH2F( "MissN_cosTvsMomentum", "Missing n cos#theta vs Momentum;cos#theta;Momentum (GeV/c)", 200, -1.0, 1.0, 2000, 0.0, 2.0);
  new TH2F( "MissN_dxdzvdydz", "Missing n dx/dz vs. dy/dz;dx/dz;dydz", 200, -1.0, 1.0, 200, -1.0, 1.0);
  new TH2F( "K0_cosTvsMomentum_CM", "K^{0} cos#theta vs Momentum;cos#theta;Momentum (GeV/c)", 200, -1.0, 1.0, 2000, 0.0, 2.0);
  new TH2F( "K0_dxdzvdydz_CM", "K^{0} dx/dz vs. dy/dz;dx/dz;dydz", 200, -1.0, 1.0, 200, -1.0, 1.0);
  new TH2F( "K0_LowSide_cosTvsMomentum_CM", "K^{0} cos#theta vs Momentum;cos#theta;Momentum (GeV/c)", 200, -1.0, 1.0, 2000, 0.0, 2.0);
  new TH2F( "K0_LowSide_dxdzvdydz_CM", "K^{0} dx/dz vs. dy/dz;dx/dz;dydz", 200, -1.0, 1.0, 200, -1.0, 1.0);
  new TH2F( "K0_UpSide_cosTvsMomentum_CM", "K^{0} cos#theta vs Momentum;cos#theta;Momentum (GeV/c)", 200, -1.0, 1.0, 2000, 0.0, 2.0);
  new TH2F( "K0_UpSide_dxdzvdydz_CM", "K^{0} dx/dz vs. dy/dz;dx/dz;dydz", 200, -1.0, 1.0, 200, -1.0, 1.0);
  new TH2F( "Lambda_cosTvsMomentum_CM", "#Lambda cos#theta vs Momentum;cos#theta;Momentum (GeV/c)", 200, -1.0, 1.0, 2000, 0.0, 2.0);
  new TH2F( "Lambda_dxdzvdydz_CM", "#Lambda dx/dz vs. dy/dz;dx/dz;dydz", 200, -1.0, 1.0, 200, -1.0, 1.0);
  new TH2F( "Lambda_LowSide_cosTvsMomentum_CM", "#Lambda cos#theta vs Momentum;cos#theta;Momentum (GeV/c)", 200, -1.0, 1.0, 2000, 0.0, 2.0);
  new TH2F( "Lambda_LowSide_dxdzvdydz_CM", "#Lambda dx/dz vs. dy/dz;dx/dz;dydz", 200, -1.0, 1.0, 200, -1.0, 1.0);
  new TH2F( "Lambda_UpSide_cosTvsMomentum_CM", "#Lambda cos#theta vs Momentum;cos#theta;Momentum (GeV/c)", 200, -1.0, 1.0, 2000, 0.0, 2.0);
  new TH2F( "Lambda_UpSide_dxdzvdydz_CM", "#Lambda dx/dz vs. dy/dz;dx/dz;dydz", 200, -1.0, 1.0, 200, -1.0, 1.0);
  new TH2F( "MissN_cosTvsMomentum_CM", "Missing n cos#theta vs Momentum;cos#theta;Momentum (GeV/c)", 200, -1.0, 1.0, 2000, 0.0, 2.0);
  new TH2F( "MissN_dxdzvdydz_CM", "Missing n dx/dz vs. dy/dz;dx/dz;dydz", 200, -1.0, 1.0, 200, -1.0, 1.0);
}
