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

void MyAnalysisCDS::Clear()
{
}

bool MyAnalysisCDS::DoAnalysis(ConfMan* conf, EventHeader* header, BeamLineHitMan* blMan, BeamLineTrackMan* bltrackMan, CDSHitMan* cdsMan, CDSTrackingMan* cdstrackMan, Particle* particle)
{
  DetectorList *dlist=DetectorList::GetInstance();
  FillHist("CDC_nGoodTrack",cdstrackMan->nGoodTrack());
  if(particle->nBeam()!=1) return false;
  pBeam* beam = particle->beam(0);

  // ###################### //
  // CDS all track analysis //
  // ###################### //
  for( int it=0; it<cdstrackMan->nGoodTrack(); it++ ){
    CDSTrack* cdc = cdstrackMan->Track(cdstrackMan->GoodTrackID(it));
    bool ELOSS = true;
    pCDS* cdstrack = TrackTools::CalcSingleAll(beam,cdstrackMan->GoodTrackID(it),cdc,cdsMan,ELOSS);
    if(cdstrack!=0){
      particle->AddCDS(*cdstrack);
    }
  }

  for(int it=0; it<particle->nCDS(); it++){
    pCDS* cdstrack = particle->cdsi(it);
    double mass2 = cdstrack->mass2();
    double beta = cdstrack->beta();
    double momentum = cdstrack->mom();
    TVector3 vtxb = cdstrack->vbeam();
    TVector3 vtxh = cdstrack->vcdc();
    FillHist("Vertex_XY", vtxb.X(), vtxb.Y());
    FillHist("Vertex_ZX", vtxb.Z(), vtxb.X());
    FillHist("Vertex_ZY", vtxb.Z(), vtxb.Y());
    FillHist("Vertex_Z", vtxb.Z());
    /* CDH */
    for(int icdh=0; icdh<cdstrack->ncdh(); icdh++){
      int seg = cdstrack->cdhseg(icdh);
      double toffs = cdstrack->dt();
      for(int i=0; i<cdsMan->nCDH(); i++){
        HodoscopeLikeHit* hit = cdsMan->CDH(i);
        if(hit->CheckRange()){
          if(hit->seg()==seg){
            double emean = hit->emean();
            double ctmean = hit->ctmean();
            const char* ud[2]={"u","d"};
            for(int iud=0; iud<2; iud++){
              FillHist(Form("CDH%s%d_TimeOffset",ud[iud],seg),toffs);
              FillHist(Form("CDH%s%d_dEMeanvsTimeOffset",ud[iud],seg),emean,toffs);
              FillHist(Form("CDH%s%d_dEMeanvsTimeOffset",ud[iud],seg),emean,toffs);
            }
            FillHist(Form("CDH%d_TimeOffset",seg),toffs);
            FillHist(Form("CDH%d_dEMeanvsTimeOffset",seg),emean,toffs);
          }
        }
      }
    }

    if(GeomTools::GetID(vtxb)==CID_Fiducial){
      FillHist("TOFvsMomentum",cdstrack->tof()-beam->t0time(),momentum);
      FillHist("Mass2vsMomentum", mass2, momentum);
      FillHist("OverbetavsMomentum", 1./beta, momentum);
      FillHist("VertexwithTRG_XY", vtxb.X(), vtxb.Y());
      FillHist("VertexwithTRG_ZX", vtxb.Z(), vtxb.X());
      FillHist("VertexwithTRG_ZY", vtxb.Z(), vtxb.Y());
      FillHist("VertexwithTRG_Z", vtxb.Z());
    }
  }

  // ###################### //
  // CDS 2 track analysis //
  // ###################### //
  for( int it1=0; it1<particle->nCDS(); it1++ ){
    for( int it2=it1+1; it2<particle->nCDS(); it2++ ){
      pCDS* cds1 = particle->cdsi(it1);
      pCDS* cds2 = particle->cdsi(it2);
      bool ELOSS = true;
      pCDS* product = TrackTools::Calc2HelixAll(beam,cds1,cds2,cdstrackMan,ELOSS);
      if(product!=0){
        particle->AddProduct(*product);
      }
    }
  }

  for(int it=0; it<particle->nProduct(); it++){
    pCDS* product = particle->product(it);
    double mass = product->mass();
    double mom = product->mom();
    TVector3 vtxb = product->vbeam();
    TVector3 vtx = product->vertex();
    int comb = product->comb();
    FillHist("PVertex_XY", vtxb.X(), vtxb.Y());
    FillHist("PVertex_ZX", vtxb.Z(), vtxb.X());
    FillHist("PVertex_ZY", vtxb.Z(), vtxb.Y());
    FillHist("PVertex_Z", vtxb.Z());
    if(GeomTools::GetID(vtxb)==CID_Fiducial){
      FillHist("PVertexwithTRG_XY", vtxb.X(), vtxb.Y());
      FillHist("PVertexwithTRG_ZX", vtxb.Z(), vtxb.X());
      FillHist("PVertexwithTRG_ZY", vtxb.Z(), vtxb.Y());
      FillHist("PVertexwithTRG_Z", vtxb.Z());

      if(comb==pow(2,CDS_PiPlus)+pow(2,CDS_PiMinus)){
        FillHist("IM_pipi",mass);
        FillHist("Momentum_pipi",mom);
      }
      if(comb==pow(2,CDS_Proton)+pow(2,CDS_PiMinus)){
        FillHist("IM_pip",mass);
        FillHist("Momentum_pip",mom);
      }
      if(comb==pow(2,CDS_PiPlus)+pow(2,CDS_Kaon)){
        FillHist("IM_piK",mass);
        FillHist("Momentum_piK",mom);
      }
      if(comb==pow(2,CDS_Proton)+pow(2,CDS_Kaon)){
        FillHist("IM_pK",mass);
        FillHist("Momentum_pK",mom);
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
  DetectorList *dlist=DetectorList::GetInstance();

  // CDH
  //std::cout << "Define Histograms for Hodoscopes" << std::endl;
  //const int nhodo = 1;
  //int hodoid[nhodo]={
  //  CID_CDH
  //};
  //for(int ihodo=0; ihodo<nhodo; ihodo++){
  //  const int cid = hodoid[ihodo];
  //  const char* name = dlist -> GetName(cid).c_str();
  //  const int nsegs = dlist -> GetNsegs(cid);
  //  for( int seg=1; seg<=nsegs; seg++ ){
  //    const char* ud[2] = {"u","d"};
  //    const char* udname[2] = {"up","down"};
  //    for(int iud=0; iud<2; iud++){
  //      new TH1F( Form("%s%s%d_TimeOffset",name,ud[iud],seg),   Form("Time offset %s-%d(%s);Time offset (ns);Counts",name,seg,udname[iud]),    2000,    -100, 100 );
  //      new TH2F( Form("%s%s%d_dEMeanvsTimeOffset",name,ud[iud],seg),   Form("dEMean CTime corr. %s-%d(%s);dE (MeV);Time offset (ns)",name,seg,udname[iud]),     200,    0, 20,  2000,    -100, 100 );
  //    }
  //    new TH1F( Form("%s%d_TimeOffset",name,seg),   Form("Time offset %s-%d;Time offset (ns);Counts",name,seg),    3000,    -50, 100 );
  //    new TH2F( Form("%s%d_dEMeanvsTimeOffset",name,seg),   Form("dEMean Time offset corr. %s-%d;dE (MeV);Time offset (ns)",name,seg),     200,    0, 20,  2000,    -100, 100 );
  //  }
  //}

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
  // CDS Particles
  std::cout << "Define Histograms for CDS particles" << std::endl;
  new TH2F( "Mass2vsMomentum", "Mass^{2} vs. Momentum;Mass^{2} (GeV/c^{2})^{2};Momentum (GeV/c)", 1500, -5, 10, 1000, -2, 18 );
  new TH2F( "OverbetavsMomentum", "1/#beta vs. Momentum;1/#beta;Momentum (GeV/c)", 1000, -0, 10, 500, -5, 5 );
  new TH2F( "TOFvsMomentum", "TOF vs. Momentum;TOF (ns);Momentum (GeV/c)", 1000, 0, 20, 500, -5, 5 );

  // Vertex
  std::cout << "Define Histograms for Vertex" << std::endl;
  new TH2F( "Vertex_XY", "Vertex XY plane;x (cm);y (cm)", 600, -15, 15, 600, -15, 15 );
  new TH2F( "Vertex_ZX", "Vertex ZX plane;z (cm);x (cm)", 1200, -30, 30, 600, -15, 15 );
  new TH2F( "Vertex_ZY", "Vertex ZY plane;z (cm);y (cm)", 1200, -30, 30, 600, -15, 15 );
  new TH1F( "Vertex_Z", "Vertex Z axis;z (cm);Counts", 1200, -30, 30 );
  new TH2F( "VertexwithTRG_XY", "Vertex XY plane;x (cm);y (cm)", 600, -15, 15, 600, -15, 15 );
  new TH2F( "VertexwithTRG_ZX", "Vertex ZX plane;z (cm);x (cm)", 1200, -30, 30, 600, -15, 15 );
  new TH2F( "VertexwithTRG_ZY", "Vertex ZY plane;z (cm);y (cm)", 1200, -30, 30, 600, -15, 15 );
  new TH1F( "VertexwithTRG_Z", "Vertex Z axis;z (cm);Counts", 1200, -30, 30 );

  new TH2F( "PVertex_XY", "Vertex XY plane;x (cm);y (cm)", 600, -15, 15, 600, -15, 15 );
  new TH2F( "PVertex_ZX", "Vertex ZX plane;z (cm);x (cm)", 1200, -30, 30, 600, -15, 15 );
  new TH2F( "PVertex_ZY", "Vertex ZY plane;z (cm);y (cm)", 1200, -30, 30, 600, -15, 15 );
  new TH1F( "PVertex_Z", "Vertex Z axis;z (cm);Counts", 1200, -30, 30 );
  new TH2F( "PVertexwithTRG_XY", "Vertex XY plane;x (cm);y (cm)", 600, -15, 15, 600, -15, 15 );
  new TH2F( "PVertexwithTRG_ZX", "Vertex ZX plane;z (cm);x (cm)", 1200, -30, 30, 600, -15, 15 );
  new TH2F( "PVertexwithTRG_ZY", "Vertex ZY plane;z (cm);y (cm)", 1200, -30, 30, 600, -15, 15 );
  new TH1F( "PVertexwithTRG_Z", "Vertex Z axis;z (cm);Counts", 1200, -30, 30 );

  // Invariant mass of product
  std::cout << "Define Histograms for product I.M." << std::endl;

  new TH1F( "IM_pipi", "Invariant mass (#pi^{-}#pi^{+});IM(#pi^{-}#pi^{+}) (GeV/c^{2});Counts", 2000, 0.0, 2.0 );
  new TH1F( "Momentum_pipi", "Momentum (#pi^{-}#pi^{+});Momentum (#pi^{-}#pi^{+}) (GeV/c);Counts", 2000, 0.0, 2.0 );
  new TH1F( "IM_pip", "Invariant mass (#pi^{-}p);IM(#pi^{-}p) (GeV/c^{2});Counts", 2000, 0.0, 2.0 );
  new TH1F( "Momentum_pip", "Momentum (#pi^{-}p);Momentum (#pi^{-}p) (GeV/c);Counts", 2000, 0.0, 2.0 );
  new TH1F( "IM_piK", "Invariant mass (#pi^{+}K^{-});IM(#pi^{+}K^{-}) (GeV/c^{2});Counts", 2000, 0.0, 2.0 );
  new TH1F( "Momentum_piK", "Momentum (#pi^{+}K^{-});Momentum (#pi^{+}K^{-}) (GeV/c);Counts", 2000, 0.0, 2.0 );
  new TH1F( "IM_pK", "Invariant mass (pK^{-});IM(pK^{-}) (GeV/c^{2});Counts", 2000, 0.0, 2.0 );
  new TH1F( "Momentum_pK", "Momentum (pK^{-});Momentum (pK^{-}) (GeV/c);Counts", 2000, 0.0, 2.0 );
}
