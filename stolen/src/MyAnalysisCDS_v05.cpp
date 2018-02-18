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
    CDSTrack* cdc = cdstrackMan->GoodTrack(cdstrackMan->GoodTrackID(it));
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
  // ##################### //
  // CDS Particle analysis //
  // ##################### //
  int ncds[6] = {0,0,0,0,0,0};
  ncds[0] = particle->nPiplus();
  ncds[1] = particle->nPiminus();
  ncds[2] = particle->nKaon();
  ncds[3] = particle->nProton();
  ncds[4] = particle->nDeuteron();
  ncds[5] = particle->nTriton();
  ncds[5] += particle->nHelium3();
  ncds[5] += particle->nOther();
  TString cdsstr[6] = {"#pi^{+}","#pi^{-}","K^{-}","p","d","Other"};
  FillHist("CDS_NumOfParticle",particle->nCDS());
  for(int it=0; it<particle->nCDS(); it++){
    pCDS* cds = particle->cdsi(it);
    if(GeomTools::GetID(cds->vbeam())==CID_Fiducial){
      double cds_vdis = cds->vdis();
      TVector3 vtxb = cds->vbeam();
      TVector3 ZeroV;
      TLorentzVector Ltgt;
      if(dlist->GetMaterial(CID_Target)=="LHydrogen"){
        Ltgt.SetVectM(ZeroV,pMass);
      }
      if(dlist->GetMaterial(CID_Target)=="LDeuterium"){
        Ltgt.SetVectM(ZeroV,dMass);
      }
      if(dlist->GetMaterial(CID_Target)=="LHelium-3"){
        Ltgt.SetVectM(ZeroV,ThreeHeMass);
      }
      TLorentzVector Lbeam = beam->GetLorentzVector(vtxb);
      TVector3 boost = (Ltgt+Lbeam).BoostVector();
      TLorentzVector cLtgt = Ltgt; cLtgt.Boost(-boost); TLorentzVector cLbeam = Lbeam; cLbeam.Boost(-boost);
      TLorentzVector Lcds = cds->GetLorentzVector();
      TLorentzVector cLcds = cds->GetLorentzVector(); cLcds.Boost(-boost);
      TString frame[2] = {"Lab","CM"};
      TString pname[8] = {"pip","p","d","t","he","pim","k","o"};
      double cds_mass[2]  = { (Lcds).M()                    , (cLcds).M()};
      double cds_mom[2]   = { (Lcds).Vect().Mag()           , (cLcds).Vect().Mag()};
      double cds_cost[2]  = { (Lcds).Vect().CosTheta()      , (cLcds).Vect().CosTheta()};
      double cds_phi[2]   = { (Lcds).Vect().Phi()           , (cLcds).Vect().Phi()};
      double cds_mmass[2]   = { (Ltgt+Lbeam-Lcds).M()              , (cLtgt+cLbeam-cLcds).M()};
      double cds_mmom[2]    = { (Ltgt+Lbeam-Lcds).Vect().Mag()     , (cLtgt+cLbeam-cLcds).Vect().Mag()};
      double cds_mcost[2]   = { (Ltgt+Lbeam-Lcds).Vect().CosTheta(), (cLtgt+cLbeam-cLcds).Vect().CosTheta()};

      FillHist(Form("%s_Vertex_XY",pname[cds->pid()].Data()), vtxb.X(), vtxb.Y());
      FillHist(Form("%s_Vertex_ZX",pname[cds->pid()].Data()), vtxb.Z(), vtxb.X());
      FillHist(Form("%s_Vertex_ZY",pname[cds->pid()].Data()), vtxb.Z(), vtxb.Y());
      FillHist(Form("%s_Vertex_Z",pname[cds->pid()].Data()), vtxb.Z());
      FillHist(Form("%s_VDIS",pname[cds->pid()].Data()),cds_vdis);
      FillHist(Form("All_VDIS"),cds_vdis);
      for(int f=1; f<2; f++){
        FillHist(Form("%s_IMvsMomentum_%s",pname[cds->pid()].Data(),frame[f].Data()),cds_mass[f],cds_mom[f]);
        FillHist(Form("%s_IMvsCosT_%s",pname[cds->pid()].Data(),frame[f].Data()),cds_mass[f],cds_cost[f]);
        FillHist(Form("%s_MMvsMomentum_%s",pname[cds->pid()].Data(),frame[f].Data()),cds_mmass[f],cds_mmom[f]);
        FillHist(Form("%s_MMvsCosT_%s",pname[cds->pid()].Data(),frame[f].Data()),cds_mmass[f],cds_mcost[f]);
      }
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
  TString pname[8]  = {"pip","p","d","t","he","pim","k","o"};
  TString pname2[8] = {"#pi^{+}","p","d","t","he","#pi^{-}","K^{-}","o"};
  for(int ip=0; ip<8; ip++){
    new TH1F( Form("%s_VDIS",pname[ip].Data()), Form("VDIS of %s;VDIS (cm);Coutns",pname2[ip].Data()), 100, 0, 10);
    new TH2F( Form("%s_Vertex_XY",pname[ip].Data()), Form("Vertex XY plane (%s);x (cm);y (cm)",pname2[ip].Data()), 600, -15, 15, 600, -15, 15 );
    new TH2F( Form("%s_Vertex_ZX",pname[ip].Data()), Form("Vertex ZX plane (%s);z (cm);x (cm)",pname2[ip].Data()), 1200, -30, 30, 600, -15, 15 );
    new TH2F( Form("%s_Vertex_ZY",pname[ip].Data()), Form("Vertex ZY plane (%s);z (cm);y (cm)",pname2[ip].Data()), 1200, -30, 30, 600, -15, 15 );
    new TH1F( Form("%s_Vertex_Z",pname[ip].Data()), Form("Vertex Z axis (%s);z (cm);Counts",pname2[ip].Data()), 1200, -30, 30 );
    new TH2F( Form("%s_IMvsMomentum_CM",pname[ip].Data()), Form("Invariant mass (%s) vs Momentum;IM(%s) (GeV/c^{2});Momentum (GeV/c)",pname2[ip].Data(),pname2[ip].Data()), 5000, 0.0, 5.0, 2000, 0.0, 2.0 );
    new TH2F( Form("%s_IMvsCosT_CM",pname[ip].Data()), Form("Invariant mass (%s) vs cos#theta;IM(%s) (GeV/c^{2});cos#theta",pname2[ip].Data(),pname2[ip].Data()), 5000, 0.0, 5.0, 2000, -1.0, 1.0 );
    new TH2F( Form("%s_MMvsMomentum_CM",pname[ip].Data()), Form("Missing mass (%s) vs Momentum;MM(%s) (GeV/c^{2});Momentum (GeV/c)",pname2[ip].Data(),pname2[ip].Data()), 5000, 0.0, 5.0, 2000, 0.0, 2.0 );
    new TH2F( Form("%s_MMvsCosT_CM",pname[ip].Data()), Form("Missing mass (%s) vs cos#theta;MM(%s) (GeV/c^{2});cos#theta",pname2[ip].Data(),pname2[ip].Data()), 5000, 0.0, 5.0, 2000, -1.0, 1.0 );
  }
  new TH1F( Form("All_VDIS"), Form("VDIS of All;VDIS (cm);Coutns"), 100, 0, 10);
  new TH1F( "CDS_NumOfParticle", "Number of CDS tracks", 10, 0, 10);
  TH2F* h2;
  h2 = new TH2F("CDS_Particle","Detected particles by CDS;Particle;Number of all tracks",6,0,6,6,0,6);
  h2->GetYaxis()->SetBinLabel(1,"1 track");
  h2->GetYaxis()->SetBinLabel(2,"2 track");
  h2->GetYaxis()->SetBinLabel(3,"3 track");
  h2->GetYaxis()->SetBinLabel(4,"4 track");
  h2->GetYaxis()->SetBinLabel(5,"5 track");
  h2->GetYaxis()->SetBinLabel(6,"6 track");
  h2->GetXaxis()->SetBinLabel(1,"#pi^{+}");
  h2->GetXaxis()->SetBinLabel(2,"#pi^{-}");
  h2->GetXaxis()->SetBinLabel(3,"K^{-}");
  h2->GetXaxis()->SetBinLabel(4,"p");
  h2->GetXaxis()->SetBinLabel(5,"d");
  h2->GetXaxis()->SetBinLabel(6,"Other");


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
