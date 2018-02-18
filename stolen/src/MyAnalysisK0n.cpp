// MyAnalysisK0n.cpp

#include "MyAnalysisK0n.h"

MyAnalysisK0n::MyAnalysisK0n(TFile* rt, ConfMan* conf)
{
  Initialize(conf);
  CutCondition();
  Clear();
}

MyAnalysisK0n::~MyAnalysisK0n()
{
  Clear();
  rtFile->cd();
  rtFile->Write();
  rtFile->Close();
}

void MyAnalysisK0n::Clear()
{
}

bool MyAnalysisK0n::DoAnalysis(ConfMan* conf, EventHeader* header, BeamLineHitMan* blMan, BeamLineTrackMan* bltrackMan, CDSHitMan* cdsMan, CDSTrackingMan* cdstrackMan, Particle* particle)
{

  rtFile->cd();

  DetectorList *dlist=DetectorList::GetInstance();
  if(particle->nBeam()!=1) return false;
  pBeam* beam = particle->beam(0);

  /* Event selection */
  //if(cdstrackMan->nGoodTrack()!=3) return true;
  //int MulCDH = 0;
  //for(int i=0; i<cdsMan->nCDH(); i++){
  //  if(cdsMan->CDH(i)->CheckRange()){
  //    MulCDH++;
  //  }
  //}
  //if(MulCDH!=3){
  //  return true;
  //}
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
  int MulPC=0;
  for(int i=0; i<blMan->nPC(); i++){
    HodoscopeLikeHit* hit = blMan->PC(i);
    if(hit->CheckRange()) MulPC++;
  }
  if(MulBVC!=0||MulCVC!=0||MulPC!=0) return true;

  // ###################### //
  // CDS 2 track analysis //
  // ###################### //
  int fk0 = -1;
  for(int it=0; it<particle->nProduct(); it++){
    pCDS* cds = particle->product(it);
    pCDS* cds1 = particle->cdsi(cds->daughter1());
    pCDS* cds2 = particle->cdsi(cds->daughter2());
    if(cds1->chi()>30||cds2->chi()>30){ continue; }
    if(cds->comb()!=pow(2,CDS_PiPlus)+pow(2,CDS_PiMinus)){ continue; }
    if(k0ll<cds->mass()&&cds->mass()<k0ul){
      fk0 = it;
    }
  }
  pCDS* k0 = 0;
  if(fk0!=-1) k0 = particle->product(fk0);
  if(k0==0) return true;

  // ############### //
  // Vertex decision //
  // ############### //
  TVector3 vertex;
  double vdis = 9999;
  bool FIDUCIAL = false;
  vdis = k0->pbdca();
  vertex = k0->vbeam();
  if(GeomTools::GetID(vertex)==CID_Fiducial){
    FIDUCIAL = true;
  }
  if(!FIDUCIAL){ return true; }

  // ############### //
  // NC hit decision //
  // ############### //
  double time = 9999;
  int fnc = -1;
  for(int it=0; it<particle->nNC(); it++){
    pNC* nc = particle->nc(it);
    nc->CalcMom(beam,vertex);
    FillHist("FWD_OverbetavsMomentum",1.0/nc->beta(),nc->mom().Mag());
    FillHist("FWD_OverbetavsEnergy",1.0/nc->beta(),nc->energy());
    FillHist("FWD_TOFvsMomentum",nc->tof(),nc->mom().Mag());
    FillHist("FWD_HitPosition",nc->hitpos().X(),nc->hitpos().Y());
    double seg = (nc->seg()-1)%16+1;
    double lay = (nc->seg()-1)/16+1;
    FillHist("FWD_HitSegment",seg,lay);
    if(nc->pid()==F_Neutron) {
      if(nc->energy()>8.0 && (time>9998||time>nc->time())){
        time = nc->time();
        fnc = it;
      }
    }
  }
  pNC*  nc =0;
  if(fnc!=-1)  nc = particle->nc(fnc);
  if(nc==0) return true;

  // Target
  TVector3 ZeroV;
  TLorentzVector Ltgt; Ltgt.SetVectM(ZeroV,ThreeHeMass);
  TLorentzVector Lbeam = beam->GetLorentzVector(vertex);
  TVector3 boost = (Ltgt+Lbeam).BoostVector();
  TLorentzVector cLtgt = Ltgt; cLtgt.Boost(-boost); TLorentzVector cLbeam = Lbeam; cLbeam.Boost(-boost);
  // K0
  TLorentzVector Lk0 = k0->GetLorentzVector();
  TLorentzVector cLk0 = k0->GetLorentzVector(); cLk0.Boost(-boost);
  // K0 + Proton
  TLorentzVector Lcds = Lk0;
  TLorentzVector cLcds = cLk0;
  // Neutron
  TLorentzVector Ln  = nc ->GetLorentzVector();
  TLorentzVector cLn  = nc ->GetLorentzVector(); cLn.Boost(-boost);

  //double k0_mass[2]  = { (Lk0).M()                    , (cLk0).M()};
  //double k0_mom[2]   = { (Lk0).Vect().Mag()           , (cLk0).Vect().Mag()};
  //double k0_cost[2]  = { (Lk0).Vect().CosTheta()      , (cLk0).Vect().CosTheta()};
  //double k0_phi[2]   = { (Lk0).Vect().Phi()           , (cLk0).Vect().Phi()};
  //double k0_mmass[2]   = { (Ltgt+Lbeam-Lk0).M()              , (cLtgt+cLbeam-cLk0).M()};
  //double k0_mmom[2]    = { (Ltgt+Lbeam-Lk0).Vect().Mag()     , (cLtgt+cLbeam-cLk0).Vect().Mag()};
  //double k0_mcost[2]   = { (Ltgt+Lbeam-Lk0).Vect().CosTheta(), (cLtgt+cLbeam-cLk0).Vect().CosTheta()};
  //double nc_mom[2]     = { Ln.Vect().Mag()                  , cLn.Vect().Mag()};
  //double nc_cost[2]    = { Ln.Vect().CosTheta()             , cLn.Vect().CosTheta()};
  //double nc_phi[2]     = { Ln.Vect().Phi()                  , cLn.Vect().Phi()};
  //double nc_mmass[2]   = { (Ltgt+Lbeam-Ln).M()              , (cLtgt+cLbeam-cLn).M()};
  //double nc_mmom[2]    = { (Ltgt+Lbeam-Ln).Vect().Mag()     , (cLtgt+cLbeam-cLn).Vect().Mag()};
  //double nc_mcost[2]   = { (Ltgt+Lbeam-Ln).Vect().CosTheta(), (cLtgt+cLbeam-cLn).Vect().CosTheta()};
  //double cds_mass[2]  = { (Lcds).M()                    , (cLcds).M()};
  //double cds_mom[2]   = { (Lcds).Vect().Mag()           , (cLcds).Vect().Mag()};
  //double cds_cost[2]  = { (Lcds).Vect().CosTheta()      , (cLcds).Vect().CosTheta()};
  //double cds_phi[2]   = { (Lcds).Vect().Phi()           , (cLcds).Vect().Phi()};
  //double cds_mmass[2]   = { (Ltgt+Lbeam-Lcds).M()              , (cLtgt+cLbeam-cLcds).M()};
  //double cds_mmom[2]    = { (Ltgt+Lbeam-Lcds).Vect().Mag()     , (cLtgt+cLbeam-cLcds).Vect().Mag()};
  //double cds_mcost[2]   = { (Ltgt+Lbeam-Lcds).Vect().CosTheta(), (cLtgt+cLbeam-cLcds).Vect().CosTheta()};
  //double ncds_mass[2]  = { (Lcds+Ln).M()                           , (cLcds+cLn).M()};
  //double ncds_mom[2]   = { (Lcds+Ln).Vect().Mag()                  , (cLcds+cLn).Vect().Mag()};
  //double ncds_cost[2]  = { (Lcds+Ln).Vect().CosTheta()             , (cLcds+cLn).Vect().CosTheta()};
  //double ncds_phi[2]   = { (Lcds+Ln).Vect().Phi()                  , (cLcds+cLn).Vect().Phi()};
  double ncds_mmass[2] = { (Ltgt+Lbeam-(Lcds+Ln)).M()              , (cLtgt+cLbeam-(cLcds+cLn)).M()};
  double ncds_mmom[2]  = { (Ltgt+Lbeam-(Lcds+Ln)).Vect().Mag()     , (cLtgt+cLbeam-(cLcds+cLn)).Vect().Mag()};
  double ncds_mcost[2] = { (Ltgt+Lbeam-(Lcds+Ln)).Vect().CosTheta(), (cLtgt+cLbeam-(cLcds+cLn)).Vect().CosTheta()};

  FillHist(Form("k0n_MMass"),ncds_mmass[1]);
  FillHist(Form("k0n_MMom"),ncds_mmom[1]);
  FillHist(Form("k0n_MCosT"),ncds_mcost[1]);

  return true;

}

bool MyAnalysisK0n::FillHist(TString name, double val1, int weight)
{
  TH1F* h1 = (TH1F*)gFile -> Get(name);
  if(h1&&weight>0){ 
    for(int i=0; i<weight; i++){
      h1 -> Fill(val1,1);
    }
    return true;
  }
  else 
    return false;
}

bool MyAnalysisK0n::FillHist(TString name, TString val1, int weight)
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

bool MyAnalysisK0n::FillHist(TString name, double val1, double val2, int weight)
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

bool MyAnalysisK0n::FillHist(TString name, TString val1, TString val2, int weight)
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

void MyAnalysisK0n::CutCondition()
{
  /* K0 mass cut */ 
  double mean  = 0.4972;
  double sigma = 0.0060;
  k0ll = mean-2*sigma; k0ul = mean+2*sigma;
  sblk0ll = mean-8*sigma; sblk0ul = mean-4*sigma;
  sbuk0ll = mean+4*sigma; sbuk0ul = mean+8*sigma;
  /* Lambda mass cut */ 
  mean  = 1.11543;
  sigma = 0.0020;
  lamll = mean-2*sigma; lamul = mean+2*sigma;
  sbllamll = mean-8*sigma; sbllamul = mean-4*sigma;
  sbulamll = mean+4*sigma; sbulamul = mean+8*sigma;
}

bool MyAnalysisK0n::Initialize(ConfMan* confMan)
{
  std::cout << "### MyAnalysisK0n::Initialize ###" << std::endl;

  std::string ofname = confMan->GetOutFileName();
  ofname.insert(ofname.find(".root"),"_anaK0n");

  rtFile =  new TFile( Form("%s",ofname.c_str()), "RECREATE");
  rtFile -> cd();

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
  new TH2F( "CDS_Mass2vsMomentum", "Mass^{2} vs. Momentum;Mass^{2} (GeV/c^{2})^{2};Momentum (GeV/c)", 3000, -5, 10, 800, -2, 2 );
  new TH2F( "CDS_Mass2vsMomentum_AfterFindMass", "Mass^{2} vs. Momentum;Mass^{2} (GeV/c^{2})^{2};Momentum (GeV/c)", 3000, -5, 10, 800, -2, 2 );
  new TH2F( "CDS_Mass2vsMomentum_PiPlus", "Mass^{2} vs. Momentum;Mass^{2} (GeV/c^{2})^{2};Momentum (GeV/c)", 3000, -5, 10, 800, -2, 2 );
  new TH2F( "CDS_Mass2vsMomentum_PiMinus", "Mass^{2} vs. Momentum;Mass^{2} (GeV/c^{2})^{2};Momentum (GeV/c)", 3000, -5, 10, 800, -2, 2 );
  new TH2F( "CDS_Mass2vsMomentum_Kaon", "Mass^{2} vs. Momentum;Mass^{2} (GeV/c^{2})^{2};Momentum (GeV/c)", 3000, -5, 10, 800, -2, 2 );
  new TH2F( "CDS_Mass2vsMomentum_Proton", "Mass^{2} vs. Momentum;Mass^{2} (GeV/c^{2})^{2};Momentum (GeV/c)", 3000, -5, 10, 800, -2, 2 );
  new TH2F( "CDS_Mass2vsMomentum_Deuteron", "Mass^{2} vs. Momentum;Mass^{2} (GeV/c^{2})^{2};Momentum (GeV/c)", 3000, -5, 10, 800, -2, 2 );
  new TH2F( "CDS_Mass2vsMomentum_Other", "Mass^{2} vs. Momentum;Mass^{2} (GeV/c^{2})^{2};Momentum (GeV/c)", 3000, -5, 10, 800, -2, 2 );
  new TH2F( "CDS_Mass2vsMomentum_All", "Mass^{2} vs. Momentum;Mass^{2} (GeV/c^{2})^{2};Momentum (GeV/c)", 3000, -5, 10, 800, -2, 2 );
  new TH2F( "CDS_Mass2vsMomentum_AfterELossCorr", "Mass^{2} vs. Momentum;Mass^{2} (GeV/c^{2})^{2};Momentum (GeV/c)", 3000, -5, 10, 800, -2, 2 );
  new TH2F( "CDS_Mass2vsMomentum_Fiducial", "Mass^{2} vs. Momentum;Mass^{2} (GeV/c^{2})^{2};Momentum (GeV/c)", 3000, -5, 10, 800, -2, 2 );
  new TH2F( "CDS_OverbetavsMomentum", "1/#beta vs. Momentum;1/#beta;Momentum (GeV/c)", 1000, 0, 10, 800, -2, 2 );
  new TH2F( "CDS_OverbetavsMomentum_AfterFindMass", "1/#beta vs. Momentum;1/#beta;Momentum (GeV/c)", 1000, 0, 10, 800, -2, 2 );
  new TH2F( "CDS_OverbetavsMomentum_PiPlus", "1/#beta vs. Momentum;1/#beta;Momentum (GeV/c)", 1000, 0, 10, 800, -2, 2 );
  new TH2F( "CDS_OverbetavsMomentum_PiMinus", "1/#beta vs. Momentum;1/#beta;Momentum (GeV/c)", 1000, 0, 10, 800, -2, 2 );
  new TH2F( "CDS_OverbetavsMomentum_Kaon", "1/#beta vs. Momentum;1/#beta;Momentum (GeV/c)", 1000, 0, 10, 800, -2, 2 );
  new TH2F( "CDS_OverbetavsMomentum_Proton", "1/#beta vs. Momentum;1/#beta;Momentum (GeV/c)", 1000, 0, 10, 800, -2, 2 );
  new TH2F( "CDS_OverbetavsMomentum_Deuteron", "1/#beta vs. Momentum;1/#beta;Momentum (GeV/c)", 1000, 0, 10, 800, -2, 2 );
  new TH2F( "CDS_OverbetavsMomentum_Other", "1/#beta vs. Momentum;1/#beta;Momentum (GeV/c)", 1000, 0, 10, 800, -2, 2 );
  new TH2F( "CDS_OverbetavsMomentum_All", "1/#beta vs. Momentum;1/#beta;Momentum (GeV/c)", 1000, 0, 10, 800, -2, 2 );
  new TH2F( "CDS_OverbetavsMomentum_AfterELossCorr", "1/#beta vs. Momentum;1/#beta;Momentum (GeV/c)", 1000, 0, 10, 800, -2, 2 );
  new TH2F( "CDS_OverbetavsMomentum_Fiducial", "1/#beta vs. Momentum;1/#beta;Momentum (GeV/c)", 1000, 0, 10, 800, -2, 2 );
  TString pname[4]  = {"k0n","k0","sblk0","sbuk0"};
  TString pname2[4] = {"K^{0}n","K^{0}","K^{0}_{lower}","K^{0}_{upper}"};
  for(int ip=0; ip<1; ip++){
    new TH1F( Form("%s_MMass",pname[ip].Data()), Form("Missing Mass (%s);MM(%s) (GeV/c^{2});Counts",pname2[ip].Data(),pname2[ip].Data()), 4000, 1.0, 3.0);
    new TH1F( Form("%s_MMomentum",pname[ip].Data()), Form("Missing Momentum (%s);Missing Momentum (%s) (GeV/c);Counts",pname2[ip].Data(),pname2[ip].Data()), 4000, 0.0, 2.0);
    new TH1F( Form("%s_MCosT",pname[ip].Data()), Form("cos#theta (%s);cos#theta(%s);Counts",pname2[ip].Data(),pname2[ip].Data()), 3000, -1.5, 1.5 );
  }

  return true;
}
