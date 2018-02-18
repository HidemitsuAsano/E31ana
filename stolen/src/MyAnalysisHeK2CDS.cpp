// MyAnalysisHeK2CDS.cpp

#include "MyAnalysisHeK2CDS.h"

MyAnalysisHeK2CDS::MyAnalysisHeK2CDS(TFile* rt, ConfMan* conf)
{
  Initialize(conf);
  CutCondition();
  Clear();
}

MyAnalysisHeK2CDS::~MyAnalysisHeK2CDS()
{
  Clear();
  rtFile->cd();
  rtFile->Write();
  rtFile->Close();
}

void MyAnalysisHeK2CDS::Clear()
{
}

bool MyAnalysisHeK2CDS::DoAnalysis(ConfMan* conf, EventHeader* header, BeamLineHitMan* blMan, BeamLineTrackMan* bltrackMan, CDSHitMan* cdsMan, CDSTrackingMan* cdstrackMan, Particle* particle)
{
  rtFile->cd();

  FillHist("EventNumber",0); /* All Events */
  DetectorList *dlist=DetectorList::GetInstance();
  if(particle->nBeam()!=1) return false;
  pBeam* beam = particle->beam(0);
  double beamtof = beam->bhdt0tof();

  /* Event selection */
  if(particle->nCDS()!=2) return true;
  FillHist("EventNumber",1); /* 2 hits in CDS*/
  if(particle->nProduct()!=1) return true;
  FillHist("EventNumber",2); /* 1 pair in CDS */

  // ############### //
  // Vertex decision //
  // ############### //
  TVector3 vertex;
  double vdis = 9999;
  int fcds = -1;
  for(int it=0; it<particle->nProduct(); it++){
    pCDS* cds = particle->product(it);
    if(cds->chi()>30) continue;
    if(vdis>9998||vdis>cds->vdis()){
      vdis = cds->vdis();
      vertex = cds->vbeam();
      if(GeomTools::GetID(vertex)==CID_Fiducial) fcds=it;
    }
  }
  if(fcds==-1) return true;
  FillHist("EventNumber",3); /* Vertex decision */
  pCDS* cds=0;
  if(fcds!=-1) cds = particle->product(fcds);

  bool lam=false; bool k0=false;
  if(cds->comb()==pow(2,CDS_PiPlus)+pow(2,CDS_PiMinus)) k0=true;
  if(cds->comb()==pow(2,CDS_Proton)+pow(2,CDS_PiMinus)) lam=true;

  if(lam) DoLambdaAna(particle,beam,cds,vertex);
  if(k0)  DoK0Ana(particle,beam,cds,vertex);

  return true;

}

bool MyAnalysisHeK2CDS::DoLambdaAna(Particle* particle, pBeam* beam, pCDS* cds, TVector3 vertex)
{

  TString frame[2] = {"Lab","CM"};

  // Target
  TVector3 ZeroV;
  TLorentzVector Ltgt; Ltgt.SetVectM(ZeroV,ThreeHeMass);
  TLorentzVector Lbeam = beam->GetLorentzVector(vertex);
  TVector3 boost = (Ltgt+Lbeam).BoostVector();
  TLorentzVector cLtgt = Ltgt; cLtgt.Boost(-boost); TLorentzVector cLbeam = Lbeam; cLbeam.Boost(-boost);
  /* CDS */
  TLorentzVector Lcds = cds->GetLorentzVector();
  TLorentzVector cLcds = cds->GetLorentzVector(); cLcds.Boost(-boost);
  /* CDS + NC */
  double cds_mass[2]  = { (Lcds).M()                           , (cLcds).M()};
  double cds_mom[2]   = { (Lcds).Vect().Mag()                  , (cLcds).Vect().Mag()};
  double cds_cost[2]  = { (Lcds).Vect().CosTheta()             , (cLcds).Vect().CosTheta()};
  double cds_phi[2]   = { (Lcds).Vect().Phi()                  , (cLcds).Vect().Phi()};
  double cds_mmass[2] = { (Ltgt+Lbeam-(Lcds)).M()              , (cLtgt+cLbeam-(cLcds)).M()};
  double cds_mmom[2]  = { (Ltgt+Lbeam-(Lcds)).Vect().Mag()     , (cLtgt+cLbeam-(cLcds)).Vect().Mag()};
  double cds_mcost[2] = { (Ltgt+Lbeam-(Lcds)).Vect().CosTheta(), (cLtgt+cLbeam-(cLcds)).Vect().CosTheta()};

  // Lambda
  for(int f=1; f<2; f++){
    FillHist(Form("lam_Mass_%s",frame[f].Data()),cds_mass[f]);
    FillHist(Form("lam_Momentum_%s",frame[f].Data()),cds_mom[f]);
    FillHist(Form("lam_CosT_%s",frame[f].Data()),cds_cost[f]);
    FillHist(Form("lam_Phi_%s",frame[f].Data()),cds_phi[f]);
    FillHist(Form("lam_MMass_%s",frame[f].Data()),cds_mmass[f]);
    FillHist(Form("lam_MMomentum_%s",frame[f].Data()),cds_mmom[f]);
    FillHist(Form("lam_MCosT_%s",frame[f].Data()),cds_mcost[f]);
  }

  // Z slice
  FillHist(Form("lam_Mass_vs_Vertex_Z"),vertex.Z(),cds_mass[1]);
  // R slice
  double r = TMath::Sqrt(vertex.X()*vertex.X()+vertex.Y()*vertex.Y());
  FillHist(Form("lam_Mass_vs_Vertex_R"),r,cds_mass[1]);
  // Phi slice
  double phi = vertex.Y()/TMath::Abs(vertex.Y())*TMath::ACos(vertex.X()/r);
  FillHist(Form("lam_Mass_vs_Vertex_Phi"),phi,cds_mass[1]);

  return true;

}

bool MyAnalysisHeK2CDS::DoK0Ana(Particle* particle, pBeam* beam, pCDS* cds, TVector3 vertex)
{

  TString frame[2] = {"Lab","CM"};

  // Target
  TVector3 ZeroV;
  TLorentzVector Ltgt; Ltgt.SetVectM(ZeroV,ThreeHeMass);
  TLorentzVector Lbeam = beam->GetLorentzVector(vertex);
  TVector3 boost = (Ltgt+Lbeam).BoostVector();
  TLorentzVector cLtgt = Ltgt; cLtgt.Boost(-boost); TLorentzVector cLbeam = Lbeam; cLbeam.Boost(-boost);
  /* CDS */
  TLorentzVector Lcds = cds->GetLorentzVector();
  TLorentzVector cLcds = cds->GetLorentzVector(); cLcds.Boost(-boost);
  /* CDS + NC */
  double cds_mass[2]  = { (Lcds).M()                           , (cLcds).M()};
  double cds_mom[2]   = { (Lcds).Vect().Mag()                  , (cLcds).Vect().Mag()};
  double cds_cost[2]  = { (Lcds).Vect().CosTheta()             , (cLcds).Vect().CosTheta()};
  double cds_phi[2]   = { (Lcds).Vect().Phi()                  , (cLcds).Vect().Phi()};
  double cds_mmass[2] = { (Ltgt+Lbeam-(Lcds)).M()              , (cLtgt+cLbeam-(cLcds)).M()};
  double cds_mmom[2]  = { (Ltgt+Lbeam-(Lcds)).Vect().Mag()     , (cLtgt+cLbeam-(cLcds)).Vect().Mag()};
  double cds_mcost[2] = { (Ltgt+Lbeam-(Lcds)).Vect().CosTheta(), (cLtgt+cLbeam-(cLcds)).Vect().CosTheta()};

  // Lambda
  for(int f=1; f<2; f++){
    FillHist(Form("k0_Mass_%s",frame[f].Data()),cds_mass[f]);
    FillHist(Form("k0_Momentum_%s",frame[f].Data()),cds_mom[f]);
    FillHist(Form("k0_CosT_%s",frame[f].Data()),cds_cost[f]);
    FillHist(Form("k0_Phi_%s",frame[f].Data()),cds_phi[f]);
    FillHist(Form("k0_MMass_%s",frame[f].Data()),cds_mmass[f]);
    FillHist(Form("k0_MMomentum_%s",frame[f].Data()),cds_mmom[f]);
    FillHist(Form("k0_MCosT_%s",frame[f].Data()),cds_mcost[f]);
  }

  // Z slice
  FillHist(Form("k0_Mass_vs_Vertex_Z"),vertex.Z(),cds_mass[1]);
  // R slice
  double r = TMath::Sqrt(vertex.X()*vertex.X()+vertex.Y()*vertex.Y());
  FillHist(Form("k0_Mass_vs_Vertex_R"),r,cds_mass[1]);
  // Phi slice
  double phi = vertex.Y()/TMath::Abs(vertex.Y())*TMath::ACos(vertex.X()/r);
  FillHist(Form("k0_Mass_vs_Vertex_Phi"),phi,cds_mass[1]);

  return true;

}


bool MyAnalysisHeK2CDS::FillHist(TString name, double val1, int weight)
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

bool MyAnalysisHeK2CDS::FillHist(TString name, TString val1, int weight)
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

bool MyAnalysisHeK2CDS::FillHist(TString name, double val1, double val2, int weight)
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

bool MyAnalysisHeK2CDS::FillHist(TString name, TString val1, TString val2, int weight)
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

void MyAnalysisHeK2CDS::CutCondition()
{
  /* K0 mass cut */ 
  double mean  = 0.4976;
  double sigma = 0.0069;
  k0ll = mean-2*sigma; k0ul = mean+2*sigma;
  sblk0ll = mean-5*sigma; sblk0ul = mean-3*sigma;
  sbuk0ll = mean+3*sigma; sbuk0ul = mean+5*sigma;
  /* Missing neutron mass cut */ 
  mean  = 0.9400;
  sigma = 0.0450;
  mnll = 0.85; mnul = 1.03;
  //mnll = mean-2*sigma; mnul = mean+2*sigma;
  sblmnll = mean-5*sigma; sblmnul = mean-3*sigma;
  sbumnll = mean+3*sigma; sbumnul = mean+5*sigma;
  /* SigmaPlus mass cut */
  mean  = 1.1888;
  sigma = 0.0044;
  spll = mean-2*sigma; spul = mean+2*sigma;
  sblspll = mean-5*sigma; sblspul = mean-3*sigma;
  sbuspll = mean+3*sigma; sbuspul = mean+5*sigma;
  /* SigmaMinus mass cut */
  mean  = 1.1970;
  sigma = 0.0054;
  smll = mean-2*sigma; smul = mean+2*sigma;
  sblsmll = mean-5*sigma; sblsmul = mean-3*sigma;
  sbusmll = mean+3*sigma; sbusmul = mean+5*sigma;
  /* Lambda mass cut */ 
  mean  = 1.1155;
  sigma = 0.0020;
  lamll = mean-2*sigma; lamul = mean+2*sigma;
  sbllamll = mean-5*sigma; sbllamul = mean-3*sigma;
  sbulamll = mean+3*sigma; sbulamul = mean+5*sigma;
}

bool MyAnalysisHeK2CDS::Initialize(ConfMan* confMan)
{
  std::cout << "### MyAnalysisHeK2CDS::Initialize ###" << std::endl;

  std::string ofname = confMan->GetOutFileName();
  ofname.insert(ofname.find(".root"),"_anaHeK2CDS");

  rtFile =  new TFile( Form("%s",ofname.c_str()), "RECREATE");
  rtFile -> cd();

  new TH1F( "EventNumber", "Number of Events", 20, 0, 20);

  std::cout << "Define Histograms for (K-,pip)X M.M." << std::endl;
  new TH1F("lam_Mass_CM","Invariant mass of #pi^{-}p;IM(#pi^{-}p) (GeV/c^{2});Counts", 5000, 1.0, 2.0);
  new TH1F("lam_Momentum_CM","Momentum of #pi^{-}p;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("lam_CosT_CM","cos#theta of #pi^{-}p;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("lam_Phi_CM","#phi of #pi^{-}p;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("lam_MMass_CM", "^{3}He(K^{-},#pi^{-}p)X missing mass;MM(#pi^{-}p) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("lam_MMomentum_CM", "^{3}He(K^{-},#pi^{-}p)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("lam_MCosT_CM", "^{3}He(K^{-},#pi^{-}p)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  /* 2D Plot */
  new TH2F("lam_Mass_vs_Vertex_Z","Invariant mass of #pi^{-}p vs. Vertex Z;z (cm);IM()#pi^{-}p (GeV/c^{2})",1000, -10, 10, 1000, 1.0, 1.2);
  new TH2F("lam_Mass_vs_Vertex_R","Invariant mass of #pi^{-}p vs. Vertex R;r (cm);IM()#pi^{-}p (GeV/c^{2})",1000, 0, 5, 1000, 1.0, 1.2);
  new TH2F("lam_Mass_vs_Vertex_Phi","Invariant mass of #pi^{-}p vs. Vertex #phi;#phi (rad.);IM()#pi^{-}p (GeV/c^{2})",1000, -1.6, 1.6, 1000, 1.0, 1.2);

  std::cout << "=== End of [Lambda::Initialize] === " << std::endl;

  std::cout << "Define Histograms for (K-,pipi)X M.M." << std::endl;
  new TH1F("k0_Mass_CM","Invariant mass of #pi^{-}#pi^{+};IM(#pi^{-}#pi^{+}) (GeV/c^{2});Counts", 5000, 0.0, 1.0);
  new TH1F("k0_Momentum_CM","Momentum of #pi^{-}#pi^{+};Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("k0_CosT_CM","cos#theta of #pi^{-}#pi^{+};cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("k0_Phi_CM","#phi of #pi^{-}#pi^{+};#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("k0_MMass_CM", "^{3}He(K^{-},#pi^{-}#pi^{+})X missing mass;MM(#pi^{-}#pi^{+}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("k0_MMomentum_CM", "^{3}He(K^{-},#pi^{-}#pi^{+})X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("k0_MCosT_CM", "^{3}He(K^{-},#pi^{-}#pi^{+})X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  /* 2D Plot */
  new TH2F("k0_Mass_vs_Vertex_Z","Invariant mass of #pi^{-}#pi^{+} vs. Vertex Z;z (cm);IM(#pi^{-}#pi^{+}) (GeV/c^{2})",1000, -10, 10, 1000, 0.4, 0.6);
  new TH2F("k0_Mass_vs_Vertex_R","Invariant mass of #pi^{-}#pi^{+} vs. Vertex R;r (cm);IM(#pi^{-}#pi^{+}) (GeV/c^{2})",1000, 0, 5, 1000, 0.4, 0.6);
  new TH2F("k0_Mass_vs_Vertex_Phi","Invariant mass of #pi^{-}#pi^{+} vs. Vertex #phi;#phi (rad.);IM(#pi^{-}#pi^{+}) (GeV/c^{2})",1000, -1.6, 1.6, 1000, 0.4, 0.6);

  std::cout << "=== End of [K0::Initialize] === " << std::endl;

  return true;
}
