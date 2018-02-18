// MyAnalysisFWDNeutral.cpp

#include "MyAnalysisFWDNeutral.h"

MyAnalysisFWDNeutral::MyAnalysisFWDNeutral(TFile* rtFile, ConfMan* conf)
{
  Initialize(rtFile, conf);
  Clear();
}

MyAnalysisFWDNeutral::~MyAnalysisFWDNeutral()
{
}

void MyAnalysisFWDNeutral::Clear()
{
}

bool MyAnalysisFWDNeutral::DoAnalysis(ConfMan* conf, EventHeader* header, BeamLineHitMan* blMan, BeamLineTrackMan* bltrackMan, CDSHitMan* cdsMan, CDSTrackingMan* cdstrackMan, Particle* particle)
{
  DetectorList *dlist=DetectorList::GetInstance();
  if(particle->nBeam()!=1) return false;
  pBeam* beam = particle->beam(0);

  // ################# //
  // Neutral selection //
  // ################# //
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
  if(MulBVC!=0||MulCVC!=0) return true;

  // ########### //
  // NC analysis //
  // ########### //
  for(int i=0; i<blMan->nNC(); i++){
    HodoscopeLikeHit* hit = blMan->NC(i);
    if(hit->CheckRange()){
      pNC* nc = new pNC();
      nc->SetSegment(hit->seg());
      nc->SetHitPosition(hit->pos().X(),hit->hitpos(),hit->pos().Z());
      nc->SetTime(hit->ctmean());
      nc->SetEnergy(hit->emean());
      particle->AddNC(*nc);
      delete nc;
    }
  }

  // ############### //
  // Vertex decision //
  // ############### //
  TVector3 vertex;
  double vdis = 9999;
  int fcds = -1;
  for(int it=0; it<particle->nCDS(); it++){
    pCDS* cdstrack = particle->cdsi(it);
    if(vdis>9998||vdis>cdstrack->vdis()){
      vdis = cdstrack->vdis();
      vertex = cdstrack->vbeam();
      if(GeomTools::GetID(vertex)==CID_Fiducial) fcds=it;
    }
  }
  double time = 9999;
  int fnc = -1;
  if(fcds>=0){
    for(int it=0; it<particle->nNC(); it++){
      pNC* nc = particle->nc(it);
      if((time>9998||time>nc->time())){
        fnc = it;
        time = nc->time();
        nc->CalcMom(beam,vertex);
      }
    }
  }
  if(fcds>=0&&fnc>=0){
    pCDS* cds = particle->cdsi(fcds);
    pNC*   nc = particle->nc(fnc);
    FillHist("TOFvsMomentum_FWD",nc->tof(),nc->mom().Mag());
    FillHist("OverbetavsMomentum_FWD",1./nc->beta(),nc->mom().Mag());
    FillHist("OverbetavsEnergy_FWD",1./nc->beta(),nc->energy());
    FillHist("NCHitPosition_XY",nc->hitpos().X(),nc->hitpos().Y());
    FillHist("NCHitSegment",nc->seg()%17,nc->seg()/17+1);
    if(nc->pid()==F_Neutron){
      TLorentzVector LB = beam->GetLorentzVector(vertex);
      TVector3 ZeroV;
      TLorentzVector TL; TL.SetVectM(ZeroV,pMass);
      TLorentzVector NfL = nc->GetLorentzVector();
      TLorentzVector CDSL = cds->GetLorentzVector();
      // Missing mass p(K-,n)X
      TLorentzVector miss = (LB+TL) - NfL;
      FillHist("MM_Nf",miss.M());
      // Invariant mass
      if(cds->pid()==CDS_PiPlus) FillHist("IM_Nfpiplus",(NfL+CDSL).M());
      if(cds->pid()==CDS_PiMinus) FillHist("IM_Nfpiminus",(NfL+CDSL).M());
      if(cds->pid()==CDS_Kaon) FillHist("IM_NfK",(NfL+CDSL).M());
    }
  }

  return true;

}

bool MyAnalysisFWDNeutral::FillHist(TString name, double val1)
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

bool MyAnalysisFWDNeutral::FillHist(TString name, double val1, double val2)
{
  TH2F* h2 = (TH2F*)gFile -> Get(name);
  if(h2){ h2 -> Fill(val1,val2);
    return true;
  }
  else {
    return false;
  }
}

bool MyAnalysisFWDNeutral::Initialize(TFile* rtFile, ConfMan* confMan)
{
  rtFile -> cd();

  // FWD neutral Particles
  std::cout << "Define Histograms for FWD neutral particles" << std::endl;
  new TH2F( "OverbetavsMomentum_FWD", "1/#beta vs. Momentum;1/#beta;Momentum (GeV/c)", 1000, -0, 10, 500, -5, 5 );
  new TH2F( "OverbetavsEnergy_FWD", "1/#beta vs. Energy deposit;1/#beta;Energy deposit (MeVee)", 1000, -0, 10, 5000, 0, 100);
  new TH2F( "TOFvsMomentum_FWD", "TOF vs. Momentum;TOF (ns);Momentum (GeV/c)", 4000, 30, 70, 500, -5, 5 );
  new TH2F( "NCHitPosition_XY", "Hit position at NC (X vs. Y);X position (cm);Y position (cm)",16,-160,160,15,-75,75);
  new TH2F( "NCHitSegment", "NC segment;Segment in layer;Layer",16,0.5,16.5,7,0.5,7.5);
  // Missing mass p(K-,n)X
  std::cout << "Define Histograms for p(K-,n)X M.M." << std::endl;
  new TH1F( "MM_Nf", "p(K^{-},n)X missing mass;MM_{p}(n_{f}) (GeV/c^{2});Counts", 2000, 0.0, 2.0 );
  // Invariant mass with FWD
  std::cout << "Define Histograms for product I.M. with FWD" << std::endl;
  new TH1F( "IM_Nfpiplus", "Invariant mass (n_{f}#pi^{+});IM(n_{f}#pi^{+}) (GeV/c^{2});Counts", 2000, 0.0, 2.0 );
  new TH1F( "IM_Nfpiminus", "Invariant mass (n_{f}#pi^{-});IM(n_{f}#pi^{-}) (GeV/c^{2});Counts", 2000, 0.0, 2.0 );
  new TH2F("pKn_MMvsMomentum_CM", "p(K^{-},n)X missing mass vs. missing momentum;MM_{p}(n) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKn_MMvsCosT_CM", "p(K^{-},n)X missing mass vs. cos#theta;MM_{p}(n) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
}
