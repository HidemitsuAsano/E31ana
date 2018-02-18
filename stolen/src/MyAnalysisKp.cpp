// MyAnalysisKp.cpp

#include "MyAnalysisKp.h"

MyAnalysisKp::MyAnalysisKp(TFile* rtFile, ConfMan* conf)
{
  Initialize(rtFile, conf);
  Clear();
}

MyAnalysisKp::~MyAnalysisKp()
{
}

void MyAnalysisKp::Clear()
{
}

bool MyAnalysisKp::DoAnalysis(ConfMan* conf, EventHeader* header, BeamLineHitMan* blMan, BeamLineTrackMan* bltrackMan, CDSHitMan* cdsMan, CDSTrackingMan* cdstrackMan, Particle* particle)
{
  DetectorList *dlist=DetectorList::GetInstance();
  if(particle->nBeam()!=1) return false;
  pBeam* beam = particle->beam(0);

  // ######################### //
  // Checking for CDS particle //
  // ######################### //
  int ncds[6] = {0,0,0,0,0,0};
  ncds[0] = particle->nPiplus();
  ncds[1] = particle->nPiminus();
  ncds[2] = particle->nKaon();
  ncds[3] = particle->nProton();
  ncds[4] = particle->nDeuteron();
  ncds[5] = particle->nTriton();
  ncds[5] += particle->nHelium3();
  ncds[5] += particle->nOther();
  TString cdsstr[6] = {"#pi^{+}","#pi^{-}","K^{-}","p","d","other"};
  for(int i=0; i<6; i++){
    FillHist("Particle_CDS",cdsstr[i],ncds[i]);
  }

  // ############### //
  // Vertex decision //
  // ############### //
  //TVector3 vertex;
  //double vdis = 9999;
  //int fcds = -1;
  //for(int it=0; it<particle->nCDS(); it++){
  //  pCDS* cdstrack = particle->cdsi(it);
  //  if(vdis>9998||vdis>cdstrack->vdis()){
  //    vdis = cdstrack->vdis();
  //    vertex = cdstrack->vbeam();
  //    if(GeomTools::GetID(vertex)==CID_Fiducial) fcds=it;
  //  }
  //}
  //if(fcds==-1) return false; /* there is no vertex in fiducial volume */

  //// ############### //
  //// NC hit decision //
  //// ############### //
  //double time = 9999;
  //int fnc = -1;
  //for(int it=0; it<particle->nNC(); it++){
  //  pNC* nc = particle->nc(it);
  //  if(nc->energy()>8.0 && (time>9998||time>nc->time())){
  //    time = nc->time();
  //    nc->CalcMom(beam,vertex);
  //    if(nc->pid()==F_Neutron) fnc = it;
  //  }
  //}
  //if(fnc==-1) return false; /* there is no NC hit */

  //pCDS* pip = particle->pip(0);
  //pCDS* pim = particle->pim(0);
  //pNC*  nf  = particle->nc(fnc);

  //TLorentzVector Lpip  =  pip->GetLorentzVector();
  //TLorentzVector Lpim  =  pim->GetLorentzVector();
  //TLorentzVector Lnf   =   nf->GetLorentzVector();
  //TLorentzVector Lbeam = beam->GetLorentzVector(vertex);
  //TVector3 ZeroV;
  //TLorentzVector Ltgt; Ltgt.SetVectM(ZeroV,pMass);

  //double im_pipi = (Lpip+Lpim).M();
  //double im_pipn = (Lpip+Lnf).M();
  //double im_pimn = (Lpim+Lnf).M();
  //double mm_pipin = (Lbeam+Ltgt-(Lpip+Lpim+Lnf)).M();

  //FillHist("TOFvsMomentum_FWD",nf->tof(),nf->mom().Mag());
  //FillHist("OverbetavsMomentum_FWD",1./nf->beta(),nf->mom().Mag());
  //FillHist("OverbetavsEnergy_FWD",1./nf->beta(),nf->energy());

  //FillHist("IM_pipi",im_pipi);
  //FillHist("IM_pipn",im_pipn);
  //FillHist("IM_pimn",im_pimn);
  //FillHist("MM_pipin",mm_pipin);

  /* Forward sigma reduction */

  return true;

}

bool MyAnalysisKp::FillHist(TString name, double val1, int weight)
{
  TH1F* h1 = (TH1F*)gFile -> Get(name);
  if(h1){
    h1 -> Fill(val1,weight);
    return true;
  }
  else {
    return false;
  }
}

bool MyAnalysisKp::FillHist(TString name, TString val1, int weight)
{
  TH1F* h1 = (TH1F*)gFile -> Get(name);
  if(h1){
    h1 -> Fill(val1,weight);
    return true;
  }
  else {
    return false;
  }
}

bool MyAnalysisKp::FillHist(TString name, double val1, double val2, int weight)
{
  TH2F* h2 = (TH2F*)gFile -> Get(name);
  if(h2){ h2 -> Fill(val1,val2,weight);
    return true;
  }
  else {
    return false;
  }
}

bool MyAnalysisKp::FillHist(TString name, TString val1, TString val2, int weight)
{
  TH2F* h2 = (TH2F*)gFile -> Get(name);
  if(h2){ h2 -> Fill(val1,val2, weight);
    return true;
  }
  else {
    return false;
  }
}

bool MyAnalysisKp::Initialize(TFile* rtFile, ConfMan* confMan)
{
  rtFile -> cd();

  // CDS Particles
  std::cout << "Define Histograms for CDS particles" << std::endl;
  TH1F* h1;
  h1 = new TH1F("Particle_CDS","Detected particles by CDS",6,0,6);
  h1->GetXaxis()->SetBinLabel(1,"#pi^{+}");
  h1->GetXaxis()->SetBinLabel(2,"#pi^{-}");
  h1->GetXaxis()->SetBinLabel(3,"K^{-}");
  h1->GetXaxis()->SetBinLabel(4,"p");
  h1->GetXaxis()->SetBinLabel(5,"d");
  h1->GetXaxis()->SetBinLabel(6,"other");
  // FWD neutral Particles
  std::cout << "Define Histograms for FWD neutral particles" << std::endl;
  new TH2F( "OverbetavsMomentum_FWD", "1/#beta vs. Momentum;1/#beta;Momentum (GeV/c)", 1000, -0, 10, 500, -5, 5 );
  new TH2F( "OverbetavsEnergy_FWD", "1/#beta vs. Energy deposit;1/#beta;Energy deposit (MeVee)", 1000, -0, 10, 5000, 0, 100);
  new TH2F( "TOFvsMomentum_FWD", "TOF vs. Momentum;TOF (ns);Momentum (GeV/c)", 4000, 30, 70, 500, -5, 5 );
  new TH2F( "NCHitPosition_XY", "Hit position at NC (X vs. Y);X position (cm);Y position (cm)",16,-160,160,15,-75,75);
  new TH2F( "NCHitSegment", "NC segment;Segment in layer;Layer",16,0.5,16.5,7,0.5,7.5);
  // Missing mass d(K-,n)X
  std::cout << "Define Histograms for p(K-,n)X M.M." << std::endl;
  new TH1F( "MM_Nf", "d(K^{-},n)X missing mass;MM_{p}(n_{f}) (GeV/c^{2});Counts", 2000, 0.0, 2.0 );
  // Invariant mass with FWD
  std::cout << "Define Histograms for product I.M. with FWD" << std::endl;
  new TH1F( "IM_pipi", "Invariant mass (#pi^{+}#pi^{-});IM(#pi^{+}#pi^{-}) (GeV/c^{2});Counts", 1000, 0.0, 1.0 );
  new TH1F( "IM_pipn", "Invariant mass (#pi^{+}n_{f});IM(#pi^{+}n_{f}) (GeV/c^{2});Counts", 1000, 1.0, 2.0 );
  new TH1F( "IM_pimn", "Invariant mass (#pi^{-}n_{f});IM(#pi^{-}n_{f}) (GeV/c^{2});Counts", 1000, 1.0, 2.0 );
  new TH1F( "MM_pipin", "d(K^{-},#pi^{+}#pi^{-}n_{f})X missing mass;MM_{d}($pi^{+}#pi^{-}n_{f}) (GeV/c^{2});Counts", 1000, 0.5, 1.5 );
}
