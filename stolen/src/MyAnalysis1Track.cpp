// MyAnalysis1Track.cpp

#include "MyAnalysis1Track.h"

MyAnalysis1Track::MyAnalysis1Track(TFile* rtFile, ConfMan* conf)
{
  Initialize(rtFile, conf);
  Clear();
}

MyAnalysis1Track::~MyAnalysis1Track()
{
}

void MyAnalysis1Track::Clear()
{
}

bool MyAnalysis1Track::DoAnalysis(ConfMan* conf, EventHeader* header, BeamLineHitMan* blMan, BeamLineTrackMan* bltrackMan, CDSHitMan* cdsMan, CDSTrackingMan* cdstrackMan, Particle* particle)
{
  DetectorList *dlist=DetectorList::GetInstance();
  if(particle->nBeam()!=1) return false;
  pBeam* beam = particle->beam(0);

  /* Only CDS 1 track analysis */
  //if(particle->nCDS()!=1) return false;

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
  if(fcds==-1) return false; /* there is no vertex in fiducial volume */

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
  TString cdsstr[6] = {"#pi^{+}","#pi^{-}","K^{-}","p","d","Other"};
  FillHist("CDS_NumOfParticle",particle->nCDS());
  for(int i=0; i<6; i++){
    FillHist("CDS_Particle",cdsstr[i],Form("%d track",particle->nCDS()),ncds[i]);
  }

  // ############### //
  // NC hit decision //
  // ############### //
  double time = 9999;
  int fnc = -1;
  for(int it=0; it<particle->nNC(); it++){
    pNC* nc = particle->nc(it);
    if(nc->energy()>8.0 && (time>9998||time>nc->time())){
      time = nc->time();
      nc->CalcMom(beam,vertex);
      if(nc->pid()==F_Neutron) fnc = it;
    }
  }
  //if(fnc==-1) return false; /* there is no NC hit */

  pCDS* cds=0;
  pNC*  nc =0;
  if(fcds!=-1) cds=particle->cdsi(fcds);
  if(fnc!=-1)  nc =particle->nc(fnc);

  /* CDS particle analysis */
  if(cds==0) return false;
  TVector3 ZeroV;
  TLorentzVector Ltgt; Ltgt.SetVectM(ZeroV,pMass);
  TLorentzVector Lbeam = beam->GetLorentzVector(vertex);
  TVector3 boost = (Ltgt+Lbeam).BoostVector();
  TLorentzVector cLtgt = Ltgt; cLtgt.Boost(-boost); TLorentzVector cLbeam = Lbeam; cLbeam.Boost(-boost);
  TLorentzVector Lcds = cds->GetLorentzVector();
  TLorentzVector cLcds = cds->GetLorentzVector(); cLcds.Boost(-boost);
  double cds_mom[2]   = { Lcds.Vect().Mag()                  , cLcds.Vect().Mag()};
  double cds_cost[2]  = { Lcds.Vect().CosTheta()             , cLcds.Vect().CosTheta()};
  double cds_phi[2]   = { Lcds.Vect().Phi()                  , cLcds.Vect().Phi()};
  double cds_mmass[2] = { (Ltgt+Lbeam-Lcds).M()              , (cLtgt+cLbeam-cLcds).M()};
  double cds_mmom[2]  = { (Ltgt+Lbeam-Lcds).Vect().Mag()     , (cLtgt+cLbeam-cLcds).Vect().Mag()};
  double cds_mcost[2] = { (Ltgt+Lbeam-Lcds).Vect().CosTheta(), (cLtgt+cLbeam-cLcds).Vect().CosTheta()};

  TString pid[8] = {"pip","p","d","t","he","pim","k","o"};
  TString frame[2] = {"Lab","CM"};
  for(int f=0; f<2; f++){
    FillHist(Form("CDS_CosTvsMomentum_%s",frame[f].Data()),cds_cost[f],cds_mom[f]);
    FillHist(Form("CDS_CosTvsPhi_%s",frame[f].Data()),cds_cost[f],cds_phi[f]);
    FillHist(Form("CDS%s_CosTvsMomentum_%s",pid[cds->pid()].Data(),frame[f].Data()),cds_cost[f],cds_mom[f]);
    FillHist(Form("CDS%s_CosTvsPhi_%s",pid[cds->pid()].Data(),frame[f].Data()),cds_cost[f],cds_phi[f]);
    FillHist(Form("pK%s_MMvsMomentum_%s",pid[cds->pid()].Data(),frame[f].Data()),cds_mmass[f],cds_mmom[f]);
    FillHist(Form("pK%s_MMvsCosT_%s",pid[cds->pid()].Data(),frame[f].Data()),cds_mmass[f],cds_mcost[f]);
  }

  /* with forward neutron */
  if(nc==0) return false;
  TLorentzVector Lnc  = nc ->GetLorentzVector();
  TLorentzVector cLnc  = nc ->GetLorentzVector(); cLnc.Boost(-boost);
  double nc_mom[2]   = { Lnc.Vect().Mag()                  , cLnc.Vect().Mag()};
  double nc_cost[2]  = { Lnc.Vect().CosTheta()             , cLnc.Vect().CosTheta()};
  double nc_phi[2]   = { Lnc.Vect().Phi()                  , cLnc.Vect().Phi()};
  double nc_mmass[2] = { (Ltgt+Lbeam-Lnc).M()              , (cLtgt+cLbeam-cLnc).M()};
  double nc_mmom[2]  = { (Ltgt+Lbeam-Lnc).Vect().Mag()     , (cLtgt+cLbeam-cLnc).Vect().Mag()};
  double nc_mcost[2] = { (Ltgt+Lbeam-Lnc).Vect().CosTheta(), (cLtgt+cLbeam-cLnc).Vect().CosTheta()};
  for(int f=0; f<2; f++){
    FillHist(Form("FWDn_CosTvsMomentum_%s",frame[f].Data()),nc_cost[f],nc_mom[f]);
    FillHist(Form("FWDn_CosTvsPhi_%s",frame[f].Data()),nc_cost[f],nc_phi[f]);
    FillHist(Form("pKn_MMvsMomentum_%s",frame[f].Data()),nc_mmass[f],nc_mmom[f]);
    FillHist(Form("pKn_MMvsCosT_%s",frame[f].Data()),nc_mmass[f],nc_mcost[f]);
  }

  /* Invariant mass CDS+NC */
  double cdsnc_imass[2] = { (Lcds+Lnc).M(), (cLcds+cLnc).M() };
  double cdsnc_imom[2] = { (Lcds+Lnc).Vect().Mag(), (cLcds+cLnc).Vect().Mag() };
  double cdsnc_icost[2] = { (Lcds+Lnc).Vect().CosTheta(), (cLcds+cLnc).Vect().CosTheta() };
  double cdsnc_mmass[2] = { (Ltgt+Lbeam-Lcds-Lnc).M(), (cLtgt+cLbeam-cLcds-cLnc).M() };
  double cdsnc_mmom[2] = { (Ltgt+Lbeam-Lcds-Lnc).Vect().Mag(), (cLtgt+cLbeam-cLcds-cLnc).Vect().Mag() };
  double cdsnc_mcost[2] = { (Ltgt+Lbeam-Lcds-Lnc).Vect().CosTheta(), (cLtgt+cLbeam-cLcds-cLnc).Vect().CosTheta() };

  for(int f=0; f<2; f++){
    FillHist(Form("n%s_IMvsMomentum_%s",pid[cds->pid()].Data(),frame[f].Data()),cdsnc_imass[f],cdsnc_imom[f]);
    FillHist(Form("n%s_IMvsCosT_%s",pid[cds->pid()].Data(),frame[f].Data()),cdsnc_imass[f],cdsnc_icost[f]);
    FillHist(Form("pKn%s_MMvsMomentum_%s",pid[cds->pid()].Data(),frame[f].Data()),cdsnc_mmass[f],cdsnc_mmom[f]);
    FillHist(Form("pKn%s_MMvsCosT_%s",pid[cds->pid()].Data(),frame[f].Data()),cdsnc_mmass[f],cdsnc_mcost[f]);
  }

  return true;

}

bool MyAnalysis1Track::FillHist(TString name, double val1, int weight)
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

bool MyAnalysis1Track::FillHist(TString name, TString val1, int weight)
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

bool MyAnalysis1Track::FillHist(TString name, double val1, double val2, int weight)
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

bool MyAnalysis1Track::FillHist(TString name, TString val1, TString val2, int weight)
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

bool MyAnalysis1Track::Initialize(TFile* rtFile, ConfMan* confMan)
{
  rtFile -> cd();

  // CDS Particles
  std::cout << "Define Histograms for CDS particles" << std::endl;
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
  new TH2F("CDS_CosTvsMomentum_Lab","cos#theta vs. momentum of CDS particle;cos#theta;Momentum (GeV/c)", 2000, -1.0, 1.0, 1000, 0.0, 1.0);
  new TH2F("CDSpip_CosTvsMomentum_Lab","cos#theta vs. momentum of CDS particle;cos#theta;Momentum (GeV/c)", 2000, -1.0, 1.0, 1000, 0.0, 1.0);
  new TH2F("CDSpim_CosTvsMomentum_Lab","cos#theta vs. momentum of CDS particle;cos#theta;Momentum (GeV/c)", 2000, -1.0, 1.0, 1000, 0.0, 1.0);
  new TH2F("CDSk_CosTvsMomentum_Lab","cos#theta vs. momentum of CDS particle;cos#theta;Momentum (GeV/c)", 2000, -1.0, 1.0, 1000, 0.0, 1.0);
  new TH2F("CDSp_CosTvsMomentum_Lab","cos#theta vs. momentum of CDS particle;cos#theta;Momentum (GeV/c)", 2000, -1.0, 1.0, 1000, 0.0, 1.0);
  new TH2F("CDSd_CosTvsMomentum_Lab","cos#theta vs. momentum of CDS particle;cos#theta;Momentum (GeV/c)", 2000, -1.0, 1.0, 1000, 0.0, 1.0);
  new TH2F("CDS_CosTvsMomentum_CM","cos#theta vs. momentum of CDS particle;cos#theta;Momentum (GeV/c)", 2000, -1.0, 1.0, 1000, 0.0, 1.0);
  new TH2F("CDSpip_CosTvsMomentum_CM","cos#theta vs. momentum of CDS particle;cos#theta;Momentum (GeV/c)", 2000, -1.0, 1.0, 1000, 0.0, 1.0);
  new TH2F("CDSpim_CosTvsMomentum_CM","cos#theta vs. momentum of CDS particle;cos#theta;Momentum (GeV/c)", 2000, -1.0, 1.0, 1000, 0.0, 1.0);
  new TH2F("CDSk_CosTvsMomentum_CM","cos#theta vs. momentum of CDS particle;cos#theta;Momentum (GeV/c)", 2000, -1.0, 1.0, 1000, 0.0, 1.0);
  new TH2F("CDSp_CosTvsMomentum_CM","cos#theta vs. momentum of CDS particle;cos#theta;Momentum (GeV/c)", 2000, -1.0, 1.0, 1000, 0.0, 1.0);
  new TH2F("CDSd_CosTvsMomentum_CM","cos#theta vs. momentum of CDS particle;cos#theta;Momentum (GeV/c)", 2000, -1.0, 1.0, 1000, 0.0, 1.0);
  new TH2F("CDS_CosTvsPhi_Lab","cos#theta vs. #phi of CDS particle;cos#theta;#phi", 2000, -1.0, 1.0, 1600, 0.0, 3.2);
  new TH2F("CDSpip_CosTvsPhi_Lab","cos#theta vs. #phi of CDS particle;cos#theta;#phi", 2000, -1.0, 1.0, 1600, 0.0, 3.2);
  new TH2F("CDSpim_CosTvsPhi_Lab","cos#theta vs. #phi of CDS particle;cos#theta;#phi", 2000, -1.0, 1.0, 1600, 0.0, 3.2);
  new TH2F("CDSk_CosTvsPhi_Lab","cos#theta vs. #phi of CDS particle;cos#theta;#phi", 2000, -1.0, 1.0, 1600, 0.0, 3.2);
  new TH2F("CDSp_CosTvsPhi_Lab","cos#theta vs. #phi of CDS particle;cos#theta;#phi", 2000, -1.0, 1.0, 1600, 0.0, 3.2);
  new TH2F("CDSd_CosTvsPhi_Lab","cos#theta vs. #phi of CDS particle;cos#theta;#phi", 2000, -1.0, 1.0, 1600, 0.0, 3.2);
  new TH2F("CDS_CosTvsPhi_CM","cos#theta vs. #phi of CDS particle;cos#theta;#phi", 2000, -1.0, 1.0, 1600, 0.0, 3.2);
  new TH2F("CDSpip_CosTvsPhi_CM","cos#theta vs. #phi of CDS particle;cos#theta;#phi", 2000, -1.0, 1.0, 1600, 0.0, 3.2);
  new TH2F("CDSpim_CosTvsPhi_CM","cos#theta vs. #phi of CDS particle;cos#theta;#phi", 2000, -1.0, 1.0, 1600, 0.0, 3.2);
  new TH2F("CDSk_CosTvsPhi_CM","cos#theta vs. #phi of CDS particle;cos#theta;#phi", 2000, -1.0, 1.0, 1600, 0.0, 3.2);
  new TH2F("CDSp_CosTvsPhi_CM","cos#theta vs. #phi of CDS particle;cos#theta;#phi", 2000, -1.0, 1.0, 1600, 0.0, 3.2);
  new TH2F("CDSd_CosTvsPhi_CM","cos#theta vs. #phi of CDS particle;cos#theta;#phi", 2000, -1.0, 1.0, 1600, 0.0, 3.2);
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
  new TH2F("pKpip_MMvsMomentum_Lab", "p(K^{-},#pi^{+})X missing mass vs. missing momentum;MM_{p}(#pi^{+}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKpim_MMvsMomentum_Lab", "p(K^{-},#pi^{-})X missing mass vs. missing momentum;MM_{p}(#pi^{-}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKk_MMvsMomentum_Lab", "p(K^{-},K^{-})X missing mass vs. missing momentum;MM_{p}(K^{-}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKp_MMvsMomentum_Lab", "p(K^{-},p)X missing mass vs. missing momentum;MM_{p}(p) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKd_MMvsMomentum_Lab", "p(K^{-},d)X missing mass vs. missing momentum;MM_{p}(d) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKpip_MMvsMomentum_CM", "p(K^{-},#pi^{+})X missing mass vs. missing momentum;MM_{p}(#pi^{+}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKpim_MMvsMomentum_CM", "p(K^{-},#pi^{-})X missing mass vs. missing momentum;MM_{p}(#pi^{-}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKk_MMvsMomentum_CM", "p(K^{-},K^{-})X missing mass vs. missing momentum;MM_{p}(K^{-}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKp_MMvsMomentum_CM", "p(K^{-},p)X missing mass vs. missing momentum;MM_{p}(p) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKd_MMvsMomentum_CM", "p(K^{-},d)X missing mass vs. missing momentum;MM_{p}(d) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKpip_MMvsCosT_Lab", "p(K^{-},#pi^{+})X missing mass vs. cos#theta;MM_{p}(#pi^{+}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKpim_MMvsCosT_Lab", "p(K^{-},#pi^{-})X missing mass vs. cos#theta;MM_{p}(#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKk_MMvsCosT_Lab", "p(K^{-},K^{-})X missing mass vs. cos#theta;MM_{p}(K^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKp_MMvsCosT_Lab", "p(K^{-},p)X missing mass vs. cos#theta;MM_{p}(p) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKd_MMvsCosT_Lab", "p(K^{-},d)X missing mass vs. cos#theta;MM_{p}(d) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKpip_MMvsCosT_CM", "p(K^{-},#pi^{+})X missing mass vs. cos#theta;MM_{p}(#pi^{+}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKpim_MMvsCosT_CM", "p(K^{-},#pi^{-})X missing mass vs. cos#theta;MM_{p}(#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKk_MMvsCosT_CM", "p(K^{-},K^{-})X missing mass vs. cos#theta;MM_{p}(K^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKp_MMvsCosT_CM", "p(K^{-},p)X missing mass vs. cos#theta;MM_{p}(p) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKd_MMvsCosT_CM", "p(K^{-},d)X missing mass vs. cos#theta;MM_{p}(d) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKnpip_MMvsMomentum_Lab", "p(K^{-},n#pi^{+})X missing mass vs. missing momentum;MM_{p}(n#pi^{+}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKnpim_MMvsMomentum_Lab", "p(K^{-},n#pi^{-})X missing mass vs. missing momentum;MM_{p}(n#pi^{-}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKnk_MMvsMomentum_Lab", "p(K^{-},nK^{-})X missing mass vs. missing momentum;MM_{p}(nK^{-}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKnp_MMvsMomentum_Lab", "p(K^{-},np)X missing mass vs. missing momentum;MM_{p}(np) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKnd_MMvsMomentum_Lab", "p(K^{-},nd)X missing mass vs. missing momentum;MM_{p}(nd) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKnpip_MMvsMomentum_CM", "p(K^{-},n#pi^{+})X missing mass vs. missing momentum;MM_{p}(n#pi^{+}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKnpim_MMvsMomentum_CM", "p(K^{-},n#pi^{-})X missing mass vs. missing momentum;MM_{p}(n#pi^{-}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKnk_MMvsMomentum_CM", "p(K^{-},nK^{-})X missing mass vs. missing momentum;MM_{p}(nK^{-}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKnp_MMvsMomentum_CM", "p(K^{-},np)X missing mass vs. missing momentum;MM_{p}(np) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKnd_MMvsMomentum_CM", "p(K^{-},nd)X missing mass vs. missing momentum;MM_{p}(nd) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKnpip_MMvsCosT_Lab", "p(K^{-},n#pi^{+})X missing mass vs. cos#theta;MM_{p}(n#pi^{+}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKnpim_MMvsCosT_Lab", "p(K^{-},n#pi^{-})X missing mass vs. cos#theta;MM_{p}(n#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKnk_MMvsCosT_Lab", "p(K^{-},nK^{-})X missing mass vs. cos#theta;MM_{p}(nK^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKnp_MMvsCosT_Lab", "p(K^{-},np)X missing mass vs. cos#theta;MM_{p}(np) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKnd_MMvsCosT_Lab", "p(K^{-},nd)X missing mass vs. cos#theta;MM_{p}(nd) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKnpip_MMvsCosT_CM", "p(K^{-},n#pi^{+})X missing mass vs. cos#theta;MM_{p}(n#pi^{+}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKnpim_MMvsCosT_CM", "p(K^{-},n#pi^{-})X missing mass vs. cos#theta;MM_{p}(n#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKnk_MMvsCosT_CM", "p(K^{-},nK^{-})X missing mass vs. cos#theta;MM_{np}(K^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKnp_MMvsCosT_CM", "p(K^{-},np)X missing mass vs. cos#theta;MM_{p}(np) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKnd_MMvsCosT_CM", "p(K^{-},nd)X missing mass vs. cos#theta;MM_{p}(nd) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  // Invariant mass with FWD
  std::cout << "Define Histograms for product I.M. with FWD" << std::endl;
  new TH2F( "npip_IMvsMomentum_Lab", "Invariant mass (n#pi^{+}) vs Momentum;IM(n#pi^{+}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0 );
  new TH2F( "npim_IMvsMomentum_Lab", "Invariant mass (n#pi^{-}) vs Momentum;IM(n#pi^{-}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0 );
  new TH2F( "nk_IMvsMomentum_Lab", "Invariant mass (nK^{-}) vs Momentum;IM(nK^{-}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0 );
  new TH2F( "np_IMvsMomentum_Lab", "Invariant mass (np) vs Momentum;IM(np) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0 );
  new TH2F( "nd_IMvsMomentum_Lab", "Invariant mass (nd) vs Momentum;IM(nd) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0 );
  new TH2F( "npip_IMvsMomentum_CM", "Invariant mass (n#pi^{+}) vs Momentum;IM(n#pi^{+}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0 );
  new TH2F( "npim_IMvsMomentum_CM", "Invariant mass (n#pi^{-}) vs Momentum;IM(n#pi^{-}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0 );
  new TH2F( "nk_IMvsMomentum_CM", "Invariant mass (nK^{-}) vs Momentum;IM(nK^{-}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0 );
  new TH2F( "np_IMvsMomentum_CM", "Invariant mass (np) vs Momentum;IM(np) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0 );
  new TH2F( "nd_IMvsMomentum_CM", "Invariant mass (nd) vs Momentum;IM(nd) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0 );
  new TH2F( "npip_IMvsCosT_Lab", "Invariant mass (n#pi^{+}) vs cos#theta;IM(n#pi^{+}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0 );
  new TH2F( "npim_IMvsCosT_Lab", "Invariant mass (n#pi^{-}) vs cos#theta;IM(n#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0 );
  new TH2F( "nk_IMvsCosT_Lab", "Invariant mass (nK^{-}) vs cos#theta;IM(nK^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0 );
  new TH2F( "np_IMvsCosT_Lab", "Invariant mass (np) vs cos#theta;IM(np) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0 );
  new TH2F( "nd_IMvsCosT_Lab", "Invariant mass (nd) vs cos#theta;IM(nd) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0 );
  new TH2F( "npip_IMvsCosT_CM", "Invariant mass (n#pi^{+}) vs cos#theta;IM(n#pi^{+}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0 );
  new TH2F( "npim_IMvsCosT_CM", "Invariant mass (n#pi^{-}) vs cos#theta;IM(n#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0 );
  new TH2F( "nk_IMvsCosT_CM", "Invariant mass (nK^{-}) vs cos#theta;IM(nK^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0 );
  new TH2F( "np_IMvsCosT_CM", "Invariant mass (np) vs cos#theta;IM(np) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0 );
  new TH2F( "nd_IMvsCosT_CM", "Invariant mass (nd) vs cos#theta;IM(nd) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0 );

  return true;
}
