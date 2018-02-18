// MyAnalysisK0.cpp

#include "MyAnalysisK0.h"

MyAnalysisK0::MyAnalysisK0(TFile* rt, ConfMan* conf)
{
  Initialize(conf);
  CutCondition();
  Clear();
}

MyAnalysisK0::~MyAnalysisK0()
{
  Clear();
  rtFile->cd();
  rtFile->Write();
  rtFile->Close();
}

void MyAnalysisK0::Clear()
{
}

bool MyAnalysisK0::DoAnalysis(ConfMan* conf, EventHeader* header, BeamLineHitMan* blMan, BeamLineTrackMan* bltrackMan, CDSHitMan* cdsMan, CDSTrackingMan* cdstrackMan, Particle* particle)
{

  rtFile->cd();

  DetectorList *dlist=DetectorList::GetInstance();
  if(particle->nBeam()!=1) return false;
  pBeam* beam = particle->beam(0);

  /* Event selection */
  if(particle->nPiplus()!=1) return false;
  if(particle->nPiminus()!=1) return false;
  if(particle->nCDS()!=2) return false;

  // ############### //
  // Vertex decision //
  // ############### //
  TVector3 vertex;
  double vdis = 9999;
  int fpro = -1;
  for(int it=0; it<particle->nProduct(); it++){
    pCDS* pro = particle->product(it);
    if(pro->comb()!=pow(2,CDS_PiPlus)+pow(2,CDS_PiMinus)) continue;
    if(vdis>9998||vdis>pro->vdis()){
      vdis = pro->vdis();
      vertex = pro->vbeam();
      if(GeomTools::GetID(vertex)==CID_Fiducial) fpro=it;
    }
  }
  if(fpro==-1) return false; /* there is no vertex in fiducial volume */

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

  pCDS* pip = particle->pip(0);
  pCDS* pim = particle->pim(0);
  pCDS* pro=0;
  if(fpro!=-1) pro=particle->product(fpro);

  /* CDS particle analysis */
  if(pro==0) return false;
  TVector3 ZeroV;
  TLorentzVector Ltgt; Ltgt.SetVectM(ZeroV,pMass);
  TLorentzVector Lbeam = beam->GetLorentzVector(vertex);
  TVector3 boost = (Ltgt+Lbeam).BoostVector();
  TLorentzVector cLtgt = Ltgt; cLtgt.Boost(-boost); TLorentzVector cLbeam = Lbeam; cLbeam.Boost(-boost);
  TLorentzVector Lpip = pip->GetLorentzVector();
  TLorentzVector cLpip = pip->GetLorentzVector(); cLpip.Boost(-boost);
  TLorentzVector Lpim = pim->GetLorentzVector();
  TLorentzVector cLpim = pim->GetLorentzVector(); cLpim.Boost(-boost);
  TLorentzVector Lpro = pro->GetLorentzVector();
  TLorentzVector cLpro = pro->GetLorentzVector(); cLpro.Boost(-boost);
  double pip_mass[2]  = { Lpip.M()                           , cLpip.M()};
  double pip_mom[2]   = { Lpip.Vect().Mag()                  , cLpip.Vect().Mag()};
  double pip_cost[2]  = { Lpip.Vect().CosTheta()             , cLpip.Vect().CosTheta()};
  double pip_phi[2]   = { Lpip.Vect().Phi()                  , cLpip.Vect().Phi()};
  double pip_mmass[2] = { (Ltgt+Lbeam-Lpip).M()              , (cLtgt+cLbeam-cLpip).M()};
  double pip_mmom[2]  = { (Ltgt+Lbeam-Lpip).Vect().Mag()     , (cLtgt+cLbeam-cLpip).Vect().Mag()};
  double pip_mcost[2] = { (Ltgt+Lbeam-Lpip).Vect().CosTheta(), (cLtgt+cLbeam-cLpip).Vect().CosTheta()};
  double pim_mass[2]  = { Lpim.M()                           , cLpim.M()};
  double pim_mom[2]   = { Lpim.Vect().Mag()                  , cLpim.Vect().Mag()};
  double pim_cost[2]  = { Lpim.Vect().CosTheta()             , cLpim.Vect().CosTheta()};
  double pim_phi[2]   = { Lpim.Vect().Phi()                  , cLpim.Vect().Phi()};
  double pim_mmass[2] = { (Ltgt+Lbeam-Lpim).M()              , (cLtgt+cLbeam-cLpim).M()};
  double pim_mmom[2]  = { (Ltgt+Lbeam-Lpim).Vect().Mag()     , (cLtgt+cLbeam-cLpim).Vect().Mag()};
  double pim_mcost[2] = { (Ltgt+Lbeam-Lpim).Vect().CosTheta(), (cLtgt+cLbeam-cLpim).Vect().CosTheta()};
  double pro_mass[2]  = { Lpro.M()                           , cLpro.M()};
  double pro_mom[2]   = { Lpro.Vect().Mag()                  , cLpro.Vect().Mag()};
  double pro_cost[2]  = { Lpro.Vect().CosTheta()             , cLpro.Vect().CosTheta()};
  double pro_phi[2]   = { Lpro.Vect().Phi()                  , cLpro.Vect().Phi()};
  double pro_mmass[2] = { (Ltgt+Lbeam-Lpro).M()              , (cLtgt+cLbeam-cLpro).M()};
  double pro_mmom[2]  = { (Ltgt+Lbeam-Lpro).Vect().Mag()     , (cLtgt+cLbeam-cLpro).Vect().Mag()};
  double pro_mcost[2] = { (Ltgt+Lbeam-Lpro).Vect().CosTheta(), (cLtgt+cLbeam-cLpro).Vect().CosTheta()};
  double pro_pbdca = pro->pbdca();
  double pro_vbdis = pro->vbdis();
  double pro_vdis  = pro->vdis();
  TString pid[8] = {"pip","p","d","t","he","pim","k","o"};
  TString frame[2] = {"Lab","CM"};

  /* PiPlus and PiMinus */
  for(int f=1; f<2; f++){
    FillHist(Form("pip_CosTvsMomentum_%s",frame[f].Data()),pip_cost[f],pip_mom[f]);
    FillHist(Form("pip_CosTvsPhi_%s",frame[f].Data()),pip_cost[f],pip_phi[f]);
    FillHist(Form("pim_CosTvsMomentum_%s",frame[f].Data()),pim_cost[f],pim_mom[f]);
    FillHist(Form("pim_CosTvsPhi_%s",frame[f].Data()),pim_cost[f],pim_phi[f]);
    FillHist(Form("pKpip_MMvsMomentum_%s",frame[f].Data()),pip_mmass[f],pip_mmom[f]);
    FillHist(Form("pKpip_MMvsCosT_%s",frame[f].Data()),pip_mmass[f],pip_mcost[f]);
    FillHist(Form("pKpim_MMvsMomentum_%s",frame[f].Data()),pim_mmass[f],pim_mmom[f]);
    FillHist(Form("pKpim_MMvsCosT_%s",frame[f].Data()),pim_mmass[f],pim_mcost[f]);
    FillHist(Form("pipMomentumvspimMomentum_%s",frame[f].Data()),pip_mom[f],pim_mom[f]);
    FillHist(Form("pipCosTvspimCosT_%s",frame[f].Data()),pip_cost[f],pim_cost[f]);
    FillHist(Form("pipPhivspimPhi_%s",frame[f].Data()),pip_phi[f],pim_phi[f]);
  }
  /* Production */
  for(int f=1; f<2; f++){
    FillHist(Form("pipi_IMvsMomentum_%s",frame[f].Data()),pro_mass[f],pro_mom[f]);
    FillHist(Form("pipi_IMvsCosT_%s",frame[f].Data()),pro_mass[f],pro_cost[f]);
    FillHist(Form("pipi_CosTvsMomentum_%s",frame[f].Data()),pro_cost[f],pro_mom[f]);
    FillHist(Form("pipi_CosTvsPhi_%s",frame[f].Data()),pro_cost[f],pro_phi[f]);
    FillHist(Form("pKpipi_MMvsMomentum_%s",frame[f].Data()),pro_mmass[f],pro_mmom[f]);
    FillHist(Form("pKpipi_MMvsCosT_%s",frame[f].Data()),pro_mmass[f],pro_mcost[f]);
  }
  FillHist("Pro_PBDCA",pro_pbdca);
  FillHist("Pro_VBDIS",pro_vbdis);
  FillHist("Pro_VDIS",pro_vdis);

  bool mnflag = false; /* missing neutron flag */
  if(mnll<pro_mmass[1]&&pro_mmass[1]<mnul) mnflag = true;
  bool k0flag = false; /* k0 flag */
  if(k0ll<pro_mass[1]&&pro_mass[1]<k0ul) k0flag = true;
  if(k0flag){
    for(int f=1; f<2; f++){
      FillHist(Form("pipi_k0_IMvsMomentum_%s",frame[f].Data()),pro_mass[f],pro_mom[f]);
      FillHist(Form("pipi_k0_IMvsCosT_%s",frame[f].Data()),pro_mass[f],pro_cost[f]);
      FillHist(Form("pKpipi_k0_MMvsMomentum_%s",frame[f].Data()),pro_mmass[f],pro_mmom[f]);
      FillHist(Form("pKpipi_k0_MMvsCosT_%s",frame[f].Data()),pro_mmass[f],pro_mcost[f]);
    }
  }
  if(mnflag){
    for(int f=1; f<2; f++){
      FillHist(Form("pip_mn_CosTvsMomentum_%s",frame[f].Data()),pip_cost[f],pip_mom[f]);
      FillHist(Form("pip_mn_CosTvsPhi_%s",frame[f].Data()),pip_cost[f],pip_phi[f]);
      FillHist(Form("pim_mn_CosTvsMomentum_%s",frame[f].Data()),pim_cost[f],pim_mom[f]);
      FillHist(Form("pim_mn_CosTvsPhi_%s",frame[f].Data()),pim_cost[f],pim_phi[f]);
      FillHist(Form("pipi_mn_IMvsMomentum_%s",frame[f].Data()),pro_mass[f],pro_mom[f]);
      FillHist(Form("pipi_mn_IMvsCosT_%s",frame[f].Data()),pro_mass[f],pro_cost[f]);
      FillHist(Form("pKpipi_mn_MMvsMomentum_%s",frame[f].Data()),pro_mmass[f],pro_mmom[f]);
      FillHist(Form("pKpipi_mn_MMvsCosT_%s",frame[f].Data()),pro_mmass[f],pro_mcost[f]);
    }
  }

  pNC*  nc =0;
  if(fnc!=-1)  nc =particle->nc(fnc);
  /* with forward neutron */
  if(nc==0) return false;
  bool nflag = true;  /* forward neutron flag */
  TLorentzVector Ln  = nc ->GetLorentzVector();
  TLorentzVector cLn  = nc ->GetLorentzVector(); cLn.Boost(-boost);
  double nc_mom[2]     = { Ln.Vect().Mag()                  , cLn.Vect().Mag()};
  double nc_cost[2]    = { Ln.Vect().CosTheta()             , cLn.Vect().CosTheta()};
  double nc_phi[2]     = { Ln.Vect().Phi()                  , cLn.Vect().Phi()};
  double nc_mmass[2]   = { (Ltgt+Lbeam-Ln).M()              , (cLtgt+cLbeam-cLn).M()};
  double nc_mmom[2]    = { (Ltgt+Lbeam-Ln).Vect().Mag()     , (cLtgt+cLbeam-cLn).Vect().Mag()};
  double nc_mcost[2]   = { (Ltgt+Lbeam-Ln).Vect().CosTheta(), (cLtgt+cLbeam-cLn).Vect().CosTheta()};
  double npip_mass[2]  = { (Lpip+Ln).M()                           , (cLpip+cLn).M()};
  double npip_mom[2]   = { (Lpip+Ln).Vect().Mag()                  , (cLpip+cLn).Vect().Mag()};
  double npip_cost[2]  = { (Lpip+Ln).Vect().CosTheta()             , (cLpip+cLn).Vect().CosTheta()};
  double npip_phi[2]   = { (Lpip+Ln).Vect().Phi()                  , (cLpip+cLn).Vect().Phi()};
  double npip_mmass[2] = { (Ltgt+Lbeam-(Lpip+Ln)).M()              , (cLtgt+cLbeam-(cLpip+cLn)).M()};
  double npip_mmom[2]  = { (Ltgt+Lbeam-(Lpip+Ln)).Vect().Mag()     , (cLtgt+cLbeam-(cLpip+cLn)).Vect().Mag()};
  double npip_mcost[2] = { (Ltgt+Lbeam-(Lpip+Ln)).Vect().CosTheta(), (cLtgt+cLbeam-(cLpip+cLn)).Vect().CosTheta()};
  double npim_mass[2]  = { (Lpim+Ln).M()                           , (cLpim+cLn).M()};
  double npim_mom[2]   = { (Lpim+Ln).Vect().Mag()                  , (cLpim+cLn).Vect().Mag()};
  double npim_cost[2]  = { (Lpim+Ln).Vect().CosTheta()             , (cLpim+cLn).Vect().CosTheta()};
  double npim_phi[2]   = { (Lpim+Ln).Vect().Phi()                  , (cLpim+cLn).Vect().Phi()};
  double npim_mmass[2] = { (Ltgt+Lbeam-(Lpim+Ln)).M()              , (cLtgt+cLbeam-(cLpim+cLn)).M()};
  double npim_mmom[2]  = { (Ltgt+Lbeam-(Lpim+Ln)).Vect().Mag()     , (cLtgt+cLbeam-(cLpim+cLn)).Vect().Mag()};
  double npim_mcost[2] = { (Ltgt+Lbeam-(Lpim+Ln)).Vect().CosTheta(), (cLtgt+cLbeam-(cLpim+cLn)).Vect().CosTheta()};
  for(int f=1; f<2; f++){
    FillHist(Form("n_CosTvsMomentum_%s",frame[f].Data()),nc_cost[f],nc_mom[f]);
    FillHist(Form("n_CosTvsPhi_%s",frame[f].Data()),nc_cost[f],nc_phi[f]);
    FillHist(Form("pKn_MMvsMomentum_%s",frame[f].Data()),nc_mmass[f],nc_mmom[f]);
    FillHist(Form("pKn_MMvsCosT_%s",frame[f].Data()),nc_mmass[f],nc_mcost[f]);
  }
  for(int f=1; f<2; f++){
    FillHist(Form("npip_IMvsMomentum_%s",frame[f].Data()),npip_mass[f],npip_mom[f]);
    FillHist(Form("npip_IMvsCosT_%s",frame[f].Data()),npip_mass[f],npip_cost[f]);
    FillHist(Form("npim_IMvsMomentum_%s",frame[f].Data()),npim_mass[f],npim_mom[f]);
    FillHist(Form("npim_IMvsCosT_%s",frame[f].Data()),npim_mass[f],npim_cost[f]);
    FillHist(Form("pKnpip_MMvsMomentum_%s",frame[f].Data()),npip_mmass[f],npip_mmom[f]);
    FillHist(Form("pKnpip_MMvsCosT_%s",frame[f].Data()),npip_mmass[f],npip_mcost[f]);
    FillHist(Form("pKnpim_MMvsMomentum_%s",frame[f].Data()),npim_mmass[f],npim_mmom[f]);
    FillHist(Form("pKnpim_MMvsCosT_%s",frame[f].Data()),npim_mmass[f],npim_mcost[f]);
  }
  if(nflag){
    for(int f=1; f<2; f++){
      FillHist(Form("pipi_n_IMvsMomentum_%s",frame[f].Data()),pro_mass[f],pro_mom[f]);
      FillHist(Form("pipi_n_IMvsCosT_%s",frame[f].Data()),pro_mass[f],pro_cost[f]);
      FillHist(Form("pKpipi_n_MMvsMomentum_%s",frame[f].Data()),pro_mmass[f],pro_mmom[f]);
      FillHist(Form("pKpipi_n_MMvsCosT_%s",frame[f].Data()),pro_mmass[f],pro_mcost[f]);
    }
  }
  if(k0flag){
    for(int f=1; f<2; f++){
      FillHist(Form("n_k0_CosTvsMomentum_%s",frame[f].Data()),nc_cost[f],nc_mom[f]);
      FillHist(Form("n_k0_CosTvsPhi_%s",frame[f].Data()),nc_cost[f],nc_phi[f]);
      FillHist(Form("pKn_k0_MMvsMomentum_%s",frame[f].Data()),nc_mmass[f],nc_mmom[f]);
      FillHist(Form("pKn_k0_MMvsCosT_%s",frame[f].Data()),nc_mmass[f],nc_mcost[f]);
    }
    for(int f=1; f<2; f++){
      FillHist(Form("npip_k0_IMvsMomentum_%s",frame[f].Data()),npip_mass[f],npip_mom[f]);
      FillHist(Form("npip_k0_IMvsCosT_%s",frame[f].Data()),npip_mass[f],npip_cost[f]);
      FillHist(Form("npim_k0_IMvsMomentum_%s",frame[f].Data()),npim_mass[f],npim_mom[f]);
      FillHist(Form("npim_k0_IMvsCosT_%s",frame[f].Data()),npim_mass[f],npim_cost[f]);
      FillHist(Form("pKnpip_k0_MMvsMomentum_%s",frame[f].Data()),npip_mmass[f],npip_mmom[f]);
      FillHist(Form("pKnpip_k0_MMvsCosT_%s",frame[f].Data()),npip_mmass[f],npip_mcost[f]);
      FillHist(Form("pKnpim_k0_MMvsMomentum_%s",frame[f].Data()),npim_mmass[f],npim_mmom[f]);
      FillHist(Form("pKnpim_k0_MMvsCosT_%s",frame[f].Data()),npim_mmass[f],npim_mcost[f]);
    }
  }

  return true;

}

bool MyAnalysisK0::FillHist(TString name, double val1, int weight)
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

bool MyAnalysisK0::FillHist(TString name, TString val1, int weight)
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

bool MyAnalysisK0::FillHist(TString name, double val1, double val2, int weight)
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

bool MyAnalysisK0::FillHist(TString name, TString val1, TString val2, int weight)
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

void MyAnalysisK0::CutCondition()
{
  /* K0 mass cut */ 
  k0ll = 0.470; k0ul = 0.525;
  /* Missing neutron mass cut */ 
  mnll = 0.910; mnul = 0.970;
}

bool MyAnalysisK0::Initialize(ConfMan* confMan)
{
  rtFile =  new TFile( Form("anaK0_%s",confMan->GetOutFileName().c_str()), "RECREATE");
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
  //new TH2F("pip_CosTvsMomentum_Lab","cos#theta vs. momentum of #pi^{+};cos#theta;Momentum (GeV/c)", 2000, -1.0, 1.0, 1000, 0.0, 1.0);
  //new TH2F("pip_CosTvsPhi_Lab","cos#theta vs. #phi of #pi^{+};cos#theta;#phi", 2000, -1.0, 1.0, 1600, 0.0, 3.2);
  //new TH2F("pim_CosTvsMomentum_Lab","cos#theta vs. momentum of #pi^{-};cos#theta;Momentum (GeV/c)", 2000, -1.0, 1.0, 1000, 0.0, 1.0);
  //new TH2F("pim_CosTvsPhi_Lab","cos#theta vs. #phi of #pi^{-};cos#theta;#phi", 2000, -1.0, 1.0, 1600, 0.0, 3.2);
  //new TH2F("Pro_CosTvsMomentum_Lab","cos#theta vs. momentum of Production;cos#theta;Momentum (GeV/c)", 2000, -1.0, 1.0, 1000, 0.0, 1.0);
  //new TH2F("Pro_CosTvsPhi_Lab","cos#theta vs. #phi of Production;cos#theta;#phi", 2000, -1.0, 1.0, 1600, 0.0, 3.2);
  //new TH2F("pipMomentumvspimMomentum_Lab","#pi^{+} momentum vs. #pi^{-} momentum;#pi^{+} momentum (GeV/c);#pi^{-} momentum (GeV/c)", 1000, 0.0, 1.0, 1000, 0.0, 1.0);
  //new TH2F("pipCosTvspimCosT_Lab","#pi^{+} cos#theta vs. #pi^{-} cos#theta;cos#theta_{#pi^{+}};cos#theta_{#pi^{-}}", 2000, -1.0, 1.0, 2000, -1.0, 1.0);
  //new TH2F("pipPhivspimPhi_Lab","#pi^{+} #phi vs. #pi^{-} #phi;#phi_#pi^{+}};#phi_{#pi^{-}}", 1600, 0.0, 3.2, 1600, 0.0, 3.2);
  //new TH2F("n_CosTvsMomentum_Lab","cos#theta vs. momentum of n;cos#theta;Momentum (GeV/c)", 2000, -1.0, 1.0, 1000, 0.0, 1.0);
  //new TH2F("n_CosTvsPhi_Lab","cos#theta vs. #phi of n;cos#theta;#phi", 2000, -1.0, 1.0, 1600, 0.0, 3.2);
  new TH2F("pipMomentumvspimMomentum_CM","#pi^{+} momentum vs. #pi^{-} momentum;#pi^{+} momentum (GeV/c);#pi^{-} momentum (GeV/c)", 1000, 0.0, 1.0, 1000, 0.0, 1.0);
  new TH2F("pipCosTvspimCosT_CM","#pi^{+} cos#theta vs. #pi^{-} cos#theta;cos#theta_{#pi^{+}};cos#theta_{#pi^{-}}", 2000, -1.0, 1.0, 2000, -1.0, 1.0);
  new TH2F("pipPhivspimPhi_CM","#pi^{+} #phi vs. #pi^{-} #phi;#phi_{#pi^{+}};#phi_{#pi^{-}}", 1600, 0.0, 3.2, 1600, 0.0, 3.2);
  new TH2F("pip_CosTvsMomentum_CM","cos#theta vs. momentum of #pi^{+};cos#theta;Momentum (GeV/c)", 2000, -1.0, 1.0, 1000, 0.0, 1.0);
  new TH2F("pip_CosTvsPhi_CM","cos#theta vs. #phi of #pi^{+};cos#theta;#phi", 2000, -1.0, 1.0, 1600, 0.0, 3.2);
  new TH2F("pim_CosTvsMomentum_CM","cos#theta vs. momentum of #pi^{-};cos#theta;Momentum (GeV/c)", 2000, -1.0, 1.0, 1000, 0.0, 1.0);
  new TH2F("pim_CosTvsPhi_CM","cos#theta vs. #phi of #pi^{-};cos#theta;#phi", 2000, -1.0, 1.0, 1600, 0.0, 3.2);
  new TH2F("pipi_CosTvsMomentum_CM","cos#theta vs. momentum of Production;cos#theta;Momentum (GeV/c)", 2000, -1.0, 1.0, 1000, 0.0, 1.0);
  new TH2F("pipi_CosTvsPhi_CM","cos#theta vs. #phi of Production;cos#theta;#phi", 2000, -1.0, 1.0, 1600, 0.0, 3.2);
  new TH1F("Pro_VDIS","Distance between tow vertexes;Distance (cm);Counts",100,0.0,10.0);
  new TH1F("Pro_VBDIS","distance between beam and vertex;Distance (cm);Counts",100,0.0,10.0);
  new TH1F("Pro_PBDCA","DCA between beam and product;DCA (cm);Counts",100,0.0,10.0);
  new TH2F("n_CosTvsMomentum_CM","cos#theta vs. momentum of n;cos#theta;Momentum (GeV/c)", 2000, -1.0, 1.0, 1000, 0.0, 1.0);
  new TH2F("n_CosTvsPhi_CM","cos#theta vs. #phi of n;cos#theta;#phi", 2000, -1.0, 1.0, 1600, 0.0, 3.2);
  // FWD neutral Particles
  std::cout << "Define Histograms for FWD neutral particles" << std::endl;
  new TH2F( "OverbetavsMomentum_FWD", "1/#beta vs. Momentum;1/#beta;Momentum (GeV/c)", 1000, -0, 10, 500, -5, 5 );
  new TH2F( "OverbetavsEnergy_FWD", "1/#beta vs. Energy deposit;1/#beta;Energy deposit (MeVee)", 1000, -0, 10, 5000, 0, 100);
  new TH2F( "TOFvsMomentum_FWD", "TOF vs. Momentum;TOF (ns);Momentum (GeV/c)", 4000, 30, 70, 500, -5, 5 );
  new TH2F( "NCHitPosition_XY", "Hit position at NC (X vs. Y);X position (cm);Y position (cm)",16,-160,160,15,-75,75);
  new TH2F( "NCHitSegment", "NC segment;Segment in layer;Layer",16,0.5,16.5,7,0.5,7.5);
  // Invariant mass of production
  std::cout << "Define Histograms for product I.M. with FWD" << std::endl;
  //new TH2F( "pipi_IMvsMomentum_Lab", "Invariant mass (#pi^{+}#pi^{-}) vs Momentum;IM(#pi^{+}#pi^{-}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0 );
  //new TH2F( "pipi_IMvsCosT_Lab", "Invariant mass (#pi^{+}#pi^{-}) vs cos#theta;IM(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0 );
  //new TH2F( "K0_IMvsMomentum_Lab", "Invariant mass (#pi^{+}#pi^{-}) vs Momentum;IM(#pi^{+}#pi^{-}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0 );
  //new TH2F( "K0_IMvsCosT_Lab", "Invariant mass (#pi^{+}#pi^{-}) vs cos#theta;IM(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0 );
  //new TH2F( "pipi_mn_IMvsMomentum_Lab", "Invariant mass (#pi^{+}#pi^{-}) vs Momentum;IM(#pi^{+}#pi^{-}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0 );
  //new TH2F( "pipi_mn_IMvsCosT_Lab", "Invariant mass (#pi^{+}#pi^{-}) vs cos#theta;IM(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0 );
  //new TH2F( "npip_IMvsMomentum_Lab", "Invariant mass (n#pi^{+}) vs Momentum;IM(n#pi^{+}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0 );
  //new TH2F( "npip_IMvsCosT_Lab", "Invariant mass (n#pi^{+}) vs cos#theta;IM(n#pi^{+}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0 );
  //new TH2F( "npim_IMvsMomentum_Lab", "Invariant mass (n#pi^{-}) vs Momentum;IM(n#pi^{-}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0 );
  //new TH2F( "npim_IMvsCosT_Lab", "Invariant mass (n#pi^{-}) vs cos#theta;IM(n#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0 );
  new TH2F( "pipi_IMvsMomentum_CM", "Invariant mass (#pi^{+}#pi^{-}) vs Momentum;IM(#pi^{+}#pi^{-}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0 );
  new TH2F( "pipi_IMvsCosT_CM", "Invariant mass (#pi^{+}#pi^{-}) vs cos#theta;IM(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0 );
  new TH2F( "pipi_k0_IMvsMomentum_CM", "Invariant mass (#pi^{+}#pi^{-}) vs Momentum;IM(#pi^{+}#pi^{-}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0 );
  new TH2F( "pipi_k0_IMvsCosT_CM", "Invariant mass (#pi^{+}#pi^{-}) vs cos#theta;IM(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0 );
  new TH2F( "pipi_mn_IMvsMomentum_CM", "Invariant mass (#pi^{+}#pi^{-}) vs Momentum;IM(#pi^{+}#pi^{-}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0 );
  new TH2F( "pipi_mn_IMvsCosT_CM", "Invariant mass (#pi^{+}#pi^{-}) vs cos#theta;IM(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0 );
  new TH2F( "pipi_n_IMvsMomentum_CM", "Invariant mass (#pi^{+}#pi^{-}) vs Momentum;IM(#pi^{+}#pi^{-}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0 );
  new TH2F( "pipi_n_IMvsCosT_CM", "Invariant mass (#pi^{+}#pi^{-}) vs cos#theta;IM(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0 );
  new TH2F( "npip_IMvsMomentum_CM", "Invariant mass (n#pi^{+}) vs Momentum;IM(n#pi^{+}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0 );
  new TH2F( "npip_IMvsCosT_CM", "Invariant mass (n#pi^{+}) vs cos#theta;IM(n#pi^{+}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0 );
  new TH2F( "npim_IMvsMomentum_CM", "Invariant mass (n#pi^{-}) vs Momentum;IM(n#pi^{-}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0 );
  new TH2F( "npim_IMvsCosT_CM", "Invariant mass (n#pi^{-}) vs cos#theta;IM(n#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0 );
  new TH2F( "npip_k0_IMvsMomentum_CM", "Invariant mass (n#pi^{+}) vs Momentum;IM(n#pi^{+}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0 );
  new TH2F( "npip_k0_IMvsCosT_CM", "Invariant mass (n#pi^{+}) vs cos#theta;IM(n#pi^{+}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0 );
  new TH2F( "npim_k0_IMvsMomentum_CM", "Invariant mass (n#pi^{-}) vs Momentum;IM(n#pi^{-}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0 );
  new TH2F( "npim_k0_IMvsCosT_CM", "Invariant mass (n#pi^{-}) vs cos#theta;IM(n#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0 );
  // Missing mass of production
  std::cout << "Define Histograms for p(K-,n)X M.M." << std::endl;
  //new TH2F("pKpip_MMvsMomentum_Lab", "p(K^{-},#pi^{+})X missing mass vs. missing momentum;MM_{p}(#pi^{+}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  //new TH2F("pKpip_MMvsCosT_Lab", "p(K^{-},#pi^{+})X missing mass vs. cos#theta;MM_{p}(#pi^{+}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  //new TH2F("pKpim_MMvsMomentum_Lab", "p(K^{-},#pi^{-})X missing mass vs. missing momentum;MM_{p}(#pi^{-}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  //new TH2F("pKpim_MMvsCosT_Lab", "p(K^{-},#pi^{-})X missing mass vs. cos#theta;MM_{p}(#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  //new TH2F("pKpipi_MMvsMomentum_Lab", "p(K^{-},#pi^{+}#pi^{-})X missing mass vs. missing momentum;MM_{p}(#pi^{+}#pi^{-}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  //new TH2F("pKpipi_MMvsCosT_Lab", "p(K^{-},#pi^{+}#pi^{-})X missing mass vs. cos#theta;MM_{p}(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  //new TH2F("pKK0_MMvsMomentum_Lab", "p(K^{-},#pi^{+}#pi^{-})X missing mass vs. missing momentum;MM_{p}(#pi^{+}#pi^{-}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  //new TH2F("pKK0_MMvsCosT_Lab", "p(K^{-},#pi^{+}#pi^{-})X missing mass vs. cos#theta;MM_{p}(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  //new TH2F("pKpipi_mn_MMvsMomentum_Lab", "p(K^{-},#pi^{+}#pi^{-})X missing mass vs. missing momentum;MM_{p}(#pi^{+}#pi^{-}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  //new TH2F("pKpipi_mn_MMvsCosT_Lab", "p(K^{-},#pi^{+}#pi^{-})X missing mass vs. cos#theta;MM_{p}(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  //new TH2F("pKn_MMvsMomentum_Lab", "p(K^{-},n)X missing mass vs. missing momentum;MM_{p}(n) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  //new TH2F("pKn_MMvsCosT_Lab", "p(K^{-},n)X missing mass vs. cos#theta;MM_{p}(n) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  //new TH2F("pKnpip_MMvsMomentum_Lab", "p(K^{-},n#pi^{+})X missing mass vs. missing momentum;MM_{p}(npi^{+}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  //new TH2F("pKnpip_MMvsCosT_Lab", "p(K^{-},n#pi^{+})X missing mass vs. cos#theta;MM_{p}(npi^{+}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  //new TH2F("pKnpim_MMvsMomentum_Lab", "p(K^{-},n#pi^{-})X missing mass vs. missing momentum;MM_{p}(npi^{-}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  //new TH2F("pKnpim_MMvsCosT_Lab", "p(K^{-},n#pi^{-})X missing mass vs. cos#theta;MM_{p}(npi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKpip_MMvsMomentum_CM", "p(K^{-},#pi^{+})X missing mass vs. missing momentum;MM_{p}(#pi^{+}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKpip_MMvsCosT_CM", "p(K^{-},#pi^{+})X missing mass vs. cos#theta;MM_{p}(#pi^{+}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKpim_MMvsMomentum_CM", "p(K^{-},#pi^{-})X missing mass vs. missing momentum;MM_{p}(#pi^{-}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKpim_MMvsCosT_CM", "p(K^{-},#pi^{-})X missing mass vs. cos#theta;MM_{p}(#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKpipi_MMvsMomentum_CM", "p(K^{-},#pi^{+}#pi^{-})X missing mass vs. missing momentum;MM_{p}(#pi^{+}#pi^{-}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKpipi_MMvsCosT_CM", "p(K^{-},#pi^{+}#pi^{-})X missing mass vs. cos#theta;MM_{p}(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKpipi_k0_MMvsMomentum_CM", "p(K^{-},#pi^{+}#pi^{-})X missing mass vs. missing momentum;MM_{p}(#pi^{+}#pi^{-}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKpipi_k0_MMvsCosT_CM", "p(K^{-},#pi^{+}#pi^{-})X missing mass vs. cos#theta;MM_{p}(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKpipi_mn_MMvsMomentum_CM", "p(K^{-},#pi^{+}#pi^{-})X missing mass vs. missing momentum;MM_{p}(#pi^{+}#pi^{-}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKpipi_mn_MMvsCosT_CM", "p(K^{-},#pi^{+}#pi^{-})X missing mass vs. cos#theta;MM_{p}(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKpipi_n_MMvsMomentum_CM", "p(K^{-},#pi^{+}#pi^{-})X missing mass vs. missing momentum;MM_{p}(#pi^{+}#pi^{-}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKpipi_n_MMvsCosT_CM", "p(K^{-},#pi^{+}#pi^{-})X missing mass vs. cos#theta;MM_{p}(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKn_MMvsMomentum_CM", "p(K^{-},n)X missing mass vs. missing momentum;MM_{p}(n) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKn_MMvsCosT_CM", "p(K^{-},n)X missing mass vs. cos#theta;MM_{p}(n) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKnpip_MMvsMomentum_CM", "p(K^{-},n#pi^{+})X missing mass vs. missing momentum;MM_{p}(npi^{+}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKnpip_MMvsCosT_CM", "p(K^{-},n#pi^{+})X missing mass vs. cos#theta;MM_{p}(npi^{+}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKnpim_MMvsMomentum_CM", "p(K^{-},n#pi^{-})X missing mass vs. missing momentum;MM_{p}(npi^{-}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKnpim_MMvsCosT_CM", "p(K^{-},n#pi^{-})X missing mass vs. cos#theta;MM_{p}(npi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKn_k0_MMvsMomentum_CM", "p(K^{-},n)X missing mass vs. missing momentum;MM_{p}(n) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKn_k0_MMvsCosT_CM", "p(K^{-},n)X missing mass vs. cos#theta;MM_{p}(n) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKnpip_k0_MMvsMomentum_CM", "p(K^{-},n#pi^{+})X missing mass vs. missing momentum;MM_{p}(npi^{+}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKnpip_k0_MMvsCosT_CM", "p(K^{-},n#pi^{+})X missing mass vs. cos#theta;MM_{p}(npi^{+}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKnpim_k0_MMvsMomentum_CM", "p(K^{-},n#pi^{-})X missing mass vs. missing momentum;MM_{p}(npi^{-}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKnpim_k0_MMvsCosT_CM", "p(K^{-},n#pi^{-})X missing mass vs. cos#theta;MM_{p}(npi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);

  return true;
}
