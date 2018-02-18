// MyAnalysis3HeKn.cpp

#include "MyAnalysis3HeKn.h"

MyAnalysis3HeKn::MyAnalysis3HeKn(TFile* rt, ConfMan* conf)
{
  Initialize(conf);
  CutCondition();
  Clear();
}

MyAnalysis3HeKn::~MyAnalysis3HeKn()
{
  Clear();
  rtFile->cd();
  rtFile->Write();
  rtFile->Close();
}

void MyAnalysis3HeKn::Clear()
{
}

bool MyAnalysis3HeKn::DoAnalysis(ConfMan* conf, EventHeader* header, BeamLineHitMan* blMan, BeamLineTrackMan* bltrackMan, CDSHitMan* cdsMan, CDSTrackingMan* cdstrackMan, Particle* particle)
{
  rtFile->cd();

  FillHist("EventNumber",0);
  DetectorList *dlist=DetectorList::GetInstance();
  if(particle->nBeam()!=1) return false;
  pBeam* beam = particle->beam(0);
  FillHist("EventNumber",1);

  /* Event selection */
  //if(particle->nCDS()!=2) return true;

  // ############### //
  // Vertex decision //
  // ############### //
  TVector3 vertex;
  double vdis = 9999;
  int fcds = -1;
  for(int it=0; it<particle->nCDS(); it++){
    pCDS* cds = particle->cdsi(it);
    if(vdis>9998||vdis>cds->vdis()){
      vdis = cds->vdis();
      vertex = cds->vbeam();
      if(GeomTools::GetID(vertex)==CID_Fiducial) fcds=it;
    }
  }
  if(fcds==-1) return true; /* there is no vertex in fiducial volume */
  FillHist("EventNumber",2);

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

  // Trigger Pattern
  FillHist("TriggerPattern",0);
  for( int i=0; i<20; i++ ){
    int val = header->pattern(i);
    if( 0<val ){
      FillHist("TriggerPattern",i);
    }
  }

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

  pCDS* cds=0;
  if(fcds!=-1) {
    cds=particle->cdsi(fcds);
  }
  if(cds==0) return false;

  TVector3 ZeroV;
  TLorentzVector Ltgt; Ltgt.SetVectM(ZeroV,ThreeHeMass);
  TLorentzVector Lbeam = beam->GetLorentzVector(vertex);
  TVector3 boost = (Ltgt+Lbeam).BoostVector();
  TLorentzVector cLtgt = Ltgt; cLtgt.Boost(-boost); TLorentzVector cLbeam = Lbeam; cLbeam.Boost(-boost);
  TLorentzVector Lcds = cds->GetLorentzVector();
  TLorentzVector cLcds = cds->GetLorentzVector(); cLcds.Boost(-boost);
  TString frame[2] = {"Lab","CM"};
  TString pname[8] = {"pip","p","d","t","he","pim","k","o"};

  pNC*  nc =0;
  if(fnc!=-1)  nc =particle->nc(fnc);
  if(nc==0) return true;
  FillHist("EventNumber",5);

  TLorentzVector Ln  = nc ->GetLorentzVector();
  TLorentzVector cLn  = nc ->GetLorentzVector(); cLn.Boost(-boost);
  double cds_mass[2]  = { (Lcds).M()                    , (cLcds).M()};
  double cds_mom[2]   = { (Lcds).Vect().Mag()           , (cLcds).Vect().Mag()};
  double cds_cost[2]  = { (Lcds).Vect().CosTheta()      , (cLcds).Vect().CosTheta()};
  double cds_phi[2]   = { (Lcds).Vect().Phi()           , (cLcds).Vect().Phi()};
  double cds_mmass[2]   = { (Ltgt+Lbeam-Lcds).M()              , (cLtgt+cLbeam-cLcds).M()};
  double cds_mmom[2]    = { (Ltgt+Lbeam-Lcds).Vect().Mag()     , (cLtgt+cLbeam-cLcds).Vect().Mag()};
  double cds_mcost[2]   = { (Ltgt+Lbeam-Lcds).Vect().CosTheta(), (cLtgt+cLbeam-cLcds).Vect().CosTheta()};
  double nc_mom[2]     = { Ln.Vect().Mag()                  , cLn.Vect().Mag()};
  double nc_cost[2]    = { Ln.Vect().CosTheta()             , cLn.Vect().CosTheta()};
  double nc_phi[2]     = { Ln.Vect().Phi()                  , cLn.Vect().Phi()};
  double nc_mmass[2]   = { (Ltgt+Lbeam-Ln).M()              , (cLtgt+cLbeam-cLn).M()};
  double nc_mmom[2]    = { (Ltgt+Lbeam-Ln).Vect().Mag()     , (cLtgt+cLbeam-cLn).Vect().Mag()};
  double nc_mcost[2]   = { (Ltgt+Lbeam-Ln).Vect().CosTheta(), (cLtgt+cLbeam-cLn).Vect().CosTheta()};
  double ncds_mass[2]  = { (Lcds+Ln).M()                           , (cLcds+cLn).M()};
  double ncds_mom[2]   = { (Lcds+Ln).Vect().Mag()                  , (cLcds+cLn).Vect().Mag()};
  double ncds_cost[2]  = { (Lcds+Ln).Vect().CosTheta()             , (cLcds+cLn).Vect().CosTheta()};
  double ncds_phi[2]   = { (Lcds+Ln).Vect().Phi()                  , (cLcds+cLn).Vect().Phi()};
  double ncds_mmass[2] = { (Ltgt+Lbeam-(Lcds+Ln)).M()              , (cLtgt+cLbeam-(cLcds+cLn)).M()};
  double ncds_mmom[2]  = { (Ltgt+Lbeam-(Lcds+Ln)).Vect().Mag()     , (cLtgt+cLbeam-(cLcds+cLn)).Vect().Mag()};
  double ncds_mcost[2] = { (Ltgt+Lbeam-(Lcds+Ln)).Vect().CosTheta(), (cLtgt+cLbeam-(cLcds+cLn)).Vect().CosTheta()};
  double seg = (nc->seg()-1)%16+1;
  double lay = (nc->seg()-1)/16+1;
  /* flag */
  bool k0flag = false; /* k0 flag */
  bool sblk0flag = false; /* lower sideband of k0 flag */
  bool sbuk0flag = false; /* upper sideband of k0 flag */
  for(int ipro=0; ipro<particle->nProduct(); ipro++){
    pCDS* pro = particle->product(ipro);
    if(pro->comb()!=pow(2,CDS_PiPlus)+pow(2,CDS_PiMinus)) continue;
    if(GeomTools::GetID(pro->vbeam())!=CID_Fiducial) continue;
    TLorentzVector Lpro = pro->GetLorentzVector();
    Lpro.Boost(-boost);
    double pro_mass = Lpro.M();
    double pro_mom  = Lpro.Vect().Mag();
    FillHist("pipi_n_IMvsMomentum_CM",pro_mass,pro_mom,1);
    if(k0ll<pro_mass&&pro_mass<k0ul) k0flag = true;
    if(sblk0ll<pro_mass&&pro_mass<sblk0ul) sblk0flag = true;
    if(sbuk0ll<pro_mass&&pro_mass<sbuk0ul) sbuk0flag = true;
  }
  bool spflag = false; /* sp flag */
  bool sblspflag = false; /* lower sideband of sp flag */
  bool sbuspflag = false; /* upper sideband of sp flag */
  for(int ipip=0; ipip<particle->nPiplus(); ipip++){
    pCDS* pip = particle->pip(ipip);
    if(GeomTools::GetID(pip->vbeam())!=CID_Fiducial) continue;
    TLorentzVector Lpip = pip->GetLorentzVector();
    TLorentzVector Lnpip = Lpip+Ln;
    Lnpip.Boost(-boost);
    double pro_mass = Lnpip.M();
    double pro_mom = Lnpip.Vect().Mag();
    if(spll<pro_mass&&pro_mass<spul) spflag = true;
    if(sblspll<pro_mass&&pro_mass<sblspul) sblspflag = true;
    if(sbuspll<pro_mass&&pro_mass<sbuspul) sbuspflag = true;
    FillHist("npip_npip_IMvsMomentum_CM",pro_mass,pro_mom,1);
  }
  bool smflag = false; /* sm flag */
  bool sblsmflag = false; /* lower sideband of sm flag */
  bool sbusmflag = false; /* upper sideband of sm flag */
  for(int ipim=0; ipim<particle->nPiminus(); ipim++){
    pCDS* pim = particle->pim(ipim);
    if(GeomTools::GetID(pim->vbeam())!=CID_Fiducial) continue;
    TLorentzVector Lpim = pim->GetLorentzVector();
    TLorentzVector Lnpim = Lpim+Ln;
    Lnpim.Boost(-boost);
    double pro_mass = Lnpim.M();
    double pro_mom = Lnpim.Vect().Mag();
    if(smll<pro_mass&&pro_mass<smul) smflag = true;
    if(sblsmll<pro_mass&&pro_mass<sblsmul) sblsmflag = true;
    if(sbusmll<pro_mass&&pro_mass<sbusmul) sbusmflag = true;
    FillHist("npim_npim_IMvsMomentum_CM",pro_mass,pro_mom,1);
  }

  FillHist("FWD_n_OverbetavsMomentum",1.0/nc->beta(),nc->mom().Mag());
  FillHist("FWD_n_OverbetavsEnergy",1.0/nc->beta(),nc->energy());
  FillHist("FWD_n_TOFvsMomentum",nc->tof(),nc->mom().Mag());
  FillHist("FWD_n_HitPosition",nc->hitpos().X(),nc->hitpos().Y());
  FillHist("FWD_n_HitSegment",seg,lay);

  // CDS particles
  for(int f=1; f<2; f++){
    FillHist(Form("%s_CosTvsMomentum_%s",pname[cds->pid()].Data(),frame[f].Data()),cds_cost[f],cds_mom[f]);
    FillHist(Form("%s_CosTvsPhi_%s",pname[cds->pid()].Data(),frame[f].Data()),cds_cost[f],cds_phi[f]);
    FillHist(Form("pK%s_MMvsMomentum_%s",pname[cds->pid()].Data(),frame[f].Data()),cds_mmass[f],cds_mmom[f]);
    FillHist(Form("pK%s_MMvsCosT_%s",pname[cds->pid()].Data(),frame[f].Data()),cds_mmass[f],cds_mcost[f]);
  }


  // Missing mass 
  for(int f=1; f<2; f++){
    FillHist(Form("n_CosTvsMomentum_%s",frame[f].Data()),nc_cost[f],nc_mom[f]);
    FillHist(Form("n_CosTvsPhi_%s",frame[f].Data()),nc_cost[f],nc_phi[f]);
    FillHist(Form("HeKn_MMvsMomentum_%s",frame[f].Data()),nc_mmass[f],nc_mmom[f]);
    FillHist(Form("HeKn_MMvsCosT_%s",frame[f].Data()),nc_mmass[f],nc_mcost[f]);
  }
  if(particle->nPiplus()>0){
    for(int f=1; f<2; f++){
      FillHist(Form("n_pip_CosTvsMomentum_%s",frame[f].Data()),nc_cost[f],nc_mom[f]);
      FillHist(Form("n_pip_CosTvsPhi_%s",frame[f].Data()),nc_cost[f],nc_phi[f]);
      FillHist(Form("HeKn_pip_MMvsMomentum_%s",frame[f].Data()),nc_mmass[f],nc_mmom[f]);
      FillHist(Form("HeKn_pip_MMvsCosT_%s",frame[f].Data()),nc_mmass[f],nc_mcost[f]);
    }
  }
  if(particle->nPiminus()>0){
    for(int f=1; f<2; f++){
      FillHist(Form("n_pim_CosTvsMomentum_%s",frame[f].Data()),nc_cost[f],nc_mom[f]);
      FillHist(Form("n_pim_CosTvsPhi_%s",frame[f].Data()),nc_cost[f],nc_phi[f]);
      FillHist(Form("HeKn_pim_MMvsMomentum_%s",frame[f].Data()),nc_mmass[f],nc_mmom[f]);
      FillHist(Form("HeKn_pim_MMvsCosT_%s",frame[f].Data()),nc_mmass[f],nc_mcost[f]);
    }
  }
  if(particle->nKaon()>0){
    for(int f=1; f<2; f++){
      FillHist(Form("n_k_CosTvsMomentum_%s",frame[f].Data()),nc_cost[f],nc_mom[f]);
      FillHist(Form("n_k_CosTvsPhi_%s",frame[f].Data()),nc_cost[f],nc_phi[f]);
      FillHist(Form("HeKn_k_MMvsMomentum_%s",frame[f].Data()),nc_mmass[f],nc_mmom[f]);
      FillHist(Form("HeKn_k_MMvsCosT_%s",frame[f].Data()),nc_mmass[f],nc_mcost[f]);
    }
  }
  if(particle->nProton()>0){
    for(int f=1; f<2; f++){
      FillHist(Form("n_p_CosTvsMomentum_%s",frame[f].Data()),nc_cost[f],nc_mom[f]);
      FillHist(Form("n_p_CosTvsPhi_%s",frame[f].Data()),nc_cost[f],nc_phi[f]);
      FillHist(Form("HeKn_p_MMvsMomentum_%s",frame[f].Data()),nc_mmass[f],nc_mmom[f]);
      FillHist(Form("HeKn_p_MMvsCosT_%s",frame[f].Data()),nc_mmass[f],nc_mcost[f]);
    }
  }
  if(particle->nPiplus()==1&&particle->nPiminus()==1){
    for(int f=1; f<2; f++){
      FillHist(Form("n_pipi_CosTvsMomentum_%s",frame[f].Data()),nc_cost[f],nc_mom[f]);
      FillHist(Form("n_pipi_CosTvsPhi_%s",frame[f].Data()),nc_cost[f],nc_phi[f]);
      FillHist(Form("HeKn_pipi_MMvsMomentum_%s",frame[f].Data()),nc_mmass[f],nc_mmom[f]);
      FillHist(Form("HeKn_pipi_MMvsCosT_%s",frame[f].Data()),nc_mmass[f],nc_mcost[f]);
    }
  }
  if(k0flag){
    for(int f=1; f<2; f++){
      FillHist(Form("n_k0_CosTvsMomentum_%s",frame[f].Data()),nc_cost[f],nc_mom[f]);
      FillHist(Form("n_k0_CosTvsPhi_%s",frame[f].Data()),nc_cost[f],nc_phi[f]);
      FillHist(Form("HeKn_k0_MMvsMomentum_%s",frame[f].Data()),nc_mmass[f],nc_mmom[f]);
      FillHist(Form("HeKn_k0_MMvsCosT_%s",frame[f].Data()),nc_mmass[f],nc_mcost[f]);
    }
  }
  if(spflag){
    for(int f=1; f<2; f++){
      FillHist(Form("n_sp_CosTvsMomentum_%s",frame[f].Data()),nc_cost[f],nc_mom[f]);
      FillHist(Form("n_sp_CosTvsPhi_%s",frame[f].Data()),nc_cost[f],nc_phi[f]);
      FillHist(Form("HeKn_sp_MMvsMomentum_%s",frame[f].Data()),nc_mmass[f],nc_mmom[f]);
      FillHist(Form("HeKn_sp_MMvsCosT_%s",frame[f].Data()),nc_mmass[f],nc_mcost[f]);
      FillHist(Form("n_s_CosTvsMomentum_%s",frame[f].Data()),nc_cost[f],nc_mom[f]);
      FillHist(Form("n_s_CosTvsPhi_%s",frame[f].Data()),nc_cost[f],nc_phi[f]);
      FillHist(Form("HeKn_s_MMvsMomentum_%s",frame[f].Data()),nc_mmass[f],nc_mmom[f]);
      FillHist(Form("HeKn_s_MMvsCosT_%s",frame[f].Data()),nc_mmass[f],nc_mcost[f]);
      if(!k0flag){
        FillHist(Form("HeKn_sp_wok0_MMvsMomentum_%s",frame[f].Data()),nc_mmass[f],nc_mmom[f]);
        FillHist(Form("HeKn_sp_wok0_MMvsCosT_%s",frame[f].Data()),nc_mmass[f],nc_mcost[f]);
        FillHist(Form("HeKn_s_wok0_MMvsMomentum_%s",frame[f].Data()),nc_mmass[f],nc_mmom[f]);
        FillHist(Form("HeKn_s_wok0_MMvsCosT_%s",frame[f].Data()),nc_mmass[f],nc_mcost[f]);
      }
    }
  }
  if(smflag){
    for(int f=1; f<2; f++){
      FillHist(Form("n_sm_CosTvsMomentum_%s",frame[f].Data()),nc_cost[f],nc_mom[f]);
      FillHist(Form("n_sm_CosTvsPhi_%s",frame[f].Data()),nc_cost[f],nc_phi[f]);
      FillHist(Form("HeKn_sm_MMvsMomentum_%s",frame[f].Data()),nc_mmass[f],nc_mmom[f]);
      FillHist(Form("HeKn_sm_MMvsCosT_%s",frame[f].Data()),nc_mmass[f],nc_mcost[f]);
      FillHist(Form("n_s_CosTvsMomentum_%s",frame[f].Data()),nc_cost[f],nc_mom[f]);
      FillHist(Form("n_s_CosTvsPhi_%s",frame[f].Data()),nc_cost[f],nc_phi[f]);
      FillHist(Form("HeKn_s_MMvsMomentum_%s",frame[f].Data()),nc_mmass[f],nc_mmom[f]);
      FillHist(Form("HeKn_s_MMvsCosT_%s",frame[f].Data()),nc_mmass[f],nc_mcost[f]);
      if(!k0flag){
        FillHist(Form("HeKn_sm_wok0_MMvsMomentum_%s",frame[f].Data()),nc_mmass[f],nc_mmom[f]);
        FillHist(Form("HeKn_sm_wok0_MMvsCosT_%s",frame[f].Data()),nc_mmass[f],nc_mcost[f]);
        FillHist(Form("HeKn_s_wok0_MMvsMomentum_%s",frame[f].Data()),nc_mmass[f],nc_mmom[f]);
        FillHist(Form("HeKn_s_wok0_MMvsCosT_%s",frame[f].Data()),nc_mmass[f],nc_mcost[f]);
      }
    }
  }
  // Invariant mass
  for(int f=1; f<2; f++){
    FillHist(Form("n%s_IMvsMomentum_%s",pname[cds->pid()].Data(),frame[f].Data()),ncds_mass[f],ncds_mom[f]);
    FillHist(Form("n%s_IMvsCosT_%s",pname[cds->pid()].Data(),frame[f].Data()),ncds_mass[f],ncds_cost[f]);
    FillHist(Form("HeKn%s_MMvsMomentum_%s",pname[cds->pid()].Data(),frame[f].Data()),ncds_mmass[f],ncds_mmom[f]);
    FillHist(Form("HeKn%s_MMvsCosT_%s",pname[cds->pid()].Data(),frame[f].Data()),ncds_mmass[f],ncds_mcost[f]);
  }
  FillHist(Form("%s_DCA",pname[cds->pid()].Data()),cds->vdis());
  FillHist(Form("all_DCA"),cds->vdis());



  return true;

}

bool MyAnalysis3HeKn::FillHist(TString name, double val1, int weight)
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

bool MyAnalysis3HeKn::FillHist(TString name, TString val1, int weight)
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

bool MyAnalysis3HeKn::FillHist(TString name, double val1, double val2, int weight)
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

bool MyAnalysis3HeKn::FillHist(TString name, TString val1, TString val2, int weight)
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

void MyAnalysis3HeKn::CutCondition()
{
  /* K0 mass cut */ 
  double mean  = 0.4976;
  double sigma = 0.0069;
  k0ll = mean-2*sigma; k0ul = mean+2*sigma;
  sblk0ll = mean-5*sigma; sblk0ul = mean-3*sigma;
  sbuk0ll = mean+3*sigma; sbuk0ul = mean+5*sigma;
  /* Missing neutron mass cut */ 
  mean  = 0.9379;
  sigma = 0.0230;
  mnll = mean-2*sigma; mnul = mean+2*sigma;
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
}

bool MyAnalysis3HeKn::Initialize(ConfMan* confMan)
{
  std::cout << "### MyAnalysis3HeKn::Initialize ###" << std::endl;

  std::string ofname = confMan->GetOutFileName();
  ofname.insert(ofname.find(".root"),"_ana3HeKn");

  rtFile =  new TFile( Form("%s",ofname.c_str()), "RECREATE");
  rtFile -> cd();

  TString pname[8]  = {"pip","p","d","t","he","pim","k","o"};
  TString pname2[8] = {"#pi^{+}","p","d","t","he","#pi^{-}","K^{-}","o"};
  TString proname[6]  = {"pipi","k0","sp","sm","s","l"};
  TString proname2[6] = {"#pi^{+}#pi^{-}","K^{0}","#Sigma^{+}","#Sigma^{-}","#Sigma^{#pm}","#Lambda"};

  new TH1F( "EventNumber", "Number of Events", 10, 0, 10);
  TH1F* h1 = new TH1F( "TriggerPattern", "Trigger Pattern", 18, -0.5, 17.5);
  h1->GetXaxis()->SetBinLabel(1,"All");
  h1->GetXaxis()->SetBinLabel(2,"Beam");
  h1->GetXaxis()->SetBinLabel(3,"Kaon");
  h1->GetXaxis()->SetBinLabel(4,"KCDH1f");
  h1->GetXaxis()->SetBinLabel(5,"Pion");
  h1->GetXaxis()->SetBinLabel(6,"Proton");
  h1->GetXaxis()->SetBinLabel(7,"KCDH1");
  h1->GetXaxis()->SetBinLabel(8,"KCDH2");
  h1->GetXaxis()->SetBinLabel(9,"PivBVC");
  h1->GetXaxis()->SetBinLabel(10,"PiCDH1");
  h1->GetXaxis()->SetBinLabel(11,"PiCDH2");
  h1->GetXaxis()->SetBinLabel(12,"Kf");
  h1->GetXaxis()->SetBinLabel(13,"1stMix");
  h1->GetXaxis()->SetBinLabel(14,"Charged");
  h1->GetXaxis()->SetBinLabel(15,"Neutral");
  h1->GetXaxis()->SetBinLabel(16,"Cosmic");
  h1->GetXaxis()->SetBinLabel(17,"Reject");
  h1->GetXaxis()->SetBinLabel(18,"SIM");

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
  // DCA
  std::cout << "Define Histograms for DCA" << std::endl;
  new TH1F( "all_DCA", "DCA;DCA (cm);Coutns", 100, 0, 10);
  new TH1F( "pip_DCA", "DCA of #pi^{+};DCA (cm);Coutns", 100, 0, 10);
  new TH1F( "pim_DCA", "DCA of #pi^{-};DCA (cm);Coutns", 100, 0, 10);
  new TH1F( "k_DCA", "DCA of K^{-};DCA (cm);Coutns", 100, 0, 10);
  new TH1F( "p_DCA", "DCA of p;DCA (cm);Coutns", 100, 0, 10);
  // FWD neutral Particles
  std::cout << "Define Histograms for FWD neutral particles" << std::endl;
  new TH2F( "FWD_OverbetavsMomentum", "1/#beta vs. Momentum;1/#beta;Momentum (GeV/c)", 1000, -0, 10, 500, -5, 5 );
  new TH2F( "FWD_OverbetavsEnergy", "1/#beta vs. Energy deposit;1/#beta;Energy deposit (MeVee)", 1000, -0, 10, 5000, 0, 100);
  new TH2F( "FWD_TOFvsMomentum", "TOF vs. Momentum;TOF (ns);Momentum (GeV/c)", 4000, 30, 70, 500, -5, 5 );
  new TH2F( "FWD_HitPosition", "Hit position at NC (X vs. Y);X position (cm);Y position (cm)",40,-400,400,25,-250,250);
  new TH2F( "FWD_HitSegment", "NC segment;Segment in layer;Layer",16,0.5,16.5,7,0.5,7.5);
  new TH2F( "FWD_n_OverbetavsMomentum", "1/#beta vs. Momentum;1/#beta;Momentum (GeV/c)", 1000, -0, 10, 500, -5, 5 );
  new TH2F( "FWD_n_OverbetavsEnergy", "1/#beta vs. Energy deposit;1/#beta;Energy deposit (MeVee)", 1000, -0, 10, 5000, 0, 100);
  new TH2F( "FWD_n_TOFvsMomentum", "TOF vs. Momentum;TOF (ns);Momentum (GeV/c)", 4000, 30, 70, 500, -5, 5 );
  new TH2F( "FWD_n_HitPosition", "Hit position at NC (X vs. Y);X position (cm);Y position (cm)",40,-400,400,25,-250,250);
  new TH2F( "FWD_n_HitSegment", "NC segment;Segment in layer;Layer",16,0.5,16.5,7,0.5,7.5);
  new TH2F( "FWD_k0_OverbetavsMomentum", "1/#beta vs. Momentum;1/#beta;Momentum (GeV/c)", 1000, -0, 10, 500, -5, 5 );
  new TH2F( "FWD_k0_OverbetavsEnergy", "1/#beta vs. Energy deposit;1/#beta;Energy deposit (MeVee)", 1000, -0, 10, 5000, 0, 100);
  new TH2F( "FWD_k0_TOFvsMomentum", "TOF vs. Momentum;TOF (ns);Momentum (GeV/c)", 4000, 30, 70, 500, -5, 5 );
  new TH2F( "FWD_k0_HitPosition", "Hit position at NC (X vs. Y);X position (cm);Y position (cm)",40,-400,400,25,-250,250);
  new TH2F( "FWD_k0_HitSegment", "NC segment;Segment in layer;Layer",16,0.5,16.5,7,0.5,7.5);
  // neutron
  new TH2F("n_CosTvsMomentum_CM","cos#theta vs. momentum of n;cos#theta;Momentum (GeV/c)", 2000, -1.0, 1.0, 1000, 0.0, 1.0);
  new TH2F("n_CosTvsPhi_CM","cos#theta vs. #phi of n;cos#theta;#phi", 2000, -1.0, 1.0, 1600, 0.0, 3.2);
  for(int ip=0; ip<8; ip++){
    new TH2F(Form("n_%s_CosTvsMomentum_CM",pname[ip].Data()),"cos#theta vs. momentum of n;cos#theta;Momentum (GeV/c)", 2000, -1.0, 1.0, 1000, 0.0, 1.0);
    new TH2F(Form("n_%s_CosTvsPhi_CM",pname[ip].Data()),"cos#theta vs. #phi of n;cos#theta;#phi", 2000, -1.0, 1.0, 1600, 0.0, 3.2);
  }
  for(int ip=0; ip<6; ip++){
    new TH2F(Form("n_%s_CosTvsMomentum_CM",proname[ip].Data()),"cos#theta vs. momentum of n;cos#theta;Momentum (GeV/c)", 2000, -1.0, 1.0, 1000, 0.0, 1.0);
    new TH2F(Form("n_%s_CosTvsPhi_CM",proname[ip].Data()),"cos#theta vs. #phi of n;cos#theta;#phi", 2000, -1.0, 1.0, 1600, 0.0, 3.2);
  }
  // Invariant mass
  std::cout << "Define Histograms for IM(nX) M.M." << std::endl;
  for(int ip=0; ip<8; ip++){
    new TH2F(Form("n%s_IMvsMomentum_CM",pname[ip].Data()), Form("Invariant mass (n%s) vs Momentum;IM(n%s) (GeV/c^{2});Momentum (GeV/c)",pname2[ip].Data(),pname2[ip].Data()), 2000, 0.0, 2.0, 1000, 0.0, 1.0 );
    new TH2F(Form("n%s_IMvsCosT_CM",pname[ip].Data()), Form("Invariant mass (n%s) vs cos#theta;IM(n%s) (GeV/c^{2});cos#theta",pname2[ip].Data(),pname2[ip].Data()), 2000, 0.0, 2.0, 2000, -1.0, 1.0 );
  }
  new TH2F(Form("pipi_n_IMvsMomentum_CM"), Form("Invariant mass (#pi^{+}#pi^{-}) vs Momentum;IM(#pi^{+}#pi^{-}) (GeV/c^{2});Momentum (GeV/c)"), 2000, 0.0, 2.0, 1000, 0.0, 2.0 );
  new TH2F(Form("npip_npip_IMvsMomentum_CM"), Form("Invariant mass (n#pi^{+}) vs Momentum;IM(n#pi^{+}) (GeV/c^{2});Momentum (GeV/c)"), 2000, 0.0, 2.0, 1000, 0.0, 2.0 );
  new TH2F(Form("npim_npim_IMvsMomentum_CM"), Form("Invariant mass (n#pi^{-}) vs Momentum;IM(n#pi^{-}) (GeV/c^{2});Momentum (GeV/c)"), 2000, 0.0, 2.0, 1000, 0.0, 2.0 );
  // Missing mass
  std::cout << "Define Histograms for p(K-,n)X M.M." << std::endl;
  new TH2F("HeKn_MMvsMomentum_CM", "p(K^{-},n)X missing mass vs. missing momentum;MM_{p}(n) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 1.0, 3.0, 1000, 0.0, 2.0);
  new TH2F("HeKn_MMvsCosT_CM", "p(K^{-},n)X missing mass vs. cos#theta;MM_{p}(n) (GeV/c^{2});cos#theta", 2000, 1.0, 3.0, 2000, -1.0, 1.0);
  for(int ip=0; ip<8; ip++){
    new TH2F(Form("HeKn_%s_MMvsMomentum_CM",pname[ip].Data()), "p(K^{-},n)X missing mass vs. missing momentum;MM_{p}(n) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 1.0, 3.0, 1000, 0.0, 2.0);
    new TH2F(Form("HeKn_%s_MMvsCosT_CM",pname[ip].Data()), "p(K^{-},n)X missing mass vs. cos#theta;MM_{p}(n) (GeV/c^{2});cos#theta", 2000, 1.0, 3.0, 2000, -1.0, 1.0);
    new TH2F(Form("HeKn%s_MMvsMomentum_CM",pname[ip].Data()), Form("p(K^{-},n%s)X missing mass vs. missing momentum;MM_{p}(n%s) (GeV/c^{2});Missing momentum (GeV/c)",pname2[ip].Data(),pname2[ip].Data()), 2000, 1.0, 3.0, 1000, 0.0, 2.0);
    new TH2F(Form("HeKn%s_MMvsCosT_CM",pname[ip].Data()), Form("p(K^{-},n%s)X missing mass vs. cos#theta;MM_{p}(n%s) (GeV/c^{2});cos#theta",pname2[ip].Data(),pname2[ip].Data()), 2000, 1.0, 3.0, 2000, -1.0, 1.0);
  }
  for(int ip=0; ip<6; ip++){
    new TH2F(Form("HeKn_%s_MMvsMomentum_CM",proname[ip].Data()), "p(K^{-},n)X missing mass vs. missing momentum;MM_{p}(n) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 1.0, 3.0, 1000, 0.0, 2.0);
    new TH2F(Form("HeKn_%s_MMvsCosT_CM",proname[ip].Data()), "p(K^{-},n)X missing mass vs. cos#theta;MM_{p}(n) (GeV/c^{2});cos#theta", 2000, 1.0, 3.0, 2000, -1.0, 1.0);
    if(2<=ip&&ip<=4){
      new TH2F(Form("HeKn_%s_wok0_MMvsMomentum_CM",proname[ip].Data()), "p(K^{-},n)X missing mass vs. missing momentum;MM_{p}(n) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 1.0, 3.0, 1000, 0.0, 2.0);
      new TH2F(Form("HeKn_%s_wok0_MMvsCosT_CM",proname[ip].Data()), "p(K^{-},n)X missing mass vs. cos#theta;MM_{p}(n) (GeV/c^{2});cos#theta", 2000, 1.0, 3.0, 2000, -1.0, 1.0);
    }
    new TH2F(Form("HeKn%s_MMvsMomentum_CM",proname[ip].Data()), Form("p(K^{-},n%s)X missing mass vs. missing momentum;MM_{p}(n%s) (GeV/c^{2});Missing momentum (GeV/c)",proname2[ip].Data(),proname2[ip].Data()), 2000, 1.0, 3.0, 1000, 0.0, 2.0);
    new TH2F(Form("HeKn%s_MMvsCosT_CM",proname[ip].Data()), Form("p(K^{-},n%s)X missing mass vs. cos#theta;MM_{p}(n%s) (GeV/c^{2});cos#theta",proname2[ip].Data(),proname2[ip].Data()), 2000, 1.0, 3.0, 2000, -1.0, 1.0);
  }

  return true;
}
