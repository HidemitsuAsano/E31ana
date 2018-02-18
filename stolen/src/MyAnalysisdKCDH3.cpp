// MyAnalysisdKCDH3.cpp

#include "MyAnalysisdKCDH3.h"

MyAnalysisdKCDH3::MyAnalysisdKCDH3(TFile* rt, ConfMan* conf)
{
  Initialize(conf);
  CutCondition();
  Clear();
}

MyAnalysisdKCDH3::~MyAnalysisdKCDH3()
{
  Clear();
  rtFile->cd();
  rtFile->Write();
  rtFile->Close();
}

void MyAnalysisdKCDH3::Clear()
{
}

bool MyAnalysisdKCDH3::DoAnalysis(ConfMan* conf, EventHeader* header, BeamLineHitMan* blMan, BeamLineTrackMan* bltrackMan, CDSHitMan* cdsMan, CDSTrackingMan* cdstrackMan, Particle* particle)
{
  rtFile->cd();

  FillHist("EventNumber",0); /* All Events */
  DetectorList *dlist=DetectorList::GetInstance();
  if(particle->nBeam()!=1) return false;
  pBeam* beam = particle->beam(0);
  double beamtof = beam->bhdt0tof();

  /* Event selection */
  if(!header->IsTrig(Trig_KCDH3)) return true;
  if(particle->nCDS()!=3) return true;
  FillHist("EventNumber",1); /* 2 hit in CDS*/
  bool neutral=false; bool charged=false;
  if(header->IsTrig(Trig_Charged)&&particle->nPC()>0) charged = true;
  if(neutral&&charged) return true;
  charged = true;
  FillHist("EventNumber",2); /* Both Neutral and Charged */
  if(charged&&(particle->nPiminus()!=2||particle->nProton()!=1)) return true;
  FillHist("EventNumber",3); /* 2 pi-  and proton in Charged */

  // ############### //
  // Vertex decision //
  // ############### //
  TVector3 vertex;
  double vdis = 9999;
  int fcds = -1;
  for(int it=0; it<particle->nCDS(); it++){
    pCDS* cds = particle->cdsi(it);
    if(cds->chi()>30) continue;
    if(vdis>9998||vdis>cds->vdis()){
      vdis = cds->vdis();
      vertex = cds->vbeam();
      if(GeomTools::GetID(vertex)==CID_Fiducial) fcds=it;
    }
  }
  if(fcds==-1) return true; /* Vertex decision */
  pCDS* cds=0;
  if(fcds!=-1) cds = particle->cdsi(fcds);

  if(neutral) DoNeutralAna(particle,beam,cds,vertex);
  if(charged) DoChargedAna(particle,beam,cds,vertex);

  return true;

}

bool MyAnalysisdKCDH3::DoNeutralAna(Particle* particle, pBeam* beam, pCDS* cds, TVector3 vertex)
{
  return true;
}

bool MyAnalysisdKCDH3::DoChargedAna(Particle* particle, pBeam* beam, pCDS* cds, TVector3 vertex)
{
  DetectorList *dlist=DetectorList::GetInstance();
  TString mode = "Charged";

  pCDS* pim1 = particle->pim(0);
  pCDS* pim2 = particle->pim(1);
  pCDS* p    = particle->proton(0);
  if(pim1->chi()>30) return true;
  if(pim2->chi()>30) return true;
  if(p->chi()>30) return true;

  // ######################### //
  // 2 Helix calculation       //
  // ######################### //
  //pCDS* pip1 = 0;
  //pCDS* pip2 = 0;
  //{
  //  for(int it=0; it<particle->nProduct(); it++){
  //    pCDS* product = particle->product(it);
  //    int comb = product->comb();
  //    if(comb==pow(2,CDS_Proton)+pow(2,CDS_PiMinus)){
  //      if(product->daughter1()==pim1->id()||product->daughter2()==pim1->id()){
  //        pip1 = product;
  //      }
  //      else if(product->daughter1()==pim2->id()||product->daughter2()==pim2->id()){
  //        pip2 = product;
  //      }
  //    }
  //  }
  //}
  //if(pip1==0||pip2==0){ return true; }
  //FillHist(Form("pip1_IMvspip2_IM_All"),pip1->mass(),pip2->mass());

  //bool pip1flag = false, pip2flag = false;
  //pCDS* lam = 0;
  //pCDS* proton = 0;
  //if(lamll<pip1->mass()&&pip1->mass()<lamul){
  //  pip1flag = true;
  //  lam = pip1;
  //  proton = p2;
  //}
  //if(lamll<pip2->mass()&&pip2->mass()<lamul){
  //  pip2flag = true;
  //  pip1flag = false;
  //  lam = pip2;
  //  proton = p1;
  //}
  //if(!pip1flag&&!pip2flag){ return true; }
  //if(lam==0){ return true; }
  //if(pip1flag){
  //  vertex = p2->vbeam();
  //  vdis = p2->vbdis();
  //  if(GeomTools::GetID(vertex)!=CID_Fiducial){ return true; }
  //}
  //if(pip2flag){
  //  vertex = p1->vbeam();
  //  vdis = p1->vbdis();
  //  if(GeomTools::GetID(vertex)!=CID_Fiducial){ return true; }
  //}

  //TString frame[2] = {"Lab","CM"};

  //// Target
  //TVector3 ZeroV;
  //TLorentzVector Ltgt;
  //if(dlist->GetMaterial(CID_Target)=="LHydrogen"){
  //  Ltgt.SetVectM(ZeroV,pMass);
  //}
  //if(dlist->GetMaterial(CID_Target)=="LDeuterium"){
  //  Ltgt.SetVectM(ZeroV,dMass);
  //}
  //if(dlist->GetMaterial(CID_Target)=="LHelium-3"){
  //  Ltgt.SetVectM(ZeroV,ThreeHeMass);
  //}
  //TLorentzVector Lbeam = beam->GetLorentzVector(vertex);
  //TVector3 boost = (Ltgt+Lbeam).BoostVector();
  //TLorentzVector cLtgt = Ltgt; cLtgt.Boost(-boost); TLorentzVector cLbeam = Lbeam; cLbeam.Boost(-boost);
  ///* Pi Minus-1 */
  //TLorentzVector Lpim1 = pim1->GetLorentzVector();
  //TLorentzVector cLpim1 = pim1->GetLorentzVector(); cLpim1.Boost(-boost);
  //double pim1_mass[2]  = { (Lpim1).M()                    , (cLpim1).M()};
  //double pim1_mom[2]   = { (Lpim1).Vect().Mag()           , (cLpim1).Vect().Mag()};
  //double pim1_cost[2]  = { (Lpim1).Vect().CosTheta()      , (cLpim1).Vect().CosTheta()};
  //double pim1_phi[2]   = { (Lpim1).Vect().Phi()                  , (cLpim1).Vect().Phi()};
  //double pim1_mmass[2]   = { (Ltgt+Lbeam-Lpim1).M()              , (cLtgt+cLbeam-cLpim1).M()};
  //double pim1_mmom[2]    = { (Ltgt+Lbeam-Lpim1).Vect().Mag()     , (cLtgt+cLbeam-cLpim1).Vect().Mag()};
  //double pim1_mcost[2]   = { (Ltgt+Lbeam-Lpim1).Vect().CosTheta(), (cLtgt+cLbeam-cLpim1).Vect().CosTheta()};
  ///* Pi Minus-2 */
  //TLorentzVector Lpim2 = pim2->GetLorentzVector();
  //TLorentzVector cLpim2 = pim2->GetLorentzVector(); cLpim2.Boost(-boost);
  //double pim2_mass[2]  = { (Lpim2).M()                    , (cLpim2).M()};
  //double pim2_mom[2]   = { (Lpim2).Vect().Mag()           , (cLpim2).Vect().Mag()};
  //double pim2_cost[2]  = { (Lpim2).Vect().CosTheta()      , (cLpim2).Vect().CosTheta()};
  //double pim2_phi[2]   = { (Lpim2).Vect().Phi()                  , (cLpim2).Vect().Phi()};
  //double pim2_mmass[2]   = { (Ltgt+Lbeam-Lpim2).M()              , (cLtgt+cLbeam-cLpim2).M()};
  //double pim2_mmom[2]    = { (Ltgt+Lbeam-Lpim2).Vect().Mag()     , (cLtgt+cLbeam-cLpim2).Vect().Mag()};
  //double pim2_mcost[2]   = { (Ltgt+Lbeam-Lpim2).Vect().CosTheta(), (cLtgt+cLbeam-cLpim2).Vect().CosTheta()};
  ///* pi + pi */
  //TLorentzVector Lpipi = Lpim1 + Lpim2;
  //TLorentzVector cLpipi = cLpim1 + cLpim2;
  //double pipi_mass[2]  = { (Lpipi).M()                    , (cLpipi).M()};
  //double pipi_mom[2]   = { (Lpipi).Vect().Mag()           , (cLpipi).Vect().Mag()};
  //double pipi_cost[2]  = { (Lpipi).Vect().CosTheta()      , (cLpipi).Vect().CosTheta()};
  //double pipi_phi[2]   = { (Lpipi).Vect().Phi()                  , (cLpipi).Vect().Phi()};
  //double pipi_mmass[2]   = { (Ltgt+Lbeam-Lpipi).M()              , (cLtgt+cLbeam-cLpipi).M()};
  //double pipi_mmom[2]    = { (Ltgt+Lbeam-Lpipi).Vect().Mag()     , (cLtgt+cLbeam-cLpipi).Vect().Mag()};
  //double pipi_mcost[2]   = { (Ltgt+Lbeam-Lpipi).Vect().CosTheta(), (cLtgt+cLbeam-cLpipi).Vect().CosTheta()};
  ///* P */
  //TLorentzVector Lp  = p ->GetLorentzVector();
  //TLorentzVector cLp  = p ->GetLorentzVector(); cLp.Boost(-boost);
  //double p_mass[2]    = { Lp.M()                           , cLp.M()};
  //double p_mom[2]     = { Lp.Vect().Mag()                  , cLp.Vect().Mag()};
  //double p_cost[2]    = { Lp.Vect().CosTheta()             , cLp.Vect().CosTheta()};
  //double p_phi[2]     = { Lp.Vect().Phi()                  , cLp.Vect().Phi()};
  //double p_mmass[2]   = { (Ltgt+Lbeam-Lp).M()              , (cLtgt+cLbeam-cLp).M()};
  //double p_mmom[2]    = { (Ltgt+Lbeam-Lp).Vect().Mag()     , (cLtgt+cLbeam-cLp).Vect().Mag()};
  //double p_mcost[2]   = { (Ltgt+Lbeam-Lp).Vect().CosTheta(), (cLtgt+cLbeam-cLp).Vect().CosTheta()};
  ///* pi1 + PC */
  //double ppim1_mass[2]  = { (Lpim1+Lp).M()                           , (cLpim1+cLp).M()};
  //double ppim1_mom[2]   = { (Lpim1+Lp).Vect().Mag()                  , (cLpim1+cLp).Vect().Mag()};
  //double ppim1_cost[2]  = { (Lpim1+Lp).Vect().CosTheta()             , (cLpim1+cLp).Vect().CosTheta()};
  //double ppim1_phi[2]   = { (Lpim1+Lp).Vect().Phi()                  , (cLpim1+cLp).Vect().Phi()};
  //double ppim1_mmass[2] = { (Ltgt+Lbeam-(Lpim1+Lp)).M()              , (cLtgt+cLbeam-(cLpim1+cLp)).M()};
  //double ppim1_mmom[2]  = { (Ltgt+Lbeam-(Lpim1+Lp)).Vect().Mag()     , (cLtgt+cLbeam-(cLpim1+cLp)).Vect().Mag()};
  //double ppim1_mcost[2] = { (Ltgt+Lbeam-(Lpim1+Lp)).Vect().CosTheta(), (cLtgt+cLbeam-(cLpim1+cLp)).Vect().CosTheta()};
  ///* pi2 + PC */
  //double ppim2_mass[2]  = { (Lpim2+Lp).M()                           , (cLpim2+cLp).M()};
  //double ppim2_mom[2]   = { (Lpim2+Lp).Vect().Mag()                  , (cLpim2+cLp).Vect().Mag()};
  //double ppim2_cost[2]  = { (Lpim2+Lp).Vect().CosTheta()             , (cLpim2+cLp).Vect().CosTheta()};
  //double ppim2_phi[2]   = { (Lpim2+Lp).Vect().Phi()                  , (cLpim2+cLp).Vect().Phi()};
  //double ppim2_mmass[2] = { (Ltgt+Lbeam-(Lpim2+Lp)).M()              , (cLtgt+cLbeam-(cLpim2+cLp)).M()};
  //double ppim2_mmom[2]  = { (Ltgt+Lbeam-(Lpim2+Lp)).Vect().Mag()     , (cLtgt+cLbeam-(cLpim2+cLp)).Vect().Mag()};
  //double ppim2_mcost[2] = { (Ltgt+Lbeam-(Lpim2+Lp)).Vect().CosTheta(), (cLtgt+cLbeam-(cLpim2+cLp)).Vect().CosTheta()};
  ///* pi1 + pi2 + PC */
  //double ppipi_mass[2]  = { (Lpipi+Lp).M()                           , (cLpipi+cLp).M()};
  //double ppipi_mom[2]   = { (Lpipi+Lp).Vect().Mag()                  , (cLpipi+cLp).Vect().Mag()};
  //double ppipi_cost[2]  = { (Lpipi+Lp).Vect().CosTheta()             , (cLpipi+cLp).Vect().CosTheta()};
  //double ppipi_phi[2]   = { (Lpipi+Lp).Vect().Phi()                  , (cLpipi+cLp).Vect().Phi()};
  //double ppipi_mmass[2] = { (Ltgt+Lbeam-(Lpipi+Lp)).M()              , (cLtgt+cLbeam-(cLpipi+cLp)).M()};
  //double ppipi_mmom[2]  = { (Ltgt+Lbeam-(Lpipi+Lp)).Vect().Mag()     , (cLtgt+cLbeam-(cLpipi+cLp)).Vect().Mag()};
  //double ppipi_mcost[2] = { (Ltgt+Lbeam-(Lpipi+Lp)).Vect().CosTheta(), (cLtgt+cLbeam-(cLpipi+cLp)).Vect().CosTheta()};

  ///* flag */
  //bool flam1flag = false; /* L1 flag */
  //bool sblflam1flag = false; /* lower sideband of L1 flag */
  //bool sbuflam1flag = false; /* upper sideband of L1 flag */
  //{
  //  double pro_mmass = ppim1_mass[1];
  //  if(flamll<pro_mmass&&pro_mmass<flamul) flam1flag = true;
  //  if(sblflamll<pro_mmass&&pro_mmass<sblflamul) sblflam1flag = true;
  //  if(sbuflamll<pro_mmass&&pro_mmass<sbuflamul) sbuflam1flag = true;
  //}
  //bool flam2flag = false; /* L2 flag */
  //bool sblflam2flag = false; /* lower sideband of L2 flag */
  //bool sbuflam2flag = false; /* upper sideband of L2 flag */
  //{
  //  double pro_mmass = ppim2_mass[1];
  //  if(flamll<pro_mmass&&pro_mmass<flamul) flam2flag = true;
  //  if(sblflamll<pro_mmass&&pro_mmass<sblflamul) sblflam2flag = true;
  //  if(sbuflamll<pro_mmass&&pro_mmass<sbuflamul) sbuflam2flag = true;
  //}
  //bool mpflag = false; /* Missing-p flag */
  //bool sblmpflag = false; /* lower sideband of Missing-p flag */
  //bool sbumpflag = false; /* upper sideband of Missing-p flag */
  //{
  //  double pro_mmass = ppipi_mmass[1];
  //  if(mpll<pro_mmass&&pro_mmass<mpul) mpflag = true;
  //  if(sblmpll<pro_mmass&&pro_mmass<sblmpul) sblmpflag = true;
  //  if(sbumpll<pro_mmass&&pro_mmass<sbumpul) sbumpflag = true;
  //}

  ////if(!mpflag) {return true;}

  //// pi1
  //for(int f=1; f<2; f++){
  //  FillHist(Form("pim1_Mass_%s_%s",frame[f].Data(),mode.Data()),pim1_mass[f]);
  //  FillHist(Form("pim1_Momentum_%s_%s",frame[f].Data(),mode.Data()),pim1_mom[f]);
  //  FillHist(Form("pim1_CosT_%s_%s",frame[f].Data(),mode.Data()),pim1_cost[f]);
  //  FillHist(Form("pim1_Phi_%s_%s",frame[f].Data(),mode.Data()),pim1_phi[f]);
  //  FillHist(Form("pim1_MMass_%s_%s",frame[f].Data(),mode.Data()),pim1_mmass[f]);
  //  FillHist(Form("pim1_MMomentum_%s_%s",frame[f].Data(),mode.Data()),pim1_mmom[f]);
  //  FillHist(Form("pim1_MCosT_%s_%s",frame[f].Data(),mode.Data()),pim1_mcost[f]);
  //}
  //// pi2
  //for(int f=1; f<2; f++){
  //  FillHist(Form("pim2_Mass_%s_%s",frame[f].Data(),mode.Data()),pim2_mass[f]);
  //  FillHist(Form("pim2_Momentum_%s_%s",frame[f].Data(),mode.Data()),pim2_mom[f]);
  //  FillHist(Form("pim2_CosT_%s_%s",frame[f].Data(),mode.Data()),pim2_cost[f]);
  //  FillHist(Form("pim2_Phi_%s_%s",frame[f].Data(),mode.Data()),pim2_phi[f]);
  //  FillHist(Form("pim2_MMass_%s_%s",frame[f].Data(),mode.Data()),pim2_mmass[f]);
  //  FillHist(Form("pim2_MMomentum_%s_%s",frame[f].Data(),mode.Data()),pim2_mmom[f]);
  //  FillHist(Form("pim2_MCosT_%s_%s",frame[f].Data(),mode.Data()),pim2_mcost[f]);
  //}
  //// pi1 + pi2
  //for(int f=1; f<2; f++){
  //  FillHist(Form("pipi_Mass_%s_%s",frame[f].Data(),mode.Data()),pipi_mass[f]);
  //  FillHist(Form("pipi_Momentum_%s_%s",frame[f].Data(),mode.Data()),pipi_mom[f]);
  //  FillHist(Form("pipi_CosT_%s_%s",frame[f].Data(),mode.Data()),pipi_cost[f]);
  //  FillHist(Form("pipi_Phi_%s_%s",frame[f].Data(),mode.Data()),pipi_phi[f]);
  //  FillHist(Form("pipi_MMass_%s_%s",frame[f].Data(),mode.Data()),pipi_mmass[f]);
  //  FillHist(Form("pipi_MMomentum_%s_%s",frame[f].Data(),mode.Data()),pipi_mmom[f]);
  //  FillHist(Form("pipi_MCosT_%s_%s",frame[f].Data(),mode.Data()),pipi_mcost[f]);
  //}
  //// p
  //for(int f=1; f<2; f++){
  //  FillHist(Form("p_Mass_%s_%s",frame[f].Data(),mode.Data()),p_mass[f]);
  //  FillHist(Form("p_Momentum_%s_%s",frame[f].Data(),mode.Data()),p_mom[f]);
  //  FillHist(Form("p_CosT_%s_%s",frame[f].Data(),mode.Data()),p_cost[f]);
  //  FillHist(Form("p_Phi_%s_%s",frame[f].Data(),mode.Data()),p_phi[f]);
  //  FillHist(Form("p_MMass_%s_%s",frame[f].Data(),mode.Data()),p_mmass[f]);
  //  FillHist(Form("p_MMomentum_%s_%s",frame[f].Data(),mode.Data()),p_mmom[f]);
  //  FillHist(Form("p_MCosT_%s_%s",frame[f].Data(),mode.Data()),p_mcost[f]);
  //}
  //// pi1 + p
  //for(int f=1; f<2; f++){
  //  FillHist(Form("ppim1_Mass_%s_%s",frame[f].Data(),mode.Data()),ppim1_mass[f]);
  //  FillHist(Form("ppim1_Momentum_%s_%s",frame[f].Data(),mode.Data()),ppim1_mom[f]);
  //  FillHist(Form("ppim1_CosT_%s_%s",frame[f].Data(),mode.Data()),ppim1_cost[f]);
  //  FillHist(Form("ppim1_Phi_%s_%s",frame[f].Data(),mode.Data()),ppim1_phi[f]);
  //  FillHist(Form("ppim1_MMass_%s_%s",frame[f].Data(),mode.Data()),ppim1_mmass[f]);
  //  FillHist(Form("ppim1_MMomentum_%s_%s",frame[f].Data(),mode.Data()),ppim1_mmom[f]);
  //  FillHist(Form("ppim1_MCosT_%s_%s",frame[f].Data(),mode.Data()),ppim1_mcost[f]);
  //}
  //// pi2 + p
  //for(int f=1; f<2; f++){
  //  FillHist(Form("ppim2_Mass_%s_%s",frame[f].Data(),mode.Data()),ppim2_mass[f]);
  //  FillHist(Form("ppim2_Momentum_%s_%s",frame[f].Data(),mode.Data()),ppim2_mom[f]);
  //  FillHist(Form("ppim2_CosT_%s_%s",frame[f].Data(),mode.Data()),ppim2_cost[f]);
  //  FillHist(Form("ppim2_Phi_%s_%s",frame[f].Data(),mode.Data()),ppim2_phi[f]);
  //  FillHist(Form("ppim2_MMass_%s_%s",frame[f].Data(),mode.Data()),ppim2_mmass[f]);
  //  FillHist(Form("ppim2_MMomentum_%s_%s",frame[f].Data(),mode.Data()),ppim2_mmom[f]);
  //  FillHist(Form("ppim2_MCosT_%s_%s",frame[f].Data(),mode.Data()),ppim2_mcost[f]);
  //}
  //// pi1 + pi2 + p
  //for(int f=1; f<2; f++){
  //  FillHist(Form("ppipi_Mass_%s_%s",frame[f].Data(),mode.Data()),ppipi_mass[f]);
  //  FillHist(Form("ppipi_Momentum_%s_%s",frame[f].Data(),mode.Data()),ppipi_mom[f]);
  //  FillHist(Form("ppipi_CosT_%s_%s",frame[f].Data(),mode.Data()),ppipi_cost[f]);
  //  FillHist(Form("ppipi_Phi_%s_%s",frame[f].Data(),mode.Data()),ppipi_phi[f]);
  //  FillHist(Form("ppipi_MMass_%s_%s",frame[f].Data(),mode.Data()),ppipi_mmass[f]);
  //  FillHist(Form("ppipi_MMomentum_%s_%s",frame[f].Data(),mode.Data()),ppipi_mmom[f]);
  //  FillHist(Form("ppipi_MCosT_%s_%s",frame[f].Data(),mode.Data()),ppipi_mcost[f]);
  //}
  //// 2D Plots 
  //for(int f=1; f<2; f++){
  //  FillHist(Form("ppim1_Mass_vsppim1_MMass_%s_%s",frame[f].Data(),mode.Data()),ppim1_mass[f],ppim1_mmass[f]);
  //  FillHist(Form("ppim2_Mass_vsppim2_MMass_%s_%s",frame[f].Data(),mode.Data()),ppim2_mass[f],ppim2_mmass[f]);
  //  FillHist(Form("ppim1_Mass_vsppim2_Mass_%s_%s",frame[f].Data(),mode.Data()),ppim1_mass[f],ppim2_mass[f]);
  //  FillHist(Form("ppim1_MMass_vsppim2_MMass_%s_%s",frame[f].Data(),mode.Data()),ppim1_mmass[f],ppim2_mmass[f]);
  //  FillHist(Form("p_MMass_vsppipi_MMass_%s_%s",frame[f].Data(),mode.Data()),p_mmass[f],ppipi_mmass[f]);
  //  FillHist(Form("p_MMass_vsppim1_Mass_%s_%s",frame[f].Data(),mode.Data()),p_mmass[f],ppim1_mass[f]);
  //  FillHist(Form("p_MMass_vsppim1_MMass_%s_%s",frame[f].Data(),mode.Data()),p_mmass[f],ppim1_mmass[f]);
  //  FillHist(Form("p_MMass_vsppim2_Mass_%s_%s",frame[f].Data(),mode.Data()),p_mmass[f],ppim2_mass[f]);
  //  FillHist(Form("p_MMass_vsppim2_MMass_%s_%s",frame[f].Data(),mode.Data()),p_mmass[f],ppim2_mmass[f]);
  //}

  return true;

}

bool MyAnalysisdKCDH3::FillHist(TString name, double val1, int weight)
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

bool MyAnalysisdKCDH3::FillHist(TString name, TString val1, int weight)
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

bool MyAnalysisdKCDH3::FillHist(TString name, double val1, double val2, int weight)
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

bool MyAnalysisdKCDH3::FillHist(TString name, TString val1, TString val2, int weight)
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

void MyAnalysisdKCDH3::CutCondition()
{
  /* K0 mass cut */ 
  double mean  = 0.49732;
  double sigma = 0.00607;
  k0ll = mean-2*sigma; k0ul = mean+2*sigma;
  sblk0ll = mean-5*sigma; sblk0ul = mean-3*sigma;
  sbuk0ll = mean+3*sigma; sbuk0ul = mean+5*sigma;
  /* Missing neutron mass cut */ 
  mean  = 0.9432;
  sigma = 0.01432;
  mnll = 0.90; mnul = 0.98;
  //mnll = mean-2*sigma; mnul = mean+2*sigma;
  sblmnll = mean-5*sigma; sblmnul = mean-3*sigma;
  sbumnll = mean+3*sigma; sbumnul = mean+5*sigma;
  /* Missing protron mass cut */ 
  mean  = 0.94909;
  sigma = 0.01359;
  //mpll = 0.90; mpul = 0.98;
  mpll = 0.99; mpul = 1.07;
  //mpll = mean-2*sigma; mpul = mean+2*sigma;
  sblmpll = mean-5*sigma; sblmpul = mean-3*sigma;
  sbumpll = mean+3*sigma; sbumpul = mean+5*sigma;
  /* SigmaPlus mass cut */
  mean  = 1.18894;
  sigma = 0.00429;
  spll = mean-2*sigma; spul = mean+2*sigma;
  sblspll = mean-5*sigma; sblspul = mean-3*sigma;
  sbuspll = mean+3*sigma; sbuspul = mean+5*sigma;
  /* SigmaMinus mass cut */
  mean  = 1.19736;
  sigma = 0.00353;
  smll = mean-2*sigma; smul = mean+2*sigma;
  sblsmll = mean-5*sigma; sblsmul = mean-3*sigma;
  sbusmll = mean+3*sigma; sbusmul = mean+5*sigma;
  /* Lambda mass cut */ 
  mean  = 1.1155;
  sigma = 0.0020;
  lamll = mean-2*sigma; lamul = mean+2*sigma;
  sbllamll = mean-5*sigma; sbllamul = mean-3*sigma;
  sbulamll = mean+3*sigma; sbulamul = mean+5*sigma;
  /* Forward Lambda mass cut */ 
  mean  = 1.11599;
  sigma = 0.00260;
  flamll = mean-2*sigma; flamul = mean+2*sigma;
  sblflamll = mean-5*sigma; sblflamul = mean-3*sigma;
  sbuflamll = mean+3*sigma; sbuflamul = mean+5*sigma;
}

bool MyAnalysisdKCDH3::Initialize(ConfMan* confMan)
{
  std::cout << "### MyAnalysisdKCDH3::Initialize ###" << std::endl;

  std::string ofname = confMan->GetOutFileName();
  ofname.insert(ofname.find(".root"),"_anadKCDH3");

  rtFile =  new TFile( Form("%s",ofname.c_str()), "RECREATE");
  rtFile -> cd();

  new TH1F( "EventNumber", "Number of Events", 20, 0, 20);
  // FWD neutral Particles
  std::cout << "Define Histograms for FWD neutral particles" << std::endl;
  new TH1F( "FWDN_Overbeta", "1/#beta;1/#beta;Counts", 5000, 0, 5);
  new TH1F( "FWDN_Overbeta_withdECut", "1/#beta;1/#beta;Counts", 5000, 0, 5);
  new TH1F( "FWDN_Overbeta_Selected", "1/#beta;1/#beta;Counts", 5000, 0, 5);
  new TH2F( "FWDN_OverbetavsEnergy", "1/#beta vs. Energy deposit;1/#beta;Energy deposit (MeVee)", 1000, -0, 10, 5000, 0, 100);
  new TH2F( "FWDN_OverbetavsEnergy_withdECut", "1/#beta vs. Energy deposit;1/#beta;Energy deposit (MeVee)", 1000, -0, 10, 5000, 0, 100);
  new TH2F( "NCHitPosition_XY", "Hit position at NC (X vs. Y);X position (cm);Y position (cm)",16,-160,160,15,-75,75);
  new TH2F( "NCHitSegment", "NC segment;Segment in layer;Layer",16,0.5,16.5,7,0.5,7.5);
  // pi+
  new TH1F("pip_Mass_CM_Neutral","Invariant mass of #pi^{+};IM(#pi^{+}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pip_Momentum_CM_Neutral","Momentum of #pi^{+};Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pip_CosT_CM_Neutral","cos#theta of #pi^{+};cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("pip_Phi_CM_Neutral","#phi of #pi^{+};#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("pip_MMass_CM_Neutral", "^{3}He(K^{-},#pi^{+})X missing mass;MM(#pi^{+}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pip_MMomentum_CM_Neutral", "^{3}He(K^{-},#pi^{+})X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pip_MCosT_CM_Neutral", "^{3}He(K^{-},#pi^{+})X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // pi-
  new TH1F("pim_Mass_CM_Neutral","Invariant mass of #pi^{-};IM(#pi^{-}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pim_Momentum_CM_Neutral","Momentum of #pi^{-};Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pim_CosT_CM_Neutral","cos#theta of #pi^{-};cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("pim_Phi_CM_Neutral","#phi of #pi^{-};#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("pim_MMass_CM_Neutral", "^{3}He(K^{-},#pi^{-})X missing mass;MM(#pi^{-}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pim_MMomentum_CM_Neutral", "^{3}He(K^{-},#pi^{-})X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pim_MCosT_CM_Neutral", "^{3}He(K^{-},#pi^{-})X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  //  pi- + pi+
  new TH1F("pipi_Mass_CM_Neutral","Invariant mass of #pi^{+}#pi^{-};IM(#pi^{+}#pi^{-}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pipi_Momentum_CM_Neutral","Momentum of #pi^{+}#pi^{-};Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pipi_CosT_CM_Neutral","cos#theta of #pi^{+}#pi^{-};cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("pipi_Phi_CM_Neutral","#phi of #pi^{+}#pi^{-};#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("pipi_MMass_CM_Neutral", "^{3}He(K^{-},#pi^{+}#pi^{-})X missing mass;MM(#pi^{+}#pi^{-}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pipi_MMomentum_CM_Neutral", "^{3}He(K^{-},#pi^{+}#pi^{-})X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pipi_MCosT_CM_Neutral", "^{3}He(K^{-},#pi^{+}#pi^{-})X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // n
  new TH1F("n_Mass_CM_Neutral","Invariant mass of n;IM(n) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("n_Momentum_CM_Neutral","Momentum of n;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("n_CosT_CM_Neutral","cos#theta of n;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("n_Phi_CM_Neutral","#phi of n;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("n_MMass_CM_Neutral", "^{3}He(K^{-},n)X missing mass;MM(n) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("n_MMomentum_CM_Neutral", "^{3}He(K^{-},n)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("n_MCosT_CM_Neutral", "^{3}He(K^{-},n)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // n + pi+
  new TH1F("npip_Mass_CM_Neutral","Invariant mass of n#pi^{+};IM(n#pi^{+}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("npip_Momentum_CM_Neutral","Momentum of n#pi^{+};Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("npip_CosT_CM_Neutral","cos#theta of n#pi^{+};cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("npip_Phi_CM_Neutral","#phi of n#pi^{+};#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("npip_MMass_CM_Neutral", "^{3}He(K^{-},n#pi^{+})X missing mass;MM(n#pi^{+}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("npip_MMomentum_CM_Neutral", "^{3}He(K^{-},n#pi^{+})X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("npip_MCosT_CM_Neutral", "^{3}He(K^{-},n#pi^{+})X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // n + pi-
  new TH1F("npim_Mass_CM_Neutral","Invariant mass of n#pi^{-};IM(n#pi^{-}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("npim_Momentum_CM_Neutral","Momentum of n#pi^{-};Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("npim_CosT_CM_Neutral","cos#theta of n#pi^{-};cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("npim_Phi_CM_Neutral","#phi of n#pi^{-};#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("npim_MMass_CM_Neutral", "^{3}He(K^{-},n#pi^{-})X missing mass;MM(n#pi^{-}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("npim_MMomentum_CM_Neutral", "^{3}He(K^{-},n#pi^{-})X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("npim_MCosT_CM_Neutral", "^{3}He(K^{-},n#pi^{-})X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // n + pi- + pi+
  new TH1F("npipi_Mass_CM_Neutral","Invariant mass of n#pi^{+}#pi^{-};IM(n#pi^{+}#pi^{-}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("npipi_Momentum_CM_Neutral","Momentum of n#pi^{+}#pi^{-};Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("npipi_CosT_CM_Neutral","cos#theta of n#pi^{+}#pi^{-};cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("npipi_Phi_CM_Neutral","#phi of n#pi^{+}#pi^{-};#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("npipi_MMass_CM_Neutral", "^{3}He(K^{-},n#pi^{+}#pi^{-})X missing mass;MM(n#pi^{+}#pi^{-}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("npipi_MMomentum_CM_Neutral", "^{3}He(K^{-},n#pi^{+}#pi^{-})X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("npipi_MCosT_CM_Neutral", "^{3}He(K^{-},n#pi^{+}#pi^{-})X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // 2D plots
  new TH2F("npip_Mass_vsnpip_MMass_CM_Neutral","Invariant mass of n#pi^{+} vs. Missing mass of n#pi^{+};IM(n#pi^{+}) (GeV/c^{2});MM(n#pi^{+}) (GeV/c^{2})", 800, 1.0, 1.8, 800, 1.0, 1.8);
  new TH2F("npim_Mass_vsnpim_MMass_CM_Neutral","Invariant mass of n#pi^{-} vs. Missing mass of n#pi^{-};IM(n#pi^{-}) (GeV/c^{2});MM(n#pi^{-}) (GeV/c^{2})", 800, 1.0, 1.8, 800, 1.0, 1.8);
  new TH2F("npip_Mass_vsnpim_Mass_CM_Neutral","Invariant mass of n#pi^{+} vs. Invariant mass of n#pi^{-};IM(n#pi^{+}) (GeV/c^{2});IM(n#pi^{-}) (GeV/c^{2})", 800, 1.0, 1.8, 800, 1.0, 1.8);
  new TH2F("npip_MMass_vsnpim_MMass_CM_Neutral","Missing mass of n#pi^{+} vs. Missing mass of n#pi^{-};MM(n#pi^{+}) (GeV/c^{2});MM(n#pi^{-}) (GeV/c^{2})", 800, 1.0, 1.8, 800, 1.0, 1.8);
  new TH2F("n_MMass_vsnpipi_MMass_CM_Neutral","Missing mass of n vs. Missing mass of n#pi^{+}#pi^{-};MM(n) (GeV/c^{2});MM(n#pi^{+}#pi^{-}) (GeV/c^{2})", 800, 1.2, 2.0, 800, 0.8, 1.6);
  new TH2F("n_MMass_vsnpip_Mass_CM_Neutral","Missing mass of n vs. Invariant mass of n#pi^{+};MM(n) (GeV/c^{2});IM(n#pi^{+}) (GeV/c^{2})", 800, 1.2, 2.0, 800, 1.0, 1.8);
  new TH2F("n_MMass_vsnpip_MMass_CM_Neutral","Missing mass of n vs. Missing mass of n#pi^{+};MM(n) (GeV/c^{2});MM(n#pi^{+}) (GeV/c^{2})", 800, 1.2, 2.0, 800, 1.0, 1.8);
  new TH2F("n_MMass_vsnpim_Mass_CM_Neutral","Missing mass of n vs. Invariant mass of n#pi^{-};MM(n) (GeV/c^{2});IM(n#pi^{-}) (GeV/c^{2})", 800, 1.2, 2.0, 800, 1.0, 1.8);
  new TH2F("n_MMass_vsnpim_MMass_CM_Neutral","Missing mass of n vs. Missing mass of n#pi^{-};MM(n) (GeV/c^{2});MM(n#pi^{-}) (GeV/c^{2})", 800, 1.2, 2.0, 800, 1.0, 1.8);

  std::cout << "=== End of [FWDNeutral::Initialize] === " << std::endl;

  std::cout << "Define Histograms for FWD Charged particles" << std::endl;
  // FDC Tracking
  new TH1F( "FDC1_nTrack", "Number of Tracks in FDC1;Num. of Tracks;Counts", 20, 0, 20 );
  new TH1F( "FDC1_Chi2", "Chi-squared of FDC1;Chi-square;Counts", 500, 0, 50);
  new TH2F( "FDC1_Position", "Position at FDC1;X (cm);Y (cm)", 1000,-50,50, 1000,-50,50);
  // FWD Charged
  new TH1F( "FWDC_Beta", "#beta of FWD Charged;#beta;Counts", 1000,0.0,1.0);
  new TH1F( "FWDC_Overbeta", "1/#beta of FWD Charged;1/#beta;Counts", 5000,0.0,5.0);
  new TH2F( "FWDC_OverbetavsEnergy", "1/#beta vs. Energy deposit;1/#beta;Energy deposit (MeVee)", 1000, -0, 5, 500, 0, 50);
  new TH1F( "FWDC_Momentum", "Momentum of FWD Charged;Momentum (GeV/c);Counts", 2000,0.0,2.0);
  new TH1F( "FWDC_MomentumR", "Momentum of FWD Charged;Momentum (GeV/c);Counts", 2000,0.0,2.0);
  new TH2F( "FWDC_MomentumvsMomentumR", "Momentum vs Mometuntum_{R} of FWD Charged;Momentum (GeV/c);Momentum_{R} (GeV/c)", 2000,0.0,2.0, 2000,0.0,2.0);
  new TH1F( "FWDC_Mass2", "Mass2 of FWD Charged;Mass^{2} (GeV^{2}/c^{4});Counts", 600,-1.0,5.0);
  new TH2F( "FWDC_Mass2vsMomentum", "Mass2 vs Momentum of FWD Charged;Mass^{2} (GeV^{2}/c^{4});Momentum", 600,-1.0,5.0, 200, 0.0, 2.0);

  // pi1
  new TH1F("pim1_Mass_CM_Charged","Invariant mass of #pi^{-}_{1};IM(#pi^{-}_{1}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pim1_Momentum_CM_Charged","Momentum of #pi^{-}_{1};Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pim1_CosT_CM_Charged","cos#theta of #pi^{-}_{1};cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("pim1_Phi_CM_Charged","#phi of #pi^{-}_{1};#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("pim1_MMass_CM_Charged", "^{3}He(K^{-},#pi^{-}_{1})X missing mass;MM(#pi^{-}_{1}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pim1_MMomentum_CM_Charged", "^{3}He(K^{-},#pi^{-}_{1})X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pim1_MCosT_CM_Charged", "^{3}He(K^{-},#pi^{-}_{1})X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // pi2
  new TH1F("pim2_Mass_CM_Charged","Invariant mass of #pi^{-}_{2};IM(#pi^{-}_{2}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pim2_Momentum_CM_Charged","Momentum of #pi^{-}_{2};Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pim2_CosT_CM_Charged","cos#theta of #pi^{-}_{2};cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("pim2_Phi_CM_Charged","#phi of #pi^{-}_{2};#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("pim2_MMass_CM_Charged", "^{3}He(K^{-},#pi^{-}_{2})X missing mass;MM(#pi^{-}_{2}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pim2_MMomentum_CM_Charged", "^{3}He(K^{-},#pi^{-}_{2})X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pim2_MCosT_CM_Charged", "^{3}He(K^{-},#pi^{-}_{2})X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // pi1 + pi2
  new TH1F("pipi_Mass_CM_Charged","Invariant mass of #pi^{-}_{1}#pi^{-}_{2};IM(#pi^{-}_{1}#pi^{-}_{2}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pipi_Momentum_CM_Charged","Momentum of #pi^{-}_{1}#pi^{-}_{2};Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pipi_CosT_CM_Charged","cos#theta of #pi^{-}_{1}#pi^{-}_{2};cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("pipi_Phi_CM_Charged","#phi of #pi^{-}_{1}#pi^{-}_{2};#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("pipi_MMass_CM_Charged", "^{3}He(K^{-},#pi^{-}_{1}#pi^{-}_{2})X missing mass;MM(#pi^{-}_{1}#pi^{-}_{2}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("pipi_MMomentum_CM_Charged", "^{3}He(K^{-},#pi^{-}_{1}#pi^{-}_{2})X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("pipi_MCosT_CM_Charged", "^{3}He(K^{-},#pi^{-}_{1}#pi^{-}_{2})X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // p
  new TH1F("p_Mass_CM_Charged","Invariant mass of n;IM(p) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("p_Momentum_CM_Charged","Momentum of p;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("p_CosT_CM_Charged","cos#theta of p;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("p_Phi_CM_Charged","#phi of p;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("p_MMass_CM_Charged", "^{3}He(K^{-},p)X missing mass;MM(p) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("p_MMomentum_CM_Charged", "^{3}He(K^{-},p)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("p_MCosT_CM_Charged", "^{3}He(K^{-},p)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // p + pi1
  new TH1F("ppim1_Mass_CM_Charged","Invariant mass of p#pi^{-}_{1};IM(p#pi^{-}_{1}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("ppim1_Momentum_CM_Charged","Momentum of p#pi^{-}_{1};Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("ppim1_CosT_CM_Charged","cos#theta of p#pi^{-}_{1};cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("ppim1_Phi_CM_Charged","#phi of p#pi^{-}_{1};#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("ppim1_MMass_CM_Charged", "^{3}He(K^{-},p#pi^{-}_{1})X missing mass;MM(p#pi^{-}_{1}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("ppim1_MMomentum_CM_Charged", "^{3}He(K^{-},p#pi^{-}_{1})X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("ppim1_MCosT_CM_Charged", "^{3}He(K^{-},p#pi^{-}_{1})X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // p + pi2
  new TH1F("ppim2_Mass_CM_Charged","Invariant mass of p#pi^{-}_{2};IM(p#pi^{-}_{2}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("ppim2_Momentum_CM_Charged","Momentum of p#pi^{-}_{2};Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("ppim2_CosT_CM_Charged","cos#theta of p#pi^{-}_{2};cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("ppim2_Phi_CM_Charged","#phi of p#pi^{-}_{2};#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("ppim2_MMass_CM_Charged", "^{3}He(K^{-},p#pi^{-}_{2})X missing mass;MM(p#pi^{-}_{2}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("ppim2_MMomentum_CM_Charged", "^{3}He(K^{-},p#pi^{-}_{2})X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("ppim2_MCosT_CM_Charged", "^{3}He(K^{-},p#pi^{-}_{2})X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // p + pi1 + pi2
  new TH1F("ppipi_Mass_CM_Charged","Invariant mass of p#pi^{-}_{1}#pi^{-}_{2};IM(p#pi^{-}_{1}#pi^{-}_{2}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("ppipi_Momentum_CM_Charged","Momentum of p#pi^{-}_{1}#pi^{-}_{2};Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("ppipi_CosT_CM_Charged","cos#theta of p#pi^{-}_{1}#pi^{-}_{2};cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F("ppipi_Phi_CM_Charged","#phi of p#pi^{-}_{1}#pi^{-}_{2};#phi;Counts", 1600, 0.0, 3.2);
  new TH1F("ppipi_MMass_CM_Charged", "^{3}He(K^{-},p#pi^{-}_{1}#pi^{-}_{2})X missing mass;MM(p#pi^{-}_{1}#pi^{-}_{2}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F("ppipi_MMomentum_CM_Charged", "^{3}He(K^{-},p#pi^{-}_{1}#pi^{-}_{2})X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F("ppipi_MCosT_CM_Charged", "^{3}He(K^{-},p#pi^{-}_{1}#pi^{-}_{2})X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  // 2D plots
  new TH2F("ppim1_Mass_vsppim1_MMass_CM_Charged","Invariant mass of p#pi^{-}_{1} vs. Missing mass of p#pi^{-}_{1};IM(p#pi^{-}_{1}) (GeV/c^{2});MM(p#pi^{-}_{1}) (GeV/c^{2})", 800, 1.0, 1.8, 800, 1.0, 1.8);
  new TH2F("ppim2_Mass_vsppim2_MMass_CM_Charged","Invariant mass of p#pi^{-}_{2} vs. Missing mass of p#pi^{-}_{2};IM(p#pi^{-}_{2}) (GeV/c^{2});MM(p#pi^{-}_{2}) (GeV/c^{2})", 800, 1.0, 1.8, 800, 1.0, 1.8);
  new TH2F("ppim1_Mass_vsppim2_Mass_CM_Charged","Invariant mass of p#pi^{-}_{1} vs. Invariant mass of p#pi^{-}_{2};IM(p#pi^{-}_{1}) (GeV/c^{2});IM(p#pi^{-}_{2}) (GeV/c^{2})", 800, 1.0, 1.8, 800, 1.0, 1.8);
  new TH2F("ppim1_MMass_vsppim2_MMass_CM_Charged","Missing mass of p#pi^{-}_{1} vs. Missing mass of p#pi^{-}_{2};MM(p#pi^{-}_{1}) (GeV/c^{2});MM(p#pi^{-}_{2}) (GeV/c^{2})", 800, 1.0, 1.8, 800, 1.0, 1.8);
  new TH2F("p_MMass_vsppipi_MMass_CM_Charged","Missing mass of p vs. Missing mass of p#pi^{-}_{1}#pi^{-}_{2};MM(p) (GeV/c^{2});MM(p#pi^{-}_{1}#pi^{-}_{2}) (GeV/c^{2})", 800, 1.2, 2.0, 800, 0.8, 1.6);
  new TH2F("p_MMass_vsppim1_Mass_CM_Charged","Missing mass of p vs. Invariant mass of p#pi^{-}_{1};MM(p) (GeV/c^{2});IM(p#pi^{-}_{1}) (GeV/c^{2})", 800, 1.2, 2.0, 800, 1.0, 1.8);
  new TH2F("p_MMass_vsppim1_MMass_CM_Charged","Missing mass of p vs. Missing mass of p#pi^{-}_{1};MM(p) (GeV/c^{2});MM(p#pi^{-}_{1}) (GeV/c^{2})", 800, 1.2, 2.0, 800, 1.0, 1.8);
  new TH2F("p_MMass_vsppim2_Mass_CM_Charged","Missing mass of p vs. Invariant mass of p#pi^{-}_{2};MM(p) (GeV/c^{2});IM(p#pi^{-}_{2}) (GeV/c^{2})", 800, 1.2, 2.0, 800, 1.0, 1.8);
  new TH2F("p_MMass_vsppim2_MMass_CM_Charged","Missing mass of p vs. Missing mass of p#pi^{-}_{2};MM(p) (GeV/c^{2});MM(p#pi^{-}_{2}) (GeV/c^{2})", 800, 1.2, 2.0, 800, 1.0, 1.8);

  std::cout << "=== End of [FWDCharged::Initialize] === " << std::endl;

  return true;
}
