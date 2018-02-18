// MyAnalysisHeKpipipp.cpp

#include "MyAnalysisHeKpipipp.h"

MyAnalysisHeKpipipp::MyAnalysisHeKpipipp(TFile* rt, ConfMan* conf)
{
  Initialize(conf);
  CutCondition();
  Clear();
}

MyAnalysisHeKpipipp::~MyAnalysisHeKpipipp()
{
  Clear();
  rtFile->cd();
  rtFile->Write();
  rtFile->Close();
}

void MyAnalysisHeKpipipp::Clear()
{
}

bool MyAnalysisHeKpipipp::DoAnalysis(ConfMan* conf, EventHeader* header, BeamLineHitMan* blMan, BeamLineTrackMan* bltrackMan, CDSHitMan* cdsMan, CDSTrackingMan* cdstrackMan, Particle* particle)
{
  rtFile->cd();
  int ievent=0;

  FillHist("EventNumber",ievent); ievent++; /* All Events */
  DetectorList *dlist=DetectorList::GetInstance();
  if(particle->nBeam()!=1) return false;
  pBeam* beam = particle->beam(0);
  double beamtof = beam->bhdt0tof();

  /* Event selection */
  if(particle->nCDS()!=4) return true;
  FillHist("EventNumber",ievent); ievent++; /* Least 1 hit in CDS*/
  for(int it=0; it<particle->nCDS(); it++){
    pCDS* cds = particle->cdsi(it);
    if(cds->chi()>30) { return true; }
  }
  FillHist("EventNumber",ievent); ievent++; /* CDC chi2 < 30 */
  bool flag_lam_pos=false; bool flag_lam_neg=false; bool flag_k0=false;
  if(particle->nProton()!=2) return true;
  FillHist("EventNumber",ievent); ievent++; /* One Proton */
  pCDS* pro = 0; 
  pCDS* one = 0; 
  pCDS* proton = 0;
  if(particle->nPiminus()==1&&particle->nPiplus()==1){
    for(int i=0; i<particle->nProduct(); i++){
      pCDS* cds = particle->product(i);
      TLorentzVector Lcds = cds->GetLorentzVector();
      int comb = cds->comb();
      double mass = Lcds.M();
      if(comb==pow(2,CDS_Proton)+pow(2,CDS_PiMinus)){
        if(lamll<mass&&mass<lamul){
          flag_lam_pos=true;
          pro = cds;
          one = particle->pip(0);
          for(int j=0; j<particle->nProton(); j++){
            pCDS* tmpcds = particle->proton(j);
            if(tmpcds->id()!=cds->daughter1()&&tmpcds->id()!=cds->daughter2()){
              proton = tmpcds;
            }
          }
        }
      }
      if(comb==pow(2,CDS_PiPlus)+pow(2,CDS_PiMinus)){
        if(k0ll<mass&&mass<k0ul){
          flag_k0=true;
          pro = cds;
          one = particle->proton(0);
        }
      }
    }
  }
  else if(particle->nPiminus()==2){
    for(int i=0; i<particle->nProduct(); i++){
      pCDS* cds = particle->product(i);
      TLorentzVector Lcds = cds->GetLorentzVector();
      double mass = Lcds.M();
      int comb = cds->comb();
      if(comb==pow(2,CDS_Proton)+pow(2,CDS_PiMinus)){
        if(lamll<mass&&mass<lamul){
          flag_lam_neg=true;
          pro = cds;
          for(int j=0; j<particle->nPiminus(); j++){
            pCDS* tmpcds = particle->pim(j);
            if(tmpcds->id()!=cds->daughter1()&&tmpcds->id()!=cds->daughter2()){
              one = tmpcds;
            }
          }
          for(int j=0; j<particle->nProton(); j++){
            pCDS* tmpcds = particle->proton(j);
            if(tmpcds->id()!=cds->daughter1()&&tmpcds->id()!=cds->daughter2()){
              proton = tmpcds;
            }
          }
        }
      }
    }
  }
  else { return true; }
  FillHist("EventNumber",ievent); ievent++; /* two pi- or pi- and pi+ event */
  if(flag_k0){ return true; }
  FillHist("EventNumber",ievent); ievent++; /* K0 reconstruction */
  if(pro==0||one==0){ return true; }
  FillHist("EventNumber",ievent); ievent++; /* No product or one */
  if(flag_lam_neg&&flag_lam_pos){ return true; }
  FillHist("EventNumber",ievent); ievent++; /* Strange flag */

  // ############### //
  // Vertex decision //
  // ############### //
  TVector3 vertex;
  vertex = pro->vbeam();
  if(GeomTools::GetID(vertex)!=CID_Fiducial) { return true; }
  vertex = one->vbeam();
  if(GeomTools::GetID(vertex)!=CID_Fiducial) { return true; }
  vertex = proton->vbeam();
  if(GeomTools::GetID(vertex)!=CID_Fiducial) { return true; }
  FillHist("EventNumber",ievent); ievent++; /* Fiducial */

  if(flag_lam_pos) DoLambdaPosAna(particle,beam,pro,one,proton,vertex);
  if(flag_lam_neg) DoLambdaNegAna(particle,beam,pro,one,proton,vertex);

  return true;

}

bool MyAnalysisHeKpipipp::DoLambdaPosAna(Particle* particle, pBeam* beam, pCDS* lam, pCDS* pip, pCDS* proton, TVector3 vertex)
{
  TString AnaName = "LambdaPos";
  TString frame[2] = {"Lab","CM"};

  // Target
  TVector3 ZeroV;
  TLorentzVector Ltgt; Ltgt.SetVectM(ZeroV,ThreeHeMass);
  TLorentzVector Lbeam = beam->GetLorentzVector(vertex);
  TVector3 boost = (Ltgt+Lbeam).BoostVector();
  TLorentzVector cLtgt = Ltgt; cLtgt.Boost(-boost); TLorentzVector cLbeam = Lbeam; cLbeam.Boost(-boost);
  /* Lambda */
  TLorentzVector Llam = lam->GetLorentzVector();
  TLorentzVector cLlam = lam->GetLorentzVector(); cLlam.Boost(-boost);
  /* pi+ */
  TLorentzVector Lpip = pip->GetLorentzVector();
  TLorentzVector cLpip = pip->GetLorentzVector(); cLpip.Boost(-boost);
  /* proton */
  TLorentzVector Lp = proton->GetLorentzVector();
  TLorentzVector cLp = proton->GetLorentzVector(); cLp.Boost(-boost);
  /* Lambda */
  double lam_mass[2]  = { (Llam).M()                           , (cLlam).M()};
  double lam_mom[2]   = { (Llam).Vect().Mag()                  , (cLlam).Vect().Mag()};
  double lam_cost[2]  = { (Llam).Vect().CosTheta()             , (cLlam).Vect().CosTheta()};
  double lam_phi[2]   = { (Llam).Vect().Phi()                  , (cLlam).Vect().Phi()};
  double lam_mmass[2] = { (Ltgt+Lbeam-(Llam)).M()              , (cLtgt+cLbeam-(cLlam)).M()};
  double lam_mmom[2]  = { (Ltgt+Lbeam-(Llam)).Vect().Mag()     , (cLtgt+cLbeam-(cLlam)).Vect().Mag()};
  double lam_mcost[2] = { (Ltgt+Lbeam-(Llam)).Vect().CosTheta(), (cLtgt+cLbeam-(cLlam)).Vect().CosTheta()};
  /* pi+ */
  double pip_mass[2]  = { (Lpip).M()                           , (cLpip).M()};
  double pip_mom[2]   = { (Lpip).Vect().Mag()                  , (cLpip).Vect().Mag()};
  double pip_cost[2]  = { (Lpip).Vect().CosTheta()             , (cLpip).Vect().CosTheta()};
  double pip_phi[2]   = { (Lpip).Vect().Phi()                  , (cLpip).Vect().Phi()};
  double pip_mmass[2] = { (Ltgt+Lbeam-(Lpip)).M()              , (cLtgt+cLbeam-(cLpip)).M()};
  double pip_mmom[2]  = { (Ltgt+Lbeam-(Lpip)).Vect().Mag()     , (cLtgt+cLbeam-(cLpip)).Vect().Mag()};
  double pip_mcost[2] = { (Ltgt+Lbeam-(Lpip)).Vect().CosTheta(), (cLtgt+cLbeam-(cLpip)).Vect().CosTheta()};
  /* proton */
  double p_mass[2]  = { (Lp).M()                           , (cLp).M()};
  double p_mom[2]   = { (Lp).Vect().Mag()                  , (cLp).Vect().Mag()};
  double p_cost[2]  = { (Lp).Vect().CosTheta()             , (cLp).Vect().CosTheta()};
  double p_phi[2]   = { (Lp).Vect().Phi()                  , (cLp).Vect().Phi()};
  double p_mmass[2] = { (Ltgt+Lbeam-(Lp)).M()              , (cLtgt+cLbeam-(cLp)).M()};
  double p_mmom[2]  = { (Ltgt+Lbeam-(Lp)).Vect().Mag()     , (cLtgt+cLbeam-(cLp)).Vect().Mag()};
  double p_mcost[2] = { (Ltgt+Lbeam-(Lp)).Vect().CosTheta(), (cLtgt+cLbeam-(cLp)).Vect().CosTheta()};
  /* pi+ + proton */
  double pipp_mass[2]  = { (Lpip+Lp).M()                           , (cLpip+cLp).M()};
  double pipp_mom[2]   = { (Lpip+Lp).Vect().Mag()                  , (cLpip+cLp).Vect().Mag()};
  double pipp_cost[2]  = { (Lpip+Lp).Vect().CosTheta()             , (cLpip+cLp).Vect().CosTheta()};
  double pipp_phi[2]   = { (Lpip+Lp).Vect().Phi()                  , (cLpip+cLp).Vect().Phi()};
  double pipp_mmass[2] = { (Ltgt+Lbeam-(Lpip+Lp)).M()              , (cLtgt+cLbeam-(cLpip+cLp)).M()};
  double pipp_mmom[2]  = { (Ltgt+Lbeam-(Lpip+Lp)).Vect().Mag()     , (cLtgt+cLbeam-(cLpip+cLp)).Vect().Mag()};
  double pipp_mcost[2] = { (Ltgt+Lbeam-(Lpip+Lp)).Vect().CosTheta(), (cLtgt+cLbeam-(cLpip+cLp)).Vect().CosTheta()};
  /* Lambda + pi+ */
  double lampi_mass[2]  = { (Llam+Lpip).M()                           , (cLlam+cLpip).M()};
  double lampi_mom[2]   = { (Llam+Lpip).Vect().Mag()                  , (cLlam+cLpip).Vect().Mag()};
  double lampi_cost[2]  = { (Llam+Lpip).Vect().CosTheta()             , (cLlam+cLpip).Vect().CosTheta()};
  double lampi_phi[2]   = { (Llam+Lpip).Vect().Phi()                  , (cLlam+cLpip).Vect().Phi()};
  double lampi_mmass[2] = { (Ltgt+Lbeam-(Llam+Lpip)).M()              , (cLtgt+cLbeam-(cLlam+Lpip)).M()};
  double lampi_mmom[2]  = { (Ltgt+Lbeam-(Llam+Lpip)).Vect().Mag()     , (cLtgt+cLbeam-(cLlam+Lpip)).Vect().Mag()};
  double lampi_mcost[2] = { (Ltgt+Lbeam-(Llam+Lpip)).Vect().CosTheta(), (cLtgt+cLbeam-(cLlam+Lpip)).Vect().CosTheta()};
  /* Lambda + p */
  double lamp_mass[2]  = { (Llam+Lp).M()                           , (cLlam+cLp).M()};
  double lamp_mom[2]   = { (Llam+Lp).Vect().Mag()                  , (cLlam+cLp).Vect().Mag()};
  double lamp_cost[2]  = { (Llam+Lp).Vect().CosTheta()             , (cLlam+cLp).Vect().CosTheta()};
  double lamp_phi[2]   = { (Llam+Lp).Vect().Phi()                  , (cLlam+cLp).Vect().Phi()};
  double lamp_mmass[2] = { (Ltgt+Lbeam-(Llam+Lp)).M()              , (cLtgt+cLbeam-(cLlam+Lp)).M()};
  double lamp_mmom[2]  = { (Ltgt+Lbeam-(Llam+Lp)).Vect().Mag()     , (cLtgt+cLbeam-(cLlam+Lp)).Vect().Mag()};
  double lamp_mcost[2] = { (Ltgt+Lbeam-(Llam+Lp)).Vect().CosTheta(), (cLtgt+cLbeam-(cLlam+Lp)).Vect().CosTheta()};
  /* Lambda + pi+ */
  double lampip_mass[2]  = { (Llam+Lpip+Lp).M()                           , (cLlam+cLpip+cLp).M()};
  double lampip_mom[2]   = { (Llam+Lpip+Lp).Vect().Mag()                  , (cLlam+cLpip+cLp).Vect().Mag()};
  double lampip_cost[2]  = { (Llam+Lpip+Lp).Vect().CosTheta()             , (cLlam+cLpip+cLp).Vect().CosTheta()};
  double lampip_phi[2]   = { (Llam+Lpip+Lp).Vect().Phi()                  , (cLlam+cLpip+cLp).Vect().Phi()};
  double lampip_mmass[2] = { (Ltgt+Lbeam-(Llam+Lpip+Lp)).M()              , (cLtgt+cLbeam-(cLlam+Lpip+Lp)).M()};
  double lampip_mmom[2]  = { (Ltgt+Lbeam-(Llam+Lpip+Lp)).Vect().Mag()     , (cLtgt+cLbeam-(cLlam+Lpip+Lp)).Vect().Mag()};
  double lampip_mcost[2] = { (Ltgt+Lbeam-(Llam+Lpip+Lp)).Vect().CosTheta(), (cLtgt+cLbeam-(cLlam+Lpip+Lp)).Vect().CosTheta()};

  // Lambda
  for(int f=1; f<2; f++){
    FillHist(Form("lam_Mass_%s_%s",frame[f].Data(),AnaName.Data()),lam_mass[f]);
    FillHist(Form("lam_Momentum_%s_%s",frame[f].Data(),AnaName.Data()),lam_mom[f]);
    FillHist(Form("lam_CosT_%s_%s",frame[f].Data(),AnaName.Data()),lam_cost[f]);
    FillHist(Form("lam_Phi_%s_%s",frame[f].Data(),AnaName.Data()),lam_phi[f]);
    FillHist(Form("lam_MMass_%s_%s",frame[f].Data(),AnaName.Data()),lam_mmass[f]);
    FillHist(Form("lam_MMomentum_%s_%s",frame[f].Data(),AnaName.Data()),lam_mmom[f]);
    FillHist(Form("lam_MCosT_%s_%s",frame[f].Data(),AnaName.Data()),lam_mcost[f]);
  }
  // pi+
  for(int f=1; f<2; f++){
    FillHist(Form("pip_Mass_%s_%s",frame[f].Data(),AnaName.Data()),pip_mass[f]);
    FillHist(Form("pip_Momentum_%s_%s",frame[f].Data(),AnaName.Data()),pip_mom[f]);
    FillHist(Form("pip_CosT_%s_%s",frame[f].Data(),AnaName.Data()),pip_cost[f]);
    FillHist(Form("pip_Phi_%s_%s",frame[f].Data(),AnaName.Data()),pip_phi[f]);
    FillHist(Form("pip_MMass_%s_%s",frame[f].Data(),AnaName.Data()),pip_mmass[f]);
    FillHist(Form("pip_MMomentum_%s_%s",frame[f].Data(),AnaName.Data()),pip_mmom[f]);
    FillHist(Form("pip_MCosT_%s_%s",frame[f].Data(),AnaName.Data()),pip_mcost[f]);
  }
  // proton
  for(int f=1; f<2; f++){
    FillHist(Form("p_Mass_%s_%s",frame[f].Data(),AnaName.Data()),p_mass[f]);
    FillHist(Form("p_Momentum_%s_%s",frame[f].Data(),AnaName.Data()),p_mom[f]);
    FillHist(Form("p_CosT_%s_%s",frame[f].Data(),AnaName.Data()),p_cost[f]);
    FillHist(Form("p_Phi_%s_%s",frame[f].Data(),AnaName.Data()),p_phi[f]);
    FillHist(Form("p_MMass_%s_%s",frame[f].Data(),AnaName.Data()),p_mmass[f]);
    FillHist(Form("p_MMomentum_%s_%s",frame[f].Data(),AnaName.Data()),p_mmom[f]);
    FillHist(Form("p_MCosT_%s_%s",frame[f].Data(),AnaName.Data()),p_mcost[f]);
  }
  // pi+ + proton
  for(int f=1; f<2; f++){
    FillHist(Form("pipp_Mass_%s_%s",frame[f].Data(),AnaName.Data()),pipp_mass[f]);
    FillHist(Form("pipp_Momentum_%s_%s",frame[f].Data(),AnaName.Data()),pipp_mom[f]);
    FillHist(Form("pipp_CosT_%s_%s",frame[f].Data(),AnaName.Data()),pipp_cost[f]);
    FillHist(Form("pipp_Phi_%s_%s",frame[f].Data(),AnaName.Data()),pipp_phi[f]);
    FillHist(Form("pipp_MMass_%s_%s",frame[f].Data(),AnaName.Data()),pipp_mmass[f]);
    FillHist(Form("pipp_MMomentum_%s_%s",frame[f].Data(),AnaName.Data()),pipp_mmom[f]);
    FillHist(Form("pipp_MCosT_%s_%s",frame[f].Data(),AnaName.Data()),pipp_mcost[f]);
  }
  // Lambda + pi+
  for(int f=1; f<2; f++){
    FillHist(Form("lampi_Mass_%s_%s",frame[f].Data(),AnaName.Data()),lampi_mass[f]);
    FillHist(Form("lampi_Momentum_%s_%s",frame[f].Data(),AnaName.Data()),lampi_mom[f]);
    FillHist(Form("lampi_CosT_%s_%s",frame[f].Data(),AnaName.Data()),lampi_cost[f]);
    FillHist(Form("lampi_Phi_%s_%s",frame[f].Data(),AnaName.Data()),lampi_phi[f]);
    FillHist(Form("lampi_MMass_%s_%s",frame[f].Data(),AnaName.Data()),lampi_mmass[f]);
    FillHist(Form("lampi_MMomentum_%s_%s",frame[f].Data(),AnaName.Data()),lampi_mmom[f]);
    FillHist(Form("lampi_MCosT_%s_%s",frame[f].Data(),AnaName.Data()),lampi_mcost[f]);
  }
  // Lambda + proton
  for(int f=1; f<2; f++){
    FillHist(Form("lamp_Mass_%s_%s",frame[f].Data(),AnaName.Data()),lamp_mass[f]);
    FillHist(Form("lamp_Momentum_%s_%s",frame[f].Data(),AnaName.Data()),lamp_mom[f]);
    FillHist(Form("lamp_CosT_%s_%s",frame[f].Data(),AnaName.Data()),lamp_cost[f]);
    FillHist(Form("lamp_Phi_%s_%s",frame[f].Data(),AnaName.Data()),lamp_phi[f]);
    FillHist(Form("lamp_MMass_%s_%s",frame[f].Data(),AnaName.Data()),lamp_mmass[f]);
    FillHist(Form("lamp_MMomentum_%s_%s",frame[f].Data(),AnaName.Data()),lamp_mmom[f]);
    FillHist(Form("lamp_MCosT_%s_%s",frame[f].Data(),AnaName.Data()),lamp_mcost[f]);
  }
  // Lambda + pi+ + proton
  for(int f=1; f<2; f++){
    FillHist(Form("lampip_Mass_%s_%s",frame[f].Data(),AnaName.Data()),lampip_mass[f]);
    FillHist(Form("lampip_Momentum_%s_%s",frame[f].Data(),AnaName.Data()),lampip_mom[f]);
    FillHist(Form("lampip_CosT_%s_%s",frame[f].Data(),AnaName.Data()),lampip_cost[f]);
    FillHist(Form("lampip_Phi_%s_%s",frame[f].Data(),AnaName.Data()),lampip_phi[f]);
    FillHist(Form("lampip_MMass_%s_%s",frame[f].Data(),AnaName.Data()),lampip_mmass[f]);
    FillHist(Form("lampip_MMomentum_%s_%s",frame[f].Data(),AnaName.Data()),lampip_mmom[f]);
    FillHist(Form("lampip_MCosT_%s_%s",frame[f].Data(),AnaName.Data()),lampip_mcost[f]);
  }

  /* 2D Plot */
  char* fill1; char* fill2;
  fill1="lampi_Mass"; /* vs */ fill2="lampi_MCosT";
  FillHist(Form("%s_vs_%s_%s_%s",fill1,fill2,frame[1].Data(),AnaName.Data()),lampi_mass[1],lampi_mcost[1],1);
  fill1="lampi_MMass"; /* vs */ fill2="lampi_MCosT";
  FillHist(Form("%s_vs_%s_%s_%s",fill1,fill2,frame[1].Data(),AnaName.Data()),lampi_mmass[1],lampi_mcost[1],1);
  fill1="lampi_Mass"; /* vs */ fill2="lampi_MMass";
  FillHist(Form("%s_vs_%s_%s_%s",fill1,fill2,frame[1].Data(),AnaName.Data()),lampi_mass[1],lampi_mmass[1],1);
  fill1="lampip_Mass"; /* vs */ fill2="lampip_MMass";
  FillHist(Form("%s_vs_%s_%s_%s",fill1,fill2,frame[1].Data(),AnaName.Data()),lampip_mass[1],lampip_mmass[1],1);
  fill1="lampi_Mass"; /* vs */ fill2="lampip_MMass";
  FillHist(Form("%s_vs_%s_%s_%s",fill1,fill2,frame[1].Data(),AnaName.Data()),lampi_mass[1],lampip_mmass[1],1);
  fill1="lampi_MMass"; /* vs */ fill2="lampip_MMass";
  FillHist(Form("%s_vs_%s_%s_%s",fill1,fill2,frame[1].Data(),AnaName.Data()),lampi_mmass[1],lampip_mmass[1],1);

  return true;
}
bool MyAnalysisHeKpipipp::DoLambdaNegAna(Particle* particle, pBeam* beam, pCDS* lam, pCDS* pim, pCDS* proton, TVector3 vertex)
{
  TString AnaName = "LambdaNeg";
  TString frame[2] = {"Lab","CM"};

  // Target
  TVector3 ZeroV;
  TLorentzVector Ltgt; Ltgt.SetVectM(ZeroV,ThreeHeMass);
  TLorentzVector Lbeam = beam->GetLorentzVector(vertex);
  TVector3 boost = (Ltgt+Lbeam).BoostVector();
  TLorentzVector cLtgt = Ltgt; cLtgt.Boost(-boost); TLorentzVector cLbeam = Lbeam; cLbeam.Boost(-boost);
  /* Lambda */
  TLorentzVector Llam = lam->GetLorentzVector();
  TLorentzVector cLlam = lam->GetLorentzVector(); cLlam.Boost(-boost);
  /* pi- */
  TLorentzVector Lpim = pim->GetLorentzVector();
  TLorentzVector cLpim = pim->GetLorentzVector(); cLpim.Boost(-boost);
  /* proton */
  TLorentzVector Lp = proton->GetLorentzVector();
  TLorentzVector cLp = proton->GetLorentzVector(); cLp.Boost(-boost);
  /* Lambda */
  double lam_mass[2]  = { (Llam).M()                           , (cLlam).M()};
  double lam_mom[2]   = { (Llam).Vect().Mag()                  , (cLlam).Vect().Mag()};
  double lam_cost[2]  = { (Llam).Vect().CosTheta()             , (cLlam).Vect().CosTheta()};
  double lam_phi[2]   = { (Llam).Vect().Phi()                  , (cLlam).Vect().Phi()};
  double lam_mmass[2] = { (Ltgt+Lbeam-(Llam)).M()              , (cLtgt+cLbeam-(cLlam)).M()};
  double lam_mmom[2]  = { (Ltgt+Lbeam-(Llam)).Vect().Mag()     , (cLtgt+cLbeam-(cLlam)).Vect().Mag()};
  double lam_mcost[2] = { (Ltgt+Lbeam-(Llam)).Vect().CosTheta(), (cLtgt+cLbeam-(cLlam)).Vect().CosTheta()};
  /* pi- */
  double pim_mass[2]  = { (Lpim).M()                           , (cLpim).M()};
  double pim_mom[2]   = { (Lpim).Vect().Mag()                  , (cLpim).Vect().Mag()};
  double pim_cost[2]  = { (Lpim).Vect().CosTheta()             , (cLpim).Vect().CosTheta()};
  double pim_phi[2]   = { (Lpim).Vect().Phi()                  , (cLpim).Vect().Phi()};
  double pim_mmass[2] = { (Ltgt+Lbeam-(Lpim)).M()              , (cLtgt+cLbeam-(cLpim)).M()};
  double pim_mmom[2]  = { (Ltgt+Lbeam-(Lpim)).Vect().Mag()     , (cLtgt+cLbeam-(cLpim)).Vect().Mag()};
  double pim_mcost[2] = { (Ltgt+Lbeam-(Lpim)).Vect().CosTheta(), (cLtgt+cLbeam-(cLpim)).Vect().CosTheta()};
  /* proton */
  double p_mass[2]  = { (Lp).M()                           , (cLp).M()};
  double p_mom[2]   = { (Lp).Vect().Mag()                  , (cLp).Vect().Mag()};
  double p_cost[2]  = { (Lp).Vect().CosTheta()             , (cLp).Vect().CosTheta()};
  double p_phi[2]   = { (Lp).Vect().Phi()                  , (cLp).Vect().Phi()};
  double p_mmass[2] = { (Ltgt+Lbeam-(Lp)).M()              , (cLtgt+cLbeam-(cLp)).M()};
  double p_mmom[2]  = { (Ltgt+Lbeam-(Lp)).Vect().Mag()     , (cLtgt+cLbeam-(cLp)).Vect().Mag()};
  double p_mcost[2] = { (Ltgt+Lbeam-(Lp)).Vect().CosTheta(), (cLtgt+cLbeam-(cLp)).Vect().CosTheta()};
  /* pi- + proton */
  double pimp_mass[2]  = { (Lpim+Lp).M()                           , (cLpim+cLp).M()};
  double pimp_mom[2]   = { (Lpim+Lp).Vect().Mag()                  , (cLpim+cLp).Vect().Mag()};
  double pimp_cost[2]  = { (Lpim+Lp).Vect().CosTheta()             , (cLpim+cLp).Vect().CosTheta()};
  double pimp_phi[2]   = { (Lpim+Lp).Vect().Phi()                  , (cLpim+cLp).Vect().Phi()};
  double pimp_mmass[2] = { (Ltgt+Lbeam-(Lpim+Lp)).M()              , (cLtgt+cLbeam-(cLpim+cLp)).M()};
  double pimp_mmom[2]  = { (Ltgt+Lbeam-(Lpim+Lp)).Vect().Mag()     , (cLtgt+cLbeam-(cLpim+cLp)).Vect().Mag()};
  double pimp_mcost[2] = { (Ltgt+Lbeam-(Lpim+Lp)).Vect().CosTheta(), (cLtgt+cLbeam-(cLpim+cLp)).Vect().CosTheta()};
  /* Lambda + pi- */
  double lampi_mass[2]  = { (Llam+Lpim).M()                           , (cLlam+cLpim).M()};
  double lampi_mom[2]   = { (Llam+Lpim).Vect().Mag()                  , (cLlam+cLpim).Vect().Mag()};
  double lampi_cost[2]  = { (Llam+Lpim).Vect().CosTheta()             , (cLlam+cLpim).Vect().CosTheta()};
  double lampi_phi[2]   = { (Llam+Lpim).Vect().Phi()                  , (cLlam+cLpim).Vect().Phi()};
  double lampi_mmass[2] = { (Ltgt+Lbeam-(Llam+Lpim)).M()              , (cLtgt+cLbeam-(cLlam+Lpim)).M()};
  double lampi_mmom[2]  = { (Ltgt+Lbeam-(Llam+Lpim)).Vect().Mag()     , (cLtgt+cLbeam-(cLlam+Lpim)).Vect().Mag()};
  double lampi_mcost[2] = { (Ltgt+Lbeam-(Llam+Lpim)).Vect().CosTheta(), (cLtgt+cLbeam-(cLlam+Lpim)).Vect().CosTheta()};
  /* Lambda + p */
  double lamp_mass[2]  = { (Llam+Lp).M()                           , (cLlam+cLp).M()};
  double lamp_mom[2]   = { (Llam+Lp).Vect().Mag()                  , (cLlam+cLp).Vect().Mag()};
  double lamp_cost[2]  = { (Llam+Lp).Vect().CosTheta()             , (cLlam+cLp).Vect().CosTheta()};
  double lamp_phi[2]   = { (Llam+Lp).Vect().Phi()                  , (cLlam+cLp).Vect().Phi()};
  double lamp_mmass[2] = { (Ltgt+Lbeam-(Llam+Lp)).M()              , (cLtgt+cLbeam-(cLlam+Lp)).M()};
  double lamp_mmom[2]  = { (Ltgt+Lbeam-(Llam+Lp)).Vect().Mag()     , (cLtgt+cLbeam-(cLlam+Lp)).Vect().Mag()};
  double lamp_mcost[2] = { (Ltgt+Lbeam-(Llam+Lp)).Vect().CosTheta(), (cLtgt+cLbeam-(cLlam+Lp)).Vect().CosTheta()};
  /* Lambda + pi+ */
  double lampip_mass[2]  = { (Llam+Lpim+Lp).M()                           , (cLlam+cLpim+cLp).M()};
  double lampip_mom[2]   = { (Llam+Lpim+Lp).Vect().Mag()                  , (cLlam+cLpim+cLp).Vect().Mag()};
  double lampip_cost[2]  = { (Llam+Lpim+Lp).Vect().CosTheta()             , (cLlam+cLpim+cLp).Vect().CosTheta()};
  double lampip_phi[2]   = { (Llam+Lpim+Lp).Vect().Phi()                  , (cLlam+cLpim+cLp).Vect().Phi()};
  double lampip_mmass[2] = { (Ltgt+Lbeam-(Llam+Lpim+Lp)).M()              , (cLtgt+cLbeam-(cLlam+Lpim+Lp)).M()};
  double lampip_mmom[2]  = { (Ltgt+Lbeam-(Llam+Lpim+Lp)).Vect().Mag()     , (cLtgt+cLbeam-(cLlam+Lpim+Lp)).Vect().Mag()};
  double lampip_mcost[2] = { (Ltgt+Lbeam-(Llam+Lpim+Lp)).Vect().CosTheta(), (cLtgt+cLbeam-(cLlam+Lpim+Lp)).Vect().CosTheta()};

  // Lambda
  for(int f=1; f<2; f++){
    FillHist(Form("lam_Mass_%s_%s",frame[f].Data(),AnaName.Data()),lam_mass[f]);
    FillHist(Form("lam_Momentum_%s_%s",frame[f].Data(),AnaName.Data()),lam_mom[f]);
    FillHist(Form("lam_CosT_%s_%s",frame[f].Data(),AnaName.Data()),lam_cost[f]);
    FillHist(Form("lam_Phi_%s_%s",frame[f].Data(),AnaName.Data()),lam_phi[f]);
    FillHist(Form("lam_MMass_%s_%s",frame[f].Data(),AnaName.Data()),lam_mmass[f]);
    FillHist(Form("lam_MMomentum_%s_%s",frame[f].Data(),AnaName.Data()),lam_mmom[f]);
    FillHist(Form("lam_MCosT_%s_%s",frame[f].Data(),AnaName.Data()),lam_mcost[f]);
  }
  // pi-
  for(int f=1; f<2; f++){
    FillHist(Form("pim_Mass_%s_%s",frame[f].Data(),AnaName.Data()),pim_mass[f]);
    FillHist(Form("pim_Momentum_%s_%s",frame[f].Data(),AnaName.Data()),pim_mom[f]);
    FillHist(Form("pim_CosT_%s_%s",frame[f].Data(),AnaName.Data()),pim_cost[f]);
    FillHist(Form("pim_Phi_%s_%s",frame[f].Data(),AnaName.Data()),pim_phi[f]);
    FillHist(Form("pim_MMass_%s_%s",frame[f].Data(),AnaName.Data()),pim_mmass[f]);
    FillHist(Form("pim_MMomentum_%s_%s",frame[f].Data(),AnaName.Data()),pim_mmom[f]);
    FillHist(Form("pim_MCosT_%s_%s",frame[f].Data(),AnaName.Data()),pim_mcost[f]);
  }
  // proton
  for(int f=1; f<2; f++){
    FillHist(Form("p_Mass_%s_%s",frame[f].Data(),AnaName.Data()),p_mass[f]);
    FillHist(Form("p_Momentum_%s_%s",frame[f].Data(),AnaName.Data()),p_mom[f]);
    FillHist(Form("p_CosT_%s_%s",frame[f].Data(),AnaName.Data()),p_cost[f]);
    FillHist(Form("p_Phi_%s_%s",frame[f].Data(),AnaName.Data()),p_phi[f]);
    FillHist(Form("p_MMass_%s_%s",frame[f].Data(),AnaName.Data()),p_mmass[f]);
    FillHist(Form("p_MMomentum_%s_%s",frame[f].Data(),AnaName.Data()),p_mmom[f]);
    FillHist(Form("p_MCosT_%s_%s",frame[f].Data(),AnaName.Data()),p_mcost[f]);
  }
  // pi- + proton
  for(int f=1; f<2; f++){
    FillHist(Form("pimp_Mass_%s_%s",frame[f].Data(),AnaName.Data()),pimp_mass[f]);
    FillHist(Form("pimp_Momentum_%s_%s",frame[f].Data(),AnaName.Data()),pimp_mom[f]);
    FillHist(Form("pimp_CosT_%s_%s",frame[f].Data(),AnaName.Data()),pimp_cost[f]);
    FillHist(Form("pimp_Phi_%s_%s",frame[f].Data(),AnaName.Data()),pimp_phi[f]);
    FillHist(Form("pimp_MMass_%s_%s",frame[f].Data(),AnaName.Data()),pimp_mmass[f]);
    FillHist(Form("pimp_MMomentum_%s_%s",frame[f].Data(),AnaName.Data()),pimp_mmom[f]);
    FillHist(Form("pimp_MCosT_%s_%s",frame[f].Data(),AnaName.Data()),pimp_mcost[f]);
  }
  // Lambda + pi+
  for(int f=1; f<2; f++){
    FillHist(Form("lampi_Mass_%s_%s",frame[f].Data(),AnaName.Data()),lampi_mass[f]);
    FillHist(Form("lampi_Momentum_%s_%s",frame[f].Data(),AnaName.Data()),lampi_mom[f]);
    FillHist(Form("lampi_CosT_%s_%s",frame[f].Data(),AnaName.Data()),lampi_cost[f]);
    FillHist(Form("lampi_Phi_%s_%s",frame[f].Data(),AnaName.Data()),lampi_phi[f]);
    FillHist(Form("lampi_MMass_%s_%s",frame[f].Data(),AnaName.Data()),lampi_mmass[f]);
    FillHist(Form("lampi_MMomentum_%s_%s",frame[f].Data(),AnaName.Data()),lampi_mmom[f]);
    FillHist(Form("lampi_MCosT_%s_%s",frame[f].Data(),AnaName.Data()),lampi_mcost[f]);
  }
  // Lambda + proton
  for(int f=1; f<2; f++){
    FillHist(Form("lamp_Mass_%s_%s",frame[f].Data(),AnaName.Data()),lamp_mass[f]);
    FillHist(Form("lamp_Momentum_%s_%s",frame[f].Data(),AnaName.Data()),lamp_mom[f]);
    FillHist(Form("lamp_CosT_%s_%s",frame[f].Data(),AnaName.Data()),lamp_cost[f]);
    FillHist(Form("lamp_Phi_%s_%s",frame[f].Data(),AnaName.Data()),lamp_phi[f]);
    FillHist(Form("lamp_MMass_%s_%s",frame[f].Data(),AnaName.Data()),lamp_mmass[f]);
    FillHist(Form("lamp_MMomentum_%s_%s",frame[f].Data(),AnaName.Data()),lamp_mmom[f]);
    FillHist(Form("lamp_MCosT_%s_%s",frame[f].Data(),AnaName.Data()),lamp_mcost[f]);
  }
  // Lambda + pi+ + proton
  for(int f=1; f<2; f++){
    FillHist(Form("lampip_Mass_%s_%s",frame[f].Data(),AnaName.Data()),lampip_mass[f]);
    FillHist(Form("lampip_Momentum_%s_%s",frame[f].Data(),AnaName.Data()),lampip_mom[f]);
    FillHist(Form("lampip_CosT_%s_%s",frame[f].Data(),AnaName.Data()),lampip_cost[f]);
    FillHist(Form("lampip_Phi_%s_%s",frame[f].Data(),AnaName.Data()),lampip_phi[f]);
    FillHist(Form("lampip_MMass_%s_%s",frame[f].Data(),AnaName.Data()),lampip_mmass[f]);
    FillHist(Form("lampip_MMomentum_%s_%s",frame[f].Data(),AnaName.Data()),lampip_mmom[f]);
    FillHist(Form("lampip_MCosT_%s_%s",frame[f].Data(),AnaName.Data()),lampip_mcost[f]);
  }

  /* 2D Plot */
  char* fill1; char* fill2;
  fill1="lampi_Mass"; /* vs */ fill2="lampi_MCosT";
  FillHist(Form("%s_vs_%s_%s_%s",fill1,fill2,frame[1].Data(),AnaName.Data()),lampi_mass[1],lampi_mcost[1],1);
  fill1="lampi_MMass"; /* vs */ fill2="lampi_MCosT";
  FillHist(Form("%s_vs_%s_%s_%s",fill1,fill2,frame[1].Data(),AnaName.Data()),lampi_mmass[1],lampi_mcost[1],1);
  fill1="lampi_Mass"; /* vs */ fill2="lampi_MMass";
  FillHist(Form("%s_vs_%s_%s_%s",fill1,fill2,frame[1].Data(),AnaName.Data()),lampi_mass[1],lampi_mmass[1],1);
  fill1="lampip_Mass"; /* vs */ fill2="lampip_MMass";
  FillHist(Form("%s_vs_%s_%s_%s",fill1,fill2,frame[1].Data(),AnaName.Data()),lampip_mass[1],lampip_mmass[1],1);
  fill1="lampi_Mass"; /* vs */ fill2="lampip_MMass";
  FillHist(Form("%s_vs_%s_%s_%s",fill1,fill2,frame[1].Data(),AnaName.Data()),lampi_mass[1],lampip_mmass[1],1);
  fill1="lampi_MMass"; /* vs */ fill2="lampip_MMass";
  FillHist(Form("%s_vs_%s_%s_%s",fill1,fill2,frame[1].Data(),AnaName.Data()),lampi_mmass[1],lampip_mmass[1],1);

  return true;
}


bool MyAnalysisHeKpipipp::DoK0Ana(Particle* particle, pBeam* beam, pCDS* k0, pCDS* proton, TVector3 vertex)
{
  TString AnaName = "K0";
  TString frame[2] = {"Lab","CM"};

  // Target
  TVector3 ZeroV;
  TLorentzVector Ltgt; Ltgt.SetVectM(ZeroV,ThreeHeMass);
  TLorentzVector Lbeam = beam->GetLorentzVector(vertex);
  TVector3 boost = (Ltgt+Lbeam).BoostVector();
  TLorentzVector cLtgt = Ltgt; cLtgt.Boost(-boost); TLorentzVector cLbeam = Lbeam; cLbeam.Boost(-boost);
  /* K0 */
  TLorentzVector Lk0 = k0->GetLorentzVector();
  TLorentzVector cLk0 = k0->GetLorentzVector(); cLk0.Boost(-boost);
  /* proton */
  TLorentzVector Lproton = proton->GetLorentzVector();
  TLorentzVector cLproton = proton->GetLorentzVector(); cLproton.Boost(-boost);
  /* K0 */
  double k0_mass[2]  = { (Lk0).M()                           , (cLk0).M()};
  double k0_mom[2]   = { (Lk0).Vect().Mag()                  , (cLk0).Vect().Mag()};
  double k0_cost[2]  = { (Lk0).Vect().CosTheta()             , (cLk0).Vect().CosTheta()};
  double k0_phi[2]   = { (Lk0).Vect().Phi()                  , (cLk0).Vect().Phi()};
  double k0_mmass[2] = { (Ltgt+Lbeam-(Lk0)).M()              , (cLtgt+cLbeam-(cLk0)).M()};
  double k0_mmom[2]  = { (Ltgt+Lbeam-(Lk0)).Vect().Mag()     , (cLtgt+cLbeam-(cLk0)).Vect().Mag()};
  double k0_mcost[2] = { (Ltgt+Lbeam-(Lk0)).Vect().CosTheta(), (cLtgt+cLbeam-(cLk0)).Vect().CosTheta()};
  /* proton */
  double p_mass[2]  = { (Lproton).M()                           , (cLproton).M()};
  double p_mom[2]   = { (Lproton).Vect().Mag()                  , (cLproton).Vect().Mag()};
  double p_cost[2]  = { (Lproton).Vect().CosTheta()             , (cLproton).Vect().CosTheta()};
  double p_phi[2]   = { (Lproton).Vect().Phi()                  , (cLproton).Vect().Phi()};
  double p_mmass[2] = { (Ltgt+Lbeam-(Lproton)).M()              , (cLtgt+cLbeam-(cLproton)).M()};
  double p_mmom[2]  = { (Ltgt+Lbeam-(Lproton)).Vect().Mag()     , (cLtgt+cLbeam-(cLproton)).Vect().Mag()};
  double p_mcost[2] = { (Ltgt+Lbeam-(Lproton)).Vect().CosTheta(), (cLtgt+cLbeam-(cLproton)).Vect().CosTheta()};
  /* Lambda + pi+ */
  double k0p_mass[2]  = { (Lk0+Lproton).M()                           , (cLk0+cLproton).M()};
  double k0p_mom[2]   = { (Lk0+Lproton).Vect().Mag()                  , (cLk0+cLproton).Vect().Mag()};
  double k0p_cost[2]  = { (Lk0+Lproton).Vect().CosTheta()             , (cLk0+cLproton).Vect().CosTheta()};
  double k0p_phi[2]   = { (Lk0+Lproton).Vect().Phi()                  , (cLk0+cLproton).Vect().Phi()};
  double k0p_mmass[2] = { (Ltgt+Lbeam-(Lk0+Lproton)).M()              , (cLtgt+cLbeam-(cLk0+Lproton)).M()};
  double k0p_mmom[2]  = { (Ltgt+Lbeam-(Lk0+Lproton)).Vect().Mag()     , (cLtgt+cLbeam-(cLk0+Lproton)).Vect().Mag()};
  double k0p_mcost[2] = { (Ltgt+Lbeam-(Lk0+Lproton)).Vect().CosTheta(), (cLtgt+cLbeam-(cLk0+Lproton)).Vect().CosTheta()};

  // K0
  for(int f=1; f<2; f++){
    FillHist(Form("k0_Mass_%s_%s",frame[f].Data(),AnaName.Data()),k0_mass[f]);
    FillHist(Form("k0_Momentum_%s_%s",frame[f].Data(),AnaName.Data()),k0_mom[f]);
    FillHist(Form("k0_CosT_%s_%s",frame[f].Data(),AnaName.Data()),k0_cost[f]);
    FillHist(Form("k0_Phi_%s_%s",frame[f].Data(),AnaName.Data()),k0_phi[f]);
    FillHist(Form("k0_MMass_%s_%s",frame[f].Data(),AnaName.Data()),k0_mmass[f]);
    FillHist(Form("k0_MMomentum_%s_%s",frame[f].Data(),AnaName.Data()),k0_mmom[f]);
    FillHist(Form("k0_MCosT_%s_%s",frame[f].Data(),AnaName.Data()),k0_mcost[f]);
  }
  // proton
  for(int f=1; f<2; f++){
    FillHist(Form("p_Mass_%s_%s",frame[f].Data(),AnaName.Data()),p_mass[f]);
    FillHist(Form("p_Momentum_%s_%s",frame[f].Data(),AnaName.Data()),p_mom[f]);
    FillHist(Form("p_CosT_%s_%s",frame[f].Data(),AnaName.Data()),p_cost[f]);
    FillHist(Form("p_Phi_%s_%s",frame[f].Data(),AnaName.Data()),p_phi[f]);
    FillHist(Form("p_MMass_%s_%s",frame[f].Data(),AnaName.Data()),p_mmass[f]);
    FillHist(Form("p_MMomentum_%s_%s",frame[f].Data(),AnaName.Data()),p_mmom[f]);
    FillHist(Form("p_MCosT_%s_%s",frame[f].Data(),AnaName.Data()),p_mcost[f]);
  }
  // Lambda + pi+
  for(int f=1; f<2; f++){
    FillHist(Form("k0p_Mass_%s_%s",frame[f].Data(),AnaName.Data()),k0p_mass[f]);
    FillHist(Form("k0p_Momentum_%s_%s",frame[f].Data(),AnaName.Data()),k0p_mom[f]);
    FillHist(Form("k0p_CosT_%s_%s",frame[f].Data(),AnaName.Data()),k0p_cost[f]);
    FillHist(Form("k0p_Phi_%s_%s",frame[f].Data(),AnaName.Data()),k0p_phi[f]);
    FillHist(Form("k0p_MMass_%s_%s",frame[f].Data(),AnaName.Data()),k0p_mmass[f]);
    FillHist(Form("k0p_MMomentum_%s_%s",frame[f].Data(),AnaName.Data()),k0p_mmom[f]);
    FillHist(Form("k0p_MCosT_%s_%s",frame[f].Data(),AnaName.Data()),k0p_mcost[f]);
  }

  /* 2D Plot */
  char* fill1; char* fill2;
  fill1="k0p_Mass"; /* vs */ fill2="k0p_MCosT";
  FillHist(Form("%s_vs_%s_%s_%s",fill1,fill2,frame[1].Data(),AnaName.Data()),k0p_mass[1],k0p_mcost[1],1);
  fill1="k0p_MMass"; /* vs */ fill2="k0p_MCosT";
  FillHist(Form("%s_vs_%s_%s_%s",fill1,fill2,frame[1].Data(),AnaName.Data()),k0p_mmass[1],k0p_mcost[1],1);
  fill1="k0p_Mass"; /* vs */ fill2="k0p_MMass";
  FillHist(Form("%s_vs_%s_%s_%s",fill1,fill2,frame[1].Data(),AnaName.Data()),k0p_mass[1],k0p_mmass[1],1);

  return true;

}


bool MyAnalysisHeKpipipp::FillHist(TString name, double val1, int weight)
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

bool MyAnalysisHeKpipipp::FillHist(TString name, TString val1, int weight)
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

bool MyAnalysisHeKpipipp::FillHist(TString name, double val1, double val2, int weight)
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

bool MyAnalysisHeKpipipp::FillHist(TString name, TString val1, TString val2, int weight)
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

void MyAnalysisHeKpipipp::CutCondition()
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

bool MyAnalysisHeKpipipp::Initialize(ConfMan* confMan)
{
  std::cout << "### MyAnalysisHeKpipipp::Initialize ###" << std::endl;

  std::string ofname = confMan->GetOutFileName();
  ofname.insert(ofname.find(".root"),"_anaHeKpipipp");

  rtFile =  new TFile( Form("%s",ofname.c_str()), "RECREATE");
  rtFile -> cd();
  char* ana = "LambdaPos";

  new TH1F( "EventNumber", "Number of Events", 20, 0, 20);
  // LambdaPos
  ana = "LambdaPos";
  std::cout << "Define Histograms for LambdaPos" << std::endl;
  new TH1F(Form("lam_Mass_CM_%s",ana),"Invariant mass of #Lambda;IM(#Lambda) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F(Form("lam_Momentum_CM_%s",ana),"Momentum of #Lambda;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F(Form("lam_CosT_CM_%s",ana),"cos#theta of #Lambda;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F(Form("lam_Phi_CM_%s",ana),"#phi of #Lambda;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F(Form("lam_MMass_CM_%s",ana), "^{3}He(K^{-},#Lambda)X missing mass;MM(#Lambda) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F(Form("lam_MMomentum_CM_%s",ana), "^{3}He(K^{-},#Lambda)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F(Form("lam_MCosT_CM_%s",ana), "^{3}He(K^{-},n)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F(Form("pip_Mass_CM_%s",ana),"Invariant mass of #pi^{+};IM(#pi^{+}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F(Form("pip_Momentum_CM_%s",ana),"Momentum of #pi^{+};Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F(Form("pip_CosT_CM_%s",ana),"cos#theta of #pi^{+};cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F(Form("pip_Phi_CM_%s",ana),"#phi of #pi^{+};#phi;Counts", 1600, 0.0, 3.2);
  new TH1F(Form("pip_MMass_CM_%s",ana), "^{3}He(K^{-},#pi^{+})X missing mass;MM(#pi^{+}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F(Form("pip_MMomentum_CM_%s",ana), "^{3}He(K^{-},#pi^{+})X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F(Form("pip_MCosT_CM_%s",ana), "^{3}He(K^{-},n)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F(Form("p_Mass_CM_%s",ana),"Invariant mass of p;IM(p) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F(Form("p_Momentum_CM_%s",ana),"Momentum of p;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F(Form("p_CosT_CM_%s",ana),"cos#theta of p;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F(Form("p_Phi_CM_%s",ana),"#phi of p;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F(Form("p_MMass_CM_%s",ana), "^{3}He(K^{-},p)X missing mass;MM(p) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F(Form("p_MMomentum_CM_%s",ana), "^{3}He(K^{-},p)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F(Form("p_MCosT_CM_%s",ana), "^{3}He(K^{-},n)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F(Form("pipp_Mass_CM_%s",ana),"Invariant mass of #pi^{+}p;IM(#pi^{+}p) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F(Form("pipp_Momentum_CM_%s",ana),"Momentum of #pi^{+}p;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F(Form("pipp_CosT_CM_%s",ana),"cos#theta of #pi^{+}p;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F(Form("pipp_Phi_CM_%s",ana),"#phi of #pi^{+}p;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F(Form("pipp_MMass_CM_%s",ana), "^{3}He(K^{-},#pi^{+}p)X missing mass;MM(#pi^{+}p) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F(Form("pipp_MMomentum_CM_%s",ana), "^{3}He(K^{-},#pi^{+}p)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F(Form("pipp_MCosT_CM_%s",ana), "^{3}He(K^{-},n)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F(Form("lampi_Mass_CM_%s",ana),"Invariant mass of #Lambda#pi^{+};IM(#Lambda#pi^{+}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F(Form("lampi_Momentum_CM_%s",ana),"Momentum of #Lambda#pi^{+};Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F(Form("lampi_CosT_CM_%s",ana),"cos#theta of #Lambda#pi^{+};cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F(Form("lampi_Phi_CM_%s",ana),"#phi of #Lambda#pi^{+};#phi;Counts", 1600, 0.0, 3.2);
  new TH1F(Form("lampi_MMass_CM_%s",ana), "^{3}He(K^{-},#Lambda#pi^{+})X missing mass;MM(#Lambda#pi^{+}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F(Form("lampi_MMomentum_CM_%s",ana), "^{3}He(K^{-},#Lambda#pi^{+})X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F(Form("lampi_MCosT_CM_%s",ana), "^{3}He(K^{-},#Lambda#pi^{+})X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F(Form("lamp_Mass_CM_%s",ana),"Invariant mass of #Lambdap;IM(#Lambdap) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F(Form("lamp_Momentum_CM_%s",ana),"Momentum of #Lambdap;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F(Form("lamp_CosT_CM_%s",ana),"cos#theta of #Lambdap;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F(Form("lamp_Phi_CM_%s",ana),"#phi of #Lambdap;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F(Form("lamp_MMass_CM_%s",ana), "^{3}He(K^{-},#Lambdap)X missing mass;MM(#Lambdap) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F(Form("lamp_MMomentum_CM_%s",ana), "^{3}He(K^{-},#Lambdap)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F(Form("lamp_MCosT_CM_%s",ana), "^{3}He(K^{-},#Lambdap)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F(Form("lampip_Mass_CM_%s",ana),"Invariant mass of #Lambda#pi^{-};IM(#Lambda#pi^{-}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F(Form("lampip_Momentum_CM_%s",ana),"Momentum of #Lambda#pi^{-};Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F(Form("lampip_CosT_CM_%s",ana),"cos#theta of #Lambda#pi^{-};cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F(Form("lampip_Phi_CM_%s",ana),"#phi of #Lambda#pi^{-};#phi;Counts", 1600, 0.0, 3.2);
  new TH1F(Form("lampip_MMass_CM_%s",ana), "^{3}He(K^{-},#Lambda#pi^{-})X missing mass;MM(#Lambda#pi^{-}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F(Form("lampip_MMomentum_CM_%s",ana), "^{3}He(K^{-},#Lambda#pi^{-})X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F(Form("lampip_MCosT_CM_%s",ana), "^{3}He(K^{-},#Lambda#pi^{-})X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  /* 2D Plot */
  new TH2F(Form("lampi_Mass_vs_lampi_MCosT_CM_%s",ana), "Invariant mass of #Lambda#pi^{+} vs. ^{3}He(K^{-},#Lambda#pi^{+})X cos#theta;IM(#Lambda#pi^{+});cos#theta;", 500, 0.0, 5.0, 300, -1.5, 1.5);
  new TH2F(Form("lampi_MMass_vs_lampi_MCosT_CM_%s",ana), "Missing mass of ^{3}He(K^{-},#Lambda#pi^{+})X vs. ^{3}He(K^{-},#Lambda#pi^{+})X cos#theta;MM(#Lambda#pi^{+});cos#theta;", 500, 0.0, 5.0, 300, -1.5, 1.5);
  new TH2F(Form("lampi_Mass_vs_lampi_MMass_CM_%s",ana), "Invariant mass of #Lambda#pi^{+} vs. Missing mass of ^{3}He(K^{-},#Lambda#pi^{+})X;IM(#Lambda#pi^{+}) (GeV/c^{2});MM(#Lambda#pi^{+}) (GeV/c^{2});", 500, 0.0, 5.0, 500, 0.0, 5.0);
  new TH2F(Form("lampip_Mass_vs_lampip_MMass_CM_%s",ana), "Invariant mass of #Lambda#pi^{+}p vs. Missing mass of ^{3}He(K^{-},#Lambda#pi^{+}p)X;IM(#Lambda#pi^{+}p) (GeV/c^{2});MM(#Lambda#pi^{+}p) (GeV/c^{2});", 500, 0.0, 5.0, 500, 0.0, 5.0);
  new TH2F(Form("lampi_Mass_vs_lampip_MMass_CM_%s",ana), "Invariant mass of #Lambda#pi^{+} vs. Missing mass of ^{3}He(K^{-},#Lambda#pi^{+}p)X;IM(#Lambda#pi^{+}) (GeV/c^{2});MM(#Lambda#pi^{+}p) (GeV/c^{2});", 500, 0.0, 5.0, 500, 0.0, 5.0);
  new TH2F(Form("lampi_MMass_vs_lampip_MMass_CM_%s",ana), "Missing mass of #Lambda#pi^{+} vs. Missing mass of ^{3}He(K^{-},#Lambda#pi^{+}p)X;IM(#Lambda#pi^{+}) (GeV/c^{2});MM(#Lambda#pi^{+}p) (GeV/c^{2});", 500, 0.0, 5.0, 500, 0.0, 5.0);
  std::cout << "=== End of [LambdaPos::Initialize] === " << std::endl;
  // LambdaNeg
  ana = "LambdaNeg";
  std::cout << "Define Histograms for LambdaPos" << std::endl;
  new TH1F(Form("lam_Mass_CM_%s",ana),"Invariant mass of #Lambda;IM(#Lambda) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F(Form("lam_Momentum_CM_%s",ana),"Momentum of #Lambda;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F(Form("lam_CosT_CM_%s",ana),"cos#theta of #Lambda;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F(Form("lam_Phi_CM_%s",ana),"#phi of #Lambda;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F(Form("lam_MMass_CM_%s",ana), "^{3}He(K^{-},#Lambda)X missing mass;MM(#Lambda) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F(Form("lam_MMomentum_CM_%s",ana), "^{3}He(K^{-},#Lambda)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F(Form("lam_MCosT_CM_%s",ana), "^{3}He(K^{-},n)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F(Form("pim_Mass_CM_%s",ana),"Invariant mass of #pi^{-};IM(#pi^{-}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F(Form("pim_Momentum_CM_%s",ana),"Momentum of #pi^{-};Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F(Form("pim_CosT_CM_%s",ana),"cos#theta of #pi^{-};cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F(Form("pim_Phi_CM_%s",ana),"#phi of #pi^{-};#phi;Counts", 1600, 0.0, 3.2);
  new TH1F(Form("pim_MMass_CM_%s",ana), "^{3}He(K^{-},#pi^{-})X missing mass;MM(#pi^{-}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F(Form("pim_MMomentum_CM_%s",ana), "^{3}He(K^{-},#pi^{-})X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F(Form("pim_MCosT_CM_%s",ana), "^{3}He(K^{-},n)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F(Form("p_Mass_CM_%s",ana),"Invariant mass of p;IM(p) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F(Form("p_Momentum_CM_%s",ana),"Momentum of p;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F(Form("p_CosT_CM_%s",ana),"cos#theta of p;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F(Form("p_Phi_CM_%s",ana),"#phi of p;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F(Form("p_MMass_CM_%s",ana), "^{3}He(K^{-},p)X missing mass;MM(p) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F(Form("p_MMomentum_CM_%s",ana), "^{3}He(K^{-},p)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F(Form("p_MCosT_CM_%s",ana), "^{3}He(K^{-},n)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F(Form("pimp_Mass_CM_%s",ana),"Invariant mass of #pi^{-}p;IM(#pi^{-}p) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F(Form("pimp_Momentum_CM_%s",ana),"Momentum of #pi^{-}p;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F(Form("pimp_CosT_CM_%s",ana),"cos#theta of #pi^{-}p;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F(Form("pimp_Phi_CM_%s",ana),"#phi of #pi^{-}p;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F(Form("pimp_MMass_CM_%s",ana), "^{3}He(K^{-},#pi^{-}p)X missing mass;MM(#pi^{-}p) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F(Form("pimp_MMomentum_CM_%s",ana), "^{3}He(K^{-},#pi^{-}p)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F(Form("pimp_MCosT_CM_%s",ana), "^{3}He(K^{-},n)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F(Form("lampi_Mass_CM_%s",ana),"Invariant mass of #Lambda#pi^{-};IM(#Lambda#pi^{-}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F(Form("lampi_Momentum_CM_%s",ana),"Momentum of #Lambda#pi^{-};Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F(Form("lampi_CosT_CM_%s",ana),"cos#theta of #Lambda#pi^{-};cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F(Form("lampi_Phi_CM_%s",ana),"#phi of #Lambda#pi^{-};#phi;Counts", 1600, 0.0, 3.2);
  new TH1F(Form("lampi_MMass_CM_%s",ana), "^{3}He(K^{-},#Lambda#pi^{-})X missing mass;MM(#Lambda#pi^{-}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F(Form("lampi_MMomentum_CM_%s",ana), "^{3}He(K^{-},#Lambda#pi^{-})X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F(Form("lampi_MCosT_CM_%s",ana), "^{3}He(K^{-},#Lambda#pi^{-})X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F(Form("lamp_Mass_CM_%s",ana),"Invariant mass of #Lambdap;IM(#Lambdap) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F(Form("lamp_Momentum_CM_%s",ana),"Momentum of #Lambdap;Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F(Form("lamp_CosT_CM_%s",ana),"cos#theta of #Lambdap;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F(Form("lamp_Phi_CM_%s",ana),"#phi of #Lambdap;#phi;Counts", 1600, 0.0, 3.2);
  new TH1F(Form("lamp_MMass_CM_%s",ana), "^{3}He(K^{-},#Lambdap)X missing mass;MM(#Lambdap) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F(Form("lamp_MMomentum_CM_%s",ana), "^{3}He(K^{-},#Lambdap)X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F(Form("lamp_MCosT_CM_%s",ana), "^{3}He(K^{-},#Lambdap)X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F(Form("lampip_Mass_CM_%s",ana),"Invariant mass of #Lambda#pi^{-};IM(#Lambda#pi^{-}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F(Form("lampip_Momentum_CM_%s",ana),"Momentum of #Lambda#pi^{-};Momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F(Form("lampip_CosT_CM_%s",ana),"cos#theta of #Lambda#pi^{-};cos#theta;Counts", 3000, -1.5, 1.5);
  new TH1F(Form("lampip_Phi_CM_%s",ana),"#phi of #Lambda#pi^{-};#phi;Counts", 1600, 0.0, 3.2);
  new TH1F(Form("lampip_MMass_CM_%s",ana), "^{3}He(K^{-},#Lambda#pi^{-})X missing mass;MM(#Lambda#pi^{-}) (GeV/c^{2});Counts", 5000, 0.0, 5.0);
  new TH1F(Form("lampip_MMomentum_CM_%s",ana), "^{3}He(K^{-},#Lambda#pi^{-})X missing momentum;Missing momentum (GeV/c);Counts", 2000, 0.0, 2.0);
  new TH1F(Form("lampip_MCosT_CM_%s",ana), "^{3}He(K^{-},#Lambda#pi^{-})X cos#theta;cos#theta;Counts", 3000, -1.5, 1.5);
  /* 2D Plot */
  new TH2F(Form("lampi_Mass_vs_lampi_MCosT_CM_%s",ana), "Invariant mass of #Lambda#pi^{-} vs. ^{3}He(K^{-},#Lambda#pi^{-})X cos#theta;IM(#Lambda#pi^{-});cos#theta;", 500, 0.0, 5.0, 300, -1.5, 1.5);
  new TH2F(Form("lampi_MMass_vs_lampi_MCosT_CM_%s",ana), "Missing mass of ^{3}He(K^{-},#Lambda#pi^{-})X vs. ^{3}He(K^{-},#Lambda#pi^{-})X cos#theta;MM(#Lambda#pi^{-});cos#theta;", 500, 0.0, 5.0, 300, -1.5, 1.5);
  new TH2F(Form("lampi_Mass_vs_lampi_MMass_CM_%s",ana), "Invariant mass of #Lambda#pi^{-} vs. Missing mass of ^{3}He(K^{-},#Lambda#pi^{-})X;IM(#Lambda#pi^{-}) (GeV/c^{2});MM(#Lambda#pi^{-}) (GeV/c^{2});", 500, 0.0, 5.0, 500, 0.0, 5.0);
  new TH2F(Form("lampip_Mass_vs_lampip_MMass_CM_%s",ana), "Invariant mass of #Lambda#pi^{-}p vs. Missing mass of ^{3}He(K^{-},#Lambda#pi^{-}p)X;IM(#Lambda#pi^{-}p) (GeV/c^{2});MM(#Lambda#pi^{-}p) (GeV/c^{2});", 500, 0.0, 5.0, 500, 0.0, 5.0);
  new TH2F(Form("lampi_Mass_vs_lampip_MMass_CM_%s",ana), "Invariant mass of #Lambda#pi^{-} vs. Missing mass of ^{3}He(K^{-},#Lambda#pi^{-}p)X;IM(#Lambda#pi^{-}) (GeV/c^{2});MM(#Lambda#pi^{-}p) (GeV/c^{2});", 500, 0.0, 5.0, 500, 0.0, 5.0);
  new TH2F(Form("lampi_MMass_vs_lampip_MMass_CM_%s",ana), "Missing mass of #Lambda#pi^{-} vs. Missing mass of ^{3}He(K^{-},#Lambda#pi^{-}p)X;IM(#Lambda#pi^{-}) (GeV/c^{2});MM(#Lambda#pi^{-}p) (GeV/c^{2});", 500, 0.0, 5.0, 500, 0.0, 5.0);
  std::cout << "=== End of [LambdaPos::Initialize] === " << std::endl;

  return true;
}
