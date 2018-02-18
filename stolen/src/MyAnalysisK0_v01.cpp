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
  //if(!header->trigmode(Mode_KCDH2)) return true;
  //if(header->IsTrig(Trig_KCDH2)&&header->IsTrig(Trig_1stMix)){
  //}
  //else {
  //  return true;
  //}

  if(particle->nPiplus()!=1) return true;
  if(particle->nPiminus()!=1) return true;
  if(particle->nCDS()!=2) return true;

  // Trigger Pattern
  FillHist("TriggerPattern",0);
  for( int i=0; i<20; i++ ){
    int val = header->pattern(i);
    if( 0<val ){
      FillHist("TriggerPattern",i);
      FillHist(Form("Trigger%02d_TDC",i),val);
    }
  }

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
  if(fpro==-1) return true; /* there is no vertex in fiducial volume */

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
    nc->CalcMom(beam,vertex);
    FillHist("FWD_OverbetavsMomentum",1.0/nc->beta(),nc->mom().Mag());
    FillHist("FWD_OverbetavsEnergy",1.0/nc->beta(),nc->energy());
    FillHist("FWD_TOFvsMomentum",nc->tof(),nc->mom().Mag());
    FillHist("FWD_HitPosition",nc->hitpos().X(),nc->hitpos().Y());
    double seg = (nc->seg()-1)%16+1;
    double lay = (nc->seg()-1)/16+1;
    FillHist("FWD_HitSegment",seg,lay);
    if(nc->pid()==F_Neutron){
      if(nc->energy()>8.0 && (time>9998||time>nc->time())){
        time = nc->time();
        fnc = it;
      }
    }
  }
  //if(fnc==-1) return true; /* there is no NC hit */

  pCDS* pip = particle->pip(0);
  pCDS* pim = particle->pim(0);
  pCDS* pro=0;
  if(fpro!=-1) pro=particle->product(fpro);

  /* CDS particle analysis */
  if(pro==0) return true;
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
  double pip_pbdca = pip->pbdca();
  double pim_pbdca = pim->pbdca();
  double pro_pbdca = pro->pbdca();
  double pro_vbdis = pro->vbdis();
  double pro_vdis  = pro->vdis();
  double m_dxdz = (Ltgt+Lbeam-Lpro).Vect().X()/(Ltgt+Lbeam-Lpro).Vect().Z();
  double m_dydz = (Ltgt+Lbeam-Lpro).Vect().Y()/(Ltgt+Lbeam-Lpro).Vect().Z();
  TVector3 mnpos;
  TVector3 ncpos;
  conf->GetGeomMapManager()->GetGPos(CID_NC,1,ncpos); ncpos.SetZ(ncpos.Z()-2.5);
  mnpos.SetX(pro->vertex().X() + m_dxdz*(ncpos.Z()-pro->vertex().Z()));
  mnpos.SetY(pro->vertex().Y() + m_dydz*(ncpos.Z()-pro->vertex().Z()));
  mnpos.SetZ(ncpos.Z());
  TString pid[8] = {"pip","p","d","t","he","pim","k","o"};
  TString frame[2] = {"Lab","CM"};

  /* flag */
  bool mnflag = false; /* missing neutron flag */
  if(mnll<pro_mmass[1]&&pro_mmass[1]<mnul) mnflag = true;
  bool sblmnflag = false; /* lower sideband of missing neutron flag */
  if(sblmnll<pro_mmass[1]&&pro_mmass[1]<sblmnul) sblmnflag = true;
  bool sbumnflag = false; /* upper sideband of missing neutron flag */
  if(sbumnll<pro_mmass[1]&&pro_mmass[1]<sbumnul) sbumnflag = true;
  bool k0flag = false; /* k0 flag */
  if(k0ll<pro_mass[1]&&pro_mass[1]<k0ul) k0flag = true;
  bool sblk0flag = false; /* lower sideband of k0 flag */
  if(sblk0ll<pro_mass[1]&&pro_mass[1]<sblk0ul) sblk0flag = true;
  bool sbuk0flag = false; /* upper sideband of k0 flag */
  if(sbuk0ll<pro_mass[1]&&pro_mass[1]<sbuk0ul) sbuk0flag = true;

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
    FillHist(Form("pipi_IMvsDCA_%s",frame[f].Data()),pro_mass[f],pro_pbdca);
    FillHist(Form("pipi_CosTvsMomentum_%s",frame[f].Data()),pro_cost[f],pro_mom[f]);
    FillHist(Form("pipi_CosTvsPhi_%s",frame[f].Data()),pro_cost[f],pro_phi[f]);
    FillHist(Form("pKpipi_MMvsMomentum_%s",frame[f].Data()),pro_mmass[f],pro_mmom[f]);
    FillHist(Form("pKpipi_MMvsCosT_%s",frame[f].Data()),pro_mmass[f],pro_mcost[f]);
  }
  FillHist("pip_PBDCA",pip_pbdca);
  FillHist("pim_PBDCA",pim_pbdca);
  FillHist("Pro_PBDCA",pro_pbdca);
  FillHist("Pro_VBDIS",pro_vbdis);
  FillHist("Pro_VDIS",pro_vdis);

  /* with k0 flag */
  if(k0flag){
    for(int f=1; f<2; f++){
      FillHist(Form("pipi_k0_IMvsMomentum_%s",frame[f].Data()),pro_mass[f],pro_mom[f]);
      FillHist(Form("pipi_k0_IMvsCosT_%s",frame[f].Data()),pro_mass[f],pro_cost[f]);
      FillHist(Form("pipi_k0_IMvsDCA_%s",frame[f].Data()),pro_mass[f],pro_pbdca);
      FillHist(Form("pKpipi_k0_MMvsMomentum_%s",frame[f].Data()),pro_mmass[f],pro_mmom[f]);
      FillHist(Form("pKpipi_k0_MMvsCosT_%s",frame[f].Data()),pro_mmass[f],pro_mcost[f]);
    }
  }
  /* with L-sidebiand k0 flag */
  if(sblk0flag){
    for(int f=1; f<2; f++){
      FillHist(Form("pipi_sblk0_IMvsMomentum_%s",frame[f].Data()),pro_mass[f],pro_mom[f]);
      FillHist(Form("pipi_sblk0_IMvsCosT_%s",frame[f].Data()),pro_mass[f],pro_cost[f]);
      FillHist(Form("pipi_sblk0_IMvsDCA_%s",frame[f].Data()),pro_mass[f],pro_pbdca);
      FillHist(Form("pKpipi_sblk0_MMvsMomentum_%s",frame[f].Data()),pro_mmass[f],pro_mmom[f]);
      FillHist(Form("pKpipi_sblk0_MMvsCosT_%s",frame[f].Data()),pro_mmass[f],pro_mcost[f]);
    }
  }
  /* with U-sidebiand k0 flag */
  if(sbuk0flag){
    for(int f=1; f<2; f++){
      FillHist(Form("pipi_sbuk0_IMvsMomentum_%s",frame[f].Data()),pro_mass[f],pro_mom[f]);
      FillHist(Form("pipi_sbuk0_IMvsCosT_%s",frame[f].Data()),pro_mass[f],pro_cost[f]);
      FillHist(Form("pipi_sbuk0_IMvsDCA_%s",frame[f].Data()),pro_mass[f],pro_pbdca);
      FillHist(Form("pKpipi_sbuk0_MMvsMomentum_%s",frame[f].Data()),pro_mmass[f],pro_mmom[f]);
      FillHist(Form("pKpipi_sbuk0_MMvsCosT_%s",frame[f].Data()),pro_mmass[f],pro_mcost[f]);
    }
  }
  /* with sidebiand k0 flag */
  if(sblk0flag||sbuk0flag){
    for(int f=1; f<2; f++){
      FillHist(Form("pipi_sbk0_IMvsMomentum_%s",frame[f].Data()),pro_mass[f],pro_mom[f]);
      FillHist(Form("pipi_sbk0_IMvsCosT_%s",frame[f].Data()),pro_mass[f],pro_cost[f]);
      FillHist(Form("pipi_sbk0_IMvsDCA_%s",frame[f].Data()),pro_mass[f],pro_pbdca);
      FillHist(Form("pKpipi_sbk0_MMvsMomentum_%s",frame[f].Data()),pro_mmass[f],pro_mmom[f]);
      FillHist(Form("pKpipi_sbk0_MMvsCosT_%s",frame[f].Data()),pro_mmass[f],pro_mcost[f]);
    }
  }
  /* without k0 flag */
  if(!k0flag){
    for(int f=1; f<2; f++){
      FillHist(Form("pipi_wok0_IMvsMomentum_%s",frame[f].Data()),pro_mass[f],pro_mom[f]);
      FillHist(Form("pipi_wok0_IMvsCosT_%s",frame[f].Data()),pro_mass[f],pro_cost[f]);
      FillHist(Form("pipi_wok0_IMvsDCA_%s",frame[f].Data()),pro_mass[f],pro_pbdca);
      FillHist(Form("pKpipi_wok0_MMvsMomentum_%s",frame[f].Data()),pro_mmass[f],pro_mmom[f]);
      FillHist(Form("pKpipi_wok0_MMvsCosT_%s",frame[f].Data()),pro_mmass[f],pro_mcost[f]);
    }
  }
  /* with missing neutron flag */
  if(mnflag){
    for(int f=1; f<2; f++){
      FillHist(Form("pipi_mn_IMvsMomentum_%s",frame[f].Data()),pro_mass[f],pro_mom[f]);
      FillHist(Form("pipi_mn_IMvsCosT_%s",frame[f].Data()),pro_mass[f],pro_cost[f]);
      FillHist(Form("pipi_mn_IMvsDCA_%s",frame[f].Data()),pro_mass[f],pro_pbdca);
      FillHist(Form("pKpipi_mn_MMvsMomentum_%s",frame[f].Data()),pro_mmass[f],pro_mmom[f]);
      FillHist(Form("pKpipi_mn_MMvsCosT_%s",frame[f].Data()),pro_mmass[f],pro_mcost[f]);
    }
    if(k0flag){
      for(int f=1; f<2; f++){
        FillHist(Form("pipi_k0mn_IMvsMomentum_%s",frame[f].Data()),pro_mass[f],pro_mom[f]);
        FillHist(Form("pipi_k0mn_IMvsCosT_%s",frame[f].Data()),pro_mass[f],pro_cost[f]);
        FillHist(Form("pipi_k0mn_IMvsDCA_%s",frame[f].Data()),pro_mass[f],pro_pbdca);
        FillHist(Form("pKpipi_k0mn_MMvsMomentum_%s",frame[f].Data()),pro_mmass[f],pro_mmom[f]);
        FillHist(Form("pKpipi_k0mn_MMvsCosT_%s",frame[f].Data()),pro_mmass[f],pro_mcost[f]);
      }
    }
  /* with L-sidebiand k0 flag */
    if(sblk0flag){
      for(int f=1; f<2; f++){
        FillHist(Form("pipi_sblk0mn_IMvsMomentum_%s",frame[f].Data()),pro_mass[f],pro_mom[f]);
        FillHist(Form("pipi_sblk0mn_IMvsCosT_%s",frame[f].Data()),pro_mass[f],pro_cost[f]);
        FillHist(Form("pipi_sblk0mn_IMvsDCA_%s",frame[f].Data()),pro_mass[f],pro_pbdca);
        FillHist(Form("pKpipi_sblk0mn_MMvsMomentum_%s",frame[f].Data()),pro_mmass[f],pro_mmom[f]);
        FillHist(Form("pKpipi_sblk0mn_MMvsCosT_%s",frame[f].Data()),pro_mmass[f],pro_mcost[f]);
      }
    }
    /* with U-sidebiand k0 flag */
    if(sbuk0flag){
      for(int f=1; f<2; f++){
        FillHist(Form("pipi_sbuk0mn_IMvsMomentum_%s",frame[f].Data()),pro_mass[f],pro_mom[f]);
        FillHist(Form("pipi_sbuk0mn_IMvsCosT_%s",frame[f].Data()),pro_mass[f],pro_cost[f]);
        FillHist(Form("pipi_sbuk0mn_IMvsDCA_%s",frame[f].Data()),pro_mass[f],pro_pbdca);
        FillHist(Form("pKpipi_sbuk0mn_MMvsMomentum_%s",frame[f].Data()),pro_mmass[f],pro_mmom[f]);
        FillHist(Form("pKpipi_sbuk0mn_MMvsCosT_%s",frame[f].Data()),pro_mmass[f],pro_mcost[f]);
      }
    }
    /* with sidebiand k0 flag */
    if(sblk0flag||sbuk0flag){
      for(int f=1; f<2; f++){
        FillHist(Form("pipi_sbk0mn_IMvsMomentum_%s",frame[f].Data()),pro_mass[f],pro_mom[f]);
        FillHist(Form("pipi_sbk0mn_IMvsCosT_%s",frame[f].Data()),pro_mass[f],pro_cost[f]);
        FillHist(Form("pipi_sbk0mn_IMvsDCA_%s",frame[f].Data()),pro_mass[f],pro_pbdca);
        FillHist(Form("pKpipi_sbk0mn_MMvsMomentum_%s",frame[f].Data()),pro_mmass[f],pro_mmom[f]);
        FillHist(Form("pKpipi_sbk0mn_MMvsCosT_%s",frame[f].Data()),pro_mmass[f],pro_mcost[f]);
      }
    }
  }
  /* without missing neutron flag */
  if(!mnflag){
    for(int f=1; f<2; f++){
      FillHist(Form("pipi_womn_IMvsMomentum_%s",frame[f].Data()),pro_mass[f],pro_mom[f]);
      FillHist(Form("pipi_womn_IMvsCosT_%s",frame[f].Data()),pro_mass[f],pro_cost[f]);
      FillHist(Form("pipi_womn_IMvsDCA_%s",frame[f].Data()),pro_mass[f],pro_pbdca);
      FillHist(Form("pKpipi_womn_MMvsMomentum_%s",frame[f].Data()),pro_mmass[f],pro_mmom[f]);
      FillHist(Form("pKpipi_womn_MMvsCosT_%s",frame[f].Data()),pro_mmass[f],pro_mcost[f]);
    }
    if(k0flag){
      for(int f=1; f<2; f++){
        FillHist(Form("pipi_k0womn_IMvsMomentum_%s",frame[f].Data()),pro_mass[f],pro_mom[f]);
        FillHist(Form("pipi_k0womn_IMvsCosT_%s",frame[f].Data()),pro_mass[f],pro_cost[f]);
        FillHist(Form("pipi_k0womn_IMvsDCA_%s",frame[f].Data()),pro_mass[f],pro_pbdca);
        FillHist(Form("pKpipi_k0womn_MMvsMomentum_%s",frame[f].Data()),pro_mmass[f],pro_mmom[f]);
        FillHist(Form("pKpipi_k0womn_MMvsCosT_%s",frame[f].Data()),pro_mmass[f],pro_mcost[f]);
      }
    }
  }
  int mnposflag[5] = {0,0,0,0,0};

  FillHist(Form("m_HitPosition"),mnpos.X(),mnpos.Y());
  if(mnflag){
    FillHist(Form("mn_HitPosition"),mnpos.X(),mnpos.Y());
    if(k0flag){
      FillHist(Form("mn_k0_HitPosition"),mnpos.X(),mnpos.Y());
      //FillHist(Form("mn_k0_Counts"),0);
      for(int i=0; i<5; i++){
        double posll = -40.0 - 10.0*(double)i;
        double posul =  40.0 + 10.0*(double)i;
        if((posll*2.0<mnpos.X()&&mnpos.X()<posul*2.0)&&(posll<mnpos.Y()&&mnpos.Y()<posul)){
          mnposflag[i] = 1;
          FillHist(Form("mn_k0_A%d_HitPosition",i+1),mnpos.X(),mnpos.Y());
          FillHist(Form("mn_k0_Counts"),i+1);
        }
      }
    }
    if(sblk0flag||sbuk0flag){
      FillHist(Form("mn_sbk0_HitPosition"),mnpos.X(),mnpos.Y());
      //FillHist(Form("mn_k0_Counts"),0);
      for(int i=0; i<5; i++){
        double posll = -40.0 - 10.0*(double)i;
        double posul =  40.0 + 10.0*(double)i;
        if((posll<mnpos.X()&&mnpos.X()<posul)&&(posll<mnpos.Y()&&mnpos.Y()<posul)){
          mnposflag[i] = 1;
          FillHist(Form("mn_sbk0_A%d_HitPosition",i+1),mnpos.X(),mnpos.Y());
          FillHist(Form("mn_sbk0_Counts"),i+1);
        }
      }
    }
  }

  pNC*  nc =0;
  if(fnc!=-1)  nc =particle->nc(fnc);
  /* with forward neutron */
  if(nc==0) return true;
  //std::cout << "###!!! NC Hit !!!###" << std::endl;
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
  double npipi_mass[2] = { (Lpro+Ln).M()                           , (cLpro+cLn).M()};
  double npipi_mom[2]  = { (Lpro+Ln).Vect().Mag()                  , (cLpro+cLn).Vect().Mag()};
  double npipi_cost[2] = { (Lpro+Ln).Vect().CosTheta()             , (cLpro+cLn).Vect().CosTheta()};
  double npipi_phi[2]  = { (Lpro+Ln).Vect().Phi()                  , (cLpro+cLn).Vect().Phi()};
  double seg = (nc->seg()-1)%16+1;
  double lay = (nc->seg()-1)/16+1;

  bool spflag = false; /* SigmaPlus flag */
  if(spll<npip_mass[1]&&npip_mass[1]<spul) spflag = true;
  bool smflag = false; /* SigmaMunus flag */
  if(smll<npim_mass[1]&&npim_mass[1]<smul) smflag = true;
  bool mk0flag = false; /* SigmaMunus flag */
  if(mk0ll<nc_mmass[1]&&nc_mmass[1]<mk0ul) mk0flag = true;

  FillHist("FWD_n_OverbetavsMomentum",1.0/nc->beta(),nc->mom().Mag());
  FillHist("FWD_n_OverbetavsEnergy",1.0/nc->beta(),nc->energy());
  FillHist("FWD_n_TOFvsMomentum",nc->tof(),nc->mom().Mag());
  FillHist("FWD_n_HitPosition",nc->hitpos().X(),nc->hitpos().Y());
  FillHist("FWD_n_HitSegment",seg,lay);
  for(int f=1; f<2; f++){
    FillHist(Form("n_CosTvsMomentum_%s",frame[f].Data()),nc_cost[f],nc_mom[f]);
    FillHist(Form("n_CosTvsPhi_%s",frame[f].Data()),nc_cost[f],nc_phi[f]);
    FillHist(Form("pKn_MMvsMomentum_%s",frame[f].Data()),nc_mmass[f],nc_mmom[f]);
    FillHist(Form("pKn_MMvsCosT_%s",frame[f].Data()),nc_mmass[f],nc_mcost[f]);
  }
  for(int f=1; f<2; f++){
    FillHist(Form("pipi_n_IMvsMomentum_%s",frame[f].Data()),pro_mass[f],pro_mom[f]);
    FillHist(Form("pipi_n_IMvsCosT_%s",frame[f].Data()),pro_mass[f],pro_cost[f]);
    FillHist(Form("pipi_n_IMvsDCA_%s",frame[f].Data()),pro_mass[f],pro_pbdca);
    FillHist(Form("pKpipi_n_MMvsMomentum_%s",frame[f].Data()),pro_mmass[f],pro_mmom[f]);
    FillHist(Form("pKpipi_n_MMvsCosT_%s",frame[f].Data()),pro_mmass[f],pro_mcost[f]);
  }
  for(int f=1; f<2; f++){
    FillHist(Form("npip_IMvsMomentum_%s",frame[f].Data()),npip_mass[f],npip_mom[f]);
    FillHist(Form("npip_IMvsCosT_%s",frame[f].Data()),npip_mass[f],npip_cost[f]);
    FillHist(Form("npim_IMvsMomentum_%s",frame[f].Data()),npim_mass[f],npim_mom[f]);
    FillHist(Form("npim_IMvsCosT_%s",frame[f].Data()),npim_mass[f],npim_cost[f]);
    FillHist(Form("npipi_IMvsMomentum_%s",frame[f].Data()),npipi_mass[f],npipi_mom[f]);
    FillHist(Form("npipi_IMvsCosT_%s",frame[f].Data()),npipi_mass[f],npipi_cost[f]);
    FillHist(Form("pKnpip_MMvsMomentum_%s",frame[f].Data()),npip_mmass[f],npip_mmom[f]);
    FillHist(Form("pKnpip_MMvsCosT_%s",frame[f].Data()),npip_mmass[f],npip_mcost[f]);
    FillHist(Form("pKnpim_MMvsMomentum_%s",frame[f].Data()),npim_mmass[f],npim_mmom[f]);
    FillHist(Form("pKnpim_MMvsCosT_%s",frame[f].Data()),npim_mmass[f],npim_mcost[f]);
    FillHist(Form("npipvsnpim_IM_%s",frame[f].Data()),npip_mass[f],npim_mass[f]);
    FillHist(Form("pipivspip_IM_%s",frame[f].Data()),pro_mass[f],npip_mass[f]);
    FillHist(Form("pipivspim_IM_%s",frame[f].Data()),pro_mass[f],npim_mass[f]);
  }
  /* with k0 flag */
  if(k0flag){
    for(int f=1; f<2; f++){
      FillHist(Form("n_k0_CosTvsMomentum_%s",frame[f].Data()),nc_cost[f],nc_mom[f]);
      FillHist(Form("n_k0_CosTvsPhi_%s",frame[f].Data()),nc_cost[f],nc_phi[f]);
      FillHist(Form("pKn_k0_MMvsMomentum_%s",frame[f].Data()),nc_mmass[f],nc_mmom[f]);
      FillHist(Form("pKn_k0_MMvsCosT_%s",frame[f].Data()),nc_mmass[f],nc_mcost[f]);
    }
    for(int f=1; f<2; f++){
      FillHist(Form("pipi_nk0_IMvsMomentum_%s",frame[f].Data()),pro_mass[f],pro_mom[f]);
      FillHist(Form("pipi_nk0_IMvsCosT_%s",frame[f].Data()),pro_mass[f],pro_cost[f]);
      FillHist(Form("pipi_nk0_IMvsDCA_%s",frame[f].Data()),pro_mass[f],pro_pbdca);
      FillHist(Form("pKpipi_nk0_MMvsMomentum_%s",frame[f].Data()),pro_mmass[f],pro_mmom[f]);
      FillHist(Form("pKpipi_nk0_MMvsCosT_%s",frame[f].Data()),pro_mmass[f],pro_mcost[f]);
    }
    for(int f=1; f<2; f++){
      FillHist(Form("npip_k0_IMvsMomentum_%s",frame[f].Data()),npip_mass[f],npip_mom[f]);
      FillHist(Form("npip_k0_IMvsCosT_%s",frame[f].Data()),npip_mass[f],npip_cost[f]);
      FillHist(Form("npim_k0_IMvsMomentum_%s",frame[f].Data()),npim_mass[f],npim_mom[f]);
      FillHist(Form("npim_k0_IMvsCosT_%s",frame[f].Data()),npim_mass[f],npim_cost[f]);
      FillHist(Form("npipi_k0_IMvsMomentum_%s",frame[f].Data()),npipi_mass[f],npipi_mom[f]);
      FillHist(Form("npipi_k0_IMvsCosT_%s",frame[f].Data()),npipi_mass[f],npipi_cost[f]);
      FillHist(Form("pKnpip_k0_MMvsMomentum_%s",frame[f].Data()),npip_mmass[f],npip_mmom[f]);
      FillHist(Form("pKnpip_k0_MMvsCosT_%s",frame[f].Data()),npip_mmass[f],npip_mcost[f]);
      FillHist(Form("pKnpim_k0_MMvsMomentum_%s",frame[f].Data()),npim_mmass[f],npim_mmom[f]);
      FillHist(Form("pKnpim_k0_MMvsCosT_%s",frame[f].Data()),npim_mmass[f],npim_mcost[f]);
    }
    FillHist("FWD_k0_OverbetavsMomentum",1.0/nc->beta(),nc->mom().Mag());
    FillHist("FWD_k0_OverbetavsEnergy",1.0/nc->beta(),nc->energy());
    FillHist("FWD_k0_TOFvsMomentum",nc->tof(),nc->mom().Mag());
    FillHist("FWD_k0_HitPosition",nc->hitpos().X(),nc->hitpos().Y());
    FillHist("FWD_k0_HitSegment",seg,lay);
  }
  /* with L-sideband k0 flag */
  if(sblk0flag){
    for(int f=1; f<2; f++){
      FillHist(Form("pKn_sblk0_MMvsMomentum_%s",frame[f].Data()),nc_mmass[f],nc_mmom[f]);
      FillHist(Form("pKn_sblk0_MMvsCosT_%s",frame[f].Data()),nc_mmass[f],nc_mcost[f]);
      FillHist(Form("npipi_sblk0_IMvsMomentum_%s",frame[f].Data()),npipi_mass[f],npipi_mom[f]);
      FillHist(Form("npipi_sblk0_IMvsCosT_%s",frame[f].Data()),npipi_mass[f],npipi_cost[f]);
    }
  }
  /* with U-sideband k0 flag */
  if(sbuk0flag){
    for(int f=1; f<2; f++){
      FillHist(Form("pKn_sbuk0_MMvsMomentum_%s",frame[f].Data()),nc_mmass[f],nc_mmom[f]);
      FillHist(Form("pKn_sbuk0_MMvsCosT_%s",frame[f].Data()),nc_mmass[f],nc_mcost[f]);
      FillHist(Form("npipi_sbuk0_IMvsMomentum_%s",frame[f].Data()),npipi_mass[f],npipi_mom[f]);
      FillHist(Form("npipi_sbuk0_IMvsCosT_%s",frame[f].Data()),npipi_mass[f],npipi_cost[f]);
    }
  }
  /* with sideband k0 flag */
  if(sblk0flag||sbuk0flag){
    for(int f=1; f<2; f++){
      FillHist(Form("pKn_sbk0_MMvsMomentum_%s",frame[f].Data()),nc_mmass[f],nc_mmom[f]);
      FillHist(Form("pKn_sbk0_MMvsCosT_%s",frame[f].Data()),nc_mmass[f],nc_mcost[f]);
      FillHist(Form("npipi_sbk0_IMvsMomentum_%s",frame[f].Data()),npipi_mass[f],npipi_mom[f]);
      FillHist(Form("npipi_sbk0_IMvsCosT_%s",frame[f].Data()),npipi_mass[f],npipi_cost[f]);
    }
  }
  /* without k0 flag */
  if(!k0flag){
    for(int f=1; f<2; f++){
      FillHist(Form("pKn_wok0_MMvsMomentum_%s",frame[f].Data()),nc_mmass[f],nc_mmom[f]);
      FillHist(Form("pKn_wok0_MMvsCosT_%s",frame[f].Data()),nc_mmass[f],nc_mcost[f]);
      FillHist(Form("pKnpip_wok0_MMvsMomentum_%s",frame[f].Data()),npip_mmass[f],npip_mmom[f]);
      FillHist(Form("pKnpip_wok0_MMvsCosT_%s",frame[f].Data()),npip_mmass[f],npip_mcost[f]);
      FillHist(Form("pKnpim_wok0_MMvsMomentum_%s",frame[f].Data()),npim_mmass[f],npim_mmom[f]);
      FillHist(Form("pKnpim_wok0_MMvsCosT_%s",frame[f].Data()),npim_mmass[f],npim_mcost[f]);
      FillHist(Form("npipi_wok0_IMvsMomentum_%s",frame[f].Data()),npipi_mass[f],npipi_mom[f]);
      FillHist(Form("npipi_wok0_IMvsCosT_%s",frame[f].Data()),npipi_mass[f],npipi_cost[f]);
    }
  }
  /* with missing neutron flag */
  if(mnflag){
    for(int f=1; f<2; f++){
      FillHist(Form("pKn_mn_MMvsMomentum_%s",frame[f].Data()),nc_mmass[f],nc_mmom[f]);
      FillHist(Form("pKn_mn_MMvsCosT_%s",frame[f].Data()),nc_mmass[f],nc_mcost[f]);
      FillHist(Form("npipi_nmn_IMvsMomentum_%s",frame[f].Data()),npipi_mass[f],npipi_mom[f]);
      FillHist(Form("npipi_nmn_IMvsCosT_%s",frame[f].Data()),npipi_mass[f],npipi_cost[f]);
    }
    if(k0flag){
      for(int f=1; f<2; f++){
        FillHist(Form("pKn_k0mn_MMvsMomentum_%s",frame[f].Data()),nc_mmass[f],nc_mmom[f]);
        FillHist(Form("pKn_k0mn_MMvsCosT_%s",frame[f].Data()),nc_mmass[f],nc_mcost[f]);
        FillHist(Form("npipi_k0mn_IMvsMomentum_%s",frame[f].Data()),npipi_mass[f],npipi_mom[f]);
        FillHist(Form("npipi_k0mn_IMvsCosT_%s",frame[f].Data()),npipi_mass[f],npipi_cost[f]);
      }
      /* without Sigma flag */
      if(!spflag&&!smflag){
        for(int f=1; f<2; f++){
          FillHist(Form("pKn_k0mnwos_MMvsMomentum_%s",frame[f].Data()),nc_mmass[f],nc_mmom[f]);
          FillHist(Form("pKn_k0mnwos_MMvsCosT_%s",frame[f].Data()),nc_mmass[f],nc_mcost[f]);
          FillHist(Form("npipi_k0mnwos_IMvsMomentum_%s",frame[f].Data()),npipi_mass[f],npipi_mom[f]);
          FillHist(Form("npipi_k0mnwos_IMvsCosT_%s",frame[f].Data()),npipi_mass[f],npipi_cost[f]);
        }
      }
    }
    if(sblk0flag||sbuk0flag){
      for(int f=1; f<2; f++){
        FillHist(Form("pKn_sbk0mn_MMvsMomentum_%s",frame[f].Data()),nc_mmass[f],nc_mmom[f]);
        FillHist(Form("pKn_sbk0mn_MMvsCosT_%s",frame[f].Data()),nc_mmass[f],nc_mcost[f]);
        FillHist(Form("npipi_sbk0mn_IMvsMomentum_%s",frame[f].Data()),npipi_mass[f],npipi_mom[f]);
        FillHist(Form("npipi_sbk0mn_IMvsCosT_%s",frame[f].Data()),npipi_mass[f],npipi_cost[f]);
      }
    }
  }
  /* without missing neutron flag */
  if(!mnflag){
    for(int f=1; f<2; f++){
      FillHist(Form("pKn_womn_MMvsMomentum_%s",frame[f].Data()),nc_mmass[f],nc_mmom[f]);
      FillHist(Form("pKn_womn_MMvsCosT_%s",frame[f].Data()),nc_mmass[f],nc_mcost[f]);
      FillHist(Form("npipi_nwomn_IMvsMomentum_%s",frame[f].Data()),npipi_mass[f],npipi_mom[f]);
      FillHist(Form("npipi_nwomn_IMvsCosT_%s",frame[f].Data()),npipi_mass[f],npipi_cost[f]);
    }
  }
  if(mnflag){
    FillHist(Form("mn_n_HitPosition"),mnpos.X(),mnpos.Y());
    /* with k0 and mn flag */
    if(k0flag){
      FillHist(Form("mn_nk0_HitPosition"),mnpos.X(),mnpos.Y());
    }
  }
  /* with Sp and mn flag */
  if(spflag){
    for(int f=1; f<2; f++){
      FillHist(Form("pKn_sp_MMvsMomentum_%s",frame[f].Data()),nc_mmass[f],nc_mmom[f]);
      FillHist(Form("pKn_sp_MMvsCosT_%s",frame[f].Data()),nc_mmass[f],nc_mcost[f]);
      FillHist(Form("pipi_sp_IMvsMomentum_%s",frame[f].Data()),pro_mass[f],pro_mom[f]);
      FillHist(Form("pipi_sp_IMvsCosT_%s",frame[f].Data()),pro_mass[f],pro_cost[f]);
      FillHist(Form("pKpipi_sp_MMvsMomentum_%s",frame[f].Data()),pro_mmass[f],pro_mmom[f]);
      FillHist(Form("pKpipi_sp_MMvsCosT_%s",frame[f].Data()),pro_mmass[f],pro_mcost[f]);
    }
    /* without k0 flag */
    if(!k0flag){
      for(int f=1; f<2; f++){
        FillHist(Form("pKn_spwok0_MMvsMomentum_%s",frame[f].Data()),nc_mmass[f],nc_mmom[f]);
        FillHist(Form("pKn_spwok0_MMvsCosT_%s",frame[f].Data()),nc_mmass[f],nc_mcost[f]);
        FillHist(Form("pipi_spwok0_IMvsMomentum_%s",frame[f].Data()),pro_mass[f],pro_mom[f]);
        FillHist(Form("pipi_spwok0_IMvsCosT_%s",frame[f].Data()),pro_mass[f],pro_cost[f]);
        FillHist(Form("pKpipi_spwok0_MMvsMomentum_%s",frame[f].Data()),pro_mmass[f],pro_mmom[f]);
        FillHist(Form("pKpipi_spwok0_MMvsCosT_%s",frame[f].Data()),pro_mmass[f],pro_mcost[f]);
      }
    }
    /* without missing neutron flag */
    if(!mnflag){
      for(int f=1; f<2; f++){
        FillHist(Form("pKn_spwomn_MMvsMomentum_%s",frame[f].Data()),nc_mmass[f],nc_mmom[f]);
        FillHist(Form("pKn_spwomn_MMvsCosT_%s",frame[f].Data()),nc_mmass[f],nc_mcost[f]);
        FillHist(Form("pipi_spwomn_IMvsMomentum_%s",frame[f].Data()),pro_mass[f],pro_mom[f]);
        FillHist(Form("pipi_spwomn_IMvsCosT_%s",frame[f].Data()),pro_mass[f],pro_cost[f]);
        FillHist(Form("pKpipi_spwomn_MMvsMomentum_%s",frame[f].Data()),pro_mmass[f],pro_mmom[f]);
        FillHist(Form("pKpipi_spwomn_MMvsCosT_%s",frame[f].Data()),pro_mmass[f],pro_mcost[f]);
      }
    }
  }
  /* without Sp and mn flag */
  if(!spflag){
    for(int f=1; f<2; f++){
      FillHist(Form("pKn_wosp_MMvsMomentum_%s",frame[f].Data()),nc_mmass[f],nc_mmom[f]);
      FillHist(Form("pKn_wosp_MMvsCosT_%s",frame[f].Data()),nc_mmass[f],nc_mcost[f]);
      FillHist(Form("pipi_wosp_IMvsMomentum_%s",frame[f].Data()),pro_mass[f],pro_mom[f]);
      FillHist(Form("pipi_wosp_IMvsCosT_%s",frame[f].Data()),pro_mass[f],pro_cost[f]);
      FillHist(Form("pKpipi_wosp_MMvsMomentum_%s",frame[f].Data()),pro_mmass[f],pro_mmom[f]);
      FillHist(Form("pKpipi_wosp_MMvsCosT_%s",frame[f].Data()),pro_mmass[f],pro_mcost[f]);
    }
  }
  /* with Sm and mn flag */
  if(smflag){
    for(int f=1; f<2; f++){
      FillHist(Form("pKn_sm_MMvsMomentum_%s",frame[f].Data()),nc_mmass[f],nc_mmom[f]);
      FillHist(Form("pKn_sm_MMvsCosT_%s",frame[f].Data()),nc_mmass[f],nc_mcost[f]);
      FillHist(Form("pipi_sm_IMvsMomentum_%s",frame[f].Data()),pro_mass[f],pro_mom[f]);
      FillHist(Form("pipi_sm_IMvsCosT_%s",frame[f].Data()),pro_mass[f],pro_cost[f]);
      FillHist(Form("pKpipi_sm_MMvsMomentum_%s",frame[f].Data()),pro_mmass[f],pro_mmom[f]);
      FillHist(Form("pKpipi_sm_MMvsCosT_%s",frame[f].Data()),pro_mmass[f],pro_mcost[f]);
    }
    /* without k0 flag */
    if(!k0flag){
      for(int f=1; f<2; f++){
        FillHist(Form("pKn_smwok0_MMvsMomentum_%s",frame[f].Data()),nc_mmass[f],nc_mmom[f]);
        FillHist(Form("pKn_smwok0_MMvsCosT_%s",frame[f].Data()),nc_mmass[f],nc_mcost[f]);
        FillHist(Form("pipi_smwok0_IMvsMomentum_%s",frame[f].Data()),pro_mass[f],pro_mom[f]);
        FillHist(Form("pipi_smwok0_IMvsCosT_%s",frame[f].Data()),pro_mass[f],pro_cost[f]);
        FillHist(Form("pKpipi_smwok0_MMvsMomentum_%s",frame[f].Data()),pro_mmass[f],pro_mmom[f]);
        FillHist(Form("pKpipi_smwok0_MMvsCosT_%s",frame[f].Data()),pro_mmass[f],pro_mcost[f]);
      }
    }
    /* without missing neutron flag */
    if(!mnflag){
      for(int f=1; f<2; f++){
        FillHist(Form("pKn_smwomn_MMvsMomentum_%s",frame[f].Data()),nc_mmass[f],nc_mmom[f]);
        FillHist(Form("pKn_smwomn_MMvsCosT_%s",frame[f].Data()),nc_mmass[f],nc_mcost[f]);
        FillHist(Form("pipi_smwomn_IMvsMomentum_%s",frame[f].Data()),pro_mass[f],pro_mom[f]);
        FillHist(Form("pipi_smwomn_IMvsCosT_%s",frame[f].Data()),pro_mass[f],pro_cost[f]);
        FillHist(Form("pKpipi_smwomn_MMvsMomentum_%s",frame[f].Data()),pro_mmass[f],pro_mmom[f]);
        FillHist(Form("pKpipi_smwomn_MMvsCosT_%s",frame[f].Data()),pro_mmass[f],pro_mcost[f]);
      }
    }
  }
  /* without Sm and mn flag */
  if(!smflag){
    for(int f=1; f<2; f++){
      FillHist(Form("pKn_wosm_MMvsMomentum_%s",frame[f].Data()),nc_mmass[f],nc_mmom[f]);
      FillHist(Form("pKn_wosm_MMvsCosT_%s",frame[f].Data()),nc_mmass[f],nc_mcost[f]);
      FillHist(Form("pipi_wosm_IMvsMomentum_%s",frame[f].Data()),pro_mass[f],pro_mom[f]);
      FillHist(Form("pipi_wosm_IMvsCosT_%s",frame[f].Data()),pro_mass[f],pro_cost[f]);
      FillHist(Form("pKpipi_wosm_MMvsMomentum_%s",frame[f].Data()),pro_mmass[f],pro_mmom[f]);
      FillHist(Form("pKpipi_wosm_MMvsCosT_%s",frame[f].Data()),pro_mmass[f],pro_mcost[f]);
    }
  }

/* NC HitPosition */
  int nposflag[5] = {0,0,0,0,0};

  if(mnflag){
    FillHist(Form("mn_n_HitPosition"),mnpos.X(),mnpos.Y());
    if(k0flag){
      FillHist(Form("mn_nk0_HitPosition"),mnpos.X(),mnpos.Y());
      //FillHist(Form("mn_k0_Counts"),0);
      if(!spflag&&!smflag){
        FillHist(Form("mn_nk0wos_HitPosition"),mnpos.X(),mnpos.Y());
      }
      for(int i=0; i<5; i++){
        double posll = -40.0 - 10.0*(double)i;
        double posul =  40.0 + 10.0*(double)i;
        if((posll*2.0<nc->hitpos().X()&&nc->hitpos().X()<posul*2.0)&&(posll<nc->hitpos().Y()&&nc->hitpos().Y()<posul)){
          nposflag[i] = 1;
          FillHist(Form("mn_nk0_A%d_HitPosition",i+1),mnpos.X(),mnpos.Y());
          FillHist(Form("mn_nk0_Counts"),i+1);
          if((-160.0<mnpos.X()&&mnpos.X()<160.0)&&(-80.0<mnpos.Y()&&mnpos.Y()<80.0)){
            FillHist(Form("mn_nk0_A%d_IRHitPosition",i+1),mnpos.X(),mnpos.Y());
            FillHist(Form("mn_nk0_IRCounts"),i+1);
          }
          else{
            FillHist(Form("mn_nk0_A%d_ORHitPosition",i+1),mnpos.X(),mnpos.Y());
            FillHist(Form("mn_nk0_ORCounts"),i+1);
          }
        }
      }
    }
    if(sblk0flag||sbuk0flag){
      FillHist(Form("mn_nsbk0_HitPosition"),mnpos.X(),mnpos.Y());
      //FillHist(Form("mn_k0_Counts"),0);
      for(int i=0; i<5; i++){
        double posll = -40.0 - 10.0*(double)i;
        double posul =  40.0 + 10.0*(double)i;
        if(nposflag[i]){
          FillHist(Form("mn_nsbk0_A%d_HitPosition",i+1),mnpos.X(),mnpos.Y());
          FillHist(Form("mn_nsbk0_Counts"),i+1);
          if((-160.0<mnpos.X()&&mnpos.X()<160.0)&&(-80.0<mnpos.Y()&&mnpos.Y()<80.0)){
            FillHist(Form("mn_nsbk0_A%d_IRHitPosition",i+1),mnpos.X(),mnpos.Y());
            FillHist(Form("mn_nsbk0_IRCounts"),i+1);
          }
          else{
            FillHist(Form("mn_nsbk0_A%d_ORHitPosition",i+1),mnpos.X(),mnpos.Y());
            FillHist(Form("mn_nsbk0_ORCounts"),i+1);
          }
        }
      }
    }
  }


  FillHist("n_HitPosition",nc->hitpos().X(),nc->hitpos().Y());
  if(k0flag){
    FillHist(Form("n_k0_HitPosition"),nc->hitpos().X(),nc->hitpos().Y());
  }
  if(mnflag){
    FillHist(Form("n_mn_HitPosition"),nc->hitpos().X(),nc->hitpos().Y());
    if(k0flag){
      FillHist(Form("n_k0mn_HitPosition"),nc->hitpos().X(),nc->hitpos().Y());
      if(!spflag&&!smflag){
        FillHist(Form("n_k0mnwos_HitPosition"),nc->hitpos().X(),nc->hitpos().Y());
        FillHist(Form("n_k0mnwos_HitPosDiff"),nc->hitpos().X()-mnpos.X(),nc->hitpos().Y()-mnpos.Y());
        for(int i=0; i<5; i++){
          if(mnposflag[i]==1){
            if(mk0flag){
              FillHist(Form("n_k0mnwos_A%d_HitPosition",i+1),nc->hitpos().X(),nc->hitpos().Y());
              FillHist(Form("n_k0mnwos_A%d_HitPosDiff",i+1),nc->hitpos().X()-mnpos.X(),nc->hitpos().Y()-mnpos.Y());
              FillHist(Form("n_k0mnwos_Counts"),i+1);
            }
          }
        } 
      }
    }
    if(sblk0flag&&sbuk0flag){
      FillHist(Form("n_sbk0mn_HitPosition"),nc->hitpos().X(),nc->hitpos().Y());
      if(!spflag&&!smflag){
        FillHist(Form("n_sbk0mnwos_HitPosition"),nc->hitpos().X(),nc->hitpos().Y());
        FillHist(Form("n_sbk0mnwos_HitPosDiff"),nc->hitpos().X()-mnpos.X(),nc->hitpos().Y()-mnpos.Y());
        for(int i=0; i<5; i++){
          if(mnposflag[i]==1){
            if(mk0flag){
              FillHist(Form("n_sbk0mnwos_A%d_HitPosition",i+1),nc->hitpos().X(),nc->hitpos().Y());
              FillHist(Form("n_sbk0mnwos_A%d_HitPosDiff",i+1),nc->hitpos().X()-mnpos.X(),nc->hitpos().Y()-mnpos.Y());
              FillHist(Form("n_sbk0mnwos_Counts"),i+1);
            }
          } 
        } 
      }
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
  /* Missing k0 mass cut */ 
  mean  = 0.4970;
  sigma = 0.0130;
  mk0ll = mean-2*sigma; mk0ul = mean+2*sigma;
  sblmk0ll = mean-5*sigma; sblmk0ul = mean-3*sigma;
  sbumk0ll = mean+3*sigma; sbumk0ul = mean+5*sigma;
}

bool MyAnalysisK0::Initialize(ConfMan* confMan)
{
  std::cout << "### MyAnalysisK0::Initialize ###" << std::endl;

  std::string ofname = confMan->GetOutFileName();
  ofname.insert(ofname.find(".root"),"_anaK0");

  rtFile =  new TFile( Form("%s",ofname.c_str()), "RECREATE");
  rtFile -> cd();

  // CDS Particles
  std::cout << "Define Histograms for CDS particles" << std::endl;
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
  for(int i=1; i<=18; i++){
    new TH1F( Form("Trigger%02d_TDC",i), Form("Trigger%02d TDC;TDC ch.;Counts",i), 4096, 0, 4096);
  }

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
  new TH2F("pipMomentumvspimMomentum_CM","#pi^{+} momentum vs. #pi^{-} momentum;#pi^{+} momentum (GeV/c);#pi^{-} momentum (GeV/c)", 1000, 0.0, 1.0, 1000, 0.0, 1.0);
  new TH2F("pipCosTvspimCosT_CM","#pi^{+} cos#theta vs. #pi^{-} cos#theta;cos#theta_{#pi^{+}};cos#theta_{#pi^{-}}", 2000, -1.0, 1.0, 2000, -1.0, 1.0);
  new TH2F("pipPhivspimPhi_CM","#pi^{+} #phi vs. #pi^{-} #phi;#phi_{#pi^{+}};#phi_{#pi^{-}}", 1600, 0.0, 3.2, 1600, 0.0, 3.2);
  new TH2F("pip_CosTvsMomentum_CM","cos#theta vs. momentum of #pi^{+};cos#theta;Momentum (GeV/c)", 2000, -1.0, 1.0, 1000, 0.0, 1.0);
  new TH2F("pip_CosTvsPhi_CM","cos#theta vs. #phi of #pi^{+};cos#theta;#phi", 2000, -1.0, 1.0, 1600, 0.0, 3.2);
  new TH2F("pim_CosTvsMomentum_CM","cos#theta vs. momentum of #pi^{-};cos#theta;Momentum (GeV/c)", 2000, -1.0, 1.0, 1000, 0.0, 1.0);
  new TH2F("pim_CosTvsPhi_CM","cos#theta vs. #phi of #pi^{-};cos#theta;#phi", 2000, -1.0, 1.0, 1600, 0.0, 3.2);
  new TH2F("pipi_CosTvsMomentum_CM","cos#theta vs. momentum of Production;cos#theta;Momentum (GeV/c)", 2000, -1.0, 1.0, 1000, 0.0, 1.0);
  new TH2F("pipi_CosTvsPhi_CM","cos#theta vs. #phi of Production;cos#theta;#phi", 2000, -1.0, 1.0, 1600, 0.0, 3.2);
  new TH1F("pip_PBDCA","DCA between beam #pi^{+};DCA (cm);Counts",100,0.0,10.0);
  new TH1F("pim_PBDCA","DCA between beam #pi^{-};DCA (cm);Counts",100,0.0,10.0);
  new TH1F("Pro_VDIS","Distance between tow vertexes;Distance (cm);Counts",100,0.0,10.0);
  new TH1F("Pro_VBDIS","distance between beam and vertex;Distance (cm);Counts",100,0.0,10.0);
  new TH1F("Pro_PBDCA","DCA between beam and product;DCA (cm);Counts",100,0.0,10.0);
  new TH2F("n_CosTvsMomentum_CM","cos#theta vs. momentum of n;cos#theta;Momentum (GeV/c)", 2000, -1.0, 1.0, 1000, 0.0, 1.0);
  new TH2F("n_CosTvsPhi_CM","cos#theta vs. #phi of n;cos#theta;#phi", 2000, -1.0, 1.0, 1600, 0.0, 3.2);
  // FWD neutral Particles
  std::cout << "Define Histograms for FWD neutral particles" << std::endl;
  new TH1F("mn_k0_Counts","Missing neutron counts;Selected region;Counts",6,-0.5,5.5);
  new TH1F("mn_sbk0_Counts","Missing neutron counts;Selected region;Counts",6,-0.5,5.5);
  new TH1F("mn_nk0_Counts","Missing neutron counts;Selected region;Counts",6,-0.5,5.5);
  new TH1F("mn_nsbk0_Counts","Missing neutron counts;Selected region;Counts",6,-0.5,5.5);
  new TH1F("n_k0mnwos_Counts","Neutron counts;Selected region;Counts",6,-0.5,5.5);
  new TH1F("n_sbk0mnwos_Counts","Neutron counts;Selected region;Counts",6,-0.5,5.5);
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
  new TH2F( "m_HitPosition", "Missing neutron position at NC (X vs. Y);X position (cm);Y position (cm)",40,-400,400,25,-250,250);
  new TH2F( "mn_HitPosition", "Missing neutron position at NC (X vs. Y);X position (cm);Y position (cm)",40,-400,400,25,-250,250);
  new TH2F( "mn_k0_HitPosition", "Missing neutron position at NC (X vs. Y);X position (cm);Y position (cm)",40,-400,400,25,-250,250);
  new TH2F( "mn_k0_A1_HitPosition", "Missing neutron position at NC (X vs. Y);X position (cm);Y position (cm)",40,-400,400,25,-250,250);
  new TH2F( "mn_k0_A2_HitPosition", "Missing neutron position at NC (X vs. Y);X position (cm);Y position (cm)",40,-400,400,25,-250,250);
  new TH2F( "mn_k0_A3_HitPosition", "Missing neutron position at NC (X vs. Y);X position (cm);Y position (cm)",40,-400,400,25,-250,250);
  new TH2F( "mn_k0_A4_HitPosition", "Missing neutron position at NC (X vs. Y);X position (cm);Y position (cm)",40,-400,400,25,-250,250);
  new TH2F( "mn_k0_A5_HitPosition", "Missing neutron position at NC (X vs. Y);X position (cm);Y position (cm)",40,-400,400,25,-250,250);
  new TH2F( "mn_n_HitPosition", "Missing neutron position at NC (X vs. Y);X position (cm);Y position (cm)",40,-400,400,25,-250,250);
  new TH2F( "mn_nk0_A1_HitPosition", "Missing neutron position at NC (X vs. Y);X position (cm);Y position (cm)",40,-400,400,25,-250,250);
  new TH2F( "mn_nk0_A2_HitPosition", "Missing neutron position at NC (X vs. Y);X position (cm);Y position (cm)",40,-400,400,25,-250,250);
  new TH2F( "mn_nk0_A3_HitPosition", "Missing neutron position at NC (X vs. Y);X position (cm);Y position (cm)",40,-400,400,25,-250,250);
  new TH2F( "mn_nk0_A4_HitPosition", "Missing neutron position at NC (X vs. Y);X position (cm);Y position (cm)",40,-400,400,25,-250,250);
  new TH2F( "mn_nk0_A5_HitPosition", "Missing neutron position at NC (X vs. Y);X position (cm);Y position (cm)",40,-400,400,25,-250,250);
  new TH2F( "mn_nk0_HitPosition", "Missing neutron position at NC (X vs. Y);X position (cm);Y position (cm)",40,-400,400,25,-250,250);
  new TH2F( "mn_nk0wos_HitPosition", "Missing neutron position at NC (X vs. Y);X position (cm);Y position (cm)",40,-400,400,25,-250,250);
  new TH2F( "n_HitPosition", "Neutron position at NC (X vs. Y);X position (cm);Y position (cm)",40,-400,400,25,-250,250);
  new TH2F( "n_k0_HitPosition", "Neutron position at NC (X vs. Y);X position (cm);Y position (cm)",40,-400,400,25,-250,250);
  new TH2F( "n_mn_HitPosition", "Neutron position at NC (X vs. Y);X position (cm);Y position (cm)",40,-400,400,25,-250,250);
  new TH2F( "n_k0mn_HitPosition", "Neutron position at NC (X vs. Y);X position (cm);Y position (cm)",40,-400,400,25,-250,250);
  new TH2F( "n_k0mnwos_HitPosition", "Neutron position at NC (X vs. Y);X position (cm);Y position (cm)",40,-400,400,25,-250,250);
  new TH2F( "n_k0mnwos_HitPosDiff", "Neutron position difference btw NC and Missng n (X vs. Y);X position (cm);Y position (cm)",40,-400,400,25,-250,250);
  new TH2F( "n_k0mnwos_A1_HitPosition", "Neutron position at NC (X vs. Y);X position (cm);Y position (cm)",40,-400,400,25,-250,250);
  new TH2F( "n_k0mnwos_A1_HitPosDiff", "Neutron position difference btw NC and Missng n (X vs. Y);X position (cm);Y position (cm)",40,-400,400,25,-250,250);
  new TH2F( "n_k0mnwos_A2_HitPosition", "Neutron position at NC (X vs. Y);X position (cm);Y position (cm)",40,-400,400,25,-250,250);
  new TH2F( "n_k0mnwos_A2_HitPosDiff", "Neutron position difference btw NC and Missng n (X vs. Y);X position (cm);Y position (cm)",40,-400,400,25,-250,250);
  new TH2F( "n_k0mnwos_A3_HitPosition", "Neutron position at NC (X vs. Y);X position (cm);Y position (cm)",40,-400,400,25,-250,250);
  new TH2F( "n_k0mnwos_A3_HitPosDiff", "Neutron position difference btw NC and Missng n (X vs. Y);X position (cm);Y position (cm)",40,-400,400,25,-250,250);
  new TH2F( "n_k0mnwos_A4_HitPosition", "Neutron position at NC (X vs. Y);X position (cm);Y position (cm)",40,-400,400,25,-250,250);
  new TH2F( "n_k0mnwos_A4_HitPosDiff", "Neutron position difference btw NC and Missng n (X vs. Y);X position (cm);Y position (cm)",40,-400,400,25,-250,250);
  new TH2F( "n_k0mnwos_A5_HitPosition", "Neutron position at NC (X vs. Y);X position (cm);Y position (cm)",40,-400,400,25,-250,250);
  new TH2F( "n_k0mnwos_A5_HitPosDiff", "Neutron position difference btw NC and Missng n (X vs. Y);X position (cm);Y position (cm)",40,-400,400,25,-250,250);
  // Invariant mass of production
  std::cout << "Define Histograms for product I.M. with FWD" << std::endl;
  new TH2F( "pipi_IMvsMomentum_CM", "Invariant mass (#pi^{+}#pi^{-}) vs Momentum;IM(#pi^{+}#pi^{-}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0 );
  new TH2F( "pipi_IMvsCosT_CM", "Invariant mass (#pi^{+}#pi^{-}) vs cos#theta;IM(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0 );
  new TH2F( "pipi_IMvsDCA_CM", "Invariant mass (#pi^{+}#pi^{-}) vs DCA;IM(#pi^{+}#pi^{-}) (GeV/c^{2});DCA (cm)", 2000, 0.0, 2.0, 100, 0.0, 10.0 );
  new TH2F( "pipi_k0_IMvsMomentum_CM", "Invariant mass (#pi^{+}#pi^{-}) vs Momentum;IM(#pi^{+}#pi^{-}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0 );
  new TH2F( "pipi_k0_IMvsCosT_CM", "Invariant mass (#pi^{+}#pi^{-}) vs cos#theta;IM(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0 );
  new TH2F( "pipi_k0_IMvsDCA_CM", "Invariant mass (#pi^{+}#pi^{-}) vs DCA;IM(#pi^{+}#pi^{-}) (GeV/c^{2});DCA (cm)", 2000, 0.0, 2.0, 100, 0.0, 10.0 );
  new TH2F( "pipi_sblk0_IMvsMomentum_CM", "Invariant mass (#pi^{+}#pi^{-}) vs Momentum;IM(#pi^{+}#pi^{-}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0 );
  new TH2F( "pipi_sblk0_IMvsCosT_CM", "Invariant mass (#pi^{+}#pi^{-}) vs cos#theta;IM(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0 );
  new TH2F( "pipi_sblk0_IMvsDCA_CM", "Invariant mass (#pi^{+}#pi^{-}) vs DCA;IM(#pi^{+}#pi^{-}) (GeV/c^{2});DCA (cm)", 2000, 0.0, 2.0, 100, 0.0, 10.0 );
  new TH2F( "pipi_sbuk0_IMvsMomentum_CM", "Invariant mass (#pi^{+}#pi^{-}) vs Momentum;IM(#pi^{+}#pi^{-}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0 );
  new TH2F( "pipi_sbuk0_IMvsCosT_CM", "Invariant mass (#pi^{+}#pi^{-}) vs cos#theta;IM(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0 );
  new TH2F( "pipi_sbuk0_IMvsDCA_CM", "Invariant mass (#pi^{+}#pi^{-}) vs DCA;IM(#pi^{+}#pi^{-}) (GeV/c^{2});DCA (cm)", 2000, 0.0, 2.0, 100, 0.0, 10.0 );
  new TH2F( "pipi_sbk0_IMvsMomentum_CM", "Invariant mass (#pi^{+}#pi^{-}) vs Momentum;IM(#pi^{+}#pi^{-}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0 );
  new TH2F( "pipi_sbk0_IMvsCosT_CM", "Invariant mass (#pi^{+}#pi^{-}) vs cos#theta;IM(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0 );
  new TH2F( "pipi_sbk0_IMvsDCA_CM", "Invariant mass (#pi^{+}#pi^{-}) vs DCA;IM(#pi^{+}#pi^{-}) (GeV/c^{2});DCA (cm)", 2000, 0.0, 2.0, 100, 0.0, 10.0 );
  new TH2F( "pipi_mn_IMvsMomentum_CM", "Invariant mass (#pi^{+}#pi^{-}) vs Momentum;IM(#pi^{+}#pi^{-}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0 );
  new TH2F( "pipi_mn_IMvsCosT_CM", "Invariant mass (#pi^{+}#pi^{-}) vs cos#theta;IM(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0 );
  new TH2F( "pipi_mn_IMvsDCA_CM", "Invariant mass (#pi^{+}#pi^{-}) vs DCA;IM(#pi^{+}#pi^{-}) (GeV/c^{2});DCA (cm)", 2000, 0.0, 2.0, 100, 0.0, 10.0 );
  new TH2F( "pipi_womn_IMvsMomentum_CM", "Invariant mass (#pi^{+}#pi^{-}) vs Momentum;IM(#pi^{+}#pi^{-}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0 );
  new TH2F( "pipi_womn_IMvsCosT_CM", "Invariant mass (#pi^{+}#pi^{-}) vs cos#theta;IM(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0 );
  new TH2F( "pipi_womn_IMvsDCA_CM", "Invariant mass (#pi^{+}#pi^{-}) vs DCA;IM(#pi^{+}#pi^{-}) (GeV/c^{2});DCA (cm)", 2000, 0.0, 2.0, 100, 0.0, 10.0 );
  new TH2F( "pipi_k0mn_IMvsMomentum_CM", "Invariant mass (#pi^{+}#pi^{-}) vs Momentum;IM(#pi^{+}#pi^{-}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0 );
  new TH2F( "pipi_k0mn_IMvsCosT_CM", "Invariant mass (#pi^{+}#pi^{-}) vs cos#theta;IM(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0 );
  new TH2F( "pipi_k0mn_IMvsDCA_CM", "Invariant mass (#pi^{+}#pi^{-}) vs DCA;IM(#pi^{+}#pi^{-}) (GeV/c^{2});DCA (cm)", 2000, 0.0, 2.0, 100, 0.0, 10.0 );
  new TH2F( "pipi_sblk0mn_IMvsMomentum_CM", "Invariant mass (#pi^{+}#pi^{-}) vs Momentum;IM(#pi^{+}#pi^{-}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0 );
  new TH2F( "pipi_sblk0mn_IMvsCosT_CM", "Invariant mass (#pi^{+}#pi^{-}) vs cos#theta;IM(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0 );
  new TH2F( "pipi_sblk0mn_IMvsDCA_CM", "Invariant mass (#pi^{+}#pi^{-}) vs DCA;IM(#pi^{+}#pi^{-}) (GeV/c^{2});DCA (cm)", 2000, 0.0, 2.0, 100, 0.0, 10.0 );
  new TH2F( "pipi_sbuk0mn_IMvsMomentum_CM", "Invariant mass (#pi^{+}#pi^{-}) vs Momentum;IM(#pi^{+}#pi^{-}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0 );
  new TH2F( "pipi_sbuk0mn_IMvsCosT_CM", "Invariant mass (#pi^{+}#pi^{-}) vs cos#theta;IM(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0 );
  new TH2F( "pipi_sbuk0mn_IMvsDCA_CM", "Invariant mass (#pi^{+}#pi^{-}) vs DCA;IM(#pi^{+}#pi^{-}) (GeV/c^{2});DCA (cm)", 2000, 0.0, 2.0, 100, 0.0, 10.0 );
  new TH2F( "pipi_sbk0mn_IMvsMomentum_CM", "Invariant mass (#pi^{+}#pi^{-}) vs Momentum;IM(#pi^{+}#pi^{-}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0 );
  new TH2F( "pipi_sbk0mn_IMvsCosT_CM", "Invariant mass (#pi^{+}#pi^{-}) vs cos#theta;IM(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0 );
  new TH2F( "pipi_sbk0mn_IMvsDCA_CM", "Invariant mass (#pi^{+}#pi^{-}) vs DCA;IM(#pi^{+}#pi^{-}) (GeV/c^{2});DCA (cm)", 2000, 0.0, 2.0, 100, 0.0, 10.0 );
  new TH2F( "pipi_k0womn_IMvsMomentum_CM", "Invariant mass (#pi^{+}#pi^{-}) vs Momentum;IM(#pi^{+}#pi^{-}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0 );
  new TH2F( "pipi_k0womn_IMvsCosT_CM", "Invariant mass (#pi^{+}#pi^{-}) vs cos#theta;IM(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0 );
  new TH2F( "pipi_k0womn_IMvsDCA_CM", "Invariant mass (#pi^{+}#pi^{-}) vs DCA;IM(#pi^{+}#pi^{-}) (GeV/c^{2});DCA (cm)", 2000, 0.0, 2.0, 100, 0.0, 10.0 );
  new TH2F( "pipi_n_IMvsMomentum_CM", "Invariant mass (#pi^{+}#pi^{-}) vs Momentum;IM(#pi^{+}#pi^{-}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0 );
  new TH2F( "pipi_n_IMvsCosT_CM", "Invariant mass (#pi^{+}#pi^{-}) vs cos#theta;IM(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0 );
  new TH2F( "pipi_n_IMvsDCA_CM", "Invariant mass (#pi^{+}#pi^{-}) vs DCA;IM(#pi^{+}#pi^{-}) (GeV/c^{2});DCA (cm)", 2000, 0.0, 2.0, 100, 0.0, 10.0 );
  new TH2F( "pipi_nk0_IMvsMomentum_CM", "Invariant mass (#pi^{+}#pi^{-}) vs Momentum;IM(#pi^{+}#pi^{-}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0 );
  new TH2F( "pipi_nk0_IMvsCosT_CM", "Invariant mass (#pi^{+}#pi^{-}) vs cos#theta;IM(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0 );
  new TH2F( "pipi_nk0_IMvsDCA_CM", "Invariant mass (#pi^{+}#pi^{-}) vs DCA;IM(#pi^{+}#pi^{-}) (GeV/c^{2});DCA (cm)", 2000, 0.0, 2.0, 100, 0.0, 10.0 );
  new TH2F( "pipi_nwok0_IMvsMomentum_CM", "Invariant mass (#pi^{+}#pi^{-}) vs Momentum;IM(#pi^{+}#pi^{-}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0 );
  new TH2F( "pipi_nwok0_IMvsCosT_CM", "Invariant mass (#pi^{+}#pi^{-}) vs cos#theta;IM(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0 );
  new TH2F( "pipi_nwok0_IMvsDCA_CM", "Invariant mass (#pi^{+}#pi^{-}) vs DCA;IM(#pi^{+}#pi^{-}) (GeV/c^{2});DCA (cm)", 2000, 0.0, 2.0, 100, 0.0, 10.0 );
  new TH2F( "pipi_nmn_IMvsMomentum_CM", "Invariant mass (#pi^{+}#pi^{-}) vs Momentum;IM(#pi^{+}#pi^{-}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0 );
  new TH2F( "pipi_nmn_IMvsCosT_CM", "Invariant mass (#pi^{+}#pi^{-}) vs cos#theta;IM(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0 );
  new TH2F( "pipi_nmn_IMvsDCA_CM", "Invariant mass (#pi^{+}#pi^{-}) vs DCA;IM(#pi^{+}#pi^{-}) (GeV/c^{2});DCA (cm)", 2000, 0.0, 2.0, 100, 0.0, 10.0 );
  new TH2F( "pipi_nwomn_IMvsMomentum_CM", "Invariant mass (#pi^{+}#pi^{-}) vs Momentum;IM(#pi^{+}#pi^{-}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0 );
  new TH2F( "pipi_nwomn_IMvsCosT_CM", "Invariant mass (#pi^{+}#pi^{-}) vs cos#theta;IM(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0 );
  new TH2F( "pipi_nwomn_IMvsDCA_CM", "Invariant mass (#pi^{+}#pi^{-}) vs DCA;IM(#pi^{+}#pi^{-}) (GeV/c^{2});DCA (cm)", 2000, 0.0, 2.0, 100, 0.0, 10.0 );
  new TH2F( "pipi_nk0mn_IMvsMomentum_CM", "Invariant mass (#pi^{+}#pi^{-}) vs Momentum;IM(#pi^{+}#pi^{-}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0 );
  new TH2F( "pipi_nk0mn_IMvsCosT_CM", "Invariant mass (#pi^{+}#pi^{-}) vs cos#theta;IM(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0 );
  new TH2F( "pipi_nk0mn_IMvsDCA_CM", "Invariant mass (#pi^{+}#pi^{-}) vs DCA;IM(#pi^{+}#pi^{-}) (GeV/c^{2});DCA (cm)", 2000, 0.0, 2.0, 100, 0.0, 10.0 );
  new TH2F( "pipi_sp_IMvsMomentum_CM", "Invariant mass (#pi^{+}#pi^{-}) vs Momentum;IM(#pi^{+}#pi^{-}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0 );
  new TH2F( "pipi_sp_IMvsCosT_CM", "Invariant mass (#pi^{+}#pi^{-}) vs cos#theta;IM(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0 );
  new TH2F( "pipi_sp_IMvsDCA_CM", "Invariant mass (#pi^{+}#pi^{-}) vs DCA;IM(#pi^{+}#pi^{-}) (GeV/c^{2});DCA (cm)", 2000, 0.0, 2.0, 100, 0.0, 10.0 );
  new TH2F( "pipi_wosp_IMvsMomentum_CM", "Invariant mass (#pi^{+}#pi^{-}) vs Momentum;IM(#pi^{+}#pi^{-}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0 );
  new TH2F( "pipi_wosp_IMvsCosT_CM", "Invariant mass (#pi^{+}#pi^{-}) vs cos#theta;IM(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0 );
  new TH2F( "pipi_wosp_IMvsDCA_CM", "Invariant mass (#pi^{+}#pi^{-}) vs DCA;IM(#pi^{+}#pi^{-}) (GeV/c^{2});DCA (cm)", 2000, 0.0, 2.0, 100, 0.0, 10.0 );
  new TH2F( "pipi_spwok0_IMvsMomentum_CM", "Invariant mass (#pi^{+}#pi^{-}) vs Momentum;IM(#pi^{+}#pi^{-}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0 );
  new TH2F( "pipi_spwok0_IMvsCosT_CM", "Invariant mass (#pi^{+}#pi^{-}) vs cos#theta;IM(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0 );
  new TH2F( "pipi_spwok0_IMvsDCA_CM", "Invariant mass (#pi^{+}#pi^{-}) vs DCA;IM(#pi^{+}#pi^{-}) (GeV/c^{2});DCA (cm)", 2000, 0.0, 2.0, 100, 0.0, 10.0 );
  new TH2F( "pipi_spwomn_IMvsMomentum_CM", "Invariant mass (#pi^{+}#pi^{-}) vs Momentum;IM(#pi^{+}#pi^{-}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0 );
  new TH2F( "pipi_spwomn_IMvsCosT_CM", "Invariant mass (#pi^{+}#pi^{-}) vs cos#theta;IM(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0 );
  new TH2F( "pipi_spwomn_IMvsDCA_CM", "Invariant mass (#pi^{+}#pi^{-}) vs DCA;IM(#pi^{+}#pi^{-}) (GeV/c^{2});DCA (cm)", 2000, 0.0, 2.0, 100, 0.0, 10.0 );
  new TH2F( "pipi_sm_IMvsMomentum_CM", "Invariant mass (#pi^{+}#pi^{-}) vs Momentum;IM(#pi^{+}#pi^{-}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0 );
  new TH2F( "pipi_sm_IMvsCosT_CM", "Invariant mass (#pi^{+}#pi^{-}) vs cos#theta;IM(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0 );
  new TH2F( "pipi_sm_IMvsDCA_CM", "Invariant mass (#pi^{+}#pi^{-}) vs DCA;IM(#pi^{+}#pi^{-}) (GeV/c^{2});DCA (cm)", 2000, 0.0, 2.0, 100, 0.0, 10.0 );
  new TH2F( "pipi_wosm_IMvsMomentum_CM", "Invariant mass (#pi^{+}#pi^{-}) vs Momentum;IM(#pi^{+}#pi^{-}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0 );
  new TH2F( "pipi_wosm_IMvsCosT_CM", "Invariant mass (#pi^{+}#pi^{-}) vs cos#theta;IM(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0 );
  new TH2F( "pipi_wosm_IMvsDCA_CM", "Invariant mass (#pi^{+}#pi^{-}) vs DCA;IM(#pi^{+}#pi^{-}) (GeV/c^{2});DCA (cm)", 2000, 0.0, 2.0, 100, 0.0, 10.0 );
  new TH2F( "pipi_smwok0_IMvsMomentum_CM", "Invariant mass (#pi^{+}#pi^{-}) vs Momentum;IM(#pi^{+}#pi^{-}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0 );
  new TH2F( "pipi_smwok0_IMvsCosT_CM", "Invariant mass (#pi^{+}#pi^{-}) vs cos#theta;IM(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0 );
  new TH2F( "pipi_smwok0_IMvsDCA_CM", "Invariant mass (#pi^{+}#pi^{-}) vs DCA;IM(#pi^{+}#pi^{-}) (GeV/c^{2});DCA (cm)", 2000, 0.0, 2.0, 100, 0.0, 10.0 );
  new TH2F( "pipi_smwomn_IMvsMomentum_CM", "Invariant mass (#pi^{+}#pi^{-}) vs Momentum;IM(#pi^{+}#pi^{-}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0 );
  new TH2F( "pipi_smwomn_IMvsCosT_CM", "Invariant mass (#pi^{+}#pi^{-}) vs cos#theta;IM(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0 );
  new TH2F( "pipi_smwomn_IMvsDCA_CM", "Invariant mass (#pi^{+}#pi^{-}) vs DCA;IM(#pi^{+}#pi^{-}) (GeV/c^{2});DCA (cm)", 2000, 0.0, 2.0, 100, 0.0, 10.0 );
  new TH2F( "pipi_nk0mnwos_IMvsMomentum_CM", "Invariant mass (#pi^{+}#pi^{-}) vs Momentum;IM(#pi^{+}#pi^{-}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0 );
  new TH2F( "pipi_nk0mnwos_IMvsCosT_CM", "Invariant mass (#pi^{+}#pi^{-}) vs cos#theta;IM(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0 );
  new TH2F( "pipi_nk0mnwos_IMvsDCA_CM", "Invariant mass (#pi^{+}#pi^{-}) vs DCA;IM(#pi^{+}#pi^{-}) (GeV/c^{2});DCA (cm)", 2000, 0.0, 2.0, 100, 0.0, 10.0 );

  new TH2F( "npip_IMvsMomentum_CM", "Invariant mass (n#pi^{+}) vs Momentum;IM(n#pi^{+}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0 );
  new TH2F( "npip_IMvsCosT_CM", "Invariant mass (n#pi^{+}) vs cos#theta;IM(n#pi^{+}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0 );
  new TH2F( "npim_IMvsMomentum_CM", "Invariant mass (n#pi^{-}) vs Momentum;IM(n#pi^{-}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0 );
  new TH2F( "npim_IMvsCosT_CM", "Invariant mass (n#pi^{-}) vs cos#theta;IM(n#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0 );
  new TH2F( "npip_k0_IMvsMomentum_CM", "Invariant mass (n#pi^{+}) vs Momentum;IM(n#pi^{+}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0 );
  new TH2F( "npip_k0_IMvsCosT_CM", "Invariant mass (n#pi^{+}) vs cos#theta;IM(n#pi^{+}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0 );
  new TH2F( "npip_wok0_IMvsMomentum_CM", "Invariant mass (n#pi^{+}) vs Momentum;IM(n#pi^{+}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0 );
  new TH2F( "npip_wok0_IMvsCosT_CM", "Invariant mass (n#pi^{+}) vs cos#theta;IM(n#pi^{+}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0 );
  new TH2F( "npim_k0_IMvsMomentum_CM", "Invariant mass (n#pi^{-}) vs Momentum;IM(n#pi^{-}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0 );
  new TH2F( "npim_k0_IMvsCosT_CM", "Invariant mass (n#pi^{-}) vs cos#theta;IM(n#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0 );
  new TH2F( "npim_wok0_IMvsMomentum_CM", "Invariant mass (n#pi^{-}) vs Momentum;IM(n#pi^{-}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0 );
  new TH2F( "npim_wok0_IMvsCosT_CM", "Invariant mass (n#pi^{-}) vs cos#theta;IM(n#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0 );
  new TH2F( "npipi_IMvsMomentum_CM", "Invariant mass (n#pi^{+}#pi^{-}) vs Momentum;IM(n#pi^{+}#pi^{-}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0 );
  new TH2F( "npipi_IMvsCosT_CM", "Invariant mass (n#pi^{+}#pi^{-}) vs cos#theta;IM(n#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0 );
  new TH2F( "npipi_k0_IMvsMomentum_CM", "Invariant mass (n#pi^{+}#pi^{-}) vs Momentum;IM(n#pi^{+}#pi^{-}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0 );
  new TH2F( "npipi_k0_IMvsCosT_CM", "Invariant mass (n#pi^{+}#pi^{-}) vs cos#theta;IM(n#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0 );
  new TH2F( "npipi_wok0_IMvsMomentum_CM", "Invariant mass (n#pi^{+}#pi^{-}) vs Momentum;IM(n#pi^{+}#pi^{-}) (GeV/c^{2});Momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0 );
  new TH2F( "npipi_wok0_IMvsCosT_CM", "Invariant mass (n#pi^{+}#pi^{-}) vs cos#theta;IM(n#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0 );
  new TH2F( "npipvsnpim_IM_CM", "IM(n#pi^{+}) vs IM(n#pi^{-});IM(n#pi^{+}) (GeV/c^{2});IM(n#pi^{-}) (GeV/c^{2})", 2000, 0.0, 2.0, 2000, 0.0, 2.0 );
  new TH2F( "pipivsnpip_IM_CM", "IM(#pi^{+}#pi^{-}) vs IM(n#pi^{+});IM(#pi^{+}#pi^{-}) (GeV/c^{2});IM(n#pi^{+}) (GeV/c^{2})", 2000, 0.0, 2.0, 2000, 0.0, 2.0 );
  new TH2F( "pipivsnpim_IM_CM", "IM(#pi^{+}#pi^{-}) vs IM(n#pi^{-});IM(#pi^{+}#pi^{-}) (GeV/c^{2});IM(n#pi^{-}) (GeV/c^{2})", 2000, 0.0, 2.0, 2000, 0.0, 2.0 );
  // Missing mass of production
  std::cout << "Define Histograms for p(K-,n)X M.M." << std::endl;
  new TH2F("pKpip_MMvsMomentum_CM", "p(K^{-},#pi^{+})X missing mass vs. missing momentum;MM_{p}(#pi^{+}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKpip_MMvsCosT_CM", "p(K^{-},#pi^{+})X missing mass vs. cos#theta;MM_{p}(#pi^{+}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKpim_MMvsMomentum_CM", "p(K^{-},#pi^{-})X missing mass vs. missing momentum;MM_{p}(#pi^{-}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKpim_MMvsCosT_CM", "p(K^{-},#pi^{-})X missing mass vs. cos#theta;MM_{p}(#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKpipi_MMvsMomentum_CM", "p(K^{-},#pi^{+}#pi^{-})X missing mass vs. missing momentum;MM_{p}(#pi^{+}#pi^{-}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKpipi_MMvsCosT_CM", "p(K^{-},#pi^{+}#pi^{-})X missing mass vs. cos#theta;MM_{p}(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKpipi_k0_MMvsMomentum_CM", "p(K^{-},#pi^{+}#pi^{-})X missing mass vs. missing momentum;MM_{p}(#pi^{+}#pi^{-}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKpipi_k0_MMvsCosT_CM", "p(K^{-},#pi^{+}#pi^{-})X missing mass vs. cos#theta;MM_{p}(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKpipi_sblk0_MMvsMomentum_CM", "p(K^{-},#pi^{+}#pi^{-})X missing mass vs. missing momentum;MM_{p}(#pi^{+}#pi^{-}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKpipi_sblk0_MMvsCosT_CM", "p(K^{-},#pi^{+}#pi^{-})X missing mass vs. cos#theta;MM_{p}(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKpipi_sbuk0_MMvsMomentum_CM", "p(K^{-},#pi^{+}#pi^{-})X missing mass vs. missing momentum;MM_{p}(#pi^{+}#pi^{-}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKpipi_sbuk0_MMvsCosT_CM", "p(K^{-},#pi^{+}#pi^{-})X missing mass vs. cos#theta;MM_{p}(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKpipi_sbk0_MMvsMomentum_CM", "p(K^{-},#pi^{+}#pi^{-})X missing mass vs. missing momentum;MM_{p}(#pi^{+}#pi^{-}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKpipi_sbk0_MMvsCosT_CM", "p(K^{-},#pi^{+}#pi^{-})X missing mass vs. cos#theta;MM_{p}(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKpipi_wok0_MMvsMomentum_CM", "p(K^{-},#pi^{+}#pi^{-})X missing mass vs. missing momentum;MM_{p}(#pi^{+}#pi^{-}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKpipi_wok0_MMvsCosT_CM", "p(K^{-},#pi^{+}#pi^{-})X missing mass vs. cos#theta;MM_{p}(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKpipi_mn_MMvsMomentum_CM", "p(K^{-},#pi^{+}#pi^{-})X missing mass vs. missing momentum;MM_{p}(#pi^{+}#pi^{-}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKpipi_mn_MMvsCosT_CM", "p(K^{-},#pi^{+}#pi^{-})X missing mass vs. cos#theta;MM_{p}(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKpipi_womn_MMvsMomentum_CM", "p(K^{-},#pi^{+}#pi^{-})X missing mass vs. missing momentum;MM_{p}(#pi^{+}#pi^{-}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKpipi_womn_MMvsCosT_CM", "p(K^{-},#pi^{+}#pi^{-})X missing mass vs. cos#theta;MM_{p}(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKpipi_k0mn_MMvsMomentum_CM", "p(K^{-},#pi^{+}#pi^{-})X missing mass vs. missing momentum;MM_{p}(#pi^{+}#pi^{-}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKpipi_k0mn_MMvsCosT_CM", "p(K^{-},#pi^{+}#pi^{-})X missing mass vs. cos#theta;MM_{p}(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKpipi_sblk0mn_MMvsMomentum_CM", "p(K^{-},#pi^{+}#pi^{-})X missing mass vs. missing momentum;MM_{p}(#pi^{+}#pi^{-}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKpipi_sblk0mn_MMvsCosT_CM", "p(K^{-},#pi^{+}#pi^{-})X missing mass vs. cos#theta;MM_{p}(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKpipi_sbuk0mn_MMvsMomentum_CM", "p(K^{-},#pi^{+}#pi^{-})X missing mass vs. missing momentum;MM_{p}(#pi^{+}#pi^{-}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKpipi_sbuk0mn_MMvsCosT_CM", "p(K^{-},#pi^{+}#pi^{-})X missing mass vs. cos#theta;MM_{p}(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKpipi_sbk0mn_MMvsMomentum_CM", "p(K^{-},#pi^{+}#pi^{-})X missing mass vs. missing momentum;MM_{p}(#pi^{+}#pi^{-}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKpipi_sbk0mn_MMvsCosT_CM", "p(K^{-},#pi^{+}#pi^{-})X missing mass vs. cos#theta;MM_{p}(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKpipi_n_MMvsMomentum_CM", "p(K^{-},#pi^{+}#pi^{-})X missing mass vs. missing momentum;MM_{p}(#pi^{+}#pi^{-}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKpipi_n_MMvsCosT_CM", "p(K^{-},#pi^{+}#pi^{-})X missing mass vs. cos#theta;MM_{p}(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKpipi_nk0_MMvsMomentum_CM", "p(K^{-},#pi^{+}#pi^{-})X missing mass vs. missing momentum;MM_{p}(#pi^{+}#pi^{-}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKpipi_nk0_MMvsCosT_CM", "p(K^{-},#pi^{+}#pi^{-})X missing mass vs. cos#theta;MM_{p}(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKpipi_nsblk0_MMvsMomentum_CM", "p(K^{-},#pi^{+}#pi^{-})X missing mass vs. missing momentum;MM_{p}(#pi^{+}#pi^{-}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKpipi_nsblk0_MMvsCosT_CM", "p(K^{-},#pi^{+}#pi^{-})X missing mass vs. cos#theta;MM_{p}(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKpipi_nsbuk0_MMvsMomentum_CM", "p(K^{-},#pi^{+}#pi^{-})X missing mass vs. missing momentum;MM_{p}(#pi^{+}#pi^{-}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKpipi_nsbuk0_MMvsCosT_CM", "p(K^{-},#pi^{+}#pi^{-})X missing mass vs. cos#theta;MM_{p}(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKpipi_nsbk0_MMvsMomentum_CM", "p(K^{-},#pi^{+}#pi^{-})X missing mass vs. missing momentum;MM_{p}(#pi^{+}#pi^{-}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKpipi_nsbk0_MMvsCosT_CM", "p(K^{-},#pi^{+}#pi^{-})X missing mass vs. cos#theta;MM_{p}(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKpipi_nwok0_MMvsMomentum_CM", "p(K^{-},#pi^{+}#pi^{-})X missing mass vs. missing momentum;MM_{p}(#pi^{+}#pi^{-}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKpipi_nwok0_MMvsCosT_CM", "p(K^{-},#pi^{+}#pi^{-})X missing mass vs. cos#theta;MM_{p}(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKpipi_nmn_MMvsMomentum_CM", "p(K^{-},#pi^{+}#pi^{-})X missing mass vs. missing momentum;MM_{p}(#pi^{+}#pi^{-}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKpipi_nmn_MMvsCosT_CM", "p(K^{-},#pi^{+}#pi^{-})X missing mass vs. cos#theta;MM_{p}(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKpipi_nwomn_MMvsMomentum_CM", "p(K^{-},#pi^{+}#pi^{-})X missing mass vs. missing momentum;MM_{p}(#pi^{+}#pi^{-}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKpipi_nwomn_MMvsCosT_CM", "p(K^{-},#pi^{+}#pi^{-})X missing mass vs. cos#theta;MM_{p}(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKpipi_nk0mn_MMvsMomentum_CM", "p(K^{-},#pi^{+}#pi^{-})X missing mass vs. missing momentum;MM_{p}(#pi^{+}#pi^{-}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKpipi_nk0mn_MMvsCosT_CM", "p(K^{-},#pi^{+}#pi^{-})X missing mass vs. cos#theta;MM_{p}(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKpipi_sp_MMvsMomentum_CM", "p(K^{-},#pi^{+}#pi^{-})X missing mass vs. missing momentum;MM_{p}(#pi^{+}#pi^{-}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKpipi_sp_MMvsCosT_CM", "p(K^{-},#pi^{+}#pi^{-})X missing mass vs. cos#theta;MM_{p}(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKpipi_wosp_MMvsMomentum_CM", "p(K^{-},#pi^{+}#pi^{-})X missing mass vs. missing momentum;MM_{p}(#pi^{+}#pi^{-}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKpipi_wosp_MMvsCosT_CM", "p(K^{-},#pi^{+}#pi^{-})X missing mass vs. cos#theta;MM_{p}(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKpipi_spwok0_MMvsMomentum_CM", "p(K^{-},#pi^{+}#pi^{-})X missing mass vs. missing momentum;MM_{p}(#pi^{+}#pi^{-}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKpipi_spwok0_MMvsCosT_CM", "p(K^{-},#pi^{+}#pi^{-})X missing mass vs. cos#theta;MM_{p}(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKpipi_spwomn_MMvsMomentum_CM", "p(K^{-},#pi^{+}#pi^{-})X missing mass vs. missing momentum;MM_{p}(#pi^{+}#pi^{-}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKpipi_spwomn_MMvsCosT_CM", "p(K^{-},#pi^{+}#pi^{-})X missing mass vs. cos#theta;MM_{p}(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKpipi_sm_MMvsMomentum_CM", "p(K^{-},#pi^{+}#pi^{-})X missing mass vs. missing momentum;MM_{p}(#pi^{+}#pi^{-}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKpipi_sm_MMvsCosT_CM", "p(K^{-},#pi^{+}#pi^{-})X missing mass vs. cos#theta;MM_{p}(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKpipi_wosm_MMvsMomentum_CM", "p(K^{-},#pi^{+}#pi^{-})X missing mass vs. missing momentum;MM_{p}(#pi^{+}#pi^{-}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKpipi_wosm_MMvsCosT_CM", "p(K^{-},#pi^{+}#pi^{-})X missing mass vs. cos#theta;MM_{p}(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKpipi_smwok0_MMvsMomentum_CM", "p(K^{-},#pi^{+}#pi^{-})X missing mass vs. missing momentum;MM_{p}(#pi^{+}#pi^{-}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKpipi_smwok0_MMvsCosT_CM", "p(K^{-},#pi^{+}#pi^{-})X missing mass vs. cos#theta;MM_{p}(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKpipi_smwomn_MMvsMomentum_CM", "p(K^{-},#pi^{+}#pi^{-})X missing mass vs. missing momentum;MM_{p}(#pi^{+}#pi^{-}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKpipi_smwomn_MMvsCosT_CM", "p(K^{-},#pi^{+}#pi^{-})X missing mass vs. cos#theta;MM_{p}(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKpipi_nk0mnwos_MMvsMomentum_CM", "p(K^{-},#pi^{+}#pi^{-})X missing mass vs. missing momentum;MM_{p}(#pi^{+}#pi^{-}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKpipi_nk0mnwos_MMvsCosT_CM", "p(K^{-},#pi^{+}#pi^{-})X missing mass vs. cos#theta;MM_{p}(#pi^{+}#pi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);

  new TH2F("pKn_MMvsMomentum_CM", "p(K^{-},n)X missing mass vs. missing momentum;MM_{p}(n) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKn_MMvsCosT_CM", "p(K^{-},n)X missing mass vs. cos#theta;MM_{p}(n) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKn_k0_MMvsMomentum_CM", "p(K^{-},n)X missing mass vs. missing momentum;MM_{p}(n) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKn_k0_MMvsCosT_CM", "p(K^{-},n)X missing mass vs. cos#theta;MM_{p}(n) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKn_sblk0_MMvsMomentum_CM", "p(K^{-},n)X missing mass vs. missing momentum;MM_{p}(n) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKn_sblk0_MMvsCosT_CM", "p(K^{-},n)X missing mass vs. cos#theta;MM_{p}(n) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKn_sbuk0_MMvsMomentum_CM", "p(K^{-},n)X missing mass vs. missing momentum;MM_{p}(n) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKn_sbuk0_MMvsCosT_CM", "p(K^{-},n)X missing mass vs. cos#theta;MM_{p}(n) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKn_sbk0_MMvsMomentum_CM", "p(K^{-},n)X missing mass vs. missing momentum;MM_{p}(n) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKn_sbk0_MMvsCosT_CM", "p(K^{-},n)X missing mass vs. cos#theta;MM_{p}(n) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKn_wok0_MMvsMomentum_CM", "p(K^{-},n)X missing mass vs. missing momentum;MM_{p}(n) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKn_wok0_MMvsCosT_CM", "p(K^{-},n)X missing mass vs. cos#theta;MM_{p}(n) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKn_mn_MMvsMomentum_CM", "p(K^{-},n)X missing mass vs. missing momentum;MM_{p}(n) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKn_mn_MMvsCosT_CM", "p(K^{-},n)X missing mass vs. cos#theta;MM_{p}(n) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKn_womn_MMvsMomentum_CM", "p(K^{-},n)X missing mass vs. missing momentum;MM_{p}(n) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKn_womn_MMvsCosT_CM", "p(K^{-},n)X missing mass vs. cos#theta;MM_{p}(n) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKn_k0mn_MMvsMomentum_CM", "p(K^{-},n)X missing mass vs. missing momentum;MM_{p}(n) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKn_k0mn_MMvsCosT_CM", "p(K^{-},n)X missing mass vs. cos#theta;MM_{p}(n) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKn_sblk0mn_MMvsMomentum_CM", "p(K^{-},n)X missing mass vs. missing momentum;MM_{p}(n) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKn_sblk0mn_MMvsCosT_CM", "p(K^{-},n)X missing mass vs. cos#theta;MM_{p}(n) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKn_sbuk0mn_MMvsMomentum_CM", "p(K^{-},n)X missing mass vs. missing momentum;MM_{p}(n) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKn_sbuk0mn_MMvsCosT_CM", "p(K^{-},n)X missing mass vs. cos#theta;MM_{p}(n) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKn_sbk0mn_MMvsMomentum_CM", "p(K^{-},n)X missing mass vs. missing momentum;MM_{p}(n) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKn_sbk0mn_MMvsCosT_CM", "p(K^{-},n)X missing mass vs. cos#theta;MM_{p}(n) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKn_sp_MMvsMomentum_CM", "p(K^{-},n)X missing mass vs. missing momentum;MM_{p}(n) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKn_sp_MMvsCosT_CM", "p(K^{-},n)X missing mass vs. cos#theta;MM_{p}(n) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKn_wosp_MMvsMomentum_CM", "p(K^{-},n)X missing mass vs. missing momentum;MM_{p}(n) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKn_wosp_MMvsCosT_CM", "p(K^{-},n)X missing mass vs. cos#theta;MM_{p}(n) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKn_spwok0_MMvsMomentum_CM", "p(K^{-},n)X missing mass vs. missing momentum;MM_{p}(n) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKn_spwok0_MMvsCosT_CM", "p(K^{-},n)X missing mass vs. cos#theta;MM_{p}(n) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKn_spwomn_MMvsMomentum_CM", "p(K^{-},n)X missing mass vs. missing momentum;MM_{p}(n) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKn_spwomn_MMvsCosT_CM", "p(K^{-},n)X missing mass vs. cos#theta;MM_{p}(n) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKn_sm_MMvsMomentum_CM", "p(K^{-},n)X missing mass vs. missing momentum;MM_{p}(n) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKn_sm_MMvsCosT_CM", "p(K^{-},n)X missing mass vs. cos#theta;MM_{p}(n) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKn_wosm_MMvsMomentum_CM", "p(K^{-},n)X missing mass vs. missing momentum;MM_{p}(n) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKn_wosm_MMvsCosT_CM", "p(K^{-},n)X missing mass vs. cos#theta;MM_{p}(n) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKn_smwok0_MMvsMomentum_CM", "p(K^{-},n)X missing mass vs. missing momentum;MM_{p}(n) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKn_smwok0_MMvsCosT_CM", "p(K^{-},n)X missing mass vs. cos#theta;MM_{p}(n) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKn_smwomn_MMvsMomentum_CM", "p(K^{-},n)X missing mass vs. missing momentum;MM_{p}(n) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKn_smwomn_MMvsCosT_CM", "p(K^{-},n)X missing mass vs. cos#theta;MM_{p}(n) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKn_k0mnwos_MMvsMomentum_CM", "p(K^{-},n)X missing mass vs. missing momentum;MM_{p}(n) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKn_k0mnwos_MMvsCosT_CM", "p(K^{-},n)X missing mass vs. cos#theta;MM_{p}(n) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKn_sblk0mnwos_MMvsMomentum_CM", "p(K^{-},n)X missing mass vs. missing momentum;MM_{p}(n) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKn_sblk0mnwos_MMvsCosT_CM", "p(K^{-},n)X missing mass vs. cos#theta;MM_{p}(n) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKn_sbuk0mnwos_MMvsMomentum_CM", "p(K^{-},n)X missing mass vs. missing momentum;MM_{p}(n) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKn_sbuk0mnwos_MMvsCosT_CM", "p(K^{-},n)X missing mass vs. cos#theta;MM_{p}(n) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKn_sbk0mnwos_MMvsMomentum_CM", "p(K^{-},n)X missing mass vs. missing momentum;MM_{p}(n) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKn_sbk0mnwos_MMvsCosT_CM", "p(K^{-},n)X missing mass vs. cos#theta;MM_{p}(n) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);

  new TH2F("pKnpip_MMvsMomentum_CM", "p(K^{-},n#pi^{+})X missing mass vs. missing momentum;MM_{p}(npi^{+}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKnpip_k0_MMvsMomentum_CM", "p(K^{-},n#pi^{+})X missing mass vs. missing momentum;MM_{p}(npi^{+}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKnpip_k0_MMvsCosT_CM", "p(K^{-},n#pi^{+})X missing mass vs. cos#theta;MM_{p}(npi^{+}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKnpip_wok0_MMvsMomentum_CM", "p(K^{-},n#pi^{+})X missing mass vs. missing momentum;MM_{p}(npi^{+}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKnpip_wok0_MMvsCosT_CM", "p(K^{-},n#pi^{+})X missing mass vs. cos#theta;MM_{p}(npi^{+}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKnpip_MMvsCosT_CM", "p(K^{-},n#pi^{+})X missing mass vs. cos#theta;MM_{p}(npi^{+}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKnpim_MMvsMomentum_CM", "p(K^{-},n#pi^{-})X missing mass vs. missing momentum;MM_{p}(npi^{-}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKnpim_MMvsCosT_CM", "p(K^{-},n#pi^{-})X missing mass vs. cos#theta;MM_{p}(npi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKnpim_k0_MMvsMomentum_CM", "p(K^{-},n#pi^{-})X missing mass vs. missing momentum;MM_{p}(npi^{-}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKnpim_k0_MMvsCosT_CM", "p(K^{-},n#pi^{-})X missing mass vs. cos#theta;MM_{p}(npi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);
  new TH2F("pKnpim_wok0_MMvsMomentum_CM", "p(K^{-},n#pi^{-})X missing mass vs. missing momentum;MM_{p}(npi^{-}) (GeV/c^{2});Missing momentum (GeV/c)", 2000, 0.0, 2.0, 1000, 0.0, 1.0);
  new TH2F("pKnpim_wok0_MMvsCosT_CM", "p(K^{-},n#pi^{-})X missing mass vs. cos#theta;MM_{p}(npi^{-}) (GeV/c^{2});cos#theta", 2000, 0.0, 2.0, 2000, -1.0, 1.0);

  return true;
}
