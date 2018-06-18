#include "MyHistReadKPpimpim.h"

using namespace std;

void initHistReadKPpimpim()
{
  new TH2F("DCA_fp_2pim",  "DCA beam #pi^{-}", 1000, 0, 25, 1000, 0, 25);
  new TH2F("DCA2_fp_2pim", "DCA bw #pi^{-}",   1000, 0, 25, 1000, 0, 25);

  new TH1F("KP_MM_pimL", "d(K^{-}, p)\"#pi^{-} #Lambda\"", 2000, 1.0, 2.0);
  new TH1F("KP_MM_pimS0", "d(K^{-}, p)\"#pi^{-} #Sigma^{0}\"", 2000, 1.0, 2.0);

  new TH1F("KPpimpim_MM", "d(K^{-}, p #pi^{-} #pi^{-})\"X\"", 2000, 0.0, 2.0);
  new TH1F("KPpimpim_MM_angle", "d(K^{-}, p #pi^{-} #pi^{-})\"X\"", 2000, 0.0, 2.0);
  new TH1F("KPpimpim_MM_RK", "d(K^{-}, p #pi^{-} #pi^{-})\"X\"", 2000, 0.0, 2.0);

  new TH1F("KPpimpim_MM_wL", "d(K^{-}, p #pi^{-} #pi^{-})\"X\"", 2000, 0.0, 2.0);
  new TH1F("KPpimpim_MM_wS0", "d(K^{-}, p #pi^{-} #pi^{-})\"X\"", 2000, 0.0, 2.0);
  new TH1F("KPpim_MM_Lmass",  "d(K^{-}, p #pi^{-})\"#Labmda\"", 2000, 0.0, 2.0);
  new TH1F("KPpim_MM_S0mass", "d(K^{-}, p #pi^{-})\"#Sigma^{0}\"", 2000, 0.0, 2.0);

  new TH1F("KPpimpim_MM_FDC1x", "FDC1 x vs d(K^{-}, p #pi^{-} #pi^{-})\"X\"", 1000,  -50, 50, 2000, 0.0, 2.0);
  new TH1F("KPpimpim_MM_FDC1y", "FDC1 y vs d(K^{-}, p #pi^{-} #pi^{-})\"X\"", 1000,  -50, 50, 2000, 0.0, 2.0);
  new TH1F("KPpimpim_MM_dx",    "dx/dz vs d(K^{-}, p #pi^{-} #pi^{-})\"X\"", 1000, -0.2, 0.2, 2000, 0.0, 2.0);
  new TH1F("KPpimpim_MM_dy",    "dy/dz x vs d(K^{-}, p #pi^{-} #pi^{-})\"X\"", 1000, -0.2, 0.2, 2000, 0.0, 2.0);
  new TH1F("KPpimpim_MM_FCseg", "CVC/PC seg x vs d(K^{-}, p #pi^{-} #pi^{-})\"X\"",   51, 0.5, 51.5, 2000, 0.0, 2.0);
  new TH1F("KPpimpim_MM_ang",   "angle vs d(K^{-}, p #pi^{-} #pi^{-})\"X\"",   51, 0.5, 51.5, 2000, 0.0, 2.0);
  new TH1F("KPpimpim_MM_dE",    "dE vs d(K^{-}, p #pi^{-} #pi^{-})\"X\"",   51, 0.5, 51.5, 2000, 0.0, 2.0);

  new TH2F("KPpim_KPpim_MM",         "d(K^{-}, p #pi^{-})\"X\" vs d(K^{-}, p #pi^{-})\"X\"", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("KPpim_KPpim_MM_mmP",      "d(K^{-}, p #pi^{-})\"X\" vs d(K^{-}, p #pi^{-})\"X\"", 2000, 0.0, 2.0, 2000, 0.0, 2.0);
  new TH2F("KPpim_KPpim_MM_mmPgamma", "d(K^{-}, p #pi^{-})\"X\" vs d(K^{-}, p #pi^{-})\"X\"", 2000, 0.0, 2.0, 2000, 0.0, 2.0);

  new TNutple("tup_kp_pimL", "", "m");
  new TNutple("tup_kp_pimS0", "", "m");
}

void fillHistReadKPpimpim(EventHeader *header, BeamLineHitMan *blMan, CDSHitMan *cdsMan, CDSTrackingMan *cdstrackMan, AnaInfo *anaInfo)
{
  if( header ){
    if( !header->IsTrig(Trig_Kaon) ) return;
    //    if( !header->trigmode(Mode_KCDH1C) ) return;
    //    if( !header->IsTrig(Trig_Charged) ) return;
  }

  if( anaInfo->nBeam()!=1 || !anaInfo->beam(0)->flag() ) return;
  if( !anaInfo->minDCA() || GeomTools::GetID(anaInfo->minDCA()->vertexBeam())!=CID_Fiducial ) return;

  //  if( anaInfo->nFCharge()!=1 || anaInfo->forwardCharge(0)->pid()!=F_Proton ) return;
  if( anaInfo->nFCharge()!=1 ) return;
  ForwardChargeInfo *fcInfo=anaInfo->forwardCharge(0);
  if( fcInfo->mass2byRK()<FC_P_MIN || FC_P_MAX<fcInfo->mass2byRK() ) return;

  if( anaInfo->nCDS(CDS_PiMinus)!=2 ) return;
  TTree *tree=(TTree*)gFile->Get("ppimpim_event");
  if( tree ){
    tree-> Fill();
  }

  BeamInfo *beam=anaInfo->beam(0);
  CDSInfo *pim0=anaInfo->CDS(CDS_PiMinus, 0);
  CDSInfo *pim1=anaInfo->CDS(CDS_PiMinus, 1);
  if( pim0->dca()>pim1->dca() ) swap(pim0, pim1);

  TLorentzVector beam_lmom=beam->lmom();
  TLorentzVector pim_lmom0=pim0->lmom();
  double dca0=pim0->dca();
  TVector3 lpos, hpos;
  pim1->GetVertex(pim0->vertexBeam(), beam->lmom.Vect()-fc->lmom.Vect()-pim0->lmom.Vect(), lpos, hpos);
  double dca2_0=(lpos-hpos).Mag();s
  TLorentzVector pim_lmom1=pim1->lmom();
  double dca1=pim_lmom1->dca();
  pim0->GetVertex(pim1->vertexBeam(), beam->lmom.Vect()-fc->lmom.Vect()-pim1->lmom.Vect(), lpos, hpos);
  double dca2_1=(lpos-hpos).Mag();
  TLorentzVector fp_lmom=fcInfo->lmom();
  TVector3 fp_mom=fp_lmom.Vect();
  TVector3 fp_mom_ang=fp_mom;
  TVector3 fp_mom_RK=fp_mom;
  fp_mom_ang.SetMag(fcInfo->momByAng());
  fp_mom_RK.SetMag(fcInfo->momByRK());
  TLorentzVector fp_lmom;
  fp_lmom_ang.SetVectM(fp_mom_ang, pMass);
  fp_lmom_RK.SetVectM(fp_mom_RK, pMass);

  TLorentzVector tgt_lmom=MyAnaTools::target_lmom();
  TVector3 FDC1pos=fcInfo->FDC1().pos();
  TVector3 dir=fp_lom.Vect().Unit();

  double mm=(beam_lmom+tgt_lmom-fp_lmom).M();
  double kppimpim_mm= (beam_lmom+tgt_lmom-fp_lmom-pim_lmom0-pim_lmom1).M();
  double kppimpim_mm_ang= (beam_lmom+tgt_lmom-fp_lmom_ang-pim_lmom0-pim_lmom1).M();
  double kppimpim_mm_RK= (beam_lmom+tgt_lmom-fp_lmom_RK-pim_lmom0-pim_lmom1).M();
  double kppim0_mm= (beam_lmom+tgt_lmom-fp_lmom-pim_lmom0).M();
  double kppim1_mm= (beam_lmom+tgt_lmom-fp_lmom-pim_lmom1).M();
  bool mmL_flag=false; bool mmS0_flag=false; bool mmP_flag=false; bool mmPgamma_flag=false;
  if( FC_mmL_MIN<kppim0_mm && kppim0_mm<FC_mmL_MAX ) mmL_flag=true;
  if( FC_mmL_MIN<kppim1_mm && kppim1_mm<FC_mmL_MAX ) mmL_flag=true;
  if( FC_mmS0_MIN<kppim0_mm && kppim0_mm<FC_mmS0_MAX ) mmS0_flag=true;
  if( FC_mmS0_MIN<kppim1_mm && kppim1_mm<FC_mmS0_MAX ) mmS0_flag=true;
  if( FC_mmP_MIN<kppimpim_mm && kppimpim_mm<FC_mmP_MAX ) mmP_flag=true;
  if( 0.99<kppimpim_mm && kppimpim_mm<1.06 ) mmPgamma_flag=true;

  MyHistTools::fillTH("KPpimpim_MM", kppimpim_mm);
  MyHistTools::fillTH("KPpimpim_MM_ang", kppimpim_mm_ang);
  MyHistTools::fillTH("KPpimpim_MM_RK", kppimpim_mm_RK);
  MyHistTools::fillTH("KPpimpim_MM_FDC1x", FDC1pos.X(), kppimpim_mm);
  MyHistTools::fillTH("KPpimpim_MM_FDC1y", FDC1pos.Y(), kppimpim_mm);
  MyHistTools::fillTH("KPpimpim_MM_dx", dir.X(), kppimpim_mm);
  MyHistTools::fillTH("KPpimpim_MM_dy", dirx.Y(), kppimpim_mm);
  MyHistTools::fillTH("DCA_fp_2pim", dca0, dca1);
  MyHistTools::fillTH("DCA2_fp_2pim", dca2_0, dca2_1);
  if( mmL_flag ) MyHistTools::fillTH("KPpimpim_MM_wL", kppimpim_mm);
  if( mmS0_flag ) MyHistTools::fillTH("KPpimpim_MM_wS0", kppimpim_mm);

  MyHistTools::fillTH("KPpim_KPpim_MM", kppim0_mm, kppim1_mm);
  if( mmP_flag ){
    MyHistTools::fillTH("KPpim_KPpim_MM_mmP", kppim0_mm, kppim1_mm);
    double L_mass=kppim1_mm;
    if( fabs(kppim0_mm-lMass)<fabs(kppim1_mm-lMass) ) L_mass=kppim0_mm;
    MyHistTools::fillTH("KPpim_MM_Lmass", L_mass);
  }
  if( mmPgamma_flag ){
    MyHistTools::fillTH("KPpim_KPpim_MM_mmPgamma", kppim0_mm, kppim1_mm);
    double S0_mass=kppim1_mm;
    if( fabs(kppim0_mm-s0Mass)<fabs(kppim1_mm-s0Mass) ) S0_mass=kppim0_mm;
    MyHistTools::fillTH("KPpim_MM_S0mass", S0_mass);
  }

  TNtuple *tup;
  if( mmP_flag && mmL_flag ){
    MyHistTools::fillTH("KP_MM_pimL", mm);
    tup=(TNtuple*)gFile-> Get("tup_kp_pimL"), tup-> fill(mm);
  }
  if( mmPgamma_flag && mmS0_flag ){
    MyHistTools::fillTH("KP_MM_pimS0", mm);
    tup=(TNtuple*)gFile-> Get("tup_kp_pimS0"), tup-> fill(mm);
  }
}
