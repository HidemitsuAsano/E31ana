#include "MyHistKN_pipi.h"

using namespace std;

void initHistKN_pipi()
{
  new TH1F("KNpipi_MM",         "d(K^{-}, n #pi^{+} #pi^{-})\"X\"", 2000, 0.0, 2.0);

  new TH1F("KN_MM_pipi",        "d(K^{-}, n)\"X\"",                 1000, 1.0, 2.0);
  new TH1F("KN_MM_pipi_wN",     "d(K^{-}, n)\"X\"",                 1000, 1.0, 2.0);
  new TH1F("KN_MM_pipi_wK0_wN",       "d(K^{-}, n)\"X\"",           1000, 1.0, 2.0);
  new TH1F("KN_MM_pipi_wSm_wN",       "d(K^{-}, n)\"X\"",           1000, 1.0, 2.0);
  new TH1F("KN_MM_pipi_wSp_wN",       "d(K^{-}, n)\"X\"",           1000, 1.0, 2.0);
  new TH1F("KN_MM_pipi_woAll_wN",     "d(K^{-}, n)\"X\"",           1000, 1.0, 2.0);

  new TH1F("CDS_IM_pipi_wN_wN", "CDS IM #pi^{+} #pi^{-}",           1000, 0.0, 1.0);
  new TH1F("Npim_IM_pip_wN",    "n #pi^{-} IM",                     1000, 1.0, 2.0);
  new TH1F("Npip_IM_pim_wN",    "n #pi^{+} IM",                     1000, 1.0, 2.0);
}

void fillHistKN_pipi(EventHeader *header, BeamLineHitMan *blMan, CDSHitMan *cdsMan, CDSTrackingMan *cdstrackMan, AnaInfo *anaInfo)
{
  if( !header->IsTrig(Trig_Neutral) ) return;
  if( !header->IsTrig(Trig_Kaon) ) return;

  if( anaInfo->nFNeutral()!=1 || anaInfo->forwardNeutral(0)->pid()!=F_Neutron ) return;
  vector<HodoscopeLikeHit*> CVChits=MyTools::getHodo(blMan, CID_CVC);
  vector<HodoscopeLikeHit*> BVChits=MyTools::getHodo(blMan, CID_BVC);
  if( CVChits.size()!=0 || BVChits.size()!=0 ) return;
  if( anaInfo->nBeam()!=1 || !anaInfo->beam(0)->flag() ) return;
  if( !anaInfo->minDCA() || GeomTools::GetID(anaInfo->minDCA()->vertexBeam())!=CID_Fiducial ) return;
  if( cdstrackMan->nGoodTrack()!=2 ) return;
  if( anaInfo->nCDS(CDS_PiMinus)!=1 || anaInfo->nCDS(CDS_PiPlus)!=1 ) return;

  BeamInfo *beam=anaInfo->beam(0);
  ForwardNeutralInfo *fnInfo=anaInfo->forwardNeutral(0);
  CDSInfo *pimInfo=anaInfo->CDS(CDS_PiMinus, 0);
  CDSInfo *pipInfo=anaInfo->CDS(CDS_PiPlus, 0);
  CDS2Info *pipiInfo=anaInfo->CDS2(CDS_PiMinus, CDS_PiPlus, 0);
  if( !pimInfo->flag() || !pipInfo->flag() || !pipiInfo->flag() ) return;
  if( pimInfo->CDHseg()==pipInfo->CDHseg() ) return;

  TLorentzVector beam_lmom=beam->lmom();
  //  cout<<" KN pipi ev : "<<header->ev()<<endl;
  //  cout<<Form(" Beam lmom m=%4.3lf (%4.3lf, %4.3lf, %4.3lf)", beam_lmom.M(), beam_lmom.X(), beam_lmom.Y(), beam_lmom.Z())<<endl;
  TLorentzVector fn_lmom=fnInfo->lmom();
  TLorentzVector mm_lmom=beam_lmom+D_LMOM-fn_lmom;
  TLorentzVector pim_lmom=pimInfo->lmom();
  TLorentzVector pip_lmom=pipInfo->lmom();

  TLorentzVector pipi_lmom=pim_lmom+pip_lmom;
  TLorentzVector npim_lmom=fn_lmom+pim_lmom;
  TLorentzVector npip_lmom=fn_lmom+pip_lmom;

  TLorentzVector knpipi_lmom=mm_lmom-pimInfo->lmom()-pipInfo->lmom();
  bool K0_flag=false, Sm_flag=false, Sp_flag=false;
  if( Npipi_K0_MIN<pipi_lmom.M() && pipi_lmom.M()<Npipi_K0_MAX ) K0_flag=true;
  if( Npipi_Sm_MIN<npim_lmom.M() && npim_lmom.M()<Npipi_Sm_MAX ) Sm_flag=true;
  if( Npipi_Sp_MIN<npip_lmom.M() && npip_lmom.M()<Npipi_Sp_MAX ) Sp_flag=true;

  MyHistTools::fillTH("KNpipi_MM", knpipi_lmom.M());
  MyHistTools::fillTH("KN_MM_pipi", mm_lmom.M());
  if( 0.9<knpipi_lmom.M() && knpipi_lmom.M()<0.98 ){
    MyHistTools::fillTH("KN_MM_pipi_wN", mm_lmom.M());

    MyHistTools::fillTH("CDS_IM_pipi_wN_wN", pipiInfo->im());
    MyHistTools::fillTH("Npim_IM_pip_wN", (pim_lmom+fn_lmom).M());
    MyHistTools::fillTH("Npip_IM_pim_wN", (pip_lmom+fn_lmom).M());

    if( K0_flag ){
      MyHistTools::fillTH("KN_MM_pipi_wK0_wN", mm_lmom.M());
    }
    if( Sm_flag ){
      MyHistTools::fillTH("KN_MM_pipi_wSm_wN", mm_lmom.M());
    }
    if( Sp_flag ){
      MyHistTools::fillTH("KN_MM_pipi_wSp_wN", mm_lmom.M());
    }
    if( !K0_flag && !Sm_flag && !Sp_flag ){
      MyHistTools::fillTH("KN_MM_pipi_woAll_wN", mm_lmom.M());
    }
  }
}
