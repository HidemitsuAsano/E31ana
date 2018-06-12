#include "MyHistReadKNpip.h"

using namespace std;

void initHistReadKNpip()
{
  new TH2F("KNpip_MM_Npip_IM", "d(K^{-}, n #pi^{+})\"X\" vs n #pi^{+} IM", 2000, 1.0, 2.0, 1000, 1.0, 2.0);
  new TH1F("Npip_IM",      "n #pi^{+} IM", 1000, 1.0, 2.0);
  new TH1F("Npip_IM_true", "n #pi^{+} IM", 1000, 1.0, 2.0);
}

void fillHistReadKNpip(EventHeader *header, BeamLineHitMan *blMan, AnaInfo *anaInfo)
{
  if( anaInfo->nFNeutral()!=1 ) return;
  if( header ){
    if( !header->IsTrig(Trig_Neutral) ) return;
  }
  vector<HodoscopeLikeHit*> CVChits=MyTools::getHodo(blMan, CID_CVC);
  vector<HodoscopeLikeHit*> BVChits=MyTools::getHodo(blMan, CID_BVC);
  if( CVChits.size()!=0 || BVChits.size()!=0 ) return;
  if( anaInfo->nBeam()!=1 || !anaInfo->beam(0)->flag() ) return;

  ForwardNeutralInfo *fnInfo=anaInfo->forwardNeutral(0);
  if( fnInfo->pid()!=F_Neutron ) return;

  TLorentzVector beam_lmom=anaInfo->beam(0)->lmom();
  TLorentzVector target_lmom=MyAnaTools::target_lmom();
  TLorentzVector fn_lmom=fnInfo->lmom();
  TLorentzVector kn_lmom=beam_lmom+target_lmom-fn_lmom;

  for( int i=0; i<anaInfo->nCDS(CDS_PiPlus); i++ ){
    CDSInfo *pip=anaInfo->CDS(CDS_PiPlus, i);
    TLorentzVector pip_lmom=pip->lmom();
    double im=(fn_lmom+pip_lmom).M();
    double mm=(kn_lmom-pip_lmom).M();

    if(GeomTools::GetID(pip->vertexBeam())==CID_Fiducial ){
      MyHistTools::fillTH("KNpip_MM_Npip_IM", mm, im);
      MyHistTools::fillTH("Npip_IM", im);
      if( Sp_MIN<im &&  im<Sp_MAX ){
	MyHistTools::fillTH("Npip_IM_true", im);
      }
    }
  }
}
