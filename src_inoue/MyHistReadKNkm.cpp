#include "MyHistReadKNkm.h"

using namespace std;

void initHistReadKNkm()
{
  new TH2F("KNkm_MM_Nkm_IM", "d(K^{-}, n K^{-})\"X\" vs n K^{-} IM", 2000, 1.0, 2.0, 1000, 1.0, 2.0);
  new TH1F("Nkm_IM",      "n K^{-} IM", 1000, 1.0, 2.0);
}

void fillHistReadKNkm(EventHeader *header, BeamLineHitMan *blMan, AnaInfo *anaInfo)
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

  for( int i=0; i<anaInfo->nCDS(CDS_Kaon); i++ ){
    CDSInfo *km=anaInfo->CDS(CDS_Kaon, i);
    TLorentzVector km_lmom=km->lmom();
    double im=(fn_lmom+km_lmom).M();
    double mm=(kn_lmom-km_lmom).M();

    if( GeomTools::GetID(km->vertexBeam())==CID_Fiducial ){
      MyHistTools::fillTH("KNkm_MM_Nkm_IM", mm, im);
      MyHistTools::fillTH("Nkm_IM", im);
    }
  }
}
