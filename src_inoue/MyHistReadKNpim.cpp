#include "MyHistReadKNpim.h"

using namespace std;

void initHistReadKNpim()
{
  new TH2F("KNpim_MM_Npim_IM", "d(K^{-}, n #pi^{-})\"X\" vs n #pi^{-} IM", 2000, 1.0, 2.0, 1000, 1.0, 2.0);
  new TH1F("Npim_IM",      "n #pi^{-} IM", 1000, 1.0, 2.0);
  new TH1F("Npim_IM_true", "n #pi^{-} IM", 1000, 1.0, 2.0);
}

void fillHistReadKNpim(EventHeader *header, BeamLineHitMan *blMan, AnaInfo *anaInfo)
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

  for( int i=0; i<anaInfo->nCDS(CDS_PiMinus); i++ ){
    CDSInfo *pim=anaInfo->CDS(CDS_PiMinus, i);
    TLorentzVector pim_lmom=pim->lmom();
    double im=(fn_lmom+pim_lmom).M();
    double mm=(kn_lmom-pim_lmom).M();

    if( GeomTools::GetID(pim->vertexBeam())==CID_Fiducial ){
      MyHistTools::fillTH("KNpim_MM_Npim_IM", mm, im);
      MyHistTools::fillTH("Npim_IM", im);
      if( Sm_MIN<im &&  im<Sm_MAX ){
	MyHistTools::fillTH("Npim_IM_true", im);
      }
    }
  }
}
