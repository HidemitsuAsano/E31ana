#include "MyHistRead_NC_ov.h"

using namespace std;

void initHistRead_NC_ov()
{
  new TH1F("NC_overbeta_wo_lay1", "NC 1/#beta", 50000, 0.0, 5.0);
  new TH1F("NC_overkill", "NC overkill check", 10, 0, 10);

  new TH1F("BVC_tof_ok", "T0-BVC tof over kill", 1000, -500, 500);
  new TH1F("CVC_tof_ok", "T0-CVC tof over kill", 1000, -500, 500);
}

void fillHistRead_NC_ov(ConfMan *conf, EventHeader *header, AnaInfo *anaInfo, BeamLineHitMan *blMan, CDSHitMan *cdsMan ,CDSTrackingMan *cdstrackMan)
{
  if( anaInfo->nBeam()!=1 || !anaInfo->beam(0)->flag() ) return;
  if( !anaInfo->minDCA() ) return;
  std::vector<HodoscopeLikeHit*> CVChits=MyTools::getHodo(blMan, CID_CVC);
  std::vector<HodoscopeLikeHit*> BVChits=MyTools::getHodo(blMan, CID_BVC);
  BeamInfo *beam=anaInfo->beam(0);

  if( GeomTools::GetID(anaInfo->minDCA()->vertexBeam())!=CID_Fiducial ) return;

  if( anaInfo->nFNeutral()==1 ){
    ForwardNeutralInfo *fnInfo=anaInfo->forwardNeutral(0);
    HodoscopeLikeHit *hit=fnInfo->NC(blMan);

    if( !header->IsTrig(Trig_Neutral) ) return;
    if( hit->seg()<17 ) return;
    MyHistTools::fillTH("NC_overbeta_wo_lay1", 1./fnInfo->beta());

    if( 1./fnInfo->beta()<1.25 || 1.3<1./fnInfo->beta() ) return;
    MyHistTools::fillTH("NC_overkill", 0);
    if( BVChits.size()>0 ){
      for( int i=0; i<BVChits.size(); i++ ){
	double tof=BVChits[i]->ctmean()-beam->T0time();
	MyHistTools::fillTH("BVC_tof_ok", tof);
      }
      MyHistTools::fillTH("NC_overkill", 1);
    }
    if( CVChits.size()>0 ){
      for( int i=0; i<CVChits.size(); i++ ){
	double tof=CVChits[i]->ctmean()-beam->T0time();
	MyHistTools::fillTH("CVC_tof_ok", tof);
      }
      MyHistTools::fillTH("NC_overkill", 2);
    }
    if( BVChits.size()>0 && CVChits.size()>0 ) MyHistTools::fillTH("NC_overkill", 3);
    if( BVChits.size()>0 || CVChits.size()>0 ) MyHistTools::fillTH("NC_overkill", 4);
  }
}
