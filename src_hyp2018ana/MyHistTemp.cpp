#include "MyHistTemp.h"

void initHistTemp()
{
  new TH1F("nCDC_fn", "nCDC_fn", 20, -0.5, 19.5);

}

void fillHistTemp(ConfMan *conf, EventHeader *header, BeamLineHitMan *blMan, CDSHitMan *cdsMan, CDSTrackingMan *cdstrackMan, AnaInfo *anaInfo)
{
  if( anaInfo->nBeam()==1 && anaInfo->beam(0)->flag() ){
    if( header->IsTrig(Trig_Neutral) ){
      if( anaInfo->nFNeutral()==1 ){
	if( anaInfo->forwardNeutral(0)->pid()==F_Neutron ){
	  MyHistTools::fillTH("nCDC_fn", anaInfo->nCDS());

	}
      }
    }
  }
}
