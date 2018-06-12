#include "MyHistCDS.h"

using namespace std;

void initHistCDS()
{
  new TH2F("CDS_mass2_mom", "CDS mass2 mom", 1000, -0.5, 9.5, 1000, -1.0, 1.0);
}

void fillCDS(CDSInfo *info, CDSHitMan *cdsMan, CDSTrackingMan *cdstrackingMan)
{
  if( info->mass2()!=DEFAULTD && info->mom()!=DEFAULTD ){
    MyHistTools::fillTH("CDS_mass2_mom", info->mass2(), info->mom());
  }
}
