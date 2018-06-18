#ifndef MYHISTREADCDS_HH
#define MYHISTREADCDS_HH

#include "AnaInfo.h"
#include "MyTools.h"
#include "MyHistTools.h"
#include "BeamLineHitMan.h"
#include "MyAnaTools.h"

#include "MyAnaParam.h"

void initHistReadCDS();
void fillHistReadCDS(AnaInfo *anaInfo, CDSHitMan *cdsMan, CDSTrackingMan *cdstrackMan);

#endif
