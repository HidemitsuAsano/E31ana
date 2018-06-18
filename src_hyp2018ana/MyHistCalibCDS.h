#ifndef MYHISTCALIBCDS_HH
#define MYHISTCALIBCDS_HH 1

#include "AnaInfo.h"
#include "MyHistTools.h"
#include "MyAnaTools.h"

void initHistCalibCDS();
void fillHistCalibCDS(EventHeader *header, ConfMan *conf, CDSHitMan *cdsMan, CDSTrackingMan *cdstrackMan, AnaInfo *anaInfo);

#endif
