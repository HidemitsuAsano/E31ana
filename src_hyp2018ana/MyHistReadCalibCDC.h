#ifndef MYHISTREADCALIBCDC_HH
#define MYHISTREADCALIBCDC_HH 1

#include "AnaInfo.h"
#include "MyHistTools.h"
#include "MyAnaTools.h"

#include "CDSTrackingMan.h"

void initHistReadCalibCDC(ConfMan *conf);
void fillHistReadCalibCDC(ConfMan *conf, BeamLineHitMan *blMan, CDSHitMan *cdsMan, CDSTrackingMan *cdstrackMan, AnaInfo *anaInfo);

#endif
