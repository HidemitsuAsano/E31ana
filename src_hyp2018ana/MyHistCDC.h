#ifndef MYHISTCDC_HH
#define MYHISTCDC_HH 1

#include "CDSHitMan.h"
#include "MyTools.h"
#include "MyHistTools.h"

void initHistCDC(ConfMan *conf);
void fillCDC(CDSHitMan *cdsMan);
void fillCDC(CDSTrackingMan *cdstrackMan, CDSHitMan *cdsMan);

#endif
