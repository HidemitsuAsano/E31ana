#ifndef MYHISTCDS_HH
#define MYHISTCDS_HH 1

#include "GlobalVariables.h"
#include "BeamLineHitMan.h"
#include "MyParam.h"
#include "MyHistTools.h"

void initHistCDS();
void fillCDS(CDSInfo *info, CDSHitMan *cdsMan, CDSTrackingMan *cdstrackMan);

#endif
