#ifndef MYHISTTEMP_HH
#define MYHISTTEMP_HH 1

#include "GlobalVariables.h"
#include "BeamLineHitMan.h"
#include "MyParam.h"
#include "MyHistTools.h"
#include "AnaInfo.h"

void initHistTemp();
void fillHistTemp(ConfMan *conf, EventHeader *header, BeamLineHitMan *blMan, CDSHitMan *cdsMan, CDSTrackingMan *cdstrackMan, AnaInfo *anaInfo);



#endif
