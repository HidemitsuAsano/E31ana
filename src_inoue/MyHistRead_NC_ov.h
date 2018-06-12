#ifndef MYHISTREAED_NC_OV
#define MYHISTREAD_NC_OV 1

#include "AnaInfo.h"
#include "MyTools.h"
#include "MyHistTools.h"

void initHistRead_NC_ov();
void fillHistRead_NC_ov(ConfMan *conf, EventHeader *header, AnaInfo *anaInfo, BeamLineHitMan *blMan, CDSHitMan *cdsMan ,CDSTrackingMan *cdstrackMan);

#endif
