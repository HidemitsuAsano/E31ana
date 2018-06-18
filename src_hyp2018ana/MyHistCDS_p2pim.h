#ifndef MYHISTCDS_P2PIM_H
#define MYHISTCDS_P2PIM_H 1

#include "AnaInfo.h"
#include "MyHistTools.h"

void initHistCDS_p2pim();
void fillCDS_p2pim(ConfMan *conf, EventHeader *header, BeamLineHitMan *blMan, CDSHitMan *cdsMan, CDSTrackingMan *cdstrackMan, 
		   AnaInfo *anaInfo, ForwardChargeInfo *fcInfo);


#endif
