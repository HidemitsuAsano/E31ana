#ifndef MYHISTREADCDS_LPIM_HH
#define MYHISTREADCDS_LPIM_HH 1

#include "AnaInfo.h"
#include "MyTools.h"
#include "MyHistTools.h"
#include "BeamLineHitMan.h"
#include "MyAnaTools.h"
#include "MyAnaParam.h"

void initHistReadCDS_Lpim();
void fillHistReadCDS_Lpim(EventHeader *header, BeamLineHitMan *blMan, AnaInfo *anaInfo);

#endif
