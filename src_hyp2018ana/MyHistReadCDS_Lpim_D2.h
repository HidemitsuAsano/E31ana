#ifndef MYHISTREADCDS_LPIM_D2_HH
#define MYHISTREADCDS_LPIM_D2_HH 1

#include "AnaInfo.h"
#include "MyTools.h"
#include "MyHistTools.h"
#include "BeamLineHitMan.h"
#include "MyAnaTools.h"
#include "MyAnaParam.h"

void initHistReadCDS_Lpim_D2();
void fillHistReadCDS_Lpim_D2(EventHeader *header, BeamLineHitMan *blMan, AnaInfo *anaInfo);

#endif
