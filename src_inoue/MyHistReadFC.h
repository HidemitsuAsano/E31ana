#ifndef MYHISTREADFC_HH
#define MYHISTREADFC_HH 1

#include "AnaInfo.h"
#include "MyHistTools.h"
#include "MyAnaTools.h"
#include "MyAnaParam.h"

void initHistReadFC();
void fillHistReadFC(EventHeader *hader, ConfMan *conf, BeamLineHitMan *blMan, AnaInfo *info);

#endif
