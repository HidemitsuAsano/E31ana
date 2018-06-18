#ifndef MYHISTREADKNPIP_HH
#define MYHISTREADKNPIP_HH

#include "AnaInfo.h"
#include "MyHistTools.h"
#include "MyAnaTools.h"

void initHistReadKNpip();
void fillHistReadKNpip(EventHeader *hader, BeamLineHitMan *blMan, AnaInfo *info);

#endif
