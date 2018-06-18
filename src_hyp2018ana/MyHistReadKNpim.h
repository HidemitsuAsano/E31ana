#ifndef MYHISTREADKNPIM_HH
#define MYHISTREADKNPIM_HH

#include "AnaInfo.h"
#include "MyHistTools.h"
#include "MyAnaTools.h"

void initHistReadKNpim();
void fillHistReadKNpim(EventHeader *hader, BeamLineHitMan *blMan, AnaInfo *info);

#endif
