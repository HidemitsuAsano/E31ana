#ifndef MYHISTREADNC_HH
#define MYHISTREADNC_HH 1

#include "AnaInfo.h"
#include "MyHistTools.h"
#include "MyAnaTools.h"

void initHistReadNC();
void fillHistReadNC(EventHeader *header, BeamLineHitMan *blMan, AnaInfo *anaInfo);

#endif
