#ifndef MYHISTREADKPPIMPIM_D2_HH
#define MYHISTREADKPPIMPIM_D2_HH 1

#include "ConfMan.h"
#include "MyTools.h"
#include "MyHistTools.h"
#include "MyAnaTools.h"
#include "MyAnaParam.h"

void initHistReadKPpimpim_D2();
void fillHistReadKPpimpim_D2(EventHeader *header, ConfMan *conf, BeamLineHitMan *blMan, CDSHitMan *cdsMan, CDSTrackingMan *cdstrackMan, AnaInfo *anaInfo); 

#endif
