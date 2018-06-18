#ifndef MYHISTREADKPPIMPIM_HH
#define MYHISTREADKPPIMPIM_HH 1

#include "ConfMan.h"
#include "MyTools.h"
#include "MyHistTools.h"
#include "MyAnaTools.h"
#include "MyAnaParam.h"

void initHistReadKPpimpim();
void fillHistReadKPpimpim(EventHeader *header, BeamLineHitMan *blMan, CDSHitMan *cdsMan, CDSTrackingMan *cdstrackMan, AnaInfo *anaInfo); 

#endif
