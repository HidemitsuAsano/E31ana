#ifndef MYHISTREADKNPIPI_HH
#define MYHISTREADKNPIPI_HH 1

#include "AnaInfo.h"
#include "MyTools.h"
#include "MyHistTools.h"
#include "BeamLineHitMan.h"
#include "MyAnaTools.h"
#include "CDSTrackingMan.h"
#include "MyParam.h"

void initHistReadKNpipi(AnaInfo *anaInfox);
void fillHistReadKNpipi(EventHeader *header, BeamLineHitMan *blMan, CDSHitMan *cdsMan, CDSTrackingMan *cdstrackMan, AnaInfo *anaInfo);

#endif
