#ifndef MYHISTREADKNPIPI_D2_HH
#define MYHISTREADKNPIPI_D2_HH 1

#include "AnaInfo.h"
#include "MyTools.h"
#include "MyHistTools.h"
#include "BeamLineHitMan.h"
#include "MyAnaTools.h"
#include "CDSTrackingMan.h"
#include "MyParam.h"

void initHistReadKNpipi_D2(EventHeader *header, AnaInfo *anaInfo);
void fillHistReadKNpipi_D2(EventHeader *header, BeamLineHitMan *blMan, CDSHitMan *cdsMan, CDSTrackingMan *cdstrackMan, AnaInfo *anaInfo);

#endif
