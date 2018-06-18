#ifndef MYHISTREADKNPIPI_H
#define MYHISTREADKNPIPI_H 1

#include "AnaInfo.h"
#include "MyTools.h"
#include "MyHistTools.h"
#include "BeamLineHitMan.h"
#include "MyAnaTools.h"
#include "CDSTrackingMan.h"
#include "EventHeader.h"
#include "BeamLineHitMan.h"

//#include "MyParam.h"
#include "MyParam_H2.h"

void initHistReadKNpipi_H2();
void fillHistReadKNpipi_H2(EventHeader *header, ConfMan *conf, BeamLineHitMan *blMan, CDSHitMan *cdsMan, CDSTrackingMan *cdstrackMan, AnaInfo *anaInfo);

#endif
