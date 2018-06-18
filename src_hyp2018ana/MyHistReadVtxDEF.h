#ifndef MYHISTREADVTXDEF_H
#define MYHISTREADVTXDEF_H 1

#include "EventHeader.h"
#include "ConfMan.h"
#include "AnaInfo.h"
#include "BeamLineTrackMan.h"
#include "CDSTrackingMan.h"

#include "MyAnaParam.h"
#include "MyHistTools.h"
#include "MyAnaTools.h"

void initHistReadVtxDEF();

void fillHistReadVtxDEF(EventHeader *header, ConfMan *conf, AnaInfo *anaInfo,
			BeamLineHitMan *blMan, BeamLineTrackMan *bltrackMan,
			CDSHitMan *cdsMan, CDSTrackingMan *cdstrackMan);

#endif
