#ifndef MYHISTMCACC_H
#define MYHISTMCACC_H 1

#include "MyMCTools.h"
#include "MyAnaTools.h"
#include "MyHistTools.h"

#include "KnuclRootData.h"
#include "AnaInfo.h"
#include "MyAnaParam.h"

void initHistMCAcc();
void fillHistMCAcc_fp(DetectorData *detData, MCData* mcData, ReactionData *reacData, AnaInfo *anaInfo,
		      BeamLineHitMan *blMan, CDSTrackingMan *cdstrackMan);

void fillHistMCAcc_fn(DetectorData *detData, MCData* mcData, ReactionData *reacData, AnaInfo *anaInfo,
		      BeamLineHitMan *blMan, CDSTrackingMan *cdstrackMan);

#endif
