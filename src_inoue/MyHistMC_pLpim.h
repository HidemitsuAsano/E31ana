#ifndef MYHISTMC_PLPIM_H
#define MYHISTMC_PLPIM_H 1

#include "AnaInfo.h"
#include "KnuclRootData.h"
#include "MyHistTools.h"
#include "MyAnaTools.h"

void initHistMC_pLpim();

void fillHistMC_pLpim(DetectorData *detData, MCData* mcData, ReactionData *reacData, AnaInfo *anaInfo,
                      BeamLineHitMan *blMan, CDSTrackingMan *cdstrackMan);
#endif
