#ifndef MYHISTMCFC_H
#define MYHISTNCFC_H 1

#include "MyHistTools.h"
#include "KnuclRootData.h"
#include "MyTools.h"
#include "MyMCTools.h"

void initHistMCFC();
void fillHistMCFC(ConfMan *conf, AnaInfo *anaInfo, CDSHitMan *cdsMan, BeamLineHitMan *blMan, 
		  CDSTrackingMan *cdstrackMan, BeamLineTrackMan *bltrackMan,
		  ReactionData *reacData, MCData *mcData, DetectorData *detData);

#endif
