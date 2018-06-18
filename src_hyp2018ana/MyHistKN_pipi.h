#ifndef MYHISTKN_PIPI_H
#define MYHISTKN_PIPI_H 1

#include "AnaInfo.h"
#include "MyHistTools.h"

void initHistKN_pipi();
void fillHistKN_pipi(EventHeader *header, BeamLineHitMan *blMan, CDSHitMan *cdsMan, CDSTrackingMan *cdstrackMan, AnaInfo *anaInfo);

#endif
