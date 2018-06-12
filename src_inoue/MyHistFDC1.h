#ifndef MYHISTFDC1_HH
#define MYHISTFDC1_HH 1

#include "GlobalVariables.h"
#include "BeamLineHitMan.h"
#include "MyParam.h"
#include "MyHistTools.h"

void initHistFDC1();
void fillFDC1(BeamLineHitMan *blMan);

void fillFDC1(LocalTrack *track);
void fillFDC1(EventHeader *header, BeamLineHitMan *blMan, BeamLineTrackMan *bltrackMan);

#endif
