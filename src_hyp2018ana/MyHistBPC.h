#ifndef MYHISTBPC_HH
#define MYHISTBPC_HH 1

#include "GlobalVariables.h"
#include "BeamLineHitMan.h"
#include "MyParam.h"
#include "MyHistTools.h"

void initHistBPC();
void fillBPC(BeamLineHitMan *blMan);

void fillBPC(LocalTrack *track);
void fillBPC(BeamLineTrackMan *bltrackMan);

#endif
