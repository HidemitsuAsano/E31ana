#ifndef MYHISTBLC1A_HH
#define MYHISTBLC1A_HH 1

#include "GlobalVariables.h"
#include "BeamLineHitMan.h"
#include "MyParam.h"
#include "MyHistTools.h"

void initHistBLC1a();
void fillBLC1a(BeamLineHitMan *blMan);

void fillBLC1a(LocalTrack *track);

#endif
