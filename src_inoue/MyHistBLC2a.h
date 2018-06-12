#ifndef MYHISTBLC2A_HH
#define MYHISTBLC2A_HH 1

#include "GlobalVariables.h"
#include "BeamLineHitMan.h"
#include "MyParam.h"
#include "MyHistTools.h"

void initHistBLC2a();
void fillBLC2a(BeamLineHitMan *blMan);

void fillBLC2a(LocalTrack *track);

#endif
