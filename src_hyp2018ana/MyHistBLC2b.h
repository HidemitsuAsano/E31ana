#ifndef MYHISTBLC2B_HH
#define MYHISTBLC2B_HH 1

#include "GlobalVariables.h"
#include "BeamLineHitMan.h"
#include "MyParam.h"
#include "MyHistTools.h"

void initHistBLC2b();
void fillBLC2b(BeamLineHitMan *blMan);

void fillBLC2b(LocalTrack *track);

#endif
