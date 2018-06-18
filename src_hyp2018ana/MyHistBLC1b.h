#ifndef MYHISTBLC1B_HH
#define MYHISTBLC1B_HH 1

#include "GlobalVariables.h"
#include "BeamLineHitMan.h"
#include "MyParam.h"
#include "MyHistTools.h"

void initHistBLC1b();
void fillBLC1b(BeamLineHitMan *blMan);

void fillBLC1b(LocalTrack *track);

#endif
