#ifndef MYHISTReduction_HH
#define MYHISTReduction_HH 1

#include "GlobalVariables.h"
#include "BeamLineHitMan.h"
#include "MyParam.h"
#include "MyHistTools.h"

void initHistReduction();
void fillReduction(BeamLineHitMan *blMan, BeamLineTrackMan *bltrackMan, BeamSpectrometer *spec, EventHeader *header);

#endif
