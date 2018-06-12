#ifndef MYHISTBEAMANA_HH
#define MYHISTBEAMANA_HH 1

#include "BeamSpectrometer.h"
#include "MyHistTools.h"

void initHistBeamAna();

void fillBeamProf(EventHeader *header, BeamLineTrackMan *bltrackMan);
void fillBLC2BPC(BeamLineTrackMan *bltrackMan);
void fillBeamSpectrometer(BeamSpectrometer *beam);

void fillBHDT0offset(BeamInfo *info, BeamLineHitMan *blMan);

#endif
