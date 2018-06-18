#ifndef USERTOOLS_HH
#define USERTOOLS_HH 1

#include "BeamLineHitMan.h"
#include "BeamLineTrackMan.h"
#include "BeamSpectrometer.h"

namespace UserTools{
  bool isGoodBeam(BeamLineHitMan *blMan, BeamLineTrackMan *bltrackMan, BeamSpectrometer *beamSpec);
  LocalTrack* getBPCTrack(BeamLineTrackMan *bltrackMan);
};

#endif
