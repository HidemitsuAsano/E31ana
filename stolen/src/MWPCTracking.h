// MWPCTracking.h

#ifndef MWPCTracking_h
#define MWPCTracking_h 1

#include "ConfMan.h"
#include "EventHeader.h"
#include "BeamLineHitMan.h"
#include "BeamLineTrackMan.h"
#include "CDSHitMan.h"
#include "CDSTrackingMan.h"
#include "Particle.h"
#include "HodoscopeLikeHit.h"

#include "TFile.h"
#include "TString.h"
#include "TH1.h"
#include "TH2.h"
#include "ELossTools.h"
#include "TrackTools.h"
#include "MathTools.h"

class MWPCTracking
{
  public:
    MWPCTracking(ConfMan*,BeamLineHitMan*);
    ~MWPCTracking();

  private:
    ConfMan* conf;
    BeamLineHitMan* blMan;

  private:
		ConfMan* confFile;
    bool Initialize(ConfMan*,BeamLineHitMan*);

  public:
    bool DoTracking(int);
    void Clear();
};

#endif


