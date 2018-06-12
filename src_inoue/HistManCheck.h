#ifndef HistManCheck_h
#define HistManCheck_h 1

#include "ConfMan.h"
#include "EventHeader.h"
#include "BeamLineHitMan.h"
#include "BeamLineTrackMan.h"
#include "BeamSpectrometer.h"
#include "CDSTrackingMan.h"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "ELossTools.h"
#include "TrackTools.h"

#include "MyTools.h"
#include "SlewingData.h"

class HistManCheck
{
 public:
  HistManCheck();
  virtual ~HistManCheck(){};

};
#endif
