// MyAnalysisFWDCharged.h

#ifndef MyAnalysisFWDCharged_h
#define MyAnalysisFWDCharged_h 1

#include "ConfMan.h"
#include "EventHeader.h"
#include "BeamLineHitMan.h"
#include "BeamLineTrackMan.h"
#include "CDSHitMan.h"
#include "CDSTrackingMan.h"
#include "Particle.h"

#include "TFile.h"
#include "TString.h"
#include "TH1.h"
#include "TH2.h"
#include "ELossTools.h"
#include "TrackTools.h"
#include "MathTools.h"

class MyAnalysisFWDCharged
{
  public:
    MyAnalysisFWDCharged(TFile*,ConfMan*,bool simflag=false);
    ~MyAnalysisFWDCharged();

  private:
    bool SIMULATION;

  private:
    bool Initialize(TFile*,ConfMan*);
    bool FillHist(TString,double);
    bool FillHist(TString,double,double);

  public:
    bool DoAnalysis(ConfMan*,EventHeader*,BeamLineHitMan*,BeamLineTrackMan*,CDSHitMan*,CDSTrackingMan*,Particle*);
    void Clear();
};

#endif


