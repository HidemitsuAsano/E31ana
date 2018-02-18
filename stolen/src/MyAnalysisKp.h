// MyAnalysisKp.h

#ifndef MyAnalysisKp_h
#define MyAnalysisKp_h 1

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

class MyAnalysisKp
{
  public:
    MyAnalysisKp(TFile*,ConfMan*);
    ~MyAnalysisKp();

  private:

  private:
    bool Initialize(TFile*,ConfMan*);
    bool FillHist(TString,double,int weight=1);
    bool FillHist(TString,TString,int weight=1);
    bool FillHist(TString,double,double,int weight=1);
    bool FillHist(TString,TString,TString,int weight=1);

  public:
    bool DoAnalysis(ConfMan*,EventHeader*,BeamLineHitMan*,BeamLineTrackMan*,CDSHitMan*,CDSTrackingMan*,Particle*);
    void Clear();
};

#endif


