// MyAnalysisDetectorCheck.h

#ifndef MyAnalysisDetectorCheck_h
#define MyAnalysisDetectorCheck_h 1

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

class MyAnalysisDetectorCheck
{
  public:
    MyAnalysisDetectorCheck(TFile*,ConfMan*);
    ~MyAnalysisDetectorCheck();

  private:
    double k0ll, k0ul;
    double sblk0ll, sblk0ul;
    double sbuk0ll, sbuk0ul;
    double mnll, mnul;
    double sblmnll, sblmnul;
    double sbumnll, sbumnul;
    double spll, spul;
    double sblspll, sblspul;
    double sbuspll, sbuspul;
    double smll, smul;
    double sblsmll, sblsmul;
    double sbusmll, sbusmul;
    double mk0ll, mk0ul;
    double sblmk0ll, sblmk0ul;
    double sbumk0ll, sbumk0ul;

  private:
    TFile* rtFile;
		ConfMan* confFile;
    void CutCondition();
    bool Initialize(ConfMan*);
    bool FillHist(TString,double,int weight=1);
    bool FillHist(TString,TString,int weight=1);
    bool FillHist(TString,double,double,int weight=1);
    bool FillHist(TString,TString,TString,int weight=1);
    double CalcTimeOffs(HodoscopeLikeHit*,Particle*);
    double CalcTimeSubOffs(HodoscopeLikeHit*,Particle*);

  public:
    bool DoAnalysis(ConfMan*,EventHeader*,BeamLineHitMan*,BeamLineTrackMan*);
    void Clear();
};

#endif


