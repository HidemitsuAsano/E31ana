// MyAnalysisCDS1Track.h

#ifndef MyAnalysisCDS1Track_h
#define MyAnalysisCDS1Track_h 1

#include "ConfMan.h"
#include "EventHeader.h"
#include "BeamLineHitMan.h"
#include "BeamLineTrackMan.h"
#include "CDSHitMan.h"
#include "CDSTrackingMan.h"

#include "TFile.h"
#include "TString.h"
#include "TH1.h"
#include "TH2.h"
#include "ELossTools.h"
#include "TrackTools.h"
#include "MathTools.h"

class MyAnalysisCDS1Track
{
  public:
    MyAnalysisCDS1Track(TFile*,ConfMan*);
    ~MyAnalysisCDS1Track();

  private:
    int BeamPID;
    double BeamMom;
    TLorentzVector BeamL;
    int T0Segment;
    double T0Timing;
    LocalTrack* trackblc1;
    LocalTrack* trackblc2;
    LocalTrack* trackbpc;

  private:
    bool Initialize(TFile*,ConfMan*);
    bool FillHist(TString,double);
    bool FillHist(TString,double,double);

  public:
    bool DoAnalysis(ConfMan*,EventHeader*,BeamLineHitMan*,BeamLineTrackMan*,CDSHitMan*,CDSTrackingMan*);
    bool Clear();
    void SetBeamPID(int pid) { BeamPID = pid; };
    void SetBeamMom(double mom) { BeamMom = mom; };
    void SetBeamL(TLorentzVector l) { BeamL = l; };
    void SetBLC1(LocalTrack* track) { trackblc1 = track; };
    void SetBLC2(LocalTrack* track) { trackblc2 = track; };
    void SetBPC(LocalTrack* track)  { trackbpc  = track; };
    void SetT0Timing(double time)   { T0Timing  = time; };
};

#endif


