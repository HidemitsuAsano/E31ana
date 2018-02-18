// MyAnalysisCDS.h

#ifndef MyAnalysisCDS_h
#define MyAnalysisCDS_h 1

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

class MyAnalysisCDS
{
  public:
    MyAnalysisCDS(TFile*,ConfMan*,bool simflag=false);
    ~MyAnalysisCDS();

  private:
    int BeamPID;
    double BeamMom;
    TLorentzVector BeamL;
    int T0Segment;
    double T0Timing;
    LocalTrack* trackblc1;
    LocalTrack* trackblc2;
    LocalTrack* trackbpc;
    bool SIMULATION;

  private:
    bool Initialize(TFile*,ConfMan*);
    bool FillHist(TString,double);
    bool FillHist(TString,double,double);

  public:
    bool DoAnalysis(ConfMan*,EventHeader*,BeamLineHitMan*,BeamLineTrackMan*,CDSHitMan*,CDSTrackingMan*,Particle*);
    bool DoAnalysisSim(ConfMan*,EventHeader*,BeamLineHitMan*,BeamLineTrackMan*,CDSHitMan*,CDSTrackingMan*,Particle*);
    void Clear();
    void SetBeamPID(int pid) { BeamPID = pid; };
    void SetBeamMom(double mom) { BeamMom = mom; };
    void SetT0Timing(double time) { T0Timing = time; };
    void SetBeamL(TLorentzVector l) { BeamL = l; };
    void SetBLC1(LocalTrack* track) { trackblc1 = track; };
    void SetBLC2(LocalTrack* track) { trackblc2 = track; };
    void SetBPC(LocalTrack* track)  { trackbpc  = track; };
};

#endif


