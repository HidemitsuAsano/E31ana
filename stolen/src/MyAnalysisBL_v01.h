// MyAnalysisBL.h

#ifndef MyAnalysisBL_h
#define MyAnalysisBL_h 1

#include "ConfMan.h"
#include "EventHeader.h"
#include "BeamLineHitMan.h"
#include "BeamLineTrackMan.h"
#include "BeamSpectrometer.h"
#include "CDSTrackingMan.h"

#include "TFile.h"
#include "TString.h"
#include "TH1.h"
#include "TH2.h"
#include "ELossTools.h"
#include "TrackTools.h"
#include "MathTools.h"

class MyAnalysisBL
{
  public:
    MyAnalysisBL(TFile*,ConfMan*);
    ~MyAnalysisBL();

  private:
    double T0Timing;
    int T0Segment;
    int BHDSegment;
    int BeamPID;
    double BeamMom;
    double BeamChi2;
    double BeamTOF;
    double BeamMass;
    TVector3 BeamP;
    TLorentzVector BeamL;
    BeamSpectrometer* D5;
    LocalTrack* trackblc1;
    LocalTrack* trackblc2;
    LocalTrack* trackbpc;

  private:
    bool Initialize(TFile*,ConfMan*);
    bool FillHist(TString,double);
    bool FillHist(TString,double,double);

  public:
    bool DoAnalysis(ConfMan*,EventHeader*,BeamLineHitMan*,BeamLineTrackMan*);
    bool Clear();
    int GetBeamPID() const { return BeamPID; };
    double GetBeamMom() const {return BeamMom; };
    double GetBeamMass() const {return BeamMass; };
    TVector3 GetBeamP() const {return BeamP; };
    TLorentzVector GetBeamL() const {return BeamL; };
    LocalTrack* GetBLC1() const {return trackblc1; };
    LocalTrack* GetBLC2() const {return trackblc2; };
    LocalTrack* GetBPC() const {return trackbpc; };
};

#endif


