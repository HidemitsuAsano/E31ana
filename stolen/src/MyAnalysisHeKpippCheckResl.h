// MyAnalysisHeKpippCheckResl.h

#ifndef MyAnalysisHeKpippCheckResl_h
#define MyAnalysisHeKpippCheckResl_h 1

#include "ConfMan.h"
#include "EventHeader.h"
#include "BeamLineHitMan.h"
#include "BeamLineTrackMan.h"
#include "CDSHitMan.h"
#include "CDSTrackingMan.h"
#include "Particle.h"
#include "KnuclRootData.h"

#include "TFile.h"
#include "TString.h"
#include "TH1.h"
#include "TH2.h"
#include "ELossTools.h"
#include "TrackTools.h"
#include "MathTools.h"
#include "/home/had/tyamaga/private/tmp/source/KinFitter/TKinFitter.h"
#include "/home/had/tyamaga/private/tmp/source/KinFitter/TFitParticlePxPyPz.h"
#include "/home/had/tyamaga/private/tmp/source/KinFitter/TFitConstraintM.h"
#include "/home/had/tyamaga/private/tmp/source/KinFitter/TFitConstraintEp.h"
#include <TDatabasePDG.h>

class MyAnalysisHeKpippCheckResl
{
  public:
    MyAnalysisHeKpippCheckResl(TFile*,ConfMan*);
    ~MyAnalysisHeKpippCheckResl();

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
    double lamll, lamul;
    double sbllamll, sbllamul;
    double sbulamll, sbulamul;


  private:
    TFile* rtFile;
    void CutCondition();
    bool Initialize(ConfMan*);
    bool FillHist(TString,double,int weight=1);
    bool FillHist(TString,TString,int weight=1);
    bool FillHist(TString,double,double,int weight=1);
    bool FillHist(TString,TString,TString,int weight=1);

  public:
    bool DoAnalysis(ConfMan*,EventHeader*,BeamLineHitMan*,BeamLineTrackMan*,CDSHitMan*,CDSTrackingMan*,Particle*, MCData*, ReactionData*, DetectorData*);
    void Clear();
};

#endif


