// MyAnalysisdKCDH3.h

#ifndef MyAnalysisdKCDH3_h
#define MyAnalysisdKCDH3_h 1

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

class MyAnalysisdKCDH3
{
  public:
    MyAnalysisdKCDH3(TFile*,ConfMan*);
    ~MyAnalysisdKCDH3();

  private:
    double k0ll, k0ul;
    double sblk0ll, sblk0ul;
    double sbuk0ll, sbuk0ul;
    double mnll, mnul;
    double sblmnll, sblmnul;
    double sbumnll, sbumnul;
    double mpll, mpul;
    double sblmpll, sblmpul;
    double sbumpll, sbumpul;
    double spll, spul;
    double sblspll, sblspul;
    double sbuspll, sbuspul;
    double smll, smul;
    double sblsmll, sblsmul;
    double sbusmll, sbusmul;
    double lamll, lamul;
    double sbllamll, sbllamul;
    double sbulamll, sbulamul;
    double flamll, flamul;
    double sblflamll, sblflamul;
    double sbuflamll, sbuflamul;


  private:
    TFile* rtFile;
    void CutCondition();
    bool Initialize(ConfMan*);
    bool FillHist(TString,double,int weight=1);
    bool FillHist(TString,TString,int weight=1);
    bool FillHist(TString,double,double,int weight=1);
    bool FillHist(TString,TString,TString,int weight=1);
    bool DoNeutralAna(Particle*,pBeam*,pCDS*,TVector3);
    bool DoChargedAna(Particle*,pBeam*,pCDS*,TVector3);

  public:
    bool DoAnalysis(ConfMan*,EventHeader*,BeamLineHitMan*,BeamLineTrackMan*,CDSHitMan*,CDSTrackingMan*,Particle*);
    void Clear();

};

#endif


