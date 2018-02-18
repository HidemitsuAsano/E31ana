// MyAnalysisBase.h
// 2014.12.22
// T. Yamaga

#ifndef MyAnalysisBase_h
#define MyAnalysisBase_h

#include <TFile.h>
#include <TTree.h>
#include "ConfMan.h"
#include "CDSHitMan.h"
#include "BeamLineHitMan.h"
#include "EventHeader.h"
#include "ScalerMan.h"
#include "BeamLineTrackMan.h"
#include "CDSTrackingMan.h"

#include "Tools.h"

class MyAnalysisBase
{
  protected:
    ConfMan* confMan;
    CDSTrackingMan* cdsTrackMan;
    BeamLineTrackMan* blTrackMan; 

    TFile* rtFile;

  public:
    MyAnalysisBase();
    virtual ~MyAnalysisBase(){};

    void SetTrackMan(CDSTrackingMan*,BeamLineTrackMan*);
    void SetTrackMan(CDSTrackingMan*);
    void SetTrackMan(BeamLineTrackMan*);
    void Initialize(ConfMan*);
    void Finalize();
    virtual bool DoScaler(ScalerMan*)=0;
    virtual bool DoAnalysis(CDSHitMan*,BeamLineHitMan*)=0;
};

#endif
