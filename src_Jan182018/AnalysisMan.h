/* AnalysisMan.h */
#ifndef AnalysisMan_h
#define AnalysisMan_h 1

#include "Common.h"

class AnalysisMan : public TObject
{
 private:
  AnalysisMan();

 public:
  AnalysisMan(ConfMan *conf);

 private:
  BeamLineChamber *BPC;

 public:
  int Execute(ConfMan *conf, BeamLineHitMan *blMan, CDSHitMan *cdsMan);
  void Clear();

  ClassDef(AnalysisMan, 1)
};

#endif
