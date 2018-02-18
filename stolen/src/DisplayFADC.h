#ifndef DISPLAY_h
#define DISPLAY_h 1

#include <vector>
#include <iostream>
#include <fstream>

#ifndef ROOT_TObject
#include "TObject.h"
#endif
#include <TSystem.h>
#include <TVirtualPad.h>
#include <TH2F.h>
#include <TGraph.h>
#include <TLine.h>

#include "ConfMan.h"
#include "SDDHitMan.h"

class DisplayFADC : public TObject
{
 public:
  DisplayFADC();
  virtual ~DisplayFADC() {};

 private:
  DisplayFADC( DisplayFADC & );

  TH2F *frameFADC;
  double XMIN,XMAX;
 public:
  bool Wait();

  void SetFADCFrame( double xmin=-80, double xmax=80, double ymin=-80, double ymax=80 );
  bool DrawFADCFrame( TVirtualPad *pad );
  bool DrawFADCFrameSep( TVirtualPad *pad );
  bool DrawFADCHit( TVirtualPad *pad, ConfMan *conf, SDDHitMan *sdd);
  bool DrawFADCHitSep( TVirtualPad *pad, ConfMan *conf, SDDHitMan *sdd);

  ClassDef( DisplayFADC, 1 );
};
#endif
