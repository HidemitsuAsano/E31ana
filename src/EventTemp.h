#ifndef EVENTTEMP_h 
#define EVENTTEMP_h 1

#include <cstdio>

#include <TFile.h>
#include <TTree.h>
#include <TGraphAsymmErrors.h>

#include "ConfMan.h"
#include "TKO.h"
#include "CDSHitMan.h"
#include "BeamLineHitMan.h"
#include "EventHeader.h"
#include "ScalerMan.h"

class EventTemp
{
 protected:
  ConfMan *confMan;

  int Event_Number;
  int VEvent_Number;
  int Block_Event_Number;

 public:
  EventTemp();
  virtual ~EventTemp(){};
  
  virtual void Initialize( ConfMan *conf )=0;
  virtual void USca( int nsca=0, unsigned int *sca=0 )=0;
  virtual void UTime( int time=0 ){};
  virtual bool UAna( TKOHitCollection *tko )=0;
  virtual void Finalize()=0;
};

#endif
