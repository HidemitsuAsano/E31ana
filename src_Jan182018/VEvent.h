// VEvent.h

#ifndef Event_h
#define Event_h 1

#include <cstdio>
#include <vector>

#include "TFile.h"
#include "TTree.h"

#include "ConfMan.h"
#include "DataForm.h"

#include "TKO.h"
#include "CDSHitMan.h"
#include "EventHeader.h"
#include "EventAlloc.h"
#include "EventTemp.h"

class VEvent
{
 private:
  //protected:
  unsigned int *Buf;
  std::vector <unsigned int *> EventHeaderPointerContainer;
  std::vector <unsigned int *> EventBufferPointerContainer;

  TKOHitCollection *tkoCol;
  EventTemp *evTemp;
  EventAlloc *evAlloc;

 public:
  VEvent( unsigned int *buf )
    : Buf(buf)
    {}
  VEvent( unsigned int *buf, int evnum )
    : Buf(buf), EventNumLast(evnum)
    {}

  VEvent()
    : Buf(0)
    {}

    ~VEvent()
      {}
    
 private:
  //protected:
  int RecordType;
  int RunNum, BlockEventNum, EventNum, Size;
  int EventNumLast, SerialNo;

  int NumEventInBlock; // set a value of first event, and used for checking the data status.

 public:
  EventTemp *EventAllocator();

  int GetRunNumber() const { return RunNum; }
  int GetEventNumber() const { return EventNum; }
  int GetSize() const { return Size; }
  int GetRecordTpye() const { return RecordType; }

  //bool Processing();
  bool Processing3(ConfMan *confMan);
  bool EventCheck();
  //virtual bool ReadEvent( FILE *fp ) { return false;}

  void SetBuffer( unsigned int *buf, int evnum );

 protected:
  bool ProcessingNormal(ConfMan *confMan);
  bool ProcessingBegin(ConfMan *confMan);
  bool ProcessingPause();
  bool ProcessingResume();
  bool ProcessingEnd(ConfMan *confMan);
  bool ProcessingUnknown();
  void PrintComment();
};

#endif
