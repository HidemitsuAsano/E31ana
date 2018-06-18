#include "EventAnalysisReadCDS.h"

using namespace std;

bool EventAnalysisReadAna::UAna()
{
  if( !MyAnaTools::goodBeam(anaInfo) ) return true;

  if( cdsFileMatching() ){
    fillHistCalibCDS(header, cdsMan, cdstrackMan, anaInfo);
  }
  return true;
}

void EventAnalysisReadAna::InitializeHistogram()
{
  initHistCalibCDS();
}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysisReadAna *event = new EventAnalysisReadAna();
  return (EventTemp*)event;
}
