#include "EventAnalysisReadNC_ov.h"

using namespace std;

bool EventAnalysisReadNC_ov::UAna()
{
  if( !MyAnaTools::goodBeam(anaInfo) ) return true;

  if( cdsFileMatching() ){

    fillHistRead_NC_ov(confMan, header, anaInfo, blMan, cdsMan, cdstrackMan);
  }
  return true;
}

void EventAnalysisReadNC_ov::InitializeHistogram()
{
  initHistRead_NC_ov();
}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysisReadNC_ov *event = new EventAnalysisReadNC_ov();
  return (EventTemp*)event;
}
