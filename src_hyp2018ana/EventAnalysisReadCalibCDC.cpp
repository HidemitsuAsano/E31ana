#include "EventAnalysisReadCalibCDC.h"

using namespace std;

bool EventAnalysisReadCalibCDC::UAna()
{
  if( !MyAnaTools::goodBeam(anaInfo) ) return true;

  if( cdsFileMatching() ){
    fillHistReadCalibCDC(confMan, blMan, cdsMan, cdstrackMan, anaInfo);
  }
  return true;
}

void EventAnalysisReadCalibCDC::InitializeHistogram()
{
  initHistReadCalibCDC(confMan);
}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysisReadCalibCDC *event = new EventAnalysisReadCalibCDC();
  return (EventTemp*)event;
}
