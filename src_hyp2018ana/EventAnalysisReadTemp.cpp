#include "EventAnalysisReadTemp.h"

using namespace std;

bool EventAnalysisReadTemp::UAna()
{
  if( !MyAnaTools::goodBeam(anaInfo) ) return true;
  // Please fill your histograms

  if( cdsFileMatching() ){

  }
  return true;
}

void EventAnalysisReadTemp::InitializeHistogram()
{
  // Please write your histogram initialization
}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysisReadTemp *event = new EventAnalysisReadTemp();
  return (EventTemp*)event;
}
