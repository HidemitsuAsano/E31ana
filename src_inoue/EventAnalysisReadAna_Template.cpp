#include "EventAnalysisReadAna.h"

using namespace std;

bool EventAnalysisReadAna::UAna()
{
  if( cdsFileMatching() ){

  }
  return true;
}

void EventAnalysisReadAna::InitializeHistogram()
{

}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysisReadAna *event = new EventAnalysisReadAna();
  return (EventTemp*)event;
}
