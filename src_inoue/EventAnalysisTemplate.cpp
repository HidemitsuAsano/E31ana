#include "EventAnalysisTemplate.h"

using namespace std;

bool EventAnalysisTemplate::UAna()
{
  TH1F *h1=NULL;
  h1=(TH1F*)gFile->Get("EventNumber"); h1->Fill(Event_Number);

  if( !UserTools::isGoodBeam(blMan, bltrackMan) ) return true;

  return true;
}

void EventAnalysisTemplate::InitializeHistogram()
{
  std::cout<<"EventAnalysisTemplate::InitializeHistogram START"<<std::endl;
  new TH1F("EventNumber", "", 1000, 0, 1.0e6);
  std::cout<<"EventAnalysisTemplate::InitializeHistogram FINISH"<<std::endl;
}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysisTemplate *event = new EventAnalysisTemplate();
  return (EventTemp*)event;
}
