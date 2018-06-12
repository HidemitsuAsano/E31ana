#include "EventAnalysisReadCDC_Template.h"

using namespace std;

bool EventAnalysisReadCDC_Template::UAna()
{
  TH1F *h1=NULL;
  h1=(TH1F*)gFile->Get("EventNumber"); h1->Fill(Event_Number);
  if( !UserTools::isGoodBeam(blMan, bltrackMan, beamSpec) ) return true;

  if( cdsFileMatching() ){
    h1=(TH1F*)gFile->Get("EventNumber_CDCmatched"); h1->Fill(Event_Number);
  }
  

  return true;
}

void EventAnalysisReadCDC_Template::InitializeHistogram()
{
  std::cout<<"EventAnalysisReadCDC_Template::InitializeHistogram START"<<std::endl;
  new TH1F("EventNumber", "", 1000, 0, 1.0e6);
  new TH1F("EventNumber_CDCmatched", "", 1000, 0, 1.0e6);
  std::cout<<"EventAnalysisReadCDC_Template::InitializeHistogram FINISH"<<std::endl;
}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysisReadCDC_Template *event = new EventAnalysisReadCDC_Template();
  return (EventTemp*)event;
}
