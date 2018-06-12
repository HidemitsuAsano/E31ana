#include "EventAnalysisReadCalib.h"

using namespace std;

bool EventAnalysisReadCalib::UAna()
{
  if( !MyAnaTools::goodBeam(anaInfo) ) return true;

  if( cdsFileMatching() ){
    for( int i=0; i<anaInfo->nCDS(); i++ ){
      anaInfo->CDS(i)->track(confMan, cdsMan, cdstrackMan);
    }

    fillHistCalibCDS(header, confMan, cdsMan, cdstrackMan, anaInfo);
  }
  return true;
}

void EventAnalysisReadCalib::InitializeHistogram()
{
  initHistCalibCDS();
}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysisReadCalib *event = new EventAnalysisReadCalib();
  return (EventTemp*)event;
}
