#include "EventAnalysisReadAna.h"

using namespace std;

bool EventAnalysisReadAna::UAna()
{
  if( !MyAnaTools::goodBeam(anaInfo) ) return true;

  if( cdsFileMatching() ){
    for( int i=0; i<anaInfo->nCDS(); i++ ){
      CDSTrack *track=anaInfo->CDS(i)->track(cdstrackMan);
      TVector3 vtxCDH=track->CDHVertex();
      if( fabs(vtxCDH.Z())>40 ) anaInfo->CDS(i)->SetFlag(false);
    }

    fillHistReadCDS(anaInfo, cdsMan, cdstrackMan);
    fillHistReadNC(header, blMan, anaInfo);
    fillHistReadFC(header, blMan, anaInfo);
    fillHistReadKNpim(header, blMan, anaInfo);
    fillHistReadKNpip(header, blMan, anaInfo);
    fillHistReadKNkm(header, blMan, anaInfo);
    fillHistReadKNpipi(header, blMan, cdsMan, cdstrackMan, anaInfo);
    fillHistReadKPpimpim(header, blMan, cdsMan, cdstrackMan, anaInfo);
    fillHistReadCDS_Lpim(header, blMan, anaInfo);
  }
  return true;
}

void EventAnalysisReadAna::InitializeHistogram()
{
  initHistReadCDS();
  initHistReadNC();
  initHistReadFC();
  initHistReadKNpim();
  initHistReadKNpip();
  initHistReadKNkm();
  initHistReadKNpipi();
  initHistReadKPpimpim();
  initHistReadCDS_Lpim();
}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysisReadAna *event = new EventAnalysisReadAna();
  return (EventTemp*)event;
}
