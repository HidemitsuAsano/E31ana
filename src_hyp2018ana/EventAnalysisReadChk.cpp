#include "EventAnalysisReadChk.h"

using namespace std;

bool EventAnalysisReadAna::UAna()
{

  fillHistTrigChk(header, confMan, blMan, cdsMan, anaInfo);

  if( !MyAnaTools::goodBeam(anaInfo) ) return true;

  if( cdsFileMatching() ){
    for( int i=0; i<anaInfo->nFCharge(); i++ ){
      if( anaInfo->forwardCharge(i)->pid()==F_Proton ){
	//	anaInfo->forwardCharge(i)->fit_forward(blMan, anaInfo->beam(0), anaInfo->minDCA(), confMan, true);
      }
    }

    // for( int i=0; i<anaInfo->nCDS(); i++ ){
    //   CDSTrack *track=anaInfo->CDS(i)->track(confMan, cdsMan, cdstrackMan);
    //   // TVector3 vtxCDH=track->CDHVertex();
    //   // if( fabs(vtxCDH.Z())>40 ) anaInfo->CDS(i)->SetFlag(false);
    // }

  }
  return true;
}

void EventAnalysisReadAna::InitializeHistogram()
{
  initHistTrigChk();
}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysisReadAna *event = new EventAnalysisReadAna();
  return (EventTemp*)event;
}
