#include "EventAnalysisReadVtxDEF.h"

using namespace std;

bool EventAnalysisReadAna::UAna()
{
  if( !MyAnaTools::goodBeam(anaInfo) ) return true;

  if( cdsFileMatching() ){
    anaInfo-> ClearFC();
    ForwardChargeInfo fcInfo=MyTools::makeFC(confMan, blMan, anaInfo);
    //    fcInfo.dump();
    if( fcInfo.step()>1 ){
      anaInfo->AddCharge(fcInfo);
    }
    // for( int i=0; i<anaInfo->nCDS(); i++ ){
    //   CDSTrack *track=anaInfo->CDS(i)->track(confMan, cdsMan, cdstrackMan);
    //   // TVector3 vtxCDH=track->CDHVertex();
    //   // if( fabs(vtxCDH.Z())>40 ) anaInfo->CDS(i)->SetFlag(false);
    // }
    
    if( MyAnaTools::targetNA()==2 ){
      fillHistReadVtxDEF(header, confMan, anaInfo, blMan, bltrackMan, cdsMan, cdstrackMan);
    }
    else{
      //      fillHistReadKNpipi(header, blMan, cdsMan, cdstrackMan, anaInfo);
    }
  }
  return true;
}

void EventAnalysisReadAna::InitializeHistogram()
{
  TTree *tree=new TTree("npipi_event", "npipi_event");
  tree-> Branch("AnaInfo", &anaInfo);

  TTree *tree2=new TTree("ppimpim_event", "ppimpim_event");
  tree2-> Branch("AnaInfo", &anaInfo);

  if( MyAnaTools::targetNA()==2 ){
    initHistReadVtxDEF();
  }
  else{
    //    initHistReadKNpipi();
  }
}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysisReadAna *event = new EventAnalysisReadAna();
  return (EventTemp*)event;
}
