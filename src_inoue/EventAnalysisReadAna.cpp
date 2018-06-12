#include "EventAnalysisReadAna.h"

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
    
    fillHistReadCDS(anaInfo, cdsMan, cdstrackMan);
    fillHistReadNC(header, blMan, anaInfo);
    fillHistReadFC(header, confMan, blMan, anaInfo);
    fillHistReadKNpim(header, blMan, anaInfo);
    fillHistReadKNpip(header, blMan, anaInfo);
    fillHistReadKNkm(header, blMan, anaInfo);
    if( MyAnaTools::targetNA()==1 ){
      fillHistReadKNpipi_H2(header, confMan, blMan, cdsMan, cdstrackMan, anaInfo);
    }
    if( MyAnaTools::targetNA()==2 ){
      fillHistReadKNpipi_D2(header, blMan, cdsMan, cdstrackMan, anaInfo);
      fillHistReadKPpimpim_D2(header, confMan, blMan, cdsMan, cdstrackMan, anaInfo);
      fillHistReadCDS_Lpim_D2(header, blMan, anaInfo);

      //      fillHistReadVtxDEF(header, confMan, anaInfo, blMan, cdsMan, cdstrackMan);
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

  initHistReadCDS();
  initHistReadNC();
  initHistReadFC();
  initHistReadKNpim();
  initHistReadKNpip();
  initHistReadKNkm();
  if( MyAnaTools::targetNA()==1 ){
    initHistReadKNpipi_H2();
  }
  if( MyAnaTools::targetNA()==2 ){
    initHistReadKPpimpim_D2();
    initHistReadKNpipi_D2(header, anaInfo);
    initHistReadCDS_Lpim_D2();

    //    initHistReadVtxDEF();
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
