#include "EventAnalysisMakeInfo.h"

using namespace std;

bool EventAnalysisMakeInfo::UAna(TKOHitCollection *tko)
{
  Event_Number++;
  if( Event_Number%1000==0 ){
    t1 = time(0);
    std::cout << " Event# : " << Event_Number << "  BlockEvent# : " << Block_Event_Number << " Time (s): "<<(t1-t0)<<std::endl;
  }

  int status = confMan->CheckEvNum( Event_Number, Block_Event_Number );
  if( status==1 ){
    return Clear();
  }
  if( status==2 ){
    return Clear(false);
  }

  header-> SetRunNumber(confMan->GetRunNumber());
  header-> SetEventNumber(Event_Number);
  header->Convert( tko, confMan );
  if( header->IsTrig(Trig_Cosmic) ){
    Num_Cosmic++;
    return Clear();
  }
  if( header->IsTrig(Trig_Reject) ){
    Num_Reject++;
    return Clear();
  }

  blMan-> Convert(tko, confMan);
  cdsMan-> Convert(tko, confMan);
  bltrackMan-> DoTracking(blMan, confMan, true, true);
  LocalTrack *trackBLC1=MyTools::trackBLC1(bltrackMan);
  LocalTrack *trackBLC2=MyTools::trackBLC2(bltrackMan);
  if( trackBLC1 && trackBLC2 ) beamSpec-> TMinuitFit(trackBLC1, trackBLC2, confMan);

  BeamInfo beam=MyTools::makeBeamInfo(confMan, header, blMan, bltrackMan, beamSpec);
  anaInfo->AddBeam(beam);

  if( cdsFileMatching() && anaInfo->beam(0)->flag() ){
    for( int i=0; i<cdstrackMan->nTrack(); i++ ){
      CDSInfo cds=MyTools::makeCDSInfo(i, cdsMan, cdstrackMan, anaInfo->beam(0), bltrackMan, confMan);
      anaInfo-> AddCDS(cds);
    }
    if( anaInfo->minDCA() ){ anaInfo->beam(0)-> SetVertex(anaInfo->minDCA()->vertexBeam()); }

    for( int i=0; i<cdstrackMan->nTrack(); i++ ){
      for( int j=i+1; j<cdstrackMan->nTrack(); j++ ){
	CDS2Info cds2=MyTools::makeCDS2Info(cdstrackMan, i, j, anaInfo);
	anaInfo-> AddCDS2(cds2);
      }
    }

    ForwardNeutralInfo fnInfo=MyTools::makeFN(blMan, anaInfo, NC_THRE, false);
    if( fnInfo.pid()==F_Gamma || fnInfo.pid()==F_Neutron ) anaInfo->AddNeutral(fnInfo);

    ForwardChargeInfo fcInfo=MyTools::makeFC(confMan, blMan, anaInfo);
    if( fcInfo.step()>1 ){
      anaInfo->AddCharge(fcInfo);
    }
  }

  /* std::cout<<" Event Number : "<<Event_Number<<std::endl; */
  /* anaInfo-> dump(); */
  /* std::string str; std::cin>>str; if( str=="q" ) return false; */
  evTree-> Fill();
  return Clear();
}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysisMakeInfo *event = new EventAnalysisMakeInfo();
  return (EventTemp*)event;
}
