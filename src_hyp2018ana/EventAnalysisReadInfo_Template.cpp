#include "EventAnalysisReadInfo.h"

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
  if( header->IsTrig(Trig_Cosmic) || header->IsTrig(Trig_Reject) ) return Clear();
  blMan-> Convert(tko, confMan);
  cdsMan-> Convert(tko, confMan);

  // cout<<" Event# "<<Event_Number<<" BlockEvent# : "<<Block_Event_Number<<endl;
  // cout<<"> Trig Mode"<<header->trigmode()<<endl;

  if( !anaFileMatching() ){
    cout<<" !!! AnaInfo file matching fault !!!"<<endl;
  }
  if( anaInfo->beam(0) && anaInfo->minDCA() ){
    anaInfo->beam(0)->SetVertex(anaInfo->minDCA()->vertexBeam());
  }
  fillReadBeam_Kf(confMan, header, blMan, anaInfo->beam(0));

  BeamInfo *beam=anaInfo->beam(0);
  if( beam->tof()<27.8588 || 29.5663<beam->tof() ) beam->SetFlag(false);

  int nBPC=0;
  BLDCTrackInfo BPCInfo;
  for( int i=0; i<beam->nBPC(); i++ ){
    if( -30<beam->BPC(i).time() && beam->BPC(i).time()<100 ){
      BPCInfo=beam->BPC(i);
      nBPC++;
    }
  }
  if( nBPC!=1 || BPCInfo.time()<-10 || 10<BPCInfo.time()) beam->SetFlag(false);
  if( BPCInfo.chi2()>10 ) beam->SetFlag(false);

  int nBLC1=0;
  BLDCTrackInfo BLC1Info;
  for( int i=0; i<beam->nBLC1(); i++ ){
    if( -30<beam->BLC1(i).time() && beam->BLC1(i).time()<100 ){
      BLC1Info=beam->BLC1(i);
      nBLC1++;
    }
  }

  int nBLC2=0;
  BLDCTrackInfo BLC2Info;
  for( int i=0; i<beam->nBLC2(); i++ ){
    if( -30<beam->BLC2(i).time() && beam->BLC2(i).time()<100 ){
      BLC2Info=beam->BLC2(i);
      nBLC2++;
    }
  }

  if( nBLC1!=1 || BLC1Info.time()<-10 || 10<BLC1Info.time() || nBLC2!=1 || BLC2Info.time()<-10 || 10<BLC2Info.time()) beam->SetFlag(false);
  if( BLC1Info.chi2()>10 || BLC2Info.chi2()>10 ) beam->SetFlag(false);

  TVector3 BLC2BPC_diff=BLC2Info.GetPosatZ(0.5*(-130-20.3))-BPCInfo.GetPosatZ(0.5*(-130-20.3));
  TVector3 BLC2BPC_dir=BLC2Info.dir()-BPCInfo.dir();
  if( BLC2BPC_diff.X()<BLC2BPC_X_MIN || BLC2BPC_diff.X()>BLC2BPC_X_MAX ) beam->SetFlag(false);
  if( BLC2BPC_diff.Y()<BLC2BPC_Y_MIN || BLC2BPC_diff.Y()>BLC2BPC_Y_MAX ) beam->SetFlag(false);
  if( BLC2BPC_dir.X()<BLC2BPC_dX_MIN || BLC2BPC_dir.X()>BLC2BPC_dX_MAX ) beam->SetFlag(false);
  if( BLC2BPC_dir.Y()<BLC2BPC_dY_MIN || BLC2BPC_dir.Y()>BLC2BPC_dY_MAX ) beam->SetFlag(false);

  if( cdsFileMatching() ){
    fillReadCDS(anaInfo);
    fillCDS_pipi(anaInfo);
    fillCDS_ppim(anaInfo);
    fillReadNC(blMan, anaInfo);
    fillHistKN_pipi(header, blMan, cdsMan, cdstrackMan, anaInfo);
  }

  return Clear();
}

void EventAnalysisMakeInfo::InitializeHistogram()
{
  initHistReadBeam();
  initHistReadCDS();
  initHistCDS_pipi();
  initHistCDS_ppim();
  initHistReadNC();
  initHistKN_pipi();
}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysisMakeInfo *event = new EventAnalysisMakeInfo();
  return (EventTemp*)event;
}
