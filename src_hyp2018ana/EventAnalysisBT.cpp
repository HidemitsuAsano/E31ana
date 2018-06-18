#include "EventAnalysisBT.h"

using namespace std;

bool EventAnalysisBT::UAnaBeam()
{
  int nT0=0;
  HodoscopeLikeHit *T0hit=0;

  MyHistTools::fillT0(blMan);
  for( int i=0; i<blMan->nT0(); i++ ){
    if( blMan->T0(i)->CheckRange() ){
      nT0++;
      T0hit=blMan->T0(i);
    }
  }
  if( nT0!=1 ) return false;

  int trig_par = header->trigparticle();
  double expected_tof;
  if( trig_par==Beam_Kaon || trig_par==Beam_Other ) expected_tof=28.6429;
  if( trig_par==Beam_Pion   ) expected_tof=25.9334;
  if( trig_par==Beam_Proton ) expected_tof=35.22;
  double tof_diff=DBL_MAX;
  HodoscopeLikeHit *BHDhit=0;
  for( int i=0; i<blMan->nBHD(); i++ ){
    HodoscopeLikeHit *hit=blMan->BHD(i);
    if( hit-> CheckRange() ){
      double tof=T0hit->ctmean()-hit->ctmean();
      if( fabs(tof-expected_tof)<tof_diff ){
	BHDhit=hit;
	tof_diff=tof;
      }
    }
  }

  bltrackMan-> DoTracking(blMan, confMan, true, true);
  int BLC1id=-1, BLC2id=-1, BPCid=-1;
  int ntrackBLC1=0, ntrackBLC2=0, ntrackBPC=0;
  for( int i=0; i<bltrackMan->ntrackBLC1(); i++ ){
    LocalTrack *track = bltrackMan->trackBLC1(i);
    if( -30<track->GetTrackTime() && track->GetTrackTime()<100 ){
      BLC1id=i;
      ntrackBLC1++;
    }
  }

  for( int i=0; i<bltrackMan->ntrackBLC2(); i++ ){
    LocalTrack *track = bltrackMan->trackBLC2(i);
    if( -30<track->GetTrackTime() && track->GetTrackTime()<100 ){
      BLC2id=i;
      ntrackBLC2++;
    }
  }

  for( int i=0; i<bltrackMan->ntrackBPC(); i++ ){
    LocalTrack *track = bltrackMan->trackBPC(i);
    if( -30<track->GetTrackTime() && track->GetTrackTime()<100 ){
      BPCid=i;
      ntrackBPC++;
    }
  }
  //#if 0
  if( trig_par==Beam_Kaon ) cout<<"beam Kaon"<<endl;
  else if( trig_par==Beam_Pion ) cout<<"beam Pion"<<endl;
  else if( trig_par==Beam_Proton ) cout<<"beam Proton"<<endl;
  else cout<<"beam id : "<<trig_par<<endl;
  cout<<" ntrack BLC1 : "<<bltrackMan->ntrackBLC1()<<endl;
  cout<<" ntrack BLC2 : "<<bltrackMan->ntrackBLC2()<<endl;
  cout<<" ntrack BPC  : "<<bltrackMan->ntrackBPC()<<endl;
  cout<<" ntrack BLC1 : "<<ntrackBLC1<<endl;
  cout<<" ntrack BLC2 : "<<ntrackBLC2<<endl;
  cout<<" ntrack BPC  : "<<ntrackBPC<<endl;
  //#endif
  if( ntrackBLC1!=1 || ntrackBLC2!=1 || ntrackBPC!=1 ) return false;

  beamSpec-> TMinuitFit(bltrackMan->trackBLC1(BLC1id), bltrackMan->trackBLC2(BLC2id), confMan);

  return true;
}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysisBT *event = new EventAnalysisBT();
  return (EventTemp*)event;
}

