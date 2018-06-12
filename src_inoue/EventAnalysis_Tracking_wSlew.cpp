#include "File.h"
#include "ConfMan.h"
#include "TKO.h"
#include "CDSHitMan.h"
#include "CDSTrackingMan.h"
#include "CircleFit.h"
#include "HelixFit.h"
#include "BeamLineHitMan.h"
#include "BeamLineTrackMan.h"
#include "BeamSpectrometer.h"
#include "EventHeader.h"
#include "ScalerMan.h"

#include "EventAlloc.h"
#include "EventTemp.h"
#include "Display.h"

#include "TrackTools.h"
#include "Tools.h"
#include "MyTools.h"

using namespace std;

class EventAnalysis_Tracking_wSlew: public EventTemp
{
public:
  EventAnalysis_Tracking_wSlew();
  ~EventAnalysis_Tracking_wSlew();
private:
  TFile *rtFile;
  TTree *evTree;
  TTree *scaTree;
  BeamLineHitMan *blMan;
  BeamLineTrackMan *bltrackMan;
  BeamSpectrometer *beamSpec;
  CDSHitMan *cdsMan;
  CDSTrackingMan *trackMan;
  EventHeader *header;
  int t0,t1;
  TString tmpname;

  int AllGoodTrack;
  int nTrack;
  int nFill;

public:
  void Initialize( ConfMan *conf );
  void USca( int nsca, unsigned int *sca );
  bool UAna( TKOHitCollection *tko );
  bool Clear();
  void Finalize();
};

EventAnalysis_Tracking_wSlew::EventAnalysis_Tracking_wSlew()
  : EventTemp()
{
}

EventAnalysis_Tracking_wSlew::~EventAnalysis_Tracking_wSlew()
{
}

void EventAnalysis_Tracking_wSlew::Initialize( ConfMan *conf )
{
  std::cout << " Enter EventAnalysis_Tracking_wSlew::Initialize " << std::endl;

  confMan = conf;
  rtFile = new TFile( confMan->GetOutFileName().c_str() , "recreate" );
  rtFile->cd();

  evTree = new TTree( "EventTree", "EventTree" );
  header = new EventHeader();
  cdsMan = new CDSHitMan();
  trackMan = new CDSTrackingMan();  
  blMan = new BeamLineHitMan();
  bltrackMan = new BeamLineTrackMan();
  evTree->Branch( "BeamLineTrackMan", &bltrackMan );
  beamSpec = new BeamSpectrometer(conf);

  t0=time(0);
  AllGoodTrack=0;
  nTrack=0;
  nFill=0;

  evTree->Branch( "EventHeader", &header );
  evTree->Branch( "CDSTrackingMan", &trackMan );
}

void EventAnalysis_Tracking_wSlew::USca( int nsca, unsigned int *sca )
{
  Block_Event_Number++;
  header->SetBlockEventNumber( Block_Event_Number );
}

bool EventAnalysis_Tracking_wSlew::UAna( TKOHitCollection *tko )
{
  Event_Number++;
  { int status = confMan->CheckEvNum( Event_Number, Block_Event_Number );
    if( status==1 ) return true;
    if( status==2 ) return false; }

  if( Event_Number%100==0)
    {
      t1=time(0);
      std::cout << " Event# : " << Event_Number << "  BlockEvent# : " << Block_Event_Number << " AllTrack# : " << nTrack <<" GoodTrack# : " << AllGoodTrack << " Time (s): " << (t1-t0) << "  nFill : "<<nFill<<std::endl;
    }
  
  header->SetRunNumber(0);
  header->SetEventNumber(Event_Number);
  
  header->Convert( tko, confMan );

  if( header->IsTrig(Trig_Cosmic) ) return Clear();

  blMan-> Convert( tko, confMan );
  cdsMan->Convert( tko, confMan );
  std::vector<HodoscopeLikeHit*> T0hits=MyTools::getHodo(blMan, CID_T0);
  std::vector<HodoscopeLikeHit*> BHDhits=MyTools::getHodo(blMan, CID_BHD);
  if( T0hits.size()!=1 || BHDhits.empty() ) return Clear();

  if( MyTools::getCDH(cdsMan).empty() ) return Clear();
  bltrackMan-> DoTracking(blMan, confMan, true, true);

  int ntrackBPC=0, ntrackBLC1=0, ntrackBLC2=0;
  LocalTrack *BPC=0, *BLC1=0, *BLC2=0;
  for( int i=0; i<bltrackMan->ntrackBLC1(); i++ ){
    LocalTrack *track = bltrackMan->trackBLC1(i);
    if( -30<track->GetTrackTime() && track->GetTrackTime()<100 ){
      ntrackBLC1++;
      BLC1=track;
    }
  }
  for( int i=0; i<bltrackMan->ntrackBLC2(); i++ ){
    LocalTrack *track = bltrackMan->trackBLC2(i);
    if( -30<track->GetTrackTime() && track->GetTrackTime()<100 ){
      ntrackBLC2++;
      BLC2=track;
    }
  }
  for( int i=0; i<bltrackMan->ntrackBPC(); i++ ){
    LocalTrack *track = bltrackMan->trackBPC(i);
    if( -30<track->GetTrackTime() && track->GetTrackTime()<100 ){
      ntrackBPC++;
      BPC=track;
    }
  }
  if( ntrackBLC1!=1 || ntrackBLC2!=1 ) return Clear();
  if( ntrackBPC!=1 || ntrackBLC1!=1 || ntrackBLC2!=1 ) return Clear();

  beamSpec-> TMinuitFit(BLC1, BLC2, confMan);
  double beam_mom = beamSpec->mom();
  int beam_pid=Beam_Kaon;
  double beam_mass = kpMass;
  if( header->IsTrig(Trig_Pion) ){
    beam_pid=Beam_Pion;
    beam_mass=piMass;
  }

  trackMan->Execute(cdsMan,confMan);

  if( ntrackBPC==1 ){
    TVector3 T0pos=BPC->GetPosatZ(-110.5);
    for( int i=0; i<trackMan->nGoodTrack(); i++ ){
      CDSTrack *track = trackMan-> GoodTrack(i);
      int CDHseg; double CDHtime;
      if( MyTools::searchCDHHit(track, cdsMan) ){
	MyTools::searchIHHit(track, cdsMan);
	track->GetCDHHit(cdsMan, CDHseg, CDHtime);
	TVector3 vtxBeam, vtxCDS;
	if( !track-> GetVertex(T0pos, BPC-> GetMomDir(), vtxBeam, vtxCDS) ) continue;
	double beam_mom_vtx, beam_tof;
	ELossTools::CalcElossBeamTGeo(T0pos, vtxBeam, beam_mom, beam_mass, beam_mom_vtx, beam_tof);

	double mass2, beta, tof;
	if( !TrackTools::FindMass2(track, BPC, CDHtime-T0hits[0]->ctmean(), beam_mom, beam_pid, beta, mass2, tof) ) continue;

	track-> Retiming(cdsMan, confMan, beta, true);
      }
    }
  }

  int nGoodTrack=trackMan->nGoodTrack();
  int nallTrack=  trackMan->nTrack();

  AllGoodTrack+=nGoodTrack;
  nTrack+=nallTrack;

  if( nallTrack<1 ) return Clear();
  evTree->Fill();

  return Clear();
}

void EventAnalysis_Tracking_wSlew::Finalize()
{
  std::cout << " Enter EventAnalysis_Tracking_wSlew::Finalize    Event_Number : "<<Event_Number<< std::endl;

  rtFile->cd();
  confMan->SaveParams();
  gFile->Write();
  gFile->Close();

  delete cdsMan;
  delete trackMan;
  delete header;
}

bool EventAnalysis_Tracking_wSlew::Clear()
{
  header->Clear();
  blMan-> Clear();
  bltrackMan->Clear();
  cdsMan->Clear();
  trackMan->Clear();
  beamSpec->Clear();

  return true;
}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysis_Tracking_wSlew *event = new EventAnalysis_Tracking_wSlew();
  return (EventTemp*)event;
}
