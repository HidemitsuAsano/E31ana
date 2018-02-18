#include "File.h"
#include "ConfMan.h"
#include "TKO.h"
#include "CDSHitMan.h"
#include "CDSTrackingMan.h"
#include "CircleFit.h"
#include "HelixFit.h"
#include "BeamLineHitMan.h"
#include "BeamLineTrackMan.h"
#include "EventHeader.h"
#include "ScalerMan.h"

#include "EventAlloc.h"
#include "EventTemp.h"
#include "Display.h"

#define CDCKilledLayer6 1

class EventAnalysisMyCDCTracking: public EventTemp
{
public:
  EventAnalysisMyCDCTracking();
  ~EventAnalysisMyCDCTracking();
private:
  TFile *rtFile;
  TTree *evTree;
  TTree *scaTree;
  BeamLineHitMan *blMan;
  CDSHitMan *cdsMan;
  CDSTrackingMan *trackMan;
  EventHeader *header;
  int t0,t1;
  TString tmpname;

  int AllGoodTrack;
  int nTrack;

public:
  void Initialize( ConfMan *conf );
  void USca( int nsca, unsigned int *sca );
  bool UAna( TKOHitCollection *tko );
  void Finalize();
};

EventAnalysisMyCDCTracking::EventAnalysisMyCDCTracking()
  : EventTemp()
{
}

EventAnalysisMyCDCTracking::~EventAnalysisMyCDCTracking()
{
}

void EventAnalysisMyCDCTracking::Initialize( ConfMan *conf )
{
#if 1
  std::cout << " Enter EventAnalysisMyCDCTracking::Initialize " << std::endl;
#endif
  confMan = conf;
  rtFile = new TFile( confMan->GetOutFileName().c_str() , "recreate" );
  rtFile->cd();

  evTree = new TTree( "EventTree", "EventTree" );
  header = new EventHeader();
  if( header==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  evTree->Branch( "EventHeader", &header );
  blMan = new BeamLineHitMan();
  if( blMan ==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  cdsMan = new CDSHitMan();
  if( cdsMan==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  trackMan = new CDSTrackingMan();  
  if( trackMan==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  evTree->Branch( "CDSTrackingMan", &trackMan );

  t0=time(0);
  AllGoodTrack=0;
  nTrack=0;
}

void EventAnalysisMyCDCTracking::USca( int nsca, unsigned int *sca )
{
#if 0
  std::cout << " Enter EventAnalysisMyCDCTracking::USca " << std::endl;
#endif
  Block_Event_Number++;
  header->SetBlockEventNumber( Block_Event_Number );
}

bool EventAnalysisMyCDCTracking::UAna( TKOHitCollection *tko )
{
#if 0
  std::cout << " Enter EventAnalysisMyCDCTracking::UAna " << std::endl;
#endif

  Event_Number++;
  { int status = confMan->CheckEvNum( Event_Number, Block_Event_Number );
    if( status==1 ) return true;
    if( status==2 ) return false; }

    if( Event_Number%1000==0)
    {
      t1=time(0);
      std::cout << " Event# : " << Event_Number << "  BlockEvent# : " << Block_Event_Number << " AllTrack# : " << nTrack <<" GoodTrack# : " << AllGoodTrack << " Time (s): " << (t1-t0) << std::endl;
    }

    header->SetRunNumber(0);
    header->SetEventNumber(Event_Number);

    header->Convert( tko, confMan );
    blMan->Convert( tko, confMan );
    cdsMan->Convert( tko, confMan );

    // Selection //
    if((header->IsTrig(Trig_Cosmic))){
      header->Clear();
      blMan->Clear();
      cdsMan->Clear();
      trackMan->Clear();
      return true;
    }

#ifdef CDCKilledLayer1
    trackMan->SetKilledLayer(1);
#endif
#ifdef CDCKilledLayer2
    trackMan->SetKilledLayer(2);
#endif
#ifdef CDCKilledLayer3
    trackMan->SetKilledLayer(3);
#endif
#ifdef CDCKilledLayer4
    trackMan->SetKilledLayer(4);
#endif
#ifdef CDCKilledLayer5
    trackMan->SetKilledLayer(5);
#endif
#ifdef CDCKilledLayer6
    trackMan->SetKilledLayer(6);
#endif
#ifdef CDCKilledLayer7
    trackMan->SetKilledLayer(7);
#endif
#ifdef CDCKilledLayer8
    trackMan->SetKilledLayer(8);
#endif
#ifdef CDCKilledLayer9
    trackMan->SetKilledLayer(9);
#endif
#ifdef CDCKilledLayer10
    trackMan->SetKilledLayer(10);
#endif
#ifdef CDCKilledLayer11
    trackMan->SetKilledLayer(11);
#endif
#ifdef CDCKilledLayer12
    trackMan->SetKilledLayer(12);
#endif
#ifdef CDCKilledLayer13
    trackMan->SetKilledLayer(13);
#endif
#ifdef CDCKilledLayer14
    trackMan->SetKilledLayer(14);
#endif
#ifdef CDCKilledLayer15
    trackMan->SetKilledLayer(15);
#endif
    trackMan->Execute(cdsMan,confMan);

    int nGoodTrack=trackMan->nGoodTrack();
    int nallTrack=  trackMan->nTrack();
    AllGoodTrack+=nGoodTrack;
    nTrack+=nallTrack;
    if( nallTrack<1 ){
      header->Clear();
      blMan->Clear();
      cdsMan->Clear();
      trackMan->Clear();
      return true;
    }

    evTree->Fill();
    cdsMan->Clear();
    blMan->Clear();
    header->Clear();
    trackMan->Clear();
    return true;
}

void EventAnalysisMyCDCTracking::Finalize()
{
  std::cout << " Enter EventAnalysisMyCDCTracking::Finalize " << std::endl;

  rtFile->cd();
  confMan->SaveCDSParam();
  gFile->Write();
  gFile->Close();
  //  gSystem->Sleep(10000);
  //  gSystem->Exec(Form("mv %s %s",tmpname.Data(),confMan->GetOutFileName().c_str()));
  delete blMan;
  delete cdsMan;
  delete trackMan;
  delete header;
}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysisMyCDCTracking *event = new EventAnalysisMyCDCTracking();
  return (EventTemp*)event;
}
