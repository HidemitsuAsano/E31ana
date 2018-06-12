#include "File.h"
#include "ConfMan.h"
#include "TKO.h"
#include "CDSHitMan.h"
#include "BeamLineHitMan.h"
#include "EventHeader.h"
#include "ScalerMan.h"

#include "EventAlloc.h"
#include "EventTemp.h"
//#include "Display.h"
#include "DisplayBLC.h"
// This class is temporary
#include "BeamLineTrackMan.h"

// This class need to Draw TCanvas and so on.
#include "TRint.h"

class EventDispBLC2: public EventTemp
{
public:
  EventDispBLC2();
  ~EventDispBLC2();
private:
  TRint *theApp;
  TFile *rtFile;
  TTree *evTree;
  TTree *scaTree;
  CDSHitMan *cdsMan;
  BeamLineHitMan *blMan;
  BeamLineTrackMan *bltrackMan;
  EventHeader *header;
  ScalerMan *scaMan;
  DisplayBLC *dispBLC;

public:
  void Initialize( ConfMan *conf );
  void USca( int nsca, unsigned int *sca );
  bool UAna( TKOHitCollection *tko );
  void Finalize();
};

EventDispBLC2::EventDispBLC2()
  : EventTemp()
{
  int argc2=0;
  char **argv2;
  theApp = new TRint("theApp", &argc2, argv2);
}

EventDispBLC2::~EventDispBLC2()
{
  delete theApp;
}

const int MaxTreeSize = 1900000000000;
void EventDispBLC2::Initialize( ConfMan *conf )
{
#if 1
  std::cout << " Enter EventDispBLC2::Initialize " << std::endl;
#endif
  confMan = conf;
  
  rtFile = new TFile( confMan->GetOutFileName().c_str(), "recreate" );
  evTree = new TTree( "EventTree", "EventTree" );
  scaTree= new TTree( "ScalerTree", "ScalerTree" );

  header = new EventHeader();
  if( header==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  cdsMan = new CDSHitMan();
  if( cdsMan==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  blMan = new BeamLineHitMan();
  if( blMan ==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  scaMan = new ScalerMan();
  if( scaMan==NULL ){ std::cerr << "!!!!" << std::endl; return; }

  bltrackMan = new BeamLineTrackMan();
  dispBLC = new DisplayBLC();
}

void EventDispBLC2::USca( int nsca, unsigned int *sca )
{
#if 0
  std::cout << " Enter EventDispBLC2::USca " << std::endl;
#endif
  Block_Event_Number++;
  header->SetBlockEventNumber( Block_Event_Number );
  scaMan->SetBlockEventNumber( Block_Event_Number );
  for( int i=0; i<nsca; i++ ){
   scaMan->AddHit( sca[i], confMan->GetScalerMapManager()->GetName(i) );
  }
#if 0
  std::cout << nsca<< std::endl;
  for( int i=0; i<nsca; i++ ){
   std::cout << "  " << sca[i];
  }
  std::cout << std::endl;
#endif
  scaMan->Clear();
}

bool EventDispBLC2::UAna( TKOHitCollection *tko )
{
#if 0
  std::cout << " Enter EventDispBLC2::UAna " << std::endl;
#endif

  Event_Number++;
  { int status = confMan->CheckEvNum( Event_Number, Block_Event_Number );
    if( status==1 ) return true;
    if( status==2 ) return false; }

  if( Event_Number%5000==0 )
    std::cout << " Event# : " << Event_Number << "  BlockEvent# : " << Block_Event_Number << std::endl;
  
  
  header->SetRunNumber(0);
  header->SetEventNumber(Event_Number);
  
  header->Convert( tko, confMan );
  blMan->Convert( tko, confMan );
  cdsMan->Convert( tko, confMan );
  //  bltrackMan-> DoTracking( blMan, confMan, true, false);
  bltrackMan-> LocalTracking( blMan, confMan, CID_BLC2a, "notiming noslope");
  bltrackMan-> LocalTracking( blMan, confMan, CID_BLC2b, "notiming noslope");

  dispBLC-> Draw(confMan, blMan, bltrackMan);
  if( !dispBLC->Wait() ) return false;

  header->Clear();
  blMan->Clear();
  bltrackMan->Clear();
  cdsMan->Clear();
  return true;
}

void EventDispBLC2::Finalize()
{
  std::cout << " Enter EventDispBLC2::Finalize " << std::endl;

  gFile->Write();
  gFile->Close();

  delete blMan;
  delete bltrackMan;
  delete cdsMan;
  delete header;
}

EventTemp *EventAlloc::EventAllocator()
{
  EventDispBLC2 *event = new EventDispBLC2();
  return (EventTemp*)event;
}
