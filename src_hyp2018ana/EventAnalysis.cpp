#include "File.h"
#include "ConfMan.h"
#include "TKO.h"
#include "CDSHitMan.h"
#include "BeamLineHitMan.h"
#include "EventHeader.h"
#include "ScalerMan.h"

#include "EventAlloc.h"
#include "EventTemp.h"
#include "Display.h"

class EventAnalysis: public EventTemp
{
public:
  EventAnalysis();
  ~EventAnalysis();
private:
  TFile *rtFile;
  TTree *evTree;
  TTree *scaTree;
  CDSHitMan *cdsMan;
  BeamLineHitMan *blMan;
  EventHeader *header;
  ScalerMan *scaMan;
public:
  void Initialize( ConfMan *conf );
  void USca( int nsca, unsigned int *sca );
  bool UAna( TKOHitCollection *tko );
  void Finalize();
};

EventAnalysis::EventAnalysis()
  : EventTemp()
{
}

EventAnalysis::~EventAnalysis()
{
}

const int MaxTreeSize = 1900000000000;
void EventAnalysis::Initialize( ConfMan *conf )
{
#if 1
  std::cout << " Enter EventAnalysis::Initialize " << std::endl;
#endif
  confMan = conf;
  
  rtFile = new TFile( confMan->GetOutFileName().c_str(), "recreate" );
  evTree = new TTree( "EventTree", "EventTree" );
  scaTree= new TTree( "ScalerTree", "ScalerTree" );

  header = new EventHeader();
  if( header==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  evTree->Branch( "EventHeader", &header );
  cdsMan = new CDSHitMan();
  if( cdsMan==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  evTree->Branch( "CDSHitMan", &cdsMan );
  blMan = new BeamLineHitMan();
  if( blMan ==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  evTree->Branch( "BeamLineHitMan", &blMan );
  scaMan = new ScalerMan();
  if( scaMan==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  scaTree->Branch( "ScalerMan", &scaMan );
}

void EventAnalysis::USca( int nsca, unsigned int *sca )
{
#if 0
  std::cout << " Enter EventAnalysis::USca " << std::endl;
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
  scaTree->Fill();
  scaMan->Clear();
}

bool EventAnalysis::UAna( TKOHitCollection *tko )
{
#if 0
  std::cout << " Enter EventAnalysis::UAna " << std::endl;
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

  evTree->Fill();

  header->Clear();
  blMan->Clear();
  cdsMan->Clear();
  return true;
}

void EventAnalysis::Finalize()
{
  std::cout << " Enter EventAnalysis::Finalize " << std::endl;

  gFile->Write();
  gFile->Close();

  delete blMan;
  delete cdsMan;
  delete header;
}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysis *event = new EventAnalysis();
  return (EventTemp*)event;
}
