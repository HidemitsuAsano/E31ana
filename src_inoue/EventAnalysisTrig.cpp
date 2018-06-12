#include "File.h"
#include "ConfMan.h"
#include "TKO.h"
#include "CDSHitMan.h"
#include "BeamLineHitMan.h"
#include "EventHeader.h"
#include "ScalerMan.h"

#include "EventAlloc.h"
#include "EventTemp.h"

class EventAnalysisTrig: public EventTemp
{
public:
  EventAnalysisTrig();
  ~EventAnalysisTrig();
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

EventAnalysisTrig::EventAnalysisTrig()
  : EventTemp()
{
}

EventAnalysisTrig::~EventAnalysisTrig()
{
}

void EventAnalysisTrig::Initialize( ConfMan *conf )
{
  std::cout << " Enter EventAnalysisTrig::Initialize " << std::endl;
  confMan = conf;

  rtFile = new TFile( confMan->GetOutFileName().c_str(), "recreate" );
  scaTree= new TTree( "ScalerTree", "ScalerTree" );

  header = new EventHeader();
  cdsMan = new CDSHitMan();
  blMan = new BeamLineHitMan();

  scaMan = new ScalerMan();
  scaTree->Branch( "ScalerMan", &scaMan );
}

void EventAnalysisTrig::USca( int nsca, unsigned int *sca )
{
  Block_Event_Number++;
  header->SetBlockEventNumber( Block_Event_Number );
  scaMan->SetBlockEventNumber( Block_Event_Number );
  for( int i=0; i<nsca; i++ ){
    scaMan->AddHit( sca[i], confMan->GetScalerMapManager()->GetName(i) );
  }

  scaTree->Fill();
  scaMan->Clear();
}

bool EventAnalysisTrig::UAna( TKOHitCollection *tko )
{
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

  header->Clear();
  blMan->Clear();
  cdsMan->Clear();
  return true;
}

void EventAnalysisTrig::Finalize()
{
  std::cout << " Enter EventAnalysisTrig::Finalize " << std::endl;

  gFile->Write();
  gFile->Close();

  delete blMan;
  delete cdsMan;
  delete header;
}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysisTrig *event = new EventAnalysisTrig();
  return (EventTemp*)event;
}
