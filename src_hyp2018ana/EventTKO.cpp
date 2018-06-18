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

class EventTKO: public EventTemp
{
public:
  EventTKO();
  ~EventTKO();
private:
  TFile *rtFile;
  TTree *evTree;
  TTree *scaTree;
  TKOHitCollection *tkoHitCol;
  //CDSHitMan *cdsMan;
  //BeamLineHitMan *blMan;
  EventHeader *header;
  ScalerMan *scaMan;
public:
  void Initialize( ConfMan *conf );
  void USca( int nsca, unsigned int *sca );
  bool UAna( TKOHitCollection *tko );
  void Finalize();
};

EventTKO::EventTKO()
  : EventTemp()
{
}

EventTKO::~EventTKO()
{
}

const long int MaxTreeSize = 1900000000;
void EventTKO::Initialize( ConfMan *conf )
{
#if 1
  std::cout << " Enter EventTKO::Initialize " << std::endl;
#endif
  confMan = conf;
  
  rtFile = new TFile( confMan->GetOutFileName().c_str(), "recreate" );
  evTree = new TTree( "EventTree", "EventTree" );
  scaTree= new TTree( "ScalerTree", "ScalerTree" );

  header = new EventHeader();
  if( header==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  evTree->Branch( "EventHeader", &header );
  tkoHitCol = new TKOHitCollection();
  if( tkoHitCol==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  evTree->Branch( "TKOHitCol", &tkoHitCol );
  /*
  cdsMan = new CDSHitMan();
  if( cdsMan==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  evTree->Branch( "CDSHitMan", &cdsMan );
  blMan = new BeamLineHitMan();
  if( blMan ==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  evTree->Branch( "BeamLineHitMan", &blMan );
  */
  scaMan = new ScalerMan();
  if( scaMan==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  scaTree->Branch( "ScalerMan", &scaMan );
  
}

void EventTKO::USca( int nsca, unsigned int *sca )
{
#if 0
  std::cout << " Enter EventTKO::USca " << std::endl;
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

bool EventTKO::UAna( TKOHitCollection *tko )
{
#if 1
  std::cout << " Enter EventTKO::UAna " << std::endl;
#endif

  Event_Number++;
  { int status = confMan->CheckEvNum( Event_Number, Block_Event_Number );
    if( status==1 ) return true;
    if( status==2 ) return false; }

  //  if( Event_Number%5000==0 )
    std::cout << " Event# : " << Event_Number << "  BlockEvent# : " << Block_Event_Number << std::endl;
  
  
  header->SetRunNumber(0);
  header->SetEventNumber(Event_Number);
  
  header->Convert( tko, confMan );
  for(int ih=0; ih<tko->entries(); ih++){  
    TKOHit *hithit = new TKOHit( tko->hit(ih)->cr(), tko->hit(ih)->sl(),tko->hit(ih)->ch(), tko->hit(ih)->data() );
  tkoHitCol->AddHit(*hithit);
  delete hithit ;
     std::cout <<"Add TKOHitCollection !!! Create:"<<tko->hit(ih)->cr()<<" Slot:"<<tko->hit(ih)->sl()<<" Channel:"<<tko->hit(ih)->ch()<<" Data:"<<tko->hit(ih)->data()<<std::endl;

  }

  evTree->Fill();
  header->Clear();
  tkoHitCol->Clear();

  return true;
}

void EventTKO::Finalize()
{
  std::cout << " Enter EventTKO::Finalize " << std::endl;

  gFile->Write();
  gFile->Close();

  //  delete blMan;
  //  delete cdsMan;
  delete header;
  delete tkoHitCol;
}

EventTemp *EventAlloc::EventAllocator()
{
  EventTKO *event = new EventTKO();
  return (EventTemp*)event;
}
