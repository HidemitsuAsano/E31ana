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

#include "MyAnalysisBase.h"
#include "MyAnalysisAlloc.h"

#define Debug 0

class EventAnalysisMyAna: public EventTemp
{
public:
  EventAnalysisMyAna();
  ~EventAnalysisMyAna();
private:
  CDSHitMan *cdsMan;
  BeamLineHitMan *blMan;
  EventHeader *header;
  ScalerMan *scaMan;
  MyAnalysisBase* myana;
  MyAnalysisAlloc* myalloc;
public:
  void Initialize( ConfMan *conf );
  void USca( int nsca, unsigned int *sca );
  bool UAna( TKOHitCollection *tko );
  void Finalize();
};

EventAnalysisMyAna::EventAnalysisMyAna()
  : EventTemp()
{
}

EventAnalysisMyAna::~EventAnalysisMyAna()
{
}

const int MaxTreeSize = 1900000000000;
void EventAnalysisMyAna::Initialize( ConfMan *conf )
{
#if 1
  std::cout << " Enter EventAnalysisMyAna::Initialize " << std::endl;
#endif
  confMan = conf;
  
  header = new EventHeader();
  if( header==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  cdsMan = new CDSHitMan();
  if( cdsMan==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  blMan = new BeamLineHitMan();
  if( blMan ==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  scaMan = new ScalerMan();
  if( scaMan==NULL ){ std::cerr << "!!!!" << std::endl; return; }

  myalloc = new MyAnalysisAlloc();
  myana = myalloc->MyAnalysisAllocator();
  myana -> Initialize(confMan);

}

void EventAnalysisMyAna::USca( int nsca, unsigned int *sca )
{
#if 0
  std::cout << " Enter EventAnalysisMyAna::USca " << std::endl;
#endif
  Block_Event_Number++;
  header->SetBlockEventNumber( Block_Event_Number );
  scaMan->SetBlockEventNumber( Block_Event_Number );
  for( int i=0; i<nsca; i++ ){
    scaMan->AddHit( sca[i], confMan->GetScalerMapManager()->GetName(i) );
  }

  myana -> DoScaler(scaMan);

  scaMan->Clear();
}

bool EventAnalysisMyAna::UAna( TKOHitCollection *tko )
{
#if 0
	std::cout << " Enter EventAnalysisMyAna::UAna " << std::endl;
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

    if(header->IsTrig(Trig_Cosmic))
      myana -> DoAnalysis(cdsMan,blMan);

    header->Clear();
    blMan->Clear();
    cdsMan->Clear();
    return true;
}

void EventAnalysisMyAna::Finalize()
{
  std::cout << " Enter EventAnalysisMyAna::Finalize " << std::endl;

  myana -> Finalize();  

  delete blMan;
  delete cdsMan;
  delete header;
  delete myalloc;
}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysisMyAna *event = new EventAnalysisMyAna();
  return (EventTemp*)event;
}
