#include "File.h"
#include "ConfMan.h"
#include "TKO.h"
#include "CDSHitMan.h"
#include "BeamLineHitMan.h"
#include "EventHeader.h"
#include "ScalerMan.h"
#include "BeamLineTrackMan.h"

#include "EventAlloc.h"
#include "EventTemp.h"
#include "Display.h"

#include "HistManBPC.h"
#include "DCEffMan.h"

class EventAnalysisBPC: public EventTemp
{
public:
  EventAnalysisBPC();
  ~EventAnalysisBPC();

private:
  TFile *rtFile;
  CDSHitMan *cdsMan;
  BeamLineHitMan *blMan;
  EventHeader *header;
  ScalerMan *scaMan;

  BeamLineTrackMan *bltrackMan;

  HistManBPC *histMan;
  DCEffMan *effMan;

public:
  void Initialize( ConfMan *conf );
  void USca( int nsca, unsigned int *sca );
  bool UAna( TKOHitCollection *tko );
  void Finalize();
};

EventAnalysisBPC::EventAnalysisBPC()
  : EventTemp()
{
}

EventAnalysisBPC::~EventAnalysisBPC()
{
}

const int MaxTreeSize = 1900000000000;
void EventAnalysisBPC::Initialize( ConfMan *conf )
{
  confMan = conf;
  
  rtFile = new TFile( confMan->GetOutFileName().c_str(), "recreate" );

  header = new EventHeader();
  cdsMan = new CDSHitMan();
  blMan = new BeamLineHitMan();
  scaMan = new ScalerMan();

  bltrackMan = new BeamLineTrackMan();

  histMan = new HistManBPC(rtFile);
  effMan = new DCEffMan(CID_BPC, confMan, blMan, bltrackMan);
}

void EventAnalysisBPC::USca( int nsca, unsigned int *sca )
{
  Block_Event_Number++;
  header->SetBlockEventNumber( Block_Event_Number );
  scaMan->SetBlockEventNumber( Block_Event_Number );
  for( int i=0; i<nsca; i++ ){
   scaMan->AddHit( sca[i], confMan->GetScalerMapManager()->GetName(i) );
  }

  scaMan->Clear();
}

bool EventAnalysisBPC::UAna( TKOHitCollection *tko )
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

  bltrackMan-> DoTracking(blMan, confMan);
  effMan-> get();

  histMan-> fill(confMan, blMan, bltrackMan, effMan);

#if 0
  effMan-> dump();
  char in;
  if( in=='q' ){
    return false;
  }
#endif

  header->Clear();
  blMan->Clear();
  cdsMan->Clear();
  bltrackMan-> Clear();
  return true;
}

void EventAnalysisBPC::Finalize()
{
  gFile->Write();
  gFile->Close();

  delete blMan;
  delete cdsMan;
  delete header;
}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysisBPC *event = new EventAnalysisBPC();
  return (EventTemp*)event;
}
