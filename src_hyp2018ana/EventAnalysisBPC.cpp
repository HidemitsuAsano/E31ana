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

  int Trig18Event;
  int TrackingEvent;

  int Trig18Event2;
  int TrackingEvent2;

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

  Trig18Event = 0;
  TrackingEvent = 0;

  Trig18Event2 = 0;
  TrackingEvent2 = 0;
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

  //  bltrackMan-> DoTracking(blMan, confMan, false, false);
  bltrackMan-> LocalTracking(blMan, confMan, CID_BPC, "notiming noslope");
  effMan-> get();
  if( effMan-> trig18() ){
    Trig18Event++;
    if( bltrackMan->ntrackBPC()>0 ){
      TrackingEvent++;
    }
  }

  if( effMan-> trig18_2() ){
    Trig18Event2++;
    if( bltrackMan->ntrackBPC()>0 ){
      TrackingEvent2++;
    }
  }

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
  std::cout<<"EventAnalysisBPC::Finalize"<<std::endl;
  std::cout<<" 1 & 8 layer hit event : "<<Trig18Event<<std::endl;
  std::cout<<" Tracking Event : "<<TrackingEvent<<std::endl;
  std::cout<<"  Efficiency = "<<100.*TrackingEvent/Trig18Event<<"[%]"<<std::endl; 
  std::cout<<"w/o 2 wire"<<std::endl;
  std::cout<<" 1 & 8 layer hit event : "<<Trig18Event2<<std::endl;
  std::cout<<" Tracking Event : "<<TrackingEvent2<<std::endl;
  std::cout<<"  Efficiency = "<<100.*TrackingEvent2/Trig18Event2<<"[%]"<<std::endl; 
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
