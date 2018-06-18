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

#include "TRint.h"
#include "DispBPC.h"
#include "BPCTrackMan.h"

class EventAnalysisBPC: public EventTemp
{
public:
  EventAnalysisBPC();
  ~EventAnalysisBPC();

private:
  TRint *theApp;

  TFile *rtFile;
  CDSHitMan *cdsMan;
  BeamLineHitMan *blMan;
  EventHeader *header;
  ScalerMan *scaMan;

  BeamLineTrackMan *bltrackMan;
  DispBPC *disp;

public:
  void Initialize( ConfMan *conf );
  void USca( int nsca, unsigned int *sca );
  bool UAna( TKOHitCollection *tko );
  void Finalize();
};

EventAnalysisBPC::EventAnalysisBPC()
  : EventTemp()
{
  int argc=0;
  char **argv;

  theApp = new TRint("theApp", &argc, argv);
}

EventAnalysisBPC::~EventAnalysisBPC()
{
  delete theApp;
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
  disp = new DispBPC(confMan);
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

  if( Event_Number%1==0 )
    std::cout << " Event# : " << Event_Number << "  BlockEvent# : " << Block_Event_Number << std::endl;
    
  header->SetRunNumber(0);
  header->SetEventNumber(Event_Number);
  
  header->Convert( tko, confMan );
  blMan->Convert( tko, confMan );
  cdsMan->Convert( tko, confMan );
  BPCTrackMan bpcTrack(confMan, blMan);
  bltrackMan-> LocalTracking(blMan, confMan, CID_BPC, "notiming noslope");
  bool flag = true;
  //  if( bltrackMan->ntrackBPC()>0 ){
  bpcTrack.DoTracking();
  bpcTrack.dumpHits();
  bpcTrack.dumpTrack();
  flag = disp-> draw(blMan, bltrackMan, bpcTrack);
    //  }

  header->Clear();
  blMan->Clear();
  cdsMan->Clear();
  bltrackMan-> Clear();

  return flag;
}

void EventAnalysisBPC::Finalize()
{
  gFile->Write();
  gFile->Close();

  delete blMan;
  delete cdsMan;
  delete header;

  delete bltrackMan;
  delete disp;
}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysisBPC *event = new EventAnalysisBPC();
  return (EventTemp*)event;
}
