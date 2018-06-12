#include "File.h"
#include "ConfMan.h"
#include "TKO.h"
#include "CDSHitMan.h"
#include "BeamLineHitMan.h"
#include "EventHeader.h"
#include "ScalerMan.h"
#include "BeamLineTrackMan.h"
#include "CDSTrackingMan.h"

#include "EventAlloc.h"
#include "EventTemp.h"

#include "AnaMan.h"
#include "HistManCheck.h"

class EventAnalysisData: public EventTemp
{
public:
  EventAnalysisData();
  ~EventAnalysisData();

private:
  int Event_Number_onSpill;
  int CDC_Event_Num;
  TFile *cdcFile;
  TTree *cdcTree;
  CDSTrackingMan *cdstrackMan;
  EventHeader *cdcHeader;

  TFile *rtFile;
  CDSHitMan *cdsMan;
  CDSTrackingMan *trackingMan;
  BeamLineHitMan *blMan;
  EventHeader *header;
  ScalerMan *scaMan;

  BeamLineTrackMan *bltrackMan;

  AnaMan *anaMan;
  HistManCheck *histMan;

public:
  void Initialize( ConfMan *conf );
  void USca( int nsca, unsigned int *sca );
  bool UAna( TKOHitCollection *tko );
  void Clear();
  void Finalize();
};

EventAnalysisData::EventAnalysisData()
  : EventTemp()
{
  std::cout<<"EventAnalysisData::Constractor Call"<<std::endl;
}

EventAnalysisData::~EventAnalysisData()
{
}

const int MaxTreeSize = 1900000000000;
void EventAnalysisData::Initialize( ConfMan *conf )
{
  std::cout<<"EventAnalysisData::Initialization START"<<std::endl;
  Event_Number_onSpill = 0;
  CDC_Event_Num = 0;
  confMan = conf;
  cdcFile = new TFile((TString)confMan->GetCDSTrackFileName());
  cdstrackMan = new CDSTrackingMan();
  cdcHeader = new EventHeader();
  cdcTree = (TTree*)cdcFile-> Get("EventTree");
  cdcTree-> SetBranchAddress("EventHeader", &cdcHeader);
  cdcTree-> SetBranchAddress("CDSTrackingMan", &cdstrackMan);
  
  rtFile = new TFile( confMan->GetOutFileName().c_str(), "recreate" );

  header = new EventHeader();
  cdsMan = new CDSHitMan();
  blMan = new BeamLineHitMan();
  scaMan = new ScalerMan();

  bltrackMan = new BeamLineTrackMan();

  anaMan = new AnaMan(confMan);
  histMan = new HistManCheck();

  std::cout<<"EventAnalysisData::Initialization FINISH"<<std::endl;
}

void EventAnalysisData::USca( int nsca, unsigned int *sca )
{
  Block_Event_Number++;
  header->SetBlockEventNumber( Block_Event_Number );
  scaMan->SetBlockEventNumber( Block_Event_Number );
  for( int i=0; i<nsca; i++ ){
   scaMan->AddHit( sca[i], confMan->GetScalerMapManager()->GetName(i) );
  }

  scaMan->Clear();
}

bool EventAnalysisData::UAna( TKOHitCollection *tko )
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
  if( header->IsTrig(Trig_Cosmic) ){
    Clear();
    return true;
  }
  Event_Number_onSpill++;

  blMan->Convert( tko, confMan );
  cdsMan->Convert( tko, confMan );
  bltrackMan-> DoTracking(blMan, confMan, true, true);

  //*** CDC Tracking file event num check ***//
  //*** if not match UAna finish          ***//
  if( CDC_Event_Num>cdcTree->GetEntries() ){
    Clear();
    return true;
  }
  while( cdcHeader->ev()<Event_Number ){
    CDC_Event_Num++;
    if( CDC_Event_Num>cdcTree->GetEntries() ){
      Clear();
      return true;
    }
    cdcTree-> GetEntry(CDC_Event_Num);
  }
  if( cdcHeader->ev()>Event_Number ){
    Clear();
    return true;
  }

  if( cdcHeader->ev()!=Event_Number ){
    std::cout<<"  !!! CDC File Event matching miss !!!"<<std::endl;
    Clear();
    return false;
  }

  Clear();
  return true;
}

void EventAnalysisData::Clear()
{
  header->Clear();
  blMan->Clear();
  cdsMan->Clear();
  bltrackMan-> Clear();
}

void EventAnalysisData::Finalize()
{
  std::cout<<"EventAnalysisData Finish   Event Number : "<<Event_Number<<"  Event Number on Spill : "<<Event_Number_onSpill<<std::endl;
  gFile->Write();
  confMan-> SaveParams();
  confMan-> SaveCode();
  gFile->Close();

  delete blMan;
  delete cdsMan;
  delete header;

  delete bltrackMan;
}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysisData *event = new EventAnalysisData();
  return (EventTemp*)event;
}
