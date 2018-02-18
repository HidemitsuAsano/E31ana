#include "File.h"
#include "ConfMan.h"
#include "TKO.h"
#include "CDSHitMan.h"
#include "BeamLineHitMan.h"
#include "BeamLineTrackMan.h"
#include "EventHeader.h"
#include "ScalerMan.h"
#include "Particle.h"

#include "EventAlloc.h"
#include "EventTemp.h"
#include "Display.h" 
#include "ELossTools.h"
#include "TrackTools.h"
#include "MyAnalysisBL.h"

#define Debug 0

class EventAnalysisMyBL: public EventTemp
{
  public:
    EventAnalysisMyBL();
    ~EventAnalysisMyBL();
  private:
    TFile *rtFile;
    TTree *evTree;
    TTree *scaTree;
    TFile *cdcFile;
    TTree *cdcTree;
    CDSHitMan *cdsMan;
    CDSTrackingMan *cdstrackMan;
    EventHeader *cdsheader;
    BeamLineHitMan *blMan;
    BeamLineTrackMan *bltrackMan;
    EventHeader *header;
    ScalerMan *scaMan;
    MyAnalysisBL* blAna;
    Particle* particle;

    int t0, t1;

  public:
    void Initialize( ConfMan *conf );
    void USca( int nsca, unsigned int *sca );
    bool UAna( TKOHitCollection *tko );
    void Finalize();
    void Clear();

};

  EventAnalysisMyBL::EventAnalysisMyBL()
: EventTemp()
{
}

EventAnalysisMyBL::~EventAnalysisMyBL()
{
}

const int MaxTreeSize = 1900000000000;
void EventAnalysisMyBL::Initialize( ConfMan *conf )
{
#if 1
  std::cout << " Enter EventAnalysisMyBL::Initialize " << std::endl;
#endif
  confMan = conf;

  rtFile = new TFile( confMan->GetOutFileName().c_str(), "recreate" );
  header = new EventHeader();
  if( header==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  cdsMan = new CDSHitMan();
  if( cdsMan==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  blMan = new BeamLineHitMan();
  if( blMan ==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  bltrackMan = new BeamLineTrackMan();
  if( bltrackMan ==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  
  scaMan = new ScalerMan();
  if( scaMan==NULL ){ std::cerr << "!!!!" << std::endl; return; }

  blAna = new MyAnalysisBL(rtFile, confMan);
  particle = new Particle();

  rtFile->cd();
  evTree = new TTree("EventTree","EventTree");
  evTree -> Branch("EventHeader",&header);
  evTree -> Branch("Particle",&particle);

  t0=time(0);
}

void EventAnalysisMyBL::USca( int nsca, unsigned int *sca )
{
  Block_Event_Number++;
  header->SetBlockEventNumber( Block_Event_Number );
  scaMan->SetBlockEventNumber( Block_Event_Number );
  for( int i=0; i<nsca; i++ ){
    scaMan->AddHit( sca[i], confMan->GetScalerMapManager()->GetName(i) );
  }

  scaMan->Clear();
}

bool EventAnalysisMyBL::UAna( TKOHitCollection *tko )
{
  Event_Number++;
  { int status = confMan->CheckEvNum( Event_Number, Block_Event_Number );
    if( status==1 ) return true;
    if( status==2 ) return false; }

    if( Event_Number%1000==0 )
    {
      t1=time(0);
      std::cout << " Event# : " << Event_Number << "  BlockEvent# : " << Block_Event_Number << " Time (s):" << (t1-t0) << std::endl;
    }

    header->SetRunNumber(0);
    header->SetEventNumber(Event_Number);

    header->Convert( tko, confMan );
    blMan->Convert( tko, confMan );
    cdsMan->Convert( tko, confMan );

    rtFile->cd();

    if(!blAna->DoAnalysis(confMan, header, blMan, bltrackMan, particle)){
      Clear();
      return true;
    }

    evTree->Fill();

    Clear();
    return true;
}

void EventAnalysisMyBL::Clear()
{
  particle->Clear();
  blAna->Clear();
  header->Clear();
  blMan->Clear();
  bltrackMan->Clear();
  cdsMan->Clear();
}

void EventAnalysisMyBL::Finalize()
{
  std::cout << " Enter EventAnalysisMyBL::Finalize " << std::endl;

  rtFile->cd();
  gFile->Write();
  gFile->Close();

  delete particle;
  delete blAna;
  delete blMan;
  delete bltrackMan;
  delete cdsMan;
  delete header;
}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysisMyBL *event = new EventAnalysisMyBL();
  return (EventTemp*)event;
}
