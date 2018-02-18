#include "File.h"
#include "ConfMan.h"
#include "TKO.h"
#include "CDSHitMan.h"
#include "BeamLineHitMan.h"
#include "BeamLineTrackMan.h"
#include "EventHeader.h"
#include "ScalerMan.h"

#include "EventAlloc.h"
#include "EventTemp.h"
#include "Display.h" 
#include "ELossTools.h"
#include "TrackTools.h"
#include "MyAnalysisBLDCCheckStop.h"
#include "MyAnalysisSDCCheck.h"
#include "MyAnalysisBLHodoCheckStop.h"

#define Debug 0

//#define BLHodoCheck
//#define BLDCCheck
#define SDCCheck

class EventAnalysisMyPreAna: public EventTemp
{
  public:
    EventAnalysisMyPreAna();
    ~EventAnalysisMyPreAna();
  private:
    TFile *rtFile;
    TTree *evTree;
    TTree *scaTree;

    BeamLineHitMan *blMan;
    BeamLineTrackMan *bltrackMan;
    EventHeader *header;
    ScalerMan *scaMan;

#ifdef BLDCCheck
    MyAnalysisBLDCCheck* bldccheckAna;
#endif
#ifdef SDCCheck
    MyAnalysisSDCCheck* bldccheckAna;
#endif
#ifdef BLHodoCheck
    MyAnalysisBLHodoCheck* hodocheckAna;
#endif

    int t0, t1;

  public:
    void Initialize( ConfMan *conf );
    void USca( int nsca, unsigned int *sca );
    bool UAna( TKOHitCollection *tko );
    void Finalize();
    void Clear();

};

  EventAnalysisMyPreAna::EventAnalysisMyPreAna()
: EventTemp()
{
}

EventAnalysisMyPreAna::~EventAnalysisMyPreAna()
{
}

const int MaxTreeSize = 1900000000000;
void EventAnalysisMyPreAna::Initialize( ConfMan *conf )
{
#if 1
  std::cout << " Enter EventAnalysisMyPreAna::Initialize " << std::endl;
#endif
  confMan = conf;

  rtFile = new TFile( confMan->GetOutFileName().c_str(), "recreate" );
  header = new EventHeader();
  if( header==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  blMan = new BeamLineHitMan();
  if( blMan ==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  bltrackMan = new BeamLineTrackMan();
  if( bltrackMan ==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  scaMan = new ScalerMan();
  if( scaMan==NULL ){ std::cerr << "!!!!" << std::endl; return; }

#ifdef BLDCCheck
  bldccheckAna  = new MyAnalysisBLDCCheck(rtFile, confMan);
#endif
#ifdef SDCCheck
  bldccheckAna  = new MyAnalysisSDCCheck(rtFile, confMan);
#endif
#ifdef BLHodoCheck
  hodocheckAna  = new MyAnalysisBLHodoCheck(rtFile, confMan);
#endif

  t0=time(0);
}

void EventAnalysisMyPreAna::USca( int nsca, unsigned int *sca )
{
  Block_Event_Number++;
  header->SetBlockEventNumber( Block_Event_Number );
  scaMan->SetBlockEventNumber( Block_Event_Number );
  for( int i=0; i<nsca; i++ ){
    scaMan->AddHit( sca[i], confMan->GetScalerMapManager()->GetName(i) );
  }

  scaMan->Clear();
}

bool EventAnalysisMyPreAna::UAna( TKOHitCollection *tko )
{
  Event_Number++;
  { int status = confMan->CheckEvNum( Event_Number, Block_Event_Number );
    if( status==1 ){
      Clear();
      return true;
    }
    if( status==2 ){
      Clear();
      return false;
    }
  }

  if( Event_Number%10000==0 )
  {
      t1=time(0);
      std::cout << " Event# : " << Event_Number << "  BlockEvent# : " << Block_Event_Number << " Time (s):" << (t1-t0) << std::endl;
  }

  header->SetRunNumber(0);
  header->SetEventNumber(Event_Number);

  header->Convert( tko, confMan );
  blMan->Convert( tko, confMan );

  rtFile->cd();

#ifdef BLDCCheck
  if(!bldccheckAna->DoAnalysis(confMan, header, blMan, bltrackMan)){
    //std::cout << " k0Ana::false " << std::endl;
    //Clear();
    //return true;
  }
#endif
#ifdef SDCCheck
  if(!bldccheckAna->DoAnalysis(confMan, header, blMan, bltrackMan)){
    //std::cout << " k0Ana::false " << std::endl;
    //Clear();
    //return true;
  }
#endif
#ifdef BLHodoCheck
  if(!hodocheckAna->DoAnalysis(confMan, header, blMan, bltrackMan)){
    //std::cout << " k0Ana::false " << std::endl;
    //Clear();
    //return true;
  }
#endif

  Clear();
  return true;
}

void EventAnalysisMyPreAna::Clear()
{
#ifdef BLDCCheck
  bldccheckAna->Clear();
#endif
#ifdef SDCCheck
  bldccheckAna->Clear();
#endif
#ifdef BLHodoCheck
  hodocheckAna->Clear();
#endif
  blMan->Clear();
  bltrackMan->Clear();
  header->Clear();
}

void EventAnalysisMyPreAna::Finalize()
{
  std::cout << " Enter EventAnalysisMyPreAna::Finalize " << std::endl;

  rtFile->cd();
  gFile->Write();
  gFile->Close();

#ifdef BLDCCheck
  delete bldccheckAna;
#endif
#ifdef SDCCheck
  delete bldccheckAna;
#endif
#ifdef BLHodoCheck
  delete hodocheckAna;
#endif
}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysisMyPreAna *event = new EventAnalysisMyPreAna();
  return (EventTemp*)event;
}
