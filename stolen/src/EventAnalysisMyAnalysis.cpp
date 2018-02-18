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

#include "MyAnalysisBL.h"

#define Debug 0

class EventAnalysisMyAnalysis: public EventTemp
{
  public:
    EventAnalysisMyAnalysis();
    ~EventAnalysisMyAnalysis();
  private:
    TFile *rtFile;
    TTree *evTree;
    TTree *scaTree;
    CDSHitMan *cdsMan;
    BeamLineHitMan *blMan;
    BeamLineTrackMan *bltrackMan;
    EventHeader *header;
    ScalerMan *scaMan;
  public:
    void Initialize( ConfMan *conf );
    void USca( int nsca, unsigned int *sca );
    bool UAna( TKOHitCollection *tko );
    void Finalize();
};

  EventAnalysisMyAnalysis::EventAnalysisMyAnalysis()
: EventTemp()
{
}

EventAnalysisMyAnalysis::~EventAnalysisMyAnalysis()
{
}

const int MaxTreeSize = 1900000000000;
void EventAnalysisMyAnalysis::Initialize( ConfMan *conf )
{
#if 1
  std::cout << " Enter EventAnalysisMyAnalysis::Initialize " << std::endl;
#endif
  confMan = conf;

  rtFile = new TFile( confMan->GetOutFileName().c_str(), "recreate" );

  header = new EventHeader();
  if( header==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  cdsMan = new CDSHitMan();
  if( cdsMan==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  blMan = new BeamLineHitMan();
  if( blMan ==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  scaMan = new ScalerMan();
  if( scaMan==NULL ){ std::cerr << "!!!!" << std::endl; return; }
}

void EventAnalysisMyAnalysis::USca( int nsca, unsigned int *sca )
{
  Block_Event_Number++;
  header->SetBlockEventNumber( Block_Event_Number );
  scaMan->SetBlockEventNumber( Block_Event_Number );
  for( int i=0; i<nsca; i++ ){
    scaMan->AddHit( sca[i], confMan->GetScalerMapManager()->GetName(i) );
  }

  TH1F *h1;
  for( int i=0; i<scaMan->nsca(); i++ ){
    int val = scaMan->sca(i)->val();
    TString name = scaMan->sca(i)->name();
    h1 = (TH1F*)gFile->Get("Scaler"); h1->Fill(i,val);
  }


  scaMan->Clear();
}

bool EventAnalysisMyAnalysis::UAna( TKOHitCollection *tko )
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

    rtFile->cd();
    TH1F *h1;
    TH2F *h2;
    DetectorList *dlist=DetectorList::GetInstance();

    // ======== //
    // Raw Data //
    // ======== //


    // Multiplicity //
    int MulBHD=0;	
    for(int i=0; i<blMan->nBHD(); i++){
      HodoscopeLikeHit* hit = blMan->BHD(i);
      if(hit->CheckRange()) MulBHD++;
    }
    int MulT0=0;	
    for(int i=0; i<blMan->nT0(); i++){
      HodoscopeLikeHit* hit = blMan->T0(i);
      if(hit->CheckRange()) MulT0++;
    }
    int MulBPD=0;	
    for(int i=0; i<blMan->nBPD(); i++){
      HodoscopeLikeHit* hit = blMan->BPD(i);
      if(hit->CheckRange()) MulBPD++;
    }
    int MulDEF=0;	
    for(int i=0; i<blMan->nDEF(); i++){
      HodoscopeLikeHit* hit = blMan->DEF(i);
      if(hit->CheckRange()) MulDEF++;
    }

    // Selection //
    if(MulT0!=1||!header->IsTrig(Trig_Pion)){
      header->Clear();
      blMan->Clear();
      cdsMan->Clear();
      return true;
    }


    header->Clear();
    blMan->Clear();
    cdsMan->Clear();
    return true;
}

void EventAnalysisMyAnalysis::Finalize()
{
  std::cout << " Enter EventAnalysisMyAnalysis::Finalize " << std::endl;

  rtFile->cd();
  gFile->Write();
  gFile->Close();

  delete blMan;
  delete cdsMan;
  delete header;
}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysisMyAnalysis *event = new EventAnalysisMyAnalysis();
  return (EventTemp*)event;
}
