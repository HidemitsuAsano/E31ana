#include "File.h"
#include "ConfMan.h"
#include "TKO.h"
#include "CDSHitMan.h"
#include "CDSTrackingMan.h"
#include "CircleFit.h"
#include "HelixFit.h"
#include "BeamLineHitMan.h"
#include "BeamLineTrackMan.h"
#include "EventHeader.h"
#include "ScalerMan.h"

#include "EventAlloc.h"
#include "EventTemp.h"
#include "Display.h"

class EventAnalysis: public EventTemp
{
public:
  EventAnalysis();
  ~EventAnalysis();
private:
  TFile *rtFile;
  TTree *evTree;
  TTree *scaTree;
  CDSHitMan *cdsMan;
  CDSTrackingMan *trackMan;
  EventHeader *header;
  int t0,t1;
  TString tmpname;

  int AllGoodTrack;
  int nTrack;

public:
  void Initialize( ConfMan *conf );
  void USca( int nsca, unsigned int *sca );
  bool UAna( TKOHitCollection *tko );
  void Finalize();
};

EventAnalysis::EventAnalysis()
  : EventTemp()
{
}

EventAnalysis::~EventAnalysis()
{
}

void EventAnalysis::Initialize( ConfMan *conf )
{
#if 1
  std::cout << " Enter EventAnalysis::Initialize " << std::endl;
#endif
  confMan = conf;
  // int i=0;
  // tmpname=Form("tmp/tmptra%d.root",i);
  // while( !gSystem->Exec( Form("test -f %s",tmpname.Data()) ))
  //   tmpname=Form("tmp/tmptra%d.root",++i);
  // std::cout<<tmpname<<std::endl;
  rtFile = new TFile( confMan->GetOutFileName().c_str() , "recreate" );
  rtFile->cd();

  evTree = new TTree( "EventTree", "EventTree" );
  header = new EventHeader();
  if( header==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  evTree->Branch( "EventHeader", &header );
  cdsMan = new CDSHitMan();
  if( cdsMan==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  //  evTree->Branch( "CDSHitMan", &cdsMan );
  trackMan = new CDSTrackingMan();  
  if( trackMan==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  evTree->Branch( "CDSTrackingMan", &trackMan );

  t0=time(0);
  AllGoodTrack=0;
  nTrack=0;
}

void EventAnalysis::USca( int nsca, unsigned int *sca )
{
#if 0
  std::cout << " Enter EventAnalysis::USca " << std::endl;
#endif
  Block_Event_Number++;
  header->SetBlockEventNumber( Block_Event_Number );
}

bool EventAnalysis::UAna( TKOHitCollection *tko )
{
#if 0
  std::cout << " Enter EventAnalysis::UAna " << std::endl;
#endif

  Event_Number++;
  { int status = confMan->CheckEvNum( Event_Number, Block_Event_Number );
    if( status==1 ) return true;
    if( status==2 ) return false; }

  if( Event_Number%1000==0)
    {
      t1=time(0);
      std::cout << " Event# : " << Event_Number << "  BlockEvent# : " << Block_Event_Number << " AllTrack# : " << nTrack <<" GoodTrack# : " << AllGoodTrack << " Time (s): " << (t1-t0) << std::endl;
    }
  
  header->SetRunNumber(0);
  header->SetEventNumber(Event_Number);
  
  header->Convert( tko, confMan );

  if( header->IsTrig(Trig_Cosmic) ){
    header->Clear();
    cdsMan->Clear();
    trackMan->Clear();
    return true;
  }
  cdsMan->Convert( tko, confMan );
  trackMan->Execute(cdsMan,confMan);

  int nGoodTrack=trackMan->nGoodTrack();
  int nallTrack=  trackMan->nTrack();
  AllGoodTrack+=nGoodTrack;
  nTrack+=nallTrack;
//  if( nGoodTrack<1 ){
  if( nallTrack<1 ){
    header->Clear();
    cdsMan->Clear();
    trackMan->Clear();
    return true;
  }
  cdsMan->Clear();
  //  trackMan->Calc(0,confMan);
  evTree->Fill();

  header->Clear();
  trackMan->Clear();
  return true;
}

void EventAnalysis::Finalize()
{
  t1=time(0);
  std::cout<<" Enter EventAnalysis::Finalize "<< std::endl;
  std::cout << " Event# : " << Event_Number << "  BlockEvent# : " << Block_Event_Number << " AllTrack# : " << nTrack <<" GoodTrack# : " << AllGoodTrack << " Time (s): " << (t1-t0) << std::endl;

  rtFile->cd();
  confMan->SaveCDSParam();
  gFile->Write();
  gFile->Close();
  //  gSystem->Sleep(10000);
  //  gSystem->Exec(Form("mv %s %s",tmpname.Data(),confMan->GetOutFileName().c_str()));
  delete cdsMan;
  delete trackMan;
  delete header;
}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysis *event = new EventAnalysis();
  return (EventTemp*)event;
}
