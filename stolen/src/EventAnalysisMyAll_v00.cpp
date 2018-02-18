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
#include "MyAnalysisBL.h"
#include "MyAnalysisCDS.h"

#define Debug 0

class EventAnalysisMyAll: public EventTemp
{
  public:
    EventAnalysisMyAll();
    ~EventAnalysisMyAll();
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
    MyAnalysisCDS* cdsAna;

    int AllGoodTrack;
    int nAllTrack;
    int CDC_Event_Number;
    int t0, t1;

  public:
    void Initialize( ConfMan *conf );
    void USca( int nsca, unsigned int *sca );
    bool UAna( TKOHitCollection *tko );
    void Finalize();

};

  EventAnalysisMyAll::EventAnalysisMyAll()
: EventTemp()
{
}

EventAnalysisMyAll::~EventAnalysisMyAll()
{
}

const int MaxTreeSize = 1900000000000;
void EventAnalysisMyAll::Initialize( ConfMan *conf )
{
#if 1
  std::cout << " Enter EventAnalysisMyAll::Initialize " << std::endl;
#endif
  confMan = conf;

  TString cdcfname=confMan->GetCDSTrackFileName();
  cdcFile = new TFile(cdcfname);
  if(!cdcFile->IsOpen()){
    std::cout<<" failed to open " <<cdcfname<< "  !!!"<<std::endl;
    exit(false);
  }
  std::cout<<" CDC Track File Successfully opend !!! " <<cdcfname<<std::endl;
  cdstrackMan   = new CDSTrackingMan();
  cdsheader    = new EventHeader();
  cdcTree=(TTree*)cdcFile->Get("EventTree");
  cdcTree->SetBranchAddress( "CDSTrackingMan", &cdstrackMan );
  cdcTree->SetBranchAddress( "EventHeader" ,&cdsheader);

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
  cdsAna = new MyAnalysisCDS(rtFile, confMan);

  AllGoodTrack = 0;
  nAllTrack = 0;
  CDC_Event_Number = 0;
  t0=time(0);
}

void EventAnalysisMyAll::USca( int nsca, unsigned int *sca )
{
  Block_Event_Number++;
  header->SetBlockEventNumber( Block_Event_Number );
  scaMan->SetBlockEventNumber( Block_Event_Number );
  for( int i=0; i<nsca; i++ ){
    scaMan->AddHit( sca[i], confMan->GetScalerMapManager()->GetName(i) );
  }

  scaMan->Clear();
}

bool EventAnalysisMyAll::UAna( TKOHitCollection *tko )
{
  Event_Number++;
  { int status = confMan->CheckEvNum( Event_Number, Block_Event_Number );
    if( status==1 ) return true;
    if( status==2 ) return false; }

    if( Event_Number%1000==0 )
    {
      t1=time(0);
      std::cout << " Event# : " << Event_Number << "  BlockEvent# : " << Block_Event_Number << " AllTrack# : " << nAllTrack <<" GoodTrack# : " << AllGoodTrack  <<"  "<<CDC_Event_Number<<" / "<<cdcTree->GetEntries() << " Time (s):" << (t1-t0) << std::endl;
    }

    if(CDC_Event_Number>=cdcTree->GetEntries()) return false;  
    cdcTree->GetEntry(CDC_Event_Number);
    if(cdsheader->ev()!=Event_Number){
      int tmpnum=CDC_Event_Number;
      while(cdsheader->ev()<Event_Number&&tmpnum<=cdcTree->GetEntries()){
        tmpnum++;    
        cdcTree->GetEntry(tmpnum);
      }
      if(cdsheader->ev()!=Event_Number)  return true;
      else CDC_Event_Number=tmpnum;
    }
    CDC_Event_Number++;

    header->SetRunNumber(0);
    header->SetEventNumber(Event_Number);

    header->Convert( tko, confMan );
    blMan->Convert( tko, confMan );
    cdsMan->Convert( tko, confMan );
    cdstrackMan -> Calc( cdsMan, confMan );
    int nGoodTrack = cdstrackMan -> nGoodTrack();
    int nAllTrack = cdstrackMan -> nTrack();
    AllGoodTrack += nGoodTrack;
    nAllTrack += nAllTrack;

    rtFile->cd();
    TH1F *h1;
    TH2F *h2;
    DetectorList *dlist=DetectorList::GetInstance();

    // ======== //
    // Raw Data //
    // ======== //

//    // Selection 1//
//    if(!header->IsTrig(Trig_Kaon)){
//      blAna->Clear();
//      cdsAna->Clear();
//      header->Clear();
//      blMan->Clear();
//      bltrackMan->Clear();
//      cdsMan->Clear();
//      return true;
//    }

//    bltrackMan -> DoTracking(blMan, confMan, true, false);

//    int MulTrackBLC1 = bltrackMan -> ntrackBLC1();
//    int MulTrackBLC2 = bltrackMan -> ntrackBLC2();

 //   // Selection 2//
 //   if(MulTrackBLC1!=1||MulTrackBLC2!=1){
 //     header->Clear();
 //     blMan->Clear();
 //     bltrackMan->Clear();
 //     cdsMan->Clear();
 //     return true;
 //   }

    if(!blAna->DoAnalysis(confMan, header, blMan, bltrackMan)){
      blAna->Clear();
      cdsAna->Clear();
      header->Clear();
      blMan->Clear();
      bltrackMan->Clear();
      cdsMan->Clear();
      return true;
    }
    cdsAna->SetBeamPID(blAna->GetBeamPID());
    cdsAna->SetBeamMom(blAna->GetBeamMom());
    cdsAna->SetBeamL(blAna->GetBeamL());
    cdsAna->SetBLC1(blAna->GetBLC1());
    cdsAna->SetBLC2(blAna->GetBLC2());
    cdsAna->SetBPC(blAna->GetBPC());
    cdsAna->DoAnalysis(confMan, header, blMan, bltrackMan, cdsMan, cdstrackMan);

    blAna->Clear();
    cdsAna->Clear();

    header->Clear();
    blMan->Clear();
    bltrackMan->Clear();
    cdsMan->Clear();
    return true;
}

void EventAnalysisMyAll::Finalize()
{
  std::cout << " Enter EventAnalysisMyAll::Finalize " << std::endl;

  rtFile->cd();
  gFile->Write();
  gFile->Close();

  delete blAna;
  delete cdsAna;
  delete blMan;
  delete bltrackMan;
  delete cdsMan;
  delete header;
}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysisMyAll *event = new EventAnalysisMyAll();
  return (EventTemp*)event;
}
