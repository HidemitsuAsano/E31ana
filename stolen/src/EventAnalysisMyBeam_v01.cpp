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

#define Debug 0

class EventAnalysisMyBeam: public EventTemp
{
  public:
    EventAnalysisMyBeam();
    ~EventAnalysisMyBeam();
  private:
    TFile *rtFile;
    TFile *rtFile2;
    TTree *evTree;
    TTree *scaTree;
    EventHeader *header;
    ScalerMan *scaMan;
    CDSHitMan *cdsMan;
    CDSTrackingMan *cdstrackMan;
    BeamLineHitMan *blMan;
    BeamLineTrackMan *bltrackMan;
    MyAnalysisBL* blAna;

    int AllGoodTrack;
    int nAllTrack;
    int nTrack;
    int CDC_Event_Number;
    bool trackflag;
    int MTDC_Event_Number;
    int t0, t1;

  public:
    void Initialize( ConfMan *conf );
    void USca( int nsca, unsigned int *sca );
    bool UAna( TKOHitCollection *tko );
    void Finalize();

    void InitializeHistogram();
};

  EventAnalysisMyBeam::EventAnalysisMyBeam()
: EventTemp()
{
}

EventAnalysisMyBeam::~EventAnalysisMyBeam()
{
}

const int MaxTreeSize = 1900000000000;
void EventAnalysisMyBeam::Initialize( ConfMan *conf )
{
#if 1
  std::cout << " Enter EventAnalysisMyBeam::Initialize " << std::endl;
#endif
  confMan = conf;

  rtFile = new TFile( confMan->GetOutFileName().c_str(), "recreate" );
  InitializeHistogram();

  header = new EventHeader();
  if( header==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  cdsMan = new CDSHitMan();
  if( cdsMan==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  blMan = new BeamLineHitMan();
  if( blMan ==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  scaMan = new ScalerMan();
  if( scaMan==NULL ){ std::cerr << "!!!!" << std::endl; return; }

  std::string ofname = confMan->GetOutFileName();
  ofname.insert(ofname.find(".root"),"_anaBeam");
  rtFile2 =  new TFile( Form("%s",ofname.c_str()), "RECREATE");
  rtFile2->cd();
  blAna = new MyAnalysisBL(rtFile2, confMan);

  AllGoodTrack = 0;
  nAllTrack = 0;
  CDC_Event_Number = 0;
  MTDC_Event_Number = 0;
  t0=time(0);
}

void EventAnalysisMyBeam::USca( int nsca, unsigned int *sca )
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

bool EventAnalysisMyBeam::UAna( TKOHitCollection *tko )
{
#if 0
  std::cout << " Enter EventAnalysisMyBeam::UAna " << std::endl;
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

    rtFile->cd();
    TH1F *h1;
    TH2F *h2;
    DetectorList *dlist=DetectorList::GetInstance();

    // ======== //
    // Raw Data //
    // ======== //


    header->Clear();
    blMan->Clear();
    cdsMan->Clear();
    return true;
}

void EventAnalysisMyBeam::Finalize()
{
  std::cout << " Enter EventAnalysisMyBeam::Finalize " << std::endl;

  rtFile->cd();
  gFile->Write();
  gFile->Close();

  delete blMan;
  delete cdsMan;
  delete header;
}

void EventAnalysisMyBeam::InitializeHistogram()
{
  rtFile->cd();

  // Scaler
  new TH1F( "Scaler", "Scaler", 50, 0, 50 );
  // Trigger
  TH1F* h1 = new TH1F( "TriggerPattern", "Trigger Pattern", 18, -0.5, 17.5);
  h1->GetXaxis()->SetBinLabel(1,"All");
  h1->GetXaxis()->SetBinLabel(2,"Beam");
  h1->GetXaxis()->SetBinLabel(3,"Kaon");
  h1->GetXaxis()->SetBinLabel(4,"KCDH1f");
  h1->GetXaxis()->SetBinLabel(5,"Pion");
  h1->GetXaxis()->SetBinLabel(6,"Proton");
  h1->GetXaxis()->SetBinLabel(7,"KCDH1");
  h1->GetXaxis()->SetBinLabel(8,"KCDH2");
  h1->GetXaxis()->SetBinLabel(9,"PivBVC");
  h1->GetXaxis()->SetBinLabel(10,"PiCDH1");
  h1->GetXaxis()->SetBinLabel(11,"PiCDH2");
  h1->GetXaxis()->SetBinLabel(12,"Kf");
  h1->GetXaxis()->SetBinLabel(13,"1stMix");
  h1->GetXaxis()->SetBinLabel(14,"Charged");
  h1->GetXaxis()->SetBinLabel(15,"Neutral");
  h1->GetXaxis()->SetBinLabel(16,"Cosmic");
  h1->GetXaxis()->SetBinLabel(17,"Reject");
  h1->GetXaxis()->SetBinLabel(18,"SIM");


}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysisMyBeam *event = new EventAnalysisMyBeam();
  return (EventTemp*)event;
}
