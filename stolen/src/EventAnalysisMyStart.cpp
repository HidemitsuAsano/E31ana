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

class EventAnalysisMyStart: public EventTemp
{
public:
  EventAnalysisMyStart();
  ~EventAnalysisMyStart();
private:
  TFile *rtFile;
  TTree *evTree;
  TTree *scaTree;
  CDSHitMan *cdsMan;
  BeamLineHitMan *blMan;
  EventHeader *header;
  ScalerMan *scaMan;
  double Spill;
public:
  void Initialize( ConfMan *conf );
  void USca( int nsca, unsigned int *sca );
  bool UAna( TKOHitCollection *tko );
  void Finalize();

  void InitializeHistogram();
};

EventAnalysisMyStart::EventAnalysisMyStart()
  : EventTemp()
{
}

EventAnalysisMyStart::~EventAnalysisMyStart()
{
}

const int MaxTreeSize = 1900000000000;
void EventAnalysisMyStart::Initialize( ConfMan *conf )
{
#if 1
  std::cout << " Enter EventAnalysisMyStart::Initialize " << std::endl;
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
}

void EventAnalysisMyStart::USca( int nsca, unsigned int *sca )
{
#if 0
  std::cout << " Enter EventAnalysisMyStart::USca " << std::endl;
#endif
  Block_Event_Number++;
  header->SetBlockEventNumber( Block_Event_Number );
  scaMan->SetBlockEventNumber( Block_Event_Number );
  for( int i=0; i<nsca; i++ ){
    scaMan->AddHit( sca[i], confMan->GetScalerMapManager()->GetName(i) );
  }
#if 0
  std::cout << nsca << std::endl;
  for( int i=0; i<nsca; i++ ){
   std::cout << "  " << sca[i];
  }
  std::cout << std::endl;
#endif

  TH1F *h1;
  for( int i=0; i<scaMan->nsca(); i++ ){
    int val = scaMan->sca(i)->val();
    TString name = scaMan->sca(i)->name();
    h1 = (TH1F*)gFile->Get("Scaler"); h1->Fill(i,val);
  }
  double val = (double)scaMan->sca(0)->val();
  Spill = val; h1=(TH1F*)gFile->Get("SpillNum"); h1->SetBinContent(1,val);

  scaMan->Clear();
}

bool EventAnalysisMyStart::UAna( TKOHitCollection *tko )
{
#if 0
	std::cout << " Enter EventAnalysisMyStart::UAna " << std::endl;
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

    // ======== //
    // Raw Data //
    // ======== //


    // Multiplicity //

    // Selection //

    // T0#3
    for( int i=0; i<blMan->nT0(); i++ ){
      HodoscopeLikeHit *hit = blMan->T0(i);
      int seg = hit->seg();
      int tu = hit->tdcu(), td = hit->tdcd();
      int ctimeu = hit->ctime(0), ctimed = hit->ctime(1);
#if Debug
      std::cout << "= trg =" << std::endl;
      std::cout << "seg = " << seg << std::endl;
#endif
      h1 = (TH1F*)gFile->Get( Form("T0u%d_TDC",seg) ); h1->Fill( tu );
      h2 = (TH2F*)gFile->Get( Form("T0u%d_TDCvEvent",seg) ); h2->Fill( Spill,tu );
      h1 = (TH1F*)gFile->Get( Form("T0d%d_TDC",seg) ); h1->Fill( td );
      h2 = (TH2F*)gFile->Get( Form("T0d%d_TDCvEvent",seg) ); h2->Fill( Spill,td );
      h1 = (TH1F*)gFile->Get( Form("T0u%d_Time",seg) ); h1->Fill( ctimeu );
      h2 = (TH2F*)gFile->Get( Form("T0u%d_TimevEvent",seg) ); h2->Fill( Spill,ctimeu );
      h1 = (TH1F*)gFile->Get( Form("T0d%d_Time",seg) ); h1->Fill( ctimed );
      h2 = (TH2F*)gFile->Get( Form("T0d%d_TimevEvent",seg) ); h2->Fill( Spill,td );
    }

    // Start
    for( int i=0; i<blMan->nStart(); i++ ){
      HodoscopeLikeHit *hit = blMan->Start(i);
      int seg = hit->seg();
      int tu = hit->tdcu();
      int timeu = hit->time(0);
#if Debug
      std::cout << "= trg =" << std::endl;
      std::cout << "seg = " << seg << std::endl;
#endif
      h1 = (TH1F*)gFile->Get( Form("Start%d_TDC",seg) ); h1->Fill( tu );
      h2 = (TH2F*)gFile->Get( Form("Start%d_TDCvEvent",seg) ); h2->Fill( Spill,tu );
      h1 = (TH1F*)gFile->Get( Form("Start%d_Time",seg) ); h1->Fill( timeu );
      h2 = (TH2F*)gFile->Get( Form("Start%d_TimevEvent",seg) ); h2->Fill( Spill,timeu );
    }

    header->Clear();
    blMan->Clear();
    cdsMan->Clear();
    return true;
}

void EventAnalysisMyStart::Finalize()
{
  std::cout << " Enter EventAnalysisMyStart::Finalize " << std::endl;

  rtFile->cd();
  gFile->Write();
  gFile->Close();

  delete blMan;
  delete cdsMan;
  delete header;
}

void EventAnalysisMyStart::InitializeHistogram()
{
  const int NumOfStartSegments = 18;

  rtFile->cd();

  // Scaler
  new TH1F( "Scaler", "Scaler", 50, 0, 50 );
  new TH1F( "SpillNum", "SpillNum", 1, 0, 1 );

  // Trigger Pattern
  new TH1F( "Pattern", "Trigger Pattern", 20, 0, 20 );

  // T0#3
  std::cout << "Define Histograms for T0" << std::endl;
  for( int seg=1; seg<=5; seg++ ){
    new TH1F( Form("T0u%d_TDC",seg),   Form("TDC T0U%d;TDC ch.;Counts",seg),    4001,    -0.5, 4000.5 );
    new TH1F( Form("T0d%d_TDC",seg),   Form("TDC T0D%d;TDC ch.;Counts",seg),    4001,    -0.5, 4000.5 );
    new TH2F( Form("T0u%d_TDCvEvent",seg),   Form("TDC T0U%d;Spill;TDC ch.",seg),    1000, 0, 1000,  4001,    -0.5, 4000.5 );
    new TH2F( Form("T0d%d_TDCvEvent",seg),   Form("TDC T0D%d;Spill;TDC ch.",seg),    1000, 0, 1000, 4001,    -0.5, 4000.5 );
    new TH1F( Form("T0u%d_Time",seg),   Form("Time T0U%d;Time (ns);Counts",seg),    4001,    -50., 50. );
    new TH1F( Form("T0d%d_Time",seg),   Form("Time T0D%d;Time (ns);Counts",seg),    4001,    -50., 50. );
    new TH2F( Form("T0u%d_TimevEvent",seg),   Form("Time T0U%d;Spill;Time (ns)",seg),    1000, 0, 1000,  4001,    -50.0125, 50.0125 );
    new TH2F( Form("T0d%d_TimevEvent",seg),   Form("Time T0D%d;Spill;Time (ns)",seg),    1000, 0, 1000, 4001,    -50.0125, 50.0125 );
  }

  // Start
  std::cout << "Define Histograms for Start" << std::endl;
  TString startname[NumOfStartSegments]={"T0 start (c2s4)","NC start (c2s7)","PC start (c2s7)",
    "CVC start (c2s7)","BVC start (c2c7)","T0 start (c2s7)","T0 start (c6s9)","T0 start (c8s1)","BVC start (c8s2)",
    "T0 start (c8s2)","T0 start (c8s4)","BVC start (c8s5)","PC start (c8s5)","T0 start (c8s15)","T0 start (c8s17)",
    "T0-OR (c2s4)","T0-retiming (c2s4)","T0-retiming after Fin/Fout (c2s4)"}; 
  new TH1F( "Start_HitPat", "Hit Pattern Start", NumOfStartSegments+1, 0, NumOfStartSegments+1 );
  for( int seg=1; seg<=NumOfStartSegments; seg++ ){
    new TH1F( Form("Start%d_TDC",seg),   Form("%s;TDC ch.;Counts",startname[seg-1].Data()),    4001,    -0.5, 4000.5 );
    new TH2F( Form("Start%d_TDCvEvent",seg),   Form("%s;Spill;TDC ch.",startname[seg-1].Data()),    1000, 0, 1000, 4001,    -0.5, 4000.5 );
    new TH1F( Form("Start%d_Time",seg),   Form("%s;Time (ns);Counts",startname[seg-1].Data()),    4001,    -0.05, 400.05  );
    new TH2F( Form("Start%d_TimevEvent",seg),   Form("%s;Spill;Time (ns)",startname[seg-1].Data()),    1000, 0, 1000, 4001,    -0.05, 400.05 );
  }

}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysisMyStart *event = new EventAnalysisMyStart();
  return (EventTemp*)event;
}
