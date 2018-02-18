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

    // Start
    for( int i=0; i<blMan->nStart(); i++ ){
      HodoscopeLikeHit *hit = blMan->Start(i);
      int seg = hit->seg();
      int tu = hit->tdcu();
      double timeu = hit->time(0);
      double etime = (double)Event_Number/60.;
#if Debug
      std::cout << "= trg =" << std::endl;
      std::cout << "seg = " << seg << std::endl;
#endif

        h1 = (TH1F*)gFile->Get( Form("Start%d_TDC",seg) ); h1->Fill( tu );
        h1 = (TH1F*)gFile->Get( Form("Start%d_Time",seg) ); h1->Fill( timeu );
        h2 = (TH2F*)gFile->Get( Form("Start%d_TDCvEvent",seg) ); h2->Fill( etime,tu );
        h2 = (TH2F*)gFile->Get( Form("Start%d_TimevEvent",seg) ); h2->Fill( etime,timeu );
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
  const int NumOfStartSegments = 16;

  rtFile->cd();

  // Scaler
  new TH1F( "Scaler", "Scaler", 50, 0, 50 );

  // Trigger Pattern
  new TH1F( "Pattern", "Trigger Pattern", 20, 0, 20 );

  // Start
  std::cout << "Define Histograms for Start" << std::endl;
  TString startname[NumOfStartSegments]=
  {"T0 start","BHD start","BVC start","BPD start","CDH start","IH start",
    "NC start","PC start","CVC start","T0-03U","T0-03D","T0-OR (Before retiming)",
    "T0-RT","T0-RT-BL (After Fin/Fout 9-1-1)","T0-RT-FWD (After Fin/Fout 11-1-4)","T0-OR (Self)"};
new TH1F( "Start_HitPat", "Hit Pattern Start", NumOfStartSegments+1, 0, NumOfStartSegments+1 );
for( int seg=1; seg<=NumOfStartSegments; seg++ ){
  new TH1F( Form("Start%d_TDC",seg),   Form("%s;TDC ch.;Counts",startname[seg-1].Data()),    4097,    -0.5, 4096.5 );
  new TH1F( Form("Start%d_Time",seg),   Form("%s;Time (ns);Counts",startname[seg-1].Data()),    4001,    -0.0125, 100.0125 );
  new TH2F( Form("Start%d_TDCvEvent",seg),   Form("%s;Elapsed time (min.);TDC ch.",startname[seg-1].Data()),    60, 0, 60, 4097,    -0.5, 4096.5 );
  new TH2F( Form("Start%d_TimevEvent",seg),   Form("%s;Elapsed time (min.);Time (ns)",startname[seg-1].Data()),    60, 0, 60, 4001,    -0.0125, 100.0125 );
}

}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysisMyStart *event = new EventAnalysisMyStart();
  return (EventTemp*)event;
}
