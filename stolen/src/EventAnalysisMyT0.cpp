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

class EventAnalysisMyT0: public EventTemp
{
public:
  EventAnalysisMyT0();
  ~EventAnalysisMyT0();
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

EventAnalysisMyT0::EventAnalysisMyT0()
  : EventTemp()
{
}

EventAnalysisMyT0::~EventAnalysisMyT0()
{
}

const int MaxTreeSize = 1900000000000;
void EventAnalysisMyT0::Initialize( ConfMan *conf )
{
#if 1
  std::cout << " Enter EventAnalysisMyT0::Initialize " << std::endl;
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

void EventAnalysisMyT0::USca( int nsca, unsigned int *sca )
{
#if 0
  std::cout << " Enter EventAnalysisMyT0::USca " << std::endl;
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

bool EventAnalysisMyT0::UAna( TKOHitCollection *tko )
{
#if 0
	std::cout << " Enter EventAnalysisMyT0::UAna " << std::endl;
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

    if(header->trigmode2(Mode_KCDH2,confMan))
      std::cout << "###!!!###" << std::endl;

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

    // Selection //
    if(MulT0!=1){
      header->Clear();
      blMan->Clear();
      cdsMan->Clear();
      return true;
    }

    // Trigger Pattern
    for( int i=0; i<20; i++ ){
      int val = header->pattern(i);
      if( 0<val ){
        h1 = (TH1F*)gFile->Get("Pattern"); h1->Fill(i);
      }
    }


    // T0
    int nT0=0;
    for( int i=0; i<blMan->nT0(); i++ ){
      HodoscopeLikeHit *hit = blMan->T0(i);
      int seg = hit->seg();
      int au = hit->adcu(), ad = hit->adcd();
      int tu = hit->tdcu(), td = hit->tdcd();
      double hitpos = hit->hitpos();
      double timeu = hit->time(0), timed = hit->time(1);
      double ctimeu = hit->ctime(0), ctimed = hit->ctime(1);
      double eneu = hit->ene(0), ened = hit->ene(1);
      double emean = hit->emean();
      double tmean = hit->tmean();
      double ctmean = hit->ctmean();
#if Debug
      std::cout << "= trg =" << std::endl;
      std::cout << "seg = " << seg << std::endl;
#endif

      h1 = (TH1F*)gFile->Get( Form("T0u%d_ADC",seg) ); h1->Fill( au );
      h1 = (TH1F*)gFile->Get( Form("T0d%d_ADC",seg) ); h1->Fill( ad );
      if(tu>0){
        h1 = (TH1F*)gFile->Get( Form("T0u%d_TDC",seg) ); h1->Fill( tu );
        h1 = (TH1F*)gFile->Get( Form("T0u%d_Time",seg) ); h1->Fill( timeu );
        h1 = (TH1F*)gFile->Get( Form("T0u%d_CTime",seg) ); h1->Fill( ctimeu );
      }
      if(td>0){
        h1 = (TH1F*)gFile->Get( Form("T0d%d_TDC",seg) ); h1->Fill( td );
        h1 = (TH1F*)gFile->Get( Form("T0d%d_Time",seg) ); h1->Fill( timed );
        h1 = (TH1F*)gFile->Get( Form("T0d%d_CTime",seg) ); h1->Fill( ctimed );
      }
      if( hit->CheckRange() ){
        h1 = (TH1F*)gFile->Get( Form("T0u%d_ADCwT",seg) ); h1->Fill( au );
        h1 = (TH1F*)gFile->Get( Form("T0d%d_ADCwT",seg) ); h1->Fill( ad );
        h1 = (TH1F*)gFile->Get( Form("T0u%d_dE",seg) ); h1->Fill( eneu );
        h1 = (TH1F*)gFile->Get( Form("T0d%d_dE",seg) ); h1->Fill( ened );
        h1 = (TH1F*)gFile->Get( Form("T0%d_dEMean",seg) ); h1->Fill( emean );
        h1 = (TH1F*)gFile->Get( Form("T0%d_TimeMean",seg) ); h1->Fill( tmean );
        h1 = (TH1F*)gFile->Get( Form("T0%d_CTimeMean",seg) ); h1->Fill( ctmean );
        h1 = (TH1F*)gFile->Get( Form("T0%d_Position",seg) ); h1->Fill( hitpos );
        if((eneu<1.5||eneu>2.5)&&(ened<1.5||ened>2.5)){
          h2 = (TH2F*)gFile->Get( Form("T0u%d_AvT",seg) );  h2->Fill( au, tu );
          h2 = (TH2F*)gFile->Get( Form("T0d%d_AvT",seg) );  h2->Fill( ad, td );
          h2 = (TH2F*)gFile->Get( Form("T0u%d_dEvTime",seg) );  h2->Fill( eneu, timeu );
          h2 = (TH2F*)gFile->Get( Form("T0d%d_dEvTime",seg) );  h2->Fill( ened, timed );
          h2 = (TH2F*)gFile->Get( Form("T0u%d_dEvCTime",seg) );  h2->Fill( eneu, ctimeu );
          h2 = (TH2F*)gFile->Get( Form("T0d%d_dEvCTime",seg) );  h2->Fill( ened, ctimed );
          h2 = (TH2F*)gFile->Get( Form("T0u%d_dEvCTimeMean",seg) );  h2->Fill( eneu, ctmean );
          h2 = (TH2F*)gFile->Get( Form("T0d%d_dEvCTimeMean",seg) );  h2->Fill( ened, ctmean );
        }
        h2 = (TH2F*)gFile->Get( Form("T0%d_dEMeanvCTimeMean",seg) );  h2->Fill( emean, ctmean );
        h1 = (TH1F*)gFile->Get( "T0_HitPat" ); h1->Fill( seg );
        nT0++;
      }
    }
    h1 = (TH1F*)gFile->Get( "T0_Mult" ); h1->Fill( nT0 );

    header->Clear();
    blMan->Clear();
    cdsMan->Clear();
    return true;
}

void EventAnalysisMyT0::Finalize()
{
  std::cout << " Enter EventAnalysisMyT0::Finalize " << std::endl;

  rtFile->cd();
  gFile->Write();
  gFile->Close();

  delete blMan;
  delete cdsMan;
  delete header;
}

void EventAnalysisMyT0::InitializeHistogram()
{
  Int_t NumOfT0Segments = 5;

  rtFile->cd();

  // Scaler
  new TH1F( "Scaler", "Scaler", 50, 0, 50 );

  // Trigger Pattern
  new TH1F( "Pattern", "Trigger Pattern", 20, 0, 20 );

  // T0
  std::cout << "Define Histograms for T0" << std::endl;
  new TH1F( "T0_HitPat", "Hit Pattern T0", NumOfT0Segments+1, 0, NumOfT0Segments+1 );
  for( int seg=1; seg<=NumOfT0Segments; seg++ ){
    new TH1F( Form("T0u%d_ADC",seg),   Form("ADC T0U%d;ADC ch.;Counts",seg),    4000,    0, 4000 );
    new TH1F( Form("T0d%d_ADC",seg),   Form("ADC T0D%d;ADC ch.;Counts",seg),    4000,    0, 4000 );
    new TH1F( Form("T0u%d_dE",seg),   Form("dE T0U%d;dE (MeV);Counts",seg),    400,    0, 40 );
    new TH1F( Form("T0d%d_dE",seg),   Form("dE T0D%d;dE (MeV);Counts",seg),    400,    0, 40 );
    new TH1F( Form("T0%d_dEMean",seg),   Form("Mean dE T0%d;dE (MeV);Counts",seg),    400,    0, 40 );
    new TH1F( Form("T0u%d_TDC",seg),   Form("TDC T0U%d;TDC ch.;Counts",seg),    4000,    0, 4000 );
    new TH1F( Form("T0d%d_TDC",seg),   Form("TDC T0D%d;TDC ch.;Counts",seg),    4000,    0, 4000 );
    new TH1F( Form("T0u%d_Time",seg),   Form("Time T0U%d;Time (ns);Counts",seg),    4000,    -100, 100 );
    new TH1F( Form("T0d%d_Time",seg),   Form("Time T0D%d;Time (ns);Counts",seg),    4000,    -100, 100 );
    new TH1F( Form("T0%d_TimeMean",seg),   Form("Mean Time T0%d;Time (ns);Counts",seg),    4000,    -100, 100 );
    new TH1F( Form("T0u%d_CTime",seg),   Form("CTime T0U%d;Time (ns);Counts",seg),    4000,    -100, 100 );
    new TH1F( Form("T0d%d_CTime",seg),   Form("CTime T0D%d;Time (ns);Counts",seg),    4000,    -100, 100 );
    new TH1F( Form("T0%d_CTimeMean",seg),   Form("Mean CTime T0%d;Time (ns);Counts",seg),    4000,    -100, 100 );
    new TH1F( Form("T0%d_Position",seg),   Form("Hit position T0%d;Position (cm);Counts",seg),    500,    -50, 50 );
    new TH1F( Form("T0u%d_ADCwT",seg), Form("ADC wTDC T0U%d;ADC ch.;Counts",seg),  4000,    0, 4000 );
    new TH1F( Form("T0d%d_ADCwT",seg), Form("ADC wTDC T0D%d;ADC ch.;Counts",seg),  4000,    0, 4000 );
    new TH2F( Form("T0u%d_AvT",seg),   Form("ADC TDC corr. T0U%d;ADC ch.;TDC ch.",seg),     200,    0, 4000,  200,    0, 4000 );
    new TH2F( Form("T0d%d_AvT",seg),   Form("ADC TDC corr. T0D%d;ADC ch.;TDC ch.",seg),     200,    0, 4000,  200,    0, 4000 );
    new TH2F( Form("T0u%d_dEvTime",seg),   Form("dE Time corr. T0U%d;dE (MeV);Time (ns)",seg),     400,    0, 40,  4000,    -100, 100 );
    new TH2F( Form("T0d%d_dEvTime",seg),   Form("dE Time corr. T0D%d;dE (MeV);Time (ns)",seg),     400,    0, 40,  4000,    -100, 100 );
    new TH2F( Form("T0u%d_dEvCTime",seg),   Form("dE CTime corr. T0U%d;dE (MeV);Time (ns)",seg),     400,    0, 40,  4000,    -100, 100 );
    new TH2F( Form("T0d%d_dEvCTime",seg),   Form("dE CTime corr. T0D%d;dE (MeV);Time (ns)",seg),     400,    0, 40,  4000,    -100, 100 );
    new TH2F( Form("T0u%d_dEvCTimeMean",seg),   Form("dE CTimeMean corr. T0U%d;dE (MeV);Time (ns)",seg),     400,    0, 40,  4000,    -100, 100 );
    new TH2F( Form("T0d%d_dEvCTimeMean",seg),   Form("dE CTimeMean corr. T0D%d;dE (MeV);Time (ns)",seg),     400,    0, 40,  4000,    -100, 100 );
    new TH2F( Form("T0%d_dEMeanvCTimeMean",seg),   Form("dEMean CTimeMean corr. T0%d;dE (MeV);Time (ns)",seg),     400,    0, 40,  4000,    -100, 100 );
  }
  new TH1F( "T0_Mult", "Multiplicity T0;Multiplicity;Counts", NumOfT0Segments+1, 0, NumOfT0Segments+1 );

}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysisMyT0 *event = new EventAnalysisMyT0();
  return (EventTemp*)event;
}
