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

class EventAnalysisMyPC: public EventTemp
{
public:
  EventAnalysisMyPC();
  ~EventAnalysisMyPC();
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

EventAnalysisMyPC::EventAnalysisMyPC()
  : EventTemp()
{
}

EventAnalysisMyPC::~EventAnalysisMyPC()
{
}

const int MaxTreeSize = 1900000000000;
void EventAnalysisMyPC::Initialize( ConfMan *conf )
{
#if 1
  std::cout << " Enter EventAnalysisMyPC::Initialize " << std::endl;
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

void EventAnalysisMyPC::USca( int nsca, unsigned int *sca )
{
#if 0
  std::cout << " Enter EventAnalysisMyPC::USca " << std::endl;
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

bool EventAnalysisMyPC::UAna( TKOHitCollection *tko )
{
#if 0
	std::cout << " Enter EventAnalysisMyPC::UAna " << std::endl;
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

		// Multihit of PC //
		int MulPC=0;	
		for(int i=0; i<blMan->nPC(); i++){
			HodoscopeLikeHit* hit = blMan->PC(i);
			if(hit->CheckRange()) MulPC++;
		}

		// Selection //
		if(MulPC==0){
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

    // PC
    int nPC=0;
    for( int i=0; i<blMan->nPC(); i++ ){
      HodoscopeLikeHit *hit = blMan->PC(i);
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

      h1 = (TH1F*)gFile->Get( Form("PCu%d_ADC",seg) ); h1->Fill( au );
      h1 = (TH1F*)gFile->Get( Form("PCd%d_ADC",seg) ); h1->Fill( ad );
      if(tu>0){
        h1 = (TH1F*)gFile->Get( Form("PCu%d_TDC",seg) ); h1->Fill( tu );
        h1 = (TH1F*)gFile->Get( Form("PCu%d_Time",seg) ); h1->Fill( timeu );
        h1 = (TH1F*)gFile->Get( Form("PCu%d_CTime",seg) ); h1->Fill( ctimeu );
      }
      if(td>0){
        h1 = (TH1F*)gFile->Get( Form("PCd%d_TDC",seg) ); h1->Fill( td );
        h1 = (TH1F*)gFile->Get( Form("PCd%d_Time",seg) ); h1->Fill( timed );
        h1 = (TH1F*)gFile->Get( Form("PCd%d_CTime",seg) ); h1->Fill( ctimed );
      }
      if( hit->CheckRange() ){
        h1 = (TH1F*)gFile->Get( Form("PCu%d_ADCwT",seg) ); h1->Fill( au );
        h1 = (TH1F*)gFile->Get( Form("PCd%d_ADCwT",seg) ); h1->Fill( ad );
        h1 = (TH1F*)gFile->Get( Form("PCu%d_dE",seg) ); h1->Fill( eneu );
        h1 = (TH1F*)gFile->Get( Form("PCd%d_dE",seg) ); h1->Fill( ened );
        h1 = (TH1F*)gFile->Get( Form("PC%d_dEMean",seg) ); h1->Fill( emean );
        h1 = (TH1F*)gFile->Get( Form("PC%d_TimeMean",seg) ); h1->Fill( tmean );
        h1 = (TH1F*)gFile->Get( Form("PC%d_CTimeMean",seg) ); h1->Fill( ctmean );
        h1 = (TH1F*)gFile->Get( Form("PC%d_Position",seg) ); h1->Fill( hitpos );
        h2 = (TH2F*)gFile->Get( Form("PCu%d_AvT",seg) );  h2->Fill( au, tu );
        h2 = (TH2F*)gFile->Get( Form("PCd%d_AvT",seg) );  h2->Fill( ad, td );
        h2 = (TH2F*)gFile->Get( Form("PCu%d_dEvTime",seg) );  h2->Fill( eneu, timeu );
        h2 = (TH2F*)gFile->Get( Form("PCd%d_dEvTime",seg) );  h2->Fill( ened, timed );
        h2 = (TH2F*)gFile->Get( Form("PCu%d_dEvCTime",seg) );  h2->Fill( eneu, ctimeu );
        h2 = (TH2F*)gFile->Get( Form("PCd%d_dEvCTime",seg) );  h2->Fill( ened, ctimed );
        h2 = (TH2F*)gFile->Get( Form("PC%d_dEMeanvCTimeMean",seg) );  h2->Fill( emean, ctmean );
        h1 = (TH1F*)gFile->Get( "PC_HitPat" ); h1->Fill( seg );
        nPC++;
      }
    }
    h1 = (TH1F*)gFile->Get( "PC_Mult" ); h1->Fill( nPC );

    header->Clear();
    blMan->Clear();
    cdsMan->Clear();
    return true;
}

void EventAnalysisMyPC::Finalize()
{
  std::cout << " Enter EventAnalysisMyPC::Finalize " << std::endl;

  rtFile->cd();
  gFile->Write();
  gFile->Close();

  delete blMan;
  delete cdsMan;
  delete header;
}

void EventAnalysisMyPC::InitializeHistogram()
{
  Int_t NumOfPCSegments = 27;

  rtFile->cd();

  // Scaler
  new TH1F( "Scaler", "Scaler", 50, 0, 50 );

  // Trigger Pattern
  new TH1F( "Pattern", "Trigger Pattern", 20, 0, 20 );

  // PC
  std::cout << "Define Histograms for PC" << std::endl;
  new TH1F( "PC_HitPat", "Hit Pattern PC", NumOfPCSegments+1, 0, NumOfPCSegments+1 );
  for( int seg=1; seg<=NumOfPCSegments; seg++ ){
    new TH1F( Form("PCu%d_ADC",seg),   Form("ADC PCU%d;ADC ch.;Counts",seg),    1000,    0, 4000 );
    new TH1F( Form("PCd%d_ADC",seg),   Form("ADC PCD%d;ADC ch.;Counts",seg),    1000,    0, 4000 );
    new TH1F( Form("PCu%d_dE",seg),   Form("dE PCU%d;dE (MeV);Counts",seg),    400,    0, 40 );
    new TH1F( Form("PCd%d_dE",seg),   Form("dE PCD%d;dE (MeV);Counts",seg),    400,    0, 40 );
    new TH1F( Form("PC%d_dEMean",seg),   Form("Mean dE PC%d;dE (MeV);Counts",seg),    400,    0, 40 );
    new TH1F( Form("PCu%d_TDC",seg),   Form("TDC PCU%d;TDC ch.;Counts",seg),    1000,    0, 4000 );
    new TH1F( Form("PCd%d_TDC",seg),   Form("TDC PCD%d;TDC ch.;Counts",seg),    1000,    0, 4000 );
    new TH1F( Form("PCu%d_Time",seg),   Form("Time PCU%d;Time (ns);Counts",seg),    4000,    -100, 100 );
    new TH1F( Form("PCd%d_Time",seg),   Form("Time PCD%d;Time (ns);Counts",seg),    4000,    -100, 100 );
    new TH1F( Form("PC%d_TimeMean",seg),   Form("Mean Time PC%d;Time (ns);Counts",seg),    4000,    -100, 100 );
    new TH1F( Form("PCu%d_CTime",seg),   Form("CTime PCU%d;Time (ns);Counts",seg),    4000,    -100, 100 );
    new TH1F( Form("PCd%d_CTime",seg),   Form("CTime PCD%d;Time (ns);Counts",seg),    4000,    -100, 100 );
    new TH1F( Form("PC%d_CTimeMean",seg),   Form("Mean CTime PC%d;Time (ns);Counts",seg),    2000,    -100, 100 );
    new TH1F( Form("PC%d_Position",seg),   Form("Hit position PC%d;Position (cm);Counts",seg),    1000,    -50, 50 );
    new TH1F( Form("PCu%d_ADCwT",seg), Form("ADC wTDC PCU%d;ADC ch.;Counts",seg),  1000,    0, 4000 );
    new TH1F( Form("PCd%d_ADCwT",seg), Form("ADC wTDC PCD%d;ADC ch.;Counts",seg),  1000,    0, 4000 );
    new TH2F( Form("PCu%d_AvT",seg),   Form("ADC TDC corr. PCU%d;ADC ch.;TDC ch.",seg),     200,    0, 4000,  200,    0, 4000 );
    new TH2F( Form("PCd%d_AvT",seg),   Form("ADC TDC corr. PCD%d;ADC ch.;TDC ch.",seg),     200,    0, 4000,  200,    0, 4000 );
    new TH2F( Form("PCu%d_dEvTime",seg),   Form("dE Time corr. PCU%d;dE (MeV);Time (ns)",seg),     400,    0, 40,  4000,    -100, 100 );
    new TH2F( Form("PCd%d_dEvTime",seg),   Form("dE Time corr. PCD%d;dE (MeV);Time (ns)",seg),     400,    0, 40,  4000,    -100, 100 );
    new TH2F( Form("PCu%d_dEvCTime",seg),   Form("dE CTime corr. PCU%d;dE (MeV);Time (ns)",seg),     400,    0, 40,  4000,    -100, 100 );
    new TH2F( Form("PCd%d_dEvCTime",seg),   Form("dE CTime corr. PCD%d;dE (MeV);Time (ns)",seg),     400,    0, 40,  4000,    -100, 100 );
    new TH2F( Form("PC%d_dEMeanvCTimeMean",seg),   Form("dEMean CTimeMean corr. PC%d;dE (MeV);Time (ns)",seg),     400,    0, 40,  4000,    -100, 100 );
  }
  new TH1F( "PC_Mult", "Multiplicity PC;Multiplicity;Counts", NumOfPCSegments+1, 0, NumOfPCSegments+1 );

}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysisMyPC *event = new EventAnalysisMyPC();
  return (EventTemp*)event;
}
