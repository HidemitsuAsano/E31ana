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

class EventAnalysisMyNC: public EventTemp
{
public:
  EventAnalysisMyNC();
  ~EventAnalysisMyNC();
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

EventAnalysisMyNC::EventAnalysisMyNC()
  : EventTemp()
{
}

EventAnalysisMyNC::~EventAnalysisMyNC()
{
}

const int MaxTreeSize = 1900000000000;
void EventAnalysisMyNC::Initialize( ConfMan *conf )
{
#if 1
  std::cout << " Enter EventAnalysisMyNC::Initialize " << std::endl;
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

void EventAnalysisMyNC::USca( int nsca, unsigned int *sca )
{
#if 0
  std::cout << " Enter EventAnalysisMyNC::USca " << std::endl;
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

bool EventAnalysisMyNC::UAna( TKOHitCollection *tko )
{
#if 0
	std::cout << " Enter EventAnalysisMyNC::UAna " << std::endl;
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

		// Multihit of NC //
		int MulNC=0;	
		for(int i=0; i<blMan->nNC(); i++){
			HodoscopeLikeHit* hit = blMan->NC(i);
			if(hit->CheckRange()) MulNC++;
		}

		// Selection //
		if(MulNC==0){
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

    // NC
    int nNC=0;
    for( int i=0; i<blMan->nNC(); i++ ){
      HodoscopeLikeHit *hit = blMan->NC(i);
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

      h1 = (TH1F*)gFile->Get( Form("NCu%d_ADC",seg) ); h1->Fill( au );
      h1 = (TH1F*)gFile->Get( Form("NCd%d_ADC",seg) ); h1->Fill( ad );
      if(tu>0){
        h1 = (TH1F*)gFile->Get( Form("NCu%d_TDC",seg) ); h1->Fill( tu );
        h1 = (TH1F*)gFile->Get( Form("NCu%d_Time",seg) ); h1->Fill( timeu );
        h1 = (TH1F*)gFile->Get( Form("NCu%d_CTime",seg) ); h1->Fill( ctimeu );
      }
      if(td>0){
        h1 = (TH1F*)gFile->Get( Form("NCd%d_TDC",seg) ); h1->Fill( td );
        h1 = (TH1F*)gFile->Get( Form("NCd%d_Time",seg) ); h1->Fill( timed );
        h1 = (TH1F*)gFile->Get( Form("NCd%d_CTime",seg) ); h1->Fill( ctimed );
      }
      if( hit->CheckRange() ){
        h1 = (TH1F*)gFile->Get( Form("NCu%d_ADCwT",seg) ); h1->Fill( au );
        h1 = (TH1F*)gFile->Get( Form("NCd%d_ADCwT",seg) ); h1->Fill( ad );
        h1 = (TH1F*)gFile->Get( Form("NCu%d_dE",seg) ); h1->Fill( eneu );
        h1 = (TH1F*)gFile->Get( Form("NCd%d_dE",seg) ); h1->Fill( ened );
        h1 = (TH1F*)gFile->Get( Form("NC%d_dEMean",seg) ); h1->Fill( emean );
        h1 = (TH1F*)gFile->Get( Form("NC%d_TimeMean",seg) ); h1->Fill( tmean );
        h1 = (TH1F*)gFile->Get( Form("NC%d_CTimeMean",seg) ); h1->Fill( ctmean );
        h1 = (TH1F*)gFile->Get( Form("NC%d_Position",seg) ); h1->Fill( hitpos );
        h2 = (TH2F*)gFile->Get( Form("NCu%d_AvT",seg) );  h2->Fill( au, tu );
        h2 = (TH2F*)gFile->Get( Form("NCd%d_AvT",seg) );  h2->Fill( ad, td );
        h2 = (TH2F*)gFile->Get( Form("NCu%d_dEvTime",seg) );  h2->Fill( eneu, timeu );
        h2 = (TH2F*)gFile->Get( Form("NCd%d_dEvTime",seg) );  h2->Fill( ened, timed );
        h2 = (TH2F*)gFile->Get( Form("NCu%d_dEvCTime",seg) );  h2->Fill( eneu, ctimeu );
        h2 = (TH2F*)gFile->Get( Form("NCd%d_dEvCTime",seg) );  h2->Fill( ened, ctimed );
        h2 = (TH2F*)gFile->Get( Form("NC%d_dEMeanvCTimeMean",seg) );  h2->Fill( emean, ctmean );
        h1 = (TH1F*)gFile->Get( "NC_HitPat" ); h1->Fill( seg );
        nNC++;
      }
    }
    h1 = (TH1F*)gFile->Get( "NC_Mult" ); h1->Fill( nNC );

    header->Clear();
    blMan->Clear();
    cdsMan->Clear();
    return true;
}

void EventAnalysisMyNC::Finalize()
{
  std::cout << " Enter EventAnalysisMyNC::Finalize " << std::endl;

  rtFile->cd();
  gFile->Write();
  gFile->Close();

  delete blMan;
  delete cdsMan;
  delete header;
}

void EventAnalysisMyNC::InitializeHistogram()
{
  Int_t NumOfNCSegments = 112;

  rtFile->cd();

  // Scaler
  new TH1F( "Scaler", "Scaler", 50, 0, 50 );

  // Trigger Pattern
  new TH1F( "Pattern", "Trigger Pattern", 20, 0, 20 );

  // NC
  std::cout << "Define Histograms for NC" << std::endl;
  new TH1F( "NC_HitPat", "Hit Pattern NC", NumOfNCSegments+1, 0, NumOfNCSegments+1 );
  for( int seg=1; seg<=NumOfNCSegments; seg++ ){
    new TH1F( Form("NCu%d_ADC",seg),   Form("ADC NCU%d;ADC ch.;Counts",seg),    1000,    0, 4000 );
    new TH1F( Form("NCd%d_ADC",seg),   Form("ADC NCD%d;ADC ch.;Counts",seg),    1000,    0, 4000 );
    new TH1F( Form("NCu%d_dE",seg),   Form("dE NCU%d;dE (MeV);Counts",seg),    400,    0, 40 );
    new TH1F( Form("NCd%d_dE",seg),   Form("dE NCD%d;dE (MeV);Counts",seg),    400,    0, 40 );
    new TH1F( Form("NC%d_dEMean",seg),   Form("Mean dE NC%d;dE (MeV);Counts",seg),    400,    0, 40 );
    new TH1F( Form("NCu%d_TDC",seg),   Form("TDC NCU%d;TDC ch.;Counts",seg),    1000,    0, 4000 );
    new TH1F( Form("NCd%d_TDC",seg),   Form("TDC NCD%d;TDC ch.;Counts",seg),    1000,    0, 4000 );
    new TH1F( Form("NCu%d_Time",seg),   Form("Time NCU%d;Time (ns);Counts",seg),    4000,    -100, 100 );
    new TH1F( Form("NCd%d_Time",seg),   Form("Time NCD%d;Time (ns);Counts",seg),    4000,    -100, 100 );
    new TH1F( Form("NC%d_TimeMean",seg),   Form("Mean Time NC%d;Time (ns);Counts",seg),    4000,    -100, 100 );
    new TH1F( Form("NCu%d_CTime",seg),   Form("CTime NCU%d;Time (ns);Counts",seg),    4000,    -100, 100 );
    new TH1F( Form("NCd%d_CTime",seg),   Form("CTime NCD%d;Time (ns);Counts",seg),    4000,    -100, 100 );
    new TH1F( Form("NC%d_CTimeMean",seg),   Form("Mean CTime NC%d;Time (ns);Counts",seg),    4000,    -100, 100 );
    new TH1F( Form("NC%d_Position",seg),   Form("Hit position NC%d;Position (cm);Counts",seg),    1000,    -50, 50 );
    new TH1F( Form("NCu%d_ADCwT",seg), Form("ADC wTDC NCU%d;ADC ch.;Counts",seg),  1000,    0, 4000 );
    new TH1F( Form("NCd%d_ADCwT",seg), Form("ADC wTDC NCD%d;ADC ch.;Counts",seg),  1000,    0, 4000 );
    new TH2F( Form("NCu%d_AvT",seg),   Form("ADC TDC corr. NCU%d;ADC ch.;TDC ch.",seg),     200,    0, 4000,  200,    0, 4000 );
    new TH2F( Form("NCd%d_AvT",seg),   Form("ADC TDC corr. NCD%d;ADC ch.;TDC ch.",seg),     200,    0, 4000,  200,    0, 4000 );
    new TH2F( Form("NCu%d_dEvTime",seg),   Form("dE Time corr. NCU%d;dE (MeV);Time (ns)",seg),     400,    0, 40,  2000,    -100, 100 );
    new TH2F( Form("NCd%d_dEvTime",seg),   Form("dE Time corr. NCD%d;dE (MeV);Time (ns)",seg),     400,    0, 40,  2000,    -100, 100 );
    new TH2F( Form("NCu%d_dEvCTime",seg),   Form("dE CTime corr. NCU%d;dE (MeV);Time (ns)",seg),     400,    0, 40,  4000,    -100, 100 );
    new TH2F( Form("NCd%d_dEvCTime",seg),   Form("dE CTime corr. NCD%d;dE (MeV);Time (ns)",seg),     400,    0, 40,  4000,    -100, 100 );
    new TH2F( Form("NC%d_dEMeanvCTimeMean",seg),   Form("dEMean CTimeMean corr. NC%d;dE (MeV);Time (ns)",seg),     400,    0, 40,  4000,    -100, 100 );
  }
  new TH1F( "NC_Mult", "Multiplicity NC;Multiplicity;Counts", NumOfNCSegments+1, 0, NumOfNCSegments+1 );

}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysisMyNC *event = new EventAnalysisMyNC();
  return (EventTemp*)event;
}
