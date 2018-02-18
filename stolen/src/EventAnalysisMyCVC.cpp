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

class EventAnalysisMyCVC: public EventTemp
{
public:
  EventAnalysisMyCVC();
  ~EventAnalysisMyCVC();
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

EventAnalysisMyCVC::EventAnalysisMyCVC()
  : EventTemp()
{
}

EventAnalysisMyCVC::~EventAnalysisMyCVC()
{
}

const int MaxTreeSize = 1900000000000;
void EventAnalysisMyCVC::Initialize( ConfMan *conf )
{
#if 1
  std::cout << " Enter EventAnalysisMyCVC::Initialize " << std::endl;
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

void EventAnalysisMyCVC::USca( int nsca, unsigned int *sca )
{
#if 0
  std::cout << " Enter EventAnalysisMyCVC::USca " << std::endl;
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

bool EventAnalysisMyCVC::UAna( TKOHitCollection *tko )
{
#if 0
	std::cout << " Enter EventAnalysisMyCVC::UAna " << std::endl;
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

		// Multihit of CVC //
		int MulCVC=0;	
		for(int i=0; i<blMan->nCVC(); i++){
			HodoscopeLikeHit* hit = blMan->CVC(i);
			if(hit->CheckRange()) MulCVC++;
		}

		// Selection //
		if(MulCVC==0){
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

    // CVC
    int nCVC=0;
    for( int i=0; i<blMan->nCVC(); i++ ){
      HodoscopeLikeHit *hit = blMan->CVC(i);
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

      h1 = (TH1F*)gFile->Get( Form("CVCu%d_ADC",seg) ); h1->Fill( au );
      h1 = (TH1F*)gFile->Get( Form("CVCd%d_ADC",seg) ); h1->Fill( ad );
      if(tu>0){
        h1 = (TH1F*)gFile->Get( Form("CVCu%d_TDC",seg) ); h1->Fill( tu );
        h1 = (TH1F*)gFile->Get( Form("CVCu%d_Time",seg) ); h1->Fill( timeu );
        h1 = (TH1F*)gFile->Get( Form("CVCu%d_CTime",seg) ); h1->Fill( ctimeu );
      }
      if(td>0){
        h1 = (TH1F*)gFile->Get( Form("CVCd%d_TDC",seg) ); h1->Fill( td );
        h1 = (TH1F*)gFile->Get( Form("CVCd%d_Time",seg) ); h1->Fill( timed );
        h1 = (TH1F*)gFile->Get( Form("CVCd%d_CTime",seg) ); h1->Fill( ctimed );
      }
      if( hit->CheckRange() ){
        h1 = (TH1F*)gFile->Get( Form("CVCu%d_ADCwT",seg) ); h1->Fill( au );
        h1 = (TH1F*)gFile->Get( Form("CVCd%d_ADCwT",seg) ); h1->Fill( ad );
        h1 = (TH1F*)gFile->Get( Form("CVCu%d_dE",seg) ); h1->Fill( eneu );
        h1 = (TH1F*)gFile->Get( Form("CVCd%d_dE",seg) ); h1->Fill( ened );
        h1 = (TH1F*)gFile->Get( Form("CVC%d_dEMean",seg) ); h1->Fill( emean );
        h1 = (TH1F*)gFile->Get( Form("CVC%d_TimeMean",seg) ); h1->Fill( tmean );
        h1 = (TH1F*)gFile->Get( Form("CVC%d_CTimeMean",seg) ); h1->Fill( ctmean );
        h1 = (TH1F*)gFile->Get( Form("CVC%d_Position",seg) ); h1->Fill( hitpos );
        h2 = (TH2F*)gFile->Get( Form("CVCu%d_AvT",seg) );  h2->Fill( au, tu );
        h2 = (TH2F*)gFile->Get( Form("CVCd%d_AvT",seg) );  h2->Fill( ad, td );
        h2 = (TH2F*)gFile->Get( Form("CVCu%d_dEvTime",seg) );  h2->Fill( eneu, timeu );
        h2 = (TH2F*)gFile->Get( Form("CVCd%d_dEvTime",seg) );  h2->Fill( ened, timed );
        h2 = (TH2F*)gFile->Get( Form("CVCu%d_dEvCTime",seg) );  h2->Fill( eneu, ctimeu );
        h2 = (TH2F*)gFile->Get( Form("CVCd%d_dEvCTime",seg) );  h2->Fill( ened, ctimed );
        h2 = (TH2F*)gFile->Get( Form("CVC%d_dEMeanvCTimeMean",seg) );  h2->Fill( emean, ctmean );
        h1 = (TH1F*)gFile->Get( "CVC_HitPat" ); h1->Fill( seg );
        nCVC++;
      }
    }
    h1 = (TH1F*)gFile->Get( "CVC_Mult" ); h1->Fill( nCVC );

    header->Clear();
    blMan->Clear();
    cdsMan->Clear();
    return true;
}

void EventAnalysisMyCVC::Finalize()
{
  std::cout << " Enter EventAnalysisMyCVC::Finalize " << std::endl;

  rtFile->cd();
  gFile->Write();
  gFile->Close();

  delete blMan;
  delete cdsMan;
  delete header;
}

void EventAnalysisMyCVC::InitializeHistogram()
{
  Int_t NumOfCVCSegments = 34;

  rtFile->cd();

  // Scaler
  new TH1F( "Scaler", "Scaler", 50, 0, 50 );

  // Trigger Pattern
  new TH1F( "Pattern", "Trigger Pattern", 20, 0, 20 );

  // CVC
  std::cout << "Define Histograms for CVC" << std::endl;
  new TH1F( "CVC_HitPat", "Hit Pattern CVC", NumOfCVCSegments+1, 0, NumOfCVCSegments+1 );
  for( int seg=1; seg<=NumOfCVCSegments; seg++ ){
    new TH1F( Form("CVCu%d_ADC",seg),   Form("ADC CVCU%d;ADC ch.;Counts",seg),    1000,    0, 4000 );
    new TH1F( Form("CVCd%d_ADC",seg),   Form("ADC CVCD%d;ADC ch.;Counts",seg),    1000,    0, 4000 );
    new TH1F( Form("CVCu%d_dE",seg),   Form("dE CVCU%d;dE (MeV);Counts",seg),    400,    0, 40 );
    new TH1F( Form("CVCd%d_dE",seg),   Form("dE CVCD%d;dE (MeV);Counts",seg),    400,    0, 40 );
    new TH1F( Form("CVC%d_dEMean",seg),   Form("Mean dE CVC%d;dE (MeV);Counts",seg),    400,    0, 40 );
    new TH1F( Form("CVCu%d_TDC",seg),   Form("TDC CVCU%d;TDC ch.;Counts",seg),    1000,    0, 4000 );
    new TH1F( Form("CVCd%d_TDC",seg),   Form("TDC CVCD%d;TDC ch.;Counts",seg),    1000,    0, 4000 );
    new TH1F( Form("CVCu%d_Time",seg),   Form("Time CVCU%d;Time (ns);Counts",seg),    4000,    -100, 100 );
    new TH1F( Form("CVCd%d_Time",seg),   Form("Time CVCD%d;Time (ns);Counts",seg),    4000,    -100, 100 );
    new TH1F( Form("CVC%d_TimeMean",seg),   Form("Mean Time CVC%d;Time (ns);Counts",seg),    4000,    -100, 100 );
    new TH1F( Form("CVCu%d_CTime",seg),   Form("CTime CVCU%d;Time (ns);Counts",seg),    4000,    -100, 100 );
    new TH1F( Form("CVCd%d_CTime",seg),   Form("CTime CVCD%d;Time (ns);Counts",seg),    4000,    -100, 100 );
    new TH1F( Form("CVC%d_CTimeMean",seg),   Form("Mean CTime CVC%d;Time (ns);Counts",seg),    4000,    -100, 100 );
    new TH1F( Form("CVC%d_Position",seg),   Form("Hit position CVC%d;Position (cm);Counts",seg),    1000,    -50, 50 );
    new TH1F( Form("CVCu%d_ADCwT",seg), Form("ADC wTDC CVCU%d;ADC ch.;Counts",seg),  1000,    0, 4000 );
    new TH1F( Form("CVCd%d_ADCwT",seg), Form("ADC wTDC CVCD%d;ADC ch.;Counts",seg),  1000,    0, 4000 );
    new TH2F( Form("CVCu%d_AvT",seg),   Form("ADC TDC corr. CVCU%d;ADC ch.;TDC ch.",seg),     200,    0, 4000,  200,    0, 4000 );
    new TH2F( Form("CVCd%d_AvT",seg),   Form("ADC TDC corr. CVCD%d;ADC ch.;TDC ch.",seg),     200,    0, 4000,  200,    0, 4000 );
    new TH2F( Form("CVCu%d_dEvTime",seg),   Form("dE Time corr. CVCU%d;dE (MeV);Time (ns)",seg),     400,    0, 40,  4000,    -100, 100 );
    new TH2F( Form("CVCd%d_dEvTime",seg),   Form("dE Time corr. CVCD%d;dE (MeV);Time (ns)",seg),     400,    0, 40,  4000,    -100, 100 );
    new TH2F( Form("CVCu%d_dEvCTime",seg),   Form("dE CTime corr. CVCU%d;dE (MeV);Time (ns)",seg),     400,    0, 40,  4000,    -100, 100 );
    new TH2F( Form("CVCd%d_dEvCTime",seg),   Form("dE CTime corr. CVCD%d;dE (MeV);Time (ns)",seg),     400,    0, 40,  4000,    -100, 100 );
    new TH2F( Form("CVC%d_dEMeanvCTimeMean",seg),   Form("dEMean CTimeMean corr. CVC%d;dE (MeV);Time (ns)",seg),     400,    0, 40,  4000,    -100, 100 );
  }
  new TH1F( "CVC_Mult", "Multiplicity CVC;Multiplicity;Counts", NumOfCVCSegments+1, 0, NumOfCVCSegments+1 );

}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysisMyCVC *event = new EventAnalysisMyCVC();
  return (EventTemp*)event;
}
