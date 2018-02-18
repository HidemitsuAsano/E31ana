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

class EventAnalysisMyBPD: public EventTemp
{
public:
  EventAnalysisMyBPD();
  ~EventAnalysisMyBPD();
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

EventAnalysisMyBPD::EventAnalysisMyBPD()
  : EventTemp()
{
}

EventAnalysisMyBPD::~EventAnalysisMyBPD()
{
}

const int MaxTreeSize = 1900000000000;
void EventAnalysisMyBPD::Initialize( ConfMan *conf )
{
#if 1
  std::cout << " Enter EventAnalysisMyBPD::Initialize " << std::endl;
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

void EventAnalysisMyBPD::USca( int nsca, unsigned int *sca )
{
#if 0
  std::cout << " Enter EventAnalysisMyBPD::USca " << std::endl;
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

bool EventAnalysisMyBPD::UAna( TKOHitCollection *tko )
{
#if 0
	std::cout << " Enter EventAnalysisMyBPD::UAna " << std::endl;
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

		// Multihit of BPD //
		int MulBPD=0;	
		for(int i=0; i<blMan->nBPD(); i++){
			HodoscopeLikeHit* hit = blMan->BPD(i);
			if(hit->CheckRange()) MulBPD++;
		}

		// Selection //
		if(MulBPD!=2){
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

    // BPD
    int nBPD=0;
    for( int i=0; i<blMan->nBPD(); i++ ){
      HodoscopeLikeHit *hit = blMan->BPD(i);
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
      h1 = (TH1F*)gFile->Get( Form("BPDu%d_ADC",seg) ); h1->Fill( au );
      h1 = (TH1F*)gFile->Get( Form("BPDd%d_ADC",seg) ); h1->Fill( ad );
      if(tu>0){
        h1 = (TH1F*)gFile->Get( Form("BPDu%d_TDC",seg) ); h1->Fill( tu );
        h1 = (TH1F*)gFile->Get( Form("BPDu%d_Time",seg) ); h1->Fill( timeu );
        h1 = (TH1F*)gFile->Get( Form("BPDu%d_CTime",seg) ); h1->Fill( ctimeu );
      }
      if(td>0){
        h1 = (TH1F*)gFile->Get( Form("BPDd%d_TDC",seg) ); h1->Fill( td );
        h1 = (TH1F*)gFile->Get( Form("BPDd%d_Time",seg) ); h1->Fill( timed );
        h1 = (TH1F*)gFile->Get( Form("BPDd%d_CTime",seg) ); h1->Fill( ctimed );
      }
      if( hit->CheckRange() ){
        h1 = (TH1F*)gFile->Get( Form("BPDu%d_ADCwT",seg) ); h1->Fill( au );
        h1 = (TH1F*)gFile->Get( Form("BPDd%d_ADCwT",seg) ); h1->Fill( ad );
        h1 = (TH1F*)gFile->Get( Form("BPDu%d_dE",seg) ); h1->Fill( eneu );
        h1 = (TH1F*)gFile->Get( Form("BPDd%d_dE",seg) ); h1->Fill( ened );
        h1 = (TH1F*)gFile->Get( Form("BPD%d_dEMean",seg) ); h1->Fill( emean );
        h1 = (TH1F*)gFile->Get( Form("BPD%d_TimeMean",seg) ); h1->Fill( tmean );
        h1 = (TH1F*)gFile->Get( Form("BPD%d_CTimeMean",seg) ); h1->Fill( ctmean );
        h1 = (TH1F*)gFile->Get( Form("BPD%d_Position",seg) ); h1->Fill( hitpos );
        h2 = (TH2F*)gFile->Get( Form("BPDu%d_AvT",seg) );  h2->Fill( au, tu );
        h2 = (TH2F*)gFile->Get( Form("BPDd%d_AvT",seg) );  h2->Fill( ad, td );
        h2 = (TH2F*)gFile->Get( Form("BPDu%d_dEvTime",seg) );  h2->Fill( eneu, timeu );
        h2 = (TH2F*)gFile->Get( Form("BPDd%d_dEvTime",seg) );  h2->Fill( ened, timed );
        h2 = (TH2F*)gFile->Get( Form("BPDu%d_dEvCTime",seg) );  h2->Fill( eneu, ctimeu );
        h2 = (TH2F*)gFile->Get( Form("BPDd%d_dEvCTime",seg) );  h2->Fill( ened, ctimed );
        h1 = (TH1F*)gFile->Get( "BPD_HitPat" ); h1->Fill( seg );
        nBPD++;
      }
    }
    h1 = (TH1F*)gFile->Get( "BPD_Mult" ); h1->Fill( nBPD );

    header->Clear();
    blMan->Clear();
    cdsMan->Clear();
    return true;
}

void EventAnalysisMyBPD::Finalize()
{
  std::cout << " Enter EventAnalysisMyBPD::Finalize " << std::endl;

  rtFile->cd();
  gFile->Write();
  gFile->Close();

  delete blMan;
  delete cdsMan;
  delete header;
}

void EventAnalysisMyBPD::InitializeHistogram()
{
  Int_t NumOfBPDSegments = 70;

  rtFile->cd();

  // Scaler
  new TH1F( "Scaler", "Scaler", 50, 0, 50 );

  // Trigger Pattern
  new TH1F( "Pattern", "Trigger Pattern", 20, 0, 20 );

  // BPD
  std::cout << "Define Histograms for BPD" << std::endl;
  new TH1F( "BPD_HitPat", "Hit Pattern BPD", NumOfBPDSegments+1, 0, NumOfBPDSegments+1 );
  for( int seg=1; seg<=NumOfBPDSegments; seg++ ){
    new TH1F( Form("BPDu%d_ADC",seg),   Form("ADC BPDU%d;ADC ch.;Counts",seg),    1000,    0, 4000 );
    new TH1F( Form("BPDd%d_ADC",seg),   Form("ADC BPDD%d;ADC ch.;Counts",seg),    1000,    0, 4000 );
    new TH1F( Form("BPDu%d_dE",seg),   Form("dE BPDU%d;dE (MeV);Counts",seg),    200,    0, 20 );
    new TH1F( Form("BPDd%d_dE",seg),   Form("dE BPDD%d;dE (MeV);Counts",seg),    200,    0, 20 );
    new TH1F( Form("BPD%d_dEMean",seg),   Form("Mean dE BPD%d;dE (MeV);Counts",seg),    200,    0, 20 );
    new TH1F( Form("BPDu%d_TDC",seg),   Form("TDC BPDU%d;TDC ch.;Counts",seg),    1000,    0, 4000 );
    new TH1F( Form("BPDd%d_TDC",seg),   Form("TDC BPDD%d;TDC ch.;Counts",seg),    1000,    0, 4000 );
    new TH1F( Form("BPDu%d_Time",seg),   Form("Time BPDU%d;Time (ns);Counts",seg),    2000,    -100, 100 );
    new TH1F( Form("BPDd%d_Time",seg),   Form("Time BPDD%d;Time (ns);Counts",seg),    2000,    -100, 100 );
    new TH1F( Form("BPD%d_TimeMean",seg),   Form("Mean Time BPD%d;Time (ns);Counts",seg),    2000,    -100, 100 );
    new TH1F( Form("BPDu%d_CTime",seg),   Form("CTime BPDU%d;Time (ns);Counts",seg),    2000,    -100, 100 );
    new TH1F( Form("BPDd%d_CTime",seg),   Form("CTime BPDD%d;Time (ns);Counts",seg),    2000,    -100, 100 );
    new TH1F( Form("BPD%d_CTimeMean",seg),   Form("Mean CTime BPD%d;Time (ns);Counts",seg),    2000,    -100, 100 );
    new TH1F( Form("BPD%d_Position",seg),   Form("Hit position BPD%d;Position (cm);Counts",seg),    1000,    -50, 50 );
    new TH1F( Form("BPDu%d_ADCwT",seg), Form("ADC wTDC BPDU%d;ADC ch.;Counts",seg),  1000,    0, 4000 );
    new TH1F( Form("BPDd%d_ADCwT",seg), Form("ADC wTDC BPDD%d;ADC ch.;Counts",seg),  1000,    0, 4000 );
    new TH2F( Form("BPDu%d_AvT",seg),   Form("ADC TDC corr. BPDU%d;ADC ch.;TDC ch.",seg),     200,    0, 4000,  200,    0, 4000 );
    new TH2F( Form("BPDd%d_AvT",seg),   Form("ADC TDC corr. BPDD%d;ADC ch.;TDC ch.",seg),     200,    0, 4000,  200,    0, 4000 );
    new TH2F( Form("BPDu%d_dEvTime",seg),   Form("dE Time corr. BPDU%d;dE (MeV);Time (ns)",seg),     200,    0, 20,  2000,    -100, 100 );
    new TH2F( Form("BPDd%d_dEvTime",seg),   Form("dE Time corr. BPDD%d;dE (MeV);Time (ns)",seg),     200,    0, 20,  2000,    -100, 100 );
    new TH2F( Form("BPDu%d_dEvCTime",seg),   Form("dE CTime corr. BPDU%d;dE (MeV);Time (ns)",seg),     200,    0, 20,  2000,    -100, 100 );
    new TH2F( Form("BPDd%d_dEvCTime",seg),   Form("dE CTime corr. BPDD%d;dE (MeV);Time (ns)",seg),     200,    0, 20,  2000,    -100, 100 );
  }
  new TH1F( "BPD_Mult", "Multiplicity BPD;Multiplicity;Counts", NumOfBPDSegments+1, 0, NumOfBPDSegments+1 );

}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysisMyBPD *event = new EventAnalysisMyBPD();
  return (EventTemp*)event;
}
