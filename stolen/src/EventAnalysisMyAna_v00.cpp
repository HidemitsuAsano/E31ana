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

class EventAnalysisMyAna: public EventTemp
{
public:
  EventAnalysisMyAna();
  ~EventAnalysisMyAna();
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

EventAnalysisMyAna::EventAnalysisMyAna()
  : EventTemp()
{
}

EventAnalysisMyAna::~EventAnalysisMyAna()
{
}

const int MaxTreeSize = 1900000000000;
void EventAnalysisMyAna::Initialize( ConfMan *conf )
{
#if 1
  std::cout << " Enter EventAnalysisMyAna::Initialize " << std::endl;
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

void EventAnalysisMyAna::USca( int nsca, unsigned int *sca )
{
#if 0
  std::cout << " Enter EventAnalysisMyAna::USca " << std::endl;
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

bool EventAnalysisMyAna::UAna( TKOHitCollection *tko )
{
#if 0
	std::cout << " Enter EventAnalysisMyAna::UAna " << std::endl;
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
		cdsMan->Convert( tko, confMan );

		rtFile->cd();
		TH1F *h1;
		TH2F *h2;

		// ======== //
		// Raw Data //
		// ======== //

		// Multihit of CDH and IH //
		int MulCDH=0;	
		for(int i=0; i<cdsMan->nCDH(); i++){
			HodoscopeLikeHit* hit = cdsMan->CDH(i);
			if(hit->CheckRange()) MulCDH++;
		}
		int MulIH=0;
		for(int i=0; i<cdsMan->nIH(); i++){
			HodoscopeLikeHit* hit = cdsMan->IH(i);
			if(hit->CheckRange()) MulIH++;
		}

		// Selection //
		if(!header->IsTrig(Trig_Cosmic) || MulCDH!=2 || MulIH!=2){
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

		// CDH
		int nCDH=0;
		for( int i=0; i<cdsMan->nCDH(); i++ ){
			HodoscopeLikeHit *hit = cdsMan->CDH(i);
			int seg = hit->seg();
			int au = hit->adcu(), ad = hit->adcd();
			int tu = hit->tdcu(), td = hit->tdcd();
#if Debug
			std::cout << "seg = " << seg << std::endl;
#endif

			h1 = (TH1F*)gFile->Get( Form("ACDHU%d",seg) ); h1->Fill( au );
			h1 = (TH1F*)gFile->Get( Form("ACDHD%d",seg) ); h1->Fill( ad );
			h1 = (TH1F*)gFile->Get( Form("TCDHU%d",seg) ); h1->Fill( tu );
			h1 = (TH1F*)gFile->Get( Form("TCDHD%d",seg) ); h1->Fill( td );
			if( hit->CheckRange() ){
				h1 = (TH1F*)gFile->Get( Form("AwTCDHU%d",seg) ); h1->Fill( au );
				h1 = (TH1F*)gFile->Get( Form("AwTCDHD%d",seg) ); h1->Fill( ad );
				h2 = (TH2F*)gFile->Get( Form("ATCDHU%d",seg) );  h2->Fill( au, tu );
				h2 = (TH2F*)gFile->Get( Form("ATCDHD%d",seg) );  h2->Fill( ad, td );
				h1 = (TH1F*)gFile->Get( "HitPatCDH" ); h1->Fill( seg );
				nCDH++;
			}
		}
		h1 = (TH1F*)gFile->Get( "MultCDH" ); h1->Fill( nCDH );

		// IH
		int nIH=0;
		for( int i=0; i<cdsMan->nIH(); i++ ){
			HodoscopeLikeHit *hit = cdsMan->IH(i);
			int seg = hit->seg();
			int au = hit->adcu(), ad = hit->adcd();
			int tu = hit->tdcu(), td = hit->tdcd();
#if Debug
			std::cout << "seg = " << seg << std::endl;
#endif

			h1 = (TH1F*)gFile->Get( Form("AIHU%d",seg) ); h1->Fill( au );
			h1 = (TH1F*)gFile->Get( Form("AIHD%d",seg) ); h1->Fill( ad );
			h1 = (TH1F*)gFile->Get( Form("TIHU%d",seg) ); h1->Fill( tu );
			h1 = (TH1F*)gFile->Get( Form("TIHD%d",seg) ); h1->Fill( td );
			if( hit->CheckRange() ){
				h1 = (TH1F*)gFile->Get( Form("AwTIHU%d",seg) ); h1->Fill( au );
				h1 = (TH1F*)gFile->Get( Form("AwTIHD%d",seg) ); h1->Fill( ad );
				h2 = (TH2F*)gFile->Get( Form("ATIHU%d",seg) );  h2->Fill( au, tu );
				h2 = (TH2F*)gFile->Get( Form("ATIHD%d",seg) );  h2->Fill( ad, td );
				h1 = (TH1F*)gFile->Get( "HitPatIH" ); h1->Fill( seg );
				nIH++;
			}
		}
		h1 = (TH1F*)gFile->Get( "MultIH" ); h1->Fill( nIH );

		header->Clear();
		blMan->Clear();
		cdsMan->Clear();
		return true;
}

void EventAnalysisMyAna::Finalize()
{
	std::cout << " Enter EventAnalysisMyAna::Finalize " << std::endl;

	rtFile->cd();
	gFile->Write();
	gFile->Close();

	delete blMan;
	delete cdsMan;
	delete header;
}

void EventAnalysisMyAna::InitializeHistogram()
{
	Int_t NumOfCDCSegments = 14;
	Int_t NumOfCDHSegments = 36;
	Int_t NumOfIHSegments = 24;

	rtFile->cd();

	// Scaler
	new TH1F( "Scaler", "Scaler", 50, 0, 50 );

	// Trigger Pattern
	new TH1F( "Pattern", "Trigger Pattern", 20, 0, 20 );

	// CDH
	std::cout << "Define Histograms for CDH" << std::endl;
	new TH1F( "HitPatCDH", "Hit Pattern CDH", NumOfCDHSegments+1, 0, NumOfCDHSegments+1 );
	for( int seg=1; seg<=NumOfCDHSegments; seg++ ){
		new TH1F( Form("ACDHU%d",seg),   Form("ADC CDHU%d",seg),    1000,    0, 4000 );
		new TH1F( Form("ACDHD%d",seg),   Form("ADC CDHD%d",seg),    1000,    0, 4000 );
		new TH1F( Form("TCDHU%d",seg),   Form("TDC CDHU%d",seg),    1000,    0, 4000 );
		new TH1F( Form("TCDHD%d",seg),   Form("TDC CDHD%d",seg),    1000,    0, 4000 );
		new TH1F( Form("AwTCDHU%d",seg), Form("ADC wTDC CDHU%d",seg),  1000,    0, 4000 );
		new TH1F( Form("AwTCDHD%d",seg), Form("ADC wTDC CDHD%d",seg),  1000,    0, 4000 );
		new TH2F( Form("ATCDHU%d",seg),   Form("ADC TDC corr. CDHU%d",seg),     200,    0, 4000,  200,    0, 4000 );
		new TH2F( Form("ATCDHD%d",seg),   Form("ADC TDC corr. CDHD%d",seg),     200,    0, 4000,  200,    0, 4000 );
	}
	new TH1F( "MultCDH", "Multipliciy CDH", NumOfCDHSegments+1, 0, NumOfCDHSegments+1 );

	// IH
	std::cout << "Define Histograms for IH" << std::endl;
	new TH1F( "HitPatIH", "Hit Pattern IH", NumOfIHSegments+1, 0, NumOfIHSegments+1 );
	for( int seg=1; seg<=NumOfIHSegments; seg++ ){
		new TH1F( Form("AIHU%d",seg),   Form("ADC IHU%d",seg),    1000,    0, 4000 );
		new TH1F( Form("AIHD%d",seg),   Form("ADC IHD%d",seg),    1000,    0, 4000 );
		new TH1F( Form("TIHU%d",seg),   Form("TDC IHU%d",seg),    1000,    0, 4000 );
		new TH1F( Form("TIHD%d",seg),   Form("TDC IHD%d",seg),    1000,    0, 4000 );
		new TH1F( Form("AwTIHU%d",seg), Form("ADC wTDC IHU%d",seg),  1000,    0, 4000 );
		new TH1F( Form("AwTIHD%d",seg), Form("ADC wTDC IHD%d",seg),  1000,    0, 4000 );
		new TH2F( Form("ATIHU%d",seg),   Form("ADC TDC corr. IHU%d",seg),     200,    0, 4000,  200,    0, 4000 );
		new TH2F( Form("ATIHD%d",seg),   Form("ADC TDC corr. IHD%d",seg),     200,    0, 4000,  200,    0, 4000 );
	}
	new TH1F( "MultIH", "Multipliciy IH", NumOfIHSegments+1, 0, NumOfIHSegments+1 );

}

EventTemp *EventAlloc::EventAllocator()
{
	EventAnalysisMyAna *event = new EventAnalysisMyAna();
	return (EventTemp*)event;
}
