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

		rtFile->cd();
		TH1F *h1;
		TH2F *h2;

		// ======== //
		// Raw Data //
		// ======== //

		// Multihit of CDH and IH //
		int MulBPD=0;	
		for(int i=0; i<blMan->nBPD(); i++){
			HodoscopeLikeHit* hit = blMan->BPD(i);
			if(hit->CheckRange()) MulBPD++;
		}
		int MulTRG=0;
		for(int i=0; i<blMan->nBHD(); i++){
			HodoscopeLikeHit* hit = blMan->BHD(i);
			if(hit->CheckRange()) MulTRG++;
		}

//		// Selection //
//		if(!header->IsTrig(Trig_Cosmic) || MulCDH!=2 || MulIH!=2){
//			header->Clear();
//			blMan->Clear();
//			cdsMan->Clear();
//			return true;
//		}

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
#if Debug
			std::cout << "= bpd =" << std::endl;
			std::cout << "seg = " << seg << std::endl;
#endif

			h1 = (TH1F*)gFile->Get( Form("ABPDU%d",seg) ); h1->Fill( au );
			h1 = (TH1F*)gFile->Get( Form("ABPDD%d",seg) ); h1->Fill( ad );
      if(tu>0){
        h1 = (TH1F*)gFile->Get( Form("TBPDU%d",seg) ); h1->Fill( tu );
      }
      if(td>0){
        h1 = (TH1F*)gFile->Get( Form("TBPDD%d",seg) ); h1->Fill( td );
      }
      if( hit->CheckRange() ){
        h1 = (TH1F*)gFile->Get( Form("AwTBPDU%d",seg) ); h1->Fill( au );
        h1 = (TH1F*)gFile->Get( Form("AwTBPDD%d",seg) ); h1->Fill( ad );
        h2 = (TH2F*)gFile->Get( Form("ATBPDU%d",seg) );  h2->Fill( au, tu );
        h2 = (TH2F*)gFile->Get( Form("ATBPDD%d",seg) );  h2->Fill( ad, td );
        h1 = (TH1F*)gFile->Get( "HitPatBPD" ); h1->Fill( seg );
        nBPD++;
      }
    }
    h1 = (TH1F*)gFile->Get( "MultBPD" ); h1->Fill( nBPD );

    // TRG
    int nTRG=0;
    for( int i=0; i<blMan->nBHD(); i++ ){
      HodoscopeLikeHit *hit = blMan->BHD(i);
      int seg = hit->seg();
      int au = hit->adcu(), ad = hit->adcd();
      int tu = hit->tdcu(), td = hit->tdcd();
#if Debug
      std::cout << "= trg =" << std::endl;
      std::cout << "seg = " << seg << std::endl;
#endif

      h1 = (TH1F*)gFile->Get( Form("ATRGU%d",seg) ); h1->Fill( au );
      h1 = (TH1F*)gFile->Get( Form("ATRGD%d",seg) ); h1->Fill( ad );
      if(tu>0){
        h1 = (TH1F*)gFile->Get( Form("TTRGU%d",seg) ); h1->Fill( tu );
      }
      if(td>0){
        h1 = (TH1F*)gFile->Get( Form("TTRGD%d",seg) ); h1->Fill( td );
      }
      if( hit->CheckRange() ){
        h1 = (TH1F*)gFile->Get( Form("AwTTRGU%d",seg) ); h1->Fill( au );
        h1 = (TH1F*)gFile->Get( Form("AwTTRGD%d",seg) ); h1->Fill( ad );
        h2 = (TH2F*)gFile->Get( Form("ATTRGU%d",seg) );  h2->Fill( au, tu );
        h2 = (TH2F*)gFile->Get( Form("ATTRGD%d",seg) );  h2->Fill( ad, td );
        h1 = (TH1F*)gFile->Get( "HitPatTRG" ); h1->Fill( seg );
        nTRG++;
      }
    }
    h1 = (TH1F*)gFile->Get( "MultTRG" ); h1->Fill( nTRG );

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
  Int_t NumOfTRGSegments = 20;
  Int_t NumOfBPDSegments = 70;

  rtFile->cd();

  // Scaler
  new TH1F( "Scaler", "Scaler", 50, 0, 50 );

  // Trigger Pattern
  new TH1F( "Pattern", "Trigger Pattern", 20, 0, 20 );

  // TRG
  std::cout << "Define Histograms for TRG" << std::endl;
  new TH1F( "HitPatTRG", "Hit Pattern TRG", NumOfTRGSegments+1, 0, NumOfTRGSegments+1 );
  for( int seg=1; seg<=NumOfTRGSegments; seg++ ){
    new TH1F( Form("ATRGU%d",seg),   Form("ADC TRGU%d",seg),    1000,    0, 4000 );
    new TH1F( Form("ATRGD%d",seg),   Form("ADC TRGD%d",seg),    1000,    0, 4000 );
    new TH1F( Form("TTRGU%d",seg),   Form("TDC TRGU%d",seg),    1000,    0, 4000 );
    new TH1F( Form("TTRGD%d",seg),   Form("TDC TRGD%d",seg),    1000,    0, 4000 );
    new TH1F( Form("AwTTRGU%d",seg), Form("ADC wTDC TRGU%d",seg),  1000,    0, 4000 );
    new TH1F( Form("AwTTRGD%d",seg), Form("ADC wTDC TRGD%d",seg),  1000,    0, 4000 );
    new TH2F( Form("ATTRGU%d",seg),   Form("ADC TDC corr. TRGU%d",seg),     200,    0, 4000,  200,    0, 4000 );
    new TH2F( Form("ATTRGD%d",seg),   Form("ADC TDC corr. TRGD%d",seg),     200,    0, 4000,  200,    0, 4000 );
  }
  new TH1F( "MultTRG", "Multipliciy TRG", NumOfTRGSegments+1, 0, NumOfTRGSegments+1 );

  // BPD
  std::cout << "Define Histograms for BPD" << std::endl;
  new TH1F( "HitPatBPD", "Hit Pattern BPD", NumOfBPDSegments+1, 0, NumOfBPDSegments+1 );
  for( int seg=1; seg<=NumOfBPDSegments; seg++ ){
    new TH1F( Form("ABPDU%d",seg),   Form("ADC BPDU%d",seg),    1000,    0, 4000 );
    new TH1F( Form("ABPDD%d",seg),   Form("ADC BPDD%d",seg),    1000,    0, 4000 );
    new TH1F( Form("TBPDU%d",seg),   Form("TDC BPDU%d",seg),    1000,    0, 4000 );
    new TH1F( Form("TBPDD%d",seg),   Form("TDC BPDD%d",seg),    1000,    0, 4000 );
    new TH1F( Form("AwTBPDU%d",seg), Form("ADC wTDC BPDU%d",seg),  1000,    0, 4000 );
    new TH1F( Form("AwTBPDD%d",seg), Form("ADC wTDC BPDD%d",seg),  1000,    0, 4000 );
    new TH2F( Form("ATBPDU%d",seg),   Form("ADC TDC corr. BPDU%d",seg),     200,    0, 4000,  200,    0, 4000 );
    new TH2F( Form("ATBPDD%d",seg),   Form("ADC TDC corr. BPDD%d",seg),     200,    0, 4000,  200,    0, 4000 );
  }
  new TH1F( "MultBPD", "Multipliciy BPD", NumOfBPDSegments+1, 0, NumOfBPDSegments+1 );


}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysisMyAna *event = new EventAnalysisMyAna();
  return (EventTemp*)event;
}
