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

class EventAnalysisMyBHD: public EventTemp
{
public:
  EventAnalysisMyBHD();
  ~EventAnalysisMyBHD();
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

EventAnalysisMyBHD::EventAnalysisMyBHD()
  : EventTemp()
{
}

EventAnalysisMyBHD::~EventAnalysisMyBHD()
{
}

const int MaxTreeSize = 1900000000000;
void EventAnalysisMyBHD::Initialize( ConfMan *conf )
{
#if 1
  std::cout << " Enter EventAnalysisMyBHD::Initialize " << std::endl;
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

void EventAnalysisMyBHD::USca( int nsca, unsigned int *sca )
{
#if 0
  std::cout << " Enter EventAnalysisMyBHD::USca " << std::endl;
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

bool EventAnalysisMyBHD::UAna( TKOHitCollection *tko )
{
#if 0
	std::cout << " Enter EventAnalysisMyBHD::UAna " << std::endl;
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
		int MulT0=0;	
    for(int i=0; i<blMan->nT0(); i++){
      HodoscopeLikeHit* hit = blMan->T0(i);
      if(hit->CheckRange()) MulT0++;
    }
		int MulBHD=0;	
    for(int i=0; i<blMan->nBHD(); i++){
      HodoscopeLikeHit* hit = blMan->BHD(i);
      if(hit->CheckRange()) MulBHD++;
    }

    // Selection //
    if(MulT0!=1 || MulBHD!=1 || !header->IsTrig(Trig_Pion)){
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


    // BHD
    int nBHD=0;
    for( int i=0; i<blMan->nBHD(); i++ ){
      HodoscopeLikeHit *hit = blMan->BHD(i);
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

      h1 = (TH1F*)gFile->Get( Form("BHDu%d_ADC",seg) ); h1->Fill( au );
      h1 = (TH1F*)gFile->Get( Form("BHDd%d_ADC",seg) ); h1->Fill( ad );
      if(tu>0){
        h1 = (TH1F*)gFile->Get( Form("BHDu%d_TDC",seg) ); h1->Fill( tu );
        h1 = (TH1F*)gFile->Get( Form("BHDu%d_Time",seg) ); h1->Fill( timeu );
        h1 = (TH1F*)gFile->Get( Form("BHDu%d_CTime",seg) ); h1->Fill( ctimeu );
      }
      if(td>0){
        h1 = (TH1F*)gFile->Get( Form("BHDd%d_TDC",seg) ); h1->Fill( td );
        h1 = (TH1F*)gFile->Get( Form("BHDd%d_Time",seg) ); h1->Fill( timed );
        h1 = (TH1F*)gFile->Get( Form("BHDd%d_CTime",seg) ); h1->Fill( ctimed );
      }
      if( hit->CheckRange() ){
        h1 = (TH1F*)gFile->Get( Form("BHDu%d_ADCwT",seg) ); h1->Fill( au );
        h1 = (TH1F*)gFile->Get( Form("BHDd%d_ADCwT",seg) ); h1->Fill( ad );
        h1 = (TH1F*)gFile->Get( Form("BHDu%d_dE",seg) ); h1->Fill( eneu );
        h1 = (TH1F*)gFile->Get( Form("BHDd%d_dE",seg) ); h1->Fill( ened );
        h1 = (TH1F*)gFile->Get( Form("BHD%d_dEMean",seg) ); h1->Fill( emean );
        h1 = (TH1F*)gFile->Get( Form("BHD%d_TimeMean",seg) ); h1->Fill( tmean );
        h1 = (TH1F*)gFile->Get( Form("BHD%d_CTimeMean",seg) ); h1->Fill( ctmean );
        h1 = (TH1F*)gFile->Get( Form("BHD%d_Position",seg) ); h1->Fill( hitpos );
        //if((eneu<1.5||eneu>2.5)&&(ened<1.5||ened>2.5)){
          h2 = (TH2F*)gFile->Get( Form("BHDu%d_AvT",seg) );  h2->Fill( au, tu );
          h2 = (TH2F*)gFile->Get( Form("BHDd%d_AvT",seg) );  h2->Fill( ad, td );
          h2 = (TH2F*)gFile->Get( Form("BHDu%d_dEvTime",seg) );  h2->Fill( eneu, timeu );
          h2 = (TH2F*)gFile->Get( Form("BHDd%d_dEvTime",seg) );  h2->Fill( ened, timed );
          h2 = (TH2F*)gFile->Get( Form("BHDu%d_dEvCTime",seg) );  h2->Fill( eneu, ctimeu );
          h2 = (TH2F*)gFile->Get( Form("BHDd%d_dEvCTime",seg) );  h2->Fill( ened, ctimed );
          h2 = (TH2F*)gFile->Get( Form("BHDu%d_dEvCTimeMean",seg) );  h2->Fill( eneu, ctmean );
          h2 = (TH2F*)gFile->Get( Form("BHDd%d_dEvCTimeMean",seg) );  h2->Fill( ened, ctmean );
        //}
        h2 = (TH2F*)gFile->Get( Form("BHD%d_dEMeanvCTimeMean",seg) );  h2->Fill( emean, ctmean );
        h1 = (TH1F*)gFile->Get( "BHD_HitPat" ); h1->Fill( seg );
        nBHD++;
      }
    }
    h1 = (TH1F*)gFile->Get( "BHD_Mult" ); h1->Fill( nBHD );

    header->Clear();
    blMan->Clear();
    cdsMan->Clear();
    return true;
}

void EventAnalysisMyBHD::Finalize()
{
  std::cout << " Enter EventAnalysisMyBHD::Finalize " << std::endl;

  rtFile->cd();
  gFile->Write();
  gFile->Close();

  delete blMan;
  delete cdsMan;
  delete header;
}

void EventAnalysisMyBHD::InitializeHistogram()
{
  Int_t NumOfBHDSegments = 20;

  rtFile->cd();

  // Scaler
  new TH1F( "Scaler", "Scaler", 50, 0, 50 );

  // Trigger Pattern
  new TH1F( "Pattern", "Trigger Pattern", 20, 0, 20 );

  // BHD
  std::cout << "Define Histograms for BHD" << std::endl;
  new TH1F( "BHD_HitPat", "Hit Pattern BHD", NumOfBHDSegments+1, 0, NumOfBHDSegments+1 );
  for( int seg=1; seg<=NumOfBHDSegments; seg++ ){
    new TH1F( Form("BHDu%d_ADC",seg),   Form("ADC BHDU%d;ADC ch.;Counts",seg),    4000,    0, 4000 );
    new TH1F( Form("BHDd%d_ADC",seg),   Form("ADC BHDD%d;ADC ch.;Counts",seg),    4000,    0, 4000 );
    new TH1F( Form("BHDu%d_dE",seg),   Form("dE BHDU%d;dE (MeV);Counts",seg),    400,    0, 40 );
    new TH1F( Form("BHDd%d_dE",seg),   Form("dE BHDD%d;dE (MeV);Counts",seg),    400,    0, 40 );
    new TH1F( Form("BHD%d_dEMean",seg),   Form("Mean dE BHD%d;dE (MeV);Counts",seg),    400,    0, 40 );
    new TH1F( Form("BHDu%d_TDC",seg),   Form("TDC BHDU%d;TDC ch.;Counts",seg),    4000,    0, 4000 );
    new TH1F( Form("BHDd%d_TDC",seg),   Form("TDC BHDD%d;TDC ch.;Counts",seg),    4000,    0, 4000 );
    new TH1F( Form("BHDu%d_Time",seg),   Form("Time BHDU%d;Time (ns);Counts",seg),    4000,    -100, 100 );
    new TH1F( Form("BHDd%d_Time",seg),   Form("Time BHDD%d;Time (ns);Counts",seg),    4000,    -100, 100 );
    new TH1F( Form("BHD%d_TimeMean",seg),   Form("Mean Time BHD%d;Time (ns);Counts",seg),    4000,    -100, 100 );
    new TH1F( Form("BHDu%d_CTime",seg),   Form("CTime BHDU%d;Time (ns);Counts",seg),    4000,    -100, 100 );
    new TH1F( Form("BHDd%d_CTime",seg),   Form("CTime BHDD%d;Time (ns);Counts",seg),    4000,    -100, 100 );
    new TH1F( Form("BHD%d_CTimeMean",seg),   Form("Mean CTime BHD%d;Time (ns);Counts",seg),    4000,    -100, 100 );
    new TH1F( Form("BHD%d_Position",seg),   Form("Hit position BHD%d;Position (cm);Counts",seg),    500,    -50, 50 );
    new TH1F( Form("BHDu%d_ADCwT",seg), Form("ADC wTDC BHDU%d;ADC ch.;Counts",seg),  4000,    0, 4000 );
    new TH1F( Form("BHDd%d_ADCwT",seg), Form("ADC wTDC BHDD%d;ADC ch.;Counts",seg),  4000,    0, 4000 );
    new TH2F( Form("BHDu%d_AvT",seg),   Form("ADC TDC corr. BHDU%d;ADC ch.;TDC ch.",seg),     200,    0, 4000,  200,    0, 4000 );
    new TH2F( Form("BHDd%d_AvT",seg),   Form("ADC TDC corr. BHDD%d;ADC ch.;TDC ch.",seg),     200,    0, 4000,  200,    0, 4000 );
    new TH2F( Form("BHDu%d_dEvTime",seg),   Form("dE Time corr. BHDU%d;dE (MeV);Time (ns)",seg),     400,    0, 40,  4000,    -100, 100 );
    new TH2F( Form("BHDd%d_dEvTime",seg),   Form("dE Time corr. BHDD%d;dE (MeV);Time (ns)",seg),     400,    0, 40,  4000,    -100, 100 );
    new TH2F( Form("BHDu%d_dEvCTime",seg),   Form("dE CTime corr. BHDU%d;dE (MeV);Time (ns)",seg),     400,    0, 40,  4000,    -100, 100 );
    new TH2F( Form("BHDd%d_dEvCTime",seg),   Form("dE CTime corr. BHDD%d;dE (MeV);Time (ns)",seg),     400,    0, 40,  4000,    -100, 100 );
    new TH2F( Form("BHDu%d_dEvCTimeMean",seg),   Form("dE CTimeMean corr. BHDU%d;dE (MeV);Time (ns)",seg),     400,    0, 40,  4000,    -100, 100 );
    new TH2F( Form("BHDd%d_dEvCTimeMean",seg),   Form("dE CTimeMean corr. BHDD%d;dE (MeV);Time (ns)",seg),     400,    0, 40,  4000,    -100, 100 );
    new TH2F( Form("BHD%d_dEMeanvCTimeMean",seg),   Form("dEMean CTimeMean corr. BHD%d;dE (MeV);Time (ns)",seg),     400,    0, 40,  4000,    -100, 100 );
  }
  new TH1F( "BHD_Mult", "Multiplicity BHD;Multiplicity;Counts", NumOfBHDSegments+1, 0, NumOfBHDSegments+1 );

}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysisMyBHD *event = new EventAnalysisMyBHD();
  return (EventTemp*)event;
}
