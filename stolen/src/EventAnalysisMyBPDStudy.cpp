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

const double TOFOffs[5] = {-0.0249496,-0.0770768,0.0312703,-0.0360381,-0.0330606};

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

		// Multihit of TRG and BPD//
    Int_t NumOfTRGSegments = 10;
    int MulTRG[3]={0,0,0};	
    for(int i=0; i<NumOfTRGSegments; i++){
      HodoscopeLikeHit* hit = blMan->BHD(i);
      if(hit->CheckRange()){
        MulTRG[0]++;
        if(hit->seg()<=5)
          MulTRG[1]++;
        else if(hit->seg()<=10)
          MulTRG[2]++;
      }
    }
    Int_t NumOfBPDSegments = 70;
    int MulBPD=0;	
    for(int i=0; i<NumOfBPDSegments; i++){
      HodoscopeLikeHit* hit = blMan->BPD(i);
      if(hit->CheckRange())
        MulBPD++;
    }

    // Selection //
    if(/*MulTRG[0]!=2||*/MulTRG[1]!=1||MulTRG[2]!=1){
      //if(MulTRG[1]!=1||MulTRG[2]!=1){
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


    // TRG
    int nTRG=0;
    int udseg[2] = {0,0};
    double CTimeMean[2] = {-999,-999};
    for( int i=0; i<NumOfTRGSegments; i++ ){
      HodoscopeLikeHit *hit = blMan->BHD(i);
      int seg = hit->seg();
      int au = hit->adcu(), ad = hit->adcd();
      int tu = hit->tdcu(), td = hit->tdcd();
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

      h1 = (TH1F*)gFile->Get( Form("TRGu%d_ADC",seg) ); h1->Fill( au );
      h1 = (TH1F*)gFile->Get( Form("TRGd%d_ADC",seg) ); h1->Fill( ad );
      if(tu>0){
        h1 = (TH1F*)gFile->Get( Form("TRGu%d_TDC",seg) ); h1->Fill( tu );
        h1 = (TH1F*)gFile->Get( Form("TRGu%d_Time",seg) ); h1->Fill( timeu );
        h1 = (TH1F*)gFile->Get( Form("TRGu%d_CTime",seg) ); h1->Fill( ctimeu );
      }
      if(td>0){
        h1 = (TH1F*)gFile->Get( Form("TRGd%d_TDC",seg) ); h1->Fill( td );
        h1 = (TH1F*)gFile->Get( Form("TRGd%d_Time",seg) ); h1->Fill( timed );
        h1 = (TH1F*)gFile->Get( Form("TRGd%d_CTime",seg) ); h1->Fill( ctimed );
      }
      if( hit->CheckRange() ){
        h1 = (TH1F*)gFile->Get( Form("TRGu%d_ADCwT",seg) ); h1->Fill( au );
        h1 = (TH1F*)gFile->Get( Form("TRGd%d_ADCwT",seg) ); h1->Fill( ad );
        h1 = (TH1F*)gFile->Get( Form("TRGu%d_dE",seg) ); h1->Fill( eneu );
        h1 = (TH1F*)gFile->Get( Form("TRGd%d_dE",seg) ); h1->Fill( ened );
        h1 = (TH1F*)gFile->Get( Form("TRG%d_dEMean",seg) ); h1->Fill( emean );
        h1 = (TH1F*)gFile->Get( Form("TRG%d_TimeMean",seg) ); h1->Fill( tmean );
        h1 = (TH1F*)gFile->Get( Form("TRG%d_CTimeMean",seg) ); h1->Fill( ctmean );
        h2 = (TH2F*)gFile->Get( Form("TRGu%d_AvT",seg) );  h2->Fill( au, tu );
        h2 = (TH2F*)gFile->Get( Form("TRGd%d_AvT",seg) );  h2->Fill( ad, td );
        h2 = (TH2F*)gFile->Get( Form("TRGu%d_dEvTime",seg) );  h2->Fill( eneu, timeu );
        h2 = (TH2F*)gFile->Get( Form("TRGd%d_dEvTime",seg) );  h2->Fill( ened, timed );
        h2 = (TH2F*)gFile->Get( Form("TRGu%d_dEvCTime",seg) );  h2->Fill( eneu, ctimeu );
        h2 = (TH2F*)gFile->Get( Form("TRGd%d_dEvCTime",seg) );  h2->Fill( ened, ctimed );
        h1 = (TH1F*)gFile->Get( "TRG_HitPat" ); h1->Fill( seg );
        if(seg <=5){
          udseg[0] = seg;
          CTimeMean[0] = ctmean;
        }
        else if(seg<=10){
          udseg[1] = seg-5;
          CTimeMean[1] = ctmean;
        }
        nTRG++;
      }
    }
    h1 = (TH1F*)gFile->Get( "TRG_Mult" ); h1->Fill( nTRG );
    h2 = (TH2F*)gFile->Get( "TRG_HitPat2D" ); h2->Fill( udseg[0],udseg[1] );
    // TOF
    double trgtof = CTimeMean[1] - CTimeMean[0];
    h1 = (TH1F*)gFile->Get( Form("TRG%dU_TRG%dD_TOF",udseg[0],udseg[1]) ); h1->Fill(trgtof);
    h1 = (TH1F*)gFile->Get( "TRGU_TRGD_TOF" ); h1->Fill(trgtof);


    // BPD
    int nBPD=0;
    double BPD_CTimeMean;
    for( int i=0; i<NumOfBPDSegments; i++ ){
      HodoscopeLikeHit *hit = blMan->BPD(i);
      int seg = hit->seg();
      int au = hit->adcu(), ad = hit->adcd();
      int tu = hit->tdcu(), td = hit->tdcd();
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
        h2 = (TH2F*)gFile->Get( Form("BPDu%d_AvT",seg) );  h2->Fill( au, tu );
        h2 = (TH2F*)gFile->Get( Form("BPDd%d_AvT",seg) );  h2->Fill( ad, td );
        h2 = (TH2F*)gFile->Get( Form("BPDu%d_dEvTime",seg) );  h2->Fill( eneu, timeu );
        h2 = (TH2F*)gFile->Get( Form("BPDd%d_dEvTime",seg) );  h2->Fill( ened, timed );
        h2 = (TH2F*)gFile->Get( Form("BPDu%d_dEvCTime",seg) );  h2->Fill( eneu, ctimeu );
        h2 = (TH2F*)gFile->Get( Form("BPDd%d_dEvCTime",seg) );  h2->Fill( ened, ctimed );
        h1 = (TH1F*)gFile->Get( "BPD_HitPat" ); h1->Fill( seg );
        BPD_CTimeMean = ctmean;
        nBPD++;

        // TOF
        double bpdtof_u = CTimeMean[0] - BPD_CTimeMean;
        double bpdtof_d = CTimeMean[1] - BPD_CTimeMean;
        if(udseg[0]==udseg[1] && udseg[0]!=0){
          h1 = (TH1F*)gFile->Get( Form("TRG%dU_BPD_TOF",udseg[0]) ); h1->Fill(bpdtof_u);
          h1 = (TH1F*)gFile->Get( Form("TRG%dD_BPD_TOF",udseg[0]) ); h1->Fill(bpdtof_d);
          h1 = (TH1F*)gFile->Get( Form("TRG%dU_BPD%02d_TOF",udseg[0],seg) ); h1->Fill(bpdtof_u);
          h1 = (TH1F*)gFile->Get( Form("TRG%dD_BPD%02d_TOF",udseg[0],seg) ); h1->Fill(bpdtof_d);
          if(seg>=25 && seg<=45){
            h1 = (TH1F*)gFile->Get( "TRGU_BPD_TOF" ); h1->Fill(bpdtof_u);
            h1 = (TH1F*)gFile->Get( "TRGD_BPD_TOF" ); h1->Fill(bpdtof_d);
          }
          h1 = (TH1F*)gFile->Get( Form("TRGU_BPD%02d_TOF",seg) ); h1->Fill(bpdtof_u);
          h1 = (TH1F*)gFile->Get( Form("TRGD_BPD%02d_TOF",seg) ); h1->Fill(bpdtof_d);
        }
      }
    }
    h1 = (TH1F*)gFile->Get( "BPD_Mult" ); h1->Fill( nBPD );


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
      Int_t NumOfTRGSegments = 10;
      Int_t NumOfBPDSegments = 70;

      rtFile->cd();

      // Scaler
      new TH1F( "Scaler", "Scaler", 50, 0, 50 );

      // Trigger Pattern
      new TH1F( "Pattern", "Trigger Pattern", 20, 0, 20 );

      // TRG
      std::cout << "Define Histograms for TRG" << std::endl;
      for( int seg=1; seg<=NumOfTRGSegments; seg++ ){
        new TH1F( Form("TRGu%d_ADC",seg),   Form("ADC TRGU%d;ADC ch.;Counts",seg),    1000,    0, 4000 );
        new TH1F( Form("TRGd%d_ADC",seg),   Form("ADC TRGD%d;ADC ch.;Counts",seg),    1000,    0, 4000 );
        new TH1F( Form("TRGu%d_dE",seg),   Form("dE TRGU%d;dE (MeV);Counts",seg),    200,    0, 20 );
        new TH1F( Form("TRGd%d_dE",seg),   Form("dE TRGD%d;dE (MeV);Counts",seg),    200,    0, 20 );
        new TH1F( Form("TRG%d_dEMean",seg),   Form("Mean dE TRG%d;dE (MeV);Counts",seg),    200,    0, 20 );
        new TH1F( Form("TRGu%d_TDC",seg),   Form("TDC TRGU%d;TDC ch.;Counts",seg),    1000,    0, 4000 );
        new TH1F( Form("TRGd%d_TDC",seg),   Form("TDC TRGD%d;TDC ch.;Counts",seg),    1000,    0, 4000 );
        new TH1F( Form("TRGu%d_Time",seg),   Form("Time TRGU%d;Time (ns);Counts",seg),    4000,    -20, 20 );
        new TH1F( Form("TRGd%d_Time",seg),   Form("Time TRGD%d;Time (ns);Counts",seg),    4000,    -20, 20 );
        new TH1F( Form("TRG%d_TimeMean",seg),   Form("Mean Time TRG%d;Time (ns);Counts",seg),    4000,    -20, 20 );
        new TH1F( Form("TRGu%d_CTime",seg),   Form("CTime TRGU%d;Time (ns);Counts",seg),    4000,    -20, 20 );
        new TH1F( Form("TRGd%d_CTime",seg),   Form("CTime TRGD%d;Time (ns);Counts",seg),    4000,    -20, 20 );
        new TH1F( Form("TRG%d_CTimeMean",seg),   Form("Mean CTime TRG%d;Time (ns);Counts",seg),    4000,    -20, 20 );
        new TH1F( Form("TRGu%d_ADCwT",seg), Form("ADC wTDC TRGU%d;ADC ch.;Counts",seg),  1000,    0, 4000 );
        new TH1F( Form("TRGd%d_ADCwT",seg), Form("ADC wTDC TRGD%d;ADC ch.;Counts",seg),  1000,    0, 4000 );
        new TH2F( Form("TRGu%d_AvT",seg),   Form("ADC TDC corr. TRGU%d;ADC ch.;TDC ch.",seg),     200,    0, 4000,  200,    0, 4000 );
        new TH2F( Form("TRGd%d_AvT",seg),   Form("ADC TDC corr. TRGD%d;ADC ch.;TDC ch.",seg),     200,    0, 4000,  200,    0, 4000 );
        new TH2F( Form("TRGu%d_dEvTime",seg),   Form("dE Time corr. TRGU%d;dE (MeV);Time (ns)",seg),     200,    0, 20,  2000,    -20, 20 );
        new TH2F( Form("TRGd%d_dEvTime",seg),   Form("dE Time corr. TRGD%d;dE (MeV);Time (ns)",seg),     200,    0, 20,  2000,    -20, 20 );
        new TH2F( Form("TRGu%d_dEvCTime",seg),   Form("dE CTime corr. TRGU%d;dE (MeV);Time (ns)",seg),     200,    0, 20,  2000,    -100, 100 );
        new TH2F( Form("TRGd%d_dEvCTime",seg),   Form("dE CTime corr. TRGD%d;dE (MeV);Time (ns)",seg),     200,    0, 20,  2000,    -100, 100 );
      }
      new TH1F( "TRG_Mult", "Multiplicity TRG;Multiplicity;Counts", NumOfTRGSegments, 1, NumOfTRGSegments+1 );
      new TH1F( "TRG_HitPat", "Hit Pattern TRG;Segment;Counts", NumOfTRGSegments, 1, NumOfTRGSegments+1 );
      new TH2F( "TRG_HitPat2D", "2D Hit Pattern TRG;Segment (up);Segment (down)", NumOfTRGSegments/2, 1, NumOfTRGSegments/2+1, NumOfTRGSegments/2, 1, NumOfTRGSegments/2+1);

      // BPD
      std::cout << "Define Histograms for BPD" << std::endl;
      for( int seg=1; seg<=NumOfBPDSegments; seg++ ){
        new TH1F( Form("BPDu%d_ADC",seg),   Form("ADC BPDU%d;ADC ch.;Counts",seg),    1000,    0, 4000 );
        new TH1F( Form("BPDd%d_ADC",seg),   Form("ADC BPDD%d;ADC ch.;Counts",seg),    1000,    0, 4000 );
        new TH1F( Form("BPDu%d_dE",seg),   Form("dE BPDU%d;dE (MeV);Counts",seg),    200,    0, 20 );
        new TH1F( Form("BPDd%d_dE",seg),   Form("dE BPDD%d;dE (MeV);Counts",seg),    200,    0, 20 );
        new TH1F( Form("BPD%d_dEMean",seg),   Form("Mean dE BPD%d;dE (MeV);Counts",seg),    200,    0, 20 );
        new TH1F( Form("BPDu%d_TDC",seg),   Form("TDC BPDU%d;TDC ch.;Counts",seg),    1000,    0, 4000 );
        new TH1F( Form("BPDd%d_TDC",seg),   Form("TDC BPDD%d;TDC ch.;Counts",seg),    1000,    0, 4000 );
        new TH1F( Form("BPDu%d_Time",seg),   Form("Time BPDU%d;Time (ns);Counts",seg),    4000,    -20, 20 );
        new TH1F( Form("BPDd%d_Time",seg),   Form("Time BPDD%d;Time (ns);Counts",seg),    4000,    -20, 20 );
        new TH1F( Form("BPD%d_TimeMean",seg),   Form("Mean Time BPD%d;Time (ns);Counts",seg),    4000,    -20, 20 );
        new TH1F( Form("BPDu%d_CTime",seg),   Form("CTime BPDU%d;Time (ns);Counts",seg),    4000,    -20, 20 );
        new TH1F( Form("BPDd%d_CTime",seg),   Form("CTime BPDD%d;Time (ns);Counts",seg),    4000,    -20, 20 );
        new TH1F( Form("BPD%d_CTimeMean",seg),   Form("Mean CTime BPD%d;Time (ns);Counts",seg),    4000,    -20, 20 );
        new TH1F( Form("BPDu%d_ADCwT",seg), Form("ADC wTDC BPDU%d;ADC ch.;Counts",seg),  1000,    0, 4000 );
        new TH1F( Form("BPDd%d_ADCwT",seg), Form("ADC wTDC BPDD%d;ADC ch.;Counts",seg),  1000,    0, 4000 );
        new TH2F( Form("BPDu%d_AvT",seg),   Form("ADC TDC corr. BPDU%d;ADC ch.;TDC ch.",seg),     200,    0, 4000,  200,    0, 4000 );
        new TH2F( Form("BPDd%d_AvT",seg),   Form("ADC TDC corr. BPDD%d;ADC ch.;TDC ch.",seg),     200,    0, 4000,  200,    0, 4000 );
        new TH2F( Form("BPDu%d_dEvTime",seg),   Form("dE Time corr. BPDU%d;dE (MeV);Time (ns)",seg),     200,    0, 20,  2000,    -20, 20 );
        new TH2F( Form("BPDd%d_dEvTime",seg),   Form("dE Time corr. BPDD%d;dE (MeV);Time (ns)",seg),     200,    0, 20,  2000,    -20, 20 );
        new TH2F( Form("BPDu%d_dEvCTime",seg),   Form("dE CTime corr. BPDU%d;dE (MeV);Time (ns)",seg),     200,    0, 20,  2000,    -100, 100 );
        new TH2F( Form("BPDd%d_dEvCTime",seg),   Form("dE CTime corr. BPDD%d;dE (MeV);Time (ns)",seg),     200,    0, 20,  2000,    -100, 100 );
      }
      new TH1F( "BPD_Mult", "Multiplicity BPD;Multiplicity;Counts", NumOfBPDSegments, 1, NumOfBPDSegments+1 );
      new TH1F( "BPD_HitPat", "Hit Pattern BPD;Segment;Counts", NumOfBPDSegments, 1, NumOfBPDSegments+1 );

      // TOF
      new TH1F( "TRGU_TRGD_TOF", "TOF TRGU and TRGD;TOF (ns);Counts", 400, -5, 5 );
      new TH1F( "TRGU_BPD_TOF", "TOF TRGU and BPD;TOF (ns);Counts", 400, -5, 5 );
      new TH1F( "TRGD_BPD_TOF", "TOF TRGD and BPD;TOF (ns);Counts", 400, -5, 5 );
      for(int j=1; j<=NumOfBPDSegments; j++){
        new TH1F( Form("TRGU_BPD%02d_TOF",j), Form("TOF TRGU and BPD%02d;TOF (ns);Counts",j), 400, -5, 5 );
        new TH1F( Form("TRGD_BPD%02d_TOF",j), Form("TOF TRGD and BPD%02d;TOF (ns);Counts",j), 400, -5, 5 );
      }
      for(int i=1; i<=NumOfTRGSegments/2; i++){
        for(int j=1; j<=NumOfTRGSegments/2; j++){
          new TH1F( Form("TRG%dU_TRG%dD_TOF",i,j), Form("TOF TRG%dU and TRG%dD;TOF (ns);Counts",i,j), 400, -5, 5 );
        }
        new TH1F( Form("TRG%dU_BPD_TOF",i), Form("TOF TRG%dU and BPD;TOF (ns);Counts",i), 400, -5, 5 );
        new TH1F( Form("TRG%dD_BPD_TOF",i), Form("TOF TRG%dD and BPD;TOF (ns);Counts",i), 400, -5, 5 );
        for(int k=1; k<=NumOfBPDSegments; k++){
          new TH1F( Form("TRG%dU_BPD%02d_TOF",i,k), Form("TOF TRG%dU and BPD%02d;TOF (ns);Counts",i,k), 400, -5, 5 );
          new TH1F( Form("TRG%dD_BPD%02d_TOF",i,k), Form("TOF TRG%dD and BPD%02d;TOF (ns);Counts",i,k), 400, -5, 5 );
        }
      }
    }

    EventTemp *EventAlloc::EventAllocator()
    {
      EventAnalysisMyBHD *event = new EventAnalysisMyBHD();
      return (EventTemp*)event;
    }
