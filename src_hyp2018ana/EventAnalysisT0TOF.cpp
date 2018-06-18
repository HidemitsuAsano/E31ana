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

class EventAnalysisRawAll: public EventTemp
{
public:
  EventAnalysisRawAll();
  ~EventAnalysisRawAll();
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

EventAnalysisRawAll::EventAnalysisRawAll()
  : EventTemp()
{
}

EventAnalysisRawAll::~EventAnalysisRawAll()
{
}

const int MaxTreeSize = 1900000000000;
void EventAnalysisRawAll::Initialize( ConfMan *conf )
{
#if 1
  std::cout << " Enter EventAnalysisRawAll::Initialize " << std::endl;
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

void EventAnalysisRawAll::USca( int nsca, unsigned int *sca )
{
#if 0
  std::cout << " Enter EventAnalysisRawAll::USca " << std::endl;
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

bool EventAnalysisRawAll::UAna( TKOHitCollection *tko )
{
#if 0
  std::cout << " Enter EventAnalysisRawAll::UAna " << std::endl;
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
//   if( !header->electron() ){
//     //if( !header->electron() && !header->pion() ){
//     header->Clear(); return true;
//   }
  blMan->Convert( tko, confMan );
  cdsMan->Convert( tko, confMan );


  rtFile->cd();
  TH1F *h1;
  TH2F *h2;

  // ======== //
  // Raw Data //
  // ======== //

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
    h1 = (TH1F*)gFile->Get( Form("AT0U%d",seg) ); h1->Fill( au );
    h1 = (TH1F*)gFile->Get( Form("AT0D%d",seg) ); h1->Fill( ad );
    h1 = (TH1F*)gFile->Get( Form("TT0U%d",seg) ); h1->Fill( tu );
    h1 = (TH1F*)gFile->Get( Form("TT0D%d",seg) ); h1->Fill( td );
    if( hit->CheckRange() ){
      h1 = (TH1F*)gFile->Get( Form("AwTT0U%d",seg) ); h1->Fill( au );
      h1 = (TH1F*)gFile->Get( Form("AwTT0D%d",seg) ); h1->Fill( ad );
      h2 = (TH2F*)gFile->Get( Form("ATT0U%d",seg) );  h2->Fill( au, tu );
      h2 = (TH2F*)gFile->Get( Form("ATT0D%d",seg) );  h2->Fill( ad, td );
      h1 = (TH1F*)gFile->Get( "HitPatT0" ); h1->Fill( seg );
      nT0++;
    }
  }
  h1 = (TH1F*)gFile->Get( "MulT0" ); h1->Fill( nT0 );
  // TOF
  int nTOF=0;
  for( int i=0; i<blMan->nTOF(); i++ ){
    HodoscopeLikeHit *hit = blMan->TOF(i);
    int seg = hit->seg();
    int au = hit->adcu(), ad = hit->adcd();
    int tu = hit->tdcu(), td = hit->tdcd();
    h1 = (TH1F*)gFile->Get( Form("ATOFU%d",seg) ); h1->Fill( au );
    h1 = (TH1F*)gFile->Get( Form("ATOFD%d",seg) ); h1->Fill( ad );
    h1 = (TH1F*)gFile->Get( Form("TTOFU%d",seg) ); h1->Fill( tu );
    h1 = (TH1F*)gFile->Get( Form("TTOFD%d",seg) ); h1->Fill( td );
    if( hit->CheckRange() ){
      h1 = (TH1F*)gFile->Get( Form("AwTTOFU%d",seg) ); h1->Fill( au );
      h1 = (TH1F*)gFile->Get( Form("AwTTOFD%d",seg) ); h1->Fill( ad );
      h2 = (TH2F*)gFile->Get( Form("ATTOFU%d",seg) );  h2->Fill( au, tu );
      h2 = (TH2F*)gFile->Get( Form("ATTOFD%d",seg) );  h2->Fill( ad, td );
      h1 = (TH1F*)gFile->Get( "HitPatTOF" ); h1->Fill( seg );
      nTOF++;
    }
  }
  h1 = (TH1F*)gFile->Get( "MulTOF" ); h1->Fill( nTOF );

  //           //
  // Corrected //
  //           //
  // T0
  for( int i=0; i<blMan->nT0(); i++ ){
    HodoscopeLikeHit *hit = blMan->T0(i);
    int seg = hit->seg();
    double eu = hit->eu(), ed = hit->ed();
    double tu = hit->tu(), td = hit->td();
    double ctu = hit->ctu(), ctd = hit->ctd();
    double ctmean = hit->ctmean(), ctsub = hit->ctsub();
    h1 = (TH1F*)gFile->Get( Form("eT0U%d",seg) ); h1->Fill( eu );
    h1 = (TH1F*)gFile->Get( Form("eT0D%d",seg) ); h1->Fill( ed );
    h1 = (TH1F*)gFile->Get( Form("tT0U%d",seg) ); h1->Fill( tu );
    h1 = (TH1F*)gFile->Get( Form("tT0D%d",seg) ); h1->Fill( td );
    h1 = (TH1F*)gFile->Get( Form("ctT0U%d",seg) ); h1->Fill( ctu );
    h1 = (TH1F*)gFile->Get( Form("ctT0D%d",seg) ); h1->Fill( ctd );
    if( hit->CheckRange() ){
      h1 = (TH1F*)gFile->Get( Form("ewtT0U%d",seg) ); h1->Fill( eu );
      h1 = (TH1F*)gFile->Get( Form("ewtT0D%d",seg) ); h1->Fill( ed );
      h2 = (TH2F*)gFile->Get( Form("etT0U%d",seg) );  h2->Fill( eu, tu );
      h2 = (TH2F*)gFile->Get( Form("etT0D%d",seg) );  h2->Fill( ed, td );
      h2 = (TH2F*)gFile->Get( Form("ectT0U%d",seg) );  h2->Fill( eu, ctu );
      h2 = (TH2F*)gFile->Get( Form("ectT0D%d",seg) );  h2->Fill( ed, ctd );
      h1 = (TH1F*)gFile->Get( Form("ctmeanT0%d",seg) ); h1->Fill( ctmean );
      h1 = (TH1F*)gFile->Get( Form("ctsubT0%d",seg) ); h1->Fill( ctsub );
    }
  }
  // TOF
  for( int i=0; i<blMan->nTOF(); i++ ){
    HodoscopeLikeHit *hit = blMan->TOF(i);
    int seg = hit->seg();
    double eu = hit->eu(), ed = hit->ed();
    double tu = hit->tu(), td = hit->td();
    double ctu = hit->ctu(), ctd = hit->ctd();
    double ctmean = hit->ctmean(), ctsub = hit->ctsub();
    h1 = (TH1F*)gFile->Get( Form("eTOFU%d",seg) ); h1->Fill( eu );
    h1 = (TH1F*)gFile->Get( Form("eTOFD%d",seg) ); h1->Fill( ed );
    h1 = (TH1F*)gFile->Get( Form("tTOFU%d",seg) ); h1->Fill( tu );
    h1 = (TH1F*)gFile->Get( Form("tTOFD%d",seg) ); h1->Fill( td );
    h1 = (TH1F*)gFile->Get( Form("ctTOFU%d",seg) ); h1->Fill( ctu );
    h1 = (TH1F*)gFile->Get( Form("ctTOFD%d",seg) ); h1->Fill( ctd );
    if( hit->CheckRange() ){
      h1 = (TH1F*)gFile->Get( Form("ewtTOFU%d",seg) ); h1->Fill( eu );
      h1 = (TH1F*)gFile->Get( Form("ewtTOFD%d",seg) ); h1->Fill( ed );
      h2 = (TH2F*)gFile->Get( Form("etTOFU%d",seg) );  h2->Fill( eu, tu );
      h2 = (TH2F*)gFile->Get( Form("etTOFD%d",seg) );  h2->Fill( ed, td );
      h2 = (TH2F*)gFile->Get( Form("ectTOFU%d",seg) );  h2->Fill( eu, ctu );
      h2 = (TH2F*)gFile->Get( Form("ectTOFD%d",seg) );  h2->Fill( ed, ctd );
      h1 = (TH1F*)gFile->Get( Form("ctmeanTOF%d",seg) ); h1->Fill( ctmean );
      h1 = (TH1F*)gFile->Get( Form("ctsubTOF%d",seg) ); h1->Fill( ctsub );
    }
  }

  // T0 and TOF
  for( int i1=0; i1<blMan->nT0(); i1++ ){
    HodoscopeLikeHit *t0 = blMan->T0(i1);
    if( !t0->CheckRange() ) continue;
    int seg1 = t0->seg();
    double eu1 = t0->eu(), ed1 = t0->ed();
    double ctmean1 = t0->ctmean();
    for( int i2=0; i2<blMan->nTOF(); i2++ ){
      HodoscopeLikeHit *tof = blMan->TOF(i2);
      if( !tof->CheckRange() ) continue;
      int seg2 = tof->seg();
      double eu2 = tof->eu(), ed2 = tof->ed();
      double ctmean2 = tof->ctmean();
      h1 = (TH1F*)gFile->Get( "TOF_T0_TOF" ); h1->Fill( ctmean2-ctmean1 );
      if( header->electron() ){ h1 = (TH1F*)gFile->Get( "TOF_T0_TOF_E" ); h1->Fill( ctmean2-ctmean1 ); }
      if( header->pion() ){     h1 = (TH1F*)gFile->Get( "TOF_T0_TOF_Pi" ); h1->Fill( ctmean2-ctmean1 );}
      if( header->kaon() ){     h1 = (TH1F*)gFile->Get( "TOF_T0_TOF_K" ); h1->Fill( ctmean2-ctmean1 ); }
      if( header->proton() ){   h1 = (TH1F*)gFile->Get( "TOF_T0_TOF_P" ); h1->Fill( ctmean2-ctmean1 ); }
      h1 = (TH1F*)gFile->Get( Form("TOF_T0%d_TOF%d",seg1,seg2) ); h1->Fill( ctmean2-ctmean1 );
      h1 = (TH1F*)gFile->Get( Form("dETOF_T0U%d_TOF%d",seg1,seg2) ); h1->Fill( eu1, ctmean2-ctmean1 );
      h1 = (TH1F*)gFile->Get( Form("dETOF_T0D%d_TOF%d",seg1,seg2) ); h1->Fill( ed1, ctmean2-ctmean1 );
      h1 = (TH1F*)gFile->Get( Form("dETOF_T0%d_TOFU%d",seg1,seg2) ); h1->Fill( eu2, ctmean2-ctmean1 );
      h1 = (TH1F*)gFile->Get( Form("dETOF_T0%d_TOFD%d",seg1,seg2) ); h1->Fill( ed2, ctmean2-ctmean1 );
    }
  }

  header->Clear();
  blMan->Clear();
  cdsMan->Clear();
  return true;
}

void EventAnalysisRawAll::Finalize()
{
  std::cout << " Enter EventAnalysisRawAll::Finalize " << std::endl;

  rtFile->cd();
  gFile->Write();
  gFile->Close();

  delete blMan;
  delete cdsMan;
  delete header;
}

void EventAnalysisRawAll::InitializeHistogram()
{
  rtFile->cd();
  
  // Scaler
  new TH1F( "Scaler", "Scaler", 50, 0, 50 );

  // Trigger Pattern
  new TH1F( "Pattern", "Trigger Pattern", 10, 0, 10 );

  // T0
  std::cout << "Define Histgram for T0" << std::endl;
  new TH1F( "MulT0", "Multiplicity T0", NumOfT0Segments+1, 0, NumOfT0Segments+1 );
  new TH1F( "HitPatT0", "Hit Pattern T0", NumOfT0Segments+1, 0, NumOfT0Segments+1 );
  for( int seg=1; seg<=NumOfT0Segments; seg++ ){
    new TH1F( Form("AT0U%d",seg),   Form("ADC T0U%d",seg),    1000,    0, 4000 );
    new TH1F( Form("AT0D%d",seg),   Form("ADC T0D%d",seg),    1000,    0, 4000 );
    new TH1F( Form("TT0U%d",seg),   Form("TDC T0U%d",seg),    1000,    0, 4000 );
    new TH1F( Form("TT0D%d",seg),   Form("TDC T0D%d",seg),    1000,    0, 4000 );
    new TH1F( Form("AwTT0U%d",seg), Form("ADC wTDC T0U%d",seg),  1000,    0, 4000 );
    new TH1F( Form("AwTT0D%d",seg), Form("ADC wTDC T0D%d",seg),  1000,    0, 4000 );
    new TH2F( Form("ATT0U%d",seg),   Form("ADC TDC corr. T0U%d",seg),     200,    0, 4000,  200,    0, 4000 );
    new TH2F( Form("ATT0D%d",seg),   Form("ADC TDC corr. T0D%d",seg),     200,    0, 4000,  200,    0, 4000 );
  }
  // TOF
  std::cout << "Define Histgram for TOF" << std::endl;
  new TH1F( "MulTOF", "Multiplicity TOF", NumOfTOFstopSegments+1, 0, NumOfTOFstopSegments+1 );
  new TH1F( "HitPatTOF", "Hit Pattern TOF", NumOfTOFstopSegments+1, 0, NumOfTOFstopSegments+1 );
  for( int seg=1; seg<=NumOfTOFstopSegments; seg++ ){
    new TH1F( Form("ATOFU%d",seg),   Form("ADC TOFU%d",seg),    1000,    0, 4000 );
    new TH1F( Form("ATOFD%d",seg),   Form("ADC TOFD%d",seg),    1000,    0, 4000 );
    new TH1F( Form("TTOFU%d",seg),   Form("TDC TOFU%d",seg),    1000,    0, 4000 );
    new TH1F( Form("TTOFD%d",seg),   Form("TDC TOFD%d",seg),    1000,    0, 4000 );
    new TH1F( Form("AwTTOFU%d",seg), Form("ADC wTDC TOFU%d",seg),  1000,    0, 4000 );
    new TH1F( Form("AwTTOFD%d",seg), Form("ADC wTDC TOFD%d",seg),  1000,    0, 4000 );
    new TH2F( Form("ATTOFU%d",seg),   Form("ADC TDC corr. TOFU%d",seg),     200,    0, 4000,  200,    0, 4000 );
    new TH2F( Form("ATTOFD%d",seg),   Form("ADC TDC corr. TOFD%d",seg),     200,    0, 4000,  200,    0, 4000 );
  }

  // Corrected
  // T0
  std::cout << "Define Histgram for T0 corrected" << std::endl;
  for( int seg=1; seg<=NumOfT0Segments; seg++ ){
    new TH1F( Form("eT0U%d",seg),   Form("dE T0U%d",seg),    200,  -0.5, 19.5 );
    new TH1F( Form("eT0D%d",seg),   Form("dE T0D%d",seg),    200,  -0.5, 19.5 );
    new TH1F( Form("tT0U%d",seg),   Form("Time T0U%d",seg),    200,  -5,  35 );
    new TH1F( Form("tT0D%d",seg),   Form("Time T0D%d",seg),    200,  -5,  35 );
    new TH1F( Form("ctT0U%d",seg),   Form("CTime T0U%d",seg),    200,  -5,  35 );
    new TH1F( Form("ctT0D%d",seg),   Form("CTime T0D%d",seg),    200,  -5,  35 );
    new TH1F( Form("ewtT0U%d",seg), Form("dE wTDC T0U%d",seg),  200, -0.5, 19.5 );
    new TH1F( Form("ewtT0D%d",seg), Form("dE wTDC T0D%d",seg),  200, -0.5, 19.5 );
    new TH2F( Form("etT0U%d",seg),   Form("dE Time corr. T0U%d",seg),     200,   -0.5, 4.5,  200,   -5,  5 );
    new TH2F( Form("etT0D%d",seg),   Form("dE Time corr. T0D%d",seg),     200,   -0.5, 4.5,  200,   -5,  5 );
    new TH2F( Form("ectT0U%d",seg),   Form("dE CTime corr. T0U%d",seg),     200, -0.5, 4.5,  200,   -5,  5 );
    new TH2F( Form("ectT0D%d",seg),   Form("dE CTime corr. T0D%d",seg),     200, -0.5, 4.5,  200,   -5,  5 );
    new TH1F( Form("ctmeanT0%d",seg),   Form("MeanCTime T0%d",seg), 200, -20, 20 );
    new TH1F( Form("ctsubT0%d",seg),   Form("DiffCTime T0%d",seg), 200, -20, 20 );
  }
  // TOF
  std::cout << "Define Histgram for TOF corrected" << std::endl;
  for( int seg=1; seg<=NumOfTOFstopSegments; seg++ ){
    new TH1F( Form("eTOFU%d",seg),   Form("dE TOFU%d",seg),    200,  -0.5, 19.5 );
    new TH1F( Form("eTOFD%d",seg),   Form("dE TOFD%d",seg),    200,  -0.5, 19.5 );
    new TH1F( Form("tTOFU%d",seg),   Form("Time TOFU%d",seg),    200,  -5,  35 );
    new TH1F( Form("tTOFD%d",seg),   Form("Time TOFD%d",seg),    200,  -5,  35 );
    new TH1F( Form("ctTOFU%d",seg),   Form("CTime TOFU%d",seg),    200,  -5,  35 );
    new TH1F( Form("ctTOFD%d",seg),   Form("CTime TOFD%d",seg),    200,  -5,  35 );
    new TH1F( Form("ewtTOFU%d",seg), Form("dE wTDC TOFU%d",seg),  200, -0.5, 19.5 );
    new TH1F( Form("ewtTOFD%d",seg), Form("dE wTDC TOFD%d",seg),  200, -0.5, 19.5 );
    new TH2F( Form("etTOFU%d",seg),   Form("dE Time corr. TOFU%d",seg),     500, -0.5, 49.5,  200,   -5,  5 );
    new TH2F( Form("etTOFD%d",seg),   Form("dE Time corr. TOFD%d",seg),     500, -0.5, 49.5,  200,   -5,  5 );
    new TH2F( Form("ectTOFU%d",seg),   Form("dE CTime corr. TOFU%d",seg),   500, -0.5, 49.5,  200,   -5,  5 );
    new TH2F( Form("ectTOFD%d",seg),   Form("dE CTime corr. TOFD%d",seg),   500, -0.5, 49.5,  200,   -5,  5 );
    new TH1F( Form("ctmeanTOF%d",seg),   Form("MeanCTime TOF%d",seg), 200, -20, 20 );
    new TH1F( Form("ctsubTOF%d",seg),   Form("DiffCTime TOF%d",seg), 200, -20, 20 );
  }

  // T0 and TOF
  std::cout << "Define Histgram for T0-TOF" << std::endl;
  new TH1F( "TOF_T0_TOF",   "TOF_T0_TOF",                  600, -20, 40 );
  new TH1F( "TOF_T0_TOF_E", "TOF_T0_TOF electron trigger", 600, -20, 40 );
  new TH1F( "TOF_T0_TOF_Pi","TOF_T0_TOF pion trigger",     600, -20, 40 );
  new TH1F( "TOF_T0_TOF_K", "TOF_T0_TOF kaon trigger",     600, -20, 40 );
  new TH1F( "TOF_T0_TOF_P", "TOF_T0_TOF proton trigger",   600, -20, 40 );
  for( int seg1=1; seg1<=NumOfT0Segments; seg1++ ){
    for( int seg2=1; seg2<=NumOfTOFstopSegments; seg2++ ){
      new TH1F( Form("TOF_T0%d_TOF%d",seg1,seg2),Form("TOF_T0%d_TOF%d",seg1,seg2), 600, -20, 40 );
      new TH2F( Form("dETOF_T0U%d_TOF%d",seg1,seg2),Form("dE(T0U) TOF corr. T0%d-TOF%d",seg1,seg2),  200,  -0.5, 4.5, 200, -5, 5 );
      new TH2F( Form("dETOF_T0D%d_TOF%d",seg1,seg2),Form("dE(T0D) TOF corr. T0%d-TOF%d",seg1,seg2),  200,  -0.5, 4.5, 200, -5, 5 );
      new TH2F( Form("dETOF_T0%d_TOFU%d",seg1,seg2),Form("dE(TOFU) TOF corr. T0%d-TOF%d",seg1,seg2), 500,  -0.5, 49.5, 200, -5, 5 );
      new TH2F( Form("dETOF_T0%d_TOFD%d",seg1,seg2),Form("dE(TOFD) TOF corr. T0%d-TOF%d",seg1,seg2), 500,  -0.5, 49.5, 200, -5, 5 );
    }
  }
}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysisRawAll *event = new EventAnalysisRawAll();
  return (EventTemp*)event;
}
