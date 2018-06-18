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

class EventAnalysisT0CDH: public EventTemp
{
public:
  EventAnalysisT0CDH();
  ~EventAnalysisT0CDH();
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

EventAnalysisT0CDH::EventAnalysisT0CDH()
  : EventTemp()
{
}

EventAnalysisT0CDH::~EventAnalysisT0CDH()
{
}

const int MaxTreeSize = 1900000000000;
void EventAnalysisT0CDH::Initialize( ConfMan *conf )
{
#if 1
  std::cout << " Enter EventAnalysisT0CDH::Initialize " << std::endl;
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

void EventAnalysisT0CDH::USca( int nsca, unsigned int *sca )
{
#if 0
  std::cout << " Enter EventAnalysisT0CDH::USca " << std::endl;
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

bool EventAnalysisT0CDH::UAna( TKOHitCollection *tko )
{
#if 0
  std::cout << " Enter EventAnalysisT0CDH::UAna " << std::endl;
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
  // CDH
  int nCDH=0;
  for( int i=0; i<cdsMan->nCDH(); i++ ){
    HodoscopeLikeHit *hit = cdsMan->CDH(i);
    int seg = hit->seg();
    int au = hit->adcu(), ad = hit->adcd();
    int tu = hit->tdcu(), td = hit->tdcd();
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
  h1 = (TH1F*)gFile->Get( "MulCDH" ); h1->Fill( nCDH );

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
  // CDH
  for( int i=0; i<cdsMan->nCDH(); i++ ){
    HodoscopeLikeHit *hit = cdsMan->CDH(i);
    int seg = hit->seg();
    double eu = hit->eu(), ed = hit->ed();
    double tu = hit->tu(), td = hit->td();
    double ctu = hit->ctu(), ctd = hit->ctd();
    double ctmean = hit->ctmean(), ctsub = hit->ctsub();
    h1 = (TH1F*)gFile->Get( Form("eCDHU%d",seg) ); h1->Fill( eu );
    h1 = (TH1F*)gFile->Get( Form("eCDHD%d",seg) ); h1->Fill( ed );
    h1 = (TH1F*)gFile->Get( Form("tCDHU%d",seg) ); h1->Fill( tu );
    h1 = (TH1F*)gFile->Get( Form("tCDHD%d",seg) ); h1->Fill( td );
    h1 = (TH1F*)gFile->Get( Form("ctCDHU%d",seg) ); h1->Fill( ctu );
    h1 = (TH1F*)gFile->Get( Form("ctCDHD%d",seg) ); h1->Fill( ctd );
    if( hit->CheckRange() ){
      h1 = (TH1F*)gFile->Get( Form("ewtCDHU%d",seg) ); h1->Fill( eu );
      h1 = (TH1F*)gFile->Get( Form("ewtCDHD%d",seg) ); h1->Fill( ed );
      h2 = (TH2F*)gFile->Get( Form("etCDHU%d",seg) );  h2->Fill( eu, tu );
      h2 = (TH2F*)gFile->Get( Form("etCDHD%d",seg) );  h2->Fill( ed, td );
      h2 = (TH2F*)gFile->Get( Form("ectCDHU%d",seg) );  h2->Fill( eu, ctu );
      h2 = (TH2F*)gFile->Get( Form("ectCDHD%d",seg) );  h2->Fill( ed, ctd );
      h1 = (TH1F*)gFile->Get( Form("ctmeanCDH%d",seg) ); h1->Fill( ctmean );
      h1 = (TH1F*)gFile->Get( Form("ctsubCDH%d",seg) ); h1->Fill( ctsub );
    }
  }

  // T0 and CDH
  for( int i1=0; i1<blMan->nT0(); i1++ ){
    HodoscopeLikeHit *t0 = blMan->T0(i1);
    if( !t0->CheckRange() ) continue;
    int seg1 = t0->seg();
    double eu1 = t0->eu(), ed1 = t0->ed();
    double ctmean1 = t0->ctmean();
    for( int i2=0; i2<cdsMan->nCDH(); i2++ ){
      HodoscopeLikeHit *cdh = cdsMan->CDH(i2);
      if( !cdh->CheckRange() ) continue;
      int seg2 = cdh->seg();
      double eu2 = cdh->eu(), ed2 = cdh->ed();
      double ctmean2 = cdh->ctmean();
      h1 = (TH1F*)gFile->Get( "TOF_T0_CDH" ); h1->Fill( ctmean2-ctmean1 );
      if( header->electron() ){ h1 = (TH1F*)gFile->Get( "TOF_T0_CDH_E" ); h1->Fill( ctmean2-ctmean1 ); }
      if( header->pion() ){     h1 = (TH1F*)gFile->Get( "TOF_T0_CDH_Pi" ); h1->Fill( ctmean2-ctmean1 );}
      if( header->kaon() ){     h1 = (TH1F*)gFile->Get( "TOF_T0_CDH_K" ); h1->Fill( ctmean2-ctmean1 ); }
      if( header->proton() ){   h1 = (TH1F*)gFile->Get( "TOF_T0_CDH_P" ); h1->Fill( ctmean2-ctmean1 ); }
      h1 = (TH1F*)gFile->Get( Form("TOF_T0%d_CDH%d",seg1,seg2) ); h1->Fill( ctmean2-ctmean1 );
      h1 = (TH1F*)gFile->Get( Form("dECDH_T0U%d_CDH%d",seg1,seg2) ); h1->Fill( eu1, ctmean2-ctmean1 );
      h1 = (TH1F*)gFile->Get( Form("dECDH_T0D%d_CDH%d",seg1,seg2) ); h1->Fill( ed1, ctmean2-ctmean1 );
      h1 = (TH1F*)gFile->Get( Form("dECDH_T0%d_CDHU%d",seg1,seg2) ); h1->Fill( eu2, ctmean2-ctmean1 );
      h1 = (TH1F*)gFile->Get( Form("dECDH_T0%d_CDHD%d",seg1,seg2) ); h1->Fill( ed2, ctmean2-ctmean1 );
    }
  }

  header->Clear();
  blMan->Clear();
  cdsMan->Clear();
  return true;
}

void EventAnalysisT0CDH::Finalize()
{
  std::cout << " Enter EventAnalysisT0CDH::Finalize " << std::endl;

  rtFile->cd();
  gFile->Write();
  gFile->Close();

  delete blMan;
  delete cdsMan;
  delete header;
}

void EventAnalysisT0CDH::InitializeHistogram()
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
  // CDH
  std::cout << "Define Histgram for CDH" << std::endl;
  new TH1F( "MulCDH", "Multiplicity CDH", NumOfCDHSegments+1, 0, NumOfCDHSegments+1 );
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
  // CDH
  std::cout << "Define Histgram for CDH corrected" << std::endl;
  for( int seg=1; seg<=NumOfCDHSegments; seg++ ){
    new TH1F( Form("eCDHU%d",seg),   Form("dE CDHU%d",seg),    200,  -0.5, 19.5 );
    new TH1F( Form("eCDHD%d",seg),   Form("dE CDHD%d",seg),    200,  -0.5, 19.5 );
    new TH1F( Form("tCDHU%d",seg),   Form("Time CDHU%d",seg),    200,  -5,  35 );
    new TH1F( Form("tCDHD%d",seg),   Form("Time CDHD%d",seg),    200,  -5,  35 );
    new TH1F( Form("ctCDHU%d",seg),   Form("CTime CDHU%d",seg),    200,  -5,  35 );
    new TH1F( Form("ctCDHD%d",seg),   Form("CTime CDHD%d",seg),    200,  -5,  35 );
    new TH1F( Form("ewtCDHU%d",seg), Form("dE wTDC CDHU%d",seg),  200, -0.5, 19.5 );
    new TH1F( Form("ewtCDHD%d",seg), Form("dE wTDC CDHD%d",seg),  200, -0.5, 19.5 );
    new TH2F( Form("etCDHU%d",seg),   Form("dE Time corr. CDHU%d",seg),     500, -0.5, 49.5,  200,   -5,  5 );
    new TH2F( Form("etCDHD%d",seg),   Form("dE Time corr. CDHD%d",seg),     500, -0.5, 49.5,  200,   -5,  5 );
    new TH2F( Form("ectCDHU%d",seg),   Form("dE CTime corr. CDHU%d",seg),   500, -0.5, 49.5,  200,   -5,  5 );
    new TH2F( Form("ectCDHD%d",seg),   Form("dE CTime corr. CDHD%d",seg),   500, -0.5, 49.5,  200,   -5,  5 );
    new TH1F( Form("ctmeanCDH%d",seg),   Form("MeanCTime CDH%d",seg), 200, -20, 20 );
    new TH1F( Form("ctsubCDH%d",seg),   Form("DiffCTime CDH%d",seg), 200, -20, 20 );
  }

  // T0 and CDH
  std::cout << "Define Histgram for T0-CDH" << std::endl;
  new TH1F( "TOF_T0_CDH",   "TOF_T0_CDH",                  600, -20, 40 );
  new TH1F( "TOF_T0_CDH_E", "TOF_T0_CDH electron trigger", 600, -20, 40 );
  new TH1F( "TOF_T0_CDH_Pi","TOF_T0_CDH pion trigger",     600, -20, 40 );
  new TH1F( "TOF_T0_CDH_K", "TOF_T0_CDH kaon trigger",     600, -20, 40 );
  new TH1F( "TOF_T0_CDH_P", "TOF_T0_CDH proton trigger",   600, -20, 40 );
  for( int seg1=1; seg1<=NumOfT0Segments; seg1++ ){
    for( int seg2=1; seg2<=NumOfCDHSegments; seg2++ ){
      new TH1F( Form("TOF_T0%d_CDH%d",seg1,seg2),Form("TOF_T0%d_CDH%d",seg1,seg2), 600, -20, 40 );
      new TH2F( Form("dECDH_T0U%d_CDH%d",seg1,seg2),Form("dE(T0U) CDH corr. T0%d-CDH%d",seg1,seg2),  200,  -0.5, 4.5, 200, -5, 5 );
      new TH2F( Form("dECDH_T0D%d_CDH%d",seg1,seg2),Form("dE(T0D) CDH corr. T0%d-CDH%d",seg1,seg2),  200,  -0.5, 4.5, 200, -5, 5 );
      new TH2F( Form("dECDH_T0%d_CDHU%d",seg1,seg2),Form("dE(CDHU) CDH corr. T0%d-CDH%d",seg1,seg2), 500,  -0.5, 49.5, 200, -5, 5 );
      new TH2F( Form("dECDH_T0%d_CDHD%d",seg1,seg2),Form("dE(CDHD) CDH corr. T0%d-CDH%d",seg1,seg2), 500,  -0.5, 49.5, 200, -5, 5 );
    }
  }
}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysisT0CDH *event = new EventAnalysisT0CDH();
  return (EventTemp*)event;
}
