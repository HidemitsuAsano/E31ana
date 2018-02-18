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

class EventAnalysisBHDT0: public EventTemp
{
public:
  EventAnalysisBHDT0();
  ~EventAnalysisBHDT0();
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

EventAnalysisBHDT0::EventAnalysisBHDT0()
  : EventTemp()
{
}

EventAnalysisBHDT0::~EventAnalysisBHDT0()
{
}

const int MaxTreeSize = 1900000000000;
void EventAnalysisBHDT0::Initialize( ConfMan *conf )
{
#if 1
  std::cout << " Enter EventAnalysisBHDT0::Initialize " << std::endl;
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

void EventAnalysisBHDT0::USca( int nsca, unsigned int *sca )
{
#if 0
  std::cout << " Enter EventAnalysisBHDT0::USca " << std::endl;
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

bool EventAnalysisBHDT0::UAna( TKOHitCollection *tko )
{
#if 0
  std::cout << " Enter EventAnalysisBHDT0::UAna " << std::endl;
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
  //     header->Clear(); return true;
  //   }
  //   if( !header->cds() ){
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

  // BeamLine
  // BHD
  int nBHD=0;
  for( int i=0; i<blMan->nBHD(); i++ ){
    HodoscopeLikeHit *hit = blMan->BHD(i);
    int seg = hit->seg();
    int au = hit->adcu(), ad = hit->adcd();
    int tu = hit->tdcu(), td = hit->tdcd();
    //std::cout << " seg:" << seg << " au:" << au << " ad:" << ad << " tu:" << tu << " td:" << td << std::endl;
    h1 = (TH1F*)gFile->Get( Form("ABHDU%d",seg) ); h1->Fill( au );
    h1 = (TH1F*)gFile->Get( Form("ABHDD%d",seg) ); h1->Fill( ad );
    h1 = (TH1F*)gFile->Get( Form("TBHDU%d",seg) ); h1->Fill( tu );
    h1 = (TH1F*)gFile->Get( Form("TBHDD%d",seg) ); h1->Fill( td );
    if( hit->CheckRange() ){
      h1 = (TH1F*)gFile->Get( Form("AwTBHDU%d",seg) ); h1->Fill( au );
      h1 = (TH1F*)gFile->Get( Form("AwTBHDD%d",seg) ); h1->Fill( ad );
      h2 = (TH2F*)gFile->Get( Form("ATBHDU%d",seg) );  h2->Fill( au, tu );
      h2 = (TH2F*)gFile->Get( Form("ATBHDD%d",seg) );  h2->Fill( ad, td );
      h1 = (TH1F*)gFile->Get( "HitPatBHD" ); h1->Fill( seg );
      nBHD++;
    }
  }
  h1 = (TH1F*)gFile->Get( "MulBHD" ); h1->Fill( nBHD );
  if( header->kaon() ){
    h1 = (TH1F*)gFile->Get( "MulBHDk" ); h1->Fill( nBHD );
  }
  if( header->pion() ){
    h1 = (TH1F*)gFile->Get( "MulBHDpi" ); h1->Fill( nBHD );
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


  // PA
  int nPA=0;
  for( int i=0; i<blMan->nPA(); i++ ){
    HodoscopeLikeHit *hit = blMan->PA(i);
    if( hit->CheckRange() ){
      nPA++;
    }
  }

  // LC2
  int nLC2=0;
  int asumLC2=0;
  for( int i=0; i<blMan->nLC2(); i++ ){
    CherenkovLikeHit *hit = blMan->LC2(i);
    int seg = hit->seg();
    int a1 = hit->adc1(), a2 = hit->adc2();
    if( 0<a1+a2 && a1+a2<2000 ){
      nLC2++;
      asumLC2 += (a1+a2);
    }
  }

  // BHD
  for( int i=0; i<blMan->nBHD(); i++ ){
    HodoscopeLikeHit *hit = blMan->BHD(i);
    int seg = hit->seg();
    double eu = hit->eu(), ed = hit->ed();
    double tu = hit->tu(), td = hit->td();
    double ctu = hit->ctu(), ctd = hit->ctd();
    double emean = hit->emean();
    double ctmean = hit->ctmean(), ctsub = hit->ctsub();
    h1 = (TH1F*)gFile->Get( Form("eBHDU%d",seg) ); h1->Fill( eu );
    h1 = (TH1F*)gFile->Get( Form("eBHDD%d",seg) ); h1->Fill( ed );
    h1 = (TH1F*)gFile->Get( Form("tBHDU%d",seg) ); h1->Fill( tu );
    h1 = (TH1F*)gFile->Get( Form("tBHDD%d",seg) ); h1->Fill( td );
    h1 = (TH1F*)gFile->Get( Form("ctBHDU%d",seg) ); h1->Fill( ctu );
    h1 = (TH1F*)gFile->Get( Form("ctBHDD%d",seg) ); h1->Fill( ctd );
    if( hit->CheckRange() ){
      h1 = (TH1F*)gFile->Get( Form("ewtBHDU%d",seg) ); h1->Fill( eu );
      h1 = (TH1F*)gFile->Get( Form("ewtBHDD%d",seg) ); h1->Fill( ed );
      h2 = (TH2F*)gFile->Get( Form("etBHDU%d",seg) );  h2->Fill( eu, tu );
      h2 = (TH2F*)gFile->Get( Form("etBHDD%d",seg) );  h2->Fill( ed, td );
      h2 = (TH2F*)gFile->Get( Form("ectBHDU%d",seg) );  h2->Fill( eu, ctu );
      h2 = (TH2F*)gFile->Get( Form("ectBHDD%d",seg) );  h2->Fill( ed, ctd );
      h1 = (TH1F*)gFile->Get( Form("ctmeanBHD%d",seg) ); h1->Fill( ctmean );
      h1 = (TH1F*)gFile->Get( Form("ctsubBHD%d",seg) ); h1->Fill( ctsub );
      h1 = (TH1F*)gFile->Get( Form("emeanBHD%d",seg) ); h1->Fill( emean );
      h2 = (TH2F*)gFile->Get( Form("eeBHD%d",seg) ); h2->Fill( eu, ed );
      h2 = (TH2F*)gFile->Get( Form("euemeanBHD%d",seg) ); h2->Fill( eu, emean );
      h2 = (TH2F*)gFile->Get( Form("edemeanBHD%d",seg) ); h2->Fill( ed, emean );
    }
  }
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

  // BHD and T0
  for( int i1=0; i1<blMan->nBHD(); i1++ ){
    HodoscopeLikeHit *bhd = blMan->BHD(i1);
    if( !bhd->CheckRange() ) continue;
    int seg1 = bhd->seg();
    double eu1 = bhd->eu(), ed1 = bhd->ed();
    double emean1 = bhd->emean();
    double tmean1 = bhd->tmean();
    double ctmean1 = bhd->ctmean();
    for( int i2=0; i2<blMan->nT0(); i2++ ){
      HodoscopeLikeHit *t0 = blMan->T0(i2);
      if( !t0->CheckRange() ) continue;
      int seg2 = t0->seg();
      double eu2 = t0->eu(), ed2 = t0->ed();
      double emean2 = t0->emean();
      double tmean2 = t0->tmean();
      double ctmean2 = t0->ctmean();
      h1 = (TH1F*)gFile->Get( "TOF_BHD_T0" ); h1->Fill( ctmean2-ctmean1 );
      if( nBHD==1 && nT0==1 ){
	h1 = (TH1F*)gFile->Get( "TOF_BHD_T0_2" ); h1->Fill( ctmean2-ctmean1 );
	if( nPA==1 ){
	  h1 = (TH1F*)gFile->Get( "TOF_BHD_T0_3" ); h1->Fill( ctmean2-ctmean1 );
	}
      }
      if( nLC2==0 ){
	h1 = (TH1F*)gFile->Get( "TOF_BHD_T0_4" ); h1->Fill( ctmean2-ctmean1 );
      }
      else if( nLC2==1 ){
	h1 = (TH1F*)gFile->Get( "TOF_BHD_T0_5" ); h1->Fill( ctmean2-ctmean1 );
      }
      else if( nLC2==2 ){
	h1 = (TH1F*)gFile->Get( "TOF_BHD_T0_6" ); h1->Fill( ctmean2-ctmean1 );
      }
      if( asumLC2<340 ){
	h1 = (TH1F*)gFile->Get( "TOF_BHD_T0_7" ); h1->Fill( ctmean2-ctmean1 );
      }
      else{
	h1 = (TH1F*)gFile->Get( "TOF_BHD_T0_8" ); h1->Fill( ctmean2-ctmean1 );
      }

      if( header->electron() ){ h1 = (TH1F*)gFile->Get( "TOF_BHD_T0_E" ); h1->Fill( ctmean2-ctmean1 ); }
      if( header->pion() ){     h1 = (TH1F*)gFile->Get( "TOF_BHD_T0_Pi" ); h1->Fill( ctmean2-ctmean1 );}
      if( header->kaon() ){     h1 = (TH1F*)gFile->Get( "TOF_BHD_T0_K" ); h1->Fill( ctmean2-ctmean1 ); }
      if( header->proton() ){   h1 = (TH1F*)gFile->Get( "TOF_BHD_T0_P" ); h1->Fill( ctmean2-ctmean1 ); }
      h1 = (TH1F*)gFile->Get( Form("TOF_BHD%d_T0%d",seg1,seg2) ); h1->Fill( ctmean2-ctmean1 );
      if( 1.2<emean1 ){
	h1 = (TH1F*)gFile->Get( Form("dETOF1_BHD%d_T0U%d",seg1,seg2) ); h1->Fill( eu2, tmean2-tmean1 );
	h1 = (TH1F*)gFile->Get( Form("dETOF1_BHD%d_T0D%d",seg1,seg2) ); h1->Fill( ed2, tmean2-tmean1 );
	h1 = (TH1F*)gFile->Get( Form("dETOF2_BHD%d_T0U%d",seg1,seg2) ); h1->Fill( eu2, ctmean2-ctmean1 );
	h1 = (TH1F*)gFile->Get( Form("dETOF2_BHD%d_T0D%d",seg1,seg2) ); h1->Fill( ed2, ctmean2-ctmean1 );
      }
      if( 1.5<emean2 ){
	h1 = (TH1F*)gFile->Get( Form("dETOF1_BHDU%d_T0%d",seg1,seg2) ); h1->Fill( eu1, tmean2-tmean1 );
	h1 = (TH1F*)gFile->Get( Form("dETOF1_BHDD%d_T0%d",seg1,seg2) ); h1->Fill( ed1, tmean2-tmean1 );
	h1 = (TH1F*)gFile->Get( Form("dETOF2_BHDU%d_T0%d",seg1,seg2) ); h1->Fill( eu1, ctmean2-ctmean1 );
	h1 = (TH1F*)gFile->Get( Form("dETOF2_BHDD%d_T0%d",seg1,seg2) ); h1->Fill( ed1, ctmean2-ctmean1 );
      }
    }
  }

  header->Clear();
  blMan->Clear();
  cdsMan->Clear();
  return true;
}

void EventAnalysisBHDT0::Finalize()
{
  std::cout << " Enter EventAnalysisBHDT0::Finalize " << std::endl;

  rtFile->cd();
  gFile->Write();
  gFile->Close();

  delete blMan;
  delete cdsMan;
  delete header;
}

void EventAnalysisBHDT0::InitializeHistogram()
{
  rtFile->cd();
  
  // Scaler
  new TH1F( "Scaler", "Scaler", 50, 0, 50 );

  // Trigger Pattern
  new TH1F( "Pattern", "Trigger Pattern", 10, 0, 10 );

  // BHD
  std::cout << "Define Histgram for BHD" << std::endl;
  new TH1F( "MulBHD", "Multiplicity BHD", NumOfBHDSegments+1, 0, NumOfBHDSegments+1 );
  new TH1F( "MulBHDk", "Multiplicity BHD at k trigger", NumOfBHDSegments+1, 0, NumOfBHDSegments+1 );
  new TH1F( "MulBHDpi", "Multiplicity BHD at pi trigger", NumOfBHDSegments+1, 0, NumOfBHDSegments+1 );
  new TH1F( "HitPatBHD", "Hit Pattern BHD", NumOfBHDSegments+1, 0, NumOfBHDSegments+1 );
  for( int seg=1; seg<=NumOfBHDSegments; seg++ ){
    new TH1F( Form("ABHDU%d",seg),   Form("ADC BHDU%d",seg),    1000,    0, 4000 );
    new TH1F( Form("ABHDD%d",seg),   Form("ADC BHDD%d",seg),    1000,    0, 4000 );
    new TH1F( Form("TBHDU%d",seg),   Form("TDC BHDU%d",seg),    1000,    0, 4000 );
    new TH1F( Form("TBHDD%d",seg),   Form("TDC BHDD%d",seg),    1000,    0, 4000 );
    new TH1F( Form("AwTBHDU%d",seg), Form("ADC wTDC BHDU%d",seg),  1000,    0, 4000 );
    new TH1F( Form("AwTBHDD%d",seg), Form("ADC wTDC BHDD%d",seg),  1000,    0, 4000 );
    new TH2F( Form("ATBHDU%d",seg),   Form("ADC TDC corr. BHDU%d",seg),     200,    0, 4000,  200,    0, 4000 );
    new TH2F( Form("ATBHDD%d",seg),   Form("ADC TDC corr. BHDD%d",seg),     200,    0, 4000,  200,    0, 4000 );
  }
  // T0
  std::cout << "Define Histgram for T0" << std::endl;
  new TH1F( "MulT0", "Multiplicity T0", NumOfT0Segments+1, 0, NumOfT0Segments );
  new TH1F( "HitPatT0", "Hit Pattern T0", NumOfT0Segments+1, 0, NumOfT0Segments );
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

  // BHD
  std::cout << "Define Histgram for BHD corrected" << std::endl;
  for( int seg=1; seg<=NumOfBHDSegments; seg++ ){
    new TH1F( Form("eBHDU%d",seg),   Form("dE BHDU%d",seg),    200,  -0.5, 19.5 );
    new TH1F( Form("eBHDD%d",seg),   Form("dE BHDD%d",seg),    200,  -0.5, 19.5 );
    new TH2F( Form("eeBHD%d",seg),   Form("dE corr. BHD%d",seg),200,-0.5, 4.5, 200,-0.5, 4.5 );
    new TH1F( Form("tBHDU%d",seg),   Form("Time BHDU%d",seg),    200,  -5,  35 );
    new TH1F( Form("tBHDD%d",seg),   Form("Time BHDD%d",seg),    200,  -5,  35 );
    new TH1F( Form("ctBHDU%d",seg),   Form("CTime BHDU%d",seg),    500,  -50,  50 );
    new TH1F( Form("ctBHDD%d",seg),   Form("CTime BHDD%d",seg),    500,  -50,  50 );
    new TH1F( Form("ewtBHDU%d",seg), Form("dE wTDC BHDU%d",seg),  200, -0.5, 19.5 );
    new TH1F( Form("ewtBHDD%d",seg), Form("dE wTDC BHDD%d",seg),  200, -0.5, 19.5 );
    new TH2F( Form("etBHDU%d",seg),   Form("dE Time corr. BHDU%d",seg),     200,   -0.5,  4.5, 200,   -5, 5 );
    new TH2F( Form("etBHDD%d",seg),   Form("dE Time corr. BHDD%d",seg),     200,   -0.5,  4.5, 200,   -5, 5 );
    new TH2F( Form("ectBHDU%d",seg),   Form("dE CTime corr. BHDU%d",seg),   200,   -0.5,  4.5, 200,   -5, 5 );
    new TH2F( Form("ectBHDD%d",seg),   Form("dE CTime corr. BHDD%d",seg),   200,   -0.5,  4.5, 200,   -5, 5 );
    new TH1F( Form("ctmeanBHD%d",seg),   Form("MeanCTime BHD%d",seg), 200, -20, 20 );
    new TH1F( Form("ctsubBHD%d",seg),   Form("DiffCTime BHD%d",seg), 200, -20, 20 );
    new TH1F( Form("emeanBHD%d",seg),   Form("MeandE BHD%d",seg), 200, -0.5, 4.5 );
    new TH2F( Form("euemeanBHD%d",seg),   Form("dEU MeandE corr. BHD%d",seg), 200, -0.5, 4.5, 200, -0.5, 4.5 );
    new TH2F( Form("edemeanBHD%d",seg),   Form("dED MeandE corr. BHD%d",seg), 200, -0.5, 4.5, 200, -0.5, 4.5 );
  }
  // T0
  std::cout << "Define Histgram for T0 corrected" << std::endl;
  for( int seg=1; seg<=NumOfT0Segments; seg++ ){
    new TH1F( Form("eT0U%d",seg),   Form("dE T0U%d",seg),    200,  -0.5, 19.5 );
    new TH1F( Form("eT0D%d",seg),   Form("dE T0D%d",seg),    200,  -0.5, 19.5 );
    new TH1F( Form("tT0U%d",seg),   Form("Time T0U%d",seg),    200,  -20,  20 );
    new TH1F( Form("tT0D%d",seg),   Form("Time T0D%d",seg),    200,  -20,  20 );
    new TH1F( Form("ctT0U%d",seg),   Form("CTime T0U%d",seg),    200,  -5,  35 );
    new TH1F( Form("ctT0D%d",seg),   Form("CTime T0D%d",seg),    200,  -5,  35 );
    new TH1F( Form("ewtT0U%d",seg), Form("dE wTDC T0U%d",seg),  200, -0.5, 19.5 );
    new TH1F( Form("ewtT0D%d",seg), Form("dE wTDC T0D%d",seg),  200, -0.5, 19.5 );
    new TH2F( Form("etT0U%d",seg),   Form("dE Time corr. T0U%d",seg),     200,   -0.5,  4.5, 200,   -5, 5 );
    new TH2F( Form("etT0D%d",seg),   Form("dE Time corr. T0D%d",seg),     200,   -0.5,  4.5, 200,   -5, 5 );
    new TH2F( Form("ectT0U%d",seg),   Form("dE CTime corr. T0U%d",seg),   200,   -0.5,  4.5, 200,   -5, 5 );
    new TH2F( Form("ectT0D%d",seg),   Form("dE CTime corr. T0D%d",seg),   200,   -0.5,  4.5, 200,   -5, 5 );
    new TH1F( Form("ctmeanT0%d",seg),   Form("MeanCTime T0%d",seg), 200, -20, 20 );
    new TH1F( Form("ctsubT0%d",seg),   Form("DiffCTime T0%d",seg), 200, -20, 20 );
  }

  // BHD and T0
  std::cout << "Define Histgram for BHD-T0" << std::endl;
  new TH1F( "TOF_BHD_T0","TOF_BHD_T0", 1000, -10, 40 );
  new TH1F( "TOF_BHD_T0_2","TOF_BHD_T0 nBHD=1 and nT0=1", 1000, -10, 40 );
  new TH1F( "TOF_BHD_T0_3","TOF_BHD_T0 nBHD=1 and nT0=1 and nPA=1", 1000, -10, 400 );
  new TH1F( "TOF_BHD_T0_4","TOF_BHD_T0 nLC2==0",  1000, -10, 40 );
  new TH1F( "TOF_BHD_T0_5","TOF_BHD_T0 nLC2==1", 1000, -10, 40 );
  new TH1F( "TOF_BHD_T0_6","TOF_BHD_T0 nLC2==2",  1000, -10, 40 );
  new TH1F( "TOF_BHD_T0_7","TOF_BHD_T0 ALC2SUM<th",  1000, -10, 40 );
  new TH1F( "TOF_BHD_T0_8","TOF_BHD_T0 ALC2SUM>th",  1000, -10, 40 );
  new TH1F( "TOF_BHD_T0_E","TOF_BHD_T0 electron trigger",  1000, -10, 40 );
  new TH1F( "TOF_BHD_T0_Pi","TOF_BHD_T0 pion trigger",  1000, -10, 40 );
  new TH1F( "TOF_BHD_T0_K","TOF_BHD_T0 kaon trigger",  1000, -10, 40 );
  new TH1F( "TOF_BHD_T0_P","TOF_BHD_T0 proton trigger", 1000, -10, 40 );
  for( int seg1=1; seg1<=NumOfBHDSegments; seg1++ ){
    for( int seg2=1; seg2<=NumOfT0Segments; seg2++ ){
      new TH1F( Form("TOF_BHD%d_T0%d",seg1,seg2),Form("TOF_BHD%d_T0%d",seg1,seg2),  1000, -10, 40 );
      new TH2F( Form("dETOF1_BHDU%d_T0%d",seg1,seg2),Form("dE(BHDU) TOF(not corrected) corr. BHD%d-T0%d",seg1,seg2), 200, -0.5, 4.5, 200, -5, 5 );
      new TH2F( Form("dETOF1_BHDD%d_T0%d",seg1,seg2),Form("dE(BHDD) TOF(not corrected) corr. BHD%d-T0%d",seg1,seg2), 200, -0.5, 4.5, 200, -5, 5 );
      new TH2F( Form("dETOF1_BHD%d_T0U%d",seg1,seg2),Form("dE(T0U) TOF(not corrected) corr. BHD%d-T0%d",seg1,seg2),  200, -0.5, 4.5, 200, -5, 5 );
      new TH2F( Form("dETOF1_BHD%d_T0D%d",seg1,seg2),Form("dE(T0D) TOF(not corrected) corr. BHD%d-T0%d",seg1,seg2),  200, -0.5, 4.5, 200, -5, 5 );
      new TH2F( Form("dETOF2_BHDU%d_T0%d",seg1,seg2),Form("dE(BHDU) TOF corr. BHD%d-T0%d",seg1,seg2), 200, -0.5, 4.5, 200, -5, 5 );
      new TH2F( Form("dETOF2_BHDD%d_T0%d",seg1,seg2),Form("dE(BHDD) TOF corr. BHD%d-T0%d",seg1,seg2), 200, -0.5, 4.5, 200, -5, 5 );
      new TH2F( Form("dETOF2_BHD%d_T0U%d",seg1,seg2),Form("dE(T0U) TOF corr. BHD%d-T0%d",seg1,seg2),  200, -0.5, 4.5, 200, -5, 5 );
      new TH2F( Form("dETOF2_BHD%d_T0D%d",seg1,seg2),Form("dE(T0D) TOF corr. BHD%d-T0%d",seg1,seg2),  200, -0.5, 4.5, 200, -5, 5 );
    }
  }
}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysisBHDT0 *event = new EventAnalysisBHDT0();
  return (EventTemp*)event;
}
