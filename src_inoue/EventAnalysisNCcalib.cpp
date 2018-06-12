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

#define TKOHIS 0
class EventAnalysisRawAll: public EventTemp
{
public:
  EventAnalysisRawAll();
  ~EventAnalysisRawAll();
private:
  TFile *rtFile;
  TTree *evTree;
  TTree *scaTree;
  //  CDSHitMan *cdsMan;
  BeamLineHitMan *blMan;
  EventHeader *header;
  ScalerMan *scaMan;
  int Layer[2];
  int sseg[2][4];
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
  Layer[0]=1;
  Layer[1]=2;
  for(int i=0;i<4;i++)
    {
      sseg[0][i]=i+1;
      sseg[1][i]=i+1;
    }
}

EventAnalysisRawAll::~EventAnalysisRawAll()
{
}

const int MaxTreeSize = 1900000000000;
void EventAnalysisRawAll::Initialize( ConfMan *conf )
{
#if 0
  std::cout << " Enter EventAnalysisRawAll::Initialize " << std::endl;
#endif
  confMan = conf;
  
  rtFile = new TFile( confMan->GetOutFileName().c_str(), "recreate" );
  InitializeHistogram();
  evTree = new TTree( "EventTree", "EventTree" );
  scaTree= new TTree( "ScalerTree", "ScalerTree" );
  header = new EventHeader();
  if( header==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  evTree->Branch( "EventHeader", &header );
  //  cdsMan = new CDSHitMan();
  //  if( cdsMan==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  //  evTree->Branch( "CDSHitMan", &cdsMan );
  blMan = new BeamLineHitMan();
  if( blMan ==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  evTree->Branch( "BeamLineHitMan", &blMan );
  scaMan = new ScalerMan();
  if( scaMan==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  scaTree->Branch( "ScalerMan", &scaMan );

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

#if 0
  for( int i=0; i<scaMan->nsca(); i++ ){
    int val = scaMan->sca(i)->val();
    TString name = scaMan->sca(i)->name();
    std::cout << std::setw(14) << name << std::setw(10) << val;
    if( (i+1)%8==0 ) std::cout << std::endl;
  }
  std::cout << std::endl;
#endif

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

  TH1F *h1;
  TH2F *h2;

  //TKO
#if TKOHIS  
  for(int ntko=0;ntko<tko->entries();ntko++)
    {
      TKOHit *tkohit=tko->hit(ntko);
      h1 = (TH1F*)gFile->Get(Form("TKOc%ds%d",tkohit->cr(),tkohit->sl() )); h1->Fill(tkohit->ch());
      h1 = (TH1F*)gFile->Get(Form("TKOc%ds%da%d",tkohit->cr(),tkohit->sl(),tkohit->ch() ) ); h1->Fill( tkohit->data() );
      
    }
#endif 
//#########

  header->Convert( tko, confMan );
  blMan->Convert( tko, confMan );
  //  cdsMan->Convert( tko, confMan );


  rtFile->cd();

  // ======== //
  // Raw Data //
  // ======== //
  // NC
  int nNC=0;
  for( int i=0; i<blMan->nNC(); i++ ){
    HodoscopeLikeHit *hit = blMan->NC(i);
    int seg = hit->seg();
    //    std::cout<<"NC seg "<<seg<<std::endl;
    int au = hit->adcu(), ad = hit->adcd();
    double eu = hit->eu(), ed = hit->ed();
    int tu = hit->tdcu(), td = hit->tdcd();
    double ctm = hit->ctmean();
    double ctsub = hit->ctsub();
    h1 = (TH1F*)gFile->Get( Form("ANCU%d",seg) ); h1->Fill( au );
    h1 = (TH1F*)gFile->Get( Form("ANCD%d",seg) ); h1->Fill( ad );
    h1 = (TH1F*)gFile->Get( Form("ENCU%d",seg) ); h1->Fill( eu );
    h1 = (TH1F*)gFile->Get( Form("ENCD%d",seg) ); h1->Fill( ed );
    h1 = (TH1F*)gFile->Get( Form("TNCU%d",seg) ); h1->Fill( tu );
    h1 = (TH1F*)gFile->Get( Form("TNCD%d",seg) ); h1->Fill( td );
    if( !hit->CheckRange() ){
      h1 = (TH1F*)gFile->Get( Form("AwoTNCU%d",seg) ); h1->Fill( au );
      h1 = (TH1F*)gFile->Get( Form("AwoTNCD%d",seg) ); h1->Fill( ad );
    }else{
      h1 = (TH1F*)gFile->Get( Form("AwTNCU%d",seg) ); h1->Fill( au );
      h1 = (TH1F*)gFile->Get( Form("AwTNCD%d",seg) ); h1->Fill( ad );
      h1 = (TH1F*)gFile->Get( Form("TsubNC%d",seg) );  h1->Fill( ctsub );
      h1 = (TH1F*)gFile->Get( Form("TmeanNC%d",seg) );  h1->Fill( ctm );
      h2 = (TH2F*)gFile->Get( Form("ATNCU%d",seg) );  h2->Fill( tu, au );
      h2 = (TH2F*)gFile->Get( Form("ATNCD%d",seg) );  h2->Fill( td, ad );
      h2 = (TH2F*)gFile->Get( Form("AudNC%d",seg) );  h2->Fill( au, ad );
      h2 = (TH2F*)gFile->Get( Form("TudNC%d",seg) );  h2->Fill( tu, td );
      //      h2 = (TH2F*)gFile->Get( Form("EvTNC%d",seg) );  h2->Fill( Event_Number, (tu-td) );
      h1 = (TH1F*)gFile->Get( "HitPatNC" ); h1->Fill( seg );
      nNC++;
    }    
  }
  h1 = (TH1F*)gFile->Get( "MulNC" ); h1->Fill( nNC );

  // PC
  int nPC=0;
  for( int i=0; i<blMan->nPC(); i++ ){
    HodoscopeLikeHit *hit = blMan->PC(i);
    int seg = hit->seg();
    int au = hit->adcu(), ad = hit->adcd();
    double eu = hit->eu(), ed = hit->ed();
    int tu = hit->tdcu(), td = hit->tdcd();
    double ctm = hit->ctmean();
    double ctsub = hit->ctsub();

    h1 = (TH1F*)gFile->Get( Form("APCU%d",seg) ); h1->Fill( au );
    h1 = (TH1F*)gFile->Get( Form("APCD%d",seg) ); h1->Fill( ad );
    h1 = (TH1F*)gFile->Get( Form("EPCU%d",seg) ); h1->Fill( eu );
    h1 = (TH1F*)gFile->Get( Form("EPCD%d",seg) ); h1->Fill( ed );
    h1 = (TH1F*)gFile->Get( Form("TPCU%d",seg) ); h1->Fill( tu );
    h1 = (TH1F*)gFile->Get( Form("TPCD%d",seg) ); h1->Fill( td );
    //    if( hit->CheckRange() ){
    if( !hit->CheckRange() ){
      h1 = (TH1F*)gFile->Get( Form("AwoTPCU%d",seg) ); h1->Fill( au );
      h1 = (TH1F*)gFile->Get( Form("AwoTPCD%d",seg) ); h1->Fill( ad );
    }else{
      h1 = (TH1F*)gFile->Get( Form("AwTPCU%d",seg) ); h1->Fill( au );
      h1 = (TH1F*)gFile->Get( Form("AwTPCD%d",seg) ); h1->Fill( ad );
      h1 = (TH1F*)gFile->Get( Form("TsubPC%d",seg) );  h1->Fill( ctsub );
      h1 = (TH1F*)gFile->Get( Form("TmeanPC%d",seg) );  h1->Fill( ctm );
      h2 = (TH2F*)gFile->Get( Form("ATPCU%d",seg) );  h2->Fill( tu, au );
      h2 = (TH2F*)gFile->Get( Form("ATPCD%d",seg) );  h2->Fill( td, ad );
      h1 = (TH1F*)gFile->Get( "HitPatPC" ); h1->Fill( seg );
      nPC++;
    }
  }
  h1 = (TH1F*)gFile->Get( "MulPC" ); h1->Fill( nPC );


  evTree->Fill();

  header->Clear();
  blMan->Clear();
  //  cdsMan->Clear();
  return true;
}

void EventAnalysisRawAll::Finalize()
{
  std::cout << " Enter EventAnalysisRawAll::Finalize() " << std::endl;

  rtFile->cd();
  gFile->Write();
  gFile->Close();

  delete blMan;
  //  delete cdsMan;
  delete header;
}

void EventAnalysisRawAll::InitializeHistogram()
{
  rtFile->cd();
  
  // Scaler
  new TH1F( "Scaler", "Scaler", 50, 0, 50 );

  // Trigger Pattern
  new TH1F( "Pattern", "Trigger Pattern", 10, 0, 10 );

#if TKOHIS
  std::cout << "Define Histgram for TKO" << std::endl;
  for( int cr=0; cr<=6; cr++ ){
    for( int sl=1; sl<=22; sl++ ){
      new TH1F( Form("TKOc%ds%d",cr,sl), Form("TKO cr%d sl%d" ,cr,sl), 32, 0, 32+1 );
      for( int ch=0; ch<=32; ch++ )
	{  
	  new TH1F( Form("TKOc%ds%da%d",cr,sl,ch), Form("TKO cr%d sl%d ch%d" ,cr,sl,ch), 4000, 0, 4000 );
	}
    }
  }
#endif



  // NC
  std::cout << "Define Histgram for NC" << std::endl;
  new TH1F( "MulNC", "Multiplicity NC", 112+1, 0, 112+1 );
  new TH1F( "HitPatNC", "Hit Pattern NC", 112+1, 0, 112+1 );
  for( int seg=1; seg<=112; seg++ ){
    new TH1F( Form("ANCU%d",seg),   Form("ADC NCU%d",seg),    4000,    0, 4000 );
    new TH1F( Form("ANCD%d",seg),   Form("ADC NCD%d",seg),    4000,    0, 4000 );
    new TH1F( Form("ENCU%d",seg),   Form("E NCU%d",seg),    4400,    -100, 1000 );
    new TH1F( Form("ENCD%d",seg),   Form("E NCD%d",seg),    4400,    -100, 1000 );
    new TH1F( Form("TNCU%d",seg),   Form("TDC NCU%d",seg),    4000,    0, 4000 );
    new TH1F( Form("TNCD%d",seg),   Form("TDC NCD%d",seg),    4000,    0, 4000 );
    new TH1F( Form("TsubNC%d",seg),   Form("Time sub NC%d",seg),    400,   -10, 10 );
    new TH1F( Form("TmeanNC%d",seg),   Form("Time mean NC%d",seg),    400,   -100, 100 );
    new TH1F( Form("AwoTNCU%d",seg),   Form("ADC woTDC NCU%d",seg),    4000,    0, 4000 );
    new TH1F( Form("AwoTNCD%d",seg),   Form("ADC woTDC NCD%d",seg),    4000,    0, 4000 );
    new TH1F( Form("AwTNCU%d",seg), Form("ADC wTDC NCU%d",seg),  4000,    0, 4000 );
    new TH1F( Form("AwTNCD%d",seg), Form("ADC wTDC NCD%d",seg),  4000,    0, 4000 );
    new TH2F( Form("ATNCU%d",seg),   Form("ADC TDC corr. NCU%d",seg),     200,    0, 4000,  200,    0, 4000 );
    new TH2F( Form("ATNCD%d",seg),   Form("ADC TDC corr. NCD%d",seg),     200,    0, 4000,  200,    0, 4000 );
    new TH2F( Form("AudNC%d",seg),   Form("ADC ud corr. NC%d",seg),     200,    0, 4000,  200,    0, 4000 );    
    new TH2F( Form("TudNC%d",seg),   Form("TDC ud corr. NC%d",seg),     200,    0, 4000,  200,    0, 4000 );    
//    new TH2F( Form("EvTNC%d",seg),   Form("Event TDCsabun corr. NC%d",seg),     1000,    0, 100000,  4000,   -2000, 2000 );
  }  
  
  // PC
  std::cout << "Define Histgram for PC" << std::endl;
  new TH1F( "MulPC", "Multiplicity PC", 32+1, 0, 32+1 );
  new TH1F( "HitPatPC", "Hit Pattern PC", 32+1, 0, 32+1 );
  for( int seg=1; seg<=32; seg++ ){
    new TH1F( Form("APCU%d",seg),   Form("ADC PCU%d",seg),    4000,    0, 4000 );
    new TH1F( Form("APCD%d",seg),   Form("ADC PCD%d",seg),    4000,    0, 4000 );
    new TH1F( Form("EPCU%d",seg),   Form("E PCU%d",seg),    4400,    -100, 1000 );
    new TH1F( Form("EPCD%d",seg),   Form("E PCD%d",seg),    4400,    -100, 1000 );
    new TH1F( Form("TPCU%d",seg),   Form("TDC PCU%d",seg),    4000,    0, 4000 );
    new TH1F( Form("TPCD%d",seg),   Form("TDC PCD%d",seg),    4000,    0, 4000 );
    new TH1F( Form("TsubPC%d",seg),   Form("Time sub PC%d",seg),    400,   -10, 10 );
    new TH1F( Form("TmeanPC%d",seg),   Form("Time mean PC%d",seg),    400,   -100, 100 );
    new TH1F( Form("AwoTPCU%d",seg),   Form("ADC woTDC PCU%d",seg),    4000,    0, 4000 );
    new TH1F( Form("AwoTPCD%d",seg),   Form("ADC woTDC PCD%d",seg),    4000,    0, 4000 );
    new TH1F( Form("AwTPCU%d",seg), Form("ADC wTDC PCU%d",seg),  4000,    0, 4000 );
    new TH1F( Form("AwTPCD%d",seg), Form("ADC wTDC PCD%d",seg),  4000,    0, 4000 );
    new TH2F( Form("ATPCU%d",seg),   Form("ADC TDC corr. PCU%d",seg),     200,    0, 4000,  200,    0, 4000 );
    new TH2F( Form("ATPCD%d",seg),   Form("ADC TDC corr. PCD%d",seg),     200,    0, 4000,  200,    0, 4000 );
  }

}



EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysisRawAll *event = new EventAnalysisRawAll();
  return (EventTemp*)event;
}
