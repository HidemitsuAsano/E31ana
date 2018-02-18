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

class EventAnalysis: public EventTemp
{
public:
  EventAnalysis();
  ~EventAnalysis();
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
};

EventAnalysis::EventAnalysis()
  : EventTemp()
{
}

EventAnalysis::~EventAnalysis()
{
}

const int MaxTreeSize = 1900000000000;
void EventAnalysis::Initialize( ConfMan *conf )
{
#if 1
  std::cout << " Enter EventAnalysis::Initialize " << std::endl;
#endif
  confMan = conf;
  
  rtFile = new TFile( confMan->GetOutFileName().c_str(), "recreate" );
  evTree = new TTree( "EventTree", "EventTree" );
  scaTree= new TTree( "ScalerTree", "ScalerTree" );

  header = new EventHeader();
  if( header==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  evTree->Branch( "EventHeader", &header );
  cdsMan = new CDSHitMan();
  if( cdsMan==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  evTree->Branch( "CDSHitMan", &cdsMan );
  blMan = new BeamLineHitMan();
  if( blMan ==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  evTree->Branch( "BeamLineHitMan", &blMan );
  scaMan = new ScalerMan();
  if( scaMan==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  scaTree->Branch( "ScalerMan", &scaMan );

  new TH1F( "Pattern", "Pattern", 10, 0, 10 );
  for( int layer=1; layer<=NumOfCDCLayers; layer++ ){
    int nwire = (int)360/(conf->GetCDCWireMapManager()->dphi(layer));
    new TH1F( Form("CDCHitPattern%d",layer), Form("CDCHitPattern%d",layer), nwire, 0, nwire );
    new TH1F( Form("CDCTDC%d",layer), Form("CDCTDC%d",layer), 1500, 0, 1500 );
  }

  new TH1F( "CDHHitPattern", "CDHHitPattern", 40, 0, 40 );
  for( int seg=1; seg<=NumOfCDHSegments; seg++ ){
    new TH1F( Form("CDHADCU%d",seg), Form("CDHADCU%d",seg), 2000, 0, 2000 );
    new TH1F( Form("CDHADCD%d",seg), Form("CDHADCD%d",seg), 2000, 0, 2000 );
    new TH1F( Form("CDHTDCU%d",seg), Form("CDHTDCU%d",seg), 2000, 0, 2000 );
    new TH1F( Form("CDHTDCD%d",seg), Form("CDHTDCD%d",seg), 2000, 0, 2000 );
    new TH1F( Form("CDHADCUWT%d",seg), Form("CDHADCUWT%d",seg), 2000, 0, 2000 );
    new TH1F( Form("CDHADCDWT%d",seg), Form("CDHADCDWT%d",seg), 2000, 0, 2000 );    
  }

}

void EventAnalysis::USca( int nsca, unsigned int *sca )
{
#if 0
  std::cout << " Enter EventAnalysis::USca " << std::endl;
#endif
  Block_Event_Number++;
  header->SetBlockEventNumber( Block_Event_Number );
  scaMan->SetBlockEventNumber( Block_Event_Number );
  for( int i=0; i<nsca; i++ ){
   scaMan->AddHit( sca[i], confMan->GetScalerMapManager()->GetName(i) );
  }
#if 0
  std::cout << nsca<< std::endl;
  for( int i=0; i<nsca; i++ ){
   std::cout << "  " << sca[i];
  }
  std::cout << std::endl;
#endif
  scaTree->Fill();
  scaMan->Clear();
}

bool EventAnalysis::UAna( TKOHitCollection *tko )
{
#if 0
  std::cout << " Enter EventAnalysis::UAna " << std::endl;
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

  TH1F *h1;
    /*** after here, we only fill histograms. ***/
    for( int i=0; i<20; i++ ){
      int val = header->pattern(i);
      if( 0<val ){ h1=(TH1F*)gFile->Get("Pattern"); h1->Fill(i); }
    }
    
    for( int layer=1; layer<=NumOfCDCLayers; layer++ ){
      //std::cout << "Layer" << layer << "  #ev=" << cdsMan->nCDC(layer) << std::endl;
      for( int i=0; i<cdsMan->nCDC(layer); i++ ){
	int wire = cdsMan->CDC(layer,i)->wire();
	int tdc = cdsMan->CDC(layer,i)->tdc();
	//std::cout << "   wire=" << wire << " tdc=" << tdc << std::endl;
	h1 = (TH1F*)gFile->Get( Form("CDCHitPattern%d",layer) ); h1->Fill(wire);
	h1 = (TH1F*)gFile->Get( Form("CDCTDC%d",layer) ); h1->Fill(tdc);
      }      
    }

    for( int i=0; i<cdsMan->nCDH(); i++ ){
      int seg  = cdsMan->CDH(i)->seg();
      int tdcu = cdsMan->CDH(i)->tdcu();
      int tdcd = cdsMan->CDH(i)->tdcd();
      int adcu = cdsMan->CDH(i)->adcu();
      int adcd = cdsMan->CDH(i)->adcd();
//       std::cout << " seg:" << seg << " tdcu:" << tdcu << " tdcd:" << tdcd
// 		<< " adcu:" << adcu << " adcd:" << adcd << std::endl;
      h1 = (TH1F*)gFile->Get( Form("CDHADCU%d",seg) ); h1->Fill(adcu);
      h1 = (TH1F*)gFile->Get( Form("CDHADCD%d",seg) ); h1->Fill(adcd);
      h1 = (TH1F*)gFile->Get( Form("CDHTDCU%d",seg) ); h1->Fill(tdcu);
      h1 = (TH1F*)gFile->Get( Form("CDHTDCD%d",seg) ); h1->Fill(tdcd);
      if( 0<tdcu && tdcu<4096 ){
	h1 = (TH1F*)gFile->Get( Form("CDHADCUWT%d",seg) ); h1->Fill(adcu);
      }
      if( 0<tdcd && tdcd<4096 ){
	h1 = (TH1F*)gFile->Get( Form("CDHADCDWT%d",seg) ); h1->Fill(adcd);
      }
      if( cdsMan->CDH(i)->CheckRange() ){
	h1 = (TH1F*)gFile->Get("CDHHitPattern"); h1->Fill(seg);
      }
    }


  evTree->Fill();

  header->Clear();
  blMan->Clear();
  cdsMan->Clear();
  return true;
}

void EventAnalysis::Finalize()
{
  std::cout << " Enter EventAnalysis::Finalize " << std::endl;


  gFile->Write();
  gFile->Close();

  delete blMan;
  delete cdsMan;
  delete header;
}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysis *event = new EventAnalysis();
  return (EventTemp*)event;
}
