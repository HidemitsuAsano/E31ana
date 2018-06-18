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
#include "GlobalVariables.h"

#define RUN43 0
#define TKOHIS 0
class EventAnalysisFDC: public EventTemp
{
public:
  EventAnalysisFDC();
  ~EventAnalysisFDC();
private:
  TFile *rtFile;
  TTree *evTree;
  TTree *scaTree;
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

EventAnalysisFDC::EventAnalysisFDC()
  : EventTemp()
{

}

EventAnalysisFDC::~EventAnalysisFDC()
{
}

//const int MaxTreeSize = 19000000000;
void EventAnalysisFDC::Initialize( ConfMan *conf )
{
#if 0
  std::cout << " Enter EventAnalysisFDC::Initialize " << std::endl;
#endif
  confMan = conf;
  
  rtFile = new TFile( confMan->GetOutFileName().c_str(), "recreate" );
  InitializeHistogram();

  header = new EventHeader();
  if( header==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  blMan = new BeamLineHitMan();
  if( blMan ==NULL ){ std::cerr << "!!!!" << std::endl; return; }

}

void EventAnalysisFDC::USca( int nsca, unsigned int *sca )
{
#if 0
  std::cout << " Enter EventAnalysisFDC::USca " << std::endl;
#endif
  Block_Event_Number++;
  header->SetBlockEventNumber( Block_Event_Number );
}


bool EventAnalysisFDC::UAna( TKOHitCollection *tko )
{
#if 0
  std::cout << " Enter EventAnalysisFDC::UAna " << std::endl;
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

  header->Convert( tko, confMan );
  blMan->Convert( tko, confMan );
  rtFile->cd();

  int nBVC=0;
  int BVCseg=-1;

  for( int i=0; i<blMan->nBVC();i++){
    HodoscopeLikeHit *hit=blMan->BVC(i);
    if(hit->CheckRange()){
      nBVC++;
      BVCseg=hit->seg();
    }
  }
  h1 = (TH1F*)gFile->Get( Form("MulBVC") ); h1->Fill( nBVC );

  int nhit[6];
  for( int layer=1; layer<=NumOfFDC1Layers; layer++ ){
    h1 = (TH1F*)gFile->Get( Form("MulFDC1_%d",layer) ); 
    h1->Fill( blMan->nFDC1(layer) );
    nhit[layer-1]=blMan->nFDC1(layer);
    for( int i=0; i<blMan->nFDC1(layer); i++ ){
      ChamberLikeHit *hit = blMan->FDC1(layer,i);
      int wire = hit->wire();
      int tdc = hit->tdc();
      h1 = (TH1F*)gFile->Get( Form("TFDC1_%d",layer) ); h1->Fill( tdc );
      h1 = (TH1F*)gFile->Get( Form("TFDC1_%d_%d",layer,wire) ); h1->Fill( tdc );
      h1 = (TH1F*)gFile->Get( Form("HitPatFDC1_%d",layer) ); h1->Fill( wire );
    }
  }

  if( nhit[0]==1 && nhit[5]==1){
    for( int layer=1; layer<=NumOfFDC1Layers; layer++ )
      if(nhit[layer-1]>=1)
	h1 = (TH1F*)gFile->Get( Form("FDC1Efficiency") ); h1->Fill( layer );
  }
    
  header->Clear();
  blMan->Clear();
  return true;
}

void EventAnalysisFDC::Finalize()
{
  std::cout << " Enter EventAnalysisFDC::Finalize " << std::endl;

  rtFile->cd();
  gFile->Write();
  gFile->Close();

  delete blMan;
  delete header;
}

void EventAnalysisFDC::InitializeHistogram()
{
  const int ndc=1;
  int dccid[ndc]={CID_FDC1};
  TString dcname[ndc]={"FDC1"};
  int NumOfLayers[ndc]={NumOfFDC1Layers};

  new TH1F("MulBVC","MulBVC",8,0.5,8.5);
  new TH1F("FDC1Efficiency","FDC1Efficiency",6,0.5,6.5);
  for(int idc=0;idc<ndc;idc++){
    std::cout << "Define Histgram for " << dcname[idc] << std::endl;
    for( int layer=1; layer<=NumOfLayers[idc]; layer++ ){
      int nwire = confMan->GetBLDCWireMapManager()->GetNWire( dccid[idc], layer );
      new TH1F( Form("Mul%s_%d",dcname[idc].Data(),layer), Form("Multiplicity %s Layer%d",dcname[idc].Data(),layer), nwire, 0, nwire+1 );
      new TH1F( Form("HitPat%s_%d",dcname[idc].Data(),layer), Form("Hit Pattern %s Layer%d",dcname[idc].Data(),layer), nwire, 0, nwire+1 );
      new TH1F( Form("T%s_%d",dcname[idc].Data(),layer), Form("TDC %s Layer%d",dcname[idc].Data(),layer), 1500, 0, 1500 );
      for( int wire=1; wire<=nwire; wire++ ){
	new TH1F( Form("T%s_%d_%d",dcname[idc].Data(),layer,wire), Form("TDC %s Layer%d Wire%d",dcname[idc].Data(),layer,wire), 1500, 0, 1500 );
      }
    }
  }
}



EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysisFDC *event = new EventAnalysisFDC();
  return (EventTemp*)event;
}
