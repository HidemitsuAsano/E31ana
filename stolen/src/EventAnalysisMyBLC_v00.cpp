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

static const int NumOfBLCLayers=8;

#define RUN43 0
#define TKOHIS 0
class EventAnalysisMyBLC: public EventTemp
{
public:
  EventAnalysisMyBLC();
  ~EventAnalysisMyBLC();
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

EventAnalysisMyBLC::EventAnalysisMyBLC()
  : EventTemp()
{

}

EventAnalysisMyBLC::~EventAnalysisMyBLC()
{
}

//const int MaxTreeSize = 19000000000;
void EventAnalysisMyBLC::Initialize( ConfMan *conf )
{
#if 0
  std::cout << " Enter EventAnalysisMyBLC::Initialize " << std::endl;
#endif
  confMan = conf;
  
  rtFile = new TFile( confMan->GetOutFileName().c_str(), "recreate" );
  InitializeHistogram();

  header = new EventHeader();
  if( header==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  blMan = new BeamLineHitMan();
  if( blMan ==NULL ){ std::cerr << "!!!!" << std::endl; return; }

}

void EventAnalysisMyBLC::USca( int nsca, unsigned int *sca )
{
#if 0
  std::cout << " Enter EventAnalysisMyBLC::USca " << std::endl;
#endif
  Block_Event_Number++;
  header->SetBlockEventNumber( Block_Event_Number );
}


bool EventAnalysisMyBLC::UAna( TKOHitCollection *tko )
{
#if 0
  std::cout << " Enter EventAnalysisMyBLC::UAna " << std::endl;
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

    int MulT0=0;
    for(int i=0; i<blMan->nT0(); i++){
      HodoscopeLikeHit* hit = blMan->T0(i);
      if(hit->CheckRange()) MulT0++;
    }

    // Selection //
    if(MulT0!=1){
      header->Clear();
      blMan->Clear();
      return true;
    }

    // BLC1a
    {    
      int nhit[NumOfBLCLayers];
      for(int i=0; i<NumOfBLCLayers; i++){
        nhit[i] = 0;
      }
      for( int layer=1; layer<=NumOfBLCLayers; layer++ ){
        h1 = (TH1F*)gFile->Get( Form("BLC1a_%d_Mult",layer) ); 
        h1->Fill( blMan->nBLC1a(layer) );
        nhit[layer-1]=blMan->nBLC1a(layer);
        for( int i=0; i<blMan->nBLC1a(layer); i++ ){
          ChamberLikeHit *hit = blMan->BLC1a(layer,i);
          int wire = hit->wire();
          int tdc = hit->tdc();
          h1 = (TH1F*)gFile->Get( Form("BLC1a_%d_TDC",layer) ); h1->Fill( tdc );
          h1 = (TH1F*)gFile->Get( Form("BLC1a_%d_%d_TDC",layer,wire) ); h1->Fill( tdc );
          h1 = (TH1F*)gFile->Get( Form("BLC1a_%d_HitPat",layer) ); h1->Fill( wire );
        }
      }
      if( nhit[0]==1 && nhit[NumOfBLCLayers-1]==1){
        for( int layer=1; layer<=NumOfBLCLayers; layer++ ){
          if(nhit[layer-1]==1){
            h1 = (TH1F*)gFile->Get( Form("BLC1a_LayerEfficiency") ); h1->Fill( layer );
          }
        }
      }
    }
    // BLC1b
    {    
      int nhit[NumOfBLCLayers];
      for(int i=0; i<NumOfBLCLayers; i++){
        nhit[i] = 0;
      }
      for( int layer=1; layer<=NumOfBLCLayers; layer++ ){
        h1 = (TH1F*)gFile->Get( Form("BLC1b_%d_Mult",layer) ); 
        h1->Fill( blMan->nBLC1b(layer) );
        nhit[layer-1]=blMan->nBLC1b(layer);
        for( int i=0; i<blMan->nBLC1b(layer); i++ ){
          ChamberLikeHit *hit = blMan->BLC1b(layer,i);
          int wire = hit->wire();
          int tdc = hit->tdc();
          h1 = (TH1F*)gFile->Get( Form("BLC1b_%d_TDC",layer) ); h1->Fill( tdc );
          h1 = (TH1F*)gFile->Get( Form("BLC1b_%d_%d_TDC",layer,wire) ); h1->Fill( tdc );
          h1 = (TH1F*)gFile->Get( Form("BLC1b_%d_HitPat",layer) ); h1->Fill( wire );
        }
      }
      if( nhit[0]==1 && nhit[NumOfBLCLayers-1]==1){
        for( int layer=1; layer<=NumOfBLCLayers; layer++ ){
          if(nhit[layer-1]==1){
            h1 = (TH1F*)gFile->Get( Form("BLC1b_LayerEfficiency") ); h1->Fill( layer );
          }
        }
      }
    }
    // BLC2a
    {    
      int nhit[NumOfBLCLayers];
      for(int i=0; i<NumOfBLCLayers; i++){
        nhit[i] = 0;
      }
      for( int layer=1; layer<=NumOfBLCLayers; layer++ ){
        h1 = (TH1F*)gFile->Get( Form("BLC2a_%d_Mult",layer) );
        h1->Fill( blMan->nBLC2a(layer) );
        nhit[layer-1]=blMan->nBLC2a(layer);
        for( int i=0; i<blMan->nBLC2a(layer); i++ ){
          ChamberLikeHit *hit = blMan->BLC2a(layer,i);
          int wire = hit->wire();
          int tdc = hit->tdc();
          h1 = (TH1F*)gFile->Get( Form("BLC2a_%d_TDC",layer) ); h1->Fill( tdc );
          h1 = (TH1F*)gFile->Get( Form("BLC2a_%d_%d_TDC",layer,wire) ); h1->Fill( tdc );
          h1 = (TH1F*)gFile->Get( Form("BLC2a_%d_HitPat",layer) ); h1->Fill( wire );
        }
      }
      if( nhit[0]==1 && nhit[NumOfBLCLayers-1]==1){
        h1 = (TH1F*)gFile->Get( Form("BLC2a_LayerEfficiency") ); h1->Fill( 0 );
        for( int layer=1; layer<=NumOfBLCLayers; layer++ ){
          if(nhit[layer-1]==1){
            h1 = (TH1F*)gFile->Get( Form("BLC2a_LayerEfficiency") ); h1->Fill( layer );
          }
        }
      }
    }
    // BLC2b
    {    
      int nhit[NumOfBLCLayers];
      for(int i=0; i<NumOfBLCLayers; i++){
        nhit[i] = 0;
      }
      for( int layer=1; layer<=NumOfBLCLayers; layer++ ){
        h1 = (TH1F*)gFile->Get( Form("BLC2b_%d_Mult",layer) ); 
        h1->Fill( blMan->nBLC2b(layer) );
        nhit[layer-1]=blMan->nBLC2b(layer);
        for( int i=0; i<blMan->nBLC2b(layer); i++ ){
          ChamberLikeHit *hit = blMan->BLC2b(layer,i);
          int wire = hit->wire();
          int tdc = hit->tdc();
          h1 = (TH1F*)gFile->Get( Form("BLC2b_%d_TDC",layer) ); h1->Fill( tdc );
          h1 = (TH1F*)gFile->Get( Form("BLC2b_%d_%d_TDC",layer,wire) ); h1->Fill( tdc );
          h1 = (TH1F*)gFile->Get( Form("BLC2b_%d_HitPat",layer) ); h1->Fill( wire );
        }
      }
      if( nhit[0]==1 && nhit[NumOfBLCLayers-1]==1){
        h1 = (TH1F*)gFile->Get( Form("BLC2b_LayerEfficiency") ); h1->Fill( 0 );
        for( int layer=1; layer<=NumOfBLCLayers; layer++ ){
          if(nhit[layer-1]==1){
            h1 = (TH1F*)gFile->Get( Form("BLC2b_LayerEfficiency") ); h1->Fill( layer );
          }
        }
      }
    }

    header->Clear();
    blMan->Clear();
    return true;
}

void EventAnalysisMyBLC::Finalize()
{
  std::cout << " Enter EventAnalysisMyBLC::Finalize " << std::endl;

  rtFile->cd();
  gFile->Write();
  gFile->Close();

  delete blMan;
  delete header;
}

void EventAnalysisMyBLC::InitializeHistogram()
{
  const int ndc=4;
  int dccid[ndc]={CID_BLC1a,CID_BLC1b,CID_BLC2a,CID_BLC2b};
  TString dcname[ndc]={"BLC1a","BLC1b","BLC2a","BLC2b"};
  int NumOfLayers[ndc]={NumOfBLCLayers,NumOfBLCLayers,NumOfBLCLayers,NumOfBLCLayers};

  for(int idc=0;idc<ndc;idc++){
    std::cout << "Define Histgram for " << dcname[idc] << std::endl;
    for( int layer=1; layer<=NumOfLayers[idc]; layer++ ){
      int nwire = confMan->GetBLDCWireMapManager()->GetNWire( dccid[idc], layer );
      new TH1F( Form("%s_%d_Mult",dcname[idc].Data(),layer), Form("Multiplicity %s Layer%d;Multiplicity;Counts",dcname[idc].Data(),layer), nwire+1, 0, nwire+1 );
      new TH1F( Form("%s_%d_HitPat",dcname[idc].Data(),layer), Form("Hit Pattern %s Layer%d;Wire;Counts",dcname[idc].Data(),layer), nwire+1, 0, nwire+1 );
      new TH1F( Form("%s_%d_TDC",dcname[idc].Data(),layer), Form("TDC %s Layer%d;TDC ch.;Counts",dcname[idc].Data(),layer), 4000, 0, 4000 );
      for( int wire=1; wire<=nwire; wire++ ){
        new TH1F( Form("%s_%d_%d_TDC",dcname[idc].Data(),layer,wire), Form("TDC %s Layer%d Wire%d;TDC ch.;Counts",dcname[idc].Data(),layer,wire), 4000, 0, 4000 );
      }
    }
    new TH1F(Form("%s_LayerEfficiency",dcname[idc].Data()),Form("%s Layer Efficiency;Layer;Counts",dcname[idc].Data()),NumOfLayers[idc]+1,0,NumOfLayers[idc]+1);
  }
}



EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysisMyBLC *event = new EventAnalysisMyBLC();
  return (EventTemp*)event;
}
