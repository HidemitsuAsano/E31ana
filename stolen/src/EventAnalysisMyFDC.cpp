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
class EventAnalysisMyFDC: public EventTemp
{
  public:
    EventAnalysisMyFDC();
    ~EventAnalysisMyFDC();
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

  EventAnalysisMyFDC::EventAnalysisMyFDC()
: EventTemp()
{

}

EventAnalysisMyFDC::~EventAnalysisMyFDC()
{
}

//const int MaxTreeSize = 19000000000;
void EventAnalysisMyFDC::Initialize( ConfMan *conf )
{
#if 0
  std::cout << " Enter EventAnalysisMyFDC::Initialize " << std::endl;
#endif
  confMan = conf;

  rtFile = new TFile( confMan->GetOutFileName().c_str(), "recreate" );
  InitializeHistogram();

  header = new EventHeader();
  if( header==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  blMan = new BeamLineHitMan();
  if( blMan ==NULL ){ std::cerr << "!!!!" << std::endl; return; }

}

void EventAnalysisMyFDC::USca( int nsca, unsigned int *sca )
{
#if 0
  std::cout << " Enter EventAnalysisMyFDC::USca " << std::endl;
#endif
  Block_Event_Number++;
  header->SetBlockEventNumber( Block_Event_Number );
}


bool EventAnalysisMyFDC::UAna( TKOHitCollection *tko )
{
#if 0
  std::cout << " Enter EventAnalysisMyFDC::UAna " << std::endl;
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

    //    int MulT0=0;
    //    for(int i=0; i<blMan->nT0(); i++){
    //      HodoscopeLikeHit* hit = blMan->T0(i);
    //      if(hit->CheckRange()) MulT0++;
    //    }
    //
    //    // Selection //
    //    if(MulT0!=1){
    //      header->Clear();
    //      blMan->Clear();
    //      return true;
    //    }

    // FDC1
    {    
      int nhit[NumOfFDCLayers];
      for(int i=0; i<NumOfFDCLayers; i++){
        nhit[i] = 0;
      }
      for( int layer=1; layer<=NumOfFDCLayers; layer++ ){
        h1 = (TH1F*)gFile->Get( Form("FDC1_Layer%d_Mult",layer) ); 
        h1->Fill( blMan->nFDC1(layer) );
        nhit[layer-1]=blMan->nFDC1(layer);
        for( int i=0; i<blMan->nFDC1(layer); i++ ){
          ChamberLikeHit *hit = blMan->FDC1(layer,i);
          int wire = hit->wire();
          int tdc = hit->tdc();
          h1 = (TH1F*)gFile->Get( Form("FDC1_Layer%d_TDC",layer) ); h1->Fill( tdc );
          h1 = (TH1F*)gFile->Get( Form("FDC1_Layer%d_Wire%d_TDC",layer,wire) ); h1->Fill( tdc );
          h1 = (TH1F*)gFile->Get( Form("FDC1_Layer%d_HitPat",layer) ); h1->Fill( wire );
          for(int mult=1; mult<=3; mult++){
            if(nhit[layer-1]==mult){
              h1 = (TH1F*)gFile->Get( Form("FDC1_Layer%d_HitPat_mult%d",layer,mult) ); h1->Fill( wire );
            }
          }
        }
      }
      for(int mult=1; mult<=9; mult++){
        h1 = (TH1F*)gFile->Get( Form("FDC1_LayerEfficiency_mult%d",mult) ); h1->Fill( 0 );
        h1 = (TH1F*)gFile->Get( Form("FDC1_LayerEfficiency_multleq%d",mult) ); h1->Fill( 0 );
        if( nhit[0]==1 && nhit[NumOfFDCLayers-1]==1){
          h1 = (TH1F*)gFile->Get( Form("FDC1_LayerEfficiency_mult%d",mult) ); h1->Fill( 1 );
          h1 = (TH1F*)gFile->Get( Form("FDC1_LayerEfficiency_multleq%d",mult) ); h1->Fill( 1 );
          h1 = (TH1F*)gFile->Get( Form("FDC1_LayerEfficiency_mult%d",mult) ); h1->Fill( NumOfFDCLayers );
          h1 = (TH1F*)gFile->Get( Form("FDC1_LayerEfficiency_multleq%d",mult) ); h1->Fill( NumOfFDCLayers );
          for( int layer=2; layer<NumOfFDCLayers; layer++ ){
            if(nhit[layer-1]==mult){
              h1 = (TH1F*)gFile->Get( Form("FDC1_LayerEfficiency_mult%d",mult) ); h1->Fill( layer );
            }
            if(nhit[layer-1]<=mult && nhit[layer-1]!=0){
              h1 = (TH1F*)gFile->Get( Form("FDC1_LayerEfficiency_multleq%d",mult) ); h1->Fill( layer );
            }
          }
        }
        bool hitflag1 = false;
        for(int layer=1; layer<=NumOfFDCLayers; layer++){
          for(int i=1; i<=NumOfFDCLayers; i++){
            if(i!=layer){
              if(nhit[i-1]==1){
                hitflag1 = true;
              }
              else{
                hitflag1 = false;
              }
              if(!hitflag1) break;
            }
          }
          h1 = (TH1F*)gFile->Get( Form("FDC1_Layer%dEfficiency_mult%d",layer,mult) ); h1->Fill( 0 );
          h1 = (TH1F*)gFile->Get( Form("FDC1_Layer%dEfficiency_multleq%d",layer,mult) ); h1->Fill( 0 );
          if(hitflag1){
            h1 = (TH1F*)gFile->Get( Form("FDC1_Layer%dEfficiency_mult%d",layer,mult) ); h1->Fill( 1 );
            if(nhit[layer-1]==mult){
              h1 = (TH1F*)gFile->Get( Form("FDC1_Layer%dEfficiency_mult%d",layer,mult) ); h1->Fill( 2 );
            }
            h1 = (TH1F*)gFile->Get( Form("FDC1_Layer%dEfficiency_multleq%d",layer,mult) ); h1->Fill( 1 );
            if(nhit[layer-1]<=mult){
              h1 = (TH1F*)gFile->Get( Form("FDC1_Layer%dEfficiency_multleq%d",layer,mult) ); h1->Fill( 2 );
            }
          }
        }
      }
    }

    header->Clear();
    blMan->Clear();
    return true;
}

void EventAnalysisMyFDC::Finalize()
{
  std::cout << " Enter EventAnalysisMyFDC::Finalize " << std::endl;

  rtFile->cd();
  gFile->Write();
  gFile->Close();

  delete blMan;
  delete header;
}

void EventAnalysisMyFDC::InitializeHistogram()
{
  const int ndc=1;
  int dccid[ndc]={CID_FDC1};
  TString dcname[ndc]={"FDC1"};
  int NumOfLayers[ndc]={NumOfFDC1Layers};

  for(int idc=0;idc<ndc;idc++){
    std::cout << "Define Histgram for " << dcname[idc] << std::endl;
    for( int layer=1; layer<=NumOfLayers[idc]; layer++ ){
      int nwire = confMan->GetBLDCWireMapManager()->GetNWire( dccid[idc], layer );
      new TH1F( Form("%s_Layer%d_Mult",dcname[idc].Data(),layer), Form("Multiplicity %s Layer%d;Multiplicity;Counts",dcname[idc].Data(),layer), nwire+1, 0, nwire+1 );
      new TH1F( Form("%s_Layer%d_HitPat",dcname[idc].Data(),layer), Form("Hit Pattern %s Layer%d;Wire;Counts",dcname[idc].Data(),layer), nwire+1, 0, nwire+1 );
      new TH1F( Form("%s_Layer%d_HitPat_mult1",dcname[idc].Data(),layer), Form("Hit Pattern %s Layer%d (mult = 1);Wire;Counts",dcname[idc].Data(),layer), nwire+1, 0, nwire+1 );
      new TH1F( Form("%s_Layer%d_HitPat_mult2",dcname[idc].Data(),layer), Form("Hit Pattern %s Layer%d (mult = 2);Wire;Counts",dcname[idc].Data(),layer), nwire+1, 0, nwire+1 );
      new TH1F( Form("%s_Layer%d_HitPat_mult3",dcname[idc].Data(),layer), Form("Hit Pattern %s Layer%d (mult = 3);Wire;Counts",dcname[idc].Data(),layer), nwire+1, 0, nwire+1 );
      new TH1F( Form("%s_Layer%d_TDC",dcname[idc].Data(),layer), Form("TDC %s Layer%d;TDC ch.;Counts",dcname[idc].Data(),layer), 4000, 0, 4000 );
      for( int wire=1; wire<=nwire; wire++ ){
        new TH1F( Form("%s_Layer%d_Wire%d_TDC",dcname[idc].Data(),layer,wire), Form("TDC %s Layer%d Wire%d;TDC ch.;Counts",dcname[idc].Data(),layer,wire), 4000, 0, 4000 );
      }
    }
    for(int mult=1; mult<=9; mult++){
      new TH1F(Form("%s_LayerEfficiency_mult%d",dcname[idc].Data(),mult),Form("%s Layer Efficiency (mult=%d);Layer;Counts",dcname[idc].Data(),mult),NumOfLayers[idc]+1,0,NumOfLayers[idc]+1);
      new TH1F(Form("%s_LayerEfficiency_multleq%d",dcname[idc].Data(),mult),Form("%s Layer Efficiency (mult#leq%d);Layer;Counts",dcname[idc].Data(),mult),NumOfLayers[idc]+1,0,NumOfLayers[idc]+1);
      for( int layer=1; layer<=NumOfLayers[idc]; layer++ ){
        new TH1F(Form("%s_Layer%dEfficiency_mult%d",dcname[idc].Data(),layer,mult),Form("%s Layer %d Efficiency (mult=%d);;Counts",dcname[idc].Data(),layer,mult),3,0,3);
        new TH1F(Form("%s_Layer%dEfficiency_multleq%d",dcname[idc].Data(),layer,mult),Form("%s Layer %d Efficiency (mult#leq%d);;Counts",dcname[idc].Data(),layer,mult),3,0,3);
      }
    }
  }
}



EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysisMyFDC *event = new EventAnalysisMyFDC();
  return (EventTemp*)event;
}
