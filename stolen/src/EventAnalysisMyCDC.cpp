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
class EventAnalysisMyCDC: public EventTemp
{
  public:
    EventAnalysisMyCDC();
    ~EventAnalysisMyCDC();
  private:
    TFile *rtFile;
    TTree *evTree;
    TTree *scaTree;
    BeamLineHitMan *blMan;
    CDSHitMan *cdsMan;
    EventHeader *header;
    ScalerMan *scaMan;

  public:
    void Initialize( ConfMan *conf );
    void USca( int nsca, unsigned int *sca );
    bool UAna( TKOHitCollection *tko );
    void Finalize();

    void InitializeHistogram();
};

  EventAnalysisMyCDC::EventAnalysisMyCDC()
: EventTemp()
{

}

EventAnalysisMyCDC::~EventAnalysisMyCDC()
{
}

//const int MaxTreeSize = 19000000000;
void EventAnalysisMyCDC::Initialize( ConfMan *conf )
{
#if 0
  std::cout << " Enter EventAnalysisMyCDC::Initialize " << std::endl;
#endif
  confMan = conf;

  rtFile = new TFile( confMan->GetOutFileName().c_str(), "recreate" );
  InitializeHistogram();

  header = new EventHeader();
  if( header==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  blMan = new BeamLineHitMan();
  if( blMan ==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  cdsMan = new CDSHitMan();
  if( cdsMan ==NULL ){ std::cerr << "!!!!" << std::endl; return; }

}

void EventAnalysisMyCDC::USca( int nsca, unsigned int *sca )
{
#if 0
  std::cout << " Enter EventAnalysisMyCDC::USca " << std::endl;
#endif
  Block_Event_Number++;
  header->SetBlockEventNumber( Block_Event_Number );
}


bool EventAnalysisMyCDC::UAna( TKOHitCollection *tko )
{
#if 0
  std::cout << " Enter EventAnalysisMyCDC::UAna " << std::endl;
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
    cdsMan->Convert( tko, confMan );
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
      cdsMan->Clear();
      return true;
    }
    // CDC
    {    
      int nhit[NumOfCDCLayers];
      for(int i=0; i<NumOfCDCLayers; i++){
        nhit[i] = 0;
      }
      for( int layer=1; layer<=NumOfCDCLayers; layer++ ){
        h1 = (TH1F*)gFile->Get( Form("CDC_Layer%d_Mult",layer) ); 
        h1->Fill( cdsMan->nCDC(layer) );
        nhit[layer-1]=cdsMan->nCDC(layer);
        for( int i=0; i<cdsMan->nCDC(layer); i++ ){
          CDCHit *hit = cdsMan->CDC(layer,i);
          int wire = hit->wire();
          int tdc = hit->tdc();
          double dt = hit->dt(), cdt = hit->cdt();
          double dl = hit->dl(), cdl = hit->cdl();
          h1 = (TH1F*)gFile->Get( Form("CDC_Layer%d_TDC",layer) ); h1->Fill( tdc );
          h1 = (TH1F*)gFile->Get( Form("CDC_Layer%d_Wire%d_TDC",layer,wire) ); h1->Fill( tdc );
          h1 = (TH1F*)gFile->Get( Form("CDC_Layer%d_Time",layer) ); h1->Fill( dt );
          h1 = (TH1F*)gFile->Get( Form("CDC_Layer%d_Wire%d_Time",layer,wire) ); h1->Fill( dt );
          //h1 = (TH1F*)gFile->Get( Form("CDC_Layer%d_CTime",layer) ); h1->Fill( cdt );
          //h1 = (TH1F*)gFile->Get( Form("CDC_Layer%d_Wire%d_CTime",layer,wire) ); h1->Fill( cdt );
          //h1 = (TH1F*)gFile->Get( Form("CDC_Layer%d_Length",layer) ); h1->Fill( dl );
          //h1 = (TH1F*)gFile->Get( Form("CDC_Layer%d_Wire%d_Length",layer,wire) ); h1->Fill( dl );
          //h1 = (TH1F*)gFile->Get( Form("CDC_Layer%d_CLength",layer) ); h1->Fill( cdl );
          //h1 = (TH1F*)gFile->Get( Form("CDC_Layer%d_Wire%d_CLength",layer,wire) ); h1->Fill( cdl );
          //h2 = (TH2F*)gFile->Get( Form("CDC_Layer%d_TimevLength",layer) ); h2->Fill( dt,dl );
          //h2 = (TH2F*)gFile->Get( Form("CDC_Layer%d_Wire%d_TimevLength",layer,wire) ); h2->Fill( dt,dl );
          //h2 = (TH2F*)gFile->Get( Form("CDC_Layer%d_CTimevLength",layer) ); h2->Fill( cdt,dl );
          //h2 = (TH2F*)gFile->Get( Form("CDC_Layer%d_Wire%d_CTimevLength",layer,wire) ); h2->Fill( cdt,dl );
          h1 = (TH1F*)gFile->Get( Form("CDC_Layer%d_HitPat",layer) ); h1->Fill( wire );
          for(int mult=1; mult<=3; mult++){
            if(nhit[layer-1]==mult){
              h1 = (TH1F*)gFile->Get( Form("CDC_Layer%d_HitPat_mult%d",layer,mult) ); h1->Fill( wire );
            }
          }
        }
      }
      for(int mult=1; mult<=9; mult++){
        h1 = (TH1F*)gFile->Get( Form("CDC_LayerEfficiency_mult%d",mult) ); h1->Fill( 0 );
        h1 = (TH1F*)gFile->Get( Form("CDC_LayerEfficiency_multleq%d",mult) ); h1->Fill( 0 );
        if( nhit[0]==1 && nhit[NumOfCDCLayers-1]==1){
          h1 = (TH1F*)gFile->Get( Form("CDC_LayerEfficiency_mult%d",mult) ); h1->Fill( 1 );
          h1 = (TH1F*)gFile->Get( Form("CDC_LayerEfficiency_multleq%d",mult) ); h1->Fill( 1 );
          h1 = (TH1F*)gFile->Get( Form("CDC_LayerEfficiency_mult%d",mult) ); h1->Fill( NumOfCDCLayers );
          h1 = (TH1F*)gFile->Get( Form("CDC_LayerEfficiency_multleq%d",mult) ); h1->Fill( NumOfCDCLayers );
          for( int layer=2; layer<NumOfCDCLayers; layer++ ){
            if(nhit[layer-1]==mult){
              h1 = (TH1F*)gFile->Get( Form("CDC_LayerEfficiency_mult%d",mult) ); h1->Fill( layer );
            }
            if(nhit[layer-1]<=mult && nhit[layer-1]!=0){
              h1 = (TH1F*)gFile->Get( Form("CDC_LayerEfficiency_multleq%d",mult) ); h1->Fill( layer );
            }
          }
        }
        bool hitflag1 = false;
        for(int layer=1; layer<=NumOfCDCLayers; layer++){
          for(int i=1; i<=NumOfCDCLayers; i++){
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
          h1 = (TH1F*)gFile->Get( Form("CDC_Layer%dEfficiency_mult%d",layer,mult) ); h1->Fill( 0 );
          h1 = (TH1F*)gFile->Get( Form("CDC_Layer%dEfficiency_multleq%d",layer,mult) ); h1->Fill( 0 );
          if(hitflag1){
            h1 = (TH1F*)gFile->Get( Form("CDC_Layer%dEfficiency_mult%d",layer,mult) ); h1->Fill( 1 );
            if(nhit[layer-1]==mult){
              h1 = (TH1F*)gFile->Get( Form("CDC_Layer%dEfficiency_mult%d",layer,mult) ); h1->Fill( 2 );
            }
            h1 = (TH1F*)gFile->Get( Form("CDC_Layer%dEfficiency_multleq%d",layer,mult) ); h1->Fill( 1 );
            if(nhit[layer-1]<=mult){
              h1 = (TH1F*)gFile->Get( Form("CDC_Layer%dEfficiency_multleq%d",layer,mult) ); h1->Fill( 2 );
            }
          }
        }
      }
    }

    header->Clear();
    blMan->Clear();
    cdsMan->Clear();
    return true;
}

void EventAnalysisMyCDC::Finalize()
{
  std::cout << " Enter EventAnalysisMyCDC::Finalize " << std::endl;

  rtFile->cd();
  gFile->Write();
  gFile->Close();

  delete cdsMan;
  delete blMan;
  delete header;
}

void EventAnalysisMyCDC::InitializeHistogram()
{
  const int ndc=1;
  int dccid[ndc]={CID_CDC};
  TString dcname[ndc]={"CDC"};
  int NumOfLayers[ndc]={NumOfCDCLayers};

  for(int idc=0;idc<ndc;idc++){
    std::cout << "Define Histgram for " << dcname[idc] << std::endl;
    for( int layer=1; layer<=NumOfLayers[idc]; layer++ ){
      int nwire = NumOfCDCWiresInLayer[layer-1];
      new TH1F( Form("%s_Layer%d_Mult",dcname[idc].Data(),layer), Form("Multiplicity %s Layer%d;Multiplicity;Counts",dcname[idc].Data(),layer), nwire+1, -0.5, nwire+0.5 );
      new TH1F( Form("%s_Layer%d_HitPat",dcname[idc].Data(),layer), Form("Hit Pattern %s Layer%d;Wire;Counts",dcname[idc].Data(),layer), nwire+1, 0, nwire+1 );
      new TH1F( Form("%s_Layer%d_HitPat_mult1",dcname[idc].Data(),layer), Form("Hit Pattern %s Layer%d (mult = 1);Wire;Counts",dcname[idc].Data(),layer), nwire+1, 0, nwire+1 );
      new TH1F( Form("%s_Layer%d_HitPat_mult2",dcname[idc].Data(),layer), Form("Hit Pattern %s Layer%d (mult = 2);Wire;Counts",dcname[idc].Data(),layer), nwire+1, 0, nwire+1 );
      new TH1F( Form("%s_Layer%d_HitPat_mult3",dcname[idc].Data(),layer), Form("Hit Pattern %s Layer%d (mult = 3);Wire;Counts",dcname[idc].Data(),layer), nwire+1, 0, nwire+1 );
      new TH1F( Form("%s_Layer%d_TDC",dcname[idc].Data(),layer), Form("TDC %s Layer%d;TDC ch.;Counts",dcname[idc].Data(),layer), 4000, 0, 4000 );
      for( int wire=1; wire<=nwire; wire++ ){
        new TH1F( Form("%s_Layer%d_Wire%d_TDC",dcname[idc].Data(),layer,wire), Form("TDC %s Layer%d Wire%d;TDC ch.;Counts",dcname[idc].Data(),layer,wire), 4000, 0, 4000 );
      }
      new TH1F( Form("%s_Layer%d_Time",dcname[idc].Data(),layer), Form("Drift time %s Layer%d;Drift time (ns);Counts",dcname[idc].Data(),layer), 4000, -500, 1500 );
      //new TH1F( Form("%s_Layer%d_CTime",dcname[idc].Data(),layer), Form("Drift time %s Layer%d;Drift time (ns);Counts",dcname[idc].Data(),layer), 4000, -500, 1500 );
      //new TH1F( Form("%s_Layer%d_Length",dcname[idc].Data(),layer), Form("Drift length %s Layer%d;Drift length (mm);Counts",dcname[idc].Data(),layer), 1100, -0.1, 1.0 );
      //new TH1F( Form("%s_Layer%d_CLength",dcname[idc].Data(),layer), Form("Drift length %s Layer%d;Drift length (mm);Counts",dcname[idc].Data(),layer), 1100, -0.1, 1.0 );
      //new TH2F( Form("%s_Layer%d_TimevLength",dcname[idc].Data(),layer), Form("Drift time vs. Drift length %s Layer%d;Drift time (ns);Drift length (mm)",dcname[idc].Data(),layer), 400, -500, 1500, 110, -0.1, 1.0 );
      //new TH2F( Form("%s_Layer%d_CTimevLength",dcname[idc].Data(),layer), Form("Drift time vs. Drift length %s Layer%d;Drift time (ns);Drift length (mm)",dcname[idc].Data(),layer), 400, -500, 1500, 110, -0.1, 1.0 );
      for( int wire=1; wire<=nwire; wire++ ){
        new TH1F( Form("%s_Layer%d_Wire%d_Time",dcname[idc].Data(),layer,wire), Form("Drift time %s Layer%d Wire%d;Drift time (ns);Counts",dcname[idc].Data(),layer,wire), 4000, -500, 1500 );
        //new TH1F( Form("%s_Layer%d_Wire%d_CTime",dcname[idc].Data(),layer,wire), Form("Drift time %s Layer%d Wire%d;Drift time (ns);Counts",dcname[idc].Data(),layer,wire), 4000, -500, 1500 );
        //new TH1F( Form("%s_Layer%d_Wire%d_Length",dcname[idc].Data(),layer,wire), Form("Drift length %s Layer%d Wire%d;Drift Length (mm);Counts",dcname[idc].Data(),layer,wire), 1100, -0.1,1.0 );
        //new TH1F( Form("%s_Layer%d_Wire%d_CLength",dcname[idc].Data(),layer,wire), Form("Drift length %s Layer%d Wire%d;Drift Length (mm);Counts",dcname[idc].Data(),layer,wire), 1100, -0.1,1.0 );
      //new TH2F( Form("%s_Layer%d_Wire%d_TimevLength",dcname[idc].Data(),layer,wire), Form("Drift time vs. Drift length %s Layer%d Wire%d;Drift time (ns);Drift length (mm)",dcname[idc].Data(),layer,wire), 400, -500, 150, 1100, -0.1, 1.0 );
      //new TH2F( Form("%s_Layer%d_Wire%d_CTimevLength",dcname[idc].Data(),layer,wire), Form("Drift time vs. Drift length %s Layer%d Wire%d;Drift time (ns);Drift length (mm)",dcname[idc].Data(),layer,wire), 400, -500, 150, 1100, -0.1, 1.0 );
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
  EventAnalysisMyCDC *event = new EventAnalysisMyCDC();
  return (EventTemp*)event;
}
