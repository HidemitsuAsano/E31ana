#include "File.h"
#include "ConfMan.h"
#include "TKO.h"
#include "CDSHitMan.h"
#include "CDSTrackingMan.h"
#include "CircleFit.h"
#include "HelixFit.h"
#include "BeamLineHitMan.h"
#include "BeamLineTrackMan.h"
#include "EventHeader.h"
#include "ScalerMan.h"

#include "EventAlloc.h"
#include "EventTemp.h"
#include "Display.h"
#include "Tools.h"

class EventAnalysis: public EventTemp
{
public:
  EventAnalysis();
  ~EventAnalysis();
private:
  TFile *rtFile;
  TFile *cdcFile;
  TTree *cdcTree;
  TTree *evTree;
  TTree *scaTree;
  CDSHitMan *cdsMan;
  CDSTrackingMan *trackMan;
  EventHeader *header;
  EventHeader *header2;
  int t0,t1;
  
  int AllGoodTrack;
  int nTrack;
  int CDC_Event_Number;
  
public:
  void Initialize( ConfMan *conf );
  void InitializeHistogram();
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

const int MaxTreeSize = 100000000000;
void EventAnalysis::Initialize( ConfMan *conf )
{
#if 1
  std::cout << " Enter EventAnalysis::Initialize " << std::endl;
#endif
  confMan = conf;
  cdcFile = new TFile(confMan->GetCDSTrackFileName().c_str());
  if(!cdcFile->IsOpen()){
    std::cout<<" failed to open " <<confMan->GetCDSTrackFileName().c_str()<< "  !!!"<<std::endl;
    exit(false);
  }

  cdcTree=(TTree*)cdcFile->Get("EventTree");
  cdcTree->SetBranchAddress( "CDSTrackingMan", &trackMan );
  cdcTree->SetBranchAddress( "EventHeader" ,&header2)
;  
  rtFile = new TFile( confMan->GetOutFileName().c_str(), "recreate" );
  rtFile->cd();
  InitializeHistogram();

  evTree = new TTree( "EventTree", "EventTree" );
  header = new EventHeader();
  if( header==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  evTree->Branch( "EventHeader", &header );
  cdsMan = new CDSHitMan();
  if( cdsMan==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  //  evTree->Branch( "CDSHitMan", &cdsMan );
  //  trackMan = new CDSTrackingMan();  
  //  if( trackMan==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  //  evTree->Branch( "CDSTrackingMan", &trackMan );

  t0=clock();
  AllGoodTrack=0;
  nTrack=0;
  CDC_Event_Number=0;
}

void EventAnalysis::USca( int nsca, unsigned int *sca )
{
#if 0
  std::cout << " Enter EventAnalysis::USca " << std::endl;
#endif
  Block_Event_Number++;
  header->SetBlockEventNumber( Block_Event_Number );
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

  if( Event_Number%1000==0)
    {
      t1=clock();
      std::cout << " Event# : " << Event_Number << "  BlockEvent# : " << Block_Event_Number << " AllTrack# : " << nTrack <<" GoodTrack# : " << AllGoodTrack << " Time (s): " << (t1-t0)*1.0e-6 << std::endl;
    }
  
  if(CDC_Event_Number>=cdcTree->GetEntries()) return false;
  cdcTree->GetEntry(CDC_Event_Number);
  //  std::cout<<header2->ev()<<"\t"<<Event_Number<<std::endl;
  if(header2->ev()!=Event_Number) return true;
  CDC_Event_Number++;
  header->SetRunNumber(0);
  header->SetEventNumber(Event_Number);
  
  header->Convert( tko, confMan );

  if( header->IsTrig(Trig_Cosmic) ){
    header->Clear();
    cdsMan->Clear();
    return true;
  }
  cdsMan->Convert( tko, confMan );
  trackMan->Calc( cdsMan, confMan);  

  int nGoodTrack=trackMan->nGoodTrack();
  int nallTrack=  trackMan->nTrack();
  AllGoodTrack+=nGoodTrack;
  nTrack+=nallTrack;
  if( nGoodTrack<1 ){
    header->Clear();
    cdsMan->Clear();
    return true;
  }

  for( int it=0; it<trackMan->nGoodTrack(); it++ ){
    CDSTrack *track = trackMan->Track( trackMan->GoodTrackID(it) );
    bool SINGLE=true;
    for(int layer=1;layer<=NumOfCDCLayers;layer++)
      if(track->nTrackHit(layer)!=1) SINGLE=false;
    double chi = track->Chi();
    double param[5];
    track->GetParameters(param);
    //      if(param[2]==0  ) continue;
    double drho=param[0], phi0=param[1], rho=1./param[2], dz=param[3], tlam=param[4];
    double mom = track->Momentum(); 
    
    //#####################
    //  Fill for CDC XT
    //#####################
    
    if(SINGLE	      
       &&fabs(param[3])<12
       &&fabs(drho)<6
       &&rho>0)
      for(int layer=1;layer<=NumOfCDCLayers;layer++)
	{
	  for(int nhit=0;nhit<track->nTrackHit(layer);nhit++)
	    {
	      CDCHit *cdc=track->TrackHit(cdsMan,layer,nhit);
	      double resi=cdc->resl();
	      double dt=cdc->dt();
	      double dl=cdc->dl();
	      double dlr=dl-resi;
	      int wire=cdc->wire();
	      
	      Tools::Fill1D(Form("CDCdt%d_%d",layer,wire), dt );
	      Tools::Fill1D(Form("CDCresid%d",layer), resi );		
	      Tools::Fill2D(Form("CDCdt_resid%d",layer), dt,resi ); 
	      Tools::Fill2D(Form("CDCdt_resid%d_%d",layer,wire),  dt,resi );
	      Tools::Fill1D(Form("CDCresid%d_%d",layer,wire),  resi );
	      Tools::Fill2D(Form("CDCdt_dl%d_%d",layer,wire),  dt,dlr );
	    }//track nhit
	} //layer
    // finish CDC XT histo filling
  }// cdc itrack
  evTree->Fill();

  header->Clear();
  cdsMan->Clear();
  return true;
}

void EventAnalysis::Finalize()
{
  std::cout << " Enter EventAnalysis::Finalize " << std::endl;

  rtFile->cd();
  //  confMan->SaveCDSParam();
  gFile->Write();
  gFile->Close();
  cdcFile->Close();

  delete cdsMan;
  delete header;
}

void EventAnalysis::InitializeHistogram()
{
  for(int layer=1;layer<=NumOfCDCLayers;layer++)
    { 
      Tools::newTH1F( Form("CDCresid%d",layer),        200,-0.1,0.1);
      Tools::newTH1F( Form("CDCresid_pi%d",layer),     200,-0.1,0.1);
      Tools::newTH1F( Form("CDCresid_proton%d",layer), 200,-0.1,0.1);
      Tools::newTH2F( Form("CDCdt_resid%d",layer),     300,-20,280,200,-0.1,0.1);
      for(int wire=1;wire<=NumOfCDCWiresInLayer[layer-1];wire++)
	{
	  //	  Tools::newTH1F( Form("CDCdt%d_%d",layer,wire),Form("CDCdt%d%d_wire",layer,wire) ,300,-60,240);
	  Tools::newTH1F( Form("CDCdt%d_%d",layer,wire),       600,-200,400);
	  Tools::newTH1F( Form("CDCresid%d_%d",layer,wire),    200,-0.2,0.2);
	  Tools::newTH2F( Form("CDCdt_resid%d_%d",layer,wire), 300,-20,280,200,-0.1,0.1);
	  Tools::newTH2F( Form("CDCdt_dl%d_%d",layer,wire),    300,-20,280,1000,-0.1,0.9);
	}
    }
}
EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysis *event = new EventAnalysis();
  return (EventTemp*)event;
}
