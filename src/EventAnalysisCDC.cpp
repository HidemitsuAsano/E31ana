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
#include "GlobalVariables.h"
#include "Tools.h"

class EventAnalsysCDC: public EventTemp
{
public:
  EventAnalsysCDC();
  ~EventAnalsysCDC();
private:
  TFile *rtFile;
  TFile *cdcFile;
  TTree *cdcTree;
  CDSHitMan *cdsMan;
  CDSTrackingMan *trackMan;
  EventHeader *header;
  Time_t Time;
  double scainit[40];
  double scaend[40];
  
  int AllGoodTrack;
  int nTrack;
  int CDC_Event_Number;
  
public:
  void Initialize( ConfMan *conf );
  void InitializeHistogram();
  void USca( int nsca, unsigned int *sca );
  bool UAna( TKOHitCollection *tko );
  void Finalize();
  bool EndOfAnEvent(bool flag=true);
  void UTime( int time ){ Time = time; };
};

EventAnalsysCDC::EventAnalsysCDC()
  : EventTemp()
{
}

EventAnalsysCDC::~EventAnalsysCDC()
{
}

void EventAnalsysCDC::Initialize( ConfMan *conf )
{
  confMan = conf;
  cdcFile = new TFile(confMan->GetCDSTrackFileName().c_str());
  if(!cdcFile->IsOpen()){
    std::cout<<" failed to open " <<confMan->GetCDSTrackFileName().c_str()<< "  !!!"<<std::endl;
    exit(false);
  }
  
  header = new EventHeader();
  if( header==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  cdsMan = new CDSHitMan();
  if( cdsMan==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  trackMan = new CDSTrackingMan();
  if( trackMan ==NULL) { std::cerr << "!!!!" << std::endl; return;}
  
  //open CDC tracking man
  cdcTree=(TTree*)cdcFile->Get("EventTree");
  cdcTree->SetBranchAddress( "EventHeader" ,&header);  
  cdcTree->SetBranchAddress( "CDSTrackingMan", &trackMan );

  rtFile = new TFile( confMan->GetOutFileName().c_str(), "recreate" );
  rtFile->cd();


  AllGoodTrack=0;
  nTrack=0;
  CDC_Event_Number=0;
}

void EventAnalsysCDC::USca( int nsca, unsigned int *sca )
{
#if 0
  std::cout << " Enter EventAnalsysCDC::USca " << std::endl;
#endif
  Block_Event_Number++;
  header->SetBlockEventNumber( Block_Event_Number );
}

bool EventAnalsysCDC::UAna( TKOHitCollection *tko )
{
#if 0
  std::cout << " Enter EventAnalsysCDC::UAna " << std::endl;
#endif

  Event_Number++;
  { int status = confMan->CheckEvNum( Event_Number, Block_Event_Number );
    if( status==1 ) return EndOfAnEvent();
    if( status==2 ) return EndOfAnEvent(false);}
  if( Event_Number%5000==0 )
    std::cout << " Event# : " << Event_Number << "  BlockEvent# : " << Block_Event_Number <<"  " <<Time<<"  "<<scaend[10]<<std::endl;
  

  if(CDC_Event_Number>=cdcTree->GetEntries()) return false;
  cdcTree->GetEntry(CDC_Event_Number);
  if(header->ev()!=Event_Number) return true;
  CDC_Event_Number++;
  header->SetRunNumber(0);
  header->SetEventNumber(Event_Number);
  
  header->Convert( tko, confMan );

  if(header->IsTrig(Trig_Cosmic)) return EndOfAnEvent();
  
  cdsMan->Convert( tko, confMan );
  
  //
  trackMan->Calc( cdsMan, confMan);  

  int nGoodTrack=trackMan->nGoodTrack();
  int nallTrack=  trackMan->nTrack();
  AllGoodTrack+=nGoodTrack;
  nTrack+=nallTrack;
  if( nGoodTrack<1 ) return EndOfAnEvent();  

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

  header->Clear();
  cdsMan->Clear();
  return true;
}

void EventAnalsysCDC::Finalize()
{
  std::cout << " Enter EventAnalsysCDC::Finalize " << std::endl;

  rtFile->cd();
  //  confMan->SaveCDSParam();
  gFile->Write();
  gFile->Close();
  cdcFile->Close();

  delete cdsMan;
  delete header;
}

void EventAnalsysCDC::InitializeHistogram()
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

bool EventAnalsysCDC::EndOfAnEvent(bool flag){
  header->Clear();
  cdsMan->Clear();
  trackMan->Clear();
  return flag;
}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalsysCDC *event = new EventAnalsysCDC();
  return (EventTemp*)event;
}
