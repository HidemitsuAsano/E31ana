#include "File.h"
#include "ConfMan.h"
#include "TKO.h"
#include "CDSHitMan.h"
#include "BeamLineHitMan.h"
#include "BeamLineTrackMan.h"
#include "EventHeader.h"
#include "ScalerMan.h"

#include "EventAlloc.h"
#include "EventTemp.h"
#include "Display.h"
#include "GlobalVariables.h"
#include "Tools.h"

class EventAnalysis: public EventTemp
{
public:
  EventAnalysis();
  ~EventAnalysis();
private:
  TFile *rtFile;
  TTree *evTree;
  BeamLineHitMan *blMan;

  BeamLineTrackMan *trackMan;
  EventHeader *header;
  int Time;

public:
  void Initialize( ConfMan *conf );
  void USca( int nsca, unsigned int *sca );
  bool UAna( TKOHitCollection *tko );
  bool EndOfAnEvent();
  void Finalize();
  void UTime( int time ){ Time = time; };
};

EventAnalysis::EventAnalysis()
  : EventTemp()
{
}

EventAnalysis::~EventAnalysis()
{
}

//const int MaxTreeSize = 19000000000;
void EventAnalysis::Initialize( ConfMan *conf )
{
#if 0
  std::cout << " Enter EventAnalysis::Initialize " << std::endl;
#endif
  confMan = conf;
  rtFile = new TFile( conf->GetOutFileName().c_str() , "recreate" );

  header = new EventHeader();
  if( header==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  blMan = new BeamLineHitMan();
  if( blMan ==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  trackMan = new BeamLineTrackMan();
  if( trackMan ==NULL ){ std::cerr << "!!!!" << std::endl; return; }
}

void EventAnalysis::USca( int nsca, unsigned int *sca )
{
  Block_Event_Number++;
  header->SetBlockEventNumber( Block_Event_Number );
}

bool EventAnalysis::EndOfAnEvent()
{
  header->Clear();
  blMan->Clear();
  trackMan->Clear();
  return true;
}
bool EventAnalysis::UAna( TKOHitCollection *tko )
{
#if 0
  std::cout << " Enter EventAnalysis::UAna " << std::endl;
#endif
  Event_Number++;
  {
    int status = confMan->CheckEvNum( Event_Number, Block_Event_Number );
    if( status==1 ) return true;
    if( status==2 ) return false; 
  }
  
  if( Event_Number%5000==0 )
    std::cout << " Event# : " << Event_Number << "  BlockEvent# : " << Block_Event_Number <<std::endl;

  header->SetRunNumber(0);
  header->SetEventNumber(Event_Number);

  DetectorList *dlist=DetectorList::GetInstance();
  //#########
  header->Convert( tko, confMan );
  blMan->Convert( tko, confMan );
  trackMan->DoTracking( blMan, confMan );
  rtFile->cd();

  if(header->IsTrig(Trig_Cosmic))  return EndOfAnEvent();

  // BeamLine Chamber

  static int dccid[]={CID_BLC1a,CID_BLC1b,
		  CID_BLC2a,CID_BLC2b,
		  CID_BPC,  CID_FDC1
  };
  static double dlmax[]={0.4,0.4,0.25,0.25,0.36,0.6};
  const int ndc=sizeof(dccid)/sizeof(dccid[0]);
  
  for(int idc=0;idc<ndc;idc++){
    const int cid= dccid[idc];
    const char* name= dlist->GetName(cid).data();
    const int nlays= dlist->GetNlayers(cid);
    // raw and converted data
    for( int layer=1; layer<=nlays; layer++ ){
      const int nwire = confMan->GetBLDCWireMapManager()->GetNWire( cid, layer );
      Tools::H1( Form("Multiplicity%s_layer%d",name,layer) , blMan->nBLDC(dccid[idc],layer), nwire+1, -0.5, nwire+0.5 );
      for( int i=0; i<blMan->nBLDC(dccid[idc],layer); i++ ){
	ChamberLikeHit *hit = blMan->BLDC(dccid[idc],layer,i);
	int wire = hit->wire();
	int tdc = hit->tdc();
	int dt = hit->dt();	
	Tools::H1( Form("h_dt_%s_layer%d",name,layer)        , dt, 3000,-200,400);
	Tools::H1( Form("h_tdc_%s_layer%d",name,layer)        , tdc, 2000, 0, 2000);
	Tools::H1( Form("h_dt_%s_layer%d_wire%d",name,layer,wire)        , dt, 3000,-200,400);
	Tools::H1( Form("h_tdc_%s_layer%d_wire%d",name,layer,wire), tdc, 2000, 0, 2000);
	Tools::H1( Form("HitProfile%s_layer%d",name,layer)   , wire, nwire, 0.5, nwire+0.5 );
	if(layer%2==1){
	  for( int i2=0; i2<blMan->nBLDC(cid,layer+1); i2++ ){
	    ChamberLikeHit *hit2 = blMan->BLDC(cid,layer+1,i2);
	    int wire2 = hit2->wire();
	    Tools::H2( Form("WireCorr%s_layer%d_%d",name,layer,layer+1), wire, wire2, 
		       nwire, 0.5, nwire+0.5, nwire, 0.5, nwire+0.5 );
	  }
	}
	if(layer<=nlays/2){
	  int layer2=layer+nlays/2;
	  for( int i2=0; i2<blMan->nBLDC(cid,layer2); i2++ ){
	    ChamberLikeHit *hit2 = blMan->BLDC(cid,layer2,i2);
	    int wire2 = hit2->wire();
	    Tools::H2( Form("WireCorr%s_layer%d_%d",name,layer,layer2), wire, wire2, 
		       nwire, 0.5, nwire+0.5, nwire, 0.5, nwire+0.5 );
	  }
	}
      }// hit
    }// layer

    // tracked data
    int ntr=trackMan->ntrackBLDC(dccid[idc]);
    Tools::H1(Form("nTrack%s",name),ntr,10,-0.5,9.5);
    if(ntr!=1) continue;
    for(int itr=0;itr<ntr;itr++){
      LocalTrack *tr=trackMan->trackBLDC(dccid[idc],itr);
      int fac[2]={1,-1};
      Tools::H1(Form("Chi2%s",name), tr->chi2all(), 1000,0,100); 
      Tools::H2(Form("AB%s",name), tr->dx(), tr->dy(), 100, -0.1, 0.1,100,-0.1,0.1);
      Tools::H2(Form("XY%s",name), tr->x(),tr->y(), 100,-15,15,100,-15,15);
      
      for( int ih=0;ih<tr->nhit();ih++ ){
	ChamberLikeHit *hit=tr->hit(ih);
	int layer=hit->layer();
	double dl=hit->dl(); // cm
	double dt=hit->dt(); // ns
	double dltrack=-999; // cm
	if(hit->xy()) dltrack=hit->y()-hit->wy();
	else          dltrack=hit->x()-hit->wx();
	double resi=hit->resl();// cm
	int lr=hit->leftright();
	Tools::H2(Form("h_dt_dltrack_%s_%d",name,layer), dt,dltrack,
		  200,-50,200,200,-dlmax[idc]*1.2,dlmax[idc]*1.2);
	Tools::H2(Form("h_dl_resi_%s_%d",name,layer), dl*fac[lr],resi,
		  200,-dlmax[idc]*1.2,dlmax[idc]*1.2,200,-0.1,0.1);
	Tools::H1(Form("h_resi_%s_%d",name,layer), resi,200,-0.1,0.1);
      } // track hit
    } // track
  }// mwdc

  return EndOfAnEvent();
}

void EventAnalysis::Finalize()
{
  std::cout << " Enter EventAnalysis::Finalize " << std::endl;
  rtFile->cd();
  gFile->Write();
  gFile->Close();
  delete trackMan;
  delete blMan;
  delete header;
}


EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysis *event = new EventAnalysis();
  return (EventTemp*)event;
}
