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
  CDSTrackingMan *trackMan;
  BeamLineHitMan *blMan;
  BeamLineTrackMan *blTrackMan;
  EventHeader *header;
  ScalerMan *scaMan;

  int AllGoodTrack;;
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

const int MaxTreeSize = 100000000000;
void EventAnalysis::Initialize( ConfMan *conf )
{
#if 1
  std::cout << " Enter EventAnalysis::Initialize " << std::endl;
#endif
  confMan = conf;
  rtFile = new TFile( confMan->GetOutFileName().c_str(), "recreate" );
  rtFile->cd();
  evTree = new TTree( "EventTree", "EventTree" );
  scaTree= new TTree( "ScalerTree", "ScalerTree" );
  header = new EventHeader();
  if( header==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  evTree->Branch( "EventHeader", &header );
  cdsMan = new CDSHitMan();
  if( cdsMan==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  evTree->Branch( "CDSHitMan", &cdsMan );
  trackMan = new CDSTrackingMan();  
  if( trackMan==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  evTree->Branch( "CDSTrackingMan", &trackMan );
  blMan = new BeamLineHitMan();
  if( blMan ==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  evTree->Branch( "BeamLineHitMan", &blMan );
  blTrackMan = new BeamLineTrackMan();
  if(blTrackMan==NULL) {std::cerr << "!!!!" << std::endl; return; }
  evTree->Branch("BeamLineTrackMan", &blTrackMan);
  scaMan = new ScalerMan();
  if( scaMan==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  scaTree->Branch( "ScalerMan", &scaMan );


  AllGoodTrack=0;

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

  if( Event_Number%1000==0 )
    {
    std::cout << " Event# : " << Event_Number << "  BlockEvent# : " << Block_Event_Number << " GoodTrack# : " << AllGoodTrack << std::endl;
    }
  
  header->SetRunNumber(0);
  header->SetEventNumber(Event_Number);
  
  header->Convert( tko, confMan );
  blMan->Convert( tko, confMan );
  cdsMan->Convert( tko, confMan );

  int nT0 = 0;
  for( int i=0; i<blMan->nT0(); i++ ){
    if( blMan->T0(i)->CheckRange() ){
      nT0++;
    }
  }
  
  /***********************/
  /**** Tracking **********/
  /***********************/
				       
  trackMan->Execute(cdsMan,confMan);
  blTrackMan->DoTracking(blMan,confMan);

    /***********************/
    /**********************/


  int nGoodTrack=trackMan->nGoodTrack();
  
  if( nGoodTrack<1 ){
    header->Clear();
    blMan->Clear();
    blTrackMan->Clear();
    cdsMan->Clear();
    trackMan->Clear();
    return true;
  }

  for(int n=0;n<nGoodTrack;n++)
    {
      int GoodTrack=trackMan->GoodTrackID(n);
      trackMan->CalcVertex_beam(GoodTrack,blTrackMan,confMan);
    }
  
  /*  
  if( blTrackMan->ntrackBLC2()<1 ){
    header->Clear();
    blMan->Clear();
    blTrackMan->Clear();
    cdsMan->Clear();
    trackMan->Clear();
    return true;
  }
  */
  //===========================
  /*
  for(int itr=0;itr<blTrackMan->ntrackBLC2();itr++)
    {
      LocalTrack *blc2=blTrackMan->trackBLC2(itr);
      std::cout<<"chi zx : zy : all "<<blc2->chi2xz()<<" : "<<blc2->chi2yz()<<" : "<<blc2->chi2all()<<std::endl;
      for(int ih=0;ih<blc2->nhit();ih++)
	{
	  ChamberLikeHit *hit=blc2->hit(ih);
	  TMarker blt_m;	
	  double hpos,wpos,dltrack;
	  if(hit->xy()==0) 
	    {
	      hpos=hit->x(); wpos=hit->wx();dltrack=fabs(hpos-wpos);
	    }  
	  else if(hit->xy()==1) 
	    {
	      hpos=hit->y(); wpos=hit->wy();dltrack=fabs(hpos-wpos);
	    }  
	  std::cout<<"cid layer wire: dt : dl : resid "
		   <<hit->cid()<<" : "<<hit->layer()
		   <<" : "<<hit->wire()<<" : "<<hit->dt()
		   <<" : "<<dltrack
		   <<" : "<<hit->resl()<<std::endl;
	}
    }
  */
  //==========================

  
  AllGoodTrack+=nGoodTrack;
  evTree->Fill();

  header->Clear();
  blMan->Clear();
  blTrackMan->Clear();
  cdsMan->Clear();
  trackMan->Clear();
  return true;
}

void EventAnalysis::Finalize()
{
  std::cout << " Enter EventAnalysis::Finalize " << std::endl;

  rtFile->cd();
  gFile->Write();
  gFile->Close();

  delete blMan;
  delete blTrackMan;
  delete cdsMan;
  delete trackMan;
  delete header;
}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysis *event = new EventAnalysis();
  return (EventTemp*)event;
}
