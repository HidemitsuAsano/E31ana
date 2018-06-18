#include "File.h"
#include "ConfMan.h"
#include "TKO.h"
#include "CDSHitMan.h"
#include "BeamLineHitMan.h"
#include "EventHeader.h"
#include "ScalerMan.h"
#include "Tools.h"
#include "DetectorList.h"

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

  BeamLineHitMan *blMan;
  EventHeader *header;

  typedef std::vector<int> SegmentContainer;
  std::map<int,SegmentContainer> Seg;

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

//const int MaxTreeSize = 1900000000000;
void EventAnalysis::Initialize( ConfMan *conf )
{
#if 1
  std::cout << " Enter EventAnalysis::Initialize " << std::endl;
#endif
  confMan = conf;
  
  rtFile = new TFile( confMan->GetOutFileName().c_str(), "recreate" );
  header = new EventHeader();
  if( header==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  blMan = new BeamLineHitMan();
  if( blMan ==NULL ){ std::cerr << "!!!!" << std::endl; return; }
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
  if( Event_Number%5000==0 )
    std::cout << " Event# : " << Event_Number << "  BlockEvent# : " << Block_Event_Number << std::endl;
  
  header->SetRunNumber(0);
  header->SetEventNumber(Event_Number);  
  header->Convert( tko, confMan );

  blMan->Convert( tko, confMan );

  rtFile->cd();

  DetectorList *dlist=DetectorList::GetInstance();

  // ======== //
  // Raw Data //
  // ======== //


  int hodoid[]={CID_BHD,CID_T0
  };
  int nhodo=sizeof(hodoid)/sizeof(int);
  bool CVC=false;
  for(int ihodo=0;ihodo<nhodo;ihodo++){
    const int cid = hodoid[ihodo];
    const char* name= dlist->GetName(cid).data();
    const int nsegs= dlist->GetNsegs(cid);
    int nHodo=0;
    for( int i=0; i<blMan->nHodo(cid); i++ ){
      HodoscopeLikeHit *hit = blMan->Hodoi(cid,i);
      int seg = hit->seg();
      int au = hit->adc(0), ad = hit->adc(1);
      int tu = hit->tdc(0), td = hit->tdc(1);
      Tools::H1( Form("A%sU%d",name,seg), au, 4096,-0.5,4095.5 );
      Tools::H1( Form("A%sD%d",name,seg), ad, 4096,-0.5,4095.5 );
      Tools::H1( Form("T%sU%d",name,seg), tu, 4096,-0.5,4095.5 );
      Tools::H1( Form("T%sD%d",name,seg), td, 4096,-0.5,4095.5 );
      if( hit->CheckRange() ){
	Seg[cid].push_back(seg);
	Tools::H1( Form("AwT%sU%d",name,seg), au, 4096,-0.5,4095.5 );
	Tools::H1( Form("AwT%sD%d",name,seg), ad, 4096,-0.5,4095.5 );
	Tools::H2( Form("AT%sU%d",name,seg), tu,au, 256,-0.5,4095.5, 256,-0.5,4095.5 );
	Tools::H2( Form("AT%sD%d",name,seg), td,ad, 256,-0.5,4095.5, 256,-0.5,4095.5 );
	Tools::H1( Form("CTMean%s%d",name,seg), hit->ctmean(), 2000,-50,150 );
	Tools::H1( Form("CTSub%s%d",name,seg), hit->ctsub(), 2000,-50,150 );
	Tools::H1( Form("HitPat%s",name), seg, nsegs,0.5, nsegs+0.5 );
      }else{
	Tools::H1( Form("AwoT%sU%d",name,seg), au, 4096,-0.5,4095.5 );
	Tools::H1( Form("AwoT%sD%d",name,seg), ad, 4096,-0.5,4095.5 );
      }
    }    
    Tools::H1( Form("Mul%s",name), nHodo, nsegs+1 , -0.5, nsegs+0.5 );
  }

  // ========== //
  // BHD-T0 TOF //
  // ========== //
  if(Seg[CID_T0].size()==1)
    {
      int t0seg=Seg[CID_T0][0];
      double time0=blMan->Hodo(CID_T0,t0seg)->ctmean();
      for( int i=0; i<(int)Seg[CID_BHD].size(); i++ ){
	int bhdseg=Seg[CID_BHD][i];
	double bhdtime=blMan->Hodo(CID_BHD,bhdseg)->ctmean();
	double tof_bhdt0=time0-bhdtime;
	Tools::H1("TOF_BHDT0",tof_bhdt0, 1000,20.,40.);
	Tools::H1(Form("TOF_BHDseg%dT0seg%d",bhdseg,t0seg),tof_bhdt0, 1000,20.,40.);
      }
    }
  header->Clear();
  blMan->Clear();
  Seg.clear();
  return true;
}

void EventAnalysis::Finalize()
{
  std::cout << " Enter EventAnalysis::Finalize " << std::endl;

  rtFile->cd();
  gFile->Write();
  gFile->Close();

  delete blMan;
  delete header;
}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysis *event = new EventAnalysis();
  return (EventTemp*)event;
}
