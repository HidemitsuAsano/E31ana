#include "File.h"
#include "ConfMan.h"
#include "TKO.h"
#include "CDSHitMan.h"
#include "BeamLineHitMan.h"
#include "EventHeader.h"
#include "ScalerMan.h"
#include "analysis.h"

#include "EventAlloc.h"
#include "EventTemp.h"
#include "Display.h"

class EventAnalysisTKO: public EventTemp
{
public:
  EventAnalysisTKO();
  ~EventAnalysisTKO();
private:
  TFile *rtFile;

public:
  void Initialize( ConfMan *conf );
  void USca( int nsca, unsigned int *sca );
  bool UAna( TKOHitCollection *tko );
  void Finalize();
  void InitializeHistogram();
};

EventAnalysisTKO::EventAnalysisTKO()
  : EventTemp()
{
}

EventAnalysisTKO::~EventAnalysisTKO()
{
}

//const int MaxTreeSize = 1900000000000;
void EventAnalysisTKO::Initialize( ConfMan *conf )
{
#if 1
  std::cout << " Enter EventAnalysisTKO::Initialize " << std::endl;
#endif
  confMan = conf;
  
  rtFile = new TFile( confMan->GetOutFileName().c_str(), "recreate" );
  InitializeHistogram();
}

void EventAnalysisTKO::USca( int nsca, unsigned int *sca )
{
#if 0
  std::cout << " Enter EventAnalysisTKO::USca " << std::endl;
#endif
  Block_Event_Number++;
}

bool EventAnalysisTKO::UAna( TKOHitCollection *tko )
{
#if 0
  std::cout << " Enter EventAnalysisTKO::UAna " << std::endl;
#endif

  Event_Number++;
  
  { int status = confMan->CheckEvNum( Event_Number, Block_Event_Number );
    if( status==1 ) return true;
    if( status==2 ) return false; }
 
  if( Event_Number%5000==0 )
  std::cout << " Event# : " << Event_Number << "  BlockEvent# : " << Block_Event_Number << std::endl;
  
  TH1F* h1;
  TKOHit* hit;
  for(int i=0;i<tko->entries();i++)
    {
      hit=tko->hit(i);
      h1=(TH1F*)gFile->Get(Form("TKOc%dn%da%d",hit->cr(),hit->sl(),hit->ch()));
      h1->Fill(hit->data());
      //      if(hit->cr()==2,hit->sl()==1)
      //	std::cout<<hit->ch()<<"\t"<<hit->data()<<std::endl;
    }
  return true;
}

void EventAnalysisTKO::InitializeHistogram(){
  rtFile->cd();  
  
  const int ncrate=MAXSMP; 
   
  int ic=0;
  for(int isl=1;isl<24;isl++)
    for(int ich=0;ich<64;ich++)
      new TH1F(Form("TKOc%dn%da%d",ic,isl,ich),Form("TKOc%dn%da%d",ic,isl,ich),4096,-0.5,4095.5);
  for(ic=1;ic<ncrate;ic++)
    for(int isl=1;isl<24;isl++)
      for(int ich=0;ich<32;ich++)
	new TH1F(Form("TKOc%dn%da%d",ic,isl,ich),Form("TKOc%dn%da%d",ic,isl,ich),4096,-0.5,4095.5);
}

void EventAnalysisTKO::Finalize()
{
  std::cout << " Enter EventAnalysisTKO::Finalize " << std::endl;
  gFile->Write();
  gFile->Close();
}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysisTKO *event = new EventAnalysisTKO();
  return (EventTemp*)event;
}
