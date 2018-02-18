#include "File.h"
#include "ConfMan.h"
#include "TKO.h"
#include "CDSHitMan.h"
#include "BeamLineHitMan.h"
#include "BeamLineTrackMan.h"
#include "EventHeader.h"
#include "ScalerMan.h"
#include "Particle.h"

#include "EventAlloc.h"
#include "EventTemp.h"
#include "Display.h" 
#include "ELossTools.h"
#include "TrackTools.h"
#include "MyAnalysisBL.h"

#define Debug 0

class EventAnalysisMyBeam: public EventTemp
{
  public:
    EventAnalysisMyBeam();
    ~EventAnalysisMyBeam();
  private:
    TFile *rtFile;
    TFile *rtFile2;
    TTree *evTree;
    TTree *scaTree;
    CDSHitMan *cdsMan;
    CDSTrackingMan *cdstrackMan;
    EventHeader *cdsheader;
    BeamLineHitMan *blMan;
    BeamLineTrackMan *bltrackMan;
    EventHeader *header;
    Particle* particle;
    ScalerMan *scaMan;
    MyAnalysisBL* blAna;

    int AllGoodTrack;
    int nAllTrack;
    int nTrack;
    int t0, t1;

    int totalkaon;
    int prekaon;

  public:
    void Initialize( ConfMan *conf );
    void USca( int nsca, unsigned int *sca );
    bool UAna( TKOHitCollection *tko );
    void Finalize();
    void Clear();

};

  EventAnalysisMyBeam::EventAnalysisMyBeam()
: EventTemp()
{
}

EventAnalysisMyBeam::~EventAnalysisMyBeam()
{
}

const int MaxTreeSize = 1900000000000;
void EventAnalysisMyBeam::Initialize( ConfMan *conf )
{
  std::cout << " Enter EventAnalysisMyBeam::Initialize " << std::endl;

  confMan = conf;

	totalkaon = 0;

  rtFile = new TFile( confMan->GetOutFileName().c_str(), "recreate" );
  header = new EventHeader();
  if( header==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  cdsMan = new CDSHitMan();
  if( cdsMan==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  blMan = new BeamLineHitMan();
  if( blMan ==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  bltrackMan = new BeamLineTrackMan();
  if( bltrackMan ==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  scaMan = new ScalerMan();
  if( scaMan==NULL ){ std::cerr << "!!!!" << std::endl; return; }

  particle = new Particle();

  std::string ofname = confMan->GetOutFileName();
  ofname.insert(ofname.find(".root"),"_anaBeam");
  rtFile2 =  new TFile( Form("%s",ofname.c_str()), "RECREATE");
  rtFile2->cd();
  blAna = new MyAnalysisBL(rtFile2, confMan);

  t0=time(0);

  std::cout << " End EventAnalysisMyBeam::Initialize " << std::endl;

}

void EventAnalysisMyBeam::USca( int nsca, unsigned int *sca )
{
  Block_Event_Number++;
  header->SetBlockEventNumber( Block_Event_Number );
  scaMan->SetBlockEventNumber( Block_Event_Number );
  for( int i=0; i<nsca; i++ ){
    scaMan->AddHit( sca[i], confMan->GetScalerMapManager()->GetName(i) );
  }

  rtFile2->cd();
  TH1F *h1;
  for( int i=0; i<scaMan->nsca(); i++ ){
    long val=0;
    if(i==0){val=scaMan->sca(i)->val();}
	
    else if(i==10){
		 val=scaMan->sca(i)->val();
		if(totalkaon<val){ totalkaon = val; }
		else if(totalkaon>10000000&&prekaon>val){std::cout << "### !!! ### " << std::endl;}
		else if(val>10000000){std::cout << "### !!! ### " << std::endl;}
		else if(totalkaon-100000000<-10000000){ totalkaon = val; }
		else if(totalkaon-100000000<val+1){ totalkaon = 100000000 + val; }
		else if(totalkaon-110000000<val+1){ totalkaon = 110000000 + val; }
		else if(totalkaon-120000000<val+1){ totalkaon = 120000000 + val; }
		else if(totalkaon-130000000<val+1){ totalkaon = 130000000 + val; }
		else if(totalkaon-140000000<val+1){ totalkaon = 140000000 + val; }
		else{ totalkaon = 150000000 + val; }

		prekaon = val;
			
	}
    else{ val=scaMan->sca(i)->val();}
    TString name = scaMan->sca(i)->name();
    h1 = (TH1F*)gFile->Get("Scaler");
    long tmpval =  h1->GetBinContent(i+1);
    if(val-tmpval<0){
      //if(tmpval>100000){
      //  val = val+10000;
      //}
    }
    h1->SetBinContent(i+1,val);
    if(i==0){
      //std::cout << "F.T. : " << val << std::endl;
    }
    if(i==10){
      h1->SetBinContent(i+1,totalkaon);
      //std::cout << "Kaon : " << totalkaon << std::endl;
      //std::cout << "kaon : " << val << std::endl;
    }
  }

  scaMan->Clear();
}

bool EventAnalysisMyBeam::UAna( TKOHitCollection *tko )
{
  Event_Number++;
  { int status = confMan->CheckEvNum( Event_Number, Block_Event_Number );
    if( status==1 ){
      Clear();
      return true;
    }
    if( status==2 ){
      Clear();
      return false;
    }
  }

  if( Event_Number%10000==0 )
  {
    t1=time(0);
    std::cout << " Event# : " << Event_Number << "  BlockEvent# : " << Block_Event_Number << " Time (s):" << (t1-t0) << std::endl;
  }

  header->SetRunNumber(0);
  header->SetEventNumber(Event_Number);

  header->Convert( tko, confMan );
  blMan->Convert( tko, confMan );
  cdsMan->Convert( tko, confMan );


  rtFile2->cd();
  if(!blAna->DoAnalysis(confMan, header, blMan, bltrackMan, particle)){
    Clear();
    return true;
  }

  Clear();
  return true;
}

void EventAnalysisMyBeam::Clear()
{
  particle->Clear();
  blAna->Clear();
  blMan->Clear();
  bltrackMan->Clear();
  cdsMan->Clear();
  header->Clear();
}

void EventAnalysisMyBeam::Finalize()
{
  std::cout << " Enter EventAnalysisMyBeam::Finalize " << std::endl;

  rtFile->cd();
  gFile->Write();
  gFile->Close();
  rtFile2->cd();
  gFile->Write();
  gFile->Close();

  delete particle;

  delete blAna;

  delete blMan;
  delete bltrackMan;
  delete cdsMan;
  delete header;
}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysisMyBeam *event = new EventAnalysisMyBeam();
  return (EventTemp*)event;
}
