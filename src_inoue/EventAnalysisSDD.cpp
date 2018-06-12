#include <fstream>

#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TVirtualPad.h>
#include <TRint.h>

#include "File.h"
#include "ConfMan.h"
#include "TKO.h"
#include "CDSHitMan.h"
#include "BeamLineHitMan.h"
#include "SDDHitMan.h"
#include "EventHeader.h"
#include "ScalerMan.h"
#include "BeamLineTrackMan.h"
#include "EventAlloc.h"
#include "EventTemp.h"
#include "Display.h"
#include "BLInfo.h"

#define TREE 1
#define META 0
#define SR 1

class EventAnalysisSDD: public EventTemp
{
public:
  EventAnalysisSDD();
  ~EventAnalysisSDD();
private:
  TFile *rtFile;
  TFile *vmeFile;
  TTree *evTree;
  TTree *scaTree;
  TTree *vtree;

  BeamLineHitMan *blMan;
  BeamLineTrackMan *trackMan;
  SDDHitMan *sddMan;
  EventHeader *header;
  ScalerMan *scaMan;
  EventStruct veve;
  BLInfo *blInfo;

  int Time;

public:
  void Initialize( ConfMan *conf );
  void USca( int nsca, unsigned int *sca );
  void UTime( int time ){ Time = time; }
  bool UAna( TKOHitCollection *tko );
  void Finalize();
  
  void InitializeHistogram();
};

EventAnalysisSDD::EventAnalysisSDD()
  : EventTemp()
{
}

EventAnalysisSDD::~EventAnalysisSDD()
{
}

const int MaxTreeSize = 1900000000000;
void EventAnalysisSDD::Initialize( ConfMan *conf )
{
#if 0
  std::cout << " Enter EventAnalysisSDD::Initialize " << std::endl;
#endif
  confMan= conf;
#if 1
  int run;
  char* str;
  if(str=getenv("RUN"))
    run=atoi(str);
  else run=4; 
  
  
  if(run<=4) run=1;
  else if(run<=7) run=2;
  else if(run<=9) run=3;
  else if(run<=12) run=4;
  else if(run<=13) run=5;
  else if(run<=15) run=6;
  else if(run<=19) run=7;
  else if(run<=23) run=8;
  else if(run<=25) run=9;
  else if(run<=26) run=8;
  else if(run<=27) run=10;
  else if(run<=31) run=11;
  else if(run<=35) run=12;
  else if(run<=37) run=13;
  
  /*
    ifstream in("param/Nov2010SDD/SDDGainMeta.param");
  ifstream in2("param/Nov2010SDD/SDDVGainMeta.param");
  ifstream in3("param/Nov2010SDD/SDDFGainMeta.param");
  */
  ifstream in("param/Nov2010SDD/tkometaparam.txt");
  ifstream in2("param/Nov2010SDD/vmemetaparam.txt");
  ifstream in3("param/Nov2010SDD/fadcmetaparam.txt");
  ifstream in4("param/Nov2010SDD/srparam.txt");
  double gain,offset;
  int temprun,tempsdd;

  while( in>>temprun>>tempsdd>>gain>>offset){
    if(tempsdd==8) continue;
    if(temprun==run){
      confMan->GetGainMapManager()->SetParam(3,15, 7+tempsdd, gain, offset);
      std::cout<<gain<<"\t"<<offset<<std::endl;
    }
  }
  in.close();
  std::cout<<"SDD TKO ADC Gain has changed"<<std::endl;

  while( in4>>temprun>>tempsdd>>gain>>offset){
    if(temprun==run){
      confMan->GetGainMapManager()->SetParam(3,15, 7+tempsdd, gain, offset);
      std::cout<<gain<<"\t"<<offset<<std::endl;
    }
  }
  in4.close();  
  std::cout<<"SDD TKO ADC Gain has changed"<<std::endl;
  

  while( in2>>temprun>>tempsdd>>gain>>offset){
    if(temprun==run){
      confMan->GetGainMapManager()->SetParam( 9, 1, tempsdd, gain, offset);
      std::cout<<gain<<"\t"<<offset<<std::endl;
    }
  }
  in2.close();
  std::cout<<"SDD VME PADC Gain has changed"<<std::endl;
  while( in3>>temprun>>tempsdd>>gain>>offset){
    if(temprun==run){
      confMan->GetGainMapManager()->SetParam( 9, 2, tempsdd, gain, offset);
      std::cout<<gain<<"\t"<<offset<<std::endl;
    }
  }
  in3.close();
  std::cout<<"SDD VME FADC Gain has changed"<<std::endl;


  ifstream in5("param/Nov2010SDD/tkoparam.txt");
  //  ifstream in6("param/Nov2010SDD/vmeparam.txt");
  //  ifstream in7("param/Nov2010SDD/fadcparam.txt");
  if(str=getenv("RUN"))
    run=atoi(str);

  while( in5>>temprun>>tempsdd>>gain>>offset){
    if(tempsdd==8){
      if(temprun==run){
	confMan->GetGainMapManager()->SetParam(3,15, 7+tempsdd, gain, offset);
	std::cout<<gain<<"\t"<<offset<<std::endl;
      }
    }
  }
  in5.close();
#endif


  rtFile = new TFile( confMan->GetOutFileName().c_str(), "recreate" );

  InitializeHistogram();
  rtFile->cd();

  evTree= new TTree( "EventTree","EventTree");
  scaTree= new TTree("ScalerTree", "ScalerTree");


  header = new EventHeader();
  if( header==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  evTree->Branch("EventHeader", &header);
  blMan = new BeamLineHitMan();
  if( blMan ==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  //  evTree->Branch("BeamLineHitMan", &blMan);
  blInfo = new BLInfo();
  if( blInfo ==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  evTree->Branch("BLInfo", &blInfo);
  trackMan = new BeamLineTrackMan();
  if( trackMan ==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  //  evTree->Branch("BeamLineTrackMan", &trackMan);
  sddMan = new SDDHitMan();
  if( sddMan ==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  evTree->Branch("SDDHitMan", &sddMan);
  scaMan = new ScalerMan();
  if( scaMan==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  scaTree->Branch("ScalerMan", &scaMan);

  if(confMan->GetVMEFileName().c_str()!=DefaultFileName){
    vmeFile = new TFile( confMan->GetVMEFileName().c_str());
    if(vmeFile->IsOpen()){
      vtree = (TTree*)vmeFile->Get("tr");
      vtree->SetBranchAddress("EventStruct",&veve);
      std::cout<<vtree->GetEntries()<<std::endl;
    }else{
      std::cout<<"------------no VME data available---------------"<<std::endl;
    }
  }else{
    std::cout<<"------------no VME data available---------------"<<std::endl;
    vmeFile=0;
  }
  rtFile->cd();
}

void EventAnalysisSDD::USca( int nsca, unsigned int *sca )
{
#if 0
  std::cout << " Enter EventAnalysisSDD::USca " << std::endl;
#endif
  Block_Event_Number++;
  header->SetBlockEventNumber( Block_Event_Number );
  scaMan->SetBlockEventNumber( Block_Event_Number );
  for( int i=0; i<nsca; i++ ){
    scaMan->AddHit( sca[i], confMan->GetScalerMapManager()->GetName(i) );
  }
#if 0
  std::cout << nsca << std::endl;
  for( int i=0; i<nsca; i++ ){
    std::cout << "  " << sca[i];
  }
  std::cout << std::endl;
#endif
  scaTree->Fill();
  scaMan->Clear();
#if 0
  std::cout << " Finish EventAnalysisSDD::USca " << std::endl;
#endif

}

bool EventAnalysisSDD::UAna( TKOHitCollection *tko )
{
#if 0
  std::cout << " Enter EventAnalysisSDD::UAna " << std::endl;
#endif

  { int status = confMan->CheckEvNum( Event_Number, Block_Event_Number );
    if( status==1 ) return true;
    if( status==2 ) return false; }
  if( Event_Number>10 ) return false;
  if( Event_Number%5000==1 )
    std::cout << " Event# : " << Event_Number << "  BlockEvent# : " << Block_Event_Number << std::endl;
  if(getenv("RUN")!=NULL)  header->SetRunNumber(atoi(getenv("RUN")));
  else header->SetRunNumber(-1);

  header->SetEventNumber( Event_Number );
  header->SetVEventNumber( VEvent_Number );
  header->SetBlockEventNumber( Block_Event_Number );
  header->Convert( tko, confMan );

  blMan->Convert( tko, confMan );
  trackMan->DoTracking( blMan, confMan );
  sddMan->Convert( tko, confMan );
  blInfo->Set(blMan);
  blInfo->Set(trackMan);

  if(vmeFile){
    if(vmeFile->IsOpen()){
      if(header->EventMatch(vtree,veve,confMan)){
#if 0
	std::cout<<"EventMatched."<<std::endl;
	std::cout<<header->evidv()<<"\t"<<header->spillv()<<std::endl;
#endif
	VEvent_Number=header->vev();
	header->SetVEventNumber(VEvent_Number);
	vtree->GetEntry(VEvent_Number);
#if 0
	std::cout<<VEvent_Number<<"\t"<<veve.evid_v<<std::endl;
#endif
	header->SetVMEinfo(veve);
	sddMan->SetData(veve);
      }
    }
  }

  if(header->utime()==-1) header->SetUtime(Time);
  sddMan->Convert( tko, confMan );
  VEvent_Number=header->vev();
  header->SetVEventNumber( VEvent_Number );

  rtFile->cd();

#if TREE
  evTree->Fill();
#endif
  header->Clear();
  blMan->Clear();
  sddMan->Clear();
  trackMan->Clear();
  blInfo->Clear();

  Event_Number++;
  VEvent_Number++;
  return true;
}

void EventAnalysisSDD::Finalize()
{
  std::cout << " Enter EventAnalysisSDD::Finalize " << std::endl;

  rtFile->cd();
  gFile->Write();
  gFile->Close();

  delete trackMan;
  delete sddMan;
  delete blMan;
  delete header;
  delete blInfo;
}

void EventAnalysisSDD::InitializeHistogram()
{
  rtFile->cd();
}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysisSDD *event = new EventAnalysisSDD();
  return (EventTemp*)event;
}
