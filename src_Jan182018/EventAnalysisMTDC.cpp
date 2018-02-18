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
#include "EventHeader.h"
#include "ScalerMan.h"
#include "BeamLineTrackMan.h"
#include "EventAlloc.h"
#include "EventTemp.h"
#include "Display.h"
#include "ana_common.h"

#define TREE 0

class EventAnalysisMTDC: public EventTemp
{
public:
  EventAnalysisMTDC();
  ~EventAnalysisMTDC();
private:
  TFile *rtFile;
  TFile *mtdcFile;
  TTree *evTree;
  TTree *scaTree;
  TTree *mtdcTree;

  BeamLineHitMan *blMan;
  BeamLineTrackMan *trackMan;
  EventHeader *header;
  ScalerMan *scaMan;
  int t0,t1;

  int spill_gpio,evid_gpio;
  int mtdc[NV1290][32][MAXHITS];

public:
  void Initialize( ConfMan *conf );
  void USca( int nsca, unsigned int *sca );
  //  void UTime( int time ){ Time = time; }
  bool UAna( TKOHitCollection *tko );
  void Finalize();
  
  void InitializeHistogram();
};

EventAnalysisMTDC::EventAnalysisMTDC()
  : EventTemp()
{
}

EventAnalysisMTDC::~EventAnalysisMTDC()
{
}

const int MaxTreeSize = 1900000000000;
void EventAnalysisMTDC::Initialize( ConfMan *conf )
{
#if 0
  std::cout << " Enter EventAnalysisMTDC::Initialize " << std::endl;
#endif
  confMan= conf;
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
  trackMan = new BeamLineTrackMan();
  if( trackMan ==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  //  evTree->Branch("BeamLineTrackMan", &trackMan);
  scaMan = new ScalerMan();
  if( scaMan==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  scaTree->Branch("ScalerMan", &scaMan);

  if(confMan->GetVMEFileName().c_str()!=DefaultFileName){
    mtdcFile = new TFile( confMan->GetVMEFileName().c_str());
    if(mtdcFile->IsOpen()){
      mtdcTree = (TTree*)mtdcFile->Get("tree_mtdc");
      mtdcTree->SetBranchAddress("spill_gpio",&spill_gpio);
      mtdcTree->SetBranchAddress("evid_gpio",&evid_gpio);
      mtdcTree->SetBranchAddress("mtdc",mtdc);
      std::cout<< confMan->GetVMEFileName().c_str() << " opened"<<std::endl;
      std::cout<<mtdcTree->GetEntries()<<std::endl;
    }else{
      std::cout<<"------------no VME data available---------------"<<std::endl;
    }
  }else{
    std::cout<<"------------no VME data available---------------"<<std::endl;
    mtdcFile=0;
  }
  rtFile->cd();
}

void EventAnalysisMTDC::USca( int nsca, unsigned int *sca )
{
#if 0
  std::cout << " Enter EventAnalysisMTDC::USca " << std::endl;
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
  //  scaTree->Fill();
  scaMan->Clear();
#if 0
  std::cout << " Finish EventAnalysisMTDC::USca " << std::endl;
#endif

}

#define debug 0
#define match 1
#define raw 0

bool EventAnalysisMTDC::UAna( TKOHitCollection *tko )
{
  //  Event_Number++;
  if(VEvent_Number==-1) VEvent_Number = 0;
  if(Event_Number==-1) Event_Number = 0;

  // { int status = confMan->CheckEvNum( Event_Number, Block_Event_Number );
  //   if( status==1 ) return true;
  //   if( status==2 ) return false; }

  if( Event_Number%5000==0 )
  std::cout << " Event# : " << Event_Number << "  BlockEvent# : " << Block_Event_Number << std::endl;
  if(getenv("RUN")!=NULL)  header->SetRunNumber(atoi(getenv("RUN")));
  else header->SetRunNumber(-1);

   // std::cout<<"NEntries = "<<mtdcTree->GetEntries()<<std::endl;
   // getchar();

  if( VEvent_Number > mtdcTree->GetEntries() ){
    std::cout<<"There are more events from TKO than MTDC."<<std::endl;
    return false;
  }

  // std::cout<<"NEntries = "<<mtdcTree->GetEntries()<<std::endl;
  // getchar();

  header->SetBlockEventNumber( Block_Event_Number );
  header->Convert( tko, confMan );
  blMan->Convert( tko, confMan );

  //  std::cout<<"-----------------A new TKO event-----------------"<<std::endl;

#if raw
  if(!mtdcFile){
    std::cout<<"no mtdc tree file found!!"<<std::endl;
    return false;
  }

  VEvent_Number++;
  mtdcTree->GetEntry(VEvent_Number);
  Event_Number++;
  header->SetEventNumber( Event_Number );
  header->SetVEventNumber(VEvent_Number);
  header->sync_set(true);
#endif


#if match
  //  std::cout<<"-----------------Match Mode-----------------"<<std::endl;
  static int error_counter_tko;
  static int error_counter_mtdc;

  if(!mtdcFile){
    std::cout<<"no mtdc tree file found!!"<<std::endl;
    return false;
  }

  if(mtdcFile->IsOpen()){
    mtdcTree->GetEntry(VEvent_Number); //to match TKO event loop
    if(header->spillt(0)==spill_gpio && header->evidt(0)==evid_gpio ){
#if debug
      std::cout<<"matched events"<<std::endl;
      std::cout<<"tko gpio = "<<header->evidt(0)<<" mtdc gpio = "<<evid_gpio<<std::endl;
      std::cout<<"Event_Number = "<<Event_Number<<" VEvent_Number = "<<VEvent_Number<<std::endl;
#endif
      Event_Number++;
      header->SetEventNumber( Event_Number );
      VEvent_Number++;
      header->SetVEventNumber(VEvent_Number);
      header->sync_set(true);
    }else{
      if(header->evidt(0)<evid_gpio ){
#if debug
	std::cout<<"mtdc gpio jumped, try to solve ..."<<std::endl;
	std::cout<<"tko gpio = "<<header->evidt(0)<<" mtdc gpio = "<<evid_gpio<<std::endl;
	std::cout<<"Event_Number = "<<Event_Number<<" VEvent_Number = "<<VEvent_Number<<std::endl;
#endif
	int gpio_holder_pre = evid_gpio;
	mtdcTree->GetEntry(VEvent_Number+1);
	int gpio_holder_aft = evid_gpio;
	if((gpio_holder_aft < gpio_holder_pre) && (gpio_holder_pre<4095)){
#if debug
	  std::cout<<"mtdc gpio jump, corredted; move on ..."<<std::endl;
	  //	  getchar();
#endif
	  mtdcTree->GetEntry(VEvent_Number);
	  evid_gpio = gpio_holder_aft - 1;
	  Event_Number++;
	  header->SetEventNumber( Event_Number );
	  VEvent_Number++;
	  header->SetVEventNumber(VEvent_Number);
	  header->sync_set(true);
	}else{
#if debug
	  std::cout<<"mtdc gpio jump, can NOT be corredted !!"<<std::endl;
	  //	  getchar();
#endif
	  Event_Number++;
	  header->SetEventNumber( Event_Number );
	  error_counter_mtdc++;
#if debug
	  std::cout<<"error_counter_mtdc = "<<error_counter_mtdc<<std::endl;
	  return true;
#endif
	}
      }
      while(header->evidt(0)>evid_gpio ){
#if debug
	std::cout<<"tko gpio jumped"<<std::endl;
	std::cout<<"tko gpio = "<<header->evidt(0)<<" mtdc gpio = "<<evid_gpio<<std::endl;
	std::cout<<"Event_Number = "<<Event_Number<<" VEvent_Number = "<<VEvent_Number<<std::endl;
#endif
	header->SetVEventNumber(VEvent_Number);
	mtdcTree->GetEntry(VEvent_Number);
	VEvent_Number++;
	error_counter_tko++;
#if debug
	std::cout<<"error_counter_tko = "<<error_counter_tko<<std::endl;
#endif
	if( VEvent_Number > mtdcTree->GetEntries() ){
	  std::cout<<"Reach the end of MTDC data."<<std::endl;
	  return false;
	}

      }
    }

  }

#endif

// #if match
//   if(mtdcFile){
//     if(mtdcFile->IsOpen()){
// #if debug
//       std::cout<<"Matching status = "<< header->EventMatch(mtdcTree,spill_gpio,evid_gpio) <<std::endl;
//       std::cout<<"----------------------------------------------- "<<std::endl;
// #endif
//       if(header->EventMatch(mtdcTree,spill_gpio,evid_gpio) == 1){
// #if debug
// 	std::cout<<"EventMatched."<<std::endl;
// 	std::cout<<header->evidv()<<"\t"<<header->spillv()<<std::endl;
// #endif
// 	mtdcTree->GetEntry(VEvent_Number);
	
//       }else if(header->EventMatch(mtdcTree,spill_gpio,evid_gpio) == 2){
// 	//for MTDC lost events
// 	std::cout<<"Matching failed: MTDC lost events"<<std::endl;
// 	std::cout<<"Increment of TKO for MTDC lost Events"<<std::endl;
// 	std::cout<<"VEvent_Number = "<<VEvent_Number<<std::endl;
// 	std::cout<<"Event_Number = "<<Event_Number<<std::endl;
// 	std::cout<<"TKO: spill = "<<header->spillt()<<" evid = "<<header->evidt()<<std::endl;
// 	std::cout<<"MTDC: spill = "<<spill_gpio<<" evid = "<<evid_gpio<<std::endl;
// 	return true;
//       }else if(header->EventMatch(mtdcTree,spill_gpio,evid_gpio) == 20){
// 	//for MTDC gpio jumping
// 	std::cout<<"Matching failed: MTDC jumping?"<<std::endl;
// 	std::cout<<"VEvent_Number = "<<VEvent_Number<<std::endl;
// 	std::cout<<"Event_Number = "<<Event_Number<<std::endl;
// 	std::cout<<"TKO: spill = "<<header->spillt()<<" evid = "<<header->evidt()<<std::endl;
// 	std::cout<<"MTDC: spill = "<<spill_gpio<<" evid = "<<evid_gpio<<std::endl;
// 	int gpio_evid_counter = evid_gpio;
// 	mtdcTree->GetEntry(VEvent_Number+1);
// 	if(gpio_evid_counter > evid_gpio){
// 	  std::cout<<"Matching failed: MTDC jumping!!"<<std::endl;
// 	  evid_gpio = evid_gpio-1;
// 	  return true;
// 	}else{return true;}
//       }else if(header->EventMatch(mtdcTree,spill_gpio,evid_gpio) == 30){
// 	//for TKO gpio jumping
// 	std::cout<<"Matching failed: TKO jumping?"<<std::endl;
// 	std::cout<<"VEvent_Number = "<<VEvent_Number<<std::endl;
// 	std::cout<<"Event_Number = "<<Event_Number<<std::endl;
// 	std::cout<<"TKO: spill = "<<header->spillt()<<" evid = "<<header->evidt()<<std::endl;
// 	std::cout<<"MTDC: spill = "<<spill_gpio<<" evid = "<<evid_gpio<<std::endl;
// 	return true;
//       }
//       while(header->EventMatch(mtdcTree,spill_gpio,evid_gpio) == 3){
// 	//for TKO lost events

// 	std::cout<<"Matching failed: TKO lost Events"<<std::endl;
// 	std::cout<<"Increment of MTDC for TKO lost Events"<<std::endl;
// 	std::cout<<"VEvent_Number = "<<VEvent_Number<<std::endl;
// 	std::cout<<"Event_Number = "<<Event_Number<<std::endl;
// 	std::cout<<"TKO: spill = "<<header->spillt()<<" evid = "<<header->evidt()<<std::endl;
// 	std::cout<<"MTDC: spill = "<<spill_gpio<<" evid = "<<evid_gpio<<std::endl;

// 	VEvent_Number++;
// 	header->SetVEventNumber(VEvent_Number);
//       }
      
//     }
//   }
// #endif
  
  //  VEvent_Number=header->vev();
  //  header->SetVEventNumber( VEvent_Number );
  rtFile->cd();

  double coeff_t0t[5]={1.05141, 1.05246, 1.05546, 1.06354, 1.06767};
  double coeff_t0b[5]={1.06854, 1.03387, 1.04729, 1.06126, 1.06415};
  double coeff_bhd[20]={2.15105, 2.18187, 2.18187, 2.16472, 2.14575,
			2.16742, 2.14401, 2.18187, 2.17285, 2.15556,
			2.18187, 2.17267, 2.16768, 2.11389, 2.16375,
			2.17758, 2.15489, 2.18187, 2.17623, 2.14666};

  TH1F* h1;
  TH2F* h2;

  if(header->sync()){

    h2=(TH2F*)gFile->Get(Form("gpio_matching")); 
    h2->Fill(header->evidt(0), evid_gpio);
    
    h2=(TH2F*)gFile->Get(Form("gpio_evidt")); 
    h2->Fill(header->ev(), header->evidt(0));

    h2=(TH2F*)gFile->Get(Form("gpio_evidv")); 
    h2->Fill(header->vev(), evid_gpio);

    int tko_cut_h=4080;
    int tko_cut_l=0000;

    for(int ii=0; ii<NV1290; ii++){
      for(int jj=0; jj<32; jj++){
	int nhits=0;
	for(int kk=0; kk<MAXHITS; kk++){
	  if( mtdc[ii][jj][kk] > 0){
	    nhits++;
	    int ref = 0;
	    if(ii == 0){
	      ref = mtdc[ii][15][0];
	    }else {
	      ref = mtdc[ii][31][0];
	    }
	    int tdc = mtdc[ii][jj][kk];
	    h1=(TH1F*)gFile->Get(Form("MTDC%dch%d",ii,jj)); 
	    h1->Fill(ref-tdc);
	    //	    h_mtdc_sub[ii][jj]->Fill(ref-tdc);
	    if(ii==0&&jj<16){
	      if(jj>0&&jj<6){
		int tkotdc=blMan->Hodo(CID_T0,jj)->tdcu();
		if(tkotdc>tko_cut_l&&tkotdc<tko_cut_h){
		  h2=(TH2F*)gFile->Get(Form("TKO_MTDC_T0t%d",jj)); 
		  h2->Fill(tkotdc,ref-mtdc[ii][jj][kk]);
		}
	      }
	      if(jj>8&&jj<14){
		int tkotdc=blMan->Hodo(CID_T0,jj-8)->tdcd();
		if(tkotdc>tko_cut_l&&tkotdc<tko_cut_h){
		  h2=(TH2F*)gFile->Get(Form("TKO_MTDC_T0b%d",jj-8)); 
		  h2->Fill(tkotdc,ref-mtdc[ii][jj][kk]);
		}
	      }
	    }
	    
	    if(ii==1&&jj<20){
	      int bhdtdc=blMan->Hodo(CID_BHD,jj+1)->tdcmean();
	      if(bhdtdc>tko_cut_l&&bhdtdc<tko_cut_h){
		h2=(TH2F*)gFile->Get(Form("TKO_MTDC_BHD%d",jj+1)); 
		h2->Fill(bhdtdc,ref-mtdc[ii][jj][kk]);
	      }
	    }	  

	  }
	}
	h1=(TH1F*)gFile->Get(Form("nhitMTDC%dch%d",ii,jj)); 
	h1->Fill(nhits);
      }
    }

    //for MTDC efficiency
    int mtdc_data[NV1290][32][MAXHITS];
    for(int ii=0; ii<NV1290; ii++){
      for(int jj=0; jj<32; jj++){
	for(int kk=0; kk<MAXHITS; kk++){
	  mtdc_data[ii][jj][kk]=-999;
	}
      }
    }

    for(int ii=0; ii<NV1290; ii++){
      for(int jj=0; jj<32; jj++){
	for(int kk=0; kk<MAXHITS; kk++){
	  if( mtdc[ii][jj][kk] > 0){
	    int ref = 0;
	    if(ii == 0){
	      ref = mtdc[ii][15][0];
	    }else {
	      ref = mtdc[ii][31][0];
	    }
	    mtdc_data[ii][jj][kk] = ref - mtdc[ii][jj][kk];

	    // std::cout<<"mtdc_data["<<ii<<"]"<<"["<<jj<<"]"<<"["<<kk<<"]="<<
	    //   (mtdc_data[ii][jj][kk])<<std::endl;
	    // std::cout<<"-------------------------------------"<<std::endl;

	  }
	}
      }
    }
    //    getchar();
    
    for(int jj=0;jj<16;jj++){
      if(jj>0&&jj<6){
    	int tkotdc=blMan->Hodo(CID_T0,jj)->tdcu();
    	if(tkotdc>tko_cut_l&&tkotdc<tko_cut_h){
    	  h1=(TH1F*)gFile->Get(Form("h_hodo_nevents")); 
    	  h1->Fill(jj);
    	  h1=(TH1F*)gFile->Get(Form("TKO_MTDC_T0t%d_eff",jj)); 
    	  bool flag=false;
    	  for(int kk=0;kk<MAXHITS;kk++){
    	    if(mtdc_data[0][jj][kk]!=-999){
    	      h1->Fill(tkotdc*coeff_t0t[jj-1]+mtdc_data[0][jj][kk]);
    	      flag=true;
	      
	      // std::cout<<"mtdc_data["<<0<<"]"<<"["<<jj<<"]"<<"["<<kk<<"]="<<
	      // 	(mtdc_data[0][jj][kk])<<std::endl;
	      // std::cout<<"tkotdc+mtdc_data["<<0<<"]"<<"["<<jj<<"]"<<"["<<kk<<"]="<<
	      // 	(tkotdc+mtdc_data[0][jj][kk])<<std::endl;

    	    }
    	  }
    	}
      }

      if(jj>8&&jj<14){
    	int tkotdc=blMan->Hodo(CID_T0,jj-8)->tdcd();
    	if(tkotdc>tko_cut_l&&tkotdc<tko_cut_h){
    	  h1=(TH1F*)gFile->Get(Form("h_hodo_nevents")); 
    	  h1->Fill(jj);
    	  h1=(TH1F*)gFile->Get(Form("TKO_MTDC_T0b%d_eff",jj-8)); 
    	  bool flag=false;
    	  for(int kk=0;kk<MAXHITS;kk++){
    	    if(mtdc_data[0][jj][kk]!=-999){
    	      h1->Fill(tkotdc*coeff_t0b[jj-9]+mtdc_data[0][jj][kk]);
    	      flag=true;
    	    }
    	  }
    	}
      }
    }

    for(int jj=0;jj<20;jj++){
      int bhdtdc=blMan->Hodo(CID_BHD,jj+1)->tdcmean();
      if(bhdtdc>tko_cut_l&&bhdtdc<tko_cut_h){
    	h1=(TH1F*)gFile->Get(Form("h_hodo_nevents")); 
    	h1->Fill(jj+100+1);
    	h1=(TH1F*)gFile->Get(Form("TKO_MTDC_BHD%d_eff",jj+1)); 
    	bool flag=false;
	//	std::cout<<"bhdtdc="<<bhdtdc<<std::endl;
    	for(int kk=0;kk<MAXHITS;kk++){
    	  if(mtdc_data[1][jj][kk]!=-999){
    	    h1->Fill(bhdtdc*coeff_bhd[jj]+mtdc_data[1][jj][kk]);
    	    flag=true;

	    // std::cout<<"mtdc_data["<<1<<"]"<<"["<<jj<<"]"<<"["<<kk<<"]="<<
	    //   (mtdc_data[1][jj][kk])<<std::endl;
	    // std::cout<<"bhdtdc+mtdc_data["<<1<<"]"<<"["<<jj<<"]"<<"["<<kk<<"]="<<
	    //   (bhdtdc+mtdc_data[1][jj][kk])<<std::endl;
	    
    	  }
    	}
      }
    }	  
    //    getchar();
  }
#if TREE
  evTree->Fill();
#endif
  header->Clear();
  blMan->Clear();
  return true;
}

void EventAnalysisMTDC::Finalize()
{
  std::cout << " Enter EventAnalysisMTDC::Finalize " << std::endl;

  rtFile->cd();
  gFile->Write();
  gFile->Close();

  delete blMan;
  delete header;
}

void EventAnalysisMTDC::InitializeHistogram()
{
  rtFile->cd();
  //  TH1F *h_hodo_nevents=new TH1F("h_hodo_nevents", "h_hodo_nevents", 20, 0, 20);
  new TH1F(Form("h_hodo_nevents"),Form("h_hodo_nevents"),150,-0.5,149.5);

  new TH2F(Form("gpio_matching"),Form("gpio_matching"),5000,-0.5,4999.5,5000,-0.5,4999.5);
  new TH2F(Form("gpio_evidv"),Form("gpio_evidv"),5000,-0.5,499999.5,5000,-0.5,4999.5);
  new TH2F(Form("gpio_evidt"),Form("gpio_evidt"),5000,-0.5,499999.5,5000,-0.5,4999.5);

  for(int imodule=0;imodule<NV1290;imodule++){
    for(int ich=0;ich<32;ich++){
      new TH1F(Form("MTDC%dch%d",imodule,ich),Form("MTDC%dch%d",imodule,ich),5000,-0.5,24999.5);
      new TH1F(Form("nhitMTDC%dch%d",imodule,ich),Form("nhitMTDC%dch%d",imodule,ich),41,-0.5,40.5);
    }
  }
  for(int seg=1;seg<=5;seg++){
    new TH2F(Form("TKO_MTDC_T0t%d",seg),Form("TKO_MTDC_T0t%d",seg),412,-0.5,4095.5,500,-0.5,24999.5);
    new TH2F(Form("TKO_MTDC_T0b%d",seg),Form("TKO_MTDC_T0b%d",seg),412,-0.5,4095.5,500,-0.5,24999.5);

    new TH1F(Form("TKO_MTDC_T0t%d_eff",seg),Form("TKO_MTDC_T0t_eff%d",seg),5000,-0.5,24999.5);
    new TH1F(Form("TKO_MTDC_T0b%d_eff",seg),Form("TKO_MTDC_T0b%d_eff",seg),5000,-0.5,24999.5);
  }

  for(int seg=1;seg<=20;seg++){
    new TH2F(Form("TKO_MTDC_BHD%d",seg),Form("TKO_MTDC_BHD%d",seg),412,-0.5,4095.5,500,-0.5,24999.5);

    new TH1F(Form("TKO_MTDC_BHD%d_eff",seg),Form("TKO_MTDC_BHD%d_eff",seg),5000,-0.5,24999.5);
  }
}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysisMTDC *event = new EventAnalysisMTDC();
  return (EventTemp*)event;
}
