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

#define TKOHIS 0
#define TREE 0
class EventAnalysisNCPC: public EventTemp
{
public:
  EventAnalysisNCPC();
  ~EventAnalysisNCPC();
private:
  TFile *rtFile;
  TTree *evTree;
  TTree *scaTree;
  BeamLineHitMan *blMan;
  EventHeader *header;
  ScalerMan *scaMan;

  int slayer[2],sseg[2];
  int run;
  double t0,t1;
public:
  void Initialize( ConfMan *conf );
  void USca( int nsca, unsigned int *sca );
  bool UAna( TKOHitCollection *tko );
  void Finalize();

  void InitializeHistogram();
};

EventAnalysisNCPC::EventAnalysisNCPC()
  : EventTemp()
{
}

EventAnalysisNCPC::~EventAnalysisNCPC()
{
}

//const int MaxTreeSize = 190000000000;
void EventAnalysisNCPC::Initialize( ConfMan *conf )
{
#if 1
  std::cout << " Enter EventAnalysisNCPC::Initialize " << std::endl;
#endif
  confMan = conf;  

  ifstream ifs("./NClayerseg.param");
  ifs>>run;
  for(int i=0;i<2;i++)
    ifs>>slayer[i]>>sseg[i];
  for(int i=0;i<2;i++)
    std::cout<<"layer,seg: "<<slayer[i]<<"  "<<sseg[i]<<std::endl;;  
  ifs.close();


  rtFile = new TFile( confMan->GetOutFileName().c_str(), "recreate" );
  InitializeHistogram();
  rtFile ->cd();

  evTree=new TTree("EventTree","EventTree");
  scaTree=new TTree("ScalerTree","ScalerTree");

  header = new EventHeader();
  if( header==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  evTree->Branch("EventHeader",&header);
  //  cdsMan = new CDSHitMan();
  //  if( cdsMan==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  blMan = new BeamLineHitMan();
  if( blMan ==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  evTree->Branch("BeamLineHitMan",&blMan);
  scaMan = new ScalerMan();
  if( scaMan==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  scaTree->Branch("ScalerMan",&scaMan);
  t0=clock();
}

void EventAnalysisNCPC::USca( int nsca, unsigned int *sca )
{
#if 0
  std::cout << " Enter EventAnalysisNCPC::USca " << std::endl;
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
#if 0
  for( int i=0; i<scaMan->nsca(); i++ ){
    int val = scaMan->sca(i)->val();
    TString name = scaMan->sca(i)->name();
    std::cout << std::setw(14) << name << std::setw(10) << val;
    if( (i+1)%8==0 ) std::cout << std::endl;
  }
  std::cout << std::endl;
#endif
  scaTree->Fill();
  scaMan->Clear();
}

bool EventAnalysisNCPC::UAna( TKOHitCollection *tko )
{
#if 0
  std::cout << " Enter EventAnalysisNCPC::UAna " << std::endl;
#endif

  Event_Number++;
  { int status = confMan->CheckEvNum( Event_Number, Block_Event_Number );
    if( status==1 ) return true;
    if( status==2 ) return false; }

  //  if( Event_Number>50000 ) return false;
  t1=clock();
  if( Event_Number%5000==0 )
    std::cout << " Event# : " << Event_Number << "  BlockEvent# : " << Block_Event_Number << " Time (s): " << (t1-t0)*1.0e-6 << std::endl;
  
  header->SetRunNumber(run);
  header->SetEventNumber(Event_Number);

  TH1F *h1;
  TH2F *h2;
  
  header->Convert( tko, confMan );
  blMan->Convert( tko, confMan );
  //  cdsMan->Convert( tko, confMan );

  rtFile->cd();

  const int ncounter=3;
  int cid[ncounter]={14,32,35};
  TString name[ncounter]={"TOF","NC","PC"};
  int NumOfSeg[ncounter]={NumOfTOFstopSegments,
			  NumOfNCSegments,
			  NumOfPCSegments};
  TString ud[2]={"u","d"};
  TString at1[5]={"A","T","Time","CTime","Ene"};
  TString at2[5]={"ADC","TDC","Time","Corrected Time","Energy"};
  
  // Fill Data
  int seg,seg2;
  int ic=0;
  double ctm[2][4],eu[2][4],ed[2][4];
  bool HIT[2]={true,true};
  for( int iset=0; iset<2; iset++ ){
    for( int iseg=0; iseg<4; iseg++ ){
      seg=(slayer[iset]-1)*16+sseg[iset]+iseg;
      HodoscopeLikeHit *hit = blMan->Hodo(cid[ic],seg);
      if(!hit)continue;
      double data[10]={};
      data[0] = hit->adcu(), data[1] = hit->adcd();
      data[2] = hit->tdcu(), data[3] = hit->tdcd();
      data[4] = hit->tu(), data[5] = hit->td();
      data[6] = hit->ctu(), data[7] = hit->ctd();
      data[8] = hit->eu(), data[9] = hit->ed();
      if(data[8]<50||data[9]<50)  HIT[iset]=false;
      for(int iud=0;iud<2;iud++){
	for(int at =0 ;at<5;at++){	  
	  h1 = (TH1F*)gFile->Get( Form("%s%s%s%d",at1[at].Data(),name[ic].Data(),ud[iud].Data(),seg) ); h1->Fill( data[at*2 + iud]  );
	}
	if( hit->CheckRange() ){
	  h1 = (TH1F*)gFile->Get( Form("AwT%s%s%d",name[ic].Data(),ud[iud].Data(),seg) ); h1->Fill( data[iud] );
	  h2 = (TH2F*)gFile->Get( Form("AT%s%s%d",name[ic].Data(),ud[iud].Data(),seg) );  h2->Fill( data[iud], data[2+iud] );
	  h2 = (TH2F*)gFile->Get( Form("ETime%s%s%d",name[ic].Data(),ud[iud].Data(),seg) );  h2->Fill( data[6+iud], data[4+iud] );
	}else{
	  HIT[iset]=false;
	}	
      }
      if(!HIT[iset]) continue;
      ctm[iset][iseg]=hit->ctmean();
      eu[iset][iseg]=hit->eu();
      ed[iset][iseg]=hit->ed();
      h1 = (TH1F*)gFile->Get( Form("TimeMean%s%d",name[ic].Data(),seg));
      h1->Fill(hit->tmean());
      h1 = (TH1F*)gFile->Get( Form("TimeSub%s%d",name[ic].Data(),seg));
      h1->Fill(hit->tsub());
      h1 = (TH1F*)gFile->Get( Form("CTimeMean%s%d",name[ic].Data(),seg));
      h1->Fill(hit->ctmean());
      h1 = (TH1F*)gFile->Get( Form("CTimeSub%s%d",name[ic].Data(),seg));
      h1->Fill(hit->ctsub());
    }
    if(HIT[iset]){
      for(int iseg=0;iseg<4;iseg++){
	seg=(slayer[iset]-1)*16+sseg[iset]+iseg;
	for(int iseg2=0;iseg2<4;iseg2++){
	  seg2=(slayer[iset]-1)*16+sseg[iset]+iseg2;
	  if(iseg==iseg2) continue;
	  h2=(TH2F*)gFile->Get( Form("E%sof%s%dTOF_%s%d_%dcorr",
				     ud[0].Data(),name[ic].Data(),
				     seg,name[ic].Data(),seg,seg2));
	  h2->Fill(eu[iset][iseg],ctm[iset][iseg]-ctm[iset][iseg2]);
	  h2=(TH2F*)gFile->Get( Form("E%sof%s%dTOF_%s%d_%dcorr",
				     ud[1].Data(),name[ic].Data(),
				     seg,name[ic].Data(),seg,seg2));
	  h2->Fill(ed[iset][iseg],ctm[iset][iseg]-ctm[iset][iseg2]);
	  
	  h1=(TH1F*)gFile->Get( Form("TOF_%s%d_%d",name[ic].Data(),seg,seg2));   
	  h1->Fill(ctm[iset][iseg]-ctm[iset][iseg2]);
	}
      }
    }
  }
#if TREE
  evTree->Fill();
#endif
  header->Clear();
  blMan->Clear();
  return true;
}

void EventAnalysisNCPC::Finalize()
{
  std::cout << " Enter EventAnalysisNCPC::Finalize " << std::endl;
  
  rtFile->cd();
  gFile->Write();
  gFile->Close();

    delete blMan;
  //  delete cdsMan;
    delete header;
  //  delete scaMan;
  //  delete rtFile;
}

void EventAnalysisNCPC::InitializeHistogram()
{
  rtFile->cd();
  

  //  int cid[3]={2,4,14};
  const int ncounter=3;
  TString name[ncounter]={"TOFS","NC","PC"};
  int NumOfSeg[ncounter]={NumOfTOFstopSegments,
			  NumOfNCSegments,
			  NumOfPCSegments};
  TString ud[2]={"u","d"};
  TString at1[5]={"A","T","Time","CTime","Ene"};
  TString at2[5]={"ADC","TDC","Time","Corrected Time","Energy"};
  for(int seg=1;seg<=112;seg++){
    new TH1F(Form("ADCwT%dU",seg),
	     Form("ADCwT%dU",seg),
	     1000,-0.5,999.5);
    new TH1F(Form("ADCwoT%dU",seg),
	     Form("ADCwoT%dU",seg),
	     1000,-0.5,999.5);
  }
  int seg,seg2;
  int ic=0;
  for( int iset=0; iset<2; iset++ ){
    for( int iseg=0; iseg<4; iseg++ ){
      seg=(slayer[iset]-1)*16+sseg[iset]+iseg;
      for(int at =0;at<2;at++){
	new TH2F( Form("%s%s%dUDCorr",at2[at].Data(),name[ic].Data(),seg),   
		  Form("%s%s%d up down correration",at2[at].Data(),name[ic].Data(),seg),
		  4000, 0, 4000, 4000, 0, 4000 );
      }
      for( int iud=0;iud<2;iud++){
	for(int at =0;at<2;at++){
	  new TH1F( Form("%s%s%s%d",at1[at].Data(),name[ic].Data(),ud[iud].Data(),seg),   
		    Form("%s %s%s%d",at2[at].Data(),name[ic].Data(),ud[iud].Data(),seg),
		    4000, 0, 4000 );
	}
	for(int at =2 ;at<4;at++){
	  new TH1F( Form("%s%s%s%d",at1[at].Data(),name[ic].Data(),ud[iud].Data(),seg),   
		    Form("%s %s%s%d",at2[at].Data(),name[ic].Data(),ud[iud].Data(),seg),
		    4000,   -100, 100 );
	}
	for(int at =4 ;at<5;at++){
	  new TH1F( Form("%s%s%s%d",at1[at].Data(),name[ic].Data(),ud[iud].Data(),seg),   
		    Form("%s %s%s%d",at2[at].Data(),name[ic].Data(),ud[iud].Data(),seg),
		    1000,   -50, 300 );
	}
	new TH1F( Form("AwT%s%s%d",name[ic].Data(),ud[iud].Data(),seg),
		  Form("ADC wTDC %s%s%d",name[ic].Data(),ud[iud].Data(),seg), 
		  4000,    0, 4000 );
	new TH2F( Form("AT%s%s%d",name[ic].Data(),ud[iud].Data(),seg),  
		  Form("ADC TDC corr. %s%s%d",name[ic].Data(),ud[iud].Data(),seg),  
		  100,    0, 1000,  100,    0, 1000 );
	new TH2F( Form("ETime%s%s%d",name[ic].Data(),ud[iud].Data(),seg),  
		  Form("Energy Time corr. %s%s%d",name[ic].Data(),ud[iud].Data(),seg),  
		  100,    -50, 300.,  100,   -50, 50 );
      }

      new TH1F( Form("TimeMean%s%d",name[ic].Data(),seg),   
		Form("Time Mean %s%d",name[ic].Data(),seg),   
		1000,   -100, 100 );
      new TH1F( Form("TimeSub%s%d",name[ic].Data(),seg),   
		Form("Time sub u-d %s%d",name[ic].Data(),seg),   
		1000,   -10, 10 );      
      new TH1F( Form("CTimeMean%s%d",name[ic].Data(),seg),   
		Form("Corrected Time Mean %s%d",name[ic].Data(),seg),   
		1000,   -50, 50 );
      new TH1F( Form("CTimeSub%s%d",name[ic].Data(),seg),   
		Form("Corrected Time sub u-d %s%d",name[ic].Data(),seg),   
		1000,   -10, 10 );      
      
      for(int iseg2=0;iseg2<4;iseg2++){
	seg2=(slayer[iset]-1)*16+sseg[iset]+iseg2;
	if(iseg==iseg2) continue;
	for( int iud=0;iud<2;iud++){
	  new TH2F( Form("E%sof%s%dTOF_%s%d_%dcorr",
			 ud[iud].Data(),name[ic].Data(),seg,name[ic].Data(),seg,seg2),  
		    Form("Energy %s of %s%d vs TOF %s%d-%d corr",
			 ud[iud].Data(),name[ic].Data(),seg,name[ic].Data(),seg,seg2),  
		    200,    -50, 300,  1000,  -50, 50 );
	}	
	new TH1F( Form("TOF_%s%d_%d",name[ic].Data(),seg,seg2),  
		  Form("TOF_%s%d-%d",name[ic].Data(),seg,seg2),  
		  2000,    -50,50 );
      }	
    }
  } 
}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysisNCPC *event = new EventAnalysisNCPC();
  return (EventTemp*)event;
}
