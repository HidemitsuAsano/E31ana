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
#include "MathTools.h"
#include "ELossTools.h"

#include "EventAlloc.h"
#include "EventTemp.h"
#include "Display.h"
#include "Tools.h"
#include "TF1.h"
#include "TGraphErrors.h"

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

		int BeginOfAnEvent();
		bool EndOfAnEvent();
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

void EventAnalysis::Initialize( ConfMan *conf )
{
#if 1
	std::cout << " Enter EventAnalysis::Initialize " << std::endl;
#endif

	confMan = conf;
	TString cdcfname=confMan->GetCDSTrackFileName();
	cdcFile = new TFile(cdcfname);
	if(!cdcFile->IsOpen()){
		std::cout<<" failed to open " <<cdcfname<< "  !!!"<<std::endl;
		exit(false);
	}
	std::cout<<" File Successfully opend !!! " <<cdcfname<<std::endl;
	trackMan   = new CDSTrackingMan();
	header2    = new EventHeader();

	cdcTree=(TTree*)cdcFile->Get("EventTree");
	cdcTree->SetBranchAddress( "CDSTrackingMan", &trackMan );
	cdcTree->SetBranchAddress( "EventHeader" ,&header2);

	rtFile = new TFile( confMan->GetOutFileName().c_str(), "recreate" );
	rtFile->cd();

	evTree     = new TTree( "EventTree", "EventTree" );
	header     = new EventHeader();
	cdsMan     = new CDSHitMan();

	t0=time(0);
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

int EventAnalysis::BeginOfAnEvent()
{
	Event_Number++;
	{ 
		int status = confMan->CheckEvNum( Event_Number, Block_Event_Number );
		if(status==1||status==2) return status;
	}
	if( Event_Number%10000==0)
	{
		t1=time(0);
		std::cout << " Event# : " << Event_Number << "  BlockEvent# : " << Block_Event_Number << " AllTrack# : " << nTrack <<" GoodTrack# : " << AllGoodTrack << " Time (s): " << (t1-t0) <<"  "<<CDC_Event_Number<<" / "<<cdcTree->GetEntries() << std::endl;

	}
	return 0;
}

bool EventAnalysis::UAna( TKOHitCollection *tko )
{
#if 0
	std::cout << " Enter EventAnalysis::UAna " << std::endl;
#endif

	int status= BeginOfAnEvent();
	if( status==1 ) return true;
	if( status==2 ) return false;

	if(CDC_Event_Number>=cdcTree->GetEntries()) return false;  
	cdcTree->GetEntry(CDC_Event_Number);
	if(header2->ev()!=Event_Number){
		int tmpnum=CDC_Event_Number;
		while(header2->ev()<Event_Number&&tmpnum<=cdcTree->GetEntries()){
			tmpnum++;    
			cdcTree->GetEntry(tmpnum);
		}
		if(header2->ev()!=Event_Number)  return true;
		else CDC_Event_Number=tmpnum;
	}
	CDC_Event_Number++;
	if( !header2->IsTrig(Trig_Cosmic) ) return EndOfAnEvent();

	header->SetRunNumber(0);
	header->SetEventNumber(Event_Number);  
	header->Convert( tko, confMan );

	cdsMan->Convert( tko, confMan );
	trackMan->Calc( cdsMan, confMan );  

	int nGoodTrack=trackMan->nGoodTrack();
	int nallTrack=  trackMan->nTrack();
	AllGoodTrack+=nGoodTrack;
	nTrack+=nallTrack;

	if( nGoodTrack<1 ) return EndOfAnEvent();

	std::vector<int> cdhseg;
	std::vector<int> ihseg;
	int nCDH=0;
	int segCDH=-1;    
	for( int i=0; i<cdsMan->nCDH(); i++ ){
		if( cdsMan->CDH(i)->CheckRange() ){
			nCDH++;
			segCDH=cdsMan->CDH(i)->seg();
			cdhseg.push_back(cdsMan->CDH(i)->seg());
		}
	}

	int nIH=0;
	int segIH=-1;    
	for( int i=0; i<cdsMan->nIH(); i++ ){
		if( cdsMan->IH(i)->CheckRange() ){
			nIH++;
			segIH=cdsMan->IH(i)->seg();
			ihseg.push_back(cdsMan->IH(i)->seg());
		}
	}
	Tools::H1("nTrack",nallTrack, 10, 0, 10 );    
	Tools::H1("nGoodTrack",nGoodTrack, 10, 0, 10 );    
	if(cdhseg.size()==2&&ihseg.size()==2){
		//if(cdhseg.size()==2){
		int tmp=TMath::Abs(cdhseg[0]-cdhseg[1]);
		if(tmp>18) tmp=36-tmp;
		if(tmp>10){
			tmp=TMath::Abs(ihseg[0]-ihseg[1]);
			if(tmp>12) tmp=24-tmp;
			if(tmp>3){
				Tools::H1("nTrackifCosmic",nallTrack, 10, 0, 10 );    
				Tools::H1("nGoodTrackifCosmic",nGoodTrack, 10, 0, 10 );    
				int cut[20];
				for(int i=0;i<20;i++) cut[i]=0;
				for( int it=0; it<nallTrack; it++ ){
					CDSTrack *track = trackMan->Track(it);
					double chi = track->Chi();
					Tools::H1("Chi2",chi, 2000, 0, 200 );    
					for(int i=0;i<20;i++)
						if(chi<(10+10*i)) cut[i]++;
				}
				Tools::H1(Form("Efficiency_Chi2"),0, 21, -5, 205);    
				for(int i=0;i<20;i++){
					if(cut[i]>2) cut[i]=2;
          //	  std::cout<<cut[i]<<std::endl;
          Tools::H1(Form("Efficiency_Chi2"),10*(i+1), 21, -5, 205, cut[i]/2. );    
        }
        //Tracking Efficiency
        Tools::H1(Form("Efficiency_track"),0, 2, 0, 2);    
        if(nGoodTrack >= cdhseg.size())
          Tools::H1(Form("Efficiency_track"),1, 2, 0, 2);    
        //Layer efficiency
        int lay[NumOfCDCLayers];
        for(int layer=1;layer<=NumOfCDCLayers;layer++){
          lay[layer-1] = 0;
        }
        for(int it=0;it<trackMan->nGoodTrack();it++){
          CDSTrack *track = trackMan->Track( trackMan->GoodTrackID(it) );
          Tools::H1(Form("Efficiency_layer"),0, 17, 0, 17);    
          Tools::H1(Form("Effic_Chi2"),0, 200, 0, 200);    
          int tmp[200];
          for(int i=1; i<=200; i++){
            if(track->Chi()<i)
              Tools::H1(Form("Effic_Chi2"),i, 200, 0, 200);    
          }
          for(int layer=1;layer<=NumOfCDCLayers;layer++){
            for(int i=0;i<track->nTrackHit(layer);i++)
              Tools::H1(Form("Efficiency_layer"),layer, 17, 0, 17);    
              lay[layer]++;
          }
        }
        //Layer hit efficiency
        Tools::H1(Form("HitEffic_layer"),0, 17, 0, 17);    
        for(int layer=1;layer<=NumOfCDCLayers;layer++){
          if(lay[layer]>=cdhseg.size())
            Tools::H1(Form("HitEffic_layer"),layer, 17, 0, 17);    
        }
      }
    }
  }

  for( int it=0; it<trackMan->nGoodTrack(); it++ ){
    //    std::cout<<it<<"  /  "<<trackMan->nGoodTrack()<<std::endl;
    CDSTrack *track = trackMan->Track( trackMan->GoodTrackID(it) );
    bool SINGLE=true;
    for(int layer=1;layer<=NumOfCDCLayers;layer++)
      if(track->nTrackHit(layer)!=1) SINGLE=false;
    //    track->HelixFitting(cdsMan);
    //    track->Calc(confMan);
    double phi = track->Phi();
    double chi = track->Chi();
    double param[5];
    track->GetParameters(param);
    double drho=param[0], phi0=param[1], rho=1./param[2], dz=param[3], tlam=param[4];
    double mom = track->Momentum(); 

    Tools::H1(Form("CDCMomentum"), mom , 300,-1.5,1.5);
    Tools::H1(Form("CDCDipAngle"), tlam, 300,-1.5,1.5);
    //#####################
    //  H for CDC XT
    //#####################

    if(SINGLE	      
        //&&fabs(param[3])<12
        //&&fabs(drho)<6
        //&&rho>0)
      )
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
            Tools::H1(Form("CDCresid%d"   ,layer), resi    , 200,-0.1,0.1);		
            Tools::H2(Form("CDCdt_resid%d",layer), dt ,resi, 300,-20,280,200,-0.1,0.1); 
            Tools::H1(Form("CDCdt%d_%d"      ,layer,wire), dt   ,600,-200,400);
            Tools::H1(Form("CDCdt%d"      ,layer), dt   ,600,-200,400);
            Tools::H1(Form("CDCHitPat%d"      ,layer), wire   ,400,0,400);
            Tools::H1(Form("CDCresid%d_%d"   ,layer,wire), resi ,200,-0.2,0.2);
            Tools::H2(Form("CDCdt_resid%d_%d",layer,wire), dt ,resi ,300,-20,280,200,-0.1,0.1);
            Tools::H2(Form("CDCdt_dlr%d_%d"   ,layer,wire), dt ,dlr  ,300,-20,280,200,-0.1,0.9);
            Tools::H2(Form("CDCdt_dlr%d"   ,layer), dt ,dlr  ,300,-20,280,200,-0.1,0.9);
            Tools::H2(Form("CDCdt_dl%d_%d"   ,layer,wire), dt ,dl  ,300,-20,280,200,-0.1,0.9);
            Tools::H2(Form("CDCdt_dl%d"   ,layer), dt ,dl  ,300,-20,280,200,-0.1,0.9);
          }//track nhit
        } //layer
  }// cdc itrack  
  return EndOfAnEvent();
  }

  bool EventAnalysis::EndOfAnEvent(){
    header->Clear();
    cdsMan->Clear();
    return true;
  }
  void EventAnalysis::Finalize()
  {
    std::cout << " Enter EventAnalysis::Finalize " << std::endl;

    rtFile->cd();

    //########  Calc Resid  ###############//  

    TF1 *gaus1=new TF1("gaus1","gaus(0)");
    gaus1->SetParameters(100,0,0.02);
    double nlayer[15]={0};
    double resid[15]={0};
    double ex[15]={0};
    double ey[15]={0};
    for(int layer=1;layer<=15;layer++)
    {
      //      cout<<"Layer "<<layer<<endl;
      TH1F *h1=(TH1F*)gFile->Get(Form("CDCresid%d",layer) );
      if(!h1) continue;
      h1->Fit("gaus1");
      double par[3]; 
      gaus1->GetParameters(par);
      nlayer[layer-1]=layer;
      resid[layer-1]=par[2];
      ey[layer-1]=gaus1->GetParError(2);      
      std::cout<<"layer "<<layer<<"   residual(cm): "<<par[2]<<std::endl; 
    }
    TGraphErrors *g1=new TGraphErrors(15,nlayer,resid,ex,ey);
    g1->SetLineColor(2);
    //  g1->Draw("ALP");
    g1->SetName("CDCResid");
    g1->Write();

    //##################

    gFile->Write();
    gFile->Close();
    cdcFile->Close();

    delete cdsMan;
    delete header;
  }

  EventTemp *EventAlloc::EventAllocator()
  {
    EventAnalysis *event = new EventAnalysis();
    return (EventTemp*)event;
  }
