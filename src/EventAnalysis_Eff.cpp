#include "File.h"
#include "ConfMan.h"
#include "TKO.h"
#include "CDSHitMan.h"
#include "CDSTrackingMan.h"
#include "CircleFit.h"
#include "HelixFit.h"
#include "BeamLineHitMan.h"
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
  EventHeader *header;
  ScalerMan *scaMan;

  int counter[5];
  int AllGoodTrack;
  double eff[15];
  double err[15];
  int ntrack[15];

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

const int MaxTreeSize = 1900000000000;
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
  scaMan = new ScalerMan();
  if( scaMan==NULL ){ std::cerr << "!!!!" << std::endl; return; }
  scaTree->Branch( "ScalerMan", &scaMan );


  rtFile->cd();

  new TH1F( "AllTrack", "AllTrack", 15, 0, 16 );
  new TH1F( "HitinLayer", "HitinLayer", 15, 0, 16 );
  //  new TH1F( "EffCDC", "EffCDC", 15, 0, 16 );


  for(int layer=1;layer<=15;layer++)
    {

      new TH1F( Form("ResidinLayer%d",layer),Form("ResidinLayer%d",layer), 3000, -0.15 , 0.15 );
     new TH1F(Form("dlinLayer%d",layer),Form("Drift length in Layer%d",layer),
2000,0,2);    
      new TH2F( Form("dt_dl%d",layer),Form("dt vs dl in Layer%d",layer), 
		700, -50 ,300,120,-0.2, 1.0 );
    }
  AllGoodTrack=0;
  for(int n=0;n<15;n++)
    { 
      eff[n]=0;
      err[n]=0;
      ntrack[n]=0;
    }
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

  rtFile->cd();
  TH1F *h1;
  TH2F *h2;

  Event_Number++;
  { int status = confMan->CheckEvNum( Event_Number, Block_Event_Number );
    if( status==1 ) return true;
    if( status==2 ) return false; }

  if( Event_Number%100==0 )
    {
    std::cout << " Event# : " << Event_Number << "  BlockEvent# : " << Block_Event_Number << " GoodTrack# : " << AllGoodTrack << std::endl;
    for(int n=0;n<15;n++)
      {
	std::cout <<"layer["<<n+1<<"] ineff=  "<<ntrack[n]-eff[n]<< " / " <<ntrack[n]  <<std::endl;
      }

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

  for(int layer=1;layer<=15;layer++)
    {

      trackMan->SetKilledLayer(layer); 
      trackMan->Execute(cdsMan,confMan);
      for(int n=0;n<trackMan->nGoodTrack();n++)
	{
	  bool eff_flag=false;
	  bool trackflag=true;
	  CDSTrack *track=trackMan->Track(trackMan->GoodTrackID(n));
	  double x_c=track->CenterX();
	  double y_c=track->CenterY();
	  double rho=track->Rho();
	  double param[5];
	  track->GetParameters(param);
	  //	  if(fabs(param[2])>0.001 || fabs(param[0])>10) continue;
	  if(fabs(param[3]) >35 || fabs(param[0])>10 ) trackflag=false;
	  if(!track->CDHFlag() )  trackflag=false;
	  CDCWireMapMan *wiremap=confMan->GetCDCWireMapManager();
	  double dphi=wiremap->dphi(layer)*TMath::DegToRad();
	  double radius=wiremap->radius(layer);
	  double cellsize=radius*dphi/2.0;
	  if(param[2]==0)
	    {
	      TVector3 Line_o,Line_d;
	      Line_o=track->LineOrigin();
	      Line_d=track->LineDirection();
	      if(fabs(Line_d.z())>0.6 && fabs(param[3])>30)
		{
		      trackflag=false;
		}
	    }
	  for(int cdc=0;cdc<cdsMan->nCDC(layer);cdc++)
	    {

	      TVector3 cdcpos,wposp;
	      CDCHit *CDC=cdsMan->CDC(layer,cdc);
	      double phi=CDC->phi();
	      double def_deg=fabs(phi-track->Theta());
	      if(def_deg>180 ) def_deg=360-def_deg;
	      if(def_deg>90) continue;

	      cdcpos.SetXYZ(CDC->wx(),CDC->wy(),CDC->wz());
	      wposp.SetXYZ(CDC->wxp(),CDC->wyp(),CDC->wzp());
	      double dl=CDC->dl();
	      double dis=0;
	      if(param[2]!=0)
		{
		  if( (1<=layer && layer<=3) || (8<=layer && layer<=9) ||(14<=layer && layer<=15))
		    {
		      double xest,yest;
		      track->PointToCircle(cdcpos.x(),cdcpos.y(),rho,x_c,y_c,dis,xest,yest);         
		    }
		  else
		    {
		      TVector3 dline,lnest,hnest;
		      dline=wposp-cdcpos;
		      
		      double lw=dline.Mag();
		      dline=dline.Unit();
		      double hphi =track->CalcHelixPhi(cdcpos.x(),cdcpos.y());
		      TVector3 hpos=track->GetPosition(hphi);
		      double k=hpos*dline-cdcpos*dline;
		      TVector3 whpos=cdcpos+k*dline;
		      if(!(track->LineToHelix(whpos,dline,param,lnest,hnest,dis) ))
			{
			  std::cout<<"Miss LtH!"<<std::endl;
			  trackflag=false;
			}
		    }
		}
	      else if(param[2]==0)
		{
		  TVector3 Line_o,Line_d;
		  Line_o=track->LineOrigin();
		  Line_d=track->LineDirection();

		  TVector3 dline,lnest,hitpos;
		  dline=wposp-cdcpos;
		  track->LineToLine(Line_o,Line_d,cdcpos,dline,dl,dis,lnest,hitpos);

		  //		  if(fabs(lnest.z()) >25 ) trackflag=false;
		  //		  if(dis>cellsize*1.5 && layer==1&&trackflag) std::cout<<"l dis "<<layer <<" "<<dis<<" "<<lnest.z()<<std::endl;
		}
	      //	      if(layer==15||layer==14)std::cout<<"layer"<<layer<<"  ev "<<Event_Number<<" dis "<<dis<<std::endl;
	      if(dis<cellsize*2)
		{
		  rtFile->cd();
		  h1 = (TH1F*)rtFile->Get(Form("ResidinLayer%d",layer));    
		  h1->Fill((dis-dl) );
		  h1 = (TH1F*)rtFile->Get(Form("dlinLayer%d",layer));    
		  h1->Fill(dis );
		  if(dis<cellsize*2) eff_flag=true;
		}
	    }//CDC
	  if(!trackflag ) continue;

	  if(eff_flag) {
	    eff[layer-1]++;
	    rtFile->cd();
	    h1 = (TH1F*)rtFile->Get("HitinLayer"); h1->Fill(layer );
   	  }
	  else  
	    {
	      //	      std::cout<<"layer"<<layer<<" NOHit  ev "<<Event_Number<<std::endl;
	    }
	  rtFile->cd();
	  h1 = (TH1F*)rtFile->Get("AllTrack"); h1->Fill(layer );
	  ntrack[layer-1]++;
	}//track

    }//layer
    /***********************/
    /**********************/


  int nGoodTrack=0;
  nGoodTrack=trackMan->nGoodTrack();
  
  if( nGoodTrack<1 ){
    header->Clear();
    blMan->Clear();
    cdsMan->Clear();
    trackMan->Clear();
    return true;
  }

  AllGoodTrack+=nGoodTrack;
  evTree->Fill();
  header->Clear();
  blMan->Clear();
  cdsMan->Clear();

  trackMan->Clear();

  return true;
}

void EventAnalysis::Finalize()
{
  std::cout << " Enter EventAnalysis::Finalize " << std::endl;


  /*
  ofstream of;
  of.open("eff_tmp.dat");
  of <<"#layer  eff  err ntrack "<< std::endl;
  for(int n=0;n<15;n++)
    {
      eff[n]=eff[n]/(double)ntrack[n];
      err[n]= sqrt(ntrack[n]*eff[n]*(1-eff[n]) )/(double)ntrack[n];
      std::cout << "eff["<<n+1<<"] "<<eff[n] <<" +-  "<<err[n]<<" "<<ntrack[n]<< std::endl;
      of << n+1<<"   "<< eff[n] <<"   "<<err[n]<<"  "<<ntrack[n]<< std::endl;
    }
  of.close();
  */

  TH1F *h1, *h2;
  rtFile->cd();
  h1 = (TH1F*)rtFile->Get("HitinLayer");
  h2 = (TH1F*)rtFile->Get("AllTrack");
  TGraphAsymmErrors *EffCDC=new  TGraphAsymmErrors(h1,h2 );
  EffCDC->SetName("EffCDC");
  EffCDC->Write();
  gFile->Write();
  gFile->Close();

  delete blMan;
  delete cdsMan;
  delete trackMan;
  delete header;
}

EventTemp *EventAlloc::EventAllocator()
{
  EventAnalysis *event = new EventAnalysis();
  return (EventTemp*)event;
}
